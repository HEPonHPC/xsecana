//
// Created by Derek Doyle on 10/9/21.
//
#include "XSecAna/TemplateFitSignalEstimator.h"

namespace xsec {
    const TH1 *
    TemplateFitSignalEstimator::
    _transpose(const TH1 * h) const {
        if(h->GetDimension()==1) return h;
        else {
            auto ret = new TH2D("", "",
                                h->GetNbinsY(),
                                h->GetYaxis()->GetXbins()->GetArray(),
                                h->GetNbinsX(),
                                h->GetXaxis()->GetXbins()->GetArray());
            for(auto i = 0; i < h->GetNbinsX()+2; i++) {
                for(auto j = 0; j < h->GetNbinsY()+2; j++) {
                    ret->SetBinContent(j, i, h->GetBinContent(i,j));
                    ret->SetBinError(j, i, h->GetBinError(i,j));
                }
            }
            return ret;
        }
    }

    TH1 *
    TemplateFitSignalEstimator::
    _mask_and_flatten(const TH1 * mask, const TH1 * templ) const {
        if(templ->GetDimension()==1) {
            std::cout << "Warning: Attempting to apply to mask 1-dimensional template fit" << std::endl;
            return (TH1*) templ->Clone();
        }
        else if(templ->GetDimension()==2) {
            assert(mask->GetDimension()==1);
            // not including under/overflow. Some may want to?
            auto nunmasked = mask->Integral(0, mask->GetNbinsX()+1);
            auto ret = new TH2D("", "",
                                templ->GetNbinsY()-2, 1, templ->GetNbinsY(),
                                nunmasked-2, 1, nunmasked-1);
            auto ii = 0;
            for(auto i = 0u; i < mask->GetNbinsX()+2; i++) {
                if(mask->GetBinContent(i)) {
                    for(auto j = 1; j <= templ->GetNbinsY(); j++) {
                        ret->SetBinContent(ii, j-1, templ->GetBinContent(i, j));
                    }
                    ii++;
                }
            }
            ret->GetXaxis()->SetTitle("Outer Bins");
            ret->GetYaxis()->SetTitle("Template Bins");
            return ret;
        }
        else {
            assert(mask->GetDimension()==2);
            auto nunmasked = ((TH2*)mask)->Integral(0, mask->GetNbinsX()+1,
                                                    0, mask->GetNbinsY()+1);

            // using under/overflow bins here in the flattened histogram
            // for easy un-packing into 1D eigen array that we'll hand to the fitter
            auto ret = new TH2D("", "",
                                templ->GetNbinsZ()-2, 1, templ->GetNbinsZ(),
                                nunmasked-2, 1, nunmasked-1);
            // Not sure why we have to loop over y first,
            // but its necessary for things to get ordered properly
            // during the un-masking
            auto ii = 0;
            for(auto j = 0u; j < mask->GetNbinsY()+2; j++) {
                for(auto i = 0u; i < mask->GetNbinsX()+2; i++) {
                    if(mask->GetBinContent(i, j)) {
                        for(auto k = 1; k <= templ->GetNbinsZ(); k++) {
                            ret->SetBinContent(k-1, ii, templ->GetBinContent(i, j, k));
                        }
                        ii++;
                    }
                }
            }
            ret->GetXaxis()->SetTitle("Template Bins");
            ret->GetYaxis()->SetTitle("Outer Bins");
            return ret;
        }
    }


    TemplateFitSignalEstimator::
    TemplateFitSignalEstimator(const TH1 * signal_template,
                               const std::map<std::string, const TH1 *> & background_templates,
                               const std::map<std::string, Systematic<TH1>> & systematics,
                               const TH1 * mask,
                               fit::TemplateFitCalculator * fit_calc)
            : fFitCalc(fit_calc),
              fMask(mask) {

        // given original template shapes, determine how to return predictions
        // set number of inner bins
        TH1 * project_predictions_like;
        std::vector<Array> templates_1d(background_templates.size()+1);
        if (signal_template->GetDimension() == 1) {
            project_predictions_like = new TH1D("", "",
                                                1, 0, 1);
        } else if (signal_template->GetDimension() == 2) {
            project_predictions_like = new TH1D("", "",
                                                signal_template->GetNbinsX(),
                                                signal_template->GetXaxis()->GetXbins()->GetArray());
        } else if (signal_template->GetDimension() == 3) {
            project_predictions_like = new TH2D("", "",
                                                signal_template->GetNbinsX(),
                                                signal_template->GetXaxis()->GetXbins()->GetArray(),
                                                signal_template->GetNbinsY(),
                                                signal_template->GetYaxis()->GetXbins()->GetArray());
        }
        else {
            throw std::runtime_error("Template fits with greater than 2 outer dimensions not supported");
        }
        // root won't draw histogram unless there are entries
        project_predictions_like->SetEntries(1);
        fProjectPredictionProps = root::TH1Props(project_predictions_like);

        // masked parameter map
        auto _map = root::MapContentsToEigen(mask);
        fParamMap = fit::detail::ParamMap(_map);

        // construct templates with masked bins removed
        fSignalTemplate = _mask_and_flatten(mask, signal_template);
        for(auto bkgd_template : background_templates) {
            fBackgroundTemplates[bkgd_template.first] = _mask_and_flatten(mask,
                                                                          bkgd_template.second);
        }
        // root stores as row-major. The calculator assumes column major
        // outer bins x inner bins
        fDims = {fSignalTemplate->GetNbinsY()+2, fSignalTemplate->GetNbinsX()+2};
        // create masked systematics
        for(auto syst : systematics) {
            std::vector<const TH1*> masked_shifts(syst.second.GetShifts().size());
            for(auto ishift = 0u; ishift < syst.second.GetShifts().size(); ishift++) {
                masked_shifts[ishift] = _mask_and_flatten(mask, syst.second.GetShifts()[ishift]);
            }
            fSystematics[syst.first] = Systematic<TH1>(syst.second.GetName(),
                                                       masked_shifts,
                                                       syst.second.GetType());
        }

        // initialize all templates as free
        for (auto temp_it: fBackgroundTemplates) {
            fIsFreeTemplate[temp_it.first] = true;
        }

        // vector of flattened templates for the fit calculator
        templates_1d[0] = root::MapContentsToEigen(fSignalTemplate); // flatten

        // calculate total in this loop for later covariance matrix calculations
        fTotalTemplate = (TH1*) fSignalTemplate->Clone();
        auto itemplate = 1u;
        for(auto background_template : background_templates) {
            fComponentLabelIdxMap[itemplate] = background_template.first;
            templates_1d[itemplate] = root::MapContentsToEigen(fBackgroundTemplates.at(background_template.first));
            fTotalTemplate->Add(fBackgroundTemplates.at(background_template.first));
            itemplate++;
        }

        // calculate and store covariance matrices
        fTotalCovariance = Matrix::Zero(fDims[0]*fDims[1], fDims[0]*fDims[1]);
        for(auto systematic : fSystematics) {
            fCovarianceMatrices[systematic.first] = systematic.second.CovarianceMatrix(fTotalTemplate);
            fTotalCovariance += fCovarianceMatrices.at(systematic.first);
        }
        fInverseCovariance = fTotalCovariance.inverse();
        fFitCalc = new fit::TemplateFitCalculator(templates_1d,
                                                  fDims,
                                                  fInverseCovariance);

    }


    void
    TemplateFitSignalEstimator::
    _eval_impl(const Array & data, const Array & error,
                                ArrayRef result, ArrayRef rerror) const {
        throw std::runtime_error(__PRETTY_FUNCTION__);
    }

    TH2D *
    TemplateFitSignalEstimator::
    GetReducedSignalTemplate() const {
        return (TH2D*) fSignalTemplate;
    }

    TH2D *
    TemplateFitSignalEstimator::
    GetReducedTotalTemplate() const {
        return (TH2D*) fTotalTemplate;
    }

    TH2D *
    TemplateFitSignalEstimator::
    GetReducedBackgroundTemplate(const std::string & bkgd_label) const {
        return (TH2D*) fBackgroundTemplates.at(bkgd_label);
    }

    const Systematic<TH1> &
    TemplateFitSignalEstimator::
    GetReducedSystematic(const std::string & systematic_label) const {
        return fSystematics.at(systematic_label);
    }

    TH2D *
    TemplateFitSignalEstimator::
    GetTotalCovariance() const {
        auto nrows = fTotalCovariance.rows();
        auto ncols = fTotalCovariance.rows();
        auto ret = new TH2D("", "Total Covariance",
                            nrows, 0, nrows,
                            ncols, 0, ncols);
        root::FillTH2Contents(ret, fTotalCovariance);
        ret->GetXaxis()->SetTitle("Template Bins");
        ret->GetYaxis()->SetTitle("Template Bins");
        return ret;
    }

    TH2D *
    TemplateFitSignalEstimator::
    GetInverseCovariance() const {
        auto nrows = fTotalCovariance.rows();
        auto ncols = fTotalCovariance.rows();
        auto ret = new TH2D("", "Inverse Total Covariance",
                            nrows, 0, nrows,
                            ncols, 0, ncols);
        root::FillTH2Contents(ret, fInverseCovariance);
        ret->GetXaxis()->SetTitle("Template Bins");
        ret->GetYaxis()->SetTitle("Template Bins");
        return ret;
    }

    TH2D *
    TemplateFitSignalEstimator::
    GetCovariance(const std::string & systematic_name) const {
        auto nrows = fTotalCovariance.rows();
        auto ncols = fTotalCovariance.rows();
        auto ret = new TH2D("", TString::Format("%s Covariance",
                                                systematic_name.c_str()).Data(),
                            nrows, 0, nrows,
                            ncols, 0, ncols);
        root::FillTH2Contents(ret, fCovarianceMatrices.at(systematic_name));
        ret->GetXaxis()->SetTitle("Template Bins");
        ret->GetYaxis()->SetTitle("Template Bins");
        return ret;
    }

    void
    TemplateFitSignalEstimator::
    FixTemplate(const std::string & template_name) {
        fIsFreeTemplate.at(template_name) = true;
    }

    TH1 *
    TemplateFitSignalEstimator::
    Background(const TH1 * data) const {
        return 0;
    }

    TH1 *
    TemplateFitSignalEstimator::
    Signal(const TH1 * data) const {
        return 0;
    }

    void
    TemplateFitSignalEstimator::
    SaveTo(TDirectory * dir, const std::string &) const {
    }

    fit::Vector
    TemplateFitSignalEstimator::
    ToCalculatorParams(const TH1 * signal_params,
                       const std::map<std::string, TH1*> & bkgd_params) const {

        Matrix calc_params(fDims[0], bkgd_params.size()+1);
        calc_params.col(0) = fParamMap.ToMinimizerParams(root::MapContentsToEigen(signal_params));
        for(auto i = 1u; i < calc_params.cols(); i++) {
            calc_params.col(i) = fParamMap.ToMinimizerParams(root::MapContentsToEigen(bkgd_params.at(fComponentLabelIdxMap.at(i))));
        }
        return calc_params.reshaped();
    }

    void
    TemplateFitSignalEstimator::
    ToUserParams(const fit::Vector & calc_params,
                 TH1 * signal_params,
                 std::map<std::string, TH1*> & bkgd_params) const {

        Matrix calc_param_m = calc_params.reshaped(fDims[0],
                                                   fBackgroundTemplates.size()+1);

        signal_params->SetContent(fParamMap.ToUserParams(calc_param_m.col(0)).data());
        for(auto i = 1u; calc_param_m.cols(); i++) {
            bkgd_params[fComponentLabelIdxMap.at(i)]->SetContent(fParamMap.ToUserParams(calc_param_m.col(i)).data());
        }
    }

    TH1 *
    TemplateFitSignalEstimator::
    Predict(const TH1 * signal_params,
            const std::map<std::string, TH1*> & bkgd_params) const {
        auto prediction = fFitCalc->Predict(ToCalculatorParams(signal_params, bkgd_params));
        return root::ToROOT(prediction,
                            fSignalTemplate);
    }

    TH1 *
    TemplateFitSignalEstimator::
    PredictProjected(const TH1 * signal_params,
                     const std::map<std::string, TH1*> & bkgd_params) const {
        Array padded_projection = Array::Zero(fProjectPredictionProps.nbins_and_uof);
        auto prediction = fFitCalc->Predict(ToCalculatorParams(signal_params, bkgd_params));
        if(fDims[0]==1) { // single analysis bin
            padded_projection(1) = prediction.sum();
        }
        else { // two-dimensional analysis binning
            padded_projection = fParamMap.ToUserParams(prediction.reshaped(fDims[1], fDims[0])
                                                               .colwise()
                                                               .sum());
        }
        return root::ToROOT(padded_projection,
                            fProjectPredictionProps);
    }

    /*
    std::pair<TH1*, std::map<std::string, TH1*>>
    TemplateFitSignalEstimator::
    PredictComponents(const TH1 * signal_params,
                      const std::map<std::string, TH1*> & bkgd_params) const {

        int nrows = fFitCalc->GetNMinimizerParams() / fFitCalc->GetNComponents();
        int ncols = fFitCalc->GetNComponents();
        Matrix calc_param_m = Matrix::Zero(nrows, ncols);
        calc_param_m.col(0) = signal_params.reshaped();
        auto signal_prediction = root::ToROOT(fFitCalc->Predict(calc_param_m.reshaped()).array(),
                                              fPredictionProps);

        std::map<std::string, TH1*> bkgd_predictions;
        for(auto i = 1u; i < ncols; i++) {
            calc_param_m = Matrix::Zero(nrows, ncols);
            calc_param_m.col(i) = bkgd_params.at(fComponentLabelIdxMap.at(i)).reshaped();
            bkgd_predictions[fComponentLabelIdxMap.at(i)] =
                    root::ToROOT(fFitCalc->Predict(calc_param_m.reshaped()).array(),
                                 fPredictionProps);
        }
        return {signal_prediction, bkgd_predictions};

    }
     */

}
