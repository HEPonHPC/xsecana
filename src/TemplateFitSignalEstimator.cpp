//
// Created by Derek Doyle on 10/9/21.
//
#include "XSecAna/TemplateFitSignalEstimator.h"

namespace xsec {
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
            // We do this instead of converting to Eigen because of limitations imposed by
            // the current implementation of Systematic
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
            project_predictions_like->GetXaxis()->SetTitle(signal_template->GetXaxis()->GetTitle());
        } else if (signal_template->GetDimension() == 3) {
            project_predictions_like = new TH2D("", "",
                                                signal_template->GetNbinsX(),
                                                signal_template->GetXaxis()->GetXbins()->GetArray(),
                                                signal_template->GetNbinsY(),
                                                signal_template->GetYaxis()->GetXbins()->GetArray());
            project_predictions_like->GetXaxis()->SetTitle(signal_template->GetXaxis()->GetTitle());
            project_predictions_like->GetYaxis()->SetTitle(signal_template->GetYaxis()->GetTitle());
        }
        else {
            throw std::runtime_error("Template fits with greater than 2 outer dimensions not supported");
        }
        // root won't draw histogram unless there are entries
        project_predictions_like->SetEntries(1);
        fProjectPredictionProps = root::TH1Props(project_predictions_like);
        fPredictionProps = root::TH1Props(signal_template);

        // construct templates with masked bins removed
        fSignalTemplate = _mask_and_flatten(mask, signal_template);
        // root stores as row-major. The calculator assumes column major, so we transpose dimensions
        // outer bins x inner bins
        fDims = {fSignalTemplate->GetNbinsY()+2, fSignalTemplate->GetNbinsX()+2};
        for(auto bkgd_template : background_templates) {
            fBackgroundTemplates[bkgd_template.first] = _mask_and_flatten(mask,
                                                                          bkgd_template.second);
        }
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

        // masked parameter map
        auto outer_map = root::MapContentsToEigen(mask);
        fOuterBinMap = fit::detail::ParamMap(outer_map);
        Array2D inner_map = Array2D::Zero(fDims[1]+2, outer_map.size());
        for(auto i = 0; i < outer_map.size(); i++) {
            if(outer_map(i)) inner_map.col(i)(Eigen::seqN(1, fDims[1])) = Array::Ones(fDims[1]);
        }
        fInnerBinMap = fit::detail::ParamMap(inner_map.reshaped());


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
            fComponentLabelIdxMap[background_template.first] = itemplate;
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

    TH1D *
    TemplateFitSignalEstimator::
    GetReducedSignalTemplate() const {
        auto n = (fSignalTemplate->GetNbinsX()+2) *
                 (fSignalTemplate->GetNbinsY()+2);
        auto ret = new TH1D("", "",
                            n, 0, n);
        Array arr(n+2);
        arr(Eigen::seqN(1, n)) = root::MapContentsToEigen(fSignalTemplate);
        ret->SetContent(arr.data());
        ret->SetEntries(1);
        ret->GetXaxis()->SetTitle("Template Bins");
        ret->GetYaxis()->SetTitle("Events");
        return ret;
    }

    TH1D *
    TemplateFitSignalEstimator::
    GetReducedTotalTemplate() const {
        auto n = (fTotalTemplate->GetNbinsX()+2) *
                 (fTotalTemplate->GetNbinsY()+2);
        auto ret = new TH1D("", "",
                            n, 0, n);
        Array arr(n+2);
        arr(Eigen::seqN(1, n)) = root::MapContentsToEigen(fTotalTemplate);
        ret->SetContent(arr.data());
        ret->SetEntries(1);
        ret->GetXaxis()->SetTitle("Template Bins");
        ret->GetYaxis()->SetTitle("Events");
        return ret;
    }

    TH1D *
    TemplateFitSignalEstimator::
    GetReducedBackgroundTemplate(const std::string & bkgd_label) const {
        auto n = (fBackgroundTemplates.at(bkgd_label)->GetNbinsX()+2) *
                 (fBackgroundTemplates.at(bkgd_label)->GetNbinsY()+2);
        auto ret = new TH1D("", "",
                            n, 0, n);
        Array arr(n+2);
        arr(Eigen::seqN(1, n)) =
                root::MapContentsToEigen(fBackgroundTemplates.at(bkgd_label));
        ret->SetContent(arr.data());
        ret->SetEntries(1);
        ret->GetXaxis()->SetTitle("Template Bins");
        ret->GetYaxis()->SetTitle("Events");
        return ret;
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
    FixComponent(const std::string & component_label, const double & val) {
        //int idx = -1;
        //for(auto comp : fComponentLabelIdxMap) {
        //    if(comp.second == template_name) {
        //        idx = comp.first;
        //    }
        //}
        //if(idx < 0) {
        //    throw std::runtime_error(("No component named " + template_name).c_str());
        //}
        // fix this component in each outer bin
        auto idx = fComponentLabelIdxMap.at(component_label);
        for(auto o = 0u; o < fDims[0]; o++) {
            fFitCalc->FixTemplate(idx*fDims[0] + o, val);
            //fFitCalc->FixTemplate(idx * fDims[0] + o, val);
        }
        fIsFreeTemplate.at(component_label) = false;
    }

    void
    TemplateFitSignalEstimator::
    ReleaseComponent(const std::string & component_label) {
        //int idx = -1;
        //for(auto comp : fComponentLabelIdxMap) {
        //    if(comp.second == template_name) {
        //        idx = comp.first;
        //    }
        //}
        //if(idx < 0) {
        //    throw std::runtime_error(("No component named " + template_name).c_str());
        //}
        auto idx = fComponentLabelIdxMap.at(component_label);
        for(auto o = 0u; o < fDims[0]; o++) {
            fFitCalc->ReleaseTemplate(idx*fDims[0] + o);
            //fFitCalc->ReleaseTemplate(idx * fDims[0] + o);
        }
        fIsFreeTemplate.at(component_label) = true;
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
        calc_params.col(0) = fOuterBinMap.ToMinimizerParams(root::MapContentsToEigen(signal_params));
        for(auto label_idx : fComponentLabelIdxMap) {
            calc_params.col(label_idx.second) = fOuterBinMap.ToMinimizerParams(
                    root::MapContentsToEigen(bkgd_params.at(label_idx.first))
            );
        }
        //for(auto i = 1u; i < calc_params.cols(); i++) {
        //    calc_params.col(i) = fOuterBinMap.ToMinimizerParams(
        //            root::MapContentsToEigen(bkgd_params.at(fComponentLabelIdxMap.at(i)))
        //    );
        //}
        return calc_params.reshaped();
    }

    void
    TemplateFitSignalEstimator::
    ToUserParams(const fit::Vector & calc_params,
                 TH1 * signal_params,
                 std::map<std::string, TH1*> & bkgd_params) const {

        Matrix calc_param_m = calc_params.reshaped(fDims[0],
                                                   fBackgroundTemplates.size()+1);

        signal_params->SetContent(fOuterBinMap.ToUserParams(calc_param_m.col(0)).data());
        for(auto label_idx : fComponentLabelIdxMap) {
            bkgd_params.at(label_idx.first)->SetContent(
                    fOuterBinMap.ToUserParams(
                            calc_param_m.col(label_idx.second)
                    ).data()
            );
        }
        //for(auto i = 1u; i < calc_param_m.cols(); i++) {
        //    bkgd_params[fComponentLabelIdxMap.at(i)]->SetContent(fOuterBinMap.ToUserParams(calc_param_m.col(i)).data());
        //}
    }

    TH1 *
    TemplateFitSignalEstimator::
    _to_template_binning(const Array & reduced_templates) const {
        return root::ToROOT(
                fInnerBinMap.ToUserParams(reduced_templates)
                        .reshaped(fDims[1] + 2, fInnerBinMap.GetMatrix().rows() / (fDims[1] + 2))
                        .transpose().reshaped(),
                fPredictionProps
        );
    }

    Array
    TemplateFitSignalEstimator::
    _predict_component(const TH1 * component_templates, const TH1 * params) const {
        return (root::MapContentsToEigen(component_templates)
                         .reshaped(fDims[1], fDims[0]) // TODO row-major
                         .matrix() *
               fOuterBinMap.ToMinimizerParams(root::MapContentsToEigen(params))
                       .asDiagonal()).reshaped();
    }

    TH1 *
    TemplateFitSignalEstimator::
    PredictTotal(const TH1 * signal_params,
                 const std::map<std::string, TH1*> & bkgd_params) const {
        return _to_template_binning(fFitCalc->Predict(ToCalculatorParams(signal_params, bkgd_params)));
    }


    TH1 *
    TemplateFitSignalEstimator::
    PredictSignal(const TH1 * signal_params) const {
        return _to_template_binning(
                fFitCalc->PredictComponent(
                        0,
                        fOuterBinMap.ToMinimizerParams(root::MapContentsToEigen(signal_params))
                )
        );
    }

    TH1 *
    TemplateFitSignalEstimator::
    PredictBackground(const std::string & background_label,
                      const TH1 * bkgd_params) const {
        return _to_template_binning(
                fFitCalc->PredictComponent(
                        fComponentLabelIdxMap.at(background_label),
                        fOuterBinMap.ToMinimizerParams(root::MapContentsToEigen(bkgd_params))
                )
        );
    }

    TH1 *
    TemplateFitSignalEstimator::
    PredictProjectedTotal(const TH1 * signal_params,
                          const std::map<std::string, TH1*> & bkgd_params) const {
        Array padded_projection = Array::Zero(fProjectPredictionProps.nbins_and_uof);
        auto prediction = fFitCalc->Predict(ToCalculatorParams(signal_params, bkgd_params));
        if(fDims[0]==1) { // single analysis bin
            padded_projection(1) = prediction.sum();
        }
        else { // two-dimensional analysis binning
            padded_projection = fOuterBinMap.ToUserParams(prediction.reshaped(fDims[1], fDims[0])
                                                               .colwise()
                                                               .sum());
        }
        return root::ToROOT(padded_projection,
                            fProjectPredictionProps);
    }

    TH1 *
    TemplateFitSignalEstimator::
    PredictProjectedSignal(const TH1 * signal_params) const {
        Array padded_projection = Array::Zero(fProjectPredictionProps.nbins_and_uof);
        auto prediction = fFitCalc->PredictComponent(
                0, fOuterBinMap.ToMinimizerParams(root::MapContentsToEigen(signal_params))
        );
        if(fDims[0]==1) {
            padded_projection(0) = prediction.sum();
        }
        else {
            padded_projection = fOuterBinMap.ToUserParams(prediction.reshaped(fDims[1], fDims[0])
                                                                  .colwise()
                                                                  .sum());
        }
        return root::ToROOT(padded_projection,
                            fProjectPredictionProps);
    }

    TH1 *
    TemplateFitSignalEstimator::
    PredictProjectedBackground(const std::string & background_label, const TH1 * bkgd_params) const {
        Array padded_projection = Array::Zero(fProjectPredictionProps.nbins_and_uof);
        auto prediction = fFitCalc->PredictComponent(
                fComponentLabelIdxMap.at(background_label),
                fOuterBinMap.ToMinimizerParams(root::MapContentsToEigen(bkgd_params))
        );
        if(fDims[0]==1) {
            padded_projection(0) = prediction.sum();
        }
        else {
            padded_projection = fOuterBinMap.ToUserParams(prediction.reshaped(fDims[1], fDims[0])
                                                                  .colwise()
                                                                  .sum());
        }
        return root::ToROOT(padded_projection,
                            fProjectPredictionProps);
    }

    void
    TemplateFitSignalEstimator::
    SetFitter(fit::IFitter * fitter) {
        fFitter = fitter;
    }

    fit::IFitter *
    TemplateFitSignalEstimator::
    GetFitter() const {
        return fFitter;
    }

    TemplateFitResult
    TemplateFitSignalEstimator::
    Fit(const TH1 * data) const {
        if(!fFitter) {
            throw std::runtime_error("This TemplateFitSignalEstimator does not have an active IFitter");
        }
        Array data_arr = root::MapContentsToEigen(_mask_and_flatten(fMask, data));
        auto result = fFitter->Fit(fFitCalc, data_arr);

        auto signal_params = (TH1*) fMask->Clone();
        auto signal_params_error_up = (TH1*) fMask->Clone();
        auto signal_params_error_down = (TH1*) fMask->Clone();
        root::CopyAxisLabels(fProjectPredictionProps, signal_params);
        root::CopyAxisLabels(fProjectPredictionProps, signal_params_error_up);
        root::CopyAxisLabels(fProjectPredictionProps, signal_params_error_down);

        signal_params->SetTitle("Normalization Parameters: Signal");
        signal_params_error_up->SetTitle("Parameters Error Up: Signal");
        signal_params_error_down->SetTitle("Parameters Error Down: Signal");

        std::map<std::string, TH1*> background_params;
        std::map<std::string, TH1*> background_params_error_up;
        std::map<std::string, TH1*> background_params_error_down;
        for(auto bkgd : fBackgroundTemplates) {
            background_params[bkgd.first] = (TH1*) fMask->Clone();
            background_params_error_up[bkgd.first] = (TH1*) fMask->Clone();
            background_params_error_down[bkgd.first] = (TH1*) fMask->Clone();

            root::CopyAxisLabels(fProjectPredictionProps, background_params.at(bkgd.first));
            root::CopyAxisLabels(fProjectPredictionProps, background_params_error_up.at(bkgd.first));
            root::CopyAxisLabels(fProjectPredictionProps, background_params_error_down.at(bkgd.first));
            background_params.at(bkgd.first)->SetTitle(("Normalization Parameters: " + bkgd.first).c_str());
            background_params_error_up.at(bkgd.first)->SetTitle(("Parameters Error Up: " + bkgd.first).c_str());
            background_params_error_down.at(bkgd.first)->SetTitle(("Parameters Error Down: " + bkgd.first).c_str());
        }
        ToUserParams(result.params, signal_params, background_params);
        ToUserParams(result.params_error_up, signal_params_error_up, background_params_error_up);
        ToUserParams(result.params_error_down, signal_params_error_down, background_params_error_down);

        return {result.fun_val,
                signal_params,
                signal_params_error_up,
                signal_params_error_down,
                background_params,
                background_params_error_up,
                background_params_error_down,
                result.covariance,
                result.fun_calls};
    }

    TemplateFitResult
    TemplateFitSignalEstimator::
    Fit(const TH1 * data, fit::IFitter * fitter) {
        SetFitter(fitter);
        return Fit(data);
    }
}
