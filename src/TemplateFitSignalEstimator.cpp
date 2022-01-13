//
// Created by Derek Doyle on 10/9/21.
//
#include "XSecAna/TemplateFitSignalEstimator.h"
#include "XSecAna/SimpleQuadSum.h"
#include "TGaxis.h"
#include "TVectorD.h"
#include "TObjString.h"
#include "TKey.h"

#include <random>
#include <iostream>
#include <iomanip>

namespace xsec {
    void
    TemplateFitResult::
    SaveTo(TDirectory * dir, const std::string & name) const {
        dir->mkdir(name.c_str());
        auto subdir = dir->GetDirectory(name.c_str());

        subdir->cd();
        TObjString("TemplateFitResult").Write("type");

        for(const auto & component : component_params) {
            subdir->mkdir(component.first.c_str());
            auto component_dir = subdir->GetDirectory(component.first.c_str());
            component_dir->cd();

            component_params.at(component.first)->Write("params");
            component_params_error_up.at(component.first)->Write("params_error_up");
            component_params_error_down.at(component.first)->Write("params_error_down");
        }
        subdir->cd();
        if(covariance) covariance->Write("param_covariance");

        TVectorD _fun_calls(1); _fun_calls[0] = fun_calls;
        TVectorD _fun_val(1); _fun_val[0] = fun_val;

        _fun_calls.Write("fun_calls");
        _fun_val.Write("fun_val");

        dir->cd();
    }

    std::unique_ptr<TemplateFitResult>
    TemplateFitResult::
    LoadFrom(TDirectory * dir, const std::string & name) {
        auto subdir = dir->GetDirectory(name.c_str());
        auto ptag = (TObjString*) subdir->Get("type");
        assert(ptag->GetString() == "TemplateFitResult" && "Type does not match TemplateFitResult");
        delete ptag;

        auto ret = std::make_unique<TemplateFitResult>();

        if(subdir->GetListOfKeys()->Contains("param_covariance")) {
            ret->covariance = (TH2D *) root::LoadTH1(subdir, "param_covariance").release();
        }
        ret->fun_calls = (*((TVectorD*)subdir->Get("fun_calls")))[0];
        ret->fun_val = (*((TVectorD*)subdir->Get("fun_val")))[0];

        for(const auto&& obj : *subdir->GetListOfKeys()) {
            if(std::string(((TKey*) obj)->GetClassName()) == "TDirectoryFile") {
                std::string component(obj->GetName());

                auto component_dir = subdir->GetDirectory(component.c_str());
                ret->component_params[component] = root::LoadTH1(component_dir, "params").release();
                ret->component_params_error_up[component] = root::LoadTH1(component_dir, "params_error_up").release();
                ret->component_params_error_down[component] = root::LoadTH1(component_dir, "params_error_down").release();
            }
        }
        return ret;
    }


    std::vector<std::string>
    TemplateFitSignalEstimator::
    GetSystematicLabels() const {
        std::vector<std::string> ret;
        for (auto syst: fSystematics) {
            ret.push_back(syst.first);
        }
        return ret;
    }

    TemplateFitSignalEstimator::
    TemplateFitSignalEstimator(const fit::TemplateFitSample & sample,
                               const TH1 * mask)
            : fReducer(mask),
              fUserComponents(sample.components) {
        fReducedComponents = fUserComponents.Reduce(fReducer);

        int icomponent = 0;
        for (const auto & component: fReducedComponents.GetComponents()) {
            fComponentIdx[component.first] = icomponent;
            icomponent++;
        }

        // given original template shapes, determine how to return predictions
        TH1 * project_predictions_like;
        auto tmp_component = fUserComponents.GetComponents().begin()->second->GetNominal();
        if (tmp_component->GetDimension() == 1) {
            project_predictions_like = new TH1D("", "",
                                                1, 0, 1);
        } else if (tmp_component->GetDimension() == 2) {
            project_predictions_like = new TH1D("", "",
                                                tmp_component->GetNbinsX(),
                                                tmp_component->GetXaxis()->GetXbins()->GetArray());
            project_predictions_like->GetXaxis()->SetTitle(tmp_component->GetXaxis()->GetTitle());
        } else if (tmp_component->GetDimension() == 3) {
            project_predictions_like = new TH2D("", "",
                                                tmp_component->GetNbinsX(),
                                                tmp_component->GetXaxis()->GetXbins()->GetArray(),
                                                tmp_component->GetNbinsY(),
                                                tmp_component->GetYaxis()->GetXbins()->GetArray());
            project_predictions_like->GetXaxis()->SetTitle(tmp_component->GetXaxis()->GetTitle());
            project_predictions_like->GetYaxis()->SetTitle(tmp_component->GetYaxis()->GetTitle());
        } else {
            throw std::runtime_error("Template fits with greater than 2 outer dimensions not supported");
        }
        // root won't draw histogram unless there are entries
        project_predictions_like->SetEntries(1);
        fProjectPredictionProps = root::TH1Props(project_predictions_like);
        fPredictionProps = root::TH1Props(tmp_component.get());

        // create reduced systematics
        for (auto syst: sample.shape_only_systematics) {
            fSystematics[syst.first] = fReducer.Reduce(syst.second);
        }

        // masked parameter map
        auto outer_map = root::MapContentsToEigen(mask);
        fOuterBinMap = fit::detail::ParamMap(outer_map);
        Array2D inner_map = Array2D::Zero(fReducedComponents.GetNInnerBins() + 2, outer_map.size());
        for (auto i = 0; i < outer_map.size(); i++) {
            if (outer_map(i))
                inner_map.col(i)(Eigen::seqN(1, fReducedComponents.GetNInnerBins())) = Array::Ones(
                        fReducedComponents.GetNInnerBins());
        }
        fInnerBinMap = fit::detail::ParamMap(inner_map.reshaped());

        // initialize all templates as free
        for (const auto & comp: fReducedComponents.GetComponents()) {
            fIsFreeTemplate[comp.first] = true;
        }

        fTotalTemplate = nullptr;
        for (const auto & component: fReducedComponents.GetComponents()) {
            if (!fTotalTemplate) {
                fTotalTemplate = (TH1 *) component.second->GetNominal()->GetHist()->Clone();
            } else {
                fTotalTemplate->Add(component.second->GetNominal()->GetHist().get());
            }
        }

        // calculate and store covariance matrices
        fTotalCovariance = 0;
        for (auto systematic: fSystematics) {
            fCovarianceMatrices[systematic.first] = systematic.second.CovarianceMatrix(fTotalTemplate);
            fCovarianceMatrices.at(systematic.first)->GetXaxis()->SetTitle("Template Bins");
            fCovarianceMatrices.at(systematic.first)->GetYaxis()->SetTitle("Template Bins");
            fCovarianceMatrices.at(systematic.first)->SetTitle((systematic.first + " Covariance").c_str());
            if (!fTotalCovariance) {
                fTotalCovariance = (TH1 *) fCovarianceMatrices.at(systematic.first)->Clone();
                fTotalCovariance->SetTitle("Total Covariance");
            } else {
                fTotalCovariance->Add(fCovarianceMatrices.at(systematic.first));
            }
        }

        Matrix total_covariance = root::MapContentsToEigenInner(fTotalCovariance)
                .matrix().reshaped(fTotalCovariance->GetNbinsX(),
                                   fTotalCovariance->GetNbinsX());
        fInvTotalCovariance = (TH1 *) fTotalCovariance->Clone();
        fInvTotalCovariance->SetTitle("Inverse Total Covariance");
        root::FillTH2Contents((TH2 *) fInvTotalCovariance, total_covariance.inverse());

        fFitCalc = new fit::TemplateFitCalculator(fReducedComponents,
                                                  {fReducedComponents.GetNOuterBins(),
                                                   fReducedComponents.GetNInnerBins()},
                                                  total_covariance);

    }

    TH2D *
    TemplateFitSignalEstimator::
    GetTotalCovariance() const {
        return (TH2D *) fTotalCovariance;
    }

    TH2D *
    TemplateFitSignalEstimator::
    GetInverseCovariance() const {
        return (TH2D *) fInvTotalCovariance;
    }

    TH2D *
    TemplateFitSignalEstimator::
    GetCovariance(const std::string & systematic_name) const {
        return (TH2D *) fCovarianceMatrices.at(systematic_name);
    }

    void
    TemplateFitSignalEstimator::
    FixComponent(const std::string & component_label, const double & val) {
        auto idx = fReducedComponents.GetComponentIdx(component_label);
        for (auto o = 0u; o < fReducedComponents.GetNOuterBins(); o++) {
            fFitCalc->FixTemplate(idx * fReducedComponents.GetNOuterBins() + o, val);
        }
        fIsFreeTemplate.at(component_label) = false;
    }

    void
    TemplateFitSignalEstimator::
    ReleaseComponent(const std::string & component_label) {
        auto idx = fReducedComponents.GetComponentIdx(component_label);
        for (auto o = 0u; o < fReducedComponents.GetNOuterBins(); o++) {
            fFitCalc->ReleaseTemplate(idx * fReducedComponents.GetNOuterBins() + o);
        }
        fIsFreeTemplate.at(component_label) = true;
    }

    void
    TemplateFitSignalEstimator::
    SaveTo(TDirectory * dir, const std::string &) const {

    }

    fit::Vector
    TemplateFitSignalEstimator::
    ToCalculatorParams(const std::map<std::string, TH1 *> & params) const {
        Matrix calc_params(fReducedComponents.GetNOuterBins(), fReducedComponents.size());
        for (auto label_idx: fComponentIdx) {
            calc_params.col(label_idx.second) = fOuterBinMap.ToMinimizerParams(
                    root::MapContentsToEigen(params.at(label_idx.first))
            );
        }
        return calc_params.reshaped();
    }

    void
    TemplateFitSignalEstimator::
    ToUserParams(const fit::Vector & calc_params,
                 std::map<std::string, TH1 *> & params) const {
        Matrix calc_param_m = calc_params.reshaped(fReducedComponents.GetNOuterBins(),
                                                   fReducedComponents.size());
        for (auto label_idx: fComponentIdx) {
            params.at(label_idx.first)->SetContent(
                    fOuterBinMap.ToUserParams(
                            calc_param_m.col(label_idx.second)
                    ).data()
            );
        }
    }

    TH1 *
    TemplateFitSignalEstimator::
    _to_template_binning(const Array & reduced_templates) const {
        Array expanded = fInnerBinMap.ToUserParams(reduced_templates);
        Array transposed = expanded
                        .reshaped(fReducedComponents.GetNInnerBins() + 2,
                                  fInnerBinMap.GetMatrix().rows() / (fReducedComponents.GetNInnerBins() + 2))
                        .transpose().reshaped();
        return root::ToROOT(
                fInnerBinMap.ToUserParams(reduced_templates)
                        .reshaped(fReducedComponents.GetNInnerBins() + 2,
                                  fInnerBinMap.GetMatrix().rows() / (fReducedComponents.GetNInnerBins() + 2))
                        .transpose().reshaped(),
                fPredictionProps
        );
    }

    TH1 *
    TemplateFitSignalEstimator::
    PredictTotal(const std::map<std::string, TH1 *> & params) const {
        return _to_template_binning(fFitCalc->Predict(ToCalculatorParams(params)));
    }

    TH1 *
    TemplateFitSignalEstimator::
    PredictProjectedTotal(const std::map<std::string, TH1 *> & params) const {
        return _project_prediction(fFitCalc->Predict(ToCalculatorParams(params)));
    }

    TH1 *
    TemplateFitSignalEstimator::
    _project_prediction(const Array & prediction) const {
        Array padded_projection = Array::Zero(fProjectPredictionProps.nbins_and_uof);
        if (fReducedComponents.GetNOuterBins() == 1) { // single analysis bin
            padded_projection(1) = prediction.sum();
        } else { // two-dimensional analysis binning
            padded_projection = fOuterBinMap.ToUserParams(prediction.reshaped(fFitCalc->GetNInnerBins(),
                                                                              fFitCalc->GetNOuterBins())
                                                                  .colwise()
                                                                  .sum().reshaped());
        }
        return root::ToROOT(padded_projection,
                            fProjectPredictionProps);
    }


    TH1 *
    TemplateFitSignalEstimator::
    NominalProjectedTotal() const {
        Array ones = Array::Ones(fFitCalc->GetNOuterBins() *
                                 fFitCalc->GetNComponents());
        return _project_prediction(fFitCalc->Predict(ones));
    }

    TH1 *
    TemplateFitSignalEstimator::
    NominalTotal() const {
        Array ones = Array::Ones(fFitCalc->GetNOuterBins() *
                                 fFitCalc->GetNComponents());
        return _to_template_binning(fFitCalc->Predict(ones));
    }


    TH1 *
    TemplateFitSignalEstimator::
    PredictComponent(const std::string & component_label, const TH1 * params) const {
        Vector mparams = fOuterBinMap.ToMinimizerParams(root::MapContentsToEigen(params));
        return _to_template_binning(fReducedComponents.PredictComponent(component_label, mparams));
    }

    TH1 *
    TemplateFitSignalEstimator::
    PredictProjectedComponent(const std::string & component_label, const TH1 * params) const {
        Vector mparams = fOuterBinMap.ToMinimizerParams(root::MapContentsToEigen(params));
        return _project_prediction(fReducedComponents.PredictComponent(component_label, mparams));
    }

    double
    TemplateFitSignalEstimator::
    Chi2(const std::shared_ptr<TH1> data,
         const std::map<std::string, TH1 *> & params) const {
        Array mparams = ToCalculatorParams(params);
        Array mdata = fReducer.Reduce(data)->GetArray();
        return fFitCalc->Chi2(ToCalculatorParams(params), fReducer.Reduce(data)->GetArray());
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

    fit::IFitCalculator *
    TemplateFitSignalEstimator::
    GetFitCalc() const {
        return fFitCalc;
    }

    TemplateFitResult
    TemplateFitSignalEstimator::
    Fit(const std::shared_ptr<TH1> data, int nrandom_seeds) const {
        if (nrandom_seeds < 0) {
            // no random seeds
            // seed at all parameters equal to 1
            return this->_fit(data, {Array::Ones(fFitCalc->GetNMinimizerParams())});
        } else {
            // throw 100 random seed configurations
            // evaluate fit calc at seed set and take the best nrandom_seeds
            // seed fit with these
            int nrandom_throws = 100;
            if (nrandom_seeds > nrandom_throws) {
                nrandom_throws = nrandom_seeds;
            }
            return this->_fit(
                    data,
                    fFitCalc->GetBestSeeds(
                            fReducer.Reduce(data)->GetArray(),
                            fFitCalc->GetRandomSeeds(nrandom_throws, -0.5, 0.5),
                            nrandom_seeds
                    )
            );
        }
    }

    TemplateFitResult
    TemplateFitSignalEstimator::
    _fit(const std::shared_ptr<TH1> data, const std::vector<Vector> & seeds) const {
        if (!fFitter) {
            throw std::runtime_error("This TemplateFitSignalEstimator does not have an active IFitter");
        }
        return this->_template_fit_result(fFitter->Fit(fFitCalc, fReducer.Reduce(data)->GetArray(), seeds));
    }

    Systematic<TH1>
    TemplateFitSignalEstimator::
    _prefit_component_uncertainty(const std::string & component_label) const {
        const auto component = fUserComponents.GetComponents().at(component_label);
        std::map<std::string, Systematic<TH1>> projected_systematics;
        for(const auto & syst : component->GetSystematics()) {
            projected_systematics[syst.first] = fit::ComponentReducer::Project(syst.second);
        }
        auto nominal = fit::ComponentReducer::Project(component->GetNominal());
        auto up = std::get<1>(SimpleQuadSum::TotalFractionalUncertainty(nominal.get(), projected_systematics)).Up();
        auto down = std::shared_ptr<TH1>((TH1*) up->Clone());
        down->Scale(-1);
        return Systematic<TH1>("", up, down);
    }

    TemplateFitResult
    TemplateFitSignalEstimator::
    Fit(const std::shared_ptr<TH1> data, fit::IFitter * fitter, int nrandom_seeds) {
        SetFitter(fitter);
        return Fit(data, nrandom_seeds);
    }

    TemplateFitResult
    TemplateFitSignalEstimator::
    _template_fit_result(const fit::FitResult & result) const {
        std::map<std::string, TH1 *> component_params;
        std::map<std::string, TH1 *> component_params_error_up;
        std::map<std::string, TH1 *> component_params_error_down;
        for (auto component: fReducedComponents.GetComponents()) {
            component_params[component.first] = (TH1 *) fReducer.GetMask()->Clone();
            component_params_error_up[component.first] = (TH1 *) fReducer.GetMask()->Clone();
            component_params_error_down[component.first] = (TH1 *) fReducer.GetMask()->Clone();
            
            component_params.at(component.first)->Reset("ICESM");
            component_params_error_up.at(component.first)->Reset("ICESM");
            component_params_error_down.at(component.first)->Reset("ICESM");

            root::CopyAxisLabels(fProjectPredictionProps, component_params.at(component.first));
            root::CopyAxisLabels(fProjectPredictionProps, component_params_error_up.at(component.first));
            root::CopyAxisLabels(fProjectPredictionProps, component_params_error_down.at(component.first));
            component_params.at(component.first)->SetTitle(("Normalization Parameters: " + component.first).c_str());
            component_params_error_up.at(component.first)->SetTitle(
                    ("Parameters Error Up: " + component.first).c_str());
            component_params_error_down.at(component.first)->SetTitle(
                    ("Parameters Error Down: " + component.first).c_str());
        }
        ToUserParams(result.params, component_params);
        ToUserParams(result.params_error_up, component_params_error_up);
        ToUserParams(result.params_error_down, component_params_error_down);
        for (auto component: fReducedComponents.GetComponents()) {
            if (!fIsFreeTemplate.at(component.first)) {
                Systematic<TH1> prefit_uncertainty = _prefit_component_uncertainty(component.first);
                component_params_error_up[component.first] = (TH1*) prefit_uncertainty.Up()->Clone();
                component_params_error_down[component.first] = (TH1*) prefit_uncertainty.Down()->Clone();
            }
        }
        auto covariance = new TH2D("", "",
                                   result.covariance.rows(), 0, result.covariance.rows(),
                                   result.covariance.rows(), 0, result.covariance.rows());
        root::FillTH2Contents(covariance, result.covariance);

        return {result.fun_val,
                component_params,
                component_params_error_up,
                component_params_error_down,
                covariance,
                result.fun_calls};
    }


    TCanvas *
    TemplateFitSignalEstimator::
    DrawParameterCorrelation(const TemplateFitResult & fit_result) const {
        auto c = new TCanvas("parameter_correlation");
        auto tmp_cor = (TH2D *) fit_result.covariance->Clone();
        for (auto i = 1; i <= tmp_cor->GetNbinsX(); i++) {
            for (auto j = 1; j <= tmp_cor->GetNbinsY(); j++) {
                double cov = fit_result.covariance->GetBinContent(i, j);
                double dii = fit_result.covariance->GetBinContent(i, i);
                double djj = fit_result.covariance->GetBinContent(j, j);
                tmp_cor->SetBinContent(i, j, cov / std::sqrt(dii * djj));
            }
        }

        this->_draw_covariance_helper(c, tmp_cor, fit_result);
        return c;
    }

    TCanvas *
    TemplateFitSignalEstimator::
    DrawParameterCovariance(const TemplateFitResult & fit_result) const {
        auto c = new TCanvas("parameter_covariance");
        auto tmp_cov = (TH2D *) fit_result.covariance->Clone();
        this->_draw_covariance_helper(c, tmp_cov, fit_result);
        return c;
    }

    void
    TemplateFitSignalEstimator::
    _draw_covariance_helper(TCanvas * c, TH1 * mat, const TemplateFitResult & fit_result) const {
        auto xmax = fit_result.covariance->GetXaxis()->GetBinLowEdge(fit_result.covariance->GetNbinsX() + 1);
        auto xaxis = new TGaxis(0, 0, xmax, 0, 0.001, fReducedComponents.size());
        auto yaxis = new TGaxis(0, 0, 0, xmax, 0.001, fReducedComponents.size());
        xaxis->SetNdivisions(fReducedComponents.size());
        yaxis->SetNdivisions(fReducedComponents.size());

        for (auto component: fComponentIdx) {
            xaxis->ChangeLabel(component.second, 40, -1, 33, -1, -1, component.first.c_str());
            yaxis->ChangeLabel(component.second, 40, -1, 33, -1, -1, component.first.c_str());
        }
        mat->GetXaxis()->SetNdivisions(fReducedComponents.size(), false);
        mat->GetYaxis()->SetNdivisions(fReducedComponents.size(), false);
        double zlim = std::max(std::abs(mat->GetMaximum()), std::abs(mat->GetMinimum()));
        mat->GetZaxis()->SetRangeUser(-zlim, zlim);
        mat->GetXaxis()->SetTickLength(0);
        mat->GetYaxis()->SetTickLength(0);
        mat->GetXaxis()->SetLabelSize(0);
        mat->GetYaxis()->SetLabelSize(0);
        mat->Draw("colz");
        xaxis->Draw("same");
        yaxis->Draw("same");
    }
}
