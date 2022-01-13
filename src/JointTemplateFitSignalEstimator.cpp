#include "XSecAna/JointTemplateFitSignalEstimator.h"
#include "XSecAna/Fit/JointTemplateFitComponent.h"
#include "XSecAna/Utils.h"
namespace xsec {
    JointTemplateFitSignalEstimator::
    JointTemplateFitSignalEstimator(const std::map<std::string, fit::TemplateFitSample> & samples,
                                    const TH1 * mask)
            : fReducer(mask) {
        fJointEstimator = new TemplateFitSignalEstimator(fit::detail::_join(samples), mask);
        for(const auto & sample : samples) {
            fSampleEstimators[sample.first] = new TemplateFitSignalEstimator(sample.second, mask);
        }
    }

    void
    JointTemplateFitSignalEstimator::
    FixComponent(const std::string & template_name, const double & val) {
        fJointEstimator->FixComponent(template_name, val);
    }

    void
    JointTemplateFitSignalEstimator::
    ReleaseComponent(const std::string & template_name) {
        fJointEstimator->ReleaseComponent(template_name);
    }

    TH1 *
    JointTemplateFitSignalEstimator::
    PredictTotal(const std::map<std::string, TH1*> & params) const {
        return fJointEstimator->PredictTotal(params);
    }

    TH1 *
    JointTemplateFitSignalEstimator::
    PredictComponent(const std::string & component_label, const TH1 * params) const {
        return fJointEstimator->PredictComponent(component_label, params);
    }

    double
    JointTemplateFitSignalEstimator::
    Chi2(const TH1 * data, const std::map<std::string, TH1*> & params) const {
        return fJointEstimator->Chi2(data, params);
    }

    TH2D *
    JointTemplateFitSignalEstimator::
    GetTotalCovariance() const {
        return fJointEstimator->GetTotalCovariance();
    }

    TH2D *
    JointTemplateFitSignalEstimator::
    GetCovariance(const std::string & systematic_name) const {
        return fJointEstimator->GetCovariance(systematic_name);
    }

    TH2D *
    JointTemplateFitSignalEstimator::
    GetInverseCovariance() const {
        return fJointEstimator->GetInverseCovariance();
    }


    std::map<std::string, TemplateFitResult>
    JointTemplateFitSignalEstimator::
    Fit(const TH1 * data, int nrandom_seeds) const {
        return _joint_fit_result(fJointEstimator->Fit(data, nrandom_seeds));
    }

    std::map<std::string, TemplateFitResult>
    JointTemplateFitSignalEstimator::
    Fit(const TH1 * data, fit::IFitter * fitter, int nrandom_seeds) {
        return _joint_fit_result(fJointEstimator->Fit(data, fitter, nrandom_seeds));
    }

    std::map<std::string, TH1*>
    JointTemplateFitSignalEstimator::
    ComplimentaryParams(const std::map<std::string, TH1*> & condi_params) const {
        std::map<std::string, TH1*> comp_params;
        for(const auto & component : condi_params) {
            comp_params[component.first] = _condi_params_to_comp_params(component.first, component.second);
        }
        return comp_params;
    }

    TH1 *
    JointTemplateFitSignalEstimator::
    _condi_params_to_comp_params(const std::string & component_name, const TH1 * condi) const {
        Array _condi = fReducer.GetMap().ToMinimizerParams(root::MapContentsToEigen(condi));
        Array _comp =
                ((fit::ReducedJointTemplateComponent*) fJointEstimator->GetReducedComponent(component_name))
                        ->ComplimentaryParams(_condi);
        return root::ToROOTLike(condi,
                                fReducer.GetMap().ToUserParams(_comp));
    }

    std::map<std::string, TemplateFitResult>
    JointTemplateFitSignalEstimator::
    _joint_fit_result(const TemplateFitResult & condi_sample_fit_result) const {
        std::map<std::string, TemplateFitResult> joint_results;

        auto sample_it = fSampleEstimators.begin();
        joint_results[sample_it->first] = condi_sample_fit_result;
        std::map<std::string, TH1*> comp_params;
        std::map<std::string, TH1*> comp_params_error_up;
        std::map<std::string, TH1*> comp_params_error_down;
        for (const auto & component: condi_sample_fit_result.component_params) {
            comp_params[component.first] =
                    _condi_params_to_comp_params(component.first,
                                                 condi_sample_fit_result.component_params.at(component.first));
            comp_params_error_up[component.first] =
                    _condi_params_to_comp_params(component.first,
                                                 condi_sample_fit_result.component_params_error_up.at(component.first));
            comp_params_error_down[component.first] =
                    _condi_params_to_comp_params(component.first,
                                                 condi_sample_fit_result.component_params_error_down.at(component.first));
        }
        sample_it++;
        joint_results[sample_it->first] = {
                condi_sample_fit_result.fun_val,
                comp_params,
                comp_params_error_up,
                comp_params_error_down,
                0, // is there a way to calculate this? Not sure yet so don't return anything misleading
                condi_sample_fit_result.fun_calls,
        };
        return joint_results;
    }
}
