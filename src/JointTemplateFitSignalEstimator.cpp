#include "XSecAna/JointTemplateFitSignalEstimator.h"
#include "XSecAna/Fit/JointTemplateFitComponent.h"
#include "XSecAna/Fit/Minuit2TemplateFitter.h"
#include "XSecAna/Utils.h"
#include <algorithm>

namespace xsec {
    template<class T>
    std::map<std::string, T>
    JointTemplateFitSignalEstimator::
    _invert_samples(const std::map<std::string, T> & samples) const {
        std::map<std::string, T> inverted;
        int i = fSampleEstimators.size()-1;
        for(auto const & sample : samples) {
            inverted[std::to_string(i) + sample.first] = sample.second;
            i--;
        }
        return inverted;
    }

    JointTemplateFitSignalEstimator::
    JointTemplateFitSignalEstimator(const std::map<std::string, fit::TemplateFitSample> & samples,
                                    const TH1 * mask)
            : fReducer(mask) {
        fJointEstimator = new TemplateFitSignalEstimator(fit::detail::_join(samples), mask);
        for(const auto & sample : samples) {
            fSampleEstimators[sample.first] = new TemplateFitSignalEstimator(sample.second, mask);
        }
        std::map<std::string, fit::TemplateFitSample> inverted_samples = _invert_samples(samples);
        fJointEstimatorInverse = new TemplateFitSignalEstimator(fit::detail::_join(inverted_samples), mask);

        auto sample_it = fSampleEstimators.begin();
        fCondiSampleLabel = sample_it->first;
        sample_it++;

        fComplimentarySampleLabel = sample_it->first;
    }

    void
    JointTemplateFitSignalEstimator::
    FixComponent(const std::string & template_name, const double & val) {
        fFixedComponentLabels.insert(template_name);
        for(auto & sample : fSampleEstimators) {
            sample.second->FixComponent(template_name, val);
        }
        fJointEstimator->FixComponent(template_name, val);
        fJointEstimatorInverse->FixComponent(template_name, val);
    }

    void
    JointTemplateFitSignalEstimator::
    ReleaseComponent(const std::string & template_name) {
        fFixedComponentLabels.erase(template_name);
        for(auto & sample : fSampleEstimators) {
            sample.second->ReleaseComponent(template_name);
        }
        fJointEstimator->ReleaseComponent(template_name);
        fJointEstimatorInverse->ReleaseComponent(template_name);
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
    Chi2(const std::shared_ptr<TH1> data, const std::map<std::string, TH1*> & params) const {
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


    //std::map<std::string, TemplateFitResult>
    //JointTemplateFitSignalEstimator::
    //Fit(const std::shared_ptr<TH1> data, int nrandom_seeds) const {
    //    //return _joint_fit_result(fJointEstimator->Fit(data, nrandom_seeds));
    //    return _joint_fit_result_run_inverse(fJointEstimator->Fit(data), data);
    //}
    //
    //std::map<std::string, TemplateFitResult>
    //JointTemplateFitSignalEstimator::
    //Fit(const std::shared_ptr<TH1> data, fit::IFitter * fitter, int nrandom_seeds) {
    //    //return _joint_fit_result(fJointEstimator->Fit(data, fitter, nrandom_seeds));
    //    //return _joint_fit_result_run_inverse(fJointEstimator->Fit(data, fitter), data);
    //}

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

    std::map<std::string, TH1*>
    JointTemplateFitSignalEstimator::
    InvertParams(const std::map<std::string, TH1*> & condi_params) const {
        std::map<std::string, TH1 *> ret;
        for (const auto & component: condi_params) {
            Array inverted_component_param =
                    ((fit::ReducedJointTemplateComponent *)
                            fJointEstimatorInverse->GetReducedComponent(component.first))->ConditionalParams(
                            fJointEstimatorInverse->ToCalculatorParamsComponent(component.second)
                    );
            ret[component.first] = fJointEstimatorInverse->ToUserParamsComponent(inverted_component_param);
        }
        return ret;
    }

    void
    JointTemplateFitSignalEstimator::
    _update_results_with_bias_uncertainty(TemplateFitResult & condi, TemplateFitResult & comp) const {

        for(auto & component : condi.component_params) {
            Array std_condi = fJointEstimator->ToCalculatorParamsComponent(component.second);
            Array inv_condi =
                    ((fit::ReducedJointTemplateComponent *)
                            fJointEstimatorInverse->GetReducedComponent(component.first))->ComplimentaryParams(
                            fJointEstimatorInverse->ToCalculatorParamsComponent(
                                    comp.component_params.at(component.first))
                    );

            Array bias_error = (std_condi - inv_condi).pow(2);

            Array fit_error_up = fJointEstimator->ToCalculatorParamsComponent(
                    condi.component_params_error_up.at(component.first)
            );
            fit_error_up = fit_error_up.pow(2);

            Array fit_error_down = fJointEstimator->ToCalculatorParamsComponent(
                    condi.component_params_error_down.at(component.first)
            );
            fit_error_down = fit_error_down.pow(2);


            Array total_error = ((fit_error_down > fit_error_up).select(fit_error_down, fit_error_up) +
                    bias_error).sqrt();

            condi.component_params_error_up.at(component.first) = fJointEstimator->ToUserParamsComponent(total_error);
            condi.component_params_error_down.at(component.first) = fJointEstimator->ToUserParamsComponent(-1 * total_error);
        }
    }

    std::map<std::string, TemplateFitResult>
    JointTemplateFitSignalEstimator::
    _joint_fit_result_run_inverse(TemplateFitResult condi_sample_fit_result,
                                  const std::shared_ptr<TH1> inverted_data,
                                  int nrandom_seeds) const {
        std::cout << "Standard ordering fit done. Moving on to inverse ordering.\n";
        std::map<std::string, TemplateFitResult> joint_results;

        std::map<std::string, TH1*> inverted_seed = this->InvertParams(condi_sample_fit_result.component_params);

        ((fit::Minuit2TemplateFitter*) fJointEstimator->GetFitter())->SetMinosErrors(false);
        auto comp_sample_fit_result = fJointEstimatorInverse->Fit(inverted_data,
                                                                  fJointEstimator->GetFitter(),
                                                                  nrandom_seeds);

        auto condi_sample_fit_result_with_bias_error = condi_sample_fit_result.Clone();
        this->_update_results_with_bias_uncertainty(condi_sample_fit_result_with_bias_error, comp_sample_fit_result);
        auto sample_it = fSampleEstimators.begin();
        joint_results[sample_it->first] = condi_sample_fit_result;
        joint_results[sample_it->first + "_with_bias_error"] = condi_sample_fit_result_with_bias_error;

        sample_it++;
        joint_results[sample_it->first] = comp_sample_fit_result;

        double chi2_inv = joint_results.at(sample_it->first).fun_val;
        double chi2_std = condi_sample_fit_result.fun_val;
        double chi2_diff = std::abs(chi2_inv - chi2_std);
        double chi2_inv_at_std = fJointEstimatorInverse->Chi2(inverted_data,
                                                              inverted_seed);
        std::cout << "Inverted chi2 = " << chi2_inv << std::endl;
        std::cout << "Standard chi2 = " << chi2_std << std::endl;
        std::cout << "Difference = " << chi2_diff << std::endl;
        std::cout << "Inverted chi2 at standard best fit parameters = " << chi2_inv_at_std << std::endl;
        return joint_results;
    }


    std::map<std::string, TemplateFitResult>
    JointTemplateFitSignalEstimator::
    _joint_fit_result(const TemplateFitResult & condi_sample_fit_result) const {
        std::map<std::string, TemplateFitResult> joint_results;


        auto sample_it = fSampleEstimators.begin();
        joint_results[sample_it->first] = condi_sample_fit_result;
        // Normalization uncertainty for fixed components need to be re-calculated
        // since the joint template fitter projects all samples into the phase space,
        // but we calculate this separately for each sample
        for(const auto & fixed_component_label : fFixedComponentLabels) {
            auto prefit_uncertainty = fSampleEstimators.at(sample_it->first)->PrefitComponentUncertainty(fixed_component_label);
            joint_results[sample_it->first].component_params_error_up.at(fixed_component_label) =
                    (TH1*) prefit_uncertainty.Up()->Clone();
            joint_results[sample_it->first].component_params_error_down.at(fixed_component_label) =
                    (TH1*) prefit_uncertainty.Down()->Clone();
            joint_results[sample_it->first].component_params_error_up.at(fixed_component_label)->Multiply(
                    joint_results[sample_it->first].component_params.at(fixed_component_label)
            );
            joint_results[sample_it->first].component_params_error_down.at(fixed_component_label)->Multiply(
                    joint_results[sample_it->first].component_params.at(fixed_component_label)
            );
        }

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

    std::shared_ptr<TH1>
    JointTemplateFitSignalEstimator::
    JoinData(const std::map<std::string, std::shared_ptr<TH1>> & data_samples) {
        std::vector<std::shared_ptr<TH1>> samples;
        for(const auto & sample : data_samples) {
            samples.push_back(sample.second);
        }
        return fit::detail::_join(samples);
    }

    std::map<std::string, TemplateFitResult>
    JointTemplateFitSignalEstimator::
    Fit(const std::map<std::string, std::shared_ptr<TH1>> data,
        int nrandom_seeds) const {
        //return _joint_fit_result(fJointEstimator->Fit(data, nrandom_seeds));
        return _joint_fit_result_run_inverse(fJointEstimator->Fit(JoinData(data), nrandom_seeds),
                                             JoinData(_invert_samples(data)),
                                             nrandom_seeds);
        //return this->Fit(JoinData(data), nrandom_seeds);
    }

    std::map<std::string, TemplateFitResult>
    JointTemplateFitSignalEstimator::
    Fit(const std::map<std::string, std::shared_ptr<TH1>> data,
        fit::IFitter * fitter,
        int nrandom_seeds) {
        return _joint_fit_result_run_inverse(fJointEstimator->Fit(JoinData(data), fitter, nrandom_seeds),
                                             JoinData(_invert_samples(data)),
                                             nrandom_seeds);
        //return this->Fit(JoinData(data), fitter, nrandom_seeds);
    }

    std::map<std::string, TemplateFitResult>
    JointTemplateFitSignalEstimator::
    Fit(const std::map<std::string, std::shared_ptr<TH1>> data,
        const std::map<std::string, TH1*> & seed) const {
        return _joint_fit_result_run_inverse(fJointEstimator->Fit(JoinData(data), seed),
                                             JoinData(_invert_samples(data)),
                                             seed);

    }

    std::map<std::string, TemplateFitResult>
    JointTemplateFitSignalEstimator::
    Fit(const std::map<std::string, std::shared_ptr<TH1>> data,
        fit::IFitter * fitter,
        const std::map<std::string, TH1*> & seed) {
        return _joint_fit_result_run_inverse(fJointEstimator->Fit(JoinData(data), fitter, seed),
                                             JoinData(_invert_samples(data)),
                                             seed);
    }
}
