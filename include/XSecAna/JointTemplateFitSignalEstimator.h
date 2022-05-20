#pragma once

#include "XSecAna/TemplateFitSignalEstimator.h"
#include "XSecAna/Fit/JointTemplateFitComponent.h"
#include <set>
namespace xsec {

    class JointTemplateFitSignalEstimator {
    public:
        JointTemplateFitSignalEstimator(const std::map<std::string, fit::TemplateFitSample> & samples,
                                        const std::map<std::string, std::string> & component_conditioning,
                                        const TH1 * mask = 0);
        
        [[nodiscard]] std::map<std::string, TemplateFitResult> Fit(const std::map<std::string, std::shared_ptr<TH1>> data,
                                                                   int nrandom_seeds=-1) const;
        [[nodiscard]] std::map<std::string, TemplateFitResult> Fit(const std::map<std::string, std::shared_ptr<TH1>> data,
                                                                   fit::IFitter * fitter,
                                                                   int nrandom_seeds=-1);

        [[nodiscard]] std::map<std::string, TemplateFitResult> Fit(const std::map<std::string, std::shared_ptr<TH1>> data,
                                                                   const std::map<std::string, TH1*> & seed) const;
        [[nodiscard]] std::map<std::string, TemplateFitResult> Fit(const std::map<std::string, std::shared_ptr<TH1>> data,
                                                                   fit::IFitter * fitter,
                                                                   const std::map<std::string, TH1*> & seed);

        TemplateFitSignalEstimator * GetSampleTemplateFit(const std::string & name) const { return fSampleEstimators.at(name); }
        TemplateFitSignalEstimator * GetJointTemplateFit() const { return fJointEstimator; }
        TemplateFitSignalEstimator * GetInvertedJointTemplateFit() const { return fJointEstimatorInverse; }

        std::map<std::string, TH1*> ComplimentaryParams(const std::map<std::string, TH1*> & condi_params) const;
        std::map<std::string, TH1*> InvertParams(const std::map<std::string, TH1*> & condi_params) const;

        void FixComponent(const std::string & template_name, const double & val=1);
        void ReleaseComponent(const std::string & template_name);

        TH1 * PredictTotal(const std::map<std::string, TH1*> & params) const;
        TH1 * PredictComponent(const std::string & component_label, const TH1 * params) const;
        TH1 * NominalTotal() const
        { return fJointEstimator->NominalTotal(); }

        double Chi2(const std::shared_ptr<TH1> data, const std::map<std::string, TH1*> & params) const;

        TH2D * GetTotalCovariance(const std::map<std::string, TH1*> & params) const;
        TH2D * GetTotalSystematicCovariance() const;
        TH2D * GetSystematicCovariance(const std::string & systematic_name) const;
        TH2D * GetInverseCovariance() const;

        static std::shared_ptr<TH1> JoinData(const std::map<std::string, std::shared_ptr<TH1>> & data_samples);
    private:
        std::map<std::string, TemplateFitResult> _joint_fit_result(const TemplateFitResult & condi_sample_fit_result) const;
        std::map<std::string, TemplateFitResult> _joint_fit_result_run_inverse(TemplateFitResult condi_sample_fit_result,
                                                                               std::shared_ptr<TH1> inverted_data,
                                                                               int nrandom_seeds) const;
        std::map<std::string, TemplateFitResult> _joint_fit_result_run_inverse(TemplateFitResult condi_sample_fit_result,
                                                                               std::shared_ptr<TH1> inverted_data,
                                                                               const std::map<std::string, TH1*> & seed) const;
        void _update_results_with_bias_uncertainty(TemplateFitResult & condi, TemplateFitResult & comp) const;

        void _update_fixed_component_fit_result(TemplateFitResult & result,
                                                std::string component_label,
                                                const Systematic<TH1> & prefit_error) const;

        template<class T>
        std::map<std::string, T> _invert_samples(const std::map<std::string, T> & samples) const;
        TH1 * _condi_params_to_comp_params(const std::string & name, const TH1 * condi) const;
        TemplateFitSignalEstimator * fJointEstimator;
        TemplateFitSignalEstimator * fJointEstimatorInverse;
        std::map<std::string, TemplateFitSignalEstimator*> fSampleEstimators;

        std::set<std::string> fFixedComponentLabels;
        fit::ComponentReducer fReducer;

        std::string fCondiSampleLabel;
        std::string fComplimentarySampleLabel;

    };

}