#pragma once

#include "XSecAna/TemplateFitSignalEstimator.h"
#include "XSecAna/Fit/JointTemplateFitComponent.h"
namespace xsec {

    class JointTemplateFitSignalEstimator {
    public:
        JointTemplateFitSignalEstimator(const std::map<std::string, fit::TemplateFitSample> & samples,
                                        const TH1 * mask = 0);

        [[nodiscard]] std::map<std::string, TemplateFitResult> Fit(const std::shared_ptr<TH1> data, int nrandom_seeds=-1) const;
        [[nodiscard]] std::map<std::string, TemplateFitResult> Fit(const std::shared_ptr<TH1> data, fit::IFitter * fitter, int nrandom_seeds=-1);

        TemplateFitSignalEstimator * GetSampleTemplateFit(const std::string & name) const { return fSampleEstimators.at(name); }
        TemplateFitSignalEstimator * GetJointTemplateFit() const { return fJointEstimator; }

        std::map<std::string, TH1*> ComplimentaryParams(const std::map<std::string, TH1*> & condi_params) const;

        void FixComponent(const std::string & template_name, const double & val=1);
        void ReleaseComponent(const std::string & template_name);

        TH1 * PredictTotal(const std::map<std::string, TH1*> & params) const;
        TH1 * PredictComponent(const std::string & component_label, const TH1 * params) const;
        TH1 * NominalTotal() const
        { return fJointEstimator->NominalTotal(); }

        double Chi2(const std::shared_ptr<TH1> data, const std::map<std::string, TH1*> & params) const;

        TH2D * GetTotalCovariance() const;
        TH2D * GetCovariance(const std::string & systematic_name) const;
        TH2D * GetInverseCovariance() const;
    private:
        std::map<std::string, TemplateFitResult> _joint_fit_result(const TemplateFitResult & condi_sample_fit_result) const;
        TH1 * _condi_params_to_comp_params(const std::string & name, const TH1 * condi) const;
        TemplateFitSignalEstimator * fJointEstimator;
        std::map<std::string, TemplateFitSignalEstimator*> fSampleEstimators;

        fit::ComponentReducer fReducer;

    };

}