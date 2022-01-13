#pragma once

#include "XSecAna/Fit/TemplateFitCalculator.h"

namespace xsec {
    namespace fit {

        class SampleCorrelatedComponent : public ITemplateComponent {
        public:
            SampleCorrelatedComponent(Matrix sample_covariance,
                                      Array sample_mean,
                                      std::vector<int> & sample_divs);
            Vector Predict(const Array & component_params) const override;
        private:
            Matrix fJointCovariance;
            Array fJointMean;
            std::vector<int> & fSampleDivs;
            std::vector<Array> fSampleMeans;

            Matrix fConditioningSampleInvCovariance;
            Matrix fCrossTermCovariance;

        };



        class JointTemplateFitCalculator : public TemplateFitCalculator {
        public:
            JointTemplateFitCalculator(const std::vector<Array> & templates,
                                       const std::vector<ITemplateComponent*> & sample_correlations,
                                       const std::vector<int> & dims,
                                       const Matrix & systematic_covariance,
                                       const std::vector<int> & divs);

            Vector Predict(const Vector & user_params) const override;
            Vector PredictComponent(const int & component_idx, const Vector & user_params) const override;

        private:
            const std::vector<int> & fSampleDivs;
            ComponentCollection fComponents;
        };
    }
}