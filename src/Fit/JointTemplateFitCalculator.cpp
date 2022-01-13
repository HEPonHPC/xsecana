#include "XSecAna/Fit/JointTemplateFitCalculator.h"

namespace xsec {
    namespace fit {

        SampleCorrelatedComponent::
        SampleCorrelatedComponent(Matrix sample_covariance,
                                  Array sample_mean,
                                  std::vector<int> & sample_divs)
                : fJointCovariance(sample_covariance),
                  fJointMean(sample_mean),
                  fSampleDivs(sample_divs) {
            fConditioningSampleInvCovariance = fJointCovariance.block(
                    0, 0, fSampleDivs[0], fSampleDivs[0]
            ).inverse();

            fCrossTermCovariance = fJointCovariance.block(
                    fSampleDivs[0], 0, fSampleDivs[0], fSampleDivs[0]
            );

            for (auto i = 0u; i < fSampleDivs.size(); i++) {
                fSampleMeans.push_back(fJointMean(Eigen::seqN(0, fSampleDivs[i])));
            }
            fSampleMeans.push_back(fJointMean(Eigen::seqN(fSampleDivs.back(), fJointMean - fSampleDivs.back())));
        }

        Vector
        SampleCorrelatedComponent::
        Predict(const Array & component_params) const {
            auto prediction = Array::Zero(fJointMean.size());
            prediction(Eigen::seqN(0, fSampleDivs[0])) =
                    component_params * fSampleMeans[0];
            prediction(Eigen::seqN(fSampleDivs[0], fSampleMeans[1].size())) =
                    fSampleMeans[1].matrix() +
                    fCrossTermCovariance * fConditioningSampleInvCovariance *
                    (fSampleMeans[0] * (component_params - 1)).matrix();
            return prediction;
        }

        JointTemplateFitCalculator::
        JointTemplateFitCalculator(const std::vector<Array> & templates,
                                   const std::vector<ITemplateComponent*> & sample_components,
                                   const std::vector<int> & dims,
                                   const Matrix & systematic_covariance,
                                   const std::vector<int> & divs)
                : TemplateFitCalculator(templates, dims, systematic_covariance),
                  fSampleDivs(divs),
                  fComponents(sample_components, dims) {
            // Due to time constraints, only support joint fits with 2 samples
            // Extending to more dimensions should be possible later on
            assert(fSampleDivs.size() == 1 && "JointTemplateFitCalculator only supports joint fit of 2 samples");

        }

        Vector
        JointTemplateFitCalculator::
        Predict(const Vector & user_params) const {
            return fComponents.Predict(user_params);
        }

        Vector
        JointTemplateFitCalculator::
        PredictComponent(const int & component_idx, const Vector & user_params) const {
            return fComponents.PredictComponent(component_idx, user_params);
        }
    }
}
