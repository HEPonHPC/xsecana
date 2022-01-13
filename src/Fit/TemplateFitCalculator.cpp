//
// Created by Derek Doyle on 10/9/21.
//

#include "XSecAna/Fit/TemplateFitCalculator.h"

namespace xsec {
    namespace fit {
        double
        TemplateFitCalculator::
        Chi2(const Vector & user_params,
             const Vector & data) const {
            Vector u = this->Predict(user_params);
            Vector v = data - u;

            Matrix total_covariance = fSystematicCovariance;
            // add statistical uncertainty of reweighted prediction
            total_covariance += u.asDiagonal();
            return v.dot(total_covariance.llt().solve(v));
        }

        double
        TemplateFitCalculator::
        fun(const Vector & minimizer_params,
            const Vector & data) const {
            fNFunCalls++;
            return this->Chi2(this->ToUserParams(minimizer_params), data);
        }

        /// \brief User-level function for returning sum
        /// of all templates given the input template normalization
        /// parameters
        Vector
        TemplateFitCalculator::
        Predict(const Vector & user_params) const {
            return fComponents.Predict(user_params);
        }

        Vector
        TemplateFitCalculator::
        PredictComponent(const std::string & component_label, const Vector & user_params) const {
            return fComponents.PredictComponent(component_label, user_params);
        }

        void
        TemplateFitCalculator::
        ReleaseTemplate(const int & template_idx) {
            fFixedParams(template_idx) = 0;
            fParamMap.UnmaskTemplate(template_idx);
        }

        void
        TemplateFitCalculator::
        FixTemplate(const int & template_idx, const double & at) {
            fFixedParams(template_idx) = at;
            fParamMap.MaskTemplate(template_idx);
        }

        Vector
        TemplateFitCalculator::
        ToUserParams(const Vector & minimizer_coords) const {
            return fParamMap.ToUserParams(minimizer_coords) + fFixedParams;
        }

        Vector
        TemplateFitCalculator::
        ToMinimizerParams(const Vector & user_coords) const {
            return fParamMap.ToMinimizerParams(user_coords);
        }

        TemplateFitCalculator::
        TemplateFitCalculator(const ReducedComponentCollection & templates,
                              const std::vector<int> & dims,
                              const Matrix & systematic_covariance)
                : fSystematicCovariance(systematic_covariance),
                  fNFunCalls(0),
                  fComponents(templates) {


            // check we have an appropriate covariance matrix
            assert(systematic_covariance.rows() == fComponents.GetNOuterBins() * fComponents.GetNInnerBins() &&
                   systematic_covariance.cols() == fComponents.GetNOuterBins() * fComponents.GetNInnerBins());

            fNUserParams = fComponents.GetNOuterBins() * fComponents.size();
            fNComponents = fComponents.size();

            fParamMap = detail::ParamMap(fNUserParams);
            fFixedParams = Eigen::RowVectorXd::Zero(fNUserParams);
        }
    }
}


