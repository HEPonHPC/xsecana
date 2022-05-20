//
// Created by Derek Doyle on 10/9/21.
//

#include "XSecAna/Fit/TemplateFitCalculator.h"
#include <cmath>

namespace xsec {
    namespace fit {

        /*
        double
        TemplateFitCalculator::
        Chi2(const Vector & user_params,
             const Vector & data) const {
            Vector u = this->Predict(user_params);
            //Vector lnu = u.array().log();

            if(!fIsBCSetup) {
                fBCData = data.array().log();
                fBCSMatrix = (1 / data.array()).matrix().asDiagonal();
                fBCLogSMatrixDeterminant = -1 * data.array().log().sum();
                fIsBCSetup = true;
            }

            Vector v = fBCData - u.array().log().matrix();

            Matrix total_covariance = fSystematicCovariance;
            // add statistical uncertainty of reweighted prediction
            total_covariance += u.asDiagonal();

            Matrix W = fBCSMatrix * total_covariance * fBCSMatrix;

            fDecomp.compute(W);
            return v.dot(fDecomp.solve(v));// + LogDetV(fDecomp);
        }
*/

        double
        TemplateFitCalculator::
        Chi2(const Vector & user_params,
             const Vector & data) const {
            Vector u = this->Predict(user_params);
            Vector v = data - u;

            Matrix total_covariance = fSystematicCovariance;
            // add statistical uncertainty of reweighted prediction
            if(!fIgnoreStatisticalUncertainty) {
                total_covariance += u.asDiagonal();
            }
            fDecomp.compute(total_covariance);
            return v.dot(fDecomp.solve(v));// + LogDetV(fDecomp);
        }

/*
        double
        TemplateFitCalculator::
        Chi2(const Vector & user_params,
             const Vector & data) const {
            Vector u = this->Predict(user_params);
            Vector v = (data.array() / u.array() - 1);

            if(!fIsNCVSetup) {
                fNCVNominal = this->Predict(Vector::Ones(user_params.size()));
                fNCVSystematicCovariance = fSystematicCovariance;
                for(int i = 0u; i < fSystematicCovariance.rows(); i++) {
                    for(int j = 0u; j < fSystematicCovariance.cols(); j++) {
                        fNCVSystematicCovariance(i,j) = fSystematicCovariance(i,j) /
                                                        (fNCVNominal(i) * fNCVNominal(j));
                    }
                }
                fIsNCVSetup = true;
            }
            Matrix total_covariance = fNCVSystematicCovariance;

            // add statistical uncertainty of reweighted prediction
            if(!fIgnoreStatisticalUncertainty) {
                total_covariance += (1 / data.array()).matrix().asDiagonal();
            }
            fDecomp.compute(total_covariance);

            return v.dot(fDecomp.solve(v));// + LogDetV(fDecomp);
        }
*/
        double
        TemplateFitCalculator::
        LogDetV(const Eigen::LLT<Matrix> & decomp) {
            Matrix L = decomp.matrixL();
            return 2 * L.diagonal().array().log().sum();
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
                              const Matrix & systematic_covariance,
                              bool ignore_statistical_uncertainty)
                : fSystematicCovariance(systematic_covariance),
                  fNFunCalls(0),
                  fComponents(templates),
                  fIgnoreStatisticalUncertainty(ignore_statistical_uncertainty) {

            // check we have an appropriate covariance matrix
            assert(systematic_covariance.rows() == fComponents.GetNOuterBins() * fComponents.GetNInnerBins() &&
                   systematic_covariance.cols() == fComponents.GetNOuterBins() * fComponents.GetNInnerBins());

            fNUserParams = fComponents.GetNOuterBins() * fComponents.size();
            fNComponents = fComponents.size();

            fParamMap = detail::ParamMap(fNUserParams);
            fFixedParams = Eigen::RowVectorXd::Zero(fNUserParams);

            WarnInversionError();
        }

        void
        TemplateFitCalculator::
        AddNoise(double noise) {
            std::cout << "Info: Adding noise to the Covariance Diagonal: " << noise << std::endl;
            Matrix epsilon = (Vector::Ones(fSystematicCovariance.rows()) * noise).asDiagonal();
            fSystematicCovariance += epsilon;
            WarnInversionError();
        }

        Matrix
        TemplateFitCalculator::
        GetTotalCovariance(const Vector & params) const {
            Vector u = this->Predict(params);
            Matrix total_covariance = fSystematicCovariance;
            // add statistical uncertainty of reweighted prediction
            if(!fIgnoreStatisticalUncertainty) {
                total_covariance += u.asDiagonal();
            }
            return total_covariance;
        }

        void
        TemplateFitCalculator::
        WarnInversionError() const {
            Matrix I = Matrix::Identity(fSystematicCovariance.rows(),
                                        fSystematicCovariance.cols());
            fDecomp.compute(fSystematicCovariance);
            Matrix x = fDecomp.solve(I);
            double error = (fSystematicCovariance * x - I).norm() / I.norm();
            std::cout << "Info: Numerical Accuracy of Covariance Matrix Inversion = " << error << std::endl;
        }
    }
}


