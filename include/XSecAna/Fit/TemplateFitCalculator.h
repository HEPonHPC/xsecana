#pragma once

#include "XSecAna/Fit/IFitter.h"
#include "XSecAna/Fit/TemplateFitComponent.h"
#include "XSecAna/Systematic.h"
#include <Eigen/Dense>
#include <cmath>
#include <vector>

namespace xsec {
    namespace fit {
        class TemplateFitCalculator : public IFitCalculator {
        public:

            TemplateFitCalculator(const ReducedComponentCollection & templates,
                                  const std::vector<int> & dims,
                                  const Matrix & systematic_covariance,
                                  bool ignore_statistical_uncertainty = false);

            void FixTemplate(const int & template_idx, const double & at = 1);

            void ReleaseTemplate(const int & template_idx);

            unsigned int GetNTemplates() const { return fNUserParams; }
            unsigned int GetNComponents() const { return fNComponents; }
            unsigned int GetNOuterBins() const { return fComponents.GetNOuterBins(); }
            unsigned int GetNInnerBins() const { return fComponents.GetNInnerBins(); }

            unsigned int GetNFunCalls() const override { return fNFunCalls; }

            unsigned int GetNMinimizerParams() const override { return fParamMap.GetNMinimizerParams(); }
            unsigned int GetNUserParams() const override { return fParamMap.GetNUserParams(); }

            const detail::ParamMap & GetParamMap() const { return fParamMap; }

            virtual Vector Predict(const Vector & user_params) const;
            virtual Vector PredictComponent(const std::string & component_label, const Vector & user_params) const;

            double Chi2(const Vector & user_params,
                        const Vector & data) const;

            Vector ToUserParams(const Vector & minimizer_coords) const override;
            Vector ToMinimizerParams(const Vector & user_coords) const override;
            double fun(const Vector & minimizer_params,
                       const Vector & data) const override;
                                                              ;

            void AddNoise(double noise);

            Matrix GetSystematicCovariance() const { return fSystematicCovariance; }
            Matrix GetTotalCovariance(const Vector & params) const;

            void SetIgnoreStatisticalUncertainty(bool ignore) { fIgnoreStatisticalUncertainty = ignore; }
            double GetIgnoreStatisticalUncertainty() const { return fIgnoreStatisticalUncertainty; }

        private:
            void SetSystematicDeterminant();
            static double LogDetV(const Eigen::LLT<Matrix> & decomp);

            void WarnInversionError() const;

            Matrix fSystematicCovariance;
            ReducedComponentCollection fComponents;
            int fNUserParams;
            int fNComponents;

            detail::ParamMap fParamMap;
            Vector fFixedParams;

            mutable unsigned int fNFunCalls;

            bool fIgnoreStatisticalUncertainty;

            mutable Eigen::LLT<Matrix> fDecomp;
            double fDetSystematicCovariance;
            //mutable Eigen::PartialPivLU<Matrix> fDecomp;

            // Box-Cox Transform Variables
            mutable bool fIsBCSetup = false;
            mutable Vector fBCData;
            mutable Matrix fBCSMatrix;
            mutable double fBCLogSMatrixDeterminant;

            // Normalized Covariance Variables
            mutable bool fIsNCVSetup = false;
            mutable Vector fNCVNominal;
            mutable Matrix fNCVSystematicCovariance;
        };
    }
}
