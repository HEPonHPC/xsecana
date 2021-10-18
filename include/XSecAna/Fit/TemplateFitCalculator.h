#pragma once

#include "XSecAna/Fit/IFitter.h"
#include <Eigen/Dense>
#include <math.h>
#include <vector>

namespace xsec {
    namespace fit {
        namespace detail {
            class ParamMap {
            public:
                ParamMap() = default;
                ParamMap(int size);
                ParamMap(const Array & mask);

                unsigned int GetNMinimizerParams() const;
                Vector ToUserParams(const Eigen::VectorXd & minimizer_params) const;
                Vector ToMinimizerParams(const Vector & user_params) const;
                void MaskTemplate(int i);
                void UnmaskTemplate(int template_idx);
                bool IsParamMasked(int i) const;
                const Matrix & GetMatrix() const;
            private:
                Matrix fM;
            };
        }

        class TemplateFitCalculator : public IFitCalculator {
        public:

            TemplateFitCalculator(const std::vector<Array> & templates,
                                  const std::vector<int> & dims,
                                  const Matrix & inverse_covariance);

            void FixTemplate(const int & template_idx, const double & at = 1);

            void ReleaseTemplate(const int & template_idx);

            unsigned int GetNTemplates() const { return fNUserParams; }
            unsigned int GetNComponents() const { return fNComponents; }
            unsigned int GetNOuterBins() const { return fNOuterBins; }
            unsigned int GetNInnerBins() const { return fNInnerBins; }

            virtual unsigned int GetNFunCalls() const override { return fNFunCalls; }

            virtual unsigned int GetNMinimizerParams() const override { return fParamMap.GetNMinimizerParams(); }

            Vector Predict(const Vector & user_params) const;

            double Chi2(const Vector & user_params,
                        const Array & data) const;

            Vector ToUserParams(const Vector & minimizer_coords) const;

            virtual double fun(const Vector & minimizer_params,
                               const Array & data) const override;

        private:

            Vector ToMinimizerParams(const Vector & user_coords) const;

            const Matrix fInvCov;
            const std::vector<int> fDims;
            Matrix fTemplates;
            int fNUserParams;
            int fNComponents;

            detail::ParamMap fParamMap;
            Vector fFixedParams;

            mutable unsigned int fNFunCalls;
            int fNOuterBins;
            int fNInnerBins;
        };


    }
}
