#pragma once

#include <Eigen/Dense>
#include "XSecAna/Fit/IFitter.h"

#include "Minuit2/FCNBase.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnPrint.h"

#include "RVersion.h"

namespace xsec {
    namespace fit {
        namespace detail {
            Eigen::Map<const Vector> STDToEigen(const std::vector<double> & v) {
                return Eigen::Map<const Vector>(v.data(), v.size(), 1);
            }
            std::vector<double> EigenToSTD(const Vector & v) {
                return std::vector<double>(v.data(), v.data() + v.size());
            }
        }
        class Minuit2TemplateFitter : public IFitter, public ROOT::Minuit2::FCNBase {
        public:
            Minuit2TemplateFitter(int strategy = 2,
                                  double up = 1,
                                  bool minos_errors = true)
                    : fMnStrategy(strategy),
                      fUp(up),
                      fMinosErrors(minos_errors)
            {}

            // IFitter overrides
            virtual FitResult Fit(IFitCalculator * fit_calc,
                                  const Vector & data,
                                  std::vector<Vector> seeds = {}) override;

            // FCNBase overrides
            virtual double operator()(const std::vector<double> & params) const override;
            virtual double Up() const override;

            void SetPrintLevel(const int & level) const;

        private:
            IFitCalculator * fFitCalc;
            int fMnStrategy;
            double fUp;
            bool fMinosErrors;
            Vector fData;
        };

    }
}