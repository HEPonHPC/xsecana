#pragma once

#include <Eigen/Dense>
#include "XSecAna/Fit/IFitter.h"

#include "Minuit2/FCNBase.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnMigrad.h"

namespace xsec {
    namespace fit {
        namespace detail {
            Eigen::Map<const Eigen::VectorXd> STDToEigen(const std::vector<double> & v) {
                return Eigen::Map<const Eigen::VectorXd>(&v[0], v.size(), 1);
            }
        }
        template<class Scalar=double,
                int Cols=Eigen::Dynamic>
        class Minuit2TemplateFitter : public IFitter<Scalar, Cols>,
                              public ROOT::Minuit2::FCNBase {
        public:
            Minuit2TemplateFitter(IFitCalculator<Scalar, Cols> * fFitCalc,
                                  int strategy = 2,
                                  double up = 1)
                          : fFitCalc(fFitCalc),
                            fMnStrategy(strategy),
                            fUp(up)
            {}

            // IFitter overrides
            virtual FitResult Fit(const Eigen::Array<Scalar, 1, Cols> & data,
                                  std::vector<Eigen::VectorXd> seeds = {}) override;

            // FCNBase overrides
            virtual double operator()(const std::vector<double> & params) const override;
            virtual double Up() const override { return fUp; }

            void SetPrintLevel(int level) { ROOT::Minuit2::MnPrint::SetLevel(level); }

        private:
            IFitCalculator<Scalar, Cols> * fFitCalc;
            int fMnStrategy;
            double fUp;
            Eigen::Array<Scalar, 1, Cols> fData;

        };

        template<class Scalar,
                int Cols>
        FitResult
        Minuit2TemplateFitter<Scalar, Cols>::
        Fit(const Eigen::Array<Scalar, 1, Cols> & data,
            std::vector<Eigen::VectorXd> seeds) {
            fData = data;
            // setup minuit fitter
            fFitCalc->GetNMinimizerParams();

            // if user hasn't provided any starting positions, we'll start at 1 for all params
            if(seeds.size() == 0) {
                seeds.push_back(Eigen::VectorXd::Ones(fFitCalc->GetNMinimizerParams()));
            }

            std::vector<ROOT::Minuit2::FunctionMinimum> mins;
            for(auto seed : seeds) {
                // minuit2 requires parameters to be initialized with errors
                // 10% seems like a reasonable starting point
                ROOT::Minuit2::MnUserParameters mn_params;
                for(auto i = 0u; i < seed.size(); i++) {
                    // can we name these better?
                    mn_params.Add(std::to_string(i),
                                  seed(i),
                                  0.01, // 10% error
                                  0,   // lower bound: 0
                                  20);  // upper bound: 20. Hopefully our MC isn't that bad
                }
                // call minimizer with this
                ROOT::Minuit2::MnMigrad minimizer(*this, mn_params, fMnStrategy);

                // fit runs here.
                // Save all the results and find the best one after.
                // We'll need the best ROOT::Minuit2::FunctionMinimum
                // to find errors with ROOT::Minuit2::MnMinos
                mins.push_back(minimizer());
            }

            // find best fit
            auto global_min_idx = 0u;
            for(auto i = 0u; i < mins.size(); i++) {
                if(mins[global_min_idx].Fval() > mins[i].Fval()) {
                    global_min_idx = i;
                }
            }

            // get errors
            ROOT::Minuit2::MnMinos minos(*this, mins[global_min_idx], fMnStrategy);
            Eigen::VectorXd plus_one_sigma_errors(fFitCalc->GetNMinimizerParams());
            Eigen::VectorXd minus_one_sigma_errors(fFitCalc->GetNMinimizerParams());
            Eigen::VectorXd best_fit_params(fFitCalc->GetNMinimizerParams());
            for(auto i = 0u; i < fFitCalc->GetNMinimizerParams(); i++) {
                // each call to minos does a fit so this part will take some time
                auto e = minos(i);
                minus_one_sigma_errors(i) = std::get<0>(e);
                plus_one_sigma_errors(i) = std::get<1>(e);
                best_fit_params(i) = mins[global_min_idx].UserState().Params()[i];
            }

            // return as FitResult
            FitResult result;
            result.plus_one_sigma_errors = plus_one_sigma_errors;
            result.minus_one_sigma_errors = minus_one_sigma_errors;
            result.fun_val = mins[global_min_idx].Fval();
            result.params = best_fit_params;
            result.fun_calls = fFitCalc->GetNFunCalls();
            return result;
        }

        template<class Scalar,
                int Cols>
        double
        Minuit2TemplateFitter<Scalar, Cols>::
        operator()(const std::vector<double> & params) const {
            return fFitCalc->fun(detail::STDToEigen(params), fData);
        }

    }
}