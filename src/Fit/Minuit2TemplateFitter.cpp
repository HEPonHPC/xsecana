//
// Created by Derek Doyle on 10/9/21.
//

#include "XSecAna/Fit/Minuit2TemplateFitter.h"

namespace xsec {
    namespace fit {

        double
        Minuit2TemplateFitter::
        Up() const {
            return fUp;
        }


        void
        Minuit2TemplateFitter::
        SetPrintLevel(const int & level) const {
#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 24, 0)
            ROOT::Minuit2::MnPrint::SetGlobalLevel(level);
#elif ROOT_VERSION_CODE >= ROOT_VERSION(6, 22, 00)
            ROOT::Minuit2::MnPrint::SetLevel(level);
#else
            std::cerr << "Attempting to set Minuit2 print level for ROOT version ";
            std::cerr << ROOT_VERSION_CODE;
            std::cerr << ", but I'm not yet sure how to do that. Provide implementation in ";
            std::cerr << __PRETTY_FUNCTION__ << std::endl;
#endif
        }

        FitResult
        Minuit2TemplateFitter::
        Fit(const Array & data,
            std::vector<Vector> seeds) {
            fData = data;
            // setup minuit fitter
            fFitCalc->GetNMinimizerParams();

            // if user hasn't provided any starting positions, we'll start at 1 for all params
            if (seeds.size() == 0) {
                seeds.push_back(Eigen::VectorXd::Ones(fFitCalc->GetNMinimizerParams()));
            }

            std::vector<ROOT::Minuit2::FunctionMinimum> mins;
            for (auto seed: seeds) {
                // minuit2 requires parameters to be initialized with errors
                // 10% seems like a reasonable starting point
                ROOT::Minuit2::MnUserParameters mn_params;
                for (auto i = 0u; i < seed.size(); i++) {
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
            for (auto i = 0u; i < mins.size(); i++) {
                if (mins[global_min_idx].Fval() > mins[i].Fval()) {
                    global_min_idx = i;
                }
            }

            // get errors
            ROOT::Minuit2::MnMinos minos(*this, mins[global_min_idx], fMnStrategy);
            Vector plus_one_sigma_errors(fFitCalc->GetNMinimizerParams());
            Vector minus_one_sigma_errors(fFitCalc->GetNMinimizerParams());
            Vector best_fit_params(fFitCalc->GetNMinimizerParams());
            for (auto i = 0u; i < fFitCalc->GetNMinimizerParams(); i++) {
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
            result.covariance = Eigen::MatrixXd::Zero(fFitCalc->GetNMinimizerParams(),
                                                      fFitCalc->GetNMinimizerParams());
            for (auto irow = 0u; irow < fFitCalc->GetNMinimizerParams(); irow++) {
                for (auto icol = 0u; icol < fFitCalc->GetNMinimizerParams(); icol++) {
                    result.covariance(irow, icol) = mins[global_min_idx].UserCovariance()(irow, icol);
                }
            }

            return result;
        }

        double
        Minuit2TemplateFitter::
        operator()(const std::vector<double> & params) const {
            return fFitCalc->fun(detail::STDToEigen(params), fData);
        }


    }
}

