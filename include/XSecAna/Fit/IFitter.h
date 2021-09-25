//
// Created by Derek Doyle on 9/24/21.
//
#pragma once
#include <Eigen/Dense>

namespace xsec {
    namespace fit {
        struct FitResult {
            double fun_val;
            Eigen::VectorXd params;
            Eigen::MatrixXd covariance;
            Eigen::VectorXd plus_one_sigma_errors;
            Eigen::VectorXd minus_one_sigma_errors;
            unsigned int fun_calls;
        };

        template<class Scalar=double,
                int Cols=Eigen::Dynamic>
        class IFitCalculator {
        public:
            virtual double fun(const Eigen::VectorXd & params,
                               const Eigen::Array<Scalar, 1, Cols> & data) const = 0;
            virtual unsigned int GetNMinimizerParams() const = 0;
            virtual unsigned int GetNFunCalls() const = 0;
        };

        template<class Scalar=double,
                int Cols=Eigen::Dynamic>
        class IFitter {
        public:
            virtual FitResult Fit(const Eigen::Array<Scalar, 1, Cols> & data,
                                  const std::vector<Eigen::VectorXd> seeds) = 0;
        };



    }
}