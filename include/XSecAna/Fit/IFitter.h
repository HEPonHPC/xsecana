//
// Created by Derek Doyle on 9/24/21.
//
#pragma once
#include <Eigen/Dense>
#include "XSecAna/IMeasurement.h"

namespace xsec {
    namespace fit {
        typedef Eigen::VectorXd Vector;
        typedef Eigen::MatrixXd Matrix;
        typedef Eigen::Map<Eigen::VectorXd> VectorMap;
        typedef Eigen::Map<Eigen::MatrixXd> MatrixMap;

        struct FitResult {
            bool is_valid;
            double fun_val;
            Vector params;
            Array2D covariance;
            Array params_error_up;
            Array params_error_down;
            unsigned int fun_calls;
        };

        class IFitCalculator {
        public:
            virtual double fun(const Vector & params,
                               const Array & data) const = 0;
            virtual unsigned int GetNMinimizerParams() const = 0;
            virtual unsigned int GetNFunCalls() const = 0;
            virtual Vector ToUserParams(const Vector & minimizer_coords) const = 0;
        };

        class IFitter {
        public:
            virtual FitResult Fit(IFitCalculator * fit_calc,
                                  const Array & data,
                                  const std::vector<Vector> seeds={}) = 0;
        };



    }
}