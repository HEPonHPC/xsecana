//
// Created by Derek Doyle on 9/24/21.
//
#pragma once
#include <Eigen/Dense>
#include "XSecAna/IMeasurement.h"
#include <iomanip>
#include <iostream>

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
                               const Vector & data) const = 0;
            virtual unsigned int GetNMinimizerParams() const = 0;
            virtual unsigned int GetNUserParams() const = 0;
            virtual unsigned int GetNFunCalls() const = 0;
            virtual Vector ToUserParams(const Vector & minimizer_coords) const = 0;
            virtual Vector ToMinimizerParams(const Vector & user_coords) const = 0;

            std::vector<Array> GetRandomSeeds(int n, double lb, double ub) const {
                std::default_random_engine generator;
                std::uniform_real_distribution<double> distribution(lb, ub);

                std::vector<Array> seeds(n);
                for(auto ithrow = 0u; ithrow < n; ithrow++) {
                    seeds[ithrow] = Array::Ones(this->GetNMinimizerParams());
                    for (auto iparam = 0u; iparam < seeds[ithrow].size(); iparam++) {
                        seeds[ithrow](iparam) += distribution(generator);
                    }
                }
                return seeds;
            }

            std::vector<Vector> GetBestSeeds(const Vector & data,
                                             const std::vector<Array> & seeds,
                                             int take) const {
                std::vector<double> chi2(seeds.size());
                std::vector<int> seed_idx(seeds.size());
                for(auto i = 0u; i < seeds.size(); i++) {
                    seed_idx[i] = i;
                    chi2[i] = this->fun(seeds[i], data);
                }

                // sort seed_idx based on ascending chi2 values
                std::sort(seed_idx.begin(), seed_idx.end(),
                          [&chi2](int a, int b) -> bool {
                              return chi2[a] < chi2[b];
                          });

                std::cout << "Seeding fit with the following parameter configurations:" << std::endl;
                std::cout << "--------------------------------------------------------" << std::endl;
                std::cout << std::left;
                std::cout << std::setw(20) << "chi-squared" << "|" << " (seed - 1) l^2 norm / N" << std::endl;
                std::cout << "--------------------------------------------------------" << std::endl;
                std::vector<Vector> best_seeds(take);
                for (auto iseed = 0u; iseed < take; iseed++) {
                    best_seeds[iseed] = seeds[seed_idx[iseed]];
                    std::cout << std::setw(20) << chi2[seed_idx[iseed]] << "| ";
                    std::cout << (best_seeds[iseed] - Vector::Ones(best_seeds[iseed].size())).lpNorm<2>() / best_seeds[iseed].size()<< std::endl;
                }
                return best_seeds;
            }
        };

        class IFitter {
        public:
            virtual FitResult Fit(IFitCalculator * fit_calc,
                                  const Vector & data,
                                  const std::vector<Vector> seeds={}) = 0;
        };



    }
}