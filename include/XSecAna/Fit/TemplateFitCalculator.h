#pragma once

#include "XSecAna/Fit/IFitter.h"
#include <Eigen/Dense>
#include <math.h>
#include <vector>

namespace xsec {
    namespace fit {
        namespace detail {
            template<class Scalar>
            struct arg_equal {
                void init(const Scalar & v, int i, int j) { return this->operator()(v, i, j); }

                void operator()(const Scalar & v, int i, int j) {
                    if (v == equal_to) {
                        row_where = i;
                        col_where = j;
                    }
                }

                Scalar equal_to;
                int row_where = -1;
                int col_where = -1;
            };

            class ParamMap {
            public:
                ParamMap() = default;

                ParamMap(int size)
                        : fM(Eigen::MatrixXd::Identity(size, size)) {}

                unsigned int GetNMinimizerParams() const { return fM.cols(); }

                Eigen::VectorXd ToUserParams(const Eigen::VectorXd & minimizer_params) const {
                    return fM * minimizer_params;
                }

                Eigen::VectorXd ToMinimizerParams(const Eigen::VectorXd & user_params) const {
                    return fM.transpose() * user_params;
                }

                void MaskTemplate(int i) {
                    // if this parameter is already being masked, don't do anything
                    if (IsParamMasked(i)) return;

                    // find column that maps param i and remove it
                    auto visitor = arg_equal<int>{1};

                    fM.row(i).visit(visitor);
                    // remove column visitor.col_where
                    auto new_nrows = fM.rows();
                    auto new_ncols = fM.cols() - 1;

                    if (visitor.col_where < new_ncols) {
                        fM.block(0, visitor.col_where, new_nrows, new_ncols - visitor.col_where) =
                                fM.block(0, visitor.col_where + 1, new_nrows, new_ncols - visitor.col_where);
                    }
                    fM.conservativeResize(new_nrows, new_ncols);
                }

                void UnmaskTemplate(int template_idx) {
                    // check bounds
                    assert(template_idx < fM.rows() &&
                           "Template index out of range");
                    // if all templates are free, do nothing
                    if(fM.rows() == fM.cols()) return;

                    // insert row to retain ordering
                    auto visitor = arg_equal<int>{1};
                    auto insert_at = -1;
                    for (auto i = 0; i < fM.cols(); i++) {
                        fM.col(i).visit(visitor);

                        if (visitor.row_where > template_idx) {
                            insert_at = i;
                            break;
                        }
                    }
                    if(insert_at < 0) insert_at = fM.cols();
                    assert(insert_at >= 0 && insert_at < fM.cols() + 1 &&
                           "could not determine where to insert column");

                    fM.conservativeResize(fM.rows(), fM.cols() + 1);
                    auto block_cols = fM.cols() - insert_at - 1;
                    fM.block(0, insert_at + 1, fM.rows(), block_cols) =
                            fM.block(0, insert_at, fM.rows(), block_cols);
                    fM.col(insert_at) = Eigen::VectorXd::Zero(fM.rows());
                    fM(template_idx, insert_at) = 1;
                }

                bool IsParamMasked(int i) const {
                    auto visitor = arg_equal<int>{1};
                    fM.row(i).visit(visitor);

                    // if -1, then visitor didn't find a match, and parameter is being masked
                    return visitor.col_where == -1;
                }

                const Eigen::MatrixXd & GetMatrix() const { return fM; }

            private:
                Eigen::MatrixXd fM;
            };
        }

        template<class Scalar,
                int Cols>
        class TemplateFitCalculator : public IFitCalculator<Scalar, Cols> {
        public:
            typedef Eigen::Matrix<Scalar, Cols, Cols> covariance_matrix_type;

            TemplateFitCalculator(const std::vector<Eigen::Array<Scalar, 1, Cols>> & templates,
                                  std::vector<int> dims,
                                  covariance_matrix_type & inverse_covariance);

            void FixTemplate(const int & template_idx, const double & at = 1);

            void ReleaseTemplate(const int & template_idx);

            unsigned int GetNTemplates() const { return fNUserParams; }

            unsigned int GetNComponents() const { return fNComponents; }
            virtual unsigned int GetNFunCalls() const override { return fNFunCalls; }

            virtual unsigned int GetNMinimizerParams() const override { return fParamMap.GetNMinimizerParams(); }

            Eigen::VectorXd Predict(const Eigen::VectorXd & user_params) const;
            double Chi2(const Eigen::VectorXd & user_params,
                        const Eigen::Array<Scalar, 1, Cols> & data) const;

            Eigen::VectorXd ToUserParams(const Eigen::VectorXd & minimizer_coords) const;

            virtual double fun(const Eigen::VectorXd & minimizer_params,
                               const Eigen::Array<Scalar, 1, Cols> & data) const override;

        private:

            Eigen::VectorXd ToMinimizerParams(const Eigen::VectorXd & user_coords) const;

            covariance_matrix_type fInvCov;
            std::vector<int> fDims;
            Eigen::MatrixXd fTemplates;
            int fNUserParams;
            int fNComponents;

            detail::ParamMap fParamMap;
            Eigen::VectorXd fFixedParams;

            mutable unsigned int fNFunCalls;
        };

        template<class Scalar,
                int Cols>
        double
        TemplateFitCalculator<Scalar, Cols>::
        Chi2(const Eigen::VectorXd & user_params,
             const Eigen::Array<Scalar, 1, Cols> & data) const {
            return (data.matrix() - this->Predict(user_params).transpose()) *
                   fInvCov *
                   (data.matrix().transpose() - this->Predict(user_params));
        }

        template<class Scalar,
                int Cols>
        double
        TemplateFitCalculator<Scalar, Cols>::
        fun(const Eigen::VectorXd & minimizer_params,
            const Eigen::Array<Scalar, 1, Cols> & data) const {
            fNFunCalls++;
            return this->Chi2(this->ToUserParams(minimizer_params), data);

        }

        /// \brief User-level function for returning sum
        /// of all templates given the input template normalization
        /// parameters
        template<class Scalar,
                int Cols>
        Eigen::VectorXd
        TemplateFitCalculator<Scalar, Cols>::
        Predict(const Eigen::VectorXd & user_params) const {
            return (fTemplates * user_params.asDiagonal())
                    .template reshaped(fInvCov.rows(), fNComponents)
                    .rowwise().sum();
        }

        template<class Scalar,
                int Cols>
        void
        TemplateFitCalculator<Scalar, Cols>::
        ReleaseTemplate(const int & template_idx) {
            fFixedParams(template_idx) = 0;
            fParamMap.UnmaskTemplate(template_idx);
        }

        template<class Scalar,
                int Cols>
        void
        TemplateFitCalculator<Scalar, Cols>::
        FixTemplate(const int & template_idx, const double & at) {
            fFixedParams(template_idx) = at;
            fParamMap.MaskTemplate(template_idx);
        }

        template<class Scalar,
                int Cols>
        Eigen::VectorXd
        TemplateFitCalculator<Scalar, Cols>::
        ToUserParams(const Eigen::VectorXd & minimizer_coords) const {
            return fParamMap.ToUserParams(minimizer_coords) + fFixedParams;
        }

        template<class Scalar,
                int Cols>
        Eigen::VectorXd
        TemplateFitCalculator<Scalar, Cols>::
        ToMinimizerParams(const Eigen::VectorXd & user_coords) const {
            return fParamMap.ToMinimizerParams(user_coords);
        }

        template<class Scalar,
                int Cols>
        TemplateFitCalculator<Scalar, Cols>::
        TemplateFitCalculator(const std::vector<Eigen::Array<Scalar, 1, Cols>> & templates,
                              std::vector<int> dims,
                              covariance_matrix_type & inverse_covariance)
                : fDims(dims),
                  fInvCov(inverse_covariance),
                  fNFunCalls(0) {
            // check all templates are consistent with one another
            bool consistent_templates = true;
            for (auto i = 0u; i < templates.size(); i++) {
                consistent_templates &= templates[0].size() == templates[i].size();
            }
            assert(consistent_templates);

            // check we have an appropriate inverse covariance matrix
            assert(inverse_covariance.rows() == templates[0].size() &&
                   inverse_covariance.cols() == templates[0].size());

            // Determine number of user parameters
            // This number may differ from number of parameters
            // used in the minimization if the user wants to hold any
            // templates fixed.
            auto n_outer_bins = 1;
            for (auto i = 0u; i < dims.size() - 1; i++) {
                n_outer_bins *= dims[i];
            }
            fNUserParams = n_outer_bins * templates.size();
            fNComponents = templates.size();

            // copy templates into an eigen matrix for fast linalg
            fTemplates = Eigen::MatrixXd(templates[0].size(),
					 templates.size());
            for (auto i = 0u; i < templates.size(); i++) {
                fTemplates.col(i) = templates[i];
            }

            // resize to flatten the outer bins
            fTemplates.resize(dims[dims.size() - 1], fNUserParams);

            fParamMap = detail::ParamMap(fNUserParams);
            fFixedParams = Eigen::RowVectorXd::Zero(fNUserParams);
        }
    }
}
