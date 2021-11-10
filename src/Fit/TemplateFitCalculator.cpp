//
// Created by Derek Doyle on 10/9/21.
//

#include "XSecAna/Fit/TemplateFitCalculator.h"

namespace xsec {
    namespace fit {
        namespace detail {
            template<typename Scalar>
            struct arg_equal {
                void init(const Scalar & v, int i, int j) { return this->operator()(v, i, j); }

                void operator()(const Scalar & v, int i, int j) {
                    if (v == equal_to) {
                        row_where = i;
                        col_where = j;
                    }
                }

                double equal_to;
                int row_where = -1;
                int col_where = -1;
            };

            ParamMap::
            ParamMap(int size)
                    : fM(Matrix::Identity(size, size)) {}
            ParamMap::
            ParamMap(const Array & mask)
                    : ParamMap(mask.size())
            {
                for(auto i = 0; i < mask.size(); i++) {
                    if(!mask(i)) this->MaskTemplate(i);
                }
            }

            unsigned int
            ParamMap::
            GetNMinimizerParams() const {
                return fM.cols();
            }

            unsigned int
            ParamMap::
            GetNUserParams() const {
                return fM.rows();
            }

            Vector
            ParamMap::
            ToUserParams(const Vector & minimizer_params) const {
                return fM * minimizer_params;
            }

            Vector
            ParamMap::
            ToMinimizerParams(const Vector & user_params) const {
                return fM.transpose() * user_params;
            }

            void
            ParamMap::
            MaskTemplate(int i) {
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

            void
            ParamMap::
            UnmaskTemplate(int template_idx) {
                // check bounds
                assert(template_idx < fM.rows() &&
                       "Template index out of range");
                // if all templates are free, do nothing
                if (fM.rows() == fM.cols()) return;

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
                if (insert_at < 0) insert_at = fM.cols();
                assert(insert_at >= 0 && insert_at < fM.cols() + 1 &&
                       "could not determine where to insert column");

                fM.conservativeResize(fM.rows(), fM.cols() + 1);
                auto block_cols = fM.cols() - insert_at - 1;
                fM.block(0, insert_at + 1, fM.rows(), block_cols) =
                        fM.block(0, insert_at, fM.rows(), block_cols);
                fM.col(insert_at) = Eigen::VectorXd::Zero(fM.rows());
                fM(template_idx, insert_at) = 1;
            }

            bool
            ParamMap::
            IsParamMasked(int i) const {
                auto visitor = arg_equal<int>{1};
                fM.row(i).visit(visitor);

                // if -1, then visitor didn't find a match, and parameter is being masked
                return visitor.col_where == -1;
            }

            const Matrix &
            ParamMap::
            GetMatrix() const {
                return fM;
            }
        }

        double
        TemplateFitCalculator::
        Chi2(const Vector & user_params,
             const Array & data) const {
            auto u = this->Predict(user_params);
            auto v = data.matrix() - u;
            Matrix total_covariance = fSystematicCovariance;
            total_covariance += u.asDiagonal();
            return v.transpose() *
                    total_covariance.llt().solve(v);
        }

        double
        TemplateFitCalculator::
        fun(const Vector & minimizer_params,
            const Array & data) const {
            fNFunCalls++;
            return this->Chi2(this->ToUserParams(minimizer_params), data);

        }

        /// \brief User-level function for returning sum
        /// of all templates given the input template normalization
        /// parameters
        Vector
        TemplateFitCalculator::
        Predict(const Vector & user_params) const {
            return (fTemplates * user_params.asDiagonal())
                    .reshaped(fSystematicCovariance.rows(), fNComponents)
                    .rowwise().sum();
        }

        Vector
        TemplateFitCalculator::
        PredictComponent(const int & component_idx, const Vector & user_params) const {
            auto component_block = fTemplates.block(0,
                                                    fNOuterBins * component_idx,
                                                    fTemplates.rows(),
                                                    fNOuterBins);
            return (component_block * user_params.asDiagonal()).reshaped();
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
        TemplateFitCalculator(const std::vector<Array> & templates,
                              const std::vector<int> & dims,
                              const Matrix & systematic_covariance)
                : fDims(dims),
                  fSystematicCovariance(systematic_covariance),
                  fNFunCalls(0) {
            // check all templates are consistent with one another
            bool consistent_templates = true;
            for (auto i = 0u; i < templates.size(); i++) {
                consistent_templates &= templates[0].size() == templates[i].size();
            }
            assert(consistent_templates);

            // check we have an appropriate inverse covariance matrix
            assert(systematic_covariance.rows() == templates[0].size() &&
                   systematic_covariance.cols() == templates[0].size());

            // Determine number of user parameters
            // This number may differ from number of parameters
            // used in the minimization if the user wants to hold any
            // templates fixed.
            fNOuterBins = 1;
            for (auto i = 0u; i < dims.size() - 1; i++) {
                fNOuterBins *= dims[i];
            }
            fNInnerBins = dims[dims.size()-1];

            fNUserParams = fNOuterBins * templates.size();
            fNComponents = templates.size();

            // copy templates into an eigen matrix for fast linalg
            fTemplates = Eigen::MatrixXd(templates[0].size(),
                                         templates.size());
            for (auto i = 0u; i < templates.size(); i++) {
                fTemplates.col(i) = templates[i];
            }

            // resize to flatten the outer bins
            fTemplates.resize(fNInnerBins, fNUserParams);

            fParamMap = detail::ParamMap(fNUserParams);
            fFixedParams = Eigen::RowVectorXd::Zero(fNUserParams);
        }
    }
}


