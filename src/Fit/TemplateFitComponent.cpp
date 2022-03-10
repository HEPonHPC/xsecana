#include "XSecAna/Fit/TemplateFitComponent.h"

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
                    : ParamMap(mask.size()) {
                for (auto i = 0; i < mask.size(); i++) {
                    if (!mask(i)) this->MaskTemplate(i);
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

        //ReducedComponent::
        //ReducedComponent(const TH1 * h)
        //        : fHist(h), fArray(root::MapContentsToEigen(h)) {
        //    assert(h->GetDimension() != 3 && "Histogram not reduced");
        //    if (h->GetDimension() == 1) {
        //        fDims = {1, h->GetNbinsX() + 2};
        //    } else if (h->GetDimension() == 2) {
        //        fDims = {h->GetNbinsY() + 2, h->GetNbinsX() + 2};
        //    }
        //}

        ReducedComponent::
        ReducedComponent(const TH1 * h, int nouter_bins, int ninner_bins)
                : fArray(root::MapContentsToEigen(h)), fNOuterBins(nouter_bins), fNInnerBins(ninner_bins) {}

        ReducedComponent::
        ReducedComponent(const Vector & a, int nouter_bins, int ninner_bins)
                : fArray(a), fNOuterBins(nouter_bins), fNInnerBins(ninner_bins) {
            fHist = std::make_shared<TH1D>("", "", a.size(), 0, a.size());
            for(int i = 0; i < a.size(); i++) {
                fHist->SetBinContent(i+1, a(i));
            }
        }

        Array
        ReducedComponent::
        Project() const {
            return fArray.reshaped(fNInnerBins, fNOuterBins).colwise().sum();
        }


        Vector
        ReducedTemplateComponent::
        Predict(const Vector & component_params) const {
            return (fMean->GetArray().reshaped(fMean->GetNInnerBins(), fMean->GetNOuterBins()) *
                   component_params.asDiagonal()).reshaped();
        }

        Vector
        ReducedTemplateComponent::
        PredictProjected(const Vector & component_params) const {
            Array nominal_projection = fMean->GetArray().reshaped(fMean->GetNInnerBins(),
                                                                  fMean->GetNOuterBins())
                    .colwise().sum().reshaped();
            return  nominal_projection * component_params.array();
        }

        IReducedTemplateComponent *
        UserTemplateComponent::
        Reduce(const ComponentReducer & reducer) const {
            std::map<std::string, Systematic<TH1>> reduced_systematics;
            for(const auto & syst : fSystematics) {
                reduced_systematics[syst.first] = reducer.Reduce(syst.second);
            }
            return new ReducedTemplateComponent(reducer.Reduce(fMean),
                                                reduced_systematics);
        }

        IReducedTemplateComponent *
        UserTemplateComponentAlternativeNominal::
        Reduce(const ComponentReducer & reducer) const {
            auto reduced_component = UserTemplateComponent::Reduce(reducer);
            return new ReducedTemplateComponentAlternativeNominal(reduced_component->GetNominal(),
                                                                  reduced_component->GetSystematics(),
                                                                  reducer.Reduce(fMeanForErrorCalc));
        }

        ComponentReducer::
        ComponentReducer(const TH1 * mask)
                : fMask(mask),
                  fMap(fit::detail::ParamMap(root::MapContentsToEigen(mask))) {}

        xsec::Systematic<TH1>
        ComponentReducer::
        Reduce(const Systematic<TH1> & syst) const {
            std::vector<std::shared_ptr<TH1>> reduced;
            for(auto i = 0u; i < syst.GetShifts().size(); i++) {
                reduced.push_back(this->Reduce(syst.GetShifts()[i])->GetHist());
            }
            return Systematic<TH1>(syst.GetName(), reduced, syst.GetType());
        }

        TH1 *
        ComponentReducer::
        Project(const std::shared_ptr<TH1> templ) {
            TH1 * ret;
            if (templ->GetDimension() == 1) {
                ret = new TH1D("", "", 1, 0, 1);
                double error;
                double integral = templ->IntegralAndError(1, templ->GetNbinsX(), error);
                ret->SetBinContent(1, integral);
                ret->SetBinError(1, error);
            } else if (templ->GetDimension() == 2) {
                ret = ((TH2 *) templ.get())->ProjectionX("",
                                                         1,
                                                         templ->GetNbinsY(),
                                                         "e");
            } else {
                // project into xy plane
                // NOF: no over flow
                // NUF: no under flow
                // e: calculate errors
                ret = ((TH3 *) templ.get())->Project3D("yx NOF NUF e");
            }
            // set new name or else ROOT will delete the histogram under our feet
            ret->SetName(root::MakeUnique("templ_project").c_str());
            return ret;
        }

        Systematic<TH1>
        ComponentReducer::
        Project(const Systematic<TH1> & syst) {
            std::vector<std::shared_ptr<TH1>> projected(syst.GetShifts().size());
            for(auto i = 0u; i < syst.GetShifts().size(); i++) {
                projected[i] = std::shared_ptr<TH1>(ComponentReducer::Project(syst.GetShifts()[i]));
            }
            return Systematic<TH1>(syst.GetName(), projected, syst.GetType());
        }

        TH1 *
        ComponentReducer::
        Compress1D(const std::shared_ptr<TH1> component) const {
            Array compressed_c = Array::Zero(fMap.GetNMinimizerParams()+2);
            compressed_c(Eigen::seqN(1, fMap.GetNMinimizerParams())) =
                    fMap.ToMinimizerParams(root::MapContentsToEigen(component.get()));
            Array compressed_e = Array::Zero(fMap.GetNMinimizerParams()+2);
            compressed_e(Eigen::seqN(1, fMap.GetNMinimizerParams())) =
                    fMap.ToMinimizerParams(root::MapErrorsToEigen(component.get()));

            auto ret = new TH1D("", "", fMap.GetNMinimizerParams(), 0, fMap.GetNMinimizerParams());
            return root::ToROOTLike(ret, compressed_c, compressed_e);
        }

        ReducedComponent *
        ComponentReducer::
        Reduce(const std::shared_ptr<TH1> component) const {
            std::vector<int> dims(2);
            if (component->GetDimension() == 1) {
                std::cout << "Warning: Attempting to apply to mask 1-dimensional template fit" << std::endl;
                return new ReducedComponent((TH1*) component->Clone(),
                                            1,
                                            component->GetNbinsX());
                        ;
            } else if (component->GetDimension() == 2) {
                assert(fMask->GetDimension() == 1);
                // not including under/overflow. Some may want to?
                int nunmasked = fMask->Integral(0, fMask->GetNbinsX() + 1);
                auto ret = new TH2D("", "",
                                    component->GetNbinsY() - 2, 1, component->GetNbinsY(),
                                    nunmasked - 2, 1, nunmasked - 1);
                auto ii = 0;
                for (auto i = 0u; i < fMask->GetNbinsX() + 2; i++) {
                    if (fMask->GetBinContent(i)) {
                        for (auto j = 1; j <= component->GetNbinsY(); j++) {
                            ret->SetBinContent(ii, j - 1, component->GetBinContent(i, j));
                            ret->SetBinError(ii, j - 1, component->GetBinError(i, j));
                        }
                        ii++;
                    }
                }
                ret->GetXaxis()->SetTitle("Outer Bins");
                ret->GetYaxis()->SetTitle("Template Bins");
                return new ReducedComponent(root::MapContentsToEigen(ret),
                                            ret->GetNbinsY()+2,
                                            ret->GetNbinsX()+2);
            } else {
                assert(fMask->GetDimension() == 2);
                int nunmasked = ((TH2 *) fMask)->Integral(0, fMask->GetNbinsX() + 1,
                                                          0, fMask->GetNbinsY() + 1);

                // using under/overflow bins here in the flattened histogram
                // for easy un-packing into 1D eigen array that we'll hand to the fitter
                // We do this instead of converting to Eigen because of limitations imposed by
                // the current implementation of Systematic
                auto ret = new TH2D("", "",
                                    component->GetNbinsZ() - 2, 1, component->GetNbinsZ(),
                                    nunmasked - 2, 1, nunmasked - 1);
                // Not sure why we have to loop over y first,
                // but its necessary for things to get ordered properly
                // during the un-masking
                auto ii = 0;
                for (auto j = 0u; j < fMask->GetNbinsY() + 2; j++) {
                    for (auto i = 0u; i < fMask->GetNbinsX() + 2; i++) {
                        if (fMask->GetBinContent(i, j)) {
                            for (auto k = 1; k <= component->GetNbinsZ(); k++) {
                                ret->SetBinContent(k - 1, ii, component->GetBinContent(i, j, k));
                                ret->SetBinError(k - 1, ii, component->GetBinError(i, j, k));
                            }
                            ii++;
                        }
                    }
                }
                ret->GetXaxis()->SetTitle("Template Bins");
                ret->GetYaxis()->SetTitle("Outer Bins");
                return new ReducedComponent(root::MapContentsToEigen(ret),
                                            ret->GetNbinsY()+2,
                                            ret->GetNbinsX()+2);
            }
        }



        ReducedComponentCollection
        UserComponentCollection::
        Reduce(const ComponentReducer & reducer) const {
            std::map<std::string, const IReducedTemplateComponent *> reduced_components;
            for (const auto & component: fComponents) {
                reduced_components[component.first] = component.second->Reduce(reducer);
            }
            return ReducedComponentCollection(reduced_components);
        }

        TH1 *
        UserComponentCollection::
        NominalTotal() const {
            TH1 * total = nullptr;
            for(const auto & component : fComponents) {
                if (!total) {
                    total = (TH1*) component.second->GetNominal()->Clone();
                }
                else {
                    total->Add(component.second->GetNominal().get());
                }
            }
            return total;
        }

        TH1 *
        UserComponentCollection::
        NominalTotalForErrorCalculation() const {
            TH1 * total = nullptr;
            for(const auto & component : fComponents) {
                if (!total) {
                    total = (TH1*) component.second->GetNominalForErrorCalculation()->Clone();
                }
                else {
                    total->Add(component.second->GetNominalForErrorCalculation().get());
                }
            }
            return total;
        }

        TH1 *
        UserComponentCollection::
        NominalProjectedTotal() const {
            TH1 * total = nullptr;
            for(const auto & component : fComponents) {
                if(!total) {
                    total = component.second->ProjectNominal();
                }
                else {
                    total->Add(component.second->ProjectNominal());
                }
            }
            return total;
        }

        TH1 *
        UserComponentCollection::
        NominalProjectedTotalForErrorCalculation() const {
            TH1 * total = nullptr;
            for(const auto & component : fComponents) {
                if(!total) {
                    total = component.second->ProjectNominalForErrorCalculation();
                }
                else {
                    total->Add(component.second->ProjectNominalForErrorCalculation());
                }
            }
            return total;
        }

        Systematic<TH1>
        UserComponentCollection::
        SystematicTotal(std::string syst_label) const {
            std::vector<std::shared_ptr<TH1>> total_shifted(
                    fComponents.begin()->second->GetSystematics().at(syst_label).GetShifts().size(), nullptr
            );
            for(auto const & component : fComponents) {
                auto component_shifts = component.second->GetSystematics().at(syst_label).GetShifts();
                for(auto i = 0u; i < total_shifted.size(); i++) {
                    if(!total_shifted[i]) {
                        total_shifted[i] = std::shared_ptr<TH1>((TH1*) component_shifts[i]->Clone());
                    }
                    else {
                        total_shifted[i]->Add(component_shifts[i].get());
                    }
                }
            }
            return Systematic<TH1>(fComponents.begin()->second->GetSystematics().at(syst_label).GetName(),
                                   total_shifted,
                                   fComponents.begin()->second->GetSystematics().at(syst_label).GetType());
        }

        Systematic<TH1>
        UserComponentCollection::
        SystematicProjectedTotal(std::string syst_label) const {
            std::vector<std::shared_ptr<TH1>> projected_total_shifted(
                    fComponents.begin()->second->GetSystematics().at(syst_label).GetShifts().size(), nullptr
            );
            for(auto const & component : fComponents) {
                auto projected_component_shifts = component.second->ProjectSystematic(syst_label).GetShifts();
                for(auto i = 0u; i < projected_total_shifted.size(); i++) {
                    if(!projected_total_shifted[i]) {
                        projected_total_shifted[i] = std::shared_ptr<TH1>((TH1*) projected_component_shifts[i]->Clone());
                    }
                    else {
                        projected_total_shifted[i]->Add(projected_component_shifts[i].get());
                    }
                }
            }
            return Systematic<TH1>(fComponents.begin()->second->GetSystematics().at(syst_label).GetName(),
                                   projected_total_shifted,
                                   fComponents.begin()->second->GetSystematics().at(syst_label).GetType());
        }

        std::map<std::string, Systematic<TH1>>
        UserComponentCollection::
        TotalSystematics() const {
            std::map<std::string, Systematic<TH1>> ret;
            for(const auto & syst : fComponents.begin()->second->GetSystematics()) {
                ret[syst.first] = this->SystematicTotal(syst.first);
            }
            return ret;
        }

        std::map<std::string, Systematic<TH1>>
        UserComponentCollection::
        ProjectedTotalSystematics() const {
            std::map<std::string, Systematic<TH1>> ret;
            for(const auto & syst : fComponents.begin()->second->GetSystematics()) {
                ret[syst.first] = this->SystematicProjectedTotal(syst.first);
            }
            return ret;
        }


        ReducedComponentCollection::
        ReducedComponentCollection(std::map<std::string, const IReducedTemplateComponent *> components)
                : fComponents(components),
                  fNInnerBins(components.begin()->second->GetNominal()->GetNInnerBins()),
                  fNOuterBins(components.begin()->second->GetNominal()->GetNOuterBins()) {
            // check all templates are consistent with one another
            bool consistent_templates = true;
            int i = 0;
            for (const auto & component : fComponents) {
                consistent_templates &= components.begin()->second->GetNominal()->GetArray().size() ==
                                        component.second->GetNominal()->GetArray().size();
                fComponentIdx[component.first] = i;
                i++;
            }
            assert(consistent_templates);
        }

        Vector
        ReducedComponentCollection::
        Predict(const Vector & user_params) const {
            Matrix prediction = Matrix::Zero(fNInnerBins * fNOuterBins,
                                             fComponents.size());
            auto user_params_mat = user_params.reshaped(fNOuterBins,
                                                        fComponents.size());
            int i = 0;
            for (const auto & component : fComponents) {
                prediction.col(i) = component.second->Predict(user_params_mat.col(i));
                i++;
            }
            return prediction.rowwise().sum();
        }

        Vector
        ReducedComponentCollection::
        PredictComponent(std::string component_label,
                         const Vector & component_params) const {
            return fComponents.at(component_label)->Predict(component_params);
        }

        ReducedComponent *
        ReducedComponentCollection::
        NominalTotal() const {
            Array total = Array::Zero(fComponents.begin()->second->GetNominal()->GetArray().size());
            for(auto component : fComponents) {
                total += component.second->GetNominal()->GetArray().array();
            }
            return new ReducedComponent(total.matrix(),
                                        fComponents.begin()->second->GetNominal()->GetNOuterBins(),
                                        fComponents.begin()->second->GetNominal()->GetNInnerBins());
        }

        TH1*
        IUserTemplateComponent::
        ProjectNominal() const {
            return ComponentReducer::Project(this->GetNominal());
        }

        TH1*
        IUserTemplateComponent::
        ProjectNominalForErrorCalculation() const {
            return ComponentReducer::Project(this->GetNominalForErrorCalculation());
        }

        std::map<std::string, Systematic<TH1>>
        IUserTemplateComponent::
        ProjectSystematics() const {
            std::map<std::string, Systematic<TH1>> projected_systs;
            for(const auto & syst : this->GetSystematics()) {
                projected_systs[syst.first] = ComponentReducer::Project(syst.second);
            }
            return projected_systs;
        }

        Systematic<TH1>
        IUserTemplateComponent::
        ProjectSystematic(std::string syst_label) const {
            return ComponentReducer::Project(this->GetSystematics().at(syst_label));
        }

        Array
        IReducedTemplateComponent::
        ProjectNominal() const {
            return this->GetNominal()->Project();
        }

        Array
        IReducedTemplateComponent::
        ProjectNominalForErrorCalculation() const {
            return this->GetNominalForErrorCalculation()->Project();
        }

        std::map<std::string, Systematic<TH1>>
        IReducedTemplateComponent::
        ProjectSystematics() const {
            std::map<std::string, Systematic<TH1>> ret;
            for(const auto & syst : this->GetSystematics()) {
                ret[syst.first] = this->ProjectSystematic(syst.first);
            }
            return ret;
        }

        Systematic<TH1>
        IReducedTemplateComponent::
        ProjectSystematic(std::string syst_label) const {
            auto project_like = new TH1D("", "",
                                         this->GetNominal()->GetNOuterBins()-2,
                                         0, this->GetNominal()->GetNOuterBins()-2);
            const Systematic<TH1> & systematic = this->GetSystematics().at(syst_label);
            std::vector<std::shared_ptr<TH1>> projected_shifts(systematic.GetShifts().size());
            for(auto i = 0u; i < systematic.GetShifts().size(); i++) {
                Array projection = root::MapContentsToEigenInner(systematic.GetShifts()[i].get())
                        .reshaped(this->GetNominal()->GetNInnerBins(), this->GetNominal()->GetNOuterBins())
                        .colwise().sum();
                projected_shifts[i] = std::shared_ptr<TH1>(root::ToROOTLike(project_like, projection));
            }
            return Systematic<TH1>(systematic.GetName(),
                                   projected_shifts,
                                   systematic.GetType());

        }

        std::vector<std::string>
        ReducedComponentCollection::
        GetSystematicLabels() const {
            std::vector<std::string> labels;
            for(const auto & syst : fComponents.begin()->second->GetSystematics()) {
                labels.push_back(syst.first);
            }
            return labels;
        }

        std::vector<std::string>
        UserComponentCollection::
        GetSystematicLabels() const {
            std::vector<std::string> labels;
            for(const auto & syst : fComponents.begin()->second->GetSystematics()) {
                labels.push_back(syst.first);
            }
            return labels;
        }
    }
}