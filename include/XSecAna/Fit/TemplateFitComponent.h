#pragma once
#include "XSecAna/Utils.h"
#include "XSecAna/Systematic.h"

namespace xsec {
    namespace fit {
        namespace detail {
            ///\brief Maps sparse arrays to/from dense arrays
            // eg. 0011100000101011 <-> 1111111
            class ParamMap {
            public:
                ParamMap() = default;
                explicit ParamMap(int size);
                explicit ParamMap(const Array & mask);
                unsigned int GetNMinimizerParams() const;
                unsigned int GetNUserParams() const;
                Vector ToUserParams(const Vector & minimizer_params) const;
                Vector ToMinimizerParams(const Vector & user_params) const;
                void MaskTemplate(int i);
                void UnmaskTemplate(int template_idx);
                bool IsParamMasked(int i) const;
                const Matrix & GetMatrix() const;
                unsigned int UserToMinimizerIdx(int user_idx) const;
                unsigned int MinimizerToUserIdx(int minimizer_idx) const;

            private:
                Matrix fM;
            };
        }

        class Reducer {
        public:
            Reducer(const TH1 * mask, bool padded)
                    : fMask(mask), fPadded(padded)
            {
                Array m = root::MapContentsToEigen(fMask);
                for(auto i = 0; i < m.size(); i++) {
                    m(i) = m(i) ? 1 : 0;
                }
                fMap = detail::ParamMap(m);

                if(fPadded) {
                    fReducedProps = root::TH1Props(new TH1D("", "",
                                                            fMap.GetNMinimizerParams(),
                                                            0,
                                                            fMap.GetNMinimizerParams()));
                }
                else {
                    fReducedProps = root::TH1Props(new TH1D("", "",
                                                            fMap.GetNMinimizerParams()-2,
                                                            0,
                                                            fMap.GetNMinimizerParams()-2));
                }
            }

            int GetReducedIndex(int xbin, int ybin, int zbin = 0) const {
                int gbin = fMask->GetBin(xbin, ybin, zbin);
                return fMap.UserToMinimizerIdx(gbin);
            }

            void GetExpandedIndex(int reduced_index, int & xbin, int & ybin, int & zbin) const {
                int expanded_gbin = fMap.MinimizerToUserIdx(reduced_index);
                fMask->GetBinXYZ(expanded_gbin, xbin, ybin, zbin);
            }

            TH1 * Reduce(const TH1 * expanded) const {
                if(fPadded) {
                    Array c(fMap.GetNMinimizerParams()+2);
                    Array e(fMap.GetNMinimizerParams()+2);
                    c(Eigen::seqN(1, fMap.GetNMinimizerParams())) =
                            fMap.ToMinimizerParams(root::MapContentsToEigen(expanded));
                    e(Eigen::seqN(1, fMap.GetNMinimizerParams())) =
                            fMap.ToMinimizerParams(root::MapErrorsToEigen(expanded));
                    return root::ToROOT(c, e, fReducedProps);
                }
                else {
                    Array r = fMap.ToMinimizerParams(root::MapContentsToEigen(expanded));
                    Array e = fMap.ToMinimizerParams(root::MapErrorsToEigen(expanded));
                    return root::ToROOT(r, fReducedProps);
                }
            }

            TH1 * Expand(const TH1 * reduced) const {
                if(fPadded) {
                    Array c = fMap.ToUserParams(root::MapContentsToEigenInner(reduced));
                    Array e = fMap.ToUserParams(root::MapErrorsToEigen(reduced, false));
                    return root::ToROOTLike(fMask, c, e);
                }
                else {
                    Array c = fMap.ToUserParams(root::MapContentsToEigen(reduced));
                    Array e = fMap.ToUserParams(root::MapErrorsToEigen(reduced));
                    return root::ToROOTLike(fMask, c, e);
                }
            }

	    const TH1 * GetMask() const { return fMask; }

        private:
            const TH1 * fMask;
            detail::ParamMap fMap;
            bool fPadded;
            root::TH1Props fReducedProps;
        };

        ///\brief representation of reduced template component
        // containing both ROOT histogram and eigen array for easy access
        // and to prevent copies/transforms during the fit
        // "Reduced" means empty bins that are not to be included in the fit have been removed
        class ReducedComponent {
        public:
            //explicit ReducedComponent(const TH1 * h);
            ReducedComponent(const Vector & a, int nouter_bins, int ninner_bins);
            ReducedComponent(const TH1 * h, int nouter_bins, int ninner_bins);
            std::shared_ptr<TH1> GetHist() const { return fHist; }
            //[[nodiscard]] TH1 * ToHist1D() const;
            const Vector & GetArray() const { return fArray; };
            int GetNInnerBins() const { return fNInnerBins; }
            int GetNOuterBins() const { return fNOuterBins; }
            Array Project() const;

        private:
            std::shared_ptr<TH1> fHist;
            const Vector fArray;
            int fNOuterBins;
            int fNInnerBins;

        };


        ///\brief Object for transforming user-level components to fitter level
        // by removing extraneous bins determined from the user-provided mask
        class ComponentReducer {
        public:
            explicit ComponentReducer(const TH1 * mask);
            [[nodiscard]] ReducedComponent * Reduce(const std::shared_ptr<TH1> component) const;
            [[nodiscard]] Systematic<TH1> Reduce(const Systematic<TH1> & syst) const;
            [[nodiscard]] TH1 * Compress1D(const std::shared_ptr<TH1> component) const;
            static TH1 * Project(const std::shared_ptr<TH1> templ);
            static Systematic<TH1> Project(const Systematic<TH1> & syst);

            const TH1 * GetMask() const { return fMask; }
            const detail::ParamMap & GetMap() const { return fMap; }
        private:
            const TH1 * fMask;
            const detail::ParamMap fMap;
        };

        ///\brief The fitter level representation of a template fit component
        // Returned predictions are reweighted by an array of weights corresponding to
        // outer bins, or analysis bins.
        class IReducedTemplateComponent {
        public:
            [[nodiscard]] virtual Vector Predict(const Vector & component_params) const = 0;
            [[nodiscard]] virtual Vector PredictProjected(const Vector & component_params) const = 0;
            virtual const ReducedComponent * GetNominal() const = 0;
            virtual const ReducedComponent * GetNominalForErrorCalculation() const { return this->GetNominal(); }
            virtual const std::map<std::string, Systematic<TH1>> & GetSystematics() const = 0;
            virtual Array ProjectNominal() const final;
            virtual Array ProjectNominalForErrorCalculation() const final;
            virtual std::map<std::string, Systematic<TH1>> ProjectSystematics() const final;
            virtual Systematic<TH1> ProjectSystematic(std::string syst_label) const final;
        };

        ///\brief The user level representation of a template fit component
        // given to template fitters
        class IUserTemplateComponent {
        public:
            [[nodiscard]] virtual IReducedTemplateComponent * Reduce(const ComponentReducer & reducer) const = 0;
            virtual std::shared_ptr<TH1> GetNominal() const = 0;
            virtual std::shared_ptr<TH1> GetNominalForErrorCalculation() const { return this->GetNominal(); }
            virtual const std::map<std::string, Systematic<TH1>> & GetSystematics() const = 0;
            virtual TH1 * ProjectNominal() const final;
            virtual TH1 * ProjectNominalForErrorCalculation() const final;
            virtual std::map<std::string, Systematic<TH1>> ProjectSystematics() const final;
            virtual Systematic<TH1> ProjectSystematic(std::string syst_label) const final;
        };

        ///\brief a collection of components that when summed together, give the total prediction
        class ReducedComponentCollection {
        public:
            ReducedComponentCollection() = default;
            explicit ReducedComponentCollection(std::map<std::string, const IReducedTemplateComponent *> components);

            [[nodiscard]] Vector Predict(const Vector & user_params) const;
            [[nodiscard]] Vector PredictComponent(std::string component_label, const Vector & user_params) const;
            int GetNOuterBins() const { return fComponents.begin()->second->GetNominal()->GetNOuterBins(); }
            int GetNInnerBins() const { return fComponents.begin()->second->GetNominal()->GetNInnerBins(); }
            size_t size() const { return fComponents.size(); }
            const std::map<std::string, const IReducedTemplateComponent*> & GetComponents() const { return fComponents; }
            const IReducedTemplateComponent * GetComponent(std::string label) const { return fComponents.at(label); }
            int GetComponentIdx(const std::string & component_label) const { return fComponentIdx.at(component_label); }
            [[nodiscard]] std::vector<std::string> GetSystematicLabels() const;
            ReducedComponent * NominalTotal() const;

        private:
            int fNInnerBins;
            int fNOuterBins;
            std::map<std::string, const IReducedTemplateComponent *> fComponents;
            std::map<std::string, int> fComponentIdx;
        };

        class UserComponentCollection {
        public:
            UserComponentCollection() = default;
            explicit UserComponentCollection(std::map<std::string, const IUserTemplateComponent*> components)
                    : fComponents(components) {}

            const std::map<std::string, const IUserTemplateComponent*> & GetComponents() const { return fComponents; }
            std::map<std::string, const IUserTemplateComponent*> & GetComponents() { return fComponents; }
            [[nodiscard]] ReducedComponentCollection Reduce(const ComponentReducer & reducer) const;
            [[nodiscard]] TH1 * NominalProjectedTotal() const;
            [[nodiscard]] TH1 * NominalTotal() const;
            [[nodiscard]] TH1 * NominalProjectedTotalForErrorCalculation() const;
            [[nodiscard]] TH1 * NominalTotalForErrorCalculation() const;
            [[nodiscard]] Systematic<TH1> SystematicProjectedTotal(std::string syst_label) const;
            [[nodiscard]] Systematic<TH1> SystematicTotal(std::string syst_label) const;
            [[nodiscard]] std::map<std::string, Systematic<TH1>> ProjectedTotalSystematics() const;
            [[nodiscard]] std::map<std::string, Systematic<TH1>> TotalSystematics() const;
            [[nodiscard]] std::vector<std::string> GetSystematicLabels() const;
            const IUserTemplateComponent * GetComponent(std::string label) const { return fComponents.at(label); }
            //[[nodiscard]] TH1 * Predict(std::map<std::string, TH1 *> params) const;
            //[[nodiscard]] TH1 * PredictComponent(std::string component_label, const TH1 * user_params) const;

        private:
            std::map<std::string, const IUserTemplateComponent*> fComponents;
        };

        ///\brief Basic user-level template component object for single-sample template fitting
        class UserTemplateComponent : public IUserTemplateComponent {
        public:
            explicit UserTemplateComponent(const std::shared_ptr<TH1> mean,
                                           std::map<std::string, Systematic<TH1>> systematics = std::map<std::string, Systematic<TH1>>())
                    : fMean(mean), fSystematics(systematics) {}

            [[nodiscard]] IReducedTemplateComponent * Reduce(const ComponentReducer & reducer) const override;
            std::shared_ptr<TH1> GetNominal() const override { return fMean; }
            const std::map<std::string, Systematic<TH1>> & GetSystematics() const override { return fSystematics; }

        private:
            const std::shared_ptr<TH1> fMean;
            const std::map<std::string, Systematic<TH1>> fSystematics;
        };

        class UserTemplateComponentAlternativeNominal : public UserTemplateComponent {
        public:
            explicit UserTemplateComponentAlternativeNominal(const std::shared_ptr<TH1> mean,
                                                             std::map<std::string, Systematic<TH1>> systematics = std::map<std::string, Systematic<TH1>>(),
                                                             const std::shared_ptr<TH1> mean_for_error_calc = 0)
            : UserTemplateComponent(mean, systematics),
              fMeanForErrorCalc(mean_for_error_calc) {}

            virtual std::shared_ptr<TH1> GetNominalForErrorCalculation() const override { return fMeanForErrorCalc; }
            [[nodiscard]] IReducedTemplateComponent * Reduce(const ComponentReducer & reducer) const override;

        private:
            std::shared_ptr<TH1> fMeanForErrorCalc;
        };

        ///\brief Basic fitter-level template component object for single-sample template fitting
        class ReducedTemplateComponent : public IReducedTemplateComponent {
        public:
            explicit ReducedTemplateComponent(const ReducedComponent * mean,
                                              const std::map<std::string, Systematic<TH1>> systematics = std::map<std::string, Systematic<TH1>>())
                    : fMean(mean), fSystematics(systematics) {}

            [[nodiscard]] Vector Predict(const Vector & component_params) const override;
            [[nodiscard]] Vector PredictProjected(const Vector & component) const override;
            const ReducedComponent * GetNominal() const override { return fMean; }
            const std::map<std::string, Systematic<TH1>> & GetSystematics() const override { return fSystematics; }

        private:
            const ReducedComponent * fMean;
            const std::map<std::string, Systematic<TH1>> fSystematics;
        };

        class ReducedTemplateComponentAlternativeNominal : public ReducedTemplateComponent {
        public:
            ReducedTemplateComponentAlternativeNominal(const ReducedComponent * mean,
                                                       const std::map<std::string, Systematic<TH1>> systematics = std::map<std::string, Systematic<TH1>>(),
                                                       const ReducedComponent * mean_for_error_calc = 0)
                    : ReducedTemplateComponent(mean, systematics),
                      fMeanForErrorCalc(mean_for_error_calc) {}
            const ReducedComponent * GetNominalForErrorCalculation() const override { return fMeanForErrorCalc; }
        private:
            const ReducedComponent * fMeanForErrorCalc;
        };

        struct TemplateFitSample {
            TemplateFitSample() = default;

            TemplateFitSample(std::map<std::string, const IUserTemplateComponent*> & _components,
                              std::map<std::string, Systematic<TH1>> _shape_only_systematics)
                    : components(_components), shape_only_systematics(_shape_only_systematics) {}

            UserComponentCollection components;
            std::map<std::string, Systematic<TH1>> shape_only_systematics;
        };
    }
}
