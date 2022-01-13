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

            private:
                Matrix fM;
            };
        }

        ///\brief representation of reduced template component
        // containing both ROOT histogram and eigen array for easy access
        // and to prevent copies/transforms during the fit
        // "Reduced" means empty bins that are not to be included in the fit have been removed
        class ReducedComponent {
        public:
            //explicit ReducedComponent(const TH1 * h);
            ReducedComponent(const Vector & a, int nouter_bins, int ninner_bins);
            ReducedComponent(const TH1 * h, int nouter_bins, int ninner_bins);
            const TH1 * GetHist() const { return fHist; }
            [[nodiscard]] TH1 * ToHist1D() const;
            const Vector & GetArray() const { return fArray; };
            int GetNInnerBins() const { return fNInnerBins; }
            int GetNOuterBins() const { return fNOuterBins; }

        private:
            const TH1 * fHist;
            const Vector fArray;
            int fNOuterBins;
            int fNInnerBins;

        };


        ///\brief Object for transforming user-level components to fitter level
        // by removing extraneous bins determined from the user-provided mask
        class ComponentReducer {
        public:
            explicit ComponentReducer(const TH1 * mask);
            [[nodiscard]] ReducedComponent * Reduce(const TH1 * component) const;
            [[nodiscard]] Systematic<TH1> Reduce(const Systematic<TH1> & syst) const;
            [[nodiscard]] TH1 * Compress1D(const TH1 * component) const;
            static TH1 * Project(const TH1 * templ);
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
            virtual const ReducedComponent * GetNominal() const = 0;
            virtual const std::map<std::string, Systematic<TH1>> & GetSystematics() const = 0;
        };

        ///\brief The user level representation of a template fit component
        // given to template fitters
        class IUserTemplateComponent {
        public:
            [[nodiscard]] virtual IReducedTemplateComponent * Reduce(const ComponentReducer & reducer) const = 0;
            virtual const TH1 * GetNominal() const = 0;
            virtual const std::map<std::string, Systematic<TH1>> & GetSystematics() const = 0;
            virtual TH1 * ProjectNominal() const final;
            virtual std::map<std::string, Systematic<TH1>> ProjectSystematics() const final;
        };

        ///\brief a collection of components that when summed together, give the total prediction
        class ReducedComponentCollection {
        public:
            ReducedComponentCollection() = default;
            explicit ReducedComponentCollection(std::map<std::string, IReducedTemplateComponent *> components);

            [[nodiscard]] Vector Predict(const Vector & user_params) const;
            [[nodiscard]] Vector PredictComponent(std::string component_label, const Vector & user_params) const;
            int GetNOuterBins() const { return fComponents.begin()->second->GetNominal()->GetNOuterBins(); }
            int GetNInnerBins() const { return fComponents.begin()->second->GetNominal()->GetNInnerBins(); }
            size_t size() const { return fComponents.size(); }
            const std::map<std::string, IReducedTemplateComponent*> & GetComponents() const { return fComponents; }
            int GetComponentIdx(const std::string & component_label) const { return fComponentIdx.at(component_label); }

        private:
            int fNInnerBins;
            int fNOuterBins;
            std::map<std::string, IReducedTemplateComponent *> fComponents;
            std::map<std::string, int> fComponentIdx;
        };

        class UserComponentCollection {
        public:
            UserComponentCollection() = default;
            explicit UserComponentCollection(std::map<std::string, IUserTemplateComponent*> components)
                    : fComponents(components) {}

            const std::map<std::string, IUserTemplateComponent*> & GetComponents() const { return fComponents; }
            ReducedComponentCollection Reduce(const ComponentReducer & reducer) const;
            [[nodiscard]] TH1 * Predict(std::map<std::string, TH1 *> params) const;
            [[nodiscard]] TH1 * PredictComponent(std::string component_label, const TH1 * user_params) const;

        private:
            std::map<std::string, IUserTemplateComponent*> fComponents;
        };

        ///\brief Basic user-level template component object for single-sample template fitting
        class UserTemplateComponent : public IUserTemplateComponent {
        public:
            explicit UserTemplateComponent(const TH1 * mean,
                                           const std::map<std::string, Systematic<TH1>> & systematics = std::map<std::string, Systematic<TH1>>())
                    : fMean(mean), fSystematics(std::move(systematics)) {}

            [[nodiscard]] IReducedTemplateComponent * Reduce(const ComponentReducer & reducer) const override;
            const TH1 * GetNominal() const override { return fMean; }
            const std::map<std::string, Systematic<TH1>> & GetSystematics() const override { return fSystematics; }

        private:
            const TH1 * fMean;
            const std::map<std::string, Systematic<TH1>> & fSystematics;
        };

        ///\brief Basic fitter-level template component object for single-sample template fitting
        class ReducedTemplateComponent : public IReducedTemplateComponent {
        public:
            explicit ReducedTemplateComponent(const ReducedComponent * mean,
                                              const std::map<std::string, Systematic<TH1>> & systematics = std::map<std::string, Systematic<TH1>>())
                    : fMean(mean), fSystematics(std::move(systematics)) {}

            [[nodiscard]] Vector Predict(const Vector & component_params) const override;
            const ReducedComponent * GetNominal() const override { return fMean; }
            const std::map<std::string, Systematic<TH1>> & GetSystematics() const override { return fSystematics; }

        private:
            const ReducedComponent * fMean;
            const std::map<std::string, Systematic<TH1>> & fSystematics;
        };

        struct TemplateFitSample {
            TemplateFitSample() = default;

            TemplateFitSample(std::map<std::string, IUserTemplateComponent*> & _components,
                              std::map<std::string, Systematic<TH1>> _shape_only_systematics)
                    : components(_components), shape_only_systematics(_shape_only_systematics) {}

            UserComponentCollection components;
            std::map<std::string, Systematic<TH1>> shape_only_systematics;
        };
    }
}

