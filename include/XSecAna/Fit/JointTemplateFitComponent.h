#pragma once

namespace xsec {
    namespace fit {
        namespace detail {
            // Functions for joining template components from multiple samples
            // All of the work is in here
            std::shared_ptr<TH1> _join(const std::shared_ptr<TH1> lhs, const std::shared_ptr<TH1> rhs);
            std::shared_ptr<TH1> _join(const std::vector<std::shared_ptr<TH1>> & samples);

            // various overloads for more complicated objects
            Systematic <TH1> _join(const Systematic <TH1> & lhs, const Systematic <TH1> & rhs);
            Systematic <TH1> _join(const std::vector <Systematic<TH1>> & samples);
            Systematic <TH1> _join(const std::vector <Systematic<TH1>> & samples);

            TemplateFitSample _join(const std::map<std::string, TemplateFitSample> & components);


            // templated function to join maps of vectors
            // requires an implementation of _join(const std::vector<T>&)->T
            template<class T>
            std::map <std::string, T>
            _join(std::vector <std::map<std::string, const T &>> samples) {
                std::map <std::string, T> joined;
                for (auto tmap: samples[0]) {
                    std::vector <T> tmp;
                    for (auto i = 0u; i < samples.size(); i++) {
                        tmp.push_back(samples[0].at(tmap.first));
                    }
                    joined[tmap.first] = _join(tmp);
                }
                return joined;
            }
        }
        /*
        class UserSingleSampleTemplateComponent : public IUserTemplateComponent {
        public:
            UserSingleSampleTemplateComponent(const TH1 * nominal,
                                              const std::map<std::string, xsec::Systematic<TH1>> & systematics)
                    : fNominal(nominal), fSystematics(systematics) {}

            const TH1 * GetNominal() const override { return fNominal; }
            const std::map<std::string, Systematic<TH1>> & GetSystematics() const { return fSystematics; }
            IReducedTemplateComponent * Reduce(const ComponentReducer & reducer) const override;
            TH1 * ProjectNominal() const;
            std::map<std::string, Systematic<TH1>> ProjectSystematics() const;

        private:
            const TH1 * fNominal;
            std::map<std::string, Systematic<TH1>> fSystematics;
        };


        class ReducedSingleSampleTemplateComponent : public ReducedTemplateComponent {
        public:
            ReducedSingleSampleTemplateComponent(const ReducedComponent * mean,
                                                 std::map<std::string, Systematic<TH1>> & systematics)
            : ReducedTemplateComponent(mean),  fSystematics(systematics) {}

            const std::map<std::string, Systematic<TH1>> & GetSystematics() const { return fSystematics; }
        private:
            const std::map<std::string, Systematic<TH1>> & fSystematics;
        };
*/
        ///\brief User-level template fit component for simultaneous fitting of
        // multiple correlated samples
        class UserJointTemplateComponent : public IUserTemplateComponent {
        public:
            UserJointTemplateComponent(const std::map<std::string, const UserTemplateComponent *> & sample_components);
            [[nodiscard]] IReducedTemplateComponent * Reduce(const ComponentReducer & reducer) const override;

            std::shared_ptr<TH1> GetNominal() const override { return fJointTemplateMean; }

            const std::map<std::string, Systematic<TH1>> & GetSystematics() const override
            { return fJointTemplateSystematics; }

            const std::shared_ptr<TH1> GetSampleNormalization(const std::string & sample_name) const
            { return fSampleNormalizations.at(sample_name); }

            const std::shared_ptr<TH1> GetJointNormalization() const
            { return fJointNormalization; }

            const std::map<std::string, const UserTemplateComponent *> & GetUserTemplateComponents() const
            { return fSamples; }

        private:
            static std::shared_ptr<TH1> _to1d(const std::shared_ptr<TH1> h);
            static std::shared_ptr<TH1> _join1d(const std::shared_ptr<TH1> lhs, const std::shared_ptr<TH1> rhs);
            static std::shared_ptr<TH1> _join1d(std::vector<std::shared_ptr<TH1>> hists);
            static Systematic<TH1> _to1d(const Systematic<TH1> & syst);

            const std::map<std::string, const UserTemplateComponent*> fUserTemplateComponents;
            std::shared_ptr<TH1> fJointTemplateMean;
            std::map<std::string, Systematic<TH1>> fJointTemplateSystematics;

            std::shared_ptr<TH1> fJointNormalization;
            std::map<std::string, std::shared_ptr<TH1>> fSampleNormalizations;
            std::map<std::string, std::map<std::string, Systematic<TH1>>> fSampleNormSystematics;
            std::map<std::string, Systematic<TH1>> fJointNormSystematics;

            TH1 * fJointNormalizationTotalCovariance;
            std::map<std::string, TH1*> fJointNormalizationCovariances;
            std::map<std::string, TH1*> fSampleNormalizationTotalCovariance;
            std::map<std::string, std::map<std::string, TH1*>> fSampleNormalizationCovariances;

            std::string fConditioningSampleLabel;
            std::map<std::string, const UserTemplateComponent*> fSamples;

        };


        ///\brief Fitter-level template fit component for simultaneous fitting of multiple
        // correlated samples
        class ReducedJointTemplateComponent : public IReducedTemplateComponent {
        public:
            ReducedJointTemplateComponent(std::map<std::string, const ReducedTemplateComponent*> samples,
                                          const ReducedComponent * joined_template_mean,
                                          std::shared_ptr<TH1> joined_normalization,
                                          std::map<std::string, std::shared_ptr<TH1>> & sample_normalizations,
                                          std::map<std::string, Systematic<TH1>> & joint_norm_systematics);

            [[nodiscard]] Vector Predict(const Vector & component_params) const override;
            [[nodiscard]] Vector ComplimentaryParams(const Vector & conditioning_params) const;

            const ReducedComponent * GetNominal() const override
            { return fJointTemplateMean; }

            const std::map<std::string, Systematic<TH1>> & GetSystematics() const override
            { return fJointTemplateSystematics; }

            TH1 * GetTotalJointNormalizationCovariance() const
            { return fJoinedNormTotalCovariance; }

            TH1 * GetJointNormalizationCovariance(std::string systematic_label) const
            { return fJoinedNormCovariances.at(systematic_label); }

            const std::shared_ptr<TH1> GetJointNormalization() const
            { return fJointNormalization; }

            const std::shared_ptr<TH1> GetSampleNormalization(std::string sample_name) const
            { return fSampleNormalizations.at(sample_name); }

            const std::map<std::string, Systematic<TH1>> & GetJointNormalizationSystematics() const
            { return fJointNormSystematics; }

            const std::shared_ptr<TH1> GetConditioningSampleNormalization() const
            { return fSampleNormalizations.at(fConditioningSampleLabel); }

            const std::shared_ptr<TH1> GetComplimentarySampleNormalization() const
            { return fSampleNormalizations.at(fComplimentarySampleLabel); }

            const ReducedTemplateComponent * GetConditioningComponent() const
            { return fConditioningSample; }

            const ReducedTemplateComponent * GetComplimentaryComponent() const
            { return fComplimentarySample; }

            bool IsSingular() const { return fIsSingular; }

        private:
            bool _is_singular();
            const ReducedComponent * fJointTemplateMean;
            std::map<std::string, Systematic<TH1>> fJointTemplateSystematics;

            std::shared_ptr<TH1> fJointNormalization;
            const std::map<std::string, std::shared_ptr<TH1>> fSampleNormalizations;
            const std::map<std::string, Systematic<TH1>> fJointNormSystematics;

            TH1 * fJoinedNormTotalCovariance;
            std::map<std::string, TH1*> fJoinedNormCovariances;
            Matrix fConditioningSampleInvCovariance;
            Matrix fCrossSampleCovariance;

            Array fConditioningSampleNormalization;
            Vector fComplimentarySampleNormalization;

            std::string fConditioningSampleLabel;
            std::string fComplimentarySampleLabel;
            int fNOuterBins;

            const ReducedTemplateComponent * fConditioningSample;
            const ReducedTemplateComponent * fComplimentarySample;
            bool fIsSingular;
        };
    }
}
