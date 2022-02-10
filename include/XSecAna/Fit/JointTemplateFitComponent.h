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

        ///\brief User-level template fit component for simultaneous fitting of
        // multiple correlated samples
        class UserJointTemplateComponent : public IUserTemplateComponent {
        public:
            UserJointTemplateComponent(const std::map<std::string, const IUserTemplateComponent *> & sample_components);
            [[nodiscard]] IReducedTemplateComponent * Reduce(const ComponentReducer & reducer) const override;

            std::shared_ptr<TH1> GetNominal() const override { return fJointTemplateMean; }

            const std::map<std::string, Systematic<TH1>> & GetSystematics() const override
            { return fJointTemplateSystematics; }

            const std::shared_ptr<TH1> GetSampleNormalization(const std::string & sample_name) const
            { return fSampleNormalizations.at(sample_name); }

            const std::shared_ptr<TH1> GetConditioningSampleNormalization() const
            { return fSampleNormalizations.at(fConditioningSampleLabel); }

            const std::shared_ptr<TH1> GetComplimentarySampleNormalization() const
            { return fSampleNormalizations.at(fComplimentarySampleLabel); }

            const std::shared_ptr<TH1> GetJointNormalization() const
            { return fJointNormalization; }

            const std::map<std::string, const IUserTemplateComponent *> & GetUserTemplateComponents() const
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
            std::shared_ptr<TH1> fJointNormalizationForErrorCalc;
            std::map<std::string, std::shared_ptr<TH1>> fSampleNormalizations;
            std::map<std::string, std::shared_ptr<TH1>> fSampleNormalizationsForErrorCalc;
            std::map<std::string, std::map<std::string, Systematic<TH1>>> fSampleNormSystematics;
            std::map<std::string, Systematic<TH1>> fJointNormSystematics;

            TH1 * fJointNormalizationTotalCovariance;
            std::map<std::string, TH1*> fJointNormalizationCovariances;
            std::map<std::string, TH1*> fSampleNormalizationTotalCovariance;
            std::map<std::string, std::map<std::string, TH1*>> fSampleNormalizationCovariances;


            std::string fConditioningSampleLabel;
            std::string fComplimentarySampleLabel;
            std::map<std::string, const IUserTemplateComponent*> fSamples;

        };


        ///\brief Fitter-level template fit component for simultaneous fitting of multiple
        // correlated samples
        class ReducedJointTemplateComponent : public IReducedTemplateComponent {
        public:
            ReducedJointTemplateComponent(std::map<std::string, const IReducedTemplateComponent*> samples,
                                          const ReducedComponent * joined_template_mean,
                                          std::shared_ptr<TH1> joined_normalization,
                                          std::shared_ptr<TH1> joined_normalization_for_error_calc,
                                          std::map<std::string, std::shared_ptr<TH1>> & sample_normalizations,
                                          std::map<std::string, std::shared_ptr<TH1>> & sample_normalizations_for_error_calc,
                                          std::map<std::string, Systematic<TH1>> & joint_norm_systematics);

            [[nodiscard]] Vector Predict(const Vector & condi_params) const override;
            [[nodiscard]] Vector PredictProjected(const Vector & condi_params) const override;
            [[nodiscard]] Vector ComplimentaryParams(const Vector & conditioning_params) const;
            [[nodiscard]] Vector ConditionalParams(const Vector & complimentary_params) const;

            const ReducedComponent * GetNominal() const override
            { return fJointTemplateMean; }

            const std::map<std::string, Systematic<TH1>> & GetSystematics() const override
            { return fJointTemplateSystematics; }

            TH1 * GetJointNormalizationTotalCovariance() const
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

            const IReducedTemplateComponent * GetConditioningComponent() const
            { return fConditioningSample; }

            const IReducedTemplateComponent * GetComplimentaryComponent() const
            { return fComplimentarySample; }

            bool IsSingular() const { return fIsSingular; }

        private:
            bool _is_singular();
            const ReducedComponent * fJointTemplateMean;
            std::map<std::string, Systematic<TH1>> fJointTemplateSystematics;

            std::shared_ptr<TH1> fJointNormalization;
            std::shared_ptr<TH1> fJointNormalizationForErrorCalc;
            const std::map<std::string, std::shared_ptr<TH1>> fSampleNormalizations;
            const std::map<std::string, std::shared_ptr<TH1>> fSampleNormalizationsForErrorCalc;
            const std::map<std::string, Systematic<TH1>> fJointNormSystematics;

            TH1 * fJoinedNormTotalCovariance;
            std::map<std::string, TH1*> fJoinedNormCovariances;
            Matrix fConditioningSampleInvCovariance;
            Matrix fCrossSampleCovariance;
            Matrix fRotationMatrix;

            Array fConditioningSampleNormalization;
            Vector fComplimentarySampleNormalization;

            std::string fConditioningSampleLabel;
            std::string fComplimentarySampleLabel;
            int fNOuterBins;

            const IReducedTemplateComponent * fConditioningSample;
            const IReducedTemplateComponent * fComplimentarySample;
            bool fIsSingular;
        };
    }
}
