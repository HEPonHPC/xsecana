#include "XSecAna/Fit/TemplateFitComponent.h"
#include "XSecAna/Fit/JointTemplateFitComponent.h"
#include "XSecAna/Utils.h"
#include "XSecAna/Systematic.h"

#include <functional>

namespace xsec {
    namespace fit {
        namespace detail {
            std::shared_ptr<TH1>
            _join(const std::shared_ptr<TH1> lhs,
                  const std::shared_ptr<TH1> rhs) {
                std::shared_ptr<TH1> ret;
                if (lhs->GetDimension() == 1) {
                    int njoined = lhs->GetNbinsX() + rhs->GetNbinsX();
                    ret = std::make_shared<TH1D>("", "",
                                                 njoined, 0, njoined);
                    for (int x = 1; x <= lhs->GetNbinsX(); x++) {
                        ret->SetBinContent(x, lhs->GetBinContent(x));
                        ret->SetBinError(x, lhs->GetBinError(x));
                    }
                    for (int x = 1; x <= rhs->GetNbinsX(); x++) {
                        ret->SetBinContent(x + lhs->GetNbinsX(), rhs->GetBinContent(x));
                        ret->SetBinError(x + lhs->GetNbinsX(), rhs->GetBinError(x));
                    }
                } else if (lhs->GetDimension() == 2) {
                    int njoined = lhs->GetNbinsY() + rhs->GetNbinsY();
                    if (lhs->GetXaxis()->IsVariableBinSize()) {
                        std::vector<double> ybins(njoined + 1);
                        for (int i = 0; i < njoined + 1; i++) ybins[i] = i;
                        ret = std::make_shared<TH2D>("", "",
                                                     lhs->GetXaxis()->GetNbins(),
                                                     lhs->GetXaxis()->GetXbins()->GetArray(),
                                                     njoined,
                                                     &ybins[0]);
                    } else {
                        ret = std::make_shared<TH2D>("", "",
                                                     lhs->GetNbinsX(),
                                                     lhs->GetXaxis()->GetBinLowEdge(1),
                                                     lhs->GetXaxis()->GetBinLowEdge(lhs->GetNbinsX() + 1),
                                                     njoined,
                                                     0,
                                                     njoined);
                    }
                    ret->GetXaxis()->SetTitle(lhs->GetXaxis()->GetTitle());
                    ret->GetYaxis()->SetTitle("Joined Template Bins");
                    for (auto x = 1; x <= lhs->GetNbinsX(); x++) {
                        for (auto y = 1; y <= lhs->GetNbinsY(); y++) {
                            ret->SetBinContent(x,
                                               y,
                                               lhs->GetBinContent(x, y));
                            ret->SetBinError(x,
                                             y,
                                             lhs->GetBinError(x, y));
                        }
                        for (auto y = 1; y <= rhs->GetNbinsY(); y++) {
                            ret->SetBinContent(x,
                                               y + lhs->GetNbinsY(),
                                               rhs->GetBinContent(x, y));
                            ret->SetBinError(x,
                                             y + lhs->GetNbinsY(),
                                             rhs->GetBinError(x, y));
                        }
                    }
                } else {
                    int njoined = lhs->GetNbinsZ() + rhs->GetNbinsZ();
                    if (lhs->GetXaxis()->IsVariableBinSize() ||
                        lhs->GetYaxis()->IsVariableBinSize()) {
                        std::vector<double> zbins(njoined + 1);
                        for (int i = 0; i < njoined + 1; i++) zbins[i] = i;
                        ret = std::make_shared<TH3D>("", "",
                                                     lhs->GetNbinsX(),
                                                     lhs->GetXaxis()->GetXbins()->GetArray(),
                                                     lhs->GetNbinsY(),
                                                     lhs->GetYaxis()->GetXbins()->GetArray(),
                                                     njoined,
                                                     &zbins[0]);
                    } else {
                        ret = std::make_shared<TH3D>("", "",
                                                     lhs->GetNbinsX(),
                                                     lhs->GetXaxis()->GetBinLowEdge(1),
                                                     lhs->GetXaxis()->GetBinLowEdge(lhs->GetNbinsX() + 1),
                                                     lhs->GetNbinsY(),
                                                     lhs->GetYaxis()->GetBinLowEdge(1),
                                                     lhs->GetYaxis()->GetBinLowEdge(lhs->GetNbinsY() + 1),
                                                     njoined,
                                                     0,
                                                     njoined);
                    }
                    for (int x = 1; x <= lhs->GetNbinsX(); x++) {
                        for (int y = 1; y <= lhs->GetNbinsY(); y++) {
                            for (auto z = 1; z <= lhs->GetNbinsZ(); z++) {
                                ret->SetBinContent(x,
                                                   y,
                                                   z,
                                                   lhs->GetBinContent(x, y, z));
                                ret->SetBinError(x,
                                                 y,
                                                 z,
                                                 lhs->GetBinError(x, y, z));
                            }
                            for (auto z = 1; z <= rhs->GetNbinsZ(); z++) {
                                ret->SetBinContent(x,
                                                   y,
                                                   z + rhs->GetNbinsZ(),
                                                   rhs->GetBinContent(x, y, z));
                                ret->SetBinError(x,
                                                 y,
                                                 z + rhs->GetNbinsZ(),
                                                 rhs->GetBinError(x, y, z));
                            }
                        }
                    }
                    ret->GetXaxis()->SetTitle(lhs->GetXaxis()->GetTitle());
                    ret->GetYaxis()->SetTitle(lhs->GetYaxis()->GetTitle());
                    ret->GetZaxis()->SetTitle("Joined Template Bins");
                }
                return ret;
            }

            TemplateFitSample
            _join(const std::map<std::string, TemplateFitSample> & samples) {
                std::map<std::string, const IUserTemplateComponent *> joined_components;
                // join components
                for (const auto & component: samples.begin()->second.components.GetComponents()) {
                    std::map<std::string, const IUserTemplateComponent*> tmp_comp;
                    for (const auto & sample: samples) {
                        tmp_comp[sample.first] =
                                (UserTemplateComponent *) samples.at(sample.first)
                                        .components.GetComponents().at(component.first);
                    }
                    joined_components[component.first] = new UserJointTemplateComponent(tmp_comp);
                }
                // join shape only systematics
                std::map<std::string, Systematic<TH1>> joined_shape_only_systematics;
                for (const auto & systematic : samples.begin()->second.shape_only_systematics) {
                    std::vector<Systematic<TH1>> tmp_systs;
                    for (const auto & sample: samples) {
                        tmp_systs.push_back(samples.at(sample.first).shape_only_systematics.at(systematic.first));
                    }
                    joined_shape_only_systematics[systematic.first] = _join(tmp_systs);
                }
                return {joined_components, joined_shape_only_systematics};
            }

            Systematic<TH1>
            _join(const Systematic<TH1> & lhs,
                  const Systematic<TH1> & rhs) {
                std::vector<std::shared_ptr<TH1>> joined;
                for (auto i = 0u; i < lhs.GetShifts().size(); i++) {
                    joined.push_back(_join(lhs.GetShifts()[i], rhs.GetShifts()[i]));
                }
                return Systematic<TH1>(lhs.GetName(), joined, lhs.GetType());
            }

            std::shared_ptr<TH1>
            _join(const std::vector<std::shared_ptr<TH1>> & samples) {
                std::shared_ptr<TH1> ret = nullptr;
                for (auto sample: samples) {
                    if (!ret) ret = std::shared_ptr<TH1>((TH1 *) sample->Clone());
                    else {
                        ret = _join(ret, sample);
                    }
                }
                return ret;
            }

            Systematic<TH1>
            _join(const std::vector<Systematic<TH1>> & samples) {
                auto ret = samples[0];
                for (auto i = 1u; i < samples.size(); i++) {
                    ret = _join(ret, samples[i]);
                }
                return ret;
            }
        }

        /*
        TH1 *
        UserSingleSampleTemplateComponent::
        ProjectNominal() const {
            return ComponentReducer::Project(fNominal);
        }

        std::map<std::string, Systematic<TH1>>
        UserSingleSampleTemplateComponent::
        ProjectSystematics() const {
            std::map<std::string, Systematic<TH1>> projected_systematics;
            for(const auto & syst : fSystematics) {
                projected_systematics[syst.first] = ComponentReducer::Project(syst.second);
            }
            return projected_systematics;
        }
       */

        std::shared_ptr<TH1>
        UserJointTemplateComponent::
        _to1d(const std::shared_ptr<TH1> h) {
            if (h->GetDimension() == 1) {
                return std::shared_ptr<TH1>((TH1 *) h->Clone());
            }
            Array contents = xsec::root::MapContentsToEigen(h.get());
            Array errors = xsec::root::MapErrorsToEigen(h.get());
            auto ret = std::make_shared<TH1D>("", "",
                                              contents.size() - 2, 0, contents.size() - 2);
            for (auto x = 0; x < contents.size(); x++) {
                ret->SetBinContent(x, contents(x));
                ret->SetBinError(x, errors(x));
            }
            return ret;
        }

        std::shared_ptr<TH1>
        UserJointTemplateComponent::
        _join1d(const std::shared_ptr<TH1> lhs, const std::shared_ptr<TH1> rhs) {
            auto lhs1d = _to1d(lhs);
            auto rhs1d = _to1d(rhs);
            return detail::_join(lhs1d, rhs1d);
        }

        std::shared_ptr<TH1>
        UserJointTemplateComponent::
        _join1d(std::vector<std::shared_ptr<TH1>> hists) {
            auto ret = std::shared_ptr<TH1>((TH1 *) hists[0]->Clone());
            for (auto i = 1u; i < hists.size(); i++) {
                ret = _join1d(ret, hists[i]);
            }
            return ret;
        }

        Systematic<TH1>
        UserJointTemplateComponent::
        _to1d(const Systematic<TH1> & syst) {
            std::vector<std::shared_ptr<TH1>> hists1d;
            for (auto i = 0u; i < syst.GetShifts().size(); i++) {
                hists1d.push_back(_to1d(syst.GetShifts()[i]));
            }
            return Systematic<TH1>(syst.GetName(), hists1d, syst.GetType());

        }

        UserJointTemplateComponent::
        UserJointTemplateComponent(const std::map<std::string, const IUserTemplateComponent *> & sample_components)
                : fSamples(sample_components) {
            // joined templates
            std::vector<std::shared_ptr<TH1>> vmeans;
            for (const auto & sample_mean: sample_components) {
                vmeans.push_back(sample_mean.second->GetNominal());
            }
            fJointTemplateMean = detail::_join(vmeans);
            for(const auto & syst : sample_components.begin()->second->GetSystematics()) {
                std::vector<Systematic<TH1>> vsysts;
                for(const auto & sample : sample_components) {
                    vsysts.push_back(sample_components.at(sample.first)->GetSystematics().at(syst.first));
                }
                fJointTemplateSystematics[syst.first] = detail::_join(vsysts);
            }


            // project sample means
            std::vector<std::shared_ptr<TH1>> projected_sample_means;
            for (const auto & sample: sample_components) {
                fSampleNormalizations[sample.first] = _to1d(
                        std::shared_ptr<TH1>(ComponentReducer::Project(sample.second->GetNominalForErrorCalculation()))
                );
                projected_sample_means.push_back(fSampleNormalizations.at(sample.first));
            }

            // join projected means
            fJointNormalization = detail::_join(projected_sample_means);

            // join projected systematics
            // loop over systematics
            for (const auto & syst: sample_components.begin()->second->GetSystematics()) {
                std::vector<Systematic<TH1>> systematics1d;
                // loop over samples
                for (const auto & sample: sample_components) {
                    fSampleNormSystematics[sample.first][syst.first] =
                            ComponentReducer::Project(sample_components.at(sample.first)->GetSystematics().at(syst.first));
                    systematics1d.push_back(_to1d(fSampleNormSystematics.at(sample.first).at(syst.first)));
                }

                fJointNormSystematics[syst.first] = detail::_join(systematics1d);
            }
            fJointNormalizationTotalCovariance = 0;
            for (auto syst: fJointNormSystematics) {
                fJointNormalizationCovariances[syst.first] = syst.second.CovarianceMatrix(fJointNormalization.get());
                if (!fJointNormalizationTotalCovariance) {
                    fJointNormalizationTotalCovariance = (TH1 *) fJointNormalizationCovariances.at(
                            syst.first)->Clone();
                } else {
                    fJointNormalizationTotalCovariance->Add(fJointNormalizationCovariances.at(syst.first));
                }
            }

            for (auto sample: fSampleNormSystematics) {
                fSampleNormalizationTotalCovariance[sample.first] = 0;
                for (auto syst: fSampleNormSystematics.at(sample.first)) {
                    fSampleNormalizationCovariances[sample.first][syst.first] =
                            syst.second.CovarianceMatrix(fSampleNormalizations.at(sample.first).get());
                    if (!fSampleNormalizationTotalCovariance.at(sample.first)) {
                        fSampleNormalizationTotalCovariance.at(sample.first) =
                                (TH1 *) fSampleNormalizationCovariances.at(sample.first).at(
                                        syst.first)->Clone();
                    } else {
                        fSampleNormalizationTotalCovariance.at(sample.first)->Add(
                                fSampleNormalizationCovariances.at(sample.first).at(syst.first)
                        );
                    }
                }
            }

            for (int x = 0; x < fJointNormalization->GetNbinsX() + 2; x++) {
                fJointNormalization->SetBinError(x,
                                                 std::sqrt(fJointNormalizationTotalCovariance->GetBinContent(x,
                                                                                                             x)));
            }
            for (auto sample: fSampleNormalizations) {
                for (int x = 0; x < fSampleNormalizations.at(sample.first)->GetNbinsX() + 2; x++) {
                    fSampleNormalizations.at(sample.first)->SetBinError(
                            x,
                            std::sqrt(fSampleNormalizationTotalCovariance.at(sample.first)->GetBinContent(x, x))
                    );
                }
            }

            fConditioningSampleLabel = fSampleNormalizations.begin()->first;
        }

        IReducedTemplateComponent *
        UserJointTemplateComponent::
        Reduce(const ComponentReducer & reducer) const {
            // loop over systematics
            std::map<std::string, Systematic<TH1>> compressed_joined_norm_systematics1d;
            for (const auto & syst: fSampleNormSystematics.begin()->second) {
                // loop over samples
                std::vector<Systematic<TH1>> systematics1d;
                for (const auto & sample: fSampleNormSystematics) {
                    std::vector<std::shared_ptr<TH1>> compressed_shifts;
                    for (auto i = 0u;
                         i < fSampleNormSystematics.at(sample.first).at(syst.first).GetShifts().size(); i++) {
                        compressed_shifts.push_back(
                                std::shared_ptr<TH1>(reducer.Compress1D(
                                        fSampleNormSystematics.at(sample.first).at(syst.first).GetShifts()[i]))
                        );
                    }
                    systematics1d.emplace_back(fSampleNormSystematics.at(sample.first).at(syst.first).GetName(),
                                               compressed_shifts,
                                               fSampleNormSystematics.at(sample.first).at(syst.first).GetType());
                }
                compressed_joined_norm_systematics1d[syst.first] = detail::_join(systematics1d);
            }

            std::vector<std::shared_ptr<TH1>> compressed_norm_means_to_join;
            std::map<std::string, std::shared_ptr<TH1>> compressed_sample_norm_means;
            for (const auto & sample: fSampleNormalizations) {
                compressed_sample_norm_means[sample.first] = std::shared_ptr<TH1>(reducer.Compress1D(sample.second));
                compressed_norm_means_to_join.push_back(compressed_sample_norm_means.at(sample.first));
            }

            std::map<std::string, const IReducedTemplateComponent *> reduced_samples;
            for (const auto & sample: fSamples) {
                reduced_samples[sample.first] = sample.second->Reduce(reducer);
            }
            return new ReducedJointTemplateComponent(reduced_samples,
                                                     reducer.Reduce(fJointTemplateMean),
                                                     detail::_join(compressed_norm_means_to_join),
                                                     compressed_sample_norm_means,
                                                     compressed_joined_norm_systematics1d);
        }

        ReducedJointTemplateComponent::
        ReducedJointTemplateComponent(std::map<std::string, const IReducedTemplateComponent *> samples,
                                      const ReducedComponent * joined_template_mean,
                                      std::shared_ptr<TH1> joined_norm_mean,
                                      std::map<std::string, std::shared_ptr<TH1>> & sample_norm_means,
                                      std::map<std::string, Systematic<TH1>> & joint_norm_systematics)
                : fJointTemplateMean(joined_template_mean),
                  fJointNormalization(joined_norm_mean),
                  fSampleNormalizations(sample_norm_means),
                  fJointNormSystematics(joint_norm_systematics),
                  fNOuterBins(joined_template_mean->GetNOuterBins()) {
            assert(sample_norm_means.size() == 2 && "Implementation only supports 2 sample joint fits");

            for(const auto & syst : samples.begin()->second->GetSystematics()) {
                std::vector<Systematic<TH1>> vsysts;
                for(const auto & sample : samples) {
                    vsysts.push_back(samples.at(sample.first)->GetSystematics().at(syst.first));
                }
                fJointTemplateSystematics[syst.first] = detail::_join(vsysts);
            }

            fJoinedNormTotalCovariance = nullptr;
            for (const auto & syst: joint_norm_systematics) {
                fJoinedNormCovariances[syst.first] = syst.second.CovarianceMatrix(fJointNormalization.get());
                if (!fJoinedNormTotalCovariance) {
                    fJoinedNormTotalCovariance = (TH1 *) fJoinedNormCovariances.at(syst.first)->Clone();
                } else {
                    fJoinedNormTotalCovariance->Add(fJoinedNormCovariances.at(syst.first));
                }
            }

            // make sure bins didn't get scrambled at some point
            // errors calculated for the reduced nominal normalization should be
            // the same as the newly created covariance matrix from reduced systematics
            for (int x = 1; x < fJointNormalization->GetNbinsX() + 2; x++) {
                assert(fJointNormalization->GetBinError(x) ==
                       std::sqrt(fJoinedNormTotalCovariance->GetBinContent(x, x)));
            }

            auto sample_it = fSampleNormalizations.begin();
            fConditioningSampleNormalization = root::MapContentsToEigenInner(sample_it->second.get());
            fConditioningSampleLabel = sample_it->first;
            fConditioningSample = samples.at(fConditioningSampleLabel);

            sample_it++;
            fComplimentarySampleNormalization = root::MapContentsToEigenInner(sample_it->second.get());
            fComplimentarySampleLabel = sample_it->first;
            fComplimentarySample = samples.at(fComplimentarySampleLabel);

            if(_is_singular()) {
                // If covariance matrix is singular, treat samples as 100% correlated
                // We can do this by fudging the covariance matrices
                Array inverse = 1. / fConditioningSampleNormalization;
                // replace infinities with 0
                inverse = (fConditioningSampleNormalization == 0).select(0, inverse);
                fConditioningSampleInvCovariance = inverse.matrix().asDiagonal();
                fCrossSampleCovariance = fComplimentarySampleNormalization.asDiagonal();
            }
            else {
                Matrix total_covariance = root::MapContentsToEigenInner(fJoinedNormTotalCovariance)
                    .reshaped(fJoinedNormTotalCovariance->GetNbinsX(),
                              fJoinedNormTotalCovariance->GetNbinsX());

                fConditioningSampleInvCovariance =
                        total_covariance.block(0,
                                               0,
                                               fNOuterBins,
                                               fNOuterBins).inverse();
                fCrossSampleCovariance =
                        total_covariance.block(fNOuterBins,
                                               0,
                                               fNOuterBins,
                                               fNOuterBins);
            }
        }

        bool
        ReducedJointTemplateComponent::
        _is_singular() {
            fIsSingular = false;
            for(int i = 0; i < fComplimentarySampleNormalization.size(); i++) {
                if(fComplimentarySampleNormalization(i) == 0) fIsSingular = true;
                if(fConditioningSampleNormalization(i) == 0) fIsSingular = true;
            }
            return fIsSingular;
        }

        Vector
        ReducedJointTemplateComponent::
        ComplimentaryParams(const Vector & conditioning_params) const {
            auto rescaled_condi_sample = fConditioningSampleNormalization * conditioning_params.array();
            Array rescaled_comp_sample =
                    fComplimentarySampleNormalization +
                    fCrossSampleCovariance * fConditioningSampleInvCovariance *
                    (rescaled_condi_sample - fConditioningSampleNormalization).matrix();
            Array complimentary_params = rescaled_comp_sample / fComplimentarySampleNormalization.array();
            return (fComplimentarySampleNormalization.array() == 0).select(0, complimentary_params);
        }

        Vector
        ReducedJointTemplateComponent::
        Predict(const Vector & condi_params) const {
            Vector comp_params = this->ComplimentaryParams(condi_params);
            Matrix condi_prediction = fConditioningSample->Predict(condi_params)
                    .reshaped(fConditioningSample->GetNominal()->GetNInnerBins(), fNOuterBins);
            Matrix comp_prediction = fComplimentarySample->Predict(comp_params)
                    .reshaped(fComplimentarySample->GetNominal()->GetNInnerBins(), fNOuterBins);

            Matrix joint_prediction(condi_prediction.rows() + comp_prediction.rows(),
                                    fNOuterBins);
            joint_prediction << condi_prediction, comp_prediction;
            return joint_prediction.reshaped();
        }

        /*
        IReducedTemplateComponent *
        UserSingleSampleTemplateComponent::
        Reduce(const ComponentReducer & reducer) const {
            std::map<std::string, Systematic<TH1>> reduced_systematics;
            for(const auto & syst : fSystematics) {
                reduced_systematics[syst.first] = reducer.Reduce(syst.second);
            }
            return new ReducedSingleSampleTemplateComponent(reducer.Reduce(fNominal),
                                                            reduced_systematics);
        }
         */
    }
}

