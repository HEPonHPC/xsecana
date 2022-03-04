//
// Created by Derek Doyle on 10/9/21.
//
#include "XSecAna/CrossSection.h"

#include "TDirectory.h"
#include "TParameter.h"

#include "XSecAna/IMeasurement.h"
#include "XSecAna/ISignalEstimator.h"
#include "XSecAna/IEfficiency.h"
#include "XSecAna/IFlux.h"
#include "XSecAna/IUnfold.h"

namespace xsec {
    ICrossSection::
    ICrossSection(IMeasurement * efficiency,
                  IMeasurement * signal_estimator,
                  IMeasurement * flux,
                  IMeasurement * unfold,
                  double ntargets)
            : fEfficiency(efficiency),
              fSignalEstimator(signal_estimator),
              fFlux(flux),
              fUnfold(unfold),
              fNTargets(ntargets) {}

    DifferentialCrossSectionEstimator::
    DifferentialCrossSectionEstimator(IMeasurement * efficiency,
                                      IMeasurement * signal_estimator,
                                      IMeasurement * flux,
                                      IMeasurement * unfold,
                                      double ntargets)
            : ICrossSection(efficiency,
                            signal_estimator,
                            flux,
                            unfold,
                            ntargets) {}


    CrossSectionEstimator::
    CrossSectionEstimator(IMeasurement * efficiency,
                          IMeasurement * signal_estimator,
                          IMeasurement * flux,
                          IMeasurement * unfold,
                          double ntargets)
            : ICrossSection(efficiency,
                            signal_estimator,
                            flux,
                            unfold,
                            ntargets) {}

    EigenDifferentialCrossSectionEstimator::
    EigenDifferentialCrossSectionEstimator(IMeasurement * efficiency,
                                           IMeasurement * signal_estimator,
                                           IMeasurement * flux,
                                           IMeasurement * unfold,
                                           double ntargets)
            : ICrossSection(efficiency,
                            signal_estimator,
                            flux,
                            unfold,
                            ntargets) {}

    EigenCrossSectionEstimator::
    EigenCrossSectionEstimator(IMeasurement * efficiency,
                               IMeasurement * signal_estimator,
                               IMeasurement * flux,
                               IMeasurement * unfold,
                               double ntargets)
            : ICrossSection(efficiency,
                            signal_estimator,
                            flux,
                            unfold,
                            ntargets) {}


    const Array
    CalculateCrossSection(const Array & unfolded_selected_signal,
                          const Array & efficiency,
                          const Array & flux,
                          const double ntargets,
                          const Array & bin_widths) {
        auto denom = efficiency * flux * ntargets / 1e4 * bin_widths;

        Array result = unfolded_selected_signal / denom;
        result = (denom == 0).select(0, result);
        return result;
    }

    const std::shared_ptr<TH1>
    CalculateCrossSection(const TH1 * unfolded_selected_signal,
                          const TH1 * efficiency,
                          const TH1 * flux,
                          const double ntargets,
                          const bool is_differential) {

        root::TH1Props prop(unfolded_selected_signal);
        Array bin_widths;
        if (is_differential) {
            bin_widths = Array::Ones(prop.nbins_and_uof);
        } else {
            bin_widths = root::MapBinWidthsToEigen(unfolded_selected_signal);
        }
        Array result = CalculateCrossSection(root::MapContentsToEigen(unfolded_selected_signal),
                                             root::MapContentsToEigen(efficiency),
                                             root::MapContentsToEigen(flux),
                                             ntargets,
                                             bin_widths);
        auto hresult = std::shared_ptr<TH1>(root::ToROOT(result, prop));
        for (auto x = 1; x <= hresult->GetNbinsX(); x++) {
            for (auto y = 1; y <= hresult->GetNbinsY(); y++) {
                for (auto z = 1; z <= hresult->GetNbinsZ(); z++) {
                    double r = hresult->GetBinContent(x, y, z);
                    if (!r) continue;

                    double deff = efficiency->GetBinError(x, y, z) / efficiency->GetBinContent(x, y, z);
                    double dflux = flux->GetBinError(x, y, z) / flux->GetBinContent(x, y, z);
                    double dsig = unfolded_selected_signal->GetBinError(x, y, z) /
                                  unfolded_selected_signal->GetBinContent(x, y, z);

                    double e = std::sqrt(deff * deff +
                                         dflux * dflux +
                                         dsig * dsig);
                    hresult->SetBinError(x, y, z, e * r);
                }
            }
        }
        return hresult;
    }


    ////////////////////////////////////////////////
    std::shared_ptr<TH1>
    DifferentialCrossSectionEstimator::
    Eval(const TH1 * data) const {

        auto signal_estimation = fSignalEstimator->Eval(data);
        auto unfolded_signal = fUnfold->Eval(signal_estimation.get());

        auto efficiency = fEfficiency->Eval(data);

        auto flux = fFlux->Eval(data);

        auto result = std::shared_ptr<TH1>(CalculateCrossSection(signal_estimation.get(),
                                                                 efficiency.get(),
                                                                 flux.get(),
                                                                 fNTargets,
                                                                 true));
        return result;
    }

    ////////////////////////////////////////////////
    std::shared_ptr<TH1>
    CrossSectionEstimator::
    Eval(const TH1 * data) const {

        auto signal_estimation = dynamic_cast<ISignalEstimator *>(fSignalEstimator)->Signal(data);
        auto unfolded_signal = dynamic_cast<IUnfoldEstimator *>(fUnfold)->Truth(signal_estimation);

        auto efficiency = fEfficiency->Eval(data);

        auto flux = fFlux->Eval(data);

        auto result = std::shared_ptr<TH1>(CalculateCrossSection(signal_estimation,
                                                                 efficiency.get(),
                                                                 flux.get(),
                                                                 fNTargets,
                                                                 false));
        return result;
    }


    ////////////////////////////////////////////////
    void
    EigenCrossSectionEstimator::
    _eval_impl(const Array & data, const Array & error,
               ArrayRef result, ArrayRef rerror) const {
        Array signal_c(data.size()), signal_e(data.size());
        dynamic_cast<IEigenSignalEstimator *>(fSignalEstimator)->_eval_impl(data, error, signal_c, signal_e);
        dynamic_cast<IEigenUnfoldEstimator *>(fUnfold)->_eval_impl(signal_c, signal_e, signal_c, signal_e);


        Array eff_c(data.size()), eff_e(data.size());
        dynamic_cast<IEigenEfficiencyEstimator *>(fEfficiency)->_eval_impl(data, error, eff_c, eff_e);

        Array flux_c(data.size()), flux_e(data.size());
        dynamic_cast<IEigenFluxEstimator *>(fFlux)->_eval_impl(data, error, flux_c, flux_e);

        result = CalculateCrossSection(signal_c,
                                       eff_c,
                                       flux_c,
                                       fNTargets,
                                       Array::Ones(data.size()));
        rerror = ((signal_e / signal_c).pow(2) +
                  (flux_e / flux_c).pow(2) +
                  (eff_e / eff_c).pow(2)).sqrt() * result;
    }

    ////////////////////////////////////////////////
    void
    EigenDifferentialCrossSectionEstimator::
    _eval_impl(const Array & data, const Array & error,
               ArrayRef result, ArrayRef rerror) const {
        Array signal_c(data.size()), signal_e(data.size());
        dynamic_cast<IEigenSignalEstimator *>(fSignalEstimator)->_eval_impl(data, error, signal_c, signal_e);
        dynamic_cast<IEigenUnfoldEstimator *>(fUnfold)->_eval_impl(signal_c, signal_e, signal_c, signal_e);


        Array eff_c(data.size()), eff_e(data.size());
        dynamic_cast<IEigenEfficiencyEstimator *>(fEfficiency)->_eval_impl(data, error, eff_c, eff_e);

        Array flux_c(data.size()), flux_e(data.size());
        dynamic_cast<IEigenFluxEstimator *>(fFlux)->_eval_impl(data, error, flux_c, flux_e);

        auto bin_widths = root::MapBinWidthsToEigen(this->GetHistProps());

        result = CalculateCrossSection(signal_c,
                                       eff_c,
                                       flux_c,
                                       fNTargets,
                                       bin_widths);

        rerror = ((signal_e / signal_c).pow(2) +
                  (flux_e / flux_c).pow(2) +
                  (eff_e / eff_c).pow(2)).sqrt() * result;
    }

    ////////////////////////////////////////////////
    void
    ICrossSection::
    SaveTo(TDirectory * dir, const std::string & subdir) const {
        TDirectory * tmp = gDirectory;
        dir = dir->mkdir(subdir.c_str());
        dir->cd();
        TObjString("ICrossSection").Write("type");

        if (fEfficiency) fEfficiency->SaveTo(dir, "fEfficiency");
        if (fFlux) fFlux->SaveTo(dir, "fFlux");
        if (fSignalEstimator) fSignalEstimator->SaveTo(dir, "fSignalEstimator");
        if (fUnfold) fUnfold->SaveTo(dir, "fUnfold");

        TParameter<double>("fNTargets", fNTargets).Write("fNTargets");

        tmp->cd();
    }

    ////////////////////////////////////////////////
    template<class DerivedCrossSection>
    std::unique_ptr<IMeasurement>
    ICrossSection::
    LoadFrom(xsec::type::LoadFunction<IMeasurement> load_efficiency,
             xsec::type::LoadFunction<IMeasurement> load_signal,
             xsec::type::LoadFunction<IMeasurement> load_flux,
             xsec::type::LoadFunction<IMeasurement> load_unfold,
             TDirectory * dir,
             const std::string & subdir) {

        dir = dir->GetDirectory(subdir.c_str());

        IMeasurement * eff = 0;
        IMeasurement * flux = 0;
        IMeasurement * sig = 0;
        IMeasurement * unfold = 0;

        if (dir->GetDirectory("fEfficiency")) eff = load_efficiency(dir, "fEfficiency").release();
        if (dir->GetDirectory("fFlux")) flux = load_flux(dir, "fFlux").release();
        if (dir->GetDirectory("fSignalEstimator"))
            sig = load_signal(dir, "fSignalEstimator").release();
        if (dir->GetDirectory("fUnfold")) unfold = load_unfold(dir, "fUnfold").release();

        auto ntargets = ((TParameter<double> *) dir->Get("fNTargets"))->GetVal();

        return std::make_unique<DerivedCrossSection>
                (eff, sig, flux, unfold, ntargets);
    }

    ////////////////////////////////////////////////
    std::unique_ptr<IMeasurement>
    EigenDifferentialCrossSectionEstimator::
    LoadFrom(xsec::type::LoadFunction<IMeasurement> load_efficiency,
             xsec::type::LoadFunction<IMeasurement> load_signal,
             xsec::type::LoadFunction<IMeasurement> load_flux,
             xsec::type::LoadFunction<IMeasurement> load_unfold,
             TDirectory * dir,
             const std::string & subdir) {
        return ICrossSection::LoadFrom<EigenDifferentialCrossSectionEstimator>(load_efficiency,
                                                                               load_signal,
                                                                               load_flux,
                                                                               load_unfold,
                                                                               dir,
                                                                               subdir);
    }


    ////////////////////////////////////////////////
    std::unique_ptr<IMeasurement>
    EigenCrossSectionEstimator::
    LoadFrom(xsec::type::LoadFunction<IMeasurement> load_efficiency,
             xsec::type::LoadFunction<IMeasurement> load_signal,
             xsec::type::LoadFunction<IMeasurement> load_flux,
             xsec::type::LoadFunction<IMeasurement> load_unfold,
             TDirectory * dir,
             const std::string & subdir) {
        return ICrossSection::LoadFrom<EigenCrossSectionEstimator>(load_efficiency,
                                                                   load_signal,
                                                                   load_flux,
                                                                   load_unfold,
                                                                   dir,
                                                                   subdir);
    }
}
