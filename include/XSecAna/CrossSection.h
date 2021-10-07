#pragma once

#include "TDirectory.h"
#include "TParameter.h"

#include "Hist.h"

#include "XSecAna/IMeasurement.h"
#include "XSecAna/IEfficiency.h"
#include "XSecAna/ISignalEstimator.h"
#include "XSecAna/IFlux.h"
#include "XSecAna/IUnfold.h"

namespace xsec {
    template<class HistType,
            bool IsDifferential=false>
    class CrossSection : public IMeasurement<HistType> {
    public:
        CrossSection() = default;

        CrossSection(IEfficiency * efficiency,
                     ISignalEstimator<HistType> * signal_estimator,
                     IFlux * flux,
                     IUnfold * unfold,
                     double ntargets = 0)
                : fEfficiency(efficiency),
                  fSignalEstimator(signal_estimator),
                  fFlux(flux),
                  fUnfold(unfold),
                  fNTargets(ntargets) {}

        void SetNTargets(double ntargets) { fNTargets = ntargets; }

        void SaveTo(TDirectory * dir, const std::string & subdir) const override;

        static std::unique_ptr<IMeasurement<HistType>> LoadFrom(xsec::type::LoadFunction<IEfficiency> load_efficiency,
                                                                xsec::type::LoadFunction<ISignalEstimator<HistType>> load_signal,
                                                                xsec::type::LoadFunction<IFlux> load_flux,
                                                                xsec::type::LoadFunction<IUnfold> load_unfold,
                                                                TDirectory * dir,
                                                                const std::string & subdir);

        HistType Eval(const HistType & data) const override;

        CrossSection<HistType,
                     true>
        ToDifferential();

    private:
        IEfficiency * fEfficiency;
        ISignalEstimator<HistType> * fSignalEstimator;
        IFlux * fFlux;
        IUnfold * fUnfold;

        double fNTargets;
    };

    ////////////////////////////////////////////////
    template<class HistType,
            bool IsDifferential>
    CrossSection<HistType,
                 true>
    CrossSection<HistType,
                 IsDifferential>::
    ToDifferential() {
        if constexpr(IsDifferential) {
            return *this;
        } else {
            return CrossSection<HistType,
                                true>(fEfficiency,
                                      fSignalEstimator,
                                      fFlux,
                                      fUnfold,
                                      fNTargets);
        }
    }

    ////////////////////////////////////////////////
    template<class HistType>
    HistType CalculateCrossSection(const HistType & unfolded_selected_signal,
                                   const HistType & efficiency,
                                   const HistType & flux,
                                   const double ntargets,
                                   const bool is_differential) {
        HistType xsec = unfolded_selected_signal;

        // don't scale efficiency by exposure
        // since exposure cancels in the ratio
        xsec = xsec.TrueDivide(efficiency);
        xsec /= flux;
        xsec *= 1e4 / ntargets; // Convert nu/m^2 to nu/cm^2
        if (is_differential) xsec = xsec.BinWidthNormalize();
        return xsec;
    }

    ////////////////////////////////////////////////
    template<class HistType,
            bool IsDifferential>
    HistType
    CrossSection<HistType,
                 IsDifferential>::
    Eval(const HistType & data) const {
        // calculate estimated signal and scale by the exposure of the data
        auto signal = fSignalEstimator->Signal(data);
        signal = signal.ScaleByExposure(data.Exposure());

        // invoke FluxType::operator* in case we want the integrated flux
        // Pass a histogram of ones through the flux parameter of
        // CalculateCrossSection to divide by one

        auto flux = fFlux->Eval(data.EdgesAndUOF());

        return CalculateCrossSection(fUnfold->Truth(signal),
                                     fEfficiency->Eval(),
                                     fFlux->Eval(data.EdgesAndUOF()),
                                     fNTargets,
                                     IsDifferential);
    }

    ////////////////////////////////////////////////
    template<class HistType,
            bool IsDifferential>
    void
    CrossSection<HistType,
                 IsDifferential>::
                 SaveTo(TDirectory * dir, const std::string & subdir) const {
        TDirectory * tmp = gDirectory;
        dir = dir->mkdir(subdir.c_str());
        dir->cd();
        TObjString("CrossSection").Write("type");

        if (fEfficiency) fEfficiency->SaveTo(dir, "fEfficiency");
        if (fFlux) fFlux->SaveTo(dir, "fFlux");
        if (fSignalEstimator) fSignalEstimator->SaveTo(dir, "fSignalEstimator");
        if (fUnfold) fUnfold->SaveTo(dir, "fUnfold");

        TParameter<double>("fNTargets", fNTargets).Write("fNTargets");

        tmp->cd();
    }

    ////////////////////////////////////////////////
    template<class HistType,
            bool IsDifferential>
    std::unique_ptr<IMeasurement<HistType>>
    CrossSection<HistType,
                 IsDifferential>::
    LoadFrom(xsec::type::LoadFunction<IEfficiency> load_efficiency,
             xsec::type::LoadFunction<ISignalEstimator<HistType>> load_signal,
             xsec::type::LoadFunction<IFlux> load_flux,
             xsec::type::LoadFunction<IUnfold> load_unfold,
             TDirectory * dir,
             const std::string & subdir) {

        dir = dir->GetDirectory(subdir.c_str());

        IEfficiency * eff = 0;
        IFlux * flux = 0;
        ISignalEstimator<HistType> * sig = 0;
        IUnfold * unfold = 0;

        if (dir->GetDirectory("fEfficiency")) eff = load_efficiency(dir, "fEfficiency").release();
        if (dir->GetDirectory("fFlux")) flux = load_flux(dir, "fFlux").release();
        if (dir->GetDirectory("fSignalEstimator"))
            sig = load_signal(dir, "fSignalEstimator").release();
        if (dir->GetDirectory("fUnfold")) unfold = load_unfold(dir, "fUnfold").release();

        auto ntargets = ((TParameter<typename HistType::scalar_type> *) dir->Get("fNTargets"))->GetVal();

        return std::make_unique<CrossSection<HistType,
                                             IsDifferential> >
                (eff, sig, flux, unfold, ntargets);
    }
}
