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

        CrossSection(IEfficiency<HistType> * efficiency,
                     ISignalEstimator<HistType> * signal_estimator,
                     IFlux<HistType> * flux,
                     IUnfold<HistType> * unfold,
                     typename HistType::scalar_type ntargets = 0)
                : fEfficiency(efficiency),
                  fSignalEstimator(signal_estimator),
                  fFlux(flux),
                  fUnfold(unfold),
                  fNTargets(ntargets) {}

        void SetNTargets(double ntargets) { fNTargets = ntargets; }

        void SaveTo(TDirectory * dir, const std::string & subdir) const override;

        /// TODO how do we serialize user's objects that inherit from the interfaces?
        static std::unique_ptr<CrossSection> LoadFrom(TDirectory * dir, const std::string & subdir);

        HistType Eval(const HistType & data) override;

        CrossSection<HistType,
                     true>
        ToDifferential();

    private:
        IEfficiency<HistType> * fEfficiency;
        ISignalEstimator<HistType> * fSignalEstimator;
        IFlux<HistType> * fFlux;
        IUnfold<HistType> * fUnfold;

        typename HistType::scalar_type fNTargets;
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
                                   const typename HistType::scalar_type ntargets,
                                   const bool is_differential) {
        HistType xsec = unfolded_selected_signal;

        // don't scale efficiency by exposure
        // since exposure cancels in the ratio
        xsec = xsec.TrueDivide(efficiency);
        xsec /= flux;
        xsec *= 1e4 / ntargets; // Convert nu/m^2 to nu/cm^2
        if (is_differential) xsec.Normalize("width");
        return xsec;
    }

    ////////////////////////////////////////////////
    template<class HistType,
            bool IsDifferential>
    HistType
    CrossSection<HistType,
                 IsDifferential>::
    Eval(const HistType & data) {
        // calculate estimated signal and scale by the exposure of the data
        auto signal = fSignalEstimator->Signal(data);
        signal = signal.ScaleByExposure(data.Exposure());

        // invoke FluxType::operator* in case we want the integrated flux
        // Pass a histogram of ones through the flux parameter of
        // CalculateCrossSection to divide by one
        return CalculateCrossSection(fUnfold->Truth(signal),
                                     (*fFlux * fEfficiency->Eval()),
                                     HistType(HistType::array_type::Ones(fEfficiency->Eval().Contents().size()),
                                              fEfficiency->Eval().Edges()),
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

        TParameter<typename HistType::scalar_type>("fNTargets", fNTargets).Write("fNTargets");

        tmp->cd();
    }

    ////////////////////////////////////////////////
    template<class HistType,
            bool IsDifferential>
    std::unique_ptr<CrossSection<HistType,
                                 IsDifferential> >
    CrossSection<HistType,
                 IsDifferential>::
    LoadFrom(TDirectory * dir, const std::string & subdir) {
        dir = dir->GetDirectory(subdir.c_str());

        IEfficiency<HistType> * eff = 0;
        IFlux<HistType> * flux = 0;
        ISignalEstimator<HistType> * sig = 0;
        IUnfold<HistType> * unfold = 0;

        if (dir->GetDirectory("fEfficiency")) eff = IEfficiency<HistType>::LoadFrom(dir, "fEfficiency").release();
        if (dir->GetDirectory("fFlux")) flux = IFlux<HistType>::LoadFrom(dir, "fFlux").release();
        if (dir->GetDirectory("fSignalEstimator"))
            sig = ISignalEstimator<HistType>::LoadFrom(dir, "fSignalEstimator").release();
        if (dir->GetDirectory("fUnfold")) unfold = IUnfold<HistType>::LoadFrom(dir, "fUnfold").release();

        auto ntargets = ((TParameter<typename HistType::scalar_type> *) dir->Get("fNTargets"))->GetVal();

        return std::make_unique<CrossSection<HistType,
                                             IsDifferential> >
                (eff, sig, flux, unfold, ntargets);
    }
}
