#pragma once

#include "TDirectory.h"
#include "TParameter.h"

#include "Hist.h"

#include "XSecAna/IMeasurement.h"

namespace xsec {
    template< class HistType,
            class SignalEstimatorType,
            class UnfoldType,
            class EfficiencyType,
            class FluxType,
            bool IsDifferential = false >
    class CrossSection : public IMeasurement< HistType > {
    public:
        CrossSection() = default;

        CrossSection(EfficiencyType * efficiency,
                     SignalEstimatorType * signal_estimator,
                     FluxType * flux,
                     UnfoldType * unfold,
                     typename HistType::scalar_type ntargets = 0)
                : fEfficiency(efficiency),
                  fSignalEstimator(signal_estimator),
                  fFlux(flux),
                  fUnfold(unfold),
                  fNTargets(ntargets) {}

        void SetNTargets(double ntargets) { fNTargets = ntargets; }

        void SaveTo(TDirectory * dir, const std::string & subdir) const override;

        static std::unique_ptr< CrossSection > LoadFrom(TDirectory * dir, const std::string& subdir);

        HistType Result(const HistType & data) override;

        CrossSection< HistType,
                SignalEstimatorType,
                UnfoldType,
                EfficiencyType,
                FluxType,
                true >
        ToDifferential();

    private:
        EfficiencyType * fEfficiency;
        SignalEstimatorType * fSignalEstimator;
        FluxType * fFlux;
        UnfoldType * fUnfold;

        typename HistType::scalar_type fNTargets;
    };

    ////////////////////////////////////////////////
    template< class HistType,
            class SignalEstimatorType,
            class UnfoldType,
            class EfficiencyType,
            class FluxType,
            bool IsDifferential >
    CrossSection< HistType,
            SignalEstimatorType,
            UnfoldType,
            EfficiencyType,
            FluxType,
            true >
    CrossSection< HistType,
            SignalEstimatorType,
            UnfoldType,
            EfficiencyType,
            FluxType,
            IsDifferential >::
    ToDifferential() {
        if constexpr(IsDifferential) {
            return *this;
        } else {
            return CrossSection< HistType,
                    SignalEstimatorType,
                    UnfoldType,
                    EfficiencyType,
                    FluxType,
                    true >(fEfficiency,
                           fSignalEstimator,
                           fFlux,
                           fUnfold,
                           fNTargets);
        }
    }

    ////////////////////////////////////////////////
    template< class HistType >
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
    template< class HistType,
            class SignalEstimatorType,
            class UnfoldType,
            class EfficiencyType,
            class FluxType,
            bool IsDifferential >
    HistType
    CrossSection< HistType,
            SignalEstimatorType,
            UnfoldType,
            EfficiencyType,
            FluxType,
            IsDifferential >::
    Result(const HistType & data) {
        // calculate estimated signal and scale by the exposure of the data
        auto signal = fSignalEstimator->Signal(data);
        signal = signal.ScaleByExposure(data.Exposure());

        // invoke FluxType::operator* in case we want the integrated flux
        // Pass a histogram of ones through the flux parameter of
        // CalculateCrossSection to divide by one
        return CalculateCrossSection(fUnfold->Truth(signal),
                                     (*fFlux * fEfficiency->ToHist()),
                                     HistType(HistType::array_type::Ones(fEfficiency->ToHist().Contents().size()),
                                              fEfficiency->ToHist().Edges()),
                                     fNTargets,
                                     IsDifferential);
    }

    ////////////////////////////////////////////////
    template< class HistType,
            class SignalEstimatorType,
            class UnfoldType,
            class EfficiencyType,
            class FluxType,
            bool IsDifferential >
    void
    CrossSection< HistType,
            SignalEstimatorType,
            UnfoldType,
            EfficiencyType,
            FluxType,
            IsDifferential >::
    SaveTo(TDirectory * dir, const std::string & subdir) const {
        TDirectory * tmp = gDirectory;
        dir = dir->mkdir(subdir.c_str());
        dir->cd();
        TObjString("CrossSection").Write("type");

        if (fEfficiency) fEfficiency->SaveTo(dir, "fEfficiency");
        if (fFlux) fFlux->SaveTo(dir, "fFlux");
        if (fSignalEstimator) fSignalEstimator->SaveTo(dir, "fSignalEstimator");
        if (fUnfold) fUnfold->SaveTo(dir, "fUnfold");

        TParameter< typename HistType::scalar_type >("fNTargets", fNTargets).Write("fNTargets");

        tmp->cd();
    }

    ////////////////////////////////////////////////
    template< class HistType,
            class SignalEstimatorType,
            class UnfoldType,
            class EfficiencyType,
            class FluxType,
            bool IsDifferential >
    std::unique_ptr< CrossSection< HistType,
            SignalEstimatorType,
            UnfoldType,
            EfficiencyType,
            FluxType,
            IsDifferential > >
    CrossSection< HistType,
            SignalEstimatorType,
            UnfoldType,
            EfficiencyType,
            FluxType,
            IsDifferential >::
    LoadFrom(TDirectory * dir, const std::string& subdir) {
        dir = dir->GetDirectory(subdir.c_str());

        EfficiencyType * eff = 0;
        FluxType * flux = 0;
        SignalEstimatorType * sig = 0;
        UnfoldType * unfold = 0;

        if (dir->GetDirectory("fEfficiency")) eff = EfficiencyType::LoadFrom(dir, "fEfficiency").release();
        if (dir->GetDirectory("fFlux")) flux = FluxType::LoadFrom(dir, "fFlux").release();
        if (dir->GetDirectory("fSignalEstimator"))
            sig = SignalEstimatorType::LoadFrom(dir, "fSignalEstimator").release();
        if (dir->GetDirectory("fUnfold")) unfold = UnfoldType::LoadFrom(dir, "fUnfold").release();

        auto ntargets = ((TParameter< typename HistType::scalar_type > *) dir->Get("fNTargets"))->GetVal();

        return std::make_unique< CrossSection< HistType,
                SignalEstimatorType,
                UnfoldType,
                EfficiencyType,
                FluxType,
                IsDifferential > >
                (eff, sig, flux, unfold, ntargets);
    }
}
