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
    class CrossSectionBase : public IMeasurement<Hist> {
    public:
        CrossSectionBase() = default;

        CrossSectionBase(IEfficiency * efficiency,
                         ISignalEstimator * signal_estimator,
                         IFlux * flux,
                         IUnfold * unfold,
                         double ntargets = 0)
                : fEfficiency(efficiency),
                  fSignalEstimator(signal_estimator),
                  fFlux(flux),
                  fUnfold(unfold),
                  fNTargets(ntargets) {}

        void SaveTo(TDirectory * dir, const std::string & subdir) const override;

        template<class DerivedCrossSection>
        static std::unique_ptr<IMeasurement<Hist>> LoadFrom(xsec::type::LoadFunction<IEfficiency> load_efficiency,
                                                            xsec::type::LoadFunction<ISignalEstimator> load_signal,
                                                            xsec::type::LoadFunction<IFlux> load_flux,
                                                            xsec::type::LoadFunction<IUnfold> load_unfold,
                                                            TDirectory * dir,
                                                            const std::string & subdir);

        virtual Hist Eval(const Hist & data) const override = 0;


    protected:
        IEfficiency * fEfficiency;
        ISignalEstimator * fSignalEstimator;
        IFlux * fFlux;
        IUnfold * fUnfold;

        double fNTargets;
    };

    class DifferentialCrossSection : public CrossSectionBase {
    public:
        DifferentialCrossSection() = default;

        DifferentialCrossSection(IEfficiency * efficiency,
                                 ISignalEstimator * signal_estimator,
                                 IFlux * flux,
                                 IUnfold * unfold,
                                 double ntargets = 0)
                : CrossSectionBase(efficiency,
                                   signal_estimator,
                                   flux,
                                   unfold,
                                   ntargets) {}

        virtual Hist Eval(const Hist & data) const override;

        static std::unique_ptr<IMeasurement<Hist>> LoadFrom(xsec::type::LoadFunction<IEfficiency> load_efficiency,
                                                            xsec::type::LoadFunction<ISignalEstimator> load_signal,
                                                            xsec::type::LoadFunction<IFlux> load_flux,
                                                            xsec::type::LoadFunction<IUnfold> load_unfold,
                                                            TDirectory * dir,
                                                            const std::string & subdir);
    };

    class CrossSection : public CrossSectionBase {
    public:
        CrossSection() = default;

        CrossSection(IEfficiency * efficiency,
                     ISignalEstimator * signal_estimator,
                     IFlux * flux,
                     IUnfold * unfold,
                     double ntargets = 0)
                : CrossSectionBase(efficiency,
                                   signal_estimator,
                                   flux,
                                   unfold,
                                   ntargets) {}

        virtual Hist Eval(const Hist & data) const override;

        static std::unique_ptr<IMeasurement<Hist>> LoadFrom(xsec::type::LoadFunction<IEfficiency> load_efficiency,
                                                            xsec::type::LoadFunction<ISignalEstimator> load_signal,
                                                            xsec::type::LoadFunction<IFlux> load_flux,
                                                            xsec::type::LoadFunction<IUnfold> load_unfold,
                                                            TDirectory * dir,
                                                            const std::string & subdir);
    };

    ////////////////////////////////////////////////
    Hist CalculateCrossSection(const Hist & unfolded_selected_signal,
                               const Hist & efficiency,
                               const Hist & flux,
                               const double ntargets,
                               const bool is_differential) {
        Hist xsec = unfolded_selected_signal;

        // don't scale efficiency by exposure
        // since exposure cancels in the ratio
        xsec = xsec.TrueDivide(efficiency);
        xsec /= flux;
        xsec *= 1e4 / ntargets; // Convert nu/m^2 to nu/cm^2
        if (is_differential) xsec = xsec.BinWidthNormalize();
        return xsec;
    }

    ////////////////////////////////////////////////
    Hist
    CrossSection::
    Eval(const Hist & data) const {
        // calculate estimated signal and scale by the exposure of the data
        auto signal = fSignalEstimator->Signal(data);
        signal = signal.ScaleByExposure(data.Exposure());

        auto flux = fFlux->Eval(data.EdgesAndUOF());

        return CalculateCrossSection(fUnfold->Truth(signal),
                                     fEfficiency->Eval(),
                                     fFlux->Eval(data.EdgesAndUOF()),
                                     fNTargets,
                                     false);
    }

    ////////////////////////////////////////////////
    Hist
    DifferentialCrossSection::
    Eval(const Hist & data) const {
        // calculate estimated signal and scale by the exposure of the data
        auto signal = fSignalEstimator->Signal(data);
        signal = signal.ScaleByExposure(data.Exposure());

        auto flux = fFlux->Eval(data.EdgesAndUOF());

        return CalculateCrossSection(fUnfold->Truth(signal),
                                     fEfficiency->Eval(),
                                     fFlux->Eval(data.EdgesAndUOF()),
                                     fNTargets,
                                     true);
    }

    ////////////////////////////////////////////////
    void
    CrossSectionBase::
    SaveTo(TDirectory * dir, const std::string & subdir) const {
        TDirectory * tmp = gDirectory;
        dir = dir->mkdir(subdir.c_str());
        dir->cd();
        TObjString("CrossSectionBase").Write("type");

        if (fEfficiency) fEfficiency->SaveTo(dir, "fEfficiency");
        if (fFlux) fFlux->SaveTo(dir, "fFlux");
        if (fSignalEstimator) fSignalEstimator->SaveTo(dir, "fSignalEstimator");
        if (fUnfold) fUnfold->SaveTo(dir, "fUnfold");

        TParameter<double>("fNTargets", fNTargets).Write("fNTargets");

        tmp->cd();
    }

    ////////////////////////////////////////////////
    template<class DerivedCrossSection>
    std::unique_ptr<IMeasurement<Hist>>
    CrossSectionBase::
    LoadFrom(xsec::type::LoadFunction<IEfficiency> load_efficiency,
             xsec::type::LoadFunction<ISignalEstimator> load_signal,
             xsec::type::LoadFunction<IFlux> load_flux,
             xsec::type::LoadFunction<IUnfold> load_unfold,
             TDirectory * dir,
             const std::string & subdir) {

        dir = dir->GetDirectory(subdir.c_str());

        IEfficiency * eff = 0;
        IFlux * flux = 0;
        ISignalEstimator * sig = 0;
        IUnfold * unfold = 0;

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
    std::unique_ptr<IMeasurement<Hist>>
    DifferentialCrossSection::
    LoadFrom(xsec::type::LoadFunction<IEfficiency> load_efficiency,
             xsec::type::LoadFunction<ISignalEstimator> load_signal,
             xsec::type::LoadFunction<IFlux> load_flux,
             xsec::type::LoadFunction<IUnfold> load_unfold,
             TDirectory * dir,
             const std::string & subdir) {
        return CrossSectionBase::LoadFrom<DifferentialCrossSection>(load_efficiency,
                                                                    load_signal,
                                                                    load_flux,
                                                                    load_unfold,
                                                                    dir,
                                                                    subdir);
    }

    ////////////////////////////////////////////////
    std::unique_ptr<IMeasurement<Hist>>
    CrossSection::
    LoadFrom(xsec::type::LoadFunction<IEfficiency> load_efficiency,
             xsec::type::LoadFunction<ISignalEstimator> load_signal,
             xsec::type::LoadFunction<IFlux> load_flux,
             xsec::type::LoadFunction<IUnfold> load_unfold,
             TDirectory * dir,
             const std::string & subdir) {
        return CrossSectionBase::LoadFrom<CrossSection>(load_efficiency,
                                                        load_signal,
                                                        load_flux,
                                                        load_unfold,
                                                        dir,
                                                        subdir);
    }
}
