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
#include "XSecAna/Math.h"

namespace xsec {
    ICrossSection::
    ICrossSection(IEfficiencyEstimator * efficiency,
                  ISignalEstimator * signal_estimator,
                  IFluxEstimator * flux,
                  IUnfoldEstimator * unfold,
                  double ntargets)
            : fEfficiency(efficiency),
              fSignalEstimator(signal_estimator),
              fFlux(flux),
              fUnfold(unfold),
              fNTargets(ntargets) {}

    EigenDifferentialCrossSectionEstimator::
    EigenDifferentialCrossSectionEstimator(IEfficiencyEstimator * efficiency,
                                           ISignalEstimator * signal_estimator,
                                           IFluxEstimator * flux,
                                           IUnfoldEstimator * unfold,
                                           double ntargets)
            : ICrossSection(efficiency,
                            signal_estimator,
                            flux,
                            unfold,
                            ntargets) {}

    EigenCrossSectionEstimator::
    EigenCrossSectionEstimator(IEfficiencyEstimator * efficiency,
                               ISignalEstimator * signal_estimator,
                               IFluxEstimator * flux,
                               IUnfoldEstimator * unfold,
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
        return unfolded_selected_signal / denom;
    }

    const TH1 *
    CalculateCrossSection(const TH1 * unfolded_selected_signal,
                          const TH1 * efficiency,
                          const TH1 * flux,
                          const double ntargets,
                          const bool is_differential) {

        root::TH1Props prop(unfolded_selected_signal, "cross_section");
        Array bin_widths;
        if (is_differential) {
            bin_widths = Array::Ones(prop.nbins_and_uof);
        } else {
            bin_widths = root::MapBinWidthsToEigen(unfolded_selected_signal);
        }
        auto result = CalculateCrossSection(root::MapContentsToEigen(unfolded_selected_signal),
                                            root::MapContentsToEigen(efficiency),
                                            root::MapContentsToEigen(flux),
                                            ntargets,
                                            bin_widths);
        return root::ToROOT(result, result, prop);
        //auto xsec = unfolded_selected_signal->Clone();
        //
        //
        //xsec = xsec->TrueDivide(efficiency);
        //xsec = xsec->Divide(flux);
        //xsec = xsec->Divide(1e4 / ntargets); // Convert nu/m^2 to nu/cm^2
        //if (is_differential) xsec = xsec->BinWidthNormalize();
        //return xsec;
    }

    ////////////////////////////////////////////////
    void
    EigenCrossSectionEstimator::
    _eval_impl(const Array & data, const Array & error,
               ArrayRef result, ArrayRef rerror) const {
        Array signal_c(data.size()), signal_e(data.size());
        dynamic_cast<IEigenSignalEstimator*>(fSignalEstimator)->_eval_impl(data, error, signal_c, signal_e);
        dynamic_cast<IEigenUnfoldEstimator*>(fUnfold)->_eval_impl(data, error, signal_c, signal_e);


        Array eff_c(data.size()), eff_e(data.size());
        dynamic_cast<IEigenEfficiencyEstimator*>(fEfficiency)->_eval_impl(data, error, eff_c, eff_e);

        Array flux_c(data.size()), flux_e(data.size());
        dynamic_cast<IEigenFluxEstimator*>(fFlux)->_eval_impl(data, error, flux_c, flux_e);

        result = CalculateCrossSection(signal_c,
                                       eff_c,
                                       flux_c,
                                       fNTargets,
                                       Array::Ones(data.size()));
        rerror = QuadSum(signal_e, eff_e, flux_e);
    }

    ////////////////////////////////////////////////
    void
    EigenDifferentialCrossSectionEstimator::
    _eval_impl(const Array & data, const Array & error,
               ArrayRef result, ArrayRef rerror) const {
        Array signal_c(data.size()), signal_e(data.size());
        dynamic_cast<IEigenSignalEstimator*>(fSignalEstimator)->_eval_impl(data, error, signal_c, signal_e);
        dynamic_cast<IEigenUnfoldEstimator*>(fUnfold)->_eval_impl(data, error, signal_c, signal_e);


        Array eff_c(data.size()), eff_e(data.size());
        dynamic_cast<IEigenEfficiencyEstimator*>(fEfficiency)->_eval_impl(data, error, eff_c, eff_e);

        Array flux_c(data.size()), flux_e(data.size());
        dynamic_cast<IEigenFluxEstimator*>(fFlux)->_eval_impl(data, error, flux_c, flux_e);

        auto bin_widths = root::MapBinWidthsToEigen(this->GetHistProps());

        result = CalculateCrossSection(signal_c,
                                       eff_c,
                                       flux_c,
                                       fNTargets,
                                       bin_widths);

        rerror = QuadSum(signal_e, eff_e, flux_e);
    }

    /*
    const _hist *
    EigenCrossSectionEstimator::
    Eval(const _hist * data) const {
        // calculate estimated signal and scale by the exposure of the data
        auto signal = fSignalEstimator->Signal(data);
        signal = signal->ScaleByExposure(data->Exposure());
        return CalculateCrossSection(fUnfold->Truth(signal),
                                     fEfficiency->Eval(),
                                     fFlux->Eval(data),
                                     fNTargets,
                                     false);
    }
      */
    /*
    ////////////////////////////////////////////////
    const _hist *
    EigenDifferentialCrossSectionEstimator::
    Eval(const _hist * data) const {
        // calculate estimated signal and scale by the exposure of the data
        auto signal = fSignalEstimator->Signal(data);
        signal = signal->ScaleByExposure(data->Exposure());

        return CalculateCrossSection(fUnfold->Truth(signal),
                                     fEfficiency->Eval(),
                                     fFlux->Eval(data),
                                     fNTargets,
                                     true);
    }
     */

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
    LoadFrom(xsec::type::LoadFunction<IEfficiencyEstimator> load_efficiency,
             xsec::type::LoadFunction<ISignalEstimator> load_signal,
             xsec::type::LoadFunction<IFluxEstimator> load_flux,
             xsec::type::LoadFunction<IUnfoldEstimator> load_unfold,
             TDirectory * dir,
             const std::string & subdir) {

        dir = dir->GetDirectory(subdir.c_str());

        IEfficiencyEstimator * eff = 0;
        IFluxEstimator * flux = 0;
        ISignalEstimator * sig = 0;
        IUnfoldEstimator * unfold = 0;

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
    LoadFrom(xsec::type::LoadFunction<IEfficiencyEstimator> load_efficiency,
             xsec::type::LoadFunction<ISignalEstimator> load_signal,
             xsec::type::LoadFunction<IFluxEstimator> load_flux,
             xsec::type::LoadFunction<IUnfoldEstimator> load_unfold,
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
    LoadFrom(xsec::type::LoadFunction<IEfficiencyEstimator> load_efficiency,
             xsec::type::LoadFunction<ISignalEstimator> load_signal,
             xsec::type::LoadFunction<IFluxEstimator> load_flux,
             xsec::type::LoadFunction<IUnfoldEstimator> load_unfold,
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
