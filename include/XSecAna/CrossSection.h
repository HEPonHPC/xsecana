#pragma once

#include "TDirectory.h"
#include "TParameter.h"

#include "XSecAna/_Hist.h"

#include "XSecAna/IMeasurement.h"
#include "XSecAna/IEfficiency.h"
#include "XSecAna/ISignalEstimator.h"
#include "XSecAna/IFlux.h"
#include "XSecAna/IUnfold.h"
#include "XSecAna/SimpleQuadSum.h"

namespace xsec {
    class ICrossSection : public virtual IMeasurement {
    public:
        ICrossSection() = default;

        ICrossSection(IEfficiencyEstimator * efficiency,
                      ISignalEstimator * signal_estimator,
                      IFluxEstimator * flux,
                      IUnfoldEstimator * unfold,
                      double ntargets = 0);

        void SaveTo(TDirectory * dir, const std::string & subdir) const override;

        template<class DerivedCrossSection>
        static std::unique_ptr<IMeasurement> LoadFrom(xsec::type::LoadFunction<IEfficiencyEstimator> load_efficiency,
                                                      xsec::type::LoadFunction<ISignalEstimator> load_signal,
                                                      xsec::type::LoadFunction<IFluxEstimator> load_flux,
                                                      xsec::type::LoadFunction<IUnfoldEstimator> load_unfold,
                                                      TDirectory * dir,
                                                      const std::string & subdir);

    protected:
        IEfficiencyEstimator * fEfficiency;
        ISignalEstimator * fSignalEstimator;
        IFluxEstimator * fFlux;
        IUnfoldEstimator * fUnfold;

        double fNTargets;
    };

    class EigenDifferentialCrossSectionEstimator : public ICrossSection,
                                                   public virtual IMeasurement,
                                                   public IEigenEval {
    public:
        EigenDifferentialCrossSectionEstimator() = default;

        EigenDifferentialCrossSectionEstimator(IEfficiencyEstimator * efficiency,
                                               ISignalEstimator * signal_estimator,
                                               IFluxEstimator * flux,
                                               IUnfoldEstimator * unfold,
                                               double ntargets = 0);

        static std::unique_ptr<IMeasurement> LoadFrom(xsec::type::LoadFunction<IEfficiencyEstimator> load_efficiency,
                                                      xsec::type::LoadFunction<ISignalEstimator> load_signal,
                                                      xsec::type::LoadFunction<IFluxEstimator> load_flux,
                                                      xsec::type::LoadFunction<IUnfoldEstimator> load_unfold,
                                                      TDirectory * dir,
                                                      const std::string & subdir);

    private:
        virtual void _eval_impl(const Array & data, const Array & error,
                                ArrayRef result, ArrayRef rerror) const override;
    };

    class EigenCrossSectionEstimator : public ICrossSection,
                                       public virtual IMeasurement,
                                       public IEigenEval {
    public:
        EigenCrossSectionEstimator() = default;

        EigenCrossSectionEstimator(IEfficiencyEstimator * efficiency,
                                   ISignalEstimator * signal_estimator,
                                   IFluxEstimator * flux,
                                   IUnfoldEstimator * unfold,
                                   double ntargets = 0);

        static std::unique_ptr<IMeasurement> LoadFrom(xsec::type::LoadFunction<IEfficiencyEstimator> load_efficiency,
                                                      xsec::type::LoadFunction<ISignalEstimator> load_signal,
                                                      xsec::type::LoadFunction<IFluxEstimator> load_flux,
                                                      xsec::type::LoadFunction<IUnfoldEstimator> load_unfold,
                                                      TDirectory * dir,
                                                      const std::string & subdir);

    private:

        virtual void _eval_impl(const Array & data, const Array & error,
                                ArrayRef result, ArrayRef rerror) const override;
    };

    const Array CalculateCrossSection(const Array * unfolded_selected_signal,
                                      const Array * efficiency,
                                      const Array flux,
                                      const double ntargets,
                                      const bool is_differential);

}
