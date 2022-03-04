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

        ICrossSection(IMeasurement * efficiency,
                      IMeasurement * signal_estimator,
                      IMeasurement * flux,
                      IMeasurement * unfold,
                      double ntargets = 0);

        void SaveTo(TDirectory * dir, const std::string & subdir) const override;

        template<class DerivedCrossSection>
        static std::unique_ptr<IMeasurement> LoadFrom(xsec::type::LoadFunction<IMeasurement> load_efficiency,
                                                      xsec::type::LoadFunction<IMeasurement> load_signal,
                                                      xsec::type::LoadFunction<IMeasurement> load_flux,
                                                      xsec::type::LoadFunction<IMeasurement> load_unfold,
                                                      TDirectory * dir,
                                                      const std::string & subdir);

    protected:
        IMeasurement * fEfficiency;
        IMeasurement * fSignalEstimator;
        IMeasurement * fFlux;
        IMeasurement * fUnfold;

        double fNTargets;
    };

    class DifferentialCrossSectionEstimator : public ICrossSection,
                                              public virtual IMeasurement {
    public:
        DifferentialCrossSectionEstimator() = default;

        DifferentialCrossSectionEstimator(IMeasurement * efficiency,
                                          IMeasurement * signal_estimator,
                                          IMeasurement * flux,
                                          IMeasurement * unfold,
                                          double ntargets = 0);
        static std::unique_ptr<IMeasurement> LoadFrom(xsec::type::LoadFunction<IMeasurement> load_efficiency,
                                                      xsec::type::LoadFunction<IMeasurement> load_signal,
                                                      xsec::type::LoadFunction<IMeasurement> load_flux,
                                                      xsec::type::LoadFunction<IMeasurement> load_unfold,
                                                      TDirectory * dir,
                                                      const std::string & subdir);

        virtual std::shared_ptr<TH1> Eval(const TH1 * data) const override;

    };

    class CrossSectionEstimator : public ICrossSection,
                                  public virtual IMeasurement {
    public:
        CrossSectionEstimator() = default;

        CrossSectionEstimator(IMeasurement * efficiency,
                              IMeasurement * signal_estimator,
                              IMeasurement * flux,
                              IMeasurement * unfold,
                              double ntargets = 0);

        static std::unique_ptr<IMeasurement> LoadFrom(xsec::type::LoadFunction<IMeasurement> load_efficiency,
                                                      xsec::type::LoadFunction<IMeasurement> load_signal,
                                                      xsec::type::LoadFunction<IMeasurement> load_flux,
                                                      xsec::type::LoadFunction<IMeasurement> load_unfold,
                                                      TDirectory * dir,
                                                      const std::string & subdir);

        virtual std::shared_ptr<TH1> Eval(const TH1 * data) const override;

    private:

    };
    class EigenDifferentialCrossSectionEstimator : public ICrossSection,
                                                   public virtual IMeasurement,
                                                   public IEigenEval {
    public:
        EigenDifferentialCrossSectionEstimator() = default;

        EigenDifferentialCrossSectionEstimator(IMeasurement * efficiency,
                                               IMeasurement * signal_estimator,
                                               IMeasurement * flux,
                                               IMeasurement * unfold,
                                               double ntargets = 0);

        static std::unique_ptr<IMeasurement> LoadFrom(xsec::type::LoadFunction<IMeasurement> load_efficiency,
                                                      xsec::type::LoadFunction<IMeasurement> load_signal,
                                                      xsec::type::LoadFunction<IMeasurement> load_flux,
                                                      xsec::type::LoadFunction<IMeasurement> load_unfold,
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

        EigenCrossSectionEstimator(IMeasurement * efficiency,
                                   IMeasurement * signal_estimator,
                                   IMeasurement * flux,
                                   IMeasurement * unfold,
                                   double ntargets = 0);

        static std::unique_ptr<IMeasurement> LoadFrom(xsec::type::LoadFunction<IMeasurement> load_efficiency,
                                                      xsec::type::LoadFunction<IMeasurement> load_signal,
                                                      xsec::type::LoadFunction<IMeasurement> load_flux,
                                                      xsec::type::LoadFunction<IMeasurement> load_unfold,
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
    const std::shared_ptr<TH1>
    CalculateCrossSection(const TH1 * unfolded_selected_signal,
                          const TH1 * efficiency,
                          const TH1 * flux,
                          const double ntargets,
                          const bool is_differential);

}
