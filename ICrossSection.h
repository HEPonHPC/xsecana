#pragma once

#include "CAFAna/Core/Spectrum.h"

#include "TH1.h"
#include "TDirectory.h"

namespace xsec {
  template<class SignalEstimatorType, 
	   class EfficiencyType,
	   class FluxType,
	   int IsDifferential>
  class CrossSection {

    CrossSection(EfficiencyType efficiency,
		 SignalEstimatorType signal_estimator,
		 FluxType flux)
      : fEfficiency(efficiency),
	fSignalEstimator(signal_estimator),
	fFlux(flux)
    {}
		    
    void SaveTo(TDirectory * dir, std::string subdir) const;
    std::unique_ptr<CrossSection> LoadFrom(TDirectory * dir, std::string subdir);

    template<class UnfoldType,
	     int UnfoldReg>
    TH1 * UnfoldedCrossSection(const ana::Spectrum * data, UnfoldType unfold, double ntargets);

  public:
    EfficiencyType      * fEfficiency;
    SignalEstimatorType * fSignalEstimator;
    FluxType            * fFlux;    
  };

  TH1 * CalculateCrossSection(TH1 * unfolded_selected_signal,
			      TH1 * efficiency,
			      TH1 * flux,
			      double ntargets,
			      bool is_differential);
}
