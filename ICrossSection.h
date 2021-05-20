#pragma once

#include "TDirectory.h"
#include "Hist.h"

namespace xsec {
  template<class SignalEstimatorType, 
	   class EfficiencyType,
	   class FluxType,
	   bool IsDifferential = false,
	   class HistType = HistXXd>
  class ICrossSection {
  public:
    ICrossSection(EfficiencyType * efficiency,
		 SignalEstimatorType * signal_estimator,
		 FluxType * flux)
      : fEfficiency(efficiency),
	fSignalEstimator(signal_estimator),
	fFlux(flux)
    {}
		    
    void SaveTo(TDirectory * dir, std::string subdir) const;
    std::unique_ptr<ICrossSection> LoadFrom(TDirectory * dir, std::string subdir);

    template<class UnfoldType>
    HistType * UnfoldedCrossSection(const HistType * data, UnfoldType unfold, double ntargets);

  public:
    EfficiencyType      * fEfficiency;
    SignalEstimatorType * fSignalEstimator;
    FluxType            * fFlux;    
  };

  template<class HistType>
  HistType * CalculateCrossSection(HistType * unfolded_selected_signal,
				   HistType * efficiency,
				   HistType * flux,
				   double ntargets,
				   bool is_differential);
}
