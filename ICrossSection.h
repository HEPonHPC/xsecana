#pragma once

#include "TDirectory.h"
#include "Hist.h"

namespace xsec {
  template<class SignalEstimatorType, 
	   class UnfoldType,
	   class EfficiencyType,
	   class FluxType,
	   bool IsDifferential = false>
  class ICrossSection {
  public:
    ICrossSection() {}

    ICrossSection(EfficiencyType * efficiency,
		  SignalEstimatorType * signal_estimator,
		  FluxType * flux,
		  UnfoldType * unfold)
      : fEfficiency(efficiency),
	fSignalEstimator(signal_estimator),
	fFlux(flux),
	fUnfold(unfold)
    {}
		    
    void SaveTo(TDirectory * dir, std::string subdir) const;
    std::unique_ptr<ICrossSection> LoadFrom(TDirectory * dir, std::string subdir);

    template<class HistType>
    HistType UnfoldedCrossSection(const HistType & data, typename HistType::scalar_type ntargets);

    template<class HistType>
    HistType CrossSection(const HistType & data, typename HistType::scalar_type ntargets);

    ICrossSection<SignalEstimatorType,
		  UnfoldType,
		  EfficiencyType,
		  FluxType,
		  true>
    ToDifferential();
    
  public:
    EfficiencyType      * fEfficiency;
    SignalEstimatorType * fSignalEstimator;
    FluxType            * fFlux;    
    UnfoldType          * fUnfold;

  };

  ////////////////////////////////////////////////
  template<class SignalEstimatorType, 
	   class UnfoldType,
	   class EfficiencyType, 
	   class FluxType,
	   bool IsDifferential>
  ICrossSection<SignalEstimatorType,
		UnfoldType,
		EfficiencyType,
		FluxType,
		true>
  ICrossSection<SignalEstimatorType,
		UnfoldType,
		EfficiencyType,
		FluxType,
		IsDifferential>::
  ToDifferential()
  {
    if constexpr(IsDifferential) {
	return *this;
      }
    else {
      return ICrossSection<SignalEstimatorType,
			   UnfoldType,
			   EfficiencyType,
			   FluxType,
			   true>(fEfficiency,
				 fSignalEstimator,
				 fFlux,
				 fUnfold);
    }
  }

  ////////////////////////////////////////////////
  template<class HistType>
  HistType CalculateCrossSection(const HistType & unfolded_selected_signal,
				 const HistType & efficiency,
				 const HistType & flux,
				 const typename HistType::scalar_type ntargets,
				 const bool is_differential)
  {
    HistType  xsec = unfolded_selected_signal;
    xsec /= efficiency;
    xsec /= flux;
    xsec *= 1e4/ntargets; // Convert nu/m^2 to nu/cm^2
    if(is_differential) xsec.Normalize("width");
    return xsec;
  }

  ////////////////////////////////////////////////
  template<class SignalEstimatorType, 
	   class UnfoldType,
	   class EfficiencyType, 
	   class FluxType,
	   bool IsDifferential>
  template<class HistType>
  HistType
  ICrossSection<SignalEstimatorType, 
		UnfoldType,
		EfficiencyType, 
		FluxType,
		IsDifferential>::
  CrossSection(const HistType & data, typename HistType::scalar_type ntargets)
  {
    // invoke FluxType::operator* in case we want the integrated flux
    // Pass a histogram of ones through the flux parameter of
    // CalculateCrossSection to divide by one
    return CalculateCrossSection(fSignalEstimator->Signal(data),
				 (*fFlux * fEfficiency->ToHist()),
				 HistType(HistType::array_type::Ones(fEfficiency->ToHist().Contents().size()),
					  fEfficiency->ToHist().Edges()),
				 ntargets,
				 IsDifferential);
  }

  ////////////////////////////////////////////////
  template<class SignalEstimatorType, 
	   class UnfoldType,
	   class EfficiencyType, 
	   class FluxType,
	   bool IsDifferential>
  template<class HistType>
  HistType
  ICrossSection<SignalEstimatorType, 
		UnfoldType,
		EfficiencyType, 
		FluxType,
		IsDifferential>::
  UnfoldedCrossSection(const HistType & data, typename HistType::scalar_type ntargets)
  {
    // invoke FluxType::operator* in case we want the integrated flux
    // Pass a histogram of ones through the flux parameter of
    // CalculateCrossSection to divide by one
    return CalculateCrossSection(fUnfold->Truth(fSignalEstimator->Signal(data)),
				 (*fFlux * fEfficiency->ToHist()),
				 HistType(HistType::array_type::Ones(fEfficiency->ToHist().Contents().size()),
					  fEfficiency->ToHist().Edges()),
				 ntargets,
				 IsDifferential);
  }

  ////////////////////////////////////////////////
  template<class SignalEstimatorType, 
	   class UnfoldType,
	   class EfficiencyType, 
	   class FluxType,
	   bool IsDifferential>
  void
  ICrossSection<SignalEstimatorType, 
		UnfoldType,
		EfficiencyType, 
		FluxType,
		IsDifferential>::
  SaveTo(TDirectory * dir, std::string subdir) const
  {
    TDirectory * tmp = gDirectory;
    dir = dir->mkdir(subdir.c_str());
    dir->cd();
    TObjString("ICrossSection").Write("type");
    
    if(fEfficiency     ) fEfficiency     ->SaveTo(dir, "fEfficiency"     );
    if(fFlux           ) fFlux           ->SaveTo(dir, "fFlux"           );
    if(fSignalEstimator) fSignalEstimator->SaveTo(dir, "fSignalEstimator");
    if(fUnfold         ) fUnfold         ->SaveTo(dir, "fUnfold"         );

    tmp->cd();
  }

  ////////////////////////////////////////////////
  template<class SignalEstimatorType, 
	   class UnfoldType,
	   class EfficiencyType, 
	   class FluxType,
	   bool IsDifferential>
  std::unique_ptr<ICrossSection<SignalEstimatorType, 
				UnfoldType,
				EfficiencyType, 
				FluxType,
				IsDifferential> >
  ICrossSection<SignalEstimatorType, 
		UnfoldType,
		EfficiencyType, 
		FluxType,
		IsDifferential>::
  LoadFrom(TDirectory * dir, std::string subdir)    
  {
    dir = dir->GetDirectory(subdir.c_str());

    EfficiencyType *      eff  = 0;
    FluxType *            flux = 0;
    SignalEstimatorType * sig  = 0;
    UnfoldType *          unfold = 0;
    
    if(dir->GetDirectory("fEfficiency"     )) eff  = EfficiencyType     ::LoadFrom(dir, "fEfficiency"     ).release();
    if(dir->GetDirectory("fFlux"           )) flux = FluxType           ::LoadFrom(dir, "fFlux"           ).release();
    if(dir->GetDirectory("fSignalEstimator")) sig  = SignalEstimatorType::LoadFrom(dir, "fSignalEstimator").release();
    if(dir->GetDirectory("fUnfold"         )) unfold  = UnfoldType      ::LoadFrom(dir, "fUnfold").release();

    return std::make_unique<ICrossSection<SignalEstimatorType, 
					  UnfoldType,
					  EfficiencyType, 
					  FluxType, 
					  IsDifferential> >
      (eff, sig, flux, unfold);
  }
}
