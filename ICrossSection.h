#pragma once

#include "TDirectory.h"
#include "Hist.h"

namespace xsec {
  template<class SignalEstimatorType, 
	   class UnfoldType,
	   class EfficiencyType,
	   class FluxType,
	   bool IsDifferential = false,
	   class HistType = HistXd>
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

    HistType * UnfoldedCrossSection(const HistType & data, double ntargets);

    HistType * CrossSection(const HistType & data, double ntargets);
    
  public:
    EfficiencyType      * fEfficiency;
    SignalEstimatorType * fSignalEstimator;
    FluxType            * fFlux;    
    UnfoldType          * fUnfold;

  };

  ////////////////////////////////////////////////
  template<typename Scalar, int Cols>
  Hist<Scalar, Cols> * CalculateCrossSection(const Hist<Scalar, Cols> * unfolded_selected_signal,
					     const Hist<Scalar, Cols> * efficiency,
					     const Hist<Scalar, Cols> * flux,
					     const Scalar ntargets,
					     const bool is_differential)
  {
    Hist<Scalar, Cols> * xsec = new Hist<Scalar, Cols>(*unfolded_selected_signal);
    *xsec /= *efficiency;
    *xsec /= *flux;
    *xsec *= 1e4/ntargets; // Convert nu/m^2 to nu/cm^2
    if(is_differential) xsec->Normalize("width");
    return xsec;
  }

  ////////////////////////////////////////////////
  template<class SignalEstimatorType, 
	   class UnfoldType,
	   class EfficiencyType, 
	   class FluxType,
	   bool IsDifferential,
	   class HistType>
  HistType *
  ICrossSection<SignalEstimatorType, 
		UnfoldType,
		EfficiencyType, 
		FluxType,
		IsDifferential,
		HistType>::
  CrossSection(const HistType & data, double ntargets)
  {
    // invoke FluxType::operator* in case we want the integrated flux
    // Pass a histogram of ones through the flux parameter of
    // CalculateCrossSection to divide by one
    return CalculateCrossSection(fSignalEstimator->Signal(data),
				 (fFlux * fEfficiency->ToHist()),
				 std::decay_t<decltype(fEfficiency.ToHist().Contents())>
				 ::Ones(fEfficiency.ToHist().Contents().size()),
				 ntargets,
				 IsDifferential);
  }

  ////////////////////////////////////////////////
  template<class SignalEstimatorType, 
	   class UnfoldType,
	   class EfficiencyType, 
	   class FluxType,
	   bool IsDifferential,
	   class HistType>
  HistType *
  ICrossSection<SignalEstimatorType, 
		UnfoldType,
		EfficiencyType, 
		FluxType,
		IsDifferential,
		HistType>::
  UnfoldedCrossSection(const HistType & data, double ntargets)
  {
    // invoke FluxType::operator* in case we want the integrated flux
    // Pass a histogram of ones through the flux parameter of
    // CalculateCrossSection to divide by one
    return CalculateCrossSection(fUnfold->Truth(fSignalEstimator->Signal(data)),
				 (fFlux * fEfficiency->ToHist()),
				 std::decay_t<decltype(fEfficiency.ToHist().Contents())>
				 ::Ones(fEfficiency.ToHist().Contents().size()),
				 ntargets,
				 IsDifferential);
  }

  ////////////////////////////////////////////////
  template<class SignalEstimatorType, 
	   class UnfoldType,
	   class EfficiencyType, 
	   class FluxType,
	   bool IsDifferential,
	   class HistType>
  void
  ICrossSection<SignalEstimatorType, 
		UnfoldType,
		EfficiencyType, 
		FluxType,
		IsDifferential,
		HistType>::
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
	   bool IsDifferential,
	   class HistType>
  std::unique_ptr<ICrossSection<SignalEstimatorType, 
				UnfoldType,
				EfficiencyType, 
				FluxType,
				IsDifferential,
				HistType> >
  ICrossSection<SignalEstimatorType, 
		UnfoldType,
		EfficiencyType, 
		FluxType,
		IsDifferential,
		HistType>::
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
					  IsDifferential,
					  HistType> >(eff, sig, flux, unfold);
  }
}
