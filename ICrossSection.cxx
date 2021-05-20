#include "XSecAna/ICrossSection.h"

#include "XSecAna/SimpleSignalEstimator.h"
#include "XSecAna/SimpleEfficiency.h"
#include "XSecAna/SimpleFlux.h"


namespace xsec {
  // make sure the interface compiles
  template class ICrossSection<SimpleSignalEstimator<HistXXd>,
			       SimpleEfficiency<HistXXd>,
			       SimpleFlux<HistXXd> >;

  ////////////////////////////////////////////////
  template<class HistType>
  HistType * CalculateCrossSection(HistType * unfolded_selected_signal,
				   HistType * efficiency,
				   HistType * flux,
				   double ntargets,
				   bool is_differential)
  {
    HistType * xsec = new HistType(*unfolded_selected_signal);
    xsec /= efficiency;
    xsec /= flux;
    xsec *= 1e4/ntargets; // Convert nu/m^2 to nu/cm^2
    if(is_differential) xsec->Normalize("width");
    return xsec;
  }

  ////////////////////////////////////////////////
  template<class SignalEstimatorType, 
	   class EfficiencyType, 
	   class FluxType,
	   bool IsDifferential,
	   class HistType>
  template<class UnfoldType>
  HistType *
  ICrossSection<SignalEstimatorType, 
	       EfficiencyType, 
	       FluxType,
	       IsDifferential,
	       HistType>::
  UnfoldedCrossSection(const HistType * data, UnfoldType unfold, double ntargets)
  {
    return CalculateCrossSection(unfold->Truth(fSignalEstimator->Signal(data)),
				 fEfficiency->Efficiency(),
				 fFlux->Flux(),
				 ntargets,
				 IsDifferential);
  }

  ////////////////////////////////////////////////
  template<class SignalEstimatorType, 
	   class EfficiencyType, 
	   class FluxType,
	   bool IsDifferential,
	   class HistType>
  void
  ICrossSection<SignalEstimatorType, 
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

    delete dir;
    tmp->cd();
  }

  ////////////////////////////////////////////////
  template<class SignalEstimatorType, 
	   class EfficiencyType, 
	   class FluxType,
	   bool IsDifferential,
	   class HistType>
  std::unique_ptr<ICrossSection<SignalEstimatorType, 
			       EfficiencyType, 
			       FluxType,
			       IsDifferential,
			       HistType> >
  ICrossSection<SignalEstimatorType, 
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
    
    if(dir->GetDirectory("fEfficiency"     )) eff  = EfficiencyType::     LoadFrom(dir, "fEfficiency"     ).release();
    if(dir->GetDirectory("fFlux"           )) flux = FluxType::           LoadFrom(dir, "fFlux"           ).release();
    if(dir->GetDirectory("fSignalEstimator")) sig  = SignalEstimatorType::LoadFrom(dir, "fSignalEstimator").release();

    return std::make_unique<ICrossSection<SignalEstimatorType, 
					  EfficiencyType, 
					  FluxType, 
					  IsDifferential,
					  HistType> >(eff, sig, flux);
  }
}
