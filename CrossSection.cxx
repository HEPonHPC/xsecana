#include "XSecAna/CrossSection.h"

namespace xsec {
  ////////////////////////////////////////////////
  TH1 * CalculateCrossSection(TH1 * unfolded_selected_signal,
			      TH1 * efficiency,
			      TH1 * flux,
			      double ntargets,
			      bool is_differential)
  {
    TH1 * xsec = (TH1*) unfolded_selected_signal->Clone(ana::UniqueName().c_str());
    xsec->Divide(efficiency);
    xsec->Divide(flux);
    xsec->Scale(1e4/ntargets, is_differential? "width" : ""); // Convert nu/m^2 to nu/cm^2
    return xsec;
  }

  ////////////////////////////////////////////////
  template<class SignalEstimatorType, 
	   class EfficiencyType, 
	   class FluxType,
	   int IsDifferential>
  template<class UnfoldType,
	   int UnfoldReg>
  TH1 *
  CrossSection<SignalEstimatorType, 
	       EfficiencyType, 
	       FluxType,
	       IsDifferential>::
  UnfoldedCrossSection(const ana::Spectrum * data, UnfoldType unfold, double ntargets)
  {
    return CalculateCrossSection(unfold->Truth(fSignalEstimator->Signal(data)),
						 fEfficiency(),
						 fFlux(),
						 ntargets,
						 IsDifferential);
  }

  ////////////////////////////////////////////////
  template<class SignalEstimatorType, 
	   class EfficiencyType, 
	   class FluxType,
	   int IsDifferential>
  void
  CrossSection<SignalEstimatorType, 
	       EfficiencyType, 
	       FluxType,
	       IsDifferential>::
  SaveTo(TDirectory * dir, std::string subdir) const
  {
    TDirectory * tmp = gDirectory;
    dir = dir->mkdir(subdir.c_str());
    dir->cd();
    TObjString("CrossSection").Write("type");
    
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
	   int IsDifferential>
  std::unique_ptr<CrossSection<SignalEstimatorType, 
			       EfficiencyType, 
			       FluxType,
			       IsDifferential> >
  CrossSection<SignalEstimatorType, 
	       EfficiencyType, 
	       FluxType,
	       IsDifferential>::
  LoadFrom(TDirectory * dir, std::string subdir)    
  {
    dir = dir->GetDirectory(subdir.c_str());

    EfficiencyType *      eff  = 0;
    FluxType *            flux = 0;
    SignalEstimatorType * sig  = 0;
    
    if(dir->GetDirectory("fEfficiency"     )) eff  = EfficiencyType::     LoadFrom(dir, "fEfficiency"     ).release();
    if(dir->GetDirectory("fFlux"           )) flux = FluxType::           LoadFrom(dir, "fFlux"           ).release();
    if(dir->GetDirectory("fSignalEstimator")) sig  = SignalEstimatorType::LoadFrom(dir, "fSignalEstimator").release();

    return std::make_unique<CrossSection<SignalEstimatorType, EfficiencyType, FluxType, IsDifferential> >(eff, sig, flux);
  }
}
