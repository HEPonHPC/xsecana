#include "XSecAna/CrossSectionAnalysis.h"

// root includes
#include "TDirectory.h"

#include <iostream>

namespace xsec {
  ///////////////////////////////////////////////////////////////////////
  template<class CrossSectionType,
	   class UnfoldType,
	   class UncertaintyPropogator,
	   class HistType>
  HistType * 
  CrossSectionAnalysis<CrossSectionType,
		       UnfoldType,
		       UncertaintyPropogator,
		       HistType>::
  AbsoluteUncertaintyUnfoldedXSec(std::string syst_name, double ntargets) 
  {
    return UncertaintyPropogator::Instance().AbsoluteUncertaintyUnfolded(fUnfold,
									 fNominalXSec,
									 fShiftedXSec.at(syst_name),
									 ntargets);
  }

  ///////////////////////////////////////////////////////////////////////
  template<class CrossSectionType,
	   class UnfoldType,
	   class UncertaintyPropogator,
	   class HistType>
  HistType * 
  CrossSectionAnalysis<CrossSectionType,
		       UnfoldType,
		       UncertaintyPropogator,
		       HistType>::
  AbsoluteUncertaintyXSec(std::string syst_name, double ntargets)
  {
    return UncertaintyPropogator::Instance().AbsoluteUncertaintyXSec(fNominalXSec,
								     fShiftedXSec.at(syst_name),
								     ntargets);
  }

  ///////////////////////////////////////////////////////////////////////
  template<class CrossSectionType,
	   class UnfoldType,
	   class UncertaintyPropogator,
	   class HistType>
  HistType * 
  CrossSectionAnalysis<CrossSectionType,
		       UnfoldType,
		       UncertaintyPropogator,
		       HistType>::
  FractionalUncertaintyUnfoldedXSec(std::string syst_name, double ntargets)
  {
    return UncertaintyPropogator::Instance().FractionalUncertaintyUnfoldedXSec(fUnfold,
									       fNominalXSec,
									       fShiftedXSec.at(syst_name),
									       ntargets);
  }

  ///////////////////////////////////////////////////////////////////////
  template<class CrossSectionType,
	   class UnfoldType,
	   class UncertaintyPropogator,
	   class HistType>
  HistType * 
  CrossSectionAnalysis<CrossSectionType,
		       UnfoldType,
		       UncertaintyPropogator,
		       HistType>::
  FractionalUncertaintyXSec(std::string syst_name, double ntargets)
  {
    return UncertaintyPropogator::Instance().FractionalUncertaintyXSec(fNominalXSec,
								       fShiftedXSec.at(syst_name),
								       ntargets);
  }

  ///////////////////////////////////////////////////////////////////////
  template<class CrossSectionType,
	   class UnfoldType,
	   class UncertaintyPropogator,
	   class HistType>
  std::pair<HistType*,HistType*> 
  CrossSectionAnalysis<CrossSectionType,
		       UnfoldType,
		       UncertaintyPropogator,
		       HistType>::
  TotalAbsoluteUncertaintyUnfoldedXSec(double ntargets)
  {
    return UncertaintyPropogator::Instance().TotalAbsoluteUncertaintyUnfoldedXSec(fUnfold,
										  fNominalXSec,
										  fShiftedXSec,
										  ntargets);
  }

  ///////////////////////////////////////////////////////////////////////
  template<class CrossSectionType,
	   class UnfoldType,
	   class UncertaintyPropogator,
	   class HistType>
  std::pair<HistType*,HistType*> 
  CrossSectionAnalysis<CrossSectionType,
		       UnfoldType,
		       UncertaintyPropogator,
		       HistType>::
  TotalAbsoluteUncertaintyXSec(double ntargets)
  {
    return UncertaintyPropogator::Instance().TotalAbsoluteUncertaintyXSec(fNominalXSec,
									  fShiftedXSec,
									  ntargets);
  }

  ///////////////////////////////////////////////////////////////////////
  template<class CrossSectionType,
	   class UnfoldType,
	   class UncertaintyPropogator,
	   class HistType>
  std::pair<HistType*,HistType*> 
  CrossSectionAnalysis<CrossSectionType,
		       UnfoldType,
		       UncertaintyPropogator,
		       HistType>::
  TotalFractionalUncertaintyUnfoldedXSec(double ntargets)
  {
    return UncertaintyPropogator::Instance().TotalFractionalUncertaintyUnfoldedXSec(fUnfold,
										    fNominalXSec,
										    fShiftedXSec,
										    ntargets);
  }

  ///////////////////////////////////////////////////////////////////////
  template<class CrossSectionType,
	   class UnfoldType,
	   class UncertaintyPropogator,
	   class HistType>
  std::pair<HistType*,HistType*> 
  CrossSectionAnalysis<CrossSectionType,
		       UnfoldType,
		       UncertaintyPropogator,
		       HistType>::
  TotalFractionalUncertaintyXSec(double ntargets)
  {
    return UncertaintyPropogator::Instance().TotalFractionalUncertaintyXSec(fNominalXSec,
									    fShiftedXSec,
									    ntargets);
  }
  

  ///////////////////////////////////////////////////////////////////////
  template<class CrossSectionType,
	   class UnfoldType,
	   class UncertaintyPropogator,
	   class HistType>
  const HistType * 
  CrossSectionAnalysis<CrossSectionType,
		       UnfoldType,
		       UncertaintyPropogator,
		       HistType>::
  UnfoldedCrossSection(std::string syst_name,
		       double ntargets)
  {    
    if(syst_name == "nominal") {
      if(!fUnfoldedNominalXSec)  {
	fUnfoldedNominalXSec = fNominalXSec->UnfoldedCrossSection(fData, 
								  fUnfold, 
								  ntargets);
      }
      return fUnfoldedNominalXSec;
    }
    else {
      if(fUnfoldedShiftedXSec.find(syst_name) == fUnfoldedShiftedXSec.end()) {
	fUnfoldedShiftedXSec[syst_name] = 
	  fShiftedXSec.at(syst_name)->Invoke(&CrossSectionType::UnfoldedCrossSection,
					     fData,
					     fUnfold,
					     ntargets);
									     
      }
    }
  }

  ///////////////////////////////////////////////////////////////////////
  template<class CrossSectionType,
	   class UnfoldType,
	   class UncertaintyPropogator,
	   class HistType>
  void
  CrossSectionAnalysis<CrossSectionType,
		       UnfoldType,
		       UncertaintyPropogator,
		       HistType>::
  SaveTo(TDirectory * dir, std::string subdir) const
  {
    TDirectory * tmp = gDirectory;

    dir = dir->mkdir(subdir.c_str()); // switch to subdir
    dir->cd();
    TObjString("CrossSectionAnalysis").Write("type");

    fNominalXSec->SaveTo(dir, "fNominalXSec");
    dir->cd();

    auto syst_dir = dir->mkdir("fShiftedXSec");
    for(auto shifted : fShiftedXSec) {
      shifted->SaveTo(syst_dir, shifted->GetName());
    }

    fUnfold->Save(dir, "fUnfold");
    
    fData->SaveTo(dir, "fData");

    delete dir;
    tmp->cd();
  }

  ///////////////////////////////////////////////////////////////////////
  template<class CrossSectionType,
	   class UnfoldType,
	   class UncertaintyPropogator,
	   class HistType>
  std::unique_ptr<CrossSectionAnalysis<CrossSectionType,
				       UnfoldType, 
				       UncertaintyPropogator,
				       HistType> > 
  CrossSectionAnalysis<CrossSectionType,
		       UnfoldType,
		       UncertaintyPropogator,
		       HistType>::
  LoadFrom(TDirectory * dir, std::string subdir)
  {
    TDirectory * tmp = gDirectory;
    dir = dir->GetDirectory(subdir.c_str());

    TObjString * ptag = (TObjString*) dir->Get("type");
    assert(ptag->GetString() == "CrossSectionAnalysis" && "Type does not match CrossSectionAnalysis");
    delete ptag;

    auto unfold = UnfoldType::LoadFrom(dir, "fUnfold");

    auto data = HistType::LoadFrom(dir, "fData");

    auto nominal_xsec = CrossSectionType::LoadFrom(dir, "fNominalXSec");
    std::map<std::string, Systematic<CrossSectionType> > shifted_xsec;
    auto syst_dir = dir->GetDirectory("fShiftedXSec");
    for(auto syst_name : *syst_dir->GetListOfKeys()) {
      shifted_xsec[syst_name->GetName()] = Systematic<CrossSectionType>::LoadFrom(syst_dir, 
										  syst_name->GetName()).release();
    }

    return std::make_unique<CrossSectionAnalysis<CrossSectionType,
						 UnfoldType, 
						 UncertaintyPropogator,
						 HistType> > 
      (nominal_xsec,
       shifted_xsec,
       unfold,
       data);
    
  }


  ///////////////////////////////////////////////////////////////////////
  template<class CrossSectionType,
	   class UnfoldType,
	   class UncertaintyPropogator,
	   class HistType>
  CrossSectionAnalysis<CrossSectionType,
		       UnfoldType,
		       UncertaintyPropogator,
		       HistType>::
  ~CrossSectionAnalysis()
  {
    delete fUnfold;
    delete fData;
    delete fNominalXSec;
    for(auto map_it = fUnfoldedShiftedXSec.begin(); map_it != fUnfoldedShiftedXSec.end(); map_it++) {
      delete map_it->second;
    }
    for(auto map_it = fShiftedXSec.begin(); map_it != fShiftedXSec.end(); map_it++) {
      delete map_it->second;
    }
    fUnfoldedShiftedXSec.clear();
  }


}
