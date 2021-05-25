#pragma once

#include <string>
#include <map>

#include "XSecAna/Systematic.h"
#include "XSecAna/Hist.h"

// root includes
#include "TDirectory.h"

namespace xsec {
  template<class CrossSectionType,
	   class UnfoldType,
	   class UncertaintyPropogator,
	   class HistType = HistXd>
  class CrossSectionAnalysis {

  public:

    CrossSectionAnalysis() {}
    
    /// \brief Forward to UncertaintyPropogator
    HistType * AbsoluteUncertaintyXSec(std::string syst_name, double ntargets);
    HistType * AbsoluteUncertaintyUnfoldedXSec(std::string syst_name, double ntargets);

    /// \brief Forward to UncertaintyPropogator
    HistType * FractionalUncertaintyXSec(std::string syst_name, double ntargets);
    HistType * FractionalUncertaintyUnfoldedXSec(std::string syst_name, double ntargets);

    /// \brief Forward to UncertaintyPropogator
    std::pair<HistType*, HistType*> TotalAbsoluteUncertaintyXSec          (double ntargets);
    std::pair<HistType*, HistType*> TotalAbsoluteUncertaintyUnfoldedXSec  (double ntargets);
    std::pair<HistType*, HistType*> TotalFractionalUncertaintyXSec        (double ntargets);
    std::pair<HistType*, HistType*> TotalFractionalUncertaintyUnfoldedXSec(double ntargets);

    /// \brief Return an unfolded cross section result for the input systematic 
    /// pass syst_name = "nominal" for the nominal result
    const HistType * UnfoldedCrossSection(std::string syst_name, double ntargets);

    //    ~CrossSectionAnalysis();

    void SaveTo(TDirectory * dir, std::string subdir) const;
    static std::unique_ptr<CrossSectionAnalysis> LoadFrom(TDirectory * dir, std::string name);

  protected:
    ///\brief constructor for loading analysis from file
    CrossSectionAnalysis(CrossSectionType * nominal_xsec,
			 std::map<std::string, Systematic<CrossSectionType> > shifted_xsec,
			 UnfoldType * unfold,
			 HistType * data)
      : fNominalXSec(nominal_xsec),
	fShiftedXSec(shifted_xsec),
	fUnfold(unfold),
	fData(data)
    {}

    // nominal is special 
    CrossSectionType * fNominalXSec;

    std::map<std::string, Systematic<CrossSectionType> > fShiftedXSec;

    UnfoldType * fUnfold = NULL;
    const HistType * fData  = NULL;

    // cache the nominal unfolded results
    HistType * fUnfoldedNominalXSec;

    // cache the unfolded shifted results
    std::map<std::string, Systematic<HistType> > fUnfoldedShiftedXSec;

    // 
    UncertaintyPropogator * fUncertaintyPropogator = 0;
  };

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
    return fUncertaintyPropogator->AbsoluteUncertaintyUnfolded(fData,
							       fNominalXSec,
							       fShiftedXSec.at(syst_name),
							       fUnfold,
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
    return fUncertaintyPropogator->AbsoluteUncertaintyXSec(fData,
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
  FractionalUncertaintyUnfoldedXSec(std::string syst_name, double ntargets)
  {
    return fUncertaintyPropogator->FractionalUncertaintyUnfoldedXSec(fData,
								     fNominalXSec,
								     fShiftedXSec.at(syst_name),
								     fUnfold,
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
    return fUncertaintyPropogator->FractionalUncertaintyXSec(fData,
							     fNominalXSec,
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
    return fUncertaintyPropogator->TotalAbsoluteUncertaintyUnfoldedXSec(fData,
									fNominalXSec,
									fShiftedXSec,
									fUnfold,
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
    return fUncertaintyPropogator->TotalAbsoluteUncertaintyXSec(fData,
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
  TotalFractionalUncertaintyUnfoldedXSec(double ntargets)
  {
    return fUncertaintyPropogator->TotalFractionalUncertaintyUnfoldedXSec(fData,
									  fNominalXSec,
									  fShiftedXSec,
									  fUnfold,
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
    return fUncertaintyPropogator->TotalFractionalUncertaintyXSec(fData,
								  fNominalXSec,
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

  /*
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
  */
}
