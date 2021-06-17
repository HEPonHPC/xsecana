#pragma once

#include <string>
#include <map>

#include "XSecAna/Systematic.h"
#include "XSecAna/Hist.h"

// root includes
#include "TDirectory.h"

namespace xsec {
  template<class CrossSectionType,
	   class UncertaintyPropagator,
	   class HistType = HistXd>
  class CrossSectionAnalysis {

  public:

    CrossSectionAnalysis() {}
    
    /// \brief Forward to UncertaintyPropagator
    HistType * AbsoluteUncertaintyXSec(std::string syst_name, double ntargets);
    HistType * AbsoluteUncertaintyUnfoldedXSec(std::string syst_name, double ntargets);

    /// \brief Forward to UncertaintyPropagator
    HistType * FractionalUncertaintyXSec(std::string syst_name, double ntargets);
    HistType * FractionalUncertaintyUnfoldedXSec(std::string syst_name, double ntargets);

    /// \brief Forward to UncertaintyPropagator
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

    ~CrossSectionAnalysis();

  protected:
    ///\brief constructor for loading analysis from file
    CrossSectionAnalysis(CrossSectionType nominal_xsec,
			 std::map<std::string, Systematic<CrossSectionType> > shifted_xsec,
			 HistType data)
      : fNominalXSec(nominal_xsec),
	fShiftedXSec(shifted_xsec),
	fData(data)
    {}

    // nominal is special 
    CrossSectionType fNominalXSec;

    std::map<std::string, Systematic<CrossSectionType> > fShiftedXSec;

    const HistType fData;

    // cache the nominal unfolded results
    HistType * fUnfoldedNominalXSec = 0;

    // cache the unfolded shifted results
    std::map<std::string, Systematic<HistType> > fUnfoldedShiftedXSec;

    // 
    UncertaintyPropagator * fUncertaintyPropagator = 0;
  };

  ///////////////////////////////////////////////////////////////////////
  template<class CrossSectionType,
	   class UncertaintyPropagator,
	   class HistType>
  HistType * 
  CrossSectionAnalysis<CrossSectionType,
		       UncertaintyPropagator,
		       HistType>::
  AbsoluteUncertaintyUnfoldedXSec(std::string syst_name, double ntargets) 
  {
    return fUncertaintyPropagator->AbsoluteUncertaintyUnfolded(fData,
							       fNominalXSec,
							       fShiftedXSec.at(syst_name),
							       ntargets);
  }


  ///////////////////////////////////////////////////////////////////////
  template<class CrossSectionType,
	   class UncertaintyPropagator,
	   class HistType>
  HistType * 
  CrossSectionAnalysis<CrossSectionType,
		       UncertaintyPropagator,
		       HistType>::
  AbsoluteUncertaintyXSec(std::string syst_name, double ntargets)
  {
    return fUncertaintyPropagator->AbsoluteUncertaintyXSec(fData,
							   fNominalXSec,
							   fShiftedXSec.at(syst_name),
							   ntargets);
  }

  ///////////////////////////////////////////////////////////////////////
  template<class CrossSectionType,
	   class UncertaintyPropagator,
	   class HistType>
  HistType * 
  CrossSectionAnalysis<CrossSectionType,
		       UncertaintyPropagator,
		       HistType>::
  FractionalUncertaintyUnfoldedXSec(std::string syst_name, double ntargets)
  {
    return fUncertaintyPropagator->FractionalUncertaintyUnfoldedXSec(fData,
								     fNominalXSec,
								     fShiftedXSec.at(syst_name),
								     ntargets);
  }

  ///////////////////////////////////////////////////////////////////////
  template<class CrossSectionType,
	   class UncertaintyPropagator,
	   class HistType>
  HistType * 
  CrossSectionAnalysis<CrossSectionType,
		       UncertaintyPropagator,
		       HistType>::
  FractionalUncertaintyXSec(std::string syst_name, double ntargets)
  {
    return fUncertaintyPropagator->FractionalUncertaintyXSec(fData,
							     fNominalXSec,
							     fShiftedXSec.at(syst_name),
							     ntargets);
  }

  ///////////////////////////////////////////////////////////////////////
  template<class CrossSectionType,
	   class UncertaintyPropagator,
	   class HistType>
  std::pair<HistType*,HistType*> 
  CrossSectionAnalysis<CrossSectionType,
		       UncertaintyPropagator,
		       HistType>::
  TotalAbsoluteUncertaintyUnfoldedXSec(double ntargets)
  {
    return fUncertaintyPropagator->TotalAbsoluteUncertaintyUnfoldedXSec(fData,
									fNominalXSec,
									fShiftedXSec,
									ntargets);
  }

  ///////////////////////////////////////////////////////////////////////
  template<class CrossSectionType,
	   class UncertaintyPropagator,
	   class HistType>
  std::pair<HistType*,HistType*> 
  CrossSectionAnalysis<CrossSectionType,
		       UncertaintyPropagator,
		       HistType>::
  TotalAbsoluteUncertaintyXSec(double ntargets)
  {
    return fUncertaintyPropagator->TotalAbsoluteUncertaintyXSec(fData,
								fNominalXSec,
								fShiftedXSec,
								ntargets);
  }

  ///////////////////////////////////////////////////////////////////////
  template<class CrossSectionType,
	   class UncertaintyPropagator,
	   class HistType>
  std::pair<HistType*,HistType*> 
  CrossSectionAnalysis<CrossSectionType,
		       UncertaintyPropagator,
		       HistType>::
  TotalFractionalUncertaintyUnfoldedXSec(double ntargets)
  {
    return fUncertaintyPropagator->TotalFractionalUncertaintyUnfoldedXSec(fData,
									  fNominalXSec,
									  fShiftedXSec,
									  ntargets);
  }

  ///////////////////////////////////////////////////////////////////////
  template<class CrossSectionType,
	   class UncertaintyPropagator,
	   class HistType>
  std::pair<HistType*,HistType*> 
  CrossSectionAnalysis<CrossSectionType,
		       UncertaintyPropagator,
		       HistType>::
  TotalFractionalUncertaintyXSec(double ntargets)
  {
    return fUncertaintyPropagator->TotalFractionalUncertaintyXSec(fData,
								  fNominalXSec,
								  fShiftedXSec,
								  ntargets);
  }
  

  ///////////////////////////////////////////////////////////////////////
  template<class CrossSectionType,
	   class UncertaintyPropagator,
	   class HistType>
  const HistType * 
  CrossSectionAnalysis<CrossSectionType,
		       UncertaintyPropagator,
		       HistType>::
  UnfoldedCrossSection(std::string syst_name,
		       double ntargets)
  {    
    if(syst_name == "nominal") {
      if(!fUnfoldedNominalXSec)  {
	fUnfoldedNominalXSec = fNominalXSec->UnfoldedCrossSection(fData, 
								  ntargets);
      }
      return fUnfoldedNominalXSec;
    }
    else {
      if(fUnfoldedShiftedXSec.find(syst_name) == fUnfoldedShiftedXSec.end()) {
	fUnfoldedShiftedXSec[syst_name] = 
	  fShiftedXSec.at(syst_name)->Invoke(&CrossSectionType::UnfoldedCrossSection,
					     fData,
					     ntargets);
									     
      }
    }
  }

  ///////////////////////////////////////////////////////////////////////
  template<class CrossSectionType,
	   class UncertaintyPropagator,
	   class HistType>
  void
  CrossSectionAnalysis<CrossSectionType,
		       UncertaintyPropagator,
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

    fData.SaveTo(dir, "fData");

    delete dir;
    tmp->cd();
  }

  ///////////////////////////////////////////////////////////////////////
  template<class CrossSectionType,
	   class UncertaintyPropagator,
	   class HistType>
  std::unique_ptr<CrossSectionAnalysis<CrossSectionType,
				       UncertaintyPropagator,
				       HistType> > 
  CrossSectionAnalysis<CrossSectionType,
		       UncertaintyPropagator,
		       HistType>::
  LoadFrom(TDirectory * dir, std::string subdir)
  {
    TDirectory * tmp = gDirectory;
    dir = dir->GetDirectory(subdir.c_str());

    TObjString * ptag = (TObjString*) dir->Get("type");
    assert(ptag->GetString() == "CrossSectionAnalysis" && "Type does not match CrossSectionAnalysis");
    delete ptag;

    auto data = HistType::LoadFrom(dir, "fData");

    auto nominal_xsec = CrossSectionType::LoadFrom(dir, "fNominalXSec");
    std::map<std::string, Systematic<CrossSectionType> > shifted_xsec;
    auto syst_dir = dir->GetDirectory("fShiftedXSec");
    for(auto syst_name : *syst_dir->GetListOfKeys()) {
      shifted_xsec[syst_name->GetName()] = Systematic<CrossSectionType>::LoadFrom(syst_dir, 
										  syst_name->GetName()).release();
    }

    return std::make_unique<CrossSectionAnalysis<CrossSectionType,
						 UncertaintyPropagator,
						 HistType> > 
      (*nominal_xsec,
       shifted_xsec,
       *data);
    
  }

  ///////////////////////////////////////////////////////////////////////
  template<class CrossSectionType,
	   class UncertaintyPropagator,
	   class HistType>
  CrossSectionAnalysis<CrossSectionType,
		       UncertaintyPropagator,
		       HistType>::
  ~CrossSectionAnalysis()
  {
    if(fUnfoldedNominalXSec) delete fUnfoldedNominalXSec;
  }

}
