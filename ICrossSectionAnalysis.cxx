#include "XSecAna/ICrossSectionAnalysis.h"

// cafana includes
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Unfold/UnfoldIterative.h"
#include "CAFAna/Unfold/UnfoldTikhonov.h"
#include "CAFAna/Unfold/UnfoldSVD.h"
#include "CAFAna/Unfold/UnfoldMaxEnt.h"

// root includes
#include "TH2.h"
#include "TH3.h"

#include <iostream>

namespace xsec {

  ///////////////////////////////////////////////////////////////////////
  template<class CrossSectionType,
	   class UnfoldType,
	   class UncertaintyPropogator>
  TH1 * 
  ICrossSectionAnalysis<CrossSectionType,
			UnfoldType,
			UncertaintyPropogator>::
  UnfoldedCrossSection(std::string syst_name)
  {    
    if(syst_name == "nominal") {
      if(!fUnfoldedNominalXSec)  {
	fUnfoldedNominalXSec = fNominalXSec->UnfoldedCrossSection(fData, 
								  fUnfold, 
								  this->NTargets());
      }
      return fUnfoldedNominalXSec;
    }
    else {
      if(fUnfoldedShiftedXSec.find(syst_name) == fUnfoldedShiftedXSec.end()) {
	fUnfoldedShiftedXSec[syst_name] = 
	  fShiftedXSec.at(syst_name)->Invoke(CrossSectionType::UnfoldedCrossSection,
					     fData,
					     fUnfold,
					     this->NTargets());
									     
      }
    }
  }

  ///////////////////////////////////////////////////////////////////////
  template<class CrossSectionType,
	   class UnfoldType,
	   class UncertaintyPropogator>
  void
  ICrossSectionAnalysis<CrossSectionType,
			UnfoldType,
			UncertaintyPropogator>::
  SaveTo(TDirectory * dir, std::string subdir) const
  {
    TDirectory * tmp = gDirectory;

    dir = dir->mkdir(subdir.c_str()); // switch to subdir
    dir->cd();
    TObjString("ICrossSectionAnalysis").Write("type");

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
	   class UncertaintyPropogator>
  std::unique_ptr<ICrossSectionAnalysis<CrossSectionType,
					UnfoldType, 
					UncertaintyPropogator> > 
  ICrossSectionAnalysis<CrossSectionType,
			UnfoldType,
			UncertaintyPropogator>::
  LoadFrom(TDirectory * dir, std::string subdir)
  {
    TDirectory * tmp = gDirectory;
    dir = dir->GetDirectory(subdir.c_str());

    TObjString * ptag = (TObjString*) dir->Get("type");
    assert(ptag->GetString() == "ICrossSectionAnalysis" && "Type does not match ICrossSectionAnalysis");
    delete ptag;

    auto unfold = UnfoldType::LoadFrom(dir, "fUnfold");

    auto data = ana::Spectrum::LoadFrom(dir, "fData");

    auto nominal_xsec = CrossSectionType::LoadFrom(dir, "fNominalXSec");
    std::map<std::string, Systematic<CrossSectionType> > shifted_xsec;
    auto syst_dir = dir->GetDirectory("fShiftedXSec");
    for(auto syst_name : *syst_dir->GetListOfKeys()) {
      shifted_xsec[syst_name->GetName()] = Systematic<CrossSectionType>::LoadFrom(syst_dir, 
										  syst_name->GetName()).release();
    }

    return std::make_unique<ICrossSectionAnalysis<CrossSectionType,
						  UnfoldType, 
						  UncertaintyPropogator> > 
      (nominal_xsec,
       shifted_xsec,
       unfold,
       data);
    
  }


  ///////////////////////////////////////////////////////////////////////
  template<class CrossSectionType,
	   class UnfoldType,
	   class UncertaintyPropogator>
  ICrossSectionAnalysis<CrossSectionType,
			UnfoldType,
			UncertaintyPropogator>::
  ~ICrossSectionAnalysis()
  {
    delete fSelectionCut;
    delete fSignalCut;
    delete fRecoAxis;
    delete fTruthAxis;
  }


}
