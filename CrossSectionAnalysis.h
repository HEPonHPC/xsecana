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
    CrossSectionAnalysis(CrossSectionType nominal_xsec,
			 std::map<std::string, Systematic<CrossSectionType> > shifted_xsec,
			 HistType data)
      : fNominalXSec(nominal_xsec),
	fShiftedXSec(shifted_xsec),
	fData(data)
    {}

    CrossSectionAnalysis(std::string name, Systematic<CrossSectionType> shifted_xsec)
      : fShiftedXSec({name, shifted_xsec})
    {}

    CrossSectionAnalysis(HistType data)
      : fData(data)
    {}


    /// \brief Forward to UncertaintyPropagator
    HistType AbsoluteUncertaintyXSec(std::string syst_name, double ntargets);
    HistType AbsoluteUncertaintyUnfoldedXSec(std::string syst_name, double ntargets);

    /// \brief Forward to UncertaintyPropagator
    HistType FractionalUncertaintyXSec(std::string syst_name, double ntargets);
    HistType FractionalUncertaintyUnfoldedXSec(std::string syst_name, double ntargets);

    /// \brief Forward to UncertaintyPropagator
    std::pair<HistType, HistType> TotalAbsoluteUncertaintyXSec          (double ntargets);
    std::pair<HistType, HistType> TotalAbsoluteUncertaintyUnfoldedXSec  (double ntargets);
    std::pair<HistType, HistType> TotalFractionalUncertaintyXSec        (double ntargets);
    std::pair<HistType, HistType> TotalFractionalUncertaintyUnfoldedXSec(double ntargets);

    /// \brief Return an unfolded cross section result for the input systematic
    const Systematic<HistType> & UnfoldedCrossSection(std::string syst_name, double ntargets);

    /// \brief Return the nominal unfolded cross section result
    const HistType & UnfoldedCrossSection(double ntargets);

    /// \brief Return an folded cross section result for the input systematic
    const Systematic<HistType> & CrossSection(std::string syst_name, double ntargets);

    /// \brief Return the nominal folded cross section result
    const HistType & CrossSection(double ntargets);

    //    ~CrossSectionAnalysis();

    void SaveTo(TDirectory * dir, std::string subdir) const;
    static std::unique_ptr<CrossSectionAnalysis> LoadFrom(TDirectory * dir, std::string name);

    ~CrossSectionAnalysis();

  protected:
    // nominal is special
    CrossSectionType fNominalXSec;

    std::map<std::string, Systematic<CrossSectionType> > fShiftedXSec;

    const HistType fData;

    // cache the nominal results
    HistType * fUnfoldedNominalXSecResult = 0;
    HistType * fNominalXSecResult = 0;

    // cache the unfolded shifted results
    std::map<std::string, Systematic<HistType> > fUnfoldedShiftedXSecResult;
    std::map<std::string, Systematic<HistType> > fShiftedXSecResult;

    //
    UncertaintyPropagator fUncertaintyPropagator;
  };

  ///////////////////////////////////////////////////////////////////////
  template<class CrossSectionType,
	   class UncertaintyPropagator,
	   class HistType>
  HistType
  CrossSectionAnalysis<CrossSectionType,
		       UncertaintyPropagator,
		       HistType>::
  AbsoluteUncertaintyUnfoldedXSec(std::string syst_name, double ntargets)
  {
    return fUncertaintyPropagator.AbsoluteUncertaintyUnfolded(fData,
							      fNominalXSec,
							      fShiftedXSec.at(syst_name),
							      ntargets);
  }


  ///////////////////////////////////////////////////////////////////////
  template<class CrossSectionType,
	   class UncertaintyPropagator,
	   class HistType>
  HistType
  CrossSectionAnalysis<CrossSectionType,
		       UncertaintyPropagator,
		       HistType>::
  AbsoluteUncertaintyXSec(std::string syst_name, double ntargets)
  {
    return fUncertaintyPropagator.AbsoluteUncertaintyXSec(fData,
							  fNominalXSec,
							  fShiftedXSec.at(syst_name),
							  ntargets);
  }

  ///////////////////////////////////////////////////////////////////////
  template<class CrossSectionType,
	   class UncertaintyPropagator,
	   class HistType>
  HistType
  CrossSectionAnalysis<CrossSectionType,
		       UncertaintyPropagator,
		       HistType>::
  FractionalUncertaintyUnfoldedXSec(std::string syst_name, double ntargets)
  {
    return fUncertaintyPropagator.FractionalUncertaintyUnfoldedXSec(fData,
								    fNominalXSec,
								    fShiftedXSec.at(syst_name),
								    ntargets);
  }

  ///////////////////////////////////////////////////////////////////////
  template<class CrossSectionType,
	   class UncertaintyPropagator,
	   class HistType>
  HistType
  CrossSectionAnalysis<CrossSectionType,
		       UncertaintyPropagator,
		       HistType>::
  FractionalUncertaintyXSec(std::string syst_name, double ntargets)
  {
    return fUncertaintyPropagator.FractionalUncertaintyXSec(fData,
							    fNominalXSec,
							    fShiftedXSec.at(syst_name),
							    ntargets);
  }

  ///////////////////////////////////////////////////////////////////////
  template<class CrossSectionType,
	   class UncertaintyPropagator,
	   class HistType>
  std::pair<HistType,HistType>
  CrossSectionAnalysis<CrossSectionType,
		       UncertaintyPropagator,
		       HistType>::
  TotalAbsoluteUncertaintyUnfoldedXSec(double ntargets)
  {
    return fUncertaintyPropagator.TotalAbsoluteUncertaintyUnfoldedXSec(fData,
								       fNominalXSec,
								       fShiftedXSec,
								       ntargets);
  }

  ///////////////////////////////////////////////////////////////////////
  template<class CrossSectionType,
	   class UncertaintyPropagator,
	   class HistType>
  std::pair<HistType,HistType>
  CrossSectionAnalysis<CrossSectionType,
		       UncertaintyPropagator,
		       HistType>::
  TotalAbsoluteUncertaintyXSec(double ntargets)
  {
    return fUncertaintyPropagator.TotalAbsoluteUncertaintyXSec(fData,
							       fNominalXSec,
							       fShiftedXSec,
							       ntargets);
  }

  ///////////////////////////////////////////////////////////////////////
  template<class CrossSectionType,
	   class UncertaintyPropagator,
	   class HistType>
  std::pair<HistType,HistType>
  CrossSectionAnalysis<CrossSectionType,
		       UncertaintyPropagator,
		       HistType>::
  TotalFractionalUncertaintyUnfoldedXSec(double ntargets)
  {
    return fUncertaintyPropagator.TotalFractionalUncertaintyUnfoldedXSec(fData,
									 fNominalXSec,
									 fShiftedXSec,
									 ntargets);
  }

  ///////////////////////////////////////////////////////////////////////
  template<class CrossSectionType,
	   class UncertaintyPropagator,
	   class HistType>
  std::pair<HistType,HistType>
  CrossSectionAnalysis<CrossSectionType,
		       UncertaintyPropagator,
		       HistType>::
  TotalFractionalUncertaintyXSec(double ntargets)
  {
    return fUncertaintyPropagator.TotalFractionalUncertaintyXSec(fData,
								 fNominalXSec,
								 fShiftedXSec,
								 ntargets);
  }

  ///////////////////////////////////////////////////////////////////////
  template<class CrossSectionType,
	   class UncertaintyPropagator,
	   class HistType>
  const Systematic<HistType> &
  CrossSectionAnalysis<CrossSectionType,
		       UncertaintyPropagator,
		       HistType>::
  UnfoldedCrossSection(std::string syst_name,
		       double ntargets)
  {
    if(fUnfoldedShiftedXSecResult.find(syst_name) == fUnfoldedShiftedXSecResult.end()) {
      fUnfoldedShiftedXSecResult[syst_name] =
	fShiftedXSec.at(syst_name).Invoke(&CrossSectionType::template UnfoldedCrossSection<HistType>,
					  fData,
					  ntargets);

    }
    return fUnfoldedShiftedXSecResult.at(syst_name);
  }

  ///////////////////////////////////////////////////////////////////////
  template<class CrossSectionType,
	   class UncertaintyPropagator,
	   class HistType>
  const Systematic<HistType> &
  CrossSectionAnalysis<CrossSectionType,
		       UncertaintyPropagator,
		       HistType>::
  CrossSection(std::string syst_name,
	       double ntargets)
  {
    if(fShiftedXSecResult.find(syst_name) == fShiftedXSecResult.end()) {
      fShiftedXSecResult[syst_name] =
	fShiftedXSec.at(syst_name).Invoke(&CrossSectionType::template CrossSection<HistType>,
					  fData,
					  ntargets);

    }
    return fShiftedXSecResult.at(syst_name);
  }

  ///////////////////////////////////////////////////////////////////////
  template<class CrossSectionType,
	   class UncertaintyPropagator,
	   class HistType>
  const HistType &
  CrossSectionAnalysis<CrossSectionType,
		       UncertaintyPropagator,
		       HistType>::
  CrossSection(double ntargets)
  {
    if(!fNominalXSecResult)  {
      fNominalXSecResult = new HistType(fNominalXSec.CrossSection(fData, ntargets));
    }
    return *fNominalXSecResult;
  }
  ///////////////////////////////////////////////////////////////////////
  template<class CrossSectionType,
	   class UncertaintyPropagator,
	   class HistType>
  const HistType &
  CrossSectionAnalysis<CrossSectionType,
		       UncertaintyPropagator,
		       HistType>::
  UnfoldedCrossSection(double ntargets)
  {
    if(!fUnfoldedNominalXSecResult)  {
      fUnfoldedNominalXSecResult = new HistType(fNominalXSec.UnfoldedCrossSection(fData, ntargets));
    }
    return *fUnfoldedNominalXSecResult;
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

    fNominalXSec.SaveTo(dir, "fNominalXSec");
    fData.SaveTo(dir, "fData");

    auto syst_dir = dir->mkdir("fShiftedXSec");
    for(auto shifted : fShiftedXSec) {
      shifted.second.SaveTo(syst_dir, shifted.first);
    }

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

    auto data = *HistType::LoadFrom(dir, "fData");

    auto nominal_xsec = *CrossSectionType::LoadFrom(dir, "fNominalXSec");


    std::map<std::string, Systematic<CrossSectionType> > shifted_xsec;
    auto syst_dir = dir->GetDirectory("fShiftedXSec");
    for(auto syst_name : *syst_dir->GetListOfKeys()) {
      shifted_xsec[syst_name->GetName()] = *Systematic<CrossSectionType>::LoadFrom(syst_dir,
										   syst_name->GetName()).release();
    }

    tmp->cd();
    return std::make_unique<CrossSectionAnalysis<CrossSectionType,
						 UncertaintyPropagator,
						 HistType> >
      (nominal_xsec,
       shifted_xsec,
       data);

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
    if(fUnfoldedNominalXSecResult) delete fUnfoldedNominalXSecResult;
    if(fNominalXSecResult) delete fNominalXSecResult;
  }

}
