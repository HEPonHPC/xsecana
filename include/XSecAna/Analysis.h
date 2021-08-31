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
  class Analysis {

  public:

    Analysis() {}
    Analysis(CrossSectionType nominal_xsec,
	     std::map<std::string, Systematic<CrossSectionType> > shifted_xsec,
	     HistType data)
      : fNominalXSec(nominal_xsec),
	fShiftedXSec(shifted_xsec),
	fData(data)
    {}

    Analysis(std::string name, Systematic<CrossSectionType> shifted_xsec)
      : fShiftedXSec({name, shifted_xsec})
    {}

    Analysis(HistType data)
      : fData(data)
    {}


    /// \brief Forward to UncertaintyPropagator
    HistType AbsoluteUncertainty(std::string syst_name, double ntargets);

    /// \brief Forward to UncertaintyPropagator
    HistType FractionalUncertainty(std::string syst_name, double ntargets);

    /// \brief Forward to UncertaintyPropagator
    std::pair<HistType, HistType> TotalAbsoluteUncertainty          (double ntargets);
    std::pair<HistType, HistType> TotalFractionalUncertainty        (double ntargets);

    /// \brief Return an folded cross section result for the input systematic
    const Systematic<HistType> & CrossSection(std::string syst_name, double ntargets);

    /// \brief Return the nominal folded cross section result
    const HistType & CrossSection(double ntargets);

    //    ~Analysis();

    void SaveTo(TDirectory * dir, std::string subdir) const;
    static std::unique_ptr<Analysis> LoadFrom(TDirectory * dir, std::string name);

    ~Analysis();

  protected:
    // nominal is special
    CrossSectionType fNominalXSec;

    std::map<std::string, Systematic<CrossSectionType> > fShiftedXSec;

    const HistType fData;

    // cache the nominal results
    HistType * fNominalXSecResult = 0;

    // cache the shifted results
    std::map<std::string, Systematic<HistType> > fShiftedXSecResult;

    //
    UncertaintyPropagator fUncertaintyPropagator;
  };

  ///////////////////////////////////////////////////////////////////////
  template<class CrossSectionType,
	   class UncertaintyPropagator,
	   class HistType>
  HistType
  Analysis<CrossSectionType,
	   UncertaintyPropagator,
	   HistType>::
  AbsoluteUncertainty(std::string syst_name, double ntargets)
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
  Analysis<CrossSectionType,
	   UncertaintyPropagator,
	   HistType>::
  FractionalUncertainty(std::string syst_name, double ntargets)
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
  Analysis<CrossSectionType,
	   UncertaintyPropagator,
	   HistType>::
  TotalAbsoluteUncertainty(double ntargets)
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
  Analysis<CrossSectionType,
	   UncertaintyPropagator,
	   HistType>::
  TotalFractionalUncertainty(double ntargets)
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
  Analysis<CrossSectionType,
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
  Analysis<CrossSectionType,
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
  void
  Analysis<CrossSectionType,
	   UncertaintyPropagator,
	   HistType>::
  SaveTo(TDirectory * dir, std::string subdir) const
  {
    TDirectory * tmp = gDirectory;

    dir = dir->mkdir(subdir.c_str()); // switch to subdir
    dir->cd();
    TObjString("Analysis").Write("type");

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
  std::unique_ptr<Analysis<CrossSectionType,
			   UncertaintyPropagator,
			   HistType> >
  Analysis<CrossSectionType,
	   UncertaintyPropagator,
	   HistType>::
  LoadFrom(TDirectory * dir, std::string subdir)
  {
    TDirectory * tmp = gDirectory;
    dir = dir->GetDirectory(subdir.c_str());

    TObjString * ptag = (TObjString*) dir->Get("type");
    assert(ptag->GetString() == "Analysis" && "Type does not match Analysis");
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
    return std::make_unique<Analysis<CrossSectionType,
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
  Analysis<CrossSectionType,
	   UncertaintyPropagator,
	   HistType>::
  ~Analysis()
  {
    if(fNominalXSecResult) delete fNominalXSecResult;
  }

}
