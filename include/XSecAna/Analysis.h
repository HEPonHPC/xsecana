#pragma once

#include <string>
#include <map>

#include "XSecAna/Systematic.h"
#include "XSecAna/Hist.h"

// root includes
#include "TDirectory.h"

namespace xsec {
  template<class MeasurementType,
	   class UncertaintyPropagator,
	   class HistType = HistXd>
  class Analysis {

  public:

    Analysis() {}
    Analysis(MeasurementType nominal_measurement,
	     std::map<std::string, Systematic<MeasurementType> > shifted_measurement,
	     HistType data)
      : fNominalMeasure(nominal_measurement),
	fShiftedMeasure(shifted_measurement),
	fData(data)
    {}

    Analysis(std::string name, Systematic<MeasurementType> shifted_measurement)
      : fShiftedMeasure({name, shifted_measurement})
    {}

    Analysis(HistType data)
      : fData(data)
    {}


    /// \brief Forward to UncertaintyPropagator
    HistType AbsoluteUncertainty(std::string syst_name);

    /// \brief Forward to UncertaintyPropagator
    HistType FractionalUncertainty(std::string syst_name);

    /// \brief Forward to UncertaintyPropagator
    std::pair<HistType, HistType> TotalAbsoluteUncertainty  ();
    std::pair<HistType, HistType> TotalFractionalUncertainty();

    /// \brief Return an folded cross section result for the input systematic
    const Systematic<HistType> & Result(std::string syst_name);

    /// \brief Return the nominal folded cross section result
    const HistType & Result();

    //    ~Analysis();

    void SaveTo(TDirectory * dir, std::string subdir) const;
    static std::unique_ptr<Analysis> LoadFrom(TDirectory * dir, std::string name);

    ~Analysis();

  protected:
    // nominal is special
    MeasurementType fNominalMeasure;

    std::map<std::string, Systematic<MeasurementType> > fShiftedMeasure;

    const HistType fData;

    // cache the nominal results
    HistType * fNominalResult = 0;

    // cache the shifted results
    std::map<std::string, Systematic<HistType> > fShiftedResult;

    //
    UncertaintyPropagator fUncertaintyPropagator;
  };

  ///////////////////////////////////////////////////////////////////////
  template<class MeasurementType,
	   class UncertaintyPropagator,
	   class HistType>
  HistType
  Analysis<MeasurementType,
	   UncertaintyPropagator,
	   HistType>::
  AbsoluteUncertainty(std::string syst_name)
  {
    return fUncertaintyPropagator.AbsoluteUncertainty(fData,
						      fNominalMeasure,
						      fShiftedMeasure.at(syst_name));
  }

  ///////////////////////////////////////////////////////////////////////
  template<class MeasurementType,
	   class UncertaintyPropagator,
	   class HistType>
  HistType
  Analysis<MeasurementType,
	   UncertaintyPropagator,
	   HistType>::
  FractionalUncertainty(std::string syst_name)
  {
    return fUncertaintyPropagator.FractionalUncertainty(fData,
							fNominalMeasure,
							fShiftedMeasure.at(syst_name));
  }

  ///////////////////////////////////////////////////////////////////////
  template<class MeasurementType,
	   class UncertaintyPropagator,
	   class HistType>
  std::pair<HistType,HistType>
  Analysis<MeasurementType,
	   UncertaintyPropagator,
	   HistType>::
  TotalAbsoluteUncertainty()
  {
    return fUncertaintyPropagator.TotalAbsoluteUncertainty(fData,
							   fNominalMeasure,
							   fShiftedMeasure);
  }

  ///////////////////////////////////////////////////////////////////////
  template<class MeasurementType,
	   class UncertaintyPropagator,
	   class HistType>
  std::pair<HistType,HistType>
  Analysis<MeasurementType,
	   UncertaintyPropagator,
	   HistType>::
  TotalFractionalUncertainty()
  {
    return fUncertaintyPropagator.TotalFractionalUncertainty(fData,
							     fNominalMeasure,
							     fShiftedMeasure);
  }

  ///////////////////////////////////////////////////////////////////////
  template<class MeasurementType,
	   class UncertaintyPropagator,
	   class HistType>
  const Systematic<HistType> &
  Analysis<MeasurementType,
	   UncertaintyPropagator,
	   HistType>::
  Result(std::string syst_name)
  {
    if(fShiftedResult.find(syst_name) == fShiftedResult.end()) {
      fShiftedResult[syst_name] =
	fShiftedMeasure.at(syst_name).Invoke(&MeasurementType::Result,
					     fData);

    }
    return fShiftedResult.at(syst_name);
  }

  ///////////////////////////////////////////////////////////////////////
  template<class MeasurementType,
	   class UncertaintyPropagator,
	   class HistType>
  const HistType &
  Analysis<MeasurementType,
	   UncertaintyPropagator,
	   HistType>::
  Result()
  {
    if(!fNominalResult)  {
      fNominalResult = new HistType(fNominalMeasure.Result(fData));
    }
    return *fNominalResult;
  }

  ///////////////////////////////////////////////////////////////////////
  template<class MeasurementType,
	   class UncertaintyPropagator,
	   class HistType>
  void
  Analysis<MeasurementType,
	   UncertaintyPropagator,
	   HistType>::
  SaveTo(TDirectory * dir, std::string subdir) const
  {
    TDirectory * tmp = gDirectory;

    dir = dir->mkdir(subdir.c_str()); // switch to subdir
    dir->cd();
    TObjString("Analysis").Write("type");

    fNominalMeasure.SaveTo(dir, "fNominalMeasure");
    fData.SaveTo(dir, "fData");

    auto syst_dir = dir->mkdir("fShiftedMeasure");
    for(auto shifted : fShiftedMeasure) {
      shifted.second.SaveTo(syst_dir, shifted.first);
    }

    tmp->cd();
  }

  ///////////////////////////////////////////////////////////////////////
  template<class MeasurementType,
	   class UncertaintyPropagator,
	   class HistType>
  std::unique_ptr<Analysis<MeasurementType,
			   UncertaintyPropagator,
			   HistType> >
  Analysis<MeasurementType,
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

    auto nominal_measurement = *MeasurementType::LoadFrom(dir, "fNominalMeasure");


    std::map<std::string, Systematic<MeasurementType> > shifted_measurement;
    auto syst_dir = dir->GetDirectory("fShiftedMeasure");
    for(auto syst_name : *syst_dir->GetListOfKeys()) {
      shifted_measurement[syst_name->GetName()] = *Systematic<MeasurementType>::LoadFrom(syst_dir,
											 syst_name->GetName()).release();
    }

    tmp->cd();
    return std::make_unique<Analysis<MeasurementType,
				     UncertaintyPropagator,
				     HistType> >
      (nominal_measurement,
       shifted_measurement,
       data);

  }

  ///////////////////////////////////////////////////////////////////////
  template<class MeasurementType,
	   class UncertaintyPropagator,
	   class HistType>
  Analysis<MeasurementType,
	   UncertaintyPropagator,
	   HistType>::
  ~Analysis()
  {
    if(fNominalResult) delete fNominalResult;
  }

}
