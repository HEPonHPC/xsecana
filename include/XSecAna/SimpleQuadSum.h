#pragma once

#include "XSecAna/Systematic.h"
#include "XSecAna/Hist.h"
#include "XSecAna/IUnfold.h"
#include "XSecAna/IUncertaintyPropagator.h"

#include <Eigen/Dense>

namespace xsec {
  
  /// \brief an UncertaintyPropagator is module
  /// that calculates uncertainty for an analysis
  ///
  /// SimpleQuadSum performs the quadrature sum of systematic shifts
  /// In the case of an asymmetric 2-sided shift,
  /// the shift is symmeterized by taking the largest shift
  template<class MeasurementType,
	   class HistType = HistXd>
  class SimpleQuadSum : IUncertaintyPropagator<MeasurementType,
					       HistType>
  {
  public:
    std::pair<HistType, HistType> 
    TotalFractionalUncertainty(const HistType & data,
			       MeasurementType & nominal_measurement,
			       std::map<std::string, Systematic<MeasurementType> > & shifted_measurement) override;

    std::pair<HistType, HistType> 
    TotalAbsoluteUncertainty(const HistType & data,
			     MeasurementType & nominal_measurement,
			     std::map<std::string, Systematic<MeasurementType> > & shifted_measurement) override;

    HistType
    FractionalUncertainty(const HistType & data,
			  MeasurementType & nominal_measurement,
			  Systematic<MeasurementType> & shifted_measurement) override;

    HistType
    AbsoluteUncertainty(const HistType & data,
			MeasurementType & nominal_measurement,
			Systematic<MeasurementType> & shifted_measurement) override;

  };

  /////////////////////////////////////////////////////////////////////////
  template<class HistType>
  inline Systematic<HistType> HandleMultiverseSystematic(const Systematic<HistType> & syst, 
							 const HistType & nominal)
  {
    // if given a multiverse systematic, return a new Systematic<HistType> that contains
    // +/- 1 sigma 
    if(syst.GetType() == SystType_t::kMultiverse) {
      return Systematic<HistType>(syst.GetName(),
				  syst.NSigmaShift( 1, nominal),
				  syst.NSigmaShift(-1, nominal));
    }
    else {
      return syst;
    }

  }

  /////////////////////////////////////////////////////////////////////////
  template<class MeasurementType,
	   class HistType>
  HistType
  SimpleQuadSum<MeasurementType, HistType>::
  AbsoluteUncertainty(const HistType & data,
		      MeasurementType & nominal_measurement,
		      Systematic<MeasurementType> & shifted_measurement)
  {
    // calculate cross sections
    auto hnominal_measurement = nominal_measurement.Result(data);
    Systematic<HistType> shifts = shifted_measurement.Invoke(&MeasurementType::Result, data);

    // convert multiverse systematic to two sided by finding 1sigma
    shifts = HandleMultiverseSystematic(shifts, hnominal_measurement);
    
    // HistType::operator- is overloaded
    // static cast to resolve
    shifts = shifts.Invoke(static_cast<HistType(HistType::*)(const HistType&) const>(&HistType::operator-),
			   hnominal_measurement);
    
    return MaxShift(shifts.Up().abs(), shifts.Down().abs());
  }

  /////////////////////////////////////////////////////////////////////////
  template<class MeasurementType,
	   class HistType>
  HistType
  SimpleQuadSum<MeasurementType, HistType>::
  FractionalUncertainty(const HistType & data,
			MeasurementType & nominal_measurement,
			Systematic<MeasurementType> & shifted_measurement)
  {
    HistType abs = AbsoluteUncertainty(data, nominal_measurement, shifted_measurement);
    return abs / nominal_measurement.Result(data);
  }

  /////////////////////////////////////////////////////////////////////////
  template<class MeasurementType,
	   class HistType>
  std::pair<HistType, HistType> 
  SimpleQuadSum<MeasurementType, HistType>::
  TotalFractionalUncertainty(const HistType & data,
			     MeasurementType & nominal_measurement,
			     std::map<std::string, Systematic<MeasurementType> > & shifted_measurement)
  {
    auto hnominal = nominal_measurement.Result(data);
    auto abs = TotalAbsoluteUncertainty(data, nominal_measurement, shifted_measurement);
    return {abs.first / hnominal, abs.second / hnominal};
  }

  /////////////////////////////////////////////////////////////////////////
  template<class MeasurementType,
	   class HistType>
  std::pair<HistType, HistType> 
  SimpleQuadSum<MeasurementType, HistType>::
  TotalAbsoluteUncertainty(const HistType & data,
			   MeasurementType & nominal_measurement,
			   std::map<std::string, Systematic<MeasurementType> > & shifted_measurement)
  {
    std::vector<HistType> shifts;
    for(auto syst_it = shifted_measurement.begin(); syst_it != shifted_measurement.end(); syst_it++) {
      shifts.push_back(AbsoluteUncertainty(data, nominal_measurement, syst_it->second));
    }
    auto result = QuadSum(shifts).sqrt();
    return {result, result};
  }
}
