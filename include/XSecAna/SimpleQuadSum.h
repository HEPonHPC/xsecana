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
  template<class CrossSectionType,
	   class HistType = HistXd>
  class SimpleQuadSum : IUncertaintyPropagator<CrossSectionType,
					       HistType>
  {
  public:
    std::pair<HistType, HistType> 
    TotalFractionalUncertainty(const HistType & data,
			       CrossSectionType & nominal_xsec,
			       std::map<std::string, Systematic<CrossSectionType> > & shifted_xsec,
			       double ntargets) override;

    std::pair<HistType, HistType> 
    TotalAbsoluteUncertainty(const HistType & data,
			     CrossSectionType & nominal_xsec,
			     std::map<std::string, Systematic<CrossSectionType> > & shifted_xsec,
			     double ntargets) override;

    HistType
    FractionalUncertainty(const HistType & data,
			  CrossSectionType & nominal_xsec,
			  Systematic<CrossSectionType> & shifted_xsec,
			  double ntargets) override;

    HistType
    AbsoluteUncertainty(const HistType & data,
			CrossSectionType & nominal_xsec,
			Systematic<CrossSectionType> & shifted_xsec,
			double ntargets) override;

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
  template<class CrossSectionType,
	   class HistType>
  HistType
  SimpleQuadSum<CrossSectionType, HistType>::
  AbsoluteUncertainty(const HistType & data,
		      CrossSectionType & nominal_xsec,
		      Systematic<CrossSectionType> & shifted_xsec,
		      double ntargets)
  {
    // calculate cross sections
    auto hnominal_xsec = nominal_xsec.CrossSection(data, ntargets);   
    Systematic<HistType> shifts = shifted_xsec.Invoke(&CrossSectionType::template CrossSection<HistType>, data, ntargets);

    // convert multiverse systematic to two sided by finding 1sigma
    shifts = HandleMultiverseSystematic(shifts, hnominal_xsec);
    
    // HistType::operator- is overloaded
    // static cast to resolve
    shifts = shifts.Invoke(static_cast<HistType(HistType::*)(const HistType&) const>(&HistType::operator-),
			   hnominal_xsec);
    
    return MaxShift(shifts.Up().abs(), shifts.Down().abs());
  }

  /////////////////////////////////////////////////////////////////////////
  template<class CrossSectionType,
	   class HistType>
  HistType
  SimpleQuadSum<CrossSectionType, HistType>::
  FractionalUncertainty(const HistType & data,
			CrossSectionType & nominal_xsec,
			Systematic<CrossSectionType> & shifted_xsec,
			double ntargets)
  {
    HistType abs = AbsoluteUncertainty(data, nominal_xsec, shifted_xsec, ntargets);
    return abs / nominal_xsec.CrossSection(data, ntargets);
  }

  /////////////////////////////////////////////////////////////////////////
  template<class CrossSectionType,
	   class HistType>
  std::pair<HistType, HistType> 
  SimpleQuadSum<CrossSectionType, HistType>::
  TotalFractionalUncertainty(const HistType & data,
			     CrossSectionType & nominal_xsec,
			     std::map<std::string, Systematic<CrossSectionType> > & shifted_xsec,
			     double ntargets)
  {
    auto hnominal = nominal_xsec.CrossSection(data, ntargets);
    auto abs = TotalAbsoluteUncertainty(data, nominal_xsec, shifted_xsec, ntargets);
    return {abs.first / hnominal, abs.second / hnominal};
  }

  /////////////////////////////////////////////////////////////////////////
  template<class CrossSectionType,
	   class HistType>
  std::pair<HistType, HistType> 
  SimpleQuadSum<CrossSectionType, HistType>::
  TotalAbsoluteUncertainty(const HistType & data,
			   CrossSectionType & nominal_xsec,
			   std::map<std::string, Systematic<CrossSectionType> > & shifted_xsec,
			   double ntargets)
  {
    std::vector<HistType> shifts;
    for(auto syst_it = shifted_xsec.begin(); syst_it != shifted_xsec.end(); syst_it++) {
      shifts.push_back(AbsoluteUncertainty(data, nominal_xsec, syst_it->second, ntargets));
    }
    auto result = QuadSum(shifts).sqrt();
    return {result, result};
  }
}
