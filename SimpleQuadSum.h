#pragma once

#include "XSecAna/Systematic.h"
#include "XSecAna/Hist.h"
#include "XSecAna/IUnfold.h"
#include "XSecAna/IUncertaintyPropogator.h"

#include <Eigen/Dense>

namespace xsec {
  /// \brief an UncertaintyPropogator is module
  /// that calculates uncertainty for an analysis
  ///
  /// SimpleQuadSum performs the quadrature sum of systematic shifts
  /// In the case of an asymmetric 2-sided shift,
  /// the shift is symmeterized by taking the largest shift
  template<class CrossSectionType,
	   class HistType = HistXd>
  class SimpleQuadSum : IUncertaintyPropogator<CrossSectionType,
					       HistType>
  {
  public:
    std::pair<HistType*, HistType*> 
    TotalFractionalUncertaintyUnfoldedXSec(const HistType & data,
					   CrossSectionType & nominal_xsec,
					   std::map<std::string, Systematic<CrossSectionType> > & shifted_xsec,
					   const IUnfold<HistType> * unfold,
					   double ntargets) override;

    std::pair<HistType*, HistType*> 
    TotalFractionalUncertaintyXSec(const HistType & data,
				   CrossSectionType & nominal_xsec,
				   std::map<std::string, Systematic<CrossSectionType> > & shifted_xsec,
				   double ntargets) override;

    std::pair<HistType*, HistType*> 
    TotalAbsoluteUncertaintyUnfoldedXSec(const HistType & data,
					 CrossSectionType & nominal_xsec,
					 std::map<std::string, Systematic<CrossSectionType> > & shifted_xsec,
					 const IUnfold<HistType> * unfold,
					 double ntargets) override;

    std::pair<HistType*, HistType*> 
    TotalAbsoluteUncertaintyXSec(const HistType & data,
				 CrossSectionType & nominal_xsec,
				 std::map<std::string, Systematic<CrossSectionType> > & shifted_xsec,
				 double ntargets) override;

    HistType *
    FractionalUncertaintyUnfoldedXSec(const HistType & data,
				      CrossSectionType & nominal_xsec,
				      Systematic<CrossSectionType> & shifted_xsec,
				      const IUnfold<HistType> * unfold,
				      double ntargets) override;

    HistType *
    FractionalUncertaintyXSec(const HistType & data,
			      CrossSectionType & nominal_xsec,
			      Systematic<CrossSectionType>  & shifted_xsec,
			      double ntargets) override;

    HistType *
    AbsoluteUncertaintyUnfoldedXSec(const HistType & data,
				    CrossSectionType & nominal_xsec,
				    Systematic<CrossSectionType> & shifted_xsec,
				    const IUnfold<HistType> * unfold,
				    double ntargets) override;
    HistType *
    AbsoluteUncertaintyXSec(const HistType & data,
			    CrossSectionType & nominal_xsec,
			    Systematic<CrossSectionType> & shifted_xsec,
			    double ntargets) override;

  };


  /////////////////////////////////////////////////////////////////////////
  template<class CrossSectionType,
	   class HistType>
  HistType *
  SimpleQuadSum<CrossSectionType, HistType>::
  AbsoluteUncertaintyXSec(const HistType & data,
			  CrossSectionType & nominal_xsec,
			  Systematic<CrossSectionType> & shifted_xsec,
			  double ntargets)
  {
    Systematic<HistType> * shifts = shifted_xsec.Invoke(&CrossSectionType::CrossSection, data, ntargets);
    shifts = shifts->Invoke(&HistType::operator-, nominal_xsec->CrossSection(data, ntargets));
    shifts = shifts->Invoke(&HistType::abs);
    
    return MaxShift(shifts->Up()->Data(),
		    shifts->Down()->Data());
  }


  /////////////////////////////////////////////////////////////////////////
  template<class CrossSectionType,
	   class HistType>
  HistType *
  SimpleQuadSum<CrossSectionType, HistType>::
  FractionalUncertaintyXSec(const HistType & data,
			    CrossSectionType & nominal_xsec,
			    Systematic<CrossSectionType> & shifted_xsec,
			    double ntargets)
  {
    HistType * abs = AbsoluteUncertaintyXSec(data, nominal_xsec, shifted_xsec, ntargets);
    return *abs / *nominal_xsec->CrossSection(data, ntargets);
  }


  /////////////////////////////////////////////////////////////////////////
  template<class CrossSectionType,
	   class HistType>
  HistType *
  SimpleQuadSum<CrossSectionType, HistType>::
  AbsoluteUncertaintyUnfoldedXSec(const HistType & data,
				  CrossSectionType & nominal_xsec,
				  Systematic<CrossSectionType> & shifted_xsec,
				  const IUnfold<HistType> * unfold,
				  double ntargets)
  {
    Systematic<HistType> * shifts = shifted_xsec.Invoke(&CrossSectionType::UnfoldedCrossSection, 
							data, 
							*unfold, 
							ntargets);

    shifts = shifts->Invoke(&HistType::operator-, nominal_xsec->UnfoldedCrossSection(data, *unfold, ntargets));
    shifts = shifts->Invoke(&HistType::abs);
    
    return MaxShift(shifts->Up()->Data(),
		    shifts->Down()->Data());
  }


  /////////////////////////////////////////////////////////////////////////
  template<class CrossSectionType,
	   class HistType>
  HistType *
  SimpleQuadSum<CrossSectionType, HistType>::
  FractionalUncertaintyUnfoldedXSec(const HistType & data,
				    CrossSectionType & nominal_xsec,
				    Systematic<CrossSectionType> & shifted_xsec,
				    const IUnfold<HistType> * unfold,
				    double ntargets)
  {
    HistType * abs = AbsoluteUncertaintyUnfoldedXSec(data, nominal_xsec, shifted_xsec, unfold, ntargets);
    return *abs / *nominal_xsec->UnfoldedCrossSection(data, *unfold, ntargets);
  }

  /////////////////////////////////////////////////////////////////////////
  template<class CrossSectionType,
	   class HistType>
  std::pair<HistType*, HistType*> 
  SimpleQuadSum<CrossSectionType, HistType>::
  TotalFractionalUncertaintyXSec(const HistType & data,
				 CrossSectionType & nominal_xsec,
				 std::map<std::string, Systematic<CrossSectionType> > & shifted_xsec,
				 double ntargets)
  {
    std::vector<HistType> shifts;
    for(auto syst_it = shifted_xsec.begin(); syst_it != shifted_xsec.end(); syst_it++) {
      shifts.push_back(FractionalUncertaintyXSec(data, nominal_xsec, syst_it->second, ntargets));
    }
    return QuadSum(shifts).sqrt();
  }

}
