#pragma once

#include "XSecAna/Systematic.h"
#include "XSecAna/Hist.h"

namespace xsec {
  /// \brief an UncertaintyPropogator is module
  /// that calculates uncertainty for an analysis
  ///
  /// IUncertaintyPropogator defines the interface assumed
  /// by the analysis object
  template<class CrossSectionType,
	   class UnfoldType,
	   class HistType = HistXXd>
  class IUncertaintyPropogator {
    virtual std::pair<HistType*, HistType*> 
    TotalFractionalUncertaintyUnfoldedXSec(const UnfoldType * unfold,
					   CrossSectionType & nominal_xsec,
					   std::map<std::string, Systematic<CrossSectionType> > shifted_xsec,
					   double ntargets) = 0;
    virtual std::pair<HistType*, HistType*> 
    TotalFractionalUncertaintyXSec(CrossSectionType & nominal_xsec,
				   std::map<std::string, Systematic<CrossSectionType> > shifted_xsec,
				   double ntargets) = 0;

    virtual std::pair<HistType*, HistType*> 
    TotalAbsoluteUncertaintyUnfoldedXSec(const UnfoldType * unfold,
					 CrossSectionType & nominal_xsec,
					 std::map<std::string, Systematic<CrossSectionType> > shifted_xsec,
					 double ntargets) = 0;
    virtual std::pair<HistType*, HistType*> 
    TotalAbsoluteUncertaintyXSec(CrossSectionType & nominal_xsec,
				 std::map<std::string, Systematic<CrossSectionType> > shifted_xsec,
				 double ntargets) = 0;

    virtual HistType *
    FractionalUncertaintyUnfoldedXSec(const UnfoldType * unfold,
				      CrossSectionType & nominal_xsec,
				      Systematic<CrossSectionType> shifted_xsec,
				      double ntargets) = 0;
    virtual HistType *
    FractionalUncertaintyXSec(CrossSectionType & nominal_xsec,
			      Systematic<CrossSectionType>  shifted_xsec,
			      double ntargets) = 0;

    virtual HistType *
    AbsoluteUncertaintyUnfoldedXSec(const UnfoldType * unfold,
				    CrossSectionType & nominal_xsec,
				    Systematic<CrossSectionType> shifted_xsec,
				    double ntargets) = 0;
    virtual HistType *
    AbsoluteUncertaintyXSec(CrossSectionType & nominal_xsec,
			    Systematic<CrossSectionType> shifted_xsec,
			    double ntargets) = 0;

  };
}
