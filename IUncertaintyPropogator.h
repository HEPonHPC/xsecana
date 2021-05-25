#pragma once

#include "XSecAna/Systematic.h"
#include "XSecAna/Hist.h"
#include "XSecAna/IUnfold.h"
namespace xsec {
  /// \brief an UncertaintyPropogator is module
  /// that calculates uncertainty for an analysis
  ///
  /// IUncertaintyPropogator defines the interface assumed
  /// by the analysis object
  template<class CrossSectionType,
	   class HistType = HistXd>
  class IUncertaintyPropogator {
  public:
    virtual std::pair<HistType*, HistType*> 
    TotalFractionalUncertaintyUnfoldedXSec(const HistType & data,
					   CrossSectionType & nominal_xsec,
					   std::map<std::string, Systematic<CrossSectionType> > & shifted_xsec,
					   const IUnfold<HistType> * unfold,
					   double ntargets) = 0;

    virtual std::pair<HistType*, HistType*> 
    TotalFractionalUncertaintyXSec(const HistType & data,
				   CrossSectionType & nominal_xsec,
				   std::map<std::string, Systematic<CrossSectionType> > & shifted_xsec,
				   double ntargets) = 0;

    virtual std::pair<HistType*, HistType*> 
    TotalAbsoluteUncertaintyUnfoldedXSec(const HistType & data,
					 CrossSectionType & nominal_xsec,
					 std::map<std::string, Systematic<CrossSectionType> > & shifted_xsec,
					 const IUnfold<HistType> * unfold,
					 double ntargets) = 0;

    virtual std::pair<HistType*, HistType*> 
    TotalAbsoluteUncertaintyXSec(const HistType & data,
				 CrossSectionType & nominal_xsec,
				 std::map<std::string, Systematic<CrossSectionType> > & shifted_xsec,
				 double ntargets) = 0;

    virtual HistType *
    FractionalUncertaintyUnfoldedXSec(const HistType & data,
				      CrossSectionType & nominal_xsec,
				      Systematic<CrossSectionType> & shifted_xsec,
				      const IUnfold<HistType> * unfold,
				      double ntargets) = 0;

    virtual HistType *
    FractionalUncertaintyXSec(const HistType & data,
			      CrossSectionType & nominal_xsec,
			      Systematic<CrossSectionType>  & shifted_xsec,
			      double ntargets) = 0;

    virtual HistType *
    AbsoluteUncertaintyUnfoldedXSec(const HistType & data,
				    CrossSectionType & nominal_xsec,
				    Systematic<CrossSectionType> & shifted_xsec,
				    const IUnfold<HistType> * unfold,
				    double ntargets) = 0;
    virtual HistType *
    AbsoluteUncertaintyXSec(const HistType & data,
			    CrossSectionType & nominal_xsec,
			    Systematic<CrossSectionType> & shifted_xsec,
			    double ntargets) = 0;

  };
  
  // inline some common functions
  template<typename Scalar,
	   int Cols>
  inline Hist<Scalar, Cols>
  MaxShift(const Hist<Scalar, Cols> & h1,
	   const Hist<Scalar, Cols> & h2)
  {
    Hist<Scalar, Cols> ret = h1; // copy binning

    // stack
    Eigen::Matrix<Scalar, 2, Cols> stack;
    stack << h1->Contents(),h2->Contents();
    ret.fContents = stack.colwise().maxCoeff();
    return ret;
  }

  template<typename Scalar,
	   int Cols>
  inline Hist<Scalar, Cols>
  QuadSum(std::vector<Hist<Scalar, Cols> > deltas)
  {
    // put deltas into an Eigen::Matrix for efficiency column-wise operations
    Eigen::Matrix<Scalar, Eigen::Dynamic, Cols> mat(deltas.size(), Cols);
    for(auto irow = 0u; irow < deltas.size(); irow++) {
      mat.row(irow) = deltas[irow].Contents();
    }
    return mat.colwise().abs2();
  }
}
