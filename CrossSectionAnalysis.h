#pragma once

#include <string>
#include <map>

#include "XSecAna/Systematic.h"
#include "XSecAna/Hist.h"

namespace xsec {
  template<class CrossSectionType,
	   class UnfoldType,
	   class UncertaintyPropogator,
	   class HistType>
  class CrossSectionAnalysis {

  public:

    /// \brief Forward to UncertaintyPropogator
    HistType * AbsoluteUncertaintyXSec(std::string syst_name, double ntargets);
    HistType * AbsoluteUncertaintyUnfoldedXSec(std::string syst_name, double ntargets);

    /// \brief Forward to UncertaintyPropogator
    HistType * FractionalUncertaintyXSec(std::string syst_name, double ntargets);
    HistType * FractionalUncertaintyUnfoldedXSec(std::string syst_name, double ntargets);

    /// \brief Forward to UncertaintyPropogator
    std::pair<HistType*, HistType*> TotalAbsoluteUncertaintyXSec          (double ntargets);
    std::pair<HistType*, HistType*> TotalAbsoluteUncertaintyUnfoldedXSec  (double ntargets);
    std::pair<HistType*, HistType*> TotalFractionalUncertaintyXSec        (double ntargets);
    std::pair<HistType*, HistType*> TotalFractionalUncertaintyUnfoldedXSec(double ntargets);

    /// \brief Return an unfolded cross section result for the input systematic 
    /// pass syst_name = "nominal" for the nominal result
    const HistType * UnfoldedCrossSection(std::string syst_name, double ntargets);

    ~CrossSectionAnalysis();

    void SaveTo(TDirectory * dir, std::string subdir) const;
    static std::unique_ptr<CrossSectionAnalysis> LoadFrom(TDirectory * dir, std::string name);

  protected:
    ///\brief constructor for loading analysis from file
    CrossSectionAnalysis(CrossSectionType * nominal_xsec,
			 std::map<std::string, Systematic<CrossSectionType> > shifted_xsec,
			 UnfoldType * unfold,
			 HistType * data)
      : fNominalXSec(nominal_xsec),
	fShiftedXSec(shifted_xsec),
	fUnfold(unfold),
	fData(data)
    {}

    // nominal is special 
    CrossSectionType * fNominalXSec;

    std::map<std::string, Systematic<CrossSectionType> > fShiftedXSec;

    UnfoldType * fUnfold = NULL;
    const HistType * fData  = NULL;

    // cache the nominal unfolded results
    HistType * fUnfoldedNominalXSec;

    // cache the unfolded shifted results
    std::map<std::string, Systematic<HistType> > fUnfoldedShiftedXSec;
  };
}
