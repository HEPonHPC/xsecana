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

    HistType * RelativeXSecUncertainty(std::string syst_name);

    /// \brief Forward to UncertaintyPropogator
    std::pair<HistType, HistType> TotalXSecUncertainty();
    
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
