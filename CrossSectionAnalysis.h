#pragma once

#include <string>
#include <map>

#include "TGraphAsymmErrors.h"
#include "TH1.h"
#include "TVector3.h"

#include "CAFAna/Core/Cut.h"
#include "CAFAna/Core/HistAxis.h"
#include "CAFAna/Core/SystShifts.h"

#include "XSecAna/Systematic.h"

namespace xsec {
  template<class CrossSectionType,
	   class UnfoldType,
	   class UncertaintyPropogator>
  class CrossSectionAnalysis {

  public:
    CrossSectionAnalysis(const ana::NuTruthCut & signal_cut,
			 const ana::Cut & selection_cut,
			 const ana::NuTruthHistAxis & truth_axis,
			 const ana::HistAxis & reco_axis)
      : fSignalCut(new ana::NuTruthCut(signal_cut)),
	fSelectionCut(new ana::Cut(selection_cut)),
	fTruthAxis(new ana::NuTruthHistAxis(truth_axis)),
	fRecoAxis(new ana::HistAxis(reco_axis))
    {}
			  

    /// \brief Forward to UncertaintyPropogator
    TH1 * RelativeXSecUncertainty(std::string syst_name);

    /// \brief Forward to UncertaintyPropogator
    std::pair<TH1*, TH1*> TotalXSecUncertainty();
    
    TH1 * UnfoldedCrossSection(std::string syst_name);

    /// \brief call Go on all loaders
    void Go();

    ~CrossSectionAnalysis();
    void SaveTo(TDirectory * dir, std::string subdir) const;
    static std::unique_ptr<CrossSectionAnalysis> LoadFrom(TDirectory * dir, std::string name);


  protected:
    ///\brief constructor for loading analysis from file
    CrossSectionAnalysis(CrossSectionType * nominal_xsec,
			 std::map<std::string, Systematic<CrossSectionType> > shifted_xsec,
			 UnfoldType * unfold,
			 ana::Spectrum * data)
      : fNominalXSec(nominal_xsec),
	fShiftedXSec(shifted_xsec),
	fUnfold(unfold),
	fData(data)
    {}

    // nominal is special 
    CrossSectionType * fNominalXSec;

    std::map<std::string, Systematic<CrossSectionType> > fShiftedXSec;

    UnfoldType * fUnfold = NULL;
    const ana::Spectrum * fData  = NULL;

    // cache the nominal unfolded results
    TH1 * fUnfoldedNominalXSec;

    // cache the unfolded shifted results
    std::map<std::string, Systematic<TH1> > fUnfoldedShiftedXSec;

    const ana::Cut * fSelectionCut          = NULL;
    const ana::HistAxis * fRecoAxis         = NULL;
    const ana::NuTruthCut * fSignalCut      = NULL;
    const ana::NuTruthHistAxis * fTruthAxis = NULL;


  };
}
