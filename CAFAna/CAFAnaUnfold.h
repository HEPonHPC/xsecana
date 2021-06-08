#pragma once
#include "XSecAna/CAFAna/CAFAnaInterface.h"

#include "XSecAna/IUnfold.h"
#include "XSecAna/Hist.h"

#include "TDirectory.h"

namespace xsec {
  template<class CAFAnaUnfoldType,
	   int UnfoldReg,
	   class HistType=HistXd>
  class CAFAnaUnfold : IUnfold<HistType> {
  public:
    HistType Truth(const HistType & reco) const;

    void SaveTo(TDirectory * dir, const std::string & subdir) const;
    static std::unique_ptr<CAFAnaUnfold> LoadFrom(TDirectory * dir, const std::string & subdir);
    
  private:
    HistXd * fMat;
    CAFAnaUnfoldType * fUnfold;
  };

  ///////////////////////////////////////////////////////
  template<class CAFAnaUnfoldType,
	   int UnfoldReg,
	   class HistType>
  HistType 
  CAFAnaUnfold<CAFAnaUnfoldType,
	       UnfoldReg, 
	       HistType>::
  Truth(const HistType & reco) const
  {
    if(!fUnfold) fUnfold = new CAFAnaUnfoldType(cafana::ToReweightableSpectrum(*fMat), UnfoldReg);
    return cafana::ToHist<HistType>(fUnfold->Truth(cafana::ToSpectrum(*reco)));
  }


  ///////////////////////////////////////////////////////
  template<class CAFAnaUnfoldType,
	   int UnfoldReg,
	   class HistType>
  void
  CAFAnaUnfold<CAFAnaUnfoldType,
	       UnfoldReg, 
	       HistType>::
  SaveTo(TDirectory * dir, const std::string & subdir) const
  {
    TDirectory * tmp = gDirectory;
    dir = dir->GetDirectory(subdir.c_str());
    dir->cd();
    
    TObjString("CAFAnaUnfold").Write("type");
    fMat->SaveTo(dir, "fMat");

    delete dir;
    tmp->cd();
  }
  
  ///////////////////////////////////////////////////////
  template<class CAFAnaUnfoldType,
	   int UnfoldReg,
	   class HistType>
  std::unique_ptr<CAFAnaUnfold<CAFAnaUnfoldType,
			       UnfoldReg, 
			       HistType> >
  CAFAnaUnfold<CAFAnaUnfoldType,
	       UnfoldReg, 
	       HistType>::
  LoadFrom(TDirectory * dir, const std::string & subdir)
  {
    dir = dir->GetDirectory(subdir.c_str());

    // make sure we're loading the right type
    TObjString * ptag = (TObjString*) dir->Get("type");
    assert(ptag->GetString() == "CAFAnaUnfold" && "Type does not match CAFAnaUnfold");
    delete ptag;
  }


}
