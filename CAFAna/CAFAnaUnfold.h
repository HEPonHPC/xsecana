#pragma once
#include "XSecAna/CAFAna/CAFAnaInterface.h"

#include "XSecAna/IUnfold.h"
#include "XSecAna/Hist.h"

#include "TDirectory.h"

namespace xsec {
  template<class CAFAnaUnfoldType,
	   int UnfoldReg,
	   class HistType=HistXXd>
  class CAFAnaUnfold : IUnfold<HistType> {
  public:
    HistType * Truth(const HistType * reco) const;

    void SaveTo(TDirectory * dir, const std::string & subdir) const;
    static std::unique_ptr<CAFAnaUnfold> LoadFrom(TDirectory * dir, const std::string & subdir);
    
  private:
    HistXXd * fMat;
    CAFAnaUnfoldType * fUnfold;
  };
}
