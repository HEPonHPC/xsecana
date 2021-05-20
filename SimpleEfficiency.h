#pragma once

#include "TDirectory.h"

#include "XSecAna/Hist.h"
#include "XSecAna/IEfficiency.h"

namespace xsec {

  template<class HistType = HistXXd>
  class SimpleEfficiency : public IEfficiency<HistType> {
  public:
    HistType * Efficiency() override;
    void SaveTo(TDirectory * dir, std::string subdir) const override;
    
    static std::unique_ptr<SimpleEfficiency> LoadFrom(TDirectory * dir, std::string name);

  private:
    HistType * fNumerator;
    HistType * fDenominator;
    HistType * fRatio;
  };
}
