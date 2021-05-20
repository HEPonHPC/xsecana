#pragma once

#include "XSecAna/Hist.h"
#include "XSecAna/ISignalEstimator.h"

namespace xsec {
  template<class HistType = HistXXd>
  class SimpleSignalEstimator : ISignalEstimator<HistType> {
  public:
    const HistType * Background(const HistType * data) const override;
    const HistType * Signal(const HistType * data) const override;

    void SaveTo(TDirectory * dir, const std::string & name) const;
    static std::unique_ptr<SimpleSignalEstimator<HistType> > LoadFrom(TDirectory * dir, const std::string & subdir);
    
  private:
    HistType * fBackground;
    
  };
}
