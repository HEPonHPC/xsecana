#pragma once

#include "TDirectory.h"

#include "XSecAna/Hist.h"
#include "XSecAna/IFlux.h"

namespace xsec {
  template<class HistType = HistXXd>
  class SimpleFlux : IFlux<HistType> {
  public:
    HistType * Flux() override;
    void SaveTo(TDirectory * dir, std::string subdir) const override;
    static std::unique_ptr<SimpleFlux<HistType> > LoadFrom(TDirectory * dir, std::string subdir);
  private:
    // hold the raw histogram
    HistType * fFlux;
  };
}
