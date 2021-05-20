#pragma once

#include "Hist.h"

namespace xsec {
  template<class HistType = HistXXd>
  class IFlux {
  public:
    virtual HistType * Flux() = 0; 
    virtual void SaveTo(TDirectory * dir, std::string subdir) const = 0;
  private:
    
  };
}
