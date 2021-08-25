#pragma once

#include "Hist.h"

namespace xsec {
  template<class HistType = HistXd>
  class IFlux {
  public:
    virtual const HistType & ToHist() = 0; 
    virtual HistType operator/(const HistType & rhs) = 0;
    virtual HistType operator*(const HistType & rhs) = 0;
    virtual void SaveTo(TDirectory * dir, std::string subdir) const = 0;
  private:
    
  };
}
