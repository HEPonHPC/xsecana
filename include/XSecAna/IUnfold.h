#pragma once

#include "XSecAna/Hist.h"

namespace xsec {
  template<class HistType = HistXd>
  class IUnfold {
  public:
    virtual HistType Truth(const HistType & reco) const = 0; 
  private:
    
  };  

}
