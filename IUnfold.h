#pragma once

namespace xsec {
  template<class HistType = HistXXd>
  class IUnfold {
  public:
    virtual HistType * Truth(const HistType * reco) const = 0; 
  private:
    
  };
}
