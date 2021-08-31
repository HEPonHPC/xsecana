#pragma once

#include "XSecAna/Hist.h"

namespace xsec {
  template<class HistType = HistXd>
  class IUnfold {
  public:
    virtual HistType Truth(const HistType & reco) const = 0; 
    virtual void SaveTo(TDirectory * dir, const std::string& name) const = 0;

    /// \brief Children must override this function
    static std::unique_ptr<IUnfold> LoadFrom(TDirectory * dir, const std::string& name )
    {
      assert(false && "IUnfold::LoadFrom not implemented");
    }

  private:
    
  };  

}
