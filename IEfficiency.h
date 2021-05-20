#pragma once

#include "TDirectory.h"

#include "XSecAna/Hist.h"

namespace xsec {
  
  ///\brief Define an interface for saving, loading, and calculating an efficiency
  template<class HistType = HistXXd>
  class IEfficiency {

  public:
    ///\brief return the calculated efficiency
    virtual HistType * Efficiency() = 0;
    virtual void SaveTo(TDirectory * dir, std::string subdir) const = 0;

    /// \brief children must override
    static std::unique_ptr<IEfficiency> * LoadFrom(TDirectory * dir, const std::string & name) 
    {
      assert(false && "IEfficiency::LoadFrom not implemented");
    }
    
  };

}
