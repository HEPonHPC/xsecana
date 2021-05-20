#pragma once

#include "TDirectory.h"

namespace xsec
{
  /// Defining interface for SignalEstimators
  template<class HistType = HistXXd>
  class ISignalEstimator
  {
  public:
    virtual const HistType * Background(const HistType * data) const = 0;
    virtual const HistType * Signal(const HistType *) const = 0;
    virtual void SaveTo(TDirectory * dir, const std::string& name) const = 0;

    /// \brief Children must override this function
    static std::unique_ptr<ISignalEstimator> LoadFrom(TDirectory * dir, const std::string& name )
    {
      assert(false && "ISignalEstimator::LoadFrom not implemented");
    }

  protected:
  };
}
