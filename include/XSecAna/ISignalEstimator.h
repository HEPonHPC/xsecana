#pragma once

#include "TDirectory.h"
#include "XSecAna/Type.h"

namespace xsec {
    /// Defining interface for SignalEstimators
    template<class HistType = HistXd>
    class ISignalEstimator {
    public:
        virtual HistType Eval(const HistType & data) = 0;

        virtual const HistType & Background(const HistType & data) = 0;

        virtual const HistType & Signal(const HistType & data) = 0;

        virtual void SaveTo(TDirectory * dir, const std::string & name) const = 0;

        static std::unique_ptr<ISignalEstimator>
        LoadFrom(xsec::type::LoadFunction<ISignalEstimator> load,
                 TDirectory * dir,
                 const std::string & name) {
            return load(dir, name);
        }

        virtual ~ISignalEstimator() = default;
    protected:
    };
}
