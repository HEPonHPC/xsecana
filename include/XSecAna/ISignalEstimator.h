#pragma once

#include "TDirectory.h"
#include "TH1.h"
#include "XSecAna/Type.h"
#include "XSecAna/IMeasurement.h"

namespace xsec {
    /// Defining interface for SignalEstimators
    class ISignal {
    public:
        virtual TH1 * Background(const TH1 * data) const = 0;

        virtual TH1 * Signal(const TH1 * data) const = 0;

        static std::unique_ptr<ISignal>
        LoadFrom(xsec::type::LoadFunction<ISignal> load,
                 TDirectory * dir,
                 const std::string & name) {
            return load(dir, name);
        }

        virtual ~ISignal() = default;

    protected:

    };

    class ISignalEstimator : public ISignal, public IMeasurement{};
    class IEigenSignalEstimator : public ISignal, public virtual IMeasurement, public IEigenEval{};
}
