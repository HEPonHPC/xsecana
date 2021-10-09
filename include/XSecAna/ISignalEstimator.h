#pragma once

#include "TDirectory.h"
#include "XSecAna/Type.h"

namespace xsec {
    /// Defining interface for SignalEstimators
    class ISignalEstimator {
    public:
        virtual const _hist * Eval(const _hist * data) = 0;

        virtual const _hist * Background(const _hist * data) = 0;

        virtual const _hist * Signal(const _hist * data) = 0;

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
