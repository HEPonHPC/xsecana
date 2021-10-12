#pragma once

#include "XSecAna/Type.h"
#include "XSecAna/IMeasurement.h"

namespace xsec {
    class IFlux {
    public:
        static std::unique_ptr<IFlux>
        LoadFrom(xsec::type::LoadFunction<IFlux> load,
                 TDirectory * dir,
                 const std::string & name) {
            return load(dir, name);
        }
    private:
    };

    class IFluxEstimator : public IFlux, public IMeasurement{};
    class IEigenFluxEstimator : public IFlux, public virtual IMeasurement, public IEigenEval{};
}
