#pragma once

#include "TDirectory.h"
#include "XSecAna/IMeasurement.h"

#include "XSecAna/Type.h"

namespace xsec {

    ///\brief Define an interface for loading, and calculating an efficiency
    class IEfficiency {
    public:
        static std::unique_ptr<IEfficiency>
        LoadFrom(xsec::type::LoadFunction<IEfficiency> load,
                 TDirectory * dir,
                 const std::string & name) {
            return load(dir, name);
        }

    protected:
    };

    class IEfficiencyEstimator : public IEfficiency, public IMeasurement{};
    class IEigenEfficiencyEstimator : public virtual IEfficiency, public virtual IMeasurement, public IEigenEval{};

}
