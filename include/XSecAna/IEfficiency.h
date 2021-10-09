#pragma once

#include "TDirectory.h"

#include "XSecAna/_Hist.h"
#include "XSecAna/Type.h"

namespace xsec {

    ///\brief Define an interface for saving, loading, and calculating an efficiency
    class IEfficiency {

    public:
        ///\brief return the calculated efficiency
        virtual _hist * Eval() = 0;

        virtual void SaveTo(TDirectory * dir, std::string subdir) const = 0;

        static std::unique_ptr<IEfficiency>
        LoadFrom(xsec::type::LoadFunction<IEfficiency> load,
                 TDirectory * dir,
                 const std::string & name) {
            return load(dir, name);
        }

        virtual ~IEfficiency() = default;

    };

}
