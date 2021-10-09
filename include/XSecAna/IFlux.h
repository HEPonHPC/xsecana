#pragma once

#include "XSecAna/_Hist.h"
#include "XSecAna/Type.h"

namespace xsec {
    class IFlux {
    public:
        virtual const _hist * Eval(const _hist * data) const = 0;

        virtual void SaveTo(TDirectory * dir, std::string subdir) const = 0;

        static std::unique_ptr<IFlux>
        LoadFrom(xsec::type::LoadFunction<IFlux> load,
                 TDirectory * dir,
                 const std::string & name) {
            return load(dir, name);
        }

        virtual ~IFlux() = default;

    private:

    };
}
