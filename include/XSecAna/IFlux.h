#pragma once

#include "XSecAna/Hist.h"
#include "XSecAna/Type.h"

namespace xsec {
    template<class HistType>
    class IFlux {
    public:
        virtual HistType Eval(const Array & edges_and_uof) const = 0;

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
