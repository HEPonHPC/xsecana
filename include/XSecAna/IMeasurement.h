#pragma once

#include "TDirectory.h"
#include "XSecAna/_Hist.h"
#include "XSecAna/Type.h"

namespace xsec {
    class IMeasurement {
    public:
        virtual const _hist * Eval(const _hist * data) const = 0;

        virtual void SaveTo(TDirectory * dir, const std::string & subdir) const = 0;

        static std::unique_ptr<IMeasurement>
        LoadFrom(xsec::type::LoadFunction<IMeasurement> load,
                 TDirectory * dir,
                 const std::string & name) {
            return load(dir, name);
        }

        virtual ~IMeasurement() = default;
    private:

    };
}
