#pragma once

#include "TDirectory.h"
#include "Hist.h"
#include "Type.h"

namespace xsec {
    class IMeasurement {
    public:
        virtual Hist Eval(const Hist & data) const = 0;

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
