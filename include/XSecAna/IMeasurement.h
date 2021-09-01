#pragma once

#include "TDirectory.h"
#include "Hist.h"

namespace xsec {
    template<class HistType>
    class IMeasurement {
    public:
        virtual HistType Result(const HistType & data) = 0;

        virtual void SaveTo(TDirectory * dir, const std::string & subdir) const = 0;

        /// \brief Children must override this function
        static std::unique_ptr<IMeasurement> LoadFrom(TDirectory * dir, const std::string & name) {
            assert(false && "IMeasurement::LoadFrom not implemented");
        }

    private:

    };
}
