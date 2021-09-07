#pragma once

#include "Hist.h"

namespace xsec {
    template<class HistType = HistXd>
    class IFlux {
    public:
        virtual HistType Eval() = 0;

        virtual HistType operator/(const HistType & rhs) = 0;

        virtual HistType operator*(const HistType & rhs) = 0;

        virtual void SaveTo(TDirectory * dir, std::string subdir) const = 0;

        /// \brief Children must override this function
        static std::unique_ptr<IFlux> LoadFrom(TDirectory * dir, const std::string & name) {
            assert(false && "IFlux::LoadFrom not implemented");
        }

    private:

    };
}
