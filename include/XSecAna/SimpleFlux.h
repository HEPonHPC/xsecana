#pragma once

#include "TDirectory.h"


#include "XSecAna/Type.h"
#include "XSecAna/IFlux.h"
#include "XSecAna/Utils.h"

namespace xsec {

    class SimpleFlux : public IEigenFluxEstimator {
    public:
        SimpleFlux() = default;

        explicit SimpleFlux(const TH1 * flux);

        virtual void SaveTo(TDirectory * dir, const std::string & subdir) const override;

        IMeasurement * IntegratedFlux() const;

        static std::unique_ptr<IMeasurement>
        LoadFrom(TDirectory * dir, const std::string & subdir);

    protected:
        virtual void _eval_impl(const Array & data, const Array & error,
                                ArrayRef result, ArrayRef rerror) const override;
        const TH1 * fFlux;
    };

    class SimpleIntegratedFlux : public IEigenFluxEstimator {
    public:
        SimpleIntegratedFlux() = default;

        explicit SimpleIntegratedFlux(const TH1 * flux);

        static std::unique_ptr<IMeasurement>
        LoadFrom(TDirectory * dir, const std::string & subdir);

        virtual void SaveTo(TDirectory * dir, const std::string & subdir) const override;

    private:
        virtual void _eval_impl(const Array & data, const Array & error,
                                ArrayRef result, ArrayRef rerror) const override;
        const TH1 * fFlux;
        const double fN=1;
        const double fdN=1;
    };
}
