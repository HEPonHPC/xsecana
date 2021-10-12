#pragma once

#include "TDirectory.h"


#include "XSecAna/Type.h"
#include "XSecAna/IFlux.h"
#include "XSecAna/Utils.h"

namespace xsec {

    class SimpleFlux : public IFluxEstimator {
    public:
        SimpleFlux() = default;

        explicit SimpleFlux(const TH1 * flux);

        virtual const TH1 * Eval(const TH1 * data) const override;

        virtual void SaveTo(TDirectory * dir, const std::string & subdir) const override;

        static std::unique_ptr<IMeasurement>
        LoadFrom(TDirectory * dir, const std::string & subdir);

    protected:
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
    };
}
