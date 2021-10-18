#pragma once

#include "XSecAna/ISignalEstimator.h"
#include "XSecAna/IMeasurement.h"

namespace xsec {
    class SimpleSignalEstimator : public IEigenSignalEstimator {
    public:
        SimpleSignalEstimator() = default;

        explicit SimpleSignalEstimator(const TH1 * bkgd);

        /// \brief background could be dependent on data
        /// in this case it isn't
        TH1 * Background(const TH1 * data) const override;

        TH1 * Signal(const TH1 * data) const override;

        void SaveTo(TDirectory * dir, const std::string & subdir) const override;

        static std::unique_ptr<IMeasurement>
        LoadFrom(TDirectory * dir, const std::string & subdir);

    private:
        void _eval_impl(const Array & data, const Array & error, ArrayRef result, ArrayRef rerror) const override;
        const TH1 * fBackground;
    };
}
