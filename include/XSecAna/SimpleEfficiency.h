#pragma once

#include "TDirectory.h"

#include "XSecAna/Hist.h"
#include "XSecAna/IEfficiency.h"

namespace xsec {

    class SimpleEfficiency : public IEigenEfficiencyEstimator {
    public:
        SimpleEfficiency(const TH1 * num,
                         const TH1 * den);

        void SaveTo(TDirectory * dir, const std::string & subdir) const override;

        static std::unique_ptr<IMeasurement> LoadFrom(TDirectory * dir, const std::string & subdir);

        const TH1 * GetNumerator() const;

        const TH1 * GetDenominator() const;

        virtual ~SimpleEfficiency() = default;

    private:
        virtual void _eval_impl(const Array & data, const Array & error,
                                ArrayRef result, ArrayRef rerror) const override;
        const TH1 * fNumerator;
        const TH1 * fDenominator;
    };
}
