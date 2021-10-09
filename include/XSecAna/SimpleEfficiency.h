#pragma once

#include "TDirectory.h"

#include "XSecAna/Hist.h"
#include "XSecAna/IEfficiency.h"

namespace xsec {

    class SimpleEfficiency : public IEfficiency {
    public:
        SimpleEfficiency(_hist * num,
                         _hist * den)
                : fNumerator(num), fDenominator(den) {}

        _hist * Eval() override;

        void SaveTo(TDirectory * dir, std::string subdir) const override;

        static std::unique_ptr<IEfficiency> LoadFrom(TDirectory * dir, const std::string & subdir);

        const _hist * GetNumerator() const { return fNumerator; }

        const _hist * GetDenominator() const { return fDenominator; }

    private:
        const _hist * fNumerator;
        const _hist * fDenominator;

        // cache the ratio
        _hist * fRatio = 0;
    };

    //////////////////////////////////////////////////////////
    _hist *
    SimpleEfficiency::
    Eval() {
        if (!fRatio) {
            fRatio = fNumerator->Divide(fDenominator);

            // binomial error
            auto e = fRatio->GetContentsAndUOF();
            auto nsig = fDenominator->GetContentsAndUOF();
            fRatio->SetErrorsAndUOF((e * (1 - e) / nsig).sqrt());
        }
        return fRatio;
    }

    //////////////////////////////////////////////////////////
    void
    SimpleEfficiency::
    SaveTo(TDirectory * dir, std::string subdir) const {
        TDirectory * tmp = gDirectory;
        dir = dir->mkdir(subdir.c_str());
        dir->cd();

        TObjString("SimpleEfficiency").Write("type");
        fNumerator->SaveTo(dir, "fNumerator");
        fDenominator->SaveTo(dir, "fDenominator");

        tmp->cd();
    }

    //////////////////////////////////////////////////////////
    std::unique_ptr<IEfficiency>
    SimpleEfficiency::
    LoadFrom(TDirectory * dir, const std::string & subdir) {
        dir = dir->GetDirectory(subdir.c_str());

        // make sure we're loading the right type
        auto ptag = (TObjString *) dir->Get("type");
        assert(ptag->GetString() == "SimpleEfficiency" && "Type does not match SimpleEfficiency");
        delete ptag;

        _hist * numerator = _hist::LoadFrom(dir, "fNumerator").release();
        _hist * denominator = _hist::LoadFrom(dir, "fDenominator").release();
        return std::make_unique<SimpleEfficiency>(numerator,
                                                  denominator);
    }
}
