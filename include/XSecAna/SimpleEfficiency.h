#pragma once

#include "TDirectory.h"

#include "XSecAna/Hist.h"
#include "XSecAna/IEfficiency.h"

namespace xsec {

    class SimpleEfficiency : public IEfficiency {
    public:
        SimpleEfficiency(Hist num,
                         Hist den)
                : fNumerator(num), fDenominator(den) {}

        Hist Eval() override;

        void SaveTo(TDirectory * dir, std::string subdir) const override;

        static std::unique_ptr<IEfficiency> LoadFrom(TDirectory * dir, const std::string & subdir);

        const Hist & GetNumerator() const { return fNumerator; }

        const Hist & GetDenominator() const { return fDenominator; }

    private:
        Hist fNumerator;
        Hist fDenominator;

        // cache the ratio
        Hist * fRatio = 0;
    };

    //////////////////////////////////////////////////////////
    Hist
    SimpleEfficiency::
    Eval() {
        if (!fRatio) {
            fRatio = new Hist(fNumerator);
            *fRatio /= fDenominator;

            // binomial error
            auto e = fRatio->ContentsAndUOF();
            auto nsig = fDenominator.ContentsAndUOF();
            fRatio->SetErrorsAndUOF((e * (1 - e) / nsig).sqrt());
        }
        return *fRatio;
    }

    //////////////////////////////////////////////////////////
    void
    SimpleEfficiency::
    SaveTo(TDirectory * dir, std::string subdir) const {
        TDirectory * tmp = gDirectory;
        dir = dir->mkdir(subdir.c_str());
        dir->cd();

        TObjString("SimpleEfficiency").Write("type");
        fNumerator.SaveTo(dir, "fNumerator");
        fDenominator.SaveTo(dir, "fDenominator");

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

        Hist numerator = *Hist::LoadFrom(dir, "fNumerator").release();
        Hist denominator = *Hist::LoadFrom(dir, "fDenominator").release();
        return std::make_unique<SimpleEfficiency>(numerator,
                                                  denominator);
    }
}
