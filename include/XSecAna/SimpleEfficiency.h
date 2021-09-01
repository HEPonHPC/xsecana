#pragma once

#include "TDirectory.h"

#include "XSecAna/Hist.h"
#include "XSecAna/IEfficiency.h"

namespace xsec {

    template< class HistType = HistXd >
    class SimpleEfficiency : public IEfficiency< HistType > {
    public:
        SimpleEfficiency(HistType num,
                         HistType den)
                : fNumerator(num), fDenominator(den) {}

        const HistType & ToHist() override;

        void SaveTo(TDirectory * dir, std::string subdir) const override;

        static std::unique_ptr< SimpleEfficiency > LoadFrom(TDirectory * dir, std::string name);

        const HistType & GetNumerator() const { return fNumerator; }

        const HistType & GetDenominator() const { return fDenominator; }

    private:
        HistType fNumerator;
        HistType fDenominator;

        // cache the ratio
        HistType * fRatio = 0;
    };

    //////////////////////////////////////////////////////////
    template< class HistType >
    const HistType &
    SimpleEfficiency< HistType >::
    ToHist() {
        if (!fRatio) {
            fRatio = new HistType(fNumerator);
            *fRatio /= fDenominator;
        }
        return *fRatio;
    }

    //////////////////////////////////////////////////////////
    template< class HistType >
    void
    SimpleEfficiency< HistType >::
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
    template< class HistType >
    std::unique_ptr< SimpleEfficiency< HistType > >
    SimpleEfficiency< HistType >::
    LoadFrom(TDirectory * dir, std::string subdir) {
        dir = dir->GetDirectory(subdir.c_str());

        // make sure we're loading the right type
        TObjString * ptag = (TObjString *) dir->Get("type");
        assert(ptag->GetString() == "SimpleEfficiency" && "Type does not match SimpleEfficiency");
        delete ptag;

        HistType numerator = *HistType::LoadFrom(dir, "fNumerator").release();
        HistType denominator = *HistType::LoadFrom(dir, "fDenominator").release();
        return std::make_unique< SimpleEfficiency< HistType > >(numerator,
                                                                denominator);
    }
}
