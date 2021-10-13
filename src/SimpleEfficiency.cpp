//
// Created by Derek Doyle on 10/9/21.
//

#include "XSecAna/SimpleEfficiency.h"

#include "TDirectory.h"

namespace xsec {
    SimpleEfficiency::
    SimpleEfficiency(const TH1 * num,
                     const TH1 * den)
            : fNumerator(num), fDenominator(den) {}

    const TH1 *
    SimpleEfficiency::
    GetNumerator() const { return fNumerator; }

    const TH1 *
    SimpleEfficiency::
    GetDenominator() const { return fDenominator; }

    //////////////////////////////////////////////////////////
    void
    SimpleEfficiency::
    _eval_impl(const Array & data, const Array & error, ArrayRef result, ArrayRef rerror) const {

        // binomial error
        result = root::MapContentsToEigen(fNumerator) / root::MapContentsToEigen(fDenominator);
        rerror = (result * (1 - result) / root::MapContentsToEigen(fDenominator)).sqrt();
    }

//////////////////////////////////////////////////////////
    void
    SimpleEfficiency::
    SaveTo(TDirectory * dir, const std::string & subdir) const {
        TDirectory * tmp = gDirectory;
        dir = dir->mkdir(subdir.c_str());
        dir->cd();

        TObjString("SimpleEfficiency").Write("type");
        fNumerator->Write("fNumerator");
        fDenominator->Write("fDenominator");

        tmp->cd();
    }

//////////////////////////////////////////////////////////
    std::unique_ptr<IMeasurement>
    SimpleEfficiency::
    LoadFrom(TDirectory * dir, const std::string & subdir) {
        dir = dir->GetDirectory(subdir.c_str());

        // make sure we're loading the right type
        auto ptag = (TObjString *) dir->Get("type");
        assert(ptag->GetString() == "SimpleEfficiency" && "Type does not match SimpleEfficiency");
        delete ptag;

        auto numerator = root::LoadTH1(dir, "fNumerator").release();
        auto denominator = root::LoadTH1(dir, "fDenominator").release();
        return std::make_unique<SimpleEfficiency>(numerator, denominator);
    }

}
