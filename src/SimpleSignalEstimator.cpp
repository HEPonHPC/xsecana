//
// Created by Derek Doyle on 10/9/21.
//

#include "XSecAna/SimpleSignalEstimator.h"
#include "XSecAna/Math.h"

#include "TDirectory.h"

namespace xsec {
    SimpleSignalEstimator::
    SimpleSignalEstimator(const TH1 * bkgd)
            : fBackground(bkgd) {}


    //////////////////////////////////////////////////////////
    TH1 *
    SimpleSignalEstimator::
    Background(const TH1 * data) const {
        return (TH1*) fBackground->Clone();
    }

    //////////////////////////////////////////////////////////
    void
    SimpleSignalEstimator::
    _eval_impl(const Array & data, const Array & error, ArrayRef result, ArrayRef rerror) const {
        result = data - root::MapContentsToEigen(fBackground);
        rerror = ((error / data).pow(2) +
                  (root::MapErrorsToEigen(fBackground) /
                   root::MapContentsToEigen(fBackground)).pow(2)).sqrt() * result;
        rerror = QuadSum(error, root::MapErrorsToEigen(fBackground));
    }

    //////////////////////////////////////////////////////////
    TH1 *
    SimpleSignalEstimator::
    Signal(const TH1 * data) const {
        return (TH1*) this->Eval(data)->Clone();
    }

    //////////////////////////////////////////////////////////
    void
    SimpleSignalEstimator::
    SaveTo(TDirectory * dir, const std::string & subdir) const {
        TDirectory * tmp = gDirectory;
        dir = dir->mkdir(subdir.c_str());
        dir->cd();

        TObjString("SimpleSignalEstimator").Write("type");
        fBackground->Write("fBackground");

        tmp->cd();
    }

    //////////////////////////////////////////////////////////
    std::unique_ptr<IMeasurement>
    SimpleSignalEstimator::
    LoadFrom(TDirectory * dir, const std::string & subdir) {
        dir = dir->GetDirectory(subdir.c_str());

        // make sure we're loading the right type
        auto ptag = (TObjString *) dir->Get("type");
        assert(ptag->GetString() == "SimpleSignalEstimator" && "Type does not match SimpleSignalEstimator");
        delete ptag;

        auto background = root::LoadTH1(dir, "fBackground").release();
        return std::make_unique<SimpleSignalEstimator>(background);
    }


}
