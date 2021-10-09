#pragma once

#include "XSecAna/_Hist.h"
#include "XSecAna/ISignalEstimator.h"

namespace xsec {
    class SimpleSignalEstimator : public ISignalEstimator {
    public:
        SimpleSignalEstimator() = default;

        explicit SimpleSignalEstimator(const _hist * bkgd)
                : fBackground(bkgd) {}

        /// \brief An implementation of Eval allows this object to
        /// interact with the systematics framework. This one just forwards a call to Signal
        /// TODO don't really like this
        const _hist * Eval(const _hist * data) override { return Signal(data); }

        /// \brief background could be dependent on data
        /// in this case it isn't
        const _hist * Background(const _hist * data) override;

        const _hist * Signal(const _hist * data) override;

        void SaveTo(TDirectory * dir, const std::string & subdir) const override;

        static std::unique_ptr<ISignalEstimator>
        LoadFrom(TDirectory * dir, const std::string & subdir);

    private:
        const _hist * fBackground;

        // cache signal _hist
        _hist * fSignal = 0;

    };

    //////////////////////////////////////////////////////////
    const _hist *
    SimpleSignalEstimator::
    Background(const _hist * data) {
        return fBackground;
    }

    //////////////////////////////////////////////////////////
    const _hist *
    SimpleSignalEstimator::
    Signal(const _hist * data) {
        if (!fSignal) fSignal = data->Subtract(fBackground);
        return fSignal;
    }

    //////////////////////////////////////////////////////////
    void
    SimpleSignalEstimator::
    SaveTo(TDirectory * dir, const std::string & subdir) const {
        TDirectory * tmp = gDirectory;
        dir = dir->mkdir(subdir.c_str());
        dir->cd();

        TObjString("SimpleSignalEstimator").Write("type");
        fBackground->SaveTo(dir, "fBackground");

        tmp->cd();
    }

    //////////////////////////////////////////////////////////
    std::unique_ptr<ISignalEstimator>
    SimpleSignalEstimator::
    LoadFrom(TDirectory * dir, const std::string & subdir) {
        dir = dir->GetDirectory(subdir.c_str());

        // make sure we're loading the right type
        auto ptag = (TObjString *) dir->Get("type");
        assert(ptag->GetString() == "SimpleSignalEstimator" && "Type does not match SimpleSignalEstimator");
        delete ptag;

        _hist * background = _hist::LoadFrom(dir, "fBackground").release();
        return std::make_unique<SimpleSignalEstimator>(background);
    }


}
