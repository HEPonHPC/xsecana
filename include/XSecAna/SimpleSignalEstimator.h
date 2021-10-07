#pragma once

#include "XSecAna/Hist.h"
#include "XSecAna/ISignalEstimator.h"

namespace xsec {
    class SimpleSignalEstimator : public ISignalEstimator {
    public:
        SimpleSignalEstimator() = default;

        explicit SimpleSignalEstimator(const Hist & bkgd)
                : fBackground(bkgd) {}

        /// \brief An implementation of Eval allows this object to
        /// interact with the systematics framework. This one just forwards a call to Signal
        /// TODO don't really like this
        Hist Eval(const Hist & data) override { return Signal(data); }

        /// \brief background could be dependent on data
        /// in this case it isn't
        const Hist & Background(const Hist & data) override;

        const Hist & Signal(const Hist & data) override;

        void SaveTo(TDirectory * dir, const std::string & subdir) const override;

        static std::unique_ptr<ISignalEstimator>
        LoadFrom(TDirectory * dir, const std::string & subdir);

    private:
        Hist fBackground;

        // cache signal hist
        Hist * fSignal = 0;

    };

    //////////////////////////////////////////////////////////
    const Hist &
    SimpleSignalEstimator::
    Background(const Hist & data) {
        return fBackground;
    }

    //////////////////////////////////////////////////////////
    const Hist &
    SimpleSignalEstimator::
    Signal(const Hist & data) {
        if (!fSignal) fSignal = new Hist(data - fBackground);
        return *fSignal;
    }

    //////////////////////////////////////////////////////////
    void
    SimpleSignalEstimator::
    SaveTo(TDirectory * dir, const std::string & subdir) const {
        TDirectory * tmp = gDirectory;
        dir = dir->mkdir(subdir.c_str());
        dir->cd();

        TObjString("SimpleSignalEstimator").Write("type");
        fBackground.SaveTo(dir, "fBackground");

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

        Hist background = *Hist::LoadFrom(dir, "fBackground");
        return std::make_unique<SimpleSignalEstimator>(background);
    }


}
