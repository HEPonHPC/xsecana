#pragma once

#include "XSecAna/Hist.h"
#include "XSecAna/ISignalEstimator.h"

namespace xsec {
    template<class HistType>
    class SimpleSignalEstimator : public ISignalEstimator<HistType> {
    public:
        SimpleSignalEstimator() = default;

        explicit SimpleSignalEstimator(const HistType & bkgd)
                : fBackground(bkgd) {}

        /// \brief An implementation of Eval allows this object to
        /// interact with the systematics framework. This one just forwards a call to Signal
        /// TODO don't really like this
        HistType Eval(const HistType & data) override { return Signal(data); }

        /// \brief background could be dependent on data
        /// in this case it isn't
        const HistType & Background(const HistType & data) override;

        const HistType & Signal(const HistType & data) override;

        void SaveTo(TDirectory * dir, const std::string & subdir) const override;

        static std::unique_ptr<ISignalEstimator<HistType> >
        LoadFrom(TDirectory * dir, const std::string & subdir);

    private:
        HistType fBackground;

        // cache signal hist
        HistType * fSignal = 0;

    };

    //////////////////////////////////////////////////////////
    template<class HistType>
    const HistType &
    SimpleSignalEstimator<HistType>::
    Background(const HistType & data) {
        return fBackground;
    }

    //////////////////////////////////////////////////////////
    template<class HistType>
    const HistType &
    SimpleSignalEstimator<HistType>::
    Signal(const HistType & data) {
        if (!fSignal) fSignal = new HistType(data - fBackground);
        return *fSignal;
    }

    //////////////////////////////////////////////////////////
    template<class HistType>
    void
    SimpleSignalEstimator<HistType>::
    SaveTo(TDirectory * dir, const std::string & subdir) const {
        TDirectory * tmp = gDirectory;
        dir = dir->mkdir(subdir.c_str());
        dir->cd();

        TObjString("SimpleSignalEstimator").Write("type");
        fBackground.SaveTo(dir, "fBackground");

        tmp->cd();
    }

    //////////////////////////////////////////////////////////
    template<class HistType>
    std::unique_ptr<ISignalEstimator<HistType> >
    SimpleSignalEstimator<HistType>::
    LoadFrom(TDirectory * dir, const std::string & subdir) {
        dir = dir->GetDirectory(subdir.c_str());

        // make sure we're loading the right type
        auto ptag = (TObjString *) dir->Get("type");
        assert(ptag->GetString() == "SimpleSignalEstimator" && "Type does not match SimpleSignalEstimator");
        delete ptag;

        HistType background = *HistType::LoadFrom(dir, "fBackground");
        return std::make_unique<SimpleSignalEstimator<HistType> >(background);
    }


}
