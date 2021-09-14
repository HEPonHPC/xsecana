#pragma once

#include "TDirectory.h"

#include "XSecAna/Hist.h"
#include "XSecAna/Type.h"
#include "XSecAna/IFlux.h"

namespace xsec {

    template<class HistType,
            bool Integrated = false>
    class SimpleFlux : public IFlux<HistType> {
    public:
        SimpleFlux() = default;

        explicit SimpleFlux(const HistType & flux)
                : fFlux(flux) {}

        // TODO this should return fResult instead
        HistType Eval() override { return fFlux; }

        HistType operator/(const HistType & rhs) override;

        HistType operator*(const HistType & rhs) override;

        void SaveTo(TDirectory * dir, std::string subdir) const override;

        static std::unique_ptr<IFlux<HistType> >
        LoadFrom(TDirectory * dir, const std::string & subdir);

    protected:
        // hold the raw histogram
        HistType fFlux;

        // cache result of (possibly) integrated flux
        // if not integrated, this will just point to fFlux
        HistType * fResult = 0;

        HistType * Result(const HistType & rhs);

    };

    // template alias deduction was introduced in c++20
    // otherwise, HistType will not be deduced.
    template<class HistType>
    using SimpleIntegratedFlux = SimpleFlux<HistType, true>;

    //////////////////////////////////////////////////////////
    template<class HistType,
            bool Integrated>
    HistType *
    SimpleFlux<HistType,
               Integrated>::
    Result(const HistType & rhs) {
        if (!fResult) {
            if constexpr(Integrated) {
                fResult = new HistType(
                        std::decay_t<decltype(rhs.Contents())>::Ones(rhs.Contents().size()) * fFlux.Integrate(),
                        rhs.Edges());
            } else {
                fResult = &fFlux;
            }
        }
        return fResult;
    }

    //////////////////////////////////////////////////////////
    template<class HistType,
            bool Integrated>
    HistType
    SimpleFlux<HistType,
               Integrated>::
    operator*(const HistType & rhs) {
        return *Result(rhs) * rhs;
    }

    //////////////////////////////////////////////////////////
    template<class HistType,
            bool Integrated>
    HistType
    SimpleFlux<HistType,
               Integrated>::
    operator/(const HistType & rhs) {
        return *Result(rhs) / rhs;
    }

    //////////////////////////////////////////////////////////
    template<class HistType,
            bool Integrated>
    void
    SimpleFlux<HistType,
               Integrated>::
    SaveTo(TDirectory * dir, std::string subdir) const {
        TDirectory * tmp = gDirectory;
        dir = dir->mkdir(subdir.c_str());
        dir->cd();

        TObjString("SimpleFlux").Write("type");
        fFlux.SaveTo(dir, "fFlux");

        tmp->cd();
    }

    //////////////////////////////////////////////////////////
    template<class HistType,
            bool Integrated>
    std::unique_ptr<IFlux<HistType> >
    SimpleFlux<HistType,
               Integrated>::
    LoadFrom(TDirectory * dir, const std::string & subdir) {
        dir = dir->GetDirectory(subdir.c_str());

        // make sure we're loading the right type
        auto ptag = (TObjString *) dir->Get("type");
        assert(ptag->GetString() == "SimpleFlux" && "Type does not match SimpleFlux");
        delete ptag;

        HistType flux = *HistType::LoadFrom(dir, "fFlux");
        return std::make_unique<SimpleFlux<HistType, Integrated>>(flux);
    }
}
