#pragma once

#include "TDirectory.h"

#include "XSecAna/Hist.h"
#include "XSecAna/Type.h"
#include "XSecAna/IFlux.h"

namespace xsec {

    template<class HistType>
    class SimpleFlux : public IFlux<HistType> {
    public:
        SimpleFlux() = default;

        explicit SimpleFlux(const HistType & flux)
                : fFlux(flux) {}

        virtual HistType Eval(const Array & edges_and_uof) const override
        { return fFlux; }

        void SaveTo(TDirectory * dir, std::string subdir) const override;

        static std::unique_ptr<IFlux<HistType>>
        LoadFrom(TDirectory * dir, const std::string & subdir);

    protected:
        HistType fFlux;
    };

    template<class HistType>
    class SimpleIntegratedFlux : public IFlux<HistType> {
    public:
        SimpleIntegratedFlux() = default;
        explicit SimpleIntegratedFlux(const Hist & flux)
                : fFlux(flux)
        {}

        virtual HistType Eval(const typename HistType::edges_type & edges_and_uof) const override;

        void SaveTo(TDirectory * dir, std::string subdir) const override;

        static std::unique_ptr<IFlux<HistType>>
        LoadFrom(TDirectory * dir, const std::string & subdir);
    private:
        Hist fFlux;
    };

    //////////////////////////////////////////////////////////
    template<class HistType>
    void
    SimpleFlux<HistType>::
    SaveTo(TDirectory * dir, std::string subdir) const {
        TDirectory * tmp = gDirectory;
        dir = dir->mkdir(subdir.c_str());
        dir->cd();

        TObjString("SimpleFlux").Write("type");
        fFlux.SaveTo(dir, "fFlux");

        tmp->cd();
    }

    //////////////////////////////////////////////////////////
    template<class HistType>
    std::unique_ptr<IFlux<HistType> >
    SimpleFlux<HistType>::
    LoadFrom(TDirectory * dir, const std::string & subdir) {
        dir = dir->GetDirectory(subdir.c_str());

        // make sure we're loading the right type
        auto ptag = (TObjString *) dir->Get("type");
        assert(ptag->GetString() == "SimpleFlux" && "Type does not match SimpleFlux");
        delete ptag;

        HistType flux = *HistType::LoadFrom(dir, "fFlux");
        return std::make_unique<SimpleFlux<HistType>>(flux);
    }

    //////////////////////////////////////////////////////////
    template<class HistType>
    HistType
    SimpleIntegratedFlux<HistType>::
    Eval(const typename HistType::edges_type & edges_and_uof) const {
        auto N = fFlux.Integrate();
        auto ones = HistType::array_and_uof_type::Ones(edges_and_uof.size() - 1);
        return HistType(ones * N,
                        edges_and_uof,
                        ones * std::sqrt(N),
                        this->fFlux.Exposure());
    }

    //////////////////////////////////////////////////////////
    template<class HistType>
    void
    SimpleIntegratedFlux<HistType>::
    SaveTo(TDirectory * dir, std::string subdir) const {
        TDirectory * tmp = gDirectory;
        dir = dir->mkdir(subdir.c_str());
        dir->cd();

        TObjString("SimpleIntegratedFlux").Write("type");
        fFlux.SaveTo(dir, "fFlux");

        tmp->cd();
    }

    //////////////////////////////////////////////////////////
    template<class HistType>
    std::unique_ptr<IFlux<HistType> >
    SimpleIntegratedFlux<HistType>::
    LoadFrom(TDirectory * dir, const std::string & subdir) {
        dir = dir->GetDirectory(subdir.c_str());

        // make sure we're loading the right type
        auto ptag = (TObjString *) dir->Get("type");
        assert(ptag->GetString() == "SimpleIntegratedFlux" && "Type does not match SimpleFlux");
        delete ptag;

        auto flux = *Hist::LoadFrom(dir, "fFlux");
        return std::make_unique<SimpleIntegratedFlux<HistType>>(flux);
    }
}
