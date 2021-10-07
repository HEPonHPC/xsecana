#pragma once

#include "TDirectory.h"

#include "XSecAna/Hist.h"
#include "XSecAna/Type.h"
#include "XSecAna/IFlux.h"

namespace xsec {

    class SimpleFlux : public IFlux {
    public:
        SimpleFlux() = default;

        explicit SimpleFlux(const Hist & flux)
                : fFlux(flux) {}

        virtual Hist Eval(const Array & edges_and_uof) const override { return fFlux; }

        void SaveTo(TDirectory * dir, std::string subdir) const override;

        static std::unique_ptr<IFlux>
        LoadFrom(TDirectory * dir, const std::string & subdir);

    protected:
        Hist fFlux;
    };

    class SimpleIntegratedFlux : public IFlux {
    public:
        SimpleIntegratedFlux() = default;

        explicit SimpleIntegratedFlux(const Hist & flux)
                : fFlux(flux) {}

        virtual Hist Eval(const Array & edges_and_uof) const override;

        void SaveTo(TDirectory * dir, std::string subdir) const override;

        static std::unique_ptr<IFlux>
        LoadFrom(TDirectory * dir, const std::string & subdir);

    private:
        Hist fFlux;
    };

    //////////////////////////////////////////////////////////
    void
    SimpleFlux::
    SaveTo(TDirectory * dir, std::string subdir) const {
        TDirectory * tmp = gDirectory;
        dir = dir->mkdir(subdir.c_str());
        dir->cd();

        TObjString("SimpleFlux").Write("type");
        fFlux.SaveTo(dir, "fFlux");

        tmp->cd();
    }

    //////////////////////////////////////////////////////////
    std::unique_ptr<IFlux>
    SimpleFlux::
    LoadFrom(TDirectory * dir, const std::string & subdir) {
        dir = dir->GetDirectory(subdir.c_str());

        // make sure we're loading the right type
        auto ptag = (TObjString *) dir->Get("type");
        assert(ptag->GetString() == "SimpleFlux" && "Type does not match SimpleFlux");
        delete ptag;

        Hist flux = *Hist::LoadFrom(dir, "fFlux");
        return std::make_unique<SimpleFlux>(flux);
    }

    //////////////////////////////////////////////////////////
    Hist
    SimpleIntegratedFlux::
    Eval(const Array & edges_and_uof) const {
        auto N = fFlux.Integrate();
        auto ones = Array::Ones(edges_and_uof.size() - 1);
        return Hist(ones * N,
                        edges_and_uof,
                        ones * std::sqrt(N),
                        this->fFlux.Exposure());
    }

    //////////////////////////////////////////////////////////
    void
    SimpleIntegratedFlux::
    SaveTo(TDirectory * dir, std::string subdir) const {
        TDirectory * tmp = gDirectory;
        dir = dir->mkdir(subdir.c_str());
        dir->cd();

        TObjString("SimpleIntegratedFlux").Write("type");
        fFlux.SaveTo(dir, "fFlux");

        tmp->cd();
    }

    //////////////////////////////////////////////////////////
    std::unique_ptr<IFlux>
    SimpleIntegratedFlux::
    LoadFrom(TDirectory * dir, const std::string & subdir) {
        dir = dir->GetDirectory(subdir.c_str());

        // make sure we're loading the right type
        auto ptag = (TObjString *) dir->Get("type");
        assert(ptag->GetString() == "SimpleIntegratedFlux" && "Type does not match SimpleFlux");
        delete ptag;

        auto flux = *Hist::LoadFrom(dir, "fFlux");
        return std::make_unique<SimpleIntegratedFlux>(flux);
    }
}
