#pragma once

#include "TDirectory.h"

#include "XSecAna/_Hist.h"
#include "XSecAna/Type.h"
#include "XSecAna/IFlux.h"

namespace xsec {

    class SimpleFlux : public IFlux {
    public:
        SimpleFlux() = default;

        explicit SimpleFlux(const _hist * flux)
                : fFlux(flux) {}

        virtual const _hist * Eval(const _hist * data) const override { return fFlux; }

        void SaveTo(TDirectory * dir, std::string subdir) const override;

        static std::unique_ptr<IFlux>
        LoadFrom(TDirectory * dir, const std::string & subdir);

    protected:
        const _hist * fFlux;
    };

    class SimpleIntegratedFlux : public IFlux {
    public:
        SimpleIntegratedFlux() = default;

        explicit SimpleIntegratedFlux(const _hist * flux)
                : fFlux(flux) {}

        virtual _hist * Eval(const _hist * data) const override;

        void SaveTo(TDirectory * dir, std::string subdir) const override;

        static std::unique_ptr<IFlux>
        LoadFrom(TDirectory * dir, const std::string & subdir);

    private:
        const _hist * fFlux;
    };

    //////////////////////////////////////////////////////////
    void
    SimpleFlux::
    SaveTo(TDirectory * dir, std::string subdir) const {
        TDirectory * tmp = gDirectory;
        dir = dir->mkdir(subdir.c_str());
        dir->cd();

        TObjString("SimpleFlux").Write("type");
        fFlux->SaveTo(dir, "fFlux");

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

        const auto * flux = _hist::LoadFrom(dir, "fFlux").release();
        return std::make_unique<SimpleFlux>(flux);
    }

    //////////////////////////////////////////////////////////
    _hist *
    SimpleIntegratedFlux::
    Eval(const _hist * data) const {
        auto N = fFlux->Integrate();
        auto ones = Array::Ones(data->GetContentsAndUOF().size());
        auto ret = data->Clone();
        ret->SetContentsAndUOF(ones * N);
        ret->SetContentsAndUOF(ones * std::sqrt(N));
        ret->SetExposure(fFlux->Exposure());
        return ret;
    }

    //////////////////////////////////////////////////////////
    void
    SimpleIntegratedFlux::
    SaveTo(TDirectory * dir, std::string subdir) const {
        TDirectory * tmp = gDirectory;
        dir = dir->mkdir(subdir.c_str());
        dir->cd();

        TObjString("SimpleIntegratedFlux").Write("type");
        fFlux->SaveTo(dir, "fFlux");

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

        const auto flux = _hist::LoadFrom(dir, "fFlux").release();
        return std::make_unique<SimpleIntegratedFlux>(flux);
    }
}
