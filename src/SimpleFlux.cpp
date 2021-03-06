//
// Created by Derek Doyle on 10/9/21.
//

#include "XSecAna/SimpleFlux.h"
#include "TDirectory.h"

namespace xsec {
    SimpleFlux::
    SimpleFlux(const TH1 * flux)
            : fFlux(flux) {}

    SimpleIntegratedFlux::
    SimpleIntegratedFlux(const TH1 * flux)
            : fFlux(flux),
              fN(root::MapContentsToEigen(fFlux).sum()),
              fdN(std::sqrt(fN)) {}

    //////////////////////////////////////////////////////////
    void
    SimpleFlux::
    SaveTo(TDirectory * dir, const std::string & subdir) const {
        TDirectory * tmp = gDirectory;
        dir = dir->mkdir(subdir.c_str());
        dir->cd();

        TObjString("SimpleFlux").Write("type");
        fFlux->Write("fFlux");

        tmp->cd();
    }

    //////////////////////////////////////////////////////////
    std::unique_ptr<IMeasurement>
    SimpleFlux::
    LoadFrom(TDirectory * dir, const std::string & subdir) {
        dir = dir->GetDirectory(subdir.c_str());

        // make sure we're loading the right type
        auto ptag = (TObjString *) dir->Get("type");
        assert(ptag->GetString() == "SimpleFlux" && "Type does not match SimpleFlux");
        delete ptag;

        auto flux = root::LoadTH1(dir, "fFlux").release();
        return std::make_unique<SimpleFlux>(flux);
        flux = 0;
    }

    void
    SimpleIntegratedFlux::
    _eval_impl(const Array & data, const Array & error, ArrayRef result, ArrayRef rerror) const {
        result = Array::Ones(result.size()) * fN;
        rerror = Array::Ones(result.size()) * fdN;
    }

    void
    SimpleFlux::
    _eval_impl(const Array & data, const Array & error, ArrayRef result, ArrayRef rerror) const {
        root::MapToEigen(this->fFlux, result, rerror);
    }

    IMeasurement *
    SimpleFlux::
    IntegratedFlux() const {
        return new SimpleIntegratedFlux(fFlux);
    }


    //////////////////////////////////////////////////////////
    void
    SimpleIntegratedFlux::
    SaveTo(TDirectory * dir, const std::string & subdir) const {
        TDirectory * tmp = gDirectory;
        dir = dir->mkdir(subdir.c_str());
        dir->cd();

        TObjString("SimpleIntegratedFlux").Write("type");
        fFlux->Write("fFlux");

        tmp->cd();
    }

    //////////////////////////////////////////////////////////
    std::unique_ptr<IMeasurement>
    SimpleIntegratedFlux::
    LoadFrom(TDirectory * dir, const std::string & subdir) {
        dir = dir->GetDirectory(subdir.c_str());

        // make sure we're loading the right type
        auto ptag = (TObjString *) dir->Get("type");
        assert(ptag->GetString() == "SimpleIntegratedFlux" && "Type does not match SimpleFlux");
        delete ptag;

        auto flux = root::LoadTH1(dir, "fFlux").release();
        return std::make_unique<SimpleIntegratedFlux>(flux);
        flux = 0;
    }
}

