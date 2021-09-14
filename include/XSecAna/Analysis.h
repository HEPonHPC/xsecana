#pragma once

#include<string>
#include <map>

#include "XSecAna/Systematic.h"
#include "XSecAna/Hist.h"
#include "XSecAna/IUncertaintyPropagator.h"
#include "XSecAna/IMeasurement.h"

// root includes
#include "TDirectory.h"

namespace xsec {
    
    template<class HistType = HistXd>
    class Analysis {

    public:
        typedef IUncertaintyPropagator<HistType,
                                       IMeasurement<HistType>,
                                       const HistType>
                UncertaintyPropagatorType;

        Analysis() = default;

        Analysis(IMeasurement<HistType>* nominal_measurement,
                 std::map<std::string, Systematic<IMeasurement<HistType>>> shifted_measurement,
                 HistType data,
                 UncertaintyPropagatorType * propagator=0)
                : fNominalMeasure(nominal_measurement),
                  fShiftedMeasure(shifted_measurement),
                  fData(data),
                  fUncertaintyPropagator(propagator) {}

        Analysis(std::string name, Systematic<IMeasurement<HistType>> shifted_measurement)
                : fShiftedMeasure({name, shifted_measurement}) {}

        explicit Analysis(HistType data)
                : fData(data) {}


        /// \brief Forward to UncertaintyPropagator
        HistType AbsoluteUncertainty(std::string syst_name);

        /// \brief Forward to UncertaintyPropagator
        HistType FractionalUncertainty(std::string syst_name);

        /// \brief Forward to UncertaintyPropagator
        std::pair<HistType, HistType> TotalAbsoluteUncertainty();

        std::pair<HistType, HistType> TotalFractionalUncertainty();

        /// \brief Return an folded cross section result for the input systematic
        const Systematic<HistType> & Result(std::string syst_name);

        /// \brief Return the nominal folded cross section result
        const HistType & Result();

        ~Analysis();

        void SaveTo(TDirectory * dir, const std::string & subdir) const;

        static std::unique_ptr<Analysis> LoadFrom(xsec::type::LoadFunction<IMeasurement<HistType> > load,
                                                  TDirectory * dir,
                                                  const std::string & subdir);

    protected:
        // nominal is special
        IMeasurement<HistType>* fNominalMeasure;

        std::map<std::string, Systematic<IMeasurement<HistType>>> fShiftedMeasure;

        const HistType fData;

        // cache the nominal results
        HistType * fNominalResult = 0;

        // cache the shifted results
        std::map<std::string, Systematic<HistType> > fShiftedResult;

        //
        UncertaintyPropagatorType * fUncertaintyPropagator;
    };

    ///////////////////////////////////////////////////////////////////////
    template<class HistType>
    HistType
    Analysis<HistType>::
    AbsoluteUncertainty(std::string syst_name) {
        return fUncertaintyPropagator->AbsoluteUncertainty(fNominalMeasure,
                                                           fShiftedMeasure.at(syst_name),
                                                           fData);
    }

    ///////////////////////////////////////////////////////////////////////
    template<class HistType>
    HistType
    Analysis<HistType>::
    FractionalUncertainty(std::string syst_name) {
        return fUncertaintyPropagator->FractionalUncertainty(fNominalMeasure,
                                                             fShiftedMeasure.at(syst_name),
                                                             fData);
    }

    ///////////////////////////////////////////////////////////////////////
    template<class HistType>
    std::pair<HistType, HistType>
    Analysis<HistType>::
    TotalAbsoluteUncertainty() {
        return fUncertaintyPropagator->TotalAbsoluteUncertainty(fNominalMeasure,
                                                                fShiftedMeasure,
                                                                fData);
    }

    ///////////////////////////////////////////////////////////////////////
    template<class HistType>
    std::pair<HistType, HistType>
    Analysis<HistType>::
    TotalFractionalUncertainty() {
        return fUncertaintyPropagator->TotalFractionalUncertainty(fNominalMeasure,
                                                                  fShiftedMeasure,
                                                                  fData);
    }

    ///////////////////////////////////////////////////////////////////////
    template<class HistType>
    const Systematic<HistType> &
    Analysis<HistType>::
    Result(std::string syst_name) {
        ForEachFunction<HistType, IMeasurement<HistType>> eval = [this](IMeasurement<HistType> * m) {
            return new HistType(m->Eval(this->fData));
        };
        if (fShiftedResult.find(syst_name) == fShiftedResult.end()) {
            fShiftedResult[syst_name] =
                    fShiftedMeasure.at(syst_name).ForEach(eval);

        }
        return fShiftedResult.at(syst_name);
    }

    ///////////////////////////////////////////////////////////////////////
    template<class HistType>
    const HistType &
    Analysis<HistType>::
    Result() {
        if (!fNominalResult) {
            fNominalResult = new HistType(fNominalMeasure->Eval(fData));
        }
        return *fNominalResult;
    }

    ///////////////////////////////////////////////////////////////////////
    template<class HistType>
    void
    Analysis<HistType>::
    SaveTo(TDirectory * dir, const std::string & subdir) const {
        TDirectory * tmp = gDirectory;

        dir = dir->mkdir(subdir.c_str()); // switch to subdir
        dir->cd();
        TObjString("Analysis").Write("type");

        fNominalMeasure->SaveTo(dir, "fNominalMeasure");
        fData.SaveTo(dir, "fData");

        auto syst_dir = dir->mkdir("fShiftedMeasure");
        for (auto shifted: fShiftedMeasure) {
            shifted.second.SaveTo(syst_dir, shifted.first);
        }

        tmp->cd();
    }

    ///////////////////////////////////////////////////////////////////////
    template<class HistType>
    std::unique_ptr<Analysis<HistType>>
    Analysis<HistType>::
    LoadFrom(xsec::type::LoadFunction<xsec::IMeasurement<HistType>> load,
             TDirectory * dir,
             const std::string & subdir) {
        TDirectory * tmp = gDirectory;
        dir = dir->GetDirectory(subdir.c_str());

        auto ptag = (TObjString *) dir->Get("type");
        assert(ptag->GetString() == "Analysis" && "Type does not match Analysis");
        delete ptag;

        auto data = *HistType::LoadFrom(dir, "fData");

        auto nominal_measurement = load(dir, "fNominalMeasure").release();


        std::map<std::string, Systematic<IMeasurement<HistType>>> shifted_measurement;
        auto syst_dir = dir->GetDirectory("fShiftedMeasure");
        for (auto syst_name: *syst_dir->GetListOfKeys()) {
            shifted_measurement[syst_name->GetName()] = *Systematic<IMeasurement<HistType>>::LoadFrom(load,
                                                                                                      syst_dir,
                                                                                                      syst_name->GetName()).release();
        }

        tmp->cd();
        return std::make_unique<Analysis<HistType> >
                (nominal_measurement,
                 shifted_measurement,
                 data);

    }

    ///////////////////////////////////////////////////////////////////////
    template<class HistType>
    Analysis<HistType>::
    ~Analysis() {
        delete fNominalMeasure;
        delete fNominalResult;
    }

}
