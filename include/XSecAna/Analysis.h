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

    template<class MeasurementType,
             class HistType = HistXd>
    class Analysis {

    public:
        typedef IUncertaintyPropagator<HistType,
                                       MeasurementType,
                                       const HistType>
                UncertaintyPropagatorType;

        Analysis() = default;

        Analysis(MeasurementType nominal_measurement,
                 std::map<std::string, Systematic<MeasurementType> > shifted_measurement,
                 HistType data,
                 UncertaintyPropagatorType * propagator)
                : fNominalMeasure(nominal_measurement),
                  fShiftedMeasure(shifted_measurement),
                  fData(data),
                  fUncertaintyPropagator(propagator) {}

        Analysis(std::string name, Systematic<MeasurementType> shifted_measurement)
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

        static std::unique_ptr<Analysis> LoadFrom(TDirectory * dir, const std::string & subdir);

    protected:
        // nominal is special
        MeasurementType fNominalMeasure;

        std::map<std::string, Systematic<MeasurementType> > fShiftedMeasure;

        const HistType fData;

        // cache the nominal results
        HistType * fNominalResult = 0;

        // cache the shifted results
        std::map<std::string, Systematic<HistType> > fShiftedResult;

        //
        UncertaintyPropagatorType * fUncertaintyPropagator;
    };

    ///////////////////////////////////////////////////////////////////////
    template<class MeasurementType,
             class HistType>
    HistType
    Analysis<MeasurementType, HistType>::
    AbsoluteUncertainty(std::string syst_name) {
        return fUncertaintyPropagator->AbsoluteUncertainty(fNominalMeasure,
                                                           fShiftedMeasure.at(syst_name),
                                                           fData);
    }

    ///////////////////////////////////////////////////////////////////////
    template<class MeasurementType,
             class HistType>
    HistType
    Analysis<MeasurementType, HistType>::
    FractionalUncertainty(std::string syst_name) {
        return fUncertaintyPropagator->FractionalUncertainty(fNominalMeasure,
                                                             fShiftedMeasure.at(syst_name),
                                                             fData);
    }

    ///////////////////////////////////////////////////////////////////////
    template<class MeasurementType,
             class HistType>
    std::pair<HistType, HistType>
    Analysis<MeasurementType, HistType>::
    TotalAbsoluteUncertainty() {
        return fUncertaintyPropagator->TotalAbsoluteUncertainty(fNominalMeasure,
                                                                fShiftedMeasure,
                                                                fData);
    }

    ///////////////////////////////////////////////////////////////////////
    template<class MeasurementType,
             class HistType>
    std::pair<HistType, HistType>
    Analysis<MeasurementType, HistType>::
    TotalFractionalUncertainty() {
        return fUncertaintyPropagator->TotalFractionalUncertainty(fNominalMeasure,
                                                                  fShiftedMeasure,
                                                                  fData);
    }

    ///////////////////////////////////////////////////////////////////////
    template<class MeasurementType,
             class HistType>
    const Systematic<HistType> &
    Analysis<MeasurementType, HistType>::
    Result(std::string syst_name) {
        if (fShiftedResult.find(syst_name) == fShiftedResult.end()) {
            fShiftedResult[syst_name] =
                    fShiftedMeasure.at(syst_name).Invoke(&std::remove_pointer<MeasurementType>::type::Eval,
                                                         fData);

        }
        return fShiftedResult.at(syst_name);
    }

    ///////////////////////////////////////////////////////////////////////
    template<class MeasurementType,
             class HistType>
    const HistType &
    Analysis<MeasurementType, HistType>::
    Result() {
        if (!fNominalResult) {
            fNominalResult = new HistType(fNominalMeasure->Eval(fData));
        }
        return *fNominalResult;
    }

    ///////////////////////////////////////////////////////////////////////
    template<class MeasurementType,
             class HistType>
    void
    Analysis<MeasurementType, HistType>::
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
    template<class MeasurementType,
             class HistType>
    std::unique_ptr<Analysis<MeasurementType, HistType> >
    Analysis<MeasurementType, HistType>::
    LoadFrom(TDirectory * dir, const std::string & subdir) {
        TDirectory * tmp = gDirectory;
        dir = dir->GetDirectory(subdir.c_str());

        auto ptag = (TObjString *) dir->Get("type");
        assert(ptag->GetString() == "Analysis" && "Type does not match Analysis");
        delete ptag;

        auto data = *HistType::LoadFrom(dir, "fData");

        auto nominal_measurement = MeasurementType::LoadFrom(dir, "fNominalMeasure").release();


        std::map<std::string, Systematic<MeasurementType> > shifted_measurement;
        auto syst_dir = dir->GetDirectory("fShiftedMeasure");
        for (auto syst_name: *syst_dir->GetListOfKeys()) {
            shifted_measurement[syst_name->GetName()] = *Systematic<MeasurementType>::LoadFrom(syst_dir,
                                                                                               syst_name->GetName()).release();
        }

        tmp->cd();
        return std::make_unique<Analysis<MeasurementType, HistType> >
                (nominal_measurement,
                 shifted_measurement,
                 data);

    }

    ///////////////////////////////////////////////////////////////////////
    template<class MeasurementType,
             class HistType>
    Analysis<MeasurementType, HistType>::
    ~Analysis() {
        delete fNominalResult;
    }

}
