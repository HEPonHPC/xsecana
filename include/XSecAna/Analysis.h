#pragma once

#include<string>
#include <map>

#include "XSecAna/Systematic.h"
#include "XSecAna/Hist.h"
#include "XSecAna/IMeasurement.h"

// root includes
#include "TDirectory.h"

namespace xsec {
    template<class HistType>
    inline
    void SaveAnalysis(xsec::IMeasurement<HistType> * nominal,
                      std::vector<xsec::Systematic<xsec::IMeasurement<HistType>>> systematics,
                      HistType & data,
                      TDirectory * dir,
                      const std::string & subdir) {
        TDirectory * tmp = gDirectory;

        dir = dir->mkdir(subdir.c_str()); // switch to subdir
        dir->cd();
        TObjString("Analysis").Write("type");

        nominal->SaveTo(dir, "Nominal");
        data.SaveTo(dir, "Data");

        TObjString(std::to_string(systematics.size()).c_str()).Write("NSystematics");
        auto syst_dir = dir->mkdir("Systematics");
        for (auto i = 0u; i < systematics.size(); i++) {
            systematics[i].SaveTo(syst_dir, std::to_string(i));
        }
        tmp->cd();
    }


    template<class HistType>
    inline
    std::tuple<xsec::IMeasurement<HistType>*,
              std::vector<xsec::Systematic<xsec::IMeasurement<HistType>>>,
              HistType>
    LoadAnalysis(xsec::type::LoadFunction<xsec::IMeasurement<HistType>> load,
                 TDirectory * dir,
                 const std::string & subdir) {
                TDirectory * tmp = gDirectory;
        dir = dir->GetDirectory(subdir.c_str());

        auto ptag = (TObjString *) dir->Get("type");
        assert(ptag->GetString() == "Analysis" && "Type does not match Analysis");
        delete ptag;

        auto data = *HistType::LoadFrom(dir, "Data");
        auto nominal = load(dir, "Nominal").release();

        auto syst_dir = dir->GetDirectory("Systematics");
        auto nsysts = std::atoi(((TObjString *) dir->Get("NSystematics"))->GetString().Data());
        std::vector<Systematic<IMeasurement<HistType>>> systematics(nsysts);
        for (auto i = 0u; i < nsysts; i++) {
            systematics[i] = *Systematic<IMeasurement<HistType>>::LoadFrom(load,
                                                                           syst_dir,
                                                                           std::to_string(i));
        }

        tmp->cd();
        return {nominal, systematics, data};
    }
}
