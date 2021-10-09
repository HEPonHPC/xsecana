#pragma once
#include<string>
#include <map>

#include "XSecAna/Systematic.h"
#include "XSecAna/Hist.h"
#include "XSecAna/IMeasurement.h"

// root includes
#include "TDirectory.h"

namespace xsec {
    inline
    void SaveAnalysis(xsec::IMeasurement * nominal,
                      std::vector<xsec::Systematic<xsec::IMeasurement>> systematics,
                      Hist & data,
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

    inline
    std::tuple<xsec::IMeasurement *,
               std::vector<xsec::Systematic<xsec::IMeasurement>>,
               Hist>
    LoadAnalysis(xsec::type::LoadFunction<xsec::IMeasurement> load,
                 TDirectory * dir,
                 const std::string & subdir) {
        TDirectory * tmp = gDirectory;
        dir = dir->GetDirectory(subdir.c_str());

        auto ptag = (TObjString *) dir->Get("type");
        assert(ptag->GetString() == "Analysis" && "Type does not match Analysis");
        delete ptag;

        auto data = *Hist::LoadFrom(dir, "Data");
        auto nominal = load(dir, "Nominal").release();

        auto syst_dir = dir->GetDirectory("Systematics");
        auto nsysts = std::atoi(((TObjString *) dir->Get("NSystematics"))->GetString().Data());
        std::vector<Systematic<IMeasurement>> systematics(nsysts);
        for (auto i = 0u; i < nsysts; i++) {
            systematics[i] = *Systematic<IMeasurement>::LoadFrom(load,
                                                                 syst_dir,
                                                                 std::to_string(i));
        }

        tmp->cd();
        return {nominal, systematics, data};
    }
}
