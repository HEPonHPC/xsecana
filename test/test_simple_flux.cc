#include "XSecAna/SimpleFlux.h"
#include "test_utils.h"
#include <iostream>

#include "TFile.h"

using namespace xsec;

int main(int argc, char ** argv) {
    bool verbose = false;
    if (argc > 1 && std::strcmp(argv[1], "-v") == 0) verbose = true;
    bool pass = true;
    bool test;

    auto flux_hist = test::utils::get_hist_of_ones();

    // arbitrary binning
    int nana_bins = 8;
    Array analysis_binning = Array::LinSpaced(nana_bins+1, 0, nana_bins);
    auto ones_analysis_binning = new TH1D("", "",
                                          nana_bins,
                                          analysis_binning.data());
    ones_analysis_binning->SetContent(Array::Ones(nana_bins+2).eval().data());

    SimpleFlux flux(flux_hist);
    SimpleIntegratedFlux integrated_flux(flux_hist);


    pass &= TEST_HIST("flux.Eval()",
                      flux.Eval(flux_hist),
                      flux_hist,
                      0,
                      verbose);

    auto target_integrated_flux = root::ToROOTLike(ones_analysis_binning,
                                                   Array::Ones(ones_analysis_binning->GetNbinsX() + 2) *
                                                   flux_hist->Integral(),
                                                   Array::Ones(ones_analysis_binning->GetNbinsX() + 2) *
                                                   std::sqrt(flux_hist->Integral()));
    TEST_HIST("integrated_flux.Eval()",
              integrated_flux.Eval(ones_analysis_binning),
              target_integrated_flux,
              0,
              verbose);


    std::string test_file_name = test::utils::test_dir() + "test_simple_flux.root";
    auto output = new TFile(test_file_name.c_str(), "recreate");
    flux.SaveTo(output, "flux");
    integrated_flux.SaveTo(output, "integrated_flux");
    output->Close();
    delete output;

    TFile * input = TFile::Open(test_file_name.c_str());
    auto loaded_flux = IMeasurement::LoadFrom(SimpleFlux::LoadFrom,
                                              input,
                                              "flux");
    auto loaded_integrated_flux = IMeasurement::LoadFrom(SimpleIntegratedFlux::LoadFrom,
                                                         input,
                                                         "integrated_flux");
    input->Close();
    delete input;

    TEST_HIST("loaded_flux",
              loaded_flux->Eval(flux_hist),
              flux_hist,
              0,
              verbose);
    TEST_HIST("loaded_integrated_flux",
              loaded_integrated_flux->Eval(ones_analysis_binning),
              target_integrated_flux,
              0,
              verbose);

    return !pass;
}
