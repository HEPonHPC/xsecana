#include "XSecAna/Hist.h"
#include "XSecAna/SimpleSignalEstimator.h"
#include "XSecAna/SimpleEfficiency.h"
#include "XSecAna/SimpleFlux.h"
#include "XSecAna/IUnfold.h"
#include "XSecAna/CrossSection.h"
#include "test_utils.h"


#include <iostream>

#include "TFile.h"

using namespace xsec;

const static double ntargets = 1e4;

/////////////////////////////////////////////////////////
// make a cross section object that evaluates out to
//  the input array when folded
IMeasurement * make_simple_xsec(const TH1 * val) {
    auto ones = root::ToROOTLike(val, Array::Ones(val->GetNbinsX() + 2));

    auto signal = (TH1 *) test::utils::get_simple_signal()->Clone();
    signal->Divide(val);
    auto efficiency = new SimpleEfficiency(signal,
                                           ones);
    auto flux = new SimpleFlux(ones);
    auto signal_estimator = new SimpleSignalEstimator(test::utils::get_simple_background());
    auto unfold = new IdentityUnfolder(ones->GetNbinsX() + 2);
    auto ret = new EigenCrossSectionEstimator(efficiency,
                                              signal_estimator,
                                              flux,
                                              unfold,
                                              ntargets);
    return ret;
}

/////////////////////////////////////////////////////////
std::vector<IMeasurement *>
make_simple_xsec_multiverse(const TH1 * hnominal, int nuniverses) {
    std::vector<IMeasurement *> xsec_universes(nuniverses);
    auto hist_universes = test::utils::make_simple_hist_multiverse(hnominal, nuniverses);

    for (auto i = 0; i < nuniverses; i++) {
        xsec_universes[i] = make_simple_xsec(hist_universes[i]);
    }
    return xsec_universes;
}


int main(int argc, char ** argv) {
    bool verbose = false;
    if (argc > 1 && std::strcmp(argv[1], "-v") == 0) verbose = true;
    bool pass = true;
    bool test;

    auto nbins = 12;
    Array bins = Array::Zero(nbins + 1);
    for (auto i = 0u; i < nbins + 1; i++) {
        bins(i) = 2 * i;
    }

    auto ones = test::utils::make_constant_hist(bins, 1);

    auto bkgd = test::utils::make_constant_hist(bins, 2);

    auto data = test::utils::make_constant_hist(bins, 4);

    auto flux_hist = test::utils::make_constant_hist(bins, 5);

    auto eff_num = test::utils::make_constant_hist(bins, 1);
    auto eff_den = test::utils::make_constant_hist(bins, 4);
    auto eff = (TH1 *) eff_num->Clone();
    eff->Divide(eff_num, eff_den, 1, 1, "B");


    auto efficiency = new SimpleEfficiency(eff_num, eff_den);         // = 1/4
    auto flux = new SimpleFlux(flux_hist);                            // = 5
    auto flux_integrated = new SimpleIntegratedFlux(flux_hist);       // = 70
    auto signal_estimator = new SimpleSignalEstimator(bkgd);          // = 2
    auto unfold = new IdentityUnfolder(bkgd->GetNbinsX() + 2); // = 1

    auto xsec_differential = new EigenDifferentialCrossSectionEstimator(efficiency,
                                                                        signal_estimator,
                                                                        flux_integrated,
                                                                        unfold,
                                                                        ntargets);
    auto xsec = new EigenCrossSectionEstimator(efficiency,
                                               signal_estimator,
                                               flux,
                                               unfold,
                                               ntargets);

    auto expected_signal = (TH1 * ) data->Clone();
    expected_signal->Add(bkgd, -1);

    auto expected_xsec = (TH1 *) data->Clone();
    expected_xsec->Add(bkgd, -1);
    expected_xsec->Divide(flux_hist);
    expected_xsec->Divide(eff);

    TEST_HIST("xsec",
              xsec->Eval(data),
              expected_xsec,
              0,
              verbose);

    auto bin_width = test::utils::make_constant_hist_like(expected_xsec, 1);
    for (auto i = 1; i <= bin_width->GetNbinsX(); i++) {
        bin_width->SetBinContent(i, bin_width->GetBinWidth(i));
        bin_width->SetBinError(i, 0);
    }
    bin_width->SetBinError(0, 0);
    bin_width->SetBinError(bin_width->GetNbinsX()+1, 0);


    auto expected_xsec_differential = (TH1 *) data->Clone();
    expected_xsec_differential->Add(bkgd, -1);
    expected_xsec_differential->Divide(eff);
    expected_xsec_differential->Divide(bin_width);

    double flux_integral_error;
    auto flux_integral = flux_hist->IntegralAndError(0, flux_hist->GetNbinsX() + 1,
                                                     flux_integral_error);
    auto integrated_flux_hist = test::utils::make_constant_hist_like(data, flux_integral);

    for (auto i = 0u; i <= integrated_flux_hist->GetNbinsX() + 1; i++) {
        integrated_flux_hist->SetBinError(i, flux_integral_error);
    }
    integrated_flux_hist->GetBinContent(0);
    expected_xsec_differential->Divide(integrated_flux_hist);
    auto result_xsec_differential = xsec_differential->Eval(data);
    TEST_HIST("xsec_differential",
              result_xsec_differential,
              expected_xsec_differential,
              1e-11,
              verbose);

    std::string test_file_name = test::utils::test_dir() + "test_simple_xsec.root";
    auto output = new TFile(test_file_name.c_str(), "recreate");
    xsec->SaveTo(output, "xsec");
    xsec_differential->SaveTo(output, "xsec_differential");
    output->Close();
    delete output;

    auto input = TFile::Open(test_file_name.c_str());
    auto loaded_xsec =
            EigenCrossSectionEstimator::LoadFrom(SimpleEfficiency::LoadFrom,
                                                 SimpleSignalEstimator::LoadFrom,
                                                 SimpleFlux::LoadFrom,
                                                 IdentityUnfolder::LoadFrom,
                                                 input,
                                                 "xsec");
    auto loaded_differential_xsec =
            EigenDifferentialCrossSectionEstimator::LoadFrom(SimpleEfficiency::LoadFrom,
                                                             SimpleSignalEstimator::LoadFrom,
                                                             SimpleIntegratedFlux::LoadFrom,
                                                             IdentityUnfolder::LoadFrom,
                                                             input,
                                                             "xsec_differential");
    input->Close();
    delete input;

    TEST_HIST("loaded xsec",
              loaded_xsec->Eval(data),
              xsec->Eval(data),
              0,
              verbose);

    TEST_HIST("loaded differential xsec",
              loaded_differential_xsec->Eval(data),
              xsec_differential->Eval(data),
              0,
              verbose);

    auto simple_data = test::utils::get_simple_data();
    auto simple_ones = (TH1 *) simple_data->Clone();
    simple_ones->Divide(simple_data);
    assert((root::MapContentsToEigen(simple_ones) -
            root::MapContentsToEigen(make_simple_xsec(simple_ones)->Eval(test::utils::get_simple_data()))
           ).isZero(0));

    return !pass;

}
