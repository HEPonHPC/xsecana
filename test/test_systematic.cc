#include <iostream>
#include <stdio.h>

#include "XSecAna/Systematic.h"
#include "XSecAna/SimpleEfficiency.h"
#include "XSecAna/SimpleFlux.h"
#include "XSecAna/SimpleSignalEstimator.h"
#include "test_utils.h"

#include <Eigen/Dense>
#include "TFile.h"

using namespace xsec;

std::string test_file_name = test::utils::test_dir() + "test_systematic.root";


bool run_tests(bool verbose, std::string dir) {

    bool pass = true;
    bool test;

    auto nominal = test::utils::get_simple_nominal_hist();
    auto up = test::utils::get_simple_up_hist();
    auto down = test::utils::get_simple_down_hist();

    ForEachFunction<TH1, TH1> subtract = [&nominal](const TH1 * h) {
        auto ret = (TH1 *) h->Clone();
        ret->Add(nominal, -1);
        return ret;
    };

    ForEachFunction<TH1, TH1> divide = [&nominal](const TH1 * h) {
        auto ret = (TH1 *) h->Clone();
        ret->Add(nominal, -1);
        ret->Divide(nominal);
        return ret;
    };


    // two sided sytematic construction
    Systematic<TH1> syst_2("syst", up, down);
    pass &= TEST_HIST("2-sided construction (up)",
                      syst_2.Up(),
                      up,
                      0,
                      verbose);
    pass &= TEST_HIST("2-sided construction (down)",
                      syst_2.Down(),
                      down,
                      0,
                      verbose);

    auto syst_diff_2 = syst_2.ForEach(subtract);
    auto _sub_up = (TH1 *) up->Clone();
    _sub_up->Add(nominal, -1);
    auto _sub_dw = (TH1 *) down->Clone();
    _sub_dw->Add(nominal, -1);
    
    pass &= TEST_HIST("two sided subtraction (up)",
                      syst_diff_2.Up(),
                      _sub_up,
                      0,
                      verbose);
    pass &= TEST_HIST("two sided subtraction (down)",
                      syst_diff_2.Down(),
                      _sub_dw,
                      0,
                      verbose);

    auto syst_div_2 = syst_2.ForEach(divide);
    auto _div_up = (TH1 *) up->Clone();
    _div_up->Add(nominal, -1);
    _div_up->Divide(nominal);

    auto _div_dw = (TH1 *) down->Clone();
    _div_dw->Add(nominal, -1);
    _div_dw->Divide(nominal);

    pass &= TEST_HIST("two sided division (up)",
                      syst_div_2.Up(),
                      _div_up,
                      0,
                      verbose);
    pass &= TEST_HIST("two sided division (down)",
                      syst_div_2.Down(),
                      _div_dw,
                      0,
                      verbose);

    // one sided systematic construction
    Systematic<TH1> syst_1("syst", up);
    pass &= TEST_HIST("construction",
                      syst_1.Up(),
                      up,
                      0,
                      verbose);


    auto syst_diff_1 = syst_1.ForEach(subtract);
    pass &= TEST_HIST("one sided subtraction (up)",
                      syst_diff_1.Up(),
                      _sub_up,
                      0,
                      verbose);

    // one sided systematic division via invoke
    auto syst_div_1 = syst_1.ForEach(divide);
    pass &= TEST_HIST("one sided division",
                      syst_div_1.Up(),
                      _div_up,
                      0,
                      verbose);
    auto expected_covariance_matrix = new TH2D("", "",
                                               nominal->GetNbinsX(), 0, nominal->GetNbinsX(),
                                               nominal->GetNbinsX(), 0, nominal->GetNbinsX());
    for(auto i = 0u; i < nominal->GetNbinsX()+2; i++) {
        for(auto j = 0u; j < nominal->GetNbinsX()+2; j++) {
            auto c_i = nominal->GetBinContent(i) - up->GetBinContent(i);
            auto c_j = nominal->GetBinContent(j) - up->GetBinContent(j);
            expected_covariance_matrix->SetBinContent(i,j, c_i * c_j);
        }
    }

    Array cov = syst_1.CovarianceMatrix(nominal).array().reshaped();
    pass &= TEST_HIST("covariance matrix",
                      root::ToROOTLike(expected_covariance_matrix, cov),
                      expected_covariance_matrix,
                      0,
                      verbose);

    auto output = new TFile(test_file_name.c_str(), "update");
    TDirectory * to = output->mkdir(dir.c_str());
    syst_2.SaveTo(to, "syst_2");
    syst_1.SaveTo(to, "syst_1");
    output->Close();
    delete output;

    TFile * input = TFile::Open(test_file_name.c_str());
    auto loaded_2 = Systematic<TH1>::LoadFrom(root::LoadTH1,
                                              input->GetDirectory(dir.c_str()),
                                              "syst_2");
    auto loaded_1 = Systematic<TH1>::LoadFrom(root::LoadTH1,
                                              input->GetDirectory(dir.c_str()),
                                              "syst_1");
    input->Close();
    delete input;

    pass &= TEST_HIST("load 2-sided (up)",
                      loaded_2->Up(),
                      up,
                      0,
                      verbose);
    pass &= TEST_HIST("load 2-sided (down)",
                      loaded_2->Down(),
                      down,
                      0,
                      verbose);

    pass &= TEST_HIST("load 1-sided (up)",
                      loaded_1->Up(),
                      up,
                      0,
                      verbose);

    // test runtime exceptions
    try {
        syst_1.Down();
        pass &= false;
    }
    catch (exceptions::SystematicTypeError & e) {
        pass &= true;
    }

    try {
        MultiverseShift(syst_1, nominal, 1);
        pass &= false;
    }
    catch (exceptions::SystematicTypeError & e) {
        pass &= true;
    }

    try {
        MultiverseShift(syst_2, nominal, 1);
        pass &= false;
    }
    catch (exceptions::SystematicTypeError & e) {
        pass &= true;
    }

    return pass;
}


template<typename Scalar, int Cols>
bool run_tests_mv(bool verbose, std::string dir) {

    bool pass = true;
    bool test;

    auto nominal = test::utils::get_simple_nominal_hist();
    int nuniverses = 50;
    std::vector<const TH1 *> universes = test::utils::make_simple_hist_multiverse(nominal, nuniverses);

    Systematic syst("test_mv", universes);

    auto plus_1sigma = MultiverseShift(syst, nominal, 1);


    auto minus_1sigma = MultiverseShift(syst, nominal, -1);

    // Systematic::NSigmaShift returns a
    // histogram representing closest universes to 1 sigma away from nominal
    // this is done bin by bin, so the returned histogram holds values that are found
    // in the multiverse
    // With a constant shift up and down, the resulting histogram will be the exact same
    // histogram from one of the universes.
    // For a maximum spread of +/- 1, the shift that gets reconstructed
    // represents the universe corresponding to the 15th and 85th percentile index, for -/+ 1sigma,
    // respectively, in a sorted array of these universes.
    int p1_idx = (0.5 - std::erf(1 / std::sqrt(2)) / 2.0) * (nuniverses - 1) + 1;
    pass &= TEST_HIST("minus 1 sigma",
                      minus_1sigma,
                      universes[p1_idx],
                      0,
                      verbose);


    // save everything for later inspection
    TFile * output = new TFile(test_file_name.c_str(), "update");
    nominal->Write("nominal");
    plus_1sigma->Write("plus_1sigma");
    minus_1sigma->Write("minus_1sigma");
    syst.SaveTo(output->mkdir(dir.c_str()), "test_mv");
    output->Close();
    delete output;

    // serialization closure
    TFile * input = TFile::Open(test_file_name.c_str());
    auto loaded = Systematic<TH1>::LoadFrom(root::LoadTH1,
                                            input->GetDirectory(dir.c_str()),
                                            "test_mv");
    input->Close();
    delete input;
    pass &= TEST_MULTIVERSE("loadfrom/saveto",
                            syst,
                            *loaded,
                            0,
                            verbose);

    // test runtime exceptions
    try {
        syst.Up();
        pass &= false;
    }
    catch (exceptions::SystematicTypeError & e) {
        pass &= true;
    }

    // test runtime exceptions
    try {
        syst.Down();
        pass &= false;
    }
    catch (exceptions::SystematicTypeError & e) {
        pass &= true;
    }

    return pass;
}


int main(int argc, char ** argv) {
    bool verbose = false;
    if (argc > 1 && std::strcmp(argv[1], "-v") == 0) verbose = true;

    bool pass = true;
    bool test;

    std::remove(test_file_name.c_str());

    pass &= run_tests(verbose, "double_dynamic");

    pass &= run_tests_mv<double, Eigen::Dynamic>(verbose, "mv_double_dynamic");

    auto den = test::utils::get_simple_nominal_hist();
    auto num = test::utils::get_simple_down_hist();

    // test of polymorphism of objects in the systematics container
    auto eff = new SimpleEfficiency(num, den);

    Systematic<IMeasurement> syst_eff("eff",
                                      eff);

    Systematic<TH1> eff_res = syst_eff.Eval(num);

    pass &= TEST_HIST("polymorphism IEfficiency",
                      eff_res.Up(),
                      eff->Eval(num),
                      0,
                      verbose);

    auto signal = new SimpleSignalEstimator(num);
    Systematic<IMeasurement> syst_signal("signal",
                                         signal);
    Systematic<TH1> signal_res = syst_signal.Eval(den);
    pass &= TEST_HIST("polymorphism ISignalEstimator",
                      signal_res.Up(),
                      signal->Eval(den),
                      0,
                      verbose);

    return !pass;
}
