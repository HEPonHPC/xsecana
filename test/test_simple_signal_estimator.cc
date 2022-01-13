#include "XSecAna/SimpleSignalEstimator.h"
#include "test_utils.h"
#include <iostream>

#include "TFile.h"

using namespace xsec;

int main(int argc, char ** argv) {
    bool verbose = false;
    if (argc > 1 && std::strcmp(argv[1], "-v") == 0) verbose = true;
    bool pass = true;
    bool test;

    auto data = test::utils::get_simple_data();
    auto expected_signal = test::utils::get_simple_signal();

    SimpleSignalEstimator signal_estimator(test::utils::get_simple_background());

    double tol = 1e-14;
    pass &= TEST_HIST("signal",
                      signal_estimator.Eval(data).get(),
                      expected_signal,
                      tol,
                      verbose);

    std::string test_file_name = test::utils::test_dir() + "test_simple_signal_estimator.root";
    auto output = new TFile(test_file_name.c_str(), "recreate");
    signal_estimator.SaveTo(output, "signal_estimator");
    output->Close();
    delete output;

    TFile * input = TFile::Open(test_file_name.c_str());
    auto loaded = IMeasurement::LoadFrom(SimpleSignalEstimator::LoadFrom,
                                         input,
                                         "signal_estimator").release();

    pass &= TEST_HIST("loadfrom",
                      loaded->Eval(data).get(),
                      expected_signal,
                      tol,
                      verbose);

    return !pass;
}
