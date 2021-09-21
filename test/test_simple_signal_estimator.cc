#include "XSecAna/Hist.h"
#include "XSecAna/SimpleSignalEstimator.h"
#include "test_utils.h"
#include <iostream>

#include "TFile.h"

using namespace xsec;

int main(int argc, char ** argv)
{
  bool verbose = false;
  if(argc > 1 && std::strcmp(argv[1], "-v") == 0)  verbose = true;
  bool pass = true;
  bool test;

  auto data = test::utils::get_simple_data<double, 10>();
  auto expected_signal = test::utils::get_simple_signal<double, 10>();

  SimpleSignalEstimator signal_estimator(test::utils::get_simple_background<double, 10>());

  TEST_HIST_AND_EDGES("signal",
                      signal_estimator.Signal(data),
                      expected_signal.ContentsAndUOF(),
                      expected_signal.EdgesAndUOF(),
                      0);

  std::string test_file_name = test::utils::test_dir() + "test_simple_signal_estimator.root";
  TFile * output = new TFile(test_file_name.c_str(), "recreate");
  signal_estimator.SaveTo(output, "signal_estimator");
  output->Close();
  delete output;

  TFile * input = TFile::Open(test_file_name.c_str());
  auto loaded = ISignalEstimator<Hist<double, 10> >::LoadFrom(SimpleSignalEstimator<Hist<double, 10>>::LoadFrom,
                                                              input,
                                                              "signal_estimator").release();

  TEST_HISTS_SAME("loadfrom",
                  loaded->Signal(data),
                  expected_signal,
                  0);

  return !pass;
}
