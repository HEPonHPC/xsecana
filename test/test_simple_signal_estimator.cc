#include "XSecAna/Hist.h"
#include "XSecAna/SimpleSignalEstimator.h"
#include "XSecAna/test/test_utils.h"
#include <iostream>

#include "TFile.h"

using namespace xsec;

int main(int argc, char ** argv)
{
  bool verbose = false;
  if(argc > 1 && std::strcmp(argv[1], "-v") == 0)  verbose = true;
  bool pass = true;
  bool test;

  Hist<double, 10> bkgd(Eigen::Array<double, 1, 10>::Ones() / 2,
			Eigen::Array<double, 1, 11>::LinSpaced(11, 0, 10));
  Hist<double, 10> data(Eigen::Array<double, 1, 10>::Ones(),
			Eigen::Array<double, 1, 11>::LinSpaced(11, 0, 10));

  SimpleSignalEstimator signal_estimator(bkgd);

  TEST_HIST("signal", signal_estimator.Signal(data), bkgd.Contents(), bkgd.Edges(), 0);

  std::string test_file_name = test::utils::test_dir() + "test_simple_signal_estimator.root";
  TFile * output = new TFile(test_file_name.c_str(), "recreate");
  signal_estimator.SaveTo(output, "signal_estimator");
  output->Close();
  delete output;

  TFile * input = TFile::Open(test_file_name.c_str());
  auto loaded = *SimpleSignalEstimator<Hist<double, 10> >::LoadFrom(input, "signal_estimator");

  TEST_HIST("loadfrom", loaded.Background(data), bkgd.Contents(), bkgd.Edges(), 0);

  return pass;
}
