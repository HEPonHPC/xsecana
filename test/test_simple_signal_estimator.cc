#include "XSecAna/Hist.h"
#include "XSecAna/SimpleSignalEstimator.h"

#include <iostream>

#include "TFile.h"

using namespace xsec;

#define TEST_ARRAY(test_name, arr1, arr2, precision)			\
  test = (arr1 - arr2).isZero(precision);					\
  if(!test || verbose) {						\
    std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << (test? ": PASSED" : ": FAILED") << std::endl; \
    std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << "\t" << arr1 << std::endl; \
    std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << "\t" << arr2 << std::endl; \
    pass = false;							\
  }									

#define TEST_HIST(test_name,HIST, target_contents, target_edges, precision) \
  test = (HIST.Contents() - target_contents).isZero(precision);			\
  if(!test || verbose) {						\
    std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << (test? ": PASSED" : ": FAILED") << std::endl; \
    std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << "\t" << HIST.Contents() << std::endl; \
    std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << "\t" << target_contents << std::endl; \
    pass = false;							\
  }									\
  test = (HIST.Edges() - target_edges).isZero(precision);			\
  if(!test || verbose) {						\
    std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << (test? ": PASSED" : ": FAILED") << std::endl; \
    std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << "\t" << HIST.Edges() << std::endl; \
    std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << "\t" << target_edges << std::endl; \
    pass = false;							\
  }								

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
  
  TFile * output = new TFile("test_simple_signal_estimator.root", "recreate");
  signal_estimator.SaveTo(output, "signal_estimator");
  output->Close();
  delete output;

  TFile * input = TFile::Open("test_simple_signal_estimator.root");
  auto loaded = *SimpleSignalEstimator<Hist<double, 10> >::LoadFrom(input, "signal_estimator");

  TEST_HIST("loadfrom", loaded.Background(data), bkgd.Contents(), bkgd.Edges(), 0);


  if(pass) std::cout << "Success!" << std::endl;
}
