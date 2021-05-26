#include <iostream>

#include "XSecAna/Hist.h"
#include "XSecAna/Systematic.h"

#include <Eigen/Dense>
#include "TFile.h"

using namespace xsec;
#define TEST_ARRAY(test_name, arr1, arr2, precision)			\
  test = (arr1 - arr2).isZero(precision);				\
  if(!test || verbose) {						\
    std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << (test? ": PASSED" : ": FAILED") << std::endl; \
    std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << "\t" << arr1 << std::endl; \
    std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << "\t" << arr2 << std::endl; \
    pass = false;							\
  }									

#define TEST_SYSTEMATIC(test_name, syst, up, down)			\
  test = *syst.GetShifts().first == (up);				\
  if(!test || verbose) {						\
    std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << (test? ": PASSED" : ": FAILED") << std::endl; \
    std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << "\t" << syst.GetShifts().first->Contents() << std::endl; \
    std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << "\t" << (up).Contents() << std::endl; \
    pass = false;							\
  }									\
  test = *syst.GetShifts().second == (down);				\
  if(!test || verbose) {				\
    std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << (test? ": PASSED" : ": FAILED") << std::endl; \
    std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << "\t" << syst.GetShifts().second->Contents() << std::endl; \
    std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << "\t" << (down).Contents() << std::endl; \
    pass = false;							\
  }									

template<typename Scalar, int Cols>
bool run_tests(bool verbose)
{
  bool pass = true;
  bool test;
  auto bins = Eigen::Array<Scalar, 1, EdgesSize(Cols)>::LinSpaced(11, 0, 10);
  
  auto vnominal = Eigen::Array<Scalar, 1, Cols>::Ones(10) * 5;
  auto vup      = Eigen::Array<Scalar, 1, Cols>::LinSpaced(10, 4, 6);
  auto vdown    = Eigen::Array<Scalar, 1, Cols>::LinSpaced(10, 4, 6).reverse();

  Hist<Scalar, Cols> nominal(vnominal, bins);
  Hist<Scalar, Cols> up     (vup     , bins);
  Hist<Scalar, Cols> down   (vdown   , bins);

  // two sided sytematic construction
  Systematic<Hist<Scalar, Cols> > syst_2("syst", up, down);
  TEST_SYSTEMATIC("construction", syst_2, up, down);

  // two sided systematic subtraction via invoke
  typedef Hist<Scalar, Cols> histtype;
  Systematic<Hist<Scalar, Cols> > syst_diff_2 = 
    syst_2.Invoke(static_cast<histtype(histtype::*)(const histtype&) const>(&histtype::operator-), nominal);
  TEST_SYSTEMATIC("two sided subtraction", syst_diff_2, up - nominal, down - nominal);

  // two sided systematic division via invoke
  Systematic<Hist<Scalar, Cols> > syst_div_2 = 
    syst_diff_2.Invoke(static_cast<histtype(histtype::*)(const histtype&) const>(&histtype::operator/), nominal);
  TEST_SYSTEMATIC("two sided division", syst_div_2, (up - nominal) / nominal, (down - nominal) / nominal);
  
  // one sided systematic construction
  Systematic<Hist<Scalar, Cols> > syst_1("syst", up);
  TEST_SYSTEMATIC("construction", syst_1, up, up);

  // one sided systematic subtraction via invoke
  Systematic<Hist<Scalar, Cols> > syst_diff_1 = 
    syst_1.Invoke(static_cast<histtype(histtype::*)(const histtype&) const>(&histtype::operator-), nominal);
  TEST_SYSTEMATIC("one sided subtraction", syst_diff_1, up - nominal, up - nominal);

  // one sided systematic division via invoke
  Systematic<Hist<Scalar, Cols> > syst_div_1 = 
    syst_diff_1.Invoke(static_cast<histtype(histtype::*)(const histtype&) const>(&histtype::operator/), nominal);
  TEST_SYSTEMATIC("one sided division", syst_div_1, (up - nominal) / nominal, (up - nominal) / nominal);

  TFile * output = new TFile("test_systematic.root", "recreate");
  syst_2.SaveTo(output, "syst_2");
  syst_1.SaveTo(output, "syst_1");
  output->Close();

  TFile * input = TFile::Open("test_systematic.root");
  auto loaded_2 = Systematic<histtype>::LoadFrom(input, "syst_2");
  auto loaded_1 = Systematic<histtype>::LoadFrom(input, "syst_1");

  TEST_SYSTEMATIC("load 2-sided", (*loaded_2), up, down);
  TEST_SYSTEMATIC("load 1-sided", (*loaded_1), up, up  );

  return pass;
}


int main(int argc, char ** argv)
{
  bool verbose = false;
  if(argc > 1 && std::strcmp(argv[1], "-v") == 0)  verbose = true;

  bool pass = true;
  pass &= run_tests<double, 10>(verbose);
  pass &= run_tests<double, Eigen::Dynamic>(verbose);
  pass &= run_tests<float, 10>(verbose);
  pass &= run_tests<float, Eigen::Dynamic>(verbose);
  if(pass) std::cout << "Success!" << std::endl;
}
