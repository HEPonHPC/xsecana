#include "XSecAna/Hist.h"

#include <Eigen/Dense>

#include <iostream>
#include <memory>
#include <iterator>

#include "TFile.h"

using namespace xsec;

#define TEST_HIST(test_name,HIST, target_contents, target_edges, verbose, precision) \
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
  
#define TEST_ARRAY(test_name, arr1, arr2, verbose, precision)			\
  test = (arr1 - arr2).isZero(precision);					\
  if(!test || verbose) {						\
    std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << (test? ": PASSED" : ": FAILED") << std::endl; \
    std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << "\t" << arr1 << std::endl; \
    std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << "\t" << arr2 << std::endl; \
    pass = false;							\
  }									

template<typename ArrayType>
ArrayType is_contiguous(const ArrayType & array)
{ 
  ArrayType a2 = array;
  auto const data = array.data();
  for (auto i = 0; i < array.size(); ++i) {
    a2(i) = *(&data[0] + i);
  }        
  return a2;
}

template<typename Scalar, int Cols>
bool run_tests(bool verbose)
{
  auto contents = Eigen::Array<Scalar, 1, Cols>::LinSpaced(10, 1, 10) - 5;
  auto bins     = Eigen::Array<Scalar, 1, xsec::EdgesSize(Cols)>::LinSpaced(11, 0, 10);

  Hist<Scalar, Cols> hist(contents, bins);

  bool pass = true;
  bool test;


  TEST_HIST("construction", hist, contents, bins, verbose, 0);
  
  TEST_HIST("constant multiply", (hist*-1), contents*-1, bins, verbose, 0);

  TEST_HIST("abs()", hist.abs(), contents.abs(), bins, verbose, 0);

  TEST_HIST("abs().sqrt()", hist.abs().sqrt(), contents.abs().sqrt(), bins, verbose, 1e-6);

  TEST_HIST("no change", hist, contents, bins, verbose, 0);

  hist = hist.abs();
  TEST_HIST("modify abs", hist, contents.abs(), bins, verbose, 0);
  
  auto ones = Eigen::Array<Scalar, 1, Cols>::Ones(10);
  TEST_ARRAY("bin width", hist.BinWidths(), ones, verbose, 0);

  TEST_ARRAY("contiguous", hist.Contents(), is_contiguous(hist.Contents()), verbose, 0);

  // saveto/loadfrom
  TFile * output = new TFile("test_hist.root", "recreate");
  hist.SaveTo(output, "hist");
  output->Close();
  
  TFile * input = TFile::Open("test_hist.root");
  auto loaded = Hist<Scalar, Cols>::LoadFrom(input, "hist");

  TEST_HIST("saveto/loadfrom", (*loaded), hist.Contents(), hist.Edges(), verbose, 0);

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

  // TODO need to test dynamic arrays
  if(pass) std::cout << "Success!" << std::endl;
}
