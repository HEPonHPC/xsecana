#include "XSecAna/Hist.h"

#include <Eigen/Dense>

#include <iostream>
#include <memory>
#include <iterator>

#include "TFile.h"

using namespace xsec;

#define TEST_HIST(test_name,HIST, target_contents, target_edges, verbose)	\
  test = (HIST.Contents() - target_contents).isZero(0);			\
  if(!test || verbose) {						\
    std::cerr << "test_hist " << test_name << (test? ": PASSED" : ": FAILED") << std::endl; \
    std::cerr << "test_hist " << test_name << "\t" << HIST.Contents() << std::endl; \
    std::cerr << "test_hist " << test_name << "\t" << target_contents << std::endl; \
    pass = false;							\
  }									\
  test = (HIST.Edges() - target_edges).isZero(0);			\
  if(!test || verbose) {						\
    std::cerr << "test_hist " << test_name << (test? ": PASSED" : ": FAILED") << std::endl; \
    std::cerr << "test_hist " << test_name << "\t" << HIST.Edges() << std::endl; \
    std::cerr << "test_hist " << test_name << "\t" << target_edges << std::endl; \
    pass = false;							\
  }								
  
#define TEST_ARRAY(test_name, arr1, arr2, verbose)				\
  test = (arr1 - arr2).isZero(0);					\
  if(!test || verbose) {						\
    std::cerr << "test_hist " << test_name << (test? ": PASSED" : ": FAILED") << std::endl; \
    std::cerr << "test_hist " << test_name << "\t" << arr1 << std::endl; \
    std::cerr << "test_hist " << test_name << "\t" << arr2 << std::endl; \
    pass = false;							\
  }									

Eigen::Array<double, 1, 10> is_contiguous(const Eigen::Array<double, 1, 10> & array)
{ 
  Eigen::Array<double, 1, 10> a2 = array;
  auto const data = array.data();
  for (auto i = 0; i < array.size(); ++i) {
    a2(i) = *(&data[0] + i);
  }        
  return a2;
}

int main(int argc, char ** argv)
{
  bool verbose = false;
  if(argc > 1 && std::strcmp(argv[1], "-v") == 0)  verbose = true;

  // TODO need to test dynamic arrays
  auto contents = Eigen::Array<double, 1, 10>::LinSpaced(10, 1, 10) - 5;
  auto bins     = Eigen::Array<double, 1, 11>::LinSpaced(11, 0, 10);

  Hist<double, 10> hist(contents, bins);

  bool pass = true;
  bool test;


  TEST_HIST("construction", hist, contents, bins, verbose);
  
  TEST_HIST("constant multiply", (hist*-1), contents*-1, bins, verbose);

  TEST_HIST("abs()", hist.abs(), contents.abs(), bins, verbose);

  TEST_HIST("abs().sqrt()", hist.abs().sqrt(), contents.abs().sqrt(), bins, verbose);

  TEST_HIST("no change", hist, contents, bins, verbose);

  hist = hist.abs();
  TEST_HIST("modify abs", hist, contents.abs(), bins, verbose);
  
  auto ones = Eigen::Array<double, 1, 10>::Ones();
  TEST_ARRAY("bin width", hist.BinWidths(), ones, verbose);

  TEST_ARRAY("contiguous", hist.Contents(), is_contiguous(hist.Contents()), verbose);

  // saveto/loadfrom
  TFile * output = new TFile("test_hist.root", "recreate");
  hist.SaveTo(output, "hist");
  output->Close();
  
  TFile * input = TFile::Open("test_hist.root");
  auto loaded = Hist<double, 10>::LoadFrom(input, "hist");

  TEST_HIST("saveto/loadfrom", loaded, hist.Contents(), hist.Edges(), verbose);

  if(pass) std::cout << "Success!" << std::endl;
}
