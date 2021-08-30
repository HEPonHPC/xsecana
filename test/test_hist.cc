#include "XSecAna/Hist.h"
#include "XSecAna/Type.h"
#include "test_utils.h"


#include <Eigen/Dense>

#include <iostream>
#include <memory>
#include <iterator>
#include <type_traits>

#include "TFile.h"

using namespace xsec;


template<typename Scalar, int Cols>
bool run_tests(bool verbose)
{
  auto contents = Eigen::Array<Scalar, 1, Cols>::LinSpaced(10, 1, 10) - 5.5;
  auto bins     = Eigen::Array<Scalar, 1, xsec::EdgesSize(Cols)>::LinSpaced(11, 0, 10);

  auto ones = Eigen::Array<Scalar, 1, Cols>::Ones(10);
  Hist<Scalar, Cols> hist (contents, bins, 1); // exposure of 1
  Hist<Scalar, Cols> hist2(contents, bins, 3); // exposure of 2
  
  double tol = 0;
  if constexpr(std::is_same<Scalar, double>::value) {
      tol = 1e-14;
    }
  else {
    tol = 1e-6;
  }

  bool pass = true;
  bool test;

  TEST_HIST("construction", hist, contents, bins, 0);

  TEST_HIST("hist multiply", hist * hist2, contents.pow(2) / 3, bins, tol);

  TEST_HIST("hist divide", hist / hist2, ones * 3, bins, tol);

  TEST_HIST("hist add", hist + hist2, contents * (4. / 3.), bins, tol);

  TEST_HIST("hist subtract", hist - hist2, contents * (2. / 3.), bins, tol);

  TEST_HIST("hist chained expression", (hist - hist2) / hist2, ones * (2.), bins, tol);
  
  TEST_HIST("constant multiply", (hist*-1), contents*-1, bins, 0);

  TEST_HIST("abs()", hist.abs(), contents.abs(), bins, 0);

  TEST_HIST("abs().sqrt()", hist.abs().sqrt(), contents.abs().sqrt(), bins, 1e-6);

  TEST_HIST("no change", hist, contents, bins, 0);

  hist = hist.abs();
  TEST_HIST("modify abs", hist, contents.abs(), bins, 0);
  
  TEST_ARRAY("bin width", hist.BinWidths(), ones, 0);

  TEST_ARRAY("contiguous", hist.Contents(), test::utils::is_contiguous(hist.Contents()), 0);

  // saveto/loadfrom
  std::string test_file_name = test::utils::test_dir() + "test_hist.root";
  TFile * output = new TFile(test_file_name.c_str(), "recreate");
  hist.SaveTo(output, "hist");
  output->Close();
  
  TFile * input = TFile::Open(test_file_name.c_str());
  auto loaded = Hist<Scalar, Cols>::LoadFrom(input, "hist");

  TEST_HIST("saveto/loadfrom", (*loaded), hist.Contents(), hist.Edges(), 0);

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

  pass &= !type::IsHist<double>();
  pass &=  type::IsHist<Hist<double, 1> >();
  
  return !pass;
}
