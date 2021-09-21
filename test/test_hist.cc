#include "XSecAna/Hist.h"
#include "XSecAna/Type.h"
#include "test_utils.h"
#include <boost/histogram.hpp>

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
  auto nbins = 10;
  auto nbins_and_uof = nbins + 2;
  auto nedges = nbins + 1;
  auto nedges_and_uof = nbins_and_uof + 1;
  auto maxx = 10;
  auto minx = 0;
  auto step = (maxx - minx) / nbins;

  auto contents = Hist<Scalar, Cols>::array_and_uof_type::LinSpaced(nbins_and_uof,
                                                                    1,
                                                                    10) - 5.5;
  auto bins     = Hist<Scalar, Cols>::edges_type::LinSpaced(nedges_and_uof,
                                                            minx - step,
                                                            maxx + step);
  auto ones = Hist<Scalar, Cols>::array_and_uof_type::Ones(nbins_and_uof);
  auto errors1 = ones * 3;
  auto errors2 = ones * 4;
  auto exposure1 = 1.;
  auto exposure2 = 3.;

  // total errors relative to hist1
  auto errors_total_rel1 = (errors1.pow(2) +
                            ((exposure1 / exposure2 ) * errors2).pow(2)).sqrt();


  Hist<Scalar, Cols> hist (contents, bins, errors1, exposure1);
  Hist<Scalar, Cols> hist2(contents, bins, errors2, exposure2);

  assert(hist.Contents().size() == nbins);
  assert(hist.Edges().size() == nbins+1);
  assert(hist.Errors().size() == nbins);

  // make sure this is a deep copy
  auto hist_cpy = hist;
  assert((hist.EdgesAndUOF() - hist_cpy.EdgesAndUOF()).isZero(0));
  assert((hist.ContentsAndUOF() - hist_cpy.ContentsAndUOF()).isZero(0));
  assert((hist.ErrorsAndUOF() - hist_cpy.ErrorsAndUOF()).isZero(0));
  assert(hist.Exposure() == hist_cpy.Exposure());
  assert(hist.EdgesAndUOF().data() != hist_cpy.EdgesAndUOF().data());
  assert(hist.ContentsAndUOF().data() != hist_cpy.ContentsAndUOF().data());
  assert(hist.ErrorsAndUOF().data() != hist_cpy.ErrorsAndUOF().data());

  double tol = 0;
  if constexpr(std::is_same<Scalar, double>::value) {
      tol = 1e-14;
    }
  else {
    tol = 1e-6;
  }

  bool pass = true;
  bool test;

  TEST_HIST_AND_EDGES("construction", hist, contents, bins, 0);

  TEST_HIST_AND_EDGES("hist multiply", hist * hist2, contents.pow(2) / 3, bins, tol);
  TEST_ARRAY_SAME("hist multiply error",
                  (hist * hist2).ErrorsAndUOF(),
                  errors_total_rel1,
                  0);

  TEST_HIST_AND_EDGES("hist divide", hist / hist2, ones * 3, bins, tol);
  TEST_ARRAY_SAME("hist divide error",
                  (hist / hist2).ErrorsAndUOF(),
                  errors_total_rel1,
                  0);

  TEST_HIST_AND_EDGES("hist add", hist + hist2, contents * (4. / 3.), bins, tol);
  TEST_ARRAY_SAME("hist add error",
                  (hist + hist2).ErrorsAndUOF(),
                  errors_total_rel1,
                  0);

  TEST_HIST_AND_EDGES("hist subtract", hist - hist2, contents * (2. / 3.), bins, tol);
  TEST_ARRAY_SAME("hist subtract error",
                  (hist - hist2).ErrorsAndUOF(),
                  errors_total_rel1,
                  0);

  // although this isn't technically correct because the expression
  // can be simplified to prevent double counting of errors from hist2,
  // we don't do the simplification when calculating errors.
  // The user needs to be made aware of this behavior so they can
  // do their own simplification.
  auto errors_total_rel1_chained_expression = (errors1.pow(2) +
                                               2*((exposure1 / exposure2 ) * errors2).pow(2)).sqrt();
  TEST_HIST_AND_EDGES("hist chained expression", (hist - hist2) / hist2, ones * (2.), bins, tol);
  TEST_ARRAY_SAME("hist chained expression error",
                  ((hist - hist2) / hist2).ErrorsAndUOF(),
                  errors_total_rel1_chained_expression,
                  1e-6);
  
  TEST_HIST_AND_EDGES("constant multiply", (hist * -1), contents * -1, bins, 0);
  TEST_ARRAY_SAME("constant multiply error",
                  (hist * -1).ErrorsAndUOF(),
                  errors1,
                  0);

  TEST_HIST_AND_EDGES("abs()", hist.abs(), contents.abs(), bins, 0);
  TEST_ARRAY_SAME("abs() error",
                  hist.abs().ErrorsAndUOF(),
                  errors1,
                  0);

  TEST_HIST_AND_EDGES("abs().sqrt()", hist.abs().sqrt(), contents.abs().sqrt(), bins, 1e-6);
  TEST_ARRAY_SAME("abs().sqrt() error",
                  hist.abs().sqrt().ErrorsAndUOF(),
                  errors1 / 2,
                  0);

  TEST_HIST_AND_EDGES("no change", hist, contents, bins, 0);

  hist = hist.abs();
  TEST_HIST_AND_EDGES("modify abs", hist, contents.abs(), bins, 0);
  TEST_ARRAY_SAME("constant multiply error",
                  hist.ErrorsAndUOF(),
                  errors1,
                  0);
  
  TEST_ARRAY_SAME("bin width",
                  hist.BinWidths(),
                  ones(Eigen::seqN(0, ones.size()-2)),
                       0);

  TEST_ARRAY_SAME("contiguous", hist.Contents(), test::utils::is_contiguous(hist.Contents()), 0);

  // saveto/loadfrom
  std::string test_file_name = test::utils::test_dir() + "test_hist.root";
  TFile * output = new TFile(test_file_name.c_str(), "recreate");
  hist.SaveTo(output, "hist");
  output->Close();
  
  TFile * input = TFile::Open(test_file_name.c_str());
  auto loaded = Hist<Scalar, Cols>::LoadFrom(input, "hist");

  TEST_HIST_AND_EDGES("saveto/loadfrom",
                      (*loaded),
                      hist.ContentsAndUOF(),
                      hist.EdgesAndUOF(),
                      0);
  TEST_ARRAY_SAME("loadfrom errors",
                  loaded->ErrorsAndUOF(),
                  hist.ErrorsAndUOF(),
                  0)

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
