#include "XSecAna/Hist.h"
#include "XSecAna/SimpleEfficiency.h"
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
  
  Hist<double, 10> num = test::utils::get_hist_of_ones<double, 10>();
  Hist<double, 10> den = test::utils::get_hist_of_ones<double, 10>() / 2;

  SimpleEfficiency eff(num, den);

  TEST_HISTS_SAME("ratio calculation",
                  eff.Eval(),
                  (num / den),
                  0);

  std::string test_file_name = test::utils::test_dir() + "test_simple_efficiency.root";
  TFile * output = new TFile(test_file_name.c_str(), "recreate");
  eff.SaveTo(output, "simple_efficiency");
  output->Close();
  delete output;

  TFile * input = TFile::Open(test_file_name.c_str());
  auto loaded = IEfficiency<Hist<double, 10> >::LoadFrom(SimpleEfficiency<Hist<double, 10>>::LoadFrom,
                                                         input,
                                                         "simple_efficiency").release();
  input->Close();
  delete input;

  TEST_HISTS_SAME("saveto/loadfrom numerator",
                  ((SimpleEfficiency<Hist<double, 10> >*)loaded)->GetNumerator(),
                  num,
                  0);
  TEST_HISTS_SAME("saveto/loadfrom denominator",
                  ((SimpleEfficiency<Hist<double, 10> > *) loaded)->GetDenominator(),
                  den,
                  0);
  TEST_HISTS_SAME("saveto/loadfrom ratio",
                  loaded->Eval(),
                  eff.Eval(),
                  0);

  return !pass;
}
