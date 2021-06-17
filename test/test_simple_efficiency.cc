#include "XSecAna/Hist.h"
#include "XSecAna/SimpleEfficiency.h"
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
  
  Hist<double, 10> num(Eigen::Array<double, 1, 10>::Ones(),
		       Eigen::Array<double, 1, 11>::LinSpaced(11, 0, 10));
  Hist<double, 10> den(Eigen::Array<double, 1, 10>::Ones() / 2,
		       Eigen::Array<double, 1, 11>::LinSpaced(11, 0, 10));

  auto expected = Eigen::Array<double, 1, 10>::Ones() * 2;

  SimpleEfficiency eff(num, den);
  
  const Hist<double, 10> * cache_hit1 = &eff.ToHist();
  const Hist<double, 10> * cache_hit2 = &eff.ToHist();

  pass &= cache_hit1 == cache_hit2;

  TEST_ARRAY("ratio calculation", eff.ToHist().Contents(), expected, 0);

  TFile * output = new TFile("test_simple_efficiency.root", "recreate");
  eff.SaveTo(output, "simple_efficiency");
  output->Close();
  delete output;

  TFile * input = TFile::Open("test_simple_efficiency.root");
  auto loaded = SimpleEfficiency<Hist<double, 10> >::LoadFrom(input, "simple_efficiency");
  input->Close();
  delete input;

  TEST_HIST("saveto/loadfrom numerator", loaded->GetNumerator(), num.Contents(), num.Edges(), 0);
  TEST_HIST("saveto/loadfrom denominator", loaded->GetDenominator(), den.Contents(), den.Edges(), 0);
  TEST_HIST("saveto/loadfrom ratio", loaded->ToHist(), eff.ToHist().Contents(), eff.ToHist().Edges(), 0);
    
  if(pass) std::cout << "Success!" << std::endl;
}
