#include <iostream>
#include <stdio.h>

#include "XSecAna/Hist.h"
#include "XSecAna/Systematic.h"
#include "test_utils.h"

#include <Eigen/Dense>
#include "TFile.h"

using namespace xsec;

std::string test_file_name = test::utils::test_dir() + "test_systematic.root";

template<typename Scalar, int Cols>
bool run_tests(bool verbose, std::string dir)
{
  bool pass = true;
  bool test;
  
  auto nominal = test::utils::get_simple_nominal_hist<Scalar, Cols>();
  auto up      = test::utils::get_simple_up_hist<Scalar, Cols>();
  auto down    = test::utils::get_simple_down_hist<Scalar, Cols>();

  // two sided sytematic construction
  Systematic<Hist<Scalar, Cols> > syst_2("syst", up, down);
  TEST_SYSTEMATIC("construction", syst_2, up, down);

  // two sided systematic subtraction via invoke
  typedef Hist<Scalar, Cols> histtype;
  auto syst_diff_2 = syst_2.Invoke(static_cast<histtype(histtype::*)(const histtype&) const>(&histtype::operator-), nominal);
  TEST_SYSTEMATIC("two sided subtraction", syst_diff_2, up - nominal, down - nominal);

  // two sided systematic division via invoke
  auto syst_div_2 = syst_diff_2.Invoke(static_cast<histtype(histtype::*)(const histtype&) const>(&histtype::operator/), nominal);
  TEST_SYSTEMATIC("two sided division", syst_div_2, (up - nominal) / nominal, (down - nominal) / nominal);
  
  // one sided systematic construction
  Systematic<Hist<Scalar, Cols> > syst_1("syst", up);
  TEST_SYSTEMATIC("construction", syst_1, up, up);

  // one sided systematic subtraction via invoke
  auto syst_diff_1 = syst_1.Invoke(static_cast<histtype(histtype::*)(const histtype&) const>(&histtype::operator-), nominal);
  TEST_SYSTEMATIC("one sided subtraction", syst_diff_1, up - nominal, up - nominal);

  // one sided systematic division via invoke
  auto syst_div_1 = syst_diff_1.Invoke(static_cast<histtype(histtype::*)(const histtype&) const>(&histtype::operator/), nominal);
  TEST_SYSTEMATIC("one sided division", syst_div_1, (up - nominal) / nominal, (up - nominal) / nominal);

  TFile * output = new TFile(test_file_name.c_str(), "update");
  TDirectory * to = output->mkdir(dir.c_str());
  syst_2.SaveTo(to, "syst_2");
  syst_1.SaveTo(to, "syst_1");
  output->Close();
  delete output;

  TFile * input = TFile::Open(test_file_name.c_str());
  auto loaded_2 = Systematic<histtype>::LoadFrom(input->GetDirectory(dir.c_str()), "syst_2");
  auto loaded_1 = Systematic<histtype>::LoadFrom(input->GetDirectory(dir.c_str()), "syst_1");
  input->Close();
  delete input;

  TEST_SYSTEMATIC("load 2-sided", (*loaded_2), up, down);
  TEST_SYSTEMATIC("load 1-sided", (*loaded_1), up, up  );

  return pass;
}


template<typename Scalar, int Cols>
bool run_tests_mv(bool verbose, std::string dir)
{
  typedef Hist<Scalar, Cols> histtype;

  bool pass = true;
  bool test;
  
  auto nominal = test::utils::get_simple_nominal_hist<Scalar, Cols>();
  int nuniverses = 50;
  std::vector<Hist<Scalar, Cols> > universes = test::utils::make_simple_hist_multiverse(nominal, nuniverses);
  
  Systematic syst("test_mv", universes);

  auto plus_1sigma = syst.NSigmaShift(1, nominal);
  auto minus_1sigma = syst.NSigmaShift(-1, nominal);

  // Systematic::NSigmaShift returns a 
  // histogram representing closest universes to 1 sigma away from nominal
  // this is done bin by bin, so the returned histogram holds values that are found
  // in the multiverse
  // With a constant shift up and down, the resulting histogram will be the exact same 
  // histogram from one of the universes.
  // For a maximum spread of +/- 1, the shift that gets reconstructed
  // represents the universe corresponding to the 15th and 85th percentile index, for -/+ 1sigma,
  // respectively, in a sorted array of these universes.
  // This 
  int p1_idx = (0.5 - std::erf(1 / std::sqrt(2)) / 2.0) * (nuniverses-1) + 1;
  TEST_ARRAY("minus 1 sigma", minus_1sigma.Contents(), universes[p1_idx].Contents(), 0);

  // Test the ability of Systematic<T> to call T::ToHist
  // No need to get too fancy here. If it compiles assume it passes
  std::vector<test::utils::Ratio<Scalar, Cols> >vratio_mv;
  for(auto i = 0u; i < syst.GetShifts().size(); i++) {
    vratio_mv.push_back(test::utils::Ratio(syst.GetShifts()[i], nominal));
  }
  
  Systematic ratio_mv("ratio_mv", vratio_mv);
  auto ratio_mv_hists = ratio_mv.Invoke(&test::utils::Ratio<Scalar, Cols>::ToHist);

  auto ratio_plus_1sigma = ratio_mv.NSigmaShift(1, test::utils::Ratio(nominal, nominal));
  auto ratio_minus_1sigma = ratio_mv.NSigmaShift(-1, test::utils::Ratio(nominal, nominal));
  
  // save everything for later inspection
  TFile * output = new TFile(test_file_name.c_str(), "update");
  syst.SaveTo(output->mkdir(dir.c_str()), "test_mv");
  nominal.SaveTo(output->GetDirectory(dir.c_str()), "nominal");

  ratio_mv_hists.SaveTo(output->GetDirectory(dir.c_str()), "ratio_mv");
  plus_1sigma.SaveTo(output->GetDirectory(dir.c_str()), "plus_1sigma");
  minus_1sigma.SaveTo(output->GetDirectory(dir.c_str()), "minus_1sigma");
  ratio_plus_1sigma.SaveTo(output->GetDirectory(dir.c_str()),  "ratio_plus_1sigma");
  ratio_minus_1sigma.SaveTo(output->GetDirectory(dir.c_str()), "ratio_minus_1sigma");

  output->Close();
  delete output;

  // serialization closure
  TFile * input = TFile::Open(test_file_name.c_str());
  auto loaded = Systematic<histtype>::LoadFrom(input->GetDirectory(dir.c_str()), "test_mv");
  input->Close();
  delete input;
  TEST_MULTIVERSE("loadfrom/saveto", syst, *loaded);

  return pass;
}



int main(int argc, char ** argv)
{
  bool verbose = false;
  if(argc > 1 && std::strcmp(argv[1], "-v") == 0)  verbose = true;

  bool pass = true;

  std::remove(test_file_name.c_str());

  pass &= run_tests<double, 10>(verbose, "double_10");
  pass &= run_tests<double, Eigen::Dynamic>(verbose, "double_dynamic");
  pass &= run_tests<float, 10>(verbose, "float_10");
  pass &= run_tests<float, Eigen::Dynamic>(verbose, "float_dynamic");

  pass &= run_tests_mv<double, 10>(verbose, "mv_double_10");
  pass &= run_tests_mv<double, Eigen::Dynamic>(verbose, "mv_double_dynamic");
  pass &= run_tests_mv<float, 10>(verbose, "mv_float_10");
  pass &= run_tests_mv<float, Eigen::Dynamic>(verbose, "mv_float_dynamic");

  

  // test of runtime
  auto nominal = test::utils::get_simple_nominal_hist<double, 10>();
  test::utils::Ratio<double, 10> up(nominal + 1, nominal);
  test::utils::Ratio<double, 10> dw(nominal - 1, nominal);

  int nuniverses = 50;
  double maxy =  1;
  double miny = -1;
  std::vector<test::utils::Ratio<double, 10> > universes(nuniverses);
  for(auto i = 0; i < nuniverses; i++) {
    universes[i] = test::utils::Ratio(nominal + (-1 + (maxy - miny) / (nuniverses-1) * i),
				      nominal);
  }
  
  std::vector<Systematic<test::utils::Ratio<double, 10> > * > systs;
  systs.push_back(new Systematic<test::utils::Ratio<double, 10> >("syst", up, dw));
  systs.push_back(new Systematic<test::utils::Ratio<double, 10> >("mv_syst", universes));

  try {
    systs[0]->NSigmaShift(0, test::utils::Ratio(nominal, nominal));
    pass &= false;
  }
  catch( exceptions::SystematicTypeError & e) {
    if(verbose) std::cout << e.what() << std::endl;
    pass &= true;
  }
  try {
    systs[0]->Up();
    pass &= true;
  }
  catch( exceptions::SystematicTypeError & e) {
    if(verbose) std::cout << e.what() << std::endl;
    pass &= false;
  }


  try {
    systs[1]->NSigmaShift(0, test::utils::Ratio(nominal, nominal));
    pass &= true;
  }
  catch( exceptions::SystematicTypeError & e) {
    if(verbose) std::cout << e.what() << std::endl;
    pass &= false;
  }
  try {
    systs[1]->Up();
    pass &= false;
  }
  catch( exceptions::SystematicTypeError & e) {
    if(verbose) std::cout << e.what() << std::endl;
    pass &= true;
  }

  return pass;
}
