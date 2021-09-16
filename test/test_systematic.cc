#include <iostream>
#include <stdio.h>

#include "XSecAna/Hist.h"
#include "XSecAna/Systematic.h"
#include "XSecAna/SimpleEfficiency.h"
#include "XSecAna/SimpleFlux.h"
#include "XSecAna/SimpleSignalEstimator.h"
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
  
  auto nominal = new Hist<Scalar, Cols>(test::utils::get_simple_nominal_hist<Scalar, Cols>());
  auto up      = new Hist<Scalar, Cols>(test::utils::get_simple_up_hist<Scalar, Cols>());
  auto down    = new Hist<Scalar, Cols>(test::utils::get_simple_down_hist<Scalar, Cols>());

  typedef Hist<Scalar, Cols> histtype;
  ForEachFunction<histtype, histtype> subtract = [&nominal](histtype * h) {
        return new histtype(*h - *nominal);
  };

  ForEachFunction<histtype, histtype> divide = [&nominal](const histtype * h) {
        return new histtype((*h - *nominal)/ *nominal);
  };


  // two sided sytematic construction
  Systematic<Hist<Scalar, Cols> > syst_2("syst", up, down);
  TEST_HISTS_SAME("2-sided construction (up)",
                  *(syst_2.Up()),
                  *(up),
                  0);
  TEST_HISTS_SAME("2-sided construction (down)",
                  *(syst_2.Down()),
                  *(down),
                  0);

  auto syst_diff_2 = syst_2.ForEach(subtract);
  TEST_HISTS_SAME("two sided subtraction (up)",
                  *(syst_diff_2.Up()),
                  (*up - *nominal),
                  0);
  TEST_HISTS_SAME("two sided subtraction (down)",
                  *(syst_diff_2.Down()),
                  (*down - *nominal),
                  0);

  auto syst_div_2 = syst_2.ForEach(divide);
  TEST_HISTS_SAME("two sided division (up)",
                  *(syst_div_2.Up()),
                  ((*up - *nominal) / *nominal),
                  0);
  TEST_HISTS_SAME("two sided division (down)",
                  *(syst_div_2.Down()),
                  ((*down - *nominal) / *nominal),
                  0);

  // one sided systematic construction
  Systematic<Hist<Scalar, Cols> > syst_1("syst", up);
  TEST_HISTS_SAME("construction",
                  *(syst_1.Up()),
                  *(up),
                  0);


  auto syst_diff_1 = syst_1.ForEach(subtract);
  TEST_HISTS_SAME("one sided subtraction (up)",
                  *(syst_diff_1.Up()),
                  (*up - *nominal),
                  0);

  // one sided systematic division via invoke
  auto syst_div_1 = syst_1.ForEach(divide);
  TEST_HISTS_SAME("one sided division",
                  *(syst_div_1.Up()),
                  ((*up - *nominal) / *nominal),
                  0);

  TFile * output = new TFile(test_file_name.c_str(), "update");
  TDirectory * to = output->mkdir(dir.c_str());
  syst_2.SaveTo(to, "syst_2");
  syst_1.SaveTo(to, "syst_1");
  output->Close();
  delete output;

  TFile * input = TFile::Open(test_file_name.c_str());
  auto loaded_2 = Systematic<histtype>::LoadFrom(histtype::LoadFrom,
                                                 input->GetDirectory(dir.c_str()),
                                                 "syst_2");
  auto loaded_1 = Systematic<histtype>::LoadFrom(histtype::LoadFrom,
                                                 input->GetDirectory(dir.c_str()),
                                                 "syst_1");
  input->Close();
  delete input;

  TEST_HISTS_SAME("load 2-sided (up)",
                  *(loaded_2->Up()),
                  *(up),
                  0);
  TEST_HISTS_SAME("load 2-sided (down)",
                  *(loaded_2->Down()),
                  *(down),
                  0);

  TEST_HISTS_SAME("load 1-sided (up)",
                  *(loaded_1->Up()),
                  *(up),
                  0);

  // test runtime exceptions
  try {
      syst_1.Down();
      pass &= false;
  }
  catch( exceptions::SystematicTypeError & e) {
    pass &= true;
  }

  try {
      MultiverseShift(syst_1, *nominal, 1);
      pass &= false;
  }
  catch( exceptions::SystematicTypeError & e) {
    pass &= true;
  }

  try {
      MultiverseShift(syst_2, *nominal, 1);
      pass &= false;
  }
  catch( exceptions::SystematicTypeError & e) {
    pass &= true;
  }

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
  std::vector<Hist<Scalar, Cols>*> universes = test::utils::make_simple_hist_multiverse(nominal, nuniverses);
  
  Systematic syst("test_mv", universes);

  auto plus_1sigma = MultiverseShift(syst, nominal, 1);
  auto minus_1sigma = MultiverseShift(syst, nominal, -1);

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
  TEST_HISTS_SAME("minus 1 sigma",
                  minus_1sigma,
                  *(universes[p1_idx]),
                  0);

  // Test ForEach
  // No need to get too fancy here. If it compiles assume it passes
  std::vector<test::utils::Ratio<Scalar, Cols>*>vratio_mv;
  for(auto i = 0u; i < syst.GetShifts().size(); i++) {
    vratio_mv.push_back(new test::utils::Ratio(*syst.GetShifts()[i], nominal));
  }
  
  Systematic ratio_mv("ratio_mv", vratio_mv);
  ForEachFunction<Hist<Scalar, Cols>, test::utils::Ratio<Scalar, Cols>> eval = [](test::utils::Ratio<Scalar, Cols> * ratio) {
      return new Hist<Scalar, Cols>(ratio->Eval());
  };
  auto ratio_mv_hists = ratio_mv.ForEach(eval);

  // save everything for later inspection
  TFile * output = new TFile(test_file_name.c_str(), "update");
  syst.SaveTo(output->mkdir(dir.c_str()), "test_mv");
  nominal.SaveTo(output->GetDirectory(dir.c_str()), "nominal");

  ratio_mv_hists.SaveTo(output->GetDirectory(dir.c_str()), "ratio_mv");
  plus_1sigma.SaveTo(output->GetDirectory(dir.c_str()), "plus_1sigma");
  minus_1sigma.SaveTo(output->GetDirectory(dir.c_str()), "minus_1sigma");

  output->Close();
  delete output;

  // serialization closure
  TFile * input = TFile::Open(test_file_name.c_str());
  auto loaded = Systematic<histtype>::LoadFrom(histtype::LoadFrom,
                                               input->GetDirectory(dir.c_str()),
                                               "test_mv");
  input->Close();
  delete input;
  TEST_MULTIVERSE("loadfrom/saveto", syst, *loaded, 0);

  // test runtime exceptions
  try {
      syst.Up();
      pass &= false;
  }
  catch( exceptions::SystematicTypeError & e) {
    pass &= true;
  }

  // test runtime exceptions
  try {
      syst.Down();
      pass &= false;
  }
  catch( exceptions::SystematicTypeError & e) {
    pass &= true;
  }

    return pass;
}



int main(int argc, char ** argv)
{
  bool verbose = false;
  if(argc > 1 && std::strcmp(argv[1], "-v") == 0)  verbose = true;

  bool pass = true;
  bool test;

  std::remove(test_file_name.c_str());

  pass &= run_tests<double, 10>(verbose, "double_10");
  pass &= run_tests<double, Eigen::Dynamic>(verbose, "double_dynamic");
  pass &= run_tests<float, 10>(verbose, "float_10");
  pass &= run_tests<float, Eigen::Dynamic>(verbose, "float_dynamic");


  pass &= run_tests_mv<double, 10>(verbose, "mv_double_10");
  pass &= run_tests_mv<double, Eigen::Dynamic>(verbose, "mv_double_dynamic");
  pass &= run_tests_mv<float, 10>(verbose, "mv_float_10");
  pass &= run_tests_mv<float, Eigen::Dynamic>(verbose, "mv_float_dynamic");


  auto nominal = test::utils::get_simple_nominal_hist<double, 10>();
  auto up = new test::utils::Ratio<double, 10>(nominal + 1, nominal);
  auto dw = new test::utils::Ratio<double, 10>(nominal - 1, nominal);

  int nuniverses = 50;
  double maxy =  1;
  double miny = -1;
  std::vector<test::utils::Ratio<double, 10>* > universes(nuniverses);
  for(auto i = 0; i < nuniverses; i++) {
    universes[i] = new test::utils::Ratio(nominal + (-1 + (maxy - miny) / (nuniverses-1) * i),
                                          nominal);
  }
  
  // test of polymorphism of objects in the systematics container
  typedef Hist<double, 10> histtype;
  auto eff = new SimpleEfficiency(nominal + 1, nominal);
  Systematic<IEfficiency<histtype>> syst_eff("eff",
                                               eff);
  ForEachFunction<histtype, IEfficiency<histtype>> eval_efficiency = [](IEfficiency<histtype> * eff) {
      return new histtype(eff->Eval());
  };
  Systematic<histtype> eff_res = syst_eff.ForEach(eval_efficiency);

  TEST_HISTS_SAME("polymorphism IEfficiency",
                  *(eff_res.Up()),
                  ((nominal + 1) / nominal),
                  0);

  auto signal = new SimpleSignalEstimator(nominal);
  auto data = nominal + 1;
  Systematic<ISignalEstimator<histtype>> syst_signal("signal",
                                                     signal);
  ForEachFunction<histtype,
                  ISignalEstimator<histtype>>
          eval_signal =
          [&data](ISignalEstimator<histtype> * sig) {
              return new histtype(sig->Eval(data));
          };
  Systematic<histtype> signal_res = syst_signal.ForEach(eval_signal);
  TEST_HISTS_SAME("polymorphism ISignalEstimator",
                  *(signal_res.Up()),
                  (data - nominal),
                  0);

  return !pass;
}
