#include <iostream>
#include <stdio.h>

#include "XSecAna/Hist.h"
#include "XSecAna/Systematic.h"

#include <Eigen/Dense>
#include "TFile.h"

using namespace xsec;

template<class Scalar, 
	 int Cols>
struct Ratio {
  Ratio(Hist<Scalar, Cols> num,
	Hist<Scalar, Cols> den)
    : numerator(num), denominator(den)
  {}

  Hist<Scalar, Cols> numerator;
  Hist<Scalar, Cols> denominator;
  
  Hist<Scalar, Cols> ToHist() const
  {
    return numerator / denominator;
  }
  
};

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

#define TEST_MULTIVERSE(test_name, mv1, mv2)	\
  test = true;								\
  for(auto imv = 0u; imv < (mv1).GetUniverses().size(); imv++) {		\
    test &= (mv1).GetUniverses()[imv] == (mv2).GetUniverses()[imv];	\
  }									\
  if(!test || verbose) {						\
    std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << (test? ": PASSED" : ": FAILED") << std::endl; \
    for(auto imv = 0u; imv < (mv1).GetUniverses().size(); imv++) {	\
      std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << "[" << imv << "]\t" << (mv1).GetUniverses()[imv].Contents() << std::endl; \
      std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << "[" << imv << "]\t" << (mv2).GetUniverses()[imv].Contents() << std::endl; \
    }									\
    pass = false;							\
  }									\


template<typename Scalar, int Cols>
bool run_tests(bool verbose, std::string dir)
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

  TFile * output = new TFile("test_systematic.root", "update");
  TDirectory * to = output->mkdir(dir.c_str());
  syst_2.SaveTo(to, "syst_2");
  syst_1.SaveTo(to, "syst_1");
  output->Close();
  delete output;

  TFile * input = TFile::Open("test_systematic.root");
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
  auto bins = Eigen::Array<Scalar, 1, EdgesSize(Cols)>::LinSpaced(11, 0, 10);
  
  auto vnominal = Eigen::Array<Scalar, 1, Cols>::LinSpaced(10, -5, 5);
  Hist<Scalar, Cols> nominal(vnominal, bins);
  int nuniverses = 50;
  Scalar maxy =  1;
  Scalar miny = -1;
  std::vector<Hist<Scalar, Cols> > universes(nuniverses);
  for(auto i = 0; i < nuniverses; i++) {
    universes[i] = nominal + (-1 + (maxy - miny) / (nuniverses-1) * i);
  }
  
  MultiverseSystematic syst("test_mv", universes);

  auto plus_1sigma = syst.NSigmaShift(1, nominal);
  auto minus_1sigma = syst.NSigmaShift(-1, nominal);

  // MultiverseSystematic::NSigmaShift returns a 
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

  // Test the ability of MultiverseSystematic<T> to call T::ToHist
  // No need to get too fancy here. If it compiles assume it passes
  std::vector<Ratio<Scalar, Cols> >vratio_mv;
  for(auto i = 0u; i < syst.GetUniverses().size(); i++) {
    vratio_mv.push_back(Ratio(syst.GetUniverses()[i], nominal));
  }
  
  MultiverseSystematic ratio_mv("ratio_mv", vratio_mv);
  MultiverseSystematic ratio_mv_hists = ratio_mv.Invoke(&Ratio<Scalar, Cols>::ToHist);

  auto ratio_plus_1sigma = ratio_mv.NSigmaShift(1, Ratio(nominal, nominal));
  auto ratio_minus_1sigma = ratio_mv.NSigmaShift(-11, Ratio(nominal, nominal));

  
  // save everything for later inspection
  TFile * output = new TFile("test_systematic.root", "update");
  syst.SaveTo(output->mkdir(dir.c_str()), "test_mv");
  ratio_mv_hists.SaveTo(output->GetDirectory(dir.c_str()), "ratio_mv");
  nominal.SaveTo(output->GetDirectory(dir.c_str()), "nominal");
  plus_1sigma.SaveTo(output->GetDirectory(dir.c_str()), "plus_1sigma");
  minus_1sigma.SaveTo(output->GetDirectory(dir.c_str()), "minus_1sigma");
  ratio_plus_1sigma.SaveTo(output->GetDirectory(dir.c_str()),  "ratio_plus_1sigma");
  ratio_minus_1sigma.SaveTo(output->GetDirectory(dir.c_str()), "ratio_minus_1sigma");


  output->Close();
  delete output;

  // serialization closure
  TFile * input = TFile::Open("test_systematic.root");
  auto loaded = MultiverseSystematic<histtype>::LoadFrom(input->GetDirectory(dir.c_str()), "test_mv");
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

  std::remove("test_systematic.root");

  pass &= run_tests<double, 10>(verbose, "double_10");
  pass &= run_tests<double, Eigen::Dynamic>(verbose, "double_dynamic");
  pass &= run_tests<float, 10>(verbose, "float_10");
  pass &= run_tests<float, Eigen::Dynamic>(verbose, "float_dynamic");

  pass &= run_tests_mv<double, 10>(verbose, "mv_double_10");
  pass &= run_tests_mv<double, Eigen::Dynamic>(verbose, "mv_double_dynamic");
  pass &= run_tests_mv<float, 10>(verbose, "mv_float_10");
  pass &= run_tests_mv<float, Eigen::Dynamic>(verbose, "mv_float_dynamic");
  if(pass) std::cout << "Success!" << std::endl;
}
