#include "XSecAna/Hist.h"
#include "XSecAna/SimpleSignalEstimator.h"
#include "XSecAna/SimpleEfficiency.h"
#include "XSecAna/SimpleFlux.h"
#include "XSecAna/IUnfold.h"
#include "XSecAna/ICrossSection.h"
#include "XSecAna/SimpleQuadSum.h"
#include <iostream>

#include "TFile.h"

using namespace xsec;

#define TEST_ARRAY(test_name, arr1, arr2, precision)			\
  test = ((arr1) - (arr2)).isZero(precision);				\
  if(!test || verbose) {						\
    std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << (test? ": PASSED" : ": FAILED") << std::endl; \
    std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << "\t" << (arr1) << std::endl; \
    std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << "\t" << (arr2) << std::endl; \
    pass &= test;							\
  }									

#define TEST_HIST(test_name,HIST, target_contents, target_edges, precision) \
  test = (HIST.Contents() - target_contents).isZero(precision);			\
  if(!test || verbose) {						\
    std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << (test? ": PASSED" : ": FAILED") << std::endl; \
    std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << "\t" << HIST.Contents() << std::endl; \
    std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << "\t" << target_contents << std::endl; \
    pass &= test;							\
  }									\
  test = (HIST.Edges() - target_edges).isZero(precision);			\
  if(!test || verbose) {						\
    std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << (test? ": PASSED" : ": FAILED") << std::endl; \
    std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << "\t" << HIST.Edges() << std::endl; \
    std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << "\t" << target_edges << std::endl; \
    pass &= test;							\
  }								

typedef Hist<double, 10> histtype;
typedef ICrossSection<SimpleSignalEstimator<histtype>,
		      DummyUnfold<double, 10>,
		      SimpleEfficiency<histtype>,
		      SimpleFlux<histtype> > SimpleCrossSection;

// make a cross section object that evaluates out to 
//  - the input array when folded
//  - 2 times the input array when unfolded
SimpleCrossSection make_simple_xsec(Hist<double, 10> val)
{
  Hist<double, 10> ones(Eigen::Array<double, 1, 10>::Ones(),
			Eigen::Array<double, 1, 11>::LinSpaced(11, 0, 10));

  auto efficiency = new SimpleEfficiency<histtype>(ones, val);           // = 1 / val
  auto flux = new SimpleFlux(ones);                                      // = 1
  auto signal_estimator = new SimpleSignalEstimator(ones);               // = 1
  auto unfold = new DummyUnfold<double, 10>(ones.Contents().size(), 2);  // = 2
  return SimpleCrossSection(efficiency,
			    signal_estimator,
			    flux,
			    unfold);
}

int main(int argc, char ** argv)
{
  bool verbose = false;
  if(argc > 1 && std::strcmp(argv[1], "-v") == 0)  verbose = true;
  bool pass = true;
  bool test;
  
  SimpleQuadSum<SimpleCrossSection, histtype> prop;
  
  Hist<double, 10> ones(Eigen::Array<double, 1, 10>::Ones(),
			Eigen::Array<double, 1, 11>::LinSpaced(11, 0, 10));
  Hist<double, 10> data(Eigen::Array<double, 1, 10>::Ones() * 2,
			Eigen::Array<double, 1, 11>::LinSpaced(11, 0, 20));
  Hist<double, 10> hnominal = ones + 0.2;
  double ntargets = 1e4;

  SimpleCrossSection nominal_xsec = make_simple_xsec(hnominal);
  SimpleCrossSection up   = make_simple_xsec(hnominal + std::sqrt(2) );  
  SimpleCrossSection down = make_simple_xsec(hnominal - std::sqrt(2) / 2);

  /* 
     multiverse example 
  */
  int nuniverses = 50;
  double maxy =  1;
  double miny = -1;
  std::vector<SimpleCrossSection> xsec_universes(nuniverses);
  std::vector<Hist<double, 10> > hist_universes(nuniverses);
  for(auto i = 0; i < nuniverses; i++) {
    xsec_universes[i] = make_simple_xsec(hnominal + (-1 + (maxy - miny) / (nuniverses-1) * i));
    hist_universes[i] = hnominal + (-1 + (maxy - miny) / (nuniverses-1) * i);
  }
  
  Systematic<SimpleCrossSection> syst_mv("mv_xsec", xsec_universes);
  Systematic<Hist<double, 10> > syst_mv_hist("mv_hist", hist_universes);

  auto abs_uncert_mv = prop.AbsoluteUncertaintyXSec(data,
						    nominal_xsec,
						    syst_mv,
						    ntargets);
  
  // index of universe representing minus 1 sigma shift
  int m1_idx = (0.5 - std::erf(1 / std::sqrt(2)) / 2.0) * (nuniverses-1) + 1;

  TEST_ARRAY("minus 1 sigma multiverse",
  	     abs_uncert_mv.Contents(),
  	     (hist_universes[m1_idx]-hnominal).Contents().abs(), 1e-14);


  /*
    Examples using 1 and 2 sided shifts
    Test each function of the propogator
  */
  Systematic<SimpleCrossSection> syst_1sided("1sided", up);
  Systematic<SimpleCrossSection> syst_2sided("2sided", up, down);

  // AbsoluteUncertaintyXSec
  auto abs_uncert_1sided = prop.AbsoluteUncertaintyXSec(data,
							nominal_xsec,
							syst_1sided,
							ntargets);
  auto abs_uncert_2sided = prop.AbsoluteUncertaintyXSec(data,
							nominal_xsec,
							syst_2sided,
							ntargets);

  TEST_ARRAY("abs_uncert",
	     abs_uncert_1sided.Contents(),
	     (ones * std::sqrt(2)).Contents(),
	     0);
  TEST_ARRAY("max shift",
	     abs_uncert_1sided.Contents(),
	     abs_uncert_2sided.Contents(),
	     0);

  // AbsoluteUncertaintyUnfoldedXSec
  auto abs_uncert_1sided_unfolded = prop.AbsoluteUncertaintyUnfoldedXSec(data,
									 nominal_xsec,
									 syst_1sided,
									 ntargets);

  TEST_ARRAY("abs_uncert unfolded",
	     abs_uncert_1sided_unfolded.Contents(),
	     (ones * (2 * std::sqrt(2))).Contents(),
	     0);

  // FractionalUncertaintyXSec
  auto frac_uncert_1sided = prop.FractionalUncertaintyXSec(data, 
							   nominal_xsec, 
							   syst_1sided,
							   ntargets);
  TEST_ARRAY("fractional uncert",
	     frac_uncert_1sided.Contents(),
	     (((hnominal + std::sqrt(2)) - hnominal) / hnominal).Contents(), 
	     0);

  // FractionalUncertaintyUnfoldedXSec
  auto frac_uncert_1sided_unfolded = prop.FractionalUncertaintyUnfoldedXSec(data, 
									    nominal_xsec, 
									    syst_1sided,
									    ntargets);
  TEST_ARRAY("fractional uncert",
	     frac_uncert_1sided_unfolded.Contents(),
	     (((hnominal + std::sqrt(2)) - hnominal) / hnominal).Contents(), 
	     0);
  

  // TotalAbsoluteUncertaintyXSec
  std::map<std::string, Systematic<SimpleCrossSection> > systs = {
    {"1sided", syst_1sided},
    {"2sided", syst_2sided},
    {"mv", syst_mv},
  };

  auto total_abs_uncert = prop.TotalAbsoluteUncertaintyXSec(data,
							    nominal_xsec,
							    systs,
							    ntargets);
  TEST_ARRAY("total absolute uncert",
	     total_abs_uncert.first.Contents(),
	     ((hist_universes[m1_idx]-hnominal).Contents().abs().pow(2) + ones.Contents() * 4).sqrt(),
	     1e-14);

  // TotalAbsoluteUncertaintyUnfoldedXSec
  auto total_abs_uncert_unfolded = prop.TotalAbsoluteUncertaintyUnfoldedXSec(data,
									     nominal_xsec,
									     systs,
									     ntargets);
  TEST_ARRAY("total absolute uncert unfolded",
	     total_abs_uncert_unfolded.first.Contents(),
	     ((hist_universes[m1_idx]-hnominal).Contents().abs().pow(2) + ones.Contents() * 4).sqrt() * 2,
	     1e-14);

  // TotalFractionalUncertaintyXSec
  auto total_frac_uncert = prop.TotalFractionalUncertaintyXSec(data,
							       nominal_xsec,
							       systs,
							       ntargets);
  TEST_ARRAY("total frac uncert",
	     total_frac_uncert.first.Contents(),
	     ((hist_universes[m1_idx]-hnominal).Contents().abs().pow(2) + ones.Contents() * 4 ).sqrt() / hnominal.Contents(),
	     1e-14);

  // TotalFractionalUncertaintyUnfoldedXSec
  auto total_frac_uncert_unfolded = prop.TotalFractionalUncertaintyUnfoldedXSec(data,
										nominal_xsec,
										systs,
										ntargets);
  TEST_ARRAY("total frac uncert unfolded",
	     total_frac_uncert_unfolded.first.Contents(),
	     ((hist_universes[m1_idx]-hnominal).Contents().abs().pow(2) + ones.Contents() * 4 ).sqrt() / hnominal.Contents(),
	     1e-14);


	     
  if(pass) std::cout << "Success!" << std::endl;
}
