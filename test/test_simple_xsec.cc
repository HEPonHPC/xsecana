#include "XSecAna/Hist.h"
#include "XSecAna/SimpleSignalEstimator.h"
#include "XSecAna/SimpleEfficiency.h"
#include "XSecAna/SimpleFlux.h"
#include "XSecAna/IUnfold.h"
#include "XSecAna/ICrossSection.h"

#include <iostream>

#include "TFile.h"

using namespace xsec;

#define TEST_ARRAY(test_name, arr1, arr2, precision)			\
  test = (arr1 - arr2).isZero(precision);					\
  if(!test || verbose) {						\
    std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << (test? ": PASSED" : ": FAILED") << std::endl; \
    std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << "\t" << arr1 << std::endl; \
    std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << "\t" << arr2 << std::endl; \
    pass = false;							\
  }									

#define TEST_HIST(test_name,HIST, target_contents, target_edges, precision) \
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



int main(int argc, char ** argv)
{
  bool verbose = false;
  if(argc > 1 && std::strcmp(argv[1], "-v") == 0)  verbose = true;
  bool pass = true;
  bool test;

  Hist<double, 10> bkgd(Eigen::Array<double, 1, 10>::Ones(),
			Eigen::Array<double, 1, 11>::LinSpaced(11, 0, 20));
  Hist<double, 10> data(Eigen::Array<double, 1, 10>::Ones() * 4,
			Eigen::Array<double, 1, 11>::LinSpaced(11, 0, 20));

  Hist<double, 10> flux_hist(Eigen::Array<double, 1, 10>::Ones() * 5,
			     Eigen::Array<double, 1, 11>::LinSpaced(11, 0, 20));

  Hist<double, 10> eff_num(Eigen::Array<double, 1, 10>::Ones() / 4,
			   Eigen::Array<double, 1, 11>::LinSpaced(11, 0, 20));
  Hist<double, 10> eff_den(Eigen::Array<double, 1, 10>::Ones(),
			   Eigen::Array<double, 1, 11>::LinSpaced(11, 0, 20));

  typedef Hist<double, 10> histtype;

  // why is this template deduction failing??
  auto efficiency = new SimpleEfficiency<histtype>(eff_num, eff_den); // = 1/4
  auto flux = new SimpleFlux(flux_hist);                              // = 5
  auto signal_estimator = new SimpleSignalEstimator(bkgd);            // = 3
  auto unfold = new DummyUnfold<double, 10>(bkgd.Contents().size());  // = 1

  ICrossSection xsec(efficiency,
		     signal_estimator,
		     flux,
		     unfold); // = 1 / 1e4

  auto xsec_differential = xsec.ToDifferential();
  
  TEST_ARRAY("xsec", 
	     xsec.CrossSection(data, (double) (12. / 5. ) * 1e4).Contents(),
	     (Eigen::Array<double, 1, 10>::Ones()),
	     0);

  TEST_ARRAY("xsec_differential", 
	     xsec_differential.CrossSection(data, (double) (12. / 5. ) * 1e4).Contents(),
	     (Eigen::Array<double, 1, 10>::Ones() / 2),
	     0);


  TFile * output = new TFile("test_simple_xsec.root", "recreate");
  xsec.SaveTo(output, "xsec");
  output->Close();
  delete output;
  
  

  if(pass) std::cout << "Success!" << std::endl;
}
