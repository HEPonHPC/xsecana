#include "XSecAna/Hist.h"
#include "XSecAna/SimpleSignalEstimator.h"
#include "XSecAna/SimpleEfficiency.h"
#include "XSecAna/SimpleFlux.h"
#include "XSecAna/IUnfold.h"
#include "XSecAna/ICrossSection.h"
#include "XSecAna/test/test_utils.h"


#include <iostream>

#include "TFile.h"

using namespace xsec;

typedef Hist<double, 10> histtype;
typedef ICrossSection<SimpleSignalEstimator<histtype>,
		      test::utils::DummyUnfold<double, 10>,
		      SimpleEfficiency<histtype>,
		      SimpleFlux<histtype> > SimpleCrossSection;

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
  auto unfold = new test::utils::DummyUnfold<double, 10>(bkgd.Contents().size(), 2);  // = 2

  SimpleCrossSection xsec(efficiency,
			  signal_estimator,
			  flux,
			  unfold); // = 1 / 1e4

  auto xsec_differential = xsec.ToDifferential();
  
  TEST_ARRAY("xsec", 
	     xsec.CrossSection(data, (double) (12. / 5. ) * 1e4).Contents(),
	     (Eigen::Array<double, 1, 10>::Ones()),
	     0);
  TEST_ARRAY("xsec unfolded", 
	     xsec.UnfoldedCrossSection(data, (double) (12. / 5. ) * 1e4).Contents(),
	     (Eigen::Array<double, 1, 10>::Ones()*2),
	     0);

  TEST_ARRAY("xsec_differential", 
	     xsec_differential.CrossSection(data, (double) (12. / 5. ) * 1e4).Contents(),
	     (Eigen::Array<double, 1, 10>::Ones() / 2),
	     0);

  std::string test_file_name = test::utils::test_dir() + "test_simple_xsec.root";
  TFile * output = new TFile(test_file_name.c_str(), "recreate");
  xsec.SaveTo(output, "xsec");
  output->Close();
  delete output;
  
  TFile * input = TFile::Open(test_file_name.c_str());
  auto loaded_xsec = *SimpleCrossSection::LoadFrom(input, "xsec");
  input->Close();
  delete input;

  TEST_ARRAY("loaded xsec", 
	     loaded_xsec.CrossSection(data, (double) (12. / 5. ) * 1e4).Contents(),
	     (Eigen::Array<double, 1, 10>::Ones()),
	     0);

  TEST_ARRAY("loaded xsec unfolded", 
	     loaded_xsec.UnfoldedCrossSection(data, (double) (12. / 5. ) * 1e4).Contents(),
	     (Eigen::Array<double, 1, 10>::Ones()*2),
	     0);
  

  return pass;
}