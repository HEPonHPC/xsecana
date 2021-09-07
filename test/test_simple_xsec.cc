#include "XSecAna/Hist.h"
#include "XSecAna/SimpleSignalEstimator.h"
#include "XSecAna/SimpleEfficiency.h"
#include "XSecAna/SimpleFlux.h"
#include "XSecAna/IUnfold.h"
#include "XSecAna/CrossSection.h"
#include "test_utils.h"


#include <iostream>

#include "TFile.h"

using namespace xsec;

typedef Hist<double, 10> histtype;
typedef CrossSection<histtype> SimpleCrossSection;

int main(int argc, char ** argv)
{
  bool verbose = false;
  if(argc > 1 && std::strcmp(argv[1], "-v") == 0)  verbose = true;
  bool pass = true;
  bool test;

  Hist<double, 10> ones(Eigen::Array<double, 1, 10>::Ones(),
                        Eigen::Array<double, 1, 11>::LinSpaced(11, 0, 20));
  Hist<double, 10> bkgd = ones * 2;

  Hist<double, 10> data(Eigen::Array<double, 1, 10>::Ones() * 4,
                        Eigen::Array<double, 1, 11>::LinSpaced(11, 0, 20),
                        test::utils::data_exposure);

  Hist<double, 10> flux_hist = ones * 5;

  Hist<double, 10> eff_num = ones / 4.;

  Hist<double, 10> eff_den = ones;


  typedef Hist<double, 10> histtype;

  // why is this template deduction failing??
  auto efficiency = new SimpleEfficiency<histtype>(eff_num, eff_den); // = 1/4 (no exposure scaling)
  auto flux = new SimpleFlux(flux_hist);                              // = 5/2 (after scaling by data exposure)
  auto signal_estimator = new SimpleSignalEstimator(bkgd);            // = 3 (after scaling by data exposure)
  auto unfold = new IdentityUnfold<double, 10>(bkgd.Contents().size()); // = 1

  SimpleCrossSection xsec(efficiency,
			  signal_estimator,
			  flux,
			  unfold); // = 1 / 1e4
  
  // test fail before targets is set

  // now set targets
  xsec.SetNTargets(test::utils::ntargets);
  
  auto xsec_differential = xsec.ToDifferential();

  TEST_ARRAY("signal",
	     (signal_estimator->Signal(data).Contents()),
	     (Eigen::Array<double, 1, 10>::Ones() * 3),
	     0);
  TEST_ARRAY("efficiency",
	     (efficiency->Eval().Contents()),
	     (Eigen::Array<double, 1, 10>::Ones() / 4.),
	     0);

  TEST_ARRAY("xsec",
             xsec.Eval(data).Contents(),
	     (Eigen::Array<double, 1, 10>::Ones() * 24. / 5.),
	     0);

  TEST_ARRAY("xsec_differential",
             xsec_differential.Eval(data).Contents(),
	     (Eigen::Array<double, 1, 10>::Ones() / 2. * 24. / 5.),
	     0);

  std::string test_file_name = test::utils::test_dir() + "test_simple_xsec.root";
  TFile * output = new TFile(test_file_name.c_str(), "recreate");
  xsec.SaveTo(output, "xsec");
  output->Close();
  delete output;

  //TFile * input = TFile::Open(test_file_name.c_str());
  //auto loaded_xsec = *SimpleCrossSection::LoadFrom(input, "xsec");
  //input->Close();
  //delete input;
    //TEST_ARRAY("loaded xsec",
    //         loaded_xsec.Eval(data).Contents(),
	//     (Eigen::Array<double, 1, 10>::Ones() * 24. / 5.),
	//     0);

  TEST_ARRAY("exposure",
	     (Eigen::Array<double, 1, 1>::Ones() * xsec.Eval(data).Exposure()),
	     (Eigen::Array<double, 1, 1>::Ones() * test::utils::data_exposure),
	     0);

  TEST_ARRAY("exposure differential",
	     (Eigen::Array<double, 1, 1>::Ones() * xsec_differential.Eval(data).Exposure()),
	     (Eigen::Array<double, 1, 1>::Ones() * test::utils::data_exposure),
	     0);


  TEST_ARRAY("test_utils::make_simple_xsec",
	     (test::utils::make_simple_xsec(ones)->Eval(test::utils::get_simple_data<double, 10>()).Contents()),
	     ones.Contents(),
	     0);

  return !pass;
}
