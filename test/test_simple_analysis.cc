#include <iostream>
#include <stdio.h>

#include "XSecAna/Hist.h"
#include "XSecAna/SimpleSignalEstimator.h"
#include "XSecAna/SimpleEfficiency.h"
#include "XSecAna/SimpleFlux.h"
#include "XSecAna/IUnfold.h"
#include "XSecAna/CrossSection.h"
#include "XSecAna/SimpleQuadSum.h"
#include "XSecAna/Analysis.h"
#include "test_utils.h"

#include <Eigen/Dense>
#include "TFile.h"

using namespace xsec;

std::string test_file_name = test::utils::test_dir() + "test_simple_analysis.root";

typedef Analysis<test::utils::SimpleCrossSection,
			     SimpleQuadSum<test::utils::SimpleCrossSection,
					   Hist<double, 10> > ,
			     Hist<double, 10> > SimpleCrossSectionAnalysis;

int main(int argc, char ** argv)
{
  bool verbose = false;
  if(argc > 1 && std::strcmp(argv[1], "-v") == 0)  verbose = true;
  bool pass = true;
  bool test;
  
  Hist<double, 10> ones(Eigen::Array<double, 1, 10>::Ones(),
			Eigen::Array<double, 1, 11>::LinSpaced(11, 0, 10));

  auto data = test::utils::get_simple_data<double, 10>();


  auto hnominal = test::utils::get_simple_nominal_hist<double, 10>();
  auto hup = test::utils::get_simple_up_hist<double, 10>();
  auto hdown = test::utils::get_simple_down_hist<double, 10>();


  test::utils::SimpleCrossSection nominal_xsec = test::utils::make_simple_xsec(hnominal);
  test::utils::SimpleCrossSection up   = test::utils::make_simple_xsec(hup);  
  test::utils::SimpleCrossSection down = test::utils::make_simple_xsec(hdown);

  auto nuniverses = 50;
  Systematic<test::utils::SimpleCrossSection> syst_mv("mv"    , test::utils::make_simple_xsec_multiverse(hnominal, nuniverses));
  Systematic<test::utils::SimpleCrossSection> syst_1 ("1sided", up);
  Systematic<test::utils::SimpleCrossSection> syst_2 ("2sided", up, down);
  
  std::map<std::string, Systematic<test::utils::SimpleCrossSection> > systs = {
    {"mv", syst_mv},
    {"1sided", syst_1},
    {"2sided", syst_2},
  };
  
  SimpleQuadSum<test::utils::SimpleCrossSection, Hist<double, 10> > prop;

  SimpleCrossSectionAnalysis analysis(nominal_xsec,
				      systs,
				      data);
  TEST_ARRAY("nominal closure", 
	     hnominal.Contents(),
	     analysis.Result().Contents(),
	     1e-14);

  TEST_ARRAY("total abs uncert",
	     prop.TotalAbsoluteUncertainty(data, nominal_xsec, systs).first.Contents(),
	     analysis.TotalAbsoluteUncertainty().first.Contents(),
	     0);

  TEST_ARRAY("total frac uncert",
	     prop.TotalFractionalUncertainty(data, nominal_xsec, systs).first.Contents(),
	     analysis.TotalFractionalUncertainty().first.Contents(),
	     0);
	     


  TFile * output = new TFile(test_file_name.c_str(), "recreate");
  analysis.SaveTo(output, "analysis");
  
  auto results_dir = output->mkdir("results");
  analysis.Result().SaveTo(results_dir,  "nominal");
  analysis.Result("mv").SaveTo(results_dir, "mv");
  analysis.Result("1sided").SaveTo(results_dir, "1sided");
  analysis.Result("2sided").SaveTo(results_dir, "2sided");
  
  auto abs_errors = prop.TotalAbsoluteUncertainty(data, nominal_xsec, systs);
  (abs_errors.first + hnominal).SaveTo(results_dir, "total_up");
  (hnominal - abs_errors.second).SaveTo(results_dir, "total_down");
  
  output->Close();
  delete output;
  
  TFile * input = TFile::Open(test_file_name.c_str());
  auto loaded_analysis = *SimpleCrossSectionAnalysis::LoadFrom(input, "analysis");
  input->Close();
  delete input;

  TEST_ARRAY("save/load",
	     analysis.Result().Contents(), 
	     loaded_analysis.Result().Contents(),
	     0);

  return !pass;
}
