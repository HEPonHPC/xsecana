#include <iostream>
#include <cstdio>

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
#include "TDirectory.h"

using namespace xsec;

auto test_file_name = test::utils::test_dir() + "test_simple_analysis.root";

typedef xsec::Hist<double, 10> histtype;

std::unique_ptr<IMeasurement<histtype>>
LoadSimpleCrossSection(TDirectory * dir,
                       const std::string & name) {
    return CrossSection<histtype>::LoadFrom(SimpleEfficiency<histtype>::LoadFrom,
                                  SimpleSignalEstimator<histtype>::LoadFrom,
                                  SimpleFlux<histtype>::LoadFrom,
                                  IdentityUnfold<double, 10>::LoadFrom,
                                  dir,
                                  name);
}

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


  auto nominal_xsec = test::utils::make_simple_xsec(hnominal);
  auto up   = test::utils::make_simple_xsec(hup);
  auto down = test::utils::make_simple_xsec(hdown);

  typedef Hist<double, 10> histtype;
  auto nuniverses = 50;
  auto xsec_universes = test::utils::make_simple_xsec_multiverse(hnominal, nuniverses);
  Systematic<IMeasurement<histtype>> syst_mv("mv", xsec_universes);
  Systematic<IMeasurement<histtype>> syst_1 ("1sided", up);
  Systematic<IMeasurement<histtype>> syst_2 ("2sided", up, down);

  std::map<std::string, Systematic<IMeasurement<histtype>> > systs = {
    {"mv", syst_mv},
    {"1sided", syst_1},
    {"2sided", syst_2},
  };
  
  SimpleQuadSum<Hist<double, 10>,
                IMeasurement<histtype>,
                const Hist<double, 10> > prop;

  Analysis<Hist<double, 10>> analysis(nominal_xsec,
                                      systs,
                                      data,
                                      &prop);
  TEST_ARRAY("nominal closure", 
             hnominal.Contents(),
             analysis.Result().Contents(),
             1e-14);

  TEST_ARRAY("total abs uncert",
             prop.TotalAbsoluteUncertainty(nominal_xsec, systs, data).first.Contents(),
             analysis.TotalAbsoluteUncertainty().first.Contents(),
             0);

  TEST_ARRAY("total frac uncert",
             prop.TotalFractionalUncertainty(nominal_xsec, systs, data).first.Contents(),
             analysis.TotalFractionalUncertainty().first.Contents(),
             0);
	     


  TFile * output = new TFile(test_file_name.c_str(), "recreate");
  analysis.SaveTo(output, "analysis");
  
  auto results_dir = output->mkdir("results");
    analysis.Result().SaveTo(results_dir, "nominal");
    analysis.Result("mv").SaveTo(results_dir, "mv");
    analysis.Result("1sided").SaveTo(results_dir, "1sided");
    analysis.Result("2sided").SaveTo(results_dir, "2sided");
  
  auto abs_errors = prop.TotalAbsoluteUncertainty(nominal_xsec, systs, data);
  (abs_errors.first + hnominal).SaveTo(results_dir, "total_up");
  (hnominal - abs_errors.second).SaveTo(results_dir, "total_down");
  
  output->Close();
  delete output;

  auto input = TFile::Open(test_file_name.c_str());
  auto loaded_analysis = Analysis<histtype>::LoadFrom(LoadSimpleCrossSection,
                                                      input, "analysis");
  input->Close();
  delete input;

  TEST_ARRAY("save/load",
             analysis.Result().Contents(),
             loaded_analysis->Result().Contents(),
             0);

  return !pass;
}
