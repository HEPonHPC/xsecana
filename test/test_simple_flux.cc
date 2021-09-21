#include "XSecAna/Hist.h"
#include "XSecAna/SimpleFlux.h"
#include "test_utils.h"
#include <iostream>

#include "TFile.h"

using namespace xsec;

int main(int argc, char ** argv)
{
  bool verbose = false;
  if(argc > 1 && std::strcmp(argv[1], "-v") == 0)  verbose = true;
  bool pass = true;
  bool test;

  Hist<double, 10> flux_hist = test::utils::get_hist_of_ones<double, 10>();

  SimpleFlux flux(flux_hist);
  SimpleIntegratedFlux<Hist<double, 10> > integrated_flux(flux_hist);

  TEST_HISTS_SAME("flux.Eval()" ,
                  flux.Eval(),
                  flux_hist,
                  0);
  TEST_HISTS_SAME("integrated_flux.Eval()",
                  integrated_flux.Eval(),
                  flux_hist,
                  0);
  
  TEST_HISTS_SAME("flux.operator/",
                  (flux / flux_hist),
                  flux_hist,
                  0);
  TEST_HISTS_SAME("integrated_flux.operator/",
                  (integrated_flux / flux.Eval()),
                  flux_hist * 12,
                  0);

  TEST_HISTS_SAME("flux.operator*",
                  (flux * flux_hist),
                  flux_hist,
                  0);
  TEST_HISTS_SAME("integrated_flux.operator*",
                  (integrated_flux * flux.Eval()),
                  flux_hist * 12,
                  0);

  std::string test_file_name = test::utils::test_dir() + "test_simple_flux.root";
  TFile * output = new TFile(test_file_name.c_str(), "recreate");
  flux.SaveTo(output, "flux");
  integrated_flux.SaveTo(output, "integrated_flux");
  output->Close();
  delete output;

  TFile * input = TFile::Open(test_file_name.c_str());
  auto loaded_flux = IFlux<Hist<double, 10> >::LoadFrom(SimpleFlux<Hist<double, 10> >::LoadFrom,
                                                        input,
                                                        "flux");
  auto loaded_integrated_flux = IFlux<Hist<double, 10> >::LoadFrom(SimpleIntegratedFlux<Hist<double, 10>>::LoadFrom,
                                                                   input,
                                                                   "integrated_flux");
  input->Close();
  delete input;

  TEST_HISTS_SAME("loaded_flux",
                  loaded_flux->Eval(),
                  flux_hist,
                  0);
  TEST_HISTS_SAME("loaded_integrated_flux",
                  loaded_integrated_flux->Eval(),
                  flux_hist,
                  0);
  
  return !pass;
}
