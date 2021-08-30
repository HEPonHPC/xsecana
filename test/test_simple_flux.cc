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

  Hist<double, 10> flux_hist(Eigen::Array<double, 1, 10>::Ones(),
			     Eigen::Array<double, 1, 11>::LinSpaced(11, 0, 10));

  SimpleFlux flux(flux_hist);
  SimpleIntegratedFlux<Hist<double, 10> > integrated_flux(flux_hist);

  TEST_HIST("flux.ToHist()"           , flux.ToHist()           , flux_hist.Contents(), flux_hist.Edges(), 0);
  TEST_HIST("integrated_flux.ToHist()", integrated_flux.ToHist(), flux_hist.Contents(), flux_hist.Edges(), 0);
  
  TEST_HIST ("flux.operator/", (flux / flux_hist), flux_hist.Contents(), flux_hist.Edges(), 0);
  TEST_HIST ("integrated_flux.operator/", 
	     (integrated_flux / flux.ToHist()),
	     flux_hist.Contents() * 10, 
	     flux_hist.Edges(), 0);

  TEST_HIST ("flux.operator*", (flux * flux_hist), flux_hist.Contents(), flux_hist.Edges(), 0);
  TEST_HIST ("integrated_flux.operator*", 
	     (integrated_flux * flux.ToHist()),
	     flux_hist.Contents() * 10, 
	     flux_hist.Edges(), 0);

  std::string test_file_name = test::utils::test_dir() + "test_simple_flux.root";
  TFile * output = new TFile(test_file_name.c_str(), "recreate");
  flux.SaveTo(output, "flux");
  integrated_flux.SaveTo(output, "integrated_flux");
  output->Close();
  delete output;

  TFile * input = TFile::Open(test_file_name.c_str());
  auto loaded_flux = *SimpleFlux<Hist<double, 10> >::LoadFrom(input, "flux");
  auto loaded_integrated_flux = *SimpleIntegratedFlux<Hist<double, 10> >::LoadFrom(input, "integrated_flux");
  input->Close();
  delete input;

  TEST_HIST("loaded_flux", loaded_flux.ToHist(), flux_hist.Contents(), flux_hist.Edges(), 0);
  TEST_HIST("loaded_integrated_flux", loaded_integrated_flux.ToHist(), flux_hist.Contents(), flux_hist.Edges(), 0);
  
  return !pass;
}
