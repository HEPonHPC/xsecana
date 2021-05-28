#include "XSecAna/Hist.h"
#include "XSecAna/SimpleFlux.h"

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

  TFile * output = new TFile("test_simple_flux.root", "recreate");
  flux.SaveTo(output, "flux");
  integrated_flux.SaveTo(output, "integrated_flux");
  output->Close();
  delete output;

  TFile * input = TFile::Open("test_simple_flux.root");
  auto loaded_flux = *SimpleFlux<Hist<double, 10> >::LoadFrom(input, "flux");
  auto loaded_integrated_flux = *SimpleIntegratedFlux<Hist<double, 10> >::LoadFrom(input, "integrated_flux");
  input->Close();
  delete input;

  TEST_HIST("loaded_flux", loaded_flux.ToHist(), flux_hist.Contents(), flux_hist.Edges(), 0);
  TEST_HIST("loaded_integrated_flux", loaded_integrated_flux.ToHist(), flux_hist.Contents(), flux_hist.Edges(), 0);
  
  if(pass) std::cout << "Success!" << std::endl;
  
}
