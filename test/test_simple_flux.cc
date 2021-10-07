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

  auto flux_hist = test::utils::get_hist_of_ones<double, -1>();
  auto flux_binning = flux_hist.EdgesAndUOF();

  HistXd flux_histxd(flux_hist.ContentsAndUOF(),
                     flux_hist.EdgesAndUOF(),
                     flux_hist.ErrorsAndUOF(),
                     flux_hist.Exposure());

  // arbitrary binning
  auto analysis_binning = Hist<double, -1>::edges_type::LinSpaced(43, 0, 41);
  Hist<double, -1> ones_analysis_binning(Hist<double, -1>::array_and_uof_type::Ones(42),
                                         analysis_binning);

  if(flux_binning.size() == analysis_binning.size()) {
      assert(!(flux_binning - analysis_binning).isZero(1e-5));
  }

  SimpleFlux flux(flux_hist);
  SimpleIntegratedFlux<Hist<double, -1> > integrated_flux(flux_histxd);


  TEST_HISTS_SAME("flux.Eval()" ,
                  flux.Eval(flux_binning),
                  flux_hist,
                  0);
  auto target_integrated_flux = ones_analysis_binning * flux_histxd.ContentsAndUOF().size();
  target_integrated_flux.SetErrorsAndUOF(ones_analysis_binning.ContentsAndUOF() *
                                         std::sqrt(flux_histxd.ContentsAndUOF().size()));
  TEST_HISTS_SAME("integrated_flux.Eval()",
                  integrated_flux.Eval(analysis_binning),
                  target_integrated_flux,
                  0);
  

  std::string test_file_name = test::utils::test_dir() + "test_simple_flux.root";
  auto output = new TFile(test_file_name.c_str(), "recreate");
  flux.SaveTo(output, "flux");
  integrated_flux.SaveTo(output, "integrated_flux");
  output->Close();
  delete output;

  TFile * input = TFile::Open(test_file_name.c_str());
  auto loaded_flux = IFlux<Hist<double, -1> >::LoadFrom(SimpleFlux<Hist<double, -1> >::LoadFrom,
                                                        input,
                                                        "flux");
  auto loaded_integrated_flux = IFlux<Hist<double, -1> >::LoadFrom(SimpleIntegratedFlux<Hist<double,-1>>::LoadFrom,
                                                                   input,
                                                                   "integrated_flux");
  input->Close();
  delete input;

  TEST_HISTS_SAME("loaded_flux",
                  loaded_flux->Eval(flux_binning),
                  flux_hist,
                  0);
  TEST_HISTS_SAME("loaded_integrated_flux",
                  loaded_integrated_flux->Eval(analysis_binning),
                  target_integrated_flux,
                  0);
  
  return !pass;
}
