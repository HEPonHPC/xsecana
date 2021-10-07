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

typedef Hist<double, -1> histtype;
typedef CrossSection<histtype> SimpleCrossSection;

int main(int argc, char ** argv)
{
  bool verbose = false;
  if(argc > 1 && std::strcmp(argv[1], "-v") == 0)  verbose = true;
  bool pass = true;
  bool test;

  auto nbins_and_uof = 12;
  Eigen::Array<double, 1, -1> bins = Eigen::Array<double, 1, -1>::Zero(nbins_and_uof+1);
  for(auto i = 0u; i < nbins_and_uof+1; i++) {
      bins(i) = 2*i;
  }

  Hist<double, -1> ones(Eigen::Array<double, 1, 12>::Ones(nbins_and_uof), bins);
  Hist<double, -1> bkgd = ones * 2;
  Hist<double, -1> data(ones.ContentsAndUOF() * 4,
                        ones.EdgesAndUOF(),
                        test::utils::data_exposure);

  Hist<double, -1> flux_hist = ones * 5;

  Hist<double, -1> eff_num = ones / 4.;

  Hist<double, -1> eff_den = ones;

  typedef Hist<double, -1> histtype;

  // why is this template deduction failing??
  auto efficiency = new SimpleEfficiency<histtype>(eff_num, eff_den); // = 1/4 (no exposure scaling)
  auto flux = new SimpleFlux(flux_hist);                              // = 5/2 (after scaling by data exposure)
  auto signal_estimator = new SimpleSignalEstimator(bkgd);            // = 3 (after scaling by data exposure)
  auto unfold = new IdentityUnfold<double, -1>(bkgd.ContentsAndUOF().size()); // = 1

  SimpleCrossSection xsec(efficiency,
                          signal_estimator,
                          flux,
                          unfold); // = 1 / 1e4

  // test fail before targets is set

  // now set targets
  xsec.SetNTargets(test::utils::ntargets);
  
  auto xsec_differential = xsec.ToDifferential();

  TEST_ARRAY_SAME("signal",
                  (signal_estimator->Signal(data).Contents()),
                  (Eigen::Array<double, 1, -1>::Ones(data.Contents().size()) * 3),
                  0);
  TEST_ARRAY_SAME("efficiency",
                  (efficiency->Eval().Contents()),
                  (Eigen::Array<double, 1, -1>::Ones(data.Contents().size()) / 4.),
                  0);

  TEST_ARRAY_SAME("xsec",
                  xsec.Eval(data).Contents(),
                  (Eigen::Array<double, 1, -1>::Ones(data.Contents().size()) * 24. / 5.),
                  0);

  auto xsec_differential_result = xsec_differential.Eval(data);
  TEST_ARRAY_SAME("xsec_differential",
                  xsec_differential.Eval(data).Contents(),
                  (Eigen::Array<double, 1, -1>::Ones(data.Contents().size()) / 2. * 24. / 5.),
                  0);

  std::string test_file_name = test::utils::test_dir() + "test_simple_xsec.root";
  TFile * output = new TFile(test_file_name.c_str(), "recreate");
  xsec.SaveTo(output, "xsec");
  output->Close();
  delete output;

  TFile * input = TFile::Open(test_file_name.c_str());
  auto loaded_xsec =
          SimpleCrossSection::LoadFrom(SimpleEfficiency<histtype>::LoadFrom,
                                       SimpleSignalEstimator<histtype>::LoadFrom,
                                       SimpleFlux<histtype>::LoadFrom,
                                       IdentityUnfold<double, -1>::LoadFrom,
                                       input,
                                       "xsec");
  input->Close();
  delete input;
  TEST_HISTS_SAME("loaded xsec",
                  loaded_xsec->Eval(data),
                  xsec.Eval(data),
                  0);

  auto simple_data = test::utils::get_simple_data<double, -1>();
  auto simple_ones = simple_data / simple_data;
  auto result = test::utils::make_simple_xsec(simple_ones)->Eval(test::utils::get_simple_data<double, -1>());
  TEST_HISTS_SAME("test_utils::make_simple_xsec",
                  simple_ones,
                  result,
                  0);

  return !pass;
}
