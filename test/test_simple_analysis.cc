#include <iostream>

#include "XSecAna/Hist.h"
#include "XSecAna/SimpleSignalEstimator.h"
#include "XSecAna/SimpleEfficiency.h"
#include "XSecAna/SimpleFlux.h"
#include "XSecAna/IUnfold.h"
#include "XSecAna/CrossSection.h"
#include "XSecAna/Analysis.h"
#include "test_utils.h"

#include <Eigen/Dense>
#include "TFile.h"
#include "TDirectory.h"

using namespace xsec;

auto test_file_name = test::utils::test_dir() + "test_simple_analysis.root";

std::unique_ptr<IMeasurement>
LoadSimpleCrossSection(TDirectory * dir,
                       const std::string & name) {
    return CrossSection::LoadFrom(SimpleEfficiency::LoadFrom,
                                  SimpleSignalEstimator::LoadFrom,
                                  SimpleFlux::LoadFrom,
                                  IdentityUnfold::LoadFrom,
                                  dir,
                                  name);
}

int main(int argc, char ** argv)
{
  bool verbose = false;
  if(argc > 1 && std::strcmp(argv[1], "-v") == 0)  verbose = true;
  bool pass = true;
  bool test;

  Hist ones = test::utils::get_hist_of_ones();

  auto data = test::utils::get_simple_data();


  auto hnominal = test::utils::get_simple_nominal_hist();
  auto hup = test::utils::get_simple_up_hist();
  auto hdown = test::utils::get_simple_down_hist();


  auto nominal_xsec = test::utils::make_simple_xsec(hnominal);
  auto up   = test::utils::make_simple_xsec(hup);
  auto down = test::utils::make_simple_xsec(hdown);

  auto nuniverses = 50;
  auto xsec_universes = test::utils::make_simple_xsec_multiverse(hnominal, nuniverses);
  Systematic<IMeasurement> syst_mv("mv", xsec_universes);
  Systematic<IMeasurement> syst_1 ("1sided", up);
  Systematic<IMeasurement> syst_2 ("2sided", up, down);

  std::vector<Systematic<IMeasurement>> systs = {
    syst_mv,
    syst_1,
    syst_2,
  };


  auto output = new TFile(test_file_name.c_str(), "recreate");
  SaveAnalysis(nominal_xsec,
               systs,
               data,
               output,
               "analysis");
  output->Close();
  delete output;

  auto input = TFile::Open(test_file_name.c_str());
  auto [loaded_nominal, loaded_systematics, loaded_data] = LoadAnalysis(LoadSimpleCrossSection,
                                                                        input,
                                                                        "analysis");
  assert(systs.size() == loaded_systematics.size());

  input->Close();
  delete input;

  TEST_HISTS_SAME("save/load",
                  loaded_nominal->Eval(loaded_data),
                  nominal_xsec->Eval(data),
                  0);

  auto eval = [](Hist & data) {
      ForEachFunction<Hist, IMeasurement> f = [&data](IMeasurement * m) {
          return new Hist(m->Eval(data));
      };
      return f;
  };

  for(auto i = 0u; i < systs.size(); i++) {
      assert(systs[i].GetType() == loaded_systematics[i].GetType());
      assert(systs[i].GetName() == loaded_systematics[i].GetName());

      auto loaded_systematic_eval = loaded_systematics[i].ForEach(eval(loaded_data));
      auto systematic_eval = systs[i].ForEach(eval(data));

      if (systs[i].GetType() == kMultiverse) {
          loaded_systematic_eval = Systematic<Hist>(loaded_systematics[i].GetName(),
                                                        new Hist(MultiverseShift(loaded_systematic_eval,
                                                                                     loaded_nominal->Eval(loaded_data),
                                                                                     1)),
                                                        new Hist(MultiverseShift(loaded_systematic_eval,
                                                                                     loaded_nominal->Eval(loaded_data),
                                                                                     -1)));
          systematic_eval = Systematic<Hist>(systs[i].GetName(),
                                                 new Hist(MultiverseShift(systematic_eval,
                                                                              loaded_nominal->Eval(loaded_data),
                                                                              1)),
                                                 new Hist(MultiverseShift(systematic_eval,
                                                                              loaded_nominal->Eval(loaded_data),
                                                                              -1)));
      }

      TEST_HISTS_SAME(("systematic " + std::to_string(i) + " up"),
                      (*loaded_systematic_eval.Up()),
                      (*systematic_eval.Up()),
                      0);
      if (systs[i].GetType() != kOneSided) {
          TEST_HISTS_SAME(("systematic " + std::to_string(i) + " down"),
                          (*loaded_systematic_eval.Down()),
                          (*systematic_eval.Down()),
                          0);
      }
  }

  return !pass;
}
