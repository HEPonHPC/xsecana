#include "XSecAna/Hist.h"
#include "XSecAna/SimpleSignalEstimator.h"
#include "XSecAna/SimpleEfficiency.h"
#include "XSecAna/SimpleFlux.h"
#include "XSecAna/IUnfold.h"
#include "XSecAna/CrossSection.h"
#include "XSecAna/SimpleQuadSum.h"
#include "test_utils.h"

#include <iostream>

#include "TFile.h"

using namespace xsec;

typedef Hist<double, 10> histtype;
typedef CrossSection<histtype> SimpleCrossSection;
typedef Systematic<IMeasurement<histtype>> CrossSectionSystematic;

int main(int argc, char ** argv)
{
  bool verbose = false;
  if(argc > 1 && std::strcmp(argv[1], "-v") == 0)  verbose = true;
  bool pass = true;
  bool test;

  Hist<double, 10> hone(Eigen::Array<double, 1, 10>::Ones(),
                        Eigen::Array<double, 1, 11>::LinSpaced(11, 0, 10));
  auto hnominal = test::utils::get_simple_nominal_hist<double, 10>();
  auto hup = test::utils::get_simple_up_hist<double, 10>();
  auto hdown = test::utils::get_simple_down_hist<double, 10>();
  auto hmax_shift = hnominal;
  for(auto i = 0u; i < hmax_shift.size(); i++) {
    hmax_shift[i] = std::max(std::abs(hnominal[i] - hup[i]),
                             std::abs(hnominal[i] - hdown[i]));
  }

  auto data = test::utils::get_simple_data<double, 10>();

  auto nominal_xsec = test::utils::make_simple_xsec(hnominal);
  auto up   = test::utils::make_simple_xsec(hup);
  auto down = test::utils::make_simple_xsec(hdown);


  // simple test of the math
  Systematic<histtype> one("1", new histtype(hone));
  Systematic<histtype> four("4", new histtype(hone * 4));
  Systematic<histtype> three("3", new histtype(hone * 3));
  auto two  = new histtype(hone * 2);


  auto one_half = SimpleQuadSum::FractionalUncertainty<histtype>(two,
                                                                 three,
                                                                 data);
  TEST_HISTS_SAME("one_half",
                  *(std::get<1>(one_half).Up()),
                  (hone / 2),
                  1e-14);


  std::map<std::string, Systematic<histtype> > syst_map = {
    {"1", one},
    {"3", three},
    {"4", four},
  };
  auto sqrt_six_halves = SimpleQuadSum::TotalFractionalUncertainty<histtype>(two,
                                                                             syst_map,
                                                                             data);
  TEST_HISTS_SAME("sqrt_six_halves",
                  *(std::get<1>(sqrt_six_halves).Up()),
                  (hone * 6 / 4).sqrt(),
                  1e-14);



  /*
    multiverse example
  */
  int nuniverses = 50;

  std::vector<IMeasurement<histtype>*> xsec_universes = test::utils::make_simple_xsec_multiverse(hnominal, nuniverses);
  std::vector<Hist<double, 10>*> hist_universes = test::utils::make_simple_hist_multiverse(hnominal, nuniverses);

  Systematic<IMeasurement<histtype>> syst_mv("mv_xsec", xsec_universes);
  Systematic<Hist<double, 10> > syst_mv_hist("mv_hist", hist_universes);


  auto abs_uncert_mv = SimpleQuadSum::AbsoluteUncertainty<histtype>(nominal_xsec,
                                                                    syst_mv,
                                                                    data);

  // index of universe representing minus 1 sigma shift
  int m1_idx = (0.5 - std::erf(1 / std::sqrt(2)) / 2.0) * (nuniverses-1) + 1;

  auto target_minus_1_sigma_multiverse = MultiverseShift(syst_mv_hist, hnominal, 1);
  TEST_HISTS_SAME("minus 1 sigma multiverse",
                  *(std::get<1>(abs_uncert_mv).Up()),
                  target_minus_1_sigma_multiverse,
                  1e-14);

  /*
    Examples using 1 and 2 sided shifts
    Test each function of the propagator
  */
  Systematic<IMeasurement<histtype>> syst_1sided("1sided", up);
  Systematic<IMeasurement<histtype>> syst_2sided("2sided", up, down);

  // AbsoluteUncertainty
  auto abs_uncert_1sided = SimpleQuadSum::AbsoluteUncertainty<histtype>(nominal_xsec,
                                                                        syst_1sided,
                                                                        data);

  auto target_abs_uncert_1_sided = hup - hnominal;
  TEST_HISTS_SAME("abs_uncert 1 sided",
                  *(std::get<1>(abs_uncert_1sided).Up()),
                  target_abs_uncert_1_sided,
                  1e-14);

  std::map<std::string, Systematic<IMeasurement<histtype>> > rmap = {
    {"1sided", syst_1sided},
  };
  auto symmetrize = SimpleQuadSum::TotalAbsoluteUncertainty<histtype>(nominal_xsec,
                                                                      rmap,
                                                                      data);
  TEST_HISTS_SAME("symmeterize",
                  *(std::get<1>(symmetrize).Up()),
                  *(std::get<1>(symmetrize).Down()),
                  0);

  auto abs_uncert_2sided = SimpleQuadSum::AbsoluteUncertainty<histtype>(nominal_xsec,
                                                                        syst_2sided,
                                                                        data);

  TEST_HISTS_SAME("abs_uncert 2 sided",
                  *(std::get<1>(abs_uncert_2sided).Up()),
                  hmax_shift,
                  1e-14);

  // FractionalUncertainty
  auto frac_uncert_1sided = SimpleQuadSum::FractionalUncertainty<histtype>(nominal_xsec,
                                                                           syst_1sided,
                                                                           data);
  TEST_HISTS_SAME("fractional uncert",
                  *(std::get<1>(frac_uncert_1sided).Up()),
                  ((hup - hnominal) / hnominal),
                  1e-14);


  // TotalAbsoluteUncertainty
  std::map<std::string, Systematic<IMeasurement<histtype>> > systs = {
    {"1sided", syst_1sided},
    {"2sided", syst_2sided},
    {"mv", syst_mv},
  };

  auto total_abs_uncert = SimpleQuadSum::TotalAbsoluteUncertainty<histtype>(nominal_xsec,
                                                                            systs,
                                                                            data);
  auto target_total_abs_uncert = (std::get<1>(abs_uncert_mv).Up()->pow(2) +
                                  std::get<1>(abs_uncert_1sided).Up()->pow(2) +
                                  std::get<1>(abs_uncert_2sided).Up()->pow(2)).sqrt();

  TEST_HISTS_SAME("total absolute uncert",
                  *(std::get<1>(total_abs_uncert).Up()),
                  target_total_abs_uncert,
                  1e-14);

  // TotalFractionalUncertainty
  auto total_frac_uncert = SimpleQuadSum::TotalFractionalUncertainty<histtype>(nominal_xsec,
                                                                               systs,
                                                                               data);
  TEST_HISTS_SAME("total frac uncert",
                  *(std::get<1>(total_frac_uncert).Up()),
                  (target_total_abs_uncert / hnominal),
                  1e-14);

  return !pass;
}
