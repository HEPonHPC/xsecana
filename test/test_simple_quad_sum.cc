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
typedef Systematic<IMeasurement<histtype>*> CrossSectionSystematic;

int main(int argc, char ** argv)
{
  bool verbose = false;
  if(argc > 1 && std::strcmp(argv[1], "-v") == 0)  verbose = true;
  bool pass = true;
  bool test;

  SimpleQuadSum<histtype, IMeasurement<histtype>*, histtype> prop;

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
  CrossSectionSystematic one("1", test::utils::make_simple_xsec(hone));
  CrossSectionSystematic four("4",
                              test::utils::make_simple_xsec(hone * 4));
  CrossSectionSystematic three("3",
                               test::utils::make_simple_xsec(hone * 3));
  auto two  = test::utils::make_simple_xsec(hone * 2);


  auto one_half = prop.FractionalUncertainty(two,
                                             three,
                                             data);
  TEST_ARRAY("one_half", one_half.Contents(), (hone / 2).Contents(), 1e-14);


  std::map<std::string, Systematic<IMeasurement<histtype>*> > syst_map = {
    {"1", one},
    {"3", three},
    {"4", four},
  };
  auto sqrt_six_halves = prop.TotalFractionalUncertainty(two,
                                                         syst_map,
                                                         data).first;
  TEST_ARRAY("sqrt_six_halves", sqrt_six_halves.Contents(), (hone * 6 / 4).sqrt().Contents(), 1e-14);



  /*
    multiverse example
  */
  int nuniverses = 50;

  std::vector<IMeasurement<histtype>*> xsec_universes = test::utils::make_simple_xsec_multiverse(hnominal, nuniverses);
  std::vector<Hist<double, 10> > hist_universes = test::utils::make_simple_hist_multiverse(hnominal, nuniverses);

  Systematic<IMeasurement<histtype>*> syst_mv("mv_xsec", xsec_universes);
  Systematic<Hist<double, 10> > syst_mv_hist("mv_hist", hist_universes);

  auto abs_uncert_mv = prop.AbsoluteUncertainty(nominal_xsec,
                                                syst_mv,
                                                data);

  // index of universe representing minus 1 sigma shift
  int m1_idx = (0.5 - std::erf(1 / std::sqrt(2)) / 2.0) * (nuniverses-1) + 1;

  TEST_ARRAY("minus 1 sigma multiverse",
  	     abs_uncert_mv.Contents(),
  	     (hist_universes[m1_idx]-hnominal).Contents().abs(), 1e-14);


  /*
    Examples using 1 and 2 sided shifts
    Test each function of the propogator
  */
  Systematic<IMeasurement<histtype>*> syst_1sided("1sided", up);
  Systematic<IMeasurement<histtype>*> syst_2sided("2sided", up, down);

  // AbsoluteUncertainty
  auto abs_uncert_1sided = prop.AbsoluteUncertainty(nominal_xsec,
                                                    syst_1sided,
                                                    data);

  TEST_ARRAY("abs_uncert 1 sided",
	     abs_uncert_1sided.Contents(),
	     (hup - hnominal).Contents(),
	     1e-14);

  std::map<std::string, Systematic<IMeasurement<histtype> *> > rmap = {
    {"1sided", syst_1sided},
  };
  auto symmetrize = prop.TotalAbsoluteUncertainty(nominal_xsec,
                                                  rmap,
                                                  data);
  TEST_ARRAY("symmeterize",
             symmetrize.first.Contents(),
             symmetrize.second.Contents(),
             0);

  auto abs_uncert_2sided = prop.AbsoluteUncertainty(nominal_xsec,
                                                    syst_2sided,
                                                    data);
  TEST_ARRAY("abs_uncert 2 sided",
             abs_uncert_2sided.Contents(),
             hmax_shift.Contents(),
             1e-14);

  // FractionalUncertainty
  auto frac_uncert_1sided = prop.FractionalUncertainty(nominal_xsec,
                                                       syst_1sided,
                                                       data);
  TEST_ARRAY("fractional uncert",
             frac_uncert_1sided.Contents(),
             ((hup - hnominal) / hnominal).Contents(),
             1e-14);


  // TotalAbsoluteUncertainty
  std::map<std::string, Systematic<IMeasurement<histtype>*> > systs = {
    {"1sided", syst_1sided},
    {"2sided", syst_2sided},
    {"mv", syst_mv},
  };

  auto total_abs_uncert = prop.TotalAbsoluteUncertainty(nominal_xsec,
                                                        systs,
                                                        data);
  auto target_total_abs_uncert = (abs_uncert_mv.pow(2) + abs_uncert_1sided.pow(2) + abs_uncert_2sided.pow(2)).sqrt();

  TEST_ARRAY("total absolute uncert",
             total_abs_uncert.first.Contents(),
             target_total_abs_uncert.Contents(),
             1e-14);

  // TotalFractionalUncertainty
  auto total_frac_uncert = prop.TotalFractionalUncertainty(nominal_xsec,
                                                           systs,
                                                           data);
  TEST_ARRAY("total frac uncert",
             total_frac_uncert.first.Contents(),
             (abs_uncert_mv.pow(2) + abs_uncert_1sided.pow(2) + abs_uncert_2sided.pow(2)).Contents().sqrt() / hnominal.Contents(),
             1e-14);

  return !pass;
}
