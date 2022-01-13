
#include "XSecAna/CrossSection.h"
#include "XSecAna/SimpleQuadSum.h"
#include "test_utils.h"
#include <iostream>

#include "TFile.h"

using namespace xsec;
const static double ntargets = 1e4;

/////////////////////////////////////////////////////////
// make a cross section object that evaluates out to
//  the input array when folded
std::shared_ptr<IMeasurement> make_simple_xsec(std::shared_ptr<const TH1> val) {
    auto ones = root::ToROOTLike(val.get(), Array::Ones(val->GetNbinsX() + 2));

    auto signal = (TH1 *) test::utils::get_simple_signal()->Clone();
    signal->Divide(val.get());
    auto efficiency = new SimpleEfficiency(signal,
                                           ones);
    auto flux = new SimpleFlux(ones);
    auto signal_estimator = new SimpleSignalEstimator(test::utils::get_simple_background());
    auto unfold = new IdentityUnfolder(ones->GetNbinsX() + 2);
    auto ret = std::shared_ptr<IMeasurement>(new EigenCrossSectionEstimator(efficiency,
                                                                            signal_estimator,
                                                                            flux,
                                                                            unfold,
                                                                            ntargets));
    return ret;
}

std::vector<std::shared_ptr<IMeasurement>> make_simple_xsec_multiverse(const TH1 * hnominal, int nuniverses) {
    auto hists = test::utils::make_simple_hist_multiverse(hnominal, nuniverses);
    std::vector<std::shared_ptr<IMeasurement>> xsecs(nuniverses);
    for (auto i = 0; i < hists.size(); i++) {
        xsecs[i] = make_simple_xsec(hists[i]);
    }
    return xsecs;
}

int main(int argc, char ** argv) {
    bool verbose = false;
    if (argc > 1 && std::strcmp(argv[1], "-v") == 0) verbose = true;
    bool pass = true;
    bool test;

    auto hone = std::shared_ptr<TH1>(test::utils::get_hist_of_ones());
    auto hnominal = std::shared_ptr<TH1>(test::utils::get_simple_nominal_hist());
    auto hup = std::shared_ptr<TH1>(test::utils::get_simple_up_hist());
    auto hdown = std::shared_ptr<const TH1>(test::utils::get_simple_down_hist());
    auto hmax_shift = std::shared_ptr<TH1>((TH1 *) hnominal->Clone());
    for (auto i = 0u; i <= hmax_shift->GetNbinsX() + 1; i++) {
        hmax_shift->SetBinContent(i, std::max(std::abs(hnominal->GetBinContent(i) - hup->GetBinContent(i)),
                                              std::abs(hnominal->GetBinContent(i) - hdown->GetBinContent(i))));
        hmax_shift->SetBinError(i, 0);
    }

    auto data = test::utils::get_simple_data();

    auto nominal_xsec = make_simple_xsec(hnominal);
    auto up = make_simple_xsec(hup);
    auto down = make_simple_xsec(hdown);


    // simple test of the math
    Systematic<TH1> one("1", hone);
    Systematic<TH1> four("4", std::shared_ptr<TH1>(test::utils::make_constant_hist_like(hone.get(), 4)));
    Systematic<TH1> three("3", std::shared_ptr<TH1>(test::utils::make_constant_hist_like(hone.get(), 3)));
    auto two = test::utils::make_constant_hist_like(hone.get(), 2);


    auto expected_one = (TH1 *) three.Up()->Clone();
    expected_one->Add(two, -1);
    expected_one->SetError(Array::Zero(data->GetNbinsX() + 2).eval().data());

    auto abs_uncert_one = SimpleQuadSum::AbsoluteUncertainty(two, three, data);
    pass &= TEST_HIST("abs_uncert_one",
                      std::get<1>(abs_uncert_one).Up().get(),
                      expected_one,
                      1e-14,
                      verbose);

    auto one_half = SimpleQuadSum::FractionalUncertainty<TH1>(two,
                                                              three,
                                                              data);

    auto expected_one_half = (TH1 *) three.Up()->Clone();
    expected_one_half->Add(two, -1);
    expected_one_half->Divide(two);
    expected_one_half->SetError(Array::Zero(data->GetNbinsX() + 2).eval().data());

    pass &= TEST_HIST("one_half",
                      std::get<1>(one_half).Up().get(),
                      expected_one_half,
                      1e-14,
                      verbose);


    std::map<std::string, Systematic<TH1> > syst_map = {
            {"1", one},
            {"3", three},
            {"4", four},
    };
    auto sqrt_six_halves = SimpleQuadSum::TotalFractionalUncertainty<TH1>(two,
                                                                          syst_map,
                                                                          data);
    auto expected_sqrt_six_halves = (TH1 *) two->Clone();
    for (auto i = 0u; i <= two->GetNbinsX() + 1; i++) {
        auto o = one.Up()->GetBinContent(i);
        auto th = three.Up()->GetBinContent(i);
        auto f = four.Up()->GetBinContent(i);
        auto tw = two->GetBinContent(i);
        expected_sqrt_six_halves->SetBinContent(i, std::sqrt(std::pow(o - tw, 2) +
                                                             std::pow(th - tw, 2) +
                                                             std::pow(f - tw, 2)) /
                                                   tw);
        expected_sqrt_six_halves->SetBinError(i, 0);
    }

    pass &= TEST_HIST("sqrt_six_halves",
                      (std::get<1>(sqrt_six_halves).Up().get()),
                      expected_sqrt_six_halves,
                      1e-14,
                      verbose);



    //
    //multiverse example
    //
    int nuniverses = 50;

    std::vector<std::shared_ptr<IMeasurement>> xsec_universes = make_simple_xsec_multiverse(hnominal.get(), nuniverses);
    std::vector<std::shared_ptr<TH1>> hist_universes = test::utils::make_simple_hist_multiverse(hnominal.get(), nuniverses);

    Systematic<IMeasurement> syst_mv("mv_xsec", xsec_universes);
    Systematic<TH1> syst_mv_hist("mv_hist", hist_universes);

    auto abs_uncert_mv = SimpleQuadSum::AbsoluteUncertainty(nominal_xsec.get(),
                                                            syst_mv,
                                                            data);

    TH1 * target_plus_1_sigma_multiverse = MultiverseShift(syst_mv_hist, hnominal.get(), 1);
    target_plus_1_sigma_multiverse->Add(hnominal.get(), -1);
    for(auto i = 0; i < target_plus_1_sigma_multiverse->GetNbinsX()+2; i++) {
        target_plus_1_sigma_multiverse->SetBinError(i, 0);
    }
    pass &= TEST_HIST("plus 1 sigma multiverse",
                      std::get<1>(abs_uncert_mv).Up().get(),
                      target_plus_1_sigma_multiverse,
                      1e-14,
                      verbose);



    //
    //  Examples using 1 and 2 sided shifts
    //  Test each function of the propagator
    //
    Systematic<IMeasurement> syst_1sided("1sided", up);
    Systematic<IMeasurement> syst_2sided("2sided", up, down);

    // AbsoluteUncertainty
    auto abs_uncert_1sided = SimpleQuadSum::AbsoluteUncertainty<IMeasurement>(nominal_xsec.get(),
                                                                              syst_1sided,
                                                                              data);

    auto target_abs_uncert_1_sided = (TH1 *) hup->Clone();
    target_abs_uncert_1_sided->Add(hnominal.get(), -1);
    target_abs_uncert_1_sided->SetError(Array::Zero(target_abs_uncert_1_sided->GetNbinsX() + 2).eval().data());
    pass &= TEST_HIST("abs_uncert 1 sided",
                      std::get<1>(abs_uncert_1sided).Up().get(),
                      target_abs_uncert_1_sided,
                      1e-14,
                      verbose);

    std::map<std::string, Systematic<IMeasurement> > rmap = {
            {"1sided", syst_1sided},
    };
    auto symmetrize = SimpleQuadSum::TotalAbsoluteUncertainty(nominal_xsec.get(),
                                                              rmap,
                                                              data);
    pass &= TEST_HIST("symmeterize",
                      std::get<1>(symmetrize).Up().get(),
                      std::get<1>(symmetrize).Down().get(),
                      0,
                      verbose);

    auto abs_uncert_2sided = SimpleQuadSum::AbsoluteUncertainty(nominal_xsec.get(),
                                                                syst_2sided,
                                                                data);

    pass &= TEST_HIST("abs_uncert 2 sided",
                      std::get<1>(abs_uncert_2sided).Up().get(),
                      hmax_shift.get(),
                      1e-14,
                      verbose);

    // FractionalUncertainty
    auto frac_uncert_1sided = SimpleQuadSum::FractionalUncertainty(nominal_xsec.get(),
                                                                   syst_1sided,
                                                                   data);
    auto expected_frac_uncert = (TH1*) hup->Clone();
    expected_frac_uncert->Add(hnominal.get(), -1);
    expected_frac_uncert->Divide(hnominal.get());
    expected_frac_uncert->SetError(Array::Zero(hup->GetNbinsX()+2).eval().data());
    pass &= TEST_HIST("fractional uncert",
                      std::get<1>(frac_uncert_1sided).Up().get(),
                      expected_frac_uncert,
                      1e-14,
                      verbose);


    // TotalAbsoluteUncertainty
    std::map<std::string, Systematic<IMeasurement> > systs = {
            {"1sided", syst_1sided},
            {"2sided", syst_2sided},
            {"mv",     syst_mv},
    };

    auto total_abs_uncert = SimpleQuadSum::TotalAbsoluteUncertainty(nominal_xsec.get(),
                                                                    systs,
                                                                    data);

    auto target_total_abs_uncert = (TH1*) hup->Clone();
    for(auto i = 0u; i <= hup->GetNbinsX()+1; i++) {
        target_total_abs_uncert->SetBinContent(i,
                                               std::sqrt(std::pow(std::get<1>(abs_uncert_mv).Up()->GetBinContent(i), 2) +
                                                         std::pow(std::get<1>(abs_uncert_1sided).Up()->GetBinContent(i), 2) +
                                                         std::pow(std::get<1>(abs_uncert_2sided).Up()->GetBinContent(i), 2)));
        target_total_abs_uncert->SetBinError(i, 0);
    }

    pass &= TEST_HIST("total absolute uncert",
                      std::get<1>(total_abs_uncert).Up().get(),
                      target_total_abs_uncert,
                      1e-14,
                      verbose);

    // TotalFractionalUncertainty
    auto total_frac_uncert = SimpleQuadSum::TotalFractionalUncertainty(nominal_xsec.get(),
                                                                       systs,
                                                                       data);
    auto target_total_frac_abs_uncert = (TH1*) target_total_abs_uncert->Clone();
    target_total_frac_abs_uncert->Divide(hnominal.get());
    target_total_frac_abs_uncert->SetError(Array::Zero(hup->GetNbinsX()+2).eval().data());
    pass &= TEST_HIST("total frac uncert",
                      std::get<1>(total_frac_uncert).Up().get(),
                      target_total_frac_abs_uncert,
                      1e-14,
                      verbose);

    return !pass;

}
