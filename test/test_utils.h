#pragma once

#include <iostream>
#include "TH1D.h"
#include "XSecAna/Utils.h"

#include "XSecAna/IUnfold.h"
#include "XSecAna/SimpleSignalEstimator.h"
#include "XSecAna/SimpleEfficiency.h"
#include "XSecAna/SimpleFlux.h"
#include "XSecAna/CrossSection.h"

[[nodiscard]] inline bool TEST_HIST(std::string test_name,
                                    const TH1 * HIST,
                                    const TH1 * target,
                                    double precision,
                                    bool verbose) {
    bool pass = true;
    xsec::Array _HIST_c = xsec::root::MapContentsToEigen(HIST);
    xsec::Array _HIST_e = xsec::root::MapErrorsToEigen(HIST);
    xsec::Array _target_c = xsec::root::MapContentsToEigen(target);
    xsec::Array _target_e = xsec::root::MapErrorsToEigen(target);

    bool test = (_HIST_c - _target_c).isZero(precision);
    if (!test || verbose) {
        std::cerr << __FUNCTION__ << "\t" << test_name + " (contents)" << (test ? ": PASSED" : ": FAILED")
                  << std::endl;
        std::cerr << __FUNCTION__ << "\t" << test_name + " (contents)" << "\t"
                  << _HIST_c.transpose() << std::endl;
        std::cerr << __FUNCTION__ << "\t" << test_name + " (contents)" << "\t" << _target_c.transpose()
                  << std::endl;
        pass &= test;
    }
    test = (_HIST_e - _target_e).isZero(precision);
    if (!test || verbose) {
        std::cerr << __FUNCTION__ << "\t" << test_name + " (errors)" << (test ? ": PASSED" : ": FAILED")
                  << std::endl;
        std::cerr << __FUNCTION__ << "\t" << test_name + " (errors)" << "\t" << _HIST_e.transpose()
                  << std::endl;
        std::cerr << __FUNCTION__ << "\t" << test_name + " (errors)" << "\t" << _target_e.transpose()
                  << std::endl;
        pass &= test;
    }
    return pass;
}

[[nodiscard]] inline bool TEST_ARRAY_SAME(const std::string & test_name,
                                          const xsec::Array & arr1,
                                          const xsec::Array & arr2,
                                          const double & precision,
                                          const bool & verbose) {
    bool test = (arr1 - arr2).isZero(precision);
    if (!test || verbose) {
        std::cerr << __FUNCTION__ << "\t" << test_name << (test ? ": PASSED" : ": FAILED") << std::endl;
        std::cerr << __FUNCTION__ << "\t" << test_name << "\t" << arr1.transpose() << std::endl;
        std::cerr << __FUNCTION__ << "\t" << test_name << "\t" << arr2.transpose() << std::endl;
    }
    return test;
}

[[nodiscard]] inline bool TEST_TENSOR_SAME(const std::string & test_name,
                                           const xsec::Array3D & arr1,
                                           const xsec::Array3D & arr2,
                                           const double & precision,
                                           const bool & verbose) {
    Eigen::ArrayXd _arr1 = xsec::ArrayMap(arr1.data(), arr1.size());
    Eigen::ArrayXd _arr2 = xsec::ArrayMap(arr2.data(), arr2.size());
    bool test = (_arr1 - _arr2).isZero(precision);
    if (!test || verbose) {
        std::cerr << test_name << (test ? ": PASSED" : ": FAILED") << std::endl;
        std::cerr << test_name << "\t" << _arr1.transpose() << std::endl;
        std::cerr << test_name << "\t" << _arr2.transpose() << std::endl;
    }
    return test;
}


[[nodiscard]] inline bool TEST_MULTIVERSE(const std::string & test_name,
                                          const xsec::Systematic<TH1> & mv1,
                                          const xsec::Systematic<TH1> & mv2,
                                          const double & precision,
                                          const bool & verbose) {
    bool test = true;
    for (auto imv = 0u; imv < (mv1).GetShifts().size(); imv++) {
        test &= (xsec::root::MapContentsToEigen(mv1.GetShifts()[imv]) -
                 xsec::root::MapContentsToEigen(mv2.GetShifts()[imv])).isZero(precision);
    }
    if (!test || verbose) {
        std::cerr << __FUNCTION__ << "\t" << test_name << (test ? ": PASSED" : ": FAILED") << std::endl;
        for (auto imv = 0u; imv < (mv1).GetShifts().size(); imv++) {
            std::cerr << __FUNCTION__ << "\t" << test_name << "[" << imv << "]\t"
                      << xsec::root::MapContentsToEigen(mv1.GetShifts()[imv]).transpose() << std::endl;
            std::cerr << __FUNCTION__ << "\t" << test_name << "[" << imv << "]\t"
                      << xsec::root::MapContentsToEigen(mv2.GetShifts()[imv]).transpose() << std::endl;
        }

    }
    return test;
}


namespace xsec {
    namespace test {
        namespace utils {
            inline TH1 * make_constant_hist_like(const TH1 * h, double c) {
                auto ret = (TH1 *) h->Clone();
                for (auto i = 0; i <= h->GetNbinsX() + 1; i++) {
                    for (auto j = 0; j <= h->GetNbinsX() + 1; j++) {
                        for (auto k = 0; k <= h->GetNbinsX() + 1; k++) {
                            ret->SetBinContent(i, j, k, c);
                        }
                    }
                }
                return ret;
            }

            inline TH1 * make_constant_hist(const Array & bins, double c) {
                auto one_a = Array::Ones(bins.size() + 1);
                auto ret = new TH1D("", "",
                                    bins.size() - 1,
                                    bins.data());
                ret->SetContent((one_a + c).eval().data());
                ret->SetEntries(ret->Integral());

                assert(std::sqrt(ret->GetBinContent(1)) == ret->GetBinError(1));
                return ret;
            }

            /////////////////////////////////////////////////////////
            inline void fill_random(TH1 * h, int n, int seed = 0) {
                std::srand(seed);
                if (h->GetDimension() == 1) {
                    auto minx = h->GetXaxis()->GetBinLowEdge(0);
                    auto maxx = h->GetXaxis()->GetBinLowEdge(h->GetNbinsX() + 2);
                    for (auto i = 0u; i < n; i++) {
                        h->Fill(((float) rand() / RAND_MAX) * (maxx - minx) + minx);
                    }
                } else if (h->GetDimension() == 2) {
                    auto h2 = dynamic_cast<TH2 *>(h);
                    auto minx = h->GetXaxis()->GetBinLowEdge(0);
                    auto maxx = h->GetXaxis()->GetBinLowEdge(h->GetNbinsX() + 2);
                    auto miny = h->GetYaxis()->GetBinLowEdge(0);
                    auto maxy = h->GetYaxis()->GetBinLowEdge(h->GetNbinsY() + 2);
                    for (auto i = 0u; i < n; i++) {
                        h2->Fill(((float) rand() / RAND_MAX) * (maxx - minx) + minx,
                                 ((float) rand() / RAND_MAX) * (maxy - miny) + miny);
                    }
                } else {
                    auto h3 = dynamic_cast<TH3 *>(h);
                    auto minx = h->GetXaxis()->GetBinLowEdge(0);
                    auto maxx = h->GetXaxis()->GetBinLowEdge(h->GetNbinsX() + 2);
                    auto miny = h->GetYaxis()->GetBinLowEdge(0);
                    auto maxy = h->GetYaxis()->GetBinLowEdge(h->GetNbinsY() + 2);
                    auto minz = h->GetZaxis()->GetBinLowEdge(0);
                    auto maxz = h->GetZaxis()->GetBinLowEdge(h->GetNbinsZ() + 2);
                    for (auto i = 0u; i < n; i++) {
                        h3->Fill(((float) rand() / RAND_MAX) * (maxx - minx) + minx,
                                 ((float) rand() / RAND_MAX) * (maxy - miny) + miny,
                                 ((float) rand() / RAND_MAX) * (maxz - minz) + minz);
                    }

                }

            }

            std::string test_dir() {
                return "./";
            }

            TH1 * make_simple_hist() {
                return new TH1D("", "",
                                10,
                                Array::LinSpaced(11, 0, 20).eval().data());
            }

            /////////////////////////////////////////////////////////
            const TH1 * get_simple_data() {
                auto ret = make_simple_hist();
                ret->SetContent((Array::LinSpaced(12, 1, 2) + 2).eval().data());
                ret->SetEntries(ret->Integral());
                return ret;
            }

            /////////////////////////////////////////////////////////
            const TH1 * get_hist_of_ones() {
                auto ret = new TH1D("", "",
                                    10,
                                    Array::LinSpaced(11, 0, 11).eval().data());
                ret->SetContent(Array::Ones(12).eval().data());
                ret->SetEntries(ret->Integral());
                return ret;
            }

            /////////////////////////////////////////////////////////
            const TH1 * get_simple_background() {
                auto ret = make_simple_hist();
                ret->SetContent(Array::Ones(12).eval().data());
                ret->SetEntries(ret->Integral());
                ret->Sumw2();
                return ret;
            }

            /////////////////////////////////////////////////////////
            const TH1 * get_simple_signal() {
                auto signal = (TH1 *) get_simple_data()->Clone();
                signal->Add(get_simple_background(), -1);
                return signal;
            }

            /////////////////////////////////////////////////////////
            const TH1 * get_simple_nominal_hist() {
                auto nominal = make_simple_hist();
                nominal->SetContent(Array::LinSpaced(12, 1, 20).eval().data());
                nominal->SetEntries(nominal->Integral());
                return nominal;
            }

            /////////////////////////////////////////////////////////
            const TH1 * get_simple_up_hist() {
                auto nominal = get_simple_nominal_hist();
                auto step = make_simple_hist();
                step->SetContent(Array::LinSpaced(12, 1.1, 1.3).pow(2).reverse().eval().data());
                auto ret = (TH1 *) nominal->Clone();
                ret->Multiply(step);
                ret->SetEntries(ret->Integral());
                return ret;
            }

            /////////////////////////////////////////////////////////
            const TH1 * get_simple_down_hist() {
                auto nominal = get_simple_nominal_hist();
                auto step = make_simple_hist();
                step->SetContent(Array::LinSpaced(12, .7, .9).pow(2).eval().data());
                auto ret = (TH1 *) nominal->Clone();
                ret->Multiply(step);
                ret->SetEntries(ret->Integral());
                return ret;
            }

            /////////////////////////////////////////////////////////
            std::vector<const TH1 *> make_simple_hist_multiverse(const TH1 * hnominal, int nuniverses) {
                double maxy = 0.1;
                double miny = -0.1;
                double step = (maxy - miny) / (nuniverses - 1);
                std::vector<const TH1 *> hist_universes(nuniverses);
                for (auto i = 0; i < nuniverses; i++) {
                    auto h = (TH1 *) hnominal->Clone();
                    h->Scale(miny + step * i);
                    h->Add(hnominal);
                    h->SetError(Array::Zero(hnominal->GetNbinsX()+2).eval().data());
                    hist_universes[i] = h;
                    h = 0;
                }
                return hist_universes;
            }
        } // utils
    } // test
} // xsec
