#pragma once

#include <iostream>
#include "XSecAna/Hist.h"
#include "XSecAna/IUnfold.h"
#include "XSecAna/SimpleSignalEstimator.h"
#include "XSecAna/SimpleEfficiency.h"
#include "XSecAna/SimpleFlux.h"
#include "XSecAna/CrossSection.h"


//----------------------------------------------------------------------
#define TEST_HIST_AND_EDGES(test_name, HIST, target_contents, target_edges, precision) \
  test = ((HIST).ContentsAndUOF() - (target_contents)).isZero(precision);    \
      if(!test || verbose) {                        \
    std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << (test? ": PASSED" : ": FAILED") << std::endl; \
    std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << "\t" << (HIST).ContentsAndUOF().transpose() << std::endl; \
    std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << "\t" << (target_contents).transpose() << std::endl; \
    pass &= test;                            \
      }                                    \
      test = ((HIST).EdgesAndUOF() - (target_edges)).isZero(precision);    \
      if(!test || verbose) {                        \
    std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << (test? ": PASSED" : ": FAILED") << std::endl; \
    std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << "\t" << (HIST).EdgesAndUOF().transpose() << std::endl; \
    std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << "\t" << (target_edges).transpose() << std::endl; \
    pass &= test;                            \
      }

//----------------------------------------------------------------------
#define TEST_HISTS_SAME(test_name, h1, h2, precision)            \
  test = ((h1) - (h2)).ContentsAndUOF().isZero(precision);                \
      if(!test || verbose) {                        \
    std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << (test? ": PASSED" : ": FAILED") << std::endl; \
    std::cerr << __PRETTY_FUNCTION__ << "\t (" << (h1).Exposure() << ") " << test_name << "\t" << (h1).ContentsAndUOF().transpose() << std::endl; \
    std::cerr << __PRETTY_FUNCTION__ << "\t (" << (h2).Exposure() << ") " << test_name << "\t" << (h2).ContentsAndUOF().transpose() << std::endl; \
    pass &= test;                            \
      }

#define TEST_ARRAY_SAME(test_name, arr1, arr2, precision)            \
  test = ((arr1) - (arr2)).isZero(precision);                \
      if(!test || verbose) {                        \
    std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << (test? ": PASSED" : ": FAILED") << std::endl; \
    std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << "\t" << (arr1).transpose() << std::endl; \
    std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << "\t" << (arr2).transpose() << std::endl; \
    pass &= test;                            \
      }


//----------------------------------------------------------------------
#define TEST_MULTIVERSE(test_name, mv1, mv2, precision)                \
  test = true;                                \
  for(auto imv = 0u; imv < (mv1).GetShifts().size(); imv++) {        \
    test &= ((mv1).GetShifts()[imv]->Contents() - (mv2).GetShifts()[imv]->Contents()).isZero(precision);        \
  }                                    \
  if(!test || verbose) {                        \
    std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << (test? ": PASSED" : ": FAILED") << std::endl; \
    for(auto imv = 0u; imv < (mv1).GetShifts().size(); imv++) {    \
      std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << "[" << imv << "]\t" << (mv1).GetShifts()[imv]->Contents().transpose() << std::endl; \
      std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << "[" << imv << "]\t" << (mv2).GetShifts()[imv]->Contents().transpose() << std::endl; \
    }                                    \
    pass &= test;                            \
  }


namespace xsec {
    namespace test {
        namespace utils {
            /////////////////////////////////////////////////////////
            struct Ratio {
                Ratio() {}

                Ratio(Hist num,
                      Hist den)
                        : numerator(num), denominator(den) {}

                Hist numerator;
                Hist denominator;

                Hist Eval() const {
                    return numerator / denominator;
                }
            };

            /////////////////////////////////////////////////////////
            template<typename ArrayType>
            ArrayType is_contiguous(const ArrayType & array) {
                ArrayType a2 = array;
                auto const data = array.data();
                for (auto i = 0; i < array.size(); ++i) {
                    a2(i) = *(&data[0] + i);
                }
                return a2;
            }

            /////////////////////////////////////////////////////////
            std::string test_dir() {
                return "./";
            }

            const static double data_exposure = 0.5;

            /////////////////////////////////////////////////////////
            Hist get_simple_data() {
                return Hist(Array::LinSpaced(12, 1, 2) + 2,
                            Array::LinSpaced(13, 0, 20),
                            data_exposure); // data will often come at different exposure from mc
            }

            /////////////////////////////////////////////////////////
            Hist get_hist_of_ones() {
                return Hist(Array::Ones(12),
                            Array::LinSpaced(13, 0, 13));
            }

            /////////////////////////////////////////////////////////
            Hist get_simple_background() {
                return Hist(Array::Ones(12),
                            Array::LinSpaced(13, 0, 20));
            }

            /////////////////////////////////////////////////////////
            Hist get_simple_signal() {
                return get_simple_data() - get_simple_background();
            }

            /////////////////////////////////////////////////////////
            Hist get_simple_nominal_hist() {
                return Hist(Array::LinSpaced(12, -0.5, 0.5) + 4,
                            Array::LinSpaced(13, 0, 20));
            }

            /////////////////////////////////////////////////////////
            Hist get_simple_up_hist() {
                auto nominal = get_simple_nominal_hist();
                auto step = Array::LinSpaced(12, 0, .3).reverse();

                return Hist(nominal.ContentsAndUOF() +
                            nominal.ContentsAndUOF() * step.pow(2),
                            nominal.EdgesAndUOF());
            }

            /////////////////////////////////////////////////////////
            Hist get_simple_down_hist() {
                auto nominal = get_simple_nominal_hist();
                auto step = Array::LinSpaced(12, 0, .3);

                return Hist(nominal.ContentsAndUOF() -
                            nominal.ContentsAndUOF() * step.pow(2),
                            nominal.EdgesAndUOF());
            }

            const static double ntargets = 1e4;

            /////////////////////////////////////////////////////////
            // make a cross section object that evaluates out to
            //  the input array when folded
            IMeasurement * make_simple_xsec(Hist val) {
                Hist ones(Array::Ones(12),
                          val.EdgesAndUOF(),
                          data_exposure);

                auto efficiency = new SimpleEfficiency(get_simple_signal() / val,
                                                       ones);
                auto flux = new SimpleFlux(Hist(ones.ContentsAndUOF(),
                                                ones.EdgesAndUOF(),
                                                get_simple_data().Exposure()));
                auto signal_estimator = new SimpleSignalEstimator(get_simple_background());
                auto unfold = new IdentityUnfold(ones.ContentsAndUOF().size());
                auto ret = new CrossSection(efficiency,
                                            signal_estimator,
                                            flux,
                                            unfold,
                                            ntargets);
                return ret;
            }


            /////////////////////////////////////////////////////////
            std::vector<Hist *> make_simple_hist_multiverse(const Hist & hnominal, int nuniverses) {
                double maxy = 0.1;
                double miny = -0.1;
                double step = (maxy - miny) / (nuniverses - 1);
                std::vector<Hist *> hist_universes(nuniverses);
                for (auto i = 0; i < nuniverses; i++) {
                    hist_universes[i] = new Hist(hnominal + hnominal * (miny + step * i));
                }
                return hist_universes;
            }

            /////////////////////////////////////////////////////////
            std::vector<IMeasurement *>
            make_simple_xsec_multiverse(const Hist & hnominal, int nuniverses) {
                std::vector<IMeasurement *> xsec_universes(nuniverses);
                auto hist_universes = make_simple_hist_multiverse(hnominal, nuniverses);

                for (auto i = 0; i < nuniverses; i++) {
                    xsec_universes[i] = make_simple_xsec(*hist_universes[i]);
                }
                return xsec_universes;
            }

        } // utils
    } // test
} // xsec
