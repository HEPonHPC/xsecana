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
  test = ((HIST).Contents() - (target_contents)).isZero(precision);    \
      if(!test || verbose) {                        \
    std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << (test? ": PASSED" : ": FAILED") << std::endl; \
    std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << "\t" << (HIST).Contents() << std::endl; \
    std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << "\t" << (target_contents) << std::endl; \
    pass &= test;                            \
      }                                    \
      test = ((HIST).Edges() - (target_edges)).isZero(precision);    \
      if(!test || verbose) {                        \
    std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << (test? ": PASSED" : ": FAILED") << std::endl; \
    std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << "\t" << (HIST).Edges() << std::endl; \
    std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << "\t" << (target_edges) << std::endl; \
    pass &= test;                            \
      }

//----------------------------------------------------------------------
#define TEST_HISTS_SAME(test_name, h1, h2, precision)            \
  test = ((h1) - (h2)).Contents().isZero(precision);                \
      if(!test || verbose) {                        \
    std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << (test? ": PASSED" : ": FAILED") << std::endl; \
    std::cerr << __PRETTY_FUNCTION__ << "\t (" << (h1).Exposure() << ") " << test_name << "\t" << (h1).Contents() << std::endl; \
    std::cerr << __PRETTY_FUNCTION__ << "\t (" << (h2).Exposure() << ") " << test_name << "\t" << (h2).Contents() << std::endl; \
    pass &= test;                            \
      }

#define TEST_ARRAY_SAME(test_name, arr1, arr2, precision)            \
  test = ((arr1) - (arr2)).isZero(precision);                \
      if(!test || verbose) {                        \
    std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << (test? ": PASSED" : ": FAILED") << std::endl; \
    std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << "\t" << (arr1) << std::endl; \
    std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << "\t" << (arr2) << std::endl; \
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
      std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << "[" << imv << "]\t" << (mv1).GetShifts()[imv]->Contents() << std::endl; \
      std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << "[" << imv << "]\t" << (mv2).GetShifts()[imv]->Contents() << std::endl; \
    }                                    \
    pass &= test;                            \
  }


namespace xsec {
    namespace test {
        namespace utils {
            /////////////////////////////////////////////////////////
            template<class Scalar,
                    int Cols>
            struct Ratio {
                Ratio() {}

                Ratio(Hist<Scalar, Cols> num,
                      Hist<Scalar, Cols> den)
                        : numerator(num), denominator(den) {}

                Hist<Scalar, Cols> numerator;
                Hist<Scalar, Cols> denominator;

                Hist<Scalar, Cols> Eval() const {
                    return numerator / denominator;
                }
            };

            /////////////////////////////////////////////////////////
            typedef Hist<double, 10> histtype;
            typedef CrossSection<histtype> SimpleCrossSection;

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
            template<class Scalar, int Cols>
            Hist<Scalar, Cols> get_simple_data() {
                return Hist<Scalar, Cols>(Eigen::Array<Scalar, 1, Cols>::LinSpaced(10, 1, 2) + 2,
                                          Eigen::Array<Scalar, 1, EdgesSize(Cols)>::LinSpaced(11, 0, 20),
                                          data_exposure); // data will often come at different exposure from mc
            }

            /////////////////////////////////////////////////////////
            template<class Scalar, int Cols>
            Hist<Scalar, Cols> get_simple_background() {
                return Hist<Scalar, Cols>(Eigen::Array<Scalar, 1, Cols>::Ones(),
                                          Eigen::Array<Scalar, 1, EdgesSize(Cols)>::LinSpaced(11, 0, 20));
            }

            /////////////////////////////////////////////////////////
            template<class Scalar, int Cols>
            Hist<Scalar, Cols> get_simple_signal() {
                return get_simple_data<Scalar, Cols>() - get_simple_background<Scalar, Cols>();
            }

            /////////////////////////////////////////////////////////
            template<class Scalar, int Cols>
            Hist<Scalar, Cols> get_simple_nominal_hist() {
                return Hist<Scalar, Cols>(Eigen::Array<Scalar, 1, Cols>::LinSpaced(10, -0.5, 0.5) + 4,
                                          Eigen::Array<Scalar, 1, EdgesSize(Cols)>::LinSpaced(11, 0, 20));
            }

            /////////////////////////////////////////////////////////
            template<class Scalar, int Cols>
            Hist<Scalar, Cols> get_simple_up_hist() {
                auto nominal = get_simple_nominal_hist<Scalar, Cols>();
                auto step = Eigen::Array<Scalar, 1, Cols>::LinSpaced(10, 0, .3).reverse();

                return Hist<Scalar, Cols>(nominal.Contents() + nominal.Contents() * step.pow(2),
                                          nominal.Edges());
            }

            /////////////////////////////////////////////////////////
            template<class Scalar, int Cols>
            Hist<Scalar, Cols> get_simple_down_hist() {
                auto nominal = get_simple_nominal_hist<Scalar, Cols>();
                auto step = Eigen::Array<Scalar, 1, Cols>::LinSpaced(10, 0, .3);

                return Hist<Scalar, Cols>(nominal.Contents() - nominal.Contents() * step.pow(2),
                                          nominal.Edges());
            }

            const static double ntargets = 1e4;

            /////////////////////////////////////////////////////////
            // make a cross section object that evaluates out to
            //  the input array when folded
            template<class Scalar, int Cols>
            IMeasurement<Hist<Scalar, Cols> > * make_simple_xsec(Hist<Scalar, Cols> val) {
                Hist<Scalar, Cols> ones(Eigen::Array<Scalar, 1, Cols>::Ones(),
                                        val.Edges());

                auto efficiency = new SimpleEfficiency<Hist<Scalar, Cols> >(get_simple_signal<Scalar, Cols>(),
                                                                            val * get_simple_data<Scalar, Cols>().Exposure() / val.Exposure());
                auto flux = new SimpleFlux(Hist<Scalar, Cols>(ones.Contents(),
                                                              ones.Edges(),
                                                              get_simple_data<Scalar, Cols>().Exposure()));
                auto signal_estimator = new SimpleSignalEstimator(get_simple_background<Scalar, Cols>());
                auto unfold = new IdentityUnfold<Scalar, Cols>(ones.Contents().size());
                auto ret = new SimpleCrossSection(efficiency,
                                              signal_estimator,
                                              flux,
                                              unfold);
                ret->SetNTargets(ntargets);
                return ret;
            }


            /////////////////////////////////////////////////////////
            template<class HistType>
            std::vector<HistType*> make_simple_hist_multiverse(const HistType & hnominal, int nuniverses) {
                double maxy = 0.1;
                double miny = -0.1;
                double step = (maxy - miny) / (nuniverses - 1);
                std::vector<HistType*> hist_universes(nuniverses);
                for (auto i = 0; i < nuniverses; i++) {
                    hist_universes[i] = new HistType(hnominal + hnominal * (miny + step * i));
                }
                return hist_universes;
            }

            /////////////////////////////////////////////////////////
            template<class HistType>
            std::vector<IMeasurement<HistType>*> make_simple_xsec_multiverse(const HistType & hnominal, int nuniverses) {
                std::vector<IMeasurement<HistType>*> xsec_universes(nuniverses);
                auto hist_universes = make_simple_hist_multiverse(hnominal, nuniverses);

                for (auto i = 0; i < nuniverses; i++) {
                    xsec_universes[i] = make_simple_xsec(*hist_universes[i]);
                }
                return xsec_universes;
            }


        } // utils
    } // test
} // xsec
