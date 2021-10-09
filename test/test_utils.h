#pragma once

#include <iostream>
#include "XSecAna/Hist.h"
#include "XSecAna/IUnfold.h"
#include "XSecAna/SimpleSignalEstimator.h"
#include "XSecAna/SimpleEfficiency.h"
#include "XSecAna/SimpleFlux.h"
#include "XSecAna/CrossSection.h"
#include "XSecAna/Array3D.h"


//----------------------------------------------------------------------
bool TEST_HIST(std::string test_name,
               const xsec::Hist * HIST,
               const xsec::Array & target_contents,
               const xsec::Array & target_edges,
               const xsec::Array & target_errors,
               double precision,
               bool verbose) {
    bool pass = true;
    bool test = (HIST->GetContentsAndUOF() - target_contents).isZero(precision);
    if (!test || verbose) {
        std::cerr << __FUNCTION__ << "\t" << test_name + " (contents)" << (test ? ": PASSED" : ": FAILED")
                  << std::endl;
        std::cerr << __FUNCTION__ << "\t" << test_name + " (contents)" << "\t"
                  << HIST->GetContentsAndUOF().transpose() << std::endl;
        std::cerr << __FUNCTION__ << "\t" << test_name + " (contents)" << "\t" << target_contents.transpose()
                  << std::endl;
        pass &= test;
    }
    test = (HIST->GetEdgesAndUOF() - target_edges).isZero(precision);
    if (!test || verbose) {
        std::cerr << __FUNCTION__ << "\t" << test_name + " (edges)" << (test ? ": PASSED" : ": FAILED")
                  << std::endl;
        std::cerr << __FUNCTION__ << "\t" << test_name + " (edges)" << "\t" << HIST->GetEdgesAndUOF().transpose()
                  << std::endl;
        std::cerr << __FUNCTION__ << "\t" << test_name + " (edges)" << "\t" << target_edges.transpose()
                  << std::endl;
        pass &= test;
    }
    test = (HIST->GetErrorsAndUOF() - target_errors).isZero(precision);
    if (!test || verbose) {
        std::cerr << __FUNCTION__ << "\t" << test_name + " (errors)" << (test ? ": PASSED" : ": FAILED")
                  << std::endl;
        std::cerr << __FUNCTION__ << "\t" << test_name + " (errors)" << "\t" << HIST->GetErrorsAndUOF().transpose()
                  << std::endl;
        std::cerr << __FUNCTION__ << "\t" << test_name + " (errors)" << "\t" << target_errors.transpose()
                  << std::endl;
        pass &= test;
    }
    return pass;
}

//----------------------------------------------------------------------
#define TEST_HISTS_SAME(test_name, h1, h2, precision)            \
  test = ((h1) - (h2)).GetGetContentsAndUOF().isZero(precision);                \
      if(!test || verbose) {                        \
    std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << (test? ": PASSED" : ": FAILED") << std::endl; \
    std::cerr << __PRETTY_FUNCTION__ << "\t (" << (h1).Exposure() << ") " << test_name << "\t" << (h1).GetContentsAndUOF().transpose() << std::endl; \
    std::cerr << __PRETTY_FUNCTION__ << "\t (" << (h2).Exposure() << ") " << test_name << "\t" << (h2).GetContentsAndUOF().transpose() << std::endl; \
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

void TEST_TENSOR_SAME(std::string test_name,
                      const xsec::Array3D & arr1,
                      const xsec::Array3D & arr2,
                      double precision,
                      bool & pass,
                      bool verbose) {
    Eigen::ArrayXd _arr1 = xsec::ToArray((arr1));
    Eigen::ArrayXd _arr2 = xsec::ToArray((arr2));
    bool test = (_arr1 - _arr2).isZero(precision);
    if (!test || verbose) {
        std::cerr << test_name << (test ? ": PASSED" : ": FAILED") << std::endl;
        std::cerr << test_name << "\t" << _arr1.transpose() << std::endl;
        std::cerr << test_name << "\t" << _arr2.transpose() << std::endl;
        pass &= test;
    }
}


//----------------------------------------------------------------------
#define TEST_MULTIVERSE(test_name, mv1, mv2, precision)                \
  test = true;                                \
  for(auto imv = 0u; imv < (mv1).GetShifts().size(); imv++) {        \
    test &= ((mv1).GetShifts()[imv]->GetContents() - (mv2).GetShifts()[imv]->GetContents()).isZero(precision);        \
  }                                    \
  if(!test || verbose) {                        \
    std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << (test? ": PASSED" : ": FAILED") << std::endl; \
    for(auto imv = 0u; imv < (mv1).GetShifts().size(); imv++) {    \
      std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << "[" << imv << "]\t" << (mv1).GetShifts()[imv]->GetContents().transpose() << std::endl; \
      std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << "[" << imv << "]\t" << (mv2).GetShifts()[imv]->GetContents().transpose() << std::endl; \
    }                                    \
    pass &= test;                            \
  }


namespace xsec {
    namespace test {
        namespace utils {
            /////////////////////////////////////////////////////////
            struct Ratio {
                Ratio() {}

                Ratio(_hist * num,
                      _hist * den)
                        : numerator(num), denominator(den) {}

                const _hist * numerator;
                const _hist * denominator;

                _hist * Eval() const {
                    return numerator->Divide(denominator);
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
            const Hist * get_simple_data() {
                return new Hist(Array::LinSpaced(12, 1, 2) + 2,
                                Array::LinSpaced(13, 0, 20),
                                data_exposure); // data will often come at different exposure from mc
            }

            /////////////////////////////////////////////////////////
            const Hist * get_hist_of_ones() {
                return new Hist(Array::Ones(12),
                                Array::LinSpaced(13, 0, 13));
            }

            /////////////////////////////////////////////////////////
            const Hist * get_simple_background() {
                return new Hist(Array::Ones(12),
                                Array::LinSpaced(13, 0, 20));
            }

            /////////////////////////////////////////////////////////
            const Hist * get_simple_signal() {
                return dynamic_cast<Hist *>(get_simple_data()->Subtract(get_simple_background()));
            }

            /////////////////////////////////////////////////////////
            const Hist * get_simple_nominal_hist() {
                return new Hist(Array::LinSpaced(12, -0.5, 0.5) + 4,
                                Array::LinSpaced(13, 0, 20));
            }

            /////////////////////////////////////////////////////////
            const Hist * get_simple_up_hist() {
                auto nominal = get_simple_nominal_hist();
                auto step = Array::LinSpaced(12, 0, .3).reverse();

                return new Hist(nominal->GetContentsAndUOF() +
                                nominal->GetContentsAndUOF() * step.pow(2),
                                nominal->GetEdgesAndUOF());
            }

            /////////////////////////////////////////////////////////
            const Hist * get_simple_down_hist() {
                auto nominal = get_simple_nominal_hist();
                auto step = Array::LinSpaced(12, 0, .3);

                return new Hist(nominal->GetContentsAndUOF() -
                                nominal->GetContentsAndUOF() * step.pow(2),
                                nominal->GetEdgesAndUOF());
            }

            const static double ntargets = 1e4;

            /////////////////////////////////////////////////////////
            // make a cross section object that evaluates out to
            //  the input array when folded
            IMeasurement * make_simple_xsec(const Hist * val) {
                auto ones = new Hist(Array::Ones(12),
                                     val->GetEdgesAndUOF(),
                                     data_exposure);

                auto efficiency = new SimpleEfficiency(get_simple_signal()->Divide(val),
                                                       ones);
                auto flux = new SimpleFlux(new Hist(ones->GetContentsAndUOF(),
                                                    ones->GetEdgesAndUOF(),
                                                    get_simple_data()->Exposure()));
                auto signal_estimator = new SimpleSignalEstimator(get_simple_background());
                auto unfold = new IdentityUnfold(ones->GetContentsAndUOF().size());
                auto ret = new CrossSection(efficiency,
                                            signal_estimator,
                                            flux,
                                            unfold,
                                            ntargets);
                return ret;
            }


            /////////////////////////////////////////////////////////
            std::vector<Hist*> make_simple_hist_multiverse(const Hist * hnominal, int nuniverses) {
                double maxy = 0.1;
                double miny = -0.1;
                double step = (maxy - miny) / (nuniverses - 1);
                std::vector<Hist *> hist_universes(nuniverses);
                for (auto i = 0; i < nuniverses; i++) {
                    hist_universes[i] = dynamic_cast<Hist*>(hnominal->Add(hnominal->Multiply(miny + step * i)));
                }
                return hist_universes;
            }

            /////////////////////////////////////////////////////////
            std::vector<IMeasurement *>
            make_simple_xsec_multiverse(const Hist * hnominal, int nuniverses) {
                std::vector<IMeasurement *> xsec_universes(nuniverses);
                auto hist_universes = make_simple_hist_multiverse(hnominal, nuniverses);

                for (auto i = 0; i < nuniverses; i++) {
                    xsec_universes[i] = make_simple_xsec(hist_universes[i]);
                }
                return xsec_universes;
            }

        } // utils
    } // test
} // xsec
