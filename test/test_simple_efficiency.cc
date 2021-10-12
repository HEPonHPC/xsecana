
#include "XSecAna/SimpleEfficiency.h"
#include "XSecAna/Utils.h"

//#include "test_utils.h"

#include <iostream>

#include "TFile.h"


using namespace xsec;

bool TEST_HIST(std::string test_name,
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

int main(int argc, char ** argv) {
    bool verbose = false;
    if (argc > 1 && std::strcmp(argv[1], "-v") == 0) verbose = true;
    bool pass = true;
    bool test;

    auto num = new TH1D("", "", 10, 0, 1);
    num->SetContent(Eigen::ArrayXd::Ones(12).eval().data());
    auto den = new TH1D("", "", 10, 0, 1);
    den->SetContent((10 * Eigen::ArrayXd::Ones(12)).eval().data());


    SimpleEfficiency eff(num, den);

    auto expected = (TH1 *) num->Clone();
    expected->Divide(num, den, 1, 1, "B");

    TEST_HIST("ratio calculation",
              eff.Eval(num), // pass num as dummy
              expected,
              0,
              verbose);

    std::string test_file_name = "test_simple_efficiency.root";
    auto output = new TFile(test_file_name.c_str(), "recreate");
    eff.SaveTo(output, "simple_efficiency");
    output->Close();
    delete output;

    TFile * input = TFile::Open(test_file_name.c_str());
    auto loaded = IMeasurement::LoadFrom(SimpleEfficiency::LoadFrom,
                                         input,
                                         "simple_efficiency").release();
    input->Close();
    delete input;

    TEST_HIST("saveto/loadfrom numerator",
              dynamic_cast<SimpleEfficiency*>(loaded)->GetNumerator(),
              num,
              0,
              verbose);
    TEST_HIST("saveto/loadfrom denominator",
              dynamic_cast<SimpleEfficiency*>(loaded)->GetDenominator(),
              den,
              0,
              verbose);
    TEST_HIST("saveto/loadfrom ratio",
              loaded->Eval(num),
              eff.Eval(num),
              0,
              verbose);


    return !pass;
}
