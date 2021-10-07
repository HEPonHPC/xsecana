#include "XSecAna/Hist.h"
#include "XSecAna/SimpleEfficiency.h"
#include "test_utils.h"

#include <iostream>

#include "TFile.h"

using namespace xsec;

int main(int argc, char ** argv) {
    bool verbose = false;
    if (argc > 1 && std::strcmp(argv[1], "-v") == 0) verbose = true;
    bool pass = true;
    bool test;

    Hist num = test::utils::get_hist_of_ones();
    Hist den = test::utils::get_hist_of_ones() / 2;


    SimpleEfficiency eff(num, den);

    TEST_HISTS_SAME("ratio calculation",
                    eff.Eval(),
                    (num / den),
                    0);

    std::string test_file_name = test::utils::test_dir() + "test_simple_efficiency.root";
    auto output = new TFile(test_file_name.c_str(), "recreate");
    eff.SaveTo(output, "simple_efficiency");
    output->Close();
    delete output;

    TFile * input = TFile::Open(test_file_name.c_str());
    auto loaded = IEfficiency<Hist>::LoadFrom(SimpleEfficiency<Hist>::LoadFrom,
                                              input,
                                              "simple_efficiency").release();
    input->Close();
    delete input;

    TEST_HISTS_SAME("saveto/loadfrom numerator",
                    ((SimpleEfficiency<Hist> *) loaded)->GetNumerator(),
                    num,
                    0);
    TEST_HISTS_SAME("saveto/loadfrom denominator",
                    ((SimpleEfficiency<Hist> *) loaded)->GetDenominator(),
                    den,
                    0);
    TEST_HISTS_SAME("saveto/loadfrom ratio",
                    loaded->Eval(),
                    eff.Eval(),
                    0);

    return !pass;
}
