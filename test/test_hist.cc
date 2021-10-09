#include "XSecAna/Hist.h"
#include "XSecAna/Array.h"
#include "test_utils.h"

#include <Eigen/Dense>

#include <iostream>
#include <memory>
#include <iterator>

#include "TFile.h"

using namespace xsec;

int main(int argc, char ** argv) {
    bool verbose = false;
    if (argc > 1 && std::strcmp(argv[1], "-v") == 0) verbose = true;
    bool pass = true;

    auto nbins = 10;
    auto nbins_and_uof = nbins + 2;
    auto nedges = nbins + 1;
    auto nedges_and_uof = nbins_and_uof + 1;
    auto maxx = 10;
    auto minx = 0;
    auto step = (maxx - minx) / nbins;

    auto exposure1 = 1.;
    auto exposure2 = 3.;

    auto contents = Array::LinSpaced(nbins_and_uof, 1, 10) - 5.5;

    auto array = WeightedArray(contents,
                               exposure1);
    auto array2 = WeightedArray(contents,
                                exposure2);
    auto bins = Array::LinSpaced(nedges_and_uof,
                                 minx - step,
                                 maxx + step);
    auto ones = Array::Ones(nbins_and_uof);
    auto errors1 = ones * 3;
    auto errors2 = ones * 4;

    // total errors relative to hist1
    auto errors_total_rel1 = (errors1.pow(2) +
                              ((exposure1 / exposure2) * errors2).pow(2)).sqrt();


    auto hist = new Hist(array, bins, errors1);
    auto hist2 = new Hist(array2, bins, errors2);

    assert((hist->GetContentsAndUOF() - array.array()).isZero(0));
    assert(hist->GetContents().size() == nbins);
    assert(hist->GetEdges().size() == nbins + 1);
    assert(hist->GetErrors().size() == nbins);

    // make sure this is a deep copy
    auto hist_cpy = dynamic_cast<Hist *>(hist->Clone());

    assert(hist != hist_cpy);
    assert((hist->GetEdgesAndUOF() - hist_cpy->GetEdgesAndUOF()).isZero(0));
    assert((hist->GetContentsAndUOF() - hist_cpy->GetContentsAndUOF()).isZero(0));
    assert((hist->GetErrorsAndUOF() - hist_cpy->GetErrorsAndUOF()).isZero(0));
    assert(hist->Exposure() == hist_cpy->Exposure());
    assert(hist->GetEdgesAndUOF().data() != hist_cpy->GetEdgesAndUOF().data());
    assert(hist->GetContentsAndUOF().data() != hist_cpy->GetContentsAndUOF().data());
    assert(hist->GetErrorsAndUOF().data() != hist_cpy->GetErrorsAndUOF().data());

    double tol = 1e-14;

    bool test;

    pass &= TEST_HIST("construction",
                      hist,
                      array.array(),
                      bins,
                      errors1,
                      0,
                      verbose);

    TEST_ARRAY_SAME("array multiply", (array * array2).array(), (contents.pow(2) / 3), tol)
    TEST_ARRAY_SAME("array divide", (array / array2).array(), ones * 3, 0);
    TEST_ARRAY_SAME("array add", (array + array2).array(), contents * (4. / 3.), tol);
    TEST_ARRAY_SAME("array subtract", (array - array2).array(), contents * (2. / 3.), tol);
    TEST_ARRAY_SAME("array chained expression", ((array - array2) / array2).array(), ones * (2.), tol);
    TEST_ARRAY_SAME("array constant multiply", (array * -1).array(), contents * -1, 0);
    TEST_ARRAY_SAME("array abs()", array.abs().array(), contents.abs(), 0);
    TEST_ARRAY_SAME("array abs().sqrt()", array.abs().sqrt().array(), contents.abs().sqrt(), 0);
    TEST_ARRAY_SAME("array no change", array.array(), contents, 0);
    array = array.abs();
    TEST_ARRAY_SAME("modify array abs", array.array(), contents.abs(), 0);

    pass &= TEST_HIST("hist multiply",
                      dynamic_cast<Hist *>(hist->Multiply(hist2)),
                      contents.pow(2) / 3,
                      bins,
                      errors_total_rel1,
                      tol,
                      verbose);
    pass &= TEST_HIST("hist divide",
                      dynamic_cast<Hist *>(hist->Divide(hist2)),
                      ones * 3,
                      bins,
                      errors_total_rel1,
                      tol,
                      verbose);


    pass &= TEST_HIST("hist add",
                      dynamic_cast<Hist *>(hist->Add(hist2)),
                      contents * (4. / 3.),
                      bins,
                      errors_total_rel1,
                      tol,
                      verbose);


    pass &= TEST_HIST("hist subtract",
                      dynamic_cast<Hist *>(hist->Subtract(hist2)),
                      contents * (2. / 3.),
                      bins,
                      errors_total_rel1,
                      tol,
                      verbose);

    // although this isn't technically correct because the expression
    // can be simplified to prevent double counting of errors from hist2,
    // we don't do the simplification when calculating errors.
    // The user needs to be made aware of this behavior so they can
    // do their own simplification.
    auto errors_total_rel1_chained_expression = (errors1.pow(2) +
                                                 2 * ((exposure1 / exposure2) * errors2).pow(2)).sqrt();

    TEST_HIST("hist chained expression",
              dynamic_cast<Hist *>(hist->Subtract(hist2)->Divide(hist2)),
              ones * (2.),
              bins,
              errors_total_rel1_chained_expression,
              tol,
              verbose);

    TEST_HIST("hist constant multiply",
              dynamic_cast<Hist*>(hist->Multiply(-1)),
              contents * -1,
              bins,
              errors1,
              verbose,
              0);


    TEST_HIST("hist abs()",
              dynamic_cast<Hist*>(hist->abs()),
              contents.abs(),
              bins,
              errors1,
              0,
              verbose);


    TEST_HIST("hist abs().sqrt()",
              dynamic_cast<Hist*>(hist->abs()->sqrt()),
              contents.abs().sqrt(),
              bins,
              errors1/ 2,
              1e-6,
              verbose);

    TEST_HIST("no change",
              hist,
              contents,
              bins,
              errors1,
              0,
              verbose);



    hist->abs(true);
    TEST_HIST("inplace hist abs",
              hist,
              contents.abs(),
              bins,
              errors1,
              0,
              verbose);

    TEST_ARRAY_SAME("bin width",
                    hist->GetBinWidths(),
                    ones(Eigen::seqN(0, ones.size() - 2)),
                    0);

    TEST_ARRAY_SAME("contiguous", hist->GetContents(), test::utils::is_contiguous(hist->GetContents()), 0);

    // saveto/loadfrom
    std::string test_file_name = test::utils::test_dir() + "test_hist->root";
    auto output = new TFile(test_file_name.c_str(), "recreate");
    hist->SaveTo(output, "hist");
    output->Close();

    auto input = TFile::Open(test_file_name.c_str());
    auto loaded = dynamic_cast<Hist*>(Hist::LoadFrom(input, "hist").release());

    TEST_HIST("saveto/loadfrom",
              loaded,
              hist->GetContentsAndUOF(),
              hist->GetEdgesAndUOF(),
              hist->GetErrorsAndUOF(),
              0,
              verbose);

    return !pass;
}