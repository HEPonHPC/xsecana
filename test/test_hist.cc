
#include <memory>
#include <iterator>

#include "XSecAna/Utils.h"

#include "test_utils.h"
#include "TFile.h"

using namespace xsec;

/*
 * Test closure of ROOT<->Eigen conversion functions
 */


inline void AreEqual(const TH1 * h1,
                     const TH1 * h2,
                     bool & same_content,
                     bool & same_error) {
    same_content &= (h1->GetDimension() == h2->GetDimension() &&
                     h1->GetNbinsX() == h2->GetNbinsX() &&
                     h1->GetNbinsY() == h2->GetNbinsY() &&
                     h1->GetNbinsZ() == h2->GetNbinsZ());
    same_error &= same_content;

    if (same_content) {
        for (auto i = 0u; i <= h1->GetNbinsX() + 1; i++) {
            for (auto j = 0u; j <= h1->GetNbinsY() + 1; j++) {
                for (auto k = 0u; k <= h1->GetNbinsZ() + 1; k++) {
                    same_content &= h1->GetBinContent(i, j, k) == h2->GetBinContent(i, j, k);
                    same_error &= h1->GetBinContent(i, j, k) == h2->GetBinContent(i, j, k);
                }
            }
        }
    }
}

int main(int argc, char ** argv) {
    bool verbose = false;
    if (argc > 1 && std::strcmp(argv[1], "-v") == 0) verbose = true;
    bool pass = true;

    auto nx = 10;
    auto ny = 3;
    auto nz = 4;

    auto th1 = new TH1D("", "",
                        nx, 0, nx);
    test::utils::fill_random(th1, 1000);

    auto th2 = new TH2D("", "",
                        nx,
                        0, nx,
                        ny,
                        0, ny);
    test::utils::fill_random(th2, 1000);

    auto th3 = new TH3D("", "",
                        nx, 0, nx,
                        ny, 0, ny,
                        nz, 0, nz);
    test::utils::fill_random(th3, 1000);

    Array arr1d_c = root::MapContentsToEigen(th1);
    Array arr1d_e = root::MapErrorsToEigen(th1);

    Array arr2d_c = root::MapContentsToEigen(th2);
    Array arr2d_e = root::MapErrorsToEigen(th2);

    Array arr3d_c = root::MapContentsToEigen(th3);
    Array arr3d_e = root::MapErrorsToEigen(th3);

    bool equal_content1 = true;
    bool equal_error1 = true;
    AreEqual(th1, root::ToROOTLike(th1, arr1d_c, arr1d_e), equal_content1, equal_error1);
    assert(equal_content1 && equal_error1);


    bool equal_content2 = true;
    bool equal_error2 = true;
    AreEqual(th2, root::ToROOTLike(th2, arr2d_c, arr2d_e), equal_content2, equal_error2);
    assert(equal_content2 && equal_error2);

    bool equal_content3 = true;
    bool equal_error3 = true;
    AreEqual(th3, root::ToROOTLike(th3, arr3d_c, arr3d_e), equal_content3, equal_error3);
    assert(equal_content3 && equal_error3);

    auto bin_width = test::utils::make_constant_hist_like(th1, 1);
    for (auto i = 1; i <= bin_width->GetNbinsX(); i++) {
        bin_width->SetBinContent(i, bin_width->GetBinWidth(i));
        bin_width->SetBinError(i, 0);
    }
    assert((root::MapContentsToEigen(bin_width) - root::MapBinWidthsToEigen(th1)).isZero(0));

    auto output = new TFile("test_hist.root", "recreate");
    th1->Write("th1");
    th2->Write("th2");
    th3->Write("th3");
    output->Close();
    delete output;

    auto input = TFile::Open("test_hist.root");
    auto lth1 = root::LoadTH1(input, "th1").release();
    auto lth2 = root::LoadTH1(input, "th2").release();
    auto lth3 = root::LoadTH1(input, "th3").release();
    input->Close();
    delete input;
    

    lth1->GetNbinsY();
    lth1->GetDimension();
    AreEqual(th1, lth1, equal_content1, equal_error1);
    assert(equal_content1 && equal_error1);

    AreEqual(th2, lth2, equal_content2, equal_error2);
    assert(equal_content2 && equal_error2);

    AreEqual(th3, lth3, equal_content3, equal_error3);
    assert(equal_content3 && equal_error3);

    return 0;
}