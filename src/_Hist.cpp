#include "XSecAna/_Hist.h"
#include "TParameter.h"
#include "XSecAna/Hist.h"

namespace xsec {
    /////////////////////////////////////////////////////////
    void
    _hist::
    SaveTo(TDirectory * dir, const std::string & subdir) const {
        TDirectory * tmp = gDirectory;
        dir->mkdir(subdir.c_str());
        dir = dir->GetDirectory(subdir.c_str());
        dir->cd();

        auto h = this->ToTH1();
        h->Write("hist");

        // exposure is saved in a histogram so it gets accumulated with ROOT's hadd
        TH1 * exposure = new TH1D("", "", 1, 0, 1);
        TParameter<int>("dims", this->GetDimensions()).Write("dims");

        exposure->SetBinContent(1, this->Exposure());
        exposure->Write("exposure");

        tmp->cd();
    }

    std::unique_ptr<_hist>
    _hist::
    LoadFrom(TDirectory * dir, const std::string & subdir) {
        dir = dir->GetDirectory(subdir.c_str());
        dir->cd();

        TH1 * h = (TH1 *) dir->Get("hist");
        if (!h) {
            std::cerr << "Object TH1 was not found in " << dir->GetPath() << std::endl;
            exit(1);
        }

        TH1 * exposure = (TH1D *) dir->Get("exposure");
        auto dims = ((TParameter<int> *) dir->Get("dims"))->GetVal();

        if(dims == 1) {
            return std::unique_ptr<_hist>(Hist::FromTH1(h, exposure->GetBinContent(1)));
        }
        else if(dims == 2) {
            return 0;
        }
        else if(dims == 3) {
            return 0;
        }
        else {
            return 0;
        }


    }

    void
    _hist::
    EnsureConsistentBinning(const _hist * rhs, const char * caller, double tol) const {
        if(typeid(this) != typeid(rhs) || !this->_is_same_binning(rhs, tol)) {
            throw exceptions::InconsistentBinningError(__PRETTY_FUNCTION__, tol);
        }
    }

    _hist *
    _hist::
    TrueDivide(const _hist * rhs) const {
        auto tmp = this->Clone();
        tmp->_true_divide(rhs);
        return tmp;
    }

    _hist *
    _hist::
    TrueDivide(const _hist * rhs, bool inplace) {
        _hist * tmp;
        if(inplace) tmp = this;
        else tmp = this->Clone();
        tmp->_true_divide(rhs);
        return tmp;
    }

    _hist *
    _hist::
    BinWidthNormalize() const {
        auto tmp = this->Clone();
        tmp->_bin_width_normalize();
        return tmp;
    }

    _hist *
    _hist::
    BinWidthNormalize(bool inplace) {
        _hist * tmp;
        if(inplace) tmp = this;
        else tmp = this->Clone();
        tmp->_bin_width_normalize();
        return tmp;
    }

    _hist *
    _hist::
    AreaNormalize() const {
        auto tmp = this->Clone();
        tmp->_area_normalize();
        return tmp;
    }

    _hist *
    _hist::
    AreaNormalize(bool inplace) {
        _hist * tmp;
        if(inplace) tmp = this;
        else tmp = this->Clone();
        tmp->_area_normalize();
        return tmp;
    }

    _hist *
    _hist::
    ScaleByExposure(const double & new_expo) const {
        auto tmp = this->Clone();
        tmp->_scale_by_exposure(new_expo);
        return tmp;
    }

    _hist *
    _hist::
    ScaleByExposure(const double & new_expo, bool inplace) {
        _hist * tmp;
        if(inplace) tmp = this;
        else tmp = this->Clone();
        tmp->_scale_by_exposure(new_expo);
        return tmp;
    }

    bool
    _hist::
    IsEqual(const _hist * rhs, const double & tol) const {
        return typeid(this) == typeid(rhs) &&
               this->_is_same_binning(rhs, tol) &&
               this->_is_same_contents(rhs, tol);
    }


    _hist *
    _hist::
    abs() const {
        auto tmp = this->Clone();
        tmp->_abs();
        return tmp;
    }

    _hist *
    _hist::
    abs(bool inplace) {
        _hist * tmp;
        if(inplace) tmp = this;
        else tmp = this->Clone();
        tmp->_abs();
        return tmp;
    }

    _hist *
    _hist::
    abs2() const {
        auto tmp = this->Clone();
        tmp->_abs2();
        return tmp;
    }

    _hist *
    _hist::
    abs2(bool inplace) {
        _hist * tmp;
        if(inplace) tmp = this;
        else tmp = this->Clone();
        tmp->_abs2();
        return tmp;
    }

    _hist *
    _hist::
    sqrt() const {
        auto tmp = this->Clone();
        tmp->_sqrt();
        return tmp;
    }

    _hist *
    _hist::
    sqrt(bool inplace) {
        _hist * tmp;
        if(inplace) tmp = this;
        else tmp = this->Clone();
        tmp->_sqrt();
        return tmp;
    }

    _hist *
    _hist::
    pow(double exp) const {
        auto tmp = this->Clone();
        tmp->_pow(exp);
        return tmp;
    }


    _hist *
    _hist::
    pow(double exp, bool inplace) {
        _hist * tmp;
        if(inplace) tmp = this;
        else tmp = this->Clone();
        tmp->_pow(exp);
        return tmp;
    }

    _hist *
    _hist::
    Subtract(const _hist * rhs, bool inplace) {
        this->EnsureConsistentBinning(rhs, __PRETTY_FUNCTION__);
        _hist * tmp;
        if(inplace) tmp = this;
        else tmp = this->Clone();
        tmp->_subtract(rhs);
        return tmp;
    }

    _hist *
    _hist::
    Subtract(const double & rhs, bool inplace) {
        _hist * tmp;
        if(inplace) tmp = this;
        else tmp = this->Clone();
        tmp->_subtract(rhs);
        return tmp;
    }

    _hist *
    _hist::
    Add(const _hist * rhs, bool inplace) {
        this->EnsureConsistentBinning(rhs, __PRETTY_FUNCTION__);
        _hist * tmp;
        if(inplace) tmp = this;
        else tmp = this->Clone();
        tmp->_add(rhs);
        return tmp;
    }

    _hist *
    _hist::
    Add(const double & rhs, bool inplace) {
        _hist * tmp;
        if(inplace) tmp = this;
        else tmp = this->Clone();
        tmp->_add(rhs);
        return tmp;
    }

    _hist *
    _hist::
    Divide(const _hist * rhs, bool inplace) {
        this->EnsureConsistentBinning(rhs, __PRETTY_FUNCTION__);
        _hist * tmp;
        if(inplace) tmp = this;
        else tmp = this->Clone();
        tmp->_divide(rhs);
        return tmp;
    }

    _hist *
    _hist::
    Divide(const double & rhs, bool inplace) {
        _hist * tmp;
        if(inplace) tmp = this;
        else tmp = this->Clone();
        tmp->_divide(rhs);
        return tmp;
    }

    _hist *
    _hist::
    Multiply(const _hist * rhs, bool inplace) {
        this->EnsureConsistentBinning(rhs, __PRETTY_FUNCTION__);
        _hist * tmp;
        if(inplace) tmp = this;
        else tmp = this->Clone();
        tmp->_multiply(rhs);
        return tmp;
    }

    _hist *
    _hist::
    Multiply(const double & rhs, bool inplace) {
        _hist * tmp;
        if(inplace) tmp = this;
        else tmp = this->Clone();
        tmp->_multiply(rhs);
        return tmp;
    }

    _hist *
    _hist::
    Subtract(const _hist * rhs) const {
        return this->Clone()->Subtract(rhs);
    }

    _hist *
    _hist::
    Subtract(const double & rhs) const {
        return this->Clone()->Subtract(rhs);
    }

    _hist *
    _hist::
    Add(const _hist * rhs) const {
        return this->Clone()->Add(rhs);
    }

    _hist *
    _hist::
    Add(const double & rhs) const {
        return this->Clone()->Add(rhs);
    }

    _hist *
    _hist::
    Divide(const _hist * rhs) const {
        return this->Clone()->Divide(rhs);
    }

    _hist *
    _hist::
    Divide(const double & rhs) const {
        return this->Clone()->Divide(rhs);
    }

    _hist *
    _hist::
    Multiply(const _hist * rhs) const {
        return this->Clone()->Multiply(rhs);
    }

    _hist *
    _hist::
    Multiply(const double & rhs) const {
        return this->Clone()->Multiply(rhs);
    }
}
