//
// Created by Derek Doyle on 10/9/21.
//
#include "XSecAna/Hist.h"
#include "TH1.h"

namespace xsec {
    Hist *
    Hist::
    FromTH1(const TH1 * h, const double & exposure) {
        const unsigned int nedges = h->GetNbinsX() + 3;
        const unsigned int nbins_and_uof = h->GetNbinsX() + 2;

        Array edges = Array::Zero(nedges);
        Array contents = Array::Zero(nbins_and_uof);
        Array errors = Array::Zero(nbins_and_uof);

        for (auto i = 0u; i < nbins_and_uof; i++) {
            edges(i) = h->GetBinLowEdge(i);
            contents(i) = h->GetBinContent(i);
            errors(i) = h->GetBinError(i);
        }
        edges(nedges - 1) = h->GetBinLowEdge(h->GetNbinsX() + 2);

        return new Hist(std::move(contents),
                        std::move(edges),
                        std::move(errors),
                        exposure);
    }

    TH1 *
    Hist::
    ToTH1(const std::string & name, const std::string & title) const {
        TH1 * h = new TH1D(name.c_str(),
                           title.c_str(),
                           this->GetEdges().size() - 1,
                           this->GetEdges().data());

        for (auto i = 0u; i < this->GetContentsAndUOF().size(); i++) {
            h->SetBinContent(i, this->GetContentsAndUOF()(i));
            h->SetBinError(i, this->GetErrorsAndUOF()(i));
        }
        return h;
    }

    /////////////////////////////////////////////////////////
    Hist::
    Hist(const int & nbins,
         const double & min,
         const double & max,
         const double & exposure) {
        fContentsAndUOF = WeightedArray(Array::Zero(nbins + 2),
                                        exposure);
        fErrorsAndUOF = Array::Zero(nbins + 2);

        auto step = (max - min) / nbins;
        fEdgesAndUOF = Array::LinSpaced(nbins + 3,
                                        min - step,
                                        max + step);
    }

    /// \brief ctor using existing arrays for contents and edges
    /// including under/overflow
    Hist::
    Hist(const Array & contents_and_uof,
         const Array & edges_and_uof,
         const Array & errors_and_uof,
         const double & exposure)
            : fContentsAndUOF(contents_and_uof, exposure),
              fEdgesAndUOF(edges_and_uof),
              fErrorsAndUOF(errors_and_uof) {
        assert(contents_and_uof.size() + 1 == edges_and_uof.size() &&
               contents_and_uof.size() == errors_and_uof.size() &&
               "Incompatible edges, contents, and/or errors");
    }

    Hist::
    Hist(const WeightedArray & contents_and_uof,
         const Array & edges_and_uof,
         const Array & errors_and_uof)
            : fContentsAndUOF(contents_and_uof),
              fEdgesAndUOF(edges_and_uof),
              fErrorsAndUOF(errors_and_uof) {
        assert(contents_and_uof.size() + 1 == edges_and_uof.size() &&
               contents_and_uof.size() == errors_and_uof.size() &&
               "Incompatible edges, contents, and/or errors");
    }


    /// \brief ctor using existing arrays for contents and edges
    /// including under/overflow
    /// Initializes errors to 0
    Hist::
    Hist(const Array & contents_and_uof,
         const Array & edges_and_uof,
         const double & exposure)
            : fContentsAndUOF(contents_and_uof, exposure),
              fEdgesAndUOF(edges_and_uof),
              fErrorsAndUOF(Array::Zero(fContentsAndUOF.size())) {
        assert(contents_and_uof.size() + 1 == edges_and_uof.size() &&
               "Incompatible edges, contents, and/or errors");
    }

    Hist::
    Hist(const WeightedArray & contents_and_uof,
         const Array & edges_and_uof)
            : fContentsAndUOF(contents_and_uof),
              fEdgesAndUOF(edges_and_uof),
              fErrorsAndUOF(Array::Zero(fContentsAndUOF.size())) {
        assert(contents_and_uof.size() + 1 == edges_and_uof.size() &&
               "Incompatible edges, contents, and/or errors");
    }

    Hist::
    Hist(const Hist & rhs)
            : fContentsAndUOF(rhs.fContentsAndUOF),
              fEdgesAndUOF(rhs.fEdgesAndUOF),
              fErrorsAndUOF(rhs.fErrorsAndUOF) {}


    /////////////////////////////////////////////////////////
    /// \brief Return an array of bin widths
    /// not including over/underflow
    Array
    Hist::
    GetBinWidths() const {
        return fEdgesAndUOF(Eigen::seqN(3, fContentsAndUOF.size() - 2)) -
               fEdgesAndUOF(Eigen::seqN(2, fContentsAndUOF.size() - 2));
    }

    const Array
    Hist::
    GetEdges() const {
        return fEdgesAndUOF(Eigen::seq(1, fEdgesAndUOF.size() - 2));
    }

    const Array
    Hist::
    GetEdgesAndUOF() const {
        return fEdgesAndUOF;
    }

    void
    Hist::
    SetContentsAndUOF(const Array & contents_and_uof) {
        assert(fContentsAndUOF.size() == contents_and_uof.size() &&
               "Incompatible contents array");
        fContentsAndUOF = WeightedArray(contents_and_uof,
                                        this->Exposure());
    }

    void
    Hist::
    SetContents(const Array & contents) {
        assert(fContentsAndUOF.size() == contents.size() &&
               "Incompatible contents array");
        fContentsAndUOF.array()(Eigen::seq(1, fContentsAndUOF.size() - 2)) = contents;
    }

    void
    Hist::
    SetErrorsAndUOF(const Array & errors_and_uof) {
        assert(fErrorsAndUOF.size() == errors_and_uof.size() &&
               "Incompatible contents array");
        fErrorsAndUOF = errors_and_uof;
    }

    void
    Hist::
    SetErrors(const Array & contents) {
        assert(fContentsAndUOF.size() == contents.size() &&
               "Incompatible contents array");
        fErrorsAndUOF.array()(Eigen::seq(1, fContentsAndUOF.size() - 2)) = contents;
    }

    void
    Hist::
    SetExposure(const double & new_exposure) {
        this->fContentsAndUOF.set_weight(new_exposure);
    }

    double &
    Hist::
    operator()(int index) {
        return fContentsAndUOF.array()(index);
    }

    double
    Hist::
    operator()(int index) const {
        return fContentsAndUOF.array()(index);
    }

    double
    Hist::
    Exposure() const {
        return fContentsAndUOF.weight();
    }

    _hist *
    Hist::
    Clone() const {
        return new Hist(*this);
    }

    double
    Hist::
    Integrate() const {
        return this->fContentsAndUOF.array().sum();
    }


    /////////////////////////////////////////////////////////
    std::unique_ptr<_hist>
    Hist::
    LoadFrom(TDirectory * dir, const std::string & subdir) {
        dir = dir->GetDirectory(subdir.c_str());
        dir->cd();

        TH1 * h = (TH1D *) dir->Get("hist");
        TH1 * exposure = (TH1D *) dir->Get("exposure");

        if (!h) {
            std::cerr << "Object TH1 was not found in " << dir->GetPath() << std::endl;
            exit(1);
        }

        return std::unique_ptr<Hist>(Hist::FromTH1(h, exposure->GetBinContent(1)));
    }

    Array
    Hist::
    GetContents() const {
        return fContentsAndUOF.array()(Eigen::seq(1, fContentsAndUOF.size() - 2));
    }

    Array
    Hist::
    GetContentsAndUOF() const {
        return fContentsAndUOF.array();
    }

    Array
    Hist::
    GetErrors() const {
        return fErrorsAndUOF.array()(Eigen::seq(1, fErrorsAndUOF.size() - 2));
    }

    Array
    Hist::
    GetErrorsAndUOF() const {
        return fErrorsAndUOF;
    }

    /////////////////////////////////////////////////////////
    /// \brief Ensure that the rhs histogram has binning
    /// consistent with this histogram.
    /// Raises
    bool
    Hist::
    _is_same_contents(const _hist * rhs, const double & tol) const {
        return (this->fContentsAndUOF - dynamic_cast<const Hist *>(rhs)->fContentsAndUOF).array().isZero(tol);
    }

    /////////////////////////////////////////////////////////
    /// \brief Ensure that the rhs histogram has binning
    /// consistent with this histogram.
    /// Raises
    bool
    Hist::
    _is_same_binning(const _hist * rhs, const double & tol) const {
        return (this->fEdgesAndUOF - dynamic_cast<const Hist *>(rhs)->fEdgesAndUOF).isZero(tol);
    }

    _hist *
    Hist::
    _scale_by_exposure(const double & new_expo) {
        fContentsAndUOF = fContentsAndUOF.ScaleByWeight(new_expo);
        return this;
    }

    _hist *
    Hist::
    _true_divide(const _hist * rhs) {
        fContentsAndUOF = WeightedArray(this->fContentsAndUOF.array() /
                                        dynamic_cast<const Hist *>(rhs)->fContentsAndUOF.array(),
                                        this->fContentsAndUOF.weight());
        return this;
    }

    _hist *
    Hist::
    _bin_width_normalize() {
        fContentsAndUOF.array()(Eigen::seqN(1, fContentsAndUOF.size() - 2)) /= this->GetBinWidths();
        return this;
    }

    _hist *
    Hist::
    _area_normalize() {
        fContentsAndUOF /= this->GetContents().sum();
        return this;
    }

    _hist *
    Hist::
    _subtract(const _hist * rhs) {
        this->fContentsAndUOF -= dynamic_cast<const Hist *>(rhs)->fContentsAndUOF;
        this->fErrorsAndUOF = (this->fErrorsAndUOF.pow(2) + (
                dynamic_cast<const Hist *>(rhs)->fErrorsAndUOF * this->Exposure() / rhs->Exposure()
        ).pow(2)).sqrt();
        return this;
    }

    _hist *
    Hist::
    _subtract(const double & rhs) {
        this->fContentsAndUOF -= rhs;
        return this;
    }

    _hist *
    Hist::
    _add(const _hist * rhs) {
        this->fContentsAndUOF += dynamic_cast<const Hist *>(rhs)->fContentsAndUOF;
        this->fErrorsAndUOF = (this->fErrorsAndUOF.pow(2) + (
                dynamic_cast<const Hist *>(rhs)->fErrorsAndUOF * this->Exposure() / rhs->Exposure()
        ).pow(2)).sqrt();
        return this;
    }

    _hist *
    Hist::
    _add(const double & rhs) {
        this->fContentsAndUOF += rhs;
        return this;
    }

    _hist *
    Hist::
    _divide(const _hist * rhs) {
        this->fContentsAndUOF /= dynamic_cast<const Hist *>(rhs)->fContentsAndUOF;
        this->fErrorsAndUOF = (this->fErrorsAndUOF.pow(2) + (
                dynamic_cast<const Hist *>(rhs)->fErrorsAndUOF * this->Exposure() / rhs->Exposure()
        ).pow(2)).sqrt();
        return this;
    }

    _hist *
    Hist::
    _divide(const double & rhs) {
        this->fContentsAndUOF /= rhs;
        return this;
    }

    _hist *
    Hist::
    _multiply(const _hist * rhs) {
        this->fContentsAndUOF *= dynamic_cast<const Hist *>(rhs)->fContentsAndUOF;
        this->fErrorsAndUOF = (this->fErrorsAndUOF.pow(2) + (
                dynamic_cast<const Hist *>(rhs)->fErrorsAndUOF * this->Exposure() / rhs->Exposure()
        ).pow(2)).sqrt();
        return this;
    }

    _hist *
    Hist::
    _multiply(const double & rhs) {
        this->fContentsAndUOF *= rhs;
        return this;
    }

    _hist *
    Hist::
    _abs() {
        fContentsAndUOF = fContentsAndUOF.abs();
        return this;
    }

    _hist *
    Hist::
    _abs2() {
        fContentsAndUOF = fContentsAndUOF.abs2();
        fErrorsAndUOF *= 2;
        return this;
    }

    _hist *
    Hist::
    _sqrt() {
        fContentsAndUOF = fContentsAndUOF.sqrt();
        fErrorsAndUOF /= 2;
        return this;
    }

    _hist *
    Hist::
    _pow(double exp) {
        fContentsAndUOF = fContentsAndUOF.pow(exp);
        fErrorsAndUOF *= exp;
        return this;
    }
}
