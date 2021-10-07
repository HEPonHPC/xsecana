#pragma once

#include <memory>
#include <Eigen/Dense>

#include "TDirectory.h"
#include "TObjString.h"
#include "TH1D.h"
#include "TH1F.h"

#include <iostream>
#include <exception>

namespace xsec {
    namespace exceptions {
        class InconsistentBinningError : public std::exception {
        public:
            InconsistentBinningError(const char * caller, double tol) {
                std::sprintf(fMsg,
                             "\nAttempted histogram operation with \ninconsistent binning in %s.\n Exceeded tolerance %1.4g",
                             caller,
                             tol);
            }

            [[nodiscard]] const char * what() const throw() override {
                return fMsg;
            }
        private:
            char fMsg[500];
        };
    }
    /// \brief Object representing a filled histogram used internally by the framework
    /// Bin contents are indexed starting at 1
    /// Underflow is index 0
    /// Overflow is index N+1 where N is number of bins
    /// Wraps Eigen arrays for fast mathematical operations especially linear algebra
    ///
    /// TODO Users can create conversion functions to/from this object
    typedef Eigen::ArrayXd Array;
    class Hist {
    public:
        Hist() {}

        /// \brief ctor using existing arrays for contents and edges
        /// including under/overflow
        Hist(const Array & contents_and_uof,
             const Array & edges_and_uof,
             const Array & errors_and_uof,
             const double & exposure = 1)
                : fContentsAndUOF(contents_and_uof),
                  fEdgesAndUOF(edges_and_uof),
                  fErrorsAndUOF(errors_and_uof),
                  fExposure(exposure)
        {
            assert(contents_and_uof.size()+1 == edges_and_uof.size() &&
                   contents_and_uof.size() == errors_and_uof.size() &&
                   "Incompatible edges, contents, and/or errors");
        }

        /// \brief ctor using existing arrays for contents and edges
        /// including under/overflow
        /// Initializes errors to 0
        Hist(const Array & contents_and_uof,
             const Array & edges_and_uof,
             const double & exposure = 1)
                : fContentsAndUOF(contents_and_uof),
                  fEdgesAndUOF(edges_and_uof),
                  fErrorsAndUOF(Array::Zero(fContentsAndUOF.size())),
                  fExposure(exposure)
        {
            assert(contents_and_uof.size()+1 == edges_and_uof.size() &&
                   "Incompatible edges, contents, and/or errors");
        }

        /// \brief convenience constructor.
        /// works for dynamic and fixed-size histograms
        /// adds two additional "bins" for under/overflow
        /// to the underlaying arrays
        Hist(const int & nbins,
             const double & min,
             const double & max);

        Hist(const Hist & rhs)
                : fContentsAndUOF(rhs.ContentsAndUOF()),
                  fEdgesAndUOF(rhs.EdgesAndUOF()),
                  fErrorsAndUOF(rhs.ErrorsAndUOF()),
                  fExposure(rhs.Exposure())
        {}

        virtual Hist BinWidthNormalize() const;
        virtual Hist & BinWidthNormalize();

        virtual Hist AreaNormalize() const;
        virtual Hist & AreaNormalize();

        virtual double Integrate() const;

        virtual Hist operator-(const Hist & rhs) const;

        virtual Hist operator+(const Hist & rhs) const;

        virtual Hist operator/(const Hist & rhs) const;

        virtual Hist operator*(const Hist & rhs) const;

        virtual bool operator==(const Hist & rhs) const;

        virtual bool operator!=(const Hist & rhs) const { return !(*this == rhs); }

        virtual Hist operator-=(const Hist & rhs);

        virtual Hist operator+=(const Hist & rhs);

        virtual Hist operator/=(const Hist & rhs);

        virtual Hist operator*=(const Hist & rhs);

        virtual Hist operator-(const double & rhs) const;

        virtual Hist operator+(const double & rhs) const;

        virtual Hist operator/(const double & rhs) const;

        virtual Hist operator*(const double & rhs) const;

        virtual Hist operator-=(const double & rhs);

        virtual Hist operator+=(const double & rhs);

        virtual Hist operator/=(const double & rhs);

        virtual Hist operator*=(const double & rhs);

        virtual Hist operator=(const Hist & rhs);

        virtual Hist operator=(Hist && rhs);

        virtual Hist ScaleByExposure(double new_expo) const;

        virtual Hist TrueDivide(const Hist & rhs) const;

        // some convenience functions
        virtual Hist abs() const;

        virtual Hist abs2() const;

        virtual Hist sqrt() const;

        virtual Hist pow(double exp) const;

        virtual const Array Contents() const
        { return fContentsAndUOF(Eigen::seq(1, fContentsAndUOF.size()-2)); }
        virtual const Array & ContentsAndUOF() const { return fContentsAndUOF; }

        virtual const Array Errors() const
        { return fErrorsAndUOF(Eigen::seq(1, fErrorsAndUOF.size()-2)); }
        virtual const Array & ErrorsAndUOF() const { return fErrorsAndUOF; }

        virtual const Array Edges() const
        { return fEdgesAndUOF(Eigen::seq(1, fEdgesAndUOF.size()-2)); }
        virtual const Array & EdgesAndUOF() const { return fEdgesAndUOF; }

        virtual void SetContentsAndUOF(const Array & contents_and_uof) {
            assert(fContentsAndUOF.size() == contents_and_uof.size() &&
                   "Incompatible contents array");
            fContentsAndUOF = contents_and_uof;
        }

        virtual void SetErrorsAndUOF(const Array & errors_and_uof) {
            assert(fErrorsAndUOF.size() == errors_and_uof.size() &&
                   "Incompatible contents array");
            fErrorsAndUOF = errors_and_uof;
        }

        virtual double Exposure() const { return fExposure; }

        virtual Array BinWidths() const;

        virtual double & operator[](int index) { return fContentsAndUOF(index); }

        virtual double operator[](int index) const { return fContentsAndUOF(index); }

        virtual void SaveTo(TDirectory * dir, const std::string & subdir) const;

        static std::unique_ptr<Hist> LoadFrom(TDirectory * dir, const std::string & subdir);

    private:
        void EnsureConsistentBinning(const Hist & rhs, const char * caller, double tol = 1e-5) const;

        Array fContentsAndUOF;
        Array fErrorsAndUOF;
        Array fEdgesAndUOF;
        double fExposure = 1;
    };
    
    // ROOT interface
    // we're still dependent enough on ROOT for this to be here
    // but eventually all ROOT things will be put into an optional interface
    namespace root {
        inline TH1 * ToTH1(const Hist & hist,
                           const std::string & name = "",
                           const std::string & title = "") {

            TH1 * h = new TH1D(name.c_str(),
                               title.c_str(),
                               hist.Edges().size() - 1,
                               hist.Edges().data());


            for (auto i = 0u; i < hist.ContentsAndUOF().size(); i++) {
                h->SetBinContent(i, hist.ContentsAndUOF()(i));
                h->SetBinError(i, hist.ErrorsAndUOF()(i));
            }

            return h;
        }

        Hist
        FromTH1(const TH1 * h, double exposure) {
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
            edges(nedges-1) = h->GetBinLowEdge(h->GetNbinsX() + 2);

            return Hist(std::move(contents),
                        std::move(edges),
                        std::move(errors),
                        exposure);
        }
    }

    /////////////////////////////////////////////////////////
    /// \brief Ensure that the rhs histogram has binning
    /// consistent with this histogram.
    /// Raises
    void
    Hist::
    EnsureConsistentBinning(const Hist & rhs, const char * caller, double tol) const {
        if(!(this->fEdgesAndUOF - rhs.fEdgesAndUOF).isZero(tol)) {
            std::cout << this->fEdgesAndUOF << std::endl;
            std::cout << rhs.fEdgesAndUOF << std::endl;
            throw exceptions::InconsistentBinningError(__PRETTY_FUNCTION__, tol);
        }
    }

    /////////////////////////////////////////////////////////
    /// \brief Divide this histogram by rhs without scaling rhs by
    /// this exposure.
    /// Returns histogram at this exposure
    Hist
    Hist::
    TrueDivide(const Hist & rhs) const {
        EnsureConsistentBinning(rhs, __PRETTY_FUNCTION__);;

        return Hist(fContentsAndUOF / rhs.ContentsAndUOF(),
                                  fEdgesAndUOF,
                                  (fErrorsAndUOF.pow(2) + rhs.fErrorsAndUOF.pow(2)).sqrt(),
                                  fExposure);
    }

    /////////////////////////////////////////////////////////
    Hist
    Hist::
    ScaleByExposure(double new_expo) const {
        return Hist(fContentsAndUOF * (new_expo / fExposure),
                                  fEdgesAndUOF,
                                  fErrorsAndUOF * (new_expo / fExposure),
                                  new_expo);
    }

    /////////////////////////////////////////////////////////
    double
    Hist::
    Integrate() const {
        return fContentsAndUOF.sum();
    }

    /////////////////////////////////////////////////////////
    Hist::
    Hist(const int & nbins,
         const double & min,
         const double & max) {
        fContentsAndUOF = Array::Zero(nbins+2);
        fErrorsAndUOF = Array::Zero(nbins+2);

        auto step = (max - min) / nbins;
        fEdgesAndUOF = Array::LinSpaced(nbins + 3,
                                        min - step,
                                        max + step);
    }

    /////////////////////////////////////////////////////////
    void
    Hist::
    SaveTo(TDirectory * dir, const std::string & subdir) const {
        TDirectory * tmp = gDirectory;
        dir->mkdir(subdir.c_str());
        dir = dir->GetDirectory(subdir.c_str());
        dir->cd();

        auto h = root::ToTH1(*this);
        h->Write("hist");

        // exposure is saved in a histogram so it gets accumulated with ROOT's hadd
        TH1 * exposure = new TH1D("", "", 1, 0, 1);

        exposure->SetBinContent(1, fExposure);
        exposure->Write("exposure");

        tmp->cd();
    }

    /////////////////////////////////////////////////////////
    std::unique_ptr<Hist>
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

        return std::make_unique<Hist >(root::FromTH1(h, exposure->GetBinContent(1)));
    }

    /////////////////////////////////////////////////////////
    bool
    Hist::
    operator==(const Hist & rhs) const {
        auto scale_exposure = this->fExposure / rhs.fExposure;
        return (this->fContentsAndUOF - rhs.fContentsAndUOF * scale_exposure).isZero(0) &&
               (this->fEdgesAndUOF - rhs.fEdgesAndUOF).isZero(0) &&
               (this->fErrorsAndUOF - rhs.fErrorsAndUOF * scale_exposure).isZero(0) &&
               (this->fExposure == rhs.fExposure);
    }

    /////////////////////////////////////////////////////////
    /// \brief Divide bin contents by corresponding bin widths in place.
    /// not including over/underflow
    Hist &
    Hist::
    BinWidthNormalize() {
        fContentsAndUOF(Eigen::seqN(1, fContentsAndUOF.size() - 2)) /= this->BinWidths();
        return *this;
    }

    /////////////////////////////////////////////////////////
    /// \brief Return new histogram with bin contents
    /// divided by corresponding bin widths.
    /// not including over/underflow
    Hist
    Hist::
    BinWidthNormalize() const {
        auto copy = *this;
        copy = copy.BinWidthNormalize();
        return copy;
    }

    /////////////////////////////////////////////////////////
    /// \brief Divide bin contents and under/overflow by
    /// sum of bin contents in place
    Hist &
    Hist::
    AreaNormalize() {
        fContentsAndUOF /= this->Contents().sum();
        return *this;
    }

    /////////////////////////////////////////////////////////
    /// \brief Divide bin contents and under/overflow by
    /// sum of bin contents in place
    Hist
    Hist::
    AreaNormalize() const {
        auto copy = *this;
        copy = copy.AreaNormalize();
        return copy;
    }

    /////////////////////////////////////////////////////////
    /// \brief Return an array of bin widths
    /// not including over/underflow
    Array
    Hist::
    BinWidths() const {
        return fEdgesAndUOF(Eigen::seqN(3, fContentsAndUOF.size()-2)) -
               fEdgesAndUOF(Eigen::seqN(2, fContentsAndUOF.size()-2));
    }

    /////////////////////////////////////////////////////////
    /// \brief Return a new histogram who's contents and errors
    /// are the absolute values of this histogram
    Hist
    Hist::
    abs() const {
        return Hist(this->fContentsAndUOF.abs(),
                                  this->fEdgesAndUOF,
                                  this->fErrorsAndUOF.abs(),
                                  this->fExposure);
    }

    /////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////
    /// \brief Return a new histogram who's contents and errors
    /// are the absolute values squared of this histogram
    Hist
    Hist::
    abs2() const {
        return Hist(this->fContentsAndUOF.abs2(),
                                  this->fEdgesAndUOF,
                                  this->fErrorsAndUOF.abs2(),
                                  this->fExposure);
    }

    /////////////////////////////////////////////////////////
    /// \brief Return a new histogram who's contents
    /// are the square root of this histogram
    /// and errors are divided by two
    Hist
    Hist::
    sqrt() const {
        return Hist(this->fContentsAndUOF.sqrt(),
                                  this->fEdgesAndUOF,
                                  this->fErrorsAndUOF / 2.,
                                  std::sqrt(this->fExposure));
    }

    /////////////////////////////////////////////////////////
    /// \brief Return a new histogram who's contents
    /// are raised to the exp power of this histogram
    /// and errors are multiplied by exp
    Hist
    Hist::
    pow(double exp) const {
        return Hist(this->fContentsAndUOF.pow(exp),
                                  this->fEdgesAndUOF,
                                  this->fErrorsAndUOF * exp,
                                  std::pow(this->fExposure, exp));
    }

    /////////////////////////////////////////////////////////
    Hist
    Hist::
    operator-(const Hist & rhs) const {
        Hist ret = *this; // copy this
        ret -= rhs;
        return ret;
    }

    /////////////////////////////////////////////////////////
    Hist
    Hist::
    operator-=(const Hist & rhs) {
        EnsureConsistentBinning(rhs, __PRETTY_FUNCTION__);;

        this->fContentsAndUOF -= rhs.fContentsAndUOF * fExposure / rhs.fExposure;
        this->fErrorsAndUOF = (this->fErrorsAndUOF.pow(2) +
                               (rhs.fErrorsAndUOF * fExposure / rhs.fExposure).pow(2)).sqrt();
        return *this;
    }

    /////////////////////////////////////////////////////////
    Hist
    Hist::
    operator+(const Hist & rhs) const {
        Hist ret = *this; // copy this
        ret += rhs;
        return ret;
    }

    /////////////////////////////////////////////////////////
    Hist
    Hist::
    operator+=(const Hist & rhs) {
        EnsureConsistentBinning(rhs, __PRETTY_FUNCTION__);;

        this->fContentsAndUOF += rhs.fContentsAndUOF * fExposure / rhs.fExposure;
        this->fErrorsAndUOF = (this->fErrorsAndUOF.pow(2) +
                               (rhs.fErrorsAndUOF * fExposure / rhs.fExposure).pow(2)).sqrt();
        return *this;
    }

    /////////////////////////////////////////////////////////
    Hist
    Hist::
    operator*(const Hist & rhs) const {
        Hist ret = *this; // copy this
        ret *= rhs;
        return ret;
    }

    /////////////////////////////////////////////////////////
    Hist
    Hist::
    operator*=(const Hist & rhs) {
        EnsureConsistentBinning(rhs, __PRETTY_FUNCTION__);;

        this->fContentsAndUOF *= rhs.fContentsAndUOF * fExposure / rhs.fExposure;
        this->fErrorsAndUOF = (this->fErrorsAndUOF.pow(2) +
                               (rhs.fErrorsAndUOF * fExposure / rhs.fExposure).pow(2)).sqrt();
        return *this;
    }

    /////////////////////////////////////////////////////////
    Hist
    Hist::
    operator/(const Hist & rhs) const {
        Hist ret = *this; // copy this
        ret /= rhs;
        return ret;
    }

    /////////////////////////////////////////////////////////
    Hist
    Hist::
    operator/=(const Hist & rhs) {
        EnsureConsistentBinning(rhs, __PRETTY_FUNCTION__);;

        this->fContentsAndUOF /= rhs.fContentsAndUOF * fExposure / rhs.fExposure;
        this->fErrorsAndUOF = (this->fErrorsAndUOF.pow(2) +
                               (rhs.fErrorsAndUOF * fExposure / rhs.fExposure).pow(2)).sqrt();
        this->fExposure = 1;
        return *this;
    }

    /////////////////////////////////////////////////////////
    Hist
    Hist::
    operator-(const double & rhs) const {
        Hist ret = *this; // copy this
        ret -= rhs;
        return ret;
    }

    /////////////////////////////////////////////////////////
    Hist
    Hist::
    operator-=(const double & rhs) {
        this->fContentsAndUOF -= rhs;
        return *this;
    }

    /////////////////////////////////////////////////////////
    Hist
    Hist::
    operator+(const double & rhs) const {
        Hist ret = *this; // copy this
        ret += rhs;
        return ret;
    }

    /////////////////////////////////////////////////////////
    Hist
    Hist::
    operator+=(const double & rhs) {
        this->fContentsAndUOF += rhs;
        return *this;
    }

    /////////////////////////////////////////////////////////
    Hist
    Hist::
    operator*(const double & rhs) const {
        Hist ret = *this; // copy this
        ret *= rhs;
        return ret;
    }

    /////////////////////////////////////////////////////////
    Hist
    Hist::
    operator/(const double & rhs) const {
        Hist ret = *this; // copy this
        ret /= rhs;
        return ret;
    }

    /////////////////////////////////////////////////////////
    Hist
    Hist::
    operator/=(const double & rhs) {
        this->fContentsAndUOF /= rhs;
        return *this;
    }


    /////////////////////////////////////////////////////////
    Hist
    Hist::
    operator*=(const double & rhs) {
        this->fContentsAndUOF *= rhs;
        return *this;
    }

    /////////////////////////////////////////////////////////
    Hist
    Hist::
    operator=(const Hist & rhs) {
        if (this == &rhs) return *this;
        fContentsAndUOF = rhs.fContentsAndUOF;
        fEdgesAndUOF = rhs.fEdgesAndUOF;
        fErrorsAndUOF = rhs.fErrorsAndUOF;
        fExposure = rhs.fExposure;
        return *this;
    }

    /////////////////////////////////////////////////////////
    Hist
    Hist::
    operator=(Hist && rhs) {
        if(this == &rhs) return *this;
        fContentsAndUOF = std::move(rhs.fContentsAndUOF);
        fEdgesAndUOF = std::move(rhs.fEdgesAndUOF);
        fErrorsAndUOF = std::move(rhs.fErrorsAndUOF);
        fExposure = rhs.fExposure;
        return *this;
    }
}
