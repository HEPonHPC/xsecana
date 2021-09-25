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
    // compile time expression for determining size of Eigen Array holding
    // N bins + 2 over/underflow
    constexpr int ContentsAndUOFSize(int Cols) { return Cols == Eigen::Dynamic ? Eigen::Dynamic : Cols + 2; }

    /// \brief Object representing a filled histogram used internally by the framework
    /// Bin contents are indexed starting at 1
    /// Underflow is index 0
    /// Overflow is index N+1 where N is number of bins
    /// Wraps Eigen arrays for fast mathematical operations especially linear algebra
    ///
    /// TODO Users can create conversion functions to/from this object
    /// template parameters are forwarded to the underlaying Eigen arrays
    /// Cols is the number of bins
    template<typename Scalar,
            int Cols = -1>
    class Hist {
    public:
        typedef Scalar scalar_type;
        typedef Eigen::Array<Scalar, 1, Cols> array_type;
        typedef Eigen::Array<Scalar, 1, ContentsAndUOFSize(Cols)> array_and_uof_type;        
        typedef Eigen::Array<Scalar, 1, Eigen::Dynamic> edges_type;

        Hist() {}

        /// \brief ctor using existing arrays for contents and edges
        /// including under/overflow
        Hist(const array_and_uof_type & contents_and_uof,
             const edges_type & edges_and_uof,
             const array_and_uof_type & errors_and_uof,
             const Scalar & exposure = 1)
                : fContentsAndUOF(contents_and_uof),
                  fEdgesAndUOF(edges_and_uof),
                  fErrorsAndUOF(errors_and_uof),
                  fExposure(exposure)
        {
            assert(contents_and_uof.size() == edges_and_uof.size()+1 &&
                   contents_and_uof.size() == edges_and_uof.size() &&
                   "Incompatible edges, contents, and/or errors");
        }

        /// \brief ctor using existing arrays for contents and edges
        /// including under/overflow
        /// Initializes errors to 0
        Hist(const array_and_uof_type & contents_and_uof,
             const edges_type & edges_and_uof,
             const Scalar & exposure = 1)
                : fContentsAndUOF(contents_and_uof),
                  fEdgesAndUOF(edges_and_uof),
                  fErrorsAndUOF(array_and_uof_type::Zero(fContentsAndUOF.size())),
                  fExposure(exposure)
        {
            assert(contents_and_uof.size() == edges_and_uof.size()+1 &&
                   "Incompatible edges, contents, and/or errors");
        }

        /// \brief convenience constructor.
        /// works for dynamic and fixed-size histograms
        /// adds two additional "bins" for under/overflow
        /// to the underlaying arrays
        Hist(const int & nbins,
             const Scalar & min,
             const Scalar & max);

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

        virtual Scalar Integrate() const;

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

        virtual Hist operator-(const Scalar & rhs) const;

        virtual Hist operator+(const Scalar & rhs) const;

        virtual Hist operator/(const Scalar & rhs) const;

        virtual Hist operator*(const Scalar & rhs) const;

        virtual Hist operator-=(const Scalar & rhs);

        virtual Hist operator+=(const Scalar & rhs);

        virtual Hist operator/=(const Scalar & rhs);

        virtual Hist operator*=(const Scalar & rhs);

        virtual Hist operator=(const Hist & rhs);

        virtual Hist operator=(Hist && rhs);

        virtual Hist ScaleByExposure(Scalar new_expo) const;

        virtual Hist TrueDivide(const Hist & rhs) const;

        // some convenience functions
        virtual Hist abs() const;

        virtual Hist abs2() const;

        virtual Hist sqrt() const;

        virtual Hist pow(Scalar exp) const;

        virtual const array_type Contents() const
        { return fContentsAndUOF(Eigen::seq(1, fContentsAndUOF.size()-2)); }
        virtual const array_and_uof_type & ContentsAndUOF() const { return fContentsAndUOF; }

        virtual const array_type Errors() const
        { return fErrorsAndUOF(Eigen::seq(1, fErrorsAndUOF.size()-2)); }
        virtual const array_and_uof_type & ErrorsAndUOF() const { return fErrorsAndUOF; }

        virtual const edges_type Edges() const
        { return fEdgesAndUOF(Eigen::seq(1, fEdgesAndUOF.size()-2)); }
        virtual const edges_type & EdgesAndUOF() const { return fEdgesAndUOF; }

        virtual void SetContentsAndUOF(const array_and_uof_type & contents_and_uof) {
            assert(fContentsAndUOF.size() == contents_and_uof.size() &&
                   "Incompatible contents array");
            fContentsAndUOF = contents_and_uof;
        }

        virtual void SetErrorsAndUOF(const array_and_uof_type & errors_and_uof) {
            assert(fErrorsAndUOF.size() == errors_and_uof.size() &&
                   "Incompatible contents array");
            fErrorsAndUOF = errors_and_uof;
        }

        virtual Scalar Exposure() const { return fExposure; }

        virtual Eigen::Array<Scalar, 1, Cols> BinWidths() const;

        virtual Scalar & operator[](int index) { return fContentsAndUOF(index); }

        virtual Scalar operator[](int index) const { return fContentsAndUOF(index); }

        virtual void SaveTo(TDirectory * dir, const std::string & subdir) const;

        static std::unique_ptr<Hist> LoadFrom(TDirectory * dir, const std::string & subdir);

    private:
        void EnsureConsistentBinning(const Hist & rhs, const char * caller, double tol = 1e-5) const;

        array_and_uof_type fContentsAndUOF;
        array_and_uof_type fErrorsAndUOF;
        edges_type fEdgesAndUOF;
        Scalar fExposure = 1;
    };

    typedef Hist<double, Eigen::Dynamic> HistXd;
    typedef Hist<float, Eigen::Dynamic> HistXf;

    // ROOT interface
    // we're still dependent enough on ROOT for this to be here
    // but eventually all ROOT things will be put into an optional interface
    namespace root {
        template<typename Scalar, int Cols>
        inline TH1 * ToTH1(const Hist<Scalar, Cols> & hist,
                           const std::string & name = "",
                           const std::string & title = "") {
            TH1 * h;
            if constexpr(std::is_same<Scalar, double>::value) {
                h = new TH1D(name.c_str(),
                             title.c_str(),
                             hist.Edges().size()-1,
                             hist.Edges().data());
            } else if constexpr(std::is_same<Scalar, float>::value) {
                h = new TH1F(name.c_str(),
                             title.c_str(),
                             hist.Edges().size()-1,
                             hist.Edges().data());
            } else {
                std::cerr << "ToTH1 not implemented for this type" << std::endl;
                exit(1);
            }

            for (auto i = 0u; i < hist.ContentsAndUOF().size(); i++) {
                h->SetBinContent(i, hist.ContentsAndUOF()(i));
                h->SetBinError(i, hist.ErrorsAndUOF()(i));
            }

            return h;
        }

        template<typename Scalar, int Cols>
        Hist<Scalar, Cols>
        FromTH1(const TH1 * h, Scalar exposure) {
            const unsigned int nedges = h->GetNbinsX() + 3;
            const unsigned int nbins_and_uof = h->GetNbinsX() + 2;

            typename Hist<Scalar, Cols>::edges_type edges = Hist<Scalar, Cols>::edges_type::Zero(nedges);
            typename Hist<Scalar, Cols>::array_and_uof_type contents = Hist<Scalar, Cols>::array_and_uof_type::Zero(nbins_and_uof);
            typename Hist<Scalar, Cols>::array_and_uof_type errors = Hist<Scalar, Cols>::array_and_uof_type::Zero(nbins_and_uof);

            for (auto i = 0u; i < nbins_and_uof; i++) {
                edges(i) = h->GetBinLowEdge(i);
                contents(i) = h->GetBinContent(i);
                errors(i) = h->GetBinError(i);
            }
            edges(nedges-1) = h->GetBinLowEdge(h->GetNbinsX() + 2);

            return Hist<Scalar, Cols>(std::move(contents),
                                      std::move(edges),
                                      std::move(errors),
                                      exposure);
        }
    }

    /////////////////////////////////////////////////////////
    /// \brief Ensure that the rhs histogram has binning
    /// consistent with this histogram.
    /// Raises
    template<typename Scalar, int Cols>
    void
    Hist<Scalar, Cols>::
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
    template<typename Scalar, int Cols>
    Hist<Scalar, Cols>
    Hist<Scalar, Cols>::
    TrueDivide(const Hist & rhs) const {
        EnsureConsistentBinning(rhs, __PRETTY_FUNCTION__);;

        return Hist<Scalar, Cols>(fContentsAndUOF / rhs.ContentsAndUOF(),
                                  fEdgesAndUOF,
                                  (fErrorsAndUOF.pow(2) + rhs.fErrorsAndUOF.pow(2)).sqrt(),
                                  fExposure);
    }

    /////////////////////////////////////////////////////////
    template<typename Scalar, int Cols>
    Hist<Scalar, Cols>
    Hist<Scalar, Cols>::
    ScaleByExposure(Scalar new_expo) const {
        return Hist<Scalar, Cols>(fContentsAndUOF * (new_expo / fExposure),
                                  fEdgesAndUOF,
                                  fErrorsAndUOF * (new_expo / fExposure),
                                  new_expo);
    }

    /////////////////////////////////////////////////////////
    template<typename Scalar, int Cols>
    Scalar
    Hist<Scalar, Cols>::
    Integrate() const {
        return fContentsAndUOF.sum();
    }

    /////////////////////////////////////////////////////////
    template<typename Scalar, int Cols>
    Hist<Scalar, Cols>::
    Hist(const int & nbins,
         const Scalar & min,
         const Scalar & max) {
        fContentsAndUOF = array_and_uof_type::Zero(nbins+2);
        fErrorsAndUOF = array_and_uof_type::Zero(nbins+2);

        auto step = (max - min) / nbins;
        fEdgesAndUOF = edges_type::LinSpaced(nbins + 3,
                                             min - step,
                                             max + step);
    }

    /////////////////////////////////////////////////////////
    template<class Scalar,
            int Cols>
    void
    Hist<Scalar, Cols>::
    SaveTo(TDirectory * dir, const std::string & subdir) const {
        TDirectory * tmp = gDirectory;
        dir->mkdir(subdir.c_str());
        dir = dir->GetDirectory(subdir.c_str());
        dir->cd();

        auto h = root::ToTH1(*this);
        h->Write("hist");

        // exposure is saved in a histogram so it gets accumulated with ROOT's hadd
        TH1 * exposure;
        if constexpr(std::is_same<Scalar, double>::value) {
            exposure = new TH1D("", "", 1, 0, 1);
        } else if constexpr(std::is_same<Scalar, float>::value) {
            exposure = new TH1F("", "", 1, 0, 1);
        }
        exposure->SetBinContent(1, fExposure);
        exposure->Write("exposure");

        tmp->cd();
    }

    /////////////////////////////////////////////////////////
    template<class Scalar,
            int Cols>
    std::unique_ptr<Hist<Scalar, Cols> >
    Hist<Scalar, Cols>::
    LoadFrom(TDirectory * dir, const std::string & subdir) {
        dir = dir->GetDirectory(subdir.c_str());
        dir->cd();

        TH1 * h;
        TH1 * exposure;
        if constexpr(std::is_same<Scalar, double>::value) {
            h = (TH1D *) dir->Get("hist");
            exposure = (TH1D *) dir->Get("exposure");
        } else if constexpr(std::is_same<Scalar, float>::value) {
            h = (TH1F *) dir->Get("hist");
            exposure = (TH1F *) dir->Get("exposure");
        } else {
            std::cerr << "Hist::LoadFrom not implemented for this type" << std::endl;
            exit(1);
        }

        if (!h) {
            std::cerr << "Object TH1 was not found in " << dir->GetPath() << std::endl;
            exit(1);
        }

        return std::make_unique<Hist<Scalar, Cols> >(root::FromTH1<Scalar, Cols>(h, exposure->GetBinContent(1)));
    }

    /////////////////////////////////////////////////////////
    template<class Scalar,
            int Cols>
    bool
    Hist<Scalar, Cols>::
    operator==(const Hist<Scalar, Cols> & rhs) const {
        auto scale_exposure = this->fExposure / rhs.fExposure;
        return (this->fContentsAndUOF - rhs.fContentsAndUOF * scale_exposure).isZero(0) &&
               (this->fEdgesAndUOF - rhs.fEdgesAndUOF).isZero(0) &&
               (this->fErrorsAndUOF - rhs.fErrorsAndUOF * scale_exposure).isZero(0) &&
               (this->fExposure == rhs.fExposure);
    }

    /////////////////////////////////////////////////////////
    /// \brief Divide bin contents by corresponding bin widths in place.
    /// not including over/underflow
    template<class Scalar,
            int Cols>
    Hist<Scalar, Cols> &
    Hist<Scalar, Cols>::
    BinWidthNormalize() {
        fContentsAndUOF(Eigen::seqN(1, fContentsAndUOF.size() - 2)) /= this->BinWidths();
        return *this;
    }

    /////////////////////////////////////////////////////////
    /// \brief Return new histogram with bin contents
    /// divided by corresponding bin widths.
    /// not including over/underflow
    template<class Scalar,
            int Cols>
    Hist<Scalar, Cols>
    Hist<Scalar, Cols>::
    BinWidthNormalize() const {
        auto copy = *this;
        copy = copy.BinWidthNormalize();
        return copy;
    }

    /////////////////////////////////////////////////////////
    /// \brief Divide bin contents and under/overflow by
    /// sum of bin contents in place
    template<class Scalar,
            int Cols>
    Hist<Scalar, Cols> &
    Hist<Scalar, Cols>::
    AreaNormalize() {
        fContentsAndUOF /= this->Contents().sum();
        return *this;
    }

    /////////////////////////////////////////////////////////
    /// \brief Divide bin contents and under/overflow by
    /// sum of bin contents in place
    template<class Scalar,
            int Cols>
    Hist<Scalar, Cols>
    Hist<Scalar, Cols>::
    AreaNormalize() const {
        auto copy = *this;
        copy = copy.AreaNormalize();
        return copy;
    }

    /////////////////////////////////////////////////////////
    /// \brief Return an array of bin widths
    /// not including over/underflow
    template<class Scalar,
            int Cols>
    Eigen::Array<Scalar, 1, Cols>
    Hist<Scalar, Cols>::
    BinWidths() const {
        return fEdgesAndUOF(Eigen::seqN(3, fContentsAndUOF.size()-2)) -
               fEdgesAndUOF(Eigen::seqN(2, fContentsAndUOF.size()-2));
    }

    /////////////////////////////////////////////////////////
    /// \brief Return a new histogram who's contents and errors
    /// are the absolute values of this histogram
    template<typename Scalar,
            int Cols>
    Hist<Scalar, Cols>
    Hist<Scalar, Cols>::
    abs() const {
        return Hist<Scalar, Cols>(this->fContentsAndUOF.abs(),
                                  this->fEdgesAndUOF,
                                  this->fErrorsAndUOF.abs(),
                                  this->fExposure);
    }

    /////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////
    /// \brief Return a new histogram who's contents and errors
    /// are the absolute values squared of this histogram
    template<typename Scalar,
            int Cols>
    Hist<Scalar, Cols>
    Hist<Scalar, Cols>::
    abs2() const {
        return Hist<Scalar, Cols>(this->fContentsAndUOF.abs2(),
                                  this->fEdgesAndUOF,
                                  this->fErrorsAndUOF.abs2(),
                                  this->fExposure);
    }

    /////////////////////////////////////////////////////////
    /// \brief Return a new histogram who's contents
    /// are the square root of this histogram
    /// and errors are divided by two
    template<typename Scalar,
            int Cols>
    Hist<Scalar, Cols>
    Hist<Scalar, Cols>::
    sqrt() const {
        return Hist<Scalar, Cols>(this->fContentsAndUOF.sqrt(),
                                  this->fEdgesAndUOF,
                                  this->fErrorsAndUOF / 2.,
                                  std::sqrt(this->fExposure));
    }

    /////////////////////////////////////////////////////////
    /// \brief Return a new histogram who's contents
    /// are raised to the exp power of this histogram
    /// and errors are multiplied by exp
    template<typename Scalar,
            int Cols>
    Hist<Scalar, Cols>
    Hist<Scalar, Cols>::
    pow(Scalar exp) const {
        return Hist<Scalar, Cols>(this->fContentsAndUOF.pow(exp),
                                  this->fEdgesAndUOF,
                                  this->fErrorsAndUOF * exp,
                                  std::pow(this->fExposure, exp));
    }

    /////////////////////////////////////////////////////////
    template<typename Scalar,
            int Cols>
    Hist<Scalar, Cols>
    Hist<Scalar, Cols>::
    operator-(const Hist & rhs) const {
        Hist<Scalar, Cols> ret = *this; // copy this
        ret -= rhs;
        return ret;
    }

    /////////////////////////////////////////////////////////
    template<typename Scalar,
            int Cols>
    Hist<Scalar, Cols>
    Hist<Scalar, Cols>::
    operator-=(const Hist & rhs) {
        EnsureConsistentBinning(rhs, __PRETTY_FUNCTION__);;

        this->fContentsAndUOF -= rhs.fContentsAndUOF * fExposure / rhs.fExposure;
        this->fErrorsAndUOF = (this->fErrorsAndUOF.pow(2) +
                               (rhs.fErrorsAndUOF * fExposure / rhs.fExposure).pow(2)).sqrt();
        return *this;
    }

    /////////////////////////////////////////////////////////
    template<typename Scalar,
            int Cols>
    Hist<Scalar, Cols>
    Hist<Scalar, Cols>::
    operator+(const Hist & rhs) const {
        Hist<Scalar, Cols> ret = *this; // copy this
        ret += rhs;
        return ret;
    }

    /////////////////////////////////////////////////////////
    template<typename Scalar,
            int Cols>
    Hist<Scalar, Cols>
    Hist<Scalar, Cols>::
    operator+=(const Hist & rhs) {
        EnsureConsistentBinning(rhs, __PRETTY_FUNCTION__);;

        this->fContentsAndUOF += rhs.fContentsAndUOF * fExposure / rhs.fExposure;
        this->fErrorsAndUOF = (this->fErrorsAndUOF.pow(2) +
                               (rhs.fErrorsAndUOF * fExposure / rhs.fExposure).pow(2)).sqrt();
        return *this;
    }

    /////////////////////////////////////////////////////////
    template<typename Scalar,
            int Cols>
    Hist<Scalar, Cols>
    Hist<Scalar, Cols>::
    operator*(const Hist & rhs) const {
        Hist<Scalar, Cols> ret = *this; // copy this
        ret *= rhs;
        return ret;
    }

    /////////////////////////////////////////////////////////
    template<typename Scalar,
            int Cols>
    Hist<Scalar, Cols>
    Hist<Scalar, Cols>::
    operator*=(const Hist & rhs) {
        EnsureConsistentBinning(rhs, __PRETTY_FUNCTION__);;

        this->fContentsAndUOF *= rhs.fContentsAndUOF * fExposure / rhs.fExposure;
        this->fErrorsAndUOF = (this->fErrorsAndUOF.pow(2) +
                               (rhs.fErrorsAndUOF * fExposure / rhs.fExposure).pow(2)).sqrt();
        return *this;
    }

    /////////////////////////////////////////////////////////
    template<typename Scalar,
            int Cols>
    Hist<Scalar, Cols>
    Hist<Scalar, Cols>::
    operator/(const Hist & rhs) const {
        Hist<Scalar, Cols> ret = *this; // copy this
        ret /= rhs;
        return ret;
    }

    /////////////////////////////////////////////////////////
    template<typename Scalar,
            int Cols>
    Hist<Scalar, Cols>
    Hist<Scalar, Cols>::
    operator/=(const Hist & rhs) {
        EnsureConsistentBinning(rhs, __PRETTY_FUNCTION__);;

        this->fContentsAndUOF /= rhs.fContentsAndUOF * fExposure / rhs.fExposure;
        this->fErrorsAndUOF = (this->fErrorsAndUOF.pow(2) +
                               (rhs.fErrorsAndUOF * fExposure / rhs.fExposure).pow(2)).sqrt();
        this->fExposure = 1;
        return *this;
    }

    /////////////////////////////////////////////////////////
    template<typename Scalar,
            int Cols>
    Hist<Scalar, Cols>
    Hist<Scalar, Cols>::
    operator-(const Scalar & rhs) const {
        Hist<Scalar, Cols> ret = *this; // copy this
        ret -= rhs;
        return ret;
    }

    /////////////////////////////////////////////////////////
    template<typename Scalar,
            int Cols>
    Hist<Scalar, Cols>
    Hist<Scalar, Cols>::
    operator-=(const Scalar & rhs) {
        this->fContentsAndUOF -= rhs;
        return *this;
    }

    /////////////////////////////////////////////////////////
    template<typename Scalar,
            int Cols>
    Hist<Scalar, Cols>
    Hist<Scalar, Cols>::
    operator+(const Scalar & rhs) const {
        Hist<Scalar, Cols> ret = *this; // copy this
        ret += rhs;
        return ret;
    }

    /////////////////////////////////////////////////////////
    template<typename Scalar,
            int Cols>
    Hist<Scalar, Cols>
    Hist<Scalar, Cols>::
    operator+=(const Scalar & rhs) {
        this->fContentsAndUOF += rhs;
        return *this;
    }

    /////////////////////////////////////////////////////////
    template<typename Scalar,
            int Cols>
    Hist<Scalar, Cols>
    Hist<Scalar, Cols>::
    operator*(const Scalar & rhs) const {
        Hist<Scalar, Cols> ret = *this; // copy this
        ret *= rhs;
        return ret;
    }

    /////////////////////////////////////////////////////////
    template<typename Scalar,
            int Cols>
    Hist<Scalar, Cols>
    Hist<Scalar, Cols>::
    operator/(const Scalar & rhs) const {
        Hist<Scalar, Cols> ret = *this; // copy this
        ret /= rhs;
        return ret;
    }

    /////////////////////////////////////////////////////////
    template<typename Scalar,
            int Cols>
    Hist<Scalar, Cols>
    Hist<Scalar, Cols>::
    operator/=(const Scalar & rhs) {
        this->fContentsAndUOF /= rhs;
        return *this;
    }


    /////////////////////////////////////////////////////////
    template<typename Scalar,
            int Cols>
    Hist<Scalar, Cols>
    Hist<Scalar, Cols>::
    operator*=(const Scalar & rhs) {
        this->fContentsAndUOF *= rhs;
        return *this;
    }

    /////////////////////////////////////////////////////////
    template<typename Scalar,
            int Cols>
    Hist<Scalar, Cols>
    Hist<Scalar, Cols>::
    operator=(const Hist<Scalar, Cols> & rhs) {
        if (this == &rhs) return *this;
        fContentsAndUOF = rhs.fContentsAndUOF;
        fEdgesAndUOF = rhs.fEdgesAndUOF;
        fErrorsAndUOF = rhs.fErrorsAndUOF;
        fExposure = rhs.fExposure;
        return *this;
    }

    /////////////////////////////////////////////////////////
    template<typename Scalar,
            int Cols>
    Hist<Scalar, Cols>
    Hist<Scalar, Cols>::
    operator=(Hist<Scalar, Cols> && rhs) {
        if(this == &rhs) return *this;
        fContentsAndUOF = std::move(rhs.fContentsAndUOF);
        fEdgesAndUOF = std::move(rhs.fEdgesAndUOF);
        fErrorsAndUOF = std::move(rhs.fErrorsAndUOF);
        fExposure = rhs.fExposure;
        return *this;
    }
}
