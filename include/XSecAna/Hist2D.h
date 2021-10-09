#pragma once

#include "XSecAna/Hist.h"
#include "XSecAna/Array2D.h"
#include "XSecAna/Hist.h"
#include <Eigen/Dense>

namespace xsec {

    class Hist2D : public Hist {
    public:
        Hist2D() = default;

        Hist2D(int nx, double xmin, double xmax,
               int ny, double ymin, double ymax,
               double exposure = 1);

        Hist2D(const Array2D contents_and_uof,
               const Array xedges_and_uof,
               const Array yedges_and_uof,
               const Array2D errors_and_uof,
               const double & exposure = 1);

        Hist2D(const WeightedArray2D contents_and_uof,
               const Array xedges_and_uof,
               const Array yedges_and_uof,
               const Array2D errors_and_uof);

        Hist2D(const Array2D & contents_and_uof,
               const Array & xedges_and_uof,
               const Array & yedges_and_uof,
               const double & exposure = 1);

        Hist2D(const WeightedArray2D & contents_and_uof,
               const Array & xedges_and_uof,
               const Array & yedges_and_uof);

        Hist2D(const Hist2D & rhs);

        Hist2D(const Hist2D && rhs);

        Hist2D BinWidthNormalize() const;

        virtual Hist2D & BinWidthNormalize();

        virtual Hist2D AreaNormalize() const;

        virtual Hist2D & AreaNormalize();

        virtual Hist2D operator-(const Hist2D & rhs) const;

        virtual Hist2D operator+(const Hist2D & rhs) const;

        virtual Hist2D operator/(const Hist2D & rhs) const;

        virtual Hist2D operator*(const Hist2D & rhs) const;

        virtual bool operator==(const Hist2D & rhs) const;

        virtual bool operator!=(const Hist2D & rhs) const { return !(*this == rhs); }

        virtual Hist2D operator-=(const Hist2D & rhs);

        virtual Hist2D operator+=(const Hist2D & rhs);

        virtual Hist2D operator/=(const Hist2D & rhs);

        virtual Hist2D operator*=(const Hist2D & rhs);

        virtual Hist2D operator-(const double & rhs) const;

        virtual Hist2D operator+(const double & rhs) const;

        virtual Hist2D operator/(const double & rhs) const;

        virtual Hist2D operator*(const double & rhs) const;

        virtual Hist2D operator-=(const double & rhs);

        virtual Hist2D operator+=(const double & rhs);

        virtual Hist2D operator/=(const double & rhs);

        virtual Hist2D operator*=(const double & rhs);

        virtual Hist2D operator=(const Hist2D & rhs);

        virtual Hist2D & operator=(Hist2D && rhs);

        virtual Hist2D ScaleByExposure(double new_expo) const;

        virtual Hist2D TrueDivide(const Hist2D & rhs) const;

        // some convenience functions
        virtual Hist2D abs() const;

        virtual Hist2D abs2() const;

        virtual Hist2D sqrt() const;

        virtual Hist2D pow(double exp) const;

        virtual const Array2D Contents() const {
            return fContentsAndUOF.array().block(1, 1,
                                                 fContentsAndUOF.rows()-1,
                                                 fContentsAndUOF.cols()-1);
        }

        virtual const Array2D ContentsAndUOF() const { return fContentsAndUOF.array(); }

        virtual const Array2D Errors() const { return fErrorsAndUOF.array()(Eigen::seq(1, fErrorsAndUOF.size() - 2)); }

        virtual const Array2D ErrorsAndUOF() const { return fErrorsAndUOF; }

        virtual const Array XEdges() const { return fXEdgesAndUOF(Eigen::seq(1, fXEdgesAndUOF.size() - 2)); }
        virtual const Array XEdges() const { return fYEdgesAndUOF(Eigen::seq(1, fYEdgesAndUOF.size() - 2)); }

        virtual const Array XEdgesAndUOF() const { return fEdgesAndUOF; }

        void SetContentsAndUOF(const Array2D & contents_and_uof) {
            assert(fContentsAndUOF.rows() == contents_and_uof.rows() &&
                   fContentsAndUOF.cols() == contents_and_uof.cols() &&
                   "Incompatible contents array");
            fContentsAndUOF = WeightedArray2D(contents_and_uof,
                                              this->Exposure());
        }

        void SetContents(const Array2D & contents) {
            assert(fContentsAndUOF.rows()-2 == contents.rows() &&
                   fContentsAndUOF.cols()-2 == contents.cols() &&
                   "Incompatible contents array");
            fContentsAndUOF.array().block(1, 1,
                                          contents.rows(),
                                          contents.cols())
                    = contents;
        }

        virtual void SetErrorsAndUOF(const Array & errors_and_uof) {
            assert(fErrorsAndUOF.rows() == errors_and_uof.rows() &&
                   fErrorsAndUOF.cols() == errors_and_uof.cols() &&
                   "Incompatible contents array");
            fErrorsAndUOF = errors_and_uof;
        }

        virtual double Exposure() const override { return fContentsAndUOF.weight(); }

        Array XBinWidths() const;
        Array YBinWidths() const;

        virtual double & operator()(int i, int j) { return fContentsAndUOF.array()(i, j); }

        virtual double operator()(int i, int j) const { return fContentsAndUOF.array()(i, j); }

        virtual void SaveTo(TDirectory * dir, const std::string & subdir) const;

        static std::unique_ptr<Hist> LoadFrom(TDirectory * dir, const std::string & subdir);

    private:
        WeightedArray2D fContentsAndUOF;
        Array2D fErrorsAndUOF;
        Array fXEdgesAndUOF;
        Array fYEdgesAndUOF;
    };

    Hist2D::
    Hist2D(int nx, xmin, xmax,
           int ny, ymin, ymax,
           double exposure) {
        fContentsAndUOF = WeightedArray2D(Array2D::Zero((nx + 2),
                                                        (ny + 2)),
                                          exposure);
        fErrorsAndUOF = Array2D::Zero(nx + 2,
                                      ny + 2);

        auto xstep = (xmax - xmin) / nx;
        fXEdgesAndUOF = Array::LinSpaced(nx + 3,
                                         xmin - xstep,
                                         xmax + xstep);
        auto ystep = (ymax - ymin) / ny;
        fYEdgesAndUOF = Array::LinSpaced(ny + 3,
                                         ymin - ystep,
                                         ymax + ystep);
    }

    Hist2D::
    Hist2D(const Array2D contents_and_uof,
           const Array xedges_and_uof,
           const Array yedges_and_uof,
           const Array2D errors_and_uof,
           const double & exposure = 1)
            : fContentsAndUOF(contents_and_uof,
                              exposure),
              fErrorsAndUOF(errors_and_uof),
              fXEdgesAndUOF(xedges_and_uof),
              fYEdgesAndUOF(yedges_and_uof) {}

    Hist2D::
    Hist2D(const WeightedArray2D contents_and_uof,
           const Array xedges_and_uof,
           const Array yedges_and_uof,
           const Array2D errors_and_uof)
            : fContentsAndUOF(contents_and_uof),

              fErrorsAndUOF(errors_and_uof),
              fXEdgesAndUOF(xedges_and_uof),
              fYEdgesAndUOF(yedges_and_uof) {
        assert(contents_and_uof.rows() + 1 == xedges_and_uof.size() &&
               contents_and_uof.cols() + 1 == yedges_and_uof.size() &&
               "Incompatible edges, contents, and/or errors");
    }

    Hist2D::
    Hist2D(const Array2D & contents_and_uof,
           const Array & xedges_and_uof,
           const Array & yedges_and_uof,
           const double & exposure = 1)
            : fContentsAndUOF(contents_and_uof, exposure),
              fErrorsAndUOF(Array2D::Zero(fContentsAndUOF.rows(),
                                          fContentsAndUOF.cols())),
              fXEdgesAndUOF(xedges_and_uof),
              fYEdgesAndUOF(yedges_and_uof) {
        assert(contents_and_uof.rows() + 1 == xedges_and_uof.size() &&
               contents_and_uof.cols() + 1 == yedges_and_uof.size() &&
               "Incompatible edges, contents, and/or errors");
    }

    Hist2D::
    Hist2D(const WeightedArray2D & contents_and_uof,
           const Array & xedges_and_uof,
           const Array & yedges_and_uof)
            : fContentsAndUOF(contents_and_uof),
              fErrorsAndUOF(Array2D::Zero(fContentsAndUOF.rows(),
                                          fContentsAndUOF.cols())),
              fXEdgesAndUOF(xedges_and_uof),
              fYEdgesAndUOF(yedges_and_uof) {
        assert(contents_and_uof.rows() + 1 == xedges_and_uof.size() &&
               contents_and_uof.cols() + 1 == yedges_and_uof.size() &&
               "Incompatible edges, contents, and/or errors");
    }

    Hist2D::
    Hist2D(const Hist2D & rhs)
            : fContentsAndUOF(rhs.fContentsAndUOF),
              fErrorsAndUOF(rhs.fContentsAndUOF),
              fXEdgesAndUOF(rhs.fXEdgesAndUOF),
              fYEdgesAndUOF(rhs.fYEdgesAndUOF) {}

    Hist2D::
    Hist2D(const Hist2D && rhs)
            : fContentsAndUOF(std::move(rhs.fContentsAndUOF)),
              fErrorsAndUOF(std::move(rhs.fContentsAndUOF)),
              fXEdgesAndUOF(std::move(rhs.fXEdgesAndUOF)),
              fYEdgesAndUOF(std::move(rhs.fYEdgesAndUOF)) {}
    Hist2D
    Hist2D::
    BinWidthNormalize() const {
        auto copy = *this;
        copy = copy.BinWidthNormalize();
        return copy;
    }

    Hist2D &
    Hist2D:: BinWidthNormalize() {
        auto temp = this->Contents();
        temp = this->BinWidthX() * temp * this->BinWidthsY();
        this->SetContents() = temp;
        return *this;
    }

    Hist2D
    Hist2D::
    AreaNormalize() const {
        auto copy = *this;
        copy = copy.AreaNormalize();
        return copy;
    }

    Hist2D &
    Hist2D::
    AreaNormalize() {
        this->SetContentsAndUOF(this->GetContentsAndUOF() / this->Integrate());
    }

    Hist2D
    Hist2D::
    operator-(const Hist2D & rhs) const {

    }

    Hist2D
    Hist2D::
    operator+(const Hist2D & rhs) const {

    }

    Hist2D
    Hist2D::
    operator/(const Hist2D & rhs) const {

    }

    Hist2D
    Hist2D::
    operator*(const Hist2D & rhs) const {

    }

    bool
    Hist2D::
    operator==(const Hist2D & rhs) const {

    }


    Hist2D
    Hist2D::
    operator-=(const Hist2D & rhs);

    Hist2D
    Hist2D::
    operator+=(const Hist2D & rhs);

    Hist2D
    Hist2D::
    operator/=(const Hist2D & rhs);

    Hist2D
    Hist2D::
    operator*=(const Hist2D & rhs);

    Hist2D
    Hist2D::
    operator-(const double & rhs) const;

    Hist2D
    Hist2D::
    operator+(const double & rhs) const;

    Hist2D
    Hist2D::
    operator/(const double & rhs) const;

    Hist2D
    Hist2D::
    operator*(const double & rhs) const;

    Hist2D
    Hist2D::
    operator-=(const double & rhs);

    Hist2D
    Hist2D::
    operator+=(const double & rhs);

    Hist2D
    Hist2D::
    operator/=(const double & rhs);

    Hist2D
    Hist2D::
    operator*=(const double & rhs);

    Hist2D
    Hist2D::
    operator=(const Hist2D & rhs);

    Hist2D &
    Hist2D::
    operator=(Hist2D && rhs);

    Hist2D
    Hist2D::
    ScaleByExposure(double new_expo) const;

    Hist2D
    Hist2D::
    TrueDivide(const Hist2D & rhs) const;

    // some convenience functions
    Hist2D
    Hist2D::
    abs() const;

    Hist2D
    Hist2D::
    abs2() const;

    Hist2D
    Hist2D::
    sqrt() const;

    Hist2D
    Hist2D::
    pow(double exp) const;

    Array
    Hist2D::
    XBinWidths() const {
        return fXEdgesAndUOF(Eigen::seqN(3, fContentsAndUOF.rows() - 2)) -
               fXEdgesAndUOF(Eigen::seqN(2, fContentsAndUOF.rows() - 2));
    }

    Array
    Hist2D::
    YBinWidths() const {
        return fYEdgesAndUOF(Eigen::seqN(3, fContentsAndUOF.cols() - 2)) -
               fYEdgesAndUOF(Eigen::seqN(2, fContentsAndUOF.cols() - 2));
    }

    void
    Hist2D::
    SaveTo(TDirectory * dir, const std::string & subdir) const;

    std::unique_ptr<Hist> LoadFrom(TDirectory * dir, const std::string & subdir);

}

