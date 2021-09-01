#pragma once

#include <memory>
#include <Eigen/Dense>
#include "TDirectory.h"
#include "TObjString.h"
#include "TH1D.h"
#include "TH1F.h"

#include <iostream>

namespace xsec {
    // compile time expression for determining size of Eigen Array holding
    // bin edges
    constexpr int EdgesSize(int Cols) { return Cols == Eigen::Dynamic ? Eigen::Dynamic : Cols + 1; }

    // Histograming object used internally by the framework
    // Wraps Eigen arrays
    // Users can create conversion functions to/from this object
    // template parameters are forwarded to the underlaying Eigen arrays
    template<typename Scalar,
            int Cols = -1>
    class Hist {
    public:

        typedef Scalar scalar_type;
        typedef Eigen::Array<Scalar, 1, Cols> array_type;
        typedef Eigen::Array<Scalar, 1, EdgesSize(Cols)> edges_type;

        Hist() {}

        Hist(const Eigen::Array<Scalar, 1, Cols> & contents,
             const Eigen::Array<Scalar, 1, EdgesSize(Cols)> & edges)
                : fContents(contents), fEdges(edges), fExposure(1) {}

        Hist(const Eigen::Array<Scalar, 1, Cols> & contents,
             const Eigen::Array<Scalar, 1, EdgesSize(Cols)> & edges,
             Scalar exposure)
                : fContents(contents), fEdges(edges), fExposure(exposure) {}


        // convenience constructor.
        // works for dynamic and fixed-size histograms
        Hist(const int & nbins,
             const Scalar & min,
             const Scalar & max);

        void Normalize(const std::string & how);

        Scalar Integrate() const;

        Hist operator-(const Hist & rhs) const;

        Hist operator+(const Hist & rhs) const;

        Hist operator/(const Hist & rhs) const;

        Hist operator*(const Hist & rhs) const;

        bool operator==(const Hist & rhs) const;

        bool operator!=(const Hist & rhs) const { return !(*this == rhs); }

        Hist operator-=(const Hist & rhs);

        Hist operator+=(const Hist & rhs);

        Hist operator/=(const Hist & rhs);

        Hist operator*=(const Hist & rhs);

        Hist operator-(const Scalar & rhs) const;

        Hist operator+(const Scalar & rhs) const;

        Hist operator/(const Scalar & rhs) const;

        Hist operator*(const Scalar & rhs) const;

        Hist operator-=(const Scalar & rhs);

        Hist operator+=(const Scalar & rhs);

        Hist operator/=(const Scalar & rhs);

        Hist operator*=(const Scalar & rhs);

        Hist ScaleByExposure(Scalar new_expo) const;

        Hist TrueDivide(const Hist & rhs) const;

        // some convenience functions
        Hist abs() const;

        Hist abs2() const;

        Hist sqrt() const;

        Hist pow(Scalar exp) const;

        const Eigen::Array<Scalar, 1, Cols> & Contents() const { return fContents; }

        const Eigen::Array<Scalar, 1, EdgesSize(Cols)> & Edges() const { return fEdges; }

        Scalar Exposure() const { return fExposure; }

        [[nodiscard]] unsigned int size() const { return fContents.size(); }

        Eigen::Array<Scalar, 1, Cols> BinWidths() const;

        Scalar & operator[](int index) { return fContents(index); }

        Scalar operator[](int index) const { return fContents(index); }

        void SaveTo(TDirectory * dir, const std::string & subdir) const;

        static std::unique_ptr<Hist> LoadFrom(TDirectory * dir, const std::string & subdir);

    private:
        Eigen::Array<Scalar, 1, Cols> fContents;
        Eigen::Array<Scalar, 1, EdgesSize(Cols)> fEdges;
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
                h = new TH1D(name.c_str(), title.c_str(), hist.size(), hist.Edges().data());
            } else if constexpr(std::is_same<Scalar, float>::value) {
                h = new TH1F(name.c_str(), title.c_str(), hist.size(), hist.Edges().data());
            } else {
                std::cerr << "ToTH1 not implemented for this type" << std::endl;
                exit(1);
            }

            for (auto i = 0u; i < hist.size(); i++) {
                h->SetBinContent(i + 1, hist.Contents()(i));
            }

            return h;
        }

        template<typename Scalar, int Cols>
        Hist<Scalar, Cols>
        FromTH1(const TH1 * h, Scalar expo) {
            Scalar edges[h->GetNbinsX() + 1];
            Scalar contents[h->GetNbinsX()];
            for (auto i = 0; i < h->GetNbinsX(); i++) {
                edges[i] = h->GetBinLowEdge(i + 1);
                contents[i] = h->GetBinContent(i + 1);
            }
            edges[h->GetNbinsX()] = h->GetBinLowEdge(h->GetNbinsX() + 1);

            // TODO this won't work for dynamic-sized arrays
            // Probably lots of work needs to be done to handle this case.
            return Hist<Scalar, Cols>(Eigen::Map<Eigen::Array<Scalar, 1, Cols> >(contents, h->GetNbinsX()),
                                      Eigen::Map<Eigen::Array<Scalar, 1, EdgesSize(Cols)> >(edges,
                                                                                            h->GetNbinsX() + 1),
                                      expo);
        }
    }

    /////////////////////////////////////////////////////////
    template<typename Scalar, int Cols>
    Hist<Scalar, Cols>
    Hist<Scalar, Cols>::
    TrueDivide(const Hist & rhs) const {
        return Hist<Scalar, Cols>(fContents / rhs.Contents(),
                                  fEdges,
                                  fExposure);
    }

    /////////////////////////////////////////////////////////
    template<typename Scalar, int Cols>
    Hist<Scalar, Cols>
    Hist<Scalar, Cols>::
    ScaleByExposure(Scalar new_expo) const {
        return Hist<Scalar, Cols>(fContents * (new_expo / fExposure),
                                  fEdges,
                                  new_expo);
    }

    /////////////////////////////////////////////////////////
    template<typename Scalar, int Cols>
    Scalar
    Hist<Scalar, Cols>::
    Integrate() const {
        return fContents.sum();
    }

    /////////////////////////////////////////////////////////
    template<typename Scalar, int Cols>
    Hist<Scalar, Cols>::
    Hist(const int & nbins,
         const Scalar & min,
         const Scalar & max) {
        fContents = Eigen::Array<Scalar, 1, Cols>::Zeros(nbins);
        fEdges = Eigen::Array<Scalar, 1, EdgesSize(Cols)>::LinSpaced(nbins + 1, min, max);
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
        return
                (this->fContents - rhs.fContents).isZero(0) &&
                (this->fEdges - rhs.fEdges).isZero(0),
                        (this->fExposure == rhs.fExposure);
    }

    /////////////////////////////////////////////////////////
    template<class Scalar,
            int Cols>
    void
    Hist<Scalar, Cols>::
    Normalize(const std::string & how) {
        if (how == "width") {
            fContents = fContents / BinWidths();
        } else if (how == "area") {
            fContents /= fContents.sum();
        }
    }

    /////////////////////////////////////////////////////////
    template<class Scalar,
            int Cols>
    Eigen::Array<Scalar, 1, Cols>
    Hist<Scalar, Cols>::
    BinWidths() const {
        return fEdges.tail(fContents.size()) - fEdges.head(fContents.size());
    }

    /////////////////////////////////////////////////////////
    template<typename Scalar,
            int Cols>
    Hist<Scalar, Cols>
    Hist<Scalar, Cols>::
    abs() const {
        return Hist<Scalar, Cols>(this->fContents.abs(),
                                  this->fEdges,
                                  this->fExposure);
    }

    /////////////////////////////////////////////////////////
    template<typename Scalar,
            int Cols>
    Hist<Scalar, Cols>
    Hist<Scalar, Cols>::
    abs2() const {
        return Hist<Scalar, Cols>(this->fContents.abs2(),
                                  this->fEdges,
                                  this->fExposure);
    }

    /////////////////////////////////////////////////////////
    template<typename Scalar,
            int Cols>
    Hist<Scalar, Cols>
    Hist<Scalar, Cols>::
    sqrt() const {
        return Hist<Scalar, Cols>(this->fContents.sqrt(),
                                  this->fEdges,
                                  std::sqrt(this->fExposure));
    }

    /////////////////////////////////////////////////////////
    template<typename Scalar,
            int Cols>
    Hist<Scalar, Cols>
    Hist<Scalar, Cols>::
    pow(Scalar exp) const {
        return Hist<Scalar, Cols>(this->fContents.pow(exp),
                                  this->fEdges,
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
        this->fContents -= rhs.ScaleByExposure(fExposure).fContents;
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
        this->fContents += rhs.ScaleByExposure(fExposure).fContents;
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
        this->fContents *= rhs.ScaleByExposure(fExposure).fContents;
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
        this->fContents /= rhs.ScaleByExposure(fExposure).fContents;
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
        this->fContents -= rhs;
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
        this->fContents += rhs;
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
        this->fContents /= rhs;
        return *this;
    }


    /////////////////////////////////////////////////////////
    template<typename Scalar,
            int Cols>
    Hist<Scalar, Cols>
    Hist<Scalar, Cols>::
    operator*=(const Scalar & rhs) {
        this->fContents *= rhs;
        return *this;
    }

}
