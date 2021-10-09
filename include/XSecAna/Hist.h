#pragma once

#include <memory>
#include <Eigen/Dense>

#include "TDirectory.h"
#include "TObjString.h"
#include "TH1D.h"
#include "TH1F.h"

#include <iostream>
#include <exception>

#include "XSecAna/Array.h"
#include "XSecAna/_Hist.h"
//#include "XSecAna/ROOT/ROOTInterface.h"

namespace xsec {

    /// \brief Object representing a filled histogram used internally by the framework
    /// Bin contents are indexed starting at 1
    /// Underflow is index 0
    /// Overflow is index N+1 where N is number of bins
    /// Wraps Eigen arrays for fast mathematical operations especially linear algebra
    ///
    /// TODO Users can create conversion functions to/from this object
    class Hist : public _hist {
    public:
        Hist() {}
        /*************** ctors ************************/

        /// \brief ctor using existing arrays for contents and edges
        /// including under/overflow
        Hist(const Array & contents_and_uof,
             const Array & edges_and_uof,
             const Array & errors_and_uof,
             const double & exposure = 1);

        Hist(const WeightedArray & contents_and_uof,
             const Array & edges_and_uof,
             const Array & errors_and_uof);

        /// \brief ctor using existing arrays for contents and edges
        /// including under/overflow
        /// Initializes errors to 0
        Hist(const Array & contents_and_uof,
             const Array & edges_and_uof,
             const double & exposure = 1);

        Hist(const WeightedArray & contents_and_uof,
             const Array & edges_and_uof);

        /// \brief convenience constructor.
        /// works for dynamic and fixed-size histograms
        /// adds two additional "bins" for under/overflow
        /// to the underlaying arrays
        Hist(const int & nbins,
             const double & min,
             const double & max,
             const double & exposure = 1);

        Hist(const Hist & rhs);


        /*************** hist1d public functions ************************/
        static Hist * FromTH1(const TH1 * h, const double & exposure);
        const Array GetEdges() const;
        const Array GetEdgesAndUOF() const;
        virtual Array GetBinWidths() const;
        double & operator()(int index);
        double operator()(int index) const;
        static std::unique_ptr<_hist> LoadFrom(TDirectory * dir, const std::string & subdir);

        /*************** virtual public functions inherited from _hist ************************/
        TH1 * ToTH1(const std::string & name = "", const std::string & title = "") const override;
        int GetDimensions() const override { return 1; }
        Array GetContentsAndUOF() const override;
        Array GetContents() const override;
        Array GetErrorsAndUOF() const override;
        Array GetErrors() const override;
        void SetContentsAndUOF(const Array & contents_and_uof) override;
        void SetErrorsAndUOF(const Array & errors_and_uof) override;
        void SetContents(const Array & contents) override;
        void SetErrors(const Array & errors) override;
        void SetExposure(const double & new_exposure) override;

        virtual double Exposure() const override;

        virtual _hist * Clone() const override;

        virtual double Integrate() const override;

    private:
        /*************** virtual private functions inherited from _hist ************************/
        _hist * _subtract(const _hist * rhs) override;

        _hist * _subtract(const double & rhs) override;

        _hist * _add(const _hist * rhs) override;

        _hist * _add(const double & rhs) override;

        _hist * _divide(const _hist * rhs) override;

        _hist * _divide(const double & rhs) override;

        _hist * _multiply(const _hist * rhs) override;

        _hist * _multiply(const double & rhs) override;

        _hist * _true_divide(const _hist * rhs) override;

        _hist * _scale_by_exposure(const double & new_expo) override;

        bool _is_same_binning(const _hist * rhs, const double & tol) const override;

        bool _is_same_contents(const _hist * rhs, const double & tol) const override;

        _hist * _abs() override;

        _hist * _abs2() override;

        _hist * _sqrt() override;

        _hist * _pow(double exp) override;

        _hist * _bin_width_normalize() override;

        _hist * _area_normalize() override;

        WeightedArray fContentsAndUOF;
        Array fErrorsAndUOF;
        Array fEdgesAndUOF;
    };
}
