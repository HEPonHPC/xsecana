#pragma once

#include <iostream>
#include <exception>

#include "TDirectory.h"
#include "TH1.h"
#include "XSecAna/Array.h"

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


    class _hist {
    public:
        // non-const final operations
        virtual _hist * Subtract(const _hist * rhs, bool inplace=false) final;
        virtual _hist * Subtract(const double & rhs, bool inplace=false) final;
        virtual _hist * Add(const _hist * rhs, bool inplace=false) final;
        virtual _hist * Add(const double & rhs, bool inplace=false) final;
        virtual _hist * Divide(const _hist * rhs, bool inplace=false) final;
        virtual _hist * Divide(const double & rhs, bool inplace=false) final;
        virtual _hist * Multiply(const _hist * rhs, bool inplace=false) final;
        virtual _hist * Multiply(const double & rhs, bool inplace=false) final;
        virtual _hist * ScaleByExposure(const double & new_expo, bool inplace=false) final;
        virtual _hist * TrueDivide(const _hist * rhs, bool inplace=false) final;
        virtual _hist * abs(bool inplace=false) final;
        virtual _hist * abs2(bool inplace=false) final;
        virtual _hist * sqrt(bool inplace=false) final;
        virtual _hist * pow(double exp, bool inplace=false) final;
        virtual _hist * BinWidthNormalize(bool inplace=false) final;
        virtual _hist * AreaNormalize(bool inplace=false) final;

        // const final operations
        virtual void SaveTo(TDirectory * dir, const std::string & subdir) const final;
        [[nodiscard]] virtual _hist * Subtract(const _hist * rhs) const final;
        [[nodiscard]] virtual _hist * Subtract(const double & rhs) const final;
        [[nodiscard]] virtual _hist * Add(const _hist * rhs) const final;
        [[nodiscard]] virtual _hist * Add(const double & rhs) const final;
        [[nodiscard]] virtual _hist * Divide(const _hist * rhs) const final;
        [[nodiscard]] virtual _hist * Divide(const double & rhs) const final;
        [[nodiscard]] virtual _hist * Multiply(const _hist * rhs) const final;
        [[nodiscard]] virtual _hist * Multiply(const double & rhs) const final;
        [[nodiscard]] virtual _hist * ScaleByExposure(const double & new_expo) const final;
        [[nodiscard]] virtual _hist * TrueDivide(const _hist * rhs) const final;
        [[nodiscard]] virtual _hist * abs() const final;
        [[nodiscard]] virtual _hist * abs2() const final;
        [[nodiscard]] virtual _hist * sqrt() const final;
        [[nodiscard]] virtual _hist * pow(double exp) const final;
        [[nodiscard]] virtual _hist * BinWidthNormalize() const final;
        [[nodiscard]] virtual _hist * AreaNormalize() const final;


        virtual bool IsEqual(const _hist * rhs, const double & tol = 0) const final;

        // must override
        virtual TH1 * ToTH1(const std::string & name = "", const std::string & title = "") const = 0;
        virtual int GetDimensions() const = 0;
        virtual Array GetContents() const = 0;
        virtual Array GetContentsAndUOF() const = 0;
        virtual Array GetErrors() const = 0;
        virtual Array GetErrorsAndUOF() const = 0;
        virtual void SetContents(const Array & contents) = 0;
        virtual void SetContentsAndUOF(const Array & contents_and_uof) = 0;
        virtual void SetErrors(const Array & errors) = 0;
        virtual void SetErrorsAndUOF(const Array & errors_and_uof) = 0;
        virtual void SetExposure(const double & new_exposure) = 0;



        virtual _hist * Clone() const = 0;
        virtual double Exposure() const = 0;
        virtual double Integrate() const = 0;


        static std::unique_ptr<_hist> LoadFrom(TDirectory * dir, const std::string & subdir);
        
        virtual ~_hist()= default;
    private:
        // must override
        virtual _hist * _subtract(const _hist * rhs) = 0;
        virtual _hist * _subtract(const double & rhs) = 0;
        virtual _hist * _add(const _hist * rhs) = 0;
        virtual _hist * _add(const double & rhs) = 0;
        virtual _hist * _divide(const _hist * rhs) = 0;
        virtual _hist * _divide(const double & rhs) = 0;
        virtual _hist * _multiply(const _hist * rhs) = 0;
        virtual _hist * _multiply(const double & rhs) = 0;
        virtual _hist * _scale_by_exposure(const double & new_expo) = 0;
        virtual _hist * _true_divide(const _hist * rhs) = 0;
        virtual bool _is_same_binning(const _hist * rhs, const double & tol) const = 0;
        virtual bool _is_same_contents(const _hist * rhs, const double & tol) const = 0;
        virtual _hist * _abs() = 0;
        virtual _hist * _abs2() = 0;
        virtual _hist * _sqrt() = 0;
        virtual _hist * _pow(double exp) = 0;
        virtual _hist * _bin_width_normalize() = 0;
        virtual _hist * _area_normalize() = 0;


        virtual void EnsureConsistentBinning(const _hist * rhs, const char * caller, double tol = 1e-5) const final;
    };
}