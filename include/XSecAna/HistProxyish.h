#pragma once

#include "XSecAna/Hist.h"
#include <exception>
#include <Eigen/Dense>

namespace xsec {
    namespace exceptions {
        class HistProxyishError : public std::exception {
        public:
            explicit HistProxyishError(const char * caller) {
                std::sprintf(fMsg,
                             "Attempted to call %s on a proxy",
                             caller);
            }

            const char * what() const throw() override {
                return &fMsg[0];
            }

        private:
            char fMsg[500];
        };
    }

#define NOTIMPLEMENTED(HistFunc) HistFunc override { throw xsec::exceptions::HistProxyishError(__FUNCTION__); }


    template<class WrappedType,
             class Scalar,
            int Cols>
    class HistProxyish : xsec::Hist<Scalar, Cols> {
    public:
        typedef xsec::Hist<Scalar, Cols> HistType;
        template<typename ... Args>
        HistProxyish(Args && ... args)
                : fWrapped(WrappedType(std::forward<Args>(args)...)) {}

        /// \brief Get the resource when it is ready
        /// This needs to be implemented for each type we are wrapping
        Hist<Scalar, Cols> GetHist();

        /// \brief First retrieve the wrapped data and save to file
        void SaveTo(TDirectory * dir, const std::string & subdir) const override;


        NOTIMPLEMENTED(void Normalize(const std::string & how));

        NOTIMPLEMENTED(Scalar Integrate() const);

        NOTIMPLEMENTED(HistType operator-(const HistType & rhs) const);

        NOTIMPLEMENTED(HistType operator+(const HistType & rhs) const);

        NOTIMPLEMENTED(HistType operator/(const HistType & rhs) const);

        NOTIMPLEMENTED(HistType operator*(const HistType & rhs) const);

        NOTIMPLEMENTED(bool operator==(const HistType & rhs) const);

        NOTIMPLEMENTED(bool operator!=(const HistType & rhs) const);

        NOTIMPLEMENTED(HistType operator-=(const HistType & rhs));

        NOTIMPLEMENTED(HistType operator+=(const HistType & rhs));

        NOTIMPLEMENTED(HistType operator/=(const HistType & rhs));

        NOTIMPLEMENTED(HistType operator*=(const HistType & rhs));

        NOTIMPLEMENTED(HistType operator-(const Scalar & rhs) const);

        NOTIMPLEMENTED(HistType operator+(const Scalar & rhs) const);

        NOTIMPLEMENTED(HistType operator/(const Scalar & rhs) const);

        NOTIMPLEMENTED(HistType operator*(const Scalar & rhs) const);

        NOTIMPLEMENTED(HistType operator-=(const Scalar & rhs));

        NOTIMPLEMENTED(HistType operator+=(const Scalar & rhs));

        NOTIMPLEMENTED(HistType operator/=(const Scalar & rhs));

        NOTIMPLEMENTED(HistType operator*=(const Scalar & rhs));

        NOTIMPLEMENTED(HistType ScaleByExposure(Scalar new_expo) const);

        NOTIMPLEMENTED(HistType TrueDivide(const HistType & rhs) const);

        // some convenience functions
        NOTIMPLEMENTED(HistType abs() const);

        NOTIMPLEMENTED(HistType abs2() const);

        NOTIMPLEMENTED(HistType sqrt() const);

        NOTIMPLEMENTED(HistType pow(Scalar exp) const);

        NOTIMPLEMENTED(Scalar Exposure() const);

        NOTIMPLEMENTED(unsigned int size() const);

        NOTIMPLEMENTED(Scalar & operator[](int index));

        NOTIMPLEMENTED(Scalar operator[](int index) const);

        // having a hard time wrapping these with the macro
        const Eigen::Array<Scalar, 1, Cols> & Contents() const override {
            throw xsec::exceptions::HistProxyishError(__FUNCTION__);
        }
        const Eigen::Array<Scalar, 1, EdgesSize(Cols)> & Edges() const override{
            throw xsec::exceptions::HistProxyishError(__FUNCTION__);
        }
        Eigen::Array<Scalar, 1, Cols> BinWidths() const override {
            throw xsec::exceptions::HistProxyishError(__FUNCTION__);
        }


    private:
        mutable WrappedType fWrapped;
    };

    template<class WrappedType,
             class Scalar,
            int Cols>
    void
    HistProxyish<WrappedType,
                 Scalar,
                 Cols>::
    SaveTo(TDirectory * dir, const std::string & subdir) const {
        fWrapped.GetHist().SaveTo(dir, subdir);
    }


}

