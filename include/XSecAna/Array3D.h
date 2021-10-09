#pragma once

#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>
#include "XSecAna/Array.h"

namespace xsec {
    typedef Eigen::Tensor<double, 3> Array3D;

    inline Array ToArray(const Array3D & array3D) {
        return Eigen::Map<const Array>(array3D.data(), array3D.size());
    }

    inline Array3D FromArray(const Array & array, const shape_t & shape) {
        return Eigen::TensorMap<const Array3D>(array.data(),
                                               shape[0],
                                               shape[1],
                                               shape[2]);
    }

    inline Array3D Array3D_Ones(int dimx, int dimy, int dimz) {
        return FromArray(Array::Ones(dimx*dimy*dimz), {dimx, dimy, dimz});
    }

    inline Array3D Array3D_Zero(int dimx, int dimy, int dimz) {
        return FromArray(Array::Zero(dimx*dimy*dimz), {dimx, dimy, dimz});
    }

    class WeightedArray3D {
    public:
        WeightedArray3D(Array3D array,
                        double weight)
                : fArray(ToArray(array),
                         weight),
                  fShape{array.dimension(0),
                         array.dimension(1),
                         array.dimension(2)}
        {}

        WeightedArray3D(const WeightedArray & array,
                        const shape_t & shape)
                : fArray(array),
                  fShape(shape)
        {}

        WeightedArray3D(const WeightedArray3D & rhs)
                : fArray(rhs.fArray), fShape(rhs.fShape)
        {}

        WeightedArray3D(const WeightedArray3D && rhs)
                : fArray(std::move(rhs.fArray)), fShape(rhs.fShape)
        {}

        double operator()(int i, int j, int k) { return this->fArray(gid(i, j, k)); }

        shape_t shape() const { return fShape; }

        Array3D array() const {
            return Eigen::TensorMap<const Array3D>(fArray.array().data(),
                                                   fShape[0],
                                                   fShape[1],
                                                   fShape[2]);
        }

        WeightedArray3D ScaleByWeight(double new_weight) const;

        virtual WeightedArray3D operator-(const WeightedArray3D & rhs) const;

        virtual WeightedArray3D operator+(const WeightedArray3D & rhs) const;

        virtual WeightedArray3D operator/(const WeightedArray3D & rhs) const;

        virtual WeightedArray3D operator*(const WeightedArray3D & rhs) const;

        virtual bool operator==(const WeightedArray3D & rhs) const;

        virtual bool operator!=(const WeightedArray3D & rhs) const { return !(*this == rhs); }

        virtual WeightedArray3D operator-=(const WeightedArray3D & rhs);

        virtual WeightedArray3D operator+=(const WeightedArray3D & rhs);

        virtual WeightedArray3D operator/=(const WeightedArray3D & rhs);

        virtual WeightedArray3D operator*=(const WeightedArray3D & rhs);

        virtual WeightedArray3D operator-(const double & rhs) const;

        virtual WeightedArray3D operator+(const double & rhs) const;

        virtual WeightedArray3D operator/(const double & rhs) const;

        virtual WeightedArray3D operator*(const double & rhs) const;

        virtual WeightedArray3D operator-=(const double & rhs);

        virtual WeightedArray3D operator+=(const double & rhs);

        virtual WeightedArray3D operator/=(const double & rhs);

        virtual WeightedArray3D operator*=(const double & rhs);

        virtual WeightedArray3D TrueDivide(const WeightedArray3D & rhs) const;

        // some convenience functions
        virtual WeightedArray3D abs() const;

        virtual WeightedArray3D abs2() const;

        virtual WeightedArray3D sqrt() const;

        virtual WeightedArray3D pow(double exp) const;

        // assignment operators
        virtual WeightedArray3D operator=(const WeightedArray3D & rhs);

        virtual WeightedArray3D & operator=(WeightedArray3D && rhs);

    private:
        int gid(const int & i, const int & j, const int & k) const {
            return (k * fShape[0] * fShape[1]) + (j * fShape[0]) + i;
        }

        WeightedArray fArray;
        shape_t fShape;
    };

    WeightedArray3D
    WeightedArray3D::
    ScaleByWeight(double new_weight) const {
        return {fArray.ScaleByWeight(new_weight),
                fShape};
    }


    WeightedArray3D
    WeightedArray3D::
    operator-(const WeightedArray3D & rhs) const {
        WeightedArray3D ret = *this;
        ret -= rhs;
        return ret;
    }

    WeightedArray3D
    WeightedArray3D::
    operator+(const WeightedArray3D & rhs) const {
        WeightedArray3D ret = *this;
        ret += rhs;
        return ret;
    }

    WeightedArray3D
    WeightedArray3D::
    operator/(const WeightedArray3D & rhs) const {
        WeightedArray3D ret = *this;
        ret /= rhs;
        return ret;
    }

    WeightedArray3D
    WeightedArray3D::
    operator*(const WeightedArray3D & rhs) const {
        WeightedArray3D ret = *this;
        ret *= rhs;
        return ret;
    }

    bool
    WeightedArray3D::
    operator==(const WeightedArray3D & rhs) const {
        if(fShape.size() != rhs.fShape.size()) return false;
        bool same = true;
        for(auto i = 0u; i < fShape.size(); i++) {
            same &= fShape[i] == rhs.fShape[i];
        }
        return same && fArray == rhs.fArray;
    }

    WeightedArray3D
    WeightedArray3D::
    operator-=(const WeightedArray3D & rhs) {
        fArray -= rhs.fArray;
        return *this;
    }

    WeightedArray3D
    WeightedArray3D::
    operator+=(const WeightedArray3D & rhs) {
        fArray += rhs.fArray;
        return *this;
    }

    WeightedArray3D
    WeightedArray3D::
    operator/=(const WeightedArray3D & rhs) {
        fArray /= rhs.fArray;
        return *this;
    }

    WeightedArray3D
    WeightedArray3D::
    operator*=(const WeightedArray3D & rhs) {
        fArray *= rhs.fArray;
        return *this;
    }

    WeightedArray3D
    WeightedArray3D::
    operator-(const double & rhs) const {
        WeightedArray3D ret = *this;
        ret -= rhs;
        return ret;
    }

    WeightedArray3D
    WeightedArray3D::
    operator+(const double & rhs) const {
        WeightedArray3D ret = *this;
        ret += rhs;
        return ret;
    }

    WeightedArray3D
    WeightedArray3D::
    operator/(const double & rhs) const {
        WeightedArray3D ret = *this;
        ret /= rhs;
        return ret;
    }

    WeightedArray3D
    WeightedArray3D::
    operator*(const double & rhs) const {
        WeightedArray3D ret = *this;
        ret *= rhs;
        return ret;
    }

    WeightedArray3D
    WeightedArray3D::
    operator-=(const double & rhs) {
        fArray -= rhs;
        return *this;
    }

    WeightedArray3D
    WeightedArray3D::
    operator+=(const double & rhs) {
        fArray += rhs;
        return *this;
    };

    WeightedArray3D
    WeightedArray3D::
    operator/=(const double & rhs) {
        fArray /= rhs;
        return *this;
    }

    WeightedArray3D
    WeightedArray3D::
    operator*=(const double & rhs) {
        fArray *= rhs;
        return *this;
    }

    WeightedArray3D
    WeightedArray3D::
    TrueDivide(const WeightedArray3D & rhs) const {
        return {fArray.TrueDivide(rhs.fArray), fShape};
    }

    // some convenience functions
    WeightedArray3D
    WeightedArray3D::
    abs() const {
        return {fArray.abs(), fShape};
    }

    WeightedArray3D
    WeightedArray3D::
    abs2() const {
        return {fArray.abs(), fShape};
    }

    WeightedArray3D
    WeightedArray3D::
    sqrt() const {
        return {fArray.sqrt(), fShape};
    }

    WeightedArray3D
    WeightedArray3D::
    pow(double exp) const {
        return {fArray.pow(exp), fShape};
    }

    WeightedArray3D
    WeightedArray3D::
    operator=(const WeightedArray3D & rhs) {
        if (this == &rhs) return *this;
        fArray = rhs.fArray;
        fShape = rhs.fShape;
        return *this;
    }

    WeightedArray3D &
    WeightedArray3D::
    operator=(WeightedArray3D && rhs) {
        if(this == &rhs) return *this;
        fArray = std::move(rhs.fArray);
        fShape = rhs.fShape;
        return *this;
    }

}