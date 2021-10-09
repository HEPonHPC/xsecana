#pragma once

#include <Eigen/Dense>
#include "XSecAna/Array.h"

namespace xsec {
    typedef Eigen::ArrayXXd Array2D;

    class WeightedArray2D {
    public:
        WeightedArray2D(Array2D array,
                        double weight)
                : fArray(array.reshaped(array.size(), 1),
                         weight),
                  fShape{array.rows(), array.cols()}
        {}

        WeightedArray2D(const WeightedArray & array,
                        const shape_t & shape)
                : fArray(array),
                  fShape(shape)
        {}

        WeightedArray2D(const WeightedArray2D & rhs)
                : fArray(rhs.fArray), fShape(rhs.fShape)
        {}

        WeightedArray2D(const WeightedArray2D && rhs)
                : fArray(std::move(rhs.fArray)), fShape(rhs.fShape)
        {}

        double operator()(int i, int j) { return this->fArray(gid(i, j)); }

        shape_t shape() const { return fShape; }

        Array2D array() const {
            return Eigen::Map<const Array2D>(fArray.array().data(),
                                             fShape[0],
                                             fShape[1]);
        }

        WeightedArray2D ScaleByWeight(double new_weight) const;

        virtual WeightedArray2D operator-(const WeightedArray2D & rhs) const;

        virtual WeightedArray2D operator+(const WeightedArray2D & rhs) const;

        virtual WeightedArray2D operator/(const WeightedArray2D & rhs) const;

        virtual WeightedArray2D operator*(const WeightedArray2D & rhs) const;

        virtual bool operator==(const WeightedArray2D & rhs) const;

        virtual bool operator!=(const WeightedArray2D & rhs) const { return !(*this == rhs); }

        virtual WeightedArray2D operator-=(const WeightedArray2D & rhs);

        virtual WeightedArray2D operator+=(const WeightedArray2D & rhs);

        virtual WeightedArray2D operator/=(const WeightedArray2D & rhs);

        virtual WeightedArray2D operator*=(const WeightedArray2D & rhs);

        virtual WeightedArray2D operator-(const double & rhs) const;

        virtual WeightedArray2D operator+(const double & rhs) const;

        virtual WeightedArray2D operator/(const double & rhs) const;

        virtual WeightedArray2D operator*(const double & rhs) const;

        virtual WeightedArray2D operator-=(const double & rhs);

        virtual WeightedArray2D operator+=(const double & rhs);

        virtual WeightedArray2D operator/=(const double & rhs);

        virtual WeightedArray2D operator*=(const double & rhs);

        virtual WeightedArray2D TrueDivide(const WeightedArray2D & rhs) const;

        // some convenience functions
        virtual WeightedArray2D abs() const;

        virtual WeightedArray2D abs2() const;

        virtual WeightedArray2D sqrt() const;

        virtual WeightedArray2D pow(double exp) const;

        // assignment operators
        virtual WeightedArray2D operator=(const WeightedArray2D & rhs);

        virtual WeightedArray2D & operator=(WeightedArray2D && rhs);

    private:
        int gid(const int & i, const int & j) const {
            return j * fShape[0] + i;
        }

        WeightedArray fArray;
        shape_t fShape;
    };

    WeightedArray2D
    WeightedArray2D::
    ScaleByWeight(double new_weight) const {
        return {fArray.ScaleByWeight(new_weight),
                fShape};
    }


    WeightedArray2D
    WeightedArray2D::
    operator-(const WeightedArray2D & rhs) const {
        WeightedArray2D ret = *this;
        ret -= rhs;
        return ret;
    }

    WeightedArray2D
    WeightedArray2D::
    operator+(const WeightedArray2D & rhs) const {
        WeightedArray2D ret = *this;
        ret += rhs;
        return ret;
    }

    WeightedArray2D
    WeightedArray2D::
    operator/(const WeightedArray2D & rhs) const {
        WeightedArray2D ret = *this;
        ret /= rhs;
        return ret;
    }

    WeightedArray2D
    WeightedArray2D::
    operator*(const WeightedArray2D & rhs) const {
        WeightedArray2D ret = *this;
        ret *= rhs;
        return ret;
    }

    bool
    WeightedArray2D::
    operator==(const WeightedArray2D & rhs) const {
        if(fShape.size() != rhs.fShape.size()) return false;
        bool same = true;
        for(auto i = 0u; i < fShape.size(); i++) {
            same &= fShape[i] == rhs.fShape[i];
        }
        return same && fArray == rhs.fArray;
    }

    WeightedArray2D
    WeightedArray2D::
    operator-=(const WeightedArray2D & rhs) {
        fArray -= rhs.fArray;
        return *this;
    }

    WeightedArray2D
    WeightedArray2D::
    operator+=(const WeightedArray2D & rhs) {
        fArray += rhs.fArray;
        return *this;
    }

    WeightedArray2D
    WeightedArray2D::
    operator/=(const WeightedArray2D & rhs) {
        fArray /= rhs.fArray;
        return *this;
    }

    WeightedArray2D
    WeightedArray2D::
    operator*=(const WeightedArray2D & rhs) {
        fArray *= rhs.fArray;
        return *this;
    }

    WeightedArray2D
    WeightedArray2D::
    operator-(const double & rhs) const {
        WeightedArray2D ret = *this;
        ret -= rhs;
        return ret;
    }

    WeightedArray2D
    WeightedArray2D::
    operator+(const double & rhs) const {
        WeightedArray2D ret = *this;
        ret += rhs;
        return ret;
    }

    WeightedArray2D
    WeightedArray2D::
    operator/(const double & rhs) const {
        WeightedArray2D ret = *this;
        ret /= rhs;
        return ret;
    }

    WeightedArray2D
    WeightedArray2D::
    operator*(const double & rhs) const {
        WeightedArray2D ret = *this;
        ret *= rhs;
        return ret;
    }

    WeightedArray2D
    WeightedArray2D::
    operator-=(const double & rhs) {
        fArray -= rhs;
        return *this;
    }

    WeightedArray2D
    WeightedArray2D::
    operator+=(const double & rhs) {
        fArray += rhs;
        return *this;
    };

    WeightedArray2D
    WeightedArray2D::
    operator/=(const double & rhs) {
        fArray /= rhs;
        return *this;
    }

    WeightedArray2D
    WeightedArray2D::
    operator*=(const double & rhs) {
        fArray *= rhs;
        return *this;
    }

    WeightedArray2D
    WeightedArray2D::
    TrueDivide(const WeightedArray2D & rhs) const {
        return {fArray.TrueDivide(rhs.fArray), fShape};
    }

    // some convenience functions
    WeightedArray2D
    WeightedArray2D::
    abs() const {
        return {fArray.abs(), fShape};
    }

    WeightedArray2D
    WeightedArray2D::
    abs2() const {
        return {fArray.abs(), fShape};
    }

    WeightedArray2D
    WeightedArray2D::
    sqrt() const {
        return {fArray.sqrt(), fShape};
    }

    WeightedArray2D
    WeightedArray2D::
    pow(double exp) const {
        return {fArray.pow(exp), fShape};
    }

    WeightedArray2D
    WeightedArray2D::
    operator=(const WeightedArray2D & rhs) {
        if (this == &rhs) return *this;
        fArray = rhs.fArray;
        fShape = rhs.fShape;
        return *this;
    }

    WeightedArray2D &
    WeightedArray2D::
    operator=(WeightedArray2D && rhs) {
        if(this == &rhs) return *this;
        fArray = std::move(rhs.fArray);
        fShape = rhs.fShape;
        return *this;
    }

}