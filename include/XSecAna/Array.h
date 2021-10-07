#pragma once

#include <Eigen/Dense>

namespace xsec {
    typedef Eigen::ArrayXd Array;

    class WeightedArray {
    public:
        WeightedArray() = default;

        WeightedArray(Array array,
                      double weight)
                : fArray(array), fWeight(weight) {}

        WeightedArray(const WeightedArray & rhs)
                : fArray(rhs.fArray),
                  fWeight(rhs.fWeight)
        {}

        WeightedArray(const WeightedArray && rhs)
                : fArray(std::move(rhs.fArray)),
                  fWeight(std::move(rhs.fWeight))
        {}

        Array array() const { return fArray; }
        Array & array() { return fArray; }
        double weight() const { return fWeight; }
        unsigned int size() const { return fArray.size(); }



        WeightedArray ScaleByWeight(double new_weight) const;

        virtual WeightedArray operator-(const WeightedArray & rhs) const;

        virtual WeightedArray operator+(const WeightedArray & rhs) const;

        virtual WeightedArray operator/(const WeightedArray & rhs) const;

        virtual WeightedArray operator*(const WeightedArray & rhs) const;

        virtual bool operator==(const WeightedArray & rhs) const;

        virtual bool operator!=(const WeightedArray & rhs) const { return !(*this == rhs); }

        virtual WeightedArray operator-=(const WeightedArray & rhs);

        virtual WeightedArray operator+=(const WeightedArray & rhs);

        virtual WeightedArray operator/=(const WeightedArray & rhs);

        virtual WeightedArray operator*=(const WeightedArray & rhs);

        virtual WeightedArray operator-(const double & rhs) const;

        virtual WeightedArray operator+(const double & rhs) const;

        virtual WeightedArray operator/(const double & rhs) const;

        virtual WeightedArray operator*(const double & rhs) const;

        virtual WeightedArray operator-=(const double & rhs);

        virtual WeightedArray operator+=(const double & rhs);

        virtual WeightedArray operator/=(const double & rhs);

        virtual WeightedArray operator*=(const double & rhs);

        virtual WeightedArray TrueDivide(const WeightedArray & rhs) const;

        // some convenience functions
        virtual WeightedArray abs() const;

        virtual WeightedArray abs2() const;

        virtual WeightedArray sqrt() const;

        virtual WeightedArray pow(double exp) const;
        
        // assignment operators
        virtual WeightedArray operator=(const WeightedArray & rhs);
        
        virtual WeightedArray & operator=(WeightedArray && rhs);

    private:
        Array fArray;
        double fWeight = 1;
    };

    WeightedArray
    WeightedArray::
    ScaleByWeight(double new_weight) const {
        return WeightedArray(fArray * (new_weight / fWeight),
                             new_weight);
    }

    WeightedArray
    WeightedArray::
    operator-(const WeightedArray & rhs) const {
        WeightedArray ret = *this; // copy this
        ret -= rhs;
        return ret;
    }

    WeightedArray
    WeightedArray::
    operator+(const WeightedArray & rhs) const {
        WeightedArray ret = *this; // copy this
        ret += rhs;
        return ret;
    }

    WeightedArray
    WeightedArray::
    operator/(const WeightedArray & rhs) const {
        WeightedArray ret = *this; // copy this
        ret /= rhs;
        return ret;
    }

    WeightedArray
    WeightedArray::
    operator*(const WeightedArray & rhs) const {
        WeightedArray ret = *this; // copy this
        ret *= rhs;
        return ret;
    }

    bool
    WeightedArray::
    operator==(const WeightedArray & rhs) const {
        return (fArray - rhs.fArray * fWeight / rhs.fWeight).isZero(0);
    }

    WeightedArray
    WeightedArray::
    operator-=(const WeightedArray & rhs) {
        fArray -= rhs.fArray * fWeight / rhs.fWeight;
        return *this;
    }

    WeightedArray
    WeightedArray::
    operator+=(const WeightedArray & rhs) {
        fArray += rhs.fArray * fWeight / rhs.fWeight;
        return *this;
    }

    WeightedArray
    WeightedArray::
    operator/=(const WeightedArray & rhs) {
        fArray /= rhs.fArray * fWeight / rhs.fWeight;
        fWeight = 1;
        return *this;
    }

    WeightedArray
    WeightedArray::
    operator*=(const WeightedArray & rhs) {
        fArray *= rhs.fArray * fWeight / rhs.fWeight;
        return *this;
    };

    WeightedArray
    WeightedArray::
    operator-(const double & rhs) const {
        WeightedArray ret = *this;
        ret -= rhs;
        return ret;
    }

    WeightedArray
    WeightedArray::
    operator+(const double & rhs) const {
        WeightedArray ret = *this;
        ret += rhs;
        return ret;
    }

    WeightedArray
    WeightedArray::
    operator/(const double & rhs) const {
        WeightedArray ret = *this;
        ret /= rhs;
        return ret;
    }

    WeightedArray
    WeightedArray::
    operator*(const double & rhs) const {
        WeightedArray ret = *this;
        ret *= rhs;
        return ret;
    }

    WeightedArray
    WeightedArray::
    operator-=(const double & rhs) {
        fArray -= rhs;
        return *this;
    }

    WeightedArray
    WeightedArray::
    operator+=(const double & rhs) {
        fArray += rhs;
        return *this;
    }

    WeightedArray
    WeightedArray::
    operator/=(const double & rhs) {
        fArray /= rhs;
        return *this;
    }

    WeightedArray
    WeightedArray::
    operator*=(const double & rhs) {
        fArray *= rhs;
        return *this;
    }

    WeightedArray
    WeightedArray::
    TrueDivide(const WeightedArray & rhs) const {
        WeightedArray ret = *this;
        ret.fArray /= rhs.fArray;
        return ret;
    }

    // some convenience functions
    WeightedArray
    WeightedArray::
    abs() const {
        return WeightedArray(fArray.abs(),
                             fWeight);
    }

    WeightedArray
    WeightedArray::
    abs2() const {
        return WeightedArray(fArray.abs2(),
                             std::pow(fWeight, 2));
    }

    WeightedArray
    WeightedArray::
    sqrt() const {
        return WeightedArray(fArray.sqrt(),
                             std::sqrt(fWeight));
    }
    WeightedArray
    WeightedArray::
    pow(double exp) const {
        return WeightedArray(fArray.pow(exp),
                             std::pow(fWeight, exp));
    }
    
    WeightedArray
    WeightedArray::
    operator=(const WeightedArray & rhs) {
        if (this == &rhs) return *this;
        fArray = rhs.fArray;
        fWeight = rhs.fWeight;
        return *this;
    }


    WeightedArray &
    WeightedArray::
    operator=(WeightedArray && rhs) {
        if(this == &rhs) return *this;
        fArray = std::move(rhs.fArray);
        fWeight = rhs.fWeight;
        return *this;
    }
}
