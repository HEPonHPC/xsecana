#pragma once

#include <Eigen/Dense>
#include <vector>
namespace xsec {
    typedef Eigen::ArrayXd Array;

    using shape_t = std::vector<long>;
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
        void set_weight(const double & new_weight) { fWeight = new_weight; }

        unsigned int size() const { return fArray.size(); }

        double operator()(const int & i) { return fArray(i); }

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

    protected:
        Array fArray;
        double fWeight = 1;
    };
}
