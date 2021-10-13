#pragma once

#include "XSecAna/Utils.h"

namespace xsec {
    inline void
    _QuadSum(ArrayRef a, const Array & b) {
        a = (a.pow(2) + b.pow(2)).sqrt();
    }

    inline Array
    QuadSum(const std::vector<Array> & deltas) {
        auto ret = deltas[0];
        for (auto i = 1u; i < deltas.size(); i++) {
            _QuadSum(ret, deltas[i]);
        }
        return ret;
    }


    inline TH1 *
    QuadSum(const std::vector<const TH1 *> & deltas) {
        std::vector<Array> _deltas_c(deltas.size());
        std::vector<Array> _deltas_e(deltas.size());
        for (auto i = 0u; i < deltas.size(); i++) {
            _deltas_c[i] = root::MapContentsToEigen(deltas[i]);
            _deltas_e[i] = root::MapErrorsToEigen(deltas[i]);

        }
        return root::ToROOT(QuadSum(_deltas_c),
                            QuadSum(_deltas_e),
                            root::TH1Props(deltas[0],
                                           root::MakeUnique(std::string(__FUNCTION__)).c_str()));
    }

    inline Array
    QuadSum(const Array & d1, const Array & d2) {
        auto ret = d1;
        _QuadSum(ret, d2);
        return ret;
    }

    inline Array
    QuadSum(const Array & d1, const Array & d2, const Array & d3) {
        auto ret = d1;
        _QuadSum(ret, d2);
        _QuadSum(ret, d3);
        return ret;
    }

    inline TH1 *
    QuadSum(const TH1 * d1, const TH1 * d2) {
        return QuadSum(std::vector<const TH1*>{d1, d2});
    }

    inline TH1 *
    QuadSum(const TH1 * d1, const TH1 * d2, const TH1 * d3) {
        return QuadSum(std::vector<const TH1*>{d1, d2, d3});
    }

    // inline some common functions
    inline Array
    MaxShift(const Array& arr1,
             const Array& arr2) {
        // stack then take max value in each column
        Eigen::Matrix<double, 2, Eigen::Dynamic> stack(2, arr1.size());
        stack.row(0) = arr1;
        stack.row(1) = arr2;
        return stack.colwise().maxCoeff();
    }
    inline TH1 *
    MaxShift(const TH1 * h1,
             const TH1 * h2) {
        return root::ToROOTLike(h1, MaxShift(root::MapContentsToEigen(h1),
                                             root::MapContentsToEigen(h2)));
    }


}