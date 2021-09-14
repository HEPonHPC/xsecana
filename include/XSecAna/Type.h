#pragma once

#include <type_traits>
#include "XSecAna/Hist.h"
#include "TDirectory.h"

namespace xsec {
    namespace type {
        template<class Scalar, int Cols>
        std::true_type _is_hist(const Hist<Scalar, Cols> *);

        std::false_type _is_hist(...);

        template<class T>
        struct IsHist : decltype(_is_hist(std::declval<T *>())) {
        };

        template<class I>
        using LoadFunction = std::unique_ptr<I> (*) (TDirectory*, const std::string&);
    }
}
