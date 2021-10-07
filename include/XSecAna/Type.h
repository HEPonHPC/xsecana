#pragma once

#include <type_traits>
#include "XSecAna/Hist.h"
#include "TDirectory.h"

namespace xsec {
    namespace type {
        template<class I>
        using LoadFunction = std::unique_ptr<I> (*) (TDirectory*, const std::string&);
    }
}
