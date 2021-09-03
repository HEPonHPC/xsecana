#pragma once

#include "XSecAna/Systematic.h"
#include "XSecAna/Hist.h"
#include "XSecAna/IUnfold.h"

namespace xsec {
    /// \brief an UncertaintyPropagator is module
    /// that calculates uncertainty for an analysis
    ///
    /// IUncertaintyPropagator defines the interface assumed
    /// by the analysis object
    template<class HistType,
             class T,
             class ... Args>
    class IUncertaintyPropagator {
    public:
        virtual std::pair<HistType, HistType>
        TotalFractionalUncertainty(T & nominal_obj,
                                   std::map<std::string,
                                            xsec::Systematic<T> > & shifted_objs,
                                   Args & ... args) = 0;

        virtual std::pair<HistType, HistType>
        TotalAbsoluteUncertainty(T & nominal_obj,
                                 std::map<std::string,
                                          xsec::Systematic<T> > & shifted_objs,
                                 Args & ... args) = 0;

        virtual HistType
        FractionalUncertainty(T & nominal_obj,
                              xsec::Systematic<T> & shifted_obj,
                              Args & ... args) = 0;
        virtual HistType
        AbsoluteUncertainty(T & nominal_obj,
                            xsec::Systematic<T> & shifted_obj,
                            Args & ... args) = 0;
    };

    // inline some common functions
    template<typename Scalar,
            int Cols>
    inline Hist <Scalar, Cols>
    MaxShift(const Hist <Scalar, Cols> & h1,
             const Hist <Scalar, Cols> & h2) {
        // stack then take max value in each column
        Eigen::Matrix<Scalar, 2, Cols> stack;
        stack << h1.Contents(), h2.Contents();

        return Hist<Scalar, Cols>(stack.colwise().maxCoeff(),
                                  h1.Edges(),
                                  h1.Exposure());
    }

    template<typename Scalar,
            int Cols>
    inline Hist <Scalar, Cols>
    QuadSum(std::vector<Hist < Scalar, Cols>

    > deltas) {
    auto exposure = deltas[0].Exposure();
    // put deltas into an Eigen::Matrix for efficiency column-wise operations
    Eigen::Matrix<Scalar, Eigen::Dynamic, Cols> mat(deltas.size(), Cols);
    for(
    auto irow = 0u;
    irow<deltas.

    size();

    irow++) {
    mat.
    row(irow) = deltas[irow].Contents() * exposure / deltas[irow].Exposure();
}
return
Hist<Scalar, Cols>(mat
.

colwise()

.

squaredNorm(),
        deltas[0]

.

Edges(),
        std::pow(deltas[0].Exposure(), 2)

);
}
}
