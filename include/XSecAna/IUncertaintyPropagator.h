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
        TotalFractionalUncertainty(T * nominal_obj,
                                   std::map<std::string,
                                            xsec::Systematic<T> > & shifted_objs,
                                   Args & ... args) = 0;

        virtual std::pair<HistType, HistType>
        TotalAbsoluteUncertainty(T * nominal_obj,
                                 std::map<std::string,
                                          xsec::Systematic<T> > & shifted_objs,
                                 Args & ... args) = 0;

        virtual HistType
        FractionalUncertainty(T * nominal_obj,
                              xsec::Systematic<T> & shifted_obj,
                              Args & ... args) = 0;

        virtual HistType
        AbsoluteUncertainty(T * nominal_obj,
                            xsec::Systematic<T> & shifted_obj,
                            Args & ... args) = 0;
    };

    // inline some common functions
    inline Hist
    MaxShift(const Hist & h1,
             const Hist & h2) {
        // stack then take max value in each column
        Eigen::Matrix<double, 2, Eigen::Dynamic> stack(2, h1.ContentsAndUOF().size());
        stack.row(0) = h1.ContentsAndUOF();
        stack.row(1) = h2.ContentsAndUOF();

        return Hist(stack.colwise().maxCoeff(),
                    h1.EdgesAndUOF(),
                    h1.Exposure());
    }

    inline Hist
    QuadSum(std::vector<const xsec::Hist *> deltas) {
        auto exposure = deltas[0]->Exposure();
        // put deltas into an Eigen::Matrix for efficiency column-wise operations
        Eigen::MatrixXd mat(deltas.size(),
                            deltas[0]->ContentsAndUOF().size());
        for (auto irow = 0u; irow < deltas.size(); irow++) {
            mat.row(irow) = deltas[irow]->ContentsAndUOF() * exposure / deltas[irow]->Exposure();
        }
        return Hist(mat.colwise().squaredNorm(),
                    deltas[0]->EdgesAndUOF(),
                    std::pow(deltas[0]->Exposure(), 2));
    }
}
