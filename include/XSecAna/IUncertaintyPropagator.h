#pragma once

#include "XSecAna/Systematic.h"
#include "XSecAna/IUnfold.h"
#include "TH1.h"
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
    inline _hist *
    MaxShift(const _hist * h1,
             const _hist * h2) {
        // stack then take max value in each column
        Eigen::Matrix<double, 2, Eigen::Dynamic> stack(2, h1->GetContentsAndUOF().size());
        stack.row(0) = h1->GetContentsAndUOF();
        stack.row(1) = h2->GetContentsAndUOF();

        auto ret = h1->Clone()
        ret->SetContentsAndUOF(stack.colwise().maxCoeff());
        return ret;
    }

    inline _hist *
    QuadSum(std::vector<const xsec::Hist *> deltas) {
        auto exposure = deltas[0]->Exposure();
        // put deltas into an Eigen::Matrix for efficiency column-wise operations
        Eigen::MatrixXd mat(deltas.size(),
                            deltas[0]->GetContentsAndUOF().size());
        for (auto irow = 0u; irow < deltas.size(); irow++) {
            mat.row(irow) = deltas[irow]->ScaleByExposure(exposure);
        }
        return Hist(mat.colwise().squaredNorm(),
                    deltas[0]->EdgesAndUOF(),
                    std::pow(deltas[0]->Exposure(), 2));
    }
}
