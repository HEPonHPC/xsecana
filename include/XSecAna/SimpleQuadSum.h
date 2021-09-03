#pragma once

#include "XSecAna/Systematic.h"
#include "XSecAna/Hist.h"
#include "XSecAna/IUnfold.h"
#include "XSecAna/IUncertaintyPropagator.h"

#include <Eigen/Dense>

namespace xsec {

    /// \brief an UncertaintyPropagator is module
    /// that calculates uncertainty for an analysis
    ///
    /// SimpleQuadSum performs the quadrature sum of systematic shifts
    /// In the case of an asymmetric 2-sided shift,
    /// the shift is symmeterized by taking the largest shift
    template<class HistType,
             class T,
             class ... Args>
    class SimpleQuadSum : IUncertaintyPropagator<HistType, T, Args...> {
    public:

        SimpleQuadSum() = default;

        std::pair<HistType, HistType>
        TotalFractionalUncertainty(T & nominal_obj,
                                   std::map<std::string,
                                            xsec::Systematic<T> > & shifted_objs,
                                   Args & ... args) override;

        std::pair<HistType, HistType>
        TotalAbsoluteUncertainty(T & nominal_obj,
                                 std::map<std::string,
                                          xsec::Systematic<T> > & shifted_objs,
                                 Args & ...  args) override;

        HistType
        FractionalUncertainty(T & nominal_obj,
                              xsec::Systematic<T> & shifted_obj,
                              Args & ...  args) override;

        HistType
        AbsoluteUncertainty(T & nominal_obj,
                            xsec::Systematic<T> & shifted_obj,
                            Args & ...  args) override;

    };

    /////////////////////////////////////////////////////////////////////////
    template<class HistType>
    inline Systematic<HistType> HandleMultiverseSystematic(const Systematic<HistType> & syst,
                                                           const HistType & nominal) {
        // if given a multiverse systematic, return a new Systematic<HistType> that contains
        // +/- 1 sigma
        if (syst.GetType() == SystType_t::kMultiverse) {
            return Systematic<HistType>(syst.GetName(),
                                        syst.NSigmaShift(1, nominal),
                                        syst.NSigmaShift(-1, nominal));
        } else {
            return syst;
        }

    }

    /////////////////////////////////////////////////////////////////////////
    template<class HistType,
             class T,
             class ... Args>
    HistType
    SimpleQuadSum<HistType, T, Args...>::
    AbsoluteUncertainty(T & nominal_obj,
                        xsec::Systematic<T> & shifted_obj,
                        Args & ...  args) {
        // calculate cross sections
        auto hnominal = nominal_obj.Eval(std::forward<Args>(args)...);
        Systematic<HistType> shifts = shifted_obj.Invoke(&T::Eval,
                                                         std::forward<Args>(args)...);

        // convert multiverse systematic to two-sided by finding 1sigma
        shifts = HandleMultiverseSystematic(shifts, hnominal);

        // HistType::operator- is overloaded
        // static cast to resolve
        shifts = shifts.Invoke(static_cast<HistType(HistType::*)(const HistType &) const>(&HistType::operator-),
                               hnominal);

        return MaxShift(shifts.Up().abs(), shifts.Down().abs());
    }

    /////////////////////////////////////////////////////////////////////////
    template<class HistType,
             class T,
             class ... Args>
    HistType
    SimpleQuadSum<HistType, T, Args...>::
    FractionalUncertainty(T & nominal_obj,
                          xsec::Systematic<T> & shifted_obj,
                          Args & ...  args) {
        HistType abs = AbsoluteUncertainty(nominal_obj,
                                           shifted_obj,
                                           args...);
        return abs / nominal_obj.Eval(std::forward<Args>(args)...);
    }

    /////////////////////////////////////////////////////////////////////////
    template<class HistType,
             class T,
             class ... Args>
    std::pair<HistType, HistType>
    SimpleQuadSum<HistType, T, Args...>::
    TotalFractionalUncertainty(T & nominal_obj,
                               std::map<std::string,
                                        xsec::Systematic<T>> & shifted_objs,
                               Args & ... args) {
        auto hnominal = nominal_obj.Eval(std::forward<Args>(args)...);
        auto abs = TotalAbsoluteUncertainty(nominal_obj,
                                            shifted_objs,
                                            args...);
        return {abs.first / hnominal, abs.second / hnominal};
    }

    /////////////////////////////////////////////////////////////////////////
    template<class HistType,
             class T,
             class ... Args>
    std::pair<HistType, HistType>
    SimpleQuadSum<HistType, T, Args...>::
    TotalAbsoluteUncertainty(T & nominal_obj,
                               std::map<std::string,
                                        xsec::Systematic<T>> & shifted_objs,
                               Args & ... args) {
        std::vector<HistType> shifts;
        for (auto syst_it = shifted_objs.begin(); syst_it != shifted_objs.end(); syst_it++) {
            shifts.push_back(AbsoluteUncertainty(nominal_obj,
                                                 syst_it->second,
                                                 args...));
        }
        auto result = QuadSum(shifts).sqrt();
        return {result, result};
    }
}
