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
    class SimpleQuadSum : public IUncertaintyPropagator<HistType, T, Args...> {
    public:

        SimpleQuadSum() = default;

        std::pair<HistType, HistType>
        TotalFractionalUncertainty(T * nominal_obj,
                                   std::map<std::string,
                                            xsec::Systematic<T> > & shifted_objs,
                                   Args & ... args) override;

        std::pair<HistType, HistType>
        TotalAbsoluteUncertainty(T * nominal_obj,
                                 std::map<std::string,
                                          xsec::Systematic<T> > & shifted_objs,
                                 Args & ...  args) override;

        HistType
        FractionalUncertainty(T * nominal_obj,
                              xsec::Systematic<T> & shifted_obj,
                              Args & ...  args) override;

        HistType
        AbsoluteUncertainty(T * nominal_obj,
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
                                        new HistType(MultiverseShift(syst, nominal, 1)),
                                        new HistType(MultiverseShift(syst, nominal, 1)));
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
    AbsoluteUncertainty(T * nominal_obj,
                        xsec::Systematic<T> & shifted_obj,
                        Args & ...  args) {

        // calculate cross sections
        HistType hnominal = nominal_obj->Eval(std::forward<Args>(args)...);
        std::vector<HistType*> shifts(shifted_obj.GetShifts().size());
        for(auto i = 0u; i < shifted_obj.GetShifts().size(); i++) {
            shifts[i] = new HistType(shifted_obj.GetShifts()[i]->Eval(std::forward<Args>(args)...));
        }

        Systematic<HistType> systematic(shifted_obj.GetName(),
                                        shifts,
                                        shifted_obj.GetType());

        // convert multiverse systematic to two-sided by finding 1sigma
        systematic = HandleMultiverseSystematic(systematic, hnominal);

        ForEachFunction<HistType, HistType> subtract = [&hnominal](HistType * h) {
            return new HistType(*h - hnominal);
        };

        systematic = systematic.ForEach(subtract);

        // if systematic is one-sided, symmeterize it
        if(systematic.GetType() == kOneSided) {
            return systematic.Up()->abs();
        }
        else{
            return MaxShift(systematic.Up()->abs(), systematic.Down()->abs());
        }
    }

    /////////////////////////////////////////////////////////////////////////
    template<class HistType,
             class T,
             class ... Args>
    HistType
    SimpleQuadSum<HistType, T, Args...>::
    FractionalUncertainty(T * nominal_obj,
                          xsec::Systematic<T> & shifted_obj,
                          Args & ...  args) {
        HistType abs = AbsoluteUncertainty(nominal_obj,
                                           shifted_obj,
                                           args...);
        return abs / nominal_obj->Eval(std::forward<Args>(args)...);
    }

    /////////////////////////////////////////////////////////////////////////
    template<class HistType,
             class T,
             class ... Args>
    std::pair<HistType, HistType>
    SimpleQuadSum<HistType, T, Args...>::
    TotalFractionalUncertainty(T * nominal_obj,
                               std::map<std::string,
                                        xsec::Systematic<T>> & shifted_objs,
                               Args & ... args) {
        auto hnominal = nominal_obj->Eval(std::forward<Args>(args)...);

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
    TotalAbsoluteUncertainty(T * nominal_obj,
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
