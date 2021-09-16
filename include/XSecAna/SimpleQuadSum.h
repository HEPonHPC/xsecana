#pragma once

#include "XSecAna/Systematic.h"
#include "XSecAna/Hist.h"
#include "XSecAna/IUnfold.h"
#include "XSecAna/IUncertaintyPropagator.h"

#include <Eigen/Dense>
#include "XSecAna/Type.h"

namespace xsec {
    /// \brief
    /// SimpleQuadSum performs the quadrature sum of systematic shifts
    /// In the case of an asymmetric 2-sided shift,
    /// the shift is symmeterized by taking the largest shift
    /// Type T must have implemented T::Eval(Args ... args)->HistType
    /// in order to transform into a HistType.
    /// The only exception is if type T is type HistType.
    namespace SimpleQuadSum {
        namespace {
            /// \brief Evaluate all the T's out to Hists
            /// Calling Eval on the T's the long way since
            /// passing parameter packs through a lambda capture
            /// is pretty opaque in < C++20
            template<class HistType,
                     class T,
                     class ... Args>
            inline
            Systematic<HistType>
            EvalSystematic(xsec::Systematic<T> & shifted_obj,
                           Args && ...  args) {
                std::vector<HistType*> shifts(shifted_obj.GetShifts().size());
                for (auto i = 0u; i < shifted_obj.GetShifts().size(); i++) {
                    shifts[i] = new HistType(shifted_obj.GetShifts()[i]->Eval(std::forward<Args>(args)...));
                }
                return Systematic<HistType>(shifted_obj.GetName(),
                                            shifts,
                                            shifted_obj.GetType());
            }

            /// \brief Internal function for calculating absolute uncertainty
            /// on Systematic<HistType>s
            template<class HistType>
            inline
            std::pair<HistType, Systematic<HistType>>
            _AbsoluteUncertainty(const HistType & nominal,
                                 const xsec::Systematic<HistType> & shifted_obj) {
                std::vector<HistType *> shifts(shifted_obj.GetShifts().size());
                HistType * up;
                HistType * down;
                // convert multiverse systematic to two-sided by finding 1sigma
                if (shifted_obj.GetType() == kMultiverse) {
                    up = new HistType(MultiverseShift(shifted_obj, nominal, 1));
                    down = new HistType(MultiverseShift(shifted_obj, nominal, -1));
                } else if (shifted_obj.GetType() == kTwoSided) {
                    up = new HistType(*shifted_obj.GetShifts()[0] - nominal);
                    down = new HistType(*shifted_obj.GetShifts()[1] - nominal);
                } else {
                    up = new HistType(*shifted_obj.GetShifts()[0] - nominal);
                    down = up;
                }

                auto diff = new HistType(MaxShift(up->abs(), down->abs()));
                return {nominal,
                        Systematic<HistType>(shifted_obj.GetName(),
                                             diff,
                                             diff)};
            }

            /// \brief Internal function for calculating fractional uncertainty
            /// on Systematic<HistType>s
            template<class HistType>
            std::pair<HistType, Systematic<HistType>>
            _FractionalUncertainty(const HistType & nominal,
                                   const xsec::Systematic<HistType> & shifted_obj) {
                auto frac = new HistType(*(std::get<1>(_AbsoluteUncertainty(nominal,
                                                                            shifted_obj)).Up()) / nominal);
                return {nominal, Systematic<HistType>(shifted_obj.GetName(), frac, frac)};

            }

            /// \brief Internal function for calculating total absolute uncertainty
            /// on Systematic<HistType>s
            template<class HistType>
            std::pair<HistType, Systematic<HistType>>
            _TotalAbsoluteUncertainty(const HistType & nominal,
                                      const std::map<std::string,
                                                     xsec::Systematic<HistType>> & shifted_objs) {
                std::vector<const HistType*> shifts;
                for (auto syst_it = shifted_objs.begin(); syst_it != shifted_objs.end(); syst_it++) {
                    shifts.push_back(std::get<1>(_AbsoluteUncertainty<HistType>(nominal,
                                                                                syst_it->second)).Up());
                }
                auto result = new HistType(QuadSum(shifts).sqrt());
                return {nominal, Systematic<HistType>("Total Absolute Uncertainty",
                                                      result,
                                                      result)};
            }

            /// \brief Internal function for calculating total fractional uncertainty
            /// on Systematic<HistType>s
            template<class HistType>
            std::pair<HistType, Systematic<HistType>>
            _TotalFractionalUncertainty(const HistType & nominal,
                                        const std::map<std::string,
                                                       xsec::Systematic<HistType>> & shifted_objs) {
                std::vector<const HistType*> shifts;
                for (auto syst_it = shifted_objs.begin(); syst_it != shifted_objs.end(); syst_it++) {
                    shifts.push_back(std::get<1>(_AbsoluteUncertainty<HistType>(nominal,
                                                                                syst_it->second)).Up());
                }
                auto result = new HistType(QuadSum(shifts).sqrt() / nominal);
                return {nominal, Systematic<HistType>("Total Fractional Uncertainty",
                                                      result,
                                                      result)};
            }
        }

        /// \brief Calculate absolute uncertainty for the given Systematic<T>.
        /// If T is not a histogram, T::Eval(args...)->HistTYpe is invoked
        /// to transform the Systematic<T> to a Systematic<HistType>,
        /// then the following calculations can be performed.
        /// If T is already a histogram, the following calculations are immediately performed.
        ///
        ///  - If type kMultiverse, the 1sigma universe is found.
        ///  - Otherwise, it is assumed the Systematic represents a 1sigma shift,
        ///    and as such, the absolute uncertainty is the difference between the Systematic
        ///    and nominal.
        ///    - If type kTwoSided, the returning Systematic is also two sided and symmeterized,
        ///      such that Up() == Down()
        ///    - If type kOneSided, the returning Systematic is also two sided and symmeterized,
        ///      such that Up() == Down()
        /// TODO can we add a switch for asymmetric errors?
        template<class HistType,
                 class T,
                 class ... Args>
        inline
        std::pair<HistType, Systematic<HistType>>
        AbsoluteUncertainty(T * nominal_obj,
                            xsec::Systematic<T> & shifted_obj,
                            Args & ...  args) {
            if constexpr(xsec::type::IsHist<T>()) {
                return _AbsoluteUncertainty(*nominal_obj, shifted_obj);
            }
            else {
                HistType hnominal = nominal_obj->Eval(std::forward<Args>(args)...);
                Systematic<HistType> hsystematic = EvalSystematic<HistType>(shifted_obj,
                                                                            std::forward<Args>(args)...);
                return _AbsoluteUncertainty(hnominal,
                                            hsystematic);
            }
        }

        /// \brief Calculate fractional uncertainty for the given Systematic<T>.
        /// If T is not a histogram, T::Eval(args...)->HistTYpe is invoked
        /// to transform the Systematic<T> to a Systematic<HistType>,
        /// then the following calculations can be performed.
        /// If T is already a histogram, the following calculations are immediately performed.
        ///
        ///  - If type kMultiverse, the 1sigma universe is found.
        ///  - Otherwise, it is assumed the Systematic represents a 1sigma shift,
        ///    and as such, the fractional is the difference between the Systematic
        ///    and nominal divided by nominal.
        ///    - If type kTwoSided, the returning Systematic is also two sided and symmeterized,
        ///      such that Up() == Down()
        ///    - If type kOneSided, the returning Systematic is also two sided and symmeterized,
        ///      such that Up() == Down()
        template<class HistType,
                 class T,
                 class ... Args>
        inline
        std::pair<HistType, Systematic<HistType>>
        FractionalUncertainty(T * nominal_obj,
                              xsec::Systematic<T> & shifted_obj,
                              Args & ...  args) {
            if constexpr(xsec::type::IsHist<T>()) {
                return _FractionalUncertainty(*nominal_obj, shifted_obj);
            }
            else {
                auto hnominal = nominal_obj->Eval(std::forward<Args>(args)...);
                auto hsystematic = EvalSystematic<HistType>(shifted_obj,
                                                            std::forward<Args>(args)...);
                return _FractionalUncertainty(hnominal, hsystematic);
            }
        }

        /// \brief Calculate the total absolute uncertainty for the given Systematic<T>'s.
        /// If T's are not histograms, T::Eval(args...)->HistTYpe is invoked
        /// to transform the Systematic<T>'s to Systematic<HistType>'s,
        /// then the following calculations can be performed.
        /// If T is already a histogram, the following calculations are immediately performed.
        ///
        /// For each Systematic<HistType>
        ///  - If type kMultiverse, the 1sigma universe is found.
        ///  - Otherwise, it is assumed the Systematic represents a 1sigma shift,
        ///    and as such, the absolute is the difference between the Systematic
        ///    and nominal.
        ///    - If type kTwoSided, the returning Systematic is also two sided and symmeterized,
        ///      such that Up() == Down()
        ///    - If type kOneSided, the returning Systematic is also two sided and symmeterized,
        ///      such that Up() == Down()
        /// then, the individual absolute uncertainties are added in quadrature.
        /// The returned Systematic contains the symmeterized total absolute uncertainty.
        template<class HistType,
                 class T,
                 class ... Args>
        inline
        std::pair<HistType, Systematic<HistType>>
        TotalAbsoluteUncertainty(T * nominal_obj,
                                 std::map<std::string,
                                          xsec::Systematic<T> > & shifted_objs,
                                 Args & ...  args) {
            if constexpr(xsec::type::IsHist<T>()) {
                return _TotalAbsoluteUncertainty(*nominal_obj, shifted_objs);
            }
            else {
                HistType hnominal = nominal_obj->Eval(std::forward<Args>(args)...);
                std::map<std::string, Systematic<HistType>> hsystematics;
                for (auto syst_it = shifted_objs.begin(); syst_it != shifted_objs.end(); syst_it++) {
                    hsystematics[syst_it->first] = EvalSystematic<HistType>(syst_it->second,
                                                                            std::forward<Args>(args)...);
                }

                return _TotalAbsoluteUncertainty(hnominal, hsystematics);
            }
        }

        /// \brief Calculate the total fractional uncertainty for the given Systematic<T>'s.
        /// If T's are not histograms, T::Eval(args...)->HistTYpe is invoked
        /// to transform the Systematic<T>'s to Systematic<HistType>'s,
        /// then the following calculations can be performed.
        /// If T is already a histogram, the following calculations are immediately performed.
        ///
        /// For each Systematic<HistType>
        ///  - If type kMultiverse, the 1sigma universe is found.
        ///  - Otherwise, it is assumed the Systematic represents a 1sigma shift,
        ///    and as such, the fractional is the difference between the Systematic
        ///    and nominal divided by nominal.
        ///    - If type kTwoSided, the returning Systematic is also two sided and symmeterized,
        ///      such that Up() == Down()
        ///    - If type kOneSided, the returning Systematic is also two sided and symmeterized,
        ///      such that Up() == Down()
        /// then, the individual fractional uncertainties are added in quadrature.
        /// The returned Systematic contains the symmeterized total fractional uncertainty.
        template<class HistType,
                 class T,
                 class ... Args>
        inline
        std::pair<HistType, Systematic<HistType>>
        TotalFractionalUncertainty(T * nominal_obj,
                                   std::map<std::string,
                                            xsec::Systematic<T> > & shifted_objs,
                                   Args & ... args) {
            if constexpr(xsec::type::IsHist<T>()) {
                return _TotalFractionalUncertainty(*nominal_obj, shifted_objs);
            }
            else {
                HistType hnominal = nominal_obj->Eval(std::forward<Args>(args)...);
                std::map<std::string, Systematic<HistType>> hsystematics;
                for (auto syst_it = shifted_objs.begin(); syst_it != shifted_objs.end(); syst_it++) {
                    hsystematics[syst_it->first] = EvalSystematic<HistType>(syst_it->second,
                                                                            std::forward<Args>(args)...);
                }

                return _TotalFractionalUncertainty(hnominal, hsystematics);
            }
        }
    }
}
