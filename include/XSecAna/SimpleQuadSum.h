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
    /// Type T must have implemented T::Eval(Args ... args)->Hist
    /// in order to transform into a Hist.
    /// The only exception is if type T is type Hist.
    namespace SimpleQuadSum {
        namespace {
            /// \brief Evaluate all the T's out to Hists
            /// Calling Eval on the T's the long way since
            /// passing parameter packs through a lambda capture
            /// is pretty opaque in < C++20
            template<class T,
                     class ... Args>
            inline
            Systematic<Hist>
            EvalSystematic(const xsec::Systematic<T> & shifted_obj,
                           Args && ...  args) {
                std::vector<Hist *> shifts(shifted_obj.GetShifts().size());
                for (auto i = 0u; i < shifted_obj.GetShifts().size(); i++) {
                    shifts[i] = new Hist(shifted_obj.GetShifts()[i]->Eval(std::forward<Args>(args)...));
                }
                return Systematic<Hist>(shifted_obj.GetName(),
                                        shifts,
                                        shifted_obj.GetType());
            }

            /// \brief Internal function for calculating absolute uncertainty
            /// on Systematic<Hist>s
            inline
            std::pair<Hist, Systematic<Hist>>
            _AbsoluteUncertainty(const Hist & nominal,
                                 const xsec::Systematic<Hist> & shifted_obj) {
                std::vector<Hist *> shifts(shifted_obj.GetShifts().size());
                Hist * up;
                Hist * down;
                // convert multiverse systematic to two-sided by finding 1sigma
                if (shifted_obj.GetType() == kMultiverse) {
                    up = new Hist(MultiverseShift(shifted_obj, nominal, 1));
                    down = new Hist(MultiverseShift(shifted_obj, nominal, -1));
                } else if (shifted_obj.GetType() == kTwoSided) {
                    up = new Hist(*shifted_obj.GetShifts()[0] - nominal);
                    down = new Hist(*shifted_obj.GetShifts()[1] - nominal);
                } else {
                    up = new Hist(*shifted_obj.GetShifts()[0] - nominal);
                    down = up;
                }

                auto diff = new Hist(MaxShift(up->abs(), down->abs()));
                return {nominal,
                        Systematic<Hist>(shifted_obj.GetName(),
                                         diff,
                                         diff)};
            }

            /// \brief Internal function for calculating fractional uncertainty
            /// on Systematic<Hist>s
            std::pair<Hist, Systematic<Hist>>
            _FractionalUncertainty(const Hist & nominal,
                                   const xsec::Systematic<Hist> & shifted_obj) {
                auto frac = new Hist(*(std::get<1>(_AbsoluteUncertainty(nominal,
                                                                        shifted_obj)).Up()) / nominal);
                return {nominal, Systematic<Hist>(shifted_obj.GetName(), frac, frac)};

            }

            /// \brief Internal function for calculating total absolute uncertainty
            /// on Systematic<Hist>s
            std::pair<Hist, Systematic<Hist>>
            _TotalAbsoluteUncertainty(const Hist & nominal,
                                      const std::map<std::string,
                                                     xsec::Systematic<Hist>> & shifted_objs) {
                std::vector<const Hist *> shifts;
                for (auto syst_it = shifted_objs.begin(); syst_it != shifted_objs.end(); syst_it++) {
                    shifts.push_back(std::get<1>(_AbsoluteUncertainty(nominal,
                                                                      syst_it->second)).Up());
                }
                auto result = new Hist(QuadSum(shifts).sqrt());
                return {nominal, Systematic<Hist>("Total Absolute Uncertainty",
                                                  result,
                                                  result)};
            }

            /// \brief Internal function for calculating total fractional uncertainty
            /// on Systematic<Hist>s
            std::pair<Hist, Systematic<Hist>>
            _TotalFractionalUncertainty(const Hist & nominal,
                                        const std::map<std::string,
                                                       xsec::Systematic<Hist>> & shifted_objs) {
                std::vector<const Hist *> shifts;
                for (auto syst_it = shifted_objs.begin(); syst_it != shifted_objs.end(); syst_it++) {
                    shifts.push_back(std::get<1>(_AbsoluteUncertainty(nominal,
                                                                      syst_it->second)).Up());
                }
                auto result = new Hist(QuadSum(shifts).sqrt() / nominal);
                return {nominal, Systematic<Hist>("Total Fractional Uncertainty",
                                                  result,
                                                  result)};
            }
        }

        /// \brief Calculate absolute uncertainty for the given Systematic<T>.
        /// If T is not a histogram, T::Eval(args...)->Hist is invoked
        /// to transform the Systematic<T> to a Systematic<Hist>,
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
        template<class T,
                 class ... Args>
        inline
        std::pair<Hist, Systematic<Hist>>
        AbsoluteUncertainty(const T * nominal_obj,
                            const xsec::Systematic<T> & shifted_obj,
                            Args & ...  args) {
            if constexpr(std::is_same<T, Hist>::value) {
                return _AbsoluteUncertainty(*nominal_obj, shifted_obj);
            } else {
                Hist hnominal = nominal_obj->Eval(std::forward<Args>(args)...);
                Systematic<Hist> hsystematic = EvalSystematic(shifted_obj,
                                                              std::forward<Args>(args)...);
                return _AbsoluteUncertainty(hnominal,
                                            hsystematic);
            }
        }

        /// \brief Calculate fractional uncertainty for the given Systematic<T>.
        /// If T is not a histogram, T::Eval(args...)->Hist is invoked
        /// to transform the Systematic<T> to a Systematic<Hist>,
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
        template<class T,
                 class ... Args>
        inline
        std::pair<Hist, Systematic<Hist>>
        FractionalUncertainty(const T * nominal_obj,
                              const xsec::Systematic<T> & shifted_obj,
                              Args & ...  args) {
            if constexpr(std::is_same<T, Hist>::value) {
                return _FractionalUncertainty(*nominal_obj, shifted_obj);
            } else {
                auto hnominal = nominal_obj->Eval(std::forward<Args>(args)...);
                auto hsystematic = EvalSystematic(shifted_obj,
                                                  std::forward<Args>(args)...);
                return _FractionalUncertainty(hnominal, hsystematic);
            }
        }

        /// \brief Calculate the total absolute uncertainty for the given Systematic<T>'s.
        /// If T's are not histograms, T::Eval(args...)->Hist is invoked
        /// to transform the Systematic<T>'s to Systematic<Hist>'s,
        /// then the following calculations can be performed.
        /// If T is already a histogram, the following calculations are immediately performed.
        ///
        /// For each Systematic<Hist>
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
        template<class T,
                 class ... Args>
        inline
        std::pair<Hist, Systematic<Hist>>
        TotalAbsoluteUncertainty(const T * nominal_obj,
                                 const std::map<std::string,
                                                xsec::Systematic<T> > & shifted_objs,
                                 Args & ...  args) {
            if constexpr(std::is_same<T, Hist>::value) {
                return _TotalAbsoluteUncertainty(*nominal_obj, shifted_objs);
            } else {
                Hist hnominal = nominal_obj->Eval(std::forward<Args>(args)...);
                std::map<std::string, Systematic<Hist>> hsystematics;
                for (auto syst_it = shifted_objs.begin(); syst_it != shifted_objs.end(); syst_it++) {
                    hsystematics[syst_it->first] = EvalSystematic(syst_it->second,
                                                                  std::forward<Args>(args)...);
                }

                return _TotalAbsoluteUncertainty(hnominal, hsystematics);
            }
        }

        /// \brief Calculate the total fractional uncertainty for the given Systematic<T>'s.
        /// If T's are not histograms, T::Eval(args...)->Hist is invoked
        /// to transform the Systematic<T>'s to Systematic<Hist>'s,
        /// then the following calculations can be performed.
        /// If T is already a histogram, the following calculations are immediately performed.
        ///
        /// For each Systematic<Hist>
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
        template<class T,
                 class ... Args>
        inline
        std::pair<Hist, Systematic<Hist>>
        TotalFractionalUncertainty(const T * nominal_obj,
                                   const std::map<std::string,
                                                  xsec::Systematic<T> > & shifted_objs,
                                   Args & ... args) {
            if constexpr(std::is_same<T, Hist>::value) {
                return _TotalFractionalUncertainty(*nominal_obj, shifted_objs);
            } else {
                Hist hnominal = nominal_obj->Eval(std::forward<Args>(args)...);
                std::map<std::string, Systematic<Hist>> hsystematics;
                for (auto syst_it = shifted_objs.begin(); syst_it != shifted_objs.end(); syst_it++) {
                    hsystematics[syst_it->first] = EvalSystematic(syst_it->second,
                                                                  std::forward<Args>(args)...);
                }

                return _TotalFractionalUncertainty(hnominal, hsystematics);
            }
        }
    }
}
