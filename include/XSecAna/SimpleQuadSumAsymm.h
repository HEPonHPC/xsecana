#pragma once

#include "XSecAna/Systematic.h"
#include "XSecAna/SimpleQuadSum.h"
#include "XSecAna/Math.h"

#include <Eigen/Dense>
#include "XSecAna/Type.h"

namespace xsec {
    /// \brief
    /// SimpleQuadSum performs the quadrature sum of systematic shifts
    /// In the case of an asymmetric 2-sided shift,
    /// the shift is symmeterized by taking the largest shift
    /// Type T must have implemented T::Eval(Args ... args)->TH1
    /// in order to transform into a TH1.
    /// The only exception is if type T is type TH1.
    namespace SimpleQuadSumAsymm {
        namespace detail {

            /// \brief Internal function for calculating absolute uncertainty
            /// on Systematic<TH1>s
            inline
            std::pair<const TH1 *, Systematic < TH1>>
            _AbsoluteUncertainty(const TH1 * nominal,
                                 const xsec::Systematic<TH1> & shifted_obj) {
                std::vector<TH1 *> shifts(shifted_obj.GetShifts().size());
                root::TH1Props props(nominal);
                Array up_c(props.nbins_and_uof);
                Array down_c(props.nbins_and_uof);

                Array nom_c = root::MapContentsToEigen(nominal);
                // convert multiverse systematic to two-sided by finding 1sigma
                if (shifted_obj.GetType() == kMultiverse) {
                    up_c = root::MapContentsToEigen(MultiverseShift(shifted_obj, nominal, 1)) - nom_c;
                    down_c = root::MapContentsToEigen(MultiverseShift(shifted_obj, nominal, -1)) - nom_c;
                } else if (shifted_obj.GetType() == kTwoSided) {
                    up_c = root::MapContentsToEigen(shifted_obj.GetShifts()[0].get());
                    down_c = root::MapContentsToEigen(shifted_obj.GetShifts()[1].get());
                    up_c = up_c - nom_c;
                    down_c = down_c - nom_c;
                } else {
                    up_c = root::MapContentsToEigen(shifted_obj.GetShifts()[0].get());
                    up_c = up_c - nom_c;
                }

                auto dup = std::shared_ptr<TH1>(root::ToROOTLike(nominal, up_c, Array::Zero(props.nbins_and_uof)));
                if(shifted_obj.GetType() == kOneSided) {
                    return {nominal,
                            Systematic<TH1>(shifted_obj.GetName(), dup)};
                }
                else {
                    auto ddw = std::shared_ptr<TH1>(
                            root::ToROOTLike(nominal, down_c, Array::Zero(props.nbins_and_uof)));
                    return {nominal,
                            Systematic<TH1>(shifted_obj.GetName(),
                                            dup,
                                            ddw)};
                }
            }

            /// \brief Internal function for calculating fractional uncertainty
            /// on Systematic<TH1>s
            inline
            std::pair<const TH1 *, Systematic < TH1>>
            _FractionalUncertainty(const TH1 * nominal,
                                   const xsec::Systematic<TH1> & shifted_obj) {
                Systematic<TH1> abs_uncert = std::get<1>(_AbsoluteUncertainty(nominal, shifted_obj));
                auto frac_up = abs_uncert.Up();
                frac_up->Divide(nominal);

                // root assumes all histograms are independent when calculating errors.
                // In this case, nominal is in numerator and denominator, and error simplifies
                // to equal the error on the absolute uncertainty
                frac_up->SetError(Array::Zero(root::TH1Props(nominal).nbins_and_uof).eval().data());
                if(abs_uncert.GetType() == kOneSided) {
                    return {nominal, Systematic<TH1>(shifted_obj.GetName(), frac_up)};
                }
                else {
                    auto frac_dw = abs_uncert.Down();
                    frac_dw->Divide(nominal);
                    frac_dw->SetError(Array::Zero(root::TH1Props(nominal).nbins_and_uof).eval().data());
                    return {nominal, Systematic<TH1>(shifted_obj.GetName(), frac_up, frac_dw)};
                }
            }

            /// \brief Internal function for calculating total absolute uncertainty
            /// on Systematic<TH1>s
            inline
            std::pair<const TH1 *, Systematic < TH1>>
            _TotalAbsoluteUncertainty(const TH1 * nominal,
                                      const std::map<std::string,
                                      xsec::Systematic<TH1>> & shifted_objs) {
                std::vector<std::shared_ptr<TH1>> shifts_up;
                std::vector<std::shared_ptr<TH1>> shifts_dw;
                for (auto syst_it = shifted_objs.begin(); syst_it != shifted_objs.end(); syst_it++) {
                    Systematic<TH1> abs_uncert = std::get<1>(_AbsoluteUncertainty(nominal, syst_it->second));
                    shifts_up.push_back(abs_uncert.Up());
                    if(abs_uncert.GetType() == kTwoSided) {
                        shifts_dw.push_back(std::get<1>(_AbsoluteUncertainty(nominal,
                                                                             syst_it->second)).Down());
                    }
                }
                auto result_up = QuadSum(shifts_up);
                result_up->SetError(Array::Zero(root::TH1Props(nominal).nbins_and_uof).eval().data());
                if(!shifts_dw.empty()) {
                    auto result_dw = QuadSum(shifts_dw);
                    result_dw->SetError(Array::Zero(root::TH1Props(nominal).nbins_and_uof).eval().data());
                    return {nominal, Systematic<TH1>("Total Absolute Uncertainty",
                                                     result_up,
                                                     result_dw)};
                }
                else {
                    return {nominal, Systematic<TH1>("Total Absolute Uncertainty",
                                                     result_up)};
                }
            }

            /// \brief Internal function for calculating total fractional uncertainty
            /// on Systematic<TH1>s
            inline
            std::pair<const TH1 *, Systematic < TH1>>
            _TotalFractionalUncertainty(const TH1 * nominal,
                                        const std::map<std::string,
                                                       xsec::Systematic<TH1>> & shifted_objs) {
                std::vector<std::shared_ptr<TH1>> shifts_up;
                std::vector<std::shared_ptr<TH1>> shifts_dw;
                for (auto syst_it = shifted_objs.begin(); syst_it != shifted_objs.end(); syst_it++) {
                    Systematic<TH1> abs_uncert = std::get<1>(_AbsoluteUncertainty(nominal, syst_it->second));
                    shifts_up.push_back(abs_uncert.Up());
                    if(abs_uncert.GetType() == kTwoSided) {
                        shifts_dw.push_back(abs_uncert.Down());
                    }
                }
                auto result_up = QuadSum(shifts_up);
                result_up->Divide(nominal);
                result_up->SetError(Array::Zero(root::TH1Props(nominal).nbins_and_uof).eval().data());
                if(!shifts_dw.empty()) {
                    auto result_dw = QuadSum(shifts_dw);
                    result_dw->Divide(nominal);
                    result_dw->SetError(Array::Zero(root::TH1Props(nominal).nbins_and_uof).eval().data());
                    return {nominal, Systematic<TH1>("Total Fractional Uncertainty",
                                                     result_up,
                                                     result_dw)};
                }
                else {
                    return {nominal, Systematic<TH1>("Total Fractional Uncertainty",
                                                     result_up)};
                }
            }
        }

        /// \brief Calculate absolute uncertainty for the given Systematic<T>.
        /// If T is not a histogram, T::Eval(args...)->TH1 is invoked
        /// to transform the Systematic<T> to a Systematic<TH1>,
        /// then the following calculations can be performed.
        /// If T is already a histogram, the following calculations are immediately performed.
        ///
        ///  - If type kMultiverse, the 1sigma universe is found.
        ///  - Otherwise, it is assumed the Systematic represents a 1sigma shift,
        ///    and as such, the absolute uncertainty is the difference between the Systematic
        ///    and nominal.
        template<class T,
                 class ... Args>
        inline
        std::pair<const TH1 *, Systematic < TH1>>
        AbsoluteUncertainty(const T * nominal_obj,
                            const xsec::Systematic<T> & shifted_obj,
                            Args & ...  args) {
            if constexpr(std::is_same<T, TH1>::value) {
                return detail::_AbsoluteUncertainty(nominal_obj, shifted_obj);
            } else {
                std::shared_ptr<TH1> hnominal = nominal_obj->Eval(std::forward<Args>(args)...);
                Systematic<TH1> hsystematic = SimpleQuadSum::detail::EvalSystematic(shifted_obj,
                                                                                    std::forward<Args>(args)...);
                return detail::_AbsoluteUncertainty(hnominal.get(),
                                                    hsystematic);
            }
        }

        /// \brief Calculate fractional uncertainty for the given Systematic<T>.
        /// If T is not a histogram, T::Eval(args...)->TH1 is invoked
        /// to transform the Systematic<T> to a Systematic<TH1>,
        /// then the following calculations can be performed.
        /// If T is already a histogram, the following calculations are immediately performed.
        ///
        ///  - If type kMultiverse, the 1sigma universe is found.
        ///  - Otherwise, it is assumed the Systematic represents a 1sigma shift,
        ///    and as such, the fractional is the difference between the Systematic
        ///    and nominal divided by nominal.
        template<class T,
                 class ... Args>
        inline
        std::pair<const TH1 *, Systematic < TH1>>

        FractionalUncertainty(const T * nominal_obj,
                              const xsec::Systematic<T> & shifted_obj,
                              Args & ...  args) {
            if constexpr(std::is_same<T, TH1>::value) {
                return detail::_FractionalUncertainty(nominal_obj, shifted_obj);
            } else {
                std::shared_ptr<TH1> hnominal = nominal_obj->Eval(std::forward<Args>(args)...);
                auto hsystematic = SimpleQuadSum::detail::EvalSystematic(shifted_obj,
                                                                         std::forward<Args>(args)...);
                return detail::_FractionalUncertainty(hnominal.get(), hsystematic);
            }
        }

        /// \brief Calculate the total absolute uncertainty for the given Systematic<T>'s.
        /// If T's are not histograms, T::Eval(args...)->TH1 is invoked
        /// to transform the Systematic<T>'s to Systematic<TH1>'s,
        /// then the following calculations can be performed.
        /// If T is already a histogram, the following calculations are immediately performed.
        ///
        /// For each Systematic<TH1>
        ///  - If type kMultiverse, the 1sigma universe is found.
        ///  - Otherwise, it is assumed the Systematic represents a 1sigma shift,
        ///    and as such, the absolute is the difference between the Systematic
        ///    and nominal.
        /// then, the individual absolute uncertainties are added in quadrature.
        /// The returned Systematic contains the asymmetrical total absolute uncertainty.
        template<class T,
                 class ... Args>
        inline
        std::pair<const TH1 *, Systematic < TH1>>

        TotalAbsoluteUncertainty(const T * nominal_obj,
                                 const std::map<std::string,
                                                xsec::Systematic<T> > & shifted_objs,
                                 Args & ...  args) {
            if constexpr(std::is_same<T, TH1>::value) {
                return detail::_TotalAbsoluteUncertainty(nominal_obj, shifted_objs);
            } else {
                std::shared_ptr<TH1> hnominal = nominal_obj->Eval(std::forward<Args>(args)...);
                std::map<std::string, Systematic<TH1>> hsystematics;
                for (auto syst_it = shifted_objs.begin(); syst_it != shifted_objs.end(); syst_it++) {
                    hsystematics[syst_it->first] = SimpleQuadSum::detail::EvalSystematic(syst_it->second,
                                                                                         std::forward<Args>(args)...);
                }

                return detail::_TotalAbsoluteUncertainty(hnominal.get(), hsystematics);
            }
        }

        /// \brief Calculate the total fractional uncertainty for the given Systematic<T>'s.
        /// If T's are not histograms, T::Eval(args...)->TH1 is invoked
        /// to transform the Systematic<T>'s to Systematic<TH1>'s,
        /// then the following calculations can be performed.
        /// If T is already a histogram, the following calculations are immediately performed.
        ///
        /// For each Systematic<TH1>
        ///  - If type kMultiverse, the 1sigma universe is found.
        ///  - Otherwise, it is assumed the Systematic represents a 1sigma shift,
        ///    and as such, the fractional is the difference between the Systematic
        ///    and nominal divided by nominal.
        /// then, the individual fractional uncertainties are added in quadrature.
        /// The returned Systematic contains the asymmetrical total fractional uncertainty.
        template<class T,
                 class ... Args>
        inline
        std::pair<const TH1 *, Systematic < TH1>>

        TotalFractionalUncertainty(const T * nominal_obj,
                                   const std::map<std::string,
                                                  xsec::Systematic<T> > & shifted_objs,
                                   Args & ... args) {
            if constexpr(std::is_same<T, TH1>::value) {
                return detail::_TotalFractionalUncertainty(nominal_obj, shifted_objs);
            } else {
                std::shared_ptr<TH1> hnominal = nominal_obj->Eval(std::forward<Args>(args)...);
                std::map<std::string, Systematic<TH1>> hsystematics;
                for (auto syst_it = shifted_objs.begin(); syst_it != shifted_objs.end(); syst_it++) {
                    hsystematics[syst_it->first] = SimpleQuadSum::detail::EvalSystematic(syst_it->second,
                                                                                         std::forward<Args>(args)...);
                }

                return detail::_TotalFractionalUncertainty(hnominal.get(), hsystematics);
            }
        }
    }
}
