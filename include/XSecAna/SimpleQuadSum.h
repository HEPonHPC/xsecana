#pragma once

#include "XSecAna/Systematic.h"
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
    namespace SimpleQuadSum {
        namespace {
            /// \brief Evaluate all the T's out to TH1s
            /// Calling Eval on the T's the long way since
            /// passing parameter packs through a lambda capture
            /// is pretty opaque in < C++20
            template<class T,
                     class ... Args>
            inline
            Systematic<TH1>
            EvalSystematic(const xsec::Systematic<T> & shifted_obj,
                           Args && ...  args) {
                std::vector<std::shared_ptr<TH1>> shifts(shifted_obj.GetShifts().size());
                for (auto i = 0u; i < shifted_obj.GetShifts().size(); i++) {
                    shifts[i] = shifted_obj.GetShifts()[i]->Eval(std::forward<Args>(args)...);
                }
                return Systematic<TH1>(shifted_obj.GetName(),
                                       shifts,
                                       shifted_obj.GetType());
            }

            /// \brief Internal function for calculating absolute uncertainty
            /// on Systematic<TH1>s
            inline
            std::pair<const TH1 *, Systematic<TH1>>
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
                    down_c = up_c;
                }
                auto max_c = MaxShift(up_c.abs(), down_c.abs());
                auto diff = std::shared_ptr<TH1>(root::ToROOTLike(nominal, max_c, Array::Zero(props.nbins_and_uof)));
                return {nominal,
                        Systematic<TH1>(shifted_obj.GetName(),
                                        diff,
                                        diff)};
            }

            /// \brief Internal function for calculating fractional uncertainty
            /// on Systematic<TH1>s
            std::pair<const TH1 *, Systematic<TH1>>
            _FractionalUncertainty(const TH1 * nominal,
                                   const xsec::Systematic<TH1> & shifted_obj) {
                auto frac = std::get<1>(_AbsoluteUncertainty(nominal, shifted_obj)).Up();
                frac->Divide(nominal);
                // root assumes all histograms are independent when calculating errors.
                // In this case, nominal is in numerator and denominator, and error simplifies
                // to equal the error on the absolute uncertainty
                frac->SetError(Array::Zero(root::TH1Props(nominal).nbins_and_uof).eval().data());
                return {nominal, Systematic<TH1>(shifted_obj.GetName(), frac, frac)};
            }

            /// \brief Internal function for calculating total absolute uncertainty
            /// on Systematic<TH1>s
            std::pair<const TH1 *, Systematic<TH1>>
            _TotalAbsoluteUncertainty(const TH1 * nominal,
                                      const std::map<std::string,
                                                     xsec::Systematic<TH1>> & shifted_objs) {
                std::vector<std::shared_ptr<TH1>> shifts;
                for (auto syst_it = shifted_objs.begin(); syst_it != shifted_objs.end(); syst_it++) {
                    shifts.push_back(std::get<1>(_AbsoluteUncertainty(nominal,
                                                                      syst_it->second)).Up());
                }
                auto result = QuadSum(shifts);
                result->SetError(Array::Zero(root::TH1Props(nominal).nbins_and_uof).eval().data());
                return {nominal, Systematic<TH1>("Total Absolute Uncertainty",
                                                 result,
                                                 result)};
            }

            /// \brief Internal function for calculating total fractional uncertainty
            /// on Systematic<TH1>s
            std::pair<const TH1 *, Systematic<TH1>>
            _TotalFractionalUncertainty(const TH1 * nominal,
                                        const std::map<std::string,
                                                       xsec::Systematic<TH1>> & shifted_objs) {
                std::vector<std::shared_ptr<TH1>> shifts;
                for (auto syst_it = shifted_objs.begin(); syst_it != shifted_objs.end(); syst_it++) {
                    shifts.push_back(std::get<1>(_AbsoluteUncertainty(nominal,
                                                                      syst_it->second)).Up());
                }
                auto result = QuadSum(shifts);
                result->Divide(nominal);
                result->SetError(Array::Zero(root::TH1Props(nominal).nbins_and_uof).eval().data());
                return {nominal, Systematic<TH1>("Total Fractional Uncertainty",
                                                 result,
                                                 result)};
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
        ///    - If type kTwoSided, the returning Systematic is also two sided and symmeterized,
        ///      such that Up() == Down()
        ///    - If type kOneSided, the returning Systematic is also two sided and symmeterized,
        ///      such that Up() == Down()
        /// TODO can we add a switch for asymmetric errors?
        template<class T,
                 class ... Args>
        inline
        std::pair<const TH1 *, Systematic<TH1>>
        AbsoluteUncertainty(const T * nominal_obj,
                            const xsec::Systematic<T> & shifted_obj,
                            Args & ...  args) {
            if constexpr(std::is_same<T, TH1>::value) {
                return _AbsoluteUncertainty(nominal_obj, shifted_obj);
            } else {
                std::shared_ptr<TH1> hnominal = nominal_obj->Eval(std::forward<Args>(args)...);
                Systematic<TH1> hsystematic = EvalSystematic(shifted_obj,
                                                             std::forward<Args>(args)...);
                return _AbsoluteUncertainty(hnominal.get(),
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
        ///    - If type kTwoSided, the returning Systematic is also two sided and symmeterized,
        ///      such that Up() == Down()
        ///    - If type kOneSided, the returning Systematic is also two sided and symmeterized,
        ///      such that Up() == Down()
        template<class T,
                 class ... Args>
        inline
        std::pair<const TH1 *, Systematic<TH1>>
        FractionalUncertainty(const T * nominal_obj,
                              const xsec::Systematic<T> & shifted_obj,
                              Args & ...  args) {
            if constexpr(std::is_same<T, TH1>::value) {
                return _FractionalUncertainty(nominal_obj, shifted_obj);
            } else {
                std::shared_ptr<TH1> hnominal = nominal_obj->Eval(std::forward<Args>(args)...);
                auto hsystematic = EvalSystematic(shifted_obj,
                                                  std::forward<Args>(args)...);
                return _FractionalUncertainty(hnominal.get(), hsystematic);
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
        ///    - If type kTwoSided, the returning Systematic is also two sided and symmeterized,
        ///      such that Up() == Down()
        ///    - If type kOneSided, the returning Systematic is also two sided and symmeterized,
        ///      such that Up() == Down()
        /// then, the individual absolute uncertainties are added in quadrature.
        /// The returned Systematic contains the symmeterized total absolute uncertainty.
        template<class T,
                 class ... Args>
        inline
        std::pair<const TH1 *, Systematic<TH1>>
        TotalAbsoluteUncertainty(const T * nominal_obj,
                                 const std::map<std::string,
                                                xsec::Systematic<T> > & shifted_objs,
                                 Args & ...  args) {
            if constexpr(std::is_same<T, TH1>::value) {
                return _TotalAbsoluteUncertainty(nominal_obj, shifted_objs);
            } else {
                std::shared_ptr<TH1> hnominal = nominal_obj->Eval(std::forward<Args>(args)...);
                std::map<std::string, Systematic<TH1>> hsystematics;
                for (auto syst_it = shifted_objs.begin(); syst_it != shifted_objs.end(); syst_it++) {
                    hsystematics[syst_it->first] = EvalSystematic(syst_it->second,
                                                                  std::forward<Args>(args)...);
                }

                return _TotalAbsoluteUncertainty(hnominal.get(), hsystematics);
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
        ///    - If type kTwoSided, the returning Systematic is also two sided and symmeterized,
        ///      such that Up() == Down()
        ///    - If type kOneSided, the returning Systematic is also two sided and symmeterized,
        ///      such that Up() == Down()
        /// then, the individual fractional uncertainties are added in quadrature.
        /// The returned Systematic contains the symmeterized total fractional uncertainty.
        template<class T,
                 class ... Args>
        inline
        std::pair<const TH1 *, Systematic<TH1>>
        TotalFractionalUncertainty(const T * nominal_obj,
                                   const std::map<std::string,
                                                  xsec::Systematic<T> > & shifted_objs,
                                   Args & ... args) {
            if constexpr(std::is_same<T, TH1>::value) {
                return _TotalFractionalUncertainty(nominal_obj, shifted_objs);
            } else {
                std::shared_ptr<TH1> hnominal = nominal_obj->Eval(std::forward<Args>(args)...);
                std::map<std::string, Systematic<TH1>> hsystematics;
                for (auto syst_it = shifted_objs.begin(); syst_it != shifted_objs.end(); syst_it++) {
                    hsystematics[syst_it->first] = EvalSystematic(syst_it->second,
                                                                  std::forward<Args>(args)...);
                }

                return _TotalFractionalUncertainty(hnominal.get(), hsystematics);
            }
        }
    }
}
