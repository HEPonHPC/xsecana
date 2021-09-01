#pragma once

#include "TDirectory.h"
#include "TObjString.h"

#include "XSecAna/Type.h"

#include <exception>
#include <cstdio>
#include <utility>

namespace xsec {
    enum SystType_t {
        kOneSided,
        kTwoSided,
        kMultiverse,
        kOneOrTwoSided,
    };

    namespace exceptions {
        class SystematicTypeError : public std::exception {
        public:
            SystematicTypeError(const char * caller,
                                SystType_t expected_type,
                                SystType_t tried_type) {
                if (expected_type == kOneSided) fExpected = "kOneSided";
                if (expected_type == kTwoSided) fExpected = "kTwoSided";
                if (expected_type == kMultiverse) fExpected = "kMultiverse";
                if (expected_type == kOneOrTwoSided) fExpected = "kOneOrTwoSided";
                if (tried_type == kOneSided) fTried = "kOneSided";
                if (tried_type == kTwoSided) fTried = "kTwoSided";
                if (tried_type == kMultiverse) fTried = "kMultiverse";
                if (tried_type == kOneOrTwoSided) fTried = "kOneOrTwoSided";
                std::sprintf(fMsg,
                             "Invalid use of %s. Expected type %s, got %s.",
                             caller,
                             fExpected.c_str(),
                             fTried.c_str());
            }

            [[nodiscard]] const char * what() const throw() override {
                return &fMsg[0];
            }

        private:
            char fMsg[500];
            std::string fExpected;
            std::string fTried;
        };
    }
    template< class T >
    class Systematic {
    public:
        Systematic() = default;

        Systematic(std::string name,
                   const T & shift)
                : fContainer({shift}),
                  fType(kOneSided),
                  fName(std::move(name)) { fContainer.shrink_to_fit(); }

        Systematic(std::string name,
                   const T & up,
                   const T & down)
                : fContainer({up, down}),
                  fType(kTwoSided),
                  fName(std::move(name)) { fContainer.shrink_to_fit(); }

        Systematic(std::string name,
                   const std::vector< T > & universes)
                : fContainer(universes),
                  fType(kMultiverse),
                  fName(std::move(name)) { fContainer.shrink_to_fit(); }

        Systematic(const std::vector< T > & container,
                   SystType_t type,
                   std::string name)
                : fContainer(container),
                  fType(type),
                  fName(std::move(name)) { fContainer.shrink_to_fit(); }

        template< class F, class... Args >
        Systematic< std::invoke_result_t< F, T, Args... > > Invoke(F && f, Args && ... args);

        template< class F, class... Args >
        Systematic< std::invoke_result_t< F, T, Args... > > Invoke(F && f, Args && ... args) const;

        void SaveTo(TDirectory * dir, const std::string& subdir) const;

        static std::unique_ptr< Systematic< T > > LoadFrom(TDirectory * dir, const std::string& subdir);

        const std::vector< T > & GetShifts() const { return fContainer; }

        const T & Up() const;

        const T & Down() const;

        [[nodiscard]] SystType_t GetType() const { return fType; }

        [[nodiscard]] std::string GetName() const { return fName; }

        ///\brief calculates nsigma shift from nominal
        /// --If T is not an instantiation of Hist
        ///   create a new MultiverseSystematic<Hist>
        ///   by calling T::ToHist on each universe and the
        ///   nominal argument, and MultiverseSystematic<Hist>::NSigmaShift
        ///   is called
        /// --If T is an instantiation of Hist, BinSigma is called
        ///
        /// In either case, return type is an instantiation
        /// of Hist, determined by return type of T::ToHist, or T itself
        auto NSigmaShift(double nsigma,
                         const T & nominal) const;


    private:
        template< class Scalar >
        Scalar BinSigma(const double & nsigma, std::vector< Scalar > & universes, const Scalar & nominal) const;


        std::vector< T > fContainer;
        SystType_t fType;
        std::string fName;
    };


    /////////////////////////////////////////////////////////////////////////
    template< class T >
    template< class F, class... Args >
    Systematic< std::invoke_result_t< F, T, Args... > >
    Systematic< T >::Invoke(F && f, Args && ... args) {
        std::vector< std::invoke_result_t< F, T, Args... > > container(this->fContainer.size());
        for (auto i = 0u; i < fContainer.size(); i++) {
            container[i] = std::invoke(std::forward< F >(f), fContainer[i], std::forward< Args >(args)...);
        }
        return Systematic< std::invoke_result_t< F, T, Args... > >(container,
                                                                   fType,
                                                                   this->fName);
    }

    /////////////////////////////////////////////////////////////////////////
    template< class T >
    template< class F, class... Args >
    Systematic< std::invoke_result_t< F, T, Args... > >
    Systematic< T >::Invoke(F && f, Args && ... args) const {
        std::vector< std::invoke_result_t< F, T, Args... > > container(this->fContainer.size());
        for (auto i = 0u; i < fContainer.size(); i++) {
            container[i] = std::invoke(std::forward< F >(f), fContainer[i], std::forward< Args >(args)...);
        }
        return Systematic< std::invoke_result_t< F, T, Args... > >(container,
                                                                   fType,
                                                                   this->fName);
    }

    /////////////////////////////////////////////////////////////////////////
    template< class T >
    template< class Scalar >
    Scalar
    Systematic< T >::
    BinSigma(const double & nsigma,
             std::vector< Scalar > & universes,
             const Scalar & nominal) const {
        int pivotbin = 0;
        std::sort(universes.begin(), universes.end());
        for (auto i = 0u; i < universes.size() - 1; i++) {
            if (nominal >= universes[i] && nominal < universes[i + 1]) {
                pivotbin = i;
                break;
            }
        }
        double count_fraction = std::erf(nsigma / std::sqrt(2));

        int nsideevents = 0;
        int lastbinindex = (int) universes.size() - 1;
        if (nsigma >= 0) nsideevents = lastbinindex - pivotbin;
        else nsideevents = pivotbin;
        int boundIdx = pivotbin + 0.5 + count_fraction * (double) nsideevents;

        int index = 0;
        if (nsigma >= 0) index = std::min(boundIdx, (int) universes.size() - 1);
        else index = std::max(boundIdx, 0);
        return universes.at(index);
    }

    /////////////////////////////////////////////////////////////////////////
    template< class T >
    auto
    Systematic< T >::
    NSigmaShift(double nsigma,
                const T & nominal) const {
        if (fType != kMultiverse) {
            throw exceptions::SystematicTypeError(__PRETTY_FUNCTION__, kMultiverse, fType);
        }
        if constexpr(type::IsHist< T >()) {
            T shift = nominal;

            for (auto ibin = 0u; ibin < nominal.size(); ibin++) {
                std::vector< typename T::scalar_type > vals;
                for (auto iuniv = 0u; iuniv < fContainer.size(); iuniv++) {
                    vals.push_back(fContainer[iuniv][ibin]);
                }
                shift[ibin] = BinSigma(nsigma, vals, nominal[ibin]);
            }
            return shift;
        } else {
            return this->Invoke(&T::ToHist).NSigmaShift(nsigma, nominal.ToHist());
        }
    }

    /////////////////////////////////////////////////////////////////////////
    template< class T >
    const T &
    Systematic< T >::
    Up() const {
        if (fType == kMultiverse) {
            throw exceptions::SystematicTypeError(__PRETTY_FUNCTION__, kOneOrTwoSided, fType);
        } else return this->fContainer[0];
    }

    /////////////////////////////////////////////////////////////////////////
    template< class T >
    const T &
    Systematic< T >::
    Down() const {
        if (fType == kMultiverse) {
            throw exceptions::SystematicTypeError(__PRETTY_FUNCTION__, kOneOrTwoSided, fType);
        }
        if (fType == kOneSided) return this->fContainer[0];
        else return this->fContainer[1];
    }

    /////////////////////////////////////////////////////////////////////////
    template< class T >
    void
    Systematic< T >::SaveTo(TDirectory * dir, const std::string& subdir) const {
        TDirectory * tmp = gDirectory;
        dir->mkdir(subdir.c_str());
        dir = dir->GetDirectory(subdir.c_str());
        dir->cd();

        TObjString("Systematic").Write("type");
        TObjString(this->fName.c_str()).Write("fName");
        TObjString(std::to_string(this->fContainer.size()).c_str()).Write("NShifts");
        TObjString(std::to_string(this->fType).c_str()).Write("fType");

        auto mv_dir = dir->mkdir("fContainer");
        for (auto i = 0u; i < fContainer.size(); i++) {
            fContainer[i].SaveTo(mv_dir, std::to_string(i));
        }

        tmp->cd();
    }

    /////////////////////////////////////////////////////////////////////////
    template< class T >
    std::unique_ptr< Systematic< T > >
    Systematic< T >::LoadFrom(TDirectory * dir, const std::string& subdir) {
        TDirectory * tmp = gDirectory;
        dir = dir->GetDirectory(subdir.c_str());
        dir->cd();

        // make sure we're loading the right type
        auto ptag = (TObjString *) dir->Get("type");
        assert(ptag->GetString() == "Systematic" && "Type does not match Systematic");
        delete ptag;

        std::string name = ((TObjString *) dir->Get("fName"))->GetString().Data();
        int nshifts = std::atoi(((TObjString *) dir->Get("NShifts"))->GetString().Data());
        int ftype = std::atoi(((TObjString *) dir->Get("fType"))->GetString().Data());

        auto container_dir = dir->GetDirectory("fContainer");
        std::vector< T > container(nshifts);
        for (auto ishift = 0; ishift < nshifts; ishift++) {
            container[ishift] = *T::LoadFrom(container_dir, std::to_string(ishift)).release();
        }

        tmp->cd();
        return std::make_unique< Systematic< T > >(container,
                                                   (SystType_t) ftype,
                                                   name);
    }
}
