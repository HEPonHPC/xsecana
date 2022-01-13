#pragma once

#include "TDirectory.h"
#include "TObjString.h"

#include "XSecAna/Type.h"
#include "XSecAna/Utils.h"

#include <exception>
#include <cstdio>
#include <utility>
#include <functional>

namespace xsec {
    enum SystType_t {
        kOneSided,
        kTwoSided,
        kMultiverse,
        kOneOrTwoSided,
    };

    template<class U, class T>
    using ForEachFunction = std::function<std::shared_ptr<U> (const std::shared_ptr<T>)>;


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
                return fMsg;
            }

        private:
            char fMsg[500];
            std::string fExpected;
            std::string fTried;
        };
    }

    template<class T>
    class Systematic {
    public:
        Systematic() = default;

        Systematic(std::string name,
                   std::shared_ptr<T> shift);

        Systematic(std::string name,
                   std::shared_ptr<T> up,
                   std::shared_ptr<T> down);

        Systematic(std::string name,
                   std::vector<std::shared_ptr<T>> & universes);

        Systematic(std::string name,
                   std::vector<std::shared_ptr<T>> & container,
                   SystType_t type);

        /*
        Systematic(const Systematic & syst);
        Systematic(Systematic && syst);
        Systematic & operator=(Systematic && rhs);
        Systematic & operator=(const Systematic & rhs);
*/
        template<class U>
        Systematic<U> ForEach(ForEachFunction<U, T> for_each, std::string new_name = "");

        Systematic<TH1> Eval(const TH1 * data, std::string new_name="") const;

        void SaveTo(TDirectory * dir, const std::string & subdir) const;

        static std::unique_ptr<Systematic<T> > LoadFrom(xsec::type::LoadFunction<T> load,
                                                        TDirectory * dir,
                                                        const std::string & subdir);

        TH1 * CovarianceMatrix(const T * nominal) const;

        const std::vector<std::shared_ptr<T>> & GetShifts() const;

        const std::shared_ptr<T> Up() const;

        const std::shared_ptr<T> Down() const;

        SystType_t GetType() const;

        std::string GetName() const;

    private:
        std::vector<std::shared_ptr<T>> fContainer;
        SystType_t fType;
        std::string fName;
    };

    TH1 * MultiverseShift(Systematic<TH1> multiverse,
                          const TH1 * nominal,
                          double nsigma = 1);
}
