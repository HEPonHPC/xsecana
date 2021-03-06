#pragma once

#include "TDirectory.h"
#include "TString.h"

#include "XSecAna/Utils.h"
#include "XSecAna/Type.h"

#include <vector>
#include <typeinfo>

namespace xsec {
    typedef Eigen::ArrayXd Array;
    typedef Eigen::Ref<Array> ArrayRef;

    class IMeasurement {
    public:
        virtual std::shared_ptr<TH1> Eval(const TH1 * data) const = 0;
        virtual void SaveTo(TDirectory * dir, const std::string & subdir) const {
            throw std::runtime_error(std::string(typeid(this).name()) + "::SaveTo not implemented");
        }
        virtual ~IMeasurement() = default;
        static std::unique_ptr<IMeasurement>
        LoadFrom(xsec::type::LoadFunction<IMeasurement> load,
                 TDirectory * dir,
                 const std::string & name) {
            return load(dir, name);
        }

        template<class T>
        T *
        cast() {
            static_assert(std::is_base_of_v<IMeasurement, T>);
            return dynamic_cast<T*>(this);
        }
    };
    typedef xsec::type::LoadFunction<IMeasurement> LoadMeasurementFunc;

    class IEigenEval : public virtual IMeasurement {
    public:
        std::shared_ptr<TH1> Eval(const TH1 * data) const final {
            fHistProps = root::TH1Props(data,
                                        root::MakeUnique("UniqueEval").c_str());

            Array _data(fHistProps.nbins_and_uof);
            Array _error(fHistProps.nbins_and_uof);
            root::MapToEigen(data, _data, _error);

            Array _result(fHistProps.nbins_and_uof);
            Array _rerror(fHistProps.nbins_and_uof);
            this->_eval_impl(_data, _error,
                             _result, _rerror);
            return std::shared_ptr<TH1>(root::ToROOT(_result, _rerror, fHistProps));

        }
        virtual ~IEigenEval()= default;
        virtual void _eval_impl(const Array & data, const Array & error,
                                ArrayRef result, ArrayRef rerror) const = 0;
        const root::TH1Props & GetHistProps() const { return fHistProps; }
    private:
        mutable root::TH1Props fHistProps;
    };
}
