//
// Created by Derek Doyle on 10/9/21.
//

#include "XSecAna/Systematic.h"
#include "XSecAna/IMeasurement.h"
#include "XSecAna/Utils.h"

#include <stdexcept>

namespace xsec {
    template<class T>
    Systematic<T>::
    Systematic(std::string name,
               const T * shift)
            : fContainer({shift}),
              fType(kOneSided),

              fName(std::move(name)) { fContainer.shrink_to_fit(); }

    template<class T>
    Systematic<T>::
    Systematic(std::string name,
               const T * up,
               const T * down)
            : fContainer({up, down}),
              fType(kTwoSided),
              fName(std::move(name)) { fContainer.shrink_to_fit(); }

    template<class T>
    Systematic<T>::
    Systematic(std::string name,
               std::vector<const T *> & universes)
            : fContainer(universes),
              fType(kMultiverse),
              fName(std::move(name)) { fContainer.shrink_to_fit(); }

    template<class T>
    Systematic<T>::
    Systematic(std::string name,
               std::vector<const T *> & container,
               SystType_t type)
            : fContainer(container),
              fType(type),
              fName(std::move(name)) { fContainer.shrink_to_fit(); }

    template<class T>
    SystType_t
    Systematic<T>::
    GetType() const { return fType; }

    template<class T>
    std::string
    Systematic<T>::
    GetName() const { return fName; }

    template<class T>
    const std::vector<const T *> &
    Systematic<T>::
    GetShifts() const { return fContainer; }

    /////////////////////////////////////////////////////////////////////////
    template<class T>
    template<class U>
    Systematic<U>
    Systematic<T>::
    ForEach(ForEachFunction<U, T> for_each, std::string new_name) {
        if (new_name == "") new_name = this->fName;
        std::vector<const U *> container(this->fContainer.size());
        for (auto i = 0u; i < this->fContainer.size(); i++) {
            container[i] = for_each(this->fContainer[i]);
        }
        return Systematic<U>(new_name, container, fType);

    }

    template<class T>
    Matrix
    Systematic<T>::
    CovarianceMatrix(const T * nominal) const {
        if constexpr (!std::is_base_of<TH1, T>::value) {
            throw std::runtime_error("Type " +
                                     std::string(typeid(T).name()) +
                                     " does not implement CovarianceMatrix. Must be of type Systematic<TH1>.");
        } else {
            Array nom_a = root::MapContentsToEigen(nominal);
            std::vector<Array> shifts_a(fContainer.size());
            for(auto i = 0u; i < fContainer.size(); i++) {
                shifts_a[i] = root::MapContentsToEigen(fContainer[i]);
            }

            Matrix cov(nom_a.size(), nom_a.size());
            if(fType == kOneSided) {
                for (auto i = 0; i < nom_a.size(); i++) {
                    for (auto j = 0; j < nom_a.size(); j++) {
                        cov(i, j) = (nom_a(i) - shifts_a[0](i)) *
                                    (nom_a(j) - shifts_a[0](j));
                    }
                }
            }
            else if(fType == kTwoSided) {
                for (auto i = 0; i < nom_a.size(); i++) {
                    for (auto j = 0; j < nom_a.size(); j++) {
                        auto di = nom_a(i) - std::max(shifts_a[0](i), shifts_a[1](i));
                        auto dj = nom_a(j) - std::max(shifts_a[0](j), shifts_a[1](j));
                        cov(i, j) = di * dj;
                    }
                }
            }
            else { // multiverse
                 for (auto i = 0; i < nom_a.size(); i++) {
                    for (auto j = 0; j < nom_a.size(); j++) {
                        double v = 0;
                        for(auto u = 0u; u < shifts_a.size(); u++) {
                            v += (shifts_a[u](i) - nom_a(i)) *
                                 (shifts_a[u](j) - nom_a(j));
                        }
                        cov(i, j) = v / shifts_a.size();
                    }
                }
            }
            return cov;
        }
    }

    template<class T>
    Systematic<TH1>
    Systematic<T>::
    Eval(const TH1 * data, std::string new_name) const {
        if constexpr (!std::is_base_of<IMeasurement, T>::value) {
            throw std::runtime_error("Type " +
                                     std::string(typeid(T).name()) +
                                     " does not implement Eval. Use ForEach instead");
        } else {
            if (new_name == "") new_name = this->fName;
            std::vector<const TH1 *> container(this->fContainer.size());
            for (auto i = 0u; i < this->fContainer.size(); i++) {
                container[i] = this->fContainer[i]->Eval(data);
            }
            return Systematic<TH1>(new_name, container, fType);
        }
    }


    /////////////////////////////////////////////////////////////////////////
    template<class T>
    const T *
    Systematic<T>::
    Up() const {
        if (fType == kMultiverse) {
            throw exceptions::SystematicTypeError(__PRETTY_FUNCTION__, kOneOrTwoSided, fType);
        } else return this->fContainer[0];
    }

    /////////////////////////////////////////////////////////////////////////
    template<class T>
    const T *
    Systematic<T>::
    Down() const {
        if (fType == kMultiverse || fType == kOneSided) {
            throw exceptions::SystematicTypeError(__PRETTY_FUNCTION__, kTwoSided, fType);
        } else return this->fContainer[1];
    }

    /////////////////////////////////////////////////////////////////////////
    template<class T>
    void
    Systematic<T>::SaveTo(TDirectory * dir, const std::string & subdir) const {
        TDirectory * tmp = gDirectory;
        dir->mkdir(subdir.c_str());
        dir = dir->GetDirectory(subdir.c_str());
        dir->cd();

        TObjString("Systematic").Write("type");
        TObjString(this->fName.c_str()).Write("fName");
        TObjString(std::to_string(this->fContainer.size()).c_str()).Write("NShifts");
        TObjString(std::to_string(this->fType).c_str()).Write("fType");

        auto mv_dir = dir->mkdir("fContainer");
        if constexpr(std::is_same<T, TH1>::value) {
            mv_dir->cd();
            for (auto i = 0u; i < fContainer.size(); i++) {
                fContainer[i]->Write(std::to_string(i).c_str());
            }
        } else if constexpr(std::is_same<T, Array>::value) {
            auto tmp = new TH1D("", "",
                                fContainer[0]->size(),
                                0, fContainer[0]->size());
            root::TH1Props props(tmp);
            for (auto i = 0u; i < fContainer.size(); i++) {
                root::ToROOT(*fContainer[i], props)->Write(std::to_string(i).c_str());
            }
        } else {
            for (auto i = 0u; i < fContainer.size(); i++) {
                fContainer[i]->SaveTo(mv_dir, std::to_string(i));
            }
        }

        tmp->cd();
    }

    /////////////////////////////////////////////////////////////////////////
    template<class T>
    std::unique_ptr<Systematic<T> >
    Systematic<T>::LoadFrom(xsec::type::LoadFunction<T> load,
                            TDirectory * dir,
                            const std::string & subdir) {
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
        std::vector<const T *> container(nshifts);
        for (auto ishift = 0; ishift < nshifts; ishift++) {
            container[ishift] = load(container_dir, std::to_string(ishift)).release();
        }

        tmp->cd();
        return std::make_unique<Systematic<T> >(name,
                                                container,
                                                (SystType_t) ftype);
    }

    /////////////////////////////////////////////////////////////////////////
    template<class Scalar>
    Scalar
    BinSigma(const double & nsigma,
             std::vector<Scalar> & universes,
             const Scalar & nominal) {
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
    TH1 *
    MultiverseShift(Systematic<TH1> multiverse,
                    const TH1 * nominal,
                    double nsigma) {
        if (multiverse.GetType() != kMultiverse) {
            throw exceptions::SystematicTypeError(__PRETTY_FUNCTION__,
                                                  kMultiverse,
                                                  multiverse.GetType());
        }
        root::TH1Props props(nominal);
        Array shift_arr(props.nbins_and_uof);
        for (auto ibin = 0u; ibin < shift_arr.size(); ibin++) {
            std::vector<double> vals_c(multiverse.GetShifts().size());
            for (auto iuniv = 0u; iuniv < multiverse.GetShifts().size(); iuniv++) {
                vals_c[iuniv] = multiverse.GetShifts()[iuniv]->GetBinContent(ibin);
            }
            shift_arr(ibin) = BinSigma(nsigma, vals_c, nominal->GetBinContent(ibin));
        }
        return root::ToROOTLike(nominal, shift_arr, Array::Zero(props.nbins_and_uof));
    }

    template
    class Systematic<Array>;

    template
    class Systematic<TH1>;

    template
    class Systematic<IMeasurement>;

    template
    Systematic<TH1>
    Systematic<TH1>::
    ForEach(ForEachFunction<TH1, TH1>, std::string);

    template
    Systematic<TH1>
    Systematic<IMeasurement>::
    ForEach(ForEachFunction<TH1, IMeasurement>, std::string);

}

