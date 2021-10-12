//
// Created by Derek Doyle on 10/9/21.
//

#include "XSecAna/Systematic.h"
#include "XSecAna/Utils.h"

namespace xsec {
    template<class T>
    Systematic<T>::
    Systematic(std::string name,
               T * shift)
            : fContainer({shift}),
              fType(kOneSided),

              fName(std::move(name)) { fContainer.shrink_to_fit(); }

    template<class T>
    Systematic<T>::
    Systematic(std::string name,
               T * up,
               T * down)
            : fContainer({up, down}),
              fType(kTwoSided),
              fName(std::move(name)) { fContainer.shrink_to_fit(); }

    template<class T>
    Systematic<T>::
    Systematic(std::string name,
               std::vector<T *> & universes)
            : fContainer(universes),
              fType(kMultiverse),
              fName(std::move(name)) { fContainer.shrink_to_fit(); }

    template<class T>
    Systematic<T>::
    Systematic(std::string name,
               std::vector<T *> & container,
               SystType_t type)
            : fContainer(container),
              fType(type),
              fName(std::move(name)) { fContainer.shrink_to_fit(); }

    template<class T>
    SystType_t
    Systematic<T>::
    GetType() const { return fType; }

    template<class T>
    const std::vector<T *> &
    Systematic<T>::
    GetShifts() const { return fContainer; }

    /////////////////////////////////////////////////////////////////////////
    template<class T>
    template<class U>
    Systematic<U>
    Systematic<T>::
    ForEach(ForEachFunction<U, T> for_each, std::string new_name) {
        if (new_name == "") new_name = this->fName;
        std::vector<U *> container(this->fContainer.size());
        for (auto i = 0u; i < this->fContainer.size(); i++) {
            container[i] = for_each(this->fContainer[i]);
        }
        return Systematic<U>(
                new_name, container,
                fType);

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
        }
        else {
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
        std::vector<T *> container(nshifts);
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

        Array shift_arr(root::TH1Props(nominal).nbins_and_uof);
        for (auto ibin = 0u; ibin < shift_arr.size(); ibin++) {
            std::vector<double> vals_c;
            for (auto iuniv = 0u; iuniv < multiverse.GetShifts().size(); iuniv++) {
                vals_c.push_back(multiverse.GetShifts()[iuniv]->GetBinContent(ibin));
            }
            shift_arr(ibin) = BinSigma(nsigma, vals_c, shift_arr(ibin));
        }
        return root::ToROOTLike(nominal, shift_arr);
    }
}

