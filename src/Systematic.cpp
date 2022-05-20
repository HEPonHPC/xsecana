//
// Created by Derek Doyle on 10/9/21.
//

#include "XSecAna/Systematic.h"
#include "XSecAna/IMeasurement.h"
#include "XSecAna/Utils.h"

#include <stdexcept>
#include <random>

namespace xsec {
    template<class T>
    Systematic<T>::
    Systematic(std::string name,
               std::shared_ptr<T> shift)
            : fContainer({shift}),
              fType(kOneSided),

              fName(std::move(name)) { fContainer.shrink_to_fit(); }

    template<class T>
    Systematic<T>::
    Systematic(std::string name,
               std::shared_ptr<T> up,
               std::shared_ptr<T> down)
            : fContainer({up, down}),
              fType(kTwoSided),
              fName(std::move(name)) { fContainer.shrink_to_fit(); }

    template<class T>
    Systematic<T>::
    Systematic(std::string name,
               std::vector<std::shared_ptr<T>> & universes)
            : fContainer(universes),
              fType(kMultiverse),
              fName(std::move(name)) { fContainer.shrink_to_fit(); }

    template<class T>
    Systematic<T>::
    Systematic(std::string name,
               std::vector<std::shared_ptr<T>> & container,
               SystType_t type)
            : fContainer(container),
              fType(type),
              fName(std::move(name)) { fContainer.shrink_to_fit(); }

    /*
    template<class T>
    Systematic<T>::
    Systematic(const Systematic<T> & syst) {
        std::cout << __PRETTY_FUNCTION__ << std::endl;
        fContainer = syst.fContainer;
        fType = syst.fType;
        fName = syst.fName;
    }

    template<class T>
    Systematic<T>::
    Systematic(Systematic<T> && syst) {
        std::cout << __PRETTY_FUNCTION__ << std::endl;
        fContainer = std::move(syst.fContainer);
        fType = syst.fType;
        fName = syst.fName;
    }

    template<class T>
    Systematic<T> &
    Systematic<T>::
    operator=(Systematic<T> && syst) {
        std::cout << __PRETTY_FUNCTION__ << std::endl;
        fContainer = std::move(syst.fContainer);
        fType = syst.fType;
        fName = syst.fName;
        return *this;
    }

    template<class T>
    Systematic<T> &
    Systematic<T>::
    operator=(const Systematic<T> & syst) {
        std::cout << __PRETTY_FUNCTION__ << std::endl;
        fContainer = syst.fContainer;
        fType = syst.fType;
        fName = syst.fName;
        return *this;
    }
     */
    template<class T>
    SystType_t
    Systematic<T>::
    GetType() const { return fType; }

    template<class T>
    std::string
    Systematic<T>::
    GetName() const { return fName; }

    template<class T>
    const std::vector<std::shared_ptr<T>> &
    Systematic<T>::
    GetShifts() const { return fContainer; }

    /////////////////////////////////////////////////////////////////////////
    template<class T>
    template<class U>
    Systematic<U>
    Systematic<T>::
    ForEach(ForEachFunction<U, T> for_each, std::string new_name) {
        if (new_name == "") new_name = this->fName;
        std::vector<std::shared_ptr<U>> container(this->fContainer.size());
        for (auto i = 0u; i < this->fContainer.size(); i++) {
            container[i] = for_each(this->fContainer[i]);
        }
        return Systematic<U>(new_name, container, fType);

    }

    ///\brief Return a random shifted distribution by sampling
    /// the covariance matrix calculated from this systematic shift
    template<class T>
    TH1 *
    Systematic<T>::
    RandomSample(const T * nominal, double seed) const {
        if constexpr(!std::is_base_of<TH1, T>::value) {
            throw std::runtime_error("Type " +
                                     std::string(typeid(T).name()) +
                                     " does not implement RandomSample. Must be of type Systematic<TH1>.");
        } else {
            return Systematic<TH1>::RandomSample(nominal, this->CovarianceMatrix(nominal), seed);
        }
    }

    ///\brief Return a random shifted distribution by sampling
    /// a covariance matrix
    /// https://juanitorduz.github.io/multivariate_normal/
    template<class T>
    TH1 *
    Systematic<T>::
    RandomSample(const T * nominal, const TH1 * covariance, double seed) {
        if constexpr(!std::is_base_of<TH1, T>::value) {
            throw std::runtime_error("Type " +
                                     std::string(typeid(T).name()) +
                                     " does not implement RandomSample. Must be of type Systematic<TH1>.");
        } else {
            Matrix cov = root::MapContentsToEigen(covariance)
                    .reshaped(nominal->GetNbinsX()+2,
                              nominal->GetNbinsX()+2);
            Matrix epsilon = (Vector::Ones(cov.rows()) * 0.0001).asDiagonal();
            cov = cov + epsilon;

            Matrix L = cov.llt().matrixL();
            Vector m = root::MapContentsToEigen(nominal);

            std::mt19937 generator(seed);
            std::normal_distribution<double> distribution(0, 1);
            Vector u(m.size());
            for(auto i = 0u; i < u.size(); i++) {
                u(i) = distribution(generator);
            }

            Vector x = m + L * u;
            return root::ToROOTLike(nominal, x);
        }
    }

    template<class T>
    TH1 *
    Systematic<T>::
    CovarianceMatrix(const T * nominal) const {
        if constexpr (!std::is_base_of<TH1, T>::value) {
            throw std::runtime_error("Type " +
                                     std::string(typeid(T).name()) +
                                     " does not implement CovarianceMatrix. Must be of type Systematic<TH1>.");
        } else {
            std::vector<Array> shifts(fContainer.size());
            for(auto i = 0u; i < fContainer.size(); i++) {
                shifts[i] = root::MapContentsToEigen(fContainer[i].get());
            }

            Array nom_a = root::MapContentsToEigen(nominal);
            Matrix cov(nom_a.size(), nom_a.size());
            if(fType == kOneSided) {
                for (auto i = 0; i < nom_a.size(); i++) {
                    for (auto j = 0; j < nom_a.size(); j++) {
                        cov(i, j) = (nom_a(i) - shifts[0](i)) *
                                    (nom_a(j) - shifts[0](j));
                    }
                }
            }
            else if(fType == kTwoSided) {
                for (auto i = 0; i < nom_a.size(); i++) {
                    for (auto j = 0; j < nom_a.size(); j++) {

                        auto di_0 = nom_a(i) - shifts[0](i);
                        auto di_1 = nom_a(i) - shifts[1](i);
                        auto dj_0 = nom_a(j) - shifts[0](j);
                        auto dj_1 = nom_a(j) - shifts[1](j);
                        auto di = (std::abs(di_0) > std::abs(di_1)) ? di_0 : di_1;
                        auto dj = (std::abs(dj_0) > std::abs(dj_1)) ? dj_0 : dj_1;

                        //cov(i, j) = di * dj;
                        cov(i, j) = di_0 * dj_0 + di_1 * dj_1;
                    }
                }
            }
            else { // multiverse
                std::vector<double> multiverse_means(nom_a.size(), 0);
                for(auto ibin = 0u; ibin < nom_a.size(); ibin++) {
                    for(auto imv = 0u; imv < shifts.size(); imv++) {
                        multiverse_means[ibin] += shifts[imv][ibin];
                    }
                    multiverse_means[ibin] /= shifts.size();
                }

                 for (auto i = 0; i < nom_a.size(); i++) {
                    for (auto j = 0; j < nom_a.size(); j++) {

                        double v = 0;
                        for(auto u = 0u; u < shifts.size(); u++) {
                            //v += (shifts[u](i) - nom_a(i)) *
                            //     (shifts[u](j) - nom_a(j));
                            v += (shifts[u](i) - multiverse_means[i]) *
                                 (shifts[u](j) - multiverse_means[j]);
                        }
                        cov(i, j) =  v / (fContainer.size()-1);
                    }
                }
            }
            auto ret = new TH2D("", "",
                                nom_a.size()-2, 0, nom_a.size()-2,
                                nom_a.size()-2, 0, nom_a.size()-2);
            ret->SetContent(cov.data());
            ret->SetEntries(cov.size());
            return ret;
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
            std::vector<std::shared_ptr<TH1>> container(this->fContainer.size());
            for (auto i = 0u; i < this->fContainer.size(); i++) {
                container[i] = this->fContainer[i]->Eval(data);
            }
            return Systematic<TH1>(new_name, container, fType);
        }
    }


    /////////////////////////////////////////////////////////////////////////
    template<class T>
    const std::shared_ptr<T>
    Systematic<T>::
    Up() const {
        if (fType == kMultiverse) {
            throw exceptions::SystematicTypeError(__PRETTY_FUNCTION__, kOneOrTwoSided, fType);
        } else return this->fContainer[0];
    }

    /////////////////////////////////////////////////////////////////////////
    template<class T>
    const std::shared_ptr<T>
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
        std::vector<std::shared_ptr<T>> container(nshifts);
        for (auto ishift = 0; ishift < nshifts; ishift++) {
            container[ishift] = std::move(load(container_dir, std::to_string(ishift)));
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

