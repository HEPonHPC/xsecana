#pragma once

#include "XSecAna/Hist.h"
#include "TParameter.h"

namespace xsec {
    template<class HistType = HistXd>
    class IUnfold {
    public:
        virtual HistType Truth(const HistType & reco) const = 0;

        virtual void SaveTo(TDirectory * dir, const std::string & name) const = 0;

        static std::unique_ptr<IUnfold> LoadFrom(xsec::type::LoadFunction<IUnfold> load,
                                                 TDirectory * dir,
                                                 const std::string & name) {
            return load(dir, name);
        }

        virtual ~IUnfold() = default;
    };

    /////////////////////////////////////////////////////////
    template<class Scalar,
            int Cols>
    class IdentityUnfold : public IUnfold<xsec::Hist<Scalar, Cols> > {
    public:

        explicit IdentityUnfold(int nbins);

        Hist <Scalar, Cols> Truth(const Hist <Scalar, Cols> & reco) const override;

        void SaveTo(TDirectory * dir, const std::string & subdir) const override;

        static std::unique_ptr<IUnfold<Hist<Scalar, Cols>>> LoadFrom(TDirectory * dir, const std::string & subdir);

    protected:
        Eigen::Matrix<Scalar,
                      ContentsAndUOFSize(Cols),
                      ContentsAndUOFSize(Cols)> fMat;
    };


/////////////////////////////////////////////////////////
    template<class Scalar,
            int Cols>
    IdentityUnfold<Scalar,
                   Cols>::
    IdentityUnfold(int nbins_and_uof) {
        fMat = Eigen::Matrix<Scalar,
                             ContentsAndUOFSize(Cols),
                             ContentsAndUOFSize(Cols)>::
        Identity(nbins_and_uof, nbins_and_uof);
    }

/////////////////////////////////////////////////////////
    template<class Scalar,
            int Cols>
    Hist <Scalar, Cols>
    IdentityUnfold<Scalar,
                   Cols>::
    Truth(const Hist <Scalar, Cols> & reco) const {
        return Hist<Scalar, Cols>(fMat * reco.ContentsAndUOF().matrix().transpose(),
                                  reco.EdgesAndUOF(),
                                  reco.Exposure());
    }

/////////////////////////////////////////////////////////
    template<class Scalar,
            int Cols>
    void
    IdentityUnfold<Scalar,
                   Cols>::
    SaveTo(TDirectory * dir, const std::string & subdir) const {
        TDirectory * tmp = gDirectory;
        dir = dir->mkdir(subdir.c_str());
        dir->cd();

        TParameter<Scalar>("cols", fMat.cols()).Write("cols");

        tmp->cd();
    }

/////////////////////////////////////////////////////////
    template<class Scalar,
            int Cols>
    std::unique_ptr<IUnfold<Hist<Scalar, Cols>>>
    IdentityUnfold<Scalar,
                   Cols>::
    LoadFrom(TDirectory
             * dir,
             const std::string & subdir) {
        TDirectory * tmp = gDirectory;
        dir = dir->GetDirectory(subdir.c_str());
        dir->cd();

        auto cols = ((TParameter<Scalar> *) dir->Get("cols"))->GetVal();
        tmp->cd();

        return std::make_unique<IdentityUnfold<Scalar, Cols> >(cols);
    }

}
