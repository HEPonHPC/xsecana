#pragma once

#include "XSecAna/Hist.h"
#include "XSecAna/Type.h"
#include "TParameter.h"


namespace xsec {
    class IUnfold {
    public:
        virtual Hist Truth(const Hist & reco) const = 0;

        virtual void SaveTo(TDirectory * dir, const std::string & name) const = 0;

        static std::unique_ptr<IUnfold> LoadFrom(xsec::type::LoadFunction<IUnfold> load,
                                                 TDirectory * dir,
                                                 const std::string & name) {
            return load(dir, name);
        }

        virtual ~IUnfold() = default;
    };

    /////////////////////////////////////////////////////////
    class IdentityUnfold : public IUnfold {
    public:

        explicit IdentityUnfold(int nbins_and_uof);

        Hist Truth(const Hist & reco) const override;

        void SaveTo(TDirectory * dir, const std::string & subdir) const override;

        static std::unique_ptr<IUnfold> LoadFrom(TDirectory * dir, const std::string & subdir);

    protected:
        Eigen::MatrixXd fMat;
    };


    /////////////////////////////////////////////////////////
    IdentityUnfold::
    IdentityUnfold(int nbins_and_uof) {
        fMat = Eigen::MatrixXd::Identity(nbins_and_uof, nbins_and_uof);
    }

    /////////////////////////////////////////////////////////
    Hist
    IdentityUnfold::
    Truth(const Hist & reco) const {
        return Hist(fMat * reco.ContentsAndUOF().matrix().transpose(),
                    reco.EdgesAndUOF(),
                    reco.Exposure());
    }

    /////////////////////////////////////////////////////////
    void
    IdentityUnfold::
    SaveTo(TDirectory * dir, const std::string & subdir) const {
        TDirectory * tmp = gDirectory;
        dir = dir->mkdir(subdir.c_str());
        dir->cd();

        TParameter<double>("cols", fMat.cols()).Write("cols");

        tmp->cd();
    }

    /////////////////////////////////////////////////////////
    std::unique_ptr<IUnfold>
    IdentityUnfold::
    LoadFrom(TDirectory
             * dir,
             const std::string & subdir) {
        TDirectory * tmp = gDirectory;
        dir = dir->GetDirectory(subdir.c_str());
        dir->cd();

        auto cols = ((TParameter<double> *) dir->Get("cols"))->GetVal();
        tmp->cd();

        return std::make_unique<IdentityUnfold>(cols);
    }

}
