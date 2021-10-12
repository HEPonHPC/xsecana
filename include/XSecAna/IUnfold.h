#pragma once

#include "XSecAna/IMeasurement.h"
#include "TParameter.h"

namespace xsec {
    class IUnfold {
    public:
        virtual const TH1 * Truth(const TH1 * reco) const = 0;
        static std::unique_ptr<IUnfold> LoadFrom(xsec::type::LoadFunction<IUnfold> load,
                                                 TDirectory * dir,
                                                 const std::string & name) {
            return load(dir, name);
        }

        virtual ~IUnfold() = default;
    };

    class IUnfoldEstimator : public IUnfold, public IMeasurement{};
    class IEigenUnfoldEstimator : public IUnfold, public virtual IMeasurement, public IEigenEval{};

    /////////////////////////////////////////////////////////
    class IdentityUnfolder : public IEigenUnfoldEstimator {
    public:

        explicit IdentityUnfolder(int nbins_and_uof);

        const TH1 * Truth(const TH1 * reco) const override;

        void SaveTo(TDirectory * dir, const std::string & subdir) const override;

        static std::unique_ptr<IMeasurement> LoadFrom(TDirectory * dir, const std::string & subdir);

    protected:
        void _eval_impl(const Array & data, const Array & error, ArrayRef result, ArrayRef rerror) const override;
        Eigen::MatrixXd fMat;
    };


    void
    IdentityUnfolder::
    _eval_impl(const Array & data, const Array & error, ArrayRef result, ArrayRef rerror) const {
        result = fMat * data.matrix();
        rerror = error;
    }

    /////////////////////////////////////////////////////////
    IdentityUnfolder::
    IdentityUnfolder(int nbins_and_uof) {
        fMat = Eigen::MatrixXd::Identity(nbins_and_uof, nbins_and_uof);
    }

    /////////////////////////////////////////////////////////
    const TH1 *
    IdentityUnfolder::
    Truth(const TH1 * reco) const {
        return this->Eval(reco);
    }

    /////////////////////////////////////////////////////////
    void
    IdentityUnfolder::
    SaveTo(TDirectory * dir, const std::string & subdir) const {
        TDirectory * tmp = gDirectory;
        dir = dir->mkdir(subdir.c_str());
        dir->cd();

        TParameter<double>("cols", fMat.cols()).Write("cols");

        tmp->cd();
    }

    /////////////////////////////////////////////////////////
    std::unique_ptr<IMeasurement>
    IdentityUnfolder::
    LoadFrom(TDirectory* dir,
             const std::string & subdir) {
        TDirectory * tmp = gDirectory;
        dir = dir->GetDirectory(subdir.c_str());
        dir->cd();

        auto cols = ((TParameter<double> *) dir->Get("cols"))->GetVal();
        tmp->cd();

        return std::make_unique<IdentityUnfolder>(cols);
    }

}
