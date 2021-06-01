#pragma once

#include "XSecAna/Hist.h"

namespace xsec {
  template<class HistType = HistXd>
  class IUnfold {
  public:
    virtual HistType Truth(const HistType & reco) const = 0; 
  private:
    
  };
  
  template<class Scalar,
	   int Cols>
  class DummyUnfold : IUnfold<Hist<Scalar, Cols> > {
  public:
    DummyUnfold(int nbins)
    {
      fMat = Eigen::Matrix<Scalar, Cols, Cols>::Identity(nbins, nbins);
    }
  
    Hist<Scalar, Cols> Truth(const Hist<Scalar, Cols> & reco) const
    { return reco*1; }
  
    void SaveTo(TDirectory * dir, std::string subdir)
    {
      TDirectory * tmp = gDirectory;
      dir = dir->mkdir(subdir.c_str());
      dir->cd();
    
      auto bins = Eigen::Array<Scalar, 1, Eigen::Dynamic>::LinSpaced(std::pow(fMat.cols(),2)+1, 
								     0,
								     std::pow(fMat.cols(),2)+1);
      Hist<Scalar, Eigen::Dynamic>(Eigen::Map<Eigen::Array<Scalar, 1, Eigen::Dynamic> >(fMat.data(),
									      1, 
									      std::pow(fMat.cols(),2)),
				   bins).SaveTo(dir, "fMat");
    
      tmp->cd();
    }

    static std::unique_ptr<DummyUnfold<Scalar, Cols> > LoadFrom(TDirectory * dir, 
								std::string subdir)
    {
      TDirectory * tmp = gDirectory;
      dir = dir->mkdir(subdir.c_str());
      dir->cd();

      auto mat = *Hist<Scalar, Eigen::Dynamic>::LoadFrom(dir, "fMat");
    
      tmp->cd();
      return std::make_unique<DummyUnfold<Scalar, Cols > >(Eigen::Map<
							   Eigen::Matrix<Scalar, Cols, Cols>
							   >
							   (mat.Contents().data(),
							    std::sqrt(mat.Contents().size()),
							    std::sqrt(mat.Contents().size())));
    
    }
    DummyUnfold(const Eigen::Matrix<Scalar, Cols, Cols> & mat)
      : fMat(mat)
    {}

  private:
    Eigen::Matrix<Scalar, Cols, Cols> fMat;
  };


}
