#pragma once

#include <iostream>
#include "XSecAna/Hist.h"
#include "XSecAna/IUnfold.h"

//----------------------------------------------------------------------
#define TEST_HIST(test_name,HIST, target_contents, target_edges, precision) \
      test = (HIST.Contents() - target_contents).isZero(precision);	\
      if(!test || verbose) {						\
	std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << (test? ": PASSED" : ": FAILED") << std::endl; \
	std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << "\t" << HIST.Contents() << std::endl; \
	std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << "\t" << target_contents << std::endl; \
	pass &= test;							\
      }									\
      test = (HIST.Edges() - target_edges).isZero(precision);		\
      if(!test || verbose) {						\
	std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << (test? ": PASSED" : ": FAILED") << std::endl; \
	std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << "\t" << HIST.Edges() << std::endl; \
	std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << "\t" << target_edges << std::endl; \
	pass &= test;							\
      }								
  
//----------------------------------------------------------------------
#define TEST_ARRAY(test_name, arr1, arr2, precision)			\
      test = (arr1 - arr2).isZero(precision);				\
      if(!test || verbose) {						\
	std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << (test? ": PASSED" : ": FAILED") << std::endl; \
	std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << "\t" << arr1 << std::endl; \
	std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << "\t" << arr2 << std::endl; \
	pass &= test;							\
      }									

//----------------------------------------------------------------------
#define TEST_SYSTEMATIC(test_name, syst, up, down)			\
  test = syst.Up() == (up);						\
  if(!test || verbose) {						\
    std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << (test? ": PASSED" : ": FAILED") << std::endl; \
    std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << "\t" << syst.Up().Contents() << std::endl; \
    std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << "\t" << (up).Contents() << std::endl; \
    pass &= test;							\
  }									\
  test = syst.Down() == (down);						\
  if(!test || verbose) {						\
    std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << (test? ": PASSED" : ": FAILED") << std::endl; \
    std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << "\t" << syst.Down().Contents() << std::endl; \
    std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << "\t" << (down).Contents() << std::endl; \
    pass &= test;							\
  }									

//----------------------------------------------------------------------
#define TEST_MULTIVERSE(test_name, mv1, mv2)				\
  test = true;								\
  for(auto imv = 0u; imv < (mv1).GetShifts().size(); imv++) {		\
    test &= (mv1).GetShifts()[imv] == (mv2).GetShifts()[imv];		\
  }									\
  if(!test || verbose) {						\
    std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << (test? ": PASSED" : ": FAILED") << std::endl; \
    for(auto imv = 0u; imv < (mv1).GetShifts().size(); imv++) {	\
      std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << "[" << imv << "]\t" << (mv1).GetShifts()[imv].Contents() << std::endl; \
      std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << "[" << imv << "]\t" << (mv2).GetShifts()[imv].Contents() << std::endl; \
    }									\
    pass &= test;							\
  }



namespace xsec {
  namespace test {
    namespace utils {
      /////////////////////////////////////////////////////////
      template<class Scalar, 
	       int Cols>
      struct Ratio {
	Ratio() {}

	Ratio(Hist<Scalar, Cols> num,
	      Hist<Scalar, Cols> den)
	  : numerator(num), denominator(den)
	{}

	Hist<Scalar, Cols> numerator;
	Hist<Scalar, Cols> denominator;
  
	Hist<Scalar, Cols> ToHist() const
	{
	  return numerator / denominator;
	}
      };

      /////////////////////////////////////////////////////////
      template<class Scalar,
	       int Cols>
      class DummyUnfold : IUnfold<Hist<Scalar, Cols> > {
      public:
	DummyUnfold(int nbins, double scale = 1)
	{
	  fMat = Eigen::Matrix<Scalar, Cols, Cols>::Identity(nbins, nbins) * scale;
	}
  
	Hist<Scalar, Cols> Truth(const Hist<Scalar, Cols> & reco) const
	{ return Hist<Scalar, Cols>(fMat * reco.Contents().matrix().transpose(), reco.Edges()); }
  
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
	  dir = dir->GetDirectory(subdir.c_str());
	  dir->cd();

	  auto mat = *Hist<Scalar, Eigen::Dynamic>::LoadFrom(dir, "fMat");
    
	  tmp->cd();

	  return std::make_unique<DummyUnfold<Scalar, Cols > >(Eigen::Map<
							       const Eigen::Matrix<Scalar, Cols, Cols>
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

      /////////////////////////////////////////////////////////
      template<typename ArrayType>
      ArrayType is_contiguous(const ArrayType & array)
      { 
	ArrayType a2 = array;
	auto const data = array.data();
	for (auto i = 0; i < array.size(); ++i) {
	  a2(i) = *(&data[0] + i);
	}        
	return a2;
      }

      /////////////////////////////////////////////////////////
      std::string test_dir() 
      {
	return std::string(std::getenv("SRT_PRIVATE_CONTEXT")) + "/XSecAna/test/";
      }

    } // utils
  } // test
} // xsec
