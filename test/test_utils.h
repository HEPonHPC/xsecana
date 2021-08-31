#pragma once

#include <iostream>
#include "XSecAna/Hist.h"
#include "XSecAna/IUnfold.h"
#include "XSecAna/SimpleSignalEstimator.h"
#include "XSecAna/SimpleEfficiency.h"
#include "XSecAna/SimpleFlux.h"
#include "XSecAna/ICrossSection.h"


//----------------------------------------------------------------------
#define TEST_HIST(test_name,HIST, target_contents, target_edges, precision) \
  test = ((HIST).Contents() - (target_contents)).isZero(precision);	\
      if(!test || verbose) {						\
	std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << (test? ": PASSED" : ": FAILED") << std::endl; \
	std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << "\t" << (HIST).Contents() << std::endl; \
	std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << "\t" << (target_contents) << std::endl; \
	pass &= test;							\
      }									\
      test = ((HIST).Edges() - (target_edges)).isZero(precision);	\
      if(!test || verbose) {						\
	std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << (test? ": PASSED" : ": FAILED") << std::endl; \
	std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << "\t" << (HIST).Edges() << std::endl; \
	std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << "\t" << (target_edges) << std::endl; \
	pass &= test;							\
      }

//----------------------------------------------------------------------
#define TEST_ARRAY(test_name, arr1, arr2, precision)			\
  test = ((arr1) - (arr2)).isZero(precision);				\
      if(!test || verbose) {						\
	std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << (test? ": PASSED" : ": FAILED") << std::endl; \
	std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << "\t" << (arr1) << std::endl; \
	std::cerr << __PRETTY_FUNCTION__ << "\t" << test_name << "\t" << (arr2) << std::endl; \
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
	{
	  return Hist<Scalar, Cols>(fMat * reco.Contents().matrix().transpose(),
				    reco.Edges(),
				    reco.Exposure());
	}

	void SaveTo(TDirectory * dir, const std::string& subdir) const override
	{
	  TDirectory * tmp = gDirectory;
	  dir = dir->mkdir(subdir.c_str());
	  dir->cd();

	  auto bins = Eigen::Array<Scalar, 1, Eigen::Dynamic>::LinSpaced(std::pow(fMat.cols(),2)+1,
									 0,
									 std::pow(fMat.cols(),2)+1);
	  Hist<Scalar, Eigen::Dynamic>(Eigen::Map<const Eigen::Array<Scalar, 1, Eigen::Dynamic> >(fMat.data(),
												  1,
												  std::pow(fMat.cols(),2)),
				       bins).SaveTo(dir, "fMat");

	  tmp->cd();
	}

	static std::unique_ptr<DummyUnfold<Scalar, Cols> > LoadFrom(TDirectory * dir,
								    const std::string& subdir)
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
      typedef Hist<double, 10> histtype;
      typedef ICrossSection<SimpleSignalEstimator<histtype>,
			    test::utils::DummyUnfold<double, 10>,
			    SimpleEfficiency<histtype>,
			    SimpleFlux<histtype> > SimpleCrossSection;

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
	return "./";
      }

      const static double data_exposure = 0.5;
      /////////////////////////////////////////////////////////
      template<class Scalar, int Cols>
      Hist<Scalar, Cols> get_simple_data()
      {
	return Hist<Scalar, Cols>(Eigen::Array<Scalar, 1, Cols>::LinSpaced(10, 1, 2) +  2,
				  Eigen::Array<Scalar, 1, EdgesSize(Cols)>::LinSpaced(11, 0, 20),
				  data_exposure); // data will often come at different exposure from mc
      }

      /////////////////////////////////////////////////////////
      template<class Scalar, int Cols>
      Hist<Scalar, Cols> get_simple_background()
      {
	return Hist<Scalar, Cols>(Eigen::Array<Scalar, 1, Cols>::Ones(),
				  Eigen::Array<Scalar, 1, EdgesSize(Cols)>::LinSpaced(11, 0, 20));
      }

      /////////////////////////////////////////////////////////
      template<class Scalar, int Cols>
      Hist<Scalar, Cols> get_simple_signal()
      {
	return get_simple_data<Scalar, Cols>() - get_simple_background<Scalar, Cols>();
      }

      /////////////////////////////////////////////////////////
      template<class Scalar, int Cols>
      Hist<Scalar, Cols> get_simple_nominal_hist()
      {
	return Hist<Scalar, Cols>(Eigen::Array<Scalar, 1, Cols>::LinSpaced(10, -0.5, 0.5) +  4,
				  Eigen::Array<Scalar, 1, EdgesSize(Cols)>::LinSpaced(11, 0, 20));
      }

      /////////////////////////////////////////////////////////
      template<class Scalar, int Cols>
      Hist<Scalar, Cols> get_simple_up_hist()
      {
	auto nominal = get_simple_nominal_hist<Scalar, Cols>();
	auto step = Eigen::Array<Scalar, 1, Cols>::LinSpaced(10, 0, .3).reverse();

	return Hist<Scalar, Cols>(nominal.Contents() + nominal.Contents() * step.pow(2),
				  nominal.Edges());
      }

      /////////////////////////////////////////////////////////
      template<class Scalar, int Cols>
      Hist<Scalar, Cols> get_simple_down_hist()
      {
	auto nominal = get_simple_nominal_hist<Scalar, Cols>();
	auto step = Eigen::Array<Scalar, 1, Cols>::LinSpaced(10, 0, .3);

	return Hist<Scalar, Cols>(nominal.Contents() - nominal.Contents() * step.pow(2),
				  nominal.Edges());
      }

      /////////////////////////////////////////////////////////
      // make a cross section object that evaluates out to
      //  - the input array when folded
      //  - 2 times the input array when unfolded
      template<class Scalar, int Cols>
      SimpleCrossSection make_simple_xsec(Hist<Scalar, Cols> val)
      {
	Hist<Scalar, Cols> ones(Eigen::Array<Scalar, 1, Cols>::Ones(),
				val.Edges());

	auto efficiency = new SimpleEfficiency<Hist<Scalar, Cols> >(get_simple_signal<Scalar, Cols>(),
								    val);
	auto flux = new SimpleFlux(Hist<Scalar, Cols>(ones.Contents(),
						      ones.Edges(),
						      get_simple_data<Scalar, Cols>().Exposure()));
	auto signal_estimator = new SimpleSignalEstimator(get_simple_background<Scalar, Cols>());
	auto unfold = new test::utils::DummyUnfold<Scalar, Cols>(ones.Contents().size(), 2);  // = 2
	return SimpleCrossSection(efficiency,
				  signal_estimator,
				  flux,
				  unfold);
      }

      const static double ntargets = 1e4;

      /////////////////////////////////////////////////////////
      template<class HistType>
      std::vector<HistType> make_simple_hist_multiverse(const HistType & hnominal, int nuniverses)
      {
	double maxy =  0.1;
	double miny = -0.1;
	double step = (maxy - miny) / (nuniverses-1);
	std::vector<HistType> hist_universes(nuniverses);
	for(auto i = 0; i < nuniverses; i++) {
	  hist_universes[i] = hnominal + hnominal * (miny + step * i);
	}
	return hist_universes;
      }

      /////////////////////////////////////////////////////////
      template<class HistType>
      std::vector<SimpleCrossSection> make_simple_xsec_multiverse(const HistType & hnominal, int nuniverses)
      {
	std::vector<SimpleCrossSection> xsec_universes(nuniverses);
	auto hist_universes = make_simple_hist_multiverse(hnominal, nuniverses);

	for(auto i = 0; i < nuniverses; i++) {
	  xsec_universes[i] = make_simple_xsec(hist_universes[i]);
	}
	return xsec_universes;
      }



    } // utils
  } // test
} // xsec
