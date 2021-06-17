
#include <iostream>
#include "XSecAna/CAFAna/CAFAnaUnfold.h"
#include "CAFAna/Unfold/UnfoldIterative.h"

#include "XSecAna/CrossSectionAnalysis.h"
#include "XSecAna/ICrossSection.h"
#include "XSecAna/SimpleQuadSum.h"
#include "XSecAna/Systematic.h"
#include "XSecAna/Hist.h"
#include "XSecAna/SimpleEfficiency.h"
#include "XSecAna/SimpleFlux.h"
#include "XSecAna/SimpleSignalEstimator.h"

#include "XSecAna/test/test_utils.h"

using namespace xsec;

int main(int argc, char ** argv)
{
  typedef ICrossSection<SimpleSignalEstimator<HistXd>,
			test::utils::DummyUnfold<double, -1>,
			SimpleEfficiency<HistXd>,
			SimpleFlux<HistXd> > SimpleCrossSection;

  HistXd hist;
  SimpleFlux<HistXd> flux;
  CrossSectionAnalysis<SimpleCrossSection,
		       SimpleQuadSum<HistXd> > analysis;

  return 1;
}
