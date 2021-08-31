#include <iostream>

#include "XSecAna/Analysis.h"
#include "XSecAna/ICrossSection.h"
#include "XSecAna/SimpleQuadSum.h"
#include "XSecAna/Systematic.h"
#include "XSecAna/Hist.h"
#include "XSecAna/SimpleEfficiency.h"
#include "XSecAna/SimpleFlux.h"
#include "XSecAna/SimpleSignalEstimator.h"

#include "test_utils.h"

using namespace xsec;

int main(int argc, char ** argv)
{
  typedef ICrossSection<HistXd,
			SimpleSignalEstimator<HistXd>,
			test::utils::DummyUnfold<double, -1>,
			SimpleEfficiency<HistXd>,
			SimpleFlux<HistXd> > SimpleCrossSection;

  HistXd hist;
  SimpleFlux<HistXd> flux;
  Analysis<SimpleCrossSection,
	   SimpleQuadSum<SimpleCrossSection> > analysis;

  return 0;
}
