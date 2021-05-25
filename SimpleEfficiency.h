#pragma once

#include "TDirectory.h"

#include "XSecAna/Hist.h"
#include "XSecAna/IEfficiency.h"

namespace xsec {

  template<class HistType = HistXd>
  class SimpleEfficiency : public IEfficiency<HistType> {
  public:
    HistType * Efficiency() override;
    void SaveTo(TDirectory * dir, std::string subdir) const override;
    
    static std::unique_ptr<SimpleEfficiency> LoadFrom(TDirectory * dir, std::string name);

  private:
    HistType * fNumerator;
    HistType * fDenominator;
    HistType * fRatio;
  };

  //////////////////////////////////////////////////////////
  template<class HistType>
  HistType * 
  SimpleEfficiency<HistType>::
  Efficiency()
  {
    if(!fRatio) {
      fRatio = new HistType(*fNumerator);
      *fRatio /= *fDenominator;
    }
    return fRatio;
  }

  //////////////////////////////////////////////////////////
  template<class HistType>
  void
  SimpleEfficiency<HistType>::
  SaveTo(TDirectory * dir, std::string subdir) const
  {
    TDirectory * tmp = gDirectory;
    dir = dir->GetDirectory(subdir.c_str());
    dir->cd();

    TObjString("SimpleEfficiency").Write("type");
    fNumerator  ->SaveTo(dir, "fNumerator"  );
    fDenominator->SaveTo(dir, "fDenominator");

    delete dir;
    tmp->cd();
  }

  //////////////////////////////////////////////////////////
  template<class HistType>
  std::unique_ptr<SimpleEfficiency<HistType> > 
  SimpleEfficiency<HistType>::
  LoadFrom(TDirectory * dir, std::string subdir)
  {
    dir = dir->GetDirectory(subdir.c_str());

    // make sure we're loading the right type
    TObjString * ptag = (TObjString*) dir->Get("type");
    assert(ptag->GetString() == "SimpleEfficiency" && "Type does not match SimpleEfficiency");
    delete ptag;
    
    HistType * numerator   = HistType::LoadFrom(dir, "fNumerator"  );
    HistType * denominator = HistType::LoadFrom(dir, "fDenominator");
    return std::make_unique<SimpleEfficiency<HistType> >(numerator,
							 denominator);
  }
}
