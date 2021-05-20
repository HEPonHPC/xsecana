#include "XSecAna/SimpleSignalEstimator.h"

namespace xsec {
  //////////////////////////////////////////////////////////
  template<class HistType>
  const HistType *
  SimpleSignalEstimator<HistType>::
  Background(const HistType * data) const
  {
    return fBackground;
  }

  //////////////////////////////////////////////////////////
  template<class HistType>
  const HistType *
  SimpleSignalEstimator<HistType>::
  Signal(const HistType * data) const
  {
    return data - fBackground;
  }
  
  //////////////////////////////////////////////////////////
  template<class HistType>
  void
  SimpleSignalEstimator<HistType>::
  SaveTo(TDirectory * dir, const std::string & subdir) const
  {
    TDirectory * tmp = gDirectory;
    dir = dir->GetDirectory(subdir.c_str());
    dir->cd();

    TObjString("SimpleSignalEstimator").Write("type");
    fBackground->SaveTo(dir, "fBackground");

    delete dir;
    tmp->cd();
  }
  
  //////////////////////////////////////////////////////////
  template<class HistType>
  std::unique_ptr<SimpleSignalEstimator<HistType> > 
  SimpleSignalEstimator<HistType>::
  LoadFrom(TDirectory * dir, const std::string & subdir)
  {
    dir = dir->GetDirectory(subdir.c_str());

    // make sure we're loading the right type
    TObjString * ptag = (TObjString*) dir->Get("type");
    assert(ptag->GetString() == "SimpleSignalEstimator" && "Type does not match SimpleSignalEstimator");
    delete ptag;
    
    HistType * background   = HistType::LoadFrom(dir, "fBackground");
    return std::make_unique<SimpleSignalEstimator<HistType> >(background);
  }
  
  
  
  
}
