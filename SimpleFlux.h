#pragma once

#include "TDirectory.h"

#include "XSecAna/Hist.h"
#include "XSecAna/IFlux.h"

namespace xsec {
  template<class HistType = HistXd>
  class SimpleFlux : IFlux<HistType> {
  public:
    SimpleFlux() {}

    HistType * Flux() override;
    void SaveTo(TDirectory * dir, std::string subdir) const override;
    static std::unique_ptr<SimpleFlux<HistType> > LoadFrom(TDirectory * dir, std::string subdir);
  private:
    // hold the raw histogram
    HistType * fFlux;
  };

  //////////////////////////////////////////////////////////
  template<class HistType>
  HistType * 
  SimpleFlux<HistType>::
  Flux()
  {
    return fFlux;
  }

  //////////////////////////////////////////////////////////
  template<class HistType>
  void
  SimpleFlux<HistType>::
  SaveTo(TDirectory * dir, std::string subdir) const
  {
    TDirectory * tmp = gDirectory;
    dir = dir->GetDirectory(subdir.c_str());
    dir->cd();

    TObjString("SimpleFlux").Write("type");
    fFlux->SaveTo(dir, "fFlux");

    delete dir;
    tmp->cd();
  }

  //////////////////////////////////////////////////////////
  template<class HistType>
  std::unique_ptr<SimpleFlux<HistType> > 
  SimpleFlux<HistType>::
  LoadFrom(TDirectory * dir, std::string subdir)
  {
    dir = dir->GetDirectory(subdir.c_str());

    // make sure we're loading the right type
    TObjString * ptag = (TObjString*) dir->Get("type");
    assert(ptag->GetString() == "SimpleFlux" && "Type does not match SimpleFlux");
    delete ptag;
    
    HistType * flux   = HistType::LoadFrom(dir, "fFlux");
    return std::make_unique<SimpleFlux<HistType> >(flux);
  }
}
