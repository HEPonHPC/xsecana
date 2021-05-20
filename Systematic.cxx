#include "XSecAna/Systematic.h"

#include "TDirectory.h"
#include "TObjString.h"
#include "TString.h"

namespace xsec {

  ////////////////////////////////////////////////////////////
  template<class T>
  template<class F, class... Args>
  OneSidedSystematic<T> *
  OneSidedSystematic<T>::Invoke(F f, Args... args) const
  {
    OneSidedSystematic<T> * ret;
    ret->fName = this->fName;
    ret->fShift = new T(*fShift);

    std::invoke(f, ret->fShift, args...);
    return ret;
  }

  ////////////////////////////////////////////////////////////
  template<class T>
  template<class F, class... Args>
  TwoSidedSystematic<T> *
  TwoSidedSystematic<T>::Invoke(F f, Args... args) const
  {
    return new TwoSidedSystematic(this->fName,
				  fUp->Invoke(f, args...),
				  fDown->Invoke(f, args...));

  }

  ////////////////////////////////////////////////////////////
  template<class T>
  template<class F, class... Args>
  MultiverseSystematic<T> *
  MultiverseSystematic<T>::Invoke(F f, Args... args) const
  {
    MultiverseSystematic<T> * ret = new MultiverseSystematic<T>(this->fName);

    ret->fUniverses = std::vector<T>(this->fUniverses.size());
    for(auto i = 0u; i < fUniverses.size(); i++) {
      ret->fUniverses[i] = new T(*fUniverses[i]);
      std::invoke(&f, ret->fUniverses[i], args...);
    }
    return ret;
  }


  ////////////////////////////////////////////////////////////
  template<class T> 
  void 
  MultiverseSystematic<T>::SaveTo(TDirectory * dir, std::string subdir) const
  {
    TDirectory * tmp = gDirectory;
    dir->mkdir(subdir.c_str());
    dir = dir->GetDirectory(subdir.c_str());
    dir->cd();

    TObjString("MultiverseSystematic").Write("type");
    TObjString(this->fName.c_str()).Write("fName");

    auto mv_dir = dir->mkdir("fUniverses");
    for(auto i = 0u; i < fUniverses.size(); i++) {
      fUniverses[i]->SaveTo(mv_dir, std::to_string(i));
    }

    delete dir;
    tmp->cd();
  }


  ////////////////////////////////////////////////////////////
  template<class T> 
  void 
  OneSidedSystematic<T>::SaveTo(TDirectory * dir, std::string subdir) const
  {
    TDirectory * tmp = gDirectory;
    dir->mkdir(subdir.c_str());
    dir = dir->GetDirectory(subdir.c_str());
    dir->cd();

    TObjString("OneSidedSystematic").Write("type");
    TObjString(this->fName.c_str()).Write("fName");

    fShift->SaveTo(dir, "fShift");


    delete dir;
    tmp->cd();
  }

  ////////////////////////////////////////////////////////////
  template<class T> 
  void 
  TwoSidedSystematic<T>::SaveTo(TDirectory * dir, std::string subdir) const
  {
    TDirectory * tmp = gDirectory;
    dir->mkdir(subdir.c_str());
    dir = dir->GetDirectory(subdir.c_str());
    dir->cd();

    TObjString("TwoSidedSystematic").Write("type");
    TObjString(this->fName.c_str()).Write("fName");
    fUp->SaveTo(dir, "fUp");
    fDown->SaveTo(dir, "fDown");

    delete dir;
    tmp->cd();
  }

  ////////////////////////////////////////////////////////////
  template<class T> 
  std::unique_ptr<MultiverseSystematic<T> >
  MultiverseSystematic<T>::LoadFrom(TDirectory * dir, std::string subdir)
  {
    dir = dir->GetDirectory(subdir.c_str());

    // make sure we're loading the right type
    TObjString * ptag = (TObjString*) dir->Get("type");
    assert(ptag->GetString() == "MultiverseSystematic" && "Type does not match MultiverseSystematic");
    delete ptag;
    
    std::string name = ((TObjString*) dir->Get("fName"))->GetString().Data();
    std::vector<T*> universes;
    auto mv_dir = dir->GetDirectory("fUniverses");
    for(auto universe : *mv_dir->GetListOfKeys()) {
      universes.push_back(T::LoadFrom(dir, universe->GetName()).release());
    }
    return std::make_unique<MultiverseSystematic<T> >(name, universes);
  }

  ////////////////////////////////////////////////////////////
  template<class T> 
  std::unique_ptr<TwoSidedSystematic<T> >
  TwoSidedSystematic<T>::LoadFrom(TDirectory * dir, std::string subdir)
  {
    dir = dir->GetDirectory(subdir.c_str());

    // make sure we're loading the right type
    TObjString * ptag = (TObjString*) dir->Get("type");
    assert(ptag->GetString() == "TwoSidedSystematic" && "Type does not match TwoSidedSystematic");
    delete ptag;
    
    std::string name = ((TObjString*) dir->Get("fName"))->GetString().Data();
    auto up = OneSidedSystematic<T>::LoadFrom(dir, "fUp").release();
    auto dw = OneSidedSystematic<T>::LoadFrom(dir, "fDown").release();
    
    return std::make_unique<TwoSidedSystematic<T> >(name, up, dw);
  }

  ////////////////////////////////////////////////////////////
  template<class T> 
  std::unique_ptr<OneSidedSystematic<T> >
  OneSidedSystematic<T>::LoadFrom(TDirectory * dir, std::string subdir)
  {
    dir = dir->GetDirectory(subdir.c_str());

    // make sure we're loading the right type
    TObjString * ptag = (TObjString*) dir->Get("type");
    assert(ptag->GetString() == "OneSidedSystematic" && "Type does not match OneSidedSystematic");
    delete ptag;
    
    std::string name = ((TObjString*) dir->Get("fName"))->GetString().Data();
    T * shift = T::LoadFrom(dir, "fShift").release();

    return std::make_unique<OneSidedSystematic<T> >(name, shift);
  }

  /*
  ////////////////////////////////////////////////////////////
  template<class T>
  T *
  MultiverseSystematic<T>::Shift(const T * nominal, double nsigma) const
  {  
    MultiverseSystematic<T> hists = this->Invoke(ToHist);
  }

  ///////////////////////////////////////////////////////////////////////
  // Calculate the universe that marks the n sigma shift
  // copied from ana::GenieMultiverse::BinSigma()
  double
  MultiverseSystematic<T>::BinSigma(std::vector<double> events,
				    double nsigma, 
				    double pivot)
  {
    int pivotbin = 0;
    double pivotbincenter = 0;
    std::sort(events.begin(), events.end());
    for(int i = 0; i < (int) events.size() - 1; i++) {
      if(pivot >= events[i] && pivot < events[i+1]) {
	pivotbin = i;
	break;
      }
    }
    pivotbincenter = pivotbin+0.5;
    double count_fraction = 2.0 * (ROOT::Math::normal_cdf(nsigma) - ROOT::Math::normal_cdf(0));

    int nsideevents = 0;
    int lastbinindex = (int) events.size() - 1;
    if(nsigma >= 0) nsideevents = lastbinindex - pivotbin;
    else nsideevents = pivotbin;
    int boundIdx = pivotbincenter + count_fraction*(double)nsideevents;

    int index = 0;
    if(nsigma >= 0) index = std::min(boundIdx, (int)events.size() - 1);
    else index = std::max(boundIdx, 0);
    return events.at(index);
  }
  */
}
