#pragma once

#include "TDirectory.h"

#include <functional>
#include "XSecAna/Type.h"
namespace xsec {

  //-------------------------------------------------------
  template<class T>
  class Systematic {
  public:
    Systematic(std::string name) 
      : fName(name)
    {}

    Systematic(std::string name,
	       T up,
	       T dw)
      : fName(name),
	fUp(new T(up)),
	fDown(new T(dw))
    {}
    Systematic(std::string name,
	       T shift)
      : fName(name),
	fUp(new T(shift)),
	fDown(NULL)
    {}

    std::pair<const T *, const T *> GetShifts() const
    { return std::pair<const T *, const T *>{fUp, fDown ? fDown : fUp}; }
    
    template<class F, class... Args>
    Systematic<std::invoke_result_t<F, T, Args...> > Invoke(F f, Args... args) const;

    void SaveTo(TDirectory * dir, std::string name) const; 
    static std::unique_ptr<Systematic> LoadFrom(TDirectory * dir, std::string name); 

    const T & Up  () const { return fUp  ; }
    const T & Down() const { return fDown ? fDown : fUp; }

  protected:
    std::string fName;
  private:
    T * fUp;
    T * fDown;
  };
    
  //-------------------------------------------------------
  template<class T>
  class MultiverseSystematic : public Systematic<T> {
  public:
    MultiverseSystematic(std::string name,
			 std::vector<T>  universes)
      : Systematic<T>(name),
	fUniverses(universes)
    {}
      
    template<class F, class... Args>
    MultiverseSystematic<std::invoke_result_t<F, T, Args...> > Invoke(F f, Args... args) const;

    void SaveTo(TDirectory * dir, std::string name) const; 
    
    std::pair<const T *, const T *> GetShifts() const;

    const std::vector<T> GetUniverses() const { return fUniverses; }
    
    static std::unique_ptr<MultiverseSystematic> LoadFrom(TDirectory * dir, std::string name); 

    auto NSigmaShift(double nsigma,
		     const T & nominal) const;
  private:

    template<class Scalar>
    Scalar BinSigma(const double & nsigma, std::vector<Scalar> & universes, const Scalar & nominal) const;
    std::vector<T> fUniverses;

  };

  ////////////////////////////////////////////////////////////
  template<class T>
  template<class F, class... Args>
  Systematic<std::invoke_result_t<F, T, Args...> >
  Systematic<T>::Invoke(F f, Args... args) const
  {
    if(fDown) {
      return Systematic<std::invoke_result_t<F, T, Args...> >(this->fName,
							      std::invoke(f, *fUp  , args...),
							      std::invoke(f, *fDown, args...));
    }
    else {
      return Systematic<std::invoke_result_t<F, T, Args...> >(this->fName,
							      std::invoke(f, *fUp  , args...));
    }

  }

  ////////////////////////////////////////////////////////////
  template<class T>
  template<class F, class... Args>
  MultiverseSystematic<std::invoke_result_t<F, T, Args...> > 
  MultiverseSystematic<T>::Invoke(F f, Args... args) const
  {
    std::vector<std::invoke_result_t<F, T, Args...> > universes(this->fUniverses.size());
    for(auto i = 0u; i < fUniverses.size(); i++) {
      universes[i] = std::invoke(f, this->fUniverses[i], args...);
    }
    return MultiverseSystematic<std::invoke_result_t<F, T, Args...> >(this->fName,
								      universes);
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
    TObjString(std::to_string(this->fUniverses.size()).c_str()).Write("NUniverses");
    
    
    auto mv_dir = dir->mkdir("fUniverses");
    for(auto i = 0u; i < fUniverses.size(); i++) {
      fUniverses[i].SaveTo(mv_dir, std::to_string(i));
    }

    tmp->cd();
  }

  ////////////////////////////////////////////////////////////
  template<class T> 
  void 
  Systematic<T>::SaveTo(TDirectory * dir, std::string subdir) const
  {
    TDirectory * tmp = gDirectory;
    dir->mkdir(subdir.c_str());
    dir = dir->GetDirectory(subdir.c_str());
    dir->cd();

    TObjString("Systematic").Write("type");
    TObjString(this->fName.c_str()).Write("fName");
    fUp->SaveTo(dir, "fUp");
    if(fDown) fDown->SaveTo(dir, "fDown");

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
    int nuniverses = std::atoi(((TObjString*) dir->Get("NUniverses"))->GetString().Data());
    std::vector<T> universes(nuniverses);
    auto mv_dir = dir->GetDirectory("fUniverses");
    for(auto imv = 0; imv < nuniverses; imv++) {
      universes[imv] = *T::LoadFrom(mv_dir, std::to_string(imv)).release();
    }
    return std::make_unique<MultiverseSystematic<T> >(name, universes);
  }

  ////////////////////////////////////////////////////////////
  template<class T> 
  std::unique_ptr<Systematic<T> >
  Systematic<T>::LoadFrom(TDirectory * dir, std::string subdir)
  {
    dir = dir->GetDirectory(subdir.c_str());

    // make sure we're loading the right type
    TObjString * ptag = (TObjString*) dir->Get("type");
    assert(ptag->GetString() == "Systematic" && "Type does not match Systematic");
    delete ptag;
    
    std::string name = ((TObjString*) dir->Get("fName"))->GetString().Data();

    if(dir->GetDirectory("fDown")) {
      std::unique_ptr<T> up = T::LoadFrom(dir, "fUp");
      std::unique_ptr<T> dw = T::LoadFrom(dir, "fDown");
      return std::make_unique<Systematic<T> >(name, *up, *dw);
    }
    else {
      std::unique_ptr<T> up = T::LoadFrom(dir, "fUp");
      return std::make_unique<Systematic<T> >(name, *up);
    }
  }

  ////////////////////////////////////////////////////////////
  /*
  template<class T>
  template<template<class, int> Hist,
	   class Scalar,
	   int Cols>  
  Hist<Scalar, Cols>
  MultiverseSystematic<T>::
  NSigmaShift(double sigma,
	      const std::vector<Hist<Scalar, Cols> > & universes,
	      const Hist<Scalar, Cols> & nominal)
  {
    Eigen::Array<Scalar, 1,  Cols> shift(nominal.NBins());

    for(auto ibin = 0u; ibin < nominal.NBins(); ibin++) {
      std::vector<Scalar> vals;
      for(auto iuniv = 0u; iuniv < universes.size(); iuniv++) {
	vals.push_back(universes.Contents()(ibin));
      }
      shift(ibin) = BinSigma(sigma, vals, nominal.Contents(ibin));
    }

    return shift;
  }
  */
  ////////////////////////////////////////////////////////////
  template<class T>
  template<class Scalar>
  Scalar
  MultiverseSystematic<T>::
  BinSigma(const double & nsigma,
	   std::vector<Scalar> & universes,
	   const Scalar & nominal) const
  {
    int pivotbin = 0;
    double pivotbincenter = 0;
    std::sort(universes.begin(), universes.end());
    for(auto i = 0u; i < universes.size() - 1; i++) {
      if(nominal >= universes[i] && nominal < universes[i+1]) {
	pivotbin = i;
	break;
      }
    }
    pivotbincenter = pivotbin+0.5;
    double count_fraction = std::erf(nsigma / std::sqrt(2));

    int nsideevents = 0;
    int lastbinindex = (int) universes.size() - 1;
    if(nsigma >= 0) nsideevents = lastbinindex - pivotbin;
    else nsideevents = pivotbin;
    int boundIdx = pivotbincenter + count_fraction*(double)nsideevents;

    int index = 0;
    if(nsigma >= 0) index = std::min(boundIdx, (int)universes.size() - 1);
    else index = std::max(boundIdx, 0);
    return universes.at(index);
  }

  ////////////////////////////////////////////////////////////
  template<class T>
  auto 
  MultiverseSystematic<T>::
  NSigmaShift(double nsigma, 
	      const T & nominal) const
  {
    if constexpr(type::IsHist<T>()) {
	T shift = nominal;

	for(auto ibin = 0u; ibin < nominal.NBins(); ibin++) {
	  std::vector<decltype(std::declval<T>().GetBinContent(1))> vals;
	  for(auto iuniv = 0u; iuniv < fUniverses.size(); iuniv++) {
	    vals.push_back(fUniverses[iuniv].GetBinContent(ibin));
	  }
	  shift.SetBinContent(ibin, BinSigma(nsigma, vals, nominal.GetBinContent(ibin)));
	}
	return shift;
      }
    else {
      return this->Invoke(&T::ToHist).NSigmaShift(nsigma, nominal.ToHist());
    }
  }

}
