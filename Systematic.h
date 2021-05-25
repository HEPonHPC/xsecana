#pragma once

#include "TDirectory.h"

#include <functional>
namespace xsec {

  //-------------------------------------------------------
  template<class T>
  class Systematic {
  public:
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
			 std::vector<T*>  universes)
      : Systematic<T>(name),
	fUniverses(universes)
    {}
      
    template<class F, class... Args>
    Systematic<std::invoke_result_t<F, T, Args...> > Invoke(F f, Args... args) const;

    void SaveTo(TDirectory * dir, std::string name) const; 
    
    std::pair<const T *, const T *> GetShifts() const;

    const std::vector<T> GetUniverses() const { return fUniverses; }
    
    static std::unique_ptr<MultiverseSystematic> LoadFrom(TDirectory * dir, std::string name); 

  private:
    double BinSigma(std::vector<double> values, double nsigma, double pivot);
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
  Systematic<std::invoke_result_t<F, T, Args...> > 
  MultiverseSystematic<T>::Invoke(F f, Args... args) const
  {
    std::vector<std::invoke_result_t<F, T, Args...> > universes(this->fUniverses.size());
    for(auto i = 0u; i < fUniverses.size(); i++) {
      universes[i] = std::invoke(&f, *this->fUniverses[i], args...);
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
  Systematic<T>::SaveTo(TDirectory * dir, std::string subdir) const
  {
    TDirectory * tmp = gDirectory;
    dir->mkdir(subdir.c_str());
    dir = dir->GetDirectory(subdir.c_str());
    dir->cd();

    TObjString("Systematic").Write("type");
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
  std::unique_ptr<Systematic<T> >
  Systematic<T>::LoadFrom(TDirectory * dir, std::string subdir)
  {
    dir = dir->GetDirectory(subdir.c_str());

    // make sure we're loading the right type
    TObjString * ptag = (TObjString*) dir->Get("type");
    assert(ptag->GetString() == "Systematic" && "Type does not match Systematic");
    delete ptag;
    
    std::string name = ((TObjString*) dir->Get("fName"))->GetString().Data();
    auto up = T::LoadFrom(dir, "fUp").release();
    auto dw = T::LoadFrom(dir, "fDown").release();
    
    return std::make_unique<Systematic<T> >(name, up->GetShift(), dw->GetShift());
  }
}
