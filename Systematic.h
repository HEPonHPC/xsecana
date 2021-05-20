#pragma once

#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/XSec/Multiverse.h"

namespace xsec {
  //-------------------------------------------------------
  template<class T>
  class Systematic {
  public:
    Systematic(std::string name) 
      : fName(name) 
    {}
    std::string GetName() const { return fName; }

  protected:
    std::string fName;
  };

  //-------------------------------------------------------
  template<class T>
  class OneSidedSystematic : public Systematic<T> {
  public:
    OneSidedSystematic(std::string name,
		       T * shift)
      : Systematic<T>(name), fShift(shift)
    {}
    
    const T * GetShift() const { return fShift; }

    template<class F, class... Args>
    OneSidedSystematic<T> * Invoke(F f, Args... args) const;

    void SaveTo(TDirectory * dir, std::string name) const; 
    static std::unique_ptr<OneSidedSystematic> LoadFrom(TDirectory * dir, std::string name); 
  private:
    T * fShift;
    
  };

  //-------------------------------------------------------
  template<class T>
  class TwoSidedSystematic : public Systematic<T> {
  public:
    TwoSidedSystematic(std::string name,
		       T * up,
		       T * dw)
      : Systematic<T>(name),
	fUp(new OneSidedSystematic<T>(name+"_up", up)),
	fDown(new OneSidedSystematic<T>(name+"_dw", dw))
    {}

    void SaveTo(TDirectory * dir, std::string name) const; 
    static std::unique_ptr<TwoSidedSystematic> LoadFrom(TDirectory * dir, std::string name); 

    template<class F, class... Args>
    TwoSidedSystematic<T> * Invoke(F f, Args... args) const;

    // LoadFrom constructor
    TwoSidedSystematic(std::string name,
		       OneSidedSystematic<T> * up,
		       OneSidedSystematic<T> * dw)
      : Systematic<T>(name),
	fUp(up),
	fDown(dw)
    {}
  private:
    OneSidedSystematic<T> * fUp;
    OneSidedSystematic<T> * fDown;
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
    MultiverseSystematic<T> * Invoke(F f, Args... args) const;

    void SaveTo(TDirectory * dir, std::string name) const; 
    static std::unique_ptr<MultiverseSystematic> LoadFrom(TDirectory * dir, std::string name); 

  private:
    std::vector<T> fUniverses;

  };

}
