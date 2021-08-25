#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/ReweightableSpectrum.h"

#include "XSecAna/Hist.h"
namespace xsec::cafana {
  template<class HistType>
  ana::ReweightableSpectrum ToReweightableSpectrum(const HistType & hist);
  
  template<class HistType>
  ana::Spectrum ToSpectrum(const HistType & hist);
  
  template<class HistType>
  HistType ToHist(const ana::Spectrum & spec);

  template<class HistType>
  HistType ToHist(const ana::ReweightableSpectrum & spec);
}
