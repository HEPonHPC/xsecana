#include "XSecAna/Hist.h"

namespace xsec {
  /////////////////////////////////////////////////////////
  template<typename Scalar, 
	   int Rows, 
	   int Cols>
  Hist<Scalar, Rows, Cols>
  Hist<Scalar, Rows, Cols>::
  operator-(const Hist & rhs) const
  {
    Hist<Scalar, Rows, Cols> ret = *this; // copy this
    ret -= rhs;
    return ret;
  }

  /////////////////////////////////////////////////////////
  template<typename Scalar, 
	   int Rows, 
	   int Cols>
  Hist<Scalar, Rows, Cols>
  Hist<Scalar, Rows, Cols>::
  operator-=(const Hist & rhs)
  {
    this->fContents -= rhs.fContents;
    return *this;
  }

  /////////////////////////////////////////////////////////
  template<typename Scalar, 
	   int Rows, 
	   int Cols>
  Hist<Scalar, Rows, Cols>
  Hist<Scalar, Rows, Cols>::
  operator+(const Hist & rhs) const
  {
    Hist<Scalar, Rows, Cols> ret = *this; // copy this
    ret += rhs;
    return ret;
  }

  /////////////////////////////////////////////////////////
  template<typename Scalar, 
	   int Rows, 
	   int Cols>
  Hist<Scalar, Rows, Cols>
  Hist<Scalar, Rows, Cols>::
  operator+=(const Hist & rhs)
  {
    this->fContents += rhs.fContents;
    return *this;
  }

  /////////////////////////////////////////////////////////
  template<typename Scalar, 
	   int Rows, 
	   int Cols>
  Hist<Scalar, Rows, Cols>
  Hist<Scalar, Rows, Cols>::
  operator*(const Hist & rhs) const
  {
    Hist<Scalar, Rows, Cols> ret = *this; // copy this
    ret *= rhs;
    return ret;
  }

  /////////////////////////////////////////////////////////
  template<typename Scalar, 
	   int Rows, 
	   int Cols>
  Hist<Scalar, Rows, Cols>
  Hist<Scalar, Rows, Cols>::
  operator*=(const Hist & rhs)
  {
    this->fContents *= rhs.fContents;
    return *this;
  }

  /////////////////////////////////////////////////////////
  template<typename Scalar, 
	   int Rows, 
	   int Cols>
  Hist<Scalar, Rows, Cols>
  Hist<Scalar, Rows, Cols>::
  operator/(const Hist & rhs) const
  {
    Hist<Scalar, Rows, Cols> ret = *this; // copy this
    ret /= rhs;
    return ret;
  }

  /////////////////////////////////////////////////////////
  template<typename Scalar, 
	   int Rows, 
	   int Cols>
  Hist<Scalar, Rows, Cols>
  Hist<Scalar, Rows, Cols>::
  operator/=(const Hist & rhs)
  {
    this->fContents /= rhs.fContents;
    return *this;
  }
















  /////////////////////////////////////////////////////////
  template<typename Scalar, 
	   int Rows, 
	   int Cols>
  Hist<Scalar, Rows, Cols> 
  Hist<Scalar, Rows, Cols>::
  operator-(const Scalar & rhs) const
  {
    Hist<Scalar, Rows, Cols> ret = *this; // copy this
    ret -= rhs;
    return ret;
  }

  /////////////////////////////////////////////////////////
  template<typename Scalar, 
	   int Rows, 
	   int Cols>
  Hist<Scalar, Rows, Cols> 
  Hist<Scalar, Rows, Cols>::
  operator-=(const Scalar & rhs)
  {
    this->fContents -= rhs;
    return *this;
  }

  /////////////////////////////////////////////////////////
  template<typename Scalar, 
	   int Rows, 
	   int Cols>
  Hist<Scalar, Rows, Cols> 
  Hist<Scalar, Rows, Cols>::
  operator+(const Scalar & rhs) const
  {
    Hist<Scalar, Rows, Cols> ret = *this; // copy this
    ret += rhs;
    return ret;
  }

  /////////////////////////////////////////////////////////
  template<typename Scalar, 
	   int Rows, 
	   int Cols>
  Hist<Scalar, Rows, Cols> 
  Hist<Scalar, Rows, Cols>::
  operator+=(const Scalar & rhs)
  {
    this->fContents += rhs;
    return *this;
  }

  /////////////////////////////////////////////////////////
  template<typename Scalar, 
	   int Rows, 
	   int Cols>
  Hist<Scalar, Rows, Cols> 
  Hist<Scalar, Rows, Cols>::
  operator*(const Scalar & rhs) const
  {
    Hist<Scalar, Rows, Cols> ret = *this; // copy this
    ret *= rhs;
    return ret;
  }

  /////////////////////////////////////////////////////////
  template<typename Scalar, 
	   int Rows, 
	   int Cols>
  Hist<Scalar, Rows, Cols> 
  Hist<Scalar, Rows, Cols>::
  operator*=(const Scalar & rhs)
  {
    this->fContents *= rhs;
    return *this;
  }

  /////////////////////////////////////////////////////////
  template<typename Scalar, 
	   int Rows, 
	   int Cols>
  Hist<Scalar, Rows, Cols> 
  Hist<Scalar, Rows, Cols>::
  operator/(const Scalar & rhs) const
  {
    Hist<Scalar, Rows, Cols> ret = *this; // copy this
    ret /= rhs;
    return ret;
  }

  /////////////////////////////////////////////////////////
  template<typename Scalar, 
	   int Rows, 
	   int Cols>
  Hist<Scalar, Rows, Cols> 
  Hist<Scalar, Rows, Cols>::
  operator/=(const Scalar & rhs)
  {
    this->fContents /= rhs;
    return *this;
  }

  
  

}
