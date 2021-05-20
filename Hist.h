#include <Eigen/Dense>

namespace xsec {
  // Histograming object used internally by the framework
  // Wraps Eigen arrays
  // Users can create conversion functions to/from this object
  // template parameters are forwarded to the underlaying Eigen arrays
  template<typename Scalar, 
	   int Rows, 
	   int Cols>
  class Hist {
  public:
    Hist(const Eigen::Array<Scalar, Rows, Cols> & contents,
	 const Eigen::Array<Scalar, Rows, Cols> & Edges);


    Hist operator-(const Hist& rhs) const;
    Hist operator+(const Hist& rhs) const;
    Hist operator/(const Hist& rhs) const;
    Hist operator*(const Hist& rhs) const;

    Hist operator-=(const Hist& rhs);
    Hist operator+=(const Hist& rhs);
    Hist operator/=(const Hist& rhs);
    Hist operator*=(const Hist& rhs);

    Hist operator-(const Scalar& rhs) const;
    Hist operator+(const Scalar& rhs) const;
    Hist operator/(const Scalar& rhs) const;
    Hist operator*(const Scalar& rhs) const;

    Hist operator-=(const Scalar& rhs);
    Hist operator+=(const Scalar& rhs);
    Hist operator/=(const Scalar& rhs);
    Hist operator*=(const Scalar& rhs);
    
  private:
    Eigen::Array<Scalar, Rows, Cols> fEdges;
    Eigen::Array<Scalar, Rows, Cols> fContents;

  };

}
