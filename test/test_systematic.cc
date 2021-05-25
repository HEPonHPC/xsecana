#include <iostream>

#include "XSecAna/Hist.h"
#include "XSecAna/Systematic.h"

#include <Eigen/Dense>

using namespace xsec;
#define TEST_ARRAY(test_name, arr1, arr2)				\
  if(!(arr1 - arr2).isZero(0)) {					\
    std::cerr << "test_systematic " << test_name << ": FAILED" << std::endl; \
    std::cerr << "test_systematic " << test_name << "\t" << arr1 << std::endl; \
    std::cerr << "test_systematic " << test_name << "\t" << arr2 << std::endl; \
    pass = false;							\
  }									

#define TEST_SYSTEMATIC(test_name, syst, up, down)			\
  if(*syst.GetShifts().first != (up)) {					\
    std::cerr << "test_systematic " << test_name << ": FAILED" << std::endl; \
    std::cerr << "test_systematic " << test_name << "\t" << syst.GetShifts().first->Contents() << std::endl; \
    std::cerr << "test_systematic " << test_name << "\t" << (up).Contents() << std::endl; \
    pass = false;							\
  }									\
  if(*syst.GetShifts().second != (down)) {				\
    std::cerr << "test_systematic " << test_name << ": FAILED" << std::endl; \
    std::cerr << "test_systematic " << test_name << "\t" << syst.GetShifts().second->Contents() << std::endl; \
    std::cerr << "test_systematic " << test_name << "\t" << (down).Contents() << std::endl; \
    pass = false;							\
  }									

#define LIFT(F) \
  ([](auto &&... xs) ->decltype(auto) { \
    return F(std::forward<decltype(xs)>(xs)...);	\
  })

int main(int argc, char ** argv)
{
  auto bins = Eigen::Array<double, 1, 11>::LinSpaced(11, 0, 10);
  
  auto vnominal = Eigen::Array<double, 1, 10>::Ones() * 5;
  auto vup      = Eigen::Array<double, 1, 10>::LinSpaced(10, 4, 6);
  auto vdown    = Eigen::Array<double, 1, 10>::LinSpaced(10, 4, 6).reverse();

  Hist<double, 10> nominal(vnominal, bins);
  Hist<double, 10> up     (vup     , bins);
  Hist<double, 10> down   (vdown   , bins);
 
  double pass = true;

  // two sided sytematic construction
  Systematic<Hist<double, 10> > syst_2("syst", up, down);
  TEST_SYSTEMATIC("construction", syst_2, up, down);

  // two sided systematic subtraction via invoke
  typedef Hist<double, 10> histtype;
  Systematic<Hist<double, 10> > syst_diff_2 = 
    syst_2.Invoke(static_cast<histtype(histtype::*)(const histtype&) const>(&histtype::operator-), nominal);
  TEST_SYSTEMATIC("two sided subtraction", syst_diff_2, up - nominal, down - nominal);

  // two sided systematic division via invoke
  Systematic<Hist<double, 10> > syst_div_2 = 
    syst_diff_2.Invoke(static_cast<histtype(histtype::*)(const histtype&) const>(&histtype::operator/), nominal);
  TEST_SYSTEMATIC("two sided division", syst_div_2, (up - nominal) / nominal, (down - nominal) / nominal);
  
  // one sided systematic construction
  Systematic<Hist<double, 10> > syst_1("syst", up);
  TEST_SYSTEMATIC("construction", syst_1, up, up);

  // one sided systematic subtraction via invoke
  Systematic<Hist<double, 10> > syst_diff_1 = 
    syst_1.Invoke(static_cast<histtype(histtype::*)(const histtype&) const>(&histtype::operator-), nominal);
  TEST_SYSTEMATIC("one sided subtraction", syst_diff_1, up - nominal, up - nominal);


  // one sided systematic division via invoke
  Systematic<Hist<double, 10> > syst_div_1 = 
    syst_diff_1.Invoke(static_cast<histtype(histtype::*)(const histtype&) const>(&histtype::operator/), nominal);
  TEST_SYSTEMATIC("one sided division", syst_div_1, (up - nominal) / nominal, (up - nominal) / nominal);

  if(pass) std::cout << "Success!" << std::endl;

  
}
