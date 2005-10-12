#ifndef CGAL_POLYNOMIAL_TOOLS_INTERVAL_ARITHMETIC_H
#define CGAL_POLYNOMIAL_TOOLS_INTERVAL_ARITHMETIC_H
#include <CGAL/Polynomial/basic.h>

/*!
  \file interval_arithmetic.h

  \todo wrapper support for boost::interval
*/

#ifdef CGAL_POLYNOMIAL_USE_CGAL
#include <CGAL/Interval_arithmetic.h>

CGAL_POLYNOMIAL_BEGIN_NAMESPACE
//! The class we use for interval arithmetic.
/*!
  This class is not protected, so the Interval_arithmetic_guard must be used when calculations
  are peformed. This is not checked. 
*/
typedef CGAL::Interval_nt_advanced Interval_nt;

//! This class sets the rounding modes for interval arithmetic.
/*!
  Create an instance of this class when you want to perform interval arithmetic.
  When the destructor is called, the mode will be automatically cleaned up. 
*/
struct Interval_arithmetic_guard {
  //! I should check that the rounding is not already upwards.
  Interval_arithmetic_guard() {
    bk_= FPU_get_and_set_cw(CGAL_FE_UPWARD);
  };
  Interval_arithmetic_guard(bool b) {
    if (b){
      bk_= FPU_get_and_set_cw(CGAL_FE_UPWARD);
    } else {
      bk_= CGAL_FE_UPWARD;
    }
  }
  ~Interval_arithmetic_guard(){
    FPU_set_cw(bk_);
  }
  bool enabled() const {
    return bk_== CGAL_FE_UPWARD;
  }
  void set_enabled(bool ft) {
    if (ft != enabled()) {
      bk_= FPU_get_and_set_cw(bk_);
    }
  }
protected:
  /*!
    Store the rounding mode to restore here. This could probably be made static. 
    Making it static might have consequences for multithreaded programs, I am not sure.
  */
  FPU_CW_t bk_;
};

template <class NT>
class To_interval: public CGAL::To_interval<NT> {
public:
  To_interval(){}
};



/*template <class NT>
std::pair<double, double> to_interval(const NT &nt){
  //bool to_interval_general;
  return CGAL::to_interval(nt);
  }*/
//using CGAL::to_interval;

CGAL_POLYNOMIAL_END_NAMESPACE

#elif POLYNOMIAL_USE_BOOST_INTERVAL

Not implemented yet.

#else

No interval arithmetic support.

#endif

CGAL_POLYNOMIAL_BEGIN_NAMESPACE

inline Extended_sign extended_sign(const Interval_nt &i) {
  if (i.inf() == 0 && i.sup() ==0) return EXTENDED_ZERO;
  else if (i.inf() <=0 && i.sup() >=0) {
    return EXTENDED_UNKNOWN;
  } else if (i.inf() >0 ){
    return EXTENDED_POSITIVE;
  } else {
    return EXTENDED_NEGATIVE;
  }
}

CGAL_POLYNOMIAL_END_NAMESPACE
#endif
