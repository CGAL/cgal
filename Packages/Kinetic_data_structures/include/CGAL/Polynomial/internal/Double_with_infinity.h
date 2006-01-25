#ifndef CGAL_POLYNOMIAL_INTERNAL_WRAPPED_DOUBLE_H
#define CGAL_POLYNOMIAL_INTERNAL_WRAPPED_DOUBLE_H
#include <CGAL/Polynomial/basic.h>

CGAL_POLYNOMIAL_BEGIN_INTERNAL_NAMESPACE

class Double_with_infinity {
  Double_with_infinity(double d): d_(d){}
  operator double() const {
    return d_;
  }
protected:
  double d_;
};

CGAL_POLYNOMIAL_END_INTERNAL_NAMESPACE

namespace std {
  template <>
  struct numeric_limits<CGAL_POLYNOMIAL_NS::internal::Double_with_infinity >: public numeric_limits<double>
  {
    static const bool is_specialized = true;
    static const bool has_infinity=true;
    static double infinity() throw() {return std::numeric_limits<double>::max();}
  };
};

#endif
