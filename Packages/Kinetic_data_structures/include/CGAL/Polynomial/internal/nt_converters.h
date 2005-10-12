#ifndef CGAL_POLYNOMIAL_NT_CONVERTERS_H
#define CGAL_POLYNOMIAL_NT_CONVERTERS_H

#include <CGAL/Polynomial/basic.h>

#ifdef CGAL_POLYNOMIAL_USE_CGAL
#include <CGAL/NT_converter.h>
#include <CGAL/number_utils_classes.h>

CGAL_POLYNOMIAL_BEGIN_NAMESPACE 

template <class NT1, class NT2>
class NT_converter: public CGAL::NT_converter<NT1, NT2> {
public:
  NT_converter(){}
};

template <class NT>
class To_double: public CGAL::To_double<NT> {
public:
  To_double(){}
};

/*template <class NT> 
double to_double(const NT &nt) {
  return CGAL::to_double(nt);
  }*/

CGAL_POLYNOMIAL_END_NAMESPACE




#else

Not implemented yet;

#endif

CGAL_POLYNOMIAL_BEGIN_NAMESPACE 

//! This does not use any CGAL code.
template <class NT>
struct Identity_converter {
  typedef NT argument_type;
  typedef NT result_type;
  Identity_converter(){}
  const NT &operator()(const NT &nt) const {
    return nt;
  }
};

CGAL_POLYNOMIAL_END_NAMESPACE
#endif
