#ifndef CGAL_POLYNOMIAL_DEFAULT_FILTERING_TRAITS_H
#define CGAL_POLYNOMIAL_DEFAULT_FILTERING_TRAITS_H
#include <CGAL/Polynomial/basic.h>
#include <CGAL/Polynomial/Polynomial.h>
#include <CGAL/Polynomial/Interval_polynomial.h>
#include <CGAL/Polynomial/polynomial_converters.h>

#ifdef POLYNOMIAL_USE_CGAL
#include <CGAL/Gmpq.h>
#define DEFAULT_FILTERING_DEFAULT_NT =CGAL::Gmpq
#endif

CGAL_POLYNOMIAL_BEGIN_NAMESPACE
template <class NT DEFAULT_FILTERING_DEFAULT_NT>
struct Default_filtering_traits {
  typedef Polynomial<NT> Exact_function;
  typedef Interval_polynomial Interval_function;
  typedef Polynomial_converter<Exact_function, 
			       Interval_function, 
			       To_interval<NT> > Exact_to_interval_converter;
};

CGAL_POLYNOMIAL_END_NAMESPACE

#ifdef DEFAULT_FILTERING_DEFAULT_NT
#undef DEFAULT_FILTERING_DEFAULT_NT
#endif
#endif
