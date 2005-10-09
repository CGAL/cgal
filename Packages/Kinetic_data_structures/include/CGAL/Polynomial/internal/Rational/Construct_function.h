#ifndef CGAL_POLYNOMIAL_INTERNAL_CONSTRUCT_FUNCTION_H
#define CGAL_POLYNOMIAL_INTERNAL_CONSTRUCT_FUNCTION_H
#include <CGAL/Polynomial/basic.h>

CGAL_POLYNOMIAL_BEGIN_INTERNAL_NAMESPACE

template <class Fn>
struct Construct_function{

  typedef Fn result_type;
  typedef typename result_type::NT argument_type;
    
  template <class It>
  result_type operator()(It b, 
			 It e) const {
    return result_type(b,e);
  }

  //! construct high degree polynomials
  result_type operator()(const argument_type &a0, 
			 const argument_type &a1=0) const {
    argument_type v[] = {a0, a1};
    return result_type(v, v+2);
  }
    
    
  //! construct high degree polynomials
  result_type operator()(const argument_type &a0, 
			 const argument_type &a1, 
			 const argument_type &a2, 
			 const argument_type &a3=0) const {
    argument_type v[] = {a0, a1, a2, a3};
    return result_type(v, v+4);
  }
    
    
  //! construct high degree polynomials
  result_type operator()(const argument_type &a0, 
			 const argument_type &a1, 
			 const argument_type &a2, 
			 const argument_type &a3,
			 const argument_type &a4, 
			 const argument_type &a5=0,
			 const argument_type &a6=0,
			 const argument_type &a7=0) const {
    argument_type v[] = {a0, a1, a2, a3, a4, a5, a6, a7};
    return result_type(v, v+8);
  }



  //! construct high degree polynomials
  result_type operator()(const argument_type &a0,
			 const argument_type &a1, 
			 const argument_type &a2, 
			 const argument_type &a3,
			 const argument_type &a4,
			 const argument_type &a5,
			 const argument_type &a6,
			 const argument_type &a7,
			 const argument_type &a8,
			 const argument_type &a9=0,
			 const argument_type &a10=0,
			 const argument_type &a11=0,
			 const argument_type &a12=0,
			 const argument_type &a13=0,
			 const argument_type &a14=0,
			 const argument_type &a15=0,
			 const argument_type &a16=0,
			 const argument_type &a17=0,
			 const argument_type &a18=0,
			 const argument_type &a19=0) const {
    argument_type v[] = {a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10,
			 a11, a12, a13, a14, a15, a16, a17, a18, a19};
    return result_type(v, v+20);
  }
};


CGAL_POLYNOMIAL_END_INTERNAL_NAMESPACE

#endif
