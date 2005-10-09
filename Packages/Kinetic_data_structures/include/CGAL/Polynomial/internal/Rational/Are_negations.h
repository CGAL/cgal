#ifndef CGAL_POLYNOMIAL_INTERNAL_ARE_NEGATIONS_H
#define CGAL_POLYNOMIAL_INTERNAL_ARE_NEGATIONS_H

#include <CGAL/Polynomial/basic.h>

CGAL_POLYNOMIAL_BEGIN_INTERNAL_NAMESPACE

//------------------------------------------------------------------
template <class P>
class Are_negations {
public:
  typedef typename P::NT NT;
  typedef P first_argument_type;
  typedef P second_argument_type;
  typedef bool result_type;
  Are_negations(){}

  result_type operator()(const first_argument_type &f0, 
			 const second_argument_type &f1) const {
    if (f0.degree() != f1.degree()) return false;
    else {
      for (int i=0; i<= f0.degree(); ++i){
	if (f0[i] != -f1[i]) return false;
      }
      return true;
    }
  }
};

CGAL_POLYNOMIAL_END_INTERNAL_NAMESPACE
#endif
