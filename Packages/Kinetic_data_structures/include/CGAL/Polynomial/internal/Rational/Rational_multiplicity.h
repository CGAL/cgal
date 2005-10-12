#ifndef CGAL_POLYNOMIAL_INTERNAL_ROOT_MULTIPLICITY_H
#define CGAL_POLYNOMIAL_INTERNAL_ROOT_MULTIPLICITY_H

#include <CGAL/Polynomial/basic.h>
#include <vector>

CGAL_POLYNOMIAL_BEGIN_INTERNAL_NAMESPACE

template <class K>
class Rational_multiplicity {
public:
  Rational_multiplicity(){}
  
  Rational_multiplicity(const typename K::Function &fh, const K &k): d_(k.differentiate_object()), 
								     k_(k){
    CGAL_Polynomial_precondition(fh.degree() != -1);
    h_.push_back(fh);
  }

  typedef unsigned int result_type;
  //typedef Bound_type argument_type;
  typedef typename K::Function::NT argument_type;
  
  result_type operator()(const argument_type &t) const {
    CGAL_Polynomial_exactness_assertion(k_.sign_at_object(h_.front())( t)== CGAL_POLYNOMIAL_NS::ZERO);
    //POLYNOMIAL_NS::Sign sn;
    //if ( k.sign_at_object(fh)(t) != POLYNOMIAL_NS::ZERO ) return 0;
    unsigned int deg=1;
    unsigned int mdegree= h_.front().degree();
    while (sign_at_of(t,deg)==CGAL_POLYNOMIAL_NS::ZERO && deg <mdegree){
      ++deg;
    }
    return deg;
  }

protected:
  CGAL_POLYNOMIAL_NS::Sign sign_at_of(const argument_type &t, unsigned int i) const {
    if (i >= h_.size()) {
      h_.push_back(d_(h_.back()));
    }
    CGAL_Polynomial_postcondition(i < h_.size());
    typename K::Sign_at sa= k_.sign_at_object(h_[i]);
    return sa(t);
  }
  
  mutable std::vector<typename K::Function> h_;
  typename K::Differentiate d_;
  K k_;
};

CGAL_POLYNOMIAL_END_INTERNAL_NAMESPACE
#endif
