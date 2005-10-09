#ifndef CGAL_POLYNOMIAL_INTERNAL_EVALUATE_H
#define CGAL_POLYNOMIAL_INTERNAL_EVALUATE_H
#include <CGAL/Polynomial/basic.h>
#include <CGAL/Polynomial/internal/numeric_solvers.h>
#include <vector>

CGAL_POLYNOMIAL_BEGIN_INTERNAL_NAMESPACE

template <class NT>
inline NT evaluate_polynomial(const std::vector<NT> &coefs, const NT &t) {
  if (coefs.empty()) return NT(0);
  typename std::vector<NT>::const_reverse_iterator rit= coefs.rbegin();
  NT result = *rit;
  ++rit;
  for (; rit != coefs.rend(); ++rit){
    result *= t;
    result += (*rit);
  }
  return result;
}

inline double evaluate_polynomial(const std::vector<double>& coefs, double t) {
  if (coefs.empty()) return 0;
  return evaluate_polynomial(&coefs.front(), &coefs.front()+coefs.size(), t);
}


template<class K, class NT >
class Evaluate_polynomial {
  typedef typename K::Function            Polynomial;

 

public:
  Evaluate_polynomial(){}
  Evaluate_polynomial(const Polynomial& p): coefs_(p.begin(), p.end()) {
    /*coefs_.resize(p.degree()+1);
    for (unsigned int i=0; i< p.degree(); ++i){
      coefs_[i]= NT(p[i]);
      }*/
  }

  typedef NT argument_type;
  typedef NT result_type;

  result_type operator()(const argument_type& x) const {
    return evaluate_polynomial(coefs_, x);
  }

protected:
  std::vector<result_type> coefs_;
};


CGAL_POLYNOMIAL_END_INTERNAL_NAMESPACE

#endif
