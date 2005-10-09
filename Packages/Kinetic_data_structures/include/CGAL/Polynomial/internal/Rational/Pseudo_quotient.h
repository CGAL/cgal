#ifndef CGAL_POLYNOMIAL_INTERNAL_PSEUDO_QUOTIENT_H
#define CGAL_POLYNOMIAL_INTERNAL_PSEUDO_QUOTIENT_H

#include <CGAL/Polynomial/basic.h>
#include <vector>

CGAL_POLYNOMIAL_BEGIN_INTERNAL_NAMESPACE
template<class Polynomial>
struct Pseudo_quotient
{
  typedef typename Polynomial::NT   NT;
  typedef Polynomial       result_type;
  typedef Polynomial       argument_type;
  typedef Polynomial       argument_type1;
  typedef Polynomial       argument_type2;

  void write(std::ostream &out) const {
    out << "pquo";
  }

  Polynomial
  operator()(const Polynomial& t, const Polynomial& v) const
  {
   Polynomial_precondition( t.degree() >= v.degree() );
   Polynomial_precondition( !v.is_zero() );
   
   int m = t.degree();
   int n = v.degree();
   int divdeg = m - n;

   std::vector<NT> r_coef(m+1);
   std::vector<NT> q_coef(divdeg+1);

    // compute the powers of v[n] from 0 to m-n
    std::vector<NT> pow(divdeg+1);
    pow[0] = NT(1);
    for (int i = 1; i <= divdeg; i++) {
      pow[i] = pow[i-1] * v[n];
    }

    // initialize r
    for (int i = 0; i < divdeg; i++) {
      r_coef[i] = t[i] * pow[divdeg - i];
    }
    for (int i = divdeg; i <= m; i++) {
      r_coef[i] = t[i];
    }

    // compute the pseudo-quotient
    for (int k = divdeg; k >= 0; k--) {
      q_coef[k] = r_coef[n + k] * pow[k];
      for (int j = n + k - 1; j >= k; j--) {
	r_coef[j] = v[n] * r_coef[j] - r_coef[n + k] * v[j - k];
      }
    }

    pow.clear();

    Polynomial q(q_coef.begin(), q_coef.end());
    return q;
  }
};


CGAL_POLYNOMIAL_END_INTERNAL_NAMESPACE
#endif
