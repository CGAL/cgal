#ifndef CGAL_POLYNOMIAL_INTERNAL_QUOTIENT_H
#define CGAL_POLYNOMIAL_INTERNAL_QUOTIENT_H

#include <CGAL/Polynomial/basic.h>
#include <CGAL/Polynomial/internal/Rational/Quotient_remainder.h>
/*!
  \file Quotient.h A class to compute quotients of polynomials.
*/

CGAL_POLYNOMIAL_BEGIN_INTERNAL_NAMESPACE

//! Compute the quotient of two polynomials.
/*!
  I pulled this out of Polynomial because I did not think polynomial should have such complicated methods. 
*/
template<class Polynomial>
class Quotient : private Quotient_remainder<Polynomial>
{
private:
  typedef Quotient_remainder<Polynomial>  Base;

public:
  typedef typename Polynomial::NT   NT;
  typedef Polynomial       result_type;
  typedef Polynomial       argument_type;
  typedef Polynomial       argument_type1;
  typedef Polynomial       argument_type2;

  void write(std::ostream &out) const {
    out << "quo";
  }

  //! compute the quotient
  result_type
  operator()(const Polynomial& t, const Polynomial& v) const
  {
    return Base::operator()(t,v).first;
  }
};


CGAL_POLYNOMIAL_END_INTERNAL_NAMESPACE

#endif // CGAL_POLYNOMIAL_INTERNAL_QUOTIENT_H
