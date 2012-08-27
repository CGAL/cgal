// Copyright (c) 2005,2006  Stanford University (USA).
// 
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$ $Date$
// 
//
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_POLYNOMIAL_POLYNOMIAL_POLYCORE_KERNEL_H
#define CGAL_POLYNOMIAL_POLYNOMIAL_POLYCORE_KERNEL_H
#include <CGAL/Polynomial/basic.h>

#include <CGAL/CORE_Expr.h>
#include <CGAL/CORE/poly/Poly.h>
#include <CGAL/CORE_BigFloat.h>
#include <CGAL/CORE_Expr.h>
#include <CGAL/CORE_BigRat.h>

namespace CGAL { namespace POLYNOMIAL { namespace internal {

struct CORE_polynomial {
  typedef CORE::Polynomial<CORE::BigRat> P;
  typedef CORE::BigRat NT;

  template <class It>
  CORE_polynomial(It b, It e): p_(std::vector<NT> (b,e)) {
    /*for (int i=0; i<= degree(); ++i){
      CGAL_precondition(P::getCoeffi(i).err()==0);
      }*/
  }
  CORE_polynomial(){
  }
  CORE_polynomial(NT n): p_(0, &n){
  }
  CORE_polynomial(int d): p_(d){}
  CORE_polynomial(const P&p): p_(p){
    /*for (int i=0; i<= degree(); ++i){
      CGAL_precondition(p_.getCoeffi(i).err()==0);
      }*/
  }

  NT operator[](unsigned int i) const {
    return p_.getCoeffi(i);
  }

  /*NT &operator[](unsigned int i) {
    return p_.getCoeffi(i);
    }*/

  NT operator()(const NT &nt) const {
    return p_.eval(nt);
  }

  bool operator==(const CORE_polynomial&o ) const {
    if (p_.getTrueDegree() != o.p_.getTrueDegree()) {
      return false;
    } else {
      for (int i=0; i<= p_.getTrueDegree(); ++i) {
	if (operator[](i) != o[i]) return false;
      }
    }
    return true;
  }
  
  void contract()  {
    p_.contract();
  }

  bool operator!=(const CORE_polynomial&o ) const {
    return !operator==(o);
  }

  CORE_polynomial operator/(const NT &nt) const {
    P ret(p_.getTrueDegree());
    for (int i=0; i<= degree(); ++i){
      ret.setCoeff(i, operator[](i)/nt);
      //CGAL_assertion(rr);
    }
    return CORE_polynomial(ret);
  }

  CORE_polynomial operator-() const {
    P ret(p_.getTrueDegree());
    for (int i=0; i<= degree(); ++i){
      ret.setCoeff(i, -operator[](i));
      //CGAL_assertion(setr);
    }
    return CORE_polynomial(ret);
  }

  CORE_polynomial operator-(const CORE_polynomial &o) const {
    CORE_polynomial r(core_polynomial());
    r.p_-=o.core_polynomial();
    return r;
  }


  CORE_polynomial operator+(const CORE_polynomial &o) const {
    CORE_polynomial r(core_polynomial());
    r.p_+=o.core_polynomial();
    return r;
  }

  CORE_polynomial operator*(const CORE_polynomial &o) const {
    CORE_polynomial r(core_polynomial());
    r.p_*=o.core_polynomial();
    return r;
  }

  CORE::Expr operator()(const CORE::Expr &nt) const {
    return p_.eval(nt);
  }

  int degree() const {
    return p_.getTrueDegree();
  }

  const P &core_polynomial() const {
    return p_;
  }

  //! write it in maple format
  template <class C, class T>
  void write(std::basic_ostream<C,T> &out) const
  {
    std::basic_ostringstream<C,T> s;
    s.flags(out.flags());
    s.imbue(out.getloc());
    s.precision(12);
    if (degree()<0) {
      s << "0";
    }
    else {
      for (int i=0; i<= degree(); ++i) {
	if (i==0) {
	  if (p_.getCoeffi(i) != 0) s << p_.getCoeffi(i);
	}
	else {
	  if ( p_.getCoeffi(i) != 0 ) {
	    if (p_.getCoeffi(i) >0) s << "+";
	    s << p_.getCoeffi(i) << "*t";
	    if (i != 1) {
	      s << "^" << i;
	    }
	  }
	}
      }
    }

    out << s.str();
  }


 //! Read very stylized input
  template <class charT, class Traits>
  void read(std::basic_istream<charT, Traits> &in) {
    std::vector<NT> coefs;
    bool pos=(in.peek()!='-');
    if (in.peek() == '+' || in.peek() == '-') {
      char c;
      in >> c;
    }
    char v='\0';
    while (true) {
      char vc, t;
      NT coef;
      // eat
      in >> coef;
      //coef.makeExact();
      if (!in.good()) return;
      unsigned int pow;
      char p= in.peek();
      if (in.peek() =='*') {
	in >> t >> vc;
	if (t != '*') {
	  in.setstate(std::ios_base::failbit);
	  return;
	  //return in;
	}
	if ( !(v=='\0' || v== vc)) {
	  in.setstate(std::ios_base::failbit);
	  return;
	  //return in;
	}
	v=vc;
	p=in.peek();
	if (in.peek() =='^') {
	  char c;
	  in  >> c >> pow;
	}
	else {
	  pow=1;
	}
      }
      else {
	pow=0;
      }

      if (coefs.size() <=pow) {
	coefs.resize(pow+1);
      }

      if (!pos) coef=-coef;
      coefs[pow]=coef;

      char n= in.peek();
      if (n=='+' || n=='-') {
	pos= (n!='-');
	char e;
	in >> e;
      } else {
	/*bool f= in.fail();
	  bool g= in.good();
	  bool e= in.eof();
	  bool b= in.bad();*/
	// This is necessary since peek could have failed if we hit the end of the buffer
	// it would better to do without peek, but that is too messy
	in.clear();
	//std::cout << f << g << e << b<<std::endl;
	break;
      }
    }
    
    p_.operator=(P(coefs));
  }

  bool is_constant() const {
    return degree() < 1;
  }

  typedef const NT* iterator;
  iterator begin() const {
    return p_.coeff;
  }
  iterator end() const {
    return p_.coeff+p_.degree+1;
  }


protected:
  P p_;
};

template < class C, class Tr>
inline std::ostream &operator<<(std::basic_ostream<C, Tr> &out,
				const CORE_polynomial &poly)
{
  poly.write(out);
  return out;
}

template <class C, class Tr>
inline std::istream &operator>>(std::basic_istream<C, Tr> &in,
				CORE_polynomial &poly)
{
  poly.read(in);
  return in;
}

CORE_polynomial operator*(const CORE_polynomial::NT &a,
			  const CORE_polynomial &p){
  //CORE_polynomial::NT ac(a);
  return CORE_polynomial(CORE_polynomial::P(0, &a)*p.core_polynomial());
}

} } } //namespace CGAL::POLYNOMIAL::internal
#endif
