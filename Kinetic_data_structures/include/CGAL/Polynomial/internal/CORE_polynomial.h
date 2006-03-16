// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: svn+ssh://scm.gforge.inria.fr/svn/cgal/trunk/Kinetic_data_structures/include/CGAL/Polynomial/Kernel.h $
// $Id: Kernel.h 28489 2006-02-14 10:08:15Z lsaboret $ $Date: 2006-02-14 02:08:15 -0800 (Tue, 14 Feb 2006) $
// 
//
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_POLYNOMIAL_POLYNOMIAL_POLYCORE_KERNEL_H
#define CGAL_POLYNOMIAL_POLYNOMIAL_POLYCORE_KERNEL_H
#include <CGAL/Polynomial/basic.h>

#include <CGAL/CORE_Expr.h>
#include <CORE/poly/Poly.h>
#include <CORE/BigFloat.h>
#include <CORE/Expr.h>
#include <CORE/BigRat.h>

CGAL_POLYNOMIAL_BEGIN_INTERNAL_NAMESPACE

struct CORE_polynomial: public CORE::Polynomial<CORE::BigRat> {
  typedef CORE::Polynomial<CORE::BigRat> P;
  typedef CORE::BigRat NT;

  template <class It>
  CORE_polynomial(It b, It e): P(std::vector<NT>(b,e)) {
    /*for (int i=0; i<= degree(); ++i){
      CGAL_precondition(P::getCoeffi(i).err()==0);
      }*/
  }
  CORE_polynomial(){
  }
  CORE_polynomial(NT n): P(0, &n){
  }
  CORE_polynomial(int d): P(d){}
  CORE_polynomial(const P&p): P(p){
    /*for (int i=0; i<= degree(); ++i){
      CGAL_precondition(P::getCoeffi(i).err()==0);
      }*/
  }

  NT operator[](unsigned int i) const {
    return P::getCoeffi(i);
  }

  /*NT &operator[](unsigned int i) {
    return P::getCoeffi(i);
    }*/

  NT operator()(const NT &nt) const {
    return P::eval(nt);
  }

  bool operator==(const CORE_polynomial&o ) const {
    if (P::getTrueDegree() != o.getTrueDegree()) {
      return false;
    } else {
      for (int i=0; i<= P::getTrueDegree(); ++i) {
	if (operator[](i) != o[i]) return false;
      }
    }
    return true;
  }

  CORE_polynomial operator/(const NT &nt) const {
    CORE_polynomial ret(P::getTrueDegree());
    for (int i=0; i<= degree(); ++i){
      ret.setCoeff(i, operator[](i)/nt);
    }
    return ret;
  }

  CORE_polynomial operator-() const {
    CORE_polynomial ret(P::getTrueDegree());
    for (int i=0; i<= degree(); ++i){
      ret.setCoeff(i, -operator[](i));
    }
    return ret;
  }

  CORE::Expr operator()(const CORE::Expr &nt) const {
    return P::eval(nt);
  }

  int degree() const {
    return P::getTrueDegree();
  }

  const P &core_polynomial() const {
    return *this;
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
	  if (P::getCoeffi(i) != 0) s << P::getCoeffi(i);
	}
	else {
	  if ( P::getCoeffi(i) != 0 ) {
	    if (P::getCoeffi(i) >0) s << "+";
	    s << P::getCoeffi(i) << "*t";
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
    
    P::operator=(P(coefs));
  }

  bool is_constant() const {
    return degree() < 1;
  }

  typedef const NT* iterator;
  iterator begin() const {
    return P::coeff;
  }
  iterator end() const {
    return P::coeff+P::degree+1;
  }

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

CGAL_POLYNOMIAL_END_INTERNAL_NAMESPACE
#endif
