// Copyright (c) 2005  Stanford University (USA).
// All rights reserved.
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
// $Id$
// 
//
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_POLYNOMIAL_INTERNAL_POLYNOMIAL_CORE_H_
#define CGAL_POLYNOMIAL_INTERNAL_POLYNOMIAL_CORE_H_

#include <CGAL/Polynomial/basic.h>
#include <vector>
#include <iostream>
#include <CGAL/Polynomial/internal/Rational/Evaluate_polynomial.h>
#include <sstream>
#include <algorithm>

#define CGAL_EXCESSIVE(x)

namespace CGAL { namespace POLYNOMIAL { namespace internal {

template <class This, class NT_t>
class Polynomial_impl;

template < class T, class NT, class C, class Tr>
inline std::ostream &operator<<(std::basic_ostream<C,Tr> &out,
				const Polynomial_impl<T, NT> &poly);

//! A basic polynomial class
/*!  This handles everything having to do with polynomials which is to
  not too numerically sensitive.  The poly is stored as a
  vector. There are no 0 top coefficients in the vector.

  I don't store a value for 0 (there is not stored offset) since the
  polynomial will in general have to be recomputed when any
  certificate is computed, and I don't see any way of benefiting from
  delaying.

  If the STRIP_ZEROS flag is true then leading zeros are stripped on
  the fly. Otherwise you have to make sure that there are no leading
  0s yourself. However, stripping leading 0s causes problems with an
  interval is used as the number type.

*/

//! todo: check resize and doubles

template <class This, class NT_t>
class Polynomial_impl
{
  typedef std::vector<NT_t>                  Coefficients;
public:
  typedef NT_t                             NT;
  typedef typename Coefficients::const_iterator  iterator;
  typedef NT                               result_type;
  typedef NT                               argument_type;

  //================
  // CONSTRUCTORS
  //================

  //! Default
  Polynomial_impl() {}

  //! Make a constant polynomial
  Polynomial_impl(const NT& c) {
    coefs_.push_back(c);
  }

  //! Get coefficients from a vector
  template <class Iterator>
  Polynomial_impl(Iterator first, Iterator beyond): coefs_(first, beyond) {}

  //========================
  // ACCESS TO COEFFICIENTS
  //========================

  //! Return the ith coefficient.
  /*!
    This must return a value for any value of i.
  */
  const NT& operator[](unsigned int i) const
  {
    CGAL_assertion( i < coefs_.size());
    //if (i < coefs_.size()) {
      return coefs_[i];
      /*}
    else {
      return zero_coef();
      }*/
  }

  //=====================================
  //     ITERATORS FOR COEFFICIENTS
  //=====================================

  iterator begin() const
  {
    return coefs_.begin();
  }

  iterator end() const
  {
    return coefs_.end();
  }

  //=========
  // DEGREE
  //=========

  //! For the more mathematical inclined (as opposed to size());
  int degree() const
  {
    return static_cast<int>(coefs_.size()) - 1;
  }

  //=============
  // PREDICATES
  //=============

  //! Returns true if the polynomial is a constant function
  bool is_constant() const
  {
    return degree() < 1;
  }

  //! Returns true if f(t)=0

  bool is_zero() const
  {
    CGAL_Polynomial_assertion( coefs_.empty() == (coefs_.size()==0) );
    return coefs_.empty();
  }

  //=======================
  //     OPERATORS
  //=======================

  //! negation
  This operator-() const
  {
    This ret;
    ret.coefs_.resize( coefs_.size() );
    for (unsigned int i=0; i < coefs_.size(); ++i) {
      ret.coefs_[i] = -coefs_[i];
    }
    return ret;
  }

  //! polynomial addition
  This operator+(const This &o) const
  {
    if (is_zero()) { return o; }
    else if (o.is_zero()) { return This(*this); }
    else {
      This ret;
      unsigned int new_deg = (std::max)(o.degree(), degree());
      ret.coefs_.resize(new_deg + 1);
      unsigned int md= (std::min)(degree(), o.degree());
      for (unsigned int i = 0; i <= md; ++i) {
	ret.coefs_[i] = operator[](i) + o[i];
      }
      for (int i=md+1; i <= degree(); ++i){
	ret.coefs_[i]= operator[](i);
      }
      for (int i=md+1; i <= o.degree(); ++i){
	ret.coefs_[i]= o[i];
      }
      CGAL_EXCESSIVE(std::cout << *this << " - " << o << " + " << ret << std::endl);
      ret.finalize();
      return ret;
    }
  }

  //! polynomial subtraction
  This operator-(const This &o) const
  {
    if (is_zero()) { return -o; }
    else if (o.is_zero()) { return This(*this); }
    else {
      This ret;
      unsigned int new_deg = (std::max)(o.degree(), degree());
      ret.coefs_.resize( new_deg + 1 );
      unsigned int md= (std::min)(degree(), o.degree());
      for (unsigned int i = 0; i <= md; ++i) {
	ret.coefs_[i] = operator[](i) - o[i];
      }
      for (int i=md+1; i <= degree(); ++i){
	ret.coefs_[i]= operator[](i);
      }
      for (int i=md+1; i <= o.degree(); ++i){
	ret.coefs_[i]= -o[i];
      }
      CGAL_EXCESSIVE(std::cout << *this << " - " << o << " = " << ret << std::endl);
      ret.finalize();
      return ret;
    }
  }

  //! polynomial multiplication
  This operator*(const This &o) const
  {
    if (is_zero()) return This(NT(0));
    else if (o.is_zero()) return o;
    This ret;
    ret.coefs_.resize( degree() + o.degree() + 1 );
    // the following for-loop makes sure that the values on ret.coefs_
    // are properly initialized with zero.
    for (int i = 0; i <= degree() + o.degree(); ++i) {
      ret.coefs_[i] = NT(0);
    }
    for (unsigned int i = 0; i < coefs_.size(); ++i) {
      for (unsigned int j = 0; j < o.coefs_.size(); ++j) {
	NT prev = ret.coefs_[i+j];
	NT result = prev + operator[](i) * o[j];
	ret.coefs_[i+j] = result;
      }
    }
    CGAL_EXCESSIVE(std::cout << *this << " * " << o << " = " << ret << std::endl);
    return ret;
  }

  //! add a scalar
  This operator+(const NT& a) const
  {
    if ( is_zero() ) { return This(a); }

    This res(*this);
    res.coefs_[0] += a;
    CGAL_EXCESSIVE(std::cout << *this << " + " << a << " = " << res << std::endl);
    return res;
  }

  //! subtract a scalar
  This operator-(const NT& a) const
  {
    return (*this) + (-a);
  }

  //! multiply with scalar
  This operator*(const NT& a) const
  {
    if ( is_zero() || CGAL_NTS sign(a)==ZERO ) { return This(); }

    This res;
    unsigned int deg = degree();
    res.coefs_.resize(deg + 1, NT(0));
    for (unsigned int i = 0; i <= deg; i++) {
      res.coefs_[i] = coefs_[i] * a;
    }
    return res;
  }

  //! divide by scalar
  This operator/(const NT& a) const
  {
    NT inv_a = NT(1) / a;
    return (*this) * inv_a;
  }

  //=======================
  // VALUE OF POLYNOMIAL
  //=======================

  //! Evaluate the polynomial at some value.
  /*!
    This is primarily for the solvers to call when they know the
    number type is exact. Note, this method performs a construction
    and is unreliable if an inexact number type is used
  */
  NT operator()(const NT &t) const
  {
    if (is_zero()) return NT(0);
    else if (degree()==0 /* || t==NT(0)*/ ) {
      return operator[](0);
    }
    else {
      return evaluate_polynomial(coefs_, t);
    }
  }


  //! The non-operator version of operator()
  /*!
   */
  NT value_at(const NT &t) const
  {
    return operator()(t);
  }


  //! !operator==
  bool operator!=(const This &o) const
  {
    return !operator==(o);
  }


  //! check if the coefficients are equal
  bool operator==(const This &o) const
  {
    if (degree() != static_cast<int>(o.degree())) return false;
    int max_size = static_cast<int>((std::max)(o.coefs_.size(), coefs_.size()));
    for (int i = 0; i < max_size; ++i) {
      if (o[i] != operator[](i)) return false;
    }
    return true;
  }


  void print()const
  {
    write(std::cout);
  }


  //! Read very stylized input
  template <class charT, class Traits>
  void read(std::basic_istream<charT, Traits> &in) {
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
      if (in.fail()) {
	//std::cout << "Read ";
	//write(std::cout);
	//std::cout << std::endl;
	return;
      }
      //std::cout << "Read " << coef << std::endl;

      //NT cp1= coef+NT(1);
      unsigned int pow;
      in.peek();
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
        in.peek();
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

      if (coefs_.size() <=pow) {
	coefs_.resize(pow+1, NT(0));
      }

      if (!pos) coef=-coef;
      //std::cout << "Read 2 " << coef << "*t^" << pow << std::endl;
      coefs_[pow]=coef;
      //std::cout << "Read " << coefs_[pow] << std::endl;

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
    };
  }


  //! write it in maple format
  template <class C, class T>
  void write(std::basic_ostream<C,T> &out) const
  {
    std::basic_ostringstream<C,T> s;
    s.flags(out.flags());
    s.imbue(out.getloc());
    s.precision(12);
    if (coefs_.size()==0) {
      s << "0";
    }
    else {
      for (unsigned int i=0; i< coefs_.size(); ++i) {
	if (i==0) {
	  if (coefs_[i] != 0) s << coefs_[i];
	}
	else {
	  if ( coefs_[i] != 0 ) {
	    if (coefs_[i] >0) s << "+";
	    s << coefs_[i] << "*t";
	    if (i != 1) {
	      s << "^" << i;
	    }
	  }
	}
      }
    }

    out << s.str();
  }


  /*
    typedef typename Coefficients::const_iterator Coefficients_iterator;
    Coefficients_iterator coefficients_begin() const {
    return coefs_.begin();

    }

    Coefficients_iterator coefficients_end() const {
    return coefs_.end();
    }*/

protected:

  std::size_t size() const {return coefs_.size();}

  static const NT & zero_coef() {
    static const NT z(0);
    return z;
  }


  //! The actual coefficients
  Coefficients coefs_;
};

template < class T, class NT,  class C, class Tr>
inline std::ostream &operator<<(std::basic_ostream<C, Tr> &out,
				const Polynomial_impl<T, NT> &poly)
{
  poly.write(out);
  return out;
}


template < class T, class NT, class C, class Tr>
inline std::istream &operator>>(std::basic_istream<C, Tr> &in,
				Polynomial_impl<T, NT> &poly)
{
  poly.read(in);
  //std::cout << "Read " << poly << std::endl;
  return in;
}


//! multiply by a constant
template <class T, class NT>
inline T operator*(const NT &a, const Polynomial_impl<T, NT> &poly)
{
  return (poly * a);
}


//! add to a constant
template <class T, class NT>
inline T operator+(const NT &a, const Polynomial_impl<T, NT> &poly)
{
  return (poly + a);
}


//! add to a constant
template < class T, class NT>
inline T  operator+(int a, const Polynomial_impl<T, NT> &poly)
{
  return (poly + NT(a));
}


//! subtract from a constant
template <class T, class NT>
inline T operator-(const NT &a, const Polynomial_impl<T, NT> &poly)
{
  return -(poly - a);
}


//! subtract from a constant
template <class T, class NT>
inline T operator-(int a, const Polynomial_impl<T, NT> &poly)
{
  return -(poly - NT(a));
}


#undef CGAL_EXCESSIVE

/*template <class NT>
class Input_rep<Polynomial<NT> , CGAL::Maple_format_tag > {
  Input_rep(Polynomial<NT>& p): p_(p){}
  std::istream & operator()(std::istream &in) const {
    
  }
  };*/

} } } //namespace CGAL::POLYNOMIAL::internal
#endif
