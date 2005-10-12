// Copyright (c) 2005  Stanford University (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_POLYNOMIAL_INTERVAL_POLYNOMIAL_H
#define CGAL_POLYNOMIAL_INTERVAL_POLYNOMIAL_H
#include <CGAL/Polynomial/basic.h>
#include <CGAL/Polynomial/internal/interval_arithmetic.h>
#include <CGAL/Polynomial/internal/Polynomial_impl.h>


CGAL_POLYNOMIAL_BEGIN_NAMESPACE

//! A polynomial specialized for interval number types
/*!
  This does not strip leading 0s. 
*/
class Interval_polynomial: public internal::Polynomial_impl<Interval_polynomial, Interval_nt> {
  typedef internal::Polynomial_impl<Interval_polynomial, Interval_nt>  Parent;
  friend class internal::Polynomial_impl<Interval_polynomial, Interval_nt>;

  void finalize(){}
  Interval_polynomial(const Parent &p): Parent(p){}
public:
  Interval_polynomial(){}
  Interval_polynomial(const Interval_nt &nt): Parent(nt){}
  template <class It>
  Interval_polynomial(It b, It e): Parent(b,e){}

  void write(std::ostream &out) const {
    if (Parent::coefs_.size()==0){
      out << "0";
    } else {
      for (unsigned int i=0; i< Parent::coefs_.size(); ++i){
	if (i==0) {
	  out << Parent::coefs_[i];
	} else {
	  if ( Parent::coefs_[i].inf() != 0 ||
	       Parent::coefs_[i].sup() != 0 ) {
	    out << "+";
	    out << Parent::coefs_[i] << "*t";
	    if (i != 1) {
	      out << "^" << i;
	    }
	  }
	}
      }
    }
  }

  void set_coef(unsigned int c, Parent::NT v) {
    CGAL_Polynomial_assertion(c < Parent::coefs_.size());
    coefs_[c]=v;
  }

};

std::ostream &operator<<(std::ostream &out, const Interval_polynomial& ip){
  ip.write(out);
  return out;
}

CGAL_POLYNOMIAL_END_NAMESPACE
#endif
