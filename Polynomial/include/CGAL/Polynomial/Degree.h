// Copyright (c) 2009 Max-Planck-Institute Saarbruecken (Germany).
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
// SPDX-License-Identifier: LGPL-3.0+
// 
//
// Author(s)     : Michael Hemmer <hemmer@mpi-inf.mpg.de> 
//
// ============================================================================


#ifndef CGAL_POLYNOMIAL_DEGREE_H
#define CGAL_POLYNOMIAL_DEGREE_H

#include <CGAL/Exponent_vector.h>
#include <CGAL/Polynomial/misc.h>


namespace CGAL {

namespace internal{

template <typename Polynomial > struct Degree;

// Polynomial musst be at least univariate ! 
template <typename Coeff_ > 
struct Degree<Polynomial<Coeff_> >
  : public CGAL::unary_function< Polynomial<Coeff_> , int  >{
  
private: 
  typedef Polynomial<Coeff_> Polynomial_d;
  
  // for univariate polynomials 
  template <typename Coeff>
  int degree(const Polynomial<Coeff>& p, int i) const {
    (void) i; 
    CGAL_assertion(i == 0);
    return p.degree();
  }
  
  // for multivariate polynomials
  template <typename Coeff>
  int degree(const Polynomial<Polynomial<Coeff> >& p, int i) const {
    typedef Polynomial<Polynomial<Coeff> > Poly;
    
    if(i == (Dimension<Poly>::value-1))
      return p.degree();
    
    int result = 0;
    for(typename Poly::const_iterator it = p.begin(); it != p.end(); it++){
      result = (std::max)(result,degree(*it,i));
    }
    return result; 
  }

public:
  int operator()(
      const Polynomial_d& p, 
      int i = (Dimension<Polynomial_d>::value-1)) const {
    CGAL_assertion(i < Dimension<Polynomial_d>::value);
    CGAL_assertion(i >= 0);
    return this->degree(p,i);
  }     

};

} // namespace internal
} //namespace CGAL

#endif //CGAL_POLYNOMIAL_DEGREE_H
