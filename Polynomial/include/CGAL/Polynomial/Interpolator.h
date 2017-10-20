// Copyright (c) 2008 Max-Planck-Institute Saarbruecken (Germany)
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
// Author(s)     : Michael Hemmer <hemmer@informatik.uni-mainz.de> 

#ifndef CGAL_POLYNOMIAL_INTERPOLATE_H
#define CGAL_POLYNOMIAL_INTERPOLATE_H

namespace CGAL {
namespace internal {

// Class for interpolation of univariate or multivariate polynomials. 
// The template argument must be a model of concept Polynomial_d
// 
//
template <class Polynomial_d_>
class Interpolator{
  typedef CGAL::Polynomial_traits_d<Polynomial_d_> PT;
    
public:
  typedef typename PT::Polynomial_d Polynomial_d; 
  typedef typename PT::Coefficient_type Coeff; 
  typedef typename PT::Innermost_coefficient_type IC;

private: 
  typedef typename CGAL::Coercion_traits<Coeff,IC>::Cast IC2Coeff;
  typedef typename PT::Construct_polynomial Construct;

  std::vector<IC> xvals; 
  std::vector<Coeff> yvals; 
  std::vector<Coeff> b; 
    
  bool valid; 
  Polynomial_d interpolant; 


  Coeff eval_newton(int n, IC z)
  {
    Coeff p(b[n]);
    for (int i = n-1; i >=0; i--){
      Coeff tmp(IC2Coeff()((z - xvals[i])));
      p = p * tmp + b[i];
    }
    return p;
  }


  Polynomial_d eval_newton_poly(int n)
  {    
    CGAL_precondition(n >=0);
    Polynomial_d p(Construct()(b[n]));
    for (int i = n-1; i >=0; i--) {
      Polynomial_d tmp = Construct()(IC2Coeff()(-xvals[i]),Coeff(1));
      p = (p * tmp) + b[i];
    }
    return p;
  }
    
public:
  Interpolator(){};
  
  // constructor from an InputIterator range with value type std::pair<IC,Coeff> 
  template<class InputIterator>
  Interpolator(InputIterator begin, InputIterator end){
    for(InputIterator it = begin; it != end; it++){
      add_interpolation_point(*it);
    }
  }
    
  /*
    Interpolator(std::vector<IC> xvals_, std::vector<Coeff> yvals_)
    : valid(false) {
    CGAL_precondition(xvals_.size() != 0);
    CGAL_precondition(xvals_.size() == yvals_.size());
    for(unsigned i = 0; i < xvals_.size();  i++){
    add_interpolation_point(xvals_[i],yvals_[i]);
    }
    }
  */
    
  // void add_interpolation_point(std::pair<IC,Coeff> point){
  //    add_interpolation_point(point.first, point.second);
  // }
    
  // void add_interpolation_point(IC xval, Coeff yval){
  void add_interpolation_point(std::pair<IC,Coeff> point){
    valid = false;
//        CGAL_precondition(0 == std::count(xval, xvals.begin(), yvals.end()));
    xvals.push_back(point.first);
    yvals.push_back(point.second);
        
    Coeff num, den;
    int k = static_cast<int>(xvals.size()) - 1; 
    if(k == 0){
      b.push_back(yvals[0]);
    }else{
      num = yvals[k] - eval_newton(k-1,xvals[k]);    
      den = Coeff(1);            
      for (int j = 0; j < k; j++) {
        // (k-j) if xvals's are sequential
        den *= (xvals[k] - xvals[j]);
      }
      b.push_back(num / den);
    }
  }
    
  Polynomial_d get_interpolant(){ 
    if (xvals.size() == 0) return Polynomial_d(0);
    // TODO: compute new interpolant from old interpolant ?
    if(!valid)
      interpolant = eval_newton_poly(static_cast<int>(xvals.size())-1);
    return interpolant; 
  }
    
};

} // namespace internal
} //namespace CGAL

#endif // CGAL_POLYNOMIAL_INTERPOLATE_H 
