// Copyright (c) 2006-2008 Inria Lorraine (France). All rights reserved.
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
// Author(s)     :  Michael Hemmer <hemmer@mpi-inf.mpg.de> 
//
// ============================================================================


namespace CGAL {

template<>
class Real_embeddable_traits< Algebraic_1  >
  : public INTERN_RET::Real_embeddable_traits_base<  Algebraic_1  , CGAL::Tag_true > {

public:
  typedef INTERN_RET::Real_embeddable_traits_base<  Algebraic_1  , CGAL::Tag_true > Base;
 
  typedef CGAL::Tag_true                 Is_real_embeddable;
  typedef bool                           Boolean;
  typedef CGAL::Sign                     Sign;
  typedef CGAL::Comparison_result        Comparison_result;


  typedef Algebraic_1  Type;
  typedef Base::Compare Compare; // todo: get a more efficient impl
   
  class Sgn
    : public std::unary_function< Type, CGAL::Sign > {
  public:
    CGAL::Sign operator()( const Type& a ) const {
      return Compare()(a, Type(0));
    }
  };

  class To_double
    : public std::unary_function< Type, double > {
  public:
    double operator()(const Type& a) const {
      return a.to_double();
    }
  };

  class To_interval
    : public std::unary_function< Type, std::pair<double, double> > {
  public:
    std::pair<double, double> operator()(const Type& a) const {        
      return a.to_interval();
    }
  };       

  class Is_zero
    : public std::unary_function< Type, Boolean> {
  public:
    bool operator()(const Type& a) const {        
      return  Sgn()(a) == CGAL::ZERO;
    }
  };
 
  class Is_finite
    :public std::unary_function< Type, Boolean> {
  public:
    bool operator()(const Type& ) const {
      return  true;
    }
  };

  class Abs
    :public  std::unary_function< Type, Type> {
  public:
    Type operator()(const Type& a) const {
      return Sgn()(a)==CGAL::NEGATIVE?-a:a;
    }
  };
};
} // namespace CGAL
