// Copyright (c) 2010  Tel-Aviv University (Israel), INRIA Sophia-Antipolis (France).
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
// $URL$
// $Id$
//
// Author(s)     : Oren Salzman <orenzalz@post.tau.ac.il >
//                 Michael Hemmer <Michael.Hemmer@sophia.inria.fr>

#ifndef CGAL_VERTICAL_SEGMENT_D_1
#define CGAL_VERTICAL_SEGMENT_D_1

#include <CGAL/tags.h>
#include <CGAL/Arr_tags.h>
#include <CGAL/assertions.h>
#include <CGAL/Object.h>

#include <CGAL/Arr_rat_arc/Base_rational_arc_ds_1.h>
#include <CGAL/Arr_rat_arc/Cache.h>
#include <CGAL/Arr_rat_arc/Rational_function.h>
#include <CGAL/Arr_rat_arc/Rational_function_canonicalized_pair.h>
#include <CGAL/Arr_rat_arc/Algebraic_point_2.h>

namespace CGAL {
namespace Arr_rational_arc {

template <  class Kernel_, 
            class Algebraic_kernel_ = Algebraic_kernel_d_1 <typename Fraction_traits <typename Kernel_::FT>::Numerator_type> >
class Vertical_segment_d_1: public Base_rational_arc_ds_1<Kernel_, Algebraic_kernel_>
{
public:
  typedef Kernel_                                           Kernel;
  typedef Algebraic_kernel_                                 Algebraic_kernel;
  typedef Base_rational_arc_ds_1<Kernel_, Algebraic_kernel> Base;
  typedef Vertical_segment_d_1<Kernel_, Algebraic_kernel>   Self;

  typedef CGAL::Arr_rational_arc::Rational_function <Kernel_, Algebraic_kernel>     Rational_function;
  typedef CGAL::Arr_rational_arc::Rational_function_pair <Kernel_,Algebraic_kernel> Rational_function_pair;
  typedef CGAL::Arr_rational_arc::Algebraic_point_2 <Kernel_,Algebraic_kernel>      Algebraic_point_2;
  typedef typename Base::Algebraic_real_1                                           Algebraic_real_1;

  typedef CGAL::Arr_rational_arc::Cache<Kernel_,Algebraic_kernel_>                  Cache;
 
public:
  ////////////////////
  //concept functions
  ////////////////////
  Vertical_segment_d_1() {}
  //construct a vertical line at p.x
  Vertical_segment_d_1(const Algebraic_point_2& p)
    :_max(p) , _min(p), _max_parameter_space(ARR_TOP_BOUNDARY), _min_parameter_space(ARR_BOTTOM_BOUNDARY)
  {}
  //construct a vertical ray at p
  Vertical_segment_d_1(const Algebraic_point_2& p, bool is_directed_up)
    :_max(p) , _min(p)
  {
    _max_parameter_space = is_directed_up ?  ARR_TOP_BOUNDARY : ARR_INTERIOR ;
    _min_parameter_space = is_directed_up ?  ARR_INTERIOR  : ARR_BOTTOM_BOUNDARY ;
  }
  //construct a vertical line at between points p1,p2
  Vertical_segment_d_1(const Algebraic_point_2& p1,const Algebraic_point_2& p2,Cache& cache)
    :_max_parameter_space(ARR_INTERIOR), _min_parameter_space(ARR_INTERIOR)
  {
    CGAL_precondition(p1.x() == p2.x());
    CGAL_precondition(p1.compare_xy_2(p2,cache) != CGAL::EQUAL);

    Rational_function_pair rat_func_pair = cache.get_rational_pair( p1.rational_function(),
                                                                    p2.rational_function());
    Comparison_result comp (rat_func_pair.compare_f_g_at(p1.x()));
    CGAL_postcondition(comp  != CGAL::EQUAL);

    _min  = (comp == CGAL::SMALLER) ? p1 : p2;
    _max  = (comp == CGAL::SMALLER) ? p2 : p1;
  }

  Vertical_segment_d_1(const Self& other)
  {
    _max = other.max();
    _min = other.min();
    _max_parameter_space = other.max_parameter_space();
    _min_parameter_space = other.min_parameter_space();
  }
  bool has_same_x(const Self& other) const
  {
    return (_min.x() == other.x()) ;
  }
  bool operator== (const Self& other) const
  {
    if ((has_same_min(other))&&
        (has_same_max(other)))
      return true;
    return false;
  }
  Self& operator= (const Self & other)
  {
    if (this != &other) // protect against invalid self-assignment
    {
        _min = other.min();
        _max = other.max();
        _min_parameter_space = other.min_parameter_space();
        _max_parameter_space = other.max_parameter_space();
    }
    return *this;
  }
  Algebraic_point_2& max() 
  {
    return _max;
  }
  const Algebraic_point_2 & max() const
  {
    return _max;
  }
  Algebraic_point_2& min() 
  {
      return _min;
  }
  const Algebraic_point_2 & min() const
  {
    return _min;
  }

  CGAL::Arr_parameter_space max_parameter_space() const
  {
    return _max_parameter_space;
  }
  CGAL::Arr_parameter_space min_parameter_space() const
  {
    return _min_parameter_space;
  }
  
  ////////////////////////
  //non concept functions
  ////////////////////////
public:
  const Algebraic_real_1& x() const 
  {
    return _min.x();
  }
  const Rational_function& min_f() const 
  {
    CGAL_precondition (_min_parameter_space == CGAL::ARR_INTERIOR );
    return _min.rational_function();
  }
  const Rational_function& max_f() const 
  {
    CGAL_precondition (_max_parameter_space == CGAL::ARR_INTERIOR );
    return _max.rational_function();
  }
  bool max_bounded() const
  {
     CGAL_precondition( (_max_parameter_space == ARR_TOP_BOUNDARY) ||
                        (_max_parameter_space == ARR_INTERIOR)) ;
     return (_max_parameter_space == ARR_TOP_BOUNDARY) ? false : true;
  }
  bool min_bounded() const
  {
   CGAL_precondition( (_min_parameter_space == ARR_BOTTOM_BOUNDARY) ||
                      (_min_parameter_space == ARR_INTERIOR)) ;
    return (_min_parameter_space == ARR_BOTTOM_BOUNDARY) ? false : true;
  }
  
private:

  bool has_same_min(const Self& other) const
  {
    if (this->_min_parameter_space != other._min_parameter_space)
      return false;
    if (this->_min_parameter_space ==  CGAL:: ARR_BOTTOM_BOUNDARY)
      return true;
    return (_min == other.min());
  }
  bool has_same_max(const Self& other) const
  {
    if (this->_max_parameter_space != other._max_parameter_space)
      return false;
    if (this->_min_parameter_space ==  CGAL::ARR_TOP_BOUNDARY)
      return true;
    return (_max == other.max());
  }
  bool is_line() const
  {
    return ((_max_parameter_space == ARR_TOP_BOUNDARY   ) && 
            (_min_parameter_space == ARR_BOTTOM_BOUNDARY) );
  }
  bool is_ray() const
  {
    CGAL_precondition( (_max_parameter_space == ARR_TOP_BOUNDARY) ||
                        (_max_parameter_space == ARR_INTERIOR)) ;
    CGAL_precondition( (_min_parameter_space == ARR_BOTTOM_BOUNDARY) ||
                      (_min_parameter_space == ARR_INTERIOR)) ;
    return (  ((_max_parameter_space == ARR_TOP_BOUNDARY) && (_min_parameter_space == ARR_INTERIOR        ))||
              ((_max_parameter_space == ARR_INTERIOR    ) && (_min_parameter_space == ARR_BOTTOM_BOUNDARY ))  );
  }
  bool is_segment() const
  {
    CGAL_precondition( (_max_parameter_space == ARR_TOP_BOUNDARY) ||
                        (_max_parameter_space == ARR_INTERIOR)) ;
    CGAL_precondition( (_min_parameter_space == ARR_BOTTOM_BOUNDARY) ||
                      (_min_parameter_space == ARR_INTERIOR)) ;
    return ((_max_parameter_space ==  ARR_INTERIOR) && 
            (_min_parameter_space ==  ARR_INTERIOR) );
  }
  
public:
  //------------------------
  // Print the vertical segment.
  std::ostream& print (std::ostream& os) const
   {    
      os << "x = " << _min.x() << " ";
      //print lower point
      if (_min_parameter_space == CGAL::ARR_BOTTOM_BOUNDARY)
      {
         os << "min y = -oo";
      }
      else
      {
         os << "y = (";    
         print_polynomial (os, _min.rational_function().numer(), 'x');
         os << ") / (";
         print_polynomial (os, _min.rational_function().denom(), 'x');
         os << ")";
      }
     
      os << " to ";
     
      //print lower point
      if (_max_parameter_space == CGAL::ARR_TOP_BOUNDARY)
      {
         os << "max y = +oo";
      }
      else
      {
         os << "y = (";    
         print_polynomial (os, _max.rational_function().numer(), 'x');
         os << ") / (";
         print_polynomial (os, _max.rational_function().denom(), 'x');
         os << ")";
      }
      return (os);
   }

private:
  Algebraic_point_2 _min;
  Algebraic_point_2 _max;
  
  CGAL::Arr_parameter_space _min_parameter_space;
  CGAL::Arr_parameter_space _max_parameter_space;
}; //Vertical_segment

//-------------------------------
//! Exporter for Vertical_segment.
template < class Kernel_,
   class Algebraic_kernel_  >
std::ostream&
operator<< (std::ostream& os, 
    const Vertical_segment_d_1<Kernel_, Algebraic_kernel_> & ver)
{
  return (ver.print (os));
}

} //namespace Arr_rational_arc {
} //namespace CGAL {


#endif  //CGAL_VERTICAL_SEGMENT_D_1
