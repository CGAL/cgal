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

#ifndef CGAL_RATIONAL_ARC_WITH_VER_D_1_H
#define CGAL_RATIONAL_ARC_WITH_VER_D_1_H

#include <CGAL/Arr_enums.h>
#include <CGAL/tags.h>
#include <CGAL/Arr_tags.h>

#include "boost/variant.hpp"
#include <CGAL/Arr_ver_support/Visitors_d_1.h>
//#include "Visitors_d_1.h"

namespace CGAL {
namespace Arr_vertical_rational_arc {

template <class Arc_,
          class Kernel_, 
          class Algebraic_kernel_>
class Rational_arc_with_ver_d_1
{
public:

  typedef Arc_                                                  Arc;
  typedef Kernel_                                               Kernel;
  typedef Algebraic_kernel_                                     Algebraic_kernel;
  
  typedef Rational_arc_with_ver_d_1<Arc,Kernel,Algebraic_kernel>   Self;
   
  typedef typename Arc::Vertical_segment_d_1                  Vertical_segment; 

  typedef typename boost::variant<Arc,Vertical_segment>       Rational_arc_with_ver_d_1_variant;
  
  typedef typename Arc::Algebraic_point_2           Algebraic_point_2;
  typedef typename Arc::Multiplicity                Multiplicity;
  typedef typename Arc::Polynomial_1                Polynomial;
  typedef typename Arc::Coefficient                 Coefficient;
  typedef typename Arc::Arithmetic_kernel           Arithmetic_kernel;
  typedef typename Arc::Rational                    Rational; 
  typedef typename Arc::Integer                     Integer;
  typedef typename Arc::Algebraic_real_1            Algebraic_real_1;
  typedef typename Arc::Algebraic_vector            Algebraic_vector;
  typedef typename Arc::Rat_vector                  Rat_vector;
  typedef typename Arc::Polynomial_1                Polynomial_1;
  typedef typename Arc::Cache                       Cache;
private:
  typedef Arr_vertical_rational_arc::Source_infinite_in_x_visitor< Arc,Vertical_segment>         Source_infinite_in_x_visitor;
  typedef Arr_vertical_rational_arc::Source_infinite_in_y_visitor< Arc,Vertical_segment>         Source_infinite_in_y_visitor;
  typedef Arr_vertical_rational_arc::Target_infinite_in_x_visitor< Arc,Vertical_segment>         Target_infinite_in_x_visitor;
  typedef Arr_vertical_rational_arc::Target_infinite_in_y_visitor< Arc,Vertical_segment>         Target_infinite_in_y_visitor;
  typedef Arr_vertical_rational_arc::Source_visitor< Arc,Vertical_segment,Algebraic_point_2>     Source_visitor;
  typedef Arr_vertical_rational_arc::Source_x_visitor< Arc,Vertical_segment,Algebraic_real_1>    Source_x_visitor;
  typedef Arr_vertical_rational_arc::Target_visitor< Arc,Vertical_segment,Algebraic_point_2>     Target_visitor;
  typedef Arr_vertical_rational_arc::Target_x_visitor< Arc,Vertical_segment,Algebraic_real_1>    Target_x_visitor;
  typedef Arr_vertical_rational_arc::Left_infinite_in_x_visitor< Arc,Vertical_segment>           Left_infinite_in_x_visitor;
  typedef Arr_vertical_rational_arc::Left_infinite_in_y_visitor< Arc,Vertical_segment>           Left_infinite_in_y_visitor;
  typedef Arr_vertical_rational_arc::Right_infinite_in_x_visitor<Arc,Vertical_segment>           Right_infinite_in_x_visitor;
  typedef Arr_vertical_rational_arc::Right_infinite_in_y_visitor<Arc,Vertical_segment>           Right_infinite_in_y_visitor;
  typedef Arr_vertical_rational_arc::Left_visitor< Arc,Vertical_segment,Algebraic_point_2>       Left_visitor;
  typedef Arr_vertical_rational_arc::Left_x_visitor< Arc,Vertical_segment,Algebraic_real_1>      Left_x_visitor;
  typedef Arr_vertical_rational_arc::Right_visitor< Arc,Vertical_segment,Algebraic_point_2>      Right_visitor;
  typedef Arr_vertical_rational_arc::Right_x_visitor< Arc,Vertical_segment,Algebraic_real_1>     Right_x_visitor;
  typedef Arr_vertical_rational_arc::Is_continuous_visitor< Arc,Vertical_segment>                Is_continuous_visitor;
  typedef Arr_vertical_rational_arc::Is_directed_right_visitor< Arc,Vertical_segment>            Is_directed_right_visitor;
private:
  //------------------------------
  //Vertical_Arc members
  //------------------------------
 
  Rational_arc_with_ver_d_1_variant _arc;    // The arc
  enum
  {
    RATIONAL_ARC_TYPE     = 0,
    VERTICAL_TYPE         = 1,
  };

public:
  const Rational_arc_with_ver_d_1_variant & variant() const
  {
    return _arc;
  }
  Rational_arc_with_ver_d_1_variant & variant() 
  {
    return _arc;
  }
  
public:
  //------------
  //Constructors
  //------------

  //---------------------------------------------------------------------------
  //default constructor
  Rational_arc_with_ver_d_1 () {}

 
  //---------------------------------------------------------------------------
  //Constructor of a whole polynomial curve defined by pcoeffs - the rational coefficients of the polynomial p(x).
  Rational_arc_with_ver_d_1 (const Rat_vector& pcoeffs,Cache& cache)
    :_arc (Arc(pcoeffs,cache)) {}
  Rational_arc_with_ver_d_1 ( const Rat_vector& pcoeffs,const Algebraic_real_1& x_s,
                              bool dir_right,Cache& cache)
    :_arc(Arc(pcoeffs,x_s,dir_right,cache)) {}
  Rational_arc_with_ver_d_1 ( const Rat_vector& pcoeffs,const Algebraic_real_1& x_s,
                              const Algebraic_real_1& x_t,Cache& cache)
    :_arc(Arc(pcoeffs,x_s,x_t,cache)) {}
  Rational_arc_with_ver_d_1 (const Rat_vector& pcoeffs, const Rat_vector& qcoeffs,Cache& cache)
    :_arc(Arc(pcoeffs,qcoeffs,cache)) {}
  Rational_arc_with_ver_d_1 ( const Rat_vector& pcoeffs, const Rat_vector& qcoeffs,
                              const Algebraic_real_1& x_s, bool dir_right,Cache& cache) 
    :_arc(Arc(pcoeffs,qcoeffs,x_s,dir_right,cache)) {}
  Rational_arc_with_ver_d_1 ( const Rat_vector& pcoeffs, const Rat_vector& qcoeffs,
                              const Algebraic_real_1& x_s, const Algebraic_real_1& x_t,Cache& cache)
    :_arc(Arc(pcoeffs,qcoeffs,x_s,x_t,cache)) {}

  Rational_arc_with_ver_d_1 (const Algebraic_point_2& p1,const Algebraic_point_2& p2) 
    :_arc (Vertical_segment(p1,p2,cache)) {}
  Rational_arc_with_ver_d_1 (const Algebraic_point_2& p,bool is_directed_up) 
    :_arc (Vertical_segment(p,is_directed_up)) {}
  Rational_arc_with_ver_d_1 (const Algebraic_point_2& p)
    :_arc (Vertical_segment(p)) {}
  Rational_arc_with_ver_d_1 (const Vertical_segment& other)
    :_arc (other) {}
  Rational_arc_with_ver_d_1 (const Arc& other) 
    :_arc (other) {}

  const Polynomial_1& numerator () const
  {
    CGAL_precondition (_arc.which() == RATIONAL_ARC_TYPE);
    if ( const Arc* pc = boost::get<Arc>( &_arc ) )
      return pc->numerator ();
    else
      CGAL_postcondition(false);
    //should not be reached
    assert (false);
    return dummy_p;
  }
  const Polynomial_1& denominator () const
  {
    CGAL_precondition (_arc.which() == RATIONAL_ARC_TYPE);
    if ( const Arc* pc = boost::get<Arc>( &_arc ) )
      return pc->denominator();
    else
      CGAL_postcondition(false);
    //should not be reached
    assert (false);
    return dummy_p;
  }

  //---------------------------------------------------------
  //Check if the x-coordinate of the source point is infinite
  Arr_parameter_space source_infinite_in_x () const
  {
    return (boost::apply_visitor(Source_infinite_in_x_visitor(),_arc));
  }

  //---------------------------------------------------------
  //Check if the y-coordinate of the source point is infinite
  Arr_parameter_space source_infinite_in_y () const
  {
    return (boost::apply_visitor(Source_infinite_in_y_visitor(),_arc));
  }

  //---------------------------------------------------------
  //Check if the x-coordinate of the target point is infinite
  Arr_parameter_space target_infinite_in_x () const
  {
    return (boost::apply_visitor(Target_infinite_in_x_visitor(),_arc));
  }

  //---------------------------------------------------------
  //Check if the y-coordinate of the target point is infinite
  Arr_parameter_space target_infinite_in_y () const
  {
    return (boost::apply_visitor(Target_infinite_in_y_visitor(),_arc));  
  }

  //--------------------
  //Get the source point
  const Algebraic_point_2& source () const
  {
    return (boost::apply_visitor(Source_visitor(),_arc));  
  }

  //----------------------------------------
  //Get the x-coordinate of the source point
  Algebraic_real_1 source_x () const
  {
    return (boost::apply_visitor(Source_x_visitor(),_arc));  
  }

  //--------------------
  //Get the target point
  const Algebraic_point_2& target () const
  {
    return (boost::apply_visitor(Target_visitor(),_arc));  
  }

  //----------------------------------------
  //Get the x-coordinate of the target point
  Algebraic_real_1 target_x () const
  {
    return (boost::apply_visitor(Target_x_visitor(),_arc));  
  }

  //-------------------------------------------------------
  //Check if the x-coordinate of the left point is infinite
  Arr_parameter_space left_infinite_in_x () const
  {
    return (boost::apply_visitor(Left_infinite_in_x_visitor(),_arc));
  }

  //-------------------------------------------------------
  //Check if the y-coordinate of the left point is infinite
  Arr_parameter_space left_infinite_in_y () const
  {
    return (boost::apply_visitor(Left_infinite_in_y_visitor(),_arc));
  }
  //--------------------------------------------------------
  //Check if the x-coordinate of the right point is infinite
  Arr_parameter_space right_infinite_in_x () const
  {
    return (boost::apply_visitor(Right_infinite_in_x_visitor(),_arc));
  }

  //--------------------------------------------------------
  //Check if the y-coordinate of the right point is infinite
  Arr_parameter_space right_infinite_in_y () const
  {
    return (boost::apply_visitor(Right_infinite_in_y_visitor(),_arc));
  }

  //------------------------------------
  //Get the x_value of the left endpoint
  const Algebraic_real_1 left_x () const
  {
     return (boost::apply_visitor(Left_x_visitor(),_arc));
  }

  //------------------------------------
  //Get the x_value of the right endpoint
  const Algebraic_real_1 right_x () const
  {
    return (boost::apply_visitor(Right_x_visitor(),_arc));
  }
  //---------------------
  //Get the left endpoint
  const Algebraic_point_2 left () const
  {
    return (boost::apply_visitor(Left_visitor(),_arc));
  }

  //----------------------
  //Get the right endpoint
  const Algebraic_point_2 right () const
  {
    return (boost::apply_visitor(Right_visitor(),_arc));
  }

  //------------------------------
  //Check if the arc is continuous
  bool is_continuous () const
  {
    return (boost::apply_visitor(Is_continuous_visitor(),_arc));
  }

  //----------------------------------
  //Check if the arc is directed right
  bool is_directed_right () const
  {
    return (boost::apply_visitor(Is_directed_right_visitor(),_arc));
  }
private:
  Polynomial_1 dummy_p;

};

//-------------------------------
//! Exporter for rational arcs.
template <class Arc_,
          class Kernel_, 
          class Algebraic_kernel_>
std::ostream&
operator<< (std::ostream& os, 
            const Rational_arc_with_ver_d_1<Arc_, Kernel_,  Algebraic_kernel_> & _arc)
{
   return (boost::apply_visitor(Print_visitor < Arc_, typename Arc_::Vertical_segment_d_1 > (os),_arc.variant()));
}

}   //name_space Arr_vertical_rational_arc
}       //namespace CGAL {
#endif //CGAL_RATIONAL_ARC_WITH_VER_D_1_H
