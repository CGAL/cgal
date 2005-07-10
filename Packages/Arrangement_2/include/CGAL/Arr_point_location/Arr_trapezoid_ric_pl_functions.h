// Copyright (c) 2005  Tel-Aviv University (Israel).
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
// Author(s)     : Idit Haran   <haranidi@post.tau.ac.il>
#ifndef CGAL_ARR_TRAPEZOID_RIC_POINT_LOCATION_FUNCTIONS_H
#define CGAL_ARR_TRAPEZOID_RIC_POINT_LOCATION_FUNCTIONS_H

/*! \file
* Member-function definitions for the Arr_trapezoid_ric_point_location<Arrangement>
* class.
*/

#define CGAL_TRAP_DEBUG

#ifdef CGAL_TRG_DEBUG
	#define TRAP_PRINT_DEBUG(expr)   std::cout << expr << std::endl
#else
	#define TRAP_PRINT_DEBUG(expr)
#endif

CGAL_BEGIN_NAMESPACE

//-----------------------------------------------------------------------------
// Locate the arrangement feature containing the given point.
//
template <class Arrangement_2>
Object Arr_trapezoid_ric_point_location<Arrangement_2>
::locate (const Point_2& p) const
{
  TRAP_PRINT_DEBUG("locate point "<<p);

  //there are different internal compiler errors if we
  // typedef the Locate_type
  typename TD::Locate_type td_lt; 

  const X_curve_plus& cv = td.locate(p,td_lt).top();

  TRAP_PRINT_DEBUG("after td.locate");

  // treat special case, where trapezoid is unbounded.
  //	for then get_parent() is not defined
  if (td_lt==TD::UNBOUNDED_TRAPEZOID)
  {
    TRAP_PRINT_DEBUG("UNBOUNDED_TRAPEZOID");
    return (CGAL::make_object (p_arr->unbounded_face()));
  }

  Halfedge_const_handle h = cv.get_parent();

  switch(td_lt)
  {
  case TD::POINT:
    {
      TRAP_PRINT_DEBUG("POINT");
      if (traits->equal_2_object()(h->target()->point(), p))
      {
        Vertex_const_handle vh = h->target();
        return (CGAL::make_object (vh));
      }
      if (traits->equal_2_object()(h->source()->point(), p))
      {
        Vertex_const_handle vh = h->source();
        return (CGAL::make_object (vh));
      }
      else
        CGAL_assertion(false);
      break;
    }

  case TD::CURVE:
    {
      TRAP_PRINT_DEBUG("CURVE");
      if ( traits->is_in_x_range_2_object()(cv,p) && 
           traits->compare_y_at_x_2_object()(p,cv) == EQUAL)
        return (CGAL::make_object(h));
      else
      {
        //ixx
        std::cerr << "curve is: "<<cv<<" point is: "<<p <<std::endl; 
        CGAL_assertion(false);
      }
      break;
    }

  case TD::TRAPEZOID:
    {
      TRAP_PRINT_DEBUG("TRAPEZOID");
      if (!(((traits->is_in_x_range_2_object()(h->curve(),p)) &&
        (traits->compare_y_at_x_2_object()(p, h->curve()) == LARGER)) ==
        (traits->compare_x_2_object()(h->source()->point(),
        h->target()->point()) == SMALLER)
        ))
        h = h->twin();
      Face_const_handle fh = h->face();
      return (CGAL::make_object(fh));
    }
  default:
    TRAP_PRINT_DEBUG("DEFAULT");
    CGAL_assertion(false);
    break;
  }

  TRAP_PRINT_DEBUG("EMPTY");
  return Object();   
}

//-----------------------------------------------------------------------------
// Locate the arrangement feature which a vertical ray emanating from the
// given point hits, considering isolated vertices.
//
template <class Arrangement>
Object Arr_trapezoid_ric_point_location<Arrangement>
::_vertical_ray_shoot (const Point_2& p, bool shoot_up) const
{
  //trying to workaround internal compiler error
  typename TD::Locate_type td_lt;
  
  X_curve_plus cv = td.vertical_ray_shoot(p, td_lt, shoot_up);

  // treat special case, where trapezoid is unbounded.
  //	for then get_parent() is not defined
  if (td_lt==TD::UNBOUNDED_TRAPEZOID)
  {
    return (CGAL::make_object (p_arr->unbounded_face()));
  }

  Halfedge_const_handle h=cv.get_parent();

  switch(td_lt)
  {
  case TD::POINT:
    if (traits->equal_2_object()(h->target()->point(), p))
    {
      Vertex_const_handle vh = h->target();
      return (CGAL::make_object (vh));
    }
    if (traits->equal_2_object()(h->source()->point(), p))
    {
      Vertex_const_handle vh = h->source();
      return (CGAL::make_object (vh));
    }
    else
      CGAL_assertion(false);
    break;

  case TD::CURVE:
    /* special case:
    h->source()->point(),p,h->target()->point() have same x
    coardinates.
    return value should be h(no h->twin()).
    */
    // orientation of h
    if (shoot_up == traits->compare_x_2_object()(h->source()->point(),
                                                 h->target()->point()))
      h=h->twin();
    return (CGAL::make_object(h));

  case TD::TRAPEZOID:
    if (!(((traits->is_in_x_range_2_object()(h->curve(),p)) &&
          (traits->compare_y_at_x_2_object()(p, h->curve()) == LARGER)) ==
          (traits->compare_x_2_object()(h->source()->point(),
                                        h->target()->point()) == SMALLER)
        ))
        h = h->twin();
    return (CGAL::make_object(h));

  default:
    CGAL_assertion(false);
    break;
  }

  return Object();  
}


CGAL_END_NAMESPACE

#endif
