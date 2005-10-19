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
    return (CGAL::make_object (this->arrangement()->unbounded_face()));
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

      //check isolated vertices
      Isolated_vertices_const_iterator   iso_verts_it;
      for (iso_verts_it = fh->isolated_vertices_begin();
          iso_verts_it != fh->isolated_vertices_end(); ++iso_verts_it)
      {
        if (traits->equal_2_object()(p, iso_verts_it->point()))
        {
          Vertex_const_handle  vh = iso_verts_it;
          return (CGAL::make_object (vh));
        }
      }

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
  Halfedge_const_handle invalid_he;
 
  X_curve_plus cv = td.vertical_ray_shoot(p, td_lt, shoot_up);

  // treat special case, where trapezoid is unbounded.
  //	for then get_parent() is not defined
  if (td_lt==TD::UNBOUNDED_TRAPEZOID)
  { 
    return (_check_isolated_for_vertical_ray_shoot(invalid_he, p, shoot_up));
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
    if ((shoot_up && h->direction() == SMALLER) ||
        (!shoot_up && h->direction() == LARGER))
      h=h->twin();

    return (CGAL::make_object(h));

  case TD::TRAPEZOID:
    if (!(((traits->is_in_x_range_2_object()(h->curve(),p)) &&
          (traits->compare_y_at_x_2_object()(p, h->curve()) == LARGER)) ==
          (traits->compare_x_2_object()(h->source()->point(),
                                        h->target()->point()) == SMALLER)
        ))
        h = h->twin();

    return (_check_isolated_for_vertical_ray_shoot(h, p, shoot_up));

  default:
    CGAL_assertion(false);
    break;
  }

  return (_check_isolated_for_vertical_ray_shoot(invalid_he, p, shoot_up));
}

//-----------------------------------------------------------------------------
// In vertical ray shoot, when the closest halfedge is found (or unbounded face)
// we check the isolated vertices inside the face to check whether there
// is an isolated vertex right above/below the query point.
// 
template <class Arrangement>
Object Arr_trapezoid_ric_point_location<Arrangement>
::_check_isolated_for_vertical_ray_shoot(Halfedge_const_handle &halfedge_found, 
                                         const Point_2& p, 
                                         bool shoot_up) const
{
  const Comparison_result point_above_under = (shoot_up ? SMALLER : LARGER);

  Isolated_vertices_const_iterator   iso_verts_it;
  Vertex_const_handle                closest_iso_v;
  const Vertex_const_handle          invalid_v;
  const Halfedge_const_handle        invalid_he;
  Face_const_handle                  face;

  if (halfedge_found == invalid_he)
    face = this->arrangement()->unbounded_face();
  else
    face = halfedge_found->face();

  for (iso_verts_it = face->isolated_vertices_begin();
       iso_verts_it != face->isolated_vertices_end(); ++iso_verts_it)
  {
    // The current isolated vertex should have the same x-coordinate as the
    // query point in order to be below or above it.
    if (compare_x (p, iso_verts_it->point()) != EQUAL)
      continue;

    // Make sure the isolated vertex is above the query point (if we shoot up)
    // or below it (if we shoot down).
    if (compare_xy (p, iso_verts_it->point()) != point_above_under)
      continue;

    // Check if the current isolated vertex lies closer to the query point than
    // the closest feature so far.
    if (closest_iso_v == invalid_v)
    {
      // Compare the current isolated vertex with the closest halfedge.
      if (closest_he == invalid_he ||
          compare_y_at_x (iso_verts_it->point(),
                          closest_he->curve()) == point_above_under)
      {
        closest_iso_v = iso_verts_it;
      }
    }
    else if (compare_xy (iso_verts_it->point(),
                         closest_iso_v->point()) == point_above_under)
    {
      closest_iso_v = iso_verts_it;
    }
  }

  if (closest_iso_v != invalid_v)
  {
    // The first object we encounter when we shoot a vertical ray from p is
    // an isolated vertex:
    return (CGAL::make_object (closest_iso_v));
  }

  if (halfedge_found == invalid_he)
    return (CGAL::make_object (face));

  //else
  return (CGAL::make_object(h));
}



CGAL_END_NAMESPACE

#endif
