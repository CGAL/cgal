// Copyright (c) 2005,2006,2007,2008,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
// Author(s)     : Ron Wein          <wein@post.tau.ac.il>
//                 Efi Fogel         <efif@post.tau.ac.il>
//                 Eric Berberich    <ericb@post.tau.ac.il>
//                 (based on old version by: Iddo Hanniel,
//                                           Eyal Flato,
//                                           Oren Nechushtan,
//                                           Ester Ezra,
//                                           Shai Hirsch,
//                                           and Eugene Lipovetsky)
//
#ifndef CGAL_ARRANGEMENT_ON_SURFACE_2_IMPL_H
#define CGAL_ARRANGEMENT_ON_SURFACE_2_IMPL_H

#ifndef CGAL_ARRANGEMENT_ON_SURFACE_INSERT_VERBOSE
#define CGAL_ARRANGEMENT_ON_SURFACE_INSERT_VERBOSE 0
#endif

/*! \file
 * Member-function definitions for the Arrangement_2<GeomTraits, TopTraits> 
 * class-template.
 */

#include <CGAL/function_objects.h> 

#if CGAL_ARRANGEMENT_ON_SURFACE_INSERT_VERBOSE
#include <CGAL/Arr_topology_traits/Sign_of_path.h>
#endif

namespace CGAL {

//-----------------------------------------------------------------------------
// Default constructor.
//
template<class GeomTraits, class TopTraits>
Arrangement_on_surface_2<GeomTraits, TopTraits>::Arrangement_on_surface_2 () :
  m_topol_traits()
{
  
  typedef has_Left_side_category<GeomTraits> Cond_left;
  typedef internal::Validate_left_side_category< GeomTraits, Cond_left::value > 
    Validate_left_side_category;
  void (Validate_left_side_category::*pleft)(void) =
    &Validate_left_side_category::template missing__Left_side_category<int>;
  (void)pleft;
  
 typedef has_Bottom_side_category<GeomTraits> Cond_bottom;
  typedef internal::Validate_bottom_side_category< GeomTraits, Cond_bottom::value > 
    Validate_bottom_side_category;
  void (Validate_bottom_side_category::*pbottom)(void) =
    &Validate_bottom_side_category::template missing__Bottom_side_category<int>;
  (void)pbottom;

 typedef has_Top_side_category<GeomTraits> Cond_top;
  typedef internal::Validate_top_side_category< GeomTraits, Cond_top::value > 
    Validate_top_side_category;
  void (Validate_top_side_category::*ptop)(void) =
    &Validate_top_side_category::template missing__Top_side_category<int>;
  (void)ptop;

 typedef has_Right_side_category<GeomTraits> Cond_right;
  typedef internal::Validate_right_side_category< GeomTraits, Cond_right::value > 
    Validate_right_side_category;
  void (Validate_right_side_category::*pright)(void) =
    &Validate_right_side_category::template missing__Right_side_category<int>;
  (void)pright;

  // Initialize the DCEL structure to represent an empty arrangement.
  m_topol_traits.init_dcel ();

  // Allocate the traits.
  m_geom_traits = new Traits_adaptor_2;
  m_own_traits = true;
}

//-----------------------------------------------------------------------------
// Copy constructor.
//
template<class GeomTraits, class TopTraits>
Arrangement_on_surface_2<GeomTraits, TopTraits>::
Arrangement_on_surface_2(const Self& arr) :
  m_geom_traits (NULL),
  m_own_traits (false)
{
  assign (arr);
}

//-----------------------------------------------------------------------------
// Constructor given a traits object.
//
template<class GeomTraits, class TopTraits>
Arrangement_on_surface_2<GeomTraits, TopTraits>::
Arrangement_on_surface_2(const Geometry_traits_2 * geom_traits) :
  m_topol_traits (geom_traits)
{
  
 typedef has_Left_side_category<GeomTraits> Cond_left;
  typedef internal::Validate_left_side_category< GeomTraits, Cond_left::value > 
    Validate_left_side_category;
  void (Validate_left_side_category::*pleft)(void) =
    &Validate_left_side_category::template missing__Left_side_category<int>;
  (void)pleft;
  
 typedef has_Bottom_side_category<GeomTraits> Cond_bottom;
  typedef internal::Validate_bottom_side_category< GeomTraits, Cond_bottom::value > 
    Validate_bottom_side_category;
  void (Validate_bottom_side_category::*pbottom)(void) =
    &Validate_bottom_side_category::template missing__Bottom_side_category<int>;
  (void)pbottom;

 typedef has_Top_side_category<GeomTraits> Cond_top;
  typedef internal::Validate_top_side_category< GeomTraits, Cond_top::value > 
    Validate_top_side_category;
  void (Validate_top_side_category::*ptop)(void) =
    &Validate_top_side_category::template missing__Top_side_category<int>;
  (void)ptop;

 typedef has_Right_side_category<GeomTraits> Cond_right;
  typedef internal::Validate_right_side_category< GeomTraits, Cond_right::value > 
    Validate_right_side_category;
  void (Validate_right_side_category::*pright)(void) =
    &Validate_right_side_category::template missing__Right_side_category<int>;
  (void)pright;

  // Initialize the DCEL structure to represent an empty arrangement.
  m_topol_traits.init_dcel ();

  // Set the traits.
  m_geom_traits = static_cast<const Traits_adaptor_2*> (geom_traits);
  m_own_traits = false;
}

//-----------------------------------------------------------------------------
// Assignment operator.
//
template<class GeomTraits, class TopTraits>
Arrangement_on_surface_2<GeomTraits, TopTraits>&
Arrangement_on_surface_2<GeomTraits, TopTraits>::operator= (const Self& arr)
{
  // Check for self-assignment.
  if (this == &arr)
    return (*this);

  assign (arr);
  return (*this);
}

//-----------------------------------------------------------------------------
// Assign an arrangement.
//
template<class GeomTraits, class TopTraits>
void Arrangement_on_surface_2<GeomTraits, TopTraits>::assign (const Self& arr)
{
  // Clear the current contents of the arrangement.
  clear();

  // Notify the observers that an assignment is to take place.
  _notify_before_assign (arr);

  // Assign the topology-traits object.
  m_topol_traits.assign (arr.m_topol_traits);

  // Go over the vertices and create duplicates of the stored points.
  typename Dcel::Vertex_iterator  vit;
  Point_2                        *dup_p;
  DVertex                        *p_v;

  for (vit = _dcel().vertices_begin(); vit != _dcel().vertices_end(); ++vit)
  {
    p_v = &(*vit);

    if (! p_v->has_null_point())
    {
      // Create the duplicate point and store it in the points container.
      dup_p = _new_point (p_v->point());

      // Associate the vertex with the duplicated point.
      p_v->set_point (dup_p);
    }
  }

  // Go over the edge and create duplicates of the stored curves.
  typename Dcel::Edge_iterator     eit;
  X_monotone_curve_2              *dup_cv;
  DHalfedge                       *p_e;

  for (eit = _dcel().edges_begin(); eit != _dcel().edges_end(); ++eit)
  {
    p_e = &(*eit);

    if (! p_e->has_null_curve())
    {
      // Create the duplicate curve and store it in the curves container.
      dup_cv = _new_curve (p_e->curve());

      // Associate the halfedge (and its twin) with the duplicated curve. 
      p_e->set_curve (dup_cv);
    }
  }

  // Take care of the traits object.
  if (m_own_traits && m_geom_traits != NULL)
    delete m_geom_traits;

  if (arr.m_own_traits)
    m_geom_traits = new Traits_adaptor_2;
  else
    m_geom_traits = arr.m_geom_traits;
  m_own_traits = arr.m_own_traits;

  // Notify the observers that the assignment has been performed.
  _notify_after_assign ();
}

//-----------------------------------------------------------------------------
// Destructor.
//
template<class GeomTraits, class TopTraits>
Arrangement_on_surface_2<GeomTraits, TopTraits>::~Arrangement_on_surface_2 ()
{  
  // Free all stored points.
  typename Dcel::Vertex_iterator       vit;

  for (vit = _dcel().vertices_begin(); vit != _dcel().vertices_end(); ++vit)
  {
    if (! vit->has_null_point())
      _delete_point (vit->point());
  }

  // Free all stores curves.
  typename Dcel::Edge_iterator         eit;

  for (eit = _dcel().edges_begin(); eit != _dcel().edges_end(); ++eit)
  {
    if (! eit->has_null_curve())
      _delete_curve (eit->curve());
  }

  // Clear the DCEL.
  _dcel().delete_all();

  // Free the traits object, if necessary.
  if (m_own_traits)
    delete m_geom_traits;

  // Detach all observers still attached to the arrangement.
  Observers_iterator  iter = m_observers.begin();
  Observers_iterator  next;
  Observers_iterator  end = m_observers.end();

  while (iter != end)
  {
    next = iter;
    ++next;
    (*iter)->detach();
    iter = next;
  }
}

//-----------------------------------------------------------------------------
// Clear the arrangement.
//
template<class GeomTraits, class TopTraits>
void Arrangement_on_surface_2<GeomTraits, TopTraits>::clear ()
{
  // Notify the observers that we are about to clear the arragement.
  _notify_before_clear ();

  // Free all stored points.
  typename Dcel::Vertex_iterator       vit;

  for (vit = _dcel().vertices_begin(); vit != _dcel().vertices_end(); ++vit)
  {
    if (! vit->has_null_point())
      _delete_point (vit->point());
  }

  // Free all stores curves.
  typename Dcel::Edge_iterator         eit;

  for (eit = _dcel().edges_begin(); eit != _dcel().edges_end(); ++eit)
  {
    if (! eit->has_null_curve())
      _delete_curve (eit->curve());
  }

  // Clear the DCEL and construct an empty arrangement.
  _dcel().delete_all();
  m_topol_traits.init_dcel ();

  // Notify the observers that we have just cleared the arragement.
  _notify_after_clear ();

  return;
}

//-----------------------------------------------------------------------------
// Insert a point as an isolated vertex in the interior of a given face.
//
template<class GeomTraits, class TopTraits>
typename Arrangement_on_surface_2<GeomTraits, TopTraits>::Vertex_handle
Arrangement_on_surface_2<GeomTraits, TopTraits>::
insert_in_face_interior(const Point_2& p, Face_handle f)
{
  DFace          *p_f = _face (f);

#if CGAL_ARRANGEMENT_ON_SURFACE_INSERT_VERBOSE
    std::cout << "Aos_2: insert_in_face_interior (interface)" << std::endl;
    std::cout << "pt   : " << p << std::endl;
    std::cout << "face : " << &(*f) << std::endl;
#endif

  // Create a new vertex associated with the given point.
  // We assume the point has no boundary conditions.
  DVertex         *v = _create_vertex (p);
  Vertex_handle   vh (v);

  // Insert v as an isolated vertex inside the given face.
  _insert_isolated_vertex (p_f, v);

  // Return a handle to the new isolated vertex.
  return (vh);
}

//-----------------------------------------------------------------------------
// Insert an x-monotone curve into the arrangement as a new hole (inner
// component) inside the given face.
//
template<class GeomTraits, class TopTraits>
typename Arrangement_on_surface_2<GeomTraits, TopTraits>::Halfedge_handle
Arrangement_on_surface_2<GeomTraits, TopTraits>::
insert_in_face_interior(const X_monotone_curve_2& cv, Face_handle f)
{

#if CGAL_ARRANGEMENT_ON_SURFACE_INSERT_VERBOSE
    std::cout << "Aos_2: insert_in_face_interior (interface)" << std::endl;
    std::cout << "cv   : " << cv << std::endl;
    std::cout << "face : " << &(*f) << std::endl;
#endif

  DFace               *p_f = _face (f);

  // Check if cv's left end has boundary conditions, and obtain a vertex v1
  // that corresponds to this end.
  const Arr_parameter_space  ps_x1 =
    m_geom_traits->parameter_space_in_x_2_object()(cv, ARR_MIN_END);
  const Arr_parameter_space  ps_y1 =
    m_geom_traits->parameter_space_in_y_2_object()(cv, ARR_MIN_END);
  DHalfedge           *fict_prev1 = NULL;

  DVertex * v1 = (ps_x1 == ARR_INTERIOR && ps_y1 == ARR_INTERIOR) ?
    // The curve has a valid left endpoint: Create a new vertex associated
    // with the curve's left endpoint.
    _create_vertex (m_geom_traits->construct_min_vertex_2_object()(cv)) :
    // Locate the DCEL features that will be used for inserting the curve's
    // left end.
    _place_and_set_curve_end (p_f, cv, ARR_MIN_END, ps_x1, ps_y1, &fict_prev1);

  // Check if cv's right end has boundary conditions, and obtain a vertex v2
  // that corresponds to this end.
  const Arr_parameter_space  ps_x2 =
    m_geom_traits->parameter_space_in_x_2_object()(cv, ARR_MAX_END);
  const Arr_parameter_space  ps_y2 =
    m_geom_traits->parameter_space_in_y_2_object()(cv, ARR_MAX_END);
  DHalfedge           *fict_prev2 = NULL;

  DVertex * v2 = (ps_x2 == ARR_INTERIOR && ps_y2 == ARR_INTERIOR) ?
    // The curve has a valid right endpoint: Create a new vertex associated
    // with the curve's right endpoint.
    _create_vertex (m_geom_traits->construct_max_vertex_2_object()(cv)) :
    // Locate the DCEL features that will be used for inserting the curve's
    // right end.
    _place_and_set_curve_end (p_f, cv, ARR_MAX_END, ps_x2, ps_y2, &fict_prev2);


  // Create the edge connecting the two vertices (note we know v1 is 
  // lexicographically smaller than v2).
  DHalfedge       *new_he;

  if (fict_prev1 == NULL && fict_prev2 == NULL)
  {
    // Both vertices represent valid points.
    new_he = _insert_in_face_interior (cv, p_f, v1, v2, SMALLER);
  }
  else if (fict_prev1 == NULL && fict_prev2 != NULL)
  {
    // v1 represents a valid point and v2 is inserted using its predecessor.
    new_he = _insert_from_vertex (cv, fict_prev2, v1, LARGER);
    new_he = new_he->opposite();
  }
  else if (fict_prev1 != NULL && fict_prev2 == NULL)
  {
    // v1 is inserted using its predecessor and v2 represents a valid point.
    new_he = _insert_from_vertex (cv, fict_prev1, v2, SMALLER);
  }
  else
  {

    // Both vertices are inserted using their predecessor halfedges.
    // Note that in this case we may create a new face.

    bool        new_face_created = false;
  
    // This is a bug fix and we could not think of a better one. 
    // In case the inserted curve has two vertical asymptotes at the top 
    // it happens that fict_prev1 is split by the max end and becomes the 
    // prev edge, which is fict_prev2. Since both pointers are equal they 
    // both point to the max end. Thus, we advance fict_prev1 by one 
    // such that it points to the min end again. 
    // Note that this only happens at the top. At the bottom everything 
    // goes fine since the insertion order is reverted with respect to the 
    // orientation of the edges. 
    if (fict_prev1 == fict_prev2)
      fict_prev1 = fict_prev1->next();

    new_he = _insert_at_vertices (cv,
                                  fict_prev1, 
                                  fict_prev2,
                                  SMALLER,
                                  new_face_created);
  
    if (new_face_created)
    {
      CGAL_assertion (! new_he->is_on_inner_ccb());
      _relocate_in_new_face (new_he);
    }
  }

  // Return a handle to the new halfedge directed from left to right.
  return (Halfedge_handle (new_he));
}

//-----------------------------------------------------------------------------
// Insert an x-monotone curve into the arrangement, such that its left 
// endpoint corresponds to a given arrangement vertex.
//
template<class GeomTraits, class TopTraits>
typename Arrangement_on_surface_2<GeomTraits, TopTraits>::Halfedge_handle
Arrangement_on_surface_2<GeomTraits, TopTraits>::
insert_from_left_vertex(const X_monotone_curve_2& cv, 
                        Vertex_handle v,
                        Face_handle f)
{
  CGAL_precondition_code (
    const bool at_obnd1 = !m_geom_traits->is_closed_2_object()(cv, ARR_MIN_END));
  CGAL_precondition_msg 
    ((! at_obnd1 &&
      m_geom_traits->equal_2_object()
      (v->point(), 
       m_geom_traits->construct_min_vertex_2_object()(cv))) ||
     (at_obnd1 && v->is_at_open_boundary()),
     "The input vertex should be the left curve end.");

  // Check if cv's right end has boundary conditions. If not, create a vertex
  // that corresponds to the right endpoint.
  const Arr_parameter_space  ps_x2 =
    m_geom_traits->parameter_space_in_x_2_object()(cv, ARR_MAX_END);
  const Arr_parameter_space  ps_y2 =
    m_geom_traits->parameter_space_in_y_2_object()(cv, ARR_MAX_END);
  DVertex             *v2 = NULL;
  DHalfedge           *fict_prev2 = NULL;

  if (ps_x2 == ARR_INTERIOR && ps_y2 == ARR_INTERIOR)
  {
    // The curve has a valid right endpoint: Create a new vertex associated
    // with the curve's right endpoint.
    v2 = _create_vertex (m_geom_traits->construct_max_vertex_2_object()(cv));
  }

  // Check if the given vertex, corresponding to the left curve end, has no
  // incident edges. 
  if (v->degree() == 0)
  {
    // The given vertex is an isolated one: We should in fact insert the curve
    // in the interior of the face containing this vertex.
    DVertex      *v1 = _vertex (v);
    DIso_vertex  *iv = NULL;
    DFace        *p_f = NULL;

    if (v->is_isolated())
    {
      // Obtain the face from the isolated vertex.
      iv = v1->isolated_vertex();
      p_f = iv->face();
    }
    else
    {
      // In this case the face that contains v should be provided by the user.
      CGAL_precondition (f != Face_handle());
      p_f = _face (f);
    }

    // If the vertex that corresponds to cv's right end has boundary
    // conditions, create it now.
    if (v2 == NULL)
    {    
      // Locate the DCEL features that will be used for inserting the curve's
      // right end.
      v2 = _place_and_set_curve_end (p_f, cv, ARR_MAX_END, ps_x2, ps_y2,
                                     &fict_prev2);
    }

    if (iv != NULL)
    {
      // Remove the isolated vertex v1, as it will not be isolated any more.
      p_f->erase_isolated_vertex (iv);
      _dcel().delete_isolated_vertex (iv);
    }

    // Create the edge connecting the two vertices (note that we know that
    // v1 is smaller than v2).
    DHalfedge  *new_he;

    if (fict_prev2 == NULL)
    {
      new_he = _insert_in_face_interior (cv, p_f, v1, v2,
                                         SMALLER);
    }
    else
    {
      new_he = _insert_from_vertex (cv,
                                    fict_prev2, v1,
                                    LARGER);
      new_he = new_he->opposite();
    }

    // Return a handle to the new halfedge directed from v1 to v2.
    return (Halfedge_handle (new_he));
  }

  // Go over the incident halfedges around v and find the halfedge after
  // which the new curve should be inserted.
  DHalfedge  *prev1 = _locate_around_vertex (_vertex (v), cv, ARR_MIN_END);
  DFace      *f1 = prev1->is_on_inner_ccb() ? prev1->inner_ccb()->face() :
                                              prev1->outer_ccb()->face();

  CGAL_assertion_msg
    (prev1 != NULL,
     "The inserted curve cannot be located in the arrangement.");

  // If the vertex that corresponds to cv's right end has boundary conditions,
  // create it now.
  if (v2 == NULL)
  {    
    // Locate the DCEL features that will be used for inserting the curve's
    // right end.
    v2 = _place_and_set_curve_end (f1, cv, ARR_MAX_END, ps_x2, ps_y2,
                                   &fict_prev2);
  }

  // Perform the insertion (note that we know that prev1->vertex is smaller
  // than v2).
  DHalfedge  *new_he;

  if (fict_prev2 == NULL)
  {
    // Insert the halfedge given the predecessor halfedge of v1.
    new_he = _insert_from_vertex (cv, prev1, v2, SMALLER);
  }
  else
  {
    // Insert the halfedge given the two predecessor halfedges.
    // Note that in this case we may create a new face.
    bool        new_face_created = false;

    new_he = _insert_at_vertices (cv,
                                  prev1, 
                                  fict_prev2,
                                  SMALLER,
                                  new_face_created);
  
    if (new_face_created)
    {
      CGAL_assertion (! new_he->is_on_inner_ccb());
      _relocate_in_new_face (new_he);
    }
  }

  // Return a handle to the halfedge directed toward the new vertex v2.
  return (Halfedge_handle (new_he));
}

//-----------------------------------------------------------------------------
// Insert an x-monotone curve into the arrangement, such that one its left
// endpoint corresponds to a given arrangement vertex, given the exact place
// for the curve in the circular list around this vertex.
//
template<class GeomTraits, class TopTraits>
typename Arrangement_on_surface_2<GeomTraits, TopTraits>::Halfedge_handle
Arrangement_on_surface_2<GeomTraits, TopTraits>::
insert_from_left_vertex(const X_monotone_curve_2& cv,
                        Halfedge_handle prev)
{

#if CGAL_ARRANGEMENT_ON_SURFACE_INSERT_VERBOSE
    std::cout << "Aos_2: insert_from_left_vertex (interface)" << std::endl;
    std::cout << "cv   : " << cv << std::endl;
    if (!prev->is_fictitious()) {
        std::cout << "prev : " << prev ->curve() << std::endl;
    } else {
      std::cout << "prev : fictitious" << std::endl;
  }
  std::cout << "dir  : " << prev ->direction() << std::endl;
#endif

  CGAL_precondition_code (
    const bool at_obnd1 = !m_geom_traits->is_closed_2_object()(cv, ARR_MIN_END));
  CGAL_precondition_msg
    ((! at_obnd1 &&
      m_geom_traits->equal_2_object()
      (prev->target()->point(), 
       m_geom_traits->construct_min_vertex_2_object()(cv))) ||
     (at_obnd1 && prev->target()->is_at_open_boundary()),
     "The target of the input halfedge should be the left curve end.");

  CGAL_precondition_msg
    (at_obnd1 || _locate_around_vertex (_vertex(prev->target()),
                                       cv, ARR_MIN_END) == _halfedge(prev),
     "In the clockwise order of curves around the vertex, "
     " cv must succeed the curve of prev.");

  // Get the predecessor halfedge for the insertion of the left curve end.
  DHalfedge           *prev1 = _halfedge (prev);
  DFace               *f1 = _face (prev->face());

  // Check if cv's right end has boundary conditions, and obtain a vertex
  // that corresponds to this end.
  const Arr_parameter_space  ps_x2 =
    m_geom_traits->parameter_space_in_x_2_object()(cv, ARR_MAX_END);
  const Arr_parameter_space  ps_y2 =
    m_geom_traits->parameter_space_in_y_2_object()(cv, ARR_MAX_END);
  DHalfedge           *fict_prev2 = NULL;

  DVertex * v2 = (ps_x2 == ARR_INTERIOR && ps_y2 == ARR_INTERIOR) ?
    // The curve has a valid right endpoint: Create a new vertex associated
    // with the curve's right endpoint.
    _create_vertex (m_geom_traits->construct_max_vertex_2_object()(cv)) :
    // Locate the DCEL features that will be used for inserting the curve's
    // right end.
    _place_and_set_curve_end (f1, cv, ARR_MAX_END, ps_x2, ps_y2, &fict_prev2);

  // Perform the insertion (note that we know that prev1->vertex is smaller
  // than v2).
  DHalfedge  *new_he;

  if (fict_prev2 == NULL)
  {
    // Insert the halfedge given the predecessor halfedge of the left vertex.
    new_he = _insert_from_vertex (cv, prev1, v2, SMALLER);
  }
  else
  {
    // Insert the halfedge given the two predecessor halfedges.
    // Note that in this case we may create a new face.
    bool        new_face_created = false;

    new_he = _insert_at_vertices (cv, prev1, fict_prev2, SMALLER,
                                  new_face_created);
  
    if (new_face_created)
    {
      CGAL_assertion (! new_he->is_on_inner_ccb());
      _relocate_in_new_face (new_he);
    }
  }

  // Return a handle to the halfedge directed toward the new vertex v2.
  return (Halfedge_handle (new_he));
}

//-----------------------------------------------------------------------------
// Insert an x-monotone curve into the arrangement, such that its right
// endpoint corresponds to a given arrangement vertex.
//
template<class GeomTraits, class TopTraits>
typename Arrangement_on_surface_2<GeomTraits, TopTraits>::Halfedge_handle
Arrangement_on_surface_2<GeomTraits, TopTraits>::
insert_from_right_vertex(const X_monotone_curve_2& cv,
                         Vertex_handle v, Face_handle f)
{
  CGAL_precondition_code
    (const bool at_obnd2 = 
     !m_geom_traits->is_closed_2_object()(cv, ARR_MAX_END));
  CGAL_precondition_msg 
    ((! at_obnd2 &&
      m_geom_traits->equal_2_object()
      (v->point(), 
       m_geom_traits->construct_max_vertex_2_object()(cv))) ||
     (at_obnd2 && v->is_at_open_boundary()),
     "The input vertex should be the right curve end.");

  // Check if cv's left end has boundary conditions. If not, create a vertex
  // that corresponds to the left endpoint.
  const Arr_parameter_space  ps_x1 =
    m_geom_traits->parameter_space_in_x_2_object()(cv, ARR_MIN_END);
  const Arr_parameter_space  ps_y1 =
    m_geom_traits->parameter_space_in_y_2_object()(cv, ARR_MIN_END);
  DVertex             *v1 = NULL;
  DHalfedge           *fict_prev1 = NULL;

  if (ps_x1 == ARR_INTERIOR && ps_y1 == ARR_INTERIOR)
  {
    // The curve has a valid left endpoint: Create a new vertex associated
    // with the curve's left endpoint.
    v1 = _create_vertex (m_geom_traits->construct_min_vertex_2_object()(cv));
  }

  // Check if the given vertex, corresponding to the right curve end, has no
  // incident edges. 
  if (v->degree() == 0)
  {
    // The given vertex is an isolated one: We should in fact insert the curve
    // in the interior of the face containing this vertex.
    DVertex      *v2 = _vertex (v);
    DIso_vertex  *iv = NULL;
    DFace        *p_f = NULL;

    if (v->is_isolated())
    {
      // Obtain the face from the isolated vertex.
      iv = v2->isolated_vertex();
      p_f = iv->face();
    }
    else
    {
      // In this case the face that contains v should be provided by the user.
      CGAL_precondition (f != Face_handle());
      p_f = _face (f);
    }

    // If the vertex that corresponds to cv's left end has boundary
    // conditions, create it now.
    if (v1 == NULL)
    {    
      // Locate the DCEL features that will be used for inserting the curve's
      // left end.
      v1 = _place_and_set_curve_end (p_f, cv, ARR_MIN_END, ps_x1, ps_y1,
                                     &fict_prev1);
    }

    if (iv != NULL)
    {
      // Remove the isolated vertex v2, as it will not be isolated any more.
      p_f->erase_isolated_vertex (iv);
      _dcel().delete_isolated_vertex (iv);
    }

    // Create the edge connecting the two vertices (note that we know that
    // v1 is smaller than v2).
    DHalfedge  *new_he;

    if (fict_prev1 == NULL)
    {
      new_he = _insert_in_face_interior (cv, p_f, v1, v2, SMALLER);
    }
    else
    {
      new_he = _insert_from_vertex (cv, fict_prev1, v2, SMALLER);
    }
    
    // Return a handle to the new halfedge whose target is the new vertex v1.
    return (Halfedge_handle (new_he->opposite()));
  }

  // Go over the incident halfedges around v and find the halfedge after
  // which the new curve should be inserted.
  DHalfedge  *prev2 = _locate_around_vertex (_vertex (v), cv, ARR_MAX_END);
  DFace      *f2 = prev2->is_on_inner_ccb() ? prev2->inner_ccb()->face() :
                                              prev2->outer_ccb()->face();

  CGAL_assertion_msg
    (prev2 != NULL,
     "The inserted curve cannot be located in the arrangement.");

  // If the vertex that corresponds to cv's left end has boundary conditions,
  // create it now.
  if (v1 == NULL)
  {    
    // Locate the DCEL features that will be used for inserting the curve's
    // left end.
    v1 = _place_and_set_curve_end (f2, cv, ARR_MIN_END, ps_x1, ps_y1,
                                   &fict_prev1);
  }

  // Perform the insertion (note that we know that prev2->vertex is larger
  // than v1).
  DHalfedge  *new_he;

  if (fict_prev1 == NULL)
  {
    // Insert the halfedge given the predecessor halfedge of v2.
    new_he = _insert_from_vertex (cv, prev2, v1, LARGER);
  }
  else
  {
    // Insert the halfedge given the two predecessor halfedges.
    // Note that in this case we may create a new face.
    bool        new_face_created = false;

    new_he = _insert_at_vertices (cv, prev2, fict_prev1, LARGER,
                                  new_face_created);
  
    if (new_face_created)
    {
      CGAL_assertion (! new_he->is_on_inner_ccb());
      _relocate_in_new_face (new_he);
    }
  }

  // Return a handle to the halfedge directed toward the new vertex v1.
  return (Halfedge_handle (new_he));
}

//-----------------------------------------------------------------------------
// Insert an x-monotone curve into the arrangement, such that its right
// endpoint corresponds to a given arrangement vertex, given the exact place
// for the curve in the circular list around this vertex.
//
template<class GeomTraits, class TopTraits>
typename Arrangement_on_surface_2<GeomTraits, TopTraits>::Halfedge_handle
Arrangement_on_surface_2<GeomTraits, TopTraits>::
insert_from_right_vertex(const X_monotone_curve_2& cv,
                         Halfedge_handle prev)
{

#if CGAL_ARRANGEMENT_ON_SURFACE_INSERT_VERBOSE
    std::cout << "Aos_2: insert_from_right_vertexs (interface)" << std::endl;
    std::cout << "cv   : " << cv << std::endl;
    if (!prev->is_fictitious()) {
      std::cout << "prev : " << prev ->curve() << std::endl;
  } else {
      std::cout << "prev : fictitious" << std::endl;
  }
  std::cout << "dir  : " << prev->direction() << std::endl;
#endif

  CGAL_precondition_code
    (const bool at_obnd2 = 
     !m_geom_traits->is_closed_2_object()(cv, ARR_MAX_END));
  CGAL_precondition_msg 
    ((! at_obnd2 &&
      m_geom_traits->equal_2_object()
      (prev->target()->point(), 
       m_geom_traits->construct_max_vertex_2_object()(cv))) ||
     (at_obnd2 && prev->target()->is_at_open_boundary()),
     "The input vertex should be the right curve end.");

  CGAL_precondition_msg
    (at_obnd2 ||
     _locate_around_vertex (_vertex(prev->target()), cv, ARR_MAX_END) ==
     _halfedge(prev),
     "In the clockwise order of curves around the vertex, "
     " cv must succeed the curve of prev.");

  // Get the predecessor halfedge for the insertion of the right curve end.
  DHalfedge           *prev2 = _halfedge (prev);
  DFace               *f2 = _face (prev->face());

  // Check if cv's left end has boundary conditions, and obtain a vertex v1
  // that corresponds to this end.
  const Arr_parameter_space  ps_x1 =
    m_geom_traits->parameter_space_in_x_2_object()(cv, ARR_MIN_END);
  const Arr_parameter_space  ps_y1 =
    m_geom_traits->parameter_space_in_y_2_object()(cv, ARR_MIN_END);
  DHalfedge           *fict_prev1 = NULL;

  DVertex * v1 = (ps_x1 == ARR_INTERIOR && ps_y1 == ARR_INTERIOR) ?
    // The curve has a valid left endpoint: Create a new vertex associated
    // with the curve's left endpoint.
    _create_vertex (m_geom_traits->construct_min_vertex_2_object()(cv)) :
    // Locate the DCEL features that will be used for inserting the curve's
    // left end.
    _place_and_set_curve_end (f2, cv, ARR_MIN_END, ps_x1, ps_y1, &fict_prev1);

  // Perform the insertion (note that we know that prev2->vertex is larger
  // than v1).
  DHalfedge  *new_he;

  if (fict_prev1 == NULL)
  {
    // Insert the halfedge given the predecessor halfedge of the right vertex.
    new_he = _insert_from_vertex (cv, prev2, v1, LARGER);
  }
  else
  {
    // Insert the halfedge given the two predecessor halfedges.
    // Note that in this case we may create a new face.
    bool        new_face_created = false;

    new_he = _insert_at_vertices (cv, prev2, fict_prev1, LARGER,
                                  new_face_created);
  
    if (new_face_created)
    {
      CGAL_assertion (! new_he->is_on_inner_ccb());
      _relocate_in_new_face (new_he);
    }
  }

  // Return a handle to the halfedge directed toward the new vertex v1.
  return (Halfedge_handle (new_he));
}

//-----------------------------------------------------------------------------
// Insert an x-monotone curve into the arrangement, such that both its
// endpoints corresponds to a given arrangement vertex.
//
template<class GeomTraits, class TopTraits>
typename Arrangement_on_surface_2<GeomTraits, TopTraits>::Halfedge_handle
Arrangement_on_surface_2<GeomTraits, TopTraits>::
insert_at_vertices(const X_monotone_curve_2& cv,
                   Vertex_handle v1, Vertex_handle v2,
                   Face_handle f)
{
  // Determine which one of the given vertices matches the left end of the
  // given curve.
  const bool at_obnd1 = !m_geom_traits->is_closed_2_object()(cv, ARR_MIN_END);
  const bool at_obnd2 = !m_geom_traits->is_closed_2_object()(cv, ARR_MAX_END);

  Arr_curve_end      ind1;
  Arr_curve_end      ind2;

  if (! at_obnd1)
  {
    CGAL_precondition_code (Vertex_handle v_right);

    if (! v1->is_at_open_boundary() &&
        m_geom_traits->equal_2_object()
        (v1->point(), 
         m_geom_traits->construct_min_vertex_2_object()(cv)))
    {
      ind1 = ARR_MIN_END;
      ind2 = ARR_MAX_END;
      CGAL_precondition_code ( v_right = v2; );
    }
    else
    {
      CGAL_precondition_msg
        (! v2->is_at_open_boundary() &&
         m_geom_traits->equal_2_object()
         (v2->point(), 
          m_geom_traits->construct_min_vertex_2_object()(cv)),
         "One of the input vertices should be the left curve end.");

      ind1 = ARR_MAX_END;
      ind2 = ARR_MIN_END;
      CGAL_precondition_code ( v_right = v1; );
    }

    CGAL_precondition_msg 
      ((! at_obnd2 &&
        m_geom_traits->equal_2_object()
        (v_right->point(), 
         m_geom_traits->construct_max_vertex_2_object()(cv))) ||
       (at_obnd2 && v_right->is_at_open_boundary()),
       "One of the input vertices should be the right curve end.");
  }
  else
  {
    if (! at_obnd2)
    {
      CGAL_precondition_code (Vertex_handle v_left);

      if (! v1->is_at_open_boundary() &&
          m_geom_traits->equal_2_object()
          (v1->point(), 
           m_geom_traits->construct_max_vertex_2_object()(cv)))
      {
        ind1 = ARR_MAX_END;
        ind2 = ARR_MIN_END;
        CGAL_precondition_code ( v_left = v2; );
      }
      else
      {
        CGAL_precondition_msg
          (! v2->is_at_open_boundary() &&
           m_geom_traits->equal_2_object()
           (v2->point(), 
            m_geom_traits->construct_max_vertex_2_object()(cv)),
           "One of the input vertices should be the right curve end.");

        ind1 = ARR_MIN_END;
        ind2 = ARR_MAX_END;
        CGAL_precondition_code ( v_left = v1; );
      }

      CGAL_precondition_msg
        (at_obnd1 && v_left->is_at_open_boundary(),
         "One of the input vertices should be the left curve end.");
    }
    else
    {
      Arr_parameter_space  ps_x1 =
        m_geom_traits->parameter_space_in_x_2_object()(cv, ARR_MIN_END);
      Arr_parameter_space  ps_y1 =
        m_geom_traits->parameter_space_in_y_2_object()(cv, ARR_MIN_END);

      // Check which vertex should be associated with the minimal curve-end
      // (so the other is associated with the maximal curve-end).
      if (m_topol_traits.are_equal (_vertex (v1), cv, ARR_MIN_END, ps_x1, ps_y1))
      {
        CGAL_assertion
          (m_topol_traits.are_equal
           (_vertex (v2), cv, ARR_MAX_END,
            m_geom_traits->parameter_space_in_x_2_object()(cv, ARR_MAX_END),
            m_geom_traits->parameter_space_in_y_2_object()(cv, ARR_MAX_END)));

        ind1 = ARR_MIN_END;
        ind2 = ARR_MAX_END;
      }
      else
      {
        CGAL_assertion (m_topol_traits.are_equal
                        (_vertex (v2), cv, ARR_MIN_END, ps_x1, ps_y1));
        CGAL_assertion
          (m_topol_traits.are_equal
           (_vertex (v1), cv, ARR_MAX_END,
            m_geom_traits->parameter_space_in_x_2_object()(cv, ARR_MAX_END),
            m_geom_traits->parameter_space_in_y_2_object()(cv, ARR_MAX_END)));

        ind1 = ARR_MAX_END;
        ind2 = ARR_MIN_END;
      }
    }
  }

  // Check whether one of the vertices has no incident halfedges.
  if (v1->degree() == 0)
  {
    // Get the face containing the isolated vertex v1.
    DVertex      *p_v1 = _vertex (v1);
    DIso_vertex  *iv1 = NULL;
    DFace        *f1 = NULL;

    if (p_v1->is_isolated())
    {
      // Obtain the containing face from the isolated vertex record.
      iv1 = p_v1->isolated_vertex();
      f1 = iv1->face();

      // Remove the isolated vertex v1, as it will not be isolated any more.
      f1->erase_isolated_vertex (iv1);
      _dcel().delete_isolated_vertex (iv1);
    }

    // Check whether the other vertex also has no incident halfedges.
    if (v2->degree() == 0)
    {
      // Both end-vertices are isolated. Make sure they are contained inside
      // the same face.
      DVertex      *p_v2 = _vertex (v2);
      DIso_vertex  *iv2 = NULL;
      DFace        *f2 = NULL;

      if (p_v2->is_isolated())
      {
        // Obtain the containing face from the isolated vertex record.
        iv2 = p_v2->isolated_vertex();
        f2 = iv2->face();

        CGAL_assertion_msg
          (f1 == NULL || f1 == f2,
           "The two isolated vertices must be located inside the same face.");

        // Remove the isolated vertex v2, as it will not be isolated any more.
        f2->erase_isolated_vertex (iv2);
        _dcel().delete_isolated_vertex (iv2);
      }
      else if(f1 == NULL)
      {
        // In this case the containing face must be given by the user.
        CGAL_precondition (f != Face_handle());
      }

      // Perform the insertion.
      Comparison_result  res = (ind1 == ARR_MIN_END) ? SMALLER : LARGER;
      DHalfedge * new_he = _insert_in_face_interior (cv, f1, p_v1, p_v2, res);

      return (Halfedge_handle (new_he));
    }

    // Go over the incident halfedges around v2 and find the halfedge after
    // which the new curve should be inserted.
    DHalfedge  *prev2 = _locate_around_vertex (_vertex (v2), cv, ind2);
    CGAL_assertion_code (
      DFace      *f2 = prev2->is_on_inner_ccb() ? prev2->inner_ccb()->face() :
                                                  prev2->outer_ccb()->face();
    );

    CGAL_assertion_msg
      (prev2 != NULL,
       "The inserted curve cannot be located in the arrangement.");

    CGAL_assertion_msg
      (f1 == NULL || f1 == f2,
       "The inserted curve should not intersect the existing arrangement.");

    // Perform the insertion. Note that the returned halfedge is directed
    // toward v1 (and not toward v2), so we return its twin.
    Comparison_result  res = (ind2 == ARR_MIN_END) ? SMALLER : LARGER;
    DHalfedge * new_he = _insert_from_vertex (cv, prev2, p_v1, res);

    return (Halfedge_handle (new_he->opposite()));
  }
  else if (v2->degree() == 0)
  {
    // Get the face containing the isolated vertex v2.
    DVertex      *p_v2 = _vertex (v2);
    DIso_vertex  *iv2 = NULL;
    DFace        *f2 = NULL;

    if (v2->is_isolated())
    {
      // Obtain the containing face from the isolated vertex record.
      iv2 = p_v2->isolated_vertex();
      f2 = iv2->face();

      // Remove the isolated vertex v2, as it will not be isolated any more.
      f2->erase_isolated_vertex (iv2);
      _dcel().delete_isolated_vertex (iv2);
    }

    // Go over the incident halfedges around v1 and find the halfedge after
    // which the new curve should be inserted.
    DHalfedge  *prev1 = _locate_around_vertex (_vertex (v1), cv, ind1);
    CGAL_assertion_code (
      DFace      *f1 = prev1->is_on_inner_ccb() ? prev1->inner_ccb()->face() :
                                                  prev1->outer_ccb()->face();
    );

    CGAL_assertion_msg
      (prev1 != NULL,
       "The inserted curve cannot be located in the arrangement.");

    CGAL_assertion_msg
      (f2 == NULL || f2 == f1,
       "The inserted curve should not intersect the existing arrangement.");

    // Perform the insertion.
    Comparison_result  res = (ind1 == ARR_MIN_END) ? SMALLER : LARGER;
    DHalfedge         *new_he = _insert_from_vertex (cv, prev1, p_v2, res);

    return (Halfedge_handle (new_he));
  }

  // Go over the incident halfedges around v1 and v2 and find the two
  // halfedges after which the new curve should be inserted, respectively.
  DHalfedge  *prev1 = _locate_around_vertex (_vertex (v1), cv, ind1);
  DHalfedge  *prev2 = _locate_around_vertex (_vertex (v2), cv, ind2);

  CGAL_assertion_msg
    (prev1 != NULL && prev2 != NULL,
     "The inserted curve cannot be located in the arrangement.");

  // Perform the insertion.
  return insert_at_vertices(cv, Halfedge_handle(prev1), Halfedge_handle(prev2));
}

//-----------------------------------------------------------------------------
// Insert an x-monotone curve into the arrangement, such that both its
// endpoints correspond to given arrangement vertices, given the exact
// place for the curve in one of the circular lists around a vertex.
//
template<class GeomTraits, class TopTraits>
typename Arrangement_on_surface_2<GeomTraits, TopTraits>::Halfedge_handle
Arrangement_on_surface_2<GeomTraits, TopTraits>::
insert_at_vertices(const X_monotone_curve_2& cv,
                   Halfedge_handle prev1,
                   Vertex_handle v2)
{
  // Determine which one of the given vertices mathces the left end of the
  // given curve.
  const bool at_obnd1 = !m_geom_traits->is_closed_2_object()(cv, ARR_MIN_END);
  const bool at_obnd2 = !m_geom_traits->is_closed_2_object()(cv, ARR_MAX_END);

  Arr_curve_end      ind2;

  if (! at_obnd1)
  {
    CGAL_precondition_code ( Vertex_handle  v_right; );

    if (! prev1->target()->is_at_open_boundary() &&
        m_geom_traits->equal_2_object()
        (prev1->target()->point(), 
         m_geom_traits->construct_min_vertex_2_object()(cv)))
    {
      ind2 = ARR_MAX_END;
      CGAL_precondition_code ( v_right = v2; );
    }
    else
    {
      CGAL_precondition_msg
        (! v2->is_at_open_boundary() &&
         m_geom_traits->equal_2_object()
         (v2->point(), 
          m_geom_traits->construct_min_vertex_2_object()(cv)),
         "One of the input vertices should be the left curve end.");

      ind2 = ARR_MIN_END;
      CGAL_precondition_code ( v_right = prev1->target(); );
    }

    CGAL_precondition_msg 
      ((! at_obnd2 &&
        m_geom_traits->equal_2_object()
        (v_right->point(), 
         m_geom_traits->construct_max_vertex_2_object()(cv))) ||
       (at_obnd2 && v_right->is_at_open_boundary()),
       "One of the input vertices should be the right curve end.");
  }
  else
  {
    if (! at_obnd2)
    {
      CGAL_precondition_code ( Vertex_handle  v_left; );

      if (! prev1->target()->is_at_open_boundary() &&
          m_geom_traits->equal_2_object()
          (prev1->target()->point(), 
           m_geom_traits->construct_max_vertex_2_object()(cv)))
      {
        ind2 = ARR_MIN_END;
        CGAL_precondition_code ( v_left = v2; );
      }
      else
      {
        CGAL_precondition_msg
          (! v2->is_at_open_boundary() &&
           m_geom_traits->equal_2_object()
           (v2->point(), 
            m_geom_traits->construct_max_vertex_2_object()(cv)),
           "One of the input vertices should be the right curve end.");

        ind2 = ARR_MAX_END;
        CGAL_precondition_code ( v_left = prev1->target(); );
      }

      CGAL_precondition_msg
        (at_obnd1 && v_left->is_at_open_boundary(),
         "One of the input vertices should be the left curve end.");
    }
    else
    {
      Arr_parameter_space  ps_x1 =
        m_geom_traits->parameter_space_in_x_2_object()(cv, ARR_MIN_END);
      Arr_parameter_space  ps_y1 =
        m_geom_traits->parameter_space_in_y_2_object()(cv, ARR_MIN_END);

      // Check which vertex should be associated with the minimal curve-end
      // (so the other is associated with the maximal curve-end).
      if (m_topol_traits.are_equal (_vertex (prev1->target()),
                                cv, ARR_MIN_END, ps_x1, ps_y1))
      {
        CGAL_assertion
          (m_topol_traits.are_equal
           (_vertex (v2), cv, ARR_MAX_END,
            m_geom_traits->parameter_space_in_x_2_object()(cv, ARR_MAX_END),
            m_geom_traits->parameter_space_in_y_2_object()(cv, ARR_MAX_END)));

        ind2 = ARR_MAX_END;
      }
      else
      {
        CGAL_assertion (m_topol_traits.are_equal
                        (_vertex (v2), cv, ARR_MIN_END, ps_x1, ps_y1));
        CGAL_assertion
          (m_topol_traits.are_equal
           (_vertex (prev1->target()), cv, ARR_MAX_END,
            m_geom_traits->parameter_space_in_x_2_object()(cv, ARR_MAX_END),
            m_geom_traits->parameter_space_in_y_2_object()(cv, ARR_MAX_END)));

        ind2 = ARR_MIN_END;
      }
    }
  }

  // Check whether v2 is has no incident edges.
  if (v2->degree() == 0)
  {
    // Get the face containing the isolated vertex v2.
    DVertex      *p_v2 = _vertex (v2);
    DIso_vertex  *iv2 = NULL;
    DFace        *f2 = NULL;

    if (v2->is_isolated())
    {
      iv2 = p_v2->isolated_vertex();
      f2 = iv2->face();

      CGAL_assertion_msg
        (f2 == _face (prev1->face()),
         "The inserted curve should not intersect the existing arrangement.");

      // Remove the isolated vertex v2, as it will not be isolated any more.
      f2->erase_isolated_vertex (iv2);
      _dcel().delete_isolated_vertex (iv2);
    }

    // Perform the insertion.
    Comparison_result  res = (ind2 == ARR_MAX_END) ? SMALLER : LARGER;
    DHalfedge         *new_he = _insert_from_vertex (cv,
                                                     _halfedge (prev1), p_v2,
                                                     res);

    return (Halfedge_handle (new_he));
  }

  // Go over the incident halfedges around v2 and find the halfedge after
  // which the new curve should be inserted.
  DHalfedge  *prev2 = _locate_around_vertex (_vertex (v2), cv, ind2);

  CGAL_assertion_msg
    (prev2 != NULL,
     "The inserted curve cannot be located in the arrangement.");

  // Perform the insertion.
  return (insert_at_vertices (cv,
                              prev1,
                              Halfedge_handle(prev2)));
}

//-----------------------------------------------------------------------------
// Insert an x-monotone curve into the arrangement, such that both its
// endpoints correspond to given arrangement vertices, given the exact
// place for the curve in both circular lists around these two vertices.
//
template<class GeomTraits, class TopTraits>
typename Arrangement_on_surface_2<GeomTraits, TopTraits>::Halfedge_handle
Arrangement_on_surface_2<GeomTraits, TopTraits>::
insert_at_vertices(const X_monotone_curve_2& cv,
                   Halfedge_handle prev1,
                   Halfedge_handle prev2)
{
#if CGAL_ARRANGEMENT_ON_SURFACE_INSERT_VERBOSE
    std::cout << "Aos_2: insert_at_vertices (interface)" << std::endl;
    std::cout << "cv   : " << cv << std::endl;
    if (!prev1->is_fictitious()) {
        std::cout << "prev1: " << prev1->curve() << std::endl;
    } else {
      std::cout << "prev1: fictitious" << std::endl;
  }
  std::cout << "dir1 : " << prev1->direction() << std::endl;
  if (!prev2->is_fictitious()) {
      std::cout << "prev2: " << prev2->curve() << std::endl;
  } else {
      std::cout << "prev2: fictitious" << std::endl;
  }
  std::cout << "dir2 : " << prev2->direction() << std::endl;
#endif
    
  // Determine which one of the given vertices (the target vertices of the
  // given halfedges) mathces the left end of the given curve.
  // Thus, we can determine the comparison result between prev1->target()
  // and prev2->target().
  const bool at_obnd1 = !m_geom_traits->is_closed_2_object()(cv, ARR_MIN_END);
  const bool at_obnd2 = !m_geom_traits->is_closed_2_object()(cv, ARR_MAX_END);
  Comparison_result  res;

  if (! at_obnd1)
  {
    CGAL_precondition_code ( Vertex_handle  v_right; );

    if (! prev1->target()->is_at_open_boundary() &&
        m_geom_traits->equal_2_object()
        (prev1->target()->point(), 
         m_geom_traits->construct_min_vertex_2_object()(cv)))
    {
      res = SMALLER;
      CGAL_precondition_code ( v_right = prev2->target(); );
    }
    else
    {
      CGAL_precondition_msg
        (! prev2->target()->is_at_open_boundary() &&
         m_geom_traits->equal_2_object()
         (prev2->target()->point(), 
          m_geom_traits->construct_min_vertex_2_object()(cv)),
         "One of the input vertices should be the left curve end.");

      res = LARGER;
      CGAL_precondition_code ( v_right = prev1->target(); );
    }

    CGAL_precondition_msg 
      ((! at_obnd2 &&
        m_geom_traits->equal_2_object()
        (v_right->point(), 
         m_geom_traits->construct_max_vertex_2_object()(cv))) ||
       (at_obnd2 && v_right->is_at_open_boundary()),
       "One of the input vertices should be the right curve end.");
  }
  else
  {
    if (! at_obnd2)
    {
      CGAL_precondition_code ( Vertex_handle  v_left; );

      if (! prev1->target()->is_at_open_boundary() &&
          m_geom_traits->equal_2_object()
          (prev1->target()->point(), 
           m_geom_traits->construct_max_vertex_2_object()(cv)))
      {
        res = LARGER;
        CGAL_precondition_code ( v_left = prev2->target(); );
      }
      else
      {
        CGAL_precondition_msg
          (! prev2->target()->is_at_open_boundary() &&
           m_geom_traits->equal_2_object()
           (prev2->target()->point(), 
            m_geom_traits->construct_max_vertex_2_object()(cv)),
           "One of the input vertices should be the right curve end.");

        res = SMALLER;
        CGAL_precondition_code ( v_left = prev1->target(); );
      }

      CGAL_precondition_msg
        (at_obnd1 && v_left->is_at_open_boundary(),
         "One of the input vertices should be the left curve end.");
    }
    else
    {
      Arr_parameter_space  ps_x1 =
        m_geom_traits->parameter_space_in_x_2_object()(cv, ARR_MIN_END);
      Arr_parameter_space  ps_y1 =
        m_geom_traits->parameter_space_in_y_2_object()(cv, ARR_MIN_END);
      // Check which vertex should be associated with the minimal curve-end
      // (so the other is associated with the maximal curve-end), and
      // determine the comparison result of the two vertices accordingly.
      if (m_topol_traits.are_equal (_vertex (prev1->target()),
                                cv, ARR_MIN_END, ps_x1, ps_y1))
      {
        CGAL_assertion
          (m_topol_traits.are_equal
           (_vertex (prev2->target()), cv, ARR_MAX_END,
            m_geom_traits->parameter_space_in_x_2_object()(cv, ARR_MAX_END),
            m_geom_traits->parameter_space_in_y_2_object()(cv, ARR_MAX_END)));

        res = SMALLER;
      }
      else
      {
        CGAL_assertion
          (m_topol_traits.are_equal
           (_vertex (prev2->target()), cv, ARR_MIN_END,
            m_geom_traits->parameter_space_in_x_2_object()(cv, ARR_MIN_END),
            m_geom_traits->parameter_space_in_y_2_object()(cv, ARR_MIN_END)));
        CGAL_assertion
          (m_topol_traits.are_equal
           (_vertex (prev1->target()), cv, ARR_MAX_END,
            m_geom_traits->parameter_space_in_x_2_object()(cv, ARR_MAX_END),
            m_geom_traits->parameter_space_in_y_2_object()(cv, ARR_MAX_END)));

        res = LARGER;
      }
    }
  }

  // Check if e1 and e2 are on the same connected component.
  DHalfedge   *p_prev1 = _halfedge (prev1);
  DHalfedge   *p_prev2 = _halfedge (prev2);
  DInner_ccb  *hole1 = (p_prev1->is_on_inner_ccb()) ? p_prev1->inner_ccb() :
                                                      NULL;
  DInner_ccb  *hole2 = (p_prev2->is_on_inner_ccb()) ? p_prev2->inner_ccb() :
                                                      NULL;
  bool         prev1_before_prev2 = true;

  if (hole1 == hole2 && hole1 != NULL)
  {
    // If prev1 and prev2 are on different components, the insertion of the
    // new curve does not generate a new face, so the way we send these
    // halfedge pointers to the auxiliary function _insert_at_vertices() does
    // not matter.
    // However, in our case, where the two halfedges are reachable from one
    // another and are located on the same hole, a new face will be created
    // and form a hole inside their current incident face. In this case we
    // have to arrange prev1 and prev2 so that the new face (hole) will be
    // incident to the correct halfedge, directed from prev1's target to
    // prev2's target.
    // To do this, we use the topology traits to determine whether prev1 lies
    // inside the new face we are about to create (or alternatively, whether
    // prev2 does not lie inside this new face).

    Comparison_result path_res = 
      _compare_induced_path_length(p_prev1, p_prev2);

    prev1_before_prev2 = (path_res == LARGER) ?
      (_is_inside_new_face (p_prev1, p_prev2, cv)) :
      (! _is_inside_new_face (p_prev2, p_prev1, cv));
  }

  // We already have the comparsion result of the target vertices of prev1
  // and prev2. If however we swap the order of the halfedges, we take the
  // opposite comparison result.
  if (! prev1_before_prev2)
    res = CGAL::opposite (res);

  // Perform the insertion.
  bool        new_face_created = false;
  DHalfedge  *new_he = (prev1_before_prev2) ?
    _insert_at_vertices (cv, p_prev1, p_prev2, res, new_face_created) :
    _insert_at_vertices (cv, p_prev2, p_prev1, res, new_face_created);

  if (new_face_created)
  {
    // In case a new face has been created (pointed by the new halfedge we
    // obtained), we have to examine the holes and isolated vertices in the
    // existing face (pointed by the twin halfedge) and move the relevant
    // holes and isolated vertices into the new face.
    _relocate_in_new_face (new_he);
  }

  // Return a handle to the new halfedge directed from prev1's target to
  // prev2's target. Note that this may be the twin halfedge of the one
  // returned by _insert_at_vertices();
  if (! prev1_before_prev2)
    new_he = new_he->opposite();

  return (Halfedge_handle (new_he));
}

//-----------------------------------------------------------------------------
// Replace the point associated with the given vertex.
//
template<class GeomTraits, class TopTraits>
typename Arrangement_on_surface_2<GeomTraits, TopTraits>::Vertex_handle
Arrangement_on_surface_2<GeomTraits, TopTraits>::
modify_vertex(Vertex_handle vh, const Point_2& p)
{
  CGAL_precondition_msg
    (! vh->is_at_open_boundary(),
     "The modified vertex must not lie on open boundary.");
  CGAL_precondition_msg (m_geom_traits->equal_2_object() (vh->point(), p),
                         "The new point is different from the current one.");

  // Modify the vertex.
  _modify_vertex (_vertex (vh), p);

  // Return a handle to the modified vertex.
  return (vh);
}

//-----------------------------------------------------------------------------
// Remove an isolated vertex from the interior of a given face.
//
template<class GeomTraits, class TopTraits>
typename Arrangement_on_surface_2<GeomTraits, TopTraits>::Face_handle
Arrangement_on_surface_2<GeomTraits, TopTraits>::
remove_isolated_vertex(Vertex_handle v)
{
  CGAL_precondition (v->is_isolated());

  // Get the face containing v.
  DVertex      *p_v = _vertex (v);
  DIso_vertex  *iv = p_v->isolated_vertex();
  DFace        *p_f = iv->face();
  Face_handle   f = Face_handle (p_f);
  
  // Notify the observers that we are abount to remove a vertex.
  _notify_before_remove_vertex (v);

  // Remove the isolated vertex from the face that contains it.
  p_f->erase_isolated_vertex (iv);
  _dcel().delete_isolated_vertex (iv);

  // Delete the vertex.
  _delete_point (p_v->point());
  _dcel().delete_vertex (p_v);

  // Notify the observers that the vertex has been removed.
  _notify_after_remove_vertex ();

  // Return a handle for the face that used to contain the deleted vertex.
  return (f);
}

//-----------------------------------------------------------------------------
// Replace the x-monotone curve associated with the given edge.
//
template<class GeomTraits, class TopTraits>
typename Arrangement_on_surface_2<GeomTraits, TopTraits>::Halfedge_handle
Arrangement_on_surface_2<GeomTraits, TopTraits>::
modify_edge(Halfedge_handle e, const X_monotone_curve_2& cv)
{
  CGAL_precondition_msg (! e->is_fictitious(),
                         "The edge must be a valid one.");
  CGAL_precondition_msg (m_geom_traits->equal_2_object() (e->curve(), cv),
                         "The new curve is different from the current one.");

  // Modify the edge.
  _modify_edge (_halfedge (e), cv);

  // Return a handle to the modified halfedge.
  return (e);
}

//-----------------------------------------------------------------------------
// Split a given edge into two, and associate the given x-monotone
// curves with the split edges.
//
template<class GeomTraits, class TopTraits>
typename Arrangement_on_surface_2<GeomTraits, TopTraits>::Halfedge_handle
Arrangement_on_surface_2<GeomTraits, TopTraits>::
split_edge(Halfedge_handle e,
           const X_monotone_curve_2& cv1, const X_monotone_curve_2& cv2)
{
  CGAL_precondition_msg (! e->is_fictitious(),
                         "The edge must be a valid one.");

  // Get the split halfedge and its twin, its source and target.
  DHalfedge * he1 = _halfedge(e);
  DHalfedge * he2 = he1->opposite();
  DVertex   * source = he2->vertex();
  CGAL_precondition_code(DVertex * target = he1->vertex(););

  // Determine the point where we split the halfedge. We also determine which
  // curve should be associated with he1 (and he2), which is the curve who
  // has an endpoint that equals e's source, and which should be associated
  // with the new pair of halfedges we are about to split (the one who has
  // an endpoint which equals e's target).
  if (m_geom_traits->parameter_space_in_x_2_object()(cv1, ARR_MAX_END) ==
      ARR_INTERIOR &&
      m_geom_traits->parameter_space_in_y_2_object()(cv1, ARR_MAX_END) ==
      ARR_INTERIOR)
  {
    const Point_2 & cv1_right = 
      m_geom_traits->construct_max_vertex_2_object() (cv1);

    if (m_geom_traits->parameter_space_in_x_2_object()(cv2, ARR_MIN_END) ==
        ARR_INTERIOR &&
        m_geom_traits->parameter_space_in_y_2_object()(cv2, ARR_MIN_END) ==
        ARR_INTERIOR &&
        m_geom_traits->equal_2_object()(m_geom_traits->
                                        construct_min_vertex_2_object()(cv2),
                                        cv1_right))
    {
      // cv1's right endpoint and cv2's left endpoint are equal, so this should
      // be the split point. Now we check whether cv1 is incident to e's source
      // and cv2 to its target, or vice versa.
      if (_are_equal (source, cv1, ARR_MIN_END)) {
        CGAL_precondition_msg
          (_are_equal (target, cv2, ARR_MAX_END),
           "The subcurve endpoints must match e's end vertices.");

        return (Halfedge_handle (_split_edge (he1, cv1_right, cv1, cv2)));
      }

      CGAL_precondition_msg
        (_are_equal (source, cv2, ARR_MAX_END) &&
         _are_equal (target, cv1, ARR_MIN_END),
         "The subcurve endpoints must match e's end vertices.");

      return (Halfedge_handle (_split_edge (he1, cv1_right, cv2, cv1)));
    }
  }
  
  if (m_geom_traits->parameter_space_in_x_2_object()(cv1, ARR_MIN_END) ==
      ARR_INTERIOR &&
      m_geom_traits->parameter_space_in_y_2_object()(cv1, ARR_MIN_END) ==
      ARR_INTERIOR)
  {
    const Point_2 & cv1_left = 
      m_geom_traits->construct_min_vertex_2_object() (cv1);

    if (m_geom_traits->parameter_space_in_x_2_object()(cv2, ARR_MAX_END) ==
        ARR_INTERIOR &&
        m_geom_traits->parameter_space_in_y_2_object()(cv2, ARR_MAX_END) ==
        ARR_INTERIOR &&
        m_geom_traits->equal_2_object()(m_geom_traits->
                                        construct_max_vertex_2_object()(cv2),
                                        cv1_left))
    {
      // cv1's left endpoint and cv2's right endpoint are equal, so this should
      // be the split point. Now we check whether cv1 is incident to e's source
      // and cv2 to its target, or vice versa.
      if (_are_equal (source, cv2, ARR_MIN_END)) {
        CGAL_precondition_msg
          (_are_equal (target, cv1, ARR_MAX_END),
           "The subcurve endpoints must match e's end vertices.");

        return (Halfedge_handle (_split_edge (he1, cv1_left, cv2, cv1)));
      }

      CGAL_precondition_msg
        (_are_equal (source, cv1, ARR_MAX_END) && 
         _are_equal (target, cv2, ARR_MIN_END),
         "The subcurve endpoints must match e's end vertices.");

      return (Halfedge_handle (_split_edge (he1, cv1_left, cv1, cv2)));
    }
  }

  CGAL_error_msg ("The two subcurves must have a common endpoint.");
  return Halfedge_handle();
}

//-----------------------------------------------------------------------------
// Merge two edges to form a single edge, and associate the given x-monotone
// curve with the merged edge.
//
template<class GeomTraits, class TopTraits>
typename Arrangement_on_surface_2<GeomTraits, TopTraits>::Halfedge_handle
Arrangement_on_surface_2<GeomTraits, TopTraits>::
merge_edge(Halfedge_handle e1, Halfedge_handle e2,
           const X_monotone_curve_2& cv)
{
  CGAL_precondition_msg (! e1->is_fictitious() && ! e2->is_fictitious(),
                         "The edges must be a valid.");

  // Assign pointers to the existing halfedges, such that we have:
  //
  //            he1      he3
  //         -------> ------->
  //       (.)      (.)v     (.)
  //         <------- <-------
  //            he2      he4
  //
  DHalfedge   *_e1 = _halfedge (e1);
  DHalfedge   *_e2 = _halfedge (e2);
  DHalfedge   *he1 = NULL;
  DHalfedge   *he2 = NULL;
  DHalfedge   *he3 = NULL;
  DHalfedge   *he4 = NULL;

  if (_e1->vertex() == _e2->opposite()->vertex())
  {
    he1 = _e1;
    he2 = he1->opposite();
    he3 = _e2;
    he4 = he3->opposite();
  }
  else if (_e1->opposite()->vertex() == _e2->opposite()->vertex())
  {
    he2 = _e1;
    he1 = he2->opposite();
    he3 = _e2;
    he4 = he3->opposite();
  }
  else if (_e1->vertex() == _e2->vertex())
  {
    he1 = _e1;
    he2 = he1->opposite();
    he4 = _e2;
    he3 = he4->opposite();
  }
  else if (_e1->opposite()->vertex() == _e2->vertex())
  {
    he2 = _e1;
    he1 = he2->opposite();
    he4 = _e2;
    he3 = he4->opposite();
  }
  else
  {
    CGAL_precondition_msg (false,
                           "The input edges do not share a common vertex.");
    return Halfedge_handle();
  }

  // The vertex we are about to delete is now he1's target vertex.
  // Make sure that he1 and he4 are the only halfedges directed to v.
  DVertex    *v = he1->vertex();

  CGAL_precondition_msg
    (! v->has_null_point(),
     "The vertex removed by the merge must not lie on open boundary.");
  CGAL_precondition_msg
    (he1->next()->opposite() == he4 &&
     he4->next()->opposite() == he1,
     "The degree of the deleted vertex is greater than 2.");

  // Make sure the curve ends match the end vertices of the merged edge.
  CGAL_precondition_msg
    ((_are_equal (he2->vertex(), cv, ARR_MIN_END) && 
      _are_equal (he3->vertex(), cv, ARR_MAX_END)) ||
     (_are_equal (he3->vertex(), cv, ARR_MIN_END) && 
      _are_equal (he2->vertex(), cv, ARR_MAX_END)),
     "The endpoints of the merged curve must match the end vertices.");

  // Keep pointers to the components that contain two halfedges he3 and he2,
  // pointing at the end vertices of the merged halfedge.
  DInner_ccb  *ic1 = (he3->is_on_inner_ccb()) ? he3->inner_ccb() : NULL;
  DOuter_ccb  *oc1 = (ic1 == NULL) ? he3->outer_ccb() : NULL;

  DInner_ccb  *ic2 = (he4->is_on_inner_ccb()) ? he4->inner_ccb() : NULL;
  DOuter_ccb  *oc2 = (ic2 == NULL) ? he4->outer_ccb() : NULL;

  // Notify the observers that we are about to merge an edge.
  _notify_before_merge_edge (e1, e2, cv);

  // As he1 and he2 will evetually represent the merged edge, while he3 and he4
  // will be deleted, check if the deleted halfedges are represantatives of a
  // the CCBs they belong to. If so, replace he3 by he1 and he4 by he2. Note
  // that as we just change the component representatives, we do not have to
  // notify the observers on the change.
  if (oc1 != NULL && oc1->halfedge() == he3)
    oc1->set_halfedge (he1);
  else if (ic1 != NULL && ic1->halfedge() == he3)
    ic1->set_halfedge (he1);

  if (oc2 != NULL && oc2->halfedge() == he4)
    oc2->set_halfedge (he2);
  else if (ic2 != NULL && ic2->halfedge() == he4)
    ic2->set_halfedge (he2);

  // If he3 is the incident halfedge to its target, replace it by he1.
  if (he3->vertex()->halfedge() == he3)
    he3->vertex()->set_halfedge (he1);

  // Disconnect he3 and he4 from the edge list.
  if (he3->next() == he4)
  {
    // he3 and he4 form an "antenna", so he1 and he2 must be connected
    // together.
    he1->set_next(he2);
  }
  else
  {
    he1->set_next (he3->next());
    he4->prev()->set_next (he2);
  }

  // Note that he1 (and its twin) is going to represent the merged edge while
  // he3 (and its twin) is going to be removed. We therefore associate the
  // merged curve with he1 and delete the curve associated with he3.
  he1->curve() = cv;
  _delete_curve (he3->curve());

  // Set the properties of the merged edge.
  he1->set_vertex (he3->vertex());

  // Notify the observers that we are about to delete a vertex.
  _notify_before_remove_vertex (Vertex_handle (v));

  // Delete the point associated with the merged vertex.
  _delete_point (v->point());

  // Delete the merged vertex.
  _dcel().delete_vertex (v);

  // Notify the observers that the vertex has been deleted.
  _notify_after_remove_vertex ();

  // Delete the redundant halfedge pair.
  _dcel().delete_edge (he3);

  // Create a handle for one of the merged halfedges.
  Halfedge_handle   hh (he1);

  // Notify the observers that the edge has been merge.
  _notify_after_merge_edge (hh);

  // Return a handle for one of the merged halfedges.
  return (hh);
}

//-----------------------------------------------------------------------------
// Remove an edge from the arrangement.
//
template<class GeomTraits, class TopTraits>
typename Arrangement_on_surface_2<GeomTraits, TopTraits>::Face_handle
Arrangement_on_surface_2<GeomTraits, TopTraits>::
remove_edge(Halfedge_handle e, bool remove_source, bool remove_target)
{
  CGAL_precondition_msg (! e->is_fictitious(),
                         "The edge must be a valid one.");

  DHalfedge   *he1 = _halfedge (e);
  DHalfedge   *he2 = he1->opposite();
  DInner_ccb  *ic1 = (he1->is_on_inner_ccb()) ? he1->inner_ccb() : NULL;
  DOuter_ccb  *oc1 = (ic1 == NULL) ? he1->outer_ccb() : NULL;
  DFace       *f1 = (ic1 != NULL) ? ic1->face() : oc1->face();
  DInner_ccb  *ic2 = (he2->is_on_inner_ccb()) ? he2->inner_ccb() : NULL;
  DOuter_ccb  *oc2 = (ic2 == NULL) ? he2->outer_ccb() : NULL;
  DFace       *f2 = (ic2 != NULL) ? ic2->face() : oc2->face();
  DFace       *f;

  if (f1 != f2 ||
      he1->next() == he2 || he2->next() == he1)
  {
    // Either the removal of he1 (and its twin halfedge) will cause the two
    // incident faces to merge, or these two halfedges form an "antenna".
    // In either case, it does not matter which halfedge we send to the
    // auxiliary function _remove_edge().
    f = _remove_edge (he1, remove_source, remove_target);
  }
  else
  {
    // In this case one of the following can happen: (a) a new hole will be
    // created by the removal of the edge (case 3.2.1 of the removal
    // procedure), or (b) an outer CCB will be split into two (case 3.2.2).
    // We begin by locating the leftmost vertex along the path from he1 to its
    // twin he2 and the leftmost vertex point along the path from the twin to
    // he1 (both paths do not include he1 and he2 themselves).
    bool                            is_perimetric1;
    bool                            at_open_bnd1;
    bool                            is_perimetric2;
    bool                            at_open_bnd2;

    std::pair<int, const DVertex*>  v_min1 =
      _find_leftmost_vertex_on_closed_loop (he1, is_perimetric1, at_open_bnd1);

    std::pair<int, const DVertex*>  v_min2 =
      _find_leftmost_vertex_on_closed_loop (he2, is_perimetric2, at_open_bnd2);

    if (! is_perimetric1 && ! is_perimetric2)
    {
      CGAL_assertion (! at_open_bnd1 || ! at_open_bnd2);

      // Both paths from he1 to he2 and back from he2 to he1 are not
      // perimetric, so case (a) occurs. As we want to determine which
      // halfedge points to the new hole that will be create (he1 or he2),
      // we have to compare the two leftmost vertices lexicographically,
      // first by the indices then by x and y. v_min2 lies to the left of
      // v_min1 if and only if he1 points at the hole we are about to create.
      //
      //         +---------------------+
      //         |                     |
      //         |   he1    +----+     |
      //         +--------->+    |     |
      //         |          +----+     |
      //         |      v_min1         |
      //         |                     |
      //  v_min2 +---------------------+
      //
      // Note that if one of the paths we have examined goes to open boundary
      // (and only of the paths may go to open boundary), 
      // then it is obvious that the other path becomes a hole in an 
      // open face.
      if (at_open_bnd2 ||
          (!at_open_bnd1 && (v_min1.first > v_min2.first ||
                             (v_min1.first == v_min2.first &&
                              m_geom_traits->compare_xy_2_object()
                              (v_min1.second->point(),
                               v_min2.second->point()) == LARGER))))
      {
        // he1 is directed to the new hole to be created.
        f = _remove_edge (he1, remove_source, remove_target);
      }
      else
      {
        CGAL_assertion (at_open_bnd1 ||
                        (v_min1.first < v_min2.first) ||
                        (v_min1.first == v_min2.first &&
                         m_geom_traits->compare_xy_2_object()
                         (v_min1.second->point(),
                          v_min2.second->point()) == SMALLER));

        // he2 is directed to the new hole to be created. As its source and
        // target are swapped with respect to the end-vertices of the given
        // halfedge e, we swap the roles of the two input flags.
        f = _remove_edge (he2, remove_target, remove_source);
      }
    }
    else if (! is_perimetric1 && is_perimetric2)
    {
      // Case (a) occurs and he1 is directed to the new hole that will be
      // created.
      f = _remove_edge (he1, remove_source, remove_target);
    }
    else if (is_perimetric1 && ! is_perimetric2)
    {
      // Case (a) occurs and he2 is directed to the new hole that will be
      // created. As its source and target are swapped with respect to the
      // end-vertices of the given halfedge e, we swap the roles of the two
      // input flags.
      f = _remove_edge (he2, remove_target, remove_source);
    }
    else
    {
      // If both paths are perimetric, then case (b) takes place:
      f = _remove_edge (he1, remove_source, remove_target);
    }
  }

  return (Face_handle (f));
}

//-----------------------------------------------------------------------------
// Protected member functions (for internal use).
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Locate the place for the given curve around the given vertex.
//
template<class GeomTraits, class TopTraits>
typename Arrangement_on_surface_2<GeomTraits, TopTraits>::DHalfedge*
Arrangement_on_surface_2<GeomTraits, TopTraits>::
_locate_around_vertex(DVertex *v,
                      const X_monotone_curve_2& cv, Arr_curve_end ind) const
{
  // Check if the given curve-end has boundary conditions.
  const Arr_parameter_space  ps_x =
    m_geom_traits->parameter_space_in_x_2_object()(cv, ind);
  const Arr_parameter_space  ps_y =
    m_geom_traits->parameter_space_in_y_2_object()(cv, ind);

  if (ps_x != ARR_INTERIOR || ps_y != ARR_INTERIOR)
  {
    // Use the topology-traits class to locate the predecessor halfedge for
    // cv around the given vertex.
    return (m_topol_traits.locate_around_boundary_vertex (v, cv, ind, ps_x, ps_y));
  }

  // In case of a non-boundary vertex, we look for the predecessor around v.
  // Get the first incident halfedge around v and the next halfedge.
  DHalfedge   *first = v->halfedge();
  DHalfedge   *curr = first;

  if (curr == NULL)
    return (NULL);

  DHalfedge   *next = curr->next()->opposite();

  // In case there is only one halfedge incident to v, return this halfedge.
  if (curr == next)
    return (curr);

  // Otherwise, we traverse the halfedges around v until we find the pair
  // of adjacent halfedges between which we should insert cv.
  typename Traits_adaptor_2::Is_between_cw_2  is_between_cw =
    m_geom_traits->is_between_cw_2_object();

  bool       eq_curr, eq_next;

  while (! is_between_cw (cv, (ind == ARR_MIN_END),
                          curr->curve(), 
                          (curr->direction() == ARR_RIGHT_TO_LEFT),
                          next->curve(), 
                          (next->direction() == ARR_RIGHT_TO_LEFT),
                          v->point(), eq_curr, eq_next))
  {
    // If cv equals one of the curves associated with the halfedges, it is
    // an illegal input curve, as it already exists in the arrangement.
    if (eq_curr || eq_next)
      return (NULL);

    // Move to the next pair of incident halfedges.
    curr = next;
    next = curr->next()->opposite();

    // If we completed a full traversal around v without locating the
    // place for cv, it follows that cv overlaps and existing curve.
    if (curr == first)
      return (NULL);
  }

  // Return the halfedge we have located.
  return (curr);
}

//-----------------------------------------------------------------------------
// Compute the distance (in halfedges) between two halfedges.
//
template<class GeomTraits, class TopTraits>
unsigned int
Arrangement_on_surface_2<GeomTraits, TopTraits>::
_halfedge_distance(const DHalfedge *e1, const DHalfedge *e2) const
{
  CGAL_precondition (e1 != e2);
  if (e1 == e2)
    return (0);

  // Traverse the halfedge chain from e1 until reaching e2.
  const DHalfedge   *curr = e1->next();
  unsigned int       dist = 1;

  while (curr != e2)
  {
    // If we have returned to e1, e2 is not reachable from e1.
    if (curr == e1)
    {
      CGAL_error();
      return (0);
    }

    curr = curr->next();
    dist++;
  }

  // We have located e2 along the boundary of e1's component - return the
  // distance (number of halfedges) between e1 and e2.
  return (dist);
}

//-----------------------------------------------------------------------------
//Compare the length of the induced paths from e1 to e2 and from e2 to e1.
// return SMALLER if e1 to e2 is shorter, EQUAL if paths lengths are equal, 
//  o/w LARGER
//
template<class GeomTraits, class TopTraits>
Comparison_result
Arrangement_on_surface_2<GeomTraits, TopTraits>::
_compare_induced_path_length(const DHalfedge *e1, const DHalfedge *e2) const
{
  CGAL_precondition (e1 != e2);
  if (e1 == e2)
    return EQUAL;

  // Traverse the halfedge chain from e1 until reaching e2.
  const DHalfedge   *curr1 = e1->next();
  // Traverse the halfedge chain from e2 until reaching e1.
  const DHalfedge   *curr2 = e2->next();

  while (curr1 != e2 && curr2 != e1)
  {
    // If we have returned to e1, e2 is not reachable from e1.
    if (curr1 == e1)
    {
      CGAL_error();
      return EQUAL;
    }
    
    // If we have returned to e2, e1 is not reachable from e2.
    if (curr2 == e2)
    {
      CGAL_error();
      return EQUAL;
    }

    curr1 = curr1->next();
    curr2 = curr2->next();
  }

  Comparison_result res;

  // Return SMALLER if e1 to e2 is shorter than e2 to e1,
  //  EQUAL if their lengths are equal, or LARGER if e2 to e1 is longer.
  if (curr1 == e2)
    res = (curr2 != e1) ? SMALLER : EQUAL;
  else
    res = LARGER;

  CGAL_postcondition_code (int dist1 = _halfedge_distance(e1,e2));
  CGAL_postcondition_code (int dist2 = _halfedge_distance(e2,e1));
  CGAL_postcondition (((dist1 < dist2) && (res == SMALLER)) ||
                      ((dist1 == dist2) && (res == EQUAL)) ||
                      ((dist1 > dist2) && (res == LARGER)));

  return res;
}

//-----------------------------------------------------------------------------
// Move a given outer CCB from one face to another.
//
template<class GeomTraits, class TopTraits>
void Arrangement_on_surface_2<GeomTraits, TopTraits>::
_move_outer_ccb (DFace *from_face, DFace *to_face, DHalfedge *he)
{
  // Get the DCEL record that represents the outer CCB.
  DOuter_ccb  *oc = he->outer_ccb();

  CGAL_assertion (oc->face() == from_face);

  // Notify the observers that we are about to move an outer CCB.
  Ccb_halfedge_circulator   circ = (Halfedge_handle(he))->ccb();

  _notify_before_move_outer_ccb (Face_handle (from_face),
                                 Face_handle (to_face),
                                 circ);

  // Remove the hole from the current face.
  from_face->erase_outer_ccb (oc);

  // Modify the component that represents the hole.
  oc->set_face (to_face);
  to_face->add_outer_ccb (oc, he);
  
  // Notify the observers that we have moved the outer CCB.
  _notify_after_move_outer_ccb (circ);

  return;
}

//-----------------------------------------------------------------------------
// Move a given inner CCB (hole) from one face to another.
//
template<class GeomTraits, class TopTraits>
void Arrangement_on_surface_2<GeomTraits, TopTraits>::
_move_inner_ccb (DFace *from_face, DFace *to_face, DHalfedge *he)
{
  // Get the DCEL record that represents the inner CCB.
  DInner_ccb  *ic = he->inner_ccb();

  CGAL_assertion (ic->face() == from_face);

  // Notify the observers that we are about to move an inner CCB.
  Ccb_halfedge_circulator   circ = (Halfedge_handle(he))->ccb();

  _notify_before_move_inner_ccb (Face_handle (from_face),
                                 Face_handle (to_face),
                                 circ);

  // Remove the hole from the current face.
  from_face->erase_inner_ccb (ic);

  // Modify the component that represents the hole.
  ic->set_face (to_face);
  to_face->add_inner_ccb (ic, he);
  
  // Notify the observers that we have moved the inner CCB.
  _notify_after_move_inner_ccb (circ);

  return;
}

//-----------------------------------------------------------------------------
// Insert the given vertex as an isolated vertex inside the given face.
//
template<class GeomTraits, class TopTraits>
void Arrangement_on_surface_2<GeomTraits, TopTraits>::
_insert_isolated_vertex (DFace *f, DVertex *v)
{
#if CGAL_ARRANGEMENT_ON_SURFACE_INSERT_VERBOSE
  std::cout << "Aos_2: _insert_isolated_verteex (internal)" << std::endl;
  if (!v->has_null_point()) {
      std::cout << "v->point: " << v->point() << std::endl;
  }
  std::cout << "face   : " << f << std::endl;
#endif
  
  Face_handle     fh (f);
  Vertex_handle   vh (v);

  // Notify the observers that we are about to insert an isolated vertex
  // inside f.
  _notify_before_add_isolated_vertex (fh, vh);

  // Create an isolated vertex-information object,
  DIso_vertex    *iv = _dcel().new_isolated_vertex();

  // Set a pointer to the face containing the vertex.
  iv->set_face (f);

  // Initiate a new hole inside the given face.
  f->add_isolated_vertex (iv, v);
  
  // Associate the information with the vertex.
  v->set_isolated_vertex (iv);

  // Notify the observers that we have formed a new isolated vertex.
  _notify_after_add_isolated_vertex (vh);

  return;
}

//-----------------------------------------------------------------------------
// Move a given isolated vertex from one face to another.
//
template<class GeomTraits, class TopTraits>
void Arrangement_on_surface_2<GeomTraits, TopTraits>::
_move_isolated_vertex(DFace *from_face, DFace *to_face, DVertex *v)
{
  // Get the DCEL isolated-vertex record.
  DIso_vertex  *iv = v->isolated_vertex();

  // Notify the observers that we are about to move an isolated vertex.
  Vertex_handle   vh (v);

  _notify_before_move_isolated_vertex (Face_handle (from_face),
                                       Face_handle (to_face),
                                       vh);

  // Set the new face is the isolated vertex-information object.
  iv->set_face (to_face);

  // Move the isolated vertex from the first face to the other.
  from_face->erase_isolated_vertex (iv);
  to_face->add_isolated_vertex (iv, v);

  // Notify the observers that we have moved the isolated vertex.
  _notify_after_move_isolated_vertex (vh);

  return;
}

//-----------------------------------------------------------------------------
// Create a new vertex and associate it with the given point.
//
template<class GeomTraits, class TopTraits>
typename Arrangement_on_surface_2<GeomTraits, TopTraits>::DVertex*
Arrangement_on_surface_2<GeomTraits, TopTraits>::
_create_vertex (const Point_2& p)
{
  // Notify the observers that we are about to create a new vertex.
  Point_2  *p_p = _new_point (p);

  _notify_before_create_vertex (*p_p);

  // Create a new vertex and associate it with the given point.
  DVertex         *v = _dcel().new_vertex();

  v->set_point (p_p);
  v->set_boundary (ARR_INTERIOR, ARR_INTERIOR);

  // Notify the observers that we have just created a new vertex.
  Vertex_handle   vh (v);
  _notify_after_create_vertex (vh);

  return (v);
}

//-----------------------------------------------------------------------------
// Create a new vertex on boundary
//
template<class GeomTraits, class TopTraits>
typename Arrangement_on_surface_2<GeomTraits, TopTraits>::DVertex*
Arrangement_on_surface_2<GeomTraits, TopTraits>::
_create_boundary_vertex(const X_monotone_curve_2& cv, Arr_curve_end ind,
                        Arr_parameter_space ps_x, Arr_parameter_space ps_y)
{
  CGAL_precondition (ps_x != ARR_INTERIOR || ps_y != ARR_INTERIOR);

  // Notify the observers that we are about to create a new boundary vertex.
  _notify_before_create_boundary_vertex (cv, ind, ps_x, ps_y);

  // Create a new vertex and set its boundary conditions.
  DVertex         *v = _dcel().new_vertex();

  v->set_boundary (ps_x, ps_y);

  // Act according to the boundary type if there is one:
  if (is_open (ps_x, ps_y))
  {
    // The curve-end lies on open boundary so the vertex is not associated with
    // a valid point.
    v->set_point (NULL);
  }
  else
  {
    // Create a boundary vertex associated with a valid point.
    Point_2 * p_p = (ind == ARR_MIN_END) ?
      _new_point (m_geom_traits->construct_min_vertex_2_object()(cv)) :
      _new_point (m_geom_traits->construct_max_vertex_2_object()(cv));

    v->set_point (p_p);
  }

  // Notify the observers that we have just created a new boundary vertex.
  Vertex_handle   vh (v);
  _notify_after_create_boundary_vertex (vh);

  return (v);
}

//-----------------------------------------------------------------------------
// Locate the DCEL features that will be used for inserting the given curve
// end, which has a boundary condition, and set the proper vertex there.
//
template<class GeomTraits, class TopTraits>
typename Arrangement_on_surface_2<GeomTraits, TopTraits>::DVertex* 
Arrangement_on_surface_2<GeomTraits, TopTraits>::
_place_and_set_curve_end(DFace *f,
                         const X_monotone_curve_2& cv, Arr_curve_end ind,
                         Arr_parameter_space ps_x, Arr_parameter_space ps_y,
                         DHalfedge **p_pred)
{
  // Use the topology traits to locate the DCEL feature that contains the
  // given curve end.
  CGAL::Object obj = m_topol_traits.place_boundary_vertex (f, cv, ind, ps_x, ps_y);
  DVertex     *v;
  DHalfedge   *fict_he;

  // Act according to the result type.
  if (CGAL::assign (fict_he, obj))
  {
    // The curve end is located on a fictitious edge. We first create a new
    // vertex that corresponds to the curve end.
    v = _create_boundary_vertex (cv, ind, ps_x, ps_y);

    // Split the fictitious halfedge at the newly created vertex.
    // The returned halfedge is the predecessor for the insertion of the curve
    // end around v.
    _notify_before_split_fictitious_edge (Halfedge_handle (fict_he),
                                          Vertex_handle (v));

    *p_pred = m_topol_traits.split_fictitious_edge (fict_he, v);

    _notify_after_split_fictitious_edge (Halfedge_handle (*p_pred),
                                         Halfedge_handle ((*p_pred)->next()));
  }
  else if (CGAL::assign (v, obj))
  {
    // In this case we are given an existing vertex that represents the curve
    // end. We now have to locate the predecessor edge for the insertion of cv
    // around this vertex.
    *p_pred = m_topol_traits.locate_around_boundary_vertex(v, cv, ind, ps_x, ps_y);
  }
  else
  {
    CGAL_assertion (obj.is_empty());

    // In this case we have to create a new vertex that reprsents the given
    // curve end.
    v = _create_boundary_vertex (cv, ind, ps_x, ps_y);

    // Notify the topology traits on the creation of the boundary vertex.
    m_topol_traits.notify_on_boundary_vertex_creation (v, cv, ind, ps_x, ps_y);

    // There are no edges incident to v, therefore no predecessor halfedge.
    *p_pred = NULL;
  }

  // Return the vertex that represents the curve end.
  return (v);
}

//-----------------------------------------------------------------------------
// Insert an x-monotone curve into the arrangement, such that both its
// endpoints correspond to free arrangement vertices (newly created vertices
// or existing isolated vertices), so a new inner CCB is formed in the face
// that contains the two vertices.
//
template<class GeomTraits, class TopTraits>
typename Arrangement_on_surface_2<GeomTraits, TopTraits>::DHalfedge*
Arrangement_on_surface_2<GeomTraits, TopTraits>::
_insert_in_face_interior(const X_monotone_curve_2& cv,
                         DFace *f, DVertex *v1, DVertex *v2,
                         Comparison_result res)
{
#if CGAL_ARRANGEMENT_ON_SURFACE_INSERT_VERBOSE
  std::cout << "Aos_2: _insert_in_face_interior (internal)" << std::endl;
  std::cout << "cv   : " << cv << std::endl;
  std::cout << "face   : " << f << std::endl;
  if (!v1->has_null_point()) {
    std::cout << "v1->point: " << v1->point() << std::endl;
  }
  if (!v2->has_null_point()) {
    std::cout << "v2->point: " << v2->point() << std::endl;
  }
#endif

  // Notify the observers that we are about to create a new edge.
  _notify_before_create_edge (cv, Vertex_handle (v1), Vertex_handle (v2));

  // Create a pair of twin halfedges connecting the two vertices,
  // and link them together to form a new connected component, a hole in f.
  DHalfedge           *he1 = _dcel().new_edge();
  DHalfedge           *he2 = he1->opposite();
  DInner_ccb          *ic = _dcel().new_inner_ccb();
  X_monotone_curve_2  *dup_cv = _new_curve (cv);

  ic->set_face (f);
  he1->set_curve (dup_cv);

  he1->set_next (he2);
  he1->set_vertex (v1);
  he1->set_inner_ccb (ic);

  he2->set_next (he1);
  he2->set_vertex (v2);
  he2->set_inner_ccb (ic);

  // Assign the incident halfedges of the two new vertices.
  v1->set_halfedge (he1);
  v2->set_halfedge (he2);

  // Set the direction of the halfedges: res indicates the direction of he2,
  // as it is the comparison result between its source (v1) and target (v2).
  const Arr_halfedge_direction   dir =
    (res == SMALLER) ? ARR_LEFT_TO_RIGHT : ARR_RIGHT_TO_LEFT;
  he2->set_direction (dir);

  // Create a handle to the new halfedge pointing at the curve target.
  Halfedge_handle   hh (he2);

  // Notify the observers that we have created a new edge.
  _notify_after_create_edge (hh);

  // Notify the observers that we are about to form a new inner CCB inside f.
  _notify_before_add_inner_ccb (Face_handle (f), hh);

  // Initiate a new inner CCB inside the given face.
  f->add_inner_ccb (ic, he2);

  // Notify the observers that we have formed a new inner CCB.
  _notify_after_add_inner_ccb (hh->ccb());

  return (he2);
}

//-----------------------------------------------------------------------------
// Insert an x-monotone curve into the arrangement, such that one of its
// endpoints corresponds to a given arrangement vertex, given the exact
// place for the curve in the circular list around this vertex. The other
// endpoint corrsponds to a free vertex (a newly created vertex or an
// isolated vertex).
//
template<class GeomTraits, class TopTraits>
typename Arrangement_on_surface_2<GeomTraits, TopTraits>::DHalfedge*
Arrangement_on_surface_2<GeomTraits, TopTraits>::
_insert_from_vertex(const X_monotone_curve_2& cv,
                    DHalfedge *prev, DVertex *v,
                    Comparison_result cmp)
{
#if CGAL_ARRANGEMENT_ON_SURFACE_INSERT_VERBOSE
  std::cout << "Aos_2: _insert_from_vertex (internal)" << std::endl;
  std::cout << "cv   : " << cv << std::endl;
  if (!prev->has_null_curve()) {
    std::cout << "prev: " << prev->curve() << std::endl;
  } else {
    std::cout << "prev: fictitious" << std::endl;
  }
  std::cout << "pref: " << (prev->is_on_inner_ccb() ? 
                            prev->inner_ccb()->face() :
                            prev->outer_ccb()->face()) << std::endl;
  if (!v->has_null_point()) {
    std::cout << "v->point: " << v->point() << std::endl;
  }
  std::cout << "cmp  : " << cmp << std::endl;
#endif

  // Get the incident face of the previous halfedge. Note that this will also
  // be the incident face of the two new halfedges we are about to create.
  DInner_ccb  *ic = (prev->is_on_inner_ccb()) ? prev->inner_ccb() : NULL;
  DOuter_ccb  *oc = (ic == NULL) ? prev->outer_ccb() : NULL;

  // The first vertex is the one that the prev halfedge points to.
  // The second vertex is given by v.
  DVertex     *v1 = prev->vertex();
  DVertex     *v2 = v;

  // Notify the observers that we are about to create a new edge.
  _notify_before_create_edge (cv, Vertex_handle (v1), Vertex_handle (v2));

  // Create a pair of twin halfedges connecting the two vertices,
  // and associate them with the given curve.
  DHalfedge           *he1 = _dcel().new_edge();
  DHalfedge           *he2 = he1->opposite();
  X_monotone_curve_2  *dup_cv = _new_curve (cv);

  he1->set_curve (dup_cv);

  he1->set_vertex (v1);
  he2->set_vertex (v2);

  // Set the component for the new halfedge pair.
  if (oc != NULL)
  {
    // On an outer component:
    he1->set_outer_ccb (oc);
    he2->set_outer_ccb (oc);
  }
  else
  {
    // On an inner component:
    he1->set_inner_ccb (ic);
    he2->set_inner_ccb (ic);
  }

  // Associate the incident halfedge of the new vertex.
  v2->set_halfedge (he2);

  // Link the new halfedges around the existing vertex v1.
  he2->set_next (he1);
  he1->set_next (prev->next());

  prev->set_next (he2);

  // Set the direction of the halfedges: res indicates the direction of he2,
  // as it is the comparison result between its source and target (v).
  const Arr_halfedge_direction   dir =
    (cmp == SMALLER) ? ARR_LEFT_TO_RIGHT : ARR_RIGHT_TO_LEFT;
  he2->set_direction (dir);

  // Notify the observers that we have created a new edge.
  _notify_after_create_edge (Halfedge_handle (he2));

  // Return a pointer to the new halfedge whose target is the free vertex v.
  return (he2);
}

//-----------------------------------------------------------------------------
// Insert an x-monotone curve into the arrangement, where the end vertices
// are given by the target points of two given halfedges.
// The two halfedges should be given such that in case a new face is formed,
// it will be the incident face of the halfedge directed from the first
// vertex to the second vertex.
//
template<class GeomTraits, class TopTraits>
typename Arrangement_on_surface_2<GeomTraits, TopTraits>::DHalfedge*
Arrangement_on_surface_2<GeomTraits, TopTraits>::
_insert_at_vertices(const X_monotone_curve_2& cv,
                    DHalfedge *prev1, DHalfedge *prev2,
                    Comparison_result cmp,
                    bool& new_face)
{
  CGAL_precondition(prev1 != NULL);
  CGAL_precondition(prev2 != NULL);
  CGAL_precondition(prev1 != prev2); 

  // Get the vertices that match cv's endpoints.
  DVertex     *v1 = prev1->vertex();
  DVertex     *v2 = prev2->vertex();

#if CGAL_ARRANGEMENT_ON_SURFACE_INSERT_VERBOSE
  typedef internal::Sign_of_path< GeomTraits, TopTraits > Sign_of_path;
  
  std::cout << "Aos_2: _insert_at_vertices (internal)" << std::endl;
  
  std::cout << "cv   : " << cv << std::endl;
  if (!prev1->has_null_curve()) {
    std::cout << "prev1: " << prev1->curve() << std::endl;
  } else {
    std::cout << "prev1: fictitious" << std::endl;
  }
  std::cout << "dir1 : " << prev1->direction() << std::endl;
  std::cout << "pref: " << (prev1->is_on_inner_ccb() ? 
                            prev1->inner_ccb()->face() :
                            prev1->outer_ccb()->face()) << std::endl;
  if (!prev2->has_null_curve()) {
      std::cout << "prev2: " << prev2->curve() << std::endl;
  } else {
      std::cout << "prev2: fictitious" << std::endl;
  }
  std::cout << "dir 2: " << prev2->direction() << std::endl;
  std::cout << "pref2: " << (prev2->is_on_inner_ccb() ? 
                             prev2->inner_ccb()->face() :
                             prev2->outer_ccb()->face()) << std::endl;
  std::cout << "cmp  : " << cmp << std::endl;
#endif

  // Get the components containing the two previous halfedges and the incident
  // face (which should be the same for the two components).
  DInner_ccb  *ic1 = (prev1->is_on_inner_ccb()) ? prev1->inner_ccb() : NULL;
  DOuter_ccb  *oc1 = (ic1 == NULL) ? prev1->outer_ccb() : NULL;
  DFace       *f = (ic1 != NULL) ? ic1->face() : oc1->face();
  DInner_ccb  *ic2 = (prev2->is_on_inner_ccb()) ? prev2->inner_ccb() : NULL;
  DOuter_ccb  *oc2 = (ic2 == NULL) ? prev2->outer_ccb() : NULL;

  CGAL_precondition_code
    (DFace       *f2 = (ic2 != NULL) ? ic2->face() : oc2->face(););

#if CGAL_ARRANGEMENT_ON_SURFACE_INSERT_VERBOSE
  std::cout << "ic1: " << ic1 << std::endl;
  std::cout << "ic2: " << ic2 << std::endl;
  std::cout << "oc1: " << oc1 << std::endl;
  std::cout << "oc2: " << oc2 << std::endl;

  std::cout << "f1: " << &(*f) << std::endl;

#if 0
  DHalfedge *curr = prev1;
  if (curr != curr->next()) {
    curr = curr->next();
    while (curr != prev1) {
      if (!curr->has_null_curve()) {
        std::cout << "curr: " << curr->curve() << std::endl;
      } else {
        std::cout << "curr: fictitious" << std::endl;
      }
      std::cout << "dir: " 
                << (curr->direction() == CGAL::ARR_LEFT_TO_RIGHT ? "L2R" : "R2L") 
                << std::endl;
      curr = curr->next();
    }
  } else {
    std::cout << "only prev1" << std::endl;
  }
#endif

  CGAL_precondition_code(std::cout << "f2: " << &(*f2) << std::endl;);

#if 0
  curr = prev2;
  if (curr != curr->next()) {
    curr = curr->next();
    while (curr != prev2) {
      if (!curr->has_null_curve()) {
        std::cout << "curr: " << curr->curve() << std::endl;
      } else {
        std::cout << "curr: fictitious" << std::endl;
      }
      std::cout << "dir: " 
                << (curr->direction() == CGAL::ARR_LEFT_TO_RIGHT ? "L2R" : "R2L") 
                << std::endl;
      curr = curr->next();
    }
  } else {
    std::cout << "only prev2" << std::endl;
  }
#endif
#endif // CGAL_ARRANGEMENT_ON_SURFACE_INSERT_VERBOSE

  CGAL_precondition_msg
    (f == f2,
     "The two halfedges must share the same incident face.");

  // In case the two previous halfedges lie on the same inner component inside
  // the face f, we use the topology-traits class to determine whether we have
  // to create a new face by splitting f, and if so - whether new face is
  // contained in the existing face or just adjacent to it.
  bool         split_new_face = true;
  bool         is_split_face_contained = false;

  if (ic1 != NULL && ic1 == ic2)
  {
    std::pair<bool, bool>   res = 
        m_topol_traits.face_split_after_edge_insertion (prev1, prev2, cv);

    split_new_face = res.first;
    is_split_face_contained = res.second;

    // The result <false, true> is illegal:
    CGAL_assertion (split_new_face || ! is_split_face_contained);
  }

  // Notify the observers that we are about to create a new edge.
  _notify_before_create_edge (cv, Vertex_handle (v1), Vertex_handle (v2));

  // Create a pair of twin halfedges connecting v1 and v2 and associate them
  // with the given curve.
  DHalfedge           *he1 = _dcel().new_edge();
  DHalfedge           *he2 = he1->opposite();
  X_monotone_curve_2  *dup_cv = _new_curve (cv);

  he1->set_curve (dup_cv);

  he1->set_vertex (v1);
  he2->set_vertex (v2);

  // Connect the new halfedges to the existing halfedges around the two
  // incident vertices.
  he1->set_next (prev1->next());
  he2->set_next (prev2->next());

  prev1->set_next (he2);
  prev2->set_next (he1);

  // Set the direction of the halfedges: res indicates the direction of he2,
  // as it is the comparison result between its source and target.
  const Arr_halfedge_direction   dir =
    (cmp == SMALLER) ? ARR_LEFT_TO_RIGHT : ARR_RIGHT_TO_LEFT;
  he2->set_direction (dir);

  // Check the various cases of insertion (in the design document: the
  // various sub-cases of case 3 in the insertion procedure).
  if ((ic1 != NULL || ic2 != NULL) && ic1 != ic2)
  {
    // In case we have to connect two disconnected components, no new face
    // is created. 
    new_face = false;

    // Check whether both halfedges are inner components (holes) in the same
    // face, or whether one is a hole and the other is on the outer boundary
    // of the face. 
    Face_handle       fh (f);

    if (ic1 != NULL && ic2 != NULL)
    {
      // In this case (3.1) we have to connect to inner CCBs (holes) inside f.
      // Notify the observers that we are about to merge two holes in the face.
      _notify_before_merge_inner_ccb (fh,
                                      (Halfedge_handle(prev1))->ccb(),
                                      (Halfedge_handle(prev2))->ccb(),
                                      Halfedge_handle (he1));

      // Remove the inner component prev2 belongs to, and unite it with the
      // inner component that prev1 belongs to.
      f->erase_inner_ccb (ic2);

      // Set the merged component for the two new halfedges.
      he1->set_inner_ccb (ic1);
      he2->set_inner_ccb (ic1);

      // Make all halfedges along ic2 to point to ic1.
      DHalfedge       *curr;

      for (curr = he2->next(); curr != he1; curr = curr->next())
        curr->set_inner_ccb (ic1);

      // Delete the redundant inner CCB.
      _dcel().delete_inner_ccb (ic2);

      // Notify the observers that we have merged the two inner CCBs.
      _notify_after_merge_inner_ccb (fh, (Halfedge_handle(he1))->ccb());
    }
    else
    {
      // In this case (3.2) we connect a hole (inner CCB) with an outer CCB
      // of the face that contains it. We remove the hole and associate the
      // pair of new halfedges with the outer boundary of the face f.
      DInner_ccb  *del_ic;
      DOuter_ccb  *oc;
      DHalfedge   *ccb_first;
      DHalfedge   *ccb_last;
      
      if (ic1 != NULL)
      {
        // We remove the inner CCB ic1 and merge in with the outer CCB oc2.
        del_ic = ic1;
        oc = oc2;
        ccb_first = he1->next();
        ccb_last = he2;
      }
      else
      {
        // We remove the inner CCB ic2 and merge in with the outer CCB oc1.
        del_ic = ic2;
        oc = oc1;
        ccb_first = he2->next();
        ccb_last = he1;
      }

      he1->set_outer_ccb (oc);
      he2->set_outer_ccb (oc);

      // Notify the observers that we are about to remove an inner CCB from
      // the face.
      _notify_before_remove_inner_ccb(fh, (Halfedge_handle(ccb_first))->ccb());

      // Remove the inner CCB from the face, as we have just connected it to
      // the outer boundary of its incident face.
      f->erase_inner_ccb (del_ic);

      // Make all halfedges along the inner CCB to point to the outer CCB of f.
      DHalfedge       *curr;

      for (curr = ccb_first; curr != ccb_last; curr = curr->next())
        curr->set_outer_ccb (oc);

      // Delete the redundant hole.
      _dcel().delete_inner_ccb (del_ic);

      // Notify the observers that we have removed an inner ccb.
      _notify_after_remove_inner_ccb (fh);
    }
  }
  else if (! split_new_face)
  {
    // RWRW: NEW!
    CGAL_assertion (ic1 == ic2 && ic1 != NULL);

    // Handle the special case where we close an inner CCB, such that
    // we form two outer CCBs of the same face.
    Face_handle       fh (f);

    // Notify the obserers we are about to remove an inner CCB from f.
    _notify_before_remove_inner_ccb (fh, (Halfedge_handle(he1))->ccb());

    // Erase the inner CCB from the incident face and delete the
    // corresponding component.
    f->erase_inner_ccb (ic1);

    _dcel().delete_inner_ccb (ic1);

    // Notify the observers that the inner CCB has been removed.
    _notify_after_remove_inner_ccb (fh);

    // Handle the first split outer CCB (the one containing he1):
    // Notify the obserers we are about to add an outer CCB to f.
    _notify_before_add_outer_ccb (fh, Halfedge_handle (he1));

    // Create a new outer CCB that for the face f, and make he1 the
    // representative halfedge of this component.
    DOuter_ccb  *f_oc1 = _dcel().new_outer_ccb();

    f->add_outer_ccb (f_oc1, he1);
    f_oc1->set_face (f);
    he1->set_outer_ccb (f_oc1);

    // Set the component of all halfedges that used to belong to he1's CCB.
    DHalfedge       *curr;

    for (curr = he1->next(); curr != he1; curr = curr->next())
      curr->set_outer_ccb (f_oc1);

    // Notify the observers that we have added an outer CCB to f.
    _notify_after_add_outer_ccb ((Halfedge_handle (he1))->ccb());

    // Handle the second split outer CCB (the one containing he2):
    // Notify the obserers we are about to add an outer CCB to f.
    _notify_before_add_outer_ccb (fh, Halfedge_handle (he2));

    // Create a new outer CCB that for the face f, and make he2 the
    // representative halfedge of this component.
    DOuter_ccb  *f_oc2 = _dcel().new_outer_ccb();

    f->add_outer_ccb (f_oc2, he2);
    f_oc2->set_face (f);
    he2->set_outer_ccb (f_oc2);

    // Set the component of all halfedges that used to belong to he2's CCB.
    for (curr = he2->next(); curr != he2; curr = curr->next())
      curr->set_outer_ccb (f_oc2);

    // Notify the observers that we have added an outer CCB to f.
    _notify_after_add_outer_ccb ((Halfedge_handle (he2))->ccb());

    // Mark that in this case no new face is created:
    new_face = false;
  }
  else if (ic1 == ic2 && oc1 == oc2)
  {
    // In this case we created a pair of halfedge that connect halfedges that
    // already belong to the same component. This means we have to cretae a
    // new face by splitting the existing face f. 
    // Notify the observers that we are about to split a face.
    Face_handle       fh (f);

    _notify_before_split_face (fh, Halfedge_handle (he1));

    // Create the new face and create a single outer component which should
    // point to he2.
    DFace       *new_f = _dcel().new_face();
    //std::cout << "New face: " << &(*new_f) << std::endl;

#if CGAL_ARRANGEMENT_ON_SURFACE_INSERT_VERBOSE
    std::cout << "new face: " << new_f << std::endl;
#endif
    
    DOuter_ccb  *new_oc = _dcel().new_outer_ccb();

    new_face = true;
    new_f->add_outer_ccb (new_oc, he2);
    new_oc->set_face (new_f);

    // Set the components of the new halfedge he2, which should be the new
    // outer comoponent of the new face.
    // Note that there are several cases for setting he1's component, so we
    // do not do it yet.
    he2->set_outer_ccb (new_oc);

    // Set the component of all halfedges that used to belong to he2's CCB.
    DHalfedge       *curr;

    for (curr = he2->next(); curr != he2; curr = curr->next())
      curr->set_outer_ccb (new_oc);

#if CGAL_ARRANGEMENT_ON_SURFACE_INSERT_VERBOSE
    std::cout << "he2 (=> prev1) defines new outer CCB" << std::endl;
    std::cout << "prev1->face(): " << (prev1->is_on_inner_ccb() ? 
                                       prev1->inner_ccb()->face() :
                                       prev1->outer_ccb()->face())
              << std::endl;
    Sign_of_path sign_of_path(topology_traits());
    //std::cout << "prev1sign: " << sign_of_path(prev1, prev1) << std::endl;
#endif

    // Check whether the two previous halfedges lie on the same innder CCB
    // or on the same outer CCB (distinguish case 3.3 and case 3.4).
    bool   is_hole;

    if (ic1 != NULL)
    {
      // In this case (3.3) we have two distinguish two sub-cases.
      if (is_split_face_contained)
      {
        // The halfedges prev1 and prev2 belong to the same inner component
        // (hole) inside the face f, such that the new edge creates a new
        // face that is contained in f (case 3.3.1).
        is_hole = true;

        // In this case, he1 lies on an inner CCB of f.
        he1->set_inner_ccb (ic1);
        
        // Note that the current representative of the inner CCB may not
        // belong to the hole any more. In this case we replace the hole
        // representative by he1.
        if (! ic1->halfedge()->is_on_inner_ccb())
          ic1->set_halfedge (he1);

#if CGAL_ARRANGEMENT_ON_SURFACE_INSERT_VERBOSE
        std::cout << "he1 (=> prev2) defines new inner CCB" << std::endl;
        std::cout << "prev2->face(): " << (prev2->is_on_inner_ccb() ? 
                                           prev2->inner_ccb()->face() :
                                           prev2->outer_ccb()->face())
                  << std::endl;
        Sign_of_path sign_of_path(topology_traits());
        //std::cout << "prev2sign: " << sign_of_path(prev2, prev2) << std::endl;
#endif

      }
      else
      {
        // The new face we have created should be adjacent to the existing
        // face (case 3.3.2).
        is_hole = false;

        // Notify the obserers we are about to add an outer CCB to f.
        _notify_before_add_outer_ccb (fh, Halfedge_handle (he1));

        // Create a new outer CCB that for the face f, and make he1 the
        // representative halfedge of this component.
        DOuter_ccb  *f_oc = _dcel().new_outer_ccb();

        f->add_outer_ccb (f_oc, he1);
        f_oc->set_face (f);
        he1->set_outer_ccb (f_oc);

        // Set the component of all halfedges that used to belong to he1's
        // CCB.
        for (curr = he1->next(); curr != he1; curr = curr->next())
          curr->set_outer_ccb (f_oc);

#if CGAL_ARRANGEMENT_ON_SURFACE_INSERT_VERBOSE
        std::cout << "he1 (=> prev2) defines adjacent outer CCB" << std::endl;
        std::cout << "prev2->face(): " << (prev2->is_on_inner_ccb() ? 
                                           prev2->inner_ccb()->face() :
                                           prev2->outer_ccb()->face())
                  << std::endl;
        Sign_of_path sign_of_path(topology_traits());
        //std::cout << "prev2sign: " << sign_of_path(prev2, prev2) << std::endl;
#endif
        
        // Notify the observers that we have added an outer CCB to f.
        _notify_after_add_outer_ccb ((Halfedge_handle (he1))->ccb());

        // Go over all other outer CCBs of f and check whether they should be
        // moved to be outer CCBs of the new face.
        DOuter_ccb_iter  oc_it = f->outer_ccbs_begin();
        DOuter_ccb_iter  oc_to_move;

        while (oc_it != f->outer_ccbs_end())
        {
          // Use the topology traits to determine whether the representative
          // of the current outer CCB should belong to the same face as he2
          // (which is on the outer boundary of the new face).
          if (*oc_it != he1 &&
              m_topol_traits.boundaries_of_same_face (*oc_it, he2))
          {
            // We increment the itrator before moving the outer CCB, because
            // this operation invalidates the iterator.
            oc_to_move = oc_it;
            ++oc_it;

            _move_outer_ccb (f, new_f, *oc_to_move);
          }
          else
          {
            ++oc_it;
          }
        }
      }
    }
    else
    {
      // In this case the face f is simply split into two (case 3.4).
      is_hole = false;

      // In this case, he1 lies on an outer CCB of f.
      he1->set_outer_ccb (oc1);

      // As the outer component of the exisitng face f may associated with
      // one of the halfedges along the boundary of the new face, we set it
      // to be he1.
      oc1->set_halfedge (he1);
    }

    // Check whether we should mark the original face and the new face as
    // bounded or as unbounded faces.
    if (! f->is_unbounded())
    {
      // The original face is bounded, so the new face split from it is
      // obviously bounded.
      new_f->set_unbounded (false);
    }
    else if (is_hole)
    {
      // The new face is a hole inside the original face, so it must be
      // bounded.
      new_f->set_unbounded (false);
    }
    else
    {
      // Use the topology traits to determine whether each of the split
      // faces is unbounded. Note that if the new face is bounded, then f
      // obviously reamins unbounded and there is no need for further checks.
      new_f->set_unbounded (m_topol_traits.is_unbounded (new_f));
      
      if (new_f->is_unbounded())
        f->set_unbounded (m_topol_traits.is_unbounded (f));
    }

    // Notify the observers that we have split the face.
    _notify_after_split_face (fh, Face_handle (new_f), is_hole);
  }
  else
  {
    CGAL_assertion (oc1 != NULL && oc2 != NULL && oc1 != oc2);

    // In case prev1 and prev2 belong to different outer CCBs of the same
    // face f (case 3.5), we have to merge this ccbs into one. Note that we
    // do not create a new face.
    new_face = false;

    // Notify the observers that we are about to merge two outer CCBs.
    Face_handle       fh (f);

    _notify_before_merge_outer_ccb (fh,
                                    (Halfedge_handle(prev1))->ccb(),
                                    (Halfedge_handle(prev2))->ccb(),
                                    Halfedge_handle (he1));

    // Remove the outer component prev2 belongs to, and unite it with the
    // outer component that prev1 belongs to.
    f->erase_outer_ccb (oc2);

    // Set the merged component for the two new halfedges.
    he1->set_outer_ccb (oc1);
    he2->set_outer_ccb (oc1);

    // Make all halfedges along oc2 to point to oc1.
    DHalfedge       *curr;

    for (curr = he2->next(); curr != he1; curr = curr->next())
      curr->set_outer_ccb (oc1);

    // Delete the redundant outer CCB.
    _dcel().delete_outer_ccb (oc2);

    // Notify the observers that we have merged the two CCBs.
    _notify_after_merge_outer_ccb (fh,
                                   (Halfedge_handle(he1))->ccb());
  }

  // Notify the observers that we have created a new edge.
  _notify_after_create_edge (Halfedge_handle (he2));

  // Return the halfedge directed from v1 to v2.
  return (he2);
}

//-----------------------------------------------------------------------------
// Relocate all inner CCBs (holes) to their proper position,
// immediately after a face has split due to the insertion of a new halfedge.
//
template<class GeomTraits, class TopTraits>
void  Arrangement_on_surface_2<GeomTraits, TopTraits>::
_relocate_inner_ccbs_in_new_face (DHalfedge *new_he)
{
  // The given halfedge points to the new face, while its twin points to the
  // old face (the one that has just been split).
  DFace        *new_face = (new_he->is_on_inner_ccb()) ?
    new_he->inner_ccb()->face() :
    new_he->outer_ccb()->face();
  DHalfedge    *opp_he = new_he->opposite();
  const bool    opp_on_inner_ccb = opp_he->is_on_inner_ccb();
  DFace        *old_face = opp_on_inner_ccb ? opp_he->inner_ccb()->face() :
    opp_he->outer_ccb()->face();
  
  CGAL_assertion (new_face != old_face);
  
  // Examine the inner CCBs inside the existing old face and move the relevant
  // ones into the new face.
  DInner_ccb_iter    ic_it = old_face->inner_ccbs_begin();
  DInner_ccb_iter    ic_to_move;

  while (ic_it != old_face->inner_ccbs_end())
  {
    // In case the new edge represents the current component in the old face
    // (note we take the opposite halfedge, as it is incident to the old face),
    // then the new face already forms a hole in the old face, and we do not
    // need to move it.
    CGAL_assertion ((*ic_it)->is_on_inner_ccb());
    
    if (opp_on_inner_ccb && (*ic_it)->inner_ccb() == opp_he->inner_ccb())
    {
      ++ic_it;
      continue;
    }
    
    // Check whether the current inner CCB is inside new face (we actually
    // check if a representative vertex is located in the new face).
    if (m_topol_traits.is_in_face (new_face,
                               (*ic_it)->vertex()->point(),
                               (*ic_it)->vertex()))
    {
      // We increment the itrator before moving the inner CCB, because this
      // operation invalidates the iterator.
      ic_to_move = ic_it;
      ++ic_it;
      
      // Move the hole.
      _move_inner_ccb (old_face, new_face, *ic_to_move);
    }
    else
    {
      ++ic_it;
    }
  }
  
  return;
}

//-----------------------------------------------------------------------------
// Relocate all isolated vertices to their proper position,
// immediately after a face has split due to the insertion of a new halfedge.
//
template<class GeomTraits, class TopTraits>
void Arrangement_on_surface_2<GeomTraits, TopTraits>::
_relocate_isolated_vertices_in_new_face (DHalfedge *new_he)
{
  // The given halfedge points to the new face, while its twin points to the
  // old face (the one that has just been split).
  DFace        *new_face = (new_he->is_on_inner_ccb()) ?
    new_he->inner_ccb()->face() :
    new_he->outer_ccb()->face();
  DHalfedge    *opp_he = new_he->opposite();
  DFace        *old_face = (opp_he->is_on_inner_ccb()) ?
    opp_he->inner_ccb()->face() :
    opp_he->outer_ccb()->face();

  CGAL_assertion (new_face != old_face);

  // Examine the isolated vertices inside the existing old face and move the
  // relevant ones into the new face.
  DIso_vertex_iter    iv_it;
  DIso_vertex_iter    iv_to_move;

  iv_it = old_face->isolated_vertices_begin();
  while (iv_it != old_face->isolated_vertices_end())
  {
    // Check whether the isolated vertex lies inside the new face.
    if (m_topol_traits.is_in_face (new_face, iv_it->point(), &(*iv_it)))
    {
      // We increment the isolated vertices itrator before moving the vertex,
      // because this operation invalidates the iterator.
      iv_to_move  = iv_it;
      ++iv_it;

      // Move the isolated vertex.
      _move_isolated_vertex (old_face, new_face, &(*iv_to_move));
    }
    else
    {
      ++iv_it;
    }
  }

  return;
}

//-----------------------------------------------------------------------------
// Relocate all holes (inner CCBs) and isolated vertices to their proper
// position, immediately after a face has split due to the insertion of a new
// halfedge.
//
template<class GeomTraits, class TopTraits>
void Arrangement_on_surface_2<GeomTraits, TopTraits>::
_relocate_in_new_face(DHalfedge *new_he)
{
  _relocate_inner_ccbs_in_new_face(new_he);
  _relocate_isolated_vertices_in_new_face(new_he);
  return;
}

//-----------------------------------------------------------------------------
// Replace the point associated with the given vertex.
//
template<class GeomTraits, class TopTraits>
void Arrangement_on_surface_2<GeomTraits, TopTraits>::
_modify_vertex (DVertex *v, const Point_2& p)
{
  // Notify the observers that we are about to modify a vertex.
  Vertex_handle  vh (v);
  _notify_before_modify_vertex (vh, p);

  // Modify the point associated with the vertex.
  v->point() = p;

  // Notify the observers that we have modified the vertex.
  _notify_after_modify_vertex (vh);

  return;
}

//-----------------------------------------------------------------------------
// Replace the x-monotone curve associated with the given edge.
//
template<class GeomTraits, class TopTraits>
void Arrangement_on_surface_2<GeomTraits, TopTraits>::
_modify_edge (DHalfedge *he, const X_monotone_curve_2& cv)
{
  // Notify the observers that we are about to modify an edge.
  Halfedge_handle  e (he);
  _notify_before_modify_edge (e, cv);

  // Modify the curve associated with the halfedge.
  he->curve() = cv;

  // Notify the observers that we have modified the edge.
  _notify_after_modify_edge (e);

  return;
}

//-----------------------------------------------------------------------------
// Check if the given vertex represents one of the ends of a given curve.
//
template<class GeomTraits, class TopTraits>
bool Arrangement_on_surface_2<GeomTraits, TopTraits>::
_are_equal (const DVertex *v,
            const X_monotone_curve_2& cv, Arr_curve_end ind) const
{
  // In case the given curve end has boundary conditions, use the topology
  // traits to determine whether it is equivalent to v.
  const Arr_parameter_space    ps_x =
    m_geom_traits->parameter_space_in_x_2_object() (cv, ind);
  const Arr_parameter_space    ps_y =
    m_geom_traits->parameter_space_in_y_2_object() (cv, ind);

  if (ps_x != ARR_INTERIOR || ps_y != ARR_INTERIOR)
  {
    return (m_topol_traits.are_equal (v, cv, ind, ps_x, ps_y));
  }

  // Otherwise, the curve end is a valid endpoint. Check that v is also
  // associated with a valid point that equals this endpoint.
  if (v->has_null_point())
    return (false);

  return (ind == ARR_MIN_END) ?
    (m_geom_traits->equal_2_object() 
     (m_geom_traits->construct_min_vertex_2_object() (cv), v->point())) :
    (m_geom_traits->equal_2_object() 
     (m_geom_traits->construct_max_vertex_2_object() (cv), v->point()));
}

//-----------------------------------------------------------------------------
// Split a given edge into two at a given point, and associate the given
// x-monotone curves with the split edges.
//
template<class GeomTraits, class TopTraits>
typename Arrangement_on_surface_2<GeomTraits, TopTraits>::DHalfedge*
Arrangement_on_surface_2<GeomTraits, TopTraits>::
_split_edge(DHalfedge *e, const Point_2& p,
            const X_monotone_curve_2& cv1, const X_monotone_curve_2& cv2)
{
  // Allocate a new vertex and associate it with the split point.
  // Note that this point must not have any boundary conditions.
  DVertex         *v = _create_vertex (p);

  // Split the edge from the given vertex.
  return (_split_edge (e, v, cv1, cv2));
}

//-----------------------------------------------------------------------------
// Split a given edge into two at a given vertex, and associate the given
// x-monotone curves with the split edges.
//
template<class GeomTraits, class TopTraits>
typename Arrangement_on_surface_2<GeomTraits, TopTraits>::DHalfedge*
Arrangement_on_surface_2<GeomTraits, TopTraits>::
_split_edge(DHalfedge *e, DVertex *v,
            const X_monotone_curve_2& cv1, const X_monotone_curve_2& cv2)
{
  // Get the split halfedge and its twin, its source and target.
  DHalfedge       *he1 = e;
  DHalfedge       *he2 = he1->opposite();
  DInner_ccb      *ic1 = (he1->is_on_inner_ccb()) ? he1->inner_ccb() : NULL;
  DOuter_ccb      *oc1 = (ic1 == NULL) ? he1->outer_ccb() : NULL;
  DInner_ccb      *ic2 = (he2->is_on_inner_ccb()) ? he2->inner_ccb() : NULL;
  DOuter_ccb      *oc2 = (ic2 == NULL) ? he2->outer_ccb() : NULL;

  // Notify the observers that we are about to split an edge.
  _notify_before_split_edge (Halfedge_handle (e), Vertex_handle (v), cv1, cv2);

  // Allocate a pair of new halfedges.
  DHalfedge   *he3 = _dcel().new_edge();
  DHalfedge   *he4 = he3->opposite();

  // Connect the new halfedges:
  //
  //            he1      he3
  //         -------> ------->
  //       (.)      (.)v     (.)
  //         <------- <-------
  //            he2      he4
  //
  v->set_halfedge (he4);

  if (he1->next() != he2)
  {
    // Connect e3 between e1 and its successor.
    he3->set_next (he1->next());

    // Insert he4 between he2 and its predecessor.
    he2->prev()->set_next (he4);
  }
  else
  {
    // he1 and he2 form an "antenna", so he4 becomes he3's successor.
    he3->set_next (he4);
  }

  if (oc1 != NULL)
    he3->set_outer_ccb (oc1);
  else
    he3->set_inner_ccb (ic1);

  he3->set_vertex (he1->vertex());
  he4->set_vertex (v);
  he4->set_next (he2);

  if (oc2 != NULL)
    he4->set_outer_ccb (oc2);
  else
    he4->set_inner_ccb (ic2);

  if (he1->vertex()->halfedge() == he1)
    // If he1 is the incident halfedge to its target, he3 replaces it.
    he1->vertex()->set_halfedge (he3);

  // Update the properties of the twin halfedges we have just split.
  he1->set_next(he3);
  he1->set_vertex(v);

  // The direction of he3 is the same as he1's (and the direction of he4 is
  // the same as he2).
  he3->set_direction (he1->direction());

  // Associate cv1 with he1 (and its twin). We allocate a new curve for cv2
  // and associate it with he3 (and its twin).
  X_monotone_curve_2  *dup_cv2 = _new_curve (cv2);

  he1->curve() = cv1;
  he3->set_curve (dup_cv2);

  // Notify the observers that we have split an edge into two.
  _notify_after_split_edge (Halfedge_handle (he1), Halfedge_handle (he3));

  // Return a handle for one of the existing halfedge that is incident to the
  // split point.
  return (he1);
}

//-----------------------------------------------------------------------------
// Compare two vertices lexicographically, while taking care of boundary
// conditions (for the special usage of _find_leftmost_vertex() alone!).
//
template<class GeomTraits, class TopTraits>
Comparison_result
Arrangement_on_surface_2<GeomTraits, TopTraits>::
_compare_vertices_xy_impl (const DVertex * v1, const DVertex * v2,
                           Arr_not_all_sides_oblivious_tag) const
{
  if (v1 == v2)
    return (EQUAL);
  
  // Check the boundary conditions in y:
  const Arr_parameter_space     ps_y1 = v1->parameter_space_in_y();
  const Arr_parameter_space     ps_y2 = v2->parameter_space_in_y();

  // In case one of the vertices is a contraction point in y, then the
  // "negative" contraction is the smallest vertex, and the "positive"
  // contraction is considered to be the largest vertex.
  if (((ps_y1 == ARR_BOTTOM_BOUNDARY) &&
      is_contracted(Bottom_side_category())) ||
      ((ps_y2 == ARR_TOP_BOUNDARY) &&
       is_contracted(Top_side_category())))
    return SMALLER;

  if (((ps_y2 == ARR_BOTTOM_BOUNDARY) &&
      is_contracted(Bottom_side_category())) ||
      ((ps_y1 == ARR_TOP_BOUNDARY) &&
       is_contracted(Top_side_category())))
    return LARGER;

  // Check the boundary conditions in x:
  const Arr_parameter_space     ps_x1 = v1->parameter_space_in_x();
  const Arr_parameter_space     ps_x2 = v2->parameter_space_in_x();

  // A more elegant code:
  // if (ps_x1 == ps_x2) || (ps_y1 == ps_y1) compare ...
  //
  
  if (ps_x1 == ARR_LEFT_BOUNDARY) {
    return (ps_x1 == ps_x2) ?
      m_geom_traits->compare_xy_2_object() (v1->point(), v2->point()) : SMALLER;
  }
  else if (ps_x1 == ARR_RIGHT_BOUNDARY) {
    return (ps_x1 == ps_x2) ?
      m_geom_traits->compare_xy_2_object() (v1->point(), v2->point()) : LARGER;
  }

  if (ps_x2 == ARR_LEFT_BOUNDARY)
    return (LARGER);
  else if (ps_x2 == ARR_RIGHT_BOUNDARY)
    return (SMALLER);

  // Check the boundary conditions in y again:
  if (ps_y1 == ARR_BOTTOM_BOUNDARY) {
    return (ps_y1 == ps_y2) ?
      m_geom_traits->compare_xy_2_object() (v1->point(), v2->point()) : SMALLER;
  }
  else if (ps_y1 == ARR_TOP_BOUNDARY) {
    return (ps_y1 == ps_y2) ?
      m_geom_traits->compare_xy_2_object() (v1->point(), v2->point()) : LARGER;
  }

  if (ps_y2 == ARR_BOTTOM_BOUNDARY) return (LARGER);
  else if (ps_y2 == ARR_TOP_BOUNDARY) return (SMALLER);
  
  // If we reached here, both vertices do not have boundary conditions, and
  // we can just compare their associated points lexicographically.
  return (m_geom_traits->compare_xy_2_object() (v1->point(), v2->point()));
}

//-----------------------------------------------------------------------------
// Locate the leftmost vertex on the a given sequence defined by two
// halfedges. This sequence is still an open loop, but it will soon be closed
// by the insertion of the given curve: The first vertex we consider is the
// target vertex of the first halfedge, and the last vertex we consider is
// the source vertex of the second halfedge, such that the new curve will
// connect these two vertices.
//
//Note that the lowest halfedge incident to the leftmost vertex (he_left_low)
//is the lowest in the open loop. In the case one needs to have the lowest incident 
//halfedge in the interior of the closed loop, further handling must be done in the
//case he_before is he_left_low (cv cannot be used instead of he_left_low->curve()
//in that case as cv may not be defined to the right of the left lowest point).
template<class GeomTraits, class TopTraits>
std::pair<const typename Arrangement_on_surface_2<GeomTraits,
                                                  TopTraits>::DVertex*,
          const typename Arrangement_on_surface_2<GeomTraits,
                                                  TopTraits>::DHalfedge*>
Arrangement_on_surface_2<GeomTraits, TopTraits>::
_find_leftmost_vertex_on_open_loop (const DHalfedge *he_before,
                                    const DHalfedge *he_after,
                                    const X_monotone_curve_2& cv,
                                    bool& is_perimetric) const
{
  // We go over the sequence of vertices, starting from he_before's target
  // vertex, until reaching he_after's source vertex, and find the leftmost
  // one. Note that we do this carefully, keeping track of the number of
  // times we crossed the identification curve in x or in y (if they exist).
  // Note that the path must not be incident to any vertex on open boundary.
  typename Traits_adaptor_2::Parameter_space_in_x_2    parameter_space_in_x =
    m_geom_traits->parameter_space_in_x_2_object(); 
  typename Traits_adaptor_2::Parameter_space_in_y_2    parameter_space_in_y =
    m_geom_traits->parameter_space_in_y_2_object(); 
  typename Traits_adaptor_2::Compare_y_at_x_right_2    compare_y_at_x_right_2 =
    m_geom_traits->compare_y_at_x_right_2_object();  
  unsigned int      x_cross_count = 0;
  unsigned int      y_cross_count = 0;
  int               index = 0;
  const DHalfedge  *he = he_before;
  const DHalfedge  *he_left_low = NULL;
  int               ind_min = 0;
  const DVertex    *v_min = he_before->vertex();
  Arr_parameter_space     ps_x, ps_y;

  is_perimetric = false;
  do
  {
    // Get the boundary conditions of the current vertex.
    ps_x = he->vertex()->parameter_space_in_x();
    ps_y = he->vertex()->parameter_space_in_y();

    CGAL_assertion(!is_open(ps_x, ps_y));

    // Get the boundary conditions of the curve ends associated with the
    // current halfedge and its next halfedge.
    if (ps_x != ARR_INTERIOR || ps_y != ARR_INTERIOR)
    {
      // In case this he is the "before" halfegde, use the boundary conditions
      // of the proper curve-end of cv. Otherwise, use the boundary conditions
      // of the curve associated with he.
      if (he == he_before)
      {
        ps_x = parameter_space_in_x (cv, ARR_MIN_END);
        ps_y = parameter_space_in_y (cv, ARR_MIN_END);
        
        if ((ps_x == ARR_INTERIOR && ps_y == ARR_INTERIOR) ||
            ! m_topol_traits.are_equal (he->vertex(), cv, ARR_MIN_END, ps_x, ps_y))
        {
          ps_x = parameter_space_in_x (cv, ARR_MAX_END);
          ps_y = parameter_space_in_y (cv, ARR_MAX_END);
          
          CGAL_assertion (m_topol_traits.are_equal (he->vertex(),
                                                cv, ARR_MAX_END, ps_x, ps_y));
        }
      }
      else
      {
        if (he->direction() == ARR_RIGHT_TO_LEFT)
        {
          ps_x = parameter_space_in_x (he->curve(), ARR_MIN_END);          ps_y = parameter_space_in_y (he->curve(), ARR_MIN_END);
        }
        else
        {
          ps_x = parameter_space_in_x (he->curve(), ARR_MAX_END);
          ps_y = parameter_space_in_y (he->curve(), ARR_MAX_END);
        }
      }

      // In case this he->next() is the "after" halfegde, use the boundary
      // conditions of the proper curve-end of cv. Otherwise, use the boundary
      // conditions of the curve associated with he->next().
      Arr_parameter_space     ps_x_next, ps_y_next;

      if (he->next() == he_after)
      {
        ps_x_next = parameter_space_in_x (cv, ARR_MIN_END);
        ps_y_next = parameter_space_in_y (cv, ARR_MIN_END);
        
        if ((ps_x_next == ARR_INTERIOR && ps_y_next == ARR_INTERIOR) ||
            ! m_topol_traits.are_equal (he->next()->opposite()->vertex(),
                                    cv, ARR_MIN_END, ps_x_next, ps_y_next))
        {
          ps_x_next = parameter_space_in_x (cv, ARR_MAX_END);
          ps_y_next = parameter_space_in_y (cv, ARR_MAX_END);
          
          CGAL_assertion (m_topol_traits.are_equal
                          (he->next()->opposite()->vertex(),
                           cv, ARR_MAX_END, ps_x_next, ps_y_next));
        }
      }
      else
      {
        if (he->next()->direction() == ARR_LEFT_TO_RIGHT)
        {
          ps_x_next = parameter_space_in_x (he->next()->curve(), ARR_MIN_END);
          ps_y_next = parameter_space_in_y (he->next()->curve(), ARR_MIN_END);
        }
        else
        {
          ps_x_next = parameter_space_in_x (he->next()->curve(), ARR_MAX_END);
          ps_y_next = parameter_space_in_y (he->next()->curve(), ARR_MAX_END);
        }
      }

      // If we cross the identification curve in x, then we must update the
      // index. Note that a crossing takes place in the following cases:
      //                .                                  .
      //                .                                  .
      //                .                                  .
      //                . v    he                   he     . v
      //       <-------(.)<---------             -------->(.)------->
      //                .                                  .
      //       (BEFORE) .    (AFTER)              (BEFORE) .  (AFTER)
      //       index-1  .      index              index    .  index+1
      //
      if ((ps_x == ARR_LEFT_BOUNDARY) && (ps_x_next == ARR_RIGHT_BOUNDARY))
      {
        CGAL_assertion(is_identified(Left_side_category()) && 
                       is_identified(Right_side_category()));
        x_cross_count++;
        index--;
      }
      else if ((ps_x == ARR_RIGHT_BOUNDARY) && (ps_x_next == ARR_LEFT_BOUNDARY))
      {
        CGAL_assertion(is_identified(Left_side_category()) && 
                       is_identified(Right_side_category()));
        x_cross_count++;
        index++;
      }

      // Check if we cross the identification curve in y.
      if (((ps_y == ARR_BOTTOM_BOUNDARY) && (ps_y_next == ARR_TOP_BOUNDARY)) ||
          ((ps_y == ARR_TOP_BOUNDARY) && (ps_y_next == ARR_BOTTOM_BOUNDARY)))
      {
        CGAL_assertion(is_identified(Bottom_side_category()) &&
                       is_identified(Top_side_category()));
        y_cross_count++;
      }
    }

    // If the halfedge is directed from right to left, its target vertex is
    // smaller than its source, so we should check whether it is also smaller
    // than the leftmost vertex so far. Note that we compare the vertices
    // lexicographically: first by the indices, then by x and y.
    if (he->direction() == ARR_RIGHT_TO_LEFT &&
        (he->next() == he_after ||
         he->next()->direction() == ARR_LEFT_TO_RIGHT))
    {
      if (v_min == he->opposite()->vertex() || v_min == he->vertex() ||
          index < ind_min ||
          (index == ind_min &&
           _compare_vertices_xy (he->vertex(), v_min) == SMALLER))
      {
        ind_min = index;
        bool v_min_updated = v_min!=he->vertex();
        v_min = he->vertex();

        if ( v_min_updated || he_left_low ==NULL ||
             compare_y_at_x_right_2(he_left_low->curve(),
                                    he->curve(),
                                    v_min->point() ) == LARGER )
        {
          // If we need to compute the lowest halfedge incident to the leftmost
          // vertex, update it now. Note that we may visit the smallest vertex
          // several times (thus the compare_y_at_x_right_2).
          he_left_low = he;
        }
      }
    }
          
    // Move to the halfedge.
    he = he->next();

  } while (he != he_after);

  // Determine if the path is perimetric, namely if there exists an identification
  // curve in x (or in y), and we have crossed it an odd number of
  // times.
  is_perimetric = (x_cross_count % 2 == 1) || (y_cross_count % 2 == 1);
  
  // Return the leftmost vertex and its index (with respect to he_before).
  return (std::make_pair (v_min, (he_left_low==he_before ? static_cast<DHalfedge*>(NULL):he_left_low) ));
}

//-----------------------------------------------------------------------------
// Locate the leftmost vertex on the a given sequence, defined by an anchor
// halfedge and its twin, which forms a closed loop (i.e., the anchor's twin
// is reachable from the anchor halfedge).
//
template<class GeomTraits, class TopTraits>
std::pair<int,
          const typename Arrangement_on_surface_2<GeomTraits,
                                                  TopTraits>::DVertex*>
Arrangement_on_surface_2<GeomTraits, TopTraits>::
_find_leftmost_vertex_on_closed_loop (const DHalfedge *he_anchor,
                                      bool& is_perimetric,
                                      bool& at_open_bnd) const
{
  // We go over the sequence of vertices, starting from he_anchor's target
  // vertex and stopping at the source vertex of its twin. As this path is a
  // closed loop, he_anchor's twin is reachable from he_anchor.
  // Note that we do this carefully, keeping track of the number of times we
  // crossed the identification curve in x or in y (if they exist). We return
  // the leftmost vertex we find, along with its index with respect of the
  // he_anchor halfedge; this index is decremented each time we cross the 
  // identification curve from right to left, and
  // incremented each time we cross it from left to right.
  typename Traits_adaptor_2::Parameter_space_in_x_2    parameter_space_in_x =
    m_geom_traits->parameter_space_in_x_2_object(); 
  typename Traits_adaptor_2::Parameter_space_in_y_2    parameter_space_in_y =
    m_geom_traits->parameter_space_in_y_2_object(); 
  unsigned int      x_cross_count = 0;
  unsigned int      y_cross_count = 0;
  int               index = 0;
  const DHalfedge  *he = he_anchor;
  int               ind_min = 0;
  const DVertex    *v_min = he->vertex();
  Arr_parameter_space     ps_x, ps_y;

  is_perimetric = false;
  at_open_bnd = false;
  do
  {
    // Get the boundary conditions of the current vertex.
    ps_x = he->vertex()->parameter_space_in_x();
    ps_y = he->vertex()->parameter_space_in_y();

    // Get the boundary conditions of the curve ends associated with the
    // current halfedge and its next halfedge.
    if (ps_x != ARR_INTERIOR || ps_y != ARR_INTERIOR) {
      // Stop here if the current vertex lies on open boundary
      if (is_open (ps_x, ps_y)) {
        at_open_bnd = true;
        index = 0;
        v_min = NULL;
        return (std::make_pair (index, v_min));
      }

      // If we are on the anchor halfedge, use the boundary conditions of
      // the curve associated with the predecessor of the anchor's twin.
      const DHalfedge  *he_curr = he;

      if (he == he_anchor)
        he = he_anchor->opposite()->prev();

      CGAL_assertion (! he_curr->has_null_curve());
      if (he_curr->direction() == ARR_RIGHT_TO_LEFT)
      {
        ps_x = parameter_space_in_x (he_curr->curve(), ARR_MIN_END);
        ps_y = parameter_space_in_y (he_curr->curve(), ARR_MIN_END);
      }
      else
      {
        ps_x = parameter_space_in_x (he_curr->curve(), ARR_MAX_END);
        ps_y = parameter_space_in_y (he_curr->curve(), ARR_MAX_END);
      }

      // Get the boundary conditions of the curve-end of the next halfedge.
      Arr_parameter_space     ps_x_next, ps_y_next;

      CGAL_assertion (! he->next()->has_null_curve());
      if (he->next()->direction() == ARR_LEFT_TO_RIGHT)
      {
        ps_x_next = parameter_space_in_x (he->next()->curve(), ARR_MIN_END);
        ps_y_next = parameter_space_in_y (he->next()->curve(), ARR_MIN_END);
      }
      else
      {
        ps_x_next = parameter_space_in_x (he->next()->curve(), ARR_MAX_END);
        ps_y_next = parameter_space_in_y (he->next()->curve(), ARR_MAX_END);
      }

      // If we cross the identification curve in x, then we must update the
      // index. Note that a crossing takes place in the following cases:
      //                .                                  .
      //                .                                  .
      //                .                                  .
      //                . v    he                   he     . v
      //       <-------(.)<---------             -------->(.)------->
      //                .                                  .
      //       (BEFORE) .    (AFTER)              (BEFORE) .  (AFTER)
      //       index-1  .      index              index    .  index+1
      //
      if ((ps_x == ARR_LEFT_BOUNDARY) && (ps_x_next == ARR_RIGHT_BOUNDARY))
      {
        CGAL_assertion(is_identified(Left_side_category()) && 
                       is_identified(Right_side_category()));
        x_cross_count++;
        index--;
      }
      else if ((ps_x == ARR_RIGHT_BOUNDARY) && (ps_x_next == ARR_LEFT_BOUNDARY))
      {
        CGAL_assertion(is_identified(Left_side_category()) && 
                       is_identified(Right_side_category()));
        x_cross_count++;
        index++;
      }

      // Check if we cross the identification curve in y.
      if (((ps_y == ARR_BOTTOM_BOUNDARY) && (ps_y_next == ARR_TOP_BOUNDARY)) ||
          ((ps_y == ARR_TOP_BOUNDARY) && (ps_y_next == ARR_BOTTOM_BOUNDARY)))
      {
        CGAL_assertion(is_identified(Bottom_side_category()) &&
                       is_identified(Top_side_category()));
        y_cross_count++;
      }
    }

    // If the halfedge is directed from right to left, its target vertex is
    // smaller than its source, so we should check whether it is also smaller
    // than the leftmost vertex so far. Note that we compare the vertices
    // lexicographically: first by the indices, then by x and y.
    if (he != he_anchor && he->direction() == ARR_RIGHT_TO_LEFT &&
        he->next()->direction() == ARR_LEFT_TO_RIGHT)
    {
      if (v_min == he->opposite()->vertex() ||
          v_min == he->vertex() ||
          index < ind_min ||
          (index == ind_min &&
           _compare_vertices_xy (he->vertex(), v_min) == SMALLER))
      {
        ind_min = index;
        v_min = he->vertex();
      }
    }
          
    // Move to the next halfedge.
    he = he->next();
    CGAL_assertion (he != he_anchor);      // Guard for infinite loops.

  } while (he->next() != he_anchor->opposite());

  // Determine if the path is perimetric, namely if there exists an identification
  // curve in x (or in y), and we have crossed it an odd number of
  // times.
  is_perimetric = (x_cross_count % 2 == 1) || (y_cross_count % 2 == 1);
    
  // Return the leftmost vertex and its index (with respect to he_anchor).
  return (std::make_pair (ind_min, v_min));
}

//-----------------------------------------------------------------------------
// Determine whether a given query halfedge lies in the interior of a new
// face we are about to create, by connecting it with another halfedge
// using a given x-monotone curve.
//
template <class GeomTraits, class TopTraits>
bool Arrangement_on_surface_2<GeomTraits, TopTraits>::
_is_inside_new_face (const DHalfedge *prev1,
                     const DHalfedge *prev2,
                     const X_monotone_curve_2& cv) const
{
  // Go over all halfedges of along boundary of the face which will eventually
  // contain prev1: As the new face is not constructed yet, this traversal
  // is simulated by going from prev2's target vertex to prev1's target vertex
  // (the source of prev1's successor) and then returning over the new curve.
  // During the traversal we locate the leftmost halfedge along the boundary
  // (i.e, the one with the lexicographically smallest target vertex, which is
  // also the lowest halfedge incident to this vertex we encountered during
  // our traversal).
  const DHalfedge              *he_last = prev1->next();
  bool                          is_perimetric;
  std::pair<const DVertex*, const DHalfedge*>   find_res =
    _find_leftmost_vertex_on_open_loop (prev2, he_last, cv, is_perimetric);

  if (is_perimetric)
  {
    // std::cout << "perimetric" << std::endl;
    // In this case the route from prev1's target to prev2's target is
    // perimetric. We use the topology traits to determine which halfedge
    // lies inside the hole (in case a hole is indeed created).
    return (m_topol_traits.is_on_new_perimetric_face_boundary (prev1, prev2, cv));
  }

  const DVertex                *v_min = find_res.first;
  const DHalfedge              *he_left_low = find_res.second;

  CGAL_assertion (! v_min->has_null_point());

  // Now note that the curves of leftmost edge and its successor are defined
  // to the right of the smallest vertex. We compare them to the right of this
  // point to determine whether prev1 is inside the hole to be created or not.
  const X_monotone_curve_2  *p_cv_curr;
  const X_monotone_curve_2  *p_cv_next;

  if (he_left_low != NULL)
  {
    p_cv_curr = &(he_left_low->curve());

    // Take special care if the next curve should really be the new curve.
    if (he_left_low->next() != he_last)
      p_cv_next = &(he_left_low->next()->curve());
    else
      p_cv_next = &cv;
  }
  else
  {
    // In this case, the leftmost edge should be the one associated with the
    // new curve (which has not been created yet).
    p_cv_curr = &cv;
    p_cv_next = &(prev2->next()->curve());
  }

  // Check if the vertex lies on the identification curve in y, in which case
  // special care must be taken.
  if (((v_min->parameter_space_in_y() == ARR_BOTTOM_BOUNDARY) &&
       is_identified(Bottom_side_category())) ||
      ((v_min->parameter_space_in_y() == ARR_TOP_BOUNDARY) &&
       is_identified(Top_side_category())))
  {
    // Both current and next curves are incident to the identification curve.
    // As v_min is the leftmost vertex, we now that their left ends must have
    // a boundary condition of type identification in y.
    Arr_parameter_space  ps_y_curr =
      m_geom_traits->parameter_space_in_y_2_object()(*p_cv_curr, ARR_MIN_END);
    Arr_parameter_space  ps_y_next =
      m_geom_traits->parameter_space_in_y_2_object()(*p_cv_next, ARR_MIN_END);

    // Check if the curves lie on opposite sides of the identification curve.
    if ((ps_y_curr == ARR_BOTTOM_BOUNDARY) && (ps_y_next == ARR_TOP_BOUNDARY))
    {
      // In this case the current curve is "above" the next one to the right
      // of v_min, in a cyclic order around the identification curve.
      return (true);
    }
    else if ((ps_y_curr == ARR_TOP_BOUNDARY) &&
             (ps_y_next == ARR_BOTTOM_BOUNDARY))
    {
      // In this case the current curve is "below" the next one to the right
      // of v_min, in a cyclic order around the identification curve.
      return (false);
    }

    // If both curves are on the same side of the identification curve, we
    // continue to compare them to the right of v_min.
    CGAL_assertion (((ps_y_curr == ARR_BOTTOM_BOUNDARY) &&
                     (ps_y_next == ARR_BOTTOM_BOUNDARY)) ||
                    ((ps_y_curr == ARR_TOP_BOUNDARY) &&
                     (ps_y_next == ARR_TOP_BOUNDARY)));
  }

  // Check if the leftmost vertex is acontraction point in y, in which case
  // special care must be taken.
  if (((v_min->parameter_space_in_y() == ARR_BOTTOM_BOUNDARY) &&
       is_contracted(Bottom_side_category())) ||
      ((v_min->parameter_space_in_y() == ARR_TOP_BOUNDARY) &&
       is_contracted(Top_side_category())))
  {
    // Get the curve-ends for cv_curr and cv_next that conincide with the
    // contraction point.
    Arr_curve_end      ind_curr;
    Arr_curve_end      ind_next;

    if (he_left_low != NULL)
    {
      // The contraction point is he_left_low's target and cv_curr is its
      // associated curve.
      ind_curr = (he_left_low->direction() == ARR_LEFT_TO_RIGHT) ?
        ARR_MAX_END : ARR_MIN_END;

      if (he_left_low->next() != he_last)
      {
        // The contraction point is he_left_low->next()'s source and cv_next
        // is its associated curve.
        ind_next = (he_left_low->next()->direction() == ARR_LEFT_TO_RIGHT) ?
          ARR_MIN_END : ARR_MAX_END;
      }
      else {
        // In this case cv_next equals cv.
        ind_next =
          (m_topol_traits.are_equal
           (v_min, cv, ARR_MIN_END,
            m_geom_traits->parameter_space_in_x_2_object() (cv, ARR_MIN_END),
            m_geom_traits->parameter_space_in_y_2_object() (cv, ARR_MIN_END))) ?
          ARR_MIN_END : ARR_MAX_END;
      }
    }
    else
    {
      // In this case cv_curr equals cv.
      ind_curr =
        (m_topol_traits.are_equal
         (v_min, cv, ARR_MIN_END,
          m_geom_traits->parameter_space_in_x_2_object() (cv, ARR_MIN_END),
          m_geom_traits->parameter_space_in_y_2_object() (cv, ARR_MIN_END))) ?
        ARR_MIN_END : ARR_MAX_END;

      // The contraction point is prev2->next()'s source and cv_next
      // is its associated curve.
      ind_next = (prev2->next()->direction() == ARR_LEFT_TO_RIGHT) ?
        ARR_MIN_END : ARR_MAX_END;
    }

    // Compare the horizontal position of the two curve-ends at the point
    // of contraction.

    Comparison_result x_res = EQUAL;

    if (ind_curr != ind_next) {
      x_res = (ind_curr == ARR_MAX_END ? SMALLER : LARGER);
    } else {
      
      CGAL_assertion(ind_curr == ind_next);
      
      x_res =
	m_geom_traits->compare_x_curve_ends_2_object() (*p_cv_curr, ind_curr,
							*p_cv_next, ind_next);
      
    }
    CGAL_assertion (x_res != EQUAL);

    return (((v_min->parameter_space_in_y() == ARR_BOTTOM_BOUNDARY) &&
             (x_res == SMALLER)) ||
            ((v_min->parameter_space_in_y() == ARR_TOP_BOUNDARY) &&
             (x_res == LARGER)));
  }

  // TODO EBEB minimal point can also be on left boundary where we have to call boundary-functors

  return (m_geom_traits->compare_y_at_x_right_2_object()
          (*p_cv_curr, *p_cv_next, v_min->point()) == LARGER);
}

//-----------------------------------------------------------------------------
// Remove a pair of twin halfedges from the arrangement.
// In case the removal causes the creation of a new hole, the given halfedge
// should point at this hole.
//
template<class GeomTraits, class TopTraits>
typename Arrangement_on_surface_2<GeomTraits, TopTraits>::DFace*
Arrangement_on_surface_2<GeomTraits, TopTraits>::
_remove_edge (DHalfedge *e, bool remove_source, bool remove_target)
{
  // Get the pair of twin edges to be removed, the connected components they
  // belong to and their incident faces.
  DHalfedge   *he1 = e;
  DHalfedge   *he2 = e->opposite();
  DInner_ccb  *ic1 = (he1->is_on_inner_ccb()) ? he1->inner_ccb() : NULL;
  DOuter_ccb  *oc1 = (ic1 == NULL) ? he1->outer_ccb() : NULL;
  DFace       *f1 = (oc1 != NULL) ? oc1->face() : ic1->face();
  DInner_ccb  *ic2 = (he2->is_on_inner_ccb()) ? he2->inner_ccb() : NULL;
  DOuter_ccb  *oc2 = (ic2 == NULL) ? he2->outer_ccb() : NULL;
  DFace       *f2 = (oc2 != NULL) ? oc2->face() : ic2->face();
  DHalfedge   *prev1 = NULL;
  DHalfedge   *prev2 = NULL;

  // Notify the observers that we are about to remove an edge.
  Halfedge_handle  hh (e);

  _notify_before_remove_edge (hh);

  // Check if the two incident faces are equal, in which case no face will be
  // merged and deleted (and a hole may be created).
  if (f1 == f2)
  {
    // Check if the two halfedges are successors along the face boundary.
    if (he1->next() == he2 && he2->next() == he1)
    {
      CGAL_assertion (ic1 != NULL && ic1 == ic2);

      // The two halfedges form a "singleton" hole inside the incident face
      // (case 1 of the removal procedure, as detailed in the design document),
      // so we simply have to remove it.
      // First notify the observers that we are about to remove this hole
      // (inner CCB).
      Face_handle               fh (f1);

      _notify_before_remove_inner_ccb (fh, (Halfedge_handle(he1))->ccb());

      // Erase the inner CCB from the incident face and delete the
      // corresponding component.      
      f1->erase_inner_ccb (ic1);

      _dcel().delete_inner_ccb (ic1);

      // Notify the observers that the inner CCB has been removed.
      _notify_after_remove_inner_ccb (fh);

      // Remove the end-vertices, if necessary.
      if (remove_target)
      {
        if (he1->vertex()->parameter_space_in_x() != ARR_INTERIOR || 
            he1->vertex()->parameter_space_in_y() != ARR_INTERIOR)
        {
          he1->vertex()->set_halfedge(NULL);    // disconnect the end vertex
          _remove_vertex_if_redundant (he1->vertex(), f1);
        }
        else
        {
          // Delete the he1's target vertex and its associated point.
          _notify_before_remove_vertex (Vertex_handle (he1->vertex()));
          
          _delete_point (he1->vertex()->point());
          _dcel().delete_vertex (he1->vertex());

          _notify_after_remove_vertex ();
        }
      }
      else
      {
        // The remaining target vertex now becomes an isolated vertex inside
        // the containing face:
        _insert_isolated_vertex (f1, he1->vertex());
      }

      if (remove_source)
      {
        if (he2->vertex()->parameter_space_in_x() != ARR_INTERIOR || 
            he2->vertex()->parameter_space_in_y() != ARR_INTERIOR)
        {
          he2->vertex()->set_halfedge(NULL);    // disconnect the end vertex
          _remove_vertex_if_redundant (he2->vertex(), f1);
        }
        else
        {
          // Delete the he1's source vertex and its associated point.
          _notify_before_remove_vertex (Vertex_handle (he2->vertex()));
          
          _delete_point (he2->vertex()->point());
          _dcel().delete_vertex (he2->vertex());
          
          _notify_after_remove_vertex ();
        }
      }
      else
      {
        // The remaining source vertex now becomes an isolated vertex inside
        // the containing face:
        _insert_isolated_vertex (f1, he2->vertex());
      }

      // Delete the curve associated with the edge to be removed.
      _delete_curve (he1->curve());
      _dcel().delete_edge (he1);

      // Notify the observers that an edge has been deleted.
      _notify_after_remove_edge();

      // Return the face that used to contain the hole.
      return (f1);
    }
    else if (he1->next() == he2 || he2->next() == he1)
    {
      CGAL_assertion (oc1 == oc2 && ic1 == ic2);

      // In this case the two halfedges form an "antenna" (case 2).
      // Make he1 point at the tip of this "antenna" (swap the pointer if
      // necessary).
      bool     remove_tip_vertex = remove_target;
 
      if (he2->next() == he1)
      {
        he1 = he2;
        he2 = he1->opposite();
        remove_tip_vertex = remove_source;
      }

      // Remove the two halfedges from the boundary chain by connecting
      // he1's predecessor with he2's successor.
      prev1 = he1->prev();
      prev1->set_next (he2->next());

      // In case the halfedges to be deleted are represantatives of their
      // CCB (note that noth should belong to the same CCB, be it an outer
      // CCB or an inner one), make prev1 the components representative.
      if (oc1 != NULL && (oc1->halfedge() == he1 ||
                          oc1->halfedge() == he2))
      {
        oc1->set_halfedge (prev1);
      }
      else if (ic1 != NULL && (ic1->halfedge() == he1 ||
                               ic1->halfedge() == he2))
      {
        ic1->set_halfedge (prev1);
      }

      // In case he2 is the representative halfedge of its target vertex,
      // replace it by prev1 (which also points at this vertex).
      if (he2->vertex()->halfedge() == he2)
        he2->vertex()->set_halfedge (prev1);

      // Try to temove the base vertex, in case it has boundary conditions.
      if (he2->vertex()->parameter_space_in_x() != ARR_INTERIOR || 
          he2->vertex()->parameter_space_in_y() != ARR_INTERIOR)
      {
        _remove_vertex_if_redundant (he2->vertex(), f1);
      }

      // Remove the redundant tip vertex, if necessary.
      if (remove_tip_vertex)
      {
        if (he1->vertex()->parameter_space_in_x() != ARR_INTERIOR || 
            he1->vertex()->parameter_space_in_y() != ARR_INTERIOR)
        {
          he1->vertex()->set_halfedge(NULL);    // disconnect the end vertex
          _remove_vertex_if_redundant (he1->vertex(), f1);
        }
        else
        {
          // Delete the vertex that forms the tip of the "antenna".
          _notify_before_remove_vertex (Vertex_handle (he1->vertex()));
          
          _delete_point (he1->vertex()->point());
          _dcel().delete_vertex (he1->vertex());

          _notify_after_remove_vertex();
        }
      }
      else
      {
        // The remaining "antenna" tip now becomes an isolated vertex inside
        // the containing face:
        _insert_isolated_vertex (f1, he1->vertex());
      }

      // Delete the curve associated with the edge to be removed.
      _delete_curve (he1->curve());
      _dcel().delete_edge (he1);

      // Notify the observers that an edge has been deleted.
      _notify_after_remove_edge();

      // Return the incident face.
      return (f1);
    }

    // In this case the degree of both end-vertices is at least 2, so we
    // can use the two predecessor halfedges of he1 and he2.
    bool        add_inner_ccb = false;

    prev1 = he1->prev();
    prev2 = he2->prev();

    if (ic1 != NULL && ic1 == ic2)
    {
      // If both halfedges lie on the same inner component (hole) inside the
      // face (case 3.1), we have to split this component into two holes.
      //
      //    +-----------------------------+
      //    |           prev1             |
      //    |   +----+ /    +----+        |
      //    |   |    +......+    |        |
      //    |   +----+      +----+        |
      //    |                             |
      //    +-----------------------------+
      //
      // Notify the observers we are about to split an inner CCB.
      _notify_before_split_inner_ccb (Face_handle (f1),
                                      (Halfedge_handle 
                                       (*(ic1->iterator())))->ccb(),
                                      Halfedge_handle (he1));

      // We first make prev1 the new representative halfedge of the first
      // inner CCB.
      ic1->set_halfedge (prev1);

      // Create a new component that represents the new hole we split.
      DInner_ccb     *new_ic = _dcel().new_inner_ccb();

      f1->add_inner_ccb (new_ic, prev2);
      new_ic->set_face (f1);

      // Associate all halfedges along the hole boundary with the new inner
      // component.
      DHalfedge       *curr;

      for (curr = he1->next(); curr != he2; curr = curr->next())
        curr->set_inner_ccb (new_ic);

      // Notify the observers that the hole has been split.
      _notify_after_split_inner_ccb (Face_handle (f1),
                                     (Halfedge_handle (prev1))->ccb(),
                                     (Halfedge_handle (prev2))->ccb());
    }
    else if (oc1 != oc2)
    {
      // RWRW: NEW!
      CGAL_assertion (oc1 != NULL && oc2 != NULL);

      // In case both halfegdes he1 and he2 are incident to the same face
      // but lie on different outer CCBs of this face, removing this pair of
      // halfedge causes the two components two merge and to become an
      // inner CCB in the face.
      // We first remove the outer CCB oc1 from f, and inform the observers
      // on doing so.
      Face_handle               fh (f1);

      _notify_before_remove_outer_ccb (fh, (Halfedge_handle(he1))->ccb());

      f1->erase_outer_ccb (oc1);
      _dcel().delete_outer_ccb (oc1);

      _notify_after_remove_outer_ccb (fh);

      // We now remove the outer CCBs oc2 from f, and inform the observers
      // on doing so.
      _notify_before_remove_outer_ccb (fh, (Halfedge_handle(he2))->ccb());

      f2->erase_outer_ccb (oc2);
      _dcel().delete_outer_ccb (oc2);

      _notify_after_remove_outer_ccb (fh);

      // Mark that we should eventually add a new inner CCB inside the face.
      add_inner_ccb = true;
    }
    else
    {
      CGAL_assertion (oc1 != NULL && oc1 == oc2);

      // If both halfedges are incident to the same outer CCB of their
      // face (case 3.2), we have to distinguish two sub-cases:
      if (m_topol_traits.hole_creation_after_edge_removal (he1))
      {
        // We have to create a new hole in the interior of the incident face
        // (case 3.2.1):
        //
        //    +-----------------------------+
        //    | prev1                       |
        //    v         +----+              |
        //    +........>+    |              |
        //    |   he1   +----+              |
        //    |                             |
        //    +-----------------------------+
        //
        // Note that it is guaranteed that he1 points at this new hole, while
        // he2 points at the boundary of the face that contains this hole.
        // First notify the observers we are about to form a new inner
        // CCB inside f1.
        _notify_before_add_inner_ccb (Face_handle (f1),
                                      Halfedge_handle (he1->next()));

        // Create a new component that represents the new hole.
        DInner_ccb   *new_ic = _dcel().new_inner_ccb();

        f1->add_inner_ccb (new_ic, he1->next());
        new_ic->set_face (f1);

        // Associate all halfedges along the hole boundary with the new inner
        // component.
        DHalfedge       *curr;

        for (curr = he1->next(); curr != he2; curr = curr->next())
          curr->set_inner_ccb (new_ic);

        // As the outer CCB of f1 may be represented by any of the
        // halfedges in between he1 -> ... -> he2 (the halfedges in between
        // represent the outer boundary of the new hole that is formed),
        // We represent the outer boundary of f1 by prev1, which definately
        // stays on the outer boundary.
        oc1->set_halfedge (prev1);

        // Notify the observers that a new hole has been formed.
        Ccb_halfedge_circulator   hccb = (Halfedge_handle(he1->next()))->ccb();

        _notify_after_add_inner_ccb (hccb);
      }
      else
      {
        // We have to split the outer CCB into two outer components
        // (case 3.2.2), such that the number of outer CCBs of the face is
        // incremented.
        // First we notify the observers that we are about to split an outer
        // component.
        _notify_before_split_outer_ccb (Face_handle (f1),
                                        Halfedge_handle (he1)->ccb(),
                                        Halfedge_handle (he1));

        // Create a new outer component.
        DOuter_ccb   *new_oc = _dcel().new_outer_ccb();

        f1->add_outer_ccb (new_oc, he1->next());
        new_oc->set_face (f1);

        // Associate all halfedges from he1 until he2 with the new CCB.
        DHalfedge       *curr;

        for (curr = he1->next(); curr != he2; curr = curr->next())
          curr->set_outer_ccb (new_oc);

        // As the outer CCB of f1 may be represented by any of the
        // halfedges in between he1 -> ... -> he2 (the halfedges in between
        // are on the new outer CCB we have just created), we represent the
        // former outer CCB by prev1, which definately stays on it.
        oc1->set_halfedge (prev1);

        // Notify the observers that a new outer CCB has been formed.
        _notify_after_split_outer_ccb (Face_handle (f1),
                                       Halfedge_handle (he1->next())->ccb(),
                                       Halfedge_handle (prev1)->ccb());
      }
    }

    // Disconnect the two halfedges we are about to delete from the edge list.
    prev1->set_next (he2->next());
    prev2->set_next (he1->next());

    // If one of these edges is an incident halfedge of its target vertex,
    // replace it by the appropriate predecessor.
    if (he1->vertex()->halfedge() == he1)
      he1->vertex()->set_halfedge (prev2);

    if (he2->vertex()->halfedge() == he2)
      he2->vertex()->set_halfedge (prev1);

    // Remove the end vertices, in case they become redundant.
    if (he1->vertex()->parameter_space_in_x() != ARR_INTERIOR || 
        he1->vertex()->parameter_space_in_y() != ARR_INTERIOR)
    {
      _remove_vertex_if_redundant (he1->vertex(), f1);
    }

    if (he2->vertex()->parameter_space_in_x() != ARR_INTERIOR || 
        he2->vertex()->parameter_space_in_y() != ARR_INTERIOR)
    {
      _remove_vertex_if_redundant (he2->vertex(), f1);
    }

    // Delete the curve associated with the edge to be removed.
    _delete_curve (he1->curve());

    // Delete the pair of halfedges.
    _dcel().delete_edge (he1);

    // RWRW: NEW!
    // In case we have to create a new inner CCB inside the face (new removal
    // case), do it now.
    if (add_inner_ccb)
    {
      // Notify the observers that we are about to create a new inner CCB
      // inside the merged face.
      Halfedge_handle   hh (prev1);

      _notify_before_add_inner_ccb (Face_handle (f1), hh);

      // Initiate a new inner CCB inside the given face.
      DInner_ccb   *new_ic = _dcel().new_inner_ccb();

      f1->add_inner_ccb (new_ic, prev1);
      new_ic->set_face (f1);

      // Set the innser CCB of the halfedges along the component boundary.
      DHalfedge   *curr = prev1;

      do
      {
        curr->set_inner_ccb (new_ic);
        curr = curr->next();
      } while (curr != prev1);

      // Notify the observers that we have formed a new inner CCB.
      _notify_after_add_inner_ccb (hh->ccb());
    }

    // Notify the observers that an edge has been deleted.
    _notify_after_remove_edge();

    // Return the incident face.
    return (f1);
  }

  // The two incident faces are not the same - in this case, the edge we are
  // about to delete separates these two faces. We therefore have to delete
  // one of these faces and merge it with the other face.
  // First notify the observers we are about to merge the two faces.
  _notify_before_merge_face (Face_handle (f1), Face_handle (f2),
                             Halfedge_handle (he1));

  // We begin by checking whether one of the face is a hole inside the other
  // face.
  DHalfedge   *curr;

  prev1 = he1->prev();
  prev2 = he2->prev();

  CGAL_assertion (ic1 == NULL || ic2 == NULL);

  if (ic1 == NULL && ic2 == NULL)
  {
    bool          add_inner_ccb = false;

    // Both halfedges lie on the outer boundary of their incident faces
    // (case 3.4). We have to distinguish two possible sub-cases.
    if (m_topol_traits.hole_creation_after_edge_removal (he1))
    {
      // We have to remove the outer CCBs of f1 and f2 that he1 and he2 lie
      // on, and create a new hole in the merged face (case 3.4.1).
      // We first remove the outer CCB oc1 from f1, and inform the observers
      // on doing so.
      _notify_before_remove_outer_ccb (Face_handle (f1),
                                       (Halfedge_handle(he1))->ccb());

      f1->erase_outer_ccb (oc1);
      _dcel().delete_outer_ccb (oc1);

      _notify_after_remove_outer_ccb (Face_handle (f1));

      // We now remove the outer CCBs oc2 from f2, and inform the observers
      // on doing so.
      _notify_before_remove_outer_ccb (Face_handle (f2),
                                       (Halfedge_handle(he2))->ccb());

      f2->erase_outer_ccb (oc2);
      _dcel().delete_outer_ccb (oc2);

      _notify_after_remove_outer_ccb (Face_handle (f2));

      // Mark that we should eventually add a new inner CCB in the merged face.
      add_inner_ccb = true;
    }
    else
    {
      // f1 and f2 are two adjacent faces (case 3.4.2), so we simply merge
      // them.
      // We first set the connected component of f2's outer-boundary halfedges
      // to be the same as f1's outer component.
      for (curr = he2->next(); curr != he2; curr = curr->next())
        curr->set_outer_ccb (oc1);
    }

    // Move the holes inside f2 to f1.
    DInner_ccb_iter    ic_it = f2->inner_ccbs_begin();
    DInner_ccb_iter    ic_to_move;
    
    while (ic_it != f2->inner_ccbs_end())
    {
      // We increment the holes itrator before moving the hole, because
      // this operation invalidates the iterator.
      ic_to_move  = ic_it;
      ++ic_it;
      
      _move_inner_ccb (f2, f1, *ic_to_move);
    }

    // In case he1, which is about to be deleted, is a representative
    // halfedge of outer component of f1, we replace it by its predecessor.
    if (oc1->halfedge() == he1)
      oc1->set_halfedge (prev1);

    // Move the isolated vertices inside f2 to f1.
    DIso_vertex_iter    iv_it = f2->isolated_vertices_begin();
    DIso_vertex_iter    iv_to_move;
    
    while (iv_it != f2->isolated_vertices_end())
    {
      // We increment the isolated vertices itrator before moving the vertex,
      // because this operation invalidates the iterator.
      iv_to_move  = iv_it;
      ++iv_it;

      _move_isolated_vertex (f2, f1, &(*iv_to_move)); 
    }

    // If he1 or he2 are the incident halfedges to their target vertices,
    // we replace them by the appropriate predecessors.
    if (he1->vertex()->halfedge() == he1)
      he1->vertex()->set_halfedge (prev2);
      
    if (he2->vertex()->halfedge() == he2)
      he2->vertex()->set_halfedge (prev1);
    
    // Disconnect the two halfedges we are about to delete from the edge
    // list.
    prev1->set_next (he2->next());
    prev2->set_next (he1->next());
      
      // Delete the curve associated with the edge to be removed.
    _delete_curve (he1->curve());

    // If the face f2 we have just merged with f1 is unbounded, then the
    // merged face is also unbounded.
    if (f2->is_unbounded())
      f1->set_unbounded(true);

    // Delete the face f2.
    _dcel().delete_face (f2);
      
    // Notify the observers that the faces have been merged.
    _notify_after_merge_face (Face_handle (f1));

    // Remove the end vertices, in case they become redundant.
    if (he1->vertex()->parameter_space_in_x() != ARR_INTERIOR || 
        he1->vertex()->parameter_space_in_y() != ARR_INTERIOR)
    {
      _remove_vertex_if_redundant (he1->vertex(), f1);
    }
    
    if (he2->vertex()->parameter_space_in_x() != ARR_INTERIOR || 
        he2->vertex()->parameter_space_in_y() != ARR_INTERIOR)
    {
      _remove_vertex_if_redundant (he2->vertex(), f1);
    }
    
    // Delete the pair of halfedges.
    _dcel().delete_edge (he1);

    // In case we have to create a new inner CCB inside the merged face
    // (case 3.4.1), do it now.
    if (add_inner_ccb)
    {
      // Notify the observers that we are about to create a new inner CCB
      // inside the merged face.
      Halfedge_handle   hh (prev1);

      _notify_before_add_inner_ccb (Face_handle (f1), hh);

      // Initiate a new inner CCB inside the given face.
      DInner_ccb   *new_ic = _dcel().new_inner_ccb();

      f1->add_inner_ccb (new_ic, prev1);
      new_ic->set_face (f1);

      // Set the innser CCB of the halfedges along the component boundary.
      curr = prev1;
      do
      {
        curr->set_inner_ccb (new_ic);
        curr = curr->next();
      } while (curr != prev1);

      // Notify the observers that we have formed a new inner CCB.
      _notify_after_add_inner_ccb (hh->ccb());
    }

    // Notify the observers that an edge has been deleted.
    _notify_after_remove_edge();
      
    // Return the merged face.
    return (f1); 
  }

  // In this case we merge a face with another face that now forms a hole
  // inside it (case 3.3). We first make sure that f1 contains the hole f2, so
  // we can merge f2 with it (we swap roles between the halfedges if
  // necessary).
  if (ic2 != NULL)
  { 
    he1 = he2;
    he2 = he1->opposite();

    ic1 = ic2;
    ic2 = NULL;

    oc2 = oc1;
    oc1 = NULL;

    DFace   *tf = f1;
    f1 = f2;
    f2 = tf;
    
    prev1 = he1->prev();
    prev2 = he2->prev();
  }

  // By removing the edge we open a closed face f2 contained in f1. By doing
  // this, the outer boundary of f2 unites with the hole boundary that ic1
  // represents. We therefore have to set the component of all halfedges
  // along the boundary of f2 to be ic1.
  for (curr = he2->next(); curr != he2; curr = curr->next())
    curr->set_inner_ccb (ic1);

  // Move the inner CCBs inside f2 to f1.
  DInner_ccb_iter    ic_it = f2->inner_ccbs_begin();
  DInner_ccb_iter    ic_to_move;

  while (ic_it != f2->inner_ccbs_end())
  {
    // We increment the holes itrator before moving the hole, because
    // this operation invalidates the iterator.
    ic_to_move  = ic_it;
    ++ic_it;

    _move_inner_ccb (f2, f1, *ic_to_move); 
  }

  // Move the isolated vertices inside f2 to f1.
  DIso_vertex_iter    iv_it = f2->isolated_vertices_begin();
  DIso_vertex_iter    iv_to_move;

  while (iv_it != f2->isolated_vertices_end())
  {
    // We increment the isolated vertices itrator before moving the vertex,
    // because this operation invalidates the iterator.
    iv_to_move  = iv_it;
    ++iv_it;

    _move_isolated_vertex (f2, f1, &(*iv_to_move)); 
  }
        
  // Notice that f2 will be merged with f1, but its boundary will still be
  // a hole inside this face. In case he1 is a represantative of this hole,
  // replace it by its predecessor.
  if (ic1->halfedge() == he1)
    ic1->set_halfedge (prev1);
  
  // If he1 or he2 are the incident halfedges to their target vertices,
  // we replace them by the appropriate predecessors.
  if (he1->vertex()->halfedge() == he1)
    he1->vertex()->set_halfedge (prev2);
        
  if (he2->vertex()->halfedge() == he2)
    he2->vertex()->set_halfedge (prev1);
      
  // Disconnect the two halfedges we are about to delete from the edge
  // list.
  prev1->set_next (he2->next());
  prev2->set_next (he1->next());

  // Delete the curve associated with the edge to be removed.
  _delete_curve (he1->curve());

  // If the face f2 we have just merged with f1 is unbounded, then the merged
  // face is also unbounded.
  if (f2->is_unbounded())
    f1->set_unbounded(true);

  // Delete the face f2.
  _dcel().delete_face (f2);
      
  // Notify the observers that the faces have been merged.
  _notify_after_merge_face (Face_handle (f1));

  // Remove the end vertices, in case they become redundant.
  if (he1->vertex()->parameter_space_in_x() != ARR_INTERIOR || 
      he1->vertex()->parameter_space_in_y() != ARR_INTERIOR)
  {
    _remove_vertex_if_redundant (he1->vertex(), f1);
  }
  
  if (he2->vertex()->parameter_space_in_x() != ARR_INTERIOR || 
      he2->vertex()->parameter_space_in_y() != ARR_INTERIOR)
  {
    _remove_vertex_if_redundant (he2->vertex(), f1);
  }

  // Delete the pair of halfedges.
  _dcel().delete_edge (he1);

  // Notify the observers that an edge has been deleted.
  _notify_after_remove_edge();

  // Return the merged face.
  return (f1);
}

//-----------------------------------------------------------------------------
// Remove a vertex in case it becomes redundant after the deletion of an
// incident edge.
//
template<class GeomTraits, class TopTraits>
void
Arrangement_on_surface_2<GeomTraits, TopTraits>::
_remove_vertex_if_redundant (DVertex *v, DFace *f)
{
  CGAL_precondition (v->parameter_space_in_x() != ARR_INTERIOR || 
                     v->parameter_space_in_y() != ARR_INTERIOR);

  // In case the vertex has no incident halfedges, remove it if it is
  // redundant. Otherwise, make it an isolated vertex.
  if (v->halfedge() == NULL)
  {
    if (m_topol_traits.is_redundant (v))
    {
      // Remove the vertex and notify the observers on doing so.
      _notify_before_remove_vertex (Vertex_handle (v));

      m_topol_traits.erase_redundant_vertex (v);

      // Note the topology traits do not free the vertex - we now do it.
      if (! v->has_null_point())
        _delete_point (v->point());
      _dcel().delete_vertex (v);

      _notify_after_remove_vertex();
    }
    else
    {
      // Keep the vertex as an isolated one.
      _insert_isolated_vertex (f, v);
    }

    return;
  }

  // Get the first two incident halfedges of v.
  DHalfedge  *he1 = v->halfedge();
  DHalfedge  *he2 = he1->next()->opposite();
  
  if (he2->next()->opposite() != he1)
  {
    // In this case there are more than two incident edges, so v obviously
    // cannot be removed.
    return;
  }

  if (! he1->has_null_curve() || ! he2->has_null_curve())
  {
    // We can only merge fictitious halfedges.
    return;
  }

  // Now check if the vertex is redundant. If it is, remove it by merging
  // its two incident fictitious halfedges.
  if (m_topol_traits.is_redundant (v))
  {
    // Use the topology traits to merge the two fictitious halfedges.
    _notify_before_merge_fictitious_edge (Halfedge_handle (he1),
                                          Halfedge_handle (he2));

    he1 = m_topol_traits.erase_redundant_vertex (v);

    _notify_after_merge_fictitious_edge (Halfedge_handle (he1));

    // Note the topology traits do not free the vertex - we now do it.
    _notify_before_remove_vertex (Vertex_handle (v));
    
    if (! v->has_null_point())
      _delete_point (v->point());
    _dcel().delete_vertex (v);

    _notify_after_remove_vertex();
  }

  return;
}

//-----------------------------------------------------------------------------
// Remove an isolated vertex from the interior of a given face (but not from
// the DCEL).
//
template<class GeomTraits, class TopTraits>
void Arrangement_on_surface_2<GeomTraits, TopTraits>::
_remove_isolated_vertex (DVertex* v)
{
  // Remove the isolated vertex from the face and delete its record.
  DIso_vertex  *iv = v->isolated_vertex();
  DFace        *f = iv->face();

  f->erase_isolated_vertex (iv);
  _dcel().delete_isolated_vertex (iv);
  
  return;
}

//---------------------------------------------------------------------------
// Check whether the arrangement is valid. In particular, check the
// validity of each vertex, halfedge, and face.
//
template<class GeomTraits, class TopTraits>
bool Arrangement_on_surface_2<GeomTraits, TopTraits>::is_valid() const
{
  Vertex_const_iterator    vit;
  bool                     is_vertex_valid;

  for (vit = vertices_begin(); vit != vertices_end(); ++vit)
  {
    is_vertex_valid = _is_valid (vit);
    if (!is_vertex_valid) 
    {
      CGAL_warning_msg (is_vertex_valid, "Invalid vertex.");
      return (false);
    }
  }
    
  Halfedge_const_iterator  heit;
  bool                     is_halfedge_valid;

  for (heit = halfedges_begin(); heit != halfedges_end(); ++heit) 
  {
    is_halfedge_valid = _is_valid (heit);
    if (! is_halfedge_valid) 
    {
      CGAL_warning_msg (is_halfedge_valid, "Invalid halfedge.");
      return (false);
    }
  }
    
  Face_const_iterator     fit;
  bool                    is_face_valid;

  for (fit = faces_begin(); fit != faces_end(); ++fit) 
  {
    is_face_valid = _is_valid (fit);
    if (! is_face_valid)
    {
      CGAL_warning_msg (is_face_valid, "Invalid face.");
      return (false);
    }
  }

  bool  are_vertices_unique = _are_vertices_unique();
  if (! are_vertices_unique)
  {
    CGAL_warning_msg (are_vertices_unique,
                      "Found two vertices with the same geometric point.");
    return (false);
  }

  // If we reached here, the arrangement is valid.
  return (true);
}

//---------------------------------------------------------------------------
// Check the validity of a vertex.
//
template<class GeomTraits, class TopTraits>
bool Arrangement_on_surface_2<GeomTraits, TopTraits>::
_is_valid (Vertex_const_handle v) const
{    
  // Do not check isolated vertices, as they have no incident halfedges.
  if (v->is_isolated())
    return (true);

  // Make sure that the vertex is the target of all its incident halfedges.
  Halfedge_around_vertex_const_circulator circ = v->incident_halfedges();
  Halfedge_around_vertex_const_circulator start = circ;

  do
  {
    if (circ->target() != v)
      return (false);
    
    ++circ;
  } while (circ != start);

  // In case of a non-boundary vertex, make sure the curves are correctly
  // ordered around this vertex.
  if (v->parameter_space_in_x() ==
      ARR_INTERIOR && v->parameter_space_in_y() == ARR_INTERIOR)
  {
    if (! _are_curves_ordered_cw_around_vertrex (v))
      return (false);
  }

  return (true);
}

//---------------------------------------------------------------------------
// Check the validity of a halfedge.
//
template<class GeomTraits, class TopTraits>
bool Arrangement_on_surface_2<GeomTraits, TopTraits>::
_is_valid (Halfedge_const_handle he) const
{
  // Check relations with the previous and the next halfedges.
  if (he->prev()->target() != he->source())
    return (false);

  if (he->target() != he->next()->source())
    return (false);

  // Check relations with the twin.
  if (he != he->twin()->twin())
    return (false);
  
  if (he->source() != he->twin()->target() || 
      he->target() != he->twin()->source())
    return (false);

  if (he->direction() == he->twin()->direction())
    return (false);

  // Stop here in case of a fictitious edge.
  if (he->is_fictitious())
    return (true);

  // Check that the end points of the curve associated with the halfedge
  // really equal the source and target vertices of this halfedge.
  const X_monotone_curve_2&  cv = he->curve();
  const DVertex             *source = _vertex (he->source());
  const DVertex             *target = _vertex (he->target());
  Comparison_result          res;

  if (source->parameter_space_in_x() == ARR_INTERIOR &&
      source->parameter_space_in_y() == ARR_INTERIOR &&
      target->parameter_space_in_x() == ARR_INTERIOR &&
      target->parameter_space_in_y() == ARR_INTERIOR)
  {
    res = m_geom_traits->compare_xy_2_object()(source->point(), target->point());
  }
  else
  {
    res = (he->direction() == ARR_LEFT_TO_RIGHT) ? SMALLER : LARGER;
  }

  if (res == SMALLER)
  {
    if (he->direction() != ARR_LEFT_TO_RIGHT)
      return (false);

    return (_are_equal (_vertex (he->source()), cv, ARR_MIN_END) &&
            _are_equal (_vertex (he->target()), cv, ARR_MAX_END));
  }
  else if (res == LARGER)
  {
    if (he->direction() != ARR_RIGHT_TO_LEFT)
      return (false);

    return (_are_equal (_vertex (he->source()), cv, ARR_MAX_END) &&
            _are_equal (_vertex (he->target()), cv, ARR_MIN_END));
  }

  // In that case, the source and target of the halfedge are equal.
  return (false);
}

//---------------------------------------------------------------------------
// Check the validity of a face.
//
template<class GeomTraits, class TopTraits>
bool Arrangement_on_surface_2<GeomTraits, TopTraits>::
_is_valid (Face_const_handle f) const
{
  // Check if all outer components of the face refer to f.
  const DFace                   *p_f = _face (f);
  DOuter_ccb_const_iter          oc_it;
  const DHalfedge               *he;
  const DOuter_ccb              *oc;

  for (oc_it = p_f->outer_ccbs_begin();
       oc_it != p_f->outer_ccbs_end(); ++oc_it)
  {
    he = *oc_it;
    if (he->is_on_inner_ccb())
      return (false);

    oc = he->outer_ccb();
    if (oc->face() != p_f)
      return (false);

    if (! _is_outer_ccb_valid (oc, he))
      return (false);
  }
 
  // Check if all inner components of the face refer to f.
  DInner_ccb_const_iter          ic_it;
  const DInner_ccb              *ic;

  for (ic_it = p_f->inner_ccbs_begin();
       ic_it != p_f->inner_ccbs_end(); ++ic_it)
  {
    he = *ic_it;
    if (! he->is_on_inner_ccb())
      return (false);

    ic = he->inner_ccb();
    if (ic->face() != p_f)
      return (false);

    if (! _is_inner_ccb_valid (ic, he))
      return (false);
  }
 
  // Check if all isolated vertices inside the face refer to f.
  DIso_vertex_const_iter          iv_it;
  const DVertex                  *v;
  const DIso_vertex              *iv;

  for (iv_it = p_f->isolated_vertices_begin();
       iv_it != p_f->isolated_vertices_end(); ++iv_it)
  {
    v = &(*iv_it);
    if (! v->is_isolated())
      return (false);

    iv = v->isolated_vertex();
    if (iv->face() != p_f)
      return (false);
  }

  // If we reached here, the face is valid.
  return (true);
}

//---------------------------------------------------------------------------
// Check the validity of an outer CCB.
//
template<class GeomTraits, class TopTraits>
bool Arrangement_on_surface_2<GeomTraits, TopTraits>::
_is_outer_ccb_valid (const DOuter_ccb *oc, const DHalfedge *first) const
{
  // Make sure that all halfedges along the CCB refer to the same component.
  const DHalfedge   *curr = first;
  bool               found_rep = false;

  do
  {
    if (curr->is_on_inner_ccb())
      return (false);

    if (oc != curr->outer_ccb())
      return (false);

    if (! found_rep && oc->halfedge() == curr)
      found_rep = true;

    curr = curr->next();
  } while (curr != first);
  
  // Return if we found the CCB representative along the outer CCB.
  return (found_rep);
}

//---------------------------------------------------------------------------
// Check the validity of an inner CCB.
//
template<class GeomTraits, class TopTraits>
bool Arrangement_on_surface_2<GeomTraits, TopTraits>::
_is_inner_ccb_valid (const DInner_ccb *ic, const DHalfedge *first) const
{
  // Make sure that all halfedges along the CCB refer to the same component.
  const DHalfedge   *curr = first;
  bool               found_rep = false;

  do
  {
    if (! curr->is_on_inner_ccb())
      return (false);

    if (ic != curr->inner_ccb())
      return (false);

    if (! found_rep && ic->halfedge() == curr)
      found_rep = true;

    curr = curr->next();
  } while (curr != first);
  
  // Return if we found the CCB representative along the inner CCB.
  return (found_rep);
}

//---------------------------------------------------------------------------
// Check that all vertices are unique (no two vertices with the same 
// geometric point).
//
template<class GeomTraits, class TopTraits>
bool Arrangement_on_surface_2<GeomTraits, TopTraits>::
_are_vertices_unique() const
{
  if (number_of_vertices() < 2)
    return (true);

  // Store all points associated with non-boundary vertices.
  std::vector<Point_2>  points_vec (number_of_vertices());
  Vertex_const_iterator vit;
  unsigned int          i = 0;
  
  for (vit = vertices_begin(); vit != vertices_end(); ++vit)
  {
    if (vit->parameter_space_in_x() == ARR_INTERIOR &&
        vit->parameter_space_in_y() == ARR_INTERIOR)
    { 
      points_vec[i] = vit->point();
      ++i;
    }
  }
  points_vec.resize (i);

  // Sort the vector of points and make sure no two adjacent points in the
  // sorted vector are equal.
  typedef typename Traits_adaptor_2::Compare_xy_2       Compare_xy_2;
  typedef typename Traits_adaptor_2::Equal_2            Equal_2;
  
  Equal_2       equal = m_geom_traits->equal_2_object();
  Compare_xy_2  compare_xy = m_geom_traits->compare_xy_2_object();
  Compare_to_less<Compare_xy_2>       cmp = compare_to_less (compare_xy);
  
  std::sort(points_vec.begin(), points_vec.end(), cmp);
  for (i = 1; i < points_vec.size(); ++i)
  {
    if (equal (points_vec[i-1], points_vec[i]))
      return (false);
  }
  
  return (true);
}

//---------------------------------------------------------------------------
// Check that the curves around a given vertex are ordered clockwise.
//
template<class GeomTraits, class TopTraits>
bool Arrangement_on_surface_2<GeomTraits, TopTraits>::
_are_curves_ordered_cw_around_vertrex (Vertex_const_handle v) const
{
  if(v->degree() < 3)
    return (true);

  typename Traits_adaptor_2::Is_between_cw_2  is_between_cw =
    m_geom_traits ->is_between_cw_2_object();

  Halfedge_around_vertex_const_circulator circ = v->incident_halfedges();
  Halfedge_around_vertex_const_circulator first = circ;
  Halfedge_around_vertex_const_circulator prev, next;
  bool                                    eq1, eq2;

  do
  {
    prev = circ; --prev;
    next = circ; ++next;

    if (!is_between_cw(circ->curve(), (circ->direction() == ARR_RIGHT_TO_LEFT),
                       prev->curve(), (prev->direction() == ARR_RIGHT_TO_LEFT),
                       next->curve(), (next->direction() == ARR_RIGHT_TO_LEFT),
                       v->point(), eq1, eq2))
      return (false);

    if (eq1 || eq2)
      return (false);
    
    ++circ;
  } while (circ != first);
   
  return (true);
}

} //namespace CGAL

#endif
