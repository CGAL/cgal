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
// $URL$
// $Id$
// 
//
// Author(s)     : Ron Wein          <wein@post.tau.ac.il>
//                 Efi Fogel         <efif@post.tau.ac.il>
//                 (based on old version by: Iddo Hanniel,
//                                           Eyal Flato,
//                                           Oren Nechushtan,
//                                           Ester Ezra,
//                                           Shai Hirsch,
//                                           and Eugene Lipovetsky)
//
#ifndef CGAL_ARRANGEMENT_2_FUNCTIONS_H
#define CGAL_ARRANGEMENT_2_FUNCTIONS_H

/*! \file
 * Member-function definitions for the Arrangement_2<Traits,Dcel> class.
 */

#include <CGAL/function_objects.h> 

CGAL_BEGIN_NAMESPACE

//-----------------------------------------------------------------------------
// Default constructor.
//
template<class Traits, class Dcel>
Arrangement_2<Traits,Dcel>::Arrangement_2 ()
{
  // Construct an empty arrangement.
  _construct_empty_arrangement ();

  // Allocate the traits.
  traits = new Traits_adaptor_2;
  own_traits = true;
}

//-----------------------------------------------------------------------------
// Copy constructor.
//
template<class Traits, class Dcel>
Arrangement_2<Traits,Dcel>::Arrangement_2 (const Self& arr) :
  v_bl (NULL),
  v_tl (NULL),
  v_br (NULL),
  v_tr (NULL),
  n_inf_verts (0),
  un_face (NULL),
  traits (NULL),
  own_traits (false)
{
  assign (arr);
}

//-----------------------------------------------------------------------------
// Constructor given a traits object.
//
template<class Traits, class Dcel>
Arrangement_2<Traits,Dcel>::Arrangement_2 (Traits_2 *tr)
{
  // Construct an empty arrangement.
  _construct_empty_arrangement ();

  // Set the traits.
  traits = static_cast<Traits_adaptor_2*>(tr);
  own_traits = false;
}

//-----------------------------------------------------------------------------
// Assignment operator.
//
template<class Traits, class Dcel>
Arrangement_2<Traits,Dcel>&
Arrangement_2<Traits,Dcel>::operator= (const Self& arr)
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
template<class Traits, class Dcel>
void Arrangement_2<Traits,Dcel>::assign (const Self& arr)
{
  // Clear the current contents of the arrangement.
  clear();

  // Notify the observers that an assignment is to take place.
  _notify_before_assign (arr);

  // Duplicate the DCEL.
  un_face = dcel.assign (arr.dcel, arr.un_face);
  n_inf_verts = arr.n_inf_verts;

  // Go over the vertices and create duplicates of the stored points.
  typename Dcel::Vertex_iterator       vit;
  Point_2                             *dup_p;
  DVertex                             *p_v;
  Boundary_type                        inf_x, inf_y;
  const DHalfedge                     *first_he, *next_he;

  for (vit = dcel.vertices_begin(); vit != dcel.vertices_end(); ++vit)
  {
    p_v = &(*vit);
    inf_x = p_v->boundary_in_x();
    inf_y = p_v->boundary_in_y();

    if (inf_x == NO_BOUNDARY && inf_y == NO_BOUNDARY)
    {
      // Normal vertex (not at inifinity): Create the duplicate point and
      // store it in the points container.
      dup_p = _new_point (p_v->point());

      // Associate the vertex with the duplicated point.
      p_v->set_point (dup_p);
    }
    else
    {
      // Check if this is one of the ficitious vertices at infinity. These
      // vertices have degree 2, as opposed to other vertices at infinity
      // which have degree 3.
      first_he = p_v->halfedge();
      next_he = first_he->next()->opposite();
      
      if (next_he->next()->opposite() == first_he)
      {
        if (inf_x == MINUS_INFINITY && inf_y == MINUS_INFINITY)
          v_bl = p_v;
        else if (inf_x == MINUS_INFINITY && inf_y == PLUS_INFINITY)
          v_tl = p_v;
        else if (inf_x == PLUS_INFINITY && inf_y == MINUS_INFINITY)
          v_br = p_v;
        else if (inf_x == PLUS_INFINITY && inf_y == PLUS_INFINITY)
          v_tr = p_v;
        else
          // We should never reach here:
          CGAL_error();
      }
    }
  }

  // Go over the edge and create duplicates of the stored curves.
  typename Dcel::Edge_iterator         eit;
  X_monotone_curve_2                  *dup_cv;
  DHalfedge                           *p_e;

  for (eit = dcel.edges_begin(); eit != dcel.edges_end(); ++eit)
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
  if (own_traits && traits != NULL)
    delete traits;

  if (arr.own_traits)
    traits = new Traits_adaptor_2;
  else
    traits = arr.traits;
  own_traits = arr.own_traits;

  // Notify the observers that the assignment has been performed.
  _notify_after_assign ();

  return;
}

//-----------------------------------------------------------------------------
// Destructor.
//
template<class Traits, class Dcel>
Arrangement_2<Traits,Dcel>::~Arrangement_2 ()
{  
  // Free all stored points.
  typename Dcel::Vertex_iterator       vit;

  for (vit = dcel.vertices_begin(); vit != dcel.vertices_end(); ++vit)
  {
    if (! vit->has_null_point())
      _delete_point (vit->point());
  }

  // Free all stores curves.
  typename Dcel::Edge_iterator         eit;

  for (eit = dcel.edges_begin(); eit != dcel.edges_end(); ++eit)
  {
    if (! eit->has_null_curve())
      _delete_curve (eit->curve());
  }

  // Clear the DCEL.
  dcel.delete_all();

  // Free the traits object, if necessary.
  if (own_traits)
    delete traits;

  // Detach all observers still attached to the arrangement.
  Observers_iterator  iter = observers.begin();
  Observers_iterator  next;
  Observers_iterator  end = observers.end();

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
template<class Traits, class Dcel>
void Arrangement_2<Traits,Dcel>::clear ()
{
  // Notify the observers that we are about to clear the arragement.
  _notify_before_clear ();

  // Free all stored points.
  typename Dcel::Vertex_iterator       vit;

  for (vit = dcel.vertices_begin(); vit != dcel.vertices_end(); ++vit)
  {
    if (! vit->has_null_point())
      _delete_point (vit->point());
  }

  // Free all stores curves.
  typename Dcel::Edge_iterator         eit;

  for (eit = dcel.edges_begin(); eit != dcel.edges_end(); ++eit)
  {
    if (! eit->has_null_curve())
      _delete_curve (eit->curve());
  }

  // Clear the DCEL and construct an empty arrangement.
  dcel.delete_all();
  _construct_empty_arrangement ();

  // Notify the observers that we have just cleared the arragement.
  _notify_after_clear (unbounded_face());

  return;
}

//-----------------------------------------------------------------------------
// Insert a point as an isolated vertex in the interior of a given face.
//
template<class Traits, class Dcel>
typename Arrangement_2<Traits,Dcel>::Vertex_handle
Arrangement_2<Traits,Dcel>::insert_in_face_interior (const Point_2& p,
                                                     Face_handle f)
{
  DFace          *p_f = _face (f);

  // Create a new vertex associated with the given point.
  DVertex         *v =_create_vertex (p);
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
template<class Traits, class Dcel>
typename Arrangement_2<Traits,Dcel>::Halfedge_handle
Arrangement_2<Traits,Dcel>::insert_in_face_interior
    (const X_monotone_curve_2& cv, 
     Face_handle f)
{
  DFace            *p_f = _face (f);

  // Check if cv's left end lies at infinity, and create a vertex v1 that
  // corresponds to this end.
  const Boundary_type  inf_x1 = traits->boundary_in_x_2_object()(cv, MIN_END);
  const Boundary_type  inf_y1 = traits->boundary_in_y_2_object()(cv, MIN_END);
  DVertex             *v1 = NULL;
  DHalfedge           *fict_prev1 = NULL;

  if (inf_x1 == NO_BOUNDARY && inf_y1 == NO_BOUNDARY)
  {
    // The curve has a valid left endpoint: Create a new vertex associated
    // with the curve's left endpoint.
    v1 = _create_vertex (traits->construct_min_vertex_2_object()(cv));
  }
  else
  {
    // Locate a halfedge of f's outer CCB such that contains cv's left end
    // in its range.
    CGAL_precondition (f->is_unbounded());

    fict_prev1 = _locate_along_ccb (p_f, cv, MIN_END,
                                    inf_x1, inf_y1);

    CGAL_assertion (fict_prev1 != NULL);

    // Split the fictitious edge and create a vertex at infinity that
    // represents cv's left end.
    fict_prev1 = _split_fictitious_edge (fict_prev1,
                                         inf_x1, inf_y1);
  }

  // Check if cv's right end lies at infinity, and create a vertex v2 that
  // corresponds to this end.
  const Boundary_type  inf_x2 = traits->boundary_in_x_2_object()(cv, MAX_END);
  const Boundary_type  inf_y2 = traits->boundary_in_y_2_object()(cv, MAX_END);
  DVertex             *v2 = NULL;
  DHalfedge           *fict_prev2 = NULL;

  if (inf_x2 == NO_BOUNDARY && inf_y2 == NO_BOUNDARY)
  {
    // The curve has a valid right endpoint: Create a new vertex associated
    // with the curve's right endpoint.
    v2 = _create_vertex (traits->construct_max_vertex_2_object()(cv));
  }
  else
  {    
    // Locate a halfedge of f's outer CCB such that contains cv's right end
    // in its range.
    CGAL_precondition (f->is_unbounded());
    
    fict_prev2 = _locate_along_ccb (p_f, cv, MAX_END,
                                    inf_x2, inf_y2);

    CGAL_assertion (fict_prev2 != NULL);

    // Split the fictitious edge and create a vertex at infinity that
    // represents cv's right end.
    fict_prev2 = _split_fictitious_edge (fict_prev2,
                                         inf_x2, inf_y2);
  }

  // Create the edge connecting the two vertices (note we know v1 is 
  // lexicographically smaller than v2).
  DHalfedge       *new_he;

  if (fict_prev1 == NULL && fict_prev2 == NULL)
  {
    // Both vertices represent valid points.
    new_he = _insert_in_face_interior (cv, p_f, v1, v2,
                                       SMALLER);
  }
  else if (fict_prev1 == NULL && fict_prev2 != NULL)
  {
    // v1 represents a valid point and v2 is at infinity.
    new_he = _insert_from_vertex (cv,
                                  fict_prev2,
                                  v1,
                                  LARGER);
    new_he = new_he->opposite();
  }
  else if (fict_prev1 != NULL && fict_prev2 == NULL)
  {
    // v1 is at infinity and v2 represents a valid point.
    new_he = _insert_from_vertex (cv,
                                  fict_prev1,
                                  v2,
                                  SMALLER);
  }
  else
  {
    // Both vertices are at infinity. Note that we may create a new unbounded
    // face.
    bool        new_face_created = false;

    new_he = _insert_at_vertices (cv,
                                  fict_prev1, 
                                  fict_prev2,
                                  SMALLER,
                                  new_face_created,
                                  true);
  
    if (new_face_created)
    {
      CGAL_assertion (! new_he->is_on_hole());
      _relocate_in_new_face (new_he);
    }
  }

  // Return a handle to the new halfedge directed from left to right.
  return (Halfedge_handle (new_he));
}

//-----------------------------------------------------------------------------
// Insert an unbounded x-monotone curve into the arrangement inside a given
// unbounded face, given a fictitious halfedge(s) that contains the unbounded
// end(s) of the curve.
//
template<class Traits, class Dcel>
typename Arrangement_2<Traits,Dcel>::Halfedge_handle 
Arrangement_2<Traits,Dcel>::insert_in_face_interior
    (const X_monotone_curve_2& cv, 
     Halfedge_handle fict_he1, Halfedge_handle fict_he2)
{
  CGAL_precondition_msg (fict_he1->is_fictitious() && 
                         (fict_he2 == Halfedge_handle() ||
                          fict_he2->is_fictitious()),
                         "The given halfedge(s) must be fictitious.");

  // Check which one of cv's ends lie at infinity (both may lie there).
  const Boundary_type   inf_x1 = traits->boundary_in_x_2_object()(cv, MIN_END);
  const Boundary_type   inf_y1 = traits->boundary_in_y_2_object()(cv, MIN_END);
  const Boundary_type   inf_x2 = traits->boundary_in_x_2_object()(cv, MAX_END);
  const Boundary_type   inf_y2 = traits->boundary_in_y_2_object()(cv, MAX_END);
  DHalfedge            *new_he;

  if (inf_x1 != NO_BOUNDARY || inf_y1 != NO_BOUNDARY)
  {
    // Split the fictitious edge and create a vertex at infinity that
    // represents cv's left end.
    DHalfedge    *fict_prev1 = _halfedge (fict_he1);

    CGAL_precondition_code (
      bool   eq_source; 
      bool   eq_target;
    );
    CGAL_precondition_msg 
      (_is_on_fictitious_edge (cv, MIN_END, inf_x1, inf_y1, fict_prev1,
                               eq_source, eq_target) &&
       ! eq_source && ! eq_target,
       "The given halfedge must contain the unbounded left end.");

    CGAL_precondition_msg
      ((inf_x2 == NO_BOUNDARY && inf_y2 == NO_BOUNDARY &&
        fict_he2 == Halfedge_handle()) ||
       (_is_on_fictitious_edge (cv, MAX_END, inf_x2, inf_y2,
                                _halfedge (fict_he2),
                                eq_source, eq_target) &&
        ! eq_source && ! eq_target),
       "The given halfedge must contain the unbounded right end.");

    fict_prev1 = _split_fictitious_edge (fict_prev1,
                                         inf_x1, inf_y1);

    if (inf_x2 != NO_BOUNDARY || inf_y2 != NO_BOUNDARY)
    {
      // Split the fictitious edge and create a vertex at infinity that
      // represents cv's right end.
      DHalfedge    *fict_prev2;
      
      if (fict_he2 != Halfedge_handle())
      {
        if (fict_he2 != fict_he1)
          fict_prev2 = _halfedge (fict_he2);
        else if (fict_prev1->direction() == SMALLER)
          fict_prev2 = fict_prev1->next();
        else
          fict_prev2 = fict_prev1;
      }
      else
      {
        if (fict_prev1->direction() == SMALLER)
          fict_prev2 = fict_prev1->next();
        else
          fict_prev2 = fict_prev1;
      }

      fict_prev2 = _split_fictitious_edge (fict_prev2,
                                           inf_x2, inf_y2);

      // As both end vertices are at infinity, we may create a new unbounded
      // face.
      bool        new_face_created = false;

      new_he = _insert_at_vertices (cv,
                                    fict_prev1, 
                                    fict_prev2,
                                    SMALLER,
                                    new_face_created,
                                    true);
  
      if (new_face_created)
      {
        CGAL_assertion (! new_he->is_on_hole());
        _relocate_in_new_face (new_he);
      }
    }
    else
    {
      // Create a new vertex that corresponds to cv's right endpoint.
      DVertex       *v2 = 
        _create_vertex (traits->construct_max_vertex_2_object()(cv));

      // Insert the curve.
      new_he = _insert_from_vertex (cv,
                                    fict_prev1,
                                    v2,
                                    SMALLER);
    }
  }
  else
  {
    // In this case only the right end of cv is unbounded.
    CGAL_precondition_msg (inf_x2 != NO_BOUNDARY || inf_y2 != NO_BOUNDARY,
                           "The inserted curve must be unbounded.");

    // Create a new vertex that corresponds to cv's left endpoint.
    DVertex       *v1 = 
      _create_vertex (traits->construct_min_vertex_2_object()(cv));
    
    // Split the fictitious edge and create a vertex at infinity that
    // represents cv's right end.
    DHalfedge    *fict_prev2;
      
    if (fict_he2 != Halfedge_handle())
      fict_prev2 = _halfedge (fict_he2);
    else
      fict_prev2 = _halfedge (fict_he1);

    CGAL_precondition_code (
      bool   eq_source; 
      bool   eq_target;
    );
    CGAL_precondition_msg 
      (_is_on_fictitious_edge (cv, MAX_END, inf_x2, inf_y2, fict_prev2,
                               eq_source, eq_target) &&
       ! eq_source && ! eq_target,
       "The given halfedge must contain the unbounded right end.");
    
    fict_prev2 = _split_fictitious_edge (fict_prev2,
                                         inf_x2, inf_y2);

    // Insert the curve.
    new_he = _insert_from_vertex (cv,
                                  fict_prev2,
                                  v1,
                                  LARGER);
    new_he = new_he->opposite();
  }

  // Return a handle to the new halfedge directed from left to right.
  return (Halfedge_handle (new_he));
}

//-----------------------------------------------------------------------------
// Insert an x-monotone curve into the arrangement, such that its left 
// endpoint corresponds to a given arrangement vertex.
//
template<class Traits, class Dcel>
typename Arrangement_2<Traits,Dcel>::Halfedge_handle
Arrangement_2<Traits,Dcel>::insert_from_left_vertex
    (const X_monotone_curve_2& cv, 
     Vertex_handle v)
{
  CGAL_precondition_msg
    (! v->is_at_infinity(),
     "The input vertex must not lie at infinity.");
  CGAL_precondition_msg 
    (traits->equal_2_object() (v->point(), 
                               traits->construct_min_vertex_2_object()(cv)),
     "The input vertex should be the left curve endpoint.");

  // Check if cv's right end lies at infinity. In case it is a regular
  // point, create a normal vertex that correspond to this point.
  const Boundary_type  inf_x2 = traits->boundary_in_x_2_object()(cv, MAX_END);
  const Boundary_type  inf_y2 = traits->boundary_in_y_2_object()(cv, MAX_END);
  DVertex             *v2 = NULL;
  DHalfedge           *fict_prev2 = NULL;

  if (inf_x2 == NO_BOUNDARY && inf_y2 == NO_BOUNDARY)
  {
    // The curve has a valid right endpoint: Create a new vertex associated
    // with the curve's right endpoint.
    v2 = _create_vertex (traits->construct_max_vertex_2_object()(cv));
  }

  // Check if the given vertex, corresponding to the left endpoint,
  // is isolated.
  if (v->is_isolated())
  {
    // The given vertex is an isolated one: We should in fact insert the curve
    // in the interior of the face containing this vertex.
    DVertex    *v1 = _vertex (v);
    DIso_vert  *iv = v1->isolated_vertex();
    DFace      *p_f = iv->face();

    // If the vertex that corresponds to cv's right end lies at infinity,
    // create it now.
    if (v2 == NULL)
    {    
      // Locate a halfedge of f's outer CCB such that contains cv's right end
      // in its range.
      CGAL_assertion (p_f->is_unbounded());

      fict_prev2 = _locate_along_ccb (p_f, cv, MAX_END,
                                      inf_x2, inf_y2);

      CGAL_assertion (fict_prev2 != NULL);

      // Split the fictitious edge and create a vertex at infinity that
      // represents cv's right end.
      fict_prev2 = _split_fictitious_edge (fict_prev2,
                                           inf_x2, inf_y2);
    }

    // Remove the isolated vertex v1, as it will not be isolated any more.
    p_f->erase_isolated_vertex (iv->iterator());
    dcel.delete_isolated_vertex (iv);

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
  DHalfedge  *prev1 = _locate_around_vertex (_vertex (v), cv, MIN_END);

  CGAL_assertion_msg
    (prev1 != NULL,
     "The inserted curve should not exist in the arrangement.");

  // If the vertex that corresponds to cv's right end lies at infinity,
  // create it now.
  if (v2 == NULL)
  {    
    // Locate a halfedge along the outer CCB of the incident face of prev1
    // that contains cv's right end in its range.
    DFace      *uf = prev1->is_on_hole() ? prev1->hole()->face() :
                                           prev1->face();
    CGAL_assertion (uf->is_unbounded());

    fict_prev2 = _locate_along_ccb (uf, cv, MAX_END,
                                    inf_x2, inf_y2);

    CGAL_assertion (fict_prev2 != NULL);

    // Split the fictitious edge and create a vertex at infinity that
    // represents cv's right end.
    fict_prev2 = _split_fictitious_edge (fict_prev2,
                                         inf_x2, inf_y2);
  }

  // Perform the insertion (note that we know that prev1->vertex is smaller
  // than v2).
  DHalfedge  *new_he;

  if (fict_prev2 == NULL)
  {
    new_he = _insert_from_vertex (cv,
                                  prev1, v2,
                                  SMALLER);
  }
  else
  {
    // v2 lies at infinity. Note that we may create a new unbounded face.
    bool        new_face_created = false;

    new_he = _insert_at_vertices (cv,
                                  prev1, 
                                  fict_prev2,
                                  SMALLER,
                                  new_face_created,
                                  true);
  
    if (new_face_created)
    {
      CGAL_assertion (! new_he->is_on_hole());
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
template<class Traits, class Dcel>
typename Arrangement_2<Traits,Dcel>::Halfedge_handle
Arrangement_2<Traits,Dcel>::insert_from_left_vertex
    (const X_monotone_curve_2& cv,
     Halfedge_handle prev)
{
  CGAL_precondition_msg
    (! prev->target()->is_at_infinity(),
     "The input vertex must not lie at infinity.");
  CGAL_precondition_msg
    (traits->equal_2_object() (prev->target()->point(),
                               traits->construct_min_vertex_2_object()(cv)),
     "The input halfedge's target should be the left curve endpoint.");
  CGAL_precondition_msg
    (_locate_around_vertex (_vertex(prev->target()),
                            cv, MIN_END) == _halfedge(prev),
     "In the clockwise order of curves around the vertex, "
     " cv must succeed the curve of prev.");

  // Check if cv's right end lies at infinity, and create a vertex v2 that
  // corresponds to this end.
  const Boundary_type  inf_x2 = traits->boundary_in_x_2_object()(cv, MAX_END);
  const Boundary_type  inf_y2 = traits->boundary_in_y_2_object()(cv, MAX_END);
  DVertex             *v2 = NULL;
  DHalfedge           *fict_prev2 = NULL;
  DHalfedge           *prev1 = _halfedge (prev);

  if (inf_x2 == NO_BOUNDARY && inf_y2 == NO_BOUNDARY)
  {
    // The curve has a valid right endpoint: Create a new vertex associated
    // with the curve's right endpoint.
    v2 = _create_vertex (traits->construct_max_vertex_2_object()(cv));
  }
  else
  {    
    // Locate a halfedge along the outer CCB of prev's incident face that
    // contains cv's right end in its range.
    DFace      *uf = prev1->is_on_hole() ? prev1->hole()->face() :
                                           prev1->face();
    CGAL_precondition (uf->is_unbounded());

    fict_prev2 = _locate_along_ccb (uf, cv, MAX_END,
                                    inf_x2, inf_y2);

    CGAL_assertion (fict_prev2 != NULL);

    // Split the fictitious edge and create a vertex at infinity that
    // represents cv's right end.
    fict_prev2 = _split_fictitious_edge (fict_prev2,
                                         inf_x2, inf_y2);
  }

  // Perform the insertion (note that we know that prev1->vertex is smaller
  // than v2).
  DHalfedge  *new_he;

  if (fict_prev2 == NULL)
  {
    new_he = _insert_from_vertex (cv,
                                  prev1, v2,
                                  SMALLER);
  }
  else
  {
    // v2 lies at infinity. Note that we may create a new unbounded face.
    bool        new_face_created = false;

    new_he = _insert_at_vertices (cv,
                                  prev1, 
                                  fict_prev2,
                                  SMALLER,
                                  new_face_created,
                                  true);
  
    if (new_face_created)
    {
      CGAL_assertion (! new_he->is_on_hole());
      _relocate_in_new_face (new_he);
    }
  }

  // Return a handle to the halfedge directed toward the new vertex v2.
  return (Halfedge_handle (new_he));
}

//-----------------------------------------------------------------------------
// Insert an x-monotone curve into the arrangement, such that its left
// endpoints corresponds to a given arrangement vertex (given the exact
// place for the curve in the circular list around this vertex), and whose
// right end is unbounded and lies on a given fictitious edge.
//
template<class Traits, class Dcel>
typename Arrangement_2<Traits,Dcel>::Halfedge_handle 
Arrangement_2<Traits,Dcel>::insert_from_left_vertex
    (const X_monotone_curve_2& cv,
     Halfedge_handle prev,
     Halfedge_handle fict_he)
{
  CGAL_precondition_msg
    (! prev->target()->is_at_infinity(),
     "The input vertex must not lie at infinity.");
  CGAL_precondition_msg
    (traits->equal_2_object() (prev->target()->point(),
                               traits->construct_min_vertex_2_object()(cv)),
     "The input halfedge's target should be the left curve endpoint.");
  CGAL_precondition_msg
    (_locate_around_vertex (_vertex(prev->target()),
                            cv, MIN_END) == _halfedge(prev),
     "In the clockwise order of curves around the vertex, "
     " cv must succeed the curve of prev.");

  DHalfedge        *prev1 = _halfedge (prev);

  // Check that cv's right end lies at infinity, and create a vertex v2 that
  // corresponds to this end.
  const Boundary_type  inf_x2 = traits->boundary_in_x_2_object()(cv, MAX_END);
  const Boundary_type  inf_y2 = traits->boundary_in_y_2_object()(cv, MAX_END);
  DHalfedge           *fict_prev2 = _halfedge (fict_he);

  CGAL_precondition_code (
    bool   eq_source; 
    bool   eq_target;
  );
  CGAL_precondition_msg 
    (_is_on_fictitious_edge (cv, MAX_END, inf_x2, inf_y2, fict_prev2,
                             eq_source, eq_target) &&
     ! eq_source && ! eq_target,
     "The given halfedge must contain the unbounded right end.");

  fict_prev2 = _split_fictitious_edge (fict_prev2,
                                       inf_x2, inf_y2);

  // Insert the curve and create an edge connecting the the two vertices.
  // Note that we may create a new unbounded face.
  bool        new_face_created = false;
  DHalfedge  *new_he = _insert_at_vertices (cv,
                                            prev1, 
                                            fict_prev2,
                                            SMALLER,
                                            new_face_created,
                                            true);
  
  if (new_face_created)
  {
    CGAL_assertion (! new_he->is_on_hole());
    _relocate_in_new_face (new_he);
  }

  // Return a handle to the halfedge directed toward the new vertex v2.
  return (Halfedge_handle (new_he));
}

//-----------------------------------------------------------------------------
// Insert an x-monotone curve into the arrangement, such that its right
// endpoint corresponds to a given arrangement vertex.
//
template<class Traits, class Dcel>
typename Arrangement_2<Traits,Dcel>::Halfedge_handle
Arrangement_2<Traits,Dcel>::insert_from_right_vertex
    (const X_monotone_curve_2& cv,
     Vertex_handle v)
{
  CGAL_precondition_msg
    (! v->is_at_infinity(),
     "The input vertex must not lie at infinity.");
  CGAL_precondition_msg
    (traits->equal_2_object() (v->point(),
                               traits->construct_max_vertex_2_object()(cv)),
     "The input vertex should be the right curve endpoint.");

  // Check if cv's left end lies at infinity. In case it is a regular
  // point, create a normal vertex that correspond to this point.
  const Boundary_type  inf_x1 = traits->boundary_in_x_2_object()(cv, MIN_END);
  const Boundary_type  inf_y1 = traits->boundary_in_y_2_object()(cv, MIN_END);
  DVertex             *v1 = NULL;
  DHalfedge           *fict_prev1 = NULL;

  if (inf_x1 == NO_BOUNDARY && inf_y1 == NO_BOUNDARY)
  {
    // The curve has a valid left endpoint: Create a new vertex associated
    // with the curve's left endpoint.
    v1 = _create_vertex (traits->construct_min_vertex_2_object()(cv));
  }

  // Check if the given vertex is isolated.
  if (v->is_isolated())
  {
    // The given vertex is an isolated one: We should in fact insert the curve
    // in the interior of the face containing this vertex.
    DVertex    *v2 = _vertex (v);
    DIso_vert  *iv = v2->isolated_vertex();
    DFace      *p_f = iv->face();

    // If the vertex that corresponds to cv's left end lies at infinity,
    // create it now.
    if (v1 == NULL)
    {    
      // Locate a halfedge of f's outer CCB such that contains cv's left end
      // in its range.
      CGAL_assertion (p_f->is_unbounded());

      fict_prev1 = _locate_along_ccb (p_f, cv, MIN_END,
                                      inf_x1, inf_y1);

      CGAL_assertion (fict_prev1 != NULL);

      // Split the fictitious edge and create a vertex at infinity that
      // represents cv's left end.
      fict_prev1 = _split_fictitious_edge (fict_prev1,
                                           inf_x1, inf_y1);
    }

    // Remove the isolated vertex v2, as it will not be isolated any more.
    p_f->erase_isolated_vertex (iv->iterator());
    dcel.delete_isolated_vertex (iv);

    // Create the edge connecting the two vertices (note that we know that
    // v1 is smaller than v2).
    DHalfedge  *new_he;

    if (fict_prev1 == NULL)
    {
      new_he = _insert_in_face_interior (cv, p_f, v1, v2,
                                         SMALLER);
    }
    else
    {
      new_he = _insert_from_vertex (cv,
                                    fict_prev1, v2,
                                    SMALLER);
    }
    
    // Return a handle to the new halfedge whose target is the new vertex v1.
    return (Halfedge_handle (new_he->opposite()));
  }

  // Go over the incident halfedges around v and find the halfedge after
  // which the new curve should be inserted.
  DHalfedge  *prev2 = _locate_around_vertex (_vertex (v), cv, MAX_END);

  CGAL_assertion_msg
    (prev2 != NULL,
     "The inserted curve should not exist in the arrangement.");
  
  // If the vertex that corresponds to cv's left end lies at infinity,
  // create it now.
  if (v1 == NULL)
  {    
    // Locate a halfedge along the outer CCB of the incident face of prev2
    // that contains cv's left end in its range.
    DFace      *uf = prev2->is_on_hole() ? prev2->hole()->face() :
                                           prev2->face();
    CGAL_assertion (uf->is_unbounded());

    fict_prev1 = _locate_along_ccb (uf, cv, MIN_END,
                                    inf_x1, inf_y1);

    CGAL_assertion (fict_prev1 != NULL);

    // Split the fictitious edge and create a vertex at infinity that
    // represents cv's left end.
    fict_prev1 = _split_fictitious_edge (fict_prev1,
                                         inf_x1, inf_y1);
  }

  // Perform the insertion (note that we know that prev2->vertex is larger
  // than v1).
  DHalfedge  *new_he;

  if (fict_prev1 == NULL)
  {
    new_he = _insert_from_vertex (cv,
                                  prev2, v1,
                                  LARGER);
  }
  else
  {
    // v1 lies at infinity. Note that we may create a new unbounded face.
    bool        new_face_created = false;

    new_he = _insert_at_vertices (cv,
                                  prev2, 
                                  fict_prev1,
                                  LARGER,
                                  new_face_created,
                                  true);
  
    if (new_face_created)
    {
      CGAL_assertion (! new_he->is_on_hole());
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
template<class Traits, class Dcel>
typename Arrangement_2<Traits,Dcel>::Halfedge_handle
Arrangement_2<Traits,Dcel>::insert_from_right_vertex
    (const X_monotone_curve_2& cv,
     Halfedge_handle prev)
{
  CGAL_precondition_msg
    (! prev->target()->is_at_infinity(),
     "The input vertex must not lie at infinity.");

  CGAL_precondition_msg
    (traits->equal_2_object() (prev->target()->point(),
                               traits->construct_max_vertex_2_object()(cv)),
     "The input halfedge's target should be the right curve endpoint.");

  CGAL_precondition_msg
    (_locate_around_vertex (_vertex(prev->target()),
                            cv, MAX_END) == _halfedge(prev),
     "In the clockwise order of curves around the vertex, "
     "cv must succeed the curve of prev.");

  // Check if cv's left end lies at infinity, and create a vertex v1 that
  // corresponds to this end.
  const Boundary_type  inf_x1 = traits->boundary_in_x_2_object()(cv, MIN_END);
  const Boundary_type  inf_y1 = traits->boundary_in_y_2_object()(cv, MIN_END);
  DVertex             *v1 = NULL;
  DHalfedge           *fict_prev1 = NULL;
  DHalfedge           *prev2 = _halfedge (prev);

  if (inf_x1 == NO_BOUNDARY && inf_y1 == NO_BOUNDARY)
  {
    // The curve has a valid left endpoint: Create a new vertex associated
    // with the curve's left endpoint.
    v1 = _create_vertex (traits->construct_min_vertex_2_object()(cv));
  }
  else
  {    
    // Locate a halfedge along the outer CCB of prev's incident face that
    // contains cv's left end in its range.
    DFace      *uf = prev2->is_on_hole() ? prev2->hole()->face() :
                                           prev2->face();
    CGAL_precondition (uf->is_unbounded());

    fict_prev1 = _locate_along_ccb (uf, cv, MIN_END,
                                    inf_x1, inf_y1);

    CGAL_assertion (fict_prev1 != NULL);

    // Split the fictitious edge and create a vertex at infinity that
    // represents cv's left end.
    fict_prev1 = _split_fictitious_edge (fict_prev1,
                                         inf_x1, inf_y1);
  }

  // Perform the insertion (note that we know that prev2->vertex is larger
  // than v1).
  DHalfedge  *new_he;

  if (fict_prev1 == NULL)
  {
    new_he = _insert_from_vertex (cv,
                                  prev2, v1,
                                  LARGER);
  }
  else
  {
    // v1 lies at infinity. Note that we may create a new unbounded face.
    bool        new_face_created = false;

    new_he = _insert_at_vertices (cv,
                                  prev2,
                                  fict_prev1, 
                                  LARGER,
                                  new_face_created,
                                  true);
  
    if (new_face_created)
    {
      CGAL_assertion (! new_he->is_on_hole());
      _relocate_in_new_face (new_he);
    }
  }

  // Return a handle to the halfedge directed toward the new vertex v1.
  return (Halfedge_handle (new_he));
}

//-----------------------------------------------------------------------------
// Insert an x-monotone curve into the arrangement, such that its right
// endpoint corresponds to a given arrangement vertex (given the exact
// place for the curve in the circular list around this vertex), and whose
// left end is unbounded and lies on a given fictitious edge.
//
template<class Traits, class Dcel>
typename Arrangement_2<Traits,Dcel>::Halfedge_handle 
Arrangement_2<Traits,Dcel>::insert_from_right_vertex
    (const X_monotone_curve_2& cv,
     Halfedge_handle prev,
     Halfedge_handle fict_he)
{
  CGAL_precondition_msg
    (! prev->target()->is_at_infinity(),
     "The input vertex must not lie at infinity.");
  CGAL_precondition_msg
    (traits->equal_2_object() (prev->target()->point(),
                               traits->construct_max_vertex_2_object()(cv)),
     "The input halfedge's target should be the right curve endpoint.");
  CGAL_precondition_msg
    (_locate_around_vertex (_vertex(prev->target()),
                            cv, MAX_END) == _halfedge(prev),
     "In the clockwise order of curves around the vertex, "
     " cv must succeed the curve of prev.");

  DHalfedge        *prev2 = _halfedge (prev);

  // Check that cv's left end lies at infinity, and create a vertex v1 that
  // corresponds to this end.
  const Boundary_type  inf_x1 = traits->boundary_in_x_2_object()(cv, MIN_END);
  const Boundary_type  inf_y1 = traits->boundary_in_y_2_object()(cv, MIN_END);
  DHalfedge           *fict_prev1 = _halfedge (fict_he);

  CGAL_precondition_code (
    bool   eq_source; 
    bool   eq_target;
  );
  CGAL_precondition_msg 
    (_is_on_fictitious_edge (cv, MIN_END, inf_x1, inf_y1, fict_prev1,
                             eq_source, eq_target) &&
     !eq_source && !eq_target,
     "The given halfedge must contain the unbounded left end.");

  fict_prev1 = _split_fictitious_edge (fict_prev1,
                                       inf_x1, inf_y1);

  // Insert the curve and create an edge connecting the the two vertices.
  // Note that we may create a new unbounded face.
  bool        new_face_created = false;
  DHalfedge  *new_he = _insert_at_vertices (cv,
                                            prev2, 
                                            fict_prev1,
                                            LARGER,
                                            new_face_created,
                                            true);
  
  if (new_face_created)
  {
    CGAL_assertion (! new_he->is_on_hole());
    _relocate_in_new_face (new_he);
  }

  // Return a handle to the halfedge directed toward the new vertex v2.
  return (Halfedge_handle (new_he));
}

//-----------------------------------------------------------------------------
// Insert an x-monotone curve into the arrangement, such that both its
// endpoints corresponds to a given arrangement vertex.
//
template<class Traits, class Dcel>
typename Arrangement_2<Traits,Dcel>::Halfedge_handle
Arrangement_2<Traits,Dcel>::insert_at_vertices (const X_monotone_curve_2& cv,
                                                Vertex_handle v1,
                                                Vertex_handle v2)
{
  CGAL_precondition_msg
    (! v1->is_at_infinity() && ! v2->is_at_infinity(),
     "The input vertices must not lie at infinity.");

  Curve_end    ind1;
  Curve_end    ind2;

  if (traits->equal_2_object() (v1->point(),
                                traits->construct_min_vertex_2_object()(cv)))
  {
    CGAL_precondition_msg 
      (traits->equal_2_object()(v2->point(),
                                traits->construct_max_vertex_2_object()(cv)),
       "The input vertices should match the curve endpoints.");

    ind1 = MIN_END;
    ind2 = MAX_END;
  }
  else
  {
    CGAL_precondition_msg
      (traits->equal_2_object()(v2->point(),
                                traits->construct_min_vertex_2_object()(cv)) &&
       traits->equal_2_object()(v1->point(),
                                traits->construct_max_vertex_2_object()(cv)),
       "The input vertices should match the curve endpoints.");

    ind1 = MAX_END;
    ind2 = MIN_END;
  }

  // Check whether one of the vertices is isolated.
  if (v1->is_isolated())
  {
    // Get the face containing the isolated vertex v1.
    DVertex    *p_v1 = _vertex (v1);
    DIso_vert  *iv1 = p_v1->isolated_vertex();
    DFace      *f1 = iv1->face();

    // Remove the isolated vertex v1, as it will not be isolated any more.
    f1->erase_isolated_vertex (iv1->iterator());
    dcel.delete_isolated_vertex (iv1);

    if (v2->is_isolated())
    {
      // Both end-vertices are isolated. Make sure they are contained inside
      // the same face.
      DVertex    *p_v2 = _vertex (v2);
      DIso_vert  *iv2 = p_v2->isolated_vertex();

      CGAL_assertion_code (
        DFace      *f2 = iv2->face();
      );

      CGAL_assertion_msg
        (f1 == f2,
         "The inserted curve should not intersect the existing arrangement.");

      // Remove the isolated vertex v2, as it will not be isolated any more.
      f1->erase_isolated_vertex (iv2->iterator());
      dcel.delete_isolated_vertex (iv2);

      // Perform the insertion.
      Comparison_result  res = traits->compare_xy_2_object() (v1->point(),
                                                              v2->point());
      CGAL_assertion (res != EQUAL);

      DHalfedge  *new_he = _insert_in_face_interior (cv, f1,
                                                     p_v1, p_v2,
                                                     res);

      return (Halfedge_handle (new_he));
    }

    // Go over the incident halfedges around v2 and find the halfedge after
    // which the new curve should be inserted.
    DHalfedge  *prev2 = _locate_around_vertex (_vertex (v2), cv, ind2);

    CGAL_assertion_msg
      (prev2 != NULL,
       "The inserted curve should not exist in the arrangement.");

    CGAL_assertion_msg
      (f1 == _face(Halfedge_handle (prev2)->face()),
       "The inserted curve should not intersect the existing arrangement.");

    // Perform the insertion. Note that the returned halfedge is directed
    // toward v1 (and not toward v2), so we return its twin.
    Comparison_result  res = traits->compare_xy_2_object() (v2->point(),
                                                            v1->point());
    CGAL_assertion (res != EQUAL);

    DHalfedge  *new_he = _insert_from_vertex (cv,
                                              prev2, p_v1,
                                              res);

    return (Halfedge_handle (new_he->opposite()));
  }
  else if (v2->is_isolated())
  {
    // Get the face containing the isolated vertex v2.
    DVertex    *p_v2 = _vertex (v2);
    DIso_vert  *iv2 = p_v2->isolated_vertex();
    DFace      *f2 = iv2->face();

    // Remove the isolated vertex v2, as it will not be isolated any more.
    f2->erase_isolated_vertex (iv2->iterator());
    dcel.delete_isolated_vertex (iv2);

    // Go over the incident halfedges around v1 and find the halfedge after
    // which the new curve should be inserted.
    DHalfedge  *prev1 = _locate_around_vertex (_vertex (v1), cv, ind1);

    CGAL_assertion_msg
      (prev1 != NULL,
       "The inserted curve should not exist in the arrangement.");

    CGAL_assertion_msg
      (f2 == _face (Halfedge_handle (prev1)->face()),
       "The inserted curve should not intersect the existing arrangement.");

    // Perform the insertion.
    Comparison_result  res = traits->compare_xy_2_object() (v1->point(),
                                                            v2->point());
    CGAL_assertion (res != EQUAL);

    DHalfedge  *new_he = _insert_from_vertex (cv,
                                              prev1, p_v2,
                                              res);

    return (Halfedge_handle (new_he));
  }

  // Go over the incident halfedges around v1 and v2 and find the two
  // halfedges after which the new curve should be inserted, respectively.
  DHalfedge  *prev1 = _locate_around_vertex (_vertex (v1), cv, ind1);
  DHalfedge  *prev2 = _locate_around_vertex (_vertex (v2), cv, ind2);

  CGAL_assertion_msg
    (prev1 != NULL && prev2 != NULL,
     "The inserted curve should not exist in the arrangement.");

  // Perform the insertion.
  return (insert_at_vertices (cv,
                              Halfedge_handle(prev1),
                              Halfedge_handle(prev2)));
}

//-----------------------------------------------------------------------------
// Insert an x-monotone curve into the arrangement, such that both its
// endpoints correspond to given arrangement vertices, given the exact
// place for the curve in one of the circular lists around a vertex.
//
template<class Traits, class Dcel>
typename Arrangement_2<Traits,Dcel>::Halfedge_handle
Arrangement_2<Traits,Dcel>::insert_at_vertices (const X_monotone_curve_2& cv,
                                                Halfedge_handle prev1,
                                                Vertex_handle v2)
{
  CGAL_precondition_msg
    (! prev1->target()->is_at_infinity() && ! v2->is_at_infinity(),
     "The input vertices must not lie at infinity.");

  Curve_end    ind2;

  if (traits->equal_2_object() (prev1->target()->point(),
                                traits->construct_min_vertex_2_object()(cv)))
  {
    CGAL_precondition_msg 
      (traits->equal_2_object()(v2->point(),
                                traits->construct_max_vertex_2_object()(cv)),
       "The input vertices should match the curve endpoints.");

    ind2 = MAX_END;
  }
  else
  {
    CGAL_precondition_msg
      (traits->equal_2_object()(v2->point(),
                                traits->construct_min_vertex_2_object()(cv)) &&
       traits->equal_2_object()(prev1->target()->point(),
                                traits->construct_max_vertex_2_object()(cv)),
       "The input vertices should match the curve endpoints.");

    ind2 = MIN_END;
  }

  // Check whether v2 is an isolated vertex.
  if (v2->is_isolated())
  {
    // Get the face containing the isolated vertex v2.
    DVertex    *p_v2 = _vertex (v2);
    DIso_vert  *iv2 = p_v2->isolated_vertex();
    DFace      *f2 = iv2->face();

    CGAL_assertion_msg
      (f2 == _face (prev1->face()),
       "The inserted curve should not intersect the existing arrangement.");

    // Remove the isolated vertex v2, as it will not be isolated any more.
    f2->erase_isolated_vertex (iv2->iterator());
    dcel.delete_isolated_vertex (iv2);

    // Perform the insertion.
    Comparison_result  res = traits->compare_xy_2_object() 
                                             (prev1->target()->point(),
                                              v2->point());
    CGAL_assertion (res != EQUAL);

    DHalfedge  *new_he = _insert_from_vertex (cv,
                                              _halfedge (prev1), p_v2,
                                              res);

    return (Halfedge_handle (new_he));
  }

  // Go over the incident halfedges around v2 and find the halfedge after
  // which the new curve should be inserted.
  DHalfedge  *prev2 = _locate_around_vertex (_vertex (v2), cv, ind2);

  CGAL_assertion_msg
    (prev2 != NULL,
     "The inserted curve should not exist in the arrangement.");

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
template<class Traits, class Dcel>
typename Arrangement_2<Traits,Dcel>::Halfedge_handle
Arrangement_2<Traits,Dcel>::insert_at_vertices (const X_monotone_curve_2& cv,
                                                Halfedge_handle prev1,
                                                Halfedge_handle prev2)
{
  CGAL_precondition_msg
    (! prev1->target()->is_at_infinity() &&
     ! prev2->target()->is_at_infinity(),
     "The input vertices must not lie at infinity.");

  CGAL_precondition_msg
    ((traits->equal_2_object()(prev1->target()->point(),
                               traits->construct_min_vertex_2_object()(cv)) &&
      traits->equal_2_object()(prev2->target()->point(),
                               traits->construct_max_vertex_2_object()(cv))) ||
     (traits->equal_2_object()(prev2->target()->point(),
                               traits->construct_min_vertex_2_object()(cv)) &&
      traits->equal_2_object()(prev1->target()->point(),
                               traits->construct_max_vertex_2_object()(cv))),
     "The input halfedges' target points should match the curve endpoints.");

  // Check if e1 and e2 are on the same connected component.
  DHalfedge  *p_prev1 = _halfedge (prev1);
  DHalfedge  *p_prev2 = _halfedge (prev2);
  DHole      *hole1 = (p_prev1->is_on_hole()) ? p_prev1->hole() : NULL;
  DHole      *hole2 = (p_prev2->is_on_hole()) ? p_prev2->hole() : NULL;
  bool        prev1_before_prev2 = true;

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
    // To do this, we check whether prev1 lies inside the new face we are
    // about to create (or alternatively, whether prev2 does not lie inside
    // this new face).
    const unsigned int  dist1 = _halfedge_distance (p_prev1, p_prev2);
    const unsigned int  dist2 = _halfedge_distance (p_prev2, p_prev1);

    prev1_before_prev2 = (dist1 > dist2) ?
      (_is_inside_new_face (p_prev1, p_prev2, cv)) :
      (! _is_inside_new_face (p_prev2, p_prev1, cv));
  }

  // Compare the points associated with the end-vertices of the edge we are
  // about to create.
  Comparison_result    res;

  if (prev1_before_prev2)
    res = traits->compare_xy_2_object() (p_prev1->vertex()->point(),
                                         p_prev2->vertex()->point());
  else
    res = traits->compare_xy_2_object() (p_prev2->vertex()->point(),
                                         p_prev1->vertex()->point());

  // Perform the insertion.
  bool        new_face_created = false;
  DHalfedge  *new_he = (prev1_before_prev2) ?
    _insert_at_vertices (cv, p_prev1, p_prev2, res, new_face_created, false) :
    _insert_at_vertices (cv, p_prev2, p_prev1, res, new_face_created, false);

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
template<class Traits, class Dcel>
typename Arrangement_2<Traits,Dcel>::Vertex_handle
Arrangement_2<Traits,Dcel>::modify_vertex (Vertex_handle vh,
                                           const Point_2& p)
{
  CGAL_precondition_msg
    (! vh->is_at_infinity(),
     "The modified vertex must not lie at infinity.");
  CGAL_precondition_msg (traits->equal_2_object() (vh->point(), p),
                         "The new point is different from the current one.");

  // Modify the vertex.
  _modify_vertex (_vertex (vh), p);

  // Return a handle to the modified vertex.
  return (vh);
}

//-----------------------------------------------------------------------------
// Remove an isolated vertex from the interior of a given face.
//
template<class Traits, class Dcel>
typename Arrangement_2<Traits,Dcel>::Face_handle
Arrangement_2<Traits,Dcel>::remove_isolated_vertex (Vertex_handle v)
{
  CGAL_precondition (v->is_isolated());

  // Get the face containing v.
  DVertex      *p_v = _vertex (v);
  DIso_vert    *iv = p_v->isolated_vertex();
  DFace        *p_f = iv->face();
  Face_handle   f = Face_handle (p_f);
  
  // Notify the observers that we are abount to remove a vertex.
  _notify_before_remove_vertex (v);

  // Remove the isolated vertex from the face that contains it.
  p_f->erase_isolated_vertex (iv->iterator());
  dcel.delete_isolated_vertex (iv);

  // Delete the vertex.
  _delete_point (p_v->point());
  dcel.delete_vertex (p_v);

  // Notify the observers that the vertex has been removed.
  _notify_after_remove_vertex ();

  // Return a handle for the face that used to contain the deleted vertex.
  return (f);
}

//-----------------------------------------------------------------------------
// Replace the x-monotone curve associated with the given edge.
//
template<class Traits, class Dcel>
typename Arrangement_2<Traits,Dcel>::Halfedge_handle
Arrangement_2<Traits,Dcel>::modify_edge (Halfedge_handle e,
                                         const X_monotone_curve_2& cv)
{
  CGAL_precondition_msg (! e->is_fictitious(),
                         "The edge must be a valid one.");
  CGAL_precondition_msg (traits->equal_2_object() (e->curve(), cv),
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
template<class Traits, class Dcel>
typename Arrangement_2<Traits,Dcel>::Halfedge_handle
Arrangement_2<Traits,Dcel>::split_edge (Halfedge_handle e,
                                        const X_monotone_curve_2& cv1,
                                        const X_monotone_curve_2& cv2)
{
  CGAL_precondition_msg (! e->is_fictitious(),
                         "The edge must be a valid one.");

  // Get the split halfedge and its twin, its source and target.
  DHalfedge       *he1 = _halfedge (e);
  DHalfedge       *he2 = he1->opposite();
  DVertex         *source = he2->vertex();
  CGAL_precondition_code (
    DVertex         *target = he1->vertex();
  );

  // Determine the point where we split the halfedge. We also determine which
  // curve should be associated with he1 (and he2), which is the curve who
  // has an endpoint that equals e's source, and which should be associated
  // with the new pair of halfedges we are about to split (the one who has
  // an endpoint which equals e's target).
  Point_2                    split_pt;
  const X_monotone_curve_2  *p_cv1 = NULL;
  const X_monotone_curve_2  *p_cv2 = NULL;

  if (traits->boundary_in_x_2_object()(cv1, MAX_END) == NO_BOUNDARY &&
      traits->boundary_in_y_2_object()(cv1, MAX_END) == NO_BOUNDARY)
  {
    const Point_2&  cv1_right = traits->construct_max_vertex_2_object() (cv1);

    if (traits->boundary_in_x_2_object()(cv2, MIN_END) == NO_BOUNDARY &&
        traits->boundary_in_y_2_object()(cv2, MIN_END) == NO_BOUNDARY &&
        traits->equal_2_object()(traits->construct_min_vertex_2_object()(cv2),
                                 cv1_right))
    {
      // cv1's right endpoint and cv2's left endpoint are equal, so this should
      // be the split point. Now we check whether cv1 is incident to e's source
      // and cv2 to its target, or vice versa.
      split_pt = cv1_right;

      if (_are_equal (source, cv1, MIN_END))
      {
        CGAL_precondition_msg
          (_are_equal (target, cv2, MAX_END),
           "The subcurve endpoints must match e's end vertices.");

        p_cv1 = &cv1;
        p_cv2 = &cv2;
      }
      else
      {
        CGAL_precondition_msg
          (_are_equal (source, cv2, MAX_END) &&
           _are_equal (target, cv1, MIN_END),
           "The subcurve endpoints must match e's end vertices.");

        p_cv1 = &cv2;
        p_cv2 = &cv1;      
      }
    }
  }

  if (p_cv1 == NULL && p_cv2 == NULL &&
      traits->boundary_in_x_2_object()(cv1, MIN_END) == NO_BOUNDARY &&
      traits->boundary_in_y_2_object()(cv1, MIN_END) == NO_BOUNDARY)
  {
    const Point_2&  cv1_left = traits->construct_min_vertex_2_object() (cv1);

    if (traits->boundary_in_x_2_object()(cv2, MAX_END) == NO_BOUNDARY &&
        traits->boundary_in_y_2_object()(cv2, MAX_END) == NO_BOUNDARY &&
        traits->equal_2_object()(traits->construct_max_vertex_2_object()(cv2),
                                 cv1_left))
    {
      // cv1's left endpoint and cv2's right endpoint are equal, so this should
      // be the split point. Now we check whether cv1 is incident to e's source
      // and cv2 to its target, or vice versa.
      split_pt = cv1_left;

      if (_are_equal (source, cv2, MIN_END))
      {
        CGAL_precondition_msg
          (_are_equal (target, cv1, MAX_END),
           "The subcurve endpoints must match e's end vertices.");

        p_cv1 = &cv2;
        p_cv2 = &cv1;
      }
      else
      {
        CGAL_precondition_msg
          (_are_equal (source, cv1, MAX_END) && 
           _are_equal (target, cv2, MIN_END),
           "The subcurve endpoints must match e's end vertices.");

        p_cv1 = &cv1;
        p_cv2 = &cv2;      
      }
    }
  }

  CGAL_precondition_msg (p_cv1 != NULL && p_cv2 != NULL,
                         "The two subcurves must have a common endpoint.");

  // Perform the split.
  return (Halfedge_handle (_split_edge (he1, split_pt, *p_cv1, *p_cv2)));
}

//-----------------------------------------------------------------------------
// Merge two edges to form a single edge, and associate the given x-monotone
// curve with the merged edge.
//
template<class Traits, class Dcel>
typename Arrangement_2<Traits,Dcel>::Halfedge_handle
Arrangement_2<Traits,Dcel>::merge_edge (Halfedge_handle e1,
                                        Halfedge_handle e2,
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
     "The vertex removed by the merge must not lie at infinity.");
  CGAL_precondition_msg
    (he1->next()->opposite() == he4 &&
     he4->next()->opposite() == he1,
     "The degree of the deleted vertex is greater than 2.");

  // Make sure the curve ends match the end vertices of the merged edge.
  CGAL_precondition_msg
    ((_are_equal (he2->vertex(), cv, MIN_END) && 
      _are_equal (he3->vertex(), cv, MAX_END)) ||
     (_are_equal (he3->vertex(), cv, MIN_END) && 
      _are_equal (he2->vertex(), cv, MAX_END)),
     "The endpoints of the merged curve must match the end vertices.");

  // Keep pointers to the components that contain two halfedges he3 and he2,
  // pointing at the end vertices of the merged halfedge.
  DHole       *hole1 = (he3->is_on_hole()) ? he3->hole() : NULL;
  DFace       *f1 = (hole1 == NULL) ? he3->face() : hole1->face();
  DHole       *hole2 = (he4->is_on_hole()) ? he4->hole() : NULL;
  DFace       *f2 = (hole2 == NULL) ? he4->face() : hole2->face();

  // Notify the observers that we are about to merge an edge.
  _notify_before_merge_edge (e1, e2, cv);

  // As he1 and he2 will evetually represent the merged edge, while he3 and he4
  // will be deleted, check if the deleted halfedges are represantatives of a
  // face boundary or a hole inside these faces. If so, replace he3 by he1 and
  // he4 by he2. Note that as we just change the hole representatives, we do
  // not have to notify the observers about the change.
  if (hole1 == NULL && f1->halfedge() == he3)
  {
    f1->set_halfedge (he1);
  }
  else if (hole1 != NULL)
  {
    if (*(hole1->iterator()) == he3)
    {
      f1->erase_hole (hole1->iterator());
      hole1->set_iterator (f1->add_hole (he1));
    }
  }

  if (hole2 == NULL && f2->halfedge() == he4)
  {
    f2->set_halfedge (he2);
  }
  else if (hole2 != NULL)
  {
    if (*(hole2->iterator()) == he4)
    {
      f2->erase_hole (hole2->iterator());
      hole2->set_iterator (f2->add_hole (he2));
    }
  }

  if (he3->vertex()->halfedge() == he3)
    // If he3 is the incident halfedge to its target, replace it by he1.
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
  dcel.delete_vertex (v);

  // Notify the observers that the vertex has been deleted.
  _notify_after_remove_vertex ();

  // Delete the redundant halfedge pair.
  dcel.delete_edge (he3);

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
template<class Traits, class Dcel>
typename Arrangement_2<Traits,Dcel>::Face_handle
Arrangement_2<Traits,Dcel>::remove_edge (Halfedge_handle e,
                                         bool remove_source,
                                         bool remove_target)
{
  CGAL_precondition_msg (! e->is_fictitious(),
                         "The edge must be a valid one.");

  DHalfedge   *he1 = _halfedge (e);
  DHalfedge   *he2 = he1->opposite();
  DHole       *hole1 = (he1->is_on_hole()) ? he1->hole() : NULL;
  DFace       *f1 = (hole1 == NULL) ? he1->face() : hole1->face();
  DHole       *hole2 = (he2->is_on_hole()) ? he2->hole() : NULL;
  DFace       *f2 = (hole2 == NULL) ? he2->face() : hole2->face();
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
    // In this case if a new hole will be created by the removal of he1 (and
    // its twin halfedge). We determine the halfedge (he1 or he2) that
    // points at the new hole in order to send it to _remove_hole().
    // We begin by locating the leftmost vertex along the path from he1 to he2
    // and the vertex point along the path from he2 to he1.
    const DHalfedge  *ccb1 = he1->next();
    const DHalfedge  *ccb2 = he2->next();
    const DVertex    *v_min1 = he1->vertex();
    const DVertex    *v_min2 = he2->vertex();

    do
    {
      if (_compare_vertices_xy (ccb1->vertex(), v_min1) == SMALLER)
        v_min1 = ccb1->vertex();

      ccb1 = ccb1->next();

    } while (ccb1 != he2);

    do
    {
      if (_compare_vertices_xy (ccb2->vertex(), v_min2) == SMALLER)
        v_min2 = ccb2->vertex();

      ccb2 = ccb2->next();

    } while (ccb2 != he1);

    // Compare the two leftmost points: v_min2 lies to the left of v_min1
    // if and only if he1 points at the hole we are about to create.
    if (_compare_vertices_xy (v_min2, v_min1) == SMALLER)
    {
      // he1 is directed to the new hole to be created.
      f = _remove_edge (he1, remove_source, remove_target);
    }
    else
    {
      // he2 is directed to the new hole to be created. As its source and
      // target are swapped with respect to the end-vertices of the given
      // halfedge e, we swap the roles of the two input flags
      f = _remove_edge (he2, remove_target, remove_source);
    }
  }

  return (Face_handle (f));
}

//-----------------------------------------------------------------------------
// Protected member functions (for internal use).
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Construct an empty arrangement (with its bounding rectangle).
//
template<class Traits, class Dcel>
void Arrangement_2<Traits,Dcel>::_construct_empty_arrangement ()
{
  // Create the fictitious unbounded face.
  un_face = dcel.new_face();
  un_face->set_halfedge(NULL);
  un_face->set_unbounded(true);

  // Create the four fictitious vertices corresponding to corners of the
  // bounding rectangle.
  v_bl = dcel.new_vertex();
  v_bl->set_boundary (MINUS_INFINITY, MINUS_INFINITY);

  v_tl = dcel.new_vertex();
  v_tl->set_boundary (MINUS_INFINITY, PLUS_INFINITY);

  v_br = dcel.new_vertex();
  v_br->set_boundary (PLUS_INFINITY, MINUS_INFINITY);

  v_tr = dcel.new_vertex();
  v_tr->set_boundary (PLUS_INFINITY, PLUS_INFINITY);

  // Create a four pairs of twin halfedges connecting the two vertices,
  // and link them together to form the bounding rectangle, forming a hole
  // in the fictitious face.
  //
  //                            he2
  //             v_tl (.) ----------------> (.) v_tr
  //                   ^ <------------------ 
  //                   ||                   ^|
  //               he1 ||        in_f       ||
  //                   ||                   || he3
  //                   |V                   ||
  //                     ------------------> V 
  //             v_bl (.) <---------------- (.) v_br
  //                             he4
  //
  DHalfedge           *he1 = dcel.new_edge();
  DHalfedge           *he1_t = he1->opposite();
  DHalfedge           *he2 = dcel.new_edge();
  DHalfedge           *he2_t = he2->opposite();
  DHalfedge           *he3 = dcel.new_edge();
  DHalfedge           *he3_t = he3->opposite();
  DHalfedge           *he4 = dcel.new_edge();
  DHalfedge           *he4_t = he4->opposite();
  DHole               *hole = dcel.new_hole();
  DFace               *in_f = dcel.new_face();

  hole->set_face (un_face);
  he1->set_curve (NULL);
  he2->set_curve (NULL);
  he3->set_curve (NULL);
  he4->set_curve (NULL);

  he1->set_next (he2);        he1_t->set_next (he4_t);
  he2->set_next (he3);        he4_t->set_next (he3_t);
  he3->set_next (he4);        he3_t->set_next (he2_t);
  he4->set_next (he1);        he2_t->set_next (he1_t);

  he1->set_vertex (v_tl);     he1_t->set_vertex (v_bl);
  he2->set_vertex (v_tr);     he2_t->set_vertex (v_tl);
  he3->set_vertex (v_br);     he3_t->set_vertex (v_tr);
  he4->set_vertex (v_bl);     he4_t->set_vertex (v_br);

  he1->set_hole (hole);       he1_t->set_face (in_f);
  he2->set_hole (hole);       he2_t->set_face (in_f);
  he3->set_hole (hole);       he3_t->set_face (in_f);
  he4->set_hole (hole);       he4_t->set_face (in_f);

  // Assign the incident halfedges of the two fictitious vertices.
  v_bl->set_halfedge (he1_t);
  v_tl->set_halfedge (he2_t);
  v_tr->set_halfedge (he3_t);
  v_br->set_halfedge (he4_t);

  // Set the direction of the halfedges:
  he1->set_direction (CGAL::SMALLER);
  he2->set_direction (CGAL::SMALLER);
  he3->set_direction (CGAL::LARGER);
  he4->set_direction (CGAL::LARGER);

  // Initiate a new hole inside the fictitious face.
  hole->set_iterator (un_face->add_hole (he1));

  // Set the real unbounded face, in the interior of the bounding rectangle.
  in_f->set_unbounded (true);
  in_f->set_halfedge (he2_t);

  // Mark that there are four vertices at infinity (the fictitious ones)
  // in the arrangement.
  n_inf_verts = 4;

  return;
}

//-----------------------------------------------------------------------------
// Get the curve associated with a vertex at infinity.
//
template<class Traits, class Dcel>
const typename Arrangement_2<Traits,Dcel>::X_monotone_curve_2&
Arrangement_2<Traits,Dcel>::_get_curve (const DVertex *v,
                                        Curve_end& ind) const
{
  // Go over the incident halfegdes of v until encountering the halfedge
  // associated with a valid curve (v should have three incident halfedges,
  // two of the are fictitious and one associated with a curve).
  const DHalfedge         *he = v->halfedge();

  while (he->has_null_curve())
  {
    he = he->next()->opposite();
    CGAL_assertion (he != v->halfedge());   // Guard for infinitive loops.
  }

  // The halfedge he is directed toward v, so if it is directed from left to
  // right, v represents the maximal end of cv, otherwise it represents its
  // minimal end.
  ind = (he->direction() == SMALLER) ? MAX_END : MIN_END;
  
  // Return the x-monotone curve.
  return (he->curve());
}

//-----------------------------------------------------------------------------
// Check whether the given infinite curve end lies on the given fictitious
// halfedge.
//
template<class Traits, class Dcel>
bool Arrangement_2<Traits,Dcel>::_is_on_fictitious_edge
    (const X_monotone_curve_2& cv, Curve_end ind,
     Boundary_type inf_x, Boundary_type inf_y,
     const DHalfedge *he,
     bool& eq_source, bool& eq_target) const
{
  eq_source = false;
  eq_target = false;

  // Get the end-vertices of the edge.
  const DVertex     *v1 = he->opposite()->vertex();
  const DVertex     *v2 = he->vertex();
  Boundary_type      he_inf;
  Comparison_result  res1, res2;
  Curve_end          v_ind;

  // Check if this is a "vertical" ficitious edge.
  if ((he_inf = v1->boundary_in_x()) != NO_BOUNDARY &&
      he_inf == v2->boundary_in_x())
  {
    // If the edge lies on x = +/- oo, the curve endpoint must also lie there.
    if (he_inf != inf_x)
      return (false);

    // Compare the y-position of the curve end to the source vertex.
    if (v1 == v_bl || v1 == v_br)
    {
      // These vertices are below any curve.
      res1 = LARGER;
    }
    else if (v1 == v_tl || v1 == v_tr)
    {
      // These vertices are above any curve.
      res1 = SMALLER;
    }
    else
    {
      const Curve_end  ind = (inf_x == MINUS_INFINITY) ? MIN_END : MAX_END;

      res1 = traits->compare_y_at_x_2_object() (cv,
                                                _get_curve (v1, v_ind),
                                                ind);
      if (res1 == EQUAL)
      {
        eq_source = true;
        return (true);
      }
    }

    // Compare the y-position of the curve end to the target vertex.
    if (v2 == v_bl || v2 == v_br)
    {
      // These vertices are below any curve.
      res2 = LARGER;
    }
    else if (v2 == v_tl || v2 == v_tr)
    {
      // These vertices are above any curve.
      res2 = SMALLER;
    }
    else
    {
      const Curve_end  ind = (inf_x == MINUS_INFINITY) ? MIN_END : MAX_END;

      res2 = traits->compare_y_at_x_2_object() (cv,
                                                _get_curve (v2, v_ind),
                                                ind);
      
      if (res2 == EQUAL)
      {
        eq_target = true;
        return (true);
      }
    }
  }
  else
  {
    // If we reched here, we have a "horizontal" fictitious halfedge.
    he_inf = v1->boundary_in_y();

    CGAL_assertion (he_inf != NO_BOUNDARY && he_inf == v2->boundary_in_y());

    // If the edge lies on y = +/- oo, the curve endpoint must also lie there
    // (and must not lies at x = +/- oo.
    if (inf_x != NO_BOUNDARY || he_inf != inf_y)
      return (false);

    // Compare the x-position of the curve end to the source vertex.
    if (v1 == v_bl || v1 == v_tl)
    {
      // These vertices are to the left of any curve.
      res1 = LARGER;
    }
    else if (v1 == v_br || v1 == v_tr)
    {
      // These vertices are to the right of any curve.
      res1 = SMALLER;
    }
    else
    {
      const X_monotone_curve_2&  cv1 = _get_curve (v1, v_ind);
      res1 = traits->compare_x_2_object() (cv, ind, cv1, v_ind);
      
      if (res1 == EQUAL)
      {
        eq_source = true;
        return (true);
      }
    }

    // Compare the x-position of the curve end to the target vertex.
    if (v2 == v_bl || v2 == v_tl)
    {
      // These vertices are to the left of any curve.
      res2 = LARGER;
    }
    else if (v2 == v_br || v2 == v_tr)
    {
      // These vertices are to the right of any curve.
      res2 = SMALLER;
    }
    else
    {
      const X_monotone_curve_2&  cv2 = _get_curve (v2, v_ind);
      res2 = traits->compare_x_2_object() (cv, ind, cv2, v_ind);

      if (res2 == EQUAL)
      {
        eq_target = true;
        return (true);
      }
    }
  }

  return (res1 != res2);
}

//-----------------------------------------------------------------------------
// Locate an ficititious halfegde on the outer CCB of a given face which
// contains the given curve end in its interior.
//
template<class Traits, class Dcel>
typename Arrangement_2<Traits,Dcel>::DHalfedge*
Arrangement_2<Traits,Dcel>::_locate_along_ccb
    (DFace *f,
     const X_monotone_curve_2& cv, Curve_end ind,
     Boundary_type inf_x, Boundary_type inf_y) const
{
  // Get a halfedge on the outer CCB of f and start traversing the CCB.
  DHalfedge         *first = f->halfedge();

  CGAL_assertion (first != NULL);

  DHalfedge         *curr = first;
  bool               eq_source, eq_target;

  do
  {
    // Note we consider only fictitious halfedges and check whether they
    // contain the relevant curve end.
    if (curr->has_null_curve() &&
        _is_on_fictitious_edge (cv, ind,
                                inf_x, inf_y,
                                curr,
                                eq_source, eq_target))
    {
      CGAL_assertion (! eq_source && ! eq_target);
      return (curr);
    }

    // Move to the next halfegde along the CCB.
    curr = curr->next();

  } while (curr != first);

  // If we reached here, we did not find a suitable halfegde.
  return (NULL);
}

//-----------------------------------------------------------------------------
// Locate the place for the given curve around the given vertex.
//
template<class Traits, class Dcel>
typename Arrangement_2<Traits,Dcel>::DHalfedge*
Arrangement_2<Traits,Dcel>::_locate_around_vertex
    (DVertex *v,
     const X_monotone_curve_2& cv,
     Curve_end ind) const
{
  // Get the first incident halfedge around v and the next halfedge.
  DHalfedge   *first = v->halfedge();
  DHalfedge   *curr = first;

  if (curr == NULL)
    return (NULL);

  DHalfedge   *next = curr->next()->opposite();

  // In case there is only one halfegde incident to v, return this halfedge.
  if (curr == next)
    return (curr);

  // Otherwise, we traverse the halfedges around v until we find the pair
  // of adjacent halfedges between which we should insert cv.
  typename Traits_adaptor_2::Is_between_cw_2  is_between_cw =
                                              traits->is_between_cw_2_object();

  bool       eq_curr, eq_next;

  while (! is_between_cw (cv, (ind == MIN_END),
                          curr->curve(), (curr->direction() == LARGER),
                          next->curve(), (next->direction() == LARGER),
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
template<class Traits, class Dcel>
unsigned int
Arrangement_2<Traits,Dcel>::_halfedge_distance (const DHalfedge *e1,
                                                const DHalfedge *e2) const
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
// Determine whether a given query halfedge lies in the interior of a new
// face we are about to create, by connecting it with another halfedge
// using a given x-monotone curve.
//
template<class Traits, class Dcel>
bool Arrangement_2<Traits,Dcel>::_is_inside_new_face
    (const DHalfedge *prev1,
     const DHalfedge *prev2,
     const X_monotone_curve_2& cv) const
{
  // Go over all halfedges of along boundary of the new face which will
  // contain prev1: As the new face is not constructed yet, this traversal
  // is simulated by going from prev2->next() to prev1 and then returning
  // over the new curve. During the traversal we locate the leftmost
  // halfedge along the boundary (i.e, the one with the lexicographically
  // smallest target vertex, which is also the first halfedge we encounter
  // if we go around this vertex in a counterclockwise direction, starting
  // from - and not including - "6 o'clock").
  // We note that this predicate is called only when closing a hole
  // comprised of finite x-monotone curves, so there are no vertices at
  // infinity or fictitious edges along the CCB.
  typename Traits_adaptor_2::Compare_xy_2  compare_xy =
                                                 traits->compare_xy_2_object();

  const DHalfedge   *left_edge = NULL;
  const DHalfedge   *first = prev2->next();
  const DHalfedge   *curr = first;
  const DHalfedge   *last = prev1->next();
  const DVertex     *v_min = prev2->vertex();
  bool               prev_is_min = true;

  do
  {
    // Note that we may visit the smallest vertex several times, for example
    // when we have:
    //
    //                          h1 /                       .
    //                           /   h2                    .
    //                        (.)-------                   .
    //                           \                         .
    //                          h3 \                       .
    //
    // If we first reach the vertex from h1, then we will encounter also h2
    // and h3, and since we take the last halfedge, we wil locate h3.
    // On the other hand, if we reach the vertex from h3, we will leave it
    // using h1's twin and will not return to it.
    if (curr->direction() == LARGER &&
        (prev_is_min ||
         curr->vertex() == v_min ||
         compare_xy (curr->vertex()->point(), v_min->point()) == SMALLER))
    {
      left_edge = curr;
      v_min = curr->vertex();
      prev_is_min = true;
    }
    else
    {
      prev_is_min = false;
    }

    curr = curr->next();
    CGAL_assertion (curr != first);

  } while (curr != last);

  // Now note that the curves of leftmost edge and its successor are defined
  // to the right of the smallest vertex. We compare them to the right of this
  // point to determine whether prev1 is inside the hole to be created or not.
  const X_monotone_curve_2  *cv_curr;
  const X_monotone_curve_2  *cv_next;

  if (left_edge != NULL)
  {
    cv_curr = &(left_edge->curve());

    // Take special care if the next curve should really be the new curve.
    if (left_edge->next() != last)
      cv_next = &(left_edge->next()->curve());
    else
      cv_next = &cv;
  }
  else
  {
    // In this case, the leftmost edge should be the one associated with the
    // new curve (which has not been created yet).
    cv_curr = &cv;
    cv_next = &(first->curve());
  }

  return (traits->compare_y_at_x_right_2_object() (*cv_curr, *cv_next,
                                                   v_min->point()) == LARGER);
}

//-----------------------------------------------------------------------------
// Compare the x-coordinates of a given vertex (which may lie at infinity) and
// the given point.
//
template<class Traits, class Dcel>
Comparison_result
Arrangement_2<Traits,Dcel>::_compare_x_imp (const Point_2& p,
                                            const DVertex* v,
                                            Tag_true) const
{
  // First check if the vertex v lies at x = -oo (then it is obviously smaller
  // than p), or at x = +oo (then it is obviously larger).
  const Boundary_type          inf_x = v->boundary_in_x();

  if (inf_x == MINUS_INFINITY)
    return (LARGER);
  else if (inf_x == PLUS_INFINITY)
    return (SMALLER);

  // Check if the vertex lies at y = +/- oo.
  const Boundary_type          inf_y = v->boundary_in_y();

  if (inf_y != NO_BOUNDARY)
  {
    // Compare the x-position of the vertical asymptote of the curve incident
    // to v with the x-coodinate of p.
    Curve_end                 ind;
    const X_monotone_curve_2& cv = _get_curve (v, ind);
    
    return (traits->compare_x_2_object() (p, cv, ind));
  }

  // In this case v represents a normal point, and we compare it with p.
  return (traits->compare_x_2_object() (p, v->point()));
}

//-----------------------------------------------------------------------------
// Compare the given vertex (which may lie at infinity) and the given point.
//
template<class Traits, class Dcel>
Comparison_result
Arrangement_2<Traits,Dcel>::_compare_xy_imp (const Point_2& p,
                                             const DVertex* v,
                                             Tag_true) const
{
  // First check if the vertex v lies at x = -oo (then it is obviously smaller
  // than p), or at x = +oo (then it is obviously larger).
  const Boundary_type          inf_x = v->boundary_in_x();

  if (inf_x == MINUS_INFINITY)
    return (LARGER);
  else if (inf_x == PLUS_INFINITY)
    return (SMALLER);

  // Check if the vertex lies at y = +/- oo.
  const Boundary_type          inf_y = v->boundary_in_y();

  if (inf_y != NO_BOUNDARY)
  {
    // Compare the x-position of the vertical asymptote of the curve incident
    // to v with the x-coodinate of p.
    Curve_end                 ind;
    const X_monotone_curve_2& cv = _get_curve (v, ind);
    Comparison_result         res =
      traits->compare_x_2_object() (p, cv, ind);

    if (res != EQUAL)
      return (res);

    // In case of equality, consider whether v lies at y = -oo or at y = +oo.
    if (inf_y == MINUS_INFINITY)
      return (LARGER);
    else
      return (SMALLER);
  }

  // In this case v represents a normal point, and we compare it with p.
  return (traits->compare_xy_2_object() (p, v->point()));
}

//-----------------------------------------------------------------------------
// Compare the relative y-position of the given point and the given edge
// (which may be fictitious if the traits class supports unbounded curves).
//
template<class Traits, class Dcel>
Comparison_result 
Arrangement_2<Traits,Dcel>::_compare_y_at_x_imp (const Point_2& p,
                                                 const DHalfedge* he,
                                                 Tag_true) const
{
  // In case of a valid edge, just compare p to its associated curve.
  if (! he->has_null_curve())
    return (traits->compare_y_at_x_2_object() (p, he->curve()));

  // Otherwise, determine on which edge of the bounding rectangle does he lie.
  // Note this can be either the top edge or the bottom edge (and not the
  // left or the right edge), as p must lie in its x-range.
  CGAL_assertion ((he->vertex()->boundary_in_x() == NO_BOUNDARY) ||
                  (he->vertex()->boundary_in_x() != 
                   he->opposite()->vertex()->boundary_in_x()));
  CGAL_assertion ((he->vertex()->boundary_in_y() != NO_BOUNDARY) &&
                  (he->vertex()->boundary_in_y() == 
                   he->opposite()->vertex()->boundary_in_y()));

  if (he->vertex()->boundary_in_y() == MINUS_INFINITY)
    // he lies on the bottom edge, so p is obviously above it.
    return (LARGER);
  else
    // he lies on the top edge, so p is obviously below it.
    return (SMALLER);
}

//-----------------------------------------------------------------------------
// Determine whether a given point lies within the region bounded by
// a boundary of a connected component.
//
template<class Traits, class Dcel>
bool Arrangement_2<Traits,Dcel>::_point_is_in (const Point_2& p,
                                               const DVertex* v,
                                               const DHalfedge *he) const
{
  CGAL_assertion (! he->is_on_hole());

  // Keep a counter of the number of x-monotone curves that intersect an upward
  // vertical emanating from p (except for some degenerate cases that are
  // explained below).
  unsigned int              n_ray_intersections = 0;

  // Go over all curves of the boundary, and count those which are above p.
  // We begin by comparing p to the source vertex of the first halfedge.
  // Note that if p coincides with this vertex, p is obviously not in the
  // interior of the component.
  const DHalfedge   *curr = he;
  Comparison_result  res_source;
  Comparison_result  res_target;
  Comparison_result  res_y_at_x;

  if (curr->opposite()->vertex() == v)
    return (false);

  res_source = _compare_xy (p, curr->opposite()->vertex());

  do
  {
    // Compare p to the target vertex of the current halfedge.
    // If the vertex v associated with p (if v is given and is not NULL)
    // on the boundary of the component, p is obviously not in the interior
    // the component.
    if (curr->vertex() == v)
      return (false);

    res_target = _compare_xy (p, curr->vertex());

    // In case the current halfedge belongs to an "antenna", namely its
    // incident face is the same as its twin's, we can simply skip it
    // (in order not to count it twice).
    if (! curr->opposite()->is_on_hole() &&
        curr->face() == curr->opposite()->face())
    {
      curr = curr->next();
      res_source = res_target;
      continue;
    }

    // Check that if we shoot a "tilted" vertical ray from p upward
    // (by "tilted" we mean the angle it forms with the x-axis is
    //  PI/2 + epsilon, where epsilon is arbitrarily small), then we hit
    // the x-monotone curve associated with curr once.
    if (res_source != res_target)
    {
      res_y_at_x = _compare_y_at_x (p, curr);

      if (res_y_at_x == SMALLER)
      {
        n_ray_intersections++;
      }        
      else if (res_y_at_x == EQUAL)
      {
        // In this case p lies on the current edge, so it is obviously not
        // contained in the interior of the component.
        return (false);
      }
    }

    // Proceed to the next halfedge along the component boundary.
    // Note that the source vertex of this halfedge is the current target.
    curr = curr->next();
    res_source = res_target;

  } while (curr != he);

  // The query point lies inside the connected components if and only if the
  // ray we shoot from it intersects the boundary an odd number of time.
  return ((n_ray_intersections % 2) != 0);
}

//-----------------------------------------------------------------------------
// Move a given hole from one face to another.
//
template<class Traits, class Dcel>
void Arrangement_2<Traits,Dcel>::_move_hole (DFace *from_face,
                                             DFace *to_face,
                                             DHoles_iter hole)
{
  // Get a halfedge lying on the boundary of the hole and the component that
  // represents the hole.
  DHalfedge   *first = *hole;
  DHole       *p_hole = first->hole();

  // Notify the observers that we are about to move a hole.
  Ccb_halfedge_circulator   circ = (Halfedge_handle(first))->ccb();

  _notify_before_move_hole (Face_handle (from_face),
                            Face_handle (to_face),
                            circ);

  // Remove the hole from the current face.
  from_face->erase_hole (hole);

  // Modify the component that represents the hole.
  p_hole->set_face (to_face);
  p_hole->set_iterator (to_face->add_hole (first));
  
  // Notify the observers that we have moved the hole.
  _notify_after_move_hole (circ);

  return;
}

//-----------------------------------------------------------------------------
// Insert the given vertex as an isolated vertex inside the given face.
//
template<class Traits, class Dcel>
void Arrangement_2<Traits,Dcel>::_insert_isolated_vertex (DFace *f,
                                                          DVertex *v)
{
  Face_handle     fh (f);
  Vertex_handle   vh (v);

  // Notify the observers that we are about to insert an isolated vertex
  // inside f.
  _notify_before_add_isolated_vertex (fh, vh);

  // Create an isolated vertex-information object,
  DIso_vert      *iv = dcel.new_isolated_vertex();

  // Set a pointer to the face containing the vertex.
  iv->set_face (f);

  // Initiate a new hole inside the given face.
  iv->set_iterator (f->add_isolated_vertex (v));
  
  // Associate the information with the vertex.
  v->set_isolated_vertex (iv);

  // Notify the observers that we have formed a new isolated vertex.
  _notify_after_add_isolated_vertex (vh);

  return;
}

//-----------------------------------------------------------------------------
// Move a given isolated vertex from one face to another.
//
template<class Traits, class Dcel>
void Arrangement_2<Traits,Dcel>::_move_isolated_vertex
    (DFace *from_face,
     DFace *to_face,
     DIsolated_vertices_iter vit)
{
  // Get the vertex the iterator points to.
  DVertex    *v = &(*vit);
  DIso_vert  *iv = v->isolated_vertex();

  // Notify the observers that we are about to move an isolated vertex.
  Vertex_handle   vh (v);

  _notify_before_move_isolated_vertex (Face_handle (from_face),
                                       Face_handle (to_face),
                                       vh);

  // Set the new face is the isolated vertex-information object.
  iv->set_face (to_face);

  // Move the isolated vertex from the first face to the other.
  from_face->erase_isolated_vertex (iv->iterator());
  iv->set_iterator (to_face->add_isolated_vertex (v));

  // Notify the observers that we have moved the isolated vertex.
  _notify_after_move_isolated_vertex (vh);

  return;
}

//-----------------------------------------------------------------------------
// Create a new vertex and associate it with the given point.
//
template<class Traits, class Dcel>
typename Arrangement_2<Traits,Dcel>::DVertex*
Arrangement_2<Traits,Dcel>::_create_vertex (const Point_2& p)
{
  // Notify the observers that we are about to create a new vertex.
  Point_2  *p_p = _new_point (p);

  _notify_before_create_vertex (*p_p);

  // Create a new vertex and associate it with the given point.
  DVertex         *v = dcel.new_vertex();

  v->set_point (p_p);

  // Notify the observers that we have just created a new vertex.
  Vertex_handle   vh (v);
  _notify_after_create_vertex (vh);

  return (v);
}

//-----------------------------------------------------------------------------
// Create a new vertex at infinity.
//
template<class Traits, class Dcel>
typename Arrangement_2<Traits,Dcel>::DVertex*
Arrangement_2<Traits,Dcel>::_create_vertex_at_infinity (Boundary_type inf_x,
                                                        Boundary_type inf_y)
{
  // Notify the observers that we are about to create a new vertex at infinity.
  _notify_before_create_vertex_at_infinity (inf_x, inf_y);

  // Create a new vertex and associate it with the given point.
  DVertex         *v = dcel.new_vertex();

  v->set_boundary (inf_x, inf_y);

  // Increment the number of vertices at infinity.
  n_inf_verts++;

  // Notify the observers that we have just created a new vertex at infinity.
  Vertex_handle   vh (v);
  _notify_after_create_vertex_at_infinity (vh);

  return (v);
}

//-----------------------------------------------------------------------------
// Insert an x-monotone curve into the arrangement, such that both its
// endpoints correspond to free arrangement vertices (newly created vertices
// or existing isolated vertices), so a new hole is formed in the face
// that contains the two vertices.
//
template<class Traits, class Dcel>
typename Arrangement_2<Traits,Dcel>::DHalfedge*
Arrangement_2<Traits,Dcel>::_insert_in_face_interior
    (const X_monotone_curve_2& cv,
     DFace *f,
     DVertex *v1, DVertex *v2,
     Comparison_result res)
{
  // Notify the observers that we are about to create a new edge.
  _notify_before_create_edge (cv, Vertex_handle (v1), Vertex_handle (v2));

  // Create a pair of twin halfedges connecting the two vertices,
  // and link them together to form a new connected component, a hole in f.
  DHalfedge           *he1 = dcel.new_edge();
  DHalfedge           *he2 = he1->opposite();
  DHole               *hole = dcel.new_hole();
  X_monotone_curve_2  *dup_cv = _new_curve (cv);

  hole->set_face (f);
  he1->set_curve (dup_cv);

  he1->set_next (he2);
  he1->set_vertex (v1);
  he1->set_hole (hole);

  he2->set_next (he1);
  he2->set_vertex (v2);
  he2->set_hole (hole);

  // Assign the incident halfedges of the two new vertices.
  v1->set_halfedge (he1);
  v2->set_halfedge (he2);

  // Set the direction of the halfedges: res indicates the direction of he2,
  // as it is the comparison result between its source (v1) and target (v2).
  he2->set_direction (res);

  // Create a handle to the new halfedge pointing at the curve target.
  Halfedge_handle   hh (he2);

  // Notify the observers that we have created a new edge.
  _notify_after_create_edge (hh);

  // Notify the observers that we are about to form a new hole inside f.
  _notify_before_add_hole (Face_handle (f), hh);

  // Initiate a new hole inside the given face.
  hole->set_iterator (f->add_hole (he2));

  // Notify the observers that we have formed a new hole.
  _notify_after_add_hole (hh->ccb());

  return (he2);
}

//-----------------------------------------------------------------------------
// Insert an x-monotone curve into the arrangement, such that one of its
// endpoints corresponds to a given arrangement vertex, given the exact
// place for the curve in the circular list around this vertex. The other
// endpoint corrsponds to a free vertex (a newly created vertex or an
// isolated vertex).
//
template<class Traits, class Dcel>
typename Arrangement_2<Traits,Dcel>::DHalfedge*
Arrangement_2<Traits,Dcel>::_insert_from_vertex (const X_monotone_curve_2& cv,
                                                 DHalfedge *prev,
                                                 DVertex *v,
                                                 Comparison_result res)
{
  // Get the incident face of the previous halfedge . Note that this will also
  // be the incident face of the two new halfedges we are about to create.
  DHole           *p_hole = (prev->is_on_hole()) ? prev->hole() : NULL;
  DFace           *p_f = (p_hole == NULL) ? prev->face() : p_hole->face();

  // The first vertex is the one that the prev halfedge points to.
  // The second vertex is given by v.
  DVertex         *v1 = prev->vertex();
  DVertex         *v2 = v;

  // Notify the observers that we are about to create a new edge.
  _notify_before_create_edge (cv, Vertex_handle (v1), Vertex_handle (v2));

  // Create a pair of twin halfedges connecting the two vertices,
  // and associate them with the given curve.
  DHalfedge           *he1 = dcel.new_edge();
  DHalfedge           *he2 = he1->opposite();
  X_monotone_curve_2  *dup_cv = _new_curve (cv);

  he1->set_curve (dup_cv);

  he1->set_vertex (v1);
  he2->set_vertex (v2);

  // Set the incident face (or hole) for the new halfedge pair.
  if (p_hole == NULL)
  {
    // On an outer component:
    he1->set_face (p_f);
    he2->set_face (p_f);
  }
  else
  {
    // On a hole:
    he1->set_hole (p_hole);
    he2->set_hole (p_hole);
  }

  // Associate the incident halfedge of the new vertex.
  v2->set_halfedge (he2);

  // Link the new halfedges around the existing vertex v1.
  he2->set_next (he1);
  he1->set_next (prev->next());

  prev->set_next (he2);

  // Set the direction of the halfedges: res indicates the direction of he2,
  // as it is the comparison result between its source and target (v).
  he2->set_direction (res);

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
template<class Traits, class Dcel>
typename Arrangement_2<Traits,Dcel>::DHalfedge*
Arrangement_2<Traits,Dcel>::_insert_at_vertices (const X_monotone_curve_2& cv,
                                                 DHalfedge *prev1,
                                                 DHalfedge *prev2,
                                                 Comparison_result res,
                                                 bool& new_face,
                                                 bool both_unbounded)
{
  // Get the vertices that match cv's endpoints.
  DVertex    *v1 = prev1->vertex();
  DVertex    *v2 = prev2->vertex();

  // Get the components containing the two previous halfedges and the incident
  // face (which should be the same for the two components).
  DHole       *hole1 = (prev1->is_on_hole()) ? prev1->hole() : NULL;
  DHole       *hole2 = (prev2->is_on_hole()) ? prev2->hole() : NULL;
  DFace       *f = (hole1 == NULL) ? prev1->face() : hole1->face();

  CGAL_precondition_code
    (DFace       *f2 = (hole2 == NULL) ? prev2->face() : hole2->face(););

  CGAL_precondition_msg
    (f == f2,
     "The two halfedges must share the same incident face.");

  // Notify the observers that we are about to create a new edge.
  _notify_before_create_edge (cv, Vertex_handle (v1), Vertex_handle (v2));

  // Create a pair of twin halfedges connecting v1 and v2 and associate them
  // with the given curve.
  DHalfedge           *he1 = dcel.new_edge();
  DHalfedge           *he2 = he1->opposite();
  X_monotone_curve_2  *dup_cv = _new_curve (cv);

  he1->set_curve (dup_cv);

  he1->set_vertex (v1);
  he2->set_vertex (v2);

  // Connect the new halfedges to the existing halfegdes around the two
  // incident vertices.
  he1->set_next (prev1->next());
  he2->set_next (prev2->next());

  prev1->set_next (he2);
  prev2->set_next (he1);

  // Set the direction of the halfedges: res indicates the direction of he2,
  // as it is the comparison result between its source and target.
  he2->set_direction (res);

  // Check whether prev1 and prev2 used to belong to the same connected
  // component or to different components.
  if (hole1 != hole2)
  {
    // In case we have connected two disconnected components (holes), no new
    // face is created. 
    new_face = false;

    // Check whether both halfedges are inner components (hole) in the same
    // face, or whether one is a hole and the other is on the outer boundary
    // of the face. 
    Face_handle       fh (f);

    if (hole1 != NULL && hole2 != NULL)
    {
      // Notify the observers that we are about to merge two holes in the face.
      _notify_before_merge_hole (fh,
                                 (Halfedge_handle(prev1))->ccb(),
                                 (Halfedge_handle(prev2))->ccb(),
                                 Halfedge_handle (he1));

      // Remove the hole prev2 belongs to, and unite it with the hole that
      // prev1 belongs to.
      f->erase_hole (hole2->iterator());

      // Set the merged component for the two new halfedges.
      he1->set_hole (hole1);
      he2->set_hole (hole1);

      // Make all halfedges along hole2 to point to hole1.
      DHalfedge       *curr;

      for (curr = he2->next(); curr != he1; curr = curr->next())
        curr->set_hole (hole1);

      // Delete the redundant hole.
      dcel.delete_hole (hole2);

      // Notify the observers that we have merged the two holes.
      _notify_after_merge_hole (fh,
                                (Halfedge_handle(he1))->ccb());
    }
    else
    {
      CGAL_assertion_msg (hole1 != NULL || hole2 != NULL,
                          "Cannot merge two outer components.");
 
      // In this case we connect a hole with the outer CCB of the face that
      // contains it. We remove the hole and associate the pair of new
      // halfedges with the outer boundary of the face f.
      DHole      *del_hole;
      DHalfedge  *ccb_first;
      DHalfedge  *ccb_last;
      
      if (hole1 != NULL)
      {
        del_hole = hole1;
        ccb_first = he1->next();
        ccb_last = he2;
      }
      else
      {
        del_hole = hole2;
        ccb_first = he2->next();
        ccb_last = he1;
      }

      he1->set_face (f);
      he2->set_face (f);

      // Notify the observers that we are about to remove a fole from the face.
      _notify_before_remove_hole (fh,
                                  (Halfedge_handle(ccb_first))->ccb());

      // Remove the hole from the face, as we have just connected this hole to
      // the outer boundary of its incident face.
      f->erase_hole (del_hole->iterator());

      // Make all halfedges along the hole to point to the outer boundary of f.
      DHalfedge       *curr;

      for (curr = ccb_first; curr != ccb_last; curr = curr->next())
        curr->set_face (f);

      // Delete the redundant hole.
      dcel.delete_hole (del_hole);

      // Notify the observers that we have removed a hole.
      _notify_after_remove_hole (fh);
    }
  }
  else
  {
    // Notify the observers that we are about to split a face.
    Face_handle       fh (f);

    _notify_before_split_face (fh, Halfedge_handle (he1));

    // The prev1 and prev2 belong to the same connected component. In this
    // case we have to create a new face, which should point to he2.
    DFace      *new_f = dcel.new_face();

    new_face = true;
    new_f->set_halfedge (he2);

    // Set the incident faces of the two new halfedges. The new face should
    // become the incident face of he2, while he1 is associated with the
    // existing face.
    he2->set_face (new_f);

    if (hole1 == NULL)
      he1->set_face (f);
    else
      he1->set_hole (hole1);

    // Set the incident face of all halfedges along he2's (new) CCB.
    DHalfedge       *curr;
    
    for (curr = he2->next(); curr != he2; curr = curr->next())
      curr->set_face (new_f);

    // Check whether the new face forms a new hole in another face.
    bool   is_hole;

    if (hole1 == NULL)
    {
      // As the outer boundary edge of the exisitng face f may be one of the
      // halfedges along the boundary of the new face, we set it to he1.
      f->set_halfedge (he1);
      is_hole = false;
    }
    else
    {
      // The prev1 and prev2 halfedges belong to a the inner boundary of their
      // component. In this case the current hole representative may not belong
      // to the hole any more, so we have to remove the hole that the used to
      // belong to the face f and replace it with the new hole, represnted
      // by he1.
      if (! (*(hole1->iterator()))->is_on_hole())
      {
        f->erase_hole (hole1->iterator());
        hole1->set_iterator (f->add_hole (he1));
      }
      is_hole = true;
    }

    // Check whether we should mark the original face and the new face as
    // bounded or as unbounded faces.
    if (both_unbounded)
    {
      // Both faces are unbounded.
      CGAL_assertion (f->is_unbounded());
      new_f->set_unbounded (true);
    }
    else
    {
      if (! f->is_unbounded())
      {
        // f is a bounded face, so the new face split from it is obviously
        // bounded.
        new_f->set_unbounded (false);
      }
      else
      {
        // We can have the two following cases: One face becomes is bounded
        // and the other is not, ot both are unbounded (but they cannot both
        // become bounded).
        std::pair<bool, bool>   unb_res = _is_face_unbounded (he1);

        if (unb_res.first || unb_res.second)
          f->set_unbounded (unb_res.first);
        new_f->set_unbounded (unb_res.second);
      }
    }

    // Notify the observers that we have split the face.
    _notify_after_split_face (fh,
                              Face_handle (new_f),
                              is_hole);
  }

  // Notify the observers that we have created a new edge.
  _notify_after_create_edge (Halfedge_handle (he2));

  // Return the halfedge directed from v1 to v2.
  return (he2);
}

//-----------------------------------------------------------------------------
// Check whether the incident face of the given halfedge and whether the
// incident face of its twin halfedge are unbounded.
//
template<class Traits, class Dcel>
std::pair<bool, bool>
Arrangement_2<Traits,Dcel>::_is_face_unbounded_imp (DHalfedge *he,
                                                    Tag_true) const
{
  int           i;
  bool          is_unbounded[2];
  DHalfedge     *first;
  DHalfedge     *curr;

  for (i = 0; i < 2; i++)
  {
    // Go over the CBB of the given halfedge and look for fictitious halfedges.
    curr = first = (i == 0) ? he : he->opposite();

    is_unbounded[i] = false;
    do
    {
      if (curr->has_null_curve())
      {
        is_unbounded[i] = true;
        break;
      }

      curr = curr->next();

    } while (curr != first);
  }

  return (std::make_pair (is_unbounded[0], is_unbounded[1]));
}

//-----------------------------------------------------------------------------
// Relocate all holes to their proper position,
// immediately after a face has split due to the insertion of a new halfedge.
//
template<class Traits, class Dcel>
void 
Arrangement_2<Traits,Dcel>::_relocate_holes_in_new_face (DHalfedge *new_he)
{
  // The given halfedge points to the new face, while its twin points to the
  // old face (the one that has just been split).
  DFace        *new_face = (new_he->is_on_hole()) ? new_he->hole()->face() :
                                                    new_he->face();
  DHalfedge    *opp_he = new_he->opposite();
  DFace        *old_face = (opp_he->is_on_hole()) ? opp_he->hole()->face() :
                                                    opp_he->face();

  CGAL_assertion (new_face != old_face);

  // Examine the holes inside the existing old face and move the relevant
  // ones into the new face.
  DHoles_iter    holes_it = old_face->holes_begin();
  DHoles_iter    hole_to_move;

  while (holes_it != old_face->holes_end())
  {  
    // In case the new edge represents the current hole in the old face
    // (note we take the opposite halfedge, as it is incident to the old face),
    // then the new face already forms a hole in the old face, and we do not
    // need to move it.
    if (*holes_it == opp_he)
    {
      ++holes_it;
      continue;
    }

    // Check whether the current hole is inside new face.
    if (_point_is_in ((*holes_it)->vertex()->point(), 
                      (*holes_it)->vertex(),
                      new_he))
    {
      // We increment the itrator before moving the hole, because this
      // operation invalidates the iterator.
      hole_to_move = holes_it;
      ++holes_it;

      // Move the hole.
      _move_hole (old_face, new_face, hole_to_move);
    }
    else
    {
      ++holes_it;
    }
  }

  return;
}

//-----------------------------------------------------------------------------
// Relocate all  vertices to their proper position,
// immediately after a face has split due to the insertion of a new halfedge.
//
template<class Traits, class Dcel>
void Arrangement_2<Traits,Dcel>::
_relocate_isolated_vertices_in_new_face (DHalfedge *new_he)
{
  // The given halfedge points to the new face, while its twin points to the
  // old face (the one that has just been split).
  DFace        *new_face = (new_he->is_on_hole()) ? new_he->hole()->face() :
                                                    new_he->face();
  DHalfedge    *opp_he = new_he->opposite();
  DFace        *old_face = (opp_he->is_on_hole()) ? opp_he->hole()->face() :
                                                    opp_he->face();

  CGAL_assertion (new_face != old_face);

  // Examine the isolated vertices inside the existing old face and move the
  // relevant ones into the new face.
  DIsolated_vertices_iter    iso_verts_it;
  DIsolated_vertices_iter    iso_vert_to_move;

  iso_verts_it = old_face->isolated_vertices_begin();
  while (iso_verts_it != old_face->isolated_vertices_end())
  {
    // Check whether the isolated vertex lies inside the new face.
    if (_point_is_in (iso_verts_it->point(), NULL, new_he))
    {
      // We increment the isolated vertices itrator before moving the vertex,
      // because this operation invalidates the iterator.
      iso_vert_to_move  = iso_verts_it;
      ++iso_verts_it;

      // Move the isolated vertex.
      _move_isolated_vertex (old_face, new_face, iso_vert_to_move);
    }
    else
    {
      ++iso_verts_it;
    }
  }

  return;
}

//-----------------------------------------------------------------------------
// Relocate all holes and isolated vertices to their proper position,
// immediately after a face has split due to the insertion of a new halfedge.
//
template<class Traits, class Dcel>
void Arrangement_2<Traits,Dcel>::_relocate_in_new_face (DHalfedge *new_he)
{
  _relocate_holes_in_new_face(new_he);
  _relocate_isolated_vertices_in_new_face(new_he);
  return;
}

//-----------------------------------------------------------------------------
// Replace the point associated with the given vertex.
//
template<class Traits, class Dcel>
void Arrangement_2<Traits,Dcel>::_modify_vertex (DVertex *v,
                                                 const Point_2& p)
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
template<class Traits, class Dcel>
void Arrangement_2<Traits,Dcel>::_modify_edge (DHalfedge *he,
                                               const X_monotone_curve_2& cv)
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
template<class Traits, class Dcel>
bool Arrangement_2<Traits,Dcel>::_are_equal (const DVertex *v,
                                             const X_monotone_curve_2& cv,
                                             Curve_end ind) const
{
  // Check if cv's end is finite and make sure that the infinity signs match
  // those of the vertex v.
  const Boundary_type    inf_x = traits->boundary_in_x_2_object() (cv, ind);
  const Boundary_type    inf_y = traits->boundary_in_y_2_object() (cv, ind);

  if (inf_x != v->boundary_in_x() || inf_y != v->boundary_in_y())
    return (false);

  if (inf_x != NO_BOUNDARY)
  {
    // The curve end lies at x = +/- oo and so does v. Make sure the curve
    // overlaps with the curve that currently induces v.
    Curve_end                  v_ind;
    const X_monotone_curve_2&  v_cv = _get_curve (v, v_ind);

    return (traits->compare_y_at_x_2_object() (cv, v_cv,
                                               v_ind) == EQUAL);
  }
  
  if (inf_y != NO_BOUNDARY)
  {
    // The curve end lies at y = +/- oo and so does v. Make sure the curve
    // overlaps with the curve that currently induces v.
    Curve_end                  v_ind;
    const X_monotone_curve_2&  v_cv = _get_curve (v, v_ind);

    return (traits->compare_x_2_object() (cv, ind,
                                          v_cv, v_ind) == EQUAL);
  }

  // If we have a valid endpoint, make sure it equals the point associated
  // with v.
  if (ind % 2 == 0)
    return (traits->equal_2_object() 
            (traits->construct_min_vertex_2_object() (cv),
             v->point()));
  else
    return (traits->equal_2_object() 
            (traits->construct_max_vertex_2_object() (cv),
             v->point()));
}

//-----------------------------------------------------------------------------
// Perform an xy-lexicographic comparison between two given vertices (which
// may lie at infinity).
//
template<class Traits, class Dcel>
Comparison_result
Arrangement_2<Traits,Dcel>::_compare_vertices_xy (const DVertex *v1,
                                                  const DVertex *v2) const
{
  // First check if one of the vertices is associated with a valid point.
  if (! v1->has_null_point())
  {
    if (! v2->has_null_point())
    {
      // If the two vertices are associated with valid points, simply compare
      // these points.
      return (traits->compare_xy_2_object() (v1->point(), v2->point()));
    }

    // Compare the valid point associated with v1 to the vertex v2.
    return (_compare_xy (v1->point(), v2));
  }
  
  if (! v2->has_null_point())
  {
    // Compare the valid point associated with v2 to the vertex v1.
    Comparison_result   res = _compare_xy (v2->point(), v1);

    if (res == EQUAL)
    {
      CGAL_assertion (v1 == v2);
      return (res);
    }

    // Swap the result.
    return ((res == LARGER) ? SMALLER : LARGER);
  }

  // In this case both vertices lie at infinity.
  const Boundary_type   inf_x1 = v1->boundary_in_x();
  const Boundary_type   inf_x2 = v2->boundary_in_x();

  if (inf_x1 != inf_x2)
  {
    // In this case, the comparison is straightforward:
    if (inf_x1 == MINUS_INFINITY || inf_x2 == PLUS_INFINITY)
      return (SMALLER);
    
    CGAL_assertion (inf_x1 == PLUS_INFINITY || inf_x2 == MINUS_INFINITY);
    return (LARGER);
  }
  else if (inf_x1 == MINUS_INFINITY)
  {
    // Both vertices lie at x = -oo, so we have to comapare their y-position.
    // We first check if one of the vertices is fictitious.
    if (v1 == v_bl || v2 == v_tl)
      return (SMALLER);
    else if (v1 == v_tl || v2 == v_bl)
      return (LARGER);

    Curve_end      ind1, ind2;

    return (traits->compare_y_at_x_2_object() (_get_curve (v1, ind1), 
                                               _get_curve (v2, ind2),
                                               MIN_END));
  }
  else if (inf_x1 == PLUS_INFINITY)
  {
    // Both vertices lie at x = +oo, so we have to comapare their y-position.
    // We first check if one of the vertices is fictitious.
    if (v1 == v_br || v2 == v_tr)
      return (SMALLER);
    else if (v1 == v_tr || v2 == v_br)
      return (LARGER);

    Curve_end      ind1, ind2;

    return (traits->compare_y_at_x_2_object() (_get_curve (v1, ind1), 
                                               _get_curve (v2, ind2),
                                               MAX_END));
  }

  // Both vertices lie at y = +/- oo and we should compare their x-coordinates.
  Curve_end                  ind1, ind2;
  const X_monotone_curve_2&  cv1 = _get_curve (v1, ind1);
  const X_monotone_curve_2&  cv2 = _get_curve (v2, ind2);
  Comparison_result          res = traits->compare_x_2_object() (cv1, ind1,
                                                                 cv2, ind2);

  if (res != EQUAL)
    return (res);

  const Boundary_type           inf_y1 = v1->boundary_in_y();
  const Boundary_type           inf_y2 = v2->boundary_in_y();

  if (inf_y1 != inf_y2)
    return ((inf_y1 == MINUS_INFINITY) ? SMALLER : LARGER);
  
  CGAL_assertion (v1 == v2);
  return (EQUAL);
}

//-----------------------------------------------------------------------------
// Split a given edge into two at a given point, and associate the given
// x-monotone curves with the split edges.
//
template<class Traits, class Dcel>
typename Arrangement_2<Traits,Dcel>::DHalfedge*
Arrangement_2<Traits,Dcel>::_split_edge (DHalfedge *e,
                                         const Point_2& p,
                                         const X_monotone_curve_2& cv1,
                                         const X_monotone_curve_2& cv2)
{
  // Allocate a new vertex and associate it with the split point.
  DVertex         *v = _create_vertex (p);

  // Split the edge from the given vertex.
  return (_split_edge (e, v, cv1, cv2));
}

//-----------------------------------------------------------------------------
// Split a given edge into two at a given vertex, and associate the given
// x-monotone curves with the split edges.
//
template<class Traits, class Dcel>
typename Arrangement_2<Traits,Dcel>::DHalfedge*
Arrangement_2<Traits,Dcel>::_split_edge (DHalfedge *e,
                                         DVertex *v,
                                         const X_monotone_curve_2& cv1,
                                         const X_monotone_curve_2& cv2)
{
  // Get the split halfedge and its twin, its source and target.
  DHalfedge       *he1 = e;
  DHalfedge       *he2 = he1->opposite();
  DHole           *hole1 = (he1->is_on_hole()) ? he1->hole() : NULL;
  DFace           *f1 = (hole1 == NULL) ? he1->face() : hole1->face();
  DHole           *hole2 = (he2->is_on_hole()) ? he2->hole() : NULL;
  DFace           *f2 = (hole2 == NULL) ? he2->face() : hole2->face();

  // Notify the observers that we are about to split an edge.
  _notify_before_split_edge (Halfedge_handle (e),
                             Vertex_handle (v),
                             cv1, cv2);

  // Allocate a pair of new halfedges.
  DHalfedge   *he3 = dcel.new_edge();
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

  if (hole1 == NULL)
    he3->set_face (f1);
  else
    he3->set_hole (hole1);

  he3->set_vertex (he1->vertex());
  he4->set_vertex (v);
  he4->set_next (he2);

  if (hole2 == NULL)
    he4->set_face (f2);
  else
    he4->set_hole (hole2);

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
// Split a given fictitious edge into two, forming a new vertex at infinity.
//
template<class Traits, class Dcel>
typename Arrangement_2<Traits,Dcel>::DHalfedge*
Arrangement_2<Traits,Dcel>::_split_fictitious_edge (DHalfedge *e,
                                                    Boundary_type inf_x,
                                                    Boundary_type inf_y)
{
  // Allocate a new vertex at infinity.
  DVertex         *v = _create_vertex_at_infinity (inf_x, inf_y);

  // Split the edge from the given vertex.
  return (_split_fictitious_edge (e, v));
}

//-----------------------------------------------------------------------------
// Split a given fictitious edge into two at a given vertex at infinity.
//
template<class Traits, class Dcel>
typename Arrangement_2<Traits,Dcel>::DHalfedge*
Arrangement_2<Traits,Dcel>::_split_fictitious_edge (DHalfedge *e,
                                                    DVertex *v)
{
  // Get the split halfedge and its twin, and their incident faces.
  // Note that he1 lies on an outer boundary of an unbounded face, while
  // its twin he2 should lie on a hole inside the fictitious face.
  DHalfedge       *he1 = e;
  DHalfedge       *he2 = he1->opposite();
  
  CGAL_assertion (! he1->is_on_hole());
  DFace           *f1 = he1->face();
  
  CGAL_assertion (f1->is_unbounded());

  CGAL_assertion (he2->is_on_hole());
  DHole           *hole2 = he2->hole();

  CGAL_assertion (hole2->face() == un_face);

  // Notify the observers that we are about to split a fictitious edge.
  _notify_before_split_fictitious_edge (Halfedge_handle (e),
                                        Vertex_handle (v));

  // Allocate a pair of new halfedges.
  DHalfedge   *he3 = dcel.new_edge();
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

  // Connect e3 between e1 and its successor.
  he3->set_next (he1->next());

  // Insert he4 between he2 and its predecessor.
  he2->prev()->set_next (he4);

  // Set the properties of the new halfedges.
  he3->set_face (f1);
  he3->set_vertex (he1->vertex());

  he4->set_vertex (v);
  he4->set_next (he2);

  he4->set_hole (hole2);

  if (he1->vertex()->halfedge() == he1)
    // If he1 is the incident halfedge to its target, he3 replaces it.
    he1->vertex()->set_halfedge (he3);

  // Update the properties of the twin halfedges we have just split.
  he1->set_next(he3);
  he1->set_vertex(v);

  // The direction of he3 is the same as he1's (and the direction of he4 is
  // the same as he2).
  he3->set_direction (he1->direction());

  // Notify the observers that we have split a fictitious edge into two.
  _notify_after_split_fictitious_edge (Halfedge_handle (he1),
                                       Halfedge_handle (he3));

  // Return a pointer to one of the existing halfedge that is incident to the
  // split vertex.
  return (he1);
}

//-----------------------------------------------------------------------------
// Remove a pair of twin halfedges from the arrangement.
// In case the removal causes the creation of a new hole, the given halfedge
// should point at this hole.
//
template<class Traits, class Dcel>
typename Arrangement_2<Traits,Dcel>::DFace*
Arrangement_2<Traits,Dcel>::_remove_edge (DHalfedge *e,
                                          bool remove_source,
                                          bool remove_target)
{
  // Get the pair of twin edges to be removed, the connected components they
  // belong to and their incident faces.
  DHalfedge   *he1 = e;
  DHalfedge   *he2 = e->opposite();
  DHole       *hole1 = (he1->is_on_hole()) ? he1->hole() : NULL;
  DFace       *f1 = (hole1 == NULL) ? he1->face() : hole1->face();
  DHole       *hole2 = (he2->is_on_hole()) ? he2->hole() : NULL;
  DFace       *f2 = (hole2 == NULL) ? he2->face() : hole2->face();
  DHalfedge   *prev1 = NULL;
  DHalfedge   *prev2 = NULL;

  // Notify the observers that we are about to remove an edge.
  Halfedge_handle  hh (e);

  _notify_before_remove_edge (hh);

  // Check if the two incident faces are equal, in which case no face will be
  // merged and deleted (and a hole may be created).
  if (f1 == f2)
  {
    CGAL_assertion (hole1 == hole2);
    
    // Check if the two halfedges are successors along the face boundary.
    if (he1->next() == he2 && he2->next() == he1)
    {
      // The two halfedges form a "singleton" hole inside the incident face -
      // remove it. First notify the observers that we are about to remove
      // this hole.
      Face_handle               fh (f1);
      Ccb_halfedge_circulator   hole = (Halfedge_handle(he1))->ccb();
      
      _notify_before_remove_hole (fh,
                                  hole);

      // Erase the hole from the incident face and delete the corresponding
      // component.      
      f1->erase_hole (hole1->iterator());

      dcel.delete_hole (hole1);

      _notify_after_remove_hole (fh);

      // Remove the end-vertices, if necessary.
      if (remove_target)
      {
        // Delete the he1's target vertex and its associated point.
        _notify_before_remove_vertex (Vertex_handle (he1->vertex()));

        _delete_point (he1->vertex()->point());
        dcel.delete_vertex (he1->vertex());

        _notify_after_remove_vertex ();
      }
      else
      {
        // The remaining target vertex now becomes an isolated vertex inside
        // the containing face:
        _insert_isolated_vertex (f1, he1->vertex());
      }

      if (remove_source)
      {
        // Delete the he1's source vertex and its associated point.
        _notify_before_remove_vertex (Vertex_handle (he2->vertex()));

        _delete_point (he2->vertex()->point());
        dcel.delete_vertex (he2->vertex());

        _notify_after_remove_vertex ();
      }
      else
      {
        // The remaining source vertex now becomes an isolated vertex inside
        // the containing face:
        _insert_isolated_vertex (f1, he2->vertex());
      }

      // Delete the curve associated with the edge to be removed.
      _delete_curve (he1->curve());
      dcel.delete_edge (he1);

      // Notify the observers that an edge has been deleted.
      _notify_after_remove_edge();

      // Return the face that used to contain the hole.
      return (f1);
    }
    else if (he1->next() == he2 || he2->next() == he1)
    {
      // In this case the two halfedges form an "antenna".
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

      // In case the halfedges to be deleted are represantatives of a face
      // boundary or a hole inside the face, replace them by prev1.
      // Note that as we just change the hole representatives, we do not have
      // to notify the observers.
      if (hole1 == NULL && (f1->halfedge() == he1 || f1->halfedge() == he2))
      {
        f1->set_halfedge (prev1);
      }
      else if (hole1 != NULL)
      {
        if (*(hole1->iterator()) == he1 || *(hole1->iterator()) == he2)
        {
          f1->erase_hole (hole1->iterator());
          hole1->set_iterator (f1->add_hole (prev1));
        }
      }

      // In case he2 is the representative halfegde of its target vertex,
      // replace it by prev1 (which also points at this vertex).
      if (he2->vertex()->halfedge() == he2)
          he2->vertex()->set_halfedge (prev1);

      // Remove the base vertex, in case it lies at infinity.
      if (he2->vertex()->has_null_point())
        _remove_vertex_at_infinity (he2->vertex());

      // Remove the redundant tip vertex, if necessary.
      if (remove_tip_vertex)
      {
        // Delete the vertex that forms the tip of the "antenna".
        _notify_before_remove_vertex (Vertex_handle (he1->vertex()));

        _delete_point (he1->vertex()->point());
        dcel.delete_vertex (he1->vertex());

        _notify_after_remove_vertex();
      }
      else
      {
        // The remaining "antenna" tip now becomes an isolated vertex inside
        // the containing face:
        _insert_isolated_vertex (f1, he1->vertex());
      }

      // Delete the curve associated with the edge to be removed.
      _delete_curve (he1->curve());
      dcel.delete_edge (he1);

      // Notify the observers that an edge has been deleted.
      _notify_after_remove_edge();

      // Return the incident face.
      return (f1);
    }

    // In this case we have to create a new hole.
    prev1 = he1->prev();
    prev2 = he2->prev();

    // Determine whether the he1's predecessor lies of the outer or the inner
    // boundary of the incident face (i.e. whether it represents a hole).
    if (hole1 == NULL)
    {      
      // We have to create a new hole in the interior of the incident face.
      //
      //    +-----------------------------+
      //    | prev1                       |
      //    |         +----+              |
      //    +.........+    |              |
      //    |         +----+              |
      //    |                             |
      //    +-----------------------------+
      //
      // Note that it is guaranteed that he1 points at this new hole, while
      // he2 points at the boundary of the face that contains this hole.
      // First notify the observers we are about to form a new hole inside f1.
      _notify_before_add_hole (Face_handle (f1),
                               Halfedge_handle (he1->next()));

      // Create a new component that represents the new hole.
      DHole   *new_hole = dcel.new_hole();

      new_hole->set_face (f1);

      // Create the new hole.
      new_hole->set_iterator (f1->add_hole (he1->next()));

      // Associate all halfedges along the hole boundary with the new hole
      // component.
      DHalfedge       *curr;

      for (curr = he1->next(); curr != he2; curr = curr->next())
        curr->set_hole (new_hole);

      // As the outer boundary of f1 may be represented by any of the
      // halfedges in between he1 -> ... -> he2 (the halfedges in between
      // represent the outer boundary of the new hole that is formed),
      // We represent the outer boundary of f1 by prev1, which definately
      // stays on the outer boundary.
      f1->set_halfedge (prev1);

      // Notify the observers that a new hole has been formed.
      Ccb_halfedge_circulator   hccb = (Halfedge_handle(he1->next()))->ccb();

      _notify_after_add_hole (hccb);
    }
    else
    {
      // We have to split an existing hole into two holes.
      //
      //    +-----------------------------+
      //    |           prev1             |
      //    |   +----+ /    +----+        |
      //    |   |    +......+    |        |
      //    |   +----+      +----+        |
      //    |                             |
      //    +-----------------------------+
      //
      // Notify the observers we are about to split a hole.
      _notify_before_split_hole (Face_handle (f1),
                                 (Halfedge_handle (*(hole1->iterator())))->
                                                                         ccb(),
                                 Halfedge_handle (he1));

      // We first remove the existing hole from the incident face, and replace
      // it with a new hole that prev1 is its representative halfedge.
      f1->erase_hole (hole1->iterator());
      hole1->set_iterator (f1->add_hole (prev1));

      // Create a new component that represents the new hole we split.
      DHole     *new_hole = dcel.new_hole();

      new_hole->set_face (f1);

      // Insert the new hole.
      new_hole->set_iterator (f1->add_hole (prev2));

      // Associate all halfedges along the hole boundary with the new hole
      // component.
      DHalfedge       *curr;

      for (curr = he1->next(); curr != he2; curr = curr->next())
        curr->set_hole (new_hole);

      // Notify the observers that the hole has been split.
      _notify_after_split_hole (Face_handle (f1),
                                (Halfedge_handle (prev1))->ccb(),
                                (Halfedge_handle (prev2))->ccb());
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

    // Remove the end vertices, in case they lie at infinity.
    if (he1->vertex()->has_null_point())
      _remove_vertex_at_infinity (he1->vertex());

    if (he2->vertex()->has_null_point())
      _remove_vertex_at_infinity (he2->vertex());

    // Delete the curve associated with the edge to be removed.
    _delete_curve (he1->curve());

    // Delete the pair of halfedges.
    dcel.delete_edge (he1);

    // Notify the observers that an edge has been deleted.
    _notify_after_remove_edge();

    // Return the incident face.
    return (f1);
  }

  // The two incident faces are not the same - in this case, the edge we are
  // about to delete separates these two faces. We therefore have to delete
  // one of these faces and merge it with the other face.
  // First notify the observers we are about to merge the two faces.
  _notify_before_merge_face (Face_handle (f1),
                             Face_handle (f2),
                             Halfedge_handle (he1));

  // We begin by checking whether one of the face is a hole inside the other
  // face.
  DHalfedge   *ccb;

  prev1 = he1->prev();
  prev2 = he2->prev();

  CGAL_assertion (hole1 == NULL || hole2 == NULL);

  if (hole1 == NULL && hole2 == NULL)
  {
    // Both halfedges lie on the outer boundary of their incident faces.
    // We first set the connected component of f2's outer-boundary halfedges
    // to be the same as f1's outer component.
    for (ccb = he2->next(); ccb != he2; ccb = ccb->next())
      ccb->set_face (f1);

    // Move the holes inside f2 to f1.
    DHoles_iter    holes_it = f2->holes_begin();
    DHoles_iter    hole_to_move;

    while (holes_it != f2->holes_end())
    {
      // We increment the holes itrator before moving the hole, because
      // this operation invalidates the iterator.
      hole_to_move  = holes_it;
      ++holes_it;

      _move_hole (f2, f1, hole_to_move);
    }

    // In case he1, which is about to be deleted, is a representative
    // halfedge of the face f1, we replace it by its predecessor.
    if (f1->halfedge() == he1)
      f1->set_halfedge (prev1);

    // Move the isolated vertices inside f2 to f1.
    DIsolated_vertices_iter    iso_verts_it;
    DIsolated_vertices_iter    iso_vert_to_move;

    iso_verts_it = f2->isolated_vertices_begin();
    while (iso_verts_it != f2->isolated_vertices_end())
    {
      // We increment the isolated vertices itrator before moving the vertex,
      // because this operation invalidates the iterator.
      iso_vert_to_move  = iso_verts_it;
      ++iso_verts_it;

      _move_isolated_vertex (f2, f1, iso_vert_to_move); 
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

    // If the face f2 we have just merged with f1 is unbounded, then the merged
    // face is also unbounded.
    if (f2->is_unbounded())
      f1->set_unbounded(true);

    // Delete the face f2.
    dcel.delete_face (f2);
      
    // Notify the observers that the faces have been merged.
    _notify_after_merge_face (Face_handle (f1));

    // Remove the end vertices, in case they lie at infinity.
    if (he1->vertex()->has_null_point())
      _remove_vertex_at_infinity (he1->vertex());

    if (he2->vertex()->has_null_point())
      _remove_vertex_at_infinity (he2->vertex());

     // Delete the pair of halfedges.
    dcel.delete_edge (he1);

    // Notify the observers that an edge has been deleted.
    _notify_after_remove_edge();

    // Return the merged face.
    return (f1); 
  }

  // In this case we merge a face with another face that now forms a hole
  // inside it. In this case we make sure that f1 contains the hole f2, so
  // we can merge f2 with it (we swap roles between the halfedges if
  // necessary).
  if (hole2 != NULL)
  { 
    he1 = he2;
    he2 = he1->opposite();

    hole1 = hole2;
    hole2 = NULL;

    DFace   *tf = f1;
    f1 = f2;
    f2 = tf;
    
    prev1 = he1->prev();
    prev2 = he2->prev();
  }

  // By removing the edge we open a closed face f2 contained in f1. By doing
  // this, the outer boundary of f2 unites with the hole boundary that hole1
  // represents. We therefore have to set component of all halfedges along the
  // boundary of f2 to be hole1.
  for (ccb = he2->next(); ccb != he2; ccb = ccb->next())
    ccb->set_hole (hole1);

  // Move the holes inside f2 to f1.
  DHoles_iter    holes_it = f2->holes_begin();
  DHoles_iter    hole_to_move;

  while (holes_it != f2->holes_end())
  {
    // We increment the holes itrator before moving the hole, because
    // this operation invalidates the iterator.
    hole_to_move  = holes_it;
    ++holes_it;

    _move_hole (f2, f1, hole_to_move); 
  }

  // Move the isolated vertices inside f2 to f1.
  DIsolated_vertices_iter    iso_verts_it;
  DIsolated_vertices_iter    iso_vert_to_move;

  iso_verts_it = f2->isolated_vertices_begin();
  while (iso_verts_it != f2->isolated_vertices_end())
  {
    // We increment the isolated vertices itrator before moving the vertex,
    // because this operation invalidates the iterator.
    iso_vert_to_move  = iso_verts_it;
    ++iso_verts_it;

    _move_isolated_vertex (f2, f1, iso_vert_to_move); 
  }
        
  // Notice that f2 will be merged with f1, but its boundary will still be
  // a hole inside this face. In case he1 is a represantative of this hole,
  // replace it by its predecessor.
  // Note that we do not have to notify the observers about this operation.
  if (*(hole1->iterator()) == he1)
  {
    f1->erase_hole (hole1->iterator());
    hole1->set_iterator (f1->add_hole (prev1));
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

  // If the face f2 we have just merged with f1 is unbounded, then the merged
  // face is also unbounded.
  if (f2->is_unbounded())
    f1->set_unbounded(true);

  // Delete the face f2.
  dcel.delete_face (f2);
      
  // Notify the observers that the faces have been merged.
  _notify_after_merge_face (Face_handle (f1));

  // Remove the end vertices, in case they lie at infinity.
  if (he1->vertex()->has_null_point())
    _remove_vertex_at_infinity (he1->vertex());
  
  if (he2->vertex()->has_null_point())
    _remove_vertex_at_infinity (he2->vertex());

  // Delete the pair of halfedges.
  dcel.delete_edge (he1);

  // Notify the observers that an edge has been deleted.
  _notify_after_remove_edge();

  // Return the merged face.
  return (f1);
}

//-----------------------------------------------------------------------------
// Remove an isolated vertex from the interior of a given face (but not from
// the dcel)
//
template<class Traits, class Dcel>
void Arrangement_2<Traits,Dcel>::_remove_isolated_vertex (DVertex* v)
{
  DIso_vert  *iv = v->isolated_vertex();
  iv->face()->erase_isolated_vertex (iv->iterator());
  dcel.delete_isolated_vertex (iv);
}

//-----------------------------------------------------------------------------
// Remove a vertex at infinity, causing its two incident fictitious edges
// to merge.
template<class Traits, class Dcel>
void Arrangement_2<Traits,Dcel>::_remove_vertex_at_infinity (DVertex* v)
{
  // Assign pointers to the halfedges incident to v (make sure that there
  // are exactly teo pairs of fictitious halfedges), such that we have:
  //
  //            he1      he3
  //         -------> ------->
  //       (.)      (.)v     (.)
  //         <------- <-------
  //            he2      he4
  //
  DHalfedge   *he1 = v->halfedge();
  DHalfedge   *he2 = he1->opposite();
  DHalfedge   *he3 = he1->next();
  DHalfedge   *he4 = he3->opposite();

  CGAL_assertion (he1->has_null_curve() && he3->has_null_curve() &&
                  he4->next() == he2);

  // Keep pointers to the components that contain two halfedges he3 and he2,
  // pointing at the end vertices of the merged halfedge.
  DHole       *hole1 = (he3->is_on_hole()) ? he3->hole() : NULL;
  DFace       *f1 = (hole1 == NULL) ? he3->face() : hole1->face();
  DHole       *hole2 = (he4->is_on_hole()) ? he4->hole() : NULL;
  DFace       *f2 = (hole2 == NULL) ? he4->face() : hole2->face();

  // Notify the observers that we are about to merge a fictitious edge.
  _notify_before_merge_fictitious_edge (Halfedge_handle (he1),
                                        Halfedge_handle (he3));

  // As he1 and he2 will evetually represent the merged edge, while he3 and he4
  // will be deleted, check if the deleted halfedges are represantatives of a
  // face boundary or a hole inside these faces. If so, replace he3 by he1 and
  // he4 by he2. Note that as we just change the hole representatives, we do
  // not have to notify the observers about the change.
  if (hole1 == NULL && f1->halfedge() == he3)
  {
    f1->set_halfedge (he1);
  }
  else if (hole1 != NULL)
  {
    if (*(hole1->iterator()) == he3)
    {
      f1->erase_hole (hole1->iterator());
      hole1->set_iterator (f1->add_hole (he1));
    }
  }

  if (hole2 == NULL && f2->halfedge() == he4)
  {
    f2->set_halfedge (he2);
  }
  else if (hole2 != NULL)
  {
    if (*(hole2->iterator()) == he4)
    {
      f2->erase_hole (hole2->iterator());
      hole2->set_iterator (f2->add_hole (he2));
    }
  }

  if (he3->vertex()->halfedge() == he3)
    // If he3 is the incident halfedge to its target, replace it by he1.
    he3->vertex()->set_halfedge (he1);

  // Disconnect he3 and he4 from the edge list.
  CGAL_assertion (he3->next() != he4);

  he1->set_next (he3->next());
  he4->prev()->set_next (he2);

  // Set the properties of the merged edge.
  he1->set_vertex (he3->vertex());

  // Notify the observers that we are about to delete a vertex.
  _notify_before_remove_vertex_at_infinity (Vertex_handle (v));

  // Delete the merged vertex.
  dcel.delete_vertex (v);

  // Notify the observers that the vertex has been deleted.
  _notify_after_remove_vertex_at_infinity ();

  // Delete the redundant halfedge pair.
  dcel.delete_edge (he3);

  //decrement the counter of vertices at infinity.
  --n_inf_verts;

  // Notify the observers that the edge has been merge.
  _notify_after_merge_fictitious_edge (Halfedge_handle (he1));

  return;
}

//---------------------------------------------------------------------------
// Check whether the arrangement is valid. In particular, check the
// validity of each vertex, halfedge, and face.
//
template<class Traits, class Dcel>
bool Arrangement_2<Traits,Dcel>::is_valid() const
{
  Vertex_const_iterator   vit;

  for (vit = vertices_begin(); vit != vertices_end(); ++vit)
  {
    bool is_vertex_valid = _is_valid (vit);
    if (!is_vertex_valid) 
    {
      CGAL_warning_msg(is_vertex_valid, "Invalid vertex.");
      return (false);
    }
  }
    
  Halfedge_const_iterator heit;

  for (heit = halfedges_begin(); heit != halfedges_end(); ++heit) 
  {
    bool is_halfedge_valid =  _is_valid (heit);
    if (! is_halfedge_valid) 
    {
      CGAL_warning_msg(is_halfedge_valid, "Invalid halfedge.");
      return (false);
    }
  }
    
  Face_const_iterator     fit;

  for (fit = faces_begin(); fit != faces_end(); ++fit) 
  {
    bool is_face_valid = _is_valid (fit);
    if (! is_face_valid)
    {
      CGAL_warning_msg(is_face_valid, "Invalid face.");
      return (false);
    }
  }

  bool  are_vertices_unique = _are_vertices_unique();
  if (!are_vertices_unique)
  {
    CGAL_warning_msg(are_vertices_unique,
                 "Found two vertices with the same geometric point.");
    return (false);
  }

  // If we reached here, the arrangement is valid.
  return (true);
}

//---------------------------------------------------------------------------
// Check the validity of a vertex.
//
template<class Traits, class Dcel>
bool Arrangement_2<Traits,Dcel>::_is_valid(Vertex_const_handle v) const
{    
  // In case of an isolated vertex, make sure it is really contained in the
  // face attached to it.
  if (v->is_isolated())
  {
    return (true);
  }

  // Make sure that the vertex is the target of all its incident halfedges.
  Halfedge_around_vertex_const_circulator circ = v->incident_halfedges();
  Halfedge_around_vertex_const_circulator start = circ;

  do
  {
    if (circ->target() != v)
      return (false);
    
    ++circ;
  } while (circ != start);

  // Make sure the curves are correctly ordered around the vertex.
  if (!_are_curves_ordered_cw_around_vertrex (v))
    return (false);
  
  return (true);
}


//---------------------------------------------------------------------------
// Check the validity of a halfedge.
//
template<class Traits, class Dcel>
bool Arrangement_2<Traits,Dcel>::_is_valid (Halfedge_const_handle he) const
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

  // Check that the end points of the curve associated with the halfedge
  // really equal the source and target vertices of this halfedge.
  const X_monotone_curve_2&  cv = he->curve();
  const DVertex             *source = _vertex (he->source());
  const DVertex             *target = _vertex (he->target());
  const Comparison_result    res = _compare_vertices_xy (source, target);

  if (he->direction() != res)
    return (false);

  if (res == SMALLER)
  {
    return (_are_equal (_vertex (he->source()), cv, MIN_END) &&
            _are_equal (_vertex (he->target()), cv, MAX_END));
  }
  else if (res == LARGER)
  {
    return (_are_equal (_vertex (he->source()), cv, MAX_END) &&
            _are_equal (_vertex (he->target()), cv, MIN_END));
  }

  // In that case, the source and target of the halfedge are equal.
  return (false);
}

//---------------------------------------------------------------------------
// Check the validity of a face.
//
template<class Traits, class Dcel>
bool Arrangement_2<Traits,Dcel>::_is_valid (Face_const_handle f) const
{   
  // check if all edges of f (on all ccb's) refer to f (as their face)
  Inner_ccb_const_iterator           iccbit;
  Ccb_halfedge_const_circulator ccb_circ;
   
  if (! f->is_unbounded())
  {
    ccb_circ = f->outer_ccb();
    if (! _is_valid (ccb_circ, f))
      return (false);
  }
    
  for (iccbit = f->holes_begin(); iccbit != f->holes_end(); ++iccbit)      
  {
    ccb_circ = *iccbit;
    if (! _is_valid (ccb_circ, f))
      return (false);
  }
   
  return (true);
}

//---------------------------------------------------------------------------
// Check the validity of a CCB of a given face.
//
template<class Traits, class Dcel>
bool
Arrangement_2<Traits,Dcel>::_is_valid (Ccb_halfedge_const_circulator start,
                                       Face_const_handle f) const
{
  // Make sure that all halfegdes along the CCB have the same incident face.
  Ccb_halfedge_const_circulator circ = start;
  
  do 
  {
    if (circ->face() != f)
      return (false);

    ++circ;
  }while (circ != start);
   
  return (true);
}

//---------------------------------------------------------------------------
// Check that all vertices are unique (no two vertices with the same 
// geometric point).
//
template<class Traits, class Dcel>
bool Arrangement_2<Traits,Dcel>::_are_vertices_unique() const
{
  if (number_of_vertices() < 2)
    return (true);

  // Store all points.
  std::vector<Point_2>  points_vec (number_of_vertices());
  Vertex_const_iterator vit;
  unsigned int          i = 0;
  
  for (vit = vertices_begin(); vit != vertices_end(); ++vit, ++i)
    points_vec[i] = vit->point();

  // Sort the vector of points and make sure no two adjacent points in the
  // sorted vector are equal.
  typedef typename Traits_adaptor_2::Compare_xy_2    Compare_xy_2; 
  
  typename Traits_adaptor_2::Equal_2  equal = traits->equal_2_object();
  Compare_xy_2                  compare_xy = traits->compare_xy_2_object();
  Compare_to_less<Compare_xy_2> cmp = compare_to_less (compare_xy);
  
  std::sort(points_vec.begin(), points_vec.end(), cmp);
  for (i = 1; i < points_vec.size(); ++i)
  {
    if(equal(points_vec[i-1], points_vec[i]))
      return (false);
  }
  
  return (true);
}

//---------------------------------------------------------------------------
// Check that the curves around a given vertex are ordered clockwise .
//
template<class Traits, class Dcel>
bool Arrangement_2<Traits,Dcel>::_are_curves_ordered_cw_around_vertrex
    (Vertex_const_handle v) const
{
  if(v->degree() < 3)
    return (true);

  typename Traits_adaptor_2::Is_between_cw_2  is_between_cw =
                                              traits->is_between_cw_2_object();

  Halfedge_around_vertex_const_circulator circ = v->incident_halfedges();
  Halfedge_around_vertex_const_circulator start = circ;
  Halfedge_around_vertex_const_circulator prev, next;
  bool                                    eq1, eq2;

  do
  {
    prev = circ; --prev;
    next = circ; ++next;

    if (!is_between_cw (circ->curve(), (circ->direction() == LARGER),
                        prev->curve(), (prev->direction() == LARGER),
                        next->curve(), (next->direction() == LARGER),
                        v->point(), eq1, eq2))
      return (false);

    ++circ;
   } while (circ != start);
   
  return (true);
}

CGAL_END_NAMESPACE

#endif
