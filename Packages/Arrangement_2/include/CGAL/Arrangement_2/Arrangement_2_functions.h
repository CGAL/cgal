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
// Author(s)     : Ron Wein          <wein@post.tau.ac.il>
//                 (based on old version by: Iddo Hanniel,
//                                           Eyal Flato,
//                                           Oren Nechushtan,
//                                           Ester Ezra,
//                                           Shai Hirsch,
//                                           and Eugene Lipovetsky)
#ifndef CGAL_ARRANGEMENT_2_FUNCTIONS_H
#define CGAL_ARRANGEMENT_2_FUNCTIONS_H

/*! \file
 * Member-function definitions for the Arrangement_2<Traits,Dcel> class.
 */

CGAL_BEGIN_NAMESPACE

//-----------------------------------------------------------------------------
// Default constructor.
//
template<class Traits, class Dcel>
Arrangement_2<Traits,Dcel>::Arrangement_2 ()
{
  // Set an empty unbounded face.
  un_face = dcel.new_face();
  un_face->set_halfedge(NULL);

  // Allocate the traits.
  traits = new Traits_wrapper_2;
  own_traits = true;
}

//-----------------------------------------------------------------------------
// Copy constructor.
//
template<class Traits, class Dcel>
Arrangement_2<Traits,Dcel>::Arrangement_2 (const Self& arr) :
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
  // Set an empty unbounded face.
  un_face = dcel.new_face();
  un_face->set_halfedge(NULL);

  // Set the traits.
  traits = static_cast<Traits_wrapper_2*>(tr);
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

  // Go over the vertices and create duplicates of the stored points.
  typename Dcel::Vertex_iterator       vit;
  Stored_point_2                      *dup_p;
  DVertex                             *p_v;

  for (vit = dcel.vertices_begin(); vit != dcel.vertices_end(); ++vit)
  {
    // Create the duplicate point and store it in the points container.
    p_v = &(*vit);
    dup_p = _new_point (p_v->point());

    // Associate the vertex with the duplicated point.
    p_v->set_point (dup_p);
  }

  // Go over the edge and create duplicates of the stored curves.
  typename Dcel::Edge_iterator         eit;
  Stored_curve_2                      *dup_cv;
  DHalfedge                           *p_e;

  for (eit = dcel.edges_begin(); eit != dcel.edges_end(); ++eit)
  {
    // Create the duplicate curve and store it in the curves container.
    p_e = &(*eit);
    dup_cv = _new_curve (p_e->curve());

    // Associate the halfedge (and its twin) with the duplicated curve. 
    p_e->set_curve (dup_cv);
  }

  // Take care of the traits object.
  if (own_traits && traits != NULL)
    delete traits;

  if (arr.own_traits)
    traits = new Traits_wrapper_2;
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
  // Clear the DCEL.
  dcel.delete_all();

  // Free all stored points.
  typename Points_container::iterator    pit = points.begin(), p_curr;

  while (pit != points.end())
  {
    p_curr = pit;
    ++pit;
    _delete_point (*p_curr);
  }
  points.destroy();

  // Free all stores curves.
  typename X_curves_container::iterator  cvit = curves.begin(), cv_curr;

  while (cvit != curves.end())
  {
    cv_curr = cvit;
    ++cvit;
    _delete_curve (*cv_curr);
  }
  curves.destroy();

  // Free the traits object, if necessary.
  if (own_traits)
    delete traits;
}

//-----------------------------------------------------------------------------
// Clear the arrangement.
template<class Traits, class Dcel>
void Arrangement_2<Traits,Dcel>::clear()
{
  // Notify the observers that we are about to clear the arragement.
  _notify_before_clear ();

  // Clear the DCEL so it contains just an empty unbounded face.
  dcel.delete_all();
  un_face = dcel.new_face();
  un_face->set_halfedge (NULL);

  // Clear the point and curve containers.
  points.destroy();
  curves.destroy();

  // Notify the observers that we have just cleared the arragement.
  _notify_after_clear (unbounded_face());

  return;
}

//-----------------------------------------------------------------------------
// Replace the point associated with the given vertex.
//
template<class Traits, class Dcel>
typename Arrangement_2<Traits,Dcel>::Vertex_handle
Arrangement_2<Traits,Dcel>::modify_vertex (Vertex_handle vh,
                                           const Point_2& p)
{
  CGAL_precondition_msg (traits->equal_2_object() (vh->point(), p),
                         "The new point is different from the current one.");

  // Modify the vertex.
  _modify_vertex (_vertex (vh), p);

  // Return a handle to the modified vertex.
  return (vh);
}

//-----------------------------------------------------------------------------
// Insert an isolated vertex in the interior of a given face.
//
template<class Traits, class Dcel>
typename Arrangement_2<Traits,Dcel>::Vertex_handle
Arrangement_2<Traits,Dcel>::insert_isolated_vertex (const Point_2& p,
                                                    Face_handle f)
{
  DFace          *p_f = _face (f);

  // Create a new vertex associated with the given point.
  Stored_point_2 *p_p = _new_point (p);

  _notify_before_create_vertex (*p_p);

  DVertex         *v = dcel.new_vertex();

  v->set_point (p_p);

  Vertex_handle   vh (v);
  _notify_after_create_vertex (vh);

  // Notify the observers that we are about to insert an isolated vertex
  // inside f.
  _notify_before_add_isolated_vertex (f, vh);

  // Associate the incident haldefge of the vertex with a halfedge on the
  // outer boundary of the face (or NULL if the face is unbounded).
  v->set_halfedge (p_f->halfedge());

  // Initiate a new hole inside the given face.
  p_f->add_isolated_vertex (v);

  // Notify the observers that we have formed a new hole.
  _notify_after_add_isolated_vertex (vh);

  // Return a handle to the new isolated vertex.
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
  DHalfedge    *p_he = p_v->halfedge();
  DFace        *p_f = (p_he != NULL) ?  p_he->face() : un_face;
  Face_handle   f = Face_handle (p_f);
  
  // Notify the observers that we are abount to remove a vertex.
  _notify_before_remove_vertex (v);

  // Locate the vertex in the isolated vertices container of the face and
  // remove it.
  DIsolated_vertices_iter     iso_verts_it;

  for (iso_verts_it = p_f->isolated_vertices_begin();
       iso_verts_it != p_f->isolated_vertices_end(); ++iso_verts_it)
  {
    if (&(*iso_verts_it) == p_v)
    {
      // We have found the isolated vertex - erase it.
      p_f->erase_isolated_vertex (iso_verts_it);
      break;
    }
  }

  CGAL_assertion_msg (iso_verts_it != p_f->isolated_vertices_end(),
                      "Failed to locate the isolated vertex inside a face.");

  // Delete the vertex.
  _delete_point (p_v->point());
  dcel.delete_vertex (p_v);

  // Notify the observers that the vertex has been removed.
  _notify_after_remove_vertex ();

  // Return a handle for the face that used to contain the deleted vertex.
  return (f);
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
  DFace           *p_f = _face (f);

  // Create a new vertex associated with the curve's left endpoints.
  // We also notify the observers on the creation of this vertex.
  Stored_point_2  *p1 =
    _new_point (traits->construct_min_vertex_2_object()(cv));

  _notify_before_create_vertex (*p1);

  DVertex         *v1 = dcel.new_vertex();

  v1->set_point (p1);

  _notify_after_create_vertex (Vertex_handle (v1));

  // Create a new vertex associated with the curve's right endpoint.
  // We also notify the observers on the creation of this vertex.
  Stored_point_2  *p2 =
    _new_point (traits->construct_max_vertex_2_object()(cv));

  _notify_before_create_vertex (*p2);

  DVertex         *v2 = dcel.new_vertex();

  v2->set_point (p2);

  _notify_after_create_vertex (Vertex_handle (v2));

  // Create the edge connecting the two vertices.
  DHalfedge       *new_he = _insert_in_face_interior (cv, p_f, v1, v2);

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
    (traits->equal_2_object() (v->point(), 
                               traits->construct_min_vertex_2_object()(cv)),
     "The input vertex should be the left curve endpoint.");

  // Create a new vertex associated with the curve's right endpoint.
  // We also notify the observers on the creation of this vertex.
  Stored_point_2  *p2 =
    _new_point (traits->construct_max_vertex_2_object()(cv));

  _notify_before_create_vertex (*p2);

  DVertex         *v2 = dcel.new_vertex();

  v2->set_point (p2);

  _notify_after_create_vertex (Vertex_handle (v2));

  // Check if the given vertex is isolated.
  if (v->is_isolated())
  {
    // The given vertex is an isolated one: We should in fact insert the curve
    // in the interior of the face containing this vertex.
    DVertex    *v1 = _vertex (v);
    DHalfedge  *p_he = v1->halfedge();
    DFace      *p_f = (p_he != NULL) ?  p_he->face() : un_face;

    // v1 will not be isolated any more - remove it from the isolated vertices
    // container of its containing face. 
    _find_and_erase_isolated_vertex (p_f, v1);

    // Create the edge connecting the two vertices.
    DHalfedge  *new_he = _insert_in_face_interior (cv, p_f, v1, v2);

    // Return a handle to the new halfedge directed from left to right.
    return (Halfedge_handle (new_he));
  }

  // Go over the incident halfedges around v and find the halfedge after
  // which the new curve should be inserted.
  DHalfedge  *prev1 = _locate_around_vertex (_vertex (v), cv);

  CGAL_assertion_msg
    (prev1 != NULL,
     "The inserted curve should not exist in the arrangement.");

  // Perform the insertion.
  DHalfedge  *new_he = _insert_from_vertex (cv,
                                           prev1, v2);

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
    (traits->equal_2_object() (prev->target()->point(),
                               traits->construct_min_vertex_2_object()(cv)),
     "The input halfedge's target should be the left curve endpoint.");
  CGAL_precondition_msg
    (_locate_around_vertex(_vertex (prev->target()), cv) == _halfedge (prev),
     "In the clockwise order of curves around the vertex, "
     " cv must succeeds the curve of prev.");

  // Create a new vertex associated with the curve's right endpoint.
  // We also notify the observers on the creation of this vertex.
  Stored_point_2  *p2 =
    _new_point (traits->construct_max_vertex_2_object()(cv));

  _notify_before_create_vertex (*p2);

  DVertex         *v2 = dcel.new_vertex();

  v2->set_point (p2);

  _notify_after_create_vertex (Vertex_handle (v2));

  // Perform the insertion.
  DHalfedge  *new_he = _insert_from_vertex (cv,
                                            _halfedge (prev), v2);

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
    (traits->equal_2_object() (v->point(),
                               traits->construct_max_vertex_2_object()(cv)),
     "The input vertex should be the right curve endpoint.");

  // Create a new vertex associated with the curve's left endpoint.
  // We also notify the observers on the creation of this vertex.
  Stored_point_2  *p1 =
    _new_point (traits->construct_min_vertex_2_object()(cv));

  _notify_before_create_vertex (*p1);

  DVertex         *v1 = dcel.new_vertex();

  v1->set_point (p1);

  _notify_after_create_vertex (Vertex_handle (v1));

  // Check if the given vertex is isolated.
  if (v->is_isolated())
  {
    // The given vertex is an isolated one: We should in fact insert the curve
    // in the interior of the face containing this vertex.
    DVertex    *v2 = _vertex (v);
    DHalfedge  *p_he = v2->halfedge();
    DFace      *p_f = (p_he != NULL) ?  p_he->face() : un_face;

    // v2 will not be isolated any more - remove it from the isolated vertices
    // container of its containing face.
    _find_and_erase_isolated_vertex (p_f, v2);

    // Create the edge connecting the two vertices.
    DHalfedge  *new_he = _insert_in_face_interior (cv, p_f, v1, v2);

    // Return a handle to the new halfedge directed from left to right.
    return (Halfedge_handle (new_he));
  }

  // Go over the incident halfedges around v and find the halfedge after
  // which the new curve should be inserted.
  DHalfedge  *prev2 = _locate_around_vertex (_vertex (v), cv);

  CGAL_assertion_msg
    (prev2 != NULL,
     "The inserted curve should not exist in the arrangement.");

  // Perform the insertion.
  DHalfedge  *new_he = _insert_from_vertex (cv,
                                           prev2, v1);

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
    (traits->equal_2_object() (prev->target()->point(),
                               traits->construct_max_vertex_2_object()(cv)),
     "The input halfedge's target should be the right curve endpoint.");

  CGAL_precondition_msg
    (_locate_around_vertex(_vertex (prev->target()), cv) == _halfedge (prev),
     "In the clockwise order of curves around the vertex, "
     "cv must succeeds the curve of prev.");

  // Create a new vertex associated with the curve's left endpoint.
  // We also notify the observers on the creation of this vertex.
  Stored_point_2  *p1 =
    _new_point (traits->construct_min_vertex_2_object()(cv));

  _notify_before_create_vertex (*p1);

  DVertex         *v1 = dcel.new_vertex();

  v1->set_point (p1);

  _notify_after_create_vertex (Vertex_handle (v1));

  // Perform the insertion.
  DHalfedge  *new_he = _insert_from_vertex (cv,
                                            _halfedge (prev), v1);

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
    ((traits->equal_2_object()(v1->point(),
                               traits->construct_min_vertex_2_object()(cv)) &&
      traits->equal_2_object()(v2->point(),
                               traits->construct_max_vertex_2_object()(cv))) ||
     (traits->equal_2_object()(v2->point(),
                               traits->construct_min_vertex_2_object()(cv)) &&
      traits->equal_2_object()(v1->point(),
                               traits->construct_max_vertex_2_object()(cv))),
     "The input vertices should match the curve endpoints.");

  // Check whether one of the vertices is isolated.
  if (v1->is_isolated())
  {
    // Get the face containing the isolated vertex v1.
    DVertex    *p_v1 = _vertex (v1);
    DHalfedge  *he1 = p_v1->halfedge();
    DFace      *f1 = (he1 != NULL) ?  he1->face() : un_face;

    // v1 will not be isolated any more - remove it from the isolated vertices
    // container of its containing face.
    _find_and_erase_isolated_vertex (f1, p_v1);

    if (v2->is_isolated())
    {
      // Both end-vertices are isolated. Make sure they are contained inside
      // the same face.
      DVertex    *p_v2 = _vertex (v2);

      CGAL_assertion_code (
        DHalfedge  *he2 = p_v2->halfedge();
        DFace      *f2 = (he2 != NULL) ?  he2->face() : un_face;
      );

      CGAL_assertion_msg
        (f1 == f2,
         "The inserted curve should not intersect the existing arrangement.");

      // v2 will not be isolated any more - remove it from the isolated
      // vertices container of its containing face.
      _find_and_erase_isolated_vertex (f1, p_v2);

      // Perform the insertion.
      DHalfedge  *new_he = _insert_in_face_interior (cv, f1,
                                                     p_v1, p_v2);

      return (Halfedge_handle (new_he));
    }

    // Go over the incident halfedges around v2 and find the halfedge after
    // which the new curve should be inserted.
    DHalfedge  *prev2 = _locate_around_vertex (_vertex (v2), cv);

    CGAL_assertion_msg
      (prev2 != NULL,
       "The inserted curve should not exist in the arrangement.");

    CGAL_assertion_msg
      (f1 == prev2->face(),
       "The inserted curve should not intersect the existing arrangement.");

    // Perform the insertion. Note that the returned halfedge is directed
    // toward v1 (and not toward v2), so we return its twin.
    DHalfedge  *new_he = _insert_from_vertex (cv,
                                              prev2, p_v1);

    return (Halfedge_handle (new_he->opposite()));
  }
  else if (v2->is_isolated())
  {
    // Get the face containing the isolated vertex v2.
    DVertex    *p_v2 = _vertex (v2);
    DHalfedge  *he2 = p_v2->halfedge();
    DFace      *f2 = (he2 != NULL) ?  he2->face() : un_face;

    // v2 will not be isolated any more - remove it from the isolated vertices
    // container of its containing face.
    _find_and_erase_isolated_vertex (f2, p_v2);

    // Go over the incident halfedges around v1 and find the halfedge after
    // which the new curve should be inserted.
    DHalfedge  *prev1 = _locate_around_vertex (_vertex (v1), cv);

    CGAL_assertion_msg
      (prev1 != NULL,
       "The inserted curve should not exist in the arrangement.");

    CGAL_assertion_msg
      (f2 == prev1->face(),
       "The inserted curve should not intersect the existing arrangement.");

    // Perform the insertion.
    DHalfedge  *new_he = _insert_from_vertex (cv,
                                              prev1, p_v2);

    return (Halfedge_handle (new_he));
  }

  // Go over the incident halfedges around v1 and v2 and find the two
  // halfedges after which the new curve should be inserted, respectively.
  DHalfedge  *prev1 = _locate_around_vertex (_vertex (v1), cv);
  DHalfedge  *prev2 = _locate_around_vertex (_vertex (v2), cv);

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
    ((traits->equal_2_object()(prev1->target()->point(),
                               traits->construct_min_vertex_2_object()(cv)) &&
      traits->equal_2_object()(v2->point(),
                               traits->construct_max_vertex_2_object()(cv))) ||
     (traits->equal_2_object()(v2->point(),
                               traits->construct_min_vertex_2_object()(cv)) &&
      traits->equal_2_object()(prev1->target()->point(),
                               traits->construct_max_vertex_2_object()(cv))),
     "The input vertex and target point should match the curve endpoints.");

  // Check whether v2 is an isolated vertex.
  if (v2->is_isolated())
  {
    // Get the face containing the isolated vertex v2.
    DVertex    *p_v2 = _vertex (v2);
    DHalfedge  *he2 = p_v2->halfedge();
    DFace      *f2 = (he2 != NULL) ?  he2->face() : un_face;

    CGAL_assertion_msg
      (f2 == (_halfedge (prev1))->face(),
       "The inserted curve should not intersect the existing arrangement.");

    // v2 will not be isolated any more - remove it from the isolated vertices
    // container of its containing face.
    _find_and_erase_isolated_vertex (f2, p_v2);

    // Perform the insertion.
    DHalfedge  *new_he = _insert_from_vertex (cv,
                                              _halfedge (prev1), p_v2);

    return (Halfedge_handle (new_he));
  }

  // Go over the incident halfedges around v2 and find the halfedge after
  // which the new curve should be inserted.
  DHalfedge  *prev2 = _locate_around_vertex (_vertex (v2), cv);

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
  const int   dist1 = _halfedge_distance (p_prev1, p_prev2);
  bool        prev1_before_prev2 = true;

  CGAL_assertion (dist1 != 0);

  if (dist1 > 0)
  {
    // If prev2 is not reachable from prev1, the insertion of the new curve
    // does not generate a new face, so the way we send these halfedge
    // pointers to the auxiliary function _insert_at_vertices() does not
    // matter. But in our case, where the two halfedges are reachable from
    // one another, a new face will be created and we have to arrange prev1
    // and prev2 so that the new face will be incident to the correct
    // halfedge. To do this, we check whether prev1 lies inside the new face
    // we are about to create (or alternatively, whether prev2 does not lie
    // inside this new face).
    const int  dist2 = _halfedge_distance (p_prev2, p_prev1);

    CGAL_assertion (dist2 > 0);

    prev1_before_prev2 = (dist1 > dist2) ?
      (_is_inside_new_face (p_prev1, p_prev2, cv)) :
      (! _is_inside_new_face (p_prev2, p_prev1, cv));
  }

  // Perform the insertion.
  bool        new_face_created = false;
  DHalfedge  *new_he = (prev1_before_prev2) ?
    _insert_at_vertices (cv, p_prev1, p_prev2, new_face_created) :
    _insert_at_vertices (cv, p_prev2, p_prev1, new_face_created);

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
// Replace the x-monotone curve associated with the given edge.
//
template<class Traits, class Dcel>
typename Arrangement_2<Traits,Dcel>::Halfedge_handle
Arrangement_2<Traits,Dcel>::modify_edge (Halfedge_handle e,
                                         const X_monotone_curve_2& cv)
{
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
  // Get the split halfedge and its twin, its source and target.
  DHalfedge       *he1 = _halfedge (e);
  DHalfedge       *he2 = he1->opposite();
  const Point_2&   source = he2->vertex()->point();
  CGAL_precondition_code (
    const Point_2&   target = he1->vertex()->point();
  );

  // Determine the point where we split the halfedge. We also determine which
  // curve should be associated with he1 (and he2), which is the curve who
  // has an endpoint that equals e's source, and which should be associated
  // with the new pair of halfedges we are about to split (the one who has
  // an endpoint which equals e's target).
  const Point_2&  cv1_left = traits->construct_min_vertex_2_object() (cv1);
  const Point_2&  cv1_right = traits->construct_max_vertex_2_object() (cv1);
  const Point_2&  cv2_left = traits->construct_min_vertex_2_object() (cv2);
  CGAL_precondition_code (
    const Point_2&  cv2_right = traits->construct_max_vertex_2_object() (cv2);
  );
  const Point_2               *split_p;
  const X_monotone_curve_2    *p_cv1;
  const X_monotone_curve_2    *p_cv2;

  CGAL_precondition_msg (traits->equal_2_object() (cv1_right, cv2_left) ||
                         traits->equal_2_object() (cv2_right, cv1_left),


                         "Curves do not share a common endpoint.");

  if (traits->equal_2_object() (cv1_right, cv2_left))
  {
    split_p = &cv1_right;

    if (traits->equal_2_object() (source, cv1_left))
    {
      CGAL_precondition (traits->equal_2_object() (target, cv2_right));

      p_cv1 = &cv1;
      p_cv2 = &cv2;
    }
    else
    {
      CGAL_precondition (traits->equal_2_object() (source, cv2_right) &&
                         traits->equal_2_object() (target, cv1_left));

      p_cv1 = &cv2;
      p_cv2 = &cv1;
    }
  }
  else
  {
    CGAL_precondition (traits->equal_2_object() (cv2_right, cv1_left));

    split_p = &cv1_left;

    if (traits->equal_2_object() (source, cv2_left))
    {
      CGAL_precondition (traits->equal_2_object() (target, cv1_right));

      p_cv1 = &cv2;
      p_cv2 = &cv1;
    }
    else
    {
      CGAL_precondition (traits->equal_2_object() (source, cv1_right) &&
                         traits->equal_2_object() (target, cv2_left));


      p_cv1 = &cv1;
      p_cv2 = &cv2;
    }
  }

  return (Halfedge_handle (_split_edge (he1, *split_p, *p_cv1, *p_cv2)));
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
  CGAL_precondition_msg
    ((traits->equal_2_object() (traits->construct_min_vertex_2_object() (cv),
                                e1->source()->point()) &&
      traits->equal_2_object() (traits->construct_max_vertex_2_object() (cv),
                                e2->target()->point())) ||
     (traits->equal_2_object() (traits->construct_min_vertex_2_object() (cv),
                                e1->source()->point()) &&
      traits->equal_2_object() (traits->construct_max_vertex_2_object() (cv),
                                e2->source()->point())) ||
     (traits->equal_2_object() (traits->construct_min_vertex_2_object() (cv),
                                e1->target()->point()) &&
      traits->equal_2_object() (traits->construct_max_vertex_2_object() (cv),
                                e2->source()->point())) ||
     (traits->equal_2_object() (traits->construct_min_vertex_2_object() (cv),
                                e1->target()->point()) &&
      traits->equal_2_object() (traits->construct_max_vertex_2_object() (cv),
                                e2->target()->point())) ||
     (traits->equal_2_object() (traits->construct_min_vertex_2_object() (cv),
                                e2->source()->point()) &&
      traits->equal_2_object() (traits->construct_max_vertex_2_object() (cv),
                                e1->source()->point())) ||
     (traits->equal_2_object() (traits->construct_min_vertex_2_object() (cv),
                                e2->source()->point()) &&
      traits->equal_2_object() (traits->construct_max_vertex_2_object() (cv),
                                e1->target()->point())) ||
     (traits->equal_2_object() (traits->construct_min_vertex_2_object() (cv),
                                e2->target()->point()) &&
      traits->equal_2_object() (traits->construct_max_vertex_2_object() (cv),
                                e1->source()->point())) ||
     (traits->equal_2_object() (traits->construct_min_vertex_2_object() (cv),
                                e2->target()->point()) &&
      traits->equal_2_object() (traits->construct_max_vertex_2_object() (cv),
                                e1->target()->point())),
     "The curve endpoints do not match the merged edges.");

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

  CGAL_precondition_msg(he1->next()->opposite() == he4 &&
                        he4->next()->opposite() == he1,
                        "The degree of the deleted vertex is greater than 2.");

  // Keep pointers to the incident faces of the two halfedges he3 and he2,
  // pointing at the end vertices of the merged halfedge.
  DFace      *f1 = he3->face();
  DFace      *f2 = he2->face();

  // Notify the observers that we are about to merge an edge.
  _notify_before_merge_edge (e1, e2, cv);

  // As he1 and he2 will evetually represent the merged edge, while he3 and he4
  // will be deleted, check if the deleted halfedges are represantatives of a
  // face boundary or a hole inside these faces. If so, replace he3 by he1 and
  // he4 by he2. Note that as we just change the hole representatives, we do
  // not have to notify the observers about the change.
  if (f1->halfedge() == he3)
  {
    f1->set_halfedge (he1);

    // Go over all isolated vertices inside f1 and update their fictitious
    // incident face (which used to be he3) to be he1.
    DIsolated_vertices_iter     iso_vert_it;

    for (iso_vert_it = f1->isolated_vertices_begin();
         iso_vert_it != f1->isolated_vertices_end(); ++iso_vert_it)
    {
      iso_vert_it->set_halfedge (he1);
    }
  }
  else
  {
    if (_find_and_erase_hole (f1, he3))
      f1->add_hole (he1);
  }

  if (f2->halfedge() == he4)

  {
    f2->set_halfedge (he2);

    // Go over all isolated vertices inside f2 and update their fictitious
    // incident face (which used to be he4) to be he2.
    DIsolated_vertices_iter     iso_vert_it;

    for (iso_vert_it = f2->isolated_vertices_begin();
         iso_vert_it != f2->isolated_vertices_end(); ++iso_vert_it)
    {
      iso_vert_it->set_halfedge (he2);
    }
  }
  else
  {
    if (_find_and_erase_hole (f2, he4))
      f2->add_hole (he2);
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
Arrangement_2<Traits,Dcel>::remove_edge (Halfedge_handle e)
{
  DHalfedge   *he1 = _halfedge (e);
  DHalfedge   *he2 = he1->opposite();
  DFace       *f;

  if (he1->face() != he2->face() ||
      he1->next() == he2 || he2->next() == he1)
  {
    // Either the removal of he1 (and its twin halfedge) will cause the two
    // incident faces to merge, or these two halfedges form an "antenna".
    // In either case, it does not matter which halfedge we send to the
    // auxiliary function _remove_edge().
    f = _remove_edge (he1);
  }
  else
  {
    // In this case if a new hole will be created by the removal of he1 (and
    // its twin halfedge). We determine the halfedge (he1 or he2) that
    // points at the new hole in order to send it to _remove_hole().
    // We begin by locating the leftmost point along the path from he1 to he2
    // and the leftmost point along the path from he2 to he1.
    typename Traits_wrapper_2::Compare_xy_2  compare_xy =
                                                 traits->compare_xy_2_object();
    DHalfedge      *ccb1 = he1->next();
    DHalfedge      *ccb2 = he2->next();
    const Point_2  *p_min1 = &(he1->vertex()->point());
    const Point_2  *p_min2 = &(he2->vertex()->point());

    do
    {
      if (compare_xy (ccb1->vertex()->point(), *p_min1) == SMALLER)
        p_min1 = &(ccb1->vertex()->point());

      ccb1 = ccb1->next();

    } while (ccb1 != he2);

    do
    {
      if (compare_xy (ccb2->vertex()->point(), *p_min2) == SMALLER)
        p_min2 = &(ccb2->vertex()->point());

      ccb2 = ccb2->next();

    } while (ccb2 != he1);

    // Compare the two leftmost points: p_min2 lies to the left of p_min1
    // if and only if he1 points at the hole we are about to create.
    if (traits->compare_x_2_object() (*p_min2, *p_min1) == SMALLER)
      f = _remove_edge (he1);
    else
      f = _remove_edge (he2);
  }

  return (Face_handle (f));
}

//-----------------------------------------------------------------------------
// Protected member functions (for internal use).
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Locate the place for the given curve around the given vertex.
//
template<class Traits, class Dcel>
typename Arrangement_2<Traits,Dcel>::DHalfedge*
Arrangement_2<Traits,Dcel>::_locate_around_vertex
    (DVertex *v,
     const X_monotone_curve_2& cv) const
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
  typename Traits_wrapper_2::Is_between_cw_2  is_between_cw =
                                              traits->is_between_cw_2_object();

  bool       eq_curr, eq_next;

  while (!is_between_cw (cv,
                         curr->curve(),
                         next->curve(),
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
int Arrangement_2<Traits,Dcel>::_halfedge_distance (const DHalfedge *e1,
                                                    const DHalfedge *e2) const
{
  // In case both given halfedges are the same:
  if (e1 == e2)
    return (0);

  // Traverse the halfedge chain from e1 and try to locate e2.
  const DHalfedge   *curr = e1->next();
  int                dist = 1;

  while (curr != e2)
  {
    // If we have returned to e1, e2 is not reachable from e1.
    if (curr == e1)
      return (-1);

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
  typename Traits_wrapper_2::Compare_xy_2  compare_xy =
                                                 traits->compare_xy_2_object();

  const DHalfedge   *left_edge = NULL;
  const DHalfedge   *first = prev2->next();
  const DHalfedge   *curr = first;
  const DHalfedge   *last = prev1->next();
  const Point_2     *p_min = &(prev2->vertex()->point());

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
    if (compare_xy (curr->vertex()->point(), *p_min) != LARGER)
    {
      left_edge = curr;
      p_min = &(curr->vertex()->point());
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
                                                   *p_min) == LARGER);
}

//-----------------------------------------------------------------------------
// Determine whether a given point lies within the region bounded by
// a boundary of a connected component.
//
template<class Traits, class Dcel>
bool Arrangement_2<Traits,Dcel>::_point_is_in (const Point_2& p,
                                               const DHalfedge *he) const
{
  // Keep a counter of the number of x-monotone curves that intersect an upward
  // vertical emanating from p (except for some degenerate cases that are
  // explained below).
  unsigned int              n_ray_intersections = 0;

  // Find the first halfedge whose associated curve is non-vertical along the
  // boundary of the given connected component.
  typename Traits_wrapper_2::Is_vertical_2    is_vertical =
                                            traits->is_vertical_2_object();
  typename Traits_wrapper_2::Equal_2          equal =
                                            traits->equal_2_object();
  typename Traits_wrapper_2::Compare_x_2      compare_x =
                                            traits->compare_x_2_object();
  typename Traits_wrapper_2::Compare_y_at_x_2 compare_y_at_x =
                                            traits->compare_y_at_x_2_object();

  const DHalfedge    *curr = he;
  Comparison_result   res_source;
  Comparison_result   res_target;

  do
  {
    curr = curr->next();
  } while (curr != he && is_vertical (curr->curve()));

  if (curr == he && is_vertical (he->curve()))
  {
    // In this case the entire component is comprised of vertical segments,
    // so it has any empty interior and p cannot lie inside it.
    return (false);
  }

  // Go over all curves of the boundary, starting from the non-vertical curve

  // we have located, and count those which are above p.
  const DHalfedge   *first = curr;

  do
  {
    // Compare the x-coordinates of the current halfedge's endpoint with the
    // query point.
    res_source = compare_x (curr->opposite()->vertex()->point(), p);
    res_target = compare_x (curr->vertex()->point(), p);

    // If the query point is a boundary vertex, it is obviuosly not in the
    // interior the component.
    if (res_target == EQUAL && equal (curr->vertex()->point(), p))
      return (false);

    // Check the query point is in the x-range (source, target] of this
    // the current curve and lies below it. Note that this condition rules
    // out vertical segment (where the x-coordinates of the source and target
    // are equal - and indeed, we do not have count intersections (actually
    // these are overlaps) of the vertical ray we shoot with vertical
    // segments along the boundary.
    if (res_source != res_target && res_source != EQUAL &&
        (compare_y_at_x (p, curr->curve()) == SMALLER))
    {
      // In the degenerate case where p lies below the target vertex of
      // the current halfedge, we have to be a bit careful:
      if (res_target == EQUAL)
      {
        // Locate the next halfedge along the boundary that does not
        // contain a vertical segment.
        const DHalfedge   *next_non_vert = curr;

        do
        {
          next_non_vert = next_non_vert->next();

          CGAL_assertion_msg (next_non_vert != curr,
                              "Infinite loop in _point_is_in().");


        } while (is_vertical (next_non_vert->curve()));

        // In case the source of the current curve and the target of
        // the next non-vertical curve lie on opposite sides of the
        // ray we shoot from p (case (a)), we have to count an
        // intersection. Otherwise, we have a "tangency" with the ray
        // (case (b)) and it is not necessary to count it.
        //
        //            +--------+              +                 .
        //            |   next                 \ next           .
        //            |                         \               .
        //            +                          +              .
        //           /                          /               .
        //     curr /                     curr /                .
        //         /                          /                 .
        //        +  (.)p                    +  (.)p            .
        //
        //          (a)                        (b)
        //
        res_target = compare_x (next_non_vert->vertex()->point(), p);

        CGAL_assertion (res_source != EQUAL && res_target != EQUAL);

        if (res_source != res_target)
          n_ray_intersections++;
      }
      else
      {
        // In this case p lies under the interior of the current x-montone
        // curve, so the vertical ray we shoot intersects it exactly once.
        n_ray_intersections++;
      }
    }

    // Proceed to the next halfedge along the component boundary.
    curr = curr->next();

  } while (curr != first);

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
  // Get a halfedge lying on the boundary of the hole.
  DHalfedge   *first = *hole;

  // Notify the observers that we are about to move a hole.
  Ccb_halfedge_circulator   circ = (Halfedge_handle(first))->ccb();

  _notify_before_move_hole (Face_handle (from_face),
                            Face_handle (to_face),
                            circ);

  // Go over all halfedges along the boundary of the hole and modify their
  // incident face.
  DHalfedge   *curr = first;

  do
  {
    curr->set_face (to_face);
    curr = curr->next();

  } while (curr != first);

  // Move the hole from the first face to the other.
  to_face->add_hole (first);
  from_face->erase_hole (hole);

  // Notify the observers that we have moved the hole.
  _notify_after_move_hole (circ);

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

  // Notify the observers that we are about to move an isolated vertex.
  Vertex_handle   vh (v);

  _notify_before_move_isolated_vertex (Face_handle (from_face),
                                       Face_handle (to_face),
                                       vh);

  // Set the new halfedge of the vertex (may be NULL if to_face is unbounded).
  v->set_halfedge (to_face->halfedge());

  // Move the isolated vertex from the first face to the other.
  to_face->add_isolated_vertex (v);
  from_face->erase_isolated_vertex (vit);

  // Notify the observers that we have moved the isolated vertex.
  _notify_after_move_isolated_vertex (vh);

  return;
}

//-----------------------------------------------------------------------------
// Check whether the given halfedge lies on the outer boundary of the given
// face.
//
template<class Traits, class Dcel>
typename Arrangement_2<Traits,Dcel>::DHalfedge*
Arrangement_2<Traits,Dcel>::_is_on_outer_boundary (DFace *f,
                                                   DHalfedge *e) const
{
  DHalfedge   *occb = f->halfedge();
  DHalfedge   *curr = occb;

  // In case f is the unbounded face, it has no outer boundary.
  if (occb == NULL)
    return (NULL);

  // Traverse the outer boundary of f and look for the given halfedge.
  do
  {
    if (curr == e)
      // We have located e, return a representative halfedge on the outer
      // boundary of the face.
      return (occb);

    curr = curr->next();
  } while (curr != occb);

  // If we reached here, we have not located e on the outer boundary of f.
  return (NULL);
}

//-----------------------------------------------------------------------------
// Check whether the given halfedge lies on the inner boundary of the given
// face.
//
template<class Traits, class Dcel>
typename Arrangement_2<Traits,Dcel>::DHalfedge*
Arrangement_2<Traits,Dcel>::_is_on_inner_boundary (DFace *f,
                                                   DHalfedge *e) const
{
  // Go over the holes of the given face.
  DHoles_iter  holes_it;
  DHalfedge   *iccb;
  DHalfedge   *curr;

  for (holes_it = f->holes_begin(); holes_it != f->holes_end(); ++holes_it)
  {
    iccb = *holes_it;
    curr = iccb;

    // Traverse the outer boundary of the current hole and look for the given
    // halfedge.
    do
    {
      if (curr == e)
        // We have located e, return a representative halfedge on the outer
        // boundary of the current hole.
        return (iccb);

      curr = curr->next();
    } while (curr != iccb);
  }

  // If we reached here, we have not located e on the inner boundary of f.
  return (NULL);
}

//-----------------------------------------------------------------------------
// Find the hole represented by a given halfedge from the holes container
// of a given face and earse this hole once it is found.
//
template<class Traits, class Dcel>
bool Arrangement_2<Traits,Dcel>::_find_and_erase_hole (DFace *f, DHalfedge* e)
{
  DHoles_iter     holes_it;

  for (holes_it = f->holes_begin(); holes_it != f->holes_end(); ++holes_it)
  {
    if (*holes_it == e)
    {
      // We have found the hole - erase it.
      f->erase_hole (holes_it);
      return (true);
    }
  }

  // We have not located the hole in the given face.
  return (false);
}

//-----------------------------------------------------------------------------
// Find the vertex in the isolated vertices container of a given face and
// earse this vertex once it is found.
//
template<class Traits, class Dcel>
bool Arrangement_2<Traits,Dcel>::_find_and_erase_isolated_vertex (DFace *f,
                                                                  DVertex *v)
{
  DIsolated_vertices_iter     iso_vert_it;

  for (iso_vert_it = f->isolated_vertices_begin();
       iso_vert_it != f->isolated_vertices_end(); ++iso_vert_it)
  {
    if (&(*iso_vert_it) == v)
    {
      // We have found the isolated vertex - erase it.
      f->erase_isolated_vertex (iso_vert_it);
      return (true);
    }
  }

  // We have not located the isolated vertex in the given face.
  return (false);
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
     DVertex *v1, DVertex *v2)
{
  // Notify the observers that we are about to create a new edge.
  _notify_before_create_edge (cv, Vertex_handle (v1), Vertex_handle (v2));

  // Create a pair of twin halfedges connecting the two vertices,
  // and link them together to form a new connected component.
  DHalfedge       *he1 = dcel.new_edge();
  DHalfedge       *he2 = he1->opposite();
  Stored_curve_2  *dup_cv = _new_curve (cv);

  he1->set_curve (dup_cv);

  he1->set_next (he2);
  he1->set_vertex (v1);
  he1->set_face (f);

  he2->set_next (he1);
  he2->set_vertex (v2);
  he2->set_face (f);

  // Assign the incident halfedges of the two new vertices.
  v1->set_halfedge (he1);
  v2->set_halfedge (he2);

  // Create a handle to the new halfedge pointing at the curve target.
  Halfedge_handle   hh (he2);

  // Notify the observers that we have created a new edge.
  _notify_after_create_edge (hh);

  // Notify the observers that we are about to form a new hole inside f.
  _notify_before_add_hole (Face_handle (f), hh);

  // Initiate a new hole inside the given face.
  f->add_hole (he2);

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
                                                 DVertex *v)
{
  // Get the incident face of the previous halfedge . Note that this will also
  // be the incident face of the two new halfedges we are about to create.
  DFace           *p_f = prev->face();

  // The first vertex is the one that the prev halfedge points to.
  // The second vertex is given by v.
  DVertex         *v1 = prev->vertex();
  DVertex         *v2 = v;

  // Notify the observers that we are about to create a new edge.
  _notify_before_create_edge (cv, Vertex_handle (v1), Vertex_handle (v2));

  // Create a pair of twin halfedges connecting the two vertices,
  // and associate them with the given curve.
  DHalfedge       *he1 = dcel.new_edge();
  DHalfedge       *he2 = he1->opposite();
  Stored_curve_2  *dup_cv = _new_curve (cv);

  he1->set_curve (dup_cv);

  he1->set_vertex (v1);
  he2->set_vertex (v2);

  he1->set_face (p_f);
  he2->set_face (p_f);

  // Associate the incident halfedge of the new vertex.
  v2->set_halfedge (he2);

  // Link the new halfedges around the existing vertex v1.
  he2->set_next (he1);
  he1->set_next (prev->next());

  prev->set_next (he2);

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
                                                 bool& new_face)
{
  // Get the vertices that match cv's endpoints.
  DVertex    *v1 = prev1->vertex();
  DVertex    *v2 = prev2->vertex();

  // Get the incident face (which should be the same for the two halfedges).
  DFace      *f = prev1->face();

  CGAL_precondition_msg
    (prev1->face() == prev2->face(),
     "The two halfedges must share the same incident face.");

  // Locate representative halfedges for the connected components of both
  // given halfedges.
  DHalfedge  *ccb1 = _is_on_outer_boundary (f, prev1);
  bool        ccb1_is_inner = false;

  if (ccb1 == NULL)
  {
    ccb1_is_inner = true;
    ccb1 = _is_on_inner_boundary (f, prev1);
  }
  CGAL_assertion (ccb1 != NULL);

  DHalfedge  *ccb2 = _is_on_outer_boundary (f, prev2);
  bool        ccb2_is_inner = false;

  if (ccb2 == NULL)
  {
    ccb2_is_inner = true;
    ccb2 = _is_on_inner_boundary (f, prev2);
  }
  CGAL_assertion (ccb2 != NULL);

  // Notify the observers that we are about to create a new edge.
  _notify_before_create_edge (cv, Vertex_handle (v1), Vertex_handle (v2));

  // Create a pair of twin halfedges connecting v1 and v2 and associate them
  // with the given curve.
  DHalfedge       *he1 = dcel.new_edge();
  DHalfedge       *he2 = he1->opposite();
  Stored_curve_2  *dup_cv = _new_curve (cv);

  he1->set_curve (dup_cv);

  he1->set_vertex (v1);
  he2->set_vertex (v2);

  // Connect the new halfedges to the existing halfegdes around the two
  // incident vertices.
  he1->set_next (prev1->next());
  he2->set_next (prev2->next());

  prev1->set_next (he2);
  prev2->set_next (he1);

  // Check whether prev1 and prev2 used to belong to the same connected
  // component or to different components.
  const bool ccbs_equal = (ccb1 == ccb2);

  if (!ccbs_equal)
  {
    // In case we have connected two disconnected components, no new face
    // is created. In this case, f is the incident face of the two new
    // halfedges.
    he1->set_face (f);
    he2->set_face (f);

    new_face = false;

    // Check whether both halfedges are inner components (hole) in the same
    // face, or whether one is a hole and the other is on the outer boundary
    // of the face. 
    Face_handle       fh (f);

    if (ccb1_is_inner && ccb2_is_inner)
    {
      // Notify the observers that we are about to merge two holes in the face.
      _notify_before_merge_hole (fh,
				 (Halfedge_handle(he1))->ccb(),
				 (Halfedge_handle(he2))->ccb(),
				 Halfedge_handle (he1));

      // Remove the hole ccb1 represents. As ccb2 is also a hole, the two
      // holes are merged into one.
      _find_and_erase_hole (f, ccb1);

      // Notify the observers that we have merged the two holes.
      _notify_after_merge_hole (fh,
				(Halfedge_handle(he1))->ccb());
    }
    else
    {
      CGAL_assertion_msg (ccb1_is_inner || ccb2_is_inner,
                          "Cannot merge two outer components.");

      DHalfedge  *ccb_to_remove = ccb1_is_inner ? ccb1 : ccb2;

      // Notify the observers that we are about to remove a fole from the face.
      _notify_before_remove_hole (fh,
				  (Halfedge_handle(ccb_to_remove))->ccb());

      // Remove the hole from the face, as we have just connected this hole to
      // the outer boundary of its incident face.
      _find_and_erase_hole (f, ccb_to_remove);

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

    // Set the incident faces of the two new halfedges.The new face should
    // become the incident face of he2, while he1 is associated with the
    // existing face.
    he2->set_face (new_f);
    he1->set_face (f);

    // Set the incident face of all halfedges along he2's (new) connected
    // component to be the new face.
    bool   is_hole;

    for (ccb2 = he2->next(); ccb2 != he2; ccb2 = ccb2->next())
      ccb2->set_face (new_f);

    if (! ccb1_is_inner)
    {
      // As the outer boundary edge of the exisitng face f may be one of the
      // halfedges along the boundary of the new face, we set it to he1.
      f->set_halfedge (he1);
      is_hole = false;
    }
    else
    {
      // The prev1 and prev2 halfedges belong to a the inner boundary of their
      // component. In this case we have to remove the hole that the used to
      // belong to the face f and replace it with the new hole, represnted
      // by he1.
      _find_and_erase_hole (f, ccb1);
      f->add_hole (he1);
      is_hole = true;
    }

    // Notify the observers that we have split the face.
    _notify_after_split_face (fh,
                              Face_handle (new_f),
                              is_hole);
  }

  // Notify the observers that we have created a new edge.
  _notify_after_create_edge (Halfedge_handle (he1));

  // Return the halfedge directed from v1 to v2.
  return (he2);
}

//-----------------------------------------------------------------------------
// Relocate all holes and isolated vertices to their proper position,
// immediately after a face has split due to the insertion of a new halfedge.
//
template<class Traits, class Dcel>
void Arrangement_2<Traits,Dcel>::_relocate_in_new_face (DHalfedge *new_he)
{
  // The given halfedge points to the new face, while its twin points to the
  // old face (the one that has just been split).
  DFace        *new_face = new_he->face();
  DFace        *old_face = new_he->opposite()->face();

  CGAL_assertion (new_face != old_face);

  // Examine the holes inside the existing old face and move the relevant
  // ones into the new face.
  DHoles_iter   holes_it = old_face->holes_begin();
  DHoles_iter   hole_to_move;

  while (holes_it != old_face->holes_end())
  {
    // Check whether the current hole is inside new face.
    if (_point_is_in((*holes_it)->vertex()->point(), new_he))
    {
      // We increment the holes itrator before moving the hole, because
      // this operation invalidates the iterator.
      hole_to_move  = holes_it;
      ++holes_it;

      // Move the hole.
      _move_hole (old_face, new_face, hole_to_move);
    }
    else
    {
      ++holes_it;
    }
  }

  // Examine the isolated vertices inside the existing old face and move the
  // relevant ones into the new face.
  DIsolated_vertices_iter    iso_verts_it;
  DIsolated_vertices_iter    iso_vert_to_move;

  iso_verts_it = old_face->isolated_vertices_begin();
  while (iso_verts_it != old_face->isolated_vertices_end())
  {
    // Check whether the isolated vertex lies inside the new face.
    if (_point_is_in (iso_verts_it->point(), new_he))
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
  // Get the split halfedge and its twin, its source and target.
  DHalfedge       *he1 = e;
  DHalfedge       *he2 = he1->opposite();

  // Allocate a new vertex and associate it with the split point.
  // We also notify the observers on the creation of this vertex.
  Stored_point_2  *dup_p = _new_point (p);

  _notify_before_create_vertex (p);

  DVertex         *v = dcel.new_vertex();

  v->set_point (dup_p);

  _notify_after_create_vertex (Vertex_handle (v));

  // Notify the observers that we are about to split an edge.
  _notify_before_split_edge (Halfedge_handle (e), cv1, cv2);

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

  he3->set_face (he1->face());
  he3->set_vertex (he1->vertex());
  he4->set_vertex (v);
  he4->set_next (he2);
  he4->set_face (he2->face());

  if (he1->vertex()->halfedge() == he1)
    // If he1 is the incident halfedge to its target, he3 replaces it.
    he1->vertex()->set_halfedge (he3);

  // Update the properties of the twin halfedges we have just split.
  he1->set_next(he3);
  he1->set_vertex(v);

  // Associate cv1 with he1 (and its twin). We allocate a new curve for cv2
  // and associate it with he3 (and its twin).
  Stored_curve_2  *dup_cv2 = _new_curve (cv2);

  he1->curve() = cv1;
  he3->set_curve (dup_cv2);

  // Notify the observers that we have split an edge into two.
  _notify_after_split_edge (Halfedge_handle (he1), Halfedge_handle (he3));

  // Return a handle for one of the existing halfedge that is incident to the
  // split point.
  return (he1);
}

//-----------------------------------------------------------------------------
// Remove a pair of twin halfedges from the arrangement.
// In case the removal causes the creation of a new hole, the given halfedge
// should point at this hole.
//
template<class Traits, class Dcel>
typename Arrangement_2<Traits,Dcel>::DFace*
Arrangement_2<Traits,Dcel>::_remove_edge (DHalfedge *e)
{
  // Get the pair of twin edges to be removed and their incident faces.
  DHalfedge   *he1 = e;
  DHalfedge   *he2 = e->opposite();
  DFace       *f1 = he1->face();
  DFace       *f2 = he2->face();
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
      // The two halfedges form a "singleton" hole inside the incident face -
      // remove it. First notify the observers that we are about to remove
      // this hole.
      Face_handle               fh (f1);
      Ccb_halfedge_circulator   hole = (Halfedge_handle(he1))->ccb();

      _notify_before_remove_hole (fh,
                                  hole);

      // Erase the hole.
      if (! _find_and_erase_hole (f1, he1))
        _find_and_erase_hole (f1, he2);

      _notify_after_remove_hole (fh);

      // Delete the first end-vertex and its associated point.
      _notify_before_remove_vertex (Vertex_handle (he1->vertex()));

      _delete_point (he1->vertex()->point());
      dcel.delete_vertex (he1->vertex());

      _notify_after_remove_vertex ();

      // Delete the second end-vertex and its associated point.
      _notify_before_remove_vertex (Vertex_handle (he2->vertex()));

      _delete_point (he2->vertex()->point());
      dcel.delete_vertex (he2->vertex());

      _notify_after_remove_vertex ();


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
      if (he2->next() == he1)
      {
        he1 = he2;
        he2 = he1->opposite();
      }

      // Remove the two halfedges from the boundary chain by connecting
      // he1's predecessor with he2's successor.
      prev1 = he1->prev();

      prev1->set_next (he2->next());

      // In case the halfedges to be deleted are represantatives of a face
      // boundary or a hole inside the face, replace them by prev1.
      // Note that as we just change the hole representatives, we do not have
      // to notify the observers.
      if (f1->halfedge() == he1 || f1->halfedge() == he2)
      {
        f1->set_halfedge (prev1);

        // Go over all isolated vertices inside f1 and update their fictitious
        // incident face (which used to be he1 or he2) to be prev1.
        DIsolated_vertices_iter     iso_vert_it;

        for (iso_vert_it = f1->isolated_vertices_begin();
             iso_vert_it != f1->isolated_vertices_end(); ++iso_vert_it)

        {
          iso_vert_it->set_halfedge (prev1);
        }
      }
      else
      {
        if (_find_and_erase_hole (f1, he1))
          f1->add_hole (prev1);
        else if (_find_and_erase_hole (f1, he2))
          f1->add_hole (prev1);
      }

      // In case he2 is the representative halfegde of its target vertex,
      // replace it by prev1 (which also points at this vertex).
      if (he2->vertex()->halfedge() == he2)
          he2->vertex()->set_halfedge (prev1);

      // Delete the vertex that forms the tip of the "antenna"
      _notify_before_remove_vertex (Vertex_handle (he1->vertex()));

      _delete_point (he1->vertex()->point());
      dcel.delete_vertex (he1->vertex());

      _notify_after_remove_vertex();

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
    // Boundary of the incident face.
    if (_is_on_outer_boundary (f1, prev1) != NULL)
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

      // Create the new hole.
      f1->add_hole (he1->next());

      // As the outer boundary of f1 may be represented by any of the
      // halfedges in between he1 -> ... -> he2 (the halfedges in between
      // represent the outer boundary of the new hole that is formed),
      // We represent the outer boundary of f1 by prev1, which definately
      // stays on the outer boundary.
      f1->set_halfedge (prev1);

      // Go over all isolated vertices inside f1 and update their fictitious
      // incident face (which used to be he1 or he2) to be prev1.
      DIsolated_vertices_iter     iso_vert_it;

      for (iso_vert_it = f1->isolated_vertices_begin();
           iso_vert_it != f1->isolated_vertices_end(); ++iso_vert_it)
      {
        iso_vert_it->set_halfedge (prev1);
      }

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
      DHalfedge   *hole =  _is_on_inner_boundary (f1, prev1);

      CGAL_assertion (hole != NULL);

      // Notify the observers we are about to split a hole.
      _notify_before_split_hole (Face_handle (f1),
				 (Halfedge_handle (hole))->ccb(),
				 Halfedge_handle (he1));

      // We first remove the existing hole from the incident face.
      if (! _find_and_erase_hole (f1, hole))
      {
        CGAL_assertion (false);
        return (NULL);
      }

      // Make he1's and he2's predecessor the representatives of the two new
      // holes.
      f1->add_hole (prev1);
      f1->add_hole (prev2);

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
  const bool   he1_is_on_hole = (_is_on_inner_boundary (f1, he1) != NULL);
  const bool   he2_is_on_hole = (_is_on_inner_boundary (f2, he2) != NULL);
  DHalfedge   *ccb;

  prev1 = he1->prev();
  prev2 = he2->prev();

  CGAL_assertion (!he1_is_on_hole || !he2_is_on_hole);

  if (!he1_is_on_hole && !he2_is_on_hole)
  {
    // Both halfedges lie on the outer boundary of their incident faces.
    // We first set the incident face of f2's boundary halfedges to be f1.
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
    {
      f1->set_halfedge (prev1);

      // Go over all isolated vertices inside f1 and update their fictitious
      // incident face (which used to be he1) to be prev1.
      DIsolated_vertices_iter     iso_vert_it;

      for (iso_vert_it = f1->isolated_vertices_begin();
           iso_vert_it != f1->isolated_vertices_end(); ++iso_vert_it)
      {
        iso_vert_it->set_halfedge (prev1);
      }
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

    // Delete the face f2 and the pair of halfdges.
    dcel.delete_face (f2);
    dcel.delete_edge (he1);
      
    // Notify the observers that the faces have been merged.
    _notify_after_merge_face (Face_handle (f1));

    // Notify the observers that an edge has been deleted.
    _notify_after_remove_edge();

    // Return the merged face.
    return (f1); 
  }

  // In this case we merge a face with another face that now forms a hole
  // inside it. In this case we make sure that f1 contains the hole f2, so
  // we can merge f2 with it (we swap roles between the halfedges if
  // necessary).
  if (he2_is_on_hole)
  { 
    he1 = he2;
    he2 = he1->opposite();
    
    f1 = he1->face();
    f2 = he2->face();
    
    prev1 = he1->prev();
    prev2 = he2->prev();
  }

  // We wish to remove f2, so we set all the incident face of all
  // halfedges along its boundary to be f1.
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
  if (_find_and_erase_hole (f1, he1))
    f1->add_hole (prev1);
  
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

  // Delete the face f2 and the pair of halfdges.
  dcel.delete_face (f2);
  dcel.delete_edge (he1);
      
  // Notify the observers that the faces have been merged.
  _notify_after_merge_face (Face_handle (f1));

  // Notify the observers that an edge has been deleted.
  _notify_after_remove_edge();

  // Return the merged face.
  return (f1);
}

CGAL_END_NAMESPACE

#endif
