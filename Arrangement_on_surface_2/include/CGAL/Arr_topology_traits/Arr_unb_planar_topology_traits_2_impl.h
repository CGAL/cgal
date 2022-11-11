// Copyright (c) 2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s): Ron Wein  <wein@post.tau.ac.il>
//            Efi Fogel <efif@post.tau.ac.il>

#ifndef CGAL_ARR_UNB_PLANAR_TOPOLOGY_TRAITS_2_IMPL_H
#define CGAL_ARR_UNB_PLANAR_TOPOLOGY_TRAITS_2_IMPL_H

#include <CGAL/license/Arrangement_on_surface_2.h>

/*! \file
 *
 * Member-function definitions for the
 * Arr_unb_planar_topology_traits_2<GeomTraits> class.
 */

namespace CGAL {

//-----------------------------------------------------------------------------
// Default constructor.
//
template <typename GeomTraits, typename Dcel_>
Arr_unb_planar_topology_traits_2<GeomTraits, Dcel_>::
Arr_unb_planar_topology_traits_2():
  Base(),
  v_bl(nullptr),
  v_tl(nullptr),
  v_br(nullptr),
  v_tr(nullptr),
  n_inf_verts(0),
  fict_face(nullptr)
{}

//-----------------------------------------------------------------------------
// Constructor with a geometry-traits class.
//
template <typename GeomTraits, typename Dcel_>
Arr_unb_planar_topology_traits_2<GeomTraits, Dcel_>::
Arr_unb_planar_topology_traits_2 (const Geometry_traits_2 * geom_traits) :
  Base (geom_traits),
  v_bl (nullptr),
  v_tl (nullptr),
  v_br (nullptr),
  v_tr (nullptr),
  n_inf_verts (0),
  fict_face (nullptr)
{}

//-----------------------------------------------------------------------------
// Assign the contents of another topology-traits class.
//
template <typename GeomTraits, typename Dcel_>
void Arr_unb_planar_topology_traits_2<GeomTraits, Dcel_>::
assign(const Self& other)
{
  // Assign the base class.
  Base::assign (other);

  // Update the topology-traits properties after the DCEL have been updated.
  dcel_updated();

  return;
}

//-----------------------------------------------------------------------------
// Make the necessary updates after the DCEL structure have been updated.
//
template <typename GeomTraits, typename Dcel_>
void Arr_unb_planar_topology_traits_2<GeomTraits, Dcel_>::dcel_updated ()
{
  // Go over the DCEL vertices and locate the four fictitious ones.
  typename Dcel::Vertex_iterator       vit;
  Halfedge                            *first_he, *next_he;

  v_bl = v_tl = v_br = v_tr = nullptr;
  n_inf_verts = 0;
  for (vit = this->m_dcel.vertices_begin();
       vit != this->m_dcel.vertices_end(); ++vit)
  {
    if (! vit->has_null_point())
      continue;

    n_inf_verts++;

    // The current vertex is not associated with a point - check if it has
    // only two incident halfedges. If so, it is one of the four fictitious
    // vertices.
    first_he = vit->halfedge();
    next_he = first_he->next()->opposite();

    if (next_he->next()->opposite() == first_he)
    {
      Arr_parameter_space ps_x = vit->parameter_space_in_x();
      Arr_parameter_space ps_y = vit->parameter_space_in_y();

      if (ps_x == ARR_LEFT_BOUNDARY && ps_y == ARR_BOTTOM_BOUNDARY)
        v_bl = &(*vit);
      else if (ps_x == ARR_LEFT_BOUNDARY && ps_y == ARR_TOP_BOUNDARY)
        v_tl = &(*vit);
      else if (ps_x == ARR_RIGHT_BOUNDARY && ps_y == ARR_BOTTOM_BOUNDARY)
        v_br = &(*vit);
      else if (ps_x == ARR_RIGHT_BOUNDARY && ps_y == ARR_TOP_BOUNDARY)
        v_tr = &(*vit);
      else
        // We should never reach here:
        CGAL_error();
    }
  }
  CGAL_assertion(v_bl != nullptr && v_tl != nullptr && v_br != nullptr && v_tr != nullptr);

  // Go over the DCEL faces and locate the fictitious face.
  typename Dcel::Face_iterator         fit;

  fict_face = nullptr;
  for (fit = this->m_dcel.faces_begin();
       fit != this->m_dcel.faces_end(); ++fit)
  {
    if (fit->is_fictitious())
    {
      CGAL_assertion (fict_face == nullptr);

      fict_face = &(*fit);
    }
  }
  CGAL_assertion (fict_face != nullptr);

  return;
}

//-----------------------------------------------------------------------------
// Initialize an empty DCEL structure.
//
template <typename GeomTraits, typename Dcel_>
void Arr_unb_planar_topology_traits_2<GeomTraits, Dcel_>::init_dcel ()
{
  // Clear the current DCEL.
  this->m_dcel.delete_all();

  // Create the fictitious unbounded face.
  fict_face = this->m_dcel.new_face();

  fict_face->set_unbounded (true);
  fict_face->set_fictitious (true);

  // Create the four fictitious vertices corresponding to corners of the
  // bounding rectangle.
  v_bl = this->m_dcel.new_vertex();
  v_bl->set_boundary (ARR_LEFT_BOUNDARY, ARR_BOTTOM_BOUNDARY);

  v_tl = this->m_dcel.new_vertex();
  v_tl->set_boundary (ARR_LEFT_BOUNDARY, ARR_TOP_BOUNDARY);

  v_br = this->m_dcel.new_vertex();
  v_br->set_boundary (ARR_RIGHT_BOUNDARY, ARR_BOTTOM_BOUNDARY);

  v_tr = this->m_dcel.new_vertex();
  v_tr->set_boundary (ARR_RIGHT_BOUNDARY, ARR_TOP_BOUNDARY);

  // Create a four pairs of twin halfedges connecting the two vertices,
  // and link them together to form the bounding rectangle, forming a hole
  // in the fictitious face.
  //
  //                            he2
  //             v_tl (.) ----------------> (.) v_tr
  //                   ^ <------------------
  //                   ||                   ^|
  //  fict_face    he1 ||        in_f       ||
  //                   ||                   || he3
  //                   |V                   ||
  //                     ------------------> V
  //             v_bl (.) <---------------- (.) v_br
  //                             he4
  //
  Halfedge           *he1 = this->m_dcel.new_edge();
  Halfedge           *he1_t = he1->opposite();
  Halfedge           *he2 = this->m_dcel.new_edge();
  Halfedge           *he2_t = he2->opposite();
  Halfedge           *he3 = this->m_dcel.new_edge();
  Halfedge           *he3_t = he3->opposite();
  Halfedge           *he4 = this->m_dcel.new_edge();
  Halfedge           *he4_t = he4->opposite();
  Outer_ccb          *oc = this->m_dcel.new_outer_ccb();
  Inner_ccb          *ic = this->m_dcel.new_inner_ccb();
  Face               *in_f = this->m_dcel.new_face();

  he1->set_curve (nullptr);
  he2->set_curve (nullptr);
  he3->set_curve (nullptr);
  he4->set_curve (nullptr);

  he1->set_next (he2);        he1_t->set_next (he4_t);
  he2->set_next (he3);        he4_t->set_next (he3_t);
  he3->set_next (he4);        he3_t->set_next (he2_t);
  he4->set_next (he1);        he2_t->set_next (he1_t);

  he1->set_vertex (v_tl);     he1_t->set_vertex (v_bl);
  he2->set_vertex (v_tr);     he2_t->set_vertex (v_tl);
  he3->set_vertex (v_br);     he3_t->set_vertex (v_tr);
  he4->set_vertex (v_bl);     he4_t->set_vertex (v_br);

  oc->set_face (in_f);
  ic->set_face (fict_face);

  he1->set_inner_ccb (ic);       he1_t->set_outer_ccb (oc);
  he2->set_inner_ccb (ic);       he2_t->set_outer_ccb (oc);
  he3->set_inner_ccb (ic);       he3_t->set_outer_ccb (oc);
  he4->set_inner_ccb (ic);       he4_t->set_outer_ccb (oc);

  // Assign the incident halfedges of the two fictitious vertices.
  v_bl->set_halfedge (he1_t);
  v_tl->set_halfedge (he2_t);
  v_tr->set_halfedge (he3_t);
  v_br->set_halfedge (he4_t);

  // Set the direction of the halfedges:
  he1->set_direction (ARR_LEFT_TO_RIGHT);
  he2->set_direction (ARR_LEFT_TO_RIGHT);
  he3->set_direction (ARR_RIGHT_TO_LEFT);
  he4->set_direction (ARR_RIGHT_TO_LEFT);

  // Set the inner component of the fictitious face.
  fict_face->add_inner_ccb (ic, he1);

  // Set the real unbounded face, in the interior of the bounding rectangle.
  in_f->add_outer_ccb (oc, he1_t);
  in_f->set_unbounded (true);

  // Mark that there are four vertices at infinity (the fictitious ones)
  // in the arrangement.
  n_inf_verts = 4;

  return;
}

//-----------------------------------------------------------------------------
// Check if the given vertex is associated with the given curve end.
//
template <typename GeomTraits, typename Dcel_>
bool Arr_unb_planar_topology_traits_2<GeomTraits, Dcel_>::
are_equal(const Vertex *v,
          const X_monotone_curve_2& cv, Arr_curve_end ind,
          Arr_parameter_space ps_x, Arr_parameter_space ps_y) const
{
  // In case the given boundary conditions do not match those of the given
  // vertex, v cannot represent the curve end.
  if (ps_x != v->parameter_space_in_x() || ps_y != v->parameter_space_in_y())
    return false;

  // Compare the curve end with the vertex.
  Comparison_result     res;

  if (ps_x != ARR_INTERIOR) {
    // The curve end lies at x = +/- oo and so does v. Check if the curve
    // overlaps with the curve that currently induces v.
    Arr_curve_end                  v_ind;
    const X_monotone_curve_2  *v_cv = _curve (v, v_ind);

    if (v_cv == nullptr)
      return (v->parameter_space_in_x() == ps_x &&
              v->parameter_space_in_y() == ps_y);

    res = this->m_geom_traits->compare_y_curve_ends_2_object()(cv,
                                                               *v_cv, v_ind);
  }
  else {
    CGAL_assertion (ps_y != ARR_INTERIOR);

    // The curve end lies at y = +/- oo and so does v. Check if the curve
    // overlaps with the curve that currently induces v.
    Arr_curve_end                  v_ind;
    const X_monotone_curve_2  *v_cv = _curve (v, v_ind);

    if (v_cv == nullptr)
      return (v->parameter_space_in_x() == ARR_INTERIOR &&
              v->parameter_space_in_y() == ps_y);

    res =
      this->m_geom_traits->compare_x_curve_ends_2_object() (cv, ind,
                                                            *v_cv, v_ind);
  }

  return (res == EQUAL);
}

//-----------------------------------------------------------------------------
// Given a curve end with boundary conditions and a face that contains the
// interior of the curve, find a place for a boundary vertex that will
// represent the curve end along the face boundary.
//
template <typename GeomTraits, typename Dcel_>
boost::optional
  <boost::variant
    <typename Arr_unb_planar_topology_traits_2<GeomTraits, Dcel_>::Vertex*,
     typename Arr_unb_planar_topology_traits_2<GeomTraits, Dcel_>::Halfedge*> >
Arr_unb_planar_topology_traits_2<GeomTraits, Dcel_>::
place_boundary_vertex(Face* f,
                      const X_monotone_curve_2& cv, Arr_curve_end ind,
                      Arr_parameter_space ps_x, Arr_parameter_space ps_y)
{
  typedef boost::variant<Vertex*, Halfedge*>    Non_optional_result;
  typedef boost::optional<Non_optional_result>  Result;

  // Get a halfedge on the outer CCB of f and start traversing the CCB.
  Halfedge* first = *(f->outer_ccbs_begin());
  Halfedge* curr = first;
  bool eq_source, eq_target;

  do {
    // Note we consider only fictitious halfedges and check whether they
    // contain the relevant curve end.
    if (curr->has_null_curve() &&
        _is_on_fictitious_edge(cv, ind, ps_x, ps_y, curr, eq_source, eq_target))
    {
      CGAL_assertion(! eq_source && ! eq_target);
      return Result(curr);
    }

    // Move to the next halfegde along the CCB.
    curr = curr->next();

  } while (curr != first);

  // If we reached here, we did not find a suitable halfegde, which should
  // never happen.
  CGAL_error();
  return boost::none;
}

//-----------------------------------------------------------------------------
// Locate a DCEL feature that contains the given unbounded curve end.
//
template <typename GeomTraits, typename Dcel_>
boost::variant
<typename Arr_unb_planar_topology_traits_2<GeomTraits, Dcel_>::Vertex*,
 typename Arr_unb_planar_topology_traits_2<GeomTraits, Dcel_>::Halfedge*,
 typename Arr_unb_planar_topology_traits_2<GeomTraits, Dcel_>::Face*>
Arr_unb_planar_topology_traits_2<GeomTraits, Dcel_>::
locate_curve_end (const X_monotone_curve_2& cv, Arr_curve_end ind,
                  Arr_parameter_space ps_x, Arr_parameter_space ps_y)
{
  typedef boost::variant<Vertex*, Halfedge*, Face*>     Result;

  // Start traversing the inner CCB of the fictitious face and try to locate
  // a feature that contains the curve end.
  Halfedge* first = *(fict_face->inner_ccbs_begin());
  Halfedge* curr = first;
  bool eq_source, eq_target;

  do {
    if (_is_on_fictitious_edge(cv, ind, ps_x, ps_y, curr, eq_source, eq_target))
    {
      if (eq_source) {
        // cv's end coincides with the source vertex of the current
        // fictitious halfedge. This means that cv overlaps the curve that
        // is associated with the only non-fictitious halfedge incident to
        // this vertex. We therefore return a pointer to this halfedge.
        Halfedge* he = curr->opposite()->next();
        CGAL_assertion(! he->has_null_curve());
        return Result(he);
      }
      else if (eq_target) {
        // cv's end coincides with the target vertex of the current
        // fictitious halfedge. This means that cv overlaps the curve that
        // is associated with the only non-fictitious halfedge incident to
        // this vertex. We therefore return a pointer to this halfedge.
        Halfedge* he = curr->opposite()->prev();
        CGAL_assertion(! he->has_null_curve());
        return Result(he);
      }

      // The current ficitious edge contains cv's end in its interior.
      // Note we use curr's twin, whose incident face is a valid
      // unbounded face (whereas the incident face of curr is the fictitious
      // face).
      Face* uf = curr->opposite()->outer_ccb()->face();
      CGAL_assertion (uf->is_unbounded() && ! uf->is_fictitious());
      return Result(uf);
    }

    curr = curr->next();

  } while (curr != first);

  // We should never reach here.
  CGAL_error();
  Vertex* v(nullptr);
  return Result(v);
}

//-----------------------------------------------------------------------------
// Split a fictitious edge using the given vertex.
//
template <typename GeomTraits, typename Dcel_>
typename Arr_unb_planar_topology_traits_2<GeomTraits, Dcel_>::Halfedge*
Arr_unb_planar_topology_traits_2<GeomTraits, Dcel_>::
split_fictitious_edge (Halfedge *e, Vertex *v)
{
  CGAL_precondition (v->parameter_space_in_x() != ARR_INTERIOR ||
                     v->parameter_space_in_y() != ARR_INTERIOR);

  // Increment the number of vertices at infinity.
  n_inf_verts++;

  // Get the split halfedge and its twin, and their incident faces.
  // Note that he1 lies on an outer boundary of an unbounded face, while
  // its twin he2 should lie on a hole (inner boundary) inside the fictitious
  // face.
  Halfedge       *he1 = e;
  Halfedge       *he2 = he1->opposite();

  CGAL_assertion (! he1->is_on_inner_ccb());
  Outer_ccb      *oc1 = he1->outer_ccb();

  CGAL_assertion (oc1->face()->is_unbounded());

  CGAL_assertion (he2->is_on_inner_ccb());
  Inner_ccb      *ic2 = he2->inner_ccb();

  CGAL_assertion (ic2->face() == fict_face);

  // Allocate a pair of new halfedges.
  Halfedge       *he3 = this->m_dcel.new_edge();
  Halfedge       *he4 = he3->opposite();

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
  he3->set_outer_ccb (oc1);
  he3->set_vertex (he1->vertex());

  he4->set_vertex (v);
  he4->set_next (he2);

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

  // Return a pointer to one of the existing halfedge that is incident to the
  // split vertex.
  return (he1);
}

//-----------------------------------------------------------------------------
// Determine whether the given face is unbounded.
//
template <typename GeomTraits, typename Dcel_>
bool Arr_unb_planar_topology_traits_2<GeomTraits, Dcel_>::
is_unbounded(const Face *f) const
{
  // Go over the outer CBB of the given face and look for fictitious halfedges.
  const Halfedge   *first = *(f->outer_ccbs_begin());
  const Halfedge   *curr = first;

  do
  {
    if (curr->has_null_curve())
      // Found a fictitious halfedge along the boundary: f is unbounded.
      return (true);

    curr = curr->next();

  } while (curr != first);

  // If we reached here, all halfedges along the face boundary are valid,
  // thus the face is bounded.
  return (false);
}

//-----------------------------------------------------------------------------
// Determine whether the given boundary vertex is redundant.
//
template <typename GeomTraits, typename Dcel_>
bool Arr_unb_planar_topology_traits_2<GeomTraits, Dcel_>::
is_redundant(const Vertex *v) const
{
  CGAL_precondition (v != v_bl && v != v_tl && v != v_br && v != v_tr);

  // A boundary vertex is redundant if there it is of degree 2 and (there
  // is no valid edge incident to it).
  const Halfedge  *first_he = v->halfedge();
  const Halfedge  *next_he = first_he->next()->opposite();

  if (next_he->next()->opposite() == first_he)
  {
    CGAL_assertion (first_he->has_null_curve() && next_he->has_null_curve());
    return (true);
  }

  return (false);
}

//-----------------------------------------------------------------------------
// Erase the given redundant vertex.
//
template <typename GeomTraits, typename Dcel_>
typename Arr_unb_planar_topology_traits_2<GeomTraits, Dcel_>::Halfedge*
Arr_unb_planar_topology_traits_2<GeomTraits, Dcel_>::
erase_redundant_vertex (Vertex *v)
{
  CGAL_precondition (is_redundant (v));

  // Assign pointers to the halfedges incident to v (make sure that there
  // are exactly teo pairs of fictitious halfedges), such that we have:
  //
  //            he1      he3
  //         -------> ------->
  //       (.)      (.)v     (.)
  //         <------- <-------
  //            he2      he4
  //
  Halfedge   *he1 = v->halfedge();
  Halfedge   *he2 = he1->opposite();
  Halfedge   *he3 = he1->next();
  Halfedge   *he4 = he3->opposite();

  CGAL_assertion (he1->has_null_curve() && he3->has_null_curve() &&
                  he4->next() == he2);

  // Keep pointers to the components that contain two halfedges he3 and he2,
  // pointing at the end vertices of the merged halfedge.
  Inner_ccb   *ic1 = (he3->is_on_inner_ccb()) ? he3->inner_ccb() : nullptr;
  Outer_ccb   *oc1 = (ic1 == nullptr) ? he3->outer_ccb() : nullptr;
  Inner_ccb   *ic2 = (he4->is_on_inner_ccb()) ? he4->inner_ccb() : nullptr;
  Outer_ccb   *oc2 = (ic2 == nullptr) ? he4->outer_ccb() : nullptr;

  // As he1 and he2 will evetually represent the merged edge, while he3 and he4
  // will be deleted, check if the deleted halfedges are represantatives of a
  // face boundary or a hole inside these faces. If so, replace he3 by he1 and
  // he4 by he2.
  if (ic1 != nullptr && ic1->halfedge() == he3)
    ic1->set_halfedge (he1);
  else if (oc1 != nullptr && oc1->halfedge() == he3)
    oc1->set_halfedge (he1);

  if (ic2 != nullptr && ic2->halfedge() == he4)
    ic2->set_halfedge (he2);
  else if (oc2 != nullptr && oc2->halfedge() == he4)
    oc2->set_halfedge (he2);

  // If he3 is the incident halfedge to its target, replace it by he1.
  if (he3->vertex()->halfedge() == he3)
    he3->vertex()->set_halfedge (he1);

  // Disconnect he3 and he4 from the edge list.
  CGAL_assertion (he3->next() != he4);

  he1->set_next (he3->next());
  he4->prev()->set_next (he2);

  // Set the properties of the merged edge.
  he1->set_vertex (he3->vertex());

  // Decrement the number of vertices at infinity (note we do not actually
  // free the vertex - the Arrangement_on_surface_2 class will do it).
  n_inf_verts--;

  // Delete the redundant halfedge pair.
  this->m_dcel.delete_edge (he3);

  return (he1);
}

//-----------------------------------------------------------------------------
// Compare the x-coordinates of a given vertex (which may lie at infinity) and
// the given point.
//
template <typename GeomTraits, typename Dcel_>
Comparison_result
Arr_unb_planar_topology_traits_2<GeomTraits, Dcel_>::
compare_x (const Point_2& p, const Vertex* v) const
{
  // First check if the vertex v lies at x = -oo (then it is obviously smaller
  // than p), or at x = +oo (then it is obviously larger).
  const Arr_parameter_space ps_x = v->parameter_space_in_x();

  if (ps_x == ARR_LEFT_BOUNDARY)
    return (LARGER);
  else if (ps_x == ARR_RIGHT_BOUNDARY)
    return (SMALLER);

  // Check if the vertex lies at y = +/- oo.
  const Arr_parameter_space ps_y = v->parameter_space_in_y();

  if (ps_y != ARR_INTERIOR)
  {
    // Compare the x-position of the vertical asymptote of the curve incident
    // to v with the x-coodinate of p.
    Arr_curve_end v_ind = ARR_MIN_END;
    const X_monotone_curve_2* v_cv = _curve (v, v_ind);

    CGAL_assertion(v_cv != nullptr);
    return
      (this->m_geom_traits->compare_x_point_curve_end_2_object()(p, *v_cv,
                                                                 v_ind));
  }

  // In this case v represents a normal point, and we compare it with p.
  return (this->m_geom_traits->compare_x_2_object() (p, v->point()));
}

//-----------------------------------------------------------------------------
// Compare the given vertex (which may lie at infinity) and the given point.
//
template <typename GeomTraits, typename Dcel_>
Comparison_result
Arr_unb_planar_topology_traits_2<GeomTraits, Dcel_>::
compare_xy (const Point_2& p, const Vertex* v) const
{
  // First check if the vertex v lies at x = -oo (then it is obviously smaller
  // than p), or at x = +oo (then it is obviously larger).
  const Arr_parameter_space ps_x = v->parameter_space_in_x();

  if (ps_x == ARR_LEFT_BOUNDARY) return (LARGER);
  else if (ps_x == ARR_RIGHT_BOUNDARY) return (SMALLER);

  // Check if the vertex lies at y = +/- oo.
  const Arr_parameter_space ps_y = v->parameter_space_in_y();

  if (ps_y != ARR_INTERIOR) {
    // Compare the x-position of the vertical asymptote of the curve incident
    // to v with the x-coodinate of p.
    Arr_curve_end v_ind = ARR_MIN_END;
    const X_monotone_curve_2* v_cv = _curve (v, v_ind);

    CGAL_assertion (v_cv != nullptr);

    Comparison_result res =
      this->m_geom_traits->compare_x_point_curve_end_2_object() (p, *v_cv,
                                                                 v_ind);

    if (res != EQUAL) return (res);

    // In case of equality, consider whether v lies at y = -oo or at y = +oo.
    return (ps_y == ARR_BOTTOM_BOUNDARY) ? LARGER : SMALLER;
  }

  // In this case v represents a normal point, and we compare it with p.
  return (this->m_geom_traits->compare_xy_2_object()(p, v->point()));
}

//-----------------------------------------------------------------------------
// Compare the relative y-position of the given point and the given edge
// (which may be fictitious).
//
template <typename GeomTraits, typename Dcel_>
Comparison_result
Arr_unb_planar_topology_traits_2<GeomTraits, Dcel_>::
compare_y_at_x(const Point_2& p, const Halfedge* he) const
{
  // In case of a valid edge, just compare p to its associated curve.
  if (! he->has_null_curve())
    return (this->m_geom_traits->compare_y_at_x_2_object()(p, he->curve()));

  // Otherwise, determine on which edge of the bounding rectangle does he lie.
  // Note this can be either the top edge or the bottom edge (and not the
  // left or the right edge), as p must lie in its x-range.
  CGAL_assertion((he->vertex()->parameter_space_in_x() == ARR_INTERIOR) ||
                 (he->vertex()->parameter_space_in_x() !=
                  he->opposite()->vertex()->parameter_space_in_x()));
  CGAL_assertion((he->vertex()->parameter_space_in_y() != ARR_INTERIOR) &&
                 (he->vertex()->parameter_space_in_y() ==
                  he->opposite()->vertex()->parameter_space_in_y()));

  // he lies on the bottom edge, so p is obviously above it.
  if (he->vertex()->parameter_space_in_y() == ARR_BOTTOM_BOUNDARY)
    return LARGER;
  // he lies on the top edge, so p is obviously below it.
  else
    return SMALLER;
}

//-----------------------------------------------------------------------------
// Get the curve associated with a boundary vertex.
//
template <typename GeomTraits, typename Dcel_>
const typename
Arr_unb_planar_topology_traits_2<GeomTraits, Dcel_>::X_monotone_curve_2*
Arr_unb_planar_topology_traits_2<GeomTraits, Dcel_>::
_curve (const Vertex* v, Arr_curve_end& ind) const
{
  // Go over the incident halfedges of v until encountering the halfedge
  // associated with a valid curve (v should have three incident halfedges,
  // two of the are fictitious and one associated with a curve).
  const Halfedge* he = v->halfedge();

  while (he->has_null_curve()) {
    he = he->next()->opposite();

    // No incident curve were found:
    if (he == v->halfedge()) return (nullptr);
  }

  // The halfedge he is directed toward v, so if it is directed from left to
  // right, v represents the maximal end of cv, otherwise it represents its
  // minimal end.
  ind = (he->direction() == ARR_LEFT_TO_RIGHT) ? ARR_MAX_END : ARR_MIN_END;

  // Return the x-monotone curve.
  return &(he->curve());
}

//-----------------------------------------------------------------------------
// Check whether the given infinite curve end lies on the given fictitious
// halfedge.
//
template <typename GeomTraits, typename Dcel_>
bool
Arr_unb_planar_topology_traits_2<GeomTraits, Dcel_>::
_is_on_fictitious_edge(const X_monotone_curve_2& cv, Arr_curve_end ind,
                       Arr_parameter_space ps_x, Arr_parameter_space ps_y,
                       const Halfedge* he,
                       bool& eq_source, bool& eq_target)
{
  eq_source = false;
  eq_target = false;

  // Get the end-vertices of the edge.
  const Vertex* v1 = he->opposite()->vertex();
  const Vertex* v2 = he->vertex();
  Comparison_result res1, res2;
  Arr_curve_end v_ind = ARR_MIN_END;

  // Check if this is a "vertical" ficitious edge.
  Arr_parameter_space he_ps_x = v1->parameter_space_in_x();
  if ((he_ps_x != ARR_INTERIOR) && (he_ps_x == v2->parameter_space_in_x())) {
    // If the edge lies on x = +/- oo, the curve endpoint must also lie there.
    CGAL_assertion((he_ps_x == ARR_LEFT_BOUNDARY) ||
                   (he_ps_x == ARR_RIGHT_BOUNDARY));

    if (he_ps_x != ps_x) return false;

    // Compare the y-position of the curve end to the source vertex.
    if ((v1 == v_bl) || (v1 == v_br)) {
      // These vertices are below any curve.
      res1 = LARGER;
    }
    else if ((v1 == v_tl) || (v1 == v_tr)) {
      // These vertices are above any curve.
      res1 = SMALLER;
    }
    else {
      const Arr_curve_end ind =
        (ps_x == ARR_LEFT_BOUNDARY) ? ARR_MIN_END : ARR_MAX_END;

      res1 =
        this->m_geom_traits->compare_y_curve_ends_2_object()(cv,
                                                             *_curve (v1, v_ind),
                                                             ind);
      if (res1 == EQUAL) {
        eq_source = true;
        return true;
      }
    }

    // Compare the y-position of the curve end to the target vertex.
    if ((v2 == v_bl) || (v2 == v_br)) {
      // These vertices are below any curve.
      res2 = LARGER;
    }
    else if ((v2 == v_tl) || (v2 == v_tr)) {
      // These vertices are above any curve.
      res2 = SMALLER;
    }
    else {
      const Arr_curve_end ind =
        (ps_x == ARR_LEFT_BOUNDARY) ? ARR_MIN_END : ARR_MAX_END;

      res2 =
        this->m_geom_traits->compare_y_curve_ends_2_object()(cv,
                                                             *_curve (v2, v_ind),
                                                             ind);

      if (res2 == EQUAL) {
        eq_target = true;
        return true;
      }
    }
  }
  else {
    // If we reched here, we have a "horizontal" fictitious halfedge.
    Arr_parameter_space he_ps_y = v1->parameter_space_in_y();

    CGAL_assertion((he_ps_y == ARR_BOTTOM_BOUNDARY ||
                    he_ps_y == ARR_TOP_BOUNDARY) &&
                   he_ps_y == v2->parameter_space_in_y());

    // If the edge lies on y = +/- oo, the curve endpoint must also lie there
    // (and must not lies at x = +/- oo.
    if ((ps_x != ARR_INTERIOR) || (he_ps_y != ps_y)) return false;

    // Compare the x-position of the curve end to the source vertex.
    if ((v1 == v_bl) || (v1 == v_tl)) {
      // These vertices are to the left of any curve.
      res1 = LARGER;
    }
    else if ((v1 == v_br) || (v1 == v_tr)) {
      // These vertices are to the right of any curve.
      res1 = SMALLER;
    }
    else {
      const X_monotone_curve_2  *v_cv1 = _curve (v1, v_ind);

      // Note that v1 is a non-fictitious vertex, therefore we expect it to
      // be associated with a valid curve end. If this is not the case, we
      // can assume that v1 has been created in order to represent the other
      // curve-end of cv, the curve that is currently being insrted into the
      // arrangement, but it hasn't been associated with a valid halfedge
      // yet, as the insertion process is still ongoing.
      // The comparison result in this case is trivial.
      if (v_cv1 != nullptr) {
        res1 =
          this->m_geom_traits->compare_x_curve_ends_2_object()(cv, ind,
                                                               *v_cv1, v_ind);

        if (res1 == EQUAL) {
          eq_source = true;
          return true;
        }
      }
      else {
        res1 = (ind == ARR_MIN_END) ? SMALLER : LARGER;
      }
    }

    // Compare the x-position of the curve end to the target vertex.
    if ((v2 == v_bl) || (v2 == v_tl)) {
      // These vertices are to the left of any curve.
      res2 = LARGER;
    }
    else if ((v2 == v_br) || (v2 == v_tr)) {
      // These vertices are to the right of any curve.
      res2 = SMALLER;
    }
    else {
      const X_monotone_curve_2* v_cv2 = _curve(v2, v_ind);

      // Note that v2 is a non-fictitious vertex, therefore we expect it to
      // be associated with a valid curve end. If this is not the case, we
      // can assume that v2 has been created in order to represent the other
      // curve-end of cv, the curve that is currently being insrted into the
      // arrangement, but it hasn't been associated with a valid halfedge
      // yet, as the insertion process is still ongoing.
      // The comparison result in this case is trivial.
      if (v_cv2 != nullptr) {
        res2 =
          this->m_geom_traits->compare_x_curve_ends_2_object()(cv, ind,
                                                               *v_cv2, v_ind);

        if (res2 == EQUAL) {
          eq_target = true;
          return true;
        }
      }
      else {
        res2 = (ind == ARR_MIN_END) ? SMALLER : LARGER;
      }
    }
  }

  return (res1 != res2);
}

} //namespace CGAL

#endif
