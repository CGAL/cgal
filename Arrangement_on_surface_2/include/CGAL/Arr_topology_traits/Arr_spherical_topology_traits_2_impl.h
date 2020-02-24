// Copyright (c) 2006,2007,2008,2009,2010,2011 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Efi Fogel         <efif@post.tau.ac.il>
//                 Ron Wein          <wein@post.tau.ac.il>

#ifndef CGAL_ARR_SPHERICAL_TOPOLOGY_TRAITS_2_IMPL_H
#define CGAL_ARR_SPHERICAL_TOPOLOGY_TRAITS_2_IMPL_H

#include <CGAL/license/Arrangement_on_surface_2.h>

/*! \file
 *
 * Member-function definitions for the
 * Arr_spherical_topology_traits_2<GeomTraits> class.
 */

namespace CGAL {

/*! \brief constructs default */
template <typename GeomTraits, typename Dcel>
Arr_spherical_topology_traits_2<GeomTraits, Dcel>::
Arr_spherical_topology_traits_2() :
  m_spherical_face(nullptr),
  m_north_pole(nullptr),
  m_south_pole(nullptr),
  m_own_geom_traits(true)
{
  m_geom_traits = new Gt_adaptor_2;
  m_boundary_vertices = Vertex_map(Vertex_key_comparer(m_geom_traits));
}

/*! \brief constructs from a geometry-traits object. */
template <typename GeomTraits, typename Dcel>
Arr_spherical_topology_traits_2<GeomTraits, Dcel>::
Arr_spherical_topology_traits_2(const Geometry_traits_2* traits) :
  m_spherical_face(nullptr),
  m_north_pole(nullptr),
  m_south_pole(nullptr),
  m_own_geom_traits(false)
{
  m_geom_traits = static_cast<const Gt_adaptor_2*>(traits);
  m_boundary_vertices = Vertex_map(Vertex_key_comparer(m_geom_traits));
}

/*! \brief destructs */
template <typename GeomTraits, typename Dcel>
Arr_spherical_topology_traits_2<GeomTraits, Dcel>::
~Arr_spherical_topology_traits_2()
{
  // Clear the DCEL.
  m_dcel.delete_all();

  if (m_own_geom_traits && (m_geom_traits != nullptr)) {
    delete m_geom_traits;
    m_geom_traits = nullptr;
  }
}

/*! \brief assigns the contents of another topology-traits class */
template <typename GeomTraits, typename Dcel>
void Arr_spherical_topology_traits_2<GeomTraits, Dcel>::
assign(const Self& other)
{
  // Clear the current DCEL and duplicate the other DCEL.
  m_dcel.delete_all();
  m_dcel.assign(other.m_dcel);

  // Take care of the traits object.
  if (m_own_geom_traits && m_geom_traits != nullptr) {
    delete m_geom_traits;
    m_geom_traits = nullptr;
  }

  if (other.m_own_geom_traits) {
    m_geom_traits = new Gt_adaptor_2;
    m_own_geom_traits = true;
  }
  else {
    m_geom_traits = other.m_geom_traits;
    m_own_geom_traits = false;
  }

  // Update the rest of the properties.
  dcel_updated();
}

/*! \brief initializes an empty DCEL structure. */
template <typename GeomTraits_, typename Dcel_>
void Arr_spherical_topology_traits_2<GeomTraits_, Dcel_>::dcel_updated()
{
  typedef Dcel_                                         Dcel;

  // Go over the DCEL vertices and locate the south and north pole (if any)
  // and any other vertex on the line of discontinuity.

  m_north_pole = nullptr;
  m_south_pole = nullptr;
  m_boundary_vertices.clear();

  typename Dcel::Vertex_iterator vit = this->m_dcel.vertices_begin();
  for (; vit != this->m_dcel.vertices_end(); ++vit) {
    Arr_parameter_space bx = vit->parameter_space_in_x();
    Arr_parameter_space by = vit->parameter_space_in_y();

    if (by == ARR_BOTTOM_BOUNDARY) m_south_pole = &(*vit);
    else if (by == ARR_TOP_BOUNDARY) m_north_pole = &(*vit);
    else if (bx != ARR_INTERIOR) {
      const Point_2& key = vit->point();
      m_boundary_vertices.insert(Vertex_value(key, &(*vit)));
    }
  }

  // Go over the DCEL faces and locate the spherical face, which is the only
  // face with no outer CCB.

  m_spherical_face = nullptr;
  typename Dcel::Face_iterator fit = this->m_dcel.faces_begin();
  for (; fit != this->m_dcel.faces_end(); ++fit) {
    if (fit->number_of_outer_ccbs() == 0) {
      CGAL_assertion(m_spherical_face == nullptr);

      m_spherical_face = &(*fit);
      break;
    }
  }
  CGAL_assertion(m_spherical_face != nullptr);
}

/*! \brief initializes an empty DCEL structure. */
template <typename GeomTraits, typename Dcel>
void Arr_spherical_topology_traits_2<GeomTraits, Dcel>::init_dcel()
{
  // std::cout << "init_dcel()" << std::endl;
  // Clear the current DCEL.
  m_dcel.delete_all();
  m_boundary_vertices.clear();

  // Create the face.
  m_spherical_face = this->m_dcel.new_face();
  m_spherical_face->set_unbounded(false);
  m_spherical_face->set_fictitious(false);

  m_north_pole = nullptr;
  m_south_pole = nullptr;
}

/*! \brief determines whether a point lies in the interior of a given face. */
template <typename GeomTraits, typename Dcel>
bool Arr_spherical_topology_traits_2<GeomTraits, Dcel>::
is_in_face(const Face* f, const Point_2& p, const Vertex* v) const
{
  // std::cout << "is_in_face()" << std::endl;
  CGAL_precondition((v == nullptr) || !v->has_null_point());
  CGAL_precondition((v == nullptr) ||
                    m_geom_traits->equal_2_object()(p, v->point()));

  /* There is always one face that contains everything else. It has no
   * outer CCB's. When a new face is constructed, we make sure that the
   * face that contains everything also contains the north pole. (In the
   * degenerate case, where a vertex coincides with the north pole, the face
   * that contains everything is incident to the north pole.)
   * If the face has no outer ccb's, it contains everything:
   */
#if 0
  std::cout << "p: " << p
            << ", f->number_of_outer_ccbs(): " << f->number_of_outer_ccbs()
            << std::endl;
#endif
  if (f->number_of_outer_ccbs() == 0) return true;
  if (((v != nullptr) && (v->parameter_space_in_y() == ARR_TOP_BOUNDARY)) ||
      (m_geom_traits->parameter_space_in_y_2_object()(p) == ARR_TOP_BOUNDARY))
    return false;

  /*! \todo a temporary test
   * if (((v != nullptr) && (v->parameter_space_in_y() == ARR_BOTTOM_BOUNDARY)) ||
   *   (p.is_min_boundary()))
   * return false;
   */

  typename Gt_adaptor_2::Parameter_space_in_x_2 ps_x_op =
    m_geom_traits->parameter_space_in_x_2_object();
  typename Gt_adaptor_2::Parameter_space_in_y_2 ps_y_op =
    m_geom_traits->parameter_space_in_y_2_object();
  typename Gt_adaptor_2::Compare_x_2 cmp_x_op =
    m_geom_traits->compare_x_2_object();
  typename Gt_adaptor_2::Compare_y_at_x_2 cmp_y_at_x_op =
    m_geom_traits->compare_y_at_x_2_object();
  typename Gt_adaptor_2::Compare_x_point_curve_end_2 cmp_x_pt_ce =
    m_geom_traits->compare_x_point_curve_end_2_object();

  // Process the input point.
  bool p_is_interior_x = !(m_geom_traits->is_on_y_identification_2_object()(p));

  /* Maintain a counter of the number of x-monotone curves that intersect an
   * upward vertical ray emanating from p. Handle degenerate cases as
   * explained below).
   */
  unsigned int num_intersections = 0;

  /* Traverse all outer CCBs of the face. For each boundary component go over
   * all its halfedges, and count those which are above p.
   */
  typename Face::Outer_ccb_const_iterator oit;
  for (oit = f->outer_ccbs_begin(); oit != f->outer_ccbs_end(); ++oit) {
    const Halfedge* first = *oit;
    const Halfedge* curr = first;

    /* Compare p to the source vertex of the first halfedge. If p coincides
     * with this vertex, p is obviously not in the interior of the face.
     */
    if (curr->opposite()->vertex() == v) return false;

    /*! We identify 2 main cases:
     * 1. The vertical ray intersects the boundary at a halfedge. In this
     * case the x-possition of p is strictly larger than the x-possition of
     * the current-curve source, and strictly smaller than x-possition of
     * the current-curve target, or vise versa.
     * 2. The vertical ray intersects the boundary at a vertex. In this case:
     * a. the x-possition of p is strictly smaller than the x-position of the
     * current-curve source, and equal to the x-position of the current-curve
     * target, and
     * b. the x-possition of p is equal to the x-position of the next-curve
     * source (not counting vertical curves in between), and strictly larger
     * than the x-possition of the next-curve target, or vise verase (that is,
     * the "smaller" and "larger" interchanged).
     */

    /* Indicates that a change between the x-position of p and the x-position
     * of the current-curve source, and the x-position of p and the x-position
     * of the current-curve target is pending. Used to handle case (2) above.
     */
    bool change_pending = false;

    /*! Indicates that the conditions described in (2.b) above are met during
     * the 1st iteration in the loop, which implies that the last curve must be
     * checked.
     */
    bool last_pending = false;

    Comparison_result res_pending = EQUAL, res_last = EQUAL,
      res_source = EQUAL, res_target;
    Arr_parameter_space ps_x_pending = ARR_INTERIOR, ps_x_last = ARR_INTERIOR,
      ps_x_source, ps_x_target = ARR_INTERIOR,
      ps_y_source, ps_y_target;

    do {
      /* Compare p to the target vertex of the current halfedge. If the
       * vertex v is on the boundary of the component, p is not in the interior
       * the face.
       */
      if (curr->vertex() == v) return false;

      // Ignore vertical curves:
      bool is_vertical = m_geom_traits->is_vertical_2_object()(curr->curve());
      if (is_vertical) {
        /* If this outer ccb chain contains the north pole, and our point
         * lies horizontaly between the two vertical curves that meet at
         * the north pole, increase the intersection counter
         */
        if (curr->direction() == ARR_LEFT_TO_RIGHT) {
          Arr_parameter_space ps_y_1 = ps_y_op(curr->curve(), ARR_MAX_END);
          Arr_parameter_space ps_y_2 = ps_y_op(curr->next()->curve(),
                                               ARR_MAX_END);
          if ((ps_y_1 == ARR_TOP_BOUNDARY) && (ps_y_2 == ARR_TOP_BOUNDARY)) {
            // Compare the x-coordinates:
            Comparison_result rc1 =
              cmp_x_pt_ce(p, curr->curve(), ARR_MAX_END);
            Comparison_result rc2 =
              cmp_x_pt_ce(p, curr->next()->curve(), ARR_MAX_END);
            if (rc1 == opposite(rc2)) ++num_intersections;
          }
        }
        curr = curr->next();
        continue;
      }

      /* If the current halfedge belongs to an "antenna". Namely, its
       * incident face is the same as its twin's, skip it to avoid counting
       * it twice.
       */
      const Face* curr_face = (curr->is_on_inner_ccb()) ?
        curr->inner_ccb()->face() : curr->outer_ccb()->face();
      const Halfedge* opp_he = curr->opposite();
      const Face* opp_curr_face = (opp_he->is_on_inner_ccb()) ?
        opp_he->inner_ccb()->face() : opp_he->outer_ccb()->face();

      if (curr_face == opp_curr_face) {
        curr = curr->next();
        continue;
      }

      Arr_curve_end ind_source, ind_target;
      if (curr->direction() == ARR_LEFT_TO_RIGHT) {
        ind_source = ARR_MIN_END;
        ind_target = ARR_MAX_END;
      }
      else {
        ind_source = ARR_MAX_END;
        ind_target = ARR_MIN_END;
      }

      ps_x_source = ps_x_op(curr->curve(), ind_source);
      ps_x_target = ps_x_op(curr->curve(), ind_target);

      ps_y_source = ps_y_op(curr->curve(), ind_source);
      ps_y_target = ps_y_op(curr->curve(), ind_target);

      if (!p_is_interior_x) {
        if (ps_x_source == ps_x_target) {
          curr = curr->next();
          continue;
        }

        if (ps_x_target != ARR_INTERIOR) {
          change_pending = true;
          ps_x_pending = (ps_x_target == ARR_LEFT_BOUNDARY) ?
            ARR_RIGHT_BOUNDARY : ARR_LEFT_BOUNDARY;
        }
        if (ps_x_source != ARR_INTERIOR) {
          if (change_pending) {
            change_pending = false;
            if (ps_x_pending == ps_x_source) {
              Comparison_result res_y_at_x = cmp_y_at_x_op(p, curr->curve());
              if (res_y_at_x == EQUAL) return false;
              if (res_y_at_x == SMALLER) num_intersections++;
            }
          } else {
            // This must be the first curve. Remember to check the last curve
            ps_x_last = (ps_x_source == ARR_LEFT_BOUNDARY) ?
              ARR_RIGHT_BOUNDARY : ARR_LEFT_BOUNDARY;
            last_pending = true;
          }
        }
        curr = curr->next();
        continue;
      }

      res_source = (ps_x_source == ARR_LEFT_BOUNDARY) ? LARGER :
        (ps_x_source == ARR_RIGHT_BOUNDARY) ? SMALLER :
        (ps_y_source == ARR_INTERIOR) ?
        cmp_x_op(p, curr->opposite()->vertex()->point()) :
        cmp_x_pt_ce(p, curr->curve(), ind_source);

      res_target = (ps_x_target == ARR_LEFT_BOUNDARY) ? LARGER :
        (ps_x_target == ARR_RIGHT_BOUNDARY) ? SMALLER :
        (ps_y_target == ARR_INTERIOR) ?
        cmp_x_op(p, curr->vertex()->point()) :
        cmp_x_pt_ce(p, curr->curve(), ind_target);

      /* If a vertical ray is shot from p upward, the x-monotone curve
       * associated with curr is hit once.
       */
      if (res_source == res_target) {
        curr = curr->next();
        continue;
      }

      if (res_source != EQUAL) {
        change_pending = true;
        res_pending = (res_source == SMALLER) ? LARGER : SMALLER;
      }
      if (res_target != EQUAL) {
        if (change_pending) {
          change_pending = false;
          if (res_pending == res_target) {
            Comparison_result res_y_at_x = cmp_y_at_x_op(p, curr->curve());
            if (res_y_at_x == EQUAL) return false;
            if (res_y_at_x == SMALLER) num_intersections++;
          }
        } else {
          // This must be the first curve. Remember to check the last curve
          res_last = (res_target == SMALLER) ? LARGER : SMALLER;
          last_pending = true;
        }
      }

      /* Proceed to the next halfedge along the component boundary.
       * Note that the source vertex of this halfedge is the current target.
       */
      curr = curr->next();
    } while (curr != first);

    if (last_pending) {
      if (!p_is_interior_x) {
        if (ps_x_last == ps_x_target) {
          Comparison_result res_y_at_x = cmp_y_at_x_op(p, curr->curve());
          if (res_y_at_x == EQUAL) return false;
          if (res_y_at_x == SMALLER) num_intersections++;
        }
        continue;
      }

      if (res_last == res_source) {
        Comparison_result res_y_at_x = cmp_y_at_x_op(p, curr->curve());
        if (res_y_at_x == EQUAL) return false;
        if (res_y_at_x == SMALLER) num_intersections++;
      }
    }
  }
  /* The query point lies inside the connected components if the face does
   * not contain the north pole, and the vertical ray intersects the
   * boundaries an odd number of times. As mentioned above, if the face does
   * contain the north pole, then it contains everything, (and has no outer
   * CCB's at all).
   */
  return (num_intersections& 0x1);
}

/*! \brief compares the relative y-position of a point and a halfedge */
template <typename GeomTraits, typename Dcel>
Comparison_result
Arr_spherical_topology_traits_2<GeomTraits, Dcel>::
compare_y_at_x(const Point_2& p, const Halfedge* he) const
{
  // std::cout << "compare_y_at_x(Point_2&,Halfedge*)" << std::endl;
  return m_geom_traits->compare_y_at_x_2_object()(p, he->curve());
}

/*! \brief determine whether a vertex is associated with a curve end */
template <typename GeomTraits, typename Dcel>
bool Arr_spherical_topology_traits_2<GeomTraits, Dcel>::
are_equal(const Vertex* v,
          const X_monotone_curve_2& xc, Arr_curve_end ind,
          Arr_parameter_space ps_x, Arr_parameter_space ps_y) const
{
#if 0
  std::cout << "are_equal"
            << ", v: " << v->point()
            << ", xc: " << xc << ", " << ind
            << std::endl;
#endif
  CGAL_precondition(ps_x == ARR_LEFT_BOUNDARY || ps_x == ARR_RIGHT_BOUNDARY ||
                    ps_y == ARR_BOTTOM_BOUNDARY || ps_y == ARR_TOP_BOUNDARY);

  // If the given boundary conditions do not match those of the given
  // vertex, v cannot represent the curve end.
  if (ps_y != v->parameter_space_in_y()) return false;

  if (ps_y != ARR_INTERIOR) return (ps_y == v->parameter_space_in_y());

  if (((ps_x == ARR_INTERIOR) && (v->parameter_space_in_x() != ARR_INTERIOR)) ||
      ((ps_x != ARR_INTERIOR) && (v->parameter_space_in_x() == ARR_INTERIOR)))
    return false;

  CGAL_assertion(ps_x != ARR_INTERIOR);
  /* Both vertices have the same x boundary conditions =>
   * comapare their y-position.
   */
  const Point_2& p1 = v->point();
  const Point_2& p2 = (ind == ARR_MIN_END) ?
    m_geom_traits->construct_min_vertex_2_object()(xc) :
    m_geom_traits->construct_max_vertex_2_object()(xc);
  return (m_geom_traits->compare_y_on_boundary_2_object()(p1, p2) == EQUAL);
}

/*! \brief receives a notification on the creation of a new boundary vertex */
template <typename GeomTraits, typename Dcel>
void
Arr_spherical_topology_traits_2<GeomTraits, Dcel>::
notify_on_boundary_vertex_creation(Vertex* v,
                                   const X_monotone_curve_2& xc,
                                   Arr_curve_end ind,
                                   Arr_parameter_space
                                     CGAL_assertion_code(ps_x),
                                   Arr_parameter_space ps_y)
{
  // std::cout << "notify_on_boundary_vertex_creation()" << std::endl;
  if (ps_y == ARR_BOTTOM_BOUNDARY) {
    m_south_pole = v;
    return;
  }
  if (ps_y == ARR_TOP_BOUNDARY) {
    m_north_pole = v;
    return;
  }
  CGAL_assertion(ps_x != ARR_INTERIOR);
  const Point_2& key = (ind == ARR_MIN_END) ?
    m_geom_traits->construct_min_vertex_2_object()(xc) :
    m_geom_traits->construct_max_vertex_2_object()(xc);
  m_boundary_vertices.insert(Vertex_value(key, v));
}

template <typename GeomTraits, typename Dcel>
bool Arr_spherical_topology_traits_2<GeomTraits, Dcel>::
let_me_decide_the_outer_ccb(std::pair< CGAL::Sign, CGAL::Sign> signs1,
                            std::pair< CGAL::Sign, CGAL::Sign> signs2,
                            bool& swap_predecessors) const
{
  // no perimetric in top-bottom for first loop
  CGAL_precondition(signs1.second == CGAL::ZERO);

  // no perimetric in top-bottom for second loop
  CGAL_precondition(signs2.second == CGAL::ZERO);

  // choose prev1 to define outer ccb of new face if it is a non-perimetric loop,
  // otherwise choose prev2
  // TODO what if both are non-zero? does it occur?
  // TODO EBEB check this!!!!
  swap_predecessors = (signs2.first != CGAL::POSITIVE);

  // but only if the at least one of the loops is perimetric, otherwise return
  // false to let leftmost-vertex decide which becomes part of the new outer ccb
  return (signs1.first != CGAL::ZERO) || (signs2.first != CGAL::ZERO);
}

/*! \brief given a curve end with boundary conditions and a face that contains
 * the interior of the curve, find a place for a boundary vertex that will
 * represent the curve end along the face boundary */
template <typename GeomTraits, typename Dcel>
CGAL::Object
Arr_spherical_topology_traits_2<GeomTraits, Dcel>::
place_boundary_vertex(Face* /* f */,
                      const X_monotone_curve_2& xc, Arr_curve_end ind,
                      Arr_parameter_space
#if !defined(CGAL_NO_ASSERTIONS)
                      ps_x
#endif
                      ,
                      Arr_parameter_space ps_y)
{
  // std::cout << "place_boundary_vertex()" << std::endl;
  if (ps_y == ARR_BOTTOM_BOUNDARY) {
    if (m_south_pole == nullptr) return Object();
    return CGAL::make_object(m_south_pole);
  }

  if (ps_y == ARR_TOP_BOUNDARY) {
    if (m_north_pole == nullptr) return Object();
    return CGAL::make_object(m_north_pole);
  }

  CGAL_assertion((ps_x == ARR_LEFT_BOUNDARY) || (ps_x == ARR_RIGHT_BOUNDARY));

  const Point_2& key = (ind == ARR_MIN_END) ?
    m_geom_traits->construct_min_vertex_2_object()(xc) :
    m_geom_traits->construct_max_vertex_2_object()(xc);
  typename Vertex_map::iterator it = m_boundary_vertices.find(key);

  if (it != m_boundary_vertices.end()) {
    Vertex* v = it->second;
    return CGAL::make_object(v);
  }

  // The vertex hasn't been created yet, return a null object:
  return Object();
}

/*! \brief locate the predecessor halfedge for the given curve around a given
 * vertex with boundary conditions. */
template <typename GeomTraits, typename Dcel>
typename Arr_spherical_topology_traits_2<GeomTraits, Dcel>::Halfedge*
Arr_spherical_topology_traits_2<GeomTraits,Dcel>::
locate_around_boundary_vertex(Vertex* v,
                              const X_monotone_curve_2& xc,
                              Arr_curve_end ind,
                              Arr_parameter_space
#if !defined(CGAL_NO_ASSERTIONS)
                              ps_x
#endif
                              ,
                              Arr_parameter_space ps_y) const
{
  // std::cout << "locate_around_boundary_vertex()" << std::endl;
  if (ps_y == ARR_BOTTOM_BOUNDARY) {
    CGAL_assertion(v == m_south_pole);
    return (_locate_around_pole(m_south_pole, xc, ind));
  }

  if (ps_y == ARR_TOP_BOUNDARY) {
    CGAL_assertion(v == m_north_pole);
    return (_locate_around_pole(m_north_pole, xc, ind));
  }

  CGAL_assertion((ps_x == ARR_LEFT_BOUNDARY) || (ps_x == ARR_RIGHT_BOUNDARY));

  return (_locate_around_vertex_on_discontinuity(v, xc, ind));
}

/*! \brief locates a DCEL feature that contains a given curve end. */
template <typename GeomTraits, typename Dcel>
CGAL::Object Arr_spherical_topology_traits_2<GeomTraits, Dcel>::
locate_curve_end(const X_monotone_curve_2& xc, Arr_curve_end ind,
                 Arr_parameter_space
#if !defined(CGAL_NO_ASSERTIONS)
                 ps_x
#endif
                 ,
                 Arr_parameter_space ps_y)
{
  // Act according to the boundary conditions.
  if (ps_y == ARR_TOP_BOUNDARY) {
    // In case the curve end coincides with the north pole, return the vertex
    // representing the north pole, if one exists. Otherwise, return the face
    // containing this pole (the spherical face).
    if (m_north_pole != nullptr) return CGAL::make_object(m_north_pole);
    return CGAL::make_object(m_spherical_face);
  }

  typename Vertex_map::iterator it;
  Vertex* v = nullptr;

  if (ps_y == ARR_BOTTOM_BOUNDARY) {
    // In case the curve end coincides with the south pole, return the vertex
    // representing the south pole, if one exists. Otherwise, search for the
    // face containing this pole.
    if (m_south_pole != nullptr) return CGAL::make_object(m_south_pole);
    it = m_boundary_vertices.begin();
  }
  else {
    CGAL_assertion((ps_x == ARR_LEFT_BOUNDARY) || (ps_x == ARR_RIGHT_BOUNDARY));

    // Check if the given curve end is incident to a vertex on the line of
    // discontinuity. If so, return this vertex. Otherwise, locate the first
    // vertex above it.
    const Point_2& key = (ind == ARR_MIN_END) ?
      m_geom_traits->construct_min_vertex_2_object()(xc) :
      m_geom_traits->construct_max_vertex_2_object()(xc);
    it = m_boundary_vertices.find(key);
    if (it != m_boundary_vertices.end()) {
      v = it->second;
      return CGAL::make_object(v);
    }

    it = m_boundary_vertices.lower_bound(key);
  }

  // At this point, the iterator it points to a vertex on the line of
  // discontinuity that is strictly above the curve end. If there is none,
  // we know the curve end is contained in the spherical face. Otherwise,
  // we return the face that lies below the vertex v.
  if (it == m_boundary_vertices.end())
    return CGAL::make_object(m_spherical_face);

  v = it->second;
  return CGAL::make_object(_face_below_vertex_on_discontinuity(v));
}

/*! \brief determines whether a given boundary vertex is redundant */
template <typename GeomTraits, typename Dcel>
bool Arr_spherical_topology_traits_2<GeomTraits, Dcel>::
is_redundant(const Vertex* v) const
{ return (v->halfedge() == nullptr); }

/* \brief erases a given redundant vertex */
template <typename GeomTraits, typename Dcel>
typename Arr_spherical_topology_traits_2<GeomTraits, Dcel>::Halfedge*
Arr_spherical_topology_traits_2<GeomTraits, Dcel>::
erase_redundant_vertex(Vertex* v)
{
  const Arr_parameter_space ps_y = v->parameter_space_in_y();
  if (ps_y == ARR_BOTTOM_BOUNDARY) {
    m_south_pole = nullptr;
    return nullptr;
  }
  if (ps_y == ARR_TOP_BOUNDARY) {
    m_north_pole = nullptr;
    return nullptr;
  }
  CGAL_assertion_code(Arr_parameter_space ps_x = v->parameter_space_in_x());
  CGAL_assertion(ps_x != ARR_INTERIOR);
  m_boundary_vertices.erase(v->point());
  return nullptr;
}

/*! \brief obtains the curve associated with a boundary vertex */
template <typename GeomTraits, typename Dcel>
const typename
Arr_spherical_topology_traits_2<GeomTraits, Dcel>::X_monotone_curve_2&
Arr_spherical_topology_traits_2<GeomTraits, Dcel>::
_curve(const Vertex* v, Arr_curve_end& ind) const
{
  // std::cout << "curve()" << std::endl;
  const Halfedge* he = v->halfedge();
  ind = (he->direction() == ARR_LEFT_TO_RIGHT) ? ARR_MAX_END : ARR_MIN_END;
  return he->curve();
}

/*! \brief returns the halfedge, the target vertex of which is given, that is
 * the predecessor of a halfedge, the curve of which is given, that is about
 * to be inserted into the dcel.
 */
template <typename GeomTraits, typename Dcel>
typename Arr_spherical_topology_traits_2<GeomTraits, Dcel>::Halfedge*
Arr_spherical_topology_traits_2<GeomTraits, Dcel>::
_locate_around_vertex_on_discontinuity(Vertex* v,
                                       const X_monotone_curve_2& xc,
                                       Arr_curve_end ind) const
{
  // If the vertex is isolated, there is no predecssor halfedge.
  if (v->is_isolated()) return nullptr;

  // Get the first incident halfedge around v and the next halfedge.
  Halfedge* first = v->halfedge();
  Halfedge* curr = first;
  CGAL_assertion(curr != nullptr);
  Halfedge* next = curr->next()->opposite();

  // If is only one halfedge incident to v, return this halfedge as xc's
  // predecessor:
  if (curr == next) return curr;

  // Otherwise, we traverse the halfedges around v until we find the pair
  // of adjacent halfedges between which we should insert xc.
  typename Gt_adaptor_2::Is_between_cw_2 is_between_cw =
    m_geom_traits->is_between_cw_2_object();
  bool eq_curr, eq_next;

  while (!is_between_cw(xc, (ind == ARR_MIN_END), curr->curve(),
                        (curr->direction() == ARR_RIGHT_TO_LEFT), next->curve(),
                        (next->direction() == ARR_RIGHT_TO_LEFT), v->point(),
                        eq_curr, eq_next))
  {
    // The curve must not be equal to one of the curves already incident to v.
    CGAL_assertion(!eq_curr && !eq_next);

    // Move to the next pair of incident halfedges.
    curr = next;
    next = curr->next()->opposite();

    // Make sure we have not completed a full traversal around v without
    // locating a place for the new curve xc.
    CGAL_assertion(curr != first);
  }

  // Return the halfedge we have located.
  return curr;
}

/*! \brief returns the halfedge, the target vertex of which is a given pole,
 * that is the predecessor of a halfedge, the curve of which is given, that
 * is about to be inserted into the dcel.
 */
template <typename GeomTraits, typename Dcel>
typename Arr_spherical_topology_traits_2<GeomTraits, Dcel>::Halfedge*
Arr_spherical_topology_traits_2<GeomTraits, Dcel>::
_locate_around_pole(Vertex* v,
                    const X_monotone_curve_2& xc, Arr_curve_end ind) const
{
  CGAL_assertion((v == m_south_pole) || (v == m_north_pole));

  // std::cout << "locate_around_pole() " << ind << std::endl;
  // If the vertex is isolated, return a null halfedge:
  if (v->is_isolated()) return nullptr;

  // Get the first incident halfedge around v and the next halfedge:
  Halfedge* first = v->halfedge();
  Halfedge* curr = first;
  CGAL_assertion(curr != nullptr);
  Halfedge* next = curr->next()->opposite();

  // If there is only one halfedge, it is the predecessor, return it:
  if (curr == next) return curr;

  // If we compare a curve and its successor around the south (resp. north)
  // pole, the result LARGER (resp. SMALLER) indicates that the line of
  // discontinuity is located in between the two curves.
  const Comparison_result cross_res = (v == m_south_pole) ? LARGER : SMALLER;

  // Traverse all other halfedges, and compare their x-positions next to the
  // pole with the query curve xc.
  typename Gt_adaptor_2::Compare_x_curve_ends_2 cmp_x_curve_ends =
    m_geom_traits->compare_x_curve_ends_2_object();
  Arr_curve_end curr_end, next_end;
  Comparison_result curr_res, next_res;
  Comparison_result curr_next_res;

  curr_end =
    (curr->direction() == ARR_RIGHT_TO_LEFT) ? ARR_MIN_END : ARR_MAX_END;
  curr_res = cmp_x_curve_ends(xc, ind, curr->curve(), curr_end);
  do {
    next_end =
      (next->direction() == ARR_RIGHT_TO_LEFT) ? ARR_MIN_END : ARR_MAX_END;
    next_res = cmp_x_curve_ends(xc, ind, next->curve(), next_end);
    curr_next_res =
      cmp_x_curve_ends(curr->curve(), curr_end, next->curve(), next_end);
    if (curr_next_res == cross_res) {
      // The line of discontinuity must lie between curr and next, so the
      // comparison result of xc with the two curves should be equal:
      if (curr_res == next_res) return curr;
    }
    else {
      // The line of discontinuity does not lie between curr and next, so the
      // comparison results must be different if xc lies in between.
      if (curr_res != next_res) return curr;
    }

    // Move to the next halfedge around the pole.
    curr = next;
    curr_end = next_end;
    curr_res = next_res;
    next = curr->next()->opposite();
  } while (curr != first);

  // We sould never reach here:
  CGAL_error();
  return nullptr;
}

/*! \brief Return the face that lies below the given vertex, which lies
 * on the line of discontinuity.
 */
template <typename GeomTraits, typename Dcel>
typename Arr_spherical_topology_traits_2<GeomTraits, Dcel>::Face*
Arr_spherical_topology_traits_2<GeomTraits, Dcel>::
_face_below_vertex_on_discontinuity(Vertex* v) const
{
  // If the vertex is isolated, just return the face that contains it.
  if (v->is_isolated()) return (v->isolated_vertex()->face());

  // Get the first incident halfedge around v and the next halfedge.
  Halfedge* first = v->halfedge();
  Halfedge* curr = first;
  CGAL_assertion(curr != nullptr);
  Halfedge* next = curr->next()->opposite();

  // If there is only one halfedge incident to v, return its incident face.
  if (curr == next)
    return ((curr->is_on_inner_ccb()) ?
            curr->inner_ccb()->face() : curr->outer_ccb()->face());

  // Otherwise, we traverse the halfedges around v and locate the first
  // halfedge we encounter if we go from "6 o'clock" clockwise.
  // First locate the lower left and the top right halfedges around v.
  typename Gt_adaptor_2::Compare_y_at_x_right_2 cmp_y_at_x_op_right =
    m_geom_traits->compare_y_at_x_right_2_object();
  typename Gt_adaptor_2::Compare_y_at_x_left_2  cmp_y_at_x_op_left =
    m_geom_traits->compare_y_at_x_left_2_object();

  Halfedge* lowest_left = nullptr;
  Halfedge* top_right = nullptr;

  do {
    // Check whether the current halfedge is defined to the left or to the
    // right of the given vertex.
    if (curr->direction() == ARR_LEFT_TO_RIGHT) {
      // The curve associated with the current halfedge is defined to the left
      // of v.
      if (lowest_left == nullptr ||
          cmp_y_at_x_op_left(curr->curve(), lowest_left->curve(), v->point())
          == SMALLER)
      {
        lowest_left = curr;
      }
    }
    else {
      // The curve associated with the current halfedge is defined to the right
      // of v.
      if (top_right == nullptr ||
          cmp_y_at_x_op_right(curr->curve(), top_right->curve(), v->point()) ==
          LARGER)
      {
        top_right = curr;
      }
    }

    // Move to the next halfedge around the vertex.
    curr = curr->next()->opposite();
  } while (curr != first);

  // The first halfedge we encounter is the lowest to the left, but if there
  // is no edge to the left, we first encounter the topmost halfedge to the
  // right. Note that as the halfedge we located has v as its target, we now
  // have to return its twin.
  first =
    (lowest_left != nullptr) ? lowest_left->opposite() : top_right->opposite();
  // std::cout << "first: " << first->opposite()->vertex()->point() << " => "
  //           << first->vertex()->point() << std::endl;

  // Face* f = (first->is_on_inner_ccb()) ?
  //   first->inner_ccb()->face() : first->outer_ccb()->face();
  // std::cout << "outer: " << f->number_of_outer_ccbs() << std::endl;
  // std::cout << "inner: " << f->number_of_inner_ccbs() << std::endl;
  // Return the incident face.
  return ((first->is_on_inner_ccb()) ?
          first->inner_ccb()->face() : first->outer_ccb()->face());
}

} //namespace CGAL

#endif
