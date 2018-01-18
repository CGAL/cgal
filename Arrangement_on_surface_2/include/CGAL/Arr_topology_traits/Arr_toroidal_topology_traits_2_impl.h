// Copyright (c) 2006,2007,2008,2009,2010,2011,2013,2014 Max-Planck-Institute Saarbruecken (Germany), Tel-Aviv University (Israel).
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
// Author(s)     : Efi Fogel         <efif@post.tau.ac.il>
//                 Ron Wein          <wein@post.tau.ac.il>
//                 Eric Berberich    <eric.berberich@cgal.org>

#ifndef CGAL_ARR_TOROIDAL_TOPOLOGY_TRAITS_2_IMPL_H
#define CGAL_ARR_TOROIDAL_TOPOLOGY_TRAITS_2_IMPL_H

/*! \file
 * Member-function definitions for the
 * Arr_toroidal_topology_traits_2<GeomTraits> class.
 */

namespace CGAL {

/*! \brief constructs default */
template <class GeomTraits, class Dcel>
Arr_toroidal_topology_traits_2<GeomTraits, Dcel>::
Arr_toroidal_topology_traits_2() :
  m_reference_face(NULL),
  m_own_traits(true)
{
  m_traits = new Traits_adaptor_2;
}

/*! \brief constructs with a geometry-traits class */
template <class GeomTraits, class Dcel>
Arr_toroidal_topology_traits_2<GeomTraits, Dcel>::
Arr_toroidal_topology_traits_2(const Geometry_traits_2* traits) :
  m_reference_face(NULL),
  m_own_traits(false)
{
  m_traits = static_cast<const Traits_adaptor_2*>(traits);
}

/*! \brief assigns the contents of another topology-traits class */
template <class GeomTraits, class Dcel>
void
Arr_toroidal_topology_traits_2<GeomTraits, Dcel>::
assign(const Self& other)
{
  CGAL_error();
  /* dummy implementation */

  // Clear the current DCEL and duplicate the other DCEL.
  m_dcel.delete_all();
  m_dcel.assign(other.m_dcel);

  // Take care of the traits object.
  if (m_own_traits && m_traits != NULL)
    delete m_traits;

  if (other.m_own_traits)
  {
    m_traits = new Traits_adaptor_2;
    m_own_traits = true;
  }
  else
  {
    m_traits = other.m_traits;
    m_own_traits = false;
  }

  // Update the rest of the properties.
  dcel_updated();

  return;
}

/*! Determine whether the DCEL reprsenets an empty structure. */
template <class GeomTraits, class Dcel>
bool
Arr_toroidal_topology_traits_2<GeomTraits, Dcel>::
is_empty_dcel() const
{
  return (m_dcel.size_of_vertices() == 0);
}

/*! \brief initializes an empty DCEL structure. */
template <class GeomTraits, class Dcel>
void
Arr_toroidal_topology_traits_2<GeomTraits, Dcel>::
dcel_updated()
{
  // Go over the DCEL vertices and locate the south and north pole (if any)
  // and any other vertex on the line of discontinuity.
  typename Dcel::Vertex_iterator       vit;
  Arr_parameter_space                  bx, by;

  CGAL_error();
  /* dummy implementation */ return;
}

/*! \brief initializes an empty DCEL structure. */
template <class GeomTraits, class Dcel>
void
Arr_toroidal_topology_traits_2<GeomTraits, Dcel>::
init_dcel()
{
#if CGAL_ARR_TOPOLOGY_TRAITS_VERBOSE
  std::cout << "TOR-TOP: init_dcel()" << std::endl;
#endif

  // Clear the current DCEL.
  m_dcel.delete_all();

  // create the reference face
  this->m_reference_face = this->m_dcel.new_face();

#if CGAL_ARR_TOPOLOGY_TRAITS_VERBOSE
  std::cout << "TOR-TOP reference face: " << this->m_reference_face << std::endl;
#endif

  // and set properties
  // bounded
  this->m_reference_face->set_unbounded (false);
  // set not fictious
  this->m_reference_face->set_fictitious (false);

  // TODO create identifications

}

/*! Determine whether the given vertex is concrete. */
template <class GeomTraits, class Dcel>
bool
Arr_toroidal_topology_traits_2<GeomTraits, Dcel>::
is_concrete_vertex(const Vertex* v) const
{
  // every vertex has geometry!
  CGAL_precondition(!v->has_null_point());
  return true;
}

/*! Obtain the number of concrete vertices. */
template <class GeomTraits, class Dcel>
typename Arr_toroidal_topology_traits_2<GeomTraits, Dcel>::Size
Arr_toroidal_topology_traits_2<GeomTraits, Dcel>::
number_of_concrete_vertices() const
{
  // every vertex has geometry
  return (m_dcel.size_of_vertices());
}

/*! Determine whether the given vertex is valid. */
template <class GeomTraits, class Dcel>
bool
Arr_toroidal_topology_traits_2<GeomTraits, Dcel>::
is_valid_vertex (const Vertex* /* v */) const {
  // all vertices are valid
  return false;
}

/*! Obtain the number of valid vertices. */
template <class GeomTraits, class Dcel>
typename Arr_toroidal_topology_traits_2<GeomTraits, Dcel>::Size
Arr_toroidal_topology_traits_2<GeomTraits, Dcel>::
number_of_valid_vertices() const
{
  // all vertices are valid
  return (m_dcel.size_of_vertices());
}

/*! Determine whether the given halfedge is valid. */
template <class GeomTraits, class Dcel>
bool
Arr_toroidal_topology_traits_2<GeomTraits, Dcel>::
is_valid_halfedge (const Halfedge* /* he */) const
{
  // all halfedges are valid
  return true;
}

/*! Obtain the number of valid halfedges. */
template <class GeomTraits, class Dcel>
typename Arr_toroidal_topology_traits_2<GeomTraits, Dcel>::Size
Arr_toroidal_topology_traits_2<GeomTraits, Dcel>::
number_of_valid_halfedges() const
{
  // all halfedges are valid
  return (m_dcel.size_of_halfedges());
}

/*! Determine whether the given face is valid. */
template <class GeomTraits, class Dcel>
bool
Arr_toroidal_topology_traits_2<GeomTraits, Dcel>::
is_valid_face (const Face* /* f */) const
{
  // each face is valid
  return true;
}

/*! Obtain the number of valid faces. */
template <class GeomTraits, class Dcel>
typename Arr_toroidal_topology_traits_2<GeomTraits, Dcel>::Size
Arr_toroidal_topology_traits_2<GeomTraits, Dcel>::
number_of_valid_faces() const
{
  // each face is valid
  return (m_dcel.size_of_faces());
}

//! \brief receives a notification on the creation of a new boundary vertex.
template <typename GeomTraits, typename Dcel>
void
Arr_toroidal_topology_traits_2<GeomTraits, Dcel>::
notify_on_boundary_vertex_creation(Vertex* v,
                                   const Point_2& p,
                                   Arr_parameter_space
                                     CGAL_assertion_code(ps_x),
                                   Arr_parameter_space ps_y)
{
  CGAL_error();
  /* dummy implementation */ return;
}

/*! \brief receives a notification on the creation of a new boundary vertex */
template <class GeomTraits, class Dcel>
void
Arr_toroidal_topology_traits_2<GeomTraits, Dcel>::
notify_on_boundary_vertex_creation(Vertex* v,
                                   const X_monotone_curve_2& xc,
                                   Arr_curve_end ind,
                                   Arr_parameter_space
                                     CGAL_assertion_code(ps_x),
                                   Arr_parameter_space ps_y)
{
  CGAL_error();
  /* dummy implementation */ return;
}

/*! Determines whether the function should decide on swapping the predecessor halfedges. */
template <class GeomTraits, class Dcel>
bool
Arr_toroidal_topology_traits_2<GeomTraits, Dcel>::
let_me_decide_the_outer_ccb(std::pair< CGAL::Sign, CGAL::Sign> signs1,
                            std::pair< CGAL::Sign, CGAL::Sign> signs2,
                            bool& swap_predecessors) const
{
  CGAL_error();
  /* dummy implementation */ return false;
}

/*! Given signs of two ccbs that show up when splitting upon insertion of
 * curve into two, determine what happens to the face(s). */
template <class GeomTraits, class Dcel>
std::pair<bool, bool>
Arr_toroidal_topology_traits_2<GeomTraits, Dcel>::
face_split_after_edge_insertion(std::pair< CGAL::Sign, CGAL::Sign > /* signs1 */,
                                std::pair< CGAL::Sign, CGAL::Sign > /* signs2 */) const
{
  CGAL_error();
  /* dummy implementation */ return std::make_pair(true, false);
}

/*! \brief determines whether a point lies in the interior of a given face. */
template <class GeomTraits, class Dcel>
bool Arr_toroidal_topology_traits_2<GeomTraits, Dcel>::
is_in_face(const Face* f, const Point_2& p, const Vertex* v) const
{
  // std::cout << "is_in_face()" << std::endl;
  CGAL_precondition(v == NULL || !v->has_null_point());
  CGAL_precondition(v == NULL || m_traits->equal_2_object()(p, v->point()));

  CGAL_error();
  /* dummy implementation */ return false;
}

/*! \brief compares the relative y-position of a point and a halfedge */
template <class GeomTraits, class Dcel>
Comparison_result
Arr_toroidal_topology_traits_2<GeomTraits, Dcel>::
compare_y_at_x(const Point_2& p, const Halfedge* he) const
{
 CGAL_error();
 /* dummy implementation */ return CGAL::LARGER;
}

/*! \brief determine whether a vertex is associated with a curve end */
template <class GeomTraits, class Dcel>
bool Arr_toroidal_topology_traits_2<GeomTraits, Dcel>::
are_equal(const Vertex* v,
          const X_monotone_curve_2& xc, Arr_curve_end ind,
          Arr_parameter_space ps_x, Arr_parameter_space ps_y) const
{
  CGAL_error();
  /* dummy implementation */ return false;
}

/*! \brief given a curve end with boundary conditions and a face that contains
 * the interior of the curve, find a place for a boundary vertex that will
 * represent the curve end along the face boundary */
template <class GeomTraits, class Dcel>
CGAL::Object
Arr_toroidal_topology_traits_2<GeomTraits, Dcel>::
place_boundary_vertex(Face* /* f */,
                      const X_monotone_curve_2& xc, Arr_curve_end ind,
                      Arr_parameter_space ps_x, Arr_parameter_space ps_y)
{
  CGAL_error();
  /* dummy implementation */ return CGAL::Object();
}

/*! \brief locate the predecessor halfedge for the given curve around a given
 * vertex with boundary conditions. */
template <class GeomTraits, class Dcel>
typename Arr_toroidal_topology_traits_2<GeomTraits, Dcel>::Halfedge*
Arr_toroidal_topology_traits_2<GeomTraits,Dcel>::
locate_around_boundary_vertex(Vertex* v,
                              const X_monotone_curve_2& xc,
                              Arr_curve_end ind,
                              Arr_parameter_space ps_x,
                              Arr_parameter_space ps_y) const
{
  CGAL_error();
  /* dummy implementation */ return NULL;
}

/*! \brief locates a DCEL feature that contains a given curve end. */
template <class GeomTraits, class Dcel>
CGAL::Object Arr_toroidal_topology_traits_2<GeomTraits, Dcel>::
locate_curve_end(const X_monotone_curve_2& xc, Arr_curve_end ind,
                 Arr_parameter_space ps_x, Arr_parameter_space ps_y)
{
  CGAL_error();
  /* dummy implementation */ return CGAL::Object();
}

template <class GeomTraits, class Dcel>
typename Arr_toroidal_topology_traits_2<GeomTraits, Dcel>::Halfedge*
Arr_toroidal_topology_traits_2<GeomTraits, Dcel>::
split_fictitious_edge(Halfedge* /* e */, Vertex* /* v */)
{
  CGAL_error();
  /* should never be called, as there are no fictitious edges */ return NULL;
}

/*! \brief determines whether a given boundary vertex is redundant */
template <class GeomTraits, class Dcel>
bool
Arr_toroidal_topology_traits_2<GeomTraits, Dcel>::
is_unbounded(const Face* f) const
{
  // no face is unbounded
  return false;
}

/*! \brief determines whether a given boundary vertex is redundant */
template <class GeomTraits, class Dcel>
bool
Arr_toroidal_topology_traits_2<GeomTraits, Dcel>::
is_redundant(const Vertex* v) const
{
  CGAL_error();
  /* dummy implementation */ return (v->halfedge() == NULL);
}

/* \brief erases a given redundant vertex */
template <class GeomTraits, class Dcel>
typename Arr_toroidal_topology_traits_2<GeomTraits, Dcel>::Halfedge*
Arr_toroidal_topology_traits_2<GeomTraits, Dcel>::
erase_redundant_vertex(Vertex* v)
{
  CGAL_error();
  /* dummy implementation */ return NULL;
}

} //namespace CGAL

#endif // CGAL_ARR_TOROIDAL_TOPOLOGY_TRAITS_2_IMPL_H
