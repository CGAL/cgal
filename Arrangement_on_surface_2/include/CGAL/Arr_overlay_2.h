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
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>
//                 Efi Fogel <efif@post.tau.ac.il>

#ifndef CGAL_ARR_OVERLAY_2_H
#define CGAL_ARR_OVERLAY_2_H

#include <CGAL/license/Arrangement_on_surface_2.h>


/*! \file
 * Definition of the global Arr_overlay_2() function.
 */

#include <boost/optional/optional.hpp>

#include <CGAL/Arrangement_on_surface_2.h>
#include <CGAL/Sweep_line_2.h>
#include <CGAL/Sweep_line_2/Arr_default_overlay_traits_base.h>

#include <vector>
#include <boost/mpl/if.hpp>
#include <boost/mpl/or.hpp>
#include <boost/type_traits.hpp>

#include <CGAL/assertions.h>

namespace CGAL {

/*!
 * Compute the overlay of two input arrangements.
 * \param arr1 The first arrangement.
 * \param arr2 The second arrangement.
 * \param arr_res Output: The resulting arrangement.
 * \param ovl_tr An overlay-traits class. As arr1, arr2 and res can be
 *               templated with different geometry-traits class and
 *               different DCELs (encapsulated in the various topology-traits
 *               classes). The geometry-traits of the result arrangement is
 *               used to construct the result arrangement. This means that all
 *               the types (e.g., Point_2, Curve_2 and X_monotone_2) of both
 *               arr1 and arr2 have to be convertible to the types
 *               in the result geometry-traits.
 *               The overlay-traits class defines the various
 *               overlay operations of pairs of DCEL features from
 *               TopTraitsA and TopTraitsB to the resulting ResDcel.
 */
template <class GeomTraitsA,
          class GeomTraitsB,
          class GeomTraitsRes,
          class TopTraitsA,
          class TopTraitsB,
          class TopTraitsRes,
          class OverlayTraits>
void overlay(const Arrangement_on_surface_2<GeomTraitsA, TopTraitsA>& arr1,
             const Arrangement_on_surface_2<GeomTraitsB, TopTraitsB>& arr2,
             Arrangement_on_surface_2<GeomTraitsRes, TopTraitsRes>& arr_res,
             OverlayTraits& ovl_tr)
{
  typedef Arrangement_on_surface_2<GeomTraitsA, TopTraitsA>       ArrA;
  typedef Arrangement_on_surface_2<GeomTraitsB, TopTraitsB>       ArrB;
  typedef Arrangement_on_surface_2<GeomTraitsRes, TopTraitsRes>   ArrRes;

  // some type assertions (not all, but better then nothing).
  CGAL_static_assertion
    ((boost::is_convertible<typename GeomTraitsA::Point_2, \
                            typename GeomTraitsRes::Point_2>::value));
  CGAL_static_assertion
    ((boost::is_convertible<typename GeomTraitsB::Point_2, \
                            typename GeomTraitsRes::Point_2 >::value));
  CGAL_static_assertion
    ((boost::is_convertible<typename GeomTraitsA::X_monotone_curve_2, \
                            typename GeomTraitsRes::X_monotone_curve_2>::value));
  CGAL_static_assertion
    ((boost::is_convertible<typename GeomTraitsB::X_monotone_curve_2, \
                            typename GeomTraitsRes::X_monotone_curve_2>::value));

  typedef typename TopTraitsRes::template
    Sweep_line_overlay_visitor<ArrA, ArrB, OverlayTraits>
                                                      Ovl_visitor;

  typedef typename Ovl_visitor::Traits_2              Ovl_traits_2;
  typedef typename Ovl_traits_2::X_monotone_curve_2   Ovl_x_monotone_curve_2;
  typedef typename Ovl_traits_2::Point_2              Ovl_point_2;

  typedef typename Ovl_traits_2::Cell_handle_red      Cell_handle_red;
  typedef typename Ovl_traits_2::Optional_cell_red    Optional_cell_red;
  typedef typename Ovl_traits_2::Cell_handle_blue     Cell_handle_blue;
  typedef typename Ovl_traits_2::Optional_cell_blue   Optional_cell_blue;

  CGAL_USE_TYPE(Optional_cell_red);
  CGAL_USE_TYPE(Optional_cell_blue);

  // The result arrangement cannot be on of the input arrangements.
  CGAL_precondition(((void *)(&arr_res) != (void *)(&arr1)) &&
                    ((void *)(&arr_res) != (void *)(&arr2)));

  // Prepare a vector of extended x-monotone curves that represent all edges
  // in both input arrangements. Each curve is associated with a halfedge
  // directed from right to left.
  typename ArrA::Halfedge_const_handle invalid_he1;
  typename ArrB::Halfedge_const_handle invalid_he2;
  std::vector<Ovl_x_monotone_curve_2>
    xcvs_vec(arr1.number_of_edges() + arr2.number_of_edges());
  unsigned int i = 0;

  typename ArrA::Edge_const_iterator eit1;
  for (eit1 = arr1.edges_begin(); eit1 != arr1.edges_end(); ++eit1, ++i) {
    typename ArrA::Halfedge_const_handle he1 = eit1;
    if (he1->direction() != ARR_RIGHT_TO_LEFT) he1 = he1->twin();
    xcvs_vec[i] = Ovl_x_monotone_curve_2(eit1->curve(), he1, invalid_he2);
  }

  typename ArrB::Edge_const_iterator eit2;
  for (eit2 = arr2.edges_begin(); eit2 != arr2.edges_end(); ++eit2, ++i) {
    typename ArrB::Halfedge_const_handle he2 = eit2;
    if (he2->direction() != ARR_RIGHT_TO_LEFT) he2 = he2->twin();
    xcvs_vec[i] = Ovl_x_monotone_curve_2(eit2->curve(), invalid_he1, he2);
  }

  // Obtain a extended traits-class object and define the sweep-line visitor.
  const typename ArrRes::Traits_adaptor_2* traits_adaptor =
    arr_res.traits_adaptor();

  /* We would like to avoid copy construction of the geometry traits class.
   * Copy construction is undesired, because it may results with data
   * duplication or even data loss.
   *
   * If the type Ovl_traits_2 is the same as the type
   * GeomTraits, use a reference to GeomTraits to avoid constructing a new one.
   * Otherwise, instantiate a local variable of the former and provide
   * the later as a single parameter to the constructor.
   *
   * Use the form 'A a(*b);' and not ''A a = b;' to handle the case where A has
   * only an implicit constructor, (which takes *b as a parameter).
   */
  typedef Arr_traits_basic_adaptor_2< GeomTraitsRes > Geom_traits_adaptor_2;
  typename boost::mpl::if_<boost::is_same< Geom_traits_adaptor_2, Ovl_traits_2>,
                           const Ovl_traits_2&, Ovl_traits_2>::type
    ex_traits(*traits_adaptor);

  Ovl_visitor visitor(&arr1, &arr2, &arr_res, &ovl_tr);
  Sweep_line_2<Ovl_traits_2, Ovl_visitor,
               typename Ovl_visitor::Subcurve, typename Ovl_visitor::Event>
    sweep_line(&ex_traits, &visitor);

  // In case both arrangement do not contain isolated vertices, go on and
  // overlay them.
  const std::size_t total_iso_verts =
    arr1.number_of_isolated_vertices() + arr2.number_of_isolated_vertices();

  if (total_iso_verts == 0) {
    // Clear the result arrangement and perform the sweep to construct it.
    arr_res.clear();
    sweep_line.sweep(xcvs_vec.begin(), xcvs_vec.end());
    xcvs_vec.clear();
    return;
  }

  // Prepare a vector of extended points that represent all isolated vertices
  // in both input arrangements.
  std::vector<Ovl_point_2> pts_vec(total_iso_verts);

  i = 0;
  typename ArrA::Vertex_const_iterator  vit1;
  for (vit1 = arr1.vertices_begin(); vit1 != arr1.vertices_end(); ++vit1) {
    if (vit1->is_isolated()) {
      typename ArrA::Vertex_const_handle v1 = vit1;
      pts_vec[i++] =
        Ovl_point_2(vit1->point(), boost::make_optional(Cell_handle_red(v1)),
                    boost::optional<Cell_handle_blue>());
    }
  }

  typename ArrB::Vertex_const_iterator  vit2;
  for (vit2 = arr2.vertices_begin(); vit2 != arr2.vertices_end(); ++vit2) {
    if (vit2->is_isolated()) {
      typename ArrB::Vertex_const_handle v2 = vit2;
      pts_vec[i++] =
        Ovl_point_2(vit2->point(), boost::optional<Cell_handle_red>(),
                    boost::make_optional(Cell_handle_blue(v2)));
    }
  }

  // Clear the result arrangement and perform the sweep to construct it.
  arr_res.clear();
  sweep_line.sweep(xcvs_vec.begin(), xcvs_vec.end(),
                   pts_vec.begin(), pts_vec.end());
  xcvs_vec.clear();
  pts_vec.clear();
}

/*!
 * Compute the (simple) overlay of two input arrangements.
 * \param arr1 The first arrangement.
 * \param arr2 The second arrangement.
 * \param arr_res Output: The resulting arrangement.
 */
template <class GeomTraitsA,
          class GeomTraitsB,
          class GeomTraitsRes,
          class TopTraitsA,
          class TopTraitsB,
          class TopTraitsRes>
void overlay(const Arrangement_on_surface_2<GeomTraitsA, TopTraitsA>& arr1,
             const Arrangement_on_surface_2<GeomTraitsB, TopTraitsB>& arr2,
             Arrangement_on_surface_2<GeomTraitsRes, TopTraitsRes>& arr_res)
{
  typedef Arrangement_on_surface_2<GeomTraitsA, TopTraitsA>     ArrA;
  typedef Arrangement_on_surface_2<GeomTraitsA, TopTraitsB>     ArrB;
  typedef Arrangement_on_surface_2<GeomTraitsRes, TopTraitsRes> ArrRes;
  _Arr_default_overlay_traits_base<ArrA, ArrB, ArrRes> ovl_traits;
  overlay(arr1, arr2, arr_res, ovl_traits);
}

} //namespace CGAL

#endif
