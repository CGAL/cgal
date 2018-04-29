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
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>
//                 Efi Fogel <efifogel@gmail.com>

#ifndef CGAL_ARR_OVERLAY_2_H
#define CGAL_ARR_OVERLAY_2_H

#include <CGAL/license/Arrangement_on_surface_2.h>

#include <CGAL/disable_warnings.h>

/*! \file
 *
 * Definition of the global Arr_overlay_2() function.
 */

#include <vector>
#include <boost/optional/optional.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/or.hpp>
#include <boost/type_traits.hpp>

#include <CGAL/Arrangement_on_surface_2.h>
#include <CGAL/Surface_sweep_2.h>
#include <CGAL/Surface_sweep_2/Arr_default_overlay_traits_base.h>
#include <CGAL/Surface_sweep_2/Arr_overlay_traits_2.h>
#include <CGAL/Surface_sweep_2/Arr_overlay_ss_visitor.h>
#include <CGAL/Surface_sweep_2/Arr_overlay_event.h>
#include <CGAL/Surface_sweep_2/Arr_overlay_subcurve.h>
#include <CGAL/assertions.h>

namespace CGAL {

/*! Compute the overlay of two input arrangements.
 * \tparam GeometryTraitsA_2 the geometry traits of the first arrangement.
 * \tparam GeometryTraitsB_2 the geometry traits of the second arrangement.
 * \tparam GeometryTraitsRes_2 the geometry traits of the resulting arrangement.
 * \tparam TopologyTraitsA the topology traits of the first arrangement.
 * \tparam TopologyTraitsB the topology traits of the second arrangement.
 * \tparam TopologyTraitsRes the topology traits of the resulting arrangement.
 * \tparam OverlayTraits An overlay-traits class. As arr1, arr2 and res can be
 *               templated with different geometry-traits class and
 *               different DCELs (encapsulated in the various topology-traits
 *               classes). The geometry-traits of the result arrangement is
 *               used to construct the result arrangement. This means that all
 *               the types (e.g., Point_2, Curve_2 and X_monotone_2) of both
 *               arr1 and arr2 have to be convertible to the types
 *               in the result geometry-traits.
 *               The overlay-traits class defines the various
 *               overlay operations of pairs of DCEL features from
 *               TopologyTraitsA and TopologyTraitsB to the resulting ResDcel.
 */
template <typename GeometryTraitsA_2,
          typename GeometryTraitsB_2,
          typename GeometryTraitsRes_2,
          typename TopologyTraitsA,
          typename TopologyTraitsB,
          typename TopologyTraitsRes,
          typename OverlayTraits>
void
overlay(const Arrangement_on_surface_2<GeometryTraitsA_2, TopologyTraitsA>& arr1,
        const Arrangement_on_surface_2<GeometryTraitsB_2, TopologyTraitsB>& arr2,
        Arrangement_on_surface_2<GeometryTraitsRes_2, TopologyTraitsRes>& arr,
        OverlayTraits& ovl_tr)
{
  typedef GeometryTraitsA_2                                     Agt2;
  typedef GeometryTraitsB_2                                     Bgt2;
  typedef GeometryTraitsRes_2                                   Rgt2;
  typedef TopologyTraitsA                                       Att;
  typedef TopologyTraitsB                                       Btt;
  typedef TopologyTraitsRes                                     Rtt;
  typedef OverlayTraits                                         Overlay_traits;

  typedef Arrangement_on_surface_2<Agt2, Att>                   Arr_a;
  typedef Arrangement_on_surface_2<Bgt2, Btt>                   Arr_b;
  typedef Arrangement_on_surface_2<Rgt2, Rtt>                   Arr_res;
  typedef typename Arr_res::Allocator                           Allocator;

  // some type assertions (not all, but better then nothing).
#if !defined(CGAL_NO_ASSERTIONS)
  typedef typename Agt2::Point_2                                A_point;
  typedef typename Bgt2::Point_2                                B_point;
  typedef typename Rgt2::Point_2                                Res_point;
#endif
  CGAL_static_assertion((boost::is_convertible<A_point, Res_point>::value));
  CGAL_static_assertion((boost::is_convertible<B_point, Res_point>::value));

#if !defined(CGAL_NO_ASSERTIONS)
  typedef typename Agt2::X_monotone_curve_2                     A_xcv;
  typedef typename Bgt2::X_monotone_curve_2                     B_xcv;
  typedef typename Rgt2::X_monotone_curve_2                     Res_xcv;
#endif
  CGAL_static_assertion((boost::is_convertible<A_xcv, Res_xcv>::value));
  CGAL_static_assertion((boost::is_convertible<B_xcv, Res_xcv>::value));

  typedef Arr_traits_basic_adaptor_2<Rgt2>              Gt_adaptor_2;
  typedef Arr_overlay_traits_2<Gt_adaptor_2, Arr_a, Arr_b>
                                                        Ovl_gt2;
  typedef Arr_overlay_event<Ovl_gt2, Arr_res, Allocator>
                                                        Ovl_event;
  typedef Arr_overlay_subcurve<Ovl_gt2, Ovl_event, Allocator>
                                                        Ovl_curve;
  typedef typename TopologyTraitsRes::template
    Overlay_helper<Ovl_gt2, Ovl_event, Ovl_curve, Arr_a, Arr_b>
                                                        Ovl_helper;
  typedef Arr_overlay_ss_visitor<Ovl_helper, Overlay_traits>
                                                        Ovl_visitor;

  typedef typename Ovl_gt2::X_monotone_curve_2          Ovl_x_monotone_curve_2;
  typedef typename Ovl_gt2::Point_2                     Ovl_point_2;
  typedef typename Ovl_gt2::Cell_handle_red             Cell_handle_red;
  typedef typename Ovl_gt2::Optional_cell_red           Optional_cell_red;
  typedef typename Ovl_gt2::Cell_handle_blue            Cell_handle_blue;
  typedef typename Ovl_gt2::Optional_cell_blue          Optional_cell_blue;

  CGAL_USE_TYPE(Optional_cell_red);
  CGAL_USE_TYPE(Optional_cell_blue);

  // The result arrangement cannot be on of the input arrangements.
  CGAL_precondition(((void*)(&arr) != (void*)(&arr1)) &&
                    ((void*)(&arr) != (void*)(&arr2)));

  // Prepare a vector of extended x-monotone curves that represent all edges
  // in both input arrangements. Each curve is associated with a halfedge
  // directed from right to left.
  typename Arr_a::Halfedge_const_handle invalid_he1;
  typename Arr_b::Halfedge_const_handle invalid_he2;
  std::vector<Ovl_x_monotone_curve_2>
    xcvs_vec(arr1.number_of_edges() + arr2.number_of_edges());
  unsigned int i = 0;

  typename Arr_a::Edge_const_iterator eit1;
  for (eit1 = arr1.edges_begin(); eit1 != arr1.edges_end(); ++eit1, ++i) {
    typename Arr_a::Halfedge_const_handle he1 = eit1;
    if (he1->direction() != ARR_RIGHT_TO_LEFT) he1 = he1->twin();
    xcvs_vec[i] = Ovl_x_monotone_curve_2(eit1->curve(), he1, invalid_he2);
  }

  typename Arr_b::Edge_const_iterator eit2;
  for (eit2 = arr2.edges_begin(); eit2 != arr2.edges_end(); ++eit2, ++i) {
    typename Arr_b::Halfedge_const_handle he2 = eit2;
    if (he2->direction() != ARR_RIGHT_TO_LEFT) he2 = he2->twin();
    xcvs_vec[i] = Ovl_x_monotone_curve_2(eit2->curve(), invalid_he1, he2);
  }

  // Obtain a extended traits-class object and define the sweep-line visitor.
  const typename Arr_res::Traits_adaptor_2* traits_adaptor =
    arr.traits_adaptor();

  /* We would like to avoid copy construction of the geometry traits class.
   * Copy construction is undesired, because it may results with data
   * duplication or even data loss.
   *
   * If the type Ovl_gt2 is the same as the type
   * GeomTraits, use a reference to GeomTraits to avoid constructing a new one.
   * Otherwise, instantiate a local variable of the former and provide
   * the later as a single parameter to the constructor.
   *
   * Use the form 'A a(*b);' and not ''A a = b;' to handle the case where A has
   * only an implicit constructor, (which takes *b as a parameter).
   */
  typename boost::mpl::if_<boost::is_same<Gt_adaptor_2, Ovl_gt2>,
                           const Ovl_gt2&, Ovl_gt2>::type
    ex_traits(*traits_adaptor);

  Ovl_visitor visitor(&arr1, &arr2, &arr, &ovl_tr);
  Ss2::Surface_sweep_2<Ovl_visitor> surface_sweep(&ex_traits, &visitor);

  // In case both arrangement do not contain isolated vertices, go on and
  // overlay them.
  const std::size_t total_iso_verts =
    arr1.number_of_isolated_vertices() + arr2.number_of_isolated_vertices();

  if (total_iso_verts == 0) {
    // Clear the result arrangement and perform the sweep to construct it.
    arr.clear();
    surface_sweep.sweep(xcvs_vec.begin(), xcvs_vec.end());
    xcvs_vec.clear();
    return;
  }

  // Prepare a vector of extended points that represent all isolated vertices
  // in both input arrangements.
  std::vector<Ovl_point_2> pts_vec(total_iso_verts);

  i = 0;
  typename Arr_a::Vertex_const_iterator  vit1;
  for (vit1 = arr1.vertices_begin(); vit1 != arr1.vertices_end(); ++vit1) {
    if (vit1->is_isolated()) {
      typename Arr_a::Vertex_const_handle v1 = vit1;
      pts_vec[i++] =
        Ovl_point_2(vit1->point(), boost::make_optional(Cell_handle_red(v1)),
                    boost::optional<Cell_handle_blue>());
    }
  }

  typename Arr_b::Vertex_const_iterator  vit2;
  for (vit2 = arr2.vertices_begin(); vit2 != arr2.vertices_end(); ++vit2) {
    if (vit2->is_isolated()) {
      typename Arr_b::Vertex_const_handle v2 = vit2;
      pts_vec[i++] =
        Ovl_point_2(vit2->point(), boost::optional<Cell_handle_red>(),
                    boost::make_optional(Cell_handle_blue(v2)));
    }
  }

  // Clear the result arrangement and perform the sweep to construct it.
  arr.clear();
  surface_sweep.sweep(xcvs_vec.begin(), xcvs_vec.end(),
                      pts_vec.begin(), pts_vec.end());
  xcvs_vec.clear();
  pts_vec.clear();
}

/*! Compute the (simple) overlay of two input arrangements.
 * \param[in] arr1 the first arrangement.
 * \param[in] arr2 the second arrangement.
 * \param[out] arr the resulting arrangement.
 */
template <typename GeometryTraitsA_2,
          typename GeometryTraitsB_2,
          typename GeometryTraitsRes_2,
          typename TopologyTraitsA,
          typename TopologyTraitsB,
          typename TopologyTraitsRes>
void
overlay(const Arrangement_on_surface_2<GeometryTraitsA_2, TopologyTraitsA>& arr1,
        const Arrangement_on_surface_2<GeometryTraitsB_2, TopologyTraitsB>& arr2,
        Arrangement_on_surface_2<GeometryTraitsRes_2, TopologyTraitsRes>& arr)
{
  typedef GeometryTraitsA_2                                     Agt2;
  typedef GeometryTraitsB_2                                     Bgt2;
  typedef GeometryTraitsRes_2                                   Rgt2;
  typedef TopologyTraitsA                                       Att;
  typedef TopologyTraitsB                                       Btt;
  typedef TopologyTraitsRes                                     Rtt;
  typedef Arrangement_on_surface_2<Agt2, Att>                   Arr_a;
  typedef Arrangement_on_surface_2<Bgt2, Btt>                   Arr_b;
  typedef Arrangement_on_surface_2<Rgt2, Rtt>                   Arr_res;

  _Arr_default_overlay_traits_base<Arr_a, Arr_b, Arr_res> ovl_traits;
  overlay(arr1, arr2, arr, ovl_traits);
}

} // namespace CGAL

#include <CGAL/enable_warnings.h>

#endif
