// Copyright (c) 2025 Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s):      Efi Fogel <efifogel@gmail.com>

#ifndef CGAL_ARR_DO_INTERSECT_OVERLAY_2_H
#define CGAL_ARR_DO_INTERSECT_OVERLAY_2_H

#include <CGAL/license/Arrangement_on_surface_2.h>

#include <CGAL/disable_warnings.h>

/*! \file
 *
 * Definition of the global do_intersect_overlay_2() function.
 */

#include <CGAL/Arrangement_on_surface_2.h>
#include <CGAL/Surface_sweep_2.h>
#include <CGAL/Surface_sweep_2/Arr_default_overlay_traits_base.h>
#include <CGAL/Surface_sweep_2/Arr_overlay_traits_2.h>
#include <CGAL/Surface_sweep_2/Arr_do_intersect_overlay_ss_visitor.h>
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
bool do_intersect_overlay(const Arrangement_on_surface_2<GeometryTraitsA_2, TopologyTraitsA>& arr1,
                          const Arrangement_on_surface_2<GeometryTraitsB_2, TopologyTraitsB>& arr2,
                          Arrangement_on_surface_2<GeometryTraitsRes_2, TopologyTraitsRes>& arr,
                          OverlayTraits& ovl_tr,
                          bool ignore_isolated_vertices = true) {
  using Agt2 = GeometryTraitsA_2;
  using Bgt2 = GeometryTraitsB_2;
  using Rgt2 = GeometryTraitsRes_2;
  using Att = TopologyTraitsA;
  using Btt = TopologyTraitsB;
  using Rtt = TopologyTraitsRes;
  using Overlay_traits = OverlayTraits;

  using Arr_a = Arrangement_on_surface_2<Agt2, Att>;
  using Arr_b = Arrangement_on_surface_2<Bgt2, Btt>;
  using Arr_res = Arrangement_on_surface_2<Rgt2, Rtt>;
  using Allocator = typename Arr_res::Allocator;

  // some type assertions (not all, but better than nothing).
  using A_point = typename Agt2::Point_2;
  using B_point = typename Bgt2::Point_2;
  using Res_point = typename Rgt2::Point_2;
  static_assert(std::is_convertible<A_point, Res_point>::value);
  static_assert(std::is_convertible<B_point, Res_point>::value);

  using A_xcv = typename Agt2::X_monotone_curve_2;
  using B_xcv = typename Bgt2::X_monotone_curve_2;
  using Res_xcv = typename Rgt2::X_monotone_curve_2;
  static_assert(std::is_convertible<A_xcv, Res_xcv>::value);
  static_assert(std::is_convertible<B_xcv, Res_xcv>::value);

  using Gt_adaptor_2 = Arr_traits_basic_adaptor_2<Rgt2>;
  using Ovl_gt2 = Arr_overlay_traits_2<Gt_adaptor_2, Arr_a, Arr_b>;
  using Ovl_event = Arr_overlay_event<Ovl_gt2, Arr_res, Allocator>;
  using Ovl_curve = Arr_overlay_subcurve<Ovl_gt2, Ovl_event, Allocator>;
  using Ovl_helper = typename TopologyTraitsRes::template Overlay_helper<Ovl_gt2, Ovl_event, Ovl_curve, Arr_a, Arr_b>;
  using Diovl_visitor = Arr_do_intersect_overlay_ss_visitor<Ovl_helper, Overlay_traits>;

  using Ovl_x_monotone_curve_2 = typename Ovl_gt2::X_monotone_curve_2;
  using Ovl_point_2 = typename Ovl_gt2::Point_2;
  using Cell_handle_red = typename Ovl_gt2::Cell_handle_red;
  using Optional_cell_red = typename Ovl_gt2::Optional_cell_red;
  using Cell_handle_blue = typename Ovl_gt2::Cell_handle_blue;
  using Optional_cell_blue = typename Ovl_gt2::Optional_cell_blue;

  CGAL_USE_TYPE(Optional_cell_red);
  CGAL_USE_TYPE(Optional_cell_blue);

  // The result arrangement cannot be on of the input arrangements.
  CGAL_precondition(((void*)(&arr) != (void*)(&arr1)) && ((void*)(&arr) != (void*)(&arr2)));

  // Prepare a vector of extended x-monotone curves that represent all edges
  // in both input arrangements. Each curve is associated with a halfedge
  // directed from right to left.
  typename Arr_a::Halfedge_const_handle invalid_he1;
  typename Arr_b::Halfedge_const_handle invalid_he2;
  std::vector<Ovl_x_monotone_curve_2> xcvs(arr1.number_of_edges() + arr2.number_of_edges());
  std::size_t i = 0;

  for (auto eit1 = arr1.edges_begin(); eit1 != arr1.edges_end(); ++eit1, ++i) {
    typename Arr_a::Halfedge_const_handle he1 = eit1;
    if (he1->direction() != ARR_RIGHT_TO_LEFT) he1 = he1->twin();
    xcvs[i] = Ovl_x_monotone_curve_2(eit1->curve(), he1, invalid_he2);
  }

  for (auto eit2 = arr2.edges_begin(); eit2 != arr2.edges_end(); ++eit2, ++i) {
    typename Arr_b::Halfedge_const_handle he2 = eit2;
    if (he2->direction() != ARR_RIGHT_TO_LEFT) he2 = he2->twin();
    xcvs[i] = Ovl_x_monotone_curve_2(eit2->curve(), invalid_he1, he2);
  }

  // Obtain an extended traits-class object and define the sweep-line visitor.
  const typename Arr_res::Traits_adaptor_2* traits_adaptor = arr.traits_adaptor();

  /* We would like to avoid copy construction of the geometry traits class.
   * Copy construction is undesired, because it may results with data
   * duplication or even data loss.
   *
   * If the type Ovl_gt2 is the same as the type
   * GeomTraits, use a reference to GeomTraits to avoid constructing a new one.
   * Otherwise, instantiate a local variable of the former and provide
   * the latter as a single parameter to the constructor.
   *
   * Use the form 'A a(*b);' and not ''A a = b;' to handle the case where A has
   * only an implicit constructor, (which takes *b as a parameter).
   */
  std::conditional_t<std::is_same_v<Gt_adaptor_2, Ovl_gt2>, const Ovl_gt2&, Ovl_gt2> ex_traits(*traits_adaptor);

  Diovl_visitor visitor(&arr1, &arr2, &arr, &ovl_tr);
  Ss2::Surface_sweep_2<Diovl_visitor> surface_sweep(&ex_traits, &visitor);

  // In case both arrangement do not contain isolated vertices, go on and overlay them.
  if (ignore_isolated_vertices ||
      ((arr1.number_of_isolated_vertices() == 0) && (arr2.number_of_isolated_vertices() == 0))) {
    // Clear the result arrangement and perform the sweep to construct it.
    arr.clear();
    if (std::is_same<typename Agt2::Bottom_side_category, Arr_contracted_side_tag>::value) {
      surface_sweep.sweep(xcvs.begin(), xcvs.end());
      xcvs.clear();
      return visitor.found_intersection();
    }
    surface_sweep.indexed_sweep(xcvs, Indexed_sweep_accessor<Arr_a, Arr_b, Ovl_x_monotone_curve_2>(arr1, arr2));
    xcvs.clear();
    return visitor.found_intersection();
  }

  // Prepare a vector of extended points that represent all isolated vertices
  // in both input arrangements.
  std::vector<Ovl_point_2> pts_vec(arr1.number_of_isolated_vertices() + arr2.number_of_isolated_vertices());

  i = 0;
  for (auto vit1 = arr1.vertices_begin(); vit1 != arr1.vertices_end(); ++vit1) {
    if (vit1->is_isolated()) {
      typename Arr_a::Vertex_const_handle v1 = vit1;
      pts_vec[i++] = Ovl_point_2(vit1->point(), std::make_optional(Cell_handle_red(v1)),
                                 std::optional<Cell_handle_blue>());
    }
  }

  for (auto vit2 = arr2.vertices_begin(); vit2 != arr2.vertices_end(); ++vit2) {
    if (vit2->is_isolated()) {
      typename Arr_b::Vertex_const_handle v2 = vit2;
      pts_vec[i++] = Ovl_point_2(vit2->point(), std::optional<Cell_handle_red>(),
                                 std::make_optional(Cell_handle_blue(v2)));
    }
  }

  // Clear the result arrangement and perform the sweep to construct it.
  arr.clear();
  if (std::is_same<typename Agt2::Bottom_side_category, Arr_contracted_side_tag>::value) {
    surface_sweep.sweep(xcvs.begin(), xcvs.end(), pts_vec.begin(), pts_vec.end());
    xcvs.clear();
    pts_vec.clear();
    return visitor.found_intersection();
  }
  surface_sweep.indexed_sweep(xcvs, Indexed_sweep_accessor<Arr_a, Arr_b, Ovl_x_monotone_curve_2>(arr1, arr2),
                              pts_vec.begin(), pts_vec.end());
  xcvs.clear();
  pts_vec.clear();
  return visitor.found_intersection();
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
bool do_intersect_overlay(const Arrangement_on_surface_2<GeometryTraitsA_2, TopologyTraitsA>& arr1,
                          const Arrangement_on_surface_2<GeometryTraitsB_2, TopologyTraitsB>& arr2,
                          Arrangement_on_surface_2<GeometryTraitsRes_2, TopologyTraitsRes>& arr) {
  using Agt2 = GeometryTraitsA_2;
  using Bgt2 = GeometryTraitsB_2;
  using Rgt2 = GeometryTraitsRes_2;
  using Att = TopologyTraitsA;
  using Btt = TopologyTraitsB;
  using Rtt = TopologyTraitsRes;
  using Arr_a = Arrangement_on_surface_2<Agt2, Att>;
  using Arr_b = Arrangement_on_surface_2<Bgt2, Btt>;
  using Arr_res = Arrangement_on_surface_2<Rgt2, Rtt>;

  _Arr_default_overlay_traits_base<Arr_a, Arr_b, Arr_res> ovl_traits;
  return do_intersect_overlay(arr1, arr2, arr, ovl_traits);
}

} // namespace CGAL

#include <CGAL/enable_warnings.h>

#endif
