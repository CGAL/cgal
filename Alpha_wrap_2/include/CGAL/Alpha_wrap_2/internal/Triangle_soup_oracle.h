// Copyright (c) 2019-2022 Google LLC (USA).
// Copyright (c) 2025 GeometryFactory (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé
//
#ifndef CGAL_ALPHA_WRAP_2_INTERNAL_TRIANGLE_SOUP_ORACLE_H
#define CGAL_ALPHA_WRAP_2_INTERNAL_TRIANGLE_SOUP_ORACLE_H

#include <CGAL/license/Alpha_wrap_2.h>

#include <CGAL/Alpha_wrap_2/internal/Alpha_wrap_AABB_geom_traits.h>
#include <CGAL/Alpha_wrap_2/internal/Oracle_base.h>
#include <CGAL/Alpha_wrap_2/internal/splitting_helper.h>

#include <CGAL/AABB_triangle_primitive_2.h>
#include <CGAL/AABB_traits_2.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/Named_function_parameters.h>

#include <boost/range/value_type.hpp>

#include <algorithm>
#include <iostream>
#include <functional>

namespace CGAL {
namespace Alpha_wraps_2 {
namespace internal {

// Just some typedefs for readability
template <typename GT_>
struct TS_oracle_traits
{
  using Geom_traits = Alpha_wrap_AABB_geom_traits<GT_>; // Wrap the kernel to add Disk_2 + custom Do_intersect_2
  using Point_2 = typename Geom_traits::Point_2;
  using AABB_traits = typename AABB_tree_splitter_traits<Point_2, Geom_traits>::AABB_traits;
  using AABB_tree = typename AABB_tree_splitter_traits<Point_2, Geom_traits>::AABB_tree;
};

template <typename GT_,
          typename BaseOracle = int,
          bool subdivide = true>
class Triangle_soup_oracle
  : // this is the base that handles calls to the AABB tree
    public AABB_tree_oracle<typename TS_oracle_traits<GT_>::Geom_traits,
                            typename TS_oracle_traits<GT_>::AABB_tree,
                            typename std::conditional<
                              /*condition*/subdivide,
                              /*true*/Splitter_traversal_traits<typename TS_oracle_traits<GT_>::AABB_traits>,
                              /*false*/Default_traversal_traits<typename TS_oracle_traits<GT_>::AABB_traits> >::type,
                            BaseOracle>,
    // this is the base that handles splitting input faces and inserting them into the AABB tree
    public AABB_tree_oracle_splitter<subdivide,
                                     typename TS_oracle_traits<GT_>::Point_2,
                                     typename TS_oracle_traits<GT_>::Geom_traits>
{
  using TSOT = TS_oracle_traits<GT_>;
  using Base_GT = GT_;

public:
  using Geom_traits = typename TSOT::Geom_traits;

private:
  using Point_2 = typename Geom_traits::Point_2;
  using Triangle_2 = typename Geom_traits::Triangle_2;

  using AABB_traits = typename TSOT::AABB_traits;
  using AABB_tree = typename TSOT::AABB_tree;
  using AABB_traversal_traits = typename std::conditional<
                                  /*condition*/subdivide,
                                  /*true*/Splitter_traversal_traits<AABB_traits>,
                                  /*false*/Default_traversal_traits<AABB_traits> >::type;

  using Oracle_base = AABB_tree_oracle<Geom_traits, AABB_tree, AABB_traversal_traits, BaseOracle>;
  using Splitter_base = AABB_tree_oracle_splitter<subdivide, Point_2, Geom_traits>;

public:
  // Constructors
  //
  // When using this constructor (and thus doing actual splitting), note that the oracle
  // will be adapted to this particular 'alpha', and so when calling again AW2(other_alpha)
  // the oracle might not have performed a split that is adapted to this other alpha value.
  Triangle_soup_oracle(const double alpha,
                       const BaseOracle& base_oracle = BaseOracle(),
                       const Base_GT& gt = Base_GT())
    : Oracle_base(base_oracle, gt), Splitter_base(alpha)
  {
    Splitter_base::initialize_tree_property_maps(this->tree());
  }

  Triangle_soup_oracle(const double alpha,
                       const Base_GT& gt,
                       const BaseOracle& base_oracle = BaseOracle())
    : Triangle_soup_oracle(alpha, base_oracle, gt)
  { }

 Triangle_soup_oracle(const BaseOracle& base_oracle,
                      const Base_GT& gt = Base_GT())
   : Triangle_soup_oracle(0. /*alpha*/, base_oracle, gt)
 { }

 Triangle_soup_oracle(const Base_GT& gt,
                      const BaseOracle& base_oracle = BaseOracle())
   : Triangle_soup_oracle(0. /*alpha*/, base_oracle, gt)
 { }

 Triangle_soup_oracle()
   : Triangle_soup_oracle(0. /*alpha*/, BaseOracle(), Base_GT())
 { }

public:
  template <typename PointRange, typename FaceRange,
            typename CGAL_NP_TEMPLATE_PARAMETERS>
  void add_triangle_soup(const PointRange& points,
                         const FaceRange& faces,
                         const CGAL_NP_CLASS& np = CGAL::parameters::default_values())
  {
    using parameters::choose_parameter;
    using parameters::get_parameter;

    using PPM = typename GetPointMap<PointRange, CGAL_NP_CLASS>::const_type;
    using Point_ref = typename boost::property_traits<PPM>::reference;

    using Face = typename boost::range_value<FaceRange>::type;

    if(points.empty() || faces.empty())
    {
#ifdef CGAL_AW2_DEBUG
      std::cout << "Warning: Input is empty (TS)" << std::endl;
#endif
      return;
    }

#ifdef CGAL_AW2_DEBUG
    std::cout << "Insert into AABB Tree (2D triangles)..." << std::endl;
#endif

    PPM pm = choose_parameter<PPM>(get_parameter(np, internal_np::point_map));
    static_assert(std::is_same<typename boost::property_traits<PPM>::value_type, Point_2>::value);

    Splitter_base::reserve(faces.size());

    typename Geom_traits::Construct_triangle_2 triangle = this->geom_traits().construct_triangle_2_object();
    typename Geom_traits::Is_degenerate_2 is_degenerate = this->geom_traits().is_degenerate_2_object();

    for(const Face& f : faces)
    {
      CGAL_precondition(std::distance(std::cbegin(f), std::cend(f)) == 3);

      auto vi = std::cbegin(f);
      CGAL_assertion(*vi < points.size());
      Point_ref p0 = get(pm, points[*vi++]);
      CGAL_assertion(*vi < points.size());
      Point_ref p1 = get(pm, points[*vi++]);
      CGAL_assertion(*vi < points.size());
      Point_ref p2 = get(pm, points[*vi]);

      const Triangle_2 tr = triangle(p0, p1, p2);
      if(is_degenerate(tr))
      {
#ifdef CGAL_AW2_DEBUG
        std::cerr << "Warning: ignoring degenerate face " << tr << std::endl;
#endif
        continue;
      }

      Splitter_base::split_and_insert_datum(tr, this->tree(), this->geom_traits());
    }

    // Manually constructing it here purely for profiling reasons: if we keep the lazy approach,
    // it will be done at the first treatment of an edge that needs a Steiner point.
    // So if one wanted to bench the flood fill runtime, it would be skewed by the time it takes
    // to accelerate the tree.
    this->tree().accelerate_distance_queries();

#ifdef CGAL_AW2_DEBUG
    std::cout << "Tree: " << this->tree().size() << " primitives (" << faces.size() << " faces in input)" << std::endl;
#endif
  }

  template <typename TriangleRange,
            typename CGAL_NP_TEMPLATE_PARAMETERS>
  void add_triangle_soup(const TriangleRange& triangles,
                         const CGAL_NP_CLASS& /*np*/ = CGAL::parameters::default_values())
  {
    if(triangles.empty())
    {
#ifdef CGAL_AW2_DEBUG
      std::cout << "Warning: Input is empty (TS)" << std::endl;
#endif
      return;
    }

#ifdef CGAL_AW2_DEBUG
    std::cout << "Insert into AABB Tree (triangles)..." << std::endl;
#endif

    typename Geom_traits::Is_degenerate_2 is_degenerate = this->geom_traits().is_degenerate_2_object();

    Splitter_base::reserve(triangles.size());

    for(const Triangle_2& tr : triangles)
    {
      if(is_degenerate(tr))
      {
#ifdef CGAL_AW2_DEBUG
        std::cerr << "Warning: ignoring degenerate triangle " << tr << std::endl;
#endif
        continue;
      }

      Splitter_base::split_and_insert_datum(tr, this->tree(), this->geom_traits());
    }

#ifdef CGAL_AW2_DEBUG
    std::cout << "Tree: " << this->tree().size() << " primitives (" << triangles.size() << " triangles in input)" << std::endl;
#endif
  }
};

} // namespace internal
} // namespace Alpha_wraps_2
} // namespace CGAL

#endif // CGAL_ALPHA_WRAP_2_INTERNAL_TRIANGLE_SOUP_ORACLE_H
