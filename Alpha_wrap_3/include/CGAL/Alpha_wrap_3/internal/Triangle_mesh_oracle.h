// Copyright (c) 2019-2022 Google LLC (USA).
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
#ifndef CGAL_ALPHA_WRAP_3_INTERNAL_TRIANGLE_MESH_ORACLE_H
#define CGAL_ALPHA_WRAP_3_INTERNAL_TRIANGLE_MESH_ORACLE_H

#include <CGAL/license/Alpha_wrap_3.h>

#include <CGAL/Alpha_wrap_3/internal/Alpha_wrap_AABB_geom_traits.h>
#include <CGAL/Alpha_wrap_3/internal/Oracle_base.h>
#include <CGAL/Alpha_wrap_3/internal/splitting_helper.h>

#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/Polygon_mesh_processing/shape_predicates.h>

#include <algorithm>
#include <iostream>
#include <functional>
#include <memory>
#include <type_traits>
#include <vector>

namespace CGAL {
namespace Alpha_wraps_3 {
namespace internal {

// Just some typedefs for readability in the main oracle class
template <typename GT_>
struct TM_oracle_traits
{
  using Geom_traits = Alpha_wrap_AABB_geom_traits<GT_>; // Wrap the kernel to add Ball_3 + custom Do_intersect_3

  using Point_3 = typename Geom_traits::Point_3;
  using AABB_traits = typename AABB_tree_splitter_traits<Point_3, Geom_traits>::AABB_traits;
  using AABB_tree = typename AABB_tree_splitter_traits<Point_3, Geom_traits>::AABB_tree;
};

// @speed could do a partial specialization 'subdivide = false' with simpler code for speed?
template <typename GT_,
          typename BaseOracle = int,
          bool subdivide = true>
class Triangle_mesh_oracle
  : // this is the base that handles calls to the AABB tree
    public AABB_tree_oracle<typename TM_oracle_traits<GT_>::Geom_traits,
                            typename TM_oracle_traits<GT_>::AABB_tree,
                            typename std::conditional<
                                      /*condition*/subdivide,
                                      /*true*/Splitter_traversal_traits<typename TM_oracle_traits<GT_>::AABB_traits>,
                                      /*false*/Default_traversal_traits<typename TM_oracle_traits<GT_>::AABB_traits> >::type,
                            BaseOracle>,
    // this is the base that handles splitting input faces and inserting them into the AABB tree
    public AABB_tree_oracle_splitter<subdivide,
                                     typename TM_oracle_traits<GT_>::Point_3,
                                     typename TM_oracle_traits<GT_>::Geom_traits>
{
  using TMOT = TM_oracle_traits<GT_>;
  using Base_GT = GT_;

public:
  using Geom_traits = typename TMOT::Geom_traits;

private:
  using Point_3 = typename Geom_traits::Point_3;
  using Triangle_3 = typename Geom_traits::Triangle_3;

  using AABB_traits = typename TMOT::AABB_traits;
  using AABB_tree = typename TMOT::AABB_tree;
  using AABB_traversal_traits = typename std::conditional<
                                  /*condition*/subdivide,
                                  /*true*/Splitter_traversal_traits<AABB_traits>,
                                  /*false*/Default_traversal_traits<AABB_traits> >::type;

  using Oracle_base = AABB_tree_oracle<Geom_traits, AABB_tree, AABB_traversal_traits, BaseOracle>;
  using Splitter_base = AABB_tree_oracle_splitter<subdivide, Point_3, Geom_traits>;

public:
  // Constructors
  //
  // When using this constructor (and thus doing actual splitting), note that the oracle
  // will be adapted to this particular 'alpha', and so when calling again AW3(other_alpha)
  // the oracle might not have performed a split that is adapted to this other alpha value.
  Triangle_mesh_oracle(const double alpha,
                       const BaseOracle& base_oracle = BaseOracle(),
                       const Base_GT& gt = Base_GT())
    : Oracle_base(base_oracle, gt), Splitter_base(alpha)
  {
    Splitter_base::initialize_tree_property_maps(this->tree());
  }

  Triangle_mesh_oracle(const double alpha,
                       const Base_GT& gt,
                       const BaseOracle& base_oracle = BaseOracle())
    : Triangle_mesh_oracle(alpha, base_oracle, gt)
  { }

 Triangle_mesh_oracle(const BaseOracle& base_oracle,
                      const Base_GT& gt = Base_GT())
   : Triangle_mesh_oracle(0. /*alpha*/, base_oracle, gt)
 { }

 Triangle_mesh_oracle(const Base_GT& gt,
                      const BaseOracle& base_oracle = BaseOracle())
   : Triangle_mesh_oracle(0. /*alpha*/, base_oracle, gt)
 { }

 Triangle_mesh_oracle()
   : Triangle_mesh_oracle(0. /*alpha*/, BaseOracle(), Base_GT())
 { }

public:
  template <typename TriangleMesh,
            typename CGAL_NP_TEMPLATE_PARAMETERS>
  void add_triangle_mesh(const TriangleMesh& tmesh,
                         const CGAL_NP_CLASS& np = CGAL::parameters::default_values())
  {
    using parameters::get_parameter;
    using parameters::choose_parameter;

    using face_descriptor = typename boost::graph_traits<TriangleMesh>::face_descriptor;

    using VPM = typename GetVertexPointMap<TriangleMesh>::const_type;
    using Point_ref = typename boost::property_traits<VPM>::reference;

    CGAL_precondition(CGAL::is_triangle_mesh(tmesh));

    if(is_empty(tmesh))
    {
#ifdef CGAL_AW3_DEBUG
      std::cout << "Warning: Input is empty (TM)" << std::endl;
#endif
      return;
    }

#ifdef CGAL_AW3_DEBUG
    std::cout << "Insert into AABB tree (faces)..." << std::endl;
#endif

    VPM vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                               get_const_property_map(vertex_point, tmesh));
    static_assert(std::is_same<typename boost::property_traits<VPM>::value_type, Point_3>::value);

    Splitter_base::reserve(num_faces(tmesh));

    for(face_descriptor f : faces(tmesh))
    {
      if(Polygon_mesh_processing::is_degenerate_triangle_face(f, tmesh, np))
      {
#ifdef CGAL_AW3_DEBUG
        std::cerr << "Warning: ignoring degenerate face " << f << std::endl;
#endif
        continue;
      }

      const Point_ref p0 = get(vpm, source(halfedge(f, tmesh), tmesh));
      const Point_ref p1 = get(vpm, target(halfedge(f, tmesh), tmesh));
      const Point_ref p2 = get(vpm, target(next(halfedge(f, tmesh), tmesh), tmesh));

      const Triangle_3 tr = this->geom_traits().construct_triangle_3_object()(p0, p1, p2);

      Splitter_base::split_and_insert_datum(tr, this->tree(), this->geom_traits());
    }

    // Manually constructing it here purely for profiling reasons: if we keep the lazy approach,
    // it will be done at the first treatment of a facet that needs a Steiner point.
    // So if one wanted to bench the flood fill runtime, it would be skewed by the time it takes
    // to accelerate the tree.
    this->tree().accelerate_distance_queries();

#ifdef CGAL_AW3_DEBUG
    std::cout << "Tree: " << this->tree().size() << " primitives (" << num_faces(tmesh) << " faces in input)" << std::endl;
#endif
  }
};

} // namespace internal
} // namespace Alpha_wraps_3
} // namespace CGAL

#endif // CGAL_ALPHA_WRAP_3_INTERNAL_TRIANGLE_MESH_ORACLE_H
