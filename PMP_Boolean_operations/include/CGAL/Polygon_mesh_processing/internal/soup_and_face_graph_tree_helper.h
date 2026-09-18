// Copyright (c) 2016 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sebastien Loriot

#ifndef CGAL_POLYGON_MESH_PROCESSING_INTERNAL_SOUP_AND_TREE_HELPER_H
#define CGAL_POLYGON_MESH_PROCESSING_INTERNAL_SOUP_AND_TREE_HELPER_H

#include <CGAL/license/Polygon_mesh_processing/corefinement.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/AABB_trees/intersection.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits_3.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_indexed_triangle_primitive_3.h>

#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/property_map.h>
#include <fstream>
#include <sstream>
#include <set>
#include <type_traits>

namespace CGAL {
namespace Polygon_mesh_processing {
namespace internal {

template<class TriangleMesh, class GT, class VPM>
struct AABB_tree_graph_helper
{
  using Primitive = AABB_face_graph_triangle_primitive<TriangleMesh, VPM>;
  using Traits = AABB_traits_3<GT, Primitive>;
  using Tree = AABB_tree<Traits>;

  using Graph_traits = boost::graph_traits<TriangleMesh>;
  using vertex_descriptor = typename Graph_traits::vertex_descriptor;
  using face_descriptor = typename Graph_traits::face_descriptor;

  template <class RPM>
  struct Split_primitives
  {
    Split_primitives(const RPM &rpm): rpm(rpm){}

    template<typename PrimitiveIterator>
    void operator()(PrimitiveIterator first,
                    PrimitiveIterator beyond,
                    const CGAL::Bbox_3& bbox) const
      {
        PrimitiveIterator middle = first + (beyond - first)/2;
        const int crd = bbox.largest_span_index();
        std::nth_element(first, middle, beyond, [this, crd](const Primitive& p1, const Primitive& p2){ return get(rpm, p1.id())[crd] < get(rpm, p2.id())[crd];});
      }
    const RPM &rpm;
  };

  // For exact side_of_triangle_mesh
  template <class BPM>
  struct Compute_bbox {
    Compute_bbox(const BPM& bpm): bpm(bpm){}

    template<typename ConstPrimitiveIterator>
    CGAL::Bbox_3 operator()(ConstPrimitiveIterator first,
                            ConstPrimitiveIterator beyond) const
    {
      CGAL::Bbox_3 bbox = get(bpm, first->id());
      for(++first; first != beyond; ++first)
        bbox += get(bpm, first->id());
      return bbox;
    }
    BPM bpm;
  };

  template<class Concurrency_tag=Sequential_tag, class FaceRange, class VertexPointMap>
  void build(const FaceRange &face_range, Tree& tree, const TriangleMesh& tm, VertexPointMap& vpm){
    using Face_bbox_tag = typename CGAL::dynamic_face_property_t<Bbox_3>;
    using Face_ref_point_tag = typename CGAL::dynamic_face_property_t<Epick::Point_3>;
    using Bbox_map = typename boost::property_map<TriangleMesh, Face_bbox_tag>::const_type;
    using Ref_point_map = typename boost::property_map<TriangleMesh, Face_ref_point_tag>::const_type;

    using VPM_kernel = typename Kernel_traits<typename boost::property_traits<VertexPointMap>::value_type>::Kernel;
    CGAL::Cartesian_converter<VPM_kernel, Epick> to_input;

    Bbox_map bb_map = get(Face_bbox_tag(), tm);
    Ref_point_map rp_map = get(Face_ref_point_tag(), tm);

#ifdef CGAL_LINKED_WITH_TBB
    if constexpr(std::is_same_v<Concurrency_tag, Parallel_tag>)
    {
      tbb::parallel_for(std::size_t(0), faces(tm).size(), [&](std::size_t i){
        face_descriptor f(i);
        put(bb_map, f, face_bbox(f, tm));
        put(rp_map, f, to_input(get(vpm, target(halfedge(f, tm), tm))) );
      });
    }
    else
#endif
    {
      for(face_descriptor f : faces(tm)){
        put(bb_map, f, face_bbox(f, tm));
        put(rp_map, f, to_input(get(vpm, target(halfedge(f, tm), tm))) );
      }
    }
    tree.insert(face_range.begin(), face_range.end(), tm, vpm);

    Compute_bbox<Bbox_map> compute_bbox(bb_map);
    Split_primitives<Ref_point_map> split_primitives(rp_map);
    tree.template custom_build<Concurrency_tag>(compute_bbox, split_primitives);
  }

  template<class ConcurrencyTag=Sequential_tag, class VertexPointMap>
  void build(Tree& tree, const TriangleMesh& tm, VertexPointMap& vpm){
    build<ConcurrencyTag>(faces(tm), tree, tm, vpm);
  }
};

#ifndef DOXYGEN_RUNNING
template <class PointRange, class VPM>
struct Property_map_for_soup
{
  using VPM_base   = VPM;
  using key_type   = std::size_t;
  using value_type = typename boost::property_traits<VPM>::value_type;
  using category   = boost::readable_property_map_tag;
  using reference  = typename boost::property_traits<VPM>::reference;

  const PointRange& points;
  VPM vpm;

  Property_map_for_soup(const PointRange& points, VPM vpm)
    : points(points)
    , vpm(vpm)
  {}

  inline friend
  reference get(const Property_map_for_soup<PointRange, VPM>& map, key_type k)
  {
    return get(map.vpm, map.points[k]);
  }
};

template<class PointRange, class TriangleRange, class GT, class VPM>
struct AABB_tree_soup_helper
{
  using Primitive = AABB_indexed_triangle_primitive_3<GT, PointRange, TriangleRange, Tag_false, typename VPM::VPM_base>;
  using Traits = AABB_traits_3<GT, Primitive>;
  using Tree = AABB_tree<Traits>;

  template <class RPM>
  struct Split_primitives
  {
    Split_primitives(const RPM &rpm): rpm(rpm){}

    template<typename PrimitiveIterator>
    void operator()(PrimitiveIterator first,
                    PrimitiveIterator beyond,
                    const CGAL::Bbox_3& bbox) const
      {
        PrimitiveIterator middle = first + (beyond - first)/2;
        const int crd = bbox.largest_span_index();
        std::nth_element(first, middle, beyond, [this, crd](const auto& p1, const auto& p2){ return get(rpm, p1.id())[crd] < get(rpm, p2.id())[crd];});
      }
    const RPM &rpm;
  };

  // For exact side_of_triangle_mesh
  template <class BPM>
  struct Compute_bbox {
    Compute_bbox(const BPM& bpm): bpm(bpm){}

    template<typename ConstPrimitiveIterator>
    CGAL::Bbox_3 operator()(ConstPrimitiveIterator first,
                            ConstPrimitiveIterator beyond) const
    {
      CGAL::Bbox_3 bbox = get(bpm, first->id());
      for(++first; first != beyond; ++first)
        bbox += get(bpm, first->id());
      return bbox;
    }
    BPM bpm;
  };

  template<class ConcurrencyTag=Sequential_tag, class FaceRange>
  void build(const FaceRange &face_range, Tree& tree, const PointRange& pts, const TriangleRange& triangles, VPM vpm){
    using PM_kernel     = typename Kernel_traits<typename boost::property_traits<VPM>::value_type>::Kernel;
    using Bbox_map      = Pointer_property_map<Bbox_3>::type;
    using Ref_point_map = Pointer_property_map<Epick::Point_3>::type;

    CGAL::Cartesian_converter<PM_kernel, Epick> to_input;

    tree.insert(face_range.begin(), face_range.end(), pts, triangles, vpm.vpm);

    std::vector<Bbox_3> bb_vector(triangles.size(), Bbox_3());
    std::vector<Epick::Point_3> rp_vector(triangles.size(), Epick::Point_3());
    Bbox_map bb_map = make_property_map(bb_vector);
    Ref_point_map rp_map = make_property_map(rp_vector);

#ifdef CGAL_LINKED_WITH_TBB
    if constexpr(std::is_same_v<ConcurrencyTag, Parallel_tag>)
    {
      tbb::parallel_for(std::size_t(0), face_range.size(), [&](std::size_t i){
        i = face_range[i];
        put(bb_map, i, get(vpm, triangles[i][0]).bbox() +
                       get(vpm, triangles[i][1]).bbox() +
                       get(vpm, triangles[i][2]).bbox());
        put(rp_map, i, to_input(get(vpm, triangles[i][0])) );
      });
    }
    else
#endif
    {
      for(std::size_t i: face_range){
        put(bb_map, i, get(vpm, triangles[i][0]).bbox() +
                       get(vpm, triangles[i][1]).bbox() +
                       get(vpm, triangles[i][2]).bbox());
        put(rp_map, i, to_input(get(vpm, triangles[i][0])) );
      }
    }
    Compute_bbox<Bbox_map> compute_bbox(bb_map);
    Split_primitives<Ref_point_map> split_primitives(rp_map);
    tree.template custom_build<ConcurrencyTag>(compute_bbox, split_primitives);
  }

  template<class ConcurrencyTag=Sequential_tag>
  void build(Tree& tree, const PointRange& pts, const TriangleRange& triangles, VPM vpm){
    build<ConcurrencyTag>(triangles, tree, pts, triangles, vpm);
  }
};

#endif

} } } // CGAL::Polygon_mesh_processing::internal

#endif // CGAL_POLYGON_MESH_PROCESSING_INTERNAL_SOUP_AND_TREE_HELPER_H
