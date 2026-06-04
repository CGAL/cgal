// Copyright (c) 2026  GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s) : Léo Valque

#ifndef CGAL_AABB_MESHES_INTERSECTIONS_H
#define CGAL_AABB_MESHES_INTERSECTIONS_H

#include <memory>

#include <CGAL/license/AABB_tree.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits_3.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

#include <CGAL/AABB_tree/internal/AABB_two_tree_traversal.h>

#include <CGAL/Polygon_mesh_processing/self_intersections.h>

#include <CGAL/box_intersection_d.h>
#include <CGAL/Box_intersection_d/Box_with_info_d.h>

#ifdef CGAL_LINKED_WITH_TBB
#include <tbb/tbb.h>
#include <CGAL/mutex.h>
#endif

namespace CGAL {

namespace Polygon_mesh_processing{

namespace experimental{

template< typename TriangleMesh1,
          typename TriangleMesh2,
          typename OutputIterator,
          typename NamedParameters1 = parameters::Default_named_parameters,
          typename NamedParameters2 = parameters::Default_named_parameters>
void mixed_meshes_intersections(const TriangleMesh1 &tm1,
                                        const TriangleMesh2 &tm2,
                                        OutputIterator out,
                                        const NamedParameters1& np1 = parameters::default_values(),
                                        const NamedParameters2& np2 = parameters::default_values())
{
  using CGAL::parameters::choose_parameter;
  using CGAL::parameters::get_parameter;

  using face_descriptor_1 = typename boost::graph_traits<TriangleMesh1>::face_descriptor;
  using halfedge_descriptor_1 = typename boost::graph_traits<TriangleMesh1>::halfedge_descriptor;

  using face_descriptor_2 = typename boost::graph_traits<TriangleMesh2>::face_descriptor;
  using halfedge_descriptor_2 = typename boost::graph_traits<TriangleMesh2>::halfedge_descriptor;

  using Concurrency_tag = typename internal_np::Lookup_named_param_def <
                                                internal_np::concurrency_tag_t,
                                                NamedParameters1,
                                                Sequential_tag
                                              > ::type;
  using GT = typename GetGeomTraits<TriangleMesh1, NamedParameters1>::type;
  GT gt = choose_parameter<GT>(get_parameter(np1, internal_np::geom_traits));

  using Face_bbox_tag = typename CGAL::dynamic_face_property_t<Bbox_3>;
  using Primitive_1 = AABB_face_graph_triangle_primitive<TriangleMesh1>;
  using Bbox_pmap_1 = typename boost::property_map<TriangleMesh1, Face_bbox_tag>::const_type;
  using Traits_1 = AABB_traits_3<GT, Primitive_1, Bbox_pmap_1>;
  using Tree_1 = AABB_tree<Traits_1>;
  using Node_1 = AABB_node<Traits_1>;
  using ConstPrimitiveIterator_1 = typename std::vector< Primitive_1 >::const_iterator;
  using Box_1 = Box_intersection_d::Box_with_info_d<double, 3, face_descriptor_1>;

  using Primitive_2 = AABB_face_graph_triangle_primitive<TriangleMesh2>;
  using Bbox_pmap_2 = typename boost::property_map<TriangleMesh2, Face_bbox_tag>::const_type;
  using Traits_2 = AABB_traits_3<GT, Primitive_2, Bbox_pmap_2>;
  using Tree_2 = AABB_tree<Traits_2>;
  using Node_2 = AABB_node<Traits_2>;
  using ConstPrimitiveIterator_2 = typename std::vector< Primitive_2 >::const_iterator;
  using Box_2 = Box_intersection_d::Box_with_info_d<double, 3, face_descriptor_2>;

  using InternOutputIterator= std::back_insert_iterator<std::vector<std::pair<Node_1*, Node_2*>>>;

  const std::size_t cutoff = 50000;

  auto vpm1 = choose_parameter(get_parameter(np1, internal_np::vertex_point),
                               get_const_property_map(vertex_point, tm1));
  auto vpm2 = choose_parameter(get_parameter(np2, internal_np::vertex_point),
                               get_const_property_map(vertex_point, tm2));

  auto triangle = [&](auto fd, const auto &vpm, const auto &tm){
    auto hd = halfedge(fd,tm);
    auto a = get(vpm, source(hd,tm));
    auto b = get(vpm, target(hd,tm));
    auto c = get(vpm, target(next(hd,tm),tm));
    return typename GT::Triangle_3(a, b, c);
  };

  auto bbox = [](auto fd, const auto &vpm, const auto &tm){
    auto hd = halfedge(fd,tm);
    Bbox_3 res = get(vpm, source(hd,tm)).bbox();
    res += get(vpm, target(hd,tm)).bbox();
    res += get(vpm, target(next(hd,tm),tm)).bbox();
    return res;
  };

  Bbox_pmap_1 bb1 = get(Face_bbox_tag(), tm1);
  Bbox_pmap_2 bb2 = get(Face_bbox_tag(), tm2);
#ifdef CGAL_LINKED_WITH_TBB
  if constexpr(std::is_same_v<Concurrency_tag, Parallel_tag>)
  {
    oneapi::tbb::parallel_for(
      oneapi::tbb::blocked_range<size_t>(0, faces(tm1).size()),
        [&](const oneapi::tbb::blocked_range<size_t>& r) {
          for (size_t i = r.begin(); i < r.end(); ++i) {
            face_descriptor_1 fd = *(faces(tm1).begin() + i);
            put(bb1, fd, bbox(fd, vpm1, tm1));
          }
        }
    );

    oneapi::tbb::parallel_for(
      oneapi::tbb::blocked_range<size_t>(0, faces(tm2).size()),
        [&](const oneapi::tbb::blocked_range<size_t>& r) {
          for (size_t i = r.begin(); i < r.end(); ++i) {
            face_descriptor_2 fd = *(faces(tm2).begin() + i);
            put(bb2, fd, bbox(fd, vpm2, tm2));
          }
        }
    );
  }
  else
#endif
  {
    for(face_descriptor_1 fd : faces(tm1))
      put(bb1, fd, bbox(fd, vpm1, tm1));
    for(face_descriptor_2 fd : faces(tm2))
      put(bb2, fd, bbox(fd, vpm2, tm2));
  }

  Traits_1 traits1(bb1);
  Tree_1 tree1(traits1);
  tree1.insert(faces(tm1).first, faces(tm1).second, tm1);

  Traits_2 traits2(bb2);
  Tree_2 tree2(traits2);
  tree2.insert(faces(tm2).first, faces(tm2).second, tm2);

#ifdef CGAL_LINKED_WITH_TBB
  if constexpr(std::is_same_v<Concurrency_tag, Parallel_tag>)
  {
    CGAL_MUTEX m;
    oneapi::tbb::task_group tg;
    tg.run([&]{ tree1.template partial_build<Concurrency_tag>(cutoff); });
    tree2.template partial_build<Concurrency_tag>(cutoff);
    tg.wait();

    tbb::concurrent_vector<std::pair<const Node_1*, const Node_2*>> inter;
    CGAL::internal::AABB_tree::two_tree_listing_intersecting_patches<Concurrency_tag>(tree1, tree2, std::back_inserter(inter), cutoff, bb1, bb2);
    oneapi::tbb::parallel_for(
      oneapi::tbb::blocked_range<size_t>(0, inter.size()),
      [&](const oneapi::tbb::blocked_range<size_t>& r) {
        for (size_t i = r.begin(); i != r.end(); ++i) {
          const auto [begin_1, end_1] = tree1.partial_node_to_primitives_iterator(*inter[i].first);
          const auto [begin_2, end_2] = tree2.partial_node_to_primitives_iterator(*inter[i].second);
          std::vector< Box_1 > boxes_1;
          std::vector< Box_2 > boxes_2;
          std::vector< Box_1* > boxes_ptr_1;
          std::vector< Box_2* > boxes_ptr_2;
          boxes_1.reserve(std::distance(begin_1, end_1));
          boxes_2.reserve(std::distance(begin_2, end_2));
          boxes_ptr_1.reserve(boxes_1.size());
          boxes_ptr_2.reserve(boxes_2.size());

          for(auto it=begin_1; it!=end_1; ++it)
            boxes_1.emplace_back(get(bb1, it->id()), it->id());
          for(auto it=begin_2; it!=end_2; ++it)
            boxes_2.emplace_back(get(bb2, it->id()), it->id());
          for(auto &b: boxes_1)
            boxes_ptr_1.push_back(std::addressof(b));
          for(auto &b: boxes_2)
            boxes_ptr_2.push_back(std::addressof(b));

          std::vector< std::pair<face_descriptor_1, face_descriptor_2>> patch_out;
          const std::ptrdiff_t cutoff = 2000;
          box_intersection_d(boxes_ptr_1.begin(), boxes_ptr_1.end(), boxes_ptr_2.begin(), boxes_ptr_2.end(),
              [&](const Box_1 *b1, const Box_2 *b2)
              {
                if(CGAL::do_intersect(triangle(b1->info(), vpm1, tm1), triangle(b2->info(), vpm2, tm2)))
                  patch_out.emplace_back(b1->info(), b2->info());
              }, cutoff);
              CGAL_SCOPED_LOCK(m);
          for(auto p: patch_out)
            *out ++ = p;
        }
      }
    );
  }
  else
#endif
  {
    tree1.template partial_build<Concurrency_tag>(cutoff);
    tree2.template partial_build<Concurrency_tag>(cutoff);

    std::vector<std::pair<const Node_1*, const Node_2*>> inter;
    CGAL::internal::AABB_tree::two_tree_listing_intersecting_patches(tree1, tree2, std::back_inserter(inter), cutoff, bb1, bb2);

    for(const auto& [n_1, n_2]: inter){
      const auto [begin_1, end_1] = tree1.partial_node_to_primitives_iterator(*n_1);
      const auto [begin_2, end_2] = tree2.partial_node_to_primitives_iterator(*n_2);
      std::vector< Box_1 > boxes_1;
      std::vector< Box_2 > boxes_2;
      std::vector< Box_1* > boxes_ptr_1;
      std::vector< Box_2* > boxes_ptr_2;
      boxes_1.reserve(std::distance(begin_1, end_1));
      boxes_2.reserve(std::distance(begin_2, end_2));
      boxes_ptr_1.reserve(boxes_1.size());
      boxes_ptr_2.reserve(boxes_2.size());

      for(auto it=begin_1; it!=end_1; ++it)
        boxes_1.emplace_back(get(bb1, it->id()), it->id());
      for(auto it=begin_2; it!=end_2; ++it)
        boxes_2.emplace_back(get(bb2, it->id()), it->id());
      for(auto &b: boxes_1)
        boxes_ptr_1.push_back(std::addressof(b));
      for(auto &b: boxes_2)
        boxes_ptr_2.push_back(std::addressof(b));

      std::vector< std::pair<face_descriptor_1, face_descriptor_2>> patch_out;
      const std::ptrdiff_t cutoff = 2000;
      box_intersection_d(boxes_ptr_1.begin(), boxes_ptr_1.end(), boxes_ptr_2.begin(), boxes_ptr_2.end(),
          [&](const Box_1 *b1, const Box_2 *b2)
          {
            if(CGAL::do_intersect(triangle(b1->info(), vpm1, tm1), triangle(b2->info(), vpm2, tm2)))
              patch_out.emplace_back(b1->info(), b2->info());
          }, cutoff);
      for(auto p: patch_out)
        *out ++ = p;
    }
  }
}

template< typename TriangleMesh1,
          typename TriangleMesh2,
          typename OutputIterator,
          typename NamedParameters1 = parameters::Default_named_parameters,
          typename NamedParameters2 = parameters::Default_named_parameters>
void box_meshes_intersections(const TriangleMesh1 &tm1,
                                        const TriangleMesh2 &tm2,
                                        OutputIterator out,
                                        const NamedParameters1& np1 = parameters::default_values(),
                                        const NamedParameters2& np2 = parameters::default_values())
{
  using CGAL::parameters::choose_parameter;
  using CGAL::parameters::get_parameter;

  using face_descriptor_1 = typename boost::graph_traits<TriangleMesh1>::face_descriptor;
  using halfedge_descriptor_1 = typename boost::graph_traits<TriangleMesh1>::halfedge_descriptor;

  using face_descriptor_2 = typename boost::graph_traits<TriangleMesh2>::face_descriptor;
  using halfedge_descriptor_2 = typename boost::graph_traits<TriangleMesh2>::halfedge_descriptor;

  using Concurrency_tag = typename internal_np::Lookup_named_param_def <
                                                internal_np::concurrency_tag_t,
                                                NamedParameters1,
                                                Sequential_tag
                                              > ::type;
  using GT = typename GetGeomTraits<TriangleMesh1, NamedParameters1>::type;
  GT gt = choose_parameter<GT>(get_parameter(np1, internal_np::geom_traits));

  using Box_1 = Box_intersection_d::Box_with_info_d<double, 3, face_descriptor_1>;
  using Box_2 = Box_intersection_d::Box_with_info_d<double, 3, face_descriptor_2>;

  auto vpm1 = choose_parameter(get_parameter(np1, internal_np::vertex_point),
                               get_const_property_map(vertex_point, tm1));
  auto vpm2 = choose_parameter(get_parameter(np2, internal_np::vertex_point),
                               get_const_property_map(vertex_point, tm2));

  auto triangle = [&](auto fd, const auto &vpm, const auto &tm){
    auto hd = halfedge(fd,tm);
    auto a = get(vpm, source(hd,tm));
    auto b = get(vpm, target(hd,tm));
    auto c = get(vpm, target(next(hd,tm),tm));
    return typename GT::Triangle_3(a, b, c);
  };

  auto bbox = [](auto fd, const auto &vpm, const auto &tm){
    auto hd = halfedge(fd,tm);
    Bbox_3 res = get(vpm, source(hd,tm)).bbox();
    res += get(vpm, target(hd,tm)).bbox();
    res += get(vpm, target(next(hd,tm),tm)).bbox();
    return res;
  };

  std::vector< Box_1 > boxes_1;
  std::vector< Box_2 > boxes_2;
  boxes_1.reserve(faces(tm1).size());
  boxes_2.reserve(faces(tm2).size());
  for(auto fd: faces(tm1))
    boxes_1.emplace_back(bbox(fd, vpm1, tm1), fd);
  for(auto fd: faces(tm2))
    boxes_2.emplace_back(bbox(fd, vpm2, tm2), fd);

  std::vector< Box_1* > boxes_ptr_1;
  std::vector< Box_2* > boxes_ptr_2;
  boxes_ptr_1.reserve(faces(tm1).size());
  boxes_ptr_2.reserve(faces(tm2).size());
  for(auto &b: boxes_1)
    boxes_ptr_1.emplace_back(std::addressof(b));
  for(auto &b: boxes_2)
    boxes_ptr_2.emplace_back(std::addressof(b));

  CGAL_MUTEX m;
  auto callback=[&](const Box_1 *b1, const Box_2 *b2){
    if(CGAL::do_intersect(triangle(b1->info(), vpm1, tm1), triangle(b2->info(), vpm2, tm2))){
      CGAL_SCOPED_LOCK(m);
      *out ++ = std::make_pair(b1->info(), b2->info());
    }
  };

  const std::ptrdiff_t cutoff = 2000;
  box_intersection_d<Concurrency_tag>(boxes_ptr_1.begin(), boxes_ptr_1.end(), boxes_ptr_2.begin(), boxes_ptr_2.end(), callback, cutoff);
}

template< typename TriangleMesh,
          typename OutputIterator,
          typename NamedParameters = parameters::Default_named_parameters>
void AABB_two_tree_self_intersections(const TriangleMesh &tm, OutputIterator out, const NamedParameters& np = parameters::default_values())
{
  using CGAL::parameters::choose_parameter;
  using CGAL::parameters::get_parameter;

  using face_descriptor = typename boost::graph_traits<TriangleMesh>::face_descriptor;
  using halfedge_descriptor = typename boost::graph_traits<TriangleMesh>::halfedge_descriptor;

  using Concurrency_tag = typename internal_np::Lookup_named_param_def <
                                                internal_np::concurrency_tag_t,
                                                NamedParameters,
                                                Sequential_tag
                                              > ::type;
  using GT = typename GetGeomTraits<TriangleMesh, NamedParameters>::type;
  GT gt = choose_parameter<GT>(get_parameter(np, internal_np::geom_traits));

  using Primitive = AABB_face_graph_triangle_primitive<TriangleMesh>;
  using Face_bbox_tag = typename CGAL::dynamic_face_property_t<Bbox_3>  ;
  using Bbox_pmap = typename boost::property_map<TriangleMesh, Face_bbox_tag>::const_type;
  using Traits = AABB_traits_3<GT, Primitive, Bbox_pmap>;
  using Tree = AABB_tree<Traits>;

  auto vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                              get_const_property_map(vertex_point, tm));

  auto bbox = [&](face_descriptor fd){
    halfedge_descriptor hd = halfedge(fd,tm);
    Bbox_3 res = get(vpm, source(hd,tm)).bbox();
    res += get(vpm, target(hd,tm)).bbox();
    res += get(vpm, target(next(hd,tm),tm)).bbox();
    return res;
  };

  Bbox_pmap bb = get(Face_bbox_tag(), tm);
  for(face_descriptor fd : faces(tm))
    put(bb, fd, bbox(fd));

  Traits traits(bb);
  Tree tree(traits);
  tree.insert(faces(tm).first, faces(tm).second, tm);
  tree.template build<Concurrency_tag>();

#ifdef CGAL_LINKED_WITH_TBB

#endif
  using InternOutputIterator= std::back_insert_iterator<std::vector<std::pair<face_descriptor, face_descriptor>>>;
  std::vector<std::pair<face_descriptor, face_descriptor>> inter;
  CGAL::internal::AABB_tree::two_tree_listing_intersecting_primitives<Concurrency_tag>(tree, tree, std::back_inserter(inter));
  for(const auto& [f_1, f_2]: inter)
    if(f_1 < f_2)
      if(Polygon_mesh_processing::internal::do_faces_intersect<GT>(f_1, f_2, tm, tm.points(), gt.construct_segment_3_object(), gt.construct_triangle_3_object(), gt.do_intersect_3_object()))
        *out ++ = std::make_pair(f_1, f_2);
}

template< typename TriangleMesh1,
          typename TriangleMesh2,
          typename OutputIterator,
          typename NamedParameters1 = parameters::Default_named_parameters,
          typename NamedParameters2 = parameters::Default_named_parameters>
void AABB_two_tree_meshes_intersections(const TriangleMesh1 &tm1,
                                        const TriangleMesh2 &tm2,
                                        OutputIterator out,
                                        const NamedParameters1& np1 = parameters::default_values(),
                                        const NamedParameters2& np2 = parameters::default_values())
{
  using CGAL::parameters::choose_parameter;
  using CGAL::parameters::get_parameter;

  using face_descriptor_1 = typename boost::graph_traits<TriangleMesh1>::face_descriptor;
  using halfedge_descriptor_1 = typename boost::graph_traits<TriangleMesh1>::halfedge_descriptor;

  using face_descriptor_2 = typename boost::graph_traits<TriangleMesh2>::face_descriptor;
  using halfedge_descriptor_2 = typename boost::graph_traits<TriangleMesh2>::halfedge_descriptor;

  using Concurrency_tag = typename internal_np::Lookup_named_param_def <
                                                internal_np::concurrency_tag_t,
                                                NamedParameters1,
                                                Sequential_tag
                                              > ::type;
  using GT = typename GetGeomTraits<TriangleMesh1, NamedParameters1>::type;
  GT gt = choose_parameter<GT>(get_parameter(np1, internal_np::geom_traits));

  using Face_bbox_tag = typename CGAL::dynamic_face_property_t<Bbox_3>;
  using Primitive_1 = AABB_face_graph_triangle_primitive<TriangleMesh1>;
  using Bbox_pmap_1 = typename boost::property_map<TriangleMesh1, Face_bbox_tag>::const_type;
  using Traits_1 = AABB_traits_3<GT, Primitive_1, Bbox_pmap_1>;
  using Tree_1 = AABB_tree<Traits_1>;

  using Primitive_2 = AABB_face_graph_triangle_primitive<TriangleMesh2>;
  using Bbox_pmap_2 = typename boost::property_map<TriangleMesh2, Face_bbox_tag>::const_type;
  using Traits_2 = AABB_traits_3<GT, Primitive_2, Bbox_pmap_2>;
  using Tree_2 = AABB_tree<Traits_2>;

  auto vpm1 = choose_parameter(get_parameter(np1, internal_np::vertex_point),
                               get_const_property_map(vertex_point, tm1));
  auto vpm2 = choose_parameter(get_parameter(np2, internal_np::vertex_point),
                               get_const_property_map(vertex_point, tm2));

  auto bbox = [](auto fd, const auto &vpm, const auto &tm){
    auto hd = halfedge(fd,tm);
    Bbox_3 res = get(vpm, source(hd,tm)).bbox();
    res += get(vpm, target(hd,tm)).bbox();
    res += get(vpm, target(next(hd,tm),tm)).bbox();
    return res;
  };

  Bbox_pmap_1 bb1 = get(Face_bbox_tag(), tm1);
  Bbox_pmap_2 bb2 = get(Face_bbox_tag(), tm2);
#ifdef CGAL_LINKED_WITH_TBB
  if constexpr(std::is_same_v<Concurrency_tag, Parallel_tag>)
  {
    oneapi::tbb::parallel_for(
      oneapi::tbb::blocked_range<size_t>(0, faces(tm1).size()),
        [&](const oneapi::tbb::blocked_range<size_t>& r) {
          for (size_t i = r.begin(); i < r.end(); ++i) {
            face_descriptor_1 fd = *(faces(tm1).begin() + i);
            put(bb1, fd, bbox(fd, vpm1, tm1));
          }
        }
    );

    oneapi::tbb::parallel_for(
      oneapi::tbb::blocked_range<size_t>(0, faces(tm2).size()),
        [&](const oneapi::tbb::blocked_range<size_t>& r) {
          for (size_t i = r.begin(); i < r.end(); ++i) {
            face_descriptor_2 fd = *(faces(tm2).begin() + i);
            put(bb2, fd, bbox(fd, vpm2, tm2));
          }
        }
    );
  }
  else
#endif
  {
    for(face_descriptor_1 fd : faces(tm1))
      put(bb1, fd, bbox(fd, vpm1, tm1));
    for(face_descriptor_2 fd : faces(tm2))
      put(bb2, fd, bbox(fd, vpm2, tm2));
  }

  Traits_1 traits1(bb1);
  Tree_1 tree1(traits1);
  tree1.insert(faces(tm1).first, faces(tm1).second, tm1);

  Traits_2 traits2(bb2);
  Tree_2 tree2(traits2);
  tree2.insert(faces(tm2).first, faces(tm2).second, tm2);

#ifdef CGAL_LINKED_WITH_TBB
  if constexpr(std::is_same_v<Concurrency_tag, Parallel_tag>)
  {
    oneapi::tbb::task_group tg;
    tg.run([&]{ tree1.template build<Concurrency_tag>(); });
    tree2.template build<Concurrency_tag>();
    tg.wait();

    using InternOutputIterator= std::back_insert_iterator<std::vector<std::pair<face_descriptor_1, face_descriptor_2>>>;
    tbb::concurrent_vector<std::pair<face_descriptor_1, face_descriptor_2>> inter;
    CGAL::internal::AABB_tree::two_tree_listing_intersecting_primitives<Concurrency_tag>(tree1, tree2, std::back_inserter(inter));
    for(const auto& p: inter)
      *out++ = p;
  }
  else
#endif
  {
    tree1.template build();
    tree2.template build();
    using InternOutputIterator= std::back_insert_iterator<std::vector<std::pair<face_descriptor_1, face_descriptor_2>>>;
    std::vector<std::pair<face_descriptor_1, face_descriptor_2>> inter;
    CGAL::internal::AABB_tree::two_tree_listing_intersecting_primitives(tree1, tree2, std::back_inserter(inter));
    for(const auto& [f_1, f_2]: inter)
      *out ++ = std::make_pair(f_1, f_2);
  }
}

} // end of namespace experimental

template< typename TriangleMesh,
          typename OutputIterator,
          typename NamedParameters = parameters::Default_named_parameters>
void AABB_self_intersections(const TriangleMesh &tm, OutputIterator out, const NamedParameters& np = parameters::default_values())
{
  using CGAL::parameters::choose_parameter;
  using CGAL::parameters::get_parameter;

  using face_descriptor = typename boost::graph_traits<TriangleMesh>::face_descriptor;
  using halfedge_descriptor = typename boost::graph_traits<TriangleMesh>::halfedge_descriptor;

  using Concurrency_tag = typename internal_np::Lookup_named_param_def <
                                                internal_np::concurrency_tag_t,
                                                NamedParameters,
                                                Sequential_tag
                                              > ::type;
  using GT = typename GetGeomTraits<TriangleMesh, NamedParameters>::type;
  GT gt = choose_parameter<GT>(get_parameter(np, internal_np::geom_traits));

  using Primitive = AABB_face_graph_triangle_primitive<TriangleMesh>;
  using Face_bbox_tag = typename CGAL::dynamic_face_property_t<Bbox_3>  ;
  using Bbox_pmap = typename boost::property_map<TriangleMesh, Face_bbox_tag>::const_type;
  using Traits = AABB_traits_3<GT, Primitive, Bbox_pmap>;
  using Tree = AABB_tree<Traits>;

  auto vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                              get_const_property_map(vertex_point, tm));

  auto bbox = [&](face_descriptor fd){
    halfedge_descriptor hd = halfedge(fd,tm);
    Bbox_3 res = get(vpm, source(hd,tm)).bbox();
    res += get(vpm, target(hd,tm)).bbox();
    res += get(vpm, target(next(hd,tm),tm)).bbox();
    return res;
  };

  Bbox_pmap bb = get(Face_bbox_tag(), tm);
  for(face_descriptor fd : faces(tm))
    put(bb, fd, bbox(fd));

  Traits traits(bb);
  Tree tree(traits);
  tree.insert(faces(tm).first, faces(tm).second, tm);
  tree.template build<Concurrency_tag>();

#ifdef CGAL_LINKED_WITH_TBB
  if constexpr(std::is_same_v<Concurrency_tag, Parallel_tag>)
  {
    CGAL_MUTEX mutex;
    tbb::concurrent_vector<std::pair<face_descriptor, face_descriptor>> inter;
    std::vector<face_descriptor> face_vec(faces(tm).begin(), faces(tm).end());
    oneapi::tbb::parallel_for(
      oneapi::tbb::blocked_range<std::size_t>(0, face_vec.size()),
      [&](const oneapi::tbb::blocked_range<std::size_t>& r) {
        for (std::size_t i = r.begin(); i != r.end(); ++i) {
          face_descriptor f_1 = face_vec[i];
          std::vector<face_descriptor> inter;
          tree.all_intersected_primitives(get(bb, f_1), std::back_inserter(inter));
          for(face_descriptor f_2: inter)
            if(f_1 < f_2)
              if(Polygon_mesh_processing::internal::do_faces_intersect<GT>(f_1, f_2, tm, tm.points(), gt.construct_segment_3_object(), gt.construct_triangle_3_object(), gt.do_intersect_3_object())){
                CGAL_SCOPED_LOCK(mutex);
                *out ++ = std::make_pair(f_1, f_2);
              }
        }
      }
    );
  }
  else
#endif
  {
    for (face_descriptor f_1: faces(tm)){
      std::vector<face_descriptor> inter;
      tree.all_intersected_primitives(get(bb, f_1), std::back_inserter(inter));
      for(face_descriptor f_2: inter)
        if(f_1 < f_2)
          if(Polygon_mesh_processing::internal::do_faces_intersect<GT>(f_1, f_2, tm, tm.points(), gt.construct_segment_3_object(), gt.construct_triangle_3_object(), gt.do_intersect_3_object()))
            *out ++ = std::make_pair(f_1, f_2);
    }
  }
}

template< typename TriangleMesh1,
          typename TriangleMesh2,
          typename OutputIterator,
          typename NamedParameters1 = parameters::Default_named_parameters,
          typename NamedParameters2 = parameters::Default_named_parameters>
void AABB_meshes_intersections(const TriangleMesh1 &tm1,
                               const TriangleMesh2 &tm2,
                               OutputIterator out,
                               const NamedParameters1& np1 = parameters::default_values(),
                               const NamedParameters2& np2 = parameters::default_values())
{
  using CGAL::parameters::choose_parameter;
  using CGAL::parameters::get_parameter;

  using face_descriptor_1 = typename boost::graph_traits<TriangleMesh1>::face_descriptor;
  using halfedge_descriptor_1 = typename boost::graph_traits<TriangleMesh1>::halfedge_descriptor;

  using face_descriptor_2 = typename boost::graph_traits<TriangleMesh1>::face_descriptor;
  using halfedge_descriptor_2 = typename boost::graph_traits<TriangleMesh1>::halfedge_descriptor;

  using Concurrency_tag = typename internal_np::Lookup_named_param_def <
                                                internal_np::concurrency_tag_t,
                                                NamedParameters1,
                                                Sequential_tag
                                              > ::type;
  using GT = typename GetGeomTraits<TriangleMesh1, NamedParameters1>::type;
  GT gt = choose_parameter<GT>(get_parameter(np1, internal_np::geom_traits));

  using Primitive = AABB_face_graph_triangle_primitive<TriangleMesh2>;
  using Face_bbox_tag = typename CGAL::dynamic_face_property_t<Bbox_3>  ;
  using Bbox_pmap = typename boost::property_map<TriangleMesh2, Face_bbox_tag>::const_type;
  using Traits = AABB_traits_3<GT, Primitive, Bbox_pmap>;
  using Tree = AABB_tree<Traits>;

  auto vpm1 = choose_parameter(get_parameter(np1, internal_np::vertex_point),
                               get_const_property_map(vertex_point, tm1));
  auto vpm2 = choose_parameter(get_parameter(np2, internal_np::vertex_point),
                               get_const_property_map(vertex_point, tm2));

  auto triangle = [&](face_descriptor_1 fd){
    halfedge_descriptor_1 hd = halfedge(fd,tm1);
    auto a = get(vpm1, source(hd,tm1));
    auto b = get(vpm1, target(hd,tm1));
    auto c = get(vpm1, target(next(hd,tm1),tm1));
    return typename GT::Triangle_3(a, b, c);
  };

  auto bbox = [&](face_descriptor_2 fd){
    halfedge_descriptor_2 hd = halfedge(fd,tm2);
    Bbox_3 res = get(vpm2, source(hd,tm2)).bbox();
    res += get(vpm2, target(hd,tm2)).bbox();
    res += get(vpm2, target(next(hd,tm2),tm2)).bbox();
    return res;
  };

  Bbox_pmap bb = get(Face_bbox_tag(), tm2);
  for(face_descriptor_2 fd : faces(tm2))
    put(bb, fd, bbox(fd));

  Traits traits(bb);
  Tree tree(traits);
  tree.insert(faces(tm2).first, faces(tm2).second, tm2);
  tree.template build<Concurrency_tag>();

#ifdef CGAL_LINKED_WITH_TBB
  if constexpr(std::is_same_v<Concurrency_tag, Parallel_tag>)
  {
    CGAL_MUTEX mutex;
    tbb::concurrent_vector<std::pair<face_descriptor_1, face_descriptor_2>> inter;
    std::vector<face_descriptor_1> face_vec(faces(tm1).begin(), faces(tm1).end());
    oneapi::tbb::parallel_for(
      oneapi::tbb::blocked_range<size_t>(0, face_vec.size()),
      [&](const oneapi::tbb::blocked_range<size_t>& r) {
        for (size_t i = r.begin(); i != r.end(); ++i) {
          face_descriptor_1 f_1 = face_vec[i];
          std::vector<face_descriptor_2> inter;
          tree.all_intersected_primitives(triangle(f_1), std::back_inserter(inter));
          CGAL_SCOPED_LOCK(mutex);
          for(face_descriptor_2 f_2: inter)
            *out ++ = std::make_pair(f_1, f_2);
        }
      }
    );
  }
  else
#endif
  {
    for (face_descriptor_1 f_1: faces(tm1)){
      std::vector<face_descriptor_2> inter;
      tree.all_intersected_primitives(triangle(f_1), std::back_inserter(inter));
      for(face_descriptor_2 f_2: inter)
        *out ++ = std::make_pair(f_1, f_2);
    }
  }
}

}
} // end namespace CGAL

#endif
