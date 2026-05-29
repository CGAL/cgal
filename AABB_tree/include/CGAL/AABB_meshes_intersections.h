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

#include <CGAL/Polygon_mesh_processing/self_intersections.h>

#ifdef CGAL_LINKED_WITH_TBB
#include <tbb/tbb.h>
#include <CGAL/mutex.h>
#endif

namespace CGAL {

namespace Polygon_mesh_processing{

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
      oneapi::tbb::blocked_range<size_t>(0, face_vec.size()),
      [&](const oneapi::tbb::blocked_range<size_t>& r) {
        for (size_t i = r.begin(); i != r.end(); ++i) {
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
