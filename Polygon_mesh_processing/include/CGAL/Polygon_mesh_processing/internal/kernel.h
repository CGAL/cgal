// Copyright (c) 2016-2025 GeometryFactory (France).
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

#ifndef CGAL_POLYGON_MESH_PROCESSING_INTERNAL_KERNEL_H
#define CGAL_POLYGON_MESH_PROCESSING_INTERNAL_KERNEL_H

#include <algorithm>
#include <random>

#include <CGAL/license/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/clip.h>
#include <CGAL/Polygon_mesh_processing/clip_convex.h>
#include <CGAL/Convex_hull_3/dual/halfspace_intersection_3.h>
#include <CGAL/Convex_hull_3/dual/halfspace_intersection_with_constructions_3.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Exact_integer.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Cartesian_converter.h>

namespace CGAL {
namespace Polygon_mesh_processing {

template <class Kernel>
struct Three_point_cut_plane_traits
{
  using FT = typename Kernel::FT;
  using Point_3 = typename Kernel::Point_3;

  struct Plane_3: public std::array<typename Kernel::Point_3, 3>{
    using Base = std::array<typename Kernel::Point_3, 3>;
    using Explicit_plane = typename Kernel::Plane_3;

    Plane_3(const Point_3 &a, const Point_3 &b, const Point_3 &c): Base({a, b, c}){}
    Plane_3(const std::array<Point_3, 3> &arr): Base(arr){}

    // Warning: it is slow (Planes are constructed each time)
    bool operator<(const Plane_3 &b) const{
      Explicit_plane pa = explicit_plane();
      Explicit_plane pb = b.explicit_plane();
      Comparison_result res = compare(pa.a(), pb.a());
      if(res == EQUAL)
        res = compare(pa.b(), pb.b());
      if(res == EQUAL)
        res = compare(pa.c(), pb.c());
      if(res == EQUAL)
        res = compare(pa.d(), pb.d());
      return res == SMALLER;
    };

    bool operator==(const Plane_3 &b) const{
      Explicit_plane pa = explicit_plane();
      Explicit_plane pb = b.explicit_plane();
      return pa==pb;
    }

    Explicit_plane explicit_plane() const{
      return  Explicit_plane((*this)[0], (*this)[1], (*this)[2]);
    }

  };
  using Vector_3 = typename Kernel::Vector_3;

  struct Does_not_support_CDT2{};

  struct Oriented_side_3
  {
    Oriented_side operator()(const Plane_3& plane, const Point_3& p)  const
    {
      return orientation(plane[0], plane[1], plane[2], p);
    }
  };

  struct Construct_plane_line_intersection_point_3
  {
    Point_3 operator()(const Plane_3& plane, const Point_3& p, const Point_3& q)
    {
      typename Kernel::Construct_plane_line_intersection_point_3 construction;
      return construction(plane[0], plane[1], plane[2], p, q);
    }
  };

  struct Construct_orthogonal_vector_3{
    Vector_3 operator()(const Plane_3& plane)
    {
      return typename Kernel::Plane_3(plane[0], plane[1], plane[2]).orthogonal_vector();
    }
  };

  struct Compute_squared_distance_3
  {
    using Compute_scalar_product_3 = typename Kernel::Compute_scalar_product_3;
    FT operator()(const Plane_3& plane, const Point_3& p)
    {
      typename Kernel::Plane_3 pl(plane[0], plane[1], plane[2]);
      return Compute_scalar_product_3()(Vector_3(ORIGIN, p), pl.orthogonal_vector())+pl.d();
    }
  };

  Oriented_side_3 oriented_side_3_object() const
  {
    return Oriented_side_3();
  }

  Construct_plane_line_intersection_point_3 construct_plane_line_intersection_point_3_object() const
  {
    return Construct_plane_line_intersection_point_3();
  }

  Construct_orthogonal_vector_3 construct_orthogonal_vector_3_object() const
  {
    return Construct_orthogonal_vector_3();
  }

  Compute_squared_distance_3 compute_squared_distance_3_object() const { return Compute_squared_distance_3(); }

#ifndef CGAL_PLANE_CLIP_DO_NOT_USE_BOX_INTERSECTION_D
// for does self-intersect
  using Segment_3 = typename Kernel::Segment_3;
  using Triangle_3 = typename Kernel::Triangle_3;
  using Construct_segment_3 = typename Kernel::Construct_segment_3;
  using Construct_triangle_3 =typename  Kernel::Construct_triangle_3;
  using Do_intersect_3 = typename Kernel::Do_intersect_3;
  Construct_segment_3 construct_segment_3_object() const { return Construct_segment_3(); }
  Construct_triangle_3 construct_triangle_3_object() const { return Construct_triangle_3(); }
  Do_intersect_3 do_intersect_3_object() const { return Do_intersect_3(); }
#endif
};

namespace internal{
template <class PolygonMesh,
          class FaceRange,
          class NamedParameters = parameters::Default_named_parameters>
PolygonMesh
kernel(const PolygonMesh& pm,
       const FaceRange& faces,
       const NamedParameters& np = parameters::default_values())
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  // graph typedefs
  using BGT = boost::graph_traits<PolygonMesh>;
  using face_descriptor = typename BGT::face_descriptor;
  // using edge_descriptor = typename BGT::edge_descriptor;
  // using halfedge_descriptor = typename BGT::halfedge_descriptor;
  using vertex_descriptor = typename BGT::vertex_descriptor;

  using GT = typename GetGeomTraits<PolygonMesh, NamedParameters>::type;
  using EK = Exact_predicates_exact_constructions_kernel;
  using K2EK = Cartesian_converter<GT, EK>;
  using EK2K = Cartesian_converter<EK, GT>;
  K2EK to_exact;
  EK2K from_exact;
  auto vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                              get_const_property_map(vertex_point, pm));

  using Point_3 = typename GT::Point_3;
  using EPoint_3 = typename EK::Point_3;
  using Plane_3 = typename Three_point_cut_plane_traits<EK>::Plane_3;

  using KernelPointMap = typename boost::property_map<PolygonMesh, dynamic_vertex_property_t<EPoint_3> >::type;

  bool bbox_filtering = choose_parameter(get_parameter(np, internal_np::use_bounding_box_filtering), true);
  bool shuffle_planes = choose_parameter(get_parameter(np, internal_np::shuffle_planes), true);

  // To speedup on stupid benchmarks
  // if (vertices(pm).size() - edges(pm).size() + faces(pm).size() != 2)
  //   return PolygonMesh();

  // Build the starting cube
  CGAL::Bbox_3 bb3 = bbox(pm, np);
  PolygonMesh kernel;
  CGAL::make_hexahedron(Point_3(bb3.xmax(),bb3.ymin(),bb3.zmin()), Point_3(bb3.xmax(),bb3.ymax(),bb3.zmin()), Point_3(bb3.xmin(),bb3.ymax(),bb3.zmin()), Point_3(bb3.xmin(),bb3.ymin(),bb3.zmin()),
                        Point_3(bb3.xmin(),bb3.ymin(),bb3.zmax()), Point_3(bb3.xmax(),bb3.ymin(),bb3.zmax()), Point_3(bb3.xmax(),bb3.ymax(),bb3.zmax()), Point_3(bb3.xmin(),bb3.ymax(),bb3.zmax()),
                        kernel);
  auto base_vpm = get_property_map(vertex_point, kernel);
  KernelPointMap kvpm = get(CGAL::dynamic_vertex_property_t<EPoint_3>(), kernel);
  for(vertex_descriptor v: vertices(kernel))
    put(kvpm, v, to_exact(get(base_vpm, v)));
  vertex_descriptor start_vertex = *vertices(pm).begin();

  std::array<vertex_descriptor, 6> bbox_vertices;
  std::array<EPoint_3, 8> corners;
  if(bbox_filtering){
    // We store the vertices that realized the bbox
    struct BBoxEntry {
      std::size_t index;
      std::function<double(const EPoint_3&)> bound;
      std::function<double(const Bbox_3&)> value;
    };
    std::array<BBoxEntry,6> entries {{
        {0, [](const EPoint_3& p){ return to_interval(p.x()).first;  }, [](const Bbox_3& b){ return b.xmin(); }},
        {1, [](const EPoint_3& p){ return to_interval(p.x()).second; }, [](const Bbox_3& b){ return b.xmax(); }},
        {2, [](const EPoint_3& p){ return to_interval(p.y()).first;  }, [](const Bbox_3& b){ return b.ymin(); }},
        {3, [](const EPoint_3& p){ return to_interval(p.y()).second; }, [](const Bbox_3& b){ return b.ymax(); }},
        {4, [](const EPoint_3& p){ return to_interval(p.z()).first;  }, [](const Bbox_3& b){ return b.zmin(); }},
        {5, [](const EPoint_3& p){ return to_interval(p.z()).second; }, [](const Bbox_3& b){ return b.zmax(); }}
    }};

    for (const auto& e : entries){
      for (vertex_descriptor v : vertices(kernel)){
        std::size_t i = e.index;
        double bound = e.bound(get(kvpm, v));
        if (bound == e.value(bb3)){
          bbox_vertices[i] = v;
          break;
        }
      }
    }
    corners = CGAL::make_array(EPoint_3(bb3.xmin(),bb3.ymin(),bb3.zmin()),
                               EPoint_3(bb3.xmin(),bb3.ymin(),bb3.zmax()),
                               EPoint_3(bb3.xmin(),bb3.ymax(),bb3.zmin()),
                               EPoint_3(bb3.xmin(),bb3.ymax(),bb3.zmax()),
                               EPoint_3(bb3.xmax(),bb3.ymin(),bb3.zmin()),
                               EPoint_3(bb3.xmax(),bb3.ymin(),bb3.zmax()),
                               EPoint_3(bb3.xmax(),bb3.ymax(),bb3.zmin()),
                               EPoint_3(bb3.xmax(),bb3.ymax(),bb3.zmax()));
  }

  Three_point_cut_plane_traits<EK> kgt;

  std::vector<face_descriptor> planes(faces.begin(), faces.end());
  if(shuffle_planes)
    std::shuffle(planes.begin(), planes.end(), std::default_random_engine());

  for(auto f: planes)
  {
    auto h = halfedge(f, pm);
    Plane_3 plane(to_exact(get(vpm,source(h, pm))),
                  to_exact(get(vpm,target(h, pm))),
                  to_exact(get(vpm,target(next(h, pm), pm))));

    if(bbox_filtering){
      // By looking the sign of the plane value, we can check only two corners
      auto pred = kgt.oriented_side_3_object();
      auto eplane = plane.explicit_plane();
      std::size_t index = (is_positive(eplane.a())?4:0) + (is_positive(eplane.b())?2:0) + (is_positive(eplane.c())?1:0);
      if(pred(plane, corners[index]) != ON_POSITIVE_SIDE)
        continue;
      if(pred(plane, corners[7-index]) == ON_POSITIVE_SIDE)
        return PolygonMesh(); // empty

      start_vertex = clip_convex(kernel, plane, CGAL::parameters::clip_volume(true).geom_traits(kgt).do_not_triangulate_faces(true).vertex_point_map(kvpm).
                                                                  bounding_box(&bbox_vertices).starting_vertex_descriptor(start_vertex));
      if (is_empty(kernel)) return kernel;

      CGAL_assertion_code(for(std::size_t i=0; i!=6; ++i))
        CGAL_assertion(kernel.is_valid(bbox_vertices[i]));

      // update bbox TODO can be smarter
      bb3 = get(kvpm, bbox_vertices[0]).bbox()+get(kvpm, bbox_vertices[1]).bbox()+get(kvpm, bbox_vertices[2]).bbox()+
            get(kvpm, bbox_vertices[3]).bbox()+get(kvpm, bbox_vertices[4]).bbox()+get(kvpm, bbox_vertices[5]).bbox();
      corners = CGAL::make_array(EPoint_3(bb3.xmin(),bb3.ymin(),bb3.zmin()),
                                 EPoint_3(bb3.xmin(),bb3.ymin(),bb3.zmax()),
                                 EPoint_3(bb3.xmin(),bb3.ymax(),bb3.zmin()),
                                 EPoint_3(bb3.xmin(),bb3.ymax(),bb3.zmax()),
                                 EPoint_3(bb3.xmax(),bb3.ymin(),bb3.zmin()),
                                 EPoint_3(bb3.xmax(),bb3.ymin(),bb3.zmax()),
                                 EPoint_3(bb3.xmax(),bb3.ymax(),bb3.zmin()),
                                 EPoint_3(bb3.xmax(),bb3.ymax(),bb3.zmax()));
    }
    else
    {
      start_vertex = clip_convex(kernel, plane, CGAL::parameters::clip_volume(true).geom_traits(kgt).do_not_triangulate_faces(true).vertex_point_map(kvpm).
                                                                  starting_vertex_descriptor(start_vertex));
      if (is_empty(kernel)) return kernel;
    }
  }

  for(vertex_descriptor v : vertices(kernel))
    put(base_vpm, v, from_exact(get(kvpm, v)));
  return kernel;
};

} // end of namespace internal

template <class PolygonMesh,
          class NamedParameters = parameters::Default_named_parameters>
PolygonMesh
kernel(const PolygonMesh& pm,
       const NamedParameters& np = parameters::default_values())
{
  return internal::kernel(pm, faces(pm), np);
}


} } // end of CGAL::Polygon_mesh_processing

#endif // CGAL_POLYGON_MESH_PROCESSING_INTERNAL_KERNEL_H