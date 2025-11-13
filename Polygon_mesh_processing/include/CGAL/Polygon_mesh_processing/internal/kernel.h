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

#include <CGAL/license/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/clip.h>
#include <CGAL/Convex_hull_3/dual/halfspace_intersection_3.h>
#include <CGAL/Convex_hull_3/dual/halfspace_intersection_with_constructions_3.h>


namespace CGAL {
namespace Polygon_mesh_processing {
namespace experimental {

template <class Kernel>
struct Three_point_cut_plane_traits
{
  using FT = typename Kernel::FT;
  using Plane_3 = std::array<typename Kernel::Point_3, 3>;
  using Point_3 = typename Kernel::Point_3;

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
      typename Kernel::Plane_3 k_plane(plane[0], plane[1], plane[2]);

      return construction(k_plane, p, q);
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


template <class TriangleMesh,
          class NamedParameters = parameters::Default_named_parameters>
TriangleMesh
kernel(const TriangleMesh& pm,
       const NamedParameters& np = parameters::default_values())
{
  // TODO: bench if with EPECK we can directly use Kernel::Plane_3
  // TODO: TriangleMesh as output is not correct since it is actually a PolygonMesh
  using parameters::choose_parameter;
  using parameters::get_parameter;

  using GT = typename GetGeomTraits<TriangleMesh, NamedParameters>::type;
  auto vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                              get_const_property_map(vertex_point, pm));
  // GT gt = choose_parameter<GT>(get_parameter(np, internal_np::geom_traits));
  using Point_3 = typename GT::Point_3;
  using parameters::choose_parameter;
  using parameters::get_parameter;

  //TODO: what do we do with a mesh that is not closed?
  //TODO: what do we do if the input is not a triangle mesh?

  if (vertices(pm).size() - edges(pm).size() + faces(pm).size() != 2)
    return TriangleMesh();



  CGAL::Bbox_3 bb3 = bbox(pm, np);
  TriangleMesh kernel;
  CGAL::make_hexahedron(Point_3(bb3.xmax(),bb3.ymin(),bb3.zmin()), Point_3(bb3.xmax(),bb3.ymax(),bb3.zmin()), Point_3(bb3.xmin(),bb3.ymax(),bb3.zmin()), Point_3(bb3.xmin(),bb3.ymin(),bb3.zmin()),
                        Point_3(bb3.xmin(),bb3.ymin(),bb3.zmax()), Point_3(bb3.xmax(),bb3.ymin(),bb3.zmax()), Point_3(bb3.xmax(),bb3.ymax(),bb3.zmax()), Point_3(bb3.xmin(),bb3.ymax(),bb3.zmax()),
                        kernel);

  Three_point_cut_plane_traits<GT> kgt;
  for (auto f : faces(pm))
  {
    auto h = halfedge(f, pm);
    auto plane = make_array( get(vpm,source(h, pm)),
                             get(vpm,target(h, pm)),
                             get(vpm,target(next(h, pm), pm)) );

#ifdef CGAL_USE_OPTI_WITH_BBOX
    auto pred = kgt.oriented_side_3_object();
    auto gbox = bbox(kernel);
    std::array<Point_3, 8> corners = CGAL::make_array(Point_3(bb3.xmax(),bb3.ymin(),bb3.zmin()),
                                                      Point_3(bb3.xmax(),bb3.ymax(),bb3.zmin()),
                                                      Point_3(bb3.xmin(),bb3.ymax(),bb3.zmin()),
                                                      Point_3(bb3.xmin(),bb3.ymin(),bb3.zmin()),
                                                      Point_3(bb3.xmin(),bb3.ymin(),bb3.zmax()),
                                                      Point_3(bb3.xmax(),bb3.ymin(),bb3.zmax()),
                                                      Point_3(bb3.xmax(),bb3.ymax(),bb3.zmax()),
                                                      Point_3(bb3.xmin(),bb3.ymax(),bb3.zmax()));
    int i=0;
    auto first_ori=pred(plane, corners[i]);
    while(++i!=8 && first_ori==ON_ORIENTED_BOUNDARY)
      first_ori=pred(plane, corners[i]);

    if (i==8) continue;
    bool all_the_same=true;
    for (;i<8;++i)
    {
      auto other_ori=pred(plane, corners[i]);
      if (other_ori!=ON_ORIENTED_BOUNDARY && other_ori!=first_ori)
      {
        all_the_same=false;
        break;
      }
    }

    if (all_the_same)
    {
      if (first_ori==ON_NEGATIVE_SIDE) continue;
      else
      {
        return TriangleMesh();
      }
    }
#endif

    clip(kernel, plane, CGAL::parameters::clip_volume(true).geom_traits(kgt).do_not_triangulate_faces(true).used_for_kernel(true));
    if (is_empty(kernel)) break;
  }

  return kernel;
};

template <class TriangleMesh,
          class NamedParameters = parameters::Default_named_parameters>
TriangleMesh
kernel_using_chull(const TriangleMesh& pm,
                   const NamedParameters& np = parameters::default_values())
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  using GT = typename GetGeomTraits<TriangleMesh, NamedParameters>::type;
  auto vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                              get_const_property_map(vertex_point, pm));
  // GT gt = choose_parameter<GT>(get_parameter(np, internal_np::geom_traits));
  using parameters::choose_parameter;
  using parameters::get_parameter;

  std::vector<typename GT::Plane_3> planes;
  planes.reserve(faces(pm).size());


  for (auto f : faces(pm))
  {
    auto h = halfedge(f, pm);
    planes.emplace_back( get(vpm,source(h, pm)),
                         get(vpm,target(h, pm)),
                         get(vpm,target(next(h, pm), pm)) );
  }

  TriangleMesh kernel;

  // if no point inside the intersection is provided, one
  // will be automatically found using linear programming
  CGAL::halfspace_intersection_3(planes.begin(),
                                 planes.end(),
                                 kernel );

  return kernel;
}

template <class TriangleMesh,
          class NamedParameters = parameters::Default_named_parameters>
TriangleMesh
kernel_using_chull_and_constructions(const TriangleMesh& pm,
                   const NamedParameters& np = parameters::default_values())
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  using GT = typename GetGeomTraits<TriangleMesh, NamedParameters>::type;
  auto vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                              get_const_property_map(vertex_point, pm));
  // GT gt = choose_parameter<GT>(get_parameter(np, internal_np::geom_traits));
  using parameters::choose_parameter;
  using parameters::get_parameter;

  std::vector<typename GT::Plane_3> planes;
  planes.reserve(faces(pm).size());


  for (auto f : faces(pm))
  {
    auto h = halfedge(f, pm);
    planes.emplace_back( get(vpm,source(h, pm)),
                         get(vpm,target(h, pm)),
                         get(vpm,target(next(h, pm), pm)) );
  }

  TriangleMesh kernel;

  // if no point inside the intersection is provided, one
  // will be automatically found using linear programming
  CGAL::halfspace_intersection_with_constructions_3(planes.begin(),
                                                    planes.end(),
                                                    kernel );

  return kernel;
}





} } } // end of CGAL::Polygon_mesh_processing::experimental

#endif // CGAL_POLYGON_MESH_PROCESSING_INTERNAL_KERNEL_H