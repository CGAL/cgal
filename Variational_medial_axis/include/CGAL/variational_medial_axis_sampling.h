// Copyright (c) 2025 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Qijia Huang

#ifndef CGAL_VARIATIONAL_MEDIAL_AXIS_SAMPLING_H
#define CGAL_VARIATIONAL_MEDIAL_AXIS_SAMPLING_H

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits_3.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Side_of_triangle_mesh.h>

namespace CGAL
{

/**
 * \ingroup PkgVMASRef
 * computes a static skeleton based on variational medial axis sampling method.
 *
 * @tparam TriangleMesh a model of `HalfedgeListGraph`, `FaceListGraph`, and `MutableFaceGraph`
 * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * @param tmesh input triangulated surface mesh
 * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{concurrency_tag}
 *     \cgalParamDescription{a tag indicating if the task should be done using one or several threads.}
 *     \cgalParamType{Either `CGAL::Sequential_tag`, or `CGAL::Parallel_tag`, or `CGAL::Parallel_if_available_tag`}
 *     \cgalParamDefault{`CGAL::Sequential_tag`}
 *   \cgalParamNEnd
 *   \cgalParamNBegin{geom_traits}
 *     \cgalParamDescription{an instance of a geometric traits class}
 *     \cgalParamType{a class model of `PMPSelfIntersectionTraits`}
 *     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
 *     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
 *   \cgalParamNEnd
 *   \cgalParamNBegin{vertex_point_map}
 *     \cgalParamDescription{a property map associating points to the vertices of `tm`}
 *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
 *                    as key type and `%Point_3` as value type}
 *     \cgalParamDefault{`boost::get(CGAL::vertex_point, tm)`}
 *     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
 *                     must be available in `TriangleMesh`.}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 */
template <class TriangleMesh, class NamedParameters = parameters::Default_named_parameters>
void
variational_medial_axis_sampling(const TriangleMesh& tmesh,
                                 const NamedParameters& np= parameters::default_values())
{
  namespace PMP = CGAL::Polygon_mesh_processing;
  using parameters::choose_parameter;
  using parameters::get_parameter;

  using GT = typename GetGeomTraits<TriangleMesh, NamedParameters>::type;
  GT traits = choose_parameter<GT>(get_parameter(np, internal_np::geom_traits));

  typedef typename internal_np::Lookup_named_param_def <
    internal_np::concurrency_tag_t,
    NamedParameters,
    Sequential_tag
  > ::type Concurrency_tag;

  constexpr bool parallel_execution = std::is_same_v<Parallel_tag, Concurrency_tag>;

  using VPM = typename CGAL::GetVertexPointMap<TriangleMesh, NamedParameters>::type;
  VPM vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                             get_const_property_map(CGAL::vertex_point, tmesh));


  using Tree = AABB_tree<AABB_traits_3<GT,AABB_face_graph_triangle_primitive<TriangleMesh,VPM>>>;

  using face_descriptor = typename boost::graph_traits<TriangleMesh>::face_descriptor;
  using halfedge_descriptor = typename boost::graph_traits<TriangleMesh>::halfedge_descriptor;
  using vertex_descriptor = typename boost::graph_traits<TriangleMesh>::vertex_descriptor;

  // build an AABB-tree of faces
  Tree tree(faces(tmesh).begin(), faces(tmesh).end(), tmesh, vpm);
  tree.accelerate_distance_queries();

  // get the closest point to the origin
  typename GT::Point_3 query(0.,0.,0.);
  typename GT::Point_3 cp = tree.closest_point(query);

  std::cout << "closest point is " << cp << "\n";

  face_descriptor closest_face = tree.closest_point_and_primitive(query).second;
  halfedge_descriptor h = halfedge(closest_face, tmesh);

  std::array<vertex_descriptor, 3> fvertices;
  fvertices[0] = target(h, tmesh);
  fvertices[1] = target(next(h, tmesh), tmesh);
  fvertices[2] = target(opposite(h, tmesh), tmesh); // or source(h, tmesh)

  std::cout << "closest triangle is " << get(vpm, fvertices[0]) << " "
                                      << get(vpm, fvertices[1]) << " "
                                      << get(vpm, fvertices[2]) << "\n";


  Side_of_triangle_mesh<TriangleMesh, GT, VPM, Tree> side_of(tree, traits);
  std::cout << "query inside bounded volume? " << (side_of(query)==CGAL::ON_BOUNDED_SIDE) << "\n";

  typename GT::Sphere_3 sphere(query, 10.);
  std::cout << "big sphere intersects triangles? " << tree.do_intersect(sphere) << "\n";

  sphere = typename GT::Sphere_3(cp, 0.01);
  std::cout << "small sphere intersects triangles? " << tree.do_intersect(sphere) << "\n";

  // add a double per face to store the area
  using Face_area_tag = CGAL::dynamic_face_property_t<double>;
  using Face_area_map = typename boost::property_map<TriangleMesh, Face_area_tag>::const_type;

  Face_area_map face_area_map = get(Face_area_tag(), tmesh, 0.); // 0. is the default value


  double total_area = 0.;
  for (face_descriptor f : faces(tmesh))
  {
    put(face_area_map, f, PMP::face_area(f, tmesh, parameters::vertex_point_map(vpm)));
    total_area += get(face_area_map, f);
  }

  std::cout << "Total area: "  << total_area << "\n";
}

} // end of CGAL namespace

#endif // CGAL_VARIATIONAL_MEDIAL_AXIS_SAMPLING_H
