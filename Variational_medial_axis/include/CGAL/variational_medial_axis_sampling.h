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
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Side_of_triangle_mesh.h>
#ifdef CGAL_EIGEN3_ENABLED
#include <CGAL/Eigen_solver_traits.h>
#endif
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
template <class TriangleMesh, class NamedParameters = parameters::Default_named_parameters >
void variational_medial_axis_sampling(const TriangleMesh& tmesh,
                                 const NamedParameters& np= parameters::default_values())
{
  namespace PMP = CGAL::Polygon_mesh_processing;
  using parameters::choose_parameter;
  using parameters::get_parameter;

  using GT                  = typename GetGeomTraits<TriangleMesh, NamedParameters>::type;
  using FT                  = GT::FT;
  using Point_3             = GT::Point_3;
  using Vector_3            = GT::Vector_3;
  using Sphere_3            = GT::Sphere_3;
  using VPM                 = typename CGAL::GetVertexPointMap<TriangleMesh, NamedParameters>::type;
  using Tree                = AABB_tree<AABB_traits_3<GT, AABB_face_graph_triangle_primitive<TriangleMesh, VPM>>>;

  using face_descriptor     = typename boost::graph_traits<TriangleMesh>::face_descriptor;
  using halfedge_descriptor = typename boost::graph_traits<TriangleMesh>::halfedge_descriptor;
  using vertex_descriptor   = typename boost::graph_traits<TriangleMesh>::vertex_descriptor;

  
#ifdef CGAL_EIGEN3_ENABLED
  using Solver_traits = CGAL::Eigen_solver_traits<Eigen::SimplicialLDLT<Eigen_matrix<FT>::EigenType>>;
#endif
  GT traits = choose_parameter<GT>(get_parameter(np, internal_np::geom_traits));

  struct MedialSphere
  {
    using Spehre_ID = std::size_t;
    Sphere_3 sphere;
    Sphere_ID id;
    std::set<Sphere_ID> neighbors;                   // IDs of neighboring spheres
    std::vector<vertex_descriptor> cluster_vertices; // vertices associated with this sphere
    FT error;                                        // error associated with this sphere
    FT cluster_area;
    MedialSphere(const Sphere_3& s, Sphere_ID i)
        : sphere(s)
        , id(i)
        , error(FT(0.0))
        , cluster_area(FT(0.0)) {}
  };

  using Vertex_normal_tag = CGAL::dynamic_vertex_property_t<Vector_3>;
  using Vertex_normal_map = typename boost::property_map<TriangleMesh, Vertex_normal_tag>::const_type;
  using Vertex_area_tag = CGAL::dynamic_vertex_property_t<FT>;
  using Vertex_area_map = typename boost::property_map<TriangleMesh, Vertex_area_tag>::const_type;
  using Vertex_error_tag = CGAL::dynamic_vertex_property_t<FT>;
  using Vertex_error_map = typename boost::property_map<TriangleMesh, Vertex_error_tag>::const_type;
  using Vertex_cluster_sphere = CGAL::dynamic_vertex_property_t<MedialSphere>;
  using Vertex_cluster_sphere_map = typename boost::property_map<TriangleMesh, Vertex_cluster_sphere>::type;

  typedef typename internal_np::Lookup_named_param_def <
    internal_np::concurrency_tag_t,
    NamedParameters,
    Sequential_tag
  > ::type Concurrency_tag;

  constexpr bool parallel_execution = std::is_same_v<Parallel_tag, Concurrency_tag>;

  
  VPM vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                             get_const_property_map(CGAL::vertex_point, tmesh));

  // build an AABB-tree of faces
  Tree tree(faces(tmesh).begin(), faces(tmesh).end(), tmesh, vpm);
  tree.accelerate_distance_queries();

  
  Vertex_normal_map vertex_normal_map = get(Vertex_normal_tag(), tmesh, Vector_3(0., 0., 0.));
  Vertex_area_map vertex_area_map = get(Vertex_area_tag(), tmesh, 0.);
  // compute vertex normal
  PMP::compute_vertex_normals(tmesh, vertex_normal_map, parameters::geom_traits(traits).vertex_point_map(vpm));

  // compute vertex areas
  for(face_descriptor f : faces(tmesh)) {
    double area = PMP::face_area(f, tmesh, parameters::vertex_point_map(vpm));
    for(vertex_descriptor v : vertices_around_face(halfedge(f, tmesh), tmesh)) {
      put(vertex_area_map, v, get(vertex_area_map, v) + area / 3.0);
    }
  }

  
  std::vector<MedialSphere> medial_spheres;
  int desired_number_of_spheres = 100; // TO DO: pass desired number of sphere by parameter
  medial_spheres.reserve(desired_number_of_spheres); 
  medial_spheres.emplace_back(Sphere_3(Point_3(0., 0., 0.), FT(1.0)), 0);

  while(medial_spheres.size() < desired_number_of_spheres) {
    //Clean data
    for(auto& sphere : medial_spheres) {
      sphere.neighbors.clear();
      sphere.cluster_vertices.clear();
      sphere.error = FT(0.0);
      sphere.cluster_area = FT(0.0);
    }

    // Compute cluster
    for(vertex_descriptor v : sphere.cluster_vertices) {
      FT area = get(vertex_area_map, v);
      Point_3 p = get(vpm, v);
      FT min_distance = std::numeric_limits<FT>::max();
      Spehre_ID closest_sphere_id;

    }
    //update sphere   
    //split sphere
    // update neighbors
  }


}

template <class GT, class TriangleMesh> 
struct MedialSphere<GT>
{   
    using Sphere_3 = typename GT::Sphere_3;
    using vertex_descriptor = typename boost::graph_traits<TriangleMesh>::vertex_descriptor;
    using Sphere_ID = std::size_t;
    Sphere_3 sphere;
    Sphere_ID id;
    std::set<Sphere_ID> neighbors;         // IDs of neighboring spheres
    std::vector<vertex_descriptor> cluster_vertices; // vertices associated with this sphere
    typename GT::FT error;                           // error associated with this sphere
};


template <class TriangleMesh, class VertexNormalMap,
          class VertexErrorMap, class SphereSet,
          class NamedParameters = parameters::Default_named_parameters>
void compute_cluster_errors(const TriangleMesh& tmesh,
                     VertexNormalMap vertex_normals,
                     VertexErrorMap vertex_error,
                     SphereSet& sphere_set,
                     const NamedParameters& np = parameters::default_values()) {
  for(auto& sphere : sphere_set) {
    
  }
}

} 
// end of CGAL namespace
// get the closest point to the origin
// Point_3 query(0.,0.,0.);
// Point_3 cp = tree.closest_point(query);

// std::cout << "closest point is " << cp << "\n";

// face_descriptor closest_face = tree.closest_point_and_primitive(query).second;
// halfedge_descriptor h = halfedge(closest_face, tmesh);

// std::array<vertex_descriptor, 3> fvertices;
// fvertices[0] = target(h, tmesh);
// fvertices[1] = target(next(h, tmesh), tmesh);
// fvertices[2] = target(opposite(h, tmesh), tmesh); // or source(h, tmesh)

// std::cout << "closest triangle is " << get(vpm, fvertices[0]) << " "
//                                     << get(vpm, fvertices[1]) << " "
//                                     << get(vpm, fvertices[2]) << "\n";

// Side_of_triangle_mesh<TriangleMesh, GT, VPM, Tree> side_of(tree, traits);
// std::cout << "query inside bounded volume? " << (side_of(query)==CGAL::ON_BOUNDED_SIDE) << "\n";

// Sphere_3 sphere(query, 10.);
// std::cout << "big sphere intersects triangles? " << tree.do_intersect(sphere) << "\n";

// sphere = Sphere_3(cp, 0.01);
// std::cout << "small sphere intersects triangles? " << tree.do_intersect(sphere) << "\n";
#endif // CGAL_VARIATIONAL_MEDIAL_AXIS_SAMPLING_H
