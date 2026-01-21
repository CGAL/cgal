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
// Author(s)     : Sebastien Loriot, Léo Valque

#ifndef CGAL_POLYGON_MESH_PROCESSING_KERNEL_H
#define CGAL_POLYGON_MESH_PROCESSING_KERNEL_H

#include <CGAL/license/Polygon_mesh_processing/corefinement.h>

#include <algorithm>
#include <random>

#include <CGAL/Polygon_mesh_processing/clip.h>
#include <CGAL/Polygon_mesh_processing/internal/clip_convex.h>
#include <CGAL/Polygon_mesh_processing/Three_point_cut_plane_traits.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Cartesian_converter.h>

#include <boost/property_map/property_map.hpp>

namespace CGAL {
namespace Polygon_mesh_processing {
namespace internal{
template <class PolygonMesh,
          class FaceRange,
          class NamedParameters = parameters::Default_named_parameters>
PolygonMesh
kernel(const PolygonMesh& pm,
       const FaceRange& face_range,
       const NamedParameters& np = parameters::default_values(),
       bool used_to_find_a_point = false,
       std::optional<typename GetGeomTraits<PolygonMesh, NamedParameters>::type::Point_3> *p = nullptr)
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  // graph typedefs
  using BGT = boost::graph_traits<PolygonMesh>;
  using face_descriptor = typename BGT::face_descriptor;
  // using edge_descriptor = typename BGT::edge_descriptor;
  using halfedge_descriptor = typename BGT::halfedge_descriptor;
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
  using EVector_3 = typename EK::Vector_3;
  using Plane_3 = typename Three_point_cut_plane_traits<EK>::Plane_3;

  using KernelPointMap = typename boost::property_map<PolygonMesh, dynamic_vertex_property_t<EPoint_3> >::type;

  bool bbox_filtering = choose_parameter(get_parameter(np, internal_np::use_bounding_box_filtering), true);
  bool shuffle_planes = choose_parameter(get_parameter(np, internal_np::shuffle_planes), true);
  // bool used_to_find_a_point = choose_parameter(get_parameter(np, internal_np::used_to_find_a_point), false);
  bool check_euler_characteristic = !choose_parameter(get_parameter(np, internal_np::allow_non_manifold_non_watertight_input), false);
  std::size_t seed = choose_parameter(get_parameter(np, internal_np::random_seed), std::random_device()());

  // Immediate exit if the input is well-formed and not of genus zero
  if(check_euler_characteristic && (vertices(pm).size() - edges(pm).size() + faces(pm).size() != 2))
    return PolygonMesh();

  // Build the starting cube
  CGAL::Bbox_3 bb3 = choose_parameter(get_parameter(np, internal_np::starting_cube), bbox(pm));
  PolygonMesh kernel;
  make_hexahedron(bb3, kernel);
  auto base_vpm = get_property_map(vertex_point, kernel);
  KernelPointMap kvpm = get(CGAL::dynamic_vertex_property_t<EPoint_3>(), kernel);
  for(vertex_descriptor v: vertices(kernel))
    put(kvpm, v, to_exact(get(base_vpm, v)));
  vertex_descriptor start_vertex = *vertices(kernel).begin();

  std::array<vertex_descriptor, 6> bbox_vertices;
  if(bbox_filtering){
    // We compute and store the vertices that realized the bbox
    struct Bbox_entry {
      std::size_t index;
      std::function<double(const EPoint_3&)> bound;
      std::function<double(const Bbox_3&)> value;
    };
    std::array<Bbox_entry,6> entries {{
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
  }


  // Get the planes and eventually shuffle them
  Three_point_cut_plane_traits<EK> kgt;
  auto oriented_side = kgt.oriented_side_3_object();
  auto orthogonal_vector = kgt.construct_orthogonal_vector_3_object();

  std::vector<face_descriptor> planes(face_range.begin(), face_range.end());
  if(shuffle_planes)
    std::shuffle(planes.begin(), planes.end(), std::default_random_engine(seed));

  // Cut iteratively the temporary kernel by halfspaces
  for(auto f: planes){
    auto h = halfedge(f, pm);
    Plane_3 plane(to_exact(get(vpm,source(h, pm))),
                  to_exact(get(vpm,target(h, pm))),
                  to_exact(get(vpm,target(next(h, pm), pm))));

    if(plane.is_degenerate())
      continue;

    if(bbox_filtering && vertices(kernel).size() >= 3 && faces(kernel).size()>1){
      // Early exit if the plane does not cut the bbox of the temporary kernel

      // By looking the sign of the plane value, we can check only two corners
      EVector_3 normal = orthogonal_vector(plane);
      // Look extreme corner according to the plane normal
      EPoint_3 corner( is_positive(normal.x())?bb3.xmax():bb3.xmin(),
                       is_positive(normal.y())?bb3.ymax():bb3.ymin(),
                       is_positive(normal.z())?bb3.zmax():bb3.zmin());
      if(oriented_side(plane, corner) != ON_POSITIVE_SIDE)
        continue;

      // Look the opposite corner
      EPoint_3 opposite_corner( is_positive(normal.x())?bb3.xmin():bb3.xmax(),
                                is_positive(normal.y())?bb3.ymin():bb3.ymax(),
                                is_positive(normal.z())?bb3.zmin():bb3.zmax());
      if(oriented_side(plane, opposite_corner) == ON_POSITIVE_SIDE)
        return PolygonMesh(); // empty

      start_vertex = clip_convex(kernel, plane, CGAL::parameters::clip_volume(true).
                                                                  geom_traits(kgt).
                                                                  do_not_triangulate_faces(true).
                                                                  vertex_point_map(kvpm).
                                                                  bounding_box(&bbox_vertices).
                                                                  starting_vertex_descriptor(start_vertex));
      if (is_empty(kernel)) return kernel;

      // update bbox, ( By looking which bbox_vertices have changed, it is possible to avoid recomputing all of them at each step )
      bb3 = get(kvpm, bbox_vertices[0]).bbox()+get(kvpm, bbox_vertices[1]).bbox()+get(kvpm, bbox_vertices[2]).bbox()+
            get(kvpm, bbox_vertices[3]).bbox()+get(kvpm, bbox_vertices[4]).bbox()+get(kvpm, bbox_vertices[5]).bbox();
    }
    else
    {
      start_vertex = clip_convex(kernel, plane, CGAL::parameters::clip_volume(true).
                                                                  geom_traits(kgt).
                                                                  do_not_triangulate_faces(true).
                                                                  vertex_point_map(kvpm).
                                                                  starting_vertex_descriptor(start_vertex));
      if (is_empty(kernel)) return kernel;
    }
  }

  if(used_to_find_a_point){
    bool require_strictly_inside = choose_parameter(get_parameter(np, internal_np::require_strictly_inside), true);

    // If we don't want a point on the surface and the kernel is degenerated, do nothing
    if(require_strictly_inside && ((vertices(kernel).size()<3) || (faces(kernel).size()==1)))
      return kernel;

    // Get the centroid
    EPoint_3 centroid(ORIGIN);
    for(auto v: vertices(kernel))
      centroid += EVector_3(ORIGIN, get(kvpm, v)) / vertices(kernel).size();

    // Approximate the centroid
    Point_3 double_centroid(to_double(centroid.x()), to_double(centroid.y()), to_double(centroid.z()));

    // Check if the approximate_centroid is inside the kernel
    bool is_valid = true;
    for(face_descriptor f: faces(kernel)){
      halfedge_descriptor h = halfedge(f, kernel);
      Plane_3 plane(get(kvpm,source(h, kernel)),
                    get(kvpm,target(h, kernel)),
                    get(kvpm,target(next(h, kernel), kernel)));
      if(oriented_side(plane, centroid) != ON_NEGATIVE_SIDE){
        is_valid = false;
        break;
      }
    }

    // If not, refine the centroid position
    if(!is_valid)
      centroid.exact();

    // Return the centroid
    *p = from_exact(centroid);
    return kernel;
  }

  // Convert points of the kernel to the type of the input mesh
  for(vertex_descriptor v : vertices(kernel))
    put(base_vpm, v, from_exact(get(kvpm, v)));
  return kernel;
};

} // end of namespace internal

 /**
  * \ingroup PMP_corefinement_grp
  *
  * \brief computes the kernel of the given mesh. The kernel is the set of all points that can see the entire surface of the mesh.
  * It is represented as a convex mesh and may be empty.
  *
  * The kernel is obtained by iteratively computing the intersection of the half-spaces defined by the faces of the mesh.
  *
  * The kernel may be degenerate. If the dimension of the kernel is 2, the output mesh consists of a single face. If the dimension of the kernel is 1, the output mesh contains only two isolated vertices. If the dimension of the kernel is 0, the output mesh contains only a single isolated vertex.
  *
  * @tparam PolygonMesh a model of `MutableFaceGraph`, `HalfedgeListGraph` and `FaceListGraph`.
  *                      An internal property map for `CGAL::vertex_point_t` must be available.
  *
  * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
  *
  * @param pm input surface mesh
  * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
  *
  * \cgalNamedParamsBegin
  *
  *   \cgalParamNBegin{vertex_point_map}
  *     \cgalParamDescription{a property map associating points to the vertices of `pm`}
  *     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
  *                    as key type and `%Point_3` as value type}
  *     \cgalParamDefault{`boost::get(CGAL::vertex_point, pm)`}
  *   \cgalParamNEnd
  *
  *   \cgalParamNBegin{use_bounding_box_filtering}
  *     \cgalParamDescription{Enables the use of the bounding box of the temporary kernel to compute the intersection of a plane with it, improving runtime in most scenarios.}
  *     \cgalParamType{bool}
  *     \cgalParamDefault{true}
  *   \cgalParamNEnd
  *
  *   \cgalParamNBegin{shuffle_planes}
  *     \cgalParamDescription{If set to `true`, the planes are considered in a random order to compute the kernel, improving runtime in most scenarios.}
  *     \cgalParamType{bool}
  *     \cgalParamDefault{true}
  *   \cgalParamNEnd
  *
  *   \cgalParamNBegin{random_seed}
  *     \cgalParamDescription{The seed use by the shuffle option}
  *     \cgalParamType{unsigned int}
  *     \cgalParamDefault{The seed of `std::random_device()`}
  *   \cgalParamNEnd
  *
  *   \cgalParamNBegin{starting_cube}
  *     \cgalParamDescription{
  *       The bounding box used to compute the kernel. The output is clipped based on this
  *       box, which must contain the kernel to guarantee a correct result.
  *     }
  *     \cgalParamType{`CGAL::Bbox_3`}
  *     \cgalParamDefault{`CGAL::Polygon_mesh_processing::bbox(pm)`}
  *   \cgalParamNEnd
  *
  *   \cgalParamNBegin{allow_non_manifold_non_watertight_input}
  *     \cgalParamDescription{
  *       If set to `true`, the input mesh is allowed to be non-manifold at vertices
  *       and/or to have boundaries. In this case, the kernel may be
  *       unbounded. The output is clipped by the bounding box provided by the `starting_cube` parameter.
  *     }
  *     \cgalParamType{bool}
  *     \cgalParamDefault{false}
  *   \cgalParamNEnd
  *
  *   \cgalParamNBegin{visitor}
  *     \cgalParamDescription{a visitor used to track the creation of new faces, edges, and faces.
  *                           Note that as there is no mesh associated with `plane`,
  *                           `boost::graph_traits<PolygonMesh>::null_halfedge()` and `boost::graph_traits<PolygonMesh>::null_face()` will be used when calling
  *                           functions of the visitor expecting a halfedge or a face from `plane`. Similarly, `pm` will be used as the mesh of `plane`.}
  *     \cgalParamType{a class model of `PMPCorefinementVisitor`}
  *     \cgalParamDefault{`Corefinement::Default_visitor<PolygonMesh>`}
  *   \cgalParamNEnd
  *
  * \cgalNamedParamsEnd
  *
  * @pre Unless the parameter `allow_non_manifold_non_watertight_input` is set to `true`. The input mesh is required to be closed, two-manifold, and free of self-intersections in order to ensure a correct result.
  *
  * @return A PolygonMesh representing the kernel of the input mesh, which may be empty.
  */
template <class PolygonMesh,
          class NamedParameters = parameters::Default_named_parameters>
PolygonMesh
kernel(const PolygonMesh& pm,
       const NamedParameters& np = parameters::default_values())
{
  return internal::kernel(pm, faces(pm), np);
}

 /**
  * \ingroup PMP_corefinement_grp
  *
  * \brief indicates whether the kernel of the given mesh is empty. The kernel is defined as the set of all points that can see the entire surface of the mesh.
  *
  * @tparam PolygonMesh a model of `MutableFaceGraph`, `HalfedgeListGraph` and `FaceListGraph`.
  *                      An internal property map for `CGAL::vertex_point_t` must be available.
  *
  * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
  *
  * @param pm input surface mesh
  * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
  *
  * \cgalNamedParamsBegin
  *
  *   \cgalParamNBegin{vertex_point_map}
  *     \cgalParamDescription{a property map associating points to the vertices of `pm`}
  *     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
  *                    as key type and `%Point_3` as value type}
  *     \cgalParamDefault{`boost::get(CGAL::vertex_point, pm)`}
  *   \cgalParamNEnd
  *
  *   \cgalParamNBegin{use_bounding_box_filtering}
  *     \cgalParamDescription{Enables the use of the bounding box of the temporary kernel to compute the intersection of a plane with it, improving runtime in most scenarios.}
  *     \cgalParamType{bool}
  *     \cgalParamDefault{true}
  *   \cgalParamNEnd
  *
  *   \cgalParamNBegin{shuffle_planes}
  *     \cgalParamDescription{If set to `true`, the planes are considered in random order to compute the kernel, improving runtime in most scenarios.}
  *     \cgalParamType{bool}
  *     \cgalParamDefault{true}
  *   \cgalParamNEnd
  *
  *   \cgalParamNBegin{random_seed}
  *     \cgalParamDescription{The seed use by the shuffle option}
  *     \cgalParamType{unsigned int}
  *     \cgalParamDefault{The seed of std::random_device()}
  *   \cgalParamNEnd
  *
  *   \cgalParamNBegin{starting_cube}
  *     \cgalParamDescription{
  *       The bounding box used to compute the kernel. The output is clipped based on this
  *       box, which must contain the kernel to guarantee a correct result.
  *     }
  *     \cgalParamType{`CGAL::Bbox_3`}
  *     \cgalParamDefault{`CGAL::Polygon_mesh_processing::bbox(pm)`}
  *   \cgalParamNEnd
  *
  *   \cgalParamNBegin{allow_non_manifold_non_watertight_input}
  *     \cgalParamDescription{
  *       If set to `true`, the input mesh is allowed to be non-manifold at vertices
  *       and/or to have boundaries. In this case, the kernel may be
  *       unbounded and the output is clipped by the bounding box provided by the `starting_cube` parameter.
  *     }
  *     \cgalParamType{bool}
  *     \cgalParamDefault{false}
  *   \cgalParamNEnd
  *
  * \cgalNamedParamsEnd
  *
  * @pre Unless the parameter `allow_non_manifold_non_watertight_input` is set to `true`. The input mesh is required to be closed, two-manifold, and free of self-intersections in order to ensure a correct result.
  *
  * @return bool
  */
template <class PolygonMesh,
          class NamedParameters = parameters::Default_named_parameters>
bool is_kernel_empty(const PolygonMesh& pm,
                     const NamedParameters& np = parameters::default_values())
{
  // TODO look if it's faster to compute with the dual instead (specifically in the none empty case)
  return is_empty(kernel(pm, np));
}

/**
  * \ingroup PMP_corefinement_grp
  *
  * \brief returns a point inside the kernel of the given mesh. The kernel is defined as the set of all points that can see the entire surface of the mesh.
  *
  * @tparam PolygonMesh a model of `MutableFaceGraph`, `HalfedgeListGraph` and `FaceListGraph`.
  *                      An internal property map for `CGAL::vertex_point_t` must be available.
  *
  * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
  *
  * @param pm input surface mesh
  * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
  *
  * \cgalNamedParamsBegin
  *
  *   \cgalParamNBegin{vertex_point_map}
  *     \cgalParamDescription{a property map associating points to the vertices of `pm`}
  *     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
  *                    as key type and `%Point_3` as value type}
  *     \cgalParamDefault{`boost::get(CGAL::vertex_point, pm)`}
  *   \cgalParamNEnd
  *
  *   \cgalParamNBegin{use_bounding_box_filtering}
  *     \cgalParamDescription{Enables the use of the bounding box of the temporary kernel to compute the intersection of a plane with it, improving runtime in most scenarios. }
  *     \cgalParamType{bool}
  *     \cgalParamDefault{true}
  *   \cgalParamNEnd
  *
  *   \cgalParamNBegin{shuffle_planes}
  *     \cgalParamDescription{If set to `true`, the planes are considered in a random order to compute the kernel, improving runtime in most scenarios }
  *     \cgalParamType{bool}
  *     \cgalParamDefault{true}
  *   \cgalParamNEnd
  *
  *   \cgalParamNBegin{starting_cube}
  *     \cgalParamDescription{
  *       The bounding box used to compute the kernel. The output is clipped based on this
  *       box, which must contain the kernel to guarantee a correct result.
  *     }
  *     \cgalParamType{`CGAL::Bbox_3`}
  *     \cgalParamDefault{`CGAL::Polygon_mesh_processing::bbox(pm)`}
  *   \cgalParamNEnd
  *
  *   \cgalParamNBegin{allow_non_manifold_non_watertight_input}
  *     \cgalParamDescription{
  *       If set to `true`, the input mesh is allowed to be non-manifold at vertices
  *       and/or to have boundaries. In this case, the kernel may be
  *       unbounded. The output is clipped by the bounding box provided by the `starting_cube` parameter.
  *     }
  *     \cgalParamType{bool}
  *     \cgalParamDefault{false}
  *   \cgalParamNEnd
  *
  *   \cgalParamNBegin{random_seed}
  *     \cgalParamDescription{The seed use by the shuffle option}
  *     \cgalParamType{unsigned int}
  *     \cgalParamDefault{The seed of std::random_device()}
  *   \cgalParamNEnd
  *
  * \cgalParamNBegin{require_strictly_inside}
  *   \cgalParamDescription{
  *     If set to `true`, the returned point is required to lie strictly inside
  *     the mesh and not on its boundary. If the mesh is degenerate, this requirement
  *     cannot be satisfied and no point is returned.}
  *   \cgalParamType{bool}
  *   \cgalParamDefault{true}
  * \cgalParamNEnd
  *
  * \cgalNamedParamsEnd
  *
  * @pre Unless the parameter `allow_non_manifold_non_watertight_input` is set to `true`. The input mesh is required to be closed, two-manifold, and free of self-intersections in order to ensure a correct result.
  *
  * @return std::optional<`%Point_3`>
  */
template <class PolygonMesh,
          class NamedParameters = parameters::Default_named_parameters>
#ifdef DOXYGEN_RUNNING
std::optional<Point_3>
#else
std::optional<typename GetGeomTraits<PolygonMesh, NamedParameters>::type::Point_3>
#endif
kernel_point(const PolygonMesh& pm,
             const NamedParameters& np = parameters::default_values())
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  bool require_strictly_inside = choose_parameter(get_parameter(np, internal_np::require_strictly_inside), true);

  std::optional<typename GetGeomTraits<PolygonMesh, NamedParameters>::type::Point_3> res;
  PolygonMesh k = internal::kernel(pm, faces(pm), np, true, &res);

  // If the kernel is empty or degenerated with strictly inside option, return empty
  if(is_empty(k) || (require_strictly_inside && ((vertices(k).size()<3) || (faces(k).size()==1))))
    return std::nullopt;

  return res;
}


} } // end of CGAL::Polygon_mesh_processing

#endif // CGAL_POLYGON_MESH_PROCESSING_KERNEL_H
