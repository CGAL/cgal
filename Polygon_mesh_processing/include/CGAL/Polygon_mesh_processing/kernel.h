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
#include <CGAL/Polygon_mesh_processing/internal/Three_point_cut_plane_traits.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Cartesian_converter.h>

#include <boost/property_map/property_map.hpp>

namespace CGAL {
namespace Polygon_mesh_processing {
namespace internal {

template <typename PolygonMesh,
          typename FaceRange,
          typename NamedParameters = parameters::Default_named_parameters,
          typename NamedParametersOut = parameters::Default_named_parameters>
void
kernel(const FaceRange& face_range,
       const PolygonMesh& pm,
       PolygonMesh& kernel,
       const NamedParameters& np = parameters::default_values(),
       const NamedParametersOut& np_out = parameters::default_values(),
       bool used_to_find_a_point = false,
       std::optional<typename GetGeomTraits<PolygonMesh, NamedParameters>::type::Point_3> *p = nullptr)
{
  using parameters::choose_parameter;
  using parameters::get_parameter;
  using parameters::is_default_parameter;

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
  auto vpm_out = choose_parameter(get_parameter(np_out, internal_np::vertex_point),
                                  get_property_map(vertex_point, kernel));

  using DefaultF2FMap = Constant_property_map<face_descriptor, face_descriptor>;

  constexpr bool is_face_to_face_map = !parameters::is_default_parameter<NamedParametersOut, internal_np::face_to_face_map_t>::value;
  auto f2f_map = choose_parameter<DefaultF2FMap>(get_parameter(np_out, internal_np::face_to_face_map));

  using Point_3 = typename GT::Point_3;
  using EPoint_3 = typename EK::Point_3;
  using EVector_3 = typename EK::Vector_3;
  using Plane_3 = typename Three_point_cut_plane_traits<EK>::Plane_3;

  using KernelPointMap = typename boost::property_map<PolygonMesh, dynamic_vertex_property_t<EPoint_3> >::type;

  bool bbox_filtering = choose_parameter(get_parameter(np, internal_np::use_bounding_box_filtering), true);
  bool shuffle_planes = choose_parameter(get_parameter(np, internal_np::shuffle_planes), true);
  bool check_euler_characteristic = !choose_parameter(get_parameter(np, internal_np::allow_open_input), false) && std::size_t(std::distance(face_range.begin(), face_range.end()))==faces(pm).size();
  unsigned seed = choose_parameter(get_parameter(np, internal_np::random_seed), unsigned(-1));
  auto rng = is_default_parameter<NamedParameters, internal_np::random_seed_t>::value ? std::default_random_engine(): std::default_random_engine(seed);

  // Immediate exit if the input is well-formed and not of genus zero
  if(check_euler_characteristic && (vertices(pm).size() - edges(pm).size() + faces(pm).size() != 2)){
    clear(kernel);
    return;
  }

  // Build the starting cube
  KernelPointMap kvpm = get(CGAL::dynamic_vertex_property_t<EPoint_3>(), kernel);
  if(is_empty(kernel))
    make_hexahedron(bbox(pm), kernel, parameters::vertex_point_map(vpm_out));
  for(vertex_descriptor v: vertices(kernel))
    put(kvpm, v, to_exact(get(vpm_out, v)));
  Bbox_3 bb3 = bbox(kernel);
  vertex_descriptor start_vertex = *vertices(kernel).begin();
  if constexpr(is_face_to_face_map)
    for(face_descriptor f: faces(kernel))
      put(f2f_map, f, BGT::null_face());

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


  // Get the planes and possibly shuffle them
  Three_point_cut_plane_traits<EK> kgt;
  auto oriented_side = kgt.oriented_side_3_object();
  auto orthogonal_vector = kgt.construct_orthogonal_vector_3_object();

  std::vector<face_descriptor> planes(face_range.begin(), face_range.end());
  if(shuffle_planes)
    std::shuffle(planes.begin(), planes.end(), rng);

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
      if(oriented_side(plane, opposite_corner) == ON_POSITIVE_SIDE){
        clear(kernel); // empty
        return;
      }

      if constexpr(is_face_to_face_map)
        start_vertex = clip_convex(kernel, plane, CGAL::parameters::clip_volume(true).
                                                                    geom_traits(kgt).
                                                                    do_not_triangulate_faces(true).
                                                                    vertex_point_map(kvpm).
                                                                    bounding_box(&bbox_vertices).
                                                                    starting_vertex_descriptor(start_vertex).
                                                                    face_to_face_map(f2f_map),
                                                                    f);
      else
        start_vertex = clip_convex(kernel, plane, CGAL::parameters::clip_volume(true).
                                                                    geom_traits(kgt).
                                                                    do_not_triangulate_faces(true).
                                                                    vertex_point_map(kvpm).
                                                                    bounding_box(&bbox_vertices).
                                                                    starting_vertex_descriptor(start_vertex));
      if (is_empty(kernel)) return;

      // update bbox, ( By looking which bbox_vertices have changed, it is possible to avoid recomputing all of them at each step )
      bb3 = get(kvpm, bbox_vertices[0]).bbox()+get(kvpm, bbox_vertices[1]).bbox()+get(kvpm, bbox_vertices[2]).bbox()+
            get(kvpm, bbox_vertices[3]).bbox()+get(kvpm, bbox_vertices[4]).bbox()+get(kvpm, bbox_vertices[5]).bbox();
    }
    else
    {
      if constexpr(is_face_to_face_map)
        start_vertex = clip_convex(kernel, plane, CGAL::parameters::clip_volume(true).
                                                                    geom_traits(kgt).
                                                                    do_not_triangulate_faces(true).
                                                                    vertex_point_map(kvpm).
                                                                    starting_vertex_descriptor(start_vertex).
                                                                    face_to_face_map(f2f_map),
                                                                    f);
      else
        start_vertex = clip_convex(kernel, plane, CGAL::parameters::clip_volume(true).
                                                                    geom_traits(kgt).
                                                                    do_not_triangulate_faces(true).
                                                                    vertex_point_map(kvpm).
                                                                    starting_vertex_descriptor(start_vertex));
      if (is_empty(kernel)) return;
    }
  }

  if(used_to_find_a_point){
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
    return;
  }

  // Convert points of the kernel to the type of the input mesh
  for(vertex_descriptor v : vertices(kernel))
    put(vpm_out, v, from_exact(get(kvpm, v)));
};

} // end of namespace internal

/**
  * \ingroup PMP_corefinement_grp
  *
  * \brief computes the kernel of the given faces of a polygon mesh.
  *
  * The kernel is defined as the convex polyhedron that is the intersection
  * of all the halfspaces on the negative side of the oriented planes defined by a range of faces
  * of the input mesh. The kernel may be empty or degenerate to a lower-dimensional convex shape.
  *
  * In the implementation, a starting shape is iteratively clipped by the faces.
  * By default, the bounding box of the input mesh is used as starting shape.
  * However, the parameter `out` may be non-empty: In this case, it must be a convex polyhedron and will be used as starting shape.
  *
  * The algorithm assumes that the faces of the input range form a closed surface as to perform a quick exit if the genus is non-zero.
  * This precondition can be relaxed using the named parameter `allow_open_input`.
  * In that case, the resulting kernel may contain faces of the starting shape.
  *
  * In case of a degenerate kernel:
  * <ul>
  *   <li>If the dimension of the kernel is `2` (i.e., the kernel is a convex polygon in 3D), the output mesh consists of a single face.</li>
  *   <li>If the dimension of the kernel is `1` (i.e., the kernel is a line segment), the output mesh consists two isolated vertices.</li>
  *   <li>If the dimension of the kernel is `0` (i.e., the kernel is a single point), the output mesh contains one isolated vertex.</li>
  * </ul>
  *
  * @tparam FaceRange a model of `ConstRange` with `boost::graph_traits<PolygonMesh>::%face_descriptor` as value type
  * @tparam PolygonMesh a model of `VertexListGraph`, `HalfedgeListGraph` and `FaceListGraph`
  * @tparam PolygonMeshOut a model of `MutableFaceGraph`, `VertexListGraph` and `FaceListGraph`
  *
  * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
  * @tparam NamedParametersOut a sequence of \ref bgl_namedparameters "Named Parameters"
  *
  * @param face_range the range of faces used
  * @param pm input surface mesh
  * @param out output surface mesh
  * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
  *
  * \cgalNamedParamsBegin
  *   \cgalParamNBegin{allow_open_input}
  *     \cgalParamDescription{If set to `true`, the input mesh is allowed to have boundaries.}
  *     \cgalParamType{Boolean}
  *     \cgalParamDefault{`false`}
  *   \cgalParamNEnd
  *
  *   \cgalParamNBegin{vertex_point_map}
  *     \cgalParamDescription{a property map associating points to the vertices of `pm`}
  *     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
  *                    as key type and `%Point_3` as value type}
  *     \cgalParamDefault{`boost::get(CGAL::vertex_point, pm)`}
  *     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t` must be available in PolygonMesh. }
  *   \cgalParamNEnd
  *
  *   \cgalParamNBegin{geom_traits}
  *     \cgalParamDescription{an instance of a geometric traits class}
  *     \cgalParamType{a class model of `Kernel`}
  *     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
  *   \cgalParamNEnd
  *
  *   \cgalParamNBegin{random_seed}
  *     \cgalParamDescription{is used to initialize the random number generator of the algorithm.}
  *     \cgalParamType{unsigned int}
  *     \cgalParamDefault{use `std::default_random_engine()`}
  *   \cgalParamNEnd
  *
  *   \cond SKIP_IN_MANUAL
  *
  *   \cgalParamNBegin{use_bounding_box_filtering}
  *     \cgalParamDescription{Enables the use of the bounding box of the temporary kernel to compute the intersection of a plane with it, improving runtime in most scenarios.}
  *     \cgalParamType{Boolean}
  *     \cgalParamDefault{`true`}
  *   \cgalParamNEnd
  *
  *   \cgalParamNBegin{shuffle_planes}
  *     \cgalParamDescription{If set to `true`, the planes are considered in a random order to compute the kernel, improving runtime in most scenarios.}
  *     \cgalParamType{Boolean}
  *     \cgalParamDefault{`true`}
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
  *   \endcond
  * \cgalNamedParamsEnd
  *
  * @param np_out an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
  *
  * \cgalNamedParamsBegin
  *   \cgalParamNBegin{vertex_point_map}
  *     \cgalParamDescription{a property map associating points to the vertices of `out`}
  *     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
  *                    as key type and `%Point_3` as value type}
  *     \cgalParamDefault{`boost::get(CGAL::vertex_point, out)`}
  *     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t` must be available in PolygonMesh. }
  *   \cgalParamNEnd
  *   \cgalParamNBegin{face_to_face_map}
  *     \cgalParamDescription{a property map storing, for each face of the output mesh, a face of the input mesh that defined the clipping plane that created it
                              (or `boost::graph_traits<PolygonMeshOut>::%null_face` if the face belongs to the starting shape)}
  *     \cgalParamType{a class model of `ReadWritePropertyMap` with
  *                   `boost::graph_traits<PolygonMeshOut>::%face_descriptor` as key type and
  *                   `boost::graph_traits<PolygonMesh>::%face_descriptor` as value type}
  *     \cgalParamDefault{unused}
  *   \cgalParamNEnd
  * \cgalNamedParamsEnd
  */
template <typename FaceRange,
          typename PolygonMesh,
          typename NamedParameters = parameters::Default_named_parameters,
          typename NamedParametersOut = parameters::Default_named_parameters>
void
kernel(const FaceRange& face_range,
       const PolygonMesh& pm,
       PolygonMesh& out,
       const NamedParameters& np = parameters::default_values(),
       const NamedParametersOut& np_out = parameters::default_values())
{
  internal::kernel(face_range, pm, out, np, np_out);
}

/**
  * \ingroup PMP_corefinement_grp
  *
  * \brief computes the kernel of the given polygon mesh.
  *
  * This is a convenience overload that calls the overload above
  * on all faces of the mesh.
  */
template <typename PolygonMesh,
          typename NamedParameters = parameters::Default_named_parameters,
          typename NamedParametersOut = parameters::Default_named_parameters>
void
kernel(const PolygonMesh& pm,
       PolygonMesh& out,
       const NamedParameters& np = parameters::default_values(),
       const NamedParametersOut& np_out = parameters::default_values())
{
  kernel(faces(pm), pm, out, np, np_out);
}

/**
  * \ingroup PMP_corefinement_grp
  *
  * \brief indicates whether the kernel of the given faces of a polygon mesh is empty.
  *
  * The kernel is defined as the convex polyhedron that is the intersection
  * of all the halfspaces on the negative side of the oriented planes defined by a range of faces
  * of the input mesh.
  *
  * See `CGAL::Polygon_mesh_processing::kernel()` for a comprehensive description of the parameters.
  */
template <typename FaceRange,
          typename PolygonMesh,
          typename CGAL_NP_TEMPLATE_PARAMETERS>
bool has_empty_kernel(const FaceRange& face_range,
                     const PolygonMesh& pm,
                     const CGAL_NP_CLASS& np = parameters::default_values())
{
  PolygonMesh k;
  kernel(face_range, pm, k, np);
  return is_empty(k);
}

/**
  * \ingroup PMP_corefinement_grp
  *
  * \brief indicates whether the kernel of the given polygon mesh is empty.
  *
  * The kernel is defined as the convex polyhedron that is the intersection
  * of all the halfspaces on the negative side of the oriented planes defined by a range of faces
  * of the input mesh.
  *
  * See `CGAL::Polygon_mesh_processing::kernel()` for a comprehensive description of the parameters.
  */
template <typename PolygonMesh,
          typename CGAL_NP_TEMPLATE_PARAMETERS>
bool has_empty_kernel(const PolygonMesh& pm,
                     const CGAL_NP_CLASS& np = parameters::default_values())
{
  return has_empty_kernel(faces(pm), pm, np);
}

/**
  * \ingroup PMP_corefinement_grp
  *
  * \brief returns a point inside the kernel of the given faces of a polygon mesh.
  *
  * The kernel is defined as the convex polyhedron that is the intersection
  * of all the halfspaces on the negative side of the oriented planes defined by a range of faces
  * of the input mesh.
  *
  * See `CGAL::Polygon_mesh_processing::kernel()` for a comprehensive description of the parameters.
  *
  * \return `std::nullopt` if and only if the kernel is empty.
  */
template <typename FaceRange,
          typename PolygonMesh,
          typename NamedParameters = parameters::Default_named_parameters>
#ifdef DOXYGEN_RUNNING
std::optional<Point_3>
#else
std::optional<typename GetGeomTraits<PolygonMesh, NamedParameters>::type::Point_3>
#endif
kernel_point(const FaceRange& face_range,
             const PolygonMesh& pm,
             const NamedParameters& np = parameters::default_values())
{
  std::optional<typename GetGeomTraits<PolygonMesh, NamedParameters>::type::Point_3> res;
  PolygonMesh k;
  internal::kernel(face_range, pm, k, np, parameters::default_values(), true, &res);

  // If the kernel is empty or degenerated with strictly inside option, return empty
  if(is_empty(k))
    return std::nullopt;

  return res;
}

/**
  * \ingroup PMP_corefinement_grp
  *
  * \brief returns a point inside the kernel of the given polygon mesh.
  *
  * The kernel is defined as the convex polyhedron that is the intersection
  * of all the halfspaces on the negative side of the oriented planes defined by a range of faces
  * of the input mesh.
  *
  * See `CGAL::Polygon_mesh_processing::kernel()` for a comprehensive description of the parameters.
  *
  * \return `std::nullopt` if and only if the kernel is empty.
  */
template <typename PolygonMesh,
          typename CGAL_NP_TEMPLATE_PARAMETERS>
#ifdef DOXYGEN_RUNNING
std::optional<Point_3>
#else
std::optional<typename GetGeomTraits<PolygonMesh, CGAL_NP_CLASS>::type::Point_3>
#endif
kernel_point(const PolygonMesh& pm,
             const CGAL_NP_CLASS& np = parameters::default_values())
{
  return kernel_point(faces(pm), pm, np);
}


} } // end of CGAL::Polygon_mesh_processing

#endif // CGAL_POLYGON_MESH_PROCESSING_KERNEL_H
