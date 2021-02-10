// Copyright (c) 2015 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sebastien Loriot,
//                 Mael Rouxel-Labbé

#ifndef CGAL_POLYGON_MESH_PROCESSING_REPAIR_H
#define CGAL_POLYGON_MESH_PROCESSING_REPAIR_H

#include <CGAL/license/Polygon_mesh_processing/repair.h>

#include <CGAL/Polygon_mesh_processing/manifoldness.h>
#include <CGAL/Polygon_mesh_processing/repair_degeneracies.h>
#include <CGAL/Polygon_mesh_processing/repair_self_intersections.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Polygon_mesh_processing/measure.h>

#include <CGAL/boost/graph/iterator.h>

#include <vector>

namespace CGAL {
namespace Polygon_mesh_processing {

/// \ingroup PMP_repairing_grp
/// removes the isolated vertices from any polygon mesh.
/// A vertex is considered isolated if it is not incident to any simplex
/// of higher dimension.
///
/// @tparam PolygonMesh a model of `FaceListGraph` and `MutableFaceGraph`
///
/// @param pmesh the polygon mesh to be repaired
///
/// @return number of removed isolated vertices
///
template <class PolygonMesh>
std::size_t remove_isolated_vertices(PolygonMesh& pmesh)
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;
  std::vector<vertex_descriptor> to_be_removed;

  for(vertex_descriptor v : vertices(pmesh))
  {
    if (CGAL::halfedges_around_target(v, pmesh).first
      == CGAL::halfedges_around_target(v, pmesh).second)
      to_be_removed.push_back(v);
  }
  std::size_t nb_removed = to_be_removed.size();
  for(vertex_descriptor v : to_be_removed)
  {
    remove_vertex(v, pmesh);
  }
  return nb_removed;
}

/// \ingroup PMP_repairing_grp
///
/// removes connected components whose area or volume is under a certain threshold value.
///
/// Thresholds are provided via \ref bgl_namedparameters "Named Parameters". (see below).
/// If thresholds are not provided by the user, default values are computed as follows:
/// - the area threshold is taken as the square of one percent of the length of the diagonal
///   of the bounding box of the mesh.
/// - the volume threshold is taken as the third power of one percent of the length of the diagonal
///   of the bounding box of the mesh.
///
/// The area and volume of a connected component will always be positive values (regardless
/// of the orientation of the mesh).
///
/// As a consequence of the last sentence, the area or volume criteria can be disabled
/// by passing zero (`0`) as threshold value.
///
/// \tparam TriangleMesh a model of `FaceListGraph` and `MutableFaceGraph`
/// \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
///
/// \param tmesh the triangulated polygon mesh
/// \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
///
/// \cgalNamedParamsBegin
///   \cgalParamNBegin{vertex_point_map}
///     \cgalParamDescription{a property map associating points to the vertices of `tmesh`}
///     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
///                    as key type and `%Point_3` as value type}
///     \cgalParamDefault{`boost::get(CGAL::vertex_point, tmesh)`}
///     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
///                     must be available in `TriangleMesh`.}
///   \cgalParamNEnd
///
///   \cgalParamNBegin{geom_traits}
///     \cgalParamDescription{an instance of a geometric traits class}
///     \cgalParamType{a class model of `Kernel`}
///     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
///     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
///     \cgalParamExtra{Exact constructions kernels are not supported by this function.}
///   \cgalParamNEnd
///
///   \cgalParamNBegin{face_index_map}
///     \cgalParamDescription{a property map associating to each face of `tmesh` a unique index between `0` and `num_faces(tmesh) - 1`}
///     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%face_descriptor`
///                    as key type and `std::size_t` as value type}
///     \cgalParamDefault{an automatically indexed internal map}
///   \cgalParamNEnd
///
///   \cgalParamNBegin{area_threshold}
///     \cgalParamDescription{a fixed value such that only connected components whose area is larger than this value are kept}
///     \cgalParamType{`geom_traits::FT`}
///     \cgalParamDefault{1\% of the length of the diagonal of the axis-aligned bounding box of the mesh, squared}
///   \cgalParamNEnd
///
///   \cgalParamNBegin{volume_threshold}
///     \cgalParamDescription{a fixed value such that only connected components whose volume is
///                           larger than this value are kept (only applies to closed connected components)}
///     \cgalParamType{`geom_traits::FT`}
///     \cgalParamDefault{1\% of the length of the diagonal of the axis-aligned bounding box of the mesh, cubed}
///     \cgalParamExtra{The mesh must be closed.}
///   \cgalParamNEnd
///
///   \cgalParamNBegin{edge_is_constrained_map}
///     \cgalParamDescription{a property map containing the constrained-or-not status of each edge of `tmesh`}
///     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%edge_descriptor`
///                    as key type and `bool` as value type. It must be default constructible.}
///     \cgalParamDefault{a default property map where no edge is constrained}
///     \cgalParamExtra{A constrained edge can be split or collapsed, but not flipped, nor its endpoints moved by smoothing.}
///   \cgalParamNEnd
///
///   \cgalParamNBegin{dry_run}
///     \cgalParamDescription{If `true`, the mesh will not be altered, but the number of components
///                           that would be removed is returned.}
///     \cgalParamType{Boolean}
///     \cgalParamDefault{`false`}
///   \cgalParamNEnd
///
///   \cgalParamNBegin{output_iterator}
///     \cgalParamDescription{An output iterator to collect the faces that would be removed by the algorithm,
///                           when using the "dry run" mode (see parameter `dry_run`)}
///     \cgalParamType{a model of `OutputIterator` with value type `face_descriptor`}
///     \cgalParamDefault{unused}
///   \cgalParamNEnd
/// \cgalNamedParamsEnd
///
/// \return the number of connected components removed (ignoring isolated vertices).
///
template <typename TriangleMesh,
          typename NamedParameters>
std::size_t remove_connected_components_of_negligible_size(TriangleMesh& tmesh,
                                                           const NamedParameters& np)
{
  using parameters::choose_parameter;
  using parameters::is_default_parameter;
  using parameters::get_parameter;

  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor          halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor              face_descriptor;

  typedef typename GetGeomTraits<TriangleMesh, NamedParameters>::type              GT;
  typedef typename GT::FT                                                          FT;
  const GT traits = choose_parameter<GT>(get_parameter(np, internal_np::vertex_point));

  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters>::const_type    VPM;
  typedef typename boost::property_traits<VPM>::value_type                         Point_3;
  const VPM vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                                   get_const_property_map(CGAL::vertex_point, tmesh));

  typedef typename GetInitializedFaceIndexMap<TriangleMesh, NamedParameters>::type FaceIndexMap;
  FaceIndexMap fim = CGAL::get_initialized_face_index_map(tmesh, np);

  FT area_threshold = choose_parameter(get_parameter(np, internal_np::area_threshold), FT(-1));
  FT volume_threshold = choose_parameter(get_parameter(np, internal_np::volume_threshold), FT(-1));

  // If no threshold is provided, compute it as a % of the bbox
  const bool is_default_area_threshold = is_default_parameter(get_parameter(np, internal_np::area_threshold));
  const bool is_default_volume_threshold = is_default_parameter(get_parameter(np, internal_np::volume_threshold));

  const bool dry_run = choose_parameter(get_parameter(np, internal_np::dry_run), false);

  typedef typename internal_np::Lookup_named_param_def<internal_np::output_iterator_t,
                                                       NamedParameters,
                                                       Emptyset_iterator>::type Output_iterator;
  Output_iterator out = choose_parameter<Output_iterator>(get_parameter(np, internal_np::output_iterator));

#ifdef CGAL_PMP_DEBUG_SMALL_CC_REMOVAL
  std::cout << "default threshold? " << is_default_area_threshold << " " << is_default_volume_threshold << std::endl;
#endif

  FT bbox_diagonal = FT(0), threshold_value = FT(0);

  if(is_default_area_threshold || is_default_volume_threshold)
  {
    if(is_empty(tmesh))
      return 0;

    const Bbox_3 bb = bbox(tmesh, np);

    bbox_diagonal = FT(CGAL::sqrt(CGAL::square(bb.xmax() - bb.xmin()) +
                                  CGAL::square(bb.ymax() - bb.ymin()) +
                                  CGAL::square(bb.zmax() - bb.zmin())));
    threshold_value = bbox_diagonal / FT(100); // default filter is 1%

#ifdef CGAL_PMP_DEBUG_SMALL_CC_REMOVAL
    std::cout << "bb xmin xmax: " << bb.xmin() << " " << bb.xmax() << std::endl;
    std::cout << "bb ymin ymax: " << bb.ymin() << " " << bb.ymax() << std::endl;
    std::cout << "bb zmin zmax: " << bb.zmin() << " " << bb.zmax() << std::endl;
    std::cout << "bbox_diagonal: " << bbox_diagonal << std::endl;
    std::cout << "threshold_value: " << threshold_value << std::endl;
#endif
  }

  if(is_default_area_threshold)
    area_threshold = CGAL::square(threshold_value);

  if(is_default_volume_threshold)
    volume_threshold = CGAL::square(threshold_value);

  const bool use_areas = (is_default_area_threshold || area_threshold > 0);
  const bool use_volumes = (is_default_volume_threshold || volume_threshold > 0);

  if(!use_areas && !use_volumes)
    return 0;

  // Compute the connected components only once
  boost::vector_property_map<std::size_t, FaceIndexMap> face_cc(static_cast<unsigned>(num_faces(tmesh)), fim);
  std::size_t num = connected_components(tmesh, face_cc, np);

#ifdef CGAL_PMP_DEBUG_SMALL_CC_REMOVAL
  std::cout << num << " different connected components" << std::endl;
#endif

  if(!dry_run)
    CGAL::Polygon_mesh_processing::remove_isolated_vertices(tmesh);

  // Compute CC-wide and total areas/volumes
  FT total_area = 0;
  std::vector<FT> component_areas(num, 0);

  if(use_areas)
  {
    for(face_descriptor f : faces(tmesh))
    {
      const FT fa = face_area(f, tmesh, np);
      component_areas[face_cc[f]] += fa;
      total_area += fa;
    }

#ifdef CGAL_PMP_DEBUG_SMALL_CC_REMOVAL
  std::cout << "area threshold: " << area_threshold << std::endl;
  std::cout << "total area: " << total_area << std::endl;
#endif
  }

  // Volumes make no sense for CCs that are not closed
  std::vector<bool> cc_closeness(num, true);
  std::vector<FT> component_volumes(num, FT(0));

  if(use_volumes)
  {
    for(halfedge_descriptor h : halfedges(tmesh))
    {
      if(is_border(h, tmesh))
        cc_closeness[face_cc[face(opposite(h, tmesh), tmesh)]] = false;
    }

    typename GT::Compute_volume_3 cv3 = traits.compute_volume_3_object();
    Point_3 origin(0, 0, 0);

    for(face_descriptor f : faces(tmesh))
    {
      const std::size_t i = face_cc[f];
      if(!cc_closeness[i])
        continue;

      const FT fv = cv3(origin,
                        get(vpm, target(halfedge(f, tmesh), tmesh)),
                        get(vpm, target(next(halfedge(f, tmesh), tmesh), tmesh)),
                        get(vpm, target(prev(halfedge(f, tmesh), tmesh), tmesh)));

      component_volumes[i] += fv;
    }

    // negative volume means the CC was oriented inward
    FT total_volume = 0;
    for(std::size_t i=0; i<num; ++i)
    {
      component_volumes[i] = CGAL::abs(component_volumes[i]);
      total_volume += component_volumes[i];
    }

#ifdef CGAL_PMP_DEBUG_SMALL_CC_REMOVAL
    std::cout << "volume threshold: " << volume_threshold << std::endl;
    std::cout << "total volume: " << total_volume << std::endl;
#endif
  }

  std::size_t res = 0;
  std::vector<bool> is_to_be_removed(num, false);

  for(std::size_t i=0; i<num; ++i)
  {
#ifdef CGAL_PMP_DEBUG_SMALL_CC_REMOVAL
    std::cout << "CC " << i << " has area: " << component_areas[i]
              << " and volume: " << component_volumes[i] << std::endl;
#endif

    if((use_volumes && cc_closeness[i] && component_volumes[i] <= volume_threshold) ||
       (use_areas && component_areas[i] <= area_threshold))
    {
      is_to_be_removed[i] = true;
      ++res;
    }
  }

#ifdef CGAL_PMP_DEBUG_SMALL_CC_REMOVAL
  std::cout << "Removing " << res << " CCs" << std::endl;
#endif

  if(dry_run)
  {
    for(face_descriptor f : faces(tmesh))
      if(is_to_be_removed[face_cc[f]])
        *out++ = f;
  }
  else
  {
    std::vector<std::size_t> ccs_to_remove;
    for(std::size_t i=0; i<num; ++i)
      if(is_to_be_removed[i])
        ccs_to_remove.push_back(i);

    remove_connected_components(tmesh, ccs_to_remove, face_cc, np);
    CGAL_expensive_postcondition(is_valid_polygon_mesh(tmesh));
  }

  return res;
}

template <typename TriangleMesh>
std::size_t remove_connected_components_of_negligible_size(TriangleMesh& tmesh)
{
  return remove_connected_components_of_negligible_size(tmesh, parameters::all_default());
}

} // namespace Polygon_mesh_processing
} // namespace CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_REPAIR_H
