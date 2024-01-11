// Copyright (c) 2019-2022 Google LLC (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Cedric Portaneri,
//                 Mael Rouxel-Labbé
//
#ifndef CGAL_ALPHA_WRAP_3_TEST_ALPHA_WRAP_VALIDATION_H
#define CGAL_ALPHA_WRAP_3_TEST_ALPHA_WRAP_VALIDATION_H

#include <CGAL/license/Alpha_wrap_3.h>

#include <CGAL/boost/graph/helpers.h>
#include <CGAL/Kernel_traits.h>
#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/Polygon_mesh_processing/intersection.h>
#include <CGAL/Polygon_mesh_processing/manifoldness.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/repair.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <CGAL/Side_of_triangle_mesh.h>

namespace CGAL {
namespace Alpha_wraps_3 {
namespace internal {

// @todo this does not detect non-manifold edges with manifold vertices
// (but PMP::self_intersections will catch it)
template <typename TriangleMesh>
bool is_combinatorially_non_manifold(const TriangleMesh& mesh)
{
  namespace PMP = CGAL::Polygon_mesh_processing;

  using halfedge_descriptor = typename boost::graph_traits<TriangleMesh>::halfedge_descriptor;

  std::vector<halfedge_descriptor> non_manifold_cones;
  PMP::non_manifold_vertices(mesh, std::back_inserter(non_manifold_cones));

  if(!non_manifold_cones.empty())
    return true;

  return false;
}

template <typename TriangleMesh,
          typename NamedParameters = parameters::Default_named_parameters>
bool has_degenerated_faces(const TriangleMesh& mesh,
                           const NamedParameters& np = CGAL::parameters::default_values())
{
  namespace PMP = CGAL::Polygon_mesh_processing;

  for(auto f : faces(mesh))
    if(PMP::is_degenerate_triangle_face(f, mesh, np))
      return true;

  return false;
}

// Edge length is bounded by twice the circumradius
template <typename TriangleMesh,
          typename NamedParameters = parameters::Default_named_parameters>
bool check_edge_length(const TriangleMesh& output_mesh,
                       const double alpha,
                       const NamedParameters& np = CGAL::parameters::default_values())
{
  const auto sq_alpha_bound = 4 * square(alpha);
  for(auto e : edges(output_mesh))
  {
    const auto sqd = Polygon_mesh_processing::squared_edge_length(e, output_mesh, np);
    if(sqd > sq_alpha_bound) // alpha is the circumradius
    {
#ifdef CGAL_AW3_DEBUG
      std::cerr << "Error: " << sqd << " greater than " << sq_alpha_bound << std::endl;
      std::cerr << get(CGAL::vertex_point, output_mesh, source(e, output_mesh)) << std::endl;
      std::cerr << get(CGAL::vertex_point, output_mesh, target(e, output_mesh)) << std::endl;
#endif
      return false;
    }
  }

  return true;
}

template <typename ConcurrencyTag = CGAL::Sequential_tag,
          typename TriangleMesh, typename FT,
          typename InputNamedParameters = parameters::Default_named_parameters,
          typename OutputNamedParameters = parameters::Default_named_parameters>
bool has_expected_Hausdorff_distance(const TriangleMesh& wrap,
                                     const TriangleMesh& input,
                                     const FT alpha, const FT offset,
                                     const InputNamedParameters& in_np = parameters::default_values(),
                                     const OutputNamedParameters& out_np = parameters::default_values())
{
  namespace PMP = CGAL::Polygon_mesh_processing;

  using face_descriptor = typename boost::graph_traits<TriangleMesh>::face_descriptor;

  std::vector<std::pair<face_descriptor, face_descriptor> > fpairs;

  const FT bound = 0.01 * (std::min)(alpha, offset);
  const FT d = PMP::bounded_error_Hausdorff_distance<ConcurrencyTag>(
                 wrap, input, bound, in_np.output_iterator(std::back_inserter(fpairs)), out_np);

#ifdef CGAL_AW3_DEBUG
  std::cout << "Alpha: " << alpha << " Offset: " << offset << " Bound: " << bound << " Hausdorff_distance: " << d << std::endl;
  std::cout << "Maximum distance on faces " << fpairs.back().first << " " << fpairs.back().second << std::endl;
#endif

  return (d < (alpha + offset + bound));
}

template <typename TriangleMesh, typename NamedParameters = parameters::Default_named_parameters>
bool is_valid_wrap(const TriangleMesh& wrap,
                   const bool check_manifoldness,
                   const NamedParameters& np = parameters::default_values())
{
  namespace PMP = CGAL::Polygon_mesh_processing;

  if(is_empty(wrap))
  {
#ifdef CGAL_AW3_DEBUG
    std::cerr << "Error: empty wrap" << std::endl;
#endif
    return false;
  }

  if(!is_valid_polygon_mesh(wrap))
  {
#ifdef CGAL_AW3_DEBUG
    std::cerr << "Error: Invalid wrap mesh" << std::endl;
#endif
    return false;
  }

  if(!is_triangle_mesh(wrap))
  {
#ifdef CGAL_AW3_DEBUG
    std::cerr << "Error: Wrap is not triangulated" << std::endl;
#endif
    return false;
  }

  if(!is_closed(wrap))
  {
    if(check_manifoldness)
    {
#ifdef CGAL_AW3_DEBUG
      std::cerr << "Error: Wrap is not closed" << std::endl;
#endif
      return false;
    }
    else
    {
#ifdef CGAL_AW3_DEBUG
      std::cerr << "W: Wrap is not closed" << std::endl;
#endif
    }
  }
  else if(!PMP::does_bound_a_volume(wrap, np))
  {
#ifdef CGAL_AW3_DEBUG
    std::cerr << "Error: Wrap does not bound a volume" << std::endl;
#endif
    return false;
  }

  if(has_degenerated_faces(wrap, np))
  {
#ifdef CGAL_AW3_DEBUG
    std::cerr << "Error: Wrap has degenerate faces" << std::endl;
#endif
    return false;
  }

  if(is_combinatorially_non_manifold(wrap))
  {
#ifdef CGAL_AW3_DEBUG
    std::cerr << "Error: Wrap is combinatorially non-manifold" << std::endl;
#endif
    return false;
  }

  if(PMP::does_self_intersect(wrap, np))
  {
    if(check_manifoldness)
    {
#ifdef CGAL_AW3_DEBUG
      std::cerr << "Error: Wrap self-intersects" << std::endl;
#endif
      return false;
    }
#ifdef CGAL_AW3_DEBUG
    std::cerr << "W: Wrap self-intersects" << std::endl;
#endif
  }

  return true;
}

template <typename TriangleMesh, typename NamedParameters = parameters::Default_named_parameters>
bool is_valid_wrap(const TriangleMesh& wrap,
                   const NamedParameters& np = parameters::default_values())
{
  return is_valid_wrap(wrap, true /*consider manifoldness*/, np);
}

template <typename InputTriangleMesh, typename OutputTriangleMesh,
          typename InputNamedParameters = parameters::Default_named_parameters,
          typename OutputNamedParameters = parameters::Default_named_parameters>
bool is_outer_wrap_of_triangle_mesh(const OutputTriangleMesh& wrap,
                                    const InputTriangleMesh& input,
                                    const OutputNamedParameters& out_np = parameters::default_values(),
                                    const InputNamedParameters& in_np = parameters::default_values())
{
  namespace PMP = CGAL::Polygon_mesh_processing;

  using parameters::get_parameter;
  using parameters::choose_parameter;

  using IVPM = typename GetVertexPointMap<InputTriangleMesh, InputNamedParameters>::const_type;
  using OVPM = typename GetVertexPointMap<OutputTriangleMesh, OutputNamedParameters>::const_type;
  using K = typename GetGeomTraits<InputTriangleMesh, InputNamedParameters>::type;

//  CGAL::Rigid_triangle_mesh_collision_detection<TriangleMesh> collision_detection;
//  collision_detection.add_mesh(input);
//  collision_detection.add_mesh(wrap);

//  auto res = collision_detection.get_all_intersections(1);
//  if(res.size() != 0)
//  {
//    std::cerr << "Error: The wrap intersects the input mesh" << std::endl;
//    return EXIT_FAILURE;
//  }

  if(PMP::do_intersect(input, wrap, in_np, out_np))
  {
#ifdef CGAL_AW3_DEBUG
    std::cerr << "Error: The wrap intersects the input mesh" << std::endl;
#endif
    return false;
  }

  IVPM in_vpm = choose_parameter(get_parameter(in_np, internal_np::vertex_point),
                                 get_const_property_map(vertex_point, input));
  OVPM out_vpm = choose_parameter(get_parameter(out_np, internal_np::vertex_point),
                                  get_const_property_map(vertex_point, wrap));

  CGAL::Side_of_triangle_mesh<OutputTriangleMesh, K, OVPM> side_of_wrap(wrap, out_vpm);

  // @speed a single vertex per CC would be sufficient
  for(auto v : vertices(input))
  {
    if(side_of_wrap(get(in_vpm, v)) != CGAL::ON_BOUNDED_SIDE)
    {
#ifdef CGAL_AW3_DEBUG
      std::cerr << "Error: Part(s) of the input mesh are outside the wrap: " << get(in_vpm, v) << std::endl;
#endif
      return false;
    }
  }

  return true;
}

template <typename OutputTriangleMesh, typename InputTriangleMesh,
          typename OutputNamedParameters = parameters::Default_named_parameters,
          typename InputNamedParameters = parameters::Default_named_parameters>
bool is_valid_wrap_of_triangle_mesh(const OutputTriangleMesh& wrap,
                                    const InputTriangleMesh& input,
                                    const OutputNamedParameters& out_np = parameters::default_values(),
                                    const InputNamedParameters& in_np = parameters::default_values())
{
  if(!is_valid_wrap(wrap, out_np))
    return false;

  if(!is_outer_wrap_of_triangle_mesh(wrap, input, out_np, in_np))
    return false;

  return true;
}

template <typename TriangleMesh, typename PointRange, typename FaceRange,
          typename OutputNamedParameters = parameters::Default_named_parameters,
          typename InputNamedParameters = parameters::Default_named_parameters>
bool is_outer_wrap_of_triangle_soup(const TriangleMesh& wrap,
                                    PointRange points, // intentional copies
                                    FaceRange faces,
                                    const OutputNamedParameters& out_np = parameters::default_values(),
                                    const InputNamedParameters& in_np = parameters::default_values())
{
  namespace PMP = CGAL::Polygon_mesh_processing;

  // Make a mesh out of the soup
  PMP::repair_polygon_soup(points, faces);
  PMP::orient_polygon_soup(points, faces);
  CGAL_assertion(PMP::is_polygon_soup_a_polygon_mesh(faces));

  TriangleMesh mesh;
  PMP::polygon_soup_to_polygon_mesh(points, faces, mesh, in_np);

  return is_outer_wrap_of_triangle_mesh(wrap, mesh, out_np);
}

template <typename TriangleMesh, typename PointRange, typename FaceRange,
          typename OutputNamedParameters = parameters::Default_named_parameters,
          typename InputNamedParameters = parameters::Default_named_parameters>
bool is_valid_wrap_of_triangle_soup(const TriangleMesh& wrap,
                                    const PointRange& points,
                                    const FaceRange& faces,
                                    const OutputNamedParameters& out_np = parameters::default_values(),
                                    const InputNamedParameters& in_np = parameters::default_values())
{
  if(!is_valid_wrap(wrap, out_np))
    return false;

  if(!is_outer_wrap_of_triangle_soup(wrap, points, faces, out_np, in_np))
    return false;

  return true;
}

template <typename TriangleMesh, typename PointRange,
          typename OutputNamedParameters = parameters::Default_named_parameters,
          typename InputNamedParameters = parameters::Default_named_parameters>
bool is_outer_wrap_of_point_set(const TriangleMesh& wrap,
                                const PointRange& points,
                                const OutputNamedParameters& out_np = parameters::default_values(),
                                const InputNamedParameters& in_np = parameters::default_values())
{
  using parameters::get_parameter;
  using parameters::choose_parameter;

  using OVPM = typename GetVertexPointMap<TriangleMesh, OutputNamedParameters>::const_type;
  using IPM = typename GetPointMap<PointRange, InputNamedParameters>::const_type;
  using K = typename Kernel_traits<typename boost::property_traits<IPM>::value_type>::Kernel;

  OVPM out_vpm = choose_parameter(get_parameter(out_np, internal_np::vertex_point),
                                  get_const_property_map(vertex_point, wrap));
  IPM in_pm = choose_parameter<IPM>(get_parameter(in_np, internal_np::point_map));

  CGAL::Side_of_triangle_mesh<TriangleMesh, K, OVPM> side_of_wrap(wrap, out_vpm);

  // @speed a single vertex per CC would be sufficient
  for(const auto& p : points)
  {
    if(side_of_wrap(get(in_pm, p)) != CGAL::ON_BOUNDED_SIDE)
    {
#ifdef CGAL_AW3_DEBUG
      std::cerr << "Part(s) of the input mesh are outside the wrap: " << get(in_pm, p) << std::endl;
#endif
      return false;
    }
  }

  return true;
}

template <typename TriangleMesh, typename PointRange,
          typename OutputNamedParameters = parameters::Default_named_parameters,
          typename InputNamedParameters = parameters::Default_named_parameters>
bool is_valid_wrap_of_point_set(const TriangleMesh& wrap,
                                const PointRange& points,
                                const OutputNamedParameters& out_np = parameters::default_values(),
                                const InputNamedParameters& in_np = parameters::default_values())
{
  if(!is_valid_wrap(wrap, out_np))
    return false;

  if(!is_outer_wrap_of_point_set(wrap, points, out_np, in_np))
    return false;

  return true;
}

} // namespace internal
} // namespace Alpha_wraps_3
} // namespace CGAL

#endif // CGAL_ALPHA_WRAP_3_TEST_ALPHA_WRAP_VALIDATION_H
