// Copyright (c) 2017-2018 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Pierre Alliez and Lingjie Zhu


#ifndef CGAL_SURFACE_MESH_APPROXIMATION_VSA_MESH_APPROXIMATION_H
#define CGAL_SURFACE_MESH_APPROXIMATION_VSA_MESH_APPROXIMATION_H

#include <CGAL/license/Surface_mesh_approximation.h>


#include <CGAL/Variational_shape_approximation.h>
#include <CGAL/boost/graph/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <CGAL/property_map.h>

#include <boost/graph/graph_traits.hpp>
#include <boost/type_traits/is_same.hpp>

#include <iostream>
#include <limits>

namespace CGAL {
namespace Surface_mesh_approximation {

/// \ingroup PkgTSMARef
/// @brief Verbose level enumeration.
enum Verbose_level {
  /// Silent
  SILENT,
  /// Main steps
  MAIN_STEPS,
  /// Verbose
  VERBOSE
};

// the named parameter header being not documented the doc is put here for now
#ifdef DOXYGEN_RUNNING
namespace parameters {

/*! \ingroup vsa_namedparameters
 * This function is used when default parameters are just fine for approximation or meshing.
 */
unspecified_type all_default();

} // namespace parameters
#endif

/*!
 * \ingroup PkgTSMARef
 * @brief approximates the input mesh with plane proxies.
 * This function uses the Variational Shape Approximation algorithm described in \cgalCite{cgal:cad-vsa-04}
 * to approximate a triangle surface mesh, with indexed triangles as output.
 *
 * @tparam TriangleMesh model of `FaceListGraph`
 * @tparam NamedParameters a sequence of \ref vsa_namedparameters
 *
 * @param tm triangle surface mesh to be approximated
 * @param np optional sequence of \ref vsa_namedparameters among the ones listed below
 * @return `true` if the indexed triangles represent a 2-manifold, oriented surface mesh, and `false` otherwise.
 *
 * \cgalNamedParamsBegin{Approximation Named Parameters}
 *  \cgalParamBegin{geom_traits} a geometric traits class instance, model of Kernel.
 *    Exact constructions kernels are not supported by this function.
 *  \cgalParamEnd
 *  \cgalParamBegin{vertex_point_map} the property map with the points associated
 *    to the vertices of `tm`. Instance of a class model of `ReadablePropertyMap`.
 *  \cgalParamEnd
 *  \cgalParamBegin{verbose_level} set verbose level.
 *  \cgalParamEnd
 *  \cgalParamBegin{seeding_method} selection of seeding method.
 *  \cgalParamEnd
 *  \cgalParamBegin{max_number_of_proxies} maximum number of proxies to approximate the input mesh.
 *  \cgalParamEnd
 *  \cgalParamBegin{min_error_drop} minimum error drop of the approximation, expressed in ratio between two iterations of proxy addition.
 *  \cgalParamEnd
 *  \cgalParamBegin{number_of_relaxations} number of relaxation iterations interleaved within seeding.
 *  \cgalParamEnd
 *  \cgalParamBegin{number_of_iterations} number of partitioning and fitting iterations after seeding.
 *  \cgalParamEnd
 * \cgalNamedParamsEnd
 * \cgalNamedParamsBegin{Meshing Named Parameters}
 *  \cgalParamBegin{subdivision_ratio} chord subdivision ratio threshold to the chord length or average edge length.
 *  \cgalParamEnd
 *  \cgalParamBegin{relative_to_chord} if `true` the `subdivision_ratio` is the ratio of the
 *    furthest vertex distance to the chord length, otherwise is the average edge length.
 *  \cgalParamEnd
 *  \cgalParamBegin{with_dihedral_angle} if set to `true` the `subdivision_ratio` is weighted by dihedral angle.
 *  \cgalParamEnd
 *  \cgalParamBegin{optimize_anchor_location} if set to `true`, optimize the anchor locations.
 *  \cgalParamEnd
 *  \cgalParamBegin{pca_plane} set `true` if use PCA plane fitting, otherwise use the default area averaged plane parameters.
 *  \cgalParamEnd
 * \cgalNamedParamsEnd
 * \cgalNamedParamsBegin{Output Named Parameters}
 *  \cgalParamBegin{face_proxy_map} a `WritablePropertyMap` with
 * `boost::graph_traits<TriangleMesh>::%face_descriptor` as key and `std::size_t` as value type.
 * A proxy is a set of connected faces which are placed under the same proxy patch (see \cgalFigureRef{iterations}).
 * The proxy-ids are contiguous in range [0, number_of_proxies - 1].
 *  \cgalParamEnd
 *  \cgalParamBegin{proxies} output iterator over proxies.
 *  \cgalParamEnd
 *  \cgalParamBegin{anchors} output iterator over anchor points.
 *  \cgalParamEnd
 *  \cgalParamBegin{triangles} output iterator over indexed triangles.
 *  \cgalParamEnd
 * \cgalNamedParamsEnd
 */
template <typename TriangleMesh, typename NamedParameters>
bool approximate_triangle_mesh(const TriangleMesh &tm, const NamedParameters &np)
{
  using parameters::get_parameter;
  using parameters::choose_parameter;
  using parameters::is_default_parameter;

  typedef typename CGAL::GetGeomTraits<TriangleMesh, NamedParameters>::type Geom_traits;
  typedef typename Geom_traits::FT FT;

  typedef typename CGAL::GetVertexPointMap<TriangleMesh, NamedParameters>::type Vertex_point_map;
  Vertex_point_map point_pmap = choose_parameter(get_parameter(np, internal_np::vertex_point),
    get_property_map(vertex_point, const_cast<TriangleMesh &>(tm)));

  typedef CGAL::Variational_shape_approximation<TriangleMesh, Vertex_point_map> L21_approx;
  typedef typename L21_approx::Error_metric L21_metric;

  const Verbose_level vl = choose_parameter(
    get_parameter(np, internal_np::verbose_level), SILENT);

  const std::size_t number_of_faces = std::distance(faces(tm).first, faces(tm).second);
  const std::size_t number_of_vertices = std::distance(vertices(tm).first, vertices(tm).second);

  if (vl == MAIN_STEPS || vl == VERBOSE) {
    std::cout << "Variational shape approximation:"
      << "\n#f " << number_of_faces
      << "\n#v " << number_of_vertices << std::endl;
  }

  L21_metric metric(tm, point_pmap);
  L21_approx approx(tm, point_pmap, metric);

  // hierarchical seeding by default
  const Seeding_method method = choose_parameter(
    get_parameter(np, internal_np::seeding_method), HIERARCHICAL);
  const std::size_t max_nb_of_proxies = choose_parameter(
    get_parameter(np, internal_np::max_number_of_proxies), 0);
  const FT min_error_drop = choose_parameter(
    get_parameter(np, internal_np::min_error_drop), FT(0.0));
  const std::size_t nb_of_relaxations = choose_parameter(
    get_parameter(np, internal_np::number_of_relaxations), 5);

  if (vl == VERBOSE) {
    std::cout << (method == RANDOM ? "Random" :
      (method == INCREMENTAL ? "Incremental" : "Hierarchical")) << " seeding.";
    std::cout << "\n#max_nb_of_proxies = " << max_nb_of_proxies
      << "\n#min_error_drop = " << min_error_drop
      << "\nnb_of_relaxations " << nb_of_relaxations << std::endl;
  }

  approx.initialize_seeds(np);

  if (vl == MAIN_STEPS || vl == VERBOSE)
    std::cout << "Seeding done." << std::endl;

  const std::size_t nb_of_iterations = choose_parameter(
    get_parameter(np, internal_np::number_of_iterations), 20);

  if (vl == VERBOSE)
    std::cout << "\n#nb_of_iterations = " << nb_of_iterations << std::endl;

  approx.run(nb_of_iterations);

  if (vl == MAIN_STEPS || vl == VERBOSE) {
    std::cout << "Approximation done."
      << "\n#proxies = " << approx.number_of_proxies() << std::endl;
  }

  // get proxy map
  approx.proxy_map( get_parameter(np, internal_np::face_proxy_map) );

  if (!parameters::is_default_parameter(get_parameter(np, internal_np::face_proxy_map))
    && (vl == MAIN_STEPS || vl == VERBOSE))
    std::cout << "Filling face proxy map done." << std::endl;

  // get proxies
  approx.proxies( get_parameter(np, internal_np::proxies) );

  if (!is_default_parameter( get_parameter(np, internal_np::proxies) )
    && (vl == MAIN_STEPS || vl == VERBOSE))
    std::cout << "Get proxies done." << std::endl;

  // meshing
  bool is_manifold = false;
  if (!is_default_parameter( get_parameter(np, internal_np::anchors))
    || !is_default_parameter( get_parameter(np, internal_np::triangles) ))
  {
    if (vl == VERBOSE) {
      const FT subdivision_ratio = choose_parameter(get_parameter(np, internal_np::subdivision_ratio), FT(5.0));
      const bool relative_to_chord = choose_parameter(get_parameter(np, internal_np::relative_to_chord), false);
      const bool with_dihedral_angle = choose_parameter(get_parameter(np, internal_np::with_dihedral_angle), false);
      const bool optimize_anchor_location = choose_parameter(get_parameter(np, internal_np::optimize_anchor_location), true);
      const bool pca_plane = choose_parameter(get_parameter(np, internal_np::pca_plane), false);
      std::cout << "Meshing: "
        << "\nchord_error = " << subdivision_ratio
        << "\nrelative_to_chord = " << relative_to_chord
        << "\nwith_dihedral_angle = " << with_dihedral_angle
        << "\noptimize_anchor_location = " << optimize_anchor_location
        << "\npca_plane = " << pca_plane << std::endl;
    }

    is_manifold = approx.extract_mesh(np);

    if (vl == MAIN_STEPS || vl == VERBOSE)
      std::cout << "Meshing done.\n"
        << (is_manifold ? "Can" : "Cannot") << " be built into 2-manifold surface." << std::endl;
  }

  // get anchor points
  approx.anchor_points( get_parameter(np, internal_np::anchors) );

  if (!is_default_parameter( get_parameter(np, internal_np::anchors) )
    && (vl == MAIN_STEPS || vl == VERBOSE))
    std::cout << "Get anchors done." << std::endl;

  // get indexed triangles
  approx.indexed_triangles( get_parameter(np, internal_np::triangles) );

  if (!is_default_parameter( get_parameter(np, internal_np::triangles) )
    && (vl == MAIN_STEPS || vl == VERBOSE))
    std::cout << "Get indexed triangles done." << std::endl;

  return is_manifold;
}

} // namespace VSA
} // end namespace CGAL

#endif // CGAL_SURFACE_MESH_APPROXIMATION_VSA_MESH_APPROXIMATION_H
