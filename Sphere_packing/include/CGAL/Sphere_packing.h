// Copyright (c) 2025 Universitaet Bremen
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s) : Rene Weller, Sven Oesau
//

#ifndef CGAL_SPHERE_PACKINGH_H_
#define CGAL_SPHERE_PACKINGH_H_

#include <CGAL/license/Sphere_packing.h>

#include <CGAL/disable_warnings.h>
#include <CGAL/Named_function_parameters.h>


namespace CGAL {

/**
 * \ingroup PkgSpherePackingAlgorithms
 * \brief Packs spheres into a closed mesh.
 * Computes a packing of spheres into a closed and self-intersection free triangle mesh. The possible range of radii of spheres can be chosen; large spheres are preferred during the packing.
 * The optimal solution is not guaranteed as it has NP-hard complexity.
 *
 * Note that the precision of the method is limited by the precision of the GPU, i.e., non-exact arithmetic with 32-bit floating point numbers is used.
 *
 * \tparam FaceGraph a model of `HalfedgeListGraph`, `FaceListGraph`, and `MutableFaceGraph`
 *
 * \tparam OutputIterator must be an output iterator accepting variables of type `Sphere_3`.
 *
 * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * \param mesh the input mesh to pack spheres into. It has to be a closed triangle mesh and may not have self-intersections.
 * \param out output iterator into which the packed spheres are written.
 * \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *
 *   \cgalParamNBegin{number_of_target_spheres}
 *     \cgalParamDescription{target number of spheres to be packed into `mesh`.}
 *     \cgalParamType{unsigned int}
 *     \cgalParamDefault{10,000}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{minimum_radius}
 *     \cgalParamDescription{minimum radius of packed spheres}
 *     \cgalParamType{float}
 *     \cgalParamDefault{0}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{maximum_radius}
 *     \cgalParamDescription{maximum radius of packed spheres}
 *     \cgalParamType{float}
 *     \cgalParamDefault{FLT_MAX}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{iteration_functor}
 *     \cgalParamDescription{a functor that is called before each split to decide whether the packing should continue.}
 *     \cgalParamType{an instance of `std::function<bool(unsigned int iteration, unsigned int performed_splits, unsigned int number_of_spheres, float sphere_volume, float object_volume, const float4 *)>`.}
 *     \cgalParamDefault{No functor is used. Packing continues until `number_of_target_spheres` or `maximum_splits` are met.}
 *     \cgalParamExtra{`float4` is a CUDA data type and has 4 members: `float x, y, z, w`.}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{initial_grid_resolution}
 *     \cgalParamDescription{number of grid cells for the longest side of the bounding box. Grid cells are cubes.}
 *     \cgalParamType{unsigned int}
 *     \cgalParamDefault{4}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{maximum_splits}
 *     \cgalParamDescription{how often the grid is subdivided at most.}
 *     \cgalParamType{unsigned int}
 *     \cgalParamDefault{6}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{iterations_between_splits}
 *     \cgalParamDescription{iterations of placing spheres between splitting the grid into a higher resolution.}
 *     \cgalParamType{unsigned int}
 *     \cgalParamDefault{30}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{vertex_point_map}
 *     \cgalParamDescription{a property map associating points to the vertices of `mesh`}
 *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<FaceGraph>::%vertex_descriptor`
 *                    as key type and `%Point_3` as value type}
 *     \cgalParamDefault{`boost::get(CGAL::vertex_point, mesh)`}
 *     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
 *                     must be available in `FaceGraph`.}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{geom_traits}
 *     \cgalParamDescription{an instance of a geometric traits class}
 *     \cgalParamType{a class model of `Kernel`}
 *     \cgalParamDefault{a \cgal Kernel deduced from the point type of `FaceGraph`, using `CGAL::Kernel_traits`}
 *     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
 *   \cgalParamNEnd
 *
 * \cgalNamedParamsEnd
 *
 * \return number of convex hulls. Can be lower than the `maximum_number_of_convex_hulls` if the specified `volume_error` is quickly met.
 */


template<typename FaceGraph, typename OutputIterator, typename NamedParameters = parameters::Default_named_parameters>
std::size_t pack_spheres(const FaceGraph& mesh, OutputIterator &out, const NamedParameters& np = parameters::default_values()) {

}

}  // end namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_SPHERE_PACKINGH_H_
