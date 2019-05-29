// Copyright (c) 2018 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Konstantinos Katrioplas (konst.katrioplas@gmail.com)
//                 Mael Rouxel-Labbé

#ifndef CGAL_POLYGON_MESH_PROCESSING_SMOOTH_MESH_H
#define CGAL_POLYGON_MESH_PROCESSING_SMOOTH_MESH_H

#include <CGAL/license/Polygon_mesh_processing/meshing_hole_filling.h>

#include <CGAL/Polygon_mesh_processing/internal/Smoothing/mesh_smoothing_impl.h>
#include <CGAL/Polygon_mesh_processing/internal/Smoothing/smoothing_evaluation.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>

#include <CGAL/boost/graph/named_function_params.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <CGAL/property_map.h>

namespace CGAL {
namespace Polygon_mesh_processing {

/*!
* \ingroup PMP_meshing_grp
* smoothes a triangulated region of a polygon mesh using angle-based criteria.
* This function improves the angles of triangle faces by iteratively moving non-constrained vertices.
* Optionally, the points are reprojected to the input surface after each iteration.
*
* @tparam TriangleMesh model of `MutableFaceGraph`.
* @tparam FaceRange range of `boost::graph_traits<TriangleMesh>::%face_descriptor`,
          model of `Range`. Its iterator type is `ForwardIterator`.
* @tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters".
*
* @param tmesh a polygon mesh with triangulated surface patches to be smoothed.
* @param faces the range of triangular faces defining one or several surface patches to be smoothed.
* @param np optional sequence of \ref pmp_namedparameters "Named Parameters" among the ones listed below.
*
* \cgalNamedParamsBegin
*  \cgalParamBegin{geom_traits} a geometric traits class instance, model of `Kernel`.
*    Exact constructions kernels are not supported by this function.
*  \cgalParamEnd
*  \cgalParamBegin{vertex_point_map} the property map with the points associated
*    to the vertices of `tmesh`. Instance of a class model of `ReadWritePropertyMap`.
*  \cgalParamEnd
*  \cgalParamBegin{number_of_iterations} the number of iterations for the
*    sequence of the smoothing iterations performed (default is 1).
*  \cgalParamEnd
*  \cgalParamBegin{vertex_is_constrained_map} a property map containing the
*    constrained-or-not status of each vertex of `tmesh`. A constrained vertex
*    cannot be modified at all during smoothing.
*  \cgalParamEnd
*  \cgalParamBegin{do_project} if `true` (default value), points are projected to the initial surface
*    after each iteration.
*  \cgalParamEnd
* \cgalNamedParamsEnd
*
* @pre `tmesh` does not contain any degenerate faces
*/
template<typename TriangleMesh, typename FaceRange, typename NamedParameters>
void smooth_angles(const FaceRange& faces,
                   TriangleMesh& tmesh,
                   const NamedParameters& np)
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor               vertex_descriptor;

  typedef typename GetGeomTraits<TriangleMesh, NamedParameters>::type                 GeomTraits;
  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters>::type             VertexPointMap;
  typedef typename boost::lookup_named_param_def<internal_np::vertex_is_constrained_t,
                                                 NamedParameters,
                                                 Constant_property_map<vertex_descriptor, bool> // default
                                                 > ::type                             VCMap;

  typedef internal::Angle_smoother<TriangleMesh, VertexPointMap, GeomTraits>          Angle_optimizer;
  typedef internal::Mesh_smoother<Angle_optimizer, TriangleMesh,
                                  VertexPointMap, VCMap, GeomTraits>                  Angle_smoother;
  typedef internal::Delaunay_edge_flipper<TriangleMesh, VertexPointMap, GeomTraits>   Delaunay_flipper;

  if(std::begin(faces) == std::end(faces))
    return;

  using boost::choose_param;
  using boost::get_param;

  // named parameters
  GeomTraits gt = choose_param(get_param(np, internal_np::geom_traits),
                               GeomTraits());
  VertexPointMap vpmap = choose_param(get_param(np, internal_np::vertex_point),
                               get_property_map(CGAL::vertex_point, tmesh));
  VCMap vcmap = choose_param(get_param(np, internal_np::vertex_is_constrained),
                             Constant_property_map<vertex_descriptor, bool>());

  std::size_t nb_iterations = choose_param(get_param(np, internal_np::number_of_iterations), 1);
  const bool do_project = choose_param(get_param(np, internal_np::do_project), true);
  const bool use_safety_constraints = choose_param(get_param(np, internal_np::use_safety_constraints), true);

  Angle_smoother smoother(tmesh, vpmap, vcmap, gt);

  smoother.init_smoothing(faces);

  for(std::size_t i=0; i<nb_iterations; ++i)
  {
    smoother.optimize(use_safety_constraints /*check for bad faces*/,
                      true /*apply all moves at once*/,
                      use_safety_constraints /*check if the min angle is improved*/);

    if(do_project)
    {
      if(use_safety_constraints && does_self_intersect(tmesh))
      {
        std::cerr << "Can't do re-projection, there are self-intersections in the mesh!\n";
        break;
      }

      smoother.project_to_surface();
    }

    // according to the paper, we're supposed to Delaunay-based edge flips!
    Delaunay_flipper delaunay_flipper(tmesh, vpmap, gt);
    delaunay_flipper(faces);
  }
}

template <typename FaceRange, typename TriangleMesh>
void smooth_angles(const FaceRange& face_range, TriangleMesh& tmesh)
{
  smooth_angles(face_range, tmesh, parameters::all_default());
}

template <typename TriangleMesh, typename CGAL_PMP_NP_TEMPLATE_PARAMETERS>
void smooth_angles(TriangleMesh& tmesh, const CGAL_PMP_NP_CLASS& np)
{
  smooth_angles(faces(tmesh), tmesh, np);
}

template<typename TriangleMesh>
void smooth_angles(TriangleMesh& tmesh)
{
  smooth_angles(faces(tmesh), tmesh, parameters::all_default());
}

/*!
* \ingroup PMP_meshing_grp
* smoothes a triangulated region of a polygon mesh using area-based criteria.
* This function tries to make the triangle area distribution as uniform as possible
* by moving non-constrained vertices.
* Optionally, the points are reprojected after each iteration.
*
* @tparam TriangleMesh model of `MutableFaceGraph`.
* @tparam FaceRange range of `boost::graph_traits<TriangleMesh>::%face_descriptor`,
          model of `Range`. Its iterator type is `ForwardIterator`.
* @tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
*
* @param tmesh a polygon mesh with triangulated surface patches to be smoothed.
* @param faces the range of triangular faces defining one or several surface patches to be smoothed.
* @param np optional sequence of \ref pmp_namedparameters "Named Parameters" among the ones listed below.
*
* \cgalNamedParamsBegin
*  \cgalParamBegin{geom_traits} a geometric traits class instance, model of `Kernel`.
*    Exact constructions kernels are not supported by this function.
*  \cgalParamEnd
*  \cgalParamBegin{vertex_point_map} the property map with the points associated
*    to the vertices of `tmesh`. Instance of a class model of `ReadWritePropertyMap`.
*  \cgalParamEnd
*  \cgalParamBegin{number_of_iterations} the number of iterations for the
*    sequence of the smoothing iterations performed (default is 1).
*  \cgalParamEnd
*  \cgalParamBegin{vertex_is_constrained_map} a property map containing the
*    constrained-or-not status of each vertex of `tmesh`. A constrained vertex
*    cannot be modified at all during smoothing.
*  \cgalParamEnd
*  \cgalParamBegin{do_project} if `true` (default value), points are projected to the initial surface after each iteration.
*  \cgalParamEnd
* \cgalNamedParamsEnd
*
* @pre `tmesh` does not contain any degenerate faces
*/
template<typename TriangleMesh, typename FaceRange, typename NamedParameters>
void smooth_areas(const FaceRange faces,
                  TriangleMesh& tmesh,
                  const NamedParameters& np)
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor               vertex_descriptor;

  typedef typename GetGeomTraits<TriangleMesh, NamedParameters>::type                 GeomTraits;
  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters>::type             VertexPointMap;
  typedef typename boost::lookup_named_param_def<internal_np::vertex_is_constrained_t,
                                                 NamedParameters,
                                                 Constant_property_map<vertex_descriptor, bool> // default
                                                 > ::type                             VCMap;

  typedef internal::Area_smoother<TriangleMesh, VertexPointMap, GeomTraits>           Area_optimizer;
  typedef internal::Mesh_smoother<Area_optimizer, TriangleMesh,
                                  VertexPointMap, VCMap, GeomTraits>                  Area_smoother;

  if(std::begin(faces) == std::end(faces))
    return;

  using boost::choose_param;
  using boost::get_param;

  // named parameters
  GeomTraits gt = choose_param(get_param(np, internal_np::geom_traits),
                               GeomTraits());
  VertexPointMap vpmap = choose_param(get_param(np, internal_np::vertex_point),
                               get_property_map(CGAL::vertex_point, tmesh));
  VCMap vcmap = choose_param(get_param(np, internal_np::vertex_is_constrained),
                             Constant_property_map<vertex_descriptor, bool>());

  std::size_t nb_iterations = choose_param(get_param(np, internal_np::number_of_iterations), 1);
  const bool do_project = choose_param(get_param(np, internal_np::do_project), true);
  const bool use_safety_constraints = choose_param(get_param(np, internal_np::use_safety_constraints), true);

  Area_smoother smoother(tmesh, vpmap, vcmap, gt);

  smoother.init_smoothing(faces);

  for(std::size_t i=0; i<nb_iterations; ++i)
  {
    smoother.optimize(use_safety_constraints /*check for bad faces*/,
                      false /*apply moves as soon as they're calculated*/,
                      false /*do not enforce a minimum angle improvement*/);

    if(do_project)
    {
      if(use_safety_constraints && does_self_intersect(tmesh))
      {
        std::cerr << "Can't do re-projection, there are self-intersections in the mesh!\n";
        break;
      }

      smoother.project_to_surface();
    }
  }
}

template <typename FaceRange, typename TriangleMesh>
void smooth_areas(const FaceRange& face_range, TriangleMesh& tmesh)
{
  smooth_areas(face_range, tmesh, parameters::all_default());
}

template <typename TriangleMesh, typename CGAL_PMP_NP_TEMPLATE_PARAMETERS>
void smooth_areas(TriangleMesh& tmesh, const CGAL_PMP_NP_CLASS& np)
{
  smooth_areas(faces(tmesh), tmesh, np);
}

template<typename TriangleMesh>
void smooth_areas(TriangleMesh& tmesh)
{
  smooth_areas(faces(tmesh), tmesh, parameters::all_default());
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///
///
///
///
///
///
///////////////////////////////////////////////////////////////////////////////////////////////////

// do both
template<typename TriangleMesh, typename FaceRange, typename NamedParameters>
void smooth(const FaceRange& faces,
            TriangleMesh& tmesh,
            const NamedParameters& np)
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor               vertex_descriptor;

  typedef typename GetGeomTraits<TriangleMesh, NamedParameters>::type                 GeomTraits;
  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters>::type             VertexPointMap;
  typedef typename boost::lookup_named_param_def<internal_np::vertex_is_constrained_t,
                                                 NamedParameters,
                                                 Constant_property_map<vertex_descriptor, bool> // default
                                                 > ::type                             VCMap;

  typedef internal::Area_smoother<TriangleMesh, VertexPointMap, GeomTraits>           Area_optimizer;
  typedef internal::Mesh_smoother<Area_optimizer, TriangleMesh,
                                  VertexPointMap, VCMap, GeomTraits>                  Area_smoother;
  typedef internal::Angle_smoother<TriangleMesh, VertexPointMap, GeomTraits>          Angle_optimizer;
  typedef internal::Mesh_smoother<Angle_optimizer, TriangleMesh,
                                  VertexPointMap, VCMap, GeomTraits>                  Angle_smoother;
  typedef internal::Delaunay_edge_flipper<TriangleMesh, VertexPointMap, GeomTraits>   Delaunay_flipper;

  if(std::begin(faces) == std::end(faces))
    return;

  using boost::choose_param;
  using boost::get_param;

  // named parameters
  GeomTraits gt = choose_param(get_param(np, internal_np::geom_traits),
                               GeomTraits());
  VertexPointMap vpmap = choose_param(get_param(np, internal_np::vertex_point),
                               get_property_map(CGAL::vertex_point, tmesh));
  VCMap vcmap = choose_param(get_param(np, internal_np::vertex_is_constrained),
                             Constant_property_map<vertex_descriptor, bool>());

  std::size_t nb_iterations = choose_param(get_param(np, internal_np::number_of_iterations), 1);
  const bool do_project = choose_param(get_param(np, internal_np::do_project), true);
  const bool use_safety_constraints = choose_param(get_param(np, internal_np::use_safety_constraints), true);

  Area_smoother area_smoother(tmesh, vpmap, vcmap, gt);
  Angle_smoother angle_smoother(tmesh, vpmap, vcmap, gt);

  area_smoother.init_smoothing(faces);
  angle_smoother.init_smoothing(faces);

  for(std::size_t i=0; i<nb_iterations; ++i)
  {
    // First apply area smoothing...
    area_smoother.optimize(use_safety_constraints /*check for bad faces*/,
                           false /*apply moves as soon as they're calculated*/,
                           false /*do not enforce a minimum angle improvement*/);
    if(do_project)
    {
      if(use_safety_constraints && does_self_intersect(tmesh))
      {
        std::cerr << "Can't do re-projection, there are self-intersections in the mesh!\n";
        break;
      }

      area_smoother.project_to_surface();
    }

    // ... then angle smoothing + Delaunay flip
    angle_smoother.optimize(use_safety_constraints /*check for bad faces*/,
                            true /*apply all moves at once*/,
                            use_safety_constraints /*check if the min angle is improved*/);

    if(do_project)
    {
      if(use_safety_constraints && does_self_intersect(tmesh))
      {
        std::cerr << "Can't do re-projection, there are self-intersections in the mesh!\n";
        break;
      }

      angle_smoother.project_to_surface();
    }

    // according to the paper, we're supposed to Delaunay-based edge flips!
    Delaunay_flipper delaunay_flipper(tmesh, vpmap, gt);
    delaunay_flipper(faces);
  }
}


template <typename FaceRange, typename TriangleMesh>
void smooth(const FaceRange& face_range, TriangleMesh& tmesh)
{
  smooth(face_range, tmesh, parameters::all_default());
}

template <typename TriangleMesh, typename CGAL_PMP_NP_TEMPLATE_PARAMETERS>
void smooth(TriangleMesh& tmesh, const CGAL_PMP_NP_CLASS& np)
{
  smooth(faces(tmesh), tmesh, np);
}

template<typename TriangleMesh>
void smooth(TriangleMesh& tmesh)
{
  smooth(faces(tmesh), tmesh, parameters::all_default());
}


///\cond SKIP_IN_MANUAL
template<typename TriangleMesh, typename GeomTraits, typename Stream>
void angles_evaluation(TriangleMesh& tmesh, GeomTraits traits, Stream& output)
{
  internal::Quality_evaluator<TriangleMesh, GeomTraits> evaluator(tmesh, traits);
  evaluator.gather_angles();
  evaluator.extract_angles(output);
}

template<typename TriangleMesh, typename GeomTraits, typename Stream>
void areas_evaluation(TriangleMesh& tmesh, GeomTraits traits, Stream& output)
{
  internal::Quality_evaluator<TriangleMesh, GeomTraits> evaluator(tmesh, traits);
  evaluator.measure_areas();
  evaluator.extract_areas(output);
}

template<typename TriangleMesh, typename GeomTraits, typename Stream>
void aspect_ratio_evaluation(TriangleMesh& tmesh, GeomTraits traits, Stream& output)
{
  internal::Quality_evaluator<TriangleMesh, GeomTraits> evaluator(tmesh, traits);
  evaluator.calc_aspect_ratios();
  evaluator.extract_aspect_ratios(output);
}
///\cond SKIP_IN_MANUAL

} // namespace Polygon_mesh_processing
} // namespace CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_SMOOTH_MESH_H
