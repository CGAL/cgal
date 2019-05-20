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

#ifndef CGAL_POLYGON_MESH_PROCESSING_SMOOTH_MESH_H
#define CGAL_POLYGON_MESH_PROCESSING_SMOOTH_MESH_H

#include <CGAL/license/Polygon_mesh_processing/meshing_hole_filling.h>

#include <CGAL/Polygon_mesh_processing/internal/Smoothing/mesh_smoothing_impl.h>
#include <CGAL/Polygon_mesh_processing/internal/Smoothing/smoothing_evaluation.h>
#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>
#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>

#include <CGAL/property_map.h>

#ifdef CGAL_PMP_SMOOTHING_VERBOSE
#include <CGAL/Timer.h>
#endif

namespace CGAL {
namespace Polygon_mesh_processing {

/*!
* \ingroup PMP_meshing_grp
* smoothes a triangulated region of a polygon mesh using angle-based criteria.
* This function improves the angles of triangle faces by iteratively moving non-constrained vertices.
* Optionally, the points are reprojected to the input surface after each iteration.
*
* @tparam PolygonMesh model of `MutableFaceGraph`.
* @tparam FaceRange range of `boost::graph_traits<PolygonMesh>::%face_descriptor`,
          model of `Range`. Its iterator type is `ForwardIterator`.
* @tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters".
*
* @param pmesh a polygon mesh with triangulated surface patches to be smoothed.
* @param faces the range of triangular faces defining one or several surface patches to be smoothed.
* @param np optional sequence of \ref pmp_namedparameters "Named Parameters" among the ones listed below.
*
* \cgalNamedParamsBegin
*  \cgalParamBegin{geom_traits} a geometric traits class instance, model of `Kernel`.
*    Exact constructions kernels are not supported by this function.
*  \cgalParamEnd
*  \cgalParamBegin{vertex_point_map} the property map with the points associated
*    to the vertices of `pmesh`. Instance of a class model of `ReadWritePropertyMap`.
*  \cgalParamEnd
*  \cgalParamBegin{number_of_iterations} the number of iterations for the
*    sequence of the smoothing iterations performed.
*  \cgalParamEnd
*  \cgalParamBegin{vertex_is_constrained_map} a property map containing the
*    constrained-or-not status of each vertex of `pmesh`. A constrained vertex
*    cannot be modified at all during smoothing.
*  \cgalParamEnd
*  \cgalParamBegin{do_project} if `true` (default value), points are projected to the initial surface after each iteration.
*  \cgalParamEnd
* \cgalNamedParamsEnd
*/
template<typename PolygonMesh, typename FaceRange, typename NamedParameters>
void smooth_angles(const FaceRange& faces,
                   PolygonMesh& pmesh,
                   const NamedParameters& np)
{
  if(boost::begin(faces)==boost::end(faces))
    return;

  using boost::choose_param;
  using boost::get_param;

  #ifdef CGAL_PMP_SMOOTHING_VERBOSE
  CGAL::Timer t;
  std::cout << "Smoothing parameters...";
  std::cout.flush();
  t.start();
  #endif

  // geom traits
  typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type GeomTraits;

  // vpmap
  typedef typename GetVertexPointMap<PolygonMesh, NamedParameters>::type   VertexPointMap;
  VertexPointMap vpmap = choose_param(get_param(np, internal_np::vertex_point),
                               get_property_map(CGAL::vertex_point, pmesh));

  // vcmap
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor     vertex_descriptor;
  typedef typename boost::lookup_named_param_def <
      internal_np::vertex_is_constrained_t,
      NamedParameters,
      Constant_property_map<vertex_descriptor, bool> // default
    > ::type                                                               VCMap;
  VCMap vcmap = choose_param(get_param(np, internal_np::vertex_is_constrained),
                             Constant_property_map<vertex_descriptor, bool>());

  // nb_iterations
  std::size_t nb_iterations = choose_param(get_param(np, internal_np::number_of_iterations), 1);
  bool do_project = choose_param(get_param(np, internal_np::do_project), true);

  internal::Compatible_smoother<PolygonMesh, VertexPointMap, VCMap, GeomTraits>
          smoother(pmesh, vpmap, vcmap);

#ifdef CGAL_PMP_SMOOTHING_VERBOSE
  t.stop();
  std::cout << " done ("<< t.time() <<" sec)." << std::endl;
  std::cout << "Removing degenerate faces..." << std::endl;
  t.reset(); t.start();
#endif

  smoother.remove_degenerate_faces();

#ifdef CGAL_PMP_SMOOTHING_VERBOSE
  t.stop();
  std::cout << " done ("<< t.time() <<" sec)." << std::endl;
  std::cout << "Initializing..." << std::endl;
  t.reset(); t.start();
#endif

  smoother.init_smoothing(faces);

#ifdef CGAL_PMP_SMOOTHING_VERBOSE
  t.stop();
  std::cout << " done ("<< t.time() <<" sec)." << std::endl;
  std::cout << "#iter = " << nb_iterations << std::endl;
  std::cout << "Smoothing ..." << std::endl;
  t.reset(); t.start();
#endif

  for(std::size_t i=0; i<nb_iterations; ++i)
  {
#ifdef CGAL_PMP_SMOOTHING_VERBOSE
    std::cout << " * Iteration " << (i + 1) << " *" << std::endl;
#endif

    smoother.angle_relaxation();
    if(do_project)
    {
      if(does_self_intersect(pmesh))
      {
        std::cerr << "Can't do re-projection, there are self-intersections in the mesh!\n";
        break;
      }
      smoother.project_to_surface();
    }
  }

#ifdef CGAL_PMP_SMOOTHING_VERBOSE
  t.stop();
  std::cout << "Smoothing done in ";
  std::cout << t.time() << " sec." << std::endl;
#endif
}

template<typename PolygonMesh, typename NamedParameters>
void smooth_angles(PolygonMesh& pmesh, const NamedParameters& np)
{
  smooth_angles(faces(pmesh), pmesh, np);
}

template<typename PolygonMesh>
void smooth_angles(PolygonMesh& pmesh)
{
  smooth_angles(faces(pmesh), pmesh, parameters::all_default());
}

/*!
* \ingroup PMP_meshing_grp
* smoothes a triangulated region of a polygon mesh using area based criteria.
* This function tries to make the triangle area distribution as uniform as possible
* by moving non constrained vertices.
* Optionally, the points are reprojected after each iteration.
*
* @tparam PolygonMesh model of `MutableFaceGraph`.
* @tparam FaceRange range of `boost::graph_traits<PolygonMesh>::%face_descriptor`,
          model of `Range`. Its iterator type is `ForwardIterator`.
* @tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
*
* @param pmesh a polygon mesh with triangulated surface patches to be smoothed.
* @param faces the range of triangular faces defining one or several surface patches to be smoothed.
* @param np optional sequence of \ref pmp_namedparameters "Named Parameters" among the ones listed below.
*
* \cgalNamedParamsBegin
*  \cgalParamBegin{geom_traits} a geometric traits class instance, model of `Kernel`.
*    Exact constructions kernels are not supported by this function.
*  \cgalParamEnd
*  \cgalParamBegin{vertex_point_map} the property map with the points associated
*    to the vertices of `pmesh`. Instance of a class model of `ReadWritePropertyMap`.
*  \cgalParamEnd
*  \cgalParamBegin{number_of_iterations} the number of iterations for the
*    sequence of the smoothing iterations performed.
*  \cgalParamEnd
*  \cgalParamBegin{vertex_is_constrained_map} a property map containing the
*    constrained-or-not status of each vertex of `pmesh`. A constrained vertex
*    cannot be modified at all during smoothing.
*  \cgalParamEnd
*  \cgalParamBegin{gradient_descent_precision} precision which is achieved by minimizing
*    the energy of each triangle. The energy is defined based on the triangle area.
*    A small value corresponds to high precision.
*  \cgalParamEnd
*  \cgalParamBegin{do_project} if `true` (default value), points are projected to the initial surface after each iteration.
*  \cgalParamEnd
* \cgalNamedParamsEnd
*/
template<typename PolygonMesh, typename FaceRange, typename NamedParameters>
void smooth_areas(const FaceRange& faces, PolygonMesh& pmesh, const NamedParameters& np)
{
  if (boost::begin(faces)==boost::end(faces))
    return;

  using boost::choose_param;
  using boost::get_param;

  #ifdef CGAL_PMP_SMOOTHING_VERBOSE
  CGAL::Timer t;
  std::cout << "Smoothing parameters...";
  std::cout.flush();
  t.start();
  #endif

  typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type GeomTraits;

  typedef typename GetVertexPointMap<PolygonMesh, NamedParameters>::type VertexPointMap;
  VertexPointMap vpmap = choose_param(get_param(np, internal_np::vertex_point),
                               get_property_map(CGAL::vertex_point, pmesh));

  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;
  typedef typename boost::lookup_named_param_def <
      internal_np::vertex_is_constrained_t,
      NamedParameters,
      Constant_property_map<vertex_descriptor, bool>
    > ::type VCMap;
  VCMap vcmap = choose_param(get_param(np, internal_np::vertex_is_constrained),
                             Constant_property_map<vertex_descriptor, bool>());

  std::size_t nb_iterations = choose_param(get_param(np, internal_np::number_of_iterations), 1);
  const double gd_precision = choose_param(get_param(np, internal_np::gradient_descent_precision), 1e-5);
  bool do_project = choose_param(get_param(np, internal_np::do_project), true);

  internal::Compatible_smoother<PolygonMesh, VertexPointMap, VCMap, GeomTraits>
          smoother(pmesh, vpmap, vcmap);

#ifdef CGAL_PMP_SMOOTHING_VERBOSE
  t.stop();
  std::cout << " done ("<< t.time() <<" sec)." << std::endl;
  std::cout << "Removing degenerate faces..." << std::endl;
  t.reset(); t.start();
#endif

  smoother.remove_degenerate_faces();

#ifdef CGAL_PMP_SMOOTHING_VERBOSE
  t.stop();
  std::cout << " done ("<< t.time() <<" sec)." << std::endl;
  std::cout << "Initializing..." << std::endl;
  t.reset(); t.start();
#endif

  smoother.init_smoothing(faces);

#ifdef CGAL_PMP_SMOOTHING_VERBOSE
  t.stop();
  std::cout << " done ("<< t.time() <<" sec)." << std::endl;
  std::cout << "#iter = " << nb_iterations << std::endl;
  std::cout << "Smoothing ..." << std::endl;
  t.reset(); t.start();
#endif

  for(std::size_t i=0; i<nb_iterations; ++i)
  {
#ifdef CGAL_PMP_SMOOTHING_VERBOSE
    std::cout << " * Iteration " << (i + 1) << " *" << std::endl;
#endif

    smoother.area_relaxation(gd_precision);
    if(do_project)
    {
      if(does_self_intersect(pmesh))
      {
        std::cerr << "Can't do re-projection, there are self-intersections in the mesh!\n";
        break;
      }
      smoother.project_to_surface();
    }
  }

#ifdef CGAL_PMP_SMOOTHING_VERBOSE
  t.stop();
  std::cout << "Smoothing done in ";
  std::cout << t.time() << " sec." << std::endl;
#endif
}

template<typename PolygonMesh, typename NamedParameters>
void smooth_areas(PolygonMesh& pmesh, const NamedParameters& np)
{
  smooth_areas(faces(pmesh), pmesh, np);
}

template<typename PolygonMesh>
void smooth_areas(PolygonMesh& pmesh)
{
  smooth_areas(faces(pmesh), pmesh, parameters::all_default());
}


///\cond SKIP_IN_MANUAL
template<typename PolygonMesh, typename GeomTraits, typename Stream>
void angles_evaluation(PolygonMesh& pmesh, GeomTraits, Stream& output)
{
  internal::Quality_evaluator<PolygonMesh, GeomTraits> evaluator(pmesh);
  evaluator.gather_angles();
  evaluator.extract_angles(output);
}

template<typename PolygonMesh, typename GeomTraits, typename Stream>
void areas_evaluation(PolygonMesh& pmesh, GeomTraits, Stream& output)
{
  internal::Quality_evaluator<PolygonMesh, GeomTraits> evaluator(pmesh);
  evaluator.measure_areas();
  evaluator.extract_areas(output);
}

template<typename PolygonMesh, typename GeomTraits, typename Stream>
void aspect_ratio_evaluation(PolygonMesh& pmesh, GeomTraits, Stream& output)
{
  internal::Quality_evaluator<PolygonMesh, GeomTraits> evaluator(pmesh);
  evaluator.calc_aspect_ratios();
  evaluator.extract_aspect_ratios(output);
}
///\cond SKIP_IN_MANUAL

} // namespace Polygon_mesh_processing
} // namespace CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_SMOOTH_MESH_H
