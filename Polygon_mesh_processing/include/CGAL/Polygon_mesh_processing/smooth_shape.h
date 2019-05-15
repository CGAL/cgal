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

#ifndef CGAL_POLYGON_MESH_PROCESSING_SMOOTH_SHAPE_H
#define CGAL_POLYGON_MESH_PROCESSING_SMOOTH_SHAPE_H

#include <CGAL/license/Polygon_mesh_processing/meshing_hole_filling.h>

#if defined(CGAL_EIGEN3_ENABLED)
#include <CGAL/Eigen_solver_traits.h>
#endif

#include <CGAL/Polygon_mesh_processing/internal/Smoothing/curvature_flow_impl.h>
#include <CGAL/Polygon_mesh_processing/internal/Smoothing/curvature_flow_explicit_impl.h>
#include <CGAL/utility.h>
#include <CGAL/property_map.h>

#ifdef CGAL_PMP_SMOOTHING_VERBOSE
#include <CGAL/Timer.h>
#endif

namespace CGAL {
namespace Polygon_mesh_processing {

// normalized explicit scheme, undocumented
template<typename PolygonMesh, typename FaceRange, typename NamedParameters>
void smooth_curvature_flow_explicit(const FaceRange& faces, PolygonMesh& pmesh, const NamedParameters& np)
{
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
                             Constant_property_map<vertex_descriptor, bool>(false));

  std::size_t nb_iterations = choose_param(get_param(np, internal_np::number_of_iterations), 1);

  #ifdef CGAL_PMP_SMOOTHING_VERBOSE
  t.stop();
  std::cout << "\rSmoothing parameters done ("<< t.time() <<" sec)" << std::endl;
  std::cout << "smoother construction...";
  std::cout.flush();
  t.reset(); t.start();
  #endif

  internal::Curvature_flow<PolygonMesh, VertexPointMap, VCMap, GeomTraits>
          curvature_smoother(pmesh, vpmap, vcmap);

  #ifdef CGAL_PMP_SMOOTHING_VERBOSE
  t.stop();
  std::cout << " done ("<< t.time() <<" sec)." << std::endl;
  std::cout << "Removing degenerate faces..." << std::endl;
  t.reset(); t.start();
  #endif

  curvature_smoother.remove_degenerate_faces();

  #ifdef CGAL_PMP_SMOOTHING_VERBOSE
  t.stop();
  std::cout << " done ("<< t.time() <<" sec)." << std::endl;
  std::cout << "Initializing..." << std::endl;
  t.reset(); t.start();
  #endif

  curvature_smoother.init_smoothing(faces);

  #ifdef CGAL_PMP_SMOOTHING_VERBOSE
  t.stop();
  std::cout << " done ("<< t.time() <<" sec)." << std::endl;
  std::cout << "#iter = " << nb_iterations << std::endl;
  std::cout << "Shape smoothing..." << std::endl;
  t.reset(); t.start();
  #endif

  for(std::size_t i=0; i<nb_iterations; ++i)
  {

    #ifdef CGAL_PMP_SMOOTHING_VERBOSE
        std::cout << " * Iteration " << (i + 1) << " *" << std::endl;
    #endif

    curvature_smoother.curvature_smoothing();

  }

  #ifdef CGAL_PMP_SMOOTHING_VERBOSE
  t.stop();
  std::cout << "Shape smoothing done in ";
  std::cout << t.time() << " sec." << std::endl;
  std::cout<<std::endl;
  #endif
}

/*!
* \ingroup PMP_meshing_grp
* smooths the overall shape of the mesh by using the mean curvature flow.
* The effect depends on the curvature of each area and on a time step which
* represents the amount by which vertices are allowed to move.
* The result conformally maps the initial surface to a sphere.
*
*
* @tparam PolygonMesh model of `MutableFaceGraph`.
* @tparam FaceRange range of `boost::graph_traits<PolygonMesh>::%face_descriptor`,
*         model of `Range`. Its iterator type is `ForwardIterator`.
* @tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
*
* @param pmesh a polygon mesh with triangulated surface patches to be smoothed.
* @param faces the range of triangular faces defining one or several surface patches to be smoothed.
* @param time a time step that corresponds to the speed by which the surface is smoothed.
*        A larger time step results in faster convergence but details may be distorted to have a larger extent
*        compared to more iterations with a smaller step. Typical values scale in the interval (1e-6, 1].
* @param np optional sequence of \ref pmp_namedparameters "Named Parameters" among the ones listed below.
*
* \cgalNamedParamsBegin
*  \cgalParamBegin{geom_traits} a geometric traits class instance, model of `Kernel`.
*    Kernels with exact constructions are not supported by this function.
*  \cgalParamEnd
*  \cgalParamBegin{vertex_point_map} the property map with the points associated
*    to the vertices of `pmesh`. Instance of a class model of `ReadWritePropertyMap`.
*  \cgalParamEnd
*  \cgalParamBegin{vertex_is_constrained_map} a property map containing the
*    constrained-or-not status of each vertex of `pmesh`. A constrained vertex
*    cannot be modified at all during smoothing.
*  \cgalParamEnd
*  \cgalParamBegin{number_of_iterations} the number of iterations for the
*    sequence of the smoothing iterations performed. Each iteration is performed
*    with the given time step.
*  \cgalParamEnd
*  \cgalParamBegin{use_explicit_scheme} a normalized explicit scheme for smoothing along the curvature flow.
*    It is not dependent on the time step parameter. Each vertex is moved sequentially without solving
*    a linear system. It can be useful for subtle noise removal and to avoid the result
*    to converge to a sphere.
* \cgalParamBegin{sparse_linear_solver} an instance of the sparse linear solver used for smoothing \cgalParamEnd
*  \cgalParamEnd
* \cgalNamedParamsEnd
*
*  @warning This function involves linear algebra, that is computed using a non-exact floating-point arithmetic.
*/
template<typename PolygonMesh, typename FaceRange, typename NamedParameters>
void smooth_along_curvature_flow(const FaceRange& faces, PolygonMesh& pmesh, const double time,
                                 const NamedParameters& np)
{
  if (boost::begin(faces)==boost::end(faces))
    return;

  using boost::choose_param;
  using boost::get_param;

  #ifdef CGAL_PMP_SMOOTHING_VERBOSE
  CGAL::Timer t;
  std::cout << "Initializing...";
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
                             Constant_property_map<vertex_descriptor, bool>(false));

  std::size_t nb_iterations = choose_param(get_param(np, internal_np::number_of_iterations), 1);

  bool use_explicit_scheme = choose_param(get_param(np, internal_np::use_explicit_scheme), false);

  if(use_explicit_scheme)
  {
    smooth_curvature_flow_explicit(faces, pmesh, parameters::number_of_iterations(nb_iterations));
  }
  else
  {
    // implicit scheme
#if defined(CGAL_EIGEN3_ENABLED)
  #if EIGEN_VERSION_AT_LEAST(3,2,0)
    typedef typename Eigen::SparseMatrix<double> Eigen_sparse_matrix;
    typedef typename Eigen::BiCGSTAB<Eigen_sparse_matrix, Eigen::IncompleteLUT<double> > Eigen_solver;
    typedef CGAL::Eigen_solver_traits<Eigen_solver> Default_solver;
  #else
    typedef bool Default_solver;//compilation should crash
      //if no solver is provided and Eigen version < 3.2
  #endif
#else
    typedef bool Default_solver;//compilation should crash
      //if no solver is provided and Eigen version < 3.2
#endif

#if defined(CGAL_EIGEN3_ENABLED)
    CGAL_static_assertion_msg(
      (!boost::is_same<typename GetSolver<NamedParameters, Default_solver>::type, bool>::value) || EIGEN_VERSION_AT_LEAST(3, 2, 0),
      "Eigen3 version 3.2 or later is required.");
#else
    CGAL_static_assertion_msg(
      (!boost::is_same<typename GetSolver<NamedParameters, Default_solver>::type, bool>::value),
      "Eigen3 version 3.2 or later is required.");
#endif

    typedef typename GetSolver<NamedParameters, Default_solver>::type Sparse_solver;
    typedef typename Sparse_solver::Matrix Eigen_matrix;
    typedef typename Sparse_solver::Vector Eigen_vector;
    Sparse_solver solver = choose_param(get_param(np, internal_np::sparse_linear_solver), Default_solver());

    std::size_t n = static_cast<int>(vertices(pmesh).size());
    Eigen_matrix A(n, n);
    Eigen_vector bx(n), by(n), bz(n), Xx(n), Xy(n), Xz(n);
    std::vector<CGAL::Triple<int, int, double> > stiffness;

    internal::Shape_smoother<PolygonMesh, VertexPointMap, VCMap, Sparse_solver, GeomTraits>
        smoother(pmesh, vpmap, vcmap);

    smoother.init_smoothing(faces);

    #ifdef CGAL_PMP_SMOOTHING_VERBOSE
    t.stop();
    std::cout << " done ("<< t.time() <<" sec)." << std::endl;
    std::cout << "calculate_stiffness_matrix_elements... ";
    t.reset(); t.start();
    #endif

    smoother.calculate_stiffness_matrix_elements(stiffness);

    #ifdef CGAL_PMP_SMOOTHING_VERBOSE
    t.stop();
    std::cout << " done ("<< t.time() <<" sec)." << std::endl;
    #endif

    for(std::size_t iter = 0; iter < nb_iterations; ++iter)
    {

      #ifdef CGAL_PMP_SMOOTHING_VERBOSE
      std::cout << "setup_system... ";
      t.reset(); t.start();
      #endif

      smoother.setup_system(A, bx, by, bz, stiffness, time);

      #ifdef CGAL_PMP_SMOOTHING_VERBOSE
      t.stop();
      std::cout << " done ("<< t.time() <<" sec)." << std::endl;
      std::cout << "solve_system... ";
      t.reset(); t.start();
      #endif

      smoother.solve_system(A, Xx, Xy, Xz, bx, by, bz, solver);

      #ifdef CGAL_PMP_SMOOTHING_VERBOSE
      t.stop();
      std::cout << " done ("<< t.time() <<" sec)." << std::endl;
      std::cout << "update_mesh... ";
      t.reset(); t.start();
      #endif

      smoother.update_mesh(Xx, Xy, Xz);

      #ifdef CGAL_PMP_SMOOTHING_VERBOSE
      t.stop();
      std::cout << " done ("<< t.time() <<" sec)." << std::endl;
      #endif
    }
  }
}

template<typename PolygonMesh, typename NamedParameters>
void smooth_along_curvature_flow(PolygonMesh& pmesh, const double time,
                                 const NamedParameters& np)
{
  smooth_along_curvature_flow(faces(pmesh), pmesh, time, np);
}

template<typename PolygonMesh>
void smooth_along_curvature_flow(PolygonMesh& pmesh, const double time)
{
  smooth_along_curvature_flow(faces(pmesh), pmesh, time,
                              parameters::all_default());
}

// API for Polyhedron demo plugin
namespace internal{

template<typename PolygonMesh, typename FaceRange, typename NamedParameters>
void solve_mcf(const FaceRange& faces, PolygonMesh& mesh, const double time,
               std::vector<CGAL::Triple<int, int, double> >& stiffness, bool compute_stiffness,
               const NamedParameters& np)
{
  using boost::choose_param;
  using boost::get_param;

  typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type GeomTraits;

  typedef typename GetVertexPointMap<PolygonMesh, NamedParameters>::type VertexPointMap;
  VertexPointMap vpmap = choose_param(get_param(np, internal_np::vertex_point),
                               get_property_map(CGAL::vertex_point, mesh));

  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;
  typedef typename boost::lookup_named_param_def <
      internal_np::vertex_is_constrained_t,
      NamedParameters,
      Constant_property_map<vertex_descriptor, bool>
    > ::type VCMap;
  VCMap vcmap = choose_param(get_param(np, internal_np::vertex_is_constrained),
                             Constant_property_map<vertex_descriptor, bool>(false));

  // nb_iterations
  std::size_t nb_iterations = choose_param(get_param(np, internal_np::number_of_iterations), 1);

#if defined(CGAL_EIGEN3_ENABLED)
  #if EIGEN_VERSION_AT_LEAST(3,2,0)
    typedef typename Eigen::SparseMatrix<double> Eigen_sparse_matrix;
    typedef typename Eigen::BiCGSTAB<Eigen_sparse_matrix, Eigen::IncompleteLUT<double> > Eigen_solver;
    typedef CGAL::Eigen_solver_traits<Eigen_solver> Default_solver;
  #else
    typedef bool Default_solver;//compilation should crash
      //if no solver is provided and Eigen version < 3.2
  #endif
#else
    typedef bool Default_solver;//compilation should crash
      //if no solver is provided and Eigen version < 3.2
#endif

#if defined(CGAL_EIGEN3_ENABLED)
    CGAL_static_assertion_msg(
      (!boost::is_same<typename GetSolver<NamedParameters, Default_solver>::type, bool>::value) || EIGEN_VERSION_AT_LEAST(3, 2, 0),
      "Eigen3 version 3.2 or later is required.");
#else
    CGAL_static_assertion_msg(
      (!boost::is_same<typename GetSolver<NamedParameters, Default_solver>::type, bool>::value),
      "Eigen3 version 3.2 or later is required.");
#endif

  typedef typename GetSolver<NamedParameters, Default_solver>::type Sparse_solver;
  typedef typename Sparse_solver::Matrix Eigen_matrix;
  typedef typename Sparse_solver::Vector Eigen_vector;
  Sparse_solver solver = choose_param(get_param(np, internal_np::sparse_linear_solver), Default_solver());

  std::size_t n = static_cast<int>(vertices(mesh).size());

  Eigen_matrix A(n, n);
  Eigen_vector bx(n), by(n), bz(n), Xx(n), Xy(n), Xz(n);

  if(compute_stiffness)
  {
    internal::Shape_smoother<PolygonMesh, VertexPointMap, VCMap, Sparse_solver, GeomTraits>
        smoother(mesh, vpmap, vcmap);
    smoother.init_smoothing(faces);
    smoother.calculate_stiffness_matrix_elements(stiffness);
  }
  else
  {
    internal::Shape_smoother<PolygonMesh, VertexPointMap, VCMap, Default_solver, GeomTraits>
        smoother(mesh, vpmap, vcmap);
    smoother.init_smoothing(faces);

    for(std::size_t i=0; i<nb_iterations; ++i)
    {
      smoother.setup_system(A, bx, by, bz, stiffness, time);
      smoother.solve_system(A, Xx, Xy, Xz, bx, by, bz, solver);
      smoother.update_mesh(Xx, Xy, Xz);
    }
  }
}

} // end internal namespace


} //Polygon_mesh_processing
} //CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_SMOOTH_SHAPE_H
