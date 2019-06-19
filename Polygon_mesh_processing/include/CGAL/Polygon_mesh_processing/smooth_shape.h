// Copyright (c) 2018-2019 GeometryFactory (France).
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
// Author(s)     : Mael Rouxel-Labbé
//                 Konstantinos Katrioplas (konst.katrioplas@gmail.com)

#ifndef CGAL_POLYGON_MESH_PROCESSING_SMOOTH_SHAPE_H
#define CGAL_POLYGON_MESH_PROCESSING_SMOOTH_SHAPE_H

#include <CGAL/license/Polygon_mesh_processing/meshing_hole_filling.h>

#if defined(CGAL_EIGEN3_ENABLED)
#include <CGAL/Eigen_solver_traits.h>
#endif

#include <CGAL/Polygon_mesh_processing/internal/Smoothing/curvature_flow_impl.h>

#include <CGAL/utility.h>
#include <CGAL/property_map.h>

#ifdef CGAL_PMP_SMOOTHING_OUTPUT_INTERMEDIARY_STEPS
#include <fstream>
#include <sstream>
#endif

namespace CGAL {
namespace Polygon_mesh_processing {

/*!
* \ingroup PMP_meshing_grp
* smooths the overall shape of the mesh by using the mean curvature flow.
* The effect depends on the curvature of each area and on a time step which
* represents the amount by which vertices are allowed to move.
* The result conformally maps the initial surface to a sphere.
*
* @tparam TriangleMesh model of `MutableFaceGraph`.
* @tparam FaceRange range of `boost::graph_traits<TriangleMesh>::%face_descriptor`,
*         model of `Range`. Its iterator type is `ForwardIterator`.
* @tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
*
* @param tmesh a polygon mesh with triangulated surface patches to be smoothed.
* @param faces the range of triangular faces defining one or several surface patches to be smoothed.
* @param time a time step that corresponds to the speed by which the surface is smoothed.
*        A larger time step results in faster convergence but details may be distorted to have a larger extent
*        compared to more iterations with a smaller step. Typical values scale in the interval (1e-6, 1].
* @param np optional sequence of \ref pmp_namedparameters "Named Parameters" among the ones listed below.
*
* \cgalNamedParamsBegin
*  \cgalParamBegin{vertex_point_map} the property map with the points associated
*    to the vertices of `tmesh`. Instance of a class model of `ReadWritePropertyMap`.
*  \cgalParamEnd
*  \cgalParamBegin{geom_traits} a geometric traits class instance, model of `Kernel`.
*  \cgalParamEnd
*  \cgalParamBegin{vertex_is_constrained_map} a property map containing the
*    constrained-or-not status of each vertex of `tmesh`. A constrained vertex
*    cannot be modified at all during smoothing.
*  \cgalParamEnd
*  \cgalParamBegin{number_of_iterations} the number of iterations for the
*    sequence of the smoothing iterations performed. Each iteration is performed
*    with the given time step.
*  \cgalParamEnd
*  \cgalParamBegin{sparse_linear_solver} an instance of the sparse linear solver used for smoothing \cgalParamEnd
*  \cgalParamEnd
* \cgalNamedParamsEnd
*
*  @warning This function involves linear algebra, that is computed using a non-exact floating-point arithmetic.
*/
template<typename TriangleMesh, typename FaceRange, typename NamedParameters>
void smooth_shape(const FaceRange& faces,
                  TriangleMesh& tmesh,
                  const double time,
                  const NamedParameters& np)
{
  if(std::begin(faces) == std::end(faces))
    return;

  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor      vertex_descriptor;
  typedef typename GetGeomTraits<TriangleMesh, NamedParameters>::type        GeomTraits;
  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters>::type    VertexPointMap;
  typedef typename boost::lookup_named_param_def<
                     internal_np::vertex_is_constrained_t,
                     NamedParameters,
                     Constant_property_map<vertex_descriptor, bool> >::type  VCMap;

  using boost::choose_param;
  using boost::get_param;

  GeomTraits gt = choose_param(get_param(np, internal_np::geom_traits), GeomTraits());
  VertexPointMap vpmap = choose_param(get_param(np, internal_np::vertex_point),
                                      get_property_map(CGAL::vertex_point, tmesh));
  VCMap vcmap = choose_param(get_param(np, internal_np::vertex_is_constrained),
                             Constant_property_map<vertex_descriptor, bool>(false));
  const std::size_t nb_iterations = choose_param(get_param(np, internal_np::number_of_iterations), 1);

#if defined(CGAL_EIGEN3_ENABLED)
#if EIGEN_VERSION_AT_LEAST(3,2,0)
  typedef typename Eigen::SparseMatrix<double>                                         Eigen_sparse_matrix;
  typedef typename Eigen::BiCGSTAB<Eigen_sparse_matrix, Eigen::IncompleteLUT<double> > Eigen_solver;
  typedef CGAL::Eigen_solver_traits<Eigen_solver>                                      Default_solver;
#else
  typedef bool Default_solver;//compilation should crash
  //if no solver is provided and Eigen version < 3.2
#endif
#else
  typedef bool Default_solver;//compilation should crash
  // if no solver is provided and Eigen version < 3.2
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

  typedef typename GetSolver<NamedParameters, Default_solver>::type               Sparse_solver;
  typedef typename Sparse_solver::Matrix                                          Eigen_matrix;
  typedef typename Sparse_solver::Vector                                          Eigen_vector;

  Sparse_solver solver = choose_param(get_param(np, internal_np::sparse_linear_solver), Default_solver());

  std::size_t n = vertices(tmesh).size();
  Eigen_matrix A(n, n);
  Eigen_vector bx(n), by(n), bz(n), Xx(n), Xy(n), Xz(n);
  std::vector<CGAL::Triple<std::size_t, std::size_t, double> > stiffness;

  internal::Shape_smoother<TriangleMesh, VertexPointMap, VCMap, Sparse_solver, GeomTraits> smoother(tmesh, vpmap, vcmap, gt);

  smoother.init_smoothing(faces);

  // For robustness reasons, the laplacian coefficients are computed only once (only the mass
  // matrix is updated at every iteration). See Kazdhan et al. "Can Mean-Curvature Flow Be Made Non-Singular?".
  smoother.calculate_stiffness_matrix_elements(stiffness);

  for(std::size_t iter=0; iter<nb_iterations; ++iter)
  {
#ifdef CGAL_PMP_SMOOTHING_VERBOSE
    std::cout << "iteration #" << iter << std::endl;
#endif

    smoother.setup_system(A, bx, by, bz, stiffness, time);

    if(smoother.solve_system(A, Xx, Xy, Xz, bx, by, bz, solver))
    {
      smoother.update_mesh(Xx, Xy, Xz);
    }
    else
    {
#ifdef CGAL_PMP_SMOOTHING_VERBOSE
      std::cerr << "Failed to solve system!" << std::endl;
#endif
      break;
    }

#ifdef CGAL_PMP_SMOOTHING_OUTPUT_INTERMEDIARY_STEPS
    std::stringstream oss;
    oss << "smoothed_mesh_step_" << iter << ".off" << std::ends;
    std::ofstream out(oss.str().c_str());
    out.precision(17);
    out << tmesh;
    out.close();
#endif
  }
}

template<typename TriangleMesh, typename FaceRange>
void smooth_shape(const FaceRange& faces,
                  TriangleMesh& tmesh,
                  const double time)
{
  smooth_shape(faces, tmesh, time, parameters::all_default());
}

template <typename TriangleMesh, typename CGAL_PMP_NP_TEMPLATE_PARAMETERS>
void smooth_shape(TriangleMesh& tmesh,
                  const double time,
                  const CGAL_PMP_NP_CLASS& np)
{
  smooth_shape(faces(tmesh), tmesh, time, np);
}

template<typename TriangleMesh>
void smooth_shape(TriangleMesh& tmesh,
                  const double time)
{
  smooth_shape(faces(tmesh), tmesh, time, parameters::all_default());
}

} // Polygon_mesh_processing
} // CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_SMOOTH_SHAPE_H
