// Copyright (c) 2018-2019 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
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

#ifdef DOXYGEN_RUNNING
#define CGAL_PMP_NP_TEMPLATE_PARAMETERS NamedParameters
#define CGAL_PMP_NP_CLASS NamedParameters
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
* @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
*
* @param tmesh a polygon mesh with triangulated surface patches to be smoothed.
* @param faces the range of triangular faces defining one or several surface patches to be smoothed.
* @param time a time step that corresponds to the speed by which the surface is smoothed.
*        A larger time step results in faster convergence but details may be distorted to have a larger extent
*        compared to more iterations with a smaller step. Typical values scale in the interval (1e-6, 1].
* @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
*
* \cgalNamedParamsBegin
*   \cgalParamNBegin{number_of_iterations}
*     \cgalParamDescription{the number of iterations for the sequence of the smoothing iterations performed}
*     \cgalParamType{unsigned int}
*     \cgalParamDefault{`1`}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{vertex_point_map}
*     \cgalParamDescription{a property map associating points to the vertices of `tmesh`}
*     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
*                    as key type and `%Point_3` as value type}
*     \cgalParamDefault{`boost::get(CGAL::vertex_point, tmesh)`}
*     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
*                     must be available in `TriangleMesh`.}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{geom_traits}
*     \cgalParamDescription{an instance of a geometric traits class}
*     \cgalParamType{a class model of `Kernel`}
*     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
*     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{vertex_is_constrained_map}
*     \cgalParamDescription{a property map containing the constrained-or-not status of each vertex of `tmesh`.}
*     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
*                    as key type and `bool` as value type. It must be default constructible.}
*     \cgalParamDefault{a default property map where no vertex is constrained}
*     \cgalParamExtra{A constrained vertex cannot be modified at all during smoothing.}
*   \cgalParamNEnd
*
*   \cgalParamNBegin{sparse_linear_solver}
*     \cgalParamDescription{an instance of the sparse linear solver used for smoothing}
*     \cgalParamType{a class model of `SparseLinearAlgebraWithFactorTraits_d`}
*     \cgalParamDefault{if \ref thirdpartyEigen "Eigen" 3.2 (or greater) is available and
*                       `CGAL_EIGEN3_ENABLED` is defined, then the following overload of `Eigen_solver_traits`
*                       is provided as default value:
*                       `CGAL::Eigen_solver_traits<Eigen::BiCGSTAB<CGAL::Eigen_sparse_matrix<double>::%EigenType, Eigen::IncompleteLUT<double> > >`}
*   \cgalParamNEnd
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
  typedef typename internal_np::Lookup_named_param_def<
                     internal_np::vertex_is_constrained_t,
                     NamedParameters,
                     Static_boolean_property_map<vertex_descriptor, false> >::type  VCMap;

  using parameters::choose_parameter;
  using parameters::get_parameter;

  GeomTraits gt = choose_parameter<GeomTraits>(get_parameter(np, internal_np::geom_traits));
  VertexPointMap vpmap = choose_parameter(get_parameter(np, internal_np::vertex_point),
                                          get_property_map(CGAL::vertex_point, tmesh));
  VCMap vcmap = choose_parameter(get_parameter(np, internal_np::vertex_is_constrained),
                                 Static_boolean_property_map<vertex_descriptor, false>());
  const unsigned int nb_iterations = choose_parameter(get_parameter(np, internal_np::number_of_iterations), 1);

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

  Sparse_solver solver = choose_parameter<Default_solver>(get_parameter(np, internal_np::sparse_linear_solver));

  std::size_t n = vertices(tmesh).size();
  Eigen_matrix A(n, n);
  Eigen_vector bx(n), by(n), bz(n), Xx(n), Xy(n), Xz(n);
  std::vector<CGAL::Triple<std::size_t, std::size_t, double> > stiffness;

  internal::Shape_smoother<TriangleMesh, VertexPointMap, VCMap, Sparse_solver, GeomTraits> smoother(tmesh, vpmap, vcmap, gt);

  smoother.init_smoothing(faces);

  // For robustness reasons, the laplacian coefficients are computed only once (only the mass
  // matrix is updated at every iteration). See Kazdhan et al. "Can Mean-Curvature Flow Be Made Non-Singular?".
  smoother.calculate_stiffness_matrix_elements(stiffness);

  for(unsigned int iter=0; iter<nb_iterations; ++iter)
  {
#ifdef CGAL_PMP_SMOOTHING_DEBUG
    std::cout << "iteration #" << iter << std::endl;
#endif

    smoother.setup_system(A, bx, by, bz, stiffness, time);

    if(smoother.solve_system(A, Xx, Xy, Xz, bx, by, bz, solver))
    {
      smoother.update_mesh(Xx, Xy, Xz);
    }
    else
    {
#ifdef CGAL_PMP_SMOOTHING_DEBUG
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
