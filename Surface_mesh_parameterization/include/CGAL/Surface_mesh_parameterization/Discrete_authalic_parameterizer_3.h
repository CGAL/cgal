// Copyright (c) 2005  INRIA (France).
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
// Author(s)     : Laurent Saboret, Pierre Alliez, Bruno Levy

#ifndef CGAL_SURFACE_MESH_PARAMETERIZATION_DISCRETE_AUTHALIC_PARAMETERIZER_3_H
#define CGAL_SURFACE_MESH_PARAMETERIZATION_DISCRETE_AUTHALIC_PARAMETERIZER_3_H

#include <CGAL/license/Surface_mesh_parameterization.h>

#include <CGAL/Surface_mesh_parameterization/internal/angles.h>
#include <CGAL/Surface_mesh_parameterization/internal/kernel_traits.h>
#include <CGAL/Surface_mesh_parameterization/Error_code.h>

#include <CGAL/Surface_mesh_parameterization/Fixed_border_parameterizer_3.h>

#include <CGAL/Default.h>

#if defined(CGAL_EIGEN3_ENABLED)
#include <CGAL/Eigen_solver_traits.h>
#endif

/// \file Discrete_authalic_parameterizer_3.h

namespace CGAL {

namespace Surface_mesh_parameterization {

/// \ingroup  PkgSurfaceParameterizationMethods
///
/// The class `Discrete_authalic_parameterizer_3`
/// implements the *Discrete Authalic Parameterization* algorithm. This method
/// is sometimes called <i>DAP</i> or just <i>Authalic parameterization</i> \cgalCite{cgal:dma-ipsm-02}.
///
/// DAP is a weak area-preserving parameterization. It is a compromise between
/// area-preserving and angle-preserving.
///
/// A one-to-one mapping is guaranteed if the surface's border is mapped onto a convex polygon.
///
/// This class is a strategy  called by the main
/// parameterization algorithm `Fixed_border_parameterizer_3::parameterize()` and it:
/// - provides the template parameters `BorderParameterizer_` and `SolverTraits_`.
/// - implements `compute_w_ij()` to compute w_ij = (i, j), coefficient of the matrix A
///   for j neighbor vertex of i, based on Discrete Authalic Parameterization algorithm.
///
/// \cgalModels `Parameterizer_3`
///
/// \tparam TriangleMesh_ must be a model of `FaceGraph`.
///
/// \tparam BorderParameterizer_ is a Strategy to parameterize the surface border
///         and must be a model of `Parameterizer_3`.<br>
///         <b>%Default:</b>
/// \code
///   Circular_border_arc_length_parameterizer_3<TriangleMesh_>
/// \endcode
///
/// \tparam SolverTraits_ must be a model of `SparseLinearAlgebraTraits_d`.<br>
///         Note that the system is *not* symmetric because `Fixed_border_parameterizer_3`
///         does not remove border vertices from the system.<br>
///         <b>%Default:</b> If \ref thirdpartyEigen "Eigen" 3.1 (or greater) is available
///         and `CGAL_EIGEN3_ENABLED` is defined, then an overload of `Eigen_solver_traits`
///         is provided as default parameter:
/// \code
///   CGAL::Eigen_solver_traits<
///           Eigen::BiCGSTAB<Eigen_sparse_matrix<double>::EigenType,
///                           Eigen::IncompleteLUT< double > > >
/// \endcode
///
/// \sa `CGAL::Surface_mesh_parameterization::Fixed_border_parameterizer_3<TriangleMesh, BorderParameterizer, SolverTraits>`
/// \sa `CGAL::Surface_mesh_parameterization::ARAP_parameterizer_3<TriangleMesh, BorderParameterizer, SolverTraits>`
/// \sa `CGAL::Surface_mesh_parameterization::Barycentric_mapping_parameterizer_3<TriangleMesh, BorderParameterizer, SolverTraits>`
/// \sa `CGAL::Surface_mesh_parameterization::Discrete_conformal_map_parameterizer_3<TriangleMesh, BorderParameterizer, SolverTraits>`
/// \sa `CGAL::Surface_mesh_parameterization::LSCM_parameterizer_3<TriangleMesh, BorderParameterizer>`
/// \sa `CGAL::Surface_mesh_parameterization::Mean_value_coordinates_parameterizer_3<TriangleMesh, BorderParameterizer, SolverTraits>`
/// \sa `CGAL::Surface_mesh_parameterization::Orbifold_Tutte_parameterizer_3<SeamMesh, SolverTraits>`
///
template < class TriangleMesh_,
           class BorderParameterizer_ = Default,
           class SolverTraits_ = Default>
class Discrete_authalic_parameterizer_3
  : public Fixed_border_parameterizer_3<
      TriangleMesh_,
      typename Default::Get<
        BorderParameterizer_,
        Circular_border_arc_length_parameterizer_3<TriangleMesh_> >::type,
      typename Default::Get<
        SolverTraits_,
#if defined(CGAL_EIGEN3_ENABLED)
        CGAL::Eigen_solver_traits<
          Eigen::BiCGSTAB<Eigen_sparse_matrix<double>::EigenType,
                          Eigen::IncompleteLUT<double> > > >::type >
#else
  #pragma message("Error: You must either provide 'SolverTraits_' or link CGAL with the Eigen library")
       SolverTraits_>::type > // no parameter provided, and Eigen is not enabled: don't compile
#endif
{
public:
#ifndef DOXYGEN_RUNNING
  typedef typename Default::Get<
    BorderParameterizer_,
    Circular_border_arc_length_parameterizer_3<TriangleMesh_> >::type  Border_parameterizer;

  typedef typename Default::Get<
    SolverTraits_,
  #if defined(CGAL_EIGEN3_ENABLED)
    CGAL::Eigen_solver_traits<
      Eigen::BiCGSTAB<Eigen_sparse_matrix<double>::EigenType,
                      Eigen::IncompleteLUT<double> > >
  #else
    SolverTraits_ // no parameter provided, and Eigen is not enabled: so don't compile!
  #endif
  >::type                                                     Solver_traits;
#else
  typedef Border_parameterizer_                               Border_parameterizer;
  typedef SolverTraits_                                       Solver_traits;
#endif

  typedef TriangleMesh_                                       TriangleMesh;

// Private types
private:
  // Superclass
  typedef Fixed_border_parameterizer_3<TriangleMesh,
                                       Border_parameterizer,
                                       Solver_traits>         Base;

// Private types
private:
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor vertex_descriptor;
  typedef CGAL::Vertex_around_target_circulator<TriangleMesh> vertex_around_target_circulator;

  // Traits subtypes:
  typedef typename Base::PPM                                   PPM;
  typedef typename Base::Kernel                                Kernel;
  typedef typename Base::NT                                    NT;
  typedef typename Base::Point_3                               Point_3;
  typedef typename Base::Vector_3                              Vector_3;

  // Solver traits subtypes:
  typedef typename Solver_traits::Vector                       Vector;
  typedef typename Solver_traits::Matrix                       Matrix;

// Public operations
public:
  /// Constructor
  Discrete_authalic_parameterizer_3(Border_parameterizer border_param = Border_parameterizer(),
                                    ///< %Object that maps the surface's border to 2D space.
                                    Solver_traits sparse_la = Solver_traits())
                                    ///< Traits object to access a sparse linear system.
    : Fixed_border_parameterizer_3<TriangleMesh,
                                   Border_parameterizer,
                                   Solver_traits>(border_param, sparse_la)
  { }

  // Default copy constructor and operator =() are fine

// Protected operations
protected:
  /// Compute w_ij = (i, j), coefficient of matrix A for j neighbor vertex of i.
  ///
  /// \param mesh a triangulated surface.
  /// \param main_vertex_v_i the vertex of `mesh` with index `i`
  /// \param neighbor_vertex_v_j the vertex of `mesh` with index `j`
  virtual NT compute_w_ij(const TriangleMesh& mesh,
                          vertex_descriptor main_vertex_v_i,
                          vertex_around_target_circulator neighbor_vertex_v_j) const
  {
    const PPM ppmap = get(vertex_point, mesh);

    const Point_3& position_v_i = get(ppmap, main_vertex_v_i);
    const Point_3& position_v_j = get(ppmap, *neighbor_vertex_v_j);

    // Compute the square norm of v_j -> v_i vector
    Vector_3 edge = position_v_i - position_v_j;
    double square_len = edge*edge;

    // Compute cotangent of (v_k,v_j,v_i) corner (i.e. cotan of v_j corner)
    // if v_k is the vertex before v_j when circulating around v_i
    vertex_around_target_circulator previous_vertex_v_k = neighbor_vertex_v_j;
    previous_vertex_v_k--;
    const Point_3& position_v_k = get(ppmap, *previous_vertex_v_k);
    NT cotg_psi_ij = internal::cotangent<Kernel>(position_v_k, position_v_j, position_v_i);

    // Compute cotangent of (v_i,v_j,v_l) corner (i.e. cotan of v_j corner)
    // if v_l is the vertex after v_j when circulating around v_i
    vertex_around_target_circulator next_vertex_v_l = neighbor_vertex_v_j;
    next_vertex_v_l++;
    const Point_3& position_v_l = get(ppmap,*next_vertex_v_l);
    NT cotg_theta_ij = internal::cotangent<Kernel>(position_v_i, position_v_j, position_v_l);

    NT weight = 0.0;
    CGAL_assertion(square_len != 0.0); // two points are identical!
    if(square_len != 0.0)
      weight = (cotg_psi_ij + cotg_theta_ij) / square_len;

    return weight;
  }
};

} // namespace Surface_mesh_parameterization

} // namespace CGAL

#endif // CGAL_SURFACE_MESH_PARAMETERIZATION_DISCRETE_AUTHALIC_PARAMETERIZER_3_H
