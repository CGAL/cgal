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
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Laurent Saboret, Pierre Alliez, Bruno Levy

#ifndef CGAL_SURFACE_MESH_PARAMETERIZATION_MEAN_VALUE_COORDINATES_PARAMETERIZER_3_H
#define CGAL_SURFACE_MESH_PARAMETERIZATION_MEAN_VALUE_COORDINATES_PARAMETERIZER_3_H

#include <CGAL/license/Surface_mesh_parameterization.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/Surface_mesh_parameterization/internal/validity.h>

#include <CGAL/Surface_mesh_parameterization/Error_code.h>
#include <CGAL/Surface_mesh_parameterization/Fixed_border_parameterizer_3.h>

#ifdef CGAL_EIGEN3_ENABLED
#include <CGAL/Eigen_solver_traits.h>
#endif

/// \file Mean_value_coordinates_parameterizer_3.h

namespace CGAL {

namespace Surface_mesh_parameterization {

/// \ingroup  PkgSurfaceParameterizationMethods
///
/// The class `Mean_value_coordinates_parameterizer_3`
/// implements *Floater Mean Value Coordinates* parameterization.
/// This method is sometimes called simply *Floater parameterization* \cgalCite{cgal:f-mvc-03}.
///
/// This is a conformal parameterization, i.e. it attempts to preserve angles.
///
/// One-to-one mapping is guaranteed if the surface's border is mapped to a convex polygon.
///
/// This class is used by the main parameterization algorithm
/// `Fixed_border_parameterizer_3::parameterize()`.
/// - It provides default `BorderParameterizer_` and `SolverTraits_` template
///   parameters.
/// - It implements `compute_w_ij()` to compute w_ij = (i, j) coefficient of matrix A
///   for j neighbor vertex of i based on Floater Mean Value Coordinates parameterization.
/// - It implements an optimized version of `is_one_to_one_mapping()`.
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
///
template < class TriangleMesh_,
           class BorderParameterizer_ = Default,
           class SolverTraits_ = Default>
class Mean_value_coordinates_parameterizer_3
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
    #pragma message("Error: You must either provide 'SolverTraits_' or link CGAL with the Eigen library");
         SolverTraits_>::type > // no parameter provided, and Eigen is not enabled: so don't compile!
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
                                      Solver_traits>          Base;

// Private types
private:
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor   vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor     face_descriptor;

  typedef typename boost::graph_traits<TriangleMesh>::vertex_iterator     vertex_iterator;
  typedef typename boost::graph_traits<TriangleMesh>::face_iterator       face_iterator;
  typedef CGAL::Vertex_around_target_circulator<TriangleMesh> vertex_around_target_circulator;

  // Mesh_TriangleMesh_3 subtypes:
  typedef typename Base::PPM              PPM;
  typedef typename Base::Kernel           Kernel;
  typedef typename Base::NT               NT;
  typedef typename Base::Point_3          Point_3;
  typedef typename Base::Vector_3         Vector_3;

  // Solver traits subtypes:
  typedef typename Solver_traits::Vector      Vector;
  typedef typename Solver_traits::Matrix      Matrix;

// Public operations
public:
  /// Constructor
  Mean_value_coordinates_parameterizer_3(Border_parameterizer border_param = Border_parameterizer(),
                                         ///< Object that maps the surface's border to 2D space.
                                         Solver_traits sparse_la = Solver_traits())
                                         ///< Traits object to access a sparse linear system.
  : Fixed_border_parameterizer_3<TriangleMesh,
                                 Border_parameterizer,
                                 Solver_traits>(border_param, sparse_la)
  { }

    // Default copy constructor and operator =() are fine

  /// Check if the 3D -> 2D mapping is one-to-one.
  template <typename VertexUVMap>
  bool is_one_to_one_mapping(const TriangleMesh& mesh,
                             halfedge_descriptor bhd,
                             const VertexUVMap uvmap) const
  {
    /// Theorem: A one-to-one mapping is guaranteed if all w_ij coefficients
    ///          are > 0 (for j vertex neighbor of i) and if the surface
    ///          border is mapped onto a 2D convex polygon.
    /// Here, all w_ij coefficients are positive (for j vertex neighbor of i), thus a
    /// valid embedding is guaranteed if the surface border is mapped
    /// onto a 2D convex polygon.
    return (Base::get_border_parameterizer().is_border_convex() ||
            internal::is_one_to_one_mapping(mesh, bhd, uvmap));
  }

// Protected operations
protected:
  /// Compute w_ij = (i, j) coefficient of matrix A for j neighbor vertex of i.
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

    // Compute the norm of v_j -> v_i vector
    Vector_3 edge = position_v_i - position_v_j;
    NT len = std::sqrt(edge * edge);

    // Compute angle of (v_j,v_i,v_k) corner (i.e. angle of v_i corner)
    // if v_k is the vertex before v_j when circulating around v_i
    vertex_around_target_circulator previous_vertex_v_k = neighbor_vertex_v_j;
    previous_vertex_v_k--;
    const Point_3& position_v_k = get(ppmap, *previous_vertex_v_k);
    NT gamma_ij = internal::compute_angle_rad<Kernel>(position_v_j, position_v_i, position_v_k);

    // Compute angle of (v_l,v_i,v_j) corner (i.e. angle of v_i corner)
    // if v_l is the vertex after v_j when circulating around v_i
    vertex_around_target_circulator next_vertex_v_l = neighbor_vertex_v_j;
    next_vertex_v_l++;
    const Point_3& position_v_l = get(ppmap, *next_vertex_v_l);
    NT delta_ij = internal::compute_angle_rad<Kernel>(position_v_l, position_v_i, position_v_j);

    NT weight = 0.0;
    CGAL_assertion(len != 0.0); // two points are identical!
    if(len != 0.0)
      weight = (std::tan(0.5*gamma_ij) + std::tan(0.5*delta_ij)) / len;
    CGAL_assertion(weight > 0);

    return weight;
  }
};

} // namespace Surface_mesh_parameterization

} // namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_SURFACE_MESH_PARAMETERIZATION_MEAN_VALUE_COORDINATES_PARAMETERIZER_3_H
