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

#ifndef CGAL_SURFACE_MESH_PARAMETERIZATION_BARYCENTRIC_MAPPING_PARAMETERIZER_3_H
#define CGAL_SURFACE_MESH_PARAMETERIZATION_BARYCENTRIC_MAPPING_PARAMETERIZER_3_H

#include <CGAL/license/Surface_mesh_parameterization.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/Surface_mesh_parameterization/internal/validity.h>
#include <CGAL/Surface_mesh_parameterization/Circular_border_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/Fixed_border_parameterizer_3.h>

#include <CGAL/Default.h>

#if defined(CGAL_EIGEN3_ENABLED)
#include <CGAL/Eigen_solver_traits.h>
#endif

/// \file Barycentric_mapping_parameterizer_3.h

namespace CGAL {

namespace Surface_mesh_parameterization {

/// \ingroup PkgSurfaceParameterizationMethods
///
/// The class `Barycentric_mapping_parameterizer_3` implements <i>Tutte Barycentric
/// Mapping algorithm</i>. This algorithm is also called
/// <i>Tutte Uniform Weights</i> by other authors \cgalCite{t-hdg-63}.
///
/// A one-to-one mapping is guaranteed if the surface's border is mapped to a convex polygon.
///
/// This class is a strategy called by the main
/// parameterization algorithm `Fixed_border_parameterizer_3::parameterize()` and it:
/// - provides the template parameters `BorderParameterizer_` and `SolverTraits_`.
/// - implements compute_w_ij() to compute `w_ij = (i,j)`, coefficient of
///   the matrix A for `j` neighbor vertex of `i`, based on Tutte Barycentric
///   Mapping method.
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
///      CGAL::Eigen_solver_traits<
///              Eigen::BiCGSTAB<Eigen_sparse_matrix<double>::EigenType,
///                              Eigen::IncompleteLUT< double > > >
/// \endcode
///
/// \sa `CGAL::Surface_mesh_parameterization::Fixed_border_parameterizer_3<TriangleMesh, BorderParameterizer, SolverTraits>`
///
template < typename TriangleMesh_,
           typename BorderParameterizer_ = Default,
           typename SolverTraits_ = Default >
class Barycentric_mapping_parameterizer_3
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
  typedef SolverTraits_                                       SolverTraits;
#endif

  typedef TriangleMesh_                                       TriangleMesh;

private:
  // Superclass
  typedef Fixed_border_parameterizer_3<TriangleMesh_,
                                       Border_parameterizer,
                                       Solver_traits>  Base;

// Private types
private:
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor    vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor  halfedge_descriptor;
  typedef CGAL::Vertex_around_target_circulator<TriangleMesh>              vertex_around_target_circulator;

  typedef typename Base::NT                       NT;

  // Solver traits subtypes:
  typedef typename Solver_traits::Vector          Vector;
  typedef typename Solver_traits::Matrix          Matrix;

// Public operations
public:
  /// Constructor
  Barycentric_mapping_parameterizer_3(Border_parameterizer border_param = Border_parameterizer(),
                                      ///< %Object that maps the surface's border to 2D space.
                                      Solver_traits sparse_la = Solver_traits())
                                      ///< Traits object to access a sparse linear system.
  : Fixed_border_parameterizer_3<TriangleMesh,
                                 Border_parameterizer,
                                 Solver_traits>(border_param, sparse_la)
  { }

  // Default copy constructor and operator =() are fine

  /// Check if the 3D -> 2D mapping is one-to-one.
  template <typename VertexUVMap,
            typename Faces_Container>
  bool is_one_to_one_mapping(const TriangleMesh& mesh,
                             halfedge_descriptor bhd,
                             const VertexUVMap uvmap) const
  {
    /// Theorem: A one-to-one mapping is guaranteed if all w_ij coefficients
    ///          are > 0 (for j vertex neighbor of i) and if the surface
    ///          border is mapped onto a 2D convex polygon.
    /// Here, all w_ij coefficients = 1 (for j vertex neighbor of i), thus a
    /// valid embedding is guaranteed if the surface border is mapped
    /// onto a 2D convex polygon.
    return (Base::get_border_parameterizer().is_border_convex() ||
            internal::is_one_to_one_mapping(mesh, bhd, uvmap));
  }

// Protected operations
protected:
  /// Compute w_ij = (i,j), coefficient of matrix A for j neighbor vertex of i.
  virtual NT compute_w_ij(const TriangleMesh& /* mesh */,
      vertex_descriptor /* main_vertex_v_i */,
      vertex_around_target_circulator /* neighbor_vertex_v_j */ ) const
  {
    /// In the Tutte Barycentric Mapping algorithm, we have w_ij = 1,
    /// for j neighbor vertex of i.
    return 1.;
  }
};

} // namespace Surface_mesh_parameterization

} // namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_SURFACE_MESH_PARAMETERIZATION_BARYCENTRIC_MAPPING_PARAMETERIZER_3_H
