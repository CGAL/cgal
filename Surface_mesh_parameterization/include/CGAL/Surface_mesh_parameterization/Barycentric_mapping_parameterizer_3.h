// Copyright (c) 2005  INRIA (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Laurent Saboret, Pierre Alliez, Bruno Levy

#ifndef CGAL_SURFACE_MESH_PARAMETERIZATION_BARYCENTRIC_MAPPING_PARAMETERIZER_3_H
#define CGAL_SURFACE_MESH_PARAMETERIZATION_BARYCENTRIC_MAPPING_PARAMETERIZER_3_H

#include <CGAL/license/Surface_mesh_parameterization.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/Surface_mesh_parameterization/internal/validity.h>
#include <CGAL/Surface_mesh_parameterization/Circular_border_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/Fixed_border_parameterizer_3.h>
#include <CGAL/Weights/uniform_weights.h>

#include <CGAL/Default.h>
#include <CGAL/iterator.h>

#if defined(CGAL_EIGEN3_ENABLED)
#include <CGAL/Eigen_solver_traits.h>
#endif

/// \file Barycentric_mapping_parameterizer_3.h

namespace CGAL {

namespace Surface_mesh_parameterization {

/// \ingroup PkgSurfaceMeshParameterizationMethods
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
/// - implements compute_w_ij() to compute `w_ij`, the `(i,j)`-coefficient of
///   the matrix `A` for `j` neighbor vertex of `i`, based on Tutte Barycentric
///   Mapping method.
///
/// \cgalModels `Parameterizer_3`
///
/// \tparam TriangleMesh_ must be a model of `FaceGraph`.
///
/// \tparam BorderParameterizer_ is a strategy to parameterize the surface border
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
  /// Border parameterizer type
  typedef Border_parameterizer_                               Border_parameterizer;

  /// Solver traits type
  typedef SolverTraits_                                       Solver_traits;
#endif

  /// Triangle mesh type
  typedef TriangleMesh_                                       Triangle_mesh;

  typedef TriangleMesh_                                       TriangleMesh;

  /// Mesh vertex type
  typedef typename boost::graph_traits<Triangle_mesh>::vertex_descriptor    vertex_descriptor;

  /// Mesh halfedge type
  typedef typename boost::graph_traits<Triangle_mesh>::halfedge_descriptor  halfedge_descriptor;

private:
  // Superclass
  typedef Fixed_border_parameterizer_3<Triangle_mesh,
                                       Border_parameterizer,
                                       Solver_traits>  Base;

// Private types
private:
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
  : Fixed_border_parameterizer_3<Triangle_mesh,
                                 Border_parameterizer,
                                 Solver_traits>(border_param, sparse_la)
  { }

  // Default copy constructor and operator =() are fine

  /// returns whether the 3D -> 2D mapping is one-to-one.
  template <typename VertexUVMap>
  bool is_one_to_one_mapping(const Triangle_mesh& mesh,
                             halfedge_descriptor bhd,
                             const VertexUVMap uvmap) const
  {
    /// Theorem: A one-to-one mapping is guaranteed if all `w_ij` coefficients
    ///          are > 0 (for `j` vertex neighbor of `i`) and if the surface
    ///          border is mapped onto a 2D convex polygon.
    /// Here, all `w_ij` coefficients are equal to `1` (for `j` vertex neighbor of `i`), thus a
    /// valid embedding is guaranteed if the surface border is mapped
    /// onto a 2D convex polygon.
    return (Base::get_border_parameterizer().is_border_convex() ||
            internal::is_one_to_one_mapping(mesh, bhd, uvmap));
  }

// Protected operations
protected:
  /// computes `w_ij`, the coefficient of matrix `A` for `j` neighbor vertex of `i`.
  virtual NT compute_w_ij(const Triangle_mesh& /* mesh */,
      vertex_descriptor /* main_vertex_v_i */,
      Vertex_around_target_circulator<Triangle_mesh> /* neighbor_vertex_v_j */ ) const
  {
    /// In the Tutte Barycentric Mapping algorithm, we have `w_ij = 1`, for `j` neighbor vertex of `i`.
    return NT(1);
  }
};

} // namespace Surface_mesh_parameterization

} // namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_SURFACE_MESH_PARAMETERIZATION_BARYCENTRIC_MAPPING_PARAMETERIZER_3_H
