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
//
//
// Author(s)     : Laurent Saboret, Pierre Alliez, Bruno Levy


#ifndef CGAL_BARYCENTRIC_MAPPING_PARAMETERIZER_3_H
#define CGAL_BARYCENTRIC_MAPPING_PARAMETERIZER_3_H

#include <CGAL/license/Surface_mesh_parameterization.h>


#include <CGAL/Fixed_border_parameterizer_3.h>
#include <CGAL/Eigen_solver_traits.h>

/// \file Barycentric_mapping_parameterizer_3.h

namespace CGAL {

/// \ingroup  PkgSurfaceParameterizationMethods
///
/// The class Barycentric_mapping_parameterizer_3 implements <i>Tutte Barycentric Mapping algorithm</i> \cgalCite{t-hdg-63}.
/// This algorithm is also called <i>Tutte Uniform Weights</i> by other authors.
///
/// One-to-one mapping is guaranteed if the surface's border is mapped to a convex polygon.
///
/// This class is used by the main
/// parameterization algorithm Fixed_border_parameterizer_3::parameterize().
/// - It provides default BorderParameterizer_3 and SparseLinearAlgebraTraits_d template
///   parameters.
/// - It implements compute_w_ij() to compute `w_ij = (i,j)` coefficient of matrix A
///   for `j` neighbor vertex of `i` based on Tutte Barycentric Mapping method.
/// - It implements an optimized version of is_one_to_one_mapping().
///
/// \cgalModels `ParameterizerTraits_3`
///
///
/// \tparam TriangleMesh       a model of `FaceGraph`
/// \tparam BorderParameterizer_3        Strategy to parameterize the surface border.
/// \tparam SparseLinearAlgebraTraits_d  Traits class to solve a sparse linear system.
///        Note: the system is *not* symmetric because `Fixed_border_parameterizer_3`
///        does not remove (yet) border vertices from the system.

/*!
\sa `CGAL::Parameterizer_traits_3<TriangleMesh>`
\sa `CGAL::Fixed_border_parameterizer_3<TriangleMesh, BorderParameterizer_3, SparseLinearAlgebraTraits_d>`
\sa `CGAL::Discrete_authalic_parameterizer_3<TriangleMesh, BorderParameterizer_3, SparseLinearAlgebraTraits_d>`
\sa `CGAL::Discrete_conformal_map_parameterizer_3<TriangleMesh, BorderParameterizer_3, SparseLinearAlgebraTraits_d>`
\sa `CGAL::LSCM_parameterizer_3<TriangleMesh, BorderParameterizer_3>`
\sa `CGAL::Mean_value_coordinates_parameterizer_3<TriangleMesh, BorderParameterizer_3, SparseLinearAlgebraTraits_d>`
 */

template
<
  class TriangleMesh,
  class BorderParameterizer_3
    = Circular_border_arc_length_parameterizer_3<TriangleMesh>,
  class SparseLinearAlgebraTraits_d
    = Eigen_solver_traits<Eigen::BiCGSTAB<Eigen_sparse_matrix<double>::EigenType,
                                          Eigen::IncompleteLUT< double > > >
>
class Barycentric_mapping_parameterizer_3
  : public Fixed_border_parameterizer_3<TriangleMesh,
                                        BorderParameterizer_3,
                                        SparseLinearAlgebraTraits_d>
{
// Private types
private:
  // Superclass
  typedef Fixed_border_parameterizer_3<TriangleMesh,
                                      BorderParameterizer_3,
                                      SparseLinearAlgebraTraits_d>  Base;

// Public types
public:
  // We have to repeat the types exported by superclass
  /// @cond SKIP_IN_MANUAL
  typedef typename Base::Error_code       Error_code;
  typedef BorderParameterizer_3           Border_param;
  typedef SparseLinearAlgebraTraits_d     Sparse_LA;
  /// @endcond

protected:
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor  vertex_descriptor;
  typedef CGAL::Vertex_around_target_circulator<TriangleMesh>  vertex_around_target_circulator;
  typedef typename Parameterizer_traits_3<TriangleMesh>::NT  NT;

  // SparseLinearAlgebraTraits_d subtypes:
  typedef typename Sparse_LA::Vector      Vector;
  typedef typename Sparse_LA::Matrix      Matrix;

// Public operations
public:
  /// Constructor
  Barycentric_mapping_parameterizer_3(Border_param border_param = Border_param(),
                                      ///< object that maps the surface's border to 2D space.
                                      Sparse_LA sparse_la = Sparse_LA())
                                      ///< Traits object to access a sparse linear system.
  : Fixed_border_parameterizer_3<TriangleMesh,
                                 Border_param,
                                 Sparse_LA>(border_param, sparse_la)
  { }

  // Default copy constructor and operator =() are fine

// Protected operations
protected:
  /// Compute w_ij = (i,j) coefficient of matrix A for j neighbor vertex of i.
  virtual NT compute_w_ij(const TriangleMesh& /* mesh */,
      vertex_descriptor /* main_vertex_v_i */,
      vertex_around_target_circulator /* neighbor_vertex_v_j */ )
  {
    /// Tutte Barycentric Mapping algorithm is the most simple one:
    /// w_ij = 1 for j neighbor vertex of i.
    return 1;
  }
};

} // namespace CGAL

#endif // CGAL_BARYCENTRIC_MAPPING_PARAMETERIZER_3_H
