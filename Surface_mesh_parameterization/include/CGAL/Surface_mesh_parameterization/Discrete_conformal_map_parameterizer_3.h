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

#ifndef CGAL_SURFACE_MESH_PARAMETERIZATION_DISCRETE_CONFORMAL_MAP_PARAMETERIZER_3_H
#define CGAL_SURFACE_MESH_PARAMETERIZATION_DISCRETE_CONFORMAL_MAP_PARAMETERIZER_3_H

#include <CGAL/license/Surface_mesh_parameterization.h>

#include <CGAL/Surface_mesh_parameterization/internal/angles.h>
#include <CGAL/Surface_mesh_parameterization/internal/kernel_traits.h>
#include <CGAL/Surface_mesh_parameterization/Error_code.h>

#include <CGAL/Surface_mesh_parameterization/Fixed_border_parameterizer_3.h>

#include <CGAL/Eigen_solver_traits.h>

/// \file Discrete_conformal_map_parameterizer_3.h

namespace CGAL {

namespace Surface_mesh_parameterization {

/// \ingroup  PkgSurfaceParameterizationMethods
///
/// The class `Discrete_conformal_map_parameterizer_3`
/// implements the <i>Discrete Conformal Map (DCM)</i> parameterization.
/// This algorithm is also called <i>Discrete Conformal Parameterization (DCP)</i>,
/// <i>Discrete Harmonic Map</i> or <i>Fixed Conformal Parameterization</i> by other authors \cgalCite{cgal:eddhls-maam-95}.
///
/// This is a conformal parameterization, i.e. it attempts to preserve angles.
///
/// A one-to-one mapping is guaranteed if the surface's border is mapped onto a convex polygon.
///
/// This class is a strategy called by the main
/// parameterization algorithm `Fixed_border_parameterizer_3::parameterize()` and it:
/// - provides the template parameters `BorderParameterizer_3` and `SparseLinearAlgebraTraits_d`.
/// - implements `compute_w_ij()` to compute w_ij = (i, j), coefficient of matrix A
///   for j neighbor vertex of i, based on Discrete Conformal Map method.
///
/// \cgalModels `Parameterizer_3`
///
/// \param TriangleMesh must be a model of `FaceGraph`.
/// \param BorderParameterizer_3 is a strategy to parameterize the surface border
///        and must be a model of `Parameterizer_3`.
/// \param SparseLinearAlgebraTraits_d is a traits class to solve a sparse linear system. <br>
///        Note: the system is *not* symmetric because `Fixed_border_parameterizer_3`
///        does not remove border vertices from the system.
///
/// \sa `CGAL::Surface_mesh_parameterization::Fixed_border_parameterizer_3<TriangleMesh, BorderParameterizer_3, SparseLinearAlgebraTraits_d>`
/// \sa `CGAL::Surface_mesh_parameterization::ARAP_parameterizer_3<TriangleMesh, BorderParameterizer_3, SparseLinearAlgebraTraits_d>`
/// \sa `CGAL::Surface_mesh_parameterization::Barycentric_mapping_parameterizer_3<TriangleMesh, BorderParameterizer_3, SparseLinearAlgebraTraits_d>`
/// \sa `CGAL::Surface_mesh_parameterization::Discrete_authalic_parameterizer_3<TriangleMesh, BorderParameterizer_3, SparseLinearAlgebraTraits_d>`
/// \sa `CGAL::Surface_mesh_parameterization::LSCM_parameterizer_3<TriangleMesh, BorderParameterizer_3>`
/// \sa `CGAL::Surface_mesh_parameterization::Mean_value_coordinates_parameterizer_3<TriangleMesh, BorderParameterizer_3, SparseLinearAlgebraTraits_d>`
///
template
<
  class TriangleMesh,
  class BorderParameterizer_3
    = Circular_border_arc_length_parameterizer_3<TriangleMesh>,
  class SparseLinearAlgebraTraits_d
    = Eigen_solver_traits<Eigen::BiCGSTAB<Eigen_sparse_matrix<double>::EigenType,
                                          Eigen::IncompleteLUT< double > > >
>
class Discrete_conformal_map_parameterizer_3
  : public Fixed_border_parameterizer_3<TriangleMesh,
                                        BorderParameterizer_3,
                                        SparseLinearAlgebraTraits_d>
{
// Private types
private:
  // Superclass
  typedef Fixed_border_parameterizer_3<TriangleMesh,
                                      BorderParameterizer_3,
                                      SparseLinearAlgebraTraits_d>    Base;

// Public types
public:
    // We have to repeat the types exported by superclass
    /// @cond SKIP_IN_MANUAL
    typedef BorderParameterizer_3           Border_param;
    typedef SparseLinearAlgebraTraits_d     Sparse_LA;
    /// @endcond

// Private types
private:
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor vertex_descriptor;
  typedef CGAL::Vertex_around_target_circulator<TriangleMesh> vertex_around_target_circulator;

  // Traits subtypes:
  typedef typename Base::Kernel             Kernel;
  typedef typename Base::NT                 NT;
  typedef typename Base::Point_3            Point_3;
  typedef typename Base::Vector_3           Vector_3;
  typedef typename Base::PPM                PPM;

  // SparseLinearAlgebraTraits_d subtypes:
  typedef typename Sparse_LA::Vector        Vector;
  typedef typename Sparse_LA::Matrix        Matrix;

// Public operations
public:
  /// Constructor
  Discrete_conformal_map_parameterizer_3(Border_param border_param = Border_param(),
                                         ///< %Object that maps the surface's border to 2D space.
                                         Sparse_LA sparse_la = Sparse_LA())
                                         ///< Traits object to access a sparse linear system.
  :   Fixed_border_parameterizer_3<TriangleMesh,
                                   Border_param,
                                   Sparse_LA>(border_param, sparse_la)
  { }

  // Default copy constructor and operator =() are fine

// Protected operations
protected:
  /// Compute w_ij = (i,j) coefficient of matrix A for j neighbor vertex of i.
  ///
  /// \param mesh a triangulated surface.
  /// \param main_vertex_v_i the vertex of `mesh` with index `i`
  /// \param neighbor_vertex_v_j the vertex of `mesh` with index `j`
  virtual NT compute_w_ij(const TriangleMesh& mesh,
                          vertex_descriptor main_vertex_v_i,
                          vertex_around_target_circulator neighbor_vertex_v_j) const // its target is main_vertex_v_i
  {
    const PPM ppmap = get(vertex_point, mesh);

    Point_3 position_v_i = get(ppmap, main_vertex_v_i);
    Point_3 position_v_j = get(ppmap, *neighbor_vertex_v_j);

    // Compute cotangent of (v_i,v_k,v_j) corner (i.e. cotan of v_k corner)
    // if v_k is the vertex before v_j when circulating around v_i
    vertex_around_target_circulator previous_vertex_v_k = neighbor_vertex_v_j;
    previous_vertex_v_k--;
    Point_3 position_v_k = get(ppmap, *previous_vertex_v_k);
    NT cotg_beta_ij = internal::cotangent<Kernel>(position_v_i, position_v_k, position_v_j);

    // Compute cotangent of (v_j,v_l,v_i) corner (i.e. cotan of v_l corner)
    // if v_l is the vertex after v_j when circulating around v_i
    vertex_around_target_circulator next_vertex_v_l = neighbor_vertex_v_j;
    next_vertex_v_l++;
    Point_3 position_v_l = get(ppmap, *next_vertex_v_l);
    NT cotg_alpha_ij = internal::cotangent<Kernel>(position_v_j, position_v_l, position_v_i);

    NT weight = cotg_beta_ij + cotg_alpha_ij;
    return weight;
  }
};

} // namespace Surface_mesh_parameterization

} // namespace CGAL

#endif // CGAL_SURFACE_MESH_PARAMETERIZATION_DISCRETE_CONFORMAL_MAP_PARAMETERIZER_3_H
