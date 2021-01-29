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

#ifndef CGAL_SURFACE_MESH_PARAMETERIZATION_DISCRETE_CONFORMAL_MAP_PARAMETERIZER_3_H
#define CGAL_SURFACE_MESH_PARAMETERIZATION_DISCRETE_CONFORMAL_MAP_PARAMETERIZER_3_H

#include <CGAL/license/Surface_mesh_parameterization.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/Surface_mesh_parameterization/internal/angles.h>
#include <CGAL/Surface_mesh_parameterization/internal/kernel_traits.h>
#include <CGAL/Surface_mesh_parameterization/Error_code.h>

#include <CGAL/Surface_mesh_parameterization/Fixed_border_parameterizer_3.h>

#if defined(CGAL_EIGEN3_ENABLED)
#include <CGAL/Eigen_solver_traits.h>
#endif

#include <CGAL/Default.h>

/// \file Discrete_conformal_map_parameterizer_3.h

namespace CGAL {

namespace Surface_mesh_parameterization {

/// \ingroup  PkgSurfaceMeshParameterizationMethods
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
/// - provides the template parameters `BorderParameterizer_` and `SolverTraits_`.
/// - implements `compute_w_ij()` to compute `w_ij`, the `(i,j)`-coefficient of matrix `A`,
///   for `j` neighbor vertex of `i`, based on Discrete Conformal Map method.
///
/// \cgalModels `Parameterizer_3`
///
/// \tparam TriangleMesh_ must be a model of `FaceGraph`
///
/// \tparam BorderParameterizer_ is a Strategy to parameterize the surface border
///         and must be a model of `Parameterizer_3`.<br>
///         <b>%Default:</b>
/// \code
///   Circular_border_arc_length_parameterizer_3<TriangleMesh_>.
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
class Discrete_conformal_map_parameterizer_3
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

  // Private types
private:
    // Superclass
  typedef Fixed_border_parameterizer_3<Triangle_mesh,
                                       Border_parameterizer,
                                       Solver_traits>    Base;

// Private types
private:
  typedef typename boost::graph_traits<Triangle_mesh>::vertex_descriptor vertex_descriptor;
  typedef CGAL::Vertex_around_target_circulator<Triangle_mesh> vertex_around_target_circulator;

  // Traits subtypes:
  typedef typename Base::Kernel             Kernel;
  typedef typename Base::NT                 NT;
  typedef typename Base::Point_3            Point_3;
  typedef typename Base::Vector_3           Vector_3;
  typedef typename Base::PPM                PPM;

  // Solver traits subtypes:
  typedef typename Solver_traits::Vector    Vector;
  typedef typename Solver_traits::Matrix    Matrix;

// Public operations
public:
  /// Constructor
  Discrete_conformal_map_parameterizer_3(Border_parameterizer border_param = Border_parameterizer(),
                                         ///< %Object that maps the surface's border to 2D space.
                                         Solver_traits sparse_la = Solver_traits())
                                         ///< Traits object to access a sparse linear system.
  :   Fixed_border_parameterizer_3<Triangle_mesh,
                                   Border_parameterizer,
                                   Solver_traits>(border_param, sparse_la)
  { }

  // Default copy constructor and operator =() are fine

// Protected operations
protected:
  /// computes `w_ij`, the `(i,j)`-coefficient of matrix `A`, for `j` neighbor vertex of `i`.
  ///
  /// \param mesh a triangulated surface.
  /// \param main_vertex_v_i the vertex of `mesh` with index `i`
  /// \param neighbor_vertex_v_j the vertex of `mesh` with index `j`
  virtual NT compute_w_ij(const Triangle_mesh& mesh,
                          vertex_descriptor main_vertex_v_i,
                          Vertex_around_target_circulator<Triangle_mesh> neighbor_vertex_v_j) const // its target is main_vertex_v_i
  {
    const PPM ppmap = get(vertex_point, mesh);

    const Point_3& position_v_i = get(ppmap, main_vertex_v_i);
    const Point_3& position_v_j = get(ppmap, *neighbor_vertex_v_j);

    // Compute cotangent of (v_i,v_k,v_j) corner (i.e. cotan of v_k corner)
    // if v_k is the vertex before v_j when circulating around v_i
    vertex_around_target_circulator previous_vertex_v_k = neighbor_vertex_v_j;
    previous_vertex_v_k--;
    const Point_3& position_v_k = get(ppmap, *previous_vertex_v_k);
    NT cotg_beta_ij = internal::cotangent<Kernel>(position_v_i, position_v_k, position_v_j);

    // Compute cotangent of (v_j,v_l,v_i) corner (i.e. cotan of v_l corner)
    // if v_l is the vertex after v_j when circulating around v_i
    vertex_around_target_circulator next_vertex_v_l = neighbor_vertex_v_j;
    next_vertex_v_l++;
    const Point_3& position_v_l = get(ppmap, *next_vertex_v_l);
    NT cotg_alpha_ij = internal::cotangent<Kernel>(position_v_j, position_v_l, position_v_i);

    NT weight = cotg_beta_ij + cotg_alpha_ij;
    return weight;
  }
};

} // namespace Surface_mesh_parameterization

} // namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_SURFACE_MESH_PARAMETERIZATION_DISCRETE_CONFORMAL_MAP_PARAMETERIZER_3_H
