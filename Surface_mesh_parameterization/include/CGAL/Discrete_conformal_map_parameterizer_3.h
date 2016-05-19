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


#ifndef CGAL_DISCRETE_CONFORMAL_MAP_PARAMETERIZER_3_H
#define CGAL_DISCRETE_CONFORMAL_MAP_PARAMETERIZER_3_H

#include <CGAL/license/Surface_mesh_parameterization.h>


#include <CGAL/Fixed_border_parameterizer_3.h>
#include <CGAL/Eigen_solver_traits.h>

/// \file Discrete_conformal_map_parameterizer_3.h

namespace CGAL {


/// \ingroup  PkgSurfaceParameterizationMethods
///
/// The class Discrete_conformal_map_parameterizer_3
/// implements the <i>Discrete Conformal Map (DCM)</i> parameterization \cgalCite{cgal:eddhls-maam-95}.
/// This algorithm is also called <i>Discrete Conformal Parameterization (DCP)</i>,
/// <i>Discrete Harmonic Map</i> or <i>Fixed Conformal Parameterization</i> by other authors.
///
/// This is a conformal parameterization, i.e. it attempts to preserve angles.
///
/// One-to-one mapping is guaranteed if surface's border is mapped onto a convex polygon.
///
/// This class is used by the main
/// parameterization algorithm `Fixed_border_parameterizer_3::parameterize()`.
/// - It provides default `BorderParameterizer_3` and `SparseLinearAlgebraTraits_d` template
///   parameters that make sense.
/// - It implements `compute_w_ij()` to compute w_ij = (i, j) coefficient of matrix A
///   for j neighbor vertex of i based on Discrete Conformal Map method.
///
/// \cgalModels `ParameterizerTraits_3`
///
///
/// \param TriangleMesh       3D surface mesh.
/// \param BorderParameterizer_3        Strategy to parameterize the surface border.
/// \param SparseLinearAlgebraTraits_d  Traits class to solve a sparse linear system.
///        Note: the system is *not* symmetric because `Fixed_border_parameterizer_3`
///        does not remove (yet) border vertices from the system.
///
/// \sa `CGAL::Parameterizer_traits_3<TriangleMesh>`
/// \sa `CGAL::Fixed_border_parameterizer_3<TriangleMesh, BorderParameterizer_3, SparseLinearAlgebraTraits_d>`
/// \sa `CGAL::Barycentric_mapping_parameterizer_3<TriangleMesh, BorderParameterizer_3, SparseLinearAlgebraTraits_d>`
/// \sa `CGAL::Discrete_authalic_parameterizer_3<TriangleMesh, BorderParameterizer_3, SparseLinearAlgebraTraits_d>`
/// \sa `CGAL::LSCM_parameterizer_3<TriangleMesh, BorderParameterizer_3>`
/// \sa `CGAL::Mean_value_coordinates_parameterizer_3<TriangleMesh, BorderParameterizer_3, SparseLinearAlgebraTraits_d>`

template
<
    class TriangleMesh,
    class BorderParameterizer_3
                = Circular_border_arc_length_parameterizer_3<TriangleMesh>,
    class SparseLinearAlgebraTraits_d
                = Eigen_solver_traits<Eigen::BiCGSTAB<Eigen_sparse_matrix<double>::EigenType, Eigen::IncompleteLUT< double > > >
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
                                        SparseLinearAlgebraTraits_d>
                                        Base;

// Public types
public:
    // We have to repeat the types exported by superclass
    /// @cond SKIP_IN_MANUAL
    typedef typename Base::Error_code       Error_code;
    typedef BorderParameterizer_3           Border_param;
    typedef SparseLinearAlgebraTraits_d     Sparse_LA;
    /// @endcond

// Private types
private:
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor vertex_descriptor;
  typedef CGAL::Halfedge_around_target_circulator<TriangleMesh> halfedge_around_target_circulator;
  typedef CGAL::Vertex_around_target_circulator<TriangleMesh> vertex_around_target_circulator;
    // Mesh_TriangleMesh_3 subtypes:
  typedef typename Parameterizer_traits_3<TriangleMesh>::NT            NT;
  typedef typename Parameterizer_traits_3<TriangleMesh>::Point_3       Point_3;
  typedef typename Parameterizer_traits_3<TriangleMesh>::Vector_3      Vector_3;

  
    // SparseLinearAlgebraTraits_d subtypes:
    typedef typename Sparse_LA::Vector      Vector;
    typedef typename Sparse_LA::Matrix      Matrix;

    using Base::cotangent;

// Public operations
public:
    /// Constructor
    Discrete_conformal_map_parameterizer_3(Border_param border_param = Border_param(),
                                            ///< Object that maps the surface's border to 2D space.
                                           Sparse_LA sparse_la = Sparse_LA())
                                            ///< Traits object to access a sparse linear system.
    :   Fixed_border_parameterizer_3<TriangleMesh,
                                     Border_param,
                                     Sparse_LA>(border_param, sparse_la)
    {}

    // Default copy constructor and operator =() are fine

// Protected operations
protected:
    /// Compute w_ij = (i,j) coefficient of matrix A for j neighbor vertex of i.
    virtual NT compute_w_ij(const TriangleMesh& mesh,
                            vertex_descriptor main_vertex_v_i,
                            halfedge_around_target_circulator neighbor_vertex_v_j) // its target is main_vertex_v_i
    {
    typedef typename Parameterizer_traits_3<TriangleMesh>::VPM PPmap;
 
        PPmap ppmap = get(vertex_point, mesh);
        Point_3 position_v_i = get(ppmap,main_vertex_v_i);
        Point_3 position_v_j = get(ppmap,source(*neighbor_vertex_v_j,mesh));

        // Compute cotangent of (v_i,v_k,v_j) corner (i.e. cotan of v_k corner)
        // if v_k is the vertex before v_j when circulating around v_i
        halfedge_around_target_circulator previous_vertex_v_k = neighbor_vertex_v_j;
        previous_vertex_v_k --;
        Point_3 position_v_k = get(ppmap,source(*previous_vertex_v_k, mesh));
        double cotg_beta_ij  = cotangent(position_v_i, position_v_k, position_v_j);

        // Compute cotangent of (v_j,v_l,v_i) corner (i.e. cotan of v_l corner)
        // if v_l is the vertex after v_j when circulating around v_i
        halfedge_around_target_circulator next_vertex_v_l = neighbor_vertex_v_j;
        next_vertex_v_l ++;
        Point_3 position_v_l = get(ppmap,source(*next_vertex_v_l,mesh));
        double cotg_alpha_ij = cotangent(position_v_j, position_v_l, position_v_i);

        double weight = cotg_beta_ij+cotg_alpha_ij;
        return weight;
    }
};


} //namespace CGAL

#endif //CGAL_DISCRETE_CONFORMAL_MAP_PARAMETERIZER_3_H
