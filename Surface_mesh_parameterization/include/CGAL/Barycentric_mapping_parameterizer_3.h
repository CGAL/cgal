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

#include <CGAL/Fixed_border_parameterizer_3.h>
#include <CGAL/surface_mesh_parameterization_assertions.h>

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
///   parameters that make sense.
/// - It implements compute_w_ij() to compute `w_ij = (i,j)` coefficient of matrix A
///   for `j` neighbor vertex of `i` based on Tutte Barycentric Mapping method.
/// - It implements an optimized version of is_one_to_one_mapping().
///
/// \cgalModels `ParameterizerTraits_3`
///
///
/// \tparam ParameterizationMesh_3       3D surface mesh.
/// \tparam BorderParameterizer_3        Strategy to parameterize the surface border.
/// \tparam SparseLinearAlgebraTraits_d  Traits class to solve a sparse linear system.
///        Note: the system is *not* symmetric because `Fixed_border_parameterizer_3`
///        does not remove (yet) border vertices from the system.

/*!
\sa `CGAL::Parameterizer_traits_3<ParameterizationMesh_3>`
\sa `CGAL::Fixed_border_parameterizer_3<ParameterizationMesh_3, BorderParameterizer_3, SparseLinearAlgebraTraits_d>`
\sa `CGAL::Discrete_authalic_parameterizer_3<ParameterizationMesh_3, BorderParameterizer_3, SparseLinearAlgebraTraits_d>`
\sa `CGAL::Discrete_conformal_map_parameterizer_3<ParameterizationMesh_3, BorderParameterizer_3, SparseLinearAlgebraTraits_d>`
\sa `CGAL::LSCM_parameterizer_3<ParameterizationMesh_3, BorderParameterizer_3, SparseLinearAlgebraTraits_d>`
\sa `CGAL::Mean_value_coordinates_parameterizer_3<ParameterizationMesh_3, BorderParameterizer_3, SparseLinearAlgebraTraits_d>`
 */

template
<
    class ParameterizationMesh_3,
    class BorderParameterizer_3
                = Circular_border_arc_length_parameterizer_3<ParameterizationMesh_3>,
    class SparseLinearAlgebraTraits_d
                = OpenNL::DefaultLinearSolverTraits<typename ParameterizationMesh_3::NT>
>
class Barycentric_mapping_parameterizer_3
    : public Fixed_border_parameterizer_3<ParameterizationMesh_3,
                                        BorderParameterizer_3,
                                        SparseLinearAlgebraTraits_d>
{
// Private types
private:
    // Superclass
    typedef Fixed_border_parameterizer_3<ParameterizationMesh_3,
                                        BorderParameterizer_3,
                                        SparseLinearAlgebraTraits_d>
                                            Base;

// Public types
public:
    // We have to repeat the types exported by superclass
    /// @cond SKIP_IN_MANUAL
    typedef typename Base::Error_code       Error_code;
    typedef ParameterizationMesh_3          Adaptor;
    typedef BorderParameterizer_3           Border_param;
    typedef SparseLinearAlgebraTraits_d     Sparse_LA;
    /// @endcond

// Private types
private:
    // Mesh_Adaptor_3 subtypes:
    typedef typename Adaptor::NT            NT;
    typedef typename Adaptor::Point_2       Point_2;
    typedef typename Adaptor::Point_3       Point_3;
    typedef typename Adaptor::Vector_2      Vector_2;
    typedef typename Adaptor::Vector_3      Vector_3;
    typedef typename Adaptor::Facet         Facet;
    typedef typename Adaptor::Facet_handle  Facet_handle;
    typedef typename Adaptor::Facet_const_handle
                                            Facet_const_handle;
    typedef typename Adaptor::Facet_iterator Facet_iterator;
    typedef typename Adaptor::Facet_const_iterator
                                            Facet_const_iterator;
    typedef typename Adaptor::Vertex        Vertex;
    typedef typename Adaptor::Vertex_handle Vertex_handle;
    typedef typename Adaptor::Vertex_const_handle
                                            Vertex_const_handle;
    typedef typename Adaptor::Vertex_iterator Vertex_iterator;
    typedef typename Adaptor::Vertex_const_iterator
                                            Vertex_const_iterator;
    typedef typename Adaptor::Border_vertex_iterator
                                            Border_vertex_iterator;
    typedef typename Adaptor::Border_vertex_const_iterator
                                            Border_vertex_const_iterator;
    typedef typename Adaptor::Vertex_around_facet_circulator
                                            Vertex_around_facet_circulator;
    typedef typename Adaptor::Vertex_around_facet_const_circulator
                                            Vertex_around_facet_const_circulator;
    typedef typename Adaptor::Vertex_around_vertex_circulator
                                            Vertex_around_vertex_circulator;
    typedef typename Adaptor::Vertex_around_vertex_const_circulator
                                            Vertex_around_vertex_const_circulator;

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
    :   Fixed_border_parameterizer_3<Adaptor,
                                   Border_param,
                                   Sparse_LA>(border_param, sparse_la)
    {}

    // Default copy constructor and operator =() are fine

// Protected operations
protected:
    /// Compute w_ij = (i,j) coefficient of matrix A for j neighbor vertex of i.
    virtual NT compute_w_ij(const Adaptor& /* mesh */,
			  Vertex_const_handle /* main_vertex_v_i */,
			  Vertex_around_vertex_const_circulator /* neighbor_vertex_v_j */ )
    {
        /// Tutte Barycentric Mapping algorithm is the most simple one:
        /// w_ij = 1 for j neighbor vertex of i.
        return 1;
    }

    /// Check if 3D -> 2D mapping is one-to-one.
    virtual bool  is_one_to_one_mapping (const Adaptor& /* mesh */,
				       const Matrix& /* A */,
				       const Vector& /* Bu */,
				       const Vector& /* Bv */)
    {
        /// Theorem: one-to-one mapping is guaranteed if all w_ij coefficients
        ///          are > 0 (for j vertex neighbor of i) and if the surface
        ///          border is mapped onto a 2D convex polygon.
        /// All w_ij coefficients = 1 (for j vertex neighbor of i), thus mapping
        /// is guaranteed if the surface border is mapped onto a 2D convex polygon.
        return Base::get_border_parameterizer().is_border_convex ();
    }
};


} //namespace CGAL

#endif //CGAL_BARYCENTRIC_MAPPING_PARAMETERIZER_3_H
