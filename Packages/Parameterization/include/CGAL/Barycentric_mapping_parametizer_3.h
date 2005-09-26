// Copyright (c) 2005  INRIA (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$
// $Name$
//
// Author(s)     : Laurent Saboret, Pierre Alliez


#ifndef CGAL_BARYCENTRIC_MAPPING_PARAMETIZER_3_H
#define CGAL_BARYCENTRIC_MAPPING_PARAMETIZER_3_H

#include <CGAL/Fixed_border_parametizer_3.h>
#include <CGAL/parameterization_assertions.h>

CGAL_BEGIN_NAMESPACE


/// Class Barycentric_mapping_parametizer_3
/// implements Tutte's barycentric mapping.
/// 1 to 1 mapping is guaranteed if surface's border is mapped to a convex polygon.
///
/// Concept: Model of the ParametizerTraits_3 concept.
///
/// Design pattern:
/// ParametizerTraits_3 models are Strategies (see [GOF95]): they implement
/// a strategy of surface parameterization for models of MeshAdaptor_3.

template
<
    class MeshAdaptor_3,              ///< 3D surface mesh
    class BorderParametizer_3         ///< Strategy to parameterize the surface border
                = Circular_border_arc_length_parametizer_3<MeshAdaptor_3>,
    class SparseLinearAlgebraTraits_d ///< Traits class to solve a sparse linear system
                = OpenNL::DefaultLinearSolverTraits<typename MeshAdaptor_3::NT>
                                      ///< Note: the sparse linear system is symmetric iff
                                      ///< Fixed_border_parametizer_3 removes fixed vertices
>
class Barycentric_mapping_parametizer_3
    : public Fixed_border_parametizer_3<MeshAdaptor_3,
                                        BorderParametizer_3,
                                        SparseLinearAlgebraTraits_d>
{
// Private types
private:

    /// Superclass
    typedef Fixed_border_parametizer_3<MeshAdaptor_3,
                                        BorderParametizer_3,
                                        SparseLinearAlgebraTraits_d>
                                            Base;

// Public types
public:
    // Export Mesh_Adaptor_3, BorderParametizer_3
    // and SparseLinearAlgebraTraits_d types
    typedef MeshAdaptor_3                   Adaptor;
    typedef typename Parametizer_traits_3<Adaptor>::Error_code
                                            Error_code;
    typedef typename Adaptor::NT            NT;
    typedef typename Adaptor::Facet_handle  Facet_handle;
    typedef typename Adaptor::Facet_const_handle
                                            Facet_const_handle;
    typedef typename Adaptor::Vertex_handle Vertex_handle;
    typedef typename Adaptor::Vertex_const_handle
                                            Vertex_const_handle;
    typedef typename Adaptor::Point_3       Point_3;
    typedef typename Adaptor::Point_2       Point_2;
    typedef typename Adaptor::Vector_3      Vector_3;
    typedef typename Adaptor::Vector_2      Vector_2;
    typedef typename Adaptor::Facet_iterator Facet_iterator;
    typedef typename Adaptor::Facet_const_iterator
                                            Facet_const_iterator;
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
    typedef BorderParametizer_3             Border_param;
    typedef SparseLinearAlgebraTraits_d     Sparse_LA;
    typedef typename Sparse_LA::Vector      Vector;
    typedef typename Sparse_LA::Matrix      Matrix;

// Public operations
public:
    /// Constructor
    Barycentric_mapping_parametizer_3(Border_param border_param = Border_param(),
                                            ///< Object that maps the surface's border to 2D space
                                      Sparse_LA sparse_la = Sparse_LA())
                                            ///< Traits object to access a sparse linear system
    :   Fixed_border_parametizer_3<Adaptor,
                                   Border_param,
                                   Sparse_LA>(border_param, sparse_la)
    {}

    // Default copy constructor and operator =() are fine

// Protected types
protected:
    typedef typename OpenNL::LinearSolver<Sparse_LA>
                                            Solver ;

// Protected operations
protected:
    /// compute wij = (i,j) coefficient of matrix A for j neighbor vertex of i
    virtual NT  compute_wij(const Adaptor& mesh,
                            Vertex_const_handle main_vertex_Vi,
                            Vertex_around_vertex_const_circulator neighbor_vertex_Vj)
    {
        /// Tutte algorithm is the most simple one: Wij = 1 for j neighbor vertex of i
        return 1;
    }

    /// Check if 3D -> 2D mapping is 1 to 1
    virtual bool  is_one_to_one_mapping (const Adaptor& mesh,
                                         const Matrix& A,
                                         const Vector& Bu,
                                         const Vector& Bv)
    {
        /// Theorem: 1 to 1 mapping is guaranteed if all Wij coefficients
        ///          are > 0 (for j vertex neighbor of i) and if the surface
        ///          boundary is mapped onto a 2D convex polygon.
        /// All Wij coefficients = 1 (for j vertex neighbor of i), thus mapping
        /// is guaranteed if the surface boundary is mapped onto a 2D convex polygon.
        return Base::get_border_parametizer().is_border_convex ();
    }
};


CGAL_END_NAMESPACE

#endif //CGAL_BARYCENTRIC_MAPPING_PARAMETIZER_3_H

