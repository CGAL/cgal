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


#ifndef CGAL_PARAMETIZER_3_H
#define CGAL_PARAMETIZER_3_H

#include <CGAL/Kernel/global_functions.h>
#include <CGAL/parameterization_assertions.h>

CGAL_BEGIN_NAMESPACE


/// Class Parametizer_traits_3
/// is the base class of all parameterization methods.
///
/// Concept:
/// Model of the ParametizerTraits_3 concept (although you cannot instanciate this class).

template<class MeshAdaptor_3>       //< 3D surface
class Parametizer_traits_3
{
// Public types
public:
    /// The various errors detected by this package
    enum Error_code {
    OK,
    ERROR_EMPTY_MESH,               ///< input mesh is empty
    ERROR_NON_TRIANGULAR_MESH,      ///< input mesh is not triangular
    ERROR_NO_SURFACE_MESH,          ///< input mesh is not a surface
    ERROR_INVALID_BOUNDARY,         ///< parameterization requires a convex border
    ERROR_BAD_MATRIX_CONDITIONING,  ///< result is mathematically unstable
    ERROR_CANNOT_SOLVE_LINEAR_SYSTEM,///< cannot solve linear system
    ERROR_NO_1_TO_1_MAPPING,        ///< parameterization does not ensure 1 to 1 mapping
    ERROR_NOT_ENOUGH_MEMORY,        ///< it's time to buy some RAM :-)
    ERROR_WRONG_PARAMETER           ///< a method received an unexpected parameter
    };

    // Export Mesh_Adaptor_3 type and subtypes
    typedef MeshAdaptor_3                   Adaptor;
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

// Public operations
public:
    /// Destructor of base class should be virtual
    virtual ~Parametizer_traits_3() {}

    // Default constructor, copy constructor and operator =() are fine

    /// Compute a 1 to 1 mapping from a 3D surface 'mesh'
    /// to a piece of the 2D space.
    /// The mapping is linear by pieces (linear in each triangle).
    /// The result is the (u,v) pair image of each vertex of the 3D surface.
    ///
    /// Preconditions:
    /// - 'mesh' must be a surface with 1 connected component
    /// - 'mesh' must be a triangular mesh
    virtual Error_code  parameterize (Adaptor* mesh) = 0;

// Protected operations
protected:
    //                                                  -> ->
    /// return cotangent of (P,Q,R) corner (ie cotan of QP,QR angle)
    double cotangent(const Point_3& P,
                     const Point_3& Q,
                     const Point_3& R)
    {
        Vector_3 u = P - Q;
        Vector_3 v = R - Q;
        // (u . v)/((u x v).len)
        double dot = (u*v);
        Vector_3 cross_vector = CGAL::cross_product(u,v);
        double cross_norm = std::sqrt(cross_vector*cross_vector);
        CGAL_parameterization_assertion(cross_norm != 0.0);
        if(cross_norm != 0.0)
            return (dot/cross_norm);
        else
            return 0.0; // undefined
    }

    //                                                  -> ->
    /// return tangent of (P,Q,R) corner (ie tangent of QP,QR angle)
    double tangent(const Point_3& P,
                   const Point_3& Q,
                   const Point_3& R)
    {
        Vector_3 u = P - Q;
        Vector_3 v = R - Q;
        // (u . v)/((u x v).len)
        double dot = (u*v);
        CGAL_parameterization_assertion(dot != 0.0);
        Vector_3 cross_vector = CGAL::cross_product(u,v);
        double cross_norm = std::sqrt(cross_vector*cross_vector);
        if(dot != 0.0)
            return (cross_norm/dot);
        else
            return 0.0; // undefined
    }

    //                                                     -> ->
    /// return angle (in radians) of of (P,Q,R) corner (ie QP,QR angle)
    static double compute_angle_rad(const Point_3& P,
                                    const Point_3& Q,
                                    const Point_3& R)
    {
        static const double PI = 3.14159265359;

        Vector_3 u = P - Q;
        Vector_3 v = R - Q;

        // check
        double product = std::sqrt(u*u) * std::sqrt(v*v);
        if(product == 0)
            return 0.0;

        // cosine
        double dot = (u*v);
        double cosine = dot / product;

        // sine
        Vector_3 w = CGAL::cross_product(u,v);
        double AbsSine = std::sqrt(w*w) / product;

        if(cosine >= 0)
            return std::asin(fix_sine(AbsSine));
        else
            return PI-std::asin(fix_sine(AbsSine));
    }

// Private operations
private:
    /// fix sine
    static double fix_sine(double sine)
    {
        if(sine >= 1)
            return 1;
        else if(sine <= -1)
            return -1;
        else
            return sine;
    }
};


CGAL_END_NAMESPACE

#endif //CGAL_PARAMETIZER_3_H
