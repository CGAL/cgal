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


#ifndef CGAL_PARAMETERIZER_3_H
#define CGAL_PARAMETERIZER_3_H

#include <CGAL/Kernel/global_functions.h>
#include <CGAL/surface_mesh_parameterization_assertions.h>

/// \file Parameterizer_traits_3.h

namespace CGAL {


/// \ingroup  PkgSurfaceParameterizationMethods
///
/// The class `Parameterizer_traits_3`
/// is the base class of all parameterization methods.
/// This class is a pure virtual class, thus cannot be instantiated.
///
/// This class doesn't do much. Its main goal is to ensure that subclasses
/// will be proper models of the `ParameterizerTraits_3` concept:
/// - `Parameterizer_traits_3` defines the Error_code list of errors detected by this package
/// - `Parameterizer_traits_3` declares a pure virtual method parameterize()
///
/// \cgalModels `ParameterizerTraits_3`
///
/// ## Design Pattern ##
/// `ParameterizerTraits_3` models are *Strategies*: they implement
/// a strategy of surface parameterization for models of `ParameterizationMesh_3`.
///
///
/// \sa `CGAL::Fixed_border_parameterizer_3<ParameterizationMesh_3, BorderParameterizer_3, SparseLinearAlgebraTraits_d>`
/// \sa `CGAL::Barycentric_mapping_parameterizer_3<ParameterizationMesh_3, BorderParameterizer_3, SparseLinearAlgebraTraits_d>`
/// \sa `CGAL::Discrete_authalic_parameterizer_3<ParameterizationMesh_3, BorderParameterizer_3, SparseLinearAlgebraTraits_d>`
/// \sa `CGAL::Discrete_conformal_map_parameterizer_3<ParameterizationMesh_3, BorderParameterizer_3, SparseLinearAlgebraTraits_d>`
/// \sa `CGAL::LSCM_parameterizer_3<ParameterizationMesh_3, BorderParameterizer_3, SparseLinearAlgebraTraits_d>`
/// \sa `CGAL::Mean_value_coordinates_parameterizer_3<ParameterizationMesh_3, BorderParameterizer_3, SparseLinearAlgebraTraits_d>`

template<class ParameterizationMesh_3>       //< 3D surface
class Parameterizer_traits_3
{
// Public types
public:
    /// List of errors detected by this package
    enum Error_code {
    OK,                             ///< Success
    ERROR_EMPTY_MESH,               ///< Input mesh is empty
    ERROR_NON_TRIANGULAR_MESH,      ///< Input mesh is not triangular
    ERROR_NO_TOPOLOGICAL_DISC,      ///< Input mesh is not a topological disc
    ERROR_BORDER_TOO_SHORT,         ///< This border parameterization requires a longer border
    ERROR_NON_CONVEX_BORDER,        ///< This parameterization method requires a convex border
    ERROR_CANNOT_SOLVE_LINEAR_SYSTEM,///< Cannot solve linear system
    ERROR_NO_1_TO_1_MAPPING,        ///< Parameterization failed: no one-to-one mapping
    ERROR_OUT_OF_MEMORY,            ///< Not enough memory
    ERROR_WRONG_PARAMETER           ///< A method received an unexpected parameter
    };

    /// Export ParameterizationMesh_3 template parameter
    typedef ParameterizationMesh_3          Adaptor;

// Protected types
protected:
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

// Public operations
public:
    /// Destructor of base class should be virtual.
    virtual ~Parameterizer_traits_3() {}

    // Default constructor, copy constructor and operator =() are fine

    /// Compute a one-to-one mapping from a 3D surface mesh
    /// to a piece of the 2D space.
    /// The mapping is linear by pieces (linear in each triangle).
    /// The result is the (u,v) pair image of each vertex of the 3D surface.
    ///
    /// \pre `mesh` must be a surface with one connected component.
    /// \pre `mesh` must be a triangular mesh.
    virtual Error_code  parameterize (Adaptor& mesh) = 0;

    /// Get message corresponding to an error code
    /// \param error_code The code returned by `parameterize()`
    /// \return           The string describing the error code
    static const char* get_error_message(int error_code)
    {
      // Messages corresponding to Error_code list above. Must be kept in sync!
      static const char* error_message[ERROR_WRONG_PARAMETER+1] = {
        "Success",
        "Input mesh is empty",
        "Input mesh is not triangular",
        "Input mesh is not a topological disc",
        "This border parameterization requires a longer border",
        "This parameterization method requires a convex border",
        "Cannot solve linear system",
        "Parameterization failed: no one-to-one mapping",
        "Not enough memory",
        "A method received an unexpected parameter"
      };

      if(error_code > ERROR_WRONG_PARAMETER || error_code < 0)
        return "Unknown error";
      else
        return error_message[error_code];
    };

// Protected operations
/// @cond SKIP_IN_MANUAL
protected:
    //                                                    -> ->
    /// Return cotangent of (P,Q,R) corner (i.e. cotan of QP,QR angle).
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
        CGAL_surface_mesh_parameterization_assertion(cross_norm != 0.0);
        if(cross_norm != 0.0)
            return (dot/cross_norm);
        else
            return 0.0; // undefined
    }

    //                                                    -> ->
    /// Return tangent of (P,Q,R) corner (i.e. tangent of QP,QR angle).
    double tangent(const Point_3& P,
                   const Point_3& Q,
                   const Point_3& R)
    {
        Vector_3 u = P - Q;
        Vector_3 v = R - Q;
        // (u . v)/((u x v).len)
        double dot = (u*v);
        CGAL_surface_mesh_parameterization_assertion(dot != 0.0);
        Vector_3 cross_vector = CGAL::cross_product(u,v);
        double cross_norm = std::sqrt(cross_vector*cross_vector);
        if(dot != 0.0)
            return (cross_norm/dot);
        else
            return 0.0; // undefined
    }

    //                                                       -> ->
    /// Return angle (in radians) of of (P,Q,R) corner (i.e. QP,QR angle).
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
/// @endcond

// Private operations
private:
    /// Fix sine.
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


} //namespace CGAL

#endif //CGAL_PARAMETERIZER_3_H
