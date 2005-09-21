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


#ifndef CGAL_CIRCULARBORDERPARAMETIZER_3_H
#define CGAL_CIRCULARBORDERPARAMETIZER_3_H

#include <CGAL/parameterization_assertions.h>
#include <CGAL/Parametizer_traits_3.h>

CGAL_BEGIN_NAMESPACE


//
// Class Circular_border_parametizer_3
//

/// Base class of strategies that parameterize the border
/// of a 3D surface onto a circle.
///
/// Implementation note:
/// To simplify the implementation, BorderParametizer_3 models know only the
/// MeshAdaptor_3 class. They don't know the parameterization algorithm
/// requirements nor the kind of sparse linear system used.
///
/// Concept: Model of the BorderParametizer_3 concept.
///
/// Design pattern:
/// BorderParametizer_3 models are Strategies (see [GOF95]): they implement
/// a strategy of boundary parameterization for models of MeshAdaptor_3

template<class MeshAdaptor_3>           //< 3D surface
class Circular_border_parametizer_3
{
// Public types
public:
    // Export Mesh_Adaptor_3 type and subtypes
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

// Public operations
public:
    /// Destructor of base class should be virtual
    virtual ~Circular_border_parametizer_3() {}

    // Default constructor, copy constructor and operator =() are fine

    /// Assign to mesh's border vertices a 2D position (ie a (u,v) pair)
    /// on border's shape. Mark them as "parameterized".
    Error_code parameterize_border (Adaptor* mesh);

    /// Indicate if border's shape is convex
    bool  is_border_convex () { return true; }

// Protected operations
protected:
    /// compute length of an edge
    virtual double compute_edge_length(const Adaptor& mesh,
                                       Vertex_const_handle source,
                                       Vertex_const_handle target) = 0;

// Private operations
private:
    /// compute total length of boundary
    double compute_boundary_length(const Adaptor& mesh);
};


/// compute  total length of boundary
template<class Adaptor>
inline
double Circular_border_parametizer_3<Adaptor>::compute_boundary_length(
                                                        const Adaptor& mesh)
{
    double len = 0.0;
    for(Border_vertex_const_iterator it = mesh.mesh_main_border_vertices_begin();
        it != mesh.mesh_main_border_vertices_end();
        it++)
    {
        CGAL_parameterization_assertion(mesh.is_vertex_on_main_border(it));

        // Get next iterator (looping)
        Border_vertex_const_iterator next = it;
        next++;
        if(next == mesh.mesh_main_border_vertices_end())
            next = mesh.mesh_main_border_vertices_begin();

        // Add 'length' of it -> next vector to 'len'
        len += compute_edge_length(mesh, it, next);
    }
    return len;
}

/// Assign to mesh's border vertices a 2D position (ie a (u,v) pair)
/// on border's shape. Mark them as "parameterized".
template<class Adaptor>
inline
typename Parametizer_traits_3<Adaptor>::Error_code
Circular_border_parametizer_3<Adaptor>::parameterize_border(Adaptor* mesh)
{
    CGAL_parameterization_assertion(mesh != NULL);

    // Nothing to do if no boundary
    if (mesh->mesh_main_border_vertices_begin() == mesh->mesh_main_border_vertices_end())
    {
        std::cerr << "  error ERROR_INVALID_BOUNDARY!" << std::endl;
        return Parametizer_traits_3<Adaptor>::ERROR_INVALID_BOUNDARY;
    }

    // compute the total boundary length
    double total_len = compute_boundary_length(*mesh);
    std::cerr << "  total boundary len: " << total_len << std::endl;
    if (total_len == 0)
    {
        std::cerr << "  error ERROR_INVALID_BOUNDARY!" << std::endl;
        return Parametizer_traits_3<Adaptor>::ERROR_INVALID_BOUNDARY;
    }

    std::cerr << "  map on a circle..." << std::endl;
    const double PI = 3.14159265359;
    const double tmp = 2*PI/total_len;
    double len = 0.0;           // current position on circle in [0, total_len]
    for(Border_vertex_iterator it = mesh->mesh_main_border_vertices_begin();
        it != mesh->mesh_main_border_vertices_end();
        it++)
    {
        CGAL_parameterization_assertion(mesh->is_vertex_on_main_border(it));

        double angle = len*tmp; // current position on the circle in radians

        // map vertex on unit circle
        Point_2 uv;
        uv = Point_2(0.5+0.5*CGAL_CLIB_STD::cos(-angle),0.5+0.5*CGAL_CLIB_STD::sin(-angle));
        mesh->set_vertex_uv(it, uv);

        // Mark vertex as "parameterized"
        mesh->set_vertex_parameterized(it, true);

        // Get next iterator (looping)
        Border_vertex_iterator next = it;
        next++;
        if(next == mesh->mesh_main_border_vertices_end())
            next = mesh->mesh_main_border_vertices_begin();

        // Add 'length' of it -> next vector to 'len'
        len += compute_edge_length(*mesh, it, next);
    }

    std::cerr << "    done" << std::endl;

    return Parametizer_traits_3<Adaptor>::OK;
}


//
// Class Circular_border_uniform_parametizer_3
//

/// This class parameterizes the border of a 3D surface onto a circle
/// on an uniform manner: points are equally spaced.
///
/// Concept: Model of the BorderParametizer_3 concept.

template<class MeshAdaptor_3>           //< 3D surface
class Circular_border_uniform_parametizer_3
    : public Circular_border_parametizer_3<MeshAdaptor_3>
{
// Public types
public:
    // Export Mesh_Adaptor_3 type and subtypes
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

// Public operations
public:
    // Default constructor, copy constructor and operator =() are fine

// Protected operations
protected:
    /// Compute length of an edge:
    virtual double compute_edge_length(const Adaptor& mesh,
                                       Vertex_const_handle source,
                                       Vertex_const_handle target)
    {
        /// Uniform boundary parameterization: points are equally spaced.
        return 1;
    }
};


//
// Class Circular_border_arc_length_parametizer_3
//

/// This class parameterizes the border of a 3D surface onto a circle
/// with an arc-length parameterization: (u,v) values are
/// proportional to the length of boundary edges.
///
/// Concept: Model of the BorderParametizer_3 concept.

template<class MeshAdaptor_3>           //< 3D surface
class Circular_border_arc_length_parametizer_3
    : public Circular_border_parametizer_3<MeshAdaptor_3>
{
// Public types
public:
    // Export Mesh_Adaptor_3 type and subtypes
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

// Public operations
public:
    // Default constructor, copy constructor and operator =() are fine

// Protected operations
protected:
    /// Compute length of an edge:
    virtual double compute_edge_length(const Adaptor& mesh,
                                       Vertex_const_handle source,
                                       Vertex_const_handle target)
    {
        /// Arc-length boundary parameterization: (u,v) values are
        /// proportional to the length of boundary edges.
        Vector_3 v = mesh.get_vertex_position(target)
                   - mesh.get_vertex_position(source);
        return std::sqrt(v*v);
    }
};


CGAL_END_NAMESPACE

#endif //CGAL_CIRCULARBORDERPARAMETIZER_3_H

