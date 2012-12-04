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


#ifndef CGAL_CIRCULARBORDERPARAMETERIZER_3_H
#define CGAL_CIRCULARBORDERPARAMETERIZER_3_H

#include <CGAL/surface_mesh_parameterization_assertions.h>
#include <CGAL/Parameterizer_traits_3.h>

/// \file Circular_border_parameterizer_3.h

namespace CGAL {


//
// Class Circular_border_parameterizer_3
//

/// \ingroup  PkgSurfaceParameterizationBorderParameterizationMethods
///
/// This is the base class of strategies that parameterize the border
/// of a 3D surface onto a circle.
/// `Circular_border_parameterizer_3` is a pure virtual class, thus
/// cannot be instantiated.
/// It implements most of the algorithm. Subclasses just
/// have to implement `compute_edge_length()` to compute a segment's length.
///
/// Implementation note:
/// To simplify the implementation, `BorderParameterizer_3` models know only the
/// `ParameterizationMesh_3` class. They do not know the parameterization algorithm
/// requirements or the kind of sparse linear system used.
///
/// \cgalModels `BorderParameterizer_3`
///
/// \sa `CGAL::Circular_border_arc_length_parameterizer_3<ParameterizationMesh_3>`
/// \sa `CGAL::Circular_border_uniform_parameterizer_3<ParameterizationMesh_3>`

template<class ParameterizationMesh_3>           //< 3D surface
class Circular_border_parameterizer_3
{
// Public types
public:
    /// Export ParameterizationMesh_3 template parameter
    typedef ParameterizationMesh_3          Adaptor;

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

// Public operations
public:
    /// Destructor of base class should be virtual.
    virtual ~Circular_border_parameterizer_3() {}

    // Default constructor, copy constructor and operator =() are fine

    /// Assign to mesh's border vertices a 2D position (i.e.\ a (u,v) pair)
    /// on border's shape. Mark them as <i>parameterized</i>.
    typename Parameterizer_traits_3<Adaptor>::Error_code
                                        parameterize_border(Adaptor& mesh);

    /// Indicate if border's shape is convex.
    bool  is_border_convex () { return true; }

// Protected operations
protected:
    /// Compute the length of an edge.
    virtual double compute_edge_length(const Adaptor& mesh,
                                       Vertex_const_handle source,
                                       Vertex_const_handle target) = 0;

// Private operations
private:
    /// Compute the total length of the border
    double compute_border_length(const Adaptor& mesh);
};


// Compute the total length of the border
template<class Adaptor>
inline
double Circular_border_parameterizer_3<Adaptor>::compute_border_length(
                                                        const Adaptor& mesh)
{
    double len = 0.0;
    for(Border_vertex_const_iterator it = mesh.mesh_main_border_vertices_begin();
        it != mesh.mesh_main_border_vertices_end();
        it++)
    {
        CGAL_surface_mesh_parameterization_assertion(mesh.is_vertex_on_main_border(it));

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

// Assign to mesh's border vertices a 2D position (i.e.\ a (u,v) pair)
// on border's shape. Mark them as "parameterized".
template<class Adaptor>
inline
typename Parameterizer_traits_3<Adaptor>::Error_code
Circular_border_parameterizer_3<Adaptor>::parameterize_border(Adaptor& mesh)
{
#ifdef DEBUG_TRACE
    std::cerr << "  map on a circle" << std::endl;
#endif

    // Nothing to do if no border
    if (mesh.mesh_main_border_vertices_begin() == mesh.mesh_main_border_vertices_end())
        return Parameterizer_traits_3<Adaptor>::ERROR_BORDER_TOO_SHORT;

    // Compute the total border length
    double total_len = compute_border_length(mesh);
    if (total_len == 0)
        return Parameterizer_traits_3<Adaptor>::ERROR_BORDER_TOO_SHORT;

    const double PI = 3.14159265359;
    const double tmp = 2*PI/total_len;
    double len = 0.0;           // current position on circle in [0, total_len]
    for(Border_vertex_iterator it = mesh.mesh_main_border_vertices_begin();
        it != mesh.mesh_main_border_vertices_end();
        it++)
    {
        CGAL_surface_mesh_parameterization_assertion(mesh.is_vertex_on_main_border(it));

        double angle = len*tmp; // current position on the circle in radians

        // map vertex on unit circle
        Point_2 uv;
        uv = Point_2(0.5+0.5*std::cos(-angle),0.5+0.5*std::sin(-angle));
        mesh.set_vertex_uv(it, uv);

        // Mark vertex as "parameterized"
        mesh.set_vertex_parameterized(it, true);

        // Get next iterator (looping)
        Border_vertex_iterator next = it;
        next++;
        if(next == mesh.mesh_main_border_vertices_end())
            next = mesh.mesh_main_border_vertices_begin();

        // Add 'length' of it -> next vector to 'len'
        len += compute_edge_length(mesh, it, next);
    }

    return Parameterizer_traits_3<Adaptor>::OK;
}


//
// Class Circular_border_uniform_parameterizer_3
//

/// \ingroup  PkgSurfaceParameterizationBorderParameterizationMethods
///
/// This class parameterizes the border of a 3D surface onto a circle
/// in a uniform manner: points are equally spaced.
/// Circular_border_parameterizer_3 implements most of the border parameterization
/// algorithm. This class implements only `compute_edge_length()` to compute a
/// segment's length.
///
/// \cgalModels `BorderParameterizer_3`
///
/// \sa `CGAL::Circular_border_arc_length_parameterizer_3<ParameterizationMesh_3>`
/// \sa `CGAL::Circular_border_parameterizer_3<ParameterizationMesh_3>`

template<class ParameterizationMesh_3>      //< 3D surface
class Circular_border_uniform_parameterizer_3
    : public Circular_border_parameterizer_3<ParameterizationMesh_3>
{
// Public types
public:
    // We have to repeat the types exported by superclass
    /// @cond SKIP_IN_MANUAL
    typedef ParameterizationMesh_3          Adaptor;
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

// Public operations
public:
    // Default constructor, copy constructor and operator =() are fine

// Protected operations
protected:
    /// Compute the length of an edge.
    virtual double compute_edge_length(const Adaptor& /* mesh */,
                                       Vertex_const_handle /* source */,
                                       Vertex_const_handle /* target */)
    {
        /// Uniform border parameterization: points are equally spaced.
        return 1;
    }
};


//
// Class Circular_border_arc_length_parameterizer_3
//

/// \ingroup  PkgSurfaceParameterizationBorderParameterizationMethods
///
/// This class parameterizes the border of a 3D surface onto a circle,
/// with an arc-length parameterization: (u,v) values are
/// proportional to the length of border edges.
/// `Circular_border_parameterizer_3` implements most of the border parameterization
/// algorithm. This class implements only `compute_edge_length()` to compute a
/// segment's length.
///
/// \cgalModels `BorderParameterizer_3`
///
/// \sa `CGAL::Circular_border_parameterizer_3<ParameterizationMesh_3>`
/// \sa `CGAL::Circular_border_uniform_parameterizer_3<ParameterizationMesh_3>`

template<class ParameterizationMesh_3>           //< 3D surface
class Circular_border_arc_length_parameterizer_3
    : public Circular_border_parameterizer_3<ParameterizationMesh_3>
{
// Public types
public:
    // We have to repeat the types exported by superclass
    /// @cond SKIP_IN_MANUAL
    typedef ParameterizationMesh_3          Adaptor;
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

// Public operations
public:
    // Default constructor, copy constructor and operator =() are fine

// Protected operations
protected:
    /// Compute the length of an edge.
    virtual double compute_edge_length(const Adaptor& mesh,
                                       Vertex_const_handle source,
                                       Vertex_const_handle target)
    {
        /// Arc-length border parameterization: (u,v) values are
        /// proportional to the length of border edges.
        Vector_3 v = mesh.get_vertex_position(target)
                   - mesh.get_vertex_position(source);
        return std::sqrt(v*v);
    }
};


} //namespace CGAL

#endif //CGAL_CIRCULARBORDERPARAMETERIZER_3_H
