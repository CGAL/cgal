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


#ifndef CGAL_TWO_VERTICES_PARAMETIZER_3_H_INCLUDED
#define CGAL_TWO_VERTICES_PARAMETIZER_3_H_INCLUDED

#include <CGAL/parameterization_assertions.h>
#include <CGAL/Parametizer_traits_3.h>

#include <cfloat>
#include <climits>

CGAL_BEGIN_NAMESPACE


//
// Declaration
//

/// Class Two_vertices_parametizer_3
/// parameterizes 2 extreme vertices of a 3D surface.
/// This kind of border parameterization is used by free border parameterizations.
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
/// a strategy of boundary parameterization for models of MeshAdaptor_3.

template<class MeshAdaptor_3>           //< 3D surface
class Two_vertices_parametizer_3
{
// Public types
public:
    /// Export MeshAdaptor_3 template parameter
    typedef MeshAdaptor_3                   Adaptor;

    /// The various errors detected by this package
    typedef typename Parametizer_traits_3<Adaptor>::Error_code
                                            Error_code;

// Private types
private:
    // Export Mesh_Adaptor_3 subtypes
    //
    /// Number type to represent coordinates
    typedef typename Adaptor::NT            NT;
    /// 2D point that represents (u,v) coordinates computed
    /// by parameterization methods. Usual methods are expected.
    typedef typename Adaptor::Point_2       Point_2;
    /// 3D point that represents vertices coordinates. Usual methods are expected.
    typedef typename Adaptor::Point_3       Point_3;
    /// 2D vector. Usual methods are expected.
    typedef typename Adaptor::Vector_2      Vector_2;
    /// 3D vector. Usual methods are expected.
    typedef typename Adaptor::Vector_3      Vector_3;
    /// Opaque type representing a facet of the 3D mesh. No methods are expected.
    typedef typename Adaptor::Facet         Facet;
    /// Handle to a facet. Model of the Handle concept.
    typedef typename Adaptor::Facet_handle  Facet_handle;
    typedef typename Adaptor::Facet_const_handle
                                            Facet_const_handle;
    /// Iterator over all mesh facets. Model of the ForwardIterator concept.
    typedef typename Adaptor::Facet_iterator Facet_iterator;
    typedef typename Adaptor::Facet_const_iterator
                                            Facet_const_iterator;
    /// Opaque type representing a vertex of the 3D mesh. No methods are expected.
    typedef typename Adaptor::Vertex        Vertex;
    /// Handle to a vertex. Model of the Handle concept.
    typedef typename Adaptor::Vertex_handle Vertex_handle;
    typedef typename Adaptor::Vertex_const_handle
                                            Vertex_const_handle;
    /// Iterator over all vertices of a mesh. Model of the ForwardIterator concept.
    typedef typename Adaptor::Vertex_iterator Vertex_iterator;
    typedef typename Adaptor::Vertex_const_iterator
                                            Vertex_const_iterator;
    /// Iterator over vertices of the mesh "main boundary".
    /// Model of the ForwardIterator concept.
    typedef typename Adaptor::Border_vertex_iterator
                                            Border_vertex_iterator;
    typedef typename Adaptor::Border_vertex_const_iterator
                                            Border_vertex_const_iterator;
    /// Counter-clockwise circulator over a facet's vertices
    /// Model of the BidirectionalCirculator concept.
    typedef typename Adaptor::Vertex_around_facet_circulator
                                            Vertex_around_facet_circulator;
    typedef typename Adaptor::Vertex_around_facet_const_circulator
                                            Vertex_around_facet_const_circulator;
    /// Clockwise circulator over the vertices incident to a vertex
    /// Model of the BidirectionalCirculator concept.
    typedef typename Adaptor::Vertex_around_vertex_circulator
                                            Vertex_around_vertex_circulator;
    typedef typename Adaptor::Vertex_around_vertex_const_circulator
                                            Vertex_around_vertex_const_circulator;

// Public operations
public:
    // Default constructor, copy constructor and operator =() are fine

    /// Map 2 extreme vertices of the 3D mesh and mark them as "parameterized"
    Error_code parameterize_border (Adaptor* mesh);

    /// Indicate if border's shape is convex.
    /// Meaningless for free border parameterization algorithms.
    bool  is_border_convex () { return false; }
};


//
// Implementation
//

/// Map 2 extreme vertices of the 3D mesh and mark them as "parameterized".
/// Return false on error.
template<class Adaptor>
inline
typename Parametizer_traits_3<Adaptor>::Error_code
Two_vertices_parametizer_3<Adaptor>::parameterize_border(Adaptor* mesh)
{
    Vertex_iterator it;

    CGAL_parameterization_assertion(mesh != NULL);

    // Nothing to do if no boundary
    if (mesh->mesh_main_border_vertices_begin() == mesh->mesh_main_border_vertices_end())
    {
        std::cerr << "  error ERROR_INVALID_BOUNDARY!" << std::endl;
        return Parametizer_traits_3<Adaptor>::ERROR_INVALID_BOUNDARY;
    }

    std::cerr << "  map 2 vertices..." << std::endl;

    // Get mesh's bounding box
    double xmin =  1e30 ;
    double ymin =  1e30 ;
    double zmin =  1e30 ;
    double xmax = -1e30 ;
    double ymax = -1e30 ;
    double zmax = -1e30 ;
    for (it = mesh->mesh_vertices_begin(); it != mesh->mesh_vertices_end(); it++)
    {
        Point_3 position = mesh->get_vertex_position(it);

        xmin = std::min(position.x(), xmin) ;
        ymin = std::min(position.y(), ymin) ;
        zmin = std::min(position.z(), zmin) ;

        xmax = std::max(position.x(), xmax) ;
        ymax = std::max(position.y(), ymax) ;
        zmax = std::max(position.z(), zmax) ;
    }

    // Find longest bounding box axes
    double dx = xmax - xmin ;
    double dy = ymax - ymin ;
    double dz = zmax - zmin ;
    enum { X_AXIS, Y_AXIS, Z_AXIS } longest_axis, second_longest_axis;
    if(dx < dy && dx < dz) {
        if(dy > dz) {
            longest_axis        = Y_AXIS;
            second_longest_axis = Z_AXIS;
        } else {
            longest_axis        = Z_AXIS;
            second_longest_axis = Y_AXIS;
        }
    } else if(dy < dx && dy < dz) {
        if(dx > dz) {
            longest_axis        = X_AXIS;
            second_longest_axis = Z_AXIS;
        } else {
            longest_axis        = Z_AXIS;
            second_longest_axis = X_AXIS;
        }
        } else { // (dz < dx && dz < dy)
        if(dx > dy) {
            longest_axis        = X_AXIS;
            second_longest_axis = Y_AXIS;
        } else {
            longest_axis        = Y_AXIS;
            second_longest_axis = X_AXIS;
        }
    }
    Vector_3 V1,                // bounding box' longest axis
             V2 ;               // bounding box' 2nd longest axis
    double V1_min=0, V1_max=0;  // bounding box' dimensions along V1
    double V2_min=0, V2_max=0;  // bounding box' dimensions along V2
    switch (longest_axis)
    {
    case X_AXIS:
        V1 = Vector_3(1,0,0) ;
        V1_min = xmin;
        V1_max = xmax;
        break;
    case Y_AXIS:
        V1 = Vector_3(0,1,0) ;
        V1_min = ymin;
        V1_max = ymax;
        break;
    case Z_AXIS:
        V1 = Vector_3(0,0,1) ;
        V1_min = zmin;
        V1_max = zmax;
        break;
    default:
        CGAL_parameterization_assertion(false);
    }
    switch (second_longest_axis)
    {
    case X_AXIS:
        V2 = Vector_3(1,0,0) ;
        V2_min = xmin;
        V2_max = xmax;
        break;
    case Y_AXIS:
        V2 = Vector_3(0,1,0) ;
        V2_min = ymin;
        V2_max = ymax;
        break;
    case Z_AXIS:
        V2 = Vector_3(0,0,1) ;
        V2_min = zmin;
        V2_max = zmax;
        break;
    default:
        CGAL_parameterization_assertion(false);
    }

    // Project onto longest bounding box axes,
    // Set extrema vertices' (u,v) in unit square and mark them as "parameterized"
    Vertex_handle vxmin = NULL ;
    double  umin  = DBL_MAX ;
    Vertex_handle vxmax = NULL ;
    double  umax  = DBL_MIN ;
    for (it = mesh->mesh_vertices_begin(); it != mesh->mesh_vertices_end(); it++)
    {
        Point_3  position = mesh->get_vertex_position(it);
        Vector_3 position_as_vector = position - Point_3(0,0,0);

        // coordinate along the bounding box' main axes
        double u = position_as_vector * V1 ;
        double v = position_as_vector * V2 ;

        // convert to unit square coordinates
        CGAL_parameterization_assertion(V1_max > V1_min);
        CGAL_parameterization_assertion(V2_max > V2_min);
        u = (u - V1_min) / (V1_max - V1_min);
        v = (v - V2_min) / (V2_max - V2_min);

        mesh->set_vertex_uv(it, Point_2(u,v)) ; // useful only for vxmin and vxmax

        if(u < umin) {
            vxmin = it ;
            umin = u ;
        }
        if(u > umax) {
            vxmax = it ;
            umax = u ;
        }
    }
    mesh->set_vertex_parameterized(vxmin, true) ;
    mesh->set_vertex_parameterized(vxmax, true) ;
    std::cerr << "    #" << mesh->get_vertex_index(vxmin) << "(" << vxmin->vertex()->index() << ") parameterized " << std::endl;
    std::cerr << "    #" << mesh->get_vertex_index(vxmax) << "(" << vxmax->vertex()->index() << ") parameterized " << std::endl;

    std::cerr << "    done" << std::endl;

    return Parametizer_traits_3<Adaptor>::OK;
}


CGAL_END_NAMESPACE

#endif //CGAL_TWO_VERTICES_PARAMETIZER_3_H_INCLUDED

