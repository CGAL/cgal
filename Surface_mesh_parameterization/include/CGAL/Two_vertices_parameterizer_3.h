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


#ifndef CGAL_TWO_VERTICES_PARAMETERIZER_3_H_INCLUDED
#define CGAL_TWO_VERTICES_PARAMETERIZER_3_H_INCLUDED

#include <CGAL/surface_mesh_parameterization_assertions.h>
#include <CGAL/Parameterizer_traits_3.h>

#include <cfloat>
#include <climits>

/// \file Two_vertices_parameterizer_3.h

namespace CGAL {


//
// Declaration
//


/// \ingroup  PkgSurfaceParameterizationBorderParameterizationMethods
///
/// The class `Two_vertices_parameterizer_3`
/// parameterizes two extreme vertices of a 3D surface.
/// This kind of border parameterization is used by free border parameterizations.
///
/// Implementation note:
/// To simplify the implementation, `BorderParameterizer_3` models know only the
/// `ParameterizationMesh_3` class. They do not know the parameterization algorithm
/// requirements or the kind of sparse linear system used.
///
/// \cgalModels `BorderParameterizer_3`
///

template<class ParameterizationMesh_3>      //< 3D surface
class Two_vertices_parameterizer_3
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
    // Default constructor, copy constructor and operator =() are fine.

    /// Map two extreme vertices of the 3D mesh and mark them as <i>parameterized</i>.
    typename Parameterizer_traits_3<Adaptor>::Error_code
                                        parameterize_border(Adaptor& mesh);

    /// Indicate if border's shape is convex.
    /// Meaningless for free border parameterization algorithms.
    bool  is_border_convex () { return false; }
};


//
// Implementation
//

// Map two extreme vertices of the 3D mesh and mark them as "parameterized".
template<class Adaptor>
inline
typename Parameterizer_traits_3<Adaptor>::Error_code
Two_vertices_parameterizer_3<Adaptor>::parameterize_border(Adaptor& mesh)
{
    Vertex_iterator it;

    // Nothing to do if no border
    if (mesh.mesh_main_border_vertices_begin() == mesh.mesh_main_border_vertices_end())
        return Parameterizer_traits_3<Adaptor>::ERROR_BORDER_TOO_SHORT;

    // Get mesh's bounding box
    double xmin = (std::numeric_limits<double>::max)() ;
    double ymin = (std::numeric_limits<double>::max)() ;
    double zmin = (std::numeric_limits<double>::max)() ;
    double xmax = (std::numeric_limits<double>::min)() ;
    double ymax = (std::numeric_limits<double>::min)() ;
    double zmax = (std::numeric_limits<double>::min)() ;
    for (it = mesh.mesh_vertices_begin(); it != mesh.mesh_vertices_end(); it++)
    {
        Point_3 position = mesh.get_vertex_position(it);

        xmin = (std::min)(position.x(), xmin) ;
        ymin = (std::min)(position.y(), ymin) ;
        zmin = (std::min)(position.z(), zmin) ;

        xmax = (std::max)(position.x(), xmax) ;
        ymax = (std::max)(position.y(), ymax) ;
        zmax = (std::max)(position.z(), zmax) ;
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
        CGAL_surface_mesh_parameterization_assertion(false);
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
        CGAL_surface_mesh_parameterization_assertion(false);
    }

    // Project onto longest bounding box axes,
    // Set extrema vertices' (u,v) in unit square and mark them as "parameterized"
    Vertex_handle vxmin = NULL ;
    double  umin  =  (std::numeric_limits<double>::max)() ;
    double vmin =  (std::numeric_limits<double>::max)(), vmax=  (std::numeric_limits<double>::min)();
    Vertex_handle vxmax = NULL ;
    double  umax  =  (std::numeric_limits<double>::min)() ;
    for (it = mesh.mesh_vertices_begin(); it != mesh.mesh_vertices_end(); it++)
    {
        Point_3  position = mesh.get_vertex_position(it);
        Vector_3 position_as_vector = position - Point_3(0,0,0);

        // coordinate along the bounding box' main axes
        double u = position_as_vector * V1 ;
        double v = position_as_vector * V2 ;

        // convert to unit square coordinates
        CGAL_surface_mesh_parameterization_assertion(V1_max > V1_min);
        CGAL_surface_mesh_parameterization_assertion(V2_max > V2_min);
        u = (u - V1_min) / (V1_max - V1_min);
        v = (v - V2_min) / (V2_max - V2_min);

        mesh.set_vertex_uv(it, Point_2(u,v)) ; // useful only for vxmin and vxmax

        if(u < umin || (u==umin && v < vmin) ) {
            vxmin = it ;
            umin = u ;
            vmin = v ;
        }
        if(u > umax || (u==umax && v > vmax) ){
            vxmax = it ;
            umax = u ;
            vmax = v ;
        }
    }
    mesh.set_vertex_parameterized(vxmin, true) ;
    mesh.set_vertex_parameterized(vxmax, true) ;

#ifdef DEBUG_TRACE
    std::cerr << "  map two vertices..." << std::endl;
    // std::cerr << "    #" << mesh.get_vertex_index(vxmin) << "(" << vxmin->vertex()->index() << ") parameterized " << std::endl;
    // std::cerr << "    #" << mesh.get_vertex_index(vxmax) << "(" << vxmax->vertex()->index() << ") parameterized " << std::endl;
#endif

    return Parameterizer_traits_3<Adaptor>::OK;
}


} //namespace CGAL

#endif //CGAL_TWO_VERTICES_PARAMETERIZER_3_H_INCLUDED
