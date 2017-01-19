// Copyright (c) 2000  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Sylvain Pion

#ifndef CGAL_IO_TRIANGULATION_GEOMVIEW_OSTREAM_2_H
#define CGAL_IO_TRIANGULATION_GEOMVIEW_OSTREAM_2_H

#include <CGAL/license/Triangulation_2.h>


#include <CGAL/IO/Geomview_stream.h>
#include <CGAL/Triangulation_2.h>

namespace CGAL {

// There are 2 drawing functions for triangulations : depending on the wired
// mode of the Geomview_stream, we draw either the edges or the faces.

// TODO :
// - Check the correctness when dimension < 2.
// - Use the current stream color instead of built-in constant.

template < class GT, class TDS >
void
show_triangulation_edges(Geomview_stream &gv, const Triangulation_2<GT,TDS> &T)
{
    // Header.
    gv.set_ascii_mode();
    gv << "(geometry " << gv.get_new_id("triangulationedge")
       << " {appearance {}{ SKEL \n"
       << T.number_of_vertices()
       << T.number_of_vertices() + T.number_of_faces()-1 << "\n";

    // Finite vertices coordinates.
    std::map<typename Triangulation_2<GT, TDS>::Vertex_handle, int> V;
    int inum = 0;
    for( typename Triangulation_2<GT, TDS>::Vertex_iterator
	  vit = T.vertices_begin(); vit != T.vertices_end(); ++vit) {
        V[vit] = inum++;
        gv << vit->point() << "\n";
    }
  
    // Finite edges indices.
    for( typename Triangulation_2<GT, TDS>::Edge_iterator
	  eit = T.edges_begin(); eit != T.edges_end(); ++eit) {
        gv << 2
           << V[(*eit).first->vertex(T.ccw((*eit).second))]
           << V[(*eit).first->vertex(T. cw((*eit).second))]
           << "\n"; // without color.
        // << 4 << drand48() << drand48() << drand48() << 1.0; // random color
    }
}

template < class GT, class TDS >
void
show_triangulation_faces(Geomview_stream &gv, const Triangulation_2<GT,TDS> &T)
{
    // Header.
    gv.set_binary_mode();
    gv << "(geometry " << gv.get_new_id("triangulation")
       << " {appearance {}{ OFF BINARY\n"
       << T.number_of_vertices() << T.number_of_faces() << 0;

    // Finite vertices coordinates.
    std::map<typename Triangulation_2<GT, TDS>::Vertex_handle, int> V;
    int inum = 0;
    for( typename Triangulation_2<GT, TDS>::Vertex_iterator
	  vit = T.vertices_begin(); vit != T.vertices_end(); ++vit) {
        V[vit] = inum++;
        gv << vit->point();
    }
  
    // Finite faces indices.
    for( typename Triangulation_2<GT, TDS>::Face_iterator
	  fit = T.faces_begin(); fit != T.faces_end(); ++fit) {
        gv << 3;
        for (int i=0; i<3; i++)
            gv << V[fit->vertex(i)];
        gv << 0; // without color.
     // gv << 4 << drand48() << drand48() << drand48() << 1.0; // random color
    }
}

template < class GT, class TDS >
Geomview_stream&
operator<<( Geomview_stream &gv, const Triangulation_2<GT,TDS> &T)
{
    bool ascii_bak = gv.get_ascii_mode();
    bool raw_bak = gv.set_raw(true);

    if (gv.get_wired())
	show_triangulation_edges(gv, T);
    else
	show_triangulation_faces(gv, T);

    // Footer.
    gv << "}})";

    gv.set_raw(raw_bak);
    gv.set_ascii_mode(ascii_bak);
    return gv;
}

} //namespace CGAL

#endif // CGAL_IO_TRIANGULATION_GEOMVIEW_OSTREAM_2_H
