// Copyright (c) 1997  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
// Author(s)     : Lutz Kettner  <kettner@mpi-sb.mpg.de>

// Cuts out a piece of a halfedge data structure defined by a predicate.


#ifndef CGAL_HALFEDGEDS_CUT_COMPONENT_H
#define CGAL_HALFEDGEDS_CUT_COMPONENT_H 1

#include <CGAL/basic.h>
#include <CGAL/HalfedgeDS_decorator.h>

namespace CGAL {


template < class HDS, class Predicate >
typename HDS::Halfedge_handle 
halfedgeDS_cut_component( HDS&                           hds,
                          typename HDS::Halfedge_handle  h,
                          Predicate                      pred,
                          typename HDS::Halfedge_handle& cut_piece)
// Cuts out a piece of a halfedge data structure for which the predicate
// `pred' is true for the vertices.
// The edge-vertex graph of the piece has to be a connected component.
// The remaining piece gets a new boundary. Returns a border halfedge of 
// the new boundary on the remaining piece. Assigns a halfedge of the 
// cut outpiece to `cut_piece'.
// The geometry for the vertices
// on the boundary and the hole have to be taken care of after this
// function call. The cut-out piece can be deleted with the member
// function erase_connected_component of the decorator class. It can
// technically happen that only an isolated vertex would remain in the
// cut out piece, in which case a dummy halfedge pair and vertex will be
// created to keep this vertex representable in the halfedge data structure.
// Precondition: pred( h->vertex()) && ! pred( h->opposite()->vertex()).
{
    typedef typename HDS::Vertex_handle    Vertex_handle;
    typedef typename HDS::Halfedge_handle  Halfedge_handle;
    typedef typename HDS::Face_handle      Face_handle;
    typedef typename HDS::Vertex           Vertex;
    typedef typename HDS::Halfedge         Halfedge;
    typedef typename HDS::Face             Face;
    CGAL::HalfedgeDS_decorator<HDS> D(hds);

    CGAL_precondition( D.is_valid(false,3));
    CGAL_precondition( pred( h->vertex()));
    CGAL_precondition( ! pred( h->opposite()->vertex()));

    Halfedge_handle start = h;
    Halfedge_handle hnew;
    Halfedge_handle hlast;
    while (true) {
        // search re-entry point
        Halfedge_handle g = h;
        while ( pred( g->next()->vertex())) {
            g = g->next();
            // create border edges around cap
            D.set_face( g, Face_handle());
        }
        if ( hnew == Halfedge_handle()) {
            // first edge, special case
            CGAL_assertion( g->next() != h && g->next()->opposite() != h);
            Halfedge_handle gnext = g->next()->opposite();
            D.remove_tip( g);
            Vertex_handle v = D.vertices_push_back( Vertex());
            D.close_tip( gnext, v);
            hnew = hds.edges_push_back( Halfedge(), Halfedge());
            hlast = hnew->opposite();
            D.insert_tip( hlast, gnext);
            D.set_face( hnew, D.get_face( gnext));
            D.set_face_halfedge( hnew);
            h = g;
            D.set_vertex_halfedge( h);                
        } else { // general case and last case
            Halfedge_handle gnext = g->next()->opposite();
            if ( gnext == start && gnext == g) {
                // last edge, special case of isolated vertex.
                // Create dummy edge and dummy vertex and attach it to g
                g = hds.edges_push_back( Halfedge(), Halfedge());
                D.insert_tip( g, gnext);
                D.close_tip(g->opposite(), D.vertices_push_back(Vertex()));
                D.set_vertex_halfedge( g);                
                D.set_vertex_halfedge( g->opposite());                
            }
            D.remove_tip( g);
            Vertex_handle v = D.vertices_push_back( Vertex());
            D.close_tip( hnew, v);
            D.insert_tip( gnext, hnew);
            hnew = hds.edges_push_back( Halfedge(), Halfedge());
            D.insert_tip( hnew->opposite(), gnext);
            D.set_face( hnew, D.get_face( gnext));
            D.set_face_halfedge( hnew);
            h = g;
            D.set_vertex_halfedge( h);                
            if ( gnext == start) {
                // last edge, special
                D.insert_tip( hnew, hlast);
                break;
            }
        }
    } // while(true)
    Face_handle fnew = D.faces_push_back( Face());
    D.set_face_in_face_loop( hlast, fnew);
    D.set_face_halfedge( hlast);
    cut_piece = h;
    CGAL_postcondition( D.is_valid(false,3));
    return hlast;
}

template < class HDS, class Predicate >
typename HDS::Halfedge_handle 
halfedgeDS_cut_component( HDS&                           hds,
                          typename HDS::Halfedge_handle  h,
                          Predicate                      pred)
// Same function as above, but deletes the cut out piece immediately.
{
    typedef typename HDS::Halfedge_handle  Halfedge_handle;
    CGAL::HalfedgeDS_decorator<HDS> D(hds);

    CGAL_precondition( D.is_valid(false,3));
    CGAL_precondition( pred( h->vertex()));
    CGAL_precondition( ! pred( h->opposite()->vertex()));
    Halfedge_handle cut_piece;
    Halfedge_handle hnew = halfedgeDS_cut_component( hds, h, pred, cut_piece);
    D.erase_connected_component( cut_piece);
    CGAL_postcondition( D.is_valid(false,3));
    return hnew;
}   

} //namespace CGAL
#endif // CGAL_HALFEDGEDS_CUT_COMPONENT_H //
// EOF //
