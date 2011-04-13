// ======================================================================
//
// Copyright (c) 2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 
//
// file          : include/CGAL/IO/Triangulation_geomview_ostream_3.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Sylvain Pion <Sylvain.Pion@sophia.inria.fr>
//
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ======================================================================

#ifndef CGAL_IO_TRIANGULATION_GEOMVIEW_OSTREAM_3_H
#define CGAL_IO_TRIANGULATION_GEOMVIEW_OSTREAM_3_H

#include <CGAL/IO/Geomview_stream.h>
#include <CGAL/Triangulation_3.h>

// TODO :
// - Check the correctness when dimension < 3.
// - Use the stream color instead of built-in constant/random.
// - If interfaces were more similar, we could think of sharing 2d and 3d ?

CGAL_BEGIN_NAMESPACE

// This one is to show the edges of a 3D triangulation.
template < class GT, class TDS >
void
show_triangulation_edges(Geomview_stream &gv, const Triangulation_3<GT,TDS> &T)
{
  // Header.
  gv.set_ascii_mode();
  gv << "(geometry " << gv.get_new_id("triangulationedge")
     << " {appearance {}{ SKEL \n"
     << T.number_of_vertices() << T.number_of_finite_edges() << "\n";

  // Finite vertices coordinates.
  std::map<Triangulation_3<GT, TDS>::Vertex_handle, int> V;
  int inum = 0;
  for( Triangulation_3<GT, TDS>::Vertex_iterator
	  vit = T.finite_vertices_begin(); vit != T.vertices_end(); ++vit) {
    V[vit] = inum++;
    gv << vit->point() << "\n";
  }
  
  // Finite edges indices.
  for( Triangulation_3<GT, TDS>::Edge_iterator
	  eit = T.finite_edges_begin(); eit != T.edges_end(); ++eit) {
      gv << 2
         << V[(*eit).first->vertex((*eit).second)]
         << V[(*eit).first->vertex((*eit).third)]
         << "\n"; // without color.
      // << 4 << drand48() << drand48() << drand48() << 1.0; // random color.
  }
}

// This one outputs the facets.
template < class GT, class TDS >
void
show_triangulation_faces(Geomview_stream &gv, const Triangulation_3<GT,TDS> &T)
{
  // Header.
  gv.set_binary_mode();
  gv << "(geometry " << gv.get_new_id("triangulation")
     << " {appearance {}{ OFF BINARY\n"
     << T.number_of_vertices() << T.number_of_finite_facets() << 0;

  // Finite vertices coordinates.
  std::map<Triangulation_3<GT, TDS>::Vertex_handle, int> V;
  int inum = 0;
  for( Triangulation_3<GT, TDS>::Vertex_iterator
	  vit = T.finite_vertices_begin(); vit != T.vertices_end(); ++vit) {
    V[vit] = inum++;
    gv << vit->point();
  }
  
  // Finite facets indices.
  for( Triangulation_3<GT, TDS>::Facet_iterator
	  fit = T.finite_facets_begin(); fit != T.facets_end(); ++fit) {
      gv << 3;
      for (int i=0; i<4; i++)
          if (i != (*fit).second)
	      gv << V[(*fit).first->vertex(i)];
      gv << 0; // without color.
      // gv << 4 << drand48() << drand48() << drand48() << 1.0; // random color
  }
}

template < class GT, class TDS >
Geomview_stream&
operator<<( Geomview_stream &gv, const Triangulation_3<GT,TDS> &T)
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

CGAL_END_NAMESPACE

#endif // CGAL_IO_TRIANGULATION_GEOMVIEW_OSTREAM_3_H
