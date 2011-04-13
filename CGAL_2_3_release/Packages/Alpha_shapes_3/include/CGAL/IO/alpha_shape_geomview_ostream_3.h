// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-2.0-I-20 $
// release_date  : $CGAL_Date: 1999/06/02 $
//
// file          : include/CGAL/IO/alpha_shape_geomview_ostream_3.h
// package       : Alpha_shapes_3(1.0)
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Tran Kai Frank DA <Frank.Da@sophia.inria.fr>
//
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ======================================================================

#ifndef CGAL_IO_ALPHA_SHAPE_GEOMVIEW_OSTREAM_3_H
#define CGAL_IO_ALPHA_SHAPE_GEOMVIEW_OSTREAM_3_H

#include <CGAL/IO/Geomview_stream.h>
#include <CGAL/Alpha_shape_3.h>

// TODO :
// - Check the correctness when dimension < 3.
// - Use the stream color instead of built-in constant/random.
// - If interfaces were more similar, we could think of sharing 2d and 3d ?

//-------------------------------------------------------------------
CGAL_BEGIN_NAMESPACE
//-------------------------------------------------------------------

// This one is to show the edges of a 3D triangulation.
// template < class GT, class TDS >
// void
// show_triangulation_edges(Geomview_stream &gv, 
//                          const Alpha_shape_3<GT,TDS> &T)
// {
//   // Header.
//   gv.set_ascii_mode();
//   gv << "(geometry " << gv.get_new_id("triangulationedge")
//      << " {appearance {}{ SKEL \n"
//      << T.number_of_vertices() << T.number_of_finite_edges() << "\n";

//   // Finite vertices coordinates.
//   std::map<Alpha_shape_3<GT, TDS>::Vertex_handle, int> V;
//   int inum = 0;
//   for( Alpha_shape_3<GT, TDS>::Vertex_iterator
// 	  vit = T.finite_vertices_begin(); vit != T.vertices_end(); ++vit) {
//     V[vit] = inum++;
//     gv << vit->point() << "\n";
//   }
  
//   // Finite edges indices.
//   for( Alpha_shape_3<GT, TDS>::Edge_iterator
// 	  eit = T.finite_edges_begin(); eit != T.edges_end(); ++eit) {
//       gv << 2
//          << V[(*eit).first->vertex((*eit).second)]
//          << V[(*eit).first->vertex((*eit).third)]
//          << "\n"; // without color.
//       // << 4 << drand48() << drand48() << drand48() << 1.0; // random color
//   }
// }


//-------------------------------------------------------------------
// This one outputs the facets.
template < class Dt >
void
Alpha_shape_3<Dt>::show_alpha_shape_faces(Geomview_stream &gv)
{
  // Finite vertices coordinates.
  typename Alpha_shape_3<Dt>::Alpha_shape_vertices_iterator Vlist_it,
    Vlist_begin = Alpha_shape_vertices_begin(),
    Vlist_end = Alpha_shape_vertices_end();

  std::map<Alpha_shape_3<Dt>::Vertex_handle, int> V;
  int number_of_vertex = 0;
  for( Vlist_it = Vlist_begin; Vlist_it != Vlist_end; Vlist_it++) {
    V[*Vlist_it] = number_of_vertex++;
  }


  typename Alpha_shape_3<Dt>::Alpha_shape_facets_iterator Flist_it,
    Flist_begin = Alpha_shape_facets_begin(),
    Flist_end = Alpha_shape_facets_end();

  std::map<Alpha_shape_3<Dt>::Facet, int> F;
  int number_of_facets = 0;
  for( Flist_it = Flist_begin; Flist_it != Flist_end; Flist_it++) {
    F[*Flist_it] = number_of_facets++;
  }

  // Header.
  gv.set_binary_mode();
  gv << "(geometry " << gv.get_new_id("alpha_shape")
     << " {appearance {}{ OFF BINARY\n"
     << number_of_vertex << number_of_facets << 0;

  for( Vlist_it = Vlist_begin; Vlist_it != Vlist_end; Vlist_it++) {
    gv << (*Vlist_it)->point();
  }
  
  // Finite facets indices.
  for( Flist_it = Flist_begin; Flist_it != Flist_end; Flist_it++){
      gv << 3;
      for (int i=0; i<4; i++)
          if (i != (*Flist_it).second)
	      gv << V[(*Flist_it).first->vertex(i)];
      gv << 0; // without color.
      // gv << 4 << drand48() << drand48() << drand48() << 1.0; // random color
  }
}

//-------------------------------------------------------------------

template < class Dt >
Geomview_stream&
operator<<( Geomview_stream &gv, Alpha_shape_3<Dt>& A)
{
    bool ascii_bak = gv.get_ascii_mode();
    bool raw_bak = gv.set_raw(true);

//     if (gv.get_wired())
//         show_alpha_shape_edges(gv, T);
//     else
        A.show_alpha_shape_faces(gv);

    // Footer.
    gv << "}})";

    gv.set_raw(raw_bak);
    gv.set_ascii_mode(ascii_bak);
    return gv;
}

//-------------------------------------------------------------------
CGAL_END_NAMESPACE
//-------------------------------------------------------------------

#endif // CGAL_IO_ALPHA_SHAPE_GEOMVIEW_OSTREAM_3_H
