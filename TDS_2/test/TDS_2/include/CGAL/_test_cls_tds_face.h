// ============================================================================
//
// Copyright (c) 1998 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       :
// release_date  :
//
// source        :
// file          : include/CGAL/_test_cls_tds_face.h
// revision      :
// revision_date :
// author(s)     : Herve Bronnimann (Herve.Bronnimann@sophia.inria.fr)
//
// coordinator   : INRIA Sophia-Antipolis
// ============================================================================


#include <cassert>

namespace CGAL {

template <class Tds>
void
_test_cls_tds_face( const Tds )
{
  std::cout << "Face Tds constructors" << std::endl;

  typedef typename Tds::Edge                 Edge;
  typedef typename Tds::Face_handle          Face_handle;
  typedef typename Tds::Vertex_handle        Vertex_handle;

  Tds tds;

  // Build a few objects
  int i;
  Vertex_handle vh1 = tds.create_vertex();
  Vertex_handle vh2 = tds.create_vertex();
  Vertex_handle vh3 = tds.create_vertex();
  Vertex_handle vh4 = tds.create_vertex();

  // Test constructors
  Face_handle fh1 = tds.create_face();
  Face_handle fh2 = tds.create_face();
  Face_handle fh3 = tds.create_face();
  Face_handle gh4 = tds.create_face(vh1,vh2,vh3);
  Face_handle gh5 = tds.create_face(vh1,vh2,vh3,fh1,fh2,fh3);
  Face_handle fh4 = tds.create_face();

  // Test vertex()
  assert( gh4->vertex(0) == vh1 );
  assert( gh4->vertex(1) == vh2 );
  assert( gh4->vertex(2) == vh3 );
  assert( gh5->vertex(0) == vh1 );
  assert( gh5->vertex(1) == vh2 );
  assert( gh5->vertex(2) == vh3 );

  // Test index(Vertex *)
  assert( gh4->index(vh1) == 0 );
  assert( gh4->index(vh2) == 1 );
  assert( gh4->index(vh3) == 2 );
  assert( gh5->index(vh1) == 0 );
  assert( gh5->index(vh2) == 1 );
  assert( gh5->index(vh3) == 2 );

  // Test has_vertex()
  assert( gh4->has_vertex(vh1) );
  assert( ! gh4->has_vertex(vh4) );
  assert( gh4->has_vertex(vh2,i) && (i==1) );
  assert( gh5->has_vertex(vh3,i) && (i==2) );

  // Test set_vertex()
  fh1->set_vertex(0,vh4);
  fh1->set_vertex(1,vh3);
  fh1->set_vertex(2,vh2);
  fh2->set_vertex(0,vh3);
  fh2->set_vertex(1,vh4);
  fh2->set_vertex(2,vh1);
  fh3->set_vertex(0,vh2);
  fh3->set_vertex(1,vh1);
  fh3->set_vertex(2,vh4);
  assert( fh1->vertex(0) == vh4 );
  assert( fh1->vertex(1) == vh3 );
  assert( fh1->vertex(2) == vh2 );
  assert( fh2->vertex(0) == vh3 );
  assert( fh2->vertex(1) == vh4 );
  assert( fh2->vertex(2) == vh1 );
  assert( fh3->vertex(0) == vh2 );
  assert( fh3->vertex(1) == vh1 );
  assert( fh3->vertex(2) == vh4 );

  // Test set_vertices()
  gh5->set_vertices();
  assert( gh5->vertex(0) == Vertex_handle());
  assert( gh5->vertex(1) == Vertex_handle());
  assert( gh5->vertex(2) == Vertex_handle());
  gh5->set_vertices( vh1, vh2, vh3 );
  assert( gh5->vertex(0) == vh1 );
  assert( gh5->vertex(1) == vh2 );
  assert( gh5->vertex(2) == vh3 );

  // Test neighbor()
  assert( gh5->neighbor(0) == fh1 );
  assert( gh5->neighbor(1) == fh2 );
  assert( gh5->neighbor(2) == fh3 );

  // Test index(Face *)
  assert( gh5->index(fh1) == 0 );
  assert( gh5->index(fh2) == 1 );
  assert( gh5->index(fh3) == 2 );

  // Test has_neighbor()
  assert( gh5->has_neighbor(fh1) );
  assert( ! gh5->has_neighbor(fh4) );
  assert( gh5->has_neighbor(fh2,i) && (i==1) );
  assert( gh5->has_neighbor(fh3,i) && (i==2) );

  // Test set_neighbor()
  fh1->set_neighbor(0,gh4);
  fh2->set_neighbor(1,gh4);
  fh3->set_neighbor(2,gh4);
  assert( fh1->neighbor(0) == gh4 );
  assert( fh2->neighbor(1) == gh4 );
  assert( fh3->neighbor(2) == gh4 );

  // Test set_neighbors()
  gh5->set_neighbors();
  assert( gh5->neighbor(0) == Face_handle());
  assert( gh5->neighbor(1) == Face_handle());
  assert( gh5->neighbor(2) == Face_handle());
  gh5->set_neighbors( fh1, fh2, fh3 );
  assert( gh5->neighbor(0) == fh1 );
  assert( gh5->neighbor(1) == fh2 );
  assert( gh5->neighbor(2) == fh3 );
  gh4->set_neighbors( fh1, fh2, fh3 );
  assert( gh4->neighbor(0) == fh1 );
  assert( gh4->neighbor(1) == fh2 );
  assert( gh4->neighbor(2) == fh3 );

  //Test mirror_vertex() mirror_index()
  assert( tds.mirror_vertex(gh4, 0) == vh4);
  assert( tds.mirror_vertex(fh1, 0) == vh1);
  assert( tds.mirror_vertex(fh2, 1) == vh2);
  assert( tds.mirror_vertex(fh3, 2) == vh3);
  assert( tds.mirror_index(gh4, 0) == 0);
  assert( tds.mirror_edge(Edge(gh4, 0)) == Edge(fh1, 0));
  assert( tds.mirror_index(fh1, 0) == 0);
  assert( tds.mirror_index(fh2, 1) == 1);
  assert( tds.mirror_index(fh3, 2) == 2);

  // Test is_valid
  assert( gh4->is_valid() );

  // Test ccw() and cw()
  assert( fh1->ccw(0) == 1 );
  assert( fh1->ccw(1) == 2 );
  assert( fh1->ccw(2) == 0 );
  assert( fh1->cw(0) == 2 );
  assert( fh1->cw(1) == 0 );
  assert( fh1->cw(2) == 1 );

  // Test dimension
  assert(gh4->dimension() == 2);

  //Test low dimensional faces
  Face_handle gh1 = tds.create_face(vh2,vh3, Vertex_handle());
  Face_handle gh2 = tds.create_face(vh3,vh1, Vertex_handle());
  Face_handle gh3 = tds.create_face(vh1,vh2, Vertex_handle());
  gh1->set_neighbors( gh2, gh3, Face_handle());
  gh2->set_neighbors( gh3, gh1, Face_handle());
  gh3->set_neighbors( gh1, gh2, Face_handle());

  assert(gh1->dimension() == 1);
  assert (gh1->is_valid());

  Face_handle hh1=tds.create_face(vh1, Vertex_handle(), Vertex_handle());
  Face_handle hh2=tds.create_face(vh2, Vertex_handle(), Vertex_handle(),
                                  hh1, Face_handle(), Face_handle());
  hh1->set_neighbor(0, hh2);
  assert (hh1->dimension() == 0);
  assert (hh1->is_valid());

  assert(tds.mirror_vertex(gh1,0) == vh1);
  assert(tds.mirror_vertex(gh1,1) == vh1);
  assert(tds.mirror_index(gh1,0) == 1);
  assert(tds.mirror_index(gh1,1) == 0);

  return;
}

} //namespace CGAL
