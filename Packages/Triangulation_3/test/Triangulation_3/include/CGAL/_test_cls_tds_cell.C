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
// file          : include/CGAL/_test_cls_tds_cell.C
// revision      : 
// revision_date : 
// author(s)     : Rebufat Francois (Francois.Rebufat@sophia.inria.fr)
//
// coordinator   : INRIA Sophia-Antipolis
// ============================================================================

#include <cassert>

template <class Cell>
void
_test_cell_tds_3( const Cell &)
{
  typedef typename Cell::Vertex            Vertex;
  typedef typename Cell::Tds               Tds;

  std::cout << "   Cells Tds Constructors " << std::endl;
  int ind;
  Tds tds;
  Vertex* v0= tds.create_vertex();
  Vertex* v1= tds.create_vertex();
  Vertex* v2= tds.create_vertex();
  Vertex* v3= tds.create_vertex();
  Cell* c1 = tds.create_cell(v0, v1, v2, v3);
  assert(c1->has_vertex(v0));
  assert(c1->has_vertex(v1));
  assert(c1->has_vertex(v2));
  assert(c1->has_vertex(v3));
  Cell* n0=tds.create_cell();
  Cell* n1=tds.create_cell();
  Cell* n2=tds.create_cell();
  Cell* n3=tds.create_cell();
  Cell* c2 = tds.create_cell(v0, v1, v2, v3, n0, n1, n2, n3);

  std::cout << "   Access cell's functions " << std::endl;
     assert(c2->has_vertex(v0));
     ind=c2->index(v0);
     assert(c2->vertex(ind)==v0);
     assert(c2->has_vertex(v0, ind));

     assert(c2->has_vertex(v1));
     ind=c2->index(v1);
     assert(c2->vertex(ind)==v1);
     assert(c2->has_vertex(v1, ind));

     assert(c2->has_vertex(v2));
     ind=c2->index(v2);
     assert(c2->vertex(ind)==v2);
     assert(c2->has_vertex(v2, ind));

     assert(c2->has_vertex(v3));
     ind=c2->index(v3);
     assert(c2->vertex(ind)==v3);
     assert(c2->has_vertex(v3, ind));

     assert(c2->has_neighbor(n0));
     ind=c2->index(n0);
     assert(c2->neighbor(ind)==n0);
     assert(c2->has_neighbor(n0,ind));

     assert(c2->has_neighbor(n1));
     ind=c2->index(n1);
     assert(c2->neighbor(ind)==n1);
     assert(c2->has_neighbor(n1,ind));

     assert(c2->has_neighbor(n2));
     ind=c2->index(n2);
     assert(c2->neighbor(ind)==n2);
     assert(c2->has_neighbor(n2,ind));

     assert(c2->has_neighbor(n3));
     ind=c2->index(n3);
     assert(c2->neighbor(ind)==n3);
     assert(c2->has_neighbor(n3,ind));

     std::cout << "   setting cell's functions " << std::endl;
   c2->set_vertex(0,v1);
   c2->set_vertex(1,v2);
   c2->set_vertex(2,v3);
   c2->set_vertex(3,v0);
   assert(c2->index(v0)==3);
   assert(c2->index(v1)==0);
   assert(c2->index(v2)==1);
   assert(c2->index(v3)==2);
   //   c2->set_vertices();
   //   assert(c2->vertex(0)==NULL);
   //   assert(c2->vertex(1)==NULL);
   //   assert(c2->vertex(2)==NULL);
   //   assert(c2->vertex(3)==NULL);
   c2->set_vertices(v0, v1, v2, v3);
   assert(c2->index(v0)==0);
   assert(c2->index(v1)==1);
   assert(c2->index(v2)==2);
   assert(c2->index(v3)==3);

   c2->set_neighbor(0,n1);
   c2->set_neighbor(1,n2);
   c2->set_neighbor(2,n3);
   c2->set_neighbor(3,n0);
   assert(c2->index(n0)==3);
   assert(c2->index(n1)==0);
   assert(c2->index(n2)==1);
   assert(c2->index(n3)==2);
   //   c2->set_neighbors();
   //   assert(c2->neighbor(0)==NULL);
   //   assert(c2->neighbor(1)==NULL);
   //   assert(c2->neighbor(2)==NULL);
   //   assert(c2->neighbor(3)==NULL);
   c2->set_neighbors(n0, n1, n2, n3);
   assert(c2->index(n0)==0);
   assert(c2->index(n1)==1);
   assert(c2->index(n2)==2);
   assert(c2->index(n3)==3);

   std::cout << "   Tds Destructors " << std::endl;
   tds.delete_vertex(v0);
   tds.delete_vertex(v1);
   tds.delete_vertex(v2);
   tds.delete_vertex(v3);
}
