// Copyright (c) 1998  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Rebufat Francois (Francois.Rebufat@sophia.inria.fr)

#include <cassert>

template <class Tds>
void
_test_cell_tds_3(const Tds &)
{
  typedef typename Tds::Vertex_handle            Vertex_handle;
  typedef typename Tds::Cell_handle              Cell_handle;

  std::cout << "   Cells Tds Constructors " << std::endl;
  int ind;
  Tds tds;
  tds.set_dimension(3);
  Vertex_handle v0= tds.create_vertex();
  Vertex_handle v1= tds.create_vertex();
  Vertex_handle v2= tds.create_vertex();
  Vertex_handle v3= tds.create_vertex();
  Cell_handle c1 = tds.create_cell(v0, v1, v2, v3);
  assert(tds.has_vertex(c1,v0));
  assert(tds.has_vertex(c1,v1));
  assert(tds.has_vertex(c1,v2));
  assert(tds.has_vertex(c1,v3));
  Cell_handle n0=tds.create_cell();
  Cell_handle n1=tds.create_cell();
  Cell_handle n2=tds.create_cell();
  Cell_handle n3=tds.create_cell();
  Cell_handle c2 = tds.create_cell(v0, v1, v2, v3, n0, n1, n2, n3);

  std::cout << "   Access cell's functions " << std::endl;
  assert(tds.has_vertex(c2, v0));
     ind=tds.index(c2,v0);
     assert(tds.vertex(c2,ind)==v0);
     assert(tds.has_vertex(c2,v0, ind));

     assert(tds.has_vertex(c2,v1));
     ind=tds.index(c2,v1);
     assert(tds.vertex(c2,ind)==v1);
     assert(tds.has_vertex(c2,v1, ind));

     assert(tds.has_vertex(c2,v2));
     ind=tds.index(c2,v2);
     assert(tds.vertex(c2,ind)==v2);
     assert(tds.has_vertex(c2,v2, ind));

     assert(tds.has_vertex(c2,v3));
     ind=tds.index(c2,v3);
     assert(tds.vertex(c2,ind)==v3);
     assert(tds.has_vertex(c2,v3, ind));

     assert(tds.has_neighbor(c2,n0));
     ind=tds.index(c2,n0);
     assert(tds.neighbor(c2,ind)==n0);
     assert(tds.has_neighbor(c2,n0,ind));

     assert(tds.has_neighbor(c2,n1));
     ind=tds.index(c2,n1);
     assert(tds.neighbor(c2,ind)==n1);
     assert(tds.has_neighbor(c2,n1,ind));

     assert(tds.has_neighbor(c2,n2));
     ind=tds.index(c2,n2);
     assert(tds.neighbor(c2,ind)==n2);
     assert(tds.has_neighbor(c2,n2,ind));

     assert(tds.has_neighbor(c2,n3));
     ind=tds.index(c2,n3);
     assert(tds.neighbor(c2,ind)==n3);
     assert(tds.has_neighbor(c2,n3,ind));

     std::cout << "   setting cell's functions " << std::endl;
<<<<<<< HEAD
   tds.set_vertex(c2,0,v1);
   tds.set_vertex(c2,1,v2);
   tds.set_vertex(c2,2,v3);
   tds.set_vertex(c2,3,v0);
   assert(tds.index(c2,v0)==3);
   assert(tds.index(c2,v1)==0);
   assert(tds.index(c2,v2)==1);
   assert(tds.index(c2,v3)==2);
   //   tds.set_vertices();
   //   assert(tds.vertex(0)==NULL);
   //   assert(tds.vertex(1)==NULL);
   //   assert(tds.vertex(2)==NULL);
   //   assert(tds.vertex(3)==NULL);
   tds.set_vertices(c2,v0, v1, v2, v3);
   assert(tds.index(c2,v0)==0);
   assert(tds.index(c2,v1)==1);
   assert(tds.index(c2,v2)==2);
   assert(tds.index(c2,v3)==3);

   tds.set_neighbor(c2,0,n1);
   tds.set_neighbor(c2,1,n2);
   tds.set_neighbor(c2,2,n3);
   tds.set_neighbor(c2,3,n0);
   assert(tds.index(c2,n0)==3);
   assert(tds.index(c2,n1)==0);
   assert(tds.index(c2,n2)==1);
   assert(tds.index(c2,n3)==2);
   //   tds.set_neighbors();
   //   assert(tds.neighbor(0)==NULL);
   //   assert(tds.neighbor(1)==NULL);
   //   assert(tds.neighbor(2)==NULL);
   //   assert(tds.neighbor(3)==NULL);
   tds.set_neighbors(c2,n0, n1, n2, n3);
   assert(tds.index(c2,n0)==0);
   assert(tds.index(c2,n1)==1);
   assert(tds.index(c2,n2)==2);
   assert(tds.index(c2,n3)==3);
=======
   c2->set_vertex(0,v1);
   c2->set_vertex(1,v2);
   c2->set_vertex(2,v3);
   c2->set_vertex(3,v0);
   assert(c2->index(v0)==3);
   assert(c2->index(v1)==0);
   assert(c2->index(v2)==1);
   assert(c2->index(v3)==2);
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
   c2->set_neighbors(n0, n1, n2, n3);
   assert(c2->index(n0)==0);
   assert(c2->index(n1)==1);
   assert(c2->index(n2)==2);
   assert(c2->index(n3)==3);
>>>>>>> cgal/master

   std::cout << "   Tds Destructors " << std::endl;
   tds.delete_vertex(v0);
   tds.delete_vertex(v1);
   tds.delete_vertex(v2);
   tds.delete_vertex(v3);
}
