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
  assert(c1->has_vertex(v0));
  assert(c1->has_vertex(v1));
  assert(c1->has_vertex(v2));
  assert(c1->has_vertex(v3));
  Cell_handle n0=tds.create_cell();
  Cell_handle n1=tds.create_cell();
  Cell_handle n2=tds.create_cell();
  Cell_handle n3=tds.create_cell();
  Cell_handle c2 = tds.create_cell(v0, v1, v2, v3, n0, n1, n2, n3);

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

   std::cout << "   Tds Destructors " << std::endl;
   tds.delete_vertex(v0);
   tds.delete_vertex(v1);
   tds.delete_vertex(v2);
   tds.delete_vertex(v3);
}
