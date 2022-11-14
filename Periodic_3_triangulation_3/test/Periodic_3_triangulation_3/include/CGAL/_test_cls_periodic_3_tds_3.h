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
// Author(s)     : Francois Rebufat
//                 Monique Teillaud <Monique.Teillaud@sophia.inria.fr>

#include <cassert>
#include <iostream>
#include <fstream>

#include <CGAL/use.h>
#include "_test_cls_periodic_3_tds_vertex.h"
#include "_test_cls_periodic_3_tds_cell.h"

template <class Tds>
void
_test_cls_periodic_3_tds_3( const Tds &)
{
  typedef typename Tds::Vertex_range      Vertex_range;
  typedef typename Tds::Cell_range        Cell_range;

  typedef typename Tds::Vertex            Vertex;
  typedef typename Tds::Cell              Cell;
  CGAL_USE_TYPE(Cell);
  typedef typename Tds::Edge              Edge;
  CGAL_USE_TYPE(typename Tds::Facet);

  typedef typename Tds::Vertex_handle     Vertex_handle;
  typedef typename Tds::Vertex_iterator   Vertex_iterator;
  CGAL_USE_TYPE(typename Tds::Facet_iterator);
  CGAL_USE_TYPE(typename Tds::Edge_iterator);
  typedef typename Tds::Cell_handle       Cell_handle;
  typedef typename Tds::Cell_iterator     Cell_iterator;

  // test Vertex and cell :
  std::cout << "    Test Vertex " << std::endl;
  _test_vertex_tds_3(Vertex());

  std::cout << "    Test Cell " << std::endl;
  _test_cell_tds_3(Tds());

  std::cout << "   Testing TDS " << std::endl;

  // Test constructors
  std::cout << "    constructors" << std::endl;
  Tds tds1;
  Tds tds2;

  // Test I/O for dimension -2
  // the other dimensions are not tested here
  // (they are implicitely tested in triangulation)
  Tds tdsfromfile;
  std::cout << "    I/O" << std::endl;
  {
    std::ofstream oFileT("Test_tds_IO_3", std::ios::out);
    oFileT << tds1;
  }
  std::ifstream iFileT("Test_tds_IO_3", std::ios::in);
  iFileT >> tdsfromfile;
  assert(tdsfromfile.is_valid());
  assert(tdsfromfile.dimension() == -2);
  assert(tdsfromfile.number_of_vertices() == 0);

  std::cout << "    copy" << std::endl;
  tds2.insert_increase_dimension();
  assert( tds2.number_of_vertices() == 1 );
  Tds tds3(tds2);

  Vertex_iterator vit;
  vit=tds3.vertices_begin();
  tds3.insert_increase_dimension(vit);
  std::cout << "ok" << std::endl;
  assert(tds3.is_valid());
  Tds tds4 = tds3;
  vit=tds4.vertices_begin();
  tds4.insert_increase_dimension(vit);
  std::cout << "ok" << std::endl;
  assert(tds4.is_valid());
  Tds tds5;
  tds5.swap(tds4);
  tds4=tds5;
  vit=tds5.vertices_begin();
  tds5.insert_increase_dimension(vit);
  std::cout << "ok" << std::endl;
  assert(tds5.is_valid());
  Tds tds6;
  tds6.swap(tds5);
  tds5=tds6;
  vit=tds6.vertices_begin();
  tds6.insert_increase_dimension(vit);
  std::cout << "ok" << std::endl;
  assert(tds6.is_valid());

  std::cout << "  Insert are tested in test_triangulation_3  " << std::endl;

  std::cout << "  Iterator and circulator are tested in test_triangulation_3  " << std::endl;

  // Access functions

  assert(tds1.dimension()==-2);
  assert(tds2.dimension()==-1);
  assert(tds3.dimension()==0);
  assert(tds4.dimension()==1);
  assert(tds5.dimension()==2);
  assert(tds6.dimension()==3);

  assert(tds3.number_of_vertices()==2);

  // Containers
  Vertex_range & vertex_c = tds3.vertices();
  Cell_range & cell_c = tds3.cells();

  assert(vertex_c.size() == 2);
  assert(cell_c.size() == 2);

  // Flips
  std::cout << "  Test flip " << std::endl;
  assert(tds6.is_valid());
  Cell_iterator cit, cdone;
  int nbflips=0;
  int i;
  cit = tds6.cells_begin();
  tds6.insert_in_cell(cit);
  cit = tds6.cells_begin();
  tds6.insert_in_cell(cit);
  cit = tds6.cells_begin();
  tds6.insert_in_cell(cit);
  assert(tds6.number_of_vertices()==8);
//   std::cout << tds6.number_of_cells()<< " cells" << std::endl;

  // We can't use the Cell_iterator while we modify the TDS.
  // However, we can store all Cell_handles beforehand,
  // since 2-3 flips do not affect the validity of existing cells.
  std::vector<Cell_handle> Cell_v;
  for (cit = tds6.cells_begin(); cit != tds6.cells_end(); ++cit)
      Cell_v.push_back(cit);

  for (typename std::vector<Cell_handle>::const_iterator ccit = Cell_v.begin();
       ccit != Cell_v.end(); ++ccit) {
    for ( i=0; i<4; i++ ) {
      assert(tds6.is_valid());
      assert(tds6.is_cell(*ccit));
      // old name (up to CGAL 3.4)
      // kept for backwards compatibility but not documented
      std::set< Vertex_handle > set_of_vertices_old;
      tds6.incident_vertices( (*ccit)->vertex(i),
                              std::inserter(set_of_vertices_old,
                                            set_of_vertices_old.begin() ) );
      if ( set_of_vertices_old.find(tds6.mirror_vertex(*ccit, i))
           == set_of_vertices_old.end() ) {
        nbflips++;
        tds6.flip_flippable( *ccit, i );
        assert(tds6.is_valid());
//         if ( tds6.flip( cit, i ) ) {
//           tds6.is_valid(true);
//           nbflips++;
//         }
      }
      // correct name
      std::set< Vertex_handle > set_of_vertices;
      tds6.adjacent_vertices( (*ccit)->vertex(i),
                              std::inserter(set_of_vertices,
                                            set_of_vertices.begin() ) );
      if ( set_of_vertices.find(tds6.mirror_vertex(*ccit, i))
           == set_of_vertices.end() ) {
        nbflips++;
        tds6.flip_flippable( *ccit, i );
        assert(tds6.is_valid());
//         if ( tds6.flip( cit, i ) ) {
//           tds6.is_valid(true);
//           nbflips++;
//         }
      }
    }
  }

  for (typename std::vector<Cell_handle>::const_iterator ccit = Cell_v.begin();
       ccit != Cell_v.end(); ++ccit) {
    for ( i=0; i<4; i++ ) {
       std::vector< Vertex_handle > vector_of_vertices_old;
     std::vector< Vertex_handle > vector_of_vertices;
      std::vector< Edge > vector_of_edges;

      // old name (up to CGAL 3.4)
      // kept for backwards compatibility but not documented
      tds6.incident_vertices
        ( (*ccit)->vertex(i), std::back_inserter(vector_of_vertices_old));
      // correct name
      tds6.adjacent_vertices
        ( (*ccit)->vertex(i), std::back_inserter(vector_of_vertices));

      tds6.incident_edges
        ( (*ccit)->vertex(i), std::back_inserter(vector_of_edges));

      assert(vector_of_edges.size() == vector_of_vertices_old.size());
      assert(vector_of_edges.size() == vector_of_vertices.size());
    }
  }


  std::cout << nbflips << " flips 2-3" << std::endl;
  assert(tds6.number_of_vertices()==8);
//  std::cout << tds6.number_of_cells()<< " cells" << std::endl;

  nbflips=0;
  bool flipped;
  int j;
  cit = tds6.cells_begin();
  cdone = tds6.cells_end();
  Cell_iterator next_cell;
  while ( cit != cdone ) {
    // NOTE : cells are deleted during loop
    // the cell_iterator is modified "by hand" (not using ++)
    flipped = false; i=0; j=1;
    next_cell = ++cit; --cit;
    while ( (! flipped) && (i<4) ) {
      if ( (i!=j) ) {
        // The Intel compiler has a bug and needs the explicit handle.
        Cell_handle ch = cit;
        flipped = tds6.flip( ch, i, j ) ;
        if (flipped) {
          nbflips++;
          assert(tds6.is_valid());
        }
      }
      if ( j==3 ) { i++; j=0; }
      else j++;
    }
    cit = next_cell;
  }
  std::cout << nbflips << " flips 3-2" << std::endl;
  assert(tds6.number_of_vertices()==8);

  // test destructor and return
  std::cout << "    test destructors and return" << std::endl;

  assert(tds1.is_valid());
  assert(tds2.is_valid());
  assert(tdsfromfile.is_valid());
  assert(tds3.is_valid());
  assert(tds4.is_valid());
  assert(tds5.is_valid());
  assert(tds6.is_valid());

//   tds1.clear();
//   tds2.clear();
//   tds3.clear();

}
