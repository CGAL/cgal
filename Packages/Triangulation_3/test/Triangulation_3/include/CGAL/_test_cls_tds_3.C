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
// file          : include/CGAL/_test_cls_tds_3.C
// revision      : 
// revision_date : 
// author(s)     : Francois Rebufat <Francois Rebufat@sophia.inria.fr>
//                 Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//
// coordinator   : INRIA Sophia-Antipolis
// ============================================================================

#include <cassert>
#include <iostream>
#include <fstream>

#include "_test_cls_tds_vertex.C"
#include "_test_cls_tds_cell.C"

template <class Tds>
void
_test_cls_tds_3( const Tds &)
{
  typedef typename Tds::Vertex            Vertex;
  typedef typename Tds::Cell              Cell;
  typedef typename Tds::Edge              Edge;
  typedef typename Tds::Facet             Facet;

  typedef typename Tds::Vertex_iterator   Vertex_iterator;
  typedef typename Tds::Facet_iterator    Facet_iterator;
  typedef typename Tds::Edge_iterator     Edge_iterator;
  typedef typename Tds::Cell_iterator     Cell_iterator;

  // test Vertex and cell :
  std::cout << "    Test Vertex " << std::endl;
  _test_vertex_tds_3(Vertex());

  std::cout << "    Test Cell " << std::endl;
  _test_cell_tds_3(Cell());

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
  std::ofstream oFileT("Test_tds_IO_3", std::ios::out);
  std::ifstream iFileT("Test_tds_IO_3", std::ios::in);
  oFileT << tds1;
  iFileT >> tdsfromfile;
  assert(tdsfromfile.is_valid());
  assert(tdsfromfile.dimension() == -2);
  assert(tdsfromfile.number_of_vertices() == 0);

  std::cout << "    copy" << std::endl;
  Vertex_iterator vit;
  tds2.insert_increase_dimension(NULL);
  assert( tds2.number_of_vertices() == 1 );
  Tds tds3(tds2);
  Vertex * v2 = tds3.create_vertex();

  vit=tds3.vertices_begin();
  tds3.insert_increase_dimension(v2,&*vit);
  std::cout << "ok" << std::endl;
  Tds tds4 = tds3;
  Vertex * v3 = tds4.create_vertex();
  vit=tds4.vertices_begin();
  tds4.insert_increase_dimension(v3,&*vit);
  std::cout << "ok" << std::endl;
  Tds tds5;
  tds5.swap(tds4);
  tds4=tds5;
  Vertex * v4 = tds5.create_vertex();
  vit=tds5.vertices_begin();
  tds5.insert_increase_dimension(v4,&*vit);
  std::cout << "ok" << std::endl;
  Tds tds6;
  tds6.swap(tds5);
  tds5=tds6;
  Vertex * v5 = tds6.create_vertex();
  vit=tds6.vertices_begin();
  tds6.insert_increase_dimension(v5,&*vit);
  std::cout << "ok" << std::endl;

  // Setting functions
  std::cout << "    setting functions" << std::endl;
  tds1.set_number_of_vertices(1);
  
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
  assert(tds1.number_of_vertices()==1);

  std::cout << "  Test flip " << std::endl;
  assert(tds6.is_valid());
  Cell_iterator cit, cdone;
  int nbflips=0;
  int i;
  Vertex * v6 = tds6.create_vertex();
  cit = tds6.cells_begin();
  tds6.insert_in_cell(v6, &(*cit));
  Vertex * v7 = tds6.create_vertex();
  cit = tds6.cells_begin();
  tds6.insert_in_cell(v7, &(*cit));
  Vertex * v8 = tds6.create_vertex();
  cit = tds6.cells_begin();
  tds6.insert_in_cell(v8, &(*cit));
  assert(tds6.number_of_vertices()==8);
//   std::cout << tds6.number_of_cells()<< " cells" << std::endl;

  cdone = tds6.cells_end();
  
  std::set< Vertex* > set_of_vertices;
  
  for ( cit = tds6.cells_begin(); cit != cdone; cit++ ) {
    // NOTE : the triangulation is modified during loop
    // --> the cell_iterator does not mean a lot
    for ( i=0; i<4; i++ ) {
      std::set< Vertex* > set_of_vertices;
      tds6.incident_vertices( (&(*cit))->vertex(i), set_of_vertices );
      if ( set_of_vertices.find
	   ( (&(*cit))->neighbor(i)->vertex
	     ( (&(*cit))->neighbor(i)->index( &(*cit) ) ) 
	     ) 
	   == set_of_vertices.end() ) {
	nbflips++;
	tds6.flip_flippable( &(*cit), i );
	assert(tds6.is_valid());
// 	if ( tds6.flip( &(*cit), i ) ) {
// 	  tds6.is_valid(true);
// 	  nbflips++;
// 	}
      }
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
	flipped = tds6.flip( &(*cit), i, j ) ;
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
