// ============================================================================
//
// Copyright (c) 1998-1999 The CGAL Consortium
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
// file          : example/Triangulation3/example_tds.C
// revision      : $Revision$
// author(s)     : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//
// coordinator   : INRIA Sophia Antipolis (Mariette Yvinec)
//
// ============================================================================

#include <CGAL/basic.h>

#include <iostream>
#include <fstream>
#include <cassert>
#include <vector>

#include <CGAL/Triangulation_cell_base_3.h>
#include <CGAL/Triangulation_vertex_base_3.h>
#include <CGAL/Triangulation_data_structure_3.h>

// We define a minimal traits class, because the Point_3 type is needed in
// order to instantiate Triangulation_vertex_base_3<>.

class empty_traits {
public:
  class Point_3 {};
};
typedef empty_traits Gt;

typedef CGAL::Triangulation_vertex_base_3<Gt>       Vb;
typedef CGAL::Triangulation_cell_base_3<Gt>         Cb;

typedef CGAL::Triangulation_data_structure_3<Vb,Cb> Tds;

typedef Tds::Cell                                   TDSCell;
typedef Tds::Vertex                                 TDSVertex;

int main()
{
  Tds T;

  assert( T.number_of_vertices() == 0 );
  assert( T.dimension() == -2 );
  assert( T.is_valid() );

  std::vector<TDSVertex*> PV(7);

  PV[0] = T.insert_increase_dimension(NULL);
  assert( T.number_of_vertices() == 1 );
  assert( T.dimension() == -1 );
  assert( T.is_valid() );

  int i;
  // each of the following insertions of vertices increases the dimension
  for ( i=1; i<5; i++ ) {
    PV[i] = T.insert_increase_dimension(NULL, PV[0]);
    assert( T.number_of_vertices() == i+1 );
    assert( T.dimension() == i-1 );
    assert( T.is_valid() );
  }
  assert( T.number_of_cells() == 5 );

  // we now have a simplex in dimension 4

  // cell incident to PV[0]
  TDSCell* c = PV[0]->cell();
  int ind;
  assert( c->has_vertex( PV[0], ind ) );
  // PV[0] is the vertex of index ind in c

  // insertion of a new vertex in the facet opposite to PV[0]
  PV[5] = T.insert_in_facet(NULL, c, ind);
  
  assert( T.number_of_vertices() == 6 );
  assert( T.dimension() == 3 );
  assert( T.is_valid() );

  // insertion of a new vertex in c
  PV[6] = T.insert_in_cell(NULL, c);

  assert( T.number_of_vertices() == 7 );
  assert( T.dimension() == 3 );
  assert( T.is_valid() );

  std::ofstream oFileT("output_tds",std::ios::out);
  // writing file output_tds; 
  oFileT << T; 

  return 0;
}
