// Copyright (c) 2005 Rijksuniversiteit Groningen (Netherlands)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
// Author(s)     : Nico Kruithof <Nico@cs.rug.nl>

#ifndef CGAL_TDS_INCREMENTAL_BUILDER_3_H
#define CGAL_TDS_INCREMENTAL_BUILDER_3_H 1

#include <CGAL/basic.h>
#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/array.h>
#include <set>
#include <list>
#include <vector>

namespace CGAL {

template < class Triangulation_3 >
class Triangulation_incremental_builder_3 {
public:
  typedef Triangulation_3                    T;
  typedef typename T::Vertex_handle          Vertex_handle;
  typedef typename T::Cell_handle            Cell_handle;
  typedef typename T::Facet                  Facet;

public:
  Triangulation_incremental_builder_3( T &t, bool verbose = false)
    : t(t), m_verbose(verbose) {
  }

  ~Triangulation_incremental_builder_3() {}

  void begin_triangulation(int dim) {
    t.clear();
    t.tds().delete_cell(t.infinite_vertex()->cell());
    // t.infinite = add_vertex();
    t.tds().set_dimension(dim); 
  }

  void end_triangulation() {
    construct_infinite_cells();
    CGAL_assertion(t.infinite_vertex()->cell() != Cell_handle());
  }

  Vertex_handle add_vertex();
  Cell_handle add_cell(Vertex_handle vh0, Vertex_handle vh1,
		       Vertex_handle vh2, Vertex_handle vh3);
	
private:
  void construct_infinite_cells();
  Cell_handle add_infinite_cell(Cell_handle ch, int i);
	
  void glue_cells(Cell_handle ch0, int ind0, Cell_handle ch1, int ind1);

  // Interior facets of the simplical cell:
  typedef std::pair < Vertex_handle, Vertex_handle > Vpair;
  typedef std::map < Vpair, Facet > MapPair;
  typedef typename MapPair::iterator   MapPairIt;
  typedef cpp11::array < Vertex_handle, 3 > Vtriple;
  typedef std::map < Vtriple, Facet > MapTriple;
  typedef typename MapTriple::iterator   MapTripleIt;
	
  Vtriple facet(Vertex_handle vh1, Vertex_handle vh2, Vertex_handle vh3) {
    if (vh1 < vh2) {
      if (vh2 < vh3) {
        return CGAL::make_array(vh1,vh2,vh3);
      } else {
        if (vh1 < vh3) {
          return CGAL::make_array(vh1,vh3,vh2);
        } else {
          return CGAL::make_array(vh3,vh1,vh2);
        }
      }
    }
    if (vh1 < vh3) {
      return CGAL::make_array(vh2,vh1,vh3);
    } else {
      if (vh2 < vh3) {
        return CGAL::make_array(vh2,vh3,vh1);
      } else {
        return CGAL::make_array(vh3,vh2,vh1);
      }
    }
  }
	
  MapTriple facets;

  T &t;
  bool m_verbose;
};



template < class TDS_>
typename Triangulation_incremental_builder_3< TDS_ >::Vertex_handle
Triangulation_incremental_builder_3< TDS_ >::add_vertex() {
  return t.tds().create_vertex();
}

template < class TDS_>
typename Triangulation_incremental_builder_3< TDS_ >::Cell_handle
Triangulation_incremental_builder_3< TDS_ >::add_cell(
  Vertex_handle vh0, Vertex_handle vh1, Vertex_handle vh2, Vertex_handle vh3) 
{
  CGAL_assertion(vh0 != NULL); CGAL_assertion(vh1 != NULL);
  CGAL_assertion(vh2 != NULL); CGAL_assertion(vh3 != NULL);
  CGAL_assertion(vh0 != vh1); CGAL_assertion(vh0 != vh2); CGAL_assertion(vh0 != vh3);
  CGAL_assertion(vh1 != vh2); CGAL_assertion(vh1 != vh3); CGAL_assertion(vh2 != vh3);
	
  Cell_handle ch =  t.tds().create_cell(vh0, vh1, vh2, vh3);
  // Neighbors are by default set to NULL
  vh0->set_cell(ch); vh1->set_cell(ch);
  vh2->set_cell(ch); vh3->set_cell(ch);

  MapTripleIt neighbIt;
  for (int i=0; i<4; i++) {
    Vtriple vtriple=facet(
		    ch->vertex((i+1)&3),
		    ch->vertex((i+2)&3),
		    ch->vertex((i+3)&3));
    neighbIt = facets.find(vtriple);
    if (neighbIt != facets.end()) {
      Facet f = (*neighbIt).second;
      glue_cells(f.first, f.second, ch, i);
      facets.erase(neighbIt);
      CGAL_assertion(f.first->neighbor(f.second) != NULL);
      CGAL_assertion(ch->neighbor(i) != NULL);
    } else {
      facets[vtriple] = Facet(ch, i);
      CGAL_assertion(ch->neighbor(i) == NULL);
    }
  }

  return ch;
}

template < class TDS_>
typename Triangulation_incremental_builder_3< TDS_ >::Cell_handle
Triangulation_incremental_builder_3< TDS_ >::add_infinite_cell(
  Cell_handle ch0, int i) 
{
  CGAL_assertion(ch0->neighbor(i) == NULL);
  Vertex_handle vh[4];
  vh[i] = t.infinite_vertex();
  vh[(i+1)&3] = ch0->vertex((i+1)&3);
  vh[(i+2)&3] = ch0->vertex((i+3)&3);
  vh[(i+3)&3] = ch0->vertex((i+2)&3);
  Cell_handle ch1 =  t.tds().create_cell(vh[0], vh[1], vh[2], vh[3]);
  // Neighbors are set to NULL
  // Do not set points to the infinite cell. All finite vertices point to
  // finite cells.
  vh[i]->set_cell(ch1);
  glue_cells(ch0, i, ch1, i);
  return ch1;
}

template < class TDS_>
void
Triangulation_incremental_builder_3< TDS_ >::glue_cells(
  Cell_handle ch0, int ind0, Cell_handle ch1, int ind1)
{
  CGAL_assertion(ch0->has_vertex(ch1->vertex((ind1+1)&3)));
  CGAL_assertion(ch0->has_vertex(ch1->vertex((ind1+2)&3)));
  CGAL_assertion(ch0->has_vertex(ch1->vertex((ind1+3)&3)));

  CGAL_assertion(ch1->has_vertex(ch0->vertex((ind0+1)&3)));
  CGAL_assertion(ch1->has_vertex(ch0->vertex((ind0+2)&3)));
  CGAL_assertion(ch1->has_vertex(ch0->vertex((ind0+3)&3)));

  CGAL_assertion(ch0->index(ch1->vertex((ind1+1)&3)) != ind0);
  CGAL_assertion(ch0->index(ch1->vertex((ind1+2)&3)) != ind0);
  CGAL_assertion(ch0->index(ch1->vertex((ind1+3)&3)) != ind0);

  CGAL_assertion(ch1->index(ch0->vertex((ind0+1)&3)) != ind1);
  CGAL_assertion(ch1->index(ch0->vertex((ind0+2)&3)) != ind1);
  CGAL_assertion(ch1->index(ch0->vertex((ind0+3)&3)) != ind1);

  ch0->set_neighbor(ind0, ch1);
  ch1->set_neighbor(ind1, ch0);
}

// Adds infinite cells to the facets on the convex hull
template < class TDS_>
void
Triangulation_incremental_builder_3< TDS_ >::construct_infinite_cells() {
  MapTripleIt ch_facet_it;
  MapPair     ch_edges;
  MapPairIt   ch_edge_it;
  Vertex_handle vh1, vh2;
	
  for (ch_facet_it = facets.begin();
       ch_facet_it != facets.end();
       ch_facet_it ++) {
    Facet ch_facet = (*ch_facet_it).second;
    Cell_handle ch0 = ch_facet.first;
    int ind0 = ch_facet.second;
    Cell_handle ch1 = add_infinite_cell(ch0, ind0);
    // Index of ch1 is also ind0
    CGAL_assertion(ch0->neighbor(ind0) != NULL);
    CGAL_assertion(ch1->neighbor(ind0) != NULL);
		
    for (int i=1; i<4; i++) {
      int i1 = (i==1?2:1);
      int i2 = (i==3?2:3);
      if (ch1->vertex((ind0+i1)&3) < ch1->vertex((ind0+i2)&3)) {
	vh1 = ch1->vertex((ind0+i1)&3);
	vh2 = ch1->vertex((ind0+i2)&3);
      } else {
	vh1 = ch1->vertex((ind0+i2)&3);
	vh2 = ch1->vertex((ind0+i1)&3);
      }
      ch_edge_it = ch_edges.find(Vpair(vh1,vh2));
      if (ch_edge_it != ch_edges.end()) {
	Facet f_opp = (*ch_edge_it).second;
	glue_cells(f_opp.first, f_opp.second, ch1, (ind0+i)&3);
	ch_edges.erase(ch_edge_it);
	CGAL_assertion(f_opp.first->neighbor(f_opp.second) != NULL);
	CGAL_assertion(ch1->neighbor((ind0+i)&3) != NULL);
      } else {
	ch_edges[Vpair(vh1,vh2)] = Facet(ch1, (ind0+i)&3);
	CGAL_assertion(ch1->neighbor((ind0+i)&3) == NULL);
      }
    }	
  }
  CGAL_assertion(ch_edges.empty());
}

} //namespace CGAL

#endif // CGAL_TDS_INCREMENTAL_BUILDER_3_H //
