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
// file          : include/CGAL/_test_cls_circulator.C
// revision      : 
// revision_date : 
// author(s)     : Francois Rebufat (Francois.Rebufat@sophia.inria.fr)
//                 Monique Teillaud
//
// coordinator   : INRIA Sophia-Antipolis
// ============================================================================

#ifndef CGAL_TEST_CLS_CIRCULATOR_C
#define CGAL_TEST_CLS_CIRCULATOR_C

#include <cassert>

template < class Triangulation >
int
_test_circulator( const Triangulation &T )
{
  typedef typename Triangulation::Finite_edges_iterator   Finite_edges_iterator;
  typedef typename Triangulation::Edge_iterator    Edge_iterator;
  typedef typename Triangulation::Vertex_iterator  Vertex_iterator;
  typedef typename Triangulation::Facet_circulator Facet_circulator;
  typedef typename Triangulation::Cell_circulator  Cell_circulator;
  typedef typename Triangulation::Cell             Cell;
  typedef typename Triangulation::Vertex           Vertex;
  typedef typename Triangulation::Vertex_handle    Vertex_handle;
  typedef typename Triangulation::Cell_handle      Cell_handle;

  int n = 0;
  Cell_circulator cc, cc0;
  Edge_iterator eit;
  Finite_edges_iterator feit;

  // testing incident_cells(edge *); 
//   for (eit=T.edges_begin(); eit!=T.edges_end();  eit++)
  eit=T.edges_begin();
  {
    cc0=cc=T.incident_cells(*eit);
	Cell_handle ch = cc0; // Test the conversion.
	(void) ch;
    do {
      assert(cc->has_vertex((*eit).first->vertex((*eit).second)));
      assert(cc->has_vertex((*eit).first->vertex((*eit).third)));
      cc++; n++;
    } while (cc != cc0);
  }
  // test of incident_cells(cellhandle,int,int) and --
  feit=T.finite_edges_begin();
  {
    cc0=cc=T.incident_cells((*feit).first,(*feit).second,(*feit).third);
      do {
	assert(cc->has_vertex((*feit).first->vertex((*feit).second)));
	assert(cc->has_vertex((*feit).first->vertex((*feit).third)));
	cc--; n++;
      } while (cc != cc0);
    }
  // testing incident_cells(edge *,cell *); and ++
//   for (eit=T.edges_begin(); eit!=T.edges_end();  eit++)
  eit=T.edges_begin();
  {
    cc0=cc=T.incident_cells(*eit,(*eit).first);
    do {
      assert(cc->has_vertex((*eit).first->vertex((*eit).second)));
      assert(cc->has_vertex((*eit).first->vertex((*eit).third)));
      cc++; n++;
    } while (cc != cc0);
  }
  // test of incident_cells(cellhandle,int,int,cellhandle) and --
  feit=T.finite_edges_begin();
  {
    cc0=cc=T.incident_cells((*feit).first,(*feit).second,(*feit).third, 
			    (*feit).first);
      do {
	assert(cc->has_vertex((*feit).first->vertex((*feit).second)));
	assert(cc->has_vertex((*feit).first->vertex((*feit).third)));
	cc--; n++;
      } while (cc != cc0);
    }
  // the following is not useful here, it tests iterators more than
  // circulators 
//   for (eit=T.finite_edges_begin(); eit!=T.edges_end();  eit++)
//     {
//      cc0=cc=T.incident_cells(*eit);
//       do {
// 	cc++; n++;
//       } while (cc != cc0);
//     }
//   for (eit=T.finite_edges_begin(); eit!=T.edges_end(); eit++)
//     {
//      cc0=cc=T.incident_cells(*eit,(*eit).first);
//       do {
// 	cc++; n++;
//       } while (cc != cc0);
//     }

  std::vector<Cell_handle> cells;
  std::vector<Vertex_handle > vertices;

  Vertex_handle vh = T.vertices_begin();

  T.incident_cells(vh,std::back_inserter(cells));
  T.incident_vertices(vh, std::back_inserter(vertices));

   Facet_circulator fc, fc0;
   int i,j;
   //   for (eit=T.edges_begin(); eit!=T.edges_end(); eit++)
   eit=T.edges_begin(); // test (edge)
   {
     fc0=fc=T.incident_facets(*eit);
      do {
	assert((*fc).first->has_vertex((*eit).first->vertex((*eit).second), i));
	assert((*fc).first->has_vertex((*eit).first->vertex((*eit).third), j));
	assert( (*fc).second == T.next_around_edge(i,j) );
	fc++; n++;
      } while (fc != fc0);
    }
   //   for (eit=T.edges_begin(); eit!=T.edges_end(); eit++)
   feit=T.finite_edges_begin(); // test (cell*,int,int)
    {
     fc0=fc=T.incident_facets((*feit).first,(*feit).second,(*feit).third);
      do {
	assert((*fc).first->has_vertex((*feit).first->vertex((*feit).second), i));
	assert((*fc).first->has_vertex((*feit).first->vertex((*feit).third), j));
	assert( (*fc).second == T.next_around_edge(i,j) );
	fc--; n++;
      } while (fc != fc0);
    }
//   int fi;
   //   for (eit=T.edges_begin(); eit!=T.edges_end(); eit++)
    eit=T.edges_begin(); // test (edge, Cell*,int)
    {
//      for (fi=0; fi!=4 ; fi++) 
//        {
// 	if (t.dimension()==2) {fi=3;}
        fc0=fc=T.incident_facets(*eit,(*eit).first,
				 T.next_around_edge((*eit).second,(*eit).third));
        do {
	assert((*fc).first->has_vertex((*eit).first->vertex((*eit).second), i));
	assert((*fc).first->has_vertex((*eit).first->vertex((*eit).third), j));
	assert( (*fc).second == T.next_around_edge(i,j) );
   	  fc++; n++;
        } while (fc != fc0);
//       }
    }
   //   for (eit=T.edges_begin(); eit!=T.edges_end(); eit++)
   feit=T.finite_edges_begin(); // test (Cell*,int,int,cell*,int)
    {
     fc0=fc=T.incident_facets((*feit).first,(*feit).second,(*feit).third,
			      (*feit).first,
			      T.next_around_edge((*feit).second,(*feit).third));
      do {
	assert((*fc).first->has_vertex((*feit).first->vertex((*feit).second), i));
	assert((*fc).first->has_vertex((*feit).first->vertex((*feit).third), j));
	assert( (*fc).second == T.next_around_edge(i,j) );
	fc--; n++;
      } while (fc != fc0);
    }

    eit=T.edges_begin(); // test (edge, Facet)
    {
//      for (fi=0; fi!=4 ; fi++) 
//        {
// 	if (t.dimension()==2) {fi=3;}
        fc0=fc=T.incident_facets(*eit,std::make_pair( (*eit).first,
					T.next_around_edge((*eit).second,
							(*eit).third)) );
        do {
	assert((*fc).first->has_vertex((*eit).first->vertex((*eit).second), i));
	assert((*fc).first->has_vertex((*eit).first->vertex((*eit).third), j));
	assert( (*fc).second == T.next_around_edge(i,j) );
   	  fc++; n++;
        } while (fc != fc0);
//       }
    }
   //   for (eit=T.edges_begin(); eit!=T.edges_end(); eit++)
   feit=T.finite_edges_begin(); // test (Cell*,int,int,Facet)
    {
     fc0=fc=T.incident_facets((*feit).first,(*feit).second,(*feit).third,
			      std::make_pair( (*feit).first,
			       T.next_around_edge((*feit).second,
					       (*feit).third)) );
      do {
	assert((*fc).first->has_vertex((*feit).first->vertex((*feit).second), i));
	assert((*fc).first->has_vertex((*feit).first->vertex((*feit).third), j));
	assert( (*fc).second == T.next_around_edge(i,j) );
	fc--; n++;
      } while (fc != fc0);
    }

  return n;
}

#endif // CGAL_TEST_CLS_CIRCULATOR_C
