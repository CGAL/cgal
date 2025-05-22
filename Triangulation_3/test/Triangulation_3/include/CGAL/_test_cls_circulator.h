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
// Author(s)     : Francois Rebufat (Francois.Rebufat@sophia.inria.fr)
//                 Monique Teillaud

#ifndef CGAL_TEST_CLS_CIRCULATOR_C
#define CGAL_TEST_CLS_CIRCULATOR_C

#include <cassert>
#include <CGAL/use.h>

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
  typedef typename Triangulation::Facet            Facet;

  CGAL_USE_TYPE(Vertex_iterator);
  CGAL_USE_TYPE(Cell);
  CGAL_USE_TYPE(Vertex);

  int n = 0;
  Cell_circulator cc, cc0;
  Edge_iterator eit;
  Finite_edges_iterator feit;

  // testing incident_cells(edge *);
//   for (eit=T.edges_begin(); eit!=T.edges_end();  eit++)
  eit=T.edges_begin();
  {
    cc0=cc=T.incident_cells(*eit);

    Cell_circulator i = cc0;
    CGAL_For_all(i, cc0) { *i; }

    Cell_handle ch = cc0; // Test the conversion.
    (void) ch;
    do {
      assert(T.tds().has_vertex(cc, T.tds().vertex(eit->first, eit->second)));
      assert(T.tds().has_vertex(cc, T.tds().vertex(eit->first, eit->third)));
      cc++; n++;
    } while (cc != cc0);
  }
  // test of incident_cells(cellhandle,int,int) and --
  feit=T.finite_edges_begin();
  {
    cc0=cc=T.incident_cells(feit->first, feit->second, feit->third);
      do {
        assert(T.tds().has_vertex(cc, T.tds().vertex(feit->first, feit->second)));
        assert(T.tds().has_vertex(cc, T.tds().vertex(feit->first, feit->third)));
        cc--; n++;
      } while (cc != cc0);
    }
  // testing incident_cells(edge *,cell *); and ++
//   for (eit=T.edges_begin(); eit!=T.edges_end();  eit++)
  eit=T.edges_begin();
  {
    cc0=cc=T.incident_cells(*eit, eit->first);
    do {
      assert(T.tds().has_vertex(cc, T.tds().vertex(eit->first, eit->second)));
      assert(T.tds().has_vertex(cc, T.tds().vertex(eit->first, eit->third)));
      cc++; n++;
    } while (cc != cc0);
  }
  // test of incident_cells(cellhandle,int,int,cellhandle) and --
  feit=T.finite_edges_begin();
  {
    cc0=cc=T.incident_cells(feit->first, feit->second, feit->third,
                            feit->first);
      do {
        assert(T.tds().has_vertex(cc, T.tds().vertex(feit->first, feit->second)));
        assert(T.tds().has_vertex(cc, T.tds().vertex(feit->first, feit->third)));
        cc--; n++;
      } while (cc != cc0);
    }
  // the following is not useful here, it tests iterators more than
  // circulators
//   for (eit=T.finite_edges_begin(); eit!=T.edges_end();  eit++)
//     {
//      cc0=cc=T.incident_cells(*eit);
//       do {
//         cc++; n++;
//       } while (cc != cc0);
//     }
//   for (eit=T.finite_edges_begin(); eit!=T.edges_end(); eit++)
//     {
//      cc0=cc=T.incident_cells(*eit, eit->first);
//       do {
//         cc++; n++;
//       } while (cc != cc0);
//     }

  std::vector<Cell_handle> cells;
  std::vector<Vertex_handle > vertices_old;
  std::vector<Vertex_handle > vertices;
  std::vector<Facet > facets;

  Vertex_handle vh = T.vertices_begin();

  T.incident_cells(vh, std::back_inserter(cells));
  // old name (up to CGAL 3.4)
  // kept for backwards compatibility but not documented
  T.incident_vertices(vh, std::back_inserter(vertices_old));
  // correct name
  T.adjacent_vertices(vh, std::back_inserter(vertices));
  T.incident_facets(vh, std::back_inserter(facets));

  for(typename std::vector<Cell_handle>::const_iterator cit = cells.begin(),
                                                        end = cells.end();
                                                        cit != end; ++cit)
    assert(T.tds().has_vertex(*cit, vh));

  for(typename std::vector<Facet>::const_iterator fit = facets.begin(),
                                                  end = facets.end();
                                                  fit != end; ++fit)
  {
    assert( T.tds().has_vertex(fit->first, vh) );
    assert( fit->second != T.tds().index(fit->first, vh) );
  }

   Facet_circulator fc, fc0, fc1;
   assert(fc1 == nullptr);
   int i,j;
   //   for (eit=T.edges_begin(); eit!=T.edges_end(); eit++)
   eit=T.edges_begin(); // test (edge)
   {
     fc0=fc=T.incident_facets(*eit);
      do {
        assert(T.tds().has_vertex(fc->first, T.tds().vertex(eit->first, eit->second), i));
        assert(T.tds().has_vertex(fc->first, T.tds().vertex(eit->first, eit->third), j));
        assert(fc->second == T.next_around_edge(i, j) );
        fc++; n++;
      } while (fc != fc0);
    }
   //   for (eit=T.edges_begin(); eit!=T.edges_end(); eit++)
   feit=T.finite_edges_begin(); // test (cell*,int,int)
    {
     fc0=fc=T.incident_facets(feit->first, feit->second, feit->third);
      do {
        assert(T.tds().has_vertex(fc->first, T.tds().vertex(feit->first, feit->second), i));
        assert(T.tds().has_vertex(fc->first, T.tds().vertex(feit->first, feit->third), j));
        assert(fc->second == T.next_around_edge(i, j) );
        fc--; n++;
      } while (fc != fc0);
    }
//   int fi;
   //   for (eit=T.edges_begin(); eit!=T.edges_end(); eit++)
    eit=T.edges_begin(); // test (edge, Cell*,int)
    {
//      for (fi=0; fi!=4 ; fi++)
//        {
//         if (t.dimension()==2) {fi=3;}
        fc0=fc=T.incident_facets(*eit, eit->first,
                                 T.next_around_edge(eit->second, eit->third));
        do {
          assert(T.tds().has_vertex(fc->first, T.tds().vertex(eit->first, eit->second), i));
          assert(T.tds().has_vertex(fc->first, T.tds().vertex(eit->first, eit->third), j));
          assert(fc->second == T.next_around_edge(i, j) );
             fc++; n++;
        } while (fc != fc0);
//       }
    }
   //   for (eit=T.edges_begin(); eit!=T.edges_end(); eit++)
   feit=T.finite_edges_begin(); // test (Cell*,int,int,cell*,int)
    {
     fc0=fc=T.incident_facets(feit->first, feit->second, feit->third,
                              feit->first,
                              T.next_around_edge(feit->second, feit->third));
      do {
        assert(T.tds().has_vertex(fc->first, T.tds().vertex(feit->first, feit->second), i));
        assert(T.tds().has_vertex(fc->first, T.tds().vertex(feit->first, feit->third), j));
        assert(fc->second == T.next_around_edge(i, j) );
        fc--; n++;
      } while (fc != fc0);
    }

    eit=T.edges_begin(); // test (edge, Facet)
    {
//      for (fi=0; fi!=4 ; fi++)
//        {
//         if (t.dimension()==2) {fi=3;}
        fc0=fc=T.incident_facets(*eit, std::make_pair( eit->first,
                                        T.next_around_edge(eit->second,
                                                           eit->third)) );
        do {
        assert(T.tds().has_vertex(fc->first, T.tds().vertex(eit->first, eit->second), i));
        assert(T.tds().has_vertex(fc->first, T.tds().vertex(eit->first, eit->third), j));
        assert(fc->second == T.next_around_edge(i, j) );
             fc++; n++;
        } while (fc != fc0);
//       }
    }
   //   for (eit=T.edges_begin(); eit!=T.edges_end(); eit++)
   feit=T.finite_edges_begin(); // test (Cell*,int,int,Facet)
    {
     fc0=fc=T.incident_facets(feit->first, feit->second, feit->third,
                              std::make_pair( feit->first,
                               T.next_around_edge(feit->second,
                                                  feit->third)) );
      do {
        assert(T.tds().has_vertex(fc->first, T.tds().vertex(feit->first, feit->second), i));
        assert(T.tds().has_vertex(fc->first, T.tds().vertex(feit->first, feit->third), j));
        assert(fc->second == T.next_around_edge(i, j) );
        fc--; n++;
      } while (fc != fc0);
    }

  return n;
}

#endif // CGAL_TEST_CLS_CIRCULATOR_C
