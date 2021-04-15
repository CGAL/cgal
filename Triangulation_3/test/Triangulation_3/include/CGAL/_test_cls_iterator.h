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

#ifndef CGAL_TEST_CLS_ITERATOR_C
#define CGAL_TEST_CLS_ITERATOR_C

template < class Triangulation >
typename Triangulation::size_type
_test_vertex_iterator( const Triangulation &T )
{
    typedef typename Triangulation::size_type       size_type;
    typedef typename Triangulation::Vertex          Vertex;
    typedef typename Triangulation::Vertex_handle   Vertex_handle;
    typedef typename Triangulation::Cell_handle     Cell_handle;
    typedef typename Triangulation::Vertex_iterator Vertex_iterator;
    typedef typename Triangulation::Finite_vertices_iterator
                                                    Finite_vertices_iterator;
    size_type n = 0;

    for (Vertex_iterator vit = T.vertices_begin(); vit != T.vertices_end(); ++vit)
        {
          Vertex_handle vh = vit; // Test the conversion.
          n++;
          const Vertex & v = *vit; // Test operator*;
          Cell_handle c = vit->cell(); // Test operator->;
          (void) vh;
          (void) v;
          (void) c;
        }
    assert( n-1 == T.number_of_vertices() );
    n=0;
    for (Finite_vertices_iterator fvit = T.finite_vertices_begin();
                    fvit != T.finite_vertices_end(); ++fvit)
    {
          Vertex_handle vh = fvit; // Test the conversion.
          const Vertex & v = *fvit; // Test operator*;
          Cell_handle c = fvit->cell(); // Test operator->;
          n++;
          (void) vh;
          (void) v;
          (void) c;
    }
    assert( n == T.number_of_vertices() );

    // Test Backward-ness of the iterators.
    n=0;
    for (Vertex_iterator vit = T.vertices_end(); vit != T.vertices_begin(); --vit)
        {
          Vertex_handle vh = vit; // Test the conversion.
          (void) vh;
          n++;
        }
    assert( n-1 == T.number_of_vertices() );
    n=0;
    for (Finite_vertices_iterator fvit = T.finite_vertices_end();
                    fvit != T.finite_vertices_begin(); --fvit)
    {
          Vertex_handle vh = fvit; // Test the conversion.
          (void) vh;
          n++;
    }
    assert( n == T.number_of_vertices() );

    return n;
}

template < class Triangulation >
int
_test_triangulation_iterator( const Triangulation &T )
{
  typedef typename Triangulation::Finite_cells_iterator  Finite_cells_iterator;
  typedef typename Triangulation::Finite_facets_iterator Finite_facets_iterator;
  typedef typename Triangulation::Finite_edges_iterator  Finite_edges_iterator;
  typedef typename Triangulation::Finite_vertices_iterator Finite_vertices_iterator;

  typedef typename Triangulation::Cell_iterator   Cell_iterator;
  typedef typename Triangulation::Facet_iterator  Facet_iterator;
  typedef typename Triangulation::Edge_iterator   Edge_iterator;
  typedef typename Triangulation::Vertex_iterator Vertex_iterator;

  typedef typename Triangulation::Vertex_handle Vertex_handle;

  typedef typename Triangulation::All_vertex_handles All_vertex_handles;
  typedef typename Triangulation::All_cell_handles All_cell_handles;
  typedef typename Triangulation::All_edges All_edges;
  typedef typename Triangulation::All_facets All_facets;
  typedef typename Triangulation::Finite_vertex_handles Finite_vertex_handles;
  typedef typename Triangulation::Finite_cell_handles Finite_cell_handles;
  typedef typename Triangulation::Finite_edges Finite_edges;
  typedef typename Triangulation::Finite_facets Finite_facets;
  typedef typename Triangulation::Points Points;

  typedef typename Triangulation::Cell            Cell;
  typedef typename Triangulation::Facet           Facet;
  typedef typename Triangulation::Edge            Edge;
  typedef typename Triangulation::Vertex          Vertex;
  typedef typename Triangulation::Cell_handle     Cell_handle;

  int n=0 , m=0 , f=0 , t=0;
  Cell_iterator Cit;
  Facet_iterator Fit;
  Edge_iterator Eit;
  Vertex_iterator Vit;
  Finite_cells_iterator FCit;
  Finite_facets_iterator FFit;
  Finite_edges_iterator FEit;
  Finite_vertices_iterator FVit;
  for (Cit = T.tds().raw_cells_begin(); Cit != T.tds().raw_cells_end(); ++Cit) {
    Cell_handle ch = Cit;
    (void) ch;
  }
  if (T.dimension()==3) {
    {
      All_vertex_handles range = T.all_vertex_handles();
      Vertex_handle vh = *(range.first);
      assert(vh == T.all_vertices_begin());
      vh = *(range.second);
      assert(vh == T.all_vertices_end());
    }
    {
      All_cell_handles range = T.all_cell_handles();
      Cell_handle vh = *(range.first);
      assert(vh == T.all_cells_begin());
      vh = *(range.second);
      assert(vh == T.all_cells_end());
    }
    {
      All_edges range = T.all_edges();
      assert(range.first == T.all_edges_begin());
      assert(range.second == T.all_edges_end());
    }
    {
      All_facets range = T.all_facets();
      assert(range.first == T.all_facets_begin());
      assert(range.second == T.all_facets_end());
    }
    {
      Finite_vertex_handles range = T.finite_vertex_handles();
      Vertex_handle vh = *(range.first);
      assert(vh == Vertex_handle(T.finite_vertices_begin()));
      vh = *(range.second);

      assert(vh == Vertex_handle(T.finite_vertices_end()));
    }
    {
      Finite_cell_handles range = T.finite_cell_handles();
      Cell_handle ch = *(range.first);
      assert(ch == Cell_handle(T.finite_cells_begin()));
      ch = *(range.second);
      assert(ch == Cell_handle(T.finite_cells_end()));
    }
    {
      Finite_edges range = T.finite_edges();
      assert(range.first == T.finite_edges_begin());
      assert(range.second == T.finite_edges_end());
    }
    {
      Finite_facets range = T.finite_facets();
      assert(range.first == T.finite_facets_begin());
      assert(range.second == T.finite_facets_end());
    }
    {
      Points range = T.points();
      assert(range.first == T.points_begin());
      assert(range.second == T.points_end());
    }


  for (FCit = T.finite_cells_begin(); FCit != T.finite_cells_end(); ++FCit)
  {
     Cell_handle ch = FCit; // Test the conversion.
     const Cell & c = *FCit; // Test operator*.
     Cell_handle ch2 = FCit->neighbor(0); // Test operator->.
     (void) ch;
     (void) c;
     (void) ch2;
     t++;
  }
  for (FFit = T.finite_facets_begin(); FFit != T.finite_facets_end(); ++FFit) {
     const Facet & f2 = *FFit; // Test operator*.
     Cell_handle ch = FFit->first; // Test operator->.
     (void) f2;
     (void) ch;
     f++;
  }
  for (FEit = T.finite_edges_begin(); FEit != T.finite_edges_end(); ++FEit) {
     const Edge & e = *FEit; // Test operator*.
     Cell_handle ch = FEit->first; // Test operator->.
     (void) e;
     (void) ch;
     m++;
  }
  for (FVit = T.finite_vertices_begin(); FVit != T.finite_vertices_end(); ++FVit) {
     const Vertex & v = *FVit; // Test operator*.
     Cell_handle ch = FVit->cell(); // Test operator->.
     (void) v;
     (void) ch;
     n++;
  }
  // Test Backward-ness of the iterators.
  for (FCit = T.finite_cells_end(); FCit != T.finite_cells_begin(); --FCit) ;
  for (FFit = T.finite_facets_end(); FFit != T.finite_facets_begin(); --FFit) ;
  for (FEit = T.finite_edges_end(); FEit != T.finite_edges_begin(); --FEit) ;
  for (FVit = T.finite_vertices_end(); FVit != T.finite_vertices_begin(); --FVit) ;

  assert((n-m+f-t)==1);
  n=0 ; m=0 ; f=0 ; t=0;
  for (Cit = T.cells_begin(); Cit != T.cells_end(); ++Cit)
  {
     Cell_handle ch = Cit; // Test the conversion.
     (void) ch;
     t++;
  }
  for (Fit = T.facets_begin(); Fit != T.facets_end(); ++Fit)
     f++;
  for (Eit = T.edges_begin(); Eit != T.edges_end(); ++Eit)
     m++;
  for (Vit = T.vertices_begin(); Vit != T.vertices_end(); ++Vit)
     n++;
  assert((n-m+f-t)==0);
  }
  if (T.dimension()==3)
    {
  Finite_cells_iterator Cit2;
  FCit = T.finite_cells_begin();
  Cit2=FCit;
  assert(T.tetrahedron(FCit)==T.tetrahedron(Cit2));
  FCit++ ; FCit-- ; ++FCit ; --FCit ;
  assert(FCit==Cit2);
  assert(T.tetrahedron(FCit)==T.tetrahedron(Cit2));
    }
  if (T.dimension() >=2)
    {
  Finite_facets_iterator Fit2;
  FFit = T.finite_facets_begin();
  Fit2=FFit;
  assert(*FFit==*Fit2);
  FFit++ ; FFit-- ; ++FFit ; --FFit ;
  assert(FFit==Fit2);
  assert(*FFit==*Fit2);
    }
  if (T.dimension() >=1)
    {
  Finite_edges_iterator Eit2;
  FEit = T.finite_edges_begin();
  Eit2=FEit;
  assert(*FEit==*Eit2);
  FEit++ ; FEit-- ; ++FEit ; --FEit ;
  assert(FEit==Eit2);
  assert(*FEit==*Eit2);
    }
  Finite_vertices_iterator Vit2;
  FVit = T.finite_vertices_begin();
  Vit2=FVit;
  assert(FVit->point()==Vit2->point());
  FVit++ ; FVit-- ; ++FVit ; --FVit ;
  assert(FVit==Vit2);
  assert(FVit->point()==Vit2->point());
  return(n-m+f-t);
}

#endif // CGAL_TEST_CLS_ITERATOR_C
