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
// file          : include/CGAL/_test_cls_iterator.C
// revision      : 
// revision_date : 
// author(s)     : Francois Rebufat (Francois.Rebufat@sophia.inria.fr)
//
// coordinator   : INRIA Sophia-Antipolis
// ============================================================================

#ifndef CGAL_TEST_CLS_ITERATOR_C
#define CGAL_TEST_CLS_ITERATOR_C

template < class Triangulation >
int
_test_vertex_iterator( const Triangulation &T )
{
    typedef typename Triangulation::Vertex_iterator Vertex_iterator;
    typedef typename Triangulation::Finite_vertex_iterator Finite_vertex_iterator;
    int n = 0;
    Vertex_iterator vit;
    for (vit = T.vertices_begin(); vit != T.vertices_end(); ++vit)
      n++;
    assert( n-1 == T.number_of_vertices() );
    n=0;
    Finite_vertex_iterator fvit;
   for (fvit = T.finite_vertices_begin(); fvit != T.finite_vertices_end(); ++fvit)
      n++;
    assert( n == T.number_of_vertices() );
    return n;
}

template < class Triangulation >
int
_test_triangulation_iterator( const Triangulation &T )
{
  typedef typename Triangulation::Finite_cell_iterator   Finite_cell_iterator;
  typedef typename Triangulation::Finite_facet_iterator  Finite_facet_iterator;
  typedef typename Triangulation::Finite_edge_iterator   Finite_edge_iterator;
  typedef typename Triangulation::Finite_vertex_iterator Finite_vertex_iterator;

  typedef typename Triangulation::Cell_iterator   Cell_iterator;
  typedef typename Triangulation::Facet_iterator  Facet_iterator;
  typedef typename Triangulation::Edge_iterator   Edge_iterator;
  typedef typename Triangulation::Vertex_iterator Vertex_iterator;

  typedef typename Triangulation::Facet           Facet;
  typedef typename Triangulation::Cell_handle            Cell_handle;
  
  int n=0 , m=0 , f=0 , t=0;
  Cell_iterator Cit;
  Facet_iterator Fit;
  Edge_iterator Eit;
  Vertex_iterator Vit;
  Finite_cell_iterator FCit;
  Finite_facet_iterator FFit;
  Finite_edge_iterator FEit;
  Finite_vertex_iterator FVit;
  if (T.dimension()==3) {
  for (FCit = T.finite_cells_begin(); FCit != T.finite_cells_end(); ++FCit)
     t++;
  for (FFit = T.finite_facets_begin(); FFit != T.finite_facets_end(); ++FFit)
     f++;
  for (FEit = T.finite_edges_begin(); FEit != T.finite_edges_end(); ++FEit)
     m++;
  for (FVit = T.finite_vertices_begin(); FVit != T.finite_vertices_end(); ++FVit)
     n++;
  assert((n-m+f-t)==1);
  n=0 ; m=0 ; f=0 ; t=0;
  for (Cit = T.cells_begin(); Cit != T.cells_end(); ++Cit)
     t++;
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
  Finite_cell_iterator Cit2;
  FCit = T.finite_cells_begin();
  Cit2=FCit;
  assert(T.tetrahedron(&*FCit)==T.tetrahedron(&*Cit2));
  FCit++ ; FCit-- ; ++FCit ; --FCit ;
  assert(FCit==Cit2);
  assert(T.tetrahedron(&*FCit)==T.tetrahedron(&*Cit2));
    }
  if (T.dimension() >=2)
    {
  Finite_facet_iterator Fit2;
  FFit = T.finite_facets_begin();
  Fit2=FFit;
  assert(*FFit==*Fit2);
  FFit++ ; FFit-- ; ++FFit ; --FFit ;
  assert(FFit==Fit2);
  assert(*FFit==*Fit2);
    }
  if (T.dimension() >=1)
    {
  Finite_edge_iterator Eit2;
  FEit = T.finite_edges_begin(); 
  Eit2=FEit;
  assert(*FEit==*Eit2);
  FEit++ ; FEit-- ; ++FEit ; --FEit ;
  assert(FEit==Eit2);
  assert(*FEit==*Eit2);
    }
  Finite_vertex_iterator Vit2;
  FVit = T.finite_vertices_begin(); 
  Vit2=FVit;
  assert(FVit->point()==Vit2->point());
  FVit++ ; FVit-- ; ++FVit ; --FVit ;
  assert(FVit==Vit2);
  assert(FVit->point()==Vit2->point());
  return(n-m+f-t);
}

#endif // CGAL_TEST_CLS_ITERATOR_C
