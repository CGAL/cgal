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
    int n = 0;
    Vertex_iterator vit;
    for (vit = T.all_vertices_begin(); vit != T.vertices_end(); ++vit)
      n++;
    assert( n-1 == T.number_of_vertices() );
    n=0;
   for (vit = T.finite_vertices_begin(); vit != T.vertices_end(); ++vit)
      n++;
    assert( n == T.number_of_vertices() );
    return n;

}

template < class Triangulation >
int
_test_triangulation_iterator( const Triangulation &T )
{
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
  if (T.dimension()==3) {
  for (Cit = T.finite_cells_begin(); Cit != T.cells_end(); ++Cit)
     t++;
  for (Fit = T.finite_facets_begin(); Fit != T.facets_end(); ++Fit)
     f++;
  for (Eit = T.finite_edges_begin(); Eit != T.edges_end(); ++Eit)
     m++;
  for (Vit = T.finite_vertices_begin(); Vit != T.vertices_end(); ++Vit)
     n++;
  assert((n-m+f-t)==1);
  n=0 ; m=0 ; f=0 ; t=0;
  for (Cit = T.all_cells_begin(); Cit != T.cells_end(); ++Cit)
     t++;
  for (Fit = T.all_facets_begin(); Fit != T.facets_end(); ++Fit)
     f++;
  for (Eit = T.all_edges_begin(); Eit != T.edges_end(); ++Eit)
     m++;
  for (Vit = T.all_vertices_begin(); Vit != T.vertices_end(); ++Vit)
     n++;
  assert((n-m+f-t)==0);
  }
  if (T.dimension()==3)
    {
  Cell_iterator Cit2;
  Cit = T.finite_cells_begin();
  Cit2=Cit;
  assert(T.tetrahedron((Cell_handle) Cit)==T.tetrahedron((Cell_handle) Cit2));
  Cit++ ; Cit-- ; ++Cit ; --Cit ;
  assert(Cit==Cit2);
  assert(T.tetrahedron((Cell_handle) Cit)==T.tetrahedron((Cell_handle) Cit2));
    }
  if (T.dimension() >=2)
    {
  Facet_iterator Fit2;
  Fit = T.finite_facets_begin();
  Fit2=Fit;
  assert(*Fit==*Fit2);
  Fit++ ; Fit-- ; ++Fit ; --Fit ;
  assert(Fit==Fit2);
  assert(*Fit==*Fit2);
    }
  if (T.dimension() >=1)
    {
  Edge_iterator Eit2;
  Eit = T.finite_edges_begin(); 
  Eit2=Eit;
  assert(*Eit==*Eit2);
  Eit++ ; Eit-- ; ++Eit ; --Eit ;
  assert(Eit==Eit2);
  assert(*Eit==*Eit2);
    }
  Vertex_iterator Vit2;
  Vit = T.finite_vertices_begin(); 
  Vit2=Vit;
  assert(Vit->point()==Vit2->point());
  Vit++ ; Vit-- ; ++Vit ; --Vit ;
  assert(Vit==Vit2);
  assert(Vit->point()==Vit2->point());
  return(n-m+f-t);
}

#endif // CGAL_TEST_CLS_ITERATOR_C
