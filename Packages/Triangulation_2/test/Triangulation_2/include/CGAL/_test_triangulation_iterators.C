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
// file          : include/CGAL/_test_triangulation_iterators.C
// revision      : 
// revision_date : 
// author(s)     : Herve Bronnimann (Herve.Bronnimann@sophia.inria.fr)
//
// coordinator   : INRIA Sophia-Antipolis
// ============================================================================


template < class Triangulation >
void
_test_iterators( const Triangulation &T )
{
  int nv = _test_cls_vertex_iterator(T);
  int ne = _test_cls_edge_iterator(T);
  int nf = _test_cls_face_iterator(T);
  int np = _test_cls_point_iterator(T);
  assert( np == nv);
  // std::cout << "Euler's relation: " << nv -ne + nf << std::endl;
  switch (T.dimension()) {
  case 0 : 
  case -1 : assert( nv == T.number_of_vertices() && ne == 0 && nf == 0);
    break;
  case 1 : assert ( nv == T.number_of_vertices() && ne == nv-1 && nf == 0);
    break;
  case 2: assert ( nv == T.number_of_vertices() && nv - ne + nf == 1);
  }
}

template < class Triangulation >
int
_test_cls_face_iterator( const Triangulation &T )
{
  typedef typename Triangulation::Finite_faces_iterator   
                                  Finite_faces_iterator;
  
  int n_finite = 0;
  Finite_faces_iterator fit;

  for (fit = T.finite_faces_begin(); 
       fit != T.finite_faces_end(); 
       ++fit)
    n_finite++;
  assert(n_finite == T.number_of_faces());

  int n=n_finite;
  for (fit = T.finite_faces_end(); 
       fit != T.finite_faces_begin();
       --fit)
    n--;
  assert(n==0);

  return n_finite;
 }


template < class Triangulation >
int
_test_cls_vertex_iterator( const Triangulation &T )
{
  typedef typename Triangulation::Finite_vertices_iterator 
                                  Finite_vertices_iterator;

  int nv = 0;
  Finite_vertices_iterator vit;

  for (vit = T.finite_vertices_begin(); 
       vit != T.finite_vertices_end(); 
       ++vit)
    nv++;
  assert( nv == T.number_of_vertices() );

  int n=nv;
  for (vit = T.finite_vertices_end(); 
       vit != T.finite_vertices_begin(); 
       --vit)
    n--;
  assert( n == 0 );

  return nv;
}

template < class Triangulation >
int
_test_cls_point_iterator( Triangulation &T )
{
  typedef typename Triangulation::Point_iterator Point_iterator;
  typedef typename Triangulation::Point          Point;

  int np = 0;
  Point_iterator pit;
  Point p;
  for (pit = T.points_begin(); 
       pit != T.points_end(); 
       ++pit) {
    np ++;
    p = *pit;
  }
  assert( np == T.number_of_vertices() );

  int n=np;
  for (pit = T.points_end(); 
       pit != T.points_begin(); 
       --pit)
    n--;
  assert( n == 0 );

  return np;  
}

template < class Triangulation >
int
_test_cls_edge_iterator( const Triangulation &T )
{
  typedef typename Triangulation::Finite_edges_iterator   
                                  Finite_edges_iterator;

  int ne = 0;
  Finite_edges_iterator eit;
  for (eit = T.finite_edges_begin(); 
       eit != T.finite_edges_end(); 
       ++eit)
    ne++;
 
  int n = ne;
  for (eit = T.finite_edges_end(); 
       eit != T.finite_edges_begin(); 
       --eit)
    n--;
  assert( n == 0 );

  return ne;
}

