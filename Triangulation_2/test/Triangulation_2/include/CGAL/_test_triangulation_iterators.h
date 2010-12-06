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
// file          : include/CGAL/_test_triangulation_iterators.h
// revision      : 
// revision_date : 
// author(s)     : Herve Bronnimann (Herve.Bronnimann@sophia.inria.fr)
//
// coordinator   : INRIA Sophia-Antipolis
// ============================================================================

#ifndef TEST_TRIANGULATION_ITERATORS_C
#define TEST_TRIANGULATION_ITERATORS_C

template < class Triangulation >
typename Triangulation::size_type
_test_cls_face_iterator( const Triangulation &T );


template < class Triangulation >
typename Triangulation::size_type
_test_cls_face_iterator( const Triangulation &T );


template < class Triangulation >
typename Triangulation::size_type
_test_cls_vertex_iterator( const Triangulation &T );


template < class Triangulation >
typename Triangulation::size_type
_test_cls_point_iterator( Triangulation &T );


template < class Triangulation >
typename Triangulation::size_type
_test_cls_edge_iterator( const Triangulation &T );


template < class Triangulation >
void
_test_iterators( const Triangulation &T )
{
  typedef typename Triangulation:: size_type size_type;
  size_type nv = _test_cls_vertex_iterator(T);
  size_type ne = _test_cls_edge_iterator(T);
  size_type nf = _test_cls_face_iterator(T);
  size_type np = _test_cls_point_iterator(T);
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
typename Triangulation::size_type
_test_cls_face_iterator( const Triangulation &T )
{
  typedef typename Triangulation::Finite_faces_iterator   
                                  Finite_faces_iterator;
  typedef typename Triangulation::Face     Face;
  typedef typename Triangulation::Vertex_handle   Vertex_handle;
  typedef typename Triangulation::Face_handle     Face_handle;
  typedef typename Triangulation::size_type       size_type;

  Face f;
  Face_handle fh;
  Vertex_handle vh;

  size_type n_finite = 0;
  Finite_faces_iterator fit;

  for (fit = T.finite_faces_begin(); 
       fit != T.finite_faces_end(); 
       ++fit) {
    f = *fit;
    fh = fit;
    vh = fit->vertex(0);
    n_finite++;
  }
  assert(n_finite == T.number_of_faces());

  size_type n=n_finite;
  for (fit = T.finite_faces_end(); 
       fit != T.finite_faces_begin();
       --fit)
    n--;
  assert(n==0);

  return n_finite;
 }


template < class Triangulation >
typename Triangulation::size_type
_test_cls_vertex_iterator( const Triangulation &T )
{
  typedef typename Triangulation::Finite_vertices_iterator 
                                  Finite_vertices_iterator;
  typedef typename Triangulation::Vertex   Vertex;
  typedef typename Triangulation::Vertex_handle   Vertex_handle;
  typedef typename Triangulation::Face_handle     Face_handle;
  typedef typename Triangulation::size_type       size_type;
  
  Vertex v;
  Face_handle fh;
  Vertex_handle vh;

  size_type nv = 0;
  Finite_vertices_iterator vit;

  for (vit = T.finite_vertices_begin(); 
       vit != T.finite_vertices_end(); 
       ++vit) {
    v = *vit;
    vh = vit;
    fh = vit->face();
    nv++;
  }
  assert( nv == T.number_of_vertices() );

  size_type n = nv;
  for (vit = T.finite_vertices_end(); 
       vit != T.finite_vertices_begin(); 
       --vit)
    n--;
  assert( n == 0 );

  return nv;
}

template < class Triangulation >
typename Triangulation::size_type
_test_cls_point_iterator( Triangulation &T )
{
  typedef typename Triangulation::Point_iterator Point_iterator;
  typedef typename Triangulation::Point          Point;
  typedef typename Triangulation::size_type      size_type;

  size_type np = 0;
  Point_iterator pit;
  Point p;
  for (pit = T.points_begin(); 
       pit != T.points_end(); 
       ++pit) {
    np ++;
    p = *pit;
  }
  assert( np == T.number_of_vertices() );

  size_type n=np;
  for (pit = T.points_end(); 
       pit != T.points_begin(); 
       --pit)
    n--;
  assert( n == 0 );

  return np;  
}

template < class Triangulation >
typename Triangulation::size_type
_test_cls_edge_iterator( const Triangulation &T )
{
  typedef typename Triangulation::Finite_edges_iterator   
                                  Finite_edges_iterator;
  typedef typename Triangulation::Edge     Edge;
  typedef typename Triangulation::Face_handle     Face_handle;
  typedef typename Triangulation::size_type      size_type;

  Edge e;
  Face_handle fh;
 
  size_type ne = 0;
  Finite_edges_iterator eit;
  for (eit = T.finite_edges_begin(); 
       eit != T.finite_edges_end(); 
       ++eit){
    e = *eit;
    fh = eit->first;
    ne++;
  }
 
  size_type n = ne;
  for (eit = T.finite_edges_end(); 
       eit != T.finite_edges_begin(); 
       --eit)
    n--;
  assert( n == 0 );

  return ne;
}

#endif //TEST_TRIANGULATION_ITERATORS_C
