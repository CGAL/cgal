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


//#include <CGAL/_test_cls_vertex_iterator.C>
//#include <CGAL/_test_cls_edge_iterator.C>
//#include <CGAL/_test_cls_face_iterator.C>

template < class Triangulation >
void
_test_iterators( const Triangulation &T )
{
  int nv = _test_cls_vertex_iterator(T);
  int ne = _test_cls_edge_iterator(T);
  int nf = _test_cls_face_iterator(T);
  // cout << "Euler's relation: " << nv -ne + nf << endl;
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
  for (fit = T.finite_faces_begin(); fit != T.finite_faces_end(); ++fit)
    n_finite++;

  assert(n_finite == T.number_of_faces());
  return n_finite;
}


template < class Triangulation >
int
_test_cls_vertex_iterator( const Triangulation &T )
{
  typedef typename Triangulation::Finite_vertices_iterator 
                                  Finite_vertices_iterator;

  int n = 0;
  Finite_vertices_iterator vit;
  for (vit = T.finite_vertices_begin(); vit != T.finite_vertices_end(); ++vit)
    n++;
  assert( n == T.number_of_vertices() );

  return n;
}

template < class Triangulation >
int
_test_cls_edge_iterator( const Triangulation &T )
{
  typedef typename Triangulation::Finite_edges_iterator   
                                  Finite_edges_iterator;

  int n = 0;
  Finite_edges_iterator eit;
  for (eit = T.finite_edges_begin(); eit != T.finite_edges_end(); ++eit)
    n++;

  return n;
}

