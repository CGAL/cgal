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
// file          : _test_cls_face_iterator.C
// revision      : 
// revision_date : 
// author(s)     : Herve Bronnimann (Herve.Bronnimann@sophia.inria.fr)
//
// coordinator   : INRIA Sophia-Antipolis
// ============================================================================


template < class Triangulation >
int
CGAL::_test_cls_face_iterator( const Triangulation &T )
{
  typedef typename Triangulation::Face_iterator   Face_iterator;

  int n_finite = 0;
  Face_iterator fit;
  for (fit = T.faces_begin(); fit != T.faces_end(); ++fit)
    n_finite++;

  typedef typename Triangulation::Face_circulator Face_circulator;
  int n_infinite = 0;
  Face_circulator fcirc = T.infinite_vertex()->incident_faces();
  if ( !fcirc.is_empty() ) 
    do {
      n_infinite++;
      ++fcirc;
    } while (fcirc != T.infinite_vertex()->incident_faces());

  assert( n_infinite + n_finite == T.number_of_faces() );

  return n_finite;
}
