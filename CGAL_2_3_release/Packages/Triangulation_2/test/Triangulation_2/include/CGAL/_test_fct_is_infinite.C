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
// file          : include/CGAL/_test_fct_is_infinite.C
// revision      : 
// revision_date : 
// author(s)     : Herve Bronnimann (Herve.Bronnimann@sophia.inria.fr)
//
// coordinator   : INRIA Sophia-Antipolis
// ============================================================================


template <class Triangulation>
void
_test_fct_is_infinite( const Triangulation &T )
{
  typedef typename Triangulation::Edge    Edge;
  typedef typename Triangulation::Finite_edges_iterator
                                  Finite_edges_iterator;

  // test infinite_face() and is_infinite(Face_handle)
  if ( !T.infinite_face().is_null() )
    assert( T.is_infinite(T.infinite_face()) );

  // test infinite_vertex() and is_infinite(Vertex_handle)
  assert( T.is_infinite(T.infinite_vertex()) );
  
  // test finite_vertex()
  // missing in the documentation is the precondition that there
  // must be at least one vertex
  if (T.number_of_vertices() != 0)
      assert( !T.is_infinite(T.finite_vertex()) );

  // test is_infinite(Edge)
  // an infinite face always has two infinite edges
  if ( !T.infinite_face().is_null() )
    { 
      int index = T.infinite_face()->index(T.infinite_vertex());
      assert( T.is_infinite( Edge(T.infinite_face(),(index+1)%3)) );
      assert( T.is_infinite( Edge(T.infinite_face(),(index+2)%3)) );
    }

  // test is_infinite(Edge_circulator)
  if ( ! T.infinite_vertex()->incident_edges().is_empty() )
    assert( T.is_infinite( T.infinite_vertex()->incident_edges()) );

  // test is_infinite(Edge_iterator)
  for(Finite_edges_iterator fei = T.finite_edges_begin();
                           fei != T.finite_edges_end(); fei++) {
    assert( ! T.is_infinite( fei));
  }
}
