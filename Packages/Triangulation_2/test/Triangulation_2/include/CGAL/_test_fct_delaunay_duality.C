// ============================================================================
//
// Copyright (c) 1998 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL
// release
// of the Computational Geometry Algorithms Library (CGAL). It is
// not
// intended for general use.
//
// ----------------------------------------------------------------------------
// 
// release       :
// release_date  :
// 
// source        : test_triangulation_tds_2.C
// file          : test_triangulation_tds_2.C
// revision      : 
// revision_date : 
// author(s)     : Herve Bronnimann
// (Herve.Bronnimann@sophia.inria.fr)
//
// coordinator   : INRIA Sophia-Antipolis
// ============================================================================


template <class Triangulation>
void
CGAL__test_delaunay_duality( const Triangulation &T )
{
  typedef Triangulation                      Cls;

  // We assume the traits class has been tested already
  // actually, any traits is good if it has been tested
  typedef typename Cls::Geom_traits          Gt;

  typedef typename Cls::Point                Point;
  typedef typename Cls::Segment              Segment;
  typedef typename Cls::Triangle             Triangle;

  typedef typename Cls::Distance             Distance;

  typedef typename Cls::Line                 Line;
  typedef typename Cls::Direction            Direction;
  typedef typename Cls::Ray                  Ray;

  typedef typename Cls::Vertex               Vertex;
  typedef typename Cls::Face                 Face;

  typedef typename Cls::Vertex_handle        Vertex_handle;
  typedef typename Cls::Face_handle          Face_handle;

  typedef pair<Face_handle,int>              Edge;

  typedef typename Cls::Vertex_iterator      Vertex_iterator;
  typedef typename Cls::Face_iterator        Face_iterator;
  typedef typename Cls::Edge_iterator        Edge_iterator;

  typedef typename Cls::Vertex_circulator    Vertex_circulator;
  typedef typename Cls::Face_circulator      Face_circulator;
  typedef typename Cls::Edge_circulator      Edge_circulator;
  typedef typename Cls::Line_face_circulator Line_face_circulator;

  typedef typename Cls::Locate_type          Locate_type;
  
  // Test dual(face iterator)
  Face_iterator fit;
  for (fit = T.faces_begin(); fit !=  T.faces_end(); ++fit)
    {
      assert( T.side_of_oriented_circle(fit, T.dual(fit)) == CGAL_ON_POSITIVE_SIDE );
    }
  
  // Test dual(edge iterator)
  Edge_iterator eit;
  for (eit =  T.edges_begin(); eit !=  T.edges_end(); ++eit)
    {
      CGAL_Object o = T.dual(eit);
      Segment s; Ray r; Line l;
      if ( CGAL_assign(s,o) ) {
        assert(  ! T.is_infinite((*eit).first) );
	assert( ! T.is_infinite(((*eit).first)->neighbor((*eit).second )) );
      } else if ( CGAL_assign(l,o) ) {
        assert( T.is_infinite((*eit).first) );
	assert( T.is_infinite(((*eit).first)->neighbor((*eit).second )) );
      } else
        assert( CGAL_assign(r,o) );
    }

  // Test dual(edge circulator)
  Edge_circulator ec=T.finite_vertex()->incident_edges(), done(ec);
  if ( !ec.is_empty() ) 
  do  
    {
      if (! T.is_infinite(ec)){
	CGAL_Object o = T.dual(ec);
	Segment s; Ray r; Line l;
	assert( CGAL_assign(s,o) || CGAL_assign(r,o) || CGAL_assign(l,o) );
      }
      ++ec;
    } while ( ec == done);
}
