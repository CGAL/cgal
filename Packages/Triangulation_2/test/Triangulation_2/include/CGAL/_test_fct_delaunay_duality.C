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
// file          : include/CGAL/_test_fct_delaunay_duality.C
// revision      : 
// revision_date : 
// author(s)     : Herve Bronnimann (Herve.Bronnimann@sophia.inria.fr)
//
// coordinator   : INRIA Sophia-Antipolis
// ============================================================================


template <class Del>
void
_test_delaunay_duality( const Del &T )
{
  //typedef Del                      Cls;

  // We assume the traits class has been tested already
  // actually, any traits is good if it has been tested
  typedef typename Del::Geom_traits          Gt;

  // typedef typename Del::Point                Point;
  //typedef typename Gt::Segment              Segment;
//   typedef typename Del::Triangle             Triangle;

//   typedef typename Del::Distance             Distance;

  //typedef typename Gt::Line                 Line;
//   typedef typename Del::Direction            Direction;
  //typedef typename Gt::Ray                  Ray;

//   typedef typename Del::Vertex               Vertex;
//   typedef typename Del::Face                 Face;

//   typedef typename Del::Vertex_handle        Vertex_handle;
//   typedef typename Del::Face_handle          Face_handle;

//   typedef std::pair<Face_handle,int>         Edge;

  //typedef typename Del::Finite_vertices_iterator      Vertex_iterator;
  typedef typename Del::Finite_faces_iterator        Face_iterator;
  typedef typename Del::Finite_edges_iterator        Edge_iterator;

 //  typedef typename Del::Vertex_circulator    Vertex_circulator;
//   typedef typename Del::Face_circulator      Face_circulator;
//   typedef typename Del::Edge_circulator      Edge_circulator;
//   typedef typename Del::Line_face_circulator Line_face_circulator;

//   typedef typename Del::Locate_type          Locate_type;
  
  // Test dual(face iterator)
  Face_iterator fit;
  for (fit = T.finite_faces_begin(); fit !=  T.finite_faces_end(); ++fit)
    {
      assert( T.side_of_oriented_circle(fit, T.dual(fit)) == 
	      CGAL::ON_POSITIVE_SIDE );
    }
  
  // Test dual(edge iterator)
  Edge_iterator eit;
  for (eit =  T.finite_edges_begin(); eit !=  T.finite_edges_end(); ++eit)
    {
         CGAL::Object o = T.dual(eit);
	 typename Gt::Ray r;
        typename Gt::Segment s;
	typename Gt::Line l;
      if ( CGAL::assign(s,o) ) {
        assert(  ! T.is_infinite((*eit).first) );
	assert( ! T.is_infinite(((*eit).first)->neighbor((*eit).second )) );
      } 
      else if ( CGAL::assign(l,o) ) {
        assert( T.dimension() == 1 );
      } 
      else {
        assert( CGAL::assign(r,o) );
      }
    }

  // Test dual(edge circulator)
 //  Edge_circulator ec=T.finite_vertex()->incident_edges(), done(ec);
//   if ( !ec.is_empty() ) 
//   do  
//     {
//       if (! T.is_infinite(ec)){
// 	// CGAL::Object o = T.dual(ec);
// // 	Segment s; Ray r; Line l;
// // 	assert( CGAL::assign(s,o) || CGAL::assign(r,o) || CGAL::assign(l,o) );
//       }
//       ++ec;
//     } while ( ec == done);
}
