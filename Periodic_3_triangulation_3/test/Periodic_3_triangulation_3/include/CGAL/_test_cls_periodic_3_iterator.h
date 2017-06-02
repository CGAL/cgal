// Copyright (c) 1998  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
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
    assert( n == T.tds().number_of_vertices() );
    n=0;
    for (Vertex_iterator vit = T.vertices_begin();
		    vit != T.vertices_end(); ++vit)
    {
	  Vertex_handle vh = vit; // Test the conversion.
	  const Vertex & v = *vit; // Test operator*;
	  Cell_handle c = vit->cell(); // Test operator->;
	  n++;
	  (void) vh;
	  (void) v;
	  (void) c;
    }
    assert( n == T.tds().number_of_vertices() );

    // Test Backward-ness of the iterators.
    n=0;
    for (Vertex_iterator vit = T.vertices_end(); vit != T.vertices_begin(); --vit)
	{
	  Vertex_handle vh = vit; // Test the conversion.
	  (void) vh;
	  n++;
	}
    assert( n == T.tds().number_of_vertices() );
    n=0;
    for (Vertex_iterator vit = T.vertices_end();
		    vit != T.vertices_begin(); --vit)
    {
	  Vertex_handle vh = vit; // Test the conversion.
	  (void) vh;
	  n++;
    }
    assert( n == T.tds().number_of_vertices() );

    return n;
}

template < class Triangulation >
typename Triangulation::size_type
_test_unique_vertex_iterator( const Triangulation &T )
{
    typedef typename Triangulation::size_type       size_type;
    typedef typename Triangulation::Vertex          Vertex;
    typedef typename Triangulation::Vertex_handle   Vertex_handle;
    typedef typename Triangulation::Cell_handle     Cell_handle;
    typedef typename Triangulation::Unique_vertex_iterator
                                                    Unique_vertex_iterator;

    size_type n = 0;

    for (Unique_vertex_iterator ovit = T.unique_vertices_begin();
	 ovit != T.unique_vertices_end(); ++ovit)
	{
	  Vertex_handle vh = ovit; // Test the conversion.
	  n++;
	  const Vertex & v = *ovit; // Test operator*;
	  Cell_handle c = ovit->cell(); // Test operator->;
	  (void) vh;
	  (void) v;
	  (void) c;
	}
    assert( n == T.number_of_vertices() );

    // Test Backward-ness of the iterators.
    n=0;
    for (Unique_vertex_iterator ovit = T.unique_vertices_end();
	 ovit != T.unique_vertices_begin(); --ovit)
	{
	  Vertex_handle vh = ovit; // Test the conversion.
	  (void) vh;
	  n++;
	}
    assert( n == T.number_of_vertices() );

    return n;
}

template < class Triangulation >
typename Triangulation::size_type
_test_triangulation_iterator( const Triangulation &T )
{
  typedef typename Triangulation::size_type       size_type;
  typedef typename Triangulation::Cell_iterator   Cell_iterator;
  typedef typename Triangulation::Facet_iterator  Facet_iterator;
  typedef typename Triangulation::Edge_iterator   Edge_iterator;
  typedef typename Triangulation::Vertex_iterator Vertex_iterator;

  typedef typename Triangulation::Cell            Cell;
  typedef typename Triangulation::Facet           Facet;
  typedef typename Triangulation::Edge            Edge;
  typedef typename Triangulation::Vertex          Vertex;
  typedef typename Triangulation::Cell_handle     Cell_handle;

  typedef typename Triangulation::Iterator_type Iterator_type;

  typedef typename Triangulation::Periodic_tetrahedron_iterator
    Periodic_tetrahedron_iterator;
  typedef typename Triangulation::Periodic_triangle_iterator
    Periodic_triangle_iterator;
  typedef typename Triangulation::Periodic_segment_iterator
    Periodic_segment_iterator;
  typedef typename Triangulation::Periodic_point_iterator
    Periodic_point_iterator;

  typedef typename Triangulation::Periodic_tetrahedron Periodic_tetrahedron;
  typedef typename Triangulation::Periodic_triangle    Periodic_triangle;
  typedef typename Triangulation::Periodic_segment     Periodic_segment;
  typedef typename Triangulation::Periodic_point       Periodic_point;

  size_type n=0 , m=0 , f=0 , t=0;
  Cell_iterator Cit;
  Facet_iterator Fit;
  Edge_iterator Eit;
  Vertex_iterator Vit;
  Periodic_tetrahedron_iterator T4it;
  Periodic_triangle_iterator T3it;
  Periodic_segment_iterator Sit;
  Periodic_point_iterator Pit;
  for (Cit = T.tds().raw_cells_begin(); Cit != T.tds().raw_cells_end(); ++Cit) {
    Cell_handle ch = Cit;
    (void) ch;
  }
  if (T.number_of_vertices()!=0) {
  for (Cit = T.cells_begin(); Cit != T.cells_end(); ++Cit)
  {
     Cell_handle ch = Cit; // Test the conversion.
     const Cell & c = *Cit; // Test operator*.
     Cell_handle ch2 = Cit->neighbor(0); // Test operator->.
     (void) ch;
     (void) c;
     (void) ch2;
     t++;
  }
  for (Fit = T.facets_begin(); Fit != T.facets_end(); ++Fit) {
     const Facet & f2 = *Fit; // Test operator*.
     Cell_handle ch = Fit->first; // Test operator->.
     (void) f2;
     (void) ch;
     f++;
  }
  for (Eit = T.edges_begin(); Eit != T.edges_end(); ++Eit) {
     const Edge & e = *Eit; // Test operator*.
     Cell_handle ch = Eit->first; // Test operator->.
     (void) e;
     (void) ch;
     m++;
  }
  for (Vit = T.vertices_begin(); Vit != T.vertices_end(); ++Vit) {
     const Vertex & v = *Vit; // Test operator*.
     Cell_handle ch = Vit->cell(); // Test operator->.
     (void) v;
     (void) ch;
     n++;
  }
  for (int i=0 ; i<4 ; i++) {
    Iterator_type it = Iterator_type(i);
    for (T4it = T.periodic_tetrahedra_begin(it);
	 T4it != T.periodic_tetrahedra_end(it); ++T4it) {
      const Periodic_tetrahedron & t = *T4it; // Test operator*.
      const Periodic_point p = T4it->at(0); // Test operator->.
      (void) t;
      (void) p;
    }
    for (T3it = T.periodic_triangles_begin(it);
	 T3it != T.periodic_triangles_end(it); ++T3it) {
      const Periodic_triangle & t = *T3it; // Test operator*.
      const Periodic_point p = T3it->at(0); // Test operator->.
      (void) t;
      (void) p;
    }
    for (Sit = T.periodic_segments_begin(it); 
	 Sit != T.periodic_segments_end(it); ++Sit) {
      const Periodic_segment & s = *Sit; // Test operator*.
      const Periodic_point p = Sit->at(0); // Test operator->.
      (void) s;
      (void) p;
    }
    for (Pit = T.periodic_points_begin(it);
	 Pit != T.periodic_points_end(it); ++Pit) {
      const Periodic_point & pp = *Pit; // Test operator*.
      const typename Triangulation::Point p = Pit->first; // Test operator->.
      (void) pp;
      (void) p;
    }
  }
  // Test Backward-ness of the iterators.
  for (Cit = T.cells_end(); Cit != T.cells_begin(); --Cit) ;
  for (Fit = T.facets_end(); Fit != T.facets_begin(); --Fit) ;
  for (Eit = T.edges_end(); Eit != T.edges_begin(); --Eit) ;
  for (Vit = T.vertices_end(); Vit != T.vertices_begin(); --Vit) ;
  for (int i=0 ; i<4 ; i++) {
    Iterator_type it = Iterator_type(i);
    for (T4it = T.periodic_tetrahedra_end(it);
	 T4it != T.periodic_tetrahedra_begin(it); --T4it) ;
    for (T3it = T.periodic_triangles_end(it);
	 T3it != T.periodic_triangles_begin(it); --T3it) ;
    for (Sit = T.periodic_segments_end(it);
	 Sit != T.periodic_segments_begin(it); --Sit) ;
    for (Pit = T.periodic_points_end(it);
	 Pit != T.periodic_points_begin(it); --Pit) ;
  }

  assert((n-m+f-t)==0);
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
  {
    Cell_iterator Cit2;
    Cit = T.cells_begin();
    Cit2=Cit;
    assert(T.construct_tetrahedron(T.periodic_tetrahedron(Cit))
     == T.construct_tetrahedron(T.periodic_tetrahedron(Cit2)));
    Cit++ ; Cit-- ; ++Cit ; --Cit ;
    assert(Cit==Cit2);
    assert(T.construct_tetrahedron(T.periodic_tetrahedron(Cit))
     == T.construct_tetrahedron(T.periodic_tetrahedron(Cit2)));
  }
  {
    Facet_iterator Fit2;
    Fit = T.facets_begin();
    Fit2=Fit;
    assert(*Fit==*Fit2);
    Fit++ ; Fit-- ; ++Fit ; --Fit ;
    assert(Fit==Fit2);
    assert(*Fit==*Fit2);
  }
  {
    Edge_iterator Eit2;
    Eit = T.edges_begin();
    Eit2=Eit;
    assert(*Eit==*Eit2);
    Eit++ ; Eit-- ; ++Eit ; --Eit ;
    assert(Eit==Eit2);
    assert(*Eit==*Eit2);
  }
  {
    Vertex_iterator Vit2;
    Vit = T.vertices_begin();
    Vit2=Vit;
    assert(Vit->point()==Vit2->point());
    Vit++ ; Vit-- ; ++Vit ; --Vit ;
    assert(Vit==Vit2);
    assert(Vit->point()==Vit2->point());
  }
  {
    Periodic_tetrahedron_iterator T4it2;
    T4it = T.periodic_tetrahedra_begin();
    T4it2=T4it;
    assert(*T4it==*T4it2);
    T4it++ ; T4it-- ; ++T4it ; --T4it ;
    assert(T4it==T4it2);
    assert(*T4it==*T4it2);
  }
  {
    Periodic_triangle_iterator T3it2;
    T3it = T.periodic_triangles_begin();
    T3it2=T3it;
    assert(*T3it==*T3it2);
    T3it++ ; T3it-- ; ++T3it ; --T3it ;
    assert(T3it==T3it2);
    assert(*T3it==*T3it2);
  }
  {
    Periodic_segment_iterator Sit2;
    Sit = T.periodic_segments_begin();
    Sit2=Sit;
    assert(*Sit==*Sit2);
    Sit++ ; Sit-- ; ++Sit ; --Sit ;
    assert(Sit==Sit2);
    assert(*Sit==*Sit2);
  }
  {
    Periodic_point_iterator Pit2;
    Pit = T.periodic_points_begin();
    Pit2=Pit;
    assert(*Pit==*Pit2);
    Pit++ ; Pit-- ; ++Pit ; --Pit ;
    assert(Pit==Pit2);
    assert(*Pit==*Pit2);
  }
  return(n-m+f-t);
}

#endif // CGAL_TEST_CLS_ITERATOR_C
