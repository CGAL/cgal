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
// file          : include/CGAL/Delaunay_triangulation_3.h
// revision      : $Revision$
//
// author(s)     : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//
// coordinator   : INRIA Sophia Antipolis (Mariette Yvinec)
//
// ============================================================================


#ifndef CGAL_DELAUNAY_TRIANGULATION_3_H
#define CGAL_DELAUNAY_TRIANGULATION_3_H

#include <set.h>
#include <CGAL/Triangulation_utils_3.h>
#include <CGAL/Triangulation_3.h>

template < class Gt, class Tds>
class CGAL_Delaunay_triangulation_3 : public CGAL_Triangulation_3<Gt,Tds>
{
public:
  
  CGAL_Delaunay_triangulation_3()
    : CGAL_Triangulation_3<Gt,Tds>() {}
  
  CGAL_Delaunay_triangulation_3(const Gt & gt)
  : CGAL_Triangulation_3<Gt,Tds>(gt) {}
  
  CGAL_Delaunay_triangulation_3(const Point & p0,
				const Point & p1,
				const Point & p2,
				const Point & p3)
    : CGAL_Triangulation_3<Gt,Tds>(p0,p1,p2,p3){} // debug

  // copy constructor duplicates vertices and cells
  CGAL_Delaunay_triangulation_3(const CGAL_Delaunay_triangulation_3<Gt,Tds> 
				& tr)
    : CGAL_Triangulation_3<Gt,Tds>(tr)
    { CGAL_triangulation_postcondition( is_valid(true) );  }
  
#ifndef CGAL_CFG_NO_MEMBER_TEMPLATES
  template < class InputIterator >
  int
  insert(InputIterator first, InputIterator last)
  {
    int n = number_of_vertices();
    while(first != last){
      insert(*first);
      ++first;
    }
    return number_of_vertices() - n;
  }
#else
#if defined(LIST_H) || defined(__SGI_STL_LIST_H)
  int
  insert(list<Point>::const_iterator first,
         list<Point>::const_iterator last)
  {
    int n = number_of_vertices();
    while(first != last){
      insert(*first);
      ++first;
    }
    return number_of_vertices() - n;
  }
#endif // LIST_H
#if defined(VECTOR_H) || defined(__SGI_STL_VECTOR_H)
  int
  insert(vector<Point>::const_iterator first,
         vector<Point>::const_iterator last)
  {
    int n = number_of_vertices();
    while(first != last){
      insert(*first);
      ++first;
    }
    return number_of_vertices() - n;
  }
#endif // VECTOR_H
#ifdef ITERATOR_H
  int
  insert(istream_iterator<Point, ptrdiff_t> first,
         istream_iterator<Point, ptrdiff_t> last)
  {
    int n = number_of_vertices();
    while(first != last){
      insert(*first);
      ++first;
    }
    return number_of_vertices() - n;
  }
#endif // ITERATOR_H
  
  int insert(Point* first,
	     Point* last)
  {
    int n = number_of_vertices();
    while(first != last){
      insert(*first);
      ++first;
    }
    return number_of_vertices() - n;
  }
#endif // CGAL_TEMPLATE_MEMBER_FUNCTIONS

  Vertex_handle insert(const Point & p )
  {
    Cell_handle start;
    if ( dimension() >= 1 ) {
      // there is at least one finite "cell" (or facet or edge)
      start = infinite_vertex()->cell()
	->neighbor( infinite_vertex()->cell()->index( infinite_vertex()) );
    }
    else {
      start = NULL;
    }
    return insert( p, start );
  }

  Vertex_handle insert(const Point & p, Cell_handle start)
  {
    switch (dimension()) {
    case 3:
      {
	Locate_type lt;
	int li, lj;
	Cell_handle c = locate( p, start, lt, li, lj);
	switch (lt) {
	case OUTSIDE_CONVEX_HULL:
	case CELL:
	case FACET:
	case EDGE:
	  {
	    Vertex_handle v = new Vertex(p);
	    set_number_of_vertices(number_of_vertices()+1);
	    set<void*, less<void*> > conflicts;
	    Cell_handle aconflict;
	    int ineighbor;
	    find_conflicts_3(conflicts,p,c,aconflict,ineighbor);
	    _tds.star_region(conflicts,&(*v),&(*aconflict),ineighbor);
	    return v;
	  }
	case VERTEX:
          return c->vertex(li);
	default :
          CGAL_triangulation_assertion(false);  // impossible
	}
	break;
      }// dim 3
    case 2:
      {
	Locate_type lt;
	int li, lj;
	Cell_handle c = locate( p, start, lt, li, lj);
	switch (lt) {
	case OUTSIDE_CONVEX_HULL:
	case CELL:
	case FACET:
	case EDGE:
	  {
	    Vertex_handle v = new Vertex(p);
	    set_number_of_vertices(number_of_vertices()+1);
	    set<void*, less<void*> > conflicts;
	    Cell_handle aconflict;
	    int ineighbor;
	    find_conflicts_2(conflicts,p,c,aconflict,ineighbor);
	    _tds.star_region(conflicts,&(*v),&(*aconflict),ineighbor);
	    return v;
	  }
	case VERTEX:
          return c->vertex(li);
	case OUTSIDE_AFFINE_HULL:
	  {
	    // if the 2d triangulation is Delaunay, the 3d
	    // triangulation will be Delaunay
	    return
	      CGAL_Triangulation_3<Gt,Tds>::insert_outside_affine_hull(p); 
	  }
	}
      }//dim 2
    default :
      // dimension <= 1
      return CGAL_Triangulation_3<Gt,Tds>::insert(p);
    }
    return CGAL_Triangulation_3<Gt,Tds>::insert(p);// to avoid warning with egcs
  }// insert(p)

private:
  void
  find_conflicts_3(set<void*, less<void*> > & conflicts, const Point & p,
		     Cell_handle c, Cell_handle & ac, int & i)
    // 3d case
    // p is in conflict with c
    // finds the set conflicts of cells in conflict with p
    // gives a cell ac having a facet on the boundary of conflicts
    // and the index i of its facet on the boundary
  {
    if ( ( conflicts.find( (void *) &(*c) ) ) != conflicts.end() )
      {
	return;   // c was already found
      }
    (void) conflicts.insert( (void *) &(*c) );

    for ( int j=0; j<4; j++ ) {
      if ( side_of_sphere( c->neighbor(j),p ) 
	   ==  CGAL_ON_BOUNDED_SIDE ) {
	find_conflicts_3(conflicts, p, c->neighbor(j), ac, i);
      }
      else {
	ac = c;
	i = j;
      }
    }
  }// find_conflicts_3
  void
  find_conflicts_2(set<void*, less<void*> > & conflicts, const Point & p,
		   Cell_handle c, Cell_handle & ac, int & i)
    // 2d case
    // p is in conflict with c
    // finds the set conflicts of cells in conflict with p
    // gives a cell ac having a facet on the boundary of conflicts
    // and the index i of its facet on the boundary
  {
    if ( ( conflicts.find( (void *) &(*c) ) ) != conflicts.end() )
      {
	return;   // c was already found
      }
    (void) conflicts.insert( (void *) &(*c) );

    for ( int j=0; j<3; j++ ) {
      if ( side_of_circle( c->neighbor(j), 3, p ) 
	   ==  CGAL_ON_BOUNDED_SIDE ) {
	find_conflicts_2(conflicts, p, c->neighbor(j), ac, i);
      }
      else {
	ac = c;
	i = j;
      }
    }
  }// find_conflicts_2
public:

  CGAL_Bounded_side
  side_of_sphere( Cell_handle c, const Point & p) const
  {
    CGAL_triangulation_precondition( dimension() == 3 );
    int i3;
    if ( ! c->has_vertex( infinite_vertex(), i3 ) ) {
      CGAL_Oriented_side
	o = geom_traits().side_of_oriented_sphere(c->vertex(0)->point(),
						  c->vertex(1)->point(),
						  c->vertex(2)->point(),
						  c->vertex(3)->point(),p);
      return ( (o == CGAL_ON_NEGATIVE_SIDE) ? CGAL_ON_UNBOUNDED_SIDE :
	       (o == CGAL_ON_POSITIVE_SIDE) ? CGAL_ON_BOUNDED_SIDE :
	       CGAL_ON_BOUNDARY );
    }
    // infinite cell :
    int i0,i1,i2;
    if ( (i3%2) == 1 ) {
      i0 = (i3+1)&3;
      i1 = (i3+2)&3;
      i2 = (i3+3)&3;
    }
    else {
      i0 = (i3+2)&3;
      i1 = (i3+1)&3;
      i2 = (i3+3)&3;
    }
    CGAL_Orientation
      o = geom_traits().orientation(c->vertex(i0)->point(),
				    c->vertex(i1)->point(),
				    c->vertex(i2)->point(),
				    p);
    switch (o) {
    case CGAL_POSITIVE:
      return CGAL_ON_BOUNDED_SIDE;
    case CGAL_NEGATIVE:
      return CGAL_ON_UNBOUNDED_SIDE;
    case CGAL_ZERO:
      {
	CGAL_Oriented_side s = 
	  geom_traits().side_of_oriented_circle
	  ( c->vertex(i0)->point(),
	    c->vertex(i1)->point(),
	    c->vertex(i2)->point(),
	    p );
	return ( (s == CGAL_ON_NEGATIVE_SIDE) ? CGAL_ON_UNBOUNDED_SIDE :
		 (s == CGAL_ON_POSITIVE_SIDE) ? CGAL_ON_BOUNDED_SIDE :
		 CGAL_ON_BOUNDARY );
// return ( (CGAL_Bounded_side) 
// 		 geom_traits().side_of_oriented_circle
// 		 ( c->vertex(i0)->point(),
// 		   c->vertex(i1)->point(),
// 		   c->vertex(i2)->point(),
// 		   p ) ); 
      }
    }
    return CGAL_ON_UNBOUNDED_SIDE;// to avoid warning with egcs
  }// end side of sphere

  CGAL_Bounded_side
  side_of_circle( const Facet & f, const Point & p) const
    {
      return side_of_circle(f.first, f.second, p);
    }

  CGAL_Bounded_side
  side_of_circle( Cell_handle c, int i, const Point & p) const
    // precondition : dimension >=2
    // in dimension 3, - for a finite facet
    // returns CGAL_ON_BOUNDARY if the point lies on the circle,
    // CGAL_ON_UNBOUNDED_SIDE when exterior, CGAL_ON_BOUNDED_SIDE
    // interior
    // for an infinite facet, considers the plane defined by the
    // adjacent finite facet of the same cell, and does the same as in 
    // dimension 2 in this plane
    // in dimension 2, for an infinite facet
    // in this case, returns CGAL_ON_BOUNDARY if the point lies on the 
    // finite edge (endpoints included) 
    // CGAL_ON_BOUNDED_SIDE for a point in the open half-plane
    // CGAL_ON_UNBOUNDED_SIDE elsewhere
  {
    CGAL_triangulation_precondition( dimension() >= 2 );
    int i3 = 5;

    if ( dimension() == 2 ) {
      CGAL_triangulation_precondition( i == 3 );
      if ( ! c->has_vertex( infinite_vertex(), i3 ) ) {
	CGAL_Oriented_side
	  o = geom_traits().side_of_oriented_circle(c->vertex(0)->point(),
						    c->vertex(1)->point(),
						    c->vertex(2)->point(),
						    p);
	// the triangulation is supposed to be valid, ie the facet
	// with vertices 0 1 2 in this order is positively oriented
	return ( (o == CGAL_ON_NEGATIVE_SIDE) ? CGAL_ON_UNBOUNDED_SIDE :
		 (o == CGAL_ON_POSITIVE_SIDE) ? CGAL_ON_BOUNDED_SIDE :
		 CGAL_ON_BOUNDARY );
      }
      // else infinite facet
      // v1, v2 finite vertices of the facet such that v1,v2,infinite
      // is positively oriented
      Vertex_handle 
	v1 = c->vertex( ccw(i3) ),
	v2 = c->vertex( cw(i3) );
      // does not work in dimension 3 :
      // we need a fourth point to look at orientation with respect
      // to v1v2
      Cell_handle n = c->neighbor(i3);
      CGAL_Orientation o =
	geom_traits().orientation_in_plane( v1->point(), 
					    v2->point(), 
					    n->vertex(n->index(c))->point(),
					    p );
      switch (o) {
      case CGAL_POSITIVE:
	// p lies on the same side of v1v2 as vn, so not in f
	{
	  return CGAL_ON_UNBOUNDED_SIDE;
	}
      case CGAL_NEGATIVE:
	// p lies in f
	{ 
	  return CGAL_ON_BOUNDED_SIDE;
	}
      case CGAL_ZERO:
	// p collinear with v1v2
	{
	  int i_e;
	  Locate_type lt;
	  CGAL_Bounded_side side = 
	    side_of_segment( p,
			     v1->point(), v2->point(),
			     lt, i_e );
	  switch (side) {
	  case CGAL_ON_UNBOUNDED_SIDE:
	    {
	      // p lies on the line defined by the finite edge, but
	      // not in edge v1v2
	      return CGAL_ON_UNBOUNDED_SIDE;
	    }
	  default :
	    {
	      // p lies in edge v1v2 (including v1 or v2)
	      return CGAL_ON_BOUNDARY;
	    }
	  }
	}
      }// switch o
    }// dim 2

    // else dimension == 3
    CGAL_triangulation_precondition( (i >= 0) && (i < 4) );
    if ( ( ! c->has_vertex(infinite_vertex(),i3) ) || ( i3 != i ) ) {
      // finite facet
      // initialization of i0 i1 i2, vertices of the facet positively 
      // oriented (if the triangulation is valid)
      int i0, i1, i2;
      switch (i) {
      case 0:
	{
	  i0 = 1;
	  i1 = 2;
	  i2 = 3;
	  break;
	}
      case 1:
	{
	  i0 = 0;
	  i1 = 2;
	  i2 = 3;
	  break;
	}
      case 2:
	{
	  i0 = 0;
	  i1 = 1;
	  i2 = 3;
	  break;
	}
      case 3:
	{
	  i0 = 0;
	  i1 = 1;
	  i2 = 2;
	  break;
	}
      }
      CGAL_triangulation_precondition( geom_traits().orientation
				       (c->vertex(i0)->point(),
					c->vertex(i1)->point(),
					c->vertex(i2)->point(),
					p) == CGAL_COPLANAR );
      CGAL_Oriented_side
	o = geom_traits().side_of_oriented_circle(c->vertex(i0)->point(),
						  c->vertex(i1)->point(),
						  c->vertex(i2)->point(),
						  p);
      return ( (o == CGAL_ON_NEGATIVE_SIDE) ? CGAL_ON_UNBOUNDED_SIDE :
	       (o == CGAL_ON_POSITIVE_SIDE) ? CGAL_ON_BOUNDED_SIDE :
	       CGAL_ON_BOUNDARY );
    }
    else {//infinite facet
      // v1, v2 finite vertices of the facet such that v1,v2,infinite
      // is positively oriented
      Vertex_handle 
	v1 = c->vertex( nextposaround(i3,i) ),
	v2 = c->vertex( nextposaround(i,i3) );
      CGAL_Orientation o =
	geom_traits().orientation_in_plane( v1->point(),
					    v2->point(),
					    c->vertex(i)->point(),
					    p );
      // then the code is duplicated from 2d case
      switch (o) {
      case CGAL_POSITIVE:
	// p lies on the same side of v1v2 as c->vertex(i), so not in f
	{
	  return CGAL_ON_UNBOUNDED_SIDE;
	}
      case CGAL_NEGATIVE:
	// p lies in f
	{ 
	  return CGAL_ON_BOUNDED_SIDE;
	}
      case CGAL_ZERO:
	// p collinear with v1v2
	{
	  int i_e;
	  Locate_type lt;
	  CGAL_Bounded_side side = 
	    side_of_segment( p,
			     v1->point(), v2->point(),
			     lt, i_e );
	  switch (side) {
	  case CGAL_ON_UNBOUNDED_SIDE:
	    {
	      // p lies on the line defined by the finite edge, but
	      // not in edge v1v2
	      return CGAL_ON_UNBOUNDED_SIDE;
	    }
	  default :
	    {
	      // p lies in edge v1v2 (including v1 or v2)
	      return CGAL_ON_BOUNDARY;
	    }
	  }
	}
      }// switch o
    }// infinite facet
    return CGAL_ON_BOUNDARY; // to avoid warning with egcs
  }// side_of_circle

  
  bool is_valid(bool verbose = false, int level = 0) const
    {
      if ( ! tds().is_valid(verbose,level) ) {
	if (verbose) { cerr << "invalid data structure" << endl; }
	CGAL_triangulation_assertion(false); return false;
      }
    
      if ( &(*infinite_vertex()) == NULL ) {
	if (verbose) { cerr << "no infinite vertex" << endl; }
	CGAL_triangulation_assertion(false); return false;
      }

      int i;

      switch ( dimension() ) {
      case 3:
	{
	  Cell_iterator it;
	  for ( it = finite_cells_begin(); it != cells_end(); ++it ) {
	    is_valid_finite((*it).handle());
	    for ( i=0; i<4; i++ ) {
	      if ( side_of_sphere
		   ( (*it).handle(), 
		     it->vertex( (it->neighbor(i))->index((*it).handle() ) )
		     ->point() )
		   == CGAL_ON_BOUNDED_SIDE ) {
		if (verbose) { 
		  cerr << "non-empty sphere " << endl;
		}
		CGAL_triangulation_assertion(false); return false;
	      }
	    }
	  }
	  break;
	}
      case 2:
	{
	  Facet_iterator it;
	  for ( it = finite_facets_begin(); it != facets_end(); ++it ) {
	    is_valid_finite((*it).first);
	    for ( i=0; i<2; i++ ) {
	      if ( side_of_circle
		   ( (*it).first, 3,
		     (*it).first
		     ->vertex( (((*it).first)->neighbor(i))->index((*it).first) )->point() )
		   == CGAL_ON_BOUNDED_SIDE ) {
		if (verbose) { 
		  cerr << "non-empty circle " << endl;
		}
		CGAL_triangulation_assertion(false); return false;
	      }
	    }
	  }
	  break;
	}
      case 1:
	{
	  Edge_iterator it;
	  for ( it = finite_edges_begin(); it != edges_end(); ++it ) {
	    is_valid_finite((*it).first);
	  }
	  break;
	}
      }
      if (verbose) { cerr << "Delaunay valid triangulation" << endl;}
      return true;
    }

};
#endif CGAL_DELAUNAY_TRIANGULATION_3_H
