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
// coordinator   : Mariette Yvinec  <Mariette.Yvinec@sophia.inria.fr>
//
// ============================================================================


#ifndef CGAL_DELAUNAY_TRIANGULATION_3_H
#define CGAL_DELAUNAY_TRIANGULATION_3_H

#include <set.h>

#include <CGAL/Triangulation_3.h>

template < class Gt, class Tds>
class CGAL_Delaunay_triangulation_3 : public CGAL_Triangulation_3<Gt,Tds>
{
public:
  
  CGAL_Delaunay_triangulation_3()
    : CGAL_Triangulation_3<Gt,Tds>() {}
  
  CGAL_Delaunay_triangulation_3(const Gt& gt)
  : CGAL_Triangulation_3<Gt,Tds>(gt) {}
  
  CGAL_Delaunay_triangulation_3(const Point & p0,
				const Point & p1,
				const Point & p2,
				const Point & p3)
    : CGAL_Triangulation_3<Gt,Tds>(p0,p1,p2,p3){} // debug

  // copy constructor duplicates vertices and cells
  CGAL_Delaunay_triangulation_3(const CGAL_Delaunay_triangulation_3<Gt,Tds> &tr)
    : CGAL_Triangulation_3<Gt,Tds>(tr)
    { CGAL_triangulation_postcondition( is_valid() );  }
  
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
  
  int
  insert(Point* first,
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

  Vertex_handle
  insert(const Point & p)
  {
    switch (dimension()) {
    case 3 :
      {
	Locate_type lt;
	int li, lj;
	Cell_handle c = locate( p, lt, li, lj);
	switch (lt) {
	case OUTSIDE_CONVEX_HULL:
	  c = c->neighbor(li);
	  // infinite cell containing p
	case CELL:
	case FACET:
	case EDGE:
	  {
	    Vertex_handle v = new Vertex(p);
	    set_number_of_vertices(number_of_vertices()+1);
	    set<void*, less<void*> > conflicts;
	    Cell_handle aconflict;
	    int ineighbor;
	    find_conflicts(conflicts,p,c,aconflict,ineighbor);
	    tds().star_region(conflicts,&(*v),&(*aconflict),ineighbor);
// 	    v->set_cell( tds().star_tetra(conflicts,v,aconflict,ineighbor) );

// 	    set<void*, less<void*> >::const_iterator it;
// 	    for(it=conflicts.begin(); it!=conflicts.end(); ++it) {
// 	      Delete( (Cell_handle) *it );
// 	    }
	    return v;
	  }
	case VERTEX:
          return c->vertex(li);
	default:
          CGAL_triangulation_assertion(false);  // impossible
	}
	break;
      }
    case 2 :
      {
	// TO BE DONE
	return CGAL_Triangulation_3<Gt,Tds>::insert(p);
	break;
      }
    default :
      return CGAL_Triangulation_3<Gt,Tds>::insert(p);
    }
  }

private:
  void
  find_conflicts(set<void*, less<void*> > & conflicts, const Point & p,
		 Cell_handle c, Cell_handle & ac, int & i)
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
    cerr << "conflit :" << endl;
    pp_cell(c);

    for ( int j=0; j<4; j++ ) {
      if ( side_of_sphere( c->neighbor(j),p ) 
	   ==  CGAL_ON_BOUNDED_SIDE ) {
	find_conflicts(conflicts, p, c->neighbor(j), ac, i);
      }
      else {
	ac = c;
	i = j;
      }
    }
  }// find_conflicts

  CGAL_Bounded_side
  side_of_sphere( Cell_handle c, const Point & p) const
  {
    int i3;
    if ( ! c->has_vertex( infinite_vertex(), i3 ) ) {
      CGAL_Oriented_side
	o = geom_traits().side_of_oriented_sphere(c->vertex(0)->point(),
						  c->vertex(1)->point(),
						  c->vertex(2)->point(),
						  c->vertex(3)->point(),p);
      return ( (o == CGAL_NEGATIVE) ? CGAL_ON_UNBOUNDED_SIDE :
	       (o == CGAL_POSITIVE) ? CGAL_ON_BOUNDED_SIDE :
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
	return ( (CGAL_Bounded_side) 
		 geom_traits().side_of_oriented_circle_in_plane
		 ( c->vertex(i0)->point(),
		   c->vertex(i1)->point(),
		   c->vertex(i2)->point(),
		   p ) ); 
      }
    default:
      // impossible
      CGAL_triangulation_assertion(false);
      return CGAL_ON_UNBOUNDED_SIDE;
    }
  }// end side of sphere
};
#endif CGAL_DELAUNAY_TRIANGULATION_3_H
