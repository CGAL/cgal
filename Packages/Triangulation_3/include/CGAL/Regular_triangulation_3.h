// ============================================================================
//
// Copyright (c) 1999 The CGAL Consortium
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
// file          : include/CGAL/Regular_triangulation_3.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//                 Sophie Balaven
//                 Sylvain Pion <Sylvain.Pion@sophia.inria.fr>
//                 
//
// coordinator   : Mariette Yvinec <Mariette.Yvinec@sophia.inria.fr>
//
// ============================================================================

#ifndef CGAL_REGULAR_TRIANGULATION_3_H
#define CGAL_REGULAR_TRIANGULATION_3_H

#include <set>
#include <CGAL/Triangulation_utils_3.h>
#include <CGAL/Triangulation_3.h>

CGAL_BEGIN_NAMESPACE

template < class Gt, class Tds>
class Regular_triangulation_3 : public Triangulation_3<Gt,Tds>
{
public:
  typedef typename Gt::Bare_point Bare_point;
  typedef typename Gt::Weighted_point Weighted_point;

  typedef Triangulation_3<Gt,Tds>::Vertex_handle Vertex_handle;
  typedef Triangulation_3<Gt,Tds>::Cell_handle Cell_handle;
  typedef Triangulation_3<Gt,Tds>::Vertex Vertex;
  typedef Triangulation_3<Gt,Tds>::Cell Cell;
  typedef Triangulation_3<Gt,Tds>::Facet Facet;
  typedef Triangulation_3<Gt,Tds>::Edge Edge;

  Regular_triangulation_3()
    : Triangulation_3<Gt,Tds>() {}
  
  Regular_triangulation_3(const Gt & gt)
  : Triangulation_3<Gt,Tds>(gt) {}
  
  // copy constructor duplicates vertices and cells
  Regular_triangulation_3(const Regular_triangulation_3<Gt,Tds> & rt)
      : Triangulation_3<Gt,Tds>(rt)
    { 
      CGAL_triangulation_postcondition( is_valid(true) );  
    }
  
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
  insert(std::list<Point>::const_iterator first,
         std::list<Point>::const_iterator last)
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
  insert(std::vector<Point>::const_iterator first,
         std::vector<Point>::const_iterator last)
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
#endif // CGAL_CFG_NO_MEMBER_TEMPLATES

  Vertex_handle insert( const Weighted_point &p );

  Vertex_handle insert( const Weighted_point & p, Cell_handle start );
    
  Bounded_side
  side_of_power_sphere( Cell_handle c, const Weighted_point &p) const;
  
  Bounded_side
  side_of_power_circle( const Facet & f, const Weighted_point & p) const
    {
      return side_of_power_circle(f.first, f.second, p);
    }

  Bounded_side
  side_of_power_circle( Cell_handle c, int i, const Weighted_point &p) const;

  Bounded_side
  side_of_power_segment( Cell_handle c, const Weighted_point &p) const;

  bool is_valid(bool verbose = false, int level = 0) const;

private:
  bool in_conflict_3(const Weighted_point & p, const Cell_handle & c)
    {
      return side_of_power_sphere(c, p) == ON_BOUNDED_SIDE;
    }

  bool in_conflict_2(const Weighted_point & p, const Cell_handle & c, int i)
    {
      return side_of_power_circle(c, i, p) == ON_BOUNDED_SIDE;
    }

  bool in_conflict_1(const Weighted_point & p, const Cell_handle & c)
    {
      return side_of_power_segment(c, p) == ON_BOUNDED_SIDE;
    }

  void find_conflicts_3(std::set<void*, std::less<void*> > &conflicts, 
			const Weighted_point & p,
			Cell_handle c, Cell_handle & ac, int & i);

  void find_conflicts_2(std::set<void*, std::less<void*> > & conflicts, 
			const Weighted_point & p,
			Cell_handle c, Cell_handle & ac, int & i);

  std::set<void*, std::less<void*> > 
  star_region_delete_points( std::set<void*, std::less<void*> > & region, 
			     Vertex* v, 
			     Cell* c, int li);
    // region is a set of connected cells
    // c belongs to region and has facet i on the boundary of region 
    // replaces the cells in region  
    // by linking v to the boundary of region
    // deleted weighted points that are not in the triangulation
    // anymore, pts will be the list of deleted points
};

template < class Gt, class Tds >
void 
Regular_triangulation_3<Gt,Tds>::
find_conflicts_3(std::set<void*, std::less<void*> > &conflicts, 
		 const Weighted_point & p, 
		 Cell_handle c, Cell_handle & ac, int & i)
  // DUPLICATED CODE !!! see Delaunay
{
  if ( ( conflicts.find( (void *) &(*c) ) ) != conflicts.end() )
    return;   // c was already found
  
  (void) conflicts.insert( (void *) &(*c) );
  
  for ( int j=0; j<4; j++ ) {
    if ( in_conflict_3( p, c->neighbor(j) ) ) {
      find_conflicts_3(conflicts, p, c->neighbor(j), ac, i);
    }
    else {
      ac = c;
      i = j;
    }
  }      
}

template < class Gt, class Tds >
void 
Regular_triangulation_3<Gt,Tds>::
find_conflicts_2(std::set<void*, std::less<void*> > & conflicts, 
		 const Weighted_point & p,
		 Cell_handle c, Cell_handle & ac, int & i)
{
  if ( ( conflicts.find( (void *) &(*c) ) ) != conflicts.end() )
      return;   // c was already found

  (void) conflicts.insert( (void *) &(*c) );

  for ( int j=0; j<3; j++ ) {
    if ( in_conflict_2( p, c->neighbor(j), 3 ) 
	 ==  ON_BOUNDED_SIDE ) {
      find_conflicts_2(conflicts, p, c->neighbor(j), ac, i);
    }
    else {
      ac = c;
      i = j;
    }
  }
}// find_conflicts_2

template < class Gt, class Tds >
Bounded_side
Regular_triangulation_3<Gt,Tds>::
side_of_power_sphere( Cell_handle c, const Weighted_point &p) const
{
  CGAL_triangulation_precondition( dimension() == 3 );
  int i3;
  if ( ! c->has_vertex( infinite_vertex(), i3 ) ) {  
    return Bounded_side( geom_traits().power_test
			 (c->vertex(0)->point(),
			  c->vertex(1)->point(),
			  c->vertex(2)->point(),
			  c->vertex(3)->point(),p) );
  }
  // else infinite cell :
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

  // general case
  Orientation
    o = geom_traits().orientation(c->vertex(i0)->point(),
				  c->vertex(i1)->point(),
				  c->vertex(i2)->point(),
				  p);
  if (o != ZERO)
    return Bounded_side(o);

  // else p coplanar with i0,i1,i2
  return Bounded_side( geom_traits().power_test
		       ( c->vertex(i0)->point(),
			 c->vertex(i1)->point(),
			 c->vertex(i2)->point(), p ) );
}// end side of power sphere

template < class Gt, class Tds >
Bounded_side
Regular_triangulation_3<Gt,Tds>::
side_of_power_circle( Cell_handle c, int i, const Weighted_point &p) const
{
  CGAL_triangulation_precondition( dimension() >= 2 );
  int i3 = 5;
  if ( dimension() == 2 ) {
    CGAL_triangulation_precondition( i == 3 );
    // the triangulation is supposed to be valid, ie the facet
    // with vertices 0 1 2 in this order is positively oriented
    if ( ! c->has_vertex( infinite_vertex(), i3 ) ) 
      return Bounded_side( geom_traits().power_test
			   (c->vertex(0)->point(),
			    c->vertex(1)->point(),
			    c->vertex(2)->point(),
			    p) );
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
    Orientation o =
      geom_traits().orientation_in_plane( v1->point(), 
					  v2->point(), 
					  n->vertex(n->index(c))->point(),
					  p );
    if ( o != ZERO ) return Bounded_side( -o );
    // case when p collinear with v1v2
    return Bounded_side( geom_traits().power_test
			 ( v1->point(), v2->point(), p ) );
  }// dim 2

  // else dimension == 3
  CGAL_triangulation_precondition( (i >= 0) && (i < 4) );
  if ( ( ! c->has_vertex(infinite_vertex(),i3) ) || ( i3 != i ) ) {
    // finite facet
    // initialization of i0 i1 i2, vertices of the facet positively 
    // oriented (if the triangulation is valid)
    int i0 = (i>0) ? 0 : 1;
    int i1 = (i>1) ? 1 : 2;
    int i2 = (i>2) ? 2 : 3;
    CGAL_triangulation_precondition( geom_traits().orientation
				     (c->vertex(i0)->point(),
				      c->vertex(i1)->point(),
				      c->vertex(i2)->point(),
				      p) == COPLANAR );
    return Bounded_side( geom_traits().power_test
			  (c->vertex(i0)->point(),
			   c->vertex(i1)->point(),
			   c->vertex(i2)->point(),
			   p) );
  }
  //else infinite facet
  // v1, v2 finite vertices of the facet such that v1,v2,infinite
  // is positively oriented
  Vertex_handle 
    v1 = c->vertex( nextposaround(i3,i) ),
    v2 = c->vertex( nextposaround(i,i3) );
  Orientation o =
    geom_traits().orientation_in_plane( v1->point(),
					v2->point(),
					c->vertex(i)->point(),
					p );
  // then the code is duplicated from 2d case
  if ( o != ZERO ) return Bounded_side( -o );
  // because p is in f iff 
  // it is not on the same side of v1v2 as c->vertex(i)
  // case when p collinear with v1v2 :
  return Bounded_side( geom_traits().power_test
		       ( v1->point(), v2->point(), p ) );
}

template < class Gt, class Tds >
Bounded_side
Regular_triangulation_3<Gt,Tds>::
side_of_power_segment( Cell_handle c, const Weighted_point &p) const
{
  CGAL_triangulation_precondition( dimension() == 1 );
  if ( ! is_infinite(c,0,1) ) 
    return Bounded_side( geom_traits().power_test
			 ( c->vertex(0)->point(), 
			   c->vertex(1)->point(), 
			   p ) );
  Locate_type lt; int i;
  return side_of_edge( p, c, lt, i );
}

template < class Gt, class Tds >
typename Regular_triangulation_3<Gt,Tds>::Vertex_handle
Regular_triangulation_3<Gt,Tds>::
insert(const Weighted_point & p )
{
  Cell_handle start;
  if ( dimension() >= 1 ) {
    // there is at least one finite "cell" (or facet or edge)
    start = infinite_vertex()->cell()
      ->neighbor(infinite_vertex()->cell()->index(infinite_vertex()));
  }
  else {
    start = NULL;
  }
  return Regular_triangulation_3<Gt,Tds>::insert( p, start );
}  

template < class Gt, class Tds >
typename Regular_triangulation_3<Gt,Tds>::Vertex_handle
Regular_triangulation_3<Gt,Tds>::
insert(const Weighted_point & p, Cell_handle start ) 
{
  switch (dimension()) {
  case 3: 
    {
      Locate_type lt;
      int li, lj;
      Cell_handle c = locate( p, start, lt, li, lj);
      
      if ( lt == VERTEX ) { return c->vertex(li); }
      // TBD : look at the weight...
      
      // else
      Vertex_handle v = new Vertex(p);
      std::set<void*, std::less<void*> > conflicts;
      std::set<void*, std::less<void*> > deleted_points;
      Cell_handle aconflict;
      int ineighbor;
      if (in_conflict_3(p, c)) {
	find_conflicts_3(conflicts, p, c, aconflict, ineighbor);
	deleted_points = 
	  star_region_delete_points(conflicts,&(*v),&(*aconflict),
				    ineighbor);
	// a voir : comment compter les sommets redondants ? (***)
	set_number_of_vertices(number_of_vertices()+1);
      }
      // else : traiter le cas des points redondants a stocker dans
      // la face associee pour le cas d'une future suppression
      return v;
    }
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
	  Vertex_handle v=NULL;
	  std::set<void*, std::less<void*> > deleted_points;
	  if (in_conflict_2(p, c, 3)) {
	    v = new Vertex(p);
	    std::set<void*, std::less<void*> > conflicts;
	    Cell_handle aconflict;
	    int ineighbor;
	    find_conflicts_2(conflicts,p,c,aconflict,ineighbor);
	    deleted_points = 
	      star_region_delete_points(conflicts,&(*v),&(*aconflict),
					ineighbor);
	    set_number_of_vertices(number_of_vertices()+1);
	  }
	  else {
	    Point* P = new Point( p );
	    deleted_points.insert( P );
	  }
	  return v;
	}
	case VERTEX:
	return c->vertex(li);
      case OUTSIDE_AFFINE_HULL:
	{
	  // if the 2d triangulation is Delaunay, the 3d
	  // triangulation will be Delaunay
	  return
	    Triangulation_3<Gt,Tds>::insert_outside_affine_hull(p); 
	}
      }
    }//dim 2
  case 1:
    {
      Locate_type lt;
      int li, lj;
      Cell_handle c = locate( p, start, lt, li, lj);
      switch (lt) {
      case OUTSIDE_CONVEX_HULL:
      case EDGE:
	{
	  Vertex_handle v=NULL;
	  std::set<void*, std::less<void*> > deleted_points;
	  Point * P;
	  if ( in_conflict_1(p, c) ) {
	    v = new Vertex(p);
	    set_number_of_vertices(number_of_vertices()+1);
	    Cell_handle bound[2];
	    Cell_handle n;
	    std::set<void*, std::less<void*> > conflicts;

	    for (int j =0; j<2; j++ ) {
	      n = c;
	      while ( ( ! is_infinite(n->vertex(1-j)) ) && 
		      in_conflict_1( p, n->neighbor(j) ) ) {
		if (n!=c) (void) conflicts.insert( (void *) &(*n) );
		P = new Point( n->vertex(1-j)->point() );
		(void) deleted_points.insert((void*) P);
		set_number_of_vertices(number_of_vertices()-1);
		n = n->neighbor(j);
	      }
	      bound[j] = n;
	    }
	    if ( bound[0] != bound[1] ) {
	      if ( (c != bound[0]) && (c != bound[1]) ) {
		(void) conflicts.insert( (void *) &(*c) );
	      }
	      bound[0]->set_vertex(0,v);
	      v->set_cell(bound[0]);
	      bound[1]->set_vertex(1,v);
	      bound[1]->set_neighbor(0,bound[0]);
	      bound[0]->set_neighbor(1,bound[1]);
	    }
	    else {
	      bound[1] = new Cell(bound[0]->vertex(0), v, NULL, NULL,
				  bound[0], bound[0]->neighbor(1), NULL, NULL);
	      _tds.add_cell(&(*bound[1]));
	      bound[0]->neighbor(1)->set_neighbor(0,bound[1]);
	      bound[0]->vertex(0)->set_cell(bound[1]);

	      bound[0]->set_neighbor(1,bound[1]);
	      bound[0]->set_vertex(0,v);
	      v->set_cell(bound[0]);
	    }

	    std::set<void*, std::less<void*> >::const_iterator it;
	    for ( it = conflicts.begin(); it != conflicts.end(); ++it) {
	      delete((Cell*)*it);
	    }
	  }
	  else {
	    Point* P = new Point( p );
	    deleted_points.insert( P );
	  }
	  return v;
	}
      case VERTEX:
	return c->vertex(li);
      case OUTSIDE_AFFINE_HULL:
	{
	  return
	    Triangulation_3<Gt,Tds>::insert_outside_affine_hull(p); 
	}
      case FACET:
      case CELL:
	// impossible in dimension 1
	return NULL;
      }
    }
  default :
    {
      // temporary : will only work for non degenerated dimensions
      // (only for the first 4 points if they form a true tetrahedron)
      return Triangulation_3<Gt,Tds>::insert(p,start);
    }
  }
}

template < class Gt, class Tds >
std::set<void*, std::less<void*> > 
Regular_triangulation_3<Gt,Tds>::
star_region_delete_points(std::set<void*, std::less<void*> > & region, 
			  Vertex* v, Cell* c, int li) 
  // region is a set of connected cells
  // c belongs to region and has facet i on the boundary of region 
  // replaces the cells in region  
  // by linking v to the boundary of region
  // deleted weighted points that are not in the triangulation
  // anymore, pts will be the list of deleted points
{
  std::set<void*, std::less<void*> > vert;
  Cell *c_tmp;
  Vertex *v_tmp;
  Vertex_handle vh;
  std::set<void*, std::less<void*> > pts;
  Point *p;
    
  // for each cell to be deleted, keep vertices
  std::set<void*, std::less<void*> >::const_iterator it;
  for( it = region.begin(); it != region.end(); ++it) {
    c_tmp = (Cell *) *it;
    // store vertices
    for (int i=0; i<=dimension() ; i++){
      vh = c_tmp->vertex(i);
      if ( (vert.find((void*) &(*(vh))) == vert.end()) &&
	   (! is_infinite(vh)) ) {
	vert.insert( (void*) &(*(vh) ) );
      }
    }
  }
    
  // Create the new faces and delete old ones
  _tds.star_region( region, v, c, li );
    
  // get the vertices incident to v
  std::set<Vertex*, std::less<Vertex*> > inc_vert;
  incident_vertices(v, inc_vert);
    
  // for each vertex, check if it is a vertex incident to v
  // if not, delete it
  for( it = vert.begin(); it != vert.end(); ++it) {
    v_tmp = (Vertex *) *it;
    if ( (inc_vert.find( v_tmp )) == inc_vert.end() ) {
      // vertex has to be deleted and point to be stored
      p = new Point( v_tmp->point() );
      pts.insert( p );
      set_number_of_vertices(number_of_vertices()-1);
      delete((Vertex *)*it);
    }
  }
    
  // returns list of deleted points
  return pts;
}

template < class Gt, class Tds >
bool 
Regular_triangulation_3<Gt,Tds>::
is_valid(bool verbose = false, int level = 0) const 
{
  if ( ! tds().is_valid(verbose,level) ) {
    if (verbose) { std::cerr << "invalid data structure" << std::endl; }
    CGAL_triangulation_assertion(false); return false;
  }

  if ( &(*infinite_vertex()) == NULL ) {
    if (verbose) { std::cerr << "no infinite vertex" << std::endl; }
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
	  if ( side_of_power_sphere
	       ( (*it).handle(), 
		 it->vertex( (it->neighbor(i))->index((*it).handle() ) )
		 ->point() )
	       == ON_BOUNDED_SIDE ) {
	    if (verbose) { 
	      std::cerr << "non-empty sphere " << std::endl;
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
	for ( i=0; i<3; i++ ) {
	  if ( side_of_power_circle
	       ( (*it).first, 3,
		 (*it).first->vertex( (((*it).first)->neighbor(i))
				      ->index((*it).first) )->point() )
	       == ON_BOUNDED_SIDE ) {
	    if (verbose) { 
	      std::cerr << "non-empty circle " << std::endl;
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
	for ( i=0; i<2; i++ ) {
	  if ( side_of_power_segment
	       ( (*it).first,
		 (*it).first->vertex( (((*it).first)->neighbor(i))
				      ->index((*it).first) )->point() )
	       == ON_BOUNDED_SIDE ) {
	    if (verbose) { 
	      std::cerr << "non-empty edge " << std::endl;
	    }
	    CGAL_triangulation_assertion(false); return false;
	  }
	}
      }
      break;
    }
  }
  if (verbose) { std::cerr << "Regular valid triangulation" << std::endl;}
  return true;
}
CGAL_END_NAMESPACE

#endif // CGAL_REGULAR_TRIANGULATION_3_H
