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
// file          : include/CGAL/Triangulation_3.h
// revision      : $Revision$
//
// author(s)     : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//
// coordinator   : Mariette Yvinec  <Mariette.Yvinec@sophia.inria.fr>
//
// ============================================================================


#ifndef CGAL_TRIANGULATION_3_H
#define CGAL_TRIANGULATION_3_H

#include <iostream.h>

#include <list.h>
#include <map.h> 
//#include <algo.h>
#include <pair.h>
#include <CGAL/triple.h>

#include <CGAL/Pointer.h>
#include <CGAL/circulator.h>

#include <CGAL/predicates_on_points_3.h>

#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_short_names_3.h>

#include <CGAL/Triangulation_data_structure_3.h>

#include <CGAL/Triangulation_cell_3.h>
#include <CGAL/Triangulation_vertex_3.h>
#include <CGAL/Triangulation_handles_3.h>
#include <CGAL/Triangulation_iterators_3.h>
#include <CGAL/Triangulation_circulators_3.h>

template < class GT, class Tds>
class CGAL_Triangulation_cell_iterator_3;

template < class GT, class Tds>
class CGAL_Triangulation_vertex_iterator_3;

template < class GT, class Tds>
class CGAL_Triangulation_edge_iterator_3;

template < class GT, class Tds>
class CGAL_Triangulation_facet_iterator_3;

template < class GT, class Tds>
class CGAL_Triangulation_cell_circulator_3;

template < class GT, class Tds>
class CGAL_T_cell;

template < class GT, class Tds >
class CGAL_Triangulation_3
{
//   friend istream& operator>> CGAL_NULL_TMPL_ARGS
//   (istream& is, CGAL_Triangulation_3<GT,Tds> &tr);
//   friend ostream& operator<< CGAL_NULL_TMPL_ARGS
//   (ostream& os, const CGAL_Triangulation_3<GT,Tds> &tr);

  friend class CGAL_T_cell<GT,Tds>;
  friend class CGAL_Triangulation_vertex_3<GT,Tds>;

  friend CGAL_Triangulation_cell_iterator_3<GT,Tds>;
  friend CGAL_Triangulation_facet_iterator_3<GT,Tds>;
  friend CGAL_Triangulation_edge_iterator_3<GT,Tds>;
  friend CGAL_Triangulation_vertex_iterator_3<GT,Tds>;
  friend CGAL_Triangulation_cell_circulator_3<GT,Tds>;

public:
  typedef CGAL_Triangulation_3<GT,Tds> Triangulation_3;

  typedef typename GT::Point Point;
  typedef typename GT::Segment Segment;
  typedef typename GT::Triangle Triangle;
  typedef typename GT::Tetrahedron Tetrahedron;

  typedef CGAL_Triangulation_cell_handle_3<GT,Tds> Cell_handle;
  typedef CGAL_Triangulation_vertex_handle_3<GT,Tds> Vertex_handle;

  typedef CGAL_T_cell<GT,Tds> Cell;
  typedef CGAL_Triangulation_vertex_3<GT,Tds> Vertex;
  typedef pair<Cell_handle, int> Facet;
  typedef CGAL_triple<Cell_handle, int, int> Edge;


  typedef CGAL_Triangulation_cell_circulator_3<GT,Tds>      Cell_circulator;
//   typedef CGAL_Triangulation_edge_circulator_3<GT,Tds>      Edge_circulator;
//   typedef CGAL_Triangulation_vertex_circulator_3<GT,Tds>    Vertex_circulator;

  typedef CGAL_Triangulation_cell_iterator_3<GT,Tds>   Cell_iterator;
  typedef CGAL_Triangulation_facet_iterator_3<GT,Tds>   Facet_iterator;
  typedef CGAL_Triangulation_edge_iterator_3<GT,Tds>   Edge_iterator;
  typedef CGAL_Triangulation_vertex_iterator_3<GT,Tds> Vertex_iterator;

  enum Locate_type {
    VERTEX=0, 
    EDGE, //1
    FACET, //2
    CELL, //3
    OUTSIDE_CONVEX_HULL, //4
    OUTSIDE_AFFINE_HULL };//5

private:
  Tds _Triangulation_data_structure_3;
  GT  _gt;
  Vertex_handle infinite; //infinite vertex
  
//   void set_finite_vertex(const Vertex_handle&  v)
//   {
//   }

//   void set_infinite_vertex(const Vertex_handle & v)
//   {
//   }

  void set_number_of_vertices(int n) 
    { _Triangulation_data_structure_3.set_number_of_vertices(n+1); }

  void init_Triangulation_data_structure_3()
    {
      infinite = new Vertex(Point(500,500,500)); // ?? debug
      _Triangulation_data_structure_3.insert_outside_affine_hull(&(*infinite));
      // forces the compiler to instanciate CGAL_debug :
      CGAL_debug( infinite ); 
      CGAL_debug( Cell_handle() );
    }
  
public:

// CONSTRUCTORS
  CGAL_Triangulation_3()
    : _Triangulation_data_structure_3(), _gt()
  {
    init_Triangulation_data_structure_3();
  }

  CGAL_Triangulation_3(const GT & gt) 
    : _Triangulation_data_structure_3(), _gt(gt)
  {
    init_Triangulation_data_structure_3();
  }

  // copy constructor duplicates vertices and cells
  CGAL_Triangulation_3(const CGAL_Triangulation_3<GT,Tds> & tr)
    : _Triangulation_data_structure_3(tr._Triangulation_data_structure_3), _gt(tr._gt), infinite(tr.infinite)
  {}

  // pour la demo
  CGAL_Triangulation_3(const Point & p0,
		     const Point & p1,
		     const Point & p2,
		     const Point & p3)
    : _Triangulation_data_structure_3(), _gt()
    {
      init_Triangulation_data_structure_3();
      insert_outside_affine_hull(p0);
      insert_outside_affine_hull(p1);
      insert_outside_affine_hull(p2);
      insert_outside_affine_hull(p3);
    }

//   CGAL_Triangulation_3(Tds Triangulation_data_structure_3, const GT & gt=GT())
//     : _Triangulation_data_structure_3(),_gt(gt)
//   {
//     _Triangulation_data_structure_3.swap(Triangulation_data_structure_3);
//   }

  // DESTRUCTOR

  ~CGAL_Triangulation_3()
  {
    clear();
  }

  void clear()
  {
    _Triangulation_data_structure_3.clear();
    infinite.Delete();
  }

 //ACCESS FUNCTIONS
  inline
  const GT & geom_traits() const 
    { return _gt;}
  
  inline
  const Tds & Triangulation_data_structure_3() const 
    { return _Triangulation_data_structure_3;}
  
  inline
  int dimension() const 
    { return _Triangulation_data_structure_3.dimension();}

  inline
  int number_of_vertices() const // number of finite vertices
    {return _Triangulation_data_structure_3.number_of_vertices()-1;}

  inline
  Vertex_handle infinite_vertex() const
  {
    return infinite;
  }
   
  inline
  Cell_handle infinite_cell() const
  {
    //    CGAL_triangulation_precondition(infinite_vertex()->cell()->
    //				    has_vertex(infinite_vertex()));
    return infinite_vertex()->cell();
  }

  
   // ASSIGNMENT
  CGAL_Triangulation_3 & operator=(const CGAL_Triangulation_3 & tr)
  {
    clear();
    _Triangulation_data_structure_3 = tr._Triangulation_data_structure_3;
    _gt = tr._gt;
    infinite = tr.infinite;
    return *this;
  }

  // HELPING FUNCTIONS
   
  void copy_triangulation(const CGAL_Triangulation_3 & tr)
  {
    clear();
    _gt = tr._gt;
    _Triangulation_data_structure_3.copy_Triangulation_data_structure_3(tr._Triangulation_data_structure_3);
    infinite = tr.infinite;
  }

  void swap(CGAL_Triangulation_3 &tr)
  {
    GT t = geom_traits();
    _gt = tr.geom_traits();
    tr._gt = t; 

    Vertex* inf = infinite_vertex();
    infinite = tr.infinite_vertex();
    tr.infinite = inf;

    _Triangulation_data_structure_3.swap(tr._Triangulation_data_structure_3);
  }

  // CHECKING
  bool is_valid(bool verbose = false, int level = 0) const
  {
    if ( ! _Triangulation_data_structure_3.is_valid(verbose,level) ) {
      if (verbose) { cerr << "invalid data structure" << endl; }
      CGAL_triangulation_assertion(false); return false;
    }

    if (dimension() == 3 ) {
      Cell_iterator it;
      for ( it = finite_cells_begin(); it != cells_end(); ++it ) {
	if ( geom_traits().orientation(it->vertex(0)->point(),
			      it->vertex(1)->point(),
			      it->vertex(2)->point(),
			      it->vertex(3)->point()) != CGAL_LEFTTURN ) {
	  if (verbose) { cerr << "badly oriented cell" << endl; }
	  CGAL_triangulation_assertion(false); return false;
	}
      }
    }
    if (verbose) { cerr << "valid triangulation" << endl;}
    return true;
  }

 // GEOMETRIC ACCESS FUNCTIONS
  
  Tetrahedron tetrahedron(const Cell_handle & c) const
  {
    CGAL_triangulation_precondition( dimension() == 3 );
    CGAL_triangulation_precondition( ! is_infinite(c) );
    return Tetrahedron(c->vertex(0)->point(),
		       c->vertex(1)->point(),
		       c->vertex(2)->point(),
		       c->vertex(3)->point());
  }

  Triangle triangle(const Cell_handle & c, int i) const
  { 
    switch ( dimension() ) {
    case 3:
      {
	CGAL_triangulation_precondition
	  ( i == 0 || i == 1 || i == 2 || i == 3 );
	break;
      }
    case 2:
      {
	CGAL_triangulation_precondition( i == 3 );
	break;
      }
    default:
	CGAL_triangulation_assertion( false );
	// return ?
    }
    CGAL_triangulation_precondition( ! is_infinite(make_pair(c,i)) );
    switch (i) {
    case 0:
      return Triangle(c->vertex(1)->point(),
		      c->vertex(2)->point(),
		      c->vertex(3)->point());
      break;
    case 1:
      return Triangle(c->vertex(0)->point(),
		      c->vertex(2)->point(),
		      c->vertex(3)->point());
    case 2:
      return Triangle(c->vertex(0)->point(),
		      c->vertex(1)->point(),
		      c->vertex(3)->point());
    case 3:
      return Triangle(c->vertex(0)->point(),
		      c->vertex(1)->point(),
		      c->vertex(2)->point());
    default:
      {
	// impossible
	CGAL_triangulation_assertion( false );
	// to avoid warning at compile time :
	return Triangle(c->vertex(1)->point(),
			c->vertex(2)->point(),
			c->vertex(3)->point());
      }
    }
  }

  Triangle triangle(const Facet & f) const
  {
    return triangle(f.first, f.second);
  }

  Segment segment(const Cell_handle & c, int i, int j) const
  {
    CGAL_triangulation_precondition( i != j );
    switch ( dimension() ) {
    case 3:
      {
	CGAL_triangulation_precondition
	  ( ( i == 0 || i == 1 || i == 2 || i == 3 ) &&
	    ( j == 0 || j == 1 || j == 2 || j == 3 ) );
	break;
      }
    case 2:
      {
	CGAL_triangulation_precondition
	  ( ( i == 0 || i == 1 || i == 2 ) &&
	    ( j == 0 || j == 1 || j == 2 ) );
	break;
      }
    case 1:
      {
	CGAL_triangulation_precondition( ( i == 0 || i == 1 ) &&
					 ( j == 0 || j == 1 ) );
	break;
      }
    default:
	CGAL_triangulation_assertion( false );
	// return ?
    }
    CGAL_triangulation_precondition( ! is_infinite(CGAL_make_triple(c,i,j)) );
    return Segment( c->vertex(i)->point(), c->vertex(j)->point() );
  }

  Segment segment(Edge e) const
  {
    return segment(e.first,e.second,e.third);
  }

//   Segment segment(const Edge_circulator& ec) const
//   {
//     return segment(*ec);
//   }
    
//   Segment segment(const Edge_iterator& ei) const
//   {
//     return segment(*ei);
//   }

  // TEST IF INFINITE FEATURES

  bool is_infinite(const Cell_handle & c) const 
    {
      CGAL_triangulation_precondition( dimension() == 3 );
      return c->has_vertex(infinite_vertex());
    }

  bool is_infinite(const Vertex_handle & v) const 
    {
      return v == infinite_vertex();
    }

  bool is_infinite(const Cell_handle c, int i) const 
    {
      switch ( dimension() ) {
      case 3:
	{
	  CGAL_triangulation_precondition
	    ( i == 0 || i == 1 || i == 2 || i == 3 );
	  break;
	}
      case 2:
	{
	  CGAL_triangulation_precondition( i == 3 );
	  break;
	}
      default:
	CGAL_triangulation_assertion( false );
	// return ?
      }
      switch (i) {
      case 0:
	return ( is_infinite(c->vertex(1)) ||
		 is_infinite(c->vertex(2)) ||
		 is_infinite(c->vertex(3)) );
	break;
      case 1:
	return ( is_infinite(c->vertex(0)) ||
		 is_infinite(c->vertex(2)) ||
		 is_infinite(c->vertex(3)) );
      case 2:
	return ( is_infinite(c->vertex(0)) ||
		 is_infinite(c->vertex(1)) ||
		 is_infinite(c->vertex(3)) );
      case 3:
	return ( is_infinite(c->vertex(0)) ||
		 is_infinite(c->vertex(1)) ||
		 is_infinite(c->vertex(2)) );
      }
      // we never get here
      CGAL_triangulation_precondition( false );
      return false;
    }

  bool is_infinite(Facet f) const 
    {
      return is_infinite(f.first,f.second);
    }

  bool is_infinite(const Cell_handle c, int i, int j) const 
    { 
      CGAL_triangulation_precondition( ! (i == j) );
      switch ( dimension() ) {
      case 3:
	{
	  CGAL_triangulation_precondition
	    ( ( i == 0 || i == 1 || i == 2 || i == 3 ) &&
	      ( j == 0 || j == 1 || j == 2 || j == 3 ) );
	  break;
	}
      case 2:
	{
	  CGAL_triangulation_precondition
	    ( ( i == 0 || i == 1 || i == 2 ) &&
	      ( j == 0 || j == 1 || j == 2 ) );
	  break;
	}
      case 1:
	{
	  CGAL_triangulation_precondition( ( i == 0 || i == 1 ) &&
					   ( j == 0 || j == 1 ) );
	  break;
	}
      default:
	CGAL_triangulation_assertion( false );
	// return 
      }
      return ( is_infinite( c->vertex(i) ) ||
	       is_infinite( c->vertex(j) ) );
    }

  bool is_infinite(Edge e) const
    {
      return is_infinite(e.first,e.second,e.third);
    }

//   bool is_infinite(const Edge_circulator& ec) const 
//     {
//       return is_infinite(*ec);
//     }

//   bool is_infinite(const Edge_iterator& ei) const 
//     {
//       return is_infinite(*ei);
//     }

  //INSERTION 

  Vertex_handle
  insert_in_cell(const Point & p, Cell_handle c)
  {
    CGAL_triangulation_precondition( dimension() == 3 );
    Locate_type lt;
    int i,j;
    CGAL_triangulation_precondition
      ( side_of_tetrahedron( p, 
			     c->vertex(0)->point(),
			     c->vertex(1)->point(),
			     c->vertex(2)->point(),
			     c->vertex(3)->point(),
			     lt,i,j ) == CGAL_ON_BOUNDED_SIDE );
    Vertex_handle v = new Vertex(p);
    _Triangulation_data_structure_3.insert_in_cell( &(*v), &(*c) );
    return v;
  }

  Vertex_handle
  insert_in_facet(const Point & p, Cell_handle c, int i)
  {
    switch ( dimension() ) {
    case 3:
      {
	CGAL_triangulation_precondition
	  ( i == 0 || i == 1 || i == 2 || i == 3 );
	break;
      }
    case 2:
      {
	CGAL_triangulation_precondition( i == 3 );
	break;
      }
    default:
      CGAL_triangulation_assertion( false );
      // return ?
    }
    Locate_type lt;
    int li,lj; 
    CGAL_triangulation_precondition
      ( (geom_traits().orientation( p, 
			   c->vertex((i+1)&3)->point(),
			   c->vertex((i+2)&3)->point(),
			   c->vertex((i+3)&3)->point() ) == CGAL_COPLANAR)
	&& 
	(side_of_triangle(p, 
			  c->vertex((i+1)&3)->point(),
			  c->vertex((i+2)&3)->point(),
			  c->vertex((i+3)&3)->point(),
			  lt,li,lj) == CGAL_ON_BOUNDED_SIDE)
	);
    Vertex_handle v = new Vertex(p);
    _Triangulation_data_structure_3.insert_in_facet(&(*v), &(*c), i);
    return v;
  }

  Vertex_handle
  insert_in_facet(const Point & p, Facet f)
  {
    return insert_in_facet(p, f.first,f.second);
  }

  Vertex_handle
  insert_in_edge(const Point & p, Cell_handle c, int i, int j)
  {
    CGAL_triangulation_precondition( ! (i == j) );
    switch ( dimension() ) {
    case 3:
      {
	CGAL_triangulation_precondition
	  ( ( i == 0 || i == 1 || i == 2 || i == 3 ) &&
	    ( j == 0 || j == 1 || j == 2 || j == 3 ) );
	CGAL_triangulation_precondition( ! is_infinite(c,i,j) );
   	break;
      }
    case 2:
      {
	CGAL_triangulation_precondition
	  ( ( i == 0 || i == 1 || i == 2 ) &&
	    ( j == 0 || j == 1 || j == 2 ) );
	CGAL_triangulation_precondition( ! is_infinite(c,i,j) );
   	break;
      }
    case 1:
      {
	CGAL_triangulation_precondition( ( i == 0 || i == 1 ) &&
					 ( j == 0 || j == 1 ) );
	break;
      }
    default:
      CGAL_triangulation_assertion( false );
      // return 
    }
    Locate_type lt;
    int li;
    CGAL_triangulation_precondition
      ( geom_traits().collinear( c->vertex(i)->point(),
			p,
			c->vertex(j)->point() )
	&&
	( side_of_segment( p,
			c->vertex(i)->point(),
			c->vertex(j)->point(),
			lt,li ) == CGAL_ON_BOUNDED_SIDE )
      );
    Vertex_handle v = new Vertex(p);
    _Triangulation_data_structure_3.insert_in_edge(&(*v), &(*c), i, j);
    return v;
  }
  
  Vertex_handle
  insert_in_edge(const Point & p, Edge e)
  {
    return insert_in_edge(p, e.first,e.second,e.third);
  }
  
  Vertex_handle
  insert_outside_convex_hull(const Point & p, Cell_handle c, 
			     int li, int lj=0)
    // c is a finite cell whose facet li lies on the convex hull boundary
    // and separates p from the triangulation
    // p is strictly outside the convex hull
    // in dimension 2, edge li,lj separates p from the triangulation
    // in dimension 1, vertex li separates p from the triangulation
    // dimension 0 not allowed, use outside-affine-hull
  {
    CGAL_triangulation_precondition( ! c->has_vertex(infinite_vertex()) );
    // is_infinite not used because of precondition on dimension
    CGAL_triangulation_precondition( dimension() > 0 );

    switch ( dimension() ) {
    case 1:
      {
	// p lies in the infinite edge neighboring c 
	// on the other side of li
	return insert_in_edge(p,c->neighbor(1-li),0,1);
      }
    case 2:
      {
	Cell_handle n = c->neighbor(3-li-lj);
	// n contains p and is infinite

	Vertex_handle v = new Vertex(p);
	set_number_of_vertices(number_of_vertices()+1);

	Locate_type loc;
	int i, j;

	// traversal in the first one of the two directions
	// of the infinite cells containing p
	// updating of the triangulation at the same time
	// by replacing the infinite vertex by v
	Cell_handle cur = n;
	Cell_handle prev = n->neighbor( ccw(n->index(infinite)) );

	while ( side_of_facet( p, cur, 3, loc, i, j ) 
		// in dimension 2, cur has only one facet numbered 3
		== CGAL_ON_BOUNDED_SIDE ) {
	    // loc must be == CELL since p supposed to be
	    // strictly outside the convex hull
	    cur->set_vertex( cur->index(infinite), v );
	    prev = cur;
	    cur = cur->neighbor( cw(cur->index(v)) ) ;
	}
	
	// creation of an infinite facet "at the end" of the sequence
	// of infinite facets containing p
	Cell_handle cnew 
	  = new Cell( _Triangulation_data_structure_3,
		      prev->vertex(ccw(prev->index(v))), v,  
		      infinite_vertex(), NULL,
		      NULL, cur, prev, NULL);
	// neighbor 0 will be set to dnew later
	prev->set_neighbor(prev->index(cur), cnew);
	cur->set_neighbor(cur->index(prev),cnew);

	// traversal in the opposite direction, and same operations
	// starts from the neighbor of n (n already modified)
	prev = n;
	cur = n->neighbor( ccw(n->index(v)) );
	
	while ( side_of_facet( p, cur, 3, loc, i, j ) 
		== CGAL_ON_BOUNDED_SIDE ) {
	  cur->set_vertex( cur->index(infinite), v );
	  prev = cur;
	  cur = cur->neighbor( ccw(cur->index(v)) ) ;
	}
	
	Cell_handle dnew 
	  = new Cell( _Triangulation_data_structure_3,
		      v, prev->vertex(cw(prev->index(v))), 
		      infinite_vertex(), NULL,
		      cur, cnew, prev, NULL);
	prev->set_neighbor(prev->index(cur), dnew);
	cur->set_neighbor(cur->index(prev),dnew);
	cnew->set_neighbor(0,dnew); // correction for cnew

	infinite->set_cell(dnew);			     
	v->set_cell( n );
	return v;
      }
    case 3:
      {
	Cell_handle n = c->neighbor(li);
	// n is an infinite cell containing p

	Vertex_handle v = new Vertex(p);
	v->set_cell( n );

	set_number_of_vertices(number_of_vertices()+1);

	link( v, hat(v,n) );
	// infinite->set_cell is done by link

	return v;
      }
    }
    // to avoid warning with eg++
    return NULL;
  }

private:
Cell_handle
hat(Vertex_handle v, Cell_handle c)
  // recursive traversal of the set of facets of the convex hull
  // that are visible from v
  // v replaces infinite_vertex in these cells
  // on the boundary, new cells with vertices v, infinite_vertex
  // and the two vertices of an edge of the boumdary are created
  // returns a cell inside the "hat", having a facet on its boundary
{
  static Cell_handle bound;

  int inf = c->index(infinite_vertex());
  c->set_vertex( inf , v );

  Cell_handle cni, cnew;
  Locate_type loc;
  int li,lj;

  int i, i1, i2;
  for ( i=0; i<4; i++ ) {
    if ( i!= inf ) {
      cni = c->neighbor(i);
      if ( ! cni->has_vertex( v ) ) {
	if ( side_of_cell( v->point(), cni, loc, li, lj )
	     == CGAL_ON_BOUNDED_SIDE ) {
	  hat( v, cni );
	}
	else { // we are on the boundary of the set of facets of the
	  // convex hull that are visible from v
	  i1 = nextposaroundij(i,inf);
	  i2 = 6-i-i1-inf;
	  cnew = new Cell( _Triangulation_data_structure_3,
			   c->vertex(i1), c->vertex(i2), 
			   v, infinite_vertex(),
			   NULL, NULL, cni, c );
	  c->set_neighbor(i,cnew);
	  cni->set_neighbor( cni->index(c), cnew );

	  bound = c;
	}
      }
    }// no else
  } // for
  return bound;
} // hat
void
link(Vertex_handle v, Cell_handle c)
  // c belongs to the hat of v anad has a fecet on its boundary
  // traverses the boundary of the hat and finds adjacencies
  // traversal is done counterclockwise as seen from v
{
  // finds a facet ib of c on the boundary of the hat
  int iv = c->index(v);
  int ib;
  for ( ib=0; ib<4; ib++ ) {
    if ( ib != iv ) {
      if ( c->neighbor(ib)->has_vertex(infinite) ) {
	break;
      }
    }
    else {
      ++ib;
    }
  }
  
  infinite->set_cell(c->neighbor(ib));

  Cell_handle bound = c;
  int i = ib;
  int next;
  Vertex_handle v1;

  do {
    iv = bound->index(v);
    // indices of the vertices != v of bound on the boundary of the hat
    // such that (i,i1,i2,iv) positive
    int i1 = nextposaroundij(i,iv);
    int i2 = 6-i-i1-iv;

    // looking for the neighbor i2 of bound :
    // we turn around i1 until we reach the boundary of the hat
    v1 = bound->vertex(i1);

    Cell_handle cur = bound;

    next = nextposaroundij(i1,iv);
    while ( ! cur->neighbor(next)->has_vertex(infinite) ) {
      cur = cur->neighbor(next);
      next = nextposaroundij(cur->index(v1),cur->index(v));
    }
    Cell_handle courant = bound->neighbor(i);
    Cell_handle trouve = cur->neighbor(next);
    courant->set_neighbor( courant->index(bound->vertex(i2)), trouve);
    trouve->set_neighbor( 6 - trouve->index(v) - 
			  trouve->index(infinite) -
			  trouve->index(v1), courant );
    bound = cur;
    i = next;
  } while ( ( bound != c ) || ( i != ib ) );
  // c may have two facets on the boundary of the hat
  // test bound != c is not enough, we must test whether
  // facet ib of c has been treated
}// end link
public:

  Vertex_handle
  insert_outside_affine_hull(const Point & p)
  {
    CGAL_triangulation_precondition( dimension() < 3 );
    bool reorient;
    switch ( dimension() ) {
    case 1:
      {
	Cell_handle c = infinite_cell();
	Cell_handle n = c->neighbor(c->index(infinite_vertex()));
	CGAL_triangulation_precondition
	  ( ! geom_traits().collinear(p,
			     n->vertex(0)->point(),
			     n->vertex(1)->point()) );
	// no reorientation : the first non-collinear point determines
	// the orientation of the plane
	reorient = false;
	break;
      }
    case 2:
      {
	Cell_handle c = infinite_cell();
	Cell_handle n = c->neighbor(c->index(infinite_vertex()));
	CGAL_triangulation_precondition
	  ( geom_traits().orientation(n->vertex(0)->point(),
			     n->vertex(1)->point(),
			     n->vertex(2)->point(),
			     p) != CGAL_COPLANAR );
	reorient = ( geom_traits().orientation( n->vertex(0)->point(),
				       n->vertex(1)->point(),
				       n->vertex(2)->point(),
				       p ) == CGAL_NEGATIVE );
	break;
      }
    default:
      reorient = false;
      break;
    }
    Vertex_handle v = new Vertex(p);
    _Triangulation_data_structure_3.insert_outside_affine_hull(&(*v), &(*infinite_vertex()), reorient);
    return v;
  }

  Vertex_handle insert(const Point & p )
  {
    Locate_type lt;
    int li, lj;
    Cell_handle c = locate( p, lt, li, lj);
    switch (lt) {
    case VERTEX:
      return c->vertex(li);
    case EDGE:
      return insert_in_edge(p, c, li, lj);
    case FACET:
      return insert_in_facet(p, c, li);
    case CELL:
      return insert_in_cell(p, c);
    case OUTSIDE_CONVEX_HULL:
      return insert_outside_convex_hull(p, c, li, lj);
    case OUTSIDE_AFFINE_HULL:
      return insert_outside_affine_hull(p);
    }
    // cannot happen, only to avoid warning with eg++
    return insert_in_edge(p, c, li, lj);
  }

#ifndef CGAL_CFG_NO_MEMBER_TEMPLATES
  template < class InputIterator >
  int insert(InputIterator first, InputIterator last)
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
  int insert(list<Point>::const_iterator first,
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
  int insert(vector<Point>::const_iterator first,
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
  int insert(istream_iterator<Point, ptrdiff_t> first,
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

  Cell_handle
  locate(const Point & p,
	 Locate_type & lt,
	 int & li,
	 int & lj) const
    // returns a cell p lies in if there is one
    // if lt == OUTSIDE_CONVEX_HULL, returns a finite cell, and li is the
    // index of a facet separating p from the rest of the triangulation
    // in dimension 2 :
    // returns a facet (Cell_handle,li) if lt == FACET
    // returns an edge (Cell_handle,li,lj) if lt == EDGE
    // returns a vertex (Cell_handle,li) if lt == VERTEX
    // if lt == OUTSIDE_CONVEX_HULL, li, lj give the edge of c
    // separating p from the rest of the triangulation
  {
    bool notfound = true;
    switch (dimension()) {
    case 3:
      {
	Cell_iterator cit = all_cells_begin();
	Cell_iterator cdone = cells_end();
	do {
	  //	  CGAL_Bounded_side side = 
	  if ( side_of_cell( p, &(*cit), lt, li, lj )
	       != CGAL_ON_UNBOUNDED_SIDE ) {
	    notfound = false;
	  }
	  else {
	    ++cit;
	  }
	} while ( cit != cdone && notfound );
	if ( notfound ) {
	  // cannot happen : there must be a cell (finite or not)
	  // containing p
	  CGAL_triangulation_assertion(false);
	  return NULL;
	}
	if ( is_infinite(&(*cit)) ) {
	  switch ( lt ) {
	  case CELL:
	    {
	      // returns the finite cell sharing the finite facet of cit
	      Cell_handle n = cit->neighbor(cit->index(infinite));
	      lt = OUTSIDE_CONVEX_HULL;
	      li = n->index(&(*cit));
	      return n;
	    }
	  case FACET:
	    {
	      // returns the finite cell sharing the finite facet of cit
	      Cell_handle n = cit->neighbor(cit->index(infinite));
	      li = n->index(&(*cit));
	      return n;
	    }
	  case EDGE:
	    {
	      // returns the finite cell sharing the finite facet of cit
	      Cell_handle n = cit->neighbor(cit->index(infinite));
	      cerr << cit->vertex(li)->point() << endl
		   << cit->vertex(lj)->point() << endl;
	      li = n->index(cit->vertex(li));
	      lj = n->index(cit->vertex(lj));
	      return n;
	    }
	  case VERTEX:
	    {
	      // returns the finite cell sharing the finite facet of cit
	      Cell_handle n = cit->neighbor(cit->index(infinite));
	      li = n->index(cit->vertex(li));
	      return n;
	    }
	  default:
	    {
	      CGAL_triangulation_assertion(false);
	      return NULL;
	    }
	  }
	}
	else { // finite cell
	  return &(*cit);
	}
	break;
      }
    case 2:
      {
	//first tests whether p is coplanar with the current triangulation
	Facet_iterator finite_fit = finite_facets_begin();
	if ( geom_traits().orientation( p, 
			       (*finite_fit).first->vertex(0)->point(),
			       (*finite_fit).first->vertex(1)->point(),
			       (*finite_fit).first->vertex(2)->point() ) 
	     != CGAL_DEGENERATE ) {
	  lt = OUTSIDE_AFFINE_HULL;
	  // li = 3; // only one facet : any cell is degenerate in dimension 2
	  return (*finite_fit).first;
	}
	// if p is coplanar, location in the triangulation
	Facet_iterator fit = all_facets_begin();
	Facet_iterator fdone = facets_end();
	do {
	  //	  CGAL_Bounded_side side = 
	  if ( side_of_facet( p, *fit, lt, li, lj )
	       != CGAL_ON_UNBOUNDED_SIDE ) {
	    notfound = false;
	  }
	  else {
	    ++fit;
	  }
	} while ( fit != fdone && notfound );
	if ( notfound ) {
	  CGAL_triangulation_assertion(false);
	  return NULL;
	}
	if ( is_infinite(*fit) ) {
	      // returns the finite facet sharing the finite edge of (*fit)
	  switch ( lt ) {
	  case FACET:
	    {
	      Cell_handle 
		n = (*fit).first->neighbor((*fit).first->index(infinite));
	      li = n->index( (*fit).first->vertex
			     (cw((*fit).first->index(infinite))) );
	      lj = n->index( (*fit).first->vertex
			     (ccw((*fit).first->index(infinite))) );
	      lt = OUTSIDE_CONVEX_HULL;
	      return n;
	    }
	  case EDGE:
	    {
	      Cell_handle 
		n = (*fit).first->neighbor((*fit).first->index(infinite));
	      li = n->index((*fit).first->vertex(li));
	      lj = n->index((*fit).first->vertex(lj));
	      return n;
	    }
	  case VERTEX:
	    {
	      Cell_handle 
		n = (*fit).first->neighbor((*fit).first->index(infinite));
	      li = n->index((*fit).first->vertex(li));
	      return n;
	    }
	  default:
	    {
	      // cannot happen, only to avoid warning with eg++
	      return (*fit).first;
	    }
	  }
	}
	else { // finite facet
	  if ( lt == FACET ) {
	    li = 3;
	  }
	  // case vertex or edge : li and lj already correct
	  // because the index of the vertices in the facet is the same as the
	  // index in the underlying degenerate cell
	  return (*fit).first;
	}
	break;
      }
    case 1:
      {
	//first tests whether p is collinear with the current triangulation
	Edge_iterator finite_eit = finite_edges_begin();
	if ( ! geom_traits().collinear(p,
			      (*finite_eit).first->vertex(0)->point(),
			      (*finite_eit).first->vertex(1)->point()) ) {
	  lt = OUTSIDE_AFFINE_HULL;
	  return (*finite_eit).first;
	}
	// if p is collinear, location :
	Edge_iterator eit = all_edges_begin();
	Edge_iterator edone = edges_end();
	do {
	  if ( side_of_edge( p, *eit, lt, li )
	       != CGAL_ON_UNBOUNDED_SIDE ) {
	    notfound = false;
	  }
	  else {
	    ++eit;
	  }
	} while ( eit != edone && notfound );
	if ( notfound ) {
	  CGAL_triangulation_assertion(false);
	  return NULL;
	}
	if ( is_infinite(*eit) ) {
	  // returns the finite edge sharing the finite vertex of *eit
	  switch ( lt ) {
	  case EDGE:
	    {
	      Cell_handle 
		n = (*eit).first->neighbor((*eit).first->index(infinite));
	      li = n->index( (*eit).first->vertex
			     ( 1 - (*eit).first->index(infinite) ) );
	      lj = 1-li;
	      lt = OUTSIDE_CONVEX_HULL;
	      return n;
	    }
	  case VERTEX:
	    {
	      Cell_handle 
		n = (*eit).first->neighbor((*eit).first->index(infinite));
	      li = n->index((*eit).first->vertex(li));
	      return n;
	    }
	  default :
	    {
	      return (*eit).first;
	    }
	  }
	}
	else { // finite edge
	  li = (*eit).second;
	  lj = (*eit).third;
	  return (*eit).first;
	}
	break;
      }
    case 0:
      {
	Vertex_iterator vit = finite_vertices_begin();
	if ( p != vit->point() ) {
	  lt = OUTSIDE_AFFINE_HULL;
	}
	else {
	  lt = VERTEX;
	  li = 0;
	}
	return vit->cell();
	break;
      }
    case -1:
      {
	lt = OUTSIDE_AFFINE_HULL;
	return NULL;
      }
    default:
      {
	CGAL_triangulation_assertion(false);
	return NULL;
      }
    }
    // impossible. only to avoid warning
    CGAL_triangulation_assertion(false);
    return NULL;
  }

  inline Cell_handle
  locate(const Point & p) const
    {
      Locate_type lt;
      int li, lj;
      return locate(p, lt, li, lj);
    }

  //TRAVERSING : ITERATORS AND CIRCULATORS
  Cell_iterator finite_cells_begin() const
  {
    CGAL_Triangulation_3<GT, Tds>* ncthis = (CGAL_Triangulation_3<GT, Tds> *)this;
    return Cell_iterator(ncthis, false); // false means "without infinite cells"
  }
  Cell_iterator all_cells_begin() const
  {
    CGAL_Triangulation_3<GT, Tds>* ncthis = (CGAL_Triangulation_3<GT, Tds> *)this;
    return Cell_iterator(ncthis, true); // true means "with infinite cells"
  }
  Cell_iterator cells_end() const
  {
    CGAL_Triangulation_3<GT, Tds>* ncthis = (CGAL_Triangulation_3<GT, Tds> *)this;
    return Cell_iterator(ncthis); // not second argument -> past-end
  }

  Vertex_iterator finite_vertices_begin() const
  {
    CGAL_Triangulation_3<GT, Tds>* ncthis = (CGAL_Triangulation_3<GT, Tds>*)this;
    return Vertex_iterator(ncthis, false);
  }
  Vertex_iterator all_vertices_begin() const
  {
    CGAL_Triangulation_3<GT, Tds>* ncthis = (CGAL_Triangulation_3<GT, Tds>*)this;
    return Vertex_iterator(ncthis, true);
  }
  Vertex_iterator vertices_end() const
  {
    CGAL_Triangulation_3<GT, Tds>* ncthis = (CGAL_Triangulation_3<GT, Tds>*)this;
    return Vertex_iterator(ncthis);
  }

  Edge_iterator finite_edges_begin() const
  {
    CGAL_Triangulation_3<GT, Tds>* ncthis = (CGAL_Triangulation_3<GT, Tds>*)this;
    return Edge_iterator(ncthis, false);
  }
  Edge_iterator all_edges_begin() const
  {
    CGAL_Triangulation_3<GT, Tds>* ncthis = (CGAL_Triangulation_3<GT, Tds>*)this;
    return Edge_iterator(ncthis, true);
  }
  Edge_iterator edges_end() const
  {
    CGAL_Triangulation_3<GT, Tds>* ncthis = (CGAL_Triangulation_3<GT, Tds>*)this;
    return Edge_iterator(ncthis);
  }

  Facet_iterator finite_facets_begin() const
  {
    CGAL_Triangulation_3<GT, Tds>* ncthis = (CGAL_Triangulation_3<GT, Tds>*)this;
    return Facet_iterator(ncthis, false);
  }
  Facet_iterator all_facets_begin() const
  {
    CGAL_Triangulation_3<GT, Tds>* ncthis = (CGAL_Triangulation_3<GT, Tds>*)this;
    return Facet_iterator(ncthis, true);
  }
  Facet_iterator facets_end() const
  {
    CGAL_Triangulation_3<GT, Tds>* ncthis = (CGAL_Triangulation_3<GT, Tds>*)this;
    return Facet_iterator(ncthis);
  }

  // PREDICATES ON POINTS ``TEMPLATED'' by the geom traits

  CGAL_Bounded_side
  side_of_segment(const Point & p,
		  const Point & p0, 
		  const Point & p1,
		  Locate_type & lt, int & i ) const
    // p0, p1 supposed to be different
    // p supposed to be collinear to p0, p1
    // returns :
    // CGAL_ON_BOUNDED_SIDE if p lies strictly inside the edge
    // CGAL_ON_BOUNDARY if p equals p0 or p1
    // CGAL_ON_UNBOUNDED_SIDE if p lies strictly outside the edge
    {
      CGAL_triangulation_precondition
	( ! geom_traits().equal(p0,p1) );
      CGAL_triangulation_precondition
	( geom_traits().collinear(p,p0,p1) );
      
      CGAL_Comparison_result c = geom_traits().compare_x(p0,p1);
      CGAL_Comparison_result c0;
      CGAL_Comparison_result c1;

      if ( c == CGAL_EQUAL ) {
	c = geom_traits().compare_y(p0,p1);
	if ( c == CGAL_EQUAL ) {
	  c = geom_traits().compare_z(p0,p1);
	  c0 = geom_traits().compare_z(p0,p);
	  c1 = geom_traits().compare_z(p,p1);
	}
	else {
	  c0 = geom_traits().compare_y(p0,p);
	  c1 = geom_traits().compare_y(p,p1);
	}
      }
      else {
	c0 = geom_traits().compare_x(p0,p);
	c1 = geom_traits().compare_x(p,p1);
      }
      
      //      if ( (c0 == CGAL_SMALLER) && (c1 == CGAL_SMALLER) ) {
      if ( c0 == c1 ) {
	lt = EDGE;
	return CGAL_ON_BOUNDED_SIDE;
      }
      if (c0 == CGAL_EQUAL) {
	lt = VERTEX;
	i = 0;
	return CGAL_ON_BOUNDARY;
      }
      if (c1 == CGAL_EQUAL) {
	lt = VERTEX;
	i = 1;
	return CGAL_ON_BOUNDARY;
      }
      lt = OUTSIDE_CONVEX_HULL;
      return CGAL_ON_UNBOUNDED_SIDE;
    }

  CGAL_Bounded_side
  side_of_edge(const Point & p,
		Cell_handle c,
		int i, int j,
		Locate_type & lt, int & li) const
    // supposes dimension 1 otherwise does not work for infinite edges
    // returns :
    // CGAL_ON_BOUNDED_SIDE if p inside the edge 
    // (for an infinite edge this means that p lies in the half line
    // defined by the vertex)
    // CGAL_ON_BOUNDARY if p equals one of the vertices
    // CGAL_ON_UNBOUNDED_SIDE if p lies outside the edge
    // (for an infinite edge this means that p lies on the other half line)
    // lt has a meaning when CGAL_ON_BOUNDED_SIDE and CGAL_ON_BOUNDARY
    // li refer to indices in the cell c 
    {//side_of_edge
      CGAL_triangulation_precondition( dimension() == 1 );
      if ( ! is_infinite(c,i,j) ) {
	return side_of_segment(p,
			       c->vertex(i)->point(),
			       c->vertex(j)->point(),
			       lt, li);
      }
      else { // infinite edge
	if ( geom_traits().equal( p, c->vertex(i)->point() ) ) {
	  lt = VERTEX;
	  li = i;
	  return CGAL_ON_BOUNDARY;
	}
	if ( geom_traits().equal( p, c->vertex(j)->point() ) ) {
	  lt = VERTEX;
	  li = j;
	  return CGAL_ON_BOUNDARY;
	}
	// does not work in dimension > 2
	int inf = c->index(infinite);
	Cell_handle n = c->neighbor(inf);
	int i_e = n->index(c);
	// we know that n is finite
	Vertex_handle
	  v0 = n->vertex(0),
	  v1 = n->vertex(1);
	CGAL_Comparison_result c = geom_traits().compare_x(v0->point(),
						  v1->point());
	CGAL_Comparison_result cp;
	if ( c == CGAL_EQUAL ) {
	  c = geom_traits().compare_y(v0->point(),
			     v1->point());
	  if ( i_e == 0 ) {
	    cp = geom_traits().compare_y( v1->point(), p );
	  }
	  else {
	    cp = geom_traits().compare_y( p, v0->point() );
	  }
	}
	else {
	  if ( i_e == 0 ) {
	    cp = geom_traits().compare_x( v1->point(), p );
	  }
	  else {
	    cp = geom_traits().compare_x( p, v0->point() );
	  }
	}
	if ( c == cp ) {
	  // p lies on the same side of n as infinite
	  lt = EDGE;
	  return CGAL_ON_BOUNDED_SIDE;
	}
	else {
	  return CGAL_ON_UNBOUNDED_SIDE;
	}
      }
    }

  CGAL_Bounded_side
  side_of_edge(const Point & p,
	       Edge e,
	       Locate_type & lt, int & li) const
    {
      return side_of_edge(p, e.first, e.second, e.third, lt, li);
    }

  CGAL_Bounded_side
  side_of_triangle(const Point & p,
		   const Point & p0, 
		   const Point & p1,
		   const Point & p2,
		   Locate_type & lt, int & i, int & j ) const
    // p0,p1,p2 supposed to define a plane
    // p supposed to lie on plane p0,p1,p2
    // triangle p0,p1,p2 defines the orientation of the plane
    // returns
    // CGAL_ON_BOUNDED_SIDE if p lies strictly inside the triangle
    // CGAL_ON_BOUNDARY if p lies on one of the edges
    // CGAL_ON_UNBOUNDED_SIDE if p lies strictly outside the triangle
    {
      CGAL_triangulation_precondition
	( ! geom_traits().collinear(p0,p1,p2) );
      CGAL_triangulation_precondition
	( geom_traits().orientation(p,p0,p1,p2) == CGAL_COPLANAR );

      // edge p0 p1 :
      CGAL_Orientation o0 = geom_traits().orientation_in_plane(p,p0,p1,p2);
      // edge p1 p2 :
      CGAL_Orientation o1 = geom_traits().orientation_in_plane(p,p1,p2,p0);
      // edge p2 p0 :
      CGAL_Orientation o2 = geom_traits().orientation_in_plane(p,p2,p0,p1);

      if ( (o0 == CGAL_NEGATIVE) ||
	   (o1 == CGAL_NEGATIVE) ||
	   (o2 == CGAL_NEGATIVE) ) {
	lt = OUTSIDE_CONVEX_HULL;
	return CGAL_ON_UNBOUNDED_SIDE;
      }

      // now all the oi's are >=0
      // sum gives the number of edges p lies on
      int sum = ( o0 == CGAL_ZERO ) + ( o1 == CGAL_ZERO ) +
	( o2 == CGAL_ZERO );

      switch (sum) {
      case 0:
	{
	  lt = FACET;
	  return CGAL_ON_BOUNDED_SIDE;
	}
      case 1:
	{
	  lt = EDGE;
	  i = ( o0 == CGAL_ZERO ) ? 0 :
	    ( o1 == CGAL_ZERO ) ? 1 :
	    2;
	  if ( i == 2 ) { j=0; }
	  else { j = i+1; }
	  return CGAL_ON_BOUNDARY;
	}
      case 2:
	{
	  lt = VERTEX;
	  i = ( o0 == CGAL_POSITIVE ) ? 2 :
	    ( o1 == CGAL_POSITIVE ) ? 0 :
	    1;
	  return CGAL_ON_BOUNDARY;
	}
      default:
	{
	  // cannot happen
	  CGAL_triangulation_assertion(false);
	  return CGAL_ON_BOUNDARY;
	}
      }
    }

  
  CGAL_Bounded_side
  side_of_facet(const Point & p,
		Cell_handle c,
		int i,
		Locate_type & lt, int & li, int & lj) const
    // supposes dimension 2 otherwise does not work for infinite facets
    // returns :
    // CGAL_ON_BOUNDED_SIDE if p inside the facet
    // (for an infinite facet this means that p lies strictly in the half plane
    // limited by its finite edge)
    // CGAL_ON_BOUNDARY if p on the boundary of the facet
    // (for an infinite facet this means that p lies on the *finite* edge)
    // CGAL_ON_UNBOUNDED_SIDE if p lies outside the facet
    // (for an infinite facet this means that p is not in the preceding two cases)
    // lt has a meaning only when CGAL_ON_BOUNDED_SIDE or CGAL_ON_BOUNDARY
    // when they mean anything, li and lj refer to indices in the cell c 
    // giving the facet (c,i)
    {//side_of_facet
      CGAL_triangulation_precondition( dimension() == 2 );
      if ( ! is_infinite(c,i) ) {
	int i0, i1, i2; // indices in the considered facet
	CGAL_Bounded_side side;
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
	default:
	  {
	    // impossible
	    CGAL_triangulation_assertion( false );
	    // to avoid warning at compile time :
	    return side_of_triangle(p,
				    c->vertex(1)->point(),
				    c->vertex(2)->point(),
				    c->vertex(3)->point(),
				    lt, li, lj);
	  }
	}
	int i_t, j_t;
	side = side_of_triangle(p,
				c->vertex(i0)->point(),
				c->vertex(i1)->point(),
				c->vertex(i2)->point(),
				lt, i_t, j_t);
	// indices in the original cell :
	li = ( i_t == 0 ) ? i0 :
	  ( i_t == 1 ) ? i1 :
	  i2;
	lj = ( j_t == 0 ) ? i0 :
	  ( j_t == 1 ) ? i1 :
	  i2;
	return side;
      }
      else { // infinite facet
	int inf = c->index(infinite);
	int i1,i2; // indices in the facet
	if ( i == (inf+1)&3 ) {
	  i1 = (inf+2)&3;
	  i2 = (inf+3)&3;
	}
	else {
	  if ( i == (inf+2)&3 ) {
	    i1 = (inf+3)&3;
	    i2 = (inf+1)&3;
	  }
	  else {
	    i1 = (inf+1)&3;
	    i2 = (inf+2)&3;
	  }
	}
	Vertex_handle 
	  v1 = c->vertex(i1),
	  v2 = c->vertex(i2);
	
	// does not work in dimension 3
	Cell_handle n = c->neighbor(inf);
	// n must be a finite cell
	CGAL_Orientation o =
	  geom_traits().orientation_in_plane( p,
				     v1->point(), 
				     v2->point(), 
				     n->vertex(n->index(c))->point() );
	switch (o) {
	case CGAL_POSITIVE:
	  // p lies on the same side of v1v2 as vn, so not in f
	  {
	    return CGAL_ON_UNBOUNDED_SIDE;
	  }
	case CGAL_NEGATIVE:
	  // p lies in f
	  { 
	    lt = FACET;
	    li = i;
	    return CGAL_ON_BOUNDED_SIDE;
	  }
	case CGAL_ZERO:
	  // p collinear with v1v2
	  {
	    int i_e;
	    CGAL_Bounded_side side = 
	      side_of_segment( p,
			       v1->point(), v2->point(),
			       lt, i_e );
	    switch (side) {
	      // computation of the indices inthe original cell
	    case CGAL_ON_BOUNDED_SIDE:
	      {
		// lt == EDGE ok
		li = i1;
		lj = i2;
		return CGAL_ON_BOUNDARY;
	      }
	    case CGAL_ON_BOUNDARY:
	      {
		// lt == VERTEX ok
		li = ( i_e == 0 ) ? i1 : i2;
		return CGAL_ON_BOUNDARY;
	      }
	    case CGAL_ON_UNBOUNDED_SIDE:
	      {
		// p lies on the line defined by the finite edge
		return CGAL_ON_UNBOUNDED_SIDE;
	      }
	    default:
	      {
		// cannot happen. only to avoid warning with eg++
		return CGAL_ON_UNBOUNDED_SIDE;
	      }
	    } 
	  }// case CGAL_ZERO
	}// switch o
      }// end infinite facet
      // cannot happen. only to avoid warning with eg++
      return CGAL_ON_UNBOUNDED_SIDE;
    }

  CGAL_Bounded_side
  side_of_facet(const Point & p,
		Facet f,
		Locate_type & lt, int & li, int & lj) const
    {
      return side_of_facet(p, f.first, f.second, lt, li, lj);
    }

  CGAL_Bounded_side
  side_of_tetrahedron(const Point & p,
		      const Point & p0, 
		      const Point & p1,
		      const Point & p2, 
		      const Point & p3,
		      Locate_type & lt, int & i, int & j ) const
    // p0,p1,p2,p3 supposed to be non coplanar
    // tetrahedron p0,p1,p2,p3 is supposed to be well oriented
    // returns :
    // CGAL_ON_BOUNDED_SIDE if p lies strictly inside the tetrahedron
    // CGAL_ON_BOUNDARY if p lies on one of the facets
    // CGAL_ON_UNBOUNDED_SIDE if p lies strictly outside the tetrahedron

    // ?? locate type...
    {
      CGAL_triangulation_precondition
	( geom_traits().orientation(p0,p1,p2,p3) == CGAL_POSITIVE );
				      
      CGAL_Orientation o0 = geom_traits().orientation(p,p1,p2,p3);
      CGAL_Orientation o1 = geom_traits().orientation(p0,p,p2,p3);
      CGAL_Orientation o2 = geom_traits().orientation(p0,p1,p,p3);
      CGAL_Orientation o3 = geom_traits().orientation(p0,p1,p2,p);

      if ( (o0 == CGAL_NEGATIVE) ||
	   (o1 == CGAL_NEGATIVE) ||
	   (o2 == CGAL_NEGATIVE) ||
	   (o3 == CGAL_NEGATIVE) ) {
	lt = OUTSIDE_CONVEX_HULL;
	return CGAL_ON_UNBOUNDED_SIDE;
      }

      // now all the oi's are >=0
      // sum gives the number of facets p lies on
      int sum = ( o0 == CGAL_ZERO ) + ( o1 == CGAL_ZERO ) 
	+ ( o2 == CGAL_ZERO ) + ( o3 == CGAL_ZERO );

      switch (sum) {
      case 0:
	{
	  lt = CELL;
	  return CGAL_ON_BOUNDED_SIDE;
	}
      case 1:
	{
	  lt = FACET;
	  // i = index such that p lies on facet(i)
	  i = ( o0 == CGAL_ZERO ) ? 0 :
	    ( o1 == CGAL_ZERO ) ? 1 :
	    ( o2 == CGAL_ZERO ) ? 2 :
	    3;
	  return CGAL_ON_BOUNDARY;
	}
      case 2:
	{
	  lt = EDGE;
	  // i = smallest index such that p does not lie on facet(i)
	  // i must be < 3 since p lies on 2 facets
	  i = ( o0 == CGAL_POSITIVE ) ? 0 :
	    ( o1 == CGAL_POSITIVE ) ? 1 :
	    2;
	  // j = larger index such that p not on facet(j)
	  // j must be > 0 since p lies on 2 facets
	  j = ( o3 == CGAL_POSITIVE ) ? 3 :
	    ( o2 == CGAL_POSITIVE ) ? 2 :
	    1;
	  return CGAL_ON_BOUNDARY;
	}
      case 3:
	{
	  lt = VERTEX;
	  // i = index such that p does not lie on facet(i)
	  i = ( o0 == CGAL_POSITIVE ) ? 0 :
	    ( o1 == CGAL_POSITIVE ) ? 1 :
	    ( o2 == CGAL_POSITIVE ) ? 2 :
	    3;
	  return CGAL_ON_BOUNDARY;
	}
      default:
	{
	  // impossible : cannot be on 4 facets for a real tetrahedron
	  CGAL_triangulation_assertion(false);
	  return CGAL_ON_BOUNDARY;
	}
      }
    }

  CGAL_Bounded_side
  side_of_cell(const Point & p, 
	       Cell_handle c,
	       Locate_type & lt, int & i, int & j) const
    // returns
    // CGAL_ON_BOUNDED_SIDE if p inside the cell
    // (for an infinite cell this means that p lies strictly in the half space
    // limited by its finite facet)
    // CGAL_ON_BOUNDARY if p on the boundary of the cell
    // (for an infinite cell this means that p lies on the *finite* facet)
    // CGAL_ON_UNBOUNDED_SIDE if p lies outside the cell
    // (for an infinite cell this means that p is not in the preceding two cases)
    // lt has a meaning only when CGAL_ON_BOUNDED_SIDE or CGAL_ON_BOUNDARY
    {
      CGAL_triangulation_precondition( dimension() == 3 );
      if ( ! is_infinite(c) ) {
	return side_of_tetrahedron(p,
				   c->vertex(0)->point(),
				   c->vertex(1)->point(),
				   c->vertex(2)->point(),
				   c->vertex(3)->point(),
				   lt, i, j);
      }
      else {
	int inf = c->index(infinite);
	CGAL_Orientation o;
	Vertex_handle 
	  v1=c->vertex((inf+1)&3), 
	  v2=c->vertex((inf+2)&3), 
	  v3=c->vertex((inf+3)&3);
	if ( inf%2 == 0 ) {
	  o = geom_traits().orientation(p,
			       v1->point(),
			       v2->point(),
			       v3->point());
	}
	else {
	  o =  geom_traits().orientation(v3->point(),
				p,
				v1->point(),
				v2->point());
	}
	switch (o) {
	case CGAL_POSITIVE:
	  {
	    lt = CELL;
	    return CGAL_ON_BOUNDED_SIDE;
	  }
	case CGAL_NEGATIVE:
	  return CGAL_ON_UNBOUNDED_SIDE;
	case CGAL_ZERO:
	  {
	    // location in the finite facet
	    int i_f, j_f;
	    CGAL_Bounded_side side = 
	      side_of_triangle(p,
			       v1->point(),
			       v2->point(),
			       v3->point(),
			       lt, i_f, j_f);
	    // lt need not be modified in most cases :
	    switch (side) {
	    case CGAL_ON_BOUNDED_SIDE:
	      {
		// lt == FACET ok
		i = inf;
		return CGAL_ON_BOUNDARY;
	      }
	    case CGAL_ON_BOUNDARY:
	      {
		// lt == VERTEX OR EDGE ok
		i = ( i_f == 0 ) ? ((inf+1)&3) :
		  ( i_f == 1 ) ? ((inf+2)&3) :
		  ((inf+3)&3);
		if ( lt == EDGE ) {
		  j = (j_f == 0 ) ? ((inf+1)&3) :
		    ( j_f == 1 ) ? ((inf+2)&3) :
		    ((inf+3)&3);
		}
		return CGAL_ON_BOUNDARY;
	      }
	    case CGAL_ON_UNBOUNDED_SIDE:
	      {
		// p lies on the plane defined by the finite facet
		// lt must be initialized
		return CGAL_ON_UNBOUNDED_SIDE;
	      }
	    default:
	      {
		CGAL_triangulation_assertion(false);
		return CGAL_ON_BOUNDARY;
	      }
	    } // switch side
	  }// case CGAL_ZERO
	default:
	  {
	    CGAL_triangulation_assertion(false);
	    return CGAL_ON_BOUNDARY;
	  }
	} // switch o
      } // else infinite cell
    } // side_of_cell
};



// template <class GT, class Tds >
// ostream& operator<<
// (ostream& os, const CGAL_Triangulation_3<GT, Tds> &tr)
// {
//     return os ;
// }

// template < class GT, class Tds >
// istream& operator>>
// (istream& is, CGAL_Triangulation_3<GT, Tds> &tr)
// {
//   return operator>>(is, tr._Triangulation_data_structure_3);
// }
    

#endif CGAL_TRIANGULATION_3_H

