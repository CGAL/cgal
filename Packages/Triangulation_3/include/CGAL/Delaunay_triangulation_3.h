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
// file          : include/CGAL/Delaunay_triangulation_3.h
// revision      : $Revision$
//
// author(s)     : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//                 Sylvain Pion <Sylvain.Pion@sophia.inria.fr>
//                 Andreas Fabri <Andreas.Fabri@sophia.inria.fr>
//
// coordinator   : INRIA Sophia Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================


#ifndef CGAL_DELAUNAY_TRIANGULATION_3_H
#define CGAL_DELAUNAY_TRIANGULATION_3_H

#include <CGAL/basic.h>

#include <set>
#include <utility>

#include <CGAL/Triangulation_short_names_3.h>
#include <CGAL/Triangulation_utils_3.h>
#include <CGAL/Triangulation_3.h>

#include <CGAL/Delaunay_remove_tds_3.h>
#include <CGAL/Triangulation_face_base_2.h>

CGAL_BEGIN_NAMESPACE

template < class Gt, 
           class Tds = Triangulation_data_structure_3 <
                                   Triangulation_vertex_base_3<Gt>,
                                   Triangulation_cell_base_3<void> > >
class Delaunay_triangulation_3 : public Triangulation_3<Gt,Tds>
{
  friend std::istream& operator >> CGAL_NULL_TMPL_ARGS
  (std::istream& is, Triangulation_3<Gt,Tds> &tr);

public:
  typedef Tds Triangulation_data_structure;
  typedef Gt  Geom_traits;
  typedef Delaunay_triangulation_3<Gt, Tds> Self;

  typedef typename Gt::Point_3       Point;
  typedef typename Gt::Segment_3     Segment;
  typedef typename Gt::Triangle_3    Triangle;
  typedef typename Gt::Tetrahedron_3 Tetrahedron;

  // types for dual:
  typedef typename Gt::Line_3        Line;
  typedef typename Gt::Ray_3         Ray;
  typedef typename Gt::Plane_3       Plane;
  typedef typename Gt::Direction_3   Direction;
  typedef typename Gt::Object_3      Object;

  typedef typename Triangulation_3<Gt,Tds>::Cell_handle   Cell_handle;
  typedef typename Triangulation_3<Gt,Tds>::Vertex_handle Vertex_handle;

  typedef typename Triangulation_3<Gt,Tds>::Cell   Cell;
  typedef typename Triangulation_3<Gt,Tds>::Vertex Vertex;
  typedef typename Triangulation_3<Gt,Tds>::Facet  Facet;
  typedef typename Triangulation_3<Gt,Tds>::Edge   Edge;

  typedef typename Triangulation_3<Gt,Tds>::Cell_circulator Cell_circulator;
  typedef typename Triangulation_3<Gt,Tds>::Cell_iterator   Cell_iterator;
  typedef typename Triangulation_3<Gt,Tds>::Facet_iterator  Facet_iterator;
  typedef typename Triangulation_3<Gt,Tds>::Edge_iterator   Edge_iterator;
  typedef typename Triangulation_3<Gt,Tds>::Vertex_iterator Vertex_iterator;

  typedef typename Triangulation_3<Gt,Tds>::Locate_type Locate_type;

protected:

  Oriented_side
  side_of_oriented_sphere(const Point &p, const Point &q, const Point &r,
	                  const Point &s, const Point &t) const
  {
      return geom_traits().side_of_oriented_sphere_3_object()(p, q, r, s, t);
  }

  Bounded_side
  coplanar_side_of_bounded_circle(const Point &p, const Point &q,
	                          const Point &r, const Point &s) const
  {
      return geom_traits().coplanar_side_of_bounded_circle_3_object()(p,q,r,s);
  }

  // for dual:
  Point
  construct_circumcenter(const Point &p, const Point &q, const Point &r) const
  {
      return geom_traits().construct_circumcenter_3_object()(p, q, r);
  }

  Point
  construct_circumcenter(const Point &p, const Point &q,
	                 const Point &r, const Point &s) const
  {
      return geom_traits().construct_circumcenter_3_object()(p, q, r, s);
  }

  Line
  construct_perpendicular_line(const Plane &pl, const Point &p) const
  {
      return geom_traits().construct_perpendicular_line_3_object()(pl, p);
  }

  Plane
  construct_plane(const Point &p, const Point &q, const Point &r) const
  {
      return geom_traits().construct_plane_3_object()(p, q, r);
  }

  Direction
  construct_direction_of_line(const Line &l) const
  {
      return geom_traits().construct_direction_of_line_3_object()(l);
  }

  Ray
  construct_ray(const Point &p, const Direction &d) const
  {
      return geom_traits().construct_ray_3_object()(p, d);
  }

  Object
  construct_object(const Point &p) const
  {
      return geom_traits().construct_object_3_object()(p);
  }

  Object
  construct_object(const Segment &s) const
  {
      return geom_traits().construct_object_3_object()(s);
  }

  Object
  construct_object(const Ray &r) const
  {
      return geom_traits().construct_object_3_object()(r);
  }

public:

  Delaunay_triangulation_3()
    : Triangulation_3<Gt,Tds>()
  {}
  
  Delaunay_triangulation_3(const Gt & gt)
    : Triangulation_3<Gt,Tds>(gt)
  {}
  
  // copy constructor duplicates vertices and cells
  Delaunay_triangulation_3(const Delaunay_triangulation_3<Gt,Tds> & tr)
    : Triangulation_3<Gt,Tds>(tr)
  { 
    CGAL_triangulation_postcondition( is_valid() );  
  }

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

  Vertex_handle insert(const Point & p, 
		       Cell_handle start = Cell_handle(), 
		       Vertex_handle v = NULL);

  Vertex_handle push_back(const Point & p)
  {
      return insert(p);
  }

  bool remove(Vertex_handle v );

  void print(Vertex_handle v) const;
  void print(Cell_handle c) const;

  Bounded_side
  side_of_sphere( Cell_handle c, const Point & p) const;

  Bounded_side
  side_of_circle( const Facet & f, const Point & p) const
    {
      return side_of_circle(f.first, f.second, p);
    }

  Bounded_side
  side_of_circle( Cell_handle c, int i, const Point & p) const;

  Point dual(Cell_handle c) const;

  Object dual(const Facet & f) const
    { return dual( f.first, f.second ); }
  Object dual(Cell_handle c, int i) const;

  bool is_valid(bool verbose = false, int level = 0) const;

  bool is_valid(Cell_handle c, bool verbose = false, int level = 0) const;

  template < class Stream> 		
  Stream& draw_dual(Stream & os)
    {
      Facet_iterator fit = finite_facets_begin();
      for (; fit != facets_end(); ++fit) {
	Object o = dual(*fit);
	Point p;
	Ray r;
	Segment s;
	if (CGAL::assign(p,o)) os << p;
	if (CGAL::assign(s,o)) os << s;
	if (CGAL::assign(r,o)) os << r; 
      }
      return os;
    }

private:

  Bounded_side
  side_of_sphere_inf_perturb(Vertex_handle v0, 
			     Vertex_handle v1, 
			     Vertex_handle v2, 
			     Vertex_handle v) const;

  Bounded_side
  side_of_sphere_finite_perturb(Vertex_handle v0, 
				Vertex_handle v1, 
				Vertex_handle v2, 
				Vertex_handle v3, 
				Vertex_handle v) const;

  int max2(int i0, int i1, int i2, int i3, int i4, int m) const;
  int maxless(int i0, int i1, int i2, int i3, int i4, int m) const;

  void delete_cells(std::vector<Cell_handle> & hole);

  void make_hole_3D_ear( Vertex_handle v, 
	                 std::vector<Facet> & boundhole,
	                 std::vector<Cell_handle> & hole);
  void undo_make_hole_3D_ear(std::vector<Facet> & boundhole,
		             std::vector<Cell_handle> & hole);
  bool fill_hole_3D_ear(std::vector<Facet> & boundhole);

private:

  class Conflict_tester_3
  {
      const Point &p;
      Self *t;

  public:

      Conflict_tester_3(const Point &pt, Self *tr)
	  : p(pt), t(tr) {}

      bool operator()(const typename Tds::Cell *c) const
      {
	  return t->side_of_sphere((Cell_handle)(Cell*)c, p)
	      == ON_BOUNDED_SIDE;
      }
  };

  class Conflict_tester_2
  {
      const Point &p;
      Self *t;

  public:

      Conflict_tester_2(const Point &pt, Self *tr)
	  : p(pt), t(tr) {}

      bool operator()(const typename Tds::Cell *c) const
      {
	  return t->side_of_circle((Cell_handle)(Cell*)c, 3, p)
	      == ON_BOUNDED_SIDE;
      }
  };
};

template < class Gt, class Tds >
Delaunay_triangulation_3<Gt,Tds>::Vertex_handle
Delaunay_triangulation_3<Gt,Tds>::
insert(const Point & p, Cell_handle start, Vertex_handle v)
{
  switch (dimension()) {
  case 3:
    {
      Locate_type lt;
      int li, lj;
      Cell_handle c = locate( p, lt, li, lj, start);
      if ( lt == VERTEX )
	  return c->vertex(li);

      set_number_of_vertices(number_of_vertices()+1);
      Conflict_tester_3 tester(p, this);
      v = (Vertex *) _tds.insert_conflict(&(*v), &(*c), tester); 
      v->set_point(p);
      return v;
    }// dim 3
  case 2:
    {
      Locate_type lt;
      int li, lj;
      Cell_handle c = locate( p, lt, li, lj, start);
      switch (lt) {
      case OUTSIDE_CONVEX_HULL:
      case CELL:
      case FACET:
      case EDGE:
	{
	  set_number_of_vertices(number_of_vertices()+1);
          Conflict_tester_2 tester(p, this);
	  v = (Vertex *) _tds.insert_conflict(&(*v), &(*c), tester); 
	  v->set_point(p);
	  return v;
	}
      case VERTEX:
	return c->vertex(li);
      case OUTSIDE_AFFINE_HULL:
	  // if the 2d triangulation is Delaunay, the 3d
	  // triangulation will be Delaunay
	return Triangulation_3<Gt,Tds>::insert_outside_affine_hull(p,v); 
      }
    }//dim 2
  default :
    // dimension <= 1
    return Triangulation_3<Gt,Tds>::insert(p,start,v);
  }
}// insert(p)

template < class Gt, class Tds >
bool
Delaunay_triangulation_3<Gt,Tds>::
remove(Vertex_handle v)
{
  CGAL_triangulation_precondition( v != Vertex_handle());
  CGAL_triangulation_precondition( !is_infinite(v));
  CGAL_triangulation_expensive_precondition( _tds.is_vertex(&(*v)) );

  if ( dimension() <3 || test_dim_down(v) ) {
    // the triangulation is rebuilt...

    Vertex_handle inf = infinite_vertex();

    std::vector< CGAL_TYPENAME_MSVC_NULL Tds::Vertex * > Vertices;
    _tds.clear_cells_only(Vertices); 

    _tds.set_number_of_vertices(0);
    _tds.set_dimension(-2);
    _tds.insert_increase_dimension( &(*inf) );
    CGAL_triangulation_assertion( inf == infinite_vertex() );

    typename std::vector<CGAL_TYPENAME_MSVC_NULL Tds::Vertex * >::iterator vit;

    for ( vit = Vertices.begin(); vit != Vertices.end(); ++vit ) 
      if ( *vit != &(*inf) && *vit != &(*v) ) 
	insert( (*vit)->point(), NULL, (Vertex *) *vit );

    _tds.delete_vertex(&(*v));

    return true;
  }

  CGAL_triangulation_assertion( dimension() == 3 );

  std::vector<Facet> boundhole; // facets on the boundary of the hole
  boundhole.reserve(64);        // 27 on average.
  std::vector<Cell_handle> hole;
  hole.reserve(64);

  make_hole_3D_ear(v, boundhole, hole);

  bool filled = fill_hole_3D_ear(boundhole);
  if(filled){
    _tds.delete_vertex(&(*v));
    delete_cells(hole);
    set_number_of_vertices(number_of_vertices()-1);
  } else {
    undo_make_hole_3D_ear(boundhole, hole);
  }

  return filled;
}

template < class Gt, class Tds >
void
Delaunay_triangulation_3<Gt,Tds>::
delete_cells(std::vector<Cell_handle> & hole)
{
  for(typename std::vector<Cell_handle>::iterator cit = hole.begin();
      cit != hole.end(); ++cit)
    _tds.delete_cell( &*(*cit) );
}

//debug
template < class Gt, class Tds >
void 
Delaunay_triangulation_3<Gt,Tds>::
print( Vertex_handle v ) const
{
  if ( is_infinite( v ) )
    std::cout << "inf" << "; ";
  else 
    std::cout << v->point() << "; ";
}

template < class Gt, class Tds >
void 
Delaunay_triangulation_3<Gt,Tds>::
print( Cell_handle c ) const
{
  for(int i = 0; i < 4; i++)
    print(c->vertex(i));
  std::cout << std::endl;
}
// end debug

template < class Gt, class Tds >
Bounded_side
Delaunay_triangulation_3<Gt,Tds>::
side_of_sphere(Cell_handle c, const Point & p) const
{
  CGAL_triangulation_precondition( dimension() == 3 );
  int i3;
  if ( ! c->has_vertex( infinite_vertex(), i3 ) ) 
    return Bounded_side( side_of_oriented_sphere( c->vertex(0)->point(),
						  c->vertex(1)->point(),
						  c->vertex(2)->point(),
						  c->vertex(3)->point(),
						  p ) );
  // else infinite cell :
  unsigned char i[3] = {(i3+1)&3, (i3+2)&3, (i3+3)&3};
  if ( (i3&1) == 0 )
      std::swap(i[0], i[1]);
  Orientation o = orientation( c->vertex(i[0])->point(),
			       c->vertex(i[1])->point(),
			       c->vertex(i[2])->point(),
			       p );
  if (o != ZERO)
    return Bounded_side(o);

  return coplanar_side_of_bounded_circle( c->vertex(i[0])->point(), 
					  c->vertex(i[1])->point(),
					  c->vertex(i[2])->point(),
					  p );
}

template < class Gt, class Tds >
Bounded_side
Delaunay_triangulation_3<Gt,Tds>::
side_of_sphere_inf_perturb(Vertex_handle v0, 
			   Vertex_handle v1, 
			   Vertex_handle v2, 
			   Vertex_handle v) const
  // v0,v1,v2 are supposed to form a facet of the triangulation
{
  CGAL_triangulation_precondition( (dimension() == 3) &&
				   (! is_infinite(v0)) &&
				   (! is_infinite(v1)) &&
				   (! is_infinite(v2)) &&
				   (! is_infinite(v)) );
  
  const Point & p0 = v0->point();
  const Point & p1 = v1->point();
  const Point & p2 = v2->point();
  const Point & p = v->point();

  Orientation o = orientation(p0,p1,p2,p);

  if (o != ZERO)
    return Bounded_side(o);

  Bounded_side bs = coplanar_side_of_bounded_circle (p0,p1,p2,p);

  if ( bs == ON_BOUNDARY ) {
    // the 4 points are coplanar, cocircular
    // the facet is supposed to exist, so:
    return ON_UNBOUNDED_SIDE;
  }
  
  return bs;
}

template < class Gt, class Tds >
int
Delaunay_triangulation_3<Gt,Tds>::
max2(int i0, int i1, int i2, int i3, int i4, int m) const
  // m supposed to be the maximum of the i's
  // could be replaced by maxless
{
  if (m == i0) return std::max(i1,std::max(i2,std::max(i3,i4)));
  if (m == i1) return std::max(i0,std::max(i2,std::max(i3,i4)));
  if (m == i2) return std::max(i0,std::max(i1,std::max(i3,i4)));
  if (m == i3) return std::max(i0,std::max(i1,std::max(i2,i4)));
  return std::max(i0,std::max(i1,std::max(i2,i3)));
}

template < class Gt, class Tds >
int
Delaunay_triangulation_3<Gt,Tds>::
maxless(int i0, int i1, int i2, int i3, int i4, int m) const
  // returns the largest i which is less than m
{
  int mtmp=-1;
  if (i0 < m) // we know that i0 > -1
    mtmp = i0;
  if ( (i1 < m) && (i1 > mtmp) )
    mtmp = i1;
  if ( (i2 < m) && (i2 > mtmp) )
    mtmp = i2;
  if ( (i3 < m) && (i3 > mtmp) )
    mtmp = i3;
  if ( (i4 < m) && (i4 > mtmp) )
    mtmp = i4;

  return mtmp;
}

template < class Gt, class Tds >
Bounded_side
Delaunay_triangulation_3<Gt,Tds>::
side_of_sphere_finite_perturb(Vertex_handle v0, 
			      Vertex_handle v1, 
			      Vertex_handle v2, 
			      Vertex_handle v3, 
			      Vertex_handle v) const
  // v0,v1,v2 and v0,v1,v3 are supposed to form two finite facets of
  // the triangulation, forming a convex angle 
  // must not be used otherwise
{
  CGAL_triangulation_precondition( (dimension() == 3) &&
				   (! is_infinite(v0)) &&
				   (! is_infinite(v1)) &&
				   (! is_infinite(v2)) &&
				   (! is_infinite(v3)) );

  const Point & p0 = v0->point();
  const Point & p1 = v1->point();
  const Point & p2 = v2->point();
  const Point & p3 = v3->point();

  CGAL_triangulation_precondition( orientation(p0,p1,p2,p3) == POSITIVE );
				   
  if (is_infinite(v))
      return ON_UNBOUNDED_SIDE;

  const Point & p = v->point();

  Bounded_side bs = Bounded_side(side_of_oriented_sphere(p0,p1,p2,p3,p));

  if ( bs != ON_BOUNDARY )
      return bs;

  int i0 = v0->get_order_of_creation();
  int i1 = v1->get_order_of_creation();
  int i2 = v2->get_order_of_creation();
  int i3 = v3->get_order_of_creation();
  int i = v->get_order_of_creation();
  Orientation o;

  int m = std::max(i0,std::max(i1,std::max(i2,std::max(i3,i))));
  // we look whether the leading monomial of the determinant has non null 
  // coefficient 
  if (m == i) 
    return ON_UNBOUNDED_SIDE;
  // since p0 p1 p2 p3 are supposed to be non coplanar and positively oriented
  if (m == i3 && (o = orientation(p0,p1,p2,p)) != ZERO )
    return Bounded_side(o);
  if (m == i2 && (o = orientation(p0,p1,p3,p)) != ZERO )
    return Bounded_side(-o);
  if (m == i1 && (o = orientation(p0,p2,p3,p)) != ZERO )
    return Bounded_side(o);
  if (m == i0 && (o = orientation(p1,p2,p3,p)) != ZERO )
    return Bounded_side(-o);

  // if not yet returned, then the leading monomial of the determinant
  // has null coefficient 
  // we look whether the 2nd monomial has non null coefficient 
  m = max2(i0,i1,i2,i3,i,m);
  if (m == i) 
    return ON_UNBOUNDED_SIDE;
  if (m == i3 && (o = orientation(p0,p1,p2,p)) != ZERO )
    return Bounded_side(o);
  if (m == i2 && (o = orientation(p0,p1,p3,p)) != ZERO )
    return Bounded_side(-o);
  if (m == i1 && (o = orientation(p0,p2,p3,p)) != ZERO )
    return Bounded_side(o);
  if (m == i0 && (o = orientation(p1,p2,p3,p)) != ZERO )
    return Bounded_side(-o);

  // we look whether the 3rd monomial also has non null coefficient 
  m = maxless(i0,i1,i2,i3,i,m);
  if (m == i) 
    return ON_UNBOUNDED_SIDE;
  if (m == i3 && (o = orientation(p0,p1,p2,p)) != ZERO )
    return Bounded_side(o);
  if (m == i2 && (o = orientation(p0,p1,p3,p)) != ZERO )
    return Bounded_side(-o);
  if (m == i1 && (o = orientation(p0,p2,p3,p)) != ZERO )
    return Bounded_side(o);
  if (m == i0 && (o = orientation(p1,p2,p3,p)) != ZERO )
    return Bounded_side(-o);

  // we look whether the 4th monomial also has non null coefficient 
  m = maxless(i0,i1,i2,i3,i,m);
  if (m == i) 
    return ON_UNBOUNDED_SIDE;
  if (m == i3 && (o = orientation(p0,p1,p2,p)) != ZERO )
    return Bounded_side(o);
  if (m == i2 && (o = orientation(p0,p1,p3,p)) != ZERO )
    return Bounded_side(-o);
  if (m == i1 && (o = orientation(p0,p2,p3,p)) != ZERO )
    return Bounded_side(o);
  if (m == i0 && (o = orientation(p1,p2,p3,p)) != ZERO )
    return Bounded_side(-o);

  // case when the first non null coefficient is the coefficient of 
  // the 5th monomial
  // moreover, the tests (m == i) were false up to here, so the
  // monomial corresponding to i is the only monomial with non-zero
  // coefficient, it is equal to orient(p0,p1,p2,p3) == positive 
  // so, no further test is required
  return ON_UNBOUNDED_SIDE; 
}

template < class Gt, class Tds >
Bounded_side
Delaunay_triangulation_3<Gt,Tds>::
side_of_circle(Cell_handle c, int i, const Point & p) const
  // precondition : dimension >=2
  // in dimension 3, - for a finite facet
  // returns ON_BOUNDARY if the point lies on the circle,
  // ON_UNBOUNDED_SIDE when exterior, ON_BOUNDED_SIDE
  // interior
  // for an infinite facet, considers the plane defined by the
  // adjacent finite facet of the same cell, and does the same as in 
  // dimension 2 in this plane
  // in dimension 2, for an infinite facet
  // in this case, returns ON_BOUNDARY if the point lies on the 
  // finite edge (endpoints included) 
  // ON_BOUNDED_SIDE for a point in the open half-plane
  // ON_UNBOUNDED_SIDE elsewhere
{
  CGAL_triangulation_precondition( dimension() >= 2 );
  int i3 = 5;

  if ( dimension() == 2 ) {
    CGAL_triangulation_precondition( i == 3 );
    // the triangulation is supposed to be valid, ie the facet
    // with vertices 0 1 2 in this order is positively oriented
    if ( ! c->has_vertex( infinite_vertex(), i3 ) ) 
      return coplanar_side_of_bounded_circle( c->vertex(0)->point(),
					      c->vertex(1)->point(),
					      c->vertex(2)->point(),
					      p );
    // else infinite facet
    // v1, v2 finite vertices of the facet such that v1,v2,infinite
    // is positively oriented
    Vertex_handle v1 = c->vertex( ccw(i3) ),
                  v2 = c->vertex( cw(i3) );
    CGAL_triangulation_assertion(coplanar_orientation(v1->point(), v2->point(),
		                 (c->mirror_vertex(i3))->point()) == NEGATIVE);
    Orientation o = coplanar_orientation(v1->point(), v2->point(), p);
    if ( o != ZERO )
	return Bounded_side( o );
    // because p is in f iff
    // is does not lie on the same side of v1v2 as vn
    int i_e;
    Locate_type lt;
    // case when p collinear with v1v2
    return side_of_segment( p,
			    v1->point(), v2->point(),
			    lt, i_e );
  }

  // else dimension == 3
  CGAL_triangulation_precondition( i >= 0 && i < 4 );
  if ( ( ! c->has_vertex(infinite_vertex(),i3) ) || ( i3 != i ) ) {
    // finite facet
    // initialization of i0 i1 i2, vertices of the facet positively 
    // oriented (if the triangulation is valid)
    int i0 = (i>0) ? 0 : 1;
    int i1 = (i>1) ? 1 : 2;
    int i2 = (i>2) ? 2 : 3;
    CGAL_triangulation_precondition( orientation( c->vertex(i0)->point(),
				                  c->vertex(i1)->point(),
				                  c->vertex(i2)->point(),
				                  p ) == COPLANAR );
    return coplanar_side_of_bounded_circle( c->vertex(i0)->point(),
					    c->vertex(i1)->point(),
					    c->vertex(i2)->point(),
					    p );
  }

  //else infinite facet
  // v1, v2 finite vertices of the facet such that v1,v2,infinite
  // is positively oriented
  Vertex_handle v1 = c->vertex( next_around_edge(i3,i) ),
                v2 = c->vertex( next_around_edge(i,i3) );
  Orientation o = (Orientation)
                  (coplanar_orientation( v1->point(), v2->point(),
			                 c->vertex(i)->point()) *
                  coplanar_orientation( v1->point(), v2->point(), p ));
  // then the code is duplicated from 2d case
  if ( o != ZERO )
      return Bounded_side( -o );
  // because p is in f iff 
  // it is not on the same side of v1v2 as c->vertex(i)
  int i_e;
  Locate_type lt;
  // case when p collinear with v1v2
  return side_of_segment( p,
			  v1->point(), v2->point(),
			  lt, i_e );
}

template < class Gt, class Tds >
Delaunay_triangulation_3<Gt,Tds>::Point
Delaunay_triangulation_3<Gt,Tds>::
dual(Cell_handle c) const
{
  CGAL_triangulation_precondition(dimension()==3);
  CGAL_triangulation_precondition( ! is_infinite(c) );
  return construct_circumcenter( c->vertex(0)->point(),
				 c->vertex(1)->point(),
				 c->vertex(2)->point(),
				 c->vertex(3)->point() );
}


template < class Gt, class Tds >
Delaunay_triangulation_3<Gt,Tds>::Object
Delaunay_triangulation_3<Gt,Tds>::
dual(Cell_handle c, int i) const
{
  CGAL_triangulation_precondition(dimension()>=2);
  CGAL_triangulation_precondition( ! is_infinite(c,i) );

  if ( dimension() == 2 ) {
    CGAL_triangulation_precondition( i == 3 );
    return construct_object( construct_circumcenter(c->vertex(0)->point(),
		                                    c->vertex(1)->point(),
					            c->vertex(2)->point()) );
  }

  // dimension() == 3
  Cell_handle n = c->neighbor(i);
  if ( ! is_infinite(c) && ! is_infinite(n) )
    return construct_object(construct_segment( dual(c), dual(n) ));

  // either n or c is infinite
  int in;
  if ( is_infinite(c) ) 
    in = n->index(c);
  else {
    n = c;
    in = i;
  }
  // n now denotes a finite cell, either c or c->neighbor(i)
  unsigned char ind[3] = {(in+1)&3,(in+2)&3,(in+3)&3};
  if ( (in&1) == 1 )
      std::swap(ind[0], ind[1]);
  const Point& p = n->vertex(ind[0])->point();
  const Point& q = n->vertex(ind[1])->point();
  const Point& r = n->vertex(ind[2])->point();
  
  Line l = construct_perpendicular_line( construct_plane(p,q,r),
					 construct_circumcenter(p,q,r) );
  return construct_object(construct_ray( dual(n),
	                                 construct_direction_of_line(l)));
}

template < class Gt, class Tds >
bool
Delaunay_triangulation_3<Gt,Tds>::
is_valid(bool verbose, int level) const
{
  if ( ! tds().is_valid(verbose,level) ) {
    if (verbose)
	std::cerr << "invalid data structure" << std::endl;
    CGAL_triangulation_assertion(false);
    return false;
  }
    
  if ( &(*infinite_vertex()) == NULL ) {
    if (verbose)
	std::cerr << "no infinite vertex" << std::endl;
    CGAL_triangulation_assertion(false);
    return false;
  }

  switch ( dimension() ) {
  case 3:
    {
      Cell_iterator it;
      for ( it = finite_cells_begin(); it != cells_end(); ++it ) {
	is_valid_finite((*it).handle());
	for (int i=0; i<4; i++ ) {
	  if ( side_of_sphere ( (*it).handle(), 
		 it->vertex( (it->neighbor(i))->index((*it).handle() ) )
		 ->point() ) == ON_BOUNDED_SIDE ) {
	    if (verbose)
	      std::cerr << "non-empty sphere " << std::endl;
	    CGAL_triangulation_assertion(false);
	    return false;
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
	for (int i=0; i<2; i++ ) {
	  if ( side_of_circle ( (*it).first, 3,
		 (*it).first->vertex( (((*it).first)->neighbor(i))
			    ->index((*it).first) )->point() )
	       == ON_BOUNDED_SIDE ) {
	    if (verbose)
	      std::cerr << "non-empty circle " << std::endl;
	    CGAL_triangulation_assertion(false);
	    return false;
	  }
	}
      }
      break;
    }
  case 1:
    {
      Edge_iterator it;
      for ( it = finite_edges_begin(); it != edges_end(); ++it )
	is_valid_finite((*it).first);
      break;
    }
  }
  if (verbose)
      std::cerr << "Delaunay valid triangulation" << std::endl;
  return true;
}

template < class Gt, class Tds >
bool
Delaunay_triangulation_3<Gt,Tds>::
is_valid(Cell_handle c, bool verbose, int level) const
{
  if ( ! (&(*c))->is_valid(dimension(),verbose,level) ) {
    if (verbose) { 
      std::cerr << "combinatorically invalid cell" ;
      for (int i=0; i <= dimension(); i++ )
	std::cerr << c->vertex(i)->point() << ", " ;
      std::cerr << std::endl;
    }
    CGAL_triangulation_assertion(false);
    return false;
  }
  switch ( dimension() ) {
  case 3:
    {
      if ( ! is_infinite(c) ) {
	is_valid_finite(c,verbose,level);
	for (int i=0; i<4; i++ ) {
	  if (side_of_sphere(c, c->vertex((c->neighbor(i))->index(c))->point())
	      == ON_BOUNDED_SIDE ) {
	    if (verbose)
		std::cerr << "non-empty sphere " << std::endl;
	    CGAL_triangulation_assertion(false);
	    return false;
	  }
	}
      }
      break;
    }
  case 2:
    {
      if ( ! is_infinite(c,3) ) {
	for (int i=0; i<2; i++ ) {
	  if (side_of_circle(c, 3, c->vertex(c->neighbor(i)->index(c))->point())
	       == ON_BOUNDED_SIDE ) {
	    if (verbose)
		std::cerr << "non-empty circle " << std::endl;
	    CGAL_triangulation_assertion(false);
	    return false;
	  }
	}
      }
      break;
    }
  }
  if (verbose)
      std::cerr << "Delaunay valid cell" << std::endl;
  return true;
}

template < class Gt, class Tds >
void
Delaunay_triangulation_3<Gt,Tds>::
make_hole_3D_ear( Vertex_handle v, 
	      std::vector<Facet> & boundhole,
	      std::vector<Cell_handle> & hole)
{
  CGAL_triangulation_expensive_precondition( ! test_dim_down(v) );

  typedef std::set<Cell_handle> Hole_cells;
  Hole_cells cells;
  incident_cells( v, cells );

  Cell_handle opp_cit;
  Vertex_handle vi;

  for (typename Hole_cells::iterator cit = cells.begin();
       cit != cells.end(); ++cit) {
    int indv = (*cit)->index(&(*v));
    opp_cit = (*cit)->neighbor( indv );
    hole.push_back(*cit);    
    boundhole.push_back( std::make_pair( opp_cit, opp_cit->index(*cit)) );

    for (int i=0; i<4; i++)
      if ( i != indv ) {
	vi = (*cit)->vertex(i);
	vi->set_cell( opp_cit );
      }
  }
}

template < class Gt, class Tds >
void
Delaunay_triangulation_3<Gt,Tds>::
undo_make_hole_3D_ear(std::vector<Facet> & boundhole,
		      std::vector<Cell_handle> & hole)
{
  typename std::vector<Cell_handle>::iterator cit = hole.begin();
  for(typename std::vector<Facet>::iterator fit = boundhole.begin();	
      fit != boundhole.end(); ++fit) {
    Cell_handle ch = (*fit).first;
    ch->set_neighbor((*fit).second, *cit);
    ++cit;
    // the vertices on the boundary of the hole still
    // point to the cells that form the boundary of the hole
  }
}

template < class Gt, class Tds >
bool
Delaunay_triangulation_3<Gt,Tds>::
fill_hole_3D_ear( std::vector<Facet> & boundhole)
{
  typedef Delaunay_remove_tds_3_2<Delaunay_triangulation_3> Surface;
  typedef typename Surface::Face_3_2 Face_3_2;
  typedef typename Surface::Vertex_3_2 Vertex_3_2;
  typedef typename Surface::Vertex_circulator Vertex_circulator_3_2;

  // The list of cells that gets created, so that we know what
  // we have to delete, in case that we cannot fill the hole
  std::vector<Cell_handle> cells;
  cells.reserve(32); // 20 on average.

  Surface surface(boundhole);

  Face_3_2 *f = &(* surface.faces_begin());
  //  Face_3_2 *last_op = NULL; // This is where the last ear was inserted
  Face_3_2 *last_op = &(* surface.faces_begin()); // This is where the last
                                                  // ear was inserted
  int k = -1;

  // This is a loop over the halfedges of the surface of the hole
  // As edges are not explicitely there, we loop over the faces instead,
  // and an index. 
  // The current face is f, the current index is k = -1, 0, 1, 2
  for(;;) {
    next_edge: ;
    k++;
    if(k == 3) {
      // The faces form a circular list. With f->n() we go to the next face.
      f = (Face_3_2*)f->n();
      if(f == last_op) {
	// We looked at all edges without doing anything, that is we are
	// in an infinite loop. ==> Panic mode, delete created cells.
	std::cerr << "\nUnable to find an ear\n" << std::endl;
	//	std::cerr <<  surface  << std::endl;
	delete_cells(cells);
	return false;
      }
      k = 0;
    }

    // The edges are marked, if they are a candidate for an ear.
    // This saves time, for example an edge gets not considered
    // from both adjacent faces.
    if(f->is_halfedge_marked(k)) {
      Vertex_3_2 *w0, *w1, *w2, *w3;
      Vertex_handle v0, v1, v2, v3;
      int i = ccw(k);
      int j = cw(k);
      Face_3_2 *n = f->neighbor(k);
      int fi = n->index(f);

      w1 = f->vertex(i);
      w2 = f->vertex(j);

      v1 = w1->info();
      v2 = w2->info();

      if( is_infinite(v1) || is_infinite(v2) ){
	// there will be another ear, so let's ignore this one,
	// because it is complicated to treat
	goto next_edge;
      }       
      w0 = f->vertex(k);
      w3 = n->vertex(fi);

      v0 = w0->info();
      v3 = w3->info();

      bool inf_0 = is_infinite(v0);
      bool inf_3 = is_infinite(v3);

      if( inf_0 || inf_3 || 
	  (orientation(v0->point(), v1->point(), 
		       v2->point(), v3->point()) == POSITIVE) ) {
	// the two faces form a concavity, in which we might plug a cell

	// we now look at all vertices that are on the boundary of the hole
	for(typename Surface::Vertex_iterator vit = surface.vertices_begin();
	    vit != surface.vertices_end();
	    ++vit) {
	  Vertex_handle v = (*vit).info();
	  if( (! is_infinite(v))
	      && (v != v0) && (v != v1) && (v != v2) && (v != v3)) {

	    Bounded_side bs;
	    
	    if (inf_0) {
 	      bs = side_of_sphere_inf_perturb(v2, v1, v3, v);
	    } else if(inf_3) {
	      bs = side_of_sphere_inf_perturb(v0, v1, v2, v);
	    } else {
	      bs = side_of_sphere_finite_perturb(v0,v1,v2,v3,v);
	    }

	    if (bs == ON_BOUNDED_SIDE)
		goto next_edge;
	  }
	}

	// we looked at all vertices

	Face_3_2 *m_i = f->neighbor(i);
	Face_3_2 *m_j = f->neighbor(j); 
	bool neighbor_i = m_i == n->neighbor(cw(fi));
	bool neighbor_j = m_j == n->neighbor(ccw(fi));

	if ((! neighbor_i) && (! neighbor_j)
	       && surface.is_edge(f->vertex(k), n->vertex(fi))) {
	  // The edge that would get introduced is on the surface
	  goto next_edge;
	}
	
	// none of the vertices violates the Delaunay property
	// We are ready to plug a new cell

        Cell_handle ch = (Cell *) _tds.create_cell(
		                       (typename Tds::Vertex *) &(*v0),
	                               (typename Tds::Vertex *) &(*v1),
			               (typename Tds::Vertex *) &(*v2),
			               (typename Tds::Vertex *) &(*v3));
	cells.push_back(ch);

	// The new cell touches the faces that form the ear
	Facet fac = n->info();
	ch->set_neighbor(0, fac.first);
	fac.first->set_neighbor(fac.second, ch);
	fac = f->info();
	ch->set_neighbor(3, fac.first);
	fac.first->set_neighbor(fac.second, ch);

	// It may touch another face, 
	// or even two other faces if it is the last cell
	if(neighbor_i) {
	  fac = m_i->info();
	  ch->set_neighbor(1, fac.first);
	  fac.first->set_neighbor(fac.second, ch);
	}
	if(neighbor_j) {
	  fac = m_j->info();
	  ch->set_neighbor(2, fac.first);
	  fac.first->set_neighbor(fac.second, ch);
	}
	
	if((! neighbor_i) && (! neighbor_j)) {
	  surface.flip(f,k);
	  int fi = n->index(f);
	  int ni = f->index(n);
	  // The flipped edge is not a concavity
	  f->unmark_edge(ni);
	  // The adjacent edges may be a concavity
	  // that is they are candidates for an ear
	  // In the list of faces they get moved behind f
	  f->mark_edge(cw(ni), f);
	  f->mark_edge(ccw(ni), f);
	  n->mark_edge(cw(fi), f);
	  n->mark_edge(ccw(fi), f);

	  f->set_info(std::make_pair(ch,2));
	  n->set_info(std::make_pair(ch,1));
	} else if (neighbor_i && (! neighbor_j)) {
	  surface.remove_degree_3(f->vertex(j), f);
	  // all three edges adjacent to f are 
	  // candidate for an ear
	  f->mark_adjacent_edges();
	  f->set_info(std::make_pair(ch,2));
	} else if ((! neighbor_i) && neighbor_j)  {
	  surface.remove_degree_3(f->vertex(i), f);
	  f->mark_adjacent_edges();
	  f->set_info(std::make_pair(ch,1));
	} else if (surface.number_of_vertices() != 4) {
	  // this should not happen at all  => panic mode, 
	  //clean up, and say that it didn't work
	  CGAL_triangulation_warning_msg(true, "panic");
	  delete_cells(cells);
	  return false;
	} else {
	  // when we leave the function the vertices and faces of the surface
	  // are deleted by the destructor
	  return true;
	}

	// we successfully inserted a cell
	last_op = f; 
	// we have to reconsider all edges incident to f
	k = -1;
      } // else if (inf_0 ...
    }// if(f->edge(k))
  } // for(;;)
}

CGAL_END_NAMESPACE

#endif // CGAL_DELAUNAY_TRIANGULATION_3_H
