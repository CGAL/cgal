// Copyright (c) 1999-2004  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//                 Sylvain Pion <Sylvain.Pion@sophia.inria.fr>

#ifndef CGAL_REGULAR_TRIANGULATION_3_H
#define CGAL_REGULAR_TRIANGULATION_3_H

#include <CGAL/basic.h>

#include <set>

#include <CGAL/Triangulation_short_names_3.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/Triangulation_cell_base_with_hidden_points_3.h>
#include <CGAL/Unique_hash_map.h>

CGAL_BEGIN_NAMESPACE

template < class Gt,
           class Tds = Triangulation_data_structure_3 <
                                   Triangulation_vertex_base_3<Gt>,
                                   Triangulation_cell_base_with_hidden_points_3<Gt> > >
class Regular_triangulation_3
  : public Triangulation_3<Gt,Tds>
{
  typedef Regular_triangulation_3<Gt, Tds>      Self;
  typedef Triangulation_3<Gt,Tds>               Tr_Base;
public:
  typedef Tds                                   Triangulation_data_structure;
  typedef Gt                                    Geom_traits;

  typedef typename Tr_Base::Vertex_handle       Vertex_handle;
  typedef typename Tr_Base::Cell_handle         Cell_handle;
  typedef typename Tr_Base::Vertex              Vertex;
  typedef typename Tr_Base::Cell                Cell;
  typedef typename Tr_Base::Facet               Facet;
  typedef typename Tr_Base::Edge                Edge;

  typedef Triple<Vertex_handle,Vertex_handle,Vertex_handle> Vertex_triple;

  typedef typename Tr_Base::Locate_type         Locate_type;
  typedef typename Tr_Base::Cell_iterator       Cell_iterator;
  typedef typename Tr_Base::Facet_iterator      Facet_iterator;
  typedef typename Tr_Base::Edge_iterator       Edge_iterator;
  typedef typename Tr_Base::Facet_circulator    Facet_circulator;

  typedef typename Tr_Base::Finite_vertices_iterator Finite_vertices_iterator;
  typedef typename Tr_Base::Finite_cells_iterator    Finite_cells_iterator;
  typedef typename Tr_Base::Finite_facets_iterator   Finite_facets_iterator;
  typedef typename Tr_Base::Finite_edges_iterator    Finite_edges_iterator;
  typedef typename Tr_Base::All_cells_iterator       All_cells_iterator;

  typedef typename Gt::Weighted_point              Weighted_point;
  typedef typename Gt::Bare_point                  Bare_point;
  typedef typename Gt::Segment_3                   Segment;
  typedef typename Gt::Triangle_3                  Triangle;
  typedef typename Gt::Tetrahedron_3               Tetrahedron;
  typedef typename Gt::Object_3                    Object;

  //Tag to distinguish Delaunay from Regular triangulations
  typedef Tag_true   Weighted_tag; 

  using Tr_Base::cw;
  using Tr_Base::ccw;
#ifndef CGAL_CFG_USING_BASE_MEMBER_BUG_2
  using Tr_Base::geom_traits;
#endif
  using Tr_Base::number_of_vertices;
  using Tr_Base::dimension;
  using Tr_Base::finite_facets_begin;
  using Tr_Base::finite_facets_end;
  using Tr_Base::finite_vertices_begin;
  using Tr_Base::finite_vertices_end;
  using Tr_Base::finite_cells_begin;
  using Tr_Base::finite_cells_end;
  using Tr_Base::finite_edges_begin;
  using Tr_Base::finite_edges_end;
  using Tr_Base::tds;
  using Tr_Base::infinite_vertex;
  using Tr_Base::next_around_edge;
  using Tr_Base::vertex_triple_index;

  Regular_triangulation_3(const Gt & gt = Gt())
    : Tr_Base(gt)
  {}

  // copy constructor duplicates vertices and cells
  Regular_triangulation_3(const Regular_triangulation_3 & rt)
      : Tr_Base(rt)
  {
      CGAL_triangulation_postcondition( is_valid() );
  }

  //insertion
  template < typename InputIterator >
  Regular_triangulation_3(InputIterator first, InputIterator last,
                          const Gt & gt = Gt())
      : Tr_Base(gt)
  {
      insert(first, last);
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

  Vertex_handle insert(const Weighted_point & p,
	               Cell_handle start = Cell_handle());

  Vertex_handle insert(const Weighted_point & p, Locate_type lt,
	               Cell_handle c, int li, int);

  void remove (Vertex_handle v);

  template < typename InputIterator >
  int remove(InputIterator first, InputIterator beyond)
  {
    int n = number_of_vertices();
    while (first != beyond) {
      remove (*first);
      ++first;
    }
    return n - number_of_vertices();
  }

private:
  void remove_2D(Vertex_handle v);
  //  void make_hole_2D(Vertex_handle v, std::list<Edge_2D> & hole);
  //  void fill_hole_delaunay_2D(std::list<Edge_2D> & hole);

  void make_canonical(Vertex_triple& t) const;

  Vertex_triple make_vertex_triple(const Facet& f) const;

#ifndef CGAL_CFG_NET2003_MATCHING_BUG
  void make_hole_3D(Vertex_handle v, 
		    std::map<Vertex_triple,Facet> &outer_map,
		    std::vector<Cell_handle> &hole);
#else
  void make_hole_3D(Vertex_handle v, 
		    std::map<Vertex_triple,Facet> &outer_map,
		    std::vector<Cell_handle> &hole)
  {
    CGAL_triangulation_expensive_precondition( ! test_dim_down(v) );

    incident_cells(v, std::back_inserter(hole));

    for (typename std::vector<Cell_handle>::iterator cit = hole.begin();
	 cit != hole.end(); ++cit) {
      int indv = (*cit)->index(v);
      Cell_handle opp_cit = (*cit)->neighbor( indv );
      Facet f(opp_cit, opp_cit->index(*cit)); 
      Vertex_triple vt = make_vertex_triple(f);
      make_canonical(vt);
      outer_map[vt] = f;
      for (int i=0; i<4; i++)
	if ( i != indv )
	  (*cit)->vertex(i)->set_cell(opp_cit);
    }
  }
#endif

  void remove_3D(Vertex_handle v);

public:

  // Queries
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
  side_of_power_segment( Cell_handle c, const Weighted_point &p)
    const;

  Vertex_handle
  nearest_power_vertex_in_cell(const Bare_point& p, 
			       const Cell_handle& c)  const;

  Vertex_handle
  nearest_power_vertex(const Bare_point& p, Cell_handle c =
		       Cell_handle()) const;

  bool is_Gabriel(Cell_handle c, int i) const;
  bool is_Gabriel(Cell_handle c, int i, int j) const;
  bool is_Gabriel(const Facet& f)const ;
  bool is_Gabriel(const Edge& e) const;
  bool is_Gabriel(Vertex_handle v) const;


  // Dual functions
  Bare_point dual(Cell_handle c) const;

//   Object dual(const Facet & f) const
//     { return dual( f.first, f.second ); }
//   Object dual(Cell_handle c, int i) const;

  template < class Stream> 		
  Stream& draw_dual(Stream & os)
    {
      typedef typename Gt::Line_3        Line;
      typedef typename Gt::Ray_3         Ray;
      Finite_facets_iterator fit = finite_facets_begin();
      for (; fit != finite_facets_end(); ++fit) {
	Object o = dual(*fit);
	Bare_point p;
	Ray r;
	Segment s;
	if (CGAL::assign(p,o)) os << p;
	if (CGAL::assign(s,o)) os << s;
	if (CGAL::assign(r,o)) os << r; 
      }
      return os;
    }

  bool is_valid(bool verbose = false, int level = 0) const;

private:
  bool
  less_power_distance(const Bare_point &p, 
		      const Weighted_point &q, 
		      const Weighted_point &r)  const 
  {
    return 
      geom_traits().compare_power_distance_3_object()(p, q, r) == SMALLER;
  }

  Bare_point
  construct_weighted_circumcenter(const Weighted_point &p, 
				  const Weighted_point &q, 
				  const Weighted_point &r, 
				  const Weighted_point &s) const
  {
     return geom_traits().construct_weighted_circumcenter_3_object()(p,q,r,s);
  }

  Vertex_handle
  nearest_power_vertex(const Bare_point &p, 
		       Vertex_handle v,
		       Vertex_handle w) const
  {
      // In case of equality, v is returned.
      CGAL_triangulation_precondition(v != w);
      if (is_infinite(v))	  return w;
      if (is_infinite(w))	  return v;
      return less_power_distance(p, w->point(), v->point()) ? w : v;
  }

  Oriented_side
  power_test(const Weighted_point &p, const Weighted_point &q) const
  {
      CGAL_precondition(equal(p, q));
      return geom_traits().power_test_3_object()(p, q);
  }

  Oriented_side
  power_test(const Weighted_point &p, const Weighted_point &q,
	     const Weighted_point &r) const
  {
      CGAL_precondition(collinear(p, q, r));
      return geom_traits().power_test_3_object()(p, q, r);
  }

  Oriented_side
  power_test(const Weighted_point &p, const Weighted_point &q,
	     const Weighted_point &r, const Weighted_point &s) const
  {
      CGAL_precondition(coplanar(p, q, r, s));
      return geom_traits().power_test_3_object()(p, q, r, s);
  }

  Oriented_side
  power_test(const Weighted_point &p, const Weighted_point &q,
	     const Weighted_point &r, const Weighted_point &s,
	     const Weighted_point &t) const
  {
      return geom_traits().power_test_3_object()(p, q, r, s, t);
  }

  bool in_conflict_3(const Weighted_point &p, const Cell_handle c) const
  {
      return side_of_power_sphere(c, p) == ON_BOUNDED_SIDE;
  }

  bool in_conflict_2(const Weighted_point &p, const Cell_handle c, int i) const
  {
      return side_of_power_circle(c, i, p) == ON_BOUNDED_SIDE;
  }

  bool in_conflict_1(const Weighted_point &p, const Cell_handle c) const
  {
      return side_of_power_segment(c, p) == ON_BOUNDED_SIDE;
  }

  bool in_conflict_0(const Weighted_point &p, const Cell_handle c) const
  {
      return power_test(c->vertex(0)->point(), p) == ON_POSITIVE_SIDE;
  }

  class Conflict_tester_3
  {
      const Weighted_point &p;
      const Self *t;
      mutable std::vector<Vertex_handle> cv;

  public:

      Conflict_tester_3(const Weighted_point &pt, const Self *tr)
	  : p(pt), t(tr) {}

      bool operator()(const Cell_handle c) const
      {
	  // We mark the vertices so that we can find the deleted ones easily.
	  if (t->in_conflict_3(p, c))
	  {
	      for (int i=0; i<4; i++)
	      {
		  Vertex_handle v = c->vertex(i);
		  if (v->cell() != Cell_handle())
		  {
		      cv.push_back(v);
		      v->set_cell(Cell_handle());
		  }
	      }
	      return true;
	  }
	  return false;
      }

      std::vector<Vertex_handle> & conflict_vector()
      {
	  return cv;
      }
  };

  class Conflict_tester_2
  {
      const Weighted_point &p;
      const Self *t;
      mutable std::vector<Vertex_handle> cv;

  public:

      Conflict_tester_2(const Weighted_point &pt, const Self *tr)
	  : p(pt), t(tr) {}

      bool operator()(const Cell_handle c) const
      {
	  if (t->in_conflict_2(p, c, 3))
	  {
	      for (int i=0; i<3; i++)
	      {
		  Vertex_handle v = c->vertex(i);
		  if (v->cell() != Cell_handle())
		  {
		      cv.push_back(v);
		      v->set_cell(Cell_handle());
		  }
	      }
	      return true;
	  }
	  return false;
      }

      std::vector<Vertex_handle> & conflict_vector()
      {
	  return cv;
      }
  };

  friend class Conflict_tester_3;
  friend class Conflict_tester_2;
};


template < class Gt, class Tds >
typename Regular_triangulation_3<Gt,Tds>::Vertex_handle
Regular_triangulation_3<Gt,Tds>::
nearest_power_vertex_in_cell(const Bare_point& p, 
			     const Cell_handle& c) const
// Returns the finite vertex of the cell c with smaller
// power distance  to p.
{
    CGAL_triangulation_precondition(dimension() >= 1);
    Vertex_handle nearest = nearest_power_vertex(p, 
						 c->vertex(0), 
						 c->vertex(1));
    if (dimension() >= 2) {
	nearest = nearest_power_vertex(p, nearest, c->vertex(2));
        if (dimension() == 3)
	    nearest = nearest_power_vertex(p, nearest, c->vertex(3));
    }
    return nearest;
}


template < class Gt, class Tds >
typename Regular_triangulation_3<Gt,Tds>::Vertex_handle
Regular_triangulation_3<Gt,Tds>::
nearest_power_vertex(const Bare_point& p, Cell_handle start) const
{
    if (number_of_vertices() == 0)  
      return Vertex_handle();

    // Use a brute-force algorithm if dimension < 3.
    if (dimension() < 3) {
	Finite_vertices_iterator vit = finite_vertices_begin();
	Vertex_handle res = vit;
	for (++vit; vit != finite_vertices_end(); ++vit)
	    res = nearest_power_vertex(p, res, vit);
	return res;
    }

    Locate_type lt;
    int li, lj;
    // I put the cast here temporarily 
    // until we solve the traits class pb of regular triangulation
    Cell_handle c = locate(static_cast<Weighted_point>(p), lt, li, lj, start);
  
    // - start with the closest vertex from the located cell.
    // - repeatedly take the nearest of its incident vertices if any
    // - if not, we're done.
    Vertex_handle nearest = nearest_power_vertex_in_cell(p, c);
    std::vector<Vertex_handle> vs;
    vs.reserve(32);
    while (true) {
	Vertex_handle tmp = nearest;
        incident_vertices(nearest, std::back_inserter(vs));
        for (typename std::vector<Vertex_handle>::const_iterator
		vsit = vs.begin(); vsit != vs.end(); ++vsit)
	    tmp = nearest_power_vertex(p, tmp, *vsit);
	if (tmp == nearest)
	    break;
	vs.clear();
	nearest = tmp;
    }
    return nearest;
}

template < class Gt, class Tds >
typename Regular_triangulation_3<Gt,Tds>::Bare_point
Regular_triangulation_3<Gt,Tds>::
dual(Cell_handle c) const
{
  CGAL_triangulation_precondition(dimension()==3);
  CGAL_triangulation_precondition( ! is_infinite(c) );
  return construct_weighted_circumcenter( c->vertex(0)->point(),
					  c->vertex(1)->point(),
					  c->vertex(2)->point(),
					  c->vertex(3)->point() );
}



template < class Gt, class Tds >
Bounded_side
Regular_triangulation_3<Gt,Tds>::
side_of_power_sphere( Cell_handle c, const Weighted_point &p) const
{
  CGAL_triangulation_precondition( dimension() == 3 );
  int i3;
  if ( ! c->has_vertex( infinite_vertex(), i3 ) ) {
    return Bounded_side( power_test (c->vertex(0)->point(),
				     c->vertex(1)->point(),
				     c->vertex(2)->point(),
				     c->vertex(3)->point(), p) );
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
  Orientation o = orientation(c->vertex(i0)->point(),
		              c->vertex(i1)->point(),
		              c->vertex(i2)->point(), p);
  if (o != ZERO)
    return Bounded_side(o);

  // else p coplanar with i0,i1,i2
  return Bounded_side( power_test( c->vertex(i0)->point(),
				   c->vertex(i1)->point(),
				   c->vertex(i2)->point(), p ) );
}

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
      return Bounded_side( power_test(c->vertex(0)->point(),
				      c->vertex(1)->point(),
				      c->vertex(2)->point(), p) );
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
    // case when p collinear with v1v2
    return Bounded_side( power_test( v1->point(), v2->point(), p ) );
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
    CGAL_triangulation_precondition( coplanar ( c->vertex(i0)->point(),
						c->vertex(i1)->point(),
						c->vertex(i2)->point(), p) );
    return Bounded_side( power_test(c->vertex(i0)->point(),
				    c->vertex(i1)->point(),
				    c->vertex(i2)->point(), p) );
  }
  //else infinite facet
  // v1, v2 finite vertices of the facet such that v1,v2,infinite
  // is positively oriented
  Vertex_handle v1 = c->vertex( next_around_edge(i3,i) ),
                v2 = c->vertex( next_around_edge(i,i3) );
  Orientation o = (Orientation)
                  (coplanar_orientation( v1->point(), v2->point(),
					c->vertex(i)->point()) *
                  coplanar_orientation( v1->point(), v2->point(), p));
  // then the code is duplicated from 2d case
  if ( o != ZERO )
      return Bounded_side( -o );
  // because p is in f iff
  // it is not on the same side of v1v2 as c->vertex(i)
  // case when p collinear with v1v2 :
  return Bounded_side( power_test( v1->point(), v2->point(), p ) );
}

template < class Gt, class Tds >
Bounded_side
Regular_triangulation_3<Gt,Tds>::
side_of_power_segment( Cell_handle c, const Weighted_point &p) const
{
  CGAL_triangulation_precondition( dimension() == 1 );
  if ( ! is_infinite(c,0,1) )
    return Bounded_side( power_test( c->vertex(0)->point(),
				     c->vertex(1)->point(), p ) );
  Locate_type lt; int i;
  Bounded_side soe = side_of_edge( p, c, lt, i );
  if (soe != ON_BOUNDARY)
    return soe;
  // Either we compare weights, or we use the finite neighboring edge
  Cell_handle finite_neighbor = c->neighbor(c->index(infinite_vertex()));
  CGAL_assertion(!is_infinite(finite_neighbor,0,1));
  return Bounded_side( power_test( finite_neighbor->vertex(0)->point(),
			           finite_neighbor->vertex(1)->point(), p ) );
}

template < class Gt, class Tds >
bool
Regular_triangulation_3<Gt,Tds>::
is_Gabriel(const Facet& f) const
{
  return is_Gabriel(f.first, f.second);
}

template < class Gt, class Tds >
bool
Regular_triangulation_3<Gt,Tds>::
is_Gabriel(Cell_handle c, int i) const
{
  CGAL_triangulation_precondition(dimension() == 3 && !is_infinite(c,i));
  typename Geom_traits::Side_of_bounded_orthogonal_sphere_3
    side_of_bounded_orthogonal_sphere = 
    geom_traits().side_of_bounded_orthogonal_sphere_3_object();

  if ((!is_infinite(c->vertex(i))) &&
      side_of_bounded_orthogonal_sphere(
	 c->vertex(vertex_triple_index(i,0))->point(),
	 c->vertex(vertex_triple_index(i,1))->point(),
	 c->vertex(vertex_triple_index(i,2))->point(),
	 c->vertex(i)->point()) == ON_BOUNDED_SIDE ) return false;

  Cell_handle neighbor = c->neighbor(i);
  int in = neighbor->index(c);

  if ((!is_infinite(neighbor->vertex(in))) &&
      side_of_bounded_orthogonal_sphere(
	 c->vertex(vertex_triple_index(i,0))->point(),
	 c->vertex(vertex_triple_index(i,1))->point(),
	 c->vertex(vertex_triple_index(i,2))->point(),	
	 neighbor->vertex(in)->point()) == ON_BOUNDED_SIDE ) return false;
 
  return true;
}

  
template < class Gt, class Tds >
bool
Regular_triangulation_3<Gt,Tds>::
is_Gabriel(const Edge& e) const
{
  return is_Gabriel(e.first, e.second, e.third);
}

template < class Gt, class Tds >
bool
Regular_triangulation_3<Gt,Tds>::
is_Gabriel(Cell_handle c, int i, int j) const
{
  CGAL_triangulation_precondition(dimension() == 3 && !is_infinite(c,i,j));
  typename Geom_traits::Side_of_bounded_orthogonal_sphere_3
    side_of_bounded_orthogonal_sphere = 
    geom_traits().side_of_bounded_orthogonal_sphere_3_object();

  Facet_circulator fcirc = incident_facets(c,i,j),
                   fdone(fcirc);
  Vertex_handle v1 = c->vertex(i);
  Vertex_handle v2 = c->vertex(j);
  do {
      // test whether the vertex of cc opposite to *fcirc
      // is inside the sphere defined by the edge e = (s, i,j)
      Cell_handle cc = (*fcirc).first;
      int ii = (*fcirc).second;
      if (!is_infinite(cc->vertex(ii)) &&
	  side_of_bounded_orthogonal_sphere( v1->point(), 
					 v2->point(),
					 cc->vertex(ii)->point())  
	  == ON_BOUNDED_SIDE ) return false;
  } while(++fcirc != fdone);
  return true;
}

template < class Gt, class Tds >
bool
Regular_triangulation_3<Gt,Tds>::
is_Gabriel(Vertex_handle v) const
{
  return nearest_power_vertex( v->point().point(), v->cell()) == v;
}

template < class Gt, class Tds >
typename Regular_triangulation_3<Gt,Tds>::Vertex_handle
Regular_triangulation_3<Gt,Tds>::
insert(const Weighted_point & p, Cell_handle start)
{
    Locate_type lt;
    int li, lj;
    Cell_handle c = locate(p, lt, li, lj, start);
    return insert(p, lt, c, li, lj);
}

template < class Gt, class Tds >
typename Regular_triangulation_3<Gt,Tds>::Vertex_handle
Regular_triangulation_3<Gt,Tds>::
insert(const Weighted_point & p, Locate_type lt, Cell_handle c, int li, int)
{
  switch (dimension()) {
  case 3:
    {
      // TODO :
      // In case the point is completely equal (including weight), then we need
      // to discard it (don't update the triangulation, nor hide it), right ?
      if (! in_conflict_3(p, c)) {  // new point is hidden
	c->hide_point(p);
	return Vertex_handle();
      }

      // Should I mark c's vertices too ?
      Conflict_tester_3 tester(p, this);
      Vertex_handle v = insert_conflict_3(c, tester);
      v->set_point(p);
      for( typename std::vector<Vertex_handle>::iterator
		it = tester.conflict_vector().begin();
		it != tester.conflict_vector().end(); ++it)
      {
	if ((*it)->cell() == Cell_handle()) {
	  // remember the hidden point
	  Cell_handle hider = locate ((*it)->point(), v->cell());
	  hider->hide_point ((*it)->point());
	  // vertex has to be deleted
	  tds().delete_vertex(*it);
	}
      }
      return v;
    }
  case 2:
    {
      switch (lt) {
      case Tr_Base::OUTSIDE_CONVEX_HULL:
      case Tr_Base::FACET:
      case Tr_Base::EDGE:
      case Tr_Base::VERTEX:
	{
          if (! in_conflict_2(p, c, 3)) {  // new point is hidden
	    c->hide_point (p); // remember the point
	    return Vertex_handle();
          }

	  Conflict_tester_2 tester(p, this);
	  Vertex_handle v = insert_conflict_2(c, tester);
	  v->set_point(p);
          for( typename std::vector<Vertex_handle>::iterator
		it = tester.conflict_vector().begin();
		it != tester.conflict_vector().end(); ++it)
	  {
            if ((*it)->cell() == Cell_handle())
	    {
	      // remember the hidden point
	      Cell_handle hider = locate ((*it)->point(), v->cell());
	      hider->hide_point ((*it)->point());
              // vertex has to be deleted
              tds().delete_vertex(*it);
	    }
	  }
	  return v;
	}
      case Tr_Base::OUTSIDE_AFFINE_HULL:
	{
	  // if the 2d triangulation is Regular, the 3d
	  // triangulation will be Regular
	  return Tr_Base::insert_outside_affine_hull(p);
	}
      default:
	CGAL_triangulation_assertion(false); // CELL cannot happen in 2D.
      }
    }//dim 2
  case 1:
    {
      switch (lt) {
      case Tr_Base::OUTSIDE_CONVEX_HULL:
      case Tr_Base::EDGE:
      case Tr_Base::VERTEX:
	{
          if (! in_conflict_1(p, c)) {  // new point is hidden
	    c->hide_point (p); // remember the point
	    return Vertex_handle();
          }

	  Cell_handle bound[2];
          // corresponding index: bound[j]->neighbor(1-j) is in conflict.
	  std::vector<Vertex_handle>  hidden_vertices;
	  std::vector<Cell_handle>    conflicts;
          conflicts.push_back(c);

          // We get all cells in conflict,
          // and remember the 2 external boundaries.

	  for (int j = 0; j<2; ++j) {
	    Cell_handle n = c->neighbor(j);
	    while ( in_conflict_1( p, n) ) {
	      conflicts.push_back(n);
              hidden_vertices.push_back(n->vertex(j));
	      n = n->neighbor(j);
	    }
	    bound[j] = n;
	  }

          // We preserve the order (like the orientation in 2D-3D).

	  Vertex_handle v = tds().create_vertex();
	  v->set_point(p);
          Cell_handle c0 = tds().create_face(v, bound[0]->vertex(0),
		                             Vertex_handle());
          Cell_handle c1 = tds().create_face(bound[1]->vertex(1), v,
		                             Vertex_handle());
          tds().set_adjacency(c0, 1, c1, 0);
          tds().set_adjacency(bound[0], 1, c0, 0);
          tds().set_adjacency(c1, 1, bound[1], 0);
          bound[0]->vertex(0)->set_cell(bound[0]);
          bound[1]->vertex(1)->set_cell(bound[1]);
          v->set_cell(c0);

	  for (typename std::vector<Vertex_handle>::iterator i = hidden_vertices.begin();
		 i != hidden_vertices.end(); ++i) {
	    Cell_handle hider = locate ((*i)->point(), c0);
	    hider->hide_point ((*i)->point());
	  }

	  tds().delete_cells(conflicts.begin(), conflicts.end());
	  tds().delete_vertices(hidden_vertices.begin(), hidden_vertices.end());
	  return v;
	}
      case Tr_Base::OUTSIDE_AFFINE_HULL:
	return Tr_Base::insert_outside_affine_hull(p);
      case Tr_Base::FACET:
      case Tr_Base::CELL:
	// impossible in dimension 1
        CGAL_assertion(false);
	return Vertex_handle();
      }
    }
  case 0:
    {
        // We need to compare the weights when the points are equal.
        if (lt == Tr_Base::VERTEX) {
	  CGAL_assertion(li == 0);
	  if (in_conflict_0(p, c)) {
	    c->hide_point (c->vertex(li)->point());
            c->vertex(li)->set_point(p); // replace by heavier point
	  } else {
	    c->hide_point (p); // hide new point
	  }
        } else
            return Tr_Base::insert(p, c);
    }
  default :
    {
      return Tr_Base::insert(p, c);
    }
  }
}

template < class Gt, class Tds >
void
Regular_triangulation_3<Gt,Tds>::
remove_2D(Vertex_handle v)
{
  // Not yet implemented
  std::cerr << "WARNING: RT3::remove() in 2D not implemented" << std::endl;
}

#ifndef CGAL_CFG_NET2003_MATCHING_BUG
template < class Gt, class Tds >
void
Regular_triangulation_3<Gt,Tds>::
make_hole_3D (Vertex_handle v, 
	      std::map<Vertex_triple,Facet>& outer_map,
	      std::vector<Cell_handle> & hole)
{
  CGAL_triangulation_expensive_precondition( ! test_dim_down(v) );

  incident_cells(v, std::back_inserter(hole));

  for (typename std::vector<Cell_handle>::iterator cit = hole.begin();
       cit != hole.end(); ++cit) {
    int indv = (*cit)->index(v);
    Cell_handle opp_cit = (*cit)->neighbor( indv );
    Facet f(opp_cit, opp_cit->index(*cit)); 
    Vertex_triple vt = make_vertex_triple(f);
    make_canonical(vt);
    outer_map[vt] = f;
    for (int i=0; i<4; i++)
      if ( i != indv )
	(*cit)->vertex(i)->set_cell(opp_cit);
  }
}
#endif

template < class Gt, class Tds >
void
Regular_triangulation_3<Gt,Tds>::
make_canonical(Vertex_triple& t) const
{
  int i = (&*(t.first) < &*(t.second))? 0 : 1;
  if(i==0) {
    i = (&*(t.first) < &*(t.third))? 0 : 2;
  } else {
    i = (&*(t.second) < &*(t.third))? 1 : 2;
  }
  Vertex_handle tmp; 
  switch(i){
  case 0: return;
  case 1:
    tmp = t.first;
    t.first = t.second;
    t.second = t.third;
    t.third = tmp;
    return;
  default:
    tmp = t.first;
    t.first = t.third;
    t.third = t.second;
    t.second = tmp;
  }
}

template < class Gt, class Tds >
typename Regular_triangulation_3<Gt,Tds>::Vertex_triple
Regular_triangulation_3<Gt,Tds>::
make_vertex_triple(const Facet& f) const
{
  // static const int vertex_triple_index[4][3] = { {1, 3, 2}, {0, 2, 3},
//                                                  {0, 3, 1}, {0, 1, 2} };
  Cell_handle ch = f.first;
  int i = f.second;
  
  return Vertex_triple(ch->vertex(vertex_triple_index(i,0)),
		       ch->vertex(vertex_triple_index(i,1)),
		       ch->vertex(vertex_triple_index(i,2))); 
}

template < class Gt, class Tds >
void
Regular_triangulation_3<Gt,Tds>::
remove_3D(Vertex_handle v)
{
  std::vector<Cell_handle> hole;
  hole.reserve(64);

  // Construct the set of vertex triples on the boundary
  // with the facet just behind
  typedef std::map<Vertex_triple,Facet> Vertex_triple_Facet_map;
  Vertex_triple_Facet_map outer_map;
  Vertex_triple_Facet_map inner_map;

  make_hole_3D (v, outer_map, hole);

  bool inf = false;
  unsigned int i;
  // collect all vertices on the boundary
  std::vector<Vertex_handle> vertices;
  vertices.reserve(64);

  incident_vertices(v, std::back_inserter(vertices));
  
  // create a Regular triangulation of the points on the boundary
  // and make a map from the vertices in aux towards the vertices in *this
  Self aux;

  Unique_hash_map<Vertex_handle,Vertex_handle> vmap;

  Cell_handle ch = Cell_handle();
  for(i=0; i < vertices.size(); i++){
    if(! is_infinite(vertices[i])){
      Vertex_handle vh = aux.insert(vertices[i]->point(), ch);
      ch = vh->cell();
      vmap[vh] = vertices[i];
    }else {
      inf = true;
    }
  }

  if(aux.dimension()==2){
    Vertex_handle fake_inf = aux.insert(v->point());
    vmap[fake_inf] = infinite_vertex();
  } else {
    vmap[aux.infinite_vertex()] = infinite_vertex();
  }

  CGAL_triangulation_assertion(aux.dimension() == 3);

  // Construct the set of vertex triples of aux
  // We reorient the vertex triple so that it matches those from outer_map
  // Also note that we use the vertices of *this, not of aux
  
  if(inf){
    for(All_cells_iterator it = aux.all_cells_begin();
	it != aux.all_cells_end();
	++it){
      for(i=0; i < 4; i++){
	Facet f = std::pair<Cell_handle,int>(it,i);
	Vertex_triple vt_aux = make_vertex_triple(f);
	Vertex_triple vt(vmap[vt_aux.first],vmap[vt_aux.third],vmap[vt_aux.second]);
	make_canonical(vt);
	inner_map[vt]= f;
      }
    }
  } else {
      for(Finite_cells_iterator it = aux.finite_cells_begin();
	it != aux.finite_cells_end();
	++it){
      for(i=0; i < 4; i++){
	Facet f = std::pair<Cell_handle,int>(it,i);
	Vertex_triple vt_aux = make_vertex_triple(f);
	Vertex_triple vt(vmap[vt_aux.first],vmap[vt_aux.third],vmap[vt_aux.second]);
	make_canonical(vt);
	inner_map[vt]= f;
      }
    }
  }
  // Grow inside the hole, by extending the surface
  while(! outer_map.empty()){
    typename Vertex_triple_Facet_map::iterator oit = outer_map.begin();
    while(is_infinite(oit->first.first) ||
	  is_infinite(oit->first.second) ||
	  is_infinite(oit->first.third)){
      ++oit;
      // otherwise the lookup in the inner_map fails
      // because the infinite vertices are different
    }
    typename Vertex_triple_Facet_map::value_type o_vt_f_pair = *oit;
    Cell_handle o_ch = o_vt_f_pair.second.first;
    unsigned int o_i = o_vt_f_pair.second.second;

    typename Vertex_triple_Facet_map::iterator iit =
             inner_map.find(o_vt_f_pair.first);
    CGAL_triangulation_assertion(iit != inner_map.end());
    typename Vertex_triple_Facet_map::value_type i_vt_f_pair = *iit;
    Cell_handle i_ch = i_vt_f_pair.second.first;
    unsigned int i_i = i_vt_f_pair.second.second;
    
    // create a new cell and glue it to the outer surface
    Cell_handle new_ch = tds().create_cell();
    new_ch->set_vertices(vmap[i_ch->vertex(0)], vmap[i_ch->vertex(1)],
			 vmap[i_ch->vertex(2)], vmap[i_ch->vertex(3)]);
    
    o_ch->set_neighbor(o_i,new_ch);
    new_ch->set_neighbor(i_i, o_ch);

    // for the other faces check, if they can also be glued
    for(i = 0; i < 4; i++){
      if(i != i_i){
	Facet f = std::pair<Cell_handle,int>(new_ch,i);
	Vertex_triple vt = make_vertex_triple(f);
	make_canonical(vt);
	std::swap(vt.second,vt.third);
	typename Vertex_triple_Facet_map::iterator oit2 = outer_map.find(vt);
	if(oit2 == outer_map.end()){
	  std::swap(vt.second,vt.third);
	  outer_map[vt]= f;
	} else {
	  // glue the faces
	  typename Vertex_triple_Facet_map::value_type o_vt_f_pair2 = *oit2;
	  Cell_handle o_ch2 = o_vt_f_pair2.second.first;
	  int o_i2 = o_vt_f_pair2.second.second;
	  o_ch2->set_neighbor(o_i2,new_ch);
	  new_ch->set_neighbor(i, o_ch2);
	  outer_map.erase(oit2);
	}
      }
    }
    outer_map.erase(oit);
  }

  // reinsert hidden points
  typename std::vector<Cell_handle>::iterator hi, hend;
  for (hi = hole.begin(), hend = hole.end(); hi != hend; ++hi) {
    int hole_i = (*hi)->index(v);
    int out_i = (*hi)->mirror_index (hole_i);
    Cell_handle out_ch = (*hi)->neighbor (hole_i);
    typename Cell::Point_iterator pi, pend;
    for (pi = (*hi)->hidden_points_begin(), pend = (*hi)->hidden_points_end();
	 pi != pend; ++pi) {
      insert (*pi, out_ch->neighbor (out_i));
    }
  }

  tds().delete_vertex(v);
  tds().delete_cells(hole.begin(), hole.end());
}

template < class Gt, class Tds >
void
Regular_triangulation_3<Gt,Tds>::
remove(Vertex_handle v)
{
  CGAL_triangulation_precondition( v != Vertex_handle());
  CGAL_triangulation_precondition( !is_infinite(v));
  CGAL_triangulation_expensive_precondition( tds().is_vertex(v) );

  if (dimension() >= 0 && test_dim_down(v)) {
      // collect all the hidden points
      std::vector<Weighted_point> hidden;
      Finite_cells_iterator ci, cend;
      for (ci = finite_cells_begin(), cend = finite_cells_end(); ci != cend; ++ci) {
	typename Cell::Point_iterator hi, hend;
	for (hi = ci->hidden_points_begin(), hend = ci->hidden_points_end();
	     hi != hend; ++hi)
	  hidden.push_back(*hi);
      }
      tds().remove_decrease_dimension(v);
      // Now try to see if we need to re-orient.
      if (dimension() == 2) {
	  Facet f = *finite_facets_begin();
          if (coplanar_orientation(f.first->vertex(0)->point(),
		                   f.first->vertex(1)->point(),
				   f.first->vertex(2)->point()) == NEGATIVE)
	      tds().reorient();
      }
      // reinsert the hidden points
      insert (hidden.begin(), hidden.end());
	     
      CGAL_triangulation_expensive_postcondition(is_valid());
      return;
  }

  if (dimension() == 1) {
      tds().remove_from_maximal_dimension_simplex(v);
      CGAL_triangulation_expensive_postcondition(is_valid());
      return;
  }

  if (dimension() == 2) {
      remove_2D(v);
      CGAL_triangulation_expensive_postcondition(is_valid());
      return;
  }

  CGAL_triangulation_assertion( dimension() == 3 );

  remove_3D(v);

  CGAL_triangulation_expensive_postcondition(is_valid());
}

template < class Gt, class Tds >
bool
Regular_triangulation_3<Gt,Tds>::
is_valid(bool verbose, int level) const
{
  if ( ! Tr_Base::is_valid(verbose,level) ) {
    if (verbose)
	std::cerr << "invalid base triangulation" << std::endl;
    CGAL_triangulation_assertion(false);
    return false;
  }

  switch ( dimension() ) {
  case 3:
    {
      Finite_cells_iterator it;
      for ( it = finite_cells_begin(); it != finite_cells_end(); ++it ) {
	is_valid_finite(it, verbose, level);
	for (int i=0; i<4; i++ ) {
	  if ( !is_infinite
	       (it->neighbor(i)->vertex(it->neighbor(i)->index(it))) ) {
	    if ( side_of_power_sphere
		 (it,
		  it->neighbor(i)->vertex(it->neighbor(i)->index(it))->point())
		  == ON_BOUNDED_SIDE ) {
	      if (verbose)
		std::cerr << "non-empty sphere " << std::endl;
	      CGAL_triangulation_assertion(false);
	      return false;
	    }
	  }
	}
      }
      break;
    }
  case 2:
    {
      Finite_facets_iterator it;
      for ( it = finite_facets_begin(); it != finite_facets_end(); ++it ) {
	is_valid_finite((*it).first, verbose, level);
	for (int i=0; i<3; i++ ) {
	  if( !is_infinite
	      ((*it).first->neighbor(i)->vertex( (((*it).first)->neighbor(i))
						 ->index((*it).first))) ) {
	    if ( side_of_power_circle
		 ( (*it).first, 3,
		   (*it).first->neighbor(i)->
		   vertex( (((*it).first)->neighbor(i))
			   ->index((*it).first) )->point() )
		 == ON_BOUNDED_SIDE ) {
	      if (verbose)
		std::cerr << "non-empty circle " << std::endl;
	      CGAL_triangulation_assertion(false);
	      return false;
	    }
	  }
	}
      }
      break;
    }
  case 1:
    {
      Finite_edges_iterator it;
      for ( it = finite_edges_begin(); it != finite_edges_end(); ++it ) {
	is_valid_finite((*it).first, verbose, level);
	for (int i=0; i<2; i++ ) {
	  if( !is_infinite
	      ((*it).first->neighbor(i)->vertex( (((*it).first)->neighbor(i))
						 ->index((*it).first))) ) {
	    if ( side_of_power_segment
		 ( (*it).first,
		   (*it).first->neighbor(i)->
		   vertex( (((*it).first)->neighbor(i))
			   ->index((*it).first) )->point() )
		 == ON_BOUNDED_SIDE ) {
	      if (verbose)
		std::cerr << "non-empty edge " << std::endl;
	      CGAL_triangulation_assertion(false);
	      return false;
	    }
	  }
	}
      }
      break;
    }
  }
  if (verbose)
      std::cerr << "valid Regular triangulation" << std::endl;
  return true;
}

CGAL_END_NAMESPACE

#endif // CGAL_REGULAR_TRIANGULATION_3_H
