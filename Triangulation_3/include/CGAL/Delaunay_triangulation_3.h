// Copyright (c) 1999-2004   INRIA Sophia-Antipolis (France).
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
// $URL$
// $Id$
//
//
// Author(s)     : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//                 Sylvain Pion
//                 Andreas Fabri <Andreas.Fabri@sophia.inria.fr>

#ifndef CGAL_DELAUNAY_TRIANGULATION_3_H
#define CGAL_DELAUNAY_TRIANGULATION_3_H

#include <CGAL/basic.h>

#include <utility>
#include <vector>

#include <CGAL/Triangulation_3.h>
#include <CGAL/internal/Delaunay_remove_tds_3.h>
#include <CGAL/iterator.h>
#include <CGAL/Location_policy.h>

CGAL_BEGIN_NAMESPACE

template < class Gt,
           class Tds_ = Default,
           class Location_policy = Default >
class Delaunay_triangulation_3;


template < class Gt, class Tds_ >
class Delaunay_triangulation_3<Gt, Tds_>
  : public Triangulation_3<Gt, Tds_>
{
  typedef Delaunay_triangulation_3<Gt, Tds_> Self;
  typedef Triangulation_3<Gt,Tds_>           Tr_Base;

public:

  typedef typename Tr_Base::Triangulation_data_structure
                                     Triangulation_data_structure;
  typedef Gt                         Geom_traits;
  typedef Compact_location           Location_policy;

  typedef typename Gt::Point_3       Point;
  typedef typename Gt::Segment_3     Segment;
  typedef typename Gt::Triangle_3    Triangle;
  typedef typename Gt::Tetrahedron_3 Tetrahedron;

  // types for dual:
  typedef typename Gt::Line_3        Line;
  typedef typename Gt::Ray_3         Ray;
  //typedef typename Gt::Plane_3       Plane;
  typedef typename Gt::Object_3      Object;

  typedef typename Tr_Base::Cell_handle   Cell_handle;
  typedef typename Tr_Base::Vertex_handle Vertex_handle;

  typedef typename Tr_Base::Cell   Cell;
  typedef typename Tr_Base::Vertex Vertex;
  typedef typename Tr_Base::Facet  Facet;
  typedef typename Tr_Base::Edge   Edge;

  typedef typename Tr_Base::Cell_circulator  Cell_circulator;
  typedef typename Tr_Base::Facet_circulator Facet_circulator;
  typedef typename Tr_Base::Cell_iterator    Cell_iterator;
  typedef typename Tr_Base::Facet_iterator   Facet_iterator;
  typedef typename Tr_Base::Edge_iterator    Edge_iterator;
  typedef typename Tr_Base::Vertex_iterator  Vertex_iterator;

  typedef typename Tr_Base::Finite_vertices_iterator Finite_vertices_iterator;
  typedef typename Tr_Base::Finite_cells_iterator    Finite_cells_iterator;
  typedef typename Tr_Base::Finite_facets_iterator   Finite_facets_iterator;
  typedef typename Tr_Base::Finite_edges_iterator    Finite_edges_iterator;

  typedef typename Tr_Base::All_cells_iterator       All_cells_iterator;

  typedef typename Tr_Base::Locate_type Locate_type;


#ifndef CGAL_CFG_USING_BASE_MEMBER_BUG_2
  using Tr_Base::cw;
  using Tr_Base::ccw;
  using Tr_Base::geom_traits;
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
  using Tr_Base::mirror_vertex;
  using Tr_Base::coplanar;
  using Tr_Base::coplanar_orientation;
  using Tr_Base::orientation;
#endif

protected:

  Oriented_side
  side_of_oriented_sphere(const Point &p0, const Point &p1, const Point &p2,
	 const Point &p3, const Point &t, bool perturb = false) const;

  Bounded_side
  coplanar_side_of_bounded_circle(const Point &p, const Point &q,
		  const Point &r, const Point &s, bool perturb = false) const;

  // for dual:
  Point
  construct_circumcenter(const Point &p, const Point &q, const Point &r) const
  {
      return geom_traits().construct_circumcenter_3_object()(p, q, r);
  }

  Line
  construct_equidistant_line(const Point &p1, const Point &p2,
                             const Point &p3) const
  {
      return geom_traits().construct_equidistant_line_3_object()(p1, p2, p3);
  }

  Ray
  construct_ray(const Point &p, const Line &l) const
  {
      return geom_traits().construct_ray_3_object()(p, l);
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

  bool
  less_distance(const Point &p, const Point &q, const Point &r) const
  {
      return geom_traits().compare_distance_3_object()(p, q, r) == SMALLER;
  }

public:

  Delaunay_triangulation_3(const Gt& gt = Gt())
    : Tr_Base(gt)
  {}

  // copy constructor duplicates vertices and cells
  Delaunay_triangulation_3(const Delaunay_triangulation_3 & tr)
    : Tr_Base(tr)
  {
    CGAL_triangulation_postcondition( is_valid() );
  }

  template < typename InputIterator >
  Delaunay_triangulation_3(InputIterator first, InputIterator last,
                           const Gt& gt = Gt())
    : Tr_Base(gt)
  {
      insert(first, last);
  }

  template < class InputIterator >
  int
  insert(InputIterator first, InputIterator last)
  {
    int n = number_of_vertices();

    std::vector<Point> points (first, last);
    std::random_shuffle (points.begin(), points.end());
    spatial_sort (points.begin(), points.end(), geom_traits());

    Vertex_handle hint;
    for (typename std::vector<Point>::const_iterator p = points.begin(), end = points.end();
            p != end; ++p)
        hint = insert(*p, hint);

    return number_of_vertices() - n;
  }

  Vertex_handle insert(const Point & p, Vertex_handle hint)
  {
    return insert(p, hint == Vertex_handle() ? this->infinite_cell() : hint->cell());
  }

  Vertex_handle insert(const Point & p, Cell_handle start = Cell_handle());

  Vertex_handle insert(const Point & p, Locate_type lt,
	               Cell_handle c, int li, int);

  Vertex_handle move_point(Vertex_handle v, const Point & p);

  template <class OutputIteratorBoundaryFacets,
            class OutputIteratorCells,
            class OutputIteratorInternalFacets>
  Triple<OutputIteratorBoundaryFacets,
         OutputIteratorCells,
         OutputIteratorInternalFacets>
  find_conflicts(const Point &p, Cell_handle c,
	         OutputIteratorBoundaryFacets bfit,
                 OutputIteratorCells cit,
		 OutputIteratorInternalFacets ifit) const
  {
      CGAL_triangulation_precondition(dimension() >= 2);

      std::vector<Cell_handle> cells;
      cells.reserve(32);
      std::vector<Facet> facets;
      facets.reserve(64);

      if (dimension() == 2) {
          Conflict_tester_2 tester(p, this);
	  ifit = Tr_Base::find_conflicts
	    (c, tester,
	     make_triple(std::back_inserter(facets),
			 std::back_inserter(cells),
			 ifit)).third;
      }
      else {
          Conflict_tester_3 tester(p, this);
	  ifit = Tr_Base::find_conflicts
	    (c, tester,
	     make_triple(std::back_inserter(facets),
			 std::back_inserter(cells),
			 ifit)).third;
      }

      // Reset the conflict flag on the boundary.
      for(typename std::vector<Facet>::iterator fit=facets.begin();
          fit != facets.end(); ++fit) {
        fit->first->neighbor(fit->second)->tds_data().clear();
	*bfit++ = *fit;
      }

      // Reset the conflict flag in the conflict cells.
      for(typename std::vector<Cell_handle>::iterator ccit=cells.begin();
        ccit != cells.end(); ++ccit) {
        (*ccit)->tds_data().clear();
	*cit++ = *ccit;
      }
      return make_triple(bfit, cit, ifit);
  }

  template <class OutputIteratorBoundaryFacets, class OutputIteratorCells>
  std::pair<OutputIteratorBoundaryFacets, OutputIteratorCells>
  find_conflicts(const Point &p, Cell_handle c,
	         OutputIteratorBoundaryFacets bfit,
                 OutputIteratorCells cit) const
  {
      Triple<OutputIteratorBoundaryFacets,
             OutputIteratorCells,
	     Emptyset_iterator> t = find_conflicts(p, c, bfit, cit,
		                                   Emptyset_iterator());
      return std::make_pair(t.first, t.second);
  }

  // Returns the vertices on the boundary of the conflict hole.
  template <class OutputIterator>
  OutputIterator
  vertices_in_conflict(const Point&p, Cell_handle c, OutputIterator res) const
  {
      CGAL_triangulation_precondition(dimension() >= 2);

      // Get the facets on the boundary of the hole.
      std::vector<Facet> facets;
      find_conflicts(p, c, std::back_inserter(facets),
	             Emptyset_iterator(), Emptyset_iterator());

      // Then extract uniquely the vertices.
      std::set<Vertex_handle> vertices;
      if (dimension() == 3) {
          for (typename std::vector<Facet>::const_iterator i = facets.begin();
	       i != facets.end(); ++i) {
	      vertices.insert(i->first->vertex((i->second+1)&3));
	      vertices.insert(i->first->vertex((i->second+2)&3));
	      vertices.insert(i->first->vertex((i->second+3)&3));
          }
      } else {
          for (typename std::vector<Facet>::const_iterator i = facets.begin();
	       i != facets.end(); ++i) {
	      vertices.insert(i->first->vertex(cw(i->second)));
	      vertices.insert(i->first->vertex(ccw(i->second)));
          }
      }

      return std::copy(vertices.begin(), vertices.end(), res);
  }

  // We return bool only for backward compatibility (it's always true).
  // The documentation mentions void.
  bool remove(Vertex_handle v);

  template < typename InputIterator >
  int remove(InputIterator first, InputIterator beyond)
  {
    int n = number_of_vertices();
    while (first != beyond) {
      remove(*first);
      ++first;
    }
    return n - number_of_vertices();
  }

private:
  void remove_3D_ear(Vertex_handle v);

  Bounded_side
  side_of_sphere(Vertex_handle v0, Vertex_handle v1,
		 Vertex_handle v2, Vertex_handle v3,
		 const Point &p, bool perturb) const;
public:

  // Queries
  Bounded_side
  side_of_sphere(Cell_handle c, const Point & p,
	         bool perturb = false) const
  {
      return side_of_sphere(c->vertex(0), c->vertex(1),
                            c->vertex(2), c->vertex(3), p, perturb);
  }

  Bounded_side
  side_of_circle( const Facet & f, const Point & p, bool perturb = false) const
  {
      return side_of_circle(f.first, f.second, p, perturb);
  }

  Bounded_side
  side_of_circle( Cell_handle c, int i, const Point & p,
	          bool perturb = false) const;

  Vertex_handle
  nearest_vertex_in_cell(const Point& p, Cell_handle c) const;

  Vertex_handle
  nearest_vertex(const Point& p, Cell_handle c = Cell_handle()) const;

  bool is_Gabriel(Cell_handle c, int i) const;
  bool is_Gabriel(Cell_handle c, int i, int j) const;
  bool is_Gabriel(const Facet& f)const ;
  bool is_Gabriel(const Edge& e) const;

// Dual functions
  Point dual(Cell_handle c) const;

  Object dual(const Facet & f) const
  { return dual( f.first, f.second ); }

  Object dual(Cell_handle c, int i) const;

  Line dual_support(Cell_handle c, int i) const;

  bool is_valid(bool verbose = false, int level = 0) const;

  bool is_valid(Cell_handle c, bool verbose = false, int level = 0) const;

  template < class Stream> 		
  Stream& draw_dual(Stream & os)
    {
      for (Finite_facets_iterator fit = finite_facets_begin(),
                                  end = finite_facets_end();
           fit != end; ++fit) {
	Object o = dual(*fit);
	if      (const Segment *s = object_cast<Segment>(&o)) os << *s;
	else if (const Ray *r     = object_cast<Ray>(&o))     os << *r;
	else if (const Point *p   = object_cast<Point>(&o))   os << *p;
      }
      return os;
    }

private:

  Vertex_handle
  nearest_vertex(const Point &p, Vertex_handle v, Vertex_handle w) const
  {
      // In case of equality, v is returned.
      CGAL_triangulation_precondition(v != w);

      if (is_infinite(v))
	  return w;
      if (is_infinite(w))
	  return v;
      return less_distance(p, w->point(), v->point()) ? w : v;
  }

#ifndef CGAL_CFG_NET2003_MATCHING_BUG
  void make_hole_3D_ear( Vertex_handle v,
	                 std::vector<Facet> & boundhole,
	                 std::vector<Cell_handle> & hole);
#else
  void make_hole_3D_ear( Vertex_handle v,
	                 std::vector<Facet> & boundhole,
                         std::vector<Cell_handle> & hole)
  {
    CGAL_triangulation_expensive_precondition( ! test_dim_down(v) );
    incident_cells(v, std::back_inserter(hole));

    for (typename std::vector<Cell_handle>::iterator cit = hole.begin();
        cit != hole.end(); ++cit) {
      int indv = (*cit)->index(v);
      Cell_handle opp_cit = (*cit)->neighbor( indv );
      boundhole.push_back(Facet( opp_cit, opp_cit->index(*cit)) );

      for (int i=0; i<4; i++)
        if ( i != indv )
	  (*cit)->vertex(i)->set_cell(opp_cit);
    }
  }
#endif

  void fill_hole_3D_ear(const std::vector<Facet> & boundhole);

  class Conflict_tester_3
  {
      const Point &p;
      const Self *t;

  public:

    Conflict_tester_3(const Point &pt, const Self *tr)
      : p(pt), t(tr) {}

    bool operator()(const Cell_handle c) const
    {
      return t->side_of_sphere(c, p, true) == ON_BOUNDED_SIDE;
    }
    Oriented_side compare_weight(const Point &, const Point &) const
    {
      return ZERO;
    }
    bool test_initial_cell(Cell_handle) const
    {
      return true;
    }
  };

  class Conflict_tester_2
  {
      const Point &p;
      const Self *t;

  public:

    Conflict_tester_2(const Point &pt, const Self *tr)
      : p(pt), t(tr) {}

    bool operator()(const Cell_handle c) const
    {
      return t->side_of_circle(c, 3, p, true) == ON_BOUNDED_SIDE;
    }
    Oriented_side compare_weight(const Point &, const Point &) const
    {
      return ZERO;
    }
    bool test_initial_cell(Cell_handle) const
    {
      return true;
    }
  };
  class Hidden_point_visitor
  {
  public:

    Hidden_point_visitor() {}

    template <class InputIterator>
    void process_cells_in_conflict(InputIterator, InputIterator) const {}
    void reinsert_vertices(Vertex_handle ) {}
    Vertex_handle replace_vertex(Cell_handle c, int index,
				 const Point &) {
      return c->vertex(index);
    }
    void hide_point(Cell_handle, const Point &) {}
  };

  class Perturbation_order {
      const Self *t;

  public:
      Perturbation_order(const Self *tr)
	  : t(tr) {}

      bool operator()(const Point *p, const Point *q) const {
	  return t->compare_xyz(*p, *q) == SMALLER;
      }
  };

  template < class DelaunayTriangulation_3 >
  class Vertex_remover;

  friend class Perturbation_order;
  friend class Conflict_tester_3;
  friend class Conflict_tester_2;

  Hidden_point_visitor hidden_point_visitor;
};

template < class Gt, class Tds >
typename Delaunay_triangulation_3<Gt,Tds>::Vertex_handle
Delaunay_triangulation_3<Gt,Tds>::
insert(const Point & p, Cell_handle start)
{
    Locate_type lt;
    int li, lj;
    Cell_handle c = locate(p, lt, li, lj, start);
    return insert(p, lt, c, li, lj);
}

template < class Gt, class Tds >
typename Delaunay_triangulation_3<Gt,Tds>::Vertex_handle
Delaunay_triangulation_3<Gt,Tds>::
insert(const Point & p, Locate_type lt, Cell_handle c, int li, int lj)
{
  switch (dimension()) {
  case 3:
    {
      Conflict_tester_3 tester(p, this);
      Vertex_handle v = insert_in_conflict(p, lt, c, li, lj,
					   tester, hidden_point_visitor);
      return v;
    }// dim 3
  case 2:
    {
      Conflict_tester_2 tester(p, this);
      return insert_in_conflict(p, lt, c, li, lj,
				tester, hidden_point_visitor);
    }//dim 2
  default :
    // dimension <= 1
    // Do not use the generic insert.
    return Tr_Base::insert(p, c);
  }
}

template < class Gt, class Tds >
typename Delaunay_triangulation_3<Gt,Tds>::Vertex_handle
Delaunay_triangulation_3<Gt,Tds>::
move_point(Vertex_handle v, const Point & p)
{
    CGAL_triangulation_precondition(! is_infinite(v));
    CGAL_triangulation_expensive_precondition(is_vertex(v));

    // Dummy implementation for a start.

    // Remember an incident vertex to restart
    // the point location after the removal.
    Cell_handle c = v->cell();
    Vertex_handle old_neighbor = c->vertex(c->index(v) == 0 ? 1 : 0);
    CGAL_triangulation_assertion(old_neighbor != v);

    remove(v);

    if (dimension() <= 0)
	return insert(p);
    return insert(p, old_neighbor->cell());
}

template <class Gt, class Tds >
template <class DelaunayTriangulation_3>
class Delaunay_triangulation_3<Gt, Tds>::Vertex_remover {
  typedef DelaunayTriangulation_3 Delaunay;
public:
  typedef Nullptr_t Hidden_points_iterator;

  Vertex_remover(Delaunay &tmp_) : tmp(tmp_) {}

  Delaunay &tmp;

  void add_hidden_points(Cell_handle) {}
  Hidden_points_iterator hidden_points_begin() { return NULL; }
  Hidden_points_iterator hidden_points_end() { return NULL; }

  Bounded_side side_of_bounded_circle(const Point &p, const Point &q,
    const Point &r, const Point &s, bool perturb = false) const {
    return tmp.coplanar_side_of_bounded_circle(p,q,r,s,perturb);
  }
};

template < class Gt, class Tds >
void
Delaunay_triangulation_3<Gt,Tds>::
remove_3D_ear(Vertex_handle v)
{
  std::vector<Facet> boundhole; // facets on the boundary of the hole
  boundhole.reserve(64);        // 27 on average.
  std::vector<Cell_handle> hole;
  hole.reserve(64);

  make_hole_3D_ear(v, boundhole, hole);

  fill_hole_3D_ear(boundhole);
  tds().delete_vertex(v);
  tds().delete_cells(hole.begin(), hole.end());
}

template < class Gt, class Tds >
bool
Delaunay_triangulation_3<Gt,Tds>::
remove(Vertex_handle v)
{
#ifdef CGAL_DELAUNAY_3_OLD_REMOVE
  if (dimension()==3 && !test_dim_down(v)) {
    remove_3D_ear(v);
    return true;
  }
#endif
  Self tmp;
  Vertex_remover<Self> remover (tmp);
  Tr_Base::remove(v,remover);

  CGAL_triangulation_expensive_postcondition(is_valid());
  return true;
}

template < class Gt, class Tds >
Oriented_side
Delaunay_triangulation_3<Gt,Tds>::
side_of_oriented_sphere(const Point &p0, const Point &p1, const Point &p2,
	                const Point &p3, const Point &p, bool perturb) const
{
    CGAL_triangulation_precondition( orientation(p0, p1, p2, p3) == POSITIVE );

    Oriented_side os =
	geom_traits().side_of_oriented_sphere_3_object()(p0, p1, p2, p3, p);

    if (os != ON_ORIENTED_BOUNDARY || !perturb)
	return os;

    // We are now in a degenerate case => we do a symbolic perturbation.

    // We sort the points lexicographically.
    const Point * points[5] = {&p0, &p1, &p2, &p3, &p};
    std::sort(points, points+5, Perturbation_order(this) );

    // We successively look whether the leading monomial, then 2nd monomial
    // of the determinant has non null coefficient.
    // 2 iterations are enough (cf paper)
    for (int i=4; i>2; --i) {
        if (points[i] == &p)
            return ON_NEGATIVE_SIDE; // since p0 p1 p2 p3 are non coplanar
	                             // and positively oriented
        Orientation o;
        if (points[i] == &p3 && (o = orientation(p0,p1,p2,p)) != COPLANAR )
            return o;
        if (points[i] == &p2 && (o = orientation(p0,p1,p,p3)) != COPLANAR )
            return o;
        if (points[i] == &p1 && (o = orientation(p0,p,p2,p3)) != COPLANAR )
            return o;
        if (points[i] == &p0 && (o = orientation(p,p1,p2,p3)) != COPLANAR )
            return o;
    }

    CGAL_triangulation_assertion(false);
    return ON_NEGATIVE_SIDE;
}

template < class Gt, class Tds >
Bounded_side
Delaunay_triangulation_3<Gt,Tds>::
coplanar_side_of_bounded_circle(const Point &p0, const Point &p1,
	       const Point &p2, const Point &p, bool perturb) const
{
    // In dim==2, we should even be able to assert orient == POSITIVE.
    CGAL_triangulation_precondition( coplanar_orientation(p0, p1, p2)
	                             != COLLINEAR );

    Bounded_side bs =
      geom_traits().coplanar_side_of_bounded_circle_3_object()(p0, p1, p2, p);

    if (bs != ON_BOUNDARY || !perturb)
	return bs;

    // We are now in a degenerate case => we do a symbolic perturbation.

    // We sort the points lexicographically.
    const Point * points[4] = {&p0, &p1, &p2, &p};
    std::sort(points, points+4, Perturbation_order(this) );

    Orientation local = coplanar_orientation(p0, p1, p2);

    // we successively look whether the leading monomial, then 2nd monimial,
    // then 3rd monomial, of the determinant which has non null coefficient
    // [syl] : TODO : Probably it can be stopped earlier like the 3D version
    for (int i=3; i>0; --i) {
        if (points[i] == &p)
            return Bounded_side(NEGATIVE); // since p0 p1 p2 are non collinear
	                           // but not necessarily positively oriented
        Orientation o;
        if (points[i] == &p2
		&& (o = coplanar_orientation(p0,p1,p)) != COLLINEAR )
	    // [syl] : TODO : I'm not sure of the signs here (nor the rest :)
            return Bounded_side(o*local);
        if (points[i] == &p1
		&& (o = coplanar_orientation(p0,p,p2)) != COLLINEAR )
            return Bounded_side(o*local);
        if (points[i] == &p0
		&& (o = coplanar_orientation(p,p1,p2)) != COLLINEAR )
            return Bounded_side(o*local);
    }

    // case when the first non null coefficient is the coefficient of
    // the 4th monomial
    // moreover, the tests (points[] == &p) were false up to here, so the
    // monomial corresponding to p is the only monomial with non-zero
    // coefficient, it is equal to coplanar_orient(p0,p1,p2) == positive
    // so, no further test is required
    return Bounded_side(-local); //ON_UNBOUNDED_SIDE;
}

template < class Gt, class Tds >
Bounded_side
Delaunay_triangulation_3<Gt,Tds>::
side_of_sphere(Vertex_handle v0, Vertex_handle v1,
	       Vertex_handle v2, Vertex_handle v3,
	       const Point &p, bool perturb) const
{
    CGAL_triangulation_precondition( dimension() == 3 );

    if (is_infinite(v0)) {
	Orientation o = orientation(v2->point(), v1->point(), v3->point(), p);
	if (o != COPLANAR)
	    return Bounded_side(o);
	return coplanar_side_of_bounded_circle(v2->point(), v1->point(), v3->point(), p, perturb);
    }

    if (is_infinite(v1)) {
	Orientation o = orientation(v2->point(), v3->point(), v0->point(), p);
	if (o != COPLANAR)
	    return Bounded_side(o);
	return coplanar_side_of_bounded_circle(v2->point(), v3->point(), v0->point(), p, perturb);
    }

    if (is_infinite(v2)) {
	Orientation o = orientation(v1->point(), v0->point(), v3->point(), p);
	if (o != COPLANAR)
	    return Bounded_side(o);
	return coplanar_side_of_bounded_circle(v1->point(), v0->point(), v3->point(), p, perturb);
    }

    if (is_infinite(v3)) {
	Orientation o = orientation(v0->point(), v1->point(), v2->point(), p);
	if (o != COPLANAR)
	    return Bounded_side(o);
	return coplanar_side_of_bounded_circle(v0->point(), v1->point(), v2->point(), p, perturb);
    }

    return (Bounded_side) side_of_oriented_sphere(v0->point(), v1->point(), v2->point(), v3->point(), p, perturb);
}

template < class Gt, class Tds >
Bounded_side
Delaunay_triangulation_3<Gt,Tds>::
side_of_circle(Cell_handle c, int i,
	       const Point & p, bool perturb) const
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
					      p, perturb);
    // else infinite facet
    // v1, v2 finite vertices of the facet such that v1,v2,infinite
    // is positively oriented
    Vertex_handle v1 = c->vertex( ccw(i3) ),
                  v2 = c->vertex( cw(i3) );
    CGAL_triangulation_assertion(coplanar_orientation(v1->point(), v2->point(),
			         mirror_vertex(c, i3)->point()) == NEGATIVE);
    Orientation o = coplanar_orientation(v1->point(), v2->point(), p);
    if ( o != COLLINEAR )
	return Bounded_side( o );
    // because p is in f iff
    // it does not lie on the same side of v1v2 as vn
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
    CGAL_triangulation_precondition( coplanar( c->vertex(i0)->point(),
				               c->vertex(i1)->point(),
				               c->vertex(i2)->point(),
					       p ) );
    return coplanar_side_of_bounded_circle( c->vertex(i0)->point(),
					    c->vertex(i1)->point(),
					    c->vertex(i2)->point(),
					    p, perturb);
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
  if ( o != COLLINEAR )
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
typename Delaunay_triangulation_3<Gt,Tds>::Vertex_handle
Delaunay_triangulation_3<Gt,Tds>::
nearest_vertex_in_cell(const Point& p, Cell_handle c) const
// Returns the finite vertex of the cell c which is the closest to p.
{
    CGAL_triangulation_precondition(dimension() >= 1);

    Vertex_handle nearest = nearest_vertex(p, c->vertex(0), c->vertex(1));
    if (dimension() >= 2) {
	nearest = nearest_vertex(p, nearest, c->vertex(2));
        if (dimension() == 3)
	    nearest = nearest_vertex(p, nearest, c->vertex(3));
    }
    return nearest;
}

template < class Gt, class Tds >
typename Delaunay_triangulation_3<Gt,Tds>::Vertex_handle
Delaunay_triangulation_3<Gt,Tds>::
nearest_vertex(const Point& p, Cell_handle start) const
{
    if (number_of_vertices() == 0)
	return Vertex_handle();

    // Use a brute-force algorithm if dimension < 3.
    if (dimension() < 3) {
	Finite_vertices_iterator vit = finite_vertices_begin();
	Vertex_handle res = vit;
	for (++vit; vit != finite_vertices_end(); ++vit)
	    res = nearest_vertex(p, res, vit);
	return res;
    }

    Locate_type lt;
    int li, lj;
    Cell_handle c = locate(p, lt, li, lj, start);
    if (lt == Tr_Base::VERTEX)
	return c->vertex(li);

    // - start with the closest vertex from the located cell.
    // - repeatedly take the nearest of its incident vertices if any
    // - if not, we're done.
    Vertex_handle nearest = nearest_vertex_in_cell(p, c);
    std::vector<Vertex_handle> vs;
    vs.reserve(32);
    while (true) {
	Vertex_handle tmp = nearest;
        adjacent_vertices(nearest, std::back_inserter(vs));
        for (typename std::vector<Vertex_handle>::const_iterator
		vsit = vs.begin(); vsit != vs.end(); ++vsit)
	    tmp = nearest_vertex(p, tmp, *vsit);
	if (tmp == nearest)
	    break;
	vs.clear();
	nearest = tmp;
    }

    return nearest;
}


template < class Gt, class Tds >
bool
Delaunay_triangulation_3<Gt,Tds>::
is_Gabriel(const Facet& f) const
{
  return is_Gabriel(f.first, f.second);
}

template < class Gt, class Tds >
bool
Delaunay_triangulation_3<Gt,Tds>::
is_Gabriel(Cell_handle c, int i) const
{
  CGAL_triangulation_precondition(dimension() == 3 && !is_infinite(c,i));
  typename Geom_traits::Side_of_bounded_sphere_3
    side_of_bounded_sphere =
    geom_traits().side_of_bounded_sphere_3_object();

  if ((!is_infinite(c->vertex(i))) &&
      side_of_bounded_sphere (
	c->vertex(vertex_triple_index(i,0))->point(),
	c->vertex(vertex_triple_index(i,1))->point(),
	c->vertex(vertex_triple_index(i,2))->point(),
	c->vertex(i)->point()) == ON_BOUNDED_SIDE ) return false;
    Cell_handle neighbor = c->neighbor(i);
  int in = neighbor->index(c);

  if ((!is_infinite(neighbor->vertex(in))) &&
      side_of_bounded_sphere(
	 c->vertex(vertex_triple_index(i,0))->point(),
	 c->vertex(vertex_triple_index(i,1))->point(),
	 c->vertex(vertex_triple_index(i,2))->point(),	
	 neighbor->vertex(in)->point()) == ON_BOUNDED_SIDE ) return false;

  return true;
}

template < class Gt, class Tds >
bool
Delaunay_triangulation_3<Gt,Tds>::
is_Gabriel(const Edge& e) const
{
  return is_Gabriel(e.first, e.second, e.third);
}

template < class Gt, class Tds >
bool
Delaunay_triangulation_3<Gt,Tds>::
is_Gabriel(Cell_handle c, int i, int j) const
{
  CGAL_triangulation_precondition(dimension() == 3 && !is_infinite(c,i,j));
  typename Geom_traits::Side_of_bounded_sphere_3
    side_of_bounded_sphere =
    geom_traits().side_of_bounded_sphere_3_object();

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
	   side_of_bounded_sphere( v1->point(),
				   v2->point(),
				   cc->vertex(ii)->point())
	  == ON_BOUNDED_SIDE ) return false;
  } while(++fcirc != fdone);
  return true;
}

template < class Gt, class Tds >
typename Delaunay_triangulation_3<Gt,Tds>::Point
Delaunay_triangulation_3<Gt,Tds>::
dual(Cell_handle c) const
{
  CGAL_triangulation_precondition(dimension()==3);
  CGAL_triangulation_precondition( ! is_infinite(c) );
  return c->circumcenter(geom_traits());
}


template < class Gt, class Tds >
typename Delaunay_triangulation_3<Gt,Tds>::Object
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

  Line l = construct_equidistant_line( p, q, r );
  return construct_object(construct_ray( dual(n), l));
}



template < class Gt, class Tds >
typename Delaunay_triangulation_3<Gt,Tds>::Line
Delaunay_triangulation_3<Gt,Tds>::
dual_support(Cell_handle c, int i) const
{
  CGAL_triangulation_precondition(dimension()>=2);
  CGAL_triangulation_precondition( ! is_infinite(c,i) );

  if ( dimension() == 2 ) {
    CGAL_triangulation_precondition( i == 3 );
    return construct_equidistant_line( c->vertex(0)->point(),
		                       c->vertex(1)->point(),
			               c->vertex(2)->point() );
  }

  return construct_equidistant_line( c->vertex((i+1)&3)->point(),
		                     c->vertex((i+2)&3)->point(),
				     c->vertex((i+3)&3)->point() );
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

  if ( infinite_vertex() == Vertex_handle() ) {
    if (verbose)
	std::cerr << "no infinite vertex" << std::endl;
    CGAL_triangulation_assertion(false);
    return false;
  }

  switch ( dimension() ) {
  case 3:
    {
      Finite_cells_iterator it;
      for ( it = finite_cells_begin(); it != finite_cells_end(); ++it ) {
	is_valid_finite(it);
	for (int i=0; i<4; i++ ) {
	  if ( !is_infinite
	       (it->neighbor(i)->vertex(it->neighbor(i)->index(it))) ) {
	    if ( side_of_sphere
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
	is_valid_finite((*it).first);
	for (int i=0; i<3; i++ ) {
	  if( !is_infinite
	      ((*it).first->neighbor(i)->vertex( (((*it).first)->neighbor(i))
						 ->index((*it).first))) ) {
	    if ( side_of_circle ( (*it).first, 3,
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
      for ( it = finite_edges_begin(); it != finite_edges_end(); ++it )
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
  if ( ! Tr_Base::is_valid(c,verbose,level) ) {
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

#ifndef CGAL_CFG_NET2003_MATCHING_BUG
template < class Gt, class Tds >
void
Delaunay_triangulation_3<Gt,Tds>::
make_hole_3D_ear( Vertex_handle v,
	          std::vector<Facet> & boundhole,
	          std::vector<Cell_handle> & hole)
{
  CGAL_triangulation_expensive_precondition( ! test_dim_down(v) );

  incident_cells(v, std::back_inserter(hole));

  for (typename std::vector<Cell_handle>::iterator cit = hole.begin();
       cit != hole.end(); ++cit) {
    int indv = (*cit)->index(v);
    Cell_handle opp_cit = (*cit)->neighbor( indv );
    boundhole.push_back(Facet( opp_cit, opp_cit->index(*cit)) );

    for (int i=0; i<4; i++)
      if ( i != indv )
	(*cit)->vertex(i)->set_cell(opp_cit);
  }
}
#endif





template < class Gt, class Tds >
void
Delaunay_triangulation_3<Gt,Tds>::
fill_hole_3D_ear(const std::vector<Facet> & boundhole)
{
  typedef internal::Delaunay_remove_tds_3_2<Delaunay_triangulation_3> Surface;
  typedef typename Surface::Face_3_2          Face_3_2;
  typedef typename Surface::Face_handle_3_2   Face_handle_3_2;
  typedef typename Surface::Vertex_handle_3_2 Vertex_handle_3_2;

  Surface surface(boundhole);

  Face_handle_3_2 f = surface.faces_begin();
  Face_handle_3_2 last_op = f; // This is where the last ear was inserted

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
      f = f->n();
      CGAL_assertion_msg(f != last_op, "Unable to find an ear");
      k = 0;
    }

    // The edges are marked, if they are a candidate for an ear.
    // This saves time, for example an edge gets not considered
    // from both adjacent faces.
    if (!f->is_halfedge_marked(k))
	continue;

    Vertex_handle_3_2 w0, w1, w2, w3;
    Vertex_handle v0, v1, v2, v3;
    int i = ccw(k);
    int j = cw(k);
    Face_handle_3_2 n = f->neighbor(k);
    int fi = n->index(f);

    w1 = f->vertex(i);
    w2 = f->vertex(j);

    v1 = w1->info();
    v2 = w2->info();

    if( is_infinite(v1) || is_infinite(v2) ){
	// there will be another ear, so let's ignore this one,
	// because it is complicated to treat
	continue;
    }
    w0 = f->vertex(k);
    w3 = n->vertex(fi);

    v0 = w0->info();
    v3 = w3->info();

    if( !is_infinite(v0) && !is_infinite(v3) &&
	  orientation(v0->point(), v1->point(),
		      v2->point(), v3->point()) != POSITIVE)
        continue;

    // the two faces form a concavity, in which we might plug a cell

    // we now look at all vertices that are on the boundary of the hole
    for(typename Surface::Vertex_iterator vit = surface.vertices_begin();
	vit != surface.vertices_end(); ++vit) {
      Vertex_handle v = vit->info();
      if (is_infinite(v) || v == v0 || v == v1 || v == v2 || v == v3)
	  continue;

      if (side_of_sphere(v0,v1,v2,v3, v->point(), true) == ON_BOUNDED_SIDE)
	  goto next_edge;
    }

    // we looked at all vertices

    Face_handle_3_2 m_i = f->neighbor(i);
    Face_handle_3_2 m_j = f->neighbor(j);
    bool neighbor_i = m_i == n->neighbor(cw(fi));
    bool neighbor_j = m_j == n->neighbor(ccw(fi));

    // Test if the edge that would get introduced is on the surface
    if ( !neighbor_i && !neighbor_j &&
	 surface.is_edge(f->vertex(k), n->vertex(fi)))
      continue;

    // none of the vertices violates the Delaunay property
    // We are ready to plug a new cell

    Cell_handle ch = tds().create_cell(v0, v1, v2, v3);

    // The new cell touches the faces that form the ear
    Facet fac = n->info();
    tds().set_adjacency(ch, 0, fac.first, fac.second);
    fac = f->info();
    tds().set_adjacency(ch, 3, fac.first, fac.second);

    // It may touch another face,
    // or even two other faces if it is the last cell
    if(neighbor_i) {
      fac = m_i->info();
      tds().set_adjacency(ch, 1, fac.first, fac.second);
    }
    if(neighbor_j) {
      fac = m_j->info();
      tds().set_adjacency(ch, 2, fac.first, fac.second);
    }

    if( !neighbor_i && !neighbor_j) {
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

      f->set_info(Facet(ch,2));
      n->set_info(Facet(ch,1));
    } else if (neighbor_i && (! neighbor_j)) {
      surface.remove_degree_3(f->vertex(j), f);
      // all three edges adjacent to f are
      // candidate for an ear
      f->mark_adjacent_edges();
      f->set_info(Facet(ch,2));
    } else if ((! neighbor_i) && neighbor_j)  {
      surface.remove_degree_3(f->vertex(i), f);
      f->mark_adjacent_edges();
      f->set_info(Facet(ch,1));
    } else {
      CGAL_assertion(surface.number_of_vertices() == 4);
      // when we leave the function the vertices and faces of the surface
      // are deleted by the destructor
      return;
    }

    // we successfully inserted a cell
    last_op = f;
    // we have to reconsider all edges incident to f
    k = -1;
  } // for(;;)
}

CGAL_END_NAMESPACE

#include <CGAL/internal/Delaunay_triangulation_hierarchy_3.h>

#endif // CGAL_DELAUNAY_TRIANGULATION_3_H
