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
//                 Sylvain Pion <Sylvain.Pion@sophia.inria.fr>
//                 
// coordinator   : INRIA Sophia Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================

#ifndef CGAL_REGULAR_TRIANGULATION_3_H
#define CGAL_REGULAR_TRIANGULATION_3_H

#include <CGAL/basic.h>

#include <set>

#include <CGAL/Triangulation_short_names_3.h>
#include <CGAL/Triangulation_3.h>

CGAL_BEGIN_NAMESPACE

template < class Gt, 
           class Tds = Triangulation_data_structure_3 <
                                   Triangulation_vertex_base_3<Gt>,
                                   Triangulation_cell_base_3<void> > >
class Regular_triangulation_3
  : public Triangulation_3<Gt,Tds>
{
  friend std::istream& operator >> CGAL_NULL_TMPL_ARGS
  (std::istream& is, Triangulation_3<Gt,Tds> &tr);

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

  typedef typename Tr_Base::Locate_type         Locate_type;
  typedef typename Tr_Base::Cell_iterator       Cell_iterator;
  typedef typename Tr_Base::Facet_iterator      Facet_iterator;
  typedef typename Tr_Base::Edge_iterator       Edge_iterator;

  typedef typename Tr_Base::Finite_vertices_iterator Finite_vertices_iterator;
  typedef typename Tr_Base::Finite_cells_iterator   Finite_cells_iterator;
  typedef typename Tr_Base::Finite_facets_iterator  Finite_facets_iterator;
  typedef typename Tr_Base::Finite_edges_iterator   Finite_edges_iterator;

  typedef typename Gt::Weighted_point              Weighted_point;

  Regular_triangulation_3(const Gt & gt = Gt())
    : Tr_Base(gt)
  {}

  // copy constructor duplicates vertices and cells
  Regular_triangulation_3(const Regular_triangulation_3 & rt)
      : Tr_Base(rt)
  { 
      CGAL_triangulation_postcondition( is_valid() );  
  }

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

  Vertex_handle insert(const Weighted_point & p, Cell_handle start = NULL);

  Vertex_handle insert(const Weighted_point & p, Locate_type lt,
	               Cell_handle c, int li, int);

  Vertex_handle push_back(const Weighted_point &p)
  {
      return insert(p);
  }

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
		  if (v->cell() != NULL)
		  {
		      cv.push_back(v);
		      v->set_cell(NULL);
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
		  if (v->cell() != NULL)
		  {
		      cv.push_back(v);
		      v->set_cell(NULL);
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
          if (lt == this->VERTEX)
              return c->vertex(li); // by coinciding point
          else
              return NULL;          // by cell
      }

      // Should I mark c's vertices too ?
      Conflict_tester_3 tester(p, this);
      Vertex_handle v = insert_conflict_3(c, tester);
      v->set_point(p);
      for( typename std::vector<Vertex_handle>::iterator
		it = tester.conflict_vector().begin();
		it != tester.conflict_vector().end(); ++it)
      {
        if ((*it)->cell() == NULL)
	{
          // vertex has to be deleted
          tds().delete_vertex(*it);
	}
      }
      // TODO : manage the hidden points.
      return v;
    }
  case 2:
    {
      switch (lt) {
      case this->OUTSIDE_CONVEX_HULL:
      case this->CELL:
      case this->FACET:
      case this->EDGE:
      case this->VERTEX:
	{
          if (! in_conflict_2(p, c, 3)) {  // new point is hidden
              if (lt == this->VERTEX)
                  return c->vertex(li); // by coinciding point
              else
                  return NULL;          // by face
          }

	  Conflict_tester_2 tester(p, this);
	  Vertex_handle v = insert_conflict_2(c, tester);
	  v->set_point(p);
          for( typename std::vector<Vertex_handle>::iterator
		it = tester.conflict_vector().begin();
		it != tester.conflict_vector().end(); ++it)
	  {
            if ((*it)->cell() == NULL)
	    {
              // vertex has to be deleted
              tds().delete_vertex(*it);
	    }
	  }
	  return v;
	}
      case this->OUTSIDE_AFFINE_HULL:
	{
	  // if the 2d triangulation is Regular, the 3d
	  // triangulation will be Regular
	  return Tr_Base::insert_outside_affine_hull(p);
	}
      }
    }//dim 2
  case 1:
    {
      switch (lt) {
      case this->OUTSIDE_CONVEX_HULL:
      case this->EDGE:
      case this->VERTEX:
	{
          if (! in_conflict_1(p, c)) {  // new point is hidden
              if (lt == this->VERTEX)
                  return c->vertex(li); // by coinciding point
              else
                  return NULL;          // by edge
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
          Cell_handle c0 = tds().create_face(v, bound[0]->vertex(0), NULL);
          Cell_handle c1 = tds().create_face(bound[1]->vertex(1), v, NULL);
          tds().set_adjacency(c0, 1, c1, 0);
          tds().set_adjacency(bound[0], 1, c0, 0);
          tds().set_adjacency(c1, 1, bound[1], 0);
          bound[0]->vertex(0)->set_cell(bound[0]);
          bound[1]->vertex(1)->set_cell(bound[1]);
          v->set_cell(c0);

	  tds().delete_cells(conflicts.begin(), conflicts.end());
	  tds().delete_vertices(hidden_vertices.begin(), hidden_vertices.end());
	  return v;
	}
      case this->OUTSIDE_AFFINE_HULL:
	return Tr_Base::insert_outside_affine_hull(p);
      case this->FACET:
      case this->CELL:
	// impossible in dimension 1
        CGAL_assertion(false);
	return NULL;
      }
    }
  case 0:
    {
        // We need to compare the weights when the points are equal.
        if (lt == this->VERTEX && in_conflict_0(p, c)) {
            CGAL_assertion(li == 0);
            c->vertex(li)->set_point(p); // replace by heavier point
        }
        else 
            return Tr_Base::insert(p, c);
    }
  default :
    {
      return Tr_Base::insert(p, c);
    }
  }
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
	  if ( side_of_power_sphere (it, 
		 it->vertex(it->neighbor(i)->index(it))->point() )
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
      Finite_facets_iterator it;
      for ( it = finite_facets_begin(); it != finite_facets_end(); ++it ) {
	is_valid_finite((*it).first, verbose, level);
	for (int i=0; i<3; i++ ) {
	  if ( side_of_power_circle
	       ( (*it).first, 3,
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
      Finite_edges_iterator it;
      for ( it = finite_edges_begin(); it != finite_edges_end(); ++it ) {
	is_valid_finite((*it).first, verbose, level);
	for (int i=0; i<2; i++ ) {
	  if ( side_of_power_segment
	       ( (*it).first,
		 (*it).first->vertex( (((*it).first)->neighbor(i))
				      ->index((*it).first) )->point() )
	       == ON_BOUNDED_SIDE ) {
	    if (verbose)
		std::cerr << "non-empty edge " << std::endl;
	    CGAL_triangulation_assertion(false);
	    return false;
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
