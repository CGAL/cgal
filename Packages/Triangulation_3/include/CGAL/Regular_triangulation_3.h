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

  Regular_triangulation_3()
    : Tr_Base()
  {}
  
  Regular_triangulation_3(const Gt & gt)
    : Tr_Base(gt)
  {}
  
  // copy constructor duplicates vertices and cells
  Regular_triangulation_3(const Regular_triangulation_3 & rt)
      : Tr_Base(rt)
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
  power_test(const Weighted_point &p, const Weighted_point &q,
	     const Weighted_point &r) const
  {
      return geom_traits().power_test_3_object()(p, q, r);
  }

  Oriented_side
  power_test(const Weighted_point &p, const Weighted_point &q,
	     const Weighted_point &r, const Weighted_point &s) const
  {
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
  return side_of_edge( p, c, lt, i );
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
      if ( lt == VERTEX )
	  return c->vertex(li);
      // choice: not to do anything
      // could reinsert the point with the new weight
      
      if (! in_conflict_3(p, c))
	  return NULL;
      // Should I mark c's vertices too ?
      Conflict_tester_3 tester(p, this);
      Vertex_handle v = insert_conflict(c, tester);
      v->set_point(p);
      for( typename std::vector<Vertex_handle>::iterator
		it = tester.conflict_vector().begin();
		it != tester.conflict_vector().end(); ++it)
      {
        if ((*it)->cell() == NULL)
	{
          // vertex has to be deleted
          _tds.delete_vertex(*it);
	}
      }
      // a voir : comment compter les sommets redondants ? (***)
      // else : traiter le cas des points redondants a stocker dans
      // la face associee pour le cas d'une future suppression
      return v;
    }
  case 2:
    {
      switch (lt) {
      case OUTSIDE_CONVEX_HULL:
      case CELL:
      case FACET:
      case EDGE:
	{
	  if (! in_conflict_2(p, c, 3))
	      return NULL;
	  Conflict_tester_2 tester(p, this);
	  Vertex_handle v = insert_conflict(c, tester);
	  v->set_point(p);
          for( typename std::vector<Vertex_handle>::iterator
		it = tester.conflict_vector().begin();
		it != tester.conflict_vector().end(); ++it)
	  {
            if ((*it)->cell() == NULL)
	    {
              // vertex has to be deleted
              _tds.delete_vertex(*it);
	    }
	  }
	  return v;
	}
      case VERTEX:
	return c->vertex(li);
	// choice: not to do anything
	// could reinsert the point with the new weight
      case OUTSIDE_AFFINE_HULL:
	{
	  // if the 2d triangulation is Delaunay, the 3d
	  // triangulation will be Delaunay
	  return Tr_Base::insert_outside_affine_hull(p);
	}
      }
    }//dim 2
  case 1:
    {
      switch (lt) {
      case OUTSIDE_CONVEX_HULL:
      case EDGE:
	{
	  if (! in_conflict_1(p, c) )
	      return NULL;
	  Cell_handle n, bound[2];
	  std::set<Cell_handle> conflicts;

	  for (int j =0; j<2; j++ ) {
	    n = c;
	    while ( ( ! is_infinite(n->vertex(1-j)) ) && 
		      in_conflict_1( p, n->neighbor(j) ) ) {
	      if (n!=c)
		  (void) conflicts.insert(n);
// 		P = new Weighted_point( n->vertex(1-j)->point() );
// 		(void) deleted_points.insert((void*) P);
	      _tds.delete_vertex( n->vertex(1-j) );
	      n = n->neighbor(j);
	    }
	    bound[j] = n;
	  }

	  Vertex_handle v = _tds.create_vertex();
	  v->set_point(p);
	  if ( bound[0] != bound[1] ) {
	    if ( (c != bound[0]) && (c != bound[1]) )
	      (void) conflicts.insert(c);
	    bound[1]->set_vertex(1,v);
	  }
	  else {
	    bound[1] = _tds.create_face(bound[0]->vertex(0), v, NULL);
	    _tds.set_adjacency(bound[1], 1, bound[0]->neighbor(1), 0);
	    bound[0]->vertex(0)->set_cell(bound[1]);
	  }
	  bound[0]->set_vertex(0,v);
	  _tds.set_adjacency(bound[1], 0, bound[0], 1);
	  v->set_cell(bound[0]);

	  _tds.delete_cells(conflicts.begin(), conflicts.end());
	  return v;
	}
      case VERTEX:
	return c->vertex(li);
	// choice: not to do anything
	// could reinsert the point with the new weight
      case OUTSIDE_AFFINE_HULL:
	return Tr_Base::insert_outside_affine_hull(p);
      case FACET:
      case CELL:
	// impossible in dimension 1
	return NULL;
      }
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
  if ( ! tds().is_valid(verbose,level) ) {
    if (verbose)
	std::cerr << "invalid data structure" << std::endl;
    CGAL_triangulation_assertion(false);
    return false;
  }

  if ( infinite_vertex() == NULL ) {
    if (verbose)
	std::cerr << "no infinite vertex" << std::endl;
    CGAL_triangulation_assertion(false);
    return false;
  }

  int i;
  switch ( dimension() ) {
  case 3:
    {
      Finite_cells_iterator it;
      for ( it = finite_cells_begin(); it != finite_cells_end(); ++it ) {
	is_valid_finite(it);
	for ( i=0; i<4; i++ ) {
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
	is_valid_finite((*it).first);
	for ( i=0; i<3; i++ ) {
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
	is_valid_finite((*it).first);
	for ( i=0; i<2; i++ ) {
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
      std::cerr << "Regular valid triangulation" << std::endl;
  return true;
}

CGAL_END_NAMESPACE

#endif // CGAL_REGULAR_TRIANGULATION_3_H
