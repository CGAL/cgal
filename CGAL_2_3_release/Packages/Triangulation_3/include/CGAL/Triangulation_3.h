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
// file          : include/CGAL/Triangulation_3.h
// revision      : $Revision$
//
// author(s)     : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//                 Sylvain Pion <Sylvain.Pion@sophia.inria.fr>
//
// coordinator   : INRIA Sophia Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================


#ifndef CGAL_TRIANGULATION_3_H
#define CGAL_TRIANGULATION_3_H

#include <CGAL/basic.h>

#include <iostream>
#include <set>
#include <map> 
#include <utility>

#include <CGAL/Triangulation_short_names_3.h>

#include <CGAL/Triangulation_utils_3.h>
#include <CGAL/triple.h>
#include <CGAL/Pointer.h>
#include <CGAL/circulator.h>

#include <CGAL/triangulation_assertions.h>

#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Triangulation_cell_base_3.h>
#include <CGAL/Triangulation_vertex_base_3.h>

#include <CGAL/Triangulation_cell_3.h>
#include <CGAL/Triangulation_vertex_3.h>
#include <CGAL/Triangulation_handles_3.h>
#include <CGAL/Triangulation_iterators_3.h>
#include <CGAL/Triangulation_circulators_3.h>

CGAL_BEGIN_NAMESPACE

template < class GT, class Tds > std::istream& operator>> 
(std::istream& is, Triangulation_3<GT,Tds> &tr);
 
template < class GT, 
           class Tds = Triangulation_data_structure_3 <
                                   Triangulation_vertex_base_3<GT>,
                                   Triangulation_cell_base_3<void> > >
class Triangulation_3
  :public Triangulation_utils_3
{
  friend std::istream& operator>> CGAL_NULL_TMPL_ARGS
  (std::istream& is, Triangulation_3<GT,Tds> &tr);

  friend class Triangulation_cell_3<GT,Tds>;
  friend class Triangulation_vertex_3<GT,Tds>;

  friend class Triangulation_cell_iterator_3<GT,Tds>;
  friend class Triangulation_facet_iterator_3<GT,Tds>;
  friend class Triangulation_edge_iterator_3<GT,Tds>;
  friend class Triangulation_vertex_iterator_3<GT,Tds>;

  friend class Triangulation_cell_circulator_3<GT,Tds>;
  friend class Triangulation_facet_circulator_3<GT,Tds>;

public:
  typedef Tds                                  Triangulation_data_structure;
  typedef GT                                   Geom_traits;
  typedef Triangulation_3<GT, Tds>             Self;

  typedef typename GT::Point_3                 Point;
  typedef typename GT::Segment_3               Segment;
  typedef typename GT::Triangle_3              Triangle;
  typedef typename GT::Tetrahedron_3           Tetrahedron;

  typedef Triangulation_cell_handle_3<GT,Tds>          Cell_handle;
  typedef Triangulation_vertex_handle_3<GT,Tds>        Vertex_handle;

  typedef Triangulation_cell_3<GT,Tds>                 Cell;
  typedef Triangulation_vertex_3<GT,Tds>               Vertex;
  typedef std::pair<Cell_handle, int>                  Facet;
  typedef triple<Cell_handle, int, int>                Edge;

  typedef Triangulation_cell_circulator_3<GT,Tds>      Cell_circulator;
  typedef Triangulation_facet_circulator_3<GT,Tds>     Facet_circulator;

  typedef Triangulation_cell_iterator_3<GT,Tds>        Cell_iterator;
  typedef Triangulation_facet_iterator_3<GT,Tds>       Facet_iterator;
  typedef Triangulation_edge_iterator_3<GT,Tds>        Edge_iterator;
  typedef Triangulation_vertex_iterator_3<GT,Tds>      Vertex_iterator;

  typedef Point                         value_type; // to have a back_inserter
  typedef const value_type&             const_reference;

  enum Locate_type {
    VERTEX=0, 
    EDGE, //1
    FACET, //2
    CELL, //3
    OUTSIDE_CONVEX_HULL, //4
    OUTSIDE_AFFINE_HULL };//5

protected:
  Tds _tds;
  GT  _gt;
  Vertex_handle infinite; //infinite vertex
 
  Comparison_result
  compare_x(const Point &p, const Point &q) const
  {
      return geom_traits().compare_x_3_object()(p, q);
  }

  Comparison_result
  compare_y(const Point &p, const Point &q) const
  {
      return geom_traits().compare_y_3_object()(p, q);
  }

  Comparison_result
  compare_z(const Point &p, const Point &q) const
  {
      return geom_traits().compare_z_3_object()(p, q);
  }

  bool
  equal(const Point &p, const Point &q) const
  {
      return geom_traits().equal_3_object()(p, q);
  }

  Orientation
  orientation(const Point &p, const Point &q,
	      const Point &r, const Point &s) const
  {
      return geom_traits().orientation_3_object()(p, q, r, s);
  }

  Orientation
  coplanar_orientation(const Point &p, const Point &q, const Point &r) const
  {
      return geom_traits().coplanar_orientation_3_object()(p, q, r);
  }

  bool
  collinear(const Point &p, const Point &q, const Point &r) const
  {
      return coplanar_orientation(p, q, r) == COLLINEAR;
  }

  Segment
  construct_segment(const Point &p, const Point &q) const
  {
      return geom_traits().construct_segment_3_object()(p, q);
  }

  Triangle
  construct_triangle(const Point &p, const Point &q, const Point &r) const
  {
      return geom_traits().construct_triangle_3_object()(p, q, r);
  }

  Tetrahedron
  construct_tetrahedron(const Point &p, const Point &q,
	                const Point &r, const Point &s) const
  {
      return geom_traits().construct_tetrahedron_3_object()(p, q, r, s);
  }

  void init_tds()
    {
      infinite = (Vertex*) _tds.insert_increase_dimension(NULL);

      // Forces the compiler to instantiate handle2pointer.
      handle2pointer( Vertex_handle() ); 
      handle2pointer( Cell_handle() );
    }
  
  bool test_dim_down(Vertex_handle v);

public:

  // CONSTRUCTORS
  Triangulation_3()
    : _tds(), _gt()
    {
      init_tds();
    }

  Triangulation_3(const GT & gt) 
    : _tds(), _gt(gt)
    {
      init_tds();
    }

  // copy constructor duplicates vertices and cells
  Triangulation_3(const Triangulation_3<GT,Tds> & tr)
    : _gt(tr._gt)
    {
      infinite = (Vertex *) _tds.copy_tds(tr._tds, &(*(tr.infinite)) );
    }

  void clear()
    {
      _tds.clear();
      init_tds();
    }

  Triangulation_3 & operator=(const Triangulation_3 & tr)
    {
      //     clear();               BUG !!
      //     infinite.Delete();
      infinite = (Vertex *) _tds.copy_tds( tr._tds, &*tr.infinite );
      _gt = tr._gt;
      return *this;
    }

  // HELPING FUNCTIONS
   
  void copy_triangulation(const Triangulation_3<GT,Tds> & tr)
    {
      //     clear();               BUG !!
      //     infinite.Delete();
      _gt = tr._gt;
      infinite = (Vertex *) _tds.copy_tds( tr._tds, &*tr.infinite );
    }

  void swap(Triangulation_3 &tr)
    {
      GT t(geom_traits());
      _gt = tr.geom_traits();
      tr._gt = t; 

      Vertex_handle inf = infinite_vertex();
      infinite = tr.infinite_vertex();
      tr.infinite = inf;

      _tds.swap(tr._tds);
    }

  //ACCESS FUNCTIONS
  const GT & geom_traits() const 
    { return _gt;}
  
  const Tds & tds() const 
    { return _tds;}
  
  int dimension() const 
    { return _tds.dimension();}

  int number_of_finite_cells() const;

  int number_of_cells() const;
 
  int number_of_finite_facets() const;

  int number_of_facets() const;

  int number_of_finite_edges() const;
 
  int number_of_edges() const;
  
  int number_of_vertices() const // number of finite vertices
    {return _tds.number_of_vertices()-1;}

  Vertex_handle infinite_vertex() const
    { return infinite; }
   
  Cell_handle infinite_cell() const
    {
      //    CGAL_triangulation_precondition(infinite_vertex()->cell()->
      //				    has_vertex(infinite_vertex()));
      return infinite_vertex()->cell();
    }

  // ASSIGNMENT
  void set_number_of_vertices(int n) 
    { _tds.set_number_of_vertices(n+1); }
   
  // GEOMETRIC ACCESS FUNCTIONS
  
  Tetrahedron tetrahedron(const Cell_handle c) const
    {
      CGAL_triangulation_precondition( dimension() == 3 );
      CGAL_triangulation_precondition( ! is_infinite(c) );
      return construct_tetrahedron(c->vertex(0)->point(),
				   c->vertex(1)->point(),
				   c->vertex(2)->point(),
				   c->vertex(3)->point());
    }

  Triangle triangle(const Cell_handle c, int i) const;

  Triangle triangle(const Facet & f) const
    { return triangle(f.first, f.second); }

  Segment segment(const Cell_handle c, int i, int j) const;

  Segment segment(const Edge & e) const
    { return segment(e.first,e.second,e.third); }

  // TEST IF INFINITE FEATURES
  bool is_infinite(const Vertex_handle v) const 
    { return v == infinite_vertex(); }

  bool is_infinite(const Cell_handle c) const 
    {
      CGAL_triangulation_precondition( dimension() == 3 );
      return c->has_vertex(infinite_vertex());
    }
  
  bool is_infinite(const Cell_handle c, int i) const;

  bool is_infinite(const Facet & f) const 
    { return is_infinite(f.first,f.second); }

  bool is_infinite(const Cell_handle c, int i, int j) const; 

  bool is_infinite(const Edge & e) const
    { return is_infinite(e.first,e.second,e.third); }


  //QUERIES

  bool is_vertex(const Point & p, Vertex_handle & v) const;

  bool is_vertex(Vertex_handle v) const;
  bool is_edge(Vertex_handle u, Vertex_handle v,
	       Cell_handle & c, int & i, int & j) const;
  bool is_facet(Vertex_handle u, Vertex_handle v, Vertex_handle w,
		Cell_handle & c, int & i, int & j, int & k) const;
  bool is_cell(Cell_handle c) const;
  bool is_cell(Vertex_handle u, Vertex_handle v, 
	       Vertex_handle w, Vertex_handle t,
	       Cell_handle & c, int & i, int & j, int & k, int & l) const;
  bool is_cell(Vertex_handle u, Vertex_handle v, 
	       Vertex_handle w, Vertex_handle t,
	       Cell_handle & c) const;

  bool has_vertex(const Facet & f, Vertex_handle v, int & j) const;
  bool has_vertex(Cell_handle c, int i, Vertex_handle v, int & j) const;
  bool has_vertex(const Facet & f, Vertex_handle v) const;
  bool has_vertex(Cell_handle c, int i, Vertex_handle v) const;

  bool are_equal(Cell_handle c, int i, Cell_handle n, int j) const;
  bool are_equal(const Facet & f, const Facet & g) const;
  bool are_equal(const Facet & f, Cell_handle n, int j) const;

  Cell_handle
  locate(const Point & p,
	 Locate_type & lt, int & li, int & lj,
	 Cell_handle start = Cell_handle()) const;

  Cell_handle
  locate(const Point & p, Cell_handle start = Cell_handle()) const
  {
      Locate_type lt;
      int li, lj;
      return locate( p, lt, li, lj, start);
  }

  // This one is for backward compatibility with CGAL 2.2.
  Cell_handle
  locate(const Point & p, Cell_handle start,
	 Locate_type & lt, int & li, int & lj) const
  {
      bool WARNING_YOU_ARE_USING_THE_DEPRECATED_VERSION_OF_LOCATE;
      return locate(p, lt, li, lj, start);
  }

  // PREDICATES ON POINTS ``TEMPLATED'' by the geom traits

  Bounded_side
  side_of_tetrahedron(const Point & p,
		      const Point & p0, 
		      const Point & p1,
		      const Point & p2, 
		      const Point & p3,
		      Locate_type & lt, int & i, int & j ) const;
  Bounded_side
  side_of_cell(const Point & p, 
	       Cell_handle c,
	       Locate_type & lt, int & i, int & j) const;
  Bounded_side
  side_of_triangle(const Point & p,
		   const Point & p0, const Point & p1, const Point & p2,
		   Locate_type & lt, int & i, int & j ) const;
  Bounded_side
  side_of_facet(const Point & p,
		Cell_handle c,
		Locate_type & lt, int & li, int & lj) const;
  Bounded_side
  side_of_facet(const Point & p,
		const Facet & f,
		Locate_type & lt, int & li, int & lj) const
    {
      CGAL_triangulation_precondition( f.second == 3 );
      return side_of_facet(p, f.first, lt, li, lj);
    }
  Bounded_side
  side_of_segment(const Point & p, 
		  const Point & p0, const Point & p1,
		  Locate_type & lt, int & i ) const;
  Bounded_side
  side_of_edge(const Point & p, 
	       Cell_handle c,
	       Locate_type & lt, int & li) const;
  Bounded_side
  side_of_edge(const Point & p,
	       const Edge & e,
	       Locate_type & lt, int & li) const
    {
      CGAL_triangulation_precondition( e.second == 0 );
      CGAL_triangulation_precondition( e.third == 1 );
      return side_of_edge(p, e.first, lt, li);
    }

  // MODIFIERS
  bool flip(Facet f);
  bool flip(Cell_handle c, int i);
  void flip_flippable(Facet f);
  void flip_flippable(Cell_handle c, int i);
  bool flip(Edge e);
  bool flip(Cell_handle c, int i, int j);
  void flip_flippable(Edge e);
  void flip_flippable(Cell_handle c, int i, int j);

  //INSERTION 

  Vertex_handle insert(const Point & p, 
		       Cell_handle start = Cell_handle(),
		       Vertex_handle v = NULL);

  Vertex_handle push_back(const Point & p)
  {
      return insert(p);
  }

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

  Vertex_handle
  insert_in_cell(const Point & p, Cell_handle c, Vertex_handle v = NULL);

  Vertex_handle
  insert_in_facet(const Point & p, Cell_handle c, int i,
		  Vertex_handle v = NULL);

  Vertex_handle
  insert_in_facet(const Point & p, const Facet & f,
		  Vertex_handle v = NULL)
    {
      return insert_in_facet(p, f.first,f.second, v);
    }

  Vertex_handle
  insert_in_edge(const Point & p, Cell_handle c, int i, int j,
		 Vertex_handle v = NULL);

  Vertex_handle
  insert_in_edge(const Point & p, const Edge & e,
		 Vertex_handle v = NULL)
    {
      return insert_in_edge(p, e.first,e.second,e.third, v);
    }
  
  Vertex_handle
  insert_outside_convex_hull(const Point & p, Cell_handle c,
			     Vertex_handle v = NULL);
  Vertex_handle
  insert_outside_affine_hull(const Point & p, Vertex_handle v = NULL );

private:
  // Here are the conflit tester function object passed to
  // _tds.insert_conflict() by insert_outside_convex_hull().
  class Conflict_tester_outside_convex_hull_3
  {
      const Point &p;
      Self *t;

  public:

      Conflict_tester_outside_convex_hull_3(const Point &pt, Self *tr)
	  : p(pt), t(tr) {}

      bool operator()(const typename Tds::Cell *c) const
      {
	  Locate_type loc;
          int i, j;
	  return t->side_of_cell( p, (Cell_handle)(Cell*) c, loc, i, j )
	      == ON_BOUNDED_SIDE;
      }
  };

  class Conflict_tester_outside_convex_hull_2
  {
      const Point &p;
      Self *t;

  public:

      Conflict_tester_outside_convex_hull_2(const Point &pt, Self *tr)
	  : p(pt), t(tr) {}

      bool operator()(const typename Tds::Cell *c) const
      {
	  Locate_type loc;
          int i, j;
	  return t->side_of_facet( p, (Cell_handle)(Cell*) c, loc, i, j )
	      == ON_BOUNDED_SIDE;
      }
  };

public:

  //TRAVERSING : ITERATORS AND CIRCULATORS
  Cell_iterator finite_cells_begin() const
    {
      if ( dimension() < 3 )
	  return cells_end();
      return Cell_iterator(this, false); // false means without infinite cells.
    }
  Cell_iterator all_cells_begin() const
    {
      if ( dimension() < 3 )
	  return cells_end();
      return Cell_iterator(this, true); // true means with infinite cells.
    }
  Cell_iterator cells_end() const
    {
      return Cell_iterator(this); // no second argument -> past-end
    }

  Vertex_iterator finite_vertices_begin() const
    {
      if ( number_of_vertices() <= 0 )
	  return vertices_end();
      return Vertex_iterator(this, false);
    }
  Vertex_iterator all_vertices_begin() const
    {
      if ( number_of_vertices() <= 0 )
	  return vertices_end();
      return Vertex_iterator(this, true);
    }
  Vertex_iterator vertices_end() const
    {
      return Vertex_iterator(this);
    }

  Edge_iterator finite_edges_begin() const
    {
      if ( dimension() < 1 )
	  return edges_end();
      return Edge_iterator(this, false);
    }
  Edge_iterator all_edges_begin() const
    {
      if ( dimension() < 1 )
	  return edges_end();
      return Edge_iterator(this, true);
    }
  Edge_iterator edges_end() const
    {
      return Edge_iterator(this);
    }

  Facet_iterator finite_facets_begin() const
    {
      if ( dimension() < 2 )
	  return facets_end();
      return Facet_iterator(this, false);
    }
  Facet_iterator all_facets_begin() const
    {
      if ( dimension() < 2 )
	  return facets_end();
      return Facet_iterator(this, true);
    }
  Facet_iterator facets_end() const
    {
      return Facet_iterator(this);
    }

  // cells around an edge
  Cell_circulator incident_cells(const Edge & e) const
    {
      CGAL_triangulation_precondition( dimension() == 3 );
      return Cell_circulator(this, e);
    }
  Cell_circulator incident_cells(Cell_handle c, int i, int j) const
    {
      CGAL_triangulation_precondition( dimension() == 3 );
      return Cell_circulator(this,c,i,j);
    }
  Cell_circulator incident_cells(const Edge & e, Cell_handle start) const
    {
      CGAL_triangulation_precondition( dimension() == 3 );
      return Cell_circulator(this, e, start);
    }
  Cell_circulator incident_cells(Cell_handle c, int i, int j, 
				 Cell_handle start) const  
    {
      CGAL_triangulation_precondition( dimension() == 3 );
      return Cell_circulator(this, c, i, j, start);
    }

  // facets around an edge
  Facet_circulator incident_facets(const Edge & e) const
    {
      CGAL_triangulation_precondition( dimension() == 3 );
      return Facet_circulator(this, e);
    }
  Facet_circulator incident_facets(Cell_handle c, int i, int j) const
    {
      CGAL_triangulation_precondition( dimension() == 3 );
      return Facet_circulator(this, c, i, j);
    }
  Facet_circulator incident_facets(const Edge & e, 
				   const Facet & start) const
    {
      CGAL_triangulation_precondition( dimension() == 3 );
      return Facet_circulator(this, e, start);
    }
  Facet_circulator incident_facets(Cell_handle c, int i, int j, 
				   const Facet & start) const  
    {
      CGAL_triangulation_precondition( dimension() == 3 );
      return Facet_circulator(this, c, i, j, start);
    }
  Facet_circulator incident_facets(const Edge & e, 
				   Cell_handle start, int f) const
    {
      CGAL_triangulation_precondition( dimension() == 3 );
      return Facet_circulator(this, e, start, f);
    }
  Facet_circulator incident_facets(Cell_handle c, int i, int j, 
				   Cell_handle start, int f) const  
    {
      CGAL_triangulation_precondition( dimension() == 3 );
      return Facet_circulator(this, c, i, j, start, f);
    }

  // around a vertex
  void
  incident_cells(Vertex_handle v,
                 std::set<Cell_handle> & cells,
		 Cell_handle c = Cell_handle() ) const;

  void
  incident_vertices(Vertex_handle v,
                    std::set<Vertex_handle> & vertices,
		    Cell_handle c = Cell_handle() ) const;

  // old methods, kept for compatibility with previous versions
  void
  incident_cells(Vertex_handle v, 
		 std::set<Cell*> & cells,
		 Cell_handle c = Cell_handle(),
		 int dummy_for_windows = 0) const;

  void
  incident_vertices(Vertex_handle v, 
		    std::set<Vertex*> & vertices,
		    Cell_handle c = Cell_handle(),
		    int dummy_for_windows = 0) const;

private:
  void 
  util_incident_vertices(Vertex_handle v, 
			 std::set<Vertex_handle> & vertices,
			 std::set<Cell_handle> & cells,
			 Cell_handle c ) const;
  void 
  util_incident_vertices(Vertex_handle v, 
			 std::set<Vertex*> & vertices,
			 std::set<Cell*> & cells,
			 Cell_handle c,
			 int dummy_for_windows = 0) const;

public:

  // CHECKING
  bool is_valid(bool verbose = false, int level = 0) const;

  bool is_valid(Cell_handle c, bool verbose = false, int level = 0) const;

  bool is_valid_finite(Cell_handle c, 
		       bool verbose = false, int level = 0) const;
};


template < class GT, class Tds >
std::istream & 
operator>> (std::istream& is, Triangulation_3<GT, Tds> &tr)
  // reads
  // the dimension
  // the number of finite vertices
  // the non combinatorial information on vertices (point, etc)
  // the number of cells
  // the cells by the indices of their vertices in the preceding list
  // of vertices, plus the non combinatorial information on each cell
  // the neighbors of each cell by their index in the preceding list of cells
  // when dimension < 3 : the same with faces of maximal dimension
{
  typedef Triangulation_3<GT, Tds>               Triangulation;
  typedef typename Triangulation::Vertex_handle  Vertex_handle;
  typedef typename Triangulation::Cell_handle    Cell_handle;
  typedef typename Triangulation::Vertex         Vertex;
  typedef typename Tds::Vertex                   TdsVertex;
  typedef typename Tds::Cell                     TdsCell;

  tr._tds.clear(); // infinite vertex deleted
  tr.infinite = (Vertex *) tr._tds.create_vertex();

  int n, d;
  is >> d >> n;
  tr._tds.set_dimension(d);
  tr.set_number_of_vertices(n);

  std::map< int, TdsVertex* > V;
  V[0] = &*(tr.infinite_vertex());
  // the infinite vertex is numbered 0

  for (int i=1; i <= n; i++) {
    V[i] = tr._tds.create_vertex();
    is >> *V[i];
  }

  std::map< int, TdsCell* > C;

  int m;
  tr._tds.read_cells(is, V, m, C);

  for (int j=0 ; j < m; j++)
    is >> *(C[j]);

  CGAL_triangulation_assertion( tr.is_valid(false) );
  return is;
}
    
template < class GT, class Tds >
std::ostream & 
operator<< (std::ostream& os, const Triangulation_3<GT, Tds> &tr)
  // writes :
  // the dimension
  // the number of finite vertices
  // the non combinatorial information on vertices (point, etc)
  // the number of cells
  // the cells by the indices of their vertices in the preceding list
  // of vertices, plus the non combinatorial information on each cell
  // the neighbors of each cell by their index in the preceding list of cells
  // when dimension < 3 : the same with faces of maximal dimension
{
  typedef Triangulation_3<GT, Tds>                 Triangulation;
  typedef typename Triangulation::Vertex_handle    Vertex_handle;
  typedef typename Triangulation::Vertex_iterator  Vertex_iterator;
  typedef typename Triangulation::Cell_iterator    Cell_iterator;
  typedef typename Triangulation::Edge_iterator    Edge_iterator;
  typedef typename Triangulation::Facet_iterator   Facet_iterator;

  // outputs dimension and number of vertices
  int n = tr.number_of_vertices();
  if (is_ascii(os))
    os << tr.dimension() << std::endl << n << std::endl;
  else
    os << tr.dimension() << n;

  if (n == 0)
    return os;
 
  std::vector<Vertex_handle> TV(n+1);
  int i = 0;

  // write the vertices
  // the vertices must be indexed by their order of creation so
  // that when reread from file, the orders of vertices are the
  // same - important for remove 

  for (Vertex_iterator it=tr.all_vertices_begin(); it!=tr.vertices_end(); ++it)
    TV[i++] = &(*it);

  CGAL_triangulation_assertion( i == n+1 ); 

  std::sort(TV.begin(), TV.end(), 
	    Vertex_tds_compare_order_of_creation<Vertex_handle>()); 

  CGAL_triangulation_assertion( tr.is_infinite(TV[0]) );

  std::map< void*, int > V;

  V[&(*tr.infinite_vertex())] = 0;
  for (i=1; i <= n; i++) {
    os << *TV[i];
    V[&(*TV[i])] = i;
    if (is_ascii(os))
	os << std::endl;
  }

  // write the non combinatorial information on the cells
  // using the << operator of Cell
  // works because the iterator of the tds traverses the cells in the
  // same order as the iterator of the triangulation
  switch ( tr.dimension() ) {
  case 3:
    {
      for(Cell_iterator it=tr.all_cells_begin(); it != tr.cells_end(); ++it)
	os << *it; // other information
      break;
    }
  case 2:
    {
      for(Facet_iterator it=tr.all_facets_begin(); it != tr.facets_end(); ++it)
	os << *((*it).first); // other information
      break;
    }
  case 1:
    {
      for(Edge_iterator it=tr.all_edges_begin(); it != tr.edges_end(); ++it)
	os << *((*it).first); // other information 
      break;
    }
  }

  // asks the tds for the combinatorial information 
  tr.tds().print_cells(os, V);
  
  return os ;
}

template < class GT, class Tds >
bool
Triangulation_3<GT,Tds>::
test_dim_down(Vertex_handle v)
  // tests whether removing v decreases the dimension of the triangulation 
  // true iff
  // v is incident to all finite cells
  // and all the other vertices are coplanar
{
  CGAL_triangulation_precondition(dimension() == 3);
  CGAL_triangulation_precondition(! is_infinite(v) );

  Cell_iterator cit = finite_cells_begin();
  Cell_iterator cdone = cells_end();

  int i, iv;
  if ( ! cit->has_vertex(v,iv) ) return false;
  Point p1=cit->vertex((iv+1)&3)->point();  
  Point p2=cit->vertex((iv+2)&3)->point();  
  Point p3=cit->vertex((iv+3)&3)->point();
  ++cit;
  
  while ( cit != cdone ) {
    if ( ! cit->has_vertex(v,iv) )
      return false;
    for ( i=1; i<4; i++ ) {
	if ( orientation
	     (p1,p2,p3,cit->vertex((iv+i)&3)->point()) != COPLANAR )
	  return false;
    }
    ++cit;
  }

  return true;
}// test_dim_down

template < class GT, class Tds >
int
Triangulation_3<GT,Tds>::
number_of_finite_cells() const 
{ 
  if ( dimension() < 3 ) return 0;
  return std::distance(finite_cells_begin(), cells_end());
}
  
template < class GT, class Tds >
int
Triangulation_3<GT,Tds>::
number_of_cells() const 
{ 
  if ( dimension() < 3 ) return 0;
  return std::distance(all_cells_begin(), cells_end());
}

template < class GT, class Tds >
int
Triangulation_3<GT,Tds>::
number_of_finite_facets() const
{
  if ( dimension() < 2 ) return 0;
  return std::distance(finite_facets_begin(), facets_end());
}

template < class GT, class Tds >
int
Triangulation_3<GT,Tds>::
number_of_facets() const
{
  if ( dimension() < 2 ) return 0;
  return std::distance(all_facets_begin(), facets_end());
}

template < class GT, class Tds >
int
Triangulation_3<GT,Tds>::
number_of_finite_edges() const
{
  if ( dimension() < 1 ) return 0;
  return std::distance(finite_edges_begin(), edges_end());
}

template < class GT, class Tds >
int
Triangulation_3<GT,Tds>::
number_of_edges() const
{
  if ( dimension() < 1 ) return 0;
  return std::distance(all_edges_begin(), edges_end());
}

template < class GT, class Tds >
Triangulation_3<GT,Tds>::Triangle
Triangulation_3<GT,Tds>::
triangle(const Cell_handle c, int i) const
{ 
  CGAL_triangulation_precondition( dimension() == 2 || dimension() == 3 );
  CGAL_triangulation_precondition( (dimension() == 2 && i == 3)
				   || (dimension() == 3 && i >= 0 && i <= 3) );
  CGAL_triangulation_precondition( ! is_infinite(std::make_pair(c, i)) );
  return construct_triangle(c->vertex(i<=0 ? 1 : 0)->point(),
			    c->vertex(i<=1 ? 2 : 1)->point(),
			    c->vertex(i<=2 ? 3 : 2)->point());
}

template < class GT, class Tds >
Triangulation_3<GT,Tds>::Segment
Triangulation_3<GT,Tds>::
segment(const Cell_handle c, int i, int j) const
{
  CGAL_triangulation_precondition( i != j );
  CGAL_triangulation_precondition( dimension() >= 1 && dimension() <= 3 );
  CGAL_triangulation_precondition( i >= 0 && i <= dimension() 
				   && j >= 0 && j <= dimension() );
  CGAL_triangulation_precondition( ! is_infinite(make_triple(c, i, j)) );
  return construct_segment( c->vertex(i)->point(), c->vertex(j)->point() );
}

template < class GT, class Tds >
inline
bool
Triangulation_3<GT,Tds>::
is_infinite(const Cell_handle c, int i) const 
{
  CGAL_triangulation_precondition( dimension() == 2 || dimension() == 3 );
  CGAL_triangulation_precondition( (dimension() == 2 && i == 3)
				   || (dimension() == 3 && i >= 0 && i <= 3) );
  return is_infinite(c->vertex(i<=0 ? 1 : 0)) ||
	 is_infinite(c->vertex(i<=1 ? 2 : 1)) ||
	 is_infinite(c->vertex(i<=2 ? 3 : 2));
}

template < class GT, class Tds >
inline
bool
Triangulation_3<GT,Tds>::
is_infinite(const Cell_handle c, int i, int j) const 
{ 
  CGAL_triangulation_precondition( i != j );
  CGAL_triangulation_precondition( dimension() >= 1 && dimension() <= 3 );
  CGAL_triangulation_precondition(
	  i >= 0 && i <= dimension() && j >= 0 && j <= dimension() );
  return is_infinite( c->vertex(i) ) || is_infinite( c->vertex(j) );
}

template < class GT, class Tds >
bool
Triangulation_3<GT,Tds>::
is_vertex(const Point & p, Vertex_handle & v) const
{
  Locate_type lt;
  int li, lj;
  Cell_handle c = locate( p, lt, li, lj );
  if ( lt != VERTEX )
    return false;
  v = c->vertex(li);
  return true;
}

template < class GT, class Tds >
inline
bool
Triangulation_3<GT,Tds>::
is_vertex(Vertex_handle v) const
{
  return _tds.is_vertex(&(*v));
}

template < class GT, class Tds >
bool
Triangulation_3<GT,Tds>::
is_edge(Vertex_handle u, Vertex_handle v,
	Cell_handle & c, int & i, int & j) const
{
  CGAL_triangulation_expensive_precondition( _tds.is_vertex(&(*u)) &&
					     _tds.is_vertex(&(*v)) );
  typename Tds::Cell* cstar;
  bool b = _tds.is_edge(&(*u), &(*v), cstar, i, j);
  if (b)
      c = ((Cell*) cstar)->handle();
  return b;
}

template < class GT, class Tds >
bool
Triangulation_3<GT,Tds>::
is_facet(Vertex_handle u, Vertex_handle v, Vertex_handle w,
	 Cell_handle & c, int & i, int & j, int & k) const
{
  CGAL_triangulation_expensive_precondition( _tds.is_vertex(&(*u)) &&
					     _tds.is_vertex(&(*v)) &&
					     _tds.is_vertex(&(*w)) );
  typename Tds::Cell* cstar;
  bool b = _tds.is_facet(&(*u), &(*v), &(*w), cstar, i, j, k);
  if (b)
      c = ((Cell*) cstar)->handle();
  return b;
}

template < class GT, class Tds >
inline
bool
Triangulation_3<GT,Tds>::
is_cell(Cell_handle c) const
{
  return _tds.is_cell(&(*c));
}

template < class GT, class Tds >
bool
Triangulation_3<GT,Tds>::
is_cell(Vertex_handle u, Vertex_handle v, 
	Vertex_handle w, Vertex_handle t,
	Cell_handle & c, int & i, int & j, int & k, int & l) const
{
  CGAL_triangulation_expensive_precondition( _tds.is_vertex(&(*u)) &&
					     _tds.is_vertex(&(*v)) &&
					     _tds.is_vertex(&(*w)) &&
					     _tds.is_vertex(&(*t)) );
  typename Tds::Cell* cstar;
  bool b = _tds.is_cell(&(*u), &(*v), &(*w), &(*t), cstar, i, j, k, l);
  if (b)
      c = ((Cell*) cstar)->handle();
  return b;
}

template < class GT, class Tds >
bool
Triangulation_3<GT,Tds>::
is_cell(Vertex_handle u, Vertex_handle v, 
	Vertex_handle w, Vertex_handle t,
	Cell_handle & c) const
{
  CGAL_triangulation_expensive_precondition( _tds.is_vertex(&(*u)) &&
					     _tds.is_vertex(&(*v)) &&
					     _tds.is_vertex(&(*w)) &&
					     _tds.is_vertex(&(*t)) );
  int i,j,k,l;
  typename Tds::Cell* cstar;
  bool b = _tds.is_cell(&(*u), &(*v), &(*w), &(*t), cstar, i, j, k, l);
  if (b)
      c = ((Cell*) cstar)->handle();
  return b;
}

template < class GT, class Tds >
inline
bool
Triangulation_3<GT,Tds>::
has_vertex(const Facet & f, Vertex_handle v, int & j) const
{
  return _tds.has_vertex(&*(f.first), f.second, &*v, j);
}

template < class GT, class Tds >
inline
bool
Triangulation_3<GT,Tds>::
has_vertex(Cell_handle c, int i, Vertex_handle v, int & j) const
{
  return _tds.has_vertex(&*c, i, &*v, j);
}

template < class GT, class Tds >
inline
bool
Triangulation_3<GT,Tds>::
has_vertex(const Facet & f, Vertex_handle v) const
{
  return _tds.has_vertex(&*(f.first), f.second, &*v);
}

template < class GT, class Tds >
inline
bool
Triangulation_3<GT,Tds>::
has_vertex(Cell_handle c, int i, Vertex_handle v) const
{
  return _tds.has_vertex(&*c, i, &*v);
}

template < class GT, class Tds >
inline
bool
Triangulation_3<GT,Tds>::
are_equal(Cell_handle c, int i, Cell_handle n, int j) const
{
  return _tds.are_equal(&*c, i, &*n, j);
}

template < class GT, class Tds >
inline
bool
Triangulation_3<GT,Tds>::
are_equal(const Facet & f, const Facet & g) const
{
  return _tds.are_equal(&*(f.first), f.second, &*(g.first), g.second);
}

template < class GT, class Tds >
inline
bool
Triangulation_3<GT,Tds>::
are_equal(const Facet & f, Cell_handle n, int j) const
{
  return _tds.are_equal(&*(f.first), f.second, &*n, j);
}

template < class GT, class Tds >
Triangulation_3<GT,Tds>::Cell_handle
Triangulation_3<GT,Tds>::
locate(const Point & p, Locate_type & lt, int & li, int & lj,
       Cell_handle start ) const
  // returns the (finite or infinite) cell p lies in
  // starts at cell "start"
  // start must be finite
  // if lt == OUTSIDE_CONVEX_HULL, li is the
  // index of a facet separating p from the rest of the triangulation
  // in dimension 2 :
  // returns a facet (Cell_handle,li) if lt == FACET
  // returns an edge (Cell_handle,li,lj) if lt == EDGE
  // returns a vertex (Cell_handle,li) if lt == VERTEX
  // if lt == OUTSIDE_CONVEX_HULL, li, lj give the edge of c
  // separating p from the rest of the triangulation
  // lt = OUTSIDE_AFFINE_HULL if p is not coplanar with the triangulation
{
  int i, inf;

  if ( dimension() >= 1 && start == Cell_handle() )
    // there is at least one finite "cell" (or facet or edge)
    start = infinite_vertex()->cell()->neighbor
            ( infinite_vertex()->cell()->index( infinite_vertex()) );

  switch (dimension()) {
  case 3:
    {
      CGAL_triangulation_precondition( start != Cell_handle() );
      Cell_handle c, previous;
      int ind_inf;
      if ( start->has_vertex(infinite, ind_inf) )
	c = start->neighbor(ind_inf);
      else
	c = start;
 
      Orientation o[4];

      // We implement the remembering visibility/stochastic walk.

      // Main locate loop
      while(1) {
	if ( c->has_vertex(infinite,li) ) {
	  // c must contain p in its interior
	  lt = OUTSIDE_CONVEX_HULL;
	  return c;
	}

	// FIXME: do more benchmarks.
	i = rand_4(); // For the (remembering) stochastic walk
	// i = 0; // For the (remembering) visibility walk. Ok for Delaunay only

        Orientation test_or = (i&1)==0 ? NEGATIVE : POSITIVE;
	const Point & p0 = c->vertex( i )->point();
	const Point & p1 = c->vertex( (i+1)&3 )->point();
	const Point & p2 = c->vertex( (i+2)&3 )->point();
	const Point & p3 = c->vertex( (i+3)&3 )->point();

	// Note : among the four Points, 3 are common with the previous cell...
	// Something can probably be done to take advantage of this, like
	// storing the four in an array and changing only one ?

	// We could make a loop of these 4 blocks, for clarity, but not speed.
	Cell_handle next = c->neighbor(i);
	if (previous != next) {
	  o[0] = orientation(p, p1, p2, p3);
	  if ( o[0] == test_or) {
	    previous = c;
	    c = next;
	    continue;
	  }
	} else
	    o[0] = (Orientation) - test_or;

	next = c->neighbor((i+1)&3);
	if (previous != next) {
	  o[1] = orientation(p0, p, p2, p3);
	  if ( o[1] == test_or) {
	    previous = c;
	    c = next;
	    continue;
	  }
	} else
	    o[1] = (Orientation) - test_or;

	next = c->neighbor((i+2)&3);
	if (previous != next) {
	  o[2] = orientation(p0, p1, p, p3);
	  if ( o[2] == test_or) {
	    previous = c;
	    c = next;
	    continue;
	  }
	} else
	    o[2] = (Orientation) - test_or;

	next = c->neighbor((i+3)&3);
	if (previous != next) {
	  o[3] = orientation(p0, p1, p2, p);
	  if ( o[3] == test_or) {
	    // previous = c; // not necessary because it's the last one.
	    c = next;
	    continue;
	  }
	} else
	    o[3] = (Orientation) - test_or;

	break;
      }
	  
	// now p is in c or on its boundary
	int sum = ( o[0] == COPLANAR )
	  + ( o[1] == COPLANAR )
	  + ( o[2] == COPLANAR )
	  + ( o[3] == COPLANAR );
	switch (sum) {
	case 0:
	  {
	    lt = CELL;
	    break;
	  }
	case 1:
	  { 
	    lt = FACET;
	    li = ( o[0] == COPLANAR ) ? i :
	         ( o[1] == COPLANAR ) ? (i+1)&3 :
	         ( o[2] == COPLANAR ) ? (i+2)&3 :
	         (i+3)&3;
	    break;
	  }
	case 2:
	  { 
	    lt = EDGE;
	    li = ( o[0] != COPLANAR ) ? i :
	         ( o[1] != COPLANAR ) ? ((i+1)&3) :
	         ((i+2)&3);
	    lj = ( o[ (li+1-i)&3 ] != COPLANAR ) ? ((li+1)&3) :
	         ( o[ (li+2-i)&3 ] != COPLANAR ) ? ((li+2)&3) :
	         ((li+3)&3);
	    CGAL_triangulation_assertion(collinear( p,
						    c->vertex( li )->point(),
						    c->vertex( lj )->point() ));
	    break;
	  }
	case 3:
	  {
	    lt = VERTEX;
	    li = ( o[0] != COPLANAR ) ? i :
	         ( o[1] != COPLANAR ) ? (i+1)&3 :
	         ( o[2] != COPLANAR ) ? (i+2)&3 :
	         (i+3)&3;
	    break;
	  }
	}
	return c;
    }
  case 2:
    {
      CGAL_triangulation_precondition( start != Cell_handle() );
      Cell_handle c;
      int ind_inf;
      if ( start->has_vertex(infinite, ind_inf) )
	c = start->neighbor(ind_inf);
      else
	c = start;

      //first tests whether p is coplanar with the current triangulation
      Facet_iterator finite_fit = finite_facets_begin();
      if ( orientation( (*finite_fit).first->vertex(0)->point(),
			(*finite_fit).first->vertex(1)->point(),
			(*finite_fit).first->vertex(2)->point(),
			p ) != DEGENERATE ) {
	lt = OUTSIDE_AFFINE_HULL;
	li = 3; // only one facet in dimension 2
	return (*finite_fit).first;
      }
      // if p is coplanar, location in the triangulation
      // only the facet numbered 3 exists in each cell
      while (1) {
	  
	if ( c->has_vertex(infinite,inf) ) {
	  // c must contain p in its interior
	  lt = OUTSIDE_CONVEX_HULL;
	  li = cw(inf);
	  lj = ccw(inf);
	  return c;
	}

	// else c is finite
	// we test its edges in a random order until we find a
	// neighbor to go further
	i = rand_3();
	const Point & p0 = c->vertex( i )->point();
	const Point & p1 = c->vertex( ccw(i) )->point();
	const Point & p2 = c->vertex( cw(i) )->point();
        Orientation o[3];
	CGAL_triangulation_assertion(coplanar_orientation(p0,p1,p2)==POSITIVE);
	o[0] = coplanar_orientation(p0,p1,p);
	if ( o[0] == NEGATIVE ) {
	  c = c->neighbor( cw(i) );
	  continue;
	}
	o[1] = coplanar_orientation(p1,p2,p);
	if ( o[1] == NEGATIVE ) {
	  c = c->neighbor( i );
	  continue;
	}
	o[2] = coplanar_orientation(p2,p0,p);
	if ( o[2] == NEGATIVE ) {
	  c = c->neighbor( ccw(i) );
	  continue;
	}

	// now p is in c or on its boundary
	int sum = ( o[0] == COLLINEAR )
	        + ( o[1] == COLLINEAR )
	        + ( o[2] == COLLINEAR );
	switch (sum) {
	case 0:
	  {
	    lt = FACET;
	    li = 3; // useless ?
	    break;
	  }
	case 1:
	  {
	    lt = EDGE;
	    li = ( o[0] == COLLINEAR ) ? i :
	         ( o[1] == COLLINEAR ) ? ccw(i) :
	         cw(i);
	    lj = ccw(li);
	    break;
	  }
	case 2:
	  {
	    lt = VERTEX;
	    li = ( o[0] != COLLINEAR ) ? cw(i) :
	         ( o[1] != COLLINEAR ) ? i :
	         ccw(i);
	    break;
	  }
	}
	return c;
      }
    }
  case 1:
    {
      CGAL_triangulation_precondition( start != Cell_handle() );
      Cell_handle c;
      int ind_inf;
      if ( start->has_vertex(infinite, ind_inf) )
	c = start->neighbor(ind_inf);
      else
	c = start;

      //first tests whether p is collinear with the current triangulation
      Edge_iterator finite_eit = finite_edges_begin();
      if ( ! collinear( p,
			(*finite_eit).first->vertex(0)->point(),
			(*finite_eit).first->vertex(1)->point()) ) {
	lt = OUTSIDE_AFFINE_HULL;
	return (*finite_eit).first;
      }
      // if p is collinear, location :
      Comparison_result o, o0, o1;
      int xyz;
      Point p0 = c->vertex(0)->point();
      Point p1 = c->vertex(1)->point();
      CGAL_triangulation_assertion( ( compare_x(p0,p1) != EQUAL ) ||
				    ( compare_y(p0,p1) != EQUAL ) ||
				    ( compare_z(p0,p1) != EQUAL ) );
      o = compare_x(p0,p1);
      if ( o == EQUAL ) {
	o = compare_y(p0,p1);
	if ( o == EQUAL ) {
	  o = compare_z(p0,p1);
	  xyz = 3;
	}
	else 
	  xyz = 2;
      }
      else 
	xyz  = 1;
      //	bool notfound = true;
      while (1) {
	if ( c->has_vertex(infinite,inf) ) {
	  // c must contain p in its interior
	  lt = OUTSIDE_CONVEX_HULL;
	  return c;
	}

	// else c is finite
	// we test on which direction to continue the traversal
	p0 = c->vertex(0)->point();
	p1 = c->vertex(1)->point();
	switch ( xyz ) {
	case 1:
	  {
	    o = compare_x(p0,p1);
	    o0 = compare_x(p0,p);
	    o1 = compare_x(p,p1);
	    break;
	  }
	case 2:
	  {
	    o = compare_y(p0,p1);
	    o0 = compare_y(p0,p);
	    o1 = compare_y(p,p1);
	    break;
	  }
	default: // case 3
	  {
	    o = compare_z(p0,p1);
	    o0 = compare_z(p0,p);
	    o1 = compare_z(p,p1);
	  }
	}
	  
	if (o0 == EQUAL) {
	  lt = VERTEX;
	  li = 0;
	  return c;
	}
	if (o1 == EQUAL) {
	  lt = VERTEX;
	  li = 1;
	  return c;
	}
	if ( o0 == o1 ) {
	  lt = EDGE;
	  li = 0;
	  lj = 1;
	  return c;
	}
	if ( o0 == o ) { 
	  c = c->neighbor(0);
	  continue;
	}
	if ( o1 == o ) { 
	  c = c->neighbor(1);
	  continue; 
	}
      }
    }
  case 0:
    {
      Vertex_iterator vit = finite_vertices_begin();
      if ( ! equal( p, vit->point() ) ) {
	lt = OUTSIDE_AFFINE_HULL;
      }
      else {
	lt = VERTEX;
	li = 0;
      }
      return vit->cell();
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
}
	  
template < class GT, class Tds >
Bounded_side
Triangulation_3<GT,Tds>::
side_of_tetrahedron(const Point & p,
		    const Point & p0, 
		    const Point & p1,
		    const Point & p2, 
		    const Point & p3,
		    Locate_type & lt, int & i, int & j ) const
  // p0,p1,p2,p3 supposed to be non coplanar
  // tetrahedron p0,p1,p2,p3 is supposed to be well oriented
  // returns :
  // ON_BOUNDED_SIDE if p lies strictly inside the tetrahedron
  // ON_BOUNDARY if p lies on one of the facets
  // ON_UNBOUNDED_SIDE if p lies strictly outside the tetrahedron
{
  CGAL_triangulation_precondition
    ( orientation(p0,p1,p2,p3) == POSITIVE );

  Orientation o0,o1,o2,o3;
  if ( ((o0 = orientation(p,p1,p2,p3)) == NEGATIVE) ||
       ((o1 = orientation(p0,p,p2,p3)) == NEGATIVE) ||
       ((o2 = orientation(p0,p1,p,p3)) == NEGATIVE) ||
       ((o3 = orientation(p0,p1,p2,p)) == NEGATIVE) ) {
    lt = OUTSIDE_CONVEX_HULL;
    return ON_UNBOUNDED_SIDE;
  }

  // now all the oi's are >=0
  // sum gives the number of facets p lies on
  int sum = ( (o0 == ZERO) ? 1 : 0 ) 
          + ( (o1 == ZERO) ? 1 : 0 ) 
          + ( (o2 == ZERO) ? 1 : 0 ) 
          + ( (o3 == ZERO) ? 1 : 0 );

  switch (sum) {
  case 0:
    {
      lt = CELL;
      return ON_BOUNDED_SIDE;
    }
  case 1:
    {
      lt = FACET;
      // i = index such that p lies on facet(i)
      i = ( o0 == ZERO ) ? 0 :
	  ( o1 == ZERO ) ? 1 :
	  ( o2 == ZERO ) ? 2 :
	  3;
      return ON_BOUNDARY;
    }
  case 2:
    {
      lt = EDGE;
      // i = smallest index such that p does not lie on facet(i)
      // i must be < 3 since p lies on 2 facets
      i = ( o0 == POSITIVE ) ? 0 :
	  ( o1 == POSITIVE ) ? 1 :
	  2;
      // j = larger index such that p not on facet(j)
      // j must be > 0 since p lies on 2 facets
      j = ( o3 == POSITIVE ) ? 3 :
	  ( o2 == POSITIVE ) ? 2 :
	  1;
      return ON_BOUNDARY;
    }
  case 3:
    {
      lt = VERTEX;
      // i = index such that p does not lie on facet(i)
      i = ( o0 == POSITIVE ) ? 0 :
	  ( o1 == POSITIVE ) ? 1 :
	  ( o2 == POSITIVE ) ? 2 :
	  3;
      return ON_BOUNDARY;
    }
  default:
    {
      // impossible : cannot be on 4 facets for a real tetrahedron
      CGAL_triangulation_assertion(false);
      return ON_BOUNDARY;
    }
  }
}

template < class GT, class Tds >
Bounded_side
Triangulation_3<GT,Tds>::
side_of_cell(const Point & p, 
	     Cell_handle c,
	     Locate_type & lt, int & i, int & j) const
  // returns
  // ON_BOUNDED_SIDE if p inside the cell
  // (for an infinite cell this means that p lies strictly in the half space
  // limited by its finite facet)
  // ON_BOUNDARY if p on the boundary of the cell
  // (for an infinite cell this means that p lies on the *finite* facet)
  // ON_UNBOUNDED_SIDE if p lies outside the cell
  // (for an infinite cell this means that p is not in the preceding
  // two cases)  
  // lt has a meaning only when ON_BOUNDED_SIDE or ON_BOUNDARY
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
    Orientation o;
    Vertex_handle 
      v1=c->vertex((inf+1)&3), 
      v2=c->vertex((inf+2)&3), 
      v3=c->vertex((inf+3)&3);
    if ( (inf&1) == 0 ) 
      o = orientation(p, v1->point(), v2->point(), v3->point());
    else 
      o =  orientation(v3->point(), p, v1->point(), v2->point());

    switch (o) {
    case POSITIVE:
      {
	lt = CELL;
	return ON_BOUNDED_SIDE;
      }
    case NEGATIVE:
      return ON_UNBOUNDED_SIDE;
    case ZERO:
      {
	// location in the finite facet
	int i_f, j_f;
	Bounded_side side = 
	  side_of_triangle(p, v1->point(), v2->point(), v3->point(),
			   lt, i_f, j_f);
	// lt need not be modified in most cases :
	switch (side) {
	case ON_BOUNDED_SIDE:
	  {
	    // lt == FACET ok
	    i = inf;
	    return ON_BOUNDARY;
	  }
	case ON_BOUNDARY:
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
	    return ON_BOUNDARY;
	  }
	case ON_UNBOUNDED_SIDE:
	  {
	    // p lies on the plane defined by the finite facet
	    // lt must be initialized
	    return ON_UNBOUNDED_SIDE;
	  }
	default:
	  {
	    CGAL_triangulation_assertion(false);
	    return ON_BOUNDARY;
	  }
	} // switch side
      }// case ZERO
    default:
      {
	CGAL_triangulation_assertion(false);
	return ON_BOUNDARY;
      }
    } // switch o
  } // else infinite cell
} // side_of_cell

template < class GT, class Tds >
Bounded_side
Triangulation_3<GT,Tds>::
side_of_triangle(const Point & p,
		 const Point & p0, 
		 const Point & p1,
		 const Point & p2,
		 Locate_type & lt, int & i, int & j ) const
  // p0,p1,p2 supposed to define a plane
  // p supposed to lie on plane p0,p1,p2
  // triangle p0,p1,p2 defines the orientation of the plane
  // returns
  // ON_BOUNDED_SIDE if p lies strictly inside the triangle
  // ON_BOUNDARY if p lies on one of the edges
  // ON_UNBOUNDED_SIDE if p lies strictly outside the triangle
{
  CGAL_triangulation_precondition( orientation(p,p0,p1,p2) == COPLANAR );

  Orientation o012 = coplanar_orientation(p0,p1,p2);
  CGAL_triangulation_precondition( o012 != COLLINEAR );

  Orientation o0; // edge p0 p1
  Orientation o1; // edge p1 p2
  Orientation o2; // edge p2 p0

  if ((o0 = coplanar_orientation(p0,p1,p)) == opposite(o012) ||
      (o1 = coplanar_orientation(p1,p2,p)) == opposite(o012) ||
      (o2 = coplanar_orientation(p2,p0,p)) == opposite(o012)) {
    lt = OUTSIDE_CONVEX_HULL;
    return ON_UNBOUNDED_SIDE;
  }

  // now all the oi's are >=0
  // sum gives the number of edges p lies on
  int sum = ( (o0 == ZERO) ? 1 : 0 ) 
          + ( (o1 == ZERO) ? 1 : 0 ) 
          + ( (o2 == ZERO) ? 1 : 0 );

  switch (sum) {
  case 0:
    {
      lt = FACET;
      return ON_BOUNDED_SIDE;
    }
  case 1:
    {
      lt = EDGE;
      i = ( o0 == ZERO ) ? 0 :
	  ( o1 == ZERO ) ? 1 :
	  2;
      if ( i == 2 ) 
	j=0;
      else 
	j = i+1;
      return ON_BOUNDARY;
    }
  case 2:
    {
      lt = VERTEX;
      i = ( o0 == o012 ) ? 2 :
	  ( o1 == o012 ) ? 0 :
	  1;
      return ON_BOUNDARY;
    }
  default:
    {
      // cannot happen
      CGAL_triangulation_assertion(false);
      return ON_BOUNDARY;
    }
  }
}

template < class GT, class Tds >
Bounded_side
Triangulation_3<GT,Tds>::
side_of_facet(const Point & p,
	      Cell_handle c,
	      Locate_type & lt, int & li, int & lj) const
  // supposes dimension 2 otherwise does not work for infinite facets
  // returns :
  // ON_BOUNDED_SIDE if p inside the facet
  // (for an infinite facet this means that p lies strictly in the half plane
  // limited by its finite edge)
  // ON_BOUNDARY if p on the boundary of the facet
  // (for an infinite facet this means that p lies on the *finite* edge)
  // ON_UNBOUNDED_SIDE if p lies outside the facet
  // (for an infinite facet this means that p is not in the
  // preceding two cases) 
  // lt has a meaning only when ON_BOUNDED_SIDE or ON_BOUNDARY
  // when they mean anything, li and lj refer to indices in the cell c 
  // giving the facet (c,i)
{
  CGAL_triangulation_precondition( dimension() == 2 );
  if ( ! is_infinite(c,3) ) {
    // The following precondition is useless because it is written
    // in side_of_facet  
    // 	CGAL_triangulation_precondition( orientation (p, 
    // 					  c->vertex(0)->point,
    // 					  c->vertex(1)->point,
    // 					  c->vertex(2)->point) == COPLANAR );
    int i_t, j_t;
    Bounded_side side = side_of_triangle(p,
			    c->vertex(0)->point(),
			    c->vertex(1)->point(),
			    c->vertex(2)->point(),
			    lt, i_t, j_t);
    // indices in the original cell :
    li = ( i_t == 0 ) ? 0 :
         ( i_t == 1 ) ? 1 : 2;
    lj = ( j_t == 0 ) ? 0 :
         ( j_t == 1 ) ? 1 : 2;
    return side;
  }
  // else infinite facet
  int inf = c->index(infinite);
    // The following precondition is useless because it is written
    // in side_of_facet  
    // 	CGAL_triangulation_precondition( orientation (p,
    // 					  c->neighbor(inf)->vertex(0)->point(),
    // 					  c->neighbor(inf)->vertex(1)->point(),
    // 					  c->neighbor(inf)->vertex(2)->point())
    // 					 == COPLANAR );
  int i2 = next_around_edge(inf,3);
  int i1 = 3-inf-i2;
  Vertex_handle v1 = c->vertex(i1),
                v2 = c->vertex(i2);

  CGAL_triangulation_assertion(coplanar_orientation(v1->point(), v2->point(),
	                       c->mirror_vertex(inf)->point()) == POSITIVE);
  Orientation o = coplanar_orientation(v1->point(), v2->point(), p);
  switch (o) {
  case POSITIVE:
    // p lies on the same side of v1v2 as vn, so not in f
    {
      return ON_UNBOUNDED_SIDE;
    }
  case NEGATIVE:
    // p lies in f
    { 
      lt = FACET;
      li = 3;
      return ON_BOUNDED_SIDE;
    }
  case ZERO:
    // p collinear with v1v2
    {
      int i_e;
      switch (side_of_segment(p, v1->point(), v2->point(), lt, i_e)) {
	// computation of the indices in the original cell
      case ON_BOUNDED_SIDE:
	{
	  // lt == EDGE ok
	  li = i1;
	  lj = i2;
	  return ON_BOUNDARY;
	}
      case ON_BOUNDARY:
	{
	  // lt == VERTEX ok
	  li = ( i_e == 0 ) ? i1 : i2;
	  return ON_BOUNDARY;
	}
      case ON_UNBOUNDED_SIDE:
	{
	  // p lies on the line defined by the finite edge
	  return ON_UNBOUNDED_SIDE;
	}
      default:
	{
	  // cannot happen. only to avoid warning with eg++
	  return ON_UNBOUNDED_SIDE;
	}
      } 
    }// case ZERO
  }// switch o
  // end infinite facet
  // cannot happen. only to avoid warning with eg++
  CGAL_triangulation_assertion(false);
  return ON_UNBOUNDED_SIDE;
}

template < class GT, class Tds >
Bounded_side
Triangulation_3<GT,Tds>::
side_of_segment(const Point & p,
		const Point & p0, 
		const Point & p1,
		Locate_type & lt, int & i ) const
  // p0, p1 supposed to be different
  // p supposed to be collinear to p0, p1
  // returns :
  // ON_BOUNDED_SIDE if p lies strictly inside the edge
  // ON_BOUNDARY if p equals p0 or p1
  // ON_UNBOUNDED_SIDE if p lies strictly outside the edge
{
  CGAL_triangulation_precondition( ! equal(p0, p1) );
  CGAL_triangulation_precondition( collinear(p, p0, p1) );
      
  Comparison_result c = compare_x(p0,p1);
  Comparison_result c0;
  Comparison_result c1;

  if ( c == EQUAL ) {
    c = compare_y(p0,p1);
    if ( c == EQUAL ) {
      c0 = compare_z(p0,p);
      c1 = compare_z(p,p1);
    }
    else {
      c0 = compare_y(p0,p);
      c1 = compare_y(p,p1);
    }
  }
  else {
    c0 = compare_x(p0,p);
    c1 = compare_x(p,p1);
  }
      
  //      if ( (c0 == SMALLER) && (c1 == SMALLER) ) {
  if ( c0 == c1 ) {
    lt = EDGE;
    return ON_BOUNDED_SIDE;
  }
  if (c0 == EQUAL) {
    lt = VERTEX;
    i = 0;
    return ON_BOUNDARY;
  }
  if (c1 == EQUAL) {
    lt = VERTEX;
    i = 1;
    return ON_BOUNDARY;
  }
  lt = OUTSIDE_CONVEX_HULL;
  return ON_UNBOUNDED_SIDE;
}

template < class GT, class Tds >
Bounded_side
Triangulation_3<GT,Tds>::
side_of_edge(const Point & p,
	     Cell_handle c,
	     Locate_type & lt, int & li) const
  // supposes dimension 1 otherwise does not work for infinite edges
  // returns :
  // ON_BOUNDED_SIDE if p inside the edge 
  // (for an infinite edge this means that p lies in the half line
  // defined by the vertex)
  // ON_BOUNDARY if p equals one of the vertices
  // ON_UNBOUNDED_SIDE if p lies outside the edge
  // (for an infinite edge this means that p lies on the other half line)
  // lt has a meaning when ON_BOUNDED_SIDE and ON_BOUNDARY  
  // li refer to indices in the cell c 
{//side_of_edge
  CGAL_triangulation_precondition( dimension() == 1 );
  if ( ! is_infinite(c,0,1) ) 
    return side_of_segment(p, c->vertex(0)->point(), c->vertex(1)->point(),
			   lt, li);
  // else infinite edge
  int inf = c->index(infinite);
  if ( equal( p, c->vertex(1-inf)->point() ) ) {
    lt = VERTEX;
    li = 1-inf;
    return ON_BOUNDARY;
  }
  // does not work in dimension > 2
  Cell_handle n = c->neighbor(inf);
  int i_e = n->index(c);
  // we know that n is finite
  Vertex_handle
    v0 = n->vertex(0),
    v1 = n->vertex(1);
  Comparison_result c01 = compare_x(v0->point(), v1->point());
  Comparison_result cp;
  if ( c01 == EQUAL ) {
    c01 = compare_y(v0->point(),v1->point());
    if ( i_e == 0 ) {
      cp = compare_y( v1->point(), p );
    }
    else {
      cp = compare_y( p, v0->point() );
    }
  }
  else {
    if ( i_e == 0 ) 
      cp = compare_x( v1->point(), p );
    else 
      cp = compare_x( p, v0->point() );
  }
  if ( c01 == cp ) {
    // p lies on the same side of n as infinite
    lt = EDGE;
    return ON_BOUNDED_SIDE;
  }
  return ON_UNBOUNDED_SIDE;
}

template < class GT, class Tds >
inline
bool
Triangulation_3<GT,Tds>::
flip( Facet f )
  // returns false if the facet is not flippable
  // true other wise and
  // flips facet i of cell c
  // c will be replaced by one of the new cells
{
  return flip( f.first, f.second);
}

template < class GT, class Tds >
bool
Triangulation_3<GT,Tds>::
flip( Cell_handle c, int i )
{
  CGAL_triangulation_precondition( (dimension() == 3) && (0<=i) && (i<4) 
				   && (number_of_vertices() > 5) );

  Cell_handle n = c->neighbor(i);
  int in = n->index(c);
  if ( is_infinite( c ) || is_infinite( n ) ) return false;
  
  if ( i%2 == 1 ) {
    if ( orientation( c->vertex((i+1)&3)->point(),
		      c->vertex((i+2)&3)->point(),
		      n->vertex(in)->point(),
		      c->vertex(i)->point() )
	 != POSITIVE ) return false;
    if ( orientation( c->vertex((i+2)&3)->point(),
		      c->vertex((i+3)&3)->point(),
		      n->vertex(in)->point(),
		      c->vertex(i)->point() )
	 != POSITIVE ) return false;
    if ( orientation( c->vertex((i+3)&3)->point(),
		      c->vertex((i+1)&3)->point(),
		      n->vertex(in)->point(),
		      c->vertex(i)->point() )
	 != POSITIVE ) return false;
  }
  else {
    if ( orientation( c->vertex((i+2)&3)->point(),
		      c->vertex((i+1)&3)->point(),
		      n->vertex(in)->point(),
		      c->vertex(i)->point() )
	 != POSITIVE ) return false;
    if ( orientation( c->vertex((i+3)&3)->point(),
		      c->vertex((i+2)&3)->point(),
		      n->vertex(in)->point(),
		      c->vertex(i)->point() )
	 != POSITIVE ) return false;
    if ( orientation( c->vertex((i+1)&3)->point(),
		      c->vertex((i+3)&3)->point(),
		      n->vertex(in)->point(),
		      c->vertex(i)->point() )
	 != POSITIVE ) return false;
  }

  _tds.flip_flippable( &(*c), i);
  return true;
}

template < class GT, class Tds >
inline
void
Triangulation_3<GT,Tds>::
flip_flippable( Facet f )
{
  return flip_flippable( f.first, f.second);
}

template < class GT, class Tds >
void
Triangulation_3<GT,Tds>::
flip_flippable( Cell_handle c, int i )
{
  CGAL_triangulation_precondition( (dimension() == 3) && (0<=i) && (i<4) 
				   && (number_of_vertices() > 5) );
  CGAL_triangulation_precondition_code( Cell_handle n = c->neighbor(i); );
  CGAL_triangulation_precondition_code( int in = n->index(c); );
  CGAL_triangulation_precondition( ( ! is_infinite( c ) ) && 
				   ( ! is_infinite( n ) ) );
  
  if ( i%2 == 1 ) {
    CGAL_triangulation_precondition( orientation( c->vertex((i+1)&3)->point(),
						  c->vertex((i+2)&3)->point(),
						  n->vertex(in)->point(),
						  c->vertex(i)->point() )
				     == POSITIVE );
    CGAL_triangulation_precondition( orientation( c->vertex((i+2)&3)->point(),
						  c->vertex((i+3)&3)->point(),
						  n->vertex(in)->point(),
						  c->vertex(i)->point() )
				     == POSITIVE );
    CGAL_triangulation_precondition( orientation( c->vertex((i+3)&3)->point(),
						  c->vertex((i+1)&3)->point(),
						  n->vertex(in)->point(),
						  c->vertex(i)->point() )
				     == POSITIVE );
  }
  else {
    CGAL_triangulation_precondition( orientation( c->vertex((i+2)&3)->point(),
						  c->vertex((i+1)&3)->point(),
						  n->vertex(in)->point(),
						  c->vertex(i)->point() )
				     == POSITIVE );
    CGAL_triangulation_precondition( orientation( c->vertex((i+3)&3)->point(),
						  c->vertex((i+2)&3)->point(),
						  n->vertex(in)->point(),
						  c->vertex(i)->point() )
				     == POSITIVE );
    CGAL_triangulation_precondition( orientation( c->vertex((i+1)&3)->point(),
						  c->vertex((i+3)&3)->point(),
						  n->vertex(in)->point(),
						  c->vertex(i)->point() )
				     == POSITIVE );
  }
  
  _tds.flip_flippable( &(*c), i);
}

template < class GT, class Tds >
inline
bool
Triangulation_3<GT,Tds>::
flip( Edge e )
  // returns false if the edge is not flippable
  // true otherwise and
  // flips edge i,j of cell c
  // c will be deleted
{
  return flip( e.first, e.second, e.third );
}

template < class GT, class Tds >
bool
Triangulation_3<GT,Tds>::
flip( Cell_handle c, int i, int j )
  // flips edge i,j of cell c
{
  CGAL_triangulation_precondition( (dimension() == 3) 
				   && (0<=i) && (i<4) 
				   && (0<=j) && (j<4)
				   && ( i != j )
				   && (number_of_vertices() > 5) );

  // checks that degree 3 and not on the convex hull
  int degree = 0;
  Cell_circulator ccir = incident_cells(c,i,j);
  Cell_circulator cdone = ccir;
  do {
    if ( is_infinite(&(*ccir)) ) return false;
    ++degree;
    ++ccir;
  } while ( ccir != cdone );

  if ( degree != 3 ) return false;

  // checks that future tetrahedra are well oriented
  Cell_handle n = c->neighbor( next_around_edge(i,j) );
  int in = n->index( c->vertex(i) );
  int jn = n->index( c->vertex(j) );
  if ( orientation( c->vertex(next_around_edge(i,j))->point(),
		    c->vertex(next_around_edge(j,i))->point(),
		    n->vertex(next_around_edge(jn,in))->point(),
		    c->vertex(j)->point() )
       != POSITIVE ) return false;
  if ( orientation( c->vertex(i)->point(),
		    c->vertex(next_around_edge(j,i))->point(),
		    n->vertex(next_around_edge(jn,in))->point(),
		    c->vertex(next_around_edge(i,j))->point() )
       != POSITIVE ) return false;

  _tds.flip_flippable( &(*c), i, j );
  return true;
}

template < class GT, class Tds >
inline
void
Triangulation_3<GT,Tds>::
flip_flippable( Edge e )
{
  return flip_flippable( e.first, e.second, e.third );
}

template < class GT, class Tds >
void
Triangulation_3<GT,Tds>::
flip_flippable( Cell_handle c, int i, int j )
  // flips edge i,j of cell c
{
#if !defined CGAL_TRIANGULATION_NO_PRECONDITIONS && \
    !defined CGAL_NO_PRECONDITIONS && !defined NDEBUG
  CGAL_triangulation_precondition( (dimension() == 3) 
				   && (0<=i) && (i<4) 
				   && (0<=j) && (j<4)
				   && ( i != j )
				   && (number_of_vertices() > 5) );
  int degree = 0;
  Cell_circulator ccir = incident_cells(c,i,j);
  Cell_circulator cdone = ccir;
  do {
    CGAL_triangulation_precondition( ! is_infinite(&(*ccir)) );
    ++degree;
    ++ccir;
  } while ( ccir != cdone );
  CGAL_triangulation_precondition( degree == 3 );

  Cell_handle n = c->neighbor( next_around_edge(i, j) );
  int in = n->index( c->vertex(i) );
  int jn = n->index( c->vertex(j) );
  CGAL_triangulation_precondition
    ( orientation( c->vertex(next_around_edge(i,j))->point(),
		   c->vertex(next_around_edge(j,i))->point(),
		   n->vertex(next_around_edge(jn,in))->point(),
		   c->vertex(j)->point() ) == POSITIVE );
  CGAL_triangulation_precondition
    ( orientation( c->vertex(i)->point(),
		   c->vertex(next_around_edge(j,i))->point(),
		   n->vertex(next_around_edge(jn,in))->point(),
		   c->vertex(next_around_edge(i,j))->point() ) == POSITIVE );
#endif
  _tds.flip_flippable( &(*c), i, j );
}

template < class GT, class Tds >
Triangulation_3<GT,Tds>::Vertex_handle
Triangulation_3<GT,Tds>::
insert(const Point & p, Cell_handle start, Vertex_handle v)
{
  Locate_type lt;
  int li, lj;
  Cell_handle c = locate( p, lt, li, lj, start);
  switch (lt) {
  case VERTEX:
    return c->vertex(li);
  case EDGE:
    return insert_in_edge(p, c, li, lj, v);
  case FACET:
    return insert_in_facet(p, c, li, v);
  case CELL:
    return insert_in_cell(p, c, v);
  case OUTSIDE_CONVEX_HULL:
    return insert_outside_convex_hull(p, c, v);
  case OUTSIDE_AFFINE_HULL:
  default:
    return insert_outside_affine_hull(p, v);
  }
}

template < class GT, class Tds >
Triangulation_3<GT,Tds>::Vertex_handle
Triangulation_3<GT,Tds>::
insert_in_cell(const Point & p, Cell_handle c, Vertex_handle v)
{
  CGAL_triangulation_precondition( dimension() == 3 );
  CGAL_triangulation_precondition_code
    ( Locate_type lt;
      int i; int j; );
  CGAL_triangulation_precondition
    ( side_of_tetrahedron( p, 
			   c->vertex(0)->point(),
			   c->vertex(1)->point(),
			   c->vertex(2)->point(),
			   c->vertex(3)->point(),
			   lt,i,j ) == ON_BOUNDED_SIDE );

    v = (Vertex*)_tds.insert_in_cell( &(*v), &(*c) );
    v->set_point(p);
    return v;
}

template < class GT, class Tds >
inline
Triangulation_3<GT,Tds>::Vertex_handle
Triangulation_3<GT,Tds>::
insert_in_facet(const Point & p, Cell_handle c, int i, Vertex_handle v)
{
  CGAL_triangulation_precondition( dimension() == 2 || dimension() == 3);
  CGAL_triangulation_precondition( (dimension() == 2 && i == 3)
	                        || (dimension() == 3 && i >= 0 && i <= 3) );
  CGAL_triangulation_precondition_code
    ( Locate_type lt;
      int li; int lj; );
  CGAL_triangulation_precondition
    ( orientation( p, 
		   c->vertex((i+1)&3)->point(),
		   c->vertex((i+2)&3)->point(),
		   c->vertex((i+3)&3)->point() ) == COPLANAR
      && 
      side_of_triangle( p, 
			c->vertex((i+1)&3)->point(),
			c->vertex((i+2)&3)->point(),
			c->vertex((i+3)&3)->point(),
			lt, li, lj) == ON_BOUNDED_SIDE );

    v = (Vertex*) _tds.insert_in_facet( &(*v), &(*c), i);
    v->set_point(p);
    return v;
}

template < class GT, class Tds >
Triangulation_3<GT,Tds>::Vertex_handle
Triangulation_3<GT,Tds>::
insert_in_edge(const Point & p, Cell_handle c, int i, int j, Vertex_handle v)
{
  CGAL_triangulation_precondition( i != j );
  CGAL_triangulation_precondition( dimension() >= 1 && dimension() <= 3 );
  CGAL_triangulation_precondition( i >= 0 && i <= dimension() 
				   && j >= 0 && j <= dimension() );
  CGAL_triangulation_precondition_code( Locate_type lt; int li; );
  switch ( dimension() ) {
  case 3:
  case 2:
    {
      CGAL_triangulation_precondition( ! is_infinite(c, i, j) );
      CGAL_triangulation_precondition( collinear( c->vertex(i)->point(),
						  p,
						  c->vertex(j)->point() )
				       &&
				       side_of_segment( p,
							c->vertex(i)->point(),
							c->vertex(j)->point(),
							lt, li ) 
				       == ON_BOUNDED_SIDE );
      break;
    }
  case 1:
    {
      CGAL_triangulation_precondition( side_of_edge(p, c, lt, li)
				       == ON_BOUNDED_SIDE );
      break;
    }
  }

  v = (Vertex*) _tds.insert_in_edge( &(*v), &(*c), i, j);
  v->set_point(p);
  return v;
}

template < class GT, class Tds >
Triangulation_3<GT,Tds>::Vertex_handle
Triangulation_3<GT,Tds>::
insert_outside_convex_hull(const Point & p, Cell_handle c, Vertex_handle v)
  // c is an infinite cell containing p
  // p is strictly outside the convex hull
  // dimension 0 not allowed, use outside-affine-hull
{
  CGAL_triangulation_precondition( dimension() > 0 );
  CGAL_triangulation_precondition( c->has_vertex(infinite) );
  // the precondition that p is in c is tested in each of the
  // insertion methods called from this method 
  switch ( dimension() ) {
  case 1:
    {
      // 	// p lies in the infinite edge neighboring c 
      // 	// on the other side of li
      // 	return insert_in_edge(p,c->neighbor(1-li),0,1);
      return insert_in_edge(p,c,0,1,v);
    }
  case 2:
    {
      set_number_of_vertices(number_of_vertices()+1);

      Conflict_tester_outside_convex_hull_2 tester(p, this);
      Vertex_handle v = (Vertex *) _tds.insert_conflict(NULL, &(*c), tester);
      v->set_point(p);
      
      return v;
    }
  case 3:
  default:
    {
      set_number_of_vertices(number_of_vertices()+1);

      Conflict_tester_outside_convex_hull_3 tester(p, this);
      Vertex_handle v = (Vertex *) _tds.insert_conflict(NULL, &(*c), tester);
      v->set_point(p);
      return v;
    }
  }
}

template < class GT, class Tds >
Triangulation_3<GT,Tds>::Vertex_handle
Triangulation_3<GT,Tds>::
insert_outside_affine_hull(const Point & p, Vertex_handle v)
{
  CGAL_triangulation_precondition( dimension() < 3 );
  bool reorient;
  switch ( dimension() ) {
  case 1:
    {
      Cell_handle c = infinite_cell();
      Cell_handle n = c->neighbor(c->index(infinite_vertex()));
      Orientation o = coplanar_orientation(n->vertex(0)->point(),
					   n->vertex(1)->point(), p);
      CGAL_triangulation_precondition ( o != COLLINEAR );
      reorient = o == NEGATIVE;
      break;
    }
  case 2:
    {
      Cell_handle c = infinite_cell();
      Cell_handle n = c->neighbor(c->index(infinite_vertex()));
      Orientation o = orientation( n->vertex(0)->point(),
			           n->vertex(1)->point(),
			           n->vertex(2)->point(), p );
      CGAL_triangulation_precondition ( o != COPLANAR );
      reorient = o == NEGATIVE;
      break;
    }
  default:
    reorient = false;
  }

  v = (Vertex*) _tds.insert_increase_dimension( &(*v), 
						&(*infinite_vertex()), 
						reorient);
  v->set_point(p);
  return v;
}

template < class GT, class Tds >
void
Triangulation_3<GT,Tds>::
incident_cells(Vertex_handle v, 
	       std::set<Cell*> & cells,
	       Cell_handle c,
	       int dummy_for_windows) const
{
  bool WARNING_THIS_FUNCTION_IS_DEPRECATED;
  CGAL_triangulation_precondition( &(*v) != NULL );
  CGAL_triangulation_expensive_precondition( _tds.is_vertex(&(*v)) );

  if ( dimension() < 3 )
      return;

  if ( &(*c) == NULL )
    c = v->cell();
  else {
    CGAL_triangulation_precondition( c->has_vertex(v) );
  }
  if ( cells.find( &(*c) ) != cells.end() )
    return; // c was already found

  cells.insert( &(*c) );
      
  for ( int j=0; j<4; j++ )
    if ( j != c->index(v) )
      incident_cells( v, cells, c->neighbor(j), dummy_for_windows);
}

template < class GT, class Tds >
void
Triangulation_3<GT,Tds>::
incident_cells(Vertex_handle v, 
	       std::set<Cell_handle> & cells,
	       Cell_handle c ) const
{
  CGAL_triangulation_precondition( &(*v) != NULL );
  CGAL_triangulation_expensive_precondition( _tds.is_vertex(&(*v)) );

  if ( dimension() < 3 )
      return;

  if ( &(*c) == NULL )
    c = v->cell();
  else
    CGAL_triangulation_precondition( c->has_vertex(v) );
  if ( cells.find( c ) != cells.end() )
    return; // c was already found

  cells.insert( c );
      
  for ( int j=0; j<4; j++ )
    if ( j != c->index(v) )
      incident_cells( v, cells, c->neighbor(j) );
}

template < class GT, class Tds >
void
Triangulation_3<GT,Tds>::
incident_vertices(Vertex_handle v, 
		  std::set<Vertex*> & vertices,
		  Cell_handle c,
		  int dummy_for_windows) const
{
  bool WARNING_THIS_FUNCTION_IS_DEPRECATED;
  CGAL_triangulation_precondition( &(*v) != NULL );
  CGAL_triangulation_expensive_precondition( _tds.is_vertex(&(*v)) );
      
  if ( number_of_vertices() < 2 )
      return;

  if ( &(*c) == NULL )
    c = v->cell();
  else 
    CGAL_triangulation_precondition( c->has_vertex(v) );

  std::set<Cell*> cells;
  util_incident_vertices(v, vertices, cells, c, dummy_for_windows);
}

template < class GT, class Tds >
void
Triangulation_3<GT,Tds>::
incident_vertices(Vertex_handle v, 
		  std::set<Vertex_handle> & vertices,
		  Cell_handle c ) const
{
  CGAL_triangulation_precondition( &(*v) != NULL );
  CGAL_triangulation_expensive_precondition( _tds.is_vertex(&(*v)) );
      
  if ( number_of_vertices() < 2 )
      return;

  if ( &(*c) == NULL )
    c = v->cell();
  else 
    CGAL_triangulation_precondition( c->has_vertex(v) );

  std::set<Cell_handle> cells;
  util_incident_vertices(v, vertices, cells, c);
}

template < class GT, class Tds >
void
Triangulation_3<GT,Tds>::
util_incident_vertices(Vertex_handle v, 
		       std::set<Vertex*> & vertices,
		       std::set<Cell*> & cells,
		       Cell_handle c,
		       int dummy_for_windows) const
{
  bool WARNING_THIS_FUNCTION_IS_DEPRECATED;
  if ( cells.find( &(*c) ) != cells.end() )
    return; // c was already visited
  cells.insert( &(*c) );

  int d = dimension();
  for (int j=0; j <= d; j++ )
    if ( j != c->index(v) ) {
      if ( vertices.find( &(*(c->vertex(j))) ) == vertices.end() )
	vertices.insert( &(*(c->vertex(j))) );
      util_incident_vertices( v, vertices, cells, c->neighbor(j), 
			      dummy_for_windows);
    }
}

template < class GT, class Tds >
void
Triangulation_3<GT,Tds>::
util_incident_vertices(Vertex_handle v, 
		       std::set<Vertex_handle> & vertices,
		       std::set<Cell_handle> & cells,
		       Cell_handle c ) const
{
  if ( cells.find( c ) != cells.end() )
    return; // c was already visited

  cells.insert( c );

  int d = dimension();
  for (int j=0; j <= d; j++ )
    if ( j != c->index(v) ) {
      if ( vertices.find( c->vertex(j) ) == vertices.end() )
	vertices.insert( c->vertex(j) );
      util_incident_vertices( v, vertices, cells, c->neighbor(j) );
    }
}

template < class GT, class Tds >
bool
Triangulation_3<GT,Tds>::
is_valid(bool verbose, int level) const
{
  if ( ! _tds.is_valid(verbose,level) ) {
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
      for ( it = finite_cells_begin(); it != cells_end(); ++it )
	is_valid_finite((*it).handle(),verbose,level);
      break;
    }
  case 2:
    {
      Facet_iterator it;
      for ( it = finite_facets_begin(); it != facets_end(); ++it )
	is_valid_finite((*it).first,verbose,level);
      break;
    }
  case 1:
    {
      Edge_iterator it;
      for ( it = finite_edges_begin(); it != edges_end(); ++it )
	is_valid_finite((*it).first,verbose,level);
      break;
    }
  }
  if (verbose)
      std::cerr << "valid triangulation" << std::endl;
  return true;
}

template < class GT, class Tds >
bool
Triangulation_3<GT,Tds>::
is_valid(Cell_handle c, bool verbose, int level) const
{
  if ( ! (&(*c))->is_valid(dimension(),verbose,level) ) {
    if (verbose) { 
      std::cerr << "combinatorially invalid cell";
      for (int i=0; i <= dimension(); i++ )
	std::cerr << c->vertex(i)->point() << ", ";
      std::cerr << std::endl;
    }
    CGAL_triangulation_assertion(false);
    return false;
  }
  if ( ! is_infinite(c) )
    is_valid_finite(c, verbose, level);
  if (verbose)
      std::cerr << "geometrically valid cell" << std::endl;
  return true;
}


template < class GT, class Tds >
bool
Triangulation_3<GT,Tds>::
is_valid_finite(Cell_handle c, bool verbose, int) const
{
  switch ( dimension() ) {
  case 3:
    {
      if ( orientation(c->vertex(0)->point(),
		       c->vertex(1)->point(),
		       c->vertex(2)->point(),
		       c->vertex(3)->point()) != POSITIVE ) {
	if (verbose)
	    std::cerr << "badly oriented cell " 
		      << c->vertex(0)->point() << ", " 
		      << c->vertex(1)->point() << ", " 
		      << c->vertex(2)->point() << ", " 
		      << c->vertex(3)->point() << std::endl; 
	CGAL_triangulation_assertion(false);
	return false;
      }
      break;
    }
  case 2:
    {
	if (coplanar_orientation(c->vertex(0)->point(),
	                         c->vertex(1)->point(),
	                         c->vertex(2)->point()) != POSITIVE) {
	  if (verbose)
	      std::cerr << "badly oriented face "
		        << c->vertex(0)->point() << ", " 
		        << c->vertex(1)->point() << ", " 
		        << c->vertex(2)->point() << std::endl;
	  CGAL_triangulation_assertion(false);
	  return false;
	}
      break;
    }
  case 1:
    {
      const Point & p0 = c->vertex(0)->point();
      const Point & p1 = c->vertex(1)->point();
	    
      if ( ! is_infinite ( c->neighbor(0)->vertex(c->neighbor(0)->index(c)) ) )
      {
	const Point & n0 =
	    c->neighbor(0)->vertex(c->neighbor(0)->index(c))->point();  
	if ( ( compare_x( p0, p1 ) != compare_x( p1, n0 ) )
	  || ( compare_y( p0, p1 ) != compare_y( p1, n0 ) )
	  || ( compare_z( p0, p1 ) != compare_z( p1, n0 ) ) ) {
	  if (verbose)
	      std::cerr << "badly oriented edge "
		        << p0 << ", " << p1 << std::endl
		        << "with neighbor 0"
		        << c->neighbor(0)->vertex(1-c->neighbor(0)->index(c))
			                 ->point() 
		        << ", " << n0 << std::endl;
	  CGAL_triangulation_assertion(false);
	  return false;
	}
      }
      if ( ! is_infinite ( c->neighbor(1)->vertex(c->neighbor(1)->index(c)) ) )
      {
	const Point & n1 = 
	  c->neighbor(1)->vertex(c->neighbor(1)->index(c))->point();
	if ( ( compare_x( p1, p0 ) != compare_x( p0, n1 ) )
	  || ( compare_y( p1, p0 ) != compare_y( p0, n1 ) )
	  || ( compare_z( p1, p0 ) != compare_z( p0, n1 ) ) ) {
	  if (verbose)
	      std::cerr << "badly oriented edge "
		        << p0 << ", " << p1 << std::endl
		        << "with neighbor 1"
		        << c->neighbor(1)->vertex(1-c->neighbor(1)->index(c))
			                 ->point() 
		        << ", " << n1 << std::endl;
	  CGAL_triangulation_assertion(false);
	  return false;
	}
      }
      break;
    }
  }
  return true;
}

CGAL_END_NAMESPACE

#endif // CGAL_TRIANGULATION_3_H
