// Copyright (c) 1997  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Tran Kai Frank DA <Frank.Da@sophia.inria.fr>
//                 Andreas Fabri <Andreas.Fabri@geometryfactory.com>

#ifndef CGAL_ALPHA_SHAPE_3_H
#define CGAL_ALPHA_SHAPE_3_H

#include <CGAL/basic.h>

#include <cassert>
#include <set>
#include <map>
#include <vector>
#include <algorithm>
#include <utility>
#include <iostream>

#include <CGAL/utility.h>
#include <CGAL/Unique_hash_map.h>
#include <CGAL/IO/Geomview_stream.h>  // TBC

//-------------------------------------------------------------------
CGAL_BEGIN_NAMESPACE
//-------------------------------------------------------------------

template < class Dt >
class Alpha_shape_3 : public Dt
{
  // DEFINITION The class Alpha_shape_3<Dt> represents the family
  // of alpha-shapes of points in a plane for all positive alpha. It
  // maintains the underlying Delaunay tetrahedralization which represents
  // connectivity and order among its simplices. Each k-dimensional simplex of
  // the Delaunay tetrahedralization is associated with an interval that
  // specifies for which values of alpha the simplex belongs to the alpha-shape
  // (sorted linear arrays resp. multimaps or interval trees). There are links
  // between the intervals and the k-dimensional simplices of the Delaunay
  // tetrahedralization (multimaps resp. hashtables).
  //

  //------------------------- TYPES ------------------------------------

public:

  typedef typename Dt::Geom_traits Gt;
  typedef typename Dt::Triangulation_data_structure Tds;

  typedef typename Gt::FT Coord_type;

  typedef typename Gt::Point_3 Point;
  typedef typename Gt::Segment_3 Segment;
  //   typedef typename Gt::Triangle_3 Triangle;
  //   typedef typename Gt::Tetrahedron_3 Tetrahedron;

  typedef typename Dt::Cell_handle Cell_handle;
  typedef typename Dt::Vertex_handle Vertex_handle;
  typedef typename Dt::Cell Cell;
  typedef typename Dt::Vertex Vertex;
  typedef typename Dt::Facet Facet;
  typedef typename Dt::Edge Edge;

  typedef typename Dt::Cell_circulator  Cell_circulator;
  typedef typename Dt::Facet_circulator Facet_circulator;

  typedef typename Dt::Cell_iterator   Cell_iterator;
  typedef typename Dt::Facet_iterator  Facet_iterator;
  typedef typename Dt::Edge_iterator   Edge_iterator;
  typedef typename Dt::Vertex_iterator Vertex_iterator;

  typedef typename Dt::Finite_cells_iterator    Finite_cells_iterator;
  typedef typename Dt::Finite_facets_iterator   Finite_facets_iterator;
  typedef typename Dt::Finite_edges_iterator    Finite_edges_iterator;
  typedef typename Dt::Finite_vertices_iterator Finite_vertices_iterator;

  typedef typename Dt::Locate_type Locate_type;

private:

  typedef std::multimap< Coord_type, Cell_handle > Interval_cell_map;
  typedef typename Interval_cell_map::value_type   Interval_cell;



  //   typedef Cell_handle const const_void;
  //   typedef std::pair< const_void, int > const_Facet;
  //   typedef std::pair< const_void, int > const_Vertex;




  typedef std::vector< Coord_type > Alpha_spectrum;
  
  typedef Unique_hash_map<Cell_handle, bool > Marked_cell_set;

public:

  //the following eight typedef were private, 
  // but operator<<(ostream) needs them
  typedef Triple<Coord_type, Coord_type, Coord_type> Interval3;
  typedef std::pair< Interval3, Edge > Interval_edge;
  typedef std::multimap< Interval3, Edge > Interval_edge_map;
  typedef std::multimap< Interval3, Facet >        Interval_facet_map;
  typedef typename Interval_facet_map::value_type  Interval_facet;
  typedef std::pair< Coord_type, Coord_type > Interval2;
  typedef std::multimap< Interval2, Vertex_handle > Interval_vertex_map; 
  typedef typename Interval_vertex_map::value_type  Interval_vertex;

  typedef typename Alpha_spectrum::const_iterator Alpha_iterator;
  // An iterator that allow to traverse the sorted sequence of
  // different alpha-values. The iterator is bidirectional and
  // non-mutable. Its value-type is Coord_type

  enum Classification_type {EXTERIOR, SINGULAR, REGULAR, INTERIOR};
  // Distinguishes the different cases for classifying a
  // k-dimensional cell of the underlying Delaunay tetrahedralization of
  // the alpha-shape.
  // 
  // `EXTERIOR' if the cell does not belong to the alpha-complex.
  // 
  // `SINGULAR' if the cell belongs to the boundary of the
  // alpha-shape, but is not incident to any higher-dimensional
  // cell of the alpha-complex
  // 
  // `REGULAR' if cell belongs to the boundary of the alpha-shape
  // and is incident to a higher-dimensional cell of the
  // alpha-complex
  // 
  // `INTERIOR' if the cell belongs to the alpha-complex, but does
  // not belong to the boundary of the alpha-shape

  enum Mode {GENERAL, REGULARIZED};
  // In general, an alpha shape is a non-connected, mixed-dimension
  // polygon. Its regularized version is formed by the set of
  // regular facets and their vertices

  typedef typename std::list< Vertex_handle >::iterator 
  Alpha_shape_vertices_iterator;
  typedef typename std::list< Facet >::iterator
  Alpha_shape_facets_iterator;
  
  class Out_alpha_shape_test
  //test if a cell belong to the alphashape
  {
    const Alpha_shape_3 * _as;
  public:
    Out_alpha_shape_test() {}
    Out_alpha_shape_test(const Alpha_shape_3 * as) {_as = as;}
    bool operator() ( const Finite_cells_iterator& fci) const {
      return _as->classify(fci) == EXTERIOR ;
    }
  };

  typedef Filter_iterator< Finite_cells_iterator, Out_alpha_shape_test>
  Alpha_shape_cells_iterator;

public:  // should be private ? --> operator should be wrappers

  // only finite facets and simplices are inserted into the maps 
  Interval_cell_map _interval_cell_map;
  Interval_facet_map _interval_facet_map;
  Interval_edge_map _interval_edge_map;
  Interval_vertex_map _interval_vertex_map;

  Alpha_spectrum _alpha_spectrum;

  Coord_type _alpha;
  Mode _mode;
 
  // should be constants
  Coord_type Infinity;
  Coord_type UNDEFINED;

  mutable std::list< Vertex_handle > Alpha_shape_vertices_list;
  mutable std::list< Facet > Alpha_shape_facets_list;

  mutable bool use_vertex_cache;
  mutable bool use_facet_cache;


  //------------------------- CONSTRUCTORS ------------------------------
public:
  // Introduces an empty alpha-shape `A' for a positive
  // alpha-value `alpha'. Precondition: `alpha' >= 0.
  Alpha_shape_3(Coord_type alpha = 0, 
		Mode m = REGULARIZED)
    : _alpha(alpha), _mode(m), Infinity(-1), UNDEFINED(-2), 
    use_vertex_cache(false), use_facet_cache(false)
    {}

  Alpha_shape_3(Dt& dt, Coord_type alpha = 0, Mode m = REGULARIZED)
    :_alpha(alpha), _mode(m), Infinity(-1), UNDEFINED(-2), 
    use_vertex_cache(false), use_facet_cache(false)
    {
      Dt::swap(dt);
       if (dimension() == 3)
	{
	  initialize_interval_cell_map();
 	  initialize_interval_facet_map();
	  initialize_interval_edge_map();
	  initialize_interval_vertex_map();
	  
	  // merge the two maps
	  initialize_alpha_spectrum();
	}
    }
 
  // Introduces an alpha-shape `A' for a positive alpha-value
  // `alpha' that is initialized with the points in the range
  // from first to last

  template < class InputIterator >  
  Alpha_shape_3(const InputIterator& first,  
		const InputIterator& last,  
		const Coord_type& alpha = 0,
		Mode m = REGULARIZED)
    : _alpha(alpha), _mode(m), Infinity(-1), UNDEFINED(-2), 
    use_vertex_cache(false), use_facet_cache(false)
    {
      Dt::insert(first, last);
      if (dimension() == 3)
	{
	  initialize_interval_cell_map();
 	  initialize_interval_facet_map();
	  initialize_interval_edge_map();
	  initialize_interval_vertex_map();
	  
	  // merge the two maps
	  initialize_alpha_spectrum();
	}
    }
 
public:

  //----------------------- OPERATIONS ---------------------------------


  template < class InputIterator >  
  int make_alpha_shape(const InputIterator& first, 
		       const InputIterator& last)
    {
      clear();
   
      int n = Dt::insert(first, last);

#ifdef DEBUG
      std::cout << "Triangulation computed" << std::endl;
#endif
      if (dimension() == 3)
	{
	  initialize_interval_cell_map();
	  initialize_interval_facet_map();
	  initialize_interval_edge_map();
 	  initialize_interval_vertex_map();

	  // merge the two maps
	  initialize_alpha_spectrum();
	}
      return n;
    }

  // Introduces an alpha-shape `A' for a positive alpha-value
  // `alpha' that is initialized with the points in the range
  // from first to last

private :

  //--------------------- INITIALIZATION OF PRIVATE MEMBERS -----------
  
  void initialize_interval_cell_map();

  void initialize_interval_facet_map();

  void initialize_interval_edge_map() {} //disabled here
  
  void initialize_interval_vertex_map();
 
  void initialize_alpha_spectrum();
  
  //---------------------------------------------------------------------

public:

  void clear()
    {
      // clears the structure
      Dt::clear();

      _interval_cell_map.clear();
      _interval_facet_map.clear();
      _interval_edge_map.clear();
      _interval_vertex_map.clear();

      _alpha_spectrum.clear();

      Alpha_shape_vertices_list.clear();

      Alpha_shape_facets_list.clear();

      set_alpha(0); 
      set_mode(REGULARIZED);
      use_vertex_cache = false;
      use_facet_cache = false;

    }

  //----------------------- PRIVATE MEMBERS --------------------

private:
   
  struct Less
  {
    bool operator()(const Interval_facet& ie, 
		    const Coord_type& alpha)
      { return ie.first.first < alpha; }

    bool operator()(const Coord_type& alpha, 
		    const Interval_facet& ie)
      { return alpha < ie.first.first; }
  };
  
  
  //----------------------- ACCESS TO PRIVATE MEMBERS -----------------

private:

  Coord_type find_interval(Cell_handle s) const
    // check whether it is faster to compute the 
    // radius directly instead of looking it up
    {
      return s->get_alpha();
    }
  
  Interval3 find_interval(const Facet& f) const
    {
      return (f.first)->get_facet_ranges(f.second);    
    }

  Interval3 find_interval(const Edge& e) const
    {
      return find_interval(e.first, e.second, e.third);
    }
  
  Interval3 find_interval(const Cell_handle& s, int i, int j) const;

  Interval2 find_interval(const Vertex_handle& v) const
    {
      return v->get_range();
    }

  //---------------------------------------------------------------------

public:

  Coord_type set_alpha(const Coord_type& alpha)
    // Sets the alpha-value to `alpha'. Precondition: `alpha' >= 0.
    // Returns the previous alpha
    {
      Coord_type previous_alpha = _alpha;
      _alpha = alpha;
      use_vertex_cache = false;
      use_facet_cache = false;
      return previous_alpha;
    }

  const Coord_type&  get_alpha() const
    // Returns the current alpha-value.
    {
      return _alpha;
    }
  

  const Coord_type&  get_nth_alpha(int n) const
    // Returns the n-th alpha-value.
    // n < size()
    {
      if (! _alpha_spectrum.empty())
	return _alpha_spectrum[n];
      else
	return UNDEFINED;
    }
  
  int number_of_alphas() const
    // Returns the number of not necessary different alpha-values
    {
      return _alpha_spectrum.size();
    }

  
  //---------------------------------------------------------------------

private:

  // the dynamic version is not yet implemented
  // desactivate the tetrahedralization member functions
  Vertex_handle insert(const Point& p) {}
  // Inserts point `p' in the alpha shape and returns the
  // corresponding vertex of the underlying Delaunay tetrahedralization.
  // If point `p' coincides with an already existing vertex, this
  // vertex is returned and the alpha shape remains unchanged.
  // Otherwise, the vertex is inserted in the underlying Delaunay
  // tetrahedralization and the associated intervals are updated.

  void remove(Vertex_handle v) {}
  // Removes the vertex from the underlying Delaunay tetrahedralization.
  // The created hole is retriangulated and the associated intervals
  // are updated.

  //---------------------------------------------------------------------

public:

  Mode set_mode(Mode mode = REGULARIZED )
    // Sets `A' to its general or regularized version. Returns the
    // previous mode.
    {
      Mode previous_mode = _mode;
      _mode = mode;
      return previous_mode;
    }

  Mode get_mode() const
    // Returns whether `A' is general or regularized.
    {
      return _mode;
    }

  //---------------------------------------------------------------------
private:

  void
  update_alpha_shape_vertex_list() const;

  //---------------------------------------------------------------------

  void
  update_alpha_shape_facet_list() const; 

  //---------------------------------------------------------------------
public:

  Alpha_shape_vertices_iterator alpha_shape_vertices_begin() const
  {
    if(!use_vertex_cache){
      update_alpha_shape_vertex_list();
    }
      return Alpha_shape_vertices_list.begin();
    }

  Alpha_shape_vertices_iterator Alpha_shape_vertices_begin() const
    {
      return alpha_shape_vertices_begin();
    }
  //---------------------------------------------------------------------

  Alpha_shape_vertices_iterator alpha_shape_vertices_end() const
    {
      return Alpha_shape_vertices_list.end();
    }

  Alpha_shape_vertices_iterator Alpha_shape_vertices_end() const
    {
      return alpha_shape_vertices_end();
    }

  //---------------------------------------------------------------------

  Alpha_shape_facets_iterator alpha_shape_facets_begin() const
    {
      if(! use_facet_cache){
	update_alpha_shape_facet_list();
      }
      return Alpha_shape_facets_list.begin();
    }

  Alpha_shape_facets_iterator Alpha_shape_facets_begin() const
    {
      return alpha_shape_facets_begin();
    }

  //---------------------------------------------------------------------

  Alpha_shape_facets_iterator alpha_shape_facets_end() const
    {
      return Alpha_shape_facets_list.end();
    }

  Alpha_shape_facets_iterator Alpha_shape_facets_end() const
    {
      return alpha_shape_facets_end();
    }

  Alpha_shape_cells_iterator alpha_shape_cells_begin() const 
    {
      return filter_iterator(finite_cells_end(),
			     Out_alpha_shape_test(this),
			     finite_cells_begin());
    }
  
  Alpha_shape_cells_iterator alpha_shape_cells_end() const
    {
      return filter_iterator(finite_cells_end(),
			     Out_alpha_shape_test(this));
    }


public: 
  
  // Traversal of the alpha-Values
  // 
  // The alpha shape class defines an iterator that allows to
  // visit the sorted sequence of alpha-values. This iterator is
  // non-mutable and bidirectional. Its value type is Coord_type.

  Alpha_iterator alpha_begin() const
    // Returns an iterator that allows to traverse the sorted sequence
    // of alpha-values of `A'.
    {
      return _alpha_spectrum.begin(); 
    }

  Alpha_iterator alpha_end() const
    // Returns the corresponding past-the-end iterator.
    {
      return _alpha_spectrum.end(); 
    }

  Alpha_iterator alpha_find(const Coord_type& alpha) const
    // Returns an iterator pointing to an element with alpha-value
    // `alpha', or the corresponding past-the-end iterator if such an
    // element is not found.
    {
      return find(_alpha_spectrum.begin(),
		  _alpha_spectrum.end(),
		  alpha);
    }

  Alpha_iterator alpha_lower_bound(const Coord_type& alpha) const
    // Returns an iterator pointing to the first element with
    // alpha-value not less than `alpha'.
    {
      return std::lower_bound(_alpha_spectrum.begin(),
			      _alpha_spectrum.end(),
			      alpha);
    }

  Alpha_iterator alpha_upper_bound(const Coord_type& alpha) const
    // Returns an iterator pointing to the first element with
    // alpha-value greater than `alpha'.
    {
      return std::upper_bound(_alpha_spectrum.begin(),
			      _alpha_spectrum.end(),
			      alpha);
    }

  //--------------------- PREDICATES -----------------------------------

  // the classification predicates take 
  //      amortized const time if STL_HASH_TABLES
  //      O(log #alpha_shape ) otherwise

  Classification_type  classify(const Point& p) const
    {
      return classify(p, get_alpha());
    }

  
  Classification_type  classify(const Point& p,   
				const Coord_type& alpha) const
    // Classifies a point `p' with respect to `A'.
    {
      Locate_type type;
      int i, j;
      Cell_handle pCell = locate(p, type, i, j);
      switch (type)
	{
	case VERTEX            : return classify(pCell->vertex(i), alpha);
	case EDGE              : return classify(pCell, i, j, alpha);
	case FACET             : return classify(pCell, i, alpha);
	case CELL              : return classify(pCell, alpha);
	case OUTSIDE_CONVEX_HULL : return EXTERIOR;
	case OUTSIDE_AFFINE_HULL : return EXTERIOR;
	default                : return EXTERIOR;
	};
    }
 
  //---------------------------------------------------------------------

  Classification_type  classify(const Cell_handle& s) const
    // Classifies the cell `f' of the underlying Delaunay
    // tetrahedralization with respect to `A'.
    {
      return classify(s, get_alpha());
    }
  
  Classification_type  classify(const Cell_handle& s, 
				const Coord_type& alpha) const
    // Classifies the cell `f' of the underlying Delaunay
    // tetrahedralization with respect to `A'.
    // we consider open spheres :
    // s->radius == alpha => f exterior
    // problem the operator [] is non-const
    {
      if (is_infinite(s)) return EXTERIOR;
    
      // the version that computes the squared radius seems to be 
      // much faster
    
      return (s->get_alpha() < alpha) ? 
	INTERIOR :
	EXTERIOR;

    }

  //---------------------------------------------------------------------
 
  Classification_type  classify(const Facet& f) const
    {  
      return classify(f.first, f.second, get_alpha());
    }

  
  Classification_type  classify(const Cell_handle& s, 
				int i) const
    {  
      return classify(s, i, get_alpha());
    }

  Classification_type  classify(const Facet& f,
				const Coord_type& alpha) const
    {  
      return classify(f.first, f.second, alpha);
    }

  Classification_type  classify(const Cell_handle& s, 
				int i,
				const Coord_type& alpha) const;
  // Classifies the face `f' of the underlying Delaunay
  // tetrahedralization with respect to `A'.

  //---------------------------------------------------------------------

  Classification_type  classify(const Edge& e) const
    {  
      return classify(e.first, e.second, e.third, get_alpha());
    }

  
  Classification_type  classify(const Cell_handle& s, 
				int i,
 				int j) const
    {  
      return classify(s, i, j, get_alpha());
    }

  Classification_type  classify(const Edge& e,
				const Coord_type& alpha ) const
    {  
      return classify(e.first, e.second, e.third, alpha);
    }

  Classification_type  classify(const Cell_handle& s, 
				int i,
				int j,
				const Coord_type& alpha) const;
  // Classifies the edge `e' of the underlying Delaunay
  // tetrahedralization with respect to `A'.
  

  //---------------------------------------------------------------------

  Classification_type  classify(const Vertex_handle& v) const
    {
      return classify(v, get_alpha());
    }

  Classification_type  classify(const Vertex_handle& v,
				const Coord_type& alpha) const;
  // Classifies the vertex `v' of the underlying Delaunay
  // tetrahedralization with respect to `A'.


  //--------------------- NB COMPONENTS ---------------------------------
  int
  number_solid_components() const
    {
      return number_of_solid_components(get_alpha());
    }

  int
  number_of_solid_components() const
    {
      return number_of_solid_components(get_alpha());
    }

  int 
  number_solid_components(const Coord_type& alpha) const
    {
      return number_of_solid_components(alpha);
    }

  int
  number_of_solid_components(const Coord_type& alpha) const;
  // Determine the number of connected solid components 
  // takes time O(#alpha_shape) amortized if STL_HASH_TABLES
  //            O(#alpha_shape log n) otherwise

private:

  void traverse(Cell_handle pCell,
		Marked_cell_set& marked_cell_set, 
		const Coord_type alpha) const;
 
  //----------------------------------------------------------------------

public:

  Alpha_iterator find_optimal_alpha(int nb_components);
  // find the minimum alpha that satisfies the properties
  // (1) nb_components solid components
  // (2) all data points on the boundary or in its interior
  
private:

  Coord_type find_alpha_solid() const;
  // compute the minumum alpha such that all data points 
  // are either on the boundary or in the interior
  // not necessarily connected
  // starting point for searching 
  // takes O(#alpha_shape) time

  //--------------------- PREDICATES ------------------------------------

private:

  bool is_attached(const Cell_handle& s, const int& i) const
    {
      int i0=(i+1)&3, i1=(i+2)&3, i2=(i+3)&3;

      Bounded_side b = 
	Gt().side_of_bounded_sphere_3_object()(s->vertex(i0)->point(),
					       s->vertex(i1)->point(),
					       s->vertex(i2)->point(),
					       s->vertex(i)->point());
 
      return (b == ON_BOUNDED_SIDE) ? true : false;
    }
  
  bool is_attached(const Cell_handle& s, const int& i0, 
		   const int& i1, const int& i) const
    {
      Bounded_side b = 
	Gt().side_of_bounded_sphere_3_object()(s->vertex(i0)->point(),
					       s->vertex(i1)->point(),
					       s->vertex(i)->point());
 
      return (b == ON_BOUNDED_SIDE) ? true : false;
    }

  bool is_attached(const Vertex_handle& vh1, 
		   const Vertex_handle& vh2,
		   const Vertex_handle& vh3) const
    {
     Bounded_side b = 
	Gt().side_of_bounded_sphere_3_object()(vh1->point(),
					       vh2->point(),
					       vh3->point());
 
      return (b == ON_BOUNDED_SIDE) ? true : false; 
    }
 
  //------------------- GEOMETRIC PRIMITIVES ----------------------------

  Coord_type squared_radius(const Cell_handle& s) const
    {
      return Gt().compute_squared_radius_3_object()(s->vertex(0)->point(),
						    s->vertex(1)->point(),
						    s->vertex(2)->point(),
						    s->vertex(3)->point());
    }

  Coord_type squared_radius(const Cell_handle& s, const int& i) const
    {
      // test which one is faster TBC

      int i0=(i+1)&3, i1=(i+2)&3, i2=(i+3)&3;

      return Gt().compute_squared_radius_3_object()(s->vertex(i0)->point(),
						    s->vertex(i1)->point(),
						    s->vertex(i2)->point());
    }

  Coord_type squared_radius(const Cell_handle& s, 
			    const int& i, const int& j) const
    {
      return Gt().compute_squared_radius_3_object()(s->vertex(i)->point(),
						    s->vertex(j)->point());
    }

  //---------------------------------------------------------------------

private:
  // prevent default copy constructor and default assigment
  
  Alpha_shape_3(const Alpha_shape_3& A)
    {}

  Alpha_shape_3& operator=(const Alpha_shape_3& A)
    {}

  //---------------------------------------------------------------------
public:  
  void show_alpha_shape_faces(Geomview_stream &gv) const;


  // to Debug
  void show_interval_facet_map() const;

}; 



//---------------------------------------------------------------------
//--------------------- MEMBER FUNCTIONS-------------------------------
//---------------------------------------------------------------------


//--------------------- INITIALIZATION OF PRIVATE MEMBERS -------------
  
template <class Dt>
void 
Alpha_shape_3<Dt>::initialize_interval_cell_map()
{  
  Finite_cells_iterator cell_it, done = finite_cells_end();
  Coord_type alpha_f;

  for( cell_it = finite_cells_begin(); cell_it != done; ++cell_it) {
      alpha_f = squared_radius(cell_it);
      _interval_cell_map.insert(Interval_cell(alpha_f, cell_it));

      // cross references
      cell_it->set_alpha(alpha_f);
    }
}


//---------------------------------------------------------------------

template <class Dt>
void 
Alpha_shape_3<Dt>::initialize_interval_facet_map()
{
  Finite_facets_iterator face_it;  
  Facet f;
  Cell_handle pCell, pNeighbor ;
  Interval3 interval;

 
  for( face_it = finite_facets_begin(); 
       face_it != finite_facets_end(); 
       ++face_it)
    {
      pCell = face_it->first;
      int i = face_it->second;

      pNeighbor = pCell->neighbor(i);
      int Neigh_i = pNeighbor->index(pCell);
      f =  Facet(pCell, i);
       
      // not on the convex hull
      if(!is_infinite(pCell) && !is_infinite(pNeighbor))
	{ 
	  Coord_type squared_radius_Cell = find_interval(pCell);
	  Coord_type squared_radius_Neighbor = find_interval(pNeighbor);
	  //f =  Facet(pCell, i);

	  if (squared_radius_Neighbor < squared_radius_Cell) 
	    {	    
	      f =  Facet(pNeighbor, Neigh_i);
	      Coord_type coord_tmp = squared_radius_Cell;
	      squared_radius_Cell = squared_radius_Neighbor;
	      squared_radius_Neighbor = coord_tmp;
	    }

	  interval = (is_attached(pCell, i) || 
		      is_attached(pNeighbor, Neigh_i)) ?
	    make_triple(UNDEFINED,
			squared_radius_Cell, 
			squared_radius_Neighbor):
	    make_triple(squared_radius(pCell, i),
			squared_radius_Cell, 
			squared_radius_Neighbor);
	}
      else  // on the convex hull
	{
	  if(is_infinite(pCell)) {
	    f =  Facet(pNeighbor, Neigh_i);
	   CGAL_triangulation_assertion(!is_infinite(pNeighbor)); 
	   interval =  (is_attached(pNeighbor, 
				    pNeighbor->index(pCell))) ?
		    make_triple(UNDEFINED,
				pNeighbor->get_alpha(),
				Infinity):
		    make_triple(squared_radius(pNeighbor, 
					       pNeighbor->index(pCell)), 
				pNeighbor->get_alpha(), 
				Infinity);
	  }
	  else  { // is_infinite(pNeighbor)
	   CGAL_triangulation_assertion(!is_infinite(pCell)); 
	   interval = is_attached(pCell, i) ?
	     make_triple(UNDEFINED,
			 find_interval(pCell),				       
			 Infinity) :
	     make_triple(squared_radius(pCell, i),
			 find_interval(pCell),
			 Infinity);
	  }
	}
      _interval_facet_map.insert(Interval_facet(interval, f));

      // cross-links
      pCell->set_facet_ranges(i, interval);
      pNeighbor->set_facet_ranges( Neigh_i,interval);
    }

  // Remark:
  // The interval_facet_map will be sorted as follows
  // first the attached faces  sorted by order of alpha_mid
  // second  the unattached faces  by order of alpha_min
  
  // 
  // if we are in regularized mode we should sort differently
  // by the second third first Key 
  // struct LessIntervalRegular
  //  {
  //    // if we are in regularized mode we should sort differently
  //    // by the second third Key
  //    bool operator()(const Interval_facet& ie1, 
  //		    const Interval_facet& ie2)
  //    { return ie1.first.second < ie2.first.second ||
  //	(ie1.first.second == ie2.first.second &&
  //	 ie1.first.third < ie2.first.third); }
  //  };


}



//---------------------------------------------------------------------
// we use a different definition than Edelsbrunner and Muecke
// singular means not incident to any 3-dimensional face
// regular means incident to a 3-dimensional face
// the interval of a vertex is defined as follows :
//   SINGULAR   first value   REGULAR   second value    INTERIOR
//---------------------------------------------------------------------
template <class Dt>
void 
Alpha_shape_3<Dt>::initialize_interval_vertex_map()
{
  Coord_type alpha_mid_v;
  Coord_type alpha_max_v;
  Coord_type alpha_s;

  Finite_vertices_iterator vertex_it;
  // only finite vertexs
  for( vertex_it = finite_vertices_begin(); 
       vertex_it != finite_vertices_end(); 
       ++vertex_it)   {
    CGAL_triangulation_assertion (! is_infinite(vertex_it));
    alpha_max_v = 0;    
    alpha_mid_v = (!_interval_cell_map.empty() ?
		   (--_interval_cell_map.end())->first :
		   0);

    std::list<Cell_handle> incidents;
    incident_cells(vertex_it, back_inserter(incidents));
    typename std::list<Cell_handle>::iterator chit = incidents.begin();
    for( ; chit != incidents.end(); ++chit){
      if (is_infinite(*chit))  alpha_max_v = Infinity;
      else {
	alpha_s = find_interval(*chit);
	alpha_mid_v  = min(alpha_mid_v, alpha_s);
	if (alpha_max_v != Infinity)
	   alpha_max_v = max(alpha_max_v, alpha_s);
      }
    }
  

    Interval2 interval = std::make_pair(alpha_mid_v, alpha_max_v);
    _interval_vertex_map.insert(Interval_vertex(interval, vertex_it));

    // cross references
    vertex_it->set_range(interval);
  }
}

//---------------------------------------------------------------------


template <class Dt>
void 
Alpha_shape_3<Dt>::initialize_alpha_spectrum()
// merges the alpha thresholds of the simplices faces (and edges)
{

  // skip the attached faces 
  // <=> _interval_facet_map.first.first == UNDEFINED
  typename Interval_facet_map::iterator
    face_it = std::upper_bound(_interval_facet_map.begin(),
			       _interval_facet_map.end(),
			       UNDEFINED,
			       Less());

  // merge the alpha values of cells
  // with the alphamin of  unattached faces.
  // eliminate duplicate values due to for example attached faces
  // merge and copy from STL since assignment should be function object
	
  typename Interval_cell_map::iterator
    cell_it = _interval_cell_map.begin();

  _alpha_spectrum.reserve(_interval_cell_map.size() +
			  _interval_facet_map.size()/ 2 );
  // should be only the number of unattached faces
  // size_type nb_unattached_facets;
  // distance(face_it, _interval_facet_map.end(), nb_unattached_facets);
  // however the distance function is expensive

  while (face_it != _interval_facet_map.end() ||
	 cell_it != _interval_cell_map.end())
    {
      if (cell_it != _interval_cell_map.end() && 
	  (face_it == _interval_facet_map.end() ||
	   (*cell_it).first < (*face_it).first.first))
	{
	  assert(cell_it != _interval_cell_map.end());

	  if (_alpha_spectrum.empty() || 
	      _alpha_spectrum.back() < (*cell_it).first)
	    _alpha_spectrum.push_back((*cell_it).first);
	  cell_it++;
	}
      else
	{
	  assert (cell_it == _interval_cell_map.end() || 
		  (face_it != _interval_facet_map.end() &&
		   (*cell_it).first >= (*face_it).first.first));

	  if (_alpha_spectrum.empty() || 
	      _alpha_spectrum.back() < (*face_it).first.first)
	    _alpha_spectrum.push_back((*face_it).first.first);
	  face_it++;
	}
    }

}

//---------------------------------------------------------------------


template <class Dt>
std::istream& operator>>(std::istream& is,  const Alpha_shape_3<Dt>& A)
  // Reads a alpha shape from stream `is' and assigns it to
  // Unknown creationvariable. Precondition: The extract operator must
  // be defined for `Point'.
{}

//---------------------------------------------------------------------

template <class Dt>
std::ostream& operator<<(std::ostream& os,  const Alpha_shape_3<Dt>& A)
  // Inserts the alpha shape into the stream `os' as an indexed face set. 
  // Precondition: The insert operator must be defined for `Point'
{ 
  typedef typename Alpha_shape_3<Dt>::Vertex_handle Vertex_handle;
  typedef typename Alpha_shape_3<Dt>::Interval_vertex_map Interval_vertex_map;
  typename Interval_vertex_map::const_iterator vertex_alpha_it;

  const typename Alpha_shape_3<Dt>::Interval2* pInterval2;

  Unique_hash_map< Vertex_handle, int > V;

  int number_of_vertices = 0;

  typedef typename Alpha_shape_3<Dt>::Interval_facet_map Interval_facet_map;
  typename Interval_facet_map::const_iterator face_alpha_it;

  const typename Alpha_shape_3<Dt>::Interval3* pInterval;

  int i0, i1, i2;

  if (A.get_mode() == Alpha_shape_3<Dt>::REGULARIZED)
    {

      Vertex_handle v;
      for (vertex_alpha_it = A._interval_vertex_map.begin(); 
	   vertex_alpha_it != A._interval_vertex_map.end() &&
	     (*vertex_alpha_it).first.first < A.get_alpha();
	   ++vertex_alpha_it)
	{
	  pInterval2 = &(*vertex_alpha_it).first;

#ifdef DEBUG
	  Alpha_shape_3<Dt>::Coord_type alpha =
	    A.get_alpha();
	  Alpha_shape_3<Dt>::Coord_type alpha_min = 
	    pInterval2->first;
	  Alpha_shape_3<Dt>::Coord_type alpha_max = 
	    pInterval2->second;
#endif // DEBUG

	  if((pInterval2->second >= A.get_alpha()
	      || pInterval2->second == A.Infinity))
	    // alpha must be larger than the min boundary
	    // and alpha is smaller than the upper boundary
	    // which might be infinity 
	    // write the vertex
	    {
	      v = (*vertex_alpha_it).second;
	      assert(A.classify(v) ==
		     Alpha_shape_3<Dt>::REGULAR);

	      V[v] = number_of_vertices++;
	      os << v->point() << std::endl;
	    }
	}
  
      // the vertices are oriented counterclockwise

      typename Alpha_shape_3<Dt>::Cell_handle s;
      int i;

      for (face_alpha_it = A._interval_facet_map.begin(); 
	   face_alpha_it != A._interval_facet_map.end() &&
	     (*face_alpha_it).first.first < A.get_alpha();
	   ++face_alpha_it)
	{
	  pInterval = &(*face_alpha_it).first;

#ifdef DEBUG
	  Alpha_shape_3<Dt>::Coord_type alpha =
	    A.get_alpha();
	  Alpha_shape_3<Dt>::Coord_type alpha_mid = 
	    pInterval->second;
	  Alpha_shape_3<Dt>::Coord_type alpha_max = 
	    pInterval->third;
#endif // DEBUG

	  assert(pInterval->second != A.Infinity);
	  // since this happens only for convex hull of dimension 2
	  // thus singular

	  if(pInterval->second < A.get_alpha() &&
	     (pInterval->third >= A.get_alpha()
	      || pInterval->third == A.Infinity))
	    // alpha must be larger than the mid boundary
	    // and alpha is smaller than the upper boundary
	    // which might be infinity 
	    // visualize the boundary
	    {

	      s = (*face_alpha_it).second.first;
	      i = (*face_alpha_it).second.second;

	      // assure that all vertices are in ccw order
	      if (A.classify(s) == Alpha_shape_3<Dt>::EXTERIOR)
		{ 
		  // take the reverse cell
		  typename Alpha_shape_3<Dt>::Cell_handle 
		    pNeighbor = s->neighbor(i);
		  i = pNeighbor->index(s);
		  s = pNeighbor;
		}
	  
	      assert(A.classify(s) == Alpha_shape_3<Dt>::INTERIOR);

	      assert(A.classify(s, i) ==
		     Alpha_shape_3<Dt>::REGULAR);

	      int i0=(i+1)&3, i1=(i+2)&3, i2=(i+3)&3;

	      os << V[s->vertex(i0)] << ' ' 
		 << V[s->vertex(i1)] << ' ' 
		 << V[s->vertex(i2)] << std::endl;
	    }
	}
    }
  else  // A.get_mode() == GENERAL -----------------------------------------
    {
 
       Vertex_handle v;
     
      // write the regular vertices

      for (vertex_alpha_it = A._interval_vertex_map.begin(); 
	   vertex_alpha_it != A._interval_vertex_map.end() &&
	     (*vertex_alpha_it).first.first < A.get_alpha();
	   ++vertex_alpha_it)
	{
	  pInterval2 = &(*vertex_alpha_it).first;

	  if((pInterval2->second >= A.get_alpha()
	      || pInterval2->second == A.Infinity))
	    // alpha must be larger than the min boundary
	    // and alpha is smaller than the upper boundary
	    // which might be infinity 
	    // write the vertex
	    {
	      v = (*vertex_alpha_it).second;
	      CGAL_triangulation_assertion(A.classify(v) ==
					   Alpha_shape_3<Dt>::REGULAR);
	      V[v] = number_of_vertices++;
	      os << v->point() << std::endl;
	    }
	}
 
      // write the singular vertices
      for (; 
	   vertex_alpha_it != A._interval_vertex_map.end();
	   ++vertex_alpha_it)
	{
	  v = (*vertex_alpha_it).second;
	  CGAL_triangulation_assertion(A.classify(v) ==
				       Alpha_shape_3<Dt>::SINGULAR);

	  V[v] = number_of_vertices++;
	  os << v->point() << std::endl;
	}
 
      // the vertices are oriented counterclockwise

      typename Alpha_shape_3<Dt>::Cell_handle s;
      int i;

      for (face_alpha_it = A._interval_facet_map.begin(); 
	   face_alpha_it != A._interval_facet_map.end() &&
	     (*face_alpha_it).first.first < A.get_alpha();
	   ++face_alpha_it)
	{
	  pInterval = &(*face_alpha_it).first;

#ifdef DEBUG
	  Alpha_shape_3<Dt>::Coord_type alpha =
	    A.get_alpha();
	  Alpha_shape_3<Dt>::Coord_type alpha_min = 
	    pInterval->first;
	  Alpha_shape_3<Dt>::Coord_type alpha_mid = 
	    pInterval->second;
	  Alpha_shape_3<Dt>::Coord_type alpha_max = 
	    pInterval->third;
#endif // DEBUG
	  
	  if(pInterval->third >= A.get_alpha()
	     || pInterval->third == A.Infinity)
	    // if alpha is smaller than the upper boundary
	    // which might be infinity 
	    // visualize the boundary
	    {
	      s = (*face_alpha_it).second.first;
	      i = (*face_alpha_it).second.second;


	      // write the regular faces
	      if (pInterval->second != A.Infinity &&
		  pInterval->second < A.get_alpha())
		{
		  CGAL_triangulation_assertion(A.classify(s, i) ==
					       Alpha_shape_3<Dt>::REGULAR);
		  // assure that all vertices are in ccw order
		  if (A.classify(s) == Alpha_shape_3<Dt>::EXTERIOR)
		    { 
		      // take the reverse cell
		      typename Alpha_shape_3<Dt>::Cell_handle 
			pNeighbor = s->neighbor(i);
		      i = pNeighbor->index(s);
		      s = pNeighbor;
		    }
	  
		  CGAL_triangulation_assertion(A.classify(s) == 
					       Alpha_shape_3<Dt>::INTERIOR);

		  i0=(i+1)&3, i1=(i+2)&3, i2=(i+3)&3;

		  os << V[s->vertex(i0)] << ' ' 
		     << V[s->vertex(i1)] << ' ' 
		     << V[s->vertex(i2)] << std::endl;
		  
		}
	      else // (pInterval->second == A.Infinity || 
		   //  pInterval->second >= A.get_alpha()))

		// pInterval->second == A.Infinity happens only for convex hull
		// of dimension 2 thus singular
		{
		  // write the singular faces
		  if (pInterval->first != A.UNDEFINED)
		    {
		      CGAL_triangulation_assertion(A.classify(s, i) ==
						   Alpha_shape_3<Dt>::SINGULAR);
		      i0=(i+1)&3, i1=(i+2)&3, i2=(i+3)&3;

		      os << V[s->vertex(i0)] << ' ' 
			 << V[s->vertex(i1)] << ' ' 
			 << V[s->vertex(i2)] << std::endl;
		      
		    }	
		}
	    }
	}
    }
  
  return os;

}

//---------------------------------------------------------------------

template <class Dt>
void
Alpha_shape_3<Dt>::update_alpha_shape_vertex_list() const
{
  Alpha_shape_vertices_list.clear();
  typedef Alpha_shape_3<Dt>::Interval_vertex_map Interval_vertex_map;
  typename Interval_vertex_map::const_iterator vertex_alpha_it;

  const Alpha_shape_3<Dt>::Interval2* pInterval2;

  Vertex_handle v;
     
  // write the regular vertices

  for (vertex_alpha_it = _interval_vertex_map.begin(); 
       vertex_alpha_it != _interval_vertex_map.end() &&
	 (*vertex_alpha_it).first.first < get_alpha();
       ++vertex_alpha_it)
    {
      pInterval2 = &(*vertex_alpha_it).first;

      if((pInterval2->second >= get_alpha()
	  || pInterval2->second == Infinity))
	// alpha must be larger than the min boundary
	// and alpha is smaller than the upper boundary
	// which might be infinity 
	// write the vertex
	{
	  v = (*vertex_alpha_it).second;
	  CGAL_triangulation_assertion(classify(v) == 
				       Alpha_shape_3<Dt>::REGULAR);
	  Alpha_shape_vertices_list.push_back(v);
	}
    }
 
  if (get_mode() == Alpha_shape_3<Dt>::GENERAL)
    {
      // write the singular vertices
      for (; 
	   vertex_alpha_it != _interval_vertex_map.end();
	   ++vertex_alpha_it)
	{
	  v = (*vertex_alpha_it).second;
	  CGAL_triangulation_assertion(classify(v) == 
				       Alpha_shape_3<Dt>::SINGULAR);

	  Alpha_shape_vertices_list.push_back(v);
	}
    }
  use_vertex_cache = true;
}

//---------------------------------------------------------------------

template <class Dt>
void
Alpha_shape_3<Dt>::update_alpha_shape_facet_list() const
{
  Alpha_shape_facets_list.clear();
  // Writes the faces of the alpha shape `A' for the current 'alpha'-value
  // to the container where 'out' refers to. Returns an output iterator 
  // which is the end of the constructed range.
  typedef  Alpha_shape_3<Dt>::Interval_facet_map Interval_facet_map;
  typename Interval_facet_map::const_iterator face_alpha_it;

  const Alpha_shape_3<Dt>::Interval3* pInterval;

  if (get_mode() == Alpha_shape_3<Dt>::REGULARIZED)
    {
      // it is much faster looking at the sorted intervals 
      // than looking at all sorted cells
      // alpha must be larger than the mid boundary
      // and alpha is smaller than the upper boundary
      for (face_alpha_it = _interval_facet_map.begin(); 
	   face_alpha_it != _interval_facet_map.end() &&
	     (*face_alpha_it).first.first < get_alpha();
	   ++face_alpha_it)
	{
	  pInterval = &(*face_alpha_it).first;

	  CGAL_triangulation_assertion(pInterval->second != Infinity);
	  // since this happens only for convex hull of dimension 2
	  // thus singular

	  if(pInterval->second < get_alpha() &&
	     (pInterval->third >= get_alpha()
	      || pInterval->third == Infinity))
	    // alpha must be larger than the mid boundary
	    // and alpha is smaller than the upper boundary
	    // which might be infinity 
	    // visualize the boundary
	    {
	      CGAL_triangulation_assertion(classify(
                                 (*face_alpha_it).second.first,
				 (*face_alpha_it).second.second) ==
			Alpha_shape_3<Dt>::REGULAR);
	      Alpha_shape_facets_list.push_back(Facet((*face_alpha_it).second.first,
						      (*face_alpha_it).second.second));
	    }
	}
    }
  else  // get_mode() == GENERAL -------------------------------------------
    {
      // draw the faces
      for (face_alpha_it = _interval_facet_map.begin(); 
	   face_alpha_it != _interval_facet_map.end() &&
	     (*face_alpha_it).first.first < get_alpha();
	   ++face_alpha_it)
	{
	  pInterval = &(*face_alpha_it).first;

	  if (pInterval->first == UNDEFINED)
	    {
	      CGAL_triangulation_assertion(pInterval->second != Infinity);
	      // since this happens only for convex hull of dimension 2
	      // thus singular

	      if(pInterval->second < get_alpha() &&
		 (pInterval->third >= get_alpha()
		  || pInterval->third == Infinity))
		// alpha must be larger than the mid boundary
		// and alpha is smaller than the upper boundary
		// which might be infinity 
		// visualize the boundary
		{
		  CGAL_triangulation_assertion(classify(
                                     (*face_alpha_it).second.first,
				     (*face_alpha_it).second.second ) ==
			    Alpha_shape_3<Dt>::REGULAR );
		  Alpha_shape_facets_list.push_back(
		    Facet((*face_alpha_it).second.first,
			  (*face_alpha_it).second.second));
		}
	    }
	  else
	    {

	      if(pInterval->third >= get_alpha()
		 || pInterval->third == Infinity)
		// if alpha is smaller than the upper boundary
		// which might be infinity 
		// visualize the boundary
		{
		  CGAL_triangulation_assertion(classify(
                                     (*face_alpha_it).second.first,
				     (*face_alpha_it).second.second) ==
			    Alpha_shape_3<Dt>::REGULAR || 
			    classify((*face_alpha_it).second.first,
				     (*face_alpha_it).second.second) ==
			    Alpha_shape_3<Dt>::SINGULAR);
		  Alpha_shape_facets_list.push_back(Facet((*face_alpha_it).second.first,
							    (*face_alpha_it).second.second));
		}
	    }

	}

    }
  use_facet_cache = true;
}

//---------------------------------------------------------------------

template < class Dt >
typename Alpha_shape_3<Dt>::Classification_type  
Alpha_shape_3<Dt>::classify(const Cell_handle& s, 
			    int i,
			    const Coord_type& alpha) const
  // Classifies the face `f' of the underlying Delaunay
  // tetrahedralization with respect to `A'.
{ 
  // the version that uses a simplified version without crosslinks
  // is much slower
  if (is_infinite(s, i))
    return EXTERIOR;
    
  Interval3 interval = find_interval(Facet (s, i));
 
  if (alpha <= interval.second)
    {
      if (get_mode() == REGULARIZED ||
	  interval.first == UNDEFINED ||
	  alpha <= interval.first)
	return EXTERIOR;
      else // alpha > interval.first
	return SINGULAR;
    }
  else    // alpha > interval.second
    {
      if (interval.third == Infinity ||
	  alpha <= interval.third)
	return REGULAR;
      else // alpha > interval.third
	return INTERIOR;
    }
   
}
 
//---------------------------------------------------------------------

template < class Dt >
typename Alpha_shape_3<Dt>::Interval3
Alpha_shape_3<Dt>::
find_interval(const Cell_handle& s, int i, int j) const
{
 // edge is said to be singular if not incident to a tetrahedra
 // in the alpha-shapes
 // there is n no use to test if edge is attached or not if the mode
 // is set to REGULAR, in this case 
 // b_attached is set to true 
// interval.first is set to UNDEFINED for any edge.

 Coord_type alpha_min_e = UNDEFINED;
 Coord_type alpha_mid_e ;
 Coord_type alpha_max_e = 0;
 bool b_attached = false;

 
  //initialisation
  Cell_circulator ccirc = incident_cells(s,i,j),
    cdone(ccirc);
  do {
    if(is_infinite(ccirc) ) alpha_max_e =  Infinity;
    else alpha_mid_e = find_interval(ccirc);
  } while(++ccirc != cdone);

  //another turn
  do {
    if (!is_infinite(ccirc) ) {
     alpha_mid_e = min ( alpha_mid_e, find_interval(ccirc));
     if (alpha_max_e != Infinity) 
       alpha_max_e = max(alpha_max_e, find_interval(ccirc));
    }
  }while(++ccirc != cdone); 
 
   //-----set alpha_min_e in the GENERAL mode--------------------
  if( get_mode() == GENERAL) {
    Facet_circulator fcirc = incident_facets(s,i,j),
      fdone(fcirc);
    do {
      // test whether the vertex of ss opposite to *fcirc
      // is inside the sphere defined by the edge e = (s, i,j)
      Cell_handle ss = (*fcirc).first;
      int ii = (*fcirc).second;
      if (!is_infinite(ss->vertex(ii)) )
	b_attached = is_attached(s->vertex(i), s->vertex(j), ss->vertex(ii));
    } while(++fcirc != fdone);

    if(! b_attached) alpha_min_e = squared_radius(s,i,j);
  }

  return make_triple(alpha_min_e,alpha_mid_e,alpha_max_e);
}

template < class Dt >
typename Alpha_shape_3<Dt>::Classification_type  
Alpha_shape_3<Dt>::classify(const Cell_handle& s, 
			    int i,
			    int j,
			    const Coord_type& alpha) const
  // Classifies the edge `e' of the underlying Delaunay
  // tetrahedralization with respect to `A'.
{ 
//   // the version that uses a simplified version without crosslinks
//   // is much slower
//   if (is_infinite(s, i, j))     return EXTERIOR;
    
//   // FIX SUGGESTED BY JUR VAN DER BERG
//   //Interval3 interval = find_interval((const Edge)(s,i,j));
//   Interval3 interval = s->get_edge_ranges(i, j);
//   //  (*(_facet_interval_map.find(const_facet(s, i)))).second;
 
//   if (alpha <= interval.second)
//     {
//       if (get_mode() == REGULARIZED ||
// 	  interval.first == UNDEFINED ||
// 	  alpha <= interval.first)
// 	return EXTERIOR;
//       else // alpha > interval.first
// 	return SINGULAR;
//     }
//   else    // alpha > interval.second
//     {
//       if (interval.third == Infinity ||
// 	  alpha <= interval.third)
// 	return REGULAR;
//       else // alpha > interval.third
// 	return INTERIOR;
//     }
  
  if (is_infinite(s, i, j))     return EXTERIOR;
  Interval3 interval = find_interval(s,i,j);

//   std::cout << interval.first << "\t"
//             << interval.second << "\t"
// 	    << interval.third << "\t"; 
 
  // edge is said to be singular if not incident to a tetrahedra
  // in the alpha-shapes
  //terahedra with circumragius=alpha are considered outside
  if (interval.third != Infinity && alpha > interval.third) return INTERIOR;
  else if ( alpha > interval.second) return REGULAR;
  else if ( get_mode() == GENERAL && 
	    interval.first != UNDEFINED &&
	    alpha >= interval.first) return SINGULAR;
  else return EXTERIOR;
}

//---------------------------------------------------------------------

template < class Dt >
typename Alpha_shape_3<Dt>::Classification_type  
Alpha_shape_3<Dt>::classify(const Vertex_handle& v,
			    const Coord_type& alpha) const
  // Classifies the vertex `v' of the underlying Delaunay
  // tetrahedralization with respect to `A'.
{
  Interval2 interval = find_interval(v);
 
  if (alpha <= interval.first)
    {
      if (get_mode() == REGULARIZED) 
	return EXTERIOR;
      else // general => vertices are never exterior
	return SINGULAR;
    }
  else    // alpha > interval.first
    {
      if (interval.second == Infinity ||
	  alpha <= interval.second)
	return REGULAR;
      else // alpha > interval.second
	return INTERIOR;
    }
}

//--------------------- NB COMPONENTS ---------------------------------

template < class Dt >
int
Alpha_shape_3<Dt>::number_of_solid_components(const Coord_type& alpha) const
    // Determine the number of connected solid components 
    // takes time O(#alpha_shape) amortized if STL_HASH_TABLES
    //            O(#alpha_shape log n) otherwise
{
  typedef typename Marked_cell_set::Data Data;
  Marked_cell_set marked_cell_set(false);
  Finite_cells_iterator cell_it, done = finite_cells_end();
  int nb_solid_components = 0;

  // only finite simplices
  for( cell_it = finite_cells_begin(); cell_it != done; ++cell_it)
    {
      Cell_handle pCell = cell_it;
      assert(pCell != NULL);
      
      if (classify(pCell, alpha) == INTERIOR){
	Data& data = marked_cell_set[pCell];
	if(data == false) { 
	  // we traverse only interior simplices
	  data = true;
	  traverse(pCell, marked_cell_set, alpha);
	  nb_solid_components++;  
	}
      }
    }
  return nb_solid_components;
}


template < class Dt >
void Alpha_shape_3<Dt>::traverse(Cell_handle pCell,
				 Marked_cell_set& marked_cell_set,
				 const Coord_type alpha) const
{
  typedef typename Marked_cell_set::Data Data;
  std::list<Cell_handle> cells;
  cells.push_back(pCell);
  Cell_handle pNeighbor;

  while(! cells.empty()){
    pCell = cells.back();
    cells.pop_back();
    for (int i=0; i<=3; i++)
      {
	pNeighbor = pCell->neighbor(i);
	assert(pNeighbor != NULL);
	if (classify(pNeighbor, alpha) == INTERIOR){
	  Data& data = marked_cell_set[pNeighbor];
	  if(data == false){
	    data = true;
	    cells.push_back(pNeighbor);
	  }
	}
      }
  } 
}

//----------------------------------------------------------------------

template <class Dt>
typename Alpha_shape_3<Dt>::Alpha_iterator 
Alpha_shape_3<Dt>::find_optimal_alpha(int nb_components)
  // find the minimum alpha that satisfies the properties
  // (1) nb_components solid components
  // (2) all data points on the boundary or in its interior
{
  Coord_type alpha = find_alpha_solid();
  // from this alpha on the alpha_solid satisfies property (2)
  
  Alpha_iterator first = alpha_lower_bound(alpha);
  if (number_of_solid_components(alpha) == nb_components)
    {
      if ((first+1) < alpha_end())
	return (first+1);
      else
	return first;
    }

  // do binary search on the alpha values
  // number_of_solid_components() is a monotone function 
  // if we start with find_alpha_solid
  
  Alpha_iterator last = alpha_end();
  Alpha_iterator middle;
  
  std::ptrdiff_t len = last - first - 1;
  std::ptrdiff_t half;

  while (len > 0)
    {
      half = len / 2;
      middle = first + half;

#ifdef DEBUG
      cerr << "first : " << *first << " last : " 
           << ((first+len != last) ? *(first+len) : *(last-1))
	   << " mid : " << *middle 
	   << " nb comps : " << number_of_solid_components(*middle) << std::endl;
#endif // DEBUG

      if (number_of_solid_components(*middle) > nb_components)
	{
	  first = middle + 1;
	  len = len - half -1; 
	} 
      else // number_of_solid_components(*middle) <= nb_components
	{
	  len = half;
	}
    }
  if ((first+1) < alpha_end())
    return (first+1);
  else
    return first;
}  	

//----------------------------------------------------------------------

template <class Dt>
typename Alpha_shape_3<Dt>::Coord_type 
Alpha_shape_3<Dt>::find_alpha_solid() const
  // compute the minumum alpha such that all data points 
  // are either on the boundary or in the interior
  // not necessarily connected
  // starting point for searching 
  // takes O(#alpha_shape) time
{
  Coord_type alpha_solid = 0;
  Finite_vertices_iterator vit, done;
  for( vit = finite_vertices_begin(); 
       vit != finite_vertices_end(); ++vit) {
    alpha_solid = max(alpha_solid, find_interval(vit).first);
  }
  return alpha_solid;
}

template <class Dt>
void
Alpha_shape_3<Dt>::
show_interval_facet_map() const
{
  typename Interval_facet_map::const_iterator 
    ifmit = _interval_facet_map.begin(),
    ifmdone = _interval_facet_map.end();
  for( ; ifmit != ifmdone; ++ifmit) {
    Interval3 interval3 = ifmit->first;
    std::cout << std::endl;
    std::cout << interval3.first << "\t"
	      << interval3.second << "\t"
	      << interval3.third << "\t" << std::endl;
    Facet facet = ifmit->second;
    interval3 = find_interval(facet);
    std::cout << interval3.first << "\t"
	      << interval3.second << "\t"
	      << interval3.third << "\t" << std::endl;    
  }
}

//-------------------------------------------------------------------
CGAL_END_NAMESPACE
//-------------------------------------------------------------------

#include <CGAL/IO/alpha_shape_geomview_ostream_3.h>

#endif //CGAL_ALPHA_SHAPE_3_H
