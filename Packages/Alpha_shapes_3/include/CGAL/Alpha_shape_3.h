// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-2.0-I-20 $
// release_date  : $CGAL_Date: 1999/06/02 $
//
// file          : include/CGAL/Alpha_shape_3.h
// package       : Alpha_shapes_3(1.0)
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Tran Kai Frank DA <Frank.Da@sophia.inria.fr>
//
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ======================================================================

#ifndef CGAL_ALPHA_SHAPE_3_H
#define CGAL_ALPHA_SHAPE_3_H

#include <assert.h>
#include <CGAL/basic.h>

#include <set>
#include <map>


#include <CGAL/triple.h>
#include <vector>

#include <algorithm>
#include <utility>
#include <iostream>

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

  typedef typename Gt::Rep Rp;
  typedef typename Rp::FT Coord_type;

  typedef typename Gt::Point_3 Point;
  typedef typename Gt::Segment_3 Segment;
  //   typedef typename Gt::Triangle Triangle;
  //   typedef typename Gt::Tetrahedron Tetrahedron;

  typedef typename Dt::Cell_handle Cell_handle;
  typedef typename Dt::Vertex_handle Vertex_handle;
  typedef typename Dt::Cell Cell;
  typedef typename Dt::Vertex Vertex;
  typedef typename Dt::Facet Facet;
  typedef typename Dt::Edge Edge;

  typedef typename Dt::Cell_circulator Cell_circulator;
  typedef typename Dt::Facet_circulator Facet_circulator;

  typedef typename Dt::Cell_iterator Cell_iterator;
  typedef typename Dt::Facet_iterator Facet_iterator;
  typedef typename Dt::Edge_iterator Edge_iterator;
  typedef typename Dt::Vertex_iterator Vertex_iterator;

  typedef typename Dt::Locate_type Locate_type;

private:

  typedef long Key;
 
  typedef std::pair< Coord_type, Cell_handle > Interval_cell;
  typedef std::multimap< Coord_type, Cell_handle, std::less<Coord_type> > 
  Interval_cell_map;

  typedef triple<Coord_type, Coord_type , Coord_type> Interval3;

  typedef std::pair< Interval3, Facet > Interval_facet;
  typedef std::multimap< Interval3, Facet, std::less<Interval3> > 
  Interval_facet_map;

  //   typedef Cell_handle const const_void;
  //   typedef std::pair< const_void, int > const_Facet;
  //   typedef std::pair< const_void, int > const_Vertex;

  typedef std::pair< Interval3, Edge > Interval_edge;
  typedef std::multimap< Interval3, Edge, std::less<Interval3> > 
  Interval_edge_map;

  typedef std::pair< Coord_type, Coord_type > Interval2;
  typedef std::pair< Interval2, Vertex_handle > Interval_vertex;
  typedef std::multimap< Interval2, Vertex_handle, std::less<Interval2> > 
  Interval_vertex_map;

  typedef std::vector< Coord_type > Alpha_spectrum;
  
  typedef std::set< Key, std::less<Key> > Marked_cell_set;

public:

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
  typedef typename std::list< Facet >::iterator Alpha_shape_facets_iterator;


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

  std::list< Vertex_handle > Alpha_shape_vertices_list;
  std::list< Facet > Alpha_shape_facets_list;


  //------------------------- CONSTRUCTORS ------------------------------
 
  // Introduces an empty alpha-shape `A' for a positive
  // alpha-value `alpha'. Precondition: `alpha' >= 0.
  Alpha_shape_3(Coord_type alpha = 0, 
		Mode m = REGULARIZED)
    : _alpha(alpha), _mode(m), Infinity(-1), UNDEFINED(-2)
    {}
 
  // Introduces an alpha-shape `A' for a positive alpha-value
  // `alpha' that is initialized with the points in the range
  // from first to last

  template < class InputIterator >  
  Alpha_shape_3(const InputIterator& first,  
		const InputIterator& last,  
		const Coord_type& alpha = 0,
		Mode = REGULARIZED)
    {
      Dt::insert(first, last);
      if (dimension() == 3)
	{
	  // Compute the associated _interval_cell_map
	  initialize_interval_cell_map();
      
	  // Compute the associated _interval_facet_map
	  initialize_interval_facet_map();
	  
	  // Compute the associated _interval_edge_map
	  initialize_interval_edge_map();
	  
	  // Compute the associated _interval_vertex_map
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
	  // Compute the associated _interval_cell_map
	  initialize_interval_cell_map();

	  // Compute the associated _interval_facet_map
	  initialize_interval_facet_map();

	  // Compute the associated _interval_edge_map
	  initialize_interval_edge_map();
  
	  // Compute the associated _interval_vertex_map
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
      return (e.first)->get_edge_ranges(e.second, e.third);
    }
  
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
      return previous_alpha;
    }

  const Coord_type&  get_alpha() const
    // Returns the current alpha-value.
    {
      return _alpha;
    }
  

  const Coord_type&  get_nth_alpha(const int& n) const
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

  Mode set_mode(Mode mode = GENERAL )
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

  std::back_insert_iterator< std::list<Vertex_handle > >
  get_alpha_shape_vertices(std::back_insert_iterator<
			   std::list<Vertex_handle > > result) const;

  //---------------------------------------------------------------------

  std::back_insert_iterator<std::list<std::pair< Cell_handle, int > > >
  get_alpha_shape_facets(std::back_insert_iterator<
			 std::list<
			 std::pair< Cell_handle, int > > > result) const;

  //---------------------------------------------------------------------
public:

  Alpha_shape_vertices_iterator alpha_shape_vertices_begin()
    {
      Alpha_shape_vertices_list.erase(Alpha_shape_vertices_list.begin(), 
				      Alpha_shape_vertices_list.end());
      std::back_insert_iterator< std::list< Vertex_handle > > 
	V_it(Alpha_shape_vertices_list);
      get_alpha_shape_vertices(V_it);
      return Alpha_shape_vertices_list.begin();
    }

  Alpha_shape_vertices_iterator Alpha_shape_vertices_begin()
    {
      return alpha_shape_vertices_begin();
    }
  //---------------------------------------------------------------------

  Alpha_shape_vertices_iterator alpha_shape_vertices_end()
    {
      return Alpha_shape_vertices_list.end();
    }

  Alpha_shape_vertices_iterator Alpha_shape_vertices_end()
    {
      return alpha_shape_vertices_end();
    }

  //---------------------------------------------------------------------

  Alpha_shape_facets_iterator alpha_shape_facets_begin()
    {
      Alpha_shape_facets_list.erase(Alpha_shape_facets_list.begin(), 
				    Alpha_shape_facets_list.end());
      std::back_insert_iterator< std::list< Facet > > 
	E_it(Alpha_shape_facets_list);
      get_alpha_shape_facets(E_it);
      return Alpha_shape_facets_list.begin();
    }

  Alpha_shape_facets_iterator Alpha_shape_facets_begin()
    {
      return alpha_shape_facets_begin();
    }

  //---------------------------------------------------------------------

  Alpha_shape_facets_iterator alpha_shape_facets_end()
    {
      return Alpha_shape_facets_list.end();
    }

  Alpha_shape_facets_iterator Alpha_shape_facets_end()
    {
      return alpha_shape_facets_end();
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
	case FACE              : return classify(pCell, i, alpha);
	case CELL              : return classify(pCell, alpha);
	case OUTSIDE           : return EXTERIOR;
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
				const int& i) const
    {  
      return classify(s, i, get_alpha());
    }

  Classification_type  classify(const Facet& f,
				const Coord_type& alpha) const
    {  
      return classify(f.first, f.second, alpha);
    }

  Classification_type  classify(const Cell_handle& s, 
				const int& i,
				const Coord_type& alpha) const;
  // Classifies the face `f' of the underlying Delaunay
  // tetrahedralization with respect to `A'.

  //---------------------------------------------------------------------

  Classification_type  classify(const Edge& e) const
    {  
      return classify(e.first, e.second, e.third, get_alpha());
    }

  
  Classification_type  classify(const Cell_handle& s, 
				const int& i,
				const int& j) const
    {  
      return classify(s, i, j, get_alpha());
    }

  Classification_type  classify(const Edge& e,
				const Coord_type& alpha ) const
    {  
      return classify(e.first, e.second, e.third, alpha);
    }

  Classification_type  classify(const Cell_handle& s, 
				const int& i,
				const int& j,
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
      return number_of_solid_components(get_alpha());
    }

  int
  number_of_solid_components(const Coord_type& alpha) const;
  // Determine the number of connected solid components 
  // takes time O(#alpha_shape) amortized if STL_HASH_TABLES
  //            O(#alpha_shape log n) otherwise

private:

  void traverse(const Cell_handle& pCell,
		Marked_cell_set& marked_cell_set, 
		const Coord_type alpha) const;
 
  //----------------------------------------------------------------------

public:

  Alpha_iterator find_optimal_alpha(const int& nb_components);
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
  void show_alpha_shape_faces(Geomview_stream &gv);

}; 



//---------------------------------------------------------------------
//--------------------- MEMBER FUNCTIONS-------------------------------
//---------------------------------------------------------------------


//--------------------- INITIALIZATION OF PRIVATE MEMBERS -------------
  
template <class Dt>
void 
Alpha_shape_3<Dt>::initialize_interval_cell_map()
{  
  Cell_iterator cell_it;
  Cell_handle pCell;
  Coord_type alpha_f;

  for( cell_it = finite_cells_begin(); 
       cell_it != cells_end(); 
       ++cell_it)
    {
      pCell = cell_it->handle();
      
      alpha_f = squared_radius(pCell);

      _interval_cell_map.insert(Interval_cell(alpha_f, (pCell)));

      // cross references
      pCell->set_alpha(alpha_f);
    }
}


//---------------------------------------------------------------------

template <class Dt>
void 
Alpha_shape_3<Dt>::initialize_interval_facet_map()
{
  Facet_iterator face_it;  // TBC
  Facet f;
  Cell_handle pCell, pNeighbor ;
 
  for( face_it = finite_facets_begin(); 
       face_it != facets_end(); 
       ++face_it)
    {
      f = *face_it;

      pCell = f.first;
      int i = f.second;

      pNeighbor = pCell->neighbor(i);
      int Neigh_i = pNeighbor->index(pCell);

      Interval3 interval;
 
      // not on the convex hull
      if(!is_infinite(pCell) && !is_infinite(pNeighbor))
	{ 
	  Coord_type squared_radius_Cell = find_interval(pCell);
	  Coord_type squared_radius_Neighbor = find_interval(pNeighbor);
	  
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
	  if(is_infinite(pCell))
	    {
	      if (!is_infinite(pNeighbor))
		{
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
	      else
		{
		  // both simplices are infinite by definition unattached
		  // the face is finite by construction
		  assert(is_infinite(pNeighbor) && is_infinite(pCell));
		  interval = make_triple(
					 squared_radius(pCell, i), 
					 Infinity,
					 Infinity);
		}
	    }
	  else // is_infinite(pNeighbor)
	    {
	      assert(is_infinite(pNeighbor) && !is_infinite(pCell));
	      if (is_attached(pCell, i))
		interval = make_triple(UNDEFINED,
				       find_interval(pCell),
				       Infinity);
	      else
		interval = make_triple(squared_radius(pCell, i),
				       find_interval(pCell),
				       Infinity);

	    }
	}
      _interval_facet_map.insert(Interval_facet(interval, f));

      // cross-links
      (f.first)->set_facet_ranges(f.second, interval);
    }

  // Remark:
  // The interval_facet_map will be sorted as follows
  // first the attached faces on the convex hull
  // second                   not on the convex hull
  // third the un-attached faces on the convex hull
  // finally                     not on the convex hull
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

// template <class Dt>
// void 
// Alpha_shape_3<Dt>::initialize_interval_edge_map()
// {
//   // TBC since our definition of singular and regular differs from
//   // Edelsbrunner and Muecke's definition
//   // verify the relation between attached edges and singular and 
//   // regular edges

//   // _interval_edge_map.reserve(number_of_vertices()); TBC
    
//   Coord_type alpha_min_e;
//   Coord_type alpha_mid_e = (!_interval_cell_map.empty() ?
// 		     _interval_cell_map.end()->first :
// 		     0);;
//   Coord_type alpha_max_e = 0;
//   bool b_attached = false;
//   Cell_handle s;

//   Edge_iterator edge_it;
//   // only finite edges
//   for( edge_it = edges_begin(); 
//        edge_it != edges_end(); 
//        ++edge_it)
//     {
//       Edge e = (*edge_it);

//       if (! is_infinite(e))              // TBC
// 	{
	  
// 	  //-------------- examine incident faces --------------------------

// 	  Edge_circulator edge_circ = opposite_edges(e),
// 	    edge_done(edge_circ);

// 	  do 
// 	    { 
// 	      // the incident face (s, i) has vertex with index j, 
// 	      s = (*edge_circ).first;
// 	      int i = (*edge_circ).second;
// 	      int j = (*edge_circ).third;

// 	      if (is_infinite(std::make_pair(s,i)))   // TBC
// 		{
// 		  alpha_max_e = Infinity;
// 		}
// 	      else
// 		{
// 		  // test whether the vertex with index j
// 		  // is inside the sphere defined by the edge e = (s, cw(i,j), 
// 		  //                                               ccw(i,j))
// 		  b_attached = is_attached(s, cw(i,j), ccw(i,j), j);

// 		  Interval3 interval3 = find_interval(const_facet(s, i));
				
// 		  alpha_mid_e = (interval3.first != UNDEFINED) ?
// 		    CGAL::min(alpha_mid_e, interval3.first): 
// 		    CGAL::min(alpha_mid_e, interval3.second); 
			
// 		  if (alpha_max_e != Infinity)
// 		    {
// 		      alpha_max_e = (interval3.third != Infinity) ?
// 			CGAL::max(alpha_max_e, interval3.third):
// 			Infinity;
// 		    }
// 		}
// 	    }
// 	  while(++edge_circ != edge_done);

// 	  alpha_min_e = (b_attached ? UNDEFINED : 
// 			 squared_radius(e.first,
// 					e.second,
// 					e.third));

// 	  Interval3 interval = make_triple(alpha_min_e, 
// 						alpha_mid_e, 
// 						alpha_max_e);
// 	  _interval_edge_map.insert(Interval_edge(interval, e));

// 	  // cross references 
// 	  // we need a canonic description since otherwise we have problem
// 	  // with our access operation
// 	  // use the two vertices

// 	  s = e.first;
// 	  const Edge const_edge(s,e.second, e.third);
// 	  _edge_interval_map[const_edge] = interval;
// 	}
//     }
// }


//---------------------------------------------------------------------

template <class Dt>
void 
Alpha_shape_3<Dt>::initialize_interval_vertex_map()
{
  Coord_type alpha_mid_v;
  Coord_type alpha_max_v;
  Coord_type alpha_s;

  Vertex_iterator vertex_it;
  // only finite vertexs
  for( vertex_it = finite_vertices_begin(); 
       vertex_it != vertices_end(); 
       ++vertex_it)
    {
      Vertex_handle v = vertex_it->handle();

      if (!is_infinite(v))         // TBC
	{
	  Cell_handle s;

	  alpha_max_v = 0;    
	  alpha_mid_v = (!_interval_cell_map.empty() ?
			 (--_interval_cell_map.end())->first :
			 0);

	  //-------------- examine incident simplices --------------------
	  // we use a different definition than Edelsbrunner and Muecke
	  // singular means not incident to any 3-dimensional face
	  // regular means incident to a 3-dimensional face

	  //--------------------------------------------------------------
// 	  Cell_circulator cell_circ = v->incident_simplices(),
// 	    done(cell_circ);

// 	  if ((*cell_circ) != NULL)
// 	    {
// 	      do
// 		{
// 		  s = (*cell_circ);
// 		  if (is_infinite(s))
// 		    {
// 		      alpha_max_v = Infinity;
// 		      // continue;
// 		    }
// 		  else
// 		    {
// 		      alpha_s = find_interval(s);
// 		      // if we define singular as not incident to a 
// 		      // 3-dimensional cell
// 		      alpha_mid_v = CGAL::min(alpha_mid_v, alpha_s);
		    
// 		      if (alpha_max_v != Infinity)
// 			alpha_max_v = CGAL::max(alpha_max_v, alpha_s);
		    
// 		    }
// 		}
// 	      while(++cell_circ != done);
// 	    }
	  //---------------------------------------------------------------

	  // TBC if cell_circulator become available
	  // at the moment takes v*s time

	  Cell_iterator cell_it;
	    for( cell_it = all_cells_begin(); 
		 cell_it != cells_end(); 
		 ++cell_it)
	      {
		s = cell_it->handle();
		if (s->has_vertex(vertex_it->handle()))
		  {
		    if (is_infinite(s))
		      {
			alpha_max_v = Infinity;
			// continue;
		      }
		    else
		      {
			alpha_s = find_interval(s);
			// if we define singular as not incident to a
                        // 3-dimensional cell
			alpha_mid_v = CGAL::min(alpha_mid_v, alpha_s);
		    
			if (alpha_max_v != Infinity)
			  alpha_max_v = CGAL::max(alpha_max_v, alpha_s);
			
		      }
		  }
	      }

	  Interval2 interval = std::make_pair(alpha_mid_v, alpha_max_v);
	  _interval_vertex_map.insert(Interval_vertex(interval, 
						      vertex_it->handle()));

	  // cross references
	  vertex_it->handle()->set_range(interval);
	}
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

  // merge the maps which is sorted and contains the alpha-values
  // of the unattached faces and the triangles.
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
  typedef Alpha_shape_3<Dt>::Interval_vertex_map Interval_vertex_map;
  typename Interval_vertex_map::const_iterator vertex_alpha_it;

  const Alpha_shape_3<Dt>::Interval2* pInterval2;

  typedef long Key;

  std::map< Key, int, std::less< Key > > V;

  int number_of_vertices = 0;

  typedef Alpha_shape_3<Dt>::Interval_facet_map Interval_facet_map;
  typename Interval_facet_map::const_iterator face_alpha_it;

  const Alpha_shape_3<Dt>::Interval3* pInterval;

  int i0, i1, i2;

  if (A.get_mode() == Alpha_shape_3<Dt>::REGULARIZED)
    {

      typename Alpha_shape_3<Dt>::Vertex_handle v;
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

	      V[(Key)&(*v)] = number_of_vertices++;
	      os << v->point() << endl;
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

	      os << V[(Key)&(*s->vertex(i0))] << ' ' 
		 << V[(Key)&(*s->vertex(i1))] << ' ' 
		 << V[(Key)&(*s->vertex(i2))] << endl;
	    }
	}
    }
  else  // A.get_mode() == GENERAL -----------------------------------------
    {
 
       typename Alpha_shape_3<Dt>::Vertex_handle v;
     
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
	      V[(Key)&(*v)] = number_of_vertices++;
	      os << v->point() << endl;
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

	  V[(Key)&(*v)] = number_of_vertices++;
	  os << v->point() << endl;
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

		  os << V[(Key)&(*s->vertex(i0))] << ' ' 
		     << V[(Key)&(*s->vertex(i1))] << ' ' 
		     << V[(Key)&(*s->vertex(i2))] << endl;
		  
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

		      os << V[(Key)&(*s->vertex(i0))] << ' ' 
			 << V[(Key)&(*s->vertex(i1))] << ' ' 
			 << V[(Key)&(*s->vertex(i2))] << endl;
		      
		    }	
		}
	    }
	}
    }
  
  return os;

}

//---------------------------------------------------------------------

template <class Dt>
std::back_insert_iterator< 
       std::list<CGAL_TYPENAME_MSVC_NULL Alpha_shape_3<Dt>::Vertex_handle > >
Alpha_shape_3<Dt>::get_alpha_shape_vertices(std::back_insert_iterator<
					    std::list<Vertex_handle > > 
					    result) const
{
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
	  *result++ = v;
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

	  *result++ = v;
	}
    }
  return result;
}

//---------------------------------------------------------------------

template <class Dt>
std::back_insert_iterator< std::list<
    std::pair<CGAL_TYPENAME_MSVC_NULL Alpha_shape_3<Dt>::Cell_handle, int > > >
Alpha_shape_3<Dt>::get_alpha_shape_facets(std::back_insert_iterator<
					  std::list<
					  std::pair< Cell_handle, int > > > 
					  result) const
{
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
	      *result++ = Facet((*face_alpha_it).second.first,
				(*face_alpha_it).second.second);
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
				     (*face_alpha_it).second.second) ==
			    Alpha_shape_3<Dt>::REGULAR);
		  *result++ = Facet((*face_alpha_it).second.first,
				    (*face_alpha_it).second.second);
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
		  *result++ = Facet((*face_alpha_it).second.first,
				    (*face_alpha_it).second.second);
		}
	    }

	}

    }
  return result;
}

//---------------------------------------------------------------------

template < class Dt >
Alpha_shape_3<Dt>::Classification_type  
Alpha_shape_3<Dt>::classify(const Cell_handle& s, 
			    const int& i,
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
Alpha_shape_3<Dt>::Classification_type  
Alpha_shape_3<Dt>::classify(const Cell_handle& s, 
			    const int& i,
			    const int& j,
			    const Coord_type& alpha) const
  // Classifies the edge `e' of the underlying Delaunay
  // tetrahedralization with respect to `A'.
{ 
  // the version that uses a simplified version without crosslinks
  // is much slower
  if (is_infinite(s, i, j))
    return EXTERIOR;
    
  Interval3 interval = find_interval((const Edge)(s,i,j));
  //  (*(_facet_interval_map.find(const_facet(s, i)))).second;
 
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
Alpha_shape_3<Dt>::Classification_type  
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
  Marked_cell_set marked_cell_set;
  Cell_iterator cell_it;
  int nb_solid_components = 0;

  // only finite simplices
  for( cell_it = finite_cells_begin(); 
       cell_it != cells_end(); 
       ++cell_it)
    {
      Cell_handle pCell = cell_it->handle();
      assert(pCell != NULL);

      if (classify(pCell, alpha) == INTERIOR &&
	  ((marked_cell_set.insert((Key)&(*pCell)))).second)
	// we traverse only interior simplices
	{
	  traverse(pCell, marked_cell_set, alpha);
	  nb_solid_components++;  
	}
    }
  return nb_solid_components;
}

//----------------------------------------------------------------------

template < class Dt >
void Alpha_shape_3<Dt>::traverse(const Cell_handle& pCell,
				 Marked_cell_set& marked_cell_set, 
				 const Coord_type alpha) const
{
  for (int i=0; i<=3; i++)
    { 
      Cell_handle pNeighbor = pCell->neighbor(i);
      assert(pNeighbor != NULL); 
      if (classify(pNeighbor, alpha) == INTERIOR &&
	  ((marked_cell_set.insert((Key)&(*pNeighbor)))).second)
	traverse(pNeighbor, marked_cell_set, alpha);
    }
}

//----------------------------------------------------------------------

template <class Dt>
Alpha_shape_3<Dt>::Alpha_iterator 
Alpha_shape_3<Dt>::find_optimal_alpha(const int& nb_components)
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
  
  ptrdiff_t len = last - first - 1;
  ptrdiff_t half;

  while (len > 0)
    {
      half = len / 2;
      middle = first + half;

#ifdef DEBUG
      cerr << "first : " << *first << " last : " 
           << ((first+len != last) ? *(first+len) : *(last-1))
	   << " mid : " << *middle 
	   << " nb comps : " << number_of_solid_components(*middle) << endl;
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
Alpha_shape_3<Dt>::Coord_type 
Alpha_shape_3<Dt>::find_alpha_solid() const
  // compute the minumum alpha such that all data points 
  // are either on the boundary or in the interior
  // not necessarily connected
  // starting point for searching 
  // takes O(#alpha_shape) time
{
  Coord_type alpha_solid = 0;

  Vertex_iterator vertex_it;
  
  // at the moment all finite + infinite vertices
  for( vertex_it = all_vertices_begin(); 
       vertex_it != vertices_end();
       ++vertex_it)
    {
      if (!is_infinite(vertex_it->handle()))
	{
	  // consider only finite vertices
	  Coord_type alpha_min_v = (--_interval_cell_map.end())->first;

	  //------------------------------------------
// 	    Cell_circulator cell_circ =
// 	    (*vertex_it)->incident_simplices(),
// 	    done(cell_circ);
// 	    do
// 	    {
// 	    Cell_handle s = (*cell_circ);
// 	    if (! is_infinite(s))
// 	    alpha_min_v = CGAL::min(find_interval(s),
// 	    alpha_min_v);
// 	    }
// 	    while (++cell_circ != done);
	  //--------------------------------------------

	  // TBC if cell_circulator become available
	  // at the moment takes v*s time
	    
	  Cell_iterator cell_it;
	  for( cell_it = all_cells_begin(); 
	       cell_it != cells_end(); 
	       ++cell_it)
	    {
	      Cell_handle s = cell_it->handle();
	      if (s->has_vertex(vertex_it->handle()))
		{
		  if (! is_infinite(s))
		    {
		      alpha_min_v = CGAL::min(find_interval(s),
					      alpha_min_v);
		    }
		}
	    }
	    
	  alpha_solid = CGAL::max(alpha_min_v, alpha_solid);
	}
    }
  return alpha_solid;
}

//-------------------------------------------------------------------
CGAL_END_NAMESPACE
//-------------------------------------------------------------------

#include <CGAL/IO/alpha_shape_geomview_ostream_3.h>

#endif //CGAL_ALPHA_SHAPE_3_H
