// Copyright (c) 1997  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the so
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
//                 Mariette Yvinec <Mariette.Yvinec@sophia.inria.fr>

#ifndef CGAL_ALPHA_SHAPE_3_H
#define CGAL_ALPHA_SHAPE_3_H

#include <CGAL/basic.h>

#include <cassert>
#include <set>
#include <map>
#include <list>
#include <vector>
#include <algorithm>
#include <utility>
#include <iostream>

#include <CGAL/utility.h>
#include <CGAL/Unique_hash_map.h>
#include <CGAL/Alpha_shape_vertex_base_3.h>
#include <CGAL/Alpha_shape_cell_base_3.h>
#ifdef CGAL_USE_GEOMVIEW
#include <CGAL/IO/Geomview_stream.h>  // TBC
#endif

//-------------------------------------------------------------------
CGAL_BEGIN_NAMESPACE
//-------------------------------------------------------------------

template < class Dt >
class Alpha_shape_3 : public Dt
{
  // DEFINITION The class Alpha_shape_3<Dt> represents the family
  // of alpha-shapes of points for all alpha. It
  // maintains the underlying Delaunay tetrahedralization which represents
  // connectivity and order among its simplices. 
  // Each k-dimensional simplex of
  // the Delaunay tetrahedralization is associated with  three 
  // alpha values specifying the status of the $k$-simplex
  // alpha_min 
  // alpha_mid
  // alpha_max
  // alpha_min exist for k<3 if the $k$-simplex is a Gabriel simplex
  // alpha_max exist for k<3 if the simplex is not on the convex hull
  // Celsl have only one value : 
  // a cell is  EXTERIOR for  alpha     <=  alpha_mid
  //            INTERIOR for  alpha_mid <   alpha
  // Edges and  facets :
  // they are EXTERIOR for  alpha     <=  alpha_min
  //          SINGULAR for  alpha_min < alpha <= alpha_mid 
  //          REGULAR  for  alpha_mid < alpha <= alpha_max
  //          INTERIOR for  alpha_max < alpha
  // Vertices  are SINGULAR for alpha <= alpha_min
  //               REGULAR for  alpha_min < alpha <=  alpha_max
  //              SUPER_REGULAR for alpha_mid < alpha <=  alpha_max
  //               INTERNAL for alpha_max < alpha
  //------------------------- TYPES ------------------------------------

public:

  typedef typename Dt::Geom_traits Gt;
  typedef typename Dt::Triangulation_data_structure Tds;

  typedef typename Gt::FT Coord_type;
  typedef Coord_type      NT;

  typedef typename Gt::Point_3 Point;
  
  typedef typename Dt::Cell_handle Cell_handle;
  typedef typename Dt::Vertex_handle Vertex_handle;
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

  using Dt::dimension;
  using Dt::finite_facets_begin;
  using Dt::finite_facets_end;
  using Dt::finite_edges_begin;
  using Dt::finite_edges_end;
  using Dt::finite_vertices_begin;
  using Dt::finite_vertices_end;
  using Dt::finite_cells_begin;
  using Dt::finite_cells_end;
  using Dt::VERTEX;
  using Dt::EDGE;
  using Dt::FACET;
  using Dt::CELL;
  using Dt::OUTSIDE_CONVEX_HULL;
  using Dt::OUTSIDE_AFFINE_HULL;
  using Dt::vertex_triple_index;

  enum Classification_type {EXTERIOR, 
			    SINGULAR, 
			    REGULAR, 
			    SUPER_REGULAR,
			    INTERIOR};
  //  a k-dimensional face of the underlying Delaunay tetrahedralization 
  //  is
  // `EXTERIOR' if the face does not belong to the alpha-complex.
  // `SINGULAR' if the face belongs to the alpha-complex
  //            but is incident to no  higher-dimensional face of the complex
  // `REGULAR' if the face  belongs to the boundary of the alpha-complex
  //            and is incident to a higher-dimensional face 
  //           of the alpha-complex
  // SUPER_REGULAR (used for vertices only)
  //              if the vertices is on the boundary of the  alpha_complex
  //              and  belongs to at least one INTERIOR cell
  // `INTERIOR' if the cell belongs to the alpha-complex, 
  //            but does not belong to the boundary of the alpha-complex

  enum Mode {GENERAL, REGULARIZED};
  // In general, an alpha_complex is a non-connected, non-regular complex
  // Its regularized version is the subcomplex formed by the INTERIOR
  // cells and their faces

  typedef CGAL::Alpha_status< NT >          Alpha_status;
  typedef std::vector< NT >                 Alpha_spectrum;

  typedef std::multimap< NT, Cell_handle >  Alpha_cell_map;
  typedef std::multimap< NT, Facet>         Alpha_facet_map;
  typedef std::multimap< NT, Edge >         Alpha_edge_map;
  typedef std::multimap< NT, Vertex_handle> Alpha_vertex_map;

  typedef std::pair<Vertex_handle, Vertex_handle> Vertex_handle_pair;
  typedef std::map<Vertex_handle_pair,Alpha_status*> Edge_alpha_map;

  typedef typename std::list< Vertex_handle >::iterator 
                                            Alpha_shape_vertices_iterator;
  typedef typename std::list< Facet >::iterator
                                            Alpha_shape_facets_iterator;
  
  //test if a cell is exterior to the alphashape
  class Exterior_cell_test{
    const Alpha_shape_3 * _as;
  public:
    Exterior_cell_test() {}
    Exterior_cell_test(const Alpha_shape_3 * as) {_as = as;}
    bool operator() ( const Finite_cells_iterator& fci) const {
      return _as->classify(fci) == EXTERIOR ;
    }
  };

  typedef Filter_iterator< Finite_cells_iterator, Exterior_cell_test>
  Alpha_shape_cells_iterator;
  typedef typename Alpha_spectrum::const_iterator Alpha_iterator;
  // An iterator that allow to traverse the sorted sequence of
  // different alpha-values. The iterator is bidirectional and
  // non-mutable. Its value-type is NT

private:
  typedef Unique_hash_map<Cell_handle, bool > Marked_cell_set;

private:
  NT _alpha;
  Mode _mode;
  mutable bool use_vertex_cache;
  mutable bool use_facet_cache;

  // only finite facets and simplices are inserted into the maps
  Alpha_cell_map     alpha_cell_map;
  //Alpha_facet_map    alpha_max_facet_map;
  Alpha_facet_map    alpha_mid_facet_map;
  Alpha_facet_map    alpha_min_facet_map;
  //Alpha_edge_map     alpha_max_edge_map;
  Alpha_edge_map     alpha_mid_edge_map;
  Alpha_edge_map     alpha_min_edge_map;
  //Alpha_vertex_map   alpha_min_vertex_map;
  Alpha_vertex_map   alpha_mid_vertex_map;
  //Alpha_vertex_map   alpha_max_vertex_map;

  Alpha_spectrum alpha_spectrum;

  Edge_alpha_map edge_alpha_map;

  mutable std::list< Vertex_handle > alpha_shape_vertices_list;
  mutable std::list< Facet > alpha_shape_facets_list;



  ///////////////////TODO :
  // uitliser un compact container ou un vcetor
  // pour stocker leas alphastatus des facets and edges
  // plutot que de faire des news
  ///////////////////////////////////////


  //------------------------- CONSTRUCTORS ------------------------------
public:
  // Introduces an empty alpha-shape `A' for a positive
  // alpha-value `alpha'. Precondition: `alpha' >= 0.
  Alpha_shape_3(NT alpha = 0, 
		Mode m = REGULARIZED)
    : _alpha(alpha), _mode(m), 
      use_vertex_cache(false), use_facet_cache(false)
    {}

  Alpha_shape_3(Dt& dt, NT alpha = 0, Mode m = REGULARIZED)
    :_alpha(alpha), _mode(m), 
    use_vertex_cache(false), use_facet_cache(false)
    {
      Dt::swap(dt);
      if (dimension() == 3) initialize_alpha();
    }
 
  // Introduces an alpha-shape `A' for a positive alpha-value
  // `alpha' that is initialized with the points in the range
  // from first to last

  template < class InputIterator >  
  Alpha_shape_3(const InputIterator& first,  
		const InputIterator& last,  
		const NT& alpha = 0,
		Mode m = REGULARIZED)
    : _alpha(alpha), _mode(m), 
    use_vertex_cache(false), use_facet_cache(false)
    {
      Dt::insert(first, last);
      if (dimension() == 3)	  initialize_alpha();
    }
 
public:

  //----------------------- OPERATIONS ---------------------------------


  template < class InputIterator >  
  int make_alpha_shape(const InputIterator& first, 
		       const InputIterator& last)
    {
      clear();
      int n = Dt::insert(first, last);
      if (dimension() == 3)	  initialize_alpha();
      return n;
    }

  // Introduces an alpha-shape `A' for a positive alpha-value
  // `alpha' that is initialized with the points in the range
  // from first to last

private :

  //--------------------- INITIALIZATION OF PRIVATE MEMBERS -----------
  // called with reinitialize=false on first initialization
  // reinitialize=true when switching the mode from REGULAR to GENERAL
  void initialize_alpha_cell_map();
  void initialize_alpha_facet_maps(bool reinitialize = false);
  void initialize_alpha_edge_maps(bool reinitialize = false);
  void initialize_alpha_vertex_maps(bool reinitialize = false);
  void initialize_alpha_spectrum();
  void initialize_alpha(bool reinitialize = false) {
    if (!reinitialize) initialize_alpha_cell_map();
    initialize_alpha_facet_maps(reinitialize);
    initialize_alpha_edge_maps(reinitialize);
    initialize_alpha_vertex_maps(reinitialize);
    initialize_alpha_spectrum();
  }

private :
  Vertex_handle_pair
  make_vertex_handle_pair( Vertex_handle v1, Vertex_handle v2) const {
    return v1 < v2 ? std::make_pair(v1,v2)
                   : std::make_pair(v2,v1);
  }

  //---------------------------------------------------------------------

public:

  void clear()
    {
      // clears the structure
      Finite_facets_iterator fit;
      for( fit = finite_facets_begin(); 
	    fit != finite_facets_end(); 
	    ++fit){
	 delete fit->first->get_facet_status(fit->second);
       }

      Finite_edges_iterator eit;
      for( eit = finite_edges_begin(); eit != finite_edges_end(); ++eit){
	Vertex_handle_pair vhp = make_vertex_handle_pair(
				  eit->first->vertex(eit->second),
				  eit->first->vertex(eit->third));
	delete  edge_alpha_map.find(vhp)->second;
      }

      Finite_vertices_iterator vit = finite_vertices_begin();
      for(  ; vit != finite_vertices_end() ; ++vit) {
	delete vit->get_alpha_status();
      }

      Dt::clear();

      alpha_cell_map.clear();
      //alpha_max_facet_map.clear();
      alpha_mid_facet_map.clear();
      alpha_min_facet_map.clear();
      //alpha_max_edge_map.clear();
      alpha_mid_edge_map.clear();
      alpha_min_edge_map.clear();
      //alpha_max_vertex_map.clear();
      alpha_mid_vertex_map.clear();
      //alpha_min_vertex_map.clear();

      alpha_spectrum.clear();

      alpha_shape_vertices_list.clear();
      alpha_shape_facets_list.clear();

      set_alpha(0); 
      set_mode(REGULARIZED);
      use_vertex_cache = false;
      use_facet_cache = false;

    }

  //---------------------------------------------------------------------

public:

  NT set_alpha(const NT& alpha)
    // Sets the alpha-value to `alpha'. Precondition: `alpha' >= 0.
    // Returns the previous alpha
    {
      NT previous_alpha = _alpha;
      _alpha = alpha;
      use_vertex_cache = false;
      use_facet_cache = false;
      return previous_alpha;
    }

  const NT&  get_alpha() const
    // Returns the current alpha-value.
    {
      return _alpha;
    }
  

  const NT&  get_nth_alpha(int n) const
    // Returns the n-th alpha-value.
    // n < size()
    {
      assert( n > 0 && n <= static_cast<int>(alpha_spectrum.size()) );
      return alpha_spectrum[n-1];
    }
  
  int number_of_alphas() const
    // Returns the number of different alpha-values
    {
      return alpha_spectrum.size();
    }

  
  //---------------------------------------------------------------------

private:

  // the dynamic version is not yet implemented
  // desactivate the tetrahedralization member functions
  void insert(const Point& p) {}
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
      if (previous_mode != _mode) {
	initialize_alpha(true);
	use_vertex_cache = false;
	use_facet_cache = false;
      }
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
    if(!use_vertex_cache)    update_alpha_shape_vertex_list();
    return alpha_shape_vertices_list.begin();
  }

  Alpha_shape_vertices_iterator Alpha_shape_vertices_begin() const
    {
      return alpha_shape_vertices_begin();
    }
  //---------------------------------------------------------------------

  Alpha_shape_vertices_iterator alpha_shape_vertices_end() const
    {
      return alpha_shape_vertices_list.end();
    }

  Alpha_shape_vertices_iterator Alpha_shape_vertices_end() const
    {
      return alpha_shape_vertices_end();
    }

  //---------------------------------------------------------------------

  Alpha_shape_facets_iterator alpha_shape_facets_begin() const
    {
      if(! use_facet_cache)	update_alpha_shape_facet_list();
      return alpha_shape_facets_list.begin();
    }

  Alpha_shape_facets_iterator Alpha_shape_facets_begin() const
    {
      return alpha_shape_facets_begin();
    }

  //---------------------------------------------------------------------

  Alpha_shape_facets_iterator alpha_shape_facets_end() const
    {
      return alpha_shape_facets_list.end();
    }

  Alpha_shape_facets_iterator Alpha_shape_facets_end() const
    {
      return alpha_shape_facets_end();
    }

  Alpha_shape_cells_iterator alpha_shape_cells_begin() const 
    {
      return filter_iterator(finite_cells_end(),
			     Exterior_cell_test(this),
			     finite_cells_begin());
    }
  
  Alpha_shape_cells_iterator alpha_shape_cells_end() const
    {
      return filter_iterator(finite_cells_end(),
			     Exterior_cell_test(this));
    }


public: 
  
  // Traversal of the alpha-Values
  // 
  // The alpha shape class defines an iterator that allows to
  // visit the sorted sequence of alpha-values. This iterator is
  // non-mutable and bidirectional. Its value type is NT.

  Alpha_iterator alpha_begin() const { return alpha_spectrum.begin(); }
  Alpha_iterator alpha_end() const {return alpha_spectrum.end();}

  Alpha_iterator alpha_find(const NT& alpha) const
    // Returns an iterator pointing to an element with alpha-value
    // `alpha', or the corresponding past-the-end iterator if such an
    // element is not found.
    {
      return std::find(alpha_spectrum.begin(),
		       alpha_spectrum.end(),
		       alpha);
    }

  Alpha_iterator alpha_lower_bound(const NT& alpha) const
    // Returns an iterator pointing to the first element with
    // alpha-value not less than `alpha'.
    {
      return std::lower_bound(alpha_spectrum.begin(),
			      alpha_spectrum.end(),
			      alpha);
    }

  Alpha_iterator alpha_upper_bound(const NT& alpha) const
    // Returns an iterator pointing to the first element with
    // alpha-value greater than `alpha'.
    {
      return std::upper_bound(alpha_spectrum.begin(),
			      alpha_spectrum.end(),
			      alpha);
    }

  //--------------------- PREDICATES -----------------------------------
private:
  Classification_type classify(const Alpha_status* as,
			       const NT& alpha) const;

public:
  Classification_type  classify(const Point& p) const
    {
      return classify(p, get_alpha());
    }

  
  Classification_type  classify(const Point& p,   
				const NT& alpha) const
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
				const NT& alpha) const
    // Classifies the cell `f' of the underlying Delaunay
    // tetrahedralization with respect to `A'.
    // we consider open spheres :
    // s->radius == alpha => f exterior
    {
      if (is_infinite(s)) return EXTERIOR;
      return (s->get_alpha() < alpha) ? INTERIOR : EXTERIOR;
    }

  //---------------------------------------------------------------------
 
  Classification_type  classify(const Facet& f) const
    {  
      return classify(f.first, f.second, get_alpha());
    }

  
  Classification_type  classify(const Cell_handle& s, 	int i) const
    {  
      return classify(s, i, get_alpha());
    }

  Classification_type  classify(const Facet& f,	const NT& alpha) const
    {  
      return classify(f.first, f.second, alpha);
    }

  Classification_type  classify(const Cell_handle& s, 
				int i,
				const NT& alpha) const;
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
				const NT& alpha ) const
    {  
      return classify(e.first, e.second, e.third, alpha);
    }

  Classification_type  classify(const Cell_handle& s, 
				int i,
				int j,
				const NT& alpha) const;
  // Classifies the edge `e' of the underlying Delaunay
  // tetrahedralization with respect to `A'.
  

  //---------------------------------------------------------------------

  Classification_type  classify(const Vertex_handle& v) const
    {
      return classify(v, get_alpha());
    }

  Classification_type  classify(const Vertex_handle& v,
				const NT& alpha) const;
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
  number_solid_components(const NT& alpha) const
    {
      return number_of_solid_components(alpha);
    }

  int
  number_of_solid_components(const NT& alpha) const;
  // Determine the number of connected solid components 
  // takes time O(#alpha_shape) amortized if STL_HASH_TABLES
  //            O(#alpha_shape log n) otherwise

private:

  void traverse(Cell_handle pCell,
		Marked_cell_set& marked_cell_set, 
		const NT alpha) const;
 
  //----------------------------------------------------------------------

public:

  Alpha_iterator find_optimal_alpha(int nb_components);
  // find the minimum alpha that satisfies the properties
  // (1) nb_components solid components
  // (2) all data points on the boundary or in its interior
  
private:

  NT find_alpha_solid() const;
  // compute the minumum alpha such that all data points 
  // are either on the boundary or in the interior
  // not necessarily connected
  // starting point for searching 
  // takes O(#alpha_shape) time


  //------------------- GEOMETRIC PRIMITIVES ----------------------------

  NT squared_radius(const Cell_handle& s) const
    {
      return Gt().compute_squared_radius_3_object()(s->vertex(0)->point(),
						    s->vertex(1)->point(),
						    s->vertex(2)->point(),
						    s->vertex(3)->point());
    }

  NT squared_radius(const Cell_handle& s, const int& i) const
    {
      return Gt().compute_squared_radius_3_object() (
		  s->vertex(vertex_triple_index(i,0))->point(),
		  s->vertex(vertex_triple_index(i,1))->point(),
		  s->vertex(vertex_triple_index(i,2))->point());
    }

  NT squared_radius(const Facet& f) const {
    return squared_radius(f.first, f.second);
  }

  NT squared_radius(const Cell_handle& s, 
			    const int& i, const int& j) const
    {
      return Gt().compute_squared_radius_3_object()(s->vertex(i)->point(),
						    s->vertex(j)->point());
    }

  NT squared_radius(const Edge& e) const {
   return  squared_radius(e.first,e.second,e.third);
  }

  //---------------------------------------------------------------------

private:
  // prevent default copy constructor and default assigment
  Alpha_shape_3(const Alpha_shape_3& A)    {}
  void operator=(const Alpha_shape_3& A)    {}

  //---------------------------------------------------------------------
public:  
#ifdef CGAL_USE_GEOMVIEW
  void show_alpha_shape_faces(Geomview_stream &gv) const;
#endif


  // to Debug
  void print_maps(); //should be const
  void print_alphas() const;
  void print_alpha_status( const Alpha_status* as) const;
  void print_alpha_status_vertex( const Alpha_status* as) const;
}; 



//---------------------------------------------------------------------
//--------------------- MEMBER FUNCTIONS-------------------------------
//---------------------------------------------------------------------


//--------------------- INITIALIZATION OF PRIVATE MEMBERS -------------
  
template <class Dt>
void 
Alpha_shape_3<Dt>::initialize_alpha_cell_map()
{ 
  Finite_cells_iterator cell_it, done = finite_cells_end();
  NT alpha ;

  for( cell_it = finite_cells_begin(); cell_it != done; ++cell_it) {
    alpha = squared_radius(cell_it);
    alpha_cell_map.insert(typename Alpha_cell_map::value_type(alpha, cell_it));

    // cross references
    cell_it->set_alpha(alpha);
  }
  return;
}


//---------------------------------------------------------------------

template <class Dt>
void 
Alpha_shape_3<Dt>::initialize_alpha_facet_maps(bool reinitialize)
{
  Finite_facets_iterator fit;  
  Cell_handle pCell, pNeighbor ;
  int i, iNeigh;
  Alpha_status* as;

  if (!reinitialize) {
    NT alpha_max, alpha_mid;
    for( fit = finite_facets_begin(); 
	 fit != finite_facets_end(); ++fit)   {

      as = new Alpha_status();
  
      pCell = fit->first;
      i = fit->second;
      pNeighbor = pCell->neighbor(i);
      iNeigh = pNeighbor->index(pCell);
  
      // not on the convex hull
      if(!is_infinite(pCell) && !is_infinite(pNeighbor))	{ 
	NT alpha_Cell = pCell->get_alpha();
	NT alpha_Neighbor = pNeighbor->get_alpha();
	if ( alpha_Cell < alpha_Neighbor) {
	  alpha_mid = alpha_Cell;
	  alpha_max = alpha_Neighbor;
	}
	else {
	  alpha_mid = alpha_Neighbor;
	  alpha_max = alpha_Cell;
	}
	as->set_is_on_chull(false);
	as->set_alpha_mid(alpha_mid);
	as->set_alpha_max(alpha_max);
	alpha_mid_facet_map.insert(typename
		                 Alpha_facet_map::value_type(alpha_mid, *fit));
	//alpha_max_facet_map.insert(typename
	//                       Alpha_facet_map::value_type(alpha_max, *fit));
      }
      else { // on the convex hull
	alpha_mid = !is_infinite(pCell) ? pCell->get_alpha() 
	                                : pNeighbor->get_alpha();
	alpha_mid_facet_map.insert(typename
		                 Alpha_facet_map::value_type(alpha_mid, *fit));
	as->set_alpha_mid(alpha_mid);
	as->set_is_on_chull(true);
      }

      //cross links
      pCell->set_facet_status(i, as);
      pNeighbor->set_facet_status(iNeigh,as);
    }
  }
    
  // initialize alpha_min if mode GENERAL 
  if(get_mode() == GENERAL) {
    alpha_min_facet_map.clear();
    NT alpha_min;
    for( fit = finite_facets_begin(); 
	 fit != finite_facets_end(); ++fit)   {
      as = fit->first->get_facet_status(fit->second);
      if (is_Gabriel(*fit)) {
	as->set_is_Gabriel(true);
	alpha_min = squared_radius(*fit);
	as->set_alpha_min(alpha_min);
	alpha_min_facet_map.insert(typename
		                 Alpha_facet_map::value_type(alpha_min, *fit));
      }
      else as->set_is_Gabriel(false);
    }
  }
  return;
 }

template <class Dt>
void 
Alpha_shape_3<Dt>::initialize_alpha_edge_maps(bool reinitialize)
{
  Finite_edges_iterator eit;
  Facet_circulator fcirc, done;
  Alpha_status *asf, *as;
  NT alpha;

  if(! reinitialize){
    for (eit = finite_edges_begin(); 
	 eit != finite_edges_end(); ++eit) {

      as = new Alpha_status();
      fcirc = incident_facets(eit->first, eit->second, eit->third);
      while (is_infinite(*fcirc) ) ++fcirc; //skip infinite incident faces
      done = fcirc;
      as->set_is_on_chull(false);
      asf = (*fcirc).first->get_facet_status((*fcirc).second);
      as->set_alpha_mid(asf->alpha_mid());  // initialise as->alpha_mid
      as->set_alpha_max(asf->alpha_mid()); // and as->alpha->max to the same
     
      do {
	if (!is_infinite(*fcirc)) {
	  asf = (*fcirc).first->get_facet_status((*fcirc).second);
	  if (get_mode() == GENERAL && asf->is_Gabriel()) 
	    alpha = asf->alpha_min();
	  else alpha = asf->alpha_mid();
	  if (alpha < as->alpha_mid())  as->set_alpha_mid(alpha) ;
	  if ( ! asf->is_on_chull()) { 
	    if( as->alpha_max() <  asf->alpha_max())
	    as->set_alpha_max( asf->alpha_max());
	  }
	  else{
	    as->set_is_on_chull(true);
	  }
	}
      } while (++fcirc != done);
 
      alpha_mid_edge_map.insert(typename Alpha_edge_map::value_type
	                        (as->alpha_mid(), *eit));
      
      //cross links
      Vertex_handle_pair vhp = make_vertex_handle_pair(
				  eit->first->vertex(eit->second),
				  eit->first->vertex(eit->third));
      edge_alpha_map.insert(std::make_pair(vhp, as));
      } 
  }

  // initialize alphamin in GENERAL mode
  if (get_mode() == GENERAL) {
    alpha_min_edge_map.clear();
    for (eit = finite_edges_begin(); 
	 eit != finite_edges_end(); ++eit) {
      Vertex_handle_pair vhp =  make_vertex_handle_pair(
				  eit->first->vertex(eit->second),
				  eit->first->vertex(eit->third));
      as = edge_alpha_map[vhp];
      if (is_Gabriel(*eit)) {
	alpha = squared_radius(*eit);
	as->set_is_Gabriel(true);
	as->set_alpha_min(alpha);
	alpha_min_edge_map.insert(typename
		                  Alpha_edge_map::value_type(alpha,*eit));
      }
      else as->set_is_Gabriel(false);
    }
  }

  // alpha_mid has to be recomputed in case of reinitialisation
  if(reinitialize) {
    alpha_mid_edge_map.clear();
    for (eit = finite_edges_begin(); 
     eit != finite_edges_end(); ++eit) {
      Vertex_handle_pair vhp =  make_vertex_handle_pair(
				  eit->first->vertex(eit->second),
				  eit->first->vertex(eit->third));
      as = edge_alpha_map[vhp];
      fcirc = incident_facets(eit->first, eit->second, eit->third);
      while (is_infinite(*fcirc) ) ++fcirc; //skip infinite incident faces
      done = fcirc;
      do {
	if (!is_infinite(*fcirc)) {
	  asf = (*fcirc).first->get_facet_status((*fcirc).second);
	  if (asf->is_Gabriel() &&  asf->alpha_min() < as->alpha_mid()) 
	    as->set_alpha_mid(asf->alpha_min());
	}
      } while (++fcirc != done);
      alpha_mid_edge_map.insert(typename Alpha_edge_map::value_type
	                        (as->alpha_mid(), *eit));
    }
  }
  return;
}

template <class Dt>
void 
Alpha_shape_3<Dt>::initialize_alpha_vertex_maps(bool reinitialize)
{
  NT alpha;

  if (! reinitialize) {
    Finite_vertices_iterator vit;

    for( vit = finite_vertices_begin(); 
	 vit != finite_vertices_end();	 ++vit) {
 
      Alpha_status * as = new Alpha_status();
      as->set_is_on_chull(false);
      
      
      // set alpha_mid and alpha_max and is_on_chull
      std::list<Cell_handle> incidents;
      incident_cells(static_cast<Vertex_handle>(vit),
	             back_inserter(incidents));
      typename std::list<Cell_handle>::iterator
	chit=incidents.begin();
      if (is_infinite(*chit)) as->set_is_on_chull(true);
      while (is_infinite(*chit)) ++chit;
      alpha = (*chit)->get_alpha();
      as->set_alpha_mid(alpha);
      as->set_alpha_max(alpha);
      for( ; chit != incidents.end(); ++chit) {
	if (is_infinite(*chit)) as->set_is_on_chull(true);
	else {
	  alpha = (*chit)->get_alpha();
	  if (alpha < as->alpha_mid()) as->set_alpha_mid(alpha);
	  if (alpha > as->alpha_max()) as->set_alpha_max(alpha);
	}
      }

      alpha_mid_vertex_map.insert(typename Alpha_vertex_map::value_type
	                          (as->alpha_mid(),vit));
  
      // cross link
      vit->set_alpha_status(as);
    }
  }

  // set alpha_min in case GENERAL - reset in case of reinitialize
  if (get_mode() == GENERAL) {
    Vertex_handle_pair vhp;
    Alpha_status *asedge, *as;
    NT alpha_min;

    for( Finite_vertices_iterator vit = finite_vertices_begin(); 
	 vit != finite_vertices_end();  ++vit) {

      as = vit->get_alpha_status();
      as->set_alpha_min(as->alpha_mid());

      std::list<Vertex_handle> incidentv;
      incident_vertices(static_cast<Vertex_handle>(vit),
	                back_inserter(incidentv));
      typename std::list<Vertex_handle>::iterator vvit=incidentv.begin();
      for( ; vvit != incidentv.end(); ++vvit){
	if (!is_infinite(*vvit)) {
	  vhp = make_vertex_handle_pair( *vvit, vit);
	  asedge = edge_alpha_map[vhp];
	  alpha_min = asedge->is_Gabriel() ? asedge->alpha_min()
	                                   : asedge->alpha_mid();
	  if ( alpha_min < as->alpha_min())   as->set_alpha_min(alpha_min);
	}
      }
    }
  }
  return;
}

//---------------------------------------------------------------------

template <class Dt>
void 
Alpha_shape_3<Dt>::initialize_alpha_spectrum()
// merges the alpha values of alpha_cell_map 
// and alpha_min_facet_map alpha_min_edge_map in GENERAL mode
// only alpha_cell_map in REGULARIZED mode
{
  typename Alpha_cell_map::iterator cit ;
  typename Alpha_facet_map::iterator fit ;
  typename Alpha_edge_map::iterator eit ;
  alpha_spectrum.clear();

  if (get_mode() == GENERAL) {
    cit = alpha_cell_map.begin();
    fit = alpha_min_facet_map.begin();
    eit = alpha_min_edge_map.begin();
    alpha_spectrum.reserve(alpha_cell_map.size() +
			   alpha_min_facet_map.size() +
			   alpha_min_edge_map.size() );
  }
  else {
    alpha_spectrum.reserve(alpha_cell_map.size());
    cit = alpha_cell_map.begin();
    fit = alpha_min_facet_map.end();
    eit = alpha_min_edge_map.end();
  }


  while (cit != alpha_cell_map.end() ||
	 fit != alpha_min_facet_map.end() ||
	 eit != alpha_min_edge_map.end() ) {

    if ( cit != alpha_cell_map.end() 
	 && ( fit == alpha_min_facet_map.end() || !(fit->first < cit->first) )
	 && ( eit == alpha_min_edge_map.end() || !(eit->first < cit->first) )
	 ) {      //advance on cit
         if (alpha_spectrum.empty() ||  alpha_spectrum.back() < cit->first){
	    alpha_spectrum.push_back(cit->first);

      }
      cit++;
      
    }

    if ( fit != alpha_min_facet_map.end() 
	 && ( cit == alpha_cell_map.end() || !(cit->first < fit->first) )
	 && ( eit == alpha_min_edge_map.end() || !(eit->first < fit->first) )
	 ) {      //advance on fit
      if (alpha_spectrum.empty() ||  alpha_spectrum.back() < fit->first){
	    alpha_spectrum.push_back(fit->first);
      }
      fit++;
    }

    if ( eit != alpha_min_edge_map.end() 
	 && ( fit == alpha_min_facet_map.end() || !(fit->first < eit->first) )
	 && ( cit == alpha_cell_map.end() || !(cit->first < eit->first) )
	 ) {      //advance on eit
         if (alpha_spectrum.empty() ||  alpha_spectrum.back() <  eit->first) {
	    alpha_spectrum.push_back(eit->first);
      }
      eit++;
    }
  }
}
  


//---------------------------------------------------------------------


#if 0
// Obviously not ready yet
template <class Dt>
std::istream& operator>>(std::istream& is,  const Alpha_shape_3<Dt>& A)
  // Reads a alpha shape from stream `is' and assigns it to
  // Unknown creationvariable. Precondition: The extract operator must
  // be defined for `Point'.
{}
#endif

//---------------------------------------------------------------------

template <class Dt>
std::ostream& operator<<(std::ostream& os,  const Alpha_shape_3<Dt>& A)
  // Inserts the alpha shape into the stream `os' as an indexed face set. 
  // Precondition: The insert operator must be defined for `Point'
{
  typedef Alpha_shape_3<Dt>                  AS;
  typedef typename AS::Vertex_handle         Vertex_handle;
  typedef typename AS::Cell_handle           Cell_handle;
  typedef typename AS::Alpha_shape_vertices_iterator 
                                             Alpha_shape_vertices_iterator;
  typedef typename AS::Alpha_shape_facets_iterator
                                             Alpha_shape_facets_iterator;

  Unique_hash_map< Vertex_handle, int > V;
  int number_of_vertices = 0;

  Alpha_shape_vertices_iterator vit;
  for( vit = A.alpha_shape_vertices_begin();
       vit != A.alpha_shape_vertices_end();
       ++vit) {
    V[Vertex_handle(vit)] = number_of_vertices++;
    os << vit->point() << std::endl;
  }

  Cell_handle c;
  int i;
  Alpha_shape_facets_iterator fit;
  for( fit = A.alpha_shape_faces_begin();
       fit != A.alpha_shape_faces_end();
       ++fit) {
    c = fit->first;
    i = fit->second;
    // the following ensures that regulat facets are output
    // in ccw order
    if (A.classify(fit) == AS::REGULAR && A.classify(c) == AS::EXTERIOR){
      c = c->neighbor(i);
      i = c->index(fit->first);
    }
    int i0=(i+1)&3, i1=(i+2)&3, i2=(i+3)&3;
    os << V[c->vertex(i0)] << ' ' 
       << V[c->vertex(i1)] << ' ' 
       << V[c->vertex(i2)] << std::endl;
  }
  return os;
}

//---------------------------------------------------------------------

template <class Dt>
void
Alpha_shape_3<Dt>::update_alpha_shape_vertex_list() const
{
  alpha_shape_vertices_list.clear();
  use_vertex_cache = true;

  // write super regular vertices
  // alpha must be stictly larger than alpha_mid (and smaller than alpha_max)
  typename Alpha_vertex_map::const_iterator vit;
  for(vit = alpha_mid_vertex_map.begin();
      vit != alpha_mid_vertex_map.lower_bound(get_alpha());
      ++vit) {
     if (classify(vit->second) == SUPER_REGULAR)
	  alpha_shape_vertices_list.push_back(vit->second);
  }

  if (get_mode() == REGULARIZED) return;
  // get REGULAR and SINGULAR vertices
   for(;
       vit != alpha_mid_vertex_map.end();
       ++vit) {
     Classification_type ct = classify(vit->second);
     if (ct == SINGULAR || ct == REGULAR)
       alpha_shape_vertices_list.push_back(vit->second);
   }
   return;
}
	 

//---------------------------------------------------------------------

template <class Dt>
void
Alpha_shape_3<Dt>::update_alpha_shape_facet_list() const
{
  alpha_shape_facets_list.clear();
  use_facet_cache = true;
  // Writes the faces of the alpha shape `A' for the current 'alpha'-value
  // to the container where 'out' refers to.

  // Get regular facets
  // alpha must be stictly larger than alpha_mid and smaller than alpha_max
  typename Alpha_facet_map::const_iterator fit;
  for(  fit = alpha_mid_facet_map.begin();
	fit != alpha_mid_facet_map.lower_bound(get_alpha());
	++fit) {
    if (classify(fit->second) == REGULAR)
	  alpha_shape_facets_list.push_back(fit->second);
  }
  if (get_mode() == REGULARIZED) return;

  // Get singular facets
  // alpha must be stricty larger than alpha_min and smaller than
  // alpha_max
  typename Alpha_facet_map::const_iterator fit2;
  for( fit2 = alpha_min_facet_map.begin();
       fit2 != alpha_min_facet_map.lower_bound(get_alpha());
       ++fit2) {
    if (classify(fit2->second) == SINGULAR)
	  alpha_shape_facets_list.push_back(fit->second);
  }
  return;
}



//---------------------------------------------------------------------

template < class Dt >
typename Alpha_shape_3<Dt>::Classification_type  
Alpha_shape_3<Dt>::classify(const Alpha_status* as,
			    const NT& alpha) const
{
 //tetrahedra with circumragius=alpha are considered outside
  if ( !as->is_on_chull() && alpha > as->alpha_max()) return INTERIOR;
  else if ( alpha > as->alpha_mid()) return REGULAR;
  else if ( get_mode() == GENERAL && 
	    as->is_Gabriel() &&
	    alpha > as->alpha_min()) return SINGULAR;
  else return EXTERIOR;
}

template < class Dt >
typename Alpha_shape_3<Dt>::Classification_type  
Alpha_shape_3<Dt>::classify(const Cell_handle& s, 
			    int i,
			    const NT& alpha) const
  // Classifies the face `f' of the underlying Delaunay
  // tetrahedralization with respect to `A'.
{ 
  if (is_infinite(s,i))   return EXTERIOR;
  return classify(s->get_facet_status(i), alpha);
}
 

template < class Dt >
typename Alpha_shape_3<Dt>::Classification_type  
Alpha_shape_3<Dt>::classify(const Cell_handle& s, 
			    int i,
			    int j,
			    const NT& alpha) const
  // Classifies the edge `e' of the underlying Delaunay
  // tetrahedralization with respect to `A'.
{ 
  if (is_infinite(s, i, j))     return EXTERIOR;
  Vertex_handle_pair
    vhp=make_vertex_handle_pair(s->vertex(i),s->vertex(j));
  return classify(edge_alpha_map.find(vhp)->second, alpha);
}

//---------------------------------------------------------------------

template < class Dt >
typename Alpha_shape_3<Dt>::Classification_type  
Alpha_shape_3<Dt>::classify(const Vertex_handle& v,
			    const NT& alpha) const
  // Classifies the vertex `v' of the underlying Delaunay
  // tetrahedralization with respect to `A'.
{
  if (is_infinite(v))     return EXTERIOR;
  const Alpha_status* as = v->get_alpha_status();

  if ( !as->is_on_chull() && alpha > as->alpha_max()) return INTERIOR;
  else if ( alpha > as->alpha_mid()) return SUPER_REGULAR;
  else if ( get_mode() == GENERAL && 
	    alpha > as->alpha_min()) return REGULAR;
    else return SINGULAR;
}

//--------------------- NB COMPONENTS ---------------------------------

template < class Dt >
int
Alpha_shape_3<Dt>::number_of_solid_components(const NT& alpha) const
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
				 const NT alpha) const
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
  // (1) nb_components solid components <= nb_components
  // (2) all data points on the boundary or in its interior
{
  NT alpha = find_alpha_solid();
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
typename Alpha_shape_3<Dt>::NT 
Alpha_shape_3<Dt>::find_alpha_solid() const
  // compute the minumum alpha such that all data points 
  // are either on the boundary or in the interior
  // not necessarily connected
{
//   NT  alpha_solid = *(alpha_spectrum.begin());
//   Finite_vertices_iterator vit, done;
//   for( vit = finite_vertices_begin(); 
//        vit != finite_vertices_end(); ++vit) {
//     alpha_solid = max(alpha_solid, vit->get_alpha_status()->alpha_mid());
//   }
//   return alpha_solid;

  return alpha_mid_vertex_map.rbegin()->first;
}

// TO  DEBUG

template <class Dt>
void 
Alpha_shape_3<Dt>::print_maps()
{
  typename Alpha_cell_map::iterator cit ;
  typename Alpha_facet_map::iterator fit ;
  typename Alpha_edge_map::iterator eit ;

  std::cerr << "size of cell_map " << alpha_cell_map.size() 
	    <<   std::endl;
  std::cerr << "size of facet_map " << alpha_min_facet_map.size() <<
    std::endl;
  std::cerr << "size of edge_map " << alpha_min_edge_map.size() <<
    std::endl;
  std::cerr << std::endl;
  std::cerr << "alpha_cell_map " << std::endl;
  for(cit = alpha_cell_map.begin();
      cit != alpha_cell_map.end(); ++cit) {
    std::cerr << cit->first << std::endl;
  }
  std::cerr << std::endl;
  std::cerr << "alpha_min_facet_map " << std::endl;
  for(fit = alpha_min_facet_map.begin();
      fit != alpha_min_facet_map.end(); ++fit) {
    std::cerr << fit->first << std::endl;
  }
  std::cerr << std::endl;
  std::cerr << "alpha_min_edge_map " << std::endl;
  for(eit = alpha_min_edge_map.begin();
      eit != alpha_min_edge_map.end(); ++eit) {
    std::cerr << eit->first << std::endl;
  }
  std::cerr << std::endl;
}


template <class Dt>
void 
Alpha_shape_3<Dt>::print_alphas() const
{
  std::cerr << std::endl;
  std::cerr << " alpha values of facets" << std::endl;
  for(Finite_facets_iterator fit = finite_facets_begin();
      fit != finite_facets_end();
      ++fit) {
    const Alpha_status* as = fit->first->get_facet_status(fit->second);
    print_alpha_status(as);
  }
  std::cerr << std::endl;
  std::cerr << " alpha values of edges " << std::endl;
  for(Finite_edges_iterator eit = finite_edges_begin();
      eit != finite_edges_end();
      ++eit) {
    Vertex_handle_pair vhp = make_vertex_handle_pair(
				  eit->first->vertex(eit->second),
				  eit->first->vertex(eit->third));
    const Alpha_status* as = edge_alpha_map.find(vhp)->second;
    print_alpha_status(as);
  }
  std::cerr << std::endl;
  std::cerr << " alpha values of vertices " << std::endl;
  for(Finite_vertices_iterator vit = finite_vertices_begin();
      vit != finite_vertices_end();
      ++vit) {
     const Alpha_status* as = vit->get_alpha_status();
    print_alpha_status_vertex(as);
  }

}

template <class Dt>
void 
Alpha_shape_3<Dt>::print_alpha_status(const Alpha_status* as) const
{
  if ( get_mode() == GENERAL &&  as->is_Gabriel())
  std::cerr << as->alpha_min() ;
  else std::cerr <<  "---   " ;
  std::cerr << "\t";
  std::cerr <<  as->alpha_mid()  << "\t";
  if(as->is_on_chull()) std::cerr <<  "---   ";
  else   std::cerr << as->alpha_max();
  std::cerr << std::endl;
}

template <class Dt>
void 
Alpha_shape_3<Dt>::print_alpha_status_vertex(const Alpha_status* as) const
{
  if ( get_mode() == GENERAL)  std::cerr << as->alpha_min() ;
  else std::cerr <<  "---   " ;
  std::cerr << "\t";
  std::cerr <<  as->alpha_mid()  << "\t";
  if(as->is_on_chull()) std::cerr <<  "---   ";
  else   std::cerr << as->alpha_max();
  std::cerr << std::endl;
}

CGAL_END_NAMESPACE

#ifdef CGAL_USE_GEOMVIEW
#include <CGAL/IO/alpha_shape_geomview_ostream_3.h>
#endif

#endif //CGAL_ALPHA_SHAPE_3_H
