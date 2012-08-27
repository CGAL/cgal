// Copyright (c) 1997, 2012  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the so
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Tran Kai Frank DA <Frank.Da@sophia.inria.fr>
//                 Andreas Fabri <Andreas.Fabri@geometryfactory.com>
//                 Mariette Yvinec <Mariette.Yvinec@sophia.inria.fr>

#ifndef CGAL_ALPHA_SHAPE_3_H
#define CGAL_ALPHA_SHAPE_3_H

#include <CGAL/basic.h>

#include <set>
#include <map>
#include <list>
#include <vector>
#include <algorithm>
#include <utility>
#include <iostream>

#include <CGAL/Triangulation_utils_3.h>
#include <CGAL/Object.h>
#include <CGAL/Unique_hash_map.h>
#include <CGAL/Compact_container.h>
#include <CGAL/Alpha_shape_vertex_base_3.h>
#include <CGAL/Alpha_shape_cell_base_3.h>
#include <CGAL/internal/Lazy_alpha_nt_3.h>
#ifdef CGAL_USE_GEOMVIEW
#include <CGAL/IO/Geomview_stream.h>  // TBC
#endif

//-------------------------------------------------------------------
namespace CGAL {
//-------------------------------------------------------------------

template < class Dt, class ExactAlphaComparisonTag = Tag_false >
class Alpha_shape_3 : public Dt
{
  // DEFINITION The class Alpha_shape_3<Dt> represents the family
  // of alpha-shapes for a set of points (or a set of weighted points)
  // for all possible values of alpha. The alphashape is defined  through
  // the Delaunay tetrahedralization of the points
  // (or the Regular tetrahedralization in case of weighted points)
  // and depends on the value of a parameter called alpha.
  // The alpha_shape is the domain of a subcomplex of this triangulation
  // called the Alpha_complex. The alpha_complex includes any simplex
  // having  a circumscribing sphere (an orthogonal sphere
  // in case of weighted points) empty of other points
  // (or suborthogonal to other sites in case of weighted points)
  // with squared radius equal or less than alpha
 
  // The alpha_shapes comes in two versions : GENERAL or REGULARIZED
  // where the REGULARIZED version is onbtaining by restricting the
  // alpha complex ti is pure 3D component.

  // The cells of the triangulation are classified as INTERIOR
  // or EXTERIOR according to the value alpha_cell of their circumsphere 
  // squared radius compared to alpha.

  // In GENERAL mode each k-dimensional simplex of the triangulation
  // for (k=0,1,2) 
  // can be classified as EXTERIOR, SINGULAR, REGULAR
  // or INTERIOR with respect to the alpha shape.
  // In GENERAL mode a $k$ simplex is REGULAR if it is on the boundary
  // of the alpha_complex and belongs to a $k+1$ simplex in the complex
  // and it is SINGULAR simplex if it is  a boundary simplex tht is not
  // included in a $k+1$ simplex of the complex.
  
  // In REGULARIZED mode each k-dimensional simplex of the triangulation
  // for (k=0,1,2) 
  // can be classified as EXTERIOR, REGULAR
  // or INTERIOR with respect to the alpha shape.
  // A $k$ simplex is REGULAR if it is on the boundary of alpha complex
  // and belong to a tetrahedral cell of the complex.

  // Roughly, the Alpha_shapes data structure computes and stores, 
  // for each simplex
  // the at most three critical value (alpha_min, alpha_mid and alpha_max)
  // which compared to the actual alpha value
  // determine the classification of the simplex.


  //------------------------- TYPES ------------------------------------

public:
  typedef Dt                                        Triangulation;
  typedef typename Dt::Geom_traits                  Gt;
  typedef typename Dt::Triangulation_data_structure Tds;

  //extra the type used for representing alpha according to ExactAlphaComparisonTag
  typedef typename internal::Alpha_nt_selector_3<Gt,ExactAlphaComparisonTag,typename Dt::Weighted_tag>::Type_of_alpha  NT;
  typedef typename internal::Alpha_nt_selector_3<Gt,ExactAlphaComparisonTag,typename Dt::Weighted_tag>::Compute_squared_radius_3 Compute_squared_radius_3;
  typedef NT      FT;
  typedef typename Gt::FT Coord_type;
  //checks whether tags are correctly set in Vertex and Cell classes
  CGAL_static_assertion( (boost::is_same<NT,typename Dt::Cell::NT>::value) );
  CGAL_static_assertion( (boost::is_same<NT,typename Dt::Vertex::Alpha_status::NT>::value) );

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

  typedef typename Dt::size_type    size_type;
  typedef typename Dt::Locate_type  Locate_type;
  typedef typename Dt::Weighted_tag Weighted_tag;

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
  using Dt::is_infinite;
  using Dt::is_Gabriel;
  using Dt::incident_cells;
  using Dt::incident_vertices;
  using Dt::incident_facets;
  using Dt::locate;

  enum Classification_type {EXTERIOR, 
			    SINGULAR, 
			    REGULAR,
			    INTERIOR};
 
  enum Mode {GENERAL, REGULARIZED};


  typedef CGAL::Alpha_status< NT >          Alpha_status;
  typedef Compact_container<Alpha_status>   Alpha_status_container;
  typedef typename Alpha_status_container::const_iterator 
                                            Alpha_status_const_iterator;
  typedef typename Alpha_status_container::iterator 
                                            Alpha_status_iterator;
  typedef std::vector< NT >                 Alpha_spectrum;

  typedef std::multimap< NT, Cell_handle >  Alpha_cell_map;
  typedef std::multimap< NT, Facet>         Alpha_facet_map;
  typedef std::multimap< NT, Edge >         Alpha_edge_map;
  typedef std::multimap< NT, Vertex_handle> Alpha_vertex_map;

  typedef std::pair<Vertex_handle, Vertex_handle> Vertex_handle_pair;
  typedef std::map<Vertex_handle_pair,Alpha_status_iterator> Edge_alpha_map;

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
  NT _alpha_solid;
  Mode _mode;
  mutable bool use_vertex_cache;
  mutable bool use_facet_cache;

  // only finite facets and simplices are inserted into the maps
  Alpha_cell_map     alpha_cell_map;
  Alpha_facet_map    alpha_min_facet_map;
  Alpha_edge_map     alpha_min_edge_map;
  Alpha_vertex_map   alpha_min_vertex_map;

  Alpha_spectrum             alpha_spectrum;
  Alpha_status_container     alpha_status_container;

  Edge_alpha_map edge_alpha_map;

  //deprecated - for backward compatibility
  mutable std::list< Vertex_handle > alpha_shape_vertices_list;
  mutable std::list< Facet > alpha_shape_facets_list;


  //------------------------- CONSTRUCTORS ------------------------------
public:
  // Introduces an empty alpha-shape `A' for a 
  // alpha-value `alpha'. 
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
 
  // Introduces an alpha-shape `A' for the alpha-value
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
  std::ptrdiff_t make_alpha_shape(const InputIterator& first, 
                                  const InputIterator& last)
    {
      clear();
      size_type n = Dt::insert(first, last);
      if (dimension() == 3){
        initialize_alpha();
      }
      return n;
    }

  // Introduces an alpha-shape `A' 
  // that is initialized with the points in the range
  // from first to last

private :

  //--------------------- INITIALIZATION OF PRIVATE MEMBERS -----------
  // called with reinitialize=false on first initialization
  // reinitialize=true when switching the mode.
  void initialize_alpha_cell_map();
  void initialize_alpha_facet_maps(bool reinitialize = false);
  void initialize_alpha_edge_maps(bool reinitialize = false);
  void initialize_alpha_vertex_maps(bool reinitialize = false);
  void initialize_alpha_spectrum();
  void initialize_alpha(bool reinitialize = false) {
    if (reinitialize == false) initialize_alpha_cell_map();
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

 // the version to be used with Tag_true is templated to avoid
 // instanciation through explicit instantiation of the whole class
  void set_alpha_min_of_vertices(Tag_false) 
  { 
    for( Finite_vertices_iterator vit = finite_vertices_begin(); 
	 vit != finite_vertices_end();  ++vit){
      Alpha_status* as = vit->get_alpha_status();
      as->set_is_Gabriel(true);  
      as->set_alpha_min(NT(0));
    }
    // insert a single vertex into the map because they all have the 
    // same alpha_min value
    alpha_min_vertex_map.insert(typename Alpha_vertex_map::value_type
				( NT(0), finite_vertices_begin()));
  }
  template <class Tag>
  void set_alpha_min_of_vertices(Tag) 
  {
    for( Finite_vertices_iterator vit = finite_vertices_begin(); 
	 vit != finite_vertices_end();  ++vit) {
      if (is_Gabriel(vit)) {
	Alpha_status* as = vit->get_alpha_status();
	as->set_is_Gabriel(true);  
	as->set_alpha_min(squared_radius(vit));      
	alpha_min_vertex_map.insert(typename Alpha_vertex_map::value_type
				    (as->alpha_min(),vit));
      }
    }
    return;
  }



  //---------------------------------------------------------------------

public:

  void clear()
    {
      // clears the structure
      alpha_status_container.clear();
      Dt::clear();

      alpha_cell_map.clear();
      alpha_min_facet_map.clear();
      alpha_min_edge_map.clear();
      alpha_min_vertex_map.clear();
   
      alpha_spectrum.clear();

      alpha_shape_vertices_list.clear();
      alpha_shape_facets_list.clear();

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
      CGAL_triangulation_assertion( n > 0 && 
		      n <= static_cast<int>(alpha_spectrum.size()) );
      return alpha_spectrum[n-1];
    }
  
  size_type number_of_alphas() const
    // Returns the number of different alpha-values
    {
      return alpha_spectrum.size();
    }

  const Edge_alpha_map* get_edge_alpha_map() const
  {
     return  &edge_alpha_map;
  }    
    
  //---------------------------------------------------------------------

private:

  // the dynamic version is not yet implemented
  // desactivate the tetrahedralization member functions
  void insert(const Point& /*p*/) {}
  // Inserts point `p' in the alpha shape and returns the
  // corresponding vertex of the underlying Delaunay tetrahedralization.
  // If point `p' coincides with an already existing vertex, this
  // vertex is returned and the alpha shape remains unchanged.
  // Otherwise, the vertex is inserted in the underlying Delaunay
  // tetrahedralization and the associated intervals are updated.

  void remove(Vertex_handle /*v*/) {}
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

  void  update_alpha_shape_vertex_list() const;
  void  update_alpha_shape_facet_list() const; 

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
      return CGAL::filter_iterator(finite_cells_end(),
			           Exterior_cell_test(this),
			           finite_cells_begin());
    }
  
  Alpha_shape_cells_iterator alpha_shape_cells_end() const
    {
      return CGAL::filter_iterator(finite_cells_end(),
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
public:
  void compute_edge_status( const Cell_handle&  c, 
			    int i, 
			    int j, 
			    Alpha_status& as) const;

  Classification_type classify(const Alpha_status& as, const NT& alpha) const;
  Classification_type classify(const Alpha_status* as, const NT& alpha) const;
  Classification_type classify(const Alpha_status_const_iterator as, 
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
      // s->radius == alpha => f interior
    {
      if (is_infinite(s)) return EXTERIOR;
      return (s->get_alpha() <=  alpha) ? INTERIOR : EXTERIOR;
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
  size_type
  number_solid_components() const
    {
      return number_of_solid_components(get_alpha());
    }

  size_type
  number_of_solid_components() const
    {
      return number_of_solid_components(get_alpha());
    }

  size_type
  number_solid_components(const NT& alpha) const
    {
      return number_of_solid_components(alpha);
    }

  size_type
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

  Alpha_iterator find_optimal_alpha(size_type nb_components) const;
  // find the minimum alpha that satisfies the properties
  // (1) all data points are on the boundary of some 3d component
  //    or in its interior
  // (2) the nb of solid components is equal or less than nb_component
  
  NT find_alpha_solid() const;
  // compute the minumum alpha such that all data points 
  // are either on the boundary or in the interior
  // not necessarily connected
  // starting point for searching 
  // takes O(#alpha_shape) time


  //------------------- GEOMETRIC PRIMITIVES ----------------------------
private:
  NT squared_radius(const Cell_handle& s) const
    {
      return Compute_squared_radius_3()(*this)(
	  this->point(s,0), this->point(s,1),
	  this->point(s,2), this->point(s,3));
    }

  NT squared_radius(const Cell_handle& s, const int& i) const
    {
      return Compute_squared_radius_3()(*this) (
	  this->point(s,vertex_triple_index(i,0)),
	  this->point(s,vertex_triple_index(i,1)),
	  this->point(s,vertex_triple_index(i,2)) );
    }

  NT squared_radius(const Facet& f) const {
    return squared_radius(f.first, f.second);
  }

  NT squared_radius(const Cell_handle& s, 
			    const int& i, const int& j) const
    {
      return Compute_squared_radius_3()(*this)(
	  this->point(s,i), this->point(s,j));
    }

  NT squared_radius(const Edge& e) const {
   return  squared_radius(e.first,e.second,e.third);
  }

  NT squared_radius(const Vertex_handle& v) const {
    return  Compute_squared_radius_3()(*this)(v->point()); 
  }


  //---------------------------------------------------------------------

private:
  // prevent default copy constructor and default assigment
  Alpha_shape_3(const Alpha_shape_3&);
  void operator=(const Alpha_shape_3&);

  //---------------------------------------------------------------------
public:  
#ifdef CGAL_USE_GEOMVIEW
  void show_alpha_shape_faces(Geomview_stream &gv) const;
#endif


  // to Debug
  void print_maps() const; 
  void print_alphas() const;
  void print_alpha_status( const Alpha_status& as) const;
  

  // To extract the alpha_shape faces for a given alpha value



  template<class OutputIterator>
  OutputIterator get_alpha_shape_cells(OutputIterator it, 
				       Classification_type type,
				       const NT& alpha) const
  {
    Finite_cells_iterator cit = finite_cells_begin();
    for( ; cit != finite_cells_end() ; ++cit){
      if (classify(cit, alpha) == type) *it++ = Cell_handle(cit);
    }
    return it;
  }

  template<class OutputIterator>
  OutputIterator get_alpha_shape_facets(OutputIterator it, 
					Classification_type type,
					const NT& alpha) const
  {
    Finite_facets_iterator fit = finite_facets_begin();
    for( ; fit != finite_facets_end() ; ++fit){
      if (classify(*fit, alpha) == type) *it++ = *fit;
    }
    return it;
  }

  template<class OutputIterator>
  OutputIterator get_alpha_shape_edges(OutputIterator it, 
				       Classification_type type,
				       const NT& alpha) const
  {
    Finite_edges_iterator eit = finite_edges_begin();
    for( ; eit != finite_edges_end() ; ++eit){
      if (classify(*eit, alpha) == type) *it++ = *eit;
    }
    return it;
  }

  template<class OutputIterator>
   OutputIterator get_alpha_shape_vertices(OutputIterator it, 
					   Classification_type type,
					   const NT& alpha) const
  {
    Finite_vertices_iterator vit = finite_vertices_begin();
    for( ; vit != finite_vertices_end() ; ++vit){
      if (classify(vit, alpha) == type) *it++ = Vertex_handle(vit);
    }
    return it;
  }

  template<class OutputIterator>
  OutputIterator get_alpha_shape_cells(OutputIterator it, 
				       Classification_type type) const
  { return get_alpha_shape_cells(it, type, get_alpha());}

  template<class OutputIterator>
  OutputIterator get_alpha_shape_facets(OutputIterator it, 
				       Classification_type type) const
  { return get_alpha_shape_facets(it, type, get_alpha());}

  template<class OutputIterator>
  OutputIterator get_alpha_shape_edges(OutputIterator it, 
				       Classification_type type) const
  { return get_alpha_shape_edges(it, type, get_alpha());}

  template<class OutputIterator>
  OutputIterator get_alpha_shape_vertices(OutputIterator it, 
				       Classification_type type) const
  { return get_alpha_shape_vertices(it, type, get_alpha());}

   template<class OutputIterator> 
   OutputIterator filtration(OutputIterator it)  const
   // scan  the  alpha_cell_map, alpha_min_facet_map,  alpha_min_edge_map  
   // and alpha_min_vertex in GENERAL mode 
   // only alpha_cell_map in REGULARIZED mode 
   // and output all the faces in order of alpha value of their appearing 
   // in the alpha complexe 
   { 
     typename Alpha_cell_map::const_iterator cit ;
     typename Alpha_facet_map::const_iterator fit ;
     typename Alpha_edge_map::const_iterator eit ;
     typename Alpha_vertex_map::const_iterator vit;
  
     if (get_mode() == GENERAL) {
       cit = alpha_cell_map.begin();
       fit = alpha_min_facet_map.begin();
       eit = alpha_min_edge_map.begin();
       vit = alpha_min_vertex_map.begin();
     }
     else { //mode==REGULARIZED do not scan maps of Gabriel elements
       cit = alpha_cell_map.begin();
       fit = alpha_min_facet_map.end();
       eit = alpha_min_edge_map.end();
       vit = alpha_min_vertex_map.end();
     }

     // sets to avoid multiple output of the same face 
     // as a regular subfaces of different faces 
     std::set<Facet>  facet_set; 
     std::set<Vertex_handle_pair>   edge_set; 
     std::set<Vertex_handle> vertex_set; 
     NT alpha_current = 0; 

     while (cit != alpha_cell_map.end()) { 

       if ( vit != alpha_min_vertex_map.end()  
 	   && (eit == alpha_min_edge_map.end() || (vit->first <= eit->first)) 
 	   && (fit == alpha_min_facet_map.end()|| (vit->first <= fit->first)) 
 	   && (cit == alpha_cell_map.end()     || (vit->first <= cit->first))) 
 	{ 
 	  //advance on vit 
 	  filtration_set_management(vit, alpha_current, 
 				    facet_set, edge_set, vertex_set); 
 	  filtration_output(vit->first, vit->second, it); 
 	  vit++; 
 	} 

       if ( eit != alpha_min_edge_map.end()  
 	 && ( fit == alpha_min_facet_map.end() || (eit->first <= fit->first) ) 
 	 && ( cit == alpha_cell_map.end()      || (eit->first <= cit->first) ) 
 	 && ( vit == alpha_min_vertex_map.end()|| (vit->first >  eit->first) ) 
 	 ) {      //advance on eit 
	 filtration_set_management(eit, alpha_current, 
 				    facet_set, edge_set, vertex_set); 
 	filtration_output(eit->first, eit->second, it, vertex_set); 
 	eit++; 
       } 

       if ( fit != alpha_min_facet_map.end()  
 	 && (cit == alpha_cell_map.end()      || (fit->first <= cit->first)) 
 	 && (eit == alpha_min_edge_map.end()  || (eit->first >  fit->first))  
 	 && (vit == alpha_min_vertex_map.end()|| (vit->first >  fit->first)) 
 	 ) {      //advance on fit 
	 filtration_set_management(fit, alpha_current, 
 				  facet_set, edge_set, vertex_set); 
	 filtration_output(fit->first, fit->second, it,  
			   edge_set, vertex_set); 
 	fit++; 
       } 

       if ( cit != alpha_cell_map.end()  
 	 && (fit == alpha_min_facet_map.end() || (fit->first > cit->first) ) 
 	 && (eit == alpha_min_edge_map.end()  || (eit->first > cit->first) ) 
 	 && (vit == alpha_min_vertex_map.end()|| (vit->first > cit->first) ) 
 	 ) {      //advance on cit 
	 filtration_set_management(cit, alpha_current, 
 				    facet_set, edge_set, vertex_set); 
	 filtration_output(cit->first, cit->second, it,
			   facet_set, edge_set, vertex_set); 
 	cit++; 
       } 
     } 
     return it; 
   } 

  private: 

   template<class Alpha_face_iterator> 
   void 
     filtration_set_management ( Alpha_face_iterator afit, 
 				NT& alpha_current, 
 				std::set<Facet>&  facet_set, 
 				std::set<Vertex_handle_pair>&   edge_set, 
 				std::set<Vertex_handle>& vertex_set)  const
   { 
     if (afit->first != alpha_current) { //new alpha_value 
       alpha_current = afit->first; 
       facet_set.clear(); 
       edge_set.clear(); 
       vertex_set.clear(); 
     } 
     return; 
   } 

   template<class OutputIterator> 
   OutputIterator   
   filtration_output( const NT & /*alpha*/,  
 		     Vertex_handle vh,  
 		     OutputIterator it,  
 		     Tag_true)   const 
   { 
     it++ = make_object(vh); 
     //std::cerr << "filtration " << alpha << " \t  VERTEX " << std::endl; 
     return it; 
   } 

   template<class OutputIterator> 
   OutputIterator   
   filtration_output( const NT& /*alpha*/,  
 		     Vertex_handle vh,  
 		     OutputIterator it,  
 		     Tag_false)     const 
   { 
     // when Delaunay, the alpha_min_vertex_map contains a single vertex 
     // because all vertices are Gabriel with the same alpha_min=0 
     // this affects only the GENERAL mode
     if (get_mode() == GENERAL){
       Finite_vertices_iterator vit=finite_vertices_begin(); 
       for( ; vit != finite_vertices_end(); vit++) { 
	 it++ = make_object( Vertex_handle(vit)); 
       } 
     }
     else {
       it++ = make_object(vh);
     }
     //std::cerr << "filtration " << alpha << " \t  VERTEX " << std::endl; 
     return it; 
   } 

   template<class OutputIterator> 
   OutputIterator   
   filtration_output( const NT& alpha,  
 		     Vertex_handle vh,  
 		     OutputIterator it) const 
   { 
     return filtration_output(alpha, vh, it, Weighted_tag()); 
   } 


  template<class OutputIterator> 
  OutputIterator   
  filtration_output( const NT& alpha,  
 		    Edge e,  
 		    OutputIterator it, 
 		    std::set<Vertex_handle>& vertex_set) const 
  { 
    Vertex_handle vh[] = {e.first->vertex(e.second),  
 			  e.first->vertex(e.third)}; 
    for(int i=0; i<2; i++) { 
      Alpha_status* as = vh[i]->get_alpha_status(); 
      if ( (get_mode()== REGULARIZED || !as->is_Gabriel())   
 	  && as->alpha_mid() == alpha  
 	  && vertex_set.find(vh[i]) == vertex_set.end() ) { 
        filtration_output( alpha, vh[i], it); 
        vertex_set.insert(vh[i]); 
      } 
    } 
    it++ = make_object(e); 
    //std::cerr << "filtration " << alpha << " \t EDGE " << std::endl; 
    return it; 
  } 
   
  template<class OutputIterator> 
  OutputIterator 
  filtration_output( const NT& alpha,  
 		    Facet f,  
 		    OutputIterator it, 
 		    std::set<Vertex_handle_pair>& edge_set, 
 		    std::set<Vertex_handle>& vertex_set ) const 
  { 
    Cell_handle c = f.first; 
    int facet_index = f.second; 

    for(int k=0; k<3; k++) { 
      int i = vertex_triple_index(facet_index, k ); 
      int j = vertex_triple_index(facet_index, this->ccw(k)); 
      Alpha_status as; 
      Vertex_handle_pair 
 	 vhp = make_vertex_handle_pair(c->vertex(i),c->vertex(j));

      if (get_mode() == GENERAL) { 
	as = *(edge_alpha_map.find(vhp)->second); 
      } 
      else{ //no edge map in REGULARIZED mode - classify on the fly 
	compute_edge_status( c, i, j, as); 
      } 
     
      if ( (get_mode()== REGULARIZED || !as.is_Gabriel())
	   && as.alpha_mid() == alpha  
	   && edge_set.find(vhp)== edge_set.end() ) {
	filtration_output( alpha, make_triple(c,i,j), it, vertex_set); 
        edge_set.insert(vhp); 
      } 
    } 

    it++ = make_object(f); 
    //std::cerr << "filtration " << alpha << " \t FACET " << std::endl; 
    return it; 
  } 

  template<class OutputIterator> 
  OutputIterator 
  filtration_output( const NT& alpha,  
 		    Cell_handle c,  
 		    OutputIterator it, 
 		    std::set<Facet>& facet_set, 
 		    std::set<Vertex_handle_pair>& edge_set, 
 		    std::set<Vertex_handle>& vertex_set) const 
  { 
    for(int i=0; i<4; i++) { 
      Alpha_status_iterator as = c->get_facet_status(i); 
      Facet f = std::make_pair(c,i); 
      if ((get_mode()== REGULARIZED || !as->is_Gabriel())
	   && as->alpha_mid() == alpha  
	   && facet_set.find(f) == facet_set.end() 
	   && facet_set.find(std::make_pair(c->neighbor(i), 
					    this->mirror_index(c, i)))
	      == facet_set.end()) { 
        filtration_output( alpha, f, it, edge_set, vertex_set); 
        facet_set.insert(f); 
      } 
    } 

    it++ = make_object(c); 
    //std::cerr << "filtration " << alpha << " \t CELL " << std::endl; 
    return it; 
  } 
 
  
};



//---------------------------------------------------------------------
//--------------------- MEMBER FUNCTIONS-------------------------------
//---------------------------------------------------------------------


//--------------------- INITIALIZATION OF PRIVATE MEMBERS -------------
  
template <class Dt,class EACT>
void 
Alpha_shape_3<Dt,EACT>::initialize_alpha_cell_map()
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

template <class Dt,class EACT>
void 
Alpha_shape_3<Dt,EACT>::initialize_alpha_facet_maps(bool reinitialize)
{
  Finite_facets_iterator fit;  
  Cell_handle pCell, pNeighbor ;
  int i, iNeigh;
  Alpha_status_iterator as;

  if (!reinitialize) {
    NT alpha_max, alpha_mid;
    for( fit = finite_facets_begin(); 
	 fit != finite_facets_end(); ++fit)   {

      as = alpha_status_container.insert(Alpha_status());
  
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
	//	alpha_mid_facet_map.insert(typename
	//	                 Alpha_facet_map::value_type(alpha_mid, *fit));
      }
      else { // on the convex hull
	alpha_mid = !is_infinite(pCell) ? pCell->get_alpha() 
	                                : pNeighbor->get_alpha();
	as->set_alpha_mid(alpha_mid);
	as->set_is_on_chull(true);
      }

      //cross links
      pCell->set_facet_status(i, as);
      pNeighbor->set_facet_status(iNeigh,as);
    }
  }
    
  // initialize alpha_min if mode GENERAL 
  if(get_mode() == GENERAL &&  alpha_min_facet_map.empty()) {
    //already done if !alpha_min_facet_map.empty()
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

template <class Dt,class EACT>
void 
Alpha_shape_3<Dt,EACT>::initialize_alpha_edge_maps(bool )
{
  // alpha_status for edges, edge_alpha_map 
  // and alpha_mid_edge and alpha_min_edge
  // are initialized only in GENERAL mode
  if(get_mode() == REGULARIZED) {return;} //no_edge_map in REGULARIZED mode
  if ( !edge_alpha_map.empty()) return; // already done

  Finite_edges_iterator eit;
  Alpha_status_iterator as;

  for (eit = finite_edges_begin(); 
       eit != finite_edges_end(); ++eit) {
    as = alpha_status_container.insert(Alpha_status());
    compute_edge_status(eit->first, eit->second, eit->third, *as);
    if ( as->is_Gabriel()) {
      alpha_min_edge_map.insert(typename
				Alpha_edge_map::value_type(as->alpha_min(),
							   *eit));
    }
     //cross links
    Vertex_handle_pair 
      vhp = make_vertex_handle_pair( eit->first->vertex(eit->second),
				     eit->first->vertex(eit->third));
    edge_alpha_map.insert(std::make_pair(vhp, as));
  }
  return;
}

template <class Dt,class EACT>
void 
Alpha_shape_3<Dt,EACT>::initialize_alpha_vertex_maps(bool reinitialize)
{
  //for a vertex 
  // alpha_max =  max of alpha values of incident cells
  // alpha_mid =  min of alpha values of incident cells in REGULAR mode
  //           =  min of alpha values of incidents faces in GENERAL mode
  // alpha_min = -squared_radius of weighted point, 
  //              if the vertex is Gabriel set only in GENERAL mode

  NT alpha, alpha_mid;
  Finite_vertices_iterator vit;

  if (reinitialize == false) _alpha_solid = alpha_cell_map.begin()->first;

  for( vit = finite_vertices_begin(); 
	 vit != finite_vertices_end();	 ++vit) {
    Alpha_status*  as = vit->get_alpha_status();

    if (reinitialize == false) {
      // set is_on_chull, compute alpha_max 
      // and alpha_mid (version REGULAR)
      // compute _alpha_solid (max of alpha_mid of vertices in REGULAR mode)
      as->set_is_on_chull(false);
      std::list<Cell_handle> incidents;
      incident_cells(static_cast<Vertex_handle>(vit),
	             back_inserter(incidents));
      typename std::list<Cell_handle>::iterator chit=incidents.begin();
      if (is_infinite(*chit)) as->set_is_on_chull(true);
      while (is_infinite(*chit)) ++chit; //skip infinte cells
      alpha = (*chit)->get_alpha();
      as->set_alpha_mid(alpha);
      as->set_alpha_max(alpha);
      ++chit;
      for( ; chit != incidents.end(); ++chit) {
	if (is_infinite(*chit)) as->set_is_on_chull(true);
	else {
	  alpha = (*chit)->get_alpha();
	  if (alpha < as->alpha_mid()) as->set_alpha_mid(alpha);
	  if (alpha > as->alpha_max()) as->set_alpha_max(alpha);
	}
      }
      if (as->alpha_mid() > _alpha_solid)  _alpha_solid = as->alpha_mid();
    }
  
    if (get_mode() == GENERAL) { //reset alpha_mid,  set alph_min
      std::list<Vertex_handle> incidentv;
      incident_vertices(static_cast<Vertex_handle>(vit),
			back_inserter(incidentv));
      typename std::list<Vertex_handle>::iterator vvit=incidentv.begin();
      for( ; vvit != incidentv.end(); ++vvit) {
	if (!is_infinite(*vvit)) {
	  Vertex_handle_pair vhp = make_vertex_handle_pair( *vvit, vit);
	  Alpha_status_iterator asedge = edge_alpha_map[vhp];
	  alpha_mid = asedge->is_Gabriel() ? asedge->alpha_min()
	    : asedge->alpha_mid();
	  if ( alpha_mid < as->alpha_mid()) as->set_alpha_mid(alpha_mid);
	}
      }
    }

    if (get_mode()== REGULARIZED && reinitialize == true) {
      // reset alpha_mid
      std::list<Cell_handle> incidents;
      incident_cells(static_cast<Vertex_handle>(vit),
		     back_inserter(incidents));
      typename std::list<Cell_handle>::iterator chit=incidents.begin();
      while (is_infinite(*chit)) ++chit; //skip infinte cells
      alpha = (*chit)->get_alpha();
      as->set_alpha_mid(alpha);
      for( ; chit != incidents.end(); ++chit) {
	if (is_infinite(*chit)) as->set_is_on_chull(true);
	else {
	  alpha = (*chit)->get_alpha();
	  if (alpha < as->alpha_mid()) as->set_alpha_mid(alpha);
	}
      }
    }
      
  }
 
  // set alpha_min in case GENERAL 
  if (get_mode() == GENERAL && alpha_min_vertex_map.empty()) {
    set_alpha_min_of_vertices(Weighted_tag());
  }
  return;
}



//---------------------------------------------------------------------

template <class Dt,class EACT>
void 
Alpha_shape_3<Dt,EACT>::initialize_alpha_spectrum()
// merges the alpha values of alpha_cell_map 
// and alpha_min_facet_map alpha_min_edge_map alpha_min_vertex in GENERAL mode
// only alpha_cell_map in REGULARIZED mode
{
  typename Alpha_cell_map::iterator cit ;
  typename Alpha_facet_map::iterator fit ;
  typename Alpha_edge_map::iterator eit ;
  typename Alpha_vertex_map::iterator vit;
  alpha_spectrum.clear();

  if (get_mode() == GENERAL) {
    cit = alpha_cell_map.begin();
    fit = alpha_min_facet_map.begin();
    eit = alpha_min_edge_map.begin();
    vit = alpha_min_vertex_map.begin();
    alpha_spectrum.reserve(alpha_cell_map.size() +
			   alpha_min_facet_map.size() +
			   alpha_min_edge_map.size() +
			   alpha_min_vertex_map.size());
  }
  else {
    alpha_spectrum.reserve(alpha_cell_map.size());
    cit = alpha_cell_map.begin();
    fit = alpha_min_facet_map.end();
    eit = alpha_min_edge_map.end();
    vit = alpha_min_vertex_map.end();
  }


  while (cit != alpha_cell_map.end() ||
	 fit != alpha_min_facet_map.end() ||
	 eit != alpha_min_edge_map.end() ) {

    if ( cit != alpha_cell_map.end() 
	 && ( fit == alpha_min_facet_map.end() || !(fit->first < cit->first) )
	 && ( eit == alpha_min_edge_map.end() || !(eit->first < cit->first) )
	 && ( vit == alpha_min_vertex_map.end() || !(vit->first < cit->first) )
	 ) {      //advance on cit
      if (alpha_spectrum.empty() ||  alpha_spectrum.back() < cit->first){
	alpha_spectrum.push_back(cit->first); 
      }
      cit++;
     }

    if ( fit != alpha_min_facet_map.end() 
	 && ( cit == alpha_cell_map.end() || !(cit->first < fit->first) )
	 && ( eit == alpha_min_edge_map.end() || !(eit->first < fit->first) ) 
	 && ( vit == alpha_min_vertex_map.end() || !(vit->first < fit->first) )
	 ) {      //advance on fit
      if (alpha_spectrum.empty() ||  alpha_spectrum.back() < fit->first){
	    alpha_spectrum.push_back(fit->first);
      }
      fit++;
    }

    if ( eit != alpha_min_edge_map.end() 
	 && ( fit == alpha_min_facet_map.end() || !(fit->first < eit->first) )
	 && ( cit == alpha_cell_map.end() || !(cit->first < eit->first) )
	 && ( vit == alpha_min_vertex_map.end() || !(vit->first < eit->first) )
	 ) {      //advance on eit
         if (alpha_spectrum.empty() ||  alpha_spectrum.back() <  eit->first) {
	    alpha_spectrum.push_back(eit->first);
      }
      eit++;
    }

    if ( vit != alpha_min_vertex_map.end() 
	 && ( fit == alpha_min_facet_map.end() || !(fit->first < vit->first) )
	 && ( cit == alpha_cell_map.end() || !(cit->first < vit->first) )
	 && ( eit == alpha_min_edge_map.end() || !(eit->first < vit->first) )
	 ) { //advance on vit
         if (alpha_spectrum.empty() ||  alpha_spectrum.back() <  vit->first) {
	    alpha_spectrum.push_back(vit->first);
      }
      vit++;
    }
  }
}
  


//---------------------------------------------------------------------


#if 0
// Obviously not ready yet
template <class Dt,class EACT>
std::istream& operator>>(std::istream& is,  const Alpha_shape_3<Dt,EACT>& A)
  // Reads a alpha shape from stream `is' and assigns it to
  // Unknown creationvariable. Precondition: The extract operator must
  // be defined for `Point'.
{}
#endif

//---------------------------------------------------------------------

template <class Dt,class EACT>
std::ostream& operator<<(std::ostream& os,  const Alpha_shape_3<Dt,EACT>& A)
  // Inserts the alpha shape into the stream `os' as an indexed face set. 
  // Precondition: The insert operator must be defined for `Point'
{
  typedef Alpha_shape_3<Dt,EACT>                  AS;
  typedef typename AS::size_type             size_type;
  typedef typename AS::Vertex_handle         Vertex_handle;
  typedef typename AS::Cell_handle           Cell_handle;
  typedef typename AS::Alpha_shape_vertices_iterator 
                                             Alpha_shape_vertices_iterator;
  typedef typename AS::Alpha_shape_facets_iterator
                                             Alpha_shape_facets_iterator;

  Unique_hash_map< Vertex_handle, size_type > V;
  size_type number_of_vertices = 0;

  Alpha_shape_vertices_iterator vit;
  for( vit = A.alpha_shape_vertices_begin();
       vit != A.alpha_shape_vertices_end();
       ++vit) {
    V[*vit] = number_of_vertices++;
    os << (*vit)->point() << std::endl;
  }

  Cell_handle c;
  int i;
  Alpha_shape_facets_iterator fit;
  for( fit = A.alpha_shape_facets_begin();
       fit != A.alpha_shape_facets_end();
       ++fit) {
    c = fit->first;
    i = fit->second;
    // the following ensures that regular facets are output
    // in ccw order
    if (A.classify(*fit) == AS::REGULAR && (A.classify(c) == AS::INTERIOR)){
      c = c->neighbor(i);
      i = c->index(fit->first);
    }
    int i0 = Triangulation_utils_3::vertex_triple_index(i,0);
    int i1 = Triangulation_utils_3::vertex_triple_index(i,1);
    int i2 = Triangulation_utils_3::vertex_triple_index(i,2);
    os << V[c->vertex(i0)] << ' ' 
       << V[c->vertex(i1)] << ' ' 
       << V[c->vertex(i2)] << std::endl;
  }
  return os;
}

//---------------------------------------------------------------------

template <class Dt,class EACT>
void
Alpha_shape_3<Dt,EACT>::update_alpha_shape_vertex_list() const
{
  alpha_shape_vertices_list.clear();
  use_vertex_cache = true;

  std::back_insert_iterator<std::list< Vertex_handle > >
    it = back_inserter(alpha_shape_vertices_list);

  get_alpha_shape_vertices(it, REGULAR);
  if (get_mode()==GENERAL) get_alpha_shape_vertices(it, SINGULAR);
  
   return;
}
	 

//---------------------------------------------------------------------

template <class Dt,class EACT>
void
Alpha_shape_3<Dt,EACT>::update_alpha_shape_facet_list() const
{
  alpha_shape_facets_list.clear();
  use_facet_cache = true;
  // Writes the faces of the alpha shape `A' for the current 'alpha'-value
  // to the container where 'out' refers to.

  std::back_insert_iterator<std::list< Facet> >
    it = back_inserter(alpha_shape_facets_list);

  get_alpha_shape_facets(it, REGULAR);
  if (get_mode()==GENERAL) get_alpha_shape_facets(it, SINGULAR);
  
  return;
}



//---------------------------------------------------------------------

template <class Dt,class EACT>
typename Alpha_shape_3<Dt,EACT>::Classification_type  
Alpha_shape_3<Dt,EACT>::classify(const Alpha_status& as,
			    const NT& alpha) const
{
 //tetrahedra with circumradius=alpha are considered inside
  if ( !as.is_on_chull() && alpha >= as.alpha_max()) return INTERIOR;
  else if ( alpha >= as.alpha_mid()) return REGULAR;
  else if ( get_mode() == GENERAL && 
	    as.is_Gabriel() &&
	    alpha >= as.alpha_min()) return SINGULAR;
  else return EXTERIOR;
}

template <class Dt,class EACT>
typename Alpha_shape_3<Dt,EACT>::Classification_type  
Alpha_shape_3<Dt,EACT>::classify(const Alpha_status* as,
			    const NT& alpha) const
{
 //tetrahedra with circumradius=alpha are considered inside
  if ( !as->is_on_chull() && alpha >= as->alpha_max()) return INTERIOR;
  else if ( alpha >= as->alpha_mid()) return REGULAR;
  else if ( get_mode() == GENERAL && 
	    as->is_Gabriel() &&
	    alpha >= as->alpha_min()) return SINGULAR;
  else return EXTERIOR;
}

template <class Dt,class EACT>
typename Alpha_shape_3<Dt,EACT>::Classification_type  
Alpha_shape_3<Dt,EACT>::classify(Alpha_status_const_iterator as,
			    const NT& alpha) const
{
  return classify(&(*as), alpha);
}

template <class Dt,class EACT>
typename Alpha_shape_3<Dt,EACT>::Classification_type  
Alpha_shape_3<Dt,EACT>::classify(const Cell_handle& s, 
			    int i,
			    const NT& alpha) const
  // Classifies the face `f' of the underlying Delaunay
  // tetrahedralization with respect to `A'.
{ 
  if (is_infinite(s,i))   return EXTERIOR;
  Alpha_status_iterator as = s->get_facet_status(i);
  return classify(as, alpha);
}
 

template <class Dt,class EACT>
typename Alpha_shape_3<Dt,EACT>::Classification_type  
Alpha_shape_3<Dt,EACT>::classify(const Cell_handle& c, 
			    int i,
			    int j,
			    const NT& alpha) const
  // Classifies the edge `e' of the underlying Delaunay
  // tetrahedralization with respect to `A'.
{ 
  if (is_infinite(c, i, j))     return EXTERIOR;
  if (get_mode() == GENERAL) {
    Alpha_status_iterator asit;
    Vertex_handle_pair
      vhp=make_vertex_handle_pair(c->vertex(i),c->vertex(j));
    asit = edge_alpha_map.find(vhp)->second;
    return classify(asit,alpha);
  }
  
  //no edge map in REGULARIZED mode - classify on the fly
  Alpha_status as;
  compute_edge_status( c, i, j, as);
  return classify(as, alpha);
}

template <class Dt,class EACT>
void
Alpha_shape_3<Dt,EACT>::
compute_edge_status( const Cell_handle& c, 
		     int i, 
		     int j,  
		     Alpha_status& as) const
{
  Facet_circulator fcirc, done;
  Alpha_status_iterator asf;
  NT alpha;
  as.set_is_on_chull(false);
  
  Cell_circulator ccirc, last;
  ccirc = incident_cells(c,i,j);
  last=ccirc;
  while (is_infinite(ccirc) ) ++ccirc; //skip infinite incident cells
  alpha = (*ccirc).get_alpha();
  as.set_alpha_mid(alpha); // initialise as.alpha_mid to alpha value of an incident cell
  as.set_alpha_max(alpha); // same for as.alpha_max 
  while (++ccirc != last) 
  {
    if (!is_infinite(ccirc)) {
      alpha = (*ccirc).get_alpha();
      if (alpha < as.alpha_mid())
        as.set_alpha_mid(alpha);
      if ( ! as.is_on_chull()) {
        if( as.alpha_max() <  alpha)
          as.set_alpha_max( alpha );
      }
    }
  }   
  
  fcirc = incident_facets(c,i,j);
  done = fcirc;  
  do {
    if (!is_infinite(*fcirc)) {
      asf = (*fcirc).first->get_facet_status((*fcirc).second);
      if (get_mode() == GENERAL && asf->is_Gabriel()){
        alpha = asf->alpha_min();
        if (alpha < as.alpha_mid())  as.set_alpha_mid(alpha);
      }
      if (asf->is_on_chull())
        as.set_is_on_chull(true);
    }
  } while (++fcirc != done);  

  // initialize alphamin
  if ( get_mode() == GENERAL){
    if (is_Gabriel(c,i,j)) {
      alpha = squared_radius(c,i,j);
      as.set_is_Gabriel(true);
      as.set_alpha_min(alpha);
    }
    else as.set_is_Gabriel(false);
  }   
}

//---------------------------------------------------------------------

template <class Dt,class EACT>
typename Alpha_shape_3<Dt,EACT>::Classification_type  
Alpha_shape_3<Dt,EACT>::classify(const Vertex_handle& v,
			    const NT& alpha) const
  // Classifies the vertex `v' of the underlying Delaunay
  // tetrahedralization with respect to `A'.
{
  if (is_infinite(v))     return EXTERIOR;
  Alpha_status* as = v->get_alpha_status();
  return classify(as, alpha);
}

//--------------------- NB COMPONENTS ---------------------------------

template <class Dt,class EACT>
typename Alpha_shape_3<Dt,EACT>::size_type
Alpha_shape_3<Dt,EACT>::number_of_solid_components(const NT& alpha) const
    // Determine the number of connected solid components 
    // takes time O(#alpha_shape) amortized if STL_HASH_TABLES
    //            O(#alpha_shape log n) otherwise
{
  typedef typename Marked_cell_set::Data Data;
  Marked_cell_set marked_cell_set(false);
  Finite_cells_iterator cell_it, done = finite_cells_end();
  size_type nb_solid_components = 0;

  // only finite simplices
  for( cell_it = finite_cells_begin(); cell_it != done; ++cell_it)
    {
      Cell_handle pCell = cell_it;
      CGAL_triangulation_assertion(pCell != NULL);
      
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


template <class Dt,class EACT>
void Alpha_shape_3<Dt,EACT>::traverse(Cell_handle pCell,
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
	CGAL_triangulation_assertion(pNeighbor != NULL);
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

template <class Dt,class EACT>
typename Alpha_shape_3<Dt,EACT>::Alpha_iterator 
Alpha_shape_3<Dt,EACT>::find_optimal_alpha(size_type nb_components) const
  // find the minimum alpha that satisfies the properties
  // (1) nb_components solid components <= nb_components
  // (2) all data points on the boundary or in its interior
{
  NT alpha = find_alpha_solid();
  // from this alpha on the alpha_solid satisfies property (2)
  
  Alpha_iterator first = alpha_lower_bound(alpha);
  if (number_of_solid_components(alpha) == nb_components)
    {
      // if ((first+1) < alpha_end()) 
      // return (first+1); 
      // else 
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

 /*      //#ifdef DEBUG */
/*       std::cerr << "first : " << *first  */
/* 		<< " last : "  */
/* 		<< ((first+len != last) ? *(first+len) : *(last-1)) */
/* 		<< " mid : " << *middle  */
/* 		<< " nb comps : " << number_of_solid_components(*middle)  */
/* 		<< std::endl; */
/*       //#endif // DEBUG */

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

 /*  std::cerr << "a la fin " << std::endl */
/* 	    << "first : " << *first  */
/* 	    << " nb comps : " << number_of_solid_components(*first) */
/* 	    << std::endl; */
/*   if ((first+1) < alpha_end())  */
/*     std::cerr << "first+1 " << *(first+1)  */
/* 	      << " nb comps : " << number_of_solid_components(*(first+1)) */
/* 	      << std::endl; */
/*   std::cerr << std::endl; */

  if (number_of_solid_components(*first) <= nb_components ) return first;
  else return first+1;
}  	

//----------------------------------------------------------------------

template <class Dt,class EACT>
typename Alpha_shape_3<Dt,EACT>::NT 
Alpha_shape_3<Dt,EACT>::find_alpha_solid() const
  // compute the minumum alpha such that all data points 
  // are either on the boundary or in the interior
  // not necessarily connected
{
  return _alpha_solid;
}

// TO  DEBUG

template <class Dt,class EACT>
void 
Alpha_shape_3<Dt,EACT>::print_maps() const
{
  typename Alpha_cell_map::const_iterator cit ;
  typename Alpha_facet_map::const_iterator fit ;
  typename Alpha_edge_map::const_iterator eit ;
  typename Alpha_vertex_map::const_iterator vit;

  std::cerr << "size of cell map " << alpha_cell_map.size() 
	    <<   std::endl;
  std::cerr << "size of facet map " << alpha_min_facet_map.size() <<
    std::endl;
  std::cerr << "size of edge map " << alpha_min_edge_map.size() <<
    std::endl;
  std::cerr << "size of vertex map " << alpha_min_vertex_map.size() <<
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
  std::cerr << "alpha_min_vertex_map " << std::endl;
  for(vit = alpha_min_vertex_map.begin();
      vit != alpha_min_vertex_map.end(); ++vit) {
    std::cerr << vit->first << std::endl;
  }
  std::cerr << std::endl;
}


template <class Dt,class EACT>
void 
Alpha_shape_3<Dt,EACT>::print_alphas() const
{
  std::cerr << std::endl;
  std::cerr << " alpha values of facets" << std::endl;
  for(Finite_facets_iterator fit = finite_facets_begin();
      fit != finite_facets_end();
      ++fit) {
    Alpha_status_iterator as = fit->first->get_facet_status(fit->second);
    print_alpha_status(*as);
  }
  std::cerr << std::endl;
  std::cerr << " alpha values of edges " << std::endl;
  if (get_mode() == GENERAL) {
    for(Finite_edges_iterator eit = finite_edges_begin();
	eit != finite_edges_end();
	++eit) {
      Vertex_handle_pair 
	vhp = make_vertex_handle_pair(eit->first->vertex(eit->second),
				      eit->first->vertex(eit->third));
      Alpha_status_iterator as = edge_alpha_map.find(vhp)->second;
      print_alpha_status(*as);
    }
  }
  std::cerr << std::endl;
  std::cerr << " alpha values of vertices " << std::endl;
  for(Finite_vertices_iterator vit = finite_vertices_begin();
      vit != finite_vertices_end();
      ++vit) {
     Alpha_status*  as = vit->get_alpha_status();
     print_alpha_status(*as);
  }

}

template <class Dt,class EACT>
void 
Alpha_shape_3<Dt,EACT>::print_alpha_status(const Alpha_status& as) const
{
  if ( get_mode() == GENERAL &&  as.is_Gabriel())
  std::cerr << as.alpha_min() ;
  else std::cerr <<  "---   " ;
  std::cerr << "\t";
  std::cerr <<  as.alpha_mid()  << "\t";
  if(as.is_on_chull()) std::cerr <<  "---   ";
  else   std::cerr << as.alpha_max();
  std::cerr << std::endl;
}

} //namespace CGAL

#ifdef CGAL_USE_GEOMVIEW
#include <CGAL/IO/alpha_shape_geomview_ostream_3.h>
#endif

#endif //CGAL_ALPHA_SHAPE_3_H
