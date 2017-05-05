// Copyright (c) 1997, 2012  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
// Author(s)     : Tran Kai Frank DA
//                 Andreas Fabri <Andreas.Fabri@geometryfactory.com>

#ifndef CGAL_ALPHA_SHAPE_2_H
#define CGAL_ALPHA_SHAPE_2_H

#include <CGAL/license/Alpha_shapes_2.h>


#include <CGAL/basic.h>

#include <list>
#include <set>
#include <map>
#include <vector>
#include <algorithm>
#include <utility>
#include <iostream>

#include <CGAL/utility.h>
#include <CGAL/Unique_hash_map.h>
#include <CGAL/Triangulation_vertex_base_2.h>
#include <CGAL/Triangulation_face_base_2.h>
#include <CGAL/Alpha_shape_vertex_base_2.h>
#include <CGAL/Alpha_shape_face_base_2.h>
#include <CGAL/internal/Lazy_alpha_nt_2.h>



namespace CGAL {

template < class Dt,class ExactAlphaComparisonTag = Tag_false>
class Alpha_shape_2 : public Dt 
{
  // DEFINITION The class Alpha_shape_2<Dt> represents the family
  // of alpha-shapes of points in a plane for all positive alpha. It
  // maintains the underlying Delaunay triangulation which represents
  // connectivity and order among its faces. Each k-dimensional face of the
  // Delaunay triangulation is associated with an interval that specifies
  // for which values of alpha the face belongs to the alpha-shape (sorted
  // linear arrays resp. multimaps or interval trees). There are links
  // between the intervals and the k-dimensional faces of the Delaunay
  // triangulation (multimaps resp. std::hashtables).
  //

  //------------------------- TYPES ------------------------------------

public:
  typedef Dt Triangulation;
  typedef typename Dt::Geom_traits Gt;
  typedef typename Dt::Triangulation_data_structure Tds;

  typedef typename internal::Alpha_nt_selector_2<
    Gt, ExactAlphaComparisonTag, typename Dt::Weighted_tag>::Type_of_alpha  Type_of_alpha;
  typedef typename internal::Alpha_nt_selector_2<
    Gt, ExactAlphaComparisonTag, typename Dt::Weighted_tag>::Compute_squared_radius_2 Compute_squared_radius_2;
  typedef typename internal::Alpha_nt_selector_2<
    Gt, ExactAlphaComparisonTag, typename Dt::Weighted_tag>::Side_of_bounded_circle_2 Side_of_bounded_circle_2;

  typedef Type_of_alpha               NT;
  typedef Type_of_alpha               FT;

  //check simplices are correctly instantiated
  CGAL_static_assertion( (boost::is_same<NT, typename Dt::Face::NT>::value) );
  CGAL_static_assertion( (boost::is_same<NT, typename Dt::Vertex::NT>::value) );

  typedef typename Dt::Point Point;

  typedef typename Gt::Point_2 Point_2;
  typedef typename Gt::Segment_2 Segment;
  typedef typename Gt::Line_2 Line;

  typedef typename Dt::Face_handle Face_handle;
  typedef typename Dt::Vertex_handle Vertex_handle;
  typedef typename Dt::Edge Edge;

  typedef typename Dt::Face_circulator Face_circulator;
  typedef typename Dt::Edge_circulator Edge_circulator;
  typedef typename Dt::Vertex_circulator Vertex_circulator;

  typedef typename Dt::Finite_faces_iterator Finite_faces_iterator;
  typedef typename Dt::Edge_iterator Edge_iterator;
  typedef typename Dt::Finite_vertices_iterator Finite_vertices_iterator;

  typedef typename Dt::Locate_type Locate_type;
  typedef typename Dt::size_type size_type;

  using Dt::finite_vertices_begin;
  using Dt::finite_vertices_end;
  using Dt::faces_begin;
  using Dt::faces_end;
  using Dt::edges_begin;
  using Dt::edges_end;
  using Dt::number_of_vertices;
  using Dt::cw;
  using Dt::ccw;
  using Dt::VERTEX;
  using Dt::EDGE;
  using Dt::FACE;
  using Dt::OUTSIDE_CONVEX_HULL;
  using Dt::OUTSIDE_AFFINE_HULL;
  using Dt::dimension;
  using Dt::is_infinite;

  // for backward compatibility
  typedef Finite_vertices_iterator Vertex_iterator;
  typedef Finite_faces_iterator Face_iterator;

private:

  typedef std::multimap< Type_of_alpha, Face_handle >  Interval_face_map;
  typedef typename Interval_face_map::value_type    Interval_face;

  typedef typename Tds::Face::Interval_3            Interval3;
  
  typedef std::multimap< Interval3, Edge >          Interval_edge_map;
  typedef typename Interval_edge_map::value_type    Interval_edge;

  typedef std::pair< Type_of_alpha, Type_of_alpha > Interval2;
  typedef std::multimap< Interval2, Vertex_handle > Interval_vertex_map;
  typedef typename Interval_vertex_map::value_type  Interval_vertex;

  typedef Face_handle const const_void;
  typedef std::pair<const_void, int> const_Edge;

  typedef std::vector< Type_of_alpha > Alpha_spectrum;
  
  typedef std::vector< Segment > Vect_seg;

  typedef Unique_hash_map< Face_handle, bool > Marked_face_set;

public:

  typedef typename std::list< Vertex_handle >::iterator 
  Alpha_shape_vertices_iterator;
  typedef typename std::list< Edge >::iterator Alpha_shape_edges_iterator;

  typedef typename Alpha_spectrum::const_iterator Alpha_iterator;
  // An iterator that allow to traverse the sorted sequence of
  // different alpha-values. The iterator is bidirectional and
  // non-mutable. Its value-type is Type_of_alpha

  enum Classification_type {EXTERIOR, SINGULAR, REGULAR, INTERIOR};
  // Distinguishes the different cases for classifying a
  // k-dimensional face of the underlying Delaunay triangulation of
  // the alpha-shape.
  // 
  // `EXTERIOR' if the face does not belong to the alpha-complex.
  // 
  // `SINGULAR' if the face belongs to the boundary of the
  // alpha-shape, but is not incident to any higher-dimensional
  // face of the alpha-complex
  // 
  // `REGULAR' if face belongs to the boundary of the alpha-shape
  // and is incident to a higher-dimensional face of the
  // alpha-complex
  // 
  // `INTERIOR' if the face belongs to the alpha-complex, but does
  // not belong to the boundary of the alpha-shape

  enum Mode {GENERAL, REGULARIZED};
  // In general, an alpha shape is a non-connected, mixed-dimension
  // polygon. Its regularized version is formed by the set of
  // regular edges and their vertices

  //------------------------ private VARIABLES -------------------------

private:

  // only finite edges and faces are inserted into the maps 
  Interval_face_map _interval_face_map;
  Interval_edge_map _interval_edge_map;
  Interval_vertex_map _interval_vertex_map;

  Alpha_spectrum _alpha_spectrum;
 
  Type_of_alpha _alpha;
  Mode _mode;
  // should be constants
  Type_of_alpha Infinity;
  Type_of_alpha UNDEFINED;
  
  mutable std::list< Vertex_handle > Alpha_shape_vertices_list;
  mutable std::list< Edge > Alpha_shape_edges_list;

  mutable bool use_vertex_cache;
  mutable bool use_edge_cache;
public:

  //------------------------- CONSTRUCTORS ------------------------------
 
  // Introduces an empty alpha-shape `A' for a positive
  // alpha-value `alpha'. Precondition: `alpha' >= 0.
  Alpha_shape_2(Type_of_alpha alpha = Type_of_alpha(0), 
		Mode m = GENERAL)
    : _alpha(alpha), _mode(m), Infinity(-1), UNDEFINED(-2),
      use_vertex_cache(false), use_edge_cache(false)
    {}
 
  // Introduces an alpha-shape `A' for a positive alpha-value
  // `alpha' that is initialized with the points in the range
  // from first to last
  template <class InputIterator>
  Alpha_shape_2(const InputIterator& first,  
		const InputIterator& last,  
		const Type_of_alpha& alpha = Type_of_alpha(0),
		Mode m = GENERAL)
    : _alpha(alpha), _mode(m), Infinity(-1), UNDEFINED(-2) ,
      use_vertex_cache(false), use_edge_cache(false)
    {
      Dt::insert(first, last);
      if (dimension() == 2)
	{
	  // Compute the associated _interval_face_map
	  initialize_interval_face_map();

	  // Compute the associated _interval_edge_map
	  initialize_interval_edge_map();
   
	  // Compute the associated _interval_vertex_map
	  initialize_interval_vertex_map();

	  // merge the two maps
	  initialize_alpha_spectrum();
	}
    }

  // Introduces an alpha-shape `A' for a positive alpha-value
  // `alpha' that is initialized with the triangulation 
  Alpha_shape_2(Dt& dt,
		const Type_of_alpha& alpha = Type_of_alpha(0),
		Mode m = GENERAL)
    : _alpha(alpha), _mode(m), Infinity(-1), UNDEFINED(-2) ,
      use_vertex_cache(false), use_edge_cache(false)
    {
      Dt::swap(dt);
      if (dimension() == 2)
	{
	  // Compute the associated _interval_face_map
	  initialize_interval_face_map();

	  // Compute the associated _interval_edge_map
	  initialize_interval_edge_map();
   
	  // Compute the associated _interval_vertex_map
	  initialize_interval_vertex_map();

	  // merge the two maps
	  initialize_alpha_spectrum();
	}
    }
 
public:

  //----------- OUTPUT POINTS CONNECTED BY PAIRS ----------------------

  std::list<Point_2> Output();
 
  std::ostream& op_ostream(std::ostream& os) const;

 
  //----------------------- OPERATIONS ---------------------------------

  // Introduces an alpha-shape `A' for a positive alpha-value
  // `alpha' that is initialized with the points in the range
  // from first to last

  template < class InputIterator >  
  std::ptrdiff_t make_alpha_shape(const InputIterator& first,  
		       const InputIterator& last) 
    {

      clear();

      size_type n =  Dt::insert(first, last);
 
      if (dimension() == 2)
	{   
	  // Compute the associated _interval_face_map
	  initialize_interval_face_map();

	  // Compute the associated _interval_edge_map
	  initialize_interval_edge_map();

	  // Compute the associated _interval_vertex_map
	  initialize_interval_vertex_map();
  
	  // merge the two maps
	  initialize_alpha_spectrum();
	}
      return n;
    }

private :

  //--------------------- INITIALIZATION OF PRIVATE MEMBERS -----------

  void initialize_interval_face_map();

  void initialize_interval_edge_map();

  void initialize_interval_vertex_map();

  void initialize_alpha_spectrum();

  //---------------------------------------------------------------------

public:

  void clear() 
    {
      // clears the structure
      Dt::clear();

      _interval_face_map.clear();
      _interval_edge_map.clear();
      _interval_vertex_map.clear();
    
      _alpha_spectrum.clear();

      Alpha_shape_vertices_list.clear();
      Alpha_shape_edges_list.clear();
    
      set_alpha(Type_of_alpha(0)); 
      set_mode(GENERAL);

    }

  //----------------------- PRIVATE MEMBERS --------------------------

private:
   
  struct Less 
  {
    bool operator()(const Interval_edge& ie, 
		    const Type_of_alpha& alpha) 
      {
	return ie.first.first < alpha; 
      }

    bool operator()( const Type_of_alpha& alpha, 
		     const Interval_edge& ie) 
      {
	return alpha < ie.first.first; 
      }

    // Needed for STL implementations of upper_bound which in debug mode 
    // check sortedness of range
    bool operator()(const Interval_edge& ie, 
		    const Interval_edge& ie2) const
      {
	return ie < ie2;
      }
  };

  
  //----------------------- ACCESS TO PRIVATE MEMBERS -----------------

private:
   
  
  Type_of_alpha find_interval(const Face_handle& f) const 
    {
      return f->get_alpha();
      // return the value Alpha f the face f
    }
 
  Interval3 find_interval(const_Edge e)  const 
    { // corriger le parametrage
      return (e.first)->get_ranges(e.second);
      // return the Interval3 for the edge n
    }


  //---------------------------------------------------------------------

public:

  Type_of_alpha set_alpha(const Type_of_alpha& alpha) 
    {
      // Sets the alpha-value to `alpha'. Precondition: `alpha' >= 0.
      // Returns the previous alpha
      Type_of_alpha previous_alpha = _alpha;
      _alpha = alpha;
      use_vertex_cache = false;
      use_edge_cache = false;
      return previous_alpha;
    }

  const Type_of_alpha&  get_alpha() const 
    {
      // Returns the current alpha-value.
      return _alpha;
    }
  
  const Type_of_alpha&  get_nth_alpha(size_type n) const 
    {
      // Returns the n-th alpha-value.
      // n < size()
      if (! _alpha_spectrum.empty())
	return _alpha_spectrum[n];
      else
	return UNDEFINED;
    }
  
  size_type number_of_alphas() const 
    {
      // Returns the number of not necessary different alpha-values
      return _alpha_spectrum.size();
    }

  //---------------------------------------------------------------------

private:

  // the dynamic version is not yet implemented
  // desactivate the triangulation member functions
  Vertex_handle insert(const Point& p);
  // Inserts point `p' in the alpha shape and returns the
  // corresponding vertex of the underlying Delaunay triangulation.
  // If point `p' coincides with an already existing vertex, this
  // vertex is returned and the alpha shape remains unchanged.
  // Otherwise, the vertex is inserted in the underlying Delaunay
  // triangulation and the associated intervals are updated.

  void remove(Vertex_handle v);
  // Removes the vertex from the underlying Delaunay triangulation.
  // The created hole is retriangulated and the associated intervals
  // are updated.

  //---------------------------------------------------------------------

public:

  Mode set_mode(Mode mode = GENERAL ) 
    {
      // Sets `A' to its general or regularized version. Returns the
      // previous mode.

      Mode previous_mode = _mode;
      _mode = mode;
      return previous_mode;
    }

  Mode get_mode() const 
    {
      // Returns whether `A' is general or regularized.
      return _mode;
    }

  //---------------------------------------------------------------------

private:


void
update_alpha_shape_vertex_list()const;

  //---------------------------------------------------------------------

  void
  update_alpha_shape_edges_list() const;
  

  //---------------------------------------------------------------------

public:

  Alpha_shape_vertices_iterator alpha_shape_vertices_begin() const
    { 
    if(!use_vertex_cache){
      update_alpha_shape_vertex_list();
    }
      return Alpha_shape_vertices_list.begin();
    }



  // for backward compatibility
  Alpha_shape_vertices_iterator Alpha_shape_vertices_begin()
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
      return Alpha_shape_vertices_list.end();
    }
  //---------------------------------------------------------------------

  Alpha_shape_edges_iterator alpha_shape_edges_begin() const
    {
      if(!use_edge_cache){
      update_alpha_shape_edges_list();
    }
      return Alpha_shape_edges_list.begin();
    }

  Alpha_shape_edges_iterator Alpha_shape_edges_begin() const
  {
    return alpha_shape_edges_begin();
  }
  //---------------------------------------------------------------------

  Alpha_shape_edges_iterator alpha_shape_edges_end() const
    {
      return Alpha_shape_edges_list.end();
    }

  Alpha_shape_edges_iterator Alpha_shape_edges_end() const
    {
      return Alpha_shape_edges_list.end();
    }
public: 
  
  // Traversal of the alpha-Values
  // 
  // The alpha shape class defines an iterator that allows to
  // visit the sorted sequence of alpha-values. This iterator is
  // non-mutable and bidirectional. Its value type is Type_of_alpha.

  Alpha_iterator alpha_begin() const 
    {
      // Returns an iterator that allows to traverse the sorted sequence
      // of alpha-values of `A'.
      return _alpha_spectrum.begin(); 
    }

  Alpha_iterator alpha_end() const 
    {
      // Returns the corresponding past-the-end iterator.
      return _alpha_spectrum.end(); 
    }
  
  Alpha_iterator alpha_find(const Type_of_alpha& alpha) const 
    {
      // Returns an iterator pointing to an element with alpha-value
      // `alpha', or the corresponding past-the-end iterator if such an
      // element is not found.
      return std::find(_alpha_spectrum.begin(),
		       _alpha_spectrum.end(),
		       alpha);
    }

  Alpha_iterator alpha_lower_bound(const Type_of_alpha& alpha) const 
    {
      // Returns an iterator pointing to the first element with
      // alpha-value not less than `alpha'.
      return std::lower_bound(_alpha_spectrum.begin(),
			      _alpha_spectrum.end(),
			      alpha);
    }

  Alpha_iterator alpha_upper_bound(const Type_of_alpha& alpha) const 
    {
      // Returns an iterator pointing to the first element with
      // alpha-value greater than `alpha'.
      return std::upper_bound(_alpha_spectrum.begin(),
			      _alpha_spectrum.end(),
			      alpha);
    }

  //--------------------- PREDICATES -----------------------------------

  // the classification predicates take 
  //      amortized const time if STL_STD::HASH_TABLES
  //      O(log #alpha_shape ) otherwise

  Classification_type  classify(const Point& p ) const 
    {
      return classify( p, get_alpha());
    }

  
  Classification_type  classify(const Point& p,   
				const Type_of_alpha& alpha) const 
    {
      // Classifies a point `p' with respect to `A'.
      Locate_type type;
      int i;
      Face_handle pFace = this->locate(p, type, i);
      switch (type) 
	{
	case VERTEX            : return classify(pFace->vertex(i), alpha);
	case EDGE              : return classify(pFace, i, alpha);
	case FACE              : return classify(pFace, alpha);
	case OUTSIDE_CONVEX_HULL : 
	case OUTSIDE_AFFINE_HULL : return EXTERIOR;
	default                : return EXTERIOR;
	}
    }

  //---------------------------------------------------------------------

  Classification_type  classify(const Face_handle& f) const 
    {
      // Classifies the face `f' of the underlying Delaunay
      // triangulation with respect to `A'.
      return classify(f, get_alpha());
    }
  
  Classification_type  classify(const Face_handle& f, 
				const Type_of_alpha& alpha) const 
    {
      // Classifies the face `f' of the underlying Delaunay
      // triangulation with respect to `A'.
      // we consider close circles :
      // f->radius < alpha <=> f exterior
      // problem the operator [] is non-const

      if (is_infinite(f)) return EXTERIOR;
    
      // the version that computes the squared radius seems to be 
      // much faster
    
      return (find_interval(f) <= alpha) ? 
	INTERIOR :
	EXTERIOR;

    }

  //---------------------------------------------------------------------
  
  Classification_type  classify(const Edge& edge) const
    {  
      return classify(edge.first, edge.second, get_alpha());
    }

  Classification_type  classify(const Face_handle& f, 
				int i) const 
    {  
      return classify(f, i, get_alpha());
    }

  Classification_type  classify(const Edge& edge,
				const Type_of_alpha& alpha) const 
    {  
      return classify(edge.first, edge.second, alpha);
    }

  Classification_type  classify(const Face_handle& f, 
				int i,
				const Type_of_alpha& alpha) const;

  
  //---------------------------------------------------------------------

  Classification_type  classify(const Vertex_handle& v) const 
    {
      return classify(v, get_alpha());
    }

  Classification_type  classify(const Vertex_handle& v,
				const Type_of_alpha& alpha) const;

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

  size_type 
  number_solid_components(const Type_of_alpha& /* alpha */) const
    {
      return number_of_solid_components(get_alpha());
    }

  size_type
  number_of_solid_components(const Type_of_alpha& alpha) const;

private:

  void traverse(const Face_handle& pFace,
		Marked_face_set& marked_face_set, 
		const Type_of_alpha alpha) const;

//  class Line_face_circulator;

  //----------------------------------------------------------------------

public:

  Alpha_iterator find_optimal_alpha(size_type nb_components);  	

private:

  Type_of_alpha find_alpha_solid() const;

  //---------------------- PREDICATES ------------------------------------

private:

  bool is_attached(const Face_handle& f, int i) const
  {
    Bounded_side b = Side_of_bounded_circle_2()(*this)(f->vertex(cw(i))->point(),
                                                       f->vertex(ccw(i))->point(),
                                                       f->vertex(i)->point());

    return (b == ON_BOUNDED_SIDE) ? true : false;
  }

  //-------------------- GEOMETRIC PRIMITIVES ----------------------------

  Type_of_alpha squared_radius(const Face_handle& f) const
  {
    return Compute_squared_radius_2()(*this)(f->vertex(0)->point(),
                                             f->vertex(1)->point(),
                                             f->vertex(2)->point());
  }

  Type_of_alpha squared_radius(const Face_handle& f, int i) const
  {
    return Compute_squared_radius_2()(*this)(f->vertex(ccw(i))->point(),
                                             f->vertex(cw(i))->point());
  }

  //---------------------------------------------------------------------

private:
  // prevent default copy constructor and default assigment
  
  Alpha_shape_2(const Alpha_shape_2& A);

  Alpha_shape_2& operator=(const Alpha_shape_2& A);


public:
  // to debug
  void print_edge_map();

};


//---------------------------------------------------------------------
//----------------------- MEMBER FUNCTIONS -----------------------------
//---------------------------------------------------------------------


template < class Dt, class EACT >
void 
Alpha_shape_2<Dt,EACT>::initialize_interval_face_map() 
{
  Type_of_alpha alpha_f;

  // only finite faces
  for(Finite_faces_iterator face_it = faces_begin(); face_it != faces_end(); ++face_it)
    {
      alpha_f = squared_radius(face_it);
      _interval_face_map.insert(Interval_face(alpha_f, face_it));

      // cross references
      face_it->set_alpha(alpha_f);
    } 
}

//-------------------------------------------------------------------------

template < class Dt, class EACT >
void 
Alpha_shape_2<Dt,EACT>::initialize_interval_edge_map() 
{  
  Edge_iterator edge_it;
  Edge edge;

  // only finite faces
  for( edge_it = edges_begin(); 
       edge_it != edges_end(); 
       ++edge_it) 
    { 
      Interval3 interval;
      edge = (*edge_it);

      Face_handle pFace = edge.first;
      int i = edge.second;
      Face_handle pNeighbor = pFace->neighbor(i);
      int Neigh_i = pNeighbor->index(pFace);
      
      // not on the convex hull
      
      if(!is_infinite(pFace) && !is_infinite(pNeighbor)) 
	{ 
	  Type_of_alpha squared_radius_Face = 
	    find_interval(pFace);
	  Type_of_alpha squared_radius_Neighbor = 
	    find_interval(pNeighbor);
	  if (squared_radius_Neighbor < squared_radius_Face) 
	    {	    
	      edge =  Edge(pNeighbor, Neigh_i);
	      Type_of_alpha coord_tmp = squared_radius_Face;
	      squared_radius_Face = squared_radius_Neighbor;
	      squared_radius_Neighbor = coord_tmp;
	    }
	  interval = (is_attached(pFace, i) || 
		      is_attached(pNeighbor, Neigh_i)) ?
	    make_triple(UNDEFINED,
			squared_radius_Face, 
			squared_radius_Neighbor):
	    make_triple(squared_radius(pFace, i),
			squared_radius_Face, 
			squared_radius_Neighbor);
	}
      else 
	{ // on the convex hull

	  if(is_infinite(pFace)) 
	    {
	      if (!is_infinite(pNeighbor)) 
		{
		  interval =  (is_attached(pNeighbor,
					   Neigh_i)) ?
		    make_triple(UNDEFINED,
				find_interval(pNeighbor), 
				Infinity):
		    make_triple(squared_radius(pNeighbor, 
					       Neigh_i), 
				find_interval(pNeighbor), 
				Infinity); 
		  edge =  Edge(pNeighbor, Neigh_i);
		}
	      else 
		{
		  // both faces are infinite by definition unattached
		  // the edge is finite by construction
		  CGAL_triangulation_precondition((is_infinite(pNeighbor) 
						   && is_infinite(pFace)));
		  interval = make_triple(squared_radius(pFace, i), 
					 Infinity,
					 Infinity);
		}
	    }
	  else 
	    { // is_infinite(pNeighbor)
	   
	      CGAL_triangulation_precondition((is_infinite(pNeighbor) 
					       && !is_infinite(pFace)));
	      if (is_attached(pFace, i))
		interval = make_triple(UNDEFINED,
				       find_interval(pFace),
				       Infinity);
	      else
		interval = make_triple(squared_radius(pFace, i),
				       find_interval(pFace),
				       Infinity);

	    }
	}
      
      _interval_edge_map.insert(Interval_edge(interval, edge));
      
      // cross-links
      (edge.first)->set_ranges(edge.second,interval);
      // MY : to fix a bug I store the interval of the edge in both faces 
      Face_handle neighbor = (edge.first)->neighbor(edge.second);
      int ni = neighbor->index(edge.first);
      neighbor->set_ranges( ni, interval);

    }

  // Remark:
  // The interval_edge_map will be sorted as follows
  // first the attached edges on the convex hull
  // second                   not on the convex hull
  // third the un-attached edges on the convex hull
  // finally                     not on the convex hull
  // 
  // if we are in regularized mode we should sort differently
  // by the second third first Key
}

//-------------------------------------------------------------------------

template < class Dt, class EACT >
void 
Alpha_shape_2<Dt,EACT>::initialize_interval_vertex_map() 
{
  Type_of_alpha alpha_mid_v;
  Type_of_alpha alpha_max_v;
  Type_of_alpha alpha_f;

  Finite_vertices_iterator vertex_it;

  for( vertex_it = finite_vertices_begin(); 
       vertex_it != finite_vertices_end(); 
       ++vertex_it) 
    {
      Vertex_handle v = vertex_it;
      Face_handle f;

      alpha_max_v = Type_of_alpha(0);    
      alpha_mid_v = (!_interval_face_map.empty() ?
		     (--_interval_face_map.end())->first :
		     Type_of_alpha(0));

      //----------------- examine incident edges --------------------------
      
// 	// if we used Edelsbrunner and Muecke's definition
// 	// singular means not incident to any higher-dimensional face
// 	// regular means incident to a higher-dimensional face
// 	Edge_circulator edge_circ = this->incident_edges(v),
// 	edge_done(edge_circ);

// 	do 
      
// 	{ 
// 	f = (*edge_circ).first;
// 	int i = (*edge_circ).second;

// 	if (is_infinite(f, i))
      
// 	{
// 	alpha_max_v = Infinity;
// 	}
// 	else
      
// 	{
// 	Interval3 interval3 = find_interval(const_Edge(f, i));
				
	
// 	alpha_mid_v = (interval3.first != UNDEFINED) ?
// 	(CGAL::min)(alpha_mid_v, interval3.first): 
// 	(CGAL::min)(alpha_mid_v, interval3.second); 
			
// 	if (alpha_max_v != Infinity)
      
// 	{
// 	alpha_max_v = (interval3.third != Infinity) ?
// 	(CGAL::max)(alpha_max_v, interval3.third):
// 	Infinity;
// 	}
// 	}
// 	}
// 	while(++edge_circ != edge_done);
     

      //-------------- examine incident faces --------------------------
      
      // we use a different definition than Edelsbrunner and Muecke
      // singular means not incident to any 2-dimensional face
      // regular means incident to a 2-dimensional face
     
      Face_circulator face_circ = this->incident_faces(v),
	done = face_circ;
   
      if (!face_circ.is_empty()) 
	{
	  do 
	    {
	      f = face_circ;
	      if (is_infinite(f)) 
		{
		  alpha_max_v = Infinity;
		  // continue;
		}
	      else 
		{
		  alpha_f = find_interval(f);
		  // if we define singular as not incident to a 2-dimensional
		  // face
		  alpha_mid_v = (CGAL::min)(alpha_mid_v, alpha_f);
		    
		  if (alpha_max_v != Infinity)
		    alpha_max_v = (CGAL::max)(alpha_max_v, alpha_f);
			    
		}
	    }
	  while(++face_circ != done);
	}
	
 
      Interval2 interval = std::make_pair(alpha_mid_v, alpha_max_v);
      _interval_vertex_map.insert(Interval_vertex(interval, vertex_it));

      // cross references
      vertex_it->set_range(interval);
    }
}

//-------------------------------------------------------------------------

template < class Dt, class EACT >
void 
Alpha_shape_2<Dt,EACT>::initialize_alpha_spectrum() 
{

  // skip the attached edges 
  // <=> _interval_edge_map.first.first == UNDEFINED
  typename Interval_edge_map::iterator 
    edge_it = std::upper_bound(_interval_edge_map.begin(),
			       _interval_edge_map.end(),
			       UNDEFINED,
			       Less());

  // merge the maps which is sorted and contains the alpha-values
  // of the unattached edges and the triangles.
  // eliminate duplicate values due to for example attached edges
  // merge and copy from STL since assignment should be function object
	
  typename Interval_face_map::iterator
    face_it = _interval_face_map.begin();

  _alpha_spectrum.reserve(_interval_face_map.size() +
			  _interval_edge_map.size()/ 2 );
  // should be only the number of unattached edges
  // size_type nb_unattached_edges;
  // distance(edge_it, _interval_edge_map.end(), nb_unattached_edges);
  // however the distance function is expensive

  while (edge_it != _interval_edge_map.end() ||
	 face_it != _interval_face_map.end()) 
    {
      if (face_it != _interval_face_map.end() &&
	  (edge_it == _interval_edge_map.end() ||
	   ((*face_it).first < (*edge_it).first.first))) 
	{
	  if (((_alpha_spectrum.empty() || 
		_alpha_spectrum.back() < (*face_it).first)) && 
	      ((*face_it).first > Type_of_alpha(0)))
	    _alpha_spectrum.push_back((*face_it).first);
	  face_it++;
	}
      else
	{
	  if (((_alpha_spectrum.empty() || 
		_alpha_spectrum.back() < (*edge_it).first.first)) &&
	      (((*edge_it).first.first) > Type_of_alpha(0)))
	    _alpha_spectrum.push_back((*edge_it).first.first);
	  edge_it++;
	}
    }
    
  while (edge_it != _interval_edge_map.end()) 
    {
      if (((_alpha_spectrum.empty() || 
	    _alpha_spectrum.back() < (*edge_it).first.first))&&
	  (((*edge_it).first.first) > Type_of_alpha(0)))
	_alpha_spectrum.push_back((*edge_it).first.first);
      edge_it++;
    }

  while (face_it != _interval_face_map.end()) 
    { 
      if (((_alpha_spectrum.empty() || 
	    _alpha_spectrum.back() < (*face_it).first))&&
	  ((*face_it).first > Type_of_alpha(0)))
	_alpha_spectrum.push_back((*face_it).first);
      face_it++;
    }

}


//-------------------------------------------------------------------------



template < class Dt, class EACT >
void
Alpha_shape_2<Dt,EACT>::update_alpha_shape_vertex_list()const {
	//typedef typename Alpha_shape_2<Dt,EACT>::Interval_vertex_map 
	//                                              Interval_vertex_map;
	typename Interval_vertex_map::const_iterator vertex_alpha_it;

	//const typename Alpha_shape_2<Dt,EACT>::Interval2* pInterval2;
	const Interval2* pInterval2;
	Vertex_handle v;
	Alpha_shape_vertices_list.clear();
	// write the regular vertices

	for (vertex_alpha_it = _interval_vertex_map.begin(); 
		vertex_alpha_it != _interval_vertex_map.end() &&
		(*vertex_alpha_it).first.first <= get_alpha();
		++vertex_alpha_it) 
		{
		pInterval2 = &(*vertex_alpha_it).first;

		if((pInterval2->second > get_alpha()
		|| pInterval2->second == Infinity)) 
		{
		// alpha must be larger than the min boundary
		// and alpha is smaller than the upper boundary
		// which might be infinity 
		// write the vertex
		v = (*vertex_alpha_it).second;
		CGAL_triangulation_assertion((classify(v) == REGULAR));
		Alpha_shape_vertices_list.push_back(v);
		}
		}
	 
	if (get_mode() == Alpha_shape_2<Dt,EACT>::GENERAL) 
		{
		// write the singular vertices
		for (; 
		vertex_alpha_it != _interval_vertex_map.end();
		++vertex_alpha_it) 
		{
		v = (*vertex_alpha_it).second;
		CGAL_triangulation_assertion((classify(v) == SINGULAR));

		Alpha_shape_vertices_list.push_back(v);
		}
		}
	use_vertex_cache = true;
  }

//-------------------------------------------------------------------------

template < class Dt, class EACT >
void
Alpha_shape_2<Dt,EACT>::update_alpha_shape_edges_list() const 
{

  // Writes the edges of the alpha shape `A' for the current $\alpha$-value
  // to the container where 'out' refers to. Returns an output iterator 
  // which is the end of the constructed range.
  //typedef  typename Alpha_shape_2<Dt,EACT>::Interval_edge_map Interval_edge_map;
  typename Interval_edge_map::const_iterator edge_alpha_it;

  //const typename Alpha_shape_2<Dt,EACT>::Interval3* pInterval;
  const Interval3* pInterval;
  Alpha_shape_edges_list.clear();
  if (get_mode() == REGULARIZED) 
    {
      // it is much faster looking at the sorted intervals 
      // than looking at all sorted faces
      // alpha must be larger than the mid boundary
      // and alpha is smaller than the upper boundary
      for (edge_alpha_it = _interval_edge_map.begin(); 
	   edge_alpha_it != _interval_edge_map.end() &&
	     (*edge_alpha_it).first.first <= get_alpha();
	   ++edge_alpha_it) 
	{
	  pInterval = &(*edge_alpha_it).first;
      
	  CGAL_triangulation_assertion(pInterval->second != Infinity);
	  // since this happens only for convex hull of dimension 2
	  // thus singular

	  if(pInterval->second <= get_alpha() &&
	     (pInterval->third > get_alpha()
	      || pInterval->third == Infinity)) 
	    {
	      // alpha must be larger than the mid boundary
	      // and alpha is smaller than the upper boundary
	      // which might be infinity 
	      // visualize the boundary
 CGAL_triangulation_assertion((classify((*edge_alpha_it).second.first,
					(*edge_alpha_it).second.second)
			       == REGULAR));
	      Alpha_shape_edges_list.push_back(Edge((*edge_alpha_it).second.first,
						    (*edge_alpha_it).second.second));
	    }
	}
    }
  else 
    {  // get_mode() == GENERAL -------------------------------------------
      // draw the edges
      for (edge_alpha_it = _interval_edge_map.begin(); 
	   edge_alpha_it != _interval_edge_map.end() &&
	     (*edge_alpha_it).first.first <= get_alpha();
	   ++edge_alpha_it) 
	{
	  pInterval = &(*edge_alpha_it).first;
      
	  if (pInterval->first == UNDEFINED) 
	    {
	      CGAL_triangulation_assertion(pInterval->second != Infinity);
	      // since this happens only for convex hull of dimension 2
	      // thus singular
	
	      if(pInterval->second <= get_alpha() &&
		 (pInterval->third > get_alpha()
		  || pInterval->third == Infinity)) 
		{
		  // alpha must be larger than the mid boundary
		  // and alpha is smaller than the upper boundary
		  // which might be infinity 
		  // visualize the boundary
 CGAL_triangulation_assertion((classify((*edge_alpha_it).second.first,
					(*edge_alpha_it).second.second) 
			       == REGULAR));
		  Alpha_shape_edges_list.push_back(Edge((*edge_alpha_it).second.first,
							(*edge_alpha_it).second.second));
		}
	    }
	  else 
	    {
	
	      if(pInterval->third > get_alpha()
		 || pInterval->third == Infinity) 
		{
		  // if alpha is smaller than the upper boundary
		  // which might be infinity 
		  // visualize the boundary
 CGAL_triangulation_assertion(((classify((*edge_alpha_it).second.first,
					 (*edge_alpha_it).second.second) 
				== REGULAR)
			       || (classify((*edge_alpha_it).second.first,
					    (*edge_alpha_it).second.second) 
				   == SINGULAR)));
		  Alpha_shape_edges_list.push_back(Edge((*edge_alpha_it).second.first,
							(*edge_alpha_it).second.second));
		}
	    }
      
	}
    
    }
  use_edge_cache = true;
}

//-------------------------------------------------------------------------

template < class Dt, class EACT >
typename Alpha_shape_2<Dt,EACT>::Classification_type  
Alpha_shape_2<Dt,EACT>::classify(const Face_handle& f, int i, 
			    const Type_of_alpha& alpha) const
{
  // Classifies the edge `e' of the underlying Delaunay
  // triangulation with respect to `A'.

  // the version that uses a simplified version without crosslinks
  // is much slower
  if (is_infinite(f, i))
    return EXTERIOR;

  // we store only finite edges in _edge_interval_map
  Interval3 interval = find_interval(const_Edge(f, i));
  //  (*(_edge_interval_map.find(const_Edge(f, i)))).second;
 
  if (alpha < interval.second) 
    {
      if (get_mode() == REGULARIZED ||
	  interval.first == UNDEFINED ||
	  alpha < interval.first)
	return EXTERIOR;
      else // alpha >= interval.first
	return SINGULAR;
    }
  else 
    {   // alpha >= interval.second
      if (interval.third == Infinity ||
	  alpha < interval.third)
	return REGULAR;
      else // alpha >= interval.third
	return INTERIOR;
    }
   
}

//-------------------------------------------------------------------------

template < class Dt, class EACT >
typename Alpha_shape_2<Dt,EACT>::Classification_type  
Alpha_shape_2<Dt,EACT>::classify(const Vertex_handle& v,
			    const Type_of_alpha& alpha) const 
{
  // Classifies the vertex `v' of the underlying Delaunay
  // triangulation with respect to `A'.
  Interval2 interval = v->get_range();
 
  if (alpha < interval.first) 
    {
      if (get_mode() == REGULARIZED) 
	return EXTERIOR;
      else // general => vertices are never exterior
	return SINGULAR;
    }
  else 
    {   // alpha >= interval.first
      if (interval.second == Infinity ||
	  alpha < interval.second)
	return REGULAR;
      else // alpha >= interval.second
	return INTERIOR;
    }
}

//-------------------------------------------------------------------------

template < class Dt, class EACT >
typename Alpha_shape_2<Dt,EACT>::size_type
Alpha_shape_2<Dt,EACT>::number_of_solid_components(const Type_of_alpha& alpha) const
{
  // Determine the number of connected solid components 
  typedef typename Marked_face_set::Data Data;
  Marked_face_set marked_face_set(false);
  Finite_faces_iterator face_it;
  size_type nb_solid_components = 0;
  
  if (number_of_vertices()==0)
    return 0;
  
  // only finite faces
  for( face_it = faces_begin(); 
       face_it != faces_end(); 
       ++face_it) 
    {
      Face_handle pFace = face_it;
      CGAL_triangulation_postcondition( pFace != NULL);
      
      if (classify(pFace, alpha) == INTERIOR){
	Data& data = marked_face_set[pFace];
	if(data == false)
	  {
	    // we traverse only interior faces
	    traverse(pFace, marked_face_set, alpha);
	    nb_solid_components++;  
	  }
      }
    }
  return nb_solid_components;
}

//-------------------------------------------------------------------------

template < class Dt, class EACT >
void
Alpha_shape_2<Dt,EACT>::traverse(const Face_handle& pFace,
			    Marked_face_set& marked_face_set, 
			    const Type_of_alpha alpha) const 
{
  typedef typename Marked_face_set::Data Data;
  std::list<Face_handle> faces;
  faces.push_back(pFace);
  Face_handle pNeighbor, fh;

  while(! faces.empty()){
    fh = faces.front();
    faces.pop_front();
    for (int i=0; i<3; i++)
      {
	pNeighbor = fh->neighbor(i);
	 CGAL_triangulation_assertion(pNeighbor != NULL);
	if (classify(pNeighbor, alpha) == INTERIOR){
	  Data& data = marked_face_set[pNeighbor];
	  if(data == false){
	    data = true;
	    faces.push_back(pNeighbor);
	  }
	}
      }
  }
}

//-------------------------------------------------------------------------

template < class Dt, class EACT >
typename Alpha_shape_2<Dt,EACT>::Alpha_iterator
Alpha_shape_2<Dt,EACT>::find_optimal_alpha(size_type nb_components) 
{
  // find the minimum alpha that satisfies the properties
  // (1) nb_components solid components
  // (2) all data points on the boundary or in its interior
  Type_of_alpha alpha = find_alpha_solid();
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

#ifdef CGAL_DEBUG_ALPHA_SHAPE_2
      std::cout << "first : " << *first << " last : " << *(first+len)
		<< " mid : " << *middle 
		<< " nb comps : " << number_of_solid_components(*middle) 
		<< std::endl;
#endif // CGAL_DEBUG_ALPHA_SHAPE_2

      if (number_of_solid_components(*middle) > nb_components) 
	{
	  first = middle + 1;
	  len = len - half - 1; 
	} 
      else 
	{ // number_of_solid_components(*middle) <= nb_components
	  len = half;
	}
    }
  if ((first+1) < alpha_end())
    return (first+1);
  else
    return first;
}

//-------------------------------------------------------------------------

template < class Dt, class EACT >
typename Alpha_shape_2<Dt,EACT>::Type_of_alpha 
Alpha_shape_2<Dt,EACT>::find_alpha_solid() const 
{
  // compute the minumum alpha such that all data points 
  // are either on the boundary or in the interior
  // not necessarily connected
  // starting point for searching 
  // takes O(#alpha_shape) time
  Type_of_alpha alpha_solid = 0;
  
  if (number_of_vertices()<3) return alpha_solid;

  Finite_vertices_iterator vertex_it;
  // only finite vertices
  for( vertex_it = finite_vertices_begin(); 
       vertex_it != finite_vertices_end();
       ++vertex_it) 
    {
      Type_of_alpha alpha_min_v = (--_interval_face_map.end())->first;

      Face_circulator face_circ = this->incident_faces(vertex_it);
      Face_circulator  done = face_circ;
      do 
	{
	  Face_handle f = face_circ;
	  if (! is_infinite(f))
	    alpha_min_v = (CGAL::min)(find_interval(f),
                                      alpha_min_v);
	}
      while (++face_circ != done);
      alpha_solid = (CGAL::max)(alpha_min_v, alpha_solid);

    }
  return alpha_solid;
}

//-------------------------------------------------------------------------

template < class Dt, class EACT >
std::ostream& 
Alpha_shape_2<Dt,EACT>::op_ostream(std::ostream& os) const
{
  typedef typename Alpha_shape_2<Dt,EACT>::Interval_vertex_map Interval_vertex_map ;
  typedef typename Alpha_shape_2<Dt,EACT>::Interval_edge_map Interval_edge_map;

  typename Interval_vertex_map::const_iterator vertex_alpha_it;
  const typename Alpha_shape_2<Dt,EACT>::Interval2* pInterval2;
  typename Interval_edge_map::const_iterator edge_alpha_it;
  const typename Alpha_shape_2<Dt,EACT>::Interval3* pInterval;

  Unique_hash_map< Vertex_handle, size_type > V;
  size_type number_of_vertices = 0;

  if (get_mode() == Alpha_shape_2<Dt,EACT>::REGULARIZED)
  {
    typename Alpha_shape_2<Dt,EACT>::Vertex_handle v;
    for (vertex_alpha_it = _interval_vertex_map.begin();
         vertex_alpha_it != _interval_vertex_map.end() &&
         (*vertex_alpha_it).first.first <= get_alpha();
         ++vertex_alpha_it)
    {
      pInterval2 = &(*vertex_alpha_it).first;

#ifdef CGAL_DEBUG_ALPHA_SHAPE_2
      typename Alpha_shape_2<Dt,EACT>::Type_of_alpha alpha =
          get_alpha();
      typename Alpha_shape_2<Dt,EACT>::Type_of_alpha alpha_mid =
          pInterval2->first;
      typename Alpha_shape_2<Dt,EACT>::Type_of_alpha alpha_max =
          pInterval2->second;
#endif // CGAL_DEBUG_ALPHA_SHAPE_2

      if((pInterval2->second > get_alpha()
          || pInterval2->second == Infinity))
      {
        // alpha must be larger than the min boundary
        // and alpha is smaller than the upper boundary
        // which might be infinity
        // write the vertex

        v = (*vertex_alpha_it).second;
        CGAL_triangulation_assertion((classify(v) ==
                                      Alpha_shape_2<Dt,EACT>::REGULAR));
        // if we used Edelsbrunner and Muecke's definition
        // regular means incident to a higher-dimensional face
        // we would write too many vertices

        V[v] = number_of_vertices++;
        os << v->point() << std::endl;
      }
    }
    // the vertices are oriented counterclockwise

    typename Alpha_shape_2<Dt,EACT>::Face_handle f;
    int i;

    for (edge_alpha_it = _interval_edge_map.begin();
         edge_alpha_it != _interval_edge_map.end() &&
         (*edge_alpha_it).first.first <= get_alpha();
         ++edge_alpha_it)
    {
      pInterval = &(*edge_alpha_it).first;

      CGAL_triangulation_assertion(pInterval->second != Infinity);
      // since this happens only for convex hull of dimension 1
      // thus singular

      if(pInterval->second <= get_alpha() &&
         (pInterval->third > get_alpha()
          || pInterval->third == Infinity))
      {
        // alpha must be larger than the mid boundary
        // and alpha is smaller than the upper boundary
        // which might be infinity
        // visualize the boundary


        f = (*edge_alpha_it).second.first;
        i = (*edge_alpha_it).second.second;

        // assure that all vertices are in ccw order
        if (classify(f) == Alpha_shape_2<Dt,EACT>::EXTERIOR)
        {
          // take the reverse face
          typename Alpha_shape_2<Dt,EACT>::Face_handle
              pNeighbor = f->neighbor(i);
          i = pNeighbor->index(f);
          f = pNeighbor;
        }

        CGAL_triangulation_assertion((classify(f) ==
                                      Alpha_shape_2<Dt,EACT>::INTERIOR));

        CGAL_triangulation_assertion((classify(f, i) ==
                                      Alpha_shape_2<Dt,EACT>::REGULAR));

        os << V[f->vertex(f->ccw(i))] << ' '
                                      << V[f->vertex(f->cw(i))] << std::endl;
      }
    }
  }
  else
  { // get_mode() == GENERAL -----------------------------------------
    typename Alpha_shape_2<Dt,EACT>::Vertex_handle v;

    // write the regular vertices
    for (vertex_alpha_it = _interval_vertex_map.begin();
         vertex_alpha_it != _interval_vertex_map.end() &&
         (*vertex_alpha_it).first.first <= get_alpha();
         ++vertex_alpha_it)
    {
      pInterval2 = &(*vertex_alpha_it).first;

      if((pInterval2->second > get_alpha()
          || pInterval2->second == Infinity))
      {
        // alpha must be larger than the min boundary
        // and alpha is smaller than the upper boundary
        // which might be infinity
        // write the vertex

        v = (*vertex_alpha_it).second;
        CGAL_triangulation_assertion((classify(v) ==
                                      Alpha_shape_2<Dt,EACT>::REGULAR));
        V[v] = number_of_vertices++;
        os << v->point() << std::endl;
      }
    }

    // write the singular vertices
    for (;
         vertex_alpha_it != _interval_vertex_map.end();
         ++vertex_alpha_it)
    {
      v = (*vertex_alpha_it).second;
      CGAL_triangulation_assertion((classify(v) ==
                                    Alpha_shape_2<Dt,EACT>::SINGULAR));

      V[v] = number_of_vertices++;
      os << v->point() << std::endl;
    }

    // the vertices are oriented counterclockwise

    typename Alpha_shape_2<Dt,EACT>::Face_handle f;
    int i;

    for (edge_alpha_it = _interval_edge_map.begin();
         edge_alpha_it != _interval_edge_map.end() &&
         (*edge_alpha_it).first.first <= get_alpha();
         ++edge_alpha_it)
    {
      pInterval = &(*edge_alpha_it).first;

#ifdef CGAL_DEBUG_ALPHA_SHAPE_2
      typename Alpha_shape_2<Dt,EACT>::Type_of_alpha alpha =
          get_alpha();
      typename Alpha_shape_2<Dt,EACT>::Type_of_alpha alpha_min =
          pInterval->first;
      typename Alpha_shape_2<Dt,EACT>::Type_of_alpha alpha_mid =
          pInterval->second;
      typename Alpha_shape_2<Dt,EACT>::Type_of_alpha alpha_max =
          pInterval->third;
#endif // CGAL_DEBUG_ALPHA_SHAPE_2

      if(pInterval->third > get_alpha()
         || pInterval->third == Infinity)
      {
        // if alpha is smaller than the upper boundary
        // which might be infinity
        // visualize the boundary

        f = (*edge_alpha_it).second.first;
        i = (*edge_alpha_it).second.second;

        // write the regular edges
        if (pInterval->second != Infinity &&
            pInterval->second <= get_alpha())
        {
          CGAL_triangulation_assertion((classify(f, i) ==
                                        Alpha_shape_2<Dt,EACT>::REGULAR));
          // assure that all vertices are in ccw order
          if (classify(f) == Alpha_shape_2<Dt,EACT>::EXTERIOR)
          {
            // take the reverse face
            typename Alpha_shape_2<Dt,EACT>::Face_handle
                pNeighbor = f->neighbor(i);

            i = pNeighbor->index(f);
            f = pNeighbor;
          }

          CGAL_triangulation_assertion((classify(f) ==
                                        Alpha_shape_2<Dt,EACT>::INTERIOR));

          os << V[f->vertex(f->ccw(i))] << ' '
                                        << V[f->vertex(f->cw(i))] << std::endl;
        }
        else
        { // pInterval->second == Infinity ||
          //                           pInterval->second >= get_alpha())
          // pInterval->second == Infinity happens only for convex hull
          // of dimension 1 thus singular

          // write the singular edges
          if (pInterval->first != UNDEFINED)
          {
            CGAL_triangulation_assertion((classify(f, i) ==
                                          Alpha_shape_2<Dt,EACT>::SINGULAR));
            os << V[f->vertex(f->ccw(i))] << ' '
                                          << V[f->vertex(f->cw(i))] << std::endl;
          }
        }
      }
    }
  }
  return os;
}

//-------------------------------------------------------------------

template < class Dt, class EACT >
std::ostream& 
operator<<(std::ostream& os, const Alpha_shape_2<Dt>& A)
{
  return A.op_ostream(os);
}

//-------------------------------------------------------------------

template < class Dt, class EACT >
std::list<typename Alpha_shape_2<Dt,EACT>::Point_2>
Alpha_shape_2<Dt,EACT>::Output () 
{
  typename Interval_edge_map::const_iterator edge_alpha_it;

  const Interval3* pInterval;
  std::list<Point_2> L;

  if (get_mode() == REGULARIZED)
  {
    // it is much faster looking at the sorted intervals
    // than looking at all sorted faces
    // alpha must be larger than the mid boundary
    // and alpha is smaller than the upper boundary
    for (edge_alpha_it = _interval_edge_map.begin();
         edge_alpha_it != _interval_edge_map.end() &&
         (*edge_alpha_it).first.first <= get_alpha();
         ++edge_alpha_it)
    {
      pInterval = &(*edge_alpha_it).first;

      if (pInterval->second != Infinity)
      {
        // since this happens only for convex hull of dimension 1
        // thus singular

        if(pInterval->second <= get_alpha() &&
           (pInterval->third > get_alpha()
            || pInterval->third == Infinity))
        {
          // alpha must be larger than the mid boundary
          // and alpha is smaller than the upper boundary
          // which might be infinity
          // visualize the boundary

          CGAL_triangulation_assertion((classify((*edge_alpha_it).second.first,
                                                 (*edge_alpha_it).second.second) ==
                                        REGULAR));

          // if we used Edelsbrunner and Muecke's definition
          // regular means incident to a higher-dimensional face
          // thus we would write to many vertices
          L.push_back((this->segment((*edge_alpha_it).second.first,
                                     (*edge_alpha_it).second.second))
                      .source());
          L.push_back((this->segment((*edge_alpha_it).second.first,
                                     (*edge_alpha_it).second.second))
                      .target());
        }
      }
    }
  }
  else
  {  // get_mode() == GENERAL
    // draw the edges
    for (edge_alpha_it = _interval_edge_map.begin();
         edge_alpha_it != _interval_edge_map.end() &&
         (*edge_alpha_it).first.first <= get_alpha();
         ++edge_alpha_it)
    {
      pInterval = &(*edge_alpha_it).first;

      if (pInterval->first == UNDEFINED)
      {

        CGAL_triangulation_assertion(pInterval->second != Infinity);
        // since this happens only for convex hull of dimension 1
        // thus singular

        if(pInterval->second <= get_alpha() &&
           (pInterval->third > get_alpha()
            || pInterval->third == Infinity))
        {
          // alpha must be larger than the mid boundary
          // and alpha is smaller than the upper boundary
          // which might be infinity
          // visualize the boundary

          CGAL_triangulation_assertion((classify((*edge_alpha_it).second.first,
                                                 (*edge_alpha_it).second.second) ==
                                        REGULAR));
          L.push_back((this->segment((*edge_alpha_it).second.first,
                                     (*edge_alpha_it).second.second))
                      .source());
          L.push_back((this->segment((*edge_alpha_it).second.first,
                                     (*edge_alpha_it).second.second))
                      .target());
        }
      }
      else
      {
        if(pInterval->third > get_alpha()
           || pInterval->third == Infinity)
        {
          // if alpha is smaller than the upper boundary
          // which might be infinity
          // visualize the boundary

          CGAL_triangulation_assertion(((classify((*edge_alpha_it).second.first,
                                                  (*edge_alpha_it).second.second) ==
                                         REGULAR) ||
                                        (classify((*edge_alpha_it).second.first,
                                                  (*edge_alpha_it).second.second) ==
                                         SINGULAR)));
          L.push_back((this->segment((*edge_alpha_it).second.first,
                                     (*edge_alpha_it).second.second))
                      .source());
          L.push_back((this->segment((*edge_alpha_it).second.first,
                                     (*edge_alpha_it).second.second))
                      .target());
        }
      }
    }
  }
  return L;
}

template < class Dt, class EACT >
void Alpha_shape_2<Dt,EACT>::print_edge_map()
{
  for (typename Interval_edge_map::iterator iemapit= _interval_edge_map.begin();
       iemapit != _interval_edge_map.end(); ++iemapit) {
    Interval3 interval = (*iemapit).first;
    Edge edge = (*iemapit).second;
    Point p1 = edge.first->vertex(cw(edge.second))->point();
    Point p2 = edge.first->vertex(ccw(edge.second))->point();
    std::cout << "[ (" <<	p1 << ") - (" << p2 << ") ] :            "
              << interval.first << " "
              << interval.second << " " << interval.third << std::endl;
  }
}

} //namespace CGAL


#endif //CGAL_ALPHA_SHAPE_2_H
