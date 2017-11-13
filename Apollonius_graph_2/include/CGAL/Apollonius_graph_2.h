// Copyright (c) 2003,2004  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>



#ifndef CGAL_APOLLONIUS_GRAPH_2_H
#define CGAL_APOLLONIUS_GRAPH_2_H

#include <CGAL/license/Apollonius_graph_2.h>


#define CGAL_APOLLONIUS_GRAPH_PSEUDO_CIRCLE_DESIGN 1


#include <iostream>
#include <vector>
#include <map>

#include <boost/tuple/tuple.hpp>

#include <CGAL/Apollonius_graph_2/basic.h>

#include <CGAL/Triangulation_2.h>
#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/Triangulation_face_base_2.h>
#include <CGAL/Apollonius_graph_vertex_base_2.h>

#include <CGAL/in_place_edge_list.h>
#include <CGAL/internal/TDS_2/edge_list.h>
#include <CGAL/Apollonius_graph_2/Traits_wrapper_2.h>

#include <CGAL/Apollonius_graph_2/Constructions_C2.h>

#include <CGAL/iterator.h>
#include <CGAL/Iterator_project.h>
#include <CGAL/Nested_iterator.h>
#include <CGAL/Concatenate_iterator.h>


namespace CGAL {


namespace internal {

  template<typename Edge, typename LTag> struct AG2_which_list;

  // use the in-place edge list
  template<typename E>
  struct AG2_which_list<E,Tag_true>
  {
    typedef E                           Edge;
    typedef In_place_edge_list<Edge>    List;
  };

  // do not use the in-place edge list
  template<typename E>
  struct AG2_which_list<E,Tag_false>
  {
    typedef E                                 Edge;
    // change the following to Tag_false in order to use
    // CGAL's Unique_hash_map
    typedef Tag_true                          Use_stl_map_tag;
    typedef Edge_list<Edge,Use_stl_map_tag>   List;
  };

  template < class Node >
  struct Project_site_2 {
    typedef Node                   argument_type;
    typedef typename Node::Site_2  Site;
    typedef Site                   result_type;
    Site&       operator()( Node& x) const { return x.site(); }
    const Site& operator()( const Node& x) const { return x.site(); }
  };

} // namespace internal


template<class Gt,class Agds,class LTag>
class Apollonius_graph_hierarchy_2;

template < class Gt,
	   class Agds = Triangulation_data_structure_2 < 
               Apollonius_graph_vertex_base_2<Gt,true>,
               Triangulation_face_base_2<Gt> >,
	   class LTag = Tag_false>
class Apollonius_graph_2
  : private Triangulation_2
  <CGAL_APOLLONIUS_GRAPH_2_NS::Apollonius_graph_traits_wrapper_2<Gt>,Agds>
{
  friend class Apollonius_graph_hierarchy_2<Gt,Agds,LTag>;
private:
  // types and access methods needed for visualization
  //--------------------------------------------------

  // types
  typedef CGAL_APOLLONIUS_GRAPH_2_NS::Construct_Apollonius_bisector_2<Gt>
  Construct_Apollonius_bisector_2;

  typedef CGAL_APOLLONIUS_GRAPH_2_NS::Construct_Apollonius_bisector_ray_2<Gt>
  Construct_Apollonius_bisector_ray_2;

  typedef
  CGAL_APOLLONIUS_GRAPH_2_NS::Construct_Apollonius_bisector_segment_2<Gt>
  Construct_Apollonius_bisector_segment_2;

  typedef CGAL_APOLLONIUS_GRAPH_2_NS::Construct_Apollonius_primal_ray_2<Gt>
  Construct_Apollonius_primal_ray_2;

  typedef CGAL_APOLLONIUS_GRAPH_2_NS::Construct_Apollonius_primal_segment_2<Gt>
  Construct_Apollonius_primal_segment_2;


  // access
  Construct_Apollonius_bisector_2
  construct_Apollonius_bisector_2_object() const {
    return Construct_Apollonius_bisector_2();
  }

  Construct_Apollonius_bisector_ray_2
  construct_Apollonius_bisector_ray_2_object() const {
    return Construct_Apollonius_bisector_ray_2();
  }

  Construct_Apollonius_bisector_segment_2
  construct_Apollonius_bisector_segment_2_object() const { 
    return Construct_Apollonius_bisector_segment_2(); 
  }

  Construct_Apollonius_primal_ray_2
  construct_Apollonius_primal_ray_2_object() const {
    return Construct_Apollonius_primal_ray_2(); 
  }

  Construct_Apollonius_primal_segment_2
  construct_Apollonius_primal_segment_2_object() const { 
    return Construct_Apollonius_primal_segment_2();
  }


protected:
  // some local types
  typedef
  CGAL_APOLLONIUS_GRAPH_2_NS::Apollonius_graph_traits_wrapper_2<Gt>
  Modified_traits;

  typedef Triangulation_2<Modified_traits,Agds>  DG;

  typedef DG                                 Delaunay_graph;

public:
  // TYPES
  //------
  typedef Agds                                   Data_structure;
  typedef Agds                                   Triangulation_data_structure;
  typedef Gt                                     Geom_traits;
  typedef typename Gt::Point_2                   Point_2;
  typedef typename Gt::Site_2                    Site_2;

  typedef typename Agds::Edge                    Edge;
  typedef typename Agds::Vertex_handle           Vertex_handle;
  typedef typename Agds::Face_handle             Face_handle;
  typedef typename Agds::Vertex                  Vertex;
  typedef typename Agds::Face                    Face;

  typedef typename Agds::Vertex_circulator       Vertex_circulator;
  typedef typename Agds::Edge_circulator         Edge_circulator;
  typedef typename Agds::Face_circulator         Face_circulator;

  typedef typename Agds::Face_iterator           All_faces_iterator;
  typedef typename Agds::Vertex_iterator         All_vertices_iterator;
  typedef typename Agds::Edge_iterator           All_edges_iterator;

  typedef typename DG::Finite_faces_iterator     Finite_faces_iterator;
  typedef typename DG::Finite_vertices_iterator  Finite_vertices_iterator;
  typedef typename DG::Finite_edges_iterator     Finite_edges_iterator;

  typedef typename Agds::size_type               size_type;

  // Auxiliary iterators for convenience
  // do not use default template argument to please VC++
  typedef internal::Project_site_2<Vertex>                   Proj_site;
  typedef Iterator_project<Finite_vertices_iterator, 
                           Proj_site>
  /*                                           */ Visible_sites_iterator;

  typedef
  Apollonius_graph_vertex_base_nested_iterator_traits<
    Finite_vertices_iterator>  Hidden_sites_nested_iterator_traits;
  

  typedef Nested_iterator<Finite_vertices_iterator,
			  Hidden_sites_nested_iterator_traits>
  /*                                            */ Hidden_sites_iterator;

  typedef Concatenate_iterator<Visible_sites_iterator,
			       Hidden_sites_iterator>     Sites_iterator;

  typedef Site_2               value_type; // to have a back_inserter
  typedef const value_type&    const_reference; 
  typedef value_type&          reference;

public:
  struct   Vertex_iterator {};
  struct   Face_iterator {};
  struct   Edge_iterator {};

protected:
  // some more local types
  //  typedef typename Agds::Vertex            Vertex;

  // point lists
  typedef std::vector<Site_2>              Site_list;
  typedef typename Site_list::iterator     Site_list_iterator;

  typedef std::map<Face_handle,bool>           Face_map;
  typedef std::map<Face_handle, Face_handle>   Face_face_map;
  typedef std::map<Vertex_handle,bool>         Vertex_map;
  typedef std::set<Edge>                       Edge_list;

  typedef std::list<Vertex_handle>         Vertex_list;
  typedef typename Vertex_list::iterator   Vertex_list_iterator;
  typedef Vertex_handle                    Vh_triple[3];

  // the edge list
  typedef typename internal::AG2_which_list<Edge,LTag>::List  List;

  typedef enum { NO_CONFLICT = -1, INTERIOR, LEFT_VERTEX,
		 RIGHT_VERTEX, BOTH_VERTICES, ENTIRE_EDGE }
  Conflict_type;

  static Conflict_type opposite(const Conflict_type& ct) {
    if ( ct == RIGHT_VERTEX ) { return LEFT_VERTEX; }
    if ( ct == LEFT_VERTEX ) { return RIGHT_VERTEX; }
    return ct;
  }

protected:
  // Less_than comparator for site weights;
  // used to sort sites by decreasing weight when a sequence of sites
  // is inserted
  class Site_less_than_comparator
  {
  private:
    const Gt& gt;
  public:
    Site_less_than_comparator(const Gt& gt) : gt(gt) {}

    bool operator ()(const Site_2& p,
		     const Site_2& q) {
      Comparison_result result = gt.compare_weight_2_object()(p, q);
      return (result == LARGER);
    }
  };

public:
  // CREATION
  //---------
  Apollonius_graph_2(const Gt& gt=Gt()) :
    DG( Modified_traits(gt) ) {}

  template< class Input_iterator >
  Apollonius_graph_2(Input_iterator first, Input_iterator beyond,
		     const Gt& gt=Gt())
    : DG( Modified_traits(gt) )
  {
    insert(first, beyond);
  }

  Apollonius_graph_2(const Apollonius_graph_2 &ag)
    : DG(ag)
  {
    CGAL_postcondition( is_valid() );
  }

  Apollonius_graph_2&
  operator=(const Apollonius_graph_2& ag)
  {
    if ( this != &ag ) {
      DG::operator=(ag);
    }
    return (*this);
  }

public:
  // ACCESS METHODS
  // --------------
  const Geom_traits& geom_traits() const {
    return DG::geom_traits();
  }

  const Data_structure& data_structure() const { return this->_tds; }
  const Triangulation_data_structure& tds() const { return this->_tds; }

  int dimension() const {
    return this->_tds.dimension();
  }

  size_type number_of_faces() const {
    return this->_tds.number_of_faces();
  }

  size_type number_of_vertices() const {
    return DG::number_of_vertices();
  }

  size_type number_of_visible_sites() const {
    return number_of_vertices();
  }

  size_type number_of_hidden_sites() const {
    //    if ( !Vertex::StoreHidden ) { return 0; }

    size_type n_hidden(0);
    for (Finite_vertices_iterator vit = finite_vertices_begin();
	 vit != finite_vertices_end(); ++vit) {
      n_hidden += vit->number_of_hidden_sites();
    }

    return n_hidden;
  }

  Vertex_handle infinite_vertex() const {
    return DG::infinite_vertex();
  }

  Face_handle infinite_face() const {
    return DG::infinite_face();
  }

  Vertex_handle finite_vertex() const {
    return DG::finite_vertex();
  }

protected:
  using Delaunay_graph::cw;
  using Delaunay_graph::ccw;

public:
  // TRAVERSAL OF THE APOLLONIUS GRAPH
  //----------------------------------
  Finite_faces_iterator finite_faces_begin() const {
    return DG::finite_faces_begin();
  }

  Finite_faces_iterator finite_faces_end() const {
    return DG::finite_faces_end();
  }

  Finite_vertices_iterator finite_vertices_begin() const {
    return DG::finite_vertices_begin();
  }

  Finite_vertices_iterator finite_vertices_end() const {
    return DG::finite_vertices_end();
  }

  Finite_edges_iterator finite_edges_begin() const {
    return DG::finite_edges_begin();    
  }
  Finite_edges_iterator finite_edges_end() const {
    return DG::finite_edges_end();    
  }


  Sites_iterator sites_begin() const {
    return Sites_iterator(visible_sites_end(),
    			  hidden_sites_begin(),
    			  visible_sites_begin());
  }

  Sites_iterator sites_end() const {
    return Sites_iterator(visible_sites_end(),
    			  hidden_sites_begin(),
    			  hidden_sites_end(),0);
  }

  Visible_sites_iterator visible_sites_begin() const {
    return Visible_sites_iterator(finite_vertices_begin());
  }

  Visible_sites_iterator visible_sites_end() const {
    return Visible_sites_iterator(finite_vertices_end());
  }

  Hidden_sites_iterator hidden_sites_begin() const {
    return Hidden_sites_iterator(finite_vertices_end(),
				 finite_vertices_begin());
  }

  Hidden_sites_iterator hidden_sites_end() const {
    return Hidden_sites_iterator(finite_vertices_end(),
				 finite_vertices_end());
  }


  All_faces_iterator all_faces_begin() const {
    return DG::all_faces_begin();
  }

  All_faces_iterator all_faces_end() const {
    return DG::all_faces_end();
  }

  All_vertices_iterator all_vertices_begin() const {
    return DG::all_vertices_begin();
  }

  All_vertices_iterator all_vertices_end() const {
    return DG::all_vertices_end();
  }

  All_edges_iterator all_edges_begin() const {
    return DG::all_edges_begin();
  }

  All_edges_iterator all_edges_end() const {
    return DG::all_edges_end();
  }

public:
  // CIRCULATORS
  //------------
  Face_circulator
  incident_faces(Vertex_handle v,
		 Face_handle f = Face_handle()) const {
    return DG::incident_faces(v, f);
  }

  Vertex_circulator
  incident_vertices(Vertex_handle v,
		    Face_handle f = Face_handle()) const { 
    return DG::incident_vertices(v, f);
  }

  Edge_circulator
  incident_edges(Vertex_handle v,
		 Face_handle f = Face_handle()) const {
    return DG::incident_edges(v, f);
  }
 
public:
  // PREDICATES
  //-----------
  bool is_infinite(const Vertex_handle& v) const {
    return DG::is_infinite(v);
  }

  bool is_infinite(const Face_handle& f) const {
    return DG::is_infinite(f);
  }

  bool is_infinite(const Face_handle& f, int i) const {
    return DG::is_infinite(f, i);
  }

  bool is_infinite(const Edge& e) const {
    return is_infinite(e.first, e.second);
  }

  bool is_infinite(const Edge_circulator& ec) const {
    return DG::is_infinite(ec);
  }

public:
  // INSERTION
  //----------
  template< class Input_iterator >
  size_type insert(Input_iterator first, Input_iterator beyond) {
    // copy to a local container
    Site_list wp_list;
    for (Input_iterator it = first; it != beyond; ++it) {
      wp_list.push_back(*it);
    }

    // sort by decreasing weight
    Site_less_than_comparator less_than(geom_traits());
    std::sort(wp_list.begin(), wp_list.end(), less_than);

    // now insert
    Site_list_iterator lit;
    for (lit = wp_list.begin(); lit != wp_list.end(); ++lit) {
      insert(*lit);
    }

    // store how many sites where in the range
    size_type num = wp_list.size();

    // clear the local container
    wp_list.clear();

    // return the number of sites in range
    return num;
  }

  Vertex_handle  insert(const Site_2& p) {
    return insert(p, Vertex_handle());
  }

  Vertex_handle  insert(const Site_2& p, Vertex_handle vnear);

public:
  // REMOVAL
  //--------
  void remove(Vertex_handle v);

public:
  // NEAREST NEIGHBOR LOCATION
  //--------------------------
  Vertex_handle  nearest_neighbor(const Point_2& p) const;
  Vertex_handle  nearest_neighbor(const Point_2& p,
				  Vertex_handle vnear) const;

public:
  // ACCESS TO THE DUAL
  //-------------------
  typename Gt::Object_2 dual(const Face_handle& f) const;

  Site_2 dual(const Finite_faces_iterator& it) const
  {
    typename Gt::Object_2 o = dual(Face_handle(it));
    Site_2 s;
    if ( assign(s, o) ) {
      return s;
    } else {
      CGAL_assertion_code( bool the_assign_statement_must_always_work(false); )
      CGAL_assertion( the_assign_statement_must_always_work );
    }
    return s;
  }


private:
  typename Gt::Object_2 dual(const Edge e) const;

  typename Gt::Object_2 dual(const Edge_circulator& ec) const {
    return dual(*ec);
  }

  typename Gt::Object_2 dual(const Finite_edges_iterator& ei) const {
    return dual(*ei);
  }

public:
  // I/O
  //----
  void file_input(std::istream&);
  void file_output(std::ostream&) const;

  template< class Stream >
  Stream& draw_primal(Stream &str) const
  {
    if ( number_of_vertices() < 2 ) {
      // do nothing
    } else if ( number_of_vertices() == 2 ) {
      Vertex_handle v1(finite_vertices_begin());
      Vertex_handle v2(++finite_vertices_begin());
      Site_2 p1 = v1->site();
      Site_2 p2 = v2->site();
      

      typename Geom_traits::Segment_2 seg =
	construct_Apollonius_primal_segment_2_object()(p1,p2);
      typename Geom_traits::Ray_2 ray1 =
	construct_Apollonius_primal_ray_2_object()(p1,p2,p2);
      typename Geom_traits::Ray_2 ray2 =
	construct_Apollonius_primal_ray_2_object()(p2,p1,p1);

      str << seg;
      str << ray1;
      str << ray2;
    } else {
      All_edges_iterator eit = all_edges_begin();
      for (; eit != all_edges_end(); ++eit) {
	draw_primal_edge< Stream >(eit, str);
      }
    }
    return str;
  }

  template < class Stream > 
  Stream& draw_dual(Stream &str) const
  {
    Finite_edges_iterator eit = finite_edges_begin();
    for (; eit != finite_edges_end(); ++eit) {
      typename Gt::Object_2 o = dual(eit);
      typename Geom_traits::Line_2     l;
      typename Geom_traits::Segment_2  s;
      typename Geom_traits::Ray_2      r;
      CGAL::Hyperbola_2<Gt>            h;
      CGAL::Hyperbola_segment_2<Gt>    hs;
      CGAL::Hyperbola_ray_2<Gt>        hr;
      if (assign(hs, o)) hs.draw(str);
      else if (assign(s, o))  str << s; 
      else if (assign(hr, o))  hr.draw(str);
      else if (assign(r, o))   str << r;
      else if (assign(h, o))  h.draw(str);
      else if (assign(l, o)) str << l;
    }
    return str;
  }

protected:
  template< class Stream >
  Stream& draw_primal_vertex(const Finite_vertices_iterator& it,
			     Stream &str) const
  {
    return str << it->site().point();
  }


  template< class Stream >
  Stream& draw_dual_vertex(const Finite_faces_iterator& it,
			   Stream &str) const
  {
    return str << dual(it);
  }

public:
  template< class Stream >
  Stream& draw_primal_edge(const Finite_edges_iterator& eit,
			   Stream &str) const
  {
    return draw_primal_edge(*eit, str);
  }

  template< class Stream >
  Stream& draw_primal_edge(const All_edges_iterator& eit,
			   Stream &str) const
  {
    return draw_primal_edge(*eit, str);
  }


  template < class Stream > 
  Stream& draw_dual_edge(const Finite_edges_iterator& eit,
			 Stream &str) const
  {
    return draw_dual_edge(*eit, str);
  }

  template< class Stream >
  Stream& draw_primal_edge(const Edge& e, Stream &str) const
  {
    typedef typename Geom_traits::Segment_2     Segment_2;
    typedef typename Geom_traits::Ray_2         Ray_2;
    typedef std::pair<Segment_2,Segment_2>      Segment_pair_2;

    typename Geom_traits::Object_2 o = primal(e);
    Segment_2       s;
    Ray_2           r;
    Segment_pair_2  s_pair;
    CGAL::Hyperbola_segment_2<Gt>    hs;
    CGAL::Parabola_segment_2<Gt>     ps;
    if (assign(hs, o))  hs.draw(str);
    if (assign(s, o))   str << s; 
    if (assign(ps, o))  ps.draw(str);
    if (assign(r, o))   str << r;
    if (assign(s_pair, o)) str << s_pair.first << s_pair.second;
    return str;
  }

  template < class Stream > 
  Stream& draw_dual_edge(const Edge& e, Stream &str) const
  {
    if ( is_infinite(e) ) { return str; }
    typename Gt::Object_2 o = dual(e);
    typename Geom_traits::Line_2     l;
    typename Geom_traits::Segment_2  s;
    typename Geom_traits::Ray_2      r;
    CGAL::Hyperbola_2<Gt>            h;
    CGAL::Hyperbola_segment_2<Gt>    hs;
    CGAL::Hyperbola_ray_2<Gt>        hr;
    if (assign(hs, o))  hs.draw(str);
    if (assign(s, o))   str << s; 
    if (assign(hr, o))  hr.draw(str);
    if (assign(r, o))   str << r;
    if (assign(h, o))   h.draw(str);
    if (assign(l, o))   str << l;

    return str;
  }

protected:
  template< class Stream >
  Stream& draw_primal_face(All_faces_iterator fit, Stream &str) const
  {
    for (int i = 0; i < 3; i++) {
      draw_primal_edge< Stream >(Edge(Face_handle(fit), i), str);
    }
    return str;
  }


  template< class Stream >
  Stream& draw_dual_face(const All_vertices_iterator& vit,
			 Stream &str) const
  {
    Edge_circulator ec_start = incident_edges(Vertex_handle(vit));
    Edge_circulator ec = ec_start;
    do {
      draw_dual_edge< Stream >(*ec, str);
      ++ec;
    } while ( ec_start != ec );
    return str;
  }

protected:
  template < class Stream > 
  Stream& draw_dual_sites(Stream &str) const
  {
    All_faces_iterator fit = all_faces_begin();
    for (; fit != all_faces_end(); ++fit) {
      Face_handle f(fit);
      if ( is_infinite(f) ) {
	if (  is_infinite(f->vertex(0))  ) {
	  str << circumcircle( f->vertex(1)->site(),
			       f->vertex(2)->site() );
	} else if (  is_infinite(f->vertex(1))  ){
	  str << circumcircle( f->vertex(2)->site(),
			       f->vertex(0)->site() );
	} else {
	  str << circumcircle( f->vertex(0)->site(),
			       f->vertex(1)->site() );	  
	}
      } else {
	Site_2 wp = circumcircle(f);
	typename Gt::Rep::Circle_2 c(wp.point(),
				     CGAL::square(wp.weight()));
	str << c;
      }
    }
    return str;
  }

public:
  // VALIDITY CHECK
  //---------------
  bool is_valid(bool verbose = false, int level = 1) const;

public:
  // MISCELLANEOUS
  //--------------
  void clear() {
    DG::clear();
  }

  void swap(Apollonius_graph_2& ag) {
    DG::swap(ag);
  }

public:
  // MK: THE FOLLOWING ARE NOT IN THE SPEC
  //======================================
  // Primal
  typename Gt::Object_2 primal(const Edge e) const;
  typename Gt::Object_2 primal(const Edge_circulator& ec) const {
    return primal(*ec);
  }
  typename Gt::Object_2 primal(const Finite_edges_iterator& ei) const {
    return primal(*ei);
  }

protected:
  // wrappers for the geometric predicates

  // checks is q is contained inside p
  bool is_hidden(const Site_2 &p, const Site_2 &q) const;

  // returns:
  //   ON_POSITIVE_SIDE if q is closer to p1
  //   ON_NEGATIVE_SIDE if q is closer to p2
  //   ON_ORIENTED_BOUNDARY if q is on the bisector of p1 and p2
  Oriented_side side_of_bisector(const Site_2 &p1,
				 const Site_2 &p2,
				 const Point_2 &q) const;

  Sign incircle(const Site_2 &p1, const Site_2 &p2,
		const Site_2 &p3, const Site_2 &q) const;

  Sign incircle(const Site_2 &p1, const Site_2 &p2,
		const Site_2 &q) const;


  Sign incircle(const Face_handle& f, const Site_2& q) const;


  Sign incircle(const Vertex_handle& v0, const Vertex_handle& v1,
		const Vertex_handle& v) const;

  Sign incircle(const Vertex_handle& v0, const Vertex_handle& v1,
		const Vertex_handle& v2, const Vertex_handle& v) const;


  
  bool finite_edge_interior(const Site_2& p1,
			    const Site_2& p2,
			    const Site_2& p3,
			    const Site_2& p4,
			    const Site_2& q,
			    bool endpoints_in_conflict) const;

  bool finite_edge_interior(const Face_handle& f, int i,
			    const Site_2& q,
			    bool endpoints_in_conflict) const;

  bool finite_edge_interior(const Vertex_handle& v1,
			    const Vertex_handle& v2,
			    const Vertex_handle& v3,
			    const Vertex_handle& v4,
			    const Vertex_handle& v,
			    bool endpoints_in_conflict) const;

  bool finite_edge_interior_degenerated(const Site_2& p1,
					const Site_2& p2,
					const Site_2& p3,
					const Site_2& q,
					bool endpoints_in_conflict) const;


  bool finite_edge_interior_degenerated(const Site_2& p1,
					const Site_2& p2,
					const Site_2& q,
					bool endpoints_in_conflict) const;

  bool finite_edge_interior_degenerated(const Face_handle& f, int i,
					const Site_2& p,
					bool endpoints_in_conflict) const;

  bool finite_edge_interior_degenerated(const Vertex_handle& v1,
					const Vertex_handle& v2,
					const Vertex_handle& v3,
					const Vertex_handle& v4,
					const Vertex_handle& v,
					bool endpoints_in_conflict) const;
  bool infinite_edge_interior(const Site_2& p2,
			      const Site_2& p3,
			      const Site_2& p4,
			      const Site_2& q,
			      bool endpoints_in_conflict) const;


  bool infinite_edge_interior(const Face_handle& f, int i,
			      const Site_2& p,
			      bool endpoints_in_conflict) const;

  bool infinite_edge_interior(const Vertex_handle& v1,
			      const Vertex_handle& v2,
			      const Vertex_handle& v3,
			      const Vertex_handle& v4,
			      const Vertex_handle& v,
			      bool endpoints_in_conflict) const;

  Conflict_type
  infinite_edge_conflict_type(const Site_2& p2,
			      const Site_2& p3,
			      const Site_2& p4,
			      const Site_2& q) const;

  Conflict_type
  finite_edge_conflict_type_degenerated(const Site_2& p1,
					const Site_2& p2,
					const Site_2& q) const;

  bool edge_interior(const Face_handle& f, int i,
		     const Site_2& p, bool b) const;


  bool edge_interior(const Edge& e,
			    const Site_2& p, bool b) const {
    return edge_interior(e.first, e.second, p, b);
  }

  bool edge_interior(const Vertex_handle& v1,
		     const Vertex_handle& v2,
		     const Vertex_handle& v3,
		     const Vertex_handle& v4,
		     const Vertex_handle& v,
		     bool endpoints_in_conflict) const;

  bool is_degenerate_edge(const Site_2& p1,
				 const Site_2& p2,
				 const Site_2& p3,
				 const Site_2& p4) const {
    return geom_traits().is_degenerate_edge_2_object()
      (p1, p2, p3, p4);
  }

  bool is_degenerate_edge(const Vertex_handle& v1,
				 const Vertex_handle& v2,
				 const Vertex_handle& v3,
				 const Vertex_handle& v4) const {
    CGAL_precondition( !is_infinite(v1) && !is_infinite(v2) &&
		       !is_infinite(v3) && !is_infinite(v4) );

    return is_degenerate_edge(v1->site(), v2->site(),
			      v3->site(), v4->site());
  }

  bool is_degenerate_edge(const Face_handle& f, int i) const {
    Vertex_handle v1 = f->vertex( ccw(i) );
    Vertex_handle v2 = f->vertex(  cw(i) );
    Vertex_handle v3 = f->vertex(     i  );
    Vertex_handle v4 = tds().mirror_vertex(f, i);

    return is_degenerate_edge(v1, v2, v3, v4);
  }

  bool is_degenerate_edge(const Edge& e) const {
    return is_degenerate_edge(e.first, e.second);
  }


protected:
  // wrappers for constructions
  Point_2 circumcenter(const Face_handle& f) const;
  Point_2 circumcenter(const Site_2& p0, 
		       const Site_2& p1, 
		       const Site_2& p2) const;

  Site_2 circumcircle(const Face_handle& f) const;
  Site_2 circumcircle(const Site_2& p0, 
		      const Site_2& p1, 
		      const Site_2& p2) const;

  typename Gt::Line_2 circumcircle(const Site_2& p0,
				   const Site_2& p1) const;

protected:
  // wrappers for combinatorial operations on the data structure

  // getting the degree of a vertex
  typename Data_structure::size_type degree(const Vertex_handle& v) {
    return this->_tds.degree(v);
  }

  // getting the symmetric edge
  Edge sym_edge(const Edge e) const {
    return sym_edge(e.first, e.second);
  }

  Edge sym_edge(const Face_handle& f, int i) const {
    Face_handle f_sym = f->neighbor(i);
    return Edge(  f_sym, f_sym->index( tds().mirror_vertex(f, i) )  );
  }

  Edge flip(Face_handle& f, int i);
  Edge flip(Edge e);

  Vertex_handle insert_in_face(Face_handle& f, const Site_2& p);

  bool          is_degree_2(const Vertex_handle& v) const;

  Vertex_handle insert_degree_2(Edge e);
  Vertex_handle insert_degree_2(Edge e, const Site_2& p);
  void          remove_degree_2(Vertex_handle v);
  void          remove_degree_3(Vertex_handle v);
  void          remove_degree_3(Vertex_handle v, Face_handle f);

  // this was defined because the hierarchy needs it
  Vertex_handle create_vertex() {
    return this->_tds.create_vertex();
  }

protected:
  // insertion of the first three sites
  Vertex_handle  insert_first(const Site_2& p);
  Vertex_handle  insert_second(const Site_2& p);
  Vertex_handle  insert_third(const Site_2& p);

  // methods for insertion
  void initialize_conflict_region(const Face_handle& f, List& l) const;
  bool check_edge_for_hidden_sites(const Face_handle& f, int i,
				   const Site_2& p, Vertex_map& vm) const;
  void expand_conflict_region(const Face_handle& f,
			      const Site_2& p,
			      List& l, Face_map& fm, Vertex_map& vm,
			      std::vector<Vh_triple*>* fe);

  Vertex_handle add_bogus_vertex(Edge e, List& l);
  Vertex_list   add_bogus_vertices(List& l);
  void          remove_bogus_vertices(Vertex_list& vl);

  void move_hidden_sites(Vertex_handle& vold, Vertex_handle& vnew);

  // MK: this is not currently used
  std::vector<Face*> get_faces_for_recycling(Face_map& fm,
					     unsigned int n_wanted);
  void remove_hidden_vertices(Vertex_map& vm);
  Vertex_handle retriangulate_conflict_region(const Site_2& p,
					      List& l,
					      Face_map& fm,
					      Vertex_map& vm);

protected:
  // methods for removal
  void  remove_first(Vertex_handle v);
  void  remove_second(Vertex_handle v);
  void  remove_third(Vertex_handle v);
  void  remove_degree_d_vertex(Vertex_handle v);
  void  minimize_degree(Vertex_handle v);

  void find_conflict_region_remove(const Vertex_handle& v,
				   const Vertex_handle& vnearest,
				   List& l, Face_map& fm,
				   Vertex_map& vm,
				   std::vector<Vh_triple*>* fe);

protected:
  // methods for I/O

  template<class T>
  bool assign(T& t2, const typename Gt::Object_2& o2) const
  {
    return geom_traits().assign_2_object()(t2, o2);
  }


protected:
  template<class OutputItFaces>
  OutputItFaces find_conflicts(const Face_handle& f,
			       const Site_2& p,
			       List& l,
			       Face_map& fm,
			       Vertex_map& vm,
			       OutputItFaces fit) const
  {
    // setting fm[f] to true means that the face has been reached and
    // that the face is available for recycling. If we do not want the
    // face to be available for recycling we must set this flag to
    // false.
    if ( fm.find(f) != fm.end() ) { return fit; }
    fm[f] = true;
    CGAL_assertion( incircle(f, p) == NEGATIVE );
    *fit++ = f;

    //  CGAL_assertion( fm.find(f) != fm.end() );
    for (int i = 0; i < 3; i++) {
      bool hidden_found = check_edge_for_hidden_sites(f, i, p, vm);

      Face_handle n = f->neighbor(i);

      if ( !hidden_found ) {
	Sign s = incircle(n, p);
	if ( s != NEGATIVE ) { continue; }

	bool interior_in_conflict = edge_interior(f, i, p, true);

	if ( !interior_in_conflict ) { continue; }
      }

      if ( fm.find(n) != fm.end() ) {
	Edge e = sym_edge(f, i);
	if ( l.is_in_list(e) ||
	     l.is_in_list(sym_edge(e)) ) {
	  l.remove(e);
	  l.remove(sym_edge(e));
	}
	continue;
      }

      Edge e = sym_edge(f, i);

      CGAL_assertion( l.is_in_list(e) );
      int j = tds().mirror_index(f, i);
      Edge e_before = sym_edge(n, ccw(j));
      Edge e_after = sym_edge(n, cw(j));
      if ( !l.is_in_list(e_before) ) {
	l.insert_before(e, e_before);
      }
      if ( !l.is_in_list(e_after) ) {
	l.insert_after(e, e_after);
      }
      l.remove(e);

      fit = find_conflicts(n, p, l, fm, vm, fit);
    } // for-loop
    return fit;
  } // find_conflicts

  bool equal(const Edge& e1, const Edge& e2) const {
    return e1.first == e2.first && e1.second == e2.second;
  }


protected:
  template<class OutputItFaces, class OutputItBoundaryEdges,
	   class OutputItHiddenVertices>
  boost::tuples::tuple<OutputItFaces, OutputItBoundaryEdges,
		       OutputItHiddenVertices>
  get_all(const Site_2& p,
	  OutputItFaces fit,
	  OutputItBoundaryEdges eit,
	  OutputItHiddenVertices vit,
	  Vertex_handle start,
	  bool find_nearest) const
  {
    CGAL_precondition( dimension() == 2 );

    // first find the nearest neighbor
    Vertex_handle vnearest = start;
    if ( find_nearest ) {
      vnearest = nearest_neighbor(p.point(), start);
      CGAL_assertion( vnearest != Vertex_handle() );
    }

    // check if it is hidden
    if ( is_hidden(vnearest->site(), p) ) {
      return boost::tuples::make_tuple(fit, eit, vit);
    }

    // find the first conflict

    // first look for conflict with vertex
    Face_circulator fc_start = incident_faces(vnearest);
    Face_circulator fc = fc_start;
    Face_handle start_f;
    Sign s;
    do {
      Face_handle f(fc);
      s = incircle(f, p);

      if ( s == NEGATIVE ) {
	start_f = f;
	break;
      }
      ++fc;
    } while ( fc != fc_start );

    // we are not in conflict with an Apollonius vertex, so we have to
    // be in conflict with the interior of an Apollonius edge
    if ( s != NEGATIVE ) {
      Edge_circulator ec_start = incident_edges(vnearest);
      Edge_circulator ec = ec_start;

      bool interior_in_conflict(false);
      Edge e;
      do {
	e = *ec;
	interior_in_conflict = edge_interior(e, p, false);
	
	if ( interior_in_conflict ) { break; }
	++ec;
      } while ( ec != ec_start );

      CGAL_assertion( interior_in_conflict );

      *eit++ = e;
      *eit++ = sym_edge(e);
      return boost::tuples::make_tuple(fit, eit, vit);
    }

    // we are in conflict with an Apollonius vertex; start from that and 
    // find the entire conflict region and then repair the diagram
    List l;
    Face_map fm;
    Vertex_map vm;

    //    *fit++ = start_f;
    initialize_conflict_region(start_f, l);
    fit = find_conflicts(start_f, p, l, fm, vm, fit);

    // output the edges on the boundary of the conflict region
    if ( l.size() > 0 ) {
      const Edge& e_front = l.front();
      // here I should be able to write: const Edge& e = l.front();
      // instead of what I have; but the compiler complains for the
      // assignment: e = l.next(e);
      Edge e = l.front();
      do {
	*eit++ = e;
	e = l.next(e);
      } while ( !equal(e, e_front) );
    }

    // output the hidden vertices
    for (typename Vertex_map::iterator it = vm.begin(); it != vm.end(); ++it) {
      *vit++ = it->first;
    }

    // clear containers
    fm.clear();
    vm.clear();
    l.clear();

    return boost::tuples::make_tuple(fit, eit, vit);
  }

public:
  template<class OutputItFaces, class OutputItBoundaryEdges,
	   class OutputItHiddenVertices>
  boost::tuples::tuple<OutputItFaces, OutputItBoundaryEdges,
		       OutputItHiddenVertices>
  get_conflicts_and_boundary_and_hidden_vertices(const Site_2& p,
						 OutputItFaces fit,
						 OutputItBoundaryEdges eit,
						 OutputItHiddenVertices vit,
						 Vertex_handle start =
						 Vertex_handle()) const
  {
    return get_all(p, fit, eit, vit, start, true);
  }

  template<class OutputItFaces, class OutputItBoundaryEdges>
  std::pair<OutputItFaces, OutputItBoundaryEdges>
  get_conflicts_and_boundary(const Site_2& p,
			     OutputItFaces fit,
			     OutputItBoundaryEdges eit,
			     Vertex_handle start =
			     Vertex_handle()) const {
    boost::tuples::tuple<OutputItFaces,OutputItBoundaryEdges,Emptyset_iterator>
      tup =
      get_conflicts_and_boundary_and_hidden_vertices(p,
						     fit,
						     eit,
						     Emptyset_iterator(),
						     start);
    return std::make_pair( boost::tuples::get<0>(tup),
			   boost::tuples::get<1>(tup) );
  }


  template<class OutputItBoundaryEdges, class OutputItHiddenVertices>
  std::pair<OutputItBoundaryEdges, OutputItHiddenVertices>
  get_boundary_of_conflicts_and_hidden_vertices(const Site_2& p,
						OutputItBoundaryEdges eit,
						OutputItHiddenVertices vit,
						Vertex_handle start =
						Vertex_handle()) const {
    boost::tuples::tuple<Emptyset_iterator,OutputItBoundaryEdges,
      OutputItHiddenVertices>
      tup =
      get_conflicts_and_boundary_and_hidden_vertices(p,
						     Emptyset_iterator(),
						     eit,
						     vit,
						     start);
    return std::make_pair( boost::tuples::get<1>(tup),
			   boost::tuples::get<2>(tup) );
  }


  template<class OutputItFaces, class OutputItHiddenVertices>
  std::pair<OutputItFaces, OutputItHiddenVertices>
  get_conflicts_and_hidden_vertices(const Site_2& p,
				    OutputItFaces fit,
				    OutputItHiddenVertices vit,
				    Vertex_handle start =
				    Vertex_handle()) const {
    boost::tuples::tuple<OutputItFaces,Emptyset_iterator,
      OutputItHiddenVertices>
      tup =
      get_conflicts_and_boundary_and_hidden_vertices(p,
						     fit,
						     Emptyset_iterator(),
						     vit,
						     start);
    return std::make_pair( boost::tuples::get<0>(tup),
			   boost::tuples::get<2>(tup) );
  }

  template<class OutputItFaces>
  OutputItFaces get_conflicts(const Site_2& p,
			      OutputItFaces fit,
			      Vertex_handle start = Vertex_handle()) const {
    boost::tuples::tuple<OutputItFaces,Emptyset_iterator,Emptyset_iterator>
      tup =
      get_conflicts_and_boundary_and_hidden_vertices(p,
						     fit,
						     Emptyset_iterator(),
						     Emptyset_iterator(),
						     start);
    return boost::tuples::get<0>(tup);
  }

  template<class OutputItBoundaryEdges>
  OutputItBoundaryEdges
  get_boundary_of_conflicts(const Site_2& p,
			    OutputItBoundaryEdges eit,
			    Vertex_handle start = Vertex_handle()) const {
    boost::tuples::tuple<Emptyset_iterator,OutputItBoundaryEdges,
      Emptyset_iterator>
      tup =
      get_conflicts_and_boundary_and_hidden_vertices(p,
						     Emptyset_iterator(),
						     eit,
						     Emptyset_iterator(),
						     start);
    return boost::tuples::get<1>(tup);
  }

  template<class OutputItHiddenVertices>
  OutputItHiddenVertices
  get_hidden_vertices(const Site_2& p,
		      OutputItHiddenVertices vit,
		      Vertex_handle start = Vertex_handle()) const {
    boost::tuples::tuple<Emptyset_iterator,Emptyset_iterator,
      OutputItHiddenVertices>
      tup =
      get_conflicts_and_boundary_and_hidden_vertices(p,
						     Emptyset_iterator(),
						     Emptyset_iterator(),
						     vit,
						     start);
    return boost::tuples::get<2>(tup);
  }


}; // Apollonius_graph_2

template<class Gt, class Agds, class LTag>
std::ostream& operator<<(std::ostream& os,
			 const Apollonius_graph_2<Gt,Agds,LTag>& ag)
{
  ag.file_output(os);
  return os;
}

template<class Gt, class Agds, class LTag>
std::istream& operator>>(std::istream& is,
			 Apollonius_graph_2<Gt,Agds,LTag>& ag)
{
  ag.file_input(is);
  return is;
}

} //namespace CGAL

#include <CGAL/Apollonius_graph_2/Apollonius_graph_2_impl.h>


#endif // CGAL_APOLLONIUS_GRAPH_2_H
