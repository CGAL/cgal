// Copyright (c) 2003,2004  INRIA Sophia-Antipolis (France) and
// Notre Dame University (U.S.A.).  All rights reserved.
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
// Author(s)     : Menelaos Karavelas <mkaravel@cse.nd.edu>



#ifndef CGAL_SEGMENT_VORONOI_DIAGRAM_2_H
#define CGAL_SEGMENT_VORONOI_DIAGRAM_2_H

#include <vector>
#include <map>
#include <algorithm>

#include <CGAL/Segment_Voronoi_diagram_short_names_2.h>

#include <CGAL/Triangulation_2.h>
#include <CGAL/Segment_Voronoi_diagram_site_2.h>
#include <CGAL/Segment_Voronoi_diagram_data_structure_2.h>
#include <CGAL/Triangulation_face_base_2.h>
#include <CGAL/Segment_Voronoi_diagram_vertex_base_2.h>

#include <CGAL/Segment_Voronoi_diagram_constructions_C2.h>

#include <CGAL/in_place_edge_list.h>
#include <CGAL/edge_list.h>
#include <CGAL/Segment_Voronoi_diagram_traits_wrapper_2.h>

#include <CGAL/Iterator_project.h>
#include <CGAL/Simple_container_wrapper.h>


/*
  Conventions:
  ------------
  1. we treat segments as open; the endpoints are separate objects
  2. a segment of length zero is treated as a point
  3. a point is deleted only if it has no segment adjacent to it
  4. when a segment is deleted it's endpoints are not deleted
  5. the user can force the deletion of endpoints; this is only done
     if condition 3 is met.
  6. when objects are written to a stream we distinguish between
     points and segments; points start by a 'p' and segments by an 's'.
*/

//#define STORE_INPUT_SITES 1
#define STORE_INPUT_SITES 0

CGAL_BEGIN_NAMESPACE

namespace CGALi {

  template<typename Edge, typename LTag> struct SVD_which_list;

  // use the in-place edge list
  template<typename E>
  struct SVD_which_list<E,Tag_true>
  {
    typedef E                           Edge;
    typedef In_place_edge_list<Edge>    List;
  };

  // do not use the in-place edge list
  template<typename E>
  struct SVD_which_list<E,Tag_false>
  {
    typedef E                                 Edge;
    // change the following to Tag_false in order to use
    // CGAL's Unique_hash_map
    typedef Tag_true                          Use_stl_map_tag;
    typedef Edge_list<Edge,Use_stl_map_tag>   List;
  };


  template < class Node >
  struct Svd_project_site_2 {
    typedef Node                   argument_type;
    typedef typename Node::Site_2  Site;
    typedef Site                   result_type;
    Site& operator()(const Node& x) const { 
      static Site s;
      s = x.site();
      return s;
    }
    //    const Site& operator()(const Node& x) const { return x.site(); }
  };

  template < class Node, class Site_t >
  struct Svd_project_input_to_site_2 {
    typedef Node                   argument_type;
    typedef Site_t                 Site;
    typedef Site                   result_type;
    Site& operator()(const Node& x) const {
      static Site s;
      if ( x.third ) { // it is a point
	s = Site::construct_site_2(x.first);
      } else {
	s = Site::construct_site_2(x.first, x.second);
      }
      return s;
    }
  };

} // namespace CGALi


template<class Gt, class STag, class DS, class LTag >
class Segment_Voronoi_diagram_hierarchy_2;



template<class Gt,
	 class DS = Segment_Voronoi_diagram_data_structure_2 < 
                Segment_Voronoi_diagram_vertex_base_2<Gt,
			    typename Gt::Intersections_tag>,
                Triangulation_face_base_2<Gt> >,
	 class LTag = Tag_false >
class Segment_Voronoi_diagram_2
  : private Triangulation_2<
          Segment_Voronoi_diagram_traits_wrapper_2<Gt>, DS >
{
  friend class Segment_Voronoi_diagram_hierarchy_2<Gt,Tag_true,DS,LTag>;
  friend class Segment_Voronoi_diagram_hierarchy_2<Gt,Tag_false,DS,LTag>;
protected:
  // LOCAL TYPES
  //------------
  typedef Segment_Voronoi_diagram_traits_wrapper_2<Gt>  Modified_traits;
  typedef Triangulation_2<Modified_traits,DS>           DG;
  typedef DG                                            Delaunay_graph;
  typedef typename DG::Vertex                           Vertex;
  typedef typename DG::Face                             Face;

  typedef LTag                                          List_tag;

public:
  // PUBLIC TYPES
  //-------------
  typedef DS                                     Data_structure;
  typedef Gt                                     Geom_traits;
  typedef typename Gt::Site_2                    Site_2;
  typedef typename Gt::Point_2                   Point_2;

  typedef typename DG::Edge                      Edge;
  typedef typename DG::Vertex_handle             Vertex_handle;
  typedef typename DG::Face_handle               Face_handle;

  typedef typename DG::Vertex_circulator         Vertex_circulator;
  typedef typename DG::Edge_circulator           Edge_circulator;
  typedef typename DG::Face_circulator           Face_circulator;

  typedef typename DG::All_faces_iterator        All_faces_iterator;
  typedef typename DG::Finite_faces_iterator     Finite_faces_iterator;
  typedef typename DG::All_vertices_iterator     All_vertices_iterator;
  typedef typename DG::Finite_vertices_iterator  Finite_vertices_iterator;
  typedef typename DG::All_edges_iterator        All_edges_iterator;
  typedef typename DG::Finite_edges_iterator     Finite_edges_iterator;

private:
  typedef typename Geom_traits::Arrangement_type_2  AT2;
  typedef typename AT2::Arrangement_type            Arrangement_type;

  typedef std::list<Point_2>                     PC;

#if STORE_INPUT_SITES
  typedef std::vector<Site_2>        Input_sites_container;
  typedef typename Input_sites_container::const_iterator
  All_inputs_iterator;
#else
  // these containers should have point handles and should replace the
  // point container...
  typedef Triple<Point_2,Point_2,bool>   Site_rep_2;
  typedef std::vector<Site_rep_2>        Input_sites_container;
  typedef typename Input_sites_container::const_iterator
  All_inputs_iterator;

  typedef CGALi::Svd_project_input_to_site_2<Site_rep_2, Site_2>
  Proj_input_to_site;
#endif

public:
  typedef Simple_container_wrapper<PC>           Point_container;
  typedef typename Point_container::iterator     Point_handle;

  typedef typename DG::size_type                 size_type;

private:
  typedef CGALi::Svd_project_site_2<Vertex>      Proj_site;

public:
#if STORE_INPUT_SITES
  typedef All_inputs_iterator Input_sites_iterator;
#else
  typedef Iterator_project<All_inputs_iterator, Proj_input_to_site>
  Input_sites_iterator;
#endif

  typedef Iterator_project<Finite_vertices_iterator, 
                           Proj_site>            Output_sites_iterator;
protected:
  // LOCAL VARIABLE(S)
  //------------------
  // the container of points
  Point_container pc_;
  Input_sites_container isc_;

protected:
  // MORE LOCAL TYPES
  //-----------------
  typedef typename Gt::Intersections_tag        Intersections_tag;

  typedef typename Data_structure::Vertex_base  Vertex_base;

  typedef std::map<Face_handle,bool>            Face_map;
  typedef std::vector<Edge>                     Edge_vector;

  typedef std::list<Vertex_handle>              Vertex_list;
  typedef typename Vertex_list::iterator        Vertex_list_iterator;

  typedef Triple<Vertex_handle,Vertex_handle,Vertex_handle>
  Vertex_triple;

  typedef std::pair<Face_handle,Face_handle>    Face_pair;

  typedef typename Vertex_base::Storage_site_2  Storage_site_2;

  // the edge list
  typedef typename CGALi::SVD_which_list<Edge,List_tag>::List  List;

public:
  // CREATION
  //---------
  Segment_Voronoi_diagram_2(const Gt& gt=Gt()) : DG(gt) {}

  template< class Input_iterator >
  Segment_Voronoi_diagram_2(Input_iterator first, Input_iterator beyond,
			    const Gt& gt=Gt())
    : DG(gt)
  {
    insert(first, beyond);
  }

  Segment_Voronoi_diagram_2(const Segment_Voronoi_diagram_2 &svd)
    : DG(svd)
  {
    CGAL_postcondition( is_valid() );
  }

  Segment_Voronoi_diagram_2&
  operator=(const Segment_Voronoi_diagram_2& svd)
  {
    DG::operator=(svd);
    return (*this);
  }

public:
  // ACCESS METHODS
  // --------------
  const Geom_traits&  geom_traits() const { return DG::geom_traits(); }

  const Data_structure&   ds() const { return this->_tds; }
  const Point_container&  point_container() const { return pc_; }

  size_type number_of_input_sites() const {
    return isc_.size();
  }

  size_type number_of_output_sites() const {
    return number_of_vertices();
  }

  size_type number_of_vertices() const {
    return DG::number_of_vertices();
  }

  size_type number_of_faces() const {
    return DG::number_of_faces();
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

  using Delaunay_graph::cw;
  using Delaunay_graph::ccw;

public:
  // TRAVERSAL OF THE DUAL GRAPH
  //----------------------------
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

  Input_sites_iterator input_sites_begin() const {
#if STORE_INPUT_SITES
    return isc_.begin();
#else
    return Input_sites_iterator(isc_.begin());
#endif
  }

  Input_sites_iterator input_sites_end() const {
#if STORE_INPUT_SITES
    return isc_.end();
#else
    return Input_sites_iterator(isc_.end());
#endif
  }

  Output_sites_iterator output_sites_begin() const {
    return Output_sites_iterator(finite_vertices_begin());
  }

  Output_sites_iterator output_sites_end() const {
    return Output_sites_iterator(finite_vertices_end());    
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
  // INSERTION METHODS
  //------------------
  template<class Input_iterator>
  size_type insert(Input_iterator first, Input_iterator beyond) {
    return insert(first, beyond, Tag_false());
  }

  template<class Input_iterator>
  size_type insert(Input_iterator first, Input_iterator beyond, Tag_true)
  {
    std::vector<Site_2> site_vec;
    for (Input_iterator it = first; it != beyond; ++it) {
      site_vec.push_back(Site_2(*it));
    }
    std::random_shuffle(site_vec.begin(), site_vec.end());
    return insert(site_vec.begin(), site_vec.end(), Tag_false());
  }

  template<class Input_iterator>
  size_type insert(Input_iterator first, Input_iterator beyond,	Tag_false)
  {
    // do it the obvious way: insert them as they come;
    // one might think though that it might be better to first insert
    // all end points and then all segments, or a variation of that.
    size_type n_before = number_of_vertices();
    for (Input_iterator it = first; it != beyond; ++it) {
      insert(*it);
    }
    size_type n_after = number_of_vertices();
    return n_after - n_before;
  }

  // insert a point
  Vertex_handle  insert(const Point_2& p) {
    // update input site container
    register_input_site(p);
    return insert_point(p, Vertex_handle());
  }

  Vertex_handle  insert(const Point_2& p, Vertex_handle vnear) {
    // update input site container
    register_input_site(p);
    return insert_point(p, vnear);
  }

protected:
  // insert a point without registering it in the input sites
  // container: useful for the hierarchy
  Vertex_handle  insert_no_register(const Point_2& p,
				    Vertex_handle vnear) {
    return insert_point(p, vnear);
  }

public:
  // insert a segment
  Vertex_handle  insert(const Point_2& p0, const Point_2& p1) {
    // update input site container
    register_input_site(p0, p1);
    return insert_segment(Site_2::construct_site_2(p0, p1), Vertex_handle());
  }

  Vertex_handle  insert(const Point_2& p0, const Point_2& p1, 
			Vertex_handle vnear) {
    // update input site container
    register_input_site(p0, p1);
    return insert_segment(Site_2::construct_site_2(p0, p1), vnear);
  }

  Vertex_handle  insert(const Site_2& t) {
    // the intended use is to unify the calls to insert(...);
    // thus the site must be an exact one; 
    CGAL_precondition( t.is_exact() );

    // update input site container
    register_input_site(t);

    if ( t.is_segment() ) {
      return insert_segment(t, Vertex_handle());
    } else if ( t.is_point() ) {
      return insert_point(t.point(), Vertex_handle());
    } else {
      CGAL_precondition ( t.is_defined() );
      return Vertex_handle(); // to avoid compiler error
    }
  }

#if 0
  Vertex_handle  insert(const Site_2& t, Vertex_handle vnear)
  {
    return insert(t, vnear, true);
  }
#endif

protected:
  void register_input_site(const Point_2& p)
  {
#if STORE_INPUT_SITES
    isc_.push_back(Site_2(p));
#else
    isc_.push_back(Site_rep_2(p, p, true));
#endif
  }

  void register_input_site(const Point_2& p0, const Point_2& p1)
  {
#if STORE_INPUT_SITES
    isc_.push_back(Site_2(p0, p1));
#else
    isc_.push_back(Site_rep_2(p0, p1, false));
#endif
  }

  void register_input_site(const Site_2& t)
  {
    CGAL_precondition( t.is_exact() );
    if ( t.is_point() ) {
      register_input_site( t.point(0) );
    } else {
      register_input_site( t.point(0), t.point(1) );
    }
  }

  Vertex_handle  insert_first(const Point_2& p);
  Vertex_handle  insert_second(const Point_2& p);
  Vertex_handle  insert_third(const Point_2& p);
  Vertex_handle  insert_third(const Site_2& t, const Storage_site_2& ss);
  Vertex_handle  insert_third(Vertex_handle v0, Vertex_handle v1);

  Vertex_handle insert_point(const Point_2& p, Vertex_handle vnear);
  Vertex_handle insert_point(const Storage_site_2& ss,
			     const Site_2& t, Vertex_handle vnear);

  Triple<Vertex_handle,Vertex_handle,Vertex_handle>
  insert_point_on_segment(const Storage_site_2& ss, const Site_2& t,
			  Vertex_handle v, const Tag_true&);

  Triple<Vertex_handle,Vertex_handle,Vertex_handle>
  insert_exact_point_on_segment(const Storage_site_2& ss, const Site_2& t,
				Vertex_handle v);

  Vertex_handle insert_segment(const Site_2& t, Vertex_handle vnear);

  Vertex_handle insert_segment_interior(const Site_2& t,
					const Storage_site_2& ss,
					Vertex_handle vnear);

  template<class ITag>
  Vertex_handle  insert_intersecting_segment(const Storage_site_2& ss,
					     const Site_2& t,
					     Vertex_handle v,
					     ITag tag)
  {
    return insert_intersecting_segment_with_tag(ss, t, v, tag);
  }

  Vertex_handle
  insert_intersecting_segment_with_tag(const Storage_site_2& ss,
				       const Site_2& t,
				       Vertex_handle v, Tag_false);

  Vertex_handle
  insert_intersecting_segment_with_tag(const Storage_site_2& ss,
				       const Site_2& t,
				       Vertex_handle v, Tag_true);



public:
  // NEAREST NEIGHBOR LOCATION
  //--------------------------
  Vertex_handle  nearest_neighbor(const Point_2& p) const {
    return nearest_neighbor(Site_2::construct_site_2(p), Vertex_handle());
  }

  Vertex_handle  nearest_neighbor(const Point_2& p,
				  Vertex_handle vnear) const {
    return nearest_neighbor(Site_2::construct_site_2(p), vnear);
  }

protected:
  Vertex_handle  nearest_neighbor(const Site_2& p,
				  Vertex_handle vnear) const;


public:
  // I/O METHODS
  //------------
  template< class Stream >
  Stream& draw_dual(Stream &str) const
  {
    Finite_edges_iterator eit = finite_edges_begin();
    for (; eit != finite_edges_end(); ++eit) {
      draw_dual_edge(*eit, str);
    }
    return str;
  }

  template < class Stream > 
  Stream& draw_skeleton(Stream &str) const
  {
    Finite_edges_iterator eit = finite_edges_begin();
    for (; eit != finite_edges_end(); ++eit) {
      Site_2 p = eit->first->vertex(  cw(eit->second) )->site();
      Site_2 q = eit->first->vertex( ccw(eit->second) )->site();

      bool is_endpoint_of_seg =
	( p.is_segment() && q.is_point() &&
	  is_endpoint_of_segment(q, p) ) ||
	( p.is_point() && q.is_segment() &&
	  is_endpoint_of_segment(p, q) );

      if ( !is_endpoint_of_seg ) {
	draw_dual_edge(*eit, str);
      }
    }
    return str;
  }

  // MK: this has to be rewritten. all the checking must be done in
  // the geometric traits class.
  template< class Stream >
  Stream& draw_dual_edge(Edge e, Stream &str) const
  {
    typename Geom_traits::Line_2              l;
    typename Geom_traits::Segment_2           s;
    typename Geom_traits::Ray_2               r;
    CGAL::Parabola_segment_2<Gt>              ps;

    Object o = primal(e);

    if (CGAL::assign(l, o))   str << l;
    if (CGAL::assign(s, o))   str << s; 
    if (CGAL::assign(r, o))   str << r;
    if (CGAL::assign(ps, o))  str << ps;

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
    pc_.clear();
    isc_.clear();
  }

  void swap(const Segment_Voronoi_diagram_2& svd) {
    DG::swap(svd);
    pc_.swap(svd.pc_);
    isc_.swap(svd.isc_);
  }

  //////////////////////////////////////////////////////////////////////
  // THE METHODS BELOW ARE LOCAL
  //////////////////////////////////////////////////////////////////////

protected:
  // HELPER METHODS FOR COMBINATORIAL OPERATIONS ON THE DATA STRUCTURE
  //------------------------------------------------------------------

  // getting the symmetric edge
  Edge sym_edge(const Edge e) const {
    return sym_edge(e.first, e.second);
  }

  Edge sym_edge(const Face_handle& f, int i) const {
    Face_handle f_sym = f->neighbor(i);
    return Edge(  f_sym, f_sym->index( f->mirror_vertex(i) )  );
  }

  Edge flip(Face_handle& f, int i) {
    CGAL_precondition ( f != Face_handle() );
    CGAL_precondition (i == 0 || i == 1 || i == 2);
    CGAL_precondition( this->dimension()==2 ); 

    CGAL_precondition( f->vertex(i) != f->mirror_vertex(i) );

    this->_tds.flip(f, i);

    return Edge(f, ccw(i));
  }

  Edge flip(Edge e) {
    return flip(e.first, e.second);
  }

  bool is_degree_2(const Vertex_handle& v) const {
    Face_circulator fc = v->incident_faces();
    Face_circulator fc1 = fc;
    ++(++fc1);
    return ( fc == fc1 );
  }

  Vertex_handle insert_degree_2(Edge e) {
    return this->_tds.insert_degree_2(e.first,e.second);
  }

  Vertex_handle insert_degree_2(Edge e, const Storage_site_2& ss) {
    Vertex_handle v = insert_degree_2(e);
    v->set_site(ss);
    return v;
  }

  void remove_degree_2(Vertex_handle v) {
    CGAL_precondition( is_degree_2(v) );
    this->_tds.remove_degree_2(v);
  }

  Vertex_handle create_vertex(const Storage_site_2& ss) {
    Vertex_handle v = this->_tds.create_vertex();
    v->set_site(ss);
    return v;
  }

  Vertex_handle create_vertex_dim_up(const Storage_site_2& ss) {
    Vertex_handle v = this->_tds.insert_dim_up(infinite_vertex());
    v->set_site(ss);
    return v;
  }

protected:
  // HELPER METHODS FOR CREATING STORAGE SITES
  //------------------------------------------
  Storage_site_2 create_storage_site(const Point_2& p) {
    Point_handle ph = pc_.insert(p);
    return Storage_site_2(ph);
  }

  Storage_site_2 create_storage_site(Vertex_handle v0,
				     Vertex_handle v1) {
    return Storage_site_2( v0->storage_site().point_handle(0),
			   v1->storage_site().point_handle(0) );
  }

  Storage_site_2 split_storage_site(const Storage_site_2& ss0,
				    const Storage_site_2& ss1,
				    unsigned int i, const Tag_false&)
  {
    // Split the first storage site which is a segment using the
    // second storage site which is a exact point
    // i denotes whether the first or second half is to be created
    CGAL_precondition( ss0.is_segment() && ss1.is_point() );
    CGAL_precondition( ss1.is_exact() );
    CGAL_precondition( i < 2 );

    if ( i == 0 ) {
      return Storage_site_2(ss0.point_handle(0), ss1.point_handle(0));
    } else {
      return Storage_site_2(ss1.point_handle(0), ss0.point_handle(1));
    }
  }

  Storage_site_2 split_storage_site(const Storage_site_2& ss0,
				    const Storage_site_2& ss1,
				    unsigned int i, const Tag_true&)
  {
    // Split the first storage site which is a segment using the
    // second storage site which is a exact point
    // i denotes whether the first or second half is to be created
    CGAL_precondition( ss0.is_segment() && ss1.is_point() );
    //    CGAL_precondition( ss1.is_exact() );
    CGAL_precondition( i < 2 );

    if ( i == 0 ) {
      if ( ss0.is_exact(0) ) {
	if ( ss1.is_exact() ) {
	  return Storage_site_2(ss0.point_handle(0), ss1.point_handle(0));
	} else {
	  Storage_site_2 supp0 = ss0.supporting_segment_site();
	  Storage_site_2 supp1 = ss1.supporting_segment_site(0);

	  if ( are_parallel(supp0.site(), supp1.site()) ) {
	    supp1 = ss1.supporting_segment_site(1);
	  }
	  return Storage_site_2(supp0.point_handle(0),
				supp0.point_handle(1),
				supp1.point_handle(0),
				supp1.point_handle(1), true);
	}
      } else {
	if ( ss1.is_exact() ) {
	  return Storage_site_2(ss0.point_handle(0), ss1.point_handle(0),
				ss0.point_handle(2), ss0.point_handle(3),
				false);
	} else {
	  Storage_site_2 supp0 = ss0.supporting_segment_site();
	  Storage_site_2 supp1 = ss1.supporting_segment_site(0);

	  if ( are_parallel(supp0.site(), supp1.site()) ) {
	    supp1 = ss1.supporting_segment_site(1);
	  }
	  return Storage_site_2(supp0.point_handle(0),
				supp0.point_handle(1),
				ss0.point_handle(2),
				ss0.point_handle(3),
				supp1.point_handle(0),
				supp1.point_handle(1) );
	}
      }
    } else { // i == 1
      if ( ss0.is_exact(1) ) {
	if ( ss1.is_exact() ) {
	  return Storage_site_2(ss1.point_handle(0), ss0.point_handle(1));
	} else {
	  Storage_site_2 supp0 = ss0.supporting_segment_site();
	  Storage_site_2 supp1 = ss1.supporting_segment_site(0);

	  if ( are_parallel(supp0.site(), supp1.site()) ) {
	    supp1 = ss1.supporting_segment_site(1);
	  }
	  return Storage_site_2(supp0.point_handle(0),
				supp0.point_handle(1),
				supp1.point_handle(0),
				supp1.point_handle(1), false);
	}
      } else {
	if ( ss1.is_exact() ) {
	  return Storage_site_2(ss1.point_handle(0), ss0.point_handle(1),
				ss0.point_handle(4), ss0.point_handle(5),
				true);
	} else {
	  Storage_site_2 supp0 = ss0.supporting_segment_site();
	  Storage_site_2 supp1 = ss1.supporting_segment_site(0);

	  if ( are_parallel(supp0.site(), supp1.site()) ) {
	    supp1 = ss1.supporting_segment_site(1);
	  }
	  return Storage_site_2(supp0.point_handle(0),
				supp0.point_handle(1),
				supp1.point_handle(0),
				supp1.point_handle(1),
				ss0.point_handle(4),
				ss0.point_handle(5) );
	}
      }
    }
  }

  Storage_site_2 create_storage_site(const Storage_site_2& ss0,
				     const Storage_site_2& ss1) {
    return Storage_site_2( ss0.point_handle(0), ss0.point_handle(1),
			   ss1.point_handle(0), ss1.point_handle(1) );
  }

  Storage_site_2 create_storage_site(const Storage_site_2& ss0,
				     const Storage_site_2& ss1,
				     bool is_first_exact) {
    return Storage_site_2( ss0.point_handle(0), ss0.point_handle(1),
			   ss1.point_handle(0), ss1.point_handle(1),
			   is_first_exact);
  }

  Storage_site_2 create_storage_site_type1(const Storage_site_2& ss0,
					   const Storage_site_2& ss1,
					   const Storage_site_2& ss2) {
    return Storage_site_2( ss0.point_handle(0), ss0.point_handle(1),
			   ss1.point_handle(2), ss1.point_handle(3),
			   ss2.point_handle(0), ss2.point_handle(1));
  }

  Storage_site_2 create_storage_site_type2(const Storage_site_2& ss0,
					   const Storage_site_2& ss1,
					   const Storage_site_2& ss2) {
    return Storage_site_2( ss0.point_handle(0), ss0.point_handle(1),
			   ss1.point_handle(0), ss1.point_handle(1),
			   ss2.point_handle(4), ss2.point_handle(5));
  }

protected:
  // MK: THE FOLLOWING ARE NOT IN THE SPEC
  //======================================
  // METHODS FOR ACCESSING THE PRIMAL GRAPH
  //---------------------------------------
  // used primarily for visualization
  Point_2  primal(const Face_handle& f) const {
    return circumcenter(f);
  }

  Object   primal(const Edge e) const;

  Object   primal(const Edge_circulator& ec) const {
    return primal(*ec); 
  }

  Object   primal(const Finite_edges_iterator& ei) const {
    return primal(*ei);
  }

protected:
  void print_error_message() const;

  void print_error_message(const Tag_false&) const
  {
    static int i = 0;

    if ( i == 0 ) {
      i++;
      std::cerr << "SVD::Insert aborted: intersecting segments found"
		<< std::endl;
    }
  }

  void print_error_message(const Tag_true&) const {}

  //protected:
public:
  // wrappers for constructions
  Point_2 circumcenter(const Face_handle& f) const {
    CGAL_precondition( this->dimension()==2 || !is_infinite(f) );
    return circumcenter(f->vertex(0)->site(),
			f->vertex(1)->site(),
			f->vertex(2)->site());
  }

  Point_2 circumcenter(const Site_2& t0, const Site_2& t1, 
		       const Site_2& t2) const {
    return
    geom_traits().construct_svd_vertex_2_object()(t0, t1, t2);
  }

protected:
  // HELPER METHODS FOR INSERTION
  //-----------------------------
  void initialize_conflict_region(const Face_handle& f, List& l);

  std::pair<Face_handle,Face_handle>
  find_faces_to_split(const Vertex_handle& v, const Site_2& t) const;

  void expand_conflict_region(const Face_handle& f, const Site_2& t,
			      const Storage_site_2& ss,
			      List& l, Face_map& fm,
			      std::map<Face_handle,Sign>& sign_map,
			      Triple<bool, Vertex_handle,
			      Arrangement_type>& vcross);

  Vertex_handle add_bogus_vertex(Edge e, List& l);
  Vertex_list   add_bogus_vertices(List& l);
  void          remove_bogus_vertices(Vertex_list& vl);

  void retriangulate_conflict_region(Vertex_handle v, List& l,
				     Face_map& fm);


protected:
  // TYPES AND ACCESS METHODS FOR VISUALIZATION
  //-------------------------------------------

  // types
  typedef CGAL::Construct_svd_circle_2<Gt,Ring_tag>
  Construct_svd_circle_2;

  typedef CGAL::Construct_svd_bisector_2<Gt,Ring_tag>
  Construct_svd_bisector_2;

  typedef CGAL::Construct_svd_bisector_ray_2<Gt,Ring_tag>
  Construct_svd_bisector_ray_2;

  typedef CGAL::Construct_svd_bisector_segment_2<Gt,Ring_tag>
  Construct_svd_bisector_segment_2;

  // access
  Construct_svd_circle_2
  construct_svd_circle_2_object() const{
    return Construct_svd_circle_2();
  }

  Construct_svd_bisector_2
  construct_svd_bisector_2_object() const {
    return Construct_svd_bisector_2();
  }

  Construct_svd_bisector_ray_2
  construct_svd_bisector_ray_2_object() const {
    return Construct_svd_bisector_ray_2();
  }

  Construct_svd_bisector_segment_2
  construct_svd_bisector_segment_2_object() const { 
    return Construct_svd_bisector_segment_2(); 
  }

protected:
  // WRAPPERS FOR GEOMETRIC PREDICATES
  //----------------------------------
  bool same_points(const Site_2& p, const Site_2& q) const {
    return geom_traits().are_same_points_2_object()(p, q);
  }

  bool same_segments(const Site_2& t, Vertex_handle v) const {
    if ( is_infinite(v) ) { return false; }
    if ( t.is_point() || v->site().is_point() ) { return false; }
    return same_segments(t, v->site());
  }

  bool same_segments(const Site_2& p, const Site_2& q) const {
    CGAL_precondition( p.is_segment() && q.is_segment() );

    return
      (same_points(p.source_site(), q.source_site()) &&
       same_points(p.target_site(), q.target_site())) ||
      (same_points(p.source_site(), q.target_site()) &&
       same_points(p.target_site(), q.source_site()));
  }

  bool is_endpoint_of_segment(const Site_2& p, const Site_2& s) const
  {
    CGAL_precondition( p.is_point() && s.is_segment() );
    return ( same_points(p, s.source_site()) ||
	     same_points(p, s.target_site()) );
  }

  bool is_degenerate_segment(const Site_2& s) const {
    CGAL_precondition( s.is_segment() );
    return same_points(s.source_site(), s.target_site());
  }

  // returns:
  //   ON_POSITIVE_SIDE if q is closer to t1
  //   ON_NEGATIVE_SIDE if q is closer to t2
  //   ON_ORIENTED_BOUNDARY if q is on the bisector of t1 and t2
  Oriented_side side_of_bisector(const Site_2 &t1, const Site_2 &t2,
				 const Site_2 &q) const {
    return geom_traits().oriented_side_of_bisector_2_object()(t1, t2, q);
  }

  Sign incircle(const Site_2 &t1, const Site_2 &t2,
		const Site_2 &t3, const Site_2 &q) const {
    return geom_traits().vertex_conflict_2_object()(t1, t2, t3, q);
  }

  Sign incircle(const Site_2 &t1, const Site_2 &t2,
		const Site_2 &q) const {
    return geom_traits().vertex_conflict_2_object()(t1, t2, q);
  }

  Sign incircle(const Face_handle& f, const Site_2& q) const;

  Sign incircle(const Vertex_handle& v0, const Vertex_handle& v1,
		const Vertex_handle& v) const {
    CGAL_precondition( !is_infinite(v0) && !is_infinite(v1)
		       && !is_infinite(v) );

    return incircle( v0->site(), v1->site(), v->site());
  }

  Sign incircle(const Vertex_handle& v0, const Vertex_handle& v1,
		const Vertex_handle& v2, const Vertex_handle& v) const;


  
  bool finite_edge_interior(const Site_2& t1, const Site_2& t2,
			    const Site_2& t3, const Site_2& t4,
			    const Site_2& q,  Sign sgn) const {
    return
      geom_traits().finite_edge_interior_conflict_2_object()
      (t1,t2,t3,t4,q,sgn);
  }

  bool finite_edge_interior(const Face_handle& f, int i,
			    const Site_2& q, Sign sgn) const {
    CGAL_precondition( !is_infinite(f) &&
		       !is_infinite(f->neighbor(i)) );
    return finite_edge_interior( f->vertex( ccw(i) )->site(),
				 f->vertex(  cw(i) )->site(),
				 f->vertex(     i  )->site(),
				 f->mirror_vertex(i)->site(), q, sgn);
  }

  bool finite_edge_interior(const Vertex_handle& v1, const Vertex_handle& v2,
			    const Vertex_handle& v3, const Vertex_handle& v4,
			    const Vertex_handle& v, Sign sgn) const {
    CGAL_precondition( !is_infinite(v1) && !is_infinite(v2) &&
		       !is_infinite(v3) && !is_infinite(v4) &&
		       !is_infinite(v) );
    return finite_edge_interior( v1->site(), v2->site(),
				 v3->site(), v4->site(),
				 v->site(), sgn);
  }

  bool finite_edge_interior(const Site_2& t1, const Site_2& t2,
			    const Site_2& t3, const Site_2& q,
			    Sign sgn) const {
    return
    geom_traits().finite_edge_interior_conflict_2_object()(t1,t2,t3,q,sgn);
  }


  bool finite_edge_interior(const Site_2& t1, const Site_2& t2,
			    const Site_2& q,  Sign sgn) const {
    return
    geom_traits().finite_edge_interior_conflict_2_object()(t1,t2,q,sgn);
  }

  bool finite_edge_interior(const Face_handle& f, int i,
			    const Site_2& p, Sign sgn, int) const;

  bool finite_edge_interior(const Vertex_handle& v1, const Vertex_handle& v2,
			    const Vertex_handle& v3, const Vertex_handle& v4,
			    const Vertex_handle& v, Sign, int) const;

  bool infinite_edge_interior(const Site_2& t2, const Site_2& t3,
			      const Site_2& t4, const Site_2& q,
			      Sign sgn) const {
    return
      geom_traits().infinite_edge_interior_conflict_2_object()
      (t2,t3,t4,q,sgn);
  }


  bool infinite_edge_interior(const Face_handle& f, int i,
			      const Site_2& q, Sign sgn) const;

  bool infinite_edge_interior(const Vertex_handle& v1,
			      const Vertex_handle& v2,
			      const Vertex_handle& v3,
			      const Vertex_handle& v4,
			      const Vertex_handle& v,
			      Sign sgn) const;

  bool edge_interior(const Face_handle& f, int i,
		     const Site_2& t, Sign sgn) const;

  bool edge_interior(const Edge& e,
		     const Site_2& t, Sign sgn) const {
    return edge_interior(e.first, e.second, t, sgn);
  }

  bool edge_interior(const Vertex_handle& v1,
		     const Vertex_handle& v2,
		     const Vertex_handle& v3,
		     const Vertex_handle& v4,
		     const Vertex_handle& v,
		     Sign sgn) const;

  Arrangement_type arrangement_type(const Site_2& t, Vertex_handle v) const;
  Arrangement_type arrangement_type(const Site_2& p, const Site_2& q) const;

  bool are_parallel(const Site_2& p, const Site_2& q) const {
    return geom_traits().are_parallel_2_object()(p, q);
  }

  Oriented_side
  oriented_side(const Site_2& s1, const Site_2& s2, const Site_2& s3,
		const Site_2& supp, const Site_2& p) const {
    CGAL_precondition( supp.is_segment() && p.is_point() );
    return geom_traits().oriented_side_2_object()(s1, s2, s3, supp, p);
  }

}; // Segment_Voronoi_diagram_2


CGAL_END_NAMESPACE


#ifdef CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION
#  include <CGAL/Segment_Voronoi_diagram_2.C>
#endif



#endif // CGAL_SEGMENT_VORONOI_DIAGRAM_2_H
