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

#include <CGAL/Triangulation_2.h>
#include <CGAL/Segment_Voronoi_diagram_site_2.h>
#include <CGAL/Segment_Voronoi_diagram_data_structure_2.h>
#include <CGAL/Triangulation_face_base_2.h>
#include <CGAL/Segment_Voronoi_diagram_vertex_base_2.h>

#include <CGAL/Segment_Voronoi_diagram_constructions_C2.h>

#include <CGAL/in_place_edge_list.h>
#include <CGAL/edge_list.h>
#include <CGAL/Segment_Voronoi_diagram_traits_wrapper_2.h>

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


} // namespace CGALi




template<class Gt, class STag, class PC, class DS, class LTag >
class Segment_Voronoi_diagram_hierarchy_2;

	 //	   typename PC = Point_container<typename Gt::Point_2>,
  //	   typename PC = Point_container<typename Gt::Point_2>,

template<class Gt,
	 class PC = std::list<typename Gt::Point_2>,
	 class DS = Segment_Voronoi_diagram_data_structure_2 < 
                Segment_Voronoi_diagram_vertex_base_2<Gt,PC,
			    typename Gt::Intersections_tag>,
                Triangulation_face_base_2<Gt> >,
	 class LTag = Tag_false >
class Segment_Voronoi_diagram_2
  : private Triangulation_2<
          Segment_Voronoi_diagram_traits_wrapper_2<Gt>, DS >
{
  friend class Segment_Voronoi_diagram_hierarchy_2<Gt,Tag_true,PC,DS,LTag>;
  friend class Segment_Voronoi_diagram_hierarchy_2<Gt,Tag_false,PC,DS,LTag>;

private:
  static const char point_descriptor;
  static const char segment_descriptor;

private:
  // types and access methods needed for visualization
  //--------------------------------------------------

  typedef typename Gt::Method_tag  Method_tag;

  // types
  typedef CGAL::Construct_svd_circle_2<Gt,Method_tag>
  Construct_svd_circle_2;

  typedef CGAL::Construct_svd_bisector_2<Gt,Method_tag>
  Construct_svd_bisector_2;

  typedef CGAL::Construct_svd_bisector_ray_2<Gt,Method_tag>
  Construct_svd_bisector_ray_2;

  typedef CGAL::Construct_svd_bisector_segment_2<Gt,Method_tag>
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
  // some local types
  typedef Segment_Voronoi_diagram_traits_wrapper_2<Gt>  Modified_traits;
  typedef Triangulation_2<Modified_traits,DS>           DG;
  typedef DG                         Delaunay_graph;
  typedef typename DG::Vertex        Vertex;
  typedef typename DG::Face          Face;

  typedef LTag            Use_in_place_edge_list_tag;

public:
  // TYPES
  //------
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

  typedef Simple_container_wrapper<PC>           Point_container;
  typedef typename Point_container::iterator     Point_handle;

  typedef typename DG::size_type                 size_type;

protected:
  // some more local types
  typedef typename Gt::Intersections_tag       Intersections_tag;

  typedef typename DS::Vertex_base             Vertex_base;

  typedef std::map<Face_handle,bool>           Face_map;
  typedef std::vector<Edge>                    Edge_vector;

  typedef std::list<Vertex_handle>         Vertex_list;
  typedef typename Vertex_list::iterator   Vertex_list_iterator;
  typedef Vertex_handle                    Vh_triple[3];

  typedef
  typename Data_structure::Vertex_base::Storage_site_2  Storage_site_2;

  // the in place edge list
  typedef typename
  CGALi::SVD_which_list<Edge,Use_in_place_edge_list_tag>::List  List;

  typedef enum { NO_CONFLICT = -1, INTERIOR, LEFT_VERTEX,
		 RIGHT_VERTEX, BOTH_VERTICES, ENTIRE_EDGE }
  Conflict_type;

  static Conflict_type opposite(const Conflict_type& ct) {
    if ( ct == RIGHT_VERTEX ) { return LEFT_VERTEX; }
    if ( ct == LEFT_VERTEX ) { return RIGHT_VERTEX; }
    return ct;
  }

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
  const Geom_traits& geom_traits() const {
    return DG::geom_traits();
  }

  size_type number_of_vertices() const {
    return DG::number_of_vertices();
  }

  size_type number_of_incident_segments(Vertex_handle v) const;


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

  //  Point_iterator points_begin() const;
  //  Point_iterator points_end() const;

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
  template<class Input_iterator>
  size_type insert(Input_iterator first, Input_iterator beyond) {
    return insert_with_tag(first, beyond, Tag_false());
  }

protected:
  template<class Input_iterator>
  size_type insert_with_tag(Input_iterator first,
			    Input_iterator beyond,
			    Tag_true)
  {
    // MK::ERROR: this changes the data I would have to copy them to
    //            a vector first and then do the suffling thing...
    std::random_shuffle(first, beyond);
    return insert_with_tag(first, beyond, Tag_false());
  }

  template<class Input_iterator>
  size_type insert_with_tag(Input_iterator first,
			    Input_iterator beyond,
			    Tag_false)
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

public:
  template<class Input_iterator, class True_false_tag>
  size_type insert(Input_iterator first, Input_iterator beyond,
		   True_false_tag tag)
  {
    return insert_with_tag(first, beyond, tag);
  }

  // insert a point
  Vertex_handle  insert(const Point_2& p) {
    return insert_point(p, Vertex_handle());
  }

  Vertex_handle  insert(const Point_2& p, Vertex_handle vnear) {
    return insert_point(p, vnear);
  }

  // insert a segment
  Vertex_handle  insert(const Point_2& p1, const Point_2& p2) {
    return
    insert_segment(Site_2(p1, p2), Vertex_handle(), true);
  }

  Vertex_handle  insert(const Point_2& p0, const Point_2& p1, 
			Vertex_handle vnear) {
    return
    insert_segment(Site_2(p0, p1), vnear, true);
  }

  // MK::ERROR: I may not want to expose this...
  // insert a site
  Vertex_handle  insert(const Site_2& t) {
    if ( t.is_segment() ) {
      return insert_segment(t, Vertex_handle(), true);
    } else if ( t.is_point() ) {
      // MK::ERROR: the following does not work if the point is not
      //            exact...
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

public:
  // REMOVAL
  //--------

  // returns the number of sites removed
  // possible answers:
  // 0 : no site was removed; this can only happen if we ask to
  //     remove a point that is the endpoint of a segment
  // 1 : a single site was removed; this can happen if
  //     (1) we ask to remove a point which is not the endpoint of
  //         a segment
  //     (2) we ask to remove an open segment, i.e., we keep its
  //         endpoints
  //     (3) we ask to remove a closed segment, but its two
  //         endpoints are also endpoints of other segments and
  //         thus cannot be removed
  // 2 : two sites were removed; this can happen if we ask to
  //     remove a closed segment, but only one of its endpoints is
  //     removed; the other is also the endpoint of another
  //     segment
  // 2 : three sites where removed; this can happen when we ask to
  //     remove a closed segment; in this case the two endpoints
  //     are not endpoints of other segment, and thus they are
  //     removed as well
  unsigned int remove(Vertex_handle v,
		      bool remove_endpoints = true);


public:
  // NEAREST NEIGHBOR LOCATION
  //--------------------------
  Vertex_handle  nearest_neighbor(const Point_2& p) const {
    return nearest_neighbor(Site_2(p), Vertex_handle());
  }

  Vertex_handle  nearest_neighbor(const Point_2& p,
				  Vertex_handle vnear) const {
    return nearest_neighbor(Site_2(p), vnear);
  }

protected:
  Vertex_handle  nearest_neighbor(const Site_2& p,
				  Vertex_handle vnear) const;


#if 0
  // MK: THE FOLLOWING ARE NOT IN THE SPEC
  //======================================
public:
  // ACCESS TO THE DUAL
  //-------------------
  Point_2 dual(const Face_handle& f) const;
  Object  dual(const Edge e) const;

  Object  dual(const Edge_circulator& ec) const {
    return dual(*ec);
  }

  Object  dual(const Finite_edges_iterator& ei) const {
    return dual(*ei);
  }
#endif

public:
  // I/O
  //----
protected:
  template < class Stream >
  Stream& write_sites(Stream& str) const
  {
    str << number_of_vertices() << std::endl;

    Finite_vertices_iterator vit = finite_vertices_begin();
    for (; vit != finite_vertices_end(); ++vit) {
      if ( vit->is_point() ) {
	str << "p " << vit->storage_site().point() << std::endl;
      } else {
	str << "s " << vit->storage_site().segment() << std::endl;
      }
    }
    return str;
  }

  template < class Stream >
  Stream& draw_sites(Stream &str) const
  {
    Finite_vertices_iterator vit = finite_vertices_begin();
    for (; vit != finite_vertices_end(); ++vit) {
      if ( vit->is_point() ) {
	str << vit->point();
      } else {
	str << vit->segment();
      }
    }
    return str;
  }

public:
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

protected:
#if 0
  template < class Stream >
  Stream& draw_Voronoi_circles(Stream& str) const
  {
    Finite_faces_iterator fit = finite_faces_begin();
    for (; fit != finite_faces_end(); ++fit) {
      typename Gt::Circle_2 c = circumcircle(Face_handle(fit));
      str << c;
      str << c.center();
    }
    return str;
  }
#endif

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
  }

  void swap(const Segment_Voronoi_diagram_2& svd) {
    DG::swap(svd);
    pc_.swap(svd.pc_);
  }

  const Data_structure&  ds() const { return this->_tds; }
  const Point_container& point_container() const { return pc_; }


protected:
  // MK: THE FOLLOWING ARE NOT IN THE SPEC
  //======================================
  // but they are needed internally
  // Primal
  Point_2  primal(const Face_handle& f) const;
  Object   primal(const Edge e) const;
  Object   primal(const Edge_circulator& ec) const {
    return primal(*ec);
  }
  Object   primal(const Finite_edges_iterator& ei) const {
    return primal(*ei);
  }

protected:
  // wrappers for the geometric predicates

  bool are_same_points(const Site_2& p, const Site_2& q) const;
  bool same_segments(const Site_2& t, Vertex_handle v) const;
  bool same_segments(const Site_2& p, const Site_2& q) const;

  bool is_endpoint_of_segment(const Site_2& p, const Site_2& s) const
  {
    CGAL_precondition( p.is_point() && s.is_segment() );
    return ( are_same_points(p, s.source_site()) ||
	     are_same_points(p, s.target_site()) );
  }

  bool is_degenerate_segment(const Site_2& s) const {
    CGAL_precondition( s.is_segment() );
    return are_same_points(s.source_site(), s.target_site());
  }

  // returns:
  //   ON_POSITIVE_SIDE if q is closer to t1
  //   ON_NEGATIVE_SIDE if q is closer to t2
  //   ON_ORIENTED_BOUNDARY if q is on the bisector of t1 and t2
  Oriented_side side_of_bisector(const Site_2 &t1,
				 const Site_2 &t2,
				 const Site_2 &q) const;

  Sign incircle(const Site_2 &t1, const Site_2 &t2,
		const Site_2 &t3, const Site_2 &q) const;

  Sign incircle(const Site_2 &t1, const Site_2 &t2,
		const Site_2 &q) const;


  Sign incircle(const Face_handle& f, const Site_2& q) const;


  Sign incircle(const Vertex_handle& v0, const Vertex_handle& v1,
		const Vertex_handle& v) const;

  Sign incircle(const Vertex_handle& v0, const Vertex_handle& v1,
		const Vertex_handle& v2, const Vertex_handle& v) const;


  
  bool finite_edge_interior(const Site_2& t1, const Site_2& t2,
			    const Site_2& t3, const Site_2& t4,
			    const Site_2& q,  Sign sgn) const;

  bool finite_edge_interior(const Face_handle& f, int i,
			    const Site_2& q, Sign sgn) const;

  bool finite_edge_interior(const Vertex_handle& v1,
			    const Vertex_handle& v2,
			    const Vertex_handle& v3,
			    const Vertex_handle& v4,
			    const Vertex_handle& v,
			    Sign sgn) const;

  bool finite_edge_interior_degenerated(const Site_2& t1, const Site_2& t2,
					const Site_2& t3, const Site_2& q,
					Sign sgn) const;


  bool finite_edge_interior_degenerated(const Site_2& t1, const Site_2& t2,
					const Site_2& q,  Sign sgn) const;

  bool finite_edge_interior_degenerated(const Face_handle& f, int i,
					const Site_2& p, Sign sgn) const;

  bool finite_edge_interior_degenerated(const Vertex_handle& v1,
					const Vertex_handle& v2,
					const Vertex_handle& v3,
					const Vertex_handle& v4,
					const Vertex_handle& v,
					Sign Sign) const;

  bool infinite_edge_interior(const Site_2& t2, const Site_2& t3,
			      const Site_2& t4, const Site_2& q,
			      Sign sgn) const;


  bool infinite_edge_interior(const Face_handle& f, int i,
			      const Site_2& q, Sign sgn) const;

  bool infinite_edge_interior(const Vertex_handle& v1,
			      const Vertex_handle& v2,
			      const Vertex_handle& v3,
			      const Vertex_handle& v4,
			      const Vertex_handle& v,
			      Sign sgn) const;

  Conflict_type
  finite_edge_conflict_type_degenerated(const Site_2& t1,
					const Site_2& t2,
					const Site_2& t) const;

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

#if 0
  bool is_degenerate_edge(const Site_2& t1,
			  const Site_2& t2,
			  const Site_2& t3,
			  const Site_2& t4) const {
    return geom_traits().is_degenerate_edge_2_object()
      (t1, t2, t3, t4);
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
    Vertex_handle v4 = f->mirror_vertex(i);

    return is_degenerate_edge(v1, v2, v3, v4);
  }

  bool is_degenerate_edge(const Edge& e) const {
    return is_degenerate_edge(e.first, e.second);
  }
#endif

  bool do_intersect(const Site_2& t, Vertex_handle v) const;
  bool do_intersect(const Site_2& p, const Site_2& q) const
  {
    std::pair<int,int> res =
      geom_traits().do_intersect_2_object()(p, q);

#if 1
    CGAL_assertion( res.first <= 4 && res.second <= 4 );

    if ( res.first == 2 ) {
      CGAL_assertion( res.second == 2 );
    } else if ( res.second == 2 ) {
      CGAL_assertion( res.first == 2 );
    }

    if ( res.first == 3 ) {
      CGAL_assertion( res.second == 3 );
    } else if ( res.second == 3 ) {
      CGAL_assertion( res.first == 3 );
    }
#endif
    if ( res.first < 2 && res.second < 2 ) { return false; }

    return (res.first != 3);
  }

  bool are_parallel(const Site_2& p, const Site_2& q) const
  {
    return geom_traits().are_parallel_2_object()(p, q);
  }

  Oriented_side
  oriented_side(const Site_2& s1, const Site_2& s2, const Site_2& s3,
		const Site_2& supp, const Site_2& p) const
  {
    CGAL_precondition( supp.is_segment() && p.is_point() );
    return geom_traits().oriented_side_2_object()(s1, s2, s3, supp, p);
  }

  void print_error_message() const
  {
    std::cerr << "SVD::Insert aborted: intersecting segments found"
	      << std::endl;
  }

  //protected:
public:
  // wrappers for constructions
  Point_2 circumcenter(const Face_handle& f) const;
  Point_2 circumcenter(const Site_2& t0, 
		       const Site_2& t1, 
		       const Site_2& t2) const;

#if 0
  typename Gt::Circle_2 circumcircle(const Face_handle& f) const;
  typename Gt::Circle_2 circumcircle(const Site_2& t0, const Site_2& t1, 
				     const Site_2& t2) const;

  typename Gt::Line_2 circumcircle(const Point_2& p0, const Point_2& p1) const;
#endif

protected:
  // wrappers for combinatorial operations on the data structure

  // getting the symmetric edge
  Edge sym_edge(const Edge e) const {
    return sym_edge(e.first, e.second);
  }

  Edge sym_edge(const Face_handle& f, int i) const {
    Face_handle f_sym = f->neighbor(i);
    return Edge(  f_sym, f_sym->index( f->mirror_vertex(i) )  );
  }

  Edge flip(Face_handle& f, int i);
  Edge flip(Edge e);

  //  Vertex_handle insert_in_face(Face_handle& f, const Weighted_point& p);

  bool          is_degree_2(const Vertex_handle& v) const;

  Vertex_handle insert_degree_2(Edge e);
  Vertex_handle insert_degree_2(Edge e, const Storage_site_2& ss);

  void          remove_degree_2(Vertex_handle v);
#if 0
  void          remove_degree_3(Vertex_handle v);
  void          remove_degree_3(Vertex_handle v, Face* f);
#endif

  // this was defined because the hierarchy needs it
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
  // insertion of the first three sites

  // the first two objects can only be points, since we always
  // add the endpoints first and then the segment.
  Storage_site_2 create_storage_site(const Point_2& p)
  {
    Point_handle ph = pc_.insert(p);
    return Storage_site_2(ph);
  }

  Storage_site_2 create_storage_site(Vertex_handle v0,
				     Vertex_handle v1)
  {
    //    typedef typename Storage_site_2::Handle_pair   Point_handle_pair;

    //    Point_handle_pair ph_pair(v0->storage_site().point_handle(),
    //			      v1->storage_site().point_handle());
    //    return Storage_site_2( ph_pair );
    return Storage_site_2( v0->storage_site().point_handle(0),
			   v1->storage_site().point_handle(0) );
  }

  Vertex_handle  insert_first(const Point_2& p);
  Vertex_handle  insert_second(const Point_2& p);
  Vertex_handle  insert_third(const Point_2& p);
  //  Vertex_handle  insert_third(const Point_2& p0, const Point_2& p1);
  Vertex_handle  insert_third(Vertex_handle v0, Vertex_handle v1);

  Vertex_handle insert_point(const Point_2& p, Vertex_handle vnear);
  Vertex_handle insert_point(const Storage_site_2& ss,
			     const Site_2& t, Vertex_handle vnear);

  Vertex_handle insert_segment(const Site_2& t, Vertex_handle vnear,
			       bool insert_endpoints);

  Vertex_handle insert_segment_interior(const Site_2& t,
					const Storage_site_2& ss,
					Vertex_handle vnear,
					bool insert_endpoints);

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



  // methods for insertion
  void initialize_conflict_region(const Face_handle& f, List& l);

  void expand_conflict_region(const Face_handle& f, const Site_2& t,
			      const Storage_site_2& ss,
			      List& l, Face_map& fm,
			      std::map<Face_handle,Sign>& sign_map,
			      std::pair<bool, Vertex_handle>& vcross,
			      std::vector<Vh_triple*>* fe);

  Vertex_handle add_bogus_vertex(Edge e, List& l);
  Vertex_list   add_bogus_vertices(List& l);
  void          remove_bogus_vertices(Vertex_list& vl);

  // MK: this is not currently used
  std::vector<Face*> get_faces_for_recycling(Face_map& fm,
					     unsigned int n_wanted);

  void retriangulate_conflict_region(Vertex_handle v, List& l,
				     Face_map& fm);

protected:
  // methods for removal
#if 0
  std::pair<Vertex_handle,Vertex_handle >
  endpoint_vertices(Vertex_handle v) const;

  bool is_endpoint_of_segment(Vertex_handle v) const;

  void  remove_first(Vertex_handle v);
  void  remove_second(Vertex_handle v);
  unsigned int remove_third(Vertex_handle v, bool remove_endpoints);
  unsigned int remove_degree_2(Vertex_handle v,
			       bool remove_endpoints);
  unsigned int remove_degree_3(Vertex_handle v,
			       bool remove_endpoints);
  unsigned int remove_degree_d(Vertex_handle v,
			       bool remove_endpoints);

  void  minimize_degree(Vertex_handle v);

  void find_conflict_region_remove(const Vertex_handle& v,
				   const Vertex_handle& vnearest,
				   List& l, Face_map& fm,
				   std::vector<Vh_triple*>* fe);
#endif
protected:
  // methods for I/O

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

private:
  Point_container pc_;

}; // Segment_Voronoi_diagram_2


CGAL_END_NAMESPACE


#ifdef CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION
#  include <CGAL/Segment_Voronoi_diagram_2.C>
#endif



#endif // CGAL_SEGMENT_VORONOI_DIAGRAM_2_H
