// Copyright (c) 1999-2016   INRIA Nancy - Grand Est (France).
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
// Author(s)     : Iordan Iordanov  <Iordan.Iordanov@loria.fr>
//                 


#ifndef CGAL_PERIODIC_4_HYPERBOLIC_TRIANGULATION_2_H
#define CGAL_PERIODIC_4_HYPERBOLIC_TRIANGULATION_2_H

#include <CGAL/basic.h>

#include <iostream>
#include <algorithm>
#include <cmath>
#include <functional>
#include <list>

#include <boost/tuple/tuple.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_smallint.hpp>
#include <boost/random/variate_generator.hpp>
#include <CGAL/Triangulation_utils_2.h>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/utility.h>
#include <CGAL/use.h>

#include <CGAL/Triangulation_2.h>
#include <CGAL/Periodic_4_hyperbolic_triangulation_ds_face_base_2.h>
#include <CGAL/Periodic_4_hyperbolic_triangulation_ds_vertex_base_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Periodic_4_hyperbolic_triangulation_iterators_2.h>
#include <CGAL/Unique_hash_map.h>
#include <CGAL/internal/Exact_type_selector.h>
#include <CGAL/NT_converter.h>
#include <CGAL/Triangulation_structural_filtering_traits.h>
#include <CGAL/Hyperbolic_octagon_word_4.h>
#include <CGAl/Hyperbolic_translation_info.h>



namespace CGAL {

#ifndef CGAL_P4T2_STRUCTURAL_FILTERING_MAX_VISITED_CELLS
	#define CGAL_P4T2_STRUCTURAL_FILTERING_MAX_VISITED_CELLS 2500
#endif // no CGAL_P4T2_STRUCTURAL_FILTERING_MAX_VISITED_CELLS


#ifndef CGAL_NO_STRUCTURAL_FILTERING
namespace internal {
// structural filtering is performed only for EPIC
struct Periodic_4_hyperbolic_structural_filtering_2_tag {};
struct No_periodic_4_hyperbolic_structural_filtering_2_tag {};

template <bool filter>
struct Periodic_4_hyperbolic_structural_filtering_selector_2 {
#ifdef FORCE_STRUCTURAL_FILTERING
  typedef Periodic_4_hyperbolic_structural_filtering_2_tag  Tag;
#else
  typedef No_periodic_4_hyperbolic_structural_filtering_2_tag  Tag;
#endif
};

template <>
struct Periodic_4_hyperbolic_structural_filtering_selector_2<true> {
  typedef Periodic_4_hyperbolic_structural_filtering_2_tag  Tag;
};
}
#endif // no CGAL_NO_STRUCTURAL_FILTERING




template < 	class GT,
			class TDS = Triangulation_data_structure_2<
				Periodic_4_hyperbolic_triangulation_ds_vertex_base_2<GT>,
				Periodic_4_hyperbolic_triangulation_ds_face_base_2<GT>
			>
		>
class Periodic_4_hyperbolic_triangulation_2 : public Triangulation_2<GT, TDS> {

	typedef Periodic_4_hyperbolic_triangulation_2<GT, TDS> 		Self;
  	typedef Triangulation_2<GT, TDS>                          	Base;

private:
	typedef typename GT::FT                      				FT;

public:

#ifndef CGAL_CFG_USING_BASE_MEMBER_BUG_2  
  using Base::cw;
  using Base::ccw;
  using Base::geom_traits;
#endif

	typedef GT 										        Geometric_traits;
	typedef TDS 									        Triangulation_data_structure;
	typedef unsigned short int 								Int;
	typedef Hyperbolic_octagon_word_4<Int, GT>				Offset;
	typedef typename GT::Circle_2 							Circle_2;
	typedef Circle_2 								      	Circle;
	typedef typename GT::Point_2 			    			Point_2;
	typedef Point_2 								      	Point;
	typedef typename GT::Segment_2 							Segment_2;
	typedef Segment_2 								    	Segment;
	typedef typename GT::Triangle_2 						Triangle_2;
	typedef Triangle_2 								    	Triangle;

	typedef std::pair<Point, Offset> 				      	Periodic_point;
	typedef array< std::pair<Point, Offset>, 2 >			Periodic_segment;
	typedef array< std::pair<Point, Offset>, 3 > 			Periodic_triangle;	

	typedef typename TDS::Vertex 					        Vertex;
	typedef typename TDS::Edge 						        Edge;
	typedef typename TDS::Face 						        Face;

	typedef typename TDS::Vertex_handle 			        Vertex_handle;
	typedef typename TDS::Face_handle 				        Face_handle;

	typedef typename TDS::size_type             	    	size_type;
	typedef typename TDS::difference_type        	  		difference_type;

	typedef typename TDS::Face_iterator         	  		Face_iterator;
	typedef typename TDS::Edge_iterator          	  		Edge_iterator;
	typedef typename TDS::Vertex_iterator        	  		Vertex_iterator;
	typedef typename TDS::Face_circulator       	  		Face_circulator;
  	typedef typename TDS::Edge_circulator         			Edge_circulator;
  	typedef typename TDS::Vertex_circulator       			Vertex_circulator;

  	typedef typename Base::Line_face_circulator     		Line_face_circulator;

	typedef Periodic_4_hyperbolic_triangulation_triangle_iterator_2<Self>
															Periodic_triangle_iterator;
	typedef Periodic_4_hyperbolic_triangulation_segment_iterator_2<Self>
															Periodic_segment_iterator;
	typedef Periodic_4_hyperbolic_triangulation_point_iterator_2<Self>
															Periodic_point_iterator;

	typedef Point                                			value_type;
	typedef const value_type&                    			const_reference;

	typedef Tag_false 								        Weighted_tag;


	enum Iterator_type {
		STORED=0,
		UNIQUE, 				// 1
		STORED_COVER_DOMAIN, 	// 2
		UNIQUE_COVER_DOMAIN 	// 3
	};

	enum Locate_type {
		VERTEX=0, 
		EDGE, 					// 1
		FACE,	   				// 2
		CELL, 					// 3
		EMPTY , 				// 4
		OUTSIDE_CONVEX_HULL, 	// unused, for compatibility with Alpha_shape_3
		OUTSIDE_AFFINE_HULL 	// unused, for compatibility with Alpha_shape_3
	}; 

protected:
	GT       _gt;
	TDS 	 _tds; 
	Circle 	 _domain;

public:
  	Periodic_4_hyperbolic_triangulation_2(Geometric_traits gt) : 
		_gt(gt), _domain( Circle_2(Point(0,0), 1*1) ), _tds() 
	{
		init_tds();
	}  

	Periodic_4_hyperbolic_triangulation_2(
		const Circle_2 domain = Circle_2(Point_2(FT(0),FT(0)), FT(1*1)), 
		const Geometric_traits &gt = Geometric_traits() ) :
		_gt(gt), _domain(domain), _tds()
	{
		init_tds();
	}

	Periodic_4_hyperbolic_triangulation_2(const Periodic_4_hyperbolic_triangulation_2& tr) :
		_gt(tr.geom_traits()),
		_domain(tr._domain)
	{
		CGAL_triangulation_expensive_postcondition(*this == tr);
	}

	Periodic_4_hyperbolic_triangulation_2& operator=(Periodic_4_hyperbolic_triangulation_2 tr) {
		swap(tr);
		return *this;
	}

	void swap(Periodic_4_hyperbolic_triangulation_2& tr) {
		_tds.swap(tr._tds);
		std::swap(tr._gt, _gt);
		std::swap(_domain, tr._domain);
	}

	void clear() {
		_tds.clear();
		init_tds();
	}


  template< class FaceIt >
  Vertex_handle insert_in_hole(const Point& p, FaceIt face_begin, FaceIt face_end) 
  {
	Vertex_handle v = _tds.insert_in_hole(face_begin, face_end);
	v->set_point(p);
	return v;
  }


  Vertex_handle insert_in_face(const Point& p, Face_handle fh) {
	Vertex_handle vh = _tds.insert_in_face(fh);
	vh->set_point(p);
	return vh;
  }


private:
	void init_tds() {
		_tds.set_dimension(-2);
	}

public:
	const Geometric_traits& geom_traits() const {
		return _gt;
	}


	/// Returns the data structure storing the triangulation.
  	const TDS & tds() const {
		return _tds;
  	}
  	
  	/// Returns the data structure storing the triangulation.
  	TDS & tds() {
		return _tds;
  	}

	const Circle_2& domain() const {
		return _domain;
	}

	void set_domain(const Circle_2 domain) {
		clear();
		_domain = domain;
		_gt.set_domain(domain);
	}


	size_type number_of_faces() const {
		return _tds.number_of_faces();
	}
  
	size_type number_of_edges() const {
		return _tds.number_of_edges();
	}
  
	size_type number_of_vertices() const {
		return _tds.number_of_vertices();
	}

	void set_offsets(Face_handle f, Offset o0, Offset o1, Offset o2) {
		f->set_offsets(o0, o1, o2);
	}

	Comparison_result compare_x(const Point &p1, const Point &p2) const {
		return geom_traits().compare_x_2_object()(p1, p2);
	}

	Comparison_result compare_x(const Point  &p1, const Point  &p2,
								const Offset &o1, const Offset &o2) const {
		return geom_traits().compare_x_2_object()(p1, p2, o1, o2);
	}


	Comparison_result compare_y(const Point &p1, const Point &p2) const {
		return geom_traits().compare_y_2_object()(p1, p2);
	}

	Comparison_result compare_y(const Point  &p1, const Point  &p2,
								const Offset &o1, const Offset &o2) const {
		return geom_traits().compare_y_2_object()(p1, p2, o1, o2);
	}


	Comparison_result compare_xy(const Point &p1, const Point &p2) const {
		return geom_traits().compare_xy_2_object()(p1, p2);
	}

	Comparison_result compare_xy(const Point  &p1, const Point  &p2,
								 const Offset &o1, const Offset &o2) const {
		return geom_traits().compare_xy_2_object()(p1, p2, o1, o2);
	}


	Orientation orientation(const Point &p1, const Point &p2, const Point &p3) const{
		return geom_traits().orientation_2_object()(p1, p2, p3);
	}

	Orientation orientation(const Point  &p1, const Point  &p2, const Point  &p3,
							const Offset &o1, const Offset &o2, const Offset &o3) const {
		return geom_traits().orientation_2_object()(p1, p2, p3, o1, o2, o3);
	}


	bool equal(const Point &p1, const Point &p2) const {
		return ( compare_xy(p1, p2) == EQUAL );
	}

	bool equal(const Point  &p1, const Point  &p2,
			   const Offset &o1, const Offset &o2) const {
		return ( compare_xy(p1, p2, o1, o2) == EQUAL );
	}


	Oriented_side side_of_oriented_circle(const Point  &p0, const Point  &p1, const Point  &p2, const Point &q) const {
		return geom_traits().side_of_oriented_circle_2_object()(p0, p1, p2, q);
	}

	Oriented_side side_of_oriented_circle(const Point  &p0, const Point  &p1, const Point  &p2, const Point &q,
										  const Offset &o0, const Offset &o1, const Offset &o2, const Offset &oq) const {
		return geom_traits().side_of_oriented_circle_2_object()(p0, p1, p2, q, o0, o1, o2, oq);
	}


	Periodic_triangle construct_periodic_3_triangle(const Point &p1, const Point &p2, const Point &p3) const {
		return make_array(	 std::make_pair(p1,Offset()),
										 std::make_pair(p2,Offset()), 
										 std::make_pair(p3,Offset()) );
	}

	Periodic_triangle construct_periodic_3_triangle(const Point  &p1, const Point  &p2, const Point  &p3,
													const Offset &o1, const Offset &o2, const Offset &o3) const {
		return make_array(	 std::make_pair(p1,o1), 
									   std::make_pair(p2,o2),
										 std::make_pair(p3,o3) );
	}

	Periodic_segment construct_periodic_3_segment(const Point &p1, const Point &p2) const {
		return make_array(	std::make_pair(p1,Offset()), 
									  std::make_pair(p2,Offset()) );
	}

	Periodic_segment construct_periodic_3_segment(	const Point  &p1, const Point  &p2,
													const Offset &o1, const Offset &o2) const {
		return make_array(	std::make_pair(p1,o1), 
									  std::make_pair(p2,o2) );
	}


	Triangle construct_triangle(const Point &p1, const Point &p2, const Point &p3) const {
		return geom_traits().construct_triangle_2_object()(p1,p2,p3);
	}

	Triangle construct_triangle(const Point  &p1, const Point  &p2, const Point  &p3,
								const Offset &o1, const Offset &o2, const Offset &o3) const {
		return geom_traits().construct_triangle_2_object()(p1,p2,p3,o1,o2,o3);
	}

	Triangle construct_triangle(const Periodic_triangle& tri) {
		return construct_triangle(	tri[0].first,  tri[1].first,  tri[2].first,
									tri[0].second, tri[1].second, tri[2].second);
	}

	Segment construct_segment(const Point &p1, const Point &p2) const {
		return geom_traits().construct_segment_2_object()(p1, p2);
	}

	Segment construct_segment(	const Point &p1, const Point &p2,
								const Offset &o1, const Offset &o2) const {
	  //std::cout << ".....making a segment between a point with offset " << o1 << " and another with offset " << o2 << std::endl; 
		return geom_traits().construct_segment_2_object()(p1,p2,o1,o2);
	}

	Segment construct_segment(const Periodic_segment& seg) const {
		return construct_segment(seg[0].first,  seg[1].first,
								 seg[0].second, seg[1].second);
	}

	Point construct_point(const Point& p, const Offset &o) const {
		return geom_traits().construct_point_2_object()(p,o);
	}

	Point construct_point(const Periodic_point& pp) const {
		return construct_point(pp.first, pp.second);
	}

	Periodic_point periodic_point( const Face_handle c, int i) const {
		return std::make_pair(c->vertex(i)->point(), c->offset(i));
	}


	Periodic_segment periodic_segment(const Face_handle c, int i, int j) const {
		CGAL_triangulation_precondition( i != j );
		CGAL_triangulation_precondition( number_of_vertices() != 0 );
		CGAL_triangulation_precondition( (i >= 0 && i <= 2) && (j >= 0 && j <= 2) );
		return make_array( 	std::make_pair(c->vertex(i)->point(), c->offset(i)),
							std::make_pair(c->vertex(j)->point(), c->offset(j)) );
	}
  

	Periodic_segment periodic_segment(const Edge & e) const {
		return periodic_segment(e.first, e.second, e.third);
	}


	Periodic_segment periodic_segment(const Edge & e, const Offset& o) const {
		Face_handle f = e.first;
		int i = e.second;
		return make_array( 	std::make_pair(f->offset(cw(i)).apply(f->vertex(cw(i))->point()),  o),
							std::make_pair(f->offset(ccw(i)).apply(f->vertex(ccw(i))->point()), o) );
	}


	Periodic_triangle periodic_triangle(const Face_handle c, int i) const;

	Periodic_triangle periodic_triangle(const Face & f) const {
		return periodic_triangle(f.first, f.second);
	}

	Point point(const Periodic_point & pp) const {
		return construct_point(	pp.first, 
								pp.second);
	}


	Segment segment(const Face_handle & fh, int idx) const {
	  return construct_segment( fh->vertex(idx)->point(),  fh->vertex(ccw(idx))->point(),
								fh->offset(idx),           fh->offset(ccw(idx)) );
	}

	Segment segment(const Point& src, const Point& tgt) {
	  return construct_segment(src, tgt);
	}

	Segment segment(const pair<Face_handle, int> & edge) {
	  return segment(edge.first, ccw(edge.second));
	}

	Segment segment(const Periodic_segment & ps) const {
		return construct_segment(	ps[0].first,  ps[1].first,
									ps[0].second, ps[1].second );
	}


	Triangle triangle(const Face_handle & fh) const {
	  return construct_triangle( fh->vertex(0)->point(), fh->vertex(1)->point(), fh->vertex(2)->point(),
								 fh->offset(0),          fh->offset(1),          fh->offset(2)        );
	}


	Triangle triangle(const Periodic_triangle & pt) const {
		return construct_triangle(	pt[0].first, pt[1].first, pt[2].first,
												pt[0].second,pt[1].second,pt[2].second );
	}


	bool is_vertex(const Point & p, Vertex_handle & v) const;

	bool is_vertex(Vertex_handle v) const {
		return _tds.is_vertex(v);
	}

	bool is_edge(Vertex_handle u, Vertex_handle v,
				 Face_handle & c, int & i, int & j) const {
		return _tds.is_edge(u, v, c, i, j);
	}

	bool is_edge(	Vertex_handle u, const Offset & off_u,
					Vertex_handle v, const Offset & off_v,
					Face_handle & c, int & i, int & j) const {
		if (!_tds.is_edge(u,v,c,i,j)) return false;
		if ((c->offset(i) == off_u) && (c->offset(j) == off_v))
			return true;
		
		// it might be that different cells containing (u,v) yield
		// different offsets, which forces us to test for all possibilities.
		else {
			Face_circulator ccirc = incident_cells(c,i,j,c);
			while (++ccirc != c) {
				i = ccirc->index(u);
				j = ccirc->index(v);
				if ((ccirc->offset(i) == off_u) && (ccirc->offset(j) == off_v)) {
					c = ccirc;
					return true;
				}
			}
			return false;
		}
	}

	bool is_face(	Vertex_handle u, Vertex_handle v, Vertex_handle w,
					Face_handle & c, int & i, int & j, int & k) const {
		return _tds.is_face(u, v, w, c, i, j, k);
	}

	bool is_face(Vertex_handle u, const Offset & off_u,
				 Vertex_handle v, const Offset & off_v,
				 Vertex_handle w, const Offset & off_w,
				 Face_handle & c, int & i, int & j, int & k) const {
		if (!_tds.is_face(u,v,w,c,i,j,k)) return false;
		if (	(get_offset(c,i) == off_u)
			 && (get_offset(c,j) == off_v)
			 && (get_offset(c,k) == off_w) )
			return true;
	
		// it might be that c and c->neighbor(l) yield different offsets
		// which forces us to test for both possibilities.
		else {
			int l = 6-i-j-k;
			c = c->neighbor(l);
			i = c->index(u);
			j = c->index(v);
			k = c->index(w);      
			return (	(c->offset(i) == off_u)
					 && (c->offset(j) == off_v)
					 && (c->offset(k) == off_w) );
		}
	}

	bool has_vertex(const Face & f, Vertex_handle v, int & j) const {
		return _tds.has_vertex(f.first, f.second, v, j);
	}

	bool has_vertex(Face_handle c, int i, Vertex_handle v, int & j) const {
		return _tds.has_vertex(c, i, v, j);
	}

	bool has_vertex(const Face & f, Vertex_handle v) const {
		return _tds.has_vertex(f.first, f.second, v);
	}

	bool has_vertex(Face_handle c, int i, Vertex_handle v) const {
		return _tds.has_vertex(c, i, v);
	}
  
	bool are_equal(Face_handle c, int i, Face_handle n, int j) const {
		return _tds.are_equal(c, i, n, j);
	}

	bool are_equal(const Face & f, const Face & g) const {
		return _tds.are_equal(f.first, f.second, g.first, g.second);
	}

	bool are_equal(const Face & f, Face_handle n, int j) const {
		return _tds.are_equal(f.first, f.second, n, j);
	}


protected:

	Orientation
	inexact_orientation(const Point &p, const Point &q, const Point &r) const {
		return geom_traits().orientation_2_object(p, q, r);	
	}

	Orientation
	inexact_orientation(const Point    &p, const Point    &q,
						const Point    &r, const Point    &s,
						const Offset& o_p, const Offset& o_q,
						const Offset& o_r, const Offset& o_s) const {
		return inexact_orientation(	construct_point(p, o_p),
									construct_point(q, o_q),
									construct_point(r, o_r),
									construct_point(s, o_s) );
	}


protected:
  template< class Conflict_tester, class OutputEdgeIterator>
  void find_in_conflict(Face_handle fh, const Offset& noff, const Conflict_tester& tester, OutputEdgeIterator it) const;

  template< class Conflict_tester>
  void find_in_conflict(Face_handle fh, const Offset& noff, const Conflict_tester& tester, vector<Face_handle>& it) const;

  template<class OutputFaceIterator>
  void 
  find_conflicts( Face_handle         d, 
				  const Point&        pt, 
				  const Offset&       current_off,
				  OutputFaceIterator  it ) const;


protected:
	Bounded_side _side_of_circle( const Face_handle &f, const Point &q,
								  const Offset &offset ) const {
	  Point p[] = { f->vertex(0)->point(),
					f->vertex(1)->point(),
					f->vertex(2)->point()};
	  
	  Offset o[]= { offset.append(f->offset(0)),
					offset.append(f->offset(1)),
					offset.append(f->offset(2))};

	  Oriented_side os = this->side_of_oriented_circle( p[0], p[1], p[2], q,
														o[0], o[1], o[2], Offset());

	  if (os == ON_NEGATIVE_SIDE) {
		return ON_UNBOUNDED_SIDE;
	  }
	  else if (os == ON_POSITIVE_SIDE) {
		return ON_BOUNDED_SIDE;
	  }
	  else {
		return ON_BOUNDARY;
	  }
   }

// Should be 'private'
public:
	Vertex_handle create_initial_triangulation(const Point &p);

public:

	Face_handle euclidean_visibility_locate(const Point& p, Offset& lo, Face_handle f = Face_handle()) const;
	Face_handle hyperbolic_locate(const Point& p, Offset& lo, Locate_type& lt, int& li, int& lj, Face_handle start = Face_handle()) const;

protected:
	// COMMON INSERTION for DELAUNAY and REGULAR TRIANGULATION
	template < class Conflict_tester, class Point_hider >
	Vertex_handle insert_in_conflict(const Point 		   &p, 		Face_handle start,
									 const Conflict_tester &tester, Point_hider &hider) {
		
		Locate_type lt = Locate_type();
		int li = 0, lj = 0;
		Face_handle c = periodic_locate(p, Offset(), lt, li, lj, start);
		return insert_in_conflict(p,lt,c,li,lj,tester,hider);
	}


public:
	Vertex_iterator vertices_begin() const {
		return _tds.vertices_begin();
	}

	Vertex_iterator vertices_end() const {
		return _tds.vertices_end();
	}

	Edge_iterator edges_begin() const {
		return _tds.edges_begin();
	}

	Edge_iterator edges_end() const {
		return _tds.edges_end();
	}

	Face_iterator faces_begin() const {
		return _tds.faces_begin();
	}

	Face_iterator faces_end() const {
		return _tds.faces_end();
	}


  // Circulators

	Face_circulator incident_faces(const Edge & e) const {
		return _tds.incident_faces(e);
	}

	Face_circulator incident_faces(Face_handle c, int i, int j) const {
		return _tds.incident_faces(c, i, j);
	}

	Face_circulator incident_faces(const Edge & e, const Face & start) const {
		return _tds.incident_faces(e, start);
	}

	Face_circulator incident_faces(Face_handle c, int i, int j, const Face & start) const {
		return _tds.incident_faces(c, i, j, start);
	}

	Face_circulator incident_faces(const Edge & e, Face_handle start, int f) const {
		return _tds.incident_faces(e, start, f);
	}

	Face_circulator incident_faces(Face_handle c, int i, int j, Face_handle start, int f) const {
		return _tds.incident_faces(c, i, j, start, f);
	}

	// around a vertex
	template <class OutputIterator>
	OutputIterator incident_faces(Vertex_handle v, OutputIterator faces) const {
		return _tds.incident_faces(v, faces);
	}

	template <class OutputIterator>
	OutputIterator incident_edges(Vertex_handle v, OutputIterator edges) const {
		return _tds.incident_edges(v, edges);
	}

	//deprecated, don't use anymore
	template <class OutputIterator>
	OutputIterator incident_vertices(Vertex_handle v, OutputIterator vertices) const {
		return _tds.incident_vertices(v, vertices);
	}

	
	template <class OutputIterator>
	OutputIterator incident_vertices(Vertex_handle v) const {
	  return _tds.incident_vertices(v, v->face());
	}


	Vertex_circulator incident_vertices(Vertex_handle v) const {
	  return _tds.incident_vertices(v, v->face());
	}


	size_type degree(Vertex_handle v) const {
		return _tds.degree(v);
	}

	// Functions forwarded from TDS.
	int mirror_index(Face_handle c, int i) const {
		return _tds.mirror_index(c, i);
	}

	Vertex_handle mirror_vertex(Face_handle c, int i) const {
		return _tds.mirror_vertex(c, i);
	}

	Face mirror_face(Face f) const {
		return _tds.mirror_face(f);
	}

private:
	bool has_self_edges() const {
		Face_iterator it;
		for ( it = faces_begin(); it != faces_end(); ++it )
			if (has_self_edges(it)) return true;
				return false;
	}
  
	bool has_self_edges(Face_handle c) const;

	bool has_cycles_length_2() const {
	  Vertex_iterator it;
	  for (it = vertices_begin(); it!= vertices_end(); ++it) {
		if (has_cycles_length_2(it)) {
		  return true;
		}
	  }
	  return false;
	}

	bool has_cycles_length_2(Vertex_handle v) const;

public:
	bool is_valid(bool verbose = false, int level = 0) const;
	bool is_valid(Face_handle c, bool verbose = false, int level = 0) const;

	// These functions give the pair (vertex, offset) that corresponds to the
	// i-th vertex of cell ch or vertex vh, respectively.
	void get_vertex(Face_handle ch, int i, Vertex_handle &vh, Offset &off) const;
	void get_vertex(Vertex_handle vh_i, Vertex_handle &vh, Offset &off) const;

protected:
	// Auxiliary functions
	Face_handle get_face(const Vertex_handle* vh) const;
  

public:

  int dimension() const {
	return _tds.dimension();
  }

  Object
  dual(const Edge_iterator& ei) const
  {
	return this->dual(*ei);
  }
  
  Object
  dual(const Edge &e) const
  {
	
	if(this->dimension() == 1) {
	  const Point& p = (e.first)->vertex(cw(e.second))->point();
	  const Point& q = (e.first)->vertex(ccw(e.second))->point();
	  
	  // hyperbolic line
	  Segment line = this->geom_traits().construct_hyperbolic_bisector_2_object()(p,q);
	  return make_object(line);
	}
	
	// incident faces
	Face_handle f1 = e.first;
	int i1 = e.second;
	
	Face_handle f2 = f1->neighbor(i1);
	int i2 = f2->index(f1);
	  
	Segment s = this->geom_traits().construct_segment_2_object()
	(dual(f1),dual(f2));

	return make_object(s);
	
  }


}; // class Periodic_4_hyperbolic_triangulation_2


/*********** FUNCTION IMPLEMENTATIONS *************/


template < class GT, class TDS >
inline typename Periodic_4_hyperbolic_triangulation_2<GT,TDS>::Periodic_triangle
Periodic_4_hyperbolic_triangulation_2<GT,TDS>::
periodic_triangle(const Face_handle c, int i) const { 
	CGAL_triangulation_precondition( number_of_vertices() != 0 );
	CGAL_triangulation_precondition( i >= 0 && i <= 2 );

	return make_array(	std::make_pair(c->vertex( (i+0)&3 )->point(), c->offset((i+0)&3)),
						std::make_pair(c->vertex( (i+1)&3 )->point(), c->offset((i+1)&3)),
						std::make_pair(c->vertex( (i+2)&3 )->point(), c->offset((i+2)&3)) );
}



/// tests if two vertices of one cell are just periodic copies of each other

template < class GT, class TDS >
inline bool Periodic_4_hyperbolic_triangulation_2<GT, TDS>::
has_self_edges(typename TDS::Face_handle c) const {
	CGAL_triangulation_assertion((c->vertex(0) != c->vertex(1)) || (c->offset(0) != c->offset(1)));
	CGAL_triangulation_assertion((c->vertex(0) != c->vertex(2)) || (c->offset(0) != c->offset(2)));
	CGAL_triangulation_assertion((c->vertex(1) != c->vertex(2)) || (c->offset(1) != c->offset(2)));

	return ( (c->vertex(0) == c->vertex(1)) || 
				 (c->vertex(0) == c->vertex(2)) ||
			   (c->vertex(1) == c->vertex(2))   );
}


template < class GT, class TDS >
inline bool Periodic_4_hyperbolic_triangulation_2<GT, TDS>::
has_cycles_length_2(typename TDS::Vertex_handle v) const {
  typename TDS::Vertex_circulator vc = incident_vertices(v), done(vc);
  std::set<Vertex_handle> check;
  do {
	std::pair<typename std::set<Vertex_handle>::iterator, bool> res = check.insert(vc);
	if (!res.second) {
	  return true;
	} 
  } while (++vc != done);
  return false;
}
	
/*! \brief Tests if the triangulation is valid.
 *
 * A triangulation is valid if
 * - A cell is not its own neighbor.
 * - A cell has no two equal neighbors
 * - A cell has no two equal vertex-offset pairs
 * - A cell is positively oriented.
 * - The point of a neighbor of cell c that does not belong to c is not inside
 *   the circumcircle of c.
 */
template < class GT, class TDS >
bool
Periodic_4_hyperbolic_triangulation_2<GT,TDS>::
is_valid(bool verbose, int level) const {

	bool error = false;
	for (Face_iterator cit = faces_begin(); cit != faces_end(); ++cit) {
		for (int i=0; i<3; i++) {
			CGAL_triangulation_assertion(cit != cit->neighbor(i));
			for (int j=i+1; j<3; j++) {
			  CGAL_triangulation_assertion(cit->neighbor(i) != cit->neighbor(j));
				  CGAL_triangulation_assertion(cit->vertex(i) != cit->vertex(j));
			}
		}

		if (orientation( cit->vertex(0)->point(), cit->vertex(1)->point(), cit->vertex(2)->point(), 
						 cit->offset(0),          cit->offset(1),          cit->offset(2) ) != POSITIVE) {
		if (verbose) {
		  std::cerr << "Orientation failed for face " << cit->get_number() << std::endl;
		  for (int j = 0; j < 3; j++) {
			std::cerr << "   v" << j << " is " << cit->vertex(j)->idx() << " with offset " << cit->offset(j) << std::endl;
		  }
		}
		error = true;
		}
	}

	if (!error) {
	  CGAL_triangulation_assertion(!has_self_edges());
	  CGAL_triangulation_assertion(!has_cycles_length_2());
	  CGAL_triangulation_assertion( _tds.number_of_vertices() + _tds.number_of_faces() + 2 == _tds.number_of_edges() );
	  return true; 
	} else {
	  return false;
	}
}


template < class GT, class TDS >
bool Periodic_4_hyperbolic_triangulation_2<GT,TDS>::
is_valid(Face_handle ch, bool verbose, int level) const {

	if ( ! _tds.is_valid(ch,verbose,level) )
		return false;
	
	bool error = false;
	const Point *p[3]; Offset off[3];
	
	for (int i=0; i<3; i++) {
		p[i] = &ch->vertex(i)->point();
		off[i] = ch->offset(i);
	}
	
	if (orientation( *p[0],  *p[1],  *p[2],  
					off[0], off[1], off[2]) != POSITIVE) {
		error = true;
	}
  
  return !error;
}



/*********************************************************************************/


template <class GT, class TDS>
template <class OutputFaceIterator>
void
Periodic_4_hyperbolic_triangulation_2<GT,TDS>::
find_conflicts( Face_handle         d, 
				const Point&        pt, 
				const Offset&       current_off,
				OutputFaceIterator  it ) const {
  if (d->tds_data().is_clear()) {
	if (_side_of_circle(d, pt, current_off) == ON_BOUNDED_SIDE) {
	  d->tds_data().mark_in_conflict();
	  d->store_offsets(current_off);
	  it++ = d;
	  for (int jj = 0; jj < 3; jj++) {
		if (d->neighbor(jj)->tds_data().is_clear()) {
		  find_conflicts(d->neighbor(jj), pt, current_off.append(d->neighbor_offset(jj)), it);
		}
	  }
	}
  }
}



/*********************************************************************************/


template <class GT, class TDS>
typename TDS::Face_handle Periodic_4_hyperbolic_triangulation_2<GT, TDS>::
euclidean_visibility_locate(const Point& p, Offset& loff, Face_handle f) const
{

  typedef typename GT::Side_of_fundamental_octagon Side_of_fundamental_octagon;
  
  Side_of_fundamental_octagon check = Side_of_fundamental_octagon();
  CGAL::Bounded_side side = check(p);
  if (side != ON_BOUNDED_SIDE) {
	return Face_handle();
  }

  // Handle the case where an initial Face_handle is not given
  if (f == Face_handle()) {
	f = _tds.faces().begin();
  }

  int curr = 0;
  int succ = ccw(curr);
  int counter = 0;

	while (true) {
	Orientation o = orientation(f->vertex(curr)->point(),       f->vertex(succ)->point(),       p,
								loff.append(f->offset(curr)),   loff.append(f->offset(succ)),   Offset());
	  
	if (o == NEGATIVE) {
	  loff = loff.append(f->neighbor_offset(cw(curr))); 
	  f = f->neighbor(cw(curr));
	  curr = ccw(curr);
	  succ = ccw(curr);
	  counter = 0;
	} else {
	  curr = succ;
	  succ = ccw(curr);
	  counter++;
	}

	if (counter > 2) {
	  break;
	}

  }

  return f;
}




template <class GT, class TDS>
typename TDS::Face_handle Periodic_4_hyperbolic_triangulation_2<GT, TDS>::
hyperbolic_locate(const Point& p, Offset& lo, Locate_type& lt, int& li, int& lj, Face_handle start) const
{
	// Get a hint of where the point is located. It's either in lf or in one of its neighbors.
	Face_handle lf = euclidean_visibility_locate(p, lo, start);
	typedef typename GT::Side_of_hyperbolic_face_2 Side_of_hyperbolic_face_2;
	Side_of_hyperbolic_face_2 sf;

	cout << "Euclidean location says face " << lf->get_number() << " with offset " << lo << endl;

	Bounded_side bs = sf(p, lf, lo);
	if (bs == ON_BOUNDED_SIDE) {
		cout << "Location was correct! The point is inside the face!" << endl;
	} else if (bs == ON_BOUNDARY) {
		cout << "Not clear if correct! The point is on the boundary of the face!" << endl;
	} else {
		Bounded_side bs1 = sf(p, lf->neighbor(0), lo);
		Bounded_side bs2 = sf(p, lf->neighbor(1), lo);
		Bounded_side bs3 = sf(p, lf->neighbor(2), lo);
		cout << "The location was not correct! The point is found inside ";
		if (bs1 == ON_BOUNDED_SIDE) 
			cout << " neighbor 0! ";
		if (bs2 == ON_BOUNDED_SIDE) 
			cout << " neighbor 1! ";
		if (bs3 == ON_BOUNDED_SIDE) 
			cout << " neighbor 2! ";
		cout << endl;
	}


 //  Orientation o1 = orientation(f->vertex(0)->point(), f->vertex(1)->point(), p,
	// 						   f->offset(0),          f->offset(1),          Offset()); //bu_encode_offset());

 //  Orientation o2 = orientation(f->vertex(1)->point(), f->vertex(2)->point(), p,
	// 						   f->offset(1),          f->offset(2),          Offset()); //bu_encode_offset());

 //  Orientation o3 = orientation(f->vertex(2)->point(), f->vertex(0)->point(), p,
	// 						   f->offset(2),          f->offset(0),          Offset()); //bu_encode_offset());
  
 //  int sum = (o1 == COLLINEAR) + (o2 == COLLINEAR) + (o3 == COLLINEAR);
 //  if (sum == 0) {
	// lt = FACE;
 //  }
 //  else if (sum == 1) {
	// lt = EDGE; 
	// li = ( o1 == COLLINEAR ? 2 : (o2 == COLLINEAR ? 0 : 1) );
 //  } else {
	// lt = VERTEX;
	// li = ( o1 != COLLINEAR ? 2 : (o2 != COLLINEAR ? 0 : 1) );
 //  }


	return lf;
}






}  // namespace CGAL


#endif   // CGAL_PERIODIC_4_HYPERBOLIC_TRIANGULATION_2_H



