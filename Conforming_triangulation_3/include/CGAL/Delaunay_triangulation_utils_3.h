//Some additional utilities for the Delaunay triangulation structure.
//Copyright (C) 2012  Utrecht University
//
//This program is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//This program is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// Author(s): Thijs van Lankveld

// Based on CGAL/Delaunay_triangulation_3.h
//
// For the love of something, why is the 3D triangulation package so badly written?

#ifndef CGAL_DELAUNAY_TRIANGULATION_3_UTILS_H
#define CGAL_DELAUNAY_TRIANGULATION_3_UTILS_H

#include <CGAL/basic.h>

#include <utility>
#include <vector>

#include <CGAL/Delaunay_triangulation_3.h>
#include "CCDT/Triangulation_segment_traverser_3.h"

CGAL_BEGIN_NAMESPACE

template < class DT > class Natural_neighbors_3;

template < class Gt, class Tds > class Triangulation_cell_traverser_3;

template < class Gt,
           class Tds = Triangulation_data_structure_3 < Triangulation_vertex_base_3<Gt>, Triangulation_cell_base_3<Gt> > >
class Delaunay_triangulation_utils_3: public Delaunay_triangulation_3<Gt, Tds> {
	typedef Triangulation_data_structure						Tds;

	typedef Delaunay_triangulation_utils_3<Gt, Tds>				Self;
	typedef Delaunay_triangulation_3<Gt, Tds>					DT;
	typedef Triangulation_3<Gt,Tds>								Tr;

	friend class Natural_neighbors_3<Self>;

public:
	typedef Tds													Triangulation_data_structure;
	typedef Gt													Geom_traits;

	typedef typename Gt::FT										FT;

	typedef typename Gt::Point_3								Point;
	typedef typename Gt::Vector_3								Vector;
	typedef typename Gt::Segment_3								Segment;
	typedef typename Gt::Line_3									Line;
	typedef typename Gt::Triangle_3								Triangle;
	typedef typename Gt::Tetrahedron_3							Tetrahedron;
	typedef typename Gt::Plane_3								Plane;

	typedef typename DT::Cell_handle							Cell_handle;
	typedef typename DT::Vertex_handle							Vertex_handle;

	typedef typename DT::Cell									Cell;
	typedef typename DT::Vertex									Vertex;
	typedef typename DT::Facet									Facet;
	typedef typename DT::Edge									Edge;

	typedef std::pair<Vertex_handle, Vertex_handle>				Bi_vertex;
	typedef Triple<Vertex_handle, Vertex_handle, Vertex_handle>	Tri_vertex;

	typedef std::pair<Point, Point>								Constraint_1;
	typedef std::vector<Point>									Constraint_2;

	typedef typename DT::Cell_circulator						Cell_circulator;
	typedef typename DT::Facet_circulator						Facet_circulator;
	typedef typename DT::Cell_iterator							Cell_iterator;
	typedef typename DT::Facet_iterator							Facet_iterator;
	typedef typename DT::Edge_iterator							Edge_iterator;
	typedef typename DT::Vertex_iterator						Vertex_iterator;

	typedef typename DT::Finite_vertices_iterator				Finite_vertices_iterator;
	typedef typename DT::Finite_cells_iterator					Finite_cells_iterator;
	typedef typename DT::Finite_facets_iterator					Finite_facets_iterator;
	typedef typename DT::Finite_edges_iterator					Finite_edges_iterator;

	typedef typename DT::All_cells_iterator						All_cells_iterator;

	typedef Triangulation_segment_traverser_3<Gt,Tds>			Cell_traverser;

	typedef typename DT::Locate_type							Locate_type;

#ifndef CGAL_CFG_USING_BASE_MEMBER_BUG_2
	using DT::cw;
	using DT::ccw;
	using DT::coplanar;
	using DT::coplanar_orientation;
	using DT::dimension;
	using DT::finite_facets_begin;
	using DT::finite_facets_end;
	using DT::finite_vertices_begin;
	using DT::finite_vertices_end;
	using DT::finite_cells_begin;
	using DT::finite_cells_end;
	using DT::finite_edges_begin;
	using DT::finite_edges_end;
	using DT::geom_traits;
	using DT::infinite_vertex;
	using DT::is_valid;
	using DT::mirror_vertex;
	using DT::next_around_edge;
	using DT::number_of_vertices;
	using DT::orientation;
	using DT::remove;
	using DT::coplanar_side_of_bounded_circle;
	using DT::side_of_oriented_sphere;
	using DT::side_of_segment;
	using DT::tds;
	using DT::vertex_triple_index;
#endif

protected:
	typedef std::list<Tri_vertex>								Tri_vertex_collection;
	typedef std::list<Bi_vertex>								Bi_vertex_collection;

protected:
	// Test whether a newly inserted point conflicts with the existing cells.
	class Conflict_tester_3 {
	protected:
		const Point& _p;
		const Self* _tr;

	public:
		Conflict_tester_3(const Point& pt, const Self* tr): _p(pt), _tr(tr) {}
		bool operator()(const Cell_handle c) const {return _tr->side_of_sphere(c, _p, true) == ON_BOUNDED_SIDE;}
		Oriented_side compare_weight(const Point&, const Point&) const {return ZERO;}
		bool test_initial_cell(Cell_handle) const {return true;}
	}; // class Conflict_tester_3

	// Test whether a newly inserted point conflicts with the existing cells.
	class Conflict_tester_2 {
	protected:
		const Point& _p;
		const Self* _tr;

	public:
		Conflict_tester_2(const Point& pt, const Self* tr): _p(pt), _tr(tr) {}
		bool operator()(const Cell_handle c) const {return _tr->side_of_circle(c, 3, _p, true) == ON_BOUNDED_SIDE;}
		Oriented_side compare_weight(const Point&, const Point&) const {return ZERO;}
		bool test_initial_cell(Cell_handle) const {return true;}
	}; // class Conflict_tester_2

	class Hidden_point_visitor {
	public:
		Hidden_point_visitor() {}

		template < class InputIterator >
		void process_cells_in_conflict(InputIterator, InputIterator) const {}
		void reinsert_vertices(Vertex_handle) {}
		Vertex_handle replace_vertex(Cell_handle c, int index, const Point&) {return c->vertex(index);}
		void hide_point(Cell_handle, const Point&) {}
	}; // class Hidden_point_visitor
	
	Hidden_point_visitor hidden_point_visitor;

public:
	// Constructors for bi-vertices.
	inline Bi_vertex to_bi_vertex(Cell_handle c, int li, int lj) const {return Bi_vertex(c->vertex(li), c->vertex(lj));}
	inline Bi_vertex to_bi_vertex(const Edge& e) const {return to_bi_vertex(e.first, e.second, e.third);}
	inline Bi_vertex sort_bi_vertex(Vertex_handle v1, Vertex_handle v2) const {if (v1 < v2) return Bi_vertex(v1, v2); return Bi_vertex(v2, v1);}
	inline Bi_vertex sort_bi_vertex(const Bi_vertex& bv) const {return sort_bi_vertex(bv.first, bv.second);}
	inline Bi_vertex sort_bi_vertex(Cell_handle c, int li, int lj) const {return sort_bi_vertex(c->vertex(li), c->vertex(lj));}
	inline Bi_vertex sort_bi_vertex(const Edge& e) const {return sort_bi_vertex(e.first, e.second, e.third);}

	// Constructors for tri-vertices.
	inline Tri_vertex to_tri_vertex(Cell_handle c, int li) const {return Tri_vertex(c->vertex((li+1)&3), c->vertex((li+2)&3), c->vertex((li+3)&3));}
	inline Tri_vertex to_tri_vertex(const Facet& f) const {return to_tri_vertex(f.first, f.second);}
	inline Tri_vertex sort_tri_vertex(Vertex_handle v1, Vertex_handle v2, Vertex_handle v3) const {
		CGAL_triangulation_precondition(v1 != v2 && v1 != v3 && v2 != v3);
		if (v1 < v2) {
			if (v2 < v3)
				return Tri_vertex(v3, v2, v1);
			if (v3 < v1)
				return Tri_vertex(v2, v1, v3);
			return Tri_vertex(v2, v3, v1);
		}
		else { // v2 < v1
			if (v3 < v2)
				return Tri_vertex(v1, v2, v3);
			if (v1 < v3)
				return Tri_vertex(v3, v1, v2);
			return Tri_vertex(v1, v3, v2);
		}
	}
	inline Tri_vertex sort_tri_vertex(const Tri_vertex& tv) const {return sort_tri_vertex(tv.first, tv.second, tv.third);}
	inline Tri_vertex sort_tri_vertex(Cell_handle c, int li) const {return sort_tri_vertex(to_tri_vertex(c, li));}
	inline Tri_vertex sort_tri_vertex(const Facet& f) const {return sort_tri_vertex(f.first, f.second);}
	inline Tri_vertex oriented_tri_vertex(Cell_handle c, int li) const {
		CGAL_triangulation_precondition(dimension() == 2 || dimension() == 3);
		CGAL_triangulation_precondition((dimension() == 2 && li == 3) || (dimension() == 3 && li >= 0 && li <= 3) );
		CGAL_triangulation_precondition(!is_infinite(c, li));
		if ((li&1) == 0)
			return Tri_vertex(c->vertex((li+2)&3),
							  c->vertex((li+1)&3),
							  c->vertex((li+3)&3));
		return Tri_vertex(c->vertex((li+1)&3),
						  c->vertex((li+2)&3),
						  c->vertex((li+3)&3));
	}
	inline Tri_vertex oriented_tri_vertex(const Facet& f) const {return oriented_tri_vertex(f.first, f.second);}

protected:
	// tri_vertex methods similar to Triangulation_3.
	Tri_vertex make_tri_vertex(const Facet& f) const {
		return Tri_vertex(f.first->vertex(vertex_triple_index(f.second, 0)),
						  f.first->vertex(vertex_triple_index(f.second, 1)),
						  f.first->vertex(vertex_triple_index(f.second, 2)));
	}
	void make_canonical(Tri_vertex& t) const {
		int i = (&*(t.first) < &*(t.second))? 0 : 1;
		if (i == 0)
			i = (&*(t.first) < &*(t.third))? 0 : 2;
		else
			i = (&*(t.second) < &*(t.third))? 1 : 2;
		Vertex_handle tmp;
		switch (i) {
			case 0: return;
			case 1:
				tmp = t.first;
				t.first = t.second;
				t.second = t.third;
				t.third = tmp;
				return;
			default:
				tmp = t.first;
				t.first = t.third;
				t.third = t.second;
				t.second = tmp;
		}
	}

	// PREDICATES
	Orientation	coplanar_orientation(const Point& p0, const Point& p1, const Point& p2, const Point& p) const {
		return Gt().coplanar_orientation_3_object()(p0, p1, p2, p);
	}

public:
	// CONSTRUCTORS
	Delaunay_triangulation_utils_3(const Gt& gt = Gt()): DT(gt) {}

	// Copy constructor duplicates vertices and cells.
	Delaunay_triangulation_utils_3(const Delaunay_triangulation_utils_3& tr): DT(tr) {
		CGAL_triangulation_postcondition(is_valid());
	}

	// Create a 3D Delaunay triangulation from a number of points.
	template < typename InputIterator >
	Delaunay_triangulation_utils_3(InputIterator first, InputIterator last, const Gt& gt = Gt()): DT(gt) {
		insert(first, last);
		CGAL_triangulation_postcondition(is_valid());
	}

public:
	// INSERT POINT

	Vertex_handle insert(const Point& p, Cell_handle start = Cell_handle());
	virtual Vertex_handle insert(const Point& p, Locate_type lt, Cell_handle c, int li, int lj);

	template < class InputIterator >
	int insert(InputIterator first, InputIterator last) {
		int n = number_of_vertices();

		std::vector<Point> points(first, last);
		std::random_shuffle(points.begin(), points.end());
		spatial_sort(points.begin(), points.end(), geom_traits());

		Cell_handle hint;
		for (std::vector<Point>::const_iterator it = points.begin(), end = points.end(); it != end; ++it)
			hint = insert(*it, hint)->cell();

		return number_of_vertices() - n;
	}

public:
	// MOVE

	Vertex_handle move_point(Vertex_handle v, const Point& p);

public:
	// Construct the edge of c opposite to ij.
	Edge opposite_edge(Cell_handle c, int li, int lj) const {
		CGAL_triangulation_precondition(li >= 0 && li < 4);
		CGAL_triangulation_precondition(lj >= 0 && lj < 4);
		CGAL_triangulation_precondition(li != lj);

		switch (6-li-lj) { // i + j + missing indices = 6.
			case 1: return Edge(c, 0, 1);
			case 2: return Edge(c, 0, 2);
			case 3: return (li == 0 || lj == 0)
						   ? Edge(c, 1, 2)
						   : Edge(c, 0, 3);
			case 4: return Edge(c, 1, 3);
			case 5: return Edge(c, 2, 3);
		}

		CGAL_triangulation_assertion(false);
		return Edge();
	}
	Edge opposite_edge(const Edge& e) const {return opposite_edge(e.first, e.second, e.third);}

	// Give the same facet as seen from the other side.
	Facet mirror_facet(Cell_handle c, int li) const {return Facet(c->neighbor(li), c->neighbor(li)->index(c));}
	Facet mirror_facet(const Facet& f) const {return mirror_facet(f.first, f.second);}

	// Construct a plane of a facet.
	Plane plane(Cell_handle c, int li) const {
		/* This should be removed to reduce orientation errors with the inexact kernel*/
		CGAL_triangulation_precondition(dimension() >= 2);
		CGAL_triangulation_precondition(dimension() == 3 || li == 3);
		CGAL_triangulation_precondition(li >= 0 && li <= 3);
		CGAL_triangulation_precondition(!is_infinite(c, li));
		if ((li&1) == 0)
			return Plane(c->vertex((li+2)&3)->point(),
						 c->vertex((li+1)&3)->point(),
						 c->vertex((li+3)&3)->point());
		return Plane(c->vertex((li+1)&3)->point(),
					 c->vertex((li+2)&3)->point(),
					 c->vertex((li+3)&3)->point());
	}
	Plane plane(const Facet& f) const {return plane(f.first, f.second);}

	// Construct a line of an edge.
	Line line(Cell_handle c, int li, int lj) const {
		CGAL_triangulation_precondition(dimension() >= 1);
		CGAL_triangulation_precondition(li >= 0 && li <= 3);
		CGAL_triangulation_precondition(lj >= 0 && lj <= 3);
		CGAL_triangulation_precondition(li != lj);
		CGAL_triangulation_precondition(li + lj < dimension() * 2);
		CGAL_triangulation_precondition(!is_infinite(c, li, lj));
		return Line(c->vertex(li)->point(),
					c->vertex(lj)->point());
	}
	Line line(const Edge& e) const {return line(e.first, e.second, e.third);}
}; // Delaunay_triangulation_utils_3

template < class Gt, class Tds >
typename Delaunay_triangulation_utils_3<Gt,Tds>::Vertex_handle
Delaunay_triangulation_utils_3<Gt,Tds>::
insert(const Point& p, Cell_handle start) {
	Locate_type lt;
	int li, lj;
	Cell_handle c = locate(p, lt, li, lj, start);
	return insert(p, lt, c, li, lj);
}

template < class Gt, class Tds >
typename Delaunay_triangulation_utils_3<Gt,Tds>::Vertex_handle
Delaunay_triangulation_utils_3<Gt,Tds>::
insert(const Point& p, Locate_type lt, Cell_handle c, int li, int lj) {
	switch (dimension()) {
		case 3: {
			Conflict_tester_3 tester(p, this);
			Vertex_handle v = insert_in_conflict(p, lt, c, li, lj, tester, hidden_point_visitor);
			return v;
		}// dim 3
		case 2: {
			Conflict_tester_2 tester(p, this);
			return insert_in_conflict(p, lt, c, li, lj, tester, hidden_point_visitor);
		}//dim 2
		default :
			// dimension <= 1
			// Do not use the generic insert.
			return Tr::insert(p, c);
	}
}

template < class Gt, class Tds >
typename Delaunay_triangulation_utils_3<Gt,Tds>::Vertex_handle
Delaunay_triangulation_utils_3<Gt,Tds>::move_point(Vertex_handle v, const Point& p) {
	CGAL_triangulation_precondition(!is_infinite(v));
	CGAL_triangulation_expensive_precondition(is_vertex(v));

	// Remember an incident vertex to restart the point location after the removal.
	Cell_handle c = v->cell();
	Vertex_handle old_neighbor = c->vertex(c->index(v) == 0 ? 1 : 0);
	CGAL_triangulation_assertion(old_neighbor != v);

	if (!remove(v))
		return v;

	if (dimension() <= 0)
		return insert(p);
	return insert(p, old_neighbor->cell());
}

CGAL_END_NAMESPACE

#endif // CGAL_DELAUNAY_TRIANGULATION_3_UTILS_H
