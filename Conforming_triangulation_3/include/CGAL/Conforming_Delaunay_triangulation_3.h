//The conforming Delaunay triangulation structure.
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

#ifndef CGAL_CONFORMING_DELAUNAY_TRIANGULATION_3_H
#define CGAL_CONFORMING_DELAUNAY_TRIANGULATION_3_H

#include <CGAL/basic.h>

#include <utility>
#include <vector>

#include "Delaunay_triangulation_utils_3.h"
#include "Conforming_triangulation_cell_base_3.h"
#include "Conforming_triangulation_vertex_base_3.h"
#include "Encroaching_collecter_3.h"

namespace CGAL {

template < class Tr > class Natural_neighbors_3;

template < class Gt,
           class Tds = Triangulation_data_structure_3 < Conforming_triangulation_vertex_base_3<Gt>,
														Conforming_triangulation_cell_base_3<Gt> >,
		   class Itag = No_intersection_tag >
class Conforming_Delaunay_triangulation_3: public Delaunay_triangulation_utils_3<Gt,Tds> {
	//typedef Triangulation_data_structure						Tds;

	typedef Conforming_Delaunay_triangulation_3<Gt,Tds,Itag>	cDT;
	typedef Delaunay_triangulation_utils_3<Gt,Tds>				DT;
	typedef Triangulation_3<Gt,Tds>								Tr;

	friend class Natural_neighbors_3<cDT>;
        typedef typename Gt::Line_3 Line;
public:
	typedef Tds													Triangulation_data_structure;
	typedef Gt													Geom_traits;
	typedef Itag												Intersection_tag;

	typedef typename Gt::FT										FT;

	typedef typename Gt::Point_3								Point;
	typedef typename Gt::Vector_3								Vector;
	typedef typename Gt::Segment_3								Segment;
	typedef typename Gt::Triangle_3								Triangle;
	typedef typename Gt::Tetrahedron_3							Tetrahedron;
	typedef typename Gt::Plane_3								Plane;

	typedef typename DT::Cell_handle							Cell_handle;
	typedef typename DT::Vertex_handle							Vertex_handle;

	typedef typename DT::Cell									Cell;
	typedef typename DT::Vertex									Vertex;
	typedef typename DT::Facet									Facet;
	typedef typename DT::Edge									Edge;

	typedef typename DT::Bi_vertex								Bi_vertex;

	typedef std::pair<Point, Point>								Constraint_1;

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

	typedef typename DT::Cell_traverser							Cell_traverser;

	typedef typename DT::Locate_type							Locate_type;

private:
	typedef std::list<Point>									Point_list;
	typedef typename Point_list::iterator						Points_iterator;
	typedef std::back_insert_iterator<Point_list>				Points_output;

public:
	typedef Encroaching_collecter_3<Gt,Tds,Itag,Points_output>	Encroaching_collecter;

#ifndef CGAL_CFG_USING_BASE_MEMBER_BUG_2
	using DT::cw;
	using DT::ccw;
	using DT::geom_traits;
	using DT::number_of_vertices;
	using DT::dimension;
	using DT::finite_facets_begin;
	using DT::finite_facets_end;
	using DT::finite_vertices_begin;
	using DT::finite_vertices_end;
	using DT::finite_cells_begin;
	using DT::finite_cells_end;
	using DT::finite_edges_begin;
	using DT::finite_edges_end;
	using DT::tds;
	using DT::infinite_vertex;
	using DT::insert;
	using DT::next_around_edge;
	using DT::vertex_triple_index;
	using DT::mirror_vertex;
	using DT::move_point;
	using DT::coplanar;
	using DT::coplanar_orientation;
	using DT::orientation;
	using DT::is_infinite;
	using DT::finite_incident_edges;
	using DT::incident_cells;
	using DT::insert_in_conflict;
	using DT::is_valid_finite;
#endif

protected:
	// Visitor that should be run before and after inserting a point
	// to maintain the conforming edges in the triangulation.
	template < class Tr >
	class Conforming_visitor {
		typedef std::pair<Bi_vertex, bool>	CBV;

		mutable std::list<CBV> _conforming;

	protected:
		Tr *_tr;

	public:
		Conforming_visitor(Tr* tr): _tr(tr) {}

		// Store and remove the conforming edges from a collection of cells.
		template <class InputIterator>
		void process_cells_in_conflict(InputIterator start, InputIterator end) const {
			for (InputIterator it = start; it != end; ++it) {
				for (int i = 0; i < _tr->dimension(); ++i) {
					for (int j = i+1; j < _tr->dimension()+1; ++j) {
						if ((*it)->is_conforming(i, j)) {
							_conforming.push_back(CBV(_tr->to_bi_vertex(*it, i, j), (*it)->is_marked(i, j)));
							_tr->set_conforming(*it, i, j, false);
						}
					}
				}
			}
		}

		// Reinsert the conforming edges into the triangulation.
		void reinsert_vertices(Vertex_handle) {
			while (!_conforming.empty()) {
				CBV b = _conforming.front();
				_conforming.pop_front();
				if (b.second)
					_tr->insert_marked(b.first);
				else
					_tr->insert_conforming(b.first);
			}
		}

		// Reinsert the conforming edges into the triangulation
		// under the assumption that no Steiner points need to be added.
		void reinsert_no_steiner() {
			Cell_handle c; int li, lj;
			while (!_conforming.empty()) {
				CBV b = _conforming.front();
				_conforming.pop_front();
				CGAL_triangulation_assertion_code(bool ok =)
					_tr->is_edge(b.first.first, b.first.second, c, li, lj);
				CGAL_triangulation_assertion(ok);

				if (b.second)
					_tr->mark_edge(c, li, lj);
				else
					_tr->set_conforming(c, li, lj, true);
			}
		}

		Vertex_handle replace_vertex(Cell_handle c, int index, const Point&) {return c->vertex(index);}
		void hide_point(Cell_handle, const Point&) {}
	}; // class Conforming_visitor

	friend class Conflict_tester_3;
	friend class Conflict_tester_2;
	friend class Encroaching_collecter_3<Gt,Tds,Itag,Points_output>;

	Conforming_visitor<cDT> conforming_visitor;

protected:
	// PREDICATES

	// Checks if point c is encroaching upon the segment ab.
	// A point is encroaching if it is inside the ball with diameter ab.
	bool is_encroaching(const Point& a, const Point& b, const Point& c) const {
		CGAL_triangulation_assertion(a != b);
		CGAL_triangulation_assertion(a != c);
		CGAL_triangulation_assertion(b != c);
		return side_of_bounded_sphere(a, b, c) != ON_UNBOUNDED_SIDE;
	}

public:
	// CONSTRUCTORS
	Conforming_Delaunay_triangulation_3(const Gt& gt = Gt()): DT(gt), conforming_visitor(this) {}

	// Copy constructor duplicates vertices and cells.
	Conforming_Delaunay_triangulation_3(const Conforming_Delaunay_triangulation_3& tr): DT(tr), conforming_visitor(this) {
		CGAL_triangulation_postcondition(is_valid());
	}

	// Create a conforming 3D Delaunay triangulation from a number of constraints.
	template < class InputIterator >
	Conforming_Delaunay_triangulation_3(InputIterator begin, InputIterator end, const Gt& gt = Gt()):
	DT(gt), conforming_visitor(this) {
		insert_conforming(begin, end);
		CGAL_triangulation_postcondition(is_valid());
	}

protected:
	// Compute the intersection of ab with a line or plane.
	bool compute_intersection(const Line& ab, const Line& l, Point& p) const {
		Object result = Gt().intersect_3_object()(ab, l);
		return assign(p, result);
	}
	bool compute_intersection(const Line& ab, const Plane& pl, Point& p) const {
		Object result = Gt().intersect_3_object()(ab, pl);
		return assign(p, result);
	}

	// Compute the intersection of ab with either a line (edge) or plane (facet).
	bool compute_intersection(Locate_type lt, Cell_handle c, int li, int lj, Vertex_handle va, Vertex_handle vb, Point& p) const {
		if (lt == DT::EDGE)
			return compute_intersection(Line(va->point(), vb->point()),
										Line(c->vertex(li)->point(), c->vertex(lj)->point()), p);
		else
			return compute_intersection(Line(va->point(), vb->point()), plane(c, li), p);
	}

	// Insert the intersection of segment ab with either a facet or an edge.
	Vertex_handle intersect(Cell_handle c, Locate_type lt, int li, int lj, Vertex_handle va, Vertex_handle vb) {
		return intersect(c, lt, li, lj, va, vb, Itag());}
	Vertex_handle intersect(Cell_handle c, Locate_type lt, int li, int lj, Vertex_handle va, Vertex_handle vb, No_intersection_tag);
	Vertex_handle intersect(Cell_handle c, Locate_type lt, int li, int lj, Vertex_handle va, Vertex_handle vb, Exact_intersections_tag);
	Vertex_handle intersect(Cell_handle c, Locate_type lt, int li, int lj, Vertex_handle va, Vertex_handle vb, Exact_predicates_tag);

	// Split segment ab based on a reference point.
	Point construct_split_point(Vertex_handle va, Vertex_handle vb, const Point& ref);
	
protected:
	// Conform an edge to important encroaching points. 
        // The boolean indicates whether the edge should be made conforming Gabriel
        template <bool make_gabriel>
	void conform_segment(Vertex_handle va, Vertex_handle vb, Points_iterator encr_begin, Points_iterator encr_end);

	// Mark a segment that it is assumed to already be in the triangulation.
	void mark_segment(Vertex_handle va, Vertex_handle vb);
	// Insert a marked segment;
	// this includes conforming the segment and marking its subsegments as conforming.
	void insert_marked(Vertex_handle va, Vertex_handle vb);
	void insert_marked(const Bi_vertex& m) {insert_marked(m.first, m.second);}

	// Mark an existing edge or facet as (un)constrained.
	void set_conforming(Cell_handle c, int li, int lj, bool C);
	void set_conforming(const Edge& e, bool C) {set_conforming(e.first, e.second, e.third, C);}
	void mark_edge(Cell_handle c, int li, int lj);
	void mark_edge(const Edge& e) {mark_edge(e.first, e.second, e.third);}

protected:
	// Restore the conforming edges of a facet based on the neighboring cell.
	void restore_conforming(Cell_handle c, int li);
	void restore_conforming(const Facet& f) {restore_conforming(f.first, f.second);}

public:
	// INSERT POINT

	virtual Vertex_handle insert(const Point& p, Locate_type lt, Cell_handle c, int li, int lj);

public:
	// INSERT CONFORMING

	// Insert a conforming segment.
	void insert_conforming(Vertex_handle va, Vertex_handle vb) {
		Point_list encr; conform_segment<false>(va, vb, encr.begin(), encr.end());}
        void insert_conforming_Gabriel(Vertex_handle va, Vertex_handle vb) {
		Point_list encr; conform_segment<true>(va, vb, encr.begin(), encr.end());}
	void insert_conforming(const Bi_vertex& c) {
		insert_conforming(c.first, c.second);}
	void insert_conforming(const Point& a, const Point& b, Cell_handle hint = Cell_handle());
	void insert_conforming(const Constraint_1& c, Cell_handle hint = Cell_handle()) {
		insert_conforming(c.first, c.second, hint);}

	// Insert a number of conforming segments, presented as a collection of pairs of Point_3.
	template < class InputIterator >
	int insert_conforming(InputIterator begin, InputIterator end, Cell_handle hint = Cell_handle());
	
	// Insert a loop of conforming segments, presented as a collection of Point_3.
	template < class InputIterator >
	int insert_conforming_loop(InputIterator begin, InputIterator end, Cell_handle hint = Cell_handle());
	template < class InputIterator >
	int insert_marked_loop(InputIterator begin, InputIterator end, Cell_handle hint = Cell_handle());
	template < class InputIterator >
	int remove_conforming_loop(InputIterator begin, InputIterator end, Cell_handle hint = Cell_handle());

protected:
	template < class InputIterator, class OutputIterator >
	void insert_loop_points(InputIterator begin, InputIterator end, OutputIterator out, Cell_handle hint = Cell_handle());

public:
	// REMOVE

	// Vertices that are incident to a conforming edge cannot be removed.
	bool remove(Vertex_handle v);

	// Remove a conforming edge or segment.
	bool remove_conforming(Cell_handle c, int li, int lj);
	bool remove_conforming(const Edge& e) {
		return remove_conforming(e.first, e.second, e.third);}
	bool remove_conforming(Vertex_handle va, Vertex_handle vb);
	bool remove_conforming(Bi_vertex c) {
		return remove_conforming(c.first, c.second);}

public:
	// QUERY CONFORMING

	// Check if an edge is conforming.
	bool is_conforming(Cell_handle c, int li, int lj) const {return c->is_conforming(li, lj);}
	bool is_conforming(const Edge& e) const {return is_conforming(e.first, e.second, e.third);}

	// Check if an edge is marked.
	bool is_marked(Cell_handle c, int li, int lj) const {return c->is_marked(li, lj);}
	bool is_marked(const Edge& e) const {return is_marked(e.first, e.second, e.third);}

	// Check if a vertex is incident to a conforming edge.
	bool are_there_incident_conforming(Vertex_handle v) const;

	// Check if the triangulation contains an edge ac on the segment ab, where c(li, lj) is that edge.
	bool includes_edge(Vertex_handle va, Vertex_handle vb, Vertex_handle& vc, Cell_handle& c, int& li, int& lj) const;

	// Check if the triangulation is valid.
	virtual bool is_valid(bool verbose = false, int level = 0) const;
	virtual bool is_valid(Cell_handle c, bool verbose = false, int level = 0) const;

private:
	// Collect the points encroaching upon the segment ac, where c is the first point on the ray ab.
	// encroaching must point to a Point.
	// Returns true if a conforming edge or constrained facet was encountered and the intersection inserted.
	template <class OutputIterator>
	bool collect_encroaching(Vertex_handle va, Vertex_handle vb, OutputIterator encroaching, Vertex_handle& vc);

	// Sort the points encroaching upon the segment ab to the from of the collection.
	// begin and end must point to Point.
	Points_iterator sort_encroaching(Points_iterator begin, Points_iterator end, const Point& a, const Point& b) const;
}; // Conforming_Delaunay_triangulation_3

template < class Gt, class Tds, class Itag >
typename Conforming_Delaunay_triangulation_3<Gt,Tds,Itag>::Vertex_handle
Conforming_Delaunay_triangulation_3<Gt,Tds,Itag>::
intersect(Cell_handle /* c */, Locate_type /* lt */, int /* li */, int /* lj */, Vertex_handle /* va */, Vertex_handle /* vb */, No_intersection_tag) {
	std::cerr << " sorry, this triangulation does not deal with" << std::endl
			  << " intersecting conforming edges" << std::endl;
	CGAL_triangulation_assertion(false);
	return Vertex_handle();
}

template < class Gt, class Tds, class Itag >
typename Conforming_Delaunay_triangulation_3<Gt,Tds,Itag>::Vertex_handle
Conforming_Delaunay_triangulation_3<Gt,Tds,Itag>::
intersect(Cell_handle c, Locate_type lt, int li, int lj, Vertex_handle va, Vertex_handle vb, Exact_intersections_tag) {
	// Compute the intersection of the segment ab and a facet or edge of cell c.
	// Note that no further check is made whether this intersection actually occurs inside the segment(s) or triangle.
	Point p;
	CGAL_triangulation_assertion_code( bool ok = )
		compute_intersection(lt, c, li, lj, va, vb, p);
	CGAL_triangulation_assertion(ok);
	Vertex_handle v = insert(p, c);
	v->steiner() = true;
	return v;
}

template < class Gt, class Tds, class Itag >
typename Conforming_Delaunay_triangulation_3<Gt,Tds,Itag>::Vertex_handle
Conforming_Delaunay_triangulation_3<Gt,Tds,Itag>::
intersect(Cell_handle c, Locate_type lt, int li, int lj, Vertex_handle va, Vertex_handle vb, Exact_predicates_tag) {
	// Compute the intersection of the segment ab and a facet or edge of cell c.
	// Note that no further check is made whether this intersection actually occurs inside the segment(s) or triangle.
	Point p;
	bool ok = compute_intersection(lt, c, li, lj, va, vb, p);
	if (ok) {
		Vertex_handle v = insert(p, c);
		v->steiner() = true;
		return v;
	}

	// Return the vertex closest to the intersection.
	Vertex_handle closest;
	switch (lt) {
		case DT::EDGE: {
			Line ab(va->point(), vb->point());
			Line ij(c->vertex(li)->point(), c->vertex(lj)->point());
			FT dist2 = squared_distance(va->point(), ij);
			FT db = squared_distance(vb->point(), ij);
			FT di = squared_distance(c->vertex(li)->point(), ab);
			FT dj = squared_distance(c->vertex(lj)->point(), ab);
			closest = va;
			if (db < dist2) {dist2 = db; closest = vb;}
			if (di < dist2) {dist2 = di; closest = c->vertex(li);}
			if (dj < dist2) {dist2 = dj; closest = c->vertex(lj);}
		}
		case DT::FACET: {
			Line ab(va->point(), vb->point());
			Plane p = plane(c, li);
			FT dist2 = squared_distance(va->point(), p);
			FT db = squared_distance(vb->point(), p);
			FT d1 = squared_distance(c->vertex((li+1)&3)->point(), ab);
			FT d2 = squared_distance(c->vertex((li+2)&3)->point(), ab);
			FT d3 = squared_distance(c->vertex((li+3)&3)->point(), ab);
			closest = va;
			if (db < dist2) {dist2 = db; closest = vb;}
			if (d1 < dist2) {dist2 = d1; closest = c->vertex((li+1)&3);}
			if (d2 < dist2) {dist2 = d2; closest = c->vertex((li+2)&3);}
			if (d3 < dist2) {dist2 = d3; closest = c->vertex((li+3)&3);}
		}
		default:
			CGAL_triangulation_assertion(false);
	}

	closest->steiner() = true;
	return closest;
}

template < class Gt, class Tds, class Itag >
typename Conforming_Delaunay_triangulation_3<Gt,Tds,Itag>::Point
Conforming_Delaunay_triangulation_3<Gt,Tds,Itag>::
construct_split_point(Vertex_handle va, Vertex_handle vb, const Point& ref) {
	// Blend of Si's and Shewchuk's segment recovery.
	// We are actually interested in constructing a layer of protecting Steiner points around each acute vertex,
	// such that no new Steiner point is inserted inside this layer.
	// Shewchuk and Si place the protecting Steiner points on a ball around the acute vertex.
	// However, determining if a vertex is acute is expensive. Additionally, the main purpose of the protecting points
	// is to place new points outside the diametral ball with their acute vertex.
	// This is always done if we place new Steiner points on the projection of the reference point onto the segment.
	// To guarantee a minimum edge length, we never place a Steiner point closer to an endpoint than to the reference.

#ifdef SPLIT_P_PROJ
	// Place the split point at the midpoint, unless the reference encroaches upon the resulting edges.
	// In that case, the split point is taken such that its projection onto ap or bp is p.
	Point s = midpoint(va->point(), vb->point());
	
	typename Geom_traits::Side_of_bounded_sphere_3 side_of_bounded_sphere = geom_traits().side_of_bounded_sphere_3_object();

	if (side_of_bounded_sphere(va->point(), mid, ref) == ON_BOUNDED_SIDE) {
		Object result = Gt().intersect_3_object()(Line(va->point(), vb->point()), Plane(ref, va->point() - ref));
		CGAL_triangulation_assertion(assign(s, result));
	}
	else if (side_of_bounded_sphere(vb->point(), mid, ref) == ON_BOUNDED_SIDE) {
		Object result = Gt().intersect_3_object()(Line(va->point(), vb->point()), Plane(ref, vb->point() - ref));
		CGAL_triangulation_assertion(assign(s, result));
	}

	return s;
#endif

	// Place the split point at the projection of the reference onto ab, unless that is too close to a or b, or it is not on the edge.
	// In that case, the split point is the midpoint between a and b.
	Line l(va->point(), vb->point());
	Point s = l.projection(ref);

	// The following will enforce a minimum edge length.
	if (has_larger_distance_to_point(s, ref, va->point()) || has_larger_distance_to_point(s, ref, vb->point()))
		s = midpoint(va->point(), vb->point());

	// Note that for inexact construction kernels, this assertion may fail for both projections and midpoints.
	//CGAL_triangulation_assertion(l.has_on(s));

	// This method doesn't actually produce a correct bounding ball, because the Steiner points may be inserted closer to the endpoint than the reference point..

	return s;

#ifdef SPLIT_SI
	// Si's point segment recovery.
	typedef Gt::Sphere_3	Sphere;
	Point s;
	Object result;
	if (c == Vertex_handle()) {
		// Rule 1.
		FT sd, r = squared_distance(a->point(), b->point());
		if (4*(sd = squared_distance(a->point(), p->point())) < r) {
			result = Gt().intersect_3_object()(Sphere(a->point(), sd), Segment(a->point(), b->point()));
		}
		else if (4*(sd = squared_distance(b->point(), p->point())) < r) {
			result = Gt().intersect_3_object()(Sphere(b->point(), sd), Segment(a->point(), b->point()));
		}
		else {
			Vertex_handle v = insert(midpoint(a->point(), b->point()), a->cell());
			return v;
		}
		if (!assign(s, result)) return Vertex_handle();
	}
	else {
		// Rule 2.
		FT r = squared_distance(c->point(), p->point());
		result = Gt().intersect_3_object()(Sphere(c->point(), r), Segment(a->point(), b->point()));
		if (!assign(s, result)) return Vertex_handle();

		FT sd = squared_distance(s, p->point());
		if (sd > squared_distance(s, b->point())) {
			//Rule 3.
			if (4*sd < squared_distance(s, a->point())) {
				result = Gt().intersect_3_object()(Sphere(s, sd), Segment(a->point(), b->point()));
				if (!assign(s, result))
					return Vertex_handle();
			}
			else {
				s = midpoint(a->point(), s);
			}
		}
	}

	return s;
#endif
}

template < class Gt, class Tds, class Itag >
template <bool make_gabriel>
void Conforming_Delaunay_triangulation_3<Gt,Tds,Itag>::
conform_segment(Vertex_handle va, Vertex_handle vb, Points_iterator encr_begin, Points_iterator encr_end) {
	CGAL_triangulation_precondition(va != vb);
	// This insertion method is an adjustment of the segment insertion
	// in the paper "Constrained Delaunay tetrahedral mesh generation and refinement" by Si (2010).
	// The points of the tetrahedra intersected by ab (encroaching points) are collected;
	// from these points, the point making the largest circumsphere with ab is the reference point.
	// The reference point indicates where ab is split and the subsegments each get a subset
	// of the encroaching points and then the process is repeated.

	// Set the initial part of the segment that is already in the triangulation to conforming.
	Vertex_handle vi;
	Cell_handle c;
	int li, lj;
	while (includes_edge(va, vb, vi, c, li, lj) && (!make_gabriel || this->is_Gabriel(c, li, lj) ) ) { 
		set_conforming(c, li, lj, true);
		va = vi;
		if (va == vb)
			return;
	}

	// If the edge is not in the triangulation, there must be encroaching points.
	// Collect the points encroaching upon and close to ab.
	if (encr_begin == encr_end) {
		Point_list encroaching;
		bool intersect = collect_encroaching(va, vb, std::back_inserter(encroaching), vi);

		if (!intersect && encroaching.empty()) {
			// There are occurrences where none of the vertices of the cells intersected by ab encroach upon ab.
			// In this case, we insert the midpoint of ab.
			vi = insert(midpoint(va->point(), vb->point()), va->cell());
			vi->steiner() = true;
			intersect = true;
		}

		// Insert ai, and if needed ib.
		conform_segment<make_gabriel>(va, vi, encroaching.begin(), encroaching.end());
		if (vi != vb) {
			conform_segment<make_gabriel>(vi, vb, encr_begin, encr_end);
		}
		return;
	}

	// Find the reference point.
	// This is the encroaching point closest to the segment.
	FT sr, radius2 = 0;
	Point ref;
	for (Points_iterator it = encr_begin; it != encr_end; ++it) {
		CGAL_triangulation_assertion(!collinear(va->point(), *it, vb->point()));
		sr = squared_radius(va->point(), *it, vb->point());
		if (radius2 < sr) {
			radius2 = sr;
			ref = *it;
		}
	}
	CGAL_triangulation_assertion(radius2 > 0);

	// Split the edge based on its reference point.
	Point s = construct_split_point(va, vb, ref);
	Vertex_handle vs = insert(s, va->cell());
	vs->steiner() = true;

	// Sort the encroaching points to encroaching upon as, sb, and only ab.
	Points_iterator after_as = sort_encroaching(encr_begin, encr_end, va->point(), vs->point());
	Points_iterator after_sb = sort_encroaching(after_as, encr_end, vs->point(), vb->point());

	// Insert the segments on both sides of s as constraints.
	conform_segment<make_gabriel>(va, vs, encr_begin, after_as);
	conform_segment<make_gabriel>(vs, vb, after_as, after_sb);
}

template < class Gt, class Tds, class Itag >
void Conforming_Delaunay_triangulation_3<Gt,Tds,Itag>::
mark_segment(Vertex_handle va, Vertex_handle vb) {
	Vertex_handle vi;
	Cell_handle c;
	int li, lj;
	while (va != vb) {
		// All segments must be in the triangulation.
		CGAL_triangulation_assertion_code(bool ok = )
			includes_edge(va, vb, vi, c, li, lj);
		CGAL_triangulation_assertion(ok);

		mark_edge(c, li, lj);

		va = vi;
	}
}

template < class Gt, class Tds, class Itag >
void Conforming_Delaunay_triangulation_3<Gt,Tds,Itag>::
insert_marked(Vertex_handle va, Vertex_handle vb) {
	// Insert the segment and then mark all parts of the segment.
	insert_conforming(va, vb);
	mark_segment(va, vb);
}

template < class Gt, class Tds, class Itag >
void Conforming_Delaunay_triangulation_3<Gt,Tds,Itag>::
set_conforming(Cell_handle c, int li, int lj, bool C) {
	// The cell is already in the triangulation, so the edge is already Delaunay.
	Vertex_handle vi = c->vertex(li), vj = c->vertex(lj);
	switch (dimension()) {
		case 3: {
			// Constrain the edge in all its incident cells.
			Cell_circulator	it = incident_cells(c, li, lj), start(it);
			do {
				it->set_conforming(it->index(vi), it->index(vj), C);
				++it;
			} while (it != start);
			break;
		}
		case 2: {
			// Constrain the edge in both cells that share it.
			Cell_handle n = c->neighbor(3-li-lj);
			n->set_conforming(n->index(vi), n->index(vj), C);
		}
		case 1: {
			// Constrain the edge in the cell itself.
			c->set_conforming(li, lj, C);
			break;
		}
	}
}

template < class Gt, class Tds, class Itag >
void Conforming_Delaunay_triangulation_3<Gt,Tds,Itag>::
mark_edge(Cell_handle c, int li, int lj) {
	// The cell is already in the triangulation, so the edge is already Delaunay.
	Vertex_handle vi = c->vertex(li), vj = c->vertex(lj);
	switch (dimension()) {
		case 3: {
			// Constrain the edge in all its incident cells.
			Cell_circulator	it = incident_cells(c, li, lj), start(it);
			do {
				it->mark(it->index(vi), it->index(vj));
				it++;
			} while (it != start);
			break;
		}
		case 2: {
			// Constrain the edge in both cells that share it.
			Cell_handle n = c->neighbor(3-li-lj);
			n->mark(n->index(vi), n->index(vj));
		}
		case 1: {
			// Constrain the edge in the cell itself.
			c->mark(li, lj);
			break;
		}
	}
}

template < class Gt, class Tds, class Itag >
template < class InputIterator >
int Conforming_Delaunay_triangulation_3<Gt,Tds,Itag>::
insert_conforming_loop(InputIterator begin, InputIterator end, Cell_handle hint) {
	// Insert the points in the loop.
	std::list<Bi_vertex> segments;
	insert_loop_points(begin, end, std::back_inserter(segments), hint);

	// Insert the segments in the loop.
	for (typename std::list<Bi_vertex>::iterator it = segments.begin(); it != segments.end(); ++it)
		insert_conforming(*it);

	return segments.size();
}

template < class Gt, class Tds, class Itag >
template < class InputIterator >
int Conforming_Delaunay_triangulation_3<Gt,Tds,Itag>::
insert_marked_loop(InputIterator begin, InputIterator end, Cell_handle hint) {
	// Insert the points in the loop.
	std::list<Bi_vertex> segments;
	insert_loop_points(begin, end, std::back_inserter(segments), hint);

	// Insert the segments in the loop.
	for (typename std::list<Bi_vertex>::iterator it = segments.begin(); it != segments.end(); ++it)
		insert_marked(*it);

	return segments.size();
}

template < class Gt, class Tds, class Itag >
template < class InputIterator >
int Conforming_Delaunay_triangulation_3<Gt,Tds,Itag>::
remove_conforming_loop(InputIterator begin, InputIterator end, Cell_handle hint) {
	// Insert the points in the loop.
	std::list<Bi_vertex> segments;
	insert_loop_points(begin, end, std::back_inserter(segments), hint);

	// Insert the segments in the loop.
	for (typename std::list<Bi_vertex>::iterator it = segments.begin(); it != segments.end(); ++it)
		remove_conforming(*it);

	return segments.size();
}

template < class Gt, class Tds, class Itag >
template < class InputIterator, class OutputIterator >
void Conforming_Delaunay_triangulation_3<Gt,Tds,Itag>::
insert_loop_points(InputIterator begin, InputIterator end, OutputIterator out, Cell_handle hint) {
	if (begin == end)
		return;

	// Insert the points in the loop.
	Vertex_handle prev = insert(*begin, hint), first(prev);
	hint = prev->cell();
	for (InputIterator it = ++begin; it != end; ++it) {
		Vertex_handle cur = insert(*it, hint);
		if (prev != cur)
			*out++ = Bi_vertex(prev, cur);
		hint = cur->cell();
		prev = cur;
	}
	if (first != prev)
		*out++ = Bi_vertex(prev, first);
}

template < class Gt, class Tds, class Itag >
void Conforming_Delaunay_triangulation_3<Gt,Tds,Itag>::
restore_conforming(Cell_handle c, int li) {
	// Restore the conforming edges of facet i of c, based on its neighboring cell.
	Cell_handle n = c->neighbor(li);
	Vertex_handle v1, v2; Edge e;
	for (int j = 0; j < dimension()+1; ++j) {
		if (j != li) {
			e = opposite_edge(c, li, (li+j)&3);
			v1 = c->vertex(e.second);
			v2 = c->vertex(e.third);
			if (n->is_marked(n->index(v1), n->index(v2)))
				e.first->mark(e.second, e.third);
			else
				e.first->set_conforming(e.second, e.third, n->is_conforming(n->index(v1), n->index(v2)));
		}
	}
}

template < class Gt, class Tds, class Itag >
typename Conforming_Delaunay_triangulation_3<Gt,Tds,Itag>::Vertex_handle
Conforming_Delaunay_triangulation_3<Gt,Tds,Itag>::
insert(const Point& p, Locate_type lt, Cell_handle c, int li, int lj) {
	switch (dimension()) {
		case 3: {
			typename DT::Conflict_tester_3 tester(p, this);
			Vertex_handle v = insert_in_conflict(p, lt, c, li, lj, tester, conforming_visitor);
			return v;
		} // dim 3
		case 2: {
			typename DT::Conflict_tester_2 tester(p, this);
			Vertex_handle v = insert_in_conflict(p, lt, c, li, lj, tester, conforming_visitor);
			return v;
		} // dim 2
		default :
			// dimension <= 1
			// Do not use the generic insert.
			return Tr::insert(p, c);
	}
}

template < class Gt, class Tds, class Itag >
void Conforming_Delaunay_triangulation_3<Gt,Tds,Itag>::
insert_conforming(const Point& a, const Point& b, Cell_handle hint) {
	// Insert the two points and then conform the segment between them.
	Vertex_handle va = insert(a, hint);
	Vertex_handle vb = insert(b, va->cell());
	if (va != vb)
		insert_conforming(va, vb);
}

template < class Gt, class Tds, class Itag >
template < class InputIterator >
int Conforming_Delaunay_triangulation_3<Gt,Tds,Itag>::
insert_conforming(InputIterator begin, InputIterator end, Cell_handle hint) {
	// Insert the points.
	std::list<Bi_vertex> segments;
	for (InputIterator it = begin; it != end; ++it) {
		Vertex_handle va = insert(it->first, hint);
		hint = va->cell();
		if (it->first != it->second) {
			Vertex_handle vb = insert(it->second, hint);
			hint = vb->cell();
			segments.push_back(Bi_vertex(va, vb));
		}
	}

	// Insert the segments
	for (typename std::list<Bi_vertex>::iterator it = segments.begin(); it != segments.end(); ++it)
		insert_conforming(*it);

	return segments.size();
}

template < class Gt, class Tds, class Itag >
bool Conforming_Delaunay_triangulation_3<Gt,Tds,Itag>::
remove(Vertex_handle v) {
	// The vertex cannot be removed if it is incident to a conforming edge.
	if (are_there_incident_conforming(v))
		return false;

	// Collect the opposite facets as seen from the outside.
	std::vector<Cell_handle> cells; cells.reserve(32);
	incident_cells(v, std::back_inserter(cells));
	std::vector<Facet> boundary; boundary.reserve(cells.size());

	for (typename std::vector<Cell_handle>::const_iterator it = cells.begin(); it != cells.end(); ++it) {
		Cell_handle n = (*it)->neighbor((*it)->index(v));
		boundary.push_back(Facet(n, n->index(*it)));
	}
	cells.clear();

	// Remove the vertex.
#ifdef CGAL_DELAUNAY_3_OLD_REMOVE
	if (dimension() == 3 && !test_dim_down(v)) {
		remove_3D_ear(v);
	} else {
#endif
	cDT tmp;
	typename DT::Vertex_remover remover(tmp);
	Tr::remove(v, remover);
#ifdef CGAL_DELAUNAY_3_OLD_REMOVE
	}
#endif

	// Reinsert any conforming edges on the boundary of the retriangulated region.
	for (typename std::vector<Facet>::iterator it = boundary.begin(); it != boundary.end(); ++it) {
		restore_conforming(mirror_facet(*it));
	}

	CGAL_triangulation_expensive_postcondition(is_valid());
	return true;
}

template < class Gt, class Tds, class Itag >
bool Conforming_Delaunay_triangulation_3<Gt,Tds,Itag>::
remove_conforming(Cell_handle c, int li, int lj) {
	// Mark the edge as non-conforming
	set_conforming(c, li, lj, false);
	return true;
}

template < class Gt, class Tds, class Itag >
bool Conforming_Delaunay_triangulation_3<Gt,Tds,Itag>::
remove_conforming(Vertex_handle va, Vertex_handle vb) {
	CGAL_triangulation_precondition(va != vb);
	// All the parts of the segment ab that are constrained and
	// not incident to a constrained facet are marked unconstrained.
	// The return is true if the whole segment ab is covered by unconstrained edges;
	// otherwise, returns false;

	// Check if first part of the segment is an edge in the triangulation.
	Vertex_handle vi;
	Cell_handle c;
	int li, lj;
	bool complete = true;
	if (includes_edge(va, vb, vi, c, li, lj)) {
		// Remove the first part.
		if (!remove_conforming(c, li, lj))
			complete = false;
		if (vi != vb)  {
			// Remove the next part.
			if (!remove_conforming(vi, vb))
				complete = false;
		}
		return complete;
	}

	return false;
}

template < class Gt, class Tds, class Itag >
bool Conforming_Delaunay_triangulation_3<Gt,Tds,Itag>::
are_there_incident_conforming(Vertex_handle v) const {
	if (dimension() > 0) {
		std::vector<Edge> edges;
		edges.reserve(32);
		finite_incident_edges(v, std::back_inserter(edges));
		for (typename std::vector<Edge>::iterator it = edges.begin(); it != edges.end(); ++it)
			if (is_conforming(*it))
				return true;
	}

	return false;
}

template < class Gt, class Tds, class Itag >
bool Conforming_Delaunay_triangulation_3<Gt,Tds,Itag>::
includes_edge(Vertex_handle va, Vertex_handle vb, Vertex_handle& vc, Cell_handle& c, int& li, int& lj) const {
	CGAL_triangulation_precondition(!is_infinite(va) && !is_infinite(vb));

	// Returns true if the line segment ab contains an edge e 
	// incident to a, false otherwise.
	// If true, vc becomes the vertex of e distinct from a,
	// c is a cell incident to e where e = (c, li, lj).
	Vertex_handle v;
	std::vector<Edge> edges; edges.reserve(64);
	finite_incident_edges(va, std::back_inserter(edges));
	for (typename std::vector<Edge>::const_iterator it = edges.begin(); it != edges.end(); ++it) {
		// Find the other vertex of the edge.
		v = it->first->vertex(it->second);
		if (v == va) v = it->first->vertex(it->third);

		if (is_infinite(v))
			continue;

		// Check if this edge is (part of) ab.
		if (v == vb) {
			vc = vb;
			c = it->first;
			li = it->second;
			lj = it->third;
			return true;
		}
		else if (sign((va->point() - v->point()) * (vb->point() - v->point())) == NEGATIVE &&
				 collinear(va->point(), v->point(), vb->point())/* &&
				 collinear_are_ordered_along_line(va->point(), v->point(), vb->point())*/) {
			vc = v;
			c = it->first;
			li = it->second;
			lj = it->third;
			return true;
		}
	}
	return false;
}

template < class Gt, class Tds, class Itag >
bool Conforming_Delaunay_triangulation_3<Gt,Tds,Itag>::
is_valid(bool verbose, int level) const {
	if (!tds().is_valid(verbose,level)) {
		if (verbose)
			std::cerr << "invalid data structure" << std::endl;
		CGAL_triangulation_assertion(false);
		return false;
	}

	if (infinite_vertex() == Vertex_handle()) {
		if (verbose)
			std::cerr << "no infinite vertex" << std::endl;
		CGAL_triangulation_assertion(false);
		return false;
	}

	int inf;
	Cell_handle n;
	Vertex_handle v1, v2;
	switch (dimension()) {
		case 3: {
			for (All_cells_iterator it = this->all_cells_begin(); it != this->all_cells_end(); ++it) {
				for (int i = 0; i < 4; ++i) {
					if (!it->neighbor(i)->has_neighbor(it)) {
						if (verbose)
							std::cerr << "inconsistent neighbors " << std::endl;
						CGAL_triangulation_assertion(false);
						return false;
					}
				}
				if (it->has_vertex(infinite_vertex(), inf)) {
					for (int i = 0; i < 4; ++i) {
						if (i != inf) {
							if (is_conforming(it, i, inf)) {
								if (verbose)
									std::cerr << "conforming infinite edge " << std::endl;
								CGAL_triangulation_assertion(false);
								return false;
							}
						}
					}
				}
				else {
					is_valid_finite(it);
					for (int i = 0; i < 4; ++i) {
						n = it->neighbor(i);
						if (!is_infinite(n->vertex(n->index(it)))) {
							if (this->side_of_sphere(it, n->vertex(n->index(it))->point()) == ON_BOUNDED_SIDE) {
								if (verbose)
									std::cerr << "non-empty sphere " << std::endl;
								CGAL_triangulation_assertion(false);
								return false;
							}
						}
					}
				}
			}
			for (Finite_edges_iterator it = finite_edges_begin(); it != finite_edges_end(); ++it) {
				v1 = it->first->vertex(it->second);
				v2 = it->first->vertex(it->third);
				bool conforming = is_conforming(*it);

				Cell_circulator	cit = incident_cells(it->first, it->second, it->third, it->first), start(cit);
				do {
					if (cit->is_marked(cit->index(v1), cit->index(v2))) {
						if (verbose)
							std::cerr << "marked edge " << std::endl;
						CGAL_triangulation_assertion(false);
						return false;
					}

					if (cit->is_conforming(cit->index(v1), cit->index(v2)) != conforming) {
						if (verbose)
							std::cerr << "inconsistent conforming edge " << std::endl;
						CGAL_triangulation_assertion(false);
						return false;
					}
					cit++;
				} while (cit != start);
			}
			break;
		}
		case 2: {
			for (typename DT::All_facets_iterator it = this->all_facets_begin(); it != this->all_facets_end(); ++it) {
				for (int i = 0; i < 3; ++i) {
					if (!it->first->neighbor(i)->has_neighbor(it->first)) {
						if (verbose)
							std::cerr << "inconsistent neighbors " << std::endl;
						CGAL_triangulation_assertion(false);
						return false;
					}
				}
				if (it->first->has_vertex(infinite_vertex(), inf) && inf < 3) {
					for (int i = 0; i < 3; ++i) {
						if (i != inf) {
							if (is_conforming(it->first, i, 3)) {
								if (verbose)
									std::cerr << "conforming infinite edge " << std::endl;
								CGAL_triangulation_assertion(false);
								return false;
							}
						}
					}
				}
				else {
					is_valid_finite(it->first);
					for (int i = 0; i < 3; ++i) {
						n = it->first->neighbor(i);
						if(!is_infinite(n->vertex(n->index(it->first)))) {
							if (this->side_of_circle(it->first, 3, n->vertex(n->index(it->first))->point()) == ON_BOUNDED_SIDE) {
								if (verbose)
									std::cerr << "non-empty circle " << std::endl;
								CGAL_triangulation_assertion(false);
								return false;
							}
						}

						Edge o = this->opposite_edge(it->first, i, 3);
						if (o.first->is_marked(o.second, o.third)) {
							if (verbose)
								std::cerr << "marked edge " << std::endl;
							CGAL_triangulation_assertion(false);
							return false;
						}

						if (is_conforming(o) != 
							is_conforming(n, n->index(o.first->vertex(o.second)), n->index(o.first->vertex(o.third)))) {
							if (verbose)
								std::cerr << "inconsistent conforming edge " << std::endl;
							CGAL_triangulation_assertion(false);
							return false;
						}
					}
				}
			}
			break;
		}
		case 1: {
			for (Finite_edges_iterator it = finite_edges_begin(); it != finite_edges_end(); ++it)
				is_valid_finite(it->first);
			break;
		}
	}
	if (verbose)
		std::cerr << "valid conforming Delaunay triangulation" << std::endl;
	return true;
}

template < class Gt, class Tds, class Itag >
bool Conforming_Delaunay_triangulation_3<Gt,Tds,Itag>::
is_valid(Cell_handle c, bool verbose, int level) const {
	if (!DT::is_valid(c, verbose, level)) {
		CGAL_triangulation_assertion(false);
		return false;
	}

	switch (dimension()) {
		case 3: {
			for (int i = 0; i < 4; ++i) {
				for (int j = i+1; j < 4; ++j) {
					Vertex_handle v1 = c->vertex(i), v2 = c->vertex(j);
					if (is_infinite(v1) || is_infinite(v2)) {
						if (is_conforming(c, i, j)) {
							if (verbose)
								std::cerr << "conforming infinite edge " << std::endl;
							CGAL_triangulation_assertion(false);
							return false;
						}
					}
					else {
						bool conforming = is_conforming(c, i, j);

						Cell_circulator	cit = incident_cells(c, i, j, c), start(cit);
						do {
							if (cit->is_marked(cit->index(v1), cit->index(v2))) {
								if (verbose)
									std::cerr << "marked edge " << std::endl;
								CGAL_triangulation_assertion(false);
								return false;
							}

							if (cit->is_conforming(cit->index(v1), cit->index(v2)) != conforming) {
								if (verbose)
									std::cerr << "inconsistent conforming edge " << std::endl;
								CGAL_triangulation_assertion(false);
								return false;
							}
							cit++;
						} while (cit != start);
					}
				}
			}
			break;
		}
		case 2: {
			for (int v=0; v<3; ++v) {
				int i = (v+1)%3,
					j = (v+2)%3;

				Vertex_handle v1 = c->vertex(i), v2 = c->vertex(j);
				if (is_infinite(v1) || is_infinite(v2)) {
					if (is_conforming(c, i, j)) {
						if (verbose)
							std::cerr << "conforming infinite edge " << std::endl;
						CGAL_triangulation_assertion(false);
						return false;
					}
				}
				else {
					if (c->is_marked(i, j)) {
						if (verbose)
							std::cerr << "marked edge " << std::endl;
						CGAL_triangulation_assertion(false);
						return false;
					}

					Cell_handle n = c->neighbor(v);
					if (c->is_conforming(i, j) != n->is_conforming(n->index(v1), n->index(v2))) {
						if (verbose)
							std::cerr << "inconsistent conforming edge " << std::endl;
						CGAL_triangulation_assertion(false);
						return false;
					}
				}
			}
			break;
		}
	}
	if (verbose)
		std::cerr << "valid Delaunay cell" << std::endl;
	return true;
}

template < class Gt, class Tds, class Itag >
template <class OutputIterator>
bool Conforming_Delaunay_triangulation_3<Gt,Tds,Itag>::
collect_encroaching(Vertex_handle va, Vertex_handle vb, OutputIterator encroaching, Vertex_handle& vc) {
	CGAL_triangulation_precondition(va != vb);
	CGAL_triangulation_precondition(!is_infinite(va) && !is_infinite(vb));
	CGAL_triangulation_precondition(dimension() >= 2);

	// Fill encroaching with the points that are encroaching upon the segment ab.
	// vc is set to the first vertex on ab starting from a;
	// if ab intersects an edge or facet the intersection is inserted and true returned.

	// This is done by traversing cells from a towards b until either b is found,
	// or any other collinear vertex, or an intersecting constraint.
	Locate_type lt; int li, lj;
	Encroaching_collecter cit(va, vb->point(), this, encroaching);
	while (cit.has_next()) {
		// Walk to the next cell.
		++cit;

		cit.traversed(lt, li, lj);

		// Stop when either a new vertex is reached, or a barrier (a conforming edge) is encountered.
		if (lt == DT::VERTEX) {
			vc = cit->vertex(li);
			if (va == vc)
				continue;
			return false;
		}

		if (cit.barrier_hit()) {
			vc = intersect(cit, lt, li, lj, va, vb);
			return true;
		}
	}

	// The cell must contain the target.
	CGAL_triangulation_assertion(cit->has_vertex(vb));
	vc = vb;
	return false;
}

template < class Gt, class Tds, class Itag >
typename Conforming_Delaunay_triangulation_3<Gt,Tds,Itag>::Points_iterator
Conforming_Delaunay_triangulation_3<Gt,Tds,Itag>::
sort_encroaching(Points_iterator begin, Points_iterator end, const Point& a, const Point& b) const {
	// Sorts the points encroaching upon ab to the front of the collection.
	// Returns an iterator to the first point not encroaching upon ab.
	Point temp;
	Points_iterator middle = end;
	for (Points_iterator it = begin; it != middle;) {
		if (collinear(a, b, *it) || !is_encroaching(a, b, *it)) {
			if (it != --middle) {
				// Swap the current point to the back.
				temp = *middle;
				*middle = *it;
				*it = temp;
			}
		}
		else
			++it;
	}
	return middle;
}

} //end of CGAL namespace

#endif // CGAL_CONFORMING_DELAUNAY_TRIANGULATION_3_H