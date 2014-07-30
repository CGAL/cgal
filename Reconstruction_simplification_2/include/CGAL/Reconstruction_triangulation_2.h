// Copyright (c) 2014  INRIA Sophia-Antipolis (France), INRIA Lorraine LORIA.
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
// Author(s)     : Fernando de Goes, Pierre Alliez, Ivo Vigan

#ifndef RECONSTRUCTION_TRIANGULATION_2_H
#define RECONSTRUCTION_TRIANGULATION_2_H

// STL
#include <map>
#include <set>
#include <list>
#include <queue>
#include <iostream>
#include <limits>

// CGAL
#include <CGAL/basic.h>
#include <CGAL/intersections.h>
#include <CGAL/Delaunay_triangulation_2.h>

// local
#include <CGAL/Sample.h>
#include <CGAL/Reconstruction_edge_2.h>
#include <CGAL/Cost.h>
#include <CGAL/Reconstruction_vertex_base_2.h>
#include <CGAL/Reconstruction_face_base_2.h>


// boost
#include <boost/multi_index/mem_fun.hpp>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/identity.hpp>
#include <boost/multi_index/member.hpp>

#define EPS   1e-15

namespace CGAL {


/// \internal
///  The Reconstruction_triangulation_2 class
///  provides the reconstruction simplex as well as the transportation plan.
/// - Each vertex stores a normal vector.
/// - A vertex a Sample which got assigned to it by the transportation plan,
///   well as the corresponding relocated Point (of type Kernel::Point_2).
/// - In order to solve a linear system over the triangulation, a vertex may be constrained
///   or not (i.e. may contribute to the right or left member of the linear system),
///   and has a unique index.
/// The vertex class must derive from Reconstruction_vertex_base_3.
///
///  @param Kernel   The underlying Kernel
///  @param Tds      Model of TriangulationDataStructure_2.
///  The vertex class must derive from Reconstruction_vertex_base_2.
///  The face   class must derive from Reconstruction_face_base_2.
template<class Kernel,
		class Tds_ = Triangulation_data_structure_2<
				Reconstruction_vertex_base_2<Kernel>,
				Reconstruction_face_base_2<Kernel> > >
class Reconstruction_triangulation_2: public Delaunay_triangulation_2<Kernel,
		Tds_> {
public:

    typedef Delaunay_triangulation_2<Kernel, Tds_> Base;

	typedef typename Kernel::FT FT;
	typedef typename Kernel::Point_2 Point;
	typedef typename Kernel::Vector_2 Vector;
	typedef typename Kernel::Ray_2 Ray;
	typedef typename Kernel::Line_2 Line;
	typedef typename Kernel::Segment_2 Segment;
	typedef typename Kernel::Triangle_2 Triangle;

	typedef typename Base::Vertex Vertex;
	typedef typename Base::Vertex_handle Vertex_handle;
	typedef typename Base::Vertex_iterator Vertex_iterator;
	typedef typename Base::Vertex_circulator Vertex_circulator;
	typedef typename Base::Finite_vertices_iterator Finite_vertices_iterator;

	typedef typename Base::Edge Edge;
	typedef typename Base::Edge_iterator Edge_iterator;
	typedef typename Base::Edge_circulator Edge_circulator;
	typedef typename Base::Finite_edges_iterator Finite_edges_iterator;

	typedef typename Base::Face Face;
	typedef typename Base::Face_handle Face_handle;
	typedef typename Base::Face_iterator Face_iterator;
	typedef typename Base::Face_circulator Face_circulator;
	typedef typename Base::Finite_faces_iterator Finite_faces_iterator;

	typedef std::map<Vertex_handle, Vertex_handle,
			less_Vertex_handle<Vertex_handle> > Vertex_handle_map;
	typedef std::map<Face_handle, Face_handle, less_Face_handle<Face_handle> >
			Face_handle_map;

	typedef std::set<Vertex_handle, less_Vertex_handle<Vertex_handle> >
			Vertex_handle_set;
	typedef std::set<Edge, less_Edge<Edge> > Edge_set;

	typedef std::list<Edge> Edge_list;

	typedef std::list<Point> Point_list;
	typedef typename Point_list::const_iterator Point_list_const_iterator;

	typedef Cost<FT> Cost;
	typedef Sample<Kernel> Sample;
	typedef std::list<Sample*> Sample_list;
	typedef typename Sample_list::const_iterator Sample_list_const_iterator;

	typedef Sample_with_priority<Sample> PSample;
	typedef std::priority_queue<PSample, std::vector<PSample>,
			greater_priority<PSample> > SQueue;

	typedef Reconstruction_edge_2<FT, Edge, Vertex_handle, Face_handle>
	Reconstruction_edge_2;

	 typedef boost::multi_index_container<
		  Reconstruction_edge_2,
		  boost::multi_index::indexed_by<
		  // sort by Reconstruction_edge_2::operator<
		  boost::multi_index::ordered_unique< boost::multi_index::identity<
		  	  Reconstruction_edge_2 > > ,
		  	// sort by Reconstruction_edge_2::priority()
		  boost::multi_index::ordered_non_unique<
		  	  boost::multi_index::const_mem_fun<
		  	  	  Reconstruction_edge_2,const FT,&Reconstruction_edge_2::priority> >
	  	  >
	  > MultiIndex;

	double m_factor; // ghost vs solid


public:
	Reconstruction_triangulation_2() {
		m_factor = 1.0;
	}

	~Reconstruction_triangulation_2() {
	}

	double& ghost_factor() {
		return m_factor;
	}

	const double& ghost_factor() const {
		return m_factor;
	}

	Edge random_finite_edge() {
		int nbf = Base::number_of_faces();
		int offset = random_int(0, nbf - 1);
		Finite_faces_iterator fi = Base::finite_faces_begin();
		for (int i = 0; i < offset; i++)
			fi++;
		Face_handle face = fi;
		int index = random_int(0, 40) % 3;
		return Edge(face, index);
	}

	// ACCESS //

	Vertex_handle source_vertex(const Edge& edge) const {
		return edge.first->vertex(Base::ccw(edge.second));
	}

	Vertex_handle target_vertex(const Edge& edge) const {
		return edge.first->vertex(Base::cw(edge.second));
	}

	Vertex_handle opposite_vertex(const Edge& edge) const {
		return edge.first->vertex(edge.second);
	}

	bool is_pinned(const Edge& edge) const {
		Vertex_handle s = source_vertex(edge);
		if (s->pinned())
			return true;
		return false;
	}

	Edge twin_edge(const Edge& edge) const {
		Face_handle f = edge.first;
		Vertex_handle v = source_vertex(edge);
		Face_handle nf = f->neighbor(edge.second);
		return Edge(nf, Base::ccw(nf->index(v)));
	}

	Edge next_edge(const Edge& edge) const {
		Face_handle f = edge.first;
		int index = Base::ccw(edge.second);
		return Edge(f, index);
	}

	Edge prev_edge(const Edge& edge) const {
		Face_handle f = edge.first;
		int index = Base::cw(edge.second);
		return Edge(f, index);
	}

	FT get_length(const Edge& edge) const {
		Segment segment = get_segment(edge);
		return std::sqrt(segment.squared_length());
	}

	Segment get_segment(const Edge& edge) const {
		const Point& ps = source_vertex(edge)->point();
		const Point& pt = target_vertex(edge)->point();
		return Segment(ps, pt);
	}

	Triangle get_triangle(Face_handle face) const {
		Vertex_handle v0 = face->vertex(0);
		Vertex_handle v1 = face->vertex(1);
		Vertex_handle v2 = face->vertex(2);
		return Triangle(v0->point(), v1->point(), v2->point());
	}

	// GET LINK //

	void get_vertices_from_edge_link(const Edge& edge,
			Vertex_handle_set& vertices) const {
		vertices.insert(opposite_vertex(edge));
		vertices.insert(opposite_vertex(twin_edge(edge)));
	}

	void get_vertices_from_vertex_link(Vertex_handle vertex,
			Vertex_handle_set& vertices) const {
		Vertex_circulator vcirc = Base::incident_vertices(vertex);
		Vertex_circulator vend = vcirc;
		CGAL_For_all(vcirc, vend)
		{
			Vertex_handle v = vcirc;
			vertices.insert(v);
		}
	}

	// boundary of star(vertex)
	// 'outward' chooses the orientation of the boundary
	void get_edges_from_star_minus_link(Vertex_handle vertex, Edge_list& hull,
			bool outward = false) const {
		Face_circulator fcirc = Base::incident_faces(vertex);
		Face_circulator fend = fcirc;
		CGAL_For_all(fcirc, fend)
		{
			Face_handle face = fcirc;
			int index = face->index(vertex);
			Edge edge(face, index);
			if (outward)
				edge = twin_edge(edge);
			hull.push_back(edge);
		}
	}

	// ATTRIBUTES //

	bool is_ghost(const Edge& edge) const {
		return edge.first->ghost(edge.second);
	}

	int get_plan(const Edge& edge) const {
		return edge.first->plan(edge.second);
	}

	void set_plan(const Edge& edge, int simplex) {
		edge.first->plan(edge.second) = simplex;
	}

	FT get_mass(const Edge& edge) const {
		return edge.first->mass(edge.second);
	}

	void set_mass(const Edge& edge, const FT mass) {
		edge.first->mass(edge.second) = mass;
	}

	const Cost& get_cost(const Edge& edge) const {
		return edge.first->cost(edge.second);
	}

	void set_vertex_cost(const Edge& edge, const Cost& cost) {
		edge.first->vertex_cost(edge.second) = cost;
	}

	void set_edge_cost(const Edge& edge, const Cost& cost) {
		edge.first->edge_cost(edge.second) = cost;
	}

	FT get_vertex_minus_edge_cost(const Edge& edge) const {
		const Cost& vcost = edge.first->vertex_cost(edge.second);
		const Cost& ecost = edge.first->edge_cost(edge.second);
		return vcost.finalize() - m_factor * ecost.finalize();
	}

	FT get_vertex_over_edge_cost(const Edge& edge) const {
		FT vvalue = edge.first->vertex_cost(edge.second).finalize();
		FT evalue = edge.first->edge_cost(edge.second).finalize();
		if (evalue == vvalue)
			return 1.0 / m_factor;
		return vvalue / (m_factor * evalue);
	}

	FT get_edge_relevance(const Edge& edge) const {
		FT M = get_mass(edge);
		if (M == 0.0)
			return 0.0;

		FT L = get_length(edge);
		FT cost = get_cost(edge).finalize();
		return M * L * L / cost;
	}

	FT get_density(const Edge& edge) const {
		FT length = get_length(edge);
		FT mass = get_mass(edge);
		return (mass / length);
	}

	unsigned int nb_samples(const Edge& edge) const {
		Edge twin = twin_edge(edge);
		return edge.first->samples(edge.second).size()
				+ twin.first->samples(twin.second).size();
	}

	void collect_samples_from_edge(const Edge& edge, Sample_list& samples) {
		const Sample_list& edge_samples = edge.first->samples(edge.second);
		samples.insert(samples.end(), edge_samples.begin(), edge_samples.end());
	}

	void collect_samples_from_vertex(Vertex_handle vertex, Sample_list& samples,
			bool cleanup) {
		Face_circulator fcirc = Base::incident_faces(vertex);
		Face_circulator fend = fcirc;
		CGAL_For_all(fcirc, fend)
		{
			Face_handle face = fcirc;
			int index = face->index(vertex);

			Edge edge(face, index);
			collect_samples_from_edge(edge, samples);

			Edge next = next_edge(edge);
			collect_samples_from_edge(next, samples);

			Edge prev = prev_edge(edge);
			collect_samples_from_edge(prev, samples);

			if (cleanup)
				face->clean_all_samples();
		}
		Sample* sample = vertex->get_sample();
		if (sample)
			samples.push_back(sample);
		if (cleanup)
			vertex->set_sample(NULL);
	}

	void collect_all_samples(Sample_list& samples) {
		for (Finite_edges_iterator ei = Base::finite_edges_begin();
				ei != Base::finite_edges_end(); ++ei) {
			Edge edge = *ei;
			Edge twin = twin_edge(edge);
			collect_samples_from_edge(edge, samples);
			collect_samples_from_edge(twin, samples);
		}
		for (Finite_vertices_iterator vi = Base::finite_vertices_begin();
				vi != Base::finite_vertices_end(); ++vi) {
			Vertex_handle v = vi;
			Sample* sample = v->get_sample();
			if (sample)
				samples.push_back(sample);
		}
	}

	void cleanup_assignments() {
		for (Finite_faces_iterator fi = Base::finite_faces_begin();
				fi != Base::finite_faces_end(); ++fi) {
			fi->clean_all_samples();
		}
		for (Finite_vertices_iterator vi = Base::finite_vertices_begin();
				vi != Base::finite_vertices_end(); ++vi) {
			vi->set_sample(NULL);
		}
	}

	// COST //

	Cost compute_total_cost() const {
		Cost sum;
		for (Finite_edges_iterator ei = Base::finite_edges_begin();
				ei != Base::finite_edges_end(); ++ei) {
			Edge edge = *ei;
			const Cost& cost = get_cost(edge);
			sum.update_max(cost);
			sum.add(cost);
		}
		return sum;
	}

	Cost compute_cost_around_vertex(Vertex_handle vertex) const {
		Cost inner;
		Cost outer;
		Face_circulator fcirc = Base::incident_faces(vertex);
		Face_circulator fend = fcirc;
		CGAL_For_all(fcirc, fend)
		{
			Face_handle face = fcirc;
			int index = face->index(vertex);

			Edge edge(face, index);
			Cost cost = get_cost(edge);
			outer.update_max(cost);
			outer.add(cost);

			edge = next_edge(edge);
			cost = get_cost(edge);
			inner.update_max(cost);
			inner.add(cost);

			edge = next_edge(edge);
			cost = get_cost(edge);
			inner.update_max(cost);
			inner.add(cost);
		}
		inner.divide(2.0);

		Cost sum;
		sum.add(inner);
		sum.add(outer);
		sum.update_max(inner);
		sum.update_max(outer);
		return sum;
	}

	void reset_all_costs() {
		for (Finite_edges_iterator ei = Base::finite_edges_begin();
				ei != Base::finite_edges_end(); ++ei) {
			Edge edge = *ei;
			update_cost(edge);
		}
	}

	void update_cost(const Edge& edge) {
		compute_mass(edge);
		compute_edge_cost(edge);
		compute_vertex_cost(edge);
		select_plan(edge);
	}

	void compute_mass(const Edge& edge) {
		FT mass = 0.0;

		typename Sample_list::const_iterator it;
		const Sample_list& samples0 = edge.first->samples(edge.second);
		for (it = samples0.begin(); it != samples0.end(); ++it) {
			Sample* sample = *it;
			mass += sample->mass();
		}

		Edge twin = twin_edge(edge);
		const Sample_list& samples1 = twin.first->samples(twin.second);
		for (it = samples1.begin(); it != samples1.end(); ++it) {
			Sample* sample = *it;
			mass += sample->mass();
		}

		set_mass(edge, mass);
		set_mass(twin, mass);
	}

	void select_plan(const Edge& edge) {
		// transport plan:
		// 0 - to vertex
		// 1 - to edge

		int plan = 0;
		FT diff = get_vertex_minus_edge_cost(edge);
		if (diff >= 0.0)
			plan = 1;

		Edge twin = twin_edge(edge);
		set_plan(edge, plan);
		set_plan(twin, plan);
	}

	void compute_edge_cost(const Edge& edge) {
		SQueue squeue;
		FT M = get_mass(edge);
		FT L = get_length(edge);
		sort_samples_from_edge(edge, squeue);
		Cost cost = compute_cost_from_squeue(squeue, M, L);

		Edge twin = twin_edge(edge);
		set_edge_cost(edge, cost);
		set_edge_cost(twin, cost);
	}

	void sort_samples_from_edge(const Edge& edge, SQueue& squeue) {
		typename Sample_list::const_iterator it;
		const Sample_list& samples0 = edge.first->samples(edge.second);
		for (it = samples0.begin(); it != samples0.end(); ++it) {
			Sample* sample = *it;
			squeue.push(PSample(sample, sample->coordinate()));
		}

		Edge twin = twin_edge(edge);
		const Sample_list& samples1 = twin.first->samples(twin.second);
		for (it = samples1.begin(); it != samples1.end(); ++it) {
			Sample* sample = *it;
			squeue.push(PSample(sample, 1.0 - sample->coordinate()));
		}
	}

	Cost compute_cost_from_squeue(SQueue& squeue, const FT M, const FT L) {
		if (squeue.empty())
			return Cost();
		if (M == 0.0)
			return Cost();

		Cost sum;
		FT start = 0.0;
		FT coef = L / M;
		while (!squeue.empty()) {
			PSample psample = squeue.top();
			squeue.pop();

			FT mass = psample.sample()->mass();
			FT coord = psample.priority() * L;
			FT bin = mass * coef;
			FT center = start + 0.5 * bin;
			FT pos = coord - center;

			FT norm2 = psample.sample()->distance2();
			FT tang2 = bin * bin / 12 + pos * pos;

			sum.add(Cost(norm2, tang2), mass);
			sum.compute_max(norm2, tang2);

			start += bin;
		}
		return sum;
	}

	void compute_vertex_cost(const Edge& edge) {
		Edge twin = twin_edge(edge);
		const Point& ps = source_vertex(edge)->point();
		const Point& pt = target_vertex(edge)->point();

		Sample_list samples;
		collect_samples_from_edge(edge, samples);
		collect_samples_from_edge(twin, samples);

		Cost sum;
		for (Sample_list_const_iterator it = samples.begin();
				it != samples.end(); ++it) {
			Sample* sample = *it;
			FT mass = sample->mass();
			const Point& query = sample->point();

			FT Ds = CGAL::squared_distance(query, ps);
			FT Dt = CGAL::squared_distance(query, pt);
			FT dist2 = ((std::min))(Ds, Dt);

			FT norm2 = sample->distance2();
			FT tang2 = dist2 - norm2;

			sum.add(Cost(norm2, tang2), mass);
			sum.compute_max(norm2, tang2);
		}
		set_vertex_cost(edge, sum);
		set_vertex_cost(twin, sum);
	}

	// SAMPLE //

	template<class Iterator> // value_type = Sample*
	void assign_samples(Iterator begin, Iterator end) {
		for (Iterator it = begin; it != end; ++it) {
			Sample* sample = *it;
			assign_sample(sample);
		}
	}

	template<class Iterator> // value_type = Sample*
	void assign_samples_brute_force(Iterator begin, Iterator end) {
		for (Iterator it = begin; it != end; ++it) {
			Sample* sample = *it;
			assign_sample_brute_force(sample);
		}
	}

	bool assign_sample(Sample* sample) {
		const Point& point = sample->point();
		Face_handle face = Base::locate(point);

		if (face == Face_handle() || Base::is_infinite(face)) {
			std::cout << "free bird" << std::endl;
			return false;
		}

		Vertex_handle vertex = find_nearest_vertex(point, face);
		if (vertex != Vertex_handle()) {
			assign_sample_to_vertex(sample, vertex);
			return true;
		}

		Edge edge = find_nearest_edge(point, face);
		assign_sample_to_edge(sample, edge);
		return true;
	}

	bool assign_sample_brute_force(Sample* sample) {
		const Point& point = sample->point();
		Face_handle nearest_face = Face_handle();
		for (Finite_faces_iterator fi = Base::finite_faces_begin();
				fi != Base::finite_faces_end(); ++fi) {
			Face_handle face = fi;
			if (face_has_point(face, point)) {
				nearest_face = face;
				break;
			}
		}

		if (nearest_face == Face_handle()) {
			std::cout << "free bird" << std::endl;
			return false;
		}

		Vertex_handle vertex = find_nearest_vertex(point, nearest_face);
		if (vertex != Vertex_handle()) {
			assign_sample_to_vertex(sample, vertex);
			return true;
		}

		Edge edge = find_nearest_edge(point, nearest_face);
		assign_sample_to_edge(sample, edge);
		return true;
	}

	bool face_has_point(Face_handle face, const Point& query) const {
		for (int i = 0; i < 3; ++i) {
			Edge edge(face, i);
			const Point& ps = source_vertex(edge)->point();
			const Point& pt = target_vertex(edge)->point();
			if (!compute_triangle_ccw(ps, pt, query))
				return false;
		}
		return true;
	}

	Vertex_handle find_nearest_vertex(const Point& point,
			Face_handle face) const {
		for (int i = 0; i < 3; ++i) {
			Vertex_handle vi = face->vertex(i);
			const Point& pi = vi->point();
			if (pi == point)
				return vi;
		}
		return Vertex_handle();
	}

	Edge find_nearest_edge(const Point& point, Face_handle face) const {
		FT min_dist2 = (std::numeric_limits<FT>::max)();
		Edge nearest(Face_handle(), 0);
		for (int i = 0; i < 3; ++i) {
			Edge edge(face, i);
			Segment segment = get_segment(edge);
			FT dist2 = compute_distance2(point, segment);
			if (dist2 < min_dist2) {
				min_dist2 = dist2;
				nearest = edge;
			}
		}

		if (nearest.first == Face_handle()) {
			std::cout << "nearest edge not found" << std::endl;
		}
		return nearest;
	}

	void assign_sample_to_vertex(Sample* sample, Vertex_handle vertex) {
		// DEBUG
		if (vertex->get_sample()) {
			std::cout << "assign to vertex: vertex already has sample"
					<< std::endl;
		}
		//

		sample->distance2() = 0.0;
		sample->coordinate() = 0.0;
		vertex->set_sample(sample);
	}

	void assign_sample_to_edge(Sample* sample, const Edge& edge) {
		Segment segment = get_segment(edge);
		const Point& query = sample->point();
		sample->distance2() = compute_distance2(query, segment);
		sample->coordinate() = compute_coordinate(query, segment);
		edge.first->add_sample(edge.second, sample);
	}

	FT compute_distance2(const Point& query, const Segment& segment) const {
		Line line = segment.supporting_line();
		if (line.has_on(query))
			return 0.0;

		Point proj = line.projection(query);
		return CGAL::squared_distance(query, proj);
	}

	FT compute_coordinate(const Point& q, const Segment& segment) const {
		const Point& p0 = segment.source();
		const Point& p1 = segment.target();
		Vector p0p1 = p1 - p0;
		Vector p0q = q - p0;
		FT t = (p0q * p0p1) / (p0p1 * p0p1);
		return t; // [0,1]
	}

	// SIGNED DISTANCE //

	// signed distance from line(a,b) to point t
	FT signed_distance(Vertex_handle a, Vertex_handle b,
			Vertex_handle t) const {
		const Point& pa = a->point();
		const Point& pb = b->point();
		const Point& pt = t->point();
		return compute_signed_distance(pa, pb, pt);
	}

	// signed distance from line(a,b) to point t
	FT compute_signed_distance(const Point& pa, const Point& pb,
			const Point& pt) const {
		if (pt == pa)
			return 0.0;
		if (pt == pb)
			return 0.0;
		if (pa == pb)
			return std::sqrt(CGAL::squared_distance(pa, pt));

		Vector vab = pb - pa;
		vab = vab / sqrt(vab * vab);
		Vector vab90(-vab.y(), vab.x());
		Vector vat = pt - pa;
		return (vat * vab90);
	}

	// signed distance from t to the intersection of line(a,b) and line(t,s)
	FT signed_distance_from_intersection(Vertex_handle a, Vertex_handle b,
			Vertex_handle t, Vertex_handle s) const {
		const Point& pa = a->point();
		const Point& pb = b->point();
		const Point& pt = t->point();
		const Point& ps = s->point();
		return compute_signed_distance_from_intersection(pa, pb, pt, ps);
	}

	// signed distance from t to the intersection of line(a,b) and line(t,s)
	FT compute_signed_distance_from_intersection(const Point& pa,
			const Point& pb, const Point& pt, const Point& ps) const {
		FT Dabt = compute_signed_distance(pa, pb, pt);
		if (Dabt == 0.0)
			return 0.0;

		Line lab(pa, pb - pa);
		Line lts(pt, ps - pt);

		FT Dqt = (std::numeric_limits<FT>::max)();
		CGAL::Object result = CGAL::intersection(lab, lts);
		const Point* iq = CGAL::object_cast<Point>(&result);
		if (iq)
			Dqt = std::sqrt(CGAL::squared_distance(*iq, pt));

		if (Dabt < 0.0)
			Dqt = -Dqt;
		return Dqt;
	}

	bool is_triangle_ccw(Vertex_handle a, Vertex_handle b,
			Vertex_handle c) const {
		const Point& pa = a->point();
		const Point& pb = b->point();
		const Point& pc = c->point();
		return compute_triangle_ccw(pa, pb, pc);
	}

	bool compute_triangle_ccw(const Point& pa, const Point& pb,
			const Point& pc) const {
		FT dist = compute_signed_distance(pa, pb, pc);
		return (dist > -EPS);
	}

	// COMBINATORIAL TESTS //

	// (a,b) is cyclic if (a,b,c) and (a,c,b) exist
	bool is_edge_cyclic(const Edge& edge) const {
		Vertex_handle f = opposite_vertex(edge);
		Vertex_handle b = opposite_vertex(twin_edge(edge));
		return (f == b);
	}

	// b from (a,b) is cyclic if (a,b,c) and (b,a,c) exist
	bool is_target_cyclic(const Edge& edge) const {
		if (!is_edge_cyclic(edge))
			return false;

		Edge twin = twin_edge(edge);
		Edge prev = prev_edge(twin);
		Face_handle fp = prev.first->neighbor(prev.second);
		Face_handle ft = twin.first->neighbor(twin.second);
		return (fp == ft);
	}

	bool is_flippable(const Edge& edge) const {
		Edge twin = twin_edge(edge);
		if (Base::is_infinite(twin.first))
			return false;
		if (Base::is_infinite(edge.first))
			return false;

		Vertex_handle vs = source_vertex(edge);
		Vertex_handle vt = target_vertex(edge);
		Vertex_handle vf = opposite_vertex(edge);
		Vertex_handle vb = opposite_vertex(twin);

		if (!is_triangle_ccw(vs, vb, vf))
			return false;
		if (!is_triangle_ccw(vt, vf, vb))
			return false;
		return true;
	}

	bool is_collapsible(const Edge& edge) const {
		if (!check_link_test(edge))
			return false;
		if (!check_kernel_test(edge))
			return false;
		return true;
	}

	bool check_link_test(const Edge& edge) const {
		Vertex_handle s = source_vertex(edge);
		Vertex_handle t = target_vertex(edge);

		if (s == t)
			return false;
		typename Vertex_handle_set::const_iterator it;

		Vertex_handle_set svertices;
		get_vertices_from_vertex_link(s, svertices);

		Vertex_handle_set tvertices;
		get_vertices_from_vertex_link(t, tvertices);

		// link(s) inter link(t)
		Vertex_handle_set ivertices;
		for (it = svertices.begin(); it != svertices.end(); ++it) {
			Vertex_handle v = *it;
			if (tvertices.find(v) != tvertices.end())
				ivertices.insert(v);
		}

		Vertex_handle_set evertices;
		get_vertices_from_edge_link(edge, evertices);

		// link(edge) =? link(s) inter link(t)
		if (evertices.size() != ivertices.size())
			return false;

		for (it = evertices.begin(); it != evertices.end(); ++it) {
			Vertex_handle v = *it;
			if (ivertices.find(v) == ivertices.end())
				return false;
		}
		return true;
	}

	bool check_kernel_test(const Edge& edge) const {
		Vertex_handle s = source_vertex(edge);
		Vertex_handle t = target_vertex(edge);

		Edge_list hull;
		get_edges_from_star_minus_link(s, hull);
		return is_in_kernel(t->point(), hull.begin(), hull.end());
	}

	template<class Iterator> // value_type = Edge
	bool is_in_kernel(const Point& query, Iterator begin, Iterator end) const {
		for (Iterator it = begin; it != end; ++it) {
			Edge edge = *it;
			const Point& pa = source_vertex(edge)->point();
			const Point& pb = target_vertex(edge)->point();
			if (!compute_triangle_ccw(pa, pb, query))
				return false;
		}
		return true;
	}

	// COLLAPSE //

	// (s,a,b) + (s,b,c) -> (s,a,c) + (a,b,c)
	// st = (source,target) from 'make_collapsible'
	// return (a,c)
	Edge flip(const Edge& sb, Edge& st, int verbose = 0) {
		Vertex_handle t = target_vertex(st);

		Edge sc = twin_edge(prev_edge(sb));
		Base::tds().flip(sb.first, sb.second);
		Edge ac = prev_edge(twin_edge(sc));

		Vertex_handle a = source_vertex(ac);
		if (a == t)
			st = prev_edge(ac);

		return ac;
	}

	void collapse(const Edge& edge, int verbose = 0) {
		if (is_edge_cyclic(edge)) {
			collapse_cyclic_edge(edge);
			return;
		}

		Edge twin = twin_edge(edge);
		Base::tds().join_vertices(twin);
	}

	// (a,b,c) + (c,b,a) + (a,c,i) + (c,a,j) ->
	// (a,c,i) + (c,a,j)
	void collapse_cyclic_edge(const Edge& bc, int verbose = 1) {
		if (verbose > 0)
			std::cout << "collapse_cyclic_edge ... ";

		Edge cb = twin_edge(bc);
		Face_handle abc = bc.first;
		Face_handle cba = cb.first;

		Vertex_handle b = source_vertex(bc);
		Vertex_handle c = target_vertex(bc);
		Vertex_handle a = opposite_vertex(bc);

		Edge ac = twin_edge(next_edge(bc));
		Edge ca = twin_edge(prev_edge(cb));

		a->set_face(ac.first);
		c->set_face(ca.first);
		ac.first->set_neighbor(ac.second, ca.first);
		ca.first->set_neighbor(ca.second, ac.first);

		this->delete_face(abc);
		this->delete_face(cba);
		this->delete_vertex(b);

		if (verbose > 0)
			std::cout << "done" << std::endl;
	}

	//TODO IV remove --------
	void print_edge(Reconstruction_edge_2 edge) {
		int i = ((edge).edge()).second;
		Point a = ((edge).edge()).first->vertex((i+1)%3)->point();
		Point b = ((edge).edge()).first->vertex((i+2)%3)->point();
		std::cout <<"( " << (edge).priority()  <<  ") ( " << a << " , " << b << " )" << std::endl;
	}
	//--------

		template <class Iterator> // value_type = Edge
		bool make_collapsible(Edge& edge, Iterator begin, Iterator end, int verbose = 0)
		{
		        Vertex_handle source = source_vertex(edge);
		        Vertex_handle target = target_vertex(edge);

		        MultiIndex multi_ind;
		        for (Iterator it = begin; it != end; ++it)
		        {
		            Edge ab = twin_edge(*it);
		            Vertex_handle a = source_vertex(ab);
		            Vertex_handle b = target_vertex(ab);
		            FT D = signed_distance_from_intersection(a, b, target, source);
		            if (D < 0.0) {
		            	multi_ind.insert(Reconstruction_edge_2(ab, D));
		            }
		        }


		        int nb_flips = 0;
		        while (!multi_ind.empty())
		        {
		        	Reconstruction_edge_2 pedge = *(multi_ind.template get<1>()).begin();
					FT Dbc = pedge.priority();
		            Edge bc = pedge.edge();
		            (multi_ind.template get<0>()).erase(pedge);

		            Edge sb = prev_edge(bc);
		            Edge ab = prev_edge(twin_edge(sb));
		            Edge sc = twin_edge(next_edge(bc));
		            Edge cd = next_edge(sc);

		            Vertex_handle a = source_vertex(ab);
		            Vertex_handle b = source_vertex(bc);
		            Vertex_handle c = target_vertex(bc);
		            Vertex_handle d = target_vertex(cd);

		            FT Dac =  std::numeric_limits<FT>::lowest();
		            if (a != c && is_triangle_ccw(a, b, c))
		                Dac = signed_distance_from_intersection(a, c, target, source);

		            FT Dbd =  std::numeric_limits<FT>::lowest();
		            if (b != d && is_triangle_ccw(b, c, d))
		                Dbd = signed_distance_from_intersection(b, d, target, source);

		            if (Dac ==  std::numeric_limits<FT>::lowest() && Dbd ==
		            		std::numeric_limits<FT>::lowest())
		            {
		                // TODO: IV comment in std::cerr << red << "---
		            	//No flips available ---" << white << std::endl;
		            	std::cerr << "--- No flips available ---"  << std::endl;
		                return false;
		            }

		            if ((std::max)(Dac, Dbd) + EPS < Dbc)
		            {
		                std::cerr.precision(10);
		                // TODO: IV comment in std::cerr << red << "--
		                //- Flip makes kernel worse ---" << white << std::endl;
		                std::cerr << "--- Flip makes kernel worse ---"
		                		<< std::endl;
		                std::cerr << Dac << " or " << Dbd << " vs "
		                		<< Dbc << std::endl;
		                std::cerr << "a: " << a->point() << std::endl;
		                std::cerr << "b: " << b->point() << std::endl;
		                std::cerr << "c: " << c->point() << std::endl;
		                std::cerr << "d: " << d->point() << std::endl;
		                std::cerr << "t: " << target->point() << std::endl;
		                std::cerr << "diff = " << Dbc - (std::max)(Dac, Dbd) << std::endl;
		                return false;
		            }

		            if (Dac > Dbd)
		            {
		                (multi_ind.template get<0>()).erase(Reconstruction_edge_2(ab));

		                Edge ac = flip(sb, edge, verbose);
		                if (Dac < 0.0) {
		                	multi_ind.insert(Reconstruction_edge_2(ac, Dac));
		                }
		            }
		            else
		            {
		                (multi_ind.template get<0>()).erase(Reconstruction_edge_2(cd));
		                Edge bd = flip(sc, edge, verbose);
		                if (Dbd < 0.0) {
		                	multi_ind.insert(Reconstruction_edge_2(bd, Dbd));
		                }
		            }
		            nb_flips++;
		        }

		        if (verbose > 1)
		        	//TODO: IV Comment in std::cerr << red << "---
		        	//Flip makes kernel worse ---" << white << std::endl;
		            std::cerr  << "Nb flips: "  << nb_flips << std::endl;

		        return true;
		    }

	int random_int(const int min, const int max) {
		int range = max - min;
		return min + int((double(rand()) / double(RAND_MAX)) * range);
	}

};
} //namespace CGAL

#endif
