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

#ifndef RECONSTRUCTION_SIMPLIFICATION_2_H_
#define RECONSTRUCTION_SIMPLIFICATION_2_H_

#include <CGAL/Reconstruction_triangulation_2.h>
#include <CGAL/Cost.h>
#include <CGAL/Reconstruction_edge_2.h>
#include <CGAL/Sample.h>
#include <CGAL/console_color.h>

#include <CGAL/property_map.h>

#include <iterator>
#include <iostream>
#include <list>
#include <algorithm>

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/mem_fun.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/identity.hpp>

namespace CGAL {

/*!
\ingroup PkgReconstructionSimplification2Classes



\brief The class `Reconstruction_simplification_2` is the main class
for executing the reconstruction and simplification tasks.
It takes an InputIterator which can be used to traverse a collection
of point-mass pairs, where the points and their masses are accessed
via the PointPMap and MassPMap `PropertyMap`s respectively.


\tparam Kernel is the geometric kernel, used throughout the reconstruction and
					simplification task.

\tparam InputIterator is the iterator type of the algorithm input.

\tparam PointPMap is a model of `ReadablePropertyMap` with a value_type = `Point_2`

\tparam MassPMap   is a model of `ReadablePropertyMap` with a value_type = `FT`

 */
template<class Kernel, class InputIterator, class PointPMap, class MassPMap>
class Reconstruction_simplification_2 {
public:


	/*!
		Number type.
	*/
	typedef typename Kernel::FT FT;

	/*!
		Point type.
	*/
	typedef typename Kernel::Point_2 Point;
	/*!
		Vector type.
	*/
	typedef typename Kernel::Vector_2 Vector;

	/*!
	The Output simplex.
	*/
	typedef Reconstruction_triangulation_2<Kernel> Triangulation;

	 /// \cond SKIP_IN_MANUAL
	typedef typename Triangulation::Vertex Vertex;
	typedef typename Triangulation::Vertex_handle Vertex_handle;
	typedef typename Triangulation::Vertex_iterator Vertex_iterator;
	typedef typename Triangulation::Vertex_circulator Vertex_circulator;
	typedef typename Triangulation::Finite_vertices_iterator
			Finite_vertices_iterator;

	typedef typename Triangulation::Edge Edge;
	typedef typename Triangulation::Edge_iterator Edge_iterator;
	typedef typename Triangulation::Edge_circulator Edge_circulator;
	typedef typename Triangulation::Finite_edges_iterator Finite_edges_iterator;

	typedef typename Triangulation::Face Face;
	typedef typename Triangulation::Face_handle Face_handle;
	typedef typename Triangulation::Face_iterator Face_iterator;
	typedef typename Triangulation::Face_circulator Face_circulator;
	typedef typename Triangulation::Finite_faces_iterator Finite_faces_iterator;

	typedef typename Triangulation::Vertex_handle_map Vertex_handle_map;
	typedef typename Triangulation::Face_handle_map Face_handle_map;

	typedef typename Triangulation::Vertex_handle_set Vertex_handle_set;
	typedef typename Triangulation::Edge_set Edge_set;

	typedef typename Triangulation::Edge_list Edge_list;

	typedef typename Triangulation::Cost Cost;
	typedef typename Triangulation::Sample Sample;
	typedef typename Triangulation::Sample_list Sample_list;
	typedef typename Triangulation::Sample_list_const_iterator
			Sample_list_const_iterator;

	typedef typename Triangulation::Point_list Point_list;
	typedef typename Triangulation::Point_list_const_iterator
			Point_list_const_iterator;

	typedef typename Triangulation::PSample PSample;
	typedef typename Triangulation::SQueue SQueue;

	typedef typename Triangulation::Reconstruction_edge_2 Reconstruction_edge_2;

	typedef typename Triangulation::MultiIndex MultiIndex;

protected:
	Triangulation m_dt;
	MultiIndex m_mindex;
	int m_ignore;
	int m_verbose;
	int m_mchoice;  // # Edges
	bool m_use_flip;
	double m_alpha; // [0, 1]
	double m_norm_tol; // [0,BBOX]
	double m_tang_tol; // [0,BBOX]
	double m_ghost; // ghost vs solid
	unsigned m_relocation; // # relocations

    // bbox
    double m_bbox_x;
    double m_bbox_y;
    double m_bbox_size;

	InputIterator start;
	InputIterator beyond;
	PointPMap point_pmap;
	MassPMap  mass_pmap;

	  /// \endcond

	// Public methods
	public:

	  /// \name Creation
	  /// @{

	/*!
	     \Instantiates a new Reconstruction_simplification_2.

	     \details Instantiates a new Reconstruction_simplification_2 object
	     	 	  for a given collection of point-mass pairs.

	     \param start_itr An InputIterator pointing the the first point-mass
	     	 	 	 	 	 pair in a collection.
   	   	 \param beyond_itr An InputIterator pointing beyond the last point-mass
	     	 	 	 	 	 pair in a collection.
	     \param in_point_pmap A `ReadablePropertyMap` used to access the input points

	     \param in_mass_pmap A `ReadablePropertyMap` used to access the input points' mass.
	*/
	Reconstruction_simplification_2(InputIterator start_itr,
									InputIterator beyond_itr,
									PointPMap in_point_pmap,
									MassPMap  in_mass_pmap) {


		start  = start_itr;
		beyond = beyond_itr;

		point_pmap = in_point_pmap;
		mass_pmap  = in_mass_pmap;

		initialize_parameters();
	}

	  /// @}


	 /// \cond SKIP_IN_MANUAL

	Reconstruction_simplification_2() {
		initialize_parameters();
	}


	~Reconstruction_simplification_2() {
		clear();
	}

	void initialize_parameters() {


		m_verbose = 0;
		m_mchoice = 0;
		m_use_flip = true;
		m_alpha = 0.5;
		m_norm_tol = 1.0;
		m_tang_tol = 1.0;
		m_ghost = 1.0;
		m_relocation = 0;

		m_bbox_x = 0.0;
        m_bbox_y = 0.0;
        m_bbox_size = 1.0;

        m_ignore = 0;
	}

	//Function if one wants to create a Reconstruction_simplification_2
	//without specifying the input yet in the constructor
	void initialize(InputIterator start_itr,
									InputIterator beyond_itr,
									PointPMap in_point_pmap,
									MassPMap  in_mass_pmap) {

		start  = start_itr;
		beyond = beyond_itr;

		point_pmap = in_point_pmap;
		mass_pmap  = in_mass_pmap;

		initialize();

	}
	 /// \endcond


	/*!
		First function to be called after instantiating a new
		Reconstruction_simplification_2 object.
		It computes an bounding box around the input points and creates a first
		(fine) output simplex as well as an initial transportation plan. This
		first output simplex
	  */
	void initialize() {

		clear();

		insert_loose_bbox(m_bbox_x, m_bbox_y, 2 * m_bbox_size);

		init(start, beyond);

		std::list<Sample*> m_samples;
		for (InputIterator it = start; it != beyond; it++) {
			Point point = get(point_pmap, *it);
			FT     mass = get( mass_pmap, *it);
			Sample* s = new Sample(point, mass);
			m_samples.push_back(s);
		}
		assign_samples(m_samples.begin(), m_samples.end());
	}

	/*!

		 Returns the solid edges present after the reconstruction process.

	 	 \details Instantiates a new Reconstruction_simplification_2 object
				  for a given collection of point-mass pairs.

		 \tparam OutputModule Concept for accessing the output

		 \param output An OutputModule in which the solid edges and vertics are
		  	  stored.
	*/
	template <class OutputModule>
	void extract_solid_elements(OutputModule& output) {
		output.store_marked_elements(m_dt, m_ignore);
	}

	 /// \cond SKIP_IN_MANUAL
	template <class Vector>
	Vector random_vec(const double scale)
	{
	    double dx = -scale + (double(rand()) / double(RAND_MAX)) * 2* scale;
	    double dy = -scale + (double(rand()) / double(RAND_MAX)) * 2* scale;
	    return Vector(dx, dy);
	}
	 /// \endcond

	/*!
			TODO COMMENT and change to perturb
	 */
	 void noise(const FT scale)
    {
        std::cerr << "noise by " << scale << "...";
         for (InputIterator it = start; it != beyond; it++)
        {
            Point point = get(point_pmap, *it);
            point = point + random_vec<Vector>(scale);
        }
        std::cerr << "done" << std::endl;
    }

	/*!
		TODO COMMENT
	 */
	void normalize_points()
    {
        noise(1e-5);
        compute_bbox(m_bbox_x, m_bbox_y, m_bbox_size);
        if (m_bbox_size == 0.0) return;

        Point center(m_bbox_x, m_bbox_y);
        for (InputIterator it = start; it != beyond; ++it)
        {
        	Point point = get(point_pmap, *it);
			Vector vec = (point - center) / m_bbox_size;
			point = CGAL::ORIGIN + vec;
        }
        m_bbox_x = m_bbox_y = 0.0;
        m_bbox_size = 1.0;
    }

	 /// \cond SKIP_IN_MANUAL
	 void compute_bbox(double &x, double &y, double &scale)
     {

        FT x_min, x_max, y_min, y_max;
        InputIterator it = start;
        Point p = get(point_pmap, *it);
        x_min = x_max = p.x();
        y_min = y_max = p.y();
        ++it;
        for ( ; it != beyond; ++it)
        {
        	p = get(point_pmap, *it);
            x_min = (std::min)(x_min, p.x());
            x_max = (std::max)(x_max, p.x());
            y_min = (std::min)(y_min, p.y());
            y_max = (std::max)(y_max, p.y());
        }

        x = 0.5 * (x_min + x_max);
        y = 0.5 * (y_min + y_max);
        scale = (std::max)(x_max - x_min, y_max - y_min);
        if (scale == 0.0) scale = 1.0;
    }

	void clear() {
		m_dt.clear();
		m_mindex.clear();
	}

	double time_duration(const double init) {
		return (clock() - init) / CLOCKS_PER_SEC;
	}

	void set_mchoice(const int mchoice) {
		m_mchoice = mchoice;
	}

	void set_verbose(const int verbose) {
		m_verbose = verbose;
	}


	void set_alpha(const double alpha) {
		m_alpha = alpha;
	}


	void set_use_flip(const bool use_flip) {
		m_use_flip = use_flip;
	}


	void set_norm_tol(const double norm_tol) {
		m_norm_tol = norm_tol;
	}


	double get_norm_tol() const {
		return m_norm_tol;
	}


	void set_tang_tol(const double tang_tol) {
		m_tang_tol = tang_tol;
	}


	double get_tang_tol() const {
		return m_tang_tol;
	}


	void set_relocation(const unsigned relocation) {
		m_relocation = relocation;
	}


	unsigned get_relocation() const {
		return m_relocation;
	}


	void set_ghost(const double g) {
		m_ghost = g;
		m_dt.ghost_factor() = m_ghost;
	}


	double get_ghost() {
		return m_ghost;
	}

	// INIT //


	void insert_loose_bbox(const double x, const double y, const double size) {
		double timer = clock();
		std::cerr << yellow << "insert loose bbox" << white << "...";

		int nb = m_dt.number_of_vertices();
		insert_point(Point(x - size, y - size), true, nb++);
		insert_point(Point(x - size, y + size), true, nb++);
		insert_point(Point(x + size, y + size), true, nb++);
		insert_point(Point(x + size, y - size), true, nb++);

		std::cerr << yellow << "done" << white << " (" << nb << " vertices, "
				<< yellow << time_duration(timer) << white << " s)"
				<< std::endl;
	}

	template<class Iterator>  // value_type = Point*
	void init(Iterator begin, Iterator end) {
		double timer = clock();
		std::cerr << yellow << "init" << white << "...";

		int nb = m_dt.number_of_vertices();
		m_dt.infinite_vertex()->pinned() = true;
		for (Iterator it = begin; it != end; it++) {
			Point point = get(point_pmap, *it);
			Vertex_handle vertex = insert_point(point, false, nb++);
		}

		std::cerr << yellow << "done" << white << " (" << nb << " vertices, "
				<< yellow << time_duration(timer) << white << " s)"
				<< std::endl;
	}

	Vertex_handle insert_point(const Point& point, const bool pinned,
			const int id) {
		Vertex_handle v = m_dt.insert(point);
		v->pinned() = pinned;
		v->id() = id;
		return v;
	}

	// ASSIGNMENT //

	void cleanup_assignments() {
		m_dt.cleanup_assignments();
	}

	template<class Iterator>  // value_type = Sample*
	void assign_samples(Iterator begin, Iterator end) {
		double timer = clock();
		std::cerr << yellow << "assign samples" << white << "...";

		m_dt.assign_samples(begin, end);
		m_dt.reset_all_costs();

		std::cerr << yellow << "done" << white << " (" << yellow
				<< time_duration(timer) << white << " s)" << std::endl;
	}

	void reassign_samples() {
		Sample_list samples;
		m_dt.collect_all_samples(samples);
		m_dt.cleanup_assignments();
		m_dt.assign_samples(samples.begin(), samples.end());
		m_dt.reset_all_costs();
	}

	void reassign_samples_around_vertex(Vertex_handle vertex) {
		Sample_list samples;
		m_dt.collect_samples_from_vertex(vertex, samples, true);
		m_dt.assign_samples(samples.begin(), samples.end());

		Edge_list hull;
		m_dt.get_edges_from_star_minus_link(vertex, hull, true);
		update_cost(hull.begin(), hull.end());
	}


	bool decimate() {
		bool ok;
		Reconstruction_edge_2 pedge;
		ok = pick_edge(m_mchoice, pedge);
		if (!ok)
			return false;

		ok = do_collapse(pedge.edge());
		if (!ok)
			return false;
		return true;
	}

	bool do_collapse(Edge edge) {
		bool ok;
		Vertex_handle s = m_dt.source_vertex(edge);
		Vertex_handle t = m_dt.target_vertex(edge);

		if (m_verbose > 0) {
			std::cerr << std::endl << green << "do collapse " << white << "("
					<< s->id() << "->" << t->id() << ") ... " << std::endl;
		}

		Sample_list samples;
		m_dt.collect_samples_from_vertex(s, samples, true);

		Edge_list hull;
		m_dt.get_edges_from_star_minus_link(s, hull, true);

		if (m_mchoice == 0)
			remove_stencil_from_pqueue(hull.begin(), hull.end());

		if (m_use_flip)
			ok = m_dt.make_collapsible(edge, hull.begin(), hull.end(),
					m_verbose);

		// debug test
		ok = m_dt.check_kernel_test(edge);
		if (!ok) {
			std::cerr << red << "do_collapse: kernel test failed: " << white
					<< std::endl;
			return false;
		}
		//

		m_dt.collapse(edge, m_verbose);

		m_dt.assign_samples(samples.begin(), samples.end());

		update_cost(hull.begin(), hull.end());

		if (m_mchoice == 0)
			push_stencil_to_pqueue(hull.begin(), hull.end());

		for (unsigned i = 0; i < m_relocation; ++i) {
			relocate_one_ring(hull.begin(), hull.end());
		}

		if (m_verbose > 0) {
			std::cerr << green << "done" << std::endl;
		}

		return true;
	}

	bool simulate_collapse(const Edge& edge, Cost& cost) {
		bool ok;
		Vertex_handle s = m_dt.source_vertex(edge);
		Vertex_handle t = m_dt.target_vertex(edge);

		if (m_verbose > 1) {
			std::cerr << green << "simulate collapse " << white << "("
					<< s->id() << "->" << t->id() << ") ... " << std::endl;
		}

		Triangulation copy;
		Edge copy_edge = copy_star(edge, copy);
		Vertex_handle copy_source = copy.source_vertex(copy_edge);

		if (m_use_flip) {
			Edge_list copy_hull;
			copy.get_edges_from_star_minus_link(copy_source, copy_hull, true);
			ok = copy.make_collapsible(copy_edge, copy_hull.begin(),
					copy_hull.end(), m_verbose);
			if (!ok) {
				std::cerr << yellow << "simulation: failed (make collapsible)"
						<< white << std::endl;
				return false;
			}
		}

		ok = copy.check_kernel_test(copy_edge);
		if (!ok) {
			std::cerr << yellow << "simulation: failed (kernel test)" << white
					<< std::endl;
			return false;
		}

		copy.collapse(copy_edge, m_verbose);

		Sample_list samples;
		m_dt.collect_samples_from_vertex(s, samples, false);

		backup_samples(samples.begin(), samples.end());
		copy.assign_samples_brute_force(samples.begin(), samples.end());
		copy.reset_all_costs();
		cost = copy.compute_total_cost();
		restore_samples(samples.begin(), samples.end());

		if (m_verbose > 1) {
			std::cerr << green << "done" << white << std::endl;
		}

		return true;
	}

	template<class Iterator> // value_type = Sample*
	void backup_samples(Iterator begin, Iterator end) {
		for (Iterator it = begin; it != end; ++it) {
			Sample* sample = *it;
			sample->backup();
		}
	}

	template<class Iterator> // value_type = Sample*
	void restore_samples(Iterator begin, Iterator end) {
		for (Iterator it = begin; it != end; ++it) {
			Sample* sample = *it;
			sample->restore();
		}
	}

	// PEDGE //

	bool create_pedge(const Edge& edge, Reconstruction_edge_2& pedge) {
		Cost after_cost;
		bool ok = simulate_collapse(edge, after_cost);
		if (!ok)
			return false;

		bool within_tol = is_within_tol(after_cost);
		if (!within_tol)
			return false;

		Vertex_handle source = m_dt.source_vertex(edge);
		Cost before_cost = m_dt.compute_cost_around_vertex(source);

		FT before = before_cost.finalize(m_alpha);
		FT after = after_cost.finalize(m_alpha);
		pedge = Reconstruction_edge_2(edge, before, after);
		return true;
	}

	bool is_within_tol(const Cost& cost) const {
		if (cost.max_norm() > m_norm_tol)
			return false;
		if (cost.max_tang() > m_tang_tol)
			return false;
		return true;
	}

	// COST //

	void init_cost() {
		m_dt.reset_all_costs();
	}

	template<class Iterator> // value_type = Edge
	void update_cost(Iterator begin, Iterator end) {
		Edge_list edges;
		collect_cost_stencil(m_dt, begin, end, edges);

		typename Edge_list::iterator ei;
		for (ei = edges.begin(); ei != edges.end(); ++ei) {
			Edge edge = *ei;
			m_dt.update_cost(edge);
		}
	}

	template<class Iterator> // value_type = Edge
	void collect_cost_stencil(const Triangulation& mesh, Iterator begin,
			Iterator end, Edge_list& edges) {
		Edge_set done;
		Edge_list fifo;
		for (Iterator it = begin; it != end; ++it) {
			Edge edge = *it;
			fifo.push_back(edge);
			done.insert(edge);
		}

		while (!fifo.empty()) {
			Edge edge = fifo.front();
			fifo.pop_front();

			edge = mesh.twin_edge(edge);
			edges.push_back(edge);

			Edge next = mesh.next_edge(edge);
			if (done.insert(next).second)
				fifo.push_back(next);

			Edge prev = mesh.prev_edge(edge);
			if (done.insert(prev).second)
				fifo.push_back(prev);
		}
	}

	// PQUEUE (MCHOICE or EXHAUSTIVE) //

	bool pick_edge(int nb, Reconstruction_edge_2& best_pedge) {
		if (m_dt.number_of_faces() < 2)
			return false;

		int ne = int(2 * m_dt.tds().number_of_edges());
		if (nb > ne)
			nb = ne;

		bool ok;
		if (nb == 0) {
			ok = pick_edge_from_pqueue(best_pedge);
			return ok;
		}
		m_mindex.clear();

		if (nb == ne) {
			ok = pick_edge_brute_force(best_pedge);
			return ok;
		}

		ok = pick_edge_randomly(nb, best_pedge);
		return ok;
	}

	bool pick_edge_from_pqueue(Reconstruction_edge_2& best_pedge) {
		if (m_mindex.empty())
			populate_pqueue();
		if (m_mindex.empty())
			return false;
		best_pedge = *(m_mindex.template get<1>()).begin();
		(m_mindex.template get<0>()).erase(best_pedge);
		return true;
	}

	bool pick_edge_brute_force(Reconstruction_edge_2& best_pedge) {
		MultiIndex mindex;
		Finite_edges_iterator ei;
		for (ei = m_dt.finite_edges_begin(); ei != m_dt.finite_edges_end();
				++ei) {
			Edge edge = *ei;
			push_to_mindex(edge, mindex);


			edge = m_dt.twin_edge(edge);
			push_to_mindex(edge, mindex);
		}
		if (mindex.empty())
			return false;
		best_pedge = *(mindex.template get<1>()).begin();
		return true;
	}

	bool pick_edge_randomly(int nb, Reconstruction_edge_2& best_pedge) {
		MultiIndex mindex;
		for (int i = 0; i < nb; ++i) {
			Reconstruction_edge_2 pedge;
			if (random_pedge(pedge))
				mindex.insert(pedge);
		}
		if (mindex.empty())
			return false;
		best_pedge = *(mindex.template get<1>()).begin();
		return true;
	}

	void populate_pqueue() {
		Finite_edges_iterator ei;
		for (ei = m_dt.finite_edges_begin(); ei != m_dt.finite_edges_end();
				++ei) {
			Edge edge = *ei;
			push_to_mindex(edge, m_mindex);

			edge = m_dt.twin_edge(edge);
			push_to_mindex(edge, m_mindex);
		}
	}


	bool push_to_mindex(const Edge& edge, MultiIndex& mindex) {
		if (m_dt.is_pinned(edge))
			return false;
		if (m_dt.is_target_cyclic(edge))
			return false;

		Reconstruction_edge_2 pedge;
		bool ok = create_pedge(edge, pedge);
		if (!ok)
			return false;
		mindex.insert(pedge);
		return true;
	}



	bool random_pedge(Reconstruction_edge_2& pedge) {
		for (unsigned i = 0; i < 10; ++i) {
			Edge edge = m_dt.random_finite_edge();
			if (m_dt.is_pinned(edge))
				continue;
			if (m_dt.is_target_cyclic(edge))
				continue;
			bool ok = create_pedge(edge, pedge);
			if (ok)
				return true;
		}
		return false;
	}

	template<class Iterator> // value_type = Edge
	void remove_stencil_from_pqueue(Iterator begin, Iterator end) {
		if (m_mindex.empty())
			return;

		Edge_list edges;
		collect_pqueue_stencil(m_dt, begin, end, edges);

		typename Edge_list::const_iterator ei;
		for (ei = edges.begin(); ei != edges.end(); ++ei) {
			Edge edge = *ei;
			(m_mindex.template get<0>()).erase(Reconstruction_edge_2(edge));
		}
	}

	template<class Iterator> // value_type = Edge
	void push_stencil_to_pqueue(Iterator begin, Iterator end) {
		Edge_list edges;
		collect_pqueue_stencil(m_dt, begin, end, edges);

		typename Edge_list::const_iterator ei;
		for (ei = edges.begin(); ei != edges.end(); ++ei) {
			Edge edge = *ei;
			push_to_mindex(edge, m_mindex);
		}
	}

	template<class Iterator> // value_type = Edge
	void collect_pqueue_stencil(const Triangulation& mesh, Iterator begin,
			Iterator end, Edge_list& edges) {
		Vertex_handle_set vertex_set;
		for (Iterator it = begin; it != end; ++it) {
			Edge edge = *it;
			Edge twin = mesh.twin_edge(edge);

			Vertex_handle s = mesh.source_vertex(edge);
			if (!s->pinned())
				vertex_set.insert(s);

			Vertex_handle t = mesh.target_vertex(edge);
			if (!t->pinned())
				vertex_set.insert(t);

			Vertex_handle f = mesh.opposite_vertex(edge);
			if (!f->pinned())
				vertex_set.insert(f);

			Vertex_handle b = mesh.opposite_vertex(twin);
			if (!b->pinned())
				vertex_set.insert(b);
		}

		typename Vertex_handle_set::const_iterator vi;
		for (vi = vertex_set.begin(); vi != vertex_set.end(); ++vi) {
			Vertex_handle v = *vi;
			Edge_circulator ecirc = mesh.incident_edges(v);
			Edge_circulator eend = ecirc;
			CGAL_For_all(ecirc, eend)
			{
				Edge edge = *ecirc;
				if (mesh.source_vertex(edge) != v)
					edge = mesh.twin_edge(edge);
				edges.push_back(edge);
			}
		}
	}

	// COPY STAR //

	// edge must not be pinned or have cyclic target
	Edge copy_star(const Edge& edge, Triangulation& copy) {
		copy.tds().set_dimension(2);
		copy.infinite_vertex()->pinned() = true;

		// copy vertices
		Vertex_handle_map cvmap;

		Vertex_handle s = m_dt.source_vertex(edge);
		Vertex_handle cs = copy.tds().create_vertex();
		cvmap[s] = copy_vertex(s, cs);

		Vertex_circulator vcirc = m_dt.incident_vertices(s);
		Vertex_circulator vend = vcirc;
		CGAL_For_all(vcirc, vend)
		{
			Vertex_handle v = vcirc;
			if (cvmap.find(v) == cvmap.end()) {
				Vertex_handle cv = copy.tds().create_vertex();
				cvmap[v] = copy_vertex(v, cv);
			}
		}

		// copy faces
		Face_handle_map cfmap;
		Face_circulator fcirc = m_dt.incident_faces(s);
		Face_circulator fend = fcirc;
		CGAL_For_all(fcirc, fend)
		{
			Face_handle f = fcirc;
			Face_handle cf = copy.tds().create_face();
			cfmap[f] = copy_face(f, cf, cvmap);
		}

		// set neighbors
		fcirc = m_dt.incident_faces(s);
		fend = fcirc;
		CGAL_For_all(fcirc, fend)
		{
			Face_handle f = fcirc;
			copy_neighbors(f, s, cvmap, cfmap);
		}

		// make copy homeomorphic to S^2
		close_copy_mesh(cs, copy);

		// copy samples surrounding star
		copy_samples(s, cs, cfmap, copy);

		// get copy of edge
		Edge copy_edge = get_copy_edge(edge, cvmap, cfmap);
		return copy_edge;
	}

	Vertex_handle copy_vertex(Vertex_handle v0, Vertex_handle v1) {
		v1->id() = v0->id();
		v1->set_point(v0->point());
		v1->pinned() = v0->pinned();
		v1->set_sample(v0->get_sample());
		return v1;
	}

	Face_handle copy_face(Face_handle f0, Face_handle f1,
			Vertex_handle_map& vmap) {
		for (unsigned i = 0; i < 3; ++i) {
			Vertex_handle v0i = f0->vertex(i);
			Vertex_handle v1i = vmap[v0i];
			f1->set_vertex(i, v1i);
			v1i->set_face(f1);
		}
		return f1;
	}

	void copy_neighbors(Face_handle f, Vertex_handle v, Vertex_handle_map& vmap,
			Face_handle_map& fmap) {
		int i = f->index(v);
		Face_handle cf = fmap[f];
		Vertex_handle cv = vmap[v];

		if (fmap.find(f->neighbor(i)) != fmap.end()) {
			Face_handle fi = f->neighbor(i);
			Face_handle cfi = fmap[fi];
			cf->set_neighbor(i, cfi);
		}

		for (unsigned j = 0; j < 2; ++j) {
			i = (i + 1) % 3;
			Face_handle fi = f->neighbor(i);
			Face_handle cfi = fmap[fi];
			cf->set_neighbor(i, cfi);
		}
	}

	void close_copy_mesh(Vertex_handle vertex, Triangulation& copy) {
		std::vector<Face_handle> outer_faces;

		Face_circulator fcirc = copy.incident_faces(vertex);
		Face_circulator fend = fcirc;
		CGAL_For_all(fcirc, fend)
		{
			Face_handle face = fcirc;
			int i = face->index(vertex);

			if (face->neighbor(i) != Face_handle())
				continue;

			Vertex_handle v1 = face->vertex((i + 1) % 3);
			Vertex_handle v2 = face->vertex((i + 2) % 3);

			Face_handle outer = copy.tds().create_face();
			outer->set_vertex(0, copy.infinite_vertex());
			outer->set_vertex(1, v2);
			outer->set_vertex(2, v1);

			face->set_neighbor(i, outer);
			outer->set_neighbor(0, face);

			outer_faces.push_back(outer);
		}

		for (unsigned i = 0; i < outer_faces.size(); ++i) {
			unsigned j = (i + 1) % outer_faces.size();
			outer_faces[i]->set_neighbor(2, outer_faces[j]);
			outer_faces[j]->set_neighbor(1, outer_faces[i]);
		}

		if (!outer_faces.empty())
			copy.infinite_vertex()->set_face(outer_faces[0]);
	}

	void copy_samples(Vertex_handle vertex, Vertex_handle copy_vertex,
			Face_handle_map& fmap, Triangulation& copy) {
		Face_circulator fcirc = m_dt.incident_faces(vertex);
		Face_circulator fend = fcirc;
		CGAL_For_all(fcirc, fend)
		{
			Face_handle face = fcirc;
			int index = face->index(vertex);
			Edge twin = m_dt.twin_edge(Edge(face, index));

			Face_handle copy_face = fmap[face];
			index = copy_face->index(copy_vertex);
			Edge copy_twin = copy.twin_edge(Edge(copy_face, index));

			Sample_list samples;
			m_dt.collect_samples_from_edge(twin, samples);
			copy_twin.first->samples(copy_twin.second) = samples;
		}
		copy_vertex->set_sample(NULL);
	}

	Edge get_copy_edge(const Edge& edge, Vertex_handle_map& vmap,
			Face_handle_map& fmap) {
		Face_handle f = edge.first;
		Vertex_handle v = f->vertex(edge.second);

		Face_handle cf = fmap[f];
		Vertex_handle cv = vmap[v];

		return Edge(cf, cf->index(cv));
	}

	// RELOCATION //

	void relocate_one_vertex(Vertex_handle vertex) {
		std::swap(vertex->point(), vertex->relocated());
		reassign_samples_around_vertex(vertex);
	}

	template<class Iterator> // value_type = Edge
	void relocate_one_ring(Iterator begin, Iterator end) {
		Vertex_handle_set vertices;
		for (Iterator it = begin; it != end; ++it) {
			Edge edge = *it;
			vertices.insert(m_dt.source_vertex(edge));
			vertices.insert(m_dt.target_vertex(edge));
		}

		typename Vertex_handle_set::const_iterator vi;
		for (vi = vertices.begin(); vi != vertices.end(); ++vi) {
			Vertex_handle v = *vi;
			if (v->pinned())
				continue;
			v->relocated() = compute_relocation(v);
		}

		for (vi = vertices.begin(); vi != vertices.end(); ++vi) {
			Vertex_handle v = *vi;
			if (v->pinned())
				continue;
			if (v->point() == v->relocated())
				continue;

			Edge_list hull;
			m_dt.get_edges_from_star_minus_link(v, hull, false);
			bool ok = m_dt.is_in_kernel(v->relocated(), hull.begin(),
					hull.end());

			if (ok) {
				// do relocation
				FT norm_bef = m_dt.compute_cost_around_vertex(v).norm();
				relocate_one_vertex(v);
				FT norm_aft = m_dt.compute_cost_around_vertex(v).norm();

				if (norm_bef < norm_aft) {
					// undo relocation
					relocate_one_vertex(v);
				} else if (m_mchoice == 0) {
					// update queue
					hull.clear();
					m_dt.get_edges_from_star_minus_link(v, hull, true);
					remove_stencil_from_pqueue(hull.begin(), hull.end());
					push_stencil_to_pqueue(hull.begin(), hull.end());
				}
			}
		}
	}
	 /// \endcond


	/*!
	Since noise and missing data result may prevent the reconstructed shape to
	have sharp corners, the algorithm offers the possibility to automatically
	relocate vertices after each edge contraction. The new location of the
	vertices is chosen such that the fitting of the output triangulation to the
	input points is improved. This is achieved by minimizing the normal component
	of the weighted $L_2$ distance. The vertices then get relocated only if the
	 resulting triangulation is still embeddable.
	  */
	void relocate_all_vertices() {
		double timer = clock();
		std::cerr << yellow << "relocate all" << white << "...";

		m_mindex.clear(); // pqueue must be recomputed

		for (Finite_vertices_iterator v = m_dt.finite_vertices_begin();
				v != m_dt.finite_vertices_end(); ++v) {
			if (v->pinned())
				continue;
			v->relocated() = compute_relocation(v);
		}

		for (Finite_vertices_iterator v = m_dt.finite_vertices_begin();
				v != m_dt.finite_vertices_end(); ++v) {
			if (v->pinned())
				continue;
			if (v->point() == v->relocated())
				continue;

			Edge_list hull;
			m_dt.get_edges_from_star_minus_link(v, hull, false);
			bool ok = m_dt.is_in_kernel(v->relocated(), hull.begin(),
					hull.end());

			if (ok) {
				// do relocation
				FT norm_bef = m_dt.compute_cost_around_vertex(v).norm();
				relocate_one_vertex(v);
				FT norm_aft = m_dt.compute_cost_around_vertex(v).norm();

				// undo relocation
				if (norm_bef < norm_aft)
					relocate_one_vertex(v);
			}
		}

		std::cerr << yellow << "done" << white << " (" << yellow
				<< time_duration(timer) << white << " s)" << std::endl;
	}

	 /// \cond SKIP_IN_MANUAL
	Vector compute_gradient(Vertex_handle vertex) {
		Vector grad(0.0, 0.0);
		Edge_circulator ecirc = m_dt.incident_edges(vertex);
		Edge_circulator eend = ecirc;
		CGAL_For_all(ecirc, eend)
		{
			Edge edge = *ecirc;
			if (m_dt.source_vertex(edge) != vertex)
				edge = m_dt.twin_edge(edge);

			if (m_dt.get_plan(edge) == 0)
				grad = grad + compute_gradient_for_plan0(edge);
			else
				grad = grad + compute_gradient_for_plan1(edge);
		}
		return grad;
	}

	Point compute_relocation(Vertex_handle vertex) {
		FT coef = 0.0;
		Vector rhs(0.0, 0.0);

		Edge_circulator ecirc = m_dt.incident_edges(vertex);
		Edge_circulator eend = ecirc;
		CGAL_For_all(ecirc, eend)
		{
			Edge edge = *ecirc;
			if (m_dt.source_vertex(edge) != vertex)
				edge = m_dt.twin_edge(edge);

			if (m_dt.get_plan(edge) == 0)
				compute_relocation_for_plan0(edge, coef, rhs);
			else
				compute_relocation_for_plan1(edge, coef, rhs);
		}
		compute_relocation_for_vertex(vertex, coef, rhs);

		if (coef == 0.0)
			return vertex->point();
		return CGAL::ORIGIN + (rhs / coef);
	}

	void compute_relocation_for_vertex(Vertex_handle vertex, FT& coef,
			Vector& rhs) {
		Sample* sample = vertex->get_sample();
		if (sample) {
			const FT m = sample->mass();
			const Point& ps = sample->point();
			rhs = rhs + m * (ps - CGAL::ORIGIN);
			coef += m;
		}
	}

	Vector compute_gradient_for_plan0(const Edge& edge) {
		Edge twin = m_dt.twin_edge(edge);
		const Point& pa = m_dt.source_vertex(edge)->point();
		const Point& pb = m_dt.target_vertex(edge)->point();

		Sample_list samples;
		m_dt.collect_samples_from_edge(edge, samples);
		m_dt.collect_samples_from_edge(twin, samples);

		Vector grad(0.0, 0.0);
		Sample_list_const_iterator it;
		for (it = samples.begin(); it != samples.end(); ++it) {
			Sample* sample = *it;
			const FT m = sample->mass();
			const Point& ps = sample->point();

			FT Da = CGAL::squared_distance(ps, pa);
			FT Db = CGAL::squared_distance(ps, pb);
			if (Da < Db)
				grad = grad + m * (pa - ps);
		}
		return grad;
	}

	void compute_relocation_for_plan0(const Edge& edge, FT& coef, Vector& rhs) {
		Edge twin = m_dt.twin_edge(edge);
		const Point& pa = m_dt.source_vertex(edge)->point();
		const Point& pb = m_dt.target_vertex(edge)->point();

		Sample_list samples;
		m_dt.collect_samples_from_edge(edge, samples);
		m_dt.collect_samples_from_edge(twin, samples);

		Sample_list_const_iterator it;
		for (it = samples.begin(); it != samples.end(); ++it) {
			Sample* sample = *it;
			const FT m = sample->mass();
			const Point& ps = sample->point();

			FT Da = CGAL::squared_distance(ps, pa);
			FT Db = CGAL::squared_distance(ps, pb);

			if (Da < Db) {
				rhs = rhs + m * (ps - CGAL::ORIGIN);
				coef += m;
			}
		}
	}

	Vector compute_gradient_for_plan1(const Edge& edge) {
		FT M = m_dt.get_mass(edge);
		const Point& pa = m_dt.source_vertex(edge)->point();
		const Point& pb = m_dt.target_vertex(edge)->point();

		SQueue queue;
		m_dt.sort_samples_from_edge(edge, queue);

		//FT start = 0.0;
		Vector grad(0.0, 0.0);
		while (!queue.empty()) {
			PSample psample = queue.top();
			queue.pop();

			const FT m = psample.sample()->mass();
			const Point& ps = psample.sample()->point();

			// normal + tangnetial
			const FT coord = psample.priority();
			Point pf = CGAL::ORIGIN + (1.0 - coord) * (pa - CGAL::ORIGIN)
					+ coord * (pb - CGAL::ORIGIN);
			grad = grad + m * (1.0 - coord) * (pf - ps);

			/*
			 // only normal
			 FT bin = m/M;
			 FT center = start + 0.5*bin;
			 Point pc = CGAL::ORIGIN + (1.0-center)*(pa - CGAL::ORIGIN) + center*(pb - CGAL::ORIGIN);
			 start += bin;
			 grad = grad + m*(bin*bin/12.0)*(pa - pb);
			 grad = grad + m*(1.0-center)*(pc - pf);
			 */
		}
		return grad;
	}

	void compute_relocation_for_plan1(const Edge& edge, FT& coef, Vector& rhs) {
		FT M = m_dt.get_mass(edge);
		const Point& pb = m_dt.target_vertex(edge)->point();

		SQueue queue;
		m_dt.sort_samples_from_edge(edge, queue);

		//FT start = 0.0;
		while (!queue.empty()) {
			PSample psample = queue.top();
			queue.pop();

			const FT m = psample.sample()->mass();
			const Point& ps = psample.sample()->point();

			const FT coord = psample.priority();
			const FT one_minus_coord = 1.0 - coord;

			// normal + tangential
			coef += m * one_minus_coord * one_minus_coord;
			rhs =
					rhs
							+ m * one_minus_coord
									* ((ps - CGAL::ORIGIN)
											- coord * (pb - CGAL::ORIGIN));

			/*
			 // only normal
			 FT bin = m/M;
			 FT center = start + 0.5*bin;
			 FT one_minus_center = 1.0 - center;
			 start += bin;

			 coef += m*bin*bin/12.0;
			 rhs = rhs + m*(bin*bin/12.0)*(pb - CGAL::ORIGIN);

			 coef += m*one_minus_center*(coord - center);
			 rhs = rhs + m*one_minus_center*(coord - center)*(pb - CGAL::ORIGIN);
			 */
		}
	}

	void print_stats_debug() const
	{
	    int nb_solid = 0;
	    int nb_ghost = 0;
	    for (Finite_edges_iterator ei = m_dt.finite_edges_begin();
	    		ei != m_dt.finite_edges_end(); ++ei)
	    {
	        Edge edge = *ei;
	        if (m_dt.is_ghost(edge)) nb_ghost++;
	        else nb_solid++;
	    }

	    std::cerr << blue << "STATS" << white << std::endl;
	    std::cerr << "# vertices : " << m_dt.number_of_vertices()-4 << std::endl;
	    std::cerr << "# triangles: " << m_dt.number_of_faces() << std::endl;
	    std::cerr << "# edges: " << m_dt.tds().number_of_edges() << std::endl;
	    std::cerr << "# solid: " << nb_solid << std::endl;
	    std::cerr << "# ghost: " << nb_ghost << std::endl;
	}

	/// \endcond

	// RECONSTRUCTION //

		 /*!
		    This function must be called after initialization().
		    It computes a shape consisting of nv vertices, reconstructing the input
		    points.

		    \param nv The number of vertices which will be present in the output.
		  */
		void reconstruct_until(const unsigned nv) {
			double timer = clock();
			std::cerr << yellow << "reconstruct until " << white << nv << " V";

			unsigned N = nv + 4;
			unsigned performed = 0;
			while (m_dt.number_of_vertices() > N) {
				bool ok = decimate();
				if (!ok)
					break;
				performed++;
			}

			std::cerr << yellow << " done" << white << " (" << performed
					<< " iters, " << m_dt.number_of_vertices() - 4 << " V "
					<< yellow << time_duration(timer) << white << " s)"
					<< std::endl;
		}

		 /*!
			This function must be called after initialization().
			It computes a shape, reconstructing the input, by performing steps many
			edge contractions on the output simplex.

			\param steps The number of edge contractions performed by the algorithm.
		  */
		void reconstruct(const unsigned steps) {
			double timer = clock();
			std::cerr << yellow << "reconstruct " << steps << white;

			unsigned performed = 0;
			for (unsigned i = 0; i < steps; ++i) {
				bool ok = decimate();
				if (!ok)
					break;
				performed++;
			}

			std::cerr << yellow << " done" << white << " (" << performed << "/"
					<< steps << " iters, " << m_dt.number_of_vertices() - 4
					<< " V, " << yellow << time_duration(timer) << white << " s)"
					<< std::endl;
		}


};
}

#endif
