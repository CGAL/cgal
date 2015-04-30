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
#include <utility>      // std::pair

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/mem_fun.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/identity.hpp>
#include <boost/property_map/property_map.hpp>

namespace CGAL {


/*!
\ingroup PkgReconstructionSimplification2Classes

This class enables to reconstruct a 1-dimensional shape from a range of points with masses.
The algorithm computes a triangulation, ....   and performs a simplification of the triangulation
by performing edge contractions. The edges are either processed in the order imposed by 
a priority queue, or in an order based on random sampling. As the priority queue guarantees
a higher quality it is the default. The user can switch to the other method, for example
for an initial simplification round, by calling `set_random_sample_size()`.

\todo @@Pierre:  In the same way discuss the other parameters and run functions.

\todo `Gt` must at least be a model of `DelaunayTriangulationTraits_2`.  @@Pierre: If more functionalty is needed we should introduce a `ReconstructionSimplificationTraits_2`.

\tparam Gt a model of the concept `Kernel`.

\tparam PointMap a model of `ReadablePropertyMap` with value type `Gt::Point_2`

\tparam MassMap a model of `ReadablePropertyMap` with value type `Gt::FT`

 */
template<class Gt,
        class PointMap = First_of_pair_property_map  <std::pair<typename Gt::Point_2 , typename Gt::FT > >,
         class MassMap  = boost::static_property_map <typename Gt::Point_2 , typename Gt::FT > >
class Reconstruction_simplification_2 {
public:

    /// \name Types 
    /// @{
    /*!
        Number type.
    */
    typedef typename Gt::FT FT;

    /*!
        Point type.
    */
    typedef typename Gt::Point_2 Point;

	 /// \cond SKIP_IN_MANUAL
	/*!
		Vector type.
	*/
	typedef typename Gt::Vector_2 Vector;


	typedef typename std::pair<Point, FT> PointMassPair;
	typedef typename std::list<PointMassPair> PointMassList;


	/*!
	The Output simplex.
	*/
	typedef Reconstruction_triangulation_2<Gt> Triangulation;


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

    /// @}

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

    PointMap point_pmap;
	MassMap  mass_pmap;

	  /// \endcond

	public:

	  /// \name Initialization
	  /// @{

	/*!
             The constructor of the reconstruction simplification class
             for a given range of point-mass pairs.
             which already builds an initial simplex.

	     \tparam InputRange is a model of `Range` with forward iterators, 
               providing input points and point mass through the following two property maps.

	     \param input_range range of input data.
	     \param point_map A `ReadablePropertyMap` used to access the input points.

	     \param mass_map A `ReadablePropertyMap` used to access the input points' mass.
             \param sample_size If `sample_size != 0`, the size of the random sample that replaces a priority queue.
             \param use_flip If `true` the flipping procedure is used for the halfedge collapse.
             \param relocation The number of point relocations that are performed between two edge collapses.
             \param verbose controls how much console output is produced by the algorithm. The values are 0,1, or >1.

	*/
	template <class InputRange>
	Reconstruction_simplification_2(const InputRange& input_range,
                                        PointMap point_map = PointMap(),
                                        MassMap  mass_map = MassMap(1),
                                        std::size_t sample_size = 0,
                                        bool use_flip = true,
                                        std::size_t relocation = 0,
                                        std::size_t verbose = 0
                                        ) {


		point_pmap = point_map;
		mass_pmap  = mass_map;

		initialize_parameters();

		initialize(input_range.begin(), input_range.end());
	}



	  /// @}

    /// \name Settting Parameters 
  /// @{
	/*!
          If `sample_size == 0`, the edge collapse is done using a priority queue.
          \todo @@Pierre:  Tell what is a good value.
          \param sample_size If `sample_size != 0`, the size of the random sample that replaces the priority queue.
	*/
  void set_random_sample_size(std::size_t sample_size) {
		m_mchoice = mchoice;
	}

	/*!
		Determines how much console output the algorithm generates.
		If set to a value larger than 0
		details about the reconstruction process are written to `std::err`.

		\param verbose The verbosity level.
	*/
  void set_verbose(std::size_t verbose) {
		m_verbose = verbose;
	}


	 /// \cond SKIP_IN_MANUAL
	void set_alpha(const double alpha) {
		m_alpha = alpha;
	}
	 /// \endcond


	/*!
		The use_flip parameter determines whether the flipping procedure
		is used for the half-edge collapse. 
	 */
	void set_use_flip(const bool use_flip) {
		m_use_flip = use_flip;
	}


	 /// \cond SKIP_IN_MANUAL
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
	/// \endcond


	/*!
		Sets the number of point relocations
		that are performed between two edge collapses.
	*/
  void set_relocation(std::size_t relocation) {
		m_relocation = relocation;
	}

	/// \cond SKIP_IN_MANUAL
	unsigned get_relocation() const {
		return m_relocation;
	}
	/// \endcond


	/*!
          \todo @@Pierre: explain what relevance means 

		\param relevance The relevance level.
	 */
	void set_relevance(const double relevance) {
		m_ghost = relevance;
		m_dt.ghost_factor() = m_ghost;
	}


	/// \cond SKIP_IN_MANUAL
	double get_ghost() {
		return m_ghost;
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
	//without yet specifying the input in the constructor.
	template <class InputIterator>
	void initialize(InputIterator start_itr,
									InputIterator beyond_itr,
									PointMap point_map,
									MassMap  mass_map) {

		point_pmap = point_map;
		mass_pmap  = mass_map;

		initialize(start_itr, beyond_itr);

	}


	template <class InputIterator>
	void initialize(InputIterator start, InputIterator beyond) {

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





	template <class Vector>
	Vector random_vec(const double scale)
	{
	    double dx = -scale + (double(rand()) / double(RAND_MAX)) * 2* scale;
	    double dy = -scale + (double(rand()) / double(RAND_MAX)) * 2* scale;
	    return Vector(dx, dy);
	}

	void clear() {
		m_dt.clear();
		m_mindex.clear();
	}

	double time_duration(const double init) {
		return (clock() - init) / CLOCKS_PER_SEC;
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
	void init(Iterator begin, Iterator beyond) {
		double timer = clock();
		std::cerr << yellow << "init" << white << "...";

		int nb = m_dt.number_of_vertices();
		m_dt.infinite_vertex()->pinned() = true;
		for (Iterator it = begin; it != beyond; it++) {
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




	 /*!
	    Returns the number of vertices present in the reconstructed triangulation.
	  */
	int number_of_vertices() {
		return m_dt.number_of_vertices()-4 ;

	}

	 /*!
	    Returns the number of (solid) edges present in the reconstructed triangulation.
	  */
	int number_of_edges() {
		int nb_solid = 0;
		for (Finite_edges_iterator ei = m_dt.finite_edges_begin();
			    		ei != m_dt.finite_edges_end(); ++ei)
		{
			Edge edge = *ei;
			if (m_dt.is_ghost(edge)) continue;
			nb_solid++;
		}
		return nb_solid;
	}


	 /*!
	    Returns the cost of the (solid) edges present in the reconstructed triangulation.
	  */
	FT total_edge_cost() {
		FT total_cost = 0;
		for (Finite_edges_iterator ei = m_dt.finite_edges_begin();
				ei != m_dt.finite_edges_end(); ++ei) {
			Edge edge = *ei;
			if (m_dt.is_ghost(edge)) continue;

			total_cost += m_dt.get_cost(edge).finalize();
		}
		return total_cost;
	}

	/// \endcond


  /// \name Simplification 
  /// You can freely mix calls of the following functions. 
  /// @{
		 /*!
		    Computes a shape consisting of `np`  points, reconstructing the input
		    points.

		    \param np The number of points which will be present in the output.
		  */
  void run_until(std::size_t np) {
			double timer = clock();
			std::cerr << yellow << "reconstruct until " << white << np << " V";

			unsigned N = np + 4;
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
			Computes a shape, reconstructing the input, by performing `steps` many
			edge contractions on the output simplex.

			\param steps The number of edge contractions performed by the algorithm.
		  */
		void run(const unsigned steps) {
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


	/*!
	Since noise and missing data may prevent the reconstructed shape to
	have sharp corners well located, the algorithm offers the possibility to automatically
	relocate points after each edge contraction. The new location of the
	points is chosen such that the fitting of the output segments to the
	input points is improved. This is achieved by minimizing the normal component
	of the weighted \f$L_2 \f$ distance. The points then get relocated only if the
	underlying triangulation is still embeddable.
	  */
	void relocate_all_points() {
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

  /// @}

    /// \name Output
    /// @{

	/*!
	Writes the points and segments of the output simplex in an indexed format into output iterators.
        \tparam  PointOutputIterator An output iterator with value type `Point`.
        \tparam IndexOutputIterator An output iterator with value type `std::size_t`
        \tparam IndexPairOutputIterator An output iterator with value type `std::pair<std::size_t,std::size_t>`

	\param points the output iterator for all points
        \param isolated_points the output iterator for the indices of isolated points
        \param segments the output iterator for the pairs of indices of segments
	*/
  template <typename PointOutputIterator,
            typename IndexOutputIterator,
            typename IndexPairOutputIterator>
  void indexed_output(PointOutputIterator points,
                      IndexOutputIterator isolated_points,
                      IndexPairOutputIterator segments) {

		typedef typename Gt::Segment_2 Segment;
		std::vector<Point> isolated_points;
		std::vector<Segment> edges;

		extract_list_output(std::back_inserter(isolated_points), std::back_inserter(edges));


		//vertices_of_edges
		std::set<Point> edge_vertices;
		for (typename std::vector<Segment>::iterator it = edges.begin();
						it != edges.end(); it++) {

			Point a = (*it).source();
			Point b = (*it).target();

			edge_vertices.insert(a);
			edge_vertices.insert(b);
		}

		os << "OFF " << isolated_points.size() + edge_vertices.size() <<
				" 0 " << edges.size()  << std::endl;

		for (typename std::vector<Point>::iterator it = isolated_points.begin();
					it != isolated_points.end(); it++) {
			os << *it << std::endl;
		}

		for (typename std::set<Point>::iterator it = edge_vertices.begin();
				it != edge_vertices.end(); it++) {

			os << *it << std::endl;
		}

		for (int i = 0; i < isolated_points.size(); i++) {
			os << "1 " <<  i << std::endl;
		}

		for (typename std::vector<Segment>::iterator it = edges.begin();
				it != edges.end(); it++) {

			//save_one_edge(os, *it,edge_vertices);

			Point a = (*it).source();
			Point b = (*it).target();

			typename std::set<Point>::iterator it_a = edge_vertices.find(a);
			typename std::set<Point>::iterator it_b = edge_vertices.find(b);

			int pos_a = std::distance(edge_vertices.begin(), it_a);
			int pos_b = std::distance(edge_vertices.begin(), it_b);

			os << "2 "  << pos_a + isolated_points.size() << " "
					<< pos_b + isolated_points.size() << std::endl;



		}
	}

	 /// \cond SKIP_IN_MANUAL

	/*!
	 Returns the solid edges and vertices present after the reconstruction
	 process finished.

	\details It takes two output iterators, one for storing the
	isolated points and one for storing the edges of the reconstructed shape.


	\tparam PointOutputIterator The output iterator type for storing the isolated points

	\tparam SegmentOutputIterator The output iterator type for storing the edges as segments.
	 */
	template<class PointOutputIterator, class SegmentOutputIterator>
	void extract_list_output(PointOutputIterator v_it, SegmentOutputIterator e_it) {

		for (Vertex_iterator vi = m_dt.vertices_begin();
						vi != m_dt.vertices_end(); ++vi)
		{

			bool incident_edges_have_sample = false;
			typename Triangulation::Edge_circulator start = m_dt.incident_edges(vi);
			typename Triangulation::Edge_circulator cur   = start;

			do {
				if (!m_dt.is_ghost(*cur)) {
					incident_edges_have_sample = true;
					break;
				}
				++cur;
			} while (cur != start);

			if (!incident_edges_have_sample) {
				if ((*vi).has_sample_assigned()) {
					Point p = (*vi).point();
					*v_it = p;
					v_it++;
				}
			}
		}

		for (Finite_edges_iterator ei = m_dt.finite_edges_begin(); ei != m_dt.finite_edges_end(); ++ei)
		{
			Edge edge = *ei;
			if (m_dt.is_ghost(edge))
				continue;

	        int index = edge.second;
	        Vertex_handle source = edge.first->vertex( (index+1)%3 );
	        Vertex_handle target = edge.first->vertex( (index+2)%3 );

	        typename Gt::Segment_2  s(source->point(), target->point());
			*e_it = s;
			e_it++;
		}


	}
	 /// \endcond



	 /// \cond SKIP_IN_MANUAL
	void extract_tds_output(Triangulation& rt2) {
		rt2 = m_dt;
		//mark vertices
		for (Vertex_iterator vi = rt2.vertices_begin();
			  vi != rt2.vertices_end(); ++vi)
		{

			bool incident_edges_have_sample = false;
			typename Triangulation::Edge_circulator start = rt2.incident_edges(vi);
			typename Triangulation::Edge_circulator cur = start;

			do {
			  if (!rt2.is_ghost(*cur)) {
				  incident_edges_have_sample = true;
				  break;
			  }
			  ++cur;
			} while (cur != start);

			if (!incident_edges_have_sample) {
			  if ((*vi).has_sample_assigned())
				  (*vi).set_relevance(1);
			}
		}


		//mark edges
		for (Finite_edges_iterator ei = rt2.finite_edges_begin(); ei != rt2.finite_edges_end(); ++ei)
		{
			Edge edge = *ei;
			FT relevance = 0;
			if (!rt2.is_ghost(edge)) {
				relevance = rt2.get_edge_relevance(edge); // >= 0
			}
			edge.first->relevance(edge.second) = relevance;
		}
	}


  /// \endcond
    /// @}

};
}

#endif
