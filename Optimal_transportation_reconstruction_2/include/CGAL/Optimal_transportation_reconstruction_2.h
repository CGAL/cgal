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
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Fernando de Goes, Pierre Alliez, Ivo Vigan, Clément Jamin

#ifndef CGAL_OPTIMAL_TRANSPORTATION_RECONSTRUCTION_2_H_
#define CGAL_OPTIMAL_TRANSPORTATION_RECONSTRUCTION_2_H_

#include <CGAL/license/Optimal_transportation_reconstruction_2.h>


#include <CGAL/OTR_2/Reconstruction_triangulation_2.h>
#include <CGAL/OTR_2/Reconstruction_edge_2.h>

#include <CGAL/property_map.h>
#include <CGAL/Real_timer.h>
#include <CGAL/Random.h>

#include <iterator>
#include <iostream>
#include <vector>
#include <list>
#include <algorithm>
#include <utility>      // std::pair

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/mem_fun.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/identity.hpp>
#include <CGAL/boost/iterator/transform_iterator.hpp>
#include <boost/type_traits/is_float.hpp>

namespace CGAL {


/*!
\ingroup PkgOptimalTransportationReconstruction2Classes

This class provides a means to reconstruct a 1-dimensional shape from a set of 2D points with masses.
The algorithm computes an initial 2D Delaunay triangulation from the input points, 
and performs a simplification of the triangulation by performing half edge collapses, edge flips and vertex relocations.

The edges are either processed in the order imposed by an priority queue, or
in an order based on random selection of edge collapse operators.
As the exhaustive priority queue guarantees a higher quality it is the default.
The user can switch to the other method, for example for an initial
simplification round, by calling `set_random_sample_size()`.

By default edge flip operators are applied to ensure that every edge of the 
triangulation are candidate to be collapsed, while preserving a valid embedding
of the triangulation. This option can be disabled by calling
\link set_use_flip() `set_use_flip(false)`\endlink to reduce the running times.

By default the vertices are not relocated after each half edge collapse.
This option can be changed by setting the number of vertex relocation steps
performed between two edge collapse operators.

The simplification is performed by calling either 
\link run_until() `run_until(n)`\endlink or \link run() `run(steps)`\endlink.
The former simplifies the triangulation until n points remain, while the latter
stops after `steps` edge collapse operators have been performed.
Furthermore, we can relocate the vertices by calling `relocate_all_points()`.

\tparam Traits a model of the concept `OptimalTransportationReconstructionTraits_2`.

\tparam PointPMap a model of `ReadablePropertyMap` with value type `Traits::Point_2`.
        Defaults to <a href="http://www.boost.org/doc/libs/release/libs/property_map/doc/identity_property_map.html">`boost::typed_identity_property_map<Traits::Point_2>`</a> 
        (for the case the input is points without mass).

\tparam MassPMap a model of `ReadablePropertyMap` with value type `Traits::FT`
        Defaults to <a href="http://www.boost.org/doc/libs/release/libs/property_map/doc/static_property_map.html">`boost::static_property_map<Traits::FT>`</a> 
        (for the case the input is points without mass).

 */
template<
  class Traits,
  class PointPMap = boost::typed_identity_property_map <typename Traits::Point_2>,
  class MassPMap  = boost::static_property_map <typename Traits::FT> >
class Optimal_transportation_reconstruction_2
{
public:

  /// \name Types
  /// @{
  /*!
        Number type.
   */
  typedef typename Traits::FT FT;

  /*!
        Point type.
   */
  typedef typename Traits::Point_2 Point;

  /*!
        Segment type.
  */
  typedef typename Traits::Segment_2 Segment;

  /// \cond SKIP_IN_MANUAL
  /*!
        Vector type.
   */
  typedef typename Traits::Vector_2 Vector;

  typedef typename std::pair<Point, FT> PointMassPair;
  typedef typename std::vector<PointMassPair> PointMassList;


  /*!
    The Output simplex.
   */
  typedef OTR_2::Reconstruction_triangulation_2<Traits>  Triangulation;

  typedef typename Triangulation::Vertex                Vertex;
  typedef typename Triangulation::Vertex_handle         Vertex_handle;
  typedef typename Triangulation::Vertex_iterator       Vertex_iterator;
  typedef typename Triangulation::Vertex_circulator     Vertex_circulator;
  typedef typename Triangulation::Finite_vertices_iterator
                                                        Finite_vertices_iterator;

  typedef typename Triangulation::Edge                  Edge;
  typedef typename Triangulation::Edge_circulator       Edge_circulator;
  typedef typename Triangulation::Finite_edges_iterator Finite_edges_iterator;

  typedef typename Triangulation::Face_handle           Face_handle;
  typedef typename Triangulation::Face_circulator       Face_circulator;
  typedef typename Triangulation::Finite_faces_iterator Finite_faces_iterator;

  typedef typename Triangulation::Vertex_handle_map     Vertex_handle_map;
  typedef typename Triangulation::Face_handle_map       Face_handle_map;

  typedef typename Triangulation::Vertex_handle_set     Vertex_handle_set;
  typedef typename Triangulation::Edge_set              Edge_set;

  typedef typename Triangulation::Edge_vector           Edge_vector;
  typedef std::list<Edge>                               Edge_list;

  typedef typename Triangulation::Cost_                 Cost_;
  typedef typename Triangulation::Sample_               Sample_;
  typedef typename Triangulation::Sample_vector         Sample_vector;
  typedef typename Triangulation::Sample_vector_const_iterator
                                                        Sample_vector_const_iterator;

  typedef typename Triangulation::PSample               PSample;
  typedef typename Triangulation::SQueue                SQueue;

  typedef typename Triangulation::Rec_edge_2            Rec_edge_2;

  typedef typename Triangulation::MultiIndex            MultiIndex;

  /// @}

protected:
  Triangulation m_dt;
  Traits const& m_traits;
  MultiIndex m_mindex;
  int m_ignore;
  int m_verbose;
  std::size_t m_mchoice;  // # Edges
  bool m_use_flip;
  FT m_alpha; // [0, 1]
  FT m_ghost; // ghost vs solid
  unsigned int m_relocation; // # relocations
  FT m_tolerance;

  PointPMap point_pmap;
  MassPMap  mass_pmap;

  /// \endcond

public:

  /// \name Initialization
  /// @{

  /*!
  Constructor of the optimal transportation reconstruction class. 
  It builds an initial simplicial complex 
  for a given range of point-mass pairs.

  \tparam InputRange is a model of `Range` with forward iterators, 
          providing input points and point masses through the 
          `PointPMap` and `MassPMap` property maps.

  \param input_range  Range of input data.
  \param point_map    A `ReadablePropertyMap` used to access the input points.
  \param mass_map     A `ReadablePropertyMap` used to access the input 
                      points' masses.
  \param sample_size  If `sample_size != 0`, the size of the random sample 
                      which replaces the exhaustive priority queue.
  \param use_flip     If `true` the edge flipping procedure is used to ensure 
                      that every edge can be made collapsible.
  \param relocation   The number of point relocations that are performed 
                      between two edge collapses.
  \param verbose      Controls how much console output is produced by 
                      the algorithm. The values are 0, 1, or > 1.
  \param traits       The traits class.
   */
  template <class InputRange>
  Optimal_transportation_reconstruction_2(
    const InputRange& input_range,
    PointPMap point_map = PointPMap(),
    MassPMap  mass_map = MassPMap(1),
    std::size_t sample_size = 0,
    bool use_flip = true,
    unsigned int relocation = 2,
    int verbose = 0,
    Traits traits = Traits())
  : m_dt(traits),
    m_traits(m_dt.geom_traits()),
    m_ignore(0), 
    m_verbose(verbose),
    m_mchoice(sample_size),
    m_use_flip(use_flip),
    m_alpha(0.5),
    m_ghost(1.0),
    m_relocation(relocation),
    m_tolerance (FT(-1.)),
    point_pmap(point_map),
    mass_pmap(mass_map)
  {
    initialize(input_range.begin(), input_range.end());
  }

  /// @}

  /// \name Settting Parameters
  /// @{
  /*!
          If `sample_size == 0`, the simplification is performed using an exhaustive priority queue.
          If `sample_size` is stricly positive the simplification is performed using a
          multiple choice approach, ie, a best-choice selection in a random sample of 
          edge collapse operators, of size `sample_size`. A typical value for the sample
          size is 15, but this value must be enlarged when targeting a very coarse simplification.
          \param sample_size If `sample_size != 0`, the size of the random sample replaces the priority queue.
   */
  void set_random_sample_size(std::size_t sample_size) {
    m_mchoice = sample_size;
  }

  /*!
        Determines how much console output the algorithm generates.
        If set to a value larger than 0
        details about the reconstruction process are written to `std::cerr`.

        \param verbose The verbosity level.
   */
  void set_verbose(int verbose) {
    m_verbose = verbose;
  }


  

  /*!
        The use_flip parameter determines whether the edge flipping procedure
        is used for the half-edge collapse. 
   */
  void set_use_flip(const bool use_flip) {
    m_use_flip = use_flip;
  }


  /*!
        Sets the number of vertex relocations
        that are performed between two edge collapses.
   */
  void set_relocation(unsigned int relocation) {
    m_relocation = relocation;
  }

  /// \cond SKIP_IN_MANUAL
  unsigned int relocation() const {
    return m_relocation;
  }
  /// \endcond


  /*!
  \param relevance The relevance threshold used for filtering the edges.
  An edge is relevant from the approximation point of view
  if it is long, covers a large mass (or equivalently the
  number of points when all masses are equal), and has a
  small transport cost. This notion is defined as 
  \f$ m(e) * |e|^2 / cost(e) \f$, where \f$ m(e) \f$ 
  denotes the mass of the points approximated by the edge,
  \f$ |e| \f$ denotes the edge length and \f$ cost(e) \f$
  its approximation error.
  As the cost is defined by mass time squared distance the
  relevance is unitless.

  The default value is 1, so that all edges receiving some mass
  are considered relevant.
  Setting a large relevance value is used to get robustness to a
  large amount of outliers.
   */
  void set_relevance(const FT relevance) {
    m_ghost = relevance;
    m_dt.ghost_factor() = m_ghost;
  }


  /// \cond SKIP_IN_MANUAL
  FT ghost() {
    return m_ghost;
  }

  FT tolerance() const { return m_tolerance; }

  /// @}

  /// \cond SKIP_IN_MANUAL

  Optimal_transportation_reconstruction_2()
  : m_traits(m_dt.geom_traits())
  {
    initialize_parameters();
  }


  ~Optimal_transportation_reconstruction_2() {
    clear();
  }

  void initialize_parameters() {
    m_verbose = 0;
    m_mchoice = 0;
    m_use_flip = true;
    m_alpha = FT(0.5);
    m_ghost = FT(1);
    m_relocation = 0;

    m_ignore = 0;
  }

  //Function if one wants to create a Optimal_transportation_reconstruction_2
  //without yet specifying the input in the constructor.
  template <class InputIterator>
  void initialize(
    InputIterator start_itr,
    InputIterator beyond_itr,
    PointPMap point_map,
    MassPMap  mass_map)
  {
    point_pmap = point_map;
    mass_pmap  = mass_map;

    initialize(start_itr, beyond_itr);
  }


  template <class InputIterator>
  void initialize(InputIterator start, InputIterator beyond) {

    clear();
    Property_map_to_unary_function<PointPMap> get_point(point_pmap);

    Bbox_2 bbox = bbox_2(
      boost::make_transform_iterator(start,get_point),
      boost::make_transform_iterator(beyond,get_point));

    insert_loose_bbox(bbox);
    init(start, beyond);

    std::vector<Sample_*> m_samples;
    for (InputIterator it = start; it != beyond; it++) {
      Point point = get(point_pmap, *it);
      FT    mass  = get( mass_pmap, *it);
      Sample_* s = new Sample_(point, mass);
      m_samples.push_back(s);
    }
    assign_samples(m_samples.begin(), m_samples.end());
  }

  template <class InputIterator>
  void initialize_with_custom_vertices(InputIterator samples_start,
                                       InputIterator samples_beyond,
                                       InputIterator vertices_start,
                                       InputIterator vertices_beyond,
                                       PointPMap point_map,
                                       MassPMap  mass_map) {
    point_pmap = point_map;
    mass_pmap  = mass_map;
    clear();
    Property_map_to_unary_function<PointPMap> get_point(point_pmap);

    Bbox_2 bbox = bbox_2(
      boost::make_transform_iterator(samples_start,get_point),
      boost::make_transform_iterator(samples_beyond,get_point));

    insert_loose_bbox(bbox);
    init(vertices_start, vertices_beyond);

    std::vector<Sample_*> m_samples;
    for (InputIterator it = samples_start; it != samples_beyond; it++) {
#ifdef CGAL_USE_PROPERTY_MAPS_API_V1
      Point point = get(point_pmap, it);
      FT    mass  = get( mass_pmap, it);
#else
      Point point = get(point_pmap, *it);
      FT    mass  = get( mass_pmap, *it);
#endif
      Sample_* s = new Sample_(point, mass);
      m_samples.push_back(s);
    }
    assign_samples(m_samples.begin(), m_samples.end());
  }


  template <class Vector>
  Vector random_vec(const FT scale) const
  {
    FT dx = -scale + get_default_random().get_double() * 2* scale;
    FT dy = -scale + get_default_random().get_double() * 2* scale;
    return m_traits.construct_vector_2_object()(dx, dy);
  }

  void clear() {
    Sample_vector samples;
    m_dt.collect_all_samples(samples);
    // Deallocate samples
    for (Sample_vector_const_iterator s_it = samples.begin();
        s_it != samples.end(); ++s_it)
    {
      delete *s_it;
    }
  }


  // INIT //
  void insert_loose_bbox(const Bbox_2& bbox) {
    CGAL::Real_timer timer;
    if (m_verbose > 0)
      std::cerr << "insert loose bbox...";

    double dl = (std::max)((bbox.xmax()-bbox.xmin()) / 2.,
                           (bbox.ymax()-bbox.ymin()) / 2.);

    timer.start();
    int nb = static_cast<int>(m_dt.number_of_vertices());
    typename Traits::Construct_point_2 point_2
      = m_traits.construct_point_2_object();
    insert_point(point_2(bbox.xmin()-dl, bbox.ymin()-dl), true, nb++);
    insert_point(point_2(bbox.xmin()-dl, bbox.ymax()+dl), true, nb++);
    insert_point(point_2(bbox.xmax()+dl, bbox.ymax()+dl), true, nb++);
    insert_point(point_2(bbox.xmax()+dl, bbox.ymin()-dl), true, nb++);

    if (m_verbose > 0)
      std::cerr << "done (" << nb << " vertices, "
                << timer.time() << " s)" << std::endl;
  }

  template<class Iterator>  // value_type = Point*
  void init(Iterator begin, Iterator beyond) {
    CGAL::Real_timer timer;
    if (m_verbose > 0)
      std::cerr << "init...";

    timer.start();
    int nb = static_cast<int>(m_dt.number_of_vertices());
    m_dt.infinite_vertex()->pinned() = true;
    for (Iterator it = begin; it != beyond; it++) {
      Point point = get(point_pmap, *it);
      insert_point(point, false, nb++);
    }

    if (m_verbose > 0)
      std::cerr << "done (" << nb << " vertices, "
                << timer.time() << " s)"
                << std::endl;
  }

private:
  Vertex_handle insert_point(
    const Point& point, const bool pinned, const int id) 
  {
    Vertex_handle v = m_dt.insert(point);
    v->pinned() = pinned;
    v->id() = id;
    return v;
  }
public:

  // ASSIGNMENT //

  void cleanup_assignments() {
    m_dt.cleanup_assignments();
  }

  template<class Iterator>  // value_type = Sample_*
  void assign_samples(Iterator begin, Iterator end) {
    CGAL::Real_timer timer;
    if (m_verbose > 0)
      std::cerr << "assign samples...";

    timer.start();
    m_dt.assign_samples(begin, end);
    m_dt.reset_all_costs();

    if (m_verbose > 0)
      std::cerr << "done (" << timer.time() << " s)" << std::endl;
  }

  void reassign_samples() {
    Sample_vector samples;
    m_dt.collect_all_samples(samples);
    m_dt.cleanup_assignments();
    m_dt.assign_samples(samples.begin(), samples.end());
    m_dt.reset_all_costs();
  }

  void reassign_samples_around_vertex(Vertex_handle vertex) {
    Sample_vector samples;
    m_dt.collect_samples_from_vertex(vertex, samples, true);
    m_dt.assign_samples(samples.begin(), samples.end());

    Edge_vector hull;
    m_dt.get_edges_from_star_minus_link(vertex, hull, true);
    update_cost(hull.begin(), hull.end());
  }


  bool do_collapse(Edge edge) {
    Vertex_handle s = m_dt.source_vertex(edge);
    Vertex_handle t = m_dt.target_vertex(edge);

    if (m_verbose > 1) {
      std::cerr << std::endl << "do collapse ("
          << s->id() << "->" << t->id() << ") ... " << std::endl;
    }

    Sample_vector samples;
    m_dt.collect_samples_from_vertex(s, samples, true);

    Edge_vector hull;
    m_dt.get_edges_from_star_minus_link(s, hull, true);

    if (m_mchoice == 0)
      remove_stencil_from_pqueue(hull.begin(), hull.end());

    if (m_use_flip)
      m_dt.make_collapsible(edge, hull.begin(), hull.end(), m_verbose);

    // debug test
    bool ok = m_dt.check_kernel_test(edge);
    if (!ok) {
      if (m_verbose > 1)
        std::cerr << "do_collapse: kernel test failed: " << std::endl;
      return false;
    }
    //

    m_dt.collapse(edge, m_verbose);

    m_dt.assign_samples(samples.begin(), samples.end());

    update_cost(hull.begin(), hull.end());

    if (m_mchoice == 0)
      push_stencil_to_pqueue(hull.begin(), hull.end());

    for (unsigned int i = 0; i < m_relocation; ++i) {
      relocate_one_ring(hull.begin(), hull.end());
    }

    if (m_verbose > 1) {
      std::cerr << "done" << std::endl;
    }

    return true;
  }

  bool simulate_collapse(const Edge& edge, Cost_& cost) {
    bool ok;
    Vertex_handle s = m_dt.source_vertex(edge);
    Vertex_handle t = m_dt.target_vertex(edge);

    if (m_verbose > 1) {
      std::cerr << "simulate collapse ("
        << s->id() << "->" << t->id() << ") ... " << std::endl;
    }

    Triangulation copy;
    Edge copy_edge = copy_star(edge, copy);
    Vertex_handle copy_source = copy.source_vertex(copy_edge);

    if (m_use_flip) {
      Edge_vector copy_hull;
      copy.get_edges_from_star_minus_link(copy_source, copy_hull, true);
      ok = copy.make_collapsible(copy_edge, copy_hull.begin(),
          copy_hull.end(), m_verbose);
      if (!ok) {
        // std::cerr << "simulation: failed (make collapsible)" << std::endl;
        return false;
      }
    }

    ok = copy.check_kernel_test(copy_edge);
    if (!ok) {
      std::cerr << "simulation: failed (kernel test)" << std::endl;
      return false;
    }

    copy.collapse(copy_edge, m_verbose);

    Sample_vector samples;
    m_dt.collect_samples_from_vertex(s, samples, false);

    backup_samples(samples.begin(), samples.end());
    copy.assign_samples_brute_force(samples.begin(), samples.end());
    copy.reset_all_costs();
    cost = copy.compute_total_cost();
    cost.set_total_weight (samples);
    restore_samples(samples.begin(), samples.end());

    if (m_verbose > 1) {
      std::cerr << "done" << std::endl;
    }

    return true;
  }

  template<class Iterator> // value_type = Sample_*
  void backup_samples(Iterator begin, Iterator end) const {
    for (Iterator it = begin; it != end; ++it) {
      Sample_* sample = *it;
      sample->backup();
    }
  }

  template<class Iterator> // value_type = Sample_*
  void restore_samples(Iterator begin, Iterator end) const {
    for (Iterator it = begin; it != end; ++it) {
      Sample_* sample = *it;
      sample->restore();
    }
  }

  // PEDGE //

  bool decimate() {
    bool ok;
    Rec_edge_2 pedge;
    ok = pick_edge(m_mchoice, pedge);
    if (!ok)
      return false;

    ok = do_collapse(pedge.edge());
    if (!ok)
      return false;
    return true;
  }

  bool is_above_tolerance (const Rec_edge_2& pedge)
  {
    if (m_tolerance == (FT)(-1.))
      return false;
    FT cost = CGAL::approximate_sqrt (pedge.after() / pedge.total_weight());
    return cost > m_tolerance;
  }

  bool create_pedge(const Edge& edge, Rec_edge_2& pedge) {
    Cost_ after_cost;
    bool ok = simulate_collapse(edge, after_cost);
    if (!ok)
      return false;

    Vertex_handle source = m_dt.source_vertex(edge);
    Cost_ before_cost = m_dt.compute_cost_around_vertex(source);

    FT before = before_cost.finalize(m_alpha);
    FT after = after_cost.finalize(m_alpha);
    pedge = Rec_edge_2(edge, before, after, after_cost.total_weight());

    if (is_above_tolerance (pedge))
      return false;
    
    return true;
  }


  // COST //

  void init_cost() {
    m_dt.reset_all_costs();
  }

  template<class Iterator> // value_type = Edge
  void update_cost(Iterator begin, Iterator end) {
    Edge_vector edges;
    collect_cost_stencil(m_dt, begin, end, edges);

    typename Edge_vector::iterator ei;
    for (ei = edges.begin(); ei != edges.end(); ++ei) {
      Edge edge = *ei;
      m_dt.update_cost(edge);
    }
  }

  template<class Iterator> // value_type = Edge
  void collect_cost_stencil(
    const Triangulation& mesh, Iterator begin, Iterator end, 
    Edge_vector& edges) const
  {
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

  bool pick_edge(std::size_t nb, Rec_edge_2& best_pedge) {
    if (m_dt.number_of_faces() < 2)
      return false;

    std::size_t ne = 2 * m_dt.tds().number_of_edges();
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

  bool pick_edge_from_pqueue(Rec_edge_2& best_pedge) {
    if (m_mindex.empty())
      populate_pqueue();
    if (m_mindex.empty())
      return false;
    best_pedge = *(m_mindex.template get<1>()).begin();
    (m_mindex.template get<0>()).erase(best_pedge);
    return true;
  }

  bool pick_edge_brute_force(Rec_edge_2& best_pedge) {
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

  bool pick_edge_randomly(std::size_t nb, Rec_edge_2& best_pedge) {
    MultiIndex mindex;
    for (std::size_t i = 0; i < nb; ++i) {
      Rec_edge_2 pedge;
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

    Rec_edge_2 pedge;
    bool ok = create_pedge(edge, pedge);
    if (!ok)
      return false;
    mindex.insert(pedge);
    return true;
  }



  bool random_pedge(Rec_edge_2& pedge) {
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
  void remove_stencil_from_pqueue(Iterator begin, Iterator end) 
  {
    if (m_mindex.empty())
      return;

    Edge_vector edges;
    collect_pqueue_stencil(m_dt, begin, end, edges);

    typename Edge_vector::const_iterator ei;
    for (ei = edges.begin(); ei != edges.end(); ++ei) {
      Edge edge = *ei;
      (m_mindex.template get<0>()).erase(Rec_edge_2(edge));
    }
  }

  template<class Iterator> // value_type = Edge
  void push_stencil_to_pqueue(Iterator begin, Iterator end) {
    Edge_vector edges;
    collect_pqueue_stencil(m_dt, begin, end, edges);

    typename Edge_vector::const_iterator ei;
    for (ei = edges.begin(); ei != edges.end(); ++ei) {
      Edge edge = *ei;
      push_to_mindex(edge, m_mindex);
    }
  }

  template<class Iterator> // value_type = Edge
  void collect_pqueue_stencil(
    const Triangulation& mesh, Iterator begin, Iterator end, 
    Edge_vector& edges) const
  {
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
    CGAL_assertion(s != m_dt.infinite_vertex() );

    Vertex_handle cs = copy.tds().create_vertex();
    cvmap[s] = copy_vertex(s, cs);

    Vertex_circulator vcirc = m_dt.incident_vertices(s);
    Vertex_circulator vend = vcirc;
    CGAL_For_all(vcirc, vend)
    {
      Vertex_handle v = vcirc;
      CGAL_assertion(v!=m_dt.infinite_vertex());
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
      copy_neighbors(f, s, cfmap);
    }

    // make copy homeomorphic to S^2
    close_copy_mesh(cs, copy);

    // copy samples surrounding star
    copy_samples(s, cs, cfmap, copy);

    // get copy of edge
    Edge copy_edge = get_copy_edge(edge, cvmap, cfmap);
    return copy_edge;
  }

  Vertex_handle copy_vertex(Vertex_handle v0, Vertex_handle v1) const {
    v1->id() = v0->id();
    v1->set_point(v0->point());
    v1->pinned() = v0->pinned();
    v1->set_sample(v0->sample());
    return v1;
  }

  Face_handle copy_face(
    Face_handle f0, Face_handle f1, Vertex_handle_map& vmap) const 
  {
    for (unsigned i = 0; i < 3; ++i) {
      Vertex_handle v0i = f0->vertex(i);
      Vertex_handle v1i = vmap[v0i];
      f1->set_vertex(i, v1i);
      v1i->set_face(f1);
    }
    return f1;
  }

  void copy_neighbors(
    Face_handle f, Vertex_handle v,
    Face_handle_map& fmap) const
  {
    int i = f->index(v);
    Face_handle cf = fmap[f];

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

  void close_copy_mesh(Vertex_handle vertex, Triangulation& copy) const {
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

  void copy_samples(
    Vertex_handle vertex, Vertex_handle copy_vertex,
    Face_handle_map& fmap, Triangulation& copy) const
  {
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

      Sample_vector samples;
      m_dt.collect_samples_from_edge(twin, samples);
      copy_twin.first->samples(copy_twin.second) = samples;
    }
    copy_vertex->set_sample(NULL);
  }

  Edge get_copy_edge(
    const Edge& edge, Vertex_handle_map& vmap, Face_handle_map& fmap) const 
  {
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

      Edge_vector hull;
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
  Vector compute_gradient(Vertex_handle vertex) const {
    Vector grad = m_traits.construct_vector_2_object()(FT(0), FT(0));
    Edge_circulator ecirc = m_dt.incident_edges(vertex);
    Edge_circulator eend = ecirc;
    CGAL_For_all(ecirc, eend)
    {
      Edge edge = *ecirc;
      if (m_dt.source_vertex(edge) != vertex)
        edge = m_dt.twin_edge(edge);

      if (m_dt.get_plan(edge) == 0)
        grad = m_traits.construct_sum_of_vectors_2_object()(
          grad, compute_gradient_for_plan0(edge));
      else
        grad = m_traits.construct_sum_of_vectors_2_object()(
          grad, compute_gradient_for_plan1(edge));
    }
    return grad;
  }

  // If the underlying number type used is not a floating point base
  // number type (like a multiprecision), the coordinates of the points
  // will increase a lot due to the relocation step. These functions
  // simply turn a relocated point to a rounded to double version.
  void relocate_on_the_double_grid(Point&, boost::true_type) const
  {}
  void relocate_on_the_double_grid(Point& p, boost::false_type) const
  {
    double x=to_double(m_traits.compute_x_2_object()(p));
    double y=to_double(m_traits.compute_y_2_object()(p));
    p=m_traits.construct_point_2_object()(FT(x),FT(y));
  }
  void relocate_on_the_double_grid(Point& p) const
  {
    relocate_on_the_double_grid(p,
      typename boost::is_float<typename Traits::FT>::type());
  }

  Point compute_relocation(Vertex_handle vertex) const {
    FT coef = FT(0);
    Vector rhs = m_traits.construct_vector_2_object()(FT(0), FT(0));

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

    if (coef == FT(0))
      return vertex->point();

    Point res = m_traits.construct_translated_point_2_object()(
      CGAL::ORIGIN,
      m_traits.construct_scaled_vector_2_object()(rhs, FT(1) / coef));
    relocate_on_the_double_grid(res);
    return res;
  }

  void compute_relocation_for_vertex(
    Vertex_handle vertex, FT& coef, Vector& rhs) const
  {
    Sample_* sample = vertex->sample();
    if (sample) {
      const FT m = sample->mass();
      const Point& ps = sample->point();
      rhs = m_traits.construct_sum_of_vectors_2_object()(rhs, 
        m_traits.construct_scaled_vector_2_object()(
          m_traits.construct_vector_2_object()(CGAL::ORIGIN, ps), m));
      coef += m;
    }
  }

  Vector compute_gradient_for_plan0(const Edge& edge) const {
    Edge twin = m_dt.twin_edge(edge);
    const Point& pa = m_dt.source_vertex(edge)->point();
    const Point& pb = m_dt.target_vertex(edge)->point();

    Sample_vector samples;
    m_dt.collect_samples_from_edge(edge, samples);
    m_dt.collect_samples_from_edge(twin, samples);

    Vector grad = m_traits.construct_vector_2_object()(FT(0), FT(0));
    Sample_vector_const_iterator it;
    for (it = samples.begin(); it != samples.end(); ++it) {
      Sample_* sample = *it;
      const FT m = sample->mass();
      const Point& ps = sample->point();

      FT Da = m_traits.compute_squared_distance_2_object()(ps, pa);
      FT Db = m_traits.compute_squared_distance_2_object()(ps, pb);
      if (Da < Db)
        grad = m_traits.construct_sum_of_vectors_2_object()(
          grad,
          m_traits.construct_scaled_vector_2_object()(
            m_traits.construct_vector_2_object()(ps, pa), m));
    }
    return grad;
  }

  void compute_relocation_for_plan0(
    const Edge& edge, FT& coef, Vector& rhs) const 
  {
    Edge twin = m_dt.twin_edge(edge);
    const Point& pa = m_dt.source_vertex(edge)->point();
    const Point& pb = m_dt.target_vertex(edge)->point();

    Sample_vector samples;
    m_dt.collect_samples_from_edge(edge, samples);
    m_dt.collect_samples_from_edge(twin, samples);

    Sample_vector_const_iterator it;
    for (it = samples.begin(); it != samples.end(); ++it) {
      Sample_* sample = *it;
      const FT m = sample->mass();
      const Point& ps = sample->point();

      FT Da = m_traits.compute_squared_distance_2_object()(ps, pa);
      FT Db = m_traits.compute_squared_distance_2_object()(ps, pb);

      if (Da < Db) {
        rhs = m_traits.construct_sum_of_vectors_2_object()(rhs, 
          m_traits.construct_scaled_vector_2_object()(
            m_traits.construct_vector_2_object()(CGAL::ORIGIN, ps), m));
        coef += m;
      }
    }
  }

  Vector compute_gradient_for_plan1(const Edge& edge) const {
    //FT M = m_dt.get_mass(edge);
    const Point& pa = m_dt.source_vertex(edge)->point();
    const Point& pb = m_dt.target_vertex(edge)->point();

    SQueue queue;
    m_dt.sort_samples_from_edge(edge, queue);

    //FT start = FT(0);
    Vector grad = m_traits.construct_vector_2_object()(FT(0), FT(0));
    while (!queue.empty()) {
      PSample psample = queue.top();
      queue.pop();

      const FT m = psample.sample()->mass();
      const Point& ps = psample.sample()->point();

      // normal + tangnetial
      const FT coord = psample.priority();
      Point pf = m_traits.construct_translated_point_2_object()(
        CGAL::ORIGIN,
        m_traits.construct_sum_of_vectors_2_object()(
          m_traits.construct_scaled_vector_2_object()(
            m_traits.construct_vector_2_object()(CGAL::ORIGIN, pa), 
            1.0 - coord),
          m_traits.construct_scaled_vector_2_object()(
            m_traits.construct_vector_2_object()(CGAL::ORIGIN, pb), 
            coord)));
      grad = m_traits.construct_sum_of_vectors_2_object()(
        grad,
        m_traits.construct_scaled_vector_2_object()(
          m_traits.construct_vector_2_object()(ps, pf), m * (1.0 - coord)));

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

  void compute_relocation_for_plan1(
    const Edge& edge, FT& coef, Vector& rhs) const 
  {
    //FT M = m_dt.get_mass(edge);
    const Point& pb = m_dt.target_vertex(edge)->point();

    SQueue queue;
    m_dt.sort_samples_from_edge(edge, queue);

    //FT start = FT(0);
    while (!queue.empty()) {
      PSample psample = queue.top();
      queue.pop();

      const FT m = psample.sample()->mass();
      const Point& ps = psample.sample()->point();

      const FT coord = psample.priority();
      const FT one_minus_coord = 1.0 - coord;

      // normal + tangential
      coef += m * one_minus_coord * one_minus_coord;
      rhs = m_traits.construct_sum_of_vectors_2_object()(
        rhs,
        m_traits.construct_scaled_vector_2_object()(
          m_traits.construct_sum_of_vectors_2_object()(
            m_traits.construct_vector_2_object()(CGAL::ORIGIN, ps),
            m_traits.construct_scaled_vector_2_object()(
              m_traits.construct_vector_2_object()(CGAL::ORIGIN, pb), -coord)),
          m * one_minus_coord));

      /*
      // only normal
      FT bin = m/M;
      FT center = start + 0.5*bin;
      Point pc = CGAL::ORIGIN + (1.0-center)*(pa - CGAL::ORIGIN) + center*(pb - CGAL::ORIGIN);
      start += bin;
      grad = grad + m*(bin*bin/12.0)*(pa - pb);
      grad = grad + m*(1.0-center)*(pc - pf);
       */

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
      if (m_dt.is_ghost(edge))
        nb_ghost++;
      else
        nb_solid++;
    }

    std::cerr << "STATS" << std::endl;
    std::cerr << "# vertices : " << m_dt.number_of_vertices()-4 << std::endl;
    std::cerr << "# isolated vertices : " << number_of_isolated_vertices() << std::endl;
    std::cerr << "# triangles: " << m_dt.number_of_faces() << std::endl;
    std::cerr << "# edges: " << m_dt.tds().number_of_edges() << std::endl;
    std::cerr << "# solid: " << nb_solid << std::endl;
    std::cerr << "# ghost: " << nb_ghost << std::endl;
  }


  /*!
    Returns the number of vertices present in the reconstructed triangulation.
   */
  std::size_t number_of_vertices() const {
    return m_dt.number_of_vertices() - 4;

  }

  /*!
    Returns the number of isolated vertices present in the reconstructed triangulation.
  */
  int number_of_isolated_vertices () const
  {
    int nb_isolated = 0;
    for (Vertex_iterator vi = m_dt.vertices_begin();
         vi != m_dt.vertices_end(); ++vi)
      {
        if (!((*vi).has_sample_assigned()))
          continue;

        typename Triangulation::Edge_circulator start = m_dt.incident_edges(vi);
        typename Triangulation::Edge_circulator cur   = start;

        do {
          if (!m_dt.is_ghost(*cur)) {
            ++nb_isolated;
            break;
          }
          ++cur;
        } while (cur != start);
      }
    return nb_isolated;
  }

  /*!
    Returns the number of (solid) edges present in the reconstructed triangulation.
   */
  int number_of_edges() const {
    int nb_solid = 0;
    for (Finite_edges_iterator ei = m_dt.finite_edges_begin();
         ei != m_dt.finite_edges_end(); ++ei)
    {
      Edge edge = *ei;
      if (m_dt.is_ghost(edge))
        continue;
      nb_solid++;
    }
    return nb_solid;
  }


  /*!
    Returns the cost of the (solid) edges present in the 
    reconstructed triangulation.
   */
  FT total_edge_cost() const {
    FT total_cost = 0;
    for (Finite_edges_iterator ei = m_dt.finite_edges_begin();
         ei != m_dt.finite_edges_end(); ++ei) 
    {
      Edge edge = *ei;
      if (m_dt.is_ghost(edge))
        continue;

      total_cost += m_dt.get_cost(edge).finalize();
    }
    return total_cost;
  }

  /// \endcond


  /// \name Simplification 
  /// You can freely mix calls of the following functions. 
  /// @{
  /*!
    Computes a shape consisting of `np` points, reconstructing the input
    points.
    \param np The number of points which will be present in the output.
    \return `true` if the number of points `np` was reached, `false`
    if the algorithm was prematurely ended because no more edge
    collapse was possible.
   */
  bool run_until(std::size_t np) {
    m_tolerance = (FT)(-1.);
    CGAL::Real_timer timer;
    if (m_verbose > 0)
      std::cerr << "reconstruct until " << np << " V";

    timer.start();
    std::size_t N = np + 4;
    std::size_t performed = 0;
    while (m_dt.number_of_vertices() > N) {
      bool ok = decimate();
      if (!ok)
        break;
      performed++;
    }

    if (m_verbose)
      std::cerr << " done" << " (" << performed
                << " iters, " << m_dt.number_of_vertices() - 4 << " V "
                << timer.time() << " s)"
                << std::endl;
    
    return (m_dt.number_of_vertices() <= N);
  }

  /*!
    Computes a shape, reconstructing the input, by performing `steps`
    edge collapse operators on the output simplex.
    \param steps The number of edge collapse operators to be performed.
    \return `true` if the required number of steps was performed,
    `false` if the algorithm was prematurely ended because no more
    edge collapse was possible.
   */
  bool run(const unsigned steps) {
    m_tolerance = (FT)(-1.);
    CGAL::Real_timer timer;
    if (m_verbose > 0)
      std::cerr << "reconstruct " << steps;

    timer.start();
    unsigned performed = 0;
    for (unsigned i = 0; i < steps; ++i) {
      bool ok = decimate();
      if (!ok)
        break;
      performed++;
    }

    if (m_verbose > 0)
      std::cerr << " done" << " (" << performed << "/"
                << steps << " iters, " << m_dt.number_of_vertices() - 4
                << " V, " << timer.time() << " s)"
                << std::endl;
    return (performed == steps);
  }


  /*!
    Computes a shape, reconstructing the input, by performing edge
    collapse operators on the output simplex until the user-defined
    tolerance is reached.

    \note The tolerance is given in the sense of the Wasserstein
    distance. It is _not_ a Hausdorff tolerance: it does not mean that
    the distance between the input samples and the output polyline is
    guaranteed to be less than `tolerance`. It means that the square
    root of transport cost per mass (homogeneous to a distance) is at
    most `tolerance`.
    
    \param tolerance Tolerance on the Wasserstein distance.
   */
  void run_under_wasserstein_tolerance (const FT tolerance) {
    m_tolerance = tolerance;
    CGAL::Real_timer timer;
    if (m_verbose > 0)
      std::cerr << "reconstruct under tolerance " << tolerance;

    timer.start();
    unsigned performed = 0;
    while (decimate ())
      performed++;

    if (m_verbose > 0)
      std::cerr << " done" << " (" << performed 
                << " iters, " << m_dt.number_of_vertices() - 4
                << " V, " << timer.time() << " s)"
                << std::endl;
  }


  /*!
    Since noise and missing data may prevent the reconstructed shape to have sharp corners well located, the algorithm offers the possibility to automatically relocate points after each edge collapse. The new location of the points is chosen such that the fitting of the output segments to the input points is improved. 
   */
  void relocate_all_points() {
    CGAL::Real_timer timer;
    if (m_verbose > 0)
      std::cerr << "relocate all points" << "...";

    timer.start();
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

      Edge_vector hull;
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

    if (m_verbose > 0)
      std::cerr << "done" << " (" << timer.time() << " s)" << std::endl;
  }

  /// @}

  /// \name Output
  /// @{

  /*!
    Writes the points and segments of the output simplex in an indexed format into output iterators.
        \tparam PointOutputIterator An output iterator with value type 
                \link Optimal_transportation_reconstruction_2::Point Point \endlink.
        \tparam IndexOutputIterator An output iterator with value type 
                `std::size_t`.
        \tparam IndexPairOutputIterator An output iterator with value type 
                `std::pair<std::size_t, std::size_t>`.

        \param points The output iterator for all points.
        \param isolated_points The output iterator for the indices of isolated points.
        \param segments The output iterator for the pairs of segment indices.
   */
  template <
    typename PointOutputIterator,
    typename IndexOutputIterator,
    typename IndexPairOutputIterator>
  CGAL::cpp11::tuple<
    PointOutputIterator, 
    IndexOutputIterator,
    IndexPairOutputIterator>
  indexed_output(
    PointOutputIterator points,
    IndexOutputIterator isolated_points,
    IndexPairOutputIterator segments) const
  {
    std::vector<Point> isolated_points_;
    std::vector<Segment> edges;

    list_output (
        std::back_inserter(isolated_points_), std::back_inserter(edges));

    // vertices_of_edges
    std::set<Point> edge_vertices;
    for (typename std::vector<Segment>::iterator it = edges.begin();
        it != edges.end(); it++) {

      Point a = (*it).source();
      Point b = (*it).target();

      edge_vertices.insert(a);
      edge_vertices.insert(b);
    }

    std::size_t count_points = 0;
    for (typename std::set<Point>::iterator it = edge_vertices.begin();
        it != edge_vertices.end(); it++) {

      *points++ = *it;
      ++count_points;
    }

    for (typename std::vector<Point>::iterator it = isolated_points_.begin();
        it != isolated_points_.end(); it++) {

      *isolated_points++ = count_points;
      *points++ = *it;
      ++count_points;
    }

    for (typename std::vector<Segment>::iterator it = edges.begin();
        it != edges.end(); it++) {

      Point const& a = it->source();
      Point const& b = it->target();

      typename std::set<Point>::iterator it_a = edge_vertices.find(a);
      typename std::set<Point>::iterator it_b = edge_vertices.find(b);

      std::size_t pos_a = std::distance(edge_vertices.begin(), it_a);
      std::size_t pos_b = std::distance(edge_vertices.begin(), it_b);

      *segments++ = std::make_pair(pos_a, pos_b);
    }

    return CGAL::cpp11::make_tuple(points, isolated_points, segments);
  }

  /*!
     Returns the solid edges and vertices present after the reconstruction
     process finished.

    \details It takes two output iterators, one for storing the
    isolated points and one for storing the edges of the reconstructed shape.

    \tparam PointOutputIterator An output iterator with value type 
            \link Optimal_transportation_reconstruction_2::Point Point \endlink.
    \tparam SegmentOutputIterator An output iterator with value type 
            \link Optimal_transportation_reconstruction_2::Segment Segment \endlink.
   */
  template<class PointOutputIterator, class SegmentOutputIterator>
  void list_output (PointOutputIterator v_it, SegmentOutputIterator e_it) const
  {
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

      Segment s = m_traits.construct_segment_2_object()(
        source->point(), target->point());
      *e_it = s;
      e_it++;
    }
  }
  /// \endcond


  /// \cond SKIP_IN_MANUAL
  const Triangulation& tds() const { return m_dt; }
  
  void extract_tds_output(Triangulation& rt2) const {
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

    // mark edges
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
} // namespace

#endif // CGAL_OPTIMAL_TRANSPORTATION_RECONSTRUCTION_2_H_
