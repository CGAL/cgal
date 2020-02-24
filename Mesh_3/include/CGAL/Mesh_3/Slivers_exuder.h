// Copyright (c) 2004-2007  INRIA Sophia-Antipolis (France).
// Copyright (c) 2008 GeometryFactory, Sophia Antipolis (France)
// Copyright (c) 2009  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Laurent Rineau, Stephane Tayeb

#ifndef CGAL_MESH_3_SLIVERS_EXUDER_H
#define CGAL_MESH_3_SLIVERS_EXUDER_H

#include <CGAL/license/Mesh_3.h>

#include <CGAL/disable_warnings.h>
#include <CGAL/Mesh_3/config.h>
#include <CGAL/Mesh_3/Concurrent_mesher_config.h>

#include <CGAL/Mesh_3/sliver_criteria.h>
#include <CGAL/Mesh_optimization_return_code.h>
#include <CGAL/Mesh_3/Null_exuder_visitor.h>
#include <CGAL/Mesh_3/Triangulation_helpers.h>

#include <CGAL/Bbox_3.h>
#include <CGAL/Double_map.h>
#include <CGAL/enum.h>
#include <CGAL/functional.h>
#include <CGAL/internal/Has_member_visited.h>
#include <CGAL/iterator.h>
#include <CGAL/Real_timer.h>

#include <CGAL/boost/iterator/transform_iterator.hpp>

#include <boost/bind.hpp>
#include <boost/format.hpp>
#include <boost/function_output_iterator.hpp>
#include <boost/optional.hpp>
#include <boost/type_traits/is_convertible.hpp>

#include <algorithm>
#include <iomanip> // std::setprecision
#include <iostream> // std::cerr/cout
#include <map>
#include <set>
#include <vector>

#ifdef CGAL_CONCURRENT_MESH_3_PROFILING
# define CGAL_PROFILE
# include <CGAL/Profile_counter.h>
#endif

#ifdef CGAL_LINKED_WITH_TBB
# include <tbb/task.h>
#endif


#ifdef CGAL_MESH_3_VERBOSE
  #define CGAL_MESH_3_EXUDER_VERBOSE
#endif


namespace CGAL {

namespace Mesh_3 {

namespace details { // various function objects

// That functor Second_of takes a pair as input (the value type of a
// map), and returns the ".second" member of that pair. It is used in
// Slivers_exuder, to constructor a transform iterator.

// It should be doable using STL bind operators, but i am not sure how
// to use them. -- Laurent Rineau, 2007/07/27
template <typename Map>
struct Second_of
  : public CGAL::cpp98::unary_function<typename Map::value_type,
                                       const typename Map::mapped_type&>
{
  typedef CGAL::cpp98::unary_function<typename Map::value_type,
                                      const typename Map::mapped_type&> Base;
  typedef typename Base::result_type                             result_type;
  typedef typename Base::argument_type                           argument_type;

  const typename Map::mapped_type&
  operator()(const typename Map::value_type& p) const
  {
    return p.second;
  }
}; // end class Second_of

} // end namespace details

/************************************************
// Class Slivers_exuder_base
// Two versions: sequential / parallel
************************************************/

// Sequential
template <typename Tr, typename Concurrency_tag>
class Slivers_exuder_base
{
protected:
  typedef typename Tr::Vertex_handle                        Vertex_handle;
  typedef typename Tr::Cell_handle                          Cell_handle;
  typedef std::vector<Cell_handle>                          Cell_vector;
  typedef typename Tr::Geom_traits                          Gt;
  typedef typename Gt::FT                                   FT;
  typedef typename std::vector<Vertex_handle>               Bad_vertices_vector;
  typedef typename Tr::Lock_data_structure                  Lock_data_structure;

  // A priority queue ordered on Tet quality (SliverCriteria)
  typedef CGAL::Double_map<Cell_handle, double>             Tet_priority_queue;
  typedef typename Tet_priority_queue::reverse_iterator     Queue_iterator;
  typedef typename Tet_priority_queue::Reverse_entry        Queue_value_type;

  Slivers_exuder_base(const Bbox_3 &, int) {}

  Lock_data_structure * get_lock_data_structure()   const { return 0; }
  void unlock_all_elements()                        const {}
  void create_root_task()                           const {}
  bool flush_work_buffers()                         const { return true; }
  void wait_for_all()                               const {}
  void destroy_root_task()                          const {}
  template <typename Func>
  void enqueue_work(Func, FT)                       const {}

protected:
  Cell_handle extract_cell_handle_from_queue_value(const Queue_value_type &qv) const
  {
    return qv.second;
  }
  FT extract_cell_quality_from_queue_value(const Queue_value_type &qv) const
  {
    return qv.first;
  }
  unsigned int extract_erase_counter_from_queue_value(const Queue_value_type &) const
  {
    return 0;
  }

  // Dummy
  unsigned int erase_counter(const Cell_handle &) const { return 0; }

  std::size_t cells_queue_size() const { return cells_queue_.size(); }
  bool cells_queue_empty()       const { return cells_queue_.empty(); }
  Queue_iterator
    cells_queue_front()                { return cells_queue_.front(); }
  void cells_queue_pop_front()         { cells_queue_.pop_front(); }
  void cells_queue_clear()             { cells_queue_.clear(); }
  void cells_queue_insert(const Cell_handle &ch, FT quality_value)
  {
    cells_queue_.insert(ch, quality_value);
  }

  /**
   * A functor to remove one \c Cell_handle from a priority queue
   */
  class Erase_from_queue
  {
  public:
    Erase_from_queue(Tet_priority_queue& queue)
    : r_queue_(queue) { }

    void operator()(const Cell_handle& cell) { r_queue_.erase(cell); }

  private:
    Tet_priority_queue& r_queue_;
  };

  /**
   * Delete cells of \c cells from \c cells_queue
   */
  void delete_cells_from_queue(const Cell_vector& cells)
  {
    std::for_each(cells.begin(), cells.end(),
                  Erase_from_queue(cells_queue_));
  }

private:

  Tet_priority_queue cells_queue_;
};


  //  The debug version of the VC++ testsuite gets a timeout when using TBB, so let's disable it
#if defined( CGAL_LINKED_WITH_TBB ) && ( !defined (BOOST_MSVC) || !defined( _DEBUG ) || !defined (CGAL_TEST_SUITE) )

// Parallel
template <typename Tr>
class Slivers_exuder_base<Tr, Parallel_tag>
{
protected:
  typedef typename Tr::Vertex_handle                        Vertex_handle;
  typedef typename Tr::Cell_handle                          Cell_handle;
  typedef std::vector<Cell_handle>                          Cell_vector;
  typedef typename Tr::Geom_traits                          Gt;
  typedef typename Gt::FT                                   FT;
  typedef typename tbb::concurrent_vector<Vertex_handle>    Bad_vertices_vector;
  typedef typename Tr::Lock_data_structure                  Lock_data_structure;

  // A priority queue ordered on Tet quality (SliverCriteria)
  typedef std::multimap<
    FT, std::pair<Cell_handle, unsigned int> >              Tet_priority_queue;
  typedef typename Tet_priority_queue::iterator             Queue_iterator;
  typedef typename Tet_priority_queue::value_type           Queue_value_type;

  Slivers_exuder_base(const Bbox_3 &bbox, int num_grid_cells_per_axis)
    : m_lock_ds(bbox, num_grid_cells_per_axis)
    , m_worksharing_ds(bbox)
  {
  }

  Lock_data_structure *get_lock_data_structure() const
  {
    return &m_lock_ds;
  }

  void unlock_all_elements() const
  {
    m_lock_ds.unlock_all_points_locked_by_this_thread();
  }

  void create_root_task() const
  {
    m_empty_root_task = new( tbb::task::allocate_root() ) tbb::empty_task;
    m_empty_root_task->set_ref_count(1);
  }

  bool flush_work_buffers() const
  {
    m_empty_root_task->set_ref_count(1);
    bool keep_flushing = m_worksharing_ds.flush_work_buffers(*m_empty_root_task);
    wait_for_all();
    return keep_flushing;
  }

  void wait_for_all() const
  {
    m_empty_root_task->wait_for_all();
  }

  void destroy_root_task() const
  {
    tbb::task::destroy(*m_empty_root_task);
    m_empty_root_task = 0;
  }

  template <typename Func>
  void enqueue_work(Func f, FT value) const
  {
    CGAL_assertion(m_empty_root_task != 0);
    m_worksharing_ds.enqueue_work(f, value, *m_empty_root_task);
  }

public:

protected:
  Cell_handle extract_cell_handle_from_queue_value(const Queue_value_type &qv) const
  {
    return qv.second.first;
  }
  FT extract_cell_quality_from_queue_value(const Queue_value_type &qv) const
  {
    return qv.first;
  }
  unsigned int extract_erase_counter_from_queue_value(const Queue_value_type &qv) const
  {
    return qv.second.second;
  }

  unsigned int erase_counter(const Cell_handle &ch) const
  {
    return ch->erase_counter();
  }

  std::size_t cells_queue_size() const { return cells_queue_.size(); }
  bool cells_queue_empty()       const { return cells_queue_.empty(); }
  Queue_iterator cells_queue_front()   { return cells_queue_.begin(); }
  void cells_queue_pop_front()         { cells_queue_.erase(cells_queue_front()); }
  void cells_queue_clear()             { cells_queue_.clear(); }
  void cells_queue_insert(const Cell_handle &ch, FT quality_value)
  {
    cells_queue_.insert(std::make_pair(
      quality_value,
      std::make_pair(ch, ch->erase_counter())));
  }

  /**
   * A functor to remove one \c Cell_handle from a priority queue
   */
  class Erase_from_queue
  {
  public:
    Erase_from_queue(Tet_priority_queue&) {}

    void operator()(const Cell_handle& cell)
    { cell->increment_erase_counter(); }
  };

  /**
   * Delete cells of \c cells from \c cells_queue
   */
  void delete_cells_from_queue(const Cell_vector& cells)
  {
    std::for_each(cells.begin(), cells.end(),
                  Erase_from_queue(cells_queue_));
  }

  mutable Lock_data_structure                 m_lock_ds;
  mutable Mesh_3::Auto_worksharing_ds         m_worksharing_ds;
  mutable tbb::task                          *m_empty_root_task;

private:
  Tet_priority_queue cells_queue_;
};
#endif // CGAL_LINKED_WITH_TBB

/************************************************
// Class Slivers_exuder
************************************************/

template <
  typename C3T3,
  typename SliverCriteria,
  typename Visitor_ = Null_exuder_visitor<C3T3>
  >
class Slivers_exuder
  : public Slivers_exuder_base<typename C3T3::Triangulation,
                               typename C3T3::Concurrency_tag>
{
public: // Types
  typedef typename C3T3::Concurrency_tag                     Concurrency_tag;
  typedef Slivers_exuder_base<
    typename C3T3::Triangulation, Concurrency_tag>           Base;

private: // Types
  typedef Slivers_exuder<C3T3, SliverCriteria, Visitor_>     Self;

  typedef typename C3T3::Triangulation                       Tr;
  typedef typename Tr::Weighted_point                        Weighted_point;
  typedef typename Tr::Bare_point                            Bare_point;
  typedef typename Tr::Cell_handle                           Cell_handle;
  typedef typename Tr::Facet                                 Facet;
  typedef typename Tr::Vertex_handle                         Vertex_handle;
  typedef typename Weighted_point::Weight                    Weight;
  typedef typename Base::Queue_value_type                    Queue_value_type;
  typedef typename Base::Cell_vector                         Cell_vector;

  typedef typename Tr::Geom_traits                           Gt;
  typedef typename Base::FT                                  FT;
  typedef typename Gt::Tetrahedron_3                         Tetrahedron_3;

  typedef typename C3T3::Cells_in_complex_iterator           Cell_iterator;
  typedef std::vector<Facet>                                 Facet_vector;

  typedef typename C3T3::Surface_patch_index                 Surface_patch_index;
  typedef typename C3T3::Subdomain_index                     Subdomain_index;
  typedef typename C3T3::Index                               Index;

  // Umbrella will store the surface_index of internal facets of a new
  // weighted point conflict zone. Such facets are represented by their edge
  // which do not contain the pumped vertex
  typedef std::pair<Vertex_handle,Vertex_handle>             Ordered_edge;
  typedef std::pair<Surface_patch_index, std::size_t>        Patch_and_counter;
  typedef std::map<Ordered_edge, Patch_and_counter>          Umbrella;

  // Boundary_facets_from_outside represents the facet of the conflict zone
  // seen from outside of it. It stores Surface_patch_index of the facet, and
  // Subdomain_index of the cell which is inside the conflict zone.
  typedef std::map<Facet,
    std::pair<Surface_patch_index, Subdomain_index> >       Boundary_facets_from_outside;

  /** Pre_star will represent the pre-star of a point. It is a (double)-map
   *  of Facet (viewed from cells inside the star), ordered by the
   *  critial_radius of the point with the cell that lies on the facet, at
   *  the exterior of the pre-star. */
  typedef CGAL::Double_map<Facet, FT>                       Pre_star;

  // Stores the value of facet for the sliver criterion functor
  typedef std::map<Facet, FT>                               Sliver_values;

  // Visitor class
  // Should define
  //  - after_cell_pumped(std::size_t cells_left_number)
  typedef Visitor_                                          Visitor;

  // Helper (for 'get_sq_distance_to_closest_vertex')
  typedef Triangulation_helpers<Tr>                         Tr_helpers;

  using Base::get_lock_data_structure;

public: // methods
  /**
   * @brief Constructor
   * @param c3t3 The mesh to exude
   * @param criteria The criterion which will be used to evaluate tet quality
   * @param d defines the maximal weight we will try:
   * max_weight(v) < d*dist(v,nearest_vertice(v))
   */
  Slivers_exuder(C3T3& c3t3,
                 const SliverCriteria& criterion,
                 const FT d = 0.45);

  /**
   * @brief pumps vertices
   * @param criterion_value_limit All vertices of tetrahedra that have a
   * quality below this bound will be pumped
   */
  Mesh_optimization_return_code operator()(Visitor visitor = Visitor())
  {
#ifdef CGAL_MESH_3_PROFILING
  WallClockTimer t;
#endif

    Mesh_optimization_return_code ret =
      pump_vertices<true>(sliver_criteria_.sliver_bound(), visitor);

#ifdef CGAL_MESH_3_PROFILING
    double exudation_time = t.elapsed();
    std::cerr << std::endl << "==== Total exudation 'wall-clock' time: "
              << exudation_time << "s ====" << std::endl;
#endif

#ifdef CGAL_MESH_3_EXPORT_PERFORMANCE_DATA
    if (ret == BOUND_REACHED)
    {
      CGAL_MESH_3_SET_PERFORMANCE_DATA("Exuder_optim_time", exudation_time);
    }
    else
    {
      CGAL_MESH_3_SET_PERFORMANCE_DATA("Exuder_optim_time",
        (ret == CANT_IMPROVE_ANYMORE ?
        "CANT_IMPROVE_ANYMORE" : "TIME_LIMIT_REACHED"));
    }
#endif

    return ret;
  }

  /// Time accessors
  void set_time_limit(double time) { time_limit_ = time; }
  double time_limit() const { return time_limit_; }

private:
  // -----------------------------------
  // Private Methods
  // -----------------------------------
  /**
   * Pumps vertices
   */
  template <bool pump_vertices_on_surfaces>
  Mesh_optimization_return_code
  pump_vertices(FT criterion_value_limit, Visitor& v);

  /**
   * Pump one vertex
   */
  template <bool pump_vertices_on_surfaces>
  bool pump_vertex(const Vertex_handle& v,
                   bool *could_lock_zone = nullptr);

  /**
   * Returns the best_weight of v
   */
  FT get_best_weight(const Vertex_handle& v,
                     bool *could_lock_zone = nullptr) const;

  /**
   * Initializes pre_star and criterion_values
   */
  void
  initialize_prestar_and_criterion_values(const Vertex_handle& v,
                                          const Cell_vector& incident_cells,
                                          Pre_star& pre_star,
                                          Sliver_values& criterion_values) const;

  /**
   * Expand pre_star with cell_to_add
   */
  bool expand_prestar(const Cell_handle& cell_to_add,
                      const Vertex_handle& pumped_vertex,
                      Pre_star& pre_star,
                      Sliver_values& criterion_values) const;

  /**
   * Returns Ordered_edge of facet which do not contains vertex
   */
  Ordered_edge get_opposite_ordered_edge(const Facet& facet,
                                         const Vertex_handle& vertex) const;

  /**
   * Returns the umbrella of internal_facets vector
   */
  boost::optional<Umbrella>
  get_umbrella(const Facet_vector& internal_facets,
               const Vertex_handle& v) const;

  /**
   * Updates the mesh with new_point
   */
  template <bool pump_vertices_on_surfaces>
  bool update_mesh(const Weighted_point& new_point,
                   const Vertex_handle& old_vertex,
                   bool *could_lock_zone = nullptr);

  /**
   * Restores cells and boundary facets of conflict zone of new_vertex in c3t3_
   */
  template <bool pump_vertices_on_surfaces>
  void restore_cells_and_boundary_facets(
    const Boundary_facets_from_outside& boundary_facets_from_outside,
    const Vertex_handle& new_vertex);

  /**
   * Restore internal facets of conflict zone of new_vertex in c3t3_
   */
  void restore_internal_facets(const Umbrella& umbrella,
                               const Vertex_handle& new_vertex);

  /**
   * Orders handles \c h1 & \c h2
   */
  template <typename Handle>
  static
  void order_two_handles(Handle& h1, Handle& h2)
  {
    if( h2 < h1 )
      std::swap(h1, h2);
  }


  template <typename Handle>
  static
  void order_three_handles(Handle& h1, Handle& h2, Handle& h3)
  {
    if(h1 > h2) std::swap(h1, h2);
    if(h2 > h3) std::swap(h2, h3);
    if(h1 > h2) std::swap(h1, h2);
  }


  /**
   * Initialization
   */
  void init(double limit_value)
  {
    if ( 0 < limit_value )
      sliver_criteria_.set_sliver_bound(limit_value);
    else
      sliver_criteria_.set_sliver_bound(sliver_criteria_.get_max_value());

    this->cells_queue_clear();
    initialize_cells_priority_queue();
    initialized_ = true;
  }


  /**
   * Initialize cells_queue w.r.t sliver_bound_
   */
  void initialize_cells_priority_queue()
  {
    for( Cell_iterator cit = c3t3_.cells_in_complex_begin() ;
        cit != c3t3_.cells_in_complex_end() ;
        ++cit)
    {
      const FT value = sliver_criteria_(cit);

      if( value < sliver_criteria_.sliver_bound() )
        this->cells_queue_insert(cit, value);
    }
  }


  /**
   * Returns the min value of second element of Ratio
   */
  FT get_min_value(const Sliver_values& criterion_values) const
  {
    using boost::make_transform_iterator;
    typedef details::Second_of<Sliver_values> Second_of;

    return *(std::min_element(
      make_transform_iterator(criterion_values.begin(), Second_of()),
      make_transform_iterator(criterion_values.end(), Second_of())));
  }


  /**
   * Returns the \c Boundary_facets_from_outside object containing mirror facets
   * of \c facets
   */
  Boundary_facets_from_outside
  get_boundary_facets_from_outside(const Facet_vector& facets) const
  {
    Boundary_facets_from_outside boundary_facets_from_outside;

    for(typename Facet_vector::const_iterator fit = facets.begin();
        fit != facets.end();
        ++fit)
    {
      boundary_facets_from_outside.insert(std::make_pair(
        tr_.mirror_facet(*fit),
        std::make_pair(c3t3_.surface_patch_index(*fit),
                       c3t3_.subdomain_index(fit->first))));
    }

    return boundary_facets_from_outside;
  }

  /**
   * Add a cell \c ch to \c cells_queue
   */
  template <bool pump_vertices_on_surfaces>
  void add_cell_to_queue(Cell_handle ch, FT criterion_value)
  {
#if defined( CGAL_LINKED_WITH_TBB ) && ( !defined (BOOST_MSVC) || !defined( _DEBUG ) || !defined (CGAL_TEST_SUITE) )
    // Parallel
    if (boost::is_convertible<Concurrency_tag, Parallel_tag>::value)
      enqueue_task<pump_vertices_on_surfaces>(
        ch, this->erase_counter(ch), criterion_value);
    // Sequential
    else
#endif
      this->cells_queue_insert(ch, criterion_value);
  }


  /**
   * A functor to remove one handle (Cell_handle/Facet_handle) from complex
   */
  class Remove_from_complex
  {
  public:
    Remove_from_complex(C3T3& c3t3)
    : c3t3_(c3t3) { }

    template <typename Handle_>
    void operator()(const Handle_& handle)
    { c3t3_.remove_from_complex(handle); }

  private:
    C3T3& c3t3_;
  };

  /**
   * Removes objects of [begin,end[ range from \c c3t3_
   */
  template<typename ForwardIterator>
  void remove_from_c3t3(ForwardIterator begin, ForwardIterator end)
  {
    std::for_each(begin, end, Remove_from_complex(c3t3_));
  }

  /**
   * Returns true if time_limit is reached
   */
  bool is_time_limit_reached() const
  {
    return ( (time_limit() > 0) && (running_time_.time() > time_limit()) );
  }

  /**
   * Returns `true` if all cells of mesh have a sliver_criteria_ value greater
   * than `sliver_bound_`
   */
  bool check_sliver_bound() const
  {
    for( Cell_iterator cit = c3t3_.cells_in_complex_begin() ;
        cit != c3t3_.cells_in_complex_end() ;
        ++cit)
    {
      const FT value = sliver_criteria_(cit);

      if( value < sliver_criteria_.sliver_bound() )
        return false;
    }

    return true;
  }


#if defined( CGAL_LINKED_WITH_TBB ) && ( !defined (BOOST_MSVC) || !defined( _DEBUG ) || !defined (CGAL_TEST_SUITE) )
  // For parallel version
  template <bool pump_vertices_on_surfaces>
  void
  enqueue_task(Cell_handle ch, unsigned int erase_counter, FT value);
#endif

private:
#if defined( CGAL_LINKED_WITH_TBB ) && ( !defined (BOOST_MSVC) || !defined( _DEBUG ) || !defined (CGAL_TEST_SUITE) )
  // Functor for enqueue_task function
  template <typename SE, bool pump_vertices_on_surfaces>
  class Pump_vertex
  {
    SE                    & m_sliver_exuder;
    const C3T3            & m_c3t3;
    Cell_handle             m_cell_handle;
    unsigned int            m_erase_counter;


  public:
    // Constructor
    Pump_vertex(SE &sliver_exuder,
                const C3T3 &c3t3,
                Cell_handle cell_handle,
                unsigned int erase_counter)
    : m_sliver_exuder(sliver_exuder),
      m_c3t3(c3t3),
      m_cell_handle(cell_handle),
      m_erase_counter(erase_counter)
    {
    }

    // Constructor
    Pump_vertex(const Pump_vertex &pvx)
    : m_sliver_exuder(pvx.m_sliver_exuder),
      m_c3t3(pvx.m_c3t3),
      m_cell_handle(pvx.m_cell_handle),
      m_erase_counter(pvx.m_erase_counter)
    {}

    // operator()
    void operator()() const
    {
#ifdef CGAL_CONCURRENT_MESH_3_PROFILING
      static Profile_branch_counter_3 bcounter(
        "early withdrawals / late withdrawals / successes [Exuder]");
#endif

      for( int i = 0; i < 4; ++i )
      {
        bool could_lock_zone;
        do
        {
          could_lock_zone = true;

          if (m_sliver_exuder.erase_counter(m_cell_handle) != m_erase_counter)
            break;

          if (!m_c3t3.triangulation().try_lock_cell(m_cell_handle))
          {
#ifdef CGAL_CONCURRENT_MESH_3_PROFILING
            bcounter.increment_branch_2(); // THIS is an early withdrawal!
#endif
            could_lock_zone = false;
            m_sliver_exuder.unlock_all_elements();
            continue;
          }

          if (m_sliver_exuder.erase_counter(m_cell_handle) != m_erase_counter)
          {
            m_sliver_exuder.unlock_all_elements();
            break;
          }

          // pump_vertices_on_surfaces is a boolean template parameter.  The
          // following condition is pruned at compiled time, if
          // pump_vertices_on_surfaces==false.
          if (pump_vertices_on_surfaces
           || m_c3t3.in_dimension(m_cell_handle->vertex(i)) > 2)
          {
            m_sliver_exuder.template pump_vertex<pump_vertices_on_surfaces>(
              m_cell_handle->vertex(i), &could_lock_zone);

#ifdef CGAL_CONCURRENT_MESH_3_PROFILING
            if (!could_lock_zone)
              bcounter.increment_branch_1(); // THIS is a late withdrawal!
            else
              ++bcounter; // Success!
#endif
          }

          m_sliver_exuder.unlock_all_elements();
        } while (!could_lock_zone);
      }

      if ( m_sliver_exuder.is_time_limit_reached() )
        tbb::task::self().cancel_group_execution();
    }
  };
#endif

  // -----------------------------------
  // Private data
  // -----------------------------------
  C3T3& c3t3_;
  Tr& tr_;
  const FT sq_delta_;

  int num_of_pumped_vertices_;
  int num_of_ignored_vertices_;
  int num_of_treated_vertices_;

  bool initialized_;
  SliverCriteria sliver_criteria_;

  // Timer
  double time_limit_;
  CGAL::Real_timer running_time_;

#ifdef CGAL_MESH_3_DEBUG_SLIVERS_EXUDER
  // -----------------------------------
  // Debug Helpers
  // -----------------------------------
private:
  /**
   * Verifies that two 'FT' are near equal
   */
  static bool near_equal(const FT d1, const FT d2)
  {
    const FT epsilon = 1e-8;
    return ( ((d1-d2) >= -1*epsilon) && ((d1-d2) <= epsilon) );
  }

  /**
   * Prints a value
   */
  static void print_FT(const FT d)
  {
    std::cerr << d << " ; ";
  }

  /** This function verifies that the pre_star contains exactly the set of
   facets given by the sequence [begin, end[.

   If v!=0, it also fills another Pre_star object, from the sequence [begin,
   end[, and checks that is in the same order as pre_star.
   */
  template <class Input_facet_it>
  bool check_pre_star(const Pre_star& pre_star,
                      Input_facet_it begin,
                      Input_facet_it end,
                      const Vertex_handle v = Vertex_handle()) const;

  /** This function verifies that the pre_star contains exactly the set of
   facets on the boundary of the conflict zone of the weighted point wp.
   The vertex handle vh is a hint for the location of wp.

   It also fills another Pre_star object, and checks that is in the same
   order as pre_star.
   */
  bool check_pre_star(const Pre_star& pre_star,
                      const Weighted_point& wp,
                      const Vertex_handle& vh) const;

  /**
   * Checks if the sliver criterion values from \c criterion_values are the same as
   * those that will be found if wp is inserted in the triangulation
   */
  bool check_ratios(const Sliver_values& criterion_values,
                    const Weighted_point& wp,
                    const Vertex_handle& vh) const;

#endif // CGAL_MESH_3_DEBUG_SLIVERS_EXUDER

}; // end class Slivers_exuder


template <typename C3T3, typename SC, typename V_>
Slivers_exuder<C3T3,SC,V_>::
Slivers_exuder(C3T3& c3t3, const SC& criteria, const FT d)
  : Base(c3t3.bbox(),
         Concurrent_mesher_config::get().locking_grid_num_cells_per_axis)
  , c3t3_(c3t3)
  , tr_(c3t3_.triangulation())
  , sq_delta_(d*d)
  , num_of_pumped_vertices_(0)
  , num_of_ignored_vertices_(0)
  , num_of_treated_vertices_(0)
  , initialized_(false)
  , sliver_criteria_(criteria)
  , time_limit_(-1)
  , running_time_()
{
  // If we're multi-thread
  tr_.set_lock_data_structure(get_lock_data_structure());
}


template <typename C3T3, typename SC, typename V_>
template <bool pump_vertices_on_surfaces>
Mesh_optimization_return_code
Slivers_exuder<C3T3,SC,V_>::
pump_vertices(FT sliver_criterion_limit,
              Visitor& visitor)
{
#ifdef CGAL_MESH_3_PROFILING
  WallClockTimer t;
#endif

  init(sliver_criterion_limit);

#ifdef CGAL_MESH_3_PROFILING
  std::cerr << std::endl << "==== Init time: "
            << t.elapsed() << "s ====" << std::endl;
#endif

#ifdef CGAL_MESH_3_EXUDER_VERBOSE
  std::cerr << "Exuding...\n";
  std::cerr << "Legend of the following line: "
            << "(#cells left,#vertices pumped,#vertices ignored)" << std::endl;

  std::cerr << "(" << this->cells_queue_size() << ",0,0)";
#endif // CGAL_MESH_3_EXUDER_VERBOSE

  running_time_.reset();
  running_time_.start();

#ifdef CGAL_MESH_3_PROFILING
  t.reset();
#endif

#if defined( CGAL_LINKED_WITH_TBB ) && ( !defined (BOOST_MSVC) || !defined( _DEBUG ) || !defined (CGAL_TEST_SUITE) )
  // Parallel
  if (boost::is_convertible<Concurrency_tag, Parallel_tag>::value)
  {
    this->create_root_task();

    while (!this->cells_queue_empty())
    {
      Queue_value_type front = *(this->cells_queue_front());
      this->cells_queue_pop_front();
      Cell_handle c = this->extract_cell_handle_from_queue_value(front);
      FT q = this->extract_cell_quality_from_queue_value(front);
      unsigned int ec = this->extract_erase_counter_from_queue_value(front);
      // Low quality first (i.e. low value of q)
      enqueue_task<pump_vertices_on_surfaces>(c, ec, q);
    }

    this->wait_for_all();

# if defined(CGAL_MESH_3_EXUDER_VERBOSE) || defined(CGAL_MESH_3_PROFILING)
    std::cerr << " Flushing";
# endif
    bool keep_flushing = true;
    while (keep_flushing)
    {
      keep_flushing = this->flush_work_buffers();
# if defined(CGAL_MESH_3_EXUDER_VERBOSE) || defined(CGAL_MESH_3_PROFILING)
      std::cerr << ".";
# endif
    }

    this->destroy_root_task();
  }
  // Sequential
  else
#endif // CGAL_LINKED_WITH_TBB
  {
    while( !this->cells_queue_empty() && !is_time_limit_reached() )
    {
      Queue_value_type front = *(this->cells_queue_front());
      Cell_handle c = this->extract_cell_handle_from_queue_value(front);

      // Low quality first (i.e. low value of cell quality)
      bool vertex_pumped = false;
      for( int i = 0; i < 4; ++i )
      {
        // pump_vertices_on_surfaces is a boolean template parameter.  The
        // following condition is pruned at compiled time, if
        // pump_vertices_on_surfaces==false.
        if( pump_vertices_on_surfaces || c3t3_.in_dimension(c->vertex(i)) > 2 )
        {
          if( pump_vertex<pump_vertices_on_surfaces>(c->vertex(i)) )
          {
            vertex_pumped = true;
            ++num_of_pumped_vertices_;
            break;
          }
          else
            ++num_of_ignored_vertices_;

          ++num_of_treated_vertices_;
        }
      }

      // if the tet could not be deleted
      if ( ! vertex_pumped )
        this->cells_queue_pop_front();

      visitor.after_cell_pumped(this->cells_queue_size());
  #ifdef CGAL_MESH_3_EXUDER_VERBOSE
      std::cerr << boost::format("\r             \r"
                                 "(%1%,%2%,%3%) (%|4$.1f| vertices/s)")
        % this->cells_queue_size()
        % num_of_pumped_vertices_
        % num_of_ignored_vertices_
        % (num_of_treated_vertices_ / running_time_.time());
  #endif // CGAL_MESH_3_EXUDER_VERBOSE
    }
  }

  running_time_.stop();

#ifdef CGAL_MESH_3_PROFILING
  std::cerr << std::endl << "==== Iterations time: "
            << t.elapsed() << "s ====" << std::endl;
#endif

#ifdef CGAL_MESH_3_EXUDER_VERBOSE
  std::cerr << std::endl;
  std::cerr << "Total exuding time: " << running_time_.time() << "s";
  std::cerr << std::endl;
#endif // CGAL_MESH_3_EXUDER_VERBOSE

  if ( is_time_limit_reached() ) {
#ifdef CGAL_MESH_3_EXUDER_VERBOSE
    std::cerr << "Exuding return code: TIME_LIMIT_REACHED\n\n";
#endif // CGAL_MESH_3_EXUDER_VERBOSE
    return TIME_LIMIT_REACHED;
  }

  if ( check_sliver_bound() ) {
#ifdef CGAL_MESH_3_EXUDER_VERBOSE
    std::cerr << "Exuding return code: BOUND_REACHED\n\n";
#endif // CGAL_MESH_3_EXUDER_VERBOSE
    return BOUND_REACHED;
  }

#ifdef CGAL_MESH_3_EXUDER_VERBOSE
    std::cerr << "Exuding return code: CANT_IMPROVE_ANYMORE\n\n";
#endif // CGAL_MESH_3_EXUDER_VERBOSE
  return CANT_IMPROVE_ANYMORE;

} // end function pump_vertices


template <typename C3T3, typename SC, typename V_>
template <bool pump_vertices_on_surfaces>
bool
Slivers_exuder<C3T3,SC,V_>::
pump_vertex(const Vertex_handle& pumped_vertex,
            bool *could_lock_zone)
{
  // Get best_weight
  FT best_weight = get_best_weight(pumped_vertex, could_lock_zone);
  if (could_lock_zone && *could_lock_zone == false)
    return false;

  typename Gt::Compare_weighted_squared_radius_3 compare_sq_radius =
    tr_.geom_traits().compare_weighted_squared_radius_3_object();

  // If best_weight <= pumped_vertex weight, nothing to do
  const Weighted_point& pumped_vertex_wp = tr_.point(pumped_vertex);
  if ( compare_sq_radius(pumped_vertex_wp, - best_weight) == CGAL::LARGER ) // best_weight > v's weight
  {
    typename Gt::Construct_point_3 cp = tr_.geom_traits().construct_point_3_object();

    const Weighted_point& pwp = tr_.point(pumped_vertex);
    Weighted_point wp(cp(pwp), best_weight);

    // Insert weighted point into mesh
    // note it can fail if the mesh is non-manifold at pumped_vertex
    return update_mesh<pump_vertices_on_surfaces>(wp,
                                                  pumped_vertex,
                                                  could_lock_zone);
  }

  return false;
}


template <typename C3T3, typename SC, typename V_>
void
Slivers_exuder<C3T3,SC,V_>::
initialize_prestar_and_criterion_values(const Vertex_handle& v,
                                        const Cell_vector& incident_cells,
                                        Pre_star& pre_star,
                                        Sliver_values& criterion_values) const
{
  for ( typename Cell_vector::const_iterator cit = incident_cells.begin() ;
       cit != incident_cells.end() ;
       ++cit )
  {
    const Cell_handle c = *cit;
    const int index = c->index(v);
    const Facet f = Facet(c, index);
    const Facet opposite_facet = tr_.mirror_facet(f);

    // Sliver criterion values initialization
    if( c3t3_.is_in_complex(c) )
    {
      criterion_values[f] = sliver_criteria_(c);
    }

    // Pre_star initialization
    // If facet is adjacent to an infinite cell, no need to put it in prestar
    // (infinite power distance radius)
    if ( tr_.is_infinite(opposite_facet.first) )
      continue;

    // Insert facet in prestar (even if it is not in complex)
    const Cell_handle opposite_cell = opposite_facet.first;
    const int index_in_opposite = opposite_cell->index(c);
    FT power_distance_to_power_sphere =
      tr_.compute_power_distance_to_power_sphere(opposite_facet.first, index_in_opposite);
    pre_star.insert(f, power_distance_to_power_sphere);
  }
}


template <typename C3T3, typename SC, typename V_>
bool
Slivers_exuder<C3T3,SC,V_>::
expand_prestar(const Cell_handle& cell_to_add,
               const Vertex_handle& pumped_vertex,
               Pre_star& pre_star,
               Sliver_values& criterion_values) const
{
  typename Gt::Compute_weight_3 cw = tr_.geom_traits().compute_weight_3_object();
  typename Gt::Construct_point_3 cp = tr_.geom_traits().construct_point_3_object();

  // Delete first facet of pre_star
  Facet start_facet = pre_star.front()->second;
  CGAL_assertion(tr_.mirror_facet(start_facet).first == cell_to_add);
#ifdef CGAL_MESH_3_DEBUG_SLIVERS_EXUDER
  FT power_distance_to_power_sphere = pre_star.front()->first;
#endif
  pre_star.pop_front();
  if ( c3t3_.is_in_complex(cell_to_add) )
  {
    criterion_values.erase(start_facet);
  }

  int start_mirror_facet_index = tr_.mirror_facet(start_facet).second;

  // For each facet of cell_to_add
  for(int i = 0; i<4 ; ++i)
  {
    // We have already treated start_facet
    if ( i == start_mirror_facet_index )
      continue;

    const Facet current_facet(cell_to_add, i);
    const Facet current_mirror_facet(tr_.mirror_facet(current_facet));
    FT new_power_distance_to_power_sphere = std::numeric_limits<FT>::infinity();

    // If current_facet_mirror is in prestar, delete it
    // (it may happen that pre_star contains two facets of the same cell)
    if ( pre_star.erase(current_mirror_facet) )
    {
      // If it is a boundary facet, stop pre_star expansion
      if ( c3t3_.is_in_complex(current_mirror_facet) )
      {
        return false;
      }

      // Update criterion_values
      if ( c3t3_.is_in_complex(cell_to_add) )
      {
        criterion_values.erase(current_mirror_facet);
      }
    }
    // If current_mirror_facet is not in prestar:
    // expand prestar & update criterion_values
    else
    {
      const Cell_handle& current_mirror_cell = current_mirror_facet.first;

      CGAL_assertion(current_mirror_cell != start_facet.first);
      CGAL_assertion(pumped_vertex != current_mirror_facet.first->vertex(0));
      CGAL_assertion(pumped_vertex != current_mirror_facet.first->vertex(1));
      CGAL_assertion(pumped_vertex != current_mirror_facet.first->vertex(2));
      CGAL_assertion(pumped_vertex != current_mirror_facet.first->vertex(3));

      // Update pre_star (we do not insert facets with infinite power distance)
      // We do insert facet of cells which are outside the complex (we just
      // don't use their sliver criterion value to get best weight)
      if ( ! tr_.is_infinite(current_mirror_cell) )
      {
        new_power_distance_to_power_sphere =
          tr_.compute_power_distance_to_power_sphere(current_mirror_cell, pumped_vertex);

        pre_star.insert(current_facet, new_power_distance_to_power_sphere);

#ifdef CGAL_MESH_3_DEBUG_SLIVERS_EXUDER
        if ( new_power_distance_to_power_sphere < power_distance_to_power_sphere )
          std::cerr << "new power distance:" << new_power_distance_to_power_sphere
                    << " / current power distance:" << power_distance_to_power_sphere
                    << std::endl;
#endif // CGAL_MESH_3_DEBUG_SLIVERS_EXUDER
      }

      // Update ratio (ratio is needed for cells of complex only)
      if ( c3t3_.is_in_complex(cell_to_add) )
      {
        const Weighted_point& pwp = tr_.point(pumped_vertex);

        // Adding the power distance to the pumped vertex's weight is a way
        // to artificially make the cell 'cell_to_add' be in conflict
        // with the pumped vertex. This is done for periodic triangulations,
        // where the function 'tetrahedron(Facet, Weighted_point)' requires
        // the cell of the facet to be in conflict with the point to determine
        // the correct offset of the weighted point (to get a correct tetrahedron).
        FT new_weight = cw(pwp);
        Facet curr_f;

        if(! tr_.is_infinite(current_mirror_cell))
        {
          // if current_mirror_cell is finite, we can re-use the value
          // 'new_power_distance_to_power_sphere'

          // Ensure that 'new_power_distance_to_power_sphere' has been initialized
          CGAL_assertion(new_power_distance_to_power_sphere != std::numeric_limits<FT>::infinity());
          new_weight += new_power_distance_to_power_sphere;

          curr_f = current_mirror_facet;
        }
        else // tr_.is_infinite(current_mirror_cell)
        {
          // We can't use 'new_power_distance_to_power_sphere' because 'current_mirror_cell'
          // is infinite. So we compute the power distance to 'cell_to_add' instead.
          CGAL_assertion(!tr_.is_infinite(cell_to_add));
          new_weight += tr_.compute_power_distance_to_power_sphere(cell_to_add, pumped_vertex);
          curr_f = current_facet;
        }

        // With 'new_weight' only, the point is only orthogonal to the power sphere,
        // so we add a little bit more weight to make sure that 'pwp' is in conflict.
        new_weight += 1e5 * std::numeric_limits<FT>::epsilon(); // 'epsilon' alone is not enough
        Weighted_point ncr_pwp(cp(pwp), new_weight);

        // The 'tetrahedron' function returns positions for which the cell is in conflict
        // with the point, which is not trivial in a periodic setting.
        Tetrahedron_3 tet = tr_.tetrahedron(curr_f, ncr_pwp);

        FT new_value = sliver_criteria_(tet);
        criterion_values.insert(std::make_pair(current_facet, new_value));
      }
    }
  }

  return true;
}


template <typename C3T3, typename SC, typename V_>
typename Slivers_exuder<C3T3,SC,V_>::FT
Slivers_exuder<C3T3,SC,V_>::
get_best_weight(const Vertex_handle& v, bool *could_lock_zone) const
{
  // incident_cells will be used both in 'initialize_prestar_and_criterion_values'
  // and to get the smallest distance between 'v' and a neighbor.
  Cell_vector incident_cells;
  incident_cells.reserve(64);
  // Parallel
  if (could_lock_zone)
  {
    if (!tr_.try_lock_and_get_incident_cells(v, incident_cells))
    {
      this->unlock_all_elements();
      *could_lock_zone = false;
      return 0.;
    }
  }
  // Sequential
  else
  {
    tr_.incident_cells(v, std::back_inserter(incident_cells));
  }

  // Get pre_star and criterion_values
  Pre_star pre_star;
  Sliver_values criterion_values;
  initialize_prestar_and_criterion_values(v, incident_cells, pre_star, criterion_values);

#ifdef CGAL_MESH_3_DEBUG_SLIVERS_EXUDER
  Pre_star pre_star_copy;
  Sliver_values ratios_copy;
#endif

  FT worst_criterion_value = get_min_value(criterion_values);
  FT best_weight = 0;

  FT sq_d_v = Tr_helpers().template get_sq_distance_to_closest_vertex
                <CGAL_NTS internal::Has_member_visited<typename Tr::Vertex> >(
                  c3t3_.triangulation(), v, incident_cells);

  // If that boolean is set to false, it means that a facet in the complex
  // is about to be flipped. In that case, the pumping is stopped.
  bool can_flip = true;

  // Main loop: find the weight which maximizes the minimum value of ratio
  while(   can_flip
        && ! pre_star.empty()
        && pre_star.front()->first < (sq_delta_ * sq_d_v)
        && ! c3t3_.is_in_complex(pre_star.front()->second) )
  {
    // Store critial radius (pre_star will be modified in expand_prestar)
    FT power_distance_to_power_sphere = pre_star.front()->first;

    // expand prestar (insert opposite_cell facets in pre_star)
    Facet link = pre_star.front()->second;
    const Cell_handle& opposite_cell = tr_.mirror_facet(link).first;

    if (could_lock_zone && !tr_.try_lock_cell(opposite_cell))
    {
      *could_lock_zone = false;
      return 0.;
    }
    can_flip = expand_prestar(opposite_cell, v, pre_star, criterion_values);

    // Update best_weight if needed
    if(can_flip)
    {
      FT min_of_pre_star = get_min_value(criterion_values);

      if( min_of_pre_star > worst_criterion_value )
      {
        // Update worst_criterion_value
        worst_criterion_value = min_of_pre_star;

        // Update best_weight
        CGAL_assertion(!pre_star.empty());
        FT next_r = pre_star.front()->first;
        best_weight = (power_distance_to_power_sphere + next_r) / 2;

#ifdef CGAL_MESH_3_DEBUG_SLIVERS_EXUDER
        pre_star_copy = pre_star;
        ratios_copy = criterion_values;
#endif // CGAL_MESH_3_DEBUG_SLIVERS_EXUDER
      }
    }
  } // end while(... can pump...)

#ifdef CGAL_MESH_3_DEBUG_SLIVERS_EXUDER
  typename Gt::Compare_weighted_squared_radius_3 compare_sq_radius =
    tr_.geom_traits().compare_weighted_squared_radius_3_object();

  const Weighted_point& vwp = tr_.point(v);
  if ( compare_sq_radius(vwp, - best_weight) == CGAL::LARGER ) // best_weight > v's weight
  {
    typename Gt::Construct_point_3 cp = tr_.geom_traits().construct_point_3_object();
    const Weighted_point& wpv = tr_.point(v);
    Weighted_point wp(cp(wpv), best_weight);
    check_pre_star(pre_star_copy, wp, v);
    check_ratios(ratios_copy, wp, v);
  }
#endif // CGAL_MESH_3_DEBUG_SLIVERS_EXUDER

  return best_weight;
}


template <typename C3T3, typename SC, typename V_>
boost::optional<typename Slivers_exuder<C3T3,SC,V_>::Umbrella >
Slivers_exuder<C3T3,SC,V_>::
get_umbrella(const Facet_vector& facets, // internal_facets of conflict zone
             const Vertex_handle& /* v, no longer used */) const
{
  Umbrella umbrella; //std::map<Ordered_edge, Patch_and_counter>

  // Insert into umbrella surface_index of facets which are on the surface
  typename Facet_vector::const_iterator fit = facets.begin();
  for ( ; fit != facets.end() ; ++fit )
  {
    if ( c3t3_.is_in_complex(*fit) )
    {
      Facet f = *fit;
      Vertex_handle v1 = f.first->vertex((f.second+1)%4);
      Vertex_handle v2 = f.first->vertex((f.second+2)%4);
      Vertex_handle v3 = f.first->vertex((f.second+3)%4);
      order_three_handles(v1, v2, v3);
      std::vector<Ordered_edge> edges;
      edges.push_back(Ordered_edge(v1, v2));
      edges.push_back(Ordered_edge(v2, v3));
      edges.push_back(Ordered_edge(v1, v3));

      for(std::size_t i = 0; i < 3; ++i)
      {
        Ordered_edge oe = edges[i];
        typename Umbrella::iterator uit = umbrella.find(oe);
        if(uit == umbrella.end()) //umbrella does not contain oe yet
        {
          umbrella.insert(std::make_pair(oe,
            std::make_pair(c3t3_.surface_patch_index(f), 1)));
        }
        else //umbrella already contains oe. Increment counter or return
        {
          std::size_t count = (*uit).second.second;
          if(count == 2) //there will be more than 3 after insertion
            return boost::none; //non-manifold configuration

          umbrella.insert(uit,
            std::make_pair(oe,
              std::make_pair(c3t3_.surface_patch_index(f), count + 1)));
        }
      }
    }
  }

  // Erase edges that have been counted twice.
  // Each Oriented_edge should appear only once.
  // Twice means that it belongs to two internal facets that are restricted.
  // Three or more corresponds to a non-manifold geometry.
  typename Umbrella::iterator uit = umbrella.begin();
  while(uit != umbrella.end())
  {
    if((*uit).second.second == 2)
    {
      typename Umbrella::iterator to_be_erased = uit++;
      umbrella.erase(to_be_erased);
    }
    else
    {
      ++uit;
    }
  }

  return umbrella;
}


template <typename C3T3, typename SC, typename V_>
template <bool pump_vertices_on_surfaces>
void
Slivers_exuder<C3T3,SC,V_>::
restore_cells_and_boundary_facets(
  const Boundary_facets_from_outside& boundary_facets_from_outside,
  const Vertex_handle& new_vertex)
{
  Cell_vector new_cells;
  new_cells.reserve(64);
  tr_.incident_cells(new_vertex, std::back_inserter(new_cells));

  // Each cell must have a facet on the boundary of the conflict zone
  CGAL_assertion(boundary_facets_from_outside.size() == new_cells.size());

  // Restore attributes of each cell
  for(typename Cell_vector::iterator cit = new_cells.begin();
      cit != new_cells.end();
      ++cit)
  {
    (*cit)->invalidate_weighted_circumcenter_cache();
    const int index = (*cit)->index(new_vertex);
    const Facet new_facet = std::make_pair(*cit, index);
    const Facet new_facet_from_outside = tr_.mirror_facet(new_facet);

    // Search new_facet_from_outside in boundary_facets_from_outside.
    // That search cannot fail.
    typename Boundary_facets_from_outside::const_iterator it =
      boundary_facets_from_outside.find(new_facet_from_outside);

    CGAL_assertion(it != boundary_facets_from_outside.end());

    // Restore facet attributes
    if ( !( it->second.first == Surface_patch_index() ) )
      c3t3_.add_to_complex(new_facet, it->second.first);

    // Restore cell attributes
    if ( !( it->second.second == Subdomain_index() ) )
      c3t3_.add_to_complex(*cit, it->second.second);

    // if the new cell is in the domain, and it criterion value is less that
    // the maximum, push it in the cells queue.
    if( c3t3_.is_in_complex(*cit) )
    {
      FT criterion_value = sliver_criteria_(*cit);

      if( criterion_value < sliver_criteria_.sliver_bound() )
        add_cell_to_queue<pump_vertices_on_surfaces>(*cit, criterion_value);
    }
  }
}


template <typename C3T3, typename SC, typename V_>
typename Slivers_exuder<C3T3,SC,V_>::Ordered_edge
Slivers_exuder<C3T3,SC,V_>::get_opposite_ordered_edge(
  const Facet& facet,
  const Vertex_handle& vertex) const
{
  CGAL_assertion(tr_.has_vertex(facet, vertex));

  Vertex_handle v1;
  Vertex_handle v2;

  // Get the two vertex of *fit which are not new_vertex
  for ( int i = 0 ; i < 4 ; ++i )
  {
    const Vertex_handle current_vertex = facet.first->vertex(i);

    if ( current_vertex != vertex && tr_.has_vertex(facet, current_vertex) )
    {
      if ( v1 == Vertex_handle() )
        v1 = current_vertex;
      else
        v2 = current_vertex;
    }
  }

  CGAL_assertion(v1 != Vertex_handle() && v2 != Vertex_handle());

  order_two_handles(v1,v2);
  return Ordered_edge(v1,v2);
}


template <typename C3T3, typename SC, typename V_>
void
Slivers_exuder<C3T3,SC,V_>::
restore_internal_facets(const Umbrella& umbrella,
                        const Vertex_handle& new_vertex)
{
  Facet_vector new_internal_facets;
  new_internal_facets.reserve(64);
  tr_.incident_facets(new_vertex, std::back_inserter(new_internal_facets));

  // Restore attributes of each facet
  for(typename Facet_vector::iterator fit = new_internal_facets.begin();
      fit != new_internal_facets.end();
      ++fit)
  {
    Ordered_edge edge = get_opposite_ordered_edge(*fit, new_vertex);

    // Search edge in umbrella.
    // If it is found, restore facet surface index from umbrella
    const typename Umbrella::const_iterator um_it = umbrella.find(edge);
    if( um_it != umbrella.end() )
    {
      c3t3_.add_to_complex(*fit, um_it->second.first);
    }
  }
}


template <typename C3T3, typename SC, typename V_>
template <bool pump_vertices_on_surfaces>
bool
Slivers_exuder<C3T3,SC,V_>::
update_mesh(const Weighted_point& new_point,
            const Vertex_handle& old_vertex,
            bool *could_lock_zone)
{
  CGAL_assertion_code(std::size_t nb_vert = tr_.number_of_vertices();)
  Cell_vector deleted_cells;
  Facet_vector internal_facets;
  Facet_vector boundary_facets;

  deleted_cells.reserve(64);
  internal_facets.reserve(64);
  boundary_facets.reserve(64);

  tr_.find_conflicts(new_point,
                     old_vertex->cell(),
                     std::back_inserter(boundary_facets),
                     std::back_inserter(deleted_cells),
                     std::back_inserter(internal_facets),
                     could_lock_zone);

  if (could_lock_zone && *could_lock_zone == false)
    return false;

  // Get some datas to restore mesh
  Boundary_facets_from_outside boundary_facets_from_outside =
    get_boundary_facets_from_outside(boundary_facets);

  boost::optional<Umbrella> umbrella = get_umbrella(internal_facets, old_vertex);
  if(umbrella == boost::none)
    return false; //abort pumping this vertex

  // Delete old cells from queue (they aren't in the triangulation anymore)
  this->delete_cells_from_queue(deleted_cells);

  // Delete old cells & facets from c3t3
  remove_from_c3t3(deleted_cells.begin(),deleted_cells.end());
  remove_from_c3t3(boundary_facets.begin(),boundary_facets.end());
  remove_from_c3t3(internal_facets.begin(),internal_facets.end());

  // Insert new point (v will be updated using a wp)
  int dimension = c3t3_.in_dimension(old_vertex);
  Index vertice_index = c3t3_.index(old_vertex);

  Vertex_handle new_vertex = tr_.insert(new_point, old_vertex->cell());
  c3t3_.set_dimension(new_vertex,dimension);
  c3t3_.set_index(new_vertex,vertice_index);

  // Only true for sequential version
  CGAL_assertion(could_lock_zone || nb_vert == tr_.number_of_vertices());

  // Restore mesh
  restore_cells_and_boundary_facets<pump_vertices_on_surfaces>(
    boundary_facets_from_outside, new_vertex);
  restore_internal_facets(*umbrella, new_vertex);

  // Only true for sequential version
  CGAL_assertion(could_lock_zone || nb_vert == tr_.number_of_vertices());

  return true; // pump was done successfully
}


#if defined( CGAL_LINKED_WITH_TBB ) && ( !defined (BOOST_MSVC) || !defined( _DEBUG ) || !defined (CGAL_TEST_SUITE) )
// For parallel version
template <typename C3T3, typename SC, typename V_>
template <bool pump_vertices_on_surfaces>
void
Slivers_exuder<C3T3,SC,V_>::
enqueue_task(Cell_handle ch, unsigned int erase_counter, FT value)
{
  this->enqueue_work(
    Pump_vertex<Self, pump_vertices_on_surfaces>(
      *this, c3t3_, ch, erase_counter),
    value);
}
#endif


#ifdef CGAL_MESH_3_DEBUG_SLIVERS_EXUDER
template <typename Pre_star>
void print_pre_stars(Pre_star pre_star_1, Pre_star pre_star_2)
{
  std::cout << "Stars: " << std::endl;
  while(!pre_star_1.empty() && !pre_star_2.empty())
  {
    std::cout << "F: " << &*((pre_star_1.front()->second).first)
                       << " " << pre_star_1.front()->second.second
              << " S: " << pre_star_1.front()->first
              << " |||||| "
              << "F: " << &*((pre_star_2.front()->second).first)
                       << " " << pre_star_2.front()->second.second
              << " S: " << pre_star_2.front()->first
              << std::endl;
    if(pre_star_1.front()->second != pre_star_2.front()->second ||
       pre_star_1.front()->first != pre_star_2.front()->first)
      std::cout << "Warning!" << std::endl;

    pre_star_1.pop_front();
    pre_star_2.pop_front();
  }
}

template <typename C3T3, typename SC, typename V_>
template <class Input_facet_it>
bool
Slivers_exuder<C3T3,SC,V_>::
check_pre_star(const Pre_star& pre_star,
               Input_facet_it begin,
               Input_facet_it end,
               const Vertex_handle v) const
{
  Pre_star pre_star_copy = pre_star;
  if(v != Vertex_handle())
  {
    Pre_star pre_star2;

    // fill pre_star2
    for(Input_facet_it fit = begin; fit != end; ++fit)
    {
      const Facet opposite_facet = tr_.mirror_facet(*fit);
      if(! tr_.is_infinite(opposite_facet.first) )
      {
        pre_star2.insert(
          *fit, tr_.compute_power_distance_to_power_sphere(opposite_facet.first, v));
      }
    }

    while(!pre_star_copy.empty() && !pre_star2.empty())
    {
      if(pre_star_copy.front()->first != pre_star2.front()->first) {
        std::cerr << "bad order\n";
        std::cerr << boost::format("pre_star.front()->first=%1%, should be %2%\n")
        % pre_star_copy.front()->first % pre_star2.front()->first;
        return false;
      }
      if ( pre_star_copy.front()->second != pre_star2.front()->second )
      {
        Facet f1 = pre_star_copy.front()->second;
        Facet f2 = pre_star2.front()->second;
        pre_star2.pop_front();
        pre_star_copy.pop_front();
        if (!pre_star_copy.empty() && !pre_star2.empty() &&
            pre_star_copy.front()->second == f2 && pre_star2.front()->second == f1 )
        {
          // It's ok
          pre_star2.pop_front();
          pre_star_copy.pop_front();
        }
        else
        {
          Facet f1 = tr_.mirror_facet(pre_star_copy.front()->second);
          Facet f2 = tr_.mirror_facet(pre_star2.front()->second);
          std::cerr << "Bad facet:" << f1.second << "/" << f2.second
          << " - " << &*f1.first << "/" << &*f2.first << std::endl;
        }
      }
      else
      {
        pre_star2.pop_front();
        pre_star_copy.pop_front();
      }
    }

    if(pre_star2.empty() && ! pre_star_copy.empty()) {
      std::cerr << "pre_star is too big!\n";
      while(!pre_star_copy.empty())
      {
        const Facet f = pre_star_copy.front()->second;
        const FT r = pre_star_copy.front()->first;
        pre_star_copy.pop_front();
        std::cerr << boost::format("extra facet (%1%,%2%) (infinite: %3%, opposite infinite: %4%), power distance: %5%\n")
        % &*f.first % f.second % tr_.is_infinite(f.first) % tr_.is_infinite(f.first->neighbor(f.second))
        % r;
      }
      return false;
    }

    if( pre_star_copy.empty() && ! pre_star2.empty() ) {
      std::cerr << "pre_star is too small!\n";
      while(!pre_star2.empty())
      {
        const Facet f = pre_star2.front()->second;
        pre_star2.pop_front();
        std::cerr << boost::format("missing facet (%1%,%2%) (infinite: %3%, opposite infinite: %4%)\n")
        % &*f.first % f.second % tr_.is_infinite(f.first) % tr_.is_infinite(f.first->neighbor(f.second));
      }
      return false;
    }
  }

  pre_star_copy = pre_star;

  for(Input_facet_it fit = begin;
      fit != end;
      ++fit)
  {
    const Facet opposite_facet = tr_.mirror_facet(*fit);
    if(!tr_.is_infinite(opposite_facet.first) && !pre_star_copy.erase(*fit))
      return false;
  }
  if( !pre_star_copy.empty() )
    return false;

  return true;
}



template <typename C3T3, typename SC, typename V_>
bool
Slivers_exuder<C3T3,SC,V_>::
check_pre_star(const Pre_star& pre_star,
               const Weighted_point& wp,
               const Vertex_handle& vh) const
{
  std::vector<Facet> boundary_facets;
  boundary_facets.reserve(64);

  tr_.find_conflicts(wp,
                     vh->cell(),
                     std::back_inserter(boundary_facets),
                     CGAL::Emptyset_iterator(),
                     CGAL::Emptyset_iterator());

  const bool result =  check_pre_star(pre_star,
                                      boundary_facets.begin(),
                                      boundary_facets.end(),
                                      vh);
  if( ! result )
  {
    std::cerr << "tested wp=" << wp << '\n';
    std::cerr << "boundary_facets.size()=" << boundary_facets.size() << std::endl;
    typename std::vector<Facet>::const_iterator bfit = boundary_facets.begin(),
                                                bfend = boundary_facets.end();
    for(; bfit != bfend; ++bfit)
      std::cerr << "Cell: " << &*(bfit->first) << " second: " << bfit->second << '\n';

    std::cerr << "\npre_star.size()=" << pre_star.size() << std::endl;
    typename Pre_star::const_iterator psit = pre_star.begin(), psend = pre_star.end();
    for(; psit != psend; ++psit)
      std::cerr << "Cell: " << &*(psit->first.first) << " second: " << psit->first.second << '\n';
    std::cerr << std::endl;
  }

  return result;
}


template <typename C3T3, typename SC, typename V_>
bool
Slivers_exuder<C3T3,SC,V_>::
check_ratios(const Sliver_values& criterion_values,
                  const Weighted_point& wp,
                  const Vertex_handle& vh) const
{
  Cell_vector deleted_cells;
  Facet_vector internal_facets;
  Facet_vector boundary_facets;

  tr_.find_conflicts(wp,
                     vh->cell(),
                     std::back_inserter(boundary_facets),
                     std::back_inserter(deleted_cells),
                     std::back_inserter(internal_facets));

  bool result = true;
  std::vector<FT> expected_ratios;
  std::vector<FT> ratio_vector;

  for ( typename Sliver_values::const_iterator rit = criterion_values.begin() ;
       rit != criterion_values.end() ;
       ++rit )
  {
    ratio_vector.push_back(rit->second);
  }

  for ( typename Facet_vector::const_iterator it = boundary_facets.begin() ;
       it != boundary_facets.end() ;
       ++ it )
  {
    if ( !c3t3_.is_in_complex((it->first)) )
      continue;

    Tetrahedron_3 tet = tr_.tetrahedron(*it, wp);
    FT ratio = sliver_criteria_(tet);
    expected_ratios.push_back(ratio);

    bool found = false;
    for ( typename Sliver_values::const_iterator rit = criterion_values.begin() ;
         rit != criterion_values.end() ;
         ++rit )
    {
      if ( near_equal(rit->second,ratio) )
      {
        found = true;
        break;
      }
    }
    if ( ! found )
    {
      result = false;
    }
  }

  if (expected_ratios.size() != criterion_values.size())
    result = false;

  if ( !result )
  {
    std::sort(expected_ratios.begin(),expected_ratios.end());
    std::sort(ratio_vector.begin(),ratio_vector.end());
    std::vector<FT> diff;
    std::set_difference(expected_ratios.begin(),expected_ratios.end(),
                        ratio_vector.begin(),ratio_vector.end(),
                        std::back_inserter(diff));

    std::cerr << "\nExpected criterion_values (size: "
              << expected_ratios.size() << ") [";
    std::for_each(expected_ratios.begin(), expected_ratios.end(), print_FT);
    std::cerr << "]\nRatios (size: " << ratio_vector.size() << ") [";
    std::for_each(ratio_vector.begin(), ratio_vector.end(), print_FT);
    std::cerr << "]\nDiff:[";
    std::for_each(diff.begin(),diff.end(), print_FT);
    std::cerr << "]\n";
  }

  return result;
}
#endif // CGAL_MESH_3_DEBUG_SLIVERS_EXUDER

} // end namespace Mesh_3

} // end namespace CGAL

#include <CGAL/enable_warnings.h>
#endif // end CGAL_MESH_3_SLIVERS_EXUDER_H
