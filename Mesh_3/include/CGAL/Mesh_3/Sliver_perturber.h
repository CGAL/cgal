// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Stephane Tayeb
//
//******************************************************************************
// File Description :
//******************************************************************************

#ifndef CGAL_MESH_3_SLIVER_PERTURBER_H
#define CGAL_MESH_3_SLIVER_PERTURBER_H

#include <CGAL/license/Mesh_3.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/Mesh_3/config.h>


#ifdef CGAL_MESH_3_VERBOSE
  #define CGAL_MESH_3_PERTURBER_VERBOSE
#endif

#ifdef CGAL_MESH_3_PERTURBER_VERBOSE
  #ifndef CGAL_MESH_3_PERTURBER_HIGH_VERBOSITY
    #define CGAL_MESH_3_PERTURBER_LOW_VERBOSITY
  #else
    #undef CGAL_MESH_3_PERTURBER_LOW_VERBOSITY
  #endif
#endif

#include <CGAL/Mesh_3/vertex_perturbation.h>
#include <CGAL/Mesh_3/C3T3_helpers.h>
#include <CGAL/Mesh_optimization_return_code.h>
#include <CGAL/Real_timer.h>
#include <CGAL/Mesh_3/Null_perturber_visitor.h>
#include <CGAL/Mesh_3/sliver_criteria.h>
#include <CGAL/Time_stamper.h>

#include <CGAL/Mesh_3/Concurrent_mesher_config.h>
#include <CGAL/Mesh_3/Worksharing_data_structures.h>

#ifdef CGAL_CONCURRENT_MESH_3_PROFILING
# define CGAL_PROFILE
# include <CGAL/Profile_counter.h>
#endif

#ifdef CGAL_LINKED_WITH_TBB
# include <tbb/task_group.h>
#endif

#include <boost/format.hpp>
#ifdef CGAL_MESH_3_USE_RELAXED_HEAP
#  error This option CGAL_MESH_3_USE_RELAXED_HEAP is no longer supported
// The reason is that the Boost relaxed heap does not ensure a strict order
// of the priority queue.
#include <boost/pending/relaxed_heap.hpp>
#else
#include <CGAL/Modifiable_priority_queue.h>
#endif //CGAL_MESH_3_USE_RELAXED_HEAP
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/type_traits/is_convertible.hpp>

#include <boost/unordered_map.hpp>

#include <memory>

namespace CGAL {

namespace Mesh_3 {

/**
* @class PVertex
* Vertex with associated perturbation datas
*/
// Sequential
template< typename FT
        , typename Vertex_handle
        , typename Point_3
        , typename SliverCriterion
        , typename Perturbation
        , typename Concurrency_tag>
class PVertex_
{
public:
typedef PVertex_<FT,
                 Vertex_handle,
                 Point_3,
                 SliverCriterion,
                 Perturbation,
                 Concurrency_tag> Self;
typedef std::size_t id_type;

/// Constructor
PVertex_()
: vertex_handle_()
, incident_sliver_nb_(0)
, min_value_((std::numeric_limits<double>::max)())
, try_nb_(0)
, p_perturbation_(nullptr)
, id_()
{ }

PVertex_(const Vertex_handle& vh, id_type id)
: vertex_handle_(vh)
, incident_sliver_nb_(0)
, min_value_((std::numeric_limits<double>::max)())
, try_nb_(0)
, p_perturbation_(nullptr)
, id_(id)
{ }

/// Associated vertex
const Vertex_handle& vertex() const { return vertex_handle_; }
void set_vertex(const Vertex_handle& vh) { vertex_handle_ = vh; }

/// Incident slivers number
unsigned int sliver_nb() const { return incident_sliver_nb_; }
void set_sliver_nb(const unsigned int n) { incident_sliver_nb_ = n; }

/// Current perturbation
const Perturbation* perturbation() const { return p_perturbation_; }
void set_perturbation(const Perturbation* p) { p_perturbation_ = p; }

/// Is perturbable
bool is_perturbable() const
{
  return (   (vertex_handle_->in_dimension() > 1)
          && (nullptr != perturbation())
          && (sliver_nb() != 0) );
}

/// Min sliver value
const FT& min_value() const { return min_value_; }
void set_min_value(const FT& min_value){ min_value_ = min_value; }

/// Try nb
const unsigned int& try_nb() const { return try_nb_; }
void set_try_nb(const unsigned int& try_nb) { try_nb_ = try_nb; }
void increment_try_nb() { ++try_nb_; }

/// Id
void set_id(const id_type& id) { id_ = id; }
id_type id() const { return id_; }

/// Operators
bool operator==(const Self& pv) const { return ( id() == pv.id() ); }

bool operator<(const Self& pv) const
{
  // vertex type (smallest-interior first)
  if ( vertex()->in_dimension() != pv.vertex()->in_dimension() )
    return vertex()->in_dimension() > pv.vertex()->in_dimension();
  // nb incident slivers (smallest first)
  else if ( sliver_nb() != pv.sliver_nb() )
    return sliver_nb() < pv.sliver_nb();
  // min angle (smallest first)
  else if ( min_value() != pv.min_value() )
    return min_value() < pv.min_value();
  // try nb (smallest first)
  else if ( try_nb() != pv.try_nb() )
    return try_nb() < pv.try_nb();
  // perturbation type (smallest first)
  else if ( perturbation() != pv.perturbation() )
    return *perturbation() < *pv.perturbation();
  return ( id() < pv.id() ); // all characteristics are the same!
}

/// Dummy functions
void update_saved_erase_counter() {}
bool is_zombie() { return false; }

private:
/// Private datas
Vertex_handle vertex_handle_;
unsigned int incident_sliver_nb_;
FT min_value_;
unsigned int try_nb_;
const Perturbation* p_perturbation_;
id_type id_;
};

#ifdef CGAL_LINKED_WITH_TBB
// Parallel
template< typename FT
        , typename Vertex_handle
        , typename Point_3
        , typename SliverCriterion
        , typename Perturbation>
class PVertex_<FT, Vertex_handle, Point_3,
               SliverCriterion, Perturbation, Parallel_tag>
{
public:
typedef PVertex_<FT,
                 Vertex_handle,
                 Point_3,
                 SliverCriterion,
                 Perturbation,
                 Parallel_tag> Self;
typedef std::size_t id_type;

/// Constructor
PVertex_()
: vertex_handle_()
, in_dimension_(-1)
, incident_sliver_nb_(0)
, min_value_((std::numeric_limits<double>::max)())
, try_nb_(0)
, p_perturbation_(nullptr)
, id_()
{ }

PVertex_(const Vertex_handle& vh, id_type id)
: vertex_handle_(vh)
, vh_erase_counter_when_added_(vh->erase_counter())
, in_dimension_(vh->in_dimension())
, incident_sliver_nb_(0)
, min_value_((std::numeric_limits<double>::max)())
, try_nb_(0)
, p_perturbation_(nullptr)
, id_(id)
{ }

/// Associated vertex
const Vertex_handle& vertex() const { return vertex_handle_; }
void set_vertex(const Vertex_handle& vh)
{
  vertex_handle_ = vh;
  update_saved_erase_counter();
}
void update_saved_erase_counter()
{
  vh_erase_counter_when_added_ = vertex_handle_->erase_counter();
}

int in_dimension() const { return in_dimension_; }

/// Incident slivers number
unsigned int sliver_nb() const { return incident_sliver_nb_; }
void set_sliver_nb(const unsigned int n) { incident_sliver_nb_ = n; }

/// Current perturbation
const Perturbation* perturbation() const { return p_perturbation_; }
void set_perturbation(const Perturbation* p) { p_perturbation_ = p; }

/// Is perturbable
bool is_perturbable() const
{
  return (   (vertex_handle_->in_dimension() > 1)
          && (nullptr != perturbation())
          && (sliver_nb() != 0) );
}

/// Min sliver value
const FT& min_value() const { return min_value_; }
void set_min_value(const FT& min_value){ min_value_ = min_value; }

/// Try nb
const unsigned int& try_nb() const { return try_nb_; }
void set_try_nb(const unsigned int& try_nb) { try_nb_ = try_nb; }
void increment_try_nb() { ++try_nb_; }

/// Id
void set_id(const id_type& id) { id_ = id; }
id_type id() const { return id_; }

/// Zombie
bool is_zombie() const
{
  return vertex_handle_->erase_counter() != vh_erase_counter_when_added_;
}

/// Operators
bool operator==(const Self& pv) const { return ( id() == pv.id() ); }

bool operator<(const Self& pv) const
{
  // vertex type (smallest-interior first)
  if ( in_dimension() != pv.in_dimension() )
    return in_dimension() > pv.in_dimension();
  // nb incident slivers (smallest first)
  else if ( sliver_nb() != pv.sliver_nb() )
    return sliver_nb() < pv.sliver_nb();
  // min angle (smallest first)
  else if ( min_value() != pv.min_value() )
    return min_value() < pv.min_value();
  // try nb (smallest first)
  else if ( try_nb() != pv.try_nb() )
    return try_nb() < pv.try_nb();
  // perturbation type (smallest first)
  else if ( perturbation() != pv.perturbation() )
    return *perturbation() < *pv.perturbation();
  return ( id() < pv.id() ); // all characteristics are the same!
}

private:
/// Private datas
Vertex_handle vertex_handle_;
unsigned int vh_erase_counter_when_added_;
int in_dimension_;
unsigned int incident_sliver_nb_;
FT min_value_;
unsigned int try_nb_;
const Perturbation* p_perturbation_;
id_type id_;
};
#endif

/************************************************
// Class Sliver_perturber_base
// Two versions: sequential / parallel
************************************************/

// Sequential
template <typename Tr, typename Concurrency_tag>
class Sliver_perturber_base
{
protected:
  typedef typename Tr::Vertex_handle                        Vertex_handle;
  typedef typename Tr::Geom_traits                          Gt;
  typedef typename Gt::FT                                   FT;
  typedef typename std::vector<Vertex_handle>               Bad_vertices_vector;
  typedef typename Tr::Lock_data_structure                  Lock_data_structure;


  Sliver_perturber_base(const Bbox_3 &, int) {}

  Lock_data_structure *
    get_lock_data_structure()                       const { return 0; }
  void unlock_all_elements()                        const {}
  void create_task_group()                          const {}
  bool flush_work_buffers()                         const { return true; }
  void wait_for_all()                               const {}
  void destroy_trask_group()                        const {}
  template <typename Func, typename PVertex>
  void enqueue_work(Func, const PVertex &)          const {}

  void increment_erase_counter(const Vertex_handle &) const {}
};

#ifdef CGAL_LINKED_WITH_TBB
// Parallel
template <typename Tr>
class Sliver_perturber_base<Tr, Parallel_tag>
{
protected:
  typedef typename Tr::Vertex_handle                        Vertex_handle;
  typedef typename Tr::Geom_traits                          Gt;
  typedef typename Gt::FT                                   FT;
  typedef typename tbb::concurrent_vector<Vertex_handle>    Bad_vertices_vector;
  typedef typename Tr::Lock_data_structure                  Lock_data_structure;

  Sliver_perturber_base(const Bbox_3 &bbox, int num_grid_cells_per_axis)
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

  void create_task_group() const
  {
    m_task_group = new tbb::task_group;
  }

  bool flush_work_buffers() const
  {
    bool keep_flushing = m_worksharing_ds.flush_work_buffers(*m_task_group);
    wait_for_all();
    return keep_flushing;
  }

  void wait_for_all() const
  {
    m_task_group->wait();
  }

  void destroy_trask_group() const
  {
    delete m_task_group;
    m_task_group = 0;
  }

  template <typename Func, typename PVertex>
  void enqueue_work(Func f, const PVertex &pv) const
  {
    CGAL_assertion(m_task_group != 0);
    m_worksharing_ds.enqueue_work(f, pv, *m_task_group);
  }

  void increment_erase_counter(const Vertex_handle &vh) const
  {
    vh->increment_erase_counter();
  }

  void cancel() const {
    return m_task_group->cancel();
  }

public:

protected:
  mutable Lock_data_structure           m_lock_ds;
  mutable Mesh_3::Auto_worksharing_ds   m_worksharing_ds;
  mutable tbb::task_group               *m_task_group;
};
#endif // CGAL_LINKED_WITH_TBB



/************************************************
// Class Sliver_perturber
************************************************/

template < typename C3T3,
           typename MeshDomain,
           typename SliverCriterion = Mesh_3::Min_dihedral_angle_criterion
                                  <typename C3T3::Triangulation::Geom_traits>,
           typename Visitor_ = Null_perturber_visitor<C3T3> >
class Sliver_perturber
: public Sliver_perturber_base<typename C3T3::Triangulation,
                               typename C3T3::Concurrency_tag>
{
  // Types
  typedef typename C3T3::Concurrency_tag Concurrency_tag;

  typedef Sliver_perturber<C3T3, MeshDomain, SliverCriterion, Visitor_> Self;
  typedef Sliver_perturber_base<
    typename C3T3::Triangulation, Concurrency_tag>                      Base;

  typedef typename C3T3::Triangulation          Tr;
  typedef typename Tr::Geom_traits              Gt;

  typedef typename Tr::Cell_handle              Cell_handle;
  typedef typename Base::Vertex_handle          Vertex_handle;
  typedef typename Tr::Vertex                   Vertex;

  typedef typename Tr::Bare_point               Bare_point;
  typedef typename Tr::Weighted_point           Weighted_point;

  typedef typename std::vector<Cell_handle>     Cell_vector;
  typedef typename std::vector<Vertex_handle>   Vertex_vector;
  typedef typename Base::Bad_vertices_vector    Bad_vertices_vector;

  typedef typename Gt::FT                       FT;

  // Helper
  typedef class C3T3_helpers<C3T3,MeshDomain> C3T3_helpers;

  using Base::get_lock_data_structure;

  // Visitor
  // Should define
  //  - bound_reached(FT bound)
  //  - end_of_perturbation_iteration(std::size_t vertices_left)
  typedef Visitor_ Visitor;

  // perturbations
public:
  typedef Abstract_perturbation<C3T3,MeshDomain,SliverCriterion> Perturbation;
  typedef boost::ptr_vector<Perturbation>                 Perturbation_vector;

private:
  // Relaxed heap

  typedef PVertex_<FT,
                   Vertex_handle,
                   Bare_point,
                   SliverCriterion,
                   Perturbation,
                   Concurrency_tag> PVertex;

  /**
   * @class PVertex_id
   * relaxed heap
   */
  class PVertex_id :
  public boost::put_get_helper<typename PVertex::id_type, PVertex_id>
  {
  public:
    typedef boost::readable_property_map_tag category;
    typedef typename PVertex::id_type value_type;
    typedef typename PVertex::id_type reference;
    typedef PVertex key_type;

    value_type operator[] (const key_type& pv) const { return pv.id(); }
  };

  typedef std::less<PVertex> less_PVertex;
  #ifdef CGAL_MESH_3_USE_RELAXED_HEAP
  typedef boost::relaxed_heap<PVertex, less_PVertex, PVertex_id> PQueue;
  #else
  typedef ::CGAL::internal::mutable_queue_with_remove<PVertex,std::vector<PVertex>, less_PVertex, PVertex_id> PQueue;
  #endif //CGAL_MESH_3_USE_RELAXED_HEAP

public:
  /**
   * Constructor
   */
  Sliver_perturber(C3T3& c3t3,
                   const MeshDomain& domain,
                   const SliverCriterion& criterion);

  /**
   * @brief Launch perturbation
   * @param sliver_bound the bound the perturber will try to achieve
   * @param delta the size of the step used by the perturber
   *
   * Runs explicit perturbation. The goal is that for each tet of the mesh,
   * SliverCriterion(tet) > sliver_bound.
   * The perturber runs step by step, using delta as step size.
   */
  Mesh_optimization_return_code
  operator()(Visitor visitor = Visitor());

  /**
   * Adds a perturbation at the end of the perturbation queue
   */
  void add_perturbation(Perturbation* perturbation);

  /// Time accessors
  void set_time_limit(double time) { time_limit_ = time; }
  double time_limit() const { return time_limit_; }

private:

  // -----------------------------------
  // Private methods
  // -----------------------------------

  /**
   * One step perturbation: tries to achieve sliver_bound quality in the mesh
   */
  bool perturb(const FT& sliver_bound, PQueue& pqueue, Visitor& v) const;

  /**
   * Builds priority queue. It will contain all vertices that have quality below
   * sliver_bound.
   * Returns priority queue size.
   *
   * precondition: pqueue.empty()
   */
  int build_priority_queue(const FT& sliver_bound, PQueue& pqueue) const;

  /**
   * Updates priority queue for all vertices of \c vertices
   */
  // Sequential
  int update_priority_queue(const Vertex_vector& vertices,
                            const FT& sliver_bound,
                            PQueue& pqueue) const;
#ifdef CGAL_LINKED_WITH_TBB
  // Parallel
  int update_priority_queue( const Vertex_vector& vertices
                           , const FT& sliver_bound
                           , Visitor& visitor
                           , Bad_vertices_vector &bad_vertices) const;
#endif

  /**
   * Updates \c pv in priority queue
   */
  int update_priority_queue(const PVertex& pv, PQueue& pqueue) const;

  // For parallel version
  void
  perturb_vertex( PVertex pv
                , const FT& sliver_bound
                , Visitor& visitor
                , Bad_vertices_vector &bad_vertices
                , bool *could_lock_zone
                ) const;

  /**
   * Returns a pvertex from a vertex handle \c vh, using id \c pv_id
   */
  PVertex
  make_pvertex(const Vertex_handle& vh,
               const FT& sliver_bound,
               const typename PVertex::id_type& pv_id) const;
  PVertex
  make_pvertex__concurrent(
               const Vertex_handle& vh,
               const FT& sliver_bound,
               const typename PVertex::id_type& pv_id) const;

  /**
   * Updates a pvertex \c pv
   */
  void update_pvertex(PVertex& pv, const FT& sliver_bound) const;
  void update_pvertex__concurrent(PVertex& pv, const FT& sliver_bound) const;

  /**
   * Returns \c vh pvertex id
   */
  typename PVertex::id_type get_pvertex_id(const Vertex_handle& vh) const
  {
    return static_cast<typename PVertex::id_type>(vh->meshing_info());
  }

  /**
   * Update bad vertices vector, wrt \c sliver_bound
   */
  // Sequential
  void update_bad_vertices(std::vector<Vertex_handle> &bad_vertices,
                           const FT& sliver_bound) const;
#ifdef CGAL_LINKED_WITH_TBB
  // Parallel
  void update_bad_vertices(tbb::concurrent_vector<Vertex_handle> &bad_vertices,
                           const FT& sliver_bound) const;
#endif

  /**
   * Initializes vertices ids
   */
  void initialize_vertices_id() const;

  /**
   * Returns true if time_limit is reached
   */
  bool is_time_limit_reached() const
  {
    return ( (time_limit() > 0) && (running_time_.time() > time_limit()) );
  }


#ifdef CGAL_MESH_3_PERTURBER_VERBOSE
  /// Verbose mode methods
  void print_perturbations_statistics() const;
  void print_final_perturbations_statistics() const;
  void reset_perturbation_counters();
#endif

#ifdef CGAL_LINKED_WITH_TBB
  // For parallel version
  void
  enqueue_task(const PVertex &pv,
               const FT& sliver_bound,
               Visitor& visitor,
               Bad_vertices_vector &bad_vertices
               ) const;
#endif

private:

#ifdef CGAL_LINKED_WITH_TBB

  // Functor for enqueue_task function
  template <typename SP, typename Bad_vertices_vector_>
  class Perturb_vertex
  {
    const SP              & m_sliver_perturber;
    PVertex                 m_pv;
    FT                      m_sliver_bound;
    Visitor               & m_visitor;
    Bad_vertices_vector_  & m_bad_vertices;

  public:
    // Constructor
    Perturb_vertex(const SP &sp,
                   const PVertex &pv,
                   FT sliver_bound,
                   Visitor& visitor,
                   Bad_vertices_vector_ &bad_vertices)
    : m_sliver_perturber(sp),
      m_pv(pv),
      m_sliver_bound(sliver_bound),
      m_visitor(visitor),
      m_bad_vertices(bad_vertices)
    {
    }

    // Constructor
    Perturb_vertex(const Perturb_vertex &pvx)
    : m_sliver_perturber(pvx.m_sliver_perturber),
      m_pv(pvx.m_pv),
      m_sliver_bound(pvx.m_sliver_bound),
      m_visitor(pvx.m_visitor),
      m_bad_vertices(pvx.m_bad_vertices)
    {}

    // operator()
    void operator()() const
    {
      bool could_lock_zone;
      do
      {
        m_sliver_perturber.perturb_vertex(
          m_pv, m_sliver_bound, m_visitor, m_bad_vertices, &could_lock_zone);
        m_sliver_perturber.unlock_all_elements();
      } while (!could_lock_zone);

      if ( m_sliver_perturber.is_time_limit_reached() )
        m_sliver_perturber.cancel();
    }
  };
#endif

  // -----------------------------------
  // Private data
  // -----------------------------------
  C3T3& c3t3_;
  Tr& tr_;
  const MeshDomain& domain_;
  SliverCriterion sliver_criterion_;
  Perturbation_vector perturbation_vector_;
  C3T3_helpers helper_;

  // Internal perturbation ordering
  int next_perturbation_order_;

  // Timer
  double time_limit_;
  CGAL::Real_timer running_time_;
};







template <typename C3T3, typename Md, typename Sc, typename V_>
Sliver_perturber<C3T3,Md,Sc,V_>::
Sliver_perturber(C3T3& c3t3,
                 const Md& domain,
                 const Sc& criterion)
  : Base(c3t3.bbox(),
         Concurrent_mesher_config::get().locking_grid_num_cells_per_axis)
  , c3t3_(c3t3)
  , tr_(c3t3_.triangulation())
  , domain_(domain)
  , sliver_criterion_(criterion)
  , helper_(c3t3_,domain_,get_lock_data_structure())
  , next_perturbation_order_(0)
  , time_limit_(-1)
  , running_time_()
{
  // If we're multi-thread
  tr_.set_lock_data_structure(get_lock_data_structure());
}



template <typename C3T3, typename Md, typename Sc, typename V_>
Mesh_optimization_return_code
Sliver_perturber<C3T3,Md,Sc,V_>::
operator()(Visitor visitor)
{
  //check criterion bound
  if ( sliver_criterion_.sliver_bound() == 0 )
    sliver_criterion_.set_sliver_bound(sliver_criterion_.get_default_value());

  // Reset sliver value cache
  helper_.reset_cache();

  // Init time counter
  if (running_time_.is_running())
    running_time_.stop();
  running_time_.reset();
  running_time_.start();

#ifdef CGAL_MESH_3_PROFILING
  WallClockTimer t;
#endif

  // Build priority queue (we use one queue for all steps)
  PQueue pqueue(tr_.number_of_vertices());

  // Initialize vertices ids
  initialize_vertices_id();


#if defined(CGAL_MESH_3_PERTURBER_VERBOSE) \
 || defined(CGAL_MESH_3_PROFILING)
  std::cerr << "Running sliver perturbation..." << std::endl;
#endif

#ifdef CGAL_MESH_3_PERTURBER_LOW_VERBOSITY
  std::cerr << "Legend of the following line: "
            << "(#vertices in pqueue, #iterations, #fails)" << std::endl;
#endif

  const FT& delta = sliver_criterion_.get_perturbation_unit();
  FT current_bound = delta;
  bool perturbation_ok = true;
  while(current_bound <= sliver_criterion_.sliver_bound() && perturbation_ok)
  {
#ifdef CGAL_MESH_3_PERTURBER_HIGH_VERBOSITY
    // reset_perturbation_counters is not const
    reset_perturbation_counters();
#endif
    perturbation_ok = perturb(current_bound, pqueue, visitor);

    visitor.bound_reached(current_bound);

    current_bound += delta;
    if ( (current_bound >= sliver_criterion_.sliver_bound())
        && (current_bound < sliver_criterion_.sliver_bound() + delta) )
    {
      current_bound = sliver_criterion_.sliver_bound();
    }
  }

#ifdef CGAL_MESH_3_PROFILING
  double perturbation_time = t.elapsed();
#endif

  running_time_.stop();
  helper_.reset_cache();//in case we re-use caches in another operation
                               // after this perturbation

#ifdef CGAL_MESH_3_PERTURBER_VERBOSE
  std::cerr << std::endl
            << "Total perturbation time: " << running_time_.time() << "s";
  std::cerr << std::endl << "Perturbation statistics:" << std::endl;
  print_final_perturbations_statistics();
#endif

#ifdef CGAL_MESH_3_PROFILING
  std::cerr << std::endl << "Total perturbation 'wall-clock' time: "
            << perturbation_time << "s" << std::endl;
#endif

  Mesh_optimization_return_code ret;

  if ( is_time_limit_reached() ) {
#if defined(CGAL_MESH_3_PERTURBER_VERBOSE) || defined(CGAL_MESH_3_PROFILING)
    std::cerr << "Perturbation return code: TIME_LIMIT_REACHED\n\n";
#endif // CGAL_MESH_3_PERTURBER_VERBOSE
    ret = TIME_LIMIT_REACHED;
  }
  else if ( !perturbation_ok ) {
#if defined(CGAL_MESH_3_PERTURBER_VERBOSE) || defined(CGAL_MESH_3_PROFILING)
    std::cerr << "Perturbation return code: CANT_IMPROVE_ANYMORE\n\n";
#endif // CGAL_MESH_3_PERTURBER_VERBOSE
    ret = CANT_IMPROVE_ANYMORE;
  }
  else
  {
#if defined(CGAL_MESH_3_PERTURBER_VERBOSE) || defined(CGAL_MESH_3_PROFILING)
    std::cerr << "Perturbation return code: BOUND_REACHED\n\n";
#endif // CGAL_MESH_3_PERTURBER_VERBOSE
    ret = BOUND_REACHED;
  }

#if defined(CGAL_MESH_3_EXPORT_PERFORMANCE_DATA) \
 && defined(CGAL_MESH_3_PROFILING)
  if (ret == BOUND_REACHED)
  {
    CGAL_MESH_3_SET_PERFORMANCE_DATA("Perturber_optim_time", perturbation_time);
  }
  else
  {
    CGAL_MESH_3_SET_PERFORMANCE_DATA("Perturber_optim_time",
      (ret == CANT_IMPROVE_ANYMORE ?
      "CANT_IMPROVE_ANYMORE" : "TIME_LIMIT_REACHED"));
  }
#endif

  return ret;
}

template <typename C3T3, typename Md, typename Sc, typename V_>
void
Sliver_perturber<C3T3,Md,Sc,V_>::
add_perturbation(Perturbation* perturbation)
{
  if ( !perturbation_vector_.empty() )
    perturbation_vector_.back().set_next(perturbation);

  if ( nullptr != perturbation )
  {
    // Set order
    perturbation->set_order(next_perturbation_order_++);

    // Add perturbation
    perturbation_vector_.push_back(perturbation);
  }
}




// -----------------------------------
// Private methods
// -----------------------------------
template <typename C3T3, typename Md, typename Sc, typename V_>
bool
Sliver_perturber<C3T3,Md,Sc,V_>::
perturb(const FT& sliver_bound, PQueue& pqueue, Visitor& visitor) const
{
#ifdef CGAL_MESH_3_PERTURBER_HIGH_VERBOSITY
  CGAL::Real_timer timer;
  timer.start();
  std::streamsize prec = std::cerr.precision(4);
  std::cerr << "Perturb sliver vertices (bound: " << sliver_bound
            << ") ..." << std::endl;
  std::cerr.precision(prec);
#endif

  // build priority queue
  int pqueue_size = build_priority_queue(sliver_bound, pqueue);

#ifdef CGAL_MESH_3_PERTURBER_HIGH_VERBOSITY
  std::cerr << "Legend of the following line: "
            << "(#vertices in pqueue, #iterations, #fails)" << std::endl;

  // Store construction time
  timer.stop();
  double construction_time = timer.time();
  timer.reset();
  timer.start();
#endif

  // Stores the vertices for which perturbation has failed
  Bad_vertices_vector bad_vertices;

#ifdef CGAL_LINKED_WITH_TBB
  // Parallel
  if (boost::is_convertible<Concurrency_tag, Parallel_tag>::value)
  {
    this->create_task_group();

    while (pqueue.size() > 0)
    {
      PVertex pv = pqueue.top();
      pqueue.pop();
      enqueue_task(pv, sliver_bound,
                   visitor, bad_vertices);
    }

    this->wait_for_all();

# if defined(CGAL_MESH_3_PERTURBER_VERBOSE) || defined(CGAL_MESH_3_PROFILING)
    std::cerr << " Flushing";
# endif
    bool keep_flushing = true;
    while (keep_flushing)
    {
      keep_flushing = this->flush_work_buffers();
# if defined(CGAL_MESH_3_PERTURBER_VERBOSE) || defined(CGAL_MESH_3_PROFILING)
      std::cerr << ".";
# endif
    }

    this->destroy_trask_group();
  }
  // Sequential
  else
#endif // CGAL_LINKED_WITH_TBB
  {
# ifdef CGAL_MESH_3_PERTURBER_VERBOSE
    int iteration_nb = 0;
# endif

    while ( !is_time_limit_reached() && !pqueue.empty() )
    {
      // Get pqueue head
      PVertex pv = pqueue.top();
      pqueue.pop();
      --pqueue_size;

      CGAL_assertion(pv.is_perturbable());

      // Get pvertex slivers list
# ifdef CGAL_NEW_INCIDENT_SLIVERS
      Cell_vector slivers;
      helper_.new_incident_slivers(pv.vertex(), sliver_criterion_, sliver_bound,
                                 std::back_inserter(slivers));
# else
      Cell_vector slivers =
        helper_.incident_slivers(pv.vertex(), sliver_criterion_, sliver_bound);
# endif

      CGAL_assertion(slivers.size() == pv.sliver_nb());

      // Perturb vertex
      Vertex_vector modified_vertices;

      // pv.perturbation() should not be nullptr if pv is in pqueue
      CGAL_assertion(pv.perturbation() != nullptr);

      std::pair<bool,Vertex_handle> perturbation_ok =
        pv.perturbation()->operator()(pv.vertex(),
                                      slivers,
                                      c3t3_,
                                      domain_,
                                      sliver_criterion_,
                                      sliver_bound,
                                      modified_vertices);

      // If vertex has changed - may happen in two cases: vertex has been moved
      // or vertex has been reverted to the same location -
      if ( perturbation_ok.second != pv.vertex() )
      {
        // Update pvertex vertex
        pv.set_vertex(perturbation_ok.second);
      }

      // If v has been moved
      if ( perturbation_ok.first )
      {
        // Update pvertex
        update_pvertex(pv,sliver_bound);

        // If pv needs to be modified again, try first perturbation
        pv.set_perturbation(&perturbation_vector_.front());
        pv.increment_try_nb();

        // update modified vertices
        pqueue_size += update_priority_queue(modified_vertices,
                                             sliver_bound,
                                             pqueue);
      }
      else
      {
        // If perturbation fails, try next one
        pv.set_perturbation(pv.perturbation()->next());

        if ( nullptr == pv.perturbation() )
        {
          bad_vertices.push_back(pv.vertex());
        }
      }

      // Update pqueue in every cases, because pv was poped
      pqueue_size += update_priority_queue(pv, pqueue);
      visitor.end_of_perturbation_iteration(pqueue_size);

# ifdef CGAL_MESH_3_PERTURBER_HIGH_VERBOSITY
      ++iteration_nb;
      std::cerr << boost::format("\r             \r"
                                 "(%1%,%2%,%4%) (%|3$.1f| iteration/s)")
      % pqueue_size
      % iteration_nb
      % (iteration_nb / timer.time())
      % bad_vertices.size();
# endif

# ifdef CGAL_MESH_3_PERTURBER_LOW_VERBOSITY
      ++iteration_nb;
      std::cerr << boost::format("\r             \r"
                                 "bound %5%: (%1%,%2%,%4%) (%|3$.1f| iteration/s)")
      % pqueue_size
      % iteration_nb
      % (iteration_nb / running_time_.time())
      % bad_vertices.size()
      % sliver_bound;
# endif
    }
  }

#ifdef CGAL_MESH_3_PERTURBER_HIGH_VERBOSITY
  std::cerr << std::endl;
  print_perturbations_statistics();
  std::cerr << "Step perturbation time: " << timer.time() + construction_time
            << "s" << std::endl << std::endl;
#endif

  if ( is_time_limit_reached() )
    return false;

  // update bad vertices list (remove those which are not bad anymore)
  update_bad_vertices(bad_vertices,sliver_bound);
  return bad_vertices.empty();
}


#ifdef CGAL_FASTER_BUILD_QUEUE

template <typename C3T3, typename Md, typename Sc, typename V_>
int
Sliver_perturber<C3T3,Md,Sc,V_>::
build_priority_queue(const FT& sliver_bound, PQueue& pqueue) const
{
  CGAL_precondition(pqueue.empty());

#ifdef CGAL_MESH_3_PERTURBER_HIGH_VERBOSITY
  CGAL::Real_timer timer;
  timer.start();
  std::cerr << "Build pqueue...";
#endif

  int pqueue_size = 0;

  typedef CGAL::Hash_handles_with_or_without_timestamps            Hash_fct;
  typedef boost::unordered_map<Vertex_handle, PVertex, Hash_fct>   M;

  M vpm;
  for ( typename Tr::Finite_cells_iterator cit = tr_.finite_cells_begin();
       cit != tr_.finite_cells_end() ;
       ++cit )
  {
    if(helper_.is_sliver(cit, sliver_criterion_, sliver_bound))
    {
      double d = cit->sliver_value();
      for(int i=0; i< 4; i++){
        Vertex_handle vh = cit->vertex(i);
        PVertex& pv = vpm[vh];
        if(pv.sliver_nb() ==0)
        {
          pv.set_vertex(vh);
          pv.set_id( get_pvertex_id(vh));
          pv.set_sliver_nb(1);
          pv.set_min_value(d);
          pv.set_perturbation(&perturbation_vector_.front());
        }
        else
        {
          pv.set_sliver_nb(pv.sliver_nb()+1);
          if(d < pv.min_value())
            pv.set_min_value(d);
        }
      }
    }
  }

  for( typename M::iterator vit = vpm.begin();
       vit != vpm.end() ;
       ++vit )
    pqueue_size += update_priority_queue(vit->second, pqueue);

#ifdef CGAL_MESH_3_PERTURBER_HIGH_VERBOSITY
  std::cerr << "done (" << pqueue_size << " vertices inserted in "
            << timer.time() << "s)\n";
#endif

  return pqueue_size;
}

#else // not CGAL_FASTER_BUILD_QUEUE

template <typename C3T3, typename Md, typename Sc, typename V_>
int
Sliver_perturber<C3T3,Md,Sc,V_>::
build_priority_queue(const FT& sliver_bound, PQueue& pqueue) const
{
  CGAL_precondition(pqueue.empty());

#ifdef CGAL_MESH_3_PERTURBER_HIGH_VERBOSITY
  CGAL::Real_timer timer;
  timer.start();
  std::cerr << "Build pqueue...";
#endif

  int pqueue_size = 0;

  for ( typename Tr::Finite_vertices_iterator vit = tr_.finite_vertices_begin();
       vit != tr_.finite_vertices_end() ;
       ++vit )
  {
    PVertex pv = make_pvertex(vit, sliver_bound, get_pvertex_id(vit));
    pqueue_size += update_priority_queue(pv, pqueue);
  }

#ifdef CGAL_MESH_3_PERTURBER_HIGH_VERBOSITY
  std::cerr << "done (" << pqueue_size << " vertices inserted in "
            << timer.time() << "s)\n";
#endif

  return pqueue_size;
}
#endif // not CGAL_FASTER_BUILD_QUEUE


template <typename C3T3, typename Md, typename Sc, typename V_>
int
Sliver_perturber<C3T3,Md,Sc,V_>::
update_priority_queue(const Vertex_vector& vertices,
                      const FT& sliver_bound,
                      PQueue& pqueue) const
{
  int modified_pv_nb = 0;
  for ( typename Vertex_vector::const_iterator vit = vertices.begin() ;
       vit != vertices.end() ;
       ++vit )
  {
    PVertex pv = make_pvertex(*vit,sliver_bound,get_pvertex_id(*vit));
    modified_pv_nb += update_priority_queue(pv, pqueue);
  }

  return modified_pv_nb;
}


#ifdef CGAL_LINKED_WITH_TBB
// For parallel version
template <typename C3T3, typename Md, typename Sc, typename V_>
int
Sliver_perturber<C3T3,Md,Sc,V_>::
update_priority_queue( const Vertex_vector& vertices
                     , const FT& sliver_bound
                     , Visitor& visitor
                     , Bad_vertices_vector &bad_vertices) const
{
  int modified_pv_nb = 0;
  for ( typename Vertex_vector::const_iterator vit = vertices.begin() ;
       vit != vertices.end() ;
       ++vit )
  {
    PVertex pv = make_pvertex__concurrent(*vit,sliver_bound,get_pvertex_id(*vit));
    if (pv.is_perturbable())
    {
      enqueue_task(pv, sliver_bound, visitor, bad_vertices);
      ++modified_pv_nb;
    }
  }

  return modified_pv_nb;
}
#endif

template <typename C3T3, typename Md, typename Sc, typename V_>
int
Sliver_perturber<C3T3,Md,Sc,V_>::
update_priority_queue(const PVertex& pv, PQueue& pqueue) const
{
  if ( pqueue.contains(pv) )
  {
    if ( pv.is_perturbable() )
    {
      pqueue.update(pv);
      return 0;
    }
    else
    {
      pqueue.remove(pv);
      return -1;
    }
  }
  else
  {
    if ( pv.is_perturbable() )
    {
      pqueue.push(pv);
      return 1;
    }
  }

  return 0;
}

#ifdef CGAL_LINKED_WITH_TBB
// For parallel version
template <typename C3T3, typename Md, typename Sc, typename V_>
void
Sliver_perturber<C3T3,Md,Sc,V_>::
perturb_vertex( PVertex pv
              , const FT& sliver_bound
              , Visitor& visitor
              , Bad_vertices_vector &bad_vertices
              , bool *could_lock_zone
              ) const
{
  typename Gt::Construct_point_3 cp = tr_.geom_traits().construct_point_3_object();

#ifdef CGAL_CONCURRENT_MESH_3_PROFILING
  static Profile_branch_counter_3 bcounter(
    "early withdrawals / late withdrawals / successes [Perturber]");
#endif

  *could_lock_zone = true;

  // Zombie?
  if (pv.is_zombie())
  {
    return;
  }

  Bare_point p = cp(pv.vertex()->point());
  if (!helper_.try_lock_point_no_spin(pv.vertex()->point()) ||
      ! tr_.geom_traits().equal_3_object()(p, cp(pv.vertex()->point())))
  {
#ifdef CGAL_CONCURRENT_MESH_3_PROFILING
    bcounter.increment_branch_2(); // THIS is an early withdrawal!
#endif
    *could_lock_zone = false;
    return;
  }

  // Zombie? (in case the vertex has changed in the meantime)
  if (pv.is_zombie())
  {
    return;
  }

  CGAL_assertion(pv.is_perturbable());

  int num_new_vertices_to_treat = 0;

  Cell_vector slivers;
  slivers.reserve(8);
  if (!helper_.try_lock_and_get_incident_slivers(
    pv.vertex(), sliver_criterion_, sliver_bound, slivers))
  {
    *could_lock_zone = false;
#ifdef CGAL_CONCURRENT_MESH_3_PROFILING
    bcounter.increment_branch_1(); // THIS is a late withdrawal!
#endif
  }
  else
  {
    // Slivers may be empty if the vertex has been modified by another thread in the meatime
    if (slivers.empty())
    {
      return;
    }

    // Perturb vertex
    Vertex_vector modified_vertices;

    // pv.perturbation() should not be nullptr if pv is in pqueue
    CGAL_assertion(pv.perturbation() != nullptr);

    std::pair<bool,Vertex_handle> perturbation_ok =
      pv.perturbation()->operator()(pv.vertex(),
                                    slivers,
                                    c3t3_,
                                    domain_,
                                    sliver_criterion_,
                                    sliver_bound,
                                    modified_vertices,
                                    could_lock_zone);

    if (*could_lock_zone)
    {
      // If vertex has changed - may happen in two cases: vertex has been moved
      // or vertex has been reverted to the same location -
      if ( perturbation_ok.second != pv.vertex() )
      {
        // Update pvertex vertex
        pv.set_vertex(perturbation_ok.second);
      }
      // If the vertex hasn't changed, we still need to "virtually" increment
      // the erase counter, because we need to invalidate the PVertex that
      // may be in other threads' queues
      else
      {
        this->increment_erase_counter(pv.vertex());
      }

      // If v has been moved
      if ( perturbation_ok.first )
      {
        // Update pvertex
        update_pvertex__concurrent(pv,sliver_bound);

        // If pv needs to be modified again, try first perturbation
        pv.set_perturbation(&perturbation_vector_.front());
        pv.increment_try_nb();

        // update modified vertices
        num_new_vertices_to_treat +=
          update_priority_queue(modified_vertices, sliver_bound, visitor, bad_vertices);
      }
      else
      {
        // If perturbation fails, try next one
        pv.set_perturbation(pv.perturbation()->next());
        pv.update_saved_erase_counter();

        if ( nullptr == pv.perturbation() )
        {
          bad_vertices.push_back(pv.vertex());
        }
      }

#ifdef CGAL_CONCURRENT_MESH_3_PROFILING
      ++bcounter;
#endif

      // Update pqueue in every cases, because pv was poped
      if (pv.is_perturbable())
      {
        enqueue_task(pv, sliver_bound, visitor, bad_vertices);
        ++num_new_vertices_to_treat;
      }
    }
    else
    {
#ifdef CGAL_CONCURRENT_MESH_3_PROFILING
      bcounter.increment_branch_1();   // THIS is a late withdrawal!
#endif
    }
  }

  visitor.end_of_perturbation_iteration(0);
}
#endif


// Sequential
template <typename C3T3, typename Md, typename Sc, typename V_>
typename Sliver_perturber<C3T3,Md,Sc,V_>::PVertex
Sliver_perturber<C3T3,Md,Sc,V_>::
make_pvertex(const Vertex_handle& vh,
             const FT& sliver_bound,
             const typename PVertex::id_type& pv_id) const
{
  CGAL_assertion(!tr_.is_infinite(vh));

  // Make pvertex in all cases
  PVertex pv(vh,pv_id);
  pv.set_perturbation(&perturbation_vector_.front());
  update_pvertex(pv, sliver_bound);

  return pv;
}

// Parallel
template <typename C3T3, typename Md, typename Sc, typename V_>
typename Sliver_perturber<C3T3,Md,Sc,V_>::PVertex
Sliver_perturber<C3T3,Md,Sc,V_>::
make_pvertex__concurrent(
             const Vertex_handle& vh,
             const FT& sliver_bound,
             const typename PVertex::id_type& pv_id) const
{
  // Make pvertex in all cases
  PVertex pv(vh,pv_id);
  pv.set_perturbation(&perturbation_vector_.front());
  update_pvertex__concurrent(pv, sliver_bound);

  return pv;
}


// Sequential
template <typename C3T3, typename Md, typename Sc, typename V_>
void
Sliver_perturber<C3T3,Md,Sc,V_>::
update_pvertex(PVertex& pv, const FT& sliver_bound) const
{
#ifdef CGAL_NEW_INCIDENT_SLIVERS
  Cell_vector slivers;
  helper_.new_incident_slivers(pv.vertex(), sliver_criterion_, sliver_bound, std::back_inserter(slivers));
#else
  Cell_vector slivers =
    helper_.incident_slivers(pv.vertex(), sliver_criterion_, sliver_bound);
#endif

  pv.set_sliver_nb(static_cast<unsigned int>(slivers.size()));
  pv.set_min_value(helper_.min_sliver_value(slivers, sliver_criterion_));
}

#ifdef CGAL_LINKED_WITH_TBB
// Parallel
template <typename C3T3, typename Md, typename Sc, typename V_>
void
Sliver_perturber<C3T3,Md,Sc,V_>::
update_pvertex__concurrent(PVertex& pv, const FT& sliver_bound) const
{
  Cell_vector slivers;
  helper_.get_incident_slivers_without_using_tds_data(
    pv.vertex(), sliver_criterion_, sliver_bound, slivers);

  pv.set_sliver_nb(static_cast<unsigned int>(slivers.size()));
  pv.set_min_value(helper_.min_sliver_value(slivers, sliver_criterion_));
}
#endif

// Sequential
template <typename C3T3, typename Md, typename Sc, typename V_>
void
Sliver_perturber<C3T3,Md,Sc,V_>::
update_bad_vertices(std::vector<Vertex_handle> &bad_vertices,
                    const FT& sliver_bound) const
{
  typename std::vector<Vertex_handle>::iterator vit = bad_vertices.begin();
  while ( vit != bad_vertices.end() )
  {
    if ( tr_.is_vertex(*vit)
        && helper_.min_incident_value(*vit,sliver_criterion_) <= sliver_bound )
    {
      ++vit;
    }
    else
    {
      vit = bad_vertices.erase(vit);
    }
  }
}

#ifdef CGAL_LINKED_WITH_TBB
// Parallel
template <typename C3T3, typename Md, typename Sc, typename V_>
void
Sliver_perturber<C3T3,Md,Sc,V_>::
update_bad_vertices(tbb::concurrent_vector<Vertex_handle> &bad_vertices,
                    const FT& sliver_bound) const
{
  tbb::concurrent_vector<Vertex_handle> tmpv;

  typename tbb::concurrent_vector<Vertex_handle>::iterator vit
    = bad_vertices.begin();
  while ( vit != bad_vertices.end() )
  {
    if ( tr_.is_vertex(*vit)
        && helper_.min_incident_value(*vit,sliver_criterion_) <= sliver_bound )
    {
      tmpv.push_back(*vit);
    }
    ++vit;
  }
  bad_vertices.swap(tmpv);
}
#endif // CGAL_LINKED_WITH_TBB


template <typename C3T3, typename Md, typename Sc, typename V_>
void
Sliver_perturber<C3T3,Md,Sc,V_>::
initialize_vertices_id() const
{
  int cur_id = 0;
  for(typename Tr::Finite_vertices_iterator it = tr_.finite_vertices_begin();
      it != tr_.finite_vertices_end(); ++it) {
    it->set_meshing_info(cur_id++);
  }
}


#ifdef CGAL_MESH_3_PERTURBER_VERBOSE
template <typename C3T3, typename Md, typename Sc, typename V_>
void
Sliver_perturber<C3T3,Md,Sc,V_>::
print_perturbations_statistics() const
{
  int total_perturbation_nb = 0;
  typename Perturbation_vector::const_iterator it = perturbation_vector_.begin();
  for ( ; it != perturbation_vector_.end() ; ++it )
  {
    total_perturbation_nb += it->counter();
  }

  if ( 0 == total_perturbation_nb )
  {
    std::cerr << "No perturbation done at this step" << std::endl;
    return;
  }

  for ( it = perturbation_vector_.begin() ;
       it != perturbation_vector_.end() ;
       ++it )
  {
    std::cerr << it->perturbation_name() << ": "
              << (double)it->counter() / (double)total_perturbation_nb * 100.
              << "% (" << it->counter() << " in " << it->time() << "s)"
              << std::endl;
  }
}



template <typename C3T3, typename Md, typename Sc, typename V_>
void
Sliver_perturber<C3T3,Md,Sc,V_>::
print_final_perturbations_statistics() const
{
  int total_perturbation_nb = 0;
  typename Perturbation_vector::const_iterator it = perturbation_vector_.begin();
  for ( ; it != perturbation_vector_.end() ; ++it )
  {
    total_perturbation_nb += it->total_counter();
  }

  if ( 0 == total_perturbation_nb )
  {
    std::cerr << "No perturbation done" << std::endl;
    return;
  }

  for ( it = perturbation_vector_.begin() ;
       it != perturbation_vector_.end() ;
       ++it )
  {
    std::cerr << it->perturbation_name() << ": "
              << (double)it->total_counter() / (double)total_perturbation_nb * 100.
              << "% (" << it->total_counter() << " in " << it->total_time() << "ms)"
              << std::endl;
  }
}



template <typename C3T3, typename Md, typename Sc, typename V_>
void
Sliver_perturber<C3T3,Md,Sc,V_>::
reset_perturbation_counters()
{
  typename Perturbation_vector::iterator it = perturbation_vector_.begin();
  for ( ; it != perturbation_vector_.end() ; ++it )
  {
    it->reset_counter();
    it->reset_timer();
  }
}
#endif


#ifdef CGAL_LINKED_WITH_TBB
// For parallel version
template <typename C3T3, typename Md, typename Sc, typename V_>
void
Sliver_perturber<C3T3,Md,Sc,V_>::
enqueue_task(const PVertex &pv,
             const FT& sliver_bound,
             Visitor& visitor,
             Bad_vertices_vector &bad_vertices
             ) const
{
  this->enqueue_work(
    Perturb_vertex<Self, Bad_vertices_vector>(
      *this, pv, sliver_bound, visitor, bad_vertices),
    pv);
}
#endif

} // end namespace Mesh_3


} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_MESH_3_SLIVERS_PERTURBER_H
