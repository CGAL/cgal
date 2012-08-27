// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
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
// Author(s)     : Stephane Tayeb
//
//******************************************************************************
// File Description : 
//******************************************************************************

#ifndef CGAL_MESH_3_SLIVER_PERTURBER_H
#define CGAL_MESH_3_SLIVER_PERTURBER_H


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
#include <CGAL/Timer.h>
#include <CGAL/Mesh_3/Null_perturber_visitor.h>

#ifdef CGAL_MESH_3_USE_RELAXED_HEAP
#include <boost/pending/relaxed_heap.hpp>
#else
#include <CGAL/Modifiable_priority_queue.h>
#endif //CGAL_MESH_3_USE_RELAXED_HEAP
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/ptr_container/ptr_vector.hpp>


namespace CGAL {

namespace Mesh_3 {
  
template < typename C3T3,
           typename MeshDomain,
           typename SliverCriterion,
           typename Visitor_ = Null_perturber_visitor<C3T3> >
class Sliver_perturber
{
  typedef typename C3T3::Triangulation  Tr;
  typedef typename Tr::Geom_traits      Gt;
  
  typedef typename Tr::Cell_handle              Cell_handle;
  typedef typename Tr::Vertex_handle            Vertex_handle;
  typedef typename Tr::Vertex                   Vertex;
  
  typedef typename std::vector<Cell_handle>     Cell_vector;
  typedef typename std::vector<Vertex_handle>   Vertex_vector;
  
  typedef typename Gt::FT                       FT;
  
  typedef Abstract_perturbation<C3T3,MeshDomain,SliverCriterion> Perturbation;
  typedef boost::ptr_vector<Perturbation>                 Perturbation_vector;

  // Helper
  typedef class C3T3_helpers<C3T3,MeshDomain> C3T3_helpers;
  
  // Visitor
  // Should define
  //  - bound_reached(FT bound)
  //  - end_of_perturbation_iteration(std::size_t vertices_left)
  typedef Visitor_ Visitor;
  
private:
  // Relaxed heap
  // -----------------------------------
  // Private classes
  // -----------------------------------
  /**
   * @class PVertex
   * Vertex with associated perturbation datas
   */
  class PVertex
  {
  public:
    typedef unsigned int id_type;
    
    /// Constructor
    PVertex(const Vertex_handle& vh, id_type id)
    : vertex_handle_(vh)
    , incident_sliver_nb_(0)
    , min_value_(SliverCriterion::max_value)
    , try_nb_(0)
    , p_perturbation_(NULL)
    , id_(id) { }
    
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
              && (NULL != perturbation()) 
              && (sliver_nb() != 0) );
    }
    
    /// Min sliver value
    const FT& min_value() const { return min_value_; }
    void set_min_value(const FT& min_value) { min_value_ = min_value; }
    
    /// Try nb
    const unsigned int& try_nb() const { return try_nb_; }
    void set_try_nb(const unsigned int& try_nb) { try_nb_ = try_nb; }
    void increment_try_nb() { ++try_nb_; }
    
    /// Id
    id_type id() const { return id_; }
    
    /// Operators
    bool operator==(const PVertex& pv) const { return ( id() == pv.id() ); }
    
    bool operator<(const PVertex& pv) const
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
      else
        return false; // all characteristics are the same!
    }
    
  private:
    /// Private datas
    Vertex_handle vertex_handle_;
    unsigned int incident_sliver_nb_;
    FT min_value_;
    unsigned int try_nb_;
    const Perturbation* p_perturbation_;
    id_type id_;
  };
  
  
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
                   const SliverCriterion& criterion = SliverCriterion());
  
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
  operator()(const FT& sliver_bound = SliverCriterion::max_value,
             const FT& delta = FT(1.),
             Visitor visitor = Visitor());
  
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
  int update_priority_queue(const Vertex_vector& vertices,
                            const FT& sliver_bound,
                            PQueue& pqueue) const;
  
  /**
   * Updates \c pv in priority queue
   */
  int update_priority_queue(const PVertex& pv, PQueue& pqueue) const;
  
  /**
   * Returns a pvertex from a vertex handle \c vh, using id \c pv_id
   */
  PVertex
  make_pvertex(const Vertex_handle& vh,
               const FT& sliver_bound,
               const typename PVertex::id_type& pv_id) const;
  
  /**
   * Updates a pvertex \c pv
   */
  void update_pvertex(PVertex& pv, const FT& sliver_bound) const;
  
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
  void update_bad_vertices(Vertex_vector& bad_vertices,
                           const FT& sliver_bound) const;
  
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
  
private:
  // -----------------------------------
  // Private data
  // -----------------------------------
  C3T3& c3t3_;
  Tr& tr_;
  const MeshDomain& domain_;
  double sliver_bound_;
  Perturbation_vector perturbation_vector_;
  SliverCriterion sliver_criterion_;
  C3T3_helpers helper_;
  
  // Internal perturbation ordering
  int next_perturbation_order_;
  
  // Timer
  double time_limit_;
  CGAL::Timer running_time_;
};
  
  
  
  
  
  
  
template <typename C3T3, typename Md, typename Sc, typename V_>
Sliver_perturber<C3T3,Md,Sc,V_>::
Sliver_perturber(C3T3& c3t3,
                 const Md& domain,
                 const Sc& criterion)
  : c3t3_(c3t3)
  , tr_(c3t3_.triangulation())
  , domain_(domain)
  , sliver_criterion_(criterion)
  , helper_(c3t3_,domain_)
  , next_perturbation_order_(0)
  , time_limit_(-1)
  , running_time_()
{
}
  
  
  
template <typename C3T3, typename Md, typename Sc, typename V_>
Mesh_optimization_return_code
Sliver_perturber<C3T3,Md,Sc,V_>::
operator()(const FT& sliver_bound, const FT& delta, Visitor visitor)
{
  // Reset sliver value cache
  helper_.reset_cache();
  
  // Init time counter
  running_time_.reset();
  running_time_.start();
  
  // Build priority queue (we use one queue for all steps)
  PQueue pqueue(tr_.number_of_vertices());

  // Initialize vertices ids
  initialize_vertices_id();
  
  
#ifdef CGAL_MESH_3_PERTURBER_VERBOSE
  std::cerr << "Running sliver perturbation..." << std::endl;
#endif
  
#ifdef CGAL_MESH_3_PERTURBER_LOW_VERBOSITY
  std::cerr << "Legend of the following line: "
            << "(#vertices in pqueue, #iterations, #fails)" << std::endl;
#endif
  
  FT current_bound = delta;
  bool perturbation_ok = true;
  while ( current_bound <= sliver_bound && perturbation_ok)
  {
#ifdef CGAL_MESH_3_PERTURBER_HIGH_VERBOSITY
    // reset_perturbation_counters is not const
    reset_perturbation_counters();
#endif     
    perturbation_ok = perturb(current_bound, pqueue, visitor);
    
    visitor.bound_reached(current_bound);
    
    current_bound += delta;
    if ( (current_bound >= sliver_bound)
        && (current_bound < sliver_bound + delta) )
    { 
      current_bound = sliver_bound;
    }
  }
  
  running_time_.stop();
  
#ifdef CGAL_MESH_3_PERTURBER_VERBOSE
  std::cerr << std::endl
            << "Total perturbation time: " << running_time_.time() << "s";
  std::cerr << std::endl << "Perturbation statistics:" << std::endl;
  print_final_perturbations_statistics();
  std::cerr << std::endl;
#endif
  
  if ( is_time_limit_reached() )
    return TIME_LIMIT_REACHED;
  
  if ( !perturbation_ok )
    return CANT_IMPROVE_ANYMORE;
  
  return BOUND_REACHED;
}

  
  
template <typename C3T3, typename Md, typename Sc, typename V_>
void
Sliver_perturber<C3T3,Md,Sc,V_>::
add_perturbation(Perturbation* perturbation)
{
  if ( !perturbation_vector_.empty() )
    perturbation_vector_.back().set_next(perturbation);
  
  if ( NULL != perturbation )
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
  CGAL::Timer timer;
  timer.start();
  std::cerr.precision(4);
  std::cerr << "Perturb sliver vertices (bound: " << sliver_bound 
            << ") ..." << std::endl;
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

#ifdef CGAL_MESH_3_PERTURBER_VERBOSE
  int iteration_nb = 0;
#endif
  
  // Stores the vertices for which perturbation has failed
  Vertex_vector bad_vertices;
  
  while ( !is_time_limit_reached() && !pqueue.empty() )
  {
    // Get pqueue head
    PVertex pv = pqueue.top();
    pqueue.pop();
    --pqueue_size;
    
    CGAL_assertion(pv.is_perturbable());
    
    // Get pvertex slivers list
    Cell_vector slivers =
      helper_.incident_slivers(pv.vertex(), sliver_criterion_, sliver_bound);
    CGAL_assertion(!slivers.empty());
    
    // Perturb vertex
    Vertex_vector modified_vertices;
    
    // pv.perturbation() should not be NULL if pv is in pqueue 
    CGAL_assertion(pv.perturbation() != NULL);
    
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
      
      if ( NULL == pv.perturbation() )
      {
        bad_vertices.push_back(pv.vertex());
      }
    }
    
    // Update pqueue in every cases, because pv was poped
    pqueue_size += update_priority_queue(pv, pqueue);
    visitor.end_of_perturbation_iteration(pqueue_size);
    
#ifdef CGAL_MESH_3_PERTURBER_HIGH_VERBOSITY
    ++iteration_nb;
    std::cerr << boost::format("\r             \r"
                               "(%1%,%2%,%4%) (%|3$.1f| iteration/s)")
    % pqueue_size
    % iteration_nb
    % (iteration_nb / timer.time())
    % bad_vertices.size();
#endif
    
#ifdef CGAL_MESH_3_PERTURBER_LOW_VERBOSITY
    ++iteration_nb;
    std::cerr << boost::format("\r             \r"
                               "bound %5%: (%1%,%2%,%4%) (%|3$.1f| iteration/s)")
    % pqueue_size
    % iteration_nb
    % (iteration_nb / running_time_.time())
    % bad_vertices.size()
    % sliver_bound;
#endif
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
  
  
  
  
template <typename C3T3, typename Md, typename Sc, typename V_>
int
Sliver_perturber<C3T3,Md,Sc,V_>::
build_priority_queue(const FT& sliver_bound, PQueue& pqueue) const
{
  CGAL_precondition(pqueue.empty());
  
#ifdef CGAL_MESH_3_PERTURBER_HIGH_VERBOSITY
  CGAL::Timer timer;
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


template <typename C3T3, typename Md, typename Sc, typename V_>
typename Sliver_perturber<C3T3,Md,Sc,V_>::PVertex
Sliver_perturber<C3T3,Md,Sc,V_>::
make_pvertex(const Vertex_handle& vh,
             const FT& sliver_bound,
             const typename PVertex::id_type& pv_id) const
{
  // Make pvertex in all cases
  PVertex pv(vh,pv_id);
  pv.set_perturbation(&perturbation_vector_.front());
  update_pvertex(pv, sliver_bound);
  
  return pv;
}

  
  
template <typename C3T3, typename Md, typename Sc, typename V_>
void
Sliver_perturber<C3T3,Md,Sc,V_>::
update_pvertex(PVertex& pv, const FT& sliver_bound) const
{
  Cell_vector slivers =
    helper_.incident_slivers(pv.vertex(), sliver_criterion_, sliver_bound);
  
  pv.set_sliver_nb(static_cast<unsigned int>(slivers.size()));
  pv.set_min_value(helper_.min_sliver_value(slivers, sliver_criterion_));
}
  
  
template <typename C3T3, typename Md, typename Sc, typename V_>
void
Sliver_perturber<C3T3,Md,Sc,V_>::
update_bad_vertices(Vertex_vector& bad_vertices,
                    const FT& sliver_bound) const
{
  typename Vertex_vector::iterator vit = bad_vertices.begin();
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
  
  
template <typename C3T3, typename Md, typename Sc, typename V_>
void
Sliver_perturber<C3T3,Md,Sc,V_>::
initialize_vertices_id() const
{
  namespace bl = boost::lambda;
  int cur_id = 0; 
  
  std::for_each(tr_.finite_vertices_begin(), tr_.finite_vertices_end(),
                bl::bind(&Vertex::set_meshing_info, &bl::_1, bl::var(cur_id)++));
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
              << "% (" << it->total_counter() << " in " << it->total_time() << "s)"
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
  
  
  
} // end namespace Mesh_3


} //namespace CGAL

#endif // CGAL_MESH_3_SLIVERS_PERTURBER_H
