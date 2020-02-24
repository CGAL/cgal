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

#ifndef CGAL_MESH_3_VERTEX_PERTURBATION_H
#define CGAL_MESH_3_VERTEX_PERTURBATION_H

#include <CGAL/license/Mesh_3.h>

#include <CGAL/disable_warnings.h>
#include <CGAL/Mesh_3/config.h>

#include <CGAL/Mesh_3/C3T3_helpers.h>
#include <CGAL/Mesh_3/Triangulation_helpers.h>
#include <CGAL/Time_stamper.h>

#ifdef CGAL_MESH_3_PERTURBER_VERBOSE
  #include <CGAL/Timer.h>
  #ifdef CGAL_LINKED_WITH_TBB
    #include <tbb/enumerable_thread_specific.h>
    #include <atomic>
  #endif
#endif

#include <boost/optional.hpp>
#include <boost/random/lagged_fibonacci.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

#include <string>
#include <ctime>

namespace CGAL {

namespace Mesh_3 {

namespace details
{
/**
   * @brief Returns the angle in radian of vectors \c u and \c v
   */
template <typename K>
typename K::FT
angle_in_radian(const typename K::Vector_3& u,
                const typename K::Vector_3& v,
                K k = K())
{
  typedef typename K::FT FT;
  typedef typename K::Vector_3 Vector_3;

  typename K::Construct_cross_product_vector_3 cross_product =
    k.construct_cross_product_vector_3_object();
  typename K::Compute_scalar_product_3 scalar_product =
    k.compute_scalar_product_3_object();
  typename K::Compute_squared_length_3 sq_length =
    k.compute_squared_length_3_object();

  // -------------------------------------
  // Angle between two vectors (in rad)
  // uv = |u||v| cos(u,v)
  // u^v  = w
  // |w| = |u||v| |sin(u,v)|
  // -------------------------------------
  FT product = CGAL::sqrt(sq_length(u) * sq_length(v));

  // Check
  if ( product == FT(0) )
    return FT(0);

  // Sine
  Vector_3 w = cross_product(u,v);
  FT abs_sin = CGAL::sqrt(sq_length(w)) / product;

  if ( abs_sin < FT(-1) ) { abs_sin = FT(-1); }
  if ( abs_sin > FT(1) ) { abs_sin = FT(1); }
  CGAL_assertion(abs_sin >= -1);
  CGAL_assertion(abs_sin <= 1);

  // We just need cosine sign
  FT cosine_sign = scalar_product(u,v);

  if ( cosine_sign >= FT(0) )
    return FT(std::asin(abs_sin));
  else
    return FT(CGAL_PI) - FT(std::asin(abs_sin));
}

/**
   * @brief Returns the angle in radian of vectors \c u and \c v
   */
template <typename Vector_3>
typename Kernel_traits<Vector_3>::Kernel::FT
angle_in_radian(const Vector_3& u, const Vector_3& v)
{
  return angle_in_radian(u,v,typename Kernel_traits<Vector_3>::Kernel());
}

/**
   * @brief Returns the squared length of edge \c e
   */
template <typename Tr>
typename Tr::Geom_traits::FT
edge_sq_length(const typename Tr::Edge& e,
               const Tr& tr)
{
  typedef typename Tr::Geom_traits     Gt;
  typedef typename Tr::Bare_point      Bare_point;
  typedef typename Tr::Weighted_point  Weighted_point;

  typename Gt::Construct_point_3 cp =
    tr.geom_traits().construct_point_3_object();
  typename Gt::Compute_squared_distance_3 sq_distance =
    tr.geom_traits().compute_squared_distance_3_object();

  const Weighted_point& wp = tr.point(e.first, e.second);
  const Weighted_point& wq = tr.point(e.first, e.third);
  const Bare_point& p = cp(wp);
  const Bare_point& q = cp(wq);

  return sq_distance(p,q);
}

/**
   * @brief Returns the minimal incident edge length of \c v
   * in triangulation \c tr
   */
template <typename Tr>
typename Tr::Geom_traits::FT
min_incident_edge_sq_length(const typename Tr::Vertex_handle& v,
                            const Tr& tr)
{
  CGAL_precondition(!tr.is_infinite(v));

  typedef typename Tr::Edge Edge;
  typedef typename Tr::Geom_traits::FT FT;

  // Get all incident edges
  std::vector<Edge> edges;
  tr.finite_incident_edges(v, std::back_inserter(edges));
  CGAL_assertion(!edges.empty());

  // Get squared min length
  typename std::vector<Edge>::iterator eit = edges.begin();
  FT min_sq_length = edge_sq_length<Tr>(*eit++, tr);

  for ( ; eit != edges.end() ; ++eit )
  {
    min_sq_length = (std::min)(min_sq_length, edge_sq_length<Tr>(*eit, tr));
  }

  return min_sq_length;
}

} // end namespace details

/**
 * @class Abstract_perturbation
 *
 * Perturbation interface. It is used by Sliver_perturber class.
 */
template <typename C3T3, typename MeshDomain, typename SliverCriterion>
class Abstract_perturbation
{
protected:
  typedef typename C3T3::Vertex_handle                    Vertex_handle;
  typedef typename C3T3::Cell_handle                      Cell_handle;

  typedef typename C3T3::Triangulation::Geom_traits::FT   FT;

public:
  /**
   * @brief constructor
   */
  Abstract_perturbation()
    : p_next_(nullptr)
    , p_previous_(nullptr)
    , order_(0)
#ifdef CGAL_MESH_3_PERTURBER_VERBOSE
    , counter_(0)
    , timer_()
#endif
  {
#ifdef CGAL_MESH_3_PERTURBER_VERBOSE
    // Initialized here in case it's some std::atomic
    total_counter_ = 0;
    total_time_ = 0;
#endif
  }

  /**
   * @brief destructor
   *
   * Note that Abstract_perturbation is not responsible of
   * p_next_ nor p_previous_ deletion
   */
  virtual ~Abstract_perturbation() {}

  /**
   * @brief This operator try to move vertex \v using the perturbation.
   * @param v the vertex to move
   * @param slivers a vector which contains incident slivers of \c v
   * @param c3t3 the c3t3
   * @param domain the domain
   * @param criterion the criterion which is used to evaluate if a cell is
   *   a sliver.
   * @param sliver_bound the bound for the above criterion.
   * @param modified_vertices an output vector which contains vertices of c3t3
   *   which may have been impacted by v relocation.
   * @return a pair containing:
   *   - a bool which is \c true if a move has been done.
   *   - a Vertex_handle which is always filled and may be the new vertex (if
   *   the move is a success), or the vertex which lies at \c v's position in
   *   the new c3t3.
   *
   * Note that this function is hill_climbing only. The min \c criterion value
   * of c3t3 could not decrease.
   */
  std::pair<bool,Vertex_handle>
  operator()(const Vertex_handle& v,
             const std::vector<Cell_handle>& slivers,
             C3T3& c3t3,
             const MeshDomain& domain,
             const SliverCriterion& criterion,
             const FT& sliver_bound,
             std::vector<Vertex_handle>& modified_vertices,
             bool *could_lock_zone = nullptr) const
  {
#ifndef CGAL_MESH_3_PERTURBER_VERBOSE
    return do_perturb(v, slivers, c3t3, domain, criterion,
                      sliver_bound, modified_vertices, could_lock_zone);
#else
    timer().start();

    // Virtual call
    std::pair<bool,Vertex_handle> perturb =
      do_perturb(v, slivers, c3t3, domain, criterion,
                 sliver_bound, modified_vertices, could_lock_zone);

    if ( perturb.first )
      ++counter_;

    timer().stop();

    return perturb;
#endif
  }

  /**
   * @brief Sets next perturbation
   */
  void set_next(Abstract_perturbation* next)
  {
    p_next_ = next;

    if ( nullptr != next )
      next->p_previous_ = this;
  }

  /**
   * Returns next perturbation
   */
  Abstract_perturbation* next() const
  {
    return p_next_;
  }

  /**
   * Returns previous perturbation
   */
  Abstract_perturbation* previous() const
  {
    return p_previous_;
  }

  /**
   * Returns the order value
   */
  int order() const
  {
    return order_;
  }

  /**
   * Sets the order value
   */
  void set_order(const int order)
  {
    order_ = order;
  }

protected:
  /**
   * Virtual function which must be implemented in children
   */
  virtual std::pair<bool,Vertex_handle>
  do_perturb(const Vertex_handle& v,
             const std::vector<Cell_handle>& slivers,
             C3T3& c3t3,
             const MeshDomain& domain,
             const SliverCriterion& criterion,
             const FT& sliver_bound,
             std::vector<Vertex_handle>& modified_vertices,
             bool *could_lock_zone = nullptr) const = 0;

  /**
   * @brief a helper function which returns the amplitude of perturbation
   */
  FT compute_perturbation_sq_amplitude(const Vertex_handle& v,
                                       const C3T3& c3t3,
                                       const FT& sq_factor) const
  {
    // We don't care if the shortest edge is inside or outside c3t3
    return   details::min_incident_edge_sq_length(v,c3t3.triangulation())
           * sq_factor;
  }

private:
  Abstract_perturbation* p_next_;
  Abstract_perturbation* p_previous_;
  // An int to have an order between Abstract_perturbation
  int order_;

#ifdef CGAL_MESH_3_PERTURBER_VERBOSE
public:
  void reset_timer() { total_time_+= 1000*time(); timer().reset(); }
  void reset_counter() { total_counter_ += counter_; counter_ = 0; }
  int counter() const { return counter_; }
  double time() const { return timer().time(); }
  int total_counter() const { return total_counter_ + counter(); }
  std::size_t total_time() const
  { return static_cast<std::size_t>(double(total_time_) + 1000*time()); }
  virtual std::string perturbation_name() const = 0;
private:
  CGAL::Timer &timer() const
  {
#ifdef CGAL_LINKED_WITH_TBB
    return timer_.local();
#else
    return timer_;
#endif
  }
  mutable int counter_;
#ifdef CGAL_LINKED_WITH_TBB
  mutable tbb::enumerable_thread_specific<CGAL::Timer> timer_;
  std::atomic<int> total_counter_;
  std::atomic<std::size_t> total_time_;
#else
  mutable CGAL::Timer timer_;
  int total_counter_;
  std::size_t total_time_;
#endif
#endif
};

template <typename C3T3, typename MD, typename SC>
inline
bool
operator<(const Abstract_perturbation<C3T3,MD,SC>& lhs,
          const Abstract_perturbation<C3T3,MD,SC>& rhs)
{
  return lhs.order() < rhs.order();
}

/**
 * @class Gradient_based_perturbation
 *
 * Base class for gradient based perturbations. The goal of these perturbations
 * is to make the sliver flip.
 */
template <typename C3T3, typename MeshDomain, typename SliverCriterion>
class Gradient_based_perturbation
  : public Abstract_perturbation<C3T3,MeshDomain,SliverCriterion>
{
protected:
  typedef Abstract_perturbation<C3T3, MeshDomain, SliverCriterion> Base;

  typedef typename Base::Vertex_handle                Vertex_handle;
  typedef typename Base::Cell_handle                  Cell_handle;

  typedef typename C3T3::Triangulation                Tr;
  typedef typename Tr::Geom_traits                    Gt;
  typedef typename Gt::FT                             FT;
  typedef typename Gt::Vector_3                       Vector_3;

  typedef typename Tr::Bare_point                     Bare_point;
  typedef typename Tr::Weighted_point                 Weighted_point;

public:
  /**
   * @brief Constructor
   */
  Gradient_based_perturbation(unsigned int max_step_nb,
                              double step_size)
    : max_step_nb_(max_step_nb)
    , sq_step_size_(step_size*step_size) { }

  /**
   * @brief destructor
   */
  virtual ~Gradient_based_perturbation() {}

protected:
  /**
   * Virtual function which must be implemented in children
   */
  virtual std::pair<bool,Vertex_handle>
  do_perturb(const Vertex_handle& v,
             const std::vector<Cell_handle>& slivers,
             C3T3& c3t3,
             const MeshDomain& domain,
             const SliverCriterion& criterion,
             const FT& sliver_bound,
             std::vector<Vertex_handle>& modified_vertices,
             bool *could_lock_zone = nullptr) const = 0;

protected:
  // -----------------------------------
  // Protected methods
  // -----------------------------------

  /**
   * Tries to apply a gradient perturbation using direction of
   * \c gradient_vector
   */
  std::pair<bool, Vertex_handle>
  apply_perturbation(const Vertex_handle& v,
                     const Vector_3& gradient_vector,
                     C3T3& c3t3,
                     const MeshDomain& domain,
                     const SliverCriterion& criterion,
                     std::vector<Vertex_handle>& modified_vertices,
                     bool *could_lock_zone = nullptr) const
  {
    typedef Triangulation_helpers<typename C3T3::Triangulation> Th;

    const Tr& tr = c3t3.triangulation();

    typename Gt::Construct_point_3 cp =
      tr.geom_traits().construct_point_3_object();
    typename Gt::Construct_weighted_point_3 cwp =
      tr.geom_traits().construct_weighted_point_3_object();
    typename Gt::Compute_squared_length_3 sq_length =
      tr.geom_traits().compute_squared_length_3_object();
    typename Gt::Construct_translated_point_3 translate =
      tr.geom_traits().construct_translated_point_3_object();
    typename Gt::Construct_vector_3 vector =
      tr.geom_traits().construct_vector_3_object();

    // create a helper
    typedef C3T3_helpers<C3T3,MeshDomain> C3T3_helpers;
    C3T3_helpers helper(c3t3, domain);

    modified_vertices.clear();

    // norm depends on the local size of the mesh
    FT sq_norm = this->compute_perturbation_sq_amplitude(v, c3t3, sq_step_size_);
    FT step_length = CGAL::sqrt(sq_norm / sq_length(gradient_vector));
    Vector_3 step_vector = step_length * gradient_vector;
    const Weighted_point& weighted_initial_loc = c3t3.triangulation().point(v);
    Bare_point initial_loc = cp(weighted_initial_loc);
    Bare_point new_loc = translate(initial_loc, step_vector);
    Bare_point final_loc = new_loc;

    if ( c3t3.in_dimension(v) < 3 )
      final_loc = helper.project_on_surface(v, new_loc);

    Vector_3 move_vector = vector(initial_loc, final_loc);

    unsigned int i = 0;
    // Concurrent-safe version
    if (could_lock_zone)
    {
      // as long as no topological change takes place
      while(Th().no_topological_change__without_set_point(c3t3.triangulation(),
                                                          v, cwp(final_loc)) &&
           ++i <= max_step_nb_ )
      {
        new_loc = translate(new_loc, step_vector);

        if ( c3t3.in_dimension(v) == 3 )
          final_loc = new_loc;
        else
          final_loc = helper.project_on_surface(v, new_loc);
      }
    }
    else
    {
      while( Th().no_topological_change(c3t3.triangulation(), v,
                                        move_vector, cwp(final_loc)) &&
             ++i <= max_step_nb_ )
      {
        new_loc = translate(new_loc, step_vector);

        if ( c3t3.in_dimension(v) == 3 )
          final_loc = new_loc;
        else
          final_loc = helper.project_on_surface(v, new_loc);

        move_vector = vector(initial_loc, final_loc);
      }
    }

    // Topology could not change moving this vertex
    if ( i > max_step_nb_ ||
         Th().inside_protecting_balls(c3t3.triangulation(), v, final_loc) )
      return std::make_pair(false, v);

    // we know that there will be a combinatorial change
    return helper.update_mesh_topo_change(v,
                                          cwp(final_loc),
                                          criterion,
                                          std::back_inserter(modified_vertices),
                                          could_lock_zone);
  }

private:
  unsigned int max_step_nb_;
  double sq_step_size_;
};

/**
 * @class Sq_radius_perturbation
 *
 * Gradient perturbation which tends to maximize cell squared radius.
 */
template <typename C3T3, typename MeshDomain, typename SliverCriterion>
class Sq_radius_perturbation
  : public Gradient_based_perturbation<C3T3,MeshDomain,SliverCriterion>
{
protected:
  typedef Gradient_based_perturbation<C3T3, MeshDomain, SliverCriterion>  Base;

  typedef typename C3T3::Triangulation                Tr;
  typedef typename C3T3::Triangulation                Triangulation;

  typedef typename Base::Vertex_handle                Vertex_handle;
  typedef typename Base::Cell_handle                  Cell_handle;

  typedef typename Tr::Geom_traits                    Gt;
  typedef typename Gt::FT                             FT;
  typedef typename Gt::Vector_3                       Vector_3;

  typedef typename Tr::Bare_point                     Bare_point;
  typedef typename Tr::Weighted_point                 Weighted_point;

public:
  /**
   * @brief Constructor
   */
  Sq_radius_perturbation(unsigned int max_step_nb,
                         double step_size)
    : Base(max_step_nb,step_size) {}

  /**
   * @brief destructor
   */
  virtual ~Sq_radius_perturbation() {}

#ifdef CGAL_MESH_3_PERTURBER_VERBOSE
  virtual std::string perturbation_name() const
  {
    return std::string("Sq radius gradient perturbation");
  }
#endif

protected:
  /**
   * do_perturb implementation
   */
  virtual std::pair<bool,Vertex_handle>
  do_perturb(const Vertex_handle& v,
             const std::vector<Cell_handle>& slivers,
             C3T3& c3t3,
             const MeshDomain& domain,
             const SliverCriterion& criterion,
             const FT&,
             std::vector<Vertex_handle>& modified_vertices,
             bool *could_lock_zone = nullptr) const
  {
    CGAL_precondition(!slivers.empty());

    Vector_3 grad_vector = compute_gradient_vector(c3t3, v, slivers);

    // Exit if grad_vector is not relevant
    if ( CGAL::NULL_VECTOR == grad_vector )
      return std::make_pair(false,v);

    return Base::apply_perturbation(v,
                                    grad_vector,
                                    c3t3,
                                    domain,
                                    criterion,
                                    modified_vertices,
                                    could_lock_zone);
  }

private:
  // -----------------------------------
  // Private methods
  // -----------------------------------

  /**
   * @brief compute the gradient vector
   */
  Vector_3
  compute_gradient_vector(const C3T3& c3t3,
                          const Vertex_handle& v,
                          const std::vector<Cell_handle>& slivers) const
  {
    switch (slivers.size())
    {
      case 1:
        return compute_gradient_vector(c3t3, slivers.front(), v);
        break;
      case 2:
      {
        Vector_3 v1 = compute_gradient_vector(c3t3, slivers.front(), v);
        Vector_3 v2 = compute_gradient_vector(c3t3, slivers.back(), v);
        if( v1 * v2 > 0 )
          // "+0.5" because sq_radius has to go up
          return 0.5*(v1 + v2);
        break;
      }
      default:
        break;
    }

    // May happen if sq_radius_gradient is not relevant for this vertex
    return CGAL::NULL_VECTOR;
  }

  /**
   * @brief compute the gradient vector
   */
  Vector_3 compute_gradient_vector(const C3T3& c3t3,
                                   const Cell_handle& cell,
                                   const Vertex_handle& v) const
  {
    const Triangulation& tr = c3t3.triangulation();

    typename Gt::Construct_point_3 cp =
      c3t3.triangulation().geom_traits().construct_point_3_object();
    typename Gt::Construct_translated_point_3 translate =
      c3t3.triangulation().geom_traits().construct_translated_point_3_object();

    unsigned int index = cell->index(v);

    const Weighted_point& wvp = tr.point(cell, index);
    const Weighted_point& wp2 = tr.point(cell, (index+1)&3);
    const Weighted_point& wp3 = tr.point(cell, (index+2)&3);
    const Weighted_point& wp4 = tr.point(cell, (index+3)&3);

    // translate the tet so that 'wp4' is the origin
    Vector_3 translate_to_origin(CGAL::ORIGIN, cp(wp4));
    const Bare_point& p1 = translate(cp(wvp), - translate_to_origin);
    const Bare_point& p2 = translate(cp(wp2), - translate_to_origin);
    const Bare_point& p3 = translate(cp(wp3), - translate_to_origin);

    // pre-compute everything
    FT sq_p1 = p1.x()*p1.x() + p1.y()*p1.y() + p1.z()*p1.z();
    FT sq_p2 = p2.x()*p2.x() + p2.y()*p2.y() + p2.z()*p2.z();
    FT sq_p3 = p3.x()*p3.x() + p3.y()*p3.y() + p3.z()*p3.z();

    // every derivative is computed w.r.t p1 (x1, y1, z1)
    FT da_dx = p2.y()*p3.z() - p3.y()*p2.z();
    FT da_dy = p2.z()*p3.x() - p2.x()*p3.z();
    FT da_dz = p2.x()*p3.y() - p3.x()*p2.y();

    FT dDx_dx = -2*p1.x()*da_dx;
    FT dDx_dy = -2*p1.y()*da_dx + sq_p2*p3.z() - sq_p3*p2.z();
    FT dDx_dz = -2*p1.z()*da_dx - sq_p2*p3.y() + sq_p3*p2.y();

    FT dDy_dx = -2*p1.x()*da_dy - sq_p2*p3.z() + sq_p3*p2.z();
    FT dDy_dy = -2*p1.y()*da_dy;
    FT dDy_dz = -2*p1.z()*da_dy + sq_p2*p3.x() - sq_p3*p2.x();

    FT dDz_dx = -2*p1.x()*da_dz + sq_p2*p3.y() - sq_p3*p2.y();
    FT dDz_dy = -2*p1.y()*da_dz - sq_p2*p3.x() + sq_p3*p2.x();
    FT dDz_dz = -2*p1.z()*da_dz;

    FT a  = p1.x()*da_dx + p1.y()*da_dy + p1.z()*da_dz;
    if ( CGAL_NTS is_zero(a) )
      return CGAL::NULL_VECTOR;

    FT Dx = -sq_p1*da_dx + p1.y()*(sq_p2*p3.z() - sq_p3*p2.z()) - p1.z()*(sq_p2*p3.y() - sq_p3*p2.y());
    FT Dy = -sq_p1*da_dy - p1.x()*(sq_p2*p3.z() - sq_p3*p2.z()) + p1.z()*(sq_p2*p3.x() - sq_p3*p2.x());
    FT Dz = -sq_p1*da_dz + p1.x()*(sq_p2*p3.y() - sq_p3*p2.y()) - p1.y()*(sq_p2*p3.x() - sq_p3*p2.x());

    // compute gradient vector
    FT sum_sqD = Dx*Dx + Dy*Dy + Dz*Dz;
    FT gx = (Dx*dDx_dx + Dy*dDy_dx + Dz*dDz_dx) / (2.0*a*a) - (da_dx * sum_sqD) / (2.0*a*a*a);
    FT gy = (Dx*dDx_dy + Dy*dDy_dy + Dz*dDz_dy) / (2.0*a*a) - (da_dy * sum_sqD) / (2.0*a*a*a);
    FT gz = (Dx*dDx_dz + Dy*dDy_dz + Dz*dDz_dz) / (2.0*a*a) - (da_dz * sum_sqD) / (2.0*a*a*a);

    return Vector_3(gx, gy, gz);
  }
};

/**
 * @class Volume_perturbation
 *
 * Gradient perturbation which tends to decrease cell volume.
 */
template <typename C3T3, typename MeshDomain, typename SliverCriterion>
class Volume_perturbation
: public Gradient_based_perturbation<C3T3,MeshDomain,SliverCriterion>
{
protected:
  typedef Gradient_based_perturbation<C3T3, MeshDomain, SliverCriterion> Base;

  typedef typename Base::Vertex_handle                Vertex_handle;
  typedef typename Base::Cell_handle                  Cell_handle;

  typedef typename C3T3::Triangulation                Tr;
  typedef typename Tr::Geom_traits                    Gt;
  typedef typename Gt::FT                             FT;
  typedef typename Gt::Vector_3                       Vector_3;

  typedef typename Tr::Bare_point                     Bare_point;
  typedef typename Tr::Weighted_point                 Weighted_point;

public:
  /**
   * @brief Constructor
   */
  Volume_perturbation(unsigned int max_step_nb,
                      double step_size)
  : Base(max_step_nb,step_size) {}

  /**
   * @brief Destructor
   */
  virtual ~Volume_perturbation() {}

#ifdef CGAL_MESH_3_PERTURBER_VERBOSE
  virtual std::string perturbation_name() const
  {
    return std::string("Volume gradient perturbation");
  }
#endif

protected:
  /**
   * @brief do_perturb implementation
   */
  virtual std::pair<bool,Vertex_handle>
  do_perturb(const Vertex_handle& v,
             const std::vector<Cell_handle>& slivers,
             C3T3& c3t3,
             const MeshDomain& domain,
             const SliverCriterion& criterion,
             const FT&,
             std::vector<Vertex_handle>& modified_vertices,
             bool *could_lock_zone = nullptr) const
  {
    CGAL_precondition(!slivers.empty());

    Vector_3 grad_vector = compute_gradient_vector(c3t3, v, slivers);

    // Exit if grad_vector is not relevant
    if ( CGAL::NULL_VECTOR == grad_vector )
      return std::make_pair(false,v);

    return Base::apply_perturbation(v,
                                    grad_vector,
                                    c3t3,
                                    domain,
                                    criterion,
                                    modified_vertices,
                                    could_lock_zone);
  }

private:
  // -----------------------------------
  // Private methods
  // -----------------------------------

  /**
   * @brief compute the gradient vector
   */
  Vector_3
  compute_gradient_vector(const C3T3& c3t3,
                          const Vertex_handle& v,
                          const std::vector<Cell_handle>& slivers) const
  {
    switch (slivers.size())
    {
      case 1:
        return -1*compute_gradient_vector(c3t3, slivers.front(), v);
        break;
      case 2:
      {
        Vector_3 v1 = compute_gradient_vector(c3t3, slivers.front(), v);
        Vector_3 v2 = compute_gradient_vector(c3t3, slivers.back(), v);
        if( v1 * v2 > 0 )
          // "-0.5" because volume has to go down
          return -0.5 * (v1 + v2);
        break;
      }
      default:
        break;
    }

    // May happen if volume_gradient is not relevant for this vertex
    return CGAL::NULL_VECTOR;
  }

  /**
   * @brief compute the gradient vector
   */
  Vector_3 compute_gradient_vector(const C3T3& c3t3,
                                   const Cell_handle& cell,
                                   const Vertex_handle& v) const
  {
    CGAL_precondition(cell->has_vertex(v));

    const typename C3T3::Triangulation& tr = c3t3.triangulation();
    typename Gt::Construct_point_3 cp = tr.geom_traits().construct_point_3_object();

    const int i = cell->index(v);

    // fixed vertices: (the ones with index != i)
    int k1 = (i+1)&3;
    int k2 = (i+2)&3;
    int k3 = (i+3)&3;

    if ( (i&1) == 0 )
      std::swap(k1,k3);

    const Weighted_point& wp1 = tr.point(cell, k1);
    const Weighted_point& wp2 = tr.point(cell, k2);
    const Weighted_point& wp3 = tr.point(cell, k3);
    const Bare_point& p1 = cp(wp1);
    const Bare_point& p2 = cp(wp2);
    const Bare_point& p3 = cp(wp3);

    FT gx =  p2.y()*p3.z() + p1.y()*(p2.z()-p3.z())
            - p3.y()*p2.z() - p1.z()*(p2.y()-p3.y());

    FT gy = -p2.x()*p3.z() - p1.x()*(p2.z()-p3.z())
            + p3.x()*p2.z() + p1.z()*(p2.x()-p3.x());

    FT gz =  p2.x()*p3.y() + p1.x()*(p2.y()-p3.y())
            - p3.x()*p2.y() - p1.y()*(p2.x()-p3.x());

    return (1.0 / 6.0 * Vector_3(gx, gy, gz));
  }
};

/**
 * @class Dihedral_angle_perturbation
 *
 * Gradient perturbation which tends to decrease cell minimal dihedral angle.
 */
template <typename C3T3, typename MeshDomain, typename SliverCriterion>
class Dihedral_angle_perturbation
  : public Gradient_based_perturbation<C3T3,MeshDomain,SliverCriterion>
{
protected:
  typedef Gradient_based_perturbation<C3T3, MeshDomain, SliverCriterion> Base;

  typedef typename Base::Vertex_handle                Vertex_handle;
  typedef typename Base::Cell_handle                  Cell_handle;

  typedef typename C3T3::Triangulation                Tr;
  typedef typename Tr::Geom_traits                    Gt;
  typedef typename Gt::FT                             FT;
  typedef typename Gt::Vector_3                       Vector_3;

  typedef typename Tr::Bare_point                     Bare_point;
  typedef typename Tr::Weighted_point                 Weighted_point;

public:
  /**
   * @brief constructor
   */
  Dihedral_angle_perturbation(unsigned int max_step_nb,
                                  double step_size)
    : Base(max_step_nb,step_size) {}

  /**
   * @brief destructor
   */
  virtual ~Dihedral_angle_perturbation() {}

#ifdef CGAL_MESH_3_PERTURBER_VERBOSE
  virtual std::string perturbation_name() const
  {
    return std::string("Dihedral angle gradient perturbation");
  }
#endif

protected:
  /**
   * @brief do_perturb implementation
   */
  virtual std::pair<bool,Vertex_handle>
  do_perturb(const Vertex_handle& v,
             const std::vector<Cell_handle>& slivers,
             C3T3& c3t3,
             const MeshDomain& domain,
             const SliverCriterion& criterion,
             const FT&,
             std::vector<Vertex_handle>& modified_vertices,
             bool *could_lock_zone = nullptr) const
  {
    CGAL_precondition(!slivers.empty());

    Vector_3 grad_vector = compute_gradient_vector(c3t3, v, slivers);

    // Exit if grad_vector is not relevant
    if ( CGAL::NULL_VECTOR == grad_vector )
      return std::make_pair(false,v);

    return Base::apply_perturbation(v,
                                    grad_vector,
                                    c3t3,
                                    domain,
                                    criterion,
                                    modified_vertices,
                                    could_lock_zone);
  }

private:
  // -----------------------------------
  // Private methods
  // -----------------------------------

  /**
   * @brief compute the gradient vector
   */
  Vector_3
  compute_gradient_vector(const C3T3& c3t3,
                          const Vertex_handle& v,
                          const std::vector<Cell_handle>& slivers) const
  {
    switch (slivers.size())
    {
      case 1:
        return -1*compute_gradient_vector(c3t3, slivers.front(), v);
        break;
      case 2:
      {
        Vector_3 v1 = compute_gradient_vector(c3t3, slivers.front(), v);
        Vector_3 v2 = compute_gradient_vector(c3t3, slivers.back(), v);
        if( v1 * v2 > 0 )
          // "-0.5" because angle has to go down
          return -0.5*(v1 + v2);
        break;
      }
      default:
        break;
    }

    // May happen if dihedral gradient is not relevant for this vertex
    return CGAL::NULL_VECTOR;
  }

  /**
   * @brief compute the gradient vector
   */
  Vector_3
  compute_gradient_vector(const C3T3& c3t3,
                          const Cell_handle& cell,
                          const Vertex_handle& v) const
  {
    CGAL_assertion(cell->has_vertex(v));
    const typename C3T3::Triangulation& tr = c3t3.triangulation();

    typename Gt::Construct_point_3 cp =
      tr.geom_traits().construct_point_3_object();
    typename Gt::Compute_squared_distance_3 sq_distance =
      tr.geom_traits().compute_squared_distance_3_object();

    const int i = cell->index(v);
    const Weighted_point& wp0 = tr.point(cell, i);
    const Bare_point& p0 = cp(wp0);

    // Other indices
    int k1 = (i+1)&3;
    int k2 = (i+2)&3;
    int k3 = (i+3)&3;

    FT angle_k1k2 = abs_dihedral_angle(c3t3, p0, k1, k2, k3, cell);
    FT angle_k2k3 = abs_dihedral_angle(c3t3, p0, k2, k3, k1, cell);
    FT angle_k3k1 = abs_dihedral_angle(c3t3, p0, k3, k1, k2, cell);

    if ( angle_k1k2 > angle_k2k3 )
      std::swap(k1,k3);

    // Here we know that min_angle_k1k2 <= min_angle_k2k3
    if ( angle_k1k2 > angle_k3k1 )
      std::swap(k2,k3);

    // Here edge k1k2 minimizes dihedral angle
    const Weighted_point& wp1 = tr.point(cell, k1);
    const Weighted_point& wp2 = tr.point(cell, k2);
    const Weighted_point& wp3 = tr.point(cell, k3);
    const Bare_point& p1 = cp(wp1);
    const Bare_point& p2 = cp(wp2);
    const Bare_point& p3 = cp(wp3);

    // grad of min dihedral angle (in cell) wrt p0
    const Vector_3 p1p0 (p1,p0);
    const Vector_3 p1p2 (p1,p2);
    const Vector_3 p1p3 (p1,p3);

    FT a_02 = CGAL::abs(details::angle_in_radian(p1p0, p1p2));
    FT a_03 = CGAL::abs(details::angle_in_radian(p1p0, p1p3));

    Vector_3 n0 = normal_estimate(c3t3, cell, k3);
    Vector_3 n1 = normal_estimate(c3t3, cell, k2);

    const FT d_p0p1 = CGAL::sqrt(sq_distance(p0,p1));
    CGAL_assertion(!is_zero(d_p0p1));

    return ( (-1./ d_p0p1) * (cotangent(a_02)*n0 + cotangent(a_03)*n1) );
  }

  /**
   * @brief returns the absolute value of dihedral angle of p,p(k1),p(k2),p(k3)
   */
  FT abs_dihedral_angle(const C3T3& c3t3,
                        const Bare_point& p,
                        const int k1,
                        const int k2,
                        const int k3,
                        const Cell_handle& cell) const
  {
    const typename C3T3::Triangulation& tr = c3t3.triangulation();

    typename Gt::Construct_point_3 cp =
      tr.geom_traits().construct_point_3_object();
    typename Gt::Compute_approximate_dihedral_angle_3 approx_dihedral_angle =
      tr.geom_traits().compute_approximate_dihedral_angle_3_object();

    const Weighted_point& wp1 = tr.point(cell, k1);
    const Weighted_point& wp2 = tr.point(cell, k2);
    const Weighted_point& wp3 = tr.point(cell, k3);

    return CGAL::abs(approx_dihedral_angle(cp(p), cp(wp1), cp(wp2), cp(wp3)));
  }

  /**
   * @brief returns the cotangent of \c value
   */
  FT cotangent(const FT& value) const
  {
    return FT(1/std::tan(CGAL::to_double(value)));
  }

  /**
   * @brief returns the normal of facet (ch,i), oriented from inside to outside
   * of \c ch
   */
  Vector_3 normal_estimate(const C3T3& c3t3, const Cell_handle& ch, const int i) const
  {
    const typename C3T3::Triangulation& tr = c3t3.triangulation();

    typename Gt::Construct_point_3 cp = tr.geom_traits().construct_point_3_object();

    int k1 = (i+1)&3;
    int k2 = (i+2)&3;
    int k3 = (i+3)&3;

    // Orient normals to the outside of cell
    if ( (i&1) == 1 )
      std::swap(k1,k2);

    const Weighted_point& wp1 = tr.point(ch, k1);
    const Weighted_point& wp2 = tr.point(ch, k2);
    const Weighted_point& wp3 = tr.point(ch, k3);
    const Bare_point& p1 = cp(wp1);
    const Bare_point& p2 = cp(wp2);
    const Bare_point& p3 = cp(wp3);

    // compute normal and return it
    return normal(p1, p2, p3);
  }
};

/**
 * @class Gradient_based_perturbation
 *
 * Base class for random based perturbations. The goal of these perturbations
 * is to improve the mesh quality (using a flip or not).
 */
template <typename C3T3, typename MeshDomain, typename SliverCriterion>
class Random_based_perturbation
  : public Abstract_perturbation<C3T3,MeshDomain,SliverCriterion>
{
protected:
  typedef Abstract_perturbation<C3T3, MeshDomain, SliverCriterion> Base;

  typedef typename Base::Vertex_handle                Vertex_handle;
  typedef typename Base::Cell_handle                  Cell_handle;

  typedef typename C3T3::Triangulation                Tr;
  typedef typename Tr::Geom_traits                    Gt;
  typedef typename Gt::FT                             FT;
  typedef typename Gt::Vector_3                       Vector_3;

  typedef typename Tr::Bare_point                     Bare_point;

public:
  /**
   * @brief constructor
   */
  Random_based_perturbation(unsigned int max_try_nb,
                            const FT& sphere_radius)
    : max_try_nb_(max_try_nb)
    , sphere_radius_(sphere_radius)
    , sphere_sq_radius_(sphere_radius*sphere_radius)
    , generator_() // initialize the random generator deterministically
    , uni_dist_(0,1)
    , random_(generator_, uni_dist_) {}

  /**
   * @brief destructor
   */
  virtual ~Random_based_perturbation() {}

protected:
  /**
   * Virtual function which must be implemented in children
   */
  virtual std::pair<bool,Vertex_handle>
  do_perturb(const Vertex_handle& v,
             const std::vector<Cell_handle>& slivers,
             C3T3& c3t3,
             const MeshDomain& domain,
             const SliverCriterion& criterion,
             const FT& sliver_bound,
             std::vector<Vertex_handle>& modified_vertices,
             bool *could_lock_zone = nullptr) const = 0;

protected:
  // -----------------------------------
  // Protected methods (children helpers)
  // -----------------------------------

  /// Access functions
  unsigned int max_try_nb() const { return max_try_nb_; }
  FT sphere_radius() const { return sphere_radius_; }
  FT sphere_sq_radius() const { return sphere_sq_radius_; }

  /**
   * @brief returns a FT between \c min and \c max
   */
  FT random_ft(const FT& min = FT(0.), const FT& max = FT(1.)) const
  {
    FT r (random_());
    return r*(max-min) + min;
  }

  /**
   * @brief returns a random vector with size \c vector_size
   */
  Vector_3 random_vector_fixed_size(const C3T3& c3t3,
                                    const FT& vector_sq_size) const
  {
    typename Gt::Compute_squared_length_3 sq_length =
      c3t3.triangulation().geom_traits().compute_squared_length_3_object();

    Vector_3 rnd_vector(random_ft(),random_ft(),random_ft());
    FT rnd_sq_size = sq_length(rnd_vector);

    if ( ! CGAL_NTS is_zero(rnd_sq_size) )
      return CGAL::sqrt(vector_sq_size / rnd_sq_size) * rnd_vector;
    else
      return CGAL::NULL_VECTOR;
  }

  /**
   * @brief returns a random vector with size between 0 and \c vector_size
   */
  Vector_3 random_vector_max_size(const C3T3& c3t3,
                                  const FT& vector_max_sq_size) const
  {
    return random_ft() * random_vector_fixed_size(c3t3, vector_max_sq_size);
  }

  /**
   * @brief returns a random vector.
   */
  Vector_3 random_vector(const C3T3& c3t3,
                         const FT& vector_sq_size,
                         bool fixed_size) const
  {
    if ( fixed_size )
      return random_vector_fixed_size(c3t3, vector_sq_size);
    else
      return random_vector_max_size(c3t3, vector_sq_size);
  }

private:
  /// The maximum random point which will be tried
  unsigned int max_try_nb_;
  // The radius of the sphere which will contain random points
  FT sphere_radius_;
  FT sphere_sq_radius_;

  // boost random generator
  typedef boost::lagged_fibonacci607 base_generator_type;
  base_generator_type generator_;
  boost::uniform_real<FT> uni_dist_;
  mutable boost::variate_generator<base_generator_type&,
                                   boost::uniform_real<FT> > random_;
};

/**
 * @class Li_random_perturbation
 *
 * Li random based perturbation.
 */
template <typename C3T3, typename MeshDomain, typename SliverCriterion>
class Li_random_perturbation
: public Random_based_perturbation<C3T3,MeshDomain,SliverCriterion>
{
protected:
  typedef Random_based_perturbation<C3T3, MeshDomain, SliverCriterion> Base;

  typedef typename Base::Vertex_handle                Vertex_handle;
  typedef typename Base::Cell_handle                  Cell_handle;

  typedef typename C3T3::Triangulation                Tr;
  typedef typename Tr::Geom_traits                    Gt;
  typedef typename Gt::FT                             FT;
  typedef typename Gt::Vector_3                       Vector_3;

  typedef typename Tr::Bare_point                     Bare_point;
  typedef typename Tr::Weighted_point                 Weighted_point;

  typedef CGAL::Hash_handles_with_or_without_timestamps Hash_fct;
  typedef boost::unordered_set<Vertex_handle, Hash_fct> Vertex_set;

public:
  /**
   * @brief constructor
   */
  Li_random_perturbation(unsigned int max_try_nb,
                         const FT& sphere_radius,
                         const bool on_sphere = false)
    : Base(max_try_nb,sphere_radius)
    , on_sphere_(on_sphere) {}

  /**
   * @brief destructor
   */
  virtual ~Li_random_perturbation() {}

#ifdef CGAL_MESH_3_PERTURBER_VERBOSE
  virtual std::string perturbation_name() const
  {
    std::stringstream name;
    name << "Li random perturbation [" << this->max_try_nb() << ", "
         << this->sphere_radius() << ", " << (on_sphere_?"on":"in") << "]";

    return name.str();
  }
#endif

protected:
  /**
   * @brief do_perturb implementation
   */
  virtual std::pair<bool,Vertex_handle>
  do_perturb(const Vertex_handle& v,
             const std::vector<Cell_handle>& slivers,
             C3T3& c3t3,
             const MeshDomain& domain,
             const SliverCriterion& criterion,
             const FT& sliver_bound,
             std::vector<Vertex_handle>& modified_vertices,
             bool *could_lock_zone = nullptr) const
  {
    CGAL_precondition(!slivers.empty());

    return apply_perturbation(v, slivers, c3t3, domain, criterion,
                              sliver_bound, modified_vertices,
                              could_lock_zone);
  }

private:
  // -----------------------------------
  // Private methods
  // -----------------------------------

  /**
   * @brief try to improve mesh using random moves of \c v
   */
  std::pair<bool,Vertex_handle>
  apply_perturbation(const Vertex_handle& v,
                     const std::vector<Cell_handle>& slivers,
                     C3T3& c3t3,
                     const MeshDomain& domain,
                     const SliverCriterion& criterion,
                     const FT& sliver_bound,
                     std::vector<Vertex_handle>& modified_vertices,
                     bool *could_lock_zone = nullptr) const
  {
    typedef Triangulation_helpers<typename C3T3::Triangulation> Th;

    typename Gt::Construct_point_3 cp =
      c3t3.triangulation().geom_traits().construct_point_3_object();
    typename Gt::Construct_weighted_point_3 cwp =
      c3t3.triangulation().geom_traits().construct_weighted_point_3_object();
    typename Gt::Construct_translated_point_3 translate =
      c3t3.triangulation().geom_traits().construct_translated_point_3_object();
    typename Gt::Construct_vector_3 vector =
      c3t3.triangulation().geom_traits().construct_vector_3_object();

    modified_vertices.clear();

    // Create a helper
    typedef C3T3_helpers<C3T3,MeshDomain> C3T3_helpers;
    C3T3_helpers helper(c3t3, domain);

    // norm depends on the local size of the mesh
    FT sq_norm = this->compute_perturbation_sq_amplitude(v, c3t3, this->sphere_sq_radius());
    const Weighted_point& weighted_initial_location = c3t3.triangulation().point(v);
    const Bare_point initial_location = cp(weighted_initial_location);

    // Initialize loop variables
    bool criterion_improved = false;
    Vertex_handle moving_vertex = v;
    Vertex_set mod_vertices;

    // TODO: use no_topo_change to compute new_angle without moving everything
    unsigned int try_nb = 0;
    while ( ++try_nb <= Base::max_try_nb() )
    {
      Vector_3 delta = this->random_vector(c3t3, sq_norm, on_sphere_);

      // always from initial_location!
      Bare_point new_location = translate(initial_location, delta);

      if ( c3t3.in_dimension(moving_vertex) < 3 )
        new_location = helper.project_on_surface(moving_vertex, new_location);

      // check that we don't insert a vertex inside a protecting ball
      if(Th().inside_protecting_balls(c3t3.triangulation(),
                                      moving_vertex, new_location))
        continue;

      // try to move vertex
      std::vector<Vertex_handle> tmp_mod_vertices;
      std::pair<bool,Vertex_handle> update =
        helper.update_mesh(moving_vertex,
                           vector(initial_location, new_location),
                           cwp(new_location),
                           criterion,
                           std::back_inserter(tmp_mod_vertices),
                           could_lock_zone);

      if (could_lock_zone && *could_lock_zone == false)
        return std::make_pair(false, Vertex_handle());

      // get new vertex
      moving_vertex = update.second;

      if ( update.first )
      {
        criterion_improved = true;

        mod_vertices.insert(tmp_mod_vertices.begin(), tmp_mod_vertices.end());

        std::size_t nois;
#ifdef CGAL_NEW_INCIDENT_SLIVERS
        nois  = helper.number_of_incident_slivers(moving_vertex, criterion, sliver_bound);
#else
        std::vector<Cell_handle> new_slivers =
          helper.incident_slivers(moving_vertex, criterion, sliver_bound);
        nois = new_slivers.size();
#endif

        // If sliver number has decreased, we won
        if(nois < slivers.size() )
        {
          std::copy(mod_vertices.begin(),
                    mod_vertices.end(),
                    std::back_inserter(modified_vertices));

          return std::make_pair(true,moving_vertex);
        }
      }
    }//end while ( ++try_nb <= Base::max_try_nb() )

    if ( criterion_improved )
    {
      std::copy(mod_vertices.begin(),
                mod_vertices.end(),
                std::back_inserter(modified_vertices));
    }

    // Moving vertex is located on the best location
    return std::make_pair(criterion_improved, moving_vertex);
  }

private:
  // If set to \c true, then random point will be generated on sphere surface.
  bool on_sphere_;
};

} // end namespace Mesh_3

} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_MESH_3_VERTEX_PERTURBATION_H
