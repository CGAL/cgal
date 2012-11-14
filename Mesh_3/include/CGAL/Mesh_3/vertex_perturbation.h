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

#ifndef CGAL_MESH_3_VERTEX_PERTURBATION_H
#define CGAL_MESH_3_VERTEX_PERTURBATION_H


#include <CGAL/Mesh_3/C3T3_helpers.h>
#include <CGAL/Mesh_3/Triangulation_helpers.h>

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
      
      typename K::Compute_scalar_product_3 scalar_product =
        k.compute_scalar_product_3_object();
      typename K::Construct_cross_product_vector_3 cross_product =
        k.construct_cross_product_vector_3_object();
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
    edge_sq_length(const typename Tr::Edge& e)
    {
      typedef typename Tr::Geom_traits Gt;
      typedef typename Gt::Point_3 Point_3;
      
      typename Gt::Compute_squared_distance_3 sq_distance 
        = Gt().compute_squared_distance_3_object();
      
      const Point_3& p = e.first->vertex(e.second)->point();
      const Point_3& q = e.first->vertex(e.third)->point();
      
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
      FT min_sq_length = edge_sq_length<Tr>(*eit++);
      
      for ( ; eit != edges.end() ; ++eit )
      {
        min_sq_length = (std::min)(min_sq_length, edge_sq_length<Tr>(*eit));
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
  typedef typename C3T3::Vertex_handle Vertex_handle;
  typedef typename C3T3::Cell_handle Cell_handle;
  typedef typename C3T3::Triangulation::Geom_traits::FT FT;
  
public:
  /**
   * @brief constructor
   */
  Abstract_perturbation() 
    : p_next_(NULL)
    , p_previous_(NULL)
    , order_(0)
#ifdef CGAL_MESH_3_PERTURBER_VERBOSE
    , counter_(0)
    , timer_()
    , total_counter_(0)
    , total_time_(0.)
#endif
  {}
  
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
             std::vector<Vertex_handle>& modified_vertices) const
  {
#ifndef CGAL_MESH_3_PERTURBER_VERBOSE
    return do_perturb(v, slivers, c3t3, domain, criterion,
                      sliver_bound, modified_vertices);
#else
    timer_.start();
    
    // Virtual call
    std::pair<bool,Vertex_handle> perturb =
      do_perturb(v, slivers, c3t3, domain, criterion,
                 sliver_bound, modified_vertices);
    
    if ( perturb.first )
      ++counter_;
    
    timer_.stop();
    
    return perturb;
#endif
  }
  
  /**
   * @brief Sets next perturbation 
   */
  void set_next(Abstract_perturbation* next)
  {
    p_next_ = next;
    
    if ( NULL != next )
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
             std::vector<Vertex_handle>& modified_vertices) const = 0;
  
  
  /**
   * @brief an helper function which returns the amplitude of perturbation
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
  void reset_timer() { total_time_+= time(); timer_.reset(); }
  void reset_counter() { total_counter_ += counter_; counter_ = 0; }
  int counter() const { return counter_; }
  double time() const { return timer_.time(); }
  int total_counter() const { return total_counter_ + counter(); }
  double total_time() const { return total_time_ + time(); }
  virtual std::string perturbation_name() const = 0;
private:
  mutable int counter_;
  mutable CGAL::Timer timer_;
  int total_counter_;
  double total_time_;
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
  typedef Abstract_perturbation<C3T3,MeshDomain,SliverCriterion> Base;
  typedef typename Base::Vertex_handle Vertex_handle;
  typedef typename Base::Cell_handle Cell_handle;
  typedef typename Base::FT FT;
  typedef typename C3T3::Triangulation::Geom_traits::Vector_3 Vector_3;
  
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
             std::vector<Vertex_handle>& modified_vertices) const = 0;
  
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
                     std::vector<Vertex_handle>& modified_vertices) const
  {
    typedef typename C3T3::Triangulation::Geom_traits Gt;
    typedef typename Gt::FT       FT;
    typedef typename Gt::Point_3  Point_3;
    typedef Triangulation_helpers<typename C3T3::Triangulation> Th;
    
    typename Gt::Compute_squared_length_3 sq_length =
      Gt().compute_squared_length_3_object();
    
    // create a helper
    typedef C3T3_helpers<C3T3,MeshDomain> C3T3_helpers;
    C3T3_helpers helper(c3t3, domain);
    
    modified_vertices.clear();
    
    // norm depends on the local size of the mesh
    FT sq_norm = this->compute_perturbation_sq_amplitude(v, c3t3, sq_step_size_);
    FT step_length = CGAL::sqrt(sq_norm/sq_length(gradient_vector));
    Point_3 new_loc = v->point() + step_length * gradient_vector;
    
    Point_3 final_loc = new_loc;
    if ( c3t3.in_dimension(v) < 3 )
      final_loc = helper.project_on_surface(new_loc, v);

    // as long as no topological change takes place
    unsigned int i = 0;
    while( Th().no_topological_change(c3t3.triangulation(), v, final_loc) 
          && (++i <= max_step_nb_) )
    {
      new_loc = new_loc + step_length * gradient_vector;
      
      if ( c3t3.in_dimension(v) == 3 )
        final_loc = new_loc;
      else 
        final_loc = helper.project_on_surface(new_loc, v);
    }
    
    // Topology could not change moving this vertex
    if ( i > max_step_nb_ )
      return std::make_pair(false,v);
    
    // we know that there will be a combinatorial change
    return helper.update_mesh(final_loc,
                              v,
                              criterion,
                              std::back_inserter(modified_vertices));
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
  typedef Gradient_based_perturbation<C3T3,MeshDomain,SliverCriterion> Base;
  typedef typename Base::Vertex_handle Vertex_handle;
  typedef typename Base::Cell_handle Cell_handle;
  typedef typename Base::FT FT;
  
  typedef typename C3T3::Triangulation::Geom_traits Gt;
  typedef typename Gt::Vector_3 Vector_3;
  typedef typename Gt::Point_3 Point_3;
  
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
             std::vector<Vertex_handle>& modified_vertices) const
  {
    CGAL_precondition(!slivers.empty());
    
    Vector_3 grad_vector = compute_gradient_vector(v,slivers);
    
    // Exit if grad_vector is not relevant
    if ( CGAL::NULL_VECTOR == grad_vector )
      return std::make_pair(false,v);
    
    return Base::apply_perturbation(v,
                                    grad_vector,
                                    c3t3,
                                    domain,
                                    criterion,
                                    modified_vertices);
  }
  
private:
  // -----------------------------------
  // Private methods
  // -----------------------------------
  
  /**
   * @brief compute the gradient vector
   */
  Vector_3
  compute_gradient_vector(const Vertex_handle& v,
                          const std::vector<Cell_handle>& slivers) const
  {
    switch (slivers.size())
    {
      case 1:
        return compute_gradient_vector(slivers.front(),v);
        break;
      case 2:
      {
        Vector_3 v1 = compute_gradient_vector(slivers.front(), v);
        Vector_3 v2 = compute_gradient_vector(slivers.back(), v);
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
  Vector_3 compute_gradient_vector(const Cell_handle& cell,
                                   const Vertex_handle& v) const
  {
    // translate the tet so that cell->vertex((i+3)&3) is 0_{R^3}
    unsigned int index = cell->index(v);
    Vector_3 translate_to_origin(CGAL::ORIGIN, cell->vertex((index+3)&3)->point()); //p4
    const Point_3& p1 = v->point() - translate_to_origin;
    const Point_3& p2 = cell->vertex((index+1)&3)->point() - translate_to_origin;
    const Point_3& p3 = cell->vertex((index+2)&3)->point() - translate_to_origin;
    
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
  typedef Gradient_based_perturbation<C3T3,MeshDomain,SliverCriterion> Base;
  typedef typename Base::Vertex_handle Vertex_handle;
  typedef typename Base::Cell_handle Cell_handle;
  typedef typename Base::FT FT;
  
  typedef typename C3T3::Triangulation::Geom_traits Gt;
  typedef typename Gt::Vector_3 Vector_3;
  typedef typename Gt::Point_3 Point_3;
  
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
             std::vector<Vertex_handle>& modified_vertices) const
  {
    CGAL_precondition(!slivers.empty());
    
    Vector_3 grad_vector = compute_gradient_vector(v,slivers);
    
    // Exit if grad_vector is not relevant
    if ( CGAL::NULL_VECTOR == grad_vector )
      return std::make_pair(false,v);
    
    return Base::apply_perturbation(v,
                                    grad_vector,
                                    c3t3,
                                    domain,
                                    criterion,
                                    modified_vertices);
  }
  
private:
  // -----------------------------------
  // Private methods
  // -----------------------------------
  
  /**
   * @brief compute the gradient vector
   */
  Vector_3
  compute_gradient_vector(const Vertex_handle& v,
                          const std::vector<Cell_handle>& slivers) const
  {
    switch (slivers.size())
    {
      case 1:
        return -1*compute_gradient_vector(slivers.front(),v);
        break;
      case 2:
      {
        Vector_3 v1 = compute_gradient_vector(slivers.front(), v);
        Vector_3 v2 = compute_gradient_vector(slivers.back(), v);
        if( v1 * v2 > 0 )
          // "-0.5" because volume has to go down
          return -0.5*(v1 + v2);
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
  Vector_3 compute_gradient_vector(const Cell_handle& cell,
                                   const Vertex_handle& v) const
  {
    CGAL_assertion(cell->has_vertex(v));
    
    const int i = cell->index(v);
    
    // fixed vertices: (the ones with index != i)
    int k1 = (i+1)&3;
    int k2 = (i+2)&3;
    int k3 = (i+3)&3;
    
    if ( (i&1) == 0 )
      std::swap(k1,k3);
    
    const Point_3& p1 = cell->vertex(k1)->point();
    const Point_3& p2 = cell->vertex(k2)->point();
    const Point_3& p3 = cell->vertex(k3)->point();
    
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
  typedef Gradient_based_perturbation<C3T3,MeshDomain,SliverCriterion> Base;
  typedef typename Base::Vertex_handle Vertex_handle;
  typedef typename Base::Cell_handle Cell_handle;
  typedef typename Base::FT FT;
  
  typedef typename C3T3::Triangulation::Geom_traits Gt;
  typedef typename Gt::Vector_3 Vector_3;
  typedef typename Gt::Point_3 Point_3;
  
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
             std::vector<Vertex_handle>& modified_vertices) const
  {
    CGAL_precondition(!slivers.empty());
    
    Vector_3 grad_vector = compute_gradient_vector(v,slivers);

    // Exit if grad_vector is not relevant
    if ( CGAL::NULL_VECTOR == grad_vector )
      return std::make_pair(false,v);
    
    return Base::apply_perturbation(v,
                                    grad_vector,
                                    c3t3,
                                    domain,
                                    criterion,
                                    modified_vertices);
  }
  
private:
  // -----------------------------------
  // Private methods
  // -----------------------------------
  
  /**
   * @brief compute the gradient vector
   */
  Vector_3
  compute_gradient_vector(const Vertex_handle& v,
                          const std::vector<Cell_handle>& slivers) const
  {
    switch (slivers.size())
    {
      case 1:
        return -1*compute_gradient_vector(slivers.front(),v);
        break;
      case 2:
      {
        Vector_3 v1 = compute_gradient_vector(slivers.front(), v);
        Vector_3 v2 = compute_gradient_vector(slivers.back(), v);
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
  compute_gradient_vector(const Cell_handle& cell,
                          const Vertex_handle& v) const
  {
    CGAL_assertion(cell->has_vertex(v));
    
    const int i = cell->index(v);
    const Point_3& p0 = v->point();
    
    // Other indices
    int k1 = (i+1)&3;
    int k2 = (i+2)&3;
    int k3 = (i+3)&3;
    
    FT angle_k1k2 = abs_dihedral_angle(p0,k1,k2,k3,cell);
    FT angle_k2k3 = abs_dihedral_angle(p0,k2,k3,k1,cell);
    FT angle_k3k1 = abs_dihedral_angle(p0,k3,k1,k2,cell);
    
    if ( angle_k1k2 > angle_k2k3 )
      std::swap(k1,k3);
    
    // Here we know that min_angle_k1k2 <= min_angle_k2k3
    if ( angle_k1k2 > angle_k3k1 )
      std::swap(k2,k3);

    // Here edge k1k2 minimizes dihedral angle
    const Point_3& p1 = cell->vertex(k1)->point();
    const Point_3& p2 = cell->vertex(k2)->point();
    const Point_3& p3 = cell->vertex(k3)->point();

    // grad of min dihedral angle (in cell) wrt p0
    const Vector_3 p1p0 (p1,p0);
    const Vector_3 p1p2 (p1,p2);
    const Vector_3 p1p3 (p1,p3);
    
    FT a_02 = CGAL::abs(details::angle_in_radian(p1p0, p1p2));
    FT a_03 = CGAL::abs(details::angle_in_radian(p1p0, p1p3));
    
    Vector_3 n0 = normal_estimate(cell,k3);
    Vector_3 n1 = normal_estimate(cell,k2);
    
    typename Gt::Compute_squared_distance_3 sq_distance
      = Gt().compute_squared_distance_3_object();
    
    const FT d_p0p1 = CGAL::sqrt(sq_distance(p0,p1));
    CGAL_assertion(!is_zero(d_p0p1));
    
    return ( (-1./ d_p0p1) * (cotangent(a_02)*n0 + cotangent(a_03)*n1) );
  }
  
  /**
   * @brief returns the absolute value of dihedral angle of p,p(k1),p(k2),p(k3)
   */
  FT abs_dihedral_angle(const Point_3& p,
                        const int k1,
                        const int k2,
                        const int k3,
                        const Cell_handle& cell) const
  { 
    return CGAL::abs(dihedral_angle(p,
                                    cell->vertex(k1)->point(),
                                    cell->vertex(k2)->point(),
                                    cell->vertex(k3)->point()));
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
  Vector_3 normal_estimate(const Cell_handle& ch, const int i) const
  {
    int k1 = (i+1)&3;
    int k2 = (i+2)&3;
    int k3 = (i+3)&3;
    
    // Orient normals to the outside of cell
    if ( (i&1) == 1 )
      std::swap(k1,k2);
    
    const Point_3& p1 = ch->vertex(k1)->point();
    const Point_3& p2 = ch->vertex(k2)->point();
    const Point_3& p3 = ch->vertex(k3)->point();
    
    // compute normal and return it
    return normal(p1,p2,p3);
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
  typedef Abstract_perturbation<C3T3,MeshDomain,SliverCriterion> Base;
  typedef typename Base::Vertex_handle Vertex_handle;
  typedef typename Base::Cell_handle Cell_handle;
  typedef typename Base::FT FT;
  typedef typename C3T3::Triangulation::Geom_traits::Vector_3 Vector_3;
  
public:
  /**
   * @brief constructor
   */
  Random_based_perturbation(unsigned int max_try_nb,
                            const FT& sphere_radius)
    : max_try_nb_(max_try_nb) 
    , sphere_radius_(sphere_radius)
    , sphere_sq_radius_(sphere_radius*sphere_radius)
    , generator_(boost::uint32_t(std::time(0)))
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
             std::vector<Vertex_handle>& modified_vertices) const = 0;
  
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
  Vector_3 random_vector_fixed_size(const FT& vector_sq_size) const
  {
    typedef typename C3T3::Triangulation::Geom_traits Gt;
    
    typename Gt::Compute_squared_length_3 sq_length =
      Gt().compute_squared_length_3_object();
    
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
  Vector_3 random_vector_max_size(const FT& vector_max_sq_size) const
  {
    return random_ft() * random_vector_fixed_size(vector_max_sq_size);
  }
  
  /**
   * @brief returns a random vector.
   */
  Vector_3 random_vector(const FT& vector_sq_size, bool fixed_size) const
  {
    if ( fixed_size )
      return random_vector_fixed_size(vector_sq_size);
    else
      return random_vector_max_size(vector_sq_size);
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
  typedef Random_based_perturbation<C3T3,MeshDomain,SliverCriterion> Base;
  typedef typename Base::Vertex_handle Vertex_handle;
  typedef typename Base::Cell_handle Cell_handle;
  typedef typename Base::FT FT;
  
  typedef typename C3T3::Triangulation::Geom_traits Gt;
  typedef typename Gt::Vector_3 Vector_3;
  typedef typename Gt::Point_3 Point_3;
  
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
             std::vector<Vertex_handle>& modified_vertices) const
  {
    CGAL_precondition(!slivers.empty());
    
    return apply_perturbation(v, slivers, c3t3, domain, criterion,
                              sliver_bound, modified_vertices);
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
                     std::vector<Vertex_handle>& modified_vertices) const
  {
    modified_vertices.clear();

    // Create an helper
    typedef C3T3_helpers<C3T3,MeshDomain> C3T3_helpers;
    C3T3_helpers helper(c3t3, domain);
    
    // norm depends on the local size of the mesh
    FT sq_norm = this->compute_perturbation_sq_amplitude(v, c3t3, this->sphere_sq_radius());
    const Point_3 initial_location = v->point();
    
    // Initialize loop variables
    bool min_angle_increased = false;
    Vertex_handle moving_vertex = v;
    Point_3 best_location = initial_location;
    std::vector<Vertex_handle> best_vertices;
    
    // TODO: use no_topo_change to compute new_angle without moving everything
    unsigned int try_nb = 0;
    while ( ++try_nb <= Base::max_try_nb() )
    {
      Vector_3 delta = this->random_vector(sq_norm,on_sphere_);
      
      // always from initial_location!
      Point_3 new_location = initial_location + delta;
      
      if ( c3t3.in_dimension(moving_vertex) < 3 )
        new_location = helper.project_on_surface(new_location, moving_vertex);
      
      // try to move vertex
      std::vector<Vertex_handle> tmp_mod_vertices;
      std::pair<bool,Vertex_handle> update =
        helper.update_mesh(new_location,
                           moving_vertex,
                           criterion,
                           std::back_inserter(tmp_mod_vertices));
      
      // get new vertex
      moving_vertex = update.second;
      
      if ( update.first )
      {
        min_angle_increased = true;
        best_location = moving_vertex->point();
        //best_vertices = tmp_mod_vertices;
        
        for ( typename std::vector<Vertex_handle>::iterator it = tmp_mod_vertices.begin() ;
             it != tmp_mod_vertices.end() ;
             ++it )
        {
          if ( std::find(best_vertices.begin(), best_vertices.end(), *it)
              == best_vertices.end() )
          { 
            best_vertices.push_back(*it);
          }
        }
        
        std::vector<Cell_handle> new_slivers =
          helper.incident_slivers(moving_vertex, criterion, sliver_bound);

        // If sliver number has decreased, we won
        if ( new_slivers.size() < slivers.size() )
        {
          std::copy(best_vertices.begin(),
                    best_vertices.end(),
                    std::back_inserter(modified_vertices));
          
          return std::make_pair(true,moving_vertex);
        }
      }
    }
    
    if ( min_angle_increased )
    {
      std::copy(best_vertices.begin(),
                best_vertices.end(),
                std::back_inserter(modified_vertices));
    }
    
    // Moving vertex is located on the best location
    return std::make_pair(min_angle_increased,moving_vertex);
  }
  
private:
  // If set to \c true, then random point will be generated on sphere surface.
  bool on_sphere_;
};
  
 
} // end namespace Mesh_3



} //namespace CGAL

#endif // CGAL_MESH_3_VERTEX_PERTURBATION_H
