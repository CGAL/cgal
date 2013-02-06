// Copyright (c) 2009-2010 INRIA Sophia-Antipolis (France).
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

#ifndef CGAL_MESH_3_PROTECT_EDGES_SIZING_FIELD_H
#define CGAL_MESH_3_PROTECT_EDGES_SIZING_FIELD_H

#include <CGAL/Delaunay_triangulation_3.h>

namespace {
  const double min_intersection_factor = .4; // (1-alpha)
  const double weight_modifier = .81; //0.9025;//0.81;
  const double distance_divisor = 2.1;
}

#include <cmath>
#include <algorithm>
#include <set>
#include <list>
#include <vector>
#include <stack>

#include <boost/iterator/transform_iterator.hpp>
#include <boost/optional.hpp>

#include <CGAL/Mesh_3/utilities.h>
#include <CGAL/enum.h>
#include <CGAL/iterator.h>
#include <CGAL/number_utils.h>

namespace CGAL {

namespace Mesh_3 {


template <typename C3T3, typename MeshDomain, typename SizingFunction>
class Protect_edges_sizing_field
{
  typedef Protect_edges_sizing_field          Self;
  
public:
  typedef typename C3T3::Triangulation        Tr;
  typedef typename Tr::Geom_traits            Gt;
  typedef typename Gt::FT                     FT;
  typedef typename Gt::Point_3                Weighted_point;
  typedef typename Weighted_point::Point      Bare_point;
  typedef typename Weighted_point::Weight     Weight;
  
  typedef typename C3T3::Cell_handle          Cell_handle;
  typedef typename C3T3::Vertex_handle        Vertex_handle;
  typedef typename C3T3::Triangulation        Triangulation;
  typedef typename C3T3::Edge                 Edge;
  
  typedef typename MeshDomain::Curve_segment_index  Curve_segment_index;
  typedef typename MeshDomain::Corner_index         Corner_index;
  typedef typename MeshDomain::Index                Index;
  
public:
  Protect_edges_sizing_field(C3T3& c3t3,
                             const MeshDomain& domain,
                             SizingFunction size=SizingFunction());
  
  void operator()(const bool refine=true);
  
private:
  typedef std::vector<std::pair<Curve_segment_index,Bare_point> >    Incident_edges;
  typedef std::vector<Vertex_handle>                                 Vertex_vector;
  typedef std::vector<std::pair<Vertex_handle,Curve_segment_index> > Incident_vertices;
  
private:
  /// Insert corners of the mesh
  void insert_corners();
  
  /// Insert balls on every edge
  void insert_balls_on_edges();
  
  /// Refine balls
  void refine_balls();

  /// Returns vertex which corresponds to corner located at point p
  Vertex_handle get_vertex_corner_from_point(const Bare_point& p,
                                             const Index& p_index) const;
  
  /// Insert point p as a curve segment point
  Vertex_handle insert_curve_point(const Bare_point& p, const Index& p_index);
  
  /// Insert point(p,w) into triangulation and set its dimension to \c dim and
  /// it's index to \c index.
  /// The newly created handle is returned
  Vertex_handle insert_point(const Bare_point& p,
                             const Weight& w,
                             int dim,
                             const Index& index);

  /**
   * Insert point(p,w) into triangulation and set its dimension to \c dim and
   * it's index to \c index.
   * The newly created handle is returned
   * This function also ensures that point(p,w) will not be inside a sphere,
   * and that no point of the triangulation will be inside its sphere.
   */
  Vertex_handle smart_insert_point(const Bare_point& p,
                                   Weight w,
                                   int dim,
                                   const Index& index);
  
  
  /// Insert balls between points which are pointed by handles \c vp and \c vq
  /// on curve identified by \c curve_index
  void insert_balls(const Vertex_handle& vp,
                    const Vertex_handle& vq,
                    const Curve_segment_index& curve_index);

  /**
   * Insert balls
   * Preconditions:
   *  - size_p < size_q
   *  - pq_geodesic > 0
   */              
  void insert_balls(const Vertex_handle& vp,
                    const Vertex_handle& vq,
                    const FT size_p,
                    const FT size_q,
                    const FT pq_geodesic,
                    const CGAL::Sign distance_sign,
                    const Curve_segment_index& curve_index);
  
  /// Returns true if balls of \c va and \c vb intersect, and (va,vb) is not
  /// an edge of the complex
  bool non_adjacent_but_intersect(const Vertex_handle& va,
                                  const Vertex_handle& vb) const;
  
  /// Change size of the ball of vertex \c v.
  Vertex_handle change_ball_size(const Vertex_handle& v, const FT size);
  
  
  /// Returns true if balls of v1 and v2 intersect "enough"
  bool is_sampling_dense_enough(const Vertex_handle& v1,
                                const Vertex_handle& v2) const;
  
  /// Takes an iterator on Vertex_handle as input and check if the sampling
  /// of those vertices is ok. If not, fix it.
  void check_and_repopulate_edges();

  /// Checks if vertex \c v is well sampled, and if its not the case, fix it.
  /// Fills out with deleted vertices during this process. out value type
  /// is Vertex_handle.
  template <typename OutputIterator>
  OutputIterator
  check_and_fix_vertex_along_edge(const Vertex_handle& v, OutputIterator out);
  
  /// Walk along edge from \c start, following the direction \c start to 
  /// \c next, and fills \c out with the vertices which do not fullfill 
  /// the sampling conditions
  template <typename OutputIterator>
  OutputIterator
  walk_along_edge(const Vertex_handle& start,
                  const Vertex_handle& next,
                  const bool test_sampling,
                  OutputIterator out) const;
  
  /// Returns next vertex along edge, i.e vertex after \c start, following
  /// the direction from \c previous to \c start
  /// \pre (previous,start) is in c3t3
  Vertex_handle next_vertex_along_edge(const Vertex_handle& start,
                                       const Vertex_handle& previous) const;
  
  /// Replace vertices between ]begin,last[ by new vertices, along curve
  /// identified by \c curve_index
  /// The value type of InputIterator is Vertex_handle.
  template <typename InputIterator>
  void repopulate(InputIterator begin,
                  InputIterator last,
                  const Curve_segment_index& index);
  
  template <typename InputIterator>
  void
  analyze_and_repopulate(InputIterator begin,
                         InputIterator last,
                         const Curve_segment_index& index);
  
  /// Checks if \c v2 size is compatible (i.e. greater) with the linear
  /// interpolation of sizes of \c v1 and \c v3 
  bool is_sizing_field_correct(const Vertex_handle& v1,
                               const Vertex_handle& v2,
                               const Vertex_handle& v3) const;
  
  /// Repopulate all incident curve around corner \c v
  /// \pre \c v is a corner of c3t3 
  template <typename OutputIterator>
  OutputIterator
  repopulate_edges_around_corner(const Vertex_handle& v, OutputIterator out);
  
  /// Returns true if edge with index \c curve_index is already treated
  bool is_treated(const Curve_segment_index& curve_index) const
  {
    return ( treated_edges_.find(curve_index) != treated_edges_.end() );
  }
  
  /// Set edge with index \c curve_index as treated
  void set_treated(const Curve_segment_index& curve_index)
  {
    treated_edges_.insert(curve_index);
  }
  
  /// Compute euclidean distance between bare points of \c va and \c vb
  FT compute_distance(const Vertex_handle& va, const Vertex_handle& vb) const
  {
    return compute_distance(va->point().point(), vb->point().point());
  }
  
  /// Compute euclidean distance between bare points \c and \c q
  FT compute_distance(const Bare_point& p, const Bare_point& q) const
  {
    return CGAL::sqrt(Gt().compute_squared_distance_3_object()(p,q));
  }
  
  /// Returns the radius of the ball of vertex \c v
  FT get_size(const Vertex_handle& v) const
  {
    return CGAL::sqrt(v->point().weight());
  }
  
private:
  C3T3& c3t3_;
  const MeshDomain& domain_;
  SizingFunction size_;
  std::set<Curve_segment_index> treated_edges_;
  std::set<Vertex_handle> unchecked_vertices_;
};


template <typename C3T3, typename MD, typename Sf>
Protect_edges_sizing_field<C3T3, MD, Sf>::
Protect_edges_sizing_field(C3T3& c3t3, const MD& domain, Sf size)
  : c3t3_(c3t3)
  , domain_(domain)
  , size_(size)
{
}


template <typename C3T3, typename MD, typename Sf>
void
Protect_edges_sizing_field<C3T3, MD, Sf>::
operator()(const bool refine)
{
#ifdef CGAL_MESH_3_VERBOSE
  std::cerr << "Inserting protection balls..." << std::endl;
#endif
  
  // Insert 0-dimensional features
  insert_corners();
#ifdef CGAL_MESH_3_VERBOSE
  std::cerr << "insert_corners() done. Nb of points in triangulation: "
            << c3t3_.triangulation().number_of_vertices() << std::endl;
#endif
  
  // Insert 1-dimensional features
  insert_balls_on_edges();
#ifdef CGAL_MESH_3_VERBOSE
  std::cerr << "insert_balls_on_edges() done. Nb of points in triangulation: "
            << c3t3_.triangulation().number_of_vertices() << std::endl;
#endif
  
  // Solve problems
  if ( refine )
  { 
    refine_balls();
#ifdef CGAL_MESH_3_VERBOSE
    std::cerr << "refine_balls() done. Nb of points in triangulation: "
              << c3t3_.triangulation().number_of_vertices() << std::endl;
#endif
    CGAL_assertion(c3t3_.is_valid());
  }
  
#ifdef CGAL_MESH_3_VERBOSE
  std::cerr << std::endl;
#endif
}


template <typename C3T3, typename MD, typename Sf>
void
Protect_edges_sizing_field<C3T3, MD, Sf>::
insert_corners()
{
  // Gt is a traits class for regular triangulations, and Gt::Kernel is its
  // CGAL kernel.
  typedef CGAL::Delaunay_triangulation_3<typename Gt::Kernel> Dt;

  // Iterate on domain corners
  typedef std::vector< std::pair<Corner_index, Bare_point> > Initial_corners;
  Initial_corners corners;
  domain_.get_corners(std::back_inserter(corners));
  
  Dt dt;
  for ( typename Initial_corners::iterator it = corners.begin(),
       end = corners.end() ; it != end ; ++it )
  {
    const Bare_point& p = it->second;
    dt.insert(p);
  }

  for ( typename Initial_corners::iterator cit = corners.begin(),
          end = corners.end() ; cit != end ; ++cit )
  {
    const Bare_point& p = cit->second;
    Index p_index = domain_.index_from_corner_index(cit->first);
    
    // Get weight (ball radius is given by size_ function)
    FT w = CGAL::square(size_(p, 0, p_index));

    // the following lines ensure that the weight w is small enough so that
    // corners balls do not intersect
    if(dt.number_of_vertices() >= 2)
    {

      typename Dt::Vertex_handle vh;
      CGAL_assertion_code( bool p_found= )
        dt.is_vertex(p, vh);
      CGAL_assertion(p_found);
      std::vector<typename Dt::Vertex_handle> vs;
      vs.reserve(32);
      dt.finite_adjacent_vertices(vh, std::back_inserter(vs));
      CGAL_assertion(!vs.empty());
      typename Dt::Point nearest = vs[0]->point();
      typename Gt::Compare_distance_3 compare_dist = 
        c3t3_.triangulation().geom_traits().compare_distance_3_object();
      for (typename std::vector<typename Dt::Vertex_handle>::const_iterator
             it = vs.begin(); it != vs.end(); ++it) 
      {
        if(compare_dist(p, (*it)->point(), nearest) == CGAL::SMALLER) {
          // 	    std::cerr << "  nearest!\n";
          nearest =  (*it)->point();
        }
      }
      typename Gt::Compute_squared_distance_3 squared_distance = 
        c3t3_.triangulation().geom_traits().compute_squared_distance_3_object();
      const FT nearest_sq_dist = squared_distance( nearest, p);
      
      w = (std::min)(w, nearest_sq_dist / FT(9));
    }
    
    // Insert corner with ball (dim is zero because p is a corner)
    Vertex_handle v = smart_insert_point(p, w, 0, p_index);
    c3t3_.add_to_complex(v,cit->first);
  }
}


template <typename C3T3, typename MD, typename Sf>
typename Protect_edges_sizing_field<C3T3, MD, Sf>::Vertex_handle
Protect_edges_sizing_field<C3T3, MD, Sf>::
insert_point(const Bare_point& p, const Weight& w, int dim, const Index& index)
{
  typedef typename Tr::size_type size_type;
  CGAL_USE_TYPE(size_type);
  
  // Insert point
  CGAL_assertion_code(size_type nb_vertices_before = c3t3_.triangulation().number_of_vertices());
  Vertex_handle v = c3t3_.triangulation().insert(Weighted_point(p,w*weight_modifier));
  
  // If point insertion created an hidden ball, fail
  CGAL_assertion ( Vertex_handle() != v );
  CGAL_assertion ( c3t3_.triangulation().number_of_vertices() == (nb_vertices_before+1) );

#ifdef PROTECTION_DEBUG
  std::cerr << "Insertion of protecting ball "
            << Weighted_point(p,w*weight_modifier)
            << " on curve #" << index << std::endl;
#endif
  
  c3t3_.set_dimension(v,dim);
  c3t3_.set_index(v,index);
  
  unchecked_vertices_.insert(v);
  
  return v;
}

  
template <typename C3T3, typename MD, typename Sf>
typename Protect_edges_sizing_field<C3T3, MD, Sf>::Vertex_handle
Protect_edges_sizing_field<C3T3, MD, Sf>::
smart_insert_point(const Bare_point& p, Weight w, int dim, const Index& index)
{
  const Tr& tr = c3t3_.triangulation();
  typename Gt::Compute_squared_distance_3 sq_distance =
    tr.geom_traits().compute_squared_distance_3_object();
  
  bool add_handle_to_unchecked = false;
  
  if ( tr.dimension() > 2 ) 
  {
    // Check that new point will not be inside a power sphere
    Cell_handle ch = tr.locate(p);
    Vertex_handle nearest_vh = tr.nearest_power_vertex(p, ch);
    FT sq_d = sq_distance(p, nearest_vh->point().point());
    CGAL_assertion( sq_d > 0 );
    
    while ( nearest_vh->point().weight() > sq_d )
    {
      // Adapt size
      Vertex_handle new_vh = change_ball_size(nearest_vh, CGAL::sqrt(sq_d));
      ch = tr.locate(p,new_vh);
      
      // Iterate
      nearest_vh = tr.nearest_power_vertex(p, ch);
      sq_d = sq_distance(p, nearest_vh->point().point());
      CGAL_assertion( sq_d > 0 );
    }

    // Change w in order to be sure that no existing point will be included
    // in (p,w)
    std::vector<Cell_handle> cells_in_conflicts;
    std::set<Vertex_handle> vertices_in_conflict_zone;
    tr.find_conflicts(Weighted_point(p, w), ch,
                      CGAL::Emptyset_iterator(),
                      std::back_inserter(cells_in_conflicts),
                      CGAL::Emptyset_iterator());

    for(typename std::vector<Cell_handle>::const_iterator 
          it = cells_in_conflicts.begin(),
          end = cells_in_conflicts.end(); it != end; ++it) 
    {
      for(int i = 0, d = tr.dimension(); i <= d; ++i) {
        vertices_in_conflict_zone.insert((*it)->vertex(i));
      }
    }
    FT min_sq_d = w;
    for(typename std::set<Vertex_handle>::const_iterator 
          it = vertices_in_conflict_zone.begin(),
          end = vertices_in_conflict_zone.end(); it != end ; ++it )
    {
      min_sq_d = (std::min)(min_sq_d, sq_distance(p, (*it)->point().point()));
    }

    if ( w > min_sq_d )
    { 
      w = min_sq_d;
      add_handle_to_unchecked = true;
    }
    CGAL_assertion_code(std::vector<Vertex_handle> hidden_vertices;);
    CGAL_assertion_code(tr.vertices_inside_conflict_zone(Weighted_point(p, w),
                                                         ch,
                                                         std::back_inserter(hidden_vertices)));
    CGAL_assertion(hidden_vertices.empty());
  }
  else // tr.dimension() <= 2
  {
    // change size of existing balls which include p
    bool restart = true;
    while ( restart )
    {
      restart = false;
      for ( typename Tr::Finite_vertices_iterator it = tr.finite_vertices_begin(),
           end = tr.finite_vertices_end() ; it != end ; ++it )
      {
        FT sq_d = sq_distance(p, it->point().point());
        if ( it->point().weight() > sq_d )
        { 
          change_ball_size(it, CGAL::sqrt(sq_d));
          restart = true;
          break;
        }
      }
    }
    
    // Change w in order to be sure that no existing point will be included
    // in (p,w)
    for ( typename Tr::Finite_vertices_iterator it = tr.finite_vertices_begin(),
         end = tr.finite_vertices_end() ; it != end ; ++it )
    {
      FT sq_d = sq_distance(p, it->point().point());
      if ( w > sq_d )
      { 
        w = sq_d;
        add_handle_to_unchecked = true;
      }
    }    
  }

  FT w_max = CGAL::square(size_(p, dim, index));
  if(w > w_max) {
    w = w_max;
    add_handle_to_unchecked = true;
  }
  
  Vertex_handle v = insert_point(p,w,dim,index);
  if ( add_handle_to_unchecked ) { unchecked_vertices_.insert(v); }
  
  return v;
}

  
template <typename C3T3, typename MD, typename Sf>
void
Protect_edges_sizing_field<C3T3, MD, Sf>::
insert_balls_on_edges()
{
  // Get features
  typedef CGAL::cpp11::tuple<Curve_segment_index,
                             std::pair<Bare_point,Index>,
                             std::pair<Bare_point,Index> >    Feature_tuple;
  typedef std::vector<Feature_tuple>                          Input_features;
  
  Input_features input_features;
  domain_.get_curve_segments(std::back_inserter(input_features));
  
  // Interate on edges
  for ( typename Input_features::iterator fit = input_features.begin(),
       end = input_features.end() ; fit != end ; ++fit )
  {
    const Curve_segment_index& curve_index = CGAL::cpp11::get<0>(*fit);
    if ( ! is_treated(curve_index) )
    {
#ifdef PROTECTION_DEBUG
      std::cerr << "** treat curve #" << curve_index << std::endl;
#endif
      const Bare_point& p = CGAL::cpp11::get<1>(*fit).first;
      const Bare_point& q = CGAL::cpp11::get<2>(*fit).first; 
      
      const Index& p_index = CGAL::cpp11::get<1>(*fit).second;
      const Index& q_index = CGAL::cpp11::get<2>(*fit).second;
      
      Vertex_handle vp,vq;
      if ( ! domain_.is_cycle(p, curve_index) )
      {
        vp = get_vertex_corner_from_point(p,p_index);
        vq = get_vertex_corner_from_point(q,q_index);
      }
      else
      {
        vp = insert_curve_point(p,p_index);
        vq = vp;
      }
      
      // Insert balls and set treated
      insert_balls(vp, vq, curve_index);
      set_treated(curve_index);
    }
  }
}


template <typename C3T3, typename MD, typename Sf>
typename Protect_edges_sizing_field<C3T3, MD, Sf>::Vertex_handle
Protect_edges_sizing_field<C3T3, MD, Sf>::
get_vertex_corner_from_point(const Bare_point& p, const Index& p_index) const
{
  // Get vertex_handle associated to corner (dim=0) point
  FT size_p = size_(p, 0, p_index);
  Vertex_handle v;
  CGAL_assertion_code( bool q_finded = )
  c3t3_.triangulation().is_vertex(Weighted_point(p,size_p), v);
  CGAL_assertion( q_finded );
  return v;
}


template <typename C3T3, typename MD, typename Sf>
typename Protect_edges_sizing_field<C3T3, MD, Sf>::Vertex_handle
Protect_edges_sizing_field<C3T3, MD, Sf>::
insert_curve_point(const Bare_point& p, const Index& p_index)
{
  const int p_dim = 1;
  FT p_weight = CGAL::square(size_(p, p_dim, p_index));
  return smart_insert_point(p, p_weight, p_dim, p_index);
}
  
  
template <typename C3T3, typename MD, typename Sf>
void
Protect_edges_sizing_field<C3T3, MD, Sf>::
insert_balls(const Vertex_handle& vp,
             const Vertex_handle& vq,
             const Curve_segment_index& curve_index)
{
  // Get size of p & q
  const Bare_point& p = vp->point().point();
  const Bare_point& q = vq->point().point();
  
  const FT sp = get_size(vp);
  const FT sq = get_size(vq);
  
  // Compute geodesic distance
  const FT pq_geo_signed = domain_.geodesic_distance(p, q, curve_index);
  const CGAL::Sign d_sign = CGAL::sign(pq_geo_signed);
  const FT pq_geo = CGAL::abs(pq_geo_signed);

  // Insert balls
  return (sp <= sq) ? insert_balls(vp, vq, sp, sq, pq_geo, d_sign, curve_index)
                    : insert_balls(vq, vp, sq, sp, pq_geo, -d_sign, curve_index);  
}



template <typename C3T3, typename MD, typename Sf>
void
Protect_edges_sizing_field<C3T3, MD, Sf>::
insert_balls(const Vertex_handle& vp,
             const Vertex_handle& vq,
             const FT sp,
             const FT sq,
             const FT d,
             const CGAL::Sign d_sign,
             const Curve_segment_index& curve_index)
{
  CGAL_precondition(d > 0);
  CGAL_precondition(sp <= sq);
  CGAL_precondition(sp > 0);

  // Notations:
  // sp = size_p,   sq = size_q,   d = pq_geodesic
  // n = nb_points,   r = delta_step_size
  // 
  // Hypothesis:
  // sp <= sq
  //
  // Let's define
  // P0 = p, Pn+1 = q, d(Pi,Pi+1) = ai
  //
  // The following constraints should be verified:
  // a0 = sp + r, an = sq,
  // ai+1 = ai + r
  // d = Sum(ai)
  //
  // The following could be obtained:
  // r = (sq - sp) / (n+1)
  // n = 2(d-sq) / (sp+sq)
  //
  // =======================
  // Calculus details:
  // ai+1 = ai + r
  // ai+1 = a0 + r*(i+1)
  //   an = a0 + r*n
  //   sq = sp + r + r*n
  //    r = (sq-sp) / (n+1)
  //
  //   d = Sum(ai)
  //   d = Sum(sp + (i+1)*r)
  //   d = (n+1)*sp + (n+1)(n+2)/2 * r
  //   d = (n+1)*sp + (n+1)(n+2)/2 * (sq-sp) / (n+1)
  // 2*d = 2(n+1)*sp + (n+2)*sq - (n+2)*sp
  // 2*d = n*sp + (n+2)*sq
  //   n = 2(d-sq) / (sp+sq)
  // =======================
  
  int n = static_cast<int>(std::floor(FT(2)*(d-sq) / (sp+sq))+.5);
  FT r = (sq - sp) / FT(n+1);
  
  // Adjust size of steps, D = covered distance
  FT D = sp*FT(n+1) + FT((n+1)*(n+2)) / FT(2) * r ;
  
  FT dleft_frac = d / D;
  CGAL_assertion(dleft_frac >= 1);
  
  // Initialize step sizes
  FT step_size = sp + r;
  FT norm_step_size = dleft_frac * step_size;
  
  // Initial distance
  FT d_signF = static_cast<FT>(d_sign);
  FT pt_dist = d_signF * norm_step_size;
  Vertex_handle prev = vp;
  const Bare_point& p = vp->point().point();

  // If there is some place to insert one point, insert it
  if ( (0 == n) && (d >= sp+sq) )
  {
    n = 1;
    step_size = sp + (d-sp-sq) / FT(2);
    pt_dist = d_signF * step_size;
    norm_step_size = step_size;
  }
  
  // Launch balls
  for ( int i = 1 ; i <= n ; ++i )
  {
    // New point position
    Bare_point new_point =
      domain_.construct_point_on_curve_segment(p, curve_index, pt_dist);
    
    // Weight (use as size the min between norm_step_size and linear interpolation)
    FT current_size = (std::min)(norm_step_size, sp + CGAL::abs(pt_dist)/d*(sq-sp));
    FT point_weight = current_size * current_size;
    
    // Index and dimension
    Index index = domain_.index_from_curve_segment_index(curve_index);
    int dim = 1; // new_point is on edge
    
    // Insert point into c3t3
    Vertex_handle new_vertex = smart_insert_point(new_point, point_weight, dim, index);
    
    // Add edge to c3t3
    c3t3_.add_to_complex(prev, new_vertex, curve_index);
    prev = new_vertex;
    
    // Step size
    step_size += r;
    norm_step_size = dleft_frac * step_size;
    
    // Increment distance
    pt_dist += d_signF * norm_step_size;
  }
  
  // Insert last edge into c3t3
  // Warning: if vp==vq (cycle) and if only 1 point was inserted,
  // then (prev,vp) == (prev,vq)
  if ( vp != vq || n > 1 )
  {
    c3t3_.add_to_complex(prev, vq, curve_index);
  }
}
  
  
template <typename C3T3, typename MD, typename Sf>
void
Protect_edges_sizing_field<C3T3, MD, Sf>::
refine_balls()
{
  Triangulation& tr = c3t3_.triangulation();
  
  // Loop
  bool restart = true;
  int nb=0;
  while ( (!unchecked_vertices_.empty() || restart) && nb<29)
  {
#ifdef PROTECTION_DEBUG
    std::cerr << "RESTART REFINE LOOP (" << nb << ")\n"
              << "\t unchecked_vertices size: " << unchecked_vertices_.size() <<"\n";
#endif
    ++nb;
    restart = false;
    std::map<Vertex_handle, FT> new_sizes;
    
    for(typename Tr::Finite_edges_iterator eit = tr.finite_edges_begin(), 
        end = tr.finite_edges_end(); eit != end; ++eit)
    {
      const Vertex_handle& va = eit->first->vertex(eit->second);
      const Vertex_handle& vb = eit->first->vertex(eit->third);
      
      // If those vertices are not adjacent 
      if( non_adjacent_but_intersect(va, vb) )
      {
        // Compute correct size of balls
        const FT ab = compute_distance(va,vb);
        FT sa_new = (std::min)(ab/distance_divisor, get_size(va));
        FT sb_new = (std::min)(ab/distance_divisor, get_size(vb));
        
        // In case of va or vb have already been in conflict, keep minimal size
        if ( new_sizes.find(va) != new_sizes.end() )
        { sa_new = (std::min)(sa_new, new_sizes[va]); }
        
        if ( new_sizes.find(vb) != new_sizes.end() )
        { sb_new = (std::min)(sb_new, new_sizes[vb]); }
        
        // Store new_sizes for va and vb
        if ( sa_new != get_size(va) )
        { new_sizes[va] = sa_new; }
        
        if ( sb_new != get_size(vb) )
        { new_sizes[vb] = sb_new; }
        
        // Loop will have to be run again
        restart = true;
      }
    }
    
    // Update size of balls
    for ( typename std::map<Vertex_handle,FT>::iterator it = new_sizes.begin(),
         end = new_sizes.end() ; it != end ; ++it)
    {
      // Set size of the ball to new value
      change_ball_size(it->first,it->second);
    }
    
    // Check edges
    check_and_repopulate_edges();
  }
}


template <typename C3T3, typename MD, typename Sf>
bool
Protect_edges_sizing_field<C3T3, MD, Sf>::
non_adjacent_but_intersect(const Vertex_handle& va, const Vertex_handle& vb) const
{
  if ( ! c3t3_.is_in_complex(va,vb) )
  {
    typename Gt::Construct_sphere_3 sphere = 
      c3t3_.triangulation().geom_traits().construct_sphere_3_object();
    
    typename Gt::Do_intersect_3 do_intersect = 
      c3t3_.triangulation().geom_traits().do_intersect_3_object();
    
    const Bare_point& a = va->point().point();
    const Bare_point& b = vb->point().point();
    
    const FT& sq_ra = va->point().weight();
    const FT& sq_rb = vb->point().weight();
    
    return do_intersect(sphere(a, sq_ra), sphere(b, sq_rb));
  }
  
  return false;
}


template <typename C3T3, typename MD, typename Sf>
typename Protect_edges_sizing_field<C3T3, MD, Sf>::Vertex_handle
Protect_edges_sizing_field<C3T3, MD, Sf>::
change_ball_size(const Vertex_handle& v, const FT size)
{
  // Check if there is something to do
  if ( get_size(v) == size )
  { return v; }
  
  // Get incident vertices along c3t3 edge
  Incident_vertices incident_vertices;
  c3t3_.adjacent_vertices_in_complex(v, std::back_inserter(incident_vertices));
  
  // Remove incident edges from complex
  for ( typename Incident_vertices::iterator vit = incident_vertices.begin(),
       vend = incident_vertices.end() ; vit != vend ; ++vit )
  {
    c3t3_.remove_from_complex(v,vit->first);
  }
  
  
  // Store point data
  Index index = v->index();
  int dim = v->in_dimension();
  Bare_point p = v->point().point();

  // Remove v from corners
  boost::optional<Corner_index> corner_index;
  if ( c3t3_.is_in_complex(v) )
  {
    corner_index = c3t3_.corner_index(v);
    c3t3_.remove_from_complex(v);
  }
  // Change v size
  c3t3_.triangulation().remove(v);
  Vertex_handle new_v = insert_point(p, size*size, dim, index);
  
  // TODO: ensure that this condition is always satisfied (Pedro's code ?)
  CGAL_assertion(v==new_v);
  //new_v->set_meshing_info(size*size);

  // Restore v in corners
  if ( corner_index )
  {
    c3t3_.add_to_complex(new_v,*corner_index);
  }
  
  // Restore c3t3 edges
  for ( typename Incident_vertices::iterator it = incident_vertices.begin(),
       end = incident_vertices.end() ; it != end ; ++it )
  {
    // Restore connectivity in c3t3
    c3t3_.add_to_complex(new_v, it->first, it->second);
  }

  // Update unchecked vertices
  unchecked_vertices_.erase(v);
  unchecked_vertices_.insert(new_v);
  return new_v;
}


template <typename C3T3, typename MD, typename Sf>
void
Protect_edges_sizing_field<C3T3, MD, Sf>::
check_and_repopulate_edges()
{
  std::set<Vertex_handle> vertices;
  std::copy( unchecked_vertices_.begin(), unchecked_vertices_.end(),
             std::inserter(vertices,vertices.begin()) );
  
  unchecked_vertices_.clear();
  
  // Fix edges
  while ( !vertices.empty() )
  {
    Vertex_handle v = *vertices.begin();
    vertices.erase(vertices.begin());
    
    Vertex_vector erased_vertices;
    check_and_fix_vertex_along_edge(v, std::back_inserter(erased_vertices));
    
    for ( typename Vertex_vector::iterator vit = erased_vertices.begin(), 
         vend = erased_vertices.end() ; vit != vend ; ++vit )
    {
      if ( c3t3_.in_dimension(*vit) > 0 )
      {
        vertices.erase(*vit);
      }
    }
  }  
}

  
template <typename C3T3, typename MD, typename Sf>
template <typename OutputIterator>
OutputIterator
Protect_edges_sizing_field<C3T3, MD, Sf>::
check_and_fix_vertex_along_edge(const Vertex_handle& v, OutputIterator out)
{
  // If v is a corner, then all incident edges have to be checked
  if ( c3t3_.is_in_complex(v) )
  {
    return repopulate_edges_around_corner(v, out);
  }
  
  // Get incident vertices along c3t3 edge
  Incident_vertices incident_vertices;
  c3t3_.adjacent_vertices_in_complex(v, std::back_inserter(incident_vertices));
  CGAL_assertion(incident_vertices.size() == 2);
  
  // Walk along edge to find the edge piece which is not correctly sampled
  typedef std::list<Vertex_handle> Vertex_list;
  Vertex_list to_repopulate;
  to_repopulate.push_front(v);
  
  const Vertex_handle& previous = incident_vertices.front().first;
  const Vertex_handle& next = incident_vertices.back().first;

  // Walk following direction (v,previous)
  walk_along_edge(v, previous, true, std::front_inserter(to_repopulate));
  
  // Check whether a complete circle has been discovered or not
  if (   to_repopulate.size() == 1
      || to_repopulate.front() != to_repopulate.back() )
  {
    // Walk in other direction (v,next)
    walk_along_edge(v, next, true, std::back_inserter(to_repopulate));    
  }
  
  // If only v is in to_repopulate, there is nothing to do
  if ( to_repopulate.size() == 1 )
  {
    *out++ = *to_repopulate.begin();
    return out;
  }
  
  // Store erased vertices
  std::copy(to_repopulate.begin(), to_repopulate.end(), out);

  // Repopulate edge
  const Curve_segment_index& curve_index = incident_vertices.front().second;
  
  analyze_and_repopulate(to_repopulate.begin(),
                         --to_repopulate.end(),
                         curve_index);
  
  return out;
}


template <typename C3T3, typename MD, typename Sf>
bool
Protect_edges_sizing_field<C3T3, MD, Sf>::
is_sampling_dense_enough(const Vertex_handle& v1, const Vertex_handle& v2) const
{
  CGAL_precondition(c3t3_.is_in_complex(v1,v2));

  // Get sizes
  FT size_v1 = get_size(v1);
  FT size_v2 = get_size(v2);
  FT distance_v1v2 = compute_distance(v1,v2);
  
  // Ensure size_v1 > size_v2
  if ( size_v1 < size_v2 ) { std::swap(size_v1, size_v2); }
  
  // Check if balls intersect
  return distance_v1v2 < (FT(min_intersection_factor) * size_v2 + size_v1);  
}


template <typename C3T3, typename MD, typename Sf>
template <typename OutputIterator>
OutputIterator
Protect_edges_sizing_field<C3T3, MD, Sf>::
walk_along_edge(const Vertex_handle& start, const Vertex_handle& next,
                bool /*test_sampling*/,
                OutputIterator out) const
{
  CGAL_precondition( c3t3_.is_in_complex(start, next) );
  
  Vertex_handle previous = start;
  Vertex_handle current = next;
  
  // Walk along edge since a corner is encountered or the balls of previous
  // and current intersects enough
  while ( ! is_sampling_dense_enough(previous, current) )
  {
    *out++ = current;
    
    // Don't go through corners
    if ( c3t3_.is_in_complex(current) || current == start )
    {
      break;
    }
    
    // Get next vertex along edge
    Vertex_handle next = next_vertex_along_edge(current,previous);
    previous = current;
    current = next;
  }
  
  return out;
}


template <typename C3T3, typename MD, typename Sf>
typename Protect_edges_sizing_field<C3T3, MD, Sf>::Vertex_handle
Protect_edges_sizing_field<C3T3, MD, Sf>::
next_vertex_along_edge(const Vertex_handle& start,
                       const Vertex_handle& previous) const
{
  CGAL_precondition( c3t3_.is_in_complex(start, previous) );
  CGAL_precondition( ! c3t3_.is_in_complex(start) );
  
  Incident_vertices incident_vertices;
  c3t3_.adjacent_vertices_in_complex(start,
                                  std::back_inserter(incident_vertices));
  CGAL_assertion(incident_vertices.size() == 2);
  
  if ( incident_vertices.front().first == previous )
  { 
    return incident_vertices.back().first;
  }
  else
  { 
    return incident_vertices.front().first;
  }
}

  
template <typename C3T3, typename MD, typename Sf>
template <typename InputIterator>
void
Protect_edges_sizing_field<C3T3, MD, Sf>::
repopulate(InputIterator begin, InputIterator last, const Curve_segment_index& index)
{
  CGAL_assertion( std::distance(begin,last) >= 0 );
  
  // May happen
  if ( begin == last ) { return; }
  
  // Valid because begin < last
  InputIterator current = begin;
  InputIterator previous = current++;
  
  // Remove edges from c3t3.
  while ( current != last )
  {
    c3t3_.remove_from_complex(*previous++, *current++);
  }
  // Keep a point between current and last
  Bare_point point_through = (*previous)->point().point();
  // Remove last edge
  c3t3_.remove_from_complex(*previous, *current);
  
  // Remove vertices (don't remove the first one and the last one)
  current = begin;
  while ( ++current != last )
  {
    c3t3_.triangulation().remove(*current);
  }
  
  // If edge is a cycle, order the iterators according to the orientation of
  // the cycle
  if (  domain_.is_cycle((*begin)->point().point(), index) 
      && domain_.distance_sign_along_cycle((*begin)->point().point(),
                                           point_through,
                                           (*last)->point().point(),
                                           index ) != CGAL::POSITIVE )
  {
    std::swap(begin,last);
  }
  
  // Repopulate edge
  insert_balls(*begin, *last, index);
}


template <typename C3T3, typename MD, typename Sf>
template <typename InputIterator>
void
Protect_edges_sizing_field<C3T3, MD, Sf>::
analyze_and_repopulate(InputIterator begin, InputIterator last, const Curve_segment_index& index)
{
  CGAL_assertion( std::distance(begin,last) >= 0 );
  
  // May happen
  if ( begin == last ) { return; }
  if ( std::distance(begin,last) == 1 )
  {
    repopulate(begin, last, index);
    return;
  }
  
  // Here std::distance(begin,last) > 1
  
  // ch_stack is the stack filled with the convex hull of element size.
  // The goal is to ensure that no ball will have its size increased
  std::stack<InputIterator> ch_stack;
  InputIterator current = begin;
  ch_stack.push(current);
  ch_stack.push(++current);

  // Compute the convex hull of the size of elements
  while ( ++current != last )
  {
    // Get last element of the stack
    InputIterator previous = ch_stack.top();
    ch_stack.pop();
    
    // If (prevprev, prev, current) is ok, then go one step forward, i.e. check
    // (prevprevprev, prevprev, current)
    while (   !ch_stack.empty() 
           && is_sizing_field_correct(*ch_stack.top(),*previous,*current) )
    {
      previous = ch_stack.top();
      ch_stack.pop();
    }
    
    // Push in the stack the furthest good element (previous)
    // and current element
    ch_stack.push(previous);   
    ch_stack.push(current);
  }
  
  // Insert last element
  ch_stack.push(last);
  
  // Repopulate edge segments
  current = ch_stack.top();
  ch_stack.pop();
  while ( !ch_stack.empty() )
  {
    InputIterator next = ch_stack.top();
    ch_stack.pop();
    // Iterators are on the reverse order in the stack, thus use [next,current]
    repopulate(next, current, index);
    current = next;
  }  
}
  
template <typename C3T3, typename MD, typename Sf>
bool
Protect_edges_sizing_field<C3T3, MD, Sf>::
is_sizing_field_correct(const Vertex_handle& v1,
                        const Vertex_handle& v2,
                        const Vertex_handle& v3) const
{
  FT s1 = get_size(v1);
  FT s2 = get_size(v2);
  FT s3 = get_size(v3);
  FT D = compute_distance(v1,v3);
  FT d = compute_distance(v1,v2);
  
  return ( s2 >= (s1 + d/D*(s3-s1)) );
}

  
template <typename C3T3, typename MD, typename Sf>
template <typename OutputIterator>
OutputIterator
Protect_edges_sizing_field<C3T3, MD, Sf>::
repopulate_edges_around_corner(const Vertex_handle& v, OutputIterator out)
{
  CGAL_precondition(c3t3_.is_in_complex(v));
  
  Incident_vertices incident_vertices;
  c3t3_.adjacent_vertices_in_complex(v, std::back_inserter(incident_vertices));
  
  for ( typename Incident_vertices::iterator vit = incident_vertices.begin(),
       vend = incident_vertices.end() ; vit != vend ; ++vit )
  {
    const Vertex_handle& next = vit->first;
    const Curve_segment_index& curve_index = vit->second;
    
    // Walk along each incident edge of the corner
    Vertex_vector to_repopulate;
    to_repopulate.push_back(v);
    walk_along_edge(v, next, true, std::back_inserter(to_repopulate));

    // Return erased vertices
    std::copy(to_repopulate.begin(), to_repopulate.end(), out);

    // Repopulate
    analyze_and_repopulate(to_repopulate.begin(), --to_repopulate.end(), curve_index);
  }
  
  return out;
}
  
} // end namespace Mesh_3


} //namespace CGAL

#endif // CGAL_MESH_3_PROTECT_EDGES_SIZING_FIELD_H
