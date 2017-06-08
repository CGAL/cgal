// Copyright (c) 2009-2010 INRIA Sophia-Antipolis (France).
// Copyright (c) 2010-2013 GeometryFactory Sarl (France)
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
// Author(s)     : Stephane Tayeb, Laurent Rineau
//
//******************************************************************************
// File Description : 
//******************************************************************************

#ifndef CGAL_MESH_3_PROTECT_EDGES_SIZING_FIELD_H
#define CGAL_MESH_3_PROTECT_EDGES_SIZING_FIELD_H

#include <CGAL/license/Mesh_3.h>


#include <CGAL/Mesh_3/config.h>

#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Mesh_3/io_signature.h>
#ifndef CGAL_NO_ASSERTIONS
#  include <boost/math/special_functions/next.hpp> // for float_prior
#endif
#include <boost/function_output_iterator.hpp>
#if ! defined(CGAL_NO_PRECONDITIONS)
#  include <sstream>
#endif
#ifdef CGAL_MESH_3_DUMP_FEATURES_PROTECTION_ITERATIONS
#include <CGAL/IO/File_binary_mesh_3.h>
#endif
#include <CGAL/Has_timestamp.h>

namespace CGAL {
namespace Mesh_3 {
namespace internal {
  const double min_intersection_factor = .4; // (1-alpha)
  const double weight_modifier = .81; //0.9025;//0.81;
  const double distance_divisor = 2.1;
  const int max_nb_vertices_to_reevaluate_size = 10;

  const int refine_balls_max_nb_of_loops = 29;
// for the origins of `refine_balls_max_nb_of_loops`, that dates from the
// very beginning of this file:
//
//     commit e9b3ff3e5730dab319a8cd581e3eb191559c98db
//     Author: Stéphane Tayeb <Stephane.Tayeb@sophia.inria.fr>
//     Date:   Tue Apr 20 14:53:11 2010 +0000
//     
//         Add draft of _with_features classes.
//
// That constant has had different values in the Git history: 9, 99, and now 29.

} // end namespace internal
} // end namespace Mesh_3
} // end namespace CGAL

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
#include <CGAL/squared_distance_3.h>

#include <fstream>
#include <sstream>

namespace CGAL {

namespace Mesh_3 {

template <typename C3t3>
void debug_dump_c3t3(const std::string filename, const C3t3& c3t3)
{
  std::cerr << "Dump current mesh to " << filename << std::endl;
  std::ofstream out(filename.c_str(), 
                    std::ios_base::out|std::ios_base::binary);
  out << "binary CGAL c3t3 " << CGAL::Get_io_signature<C3t3>()() << "\n";
  CGAL::set_binary_mode(out);
  out << c3t3;
}


template <typename C3T3, typename MeshDomain, typename SizingFunction>
class Protect_edges_sizing_field
{
  typedef Protect_edges_sizing_field          Self;
  
public:
  typedef typename C3T3::Triangulation        Tr;
  typedef typename Tr::Bare_point             Bare_point;
  typedef typename Tr::Weighted_point         Weighted_point;
  typedef typename Weighted_point::Weight     Weight;

  typedef typename Tr::Geom_traits            Gt;
  typedef typename Gt::FT                     FT;

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
                             SizingFunction size=SizingFunction(),
                             const FT minimal_size = FT());
  
  void operator()(const bool refine=true);

  void set_nonlinear_growth_of_balls(bool b = true) {
    nonlinear_growth_of_balls = b;
  }
  
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
  
  /// Insert point(p,w) into triangulation and set its dimension to \c dim and
  /// it's index to \c index.
  /// The handle of the newly created vertex is returned.
  Vertex_handle insert_point(const Bare_point& p,
                             const Weight& w,
                             int dim,
                             const Index& index,
                             const bool special_ball = false);

  /**
   * Insert point(p,w) into triangulation and set its dimension to \c dim and
   * it's index to \c index.
   * The handle of the newly created vertex is returned.
   * 
   * This function also ensures that point(p,w) will not be inside a
   * sphere, by decreasing the radius of any sphere that contains it.
   * It also ensures that no point of the triangulation will be inside its
   * sphere, by decreasing w.
   */
  template <typename ErasedVeOutIt>
  std::pair<Vertex_handle, ErasedVeOutIt>
  smart_insert_point(const Bare_point& p,
		     Weight w,
		     int dim,
		     const Index& index,
		     ErasedVeOutIt out);
  
  
  /// Insert balls between points which are pointed by handles \c vp and \c vq
  /// on curve identified by \c curve_index
  template <typename ErasedVeOutIt>
  ErasedVeOutIt insert_balls(const Vertex_handle& vp,
			     const Vertex_handle& vq,
			     const Curve_segment_index& curve_index,
			     ErasedVeOutIt out);

  /**
   * Insert balls
   * Preconditions:
   *  - size_p < size_q
   *  - pq_geodesic > 0
   */
  template <typename ErasedVeOutIt>
  ErasedVeOutIt insert_balls(const Vertex_handle& vp,
			     const Vertex_handle& vq,
			     const FT size_p,
			     const FT size_q,
			     const FT pq_geodesic,
			     const CGAL::Sign distance_sign,
			     const Curve_segment_index& curve_index,
			     ErasedVeOutIt out);
  
  /// Returns true if balls of \c va and \c vb intersect, and (va,vb) is not
  /// an edge of the complex
  bool non_adjacent_but_intersect(const Vertex_handle& va,
                                  const Vertex_handle& vb) const;
  
  /// Returns true if balls of \c va and \c vb intersect
  bool do_balls_intersect(const Vertex_handle& va,
                          const Vertex_handle& vb) const;
  
  /// Change size of the ball of vertex \c v.
  Vertex_handle change_ball_size(const Vertex_handle& v, const FT size,
                                 const bool special_ball = false);
  
  
  /// Returns true if balls of v1 and v2 intersect "enough"
  bool is_sampling_dense_enough(const Vertex_handle& v1,
                                const Vertex_handle& v2) const;
  
  /// Takes an iterator on Vertex_handle as input and check if the sampling
  /// of those vertices is ok. If not, fix it.
  void check_and_repopulate_edges();

  /// Checks if vertex \c v is well sampled, and if its not the case, fix it.
  /// Fills out with deleted vertices during this process. out value type
  /// is Vertex_handle.
  template <typename ErasedVeOutIt>
  ErasedVeOutIt
  check_and_fix_vertex_along_edge(const Vertex_handle& v, ErasedVeOutIt out);
  
  /// Walk along edge from \c start, following the direction \c start to 
  /// \c next, and fills \c out with the vertices which do not fullfill 
  /// the sampling conditions
  template <typename ErasedVeOutIt>
  ErasedVeOutIt
  walk_along_edge(const Vertex_handle& start,
                  const Vertex_handle& next,
                  const bool test_sampling,
                  ErasedVeOutIt out) const;
  
  /// Returns next vertex along edge, i.e vertex after \c start, following
  /// the direction from \c previous to \c start
  /// \pre (previous,start) is in c3t3
  Vertex_handle next_vertex_along_edge(const Vertex_handle& start,
                                       const Vertex_handle& previous) const;
  
  /// Replace vertices between ]begin,last[ by new vertices, along curve
  /// identified by \c curve_index
  /// The value type of InputIterator is Vertex_handle.
  template <typename InputIterator, typename ErasedVeOutIt>
  ErasedVeOutIt repopulate(InputIterator begin,
			   InputIterator last,
			   const Curve_segment_index& index,
			   ErasedVeOutIt out);

  template <typename InputIterator, typename ErasedVeOutIt>
  ErasedVeOutIt
  analyze_and_repopulate(InputIterator begin,
                         InputIterator last,
                         const Curve_segment_index& index,
			 ErasedVeOutIt out);
  
  /// Checks if \c v2 size is compatible (i.e. greater) with the linear
  /// interpolation of sizes of \c v1 and \c v3 
  bool is_sizing_field_correct(const Vertex_handle& v1,
                               const Vertex_handle& v2,
                               const Vertex_handle& v3,
                               const Curve_segment_index& index) const;

  /// Repopulate all incident curve around corner \c v
  /// \pre \c v is a corner of c3t3 
  template <typename ErasedVeOutIt>
  ErasedVeOutIt
  repopulate_edges_around_corner(const Vertex_handle& v, ErasedVeOutIt out);
  
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
    typename C3T3::Triangulation::Geom_traits::Construct_point_3 wp2p =
      c3t3_.triangulation().geom_traits().construct_point_3_object();

    return compute_distance(wp2p(va->point()), wp2p(vb->point()));
  }
  
  /// Compute euclidean distance between bare points \c and \c q
  FT compute_distance(const Bare_point& p, const Bare_point& q) const
  {
    return CGAL::sqrt(c3t3_.triangulation().geom_traits().
                        compute_squared_distance_3_object()(p,q));
  }
  
  /// Returns the radius of the ball of vertex \c v
  FT get_radius(const Vertex_handle& v) const
  {
    return CGAL::sqrt(v->point().weight());
  }

  /// Test if a given vertex is a special protecting ball
  /// A special ball is a protecting ball whose radius is set to the
  /// minimal radius. Such a ball can exceptionally intersect a ball that
  /// is on a different curve. Special balls are marked with special values
  /// of 'in_dimension'.
  bool is_special(const Vertex_handle&v) const
  {
    return v->is_special();
  }

  /// Set to a negative dimension to mark this ball as a special one.
  void set_special(const Vertex_handle&v) const
  {
    v->set_special();
  }

  /// Returns the real dimension of the vertex
  int get_dimension(const Vertex_handle& v) const
  {
    return c3t3_.in_dimension(v);
  }

  /// Query the sizing field and returns its value at the point p, or
  /// minimal_size if the later is greater.
  FT query_size(const Bare_point& p, int dim, const Index& index) const
  {
    if(dim < 0) dim = -1 - dim; // Convert the dimension if it was set to a
                                // negative value (marker for special balls).
    const FT s = size_(p, dim, index);
    // if(minimal_size_ != FT() && s < minimal_size_) 
    //   return minimal_size_;
    // else
    if(s <= FT(0)) {
      std::stringstream msg;
      msg.precision(17);
      msg << "Error: the sizing field is null at ";
      if(dim == 0) msg << "corner (";
      else msg << "point (";
      msg << p << ")";
      CGAL_error_msg(msg.str().c_str());
    }
    return s;
  }

  template <typename Vertex_handle>
  static std::string disp_vert(Vertex_handle v, Tag_true) {
    std::stringstream ss;
    ss << (void*)(&*v) << "[ts=" << v->time_stamp() << "]"
       << "(" << v->point() <<")";
    return ss.str();
  }

  template <typename Vertex_handle>
  static std::string disp_vert(Vertex_handle v, Tag_false) {
    std::stringstream ss;
    ss << (void*)(&*v) << "(" << v->point() <<")";
    return ss.str();
  }

  template <typename Vertex_handle>
  static std::string disp_vert(Vertex_handle v)
  {
    typedef typename std::iterator_traits<Vertex_handle>::value_type Vertex;
    return disp_vert(v, CGAL::internal::Has_timestamp<Vertex>());
  }

private:
  C3T3& c3t3_;
  const MeshDomain& domain_;
  SizingFunction size_;
  FT minimal_size_;
  Weight minimal_weight_;
  std::set<Curve_segment_index> treated_edges_;
  std::set<Vertex_handle> unchecked_vertices_;
  int refine_balls_iteration_nb;
  bool nonlinear_growth_of_balls;
};


template <typename C3T3, typename MD, typename Sf>
Protect_edges_sizing_field<C3T3, MD, Sf>::
Protect_edges_sizing_field(C3T3& c3t3, const MD& domain, 
                           Sf size, const FT minimal_size)
  : c3t3_(c3t3)
  , domain_(domain)
  , size_(size)
  , minimal_size_(minimal_size)
  , minimal_weight_(CGAL::square(minimal_size))
  , refine_balls_iteration_nb(0)
  , nonlinear_growth_of_balls(false)
{
#ifndef CGAL_MESH_3_NO_PROTECTION_NON_LINEAR
  set_nonlinear_growth_of_balls();
#endif
}


template <typename C3T3, typename MD, typename Sf>
void
Protect_edges_sizing_field<C3T3, MD, Sf>::
operator()(const bool refine)
{
#ifdef CGAL_MESH_3_VERBOSE
  std::cerr << "Inserting protection balls..." << std::endl
            << "  refine_balls = " << std::boolalpha << refine << std::endl
            << "  min_balls_radius = " << minimal_size_ << std::endl
            << "  min_balls_weight = " << minimal_weight_ << std::endl;
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
    CGAL_assertion(minimal_size_ > 0 || c3t3_.is_valid());
 }
  
  // debug_dump_c3t3("dump-mesh-after-protect_edges.binary.cgal", c3t3_);

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
    FT w = CGAL::square(query_size(p, 0, p_index));

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
    Vertex_handle v = smart_insert_point(p, w, 0, p_index,
					 CGAL::Emptyset_iterator()).first;
    CGAL_assertion(v != Vertex_handle());

    // As C3t3::add_to_complex modifies the 'in_dimension' of the vertex,
    // we need to backup and re-set the 'is_special' marker after.
    const bool special_ball = is_special(v);
    c3t3_.add_to_complex(v,cit->first);
    if(special_ball) {
      set_special(v);
    }
  }
} //end insert_corners()


template <typename C3T3, typename MD, typename Sf>
typename Protect_edges_sizing_field<C3T3, MD, Sf>::Vertex_handle
Protect_edges_sizing_field<C3T3, MD, Sf>::
insert_point(const Bare_point& p, const Weight& w, int dim, const Index& index,
             const bool special_ball /* = false */)
{
  using CGAL::Mesh_3::internal::weight_modifier;

  if(dim < 0) dim = -1 - dim; // Convert the dimension if it was set to a
                              // negative value (marker for special balls).

  typedef typename Tr::size_type size_type;
  CGAL_USE_TYPE(size_type);
  
  // Insert point
  CGAL_assertion_code(size_type nb_vertices_before = c3t3_.triangulation().number_of_vertices());

  typename Gt::Construct_weighted_point_3 cwp =
    c3t3_.triangulation().geom_traits().construct_weighted_point_3_object();

  const Weighted_point wp = cwp(p,w*weight_modifier);

  typename Tr::Locate_type lt;
  int li, lj;
  const typename Tr::Cell_handle ch = c3t3_.triangulation().locate(wp, lt, li, lj);
  Vertex_handle v = c3t3_.triangulation().insert(wp, lt, ch, li, lj);

  // If point insertion created an hidden ball, fail
  CGAL_assertion ( Vertex_handle() != v );
  CGAL_assertion ( lt == Tr::VERTEX ||
                   c3t3_.triangulation().number_of_vertices() == (nb_vertices_before+1) );

#if CGAL_MESH_3_PROTECTION_DEBUG & 1
  std::cerr << "Insertion of ";
  if(special_ball) std::cerr << "SPECIAL ";
  std::cerr << "protecting ball " << cwp(p,w*weight_modifier);
  switch(dim) {
  case 0:
    std::cerr << " on corner #";
    break;
  case 1:
    std::cerr << " on curve #";
    break;
  default:
    std::cerr << " ERROR dim=" << dim << " index=";
  }
  std::cerr  << CGAL::oformat(index) << std::endl;
  if(v == Vertex_handle()) std::cerr << "  HIDDEN!\n";
  std::cerr << "The weight was " << w << std::endl;
#endif // CGAL_MESH_3_PROTECTION_DEBUG

  c3t3_.set_dimension(v, dim);
  if(special_ball)
    set_special(v);
  c3t3_.set_index(v,index);
  
  unchecked_vertices_.insert(v);
  
  return v;
}

  
template <typename C3T3, typename MD, typename Sf>
template <typename ErasedVeOutIt>
std::pair<typename Protect_edges_sizing_field<C3T3, MD, Sf>::Vertex_handle,
	  ErasedVeOutIt>
Protect_edges_sizing_field<C3T3, MD, Sf>::
smart_insert_point(const Bare_point& p, Weight w, int dim, const Index& index,
		   ErasedVeOutIt out)
{
#if CGAL_MESH_3_PROTECTION_DEBUG & 1
  std::cerr << "smart_insert_point( (" << p 
            << "), w=" << w
            << ", dim=" << dim
            << ", index=" << CGAL::oformat(index) << ")\n";
#endif
  const Tr& tr = c3t3_.triangulation();
  typename Gt::Compute_squared_distance_3 sq_distance =
    tr.geom_traits().compute_squared_distance_3_object();
  typename Gt::Construct_point_3 wp2p =
    tr.geom_traits().construct_point_3_object();
  typename Gt::Construct_weighted_point_3 cwp =
    tr.geom_traits().construct_weighted_point_3_object();

  bool add_handle_to_unchecked = false; /// add or not the new vertex to
                                        /// the set 'unchecked_vertices'
  bool insert_a_special_ball = false; /// will be passed to the function
                                      /// this->insert_point

  const Weighted_point wp0 = cwp(p); // with weight 0, used for locate()

  if ( tr.dimension() > 2 )
  {
    // Check that new point will not be inside a power sphere

    typename Tr::Locate_type lt;
    int li, lj;
    Cell_handle ch = tr.locate(wp0, lt, li, lj);
    Vertex_handle nearest_vh = tr.nearest_power_vertex(p, ch);
    FT sq_d = sq_distance(p, wp2p(nearest_vh->point()));
    while ( nearest_vh->point().weight() > sq_d && ! is_special(nearest_vh) )
    {
      CGAL_assertion( minimal_size_ >0 || sq_d > 0 );

      bool special_ball = false;
      if(minimal_weight_ != Weight() && sq_d < minimal_weight_) {
        sq_d = minimal_weight_;
        w = minimal_weight_;
        special_ball = true;
        insert_a_special_ball = true;
      }

      // Adapt size
      *out++ = nearest_vh;
      Vertex_handle new_vh = change_ball_size(nearest_vh, CGAL::sqrt(sq_d),
                                              special_ball);
      ch = tr.locate(wp0, lt, li, lj, new_vh);

      // Iterate
      nearest_vh = tr.nearest_power_vertex(p, ch);
      sq_d = sq_distance(p, wp2p(nearest_vh->point()));
    }

    if( is_special(nearest_vh) && nearest_vh->point().weight() > sq_d )
    {
      w = minimal_weight_;
      insert_a_special_ball = true;
    }

    // Change w in order to be sure that no existing point will be included
    // in (p,w)
    std::vector<Vertex_handle> vertices_in_conflict_zone;
    { // fill vertices_in_conflict_zone
      std::set<Vertex_handle> vertices_in_conflict_zone_set;
      std::vector<Cell_handle> cells_in_conflicts;
      Weighted_point wp = cwp(p, w);
      tr.find_conflicts(wp, ch,
                        CGAL::Emptyset_iterator(),
                        std::back_inserter(cells_in_conflicts),
                        CGAL::Emptyset_iterator());

      for(typename std::vector<Cell_handle>::const_iterator 
            it = cells_in_conflicts.begin(),
            end = cells_in_conflicts.end(); it != end; ++it) 
      {
        for(int i = 0, d = tr.dimension(); i <= d; ++i) {
          const Vertex_handle v = (*it)->vertex(i);
          if( ! c3t3_.triangulation().is_infinite(v) ) {
            vertices_in_conflict_zone_set.insert(v);
          }
        }
      }
      vertices_in_conflict_zone.insert(vertices_in_conflict_zone.end(),
				       vertices_in_conflict_zone_set.begin(),
				       vertices_in_conflict_zone_set.end());
    }
    FT min_sq_d = w;
#if CGAL_MESH_3_PROTECTION_DEBUG & 1
    typename Tr::Point nearest_point;
#endif
    for(typename std::vector<Vertex_handle>::const_iterator
          it = vertices_in_conflict_zone.begin(),
          end = vertices_in_conflict_zone.end(); it != end ; ++it )
    {
      const FT sq_d = sq_distance(p, wp2p((*it)->point()));
      if(minimal_weight_ != Weight() && sq_d < minimal_weight_) {
        insert_a_special_ball = true;
#if CGAL_MESH_3_PROTECTION_DEBUG & 1
        nearest_point = (*it)->point();
#endif
        min_sq_d = minimal_weight_;
        if( ! is_special(*it) ) {
	  *out++ = *it;
          ch = change_ball_size(*it, minimal_size_, true)->cell(); // special ball
        }
      } else {
        if(sq_d < min_sq_d) {
#if CGAL_MESH_3_PROTECTION_DEBUG & 1
          nearest_point = (*it)->point();
#endif
          min_sq_d = sq_d;
        }
      }
    }

    if ( w > min_sq_d )
    { 
#if CGAL_MESH_3_PROTECTION_DEBUG & 1
      std::cerr << "smart_insert_point: weight " << w
                << " reduced to " << min_sq_d
                << "\n (near existing point: " << nearest_point << " )\n";
#endif
      w = min_sq_d;
      add_handle_to_unchecked = true;
    }

    if(lt != Tr::VERTEX) {
      using CGAL::Mesh_3::internal::weight_modifier;
      CGAL_assertion_code(std::vector<Vertex_handle> hidden_vertices;);
      CGAL_assertion_code(ch = tr.locate(wp0, lt, li, lj, ch););
      CGAL_assertion_code(Weighted_point wpp = cwp(p, w*weight_modifier);)
      CGAL_assertion_code(tr.vertices_inside_conflict_zone(wpp,
                                                           ch,
                                                           std::back_inserter(hidden_vertices)));

      CGAL_assertion(hidden_vertices.empty());
    }
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
        FT sq_d = sq_distance(p, wp2p(it->point()));
        if ( it->point().weight() > sq_d )
        { 
          bool special_ball = false;
          if(minimal_weight_ != Weight() && sq_d > minimal_weight_) {
            sq_d = minimal_weight_;
            w = minimal_weight_;
            special_ball = true;
            insert_a_special_ball = true;
          }
          if( ! is_special(it) ) { 
	    *out++ = it;
            change_ball_size(it, CGAL::sqrt(sq_d), special_ball);
            restart = true;
          }
          break;
        }
      }
    }
    
    FT min_sq_d = w;
    typename Tr::Point nearest_point;
    // Change w in order to be sure that no existing point will be included
    // in (p,w)
    for ( typename Tr::Finite_vertices_iterator it = tr.finite_vertices_begin(),
         end = tr.finite_vertices_end() ; it != end ; ++it )
    {
      FT sq_d = sq_distance(p, wp2p(it->point()));
      if(sq_d < min_sq_d) {
        min_sq_d = sq_d;
        nearest_point = it->point();
      }
    }    
    if ( w > min_sq_d )
    { 
#if CGAL_MESH_3_PROTECTION_DEBUG & 1
      std::cerr << "smart_insert_point: weight " << w
                << " reduced to " << min_sq_d
                << "\n (near existing point: " << nearest_point << " )\n";
#endif
      w = min_sq_d;
      add_handle_to_unchecked = true;
    }
  }

  const FT w_max = CGAL::square(query_size(p, dim, index));

  if(w > w_max) {
#if CGAL_MESH_3_PROTECTION_DEBUG & 1
    std::cerr << "smart_insert_point: weight " << w
              << " reduced to " << w_max << " (sizing field)\n";
#endif
    w = w_max;
    add_handle_to_unchecked = true;
  }
  if( w < minimal_weight_) {
    w = minimal_weight_;
    insert_a_special_ball = true;
  }
  Vertex_handle v = insert_point(p,w,dim,index, insert_a_special_ball);
  /// @TODO `insert_point` does insert in unchecked_vertices anyway!
  if ( add_handle_to_unchecked ) { unchecked_vertices_.insert(v); }
  
  return std::pair<Vertex_handle, ErasedVeOutIt>(v, out);
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
#if CGAL_MESH_3_PROTECTION_DEBUG & 1
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
        typename Gt::Construct_weighted_point_3 cwp =
          c3t3_.triangulation().geom_traits().construct_weighted_point_3_object();

        // Even if the curve is a cycle, it can intersect other curves at
        // its first point (here 'p'). In that case, 'p' is a corner, even
        // if the curve is a cycle.
        if(!c3t3_.triangulation().is_vertex(cwp(p), vp))
        {
          // if 'p' is not a corner, find out a second point 'q' on the
          // curve, "far" from 'p', and limit the radius of the ball of 'p'
          // with the third of the distance from 'p' to 'q'.
          FT p_size = query_size(p, 1, p_index);

          FT curve_lenght = domain_.geodesic_distance(p, p, curve_index);

          Bare_point other_point =
            domain_.construct_point_on_curve_segment(p,
                                                     curve_index,
                                                     curve_lenght / 2);
          p_size = (std::min)(p_size,
                              compute_distance(p, other_point) / 3);
          vp = smart_insert_point(p,
                                  CGAL::square(p_size),
                                  1,
                                  p_index,
				  CGAL::Emptyset_iterator()).first;
        }
        // No 'else' because in that case 'is_vertex(..)' already filled
        // the variable 'vp'.
        vq = vp;
      }
      
      // Insert balls and set treated
      // if(do_balls_intersect(vp, vq)) {
      //   CGAL_assertion(is_special(vp) || is_special(vq));
      // }
      // else
      {
        insert_balls(vp, vq, curve_index, Emptyset_iterator());
      }
      set_treated(curve_index);
    }
    // std::stringstream s;
    // s << "dump-mesh-curve-" << curve_index << ".binary.cgal";
    // debug_dump_c3t3(s.str(), c3t3_);
  }
} //end insert_balls_on_edges()


template <typename C3T3, typename MD, typename Sf>
typename Protect_edges_sizing_field<C3T3, MD, Sf>::Vertex_handle
Protect_edges_sizing_field<C3T3, MD, Sf>::
get_vertex_corner_from_point(const Bare_point& p, const Index&) const
{
  typename Gt::Construct_weighted_point_3 cwp =
    c3t3_.triangulation().geom_traits().construct_weighted_point_3_object();

  // Get vertex_handle associated to corner (dim=0) point
  Vertex_handle v;
  CGAL_assertion_code( bool q_finded = )
  c3t3_.triangulation().is_vertex(cwp(p), v);
  // Let the weight be 0, because is_vertex only locates the point, and
  // check that the location type is VERTEX.
  CGAL_assertion( q_finded );
  return v;
}


template <typename C3T3, typename MD, typename Sf>
template <typename ErasedVeOutIt>
ErasedVeOutIt
Protect_edges_sizing_field<C3T3, MD, Sf>::
insert_balls(const Vertex_handle& vp,
             const Vertex_handle& vq,
             const Curve_segment_index& curve_index,
	     ErasedVeOutIt out)
{
  typename C3T3::Triangulation::Geom_traits::Construct_point_3 wp2p =
    c3t3_.triangulation().geom_traits().construct_point_3_object();

  // Get size of p & q
  const Bare_point& p = wp2p(vp->point());
  const Bare_point& q = wp2p(vq->point());

  const FT sp = get_radius(vp);
  const FT sq = get_radius(vq);

  // Compute geodesic distance
  const FT pq_geo_signed = domain_.geodesic_distance(p, q, curve_index);
  const CGAL::Sign d_sign = CGAL::sign(pq_geo_signed);
  const FT pq_geo = CGAL::abs(pq_geo_signed);

  // Insert balls
  return
    (sp <= sq) ?
    insert_balls(vp, vq, sp, sq, pq_geo, d_sign, curve_index, out) :
    insert_balls(vq, vp, sq, sp, pq_geo, -d_sign, curve_index, out);
}



template <typename C3T3, typename MD, typename Sf>
template <typename ErasedVeOutIt>
ErasedVeOutIt
Protect_edges_sizing_field<C3T3, MD, Sf>::
insert_balls(const Vertex_handle& vp,
             const Vertex_handle& vq,
             const FT sp,
             const FT sq,
             const FT d,
             const CGAL::Sign d_sign,
             const Curve_segment_index& curve_index,
	     ErasedVeOutIt out)
{
  typename C3T3::Triangulation::Geom_traits::Construct_point_3 wp2p =
    c3t3_.triangulation().geom_traits().construct_point_3_object();

#if CGAL_MESH_3_PROTECTION_DEBUG & 1
  std::cerr << "insert_balls(vp=" << disp_vert(vp) << ",\n"
            << "             vq=" << disp_vert(vq) << ",\n"
            << "             sp=" << sp << ",\n"
            << "             sq=" << sq << ",\n"
            << "              d=" << d << ",\n"
            << "             d_sign=" << d_sign << ",\n"
            << "             curve_index=" << curve_index << ")\n";
#endif
  CGAL_precondition(d > 0);
  CGAL_precondition(sp <= sq);
#if ! defined(CGAL_NO_PRECONDITIONS)
  if(sp <= 0) {
    std::stringstream msg;;
    msg.precision(17);
    msg << "Error: the mesh sizing field is null at point (";
    msg << wp2p(vp->point()) << ")!";
    CGAL_precondition_msg(sp > 0, msg.str().c_str());
  }
#endif // ! CGAL_NO_PRECONDITIONS

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
  // if( minimal_weight_ != 0 && n == 0 ) return;

  if(nonlinear_growth_of_balls && refine_balls_iteration_nb < 3)
  {
    // This block tries not to apply the general rule that the size of
    // protecting balls is a linear interpolation of the size of protecting
    // balls at corner. When the curve segment is long enough, pick a point
    // at the middle and choose a new size.
    if(n >= internal::max_nb_vertices_to_reevaluate_size &&
       d >= (internal::max_nb_vertices_to_reevaluate_size * minimal_weight_)) {
#if CGAL_MESH_3_PROTECTION_DEBUG & 1
      std::cerr << "Number of to-be-inserted balls is: "
                << n << "\n  between points ("
                << vp->point() << ") and (" << vq->point()
                << ") (geodesic distance: "
                << domain_.geodesic_distance(wp2p(vp->point()), wp2p(vq->point()),
                                             curve_index)
                << ")\n";
#endif
      const Bare_point new_point =
        domain_.construct_point_on_curve_segment(wp2p(vp->point()),
                                                 curve_index,
                                                 d_sign * d / 2);
      const int dim = 1; // new_point is on edge
      const Index index = domain_.index_from_curve_segment_index(curve_index);
      const FT point_weight = CGAL::square(size_(new_point, dim, index));
#if CGAL_MESH_3_PROTECTION_DEBUG & 1
      std::cerr << "  middle point: " << new_point << std::endl;
      std::cerr << "  new weight: " << point_weight << std::endl;
#endif
      std::pair<Vertex_handle, ErasedVeOutIt> pair =
        smart_insert_point(new_point,
                           point_weight,
                           dim,
                           index,
                           out);
      const Vertex_handle new_vertex = pair.first;
      out = pair.second;
      const FT sn = get_radius(new_vertex);
      if(sp <= sn) {
        out=insert_balls(vp, new_vertex, sp, sn, d/2, d_sign, curve_index, out);
      } else {
        out=insert_balls(new_vertex, vp, sn, sp, d/2, -d_sign, curve_index, out);
      }
      if(sn <= sq) {
        out=insert_balls(new_vertex, vq, sn, sq, d/2, d_sign, curve_index, out);
      } else {
        out=insert_balls(vq, new_vertex, sq, sn, d/2, -d_sign, curve_index, out);
      }
      return out;
    }
  } // nonlinear_growth_of_balls

  FT r = (sq - sp) / FT(n+1);

#if CGAL_MESH_3_PROTECTION_DEBUG & 1
  std::cerr << "  n=" << n
            << "\n  r=" << r << std::endl;
#endif


  // Adjust size of steps, D = covered distance
  FT D = sp*FT(n+1) + FT((n+1)*(n+2)) / FT(2) * r ;
  
  FT dleft_frac = d / D;

  // Initialize step sizes
  FT step_size = sp + r;
  FT norm_step_size = dleft_frac * step_size;
  
  // Initial distance
  FT d_signF = static_cast<FT>(d_sign);
  FT pt_dist = d_signF * norm_step_size;
  Vertex_handle prev = vp;
  const Bare_point& p = wp2p(vp->point());

  // if ( (0 == n) && 
  //      ( (d >= sp+sq) || !is_sampling_dense_enough(vp, vq) ) )

  // If there is some place to insert one point, insert it
  if ( (0 == n) && (d >= sp+sq) )
  {
    n = 1;
    step_size = sp + (d-sp-sq) / FT(2);
    pt_dist = d_signF * step_size;
    norm_step_size = step_size;
  } else if(vp == vq && n == 1) {
    // In case we sample a full cycle, we want to add at least 2
    // balls, equally distributed.
    n = 2;
    step_size = d / FT(n+1);
    pt_dist = d_signF * step_size;
    norm_step_size = step_size;
  } else {
    CGAL_assertion_code(using boost::math::float_prior);
    CGAL_assertion(n==0 ||
                   dleft_frac >= float_prior(float_prior(1.)));
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
    std::pair<Vertex_handle, ErasedVeOutIt> pair =
      smart_insert_point(new_point, point_weight, dim, index, out);
    Vertex_handle new_vertex = pair.first;
    out = pair.second;
    
    // Add edge to c3t3
    if(!c3t3_.is_in_complex(prev, new_vertex)) {
      c3t3_.add_to_complex(prev, new_vertex, curve_index);
    }
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
    if(!c3t3_.is_in_complex(prev, vq)) {
      c3t3_.add_to_complex(prev, vq, curve_index);
    }
  }
  return out;
}
  
  
template <typename C3T3, typename MD, typename Sf>
void
Protect_edges_sizing_field<C3T3, MD, Sf>::
refine_balls()
{
#if CGAL_MESH_3_PROTECTION_DEBUG & 4
  dump_c3t3(c3t3_, "dump-before-refine_balls");
  dump_c3t3_edges(c3t3_, "dump-before-refine_balls");
#endif
  Triangulation& tr = c3t3_.triangulation();
  
  // Loop
  bool restart = true;
  using CGAL::Mesh_3::internal::refine_balls_max_nb_of_loops;
  this->refine_balls_iteration_nb = 0;
  while ( (!unchecked_vertices_.empty() || restart) &&
          this->refine_balls_iteration_nb < refine_balls_max_nb_of_loops)
  {
#ifdef CGAL_MESH_3_DUMP_FEATURES_PROTECTION_ITERATIONS
    std::ostringstream oss;
    oss << "dump_protecting_balls_" << refine_balls_iteration_nb << ".cgal";
    std::ofstream outfile(oss.str(), std::ios_base::binary | std::ios_base::out);
    CGAL::Mesh_3::save_binary_file(outfile, c3t3_, true);
#endif //CGAL_MESH_3_DUMP_FEATURES_PROTECTION_ITERATIONS

#if CGAL_MESH_3_PROTECTION_DEBUG & 1
    std::cerr << "RESTART REFINE LOOP (" << refine_balls_iteration_nb << ")\n"
              << "\t unchecked_vertices size: " << unchecked_vertices_.size() <<"\n";
#endif
    ++refine_balls_iteration_nb;
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
        using CGAL::Mesh_3::internal::distance_divisor;

        // Compute correct size of balls
        const FT ab = compute_distance(va,vb);
	/// @TOTO pb: get_radius(va) is called several times
        FT sa_new = (std::min)(ab/distance_divisor, get_radius(va));
        FT sb_new = (std::min)(ab/distance_divisor, get_radius(vb));
        
        // In case of va or vb have already been in conflict, keep minimal size
        if ( new_sizes.find(va) != new_sizes.end() )
        { sa_new = (std::min)(sa_new, new_sizes[va]); }
        
        if ( new_sizes.find(vb) != new_sizes.end() )
        { sb_new = (std::min)(sb_new, new_sizes[vb]); }
        
        // Store new_sizes for va and vb
        if ( sa_new != get_radius(va) )
        { new_sizes[va] = sa_new; }
        
        if ( sb_new != get_radius(vb) )
        { new_sizes[vb] = sb_new; }
      }
    }

    // The std::map with Vertex_handle as the key is not robust, because
    // the time stamp of vertices can change during the following loop. The
    // solution is to copy it in a vector.
    std::vector<std::pair<Vertex_handle, FT> >
      new_sizes_copy(new_sizes.begin(), new_sizes.end());
    new_sizes.clear();
    
    // Update size of balls
    for ( typename std::vector<std::pair<Vertex_handle,FT> >::iterator
	    it = new_sizes_copy.begin(),
	    end = new_sizes_copy.end();
	  it != end ; ++it )
    {
      const Vertex_handle v = it->first;
      const FT new_size = it->second;
      // Set size of the ball to new value
      if(minimal_size_ != FT() && new_size < minimal_size_) {
        if(!is_special(v)) {
          change_ball_size(v,minimal_size_,true); // special ball
                  
          // Loop will have to be run again
          restart = true;
        }
      } else {
        change_ball_size(v,new_size);
                
        // Loop will have to be run again
        restart = true;
      }
    }
    
#if CGAL_MESH_3_PROTECTION_DEBUG & 4
    dump_c3t3(c3t3_, "dump-before-check_and_repopulate_edges");
    dump_c3t3_edges(c3t3_, "dump-before-check_and_repopulate_edges");
#endif
    // Check edges
    check_and_repopulate_edges();
  }

  if(this->refine_balls_iteration_nb == refine_balls_max_nb_of_loops)
    std::cerr << "Warning : features protection has reached maximal "
              << " number of loops." << std::endl
              << "          It might result in a crash." << std::endl;

} // end refine_balls()


template <typename C3T3, typename MD, typename Sf>
bool
Protect_edges_sizing_field<C3T3, MD, Sf>::
non_adjacent_but_intersect(const Vertex_handle& va, const Vertex_handle& vb) const
{
  if ( ! c3t3_.is_in_complex(va,vb) )
  {
    return do_balls_intersect(va, vb);
  }
  
  return false;
}


template <typename C3T3, typename MD, typename Sf>
bool
Protect_edges_sizing_field<C3T3, MD, Sf>::
do_balls_intersect(const Vertex_handle& va, const Vertex_handle& vb) const
{
  typename Gt::Construct_sphere_3 sphere = 
    c3t3_.triangulation().geom_traits().construct_sphere_3_object();
    
  typename Gt::Do_intersect_3 do_intersect = 
    c3t3_.triangulation().geom_traits().do_intersect_3_object();
  typename C3T3::Triangulation::Geom_traits::Construct_point_3 wp2p =
    c3t3_.triangulation().geom_traits().construct_point_3_object();

  const Bare_point& a = wp2p(va->point());
  const Bare_point& b = wp2p(vb->point());

  const FT& sq_ra = va->point().weight();
  const FT& sq_rb = vb->point().weight();
    
  return do_intersect(sphere(a, sq_ra), sphere(b, sq_rb));
}

template <typename C3T3, typename MD, typename Sf>
typename Protect_edges_sizing_field<C3T3, MD, Sf>::Vertex_handle
Protect_edges_sizing_field<C3T3, MD, Sf>::
change_ball_size(const Vertex_handle& v, const FT size, const bool special_ball)
{
  // Check if there is something to do
  Weight w = CGAL::square(size);
  // if(v->point().weight() == w)
  // { return v; }

#if CGAL_MESH_3_PROTECTION_DEBUG & 1
  std::cerr << "change_ball_size(v=" << disp_vert(v)
            << " dim=" << c3t3_.in_dimension(v)
            << " index=" << CGAL::oformat(c3t3_.index(v))
            << " ,\n"
            << "                 size=" << size 
            << ", special_ball=" << std::boolalpha << special_ball << std::endl;
#endif
  
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
  typename Gt::Construct_point_3 wp2p =
    c3t3_.triangulation().geom_traits().construct_point_3_object();

  Index index = c3t3_.index(v);
  int dim = get_dimension(v);
  Bare_point p = wp2p(v->point());

  // Remove v from corners
  boost::optional<Corner_index> corner_index =
    boost::make_optional(false, Corner_index());
  if ( c3t3_.is_in_complex(v) )
  {
    corner_index = c3t3_.corner_index(v);
    c3t3_.remove_from_complex(v);
  }

  unchecked_vertices_.erase(v);
  // Change v size
  c3t3_.triangulation().remove(v);

  CGAL_assertion_code(typename Gt::Construct_weighted_point_3 cwp =
                        c3t3_.triangulation().geom_traits().construct_weighted_point_3_object();)
  CGAL_assertion_code(const Weighted_point wp = cwp(p,w);)
  CGAL_assertion_code(Tr& tr = c3t3_.triangulation());
  CGAL_assertion_code(Cell_handle ch = tr.locate(wp));
  CGAL_assertion_code(std::vector<Vertex_handle> hidden_vertices);
  CGAL_assertion_code(if(tr.dimension() > 2)
                      tr.vertices_inside_conflict_zone(wp,
                                                       ch,
                                                       std::back_inserter(hidden_vertices)));

  Vertex_handle new_v = insert_point(p, w , dim, index, special_ball);
  CGAL_assertion(hidden_vertices.empty());

  CGAL_assertion( (! special_ball) || is_special(new_v) ); 
  
  // TODO: ensure that this condition is always satisfied (Pedro's code ?)
  CGAL_assertion(v==new_v);
  //new_v->set_meshing_info(size*size);

  // Restore v in corners
  if ( corner_index )
  {
    // As C3t3::add_to_complex modifies the 'in_dimension' of the vertex,
    // we need to backup and re-set the 'is_special' marker after.
    const bool special_ball = is_special(new_v);
    c3t3_.add_to_complex(new_v,*corner_index);
    if(special_ball) {
        set_special(new_v);
    }
  }
  
  // Restore c3t3 edges
  for ( typename Incident_vertices::iterator it = incident_vertices.begin(),
       end = incident_vertices.end() ; it != end ; ++it )
  {
    // Restore connectivity in c3t3
    c3t3_.add_to_complex(new_v, it->first, it->second);
  }

  // Update unchecked vertices
  unchecked_vertices_.insert(new_v);
  return new_v;
}
namespace details {

  // Functor used by Protect_edges_sizing_field::check_and_repopulate_edges, below
  template <typename Set>
  class Erase_element_from_set {
    Set* set_ptr;
  public:
    Erase_element_from_set(Set& set) : set_ptr(&set) {}

    void operator()(const typename Set::key_type& x)
    {
      set_ptr->erase(x);
    }
  }; // end class Erase_element_from_set

} // end namespace details

template <typename C3T3, typename MD, typename Sf>
void
Protect_edges_sizing_field<C3T3, MD, Sf>::
check_and_repopulate_edges()
{
#if CGAL_MESH_3_PROTECTION_DEBUG & 1
  std::cerr << "check_and_repopulate_edges()\n";
#endif
  typedef std::set<Vertex_handle> Vertices;
  Vertices vertices;
  std::copy( unchecked_vertices_.begin(), unchecked_vertices_.end(),
             std::inserter(vertices,vertices.begin()) );
  
  unchecked_vertices_.clear();
  
  // Fix edges
  while ( !vertices.empty() )
  {
    Vertex_handle v = *vertices.begin();
    vertices.erase(vertices.begin());

    details::Erase_element_from_set<Vertices> erase_from_vertices(vertices);

    check_and_fix_vertex_along_edge(v,
      boost::make_function_output_iterator(erase_from_vertices));
  }
}

  
template <typename C3T3, typename MD, typename Sf>
template <typename ErasedVeOutIt>
ErasedVeOutIt
Protect_edges_sizing_field<C3T3, MD, Sf>::
check_and_fix_vertex_along_edge(const Vertex_handle& v, ErasedVeOutIt out)
{
#if CGAL_MESH_3_PROTECTION_DEBUG & 1
  std::cerr << "check_and_fix_vertex_along_edge(" 
            << disp_vert(v)
            << " dim=" << get_dimension(v)
            << " index=" << CGAL::oformat(c3t3_.index(v))
            << " special=" << std::boolalpha << is_special(v)
            << ")\n";
#endif
  // If v is a corner, then all incident edges have to be checked
  if ( c3t3_.is_in_complex(v) )
  {
    return repopulate_edges_around_corner(v, out);
  }
  
  // Get incident vertices along c3t3 edge
  Incident_vertices incident_vertices;
  c3t3_.adjacent_vertices_in_complex(v, std::back_inserter(incident_vertices));
  CGAL_assertion(v->is_special()
                 || incident_vertices.size() == 0
                 || incident_vertices.size() == 1
                 || incident_vertices.size() == 2);
  if(incident_vertices.size() == 0) return out;
  // The size of 'incident_vertices' can be 0 if v is a ball that covers
  // entirely a closed curve. 
  // The size can also be 1 if the curve is a cycle, and the temporary
  // mesh is only two balls on the cycle: then each ball has only one
  // neighbor.

  // Walk along edge to find the edge piece which is not correctly sampled
  typedef std::list<Vertex_handle> Vertex_list;
  Vertex_list to_repopulate;
  to_repopulate.push_front(v);
  
  const Vertex_handle& previous = incident_vertices.front().first;
  const Vertex_handle& next = incident_vertices.back().first;

  // Walk following direction (v,previous)
  walk_along_edge(v, previous, true, std::front_inserter(to_repopulate));
#if CGAL_MESH_3_PROTECTION_DEBUG & 1
  std::cerr <<  "to_repopulate.size()=" << to_repopulate.size() << "\n";
#endif // CGAL_MESH_3_PROTECTION_DEBUG  

  // Check whether a complete circle has been discovered or not
  if (   to_repopulate.size() == 1
      || to_repopulate.front() != to_repopulate.back() )
  {
    // Walk in other direction (v,next)
    walk_along_edge(v, next, true, std::back_inserter(to_repopulate));    
#if CGAL_MESH_3_PROTECTION_DEBUG & 1
    std::cerr <<  "to_repopulate.size()=" << to_repopulate.size() << "\n";
#endif // CGAL_MESH_3_PROTECTION_DEBUG  
  }
  
  // If only v is in to_repopulate, there is nothing to do
  if ( to_repopulate.size() == 1 )
  {
    *out++ = *to_repopulate.begin();
    return out;
  }
  
  // Store erased vertices
  // out = std::copy(to_repopulate.begin(), to_repopulate.end(), out);

  // Repopulate edge
  const Curve_segment_index& curve_index = incident_vertices.front().second;
  
  out = analyze_and_repopulate(to_repopulate.begin(),
			       --to_repopulate.end(),
			       curve_index,
			       out);
  
  return out;
}


template <typename C3T3, typename MD, typename Sf>
bool
Protect_edges_sizing_field<C3T3, MD, Sf>::
is_sampling_dense_enough(const Vertex_handle& v1, const Vertex_handle& v2) const
{
  using CGAL::Mesh_3::internal::min_intersection_factor;
  CGAL_precondition(c3t3_.is_in_complex(v1,v2));

  typename C3T3::Triangulation::Geom_traits::Construct_point_3 wp2p =
    c3t3_.triangulation().geom_traits().construct_point_3_object();

  // Get sizes
  FT size_v1 = get_radius(v1);
  FT size_v2 = get_radius(v2);

  Curve_segment_index curve_index = Curve_segment_index();
  if(get_dimension(v1) == 1) {
    curve_index = domain_.curve_segment_index(v1->index());
  } else if(get_dimension(v2) == 1) {
    curve_index = domain_.curve_segment_index(v2->index());
  }

  if(curve_index != Curve_segment_index()) {
    FT geodesic_distance = domain_.geodesic_distance(wp2p(v1->point()),
                                                     wp2p(v2->point()),
                                                     curve_index);
    if(domain_.is_cycle(wp2p(v1->point()), curve_index)) {
      geodesic_distance =
        (std::min)(CGAL::abs(geodesic_distance),
                   CGAL::abs(domain_.geodesic_distance(wp2p(v2->point()),
                                                       wp2p(v1->point()),
                                                       curve_index)));
    } else {
      geodesic_distance = CGAL::abs(geodesic_distance);
    }
    // Sufficient condition so that the curve portion between v1 and v2 is
    // inside the union of the two balls.
    if(geodesic_distance > (size_v1 + size_v2)) {
#if CGAL_MESH_3_PROTECTION_DEBUG & 1
      std::cerr << "Note: on curve #" << curve_index << ", between "
                << disp_vert(v1) << " and " << disp_vert(v2) << ", the "
                << "geodesic distance is " << geodesic_distance << " and the"
                << " sum of radii is " << size_v1 + size_v2 << std::endl;
#endif
      return false;
    }
  }
  
  const FT distance_v1v2 = compute_distance(v1,v2);
  // Ensure size_v1 > size_v2
  if ( size_v1 < size_v2 ) { std::swap(size_v1, size_v2); }
  
  // Check if balls intersect
  return distance_v1v2 < (FT(min_intersection_factor) * size_v2 + size_v1);  
}


template <typename C3T3, typename MD, typename Sf>
template <typename ErasedVeOutIt>
ErasedVeOutIt
Protect_edges_sizing_field<C3T3, MD, Sf>::
walk_along_edge(const Vertex_handle& start, const Vertex_handle& next,
                bool /*test_sampling*/,
                ErasedVeOutIt out) const
{
#if CGAL_MESH_3_PROTECTION_DEBUG & 4
  if(!c3t3_.is_in_complex(start, next)) {
    std::cerr << "ERROR: the edge ( " << start->point() << " , "
              << next->point() << " ) is not in complex!\n";
    dump_c3t3(c3t3_, "dump-bug");
    dump_c3t3_edges(c3t3_, "dump-bug-c3t3");
  }
#endif
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
template <typename InputIterator, typename ErasedVeOutIt>
ErasedVeOutIt
Protect_edges_sizing_field<C3T3, MD, Sf>::
repopulate(InputIterator begin, InputIterator last,
	   const Curve_segment_index& index, ErasedVeOutIt out)
{
#if CGAL_MESH_3_PROTECTION_DEBUG & 1
  std::cerr << "repopulate(begin=" << disp_vert(*begin) << "\n"
            << "            last=" << disp_vert(*last)  << "\n"
            << "                  distance(begin, last)=" 
            << std::distance(begin, last) << ",\n"
            << "           index=" << CGAL::oformat(index) << ")\n";
#endif
  CGAL_assertion( std::distance(begin,last) >= 0 );
  
  // May happen
  if ( begin == last ) { return out; }
  
  // Valid because begin < last
  InputIterator current = begin;
  InputIterator previous = current++;
  
  // Remove edges from c3t3.
  while ( current != last )
  {
    c3t3_.remove_from_complex(*previous++, *current++);
  }
  // Keep a point between current and last
  typename C3T3::Triangulation::Geom_traits::Construct_point_3 wp2p =
    c3t3_.triangulation().geom_traits().construct_point_3_object();
  Bare_point point_through = wp2p((*previous)->point());

  // Remove last edge
  c3t3_.remove_from_complex(*previous, *current);
  
  // Remove vertices (don't remove the first one and the last one)
  current = begin;
  while ( ++current != last )
  {
#if CGAL_MESH_3_PROTECTION_DEBUG & 1
    std::cerr << "Removal of ";
    if(is_special(*current)) std::cerr << "SPECIAL ";
    std::cerr << "protecting ball "
              << (*current)->point();
    switch(get_dimension(*current)) {
    case 0:
      std::cerr << " on corner #";
      break;
    case 1:
      std::cerr << " on curve #";
      break;
    default:
      std::cerr << " ERROR dim=" << get_dimension(*current)  << " index=";
    }
    std::cerr  << CGAL::oformat(c3t3_.index(*current)) << std::endl;
#endif // CGAL_MESH_3_PROTECTION_DEBUG
    *out++ = *current;
    c3t3_.triangulation().remove(*current);
  }
  
  // If edge is a cycle, order the iterators according to the orientation of
  // the cycle
  if (  domain_.is_cycle(wp2p((*begin)->point()), index)
      && domain_.distance_sign_along_cycle(wp2p((*begin)->point()),
                                           point_through,
                                           wp2p((*last)->point()),
                                           index ) != CGAL::POSITIVE )
  {
    std::swap(begin,last);
  }
  
  // Repopulate edge
  return insert_balls(*begin, *last, index, out);
}


template <typename C3T3, typename MD, typename Sf>
template <typename InputIterator, typename ErasedVeOutIt>
ErasedVeOutIt
Protect_edges_sizing_field<C3T3, MD, Sf>::
analyze_and_repopulate(InputIterator begin, InputIterator last,
		       const Curve_segment_index& index, ErasedVeOutIt out)
{
#if CGAL_MESH_3_PROTECTION_DEBUG & 1
  std::cerr << "analyze_and_repopulate(begin=" << disp_vert(*begin) << "\n"
            << "                       last=" << disp_vert(*last) << "\n"
            << "                              distance(begin, last)=" 
            << std::distance(begin, last) << ",\n"
            << "                       index=" << CGAL::oformat(index) << ")\n";
#endif
  CGAL_assertion( std::distance(begin,last) >= 0 );
  
  // May happen
  if ( begin == last ) { return out; }
  if ( std::distance(begin,last) == 1 )
  {
    out = repopulate(begin, last, index, out);
    return out;
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
           && is_sizing_field_correct(*ch_stack.top(),*previous,*current, index) )
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
    out = repopulate(next, current, index, out);
    current = next;
  }
  return out;
}
  
template <typename C3T3, typename MD, typename Sf>
bool
Protect_edges_sizing_field<C3T3, MD, Sf>::
is_sizing_field_correct(const Vertex_handle& v1,
                        const Vertex_handle& v2,
                        const Vertex_handle& v3,
                        const Curve_segment_index& curve_index) const
{
  typename C3T3::Triangulation::Geom_traits::Construct_point_3 wp2p =
    c3t3_.triangulation().geom_traits().construct_point_3_object();

  FT s1 = get_radius(v1);
  FT s2 = get_radius(v2);
  FT s3 = get_radius(v3);
  FT D = CGAL::abs(domain_.geodesic_distance(wp2p(v1->point()), wp2p(v3->point()), curve_index));
  FT d = CGAL::abs(domain_.geodesic_distance(wp2p(v1->point()), wp2p(v2->point()), curve_index));
  if(domain_.is_cycle(v1->point().point(), curve_index)) {
    D = (std::min)(D, CGAL::abs(domain_.geodesic_distance(wp2p(v3->point()), wp2p(v1->point()), curve_index)));
    d = (std::min)(d, CGAL::abs(domain_.geodesic_distance(wp2p(v2->point()), wp2p(v1->point()), curve_index)));
  }
  
  return ( s2 >= (s1 + d/D*(s3-s1)) );
}

  
template <typename C3T3, typename MD, typename Sf>
template <typename ErasedVeOutIt>
ErasedVeOutIt
Protect_edges_sizing_field<C3T3, MD, Sf>::
repopulate_edges_around_corner(const Vertex_handle& v, ErasedVeOutIt out)
{
  CGAL_precondition(c3t3_.is_in_complex(v));

  typename C3T3::Triangulation::Geom_traits::Construct_point_3 wp2p =
    c3t3_.triangulation().geom_traits().construct_point_3_object();

  Incident_vertices incident_vertices;
  c3t3_.adjacent_vertices_in_complex(v, std::back_inserter(incident_vertices));
  
  for ( typename Incident_vertices::iterator vit = incident_vertices.begin(),
       vend = incident_vertices.end() ; vit != vend ; ++vit )
  {
    const Vertex_handle& next = vit->first;
    const Curve_segment_index& curve_index = vit->second;
    
    // if `v` is incident to a cycle, it might be that the full cycle,
    // including the edge `[next, v]`, has already been processed by
    // `analyze_and_repopulate()` walking in the other direction.
    if(domain_.is_cycle(wp2p(v->point()), curve_index) &&
       !c3t3_.is_in_complex(v, next)) continue;

    // Walk along each incident edge of the corner
    Vertex_vector to_repopulate;
    to_repopulate.push_back(v);
    walk_along_edge(v, next, true, std::back_inserter(to_repopulate));

    // Return erased vertices
    std::copy(to_repopulate.begin(), to_repopulate.end(), out);

    // Repopulate
    out = analyze_and_repopulate(to_repopulate.begin(), --to_repopulate.end(),
				 curve_index, out);
  }
  
  return out;
}
  
} // end namespace Mesh_3


} //namespace CGAL

#endif // CGAL_MESH_3_PROTECT_EDGES_SIZING_FIELD_H
