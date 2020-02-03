// Copyright (c) 2009-2010 INRIA Sophia-Antipolis (France).
// Copyright (c) 2010-2017 GeometryFactory Sarl (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Stephane Tayeb, Laurent Rineau, Mael Rouxel-Labbé
//
//******************************************************************************
// File Description :
//******************************************************************************

#ifndef CGAL_PERIODIC_3_MESH_3_PROTECT_EDGES_SIZING_FIELD_H
#define CGAL_PERIODIC_3_MESH_3_PROTECT_EDGES_SIZING_FIELD_H

// Main differences with Mesh_3's version:
// - A map [Vertex_handle] --> [position in full space] because:
//   * This position might not be (in full space coordinates) in the canonical instance
//   * Different curves might represent the same corner at different positions
// - change_ball_size()'s signature is different: it can return false if we failed
//   to change the ball size (if the triangulation cover changes during the process).
// - Dummy point processing to avoid overrefinement when dummy points are close
//   to sharp features.
// - Various minor improvements over the original code, which should be brought back
//   (unordered sets, improving 'unchecked_vertices', etc.)
// - No wonky indentation :)

#include <CGAL/license/Periodic_3_mesh_3.h>

#include <CGAL/Mesh_3/io_signature.h>
#ifdef CGAL_MESH_3_DUMP_FEATURES_PROTECTION_ITERATIONS
#include <CGAL/IO/File_binary_mesh_3.h>
#endif
#include <CGAL/Mesh_3/Protect_edges_sizing_field.h>
#include <CGAL/Mesh_3/utilities.h>
#include <CGAL/Mesh_3/Triangulation_helpers.h>

#include <CGAL/enum.h>
#include <CGAL/internal/Has_member_visited.h>
#include <CGAL/iterator.h>
#include <CGAL/number_utils.h>
#include <CGAL/Periodic_3_Delaunay_triangulation_traits_3.h>
#include <CGAL/Periodic_3_Delaunay_triangulation_3.h>
#include <CGAL/Time_stamper.h>

#include <CGAL/boost/iterator/transform_iterator.hpp>

#include <boost/bind.hpp>
#include <boost/function_output_iterator.hpp>
#ifndef CGAL_NO_ASSERTIONS
#  include <boost/math/special_functions/next.hpp> // for float_prior
#endif
#include <boost/optional.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iterator>
#include <list>
#include <sstream>
#include <set>
#include <stack>
#include <string>
#include <utility>
#include <vector>

namespace CGAL {
namespace Periodic_3_mesh_3 {

template <typename C3T3, typename MeshDomain, typename SizingFunction>
class Protect_edges_sizing_field
  : public CGAL::Mesh_3::internal::Debug_messages_tools
{
  typedef Protect_edges_sizing_field          Self;

public:
  typedef typename C3T3::Triangulation        Tr;
  typedef typename Tr::Bare_point             Bare_point;
  typedef typename Tr::Weighted_point         Weighted_point;
  typedef typename Weighted_point::Weight     Weight;

  typedef typename Tr::Geom_traits            Gt;
  typedef typename Gt::FT                     FT;
  typedef typename Gt::Vector_3               Vector_3;

  typedef typename Tr::Vertex                 Vertex;
  typedef typename C3T3::Vertex_handle        Vertex_handle;
  typedef typename C3T3::Edge                 Edge;
  typedef typename C3T3::Cell_handle          Cell_handle;

  typedef typename MeshDomain::Curve_index          Curve_index;
  typedef typename MeshDomain::Corner_index         Corner_index;
  typedef typename MeshDomain::Index                Index;

private:
  typedef typename CGAL::Kernel_traits<MeshDomain>::Kernel   Kernel;

  typedef Periodic_3_Delaunay_triangulation_traits_3<Kernel> DtGt;
  typedef Periodic_3_Delaunay_triangulation_3<DtGt>          Dt;
  typedef Mesh_3::Triangulation_helpers<Dt>                  Dt_helpers;
  typedef typename Dt::Vertex_handle                         Dt_Vertex_handle;

  typedef Mesh_3::Triangulation_helpers<Tr>                  Tr_helpers;

  typedef std::vector<std::pair<Curve_index,Bare_point> >    Incident_edges;
  typedef std::vector<Vertex_handle>                         Vertex_vector;
  typedef std::vector<std::pair<Vertex_handle,Curve_index> > Adjacent_vertices;

  typedef CGAL::Hash_handles_with_or_without_timestamps      Hash_fct;
  typedef boost::unordered_set<Vertex_handle, Hash_fct>      Vertex_set;

private:
  // below could be a set of point pointers
  typedef std::set<Bare_point>                               Domain_pts_container;
  typedef boost::unordered_map<Vertex_handle,
                               Domain_pts_container,
                               Hash_fct>                     Curve_correspondence_map;
  typedef std::map<Curve_index, Curve_correspondence_map>    Correspondence_map;

  typedef typename Curve_correspondence_map::iterator        CCMit;
  typedef typename Curve_correspondence_map::const_iterator  CCMcit;

public:
  Protect_edges_sizing_field(C3T3& c3t3,
                             const MeshDomain& domain,
                             SizingFunction size=SizingFunction(),
                             const FT minimal_size = FT());

  void set_nonlinear_growth_of_balls(bool b = true) { nonlinear_growth_of_balls = b; }

  void operator()(const bool refine = true);

private:
  /// Dump the dummy vertices to a file.
  void dump_dummy_points(const std::string filename) const;

  /// Attempt to remove a vertex that corresponds to a dummy point from the triangulation.
  /// Failure happens when the removal of that vertex would cause the triangulation
  /// to switch to 27-sheets.
  bool try_to_remove_dummy_vertex(const Vertex_handle dummy_vertex) const;

  /// Given a maximal weight 'intended_weight', compute the maximal squared size (weight)
  /// of the protection ball at 'protection_vertex' such that the ball does not contain
  /// any foreign point.
  FT get_maximum_weight(const Vertex_handle protection_vertex, const FT intended_weight) const;

  /// Attempt to change the ball size of 'v'. The difference with 'change_ball_size()' is
  /// that the dummy point 'dummy_point' is inserted in the triangulation before
  /// size change and removed after. This is done to try to keep the triangulation
  /// in a 1-cover mode.
  bool change_ball_size_with_dummy_scaffholding(Vertex_handle& v, const FT w, bool is_special,
                                                const Weighted_point& dummy_point);

  /// Insert the point 'dummy_point' in the triangulation and sets its index to '0'
  /// to mark it as a dummy point.
  Vertex_handle insert_dummy_point(const Weighted_point& dummy_point);

  /// Try to grow the ball at 'protection_vertex' by removing a close dummy vertex.
  bool try_to_remove_close_dummy_vertex(Vertex_handle& protection_vertex,
                                        const Vertex_handle dummy_vertex,
                                        const FT intended_weight);

  /// Attempt to move a dummy vertex to the prescribed position 'new_position'.
  bool try_to_move_dummy_vertex(const Vertex_handle dummy_vertex,
                                const Weighted_point& new_position);

  /// Attempt to move a dummy vertex as far as possible from the protection vertex,
  /// up a to a maximum distance. If successful, try to grow the protection ball
  /// at 'protection_vertex' to its intended size.
  bool move_dummy_vertex_away_from_corner_with_direction(Vertex_handle& protection_vertex,
                                                         const Vertex_handle dummy_vertex,
                                                         const FT intended_weight,
                                                         const Vector_3& dir,
                                                         const FT min_sq_dist);

  /// Attempt to grow the ball at 'protection_vertex' by moving a close dummy vertex.
  bool try_to_move_close_dummy_vertex(Vertex_handle& protection_vertex,
                                      const Vertex_handle dummy_vertex,
                                      const FT intended_weight);

  /// Attempt to grow the ball at 'protection_vertex' by removing or moving a close dummy vertex.
  bool try_to_solve_close_dummy_point(Vertex_handle& protection_vertex,
                                      Vertex_handle dummy_vertex,
                                      FT intended_weight);

private:
  /// Insert corners of the mesh.
  void insert_corners();

  /// Insert balls on every edge.
  void insert_balls_on_edges();

  /// Refine balls.
  void refine_balls();

  /// Return the vertex which corresponds to the corner located at point `p`.
  Vertex_handle get_vertex_corner_from_point(const Bare_point& p,
                                             const Index& p_index) const;

  /// Insert point(p,w) into the triangulation and set its dimension to \c dim
  /// and its index to \c index.
  /// The handle of the newly created vertex is returned.
  template <typename Curve_index_container>
  Vertex_handle insert_point(const Bare_point& p,
                             const Weight& w,
                             int dim,
                             const Index& index,
                             const Curve_index_container& curve_indices,
                             const bool special_ball);

  /**
   * Insert point(p,w) into the triangulation and set its dimension to \c dim and
   * its index to \c index.
   * The handle of the newly created vertex is returned.
   *
   * This function also ensures that `point(p,w)` will not be inside a
   * sphere by decreasing the radius of any sphere that contains it.
   * It also ensures that no point of the triangulation will be inside its
   * sphere by decreasing w.
   */
  template <typename ErasedVeOutIt, typename Curve_index_container>
  std::pair<Vertex_handle, ErasedVeOutIt>
  smart_insert_point(const Bare_point& p,
                     Weight w,
                     int dim,
                     const Index& index,
                     const Curve_index_container& curve_indices,
                     ErasedVeOutIt out);

  /// Insert balls between the points identified by the handles \c vp and \c vq
  /// on the curve identified by \c curve_index.
  ///
  /// \param orientation Orientation of the curve segment between \c vp and
  ///        \c vq, given the orientation of the curve of index \c curve_index
  template <typename ErasedVeOutIt>
  ErasedVeOutIt insert_balls(const Vertex_handle& vp,
                             const Vertex_handle& vq,
                             const Curve_index& curve_index,
                             const CGAL::Orientation orientation,
                             ErasedVeOutIt out);

  /**
   * Insert balls.
   * \pre size_p < size_q
   * \pre pq_geodesic > 0
   */
  template <typename ErasedVeOutIt>
  ErasedVeOutIt insert_balls(const Vertex_handle& vp,
                             const Vertex_handle& vq,
                             const FT size_p,
                             const FT size_q,
                             const FT pq_length,
                             const CGAL::Orientation orientation,
                             const Curve_index& curve_index,
                             ErasedVeOutIt out);

  /// Return `true` if the balls of \c va and \c vb intersect, and (va,vb) is not
  /// an edge of the complex.
  bool non_adjacent_but_intersect(const Vertex_handle& va,
                                  const Vertex_handle& vb) const;

  /// Return `true` if the balls of \c va and \c vb intersect.
  bool do_balls_intersect(const Vertex_handle& va,
                          const Vertex_handle& vb) const;

  /// Change the size of the ball of the vertex \c v.
  bool change_ball_size(Vertex_handle& v,
                        const FT squared_size,
                        const bool special_ball = false);


  /// Return `true` if balls of v1 and v2 intersect "enough".
  /// \param orientation Orientation of the curve segment between \c v1 and
  ///        \c v2, given the orientation of the curve of index
  ///        \c curve_index
  /// \pre `c3t3.curve_index(v1, v2) == curve_index`
  bool is_sampling_dense_enough(const Vertex_handle& v1,
                                const Vertex_handle& v2,
                                const Curve_index& curve_index,
                                const CGAL::Orientation orientation) const;

  /// Take an iterator on Vertex_handle as input and check if the sampling
  /// of those vertices is ok. If not, fix it.
  void check_and_repopulate_edges();

  /// Check if the vertex \c v is well sampled, and if its not the case, fix it.
  /// Fill `out` with deleted vertices during this process. The value type of `out`
  /// is `Vertex_handle`.
  template <typename ErasedVeOutIt>
  ErasedVeOutIt
  check_and_fix_vertex_along_edge(const Vertex_handle& v, ErasedVeOutIt out);

  /// Given two vertices \c start and \c next inserted on the curve with
  /// index \c curve_index, return `CGAL::POSITIVE` if the curve arc from
  /// \c start to \c next is oriented in the same orientation as the curve
  /// segment with index \c curve_index, or `CGAL::NEGATIVE` otherwise.
  ///
  /// \pre `c3t3.curve_index(v1, v2) == curve_index`
  CGAL::Orientation
  orientation_of_walk(const Vertex_handle& start,
                      const Vertex_handle& next,
                      Curve_index curve_index) const;

  /// Walk along the edge from \c start, following the direction \c start to
  /// \c next, and fills \c out with the vertices which do not fulfill
  /// the sampling conditions.
  ///
  /// \param orientation Orientation of the curve segment between \c v1 and
  ///        \c v2, given the orientation of the curve of index
  ///        \c curve_index
  ///
  /// \pre `c3t3.curve_index(v1, v2) == curve_index`
  template <typename ErasedVeOutIt>
  ErasedVeOutIt
  walk_along_edge(const Vertex_handle& start,
                  const Vertex_handle& next,
                  Curve_index curve_index,
                  const CGAL::Orientation orientation,
                  ErasedVeOutIt out) const;

  /// Return the next vertex along edge, i.e the vertex after \c start, following
  /// the direction from \c previous to \c start.
  /// \pre (previous,start) is in c3t3
  /// \pre `c3t3.curve_index(start, previous) == curve_index`
  Vertex_handle
  next_vertex_along_curve(const Vertex_handle& start,
                          const Vertex_handle& previous,
                          const Curve_index& curve_index) const;

  /// Replace the vertices between ]begin,last[ with new vertices, along the curve
  /// identified by \c curve_index.
  /// The value type of `InputIterator` is `Vertex_handle`.
  ///
  /// \param orientation Orientation of the curve segment between \c begin and
  ///        \c last, given the orientation of the curve of index
  ///        \c curve_index
  ///
  template <typename InputIterator, typename ErasedVeOutIt>
  ErasedVeOutIt repopulate(InputIterator begin,
                           InputIterator last,
                           const Curve_index& curve_index,
                           const CGAL::Orientation orientation,
                           ErasedVeOutIt out);

  template <typename InputIterator, typename ErasedVeOutIt>
  ErasedVeOutIt analyze_and_repopulate(InputIterator begin,
                                       InputIterator last,
                                       const Curve_index& index,
                                       const CGAL::Orientation orientation,
                                       ErasedVeOutIt out);

  /// Check if the size of \c v2 is compatible (i.e. greater) with the linear
  /// interpolation of the sizes of \c v1 and \c v3.
  bool is_sizing_field_correct(const Vertex_handle& v1,
                               const Vertex_handle& v2,
                               const Vertex_handle& v3,
                               const Curve_index& index,
                               const CGAL::Orientation orientation) const;

  /// Repopulate all incident curves around the corner \c v.
  /// \pre \c v is a corner of c3t3
  template <typename ErasedVeOutIt>
  ErasedVeOutIt
  repopulate_edges_around_corner(const Vertex_handle& v, ErasedVeOutIt out);

  /// Return `true` if the edge with index \c curve_index is already treated.
  bool is_treated(const Curve_index& curve_index) const
  {
    return (treated_edges_.find(curve_index) != treated_edges_.end());
  }

  /// Set the edge with index \c curve_index as treated.
  void set_treated(const Curve_index& curve_index)
  {
    treated_edges_.insert(curve_index);
  }

  /// Compute the Euclidean distance between the bare points \c and \c q.
  FT compute_distance(const Bare_point& p, const Bare_point& q) const
  {
    // No need to call min_squared_distance() because 'p' and 'q' have been
    // built to have proper respective positions in space
    return CGAL::sqrt(c3t3_.triangulation().geom_traits().
                        compute_squared_distance_3_object()(p,q));
  }

  FT compute_distance(const Vertex_handle& va, const Vertex_handle& vb,
                      const Curve_index& curve_index, const CGAL::Orientation orientation) const
  {
    Bare_point pa, pb;
    boost::tie(pa, pb) = get_positions(va, vb, curve_index, orientation);
    return compute_distance(pa, pb);
  }

  /// Compute the Euclidean distance between the bare points of \c va and \c vb.
  FT compute_distance(const Vertex_handle& va, const Vertex_handle& vb) const
  {
    typename Gt::Construct_point_3 cp =
      c3t3_.triangulation().geom_traits().construct_point_3_object();

    const Weighted_point& wpa = c3t3_.triangulation().point(va);
    const Weighted_point& wpb = c3t3_.triangulation().point(vb);

    return CGAL::sqrt(c3t3_.triangulation().min_squared_distance(cp(wpa), cp(wpb)));
  }

  /// Return the radius of the ball of vertex \c v.
  FT get_radius(const Vertex_handle& v) const
  {
    typename Gt::Compute_weight_3 cw =
      c3t3_.triangulation().geom_traits().compute_weight_3_object();

    const Weighted_point& vwp = c3t3_.triangulation().point(v);
    return CGAL::sqrt(cw(vwp));
  }

  /// Test if a given vertex is a special protecting ball.
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

  /// Return the real dimension of the vertex `v`.
  int get_dimension(const Vertex_handle& v) const
  {
    return c3t3_.in_dimension(v);
  }

  /// Query the sizing field and return its value at the point `p`, or
  /// `minimal_size` if the latter is greater.
  FT query_size(const Bare_point& p, int dim, const Index& index) const
  {
    // Convert the dimension if it was set to a
    // negative value (marker for special balls).
    if(dim < 0)
      dim = -1 - dim;

    const FT s = size_(p, dim, index);
//     if(minimal_size_ != FT() && s < minimal_size_)
//       return minimal_size_;
//     else
    if(s <= FT(0))
    {
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

  /// Output the correspondence map (see definition near class member).
  void dump_correspondence_map() const;

  template<typename Curve_index_container>
  void insert_in_correspondence_map(const Vertex_handle v,
                                    const Bare_point& p,
                                    const Curve_index_container& curve_indices,
                                    typename Curve_index_container::const_iterator* /*sfinae*/ = nullptr);

  void insert_in_correspondence_map(const Vertex_handle v,
                                    const Bare_point& p,
                                    const Curve_index& curve_index);

  void remove_from_correspondence_map(const Vertex_handle v,
                                      const Curve_index& curve_index);

  Bare_point get_position(const Vertex_handle v1,
                          const Curve_index& curve_index) const;

  std::pair<Bare_point, Bare_point>
  get_positions_with_vertex_at_extremity(const Bare_point& known_point,
                                         const Domain_pts_container& extremity_points,
                                         const bool inverted_return_order,
                                         const Curve_index& curve_index,
                                         const CGAL::Orientation orientation) const;

  std::pair<Bare_point, Bare_point> get_positions(const Vertex_handle v1,
                                                  const Vertex_handle v2,
                                                  const Curve_index& curve_index,
                                                  const CGAL::Orientation orientation) const;

  std::pair<Bare_point, Bare_point> get_positions_with_unknown_orientation(const Vertex_handle v1,
                                                                           const Vertex_handle v2,
                                                                           const Curve_index& curve_index) const;

  boost::tuple<Bare_point, Bare_point, Bare_point> get_positions(const Vertex_handle v1,
                                                                 const Vertex_handle v2,
                                                                 const Vertex_handle v3,
                                                                 const Curve_index& curve_index,
                                                                 const CGAL::Orientation orientation) const;

private:
  C3T3& c3t3_;
  const MeshDomain& domain_;
  SizingFunction size_;
  FT minimal_size_;
  Weight minimal_weight_;
  std::set<Curve_index> treated_edges_;
  Vertex_set unchecked_vertices_;
  int refine_balls_iteration_nb;
  bool nonlinear_growth_of_balls;

  // The correspondence map is a map that keeps track of the position of the points
  // in the domain class, which might differ due to periodicity. For example,
  // someone might be using the unit cube as fundamental domain and add the vertex(1,0,0)
  // as corner. If 'v' is the vertex handle of that corner, we have v->point() = (0,0,0).
  //
  // This map must be kept perfectly up to date: removing a vertex from the triangulation
  // means that the correspondence must be removed in the map too.
  Correspondence_map vertex_handle_to_full_space_position;
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
  , vertex_handle_to_full_space_position()
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
  // This class is only meant to be used with periodic triangulations
  CGAL_assertion((boost::is_same<typename Tr::Periodic_tag, CGAL::Tag_true>::value));

#ifdef CGAL_MESH_3_VERBOSE
  std::cerr << "Inserting protection balls..." << std::endl
            << "  refine_balls = " << std::boolalpha << refine << std::endl
            << "  min_balls_radius = " << minimal_size_ << std::endl
            << "  min_balls_weight = " << minimal_weight_ << std::endl;
#endif

#if defined(CGAL_PERIODIC_PROTECTION_ATTEMPT_TO_REMOVE_DUMMY_PTS) && defined(CGAL_MESH_3_VERBOSE)
  std::cerr << "Dumping dummy vertices pre protection treatment" << std::endl;
  dump_dummy_points("dummies_pre_treatement.xyz");
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
  if(refine)
  {
    refine_balls();
#ifdef CGAL_MESH_3_VERBOSE
    std::cerr << "refine_balls() done. Nb of points in triangulation: "
              << c3t3_.triangulation().number_of_vertices() << std::endl;
#endif

    // The assertion below is disabled because c3t3.is_valid() may return false
    // positives for periodic triangulations, specifically because of the line:
    // ' if(! do_intersect(sphere(p, sq_rp), sphere(q, sq_rq))) ... '
    // which does not take periodicity into account
//    CGAL_assertion(minimal_size_ > 0 || c3t3_.is_valid());
 }

  // debug_dump_c3t3("dump-mesh-after-protect_edges.binary.cgal", c3t3_);

#if defined(CGAL_PERIODIC_PROTECTION_ATTEMPT_TO_REMOVE_DUMMY_PTS) && defined(CGAL_MESH_3_VERBOSE)
  std::cerr << "Dumping dummy vertices post protection treatment" << std::endl;
  dump_dummy_points("dummies_post_treatement.xyz");
#endif

  CGAL_expensive_postcondition(c3t3_.triangulation().tds().is_valid());
  CGAL_expensive_postcondition(c3t3_.triangulation().is_valid());
  CGAL_expensive_postcondition(c3t3_.is_valid());
}

template <typename C3T3, typename MD, typename Sf>
void
Protect_edges_sizing_field<C3T3, MD, Sf>::
dump_correspondence_map() const
{
  typename Correspondence_map::const_iterator cmcit = vertex_handle_to_full_space_position.begin(),
                                              cmcend = vertex_handle_to_full_space_position.end();
  for(; cmcit!= cmcend; ++cmcit)
  {
    std::cerr << "Curve " << cmcit->first << std::endl;
    const Curve_correspondence_map& cmap = cmcit->second;
    typename Curve_correspondence_map::const_iterator ccmcit = cmap.begin(),
                                                      ccmcend = cmap.end();
    for(; ccmcit!=ccmcend; ++ccmcit)
    {
      std::cerr << "  Vertex: " << &*(ccmcit->first)
                << " canonical pos: " << c3t3_.triangulation().point(ccmcit->first) << std::endl;
      const Domain_pts_container& vpts = ccmcit->second;
      typename Domain_pts_container::const_iterator ptit = vpts.begin();
      for(; ptit!=vpts.end(); ++ptit)
        std::cerr << "    Point: " << *ptit << std::endl;
    }
  }
}


template <typename C3T3, typename MD, typename Sf>
template<typename Curve_index_container>
void
Protect_edges_sizing_field<C3T3, MD, Sf>::
insert_in_correspondence_map(const Vertex_handle v,
                             const Bare_point& p,
                             const Curve_index_container& curve_indices,
                             typename Curve_index_container::const_iterator* /*sfinae*/)
{
  CGAL_precondition(c3t3_.triangulation().geom_traits().construct_point_3_object()(
                      c3t3_.triangulation().point(v)) ==
                        c3t3_.triangulation().canonicalize_point(p));

  typename Curve_index_container::const_iterator cicit = curve_indices.begin(),
                                                 ciend = curve_indices.end();
  for(; cicit!=ciend; ++cicit)
    insert_in_correspondence_map(v, p, *cicit);
}


template <typename C3T3, typename MD, typename Sf>
void
Protect_edges_sizing_field<C3T3, MD, Sf>::
insert_in_correspondence_map(const Vertex_handle v,
                             const Bare_point& p,
                             const Curve_index& curve_index)
{
#if CGAL_MESH_3_PROTECTION_DEBUG & 2
  std::cerr << "insert_in_correspondence_map: " << &*v
            << " pos: " << p
            << " index: " << curve_index << std::endl;
#endif

  // Get the correspondence map of the correct index
  Curve_correspondence_map& cmap =
    vertex_handle_to_full_space_position[curve_index];

  // Try to insert
  Domain_pts_container bps;
  bps.insert(p);
  std::pair<CCMit, bool> is_insert_successful = cmap.insert(std::make_pair(v, bps));

  if(!is_insert_successful.second) // already existed in the cmap
  {
    Domain_pts_container& v_bps = is_insert_successful.first->second;

#if CGAL_MESH_3_PROTECTION_DEBUG & 4
    std::cerr << "Point already existed for this curve, "
              << "at position: " << *(v_bps.begin()) << std::endl;
#endif

    CGAL_assertion(v_bps.size() == 1);
    v_bps.insert(p);
  }

  // points of a loop can only have one full-space position
  CGAL_assertion_code(const Domain_pts_container& v_bps = is_insert_successful.first->second;)
  CGAL_assertion(v_bps.size() == 1 ||
                 (!domain_.is_loop(curve_index) && v_bps.size() == 2));
}

template <typename C3T3, typename MD, typename Sf>
void
Protect_edges_sizing_field<C3T3, MD, Sf>::
remove_from_correspondence_map(const Vertex_handle v,
                               const Curve_index& curve_index)
{
  typename Correspondence_map::iterator cmap_it =
    vertex_handle_to_full_space_position.find(curve_index);
  CGAL_precondition(cmap_it != vertex_handle_to_full_space_position.end());

  Curve_correspondence_map& cmap = cmap_it->second;
  CCMcit v_in_cmap = cmap.find(v);
  CGAL_precondition(v_in_cmap != cmap.end());

  // Only corners can have 2 different full-space positions and we should never be removing corners.
  CGAL_assertion_code(const Domain_pts_container& v_bps = v_in_cmap->second;)
  CGAL_assertion(v_bps.size() == 1);

  cmap.erase(v_in_cmap);
}


template <typename C3T3, typename MD, typename Sf>
typename Protect_edges_sizing_field<C3T3, MD, Sf>::Bare_point
Protect_edges_sizing_field<C3T3, MD, Sf>::
get_position(const Vertex_handle v1,
             const Curve_index& curve_index) const
{
  typename Correspondence_map::const_iterator cmap_it =
    vertex_handle_to_full_space_position.find(curve_index);
  CGAL_precondition(cmap_it != vertex_handle_to_full_space_position.end());
  const Curve_correspondence_map& cmap = cmap_it->second;

  CCMcit v1_in_cmap = cmap.find(v1);
  CGAL_precondition(v1_in_cmap != cmap.end());

  // this is the 'get_position' for a loop, so there can only be one position
  // per vertex.
  const Domain_pts_container& v1_bps = v1_in_cmap->second;
  CGAL_assertion(v1_bps.size() == 1);

  return *(v1_bps.begin());
}


template <typename C3T3, typename MD, typename Sf>
std::pair<typename Protect_edges_sizing_field<C3T3, MD, Sf>::Bare_point,
          typename Protect_edges_sizing_field<C3T3, MD, Sf>::Bare_point>
Protect_edges_sizing_field<C3T3, MD, Sf>::
get_positions_with_vertex_at_extremity(const Bare_point& known_point,
                                       const Domain_pts_container& extremity_points,
                                       const bool inverted_return_order,
                                       const Curve_index& curve_index,
                                       const CGAL::Orientation orientation) const
{
#if CGAL_MESH_3_PROTECTION_DEBUG & 2
  std::cerr << "get_positions_with_vertex_at_extremity()" << std::endl
            << "known_point: " << known_point << " on curve " << CGAL::oformat(curve_index)
            << " orientation: " << orientation
            << " inverted order ? " << std::boolalpha << inverted_return_order << std::endl;
#endif

  CGAL_precondition(!domain_.is_loop(curve_index));
  CGAL_precondition(extremity_points.size() == 2);

  // 'known_point' is a point on the curve, 'extremity_vertex' is the vertex that
  // corresponds to the corner that is both the beginning and (periodically)
  // the end of the curve with index 'curve_index'.

  const Bare_point& extr_pt = *(extremity_points.begin());
  const Bare_point& other_extr_pt = *(--extremity_points.end());
  CGAL_assertion(known_point != extr_pt);

#if CGAL_MESH_3_PROTECTION_DEBUG & 4
  std::cerr << "extr/other: " << extr_pt << " // " << other_extr_pt << std::endl;
#endif

  // Check for unknown orientation, which is indicated by a 'Zero'
  if(orientation == CGAL::ZERO)
  {
    // Three things:
    // - The only moment where we need to know positions without knowing orientation
    //   is when we are walking on edges of the mesh.
    // - We are in 1-cover and thus the edges of the mesh are small (no cycles, etc).
    // - We are in this function because one of the vertex 'v' corresponds to
    //   two points 'p1' and 'p2' (it's an extremity of a non-loop) on the curve
    //   of index 'curve_index'.
    // ==> the known point is necessarily much closer to one of the points 'p1' or 'p2'
    //     (because the edges are small).
    // ==> the point we seek is the extremity that is closest to the known point.

    const FT distance_to_extr_pt =
      domain_.curve_segment_length(known_point, extr_pt, curve_index,
                                   domain_.distance_sign(known_point, extr_pt, curve_index));
    const FT distance_to_other_extr_pt =
      domain_.curve_segment_length(known_point, other_extr_pt, curve_index,
                                   domain_.distance_sign(known_point, other_extr_pt, curve_index));

#if CGAL_MESH_3_PROTECTION_DEBUG & 4
    std::cerr << "distances: " << distance_to_extr_pt << " // "
              << distance_to_other_extr_pt << std::endl;
#endif

    CGAL_assertion(distance_to_extr_pt != distance_to_other_extr_pt);

    if(CGAL::abs(distance_to_extr_pt) > CGAL::abs(distance_to_other_extr_pt))
    {
      // 'other_extr_pt' is closer
      if(inverted_return_order)
        return std::make_pair(other_extr_pt, known_point);
      else
        return std::make_pair(known_point, other_extr_pt);
    }
    else
    {
      // 'extr_pt' is closer
      if(inverted_return_order)
        return std::make_pair(extr_pt, known_point);
      else
        return std::make_pair(known_point, extr_pt);
    }
  }

  // Case where the orientation is known
  CGAL::Orientation distance_sign;
  if(inverted_return_order)
    distance_sign = domain_.distance_sign(extr_pt, known_point, curve_index);
  else
    distance_sign = domain_.distance_sign(known_point, extr_pt, curve_index);

  if(distance_sign == orientation)
  {
    if(inverted_return_order)
      return std::make_pair(extr_pt, known_point);
    else
      return std::make_pair(known_point, extr_pt);
  }
  else // it's the other extremity
  {
    if(inverted_return_order)
      return std::make_pair(other_extr_pt, known_point);
    else
      return std::make_pair(known_point, other_extr_pt);
  }
}


template <typename C3T3, typename MD, typename Sf>
std::pair<typename Protect_edges_sizing_field<C3T3, MD, Sf>::Bare_point,
          typename Protect_edges_sizing_field<C3T3, MD, Sf>::Bare_point>
Protect_edges_sizing_field<C3T3, MD, Sf>::
get_positions(const Vertex_handle v1,
              const Vertex_handle v2,
              const Curve_index& curve_index,
              const CGAL::Orientation orientation) const
{
#if CGAL_MESH_3_PROTECTION_DEBUG & 2
  std::cerr << "get positions: " << &*v1 << " " << &*v2
            << " in curve " << curve_index << " and orientation: " << orientation << std::endl;
#endif

  typename Correspondence_map::const_iterator cmap_it =
    vertex_handle_to_full_space_position.find(curve_index);
  CGAL_precondition(cmap_it != vertex_handle_to_full_space_position.end());
  const Curve_correspondence_map& cmap = cmap_it->second;

  CCMcit v1_in_cmap = cmap.find(v1);
  CCMcit v2_in_cmap = cmap.find(v2);

  CGAL_precondition(v1_in_cmap != cmap.end());
  CGAL_precondition(v2_in_cmap != cmap.end());

  const Domain_pts_container& v1_bps = v1_in_cmap->second;
  const Domain_pts_container& v2_bps = v2_in_cmap->second;
  CGAL_precondition(v1_bps.size() == 1 || v1_bps.size() == 2);
  CGAL_precondition(v2_bps.size() == 1 || v2_bps.size() == 2);

  bool found_p1 = false, found_p2 = false;
  Bare_point p1, p2;

  if(v1_bps.size() == 1)
  {
    found_p1 = true;
    p1 = *(v1_bps.begin());
  }

  if(v2_bps.size() == 1)
  {
    found_p2 = true;
    p2 = *(v2_bps.begin());
  }

  if(found_p1 && found_p2)
    return std::make_pair(p1, p2);

  if(!found_p1 && !found_p2)
  {
    // A vertex can only correspond to two points if it is an extremity
    // that is the same periodically. Additionally, a curve only has two extremities.
    CGAL_assertion(v1 == v2);

    // Return two different points, with the correct orientation.
    p1 = *(v1_bps.begin());
    p2 = *(--v1_bps.end());
    CGAL_assertion(p1 != p2);

    CGAL::Orientation distance_sign = domain_.distance_sign(p1, p2, curve_index);
    if(distance_sign == orientation)
      return std::make_pair(p1, p2);
    else
      return std::make_pair(p2, p1);
  }

  // If we're here, only one point has been found, and the other corresponds
  // to an extremity of a curve that has the same beginning and end periodically
  // (seen by the fact that the vertex corresponds to two points).
  if(!found_p1)
  {
    return get_positions_with_vertex_at_extremity(p2, v1_bps, true /*inverted return order*/,
                                                  curve_index, orientation);
  }
  else
  {
    CGAL_assertion(!found_p2);
    return get_positions_with_vertex_at_extremity(p1, v2_bps, false /*inverted return order*/,
                                                  curve_index, orientation);
  }
}


template <typename C3T3, typename MD, typename Sf>
std::pair<typename Protect_edges_sizing_field<C3T3, MD, Sf>::Bare_point,
          typename Protect_edges_sizing_field<C3T3, MD, Sf>::Bare_point>
Protect_edges_sizing_field<C3T3, MD, Sf>::
get_positions_with_unknown_orientation(const Vertex_handle v1,
                                       const Vertex_handle v2,
                                       const Curve_index& curve_index) const
{
  // We don't know the orientation! This is problematic if a vertex corresponds
  // to two points (extremity of a periodic loop).
  CGAL::Orientation orientation = CGAL::ZERO;
  return get_positions(v1, v2, curve_index, orientation);
}


template <typename C3T3, typename MD, typename Sf>
boost::tuple<typename Protect_edges_sizing_field<C3T3, MD, Sf>::Bare_point,
             typename Protect_edges_sizing_field<C3T3, MD, Sf>::Bare_point,
             typename Protect_edges_sizing_field<C3T3, MD, Sf>::Bare_point>
Protect_edges_sizing_field<C3T3, MD, Sf>::
get_positions(const Vertex_handle v1,
              const Vertex_handle v2,
              const Vertex_handle v3,
              const Curve_index& curve_index,
              const CGAL::Orientation orientation) const
{
  Bare_point p1, p2_check, p2, p3;

  boost::tie(p1, p2) = get_positions(v1, v2, curve_index, orientation);
  boost::tie(p2_check, p3) = get_positions(v2, v3, curve_index, orientation);
  CGAL_assertion(p2_check == p2);

  return boost::make_tuple(p1, p2, p3);
}


////////////////////////////////////////////////////////////////////////////////
template <typename C3T3, typename MD, typename Sf>
void
Protect_edges_sizing_field<C3T3, MD, Sf>::
dump_dummy_points(const std::string filename) const
{
  std::ofstream dummy_out(filename.c_str());

  typename Gt::Construct_point_3 cp =
    c3t3_.triangulation().geom_traits().construct_point_3_object();
  typename Tr::All_vertices_iterator vit = c3t3_.triangulation().all_vertices_begin(),
                                     vend = c3t3_.triangulation().all_vertices_end();
  for(; vit!=vend; ++vit)
  {
    Index index = c3t3_.index(vit);
    const int* i = boost::get<int>(&index);
    if(i && *i == 0)
      dummy_out << cp(c3t3_.triangulation().point(vit)) << std::endl;
  }
}


template <typename C3T3, typename MD, typename Sf>
bool
Protect_edges_sizing_field<C3T3, MD, Sf>::
try_to_remove_dummy_vertex(const Vertex_handle dummy_vertex) const
{
  // 'dummy_vertex' must correspond to a dummy point
  CGAL_precondition_code(Index index = c3t3_.index(dummy_vertex);)
  CGAL_precondition_code(if(const int* i = boost::get<int>(&index)) {)
  CGAL_precondition(*i == 0);
  CGAL_precondition_code(})

  return c3t3_.triangulation().remove(dummy_vertex);
}


template <typename C3T3, typename MD, typename Sf>
typename Protect_edges_sizing_field<C3T3, MD, Sf>::FT
Protect_edges_sizing_field<C3T3, MD, Sf>::
get_maximum_weight(const Vertex_handle protection_vertex, const FT intended_weight) const
{
#if CGAL_MESH_3_PROTECTION_DEBUG & 4
  std::cerr << "Compute max weight(), intended weight: " << intended_weight << std::endl;
#endif

  FT max_possible_weight = intended_weight;

  // Note that incident_cells = finite_incident_cells for P3DT3
  std::vector<Cell_handle> finite_incident_cells;
  finite_incident_cells.reserve(64);
  c3t3_.triangulation().incident_cells(protection_vertex, std::back_inserter(finite_incident_cells));

  Tr_helpers tr_helpers;
  const FT nearest_sq_dist = tr_helpers.template get_sq_distance_to_closest_vertex
                               <CGAL_NTS internal::Has_member_visited<typename Dt::Vertex> >(
                                 c3t3_.triangulation(), protection_vertex, finite_incident_cells);

  if(nearest_sq_dist < intended_weight)
    max_possible_weight = nearest_sq_dist;

  if(max_possible_weight < minimal_weight_)
    max_possible_weight = minimal_weight_;

  CGAL_assertion_code(const Weighted_point& pvwp = c3t3_.triangulation().point(protection_vertex);)
  CGAL_assertion_code(const Bare_point& pvp =
    c3t3_.triangulation().geom_traits().construct_point_3_object()(pvwp);)
  CGAL_assertion_code(const int dim = get_dimension(protection_vertex);)
  CGAL_assertion_code(const Index index = c3t3_.index(protection_vertex);)
  CGAL_assertion_code(const FT w_max = CGAL::square(query_size(pvp, dim, index));)
  CGAL_assertion(max_possible_weight <= w_max);

  return max_possible_weight;
}


template <typename C3T3, typename MD, typename Sf>
bool
Protect_edges_sizing_field<C3T3, MD, Sf>::
change_ball_size_with_dummy_scaffholding(Vertex_handle& v, const FT w, bool is_special,
                                         const Weighted_point& dummy_point)
{
#if CGAL_MESH_3_PROTECTION_DEBUG & 2
  std::cerr << "Changing ball size with scaffholding dummy pt: " << dummy_point << std::endl;
#endif

  // If we are here, we have removed the dummy vertex and are still in 1-cover,
  // and we've computed the best weight for the vertex v. However, since "change_ball_size"
  // is basically a call to "tr.remove() + tr.insert()", it might itself change the cover !
  CGAL_precondition(c3t3_.triangulation().is_1_cover());

  // The dummy point should have been removed previously
  CGAL_precondition_code(Vertex_handle useless;)
  CGAL_precondition(!c3t3_.triangulation().is_vertex(dummy_point, useless));

  typename Gt::Compute_weight_3 cw =
    c3t3_.triangulation().geom_traits().compute_weight_3_object();

  // The trick is to re-insert the dummy point, change the ball size, and once again remove
  // the dummy vertex!
  const Weighted_point& v_wp = c3t3_.triangulation().point(v);
  const FT previous_weight = cw(v_wp);

  Vertex_handle dummy_vertex = insert_dummy_point(dummy_point);
  CGAL_assertion(dummy_vertex != Vertex_handle());

  if(change_ball_size(v, w, is_special))
  {
    // The ball changing size mide hide the dummy vertex, so check for existence first
    if(!c3t3_.triangulation().is_vertex(dummy_vertex) ||
       try_to_remove_dummy_vertex(dummy_vertex))
    {
#if CGAL_MESH_3_PROTECTION_DEBUG & 4
      std::cerr << "All good, changed weight and re-removed the dummy!" << std::endl;
#endif
      return true;
    }
    else
    {
#if CGAL_MESH_3_PROTECTION_DEBUG & 4
      std::cerr << "The dummy point now refuses to be removed... ?" << std::endl;
#endif
      if(!change_ball_size(v, previous_weight, is_special)) // restore the weight
      {
        std::cerr << "Error: failed to remove scaffholding dummy point" << std::endl;
        CGAL_assertion(false); // do not allow the weight restoration to fail
      }
    }
  }
  else
  {
#if CGAL_MESH_3_PROTECTION_DEBUG & 4
    std::cerr << "Couldn't change the weight... ?" << std::endl;
#endif
  }
  return false;
}


template <typename C3T3, typename MD, typename Sf>
typename Protect_edges_sizing_field<C3T3, MD, Sf>::Vertex_handle
Protect_edges_sizing_field<C3T3, MD, Sf>::
insert_dummy_point(const Weighted_point& dummy_point)
{
#if CGAL_MESH_3_PROTECTION_DEBUG & 4
  std::cerr << "(Re-)inserting dummy point: " << dummy_point << std::endl;
#endif
  typename Tr::Locate_type lt;
  int li, lj;
  Cell_handle ch = c3t3_.triangulation().locate(dummy_point, lt, li, lj);
  Vertex_handle new_v = c3t3_.triangulation().insert(dummy_point, lt, ch, li, lj);
  c3t3_.set_index(new_v, 0);
  return new_v;
}


template <typename C3T3, typename MD, typename Sf>
bool
Protect_edges_sizing_field<C3T3, MD, Sf>::
try_to_remove_close_dummy_vertex(Vertex_handle& protection_vertex,
                                 const Vertex_handle dummy_vertex,
                                 const FT intended_weight)
{
  CGAL_precondition(dummy_vertex != Vertex_handle());
  const Weighted_point dummy_point = c3t3_.triangulation().point(dummy_vertex); // intentional copy

  CGAL_precondition(c3t3_.triangulation().is_vertex(dummy_vertex));
  if(try_to_remove_dummy_vertex(dummy_vertex))
  {
#if CGAL_MESH_3_PROTECTION_DEBUG & 4
    std::cerr << "Successfully removed the dummy vertex" << std::endl;
#endif
    CGAL_expensive_assertion(c3t3_.triangulation().tds().is_valid());
    CGAL_expensive_assertion(c3t3_.triangulation().is_valid());
    CGAL_expensive_assertion(c3t3_.is_valid());
    CGAL_assertion(c3t3_.triangulation().is_1_cover());

    // 'intended_weight' might still be too much, reduce it as little as possible
    FT weight = get_maximum_weight(protection_vertex, intended_weight);
#if CGAL_MESH_3_PROTECTION_DEBUG & 4
    std::cerr << "maximum possible weight: " << weight << std::endl;
#endif

    // @todo (?) : if the weight is small compared to 'intended_weight' and the nearest
    // point is again a dummy, also try to remove that one ?

    if(change_ball_size_with_dummy_scaffholding(
         protection_vertex, weight, is_special(protection_vertex), dummy_point))
    {
#if CGAL_MESH_3_PROTECTION_DEBUG & 4
      std::cerr << "Successfuly removed the dummy point and changed the vertex weight" << std::endl;
#endif
      return true;
    }
    else
    {
      // 'change_ball_size_with_dummy_scaffholding()' re-inserts the dummy point on failure
#if CGAL_MESH_3_PROTECTION_DEBUG & 4
      std::cerr << "Failed to change the weight..." << std::endl;
#endif
    }
  }
  else
  {
#if CGAL_MESH_3_PROTECTION_DEBUG & 4
    std::cerr << "Failed to remove dummy point" << std::endl;
#endif
  }
  return false;
}


template <typename C3T3, typename MD, typename Sf>
bool
Protect_edges_sizing_field<C3T3, MD, Sf>::
try_to_move_dummy_vertex(const Vertex_handle dummy_vertex,
                         const Weighted_point& new_position)
{
  // Insert first to maximise likeliness of success
  Vertex_handle new_dummy = insert_dummy_point(new_position);

  if(!try_to_remove_dummy_vertex(dummy_vertex))
  {
    // If we didn't manage to remove the old position, remove the new position
    if(!try_to_remove_dummy_vertex(new_dummy))
    {
      std::cerr << "Error: failed to remove scaffholding dummy point" << std::endl;
      CGAL_assertion(false); // Do not allow failure in removing the new point
    }
    return false;
  }

#if CGAL_MESH_3_PROTECTION_DEBUG & 2
  std::cerr << "Moved dummy vertex to new position: " << new_position << std::endl;
#endif
  return true;
}


template <typename C3T3, typename MD, typename Sf>
bool
Protect_edges_sizing_field<C3T3, MD, Sf>::
move_dummy_vertex_away_from_corner_with_direction(Vertex_handle& protection_vertex,
                                                  const Vertex_handle dummy_vertex,
                                                  const FT intended_weight,
                                                  const Vector_3& dir,
                                                  const FT min_sq_dist)
{
  typename Gt::Construct_point_3 cp =
    c3t3_.triangulation().geom_traits().construct_point_3_object();

  int number_of_attempts = 10;
  FT max_sq_length = 2 * intended_weight;
  FT step_sq_length = max_sq_length / number_of_attempts;

  const Weighted_point& protection_position = c3t3_.triangulation().point(protection_vertex);
  const Weighted_point& dummy_point = c3t3_.triangulation().point(dummy_vertex);

  // First, try to find the farthest position that we can push the dummy point away
  // in the direction, with a limit of 2*intended_weight
  for(int i=0; i<number_of_attempts; ++i)
  {
    FT sq_distance = max_sq_length - i * step_sq_length;

    if(sq_distance < min_sq_dist)
      break; // since the next tries would be even closer, there's no point trying more

    const Bare_point moved_bp = cp(protection_position) + CGAL::sqrt(sq_distance) * dir;
#ifdef CGAL_PERIODIC_3_MESH_3_DEBUG_DUMMY_TREATMENT
    std::cerr << "Trying to move " << cp(dummy_point) << " to " << moved_bp << std::endl;
#endif

    // The modified position must not be in an existing ball
    Tr_helpers tr_helpers;
    if(tr_helpers.inside_protecting_balls(c3t3_.triangulation(), dummy_vertex, moved_bp))
      continue;

    // Insert the new dummy point
    const Weighted_point moved_dummy_point(moved_bp, 0.);
    Vertex_handle moved_dummy_vertex = insert_dummy_point(moved_dummy_point);
    if(moved_dummy_vertex == Vertex_handle()) // hidden point
      continue;

    // Now, try to remove the problematic dummy
    if(try_to_remove_dummy_vertex(dummy_vertex))
    {
#ifdef CGAL_PERIODIC_3_MESH_3_DEBUG_DUMMY_TREATMENT
      std::cerr << "With this move, the initial dummy vertex can be removed" << std::endl;
#endif
      FT weight = get_maximum_weight(protection_vertex, intended_weight);

      // @todo (?) : if the weight is small compared to 'intended_weight' and the nearest
      // point is again a dummy, also try to remove that one ?

      if(change_ball_size_with_dummy_scaffholding(
           protection_vertex, weight, is_special(protection_vertex), dummy_point))
      {
#ifdef CGAL_PERIODIC_3_MESH_3_DEBUG_DUMMY_TREATMENT
        std::cerr << "Successfully moved the dummy point and changed the vertex weight" << std::endl;
#endif
        return true;
      }
      else
      {
#ifdef CGAL_PERIODIC_3_MESH_3_DEBUG_DUMMY_TREATMENT
        std::cerr << "Failed to change the weight..." << std::endl;
#endif

        // 'change_ball_size_with_dummy_scaffholding()' re-inserts the initial dummy point on failure
        if(!try_to_remove_dummy_vertex(moved_dummy_vertex))
        {
          std::cerr << "Error: Failed to remove moved dummy vertex" << std::endl;
          CGAL_assertion(false); // do not allow the removal of the moved vertex to fail
        }
        continue;
      }
    }
    else // failed to remove the dummy point
    {
#ifdef CGAL_PERIODIC_3_MESH_3_DEBUG_DUMMY_TREATMENT
      std::cerr << "Can't remove the initial dummy point" << std::endl;
#endif
      if(try_to_remove_dummy_vertex(moved_dummy_vertex))
      {
        std::cerr << "Error: Failed to remove moved dummy vertex" << std::endl;
        CGAL_assertion(false); // do not allow the removal of the moved vertex to fail
      }
      continue;
    }
  }
  return false;
}


template <typename C3T3, typename MD, typename Sf>
bool
Protect_edges_sizing_field<C3T3, MD, Sf>::
try_to_move_close_dummy_vertex(Vertex_handle& protection_vertex,
                               const Vertex_handle dummy_vertex,
                               const FT intended_weight)
{
  typename Gt::Construct_point_3 cp =
    c3t3_.triangulation().geom_traits().construct_point_3_object();
  typename Gt::Compute_squared_distance_3 csd =
    c3t3_.triangulation().geom_traits().compute_squared_distance_3_object();

  const Weighted_point& pv_wp = c3t3_.triangulation().point(protection_vertex);
  const Weighted_point& dv_wp = c3t3_.triangulation().point(dummy_vertex);
  const Weighted_point canonical_dv_wp = c3t3_.triangulation().get_closest_point(pv_wp, dv_wp);
  CGAL_precondition(pv_wp != canonical_dv_wp);

  Vector_3 dir(cp(pv_wp), cp(canonical_dv_wp));
  dir /= CGAL::sqrt(dir * dir);
  FT sq_dist = csd(cp(pv_wp), cp(canonical_dv_wp));

  // First direction to move is the corner-dummy line.
  if(move_dummy_vertex_away_from_corner_with_direction(protection_vertex, dummy_vertex,
                                                       intended_weight, dir, sq_dist))
    return true;

  // @todo (?) If failure, try to move it away from the corner, along another (random) line

  return false;
}


template <typename C3T3, typename MD, typename Sf>
bool
Protect_edges_sizing_field<C3T3, MD, Sf>::
try_to_solve_close_dummy_point(Vertex_handle& protection_vertex,
                                      Vertex_handle dummy_vertex,
                                      FT intended_weight)
{
  // dummy_vertex must be a dummy point
  CGAL_precondition_code(Index index = c3t3_.index(dummy_vertex);)
  CGAL_precondition_code(if(const int* i = boost::get<int>(&index)) {)
  CGAL_precondition(*i == 0);
  CGAL_precondition_code(})

#if CGAL_MESH_3_PROTECTION_DEBUG & 2
      std::cerr << "Dummy vertex: " << c3t3_.triangulation().point(dummy_vertex) << " is too close "
            << "to vertex: " << c3t3_.triangulation().point(protection_vertex) << std::endl;
#endif

  if(try_to_remove_close_dummy_vertex(protection_vertex, dummy_vertex, intended_weight))
    return true;

  if(try_to_move_close_dummy_vertex(protection_vertex, dummy_vertex, intended_weight))
     return true;

  // If failed at both removing and moving... Well, that is one stinky dummy point...
  return false;
}


template <typename C3T3, typename MD, typename Sf>
void
Protect_edges_sizing_field<C3T3, MD, Sf>::
insert_corners()
{
  // Iterate on domain corners
  typedef std::vector< std::pair<Corner_index, Bare_point> > Initial_corners;

  Initial_corners corners;
  domain_.get_corners(std::back_inserter(corners));

  Dt dt(c3t3_.triangulation().domain());
  std::map<Bare_point, Dt_Vertex_handle> dt_point_vh_map;

  for(typename Initial_corners::iterator it = corners.begin(),
        end = corners.end() ; it != end ; ++it)
  {
    const Bare_point& p = it->second;
    Dt_Vertex_handle dtvh = dt.insert(c3t3_.triangulation().canonicalize_point(p));
    CGAL_assertion(dtvh != Dt_Vertex_handle());
    dt_point_vh_map[p] = dtvh; // note that p is not canonicalized here (there's no point to)
  }

  for(typename Initial_corners::iterator cit = corners.begin(),
        end = corners.end() ; cit != end ; ++cit)
  {
    const Bare_point& p = cit->second;
    const Corner_index& corner_index = cit->first;
    Index p_index = domain_.index_from_corner_index(corner_index);

#if CGAL_MESH_3_PROTECTION_DEBUG & 1
    std::cerr << "** treat corner #" << CGAL::oformat(p_index) << std::endl;
#endif

    // Get weight (the ball radius is given by the 'query_size' function)
    FT w = CGAL::square(query_size(p, 0, p_index));

#if CGAL_MESH_3_PROTECTION_DEBUG & 1
    std::cerr << "Point: " << p << std::endl;
    std::cerr << "Weight from sizing field: " << w << std::endl;
#endif

    //--------------------------------------------------------------------------
    // The following lines ensure that the weight w is small enough so that
    // corners balls do not intersect
    if(dt.number_of_vertices() >= 2)
    {
      Dt_Vertex_handle vh = dt_point_vh_map.at(p); // 'at' because it must be there

      // Note that incident_cells = finite_incident_cells for P3DT3
      std::vector<typename Dt::Cell_handle> finite_incident_cells;
      finite_incident_cells.reserve(64);
      dt.incident_cells(vh, std::back_inserter(finite_incident_cells));

      Dt_helpers helpers;
      const FT nearest_sq_dist = helpers.template get_sq_distance_to_closest_vertex
                                   <CGAL_NTS internal::Has_member_visited<typename Dt::Vertex> >(
                                     dt, vh, finite_incident_cells);

      w = (std::min)(w, nearest_sq_dist / FT(9));
#if CGAL_MESH_3_PROTECTION_DEBUG & 1
      std::cerr << "sq distance to nearest neighbor: " << nearest_sq_dist << std::endl;
#endif
    }

#if CGAL_MESH_3_PROTECTION_DEBUG & 1
    std::cerr << "Final weight: " << w << std::endl;
#endif

    // get curve indices here
    std::vector<Curve_index> incident_curves_indices;
    domain_.get_corner_incident_curves(corner_index,
                                       std::back_inserter(incident_curves_indices));
    if(incident_curves_indices.size() == 0)
    {
      // This corner won't be added to the correspondence map, but it's okay
      // because we will never need it.
#if CGAL_MESH_3_PROTECTION_DEBUG & 1
      std::cerr << "Warning: corner without any incident curve" << std::endl;
#endif
    }

    // Insert corner with ball (dim is zero because p is a corner)
    Vertex_handle v = smart_insert_point(p, w, 0, p_index, incident_curves_indices,
                                         CGAL::Emptyset_iterator()).first;
    CGAL_postcondition(v != Vertex_handle());

    // As C3t3::add_to_complex modifies the 'in_dimension' of the vertex,
    // we need to backup and re-set the 'is_special' marker after.
#if CGAL_MESH_3_PROTECTION_DEBUG & 2
    std::cerr << "Adding: " << c3t3_.triangulation().point(v) << " to complex" << std::endl;
#endif
    c3t3_.add_to_complex(v,cit->first);

    const bool special_ball = is_special(v);
    if(special_ball)
      set_special(v);
  }

#if CGAL_MESH_3_PROTECTION_DEBUG & 4
  std::cerr << "Correspondence map after insert corners: " << std::endl;
  dump_correspondence_map();
#endif
}


template <typename C3T3, typename MD, typename Sf>
template <typename Curve_index_container>
typename Protect_edges_sizing_field<C3T3, MD, Sf>::Vertex_handle
Protect_edges_sizing_field<C3T3, MD, Sf>::
insert_point(const Bare_point& p, const Weight& w, int dim, const Index& index,
             const Curve_index_container& curve_indices, const bool special_ball)
{
#if CGAL_MESH_3_PROTECTION_DEBUG & 1
  std::cerr << "insert_point()" << std::endl;
  std::cerr << "pos: " << p << " weight: " << w
            << " dim: " << dim << " index: " << CGAL::oformat(index) << std::endl;
#endif

  using CGAL::Mesh_3::internal::weight_modifier;

  // Convert the dimension if it was set to a negative value (marker for special balls).
  if(dim < 0)
    dim = -1 - dim;

  typedef typename Tr::size_type size_type;
  CGAL_USE_TYPE(size_type);

  // Insert point
#ifndef CGAL_PERIODIC_PROTECTION_ATTEMPT_TO_REMOVE_DUMMY_PTS
  CGAL_assertion_code(size_type nb_vertices_before = c3t3_.triangulation().number_of_vertices());
#endif

  typename Gt::Construct_weighted_point_3 cwp =
    c3t3_.triangulation().geom_traits().construct_weighted_point_3_object();

  const Weighted_point wp = cwp(p, w*weight_modifier);

  typename Tr::Locate_type lt;
  int li, lj;
  const typename Tr::Cell_handle ch = c3t3_.triangulation().locate(wp, lt, li, lj);

  CGAL_assertion(lt != Tr::VERTEX); // that case is handled before

  Vertex_handle v = c3t3_.triangulation().insert(wp, lt, ch, li, lj);
  CGAL_postcondition(Vertex_handle() != v);
#if CGAL_MESH_3_PROTECTION_DEBUG & 1
  std::cerr << "Inserted point: " << wp << " vertex_handle: " << &*v << std::endl;
#endif

#ifndef CGAL_PERIODIC_PROTECTION_ATTEMPT_TO_REMOVE_DUMMY_PTS
  // Contrary to Mesh_3, we are okay with hiding (dummy) points while changing weights
  CGAL_assertion (lt == Tr::VERTEX ||
                   c3t3_.triangulation().number_of_vertices() == (nb_vertices_before+1));
#endif

  // add the full space position
  insert_in_correspondence_map(v, p, curve_indices);

#if CGAL_MESH_3_PROTECTION_DEBUG & 1
  std::cerr << "Insertion of ";
  if(special_ball)
    std::cerr << "SPECIAL ";
  std::cerr << "protecting ball ";
  if(v == Vertex_handle())
    std::cerr << cwp(p,w*weight_modifier);
  else
    std::cerr << disp_vert(v);

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

  std::cerr << CGAL::oformat(index) << std::endl;
  if(v == Vertex_handle())
    std::cerr << "  HIDDEN!\n";
  std::cerr << "The weight was " << w << std::endl;
#endif // CGAL_MESH_3_PROTECTION_DEBUG

  c3t3_.set_dimension(v, dim);
  if(special_ball)
    set_special(v);
  c3t3_.set_index(v, index);

  unchecked_vertices_.insert(v);

  return v;
}


template <typename C3T3, typename MD, typename Sf>
template <typename ErasedVeOutIt, typename Curve_index_container>
std::pair<typename Protect_edges_sizing_field<C3T3, MD, Sf>::Vertex_handle,
          ErasedVeOutIt>
Protect_edges_sizing_field<C3T3, MD, Sf>::
smart_insert_point(const Bare_point& p, Weight w, int dim, const Index& index,
                   const Curve_index_container& curve_indices, ErasedVeOutIt out)
{
#if CGAL_MESH_3_PROTECTION_DEBUG & 1
  std::cerr << "smart_insert_point((" << p
            << "), w=" << w
            << ", dim=" << dim
            << ", index=" << CGAL::oformat(index) << ")\n";
#endif
  const Tr& tr = c3t3_.triangulation();

#if CGAL_MESH_3_PROTECTION_DEBUG & 2
  std::cerr << " Triangulation has " << tr.number_of_vertices() << " vertices" << std::endl;
#endif

  typename Gt::Compare_weighted_squared_radius_3 cwsr =
    tr.geom_traits().compare_weighted_squared_radius_3_object();
  typename Gt::Construct_point_3 cp =
    tr.geom_traits().construct_point_3_object();
  typename Gt::Construct_weighted_point_3 cwp =
    tr.geom_traits().construct_weighted_point_3_object();

  bool insert_a_special_ball = false; // will be passed to the function this->insert_point

  const Weighted_point wp0 = cwp(p); // with weight 0, used for locate()

  // ---------------------------------------------------------------------------
  // Check that new point will not be inside a power sphere
  typename Tr::Locate_type lt;
  int li, lj;
  Cell_handle ch = tr.locate(wp0, lt, li, lj);

#if CGAL_MESH_3_PROTECTION_DEBUG & 2
  std::cerr << "locate returns: " << &*ch << " lt: " << lt << " lilj: " << li << " " << lj << std::endl;
  for(int i=0; i<4; ++i)
    std::cerr << " -> " << c3t3_.triangulation().point(ch, i) << std::endl;
#endif

  if(lt == Tr::VERTEX)
  {
#if CGAL_MESH_3_PROTECTION_DEBUG & 1
    std::cerr << "Warning: point " << p << " already inserted" << std::endl;
#endif
    Vertex_handle v = ch->vertex(li);

    Index existing_vertex_index = c3t3_.index(v);
    const int* i = boost::get<int>(&existing_vertex_index);

    if(i && *i == 0)
    {
#if CGAL_MESH_3_PROTECTION_DEBUG & 1
      std::cerr << "Trying to insert on a dummy point" << std::endl;
#endif

      // The point we want to insert is unfortunately already in as a dummy...
      // Slightly perturb the dummy and call the same function again.

      // Note that incident_cells = finite_incident_cells for P3DT3
      std::vector<Cell_handle> finite_incident_cells;
      finite_incident_cells.reserve(64);
      c3t3_.triangulation().incident_cells(v, std::back_inserter(finite_incident_cells));

      Tr_helpers tr_helpers;
      const FT nearest_sq_dist = tr_helpers.template get_sq_distance_to_closest_vertex
                                   <CGAL_NTS internal::Has_member_visited<typename Dt::Vertex> >(
                                     c3t3_.triangulation(), v, finite_incident_cells);

      CGAL_assertion(nearest_sq_dist > 0.);

      // Use some scaffholding to safely to remove the dummy.
      Weighted_point close_pt(tr.point(v));
      Vector_3 rnd_direction(1., 0., 0.); // arbitrary direction

      typename Gt::Construct_translated_point_3 translate =
        c3t3_.triangulation().geom_traits().construct_translated_point_3_object();
      typename Gt::Construct_point_3 cp =
        c3t3_.triangulation().geom_traits().construct_point_3_object();

      // small perturbation
      close_pt = cwp(translate(cp(close_pt), 1e-5 * CGAL::sqrt(nearest_sq_dist) * rnd_direction));

      if(!try_to_move_dummy_vertex(v, close_pt))
      {
        // If this ever fails, probably need to pick a closer position than the
        // '1e-5' value above. Or some adaptive algorithm.
        CGAL_assertion(false);
      }

#if CGAL_MESH_3_PROTECTION_DEBUG & 1
      std::cerr << "Moved dummy point, calling smart_insert() again" << std::endl;
#endif
      return smart_insert_point(p, w, dim, index, curve_indices, out);
    }
    else
    {
      // The corner has already been inserted and necessary adjustements to its weight
      // have already been performed during its insertion and during the insertion
      // of other points. Only thing missing is to add it the correspondence map.
      insert_in_correspondence_map(v, p, curve_indices);
      return std::make_pair(v, out);
    }
  }

  Vertex_handle nearest_vh;
  FT sq_d;
  boost::tie(nearest_vh, sq_d) = tr.nearest_power_vertex_with_sq_distance(p, ch);
  CGAL_assertion(nearest_vh != Vertex_handle());
  CGAL_assertion(tr.point(nearest_vh) != cwp(tr.canonicalize_point(p)));

#if CGAL_MESH_3_PROTECTION_DEBUG & 2
  std::cerr << "Nearest power vertex of (" << p << ") is "
            << &*nearest_vh << " (" << c3t3_.triangulation().point(nearest_vh) << ") "
            << "at distance: " << sq_d << std::endl;

  Index nearest_vh_index = c3t3_.index(nearest_vh);
  int* i = boost::get<int>(&nearest_vh_index);
  if(i && *i == 0)
    std::cerr << "Nearest power vertex is a dummy point" << std::endl;
#endif

  // This will never happen for a dummy point
  while(cwsr(c3t3_.triangulation().point(nearest_vh), - sq_d) == CGAL::SMALLER &&
        ! is_special(nearest_vh))
  {
    CGAL_assertion(minimal_size_ > 0 || sq_d > 0);

    bool special_ball = false;
    if(minimal_weight_ != Weight() && sq_d < minimal_weight_)
    {
      sq_d = minimal_weight_;
      w = minimal_weight_;
      special_ball = true;
      insert_a_special_ball = true;
    }

    // Adapt size
    *out++ = nearest_vh;
    change_ball_size(nearest_vh, sq_d, special_ball);

    // Iterate
    ch = tr.locate(wp0, lt, li, lj, nearest_vh);
    boost::tie(nearest_vh, sq_d) = tr.nearest_power_vertex_with_sq_distance(p, ch);
    CGAL_assertion(nearest_vh != Vertex_handle());
  }

  if(is_special(nearest_vh) &&
     cwsr(c3t3_.triangulation().point(nearest_vh), - sq_d) == CGAL::SMALLER)
  {
    w = minimal_weight_;
    insert_a_special_ball = true;
  }

  // ---------------------------------------------------------------------------
  // Change w in order to be sure that no existing point will be included in (p,w)
  Vertex_handle nearest_vertex;
  FT min_sq_distance = w;

#ifdef CGAL_PERIODIC_PROTECTION_ATTEMPT_TO_REMOVE_DUMMY_PTS
  FT initial_weight = w;
  bool was_weight_reduced_due_to_close_dummy_vertex = false;
#endif

  // fill vertices_in_conflict_zone
  Vertex_set vertices_in_conflict_zone_set;
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
    for(int i=0; i<4; ++i)
    {
      Vertex_handle v = (*it)->vertex(i);
      if(!vertices_in_conflict_zone_set.insert(v).second)
        continue;

#ifdef CGAL_PERIODIC_PROTECTION_ATTEMPT_TO_REMOVE_DUMMY_PTS
      Index v_index = c3t3_.index(v);
      const int* id = boost::get<int>(&v_index);
      bool is_v_dummy_vertex(id && *id == 0);
#endif

      const FT sq_d = tr.min_squared_distance(p, cp(c3t3_.triangulation().point(v)));

      if(minimal_weight_ != Weight() && sq_d < minimal_weight_)
      {
        insert_a_special_ball = true;
        nearest_vertex = v;
        min_sq_distance = minimal_weight_;

#ifdef CGAL_PERIODIC_PROTECTION_ATTEMPT_TO_REMOVE_DUMMY_PTS
        was_weight_reduced_due_to_close_dummy_vertex = is_v_dummy_vertex;
#endif

        if(! is_special(v))
        {
          *out++ = v;
          change_ball_size(v, minimal_weight_, true /*special ball*/);
          ch = v->cell();
        }
      }
      else
      {
        if(sq_d < min_sq_distance)
        {
          nearest_vertex = v;
          min_sq_distance = sq_d;

#ifdef CGAL_PERIODIC_PROTECTION_ATTEMPT_TO_REMOVE_DUMMY_PTS
          was_weight_reduced_due_to_close_dummy_vertex = is_v_dummy_vertex;
#endif
        }
      }
    } // end of "for(int i=0; i<4; ++i)"
  }

  if(w > min_sq_distance)
  {
#if CGAL_MESH_3_PROTECTION_DEBUG & 1
    std::cerr << "smart_insert_point: weight " << w
              << " reduced to " << min_sq_distance
              << "\n (near existing point: "
              << c3t3_.triangulation().point(nearest_vertex) << ")\n";

    Index nearest_vertex_index = c3t3_.index(nearest_vertex);
    i = boost::get<int>(&nearest_vertex_index);
    if(i && *i == 0)
      std::cerr << "reduced due to dummy" << std::endl;
#endif
    w = min_sq_distance;
  }
  else
  {
#ifdef CGAL_PERIODIC_PROTECTION_ATTEMPT_TO_REMOVE_DUMMY_PTS
    was_weight_reduced_due_to_close_dummy_vertex = false;
#endif
  }

  if(lt != Tr::VERTEX)
  {
    CGAL_assertion_code(using CGAL::Mesh_3::internal::weight_modifier;)
    CGAL_assertion_code(std::vector<Vertex_handle> hidden_vertices;)
    CGAL_assertion_code(ch = tr.locate(wp0, lt, li, lj, ch);)
    CGAL_assertion_code(Weighted_point wpp = cwp(p, w*weight_modifier);)
    CGAL_assertion_code(tr.vertices_inside_conflict_zone(wpp,
                                                         ch,
                                                         std::back_inserter(hidden_vertices));)
    CGAL_assertion(hidden_vertices.empty());
  }

  const FT w_max = CGAL::square(query_size(p, dim, index));
  if(w > w_max)
  {
#if CGAL_MESH_3_PROTECTION_DEBUG & 1
    std::cerr << "smart_insert_point: weight " << w
              << " reduced to " << w_max << " (sizing field)\n";
#endif
    w = w_max;
#ifdef CGAL_PERIODIC_PROTECTION_ATTEMPT_TO_REMOVE_DUMMY_PTS
    was_weight_reduced_due_to_close_dummy_vertex = false;
#endif
  }

  if(w < minimal_weight_)
  {
    w = minimal_weight_;
    insert_a_special_ball = true;
#ifdef CGAL_PERIODIC_PROTECTION_ATTEMPT_TO_REMOVE_DUMMY_PTS
    was_weight_reduced_due_to_close_dummy_vertex = false;
#endif
  }

  // Insert the corner
  Vertex_handle new_vertex = insert_point(p, w, dim, index, curve_indices, insert_a_special_ball);

#ifdef CGAL_PERIODIC_PROTECTION_ATTEMPT_TO_REMOVE_DUMMY_PTS
  if(was_weight_reduced_due_to_close_dummy_vertex)
  {
    // The weight was reduced due to a close dummy point, try to remove that dummy
    // point and restore the weight to its former glory.
    if(try_to_solve_close_dummy_point(new_vertex, nearest_vertex, initial_weight))
    {
#ifdef CGAL_PERIODIC_3_MESH_3_DEBUG_DUMMY_TREATMENT
      std::cerr << "The dummy vertex was (re)moved and the weight was increased!" << std::endl;
#endif
    }
    else
    {
#ifdef CGAL_PERIODIC_3_MESH_3_DEBUG_DUMMY_TREATMENT
      std::cerr << "The dummy vertex is too tough to handle..." << std::endl;
#endif
    }

    CGAL_expensive_postcondition(c3t3_.triangulation().tds().is_valid());
    CGAL_expensive_postcondition(c3t3_.triangulation().is_valid());
    CGAL_expensive_postcondition(c3t3_.is_valid());
  }
#endif

  return std::make_pair(new_vertex, out);
}


template <typename C3T3, typename MD, typename Sf>
void
Protect_edges_sizing_field<C3T3, MD, Sf>::
insert_balls_on_edges()
{
  // Get features
  typedef std::tuple<Curve_index,
                             std::pair<Bare_point,Index>,
                             std::pair<Bare_point,Index> >    Feature_tuple;
  typedef std::vector<Feature_tuple>                          Input_features;

  Input_features input_features;
  domain_.get_curves(std::back_inserter(input_features));

  // Interate on edges
  for(typename Input_features::iterator fit = input_features.begin(),
       end = input_features.end() ; fit != end ; ++fit)
  {
    const Curve_index& curve_index = std::get<0>(*fit);
    if(! is_treated(curve_index))
    {
#if CGAL_MESH_3_PROTECTION_DEBUG & 1
      std::cerr << "** treat curve #" << curve_index << std::endl;
      std::cerr << "is it a loop? " << domain_.is_loop(curve_index) << std::endl;
#endif
      const Bare_point& p = std::get<1>(*fit).first;
      const Bare_point& q = std::get<2>(*fit).first;

      const Index& p_index = std::get<1>(*fit).second;
      const Index& q_index = std::get<2>(*fit).second;

      Vertex_handle vp,vq;
      if(! domain_.is_loop(curve_index))
      {
        vp = get_vertex_corner_from_point(p, p_index);
        vq = get_vertex_corner_from_point(q, q_index);
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

          FT curve_length = domain_.curve_length(curve_index);

          Bare_point other_point = domain_.construct_point_on_curve(p,
                                                                    curve_index,
                                                                    curve_length / 2);
          const FT po_distance = compute_distance(p, other_point);
          p_size = (std::min)(p_size, po_distance / 3);
          vp = smart_insert_point(p,
                                  CGAL::square(p_size),
                                  1 /*dim*/,
                                  p_index,
                                  curve_index,
                                  CGAL::Emptyset_iterator()).first;
        }
        // No 'else' because in that case 'is_vertex(..)' already filled
        // the variable 'vp'.
        vq = vp;
      }

      // Insert balls and set treated
//      if(do_balls_intersect(vp, vq))
//      {
//        CGAL_assertion(is_special(vp) || is_special(vq));
//      }
//      else
      {
        insert_balls(vp, vq, curve_index, CGAL::POSITIVE, Emptyset_iterator());
      }
      set_treated(curve_index);
    }

    // std::stringstream s;
    // s << "dump-mesh-curve-" << curve_index << ".binary.cgal";
    // debug_dump_c3t3(s.str(), c3t3_);
  }
}


template <typename C3T3, typename MD, typename Sf>
typename Protect_edges_sizing_field<C3T3, MD, Sf>::Vertex_handle
Protect_edges_sizing_field<C3T3, MD, Sf>::
get_vertex_corner_from_point(const Bare_point& p, const Index&) const
{
  typename Gt::Construct_weighted_point_3 cwp =
    c3t3_.triangulation().geom_traits().construct_weighted_point_3_object();

  // Get vertex_handle associated to corner (dim=0) point
  Vertex_handle v;
  CGAL_assertion_code(bool q_found =)

  // Let the weight be 0, because is_vertex only locates the point, and
  // check that the location type is VERTEX.
  c3t3_.triangulation().is_vertex(cwp(p), v);

  CGAL_assertion(q_found);
  return v;
}


template <typename C3T3, typename MD, typename Sf>
template <typename ErasedVeOutIt>
ErasedVeOutIt
Protect_edges_sizing_field<C3T3, MD, Sf>::
insert_balls(const Vertex_handle& vp,
             const Vertex_handle& vq,
             const Curve_index& curve_index,
             const CGAL::Orientation orientation,
             ErasedVeOutIt out)
{
  // Get size of p & q
  Bare_point vpp, vqp;
  boost::tie(vpp, vqp) = get_positions(vp, vq, curve_index, orientation);

  const FT sp = get_radius(vp);
  const FT sq = get_radius(vq);

  CGAL_assertion(vpp != vqp || domain_.is_loop(curve_index));

  // Compute geodesic distance
  const FT pq_length = (vp == vq) ?
                         domain_.curve_length(curve_index) :
                         domain_.curve_segment_length(vpp, vqp, curve_index, orientation);

  // Insert balls
  return
    (sp <= sq) ?
      insert_balls(vp, vq, sp, sq, pq_length, orientation, curve_index, out) :
      insert_balls(vq, vp, sq, sp, pq_length, -orientation, curve_index, out);
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
             const CGAL::Orientation d_sign,
             const Curve_index& curve_index,
             ErasedVeOutIt out)
{
#if CGAL_MESH_3_PROTECTION_DEBUG & 1
  std::cerr << "insert_balls(vp=" << disp_vert(vp) << ",\n"
            << "             vq=" << disp_vert(vq) << ",\n"
            << "             sp=" << sp << ",\n"
            << "             sq=" << sq << ",\n"
            << "             d=" << d << ",\n"
            << "             d_sign=" << d_sign << ",\n"
            << "             curve_index=" << curve_index << ")\n";
  std::cerr << "nonlinear growth? " << std::boolalpha << nonlinear_growth_of_balls << '\n'
            << "refine_balls_iteration_nb: " << refine_balls_iteration_nb << '\n'
            << "max_nb_vertices_to_reevaluate_size" << Mesh_3::internal::max_nb_vertices_to_reevaluate_size << std::endl;
#endif
  CGAL_precondition(d > 0);
  CGAL_precondition(sp <= sq);

  Bare_point vpp, vqp;
  boost::tie(vpp, vqp) = get_positions(vp, vq, curve_index, d_sign);

#if ! defined(CGAL_NO_PRECONDITIONS)
  if(sp <= 0)
  {
    std::stringstream msg;;
    msg.precision(17);
    msg << "Error: the mesh sizing field is null at point (";
    msg << vpp << ")!";
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
  // if(minimal_weight_ != 0 && n == 0) return;

  if(nonlinear_growth_of_balls && refine_balls_iteration_nb < 3)
  {
    // This block tries not to apply the general rule that the size of
    // protecting balls is a linear interpolation of the size of protecting
    // balls at corner. When the curve segment is long enough, pick a point
    // at the middle and choose a new size.
    if(n >= Mesh_3::internal::max_nb_vertices_to_reevaluate_size &&
       d >= (Mesh_3::internal::max_nb_vertices_to_reevaluate_size * minimal_weight_))
    {
#if CGAL_MESH_3_PROTECTION_DEBUG & 1
      std::cerr << "Number of to-be-inserted balls is: "
                << n << "\n  between points ("
                << vpp << ") and (" << vqp
                << ") (arc length: "
                << domain_.curve_segment_length(vpp,
                                                vqp,
                                                curve_index, d_sign)
                << ")\n";
#endif
      const Bare_point new_point =
        domain_.construct_point_on_curve(vpp, curve_index, d_sign * d / 2);
      const int dim = 1; // new_point is on edge
      const Index index = domain_.index_from_curve_index(curve_index);
      const FT point_weight = CGAL::square(size_(new_point, dim, index));

#if CGAL_MESH_3_PROTECTION_DEBUG & 1
      std::cerr << "  middle point: " << new_point << std::endl;
      std::cerr << "  new weight: " << point_weight << std::endl;
#endif
      std::pair<Vertex_handle, ErasedVeOutIt> pair =
        smart_insert_point(new_point, point_weight, dim, index, curve_index, out);
      Vertex_handle new_vertex = pair.first;

      out = pair.second;
      const FT sn = get_radius(new_vertex);
      if(sp <= sn)
        out = insert_balls(vp, new_vertex, sp, sn, d/2, d_sign, curve_index, out);
      else
        out = insert_balls(new_vertex, vp, sn, sp, d/2, -d_sign, curve_index, out);

      if(sn <= sq)
        out = insert_balls(new_vertex, vq, sn, sq, d/2, d_sign, curve_index, out);
      else
        out = insert_balls(vq, new_vertex, sq, sn, d/2, -d_sign, curve_index, out);

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

  // if((0 == n) &&
  //      ((d >= sp+sq) || !is_sampling_dense_enough(vp, vq)))

  // If there is some place to insert one point, insert it
  if((0 == n) && (d >= sp+sq))
  {
    n = 1;
    step_size = sp + (d-sp-sq) / FT(2);
    pt_dist = d_signF * step_size;
    norm_step_size = step_size;
  }
  else if(vp == vq && n == 1)
  {
    // In case we sample a full cycle, we want to add at least 2
    // balls, equally distributed.
    n = 2;
    step_size = d / FT(n+1);
    pt_dist = d_signF * step_size;
    norm_step_size = step_size;
  }
  else
  {
    CGAL_assertion_code(using boost::math::float_prior);
    CGAL_assertion(n==0 ||
                   dleft_frac >= float_prior(float_prior(1.)));
  }

  // Launch balls
  for(int i = 1 ; i <= n ; ++i)
  {
    // New point position
    Bare_point new_point =
      domain_.construct_point_on_curve(vpp, curve_index, pt_dist);

    // Weight (use as size the min between norm_step_size and linear interpolation)
    FT current_size = (std::min)(norm_step_size, sp + CGAL::abs(pt_dist)/d*(sq-sp));
    FT point_weight = current_size * current_size;

    // Index and dimension
    Index index = domain_.index_from_curve_index(curve_index);
    int dim = 1; // new_point is on edge

    // Insert point into c3t3
    std::pair<Vertex_handle, ErasedVeOutIt> pair =
      smart_insert_point(new_point, point_weight, dim, index, curve_index, out);
    Vertex_handle new_vertex = pair.first;
    out = pair.second;

    // Add edge to c3t3
    if(!c3t3_.is_in_complex(prev, new_vertex))
    {
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
  if(vp != vq || n > 1)
  {
    if(!c3t3_.is_in_complex(prev, vq))
    {
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
  Tr& tr = c3t3_.triangulation();

#if CGAL_MESH_3_PROTECTION_DEBUG & 1
  std::cerr << "Refine balls() " << std::endl;
#endif

  // Loop
  bool restart = true;
  using CGAL::Mesh_3::internal::refine_balls_max_nb_of_loops;
  this->refine_balls_iteration_nb = 0;
  while((!unchecked_vertices_.empty() || restart) &&
          this->refine_balls_iteration_nb < refine_balls_max_nb_of_loops)
  {
#ifdef CGAL_MESH_3_DUMP_FEATURES_PROTECTION_ITERATIONS
    std::ostringstream oss;
    oss << "dump_protecting_balls_" << refine_balls_iteration_nb << ".cgal";
    std::ofstream outfile(oss.str().c_str(), std::ios_base::binary | std::ios_base::out);
    CGAL::Mesh_3::save_binary_file(outfile, c3t3_, true);
    outfile.close();
#endif //CGAL_MESH_3_DUMP_FEATURES_PROTECTION_ITERATIONS

#if CGAL_MESH_3_PROTECTION_DEBUG & 1
    std::cerr << "RESTART REFINE LOOP (" << refine_balls_iteration_nb << ")\n"
              << "\t unchecked_vertices size: " << unchecked_vertices_.size() <<"\n";
#endif
    ++refine_balls_iteration_nb;
    restart = false;
    boost::unordered_map<Vertex_handle, FT, Hash_fct> new_sizes;

    for(typename Tr::Finite_edges_iterator eit = tr.finite_edges_begin(),
        end = tr.finite_edges_end(); eit != end; ++eit)
    {
      Vertex_handle va = eit->first->vertex(eit->second);
      Vertex_handle vb = eit->first->vertex(eit->third);

#if CGAL_MESH_3_PROTECTION_DEBUG & 16
      std::cerr << "Treat edge: " << c3t3_.triangulation().point(va) << " || "
                                  << c3t3_.triangulation().point(vb) << std::endl;
      std::cerr << "Intersecting edges ? " << do_balls_intersect(va, vb) << std::endl;
#endif

      // If those vertices are not adjacent in the complex but still intersect,
      // must reduce the balls size
      if(non_adjacent_but_intersect(va, vb))
      {
        using CGAL::Mesh_3::internal::distance_divisor;

        // Compute correct size of balls
        const FT ab = compute_distance(va, vb);

        FT ra = get_radius(va);
        FT rb = get_radius(vb);
        FT sa_new = (std::min)(ab/distance_divisor, ra);
        FT sb_new = (std::min)(ab/distance_divisor, rb);

#if CGAL_MESH_3_PROTECTION_DEBUG & 16
        std::cerr << "ab: " << ab << std::endl;
        std::cerr << "distance divisor: " << distance_divisor << std::endl;
        std::cerr << "ra | rb: " << ra << " " << rb << std::endl;
        std::cerr << "sa_new | sb_new: " << sa_new << " " << sb_new << std::endl;
#endif

        // In case of va or vb have already been in conflict, keep minimal size
        if(new_sizes.find(va) != new_sizes.end())
        {
          sa_new = (std::min)(sa_new, new_sizes[va]);
        }

        if(new_sizes.find(vb) != new_sizes.end())
        {
          sb_new = (std::min)(sb_new, new_sizes[vb]);
        }

#if CGAL_MESH_3_PROTECTION_DEBUG & 16
        std::cerr << "refine_balls: " << disp_vert(va) << " and "
                  << disp_vert(vb) << " are non-adjacent but do intersect\n";
#endif // CGAL_MESH_3_PROTECTION_DEBUG & 16

        // Store new_sizes for va and vb
        if(sa_new != ra)
        {
#if CGAL_MESH_3_PROTECTION_DEBUG & 16
          std::cerr << "  new_sizes[" << disp_vert(va) << ":"
                    << new_sizes[va] << "\n";
#endif // CGAL_MESH_3_PROTECTION_DEBUG & 16
          new_sizes[va] = sa_new;
        }

        if(sb_new != rb)
        {
#if CGAL_MESH_3_PROTECTION_DEBUG & 16
          std::cerr << "  new_sizes[" << disp_vert(vb) << ":"
                    << new_sizes[vb] << "\n";
#endif // CGAL_MESH_3_PROTECTION_DEBUG & 16
          new_sizes[vb] = sb_new;
        }
      }
    }

    // The std::map with Vertex_handle as the key is not robust, because
    // the time stamp of vertices can change during the following loop. The
    // solution is to copy it in a vector.
    std::vector<std::pair<Vertex_handle, FT> >
      new_sizes_copy(new_sizes.begin(), new_sizes.end());
    new_sizes.clear();

    // Update size of balls
    for(typename std::vector<std::pair<Vertex_handle,FT> >::iterator
          it = new_sizes_copy.begin(), end = new_sizes_copy.end(); it != end ; ++it)
    {
      Vertex_handle v = it->first;
      const FT new_size = it->second;
      // Set size of the ball to new value
      if(minimal_size_ != FT() && new_size < minimal_size_)
      {
        if(!is_special(v))
        {
          change_ball_size(v, minimal_weight_, true); // special ball

          // Loop will have to be run again
          restart = true;
        }
      }
      else
      {
        change_ball_size(v, CGAL::square(new_size));

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
  {
#if CGAL_MESH_3_PROTECTION_DEBUG & 1
    std::cerr << "Warning : features protection has reached maximal "
              << " number of loops." << std::endl
              << "          It might result in a crash." << std::endl;
#endif
  }
} // end refine_balls()


template <typename C3T3, typename MD, typename Sf>
bool
Protect_edges_sizing_field<C3T3, MD, Sf>::
non_adjacent_but_intersect(const Vertex_handle& va, const Vertex_handle& vb) const
{
  CGAL_precondition(va != vb);
  CGAL_precondition(c3t3_.triangulation().point(va) != c3t3_.triangulation().point(vb));

  if(! c3t3_.is_in_complex(va, vb))
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
  typename Gt::Compute_weight_3 cw =
    c3t3_.triangulation().geom_traits().compute_weight_3_object();
  typename Gt::Construct_point_3 cp =
    c3t3_.triangulation().geom_traits().construct_point_3_object();

  // Using the vertex->full space position map is not enough: a sphere at 0.9
  // with radius 0.1 does intersect a sphere at 1.1 (i.e. 0.1) with radius 0.1.
  // --> need min_squared_distance
  const Weighted_point& wpa = c3t3_.triangulation().point(va);
  const Weighted_point& wpb = c3t3_.triangulation().point(vb);
  CGAL_precondition(va != vb && wpa != wpb);

  return (c3t3_.triangulation().min_squared_distance(cp(wpa), cp(wpb)) <=
            CGAL::square(CGAL::sqrt(cw(wpa)) + CGAL::sqrt(cw(wpa))));
}

template <typename C3T3, typename MD, typename Sf>
bool
Protect_edges_sizing_field<C3T3, MD, Sf>::
change_ball_size(Vertex_handle& v, const FT squared_size, const bool special_ball)
{
  Weight w(squared_size);

  // Check if there is something to do
  // if(c3t3_.triangulation().point(v).weight() == w)
  //   return true;

#if CGAL_MESH_3_PROTECTION_DEBUG & 1
  std::cerr << "change_ball_size(v=" << disp_vert(v)
            << " dim=" << c3t3_.in_dimension(v)
            << " index=" << CGAL::oformat(c3t3_.index(v))
            << " ,\n"
            << "                 (squared) size=" << w
            << ", special_ball=" << std::boolalpha << special_ball << std::endl;
#endif

  // Get incident vertices along c3t3 edge
  Adjacent_vertices adjacent_vertices;
  c3t3_.adjacent_vertices_in_complex(v, std::back_inserter(adjacent_vertices));

  // Remove incident edges from complex
  for(typename Adjacent_vertices::iterator vit = adjacent_vertices.begin(),
       vend = adjacent_vertices.end() ; vit != vend ; ++vit)
  {
    c3t3_.remove_from_complex(v, vit->first);
  }

  typename Gt::Construct_point_3 cp =
    c3t3_.triangulation().geom_traits().construct_point_3_object();

  // Store point data
  int dim = get_dimension(v);
  const Bare_point p = cp(c3t3_.triangulation().point(v)); // intentional copy

  // Remove v from the set of corners
  boost::optional<Corner_index> corner_index = boost::make_optional(false, Corner_index());
  if(c3t3_.is_in_complex(v))
  {
    corner_index = c3t3_.corner_index(v);
    c3t3_.remove_from_complex(v);
  }

  unchecked_vertices_.erase(v);

  // Attempt to change weight, abort if the remove returns 'false', which means
  // that the triangulation would switch to 27-sheets if the point were removed.
  //
  // Note that there is no need to change the correspondence map here because
  // we are directly re-inserting at the same position and with the same vertex handle.
  if(!c3t3_.triangulation().remove(v))
  {
    // This should only happen if we're fiddling with dummy points
#if CGAL_MESH_3_PROTECTION_DEBUG & 4
    std::cerr << "Warning: change_ball_size() can't remove vertex "
              << "as it would switch the triangulation's cover" << std::endl;
#endif

    // Restore the different properties
    unchecked_vertices_.insert(v);

    // Restore corner, if needed
    if(corner_index)
      c3t3_.add_to_complex(v, *corner_index);

    // Restore connectivity in c3t3
    for(typename Adjacent_vertices::iterator it = adjacent_vertices.begin(),
         end = adjacent_vertices.end() ; it != end ; ++it)
    {
      c3t3_.add_to_complex(v, it->first, it->second);
    }
    return false;
  }

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

  // We do not need to change anything in the correspondence map, so leave the vector empty
  Index index = c3t3_.index(v);
  std::vector<Curve_index> incident_curves_indices;

  Vertex_handle new_v = insert_point(p, w , dim, index, incident_curves_indices, special_ball);
  CGAL_assertion(hidden_vertices.empty());

  CGAL_assertion((! special_ball) || is_special(new_v));

  // TODO: ensure that this condition is always satisfied (Pedro's code ?)
  CGAL_assertion(v == new_v);
  //  new_v->set_meshing_info(w);

  // Restore v in corners
  if(corner_index)
  {
    // As C3t3::add_to_complex modifies the 'in_dimension' of the vertex,
    // we need to backup and re-set the 'is_special' marker after.
    c3t3_.add_to_complex(new_v, *corner_index);

    const bool special_ball = is_special(new_v);
    if(special_ball)
        set_special(new_v);
  }

  // Restore c3t3 edges
  for(typename Adjacent_vertices::iterator it = adjacent_vertices.begin(),
       end = adjacent_vertices.end() ; it != end ; ++it)
  {
    // Restore connectivity in c3t3
    c3t3_.add_to_complex(new_v, it->first, it->second);
  }

  // Update unchecked vertices
  unchecked_vertices_.insert(new_v);
  return true;
}

namespace details {

// Functor used by Protect_edges_sizing_field::check_and_repopulate_edges, below
template <typename Set>
class Erase_element_from_set
{
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
  std::cerr << c3t3_.number_of_vertices_in_complex() << " corners" << std::endl;
  std::cerr << c3t3_.triangulation().number_of_vertices() << " vertices" << std::endl;
#endif

  Vertex_set vertices;
  std::copy(unchecked_vertices_.begin(), unchecked_vertices_.end(),
            std::inserter(vertices,vertices.begin()));
  unchecked_vertices_.clear();

  // Fix edges
  while(!vertices.empty())
  {
    Vertex_handle v = *vertices.begin();
    vertices.erase(vertices.begin());

    details::Erase_element_from_set<Vertex_set> erase_from_vertices(vertices);

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
  if(c3t3_.is_in_complex(v))
  {
    return repopulate_edges_around_corner(v, out);
  }

  // Get incident vertices along c3t3 edge
  Adjacent_vertices adjacent_vertices;
  c3t3_.adjacent_vertices_in_complex(v, std::back_inserter(adjacent_vertices));

  // The size of 'adjacent_vertices' can be 0 if v is a ball that covers
  // entirely a closed curve.
  // The size can also be 1 if the curve is a cycle, and the temporary
  // mesh is only two balls on the cycle: then each ball has only one
  // neighbor.
  CGAL_assertion(v->is_special() || adjacent_vertices.size() < 3);

  if(adjacent_vertices.size() == 0)
    return out;

  typename Adjacent_vertices::const_iterator ivb = adjacent_vertices.begin();
  typename Adjacent_vertices::const_iterator ivl = --(adjacent_vertices.end());

#if CGAL_MESH_3_PROTECTION_DEBUG & 4
  std::cerr << adjacent_vertices.size() << " incident vertices" << std::endl;
  std::cerr << "ivb: " << *(ivb->first) << " // " << ivb->second << std::endl;
  std::cerr << "ivl: " << *(ivl->first) << " // " << ivl->second << std::endl;
#endif

  const Curve_index& curve_index = ivb->second;
  CGAL_assertion(ivl->second == curve_index);

  // Walk along edge to find the edge piece which is not correctly sampled
  typedef std::list<Vertex_handle> Vertex_list;
  Vertex_list to_repopulate;
  to_repopulate.push_front(v);

  const Vertex_handle& previous = ivb->first;
  const Vertex_handle& next = ivl->first;

  const CGAL::Orientation orientation = orientation_of_walk(v, next, curve_index);

  // Walk following direction (v,previous)
  walk_along_edge(v, previous, curve_index, -orientation, std::front_inserter(to_repopulate));

#if CGAL_MESH_3_PROTECTION_DEBUG & 1
  std::cerr <<  "to_repopulate.size()=" << to_repopulate.size() << "\n";
#endif // CGAL_MESH_3_PROTECTION_DEBUG

  // Check whether a complete circle has been discovered or not
  if(to_repopulate.size() == 1
      || to_repopulate.front() != to_repopulate.back())
  {
    // Walk in other direction (v,next)
    walk_along_edge(v, next, curve_index, orientation,
                    std::back_inserter(to_repopulate));

#if CGAL_MESH_3_PROTECTION_DEBUG & 1
    std::cerr <<  "to_repopulate.size()=" << to_repopulate.size() << "\n";
#endif // CGAL_MESH_3_PROTECTION_DEBUG
  }

  // If only v is in to_repopulate, there is nothing to do
  if(to_repopulate.size() == 1)
  {
    *out++ = *to_repopulate.begin();
    return out;
  }

  // Store erased vertices
  // out = std::copy(to_repopulate.begin(), to_repopulate.end(), out);

  // Repopulate edge
  out = analyze_and_repopulate(to_repopulate.begin(),
                               --to_repopulate.end(),
                               curve_index,
                               orientation,
                               out);

  return out;
}


template <typename C3T3, typename MD, typename Sf>
bool
Protect_edges_sizing_field<C3T3, MD, Sf>::
is_sampling_dense_enough(const Vertex_handle& v1, const Vertex_handle& v2,
                         const Curve_index& curve_index,
                         const CGAL::Orientation orientation) const
{
  using CGAL::Mesh_3::internal::min_intersection_factor;
  CGAL_precondition(c3t3_.curve_index(v1,v2) == curve_index);

  typename Gt::Compute_weight_3 cw =
    c3t3_.triangulation().geom_traits().compute_weight_3_object();

  // Get sizes
  FT size_v1 = get_radius(v1);
  FT size_v2 = get_radius(v2);

  CGAL_assertion(get_dimension(v1) != 1 ||
                 curve_index == domain_.curve_index(v1->index()));
  CGAL_assertion(get_dimension(v2) != 1 ||
                 curve_index == domain_.curve_index(v2->index()));

  Bare_point v1p, v2p;
  boost::tie(v1p, v2p) = get_positions(v1, v2, curve_index, orientation);

  FT arc_length = domain_.curve_segment_length(v1p,
                                               v2p,
                                               curve_index,
                                               orientation);

  // Sufficient condition so that the curve portion between v1 and v2 is
  // inside the union of the two balls.
  if(arc_length > (size_v1 + size_v2))
  {
#if CGAL_MESH_3_PROTECTION_DEBUG & 1
    std::cerr << "Note: on curve #" << curve_index << ", between "
              << disp_vert(v1) << " and " << disp_vert(v2) << ", the "
              << "arc length is " << arc_length << " and the"
              << " sum of radii is " << size_v1 + size_v2 << std::endl;
#endif // CGAL_MESH_3_PROTECTION_DEBUG & 1

    if(!do_balls_intersect(v2, v1))
    {
#if CGAL_MESH_3_PROTECTION_DEBUG & 1
      std::cerr << "     And spheres do NOT intersect!\n";
#endif // CGAL_MESH_3_PROTECTION_DEBUG & 1
      return false;
    }

    const Weighted_point& v1_wp = c3t3_.triangulation().point(v1);
    const Weighted_point& v2_wp = c3t3_.triangulation().point(v2);

    const bool cov = domain_.is_curve_segment_covered(curve_index,
                                                      orientation,
                                                      v1p, v2p,
                                                      cw(v1_wp), cw(v2_wp));
#if CGAL_MESH_3_PROTECTION_DEBUG & 1
    if(cov)
      std::cerr << "      But the curve is locally covered\n";
    else
      std::cerr << "      And the curve is NOT locally covered\n";
#endif
    return cov;
  }

  const FT distance_v1v2 = compute_distance(v1, v2, curve_index, orientation);

  // Ensure size_v1 > size_v2
  if(size_v1 < size_v2)
    std::swap(size_v1, size_v2);

  // Check if balls intersect
  return distance_v1v2 < (FT(min_intersection_factor) * size_v2 + size_v1);
}


template <typename C3T3, typename MD, typename Sf>
CGAL::Orientation
Protect_edges_sizing_field<C3T3, MD, Sf>::
orientation_of_walk(const Vertex_handle& start,
                    const Vertex_handle& next,
                    Curve_index curve_index) const
{
#if CGAL_MESH_3_PROTECTION_DEBUG & 4
  std::cerr << "orientation_of_walk("
            << c3t3_.triangulation().point(start) << ", "
            << c3t3_.triangulation().point(next) << ")" << std::endl;
#endif

  Bare_point start_p, next_p;
  boost::tie(start_p, next_p) = get_positions_with_unknown_orientation(start, next, curve_index);
#if CGAL_MESH_3_PROTECTION_DEBUG & 4
  std::cerr << "positions to determine orientation: " << start_p << " " << next_p << std::endl;
#endif

  if(domain_.is_loop(curve_index))
  {
    // if the curve is a cycle, the direction is the direction passing
    // through the next vertex, and the next-next vertex
    Vertex_handle next_along_curve = next_vertex_along_curve(next,start,curve_index);
    const Bare_point& next_along_curve_wp = get_position(next_along_curve, curve_index);

    return domain_.distance_sign_along_loop(
             start_p, next_p, next_along_curve_wp, curve_index);
  }
  else
  {
    // otherwise, the sign is just the sign of the geodesic distance
    return domain_.distance_sign(start_p, next_p, curve_index);
  }
}

template <typename C3T3, typename MD, typename Sf>
template <typename ErasedVeOutIt>
ErasedVeOutIt
Protect_edges_sizing_field<C3T3, MD, Sf>::
walk_along_edge(const Vertex_handle& start, const Vertex_handle& next,
                Curve_index curve_index,
                const CGAL::Orientation orientation,
                ErasedVeOutIt out) const
{
#if CGAL_MESH_3_PROTECTION_DEBUG & 4
  std::cerr << "walk_along_edge("
            << c3t3_.triangulation().point(start) << ", "
            << c3t3_.triangulation().point(next)
            << ", orientation: " << orientation << ")" << std::endl;
#endif

#if CGAL_MESH_3_PROTECTION_DEBUG & 4
  if(!c3t3_.is_in_complex(start, next))
  {
    std::cerr << "ERROR: the edge (" << c3t3_.triangulation().point(start) << " , "
              << c3t3_.triangulation().point(next) << ") is not in complex!\n";
    dump_c3t3(c3t3_, "dump-bug");
    dump_c3t3_edges(c3t3_, "dump-bug-c3t3");
  }
#endif
  CGAL_precondition(c3t3_.is_in_complex(start, next));

  Vertex_handle previous = start;
  Vertex_handle current = next;

  // Walk along edge since a corner is encountered or the balls of previous
  // and current intersects enough
  while(! is_sampling_dense_enough(previous, current, curve_index, orientation))
  {
    *out++ = current;

    // Don't go through corners
    if(c3t3_.is_in_complex(current) || current == start)
    {
      break;
    }

    // Get next vertex along edge
    Vertex_handle next = next_vertex_along_curve(current,previous,curve_index);
    previous = current;
    current = next;
  }

  return out;
}


template <typename C3T3, typename MD, typename Sf>
typename Protect_edges_sizing_field<C3T3, MD, Sf>::Vertex_handle
Protect_edges_sizing_field<C3T3, MD, Sf>::
next_vertex_along_curve(const Vertex_handle& start,
                        const Vertex_handle& previous,
                        const Curve_index& curve_index) const
{
  CGAL_USE(curve_index);
  CGAL_precondition(c3t3_.curve_index(start, previous) == curve_index);
  CGAL_precondition(domain_.is_loop(curve_index) ||
                     (! c3t3_.is_in_complex(start)));

  Adjacent_vertices adjacent_vertices;
  c3t3_.adjacent_vertices_in_complex(start, std::back_inserter(adjacent_vertices));

  adjacent_vertices.erase
    (std::remove_if(adjacent_vertices.begin(),
                    adjacent_vertices.end(),
                    boost::bind(&Adjacent_vertices::value_type::second, _1) != curve_index),
     adjacent_vertices.end());

//  typename Adjacent_vertices::const_iterator iv = adjacent_vertices.begin();
//  while(iv!=adjacent_vertices.end())
//  {
//    typename Adjacent_vertices::const_iterator iv2 = iv++;
//    if(iv2->second != curve_index)
//      adjacent_vertices.erase(iv2);
//  }

  CGAL_assertion(adjacent_vertices.size() == 2);

  typename Adjacent_vertices::const_iterator ivb = adjacent_vertices.begin();
  if(ivb->first == previous)
  {
    typename Adjacent_vertices::const_iterator ivl = --(adjacent_vertices.end());
    return ivl->first;
  }
  else
  {
    return ivb->first;
  }
}


template <typename C3T3, typename MD, typename Sf>
template <typename InputIterator, typename ErasedVeOutIt>
ErasedVeOutIt
Protect_edges_sizing_field<C3T3, MD, Sf>::
repopulate(InputIterator begin, InputIterator last,
           const Curve_index& curve_index,
           const CGAL::Orientation orientation,
           ErasedVeOutIt out)
{
#if CGAL_MESH_3_PROTECTION_DEBUG & 1
  std::cerr << "repopulate(begin=" << disp_vert(*begin) << "\n"
            << "           last=" << disp_vert(*last)  << "\n"
            << "           distance(begin, last)=" << std::distance(begin, last) << ",\n"
            << "           curve_index=" << CGAL::oformat(curve_index) << ",\n"
            << "           orientation=" << orientation << ")\n";
#endif
  CGAL_assertion(std::distance(begin,last) >= 0);

  // May happen
  if(begin == last)
    return out;

  // Valid because begin < last
  InputIterator current = begin;
  InputIterator previous = current++;

  // Remove edges from c3t3.
  while(current != last)
  {
    c3t3_.remove_from_complex(*previous++, *current++);
  }

  // Remove last edge
  c3t3_.remove_from_complex(*previous, *current);

  // Remove vertices (don't remove the first one and the last one)
  current = begin;
  while(++current != last)
  {
#if CGAL_MESH_3_PROTECTION_DEBUG & 1
    std::cerr << "Removal of ";
    if(is_special(*current)) std::cerr << "SPECIAL ";
    std::cerr << "protecting ball "
              << c3t3_.triangulation().point((*current));
    switch(get_dimension(*current))
    {
    case 0:
      std::cerr << " on corner #";
      break;
    case 1:
      std::cerr << " on curve #";
      break;
    default:
      std::cerr << " ERROR dim=" << get_dimension(*current)  << " curve_index=";
    }
    std::cerr  << CGAL::oformat(c3t3_.index(*current)) << std::endl;
#endif // CGAL_MESH_3_PROTECTION_DEBUG
    *out++ = *current;
    remove_from_correspondence_map(*current, curve_index);
    c3t3_.triangulation().remove(*current);
  }

  // Repopulate edge
  return insert_balls(*begin, *last, curve_index, orientation, out);
}


template <typename C3T3, typename MD, typename Sf>
template <typename InputIterator, typename ErasedVeOutIt>
ErasedVeOutIt
Protect_edges_sizing_field<C3T3, MD, Sf>::
analyze_and_repopulate(InputIterator begin, InputIterator last,
                       const Curve_index& curve_index,
                       const CGAL::Orientation orientation,
                       ErasedVeOutIt out)
{
#if CGAL_MESH_3_PROTECTION_DEBUG & 1
  std::cerr << "analyze_and_repopulate(begin=" << disp_vert(*begin) << "\n"
            << "                       last=" << disp_vert(*last) << "\n"
            << "                       distance(begin, last)=" << std::distance(begin, last) << ",\n"
            << "                       curve_index=" << CGAL::oformat(curve_index) << ",\n"
            << "                       orientation=" << orientation << ")\n";
#endif
  CGAL_assertion(std::distance(begin,last) >= 0);

  // May happen
  if(begin == last)
    return out;

  if(std::distance(begin,last) == 1)
  {
    out = repopulate(begin, last, curve_index, orientation, out);
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
  while(++current != last)
  {
    // Get last element of the stack
    InputIterator previous = ch_stack.top();
    ch_stack.pop();

    // if(prevprev, prev, current) is ok, then go one step forward, i.e. check
    // (prevprevprev, prevprev, current)
    while(!ch_stack.empty() &&
          is_sizing_field_correct(*ch_stack.top(),*previous,*current, curve_index, orientation))
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
  while(!ch_stack.empty())
  {
    InputIterator next = ch_stack.top();
    ch_stack.pop();
    // Iterators are on the reverse order in the stack, thus use [next,current]
    out = repopulate(next, current, curve_index, orientation, out);
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
                        const Curve_index& curve_index,
                        const CGAL::Orientation orientation) const
{
  FT s1 = get_radius(v1);
  FT s2 = get_radius(v2);
  FT s3 = get_radius(v3);

  Bare_point p1, p2, p3;
  boost::tie(p1, p2, p3) = get_positions(v1, v2, v3, curve_index, orientation);

  FT D = domain_.curve_segment_length(p1, p3, curve_index, orientation);
  FT d = domain_.curve_segment_length(p1, p2, curve_index, orientation);

  return (s2 >= (s1 + d/D*(s3-s1)));
}


template <typename C3T3, typename MD, typename Sf>
template <typename ErasedVeOutIt>
ErasedVeOutIt
Protect_edges_sizing_field<C3T3, MD, Sf>::
repopulate_edges_around_corner(const Vertex_handle& v, ErasedVeOutIt out)
{
  CGAL_precondition(c3t3_.is_in_complex(v));

#if CGAL_MESH_3_PROTECTION_DEBUG & 2
  std::cerr << "repopulate_edges_around_corner: " << c3t3_.triangulation().point(v) << std::endl;
#endif

  Adjacent_vertices adjacent_vertices;
  c3t3_.adjacent_vertices_in_complex(v, std::back_inserter(adjacent_vertices));

  for(typename Adjacent_vertices::iterator vit = adjacent_vertices.begin(),
       vend = adjacent_vertices.end() ; vit != vend ; ++vit)
  {
    Vertex_handle next = vit->first;
    const Curve_index curve_index = vit->second;

    // if `v` is incident to a cycle, it might be that the full cycle,
    // including the edge `[next, v]`, has already been processed by
    // `analyze_and_repopulate()` walking in the other direction.
    if(/*domain_.is_loop(curve_index) &&*/ !c3t3_.is_in_complex(v, next))
      continue;

    const CGAL::Orientation orientation = orientation_of_walk(v, next, curve_index);

    // Walk along each incident edge of the corner
    Vertex_vector to_repopulate;
    to_repopulate.push_back(v);
    walk_along_edge(v, next, curve_index, orientation,
                    std::back_inserter(to_repopulate));

    // Return erased vertices, with an exception: if `to_repopulate.back()`
    // is a second corner, it must not be inserted in the output iterator,
    // because otherwise that second corner would be removed from
    // `check_and_repopulate_edges()::vertices` before it is passed itself
    // to `repopulate_edges_around_corner()`.
    if(c3t3_.is_in_complex(to_repopulate.back()))
      std::copy(to_repopulate.begin(), boost::prior(to_repopulate.end()), out);
    else
      std::copy(to_repopulate.begin(), to_repopulate.end(), out);

    // Repopulate
    out = analyze_and_repopulate(to_repopulate.begin(), --to_repopulate.end(),
                                 curve_index, orientation,
                                 out);
  }

  return out;
}

} // namespace Periodic_3_mesh_3

} // namespace CGAL

#endif // CGAL_PERIODIC_3_MESH_3_PROTECT_EDGES_SIZING_FIELD_H
