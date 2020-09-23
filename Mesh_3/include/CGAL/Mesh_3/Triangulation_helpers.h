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

#ifndef CGAL_MESH_3_TRIANGULATION_HELPERS_H
#define CGAL_MESH_3_TRIANGULATION_HELPERS_H

#include <CGAL/license/Mesh_3.h>

#include <CGAL/enum.h>
#include <CGAL/internal/Has_nested_type_Bare_point.h>
#include <CGAL/Time_stamper.h>

#include <boost/mpl/if.hpp>
#include <boost/mpl/identity.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/unordered_set.hpp>
#include <boost/utility/enable_if.hpp>

#include <algorithm>
#include <iostream>
#include <iterator>
#include <limits>
#include <utility>
#include <vector>

namespace CGAL {

namespace Mesh_3 {

template<typename Tr>
class Triangulation_helpers
{
  typedef typename Tr::Geom_traits              Gt;

  typedef typename Gt::FT                       FT;
  typedef typename Gt::Vector_3                 Vector_3;

  // If `Tr` is not a triangulation that has defined Bare_point,
  // use Point_3 as defined in the traits class.
  typedef typename boost::mpl::eval_if_c<
    CGAL::internal::Has_nested_type_Bare_point<Tr>::value,
    typename CGAL::internal::Bare_point_type<Tr>,
    boost::mpl::identity<typename Gt::Point_3>
  >::type                                       Bare_point;

  // 'Point' is either a bare point or a weighted point, depending on the triangulation.
  // Since 'Triangulation_helpers' can be templated by an unweighted triangulation,
  // this is one of the rare cases where we use 'Point'.
  typedef typename Tr::Point                    Point;

  typedef typename Tr::Vertex                   Vertex;
  typedef typename Tr::Vertex_handle            Vertex_handle;
  typedef typename Tr::Cell                     Cell;
  typedef typename Tr::Cell_handle              Cell_handle;
  typedef typename Tr::Cell_iterator            Cell_iterator;

  typedef std::vector<Cell_handle>              Cell_vector;

  /**
   * A functor to reset visited flag of each facet of a cell
   */
  struct Reset_facet_visited
  {
    void operator()(const Cell_handle& c) const
    {
      for ( int i=0; i<4 ; ++i )
        c->reset_visited(i);
    }
  };

  /**
   * A functor to get the point of a vertex vh, but replacing
   * it by m_p when vh == m_vh
   */
  struct Point_getter
  {
    /// When the requested will be about vh, the returned point will be p
    /// instead of vh->point()
    Point_getter(const Vertex_handle &vh, const Point&p)
      : m_vh(vh), m_p(p)
    {}

    const Point& operator()(const Vertex_handle &vh) const
    {
      return (vh == m_vh ? m_p : vh->point());
    }

  private:
    const Vertex_handle m_vh;
    const Point &m_p;
  };

public:
  /// Constructor / Destructor
  Triangulation_helpers() {}
  ~Triangulation_helpers() {}

  /**
   * Returns true if moving \c v to \c p makes no topological
   * change in \c tr
   */
  bool no_topological_change(Tr& tr,
                             const Vertex_handle v,
                             const Vector_3& move,
                             const Point& p,
                             Cell_vector& cells_tos) const;
  bool no_topological_change__without_set_point(
                             const Tr& tr,
                             const Vertex_handle v,
                             const Point& p,
                             Cell_vector& cells_tos) const;

  bool no_topological_change(Tr& tr,
                             const Vertex_handle v,
                             const Vector_3& move,
                             const Point& p) const;
  bool no_topological_change__without_set_point(
                             const Tr& tr,
                             const Vertex_handle v,
                             const Point& p) const;

  bool inside_protecting_balls(const Tr& tr,
                               const Vertex_handle v,
                               const Bare_point& p) const;

  /**
   * Returns the squared distance from \c vh to its closest vertex
   *
   * \pre `vh` is not the infinite vertex
   */
  template<typename Tag> // Two versions to distinguish using 'Has_visited_for_vertex_extractor'
  FT get_sq_distance_to_closest_vertex(const Tr& tr,
                                       const Vertex_handle& vh,
                                       const Cell_vector& incident_cells,
                                       typename boost::enable_if_c<Tag::value>::type* = nullptr) const;

  // @todo are the two versions really worth it, I can't tell the difference from a time POV...
  template<typename Tag>
  FT get_sq_distance_to_closest_vertex(const Tr& tr,
                                       const Vertex_handle& vh,
                                       const Cell_vector& incident_cells,
                                       typename boost::disable_if_c<Tag::value>::type* = nullptr) const;

private:
  /**
   * Returns true if \c v is well_oriented on each cell of \c cell_tos
   */
  // For sequential version
  bool well_oriented(const Tr& tr,
                     const Cell_vector& cell_tos) const;
  // For parallel version
  bool well_oriented(const Tr& tr,
                     const Cell_vector& cell_tos,
                     const Point_getter& pg) const;
};

template<typename Tr>
bool
Triangulation_helpers<Tr>::
no_topological_change(Tr& tr,
                      const Vertex_handle v0,
                      const Vector_3& move,
                      const Point& p,
                      Cell_vector& cells_tos) const
{
  if(tr.point(v0) == p)
    return true;

  // For periodic triangulations, calling this function actually takes longer than
  // just removing and re-inserting the point directly because the periodic mesh
  // triangulation's side_of_power_sphere() function is somewhat brute-forcy:
  // it has to check for 27 copies of the query point to determine correctly the side.
  // Thus, we simply return 'false' if the triangulation is a periodic triangulation.
  //
  // Note that the function was nevertheless adapted to work with periodic triangulation
  // so this hack can be disabled if one day 'side_of_power_sphere()' is improved.
  if(boost::is_same<typename Tr::Periodic_tag, Tag_true>::value)
    return false;

  typename Gt::Construct_opposite_vector_3 cov =
      tr.geom_traits().construct_opposite_vector_3_object();

  bool np = true;
  const Point fp = tr.point(v0);

  // move the point
  tr.set_point(v0, move, p);

  if(!well_oriented(tr, cells_tos))
  {
    // Reset (restore) v0
    tr.set_point(v0, cov(move), fp);
    return false;
  }

  // Reset visited tags of facets
  std::for_each(cells_tos.begin(), cells_tos.end(), Reset_facet_visited());

  for ( typename Cell_vector::iterator cit = cells_tos.begin() ;
                                       cit != cells_tos.end() ;
                                       ++cit )
  {
    Cell_handle c = *cit;
    for(int j=0; j<4; j++)
    {
      // Treat each facet only once
      if(c->is_facet_visited(j))
        continue;

      // Set facet and its mirrored version as visited
      Cell_handle cj = c->neighbor(j);
      int mj = tr.mirror_index(c, j);
      c->set_facet_visited(j);
      cj->set_facet_visited(mj);

      if(tr.is_infinite(c->vertex(j)))
      {
        if(tr.side_of_power_sphere(c, tr.point(cj->vertex(mj)), false)
           != CGAL::ON_UNBOUNDED_SIDE)
        {
          np = false;
          break;
        }
      }
      else
      {
        if(tr.side_of_power_sphere(cj, tr.point(c->vertex(j)), false)
           != CGAL::ON_UNBOUNDED_SIDE)
        {
          np = false;
          break;
        }
      }
    }
  }

  // Reset (restore) v0
  tr.set_point(v0, cov(move), fp);

  return np;
}

template<typename Tr>
bool
Triangulation_helpers<Tr>::
no_topological_change__without_set_point(
  const Tr& tr,
  const Vertex_handle v0,
  const Point& p,
  Cell_vector& cells_tos) const
{
  bool np = true;

  Point_getter pg(v0, p);

  if(!well_oriented(tr, cells_tos, pg))
  {
    return false;
  }

  // Reset visited tags of facets
  std::for_each(cells_tos.begin(), cells_tos.end(), Reset_facet_visited());

  for ( typename Cell_vector::iterator cit = cells_tos.begin() ;
        cit != cells_tos.end() ;
        ++cit )
  {
    Cell_handle c = *cit;
    for(int j=0; j<4; j++)
    {
      // Treat each facet only once
      if(c->is_facet_visited(j))
        continue;

      // Set facet and it's mirror's one visited
      Cell_handle cj = c->neighbor(j);
      int mj = tr.mirror_index(c, j);
      c->set_facet_visited(j);
      cj->set_facet_visited(mj);

      Vertex_handle v1 = c->vertex(j);
      typedef typename Tr::Triangulation_data_structure TDS;
      typedef typename TDS::Cell_range Cell_range;
      typedef typename TDS::Vertex_range Vertex_range;
      if(tr.is_infinite(v1))
      {
        // Build a copy of c, and replace V0 by a temporary vertex (position "p")
        typename Cell_handle::value_type c_copy (*c);
        Cell_range::Time_stamper_impl::initialize_time_stamp(&c_copy);
        int i_v0;
        typename Vertex_handle::value_type v;
        Vertex_range::Time_stamper_impl::initialize_time_stamp(&v);
        if (c_copy.has_vertex(v0, i_v0))
        {
          v.set_point(p);
          c_copy.set_vertex(i_v0, Vertex_range::s_iterator_to(v));
        }

        Cell_handle c_copy_h = Cell_range::s_iterator_to(c_copy);
        if(tr.side_of_power_sphere(c_copy_h, pg(cj->vertex(mj)), false)
           != CGAL::ON_UNBOUNDED_SIDE)
        {
          np = false;
          break;
        }
      }
      else
      {
        // Build a copy of cj, and replace V0 by a temporary vertex (position "p")
        typename Cell_handle::value_type cj_copy (*cj);
        Cell_range::Time_stamper_impl::initialize_time_stamp(&cj_copy);
        int i_v0;
        typename Vertex_handle::value_type v;
        Vertex_range::Time_stamper_impl::initialize_time_stamp(&v);
        if (cj_copy.has_vertex(v0, i_v0))
        {
          v.set_point(p);
          cj_copy.set_vertex(i_v0, Vertex_range::s_iterator_to(v));
        }

        Cell_handle cj_copy_h = Cell_range::s_iterator_to(cj_copy);
        if(tr.side_of_power_sphere(cj_copy_h, pg(v1), false)
           != CGAL::ON_UNBOUNDED_SIDE)
        {
          np = false;
          break;
        }
      }
    }
  }

  return np;
}


template<typename Tr>
bool
Triangulation_helpers<Tr>::
no_topological_change(Tr& tr,
                      const Vertex_handle v0,
                      const Vector_3& move,
                      const Point& p) const
{
  Cell_vector cells_tos;
  cells_tos.reserve(64);
  tr.incident_cells(v0, std::back_inserter(cells_tos));
  return no_topological_change(tr, v0, move, p, cells_tos);
}

template<typename Tr>
bool
Triangulation_helpers<Tr>::
no_topological_change__without_set_point(
                      const Tr& tr,
                      const Vertex_handle v0,
                      const Point& p) const
{
  Cell_vector cells_tos;
  cells_tos.reserve(64);
  tr.incident_cells(v0, std::back_inserter(cells_tos));
  return no_topological_change__without_set_point(tr, v0, p, cells_tos);
}


template<typename Tr>
bool
Triangulation_helpers<Tr>::
inside_protecting_balls(const Tr& tr,
                        const Vertex_handle v,
                        const Bare_point& p) const
{
  typename Gt::Compare_weighted_squared_radius_3 cwsr =
    tr.geom_traits().compare_weighted_squared_radius_3_object();

  Vertex_handle nv = tr.nearest_power_vertex(p, v->cell());
  const Point& nvwp = tr.point(nv);
  if(cwsr(nvwp, FT(0)) == CGAL::SMALLER)
  {
    typename Tr::Geom_traits::Construct_point_3 cp = tr.geom_traits().construct_point_3_object();
    const Point& nvwp = tr.point(nv);
    // 'true' if the distance between 'p' and 'nv' is smaller or equal than the weight of 'nv'
    return (cwsr(nvwp , - tr.min_squared_distance(p, cp(nvwp))) != CGAL::LARGER);
  }

  return false;
}

/// Return the squared distance from vh to its closest vertex
/// if `Has_visited_for_vertex_extractor` is `true`
template<typename Tr>
template<typename Tag>
typename Triangulation_helpers<Tr>::FT
Triangulation_helpers<Tr>::
get_sq_distance_to_closest_vertex(const Tr& tr,
                                  const Vertex_handle& vh,
                                  const Cell_vector& incident_cells,
                                  typename boost::enable_if_c<Tag::value>::type*) const
{
  CGAL_precondition(!tr.is_infinite(vh));

  typedef std::vector<Vertex_handle>              Vertex_container;

  // There is no need to use tr.min_squared_distance() here because we are computing
  // distances between 'v' and a neighbor within their common cell, which means
  // that even if we are using a periodic triangulation, the distance is correctly computed.
  typename Gt::Compute_squared_distance_3 csqd = tr.geom_traits().compute_squared_distance_3_object();
  typename Gt::Construct_point_3 cp = tr.geom_traits().construct_point_3_object();

  Vertex_container treated_vertices;
  FT min_sq_dist = std::numeric_limits<FT>::infinity();

  for(typename Cell_vector::const_iterator cit = incident_cells.begin();
                                           cit != incident_cells.end(); ++cit)
  {
    const Cell_handle c = (*cit);
    const int k = (*cit)->index(vh);
    const Point& wpvh = tr.point(c, k);

    // For each vertex of the cell
    for(int i=1; i<4; ++i)
    {
      const int n = (k+i)&3;
      const Vertex_handle& vn = c->vertex(n);

      if(vn == Vertex_handle() ||
         tr.is_infinite(vn) ||
         vn->visited_for_vertex_extractor)
        continue;

      vn->visited_for_vertex_extractor = true;
      treated_vertices.push_back(vn);

      const Point& wpvn = tr.point(c, n);
      const FT sq_d = csqd(cp(wpvh), cp(wpvn));

      if(sq_d < min_sq_dist)
        min_sq_dist = sq_d;
    }
  }

  for(std::size_t i=0; i < treated_vertices.size(); ++i)
    treated_vertices[i]->visited_for_vertex_extractor = false;

  return min_sq_dist;
}

/// Return the squared distance from vh to its closest vertex
/// if `Has_visited_for_vertex_extractor` is `false`
template<typename Tr>
template<typename Tag>
typename Triangulation_helpers<Tr>::FT
Triangulation_helpers<Tr>::
get_sq_distance_to_closest_vertex(const Tr& tr,
                                  const Vertex_handle& vh,
                                  const Cell_vector& incident_cells,
                                  typename boost::disable_if_c<Tag::value>::type*) const
{
  CGAL_precondition(!tr.is_infinite(vh));

  typedef CGAL::Hash_handles_with_or_without_timestamps      Hash_fct;
  typedef boost::unordered_set<Vertex_handle, Hash_fct>      Vertex_container;
  typedef typename Vertex_container::iterator                VC_it;

  // There is no need to use tr.min_squared_distance() here because we are computing
  // distances between 'v' and a neighbor within their common cell, which means
  // that even if we are using a periodic triangulation, the distance is correctly computed.
  typename Gt::Compute_squared_distance_3 csqd = tr.geom_traits().compute_squared_distance_3_object();
  typename Gt::Construct_point_3 cp = tr.geom_traits().construct_point_3_object();

  Vertex_container treated_vertices;
  FT min_sq_dist = std::numeric_limits<FT>::infinity();

  for(typename Cell_vector::const_iterator cit = incident_cells.begin();
                                           cit != incident_cells.end(); ++cit)
  {
    const Cell_handle c = (*cit);
    const int k = (*cit)->index(vh);
    const Point& wpvh = tr.point(c, k);

    // For each vertex of the cell
    for(int i=1; i<4; ++i)
    {
      const int n = (k+i)&3;
      const Vertex_handle& vn = c->vertex(n);

      if(vn == Vertex_handle() ||
         tr.is_infinite(vn))
        continue;

      std::pair<VC_it, bool> is_insert_successful = treated_vertices.insert(vn);
      if(! is_insert_successful.second) // vertex has already been treated
        continue;

      const Point& wpvn = tr.point(c, n);
      const FT sq_d = csqd(cp(wpvh), cp(wpvn));

      if(sq_d < min_sq_dist)
        min_sq_dist = sq_d;
    }
  }

  return min_sq_dist;
}


/// This function well_oriented is called by no_topological_change after the
/// position of the vertex has been (tentatively) modified.
template<typename Tr>
bool
Triangulation_helpers<Tr>::
well_oriented(const Tr& tr,
              const Cell_vector& cells_tos) const
{
  typedef typename Tr::Geom_traits Gt;
  typename Gt::Orientation_3 orientation = tr.geom_traits().orientation_3_object();
  typename Gt::Construct_point_3 cp = tr.geom_traits().construct_point_3_object();

  typename Cell_vector::const_iterator it = cells_tos.begin();
  for( ; it != cells_tos.end() ; ++it)
  {
    Cell_handle c = *it;
    if( tr.is_infinite(c) )
    {
      int iv = c->index(tr.infinite_vertex());
      Cell_handle cj = c->neighbor(iv);
      int mj = tr.mirror_index(c, iv);

      const Point& mjwp = tr.point(cj, mj);
      const Point& fwp1 = tr.point(c, (iv+1)&3);
      const Point& fwp2 = tr.point(c, (iv+2)&3);
      const Point& fwp3 = tr.point(c, (iv+3)&3);

      if(orientation(cp(mjwp), cp(fwp1), cp(fwp2), cp(fwp3)) != CGAL::NEGATIVE)
        return false;
    }
    else
    {
      const Point& cwp0 = tr.point(c, 0);
      const Point& cwp1 = tr.point(c, 1);
      const Point& cwp2 = tr.point(c, 2);
      const Point& cwp3 = tr.point(c, 3);

      if(orientation(cp(cwp0), cp(cwp1), cp(cwp2), cp(cwp3)) != CGAL::POSITIVE)
      return false;
    }
  }
  return true;
}

/// Another version for the parallel version
/// Here, the set_point is not done before, but we use a Point_getter instance
/// to get the point of a vertex.
template<typename Tr>
bool
Triangulation_helpers<Tr>::
well_oriented(const Tr& tr,
              const Cell_vector& cells_tos,
              const Point_getter& pg) const
{
  typedef typename Tr::Geom_traits Gt;
  typename Gt::Orientation_3 orientation = tr.geom_traits().orientation_3_object();
  typename Gt::Construct_point_3 cp = tr.geom_traits().construct_point_3_object();

  typename Cell_vector::const_iterator it = cells_tos.begin();
  for( ; it != cells_tos.end() ; ++it)
  {
    Cell_handle c = *it;
    if( tr.is_infinite(c) )
    {
      int iv = c->index(tr.infinite_vertex());
      Cell_handle cj = c->neighbor(iv);
      int mj = tr.mirror_index(c, iv);
      if(orientation(cp(pg(cj->vertex(mj))),
                     cp(pg(c->vertex((iv+1)&3))),
                     cp(pg(c->vertex((iv+2)&3))),
                     cp(pg(c->vertex((iv+3)&3)))) != CGAL::NEGATIVE)
        return false;
    }
    else if(orientation(cp(pg(c->vertex(0))),
                        cp(pg(c->vertex(1))),
                        cp(pg(c->vertex(2))),
                        cp(pg(c->vertex(3)))) != CGAL::POSITIVE)
      return false;
  }
  return true;
}



} // end namespace Mesh_3

} //namespace CGAL

#endif // CGAL_MESH_3_TRIANGULATION_HELPERS_H
