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

#ifndef CGAL_MESH_3_TRIANGULATION_HELPERS_H
#define CGAL_MESH_3_TRIANGULATION_HELPERS_H

#include <CGAL/license/Mesh_3.h>


#include <vector>
#include <CGAL/squared_distance_3.h>

namespace CGAL {

namespace Mesh_3 {


template<typename Tr>
class Triangulation_helpers
{
  typedef typename Tr::Geom_traits              Gt;

  typedef typename Tr::Bare_point               Bare_point;
  typedef typename Tr::Weighted_point           Weighted_point;
  typedef typename Gt::Vector_3                 Vector_3;

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
    Point_getter(const Vertex_handle &vh, const Weighted_point&p)
      : m_vh(vh), m_p(p)
    {}

    const Weighted_point& operator()(const Vertex_handle &vh) const
    {
      return (vh == m_vh ? m_p : vh->point());
    }

  private:
    const Vertex_handle m_vh;
    const Weighted_point &m_p;
  };

public:
  /// Constructor / Destructor
  Triangulation_helpers() {}
  ~Triangulation_helpers() {}

  /**
   * Moves point from \c v to \c p.
   */
  void move_point(Tr& tr,
                  const Vertex_handle v,
                  const Vector_3& move,
                  const Weighted_point& p) const;

  /**
   * Returns true if moving \c v to \c p makes no topological
   * change in \c tr
   */
  bool no_topological_change(Tr& tr,
                             const Vertex_handle v,
                             const Vector_3& move,
                             const Weighted_point& p,
                             Cell_vector& cells_tos) const;
  bool no_topological_change__without_set_point(
                             const Tr& tr,
                             const Vertex_handle v,
                             const Weighted_point& p,
                             Cell_vector& cells_tos) const;

  bool no_topological_change(Tr& tr,
                             const Vertex_handle v,
                             const Vector_3& move,
                             const Weighted_point& p) const;
  bool no_topological_change__without_set_point(
                             const Tr& tr,
                             const Vertex_handle v,
                             const Weighted_point& p) const;

  bool inside_protecting_balls(const Tr& tr,
                               const Vertex_handle v,
                               const Bare_point& p) const;

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
void
Triangulation_helpers<Tr>::
move_point(Tr& tr,
           const Vertex_handle v,
           const Vector_3& move,
           const Weighted_point& p) const
{
  if ( no_topological_change(tr, v, move, p) )
    tr.set_point(v, move, p);
  else
  {
    tr.insert(p);
    tr.remove(v);
  }
}

template<typename Tr>
bool
Triangulation_helpers<Tr>::
no_topological_change(Tr& tr,
                      const Vertex_handle v0,
                      const Vector_3& move,
                      const Weighted_point& p,
                      Cell_vector& cells_tos) const
{
  typename Gt::Construct_opposite_vector_3 cov =
      tr.geom_traits().construct_opposite_vector_3_object();

  bool np = true;
  const Weighted_point fp = tr.point(v0);

//#define CGAL_PERIODIC_BACKUP_CHECK
//#define CGAL_PERIODIC_SIDE_OF_DEBUG
//#define CGAL_PERIODIC_CELL_VERBOSE

#ifdef CGAL_PERIODIC_BACKUP_CHECK
  std::vector<Cell> cells_backup;
  for(Cell_iterator cit=tr.all_cells_begin(); cit!=tr.all_cells_end(); ++cit)
    cells_backup.push_back(*cit);

  std::vector<Vertex> vertex_backup;
  for(typename Tr::Vertex_iterator vit=tr.all_vertices_begin();
                                   vit!=tr.all_vertices_end(); ++vit)
    vertex_backup.push_back(*vit);
#endif

#ifdef CGAL_PERIODIC_SIDE_OF_DEBUG
  CGAL_assertion(well_oriented(tr, cells_tos));

  for(Cell_iterator cit=tr.all_cells_begin(); cit!=tr.all_cells_end(); ++cit)
  {
    for(int j=0; j<4; j++)
      cit->reset_visited(j);
  }

  for(Cell_iterator cit=tr.all_cells_begin(); cit!=tr.all_cells_end(); ++cit)
  {
    Cell_handle c = cit;
    for(int j=0; j<4; j++)
    {
      // Treat each facet only once
      if(c->is_facet_visited(j))
        continue;

      // Set facet and its mirrored version as visited
      Cell_handle cj = c->neighbor(j);
      CGAL_assertion(c!=cj);
      int mj = tr.mirror_index(c, j);
      c->set_facet_visited(j);
      cj->set_facet_visited(mj);

      if(tr.is_infinite(c->vertex(j)))
      {
        CGAL_assertion(tr.side_of_power_sphere(c, tr.point(cj, mj), false)
                         != CGAL::ON_BOUNDED_SIDE);
      }
      else
      {
#ifdef CGAL_PERIODIC_CELL_VERBOSE
        std::cout << "--- pre" << std::endl;
        std::cout << "Cj: " << &*cj << std::endl
                  << cj->vertex(0)->point() << " Off: " << cj->offset(0) << " tra: " << tr.point(cj, 0) << std::endl
                  << cj->vertex(1)->point() << " Off: " << cj->offset(1) << " tra: " << tr.point(cj, 1) << std::endl
                  << cj->vertex(2)->point() << " Off: " << cj->offset(2) << " tra: " << tr.point(cj, 2) << std::endl
                  << cj->vertex(3)->point() << " Off: " << cj->offset(3) << " tra: " << tr.point(cj, 3) << std::endl
                  << "fifth: " << tr.point(c, j) << std::endl;
#endif
        if(tr.side_of_power_sphere(cj, c->vertex(j)->point(), false) == CGAL::ON_BOUNDARY)
        {
//          std::cerr << "warning: on boundary" << std::endl;
        }

        CGAL_assertion(tr.side_of_power_sphere(cj, c->vertex(j)->point(), false)
                        != CGAL::ON_BOUNDED_SIDE);
      }
    }
  }
#endif

//  CGAL_assertion(tr.is_valid(true, 1000));

  // move point
  tr.set_point(v0, move, p);

  if(!well_oriented(tr, cells_tos))
  {
    // Reset (restore) v0
    tr.set_point(v0, cov(move), fp);
    return false;
  }

#if 1
  // Reset visited tags of facets
  std::for_each(cells_tos.begin(), cells_tos.end(), Reset_facet_visited());

  for ( typename Cell_vector::iterator cit = cells_tos.begin() ;
                                       cit != cells_tos.end() ;
                                       ++cit )
  {
    Cell_handle c = *cit;
#else
  // Reset visited tags of facets
  for(Cell_iterator cit=tr.all_cells_begin(); cit!=tr.all_cells_end(); ++cit)
  {
    for(int j=0; j<4; j++)
      cit->reset_visited(j);
  }

  for(Cell_iterator cit=tr.all_cells_begin(); cit!=tr.all_cells_end(); ++cit)
  {
    Cell_handle c = cit;
#endif
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
        if(tr.side_of_power_sphere(c, cj->vertex(mj)->point(), false)
           != CGAL::ON_UNBOUNDED_SIDE)
        {
          np = false;
          break;
        }
      }
      else
      {
#ifdef CGAL_PERIODIC_CELL_VERBOSE
        std::cout << "--- inprogress" << std::endl;
        std::cout << "Cj: " << &*cj << std::endl
                  << cj->vertex(0)->point() << " Off: " << cj->offset(0) << " tra: " << tr.point(cj, 0) << std::endl
                  << cj->vertex(1)->point() << " Off: " << cj->offset(1) << " tra: " << tr.point(cj, 1) << std::endl
                  << cj->vertex(2)->point() << " Off: " << cj->offset(2) << " tra: " << tr.point(cj, 2) << std::endl
                  << cj->vertex(3)->point() << " Off: " << cj->offset(3) << " tra: " << tr.point(cj, 3) << std::endl
                  << "fifth: " << tr.point(c, j) << std::endl;
#endif
        if(tr.side_of_power_sphere(cj, c->vertex(j)->point(), false)
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

#ifdef CGAL_PERIODIC_SIDE_OF_DEBUG
  CGAL_assertion(well_oriented(tr, cells_tos));

  for(Cell_iterator cit=tr.all_cells_begin(); cit!=tr.all_cells_end(); ++cit)
  {
    for(int j=0; j<4; j++)
      cit->reset_visited(j);
  }

  for(Cell_iterator cit=tr.all_cells_begin(); cit!=tr.all_cells_end(); ++cit)
  {
    Cell_handle c = cit;
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
        CGAL_assertion(tr.side_of_power_sphere(c, tr.point(cj, mj), false)
                         != CGAL::ON_BOUNDED_SIDE);
      }
      else
      {
#ifdef CGAL_PERIODIC_CELL_VERBOSE
        std::cout << "--- post" << std::endl;
        std::cout << "Cj: " << &*cj << std::endl
                  << cj->vertex(0)->point() << " Off: " << cj->offset(0) << " tra: " << tr.point(cj, 0) << std::endl
                  << cj->vertex(1)->point() << " Off: " << cj->offset(1) << " tra: " << tr.point(cj, 1) << std::endl
                  << cj->vertex(2)->point() << " Off: " << cj->offset(2) << " tra: " << tr.point(cj, 2) << std::endl
                  << cj->vertex(3)->point() << " Off: " << cj->offset(3) << " tra: " << tr.point(cj, 3) << std::endl
                  << "fifth: " << tr.point(c, j) << std::endl;
#endif
        CGAL_assertion(tr.side_of_power_sphere(cj, c->vertex(j)->point(), false)
                        != CGAL::ON_BOUNDED_SIDE);
      }
    }
  }
#endif

  // DEBUG: check that the triangulation is the same as before...
#ifdef CGAL_PERIODIC_BACKUP_CHECK
  typename std::vector<Cell>::const_iterator ccit = cells_backup.begin();
  for(Cell_iterator cit=tr.all_cells_begin(); cit!=tr.all_cells_end(); ++cit)
  {
    const Cell& c1 = *ccit++;
    const Cell& c2 = *cit;

    CGAL_assertion(c1.vertex(0) == c2.vertex(0));
    CGAL_assertion(c1.vertex(1) == c2.vertex(1));
    CGAL_assertion(c1.vertex(2) == c2.vertex(2));
    CGAL_assertion(c1.vertex(3) == c2.vertex(3));

    CGAL_assertion(c1.offset(0) == c2.offset(0));
    CGAL_assertion(c1.offset(1) == c2.offset(1));
    CGAL_assertion(c1.offset(2) == c2.offset(2));
    CGAL_assertion(c1.offset(3) == c2.offset(3));
  }

  typename std::vector<Vertex>::iterator vvit = vertex_backup.begin();
  for(typename Tr::Vertex_iterator vit=tr.all_vertices_begin();
                                   vit!=tr.all_vertices_end(); ++vit)
  {
    CGAL_assertion(vvit->point() == vit->point());
    ++vvit;
  }
#endif

//  CGAL_assertion(tr.is_valid(true, 1000));

  return np;
}

template<typename Tr>
bool
Triangulation_helpers<Tr>::
no_topological_change__without_set_point(
  const Tr& tr,
  const Vertex_handle v0,
  const Weighted_point& p,
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
      if(tr.is_infinite(v1))
      {
        // Build a copy of c, and replace V0 by a temporary vertex (position "p")
        typename Cell_handle::value_type c_copy (*c);
        int i_v0;
        typename Vertex_handle::value_type v;
        if (c_copy.has_vertex(v0, i_v0))
        {
          v.set_point(p);
          c_copy.set_vertex(i_v0,
            Tr::Triangulation_data_structure::Vertex_range::s_iterator_to(v));
        }

        Cell_handle c_copy_h =
          Tr::Triangulation_data_structure::Cell_range::s_iterator_to(c_copy);
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
        int i_v0;
        typename Vertex_handle::value_type v;
        if (cj_copy.has_vertex(v0, i_v0))
        {
          v.set_point(p);
          cj_copy.set_vertex(i_v0,
            Tr::Triangulation_data_structure::Vertex_range::s_iterator_to(v));
        }

        Cell_handle cj_copy_h =
          Tr::Triangulation_data_structure::Cell_range::s_iterator_to(cj_copy);
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
                      const Weighted_point& p) const
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
                      const Weighted_point& p) const
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
  Vertex_handle nv = tr.nearest_power_vertex(p, v->cell());
  if(nv->point().weight() > 0)
  {
    typename Tr::Geom_traits::Construct_point_3 cp = tr.geom_traits().construct_point_3_object();
    const Weighted_point& nvwp = tr.point(nv);
    return (tr.min_squared_distance(p, cp(nvwp)) <= nv->point().weight());
  }

  return false;
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

      const Weighted_point& mjwp = tr.point(cj, mj);
      const Weighted_point& fwp1 = tr.point(c, (iv+1)&3);
      const Weighted_point& fwp2 = tr.point(c, (iv+2)&3);
      const Weighted_point& fwp3 = tr.point(c, (iv+3)&3);

      if(orientation(cp(mjwp), cp(fwp1), cp(fwp2), cp(fwp3)) != CGAL::NEGATIVE)
        return false;
    }
    else
    {
      const Weighted_point& cwp0 = tr.point(c, 0);
      const Weighted_point& cwp1 = tr.point(c, 1);
      const Weighted_point& cwp2 = tr.point(c, 2);
      const Weighted_point& cwp3 = tr.point(c, 3);

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
