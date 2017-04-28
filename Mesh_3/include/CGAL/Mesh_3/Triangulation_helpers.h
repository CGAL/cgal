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
  typedef typename Tr::Vertex_handle            Vertex_handle;
  typedef typename Tr::Cell_handle              Cell_handle;
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
                  const Vertex_handle& v,
                  const Weighted_point& p) const;

  /**
   * Returns true if moving \c v to \c p makes no topological
   * change in \c tr
   */
  bool no_topological_change(const Tr& tr,
                             const Vertex_handle& v,
                             const Weighted_point& p,
                             Cell_vector& cells_tos) const;
  bool no_topological_change__without_set_point(
                             const Tr& tr,
                             const Vertex_handle& v,
                             const Weighted_point& p,
                             Cell_vector& cells_tos) const;

  bool no_topological_change(const Tr& tr,
                             const Vertex_handle& v,
                             const Weighted_point& p) const;
  bool no_topological_change__without_set_point(
                             const Tr& tr,
                             const Vertex_handle& v,
                             const Weighted_point& p) const;


  bool inside_protecting_balls(const Tr& tr,
                               const Vertex_handle& v,
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
           const Vertex_handle& v,
           const Weighted_point& p) const
{
  if ( no_topological_change(tr, v, p) )
    v->set_point(p);
  else
  {
    tr.insert(p);
    tr.remove(v);
  }
}

template<typename Tr>
bool
Triangulation_helpers<Tr>::
no_topological_change(const Tr& tr,
                      const Vertex_handle& v0,
                      const Weighted_point& p,
                      Cell_vector& cells_tos) const
{
  bool np = true;
  const Weighted_point fp = v0->point();
  v0->set_point(p);

  if(!well_oriented(tr, cells_tos))
  {
    // Reset (restore) v0
    v0->set_point(fp);
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
        if(tr.side_of_power_sphere(c, cj->vertex(mj)->point(), false)
           != CGAL::ON_UNBOUNDED_SIDE)
        {
          np = false;
          break;
        }
      }
      else
      {
        if(tr.side_of_power_sphere(cj, v1->point(), false)
           != CGAL::ON_UNBOUNDED_SIDE)
        {
          np = false;
          break;
        }
      }
    }
  }

  // Reset (restore) v0
  v0->set_point(fp);

  return np;
}

template<typename Tr>
bool
Triangulation_helpers<Tr>::
no_topological_change__without_set_point(
  const Tr& tr,
  const Vertex_handle& v0,
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
no_topological_change(const Tr& tr,
                      const Vertex_handle& v0,
                      const Weighted_point& p) const
{
  Cell_vector cells_tos;
  cells_tos.reserve(64);
  tr.incident_cells(v0, std::back_inserter(cells_tos));
  return no_topological_change(tr, v0, p, cells_tos);
}

template<typename Tr>
bool
Triangulation_helpers<Tr>::
no_topological_change__without_set_point(
                      const Tr& tr,
                      const Vertex_handle& v0,
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
                        const Vertex_handle& v,
                        const Bare_point& p) const
{
  Vertex_handle nv = tr.nearest_power_vertex(p, v->cell());
  if(nv->point().weight() > 0)
    return tr.geom_traits().compare_squared_distance_3_object()(
          p, nv->point(), nv->point().weight()) != CGAL::LARGER;
  return false;
}

  
/// This function well_oriented is called by no_topological_change after a
/// v->set_point(p)
template<typename Tr>
bool
Triangulation_helpers<Tr>::
well_oriented(const Tr& tr,
              const Cell_vector& cells_tos) const
{
  typedef typename Tr::Geom_traits Gt;
  typename Gt::Orientation_3 orientation = tr.geom_traits().orientation_3_object();
  typename Gt::Construct_point_3 wp2p = tr.geom_traits().construct_point_3_object();

  typename Cell_vector::const_iterator it = cells_tos.begin();
  for( ; it != cells_tos.end() ; ++it)
  {
    Cell_handle c = *it;
    if( tr.is_infinite(c) )
    {
      int iv = c->index(tr.infinite_vertex());
      Cell_handle cj = c->neighbor(iv);
      int mj = tr.mirror_index(c, iv);
      if(orientation(wp2p(cj->vertex(mj)->point()),
                     wp2p(c->vertex((iv+1)&3)->point()),
                     wp2p(c->vertex((iv+2)&3)->point()),
                     wp2p(c->vertex((iv+3)&3)->point())) != CGAL::NEGATIVE)
        return false;
    }
    else if(orientation(wp2p(c->vertex(0)->point()),
                        wp2p(c->vertex(1)->point()),
                        wp2p(c->vertex(2)->point()),
                        wp2p(c->vertex(3)->point())) != CGAL::POSITIVE)
      return false;
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
  typename Gt::Construct_point_3 wp2p = tr.geom_traits().construct_point_3_object();

  typename Cell_vector::const_iterator it = cells_tos.begin();
  for( ; it != cells_tos.end() ; ++it)
  {
    Cell_handle c = *it;
    if( tr.is_infinite(c) )
    {
      int iv = c->index(tr.infinite_vertex());
      Cell_handle cj = c->neighbor(iv);
      int mj = tr.mirror_index(c, iv);
      if(orientation(wp2p(pg(cj->vertex(mj))),
                     wp2p(pg(c->vertex((iv+1)&3))),
                     wp2p(pg(c->vertex((iv+2)&3))),
                     wp2p(pg(c->vertex((iv+3)&3)))) != CGAL::NEGATIVE)
        return false;
    }
    else if(orientation(wp2p(pg(c->vertex(0))),
                        wp2p(pg(c->vertex(1))),
                        wp2p(pg(c->vertex(2))),
                        wp2p(pg(c->vertex(3)))) != CGAL::POSITIVE)
      return false;
  }
  return true;
}



} // end namespace Mesh_3

} //namespace CGAL

#endif // CGAL_MESH_3_TRIANGULATION_HELPERS_H
