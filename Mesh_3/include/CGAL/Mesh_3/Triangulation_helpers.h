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

#include <vector>

namespace CGAL {

namespace Mesh_3 {
  
  
template<typename Tr>
class Triangulation_helpers
{
  typedef typename Tr::Geom_traits              Gt;
  typedef typename Gt::Point_3                  Point_3;
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
  
public:
  /// Constructor / Destructor
  Triangulation_helpers() {}
  ~Triangulation_helpers() {}
  
  /**
   * Moves point from \c v to \c p.
   */
  void move_point(Tr& tr,
                  const Vertex_handle& v,
                  const Point_3& p) const;
  
  /**
   * Returns true if moving \c v to \c p makes no topological
   * change in \c tr
   */
  bool no_topological_change(const Tr& tr,
                             const Vertex_handle& v,
                             const Point_3& p) const;
  
private:
  /**
   * Returns true if \c v is well_oriented on each cell of \c cell_tos
   */
  bool well_oriented(const Tr& tr,
                     const Cell_vector& cell_tos) const;
};
  
  
template<typename Tr>
void
Triangulation_helpers<Tr>::
move_point(Tr& tr,
           const Vertex_handle& v,
           const Point_3& p) const
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
                      const Point_3& p) const
{
  bool np = true;
  Point_3 fp = v0->point();
  v0->set_point(p);
  
  Cell_vector cells_tos;
  tr.incident_cells(v0, std::back_inserter(cells_tos));
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
well_oriented(const Tr& tr,
              const Cell_vector& cells_tos) const
{
  typename Cell_vector::const_iterator it = cells_tos.begin();
  for( ; it != cells_tos.end() ; ++it)
  {
    Cell_handle c = *it;
    if( tr.is_infinite(c) ) 
    {
      int iv = c->index(tr.infinite_vertex());
      Cell_handle cj = c->neighbor(iv);
      int mj = tr.mirror_index(c, iv);
      if(CGAL::orientation(cj->vertex(mj)->point(),
                           c->vertex((iv+1)&3)->point(),
                           c->vertex((iv+2)&3)->point(),
                           c->vertex((iv+3)&3)->point()) != CGAL::NEGATIVE)
        return false;
    }
    else if(CGAL::orientation(c->vertex(0)->point(),
                              c->vertex(1)->point(),
                              c->vertex(2)->point(),
                              c->vertex(3)->point()) != CGAL::POSITIVE)
      return false;
  }
  return true;
} 



} // end namespace Mesh_3 
  
} //namespace CGAL

#endif // CGAL_MESH_3_TRIANGULATION_HELPERS_H
