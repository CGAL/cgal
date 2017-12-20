// Copyright (c) 2017 GeometryFactory (France).
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
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Sebastien Loriot


#ifndef CGAL_POLYGON_MESH_PROCESSING_INTERNAL_REPAIR_EXTRA_H
#define CGAL_POLYGON_MESH_PROCESSING_INTERNAL_REPAIR_EXTRA_H

#include <CGAL/license/Polygon_mesh_processing/repair.h>

#include <CGAL/Polygon_mesh_processing/stitch_borders.h>
#include <CGAL/box_intersection_d.h>
#include <CGAL/Union_find.h>

#include <boost/unordered_map.hpp>

#ifndef DOXYGEN_RUNNING

namespace CGAL {

namespace Polygon_mesh_processing {

namespace internal {


template <class PM, class Vpm, class Halfedge_multiplicity>
struct Edges_proximity_report{
  typedef typename boost::graph_traits<PM>::halfedge_descriptor halfedge_descriptor;
  typedef std::vector< std::pair<halfedge_descriptor, halfedge_descriptor> > Halfedge_pairs;

  double m_epsilon;
  Vpm m_vpm;
  PM& m_pm;
  Halfedge_multiplicity& m_multiplicity;

  Halfedge_pairs& m_matching_hedges;

  Edges_proximity_report(double epsilon, Vpm vpm, PM& pm,
                         Halfedge_multiplicity& multiplicity,
                         Halfedge_pairs& matching_hedges)
    : m_epsilon( epsilon )
    , m_vpm(vpm)
    , m_pm(pm)
    , m_multiplicity(multiplicity)
    , m_matching_hedges(matching_hedges)
  {}

  // callback functor that reports all truly intersecting triangles
  template<class Box>
  void operator()(const Box* b1_ptr, const Box* b2_ptr) const
  {
    typedef typename boost::property_traits<Vpm>::reference Point_ref;

    const Box& b1 = *b1_ptr;
    const Box& b2 = *b2_ptr;

    halfedge_descriptor h1 = halfedge(b1.info(), m_pm);
    halfedge_descriptor h2 = halfedge(b2.info(), m_pm);
    if ( !is_border(h1, m_pm) ) h1=opposite(h1, m_pm);
    if ( !is_border(h2, m_pm) ) h2=opposite(h2, m_pm);

    Point_ref src1 = get(m_vpm, source(h1, m_pm));
    Point_ref tgt1 = get(m_vpm, target(h1, m_pm));
    Point_ref src2 = get(m_vpm, source(h2, m_pm));
    Point_ref tgt2 = get(m_vpm, target(h2, m_pm));

    if ( squared_distance(src1,tgt2) < m_epsilon * m_epsilon &&
         squared_distance(tgt1,src2) < m_epsilon * m_epsilon &&
         angle(src1, tgt1, tgt2, src2) == ACUTE )
    {
      // candidate for stitching
      m_matching_hedges.push_back( std::make_pair(h1,h2) );
      ++(m_multiplicity.insert(std::make_pair(h1,0)).first->second);
      ++(m_multiplicity.insert(std::make_pair(h2,0)).first->second);
    }
  }
};

} // end of internal namespace

// \todo document me and move me into another file
template <class PM, class Vpm>
bool is_null_edge(typename boost::graph_traits<PM>::halfedge_descriptor h,
                  PM& pm, Vpm vpm)
{
  return get(vpm, source(h, pm)) ==  get(vpm, target(h,pm));
}

template <class PM, class Vpm>
void collect_close_stitchable_boundary_edges(PM& pm,
                                             double epsilon,
                                             Vpm vpm,
                                             std::vector< std::pair<
                                                typename boost::graph_traits<PM>::halfedge_descriptor,
                                                typename boost::graph_traits<PM>::halfedge_descriptor> >& halfedges_to_stitch,
                                             bool merge_points=true)
{
  typedef typename boost::graph_traits<PM>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<PM>::edge_descriptor edge_descriptor;
  typedef typename boost::graph_traits<PM>::vertex_descriptor vertex_descriptor;

  typedef boost::unordered_map<halfedge_descriptor, int> Halfedge_multiplicity;
  typedef std::vector<std::pair<halfedge_descriptor, halfedge_descriptor> > Halfedge_pairs;

  typedef typename Box_intersection_d::Box_with_info_d<double, 3, edge_descriptor> Box;

  typedef Union_find<vertex_descriptor> UF_vertices;
  typedef std::map<vertex_descriptor, typename UF_vertices::handle> Handle_map;

  typedef typename boost::property_traits<Vpm>::reference Point_ref;

  std::vector<Box> boxes;
  BOOST_FOREACH(edge_descriptor ed, edges(pm))
  {
    if (is_border(ed, pm))
    {
      Point_ref src = get(vpm, source(ed, pm));
      Point_ref tgt = get(vpm, target(ed, pm));

      boxes.push_back( Box(
        Bbox_3( src.x()-epsilon/2, src.y()-epsilon/2, src.z()-epsilon/2,
                      src.x()+epsilon/2, src.y()+epsilon/2, src.z()+epsilon/2 )
          +
        Bbox_3( tgt.x()-epsilon/2, tgt.y()-epsilon/2, tgt.z()-epsilon/2,
                      tgt.x()+epsilon/2, tgt.y()+epsilon/2, tgt.z()+epsilon/2 ),
        ed )
      );
    }
  }

  std::vector<Box*> box_ptrs;
  box_ptrs.reserve(boxes.size());
  BOOST_FOREACH(Box& b, boxes)
    box_ptrs.push_back(&b);


  Halfedge_multiplicity multiplicity;
  Halfedge_pairs matching_hedges;

  box_self_intersection_d(box_ptrs.begin(), box_ptrs.end(),
                          internal::Edges_proximity_report<PM, Vpm, Halfedge_multiplicity>(
                          epsilon, vpm, pm, multiplicity, matching_hedges));


  UF_vertices uf_vertices;
  Handle_map handles;

  typedef std::pair<halfedge_descriptor, halfedge_descriptor> Halfedge_pair;
  BOOST_FOREACH(const Halfedge_pair& p, matching_hedges)
  {
    CGAL_assertion(multiplicity.count(p.first)==1 && multiplicity.count(p.second)==1);
    if (multiplicity[p.first]==1 && multiplicity[p.second]==1)
    {
      bool skip=false;
      BOOST_FOREACH(halfedge_descriptor h, halfedges_around_source(p.first, pm))
        if ( get(vpm, target(h, pm)) == get(vpm, source(h, pm)) )
        {
          // ignore that edge
          skip=true;
          break;
        }
      if (skip) continue;

      BOOST_FOREACH(halfedge_descriptor h, halfedges_around_target(p.first, pm))
        if ( get(vpm, target(h, pm)) == get(vpm, source(h, pm)) )
        {
          // ignore that edge
          skip=true;
          break;
        }
      if (skip) continue;

      // put the opposite vertices in the same set
      if (merge_points)
      {
        internal::uf_join_vertices( source(p.first, pm), target(p.second, pm), uf_vertices, handles );
        internal::uf_join_vertices( target(p.first, pm), source(p.second, pm), uf_vertices, handles );
      }

      // set the halfedges for stitching
      halfedges_to_stitch.push_back(p);
    }
  }

  // update location of vertices in the same set (use the position of the master vertex)
  for(typename UF_vertices::iterator it=uf_vertices.begin(),
                                 it_end=uf_vertices.end(); it!=it_end; ++it)
  {
    typename UF_vertices::handle master = uf_vertices.find(it);
    if (master!=it)
    {
      put(vpm, *it, get(vpm, *master));
    }
  }

  //filter pairs for incident degenerate borders
  std::vector<bool> pair_to_remove(halfedges_to_stitch.size(), false);
  std::size_t i=0;
  std::size_t nb_pairs_to_remove=0;
  BOOST_FOREACH(const Halfedge_pair& p, halfedges_to_stitch)
  {
    cpp11::array<halfedge_descriptor, 4> hedges = {{ p.first, p.second, opposite(p.first, pm), opposite(p.second, pm) }};
    bool null_edge_found=false;

    BOOST_FOREACH( halfedge_descriptor h, hedges)
    {
      if ( is_null_edge(next(h, pm), pm, vpm) ||
           is_null_edge(prev(h, pm), pm, vpm) )
      {
        null_edge_found=true;
        break;
      }
    }
    if (null_edge_found)
    {
      pair_to_remove[i]=true;
      ++nb_pairs_to_remove;
    }
    ++i;
  }

  if (nb_pairs_to_remove!=0)
  {
    std::vector<Halfedge_pair> buffer;
    buffer.reserve(halfedges_to_stitch.size()-nb_pairs_to_remove);
    i=0;
    BOOST_FOREACH(const Halfedge_pair& p, halfedges_to_stitch)
    {
      if (!pair_to_remove[i])
        buffer.push_back(p);
      ++i;
    }
    halfedges_to_stitch.swap(buffer);
  }
}

} } // end of CGAL::Polygon_mesh_processing

#endif // !DOXYGEN_RUNNING

#endif // CGAL_POLYGON_MESH_PROCESSING_INTERNAL_REPAIR_EXTRA_H
