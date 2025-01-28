// Copyright (c) 2017 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sebastien Loriot


#ifndef CGAL_POLYGON_MESH_PROCESSING_INTERNAL_REPAIR_EXTRA_H
#define CGAL_POLYGON_MESH_PROCESSING_INTERNAL_REPAIR_EXTRA_H

#include <CGAL/license/Polygon_mesh_processing/geometric_repair.h>

#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/box_intersection_d.h>
#include <CGAL/Kernel/global_functions_3.h>
#include <CGAL/Polygon_mesh_processing/stitch_borders.h>
#include <CGAL/Union_find.h>

#include <unordered_map>

#include <array>
#include <utility>
#include <vector>

#ifndef DOXYGEN_RUNNING

namespace CGAL {

namespace Polygon_mesh_processing {

namespace internal {


template <class PM, class Vpm, class Halfedge_multiplicity>
struct Edges_proximity_report
{
  typedef typename boost::graph_traits<PM>::halfedge_descriptor halfedge_descriptor;
  typedef std::vector< std::pair<halfedge_descriptor, halfedge_descriptor> > Halfedge_pairs;

  typedef typename GetGeomTraits<PM>::type::FT FT;

  FT m_sq_epsilon;
  Vpm m_vpm;
  PM& m_pm;
  Halfedge_multiplicity& m_multiplicity;

  Halfedge_pairs& m_matching_hedges;

  Edges_proximity_report(double epsilon, Vpm vpm, PM& pm,
                         Halfedge_multiplicity& multiplicity,
                         Halfedge_pairs& matching_hedges)
    : m_sq_epsilon(square(epsilon))
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

    if ( compare_squared_distance(src1, tgt2, m_sq_epsilon) == SMALLER &&
         compare_squared_distance(tgt1, src2, m_sq_epsilon) == SMALLER &&
         angle(src1, tgt1, tgt2, src2) == ACUTE )
    {
      // candidate for stitching
      m_matching_hedges.emplace_back(h1,h2);
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

  typedef std::unordered_map<halfedge_descriptor, int> Halfedge_multiplicity;
  typedef std::vector<std::pair<halfedge_descriptor, halfedge_descriptor> > Halfedge_pairs;

  typedef CGAL::Box_intersection_d::ID_FROM_BOX_ADDRESS Box_policy;
  typedef CGAL::Box_intersection_d::Box_with_info_d<double, 3, edge_descriptor, Box_policy> Box;

  typedef Union_find<vertex_descriptor> UF_vertices;
  typedef std::map<vertex_descriptor, typename UF_vertices::handle> Handle_map;

  typedef typename boost::property_traits<Vpm>::reference Point_ref;

  const double half_eps = 0.5 * epsilon;

  std::vector<Box> boxes;
  for(edge_descriptor ed : edges(pm))
  {
    if (is_border(ed, pm))
    {
      Point_ref src = get(vpm, source(ed, pm));
      Point_ref tgt = get(vpm, target(ed, pm));

      const double sx = to_double(src.x());
      const double sy = to_double(src.y());
      const double sz = to_double(src.z());
      const double tx = to_double(tgt.x());
      const double ty = to_double(tgt.y());
      const double tz = to_double(tgt.z());

      boxes.emplace_back(Bbox_3(sx - half_eps, sy - half_eps, sz - half_eps,
                                sx + half_eps, sy + half_eps, sz + half_eps) +
                         Bbox_3(tx - half_eps, ty - half_eps, tz - half_eps,
                                tx + half_eps, ty + half_eps, tz + half_eps),
                         ed);
    }
  }

  std::vector<Box*> box_ptrs;
  box_ptrs.reserve(boxes.size());
  for(Box& b : boxes)
    box_ptrs.push_back(&b);

  Halfedge_multiplicity multiplicity;
  Halfedge_pairs matching_hedges;

  box_self_intersection_d(box_ptrs.begin(), box_ptrs.end(),
                          internal::Edges_proximity_report<PM, Vpm, Halfedge_multiplicity>(
                          epsilon, vpm, pm, multiplicity, matching_hedges));


  UF_vertices uf_vertices;
  Handle_map handles;

  typedef std::pair<halfedge_descriptor, halfedge_descriptor> Halfedge_pair;
  for(const Halfedge_pair& p : matching_hedges)
  {
    CGAL_assertion(multiplicity.count(p.first)==1 && multiplicity.count(p.second)==1);
    if (multiplicity[p.first]==1 && multiplicity[p.second]==1)
    {
      bool skip=false;
      for(halfedge_descriptor h : halfedges_around_source(p.first, pm))
        if ( get(vpm, target(h, pm)) == get(vpm, source(h, pm)) )
        {
          // ignore that edge
          skip=true;
          break;
        }
      if (skip) continue;

      for(halfedge_descriptor h : halfedges_around_target(p.first, pm))
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
  for(const Halfedge_pair& p : halfedges_to_stitch)
  {
    std::array<halfedge_descriptor, 4> hedges = {{ p.first, p.second, opposite(p.first, pm), opposite(p.second, pm) }};
    bool null_edge_found=false;

    for(halfedge_descriptor h : hedges)
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
    for(const Halfedge_pair& p : halfedges_to_stitch)
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
