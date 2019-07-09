// Copyright (c) 2019 GeometryFactory (France).
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
// Author(s)     : Sebastien Loriot,
//                 Mael Rouxel-Labbé

#ifndef CGAL_POLYGON_MESH_PROCESSING_REMOVE_DEGENERACIES_H
#define CGAL_POLYGON_MESH_PROCESSING_REMOVE_DEGENERACIES_H

#include <CGAL/license/Polygon_mesh_processing/repair.h>

#include <CGAL/Polygon_mesh_processing/shape_predicates.h>

#include <CGAL/boost/graph/Euler_operations.h>

#include <boost/graph/graph_traits.hpp>

#include <iostream>
#include <fstream> // @tmp
#include <set>
#include <sstream> // @tmp

namespace CGAL {
namespace Polygon_mesh_processing {
namespace internal {

template <typename TriangleMesh>
bool is_face_incident_to_border(const typename boost::graph_traits<TriangleMesh>::face_descriptor f,
                                const TriangleMesh& tmesh)
{
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor        halfedge_descriptor;

  for(halfedge_descriptor h : CGAL::halfedges_around_face(halfedge(f, tmesh), tmesh))
  {
    if(is_border_edge(h, tmesh))
      return true;
  }

  return false;
}

template <typename TriangleMesh, typename EdgeContainer, typename NamedParameters>
void add_if_badly_shaped(const typename boost::graph_traits<TriangleMesh>::face_descriptor f,
                         TriangleMesh& tmesh,
                         EdgeContainer& edges_to_collapse,
                         EdgeContainer& edges_to_flip,
                         const NamedParameters& np)
{
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor       halfedge_descriptor;

  // @todo parameters
  const double needle_threshold = 4; // longest edge / shortest edge over this ratio ==> needle
  const double cap_threshold = std::cos(160. / 180 * CGAL_PI); // angle over 120° ==> cap

  if(is_face_incident_to_border(f, tmesh))
    return;

  halfedge_descriptor res = CGAL::Polygon_mesh_processing::is_needle_triangle_face(f, tmesh, needle_threshold, np);
  if(res != boost::graph_traits<TriangleMesh>::null_halfedge())
  {
    std::cout << "add new needle: " << res << std::endl;
    edges_to_collapse.insert(edge(res, tmesh));
  }
  else // let's not make it possible to have a face be both a cap and a needle fo(for now)
  {
    res = CGAL::Polygon_mesh_processing::is_cap_triangle_face(f, tmesh, cap_threshold, np);
    if(res != boost::graph_traits<TriangleMesh>::null_halfedge())
    {
      std::cout << "add new cap: " << res << std::endl;
      edges_to_flip.insert(edge(res, tmesh));
    }
  }
}

} // namespace internal

template <typename FaceRange, typename TriangleMesh, typename NamedParameters>
bool remove_almost_degenerate_faces(const FaceRange& face_range,
                                    TriangleMesh& tmesh,
                                    const NamedParameters& np)
{
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor         vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor       halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::edge_descriptor           edge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor           face_descriptor;

  std::set<edge_descriptor> edges_to_collapse;
  std::set<edge_descriptor> edges_to_flip;

  // @todo could probably do something a bit better by looping edges, consider the incident faces
  // f1 / f2 and look at f1 if f1<f2, and the edge is smaller than the two other edges...
  for(face_descriptor f : face_range)
    internal::add_if_badly_shaped(f, tmesh, edges_to_collapse, edges_to_flip, np);

  int iter = 0;
  for(;;)
  {
    std::cout << edges_to_collapse.size() << " needles and " << edges_to_flip.size() << " caps" << std::endl;

    std::ostringstream oss;
    oss << "degen_cleaning_iter_" << iter++ << ".off";
    std::ofstream out(oss.str().c_str());
    out << std::setprecision(17);
    out << tmesh;
    out.close();

    if(edges_to_collapse.empty() && edges_to_flip.empty())
      return true;

    std::set<edge_descriptor> next_edges_to_collapse;
    std::set<edge_descriptor> next_edges_to_flip;

    // treat needles
    for(edge_descriptor e : edges_to_collapse)
    {
      std::cout << "treat needle: " << e << " (" << tmesh.point(source (e, tmesh)) << " --- " << tmesh.point(target(e, tmesh)) << ")" << std::endl;
      if(CGAL::Euler::does_satisfy_link_condition(e, tmesh))
      {
        vertex_descriptor v = Euler::collapse_edge(e, tmesh); // @todo move 'v' to the midpoint?

        edges_to_flip.erase(e);

        // The geometry of all the faces incident to 'v' has changed and so we recompute their badness
        // @fixme nasty complexity, use tags or something...
        for(halfedge_descriptor inc_h : CGAL::halfedges_around_target(v, tmesh))
        {
          if(is_border_edge(inc_h, tmesh))
            continue;

          // since the 'bad' edge of a face incident to 'v' might not be an incident edge
          for(halfedge_descriptor other_h : CGAL::halfedges_around_face(inc_h, tmesh))
          {
            edge_descriptor other_e = edge(other_h, tmesh);

            if(other_e != e) // erasing the current position while looping is undefined behavior
              edges_to_collapse.erase(other_e);
            edges_to_flip.erase(other_e);
          }

          // adding directly to 'edges_to_flip'
          internal::add_if_badly_shaped(face(inc_h, tmesh), tmesh, next_edges_to_collapse, edges_to_flip, np);
        }
      }
      else
      {
        std::cerr << "Warning: uncollapsable edge! " << tmesh.point(source(e, tmesh)) << " --- "
                                                     << tmesh.point(target(e, tmesh)) << std::endl;
      }
    }

    // treat caps
    for(edge_descriptor e : edges_to_flip)
    {
      std::cout << "treat cap: " << e << " (" << tmesh.point(source (e, tmesh)) << " --- " << tmesh.point(target(e, tmesh)) << ")" << std::endl;
      halfedge_descriptor h = halfedge(e, tmesh);

      // condition for the flip to be valid (the edge to be created does not already exist)
      if(!halfedge(target(next(h, tmesh), tmesh),
                   target(next(opposite(h, tmesh), tmesh), tmesh), tmesh).second)
      {
        std::cout << "Flippin!" << std::endl;
        Euler::flip_edge(h, tmesh);

        internal::add_if_badly_shaped(face(h, tmesh), tmesh, next_edges_to_collapse, next_edges_to_flip, np);
        internal::add_if_badly_shaped(face(opposite(h, tmesh), tmesh), tmesh, next_edges_to_collapse, next_edges_to_flip, np);
      }
      else
      {
        std::cerr << "Warning: unflippable edge! " << tmesh.point(source(h, tmesh)) << " --- "
                                                   << tmesh.point(target(h, tmesh)) << std::endl;
      }
    }

    std::swap(edges_to_collapse, next_edges_to_collapse);
    std::swap(edges_to_flip, next_edges_to_flip);
  }

  return false;
}

template <typename FaceRange, typename TriangleMesh>
bool remove_almost_degenerate_faces(const FaceRange& face_range,
                                    TriangleMesh& tmesh)
{
  return remove_almost_degenerate_faces(face_range, tmesh, parameters::all_default());
}

template <typename TriangleMesh, typename CGAL_PMP_NP_TEMPLATE_PARAMETERS>
bool remove_almost_degenerate_faces(TriangleMesh& tmesh,
                                    const CGAL_PMP_NP_CLASS& np)
{
  return remove_almost_degenerate_faces(faces(tmesh), tmesh, np);
}

template<class TriangleMesh>
bool remove_almost_degenerate_faces(TriangleMesh& tmesh)
{
  return remove_almost_degenerate_faces(tmesh, CGAL::parameters::all_default());
}

} // namespace Polygon_mesh_processing
} // namespace CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_REMOVE_DEGENERACIES_H
