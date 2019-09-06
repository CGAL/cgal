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
#include <set>
#ifdef CGAL_PMP_DEBUG_REMOVE_DEGENERACIES
#include <sstream>
#include <fstream>
#endif

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

template <typename TriangleMesh, typename NamedParameters>
std::array<typename boost::graph_traits<TriangleMesh>::halfedge_descriptor, 2>
is_badly_shaped(const typename boost::graph_traits<TriangleMesh>::face_descriptor f,
                      TriangleMesh& tmesh,
                      const NamedParameters& np)
{
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor       halfedge_descriptor;
  halfedge_descriptor null_halfedge = boost::graph_traits<TriangleMesh>::null_halfedge();

  // @todo parameters
  const double needle_threshold = 4; // longest edge / shortest edge over this ratio ==> needle
  const double cap_threshold = std::cos(160. / 180 * CGAL_PI); // angle over 120° ==> cap

  halfedge_descriptor res = CGAL::Polygon_mesh_processing::is_needle_triangle_face(f, tmesh, needle_threshold, np);
  if(res != boost::graph_traits<TriangleMesh>::null_halfedge())
  {
    return make_array(res, null_halfedge);
  }
  else // let's not make it possible to have a face be both a cap and a needle (for now)
  {
    res = CGAL::Polygon_mesh_processing::is_cap_triangle_face(f, tmesh, cap_threshold, np);
    if(res != boost::graph_traits<TriangleMesh>::null_halfedge())
    {
      return make_array(null_halfedge, res);
    }
  }
  return make_array(null_halfedge, null_halfedge);
}

template <typename TriangleMesh, typename EdgeContainer, typename NamedParameters>
void add_if_badly_shaped(const typename boost::graph_traits<TriangleMesh>::face_descriptor f,
                         TriangleMesh& tmesh,
                         EdgeContainer& edges_to_collapse,
                         EdgeContainer& edges_to_flip,
                         const NamedParameters& np)
{
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor       halfedge_descriptor;

  std::array<halfedge_descriptor, 2> res = is_badly_shaped(f, tmesh, np);

  if(res[0] != boost::graph_traits<TriangleMesh>::null_halfedge())
  {
#ifdef  CGAL_PMP_DEBUG_REMOVE_DEGENERACIES
    std::cout << "add new needle: " << res << std::endl;
#endif
    edges_to_collapse.insert(edge(res[0], tmesh));
  }
  else // let's not make it possible to have a face be both a cap and a needle (for now)
  {
    if(res[1] != boost::graph_traits<TriangleMesh>::null_halfedge())
    {
#ifdef  CGAL_PMP_DEBUG_REMOVE_DEGENERACIES
      std::cout << "add new cap: " << res << std::endl;
#endif
      edges_to_flip.insert(edge(res[1], tmesh));
    }
  }
}

} // namespace internal

template <typename FaceRange, typename TriangleMesh, typename NamedParameters>
bool remove_almost_degenerate_faces(const FaceRange& face_range,
                                    TriangleMesh& tmesh,
                                    const NamedParameters& np)
{
//  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor         vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor       halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::edge_descriptor           edge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor           face_descriptor;

  std::set<edge_descriptor> edges_to_collapse;
  std::set<edge_descriptor> edges_to_flip;

  // @todo could probably do something a bit better by looping edges, consider the incident faces
  // f1 / f2 and look at f1 if f1<f2, and the edge is smaller than the two other edges...
  for(face_descriptor f : face_range)
    internal::add_if_badly_shaped(f, tmesh, edges_to_collapse, edges_to_flip, np);

#ifdef  CGAL_PMP_DEBUG_REMOVE_DEGENERACIES
  int iter = 0;
#endif
  for(;;)
  {
#ifdef  CGAL_PMP_DEBUG_REMOVE_DEGENERACIES
    std::cout << edges_to_collapse.size() << " needles and " << edges_to_flip.size() << " caps" << std::endl;
    std::ostringstream oss;
    oss << "degen_cleaning_iter_" << iter++ << ".off";
    std::ofstream out(oss.str().c_str());
    out << std::setprecision(17);
    out << tmesh;
    out.close();
#endif

    if(edges_to_collapse.empty() && edges_to_flip.empty())
      return true;

    std::set<edge_descriptor> next_edges_to_collapse;
    std::set<edge_descriptor> next_edges_to_flip;

    // treat needles
    for(edge_descriptor e : edges_to_collapse)
    {
#ifdef  CGAL_PMP_DEBUG_REMOVE_DEGENERACIES
      std::cout << "treat needle: " << e << " (" << tmesh.point(source (e, tmesh)) << " --- " << tmesh.point(target(e, tmesh)) << ")" << std::endl;
#endif
      if(CGAL::Euler::does_satisfy_link_condition(e, tmesh))
      {
        // the following edges are removed by the collapse
        halfedge_descriptor h = halfedge(e, tmesh);
        if (!is_border(h, tmesh))
        {
          edge_descriptor pe=edge(prev(h,tmesh), tmesh);
          edges_to_flip.erase(pe);
          edges_to_collapse.erase(pe);
        }
        h = opposite(h, tmesh);
        if (!is_border(h, tmesh))
        {
          edge_descriptor pe=edge(prev(h,tmesh), tmesh);
          edges_to_flip.erase(pe);
          edges_to_collapse.erase(pe);
        }
        edges_to_flip.erase(e);

        /* vertex_descriptor v = */ Euler::collapse_edge(e, tmesh); // @todo move 'v' to the midpoint?
      }
#ifdef  CGAL_PMP_DEBUG_REMOVE_DEGENERACIES
      else
      {
        std::cerr << "Warning: uncollapsable edge! " << tmesh.point(source(e, tmesh)) << " --- "
                                                     << tmesh.point(target(e, tmesh)) << std::endl;
      }
#endif
    }

    // treat caps
    for(edge_descriptor e : edges_to_flip)
    {
#ifdef  CGAL_PMP_DEBUG_REMOVE_DEGENERACIES
      std::cout << "treat cap: " << e << " (" << tmesh.point(source (e, tmesh)) << " --- " << tmesh.point(target(e, tmesh)) << ")" << std::endl;
#endif

      // special case on the border
      if ( is_border(e, tmesh) )
      {
        // remove the triangle
        halfedge_descriptor h = halfedge(e, tmesh);
        if (is_border(h, tmesh))
          h = opposite(h, tmesh);
        edges_to_flip.erase(edge(prev(h, tmesh), tmesh));
        edges_to_flip.erase(edge(next(h, tmesh), tmesh));
        Euler::remove_face(h, tmesh);
        continue;
      }

      halfedge_descriptor h = halfedge(e, tmesh);

      // condition for the flip to be valid (the edge to be created does not already exist)
      if(!halfedge(target(next(h, tmesh), tmesh),
                   target(next(opposite(h, tmesh), tmesh), tmesh), tmesh).second)
      {
#ifdef  CGAL_PMP_DEBUG_REMOVE_DEGENERACIES
        std::cout << "Flipping" << std::endl;
#endif
        Euler::flip_edge(h, tmesh);

        // remove edges from the new faces
        edges_to_flip.erase(edge(prev(h, tmesh), tmesh));
        edges_to_flip.erase(edge(next(h, tmesh), tmesh));
        edges_to_flip.erase(edge(prev(opposite(h, tmesh), tmesh), tmesh));
        edges_to_flip.erase(edge(next(opposite(h, tmesh), tmesh), tmesh));

        // make sure we do not enter in an infinite loop
        if (!is_border(h, tmesh))
        {
          std::array<halfedge_descriptor,2> hp = internal::is_badly_shaped(face(h, tmesh), tmesh, np);
          if (hp[1]!=boost::graph_traits<TriangleMesh>::null_halfedge() &&
              edge(hp[1], tmesh) != edge(h, tmesh)) // avoid infinite loop
          {
            next_edges_to_flip.insert( edge(hp[1], tmesh) );
          }
        }
        h=opposite(h, tmesh);
        if (!is_border(h, tmesh))
        {
          std::array<halfedge_descriptor,2> hp = internal::is_badly_shaped(face(h, tmesh), tmesh, np);
          if (hp[1]!=boost::graph_traits<TriangleMesh>::null_halfedge() &&
              edge(hp[1], tmesh) != edge(h, tmesh)) // avoid infinite loop
          {
            next_edges_to_flip.insert( edge(hp[1], tmesh) );
          }
        }
      }
#ifdef  CGAL_PMP_DEBUG_REMOVE_DEGENERACIES
      else
      {
        std::cerr << "Warning: unflippable edge! " << tmesh.point(source(h, tmesh)) << " --- "
                                                   << tmesh.point(target(h, tmesh)) << std::endl;
      }
#endif
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
