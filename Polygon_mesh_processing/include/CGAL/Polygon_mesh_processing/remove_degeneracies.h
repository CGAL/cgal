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
#include <CGAL/Polygon_mesh_processing/measure.h>

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
  const halfedge_descriptor null_halfedge = boost::graph_traits<TriangleMesh>::null_halfedge();

  // @todo parameters
  const double needle_threshold = 4; // longest edge / shortest edge over this ratio ==> needle
  const double cap_threshold = std::cos(160. / 180 * CGAL_PI); // angle over 160° ==> cap
  const double collapse_threshold = 0.2; // max length of edges allowed to be collapsed

  halfedge_descriptor res = CGAL::Polygon_mesh_processing::is_needle_triangle_face(f, tmesh, needle_threshold, np);
  if(res != boost::graph_traits<TriangleMesh>::null_halfedge())
  {
    if ( edge_length(res, tmesh, np) <= collapse_threshold )
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
void collect_badly_shaped_triangles(const typename boost::graph_traits<TriangleMesh>::face_descriptor f,
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
    std::cout << "add new needle: " << edge(res[0], tmesh) << std::endl;
#endif
    edges_to_collapse.insert(edge(res[0], tmesh));
  }
  else // let's not make it possible to have a face be both a cap and a needle (for now)
  {
    if(res[1] != boost::graph_traits<TriangleMesh>::null_halfedge())
    {
#ifdef  CGAL_PMP_DEBUG_REMOVE_DEGENERACIES
      std::cout << "add new cap: " << edge(res[1],tmesh) << std::endl;
#endif
      edges_to_flip.insert(edge(res[1], tmesh));
    }
  }
}

// Following Ronfard et al. 96 we look at variation of the normal after the collapse
// the collapse must be topologically valid
template <class TriangleMesh, class NamedParameters>
bool is_collapse_geometrically_valid(typename boost::graph_traits<TriangleMesh>::halfedge_descriptor h,
                                     const TriangleMesh& tmesh,
                                     const NamedParameters& np)
{
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters>::const_type VPM;
  typedef typename boost::property_traits<VPM>::reference Point_ref;
  typedef typename GetGeomTraits<TriangleMesh, NamedParameters>::type Traits;

  VPM vpm = choose_param(get_param(np, internal_np::vertex_point),
                         get_const_property_map(vertex_point, tmesh));


/// @todo handle boundary edges

  h = opposite(h, tmesh); // Euler::collapse edge keeps the target and removes the source

  // source is kept, target is removed
  CGAL_assertion( target(h, tmesh)==vertex_removed );
  Point_ref kept = get(vpm, source(h, tmesh));
  Point_ref removed= get(vpm, target(h, tmesh));

  // consider triangles incident to the vertex removed
  halfedge_descriptor stop = prev(opposite(h, tmesh), tmesh);
  halfedge_descriptor hi = opposite(next(h, tmesh), tmesh);

  std::vector<halfedge_descriptor> triangles;
  while (hi!=stop)
  {
    if (!is_border(hi, tmesh))
    {
      Point_ref a = get(vpm, target(next(hi, tmesh), tmesh));
      Point_ref b = get(vpm, source(hi, tmesh));

      //ack a-b-point_remove and a-b-point_kept has a compatible orientation
      /// @todo use a predicate
      typename Traits::Vector_3 n1 = CGAL::cross_product(removed-a, b-a);
      typename Traits::Vector_3 n2 = CGAL::cross_product(kept-a, b-a);
      if ( n1*n2 <=0 )
        return false;
    }
    hi = opposite(next(hi, tmesh), tmesh);
  }

  return true;
}

} // namespace internal

// @todo check what to use as priority queue with removable elements, set might not be optimal
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
    internal::collect_badly_shaped_triangles(f, tmesh, edges_to_collapse, edges_to_flip, np);

#ifdef  CGAL_PMP_DEBUG_REMOVE_DEGENERACIES
  int iter = 0;
#endif

  for(;;)
  {
    bool something_was_done = false;

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

    /// @todo maybe using a priority queue handling the more almost degenerate elements should be used
    std::set<edge_descriptor> next_edges_to_collapse;
    std::set<edge_descriptor> next_edges_to_flip;

    // treat needles
#ifdef  CGAL_PMP_DEBUG_REMOVE_DEGENERACIES_EXTRA
    int kk=0;
    std::ofstream(std::string("tmp/n-00000.off")) << tmesh;
#endif
    while (!edges_to_collapse.empty())
    {
      edge_descriptor e = *edges_to_collapse.begin();
      edges_to_collapse.erase(edges_to_collapse.begin());

#ifdef  CGAL_PMP_DEBUG_REMOVE_DEGENERACIES
      std::cout << "  treat needle: " << e << " (" << tmesh.point(source (e, tmesh)) << " --- " << tmesh.point(target(e, tmesh)) << ")" << std::endl;
#endif
      if(CGAL::Euler::does_satisfy_link_condition(e, tmesh))
      {
        // the following edges are removed by the collapse
        halfedge_descriptor h = halfedge(e, tmesh);
        CGAL_assertion(!is_border(h, tmesh)); // because extracted from a face

        std::array<halfedge_descriptor,2> nc = internal::is_badly_shaped(face(h, tmesh), tmesh, np);
        if (nc[0]!=h)
        {
#ifdef  CGAL_PMP_DEBUG_REMOVE_DEGENERACIES
          std::cerr << "Warning: Needle criteria no longer verified " << tmesh.point(source(e, tmesh)) << " "
                                                                      << tmesh.point(target(e, tmesh)) << std::endl;
#endif
          // the opposite edge might also have been inserted in the set and might still be a needle
          h = opposite(h, tmesh);
          if (is_border(h, tmesh) ) continue;
          nc = internal::is_badly_shaped(face(h, tmesh), tmesh, np);
          if (nc[0] != h)
            continue;
        }

        for (int i=0; i<2; ++i)
        {
          if (!is_border(h, tmesh))
          {
            edge_descriptor pe=edge(prev(h,tmesh), tmesh);
            edges_to_flip.erase(pe);
            next_edges_to_collapse.erase(pe);
            edges_to_collapse.erase(pe);
          }
          h = opposite(h, tmesh);
        }

        if ( !internal::is_collapse_geometrically_valid(h, tmesh, np) )
        {
          h = opposite(h, tmesh);
          if ( !internal::is_collapse_geometrically_valid(h, tmesh, np) )
          {
#ifdef  CGAL_PMP_DEBUG_REMOVE_DEGENERACIES
            std::cerr << "Warning: geometrically invalid edge collapse! "
                      << tmesh.point(source(e, tmesh)) << " "
                      << tmesh.point(target(e, tmesh)) << std::endl;
#endif
            next_edges_to_collapse.insert(e);
            continue;
          }
          e = edge(h, tmesh); // update the edge
        }

        edges_to_flip.erase(e);
        next_edges_to_collapse.erase(e); // for edges added in faces incident to a vertex kept after a collapse
#ifdef  CGAL_PMP_DEBUG_REMOVE_DEGENERACIES_EXTRA
        std::cerr << "  " << kk << " -- Collapsing " << tmesh.point(source(h, tmesh)) << "  "
                                     << tmesh.point(target(h, tmesh)) << std::endl;
#endif
        // moving to the midpoint is not a good idea. On a circle for example you might endpoint with
        // a bad geometry because you iteratively move one point
        // auto mp = midpoint(tmesh.point(source(h, tmesh)), tmesh.point(target(h, tmesh)));
        vertex_descriptor v = Euler::collapse_edge(e, tmesh);

        //tmesh.point(v) = mp;
        // examine all faces incident to the vertex kept
        for (halfedge_descriptor hv : halfedges_around_target(v, tmesh))
          if (!is_border(hv, tmesh))
            internal::collect_badly_shaped_triangles(face(hv, tmesh), tmesh, edges_to_collapse, edges_to_flip, np);

#ifdef  CGAL_PMP_DEBUG_REMOVE_DEGENERACIES_EXTRA
        std::string nb = std::to_string(++kk);
        if (kk<10) nb = std::string("0")+nb;
        if (kk<100) nb = std::string("0")+nb;
        if (kk<1000) nb = std::string("0")+nb;
        if (kk<10000) nb = std::string("0")+nb;
        std::ofstream(std::string("tmp/n-")+nb+std::string(".off")) << tmesh;
#endif
        something_was_done = true;
      }
      else
      {
#ifdef  CGAL_PMP_DEBUG_REMOVE_DEGENERACIES
        std::cerr << "Warning: uncollapsable edge! " << tmesh.point(source(e, tmesh)) << " "
                                                     << tmesh.point(target(e, tmesh)) << std::endl;
#endif
        next_edges_to_collapse.insert(e);
      }
    }

    // treat caps
#ifdef  CGAL_PMP_DEBUG_REMOVE_DEGENERACIES_EXTRA
    kk=0;
    std::ofstream(std::string("tmp/c-000.off")) << tmesh;
#endif
    while (!edges_to_flip.empty())
    {
      edge_descriptor e = *edges_to_flip.begin();
      edges_to_flip.erase(edges_to_flip.begin());

#ifdef  CGAL_PMP_DEBUG_REMOVE_DEGENERACIES
      std::cout << "treat cap: " << e << " (" << tmesh.point(source (e, tmesh)) << " --- " << tmesh.point(target(e, tmesh)) << ")" << std::endl;
#endif

      halfedge_descriptor h = halfedge(e, tmesh);
      std::array<halfedge_descriptor,2> nc = internal::is_badly_shaped(face(h, tmesh), tmesh, np);
      // First check the triangle is still a cap
      if (nc[1]!=h)
      {
#ifdef  CGAL_PMP_DEBUG_REMOVE_DEGENERACIES
        std::cerr << "Warning: Cap criteria no longer verified " << tmesh.point(source(e, tmesh)) << " --- "
                                                                  << tmesh.point(target(e, tmesh)) << std::endl;
#endif
        // the opposite edge might also have been inserted in the set and might still be a cap
        h = opposite(h, tmesh);
        if (is_border(h, tmesh) ) continue;
        nc = internal::is_badly_shaped(face(h, tmesh), tmesh, np);
        if (nc[1] != h)
          continue;
      }

      // special case on the border
      if ( is_border(opposite(h, tmesh), tmesh) )
      {
        // remove the triangle
        edges_to_flip.erase(edge(prev(h, tmesh), tmesh));
        edges_to_flip.erase(edge(next(h, tmesh), tmesh));
        next_edges_to_collapse.erase(edge(prev(h, tmesh), tmesh));
        next_edges_to_collapse.erase(edge(next(h, tmesh), tmesh));
        Euler::remove_face(h, tmesh);
        something_was_done = true;
        continue;
      }

      // condition for the flip to be valid (the edge to be created does not already exist)
      if(!halfedge(target(next(h, tmesh), tmesh),
                   target(next(opposite(h, tmesh), tmesh), tmesh), tmesh).second)
      {
#ifdef  CGAL_PMP_DEBUG_REMOVE_DEGENERACIES
        std::cout << "Flipping" << std::endl;
#endif
#ifdef  CGAL_PMP_DEBUG_REMOVE_DEGENERACIES_EXTRA
        std::cerr << "step " << kk << "\n";
        std::cerr << "  Flipping " << tmesh.point(source(h, tmesh)) << "  "
                                                   << tmesh.point(target(h, tmesh)) << std::endl;
#endif
        Euler::flip_edge(h, tmesh);
        CGAL_assertion( edge(h, tmesh) == e );

        // handle face updates
        bool undo = false;
        for (int i=0; i<2; ++i)
        {
          CGAL_assertion(!is_border(h, tmesh));
          std::array<halfedge_descriptor,2> nc = internal::is_badly_shaped(face(h, tmesh), tmesh, np);
          if (nc[1]!=boost::graph_traits<TriangleMesh>::null_halfedge())
          {
            if (edge(nc[1], tmesh) != e)
              next_edges_to_flip.insert( edge(nc[1], tmesh) );
            else
              undo = true;
          }
          else
          {
            if (nc[0]!=boost::graph_traits<TriangleMesh>::null_halfedge())
            {
              next_edges_to_collapse.insert(edge(nc[0], tmesh));
            }
          }
          h=opposite(h, tmesh);
        }

        if (undo)
        {
          /// @todo we could decide on a metric and allow the filp only if it is improving the mesh
          Euler::flip_edge(h, tmesh); // avoid infinite loop: flip only if a cap is removed
          // even if we inserted some edges above these will be tested before flipping (not optimal but it does not harm)
          next_edges_to_flip.insert(e); // we add it for the next round in case some other flip  fix the loop
          continue;
        }

        something_was_done = true;
      }
#ifdef  CGAL_PMP_DEBUG_REMOVE_DEGENERACIES
      else
      {
        std::cerr << "Warning: unflippable edge! " << tmesh.point(source(h, tmesh)) << " --- "
                                                   << tmesh.point(target(h, tmesh)) << std::endl;
        next_edges_to_flip.insert(e);
      }
#endif
#ifdef  CGAL_PMP_DEBUG_REMOVE_DEGENERACIES_EXTRA
      std::string nb = std::to_string(++kk);
      if (kk<10) nb = std::string("0")+nb;
      if (kk<100) nb = std::string("0")+nb;
      if (kk<1000) nb = std::string("0")+nb;
      if (kk<10000) nb = std::string("0")+nb;
      std::ofstream(std::string("tmp/c-")+nb+std::string(".off")) << tmesh;
#endif
    }

    std::swap(edges_to_collapse, next_edges_to_collapse);
    std::swap(edges_to_flip, next_edges_to_flip);

    if (!something_was_done)
      return false;
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
