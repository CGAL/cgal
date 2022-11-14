// Copyright (c) 2015 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Ilker O. Yaz

#ifndef CGAL_HOLE_FILLING_TRIANGULATE_HOLE_POLYHEDRON_3_H
#define CGAL_HOLE_FILLING_TRIANGULATE_HOLE_POLYHEDRON_3_H

#include <CGAL/license/Polygon_mesh_processing/meshing_hole_filling.h>


#include <CGAL/Polygon_mesh_processing/internal/Hole_filling/Triangulate_hole_polyline.h>
#ifdef CGAL_PMP_HOLE_FILLING_DEBUG
#include <CGAL/Timer.h>
#endif
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/boost/graph/Euler_operations.h>
#include <vector>

namespace CGAL {
namespace Polygon_mesh_processing {
namespace internal {


template<class PolygonMesh, class OutputIterator>
struct Tracer_polyhedron
{
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;

  Tracer_polyhedron(OutputIterator out,
                    PolygonMesh& pmesh,
                    std::vector<halfedge_descriptor>& P)
    : out(out), pmesh(pmesh), P(P)
  { }

  template <class LookupTable>
  halfedge_descriptor
  operator()(const LookupTable& lambda,
             int i, int k,
             bool last = true)
  {
    if(i + 1 == k)
      return P[i+1];

    halfedge_descriptor h, g;
    if(i+2 == k)
    {
      if(last)
      {
        h = P[i + 1];
        Euler::fill_hole(h, pmesh);
      }
      else
      {
        h = Euler::add_face_to_border(prev(P[i + 1], pmesh), P[i + 2 /*k*/], pmesh);
      }

      CGAL_assertion(face(h,pmesh) != boost::graph_traits<PolygonMesh>::null_face());
      *out++ = face(h,pmesh);
      return opposite(h,pmesh);
    }
    else
    {
      int la = lambda.get(i, k);
      h = operator()(lambda, i, la, false);
      g = operator()(lambda, la, k, false);

      if(last)
      {
        h = g;
        Euler::fill_hole(g, pmesh);
      }
      else
      {
        h = Euler::add_face_to_border(prev(h, pmesh), g, pmesh);
      }

      CGAL_assertion(face(h,pmesh) != boost::graph_traits<PolygonMesh>::null_face());
      *out++ = face(h,pmesh);
      return opposite(h,pmesh);
    }
  }

  OutputIterator out;
  PolygonMesh& pmesh;
  std::vector<halfedge_descriptor>& P;
};

// This function is used in test cases (since it returns not just OutputIterator but also Weight)
template<class PolygonMesh, class OutputIterator, class VertexPointMap, class Kernel, class Visitor>
std::pair<OutputIterator, CGAL::internal::Weight_min_max_dihedral_and_area>
triangulate_hole_polygon_mesh(PolygonMesh& pmesh,
            typename boost::graph_traits<PolygonMesh>::halfedge_descriptor border_halfedge,
            OutputIterator out,
            VertexPointMap vpmap,
            bool use_delaunay_triangulation,
            const Kernel& k,
            const bool use_cdt,
            const bool skip_cubic_algorithm,
            Visitor& visitor,
            const typename Kernel::FT max_squared_distance)
{
  typedef Halfedge_around_face_circulator<PolygonMesh>   Hedge_around_face_circulator;
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename Kernel::Point_3 Point_3;

  typedef std::map<vertex_descriptor, int>    Vertex_map;
  typedef typename Vertex_map::iterator       Vertex_map_it;

#ifdef CGAL_PMP_HOLE_FILLING_DEBUG
  CGAL::Timer timer; timer.start();
#endif

  std::vector<Point_3> P, Q;
  std::vector<halfedge_descriptor> P_edges;
  Vertex_map vertex_map;

  int id = 0;
  Hedge_around_face_circulator circ(border_halfedge,pmesh), done(circ);
  do
  {
    P.push_back(get(vpmap, target(*circ, pmesh)));
    Q.push_back(get(vpmap, target(next(opposite(next(*circ,pmesh),pmesh),pmesh),pmesh)));
    P_edges.push_back(*circ);
    if(!vertex_map.insert(std::make_pair(target(*circ,pmesh), id++)).second)
    {
#ifndef CGAL_TEST_SUITE
      CGAL_warning_msg(false, "Returning no output. Non-manifold vertex is found on boundary!");
#else
      std::cerr << "W: Returning no output. Non-manifold vertex is found on boundary!\n";
#endif
      return std::make_pair(out, CGAL::internal::Weight_min_max_dihedral_and_area::NOT_VALID());
    }
  } while (++circ != done);

  // existing_edges contains neighborhood information between boundary vertices
  // more precisely if v_i is neighbor to any other vertex than v_(i-1) and v_(i+1),
  // this edge is put into existing_edges
  std::vector<std::pair<int, int> > existing_edges;
  for(Vertex_map_it v_it = vertex_map.begin(); v_it != vertex_map.end(); ++v_it)
  {
    int v_it_id = v_it->second;
    int v_it_prev = (v_it_id == 0)   ? (id-1) : (v_it_id-1);
    int v_it_next = (v_it_id == id-1) ?  0    : (v_it_id+1);

    Halfedge_around_target_circulator<PolygonMesh>
      circ_vertex(halfedge(v_it->first,pmesh),pmesh),
      done_vertex(circ_vertex);
    do
    {
      Vertex_map_it v_it_neigh_it = vertex_map.find(source(*circ_vertex,pmesh));

      if(v_it_neigh_it != vertex_map.end()) //other endpoint found in the map
      {
        int v_it_neigh_id = v_it_neigh_it->second;
        if(v_it_neigh_id != v_it_prev && v_it_neigh_id != v_it_next)
        {
          //there is an edge incident to v_it, which is not next or previous
          //from vertex_map (checked by comparing IDs)
          if(v_it_id < v_it_neigh_id) // to include each edge only once
            existing_edges.push_back(std::make_pair(v_it_id, v_it_neigh_id));
        }
      }
    } while(++circ_vertex != done_vertex);
  }

//#define CGAL_USE_WEIGHT_INCOMPLETE
#ifdef CGAL_USE_WEIGHT_INCOMPLETE
  typedef CGAL::internal::Weight_calculator<Weight_incomplete<CGAL::internal::Weight_min_max_dihedral_and_area>,
        CGAL::internal::Is_valid_existing_edges_and_degenerate_triangle> WC;
#else
  typedef CGAL::internal::Weight_calculator<CGAL::internal::Weight_min_max_dihedral_and_area,
        CGAL::internal::Is_valid_existing_edges_and_degenerate_triangle> WC;
#endif

  CGAL::internal::Is_valid_existing_edges_and_degenerate_triangle is_valid(existing_edges);

  // fill hole using polyline function, with custom tracer for PolygonMesh
  Tracer_polyhedron<PolygonMesh, OutputIterator> tracer(out, pmesh, P_edges);

#ifndef CGAL_HOLE_FILLING_DO_NOT_USE_CDT2
  if(use_cdt && triangulate_hole_polyline_with_cdt(P, tracer, visitor, is_valid, k, max_squared_distance))
    return std::make_pair(tracer.out, CGAL::internal::Weight_min_max_dihedral_and_area(0,0));
#endif
  CGAL::internal::Weight_min_max_dihedral_and_area weight =
#ifndef CGAL_USE_WEIGHT_INCOMPLETE
  triangulate_hole_polyline(P, Q, tracer, WC(is_valid), visitor, use_delaunay_triangulation, skip_cubic_algorithm, k);
#else
  // get actual weight in Weight_incomplete
  triangulate_hole_polyline(P, Q, tracer, WC(is_valid), visitor, use_delaunay_triangulation, k).weight;
#endif

#ifdef CGAL_PMP_HOLE_FILLING_DEBUG
  std::cerr << "Hole filling: " << timer.time() << " sc." << std::endl; timer.reset();
#endif

  return std::make_pair(tracer.out, weight);
}

} // namespace internal
} // namespace Polygon_mesh_processing
} // namespace CGAL

#endif //CGAL_HOLE_FILLING_TRIANGULATE_HOLE_POLYHEDRON_3_H
