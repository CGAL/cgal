// Copyright (c) 2023-2026 GeometryFactory and Claudio Mancinelli.
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Claudio Mancinelli and Sébastien Loriot

#ifndef CGAL_POLYGON_MESH_PROCESSING_BSURF_RECURSIVE_DE_CASTELJAU_H
#define CGAL_POLYGON_MESH_PROCESSING_BSURF_RECURSIVE_DE_CASTELJAU_H

#include <CGAL/license/Vector_graphics_on_surfaces.h>

#include <CGAL/Vector_graphics_on_surfaces/locally_shortest_path.h>


namespace CGAL {
namespace Vector_graphics_on_surfaces {

namespace internal {

template <class K, class TriangleMesh, class VertexPointMap>
struct Bezier_tracing_impl
{
  using face_descriptor =
      typename boost::graph_traits<TriangleMesh>::face_descriptor;
  using vertex_descriptor =
      typename boost::graph_traits<TriangleMesh>::vertex_descriptor;
  using halfedge_descriptor =
      typename boost::graph_traits<TriangleMesh>::halfedge_descriptor;

  using Point_2 = typename K::Point_2;
  using Point_3 = typename K::Point_3;
  using Vector_2 = typename K::Vector_2;
  using Vector_3 = typename K::Vector_3;
  using FT = typename K::FT;

  using Face_location = Polygon_mesh_processing::Face_location<TriangleMesh, FT>;
  using Edge_location = Polygon_mesh_processing::Edge_location<TriangleMesh, FT>;

  template <class EdgeLocationRange>
  static
  std::vector<Point_3>
  get_positions(const EdgeLocationRange& edge_locations,
                const TriangleMesh& mesh,
                const Face_location& src,
                const Face_location& tgt)
  {
    namespace PMP = CGAL::Polygon_mesh_processing;

    std::vector<Point_3> result;
    result.reserve(edge_locations.size()+2);
    result.push_back(PMP::construct_point(src,mesh));
    for(auto& e: edge_locations)
        result.push_back(PMP::construct_point(e,mesh));

    result.push_back(PMP::construct_point(tgt,mesh));
//TODO: we must guarantee that result is sorted and unique (rounding issue?)
    return result;
  }

  template <class EdgeLocationRange>
  static
  std::vector<FT>
  path_parameters(const EdgeLocationRange& edge_locations,
                  const TriangleMesh& mesh,
                  const Face_location& src,
                  const Face_location& tgt)
  {
    std::vector<Point_3> pos=get_positions(edge_locations,mesh,src,tgt);
    FT L=0.;
    std::vector<FT> result(pos.size());
    for(std::size_t i=0;i<pos.size();++i)
    {
      if(i) L+=sqrt(squared_distance(pos[i],pos[i-1]));
      result[i]=L;
    }

    for(auto& t:result) t/=L;

    return result;
  }

  template <class EdgeLocationRange>
  static
  Face_location
  eval_point_on_geodesic(const EdgeLocationRange& edge_locations,
                         const TriangleMesh& mesh,
                         const Face_location& src,
                         const Face_location& tgt,
                         const std::vector<FT>& parameters,/// edge length parameterization of the path from src to tgt through edge_locations
                         const FT& t)
  {
    if (t==0) return src;
    if (t==1) return tgt;

    if(src.first==tgt.first)
    {
      std::array<FT,3> bary;
      bary[0]=(1-t)*src.second[0]+t*tgt.second[0];
      bary[1]=(1-t)*src.second[1]+t*tgt.second[1];
      bary[2]=(1-t)*src.second[2]+t*tgt.second[2];
      return {src.first,bary};
    }

    std::size_t i = 0;
    for (; i < parameters.size() - 1; i++)
    {
      if (parameters[i + 1] >= t) break;
    }
    FT t_low = parameters[i];
    FT t_high = parameters[i + 1];
    CGAL_assertion(t_high!=t_low);
    FT alpha = (t - t_low) / (t_high - t_low);
    std::array<FT,3> bary_low;
    std::array<FT,3> bary_high;

    // warning there is an offset of the index: parameters contains one extra element (src) at 0
    // while edge_locations does not
    face_descriptor curr_tid = i==0?src.first:face(halfedge(edge_locations[i-1].first,mesh),mesh);
    halfedge_descriptor h_face = halfedge(curr_tid, mesh);
    auto edge_barycentric_coordinate =
      [&mesh, h_face](halfedge_descriptor h_edge,
                     const std::array<FT,2>& bary_edge)
    {
      std::array<FT,3> bary_edge_in_face;
      if (h_face!=h_edge)
      {
        if (h_face==next(h_edge, mesh))
        {
          bary_edge_in_face[0]=bary_edge[1];
          bary_edge_in_face[1]=0;
          bary_edge_in_face[2]=bary_edge[0];
        }
        else
        {
          bary_edge_in_face[0]=0;
          bary_edge_in_face[1]=bary_edge[0];
          bary_edge_in_face[2]=bary_edge[1];
        }
      }
      else
      {
        bary_edge_in_face[0]=bary_edge[0];
        bary_edge_in_face[1]=bary_edge[1];
        bary_edge_in_face[2]=0;
      }

      return bary_edge_in_face;
    };

    if(i==0)
      bary_low=src.second;
    else
    {
      halfedge_descriptor h_low = halfedge(edge_locations[i-1].first, mesh);
      bary_low = edge_barycentric_coordinate(h_low, edge_locations[i-1].second);
    }

    if(i==parameters.size()-2)
      bary_high=tgt.second;
    else
    {
      halfedge_descriptor h_high = opposite(halfedge(edge_locations[i].first, mesh), mesh);
      CGAL_assertion(face(h_high,mesh)==curr_tid);
      std::array<FT,2> edge_bary_high=edge_locations[i].second;
      std::swap(edge_bary_high[0],edge_bary_high[1]);
      bary_high = edge_barycentric_coordinate(h_high, edge_bary_high);
    }

    std::array<FT,3> bary;
    bary[0]=(1-alpha)*bary_low[0]+alpha*bary_high[0];
    bary[1]=(1-alpha)*bary_low[1]+alpha*bary_high[1];
    bary[2]=(1-alpha)*bary_low[2]+alpha*bary_high[2];

    return {curr_tid,bary};
  }

  static
  Face_location
  geodesic_lerp(const TriangleMesh &mesh,
                const Face_location& src,
                const Face_location& tgt,const FT& t
#ifndef CGAL_BSURF_USE_DIJKSTRA_SP
               , const Dual_geodesic_solver<FT>& solver
#endif
  )
  {
    std::vector<Edge_location> edge_locations;
    locally_shortest_path<FT>(src,tgt,mesh, edge_locations, solver);
    std::vector<FT> parameters=path_parameters(edge_locations,mesh,src,tgt);
    Face_location point = eval_point_on_geodesic(edge_locations,mesh,src,tgt,parameters,t);
    return point;
  }


  static
  std::pair<Bezier_segment<TriangleMesh, FT>,Bezier_segment<TriangleMesh,FT>>
  subdivide_Bezier_polygon(const TriangleMesh& mesh,
                           const Bezier_segment<TriangleMesh,FT>& polygon,
                           const FT& t
#ifndef CGAL_BSURF_USE_DIJKSTRA_SP
                           , const Dual_geodesic_solver<FT>& solver
#endif
  )
  {
#ifndef CGAL_BSURF_USE_DIJKSTRA_SP
     Face_location Q0 = geodesic_lerp(mesh, polygon[0], polygon[1], t, solver);
     Face_location Q1 = geodesic_lerp(mesh, polygon[1], polygon[2], t, solver);
     Face_location Q2 = geodesic_lerp(mesh, polygon[2], polygon[3], t, solver);
     Face_location R0 = geodesic_lerp(mesh, Q0, Q1, t, solver);
     Face_location R1 = geodesic_lerp(mesh, Q1, Q2, t, solver);
     Face_location S  = geodesic_lerp(mesh, R0, R1, t, solver);
#else
     Face_location Q0 = geodesic_lerp(mesh, polygon[0], polygon[1], t);
     Face_location Q1 = geodesic_lerp(mesh, polygon[1], polygon[2], t);
     Face_location Q2 = geodesic_lerp(mesh, polygon[2], polygon[3], t);
     Face_location R0 = geodesic_lerp(mesh, Q0, Q1, t);
     Face_location R1 = geodesic_lerp(mesh, Q1, Q2, t);
     Face_location S  = geodesic_lerp(mesh, R0, R1, t);
#endif
    return {{polygon[0], Q0, R0, S}, {S, R1, Q2, polygon[3]}};
  }
};

} // end of internal namespace

/*!
 * \ingroup VGSFunctions
 * computes a discretization of a Bézier segment defined by the location of four control points on `tmesh`.
 * All control points must be on the same connected component. This functions applies several iterations of
 * the de Casteljau algorithm, and geodesic shortest paths are drawn between the control points.
 * \tparam TriangleMesh a model of `FaceListGraph` and `EdgeListGraph`
 * \tparam FT floating point number type (float or double)
 * \param mesh input triangle mesh to compute the path on
 * \param control_points control points of the Bézier segment
 * \param num_subdiv the number of iterations of the subdivision algorithm
 * \param solver container for the precomputed information. If not initialized, it will be initialized internally.
 * \return discretization of the Bézier segment as face locations
 * \todo add named parameters
 * \todo do we want to also have a way to return Bézier segments? The output is actually Bézier segments subdivided.
 */
template <class TriangleMesh, class FT>
std::vector<CGAL::Polygon_mesh_processing::Face_location<TriangleMesh, FT>>
recursive_de_Casteljau(const TriangleMesh& mesh,
                       const Bezier_segment<TriangleMesh, FT>& control_points,
                       const int num_subdiv
#ifndef CGAL_BSURF_USE_DIJKSTRA_SP
                       , const Dual_geodesic_solver<FT>& solver = Dual_geodesic_solver<FT>()
#endif
                       )
{
  namespace PMP = CGAL::Polygon_mesh_processing;

  //TODO replace with named parameter
  using VPM = typename boost::property_map<TriangleMesh, CGAL::vertex_point_t>::const_type;
  using K =  typename Kernel_traits<typename boost::property_traits<VPM>::value_type>::type;
  using Impl = internal::Bezier_tracing_impl<K, TriangleMesh, VPM>;

  // init solver if empty
  const Dual_geodesic_solver<FT>* solver_ptr=&solver;
  Dual_geodesic_solver<FT> local_solver;
  if (solver.graph.empty())
  {
    solver_ptr = &local_solver;
    init_geodesic_dual_solver(local_solver, mesh);
  }

  std::vector<Bezier_segment<TriangleMesh, FT>> segments(1,control_points);
  std::vector<Bezier_segment<TriangleMesh, FT>> result;
  for (auto subdivision = 0; subdivision < num_subdiv; ++subdivision)
  {
    result.clear();
    result.reserve(segments.size() * 2);
    for (std::size_t i = 0; i < segments.size(); ++i)
    {
      auto [split0, split1] = Impl::subdivide_Bezier_polygon(mesh, segments[i], 0.5, *solver_ptr);
      result.push_back(split0);
      result.push_back(split1);
    }

    CGAL_assertion( 2*segments.size()==result.size() );

    std::swap(segments, result);
  }

  std::vector<PMP::Face_location<TriangleMesh, FT>> final_result;
  final_result.reserve(result.size() * 3 + 1);
  final_result.push_back(std::move(result.front()[0]));
  for (Bezier_segment<TriangleMesh, FT>& bs : result)
  {
    for (int i=1;i<4;++i)
      final_result.push_back(std::move(bs[i]));
  }

#if 0
  // nasty trick to build the vector from a pair of iterators
  // using the fact that data in array and vector are contiguous
  return {(PMP::Face_location<TriangleMesh, FT>*)segments.data(),
          (PMP::Face_location<TriangleMesh, FT>*)segments.data() + segments.size() * 4};
#endif
  return final_result;
}

} } // end of CGAL::Vector_graphics_on_surfaces namespace

#endif // CGAL_POLYGON_MESH_PROCESSING_BSURF_RECURSIVE_DE_CASTELJAU_H
