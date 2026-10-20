// Copyright (c) 2026 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sébastien Loriot

#ifndef CGAL_POLYGON_MESH_PROCESSING_TRIANGLE_SOUP_BOOLEAN_OPERATIONS_3_H
#define CGAL_POLYGON_MESH_PROCESSING_TRIANGLE_SOUP_BOOLEAN_OPERATIONS_3_H

#include <CGAL/license/Polygon_mesh_processing/autorefinement.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/autorefinement.h>
#include <CGAL/Polygon_mesh_processing/internal/Corefinement/predicates.h>
#include <CGAL/Union_find.h>

#include <fstream>

#ifdef CGAL_LINKED_WITH_TBB
#include <oneapi/tbb/parallel_sort.h>
#endif

#if defined(CGAL_BO3_TIMERS) || defined(CGAL_BO3_SELECT_TIMERS)
#include <CGAL/Real_timer.h>
#endif

namespace CGAL {
namespace boolops_3 {

//TODO: update nm_edges for identical coplanar status (less test below then)
template <class EdgeConnectivity, class FaceConnectivity>
std::vector<std::size_t>
connected_components(std::size_t begin,
                     std::size_t end,
                     int offset,
                     int increment,
                     const EdgeConnectivity& adjacent_edges,
                     const FaceConnectivity& adjacent_faces,
                     std::vector<std::size_t>& ccids)
{
  auto is_handled = [&](std::size_t tid){ return ccids[tid]!=std::size_t(-1); };
  std::vector<std::size_t> one_triangle_per_cc;

  std::size_t i=offset;

  std::vector<std::size_t> queue;
  for (std::size_t tid=begin; tid<end; ++tid)
  {
    if (is_handled(tid)) continue;
    queue.push_back(tid);
    ccids[tid]=i;
    while(!queue.empty())
    {
      std::size_t q_tid = queue.back();
      queue.pop_back();

      for (int k=0; k<3; ++k) // TODO: change that if you want more than triangles
      {
        auto ei = adjacent_edges[q_tid][k];
        const auto &tids = adjacent_faces[ei];
        CGAL_assertion(tids.size() > 0);
        switch( tids.size() )
        {
          case 1:
          // boundary edge
          break;
          case 2:
          {
            std::size_t n_tid = tids[0]!=q_tid ? tids[0] : tids[1];
            if (!is_handled(n_tid)){
              ccids[n_tid]=i;
              queue.push_back(n_tid);
            }
          }
          break;
          default:
          // non-manifold edge
          break;
        }
      }
    }
    one_triangle_per_cc.push_back(tid);
    i += increment;
  }
  return one_triangle_per_cc;
}

template <class EdgeConnectivity, class FaceConnectivity>
std::vector<std::size_t>
connected_components(const EdgeConnectivity& adjacent_edges,
                     const FaceConnectivity& adjacent_faces,
                     std::vector<std::size_t>& ccids)
{
  std::size_t nb_poly = adjacent_edges.size();
  ccids.resize(nb_poly, std::size_t(-1));
  return connected_components(0, nb_poly, 0, 1, adjacent_edges, adjacent_faces, ccids);
}

#ifdef CGAL_LINKED_WITH_TBB
// we actually only parallelise the exploraton in the triangles of the two (for now)
// input soups. Since autorefine only creates new triangles and subdivided one are bounded
// by non-manifold edges we have 2 non connected ranges of triangles that can be processed
// separately.
template <class EdgeConnectivity, class FaceConnectivity>
std::vector<std::size_t>
parallel_connected_components(const size_t polyids_end1,
                              const size_t polyids_end2,
                              const EdgeConnectivity& adjacent_edges,
                              const FaceConnectivity& adjacent_faces,
                              std::vector<std::size_t>& ccids)
{
  std::size_t nb_poly = adjacent_edges.size();
  ccids.resize(nb_poly, std::size_t(-1));
  std::vector<std::size_t> trs_per_cc_1, trs_per_cc_2, trs_per_cc_end, merged_trs_per_cc;

  tbb::task_group task_group;
  task_group.run([&](){
    trs_per_cc_1 = connected_components(0, polyids_end1, 0, 2, adjacent_edges, adjacent_faces, ccids);
  });
  trs_per_cc_2 = connected_components(polyids_end1, polyids_end2, 1, 2, adjacent_edges, adjacent_faces, ccids);
  task_group.wait();
  std::size_t max_size = (std::max)(trs_per_cc_1.size(), trs_per_cc_2.size());
  trs_per_cc_end = connected_components(polyids_end2, nb_poly, max_size, 1, adjacent_edges, adjacent_faces, ccids);

  merged_trs_per_cc.resize(2 * max_size + trs_per_cc_end.size(), std::size_t(-1));
  for(std::size_t i=0; i<trs_per_cc_1.size(); ++i)
    merged_trs_per_cc[2*i] = trs_per_cc_1[i];
  for(std::size_t i=0; i<trs_per_cc_2.size(); ++i)
    merged_trs_per_cc[2*i+1] = trs_per_cc_2[i];
  for(std::size_t i=0; i<trs_per_cc_end.size(); ++i)
    merged_trs_per_cc[max_size+i] = trs_per_cc_end[i];
  return merged_trs_per_cc;
}
#endif

struct Selection_data
{
  std::vector<std::size_t> ccids;
  std::vector<bool> coplanar_patches;
  std::vector<bool> coplanar_patches_for_union_and_intersection;
  std::vector< std::vector<bool> > cc_is_inside_per_mesh;
  std::vector<unsigned int> input_ids;
};

template <class Concurrency_tag, class Exact_point, class Triangle>
Selection_data
select(const std::vector<Exact_point>& points,
       const std::vector<Triangle>& triangles,
             std::vector<unsigned int>& input_ids,
#ifdef CGAL_LINKED_WITH_TBB
             std::size_t polyids_end1,
             std::size_t polyids_end2,
#else
             std::size_t /* polyids_end1 */,
             std::size_t /* polyids_end2 */,
#endif
       const unsigned int nb_meshes)
{
#ifdef CGAL_BO3_SELECT_TIMERS
  Real_timer t;
  t.start();
#endif

  namespace PMP = CGAL::Polygon_mesh_processing;
  namespace pred = PMP::Corefinement;
  Selection_data data_out;

  // Create connectivity information
  using PID = std::size_t;
  using EID = std::size_t;
  using TID = std::size_t;

  // create a temporary vector to store all the edges of the triangles.
  // We will then sort those triangle edges so that identical edges are consecutive in the container.
  // We will then use this adjacency to build a data structure providing the adjacency between triangles.
  using Edge_record = std::tuple<PID, PID, TID>;
  std::vector<Edge_record> all_triangle_edges;
  all_triangle_edges.resize(3 * triangles.size());
#ifdef CGAL_LINKED_WITH_TBB
  if constexpr (Concurrency_tag::is_parallel)
  {
    tbb::parallel_for(std::size_t(0), triangles.size(), [&](std::size_t tid)
    {
      auto add_edge = [&](std::size_t idx, PID a, PID b){
        if (b < a) std::swap(a, b);
        all_triangle_edges[idx] = Edge_record{a, b, tid};
      };
      add_edge(3 * tid + 0, triangles[tid][0], triangles[tid][1]);
      add_edge(3 * tid + 1, triangles[tid][1], triangles[tid][2]);
      add_edge(3 * tid + 2, triangles[tid][2], triangles[tid][0]);
    });
  }
  else
#endif
  {
    for (TID tid = 0; tid < triangles.size(); ++tid)
    {
      auto add_edge = [&](std::size_t idx, PID a, PID b){
        if (b < a) std::swap(a, b);
        all_triangle_edges[idx] = Edge_record{a, b, tid};
      };
      add_edge(3 * tid + 0, triangles[tid][0], triangles[tid][1]);
      add_edge(3 * tid + 1, triangles[tid][1], triangles[tid][2]);
      add_edge(3 * tid + 2, triangles[tid][2], triangles[tid][0]);
    }
  }

#ifdef CGAL_BO3_SELECT_TIMERS
  std::cout << t.time() << " sec. to collect edge records\n";
#endif

#ifdef CGAL_LINKED_WITH_TBB
  if constexpr (Concurrency_tag::is_parallel)
  {
    tbb::parallel_sort(all_triangle_edges.begin(), all_triangle_edges.end());
  }
  else
#endif
  {
    std::sort(all_triangle_edges.begin(), all_triangle_edges.end());
  }

#ifdef CGAL_BO3_SELECT_TIMERS
  std::cout << t.time() << " sec. to sort edge records\n";
#endif

  // Build adjacencies between edges and faces
  std::vector<std::array<EID, 3>> adjacent_edges(triangles.size());
  std::vector<boost::container::small_vector<TID, 2>> adjacent_faces;
#ifdef CGAL_LINKED_WITH_TBB
  if constexpr (Concurrency_tag::is_parallel)
  {
    // split all_triangle_edges in blocks and compute the number of times we change edges for each block.
    // for the first block we have an extra +1 making first_edge_id_per_block[nb_blocks] the number of edges
    const std::size_t block_size = 4096;
    const std::size_t nb_blocks = (all_triangle_edges.size() + block_size - 1) / block_size;

    std::vector<std::size_t> nb_edge_changes_per_block(nb_blocks, 0);

    // count how many times we change from an edge to the next one into each block
    tbb::parallel_for(tbb::blocked_range<std::size_t>(0, nb_blocks, 1), [&](const tbb::blocked_range<std::size_t>& r)
    {
      for (std::size_t bi = r.begin(); bi != r.end(); ++bi){
        const std::size_t begin = bi * block_size;
        const std::size_t end = (std::min)(begin + block_size, all_triangle_edges.size());

        for (std::size_t i = begin; i < end; ++i)
          if (i == 0 || std::get<0>(all_triangle_edges[i]) != std::get<0>(all_triangle_edges[i-1]) ||
                        std::get<1>(all_triangle_edges[i]) != std::get<1>(all_triangle_edges[i-1]))
            ++nb_edge_changes_per_block[bi];
      }
    });

    // Compute offset for each block
    std::vector<std::size_t> first_edge_id_per_block(nb_blocks+1, 0);
    for (std::size_t i=0; i < nb_blocks; ++i)
      first_edge_id_per_block[i+1] = first_edge_id_per_block[i] + nb_edge_changes_per_block[i];

    adjacent_faces.resize(first_edge_id_per_block.back());

    tbb::parallel_for(tbb::blocked_range<std::size_t>(0, nb_blocks, 1), [&](const tbb::blocked_range<std::size_t>& r)
    {
      for (std::size_t bi = r.begin(); bi != r.end(); ++bi)
      {
        std::size_t begin = bi * block_size;
        const std::size_t end = (std::min)(begin + block_size, all_triangle_edges.size());
        std::size_t idx = first_edge_id_per_block[bi];

        // Skip a group that started in the previous block.
        if(begin != 0)
          while(std::get<0>(all_triangle_edges[begin-1]) == std::get<0>(all_triangle_edges[begin]) &&
                std::get<1>(all_triangle_edges[begin-1]) == std::get<1>(all_triangle_edges[begin]))
            ++begin;

        for (std::size_t i = begin; i < end; )
        {
          const PID a = std::get<0>(all_triangle_edges[i]);
          const PID b = std::get<1>(all_triangle_edges[i]);

          auto& faces = adjacent_faces[idx];
          CGAL_assertion(faces.empty());

          // may extend beyond the block's end, but is processed only once.
          do
          {
            const TID tid = std::get<2>(all_triangle_edges[i]);
            const auto& t = triangles[tid];
            faces.push_back(tid);

            const std::size_t k = (t[0] == a || t[0] == b) ? ((t[1] == a || t[1] == b) ? 2 : 1) : 0;
            CGAL_assertion((t[(k+1)%3] == a && t[(k+2)%3] == b) || (t[(k+1)%3] == b && t[(k+2)%3] == a));
            adjacent_edges[tid][k] = idx;
            ++i;
          } while (i < all_triangle_edges.size() && std::get<0>(all_triangle_edges[i]) == a &&
                                                    std::get<1>(all_triangle_edges[i]) == b);
          ++idx;
        }
        CGAL_assertion(idx == first_edge_id_per_block[bi+1]);
      }
    });
  }
  else
#endif
  {
    adjacent_faces.reserve(all_triangle_edges.size());
    for (std::size_t i = 0; i < all_triangle_edges.size(); )
    {
      const PID a = std::get<0>(all_triangle_edges[i]);
      const PID b = std::get<1>(all_triangle_edges[i]);
      const std::size_t idx = adjacent_faces.size();
      adjacent_faces.emplace_back();
      CGAL_assertion(adjacent_faces.back().empty());
      do
      {
        const TID tid = std::get<2>(all_triangle_edges[i]);
        const auto& t = triangles[tid];
        adjacent_faces.back().push_back(tid);
        ++i;

        std::size_t k = (t[0] == a || t[0] == b) ? ((t[1] == a || t[1] == b) ? 2 : 1) : 0;
        CGAL_assertion((t[(k+1)%3] == a && t[(k+2)%3] == b) || (t[(k+1)%3] == b && t[(k+2)%3] == a));
        adjacent_edges[tid][k] = idx;
      } while (i < all_triangle_edges.size() && std::get<0>(all_triangle_edges[i]) == a && std::get<1>(all_triangle_edges[i]) == b);
      CGAL_assertion(!adjacent_faces.back().empty());
    }
  }
  all_triangle_edges.clear(); // no longer needed

#ifdef CGAL_BO3_SELECT_TIMERS
  std::cout << t.time() << " sec. to build adjacencies\n";
#endif

  // extract connected_components bounded by non-manifold edges
  std::vector<std::size_t>& ccids = data_out.ccids;
  std::vector<std::size_t> one_triangle_per_cc;
#ifdef CGAL_LINKED_WITH_TBB
  if constexpr(Concurrency_tag::is_parallel)
  {
    one_triangle_per_cc = boolops_3::parallel_connected_components(polyids_end1, polyids_end2, adjacent_edges, adjacent_faces, ccids);
  }
  else
#endif
  {
    one_triangle_per_cc = boolops_3::connected_components(adjacent_edges, adjacent_faces, ccids);
  }
  std::size_t nb_cc = one_triangle_per_cc.size();

#ifdef CGAL_BO3_SELECT_TIMERS
  std::cout << t.time() << " sec. to compute connected components\n";
#endif

  // collect ids of non-manifold edges
  std::vector< EID > nm_edges;
  for (EID eid=0; eid<adjacent_faces.size(); ++eid)
    if (adjacent_faces[eid].size()>2)
      nm_edges.push_back(eid);

#ifdef CGAL_BO_AUTOREF_DEBUG
  std::cout << "nb_cc " << nb_cc << "\n";
  std::vector< std::vector< boost::container::small_vector<std::size_t, 3> > > triangles_per_cc(nb_cc);
  for (std::size_t tid=0; tid<triangles.size(); ++tid)
    triangles_per_cc[ccids[tid]].push_back(triangles[tid]);
  for (std::size_t i=0; i<nb_cc; ++i)
  {
    std::vector<Exact_point> cc_points = points;
    PMP::remove_isolated_points_in_polygon_soup(cc_points, triangles_per_cc[i]);
    CGAL::IO::write_OFF("cc_"+std::to_string(i)+".off", cc_points, triangles_per_cc[i], CGAL::parameters::stream_precision(17));
  }
#endif

  std::size_t nb_handled = 0;
  std::vector<bool> cc_handled(nb_cc, false);
  std::vector<bool>& coplanar_patches=data_out.coplanar_patches;
  std::vector<bool>& coplanar_patches_for_union_and_intersection=data_out.coplanar_patches_for_union_and_intersection;
  std::vector< std::vector<bool> >& cc_is_inside_per_mesh=data_out.cc_is_inside_per_mesh;

  for(std::size_t cc=0; cc<nb_cc; ++cc)
    if( one_triangle_per_cc[cc] == std::size_t(-1)){
      ++nb_handled;
      cc_handled[cc] = true;
    }

  coplanar_patches.assign(nb_cc, false);
  coplanar_patches_for_union_and_intersection.assign(nb_cc, false);
#ifdef CGAL_BO_AUTOREF_DEBUG
  // Note: if we want to avoid this debug vector, we can use an optional in cc_is_inside_per_mesh.
  // However, we would loose the compactness of std::vector<bool>, even if nb_cc is not that big.
  std::vector< std::vector<bool> > cc_is_inside_per_mesh_init(nb_meshes, std::vector<bool>(nb_cc, false));
  int case_id=-1;
#endif
  cc_is_inside_per_mesh.assign(nb_meshes, std::vector<bool>(nb_cc, false));

  for (auto e : nm_edges)
  {
    if (nb_handled==nb_cc) break;

    bool all_handled=true;
    auto &tids = adjacent_faces[e];
    for (TID tid : tids)
    {
      if (!cc_handled[ccids[tid]])
      {
        all_handled=false;
        break;
      }
    }

    if (all_handled) continue;

    bool all_from_same_mesh = true;

    for (TID tid : tids)
    {
      if (input_ids[tid]!=input_ids[tids.front()])
      {
        all_from_same_mesh=false;
        break;
      }
    }
    if (all_from_same_mesh) continue;
#ifdef CGAL_BO_AUTOREF_DEBUG
    std::cout << "case #" << ++case_id << " looking at edge " << pid << " " << e << "\n";
    std::cout << "  incident CCs:";
    for (TID tid : tids) std::cout << " " << ccids[tid] << "("<< input_ids[tid] <<")";
    std::cout << "\n";
#endif
    for (TID tid : tids)
    {
      if (!cc_handled[ccids[tid]]) ++nb_handled;
      cc_handled[ccids[tid]] = true;
    }

    std::vector<std::pair<int, TID>> to_sort;
    to_sort.reserve(tids.size());

    for (TID tid : tids)
    {
      int i = adjacent_edges[tid][0] == e ? 0 : (adjacent_edges[tid][1] == e ? 1 : 2);
      CGAL_assertion( adjacent_edges[tid][i] == e );
      to_sort.emplace_back(i, tid);
    }

    // we arbitrarily take the last triangle
    std::pair<int, TID> tref=to_sort.front();
    std::size_t i0=triangles[tref.second][(tref.first+1)%3],
                i1=triangles[tref.second][(tref.first+2)%3];

    const Exact_point& pref=points[triangles[tref.second][tref.first]];
    const Exact_point& p0=points[i0];
    const Exact_point& p1=points[i1];

    //Considering the plane with normal vector [o_prime,o] and containing o.
    //We define the counterclockwise order around o when looking from
    //the side of the plane containing o_prime.
    //We consider the portion of the plane defined by rotating a ray starting at o
    //from the planar projection of p1 to the planar projection of p2 in
    //counterclockwise order.
    //The predicates indicates whether the planar projection of point q lies in this
    //portion of the plane.
    //Preconditions:
    //  o_prime,o,p1 are not collinear
    //  o_prime,o,p2 are not collinear
    //  o_prime,o,q are not collinear
    //  o_prime,o,p1,q are not coplanar or coplanar_orientation(o,o_prime,p1,q)==NEGATIVE
    //  o_prime,o,p2,q are not coplanar or coplanar_orientation(o,o_prime,p2,q)==NEGATIVE
    // template <class Kernel>
    // bool  sorted_around_edge(
    //  const typename Kernel::Point_3& o_prime, const typename Kernel::Point_3& o,
    //  const typename Kernel::Point_3& p1, const typename Kernel::Point_3& p2,
    //  const typename Kernel::Point_3& q)
    auto less = [&tref, &pref, &points, &triangles, &p0, &p1](const std::pair<int, TID>& e1,
                                                              const std::pair<int, TID>& e2)
    {
      using EPECK = CGAL::Exact_predicates_exact_constructions_kernel;
      // handling coplanar
      //TODO: we will make several time the comparision with the same point, not sure if it is coslty. Maybe group them?
      if (triangles[e1.second][e1.first]==triangles[e2.second][e2.first])
        return e1.second < e2.second;
      if (triangles[e1.second][e1.first]==triangles[tref.second][tref.first])
        return true;
      if (triangles[e2.second][e2.first]==triangles[tref.second][tref.first])
        return false;

      return pred::sorted_around_edge<EPECK>(p0, p1, pref,
                                             points[triangles[e2.second][e2.first]],
                                             points[triangles[e1.second][e1.first]]);
    };

#ifdef CGAL_BO_AUTOREF_DEBUG
    std::set<std::size_t> third_point;
    for (const std::pair<int, TID>& p : to_sort)
      third_point.insert(triangles[p.second][p.first]);

    if (third_point.size()!=to_sort.size()+1)
      std::cerr << "  presence of coplanar triangles\n";
#endif

    std::sort(std::next(to_sort.begin()), to_sort.end(), less);
    boost::dynamic_bitset<> is_inside_mesh(nb_meshes, 0); //indicate the status while turning around the edge

    std::vector<bool> init(nb_meshes, false);
    init[input_ids[tref.second]]=true;
    is_inside_mesh[input_ids[tref.second]]=true;
    for (std::size_t pos=1; pos<to_sort.size(); ++pos)
    {
      const std::pair<int, TID>& p = to_sort[pos];
      if (init[input_ids[p.second]]) continue;
      init[input_ids[p.second]]=true;
      CGAL_assertion(pos!=to_sort.size()-1); // not possible since there must be at least 2 triangles per mesh

      //check if collinear with ref
      if (triangles[p.second][p.first]==triangles[tref.second][tref.first])
      {
        // incompatible orientation -> flip the bit
        if (triangles[p.second][(p.first+1)%3]==triangles[tref.second][(tref.first+1)%3])
          is_inside_mesh.set(input_ids[p.second]);
      }
      else
        // incompatible orientation -> flip the bit
        if (triangles[p.second][(p.first+1)%3]==triangles[tref.second][(tref.first+1)%3])
          is_inside_mesh.set(input_ids[p.second]);
    }
#ifdef CGAL_BO_AUTOREF_DEBUG
    auto print_t_info = [&](std::size_t tid)
    {
      std::cout << "         " << tid << " from mesh " << input_ids[tid] << " part of CC " << ccids[tid] << " "
                << points[triangles[tid][0]] << " " << points[triangles[tid][1]] << " " << points[triangles[tid][2]] << "\n";
    };

    int nb=-1;
    for (const std::pair<int, TID>& p : to_sort)
    {
      std::ofstream debug("case_"+ std::to_string(case_id)+"_"+std::to_string(++nb)+".polylines.txt");
      debug << 4 << " " << points[triangles[p.second][0]] << " "<< " " << points[triangles[p.second][1]] << " " << points[triangles[p.second][2]] << " " << points[triangles[p.second][0]] << "\n";
      debug.close();
      print_t_info(p.second);
    }
    std::cout << "  is_inside_mesh = ";
    for (unsigned int i=0; i<nb_meshes; ++i) std::cout << is_inside_mesh[i];
    std::cout <<"\n";
#endif

    CGAL_assertion(triangles[tref.second][tref.first]!=triangles[to_sort.back().second][to_sort.back().first]);// coplanar with tref are smaller than the others
    for (std::size_t pos=0; pos<to_sort.size(); )
    {
      std::size_t start_pos = pos;
      const std::pair<int, TID>& p = to_sort[pos];
      ++pos;
      std::vector<unsigned int> mesh_ids;
      mesh_ids.push_back(input_ids[p.second]);

      while (pos<to_sort.size() && triangles[p.second][p.first]==triangles[to_sort[pos].second][to_sort[pos].first])
      {
        mesh_ids.push_back(input_ids[to_sort[pos].second]);
        coplanar_patches[ccids[to_sort[pos].second]]=true;
        ++pos;
      }
      if (mesh_ids.size()>1)
      {
        coplanar_patches[ccids[p.second]]=true;
        std::sort(mesh_ids.begin(), mesh_ids.end());

        //check orientation -> mark only if all the orientations are identical
        bool same_orientation=true;
        for (std::size_t k=start_pos+1;k<pos; ++k)
        {
          if (triangles[to_sort[start_pos].second][(to_sort[start_pos].first+1)%3] !=
              triangles[to_sort[k].second][(to_sort[k].first+1)%3])
          {
            same_orientation=false;
            break;
          }
        }
        if (same_orientation)
        {
          for (std::size_t k=start_pos;k<pos; ++k)
          {
            coplanar_patches_for_union_and_intersection[ccids[to_sort[k].second]]=true;
          }
        }
      }
      for (auto mid = is_inside_mesh.find_first();
                mid < is_inside_mesh.npos;
                mid = is_inside_mesh.find_next(mid))
      {
        if (!std::binary_search(mesh_ids.begin(), mesh_ids.end(), mid))
          if (is_inside_mesh[mid])
          {
#ifdef CGAL_BO_AUTOREF_DEBUG
            CGAL_assertion( !cc_is_inside_per_mesh_init[mid][ccids[p.second]] || cc_is_inside_per_mesh[mid][ccids[p.second]]==true);
            cc_is_inside_per_mesh_init[mid][ccids[p.second]]=true;
#endif
            cc_is_inside_per_mesh[mid][ccids[p.second]]=true;
#ifdef CGAL_BO_AUTOREF_DEBUG
            std::cout << "CC " << ccids[p.second] << "(" << input_ids[p.second] << ")" << " is inside mesh " << mid << "\n";
#endif
          }
      }
      for (unsigned int mesh_id : mesh_ids)
        is_inside_mesh.flip(mesh_id);
    }
  }

  if (nb_cc!=nb_handled)
  {
    std::vector< std::vector<std::size_t> > triangles_of_intersect_per_cc(nb_cc);
    std::vector< std::pair<Exact_point, std::size_t> > points_and_cc;
    points_and_cc.reserve(nb_cc-nb_handled);
    for (std::size_t cc=0; cc<nb_cc; ++cc)
    {
      if (!cc_handled[cc])
      {
        points_and_cc.emplace_back(points[triangles[one_triangle_per_cc[cc]][0]], cc);
      }
    }

    //TODO: we can also filter using
    std::vector<CGAL::Bbox_3> bb_per_cc(nb_cc);
    for (std::size_t tid=0; tid<triangles.size(); ++tid)
    {
      CGAL::Bbox_3 t_bbox;
      t_bbox+=points[triangles[tid][0]].bbox();
      t_bbox+=points[triangles[tid][1]].bbox();
      t_bbox+=points[triangles[tid][2]].bbox();
      bb_per_cc[ccids[tid]]+=t_bbox;

      for (auto [pt, cc] : points_and_cc)
      {
        if (input_ids[tid]==input_ids[one_triangle_per_cc[cc]]) continue;

        //TODO: we might choose a direction depending on the bbox extents (true choose should be based on the number of triangles)
        //      for now we will consider a ray from the point toward z = +infinity
        if (t_bbox.xmin() > pt.x() || t_bbox.xmax() < pt.x() ||
            t_bbox.ymin() > pt.y() || t_bbox.ymax() < pt.y() ||
            t_bbox.zmax() < pt.z() ) continue;
        triangles_of_intersect_per_cc[cc].push_back(tid);
      }
    }

    //bbox of 2 CCs that are not nested --> cc are disjoint (because points are shared)
    //TODO: this part works as is only for 2 meshes
    for (std::size_t cc=0; cc<nb_cc; ++cc)
    {
      if (!cc_handled[cc])
      {
        bool all_disjoint = true;
        for(std::size_t other_cc=0; other_cc<nb_cc; ++other_cc)
        {
          if (cc!=other_cc && input_ids[one_triangle_per_cc[cc]]!=input_ids[one_triangle_per_cc[other_cc]])
          {
            if (do_overlap(bb_per_cc[cc], bb_per_cc[other_cc]))
            {
              all_disjoint=false;
              break;
            }
          }
        }
        if (all_disjoint) cc_handled[cc]=true;
      }
    }


    for (std::size_t cc=0; cc<nb_cc; ++cc)
    {
      if (!cc_handled[cc])
      {
        // TODO: this is a bit naive but works. To begin with, we can use a better do-intersect + use bbox to filter distances
        std::vector<std::size_t>& subtri=triangles_of_intersect_per_cc[cc];
        const Exact_point& pt = points[triangles[one_triangle_per_cc[cc]][0]];
        Exact_predicates_exact_constructions_kernel::Ray_3 ray(pt, Exact_predicates_exact_constructions_kernel::Vector_3(0.,0.,1.));
        using T3 = Exact_predicates_exact_constructions_kernel::Triangle_3;
        std::vector<std::pair<std::size_t,T3>> to_sort;
        for (std::size_t tid : subtri)
        {
          T3 tr(points[triangles[tid][0]],points[triangles[tid][1]],points[triangles[tid][2]]);
          if (do_intersect(ray, tr))
            to_sort.emplace_back(tid, tr);
        }

        std::sort(to_sort.begin(), to_sort.end(),
                  [&pt](const std::pair<std::size_t, T3>& p1, const std::pair<std::size_t, T3>& p2)
                  {
                    return compare_distance(pt, p1.second, p2.second) == CGAL::SMALLER;
                  }
        );

        bool is_outside = true;
        for (const std::pair<std::size_t, T3>& p : to_sort)
        {
          CGAL::Orientation ori = orientation(p.second[0], p.second[1], p.second[2], pt);
          if (ori == CGAL::COPLANAR) continue;
          is_outside = (ori == POSITIVE);

          break;
        }

        cc_handled[cc]=true;
        if (!is_outside)
          // TODO: clearly working only for 2 meshes
          cc_is_inside_per_mesh[ (input_ids[one_triangle_per_cc[cc]]+1)%2 ][cc]=true;
      }
    }
  }


  data_out.input_ids.swap(input_ids);
#ifdef CGAL_BO3_SELECT_TIMERS
  std::cout << "select done in " << t.time() << "\n";
#endif
  return data_out;
}

struct Autoref_visitor
{
  Autoref_visitor(std::vector<unsigned int>& input_ids,
                  std::size_t nbt1, std::size_t nbt2)
    : input_ids(input_ids)
    , nbt1(nbt1)
    , nbt2(nbt2)
  {}

  inline void number_of_output_triangles(std::size_t nbt)
  {
    input_ids.resize(nbt, -1);
  }
  inline void verbatim_triangle_copy(std::size_t tgt_id, std::size_t src_id)
  {
    input_ids[tgt_id] = src_id < nbt1? 0 : 1;
  }
  inline void new_subtriangle(std::size_t tgt_id, std::size_t src_id)
  {
    input_ids[tgt_id] = src_id < nbt1? 0 : 1;
  }
  inline void delete_triangle(std::size_t /*src_id*/) {}

  std::vector<unsigned int>& input_ids;
  std::size_t nbt1;
  std::size_t nbt2;
};

template <class Concurrency_tag = Sequential_tag, class PointRange, class TriangleRange>
Selection_data
prepare_data(const PointRange& points_1, const TriangleRange& triangles_1,
             const PointRange& points_2, const TriangleRange& triangles_2,
                   PointRange& out_points,     TriangleRange& out_triangles)
{
  namespace PMP = CGAL::Polygon_mesh_processing;

#ifdef CGAL_BO3_TIMERS
  Real_timer timer;
  timer.start();
#endif
  out_points.reserve(points_1.size()+points_2.size());
  out_points.insert(out_points.end(), points_1.begin(), points_1.end());
  out_points.insert(out_points.end(), points_2.begin(), points_2.end());
  out_triangles.reserve(triangles_1.size()+triangles_2.size());
  out_triangles.insert(out_triangles.end(), triangles_1.begin(), triangles_1.end());
  std::size_t nbt_1=triangles_1.size(), nbt_2=triangles_2.size();

  for (auto t : triangles_2)
  {
    for (int i=0;i<3;++i)
      t[i]+=points_1.size();
    out_triangles.push_back(t);
  }

  std::vector<unsigned int> input_ids(nbt_1+nbt_2, 0);
  for (std::size_t i=nbt_1; i<nbt_1+nbt_2; ++i)
    input_ids[i]=1;
  Autoref_visitor visitor(input_ids, nbt_1, nbt_2);

#ifdef CGAL_BO3_TIMERS
  timer.stop();
  std::cout << "merged input done in " << timer.time() << "\n";
  timer.reset();
  timer.start();
#endif

  //TODO:  investigate if a detection of identical triangles in the input could be interesting
  //TODO: import optimisation on the triangle ranges with bboxes meshes intersection to restrict the number of intersection tests
  //TODO: import identical patches opti
  PMP::autorefine_triangle_soup(out_points, out_triangles, CGAL::parameters::concurrency_tag(Concurrency_tag())
                                                                            .visitor(visitor).split_triangle_range_at(nbt_1));
#ifdef CGAL_BO3_TIMERS
  timer.stop();
  std::cout << "autorefine done in " << timer.time() << "\n";
  timer.reset();
  timer.start();
#endif
#ifdef CGAL_BO_AUTOREF_DEBUG
  CGAL::IO::write_OFF("autorefined.off", out_points, out_triangles, CGAL::parameters::stream_precision(17));
#endif
  auto res =select<Concurrency_tag>(out_points, out_triangles, input_ids, nbt_1, nbt_1+nbt_2, 2);
#ifdef CGAL_BO3_TIMERS
  timer.stop();
  std::cout << "select done in " << timer.time() << "\n";
#endif
  return res;
}

template <class Concurrency_tag = Sequential_tag, class PointRange, class TriangleRange>
void
update_output(PointRange& out_points,  TriangleRange& out_triangles, std::size_t nbt, const std::vector<bool>& to_keep, const std::vector<bool>& to_flip=std::vector<bool>())
{
  if (nbt==out_triangles.size() && to_flip.empty()) return;

  TriangleRange triangles;
  triangles.reserve(nbt);

  PointRange points;
  points.reserve(nbt<3?3:((nbt-2)/2+1));
  std::vector<std::size_t> point_map(out_points.size(), out_points.size());

  auto point_id = [&](std::size_t i)
  {
    if (point_map[i]==out_points.size())
    {
      point_map[i]=points.size();
      points.push_back(out_points[i]);
    }
    return point_map[i];
  };

  for (std::size_t ti=0; ti<out_triangles.size(); ++ti)
  {
    if (to_keep[ti])
    {
      auto t = out_triangles[ti];
      t[0]=point_id(t[0]);
      t[1]=point_id(t[1]);
      t[2]=point_id(t[2]);
      if (to_flip.empty() || !to_flip[ti])
        triangles.push_back(t);
      else
      {
        std::swap(t[0], t[1]);
        triangles.push_back(t);
      }
    }
  }

  out_triangles.swap(triangles);
  out_points.swap(points);
}

} // end of boolops_3 namespace

namespace Polygon_mesh_processing
{

#ifndef DOXYGEN_RUNNING
template <class PointRange, class TriangleRange, class NamedParameters = parameters::Default_named_parameters>
void compute_self_union(PointRange& points, TriangleRange& triangles,
                        const NamedParameters& np = parameters::default_values())
{
  CGAL_USE(np);

  typedef typename internal_np::Lookup_named_param_def <
    internal_np::concurrency_tag_t,
    NamedParameters,
    Sequential_tag
  > ::type Concurrency_tag;

  namespace PMP = CGAL::Polygon_mesh_processing;
  using namespace boolops_3;

#ifdef CGAL_BO3_TIMERS
  Real_timer timer;
  timer.start();
#endif
  // TODO: convert points to EPECK if not already the case
  using Exact_point = typename std::iterator_traits<typename PointRange::const_iterator>::value_type;
  PMP::autorefine_triangle_soup(points, triangles, CGAL::parameters::concurrency_tag(Concurrency_tag()));
#ifdef CGAL_BO3_TIMERS
  timer.stop();
  std::cout << "autorefine done in " << timer.time() << "\n";
  timer.reset();
#endif
#ifdef CGAL_BO_AUTOREF_DEBUG
  CGAL::IO::write_OFF("autorefined.off", out_points, out_triangles, CGAL::parameters::stream_precision(17));
#endif

  //TODO: use Concurrency_tagf
  namespace PMP = CGAL::Polygon_mesh_processing;
  namespace pred = PMP::Corefinement;
  //~ Selection_data data_out;

  // fill edge map (TODO: move to a function)
  using PID = std::size_t;
  using EID = std::size_t;
  using TID = std::size_t;
  using Edge_record = std::tuple<PID, PID, TID>;

  std::vector<Edge_record> all_triangle_edges(3 * triangles.size());
  for (TID tid = 0; tid < triangles.size(); ++tid)
  {
    auto add_edge = [&](std::size_t idx, PID a, PID b)
    {
      if (b < a) std::swap(a, b);
      all_triangle_edges[idx] = Edge_record{a, b, tid};
    };
    add_edge(3 * tid, triangles[tid][0], triangles[tid][1]);
    add_edge(3 * tid + 1, triangles[tid][1], triangles[tid][2]);
    add_edge(3 * tid + 2, triangles[tid][2], triangles[tid][0]);
  }
  std::sort(all_triangle_edges.begin(), all_triangle_edges.end());

  std::vector<std::array<EID, 3>> adjacent_edges(triangles.size());
  std::vector<boost::container::small_vector<TID, 2>> adjacent_faces;
  adjacent_faces.reserve(all_triangle_edges.size());
  for (std::size_t i = 0; i < all_triangle_edges.size(); )
  {
    const PID a = std::get<0>(all_triangle_edges[i]);
    const PID b = std::get<1>(all_triangle_edges[i]);
    const EID eid = adjacent_faces.size();
    adjacent_faces.emplace_back();
    do
    {
      const TID tid = std::get<2>(all_triangle_edges[i++]);
      adjacent_faces.back().push_back(tid);
      for (std::size_t j = 0; j < 3; ++j)
        if ((triangles[tid][j] == a && triangles[tid][(j+1)%3] == b) ||
            (triangles[tid][j] == b && triangles[tid][(j+1)%3] == a))
        {
          adjacent_edges[tid][(j+2)%3] = eid;
          break;
        }
    } while (i < all_triangle_edges.size() && std::get<0>(all_triangle_edges[i]) == a && std::get<1>(all_triangle_edges[i]) == b);
  }

  // extract connected_components bounded by non-manifold edges
  std::vector<std::size_t> ccids;
  std::vector<std::size_t> one_triangle_per_cc =
    boolops_3::connected_components(adjacent_edges, adjacent_faces, ccids);
  std::size_t nb_cc = one_triangle_per_cc.size();

  std::vector<EID> nm_edges;
  for (EID eid=0; eid<adjacent_faces.size(); ++eid)
    if (adjacent_faces[eid].size()>2)
      nm_edges.push_back(eid);

  // Get the vertices associate to an edge
  auto edge_vertices = [&](EID eid) {
    std::size_t tid = adjacent_faces[eid].front();
    std::size_t k = adjacent_edges[tid][0] == eid ? 0 :adjacent_edges[tid][1] == eid ? 1 : 2;
    CGAL_assertion(adjacent_edges[tid][k] == eid);
    return std::make_pair(triangles[tid][(k+1)%3], triangles[tid][(k+2)%3]);
  };


  using UF_CC = CGAL::Union_find<std::size_t>;
  using UF_handle = UF_CC::handle;

  UF_CC uf_cc;
  std::vector<UF_handle> uf_handles(nb_cc);
  for (std::size_t i=0; i<nb_cc; ++i)
    uf_handles[i]=uf_cc.make_set(i);


  // cc -> edges
  // TODO we need split graph into polylines I guess
  std::vector<std::vector<EID> > cc_to_edges(nb_cc);
  std::vector<bool> pts_on_edge(points.size(), false); // TODO: fill it in autoref visitor?

  for (EID eid : nm_edges)
  {
    std::set<std::size_t> ccs;
    const auto [pid, other_pid] = edge_vertices(eid);
    pts_on_edge[pid]=true;
    pts_on_edge[other_pid]=true;
    for (TID tid : adjacent_faces[eid])
      ccs.insert(ccids[tid]);
    for (std::size_t ccid : ccs)
      cc_to_edges[ccid].push_back(eid);
  }

  std::size_t min_pt_id=0;
  for (std::size_t i=1; i<points.size(); ++i)
  {
    if (points[min_pt_id].x()==points[i].x())
    {
      if (points[min_pt_id].y()==points[i].y())
      {
        CGAL_assertion(points[min_pt_id].z()!=points[i].z());
        if (points[min_pt_id].z()>points[i].z()) min_pt_id=i;
      }
      else
        if (points[min_pt_id].y()>points[i].y()) min_pt_id=i;
    }
    else
      if (points[min_pt_id].x()>points[i].x()) min_pt_id=i;
  }

  //TODO: temporary condition assuming extreme point is manifold, need another selection method otherwise
  CGAL_precondition(!pts_on_edge[min_pt_id]);

  std::size_t selected_tid=triangles.size();
  for (TID tid=0; tid<triangles.size(); ++tid)
  {
    for (int i=0;i<3;++i)
    if (triangles[tid][i]==min_pt_id)
    {
      selected_tid=tid;
      break;
    }
  }

  //TODO: temporary condition, need another selection method is we have isolated vertices
  CGAL_precondition(selected_tid!=triangles.size());

  std::vector<bool> ccs_handled(nb_cc, false);
  std::vector<std::size_t> cc_stack;
  cc_stack.push_back(ccids[selected_tid]);

  // main loop to classify cc
  while(!cc_stack.empty())
  {
    std::size_t ccid = cc_stack.back();
    cc_stack.pop_back();

    if (ccs_handled[ccid]) continue;
    ccs_handled[ccid]=true;

    // iterate on all non-manifold edges of the current cc
    for (EID eid : cc_to_edges[ccid])
    {
      const auto [pid, other_pid] = edge_vertices(eid);
#ifdef CGAL_BO_AUTOREF_DEBUG
      std::cout << "case #" << ++case_id << " looking at edge " << pid << " " << other_pid << "\n";
      std::cout << "  incident CCs:";
      for (TID tid : adjacent_faces[eid]) std::cout << " " << ccids[tid] << "("<< input_ids[tid] <<")";
      std::cout << "\n";
#endif

      std::vector<std::pair<int, TID>> to_sort;
      to_sort.reserve(adjacent_faces[eid].size());

      std::size_t i0=pid, i1=other_pid;
      for (TID tid : adjacent_faces[eid])
      {
        for (int i=0; i<3; ++i)
        {
          if (triangles[tid][i]!=i0 && triangles[tid][i]!=i1)
            to_sort.emplace_back(i, tid);
        }
        if (ccids[tid]==ccid && to_sort.size()>1)
          std::swap(to_sort.front(), to_sort.back());
      }

      std::pair<int, TID> tref=to_sort.front(); // tid is from ccid
      if (triangles[tref.second][(tref.first+1)%3]!=i0)
        std::swap(i0,i1);

      const Exact_point& pref=points[triangles[tref.second][tref.first]];
      const Exact_point& p0=points[i0];
      const Exact_point& p1=points[i1];

      //Considering the plane with normal vector [o_prime,o] and containing o.
      //We define the counterclockwise order around o when looking from
      //the side of the plane containing o_prime.
      //We consider the portion of the plane defined by rotating a ray starting at o
      //from the planar projection of p1 to the planar projection of p2 in
      //counterclockwise order.
      //The predicates indicates whether the planar projection of point q lies in this
      //portion of the plane.
      //Preconditions:
      //  o_prime,o,p1 are not collinear
      //  o_prime,o,p2 are not collinear
      //  o_prime,o,q are not collinear
      //  o_prime,o,p1,q are not coplanar or coplanar_orientation(o,o_prime,p1,q)==NEGATIVE
      //  o_prime,o,p2,q are not coplanar or coplanar_orientation(o,o_prime,p2,q)==NEGATIVE
      // template <class Kernel>
      // bool  sorted_around_edge(
      //  const typename Kernel::Point_3& o_prime, const typename Kernel::Point_3& o,
      //  const typename Kernel::Point_3& p1, const typename Kernel::Point_3& p2,
      //  const typename Kernel::Point_3& q)
      auto less = [&tref, &pref, &points, &triangles, &p0, &p1](const std::pair<int, TID>& e1,
                                                                const std::pair<int, TID>& e2)
      {
        using EPECK = CGAL::Exact_predicates_exact_constructions_kernel;
        // handling coplanar
        //TODO: we will make several time the comparision with the same point, not sure if it is coslty. Maybe group them?
        if (triangles[e1.second][e1.first]==triangles[e2.second][e2.first])
          return e1.second < e2.second;
        if (triangles[e1.second][e1.first]==triangles[tref.second][tref.first])
          return true;
        if (triangles[e2.second][e2.first]==triangles[tref.second][tref.first])
          return false;

        return pred::sorted_around_edge<EPECK>(p0, p1, pref,
                                               points[triangles[e2.second][e2.first]],
                                               points[triangles[e1.second][e1.first]]);
      };

  #ifdef CGAL_BO_AUTOREF_DEBUG
      std::set<std::size_t> third_point;
      for (const std::pair<int, TID>& p : to_sort)
        third_point.insert(triangles[p.second][p.first]);

      if (third_point.size()!=to_sort.size()+1)
        std::cerr << "  presence of coplanar triangles\n";
  #endif

      std::sort(std::next(to_sort.begin()), to_sort.end(), less);

      // TODO: we need to make sure orentations are compatible

      std::size_t n_ccid = ccids[to_sort[1].second];
      uf_handles[n_ccid] = uf_cc.find(uf_handles[n_ccid]);
      uf_handles[ccid] = uf_cc.find(uf_handles[ccid]);
      uf_cc.unify_sets(uf_handles[n_ccid], uf_handles[ccid]);
      if (!ccs_handled[n_ccid])
        cc_stack.push_back(n_ccid);
    }
  }

  // mark ccids to keep
  std::size_t ref_cc = ccids[selected_tid];
  UF_handle ref_handle = uf_cc.find(uf_handles[ref_cc]);
  std::vector<bool> selected_ccs(nb_cc);
  for (std::size_t ccid=0; ccid<nb_cc; ++ccid)
  {
    UF_handle cc_h = uf_cc.find(uf_handles[ccid]);
    selected_ccs[ccid] = (cc_h==ref_handle);
  }

  std::vector<bool> to_keep(triangles.size(), false);
  std::size_t nb_out_triangles=0;
  for (std::size_t tid=0; tid<triangles.size(); ++tid)
  {
    if (selected_ccs[ccids[tid]])
    {
      to_keep[tid]=true;
      ++nb_out_triangles;
    }
  }
  boolops_3::update_output<Concurrency_tag>(points, triangles, nb_out_triangles, to_keep);
}
#endif

/**
  * \ingroup PMP_boolop_soup_grp
  *
  * computes the union of two volumes, each bounded by a triangle soup.
  *
  * \pre \link CGAL::Polygon_mesh_processing::does_triangle_soup_self_intersect() `!CGAL::Polygon_mesh_processing::does_triangle_soup_self_intersect(points_1, triangles_1)` \endlink
  * \pre \link CGAL::Polygon_mesh_processing::does_triangle_soup_self_intersect() `!CGAL::Polygon_mesh_processing::does_triangle_soup_self_intersect(points_2, triangles_2)` \endlink
  *
  * @tparam PointRange a model of the concepts `RandomAccessContainer` whose value type is `GeomTraits::Point_3` as value type, `GeomTraits` being the type of the named parameter `geom_traits`.
  * @tparam TriangleRange a model of the concept `RandomAccessContainer` whose `value_type` is a model of the concept `RandomAccessContainer` whose `value_type` is `std::size_t`.
  *
  * @param points_1 points of the first triangle soup.
  * @param triangles_1 triangles of the first triangle soup, each element in the vector describes a triangle using the indices of the points in `points_1`.
  *
  * @param points_2 points of the second triangle soup.
  * @param triangles_2 triangles of the second triangle soup, each element in the vector describes a triangle using the indices of the points in `points_2`.
  *
  * @param out_points points of the output triangle soup bounding the volume result of the Boolean operation.
  * @param out_triangles triangles of the ouput triangle soup, each element in the vector describes a triangle using the indices of the points in `out_points`.
  *
  * @param np optional sequence of named parameters among the ones listed below
  *
  * \cgalNamedParamsBegin
  *   \cgalParamNBegin{concurrency_tag}
  *     \cgalParamDescription{a tag indicating if the task should be done using one or several threads.}
  *     \cgalParamType{Either `CGAL::Sequential_tag`, or `CGAL::Parallel_tag`, or `CGAL::Parallel_if_available_tag`}
  *     \cgalParamDefault{`CGAL::Sequential_tag`}
  *   \cgalParamNEnd
  *   \cgalParamNBegin{geom_traits}
  *     \cgalParamDescription{an instance of a geometric traits class}
  *     \cgalParamType{a class model of `Kernel`}
  *     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
  *     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
  *   \cgalParamNEnd
  * \cgalNamedParamsEnd
  */
template <class PointRange, class TriangleRange, class NamedParameters = parameters::Default_named_parameters>
void compute_union(const PointRange& points_1, const TriangleRange& triangles_1,
                   const PointRange& points_2, const TriangleRange& triangles_2,
                         PointRange& out_points,     TriangleRange& out_triangles,
                   const NamedParameters& np = parameters::default_values())
{
  CGAL_USE(np);

  typedef typename internal_np::Lookup_named_param_def <
    internal_np::concurrency_tag_t,
    NamedParameters,
    Sequential_tag
  > ::type Concurrency_tag;

#ifdef CGAL_BO3_TIMERS
  Real_timer t;
  t.start();
#endif

  using namespace boolops_3;
  Selection_data data_out=prepare_data<Concurrency_tag>(points_1, triangles_1, points_2, triangles_2, out_points, out_triangles);
  const unsigned int nb_meshes=2;

#ifdef CGAL_BO3_TIMERS
  std::cout << "Prepare_data + autoref: " << t.time() << "\n";
#endif
  std::vector<bool> to_keep(out_triangles.size(), false);
  std::size_t nb_out_triangles=0;

  for (std::size_t ti=0; ti<out_triangles.size(); ++ti)
  {
    unsigned int t_mesh_id=data_out.input_ids[ti];
    bool keep=true;
    //TODO works only for 2 meshes
    if (data_out.coplanar_patches[data_out.ccids[ti]])
    {
      if (!data_out.coplanar_patches_for_union_and_intersection[data_out.ccids[ti]] || t_mesh_id==1)
      {
        keep=false;
      }
    }
    else
      for(unsigned int mi=0; mi<nb_meshes; ++mi)
      {
        if (mi!=t_mesh_id && data_out.cc_is_inside_per_mesh[mi][data_out.ccids[ti]])
        {
          keep=false;
          break;
        }
      }
    if (keep)
    {
      ++nb_out_triangles;
      to_keep[ti]=true;
    }
  }

#ifdef CGAL_BO3_TIMERS
  std::cout << "Collect patches: " << t.time() << "\n";
#endif
  boolops_3::update_output<Concurrency_tag>(out_points, out_triangles, nb_out_triangles, to_keep);
#ifdef CGAL_BO3_TIMERS
  std::cout << "BO3 Building output: " << t.time() << "\n";
#endif
}

/**
  * \ingroup PMP_boolop_soup_grp
  * computes the intersection of two volumes, each bounded by a triangle soup.
  * \copydetails CGAL::Polygon_mesh_processing::compute_union()
  */
template <class PointRange, class TriangleRange, class NamedParameters = parameters::Default_named_parameters>
void compute_intersection(const PointRange& points_1, const TriangleRange& triangles_1,
                          const PointRange& points_2, const TriangleRange& triangles_2,
                                PointRange& out_points,     TriangleRange& out_triangles,
                          const NamedParameters& np = parameters::default_values())
{
  CGAL_USE(np);

  typedef typename internal_np::Lookup_named_param_def <
    internal_np::concurrency_tag_t,
    NamedParameters,
    Sequential_tag
  > ::type Concurrency_tag;

  using namespace boolops_3;
  Selection_data data_out=prepare_data<Concurrency_tag>(points_1, triangles_1, points_2, triangles_2, out_points, out_triangles);
  const unsigned int nb_meshes=2;

  std::vector<bool> to_keep(out_triangles.size(), false);
  std::size_t nb_out_triangles=0;

  for (std::size_t ti=0; ti<out_triangles.size(); ++ti)
  {
    unsigned int t_mesh_id=data_out.input_ids[ti];
    bool keep=true;
    //TODO works only for 2 meshes
    if (data_out.coplanar_patches[data_out.ccids[ti]])
    {
      if (!data_out.coplanar_patches_for_union_and_intersection[data_out.ccids[ti]] || t_mesh_id==1)
      {
        keep=false;
      }
    }
    else
      for(unsigned int mi=0; mi<nb_meshes; ++mi)
      {
        if (mi!=t_mesh_id && !data_out.cc_is_inside_per_mesh[mi][data_out.ccids[ti]])
        {
          keep=false;
          break;
        }
      }

    if (keep)
    {
      ++nb_out_triangles;
      to_keep[ti]=true;
    }
  }

  boolops_3::update_output<Concurrency_tag>(out_points, out_triangles, nb_out_triangles, to_keep);
}

/**
  * \ingroup PMP_boolop_soup_grp
  * computes the difference of two volumes, each bounded by a triangle soup.
  * \copydetails CGAL::Polygon_mesh_processing::compute_union()
  */
template <class PointRange, class TriangleRange, class NamedParameters = parameters::Default_named_parameters>
void compute_difference(const PointRange& points_1, const TriangleRange& triangles_1,
                        const PointRange& points_2, const TriangleRange& triangles_2,
                              PointRange& out_points,     TriangleRange& out_triangles,
                        const NamedParameters& np = parameters::default_values())
{
  CGAL_USE(np);

  typedef typename internal_np::Lookup_named_param_def <
    internal_np::concurrency_tag_t,
    NamedParameters,
    Sequential_tag
  > ::type Concurrency_tag;

  using namespace boolops_3;
  Selection_data data_out=prepare_data<Concurrency_tag>(points_1, triangles_1, points_2, triangles_2, out_points, out_triangles);
  // const unsigned int nb_meshes=2;

  std::vector<bool> to_keep(out_triangles.size(), false);
  std::vector<bool> to_flip(out_triangles.size(), false);
  std::size_t nb_out_triangles=0;

  for (std::size_t ti=0; ti<out_triangles.size(); ++ti)
  {

    if (data_out.input_ids[ti]==0)
    {
      if (data_out.coplanar_patches[data_out.ccids[ti]] && data_out.coplanar_patches_for_union_and_intersection[data_out.ccids[ti]]) continue;
      if (!data_out.cc_is_inside_per_mesh[1][data_out.ccids[ti]])
      {
        ++nb_out_triangles;
        to_keep[ti]=true;
        continue;
      }
    }
    else
    {
      if (data_out.coplanar_patches[data_out.ccids[ti]]) continue;
      if (data_out.cc_is_inside_per_mesh[0][data_out.ccids[ti]])
      {
        ++nb_out_triangles;
        to_keep[ti]=true;
        to_flip[ti]=true;
      }
    }
  }

  boolops_3::update_output<Concurrency_tag>(out_points, out_triangles, nb_out_triangles, to_keep, to_flip);
}

} } // CGAL::Polygon_mesh_processing namespace

#endif // CGAL_POLYGON_MESH_PROCESSING_TRIANGLE_SOUP_BOOLEAN_OPERATIONS_3_H
