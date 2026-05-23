// Copyright (c) 2026 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: LicenseRef-Commercial
//
//
// Author(s)     : Sébastien Loriot

#ifndef CGAL_POLYGON_MESH_PROCESSING_TRIANGLE_SOUP_BOOLEAN_OPERATIONS_3_H
#define CGAL_POLYGON_MESH_PROCESSING_TRIANGLE_SOUP_BOOLEAN_OPERATIONS_3_H

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/autorefinement.h>
#include <CGAL/Polygon_mesh_processing/internal/Corefinement/predicates.h>
#include <CGAL/Union_find.h>

#include <fstream>

#ifdef CGAL_BO3_TIMERS
#include <CGAL/Real_timer.h>
#endif

namespace CGAL {
namespace boolops_3 {

//TODO: update nm_edges for identical coplanar status (less test below then)
template <class Polygons, class EdgeMap, class NM_edges>
std::vector<std::size_t>
connected_components(const Polygons& triangles,
                     const EdgeMap& edge_map,
                     std::vector<std::size_t>& ccids,
                     NM_edges& nm_edges)
{
  std::vector<std::size_t> one_triangle_per_cc;

  std::size_t nb_poly = triangles.size();
  ccids.resize(nb_poly);
  std::size_t i=0;

  std::vector<bool> handled(nb_poly, false);
  for (std::size_t tid=0; tid<nb_poly; ++tid)
  {
    if (handled[tid]) continue;
    std::vector<std::size_t> queue;
    queue.push_back(tid);
    while(!queue.empty())
    {
      std::size_t q_tid = queue.back();
      queue.pop_back();
      if ( handled[q_tid]) continue;
      handled[q_tid]=true;
      ccids[q_tid]=i;


      for (int k=0; k<3; ++k) // TODO: change that if you want more than triangles
      {
        auto p = CGAL::make_sorted_pair(triangles[q_tid][k], triangles[q_tid][(k+1)%3]);
        auto it = edge_map[p.first].find(p.second);
        CGAL_assertion(it!=edge_map[p.first].end());
        const auto& tids =  it->second;
        switch( tids.size() )
        {
          case 1:
          break;
          case 2:
          {
            std::size_t n_tid = tids[0]!=q_tid ? tids[0] : tids[1];
            if (!handled[n_tid])
              queue.push_back(n_tid);
          }
          break;
          default:
            nm_edges.emplace_back(p.first, it);
          //nm edge
          break;
        }
      }
    }
    one_triangle_per_cc.push_back(tid);
    ++i;
  }
  return one_triangle_per_cc;
}

struct Selection_data
{
  std::vector<std::size_t> ccids;
  std::vector<bool> coplanar_patches;
  std::vector<bool> coplanar_patches_for_union_and_intersection;
  std::vector< std::vector<bool> > cc_is_inside_per_mesh; // TODO: use an optional
  std::vector<unsigned int> input_ids;
};

template <class Concurrency_tag, class Exact_point, class Triangle>
Selection_data
select(const std::vector<Exact_point>& points,
       const std::vector<Triangle>& triangles,
             std::vector<unsigned int>& input_ids,
       const unsigned int nb_meshes)
{
  //TODO: use Concurrency_tag
  namespace PMP = CGAL::Polygon_mesh_processing;
  namespace pred = PMP::Corefinement;
  Selection_data data_out;

  // fill edge map
  typedef std::size_t PID;
  typedef std::size_t TID;
  std::vector< std::unordered_map<PID, std::vector<TID> > > edge_map(points.size());
  using Boundary_edge_location = std::pair<PID, std::unordered_map<PID, std::vector<TID> >::const_iterator>;

  for (TID tid=0; tid<triangles.size(); ++tid)
  {
    auto p = CGAL::make_sorted_pair(triangles[tid][0], triangles[tid][1]);
    edge_map[p.first][p.second].push_back(tid);
    p = CGAL::make_sorted_pair(triangles[tid][1], triangles[tid][2]);
    edge_map[p.first][p.second].push_back(tid);
    p = CGAL::make_sorted_pair(triangles[tid][2], triangles[tid][0]);
    edge_map[p.first][p.second].push_back(tid);
  }

  // extract connected_components bounded by non-manifold edges
  std::vector<std::size_t>& ccids = data_out.ccids;
  std::vector< Boundary_edge_location > nm_edges;
  std::vector<std::size_t> one_triangle_per_cc =
    boolops_3::connected_components(triangles, edge_map, ccids, nm_edges);
  std::size_t nb_cc = one_triangle_per_cc.size();

  // TODO: direct filling rather that sort  + unique?
  //~ std::sort(nm_edges.begin(), nm_edges.end(),
            //~ [](const std::pair<PID, PID>& p1, const std::pair<PID, PID>& p2)
            //~ {
              //~ if (p1.first!=p2.first) return p1.first<p2.first;
              //~ return p1.second < p2.second;
            //~ });
  //~ nm_edges.erase(std::unique(nm_edges.begin(), nm_edges.end()), nm_edges.end());

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

  coplanar_patches.assign(nb_cc, false);
  coplanar_patches_for_union_and_intersection.assign(nb_cc, false);
#ifdef CGAL_BO_AUTOREF_DEBUG
  std::vector< std::vector<bool> > cc_is_inside_per_mesh_init(nb_meshes, std::vector<bool>(nb_cc, false));
  int case_id=-1;
#endif
  cc_is_inside_per_mesh.assign(nb_meshes, std::vector<bool>(nb_cc, false)); // TODO: use an optional

  for (auto [pid, it] : nm_edges)
  {
    if (nb_handled==nb_cc) break;

    bool all_handled=true;
    for (TID tid : it->second)
    {
      if (!cc_handled[ccids[tid]])
      {
        all_handled=false;
        break;
      }
    }

    if (all_handled) continue;

    bool all_from_same_mesh = true;

    for (TID tid : it->second)
    {
      if (input_ids[tid]!=input_ids[it->second.front()])
      {
        all_from_same_mesh=false;
        break;
      }
    }
    if (all_from_same_mesh) continue;
#ifdef CGAL_BO_AUTOREF_DEBUG
    std::cout << "case #" << ++case_id << " looking at edge " << pid << " " << it->first << "\n";
    std::cout << "  incident CCs:";
    for (TID tid : it->second) std::cout << " " << ccids[tid] << "("<< input_ids[tid] <<")";
    std::cout << "\n";
#endif
    for (TID tid : it->second)
    {
      if (!cc_handled[ccids[tid]]) ++nb_handled;
      cc_handled[ccids[tid]] = true;
    }

    std::vector<std::pair<int, TID>> to_sort;
    to_sort.reserve(it->second.size());

    std::size_t i0=pid, i1=it->first;
    for (TID tid : it->second)
    {
      for (int i=0; i<3; ++i)
      {
        if (triangles[tid][i]!=i0 && triangles[tid][i]!=i1)
          to_sort.emplace_back(i, tid);
      }
    }

    // we arbitrarily take the last triangle
    std::pair<int, TID> tref=to_sort.front();
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

#ifdef CGAL_BO3_TIMERS
  Real_timer timer;
  timer.start();
#endif
  std::vector<unsigned int> input_ids(nbt_1+nbt_2, 0);
  for (std::size_t i=nbt_1; i<nbt_1+nbt_2; ++i)
    input_ids[i]=1;
  Autoref_visitor visitor(input_ids, nbt_1, nbt_2);
  //TODO:  investigate if a detection of identical triangles in the input could be interesting
  //TODO: import optimisation on the triangle ranges with bboxes meshes intersection to restrict the number of intersection tests
  //TODO: import identical patches opti
  PMP::autorefine_triangle_soup(out_points, out_triangles, CGAL::parameters::concurrency_tag(Concurrency_tag())
                                                                            .visitor(visitor).split_triangle_range_at(nbt_1));
#ifdef CGAL_BO3_TIMERS
  timer.stop();
  std::cout << "autorefine done in " << timer.time() << "\n";
  timer.reset();
#endif
#ifdef CGAL_BO_AUTOREF_DEBUG
  CGAL::IO::write_OFF("autorefined.off", out_points, out_triangles, CGAL::parameters::stream_precision(17));
#endif
  return select<Concurrency_tag>(out_points, out_triangles, input_ids, 2);
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

template <class Concurrency_tag = Sequential_tag, class PointRange, class TriangleRange>
void compute_self_union(PointRange& points, TriangleRange& triangles)
{
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
  typedef std::size_t PID;
  typedef std::size_t TID;
  std::vector< std::unordered_map<PID, std::vector<TID> > > edge_map(points.size());
  using Boundary_edge_location = std::pair<PID, std::unordered_map<PID, std::vector<TID> >::const_iterator>;

  for (TID tid=0; tid<triangles.size(); ++tid)
  {
    auto p = CGAL::make_sorted_pair(triangles[tid][0], triangles[tid][1]);
    edge_map[p.first][p.second].push_back(tid);
    p = CGAL::make_sorted_pair(triangles[tid][1], triangles[tid][2]);
    edge_map[p.first][p.second].push_back(tid);
    p = CGAL::make_sorted_pair(triangles[tid][2], triangles[tid][0]);
    edge_map[p.first][p.second].push_back(tid);
  }

  // extract connected_components bounded by non-manifold edges  (TODO: move to a function)
  std::vector<std::size_t> ccids;
  std::vector< Boundary_edge_location > nm_edges;
  std::vector<std::size_t> one_triangle_per_cc =
    boolops_3::connected_components(triangles, edge_map, ccids, nm_edges);
  std::size_t nb_cc = one_triangle_per_cc.size();


  using UF_CC = CGAL::Union_find<std::size_t>;
  using UF_handle = UF_CC::handle;

  UF_CC uf_cc;
  std::vector<UF_handle> uf_handles(nb_cc);
  for (std::size_t i=0; i<nb_cc; ++i)
    uf_handles[i]=uf_cc.make_set(i);


  // cc -> edges
  // TODO we need split graph into polylines I guess
  std::vector<std::vector<Boundary_edge_location> > cc_to_edges(nb_cc);
  std::vector<bool> pts_on_edge(points.size(), false); // TODO: fill it in autoref visitor?

  for (auto [pid, it] : nm_edges)
  {
    std::set<std::size_t> ccs;
    pts_on_edge[pid]=true;
    pts_on_edge[it->first]=true;
    for (TID tid : it->second)
      ccs.insert(ccids[tid]);
    for (std::size_t ccid : ccs)
      cc_to_edges[ccid].emplace_back(pid, it);
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
    for (auto [pid, it] : cc_to_edges[ccid])
    {
#ifdef CGAL_BO_AUTOREF_DEBUG
      std::cout << "case #" << ++case_id << " looking at edge " << pid << " " << it->first << "\n";
      std::cout << "  incident CCs:";
      for (TID tid : it->second) std::cout << " " << ccids[tid] << "("<< input_ids[tid] <<")";
      std::cout << "\n";
#endif

      std::vector<std::pair<int, TID>> to_sort;
      to_sort.reserve(it->second.size());

      std::size_t i0=pid, i1=it->first;
      for (TID tid : it->second)
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

template <class Concurrency_tag = Sequential_tag, class PointRange, class TriangleRange>
void compute_union(const PointRange& points_1, const TriangleRange& triangles_1,
                   const PointRange& points_2, const TriangleRange& triangles_2,
                         PointRange& out_points,     TriangleRange& out_triangles)
{
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

  boolops_3::update_output<Concurrency_tag>(out_points, out_triangles, nb_out_triangles, to_keep);
}

template <class Concurrency_tag = Sequential_tag, class PointRange, class TriangleRange>
void compute_intersection(const PointRange& points_1, const TriangleRange& triangles_1,
                          const PointRange& points_2, const TriangleRange& triangles_2,
                                PointRange& out_points,     TriangleRange& out_triangles)
{
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

template <class Concurrency_tag = Sequential_tag, class PointRange, class TriangleRange>
void compute_difference(const PointRange& points_1, const TriangleRange& triangles_1,
                        const PointRange& points_2, const TriangleRange& triangles_2,
                              PointRange& out_points,     TriangleRange& out_triangles)
{
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
