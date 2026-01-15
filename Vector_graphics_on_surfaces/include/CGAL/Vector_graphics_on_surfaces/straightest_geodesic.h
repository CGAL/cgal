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

#ifndef CGAL_POLYGON_MESH_PROCESSING_BSURF_STRAIGHTEST_GEODESIC_H
#define CGAL_POLYGON_MESH_PROCESSING_BSURF_STRAIGHTEST_GEODESIC_H

#include <CGAL/license/Vector_graphics_on_surfaces.h>

#include <CGAL/Vector_graphics_on_surfaces/internal/utils.h>

namespace CGAL {
namespace Vector_graphics_on_surfaces {

namespace internal {

template <class K, class TriangleMesh, class VertexPointMap>
struct Straightest_geodesic_imp
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

  // dir is defined wrt the halfedge of start_loc.first that is used as y axis
  static
  std::vector<CGAL::Polygon_mesh_processing::Face_location<TriangleMesh, FT>>
  straightest_geodesic(const CGAL::Polygon_mesh_processing::Face_location<TriangleMesh, FT>& start_loc,
                       const TriangleMesh& mesh,
                       const VertexPointMap &vpm,
                       Vector_2 dir,const FT& len)
  {
    namespace PMP = CGAL::Polygon_mesh_processing;

#ifdef CGAL_DEBUG_BSURF
    {
      std::ofstream log("/tmp/log.txt");
      log << std::setprecision(17);
      log << start_loc.first << " " << start_loc.second[0]  << " " << start_loc.second[1]  << " " << start_loc.second[2] << "\n";
      log << "dir = " << dir << "\n";
      log << "len = " << len << "\n";
    }
#endif

    if (len == 0) return {start_loc};

    auto get_halfedge_offset=[&mesh](const halfedge_descriptor& h_ref,const halfedge_descriptor& h_curr)
    {
      if( h_ref == h_curr) return 0;
      if( next(h_ref, mesh) == h_curr ) return 1;
      if( prev(h_ref, mesh) == h_curr ) return 2;

      std::cout<<"Error! Halfedges are in different faces"<<std::endl;

      CGAL_assertion(false);
      return -1;
    };

#ifdef CGAL_DEBUG_BSURF
    auto fn = compute_face_normal(start_loc.first, mesh);
    double theta= std::acos(dir.y()/sqrt(dir*dir));
    if(dir.x()<0)
      theta*=-1;
    auto dir3 = rotate_vector(mesh.point(target(halfedge(start_loc.first, mesh), mesh)) -
                              mesh.point(source(halfedge(start_loc.first, mesh), mesh)),
                              fn, theta);

    std::cout << "direction " << PMP::construct_point(start_loc, mesh) << " " << PMP::construct_point(start_loc, mesh)+dir3 << "\n";
#endif

    auto get_vid_offset=[&mesh](const halfedge_descriptor& h_ref,const vertex_descriptor& vid)
    {
      if(source(h_ref,mesh)==vid)
        return 0;

      if(target(h_ref,mesh)==vid)
        return 1;

      if(target(next(h_ref,mesh),mesh)==vid)
        return 2;

      std::cerr<<"Error! Halfedges are in different faces"<<std::endl;

      CGAL_assertion(false);

      return -1;
    };

    auto get_halfedge=[&mesh](const int k,const halfedge_descriptor& h_ref)
    {
      switch(k)
      {
        case 0:
         return h_ref;

        case 1:
         return next(h_ref,mesh);

        default:
         return prev(h_ref,mesh);
      }
    };

    auto edge_barycentric_coordinate =
      [&mesh](const halfedge_descriptor h_edge,
              const halfedge_descriptor h_face,
              const std::array<FT,2>& bary_edge)
    {
      std::array<FT,3> bary_edge_in_face=make_array(0.,0.,0.);
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
          assert(h_face==prev(h_edge,mesh));
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

    std::vector<PMP::Face_location<TriangleMesh, FT>> result;
    FT accumulated=0.;
    face_descriptor curr_tid=start_loc.first;
    std::array<Vector_2, 3> curr_flat_tid=init_flat_triangle<K>(halfedge(curr_tid,mesh),vpm,mesh);
    Vector_2 flat_p= start_loc.second[0]*curr_flat_tid[0]+start_loc.second[1]*curr_flat_tid[1]+start_loc.second[2]*curr_flat_tid[2];
    PMP::Face_location<TriangleMesh, FT> curr_p=start_loc;

    PMP::descriptor_variant<TriangleMesh> var_des = PMP::get_descriptor_from_location(start_loc,mesh);
    switch(var_des.index())
    {
      case 0:
      {
        // TODO: handle is src_vertex is a border vertex: does not necessarily mean early exist as it depends on the direction
        vertex_descriptor src_vertex = std::get<vertex_descriptor>(var_des);
        halfedge_descriptor href = halfedge(curr_tid, mesh);
        int src_index = source(href, mesh)==src_vertex?0:(target(href,mesh)==src_vertex?1:2);
        // first get the direction in the basis of the triangle to get its angle with an edge incident to the vertex

        Vector_2 v02 = curr_flat_tid[(src_index+2)%3] - curr_flat_tid[src_index];
        Vector_2 v01 = curr_flat_tid[(src_index+1)%3] - curr_flat_tid[src_index];
        double unnorm_cos_theta_1 = v02 * dir;
        double unnorm_cos_theta_2 = v02 * v01;

        // 2 orientation tests to check if the angle is larger or smaller than Pi
        // 2 orientation tests are needed as curr_flat_tid orientation could be inverted during flattening
        Orientation dir_ori=orientation(v02,dir), v01_ori=orientation(v02,v01);

        // check if dir is  the cone centered at curr_flat_tid[src_vertex] (TODO: should be a predicate!)
        if (orientation(dir, v01)==LEFT_TURN || orientation(dir, v02)==RIGHT_TURN)
        {
          //TODO: should be made robust with snapping and predicates or something à la robust construction traits

          // compute all the angles around the vertex from where the path starts
          std::vector<double> angles;
          double acc_angle=0;
          halfedge_descriptor hloop=href;
          while (target(hloop, mesh)!=src_vertex) hloop=next(hloop, mesh);
          halfedge_descriptor hstart=hloop;
          Vector_3 n(get(vpm,target(hloop, mesh)), get(vpm,source(hloop, mesh)));
          n /= std::sqrt(n*n);
          do{
            hloop=opposite(next(hloop, mesh), mesh);
            Vector_3 nv(get(vpm,target(hloop, mesh)), get(vpm,source(hloop, mesh)));
            nv /= std::sqrt(nv*nv);
            angles.push_back( std::acos(n * nv) );
            acc_angle+=angles.back();
            n = nv;
          }while(hloop!=hstart);

          CGAL_assertion( std::acos(unnorm_cos_theta_2 / std::sqrt((v01*v01) * (v02*v02))) == angles[0]); // will not be true because of rounding

          // normalize the angle to bring them in [0,1]
          for (double& angle : angles)
            angle /= acc_angle;

          // compute and normalize the target angle
          double target_theta = std::acos(unnorm_cos_theta_1 / std::sqrt((dir*dir) * (v02*v02)));
          if ( dir_ori!=v01_ori ) target_theta = 2 * CGAL_PI - target_theta;
          target_theta /= acc_angle;

          double curr_angle = 0;
          int ia = 0;
          do{
            curr_angle+=angles[ia];
            if (curr_angle >= target_theta) break;
            hloop=opposite(next(hloop, mesh), mesh);
            ++ia;
          }while(hloop!=hstart);

          CGAL_assertion(hloop!=hstart);

          // angle in target face wrt hloop
          double delta = (target_theta - curr_angle + angles[ia]) * acc_angle;

          // using the law of the sinus we have:
          CGAL_assertion(target(hloop, mesh) == src_vertex);
          Vector_3 vloop(get(vpm, source(hloop, mesh)),get(vpm, target(hloop, mesh)));
          Vector_3 vprev(get(vpm, source(hloop, mesh)),get(vpm, target(next(hloop, mesh), mesh)));
          double dvloop = std::sqrt(vloop*vloop);
          double dvprev = std::sqrt(vprev*vprev);
          double theta = std::acos( vloop*vprev / dvloop / dvprev );
          double d_opp_delta =  dvloop * sin(delta) / sin(CGAL_PI-delta-theta);

          double alpha = d_opp_delta/dvprev;

          if (is_border(hloop,mesh))
          {
            // border case falling outside of the mesh TODO: check with Claudio
            result.push_back(curr_p);
            return result;
          }

          halfedge_descriptor href = halfedge(face(hloop,mesh),mesh);
          curr_flat_tid=init_flat_triangle<K>(href,vpm,mesh);
          int src_id = source(href, mesh)==src_vertex?0:(target(href,mesh)==src_vertex?1:2);

#ifdef CGAL_DEBUG_BSURF
          std::cout << "new_face= " << face(href, mesh) << "\n";
          std::array<Vector_3, 3> debug_pt;
          debug_pt[0]=mesh.point(source(href,mesh)) - Point_3(0.,0.,0.);
          debug_pt[1]=mesh.point(target(href,mesh)) - Point_3(0.,0.,0.);
          debug_pt[2]=mesh.point(target(next(href,mesh),mesh)) - Point_3(0.,0.,0.);
          std::cout << "direction inter pt: " << (alpha) * debug_pt[(src_id+1)%3] + (1-alpha) * debug_pt[(src_id+2)%3] << "\n";
#endif
          // point on the opposite edge of src_vertex in face(href, mesh)
          Vector_2 ip = (alpha) * curr_flat_tid[(src_id+1)%3] + (1-alpha) * curr_flat_tid[(src_id+2)%3];
          dir = ip - curr_flat_tid[src_id];
          curr_tid=face(href, mesh);
          flat_p = curr_flat_tid[src_id];
          curr_p.second=CGAL::make_array(0.,0.,0.);
          curr_p.second[src_id]=1.;
          curr_p.first=curr_tid;
        }
      }
      break;
      case 1:
      {
        halfedge_descriptor h = std::get<halfedge_descriptor>(var_des);
        halfedge_descriptor href = halfedge(face(h,mesh), mesh);
        int opp_id = h==href?2:(next(href,mesh)==h?0:1);
        Vector_2 h_2d = curr_flat_tid[(opp_id+2)%3] - curr_flat_tid[(opp_id+1)%3];

        // check if the line starts in the current face (TODO should be a predicate)
        if ( orientation(h_2d, dir) == RIGHT_TURN )
        {
          if ( is_border(opposite(h, mesh), mesh) )
          {
            result.push_back(curr_p);
            return result;
          }

          // take the same point but in the opposite face
          halfedge_descriptor h_start=h;
          h=opposite(h, mesh);
          curr_tid=face(h, mesh);
          curr_p.first=curr_tid;
          href=halfedge(curr_tid, mesh);
          CGAL_assertion(start_loc.second[opp_id]==0);
          std::array<FT, 2> bary2 = CGAL::make_array(start_loc.second[(opp_id+2 )%3], start_loc.second[(opp_id+1)%3]); // swap done here
          opp_id = h==href?2:(next(href,mesh)==h?0:1);
          curr_p.second[opp_id]=FT(0);
          curr_p.second[(opp_id+1)%3]=bary2[0];
          curr_p.second[(opp_id+2)%3]=bary2[1];


          // dir must also be updated:
          // unfold the new triangle into the basis of the input one
          // compute the angle between the new triangle ref halfedge and dir
          // compute the new dir thanks to that angle.
          halfedge_descriptor start_href=halfedge(start_loc.first, mesh);
          int k=h_start==prev(start_href,mesh)?2 :( h_start==next(start_href,mesh)?1:0);
          std::array<Vector_2,3> new_flat_tid_in_curr_basis=unfold_face<K>(h_start,vpm, mesh,curr_flat_tid, k);

          Vector_2 ybase = new_flat_tid_in_curr_basis[1]-new_flat_tid_in_curr_basis[0];
          double cos_theta = ybase * dir / std::sqrt(ybase*ybase) / std::sqrt(dir*dir);
          double theta = std::acos(cos_theta);
          dir = Vector_2(std::sin(theta), std::cos(theta));

          curr_flat_tid=init_flat_triangle<K>(halfedge(curr_tid,mesh),vpm,mesh);
          flat_p= curr_p.second[0]*curr_flat_tid[0]+curr_p.second[1]*curr_flat_tid[1]+curr_p.second[2]*curr_flat_tid[2];
        }
      }
      break;
      default:
      break;
    }

    PMP::Face_location<TriangleMesh, FT> prev_p;
    Vector_2 curr_dir=dir;
    halfedge_descriptor h_ref=halfedge(curr_tid,mesh);
    halfedge_descriptor h_curr=h_ref;



    result.push_back(curr_p);
#ifdef CGAL_DEBUG_BSURF
  std::cout << "p= " << PMP::construct_point(curr_p,mesh) << ")\n";
#endif

    //TODO not sure why we need that
    auto [is_vert, kv] = point_is_vert<K>(curr_p);
    auto [is_edge, ke] = point_is_edge<K>(curr_p);
    if (is_vert)
      h_curr=get_halfedge(kv,h_ref);
    else if (is_edge)
      h_curr=get_halfedge(ke,h_ref);

#ifdef CGAL_DEBUG_BSURF
    std::cout << "Accumulated loop starts\n";
#endif
    while (accumulated < len)
    {
#ifdef CGAL_DEBUG_BSURF
      std::cout << "--->" << accumulated << " vs " << len << "\n";
#endif
      int curr_offset=get_halfedge_offset(h_ref,h_curr);

      auto [k, t1] = segment_in_tri<K>(flat_p, curr_flat_tid, curr_dir,curr_offset);
#ifdef CGAL_DEBUG_BSURF
      std::cout << "  t1 = " << t1 << std::endl;
#endif
      CGAL_assertion(k!=-1);
      std::array<FT,3> new_bary=make_array(0.,0.,0.);
      PMP::Face_location<TriangleMesh, FT> point_on_edge;

      new_bary[k] = 1 - t1;
      new_bary[(k + 1) % 3] = t1;
      point_on_edge.first=curr_tid;
      point_on_edge.second=new_bary;
      std::tie(is_vert, kv) = point_is_vert<K>(point_on_edge);
#ifdef CGAL_DEBUG_BSURF
      std::cout << "  is_vert? " << is_vert << "\n";
#endif
      if (is_vert)
      {
        point_on_edge.second = make_array(0.,0.,0.);
        point_on_edge.second[kv]=1;
        vertex_descriptor vid = target(get_halfedge(kv, prev(h_ref,mesh)), mesh);
#ifdef CGAL_DEBUG_BSURF
        std::cout<< "hit vertex "<< vid <<std::endl;
        std::cout<< "href "<< edge(h_ref,mesh) <<std::endl;
        std::cout<< "kv "<< kv <<std::endl;
#endif
        accumulated +=
            sqrt(squared_distance(PMP::construct_point(curr_p,mesh), get(vpm,vid)));

        // stop if hitting the boundary
        bool border_reached = false;
        for (halfedge_descriptor h : halfedges_around_target(vid, mesh))
          if (is_border(h, mesh))
          {
            //TODO: that's not necessarily a break! depends on dir
            result.push_back(point_on_edge);
            border_reached=true;
            prev_p=point_on_edge; //useful if accumulated > len TODO: what about below for edges?
            break;
          }
        if (border_reached) break;


  //   Point_3 vert = get(vpm,vid);
  //   Point_3 vert_adj=get(vpm,source(h,mesh));
  //   Vector_3 v = vert - vert_adj;

        //TODO add a 2D version of approximate_angle in CGAL
        // FT init_angle = approximate_angle(curr_flat_tid[(kv+2)%3]-curr_flat_tid[kv], curr_dir);
        auto tmp =curr_flat_tid[kv]-curr_flat_tid[(kv+2)%3];
        FT init_angle = approximate_angle(Vector_3(tmp.x(), tmp.y(), 0), Vector_3(curr_dir.x(), curr_dir.y(), 0));
        std::tie(curr_dir, curr_tid, h_curr) =
            polthier_condition_at_vert<K>(mesh,vpm,vid,curr_tid, init_angle);

        h_ref=halfedge(curr_tid,mesh);
        kv=get_vid_offset(h_ref,vid);
        CGAL_assertion(kv!=-1);
        std::array<FT,3> tmp_bary=make_array(0.,0.,0.);
        tmp_bary[kv]=1;

        curr_flat_tid=init_flat_triangle<K>(h_ref,vpm,mesh);
        prev_p = curr_p;

        curr_p=PMP::Face_location<TriangleMesh,FT>(curr_tid,tmp_bary);
        int k=get_vid_offset(h_ref,target(h_curr,mesh));
        flat_p=curr_flat_tid[k];
      }
      else
      {
        h_curr=opposite(get_halfedge(k, h_ref),mesh);

        // break if hitting the boundary
        if (is_border(h_curr, mesh))
        {
          result.push_back(point_on_edge);
          accumulated += sqrt(squared_distance(PMP::construct_point(point_on_edge,mesh),
                                               PMP::construct_point(curr_p,mesh)));
          prev_p=point_on_edge; //if accumulated > len we need to be in the common face
          break;
        }

        face_descriptor adj = face(h_curr,mesh);
        std::array<FT,2> curr_alpha=make_array(t1,1-t1); //reversed because will switch face
        new_bary=edge_barycentric_coordinate(h_curr,halfedge(adj,mesh),curr_alpha);
        prev_p = curr_p;
        curr_p.first=adj;
        curr_p.second= new_bary;
        accumulated += sqrt(squared_distance(PMP::construct_point(curr_p,mesh),
                                             PMP::construct_point(prev_p,mesh)));
        curr_tid = adj;
        h_ref=halfedge(curr_tid,mesh);
        //TODO curr_dir should be normalized every time (to try to avoid numerical errors)
        curr_dir = compute_new_dir<K>(h_ref,h_curr,prev_p.second,curr_p.second,vpm,mesh);
        curr_flat_tid=init_flat_triangle<K>(h_ref,vpm,mesh);
        flat_p= curr_p.second[0]*curr_flat_tid[0]+curr_p.second[1]*curr_flat_tid[1]+curr_p.second[2]*curr_flat_tid[2];
#ifdef CGAL_DEBUG_BSURF
        std::cout << "  h_curr " << edge(h_curr, mesh)<<"\n";
        std::cout << "  adj " << adj<<"\n";
        Vector_2 intersection_point=new_bary[0]*curr_flat_tid[0]+new_bary[1]*curr_flat_tid[1]+new_bary[2]*curr_flat_tid[2];
        std::cout << "  New intersection point is "<< intersection_point<<std::endl;
        std::cout << "  barycentric coordinates"<<new_bary[0]<<" "<<new_bary[1]<< " "<<new_bary[2]<<"\n";
        std::cout << "  curr_flat_tid"<<std::endl;
        std::cout << "  " << curr_flat_tid[0] << " " << curr_flat_tid[1] << " " << curr_flat_tid[2] << "\n";
        std::cout << "  New_Flat_p "<<flat_p.x()<<" "<< flat_p.y()<<" )\n";
        std::cout << "  New Dir"<<flat_p.x()<<" "<< flat_p.y()<<" "<<flat_p.x()+curr_dir.x()<<" "<<flat_p.y()+ curr_dir.y()<<"\n";
#endif
      }
      result.push_back(curr_p);
#ifdef CGAL_DEBUG_BSURF
      std::cout << "  inter pt: " << PMP::construct_point(curr_p, mesh) << "\n";
#endif
    }

    double excess = accumulated - len;

    if (excess <= 0)
    {
#ifdef CGAL_DEBUG_BSURF
      std::cout << "excess " << excess << std::endl;
#endif
      return result;
    }

    Point_3 prev_pos = PMP::construct_point(*std::next(result.rbegin()),mesh);
    Point_3 last_pos = PMP::construct_point(result.back(),mesh);
    double alpha = excess / sqrt((last_pos - prev_pos).squared_length());
    Point_3 pos = barycenter(prev_pos, alpha, last_pos, 1-alpha);
#ifdef CGAL_DEBUG_BSURF
    std::cout << "excess " << excess << std::endl;
    std::cout << "prev_pos " << prev_pos << std::endl;
    std::cout << "last_pos " << last_pos << std::endl;

    std::cout << "prev_p "<< prev_p.first << std::endl;
    std::cout << "pos " << pos << std::endl;

#endif

    auto [inside, bary] =
        point_in_triangle<K>(vpm,mesh,prev_p.first,pos); // TODO replace with function in PMP/locate.h
    if (!inside)
    {
      std::cout << "prev_pos " << prev_pos << "\n";
      std::cout << "last_pos " << last_pos << "\n";
      std::cout << "pos " << pos << "\n";
      std::cout << "error!This point should be in the triangle" << std::endl; //TODO this is a debug
    }

    result.pop_back();
    prev_p.second=bary;
    result.push_back(prev_p);

    return result;
  }

  // static
  // std::vector<mesh_point>
  // polthier_straightest_geodesic(const vector<vec3i> &triangles, const vector<vec3f> &positions,
  //                               const vector<vec3i> &adjacencies, const vector<vector<int>> &v2t,
  //                               const vector<vector<float>> &angles, const vector<float> &total_angles,
  //                               const mesh_point &p, const vec2f &dir, const float &len)
  // {
  //   auto result = vector<mesh_point>{};
  //   auto accumulated = 0.f;
  //   auto curr_tid = p.face;
  //   auto curr_p = p;
  //   auto prev_p = mesh_point{};
  //   auto curr_dir = dir;

  //   result.push_back(p);

  //   auto k_start = 0;
  //   auto [is_vert, kv] = point_is_vert(p);
  //   auto [is_edge, ke] = point_is_edge(p);
  //   if (is_vert)
  //     k_start = kv;
  //   else if (is_edge)
  //     k_start = ke;

  //   auto count = 0;
  //   while (accumulated < len)
  //   {
  //     ++count;
  //     auto [k, t1] = straightest_path_in_tri(triangles, positions, curr_p,
  //                                           curr_dir, k_start);

  //     auto new_bary = zero3f;
  //     auto point_on_edge = mesh_point{};
  //     if (k != -1) {
  //       new_bary[k] = 1 - t1;
  //       new_bary[(k + 1) % 3] = t1;
  //       point_on_edge = mesh_point{curr_tid, vec2f{new_bary.y, new_bary.z}};
  //     } else {
  //       std::tie(is_edge, ke) = point_is_edge(curr_p, 5e-3);
  //       std::tie(is_vert, kv) = point_is_vert(curr_p, 5e-3);
  //       auto bary = get_bary(curr_p.uv);
  //       if (is_edge) {
  //         k = ke;
  //         t1 = bary[(k + 1) % 3];
  //         point_on_edge = curr_p;
  //       } else if (is_vert) {
  //         auto bary3 = zero3f;
  //         bary3[kv] = 1;
  //         point_on_edge = {curr_p.face, {bary3.y, bary3.z}};
  //       } else {
  //         std::cout << "Error!This should not happen" << std::endl;
  //         return result;
  //       }
  //     }

  //     std::tie(is_vert, kv) = point_is_vert(point_on_edge);

  //     if (is_vert) {
  //       auto vid = triangles[curr_tid][kv];
  //       if (angles[vid].size() == 0)
  //         return result;

  //       accumulated +=
  //           length(eval_position(triangles, positions, curr_p) - positions[vid]);
  //       auto dir3d = normalize(eval_position(triangles, positions, curr_p) -
  //                             positions[vid]);
  //       prev_p = curr_p;
  //       std::tie(curr_dir, curr_p, k_start) =
  //           polthier_condition_at_vert(triangles, positions, adjacencies,
  //                                     total_angles, vid, curr_tid, dir3d);
  //       curr_tid = curr_p.face;
  //       if (curr_tid == -1)
  //         return result;

  //     } else {
  //       auto adj = adjacencies[curr_tid][k];
  //       if (adj == -1)
  //         return result;
  //       auto h = find(adjacencies[adj], curr_tid);

  //       new_bary = zero3f;
  //       new_bary[h] = t1;
  //       new_bary[(h + 1) % 3] = 1 - t1;

  //       prev_p = curr_p;
  //       curr_p = mesh_point{adj, vec2f{new_bary.y, new_bary.z}};
  //       accumulated += length(eval_position(triangles, positions, curr_p) -
  //                             eval_position(triangles, positions, prev_p));

  //       auto T = switch_reference_frame(triangles, positions, adj, curr_tid);
  //       curr_dir = switch_reference_frame_vector(T, curr_dir);

  //       curr_tid = adj;
  //       k_start = h;
  //     }

  //     result.push_back(curr_p);
  //   }

  //   auto excess = accumulated - len;
  //   auto prev_pos = eval_position(triangles, positions, result.rbegin()[1]);
  //   auto last_pos = eval_position(triangles, positions, result.back());
  //   auto alpha = excess / length(last_pos - prev_pos);
  //   auto pos = alpha * prev_pos + (1 - alpha) * last_pos;

  //   auto [inside, bary] =
  //       point_in_triangle(triangles, positions, prev_p.face, pos);
  //   if (!inside)
  //     std::cout << "error!This point should be in the triangle" << std::endl;

  //   result.pop_back();
  //   result.push_back(mesh_point{prev_p.face, bary});

  //   return result;
  // }

}; // end of Straightest_geodesic_imp
} //end of internal namespace


/*!
 * \ingroup VGSFunctions
 * computes a path on a triangle mesh that is computed by starting a walk on `tmesh`
 * given a direction and a maximum distance. The distance will not be achieved if a border edge
 * is reached before.
 * \tparam TriangleMesh a model of `FaceListGraph` and `EdgeListGraph`
 * \tparam FT floating point number type (float or double)
 * \param tmesh input triangle mesh to compute the path on
 * \param src the source of the path
 * \param len the distance to walk along the straightest
 * \param dir the initial direction of the walk, given as a 2D vector in the face of src, the halfedge of the face being the y-axis.
 * \return the straightest path (not containing `src`)
 * \todo add named parameters
 * \todo do we want to also have a way to return Bézier segments? The output is actually Bézier segments subdivided.
 * \todo offer something better than a 2D vector for the direction
 */
template <class K, class TriangleMesh>
std::vector<CGAL::Polygon_mesh_processing::Face_location<TriangleMesh, typename K::FT>>
straightest_geodesic(const CGAL::Polygon_mesh_processing::Face_location<TriangleMesh, typename K::FT> &src,
                     const typename K::Vector_2& dir,
                     const typename K::FT len,
                     const TriangleMesh &tmesh)
{
  //TODO replace with named parameter
  using VPM = typename boost::property_map<TriangleMesh, CGAL::vertex_point_t>::const_type;
  using Impl = internal::Straightest_geodesic_imp<K, TriangleMesh, VPM>;
  VPM vpm = get(CGAL::vertex_point, tmesh);


  return Impl::straightest_geodesic(src, tmesh, vpm, dir, len);
}

} } // end of CGAL::Vector_graphics_on_surfaces

#endif // CGAL_POLYGON_MESH_PROCESSING_BSURF_STRAIGHTEST_GEODESIC_H
