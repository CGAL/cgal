// Copyright (c) 2019 GeometryFactory (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Andreas Fabri

#ifndef CGAL_UNION_OF_POLYGONS_POYGONS_H
#define CGAL_UNION_OF_POLYGONS_POYGONS_H

#include <iostream>
#include <fstream>

#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Constrained_triangulation_face_base_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Projection_traits_xy_3.h>

#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Interval_skip_list.h>
#include <CGAL/Interval_skip_list_interval.h>

#include <vector>


namespace CGAL {

namespace union_of_polygons_impl
{
  template <class K>
  struct Local_construct_point_2
  {
    typename K::Point_2
    operator()(typename K::FT x, typename K::FT y) const
    {
      return typename K::Construct_point_2()(x,y);
    }
  };

  template <class K >
  struct Local_construct_point_2<Projection_traits_xy_3<K> >
  {
    typename K::Point_3
    operator()(typename K::FT x, typename K::FT y) const
    {
      return typename K::Construct_point_3()(x,y,0);
    }
  };
}

template <typename K>
struct Union_of_polygons
{

  struct V {
    int polylines;
    bool non_manifold;
    std::size_t id;
    V()
      : polylines(0), non_manifold(false), id(-1)
    {}
  };

  struct F {
    bool mark;
    int cc;
    F()
      :mark(false), cc(-1)
    {}
  };

  typedef typename K::FT                                              FT;
  typedef typename K::Point_2                                         Point_2;
  typedef typename K::Segment_2                                       Segment_2;
  typedef Exact_intersections_tag                      Itag;
  typedef Triangulation_vertex_base_with_info_2<V,K>   Vb;
  typedef Triangulation_face_base_with_info_2<F,K>     Fbb;
  typedef Constrained_triangulation_face_base_2<K,Fbb> Fb;
  typedef Triangulation_data_structure_2<Vb,Fb>        Tds;
  typedef Constrained_Delaunay_triangulation_2<K,Tds,Itag> CDT;
  typedef Constrained_triangulation_plus_2<CDT>          CDTP;
  typedef typename CDTP::Edge                                         Edge;
  typedef typename CDTP::Constraint_id                                Cid;
  typedef typename CDTP::Vertex_handle                                Vertex_handle;
  typedef typename CDTP::Face_handle                                  Face_handle;
  typedef typename CDTP::Finite_faces_iterator                        Finite_faces_iterator;
  typedef typename CDTP::Finite_vertices_iterator                     Finite_vertices_iterator;
  typedef typename CDTP::Vertices_in_constraint_iterator              Vertices_in_constraint_iterator;
  typedef typename CDTP::Edge_circulator                              Edge_circulator;
  typedef Polygon_with_holes_2<K>                                     Polygon_with_holes;
  typedef typename Polygon_with_holes::General_polygon_2            Polygon_2;
  typedef typename Polygon_with_holes::Hole_const_iterator          Hole_const_iterator;

  CDTP cdt;
  std::vector<Cid> cids;
  bool dilate = false;

  struct Interval_skip_list_interval_with_uv_and_sign
    : public Interval_skip_list_interval<FT>
  {
    typedef Interval_skip_list_interval<FT> Base;
    typedef typename Base::Value Value;
    Vertex_handle low, high;
    int sign;

    Interval_skip_list_interval_with_uv_and_sign(const Value& inf,
                                                 const Value& sup,
                                                 Vertex_handle low,
                                                 Vertex_handle high,
                                                 int sign,
                                                 bool lb,
                                                 bool hb)
      : Base(inf, sup,lb,hb)
      , low(low)
      , high(high)
      , sign(sign)
    {}

    bool operator==(const Interval_skip_list_interval_with_uv_and_sign& I) const
    {
      return ( (low == I.low) && (high == I.high) && (sign == I.sign) &&
               Base::operator==(I) );
    }

  };

  void mark_inside(Face_handle fh, int cc)
  {
    std::list<Face_handle> queue;
    queue.push_back(fh);
    fh->info().cc = cc;
    while(! queue.empty()){
      Face_handle fh = queue.front();
      queue.pop_front();
      for(int i = 0; i < 3; i++){
        typename CDT::Edge e(fh,i);
        if(! cdt.is_constrained(e)){
          Face_handle nh = fh->neighbor(i);
          if(nh->info().cc == -1){
            queue.push_back(nh);
            nh->info().cc = cc;
          }
        }
      }
    }
  }

public:

  Union_of_polygons(const std::list<Polygon_2>& polys, std::list<Polygon_2>&result)
  {
    for(const Polygon_2& p : polys){
      if(p.size() > 2){
        Cid cid = cdt.insert_constraint(p.vertices_begin(), p.vertices_end(), true);
        cids.push_back(cid);
      }
    }
    CGAL_assertion(cdt.is_valid(true));

    bool rerun = run();

    if(rerun){
      dilate_non_manifold_vertices();
      reset_markers();
      run();
    }

    write_result(result);
  }

  template <class SurfaceMesh>
  Union_of_polygons(const std::vector< std::vector<Point_2> >& polys, SurfaceMesh& result)
  {
    for(const std::vector<Point_2>& p : polys){
      if(p.size() > 2){
        Cid cid = cdt.insert_constraint(p.begin(), p.end(), true);
        cids.push_back(cid);
      }
    }
    CGAL_assertion(cdt.is_valid(true));

    bool rerun = run();

    if(rerun){
      dilate_non_manifold_vertices();
      reset_markers();
      run();
    }

    write_result_as_a_surface_mesh(result);
  }

  void reset_markers()
  {
    for(Finite_faces_iterator fh = cdt.finite_faces_begin(); fh != cdt.finite_faces_end(); fh++){
      fh->info() = F();
    }
    for(Finite_vertices_iterator vh = cdt.finite_vertices_begin(); vh != cdt.finite_vertices_end(); vh++){
      vh->info() = V();
    }
  }

  bool run()
  {
    std::vector<Face_handle> face_for_cc;
    typename K::Left_turn_2 left_turn;
    int cc = 0;
    for (Finite_faces_iterator fh = cdt.finite_faces_begin(); fh != cdt.finite_faces_end(); fh++) {
      if(fh->info().cc == -1){
        mark_inside(fh, cc);
        face_for_cc.push_back(fh);
        ++cc;
      }
    }

    std::vector<int> sums(cc);

    // put non-horizontal edges in ISL
    typedef Interval_skip_list_interval_with_uv_and_sign Interval;
    std::vector<Interval> intervals;
    // todo  interval.reserve();

    for(Cid cid : cids){
      for(Vertices_in_constraint_iterator it = cdt.vertices_in_constraint_begin(cid); ; ){
        Vertex_handle u = *it;
        ++it;
        if(it  == cdt.vertices_in_constraint_end(cid)){
          break;
        }
        Vertex_handle v = *it;
        if(u->point().y() == v->point().y()){
          continue;
        }
        if(u->point().y() < v->point().y()) { // upwards
          intervals.push_back(Interval(u->point().y(),v->point().y(), u, v, 1, true, false));
        } else {
          intervals.push_back(Interval(v->point().y(),u->point().y(), v, u, -1, true, false));
        }
      }
    }
    Interval_skip_list<Interval> isl(intervals.begin(),intervals.end());

    for(int i = 0; i < cc; i++){
      Face_handle fh = face_for_cc[i];
      Point_2 p = centroid(fh->vertex(0)->point(),
                                 fh->vertex(1)->point(),
                                 fh->vertex(2)->point());
      std::vector<Interval> res;

      int sum = 0;
      isl.find_intervals(p.y(), std::back_inserter(res));
      for(const Interval& in : res){
        if(left_turn(in.low->point(),in.high->point(), p)){
          sum += in.sign;
        }
      }
      sums[i] = sum;
    }

    for (Finite_faces_iterator fh = cdt.finite_faces_begin(); fh != cdt.finite_faces_end(); ++fh) {
      fh->info().mark = sums[fh->info().cc] > 0;
    }

    // For each vertex determine the number of polylines passing through
    for(Cid cid : cids){
      for(Vertices_in_constraint_iterator it = cdt.vertices_in_constraint_begin(cid);
          it  != cdt.vertices_in_constraint_end(cid); ++it){
        Vertex_handle u = *it;
        ++(u->info().polylines);
        if(u->info().polylines > 1){
          dilate = true;
        }
      }
    }
    return dilate;

  }

  Point_2 construct_point_2(FT x, FT y)
  {
    return
      typename union_of_polygons_impl::Local_construct_point_2<K>()(x,y);
  }

  void dilate_non_manifold_vertices()
  {
    std::map<Vertex_handle,double> non_manifold_vertices;

    // For each vertex determine if it is non-manifold and compute the distance to the one-ring
    for(Cid cid : cids){
      for(Vertices_in_constraint_iterator it = cdt.vertices_in_constraint_begin(cid);
          it  != cdt.vertices_in_constraint_end(cid); ++it){
        Vertex_handle u = *it;
        if(non_manifold_vertices.find(u) == non_manifold_vertices.end()){
          Point_2 up =u->point();
          if(u->info().polylines > 1){
            Edge_circulator ec = cdt.incident_edges(u), done(ec);
            int in_out = 0;
            double sdist = (std::numeric_limits<double>::max)();
            do {
              Edge e = *ec;
              Face_handle fh = e.first;
              int fi = e.second;
              CGAL_assertion(fh->vertex(CDTP::cw(fi)) == u);
              Vertex_handle v = fh->vertex(CDTP::ccw(fi));
              Vertex_handle w = fh->vertex(fi);
              if(!cdt.is_infinite(v)){
                double sqd = to_double(squared_distance(up, v->point()));
                sdist = (std::min)(sqd,sdist);
              }
              if( (!cdt.is_infinite(v)) &&  (!cdt.is_infinite(w))){
                Segment_2 seg(v->point(),w->point());
                double sqd = to_double(squared_distance(up, seg));
                sdist = (std::min)(sqd,sdist);
              }

              Face_handle nh = fh->neighbor(e.second);
              if(fh->info().mark != nh->info().mark){
                ++in_out;
              }
              ++ec;
            }while (ec != done);
            if(in_out > 2){
              if(! u->info().non_manifold){
                u->info().non_manifold = true;
                non_manifold_vertices[u] = std::sqrt(sdist)/10.0;
              }
            }
          }
        }
      }
    }

    for(const std::pair<Vertex_handle,double>& vh_dist : non_manifold_vertices){
      Polygon_2 poly;
      Point_2 p = vh_dist.first->point();
      FT d = vh_dist.second;
      poly.push_back(construct_point_2(p.x() - d,  p.y() - d));
      poly.push_back(construct_point_2(p.x() + d,  p.y() - d));
      poly.push_back(construct_point_2(p.x() + d,  p.y() + d));
      poly.push_back(construct_point_2(p.x() - d,  p.y() + d));
      Cid cid = cdt.insert_constraint(poly.vertices_begin(), poly.vertices_end(), true);
      cids.push_back(cid);
    }
  }

  void write_result(std::list<Polygon_2>& result)
  {
    // Collect all edges with a marked face to their left
    std::map<Vertex_handle, Vertex_handle> boundary_edges;

    for(Cid cid : cids){
      for(Vertices_in_constraint_iterator it = cdt.vertices_in_constraint_begin(cid); ; ){
        Vertex_handle u = *it;
        ++it;
        if(it  == cdt.vertices_in_constraint_end(cid)){
          break;
        }
        Vertex_handle v = *it;

        Face_handle fh;
        int fi;
        cdt.is_edge(u,v,fh,fi);
        Face_handle nh = fh->neighbor(fi);

        if(fh->info().mark != nh->info().mark){
          if(fh->info().mark){
            boundary_edges[v] = u;
          }else{
            boundary_edges[u] = v;
          }
        }
      }
    }

    // Extract polylines which are simple by construction
    while(! boundary_edges.empty()){
      typename std::map<Vertex_handle, Vertex_handle>::iterator it = boundary_edges.begin();
      Vertex_handle s = it->first;
      Polygon_2 poly;
      while((it = boundary_edges.find(s))!= boundary_edges.end()){
        s = it->second;
        poly.push_back(s->point());
        boundary_edges.erase(it);
      }
      result.push_back(poly);
    }
  }

  template <class SurfaceMesh>
  void write_result_as_a_surface_mesh(SurfaceMesh& sm)
  {
    std::vector<typename SurfaceMesh::Vertex_index> sm_verts;
    for (Finite_faces_iterator fh = cdt.finite_faces_begin(); fh != cdt.finite_faces_end(); ++fh)
    {
      if(fh->info().mark)
      {
        for (int i=0;i<3;++i)
        {
          if (fh->vertex(i)->info().id == std::size_t(-1))
          {
            fh->vertex(i)->info().id = sm_verts.size();
            sm_verts.push_back( sm.add_vertex( fh->vertex(i)->point() ) );
          }
        }
      }
    }

    for (typename CDT::Finite_faces_iterator fit=cdt.finite_faces_begin(), end=cdt.finite_faces_end(); fit!=end; ++fit)
    {
      if (fit->info().mark)
      {
        std::array<typename SurfaceMesh::Vertex_index, 3> fverts = { sm_verts[fit->vertex(0)->info().id],
                                                                     sm_verts[fit->vertex(1)->info().id],
                                                                     sm_verts[fit->vertex(2)->info().id] };
        sm.add_face(fverts);
      }
    }
  }
};// struct Union_of_polygons

} // namespace CGAL

#endif // CGAL_UNION_OF_POLYGONS_POYGONS_H
