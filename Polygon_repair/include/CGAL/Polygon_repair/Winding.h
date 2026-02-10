// Copyright (c) 2024  GeometryFactory SARL (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Andreas Fabri

#ifndef CGAL_POLYGON_REPAIR_WINDING_H
#define CGAL_POLYGON_REPAIR_WINDING_H

#include <CGAL/license/Polygon_repair.h>

#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Multipolygon_with_holes_2.h>

#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>

#include <vector>
#include <iostream>

#include <boost/property_map/property_map.hpp>

namespace CGAL {
namespace Polygon_repair {

#ifndef DOXYGEN_RUNNING
/*
\ingroup PkgPolygonRepairFunctions

\tparam Kernel must be
*/

template <typename Kernel>
class Winding {

private:
  struct FaceInfo {

    FaceInfo()
    {}

    int wind;
    int label;
    int nesting_level[2];
    bool processed;

    bool
    in_domain(int i) const
    {
      return nesting_level[i] % 2 == 1;
    }

    template <typename Fct>
    bool
    in_domain(const Fct& fct) const
    {
      return fct(in_domain(0), in_domain(1));
    }

  };

  using K = Kernel;
  using Point_2 = typename K::Point_2;
  using Polygon_2 = CGAL::Polygon_2<K>;
  using Polygon_with_holes_2 = CGAL::Polygon_with_holes_2<K>;
  using Multipolygon_with_holes_2 = CGAL::Multipolygon_with_holes_2<K>;

  using Itag = std::conditional_t<std::is_floating_point_v<typename K::FT>, Exact_predicates_tag, Exact_intersections_tag>;
  using Vb = CGAL::Triangulation_vertex_base_2<K>;
  using Fbb = CGAL::Triangulation_face_base_with_info_2<FaceInfo,K>;
  using Fb = CGAL::Constrained_triangulation_face_base_2<K,Fbb>;
  using Tds = CGAL::Triangulation_data_structure_2<Vb,Fb>;
  using CDT = CGAL::Constrained_Delaunay_triangulation_2<K,Tds,Itag>;
  using CDTplus = CGAL::Constrained_triangulation_plus_2<CDT>;
  using Face_handle = typename CDTplus::Face_handle;
  using Face_circulator = typename CDTplus::Face_circulator;
  using Vertex_handle = typename CDTplus::Vertex_handle;
  using Constraint_id = typename CDTplus::Constraint_id;
  using Edge = typename CDTplus::Edge;
  using Context = typename CDTplus::Context;

  CDTplus cdt;
  constexpr static int uninitialized = std::numeric_limits<int>::lowest();


  struct Polygon_less {
    bool operator()(const Polygon_2& pa, const Polygon_2& pb) const {
      typename Polygon_2::Vertex_iterator va = pa.vertices_begin();
      typename Polygon_2::Vertex_iterator vb = pb.vertices_begin();
      while (va != pa.vertices_end() && vb != pb.vertices_end()) {
        if (*va != *vb) return *va < *vb;
        ++va;
        ++vb;
      }
      if (vb == pb.vertices_end()) return false;
      return true;
    }
  };

  struct Polygon_with_holes_less {
    Polygon_less pl;
    bool operator()(const Polygon_with_holes_2& pa, const Polygon_with_holes_2& pb) const {
      if (pl(pa.outer_boundary(), pb.outer_boundary())) return true;
      if (pl(pb.outer_boundary(), pa.outer_boundary())) return false;
      typename Polygon_with_holes_2::Hole_const_iterator ha = pa.holes_begin();
      typename Polygon_with_holes_2::Hole_const_iterator hb = pb.holes_begin();
      while (ha != pa.holes_end() && hb != pb.holes_end()) {
        if (pl(*ha, *hb)) return true;
        if (pl(*hb, *ha)) return false;
      }
      if (hb == pb.holes_end()) return false;
      return true;
    }
  };

public:

/*!
default constructor.
*/
  Winding() = default;


/*!
sets the polygon as input of the winding number computation.
*/
  void
  insert(const Polygon_2& p)
  {
      cdt.insert_constraint(p.vertices_begin(), p.vertices_end(), true);
  }

/*!
sets the polygon as input of the winding number computation.
*/
  void
  insert(const Polygon_with_holes_2& pwh)
  {
      cdt.insert_constraint(pwh.outer_boundary().vertices_begin(), pwh.outer_boundary().vertices_end(), true);
      for(auto const& hole : pwh.holes()){
        cdt.insert_constraint(hole.vertices_begin(), hole.vertices_end(), true);
      }
  }

/*!
sets the polygon as input of the winding number computation.
*/
  void
  insert(const Multipolygon_with_holes_2& mpwh)
  {
    for(const auto& pwh : mpwh.polygons_with_holes()){
      cdt.insert_constraint(pwh.outer_boundary().vertices_begin(), pwh.outer_boundary().vertices_end(), true);
      for(auto const& hole : pwh.holes()){
        cdt.insert_constraint(hole.vertices_begin(), hole.vertices_end(), true);
      }
    }
  }

  void label(Face_handle f, int index, std::list<std::pair<Edge,int>>& border)
  {
    std::list<Face_handle> queue;
    queue.push_back(f);
    while(! queue.empty()){
      Face_handle fh = queue.front();
      queue.pop_front();
      fh->info().wind = index;
      for(int i = 0; i < 3; i++){
        Edge e(fh,i);
        Face_handle n = fh->neighbor(i);
        if(! cdt.is_constrained(e)){
          if(n->info().wind != index){
            queue.push_back(n);
          }
        }else{
          if(n->info().wind == uninitialized){
            Vertex_handle u = e.first->vertex(cdt.cw(e.second));
            Vertex_handle v = e.first->vertex(cdt.ccw(e.second));
            int delta = 0;
            for(Context c : cdt.contexts(u,v)){
              if(*c.current() ==u && *std::next(c.current()) == v){
                ++delta;
              }else if(*c.current() ==v && *std::next(c.current()) == u){
                --delta;
              }else{
                CGAL_assertion(false);
              }
            }
            border.push_back(std::make_pair(Edge(n,n->index(fh)),index+delta));
          }
        }
     }
   }
  }

  void label()
  {
    std::list<std::pair<Edge,int>> border;
    for(Face_handle f : cdt.all_face_handles()){
      f->info().wind = uninitialized;
    }
    int index = 0;
    label(cdt.infinite_face(),index++, border);
    while(! border.empty()){
      Edge e;
      int wind;
      std::tie(e,wind) = border.front();
      border.pop_front();
      if(e.first->info().wind == uninitialized){
        label(e.first, wind,border);
      }else{
        CGAL_assertion(e.first->info().wind == wind);
      }
    }
  }

  void
  label_domain(Face_handle start, int label)
  {
    std::list<Face_handle> queue;
    queue.push_back(start);

    while(! queue.empty()){
      Face_handle fh = queue.front();
      queue.pop_front();
      fh->info().label = label;

      for(int i = 0; i < 3; i++){
        Face_handle n = fh->neighbor(i);
        if(n->info().wind != 0 && n->info().label == 0){
            queue.push_back(n);
          }
      }
    }
  }

  int
  label_domains()
  {
    for (auto const face: cdt.all_face_handles()) {
      face->info().processed = false;
      face->info().label = 0;
    }
    int label = 1;
    for (auto const face: cdt.all_face_handles()) {
      if(face->info().wind != 0 && face->info().label == 0){
        label_domain(face, label);
        ++label;
      }
    }
    return label;
  }

  void
  reconstruct_ring(std::list<Point_2>& ring,
                   Face_handle face_adjacent_to_boundary,
                   int opposite_vertex)
  {
    // Create ring
    Face_handle current_face = face_adjacent_to_boundary;
    int current_opposite_vertex = opposite_vertex;
    CGAL_assertion(face_adjacent_to_boundary->info().wind != 0);
    do {
      current_face->info().processed = true;
      Vertex_handle pivot_vertex = current_face->vertex(current_face->cw(current_opposite_vertex));
      // std::cout << "\tAdding point " << pivot_vertex->point() << std::endl;
      ring.push_back(pivot_vertex->point());
      Face_circulator fc = cdt.incident_faces(pivot_vertex, current_face);
      do {
        ++fc;
      } while (fc->info().label != current_face->info().label);
      current_face = fc;
      current_opposite_vertex = fc->cw(fc->index(pivot_vertex));
    } while (current_face != face_adjacent_to_boundary ||
             current_opposite_vertex != opposite_vertex);

    // Start at lexicographically smallest vertex
    typename std::list<Point_2>::iterator smallest_vertex = ring.begin();
    for (typename std::list<Point_2>::iterator current_vertex = ring.begin();
         current_vertex != ring.end(); ++current_vertex) {
      if (*current_vertex < *smallest_vertex) smallest_vertex = current_vertex;
    }
    if (ring.front() != *smallest_vertex) {
      ring.splice(ring.begin(), ring, smallest_vertex, ring.end());
    }
  }

  Multipolygon_with_holes_2
  operator()()
  {
    int number_of_polygons = label_domains() - 1;

    Multipolygon_with_holes_2 mp;
    std::vector<Polygon_2> polygons; // outer boundaries
    std::vector<std::set<Polygon_2, Polygon_less>> holes; // holes are ordered (per polygon)
    polygons.resize(number_of_polygons);
    holes.resize(number_of_polygons);

    for (auto const face: cdt.all_face_handles()) {
      face->info().processed = false;
    }

    for (auto const &face: cdt.finite_face_handles()) {
      if (face->info().wind==0) continue; // exterior triangle
      if (face->info().processed) continue; // already reconstructed
      for (int opposite_vertex = 0; opposite_vertex < 3; ++opposite_vertex) {

        if ((face->info().wind != 0) == (face->neighbor(opposite_vertex)->info().wind != 0)) continue; // not adjacent to boundary

        // Reconstruct ring
        std::list<Point_2> ring;
        reconstruct_ring(ring, face, opposite_vertex);

        // Put ring in polygons
        Polygon_2 polygon;
        polygon.reserve(ring.size());
        polygon.insert(polygon.vertices_end(), ring.begin(), ring.end());
        if (polygon.orientation() == CGAL::COUNTERCLOCKWISE) {
          polygons[face->info().label-1] = std::move(polygon);
        } else {
          holes[face->info().label-1].insert(std::move(polygon));
        } break;
      }
    }

    // Create polygons with holes and put in multipolygon
    std::set<Polygon_with_holes_2, Polygon_with_holes_less> ordered_polygons;
    for (std::size_t i = 0; i < polygons.size(); ++i) {
      ordered_polygons.insert(Polygon_with_holes_2(std::move(polygons[i]),
                                                   std::make_move_iterator(holes[i].begin()),
                                                   std::make_move_iterator(holes[i].end())));
    }
    for (auto const& polygon: ordered_polygons) {
      mp.add_polygon_with_holes(std::move(polygon));
    }
    return mp;
  }

};

#endif

} // namespace Polygon_repair

} // namespace CGAL

#endif // CGAL_POLYGON_REPAIR_WINDING_H
