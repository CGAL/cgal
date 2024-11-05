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

#ifndef CGAL_POLYGON_REPAIR_BOOLEAN_H
#define CGAL_POLYGON_REPAIR_BOOLEAN_H

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

/*!
\ingroup PkgTriangulation2Miscellaneous

\tparam Kernel must be
*/

template <typename Kernel>
class Boolean {

private:
  struct FaceInfo {

    FaceInfo()
    {}

    int label;
    int layers;
    bool processed;

    template <typename Fct>
    bool
    in_domain(const Fct& fct) const
    {
      return fct(layers);
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


  // @todo taken from Polygon_repair should be factorized
  struct Polygon_less {

    bool
    operator()(const Polygon_2& pa, const Polygon_2& pb) const
    {
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


  // @todo taken from Polygon_repair should be factorized
  struct Polygon_with_holes_less {
    Polygon_less pl;

    bool
    operator()(const Polygon_with_holes_2& pa, const Polygon_with_holes_2& pb) const
    {
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


  CDTplus cdt;

public:

/*!
default constructor.
*/
  Boolean() = default;


/*!
sets the polygons as input of the %Boolean operation.
*/
  void
  insert(const Multipolygon_with_holes_2& pA)
  {
    for(const auto& pwh : pA.polygons_with_holes()){
      cdt.insert_constraint(pwh.outer_boundary().vertices_begin(), pwh.outer_boundary().vertices_end(), true);
      for(auto const& hole : pwh.holes()){
        cdt.insert_constraint(hole.vertices_begin(), hole.vertices_end(), true);
      }
    }

    mark_domains();
  }

private:

  void
  mark_domains(Face_handle start,
               int index,
               std::list<Edge>& border)
  {
    if(start->info().layers != -1){
      return;
    }
    std::list<Face_handle> queue;
    queue.push_back(start);

    while(! queue.empty()){
      Face_handle fh = queue.front();
      queue.pop_front();
      if(fh->info().layers == -1){
        fh->info().layers = index;
        for(int i = 0; i < 3; i++){
          Edge e(fh,i);
          Face_handle n = fh->neighbor(i);
          if(n->info().layers == -1){
            if(cdt.is_constrained(e)){
              border.push_back(e);
            }else{
              queue.push_back(n);
            }
          }
        }
      }
    }
  }


  // this marks how many multipolygon interiors overlap a cell of the arrangement of mutipolygons
  void
  mark_domains()
  {
    for(Face_handle f : cdt.all_face_handles()){
      f->info().layers = -1;
    }

    int overlays = 0;
    std::list<Edge> border;
    mark_domains(cdt.infinite_face(), overlays, border);

    while(! border.empty()){
      Edge e = border.front();
      border.pop_front();
      Face_handle fh = e.first;
      int  fi = e.second;
      Face_handle n = fh->neighbor(fi);
      if(n->info().layers == -1){
        Vertex_handle u = fh->vertex(cdt.cw(fi)), v = fh->vertex(cdt.ccw(fi));
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
        mark_domains(n, fh->info().layers+delta, border);
      }
    }
  }


  template <typename Fct>
  void
  label_domain(Face_handle start, int label, const Fct& fct)
  {
    std::list<Face_handle> queue;
    start->info().label = label;
    queue.push_back(start);

    while(! queue.empty()){
      Face_handle fh = queue.front();
      queue.pop_front();

      for(int i = 0; i < 3; i++){
        Face_handle n = fh->neighbor(i);
        if(n->info().in_domain(fct)){
          if(n->info().label == 0){
            n->info().label = label;
            queue.push_back(n);
          }
        }
      }
    }
  }

  // this marks the domain for the Boolean operation applied on the two multipolygons
  template <typename Fct>
  int
  label_domains(const Fct& fct)
  {
    for (auto const face: cdt.all_face_handles()) {
      face->info().processed = false;
      face->info().label = 0;
    }
    int label = 1;
    for (auto const face: cdt.all_face_handles()) {
      if(face->info().in_domain(fct) && face->info().label == 0){
        label_domain(face, label, fct);
        ++label;
      }
    }
    return label;
  }



  // @todo taken from Polygon_repair and adapted; might be factorized
  // Reconstruct ring boundary starting from an edge (face + opposite vertex) that is part of it
  template <typename Fct>
  void
  reconstruct_ring(std::list<Point_2>& ring,
                   Face_handle face_adjacent_to_boundary,
                   int opposite_vertex,
                   const Fct& fct)
  {
    // Create ring
    Face_handle current_face = face_adjacent_to_boundary;
    int current_opposite_vertex = opposite_vertex;
    CGAL_assertion(face_adjacent_to_boundary->info().in_domain(fct));
    do {
      CGAL_assertion(current_face->info().in_domain(fct) == face_adjacent_to_boundary->info().in_domain(fct));
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


public:

  // @todo taken from Polygon_repair and adapted; might be factorized

/*!
performs the Boolean operation applying `fct` and returns the result as a multipolygon with holes.

\tparam Fct must have the operator `bool operator()(bool, bool)`.
*/
  template <typename Fct>
  Multipolygon_with_holes_2
  operator()(const Fct& fct)
  {
    int number_of_polygons = label_domains(fct) - 1;
    Multipolygon_with_holes_2 mp;
    std::vector<Polygon_2> polygons; // outer boundaries
    std::vector<std::set<Polygon_2, Polygon_less>> holes; // holes are ordered (per polygon)
    polygons.resize(number_of_polygons);
    holes.resize(number_of_polygons);

    for (auto const face: cdt.all_face_handles()) {
      face->info().processed = false;
    }

    for (auto const &face: cdt.finite_face_handles()) {
      if (! face->info().in_domain(fct)) continue; // exterior triangle
      if (face->info().processed) continue; // already reconstructed
      for (int opposite_vertex = 0; opposite_vertex < 3; ++opposite_vertex) {
        if (face->info().in_domain(fct) == face->neighbor(opposite_vertex)->info().in_domain(fct)) continue; // not adjacent to boundary

        // Reconstruct ring
        std::list<Point_2> ring;
        reconstruct_ring(ring, face, opposite_vertex, fct);

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


/*!
access to the underlying constrained triangulation.
*/
  const CDTplus&
  triangulation() const
  {
    return cdt;
  }

};

template <typename K>
Multipolygon_with_holes_2<K>
join(const Multipolygon_with_holes_2<K>& pA)
{
  CGAL::Polygon_repair::Boolean<K> bops;
  bops.insert(pA);
  struct Larger_than_zero {
    bool operator()(int i) const
    {
      return i > 0;
    }
  };
  Larger_than_zero ltz;
  return bops(ltz);
}



template <typename K>
Multipolygon_with_holes_2<K>
intersection(const Multipolygon_with_holes_2<K>& pA)
{
  CGAL::Polygon_repair::Boolean<K> bops;
  bops.insert(pA);
  struct Equal  {
    int val;
    Equal(int val)
    : val(val)
    {}

    bool operator()(int i) const
    {
      return i == val;
    }
  };
  Equal equal(pA.number_of_polygons_with_holes());
  return bops(equal);
}

} // namespace Polygon_repair
} //namespace CGAL

#endif CGAL_POLYGON_REPAIR_BOOLEAN_H
