// Copyright (c) 2005 Rijksuniversiteit Groningen (Netherlands)
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
//
//
// Author(s)     : Nico Kruithof <Nico@cs.rug.nl>

#ifndef CGAL_TRIANGULATE_POWER_DIAGRAM_3_H
#define CGAL_TRIANGULATE_POWER_DIAGRAM_3_H

#include <CGAL/license/Skin_surface_3.h>


#include <CGAL/Compute_anchor_3.h>

#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Triangulated_mixed_complex_observer_3.h>
#include <CGAL/Triangulation_incremental_builder_3.h>
// NGHK: move this one to SkinSurfaceTraits
#include <CGAL/Compute_anchor_3.h>

#include <CGAL/Union_find.h>

namespace CGAL {

template <
  class RegularTriangulation_3,
  class TriangulatedMixedComplex_3,
  class TriangulatedMixedComplexObserver_3 =
    Triangulated_mixed_complex_observer_3<TriangulatedMixedComplex_3,
                                          RegularTriangulation_3> >
class Power_diagram_triangulator_3
{
public:
  typedef typename RegularTriangulation_3::Geom_traits       Regular_traits;
  typedef typename TriangulatedMixedComplex_3::Geom_traits   Triangulated_mixed_complex_traits;

  typedef RegularTriangulation_3                     Regular;
  typedef TriangulatedMixedComplex_3                 Triangulated_mixed_complex;
  typedef TriangulatedMixedComplexObserver_3
  Triangulated_mixed_complex_observer;

private:
  typedef typename Regular::Vertex_handle            Rt_Vertex_handle;
  typedef typename Regular::Edge                     Rt_Edge;
  typedef typename Regular::Facet                    Rt_Facet;
  typedef typename Regular::Cell_handle              Rt_Cell_handle;

  typedef typename Regular::Finite_vertices_iterator Rt_Finite_vertices_iterator;
  typedef typename Regular::Finite_edges_iterator    Rt_Finite_edges_iterator;
  typedef typename Regular::Finite_facets_iterator   Rt_Finite_facets_iterator;
  typedef typename Regular::All_cells_iterator       Rt_All_cells_iterator;
  typedef typename Regular::Finite_cells_iterator    Rt_Finite_cells_iterator;

  typedef typename Regular::Cell_circulator          Rt_Cell_circulator;

  typedef Triangulation_simplex_3<Regular>           Rt_Simplex;
  typedef typename Regular::Bare_point               Rt_Bare_point;
  typedef typename Regular_traits::FT                Rt_FT;
  typedef typename Regular::Weighted_point           Rt_Weighted_point;

  typedef typename Triangulated_mixed_complex::Vertex_handle    Tmc_Vertex_handle;
  typedef typename Triangulated_mixed_complex::Edge             Tmc_Edge;
  typedef typename Triangulated_mixed_complex::Facet            Tmc_Facet;
  typedef typename Triangulated_mixed_complex::Cell_handle      Tmc_Cell_handle;

  typedef typename Triangulated_mixed_complex::Finite_vertices_iterator
  Tmc_Finite_vertices_iterator;
  typedef typename Triangulated_mixed_complex::Finite_edges_iterator
  Tmc_Finite_edges_iterator;
  typedef typename Triangulated_mixed_complex::Finite_facets_iterator
  Tmc_Finite_facets_iterator;
  typedef typename Triangulated_mixed_complex::All_cells_iterator
  Tmc_All_cells_iterator;
  typedef typename Triangulated_mixed_complex::Finite_cells_iterator
  Tmc_Finite_cells_iterator;
  typedef typename Triangulated_mixed_complex::Cell_circulator
  Tmc_Cell_circulator;

  typedef typename TriangulatedMixedComplex_3::Geom_traits  Tmc_traits;
  typedef typename Tmc_traits::Point_3                      Tmc_Point;
  typedef typename Tmc_traits::RT                           Tmc_RT;

  typedef Triangulation_incremental_builder_3<Triangulated_mixed_complex>
                                                          Triangulation_incremental_builder;

  typedef Compute_anchor_3<Regular>                       Compute_anchor;
  typedef std::pair<Rt_Simplex,Rt_Simplex>                Symb_anchor;

  // You might get type differences here:
//  struct Anchor_map_iterator_tmp;
//  typedef std::map<Rt_Simplex, Anchor_map_iterator_tmp>     Anchor_map;
//  struct Anchor_map_iterator_tmp : Anchor_map::iterator {
//    Anchor_map_iterator_tmp()
//      : Anchor_map::iterator() {}
//    Anchor_map_iterator_tmp(typename Anchor_map::iterator const &it)
//      : Anchor_map::iterator(it) {}
//  };
//  typedef typename Anchor_map::iterator                     Anchor_map_iterator;

  typedef Union_find<Rt_Simplex>                            Union_find_anchor;
  typedef std::map<Rt_Simplex,
                   typename Union_find_anchor::handle>      Simplex_UF_map;

public:
  Power_diagram_triangulator_3(Regular &regular,
                               Triangulated_mixed_complex &triangulated_mixed_complex,
                               Triangulated_mixed_complex_observer &observer,
                               bool verbose)
    : regular(regular),
      _tmc(triangulated_mixed_complex),
      observer(observer),
      triangulation_incr_builder(triangulated_mixed_complex),
      compute_anchor_obj(regular),
      verbose(verbose)
  {
    build();
  }

private:
  void build()
  {
    triangulation_incr_builder.begin_triangulation(3);

    if (verbose)
      std::cout << "Construct vertices" << std::endl;
    construct_vertices();

    if (verbose)
      std::cout << "Construct cells" << std::endl;
    construct_cells(); // mixed cells corresponding to regular vertices

    triangulation_incr_builder.end_triangulation();

    anchors.clear();

    CGAL_assertion(_tmc.is_valid());

    //remove_small_edges();

//     { // NGHK: debug code:
//       CGAL_assertion(_tmc.is_valid());
//       std::vector<Tmc_Vertex_handle> ch_vertices;
//       _tmc.incident_vertices(_tmc.infinite_vertex(),
//                              std::back_inserter(ch_vertices));
//       for (typename std::vector<Tmc_Vertex_handle>::iterator
//            vit = ch_vertices.begin(); vit != ch_vertices.end(); vit++) {
//         CGAL_assertion((*vit)->sign() == POSITIVE);
//       }
//     }
  }

  Tmc_Vertex_handle add_vertex(Rt_Simplex const &anchor);

  Tmc_Cell_handle add_cell(Tmc_Vertex_handle vh[], int orient, Rt_Simplex s);

  Tmc_Vertex_handle get_vertex(Rt_Simplex &sVor);

  void construct_anchor_vor(Rt_Simplex const &sVor);

  void construct_anchors();

  Rt_Simplex get_anchor_vor(Rt_Simplex const &sVor)
  {
    typename Simplex_UF_map::iterator it = anchor_vor_map.find(sVor);
    CGAL_assertion(it != anchor_vor_map.end());
    return *anchor_vor_uf.find(it->second);
  }

  void construct_vertices();

  Tmc_Point get_orthocenter(Rt_Simplex const &s);

  Tmc_Point get_anchor(Rt_Simplex const &sVor);

  template <class Point>
  Point construct_anchor_point(const Point &center_vor) {
    return center_vor;
  }

  void construct_cells();

  void remove_small_edges();

  bool is_collapsible(Tmc_Vertex_handle vh,
                      Tmc_Vertex_handle &vh_collapse_to,
                      Tmc_RT sq_length);
  void do_collapse(Tmc_Vertex_handle vh, Tmc_Vertex_handle vh_collapse_to);

private:
  Regular const &regular;
  Triangulated_mixed_complex &_tmc;
  Triangulated_mixed_complex_observer &observer;
  Triangulation_incremental_builder triangulation_incr_builder;
  typename Tmc_traits::Construct_weighted_circumcenter_3 orthocenter_obj;
  typename Tmc_traits::Compute_squared_radius_smallest_orthogonal_sphere_3 orthoweight_obj;
  Compute_anchor_3<Regular> compute_anchor_obj;
  bool verbose;

  Cartesian_converter<typename Rt_Bare_point::R,
                      Triangulated_mixed_complex_traits > r2t_converter_object;

  static const int edge_index[4][4];
  struct Index_c4 { Tmc_Vertex_handle V[4]; };
  struct Index_c6 { Tmc_Vertex_handle V[6]; };
  struct Index_c44 { Tmc_Vertex_handle V[4][4]; };
  struct Index_v {
    Unique_hash_map < Rt_Vertex_handle, Tmc_Vertex_handle > V;
  };

  // index to vertex
  Unique_hash_map < Rt_Cell_handle, Index_c4 > index_03;

  Union_find_anchor                            anchor_vor_uf;
  Simplex_UF_map                               anchor_vor_map;
//  Anchor_map                                   anchor_vor2;
  std::map<Rt_Simplex, Tmc_Vertex_handle>      anchors;
};

template <class RegularTriangulation_3,
          class TriangulatedMixedComplex_3,
          class TriangulatedMixedComplexObserver_3>
const int Power_diagram_triangulator_3<RegularTriangulation_3,
                                       TriangulatedMixedComplex_3,
                                       TriangulatedMixedComplexObserver_3>::
edge_index[4][4] = {{-1,0,1,2},{0,-1,3,4},{1,3,-1,5},{2,4,5,-1}};

template <class RegularTriangulation_3,
          class TriangulatedMixedComplex_3,
          class TriangulatedMixedComplexObserver_3>
void
Power_diagram_triangulator_3<RegularTriangulation_3,
                             TriangulatedMixedComplex_3,
                             TriangulatedMixedComplexObserver_3>::
construct_anchor_vor(Rt_Simplex const &sVor)
{
  Rt_Simplex s = compute_anchor_obj.anchor_vor(sVor);

  typename Union_find_anchor::handle sVor_handle, s_handle;
  sVor_handle = anchor_vor_uf.make_set(sVor);
  anchor_vor_map[sVor] = sVor_handle;

  typename Simplex_UF_map::iterator s_it = anchor_vor_map.find(s);
  CGAL_assertion(s_it != anchor_vor_map.end());

  anchor_vor_uf.unify_sets(sVor_handle, s_it->second);

  // degenerate simplices:
  if (compute_anchor_obj.is_degenerate()) {
    typename Compute_anchor::Simplex_iterator degenerate_it;
    typename Simplex_UF_map::iterator deg_map_it;
    for (degenerate_it = compute_anchor_obj.equivalent_anchors_begin();
         degenerate_it != compute_anchor_obj.equivalent_anchors_end();
         degenerate_it++) {
      deg_map_it = anchor_vor_map.find(*degenerate_it);

      // Possibly not found for 2 Voronoi vertices with the same center,
      // If the first vertex is inserted and the second is already found.
      // see compute_anchor_obj.anchor_vor(Cell_handle)
      if (deg_map_it != anchor_vor_map.end()) {
        anchor_vor_uf.unify_sets(sVor_handle, deg_map_it->second);
      }
    }
  }
}

template <class RegularTriangulation_3,
          class TriangulatedMixedComplex_3,
          class TriangulatedMixedComplexObserver_3>
void
Power_diagram_triangulator_3<RegularTriangulation_3,
                             TriangulatedMixedComplex_3,
                             TriangulatedMixedComplexObserver_3>::
construct_anchors()
{
  Rt_Finite_vertices_iterator vit;
  Rt_Finite_edges_iterator eit;
  Rt_Finite_facets_iterator fit;
  Rt_Finite_cells_iterator cit;
  Rt_Simplex s;

  // Compute anchor points:
  for (cit=regular.finite_cells_begin();
       cit!=regular.finite_cells_end(); cit++) {
    s = Rt_Simplex(cit);
    construct_anchor_vor(s);
    CGAL_assertion(s.dimension() == 3);
  }
  for (fit=regular.finite_facets_begin();
       fit!=regular.finite_facets_end(); fit++) {
    s = Rt_Simplex(*fit);
    construct_anchor_vor(s);
    CGAL_assertion(s.dimension() == 2);
  }
  for (eit=regular.finite_edges_begin();
       eit!=regular.finite_edges_end(); eit++) {
    s = Rt_Simplex(*eit);
    construct_anchor_vor(s);
    CGAL_assertion(s.dimension() == 1);
  }
  for (vit=regular.finite_vertices_begin();
       vit!=regular.finite_vertices_end(); vit++) {
    CGAL_assertion(vit->cell() != Rt_Cell_handle());
    s = Rt_Simplex(vit);
    construct_anchor_vor(s);
    CGAL_assertion(s.dimension() == 0);
  }
}

// Constructs the vertices of the simplicial complex
template <class RegularTriangulation_3,
          class TriangulatedMixedComplex_3,
          class TriangulatedMixedComplexObserver_3>
void
Power_diagram_triangulator_3<RegularTriangulation_3,
                             TriangulatedMixedComplex_3,
                             TriangulatedMixedComplexObserver_3>::
construct_vertices()
{
  Rt_All_cells_iterator acit;
  Rt_Finite_cells_iterator cit;
  Rt_Finite_facets_iterator fit;
  Rt_Finite_edges_iterator eit;
  Rt_Finite_vertices_iterator vit;
  Rt_Cell_circulator ccir, cstart;
  Rt_Vertex_handle v1, v2, v3;
  Rt_Edge e;
  Rt_Cell_handle c1, c2;
  Rt_Simplex sVor;
  Tmc_Vertex_handle vh;

  if (verbose) std::cout << "construct_anchors" << std::endl;
  construct_anchors();

  if (verbose) std::cout << "4 ";
  // anchor dimDel=0, dimVor=3
  for (cit=regular.finite_cells_begin();
       cit!=regular.finite_cells_end(); cit++) {
    sVor = get_anchor_vor(Rt_Simplex(cit));
    if (anchors.find(sVor) == anchors.end()) {
      vh = add_vertex(sVor);
      anchors[sVor] = vh;
      CGAL_assertion(vh == get_vertex(sVor));
    }
  }

  if (verbose) std::cout << "3 ";
  // anchor dimDel=2, dimVor=3 and dimDel=0, dimVor=2
  for (fit=regular.finite_facets_begin(); fit!=regular.finite_facets_end(); fit++) {
    // anchor dimDel=0, dimVor=2
    sVor = get_anchor_vor(*fit);
    if (anchors.find(sVor) == anchors.end()) {
      vh = add_vertex(sVor);
      anchors[sVor] = vh;
      CGAL_assertion(vh == get_vertex(sVor));
    }
  }

  if (verbose) std::cout << "2 ";
  // anchor dimDel=0, dimVor=1
  for (eit=regular.finite_edges_begin(); eit!=regular.finite_edges_end(); eit++) {
    sVor = get_anchor_vor(*eit);
    if (anchors.find(sVor) == anchors.end()) {
      vh = add_vertex(sVor);
      anchors[sVor] = vh;
      CGAL_assertion(vh == get_vertex(sVor));
    }
  }

  if (verbose) std::cout << "1 ";
  // anchor dimDel=0, dimVor=0
  for (vit=regular.finite_vertices_begin(); vit!=regular.finite_vertices_end(); vit++) {
    sVor = get_anchor_vor(Rt_Simplex(vit));
    if (anchors.find(sVor) == anchors.end()) {
      vh = add_vertex(sVor);
      anchors[sVor] = vh;
      CGAL_assertion(vh == get_vertex(sVor));
    }
  }
}

// Constructs the cells of the mixed complex corresponding
// to Regular vertices
template <class RegularTriangulation_3,
          class TriangulatedMixedComplex_3,
          class TriangulatedMixedComplexObserver_3>
void
Power_diagram_triangulator_3<RegularTriangulation_3,
                             TriangulatedMixedComplex_3,
                             TriangulatedMixedComplexObserver_3>::
construct_cells()
{
  Rt_Simplex sVor_v, sVor_e, sVor_f, sVor_c;
  Tmc_Vertex_handle vh[4];

  for (Rt_Finite_vertices_iterator vit=regular.finite_vertices_begin();
       vit!=regular.finite_vertices_end(); vit++) {

    Rt_Simplex simplex(vit);
    sVor_v = get_anchor_vor(Rt_Simplex(vit));
    vh[0] = get_vertex(sVor_v);

    std::list<Rt_Cell_handle> adj_cells;
    typename std::list<Rt_Cell_handle>::iterator adj_cell;
    regular.incident_cells(vit, std::back_inserter(adj_cells));

    // Construct cells:
    for (adj_cell = adj_cells.begin();
         adj_cell != adj_cells.end();
         adj_cell ++) {
      if (!regular.is_infinite(*adj_cell)) {
        sVor_c = get_anchor_vor(Rt_Simplex(*adj_cell));
        vh[3] = get_vertex(sVor_c);
        int index = (*adj_cell)->index(vit);
        for (int i=1; i<4; i++) {
          sVor_f = get_anchor_vor(
                     Rt_Simplex(Rt_Facet(*adj_cell,(index+i)&3)));
          vh[2] = get_vertex(sVor_f);

          for (int j=1; j<4; j++) {
            if (j!=i) {
              sVor_e = get_anchor_vor(
                         Rt_Simplex(Rt_Edge(*adj_cell,index,(index+j)&3)));
              vh[1] = get_vertex(sVor_e);
              if ((vh[0] != vh[1]) && (vh[1] != vh[2]) && (vh[2] != vh[3])) {
                CGAL_assertion(sVor_v != sVor_e);
                CGAL_assertion(sVor_e != sVor_f);
                CGAL_assertion(sVor_f != sVor_c);
                // Tmc_Cell_handle ch =
                add_cell(vh,(index + (j==(i%3+1)? 1:0))&1,simplex);
              }
            }
          }
        }
      }
    }
  }
}

// Adds a vertex to the simplicial complex
template <class RegularTriangulation_3,
          class TriangulatedMixedComplex_3,
          class TriangulatedMixedComplexObserver_3>
typename Power_diagram_triangulator_3<
  RegularTriangulation_3,
  TriangulatedMixedComplex_3,
  TriangulatedMixedComplexObserver_3>::Tmc_Vertex_handle
Power_diagram_triangulator_3<RegularTriangulation_3,
                             TriangulatedMixedComplex_3,
                             TriangulatedMixedComplexObserver_3>::
add_vertex(Rt_Simplex const &anchor)
{
  Tmc_Vertex_handle vh;
  vh = triangulation_incr_builder.add_vertex();
  vh->point() = get_anchor(anchor);
  observer.after_vertex_insertion(anchor, anchor, vh);

  return vh;
}

// Gets a vertex from the simplicial complex based on the anchors
template <class RegularTriangulation_3,
          class TriangulatedMixedComplex_3,
          class TriangulatedMixedComplexObserver_3>
typename Power_diagram_triangulator_3<
  RegularTriangulation_3,
  TriangulatedMixedComplex_3,
  TriangulatedMixedComplexObserver_3>::Tmc_Vertex_handle
Power_diagram_triangulator_3<RegularTriangulation_3,
                             TriangulatedMixedComplex_3,
                             TriangulatedMixedComplexObserver_3>::
get_vertex (Rt_Simplex &sVor)
{
  Rt_Simplex sVor2 = get_anchor_vor(sVor);
  CGAL_assertion(sVor == sVor2);
  Tmc_Vertex_handle vh = anchors[sVor2];
  CGAL_assertion(vh != Tmc_Vertex_handle());
  return vh;
}

// Adds a cell to the simplicial complex
template <class RegularTriangulation_3,
          class TriangulatedMixedComplex_3,
          class TriangulatedMixedComplexObserver_3>
typename Power_diagram_triangulator_3<
  RegularTriangulation_3,
  TriangulatedMixedComplex_3,
  TriangulatedMixedComplexObserver_3>::Tmc_Cell_handle
Power_diagram_triangulator_3<RegularTriangulation_3,
                             TriangulatedMixedComplex_3,
                             TriangulatedMixedComplexObserver_3>::
add_cell(Tmc_Vertex_handle vh[], int orient, Rt_Simplex s)
{
  CGAL_assertion((orient==0) || (orient==1));
  CGAL_assertion(vh[0] != Tmc_Vertex_handle()); CGAL_assertion(vh[1] != Tmc_Vertex_handle());
  CGAL_assertion(vh[2] != Tmc_Vertex_handle()); CGAL_assertion(vh[3] != Tmc_Vertex_handle());
  CGAL_assertion(vh[0] != vh[1]); CGAL_assertion(vh[0] != vh[2]); CGAL_assertion(vh[0] != vh[3]);
  CGAL_assertion(vh[1] != vh[2]); CGAL_assertion(vh[1] != vh[3]); CGAL_assertion(vh[2] != vh[3]);

  Tmc_Cell_handle ch;

  if (orient) {
    if (orientation(vh[0]->point(), vh[1]->point(),
                    vh[2]->point(), vh[3]->point()) != POSITIVE) {
      std::cout << orientation(vh[0]->point(), vh[1]->point(),
                               vh[2]->point(), vh[3]->point())<< std::endl;
    }
    CGAL_assertion(orientation(vh[0]->point(), vh[1]->point(),
                               vh[2]->point(), vh[3]->point()) == POSITIVE);
    ch = triangulation_incr_builder.add_cell(vh[0], vh[1], vh[2], vh[3]);
  } else {
    CGAL_assertion(orientation(vh[0]->point(), vh[1]->point(),
                               vh[3]->point(), vh[2]->point()) == POSITIVE);
    ch = triangulation_incr_builder.add_cell(vh[0], vh[1], vh[3], vh[2]);
  }
  observer.after_cell_insertion(s, ch);
  return ch;
}

template <class RegularTriangulation_3,
          class TriangulatedMixedComplex_3,
          class TriangulatedMixedComplexObserver_3>
typename TriangulatedMixedComplex_3::Geom_traits::Point_3
Power_diagram_triangulator_3<RegularTriangulation_3,
                             TriangulatedMixedComplex_3,
                             TriangulatedMixedComplexObserver_3>::
get_orthocenter(Rt_Simplex const &s)
{
  Rt_Vertex_handle vh;
  Rt_Edge           e;
  Rt_Facet          f;
  Rt_Cell_handle   ch;

  Tmc_Point result;
  switch (s.dimension()) {
  case 0:
    vh=s;
    result = Tmc_traits().construct_point_3_object()(r2t_converter_object(vh->point()));
    break;
  case 1:
    e=s;
    result = orthocenter_obj(
           r2t_converter_object(e.first->vertex(e.second)->point()),
           r2t_converter_object(e.first->vertex(e.third)->point()));
    break;
  case 2:
    f=s;
    result = orthocenter_obj(
           r2t_converter_object(
              f.first->vertex((f.second+1)&3)->point()),
           r2t_converter_object(
              f.first->vertex((f.second+2)&3)->point()),
           r2t_converter_object(
              f.first->vertex((f.second+3)&3)->point()));
    break;
  case 3:
    ch=s;
    result = orthocenter_obj(
           r2t_converter_object(ch->vertex(0)->point()),
           r2t_converter_object(ch->vertex(1)->point()),
           r2t_converter_object(ch->vertex(2)->point()),
           r2t_converter_object(ch->vertex(3)->point()));
    break;
  }
  return result;
}

template <class RegularTriangulation_3,
          class TriangulatedMixedComplex_3,
          class TriangulatedMixedComplexObserver_3>
typename TriangulatedMixedComplex_3::Geom_traits::Point_3
Power_diagram_triangulator_3<
  RegularTriangulation_3,
  TriangulatedMixedComplex_3,
  TriangulatedMixedComplexObserver_3>::
get_anchor(Rt_Simplex const &sVor)
{
  return get_orthocenter(sVor);
}

template <class RegularTriangulation_3,
          class TriangulatedMixedComplex_3,
          class TriangulatedMixedComplexObserver_3>
void
Power_diagram_triangulator_3<RegularTriangulation_3,
                             TriangulatedMixedComplex_3,
                             TriangulatedMixedComplexObserver_3>::
remove_small_edges()
{
  Bbox_3 bbox;
  for (Tmc_Finite_vertices_iterator vit = _tmc.finite_vertices_begin();
       vit != _tmc.finite_vertices_end(); vit++) {
    bbox = bbox+vit->point().bbox();
  }
  // Tmc_RT sq_length = ((bbox.xmax()-bbox.xmin())*(bbox.xmax()-bbox.xmin()) +
  //    (bbox.ymax()-bbox.ymin())*(bbox.ymax()-bbox.ymin()) +
  //    (bbox.zmax()-bbox.zmin())*(bbox.zmax()-bbox.zmin()))/100000000;

  Tmc_RT sq_length = 1e-6;
  // NGHK: This may intrudoce rounding errors, since the quadratic surface
  // may change:
  Tmc_Vertex_handle vh, vh_collapse_to;
  Tmc_Finite_vertices_iterator vit = _tmc.finite_vertices_begin();
  int nCollapsed=0;
  while (vit != _tmc.finite_vertices_end()) {
    vh = vit;
    vit++;
    if (is_collapsible(vh, vh_collapse_to,sq_length)) {
      nCollapsed ++;
      do_collapse(vh,vh_collapse_to);
    }
  }
  std::cout << "Collapsed: " << nCollapsed << std::endl;
}

template <class RegularTriangulation_3,
          class TriangulatedMixedComplex_3,
          class TriangulatedMixedComplexObserver_3>
bool
Power_diagram_triangulator_3<RegularTriangulation_3,
                             TriangulatedMixedComplex_3,
                             TriangulatedMixedComplexObserver_3>::
is_collapsible(Tmc_Vertex_handle vh,
               Tmc_Vertex_handle &vh_collapse_to,
               Tmc_RT sq_length)
{
  std::vector<Tmc_Cell_handle> incident_cells;
  CGAL_assertion(_tmc.is_vertex(vh));
  incident_cells.reserve(64);
  _tmc.incident_cells(vh, std::back_inserter(incident_cells));

  std::set<Tmc_Vertex_handle> incident_vertices;
  for(typename std::vector<Tmc_Cell_handle>::iterator
      cit = incident_cells.begin();
      cit != incident_cells.end(); ++cit) {
    // Put all incident vertices in incident_vertices.
    for (int j=0; j<4; ++j)
      if ((*cit)->vertex(j) != vh)
        incident_vertices.insert((*cit)->vertex(j));
  }

  for (typename std::set<Tmc_Vertex_handle>::iterator
       it = incident_vertices.begin();
       it != incident_vertices.end(); it++) {
    if ((_tmc.geom_traits().compute_squared_distance_3_object()(vh->point(),
                                                                (*it)->point())
         < sq_length) &&
        (vh->cell()->surf == (*it)->cell()->surf) &&
        (vh->sign() == (*it)->sign())) {
      bool ok = true;
      for (typename std::vector<Tmc_Cell_handle>::iterator
           cit = incident_cells.begin();
           ok && (cit != incident_cells.end()); cit++) {
        if (!(*cit)->has_vertex(*it)) {
          const Tmc_Point* pts[4] = { &((*cit)->vertex(0)->point()),
                                      &((*cit)->vertex(1)->point()),
                                      &((*cit)->vertex(2)->point()),
                                      &((*cit)->vertex(3)->point()) };
          pts[(*cit)->index(vh)] = &(*it)->point();

          ok = (_tmc.geom_traits().orientation_3_object()
                (*pts[0],*pts[1],*pts[2],*pts[3]) == CGAL::POSITIVE);
        }
      }
      if (ok) {
        vh_collapse_to = *it;
        return true;
      }
    }
  }
  return false;
}

template <class RegularTriangulation_3,
          class TriangulatedMixedComplex_3,
          class TriangulatedMixedComplexObserver_3>
void
Power_diagram_triangulator_3<RegularTriangulation_3,
TriangulatedMixedComplex_3,
TriangulatedMixedComplexObserver_3>::
do_collapse(Tmc_Vertex_handle vh, Tmc_Vertex_handle vh_collapse_to)
{
  std::vector<Tmc_Cell_handle> incident_cells;
  incident_cells.reserve(32);
  _tmc.incident_cells(vh, std::back_inserter(incident_cells));
  int i,i2;
  for (typename std::vector<Tmc_Cell_handle>::iterator
       it = incident_cells.begin(); it != incident_cells.end(); it++) {
    i = (*it)->index(vh);
    if ((*it)->has_vertex(vh_collapse_to,i2)) {
      // This cell is collapsed, set neighbor information of the new facet
      // and set the cell-pointer of the incident vertices.
      Tmc_Cell_handle ch1 = (*it)->neighbor(i);
      Tmc_Cell_handle ch2 = (*it)->neighbor(i2);
      ch1->set_neighbor(ch1->index((*it)), ch2);
      ch2->set_neighbor(ch2->index((*it)), ch1);
      for (int i=0; i<4; i++) {
        // Try to point to a cell with the same surface:
        if ((*it)->vertex(i)->cell() == (*it)) {
          if ((*it)->surf == ch1->surf) {
            (*it)->vertex(i)->set_cell(ch1);
          } else {
            (*it)->vertex(i)->set_cell(ch2);
          }
        }
      }
      _tmc.tds().delete_cell((*it));
    } else {
      // This cell is changed, set pointer to the new vertex
      (*it)->set_vertex(i,vh_collapse_to);
    }
  }
  _tmc.tds().delete_vertex(vh);
}

template <class RegularTriangulation_3,
          class TriangulatedMixedComplex_3,
          class TriangulatedMixedComplexObserver_3>
void
triangulate_power_diagram_3(RegularTriangulation_3 &rt,
                            TriangulatedMixedComplex_3 &tmc,
                            TriangulatedMixedComplexObserver_3 &observer,
                            bool verbose)
{
  typedef Power_diagram_triangulator_3<
    RegularTriangulation_3,
    TriangulatedMixedComplex_3,
    TriangulatedMixedComplexObserver_3>  Power_diagram_triangulator;
  Power_diagram_triangulator(rt, tmc, observer, verbose);
}

template <class RegularTriangulation_3,
          class TriangulatedMixedComplex_3>
void
triangulate_power_diagram_3(RegularTriangulation_3 const &regular,
                            TriangulatedMixedComplex_3 &tmc,
                            bool verbose)
{
  Triangulated_mixed_complex_observer_3<
    TriangulatedMixedComplex_3, const RegularTriangulation_3> observer(1);

  triangulate_power_diagram_3(regular, tmc, observer, verbose);
}

} //namespace CGAL

#endif // CGAL_TRIANGULATE_POWER_DIAGRAM_3_H
