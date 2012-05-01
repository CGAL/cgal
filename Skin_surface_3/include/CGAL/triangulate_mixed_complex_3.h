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

#ifndef CGAL_TRIANGULATE_MIXED_COMPLEX_3
#define CGAL_TRIANGULATE_MIXED_COMPLEX_3

// #include <CGAL/Unique_hash_map.h>
#include <CGAL/Compute_anchor_3.h>

#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Triangulated_mixed_complex_observer_3.h>
#include <CGAL/Triangulation_incremental_builder_3.h>

#include <CGAL/Skin_surface_base_3.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

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
class Mixed_complex_triangulator_3 {
public:
  typedef typename RegularTriangulation_3::Geom_traits
  Regular_traits;
  typedef typename TriangulatedMixedComplex_3::Geom_traits
  Triangulated_mixed_complex_traits;

  typedef RegularTriangulation_3                   Regular;
  typedef TriangulatedMixedComplex_3               Triangulated_mixed_complex;
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
  typedef typename Regular::Bare_point               Rt_Point;
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
  // The map that maps a Rt_Simplex to an iterator of the map 
  // (used as union_find_structure)
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
                   typename Union_find_anchor::handle> Simplex_UF_map;


public:
  Mixed_complex_triangulator_3(Regular &regular,
			       Rt_FT const &shrink,
			       Triangulated_mixed_complex
			       &triangulated_mixed_complex,
			       Triangulated_mixed_complex_observer &observer,
			       bool verbose)
    : regular(regular),
      shrink(shrink),
      _tmc(triangulated_mixed_complex),
      observer(observer),
      triangulation_incr_builder(triangulated_mixed_complex), 
      construct_anchor_point_3_obj(r2t_converter_object(shrink)),
      compute_anchor_obj(regular),
      verbose(verbose)  {
    
    build();
  }

private:
  void build() {

    triangulation_incr_builder.begin_triangulation(3);

    if (verbose) std::cout << "Construct vertices" << std::endl;
    construct_vertices();

    // mixed cells corresponding to regular vertices
    if (verbose) std::cout << "Construct 0 cells" << std::endl;
    for (Rt_Finite_vertices_iterator vit = regular.finite_vertices_begin();
	 vit != regular.finite_vertices_end(); vit ++) {
      construct_0_cell(vit);
    }

    // mixed cells corresponding to regular edges
    if (verbose) std::cout << "Construct 1 cells" << std::endl;
    for (Rt_Finite_edges_iterator eit = regular.finite_edges_begin();
	 eit != regular.finite_edges_end(); eit ++) {
      construct_1_cell(eit);
    }

    // mixed cells corresponding to regular facets
    if (verbose) std::cout << "Construct 2 cells" << std::endl;
    for (Rt_Finite_facets_iterator fit = regular.finite_facets_begin();
	 fit != regular.finite_facets_end(); fit ++) {
      construct_2_cell(fit);
    }
    
    // mixed cells corresponding to regular cells
    if (verbose) std::cout << "Construct 3 cells" << std::endl;
    for (Rt_Finite_cells_iterator cit = regular.finite_cells_begin();
	 cit != regular.finite_cells_end();
	 cit++) {
      construct_3_cell(cit);
    }

    triangulation_incr_builder.end_triangulation();
    
    anchors.clear();

  }

  Tmc_Vertex_handle add_vertex(Symb_anchor const &anchor); 
  Tmc_Cell_handle add_cell(Tmc_Vertex_handle vh[], int orient, Rt_Simplex s);
	
  Tmc_Vertex_handle get_vertex(Rt_Simplex &sDel, Rt_Simplex &sVor);


  void construct_anchor_del(Rt_Simplex const &sDel);
  void construct_anchor_vor(Rt_Simplex const &sVor);
  void construct_anchors();
  Rt_Simplex &get_anchor_del(Rt_Simplex const &sDel) {
    typename Simplex_UF_map::iterator it = anchor_del_map.find(sDel);
    CGAL_assertion(it != anchor_del_map.end());
    return *anchor_del_uf.find(it->second);
  }
  Rt_Simplex &get_anchor_vor(Rt_Simplex const &sVor) {
    typename Simplex_UF_map::iterator it = anchor_vor_map.find(sVor);
    CGAL_assertion(it != anchor_vor_map.end());
    return *anchor_vor_uf.find(it->second);
  }  
//  Anchor_map_iterator find_anchor(Anchor_map &a_map, Rt_Simplex const&s) {
//    return find_anchor(a_map, a_map.find(s));
//  }
//  Anchor_map_iterator find_anchor(Anchor_map &a_map,
//				  Anchor_map_iterator const&it) {
//    CGAL_assertion(it != a_map.end());
//    Anchor_map_iterator it2 = it->second;
//    while (it2 != it2->second) {
//      it->second = it2->second;
//      // NGHK: changed the type for the map-iterator-hack
//      it2->second = it;
//      it2 = it->second;
//    }
//    return it2;
//  }
  void construct_vertices();
  
  Tmc_Point get_weighted_circumcenter(Rt_Simplex const &s);
  Tmc_Point get_anchor(Rt_Simplex const &sDel, Rt_Simplex const &sVor);
  template <class Point>
  Point construct_anchor_point(const Point &center_del, 
                               const Point &center_vor) {
    return construct_anchor_point_3_obj(center_del,center_vor);
  }
 
  void construct_0_cell(Rt_Vertex_handle rt_vh);
  void construct_1_cell(const Rt_Finite_edges_iterator &eit);
  void construct_2_cell(const Rt_Finite_facets_iterator &fit);
  void construct_3_cell(Rt_Cell_handle rt_ch);
	
  void remove_small_edges();
  bool is_collapsible(Tmc_Vertex_handle vh, 
		      Tmc_Vertex_handle &vh_collapse_to,
		      Tmc_RT sq_length);
  void do_collapse(Tmc_Vertex_handle vh, Tmc_Vertex_handle vh_collapse_to);
  

  Sign orientation(Tmc_Cell_handle ch);



private:

  Regular const &regular;
  Rt_FT const &shrink;
  Triangulated_mixed_complex &_tmc;
  Triangulated_mixed_complex_observer &observer;

  Triangulation_incremental_builder triangulation_incr_builder;

  Construct_weighted_circumcenter_3<
    Regular_triangulation_euclidean_traits_3<
    Triangulated_mixed_complex_traits> >                weighted_circumcenter_obj;

  Construct_anchor_point_3<
    Regular_triangulation_euclidean_traits_3<
    Triangulated_mixed_complex_traits> >                construct_anchor_point_3_obj;

  Compute_squared_radius_smallest_orthogonal_sphere_3<
    Regular_triangulation_euclidean_traits_3<
    Triangulated_mixed_complex_traits> >       orthoweight_obj;
  Compute_anchor_3<Regular> compute_anchor_obj;
  bool verbose;

  Weighted_converter_3<
    Cartesian_converter<typename Regular_traits::Bare_point::R, 
			Triangulated_mixed_complex_traits > >
  r2t_converter_object;
    

  static const int edge_index[4][4];
  struct Index_c4 { Tmc_Vertex_handle V[4]; };
  struct Index_c6 { Tmc_Vertex_handle V[6]; };
  struct Index_c44 { Tmc_Vertex_handle V[4][4]; };
  struct Index_v {
    Unique_hash_map < Rt_Vertex_handle, Tmc_Vertex_handle > V;
  };

  // index to vertex
  Unique_hash_map < Rt_Cell_handle, Index_c4 > index_03;
  
  Union_find_anchor                            anchor_del_uf, anchor_vor_uf;
  Simplex_UF_map                               anchor_del_map, anchor_vor_map;
                   
//  Anchor_map                                     anchor_del2, anchor_vor2;
  std::map<Symb_anchor, Tmc_Vertex_handle>        anchors;
};

template < 
  class RegularTriangulation_3,
  class TriangulatedMixedComplex_3,
  class TriangulatedMixedComplexObserver_3>
const int Mixed_complex_triangulator_3<
  RegularTriangulation_3,
  TriangulatedMixedComplex_3,
  TriangulatedMixedComplexObserver_3>::
edge_index[4][4] = {{-1,0,1,2},{0,-1,3,4},{1,3,-1,5},{2,4,5,-1}};


template < 
  class RegularTriangulation_3,
  class TriangulatedMixedComplex_3,
  class TriangulatedMixedComplexObserver_3>
void
Mixed_complex_triangulator_3<
  RegularTriangulation_3,
  TriangulatedMixedComplex_3,
  TriangulatedMixedComplexObserver_3>::
construct_anchor_del(Rt_Simplex const &sDel) {
  Rt_Simplex s = compute_anchor_obj.anchor_del(sDel);
  
  typename Union_find_anchor::handle sDel_handle, s_handle;
  sDel_handle = anchor_del_uf.make_set(sDel);
  anchor_del_map[sDel] = sDel_handle;
  
  typename Simplex_UF_map::iterator s_it = anchor_del_map.find(s); 
  CGAL_assertion(s_it != anchor_del_map.end());

  anchor_del_uf.unify_sets(sDel_handle, s_it->second);

  // degenerate simplices:
  if (compute_anchor_obj.is_degenerate()) {
    typename Compute_anchor::Simplex_iterator degenerate_it;
    typename Simplex_UF_map::iterator deg_map_it;
    for (degenerate_it = compute_anchor_obj.equivalent_anchors_begin();
         degenerate_it != compute_anchor_obj.equivalent_anchors_end(); 
         degenerate_it++) {
      deg_map_it = anchor_del_map.find(*degenerate_it); 
      CGAL_assertion(deg_map_it != anchor_del_map.end());
      
      anchor_del_uf.unify_sets(sDel_handle, deg_map_it->second);
    }
  }
}

template < 
  class RegularTriangulation_3,
  class TriangulatedMixedComplex_3,
  class TriangulatedMixedComplexObserver_3>
void
Mixed_complex_triangulator_3<
  RegularTriangulation_3,
  TriangulatedMixedComplex_3,
  TriangulatedMixedComplexObserver_3>::
construct_anchor_vor(Rt_Simplex const &sVor) {
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
//  Rt_Simplex s = compute_anchor_obj.anchor_vor(sVor);
//  anchor_vor2[sVor] = Anchor_map_iterator();
//
//  Anchor_map_iterator it = anchor_vor2.find(sVor);
//  Anchor_map_iterator it2 = anchor_vor2.find(s);
//  CGAL_assertion(it != anchor_vor2.end());
//  CGAL_assertion(it2 != anchor_vor2.end());
//  it->second = it2;
//
//  // degenerate simplices:
//  if (compute_anchor_obj.is_degenerate()) {
//    it = find_anchor(anchor_vor2, it);
//    typename Compute_anchor::Simplex_iterator degenerate_it;
//    for (degenerate_it = compute_anchor_obj.equivalent_anchors_begin();
//	 degenerate_it != compute_anchor_obj.equivalent_anchors_end(); 
//	 degenerate_it++) {
//      Anchor_map_iterator tmp;
//      it2 = anchor_vor2.find(*degenerate_it);
//      // Possibly not found for 2 Voronoi vertices with the same center,
//      // If the first vertex is inserted and the second is already found.
//      // see compute_anchor_obj.anchor_vor(Cell_handle)
//      if (it2 != anchor_vor2.end()) {
//	CGAL_assertion(it2 != anchor_vor2.end());
//	// Merge sets:
//	while (it2 != it2->second) {
//	  tmp = it2->second;
//	  it2->second = it->second;
//	  it2 = tmp;
//	  CGAL_assertion(it2 != anchor_vor2.end());
//	}
//	it2->second = it->second;
//      }
//    }
//  }
}

template < 
  class RegularTriangulation_3,
  class TriangulatedMixedComplex_3,
  class TriangulatedMixedComplexObserver_3>
void
Mixed_complex_triangulator_3<
  RegularTriangulation_3,
  TriangulatedMixedComplex_3,
  TriangulatedMixedComplexObserver_3>::
construct_anchors() {
  Rt_Finite_vertices_iterator vit;
  Rt_Finite_edges_iterator eit;
  Rt_Finite_facets_iterator fit;
  Rt_Finite_cells_iterator cit;
  Rt_Simplex s;
  
  // Compute anchor points:
  for (vit=regular.finite_vertices_begin();
       vit!=regular.finite_vertices_end(); vit++) {
    construct_anchor_del(Rt_Simplex(vit));
  }  
  for (eit=regular.finite_edges_begin();
       eit!=regular.finite_edges_end(); eit++) {
    s = Rt_Simplex(*eit);
    construct_anchor_del(s);
    CGAL_assertion(s.dimension() == 1);
  }
  for (fit=regular.finite_facets_begin();
       fit!=regular.finite_facets_end(); fit++) {
    s = Rt_Simplex(*fit);
    construct_anchor_del(s);
    CGAL_assertion(s.dimension() == 2);
  }
  for (cit=regular.finite_cells_begin();
       cit!=regular.finite_cells_end(); cit++) {
    s = Rt_Simplex(cit);
    construct_anchor_del(s);
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
template < 
  class RegularTriangulation_3,
  class TriangulatedMixedComplex_3,
  class TriangulatedMixedComplexObserver_3>
void
Mixed_complex_triangulator_3<
  RegularTriangulation_3,
  TriangulatedMixedComplex_3,
  TriangulatedMixedComplexObserver_3>::
construct_vertices() {
  Rt_All_cells_iterator acit;
  Rt_Finite_cells_iterator cit;
  Rt_Finite_facets_iterator fit;
  Rt_Finite_edges_iterator eit;
  Rt_Finite_vertices_iterator vit;
  Rt_Cell_circulator ccir, cstart;
  Rt_Vertex_handle v1, v2, v3;
  Rt_Edge e;
  Rt_Cell_handle c1, c2;
  Rt_Simplex sDel, sVor;
  Tmc_Vertex_handle vh;

  if (verbose) std::cout << "construct_anchors" << std::endl;
  construct_anchors();

  if (verbose) std::cout << "9 ";
  // anchor dimDel=0, dimVor=3
  for (cit=regular.finite_cells_begin();
       cit!=regular.finite_cells_end(); cit++) {
    sVor = get_anchor_vor(Rt_Simplex(cit));
    for (int i=0; i<4; i++) {
      sDel = get_anchor_del(Rt_Simplex(cit->vertex(i)));
      if (anchors.find(Symb_anchor(sDel,sVor)) == anchors.end()) {
	vh = add_vertex(Symb_anchor(sDel,sVor));
	anchors[Symb_anchor(sDel,sVor)] = vh;
	CGAL_assertion(vh == get_vertex(sDel, sVor));
      }
    }
  }

  if (verbose) std::cout << "8 ";
  // anchor dimDel=1, dimVor=3
  for (cit=regular.finite_cells_begin(); cit!=regular.finite_cells_end(); cit++) {
    sVor = get_anchor_vor(Rt_Simplex(cit));
    for (int i=0; i<3; i++) {
      for (int j=i+1; j<4; j++) {
	sDel = get_anchor_del(Rt_Simplex(Rt_Edge(cit,i,j)));
	if (anchors.find(Symb_anchor(sDel,sVor)) == anchors.end()) {
	  vh = add_vertex(Symb_anchor(sDel,sVor));
	  anchors[Symb_anchor(sDel,sVor)] = vh;
	  CGAL_assertion(vh == get_vertex(sDel, sVor));
	}
      }
    }
  }

  if (verbose) std::cout << "7 ";
  // anchor dimDel=2, dimVor=3 and dimDel=0, dimVor=2
  for (fit=regular.finite_facets_begin(); fit!=regular.finite_facets_end(); fit++) {
    // anchor dimDel=2, dimVor=3
    c1 = fit->first;
    c2 = c1->neighbor(fit->second);

    sDel = get_anchor_del(*fit);
    if (!regular.is_infinite(c1)) {
      sVor = get_anchor_vor(c1);
      if (anchors.find(Symb_anchor(sDel,sVor)) == anchors.end()) {
	vh = add_vertex(Symb_anchor(sDel,sVor));
	anchors[Symb_anchor(sDel,sVor)] = vh;
	CGAL_assertion(vh == get_vertex(sDel, sVor));
      }
    }
    if (!regular.is_infinite(c2)) {
      sVor = get_anchor_vor(c2);
      if (anchors.find(Symb_anchor(sDel,sVor)) == anchors.end()) {
	vh = add_vertex(Symb_anchor(sDel,sVor));
	anchors[Symb_anchor(sDel,sVor)] = vh;
	CGAL_assertion(vh == get_vertex(sDel, sVor));
      }
    }
    // anchor dimDel=0, dimVor=2
    sVor = get_anchor_vor(*fit);
    for (int i=1; i<4; i++) {
      sDel = get_anchor_del(Rt_Simplex(c1->vertex((fit->second+i)&3)));
      if (anchors.find(Symb_anchor(sDel,sVor)) == anchors.end()) {
	vh = add_vertex(Symb_anchor(sDel,sVor));
	anchors[Symb_anchor(sDel,sVor)] = vh;
	CGAL_assertion(vh == get_vertex(sDel, sVor));
      } else {
        vh = get_vertex(sDel, sVor);
      }
    }
  }
	
  if (verbose) std::cout << "6 ";
  // anchor dimDel=0, dimVor=1
  for (eit=regular.finite_edges_begin(); eit!=regular.finite_edges_end(); eit++) {
    sVor = get_anchor_vor(*eit);

    v1 = eit->first->vertex(eit->second);
    v2 = eit->first->vertex(eit->third);
    sDel = get_anchor_del(v1);
    if (anchors.find(Symb_anchor(sDel,sVor)) == anchors.end()) {
      vh = add_vertex(Symb_anchor(sDel,sVor));
      anchors[Symb_anchor(sDel,sVor)] = vh;
      CGAL_assertion(vh == get_vertex(sDel, sVor));
    }
			
    sDel = get_anchor_del(v2);
    if (anchors.find(Symb_anchor(sDel,sVor)) == anchors.end()) {
      vh = add_vertex(Symb_anchor(sDel,sVor));
      anchors[Symb_anchor(sDel,sVor)] = vh;
      CGAL_assertion(vh == get_vertex(sDel, sVor));
    }
  }
	
  if (verbose) std::cout << "5 ";
  // anchor dimDel=3, dimVor=3
  for (cit=regular.finite_cells_begin(); cit!=regular.finite_cells_end(); cit++) {
    sDel = get_anchor_del(Rt_Simplex(cit));
    sVor = get_anchor_vor(Rt_Simplex(cit));
    if (anchors.find(Symb_anchor(sDel,sVor)) == anchors.end()) {
      vh = add_vertex(Symb_anchor(sDel,sVor));
      anchors[Symb_anchor(sDel,sVor)] = vh;
      CGAL_assertion(vh == get_vertex(sDel, sVor));
    }
  }


  if (verbose) std::cout << "4 ";
  // anchor dimDel=0, dimVor=0
  for (vit=regular.finite_vertices_begin(); vit!=regular.finite_vertices_end(); vit++) {
    sDel = get_anchor_del(Rt_Simplex(vit));
    sVor = get_anchor_vor(Rt_Simplex(vit));
    if (anchors.find(Symb_anchor(sDel,sVor)) == anchors.end()) {
      vh = add_vertex(Symb_anchor(sDel,sVor));
      anchors[Symb_anchor(sDel,sVor)] = vh;
      CGAL_assertion(vh == get_vertex(sDel, sVor));
    }
  }
	
  if (verbose) std::cout << "3 ";
  // anchor dimDel=1, dimVor=2
  for (fit=regular.finite_facets_begin(); fit!=regular.finite_facets_end(); fit++) {
    c1 = fit->first;
    c2 = c1->neighbor(fit->second);

    sVor = get_anchor_vor(Rt_Simplex(*fit));
    for (int i=1; i<3; i++) {
      for (int j=i+1; j<4; j++) {
        e.first = c1;
        e.second = (fit->second+i)&3;
        e.third = (fit->second+j)&3;
        sDel = get_anchor_del(Rt_Simplex(e));
	if (anchors.find(Symb_anchor(sDel,sVor)) == anchors.end()) {
	  vh = add_vertex(Symb_anchor(sDel,sVor));
	  anchors[Symb_anchor(sDel,sVor)] = vh;
	  CGAL_assertion(vh == get_vertex(sDel, sVor));
	}
      }
    }
  }
	
  if (verbose) std::cout << "2 ";
  // anchor dimDel=2, dimVor=2
  for (fit=regular.finite_facets_begin(); fit!=regular.finite_facets_end(); fit++) {
    c1 = fit->first;
    c2 = c1->neighbor(fit->second);

    sVor = get_anchor_vor(Rt_Simplex(*fit));
    sDel = get_anchor_del(Rt_Simplex(*fit));
    if (anchors.find(Symb_anchor(sDel,sVor)) == anchors.end()) {
      vh = add_vertex(Symb_anchor(sDel,sVor));
      anchors[Symb_anchor(sDel,sVor)] = vh;
      CGAL_assertion(vh == get_vertex(sDel, sVor));
    }
  }
	
  if (verbose) std::cout << "1" << std::endl;
  // anchor dimDel=1, dimVor=1
  for (eit=regular.finite_edges_begin(); eit!=regular.finite_edges_end(); eit++) {
    v1 = eit->first->vertex(eit->second);
    v2 = eit->first->vertex(eit->third);

    sVor = get_anchor_vor(Rt_Simplex(*eit));
    sDel = get_anchor_del(Rt_Simplex(*eit));

    if (anchors.find(Symb_anchor(sDel,sVor)) == anchors.end()) {
      vh = add_vertex(Symb_anchor(sDel,sVor));
      anchors[Symb_anchor(sDel,sVor)] = vh;
      CGAL_assertion(vh == get_vertex(sDel, sVor));
    }
  }
}

// Constructs the cells of the mixed complex corresponding
// to Regular vertices
template < 
  class RegularTriangulation_3,
  class TriangulatedMixedComplex_3,
  class TriangulatedMixedComplexObserver_3>
void
Mixed_complex_triangulator_3<
  RegularTriangulation_3,
  TriangulatedMixedComplex_3,
  TriangulatedMixedComplexObserver_3>::
construct_0_cell(Rt_Vertex_handle rt_vh) {
  Rt_Simplex sDel_v, sVor_v, sVor_e, sVor_f, sVor_c;
  Tmc_Vertex_handle vh[4];
  
  Rt_Simplex simplex(rt_vh);
  sDel_v = get_anchor_del(Rt_Simplex(rt_vh));
  sVor_v = get_anchor_vor(Rt_Simplex(rt_vh));
  vh[0] = get_vertex(sDel_v,sVor_v);
    
  std::list<Rt_Cell_handle> adj_cells;
  typename std::list<Rt_Cell_handle>::iterator adj_cell;
  regular.incident_cells(rt_vh, std::back_inserter(adj_cells));
    
  // Construct cells:
  for (adj_cell = adj_cells.begin();
       adj_cell != adj_cells.end();
       adj_cell ++) {
    if (!regular.is_infinite(*adj_cell)) {
      sVor_c = get_anchor_vor(Rt_Simplex(*adj_cell));
      vh[3] = get_vertex(sDel_v,sVor_c);
      int index = (*adj_cell)->index(rt_vh);
      for (int i=1; i<4; i++) {
	sVor_f = get_anchor_vor(Rt_Simplex(Rt_Facet(*adj_cell,(index+i)&3)));
	vh[2] = get_vertex(sDel_v,sVor_f);
	  
	for (int j=1; j<4; j++) {
	  if (j!=i) {
	    sVor_e = get_anchor_vor(
				    Rt_Simplex(Rt_Edge(*adj_cell,index,(index+j)&3)));
	    vh[1] = get_vertex(sDel_v,sVor_e);
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

// Constructs 1-cells of the mixed complex corresponding to edges
// of the regular triangulation
template < 
  class RegularTriangulation_3,
  class TriangulatedMixedComplex_3,
  class TriangulatedMixedComplexObserver_3>
void
Mixed_complex_triangulator_3<
  RegularTriangulation_3,
  TriangulatedMixedComplex_3,
  TriangulatedMixedComplexObserver_3>::
construct_1_cell(const Rt_Finite_edges_iterator &e) {
  Rt_Simplex sDel_v, sDel_e, sVor_e, sVor_f, sVor_c;
  Tmc_Vertex_handle vh[4];
  Rt_Vertex_handle v[2];
  Tmc_Cell_handle ch;
  
  Rt_Simplex mixed_cell_simplex(*e);
  sDel_e = get_anchor_del(Rt_Simplex(*e));
  sVor_e = get_anchor_vor(Rt_Simplex(*e));
    
  v[0] = e->first->vertex(e->second);
  v[1] = e->first->vertex(e->third);
    
  // Construct cells on the side of v[vi]:
  for (int vi=0; vi<2; vi++) {
    sDel_v = get_anchor_del(Rt_Simplex(v[vi]));
    if (!(sDel_v == sDel_e)) {
      Rt_Cell_circulator ccir, cstart;
      ccir = cstart = regular.incident_cells(*e);
      do {
	if (!regular.is_infinite(ccir)) {
	  int index0 = ccir->index(v[vi]);
	  int index1 = ccir->index(v[1-vi]);

	  sVor_c = get_anchor_vor(Rt_Simplex(ccir));

	  for (int fi=1; fi<4; fi++) {
	    if (((index0+fi)&3) != index1) {
	      sVor_f =
		get_anchor_vor(Rt_Simplex(Rt_Facet(ccir,(index0+fi)&3)));
	      if ((sVor_c != sVor_f) && (sVor_f != sVor_e)) {
		vh[0] = get_vertex(sDel_v, sVor_e);
		vh[1] = get_vertex(sDel_e, sVor_e);
		vh[2] = get_vertex(sDel_e, sVor_f);
		vh[3] = get_vertex(sDel_e, sVor_c);
		int orient;
		if (((4+index1-index0)&3) == 1) {
		  orient = (index1 + (fi==2))&1;
		} else {
		  orient = (index1 + (fi==1))&1;
		}
		// vh: dimension are (01,11,12,13)
		ch = add_cell(vh,orient,mixed_cell_simplex);
									
		vh[1] = get_vertex(sDel_v, sVor_f);
		// vh: dimension are (01,02,12,13)
		ch = add_cell(vh,1-orient,mixed_cell_simplex);
									
		vh[2] = get_vertex(sDel_v, sVor_c);
		// vh: dimension are (01,02,03,13)
		ch = add_cell(vh,orient,mixed_cell_simplex);
	      }
	    }
	  }
	}
	ccir ++;
      } while (ccir != cstart);
    }
  }
}


// Constructs 2-cells of the mixed complex corresponding to facets
// of the regular triangulation
template < 
  class RegularTriangulation_3,
  class TriangulatedMixedComplex_3,
  class TriangulatedMixedComplexObserver_3>
void
Mixed_complex_triangulator_3<
  RegularTriangulation_3,
  TriangulatedMixedComplex_3,
  TriangulatedMixedComplexObserver_3>::
construct_2_cell(const Rt_Finite_facets_iterator &fit) {
  Rt_Simplex sDel_v, sDel_e, sDel_f, sVor_f, sVor_c;
  Tmc_Vertex_handle vh[4]; // Implicit function over vLabels is increasing ...
  Rt_Cell_handle rt_ch;
  int index;
	
  rt_ch = fit->first;
  index = fit->second;
  Rt_Simplex simplex(*fit);
  sDel_f = get_anchor_del(Rt_Simplex(*fit));
  sVor_f = get_anchor_vor(Rt_Simplex(*fit));
		
  for (int i=0; i<2; i++) { // Do this twice
    if (!regular.is_infinite(rt_ch)) {
      sVor_c = get_anchor_vor(Rt_Simplex(rt_ch));
	
      vh[3] = get_vertex(sDel_f, sVor_c);
      Tmc_Vertex_handle vh2 = get_vertex(sDel_f, sVor_f);
      if (vh2 != vh[3]) {
	// Facet and cell do not coincide ..
	for (int vi=1; vi<4; vi++) {
	  sDel_v = get_anchor_del(Rt_Simplex(rt_ch->vertex((index+vi)&3)));
	  //index_02[rt_ch].V[index][(index+vi)&3];
	  vh[0] = get_vertex(sDel_v, sVor_f);
	  for (int ei=1; ei<4; ei++) {
	    if (vi != ei) {
	      vh[2] = vh2;
	      int index0 = (index+vi)&3;
	      int index1 = (index+ei)&3;
	      int fi = (6+index-vi-ei)&3;//6-index-index0-index1;
	      sDel_e =
		get_anchor_del(Rt_Simplex(Rt_Edge(rt_ch, index0, index1)));
	      vh[1] = get_vertex(sDel_e, sVor_f);
	      //index_12[rt_ch].V[index][(6+index-vi-ei)&3];
	      if ((vh[0] != vh[1]) && (vh[1] != vh[2])) {
		// index0: v0
		// index1: v1
		// index0+fi&3 == facet
		int orient;
									
		if (((4+index1-index0)&3) == 3) {
		  orient = (index1 + (((4+index0-fi)&3)==2))&1;
		} else {
		  orient = (index1 + (((4+index0-fi)&3)==1))&1;
		}

		add_cell(vh,orient,simplex);
									
		vh[2] = get_vertex(sDel_e, sVor_c);
		add_cell(vh,1-orient,simplex);
									
		vh[1] = get_vertex(sDel_v, sVor_c);
		add_cell(vh,orient,simplex);
	      } 
	    }
	  }
	}
      }
    }
    // swap to the other cell
    Rt_Cell_handle ch_old = rt_ch;
    rt_ch = rt_ch->neighbor(index);
    index = rt_ch->index(ch_old);
  }

  CGAL_assertion(rt_ch == fit->first);
  CGAL_assertion(index == fit->second);
}


// Constructs 3-cells of the mixed complex corresponding to cells
// of the regular triangulation
template < 
  class RegularTriangulation_3,
  class TriangulatedMixedComplex_3,
  class TriangulatedMixedComplexObserver_3>
void
Mixed_complex_triangulator_3<
  RegularTriangulation_3,
  TriangulatedMixedComplex_3,
  TriangulatedMixedComplexObserver_3>::
construct_3_cell(Rt_Cell_handle rt_ch) {
  Rt_Simplex sDel_v, sDel_e, sDel_f, sDel_c, sVor_c;
  Tmc_Vertex_handle vh[4];
  Tmc_Cell_handle ch;

  // construct the tetrahedron:
  //   C[ch], C[Facet(ch,fi)], C[Edge(ch,ei,vi)], C[ch->vertex(vi)]
  sDel_c = get_anchor_del(Rt_Simplex(rt_ch));
  sVor_c = get_anchor_vor(Rt_Simplex(rt_ch));
  Rt_Simplex simplex = Rt_Simplex(rt_ch);
  vh[0] = get_vertex(sDel_c, sVor_c); 
  for (int fi=0; fi<4; fi++) {
    sDel_f = get_anchor_del(Rt_Simplex(Rt_Facet(rt_ch, fi)));
    vh[1] = get_vertex(sDel_f, sVor_c);
    if (vh[0] != vh[1]) {
      for (int vi=1; vi<4; vi++) {
	int index0 = (fi+vi)&3;
	sDel_v = get_anchor_del(Rt_Simplex(rt_ch->vertex(index0)));
	for (int ei=1; ei<4; ei++) {
	  int index1 = (fi+ei)&3;
	  if (vi != ei) {
	    sDel_e = get_anchor_del(Rt_Simplex(Rt_Edge(rt_ch, index0, index1)));
	    vh[2] = get_vertex(sDel_e, sVor_c);
	    // index_13[rt_ch].V[edge_index[index0][index1]];
	    vh[3] = get_vertex(sDel_v, sVor_c);
	    // index_03[rt_cit].V[index0];
	    if ((vh[1] != vh[2]) && (vh[2] != vh[3])) {
	      int orient;
								
	      if (((4+index1-index0)&3) == 1) {
		orient = (index1 + (vi==2))&1;
	      } else {
		orient = (index1 + (vi==3))&1;
	      }
	      ch = add_cell(vh, orient, simplex);
	    }
	  }
	}
      }
    }
  }
}

// Adds a vertex to the simplicial complex
template < 
  class RegularTriangulation_3,
  class TriangulatedMixedComplex_3,
  class TriangulatedMixedComplexObserver_3>
typename Mixed_complex_triangulator_3<
  RegularTriangulation_3,
  TriangulatedMixedComplex_3,
  TriangulatedMixedComplexObserver_3>::Tmc_Vertex_handle
Mixed_complex_triangulator_3<
  RegularTriangulation_3,
  TriangulatedMixedComplex_3,
  TriangulatedMixedComplexObserver_3>::
add_vertex (Symb_anchor const &anchor)
{
  Tmc_Vertex_handle vh;
  vh = triangulation_incr_builder.add_vertex();
  observer.after_vertex_insertion(anchor.first, anchor.second, vh); 
  
  Protect_FPU_rounding<true> P;
  vh->point() = get_anchor(anchor.first, anchor.second);

//   std::cout << "@ [" 
//             << vh->info().first << " - " 
//             << vh->info().second << "] -- [" 
//             << vh->point() << "] -- ["
//             << get_weighted_circumcenter(vh->info().first) << " - "
//             << get_weighted_circumcenter(vh->info().second) 
//             << "]" << std::endl;

  return vh;
}

// Gets a vertex from the simplicial complex based on the anchors
template < 
  class RegularTriangulation_3,
  class TriangulatedMixedComplex_3,
  class TriangulatedMixedComplexObserver_3>
typename Mixed_complex_triangulator_3<
  RegularTriangulation_3,
  TriangulatedMixedComplex_3,
  TriangulatedMixedComplexObserver_3>::Tmc_Vertex_handle
Mixed_complex_triangulator_3<
  RegularTriangulation_3,
  TriangulatedMixedComplex_3,
  TriangulatedMixedComplexObserver_3>::get_vertex (
						   Rt_Simplex &sDel, Rt_Simplex &sVor)
{
  Rt_Simplex sDel2 = get_anchor_del(sDel);
  Rt_Simplex sVor2 = get_anchor_vor(sVor);
  CGAL_assertion(sDel == sDel2);
  CGAL_assertion(sVor == sVor2);
  Tmc_Vertex_handle vh = anchors[Symb_anchor(sDel2,sVor2)];
  CGAL_assertion(vh != Tmc_Vertex_handle());
  return vh;
}

// Adds a cell to the simplicial complex
template < 
  class RegularTriangulation_3,
  class TriangulatedMixedComplex_3,
  class TriangulatedMixedComplexObserver_3>
typename Mixed_complex_triangulator_3<
  RegularTriangulation_3,
  TriangulatedMixedComplex_3,
  TriangulatedMixedComplexObserver_3>::Tmc_Cell_handle
Mixed_complex_triangulator_3<
  RegularTriangulation_3,
  TriangulatedMixedComplex_3,
  TriangulatedMixedComplexObserver_3>::
add_cell(Tmc_Vertex_handle vh[], int orient, Rt_Simplex s) {
  CGAL_assertion((orient==0) || (orient==1));
  CGAL_assertion(vh[0] != Tmc_Vertex_handle()); 
  CGAL_assertion(vh[1] != Tmc_Vertex_handle());
  CGAL_assertion(vh[2] != Tmc_Vertex_handle()); 
  CGAL_assertion(vh[3] != Tmc_Vertex_handle());
  CGAL_assertion(vh[0] != vh[1]); 
  CGAL_assertion(vh[0] != vh[2]); 
  CGAL_assertion(vh[0] != vh[3]); 
  CGAL_assertion(vh[1] != vh[2]); 
  CGAL_assertion(vh[1] != vh[3]); 
  CGAL_assertion(vh[2] != vh[3]);

  Tmc_Cell_handle ch;

  if (orient) {
    ch = triangulation_incr_builder.add_cell(vh[0], vh[1], vh[2], vh[3]);
  } else {
    ch = triangulation_incr_builder.add_cell(vh[0], vh[1], vh[3], vh[2]);
  }
  CGAL_assertion(orientation(ch) == POSITIVE);
  observer.after_cell_insertion(s, ch);
  return ch;
}

template < 
  class RegularTriangulation_3,
  class TriangulatedMixedComplex_3,
  class TriangulatedMixedComplexObserver_3>
typename Mixed_complex_triangulator_3<
  RegularTriangulation_3,
  TriangulatedMixedComplex_3,
  TriangulatedMixedComplexObserver_3>::Tmc_Point
Mixed_complex_triangulator_3<
  RegularTriangulation_3,
  TriangulatedMixedComplex_3,
  TriangulatedMixedComplexObserver_3>::
get_weighted_circumcenter(Rt_Simplex const &s) {
  Rt_Vertex_handle vh;
  Rt_Edge           e;
  Rt_Facet          f;
  Rt_Cell_handle   ch;
  
  Tmc_Point result;
  switch (s.dimension()) {
  case 0:
    vh=s;
    result = r2t_converter_object(vh->point());
    break;
  case 1:
    e=s;
    result = weighted_circumcenter_obj(
			     r2t_converter_object(e.first->vertex(e.second)->point()),
			     r2t_converter_object(e.first->vertex(e.third)->point()));
    break;
  case 2:
    f=s;
    result = weighted_circumcenter_obj(
			     r2t_converter_object(f.first->vertex((f.second+1)&3)->point()),
			     r2t_converter_object(f.first->vertex((f.second+2)&3)->point()),
			     r2t_converter_object(f.first->vertex((f.second+3)&3)->point()));
    break;
  case 3:
    ch=s;
    result = weighted_circumcenter_obj(
           r2t_converter_object(ch->vertex(0)->point()),
           r2t_converter_object(ch->vertex(1)->point()),
           r2t_converter_object(ch->vertex(2)->point()),
           r2t_converter_object(ch->vertex(3)->point()));
    break;
  default:
    CGAL_error();
  }
  return result;
}

template < 
  class RegularTriangulation_3,
  class TriangulatedMixedComplex_3,
  class TriangulatedMixedComplexObserver_3>
typename Mixed_complex_triangulator_3<
  RegularTriangulation_3,
  TriangulatedMixedComplex_3,
  TriangulatedMixedComplexObserver_3>::Tmc_Point
Mixed_complex_triangulator_3<
  RegularTriangulation_3,
  TriangulatedMixedComplex_3,
  TriangulatedMixedComplexObserver_3>::
get_anchor(Rt_Simplex const &sDel, Rt_Simplex const &sVor)
{
  Protect_FPU_rounding<true> P;

  Tmc_Point dfoc = get_weighted_circumcenter(sDel);
  Tmc_Point vfoc = get_weighted_circumcenter(sVor);
	
  return construct_anchor_point(dfoc, vfoc);
}

template < 
  class RegularTriangulation_3,
  class TriangulatedMixedComplex_3,
  class TriangulatedMixedComplexObserver_3>
void
Mixed_complex_triangulator_3<
  RegularTriangulation_3,
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
  for (Tmc_Finite_vertices_iterator vit = _tmc.finite_vertices_begin();
       vit != _tmc.finite_vertices_end(); ) {
    vh = vit;
    vit++;
    if (is_collapsible(vh, vh_collapse_to,sq_length)) {
      do_collapse(vh,vh_collapse_to);
    }
  }
}

template < 
  class RegularTriangulation_3,
  class TriangulatedMixedComplex_3,
  class TriangulatedMixedComplexObserver_3>
bool
Mixed_complex_triangulator_3<
  RegularTriangulation_3,
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

template < 
  class RegularTriangulation_3,
  class TriangulatedMixedComplex_3,
  class TriangulatedMixedComplexObserver_3>
void
Mixed_complex_triangulator_3<
  RegularTriangulation_3,
  TriangulatedMixedComplex_3,
  TriangulatedMixedComplexObserver_3>::
do_collapse(Tmc_Vertex_handle vh, 
	    Tmc_Vertex_handle vh_collapse_to)
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

template < 
  class RegularTriangulation_3,
  class TriangulatedMixedComplex_3,
  class TriangulatedMixedComplexObserver_3>
Sign 
Mixed_complex_triangulator_3<
  RegularTriangulation_3,
  TriangulatedMixedComplex_3,
  TriangulatedMixedComplexObserver_3>::
orientation(Tmc_Cell_handle ch) {
    Orientation o;
    // Protection is outside the try block as VC8 has the CGAL_CFG_FPU_ROUNDING_MODE_UNWINDING_VC_BUG
    Protect_FPU_rounding<true> P;
    try {
    Tmc_Point pts[4];
    for (int i=0; i<4; i++) pts[i] = ch->vertex(i)->point();

    
    // filtered kernel
    o = _tmc.geom_traits().orientation_3_object()(pts[0], pts[1], 
                                                  pts[2], pts[3]);
  } catch (Uncertain_conversion_exception) {
    Protect_FPU_rounding<false> P(CGAL_FE_TONEAREST);
    typedef Exact_predicates_exact_constructions_kernel EK;
    typedef Cartesian_converter<EK, Tmc_traits>         Exact_converter;
    typedef Skin_surface_traits_3<EK>                   Exact_traits;

    Exact_converter converter;
    Exact_traits    exact_traits(shrink);
    typename EK::Point_3 e_pts[4];

    for (int k=0; k<4; k++) {
      e_pts[k] = 
        Skin_surface_base_3<Regular_traits>::
        get_anchor_point(ch->vertex(k)->info(), exact_traits);
      // Store the more precise point
      ch->vertex(k)->point() = converter(e_pts[k]);
    }
    o = exact_traits.orientation_3_object()(e_pts[0], e_pts[1], 
                                            e_pts[2], e_pts[3]);
  }
  return o;
}

template < 
  class RegularTriangulation_3,
  class TriangulatedMixedComplex_3,
  class TriangulatedMixedComplexObserver_3>
void 
triangulate_mixed_complex_3(RegularTriangulation_3 &rt,
			    typename RegularTriangulation_3::Geom_traits::FT
			    const & shrink_factor,
			    TriangulatedMixedComplex_3 &tmc,
			    TriangulatedMixedComplexObserver_3 &observer,
			    bool verbose) 
{
  typedef Mixed_complex_triangulator_3<
    RegularTriangulation_3,
    TriangulatedMixedComplex_3,
    TriangulatedMixedComplexObserver_3>  Mixed_complex_triangulator;
  Mixed_complex_triangulator(rt, shrink_factor, tmc, observer, verbose);
}


template < 
  class RegularTriangulation_3,
  class TriangulatedMixedComplex_3>
void 
triangulate_mixed_complex_3(RegularTriangulation_3 const &regular, 
			    typename RegularTriangulation_3::Geom_traits::FT
			    const &shrink_factor,
			    TriangulatedMixedComplex_3 &tmc,
			    bool verbose)
{
  Triangulated_mixed_complex_observer_3<
    TriangulatedMixedComplex_3, const RegularTriangulation_3> 
    observer(shrink_factor);
  triangulate_mixed_complex_3(regular, shrink_factor, tmc, observer, verbose);
}

} //namespace CGAL

#endif // CGAL_TRIANGULATE_MIXED_COMPLEX_H
