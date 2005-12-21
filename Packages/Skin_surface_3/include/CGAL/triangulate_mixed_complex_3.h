// Copyright (c) 1999-2003  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Nico Kruithof <Nico@cs.rug.nl>

#ifndef TRIANGULATE_MIXED_COMPLEX_3
#define TRIANGULATE_MIXED_COMPLEX_3

// #include <CGAL/Unique_hash_map.h>
#include <CGAL/Compute_anchor_3.h>

#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Triangulated_mixed_complex_observer_3.h>
#include <CGAL/Triangulation_incremental_builder_3.h>
// NGHK: move this one to SkinSurfaceTraits
#include <CGAL/Compute_anchor_3.h>

#include <CGAL/Union_find.h>

CGAL_BEGIN_NAMESPACE


template < 
  class SkinSurfaceTraits_3,
  class RegularTriangulation_3,
  class TriangulatedMixedComplex_3,
  class TriangulatedMixedComplexObserver_3 =
  Triangulated_mixed_complex_observer_3<SkinSurfaceTraits_3,
					TriangulatedMixedComplex_3,
					RegularTriangulation_3> >
class Mixed_complex_triangulator_3 {
public:
  typedef SkinSurfaceTraits_3                 Skin_traits_3;
  typedef typename SkinSurfaceTraits_3::Regular_traits
                                              Regular_traits;
  typedef typename SkinSurfaceTraits_3::Triangulated_mixed_complex_traits
                                              Triangulated_mixed_complex_kernel;

  typedef RegularTriangulation_3                   Regular;
  typedef TriangulatedMixedComplex_3               Triangulated_mixed_complex;
  typedef TriangulatedMixedComplexObserver_3
  Triangulated_mixed_complex_observer;

  typedef typename Skin_traits_3::R2T_converter R2T_converter;
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
  typedef typename Skin_traits_3::Regular_kernel          Rt_Geom_traits;
  typedef typename Regular::Bare_point               Rt_Point;
  typedef typename Regular::Geom_traits::RT          Rt_RT;
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
	
  typedef typename Skin_traits_3::Triangulated_mixed_complex_traits
    Tmc_Geom_traits;
  typedef typename Tmc_Geom_traits::Point_3              Tmc_Point;
  typedef typename Tmc_Geom_traits::RT                   Tmc_RT;

  typedef Triangulation_incremental_builder_3<Triangulated_mixed_complex>
  Triangulation_incremental_builder;

  typedef Compute_anchor_3<Regular>                       Compute_anchor;
  typedef std::pair<Rt_Simplex,Rt_Simplex>                Symb_anchor;
public:

  Mixed_complex_triangulator_3(
    Regular &regular,
    Triangulated_mixed_complex &triangulated_mixed_complex,
    Skin_traits_3 &skin_traits)
    : regular(regular),
      triangulated_mixed_complex(triangulated_mixed_complex),
      skin_traits(skin_traits),
      observer(skin_traits.shrink_factor()), 
      triangulation_incr_builder(triangulated_mixed_complex), 
      compute_anchor_obj(regular) {

    build();
  }

  Mixed_complex_triangulator_3(
    Regular &regular,
    Triangulated_mixed_complex &triangulated_mixed_complex,
    Skin_traits_3 &skin_traits,
    Triangulated_mixed_complex_observer &observer)
    : regular(regular),
      triangulated_mixed_complex(triangulated_mixed_complex),
      skin_traits(skin_traits),
      observer(observer),
      triangulation_incr_builder(triangulated_mixed_complex), 
      compute_anchor_obj(regular) {
	  
    build();
  }


private:
  void build() {
    triangulation_incr_builder.begin_triangulation(3);

    construct_vertices();

    construct_0_cells(); // mixed cells corresponding to regular vertices
    construct_1_cells(); // mixed cells corresponding to regular edges
    construct_2_cells(); // mixed cells corresponding to regular facets
    construct_3_cells(); // mixed cells corresponding to regular cells

    triangulation_incr_builder.end_triangulation();
  }

  Tmc_Vertex_handle add_vertex(Symb_anchor const &anchor); 
  Tmc_Cell_handle add_cell(Tmc_Vertex_handle vh[], int orient, Rt_Simplex s);
	
  Tmc_Vertex_handle get_vertex(Rt_Simplex &sDel, Rt_Simplex &sVor);


  void construct_anchor_del(Rt_Simplex const &sDel);
  void construct_anchor_vor(Rt_Simplex const &sVor);
  void construct_anchors();
  Rt_Simplex get_anchor_del(Rt_Simplex const &sDel) {
    return *anchor_del.find(map_del[sDel]);
  }
  Rt_Simplex get_anchor_vor(Rt_Simplex const &sVor) {
    return *anchor_vor.find(map_vor[sVor]);
  }  
  void construct_vertices();
  
  Tmc_Point get_orthocenter(Rt_Simplex const &s);
  Tmc_Point get_anchor(Rt_Simplex const &sDel, Rt_Simplex const &sVor);

  void construct_0_cells();
  void construct_1_cells();
  void construct_2_cells();
  void construct_3_cells();
	
private:
  Regular &regular;
  Triangulated_mixed_complex &triangulated_mixed_complex;
  Skin_traits_3 &skin_traits;
  Triangulated_mixed_complex_observer &observer;

  Triangulation_incremental_builder triangulation_incr_builder;

  Construct_weighted_circumcenter_3<
    Regular_triangulation_euclidean_traits_3<
    Tmc_Geom_traits> >                                orthocenter_obj;
  Compute_squared_radius_smallest_orthogonal_sphere_3<
    Regular_triangulation_euclidean_traits_3<
    typename Triangulated_mixed_complex::Geom_traits> >              orthoweight_obj;
  Compute_anchor_3<Regular> compute_anchor_obj;

  const static int edge_index[4][4];
  struct Index_c4 { Tmc_Vertex_handle V[4]; };
  struct Index_c6 { Tmc_Vertex_handle V[6]; };
  struct Index_c44 { Tmc_Vertex_handle V[4][4]; };
  struct Index_v {
    Unique_hash_map < Rt_Vertex_handle, Tmc_Vertex_handle > V;
  };
  // Facets on the border of the simplicial complex:
  // name is given by (dim del,dim vor)

  // index to vertex
  Unique_hash_map < Rt_Cell_handle, Index_c4 > index_03;
  
  typedef Union_find<Rt_Simplex>                 Union_find_anchor;
  typedef typename Union_find_anchor::handle     Union_find_anchor_handle;
  typedef typename Union_find_anchor::iterator   Union_find_anchor_iterator;

  Union_find_anchor                              anchor_del, anchor_vor;
  std::map<Rt_Simplex, Union_find_anchor_handle> map_del, map_vor;
  std::map<Symb_anchor, Tmc_Vertex_handle>        anchors;
};

template < class SkinSurfaceTraits_3, class TriangulatedMixedComplex,
	   class Regular_TDS, class Mixed_complex_observer_>
const int Mixed_complex_triangulator_3<
  SkinSurfaceTraits_3, TriangulatedMixedComplex,
  Regular_TDS, Mixed_complex_observer_>::edge_index[4][4] =
  {{-1,0,1,2},{0,-1,3,4},{1,3,-1,5},{2,4,5,-1}};


template < class SkinSurfaceTraits_3, class TriangulatedMixedComplex,
	   class Regular_TDS, class Mixed_complex_observer_>
void
Mixed_complex_triangulator_3<SkinSurfaceTraits_3, TriangulatedMixedComplex,
			     Regular_TDS, Mixed_complex_observer_>::
construct_anchor_del(Rt_Simplex const &sDel) {
  Rt_Simplex sim = sDel;
  Union_find_anchor_handle handle = anchor_del.make_set(sDel);
  map_del[sDel] = handle;
  CGAL_assertion(sim == sDel);
  
  Rt_Simplex s = compute_anchor_obj.anchor_del(sDel);
  CGAL_assertion(sim == sDel);
  if (sDel != s) {
    CGAL_assertion(s != sim);
    anchor_del.unify_sets(handle, map_del[s]);
  }

  // degenerate simplices:
  if (compute_anchor_obj.is_degenerate()) {
    typename Compute_anchor::Simplex_iterator it;
    for (it = compute_anchor_obj.equivalent_anchors_begin();
	 it != compute_anchor_obj.equivalent_anchors_end(); it++) {
      anchor_del.unify_sets(handle, map_del[*it]);
    }
  }
}

template < class SkinSurfaceTraits_3, class TriangulatedMixedComplex,
	   class Regular_TDS, class Mixed_complex_observer_>
void
Mixed_complex_triangulator_3<SkinSurfaceTraits_3, TriangulatedMixedComplex,
			     Regular_TDS, Mixed_complex_observer_>::
construct_anchor_vor(Rt_Simplex const &sVor) {
  Union_find_anchor_handle handle = anchor_vor.make_set(sVor);
  map_vor[sVor] = handle;

  Rt_Simplex s = compute_anchor_obj.anchor_vor(sVor);
  if (sVor != s) {
    anchor_vor.unify_sets(handle, map_vor[s]);
  }

  // degenerate simplices:
  if (compute_anchor_obj.is_degenerate()) {
    typename Compute_anchor::Simplex_iterator it;
    for (it = compute_anchor_obj.equivalent_anchors_begin();
	 it != compute_anchor_obj.equivalent_anchors_end(); it++) {
      typename std::map<Rt_Simplex, Union_find_anchor_handle>::iterator h_it;
      h_it = map_vor.find(*it);
      // Possibly not found for 2 Voronoi vertices with the same center,
      // If the first vertex is inserted and the second is already found.
      if (h_it != map_vor.end()) {
	anchor_vor.unify_sets(handle, (*h_it).second);
      } else {
	CGAL_assertion(s.dimension() == 3);
      }
    }
  }
}

template < class SkinSurfaceTraits_3, class TriangulatedMixedComplex,
	   class Regular_TDS, class Mixed_complex_observer_>
void
Mixed_complex_triangulator_3<SkinSurfaceTraits_3, TriangulatedMixedComplex,
			     Regular_TDS, Mixed_complex_observer_>::
construct_anchors() {
  Rt_Finite_vertices_iterator vit;
  Rt_Finite_edges_iterator eit;
  Rt_Finite_facets_iterator fit;
  Rt_Finite_cells_iterator cit;
  Rt_Simplex s;
  
  // Compute anchor points:
  for (vit=regular.finite_vertices_begin();
       vit!=regular.finite_vertices_end(); vit++) {
    s = Rt_Simplex(vit);
    construct_anchor_del(s);
    CGAL_assertion(s.dimension() == 0);
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
    s = Rt_Simplex(vit);
    construct_anchor_vor(s);
    CGAL_assertion(s.dimension() == 0);
  }
}


// Constructs the vertices of the simplicial complex
template < class SkinSurfaceTraits_3, class TriangulatedMixedComplex,
	   class Regular_TDS, class Mixed_complex_observer_>
void
Mixed_complex_triangulator_3<SkinSurfaceTraits_3, TriangulatedMixedComplex,
			     Regular_TDS, Mixed_complex_observer_>::
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

  construct_anchors();

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

  // anchor dimDel=1, dimVor=3
  for (cit=regular.finite_cells_begin(); cit!=regular.finite_cells_end(); cit++) {
    sVor = get_anchor_vor(Rt_Simplex(cit));
    for (int i=0; i<3; i++) {
      for (int j=i+1; j<4; j++) {
	sDel = get_anchor_del(Rt_Simplex(Rt_Edge(cit,i,j)));
	if (anchors.find(Symb_anchor(sDel,sVor)) == anchors.end()) {
	  vh = add_vertex(Symb_anchor(sDel,sVor));
	  anchors[Symb_anchor(sDel,sVor)] = vh;
	  assert(vh == get_vertex(sDel, sVor));
	}
      }
    }
  }

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
	assert(vh == get_vertex(sDel, sVor));
      }
    }
    if (!regular.is_infinite(c2)) {
      sVor = get_anchor_vor(c2);
      if (anchors.find(Symb_anchor(sDel,sVor)) == anchors.end()) {
	vh = add_vertex(Symb_anchor(sDel,sVor));
	anchors[Symb_anchor(sDel,sVor)] = vh;
	assert(vh == get_vertex(sDel, sVor));
      }
    }
    // anchor dimDel=0, dimVor=2
    sVor = get_anchor_vor(*fit);
    for (int i=1; i<4; i++) {
      sDel = get_anchor_del(Rt_Simplex(c1->vertex((fit->second+i)&3)));
      if (anchors.find(Symb_anchor(sDel,sVor)) == anchors.end()) {
	vh = add_vertex(Symb_anchor(sDel,sVor));
	anchors[Symb_anchor(sDel,sVor)] = vh;
	assert(vh == get_vertex(sDel, sVor));
      } else {
        vh = get_vertex(sDel, sVor);
      }
    }
  }
	
  // anchor dimDel=0, dimVor=1
  for (eit=regular.finite_edges_begin(); eit!=regular.finite_edges_end(); eit++) {
    sVor = get_anchor_vor(*eit);

    v1 = eit->first->vertex(eit->second);
    v2 = eit->first->vertex(eit->third);
    sDel = get_anchor_del(v1);
    if (anchors.find(Symb_anchor(sDel,sVor)) == anchors.end()) {
      vh = add_vertex(Symb_anchor(sDel,sVor));
      anchors[Symb_anchor(sDel,sVor)] = vh;
      assert(vh == get_vertex(sDel, sVor));
    }
			
    sDel = get_anchor_del(v2);
    if (anchors.find(Symb_anchor(sDel,sVor)) == anchors.end()) {
      vh = add_vertex(Symb_anchor(sDel,sVor));
      anchors[Symb_anchor(sDel,sVor)] = vh;
      assert(vh == get_vertex(sDel, sVor));
    }
  }
	
  // anchor dimDel=3, dimVor=3
  for (cit=regular.finite_cells_begin(); cit!=regular.finite_cells_end(); cit++) {
    sDel = get_anchor_del(Rt_Simplex(cit));
    sVor = get_anchor_vor(Rt_Simplex(cit));
    if (anchors.find(Symb_anchor(sDel,sVor)) == anchors.end()) {
      vh = add_vertex(Symb_anchor(sDel,sVor));
      anchors[Symb_anchor(sDel,sVor)] = vh;
      assert(vh == get_vertex(sDel, sVor));
    }
  }


  // anchor dimDel=0, dimVor=0
  for (vit=regular.finite_vertices_begin(); vit!=regular.finite_vertices_end(); vit++) {
    sDel = get_anchor_del(Rt_Simplex(vit));
    sVor = get_anchor_vor(Rt_Simplex(vit));
    if (anchors.find(Symb_anchor(sDel,sVor)) == anchors.end()) {
      vh = add_vertex(Symb_anchor(sDel,sVor));
      anchors[Symb_anchor(sDel,sVor)] = vh;
      assert(vh == get_vertex(sDel, sVor));
    }
  }
	
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
	  assert(vh == get_vertex(sDel, sVor));
	}
      }
    }
  }
	
  // anchor dimDel=2, dimVor=2
  for (fit=regular.finite_facets_begin(); fit!=regular.finite_facets_end(); fit++) {
    c1 = fit->first;
    c2 = c1->neighbor(fit->second);

    sVor = get_anchor_vor(Rt_Simplex(*fit));
    sDel = get_anchor_del(Rt_Simplex(*fit));
    if (anchors.find(Symb_anchor(sDel,sVor)) == anchors.end()) {
      vh = add_vertex(Symb_anchor(sDel,sVor));
      anchors[Symb_anchor(sDel,sVor)] = vh;
      assert(vh == get_vertex(sDel, sVor));
    }
  }
	
  // anchor dimDel=1, dimVor=1
  for (eit=regular.finite_edges_begin(); eit!=regular.finite_edges_end(); eit++) {
    v1 = eit->first->vertex(eit->second);
    v2 = eit->first->vertex(eit->third);

    sVor = get_anchor_vor(Rt_Simplex(*eit));
    sDel = get_anchor_del(Rt_Simplex(*eit));

    if (anchors.find(Symb_anchor(sDel,sVor)) == anchors.end()) {
      vh = add_vertex(Symb_anchor(sDel,sVor));
      anchors[Symb_anchor(sDel,sVor)] = vh;
      assert(vh == get_vertex(sDel, sVor));
    }
  }
}

// Constructs the cells of the mixed complex corresponding
// to Regular vertices
template < class SkinSurfaceTraits_3, class TriangulatedMixedComplex,
	   class Regular_TDS, class Mixed_complex_observer_>
void
Mixed_complex_triangulator_3<SkinSurfaceTraits_3, TriangulatedMixedComplex,
			     Regular_TDS, Mixed_complex_observer_>::
construct_0_cells() {
  Rt_Simplex sDel_v, sVor_v, sVor_e, sVor_f, sVor_c;
  Tmc_Vertex_handle vh[4];
  
  for (Rt_Finite_vertices_iterator vit=regular.finite_vertices_begin();
       vit!=regular.finite_vertices_end(); vit++) {
    
    Rt_Simplex simplex(vit);
    sDel_v = get_anchor_del(Rt_Simplex(vit));
    sVor_v = get_anchor_vor(Rt_Simplex(vit));
    vh[0] = get_vertex(sDel_v,sVor_v);
    
    std::list<Rt_Cell_handle> adj_cells;
    typename std::list<Rt_Cell_handle>::iterator adj_cell;
    regular.incident_cells(vit, std::back_inserter(adj_cells));
    
    // Construct cells:
    for (adj_cell = adj_cells.begin();
	 adj_cell != adj_cells.end();
	 adj_cell ++) {
      if (!regular.is_infinite(*adj_cell)) {
	sVor_c = get_anchor_vor(Rt_Simplex(*adj_cell));
	vh[3] = get_vertex(sDel_v,sVor_c);
	int index = (*adj_cell)->index(vit);
	for (int i=1; i<4; i++) {
	  sVor_f = get_anchor_vor(
	    Rt_Simplex(Rt_Facet(*adj_cell,(index+i)&3)));
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
		Tmc_Cell_handle ch =
		  add_cell(vh,(index + (j==(i%3+1)? 1:0))&1,simplex);
	      }
	    }
	  }
	}
      }
    }
  }
}

// Constructs 1-cells of the mixed complex corresponding to edges
// of the regular triangulation
template < class SkinSurfaceTraits_3, class TriangulatedMixedComplex,
	   class Regular_TDS, class Mixed_complex_observer_>
void
Mixed_complex_triangulator_3<SkinSurfaceTraits_3, TriangulatedMixedComplex,
			     Regular_TDS, Mixed_complex_observer_>::construct_1_cells() {
  Rt_Simplex sDel_v, sDel_e, sVor_e, sVor_f, sVor_c;
  Tmc_Vertex_handle vh[4];
  Rt_Vertex_handle v[2];
  Tmc_Cell_handle ch;
  
  for (Rt_Finite_edges_iterator eit=regular.finite_edges_begin();
       eit!=regular.finite_edges_end(); eit++) {
    Rt_Simplex mixed_cell_simplex(*eit);
    sDel_e = get_anchor_del(Rt_Simplex(*eit));
    sVor_e = get_anchor_vor(Rt_Simplex(*eit));
    
    v[0] = eit->first->vertex(eit->second);
    v[1] = eit->first->vertex(eit->third);
    
    // Construct cells on the side of v[vi]:
    for (int vi=0; vi<2; vi++) {
      sDel_v = get_anchor_del(Rt_Simplex(v[vi]));
      if (!(sDel_v == sDel_e)) {
	Rt_Cell_circulator ccir, cstart;
	ccir = cstart = regular.incident_cells(*eit);
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
}


// Constructs 2-cells of the mixed complex corresponding to facets
// of the regular triangulation
template < class SkinSurfaceTraits_3, class TriangulatedMixedComplex,
	   class Regular_TDS, class Mixed_complex_observer_>
void
Mixed_complex_triangulator_3<SkinSurfaceTraits_3, TriangulatedMixedComplex,
			     Regular_TDS, Mixed_complex_observer_>::
construct_2_cells() {
  Rt_Simplex sDel_v, sDel_e, sDel_f, sVor_f, sVor_c;
  Tmc_Vertex_handle vh[4]; // Implicit function over vLabels is increasing ...
  Rt_Cell_handle rt_ch;
  int index;
	
  for (Rt_Finite_facets_iterator fit = regular.finite_facets_begin();
       fit != regular.finite_facets_end();
       fit ++) {

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
}


// Constructs 3-cells of the mixed complex corresponding to cells
// of the regular triangulation
template < class SkinSurfaceTraits_3, class TriangulatedMixedComplex,
	   class Regular_TDS, class Mixed_complex_observer_>
void
Mixed_complex_triangulator_3<SkinSurfaceTraits_3, TriangulatedMixedComplex,
			     Regular_TDS, Mixed_complex_observer_>::
construct_3_cells() {
  Rt_Simplex sDel_v, sDel_e, sDel_f, sDel_c, sVor_c;
  Tmc_Vertex_handle vh[4];
  Tmc_Cell_handle ch;

  for (Rt_Finite_cells_iterator cit = regular.finite_cells_begin();
       cit != regular.finite_cells_end();
       cit++) {
    // construct the tetrahedron:
    //   C[ch], C[Facet(ch,fi)], C[Edge(ch,ei,vi)], C[ch->vertex(vi)]
    sDel_c = get_anchor_del(Rt_Simplex(cit));
    sVor_c = get_anchor_vor(Rt_Simplex(cit));
    Rt_Simplex simplex = Rt_Simplex(cit);
    vh[0] = get_vertex(sDel_c, sVor_c); 
    for (int fi=0; fi<4; fi++) {
      sDel_f = get_anchor_del(Rt_Simplex(Rt_Facet(cit, fi)));
      vh[1] = get_vertex(sDel_f, sVor_c);
      if (vh[0] != vh[1]) {
	for (int vi=1; vi<4; vi++) {
	  int index0 = (fi+vi)&3;
	  sDel_v = get_anchor_del(Rt_Simplex(cit->vertex(index0)));
	  for (int ei=1; ei<4; ei++) {
	    int index1 = (fi+ei)&3;
	    if (vi != ei) {
	      sDel_e = get_anchor_del(Rt_Simplex(Rt_Edge(cit, index0, index1)));
	      vh[2] = get_vertex(sDel_e, sVor_c);
	      // index_13[cit].V[edge_index[index0][index1]];
	      vh[3] = get_vertex(sDel_v, sVor_c);
	      // index_03[cit].V[index0];
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
}

// Adds a vertex to the simplicial complex
template < class SkinSurfaceTraits_3, class TriangulatedMixedComplex,
	   class Regular_TDS, class Mixed_complex_observer_>
typename Mixed_complex_triangulator_3<SkinSurfaceTraits_3, TriangulatedMixedComplex,
				      Regular_TDS, Mixed_complex_observer_>::Tmc_Vertex_handle
Mixed_complex_triangulator_3<SkinSurfaceTraits_3, TriangulatedMixedComplex,
			     Regular_TDS, Mixed_complex_observer_>::
add_vertex (Symb_anchor const &anchor)
{
  Tmc_Vertex_handle vh;
  vh = triangulation_incr_builder.add_vertex();
  vh->point() = get_anchor(anchor.first, anchor.second);
  observer.after_vertex_insertion(anchor.first, anchor.second, vh); 

  return vh;
}

// Gets a vertex from the simplicial complex based on the anchors
template < class SkinSurfaceTraits_3, class TriangulatedMixedComplex,
	   class Regular_TDS, class Mixed_complex_observer_>
typename Mixed_complex_triangulator_3<SkinSurfaceTraits_3, TriangulatedMixedComplex,
				      Regular_TDS, Mixed_complex_observer_>::Tmc_Vertex_handle
Mixed_complex_triangulator_3<SkinSurfaceTraits_3, TriangulatedMixedComplex,
			     Regular_TDS, Mixed_complex_observer_>::get_vertex (
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
template < class SkinSurfaceTraits_3, class TriangulatedMixedComplex,
	   class Regular_TDS, class Mixed_complex_observer_>
typename Mixed_complex_triangulator_3<SkinSurfaceTraits_3, TriangulatedMixedComplex,
				      Regular_TDS, Mixed_complex_observer_>::Tmc_Cell_handle
Mixed_complex_triangulator_3<SkinSurfaceTraits_3, TriangulatedMixedComplex,
			     Regular_TDS, Mixed_complex_observer_>::
add_cell(Tmc_Vertex_handle vh[], int orient, Rt_Simplex s) {
  assert((orient==0) || (orient==1));
  assert(vh[0] != Tmc_Vertex_handle()); assert(vh[1] != Tmc_Vertex_handle());
  assert(vh[2] != Tmc_Vertex_handle()); assert(vh[3] != Tmc_Vertex_handle());
  assert(vh[1] != vh[2]); assert(vh[1] != vh[3]); assert(vh[1] != vh[4]);
  assert(vh[2] != vh[3]); assert(vh[2] != vh[4]); assert(vh[3] != vh[4]);
  Tmc_Cell_handle ch;

  if (orient) {
    CGAL_assertion(orientation(
		     vh[0]->point(), vh[1]->point(),
		     vh[2]->point(), vh[3]->point()) == POSITIVE);
    ch = triangulation_incr_builder.add_cell(vh[0], vh[1], vh[2], vh[3]);
  } else {
    CGAL_assertion(orientation(
		     vh[0]->point(), vh[1]->point(),
		     vh[3]->point(), vh[2]->point()) == POSITIVE);
    ch = triangulation_incr_builder.add_cell(vh[0], vh[1], vh[3], vh[2]);
  }
  observer.after_cell_insertion(s, ch);
  return ch;
}

template < class SkinSurfaceTraits_3, class TriangulatedMixedComplex,
	   class Regular_TDS, class Mixed_complex_observer_>
typename SkinSurfaceTraits_3::Triangulated_mixed_complex_traits::Point_3
Mixed_complex_triangulator_3<SkinSurfaceTraits_3, TriangulatedMixedComplex,
			     Regular_TDS, Mixed_complex_observer_>::
get_orthocenter(Rt_Simplex const &s) {
  Rt_Vertex_handle vh;
  Rt_Edge           e;
  Rt_Facet          f;
  Rt_Cell_handle   ch;
  
  Tmc_Point result;
  switch (s.dimension()) {
    case 0:
      vh=s;
      result = skin_traits.r2t_converter_object()(vh->point());
      break;
    case 1:
      e=s;
      result = orthocenter_obj(
	skin_traits.r2t_converter_object()(e.first->vertex(e.second)->point()),
	skin_traits.r2t_converter_object()(e.first->vertex(e.third)->point()));
      break;
    case 2:
      f=s;
      result = orthocenter_obj(
	skin_traits.r2t_converter_object()(
	  f.first->vertex((f.second+1)&3)->point()),
	skin_traits.r2t_converter_object()(
	  f.first->vertex((f.second+2)&3)->point()),
	skin_traits.r2t_converter_object()(
	  f.first->vertex((f.second+3)&3)->point()));
      break;
    case 3:
      ch=s;
      result = orthocenter_obj(
	skin_traits.r2t_converter_object()(ch->vertex(0)->point()),
	skin_traits.r2t_converter_object()(ch->vertex(1)->point()),
	skin_traits.r2t_converter_object()(ch->vertex(2)->point()),
	skin_traits.r2t_converter_object()(ch->vertex(3)->point()));
      break;
  }
  return result;
}

template < class SkinSurfaceTraits_3, class TriangulatedMixedComplex,
	   class Regular_TDS, class Mixed_complex_observer_>
typename SkinSurfaceTraits_3::Triangulated_mixed_complex_traits::Point_3
Mixed_complex_triangulator_3<SkinSurfaceTraits_3, TriangulatedMixedComplex,
			     Regular_TDS, Mixed_complex_observer_>::
get_anchor(Rt_Simplex const &sDel, Rt_Simplex const &sVor)
{
  Tmc_Point dfoc = get_orthocenter(sDel);
  Tmc_Point vfoc = get_orthocenter(sVor);
	
  return skin_traits.construct_anchor_point_3_object()(dfoc, vfoc);
}

template < 
  class SkinSurfaceTraits_3,
  class Regular_3,
  class TriangulatedMixedComplex_3,
  class MixedComplexObserver_3 >
void triangulate_mixed_complex_3(
  Regular_3 &rt, TriangulatedMixedComplex_3 &sc,
  SkinSurfaceTraits_3 &traits,
  MixedComplexObserver_3 &observer) {

  typedef Mixed_complex_triangulator_3<
    SkinSurfaceTraits_3,
    Regular_3,
    TriangulatedMixedComplex_3,
    MixedComplexObserver_3>              Mixed_complex_triangulator;
    
  Mixed_complex_triangulator(rt, sc, traits, observer);
}


template < 
  class SkinSurfaceTraits_3,
  class Regular_3,
  class TriangulatedMixedComplex_3>
void triangulate_mixed_complex_3(
  Regular_3 &regular, TriangulatedMixedComplex_3 &triangulated_mixed_complex,
  SkinSurfaceTraits_3 &traits) {
  // NGHK: Make the shrink factor the right type:
  Triangulated_mixed_complex_observer_3<
    SkinSurfaceTraits_3, TriangulatedMixedComplex_3, Regular_3>
    observer(CGAL::to_double(traits.shrink_factor()));
  triangulate_mixed_complex_3(
    regular, triangulated_mixed_complex, traits, observer);
}

CGAL_END_NAMESPACE

#endif // TRIANGULATE_MIXED_COMPLEX_H
