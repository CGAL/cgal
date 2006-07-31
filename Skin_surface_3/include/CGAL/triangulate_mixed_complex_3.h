// Copyright (c) 2005 Rijksuniversiteit Groningen (Netherlands)
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
// NGHK: move this one to SkinSurfaceTraits
#include <CGAL/Mixed_complex_traits_3.h>

CGAL_BEGIN_NAMESPACE


template <class RegularTriangulation_3>
class Combinatorial_mixed_complex_triangulator_3 {
public:
  typedef RegularTriangulation_3                   Regular;
  typedef typename Regular::Geom_traits            Regular_traits;

  // NGHK: Make private again: private:
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

  typedef Compute_anchor_3<Regular>                  Compute_anchor;
  typedef std::pair<Rt_Simplex,Rt_Simplex>           Symb_anchor;

  typedef Symb_anchor                                Vertex;
  class Cell {
  public:
    Cell(Symb_anchor vs[]) {
      for (int i=0; i<4; i++) _vs[i] = vs[i];
    }
    Cell(const Symb_anchor &vh0,
	 const Symb_anchor &vh1, 
	 const Symb_anchor &vh2, 
	 const Symb_anchor &vh3) 
    {
      _vs[0] = vh0;
      _vs[1] = vh1;
      _vs[2] = vh2;
      _vs[3] = vh3;
    }
    Vertex operator[](int i) {
      CGAL_assertion((0 <= i) && (i<4));
      return _vs[i];
    }
  private:
    Symb_anchor _vs[4];
  };

  // You might get type differences here:
  // The map that maps a Rt_Simplex to an iterator of the map 
  // (used as union_find_structure)
  struct Anchor_map_iterator_tmp;
  typedef std::map<Rt_Simplex, Anchor_map_iterator_tmp>     Anchor_map;
  struct Anchor_map_iterator_tmp : Anchor_map::iterator {
    Anchor_map_iterator_tmp() 
      : Anchor_map::iterator() {}
    Anchor_map_iterator_tmp(typename Anchor_map::iterator const &it) 
      : Anchor_map::iterator(it) {}
  };
  typedef typename Anchor_map::iterator                     Anchor_map_iterator;

public:
  Combinatorial_mixed_complex_triangulator_3(Regular const &regular,
			       Rt_FT const &shrink, 
			       bool verbose)
    : regular(regular), shrink(shrink), verbose(verbose),
      compute_anchor_obj(regular) {
  }


public:
  template <class OutputIteratorVertices>
  void construct_vertices(OutputIteratorVertices vertices);

  template <class OutputIteratorCells>
  void construct_0_cell(Rt_Vertex_handle rt_vh, 
			OutputIteratorCells out);
  template <class OutputIteratorCells>
  void construct_1_cell(const Rt_Finite_edges_iterator &eit, 
			OutputIteratorCells out);
  template <class OutputIteratorCells>
  void construct_2_cell(const Rt_Finite_facets_iterator &fit, 
			OutputIteratorCells out);
  template <class OutputIteratorCells>
  void construct_3_cell(Rt_Cell_handle rt_ch, 
			OutputIteratorCells out);

  template <class MixedComplexTraits_3> 
  typename MixedComplexTraits_3::Bare_point
  location(const Vertex &v,
	   const MixedComplexTraits_3 &traits) const {
    typename MixedComplexTraits_3::Bare_point p_del = 
      orthocenter(v.first, traits);
    typename MixedComplexTraits_3::Bare_point p_vor = 
      orthocenter(v.second, traits);

    return traits.construct_anchor_point_3_object()(p_del, p_vor);
  }
  void location(const Vertex &v) const {
  }

private:
  template <class MixedComplexTraits_3> 
  typename MixedComplexTraits_3::Bare_point
  orthocenter(const Rt_Simplex &s,
	      const MixedComplexTraits_3 &traits) const {
    Weighted_converter_3
      <Cartesian_converter<typename Regular_traits::Bare_point::R, 
                           typename MixedComplexTraits_3::K> >
      converter;
    switch(s.dimension()) {
    case 0: 
      {
	Rt_Vertex_handle vh = s;
	return converter(vh->point());
      }
    case 1:
      {
	Rt_Edge e = s;
	return traits.construct_weighted_circumcenter_3_object()
	  (converter(e.first->vertex(e.second)->point()),
	   converter(e.first->vertex(e.third)->point()));
      }
    case 2: 
      {
	Rt_Facet f = s;
	return traits.construct_weighted_circumcenter_3_object()
	  (converter(f.first->vertex((f.second+1)&3)->point()),
	   converter(f.first->vertex((f.second+2)&3)->point()),
	   converter(f.first->vertex((f.second+3)&3)->point()));
      }
    case 3: 
      {
	Rt_Cell_handle ch = s;
	return traits.construct_weighted_circumcenter_3_object()
	  (converter(ch->vertex(0)->point()),
	   converter(ch->vertex(1)->point()),
	   converter(ch->vertex(2)->point()),
	   converter(ch->vertex(3)->point()));
      }
    }
    CGAL_assertion(false);
    return typename MixedComplexTraits_3::Weighted_point_3();
  }

  template <class OutputIteratorVertices>
  Vertex add_vertex(Symb_anchor const &anchor, 
			OutputIteratorVertices vertices); 
  template <class OutputIteratorCells>
  void add_cell(Vertex vh[], 
		int orient, 
		OutputIteratorCells cells);
	
  Vertex &get_vertex(Rt_Simplex &sDel, Rt_Simplex &sVor);


  void construct_anchor_del(Rt_Simplex const &sDel);
  void construct_anchor_vor(Rt_Simplex const &sVor);
  void construct_anchors();
  Rt_Simplex get_anchor_del(Rt_Simplex const &sDel) {
    return find_anchor(anchor_del2, sDel)->first;
  }
  Rt_Simplex get_anchor_vor(Rt_Simplex const &sVor) {
    return find_anchor(anchor_vor2, sVor)->first;
  }  
  Anchor_map_iterator find_anchor(Anchor_map &a_map, Rt_Simplex const&s) {
    return find_anchor(a_map, a_map.find(s));
  }
  Anchor_map_iterator find_anchor(Anchor_map &a_map,
				  Anchor_map_iterator const&it) {
    CGAL_assertion(it != a_map.end());
    Anchor_map_iterator it2 = it->second;
    while (it2 != it2->second) {
      it->second = it2->second;
      // NGHK: changed the type for the map-iterator-hack
      it2->second = it;
      it2 = it->second;
    }
    return it2;
  }
  
private:
  Regular const &regular;
  Rt_FT const &shrink;
  bool verbose;

  Compute_anchor_3<Regular> compute_anchor_obj;

  Anchor_map                                     anchor_del2, anchor_vor2;
  std::map<Symb_anchor, Vertex>        anchors;
};

template <class RegularTriangulation_3>
void
Combinatorial_mixed_complex_triangulator_3<RegularTriangulation_3>::
construct_anchor_del(Rt_Simplex const &sDel) {
  Rt_Simplex s = compute_anchor_obj.anchor_del(sDel);
  anchor_del2[sDel] = Anchor_map_iterator();

  Anchor_map_iterator it = anchor_del2.find(sDel);
  Anchor_map_iterator it2 = anchor_del2.find(s);
  CGAL_assertion(it != anchor_del2.end());
  CGAL_assertion(it2 != anchor_del2.end());
  it->second = it2;

  // degenerate simplices:
  if (compute_anchor_obj.is_degenerate()) {
    it = find_anchor(anchor_del2, it);
    typename Compute_anchor::Simplex_iterator degenerate_it;
    for (degenerate_it = compute_anchor_obj.equivalent_anchors_begin();
	 degenerate_it != compute_anchor_obj.equivalent_anchors_end(); 
	 degenerate_it++) {
      Anchor_map_iterator tmp;
      it2 = anchor_del2.find(*degenerate_it);
      CGAL_assertion(it2 != anchor_del2.end());
      // Merge sets:
      while (it2 != it2->second) {
	tmp = it2->second;
	it2->second = it->second;
	it2 = tmp;
	CGAL_assertion(it2 != anchor_del2.end());
      }
      it2->second = it->second;
    }
  }
}

template <class RegularTriangulation_3>
void
Combinatorial_mixed_complex_triangulator_3<RegularTriangulation_3>::
construct_anchor_vor(Rt_Simplex const &sVor) {
  Rt_Simplex s = compute_anchor_obj.anchor_vor(sVor);
  anchor_vor2[sVor] = Anchor_map_iterator();

  Anchor_map_iterator it = anchor_vor2.find(sVor);
  Anchor_map_iterator it2 = anchor_vor2.find(s);
  CGAL_assertion(it != anchor_vor2.end());
  CGAL_assertion(it2 != anchor_vor2.end());
  it->second = it2;

  // degenerate simplices:
  if (compute_anchor_obj.is_degenerate()) {
    it = find_anchor(anchor_vor2, it);
    typename Compute_anchor::Simplex_iterator degenerate_it;
    for (degenerate_it = compute_anchor_obj.equivalent_anchors_begin();
	 degenerate_it != compute_anchor_obj.equivalent_anchors_end(); 
	 degenerate_it++) {
      Anchor_map_iterator tmp;
      it2 = anchor_vor2.find(*degenerate_it);
      // Possibly not found for 2 Voronoi vertices with the same center,
      // If the first vertex is inserted and the second is already found.
      // see compute_anchor_obj.anchor_vor(Cell_handle)
      if (it2 != anchor_vor2.end()) {
	CGAL_assertion(it2 != anchor_vor2.end());
	// Merge sets:
	while (it2 != it2->second) {
	  tmp = it2->second;
	  it2->second = it->second;
	  it2 = tmp;
	  CGAL_assertion(it2 != anchor_vor2.end());
	}
	it2->second = it->second;
      }
    }
  }
}

template <class RegularTriangulation_3>
void
Combinatorial_mixed_complex_triangulator_3<RegularTriangulation_3>::
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
template <class RegularTriangulation_3>
  template <class OutputIteratorVertices>
void
Combinatorial_mixed_complex_triangulator_3<RegularTriangulation_3>::
construct_vertices(OutputIteratorVertices vertices)
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
  Rt_Simplex sDel, sVor;
  Vertex vh;

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
	vh = add_vertex(Symb_anchor(sDel,sVor), vertices);
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
	  vh = add_vertex(Symb_anchor(sDel,sVor), vertices);
	  anchors[Symb_anchor(sDel,sVor)] = vh;
	  assert(vh == get_vertex(sDel, sVor));
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
	vh = add_vertex(Symb_anchor(sDel,sVor), vertices);
	anchors[Symb_anchor(sDel,sVor)] = vh;
	assert(vh == get_vertex(sDel, sVor));
      }
    }
    if (!regular.is_infinite(c2)) {
      sVor = get_anchor_vor(c2);
      if (anchors.find(Symb_anchor(sDel,sVor)) == anchors.end()) {
	vh = add_vertex(Symb_anchor(sDel,sVor), vertices);
	anchors[Symb_anchor(sDel,sVor)] = vh;
	assert(vh == get_vertex(sDel, sVor));
      }
    }
    // anchor dimDel=0, dimVor=2
    sVor = get_anchor_vor(*fit);
    for (int i=1; i<4; i++) {
      sDel = get_anchor_del(Rt_Simplex(c1->vertex((fit->second+i)&3)));
      if (anchors.find(Symb_anchor(sDel,sVor)) == anchors.end()) {
	vh = add_vertex(Symb_anchor(sDel,sVor), vertices);
	anchors[Symb_anchor(sDel,sVor)] = vh;
	assert(vh == get_vertex(sDel, sVor));
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
      vh = add_vertex(Symb_anchor(sDel,sVor), vertices);
      anchors[Symb_anchor(sDel,sVor)] = vh;
      assert(vh == get_vertex(sDel, sVor));
    }
			
    sDel = get_anchor_del(v2);
    if (anchors.find(Symb_anchor(sDel,sVor)) == anchors.end()) {
      vh = add_vertex(Symb_anchor(sDel,sVor), vertices);
      anchors[Symb_anchor(sDel,sVor)] = vh;
      assert(vh == get_vertex(sDel, sVor));
    }
  }
	
  if (verbose) std::cout << "5 ";
  // anchor dimDel=3, dimVor=3
  for (cit=regular.finite_cells_begin(); cit!=regular.finite_cells_end(); cit++) {
    sDel = get_anchor_del(Rt_Simplex(cit));
    sVor = get_anchor_vor(Rt_Simplex(cit));
    if (anchors.find(Symb_anchor(sDel,sVor)) == anchors.end()) {
      vh = add_vertex(Symb_anchor(sDel,sVor), vertices);
      anchors[Symb_anchor(sDel,sVor)] = vh;
      assert(vh == get_vertex(sDel, sVor));
    }
  }


  if (verbose) std::cout << "4 ";
  // anchor dimDel=0, dimVor=0
  for (vit=regular.finite_vertices_begin(); vit!=regular.finite_vertices_end(); vit++) {
    sDel = get_anchor_del(Rt_Simplex(vit));
    sVor = get_anchor_vor(Rt_Simplex(vit));
    if (anchors.find(Symb_anchor(sDel,sVor)) == anchors.end()) {
      vh = add_vertex(Symb_anchor(sDel,sVor), vertices);
      anchors[Symb_anchor(sDel,sVor)] = vh;
      assert(vh == get_vertex(sDel, sVor));
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
	  vh = add_vertex(Symb_anchor(sDel,sVor), vertices);
	  anchors[Symb_anchor(sDel,sVor)] = vh;
	  assert(vh == get_vertex(sDel, sVor));
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
      vh = add_vertex(Symb_anchor(sDel,sVor), vertices);
      anchors[Symb_anchor(sDel,sVor)] = vh;
      assert(vh == get_vertex(sDel, sVor));
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
      vh = add_vertex(Symb_anchor(sDel,sVor), vertices);
      anchors[Symb_anchor(sDel,sVor)] = vh;
      assert(vh == get_vertex(sDel, sVor));
    }
  }
}

// Constructs the cells of the mixed complex corresponding
// to Regular vertices
template <class RegularTriangulation_3>
template <class OutputIteratorCells>
void
Combinatorial_mixed_complex_triangulator_3<RegularTriangulation_3>::
construct_0_cell(Rt_Vertex_handle rt_vh, OutputIteratorCells cells)
{
  Rt_Simplex sDel_v, sVor_v, sVor_e, sVor_f, sVor_c;
  Vertex vh[4];
  
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
	      
	      add_cell(vh,(index + (j==(i%3+1)? 1:0))&1, cells);
	    }
	  }
	}
      }
    }
  }
}

// Constructs 1-cells of the mixed complex corresponding to edges
// of the regular triangulation
template <class RegularTriangulation_3>
template <class OutputIteratorCells>
void
Combinatorial_mixed_complex_triangulator_3<RegularTriangulation_3>::
construct_1_cell(const Rt_Finite_edges_iterator &e,
		 OutputIteratorCells cells) {
  Rt_Simplex sDel_v, sDel_e, sVor_e, sVor_f, sVor_c;
  Vertex vh[4];
  Rt_Vertex_handle v[2];
  
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
		add_cell(vh,orient,cells);
									
		vh[1] = get_vertex(sDel_v, sVor_f);
		// vh: dimension are (01,02,12,13)
		add_cell(vh,1-orient,cells);
									
		vh[2] = get_vertex(sDel_v, sVor_c);
		// vh: dimension are (01,02,03,13)
		add_cell(vh,orient,cells);
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
template <class RegularTriangulation_3>
template <class OutputIteratorCells>
void
Combinatorial_mixed_complex_triangulator_3<RegularTriangulation_3>::
construct_2_cell(const Rt_Finite_facets_iterator &fit,
		 OutputIteratorCells cells) {
  Rt_Simplex sDel_v, sDel_e, sDel_f, sVor_f, sVor_c;
  Vertex vh[4]; // Implicit function over vLabels is increasing ...
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
      Vertex vh2 = get_vertex(sDel_f, sVor_f);
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

		add_cell(vh,orient,cells);
									
		vh[2] = get_vertex(sDel_e, sVor_c);
		add_cell(vh,1-orient,cells);
									
		vh[1] = get_vertex(sDel_v, sVor_c);
		add_cell(vh,orient,cells);
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
template <class RegularTriangulation_3>
template <class OutputIteratorCells>
void
Combinatorial_mixed_complex_triangulator_3<RegularTriangulation_3>::
construct_3_cell(Rt_Cell_handle rt_ch,
		 OutputIteratorCells cells) {
  Rt_Simplex sDel_v, sDel_e, sDel_f, sDel_c, sVor_c;
  Vertex vh[4];

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
	      add_cell(vh, orient, cells);
	    }
	  }
	}
      }
    }
  }
}

// Adds a vertex to the simplicial complex
template <class RegularTriangulation_3>
template <class OutputIteratorVertices>
typename Combinatorial_mixed_complex_triangulator_3<RegularTriangulation_3>::
Vertex
Combinatorial_mixed_complex_triangulator_3<RegularTriangulation_3>::
add_vertex (Symb_anchor const &anchor,
	    OutputIteratorVertices vertices)
{
  *vertices++ = anchor;

  return anchor;
}

// Gets a vertex from the simplicial complex based on the anchors
template <class RegularTriangulation_3>
typename Combinatorial_mixed_complex_triangulator_3<RegularTriangulation_3>::
Vertex &
Combinatorial_mixed_complex_triangulator_3<RegularTriangulation_3>::
get_vertex (Rt_Simplex &sDel, Rt_Simplex &sVor)
{
  Rt_Simplex sDel2 = get_anchor_del(sDel);
  Rt_Simplex sVor2 = get_anchor_vor(sVor);
  CGAL_assertion(sDel == sDel2);
  CGAL_assertion(sVor == sVor2);
  Vertex &vh = anchors[Symb_anchor(sDel2,sVor2)];
  CGAL_assertion(vh != Vertex());
  return vh;
}

// Adds a cell to the simplicial complex
template <class RegularTriangulation_3>
template <class OutputIteratorCells>
void
Combinatorial_mixed_complex_triangulator_3<RegularTriangulation_3>::
add_cell(Vertex vh[], int orient, OutputIteratorCells cells) {
  assert((orient==0) || (orient==1));
  assert(vh[0] != Vertex()); assert(vh[1] != Vertex());
  assert(vh[2] != Vertex()); assert(vh[3] != Vertex());
  assert(vh[0] != vh[1]); assert(vh[0] != vh[2]); assert(vh[0] != vh[3]);
  assert(vh[1] != vh[2]); assert(vh[1] != vh[3]); assert(vh[2] != vh[3]);

  if (orient) {
    *cells++ = Cell(vh[0], vh[1], vh[2], vh[3]);
  } else {
    *cells++ = Cell(vh[0], vh[1], vh[3], vh[2]);
  }
}

template <class RegularTriangulation_3,
	  class OutputTriangulation_3,
	  class TriangulatedMixedComplexObserver_3>
void 
triangulate_mixed_complex_3(RegularTriangulation_3 &rt,
			    typename RegularTriangulation_3::Geom_traits::FT
			    const & shrink_factor,
			    OutputTriangulation_3 &tmc,
			    TriangulatedMixedComplexObserver_3 &observer,
			    bool verbose) 
{
  typedef RegularTriangulation_3            Regular;
  typedef typename Regular::Finite_vertices_iterator Rt_Vertices_iterator;
  typedef typename Regular::Finite_edges_iterator    Rt_Edges_iterator;
  typedef typename Regular::Finite_facets_iterator   Rt_Facets_iterator;
  typedef typename Regular::Finite_cells_iterator    Rt_Cells_iterator;
  typedef Triangulation_simplex_3<Regular>           Rt_Simplex;

  typedef Combinatorial_mixed_complex_triangulator_3<RegularTriangulation_3>  
    CMCT;
  typedef typename CMCT::Vertex Cmct_Vertex;
  typedef typename CMCT::Cell   Cmct_Cell;
  
  typedef Triangulation_incremental_builder_3<OutputTriangulation_3>
    Triangulation_incremental_builder;
  typedef Mixed_complex_traits_3<typename OutputTriangulation_3::Geom_traits>
    Mc_traits;

  typedef typename OutputTriangulation_3::Vertex_handle Out_Vertex_handle;
  typedef typename OutputTriangulation_3::Cell_handle   Out_Cell_handle;

  CMCT mc_triangulator(rt, shrink_factor, verbose);
  Triangulation_incremental_builder triangulation_incr_builder(tmc);

  std::map <Cmct_Vertex, Out_Vertex_handle> vertex_map;

  triangulation_incr_builder.begin_triangulation(3);

  { // Vertices
    if (verbose) std::cout << "Construct vertices" << std::endl;
    std::list<Cmct_Vertex> vertices;
    mc_triangulator.construct_vertices(std::back_inserter(vertices));
    
    for (typename std::list<Cmct_Vertex>::iterator vit = vertices.begin();
	   vit != vertices.end(); vit++) {
      Out_Vertex_handle vh = triangulation_incr_builder.add_vertex();
      vh->point() = mc_triangulator.location(*vit, Mc_traits(shrink_factor));
      vertex_map[*vit] = vh;
      observer.after_vertex_insertion(vit->first, vit->second, vh);
    }
  }

  { // Cells
    std::vector<Cmct_Cell> cells;

    // mixed cells corresponding to regular vertices
    if (verbose) std::cout << "Construct 0 cells" << std::endl;
    for (Rt_Vertices_iterator vit = rt.finite_vertices_begin();
	 vit != rt.finite_vertices_end(); vit ++) {
      mc_triangulator.construct_0_cell(vit, std::back_inserter(cells));
      Rt_Simplex s(vit);
      for (typename std::vector<Cmct_Cell>::iterator it = cells.begin();
	   it != cells.end(); it++) {
	Out_Cell_handle ch = 
	  triangulation_incr_builder.add_cell(vertex_map[((*it)[0])],
					      vertex_map[((*it)[1])],
					      vertex_map[((*it)[2])],
					      vertex_map[((*it)[3])]);      
	observer.after_cell_insertion(s, ch);
      }
      cells.clear();
    }

    // mixed cells corresponding to regular edges
    if (verbose) std::cout << "Construct 1 cells" << std::endl;
    for (Rt_Edges_iterator eit = rt.finite_edges_begin();
	 eit != rt.finite_edges_end(); eit ++) {
      mc_triangulator.construct_1_cell(eit, std::back_inserter(cells));
      Rt_Simplex s(*eit);
      for (typename std::vector<Cmct_Cell>::iterator it = cells.begin();
	   it != cells.end(); it++) {
	Out_Cell_handle ch = 
	  triangulation_incr_builder.add_cell(vertex_map[((*it)[0])],
					      vertex_map[((*it)[1])],
					      vertex_map[((*it)[2])],
					      vertex_map[((*it)[3])]);      
	observer.after_cell_insertion(s, ch);
      }
      cells.clear();
    }

    // mixed cells corresponding to regular facets
    if (verbose) std::cout << "Construct 2 cells" << std::endl;
    for (Rt_Facets_iterator fit = rt.finite_facets_begin();
	 fit != rt.finite_facets_end(); fit ++) {
      mc_triangulator.construct_2_cell(fit, std::back_inserter(cells));
      Rt_Simplex s(*fit);
      for (typename std::vector<Cmct_Cell>::iterator it = cells.begin();
	   it != cells.end(); it++) {
	Out_Cell_handle ch = 
	  triangulation_incr_builder.add_cell(vertex_map[((*it)[0])],
					      vertex_map[((*it)[1])],
					      vertex_map[((*it)[2])],
					      vertex_map[((*it)[3])]);      
	observer.after_cell_insertion(s, ch);
      }
      cells.clear();
    }

    // mixed cells corresponding to regular cells
    if (verbose) std::cout << "Construct 3 cells" << std::endl;
    for (Rt_Cells_iterator cit = rt.finite_cells_begin();
	 cit != rt.finite_cells_end(); cit ++) {
      mc_triangulator.construct_3_cell(cit, std::back_inserter(cells));
      Rt_Simplex s(cit);
      for (typename std::vector<Cmct_Cell>::iterator it = cells.begin();
	   it != cells.end(); it++) {
	Out_Cell_handle ch = 
	  triangulation_incr_builder.add_cell(vertex_map[((*it)[0])],
					      vertex_map[((*it)[1])],
					      vertex_map[((*it)[2])],
					      vertex_map[((*it)[3])]);      
	observer.after_cell_insertion(s, ch);
      }
      cells.clear();
    }
  }
  
  triangulation_incr_builder.end_triangulation();

//   // mixed cells corresponding to regular edges
//   if (verbose) std::cout << "Construct 1 cells" << std::endl;
//   for (Rt_Finite_edges_iterator eit = regular.finite_edges_begin();
//        eit != regular.finite_edges_end(); eit ++) {
//     construct_1_cell(eit, std::back_inserter(cells));
//   }

//   // mixed cells corresponding to regular facets
//   if (verbose) std::cout << "Construct 2 cells" << std::endl;
//   for (Rt_Finite_facets_iterator fit = regular.finite_facets_begin();
//        fit != regular.finite_facets_end(); fit ++) {
//     construct_2_cell(fit, std::back_inserter(cells));
//   }
    
//   // mixed cells corresponding to regular cells
//   if (verbose) std::cout << "Construct 3 cells" << std::endl;
//   for (Rt_Finite_cells_iterator cit = regular.finite_cells_begin();
//        cit != regular.finite_cells_end();
//        cit++) {
//     construct_3_cell(cit, std::back_inserter(cells));
//   }

//   triangulation_incr_builder.end_triangulation();
    
//   anchors.clear();

//   //remove_small_edges();
//   { // NGHK: debug code:
//     CGAL_assertion(_tmc.is_valid());
//     std::vector<Vertex> ch_vertices;
//     _tmc.incident_vertices(_tmc.infinite_vertex(), 
// 			   std::back_inserter(ch_vertices));
//     for (typename std::vector<Vertex>::iterator
// 	   vit = ch_vertices.begin(); vit != ch_vertices.end(); vit++) {
//       CGAL_assertion((*vit)->sign() == POSITIVE);
//     }
//   }
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

CGAL_END_NAMESPACE

#endif // CGAL_TRIANGULATE_MIXED_COMPLEX_H
