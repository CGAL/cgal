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

#include <CGAL/Compute_anchor_3.h>

#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Triangulated_mixed_complex_observer_3.h>
#include <CGAL/Triangulation_incremental_builder_3.h>

#include <CGAL/Mixed_complex_traits_3.h>

#include <CGAL/iterator.h>

CGAL_BEGIN_NAMESPACE

template < class CMCT >
class CMCT_Cell {
public:
  typedef CMCT_Cell<CMCT>                   Self;
  typedef typename CMCT::Vertex_handle Vertex_handle;

  CMCT_Cell() {
  }
  CMCT_Cell(Vertex_handle vs[]) {
    for (int i=0; i<4; i++) _vs[i] = vs[i];
  }
  CMCT_Cell(const Vertex_handle &vh0,
       const Vertex_handle &vh1, 
       const Vertex_handle &vh2, 
       const Vertex_handle &vh3) 
  {
    _vs[0] = vh0;
    _vs[1] = vh1;
    _vs[2] = vh2;
    _vs[3] = vh3;
  }
  Vertex_handle operator[](int i) const {
    CGAL_assertion((0 <= i) && (i<4));
    return _vs[i];
  }
  Vertex_handle vertex(int i) const {
    CGAL_assertion((0 <= i) && (i<4));
    return _vs[i];
  }
  bool operator==(const Self &other) const {
    return ((_vs[0] == other._vs[0]) &&
	    (_vs[1] == other._vs[1]) &&
	    (_vs[2] == other._vs[2]) &&
	    (_vs[3] == other._vs[3]));
  }
  bool operator!=(const Self &other) const { return !operator==(other); }
private:
  Vertex_handle _vs[4];
};

template < class CMCT >
class CMCT_Cell_iterator {
public:
  typedef CMCT_Cell_iterator<CMCT> Self;
  typedef typename CMCT::Cell Cell;
  typedef typename CMCT::Rt_Finite_vertices_iterator 
    Rt_Finite_vertices_iterator;
  typedef typename CMCT::Rt_Finite_edges_iterator 
    Rt_Finite_edges_iterator;
  typedef typename CMCT::Rt_Finite_facets_iterator 
    Rt_Finite_facets_iterator;
  typedef typename CMCT::Rt_Finite_cells_iterator 
    Rt_Finite_cells_iterator;


  // types of an iterator:
  typedef std::forward_iterator_tag iterator_category;
  typedef Cell                      value_type;
  typedef std::ptrdiff_t            difference_type;
  typedef value_type*               pointer;
  typedef value_type&               reference;

  CMCT_Cell_iterator(const CMCT *cmct_, 
		int d=0) : dim(d), cmct(cmct_) {
    if (cmct->regular.number_of_vertices()==0) {
      dim = 4;
      return;
    }
    vit = cmct->regular.finite_vertices_begin();
    eit = cmct->regular.finite_edges_begin();
    fit = cmct->regular.finite_facets_begin();
    cit = cmct->regular.finite_cells_begin();

    cmct->construct_0_cell(vit, std::back_inserter(cells));
    vit++;
    curr = cells.begin();
  }
  CMCT_Cell_iterator(const CMCT_Cell_iterator &copy) {
    dim=copy.dim; cmct = copy.cmct;
    vit=copy.vit; eit=copy.eit; fit=copy.fit; cit=copy.cit;
    for (typename std::vector<Cell>::const_iterator it = copy.cells.begin();
	 it != copy.cells.end(); it++) {
      cells.push_back(*it);
    }
    CGAL_assertion(std::equal(copy.cells.begin(), copy.cells.end(), cells.begin()));
    int ncells = cells.size()-std::distance
      ((typename std::vector<Cell>::const_iterator)copy.curr, 
       (typename std::vector<Cell>::const_iterator)copy.cells.end());
    CGAL_assertion(ncells >= 0);
    curr = cells.begin();
    while (ncells--!=0) curr++;
    CGAL_assertion
      (std::distance((typename std::vector<Cell>::const_iterator)curr, 
		     (typename std::vector<Cell>::const_iterator)cells.end())
       ==
       std::distance((typename std::vector<Cell>::const_iterator)copy.curr, 
		     (typename std::vector<Cell>::const_iterator)copy.cells.end()));
  }
  Self &operator=(const CMCT_Cell_iterator &other) {
    if (&other != this) {
      dim=other.dim; cmct = other.cmct;
      vit=other.vit; eit=other.eit; fit=other.fit; cit=other.cit;
      for (typename std::vector<Cell>::const_iterator 
	     it = other.cells.begin();
	   it != other.cells.end(); it++) {
	cells.push_back(*it);
      }
      int ncells = cells.size() - std::distance(other.curr, other.cells.end());
      CGAL_assertion(ncells > 0);
      curr = cells.begin();
      while (ncells--!=0) curr++;
      CGAL_assertion(std::distance(curr, (typename std::vector<Cell>::const_iterator) cells.end()) ==
		     std::distance(other.curr, other.cells.end()));
    }
    return *this;
  }

  Self & operator++() { 
    CGAL_assertion(curr != cells.end());
    curr++;
    if (curr != cells.end()) return *this;

    cells.clear();
    while (cells.empty() && (dim < 4)) {
      switch (dim) {
      case 0:
	{
	  if (vit == cmct->regular.finite_vertices_end()) {
	    dim ++;
	  } else {
	    cmct->construct_0_cell(vit, std::back_inserter(cells));
	    vit ++;
	  }	    
	  break;
	}
      case 1:
	{
	  if (eit == cmct->regular.finite_edges_end()) {
	    dim ++;
	  } else {
	    cmct->construct_1_cell(eit, std::back_inserter(cells));
	    eit ++;
	  }	    
	  break;
	}
      case 2:
	{
	  if (fit == cmct->regular.finite_facets_end()) {
	    dim ++;
	  } else {
	    cmct->construct_2_cell(fit, std::back_inserter(cells));
	    fit ++;
	  }	    
	  break;
	}
      case 3:
	{
	  if (cit == cmct->regular.finite_cells_end()) {
	    dim ++;
	  } else {
	    cmct->construct_3_cell(cit, std::back_inserter(cells));
	    cit ++;
	  }	    
	  break;
	}
      };
    }
    curr = cells.begin();
    return *this; 
  }
  CMCT_Cell_iterator & operator++(int i) { 
    CMCT_Cell_iterator temp = *this;
    --*this;
    return temp;
  }

  Cell & operator*() const { 
    CGAL_assertion(dim < 4);
    CGAL_assertion(curr != cells.end());
    return *curr; 
  }
  pointer operator->() const { 
    CGAL_assertion(dim < 4);
    CGAL_assertion(curr != cells.end());
    return &*curr; 
  }


  operator pointer() { return *curr; }

  bool operator==(const CMCT_Cell_iterator &other) const {
    if ((dim == 4) || (other.dim == 4) )
      return (dim==other.dim);
	  
    return (*curr == *other.curr);
  }
  bool operator!=(const CMCT_Cell_iterator &other) const {
    return !(*this==other);
  }
private:
  int dim;
  const CMCT *cmct;
  Rt_Finite_vertices_iterator vit;
  Rt_Finite_edges_iterator eit;
  Rt_Finite_facets_iterator fit;
  Rt_Finite_cells_iterator cit;
  std::vector<Cell> cells; // Cells in the current mixed cell;
  //typename std::vector<Cell>::const_iterator curr;
  typename std::vector<Cell>::iterator curr;
};

template <class MixedComplexTraits_3>
class Combinatorial_mixed_complex_triangulator_3 {
  typedef Combinatorial_mixed_complex_triangulator_3<MixedComplexTraits_3>
  Self;
public:
  typedef MixedComplexTraits_3                   Regular;
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
  typedef const Vertex *                             Vertex_handle;

  typedef CMCT_Cell<Self>                            Cell;
  
  friend class CMCT_Cell_iterator<Self>;

  typedef CMCT_Cell_iterator<Self>                   Cell_iterator;
  

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
  typedef typename Anchor_map::iterator              Anchor_map_iterator;
  typedef typename Anchor_map::const_iterator        Anchor_map_const_iterator;



  typedef std::set<Vertex>                           Vertex_container;
  typedef typename Vertex_container::iterator        Vertex_container_it;
  typedef Vertex_container_it                        Vertex_iterator;  
  
  typedef std::map<Vertex_handle, Vertex_handle>     Symb_vertex_map;
  typedef typename Symb_vertex_map::const_iterator   Symb_vertex_map_const_it;

  class Is_no_vertex {
    const Self *ctmc;
  public:
    Is_no_vertex(const Self *ctmc_) : ctmc(ctmc_) {}
    bool operator()(const Symb_vertex_map_const_it & it) const 
    {
      return (it->first != it->second);
    }
  };

  //typedef Filter_iterator<Symb_vertex_map_const_it, Is_no_vertex> Vertex_iterator;
  

//   class Vertex_iterator {
//     typedef Self Combinatorial_mixed_complex_triangulator_3;
//   public:
//     // types of an iterator:
//     typedef std::forward_iterator_tag iterator_category;
//     typedef Vertex                    value_type;
//     typedef std::ptrdiff_t            difference_type;
//     typedef value_type*               pointer;
//     typedef value_type&               reference;

//     Vertex_iterator(Symb_vertex_map_const_it it, Symb_vertex_map_const_it end) :
//       it(it), end(end) {}
//     Vertex_iterator(Symb_vertex_map_const_it end) :
//       it(end), end(end) {}
    
//     Vertex_iterator& operator++() {
//       do { ++it; } while (it != end && (it->first != it->second));
//       return *this;
//     }
    
//     Vertex_iterator operator++(int) {
//       Vertex_iterator tmp(*this);
//       ++(*this);
//       return tmp;
//     }
    
//     reference operator*() const { return (it->first);  }
//     pointer operator->() const  { return (it->first); }

//     bool operator==(const Vertex_iterator other) {
//       return (other.it == it);
//     }
//     bool operator!=(const Vertex_iterator other) {
//       return !(other.it == it);
//     }
//   private:
//     Symb_vertex_map_const_it it, end;
//   };

public:
  Combinatorial_mixed_complex_triangulator_3(Regular const &regular,
					     bool verbose)
    : regular(regular), verbose(verbose),
      compute_anchor_obj(regular) {
    construct_vertices();
  }

  void construct_vertices();

  template <class OutputIteratorVertices>
  void construct_vertices(OutputIteratorVertices vertices);

  Vertex_iterator vertices_begin() {
    return vertices.begin();
  }
  Vertex_iterator vertices_end() {
    return vertices.end();
  }

  Cell_iterator cells_begin() {
    return Cell_iterator(this);
  }
  Cell_iterator cells_end() {
    return Cell_iterator(this, 4);
  }

  template <class OutputIteratorCells>
  void construct_0_cell(Rt_Vertex_handle rt_vh, 
			OutputIteratorCells out) const;
  template <class OutputIteratorCells>
  void construct_1_cell(const Rt_Edge &e, 
			OutputIteratorCells out) const;
  template <class OutputIteratorCells>
  void construct_1_cell(Rt_Finite_edges_iterator eit, 
			OutputIteratorCells out) const {
    construct_1_cell(*eit, out);
  }
  template <class OutputIteratorCells>
  void construct_2_cell(const Rt_Facet &f, 
			OutputIteratorCells out) const;
  template <class OutputIteratorCells>
  void construct_2_cell(const Rt_Finite_facets_iterator &fit, 
			OutputIteratorCells out) const {
    construct_2_cell(*fit, out);
  }
  template <class OutputIteratorCells>
  void construct_3_cell(Rt_Cell_handle rt_ch, 
			OutputIteratorCells out) const;

  template <class Other_MixedComplexTraits_3> 
  typename Other_MixedComplexTraits_3::Bare_point
  location(const Vertex &v,
	   const Other_MixedComplexTraits_3 &traits) const {
    typename Other_MixedComplexTraits_3::Bare_point p_del = 
      orthocenter(v.first, traits);
    typename Other_MixedComplexTraits_3::Bare_point p_vor = 
      orthocenter(v.second, traits);

    return traits.construct_anchor_point_3_object()(p_del, p_vor);
  }
  template <class Other_MixedComplexTraits_3> 
  typename Other_MixedComplexTraits_3::Bare_point
  location(const Vertex_handle vh,
	   const Other_MixedComplexTraits_3 &traits) const {
    typename Other_MixedComplexTraits_3::Bare_point p_del = 
      orthocenter(vh->first, traits);
    typename Other_MixedComplexTraits_3::Bare_point p_vor = 
      orthocenter(vh->second, traits);

    typename Other_MixedComplexTraits_3::Bare_point result = 
      traits.construct_anchor_point_3_object()(p_del, p_vor);

    return result;
  }

  template <class Other_MixedComplexTraits_3>
  Bounded_side bounded_side(const typename Regular_traits::Bare_point &p,
			    const Cell &c,
			    const Other_MixedComplexTraits_3 &traits) const {
    typedef Other_MixedComplexTraits_3                          Traits;
    typedef typename Traits::Bare_point::R::Tetrahedron_3 Tetrahedron;
    typedef Cartesian_converter<
      typename Regular_traits::Bare_point::R, 
      typename Traits::Bare_point::R>                     Converter;


    typename Traits::Bare_point pts[5];
    for (int i=0; i<4; i++) pts[i] = location(c[i], traits);
    pts[4] = Converter()(p);

    return 
      Tetrahedron(pts[0],pts[1],pts[2],pts[3]).bounded_side(pts[4]);
  }

private:
  template <class Other_MixedComplexTraits_3> 
  typename Other_MixedComplexTraits_3::Bare_point
  orthocenter(const Rt_Simplex &s,
	      const Other_MixedComplexTraits_3 &traits) const {
    Weighted_converter_3
      <Cartesian_converter<typename Regular_traits::Bare_point::R, 
                           typename Other_MixedComplexTraits_3::K> >
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
    return typename Other_MixedComplexTraits_3::Weighted_point();
  }

  Vertex_handle add_vertex(Symb_anchor anchor); 
  template <class OutputIteratorCells>
  void add_cell(Vertex_handle vh[], 
		int orient, 
		OutputIteratorCells cells) const;
	
  Vertex_handle get_vertex(Rt_Simplex &sDel, Rt_Simplex &sVor) const;


  void construct_anchor_del(Rt_Simplex const &sDel);
  void construct_anchor_vor(Rt_Simplex const &sVor);
  void construct_anchors();
  const Rt_Simplex get_anchor_del(Rt_Simplex const &sDel) const {
    return find_anchor(anchor_del2, sDel)->first;
  }
  const Rt_Simplex get_anchor_vor(Rt_Simplex const &sVor) const {
    Anchor_map_const_iterator it = find_anchor(anchor_vor2, sVor);
    return it->first;
  }  
  Anchor_map_const_iterator find_anchor(const Anchor_map &a_map, 
					Rt_Simplex const&s) const {
    Anchor_map_const_iterator it = a_map.find(s);
    return find_anchor(a_map, it);
  }
  Anchor_map_const_iterator find_anchor(const Anchor_map &a_map,
				  Anchor_map_const_iterator &start) const {
    CGAL_assertion(start != a_map.end());
    Anchor_map_const_iterator it = start;
    while (it != it->second) it = it->second;

    return it;
  }
  Anchor_map_iterator find_anchor(Anchor_map &a_map,
				  Anchor_map_iterator &start) {
    // do the union-find-trick:
    CGAL_assertion(start != a_map.end());
    Anchor_map_iterator it1 = start;
    Anchor_map_iterator it2 = it1->second;
    while (it2 != it2->second) {
      it1->second = it2->second;
      // NGHK: changed the type for the map-iterator-hack
      it2->second = it1;
      it2 = it1->second;
    }
    return it2;
  }
  
private:
  Regular const &regular;
  bool verbose;

  Compute_anchor compute_anchor_obj;

  Anchor_map             anchor_del2, anchor_vor2;
  Vertex_container       vertices;
};

template <class MixedComplexTraits_3>
void
Combinatorial_mixed_complex_triangulator_3<MixedComplexTraits_3>::
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

template <class MixedComplexTraits_3>
void
Combinatorial_mixed_complex_triangulator_3<MixedComplexTraits_3>::
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

template <class MixedComplexTraits_3>
void
Combinatorial_mixed_complex_triangulator_3<MixedComplexTraits_3>::
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
template <class MixedComplexTraits_3>
void
Combinatorial_mixed_complex_triangulator_3<MixedComplexTraits_3>::
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
  Rt_Simplex sDel, sVor;

  if (verbose) std::cout << "construct_anchors" << std::endl;
  construct_anchors();

  if (verbose) std::cout << "9 ";
  // anchor dimDel=0, dimVor=3
  for (cit=regular.finite_cells_begin();
       cit!=regular.finite_cells_end(); cit++) {
    sVor = get_anchor_vor(Rt_Simplex(cit));
    for (int i=0; i<4; i++) {
      sDel = get_anchor_del(Rt_Simplex(cit->vertex(i)));
      add_vertex(Symb_anchor(sDel,sVor));
    }
  }
  
  if (verbose) std::cout << "8 ";
  // anchor dimDel=1, dimVor=3
  for (cit = regular.finite_cells_begin(); 
       cit != regular.finite_cells_end(); cit++) {
    sVor = get_anchor_vor(Rt_Simplex(cit));
    for (int i=0; i<3; i++) {
      for (int j=i+1; j<4; j++) {
	sDel = get_anchor_del(Rt_Simplex(Rt_Edge(cit,i,j)));
	add_vertex(Symb_anchor(sDel,sVor));
      }
    }
  }
  
  if (verbose) std::cout << "7 ";
  // anchor dimDel=2, dimVor=3 and dimDel=0, dimVor=2
  for (fit = regular.finite_facets_begin(); 
       fit != regular.finite_facets_end(); fit++) {
    // anchor dimDel=2, dimVor=3
    c1 = fit->first;
    c2 = c1->neighbor(fit->second);
    
    sDel = get_anchor_del(*fit);
    if (!regular.is_infinite(c1)) {
      sVor = get_anchor_vor(c1);
      add_vertex(Symb_anchor(sDel,sVor));
    }
    
    if (!regular.is_infinite(c2)) {
      sVor = get_anchor_vor(c2);
      add_vertex(Symb_anchor(sDel,sVor));
    }
    // anchor dimDel=0, dimVor=2
    sVor = get_anchor_vor(*fit);
    for (int i=1; i<4; i++) {
      sDel = get_anchor_del(Rt_Simplex(c1->vertex((fit->second+i)&3)));
      add_vertex(Symb_anchor(sDel,sVor));
    }
  }
	
  if (verbose) std::cout << "6 ";
  // anchor dimDel=0, dimVor=1
  for (eit=regular.finite_edges_begin(); eit!=regular.finite_edges_end(); eit++) {
    sVor = get_anchor_vor(*eit);
    v1 = eit->first->vertex(eit->second);
    v2 = eit->first->vertex(eit->third);
    sDel = get_anchor_del(v1);
    add_vertex(Symb_anchor(sDel,sVor));
			
    sDel = get_anchor_del(v2);
    add_vertex(Symb_anchor(sDel,sVor));
  }
	
  if (verbose) std::cout << "5 ";
  // anchor dimDel=3, dimVor=3
  for (cit=regular.finite_cells_begin(); 
       cit!=regular.finite_cells_end(); cit++) {
    sDel = get_anchor_del(Rt_Simplex(cit));
    sVor = get_anchor_vor(Rt_Simplex(cit));
    add_vertex(Symb_anchor(sDel,sVor));
  }


  if (verbose) std::cout << "4 ";
  // anchor dimDel=0, dimVor=0
  for (vit=regular.finite_vertices_begin(); 
       vit!=regular.finite_vertices_end(); vit++) {
    sDel = get_anchor_del(Rt_Simplex(vit));
    sVor = get_anchor_vor(Rt_Simplex(vit));
    add_vertex(Symb_anchor(sDel,sVor));
  }
	
  if (verbose) std::cout << "3 ";
  // anchor dimDel=1, dimVor=2
  for (fit = regular.finite_facets_begin(); 
       fit != regular.finite_facets_end(); fit++) {
    c1 = fit->first;
    c2 = c1->neighbor(fit->second);

    sVor = get_anchor_vor(Rt_Simplex(*fit));
    for (int i=1; i<3; i++) {
      for (int j=i+1; j<4; j++) {
	e.first = c1;
	e.second = (fit->second+i)&3;
	e.third = (fit->second+j)&3;
	sDel = get_anchor_del(Rt_Simplex(e));
	add_vertex(Symb_anchor(sDel,sVor));
      }
    }
  }
	
  if (verbose) std::cout << "2 ";
  // anchor dimDel=2, dimVor=2
  for (fit=regular.finite_facets_begin(); 
       fit!=regular.finite_facets_end(); fit++) {
    c1 = fit->first;
    c2 = c1->neighbor(fit->second);

    sVor = get_anchor_vor(Rt_Simplex(*fit));
    sDel = get_anchor_del(Rt_Simplex(*fit));
    add_vertex(Symb_anchor(sDel,sVor));
  }
	
  if (verbose) std::cout << "1" << std::endl;
  // anchor dimDel=1, dimVor=1
  for (eit=regular.finite_edges_begin(); 
       eit!=regular.finite_edges_end(); eit++) {
    v1 = eit->first->vertex(eit->second);
    v2 = eit->first->vertex(eit->third);

    sVor = get_anchor_vor(Rt_Simplex(*eit));
    sDel = get_anchor_del(Rt_Simplex(*eit));
    add_vertex(Symb_anchor(sDel,sVor));
  }
}

// Constructs the vertices of the simplicial complex
template <class MixedComplexTraits_3>
  template <class OutputIteratorVertices>
void
Combinatorial_mixed_complex_triangulator_3<MixedComplexTraits_3>::
construct_vertices(OutputIteratorVertices out)
{
  for (Vertex_iterator it = vertices.begin();
       it != vertices.end(); it++) {
    if (it->first == it->second) *out++ = *it;
  }
}


// Constructs the cells of the mixed complex corresponding
// to Regular vertices
template <class MixedComplexTraits_3>
template <class OutputIteratorCells>
void
Combinatorial_mixed_complex_triangulator_3<MixedComplexTraits_3>::
construct_0_cell(Rt_Vertex_handle rt_vh, OutputIteratorCells cells) const
{
  Rt_Simplex sDel_v, sVor_v, sVor_e, sVor_f, sVor_c;
  Vertex_handle vh[4];
  
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
	    sVor_e = 
	      get_anchor_vor(Rt_Simplex(Rt_Edge(*adj_cell,index,(index+j)&3)));
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
template <class MixedComplexTraits_3>
template <class OutputIteratorCells>
void
Combinatorial_mixed_complex_triangulator_3<MixedComplexTraits_3>::
construct_1_cell(const Rt_Edge &e,
		 OutputIteratorCells cells) const {
  Rt_Simplex sDel_v, sDel_e, sVor_e, sVor_f, sVor_c;
  Vertex_handle vh[4];
  Rt_Vertex_handle v[2];
  
  Rt_Simplex mixed_cell_simplex(e);
  sDel_e = get_anchor_del(Rt_Simplex(e));
  sVor_e = get_anchor_vor(Rt_Simplex(e));
    
  v[0] = e.first->vertex(e.second);
  v[1] = e.first->vertex(e.third);
    
  // Construct cells on the side of v[vi]:
  for (int vi=0; vi<2; vi++) {
    sDel_v = get_anchor_del(Rt_Simplex(v[vi]));
    if (!(sDel_v == sDel_e)) {
      Rt_Cell_circulator ccir, cstart;
      ccir = cstart = regular.incident_cells(e);
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
template <class MixedComplexTraits_3>
template <class OutputIteratorCells>
void
Combinatorial_mixed_complex_triangulator_3<MixedComplexTraits_3>::
construct_2_cell(const Rt_Facet &f,
		 OutputIteratorCells cells) const {
  Rt_Simplex sDel_v, sDel_e, sDel_f, sVor_f, sVor_c;
  Vertex_handle vh[4]; // Implicit function over vLabels is increasing ...
  Rt_Cell_handle rt_ch;
  int index;
	
  rt_ch = f.first;
  index = f.second;
  Rt_Simplex simplex(f);
  sDel_f = get_anchor_del(Rt_Simplex(f));
  sVor_f = get_anchor_vor(Rt_Simplex(f));
		
  for (int i=0; i<2; i++) { // Do this twice
    if (!regular.is_infinite(rt_ch)) {
      sVor_c = get_anchor_vor(Rt_Simplex(rt_ch));
	
      vh[3] = get_vertex(sDel_f, sVor_c);
      Vertex_handle vh2 = get_vertex(sDel_f, sVor_f);
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

  CGAL_assertion(rt_ch == f.first);
  CGAL_assertion(index == f.second);
}


// Constructs 3-cells of the mixed complex corresponding to cells
// of the regular triangulation
template <class MixedComplexTraits_3>
template <class OutputIteratorCells>
void
Combinatorial_mixed_complex_triangulator_3<MixedComplexTraits_3>::
construct_3_cell(Rt_Cell_handle rt_ch,
		 OutputIteratorCells cells) const {
  Rt_Simplex sDel_v, sDel_e, sDel_f, sDel_c, sVor_c;
  Vertex_handle vh[4];

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
template <class MixedComplexTraits_3>
typename Combinatorial_mixed_complex_triangulator_3<MixedComplexTraits_3>::
Vertex_handle
Combinatorial_mixed_complex_triangulator_3<MixedComplexTraits_3>::
add_vertex (Symb_anchor anchor)
{
  Vertex_iterator vit = vertices.find(anchor);
  if (vit == vertices.end()) 
    return &*(vertices.insert(anchor).first);
  return &*vit;
}

// Gets a vertex from the simplicial complex based on the anchors
template <class MixedComplexTraits_3>
typename Combinatorial_mixed_complex_triangulator_3<MixedComplexTraits_3>::
Vertex_handle 
Combinatorial_mixed_complex_triangulator_3<MixedComplexTraits_3>::
get_vertex (Rt_Simplex &sDel, Rt_Simplex &sVor) const
{
  Rt_Simplex sDel2 = get_anchor_del(sDel);
  Rt_Simplex sVor2 = get_anchor_vor(sVor);
  CGAL_assertion(sDel == sDel2);
  CGAL_assertion(sVor == sVor2);
  Vertex_container_it it = vertices.find(Symb_anchor(sDel2,sVor2));
//   Symb_vertex_map_const_it it = anchors.find(Symb_anchor(sDel2,sVor2));
  CGAL_assertion(it != vertices.end());
  Vertex_handle vh = &*it;
  CGAL_assertion(*vh != Vertex());
  return vh;
}

// Adds a cell to the simplicial complex
template <class MixedComplexTraits_3>
template <class OutputIteratorCells>
void
Combinatorial_mixed_complex_triangulator_3<MixedComplexTraits_3>::
add_cell(Vertex_handle vh[], int orient, OutputIteratorCells cells) const {
  assert((orient==0) || (orient==1));
  assert(*vh[0] != Vertex()); assert(*vh[1] != Vertex());
  assert(*vh[2] != Vertex()); assert(*vh[3] != Vertex());
  assert(*vh[0] != *vh[1]); assert(*vh[0] != *vh[2]); assert(*vh[0] != *vh[3]);
  assert(*vh[1] != *vh[2]); assert(*vh[1] != *vh[3]); assert(*vh[2] != *vh[3]);

  if (orient) {
    *cells++ = Cell(vh[0], vh[1], vh[2], vh[3]);
  } else {
    *cells++ = Cell(vh[0], vh[1], vh[3], vh[2]);
  }
}

template <class MixedComplexTraits_3,
	  class OutputTriangulation_3,
	  class TriangulatedMixedComplexObserver_3>
void 
triangulate_mixed_complex_3(MixedComplexTraits_3 &rt,
			    typename MixedComplexTraits_3::Geom_traits::FT
			    const & shrink_factor,
			    OutputTriangulation_3 &tmc,
			    TriangulatedMixedComplexObserver_3 &observer,
			    bool verbose) 
{
  typedef MixedComplexTraits_3            Regular;
  typedef typename Regular::Finite_vertices_iterator Rt_Vertices_iterator;
  typedef typename Regular::Finite_edges_iterator    Rt_Edges_iterator;
  typedef typename Regular::Finite_facets_iterator   Rt_Facets_iterator;
  typedef typename Regular::Finite_cells_iterator    Rt_Cells_iterator;
  typedef Triangulation_simplex_3<Regular>           Rt_Simplex;

  typedef Combinatorial_mixed_complex_triangulator_3<MixedComplexTraits_3>  
    CMCT;
  typedef typename CMCT::Vertex          Cmct_Vertex;
  typedef typename CMCT::Vertex_handle   Cmct_Vertex_handle;
  typedef typename CMCT::Cell            Cmct_Cell;
  typedef typename CMCT::Vertex_iterator Cmct_Vertex_iterator;

  typedef Triangulation_incremental_builder_3<OutputTriangulation_3>
    Triangulation_incremental_builder;
  typedef Mixed_complex_traits_3<typename OutputTriangulation_3::Geom_traits>
    Mc_traits;

  typedef typename OutputTriangulation_3::Vertex_handle Out_Vertex_handle;
  typedef typename OutputTriangulation_3::Cell_handle   Out_Cell_handle;

  CMCT mc_triangulator(rt, verbose);
  Triangulation_incremental_builder triangulation_incr_builder(tmc);

  std::map <Cmct_Vertex_handle, Out_Vertex_handle> vertex_map;

  triangulation_incr_builder.begin_triangulation(3);

  { // Vertices
    if (verbose) std::cout << "Construct vertices" << std::endl;
//     std::list<const Cmct_Vertex_handle> vertices;
    //std::list<Cmct_Vertex_iterator> v2;
    //mc_triangulator.construct_vertices(std::back_inserter(vertices));
//     { // NGHK: DEBUG CODE
//       Cmct_Vertex_iterator vit = mc_triangulator.vertices_begin();
//       CGAL_assertion((int)std::distance(mc_triangulator.vertices_begin(), 
// 					mc_triangulator.vertices_end())==
// 		     (int)vertices.size());
//     }
    for (Cmct_Vertex_iterator vit = mc_triangulator.vertices_begin();
	 vit != mc_triangulator.vertices_end(); vit++) {
      Out_Vertex_handle vh = triangulation_incr_builder.add_vertex();
      vh->point() = mc_triangulator.location(*vit, Mc_traits(shrink_factor));
      Cmct_Vertex_handle cmct_vh = &*vit;
      vertex_map[cmct_vh] = vh;
      observer.after_vertex_insertion((*vit).first, (*vit).second, vh);
    }
  }

  { // Cells
    std::vector<Cmct_Cell> cells;
    int nCells=0, nSimplices=0;

    // mixed cells corresponding to regular vertices
    if (verbose) std::cout << "Construct 0 cells" << std::endl;
    for (Rt_Vertices_iterator vit = rt.finite_vertices_begin();
	 vit != rt.finite_vertices_end(); vit ++) {
      mc_triangulator.construct_0_cell(vit, std::back_inserter(cells));
      nCells += cells.size();nSimplices++;
      Rt_Simplex s(vit);
      for (typename std::vector<Cmct_Cell>::iterator it = cells.begin();
	   it != cells.end(); it++) {
	Out_Cell_handle ch = 
	  triangulation_incr_builder.add_cell(vertex_map[(it->vertex(0))],
					      vertex_map[(it->vertex(1))],
					      vertex_map[(it->vertex(2))],
					      vertex_map[(it->vertex(3))]);      
	observer.after_cell_insertion(s, ch);
      }
      cells.clear();
    }

    // mixed cells corresponding to regular edges
    if (verbose) std::cout << "Construct 1 cells" << std::endl;
    for (Rt_Edges_iterator eit = rt.finite_edges_begin();
	 eit != rt.finite_edges_end(); eit ++) {
      mc_triangulator.construct_1_cell(eit, std::back_inserter(cells));
      nCells += cells.size();nSimplices++;
      Rt_Simplex s(*eit);
      for (typename std::vector<Cmct_Cell>::iterator it = cells.begin();
	   it != cells.end(); it++) {
	Out_Cell_handle ch = 
	  triangulation_incr_builder.add_cell(vertex_map[(it->vertex(0))],
					      vertex_map[(it->vertex(1))],
					      vertex_map[(it->vertex(2))],
					      vertex_map[(it->vertex(3))]);      
	observer.after_cell_insertion(s, ch);
      }
      cells.clear();
    }

    // mixed cells corresponding to regular facets
    if (verbose) std::cout << "Construct 2 cells" << std::endl;
    for (Rt_Facets_iterator fit = rt.finite_facets_begin();
	 fit != rt.finite_facets_end(); fit ++) {
      mc_triangulator.construct_2_cell(fit, std::back_inserter(cells));
      nCells += cells.size();nSimplices++;
      Rt_Simplex s(*fit);
      for (typename std::vector<Cmct_Cell>::iterator it = cells.begin();
	   it != cells.end(); it++) {
	Out_Cell_handle ch = 
	  triangulation_incr_builder.add_cell(vertex_map[(it->vertex(0))],
					      vertex_map[(it->vertex(1))],
					      vertex_map[(it->vertex(2))],
					      vertex_map[(it->vertex(3))]);      
	observer.after_cell_insertion(s, ch);
      }
      cells.clear();
    }

    // mixed cells corresponding to regular cells
    if (verbose) std::cout << "Construct 3 cells" << std::endl;
    for (Rt_Cells_iterator cit = rt.finite_cells_begin();
	 cit != rt.finite_cells_end(); cit ++) {
      mc_triangulator.construct_3_cell(cit, std::back_inserter(cells));
      nCells += cells.size();nSimplices++;
      Rt_Simplex s(cit);
      for (typename std::vector<Cmct_Cell>::iterator it = cells.begin();
	   it != cells.end(); it++) {
	Out_Cell_handle ch = 
	  triangulation_incr_builder.add_cell(vertex_map[(it->vertex(0))],
					      vertex_map[(it->vertex(1))],
					      vertex_map[(it->vertex(2))],
					      vertex_map[(it->vertex(3))]);      
	observer.after_cell_insertion(s, ch);
      }
      cells.clear();
    }
    CGAL_assertion(std::distance(mc_triangulator.cells_begin(), 
				 mc_triangulator.cells_end()) == nCells);
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
  class MixedComplexTraits_3,
  class TriangulatedMixedComplex_3>
void 
triangulate_mixed_complex_3(MixedComplexTraits_3 const &regular, 
			    typename MixedComplexTraits_3::Geom_traits::FT
			    const &shrink_factor,
			    TriangulatedMixedComplex_3 &tmc,
			    bool verbose)
{
  Triangulated_mixed_complex_observer_3<
    TriangulatedMixedComplex_3, const MixedComplexTraits_3> 
    observer(shrink_factor);
  triangulate_mixed_complex_3(regular, shrink_factor, tmc, observer, verbose);
}

CGAL_END_NAMESPACE

#endif // CGAL_TRIANGULATE_MIXED_COMPLEX_H
