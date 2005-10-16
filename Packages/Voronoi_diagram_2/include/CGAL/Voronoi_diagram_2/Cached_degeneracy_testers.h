// Copyright (c) 2005 Foundation for Research and Technology-Hellas (Greece).
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
// Author(s)     : Menelaos Karavelas <mkaravel@tem.uoc.gr>

#ifndef CGAL_VORONOI_DIAGRAM_2_CACHED_DEGENERACY_TESTERS_H
#define CGAL_VORONOI_DIAGRAM_2_CACHED_DEGENERACY_TESTERS_H 1

#include <CGAL/Voronoi_diagram_2/basic.h>
#ifdef USE_STD_MAP
#  include <map>
#else
#  include <CGAL/Unique_hash_map.h>
#endif
#include <CGAL/edge_list.h>
#include <cstdlib>


CGAL_BEGIN_NAMESPACE

CGAL_VORONOI_DIAGRAM_2_BEGIN_NAMESPACE

//=========================================================================
//=========================================================================

#ifdef USE_STD_MAP
template<class Edge_t>
struct Edge_less
{
  typedef Edge_t        Edge;
  typedef bool          result_type;
  typedef Arity_tag<2>  Arity;

  bool operator()(const Edge& e1, const Edge& e2) const {
    if ( e1.first != e2.first ) { return e1.first < e2.first; }
    return e1.second < e2.second;
  }
};
#endif

//=========================================================================


template<class Edge_tester_t>
class Cached_edge_degeneracy_tester
{
 public:
  typedef Edge_tester_t    Edge_degeneracy_tester;
  typedef typename Edge_degeneracy_tester::Delaunay_graph Delaunay_graph;
  typedef typename Edge_degeneracy_tester::Edge Edge;
  typedef typename Edge_degeneracy_tester::Face_handle Face_handle;
  typedef typename Edge_degeneracy_tester::Edge_circulator Edge_circulator;

  typedef typename Edge_degeneracy_tester::All_edges_iterator
  All_edges_iterator;

  typedef typename Edge_degeneracy_tester::Finite_edges_iterator
  Finite_edges_iterator;

  typedef typename Edge_degeneracy_tester::result_type  result_type;
  typedef typename Edge_degeneracy_tester::Arity        Arity;

 private:
  typedef Cached_edge_degeneracy_tester<Edge_degeneracy_tester> Self;

  // true if degenerate, false otherwise
#ifdef USE_STD_MAP
  typedef std::map<Edge,bool,Edge_less<Edge> > Edge_map;
#else
  enum Three_valued { UNDEFINED = -1, False, True };
  typedef Unique_hash_map<Edge,Three_valued,CGALi::Edge_hash_function>
  Edge_map;
#endif

  Edge opposite(const Delaunay_graph& dual, const Edge& e) const {
    int i_mirror = dual.tds().mirror_index(e.first, e.second);
    return Edge( e.first->neighbor(e.second), i_mirror );
  }

 public:
  bool operator()(const Delaunay_graph& dual, const Edge& e) const {
    if ( dual.dimension() < 2 ) { return false; }
#ifdef USE_STD_MAP
    typename Edge_map::iterator it = emap.find(e);
    if ( it != emap.end() ) { return it->second; }
#else
    if ( emap.is_defined(e) && emap[e] != UNDEFINED ) { return emap[e]; }
#endif
    bool b = e_tester(dual, e);
#ifdef USE_STD_MAP
    Edge e_opp = opposite(dual, e);
    std::pair<Edge,bool> p = std::make_pair(e,b);
    emap.insert(p);
    p.first = e_opp;
    emap.insert(p);
#else
    Three_valued b3 = (b ? True : False);
    emap[e] = b3;
    emap[opposite(dual, e)] = b3;
#endif
    return b;
  }

  bool operator()(const Delaunay_graph& dual,
		  const Face_handle& f,	int i) const {
    return operator()(dual, Edge(f,i));
  }

  bool operator()(const Delaunay_graph& dual,
		  const Edge_circulator& ec) const {
    return operator()(dual, *ec);
  }

  bool operator()(const Delaunay_graph& dual,
		  const All_edges_iterator& eit) const {
    return operator()(dual, *eit);
  }

  bool operator()(const Delaunay_graph& dual,
		  const Finite_edges_iterator& eit) const {
    return operator()(dual, *eit);
  }

  bool erase(const Edge& e) const {
#ifdef USE_STD_MAP
    typename Edge_map::iterator it = emap.find(e);
    if ( it == emap.end() ) { return false; }

    // erase the edge from the map
    emap.erase(it);
    return true;
#else
    if ( emap.is_defined(e) ) { emap[e] = UNDEFINED; }
    return true;
#endif
  }

  void clear() {
#ifdef USE_STD_MAP
    emap.clear();
#else
    emap.clear(UNDEFINED);
#endif
  }

  void swap(Self& other) {
    e_tester.swap(other.e_tester);
    std::swap(emap, other.emap);
  }

  bool is_valid() const { return true; }

  bool is_valid(const Delaunay_graph& dual) const {
#ifdef USE_STD_MAP
    bool valid = true;
    typename Edge_map::iterator it;
    for (it = emap.begin(); it != emap.end(); ++it) {
      valid = valid && dual.tds().is_edge(it->first.first,
					  it->first.second);
    }
    return valid;
#else
    bool valid = true;
    typename Delaunay_graph::All_edges_iterator eit;
    for (eit = dual.all_edges_begin(); eit != dual.all_edges_end(); ++eit) {
      Edge e = *eit;
      bool b = !emap.is_defined(e) || (emap[e] != UNDEFINED);
      valid = valid && b;
    }
    return valid;
#endif
  }

 private:
  Edge_degeneracy_tester e_tester;
  mutable Edge_map emap;
};


//=========================================================================
//=========================================================================


template<class Face_degeneracy_t>
class Cached_face_degeneracy_tester
{
  // tests whether a face has zero area
 public:
  typedef Face_degeneracy_t                     Face_degeneracy_tester;

  typedef typename Face_degeneracy_tester::Delaunay_graph  Delaunay_graph;
  typedef typename Face_degeneracy_tester::Vertex_handle   Vertex_handle;

  typedef typename Face_degeneracy_tester::result_type  result_type;
  typedef typename Face_degeneracy_tester::Arity        Arity;

 private:
  typedef Cached_face_degeneracy_tester<Face_degeneracy_tester>   Self;

#ifdef USE_STD_MAP
  typedef std::map<Vertex_handle,bool>  Vertex_map;
#else
  enum Three_valued { UNDEFINED = -1, False, True };
  typedef Unique_hash_map<Vertex_handle,Three_valued>   Vertex_map;
#endif

 public:
  bool operator()(const Delaunay_graph& dual, const Vertex_handle& v) const {
    if ( dual.dimension() < 2 ) { return false; }
#ifdef USE_STD_MAP
    typename Vertex_map::iterator it = vmap.find(v);
    if ( it != vmap.end() ) { return it->second; }
#else
    if ( vmap.is_defined(v) && vmap[v] != UNDEFINED ) { return vmap[v]; }
#endif
    bool b = f_tester(dual, v);
#ifdef USE_STD_MAP
    vmap.insert( std::make_pair(v,b) );
#else
    Three_valued b3 = (b ? True : False);
    vmap[v] = b3;
#endif
    return b;
  }

  bool erase(const Vertex_handle& v) const {
#ifdef USE_STD_MAP
    typename Vertex_map::iterator it = vmap.find(v);
    if ( it == vmap.end() ) { return false; }

    // erase the edge from the map
    vmap.erase(it);
    return true;
#else
    if ( vmap.is_defined(v) ) { vmap[v] = UNDEFINED; }
    return true;
#endif
  }


  void clear() {
    vmap.clear();
  }

  void swap(Self& other) {
    f_tester.swap(other.f_tester);
    std::swap(vmap, other.vmap);
  }

  bool is_valid() const { return true; }

  bool is_valid(const Delaunay_graph& dual) const {
#ifdef USE_STD_MAP
    bool valid = true;
    typename Vertex_map::iterator it;
    for (it = vmap.begin(); it != vmap.end(); ++it) {
      valid = valid && dual.tds().is_vertex(it->first);
    }
    return valid;
#else
    bool valid = true;
    typename Delaunay_graph::All_vertices_iterator vit;
    for (vit = dual.all_vertices_begin();
	 vit != dual.all_vertices_end(); ++vit) {
      bool b = !vmap.is_defined(vit) || (vmap[vit] != UNDEFINED);
      valid = valid && b;
    }
    return valid;
#endif
  }

 private:
  Face_degeneracy_tester f_tester;
  mutable Vertex_map vmap;
};

//=========================================================================

template<class DG> class Identity_face_degeneracy_tester;

template<class DG>
class Cached_face_degeneracy_tester<Identity_face_degeneracy_tester<DG> >
  : public Identity_face_degeneracy_tester<DG>
{
 private:
  typedef Identity_face_degeneracy_tester<DG>   Base;
  typedef Cached_face_degeneracy_tester<Base>   Self;

 public:
  bool erase(const typename Base::Vertex_handle&) const { return true; }

  void clear() {}
  void swap(Self&) {}

  bool is_valid() const { return true; }

  bool is_valid(const typename Base::Delaunay_graph&) const { return true; }
};

//=========================================================================
//=========================================================================


CGAL_VORONOI_DIAGRAM_2_END_NAMESPACE

CGAL_END_NAMESPACE


#endif // CGAL_VORONOI_DIAGRAM_2_CACHED_DEGENERACY_TESTERS_H
