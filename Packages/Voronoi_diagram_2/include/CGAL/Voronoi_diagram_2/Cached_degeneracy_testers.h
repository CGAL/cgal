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

#include <CGAL/Voronoi_diagram_2/Edge_less.h>
#include <map>

#include <CGAL/Unique_hash_map.h>
#include <CGAL/edge_list.h>


CGAL_BEGIN_NAMESPACE

CGAL_VORONOI_DIAGRAM_2_BEGIN_NAMESPACE

//=========================================================================
//=========================================================================

template<class Edge_tester_t, class Use_std_map = Tag_false>
class Cached_edge_degeneracy_tester;

template<class Face_degeneracy_t, class Use_std_tag = Tag_false>
class Cached_face_degeneracy_tester;

//=========================================================================
//=========================================================================

template<class Edge_tester_t>
class Cached_edge_degeneracy_tester<Edge_tester_t,Tag_true>
{
private:
  typedef Cached_edge_degeneracy_tester<Edge_tester_t,Tag_true>  Self;

public:
  typedef Edge_tester_t    Edge_degeneracy_tester;

  typedef typename Edge_degeneracy_tester::Delaunay_graph  Delaunay_graph;
  typedef typename Edge_degeneracy_tester::Edge            Edge;
  typedef typename Edge_degeneracy_tester::Face_handle     Face_handle;
  typedef typename Edge_degeneracy_tester::Edge_circulator Edge_circulator;

  typedef typename Edge_degeneracy_tester::All_edges_iterator
  All_edges_iterator;

  typedef typename Edge_degeneracy_tester::Finite_edges_iterator
  Finite_edges_iterator;

  typedef typename Edge_degeneracy_tester::result_type  result_type;
  typedef typename Edge_degeneracy_tester::Arity        Arity;

private:
  typedef std::map<Edge,bool,Edge_less<Edge> >   Edge_map;

  Edge opposite(const Delaunay_graph& dual, const Edge& e) const {
    int i_mirror = dual.tds().mirror_index(e.first, e.second);
    return Edge( e.first->neighbor(e.second), i_mirror );
  }

public:
  bool operator()(const Delaunay_graph& dual, const Edge& e) const {
    if ( dual.dimension() < 2 ) { return false; }

    typename Edge_map::iterator it = emap.find(e);
    if ( it != emap.end() ) { return it->second; }

    bool b = e_tester(dual, e);

    Edge e_opp = opposite(dual, e);
    std::pair<Edge,bool> p = std::make_pair(e,b);
    emap.insert(p);
    p.first = e_opp;
    emap.insert(p);

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
    typename Edge_map::iterator it = emap.find(e);
    if ( it == emap.end() ) { return false; }

    // erase the edge from the map
    emap.erase(it);
    return true;
  }

  void clear() {
    emap.clear();
  }

  void swap(Self& other) {
    e_tester.swap(other.e_tester);
    std::swap(emap, other.emap);
  }

  bool is_valid() const { return true; }

  bool is_valid(const Delaunay_graph& dual) const {
    bool valid = true;
    typename Edge_map::iterator it;
    for (it = emap.begin(); it != emap.end(); ++it) {
      valid = valid && dual.tds().is_edge(it->first.first,
					  it->first.second);
    }
    return valid;
  }

private:
  Edge_degeneracy_tester e_tester;
  mutable Edge_map emap;
};

//=========================================================================

template<class Edge_tester_t>
class Cached_edge_degeneracy_tester<Edge_tester_t,Tag_false>
{
private:
  typedef Cached_edge_degeneracy_tester<Edge_tester_t,Tag_false>  Self;

public:
  typedef Edge_tester_t    Edge_degeneracy_tester;

  typedef typename Edge_degeneracy_tester::Delaunay_graph  Delaunay_graph;
  typedef typename Edge_degeneracy_tester::Edge            Edge;
  typedef typename Edge_degeneracy_tester::Face_handle     Face_handle;
  typedef typename Edge_degeneracy_tester::Edge_circulator Edge_circulator;

  typedef typename Edge_degeneracy_tester::All_edges_iterator
  All_edges_iterator;

  typedef typename Edge_degeneracy_tester::Finite_edges_iterator
  Finite_edges_iterator;

  typedef typename Edge_degeneracy_tester::result_type  result_type;
  typedef typename Edge_degeneracy_tester::Arity        Arity;

private:
  enum Three_valued { UNDEFINED = -1, False, True };

  typedef Unique_hash_map<Edge,Three_valued,CGALi::Edge_hash_function>
  Edge_map;

  Edge opposite(const Delaunay_graph& dual, const Edge& e) const {
    int i_mirror = dual.tds().mirror_index(e.first, e.second);
    return Edge( e.first->neighbor(e.second), i_mirror );
  }

public:
  bool operator()(const Delaunay_graph& dual, const Edge& e) const {
    if ( dual.dimension() < 2 ) { return false; }
    if ( emap.is_defined(e) && emap[e] != UNDEFINED ) {
      return emap[e];
    }

    bool b = e_tester(dual, e);
    Three_valued b3 = (b ? True : False);
    emap[e] = b3;
    emap[opposite(dual, e)] = b3;

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
    if ( emap.is_defined(e) ) { emap[e] = UNDEFINED; }
    return true;
  }

  void clear() {
    emap.clear(UNDEFINED);
  }

  void swap(Self& other) {
    e_tester.swap(other.e_tester);
    std::swap(emap, other.emap);
  }

  bool is_valid() const { return true; }

  bool is_valid(const Delaunay_graph& dual) const {
    bool valid = true;
    All_edges_iterator eit;
    for (eit = dual.all_edges_begin(); eit != dual.all_edges_end(); ++eit) {
      Edge e = *eit;
      bool b = !emap.is_defined(e) || (emap[e] != UNDEFINED);
      valid = valid && b;
    }
    return valid;
  }

private:
  Edge_degeneracy_tester e_tester;
  mutable Edge_map emap;
};


//=========================================================================
//=========================================================================


template<class Face_degeneracy_t>
class Cached_face_degeneracy_tester<Face_degeneracy_t,Tag_true>
{
  // tests whether a face has zero area
public:
  typedef Face_degeneracy_t                     Face_degeneracy_tester;

  typedef typename Face_degeneracy_tester::Delaunay_graph  Delaunay_graph;
  typedef typename Face_degeneracy_tester::Vertex_handle   Vertex_handle;

  typedef typename Face_degeneracy_tester::result_type  result_type;
  typedef typename Face_degeneracy_tester::Arity        Arity;

private:
  typedef Cached_face_degeneracy_tester<Face_degeneracy_tester,Tag_true>  Self;

  typedef std::map<Vertex_handle,bool>  Vertex_map;

public:
  bool operator()(const Delaunay_graph& dual, const Vertex_handle& v) const {
    if ( dual.dimension() < 2 ) { return false; }

    typename Vertex_map::iterator it = vmap.find(v);
    if ( it != vmap.end() ) { return it->second; }

    bool b = f_tester(dual, v);
    vmap.insert( std::make_pair(v,b) );
    return b;
  }

  bool erase(const Vertex_handle& v) const {
    typename Vertex_map::iterator it = vmap.find(v);
    if ( it == vmap.end() ) { return false; }

    // erase the edge from the map
    vmap.erase(it);
    return true;
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
    bool valid = true;
    typename Vertex_map::iterator it;
    for (it = vmap.begin(); it != vmap.end(); ++it) {
      valid = valid && dual.tds().is_vertex(it->first);
    }
    return valid;
  }

private:
  Face_degeneracy_tester f_tester;
  mutable Vertex_map vmap;
};

//=========================================================================

template<class Face_degeneracy_t>
class Cached_face_degeneracy_tester<Face_degeneracy_t,Tag_false>
{
  // tests whether a face has zero area
public:
  typedef Face_degeneracy_t                     Face_degeneracy_tester;

  typedef typename Face_degeneracy_tester::Delaunay_graph  Delaunay_graph;
  typedef typename Face_degeneracy_tester::Vertex_handle   Vertex_handle;

  typedef typename Face_degeneracy_tester::result_type  result_type;
  typedef typename Face_degeneracy_tester::Arity        Arity;

private:
  typedef Cached_face_degeneracy_tester<Face_degeneracy_tester,Tag_false> Self;

  enum Three_valued { UNDEFINED = -1, False, True };
  typedef Unique_hash_map<Vertex_handle,Three_valued>   Vertex_map;

public:
  bool operator()(const Delaunay_graph& dual, const Vertex_handle& v) const {
    if ( dual.dimension() < 2 ) { return false; }
    if ( vmap.is_defined(v) && vmap[v] != UNDEFINED ) { return vmap[v]; }

    bool b = f_tester(dual, v);
    Three_valued b3 = (b ? True : False);
    vmap[v] = b3;
    return b;
  }

  bool erase(const Vertex_handle& v) const {
    if ( vmap.is_defined(v) ) { vmap[v] = UNDEFINED; }
    return true;
  }


  void clear() {
    vmap.clear(UNDEFINED);
  }

  void swap(Self& other) {
    f_tester.swap(other.f_tester);
    std::swap(vmap, other.vmap);
  }

  bool is_valid() const { return true; }

  bool is_valid(const Delaunay_graph& dual) const {
    bool valid = true;
    typename Delaunay_graph::All_vertices_iterator vit;
    for (vit = dual.all_vertices_begin();
	 vit != dual.all_vertices_end(); ++vit) {
      bool b = !vmap.is_defined(vit) || (vmap[vit] != UNDEFINED);
      valid = valid && b;
    }
    return valid;
  }

private:
  Face_degeneracy_tester f_tester;
  mutable Vertex_map vmap;
};

//=========================================================================
//=========================================================================

// Specialization for the identity face degeneracy tester

template<class DG> class Identity_face_degeneracy_tester;

template<class DG>
class Cached_face_degeneracy_tester<Identity_face_degeneracy_tester<DG>,
				    Tag_false>
  : public Identity_face_degeneracy_tester<DG>
{
 private:
  typedef Identity_face_degeneracy_tester<DG>            Base;
  typedef Cached_face_degeneracy_tester<Base,Tag_false>  Self;

 public:
  bool erase(const typename Base::Vertex_handle&) const { return true; }

  bool is_valid() const { return true; }

  bool is_valid(const typename Base::Delaunay_graph&) const { return true; }
};

template<class DG>
class Cached_face_degeneracy_tester<Identity_face_degeneracy_tester<DG>,
				    Tag_true>
  : public Cached_face_degeneracy_tester<Identity_face_degeneracy_tester<DG>,
					 Tag_false>
{};

//=========================================================================
//=========================================================================

CGAL_VORONOI_DIAGRAM_2_END_NAMESPACE

CGAL_END_NAMESPACE

#endif // CGAL_VORONOI_DIAGRAM_2_CACHED_DEGENERACY_TESTERS_H
