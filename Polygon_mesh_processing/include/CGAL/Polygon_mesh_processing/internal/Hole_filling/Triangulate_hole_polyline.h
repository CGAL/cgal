// Copyright (c) 2015 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Ilker O. Yaz

#ifndef CGAL_HOLE_FILLING_TRIANGULATE_HOLE_POLYLINE_H
#define CGAL_HOLE_FILLING_TRIANGULATE_HOLE_POLYLINE_H

#include <CGAL/license/Polygon_mesh_processing/meshing_hole_filling.h>


#include <CGAL/value_type_traits.h>
#ifndef CGAL_HOLE_FILLING_DO_NOT_USE_DT3
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Delaunay_triangulation_cell_base_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#endif

#ifndef CGAL_HOLE_FILLING_DO_NOT_USE_CDT2
#include <CGAL/centroid.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Projection_traits_3.h>
#include <queue>
#include <CGAL/Polygon_2_algorithms.h>
#endif

#include <CGAL/utility.h>
#include <CGAL/iterator.h>
#include <CGAL/use.h>
#include <CGAL/Kernel/global_functions_3.h>
#include <CGAL/squared_distance_3.h>

#include <vector>
#include <stack>
#include <map>
#include <unordered_set>

#include <CGAL/boost/iterator/transform_iterator.hpp>

namespace CGAL {
namespace internal {

/************************************************************************
 * Lookup tables
 ************************************************************************/
// Wrapper around vector
template<class T>
class Lookup_table {
public:
  Lookup_table(int n, const T& t) : n(n), table(n*n, t) { }
  void put(int i, int j, const T& t) {
    CGAL_assertion(bound_check(i,j));
    table[i*n + j] = t;
  }
  const T& get(int i, int j) const {
    CGAL_assertion(bound_check(i,j));
    return table[i*n + j];
  }

  int n;
private:
  bool bound_check(int i, int j) const {
    CGAL_assertion(i >= 0 && i < n);
    CGAL_assertion(j >= 0 && j < n);
    CGAL_assertion(i < j);
    CGAL_USE(i);
    CGAL_USE(j);
    // previous implementation was based on directly vector and i supposed to be always smaller than j.
    // this check actually can be removed and i =min(i,j) j = max(i,j) can be used for reflexive access
    return true;
  }
  std::vector<T> table;
};

// Wrapper around map, where if i,j is not found a default value is returned,
// and if default value inserted i,j erased.
template<class T>
class Lookup_table_map {
public:
  Lookup_table_map(int n, const T& default_) : n(n), default_(default_) { }

  void put(int i, int j, const T& t) {
    CGAL_assertion(bound_check(i,j));

    if(t == default_) {
      table.erase(std::make_pair(i,j));
      return;
    }

    std::pair<typename Map::iterator, bool> inserted = table.insert(std::make_pair(std::make_pair(i,j), t));
    if(!inserted.second) { inserted.first->second = t;}
  }
  const T& get(int i, int j) const {
    CGAL_assertion(bound_check(i,j));
    typename Map::const_iterator ij = table.find(std::make_pair(i,j));
    if(ij != table.end()) {
      return ij->second;
    }
    return default_;
  }

  void set_range_to_default(int b, int e/*inclusive*/) {
    // given a range b-e erase each entity which falls into b-e
    // example: given range 2-6, entries need to be deleted  2-3, 2-4, 2-5, 2-6
    //                                                            3-4, 3-5, 3-6
    //                                                                 4-5, 4-6
    //                                                                      5-6
    typename Map::iterator it;
    if(b == 0) { it = table.begin(); }
    else       {
      it = table.upper_bound(std::make_pair(b-1, n)); // to find first entry where entry.first == b
    }

    while(it != table.end() && it->first.first != e) {
      if(it->first.second <= e)
      { table.erase(it++); }
      else
      { ++it; }
    }
  }

  int n;
private:
  typedef std::map<std::pair<int,int>, T> Map;
  bool bound_check(int i, int j) const {
    CGAL_assertion(i >= 0 && i < n);
    CGAL_assertion(j >= 0 && j < n);
    CGAL_assertion(i < j);
    CGAL_USE(i);
    CGAL_USE(j);
    return true;
  }
  Map table;
  T default_;
};

/************************************************************************
 * Is_valid classes (to be used in Weight_calculator)
 ************************************************************************/
struct Is_valid_existing_edges
{
  typedef std::vector<std::pair<int, int> > Edge_container;

  Is_valid_existing_edges(Edge_container& existing_edges)
    : existing_edges(existing_edges)
  {
    std::sort(existing_edges.begin(), existing_edges.end());
#ifndef CGAL_NDEBUG
    // all pairs need to be satisfy pair.first < pair.second
    for(Edge_container::iterator it = existing_edges.begin(); it != existing_edges.end(); ++it) {
      CGAL_assertion(it->first < it->second);
    }
#endif
  }

  template<class Point_3>
  bool operator()(const std::vector<Point_3>&,
                  int v0, int v1, int v2) const
  {
    CGAL_assertion(v0 < v1 && v1 < v2);

    if(v0 + 1 != v1 && // border edges can not be inside existing_edges, so no need to check
       std::binary_search(existing_edges.begin(), existing_edges.end(), std::make_pair(v0,v1)) )
    { return false; }

    if(v1 + 1 != v2 &&
      std::binary_search(existing_edges.begin(), existing_edges.end(), std::make_pair(v1,v2)) )
    { return false; }

    if(std::binary_search(existing_edges.begin(), existing_edges.end(), std::make_pair(v2,v0)) )
    { return false; }

    return true;
  }

  Edge_container& existing_edges;
};

struct Is_not_degenerate_triangle
{
  template<class Point_3>
  bool operator()(const std::vector<Point_3>& P,
                  int v0, int v1, int v2) const
  {
    return !CGAL::collinear(P[v0], P[v1], P[v2]);
  }
};

// Combine above two
struct Is_valid_existing_edges_and_degenerate_triangle
{
  Is_valid_existing_edges_and_degenerate_triangle(Is_valid_existing_edges::Edge_container& edges)
    : is_valid_edges(edges) { }

  template<class Point_3>
  bool operator()(const std::vector<Point_3>& P,
                  int v0, int v1, int v2) const
  {
    return Is_not_degenerate_triangle()(P,v0,v1,v2)
                      && is_valid_edges(P,v0,v1,v2);
  }

  Is_valid_existing_edges is_valid_edges;
};

/************************************************************************
 * Weights
 ************************************************************************/

// First minimizes the worst dihedral angle between patch triangles, then the total surface area as a tiebreaker.
class Weight_min_max_dihedral_and_area
{
  template<class Weight_, class IsValid>
  friend struct Weight_calculator;

  template<class Weight_>
  friend class Weight_incomplete;

public:
  // these two should not be used (used in test code)
  std::pair<double,double> w;
  Weight_min_max_dihedral_and_area(double angle, double area) : w(angle, area) { }

// below required by Weight concept
private:
  template<class Point_3, class LookupTable>
  Weight_min_max_dihedral_and_area(const std::vector<Point_3>& P,
                                   const std::vector<Point_3>& Q,
                                   int i, int j, int k,
                                   const LookupTable& lambda)
  {
    CGAL_assertion(i < j);
    CGAL_assertion(j < k);
    int n = static_cast<int>(P.size()) -1; // because the first and last point are equal

    // The CGAL::dihedral angle is measured between the oriented triangles, that is it goes from [-pi, pi]
    // What we need is the angle between the normals of the triangles between [0, pi]
    double ang_max = 0;

    // Test each edge
    int vertices[] = {i, j, k};
    for(int e = 0; e < 3; ++e)
    {
      int v0      = vertices[e];
      int v1      = vertices[(e+1)%3];
      int v_other = vertices[(e+2)%3];
      double angle = 0;
      // check whether the edge is border
      if( (v0 + 1 == v1 || (v0 == n-1 && v1 == 0) ) && !Q.empty() ) {
        angle = 180 - CGAL::abs(
           to_double(CGAL::approximate_dihedral_angle(P[v0],P[v1],P[v_other],Q[v0])) );
      }
      else {
        if(e == 2) { continue; }
        if(lambda.get(v0, v1) != -1){
          const Point_3& p01 = P[lambda.get(v0, v1)];
          angle = 180 - CGAL::abs(
            to_double(CGAL::approximate_dihedral_angle(P[v0],P[v1],P[v_other],p01)) );
        }
      }
      ang_max = (std::max)(ang_max, angle);
    }

    w = std::make_pair(ang_max, to_double(CGAL::approximate_sqrt(CGAL::squared_area(P[i],P[j],P[k]))));
  }

public:
  Weight_min_max_dihedral_and_area operator+(const Weight_min_max_dihedral_and_area& w2) const
  {
    CGAL_assertion((*this) != NOT_VALID());
    CGAL_assertion(w2 != NOT_VALID());

    return Weight_min_max_dihedral_and_area((std::max)(w.first, w2.w.first), w.second + w2.w.second);
  }

  bool operator<(const Weight_min_max_dihedral_and_area& w2) const
  {
    CGAL_assertion((*this) != NOT_VALID());
    CGAL_assertion(w2 != NOT_VALID());

    if(w.first == w2.w.first)
    { return w.second < w2.w.second; }
    return w.first < w2.w.first;
  }

  bool operator==(const Weight_min_max_dihedral_and_area& w2) const
  { return w.first == w2.w.first && w.second == w2.w.second; }

  bool operator!=(const Weight_min_max_dihedral_and_area& w2) const
  { return !(*this == w2); }

  static const Weight_min_max_dihedral_and_area DEFAULT() // rule: x + DEFAULT() == x
  { return Weight_min_max_dihedral_and_area(0,0); }
  static const Weight_min_max_dihedral_and_area NOT_VALID()
  { return Weight_min_max_dihedral_and_area(-1,-1); }

  friend std::ostream& operator<<(std::ostream& out, const Weight_min_max_dihedral_and_area& w) {
    out << "Max dihedral: " << w.w.first << ", Total area: " << w.w.second;
    return out;
  }
};

// For proof of concept. Tested weakly.
class Weight_total_edge {
  template<class Weight_, class IsValid>
  friend struct Weight_calculator;

private:
  double total_length;

  Weight_total_edge(double total_length = 0) : total_length(total_length) { }

  template<class Point_3, class LookupTable>
  Weight_total_edge(const std::vector<Point_3>& P,
                    const std::vector<Point_3>&,
                    int i, int j, int k,
                    const LookupTable&)
    : total_length(0)
  {
    CGAL_assertion(i < j);
    CGAL_assertion(j < k);
    int n = P.size() -1; // because the first and last point are equal

    // Test each edge
    int vertices[] = {i, j, k};
    for(int e = 0; e < 3; ++e)
    {
      int v0      = vertices[e];
      int v1      = vertices[(e+1)%3];

      // check whether the edge is border
      bool border = (v0 + 1 == v1) || (v0 == n-1 && v1 == 0);
      if(!border) {
        total_length += CGAL::sqrt(CGAL::squared_distance(P[v0],P[v1]));
      }
    }
  }

public:
  Weight_total_edge operator+(const Weight_total_edge& w2) const
  {
    CGAL_assertion((*this) != NOT_VALID());
    CGAL_assertion(w2 != NOT_VALID());
    return Weight_total_edge(total_length + w2.total_length);
  }

  bool operator<(const Weight_total_edge& w2) const
  {
    CGAL_assertion((*this) != NOT_VALID());
    CGAL_assertion(w2 != NOT_VALID());
    return total_length < w2.total_length;
  }

  bool operator==(const Weight_total_edge& w2) const
  { return total_length == w2.total_length; }
  bool operator!=(const Weight_total_edge& w2) const
  { return !(*this == w2); }

  static const Weight_total_edge DEFAULT() { return Weight_total_edge(0); } // rule: x + DEFAULT() == x
  static const Weight_total_edge NOT_VALID() { return Weight_total_edge(-1); }
  friend std::ostream& operator<<(std::ostream& out, const Weight_total_edge& w) {
    out << "Total edge length : " << w.total_length;
    return out;
  }
};

// Weights for incomplete patches. It maximize patch size, then ActualWeight as a tiebreaker.
template<class ActualWeight>
class Weight_incomplete
{
  template<class Weight_, class IsValid>
  friend struct Weight_calculator;

private:
  template<class Point_3, class LookupTable>
  Weight_incomplete(const std::vector<Point_3>& P,
                    const std::vector<Point_3>& Q,
                    int i, int j, int k,
                    const LookupTable& lambda)
    : weight(P,Q,i,j,k,lambda), patch_size(1)
  { }

  Weight_incomplete(const ActualWeight& weight, int patch_size)
    : weight(weight), patch_size(patch_size)
  { }

public:
  Weight_incomplete operator+(const Weight_incomplete& w2) const
  {
    CGAL_assertion((*this) != NOT_VALID());
    CGAL_assertion(w2 != NOT_VALID());

    return Weight_incomplete(weight + w2.weight, patch_size + w2.patch_size);
  }

  bool operator<(const Weight_incomplete& w2) const
  {
    CGAL_assertion((*this) != NOT_VALID());
    CGAL_assertion(w2 != NOT_VALID());

    if(patch_size == w2.patch_size) {
      return weight < w2.weight;
    }
    return patch_size > w2.patch_size; // if patch size is larger, then the weight is smaller
  }

  bool operator==(const Weight_incomplete& w2) const
  { return weight == w2.weight && patch_size == w2.patch_size; }

  bool operator!=(const Weight_incomplete& w2) const
  { return !(*this == w2); }

  static const Weight_incomplete DEFAULT() // rule: x + DEFAULT() == x
  { return Weight_incomplete(ActualWeight::DEFAULT(), 0); }
  static const Weight_incomplete NOT_VALID()
  { return Weight_incomplete(ActualWeight::NOT_VALID(), 0); }

  friend std::ostream& operator<<(std::ostream& out, const Weight_incomplete& w) {
    out << "Patch size: " << w.patch_size << ", Actual weight: " << w.weight;
    return out;
  }

  ActualWeight weight;
  int patch_size;
};

// Weight calculator class is both responsible from calculating weights, and checking validity of triangle
template<class Weight_, class IsValid>
struct Weight_calculator
{
  typedef Weight_ Weight;
  Weight_calculator(const IsValid& is_valid = IsValid()) : is_valid(is_valid) { }

  template<class Point_3, class LookupTable>
  Weight operator()(const std::vector<Point_3>& P,
                    const std::vector<Point_3>& Q,
                    int i, int j, int k,
                    const LookupTable& lambda) const
  {
    if( !is_valid(P,i,j,k) )
    { return Weight::NOT_VALID(); }
    return Weight(P, Q, i,j,k, lambda);
  }

  IsValid is_valid;
};
/************************************************************************
 * Tracer
 ************************************************************************/
// It can produce a patch from both complete and incomplete lambda
template<class OutputIteratorValueType, class OutputIteratorPatch, class OutputIteratorHole>
struct Tracer_polyline_incomplete {
  Tracer_polyline_incomplete(OutputIteratorPatch out, OutputIteratorHole out_hole)
    : out(out), out_hole(out_hole)
  { }

  template <class LookupTable>
  void
  operator()(const LookupTable& lambda, int v0, int v1)
  {
    CGAL_assertion_code( const int n = lambda.n; )
    std::stack<std::pair<int, int> > ranges;
    ranges.push(std::make_pair(v0, v1));

    while(!ranges.empty()) {
      std::pair<int, int> r = ranges.top();
      ranges.pop();
      CGAL_assertion(r.first >= 0 && r.first < n);
      CGAL_assertion(r.second >= 0 && r.second < n);

      if(r.first + 1 == r.second) { continue; }

      int la = lambda.get(r.first, r.second);
      if(la == -1) {
        *out_hole++ = std::make_pair(r.first, r.second);
        continue;
      }

      CGAL_assertion(la >= 0 && la < n);
      CGAL_assertion( (r.first < la) && (r.second > la) );
      *out++ = OutputIteratorValueType(r.first, la, r.second);

      ranges.push(std::make_pair(r.first, la));
      ranges.push(std::make_pair(la, r.second));
    }
  }

  OutputIteratorPatch out;
  OutputIteratorHole  out_hole;
};

/************************************************************************
 * Triangulate hole with support of 3D Triangulation
 ************************************************************************/

// to support incident_facets(Edge e) function for both dimension 2 and 3
template<unsigned int Dimension, class Triangulator>
struct Incident_facet_circulator;

template<class Triangulator>
struct Incident_facet_circulator_base
{
  typedef typename Triangulator::Facet         Facet;
  typedef typename Triangulator::Edge          Edge;
  typedef typename Triangulator::Cell_handle   Cell_handle;

  struct Edge_wrapper {
    Edge_wrapper(Edge e) : e(e) { }

    int vertex_first()
    { return (std::min)(e.first->vertex(e.second)->info(), e.first->vertex(e.third)->info()); }
    int vertex_second()
    { return (std::max)(e.first->vertex(e.second)->info(), e.first->vertex(e.third)->info()); }

    Edge e;
  };

  // Finds the other vertex than e.v0 and e.v1 in facet f
  // Note that this may return infinite vertex
  int get_third(Facet f, Edge e) {
    int v0_info = e.first->vertex(e.second)->info();
    int v1_info = e.first->vertex(e.third )->info();
    // warning: it should be designed to handle dimension 2 (e.g. f.first->vertex(3)->info() will crash)
    for(int i = 0; i < 4; ++i) {
      if(i == f.second) { continue; } // skip the vertex which is not on `f`
      int f3 = f.first->vertex(i)->info();
      if(f3 != v0_info && f3 != v1_info) {
        return f3;
      }
    }
    CGAL_assertion(false);
    return -1;
  }

  Edge_wrapper edge_first(Facet f, Edge e) {
    int v0_info = (std::min)(e.first->vertex(e.second)->info(),
                             e.first->vertex(e.third)->info());
    return Edge(f.first,
                get_vertex_index(f.first, v0_info) ,
                get_vertex_index(f.first, get_third(f,e)));
  }

  Edge_wrapper edge_second(Facet f, Edge e) {
    int v1_info = (std::max)(e.first->vertex(e.second)->info(),
                             e.first->vertex(e.third)->info());
    return Edge(f.first,
                get_vertex_index(f.first, v1_info),
                get_vertex_index(f.first, get_third(f,e)));
  }

  int get_vertex_index(Cell_handle ch, int info) {
    // warning: it should be designed to handle dimension 2 (e.g. f.first->vertex(3)->info() will crash)
    for(int i = 0; i < 4; ++i) {
      int v = ch->vertex(i)->info();
      if(v == info) { return i; }
    }
    CGAL_assertion(false);
    return -1;
  }
};

// Use the fact that an edge can be incident to 2 facets in dimension 2
// and all valid facets (which contains finite + infinite vertices but not the NULL vertex) are
// pointed by index 3 in cells
template<class Triangulator>
struct Incident_facet_circulator<2, Triangulator>
  : Incident_facet_circulator_base<Triangulator>
{
  typedef typename Triangulator::Facet         Facet;
  typedef typename Triangulator::Edge          Edge;
  typedef typename Triangulator::Triangulation Triangulation;
  typedef typename Incident_facet_circulator_base<Triangulator>::Edge_wrapper Edge_wrapper;

  Incident_facet_circulator(Edge_wrapper ew, const Triangulation&)
    : f1( Facet(ew.e.first, 3) ),
      f2( Facet(ew.e.first->neighbor(3 - ew.e.second - ew.e.third), 3) ),
      it(f1), e(ew.e)
  {
     CGAL_assertion(f1 != f2);
     CGAL_assertion(e.second < 3 && e.third < 3);
  }
  Incident_facet_circulator& operator++() {
    it = it == f1 ? f2 : f1;
    return *this;
  }
  operator bool() const { return it != f1; }

  int get_third()
  { return Incident_facet_circulator_base<Triangulator>::get_third(it, e); }
  Edge_wrapper edge_first()
  { return Incident_facet_circulator_base<Triangulator>::edge_first(it, e); }
  Edge_wrapper edge_second()
  { return Incident_facet_circulator_base<Triangulator>::edge_second(it, e); }

  Facet f1, f2, it;
  Edge e;
};

// Just a wrapper around Facet_circulator
template<class Triangulator>
struct Incident_facet_circulator<3, Triangulator>
  : Incident_facet_circulator_base<Triangulator>
{
  typedef typename Triangulator::Facet            Facet;
  typedef typename Triangulator::Edge             Edge;
  typedef typename Triangulator::Triangulation    Triangulation;
  typedef typename Triangulator::Facet_circulator Facet_circulator;
  typedef typename Incident_facet_circulator_base<Triangulator>::Edge_wrapper Edge_wrapper;

  Incident_facet_circulator(Edge_wrapper ew, const Triangulation& tr)
    : it(tr.incident_facets(ew.e)), end(it), e(ew.e)
  { }
  Incident_facet_circulator& operator++() {
    ++it;
    return *this;
  }
  operator bool() const { return it != end; }

  int get_third()
  { return Incident_facet_circulator_base<Triangulator>::get_third(*it, e); }
  Edge_wrapper edge_first()
  { return Incident_facet_circulator_base<Triangulator>::edge_first(*it, e); }
  Edge_wrapper edge_second()
  { return Incident_facet_circulator_base<Triangulator>::edge_second(*it, e); }

  Facet_circulator it;
  Facet_circulator end;
  Edge e;
};

// Another DS for search space, which can be used in triangulate_DT
// It is useful for extending the search space of 3D Triangulation by appending new triangles
struct Edge_graph
{
  struct Edge_comp {
    bool operator()(std::pair<int, int> p0, std::pair<int, int> p1) const {
      if(p0.first > p0.second) { std::swap(p0.first, p0.second); }
      if(p1.first > p1.second) { std::swap(p1.first, p1.second); }
      return p0 < p1;
    }
  };

  typedef std::unordered_set<int> Vertex_container;
  // contains edges as key, and each edge contains set of third vertices which denote neighbor facets to that edge
  typedef std::map<std::pair<int, int>, Vertex_container, Edge_comp> Graph;

  struct Edge_wrapper {
    Edge_wrapper(std::pair<int, int> e) : e(e) { }

    int vertex_first()  const { return (std::min)(e.first, e.second); }
    int vertex_second() const { return (std::max)(e.first, e.second); }
    std::pair<int, int> e;
  };

  struct Incident_facet_circulator
  {
    Incident_facet_circulator(Edge_wrapper ew, const Edge_graph& tr)
    {
      Graph::const_iterator it_e = tr.graph.find(ew.e);
      CGAL_assertion(it_e != tr.graph.end());
      it = it_e->second.begin();
      end = it_e->second.end();
      e = ew.e;
    }
    Incident_facet_circulator& operator++() {
      ++it;
      return *this;
    }
    operator bool() const { return it != end; }

    int get_third() const { return *it; }
    Edge_wrapper edge_first()  const { return std::make_pair(Edge_wrapper(e).vertex_first(), *it); }
    Edge_wrapper edge_second() const { return std::make_pair(Edge_wrapper(e).vertex_second(), *it); }

    Vertex_container::const_iterator it;
    Vertex_container::const_iterator end;
    std::pair<int, int> e;
  };

  template<class IncidentFacetCirculator, class Triangulation>
  void init(const Triangulation& tr, const std::vector<bool>& edge_exist)
  {
    typedef typename Triangulation::Finite_edges_iterator Finite_edges_iterator;

    n = static_cast<int>(edge_exist.size());
    for(Finite_edges_iterator eb = tr.finite_edges_begin(); eb != tr.finite_edges_end(); ++eb)
    {
      int v0 = eb->first->vertex(eb->second)->info();
      int v1 = eb->first->vertex(eb->third )->info();
      Vertex_container& e_neighs = graph[std::make_pair(v0, v1)];

      IncidentFacetCirculator fb(*eb, tr);
      do {
        int v2 = fb.get_third();
        if(v2 == -1) { continue; }
        e_neighs.insert(v2);
      } while(++fb);
    }

    for(int i = 0; i < n; ++i) {
      if(edge_exist[i]) { continue; }
      int v0 = i == n-1 ? 0   : i;
      int v1 = i == n-1 ? n-1 : i+1;
      add_all_possible_to_edge(std::make_pair(v0, v1));
    }
  }

  void add_all_possible_to_edge(std::pair<int, int> e) {
    Vertex_container& e_neighs = graph[e];
    for(int i = 0; i < n; ++i) {
      if(i == e.first || i == e.second) { continue; }
      e_neighs.insert(i);
      graph[std::make_pair(i, e.first)].insert(e.second);
      graph[std::make_pair(i, e.second)].insert(e.first);
    }
  }

  Graph graph;
  int n;
};

template<
  class Kernel,
  class Tracer,
  class WeightCalculator,
  class Visitor,
  template <class> class LookupTable = Lookup_table
>
class Triangulate_hole_polyline;

#ifndef CGAL_HOLE_FILLING_DO_NOT_USE_DT3
// By default Lookup_table_map is used, since Lookup_table requires n*n mem.
// Performance decrease is nearly 2x (for n = 10,000, for larger n Lookup_table just goes out of mem)
template<
  class Kernel,
  class Tracer,
  class WeightCalculator,
  class Visitor,
  template <class> class LookupTable = Lookup_table_map
>
class Triangulate_hole_polyline_DT
{
  struct Auto_count {
    typedef std::pair<typename Kernel::Point_3, int> result_type;

    Auto_count(int count = 0) : count(count)  { }
    result_type operator()(const typename Kernel::Point_3& p) const
    { return std::make_pair(p, count++); }
    mutable int count;
  };

public:
  typedef typename WeightCalculator::Weight                   Weight;
  typedef typename Kernel::Point_3                            Point_3;
  typedef std::vector<Point_3>                                Polyline_3;

  typedef Triangulation_vertex_base_with_info_3<int, Kernel>  VB_with_id;
  typedef Delaunay_triangulation_cell_base_3<Kernel>          CB;
  typedef Triangulation_data_structure_3<VB_with_id, CB>      TDS;
  typedef Delaunay_triangulation_3<Kernel, TDS>               Triangulation;

  typedef typename Triangulation::Finite_edges_iterator       Finite_edges_iterator;
  typedef typename Triangulation::Facet_circulator            Facet_circulator;
  typedef typename Triangulation::Cell_handle                 Cell_handle;
  typedef typename Triangulation::Vertex_handle               Vertex_handle;
  typedef typename Triangulation::Edge                        Edge;
  typedef typename Triangulation::Facet                       Facet;

  typedef Incident_facet_circulator<2, Triangulate_hole_polyline_DT> IFC_2;
  typedef Incident_facet_circulator<3, Triangulate_hole_polyline_DT> IFC_3;

  Weight operator()(const Polyline_3& P,
                    const Polyline_3& Q,
                    Tracer& tracer,
                    const WeightCalculator& WC,
                    Visitor& visitor) const
  {
    CGAL_assertion(P.front() == P.back());
    CGAL_assertion(Q.empty() || (Q.front() == Q.back()));
    CGAL_assertion(Q.empty() || (P.size() == Q.size()));

    int n = static_cast<int>(P.size())-1; // because the first and last point are equal
    Triangulation tr;
    std::vector<bool> edge_exist;
    std::pair<int, int> range(0, n-1);
    std::tuple<boost::optional<Edge>, bool, bool> res = construct_3D_triangulation(P, range, tr, edge_exist);
    if(!std::get<2>(res)) {
#ifdef CGAL_HOLE_FILLING_VERBOSE
      #ifndef CGAL_TEST_SUITE
      CGAL_warning_msg(false, "Returning no output. Dimension of 3D Triangulation is below 2!");
      #else
      std::cerr << "W: Returning no output. Dimension of 3D Triangulation is below 2!\n";
      #endif
#endif
      return Weight::NOT_VALID();
    }

    // all border edges inside 3D Triangulation
    if(std::get<1>(res)) {
      LookupTable<Weight> W(n, Weight::DEFAULT()); // do not forget that these default values are not changed for [i, i+1]
      LookupTable<int>    lambda(n,-1);

      typename Incident_facet_circulator_base<Triangulate_hole_polyline_DT>::Edge_wrapper
        e_start(*std::get<0>(res));
      if(tr.dimension() == 3) {
        visitor.start_quadratic_phase(tr.number_of_finite_facets());
        triangulate_DT<IFC_3>(P, Q, W, lambda, e_start, tr, WC, visitor, false);
      }
      else {
        CGAL_assertion(tr.dimension() == 2);
        triangulate_DT<IFC_2>(P, Q, W, lambda, e_start, tr, WC, visitor, false);
      }

      if(W.get(0, n-1) == Weight::NOT_VALID()) {
#ifdef CGAL_HOLE_FILLING_VERBOSE
        #ifndef CGAL_TEST_SUITE
        CGAL_warning_msg(false, "Returning no output. No possible triangulation is found!");
        #else
        std::cerr << "W: Returning no output. No possible triangulation is found!\n";
        #endif
#endif
        visitor.end_quadratic_phase(false);
        return Weight::NOT_VALID();
      }

      tracer(lambda, 0, n-1);
      visitor.end_quadratic_phase(true);
      return W.get(0,n-1);
    }

    // How to handle missing border edges
    #if 1
    visitor.start_quadratic_phase(tr.number_of_finite_facets());
    Weight w = fill_by_extra_triangles(tr, edge_exist, P, Q, tracer, WC, visitor);
    #else
    // This approach produce better patches when used with Weight_incomplete
    // (which should be arranged in internal::triangulate_hole_Polyhedron, triangulate_polyline)
    Weight w = fill_by_incomplete_patches(tr, std::get<0>(res), edge_exist, P, Q, tracer, WC, visitor);
    #endif

    visitor.end_quadratic_phase(w != Weight::NOT_VALID());

    return w;
  }

private:
  /************************************************************************
  * Main algorithm which construct a minimum patch top-down searching through the space of tr
  *
  * + Edge_DT should have:
  *   - vertex_first() vertex_second() functions where vertex_first() always returns vertex with smaller id
  * + IncidentFacetCirculator should have:
  *   - constructor with Edge_DT and Triangulation_DT
  *   - pre-increment, conversion to bool
  *   - get_third() returning the third vertex of facet pointed by circulator
  *   - edge_first() and edge_second() neighbor edges to get_third() vertex
  ************************************************************************/
  template<class IncidentFacetCirculator, class Edge_DT, class Triangulation_DT>
  void triangulate_DT(const Polyline_3& P,
                      const Polyline_3& Q,
                      LookupTable<Weight>& W,
                      LookupTable<int>& lambda,
                      Edge_DT e,
                      const Triangulation_DT& tr,
                      const WeightCalculator& WC,
                      Visitor& visitor,
                      const bool produce_incomplete) const
  {
    /**********************************************************************
     *  + Default W value is Weight::DEFAULT(), default lambda value is -1.
     *  + DEFAULT() is used to check whether the region (v0-v1) is processed.
     *  + If a range v0-v1 does not contains any possible triangulation, then W[v0,v1] = NOT_VALID() and lambda[v0,v1] = -1
     *  + Note that w + DEFAULT() == w must hold
     */
    int v0 = e.vertex_first();
    int v1 = e.vertex_second();
    CGAL_assertion(v0 < v1);  // vertex_first() should always return vertex with smaller index
    CGAL_assertion(v0 != -1); // edge can not be incident to infinite vertex

    if( v0 + 1 == v1 || // border edge - should not check v0 = 0, v1 = n-1, because it is the initial edge where the algorithm starts
        W.get(v0, v1) != Weight::DEFAULT() ) // the range is previously processed
    { return; }

    visitor.quadratic_step();

    int m_min = -1;
    Weight w_min = Weight::NOT_VALID();

    IncidentFacetCirculator fb(e, tr);
    do {
      int v2 = fb.get_third();
      if(v2 < v0 || v2 > v1) { continue; } // this will also skip infinite vertex

      if(WC(P,Q, v0,v2,v1, lambda) == Weight::NOT_VALID())
      { continue; } // computed weight in here is not correct weight
                    // since max dih angle requires neighbor ranges to be already computed. It is just for checking validity.

      Weight w = Weight::DEFAULT();

      Edge_DT e0 = fb.edge_first(); // edge v0-v2
      CGAL_assertion(e0.vertex_first() == v0 && e0.vertex_second() == v2);
      triangulate_DT<IncidentFacetCirculator>(P, Q, W, lambda, e0, tr, WC, visitor, produce_incomplete); // region v0-v2
      const Weight& we0 = W.get(v0, v2);

      if(!produce_incomplete && we0 == Weight::NOT_VALID())
      { continue; } // not producing incomplete patches and failed to fill sub-range v0-v2, so no reason to proceed
      if(we0 != Weight::NOT_VALID()) // to not consider we0 if it is NOT_VALID (it is valid when produce_incomplete = true)
      { w = w + we0; }

      Edge_DT e1 = fb.edge_second(); // edge v2-v1
      CGAL_assertion(e1.vertex_first() == v2 && e1.vertex_second() == v1);
      triangulate_DT<IncidentFacetCirculator>(P, Q, W, lambda, e1, tr, WC, visitor, produce_incomplete); // region v2-v1
      const Weight& we1 = W.get(v2, v1);

      if(!produce_incomplete && we1 == Weight::NOT_VALID())
      { continue; }
      if(we1 != Weight::NOT_VALID()) // to not consider we1 if it is NOT_VALID (it is valid when produce_incomplete = true)
      { w = w + we1; }

      w = w + WC(P,Q, v0,v2,v1, lambda);
      if(m_min == -1 || w < w_min){
        w_min = w;
        m_min = v2;
      }
    } while(++fb);

    // can be m_min = -1 and w_min = NOT_VALID which means no possible triangulation between v0-v1
    W.put(v0,v1, w_min);
    lambda.put(v0,v1, m_min);
  }

  // returns [h.first-h.second edge, true if all edges inside 3D triangulation, true if tr.dimension() >= 2]
  std::tuple<boost::optional<Edge>, bool, bool>
  construct_3D_triangulation(const Polyline_3& P,
                             std::pair<int,int> h,
                             Triangulation& tr,
                             std::vector<bool>& edge_exist) const
  {
    // construct 3D tr with P[h.first], P[h.second] also assign ids from h.first to h.second
    boost::optional<Edge> e;
    int n_border = h.second - h.first + 1;
    tr.insert(boost::make_transform_iterator(std::next(P.begin(), h.first), Auto_count(h.first)),
              boost::make_transform_iterator(std::next(P.begin(), h.second +1), Auto_count(h.first)));
    tr.infinite_vertex()->info() = -1;

    if(tr.dimension() < 2) { return std::make_tuple(e, false, false); }

    // check whether all edges are included in DT, and get v0-vn-1 edge
    edge_exist.assign(n_border, false);
    int nb_exists = 0;
    Finite_edges_iterator v_first_v_second_edge; // range.first - range.second edge
    for(Finite_edges_iterator eb = tr.finite_edges_begin();
        eb != tr.finite_edges_end();
        ++eb)
    {
      int v0_id = eb->first->vertex(eb->second)->info();
      int v1_id = eb->first->vertex(eb->third )->info();
      if(v0_id > v1_id) { std::swap(v0_id, v1_id); }

      // to start from v0 vn-1 edge
      if(v0_id == h.first && v1_id == h.second) { v_first_v_second_edge = eb; }

      // check whether the edge is border edge
      int border_id = -1;
      if(v0_id + 1 == v1_id)                         { border_id = v0_id; }
      else if(v0_id == h.first && v1_id == h.second) { border_id = v1_id; }

      if(border_id != -1 && !edge_exist[border_id - h.first]) {
        ++nb_exists;
        edge_exist[border_id - h.first] = true;
      }
    }
    CGAL_assertion(n_border >= nb_exists);

    bool is_3D_T_complete = (nb_exists == n_border);
    if(edge_exist[n_border-1]) { e = *v_first_v_second_edge; }
    return std::make_tuple(e, is_3D_T_complete, true);
  }

  /************************************************************************
  * Try to construct hole part by part.
  *
  * What need to be improved:
  *  + if 3D triangulation does not contain the start-edge (edge between v0 vn-1) we directly switch to all space.
  *  + when switched to all-space, we use map based lookup tables.
  ************************************************************************/
  Weight fill_by_incomplete_patches(Triangulation& tr,
                                    boost::optional<Edge> start_edge,
                                    std::vector<bool>& edge_exist,
                                    const Polyline_3& P,
                                    const Polyline_3& Q,
                                    Tracer& tracer,
                                    const WeightCalculator& WC,
                                    Visitor& visitor) const
  {
    typedef std::pair<int, int> Range;
    typedef std::back_insert_iterator<std::vector<Range> > Output_hole_iterator;
    typedef Tracer_polyline_incomplete<std::tuple<int, int, int>, Emptyset_iterator, Output_hole_iterator> Remaining_holes_tracer;

    std::vector<Range> remaining_holes;

    int n_all = P.size()-1;// because the first point and last point are equal
    remaining_holes.push_back(Range(0, n_all-1)); // corresponds to start_edge

    LookupTable<Weight> W(n_all, Weight::DEFAULT()); // do not forget that these default values are not changed for [i, i+1]
    LookupTable<int>    lambda(n_all,-1);

    while(true) {
      Range h = remaining_holes.back();
      remaining_holes.pop_back();

      if(start_edge) {
        typename IFC_3::Edge_wrapper e(*start_edge);
        CGAL_assertion(h.first == e.vertex_first() &&
                       h.second == e.vertex_second());
      }

      if(!start_edge) {
        // switch to brute force
        Triangulate_hole_polyline<Kernel, Tracer, WeightCalculator, Visitor, LookupTable> all_space;
        all_space.triangulate_all(P, Q, WC, visitor, std::make_pair(h.first,  h.second), W, lambda);
        if(W.get(h.first, h.second) == Weight::NOT_VALID()) {
#ifdef CGAL_HOLE_FILLING_VERBOSE
          CGAL_warning_msg(false, "Returning no output. Filling hole with incomplete patches is not successful!");
#endif
          return Weight::NOT_VALID();
        }
      }
      else {
        // run the algorithm
        typename IFC_3::Edge_wrapper e(*start_edge);
        if(tr.dimension() == 3)
        {
          triangulate_DT<IFC_3>(P, Q, W, lambda, e, tr, WC, true);
        }
        else {
          CGAL_assertion(tr.dimension() == 2);
          triangulate_DT<IFC_2>(P, Q, W, lambda, e, tr, WC, true);
        }
        // check whether there is any improvement (at least we should construct one triangle)
        if(W.get(h.first, h.second) == Weight::NOT_VALID()) {
          // switch to brute force
          Triangulate_hole_polyline<Kernel, Tracer, WeightCalculator, Visitor, LookupTable> all_space;
          all_space.triangulate_all(P, Q, WC, visitor, std::make_pair(h.first, h.second), W, lambda);
          if(W.get(h.first, h.second) == Weight::NOT_VALID()) {
#ifdef CGAL_HOLE_FILLING_VERBOSE
            CGAL_warning_msg(false, "Returning no output. Filling hole with incomplete patches is not successful!");
#endif
            return Weight::NOT_VALID();
          }
        }
        // gather remaining holes
        Remaining_holes_tracer hole_tracer((Emptyset_iterator()), (Output_hole_iterator(remaining_holes)));
        hole_tracer(lambda, e.vertex_first(), e.vertex_second());
      }

      if(remaining_holes.empty()) { break; }
      // construct tr for next coming hole
      h = remaining_holes.back();
      tr.clear();
      std::tuple<boost::optional<Edge>, bool, bool> res = construct_3D_triangulation(P, h, tr, edge_exist);
      if(!std::get<0>(res)) {
#ifdef CGAL_HOLE_FILLING_VERBOSE
        CGAL_warning_msg(false, "Returning no output. Filling hole with incomplete patches is not successful!");
#endif
        return Weight::NOT_VALID();
      }
      start_edge = *std::get<0>(res);
      // clear related regions in W, lambda for next coming hole
      W.set_range_to_default(h.first, h.second);
      lambda.set_range_to_default(h.first, h.second);
    }
    tracer(lambda, 0, n_all-1);

    // W.get(0, n_all -1) is not correct weight (since we do not update weights while we are filling remaining holes),
    // we need to recalculate it
    std::stack<std::pair<int, int> > ranges;
    ranges.push(std::make_pair(0, n_all-1));
    Weight total_weight = Weight::DEFAULT();
    while(!ranges.empty()) {
      std::pair<int, int> r = ranges.top();
      ranges.pop();
      if(r.first + 1 == r.second) { continue; }
      int la = lambda.get(r.first, r.second);
      total_weight = total_weight + WC(P,Q, r.first,la,r.second, lambda);
      ranges.push(std::make_pair(r.first, la));
      ranges.push(std::make_pair(la, r.second));
    }
    return total_weight;
  }

  /************************************************************************
   * This approach extends the search space by adding extra triangles.
   *
   * + Initial search space is 3D Triangulation.
   * + For each border edge which is not inside 3DT, we add all possible triangles containing that edge to the search space.
   *   Example: say border edge [4-5] is not found in 3DT, then triangles = { [0,4,5] [1,4,5] [2,4,5] ... [ n-1,4,5] }
   *   are added to the search space.
   *
   * I guess this approach does not make the search space complete, since there are some cases that it still returns no patch.
   ************************************************************************/
  Weight fill_by_extra_triangles(const Triangulation& tr,
                                 const std::vector<bool>& edge_exist,
                                 const Polyline_3& P,
                                 const Polyline_3& Q,
                                 Tracer& tracer,
                                 const WeightCalculator& WC,
                                 Visitor& visitor) const
  {
    int n = static_cast<int>(edge_exist.size());
    LookupTable<Weight> W(n, Weight::DEFAULT()); // do not forget that these default values are not changed for [i, i+1]
    LookupTable<int>    lambda(n,-1);

    Edge_graph edge_graph;

    if(tr.dimension() == 3)
    { edge_graph.init<IFC_3>(tr, edge_exist); }
    else
    { edge_graph.init<IFC_2>(tr, edge_exist); }

    Edge_graph::Edge_wrapper e_start(std::make_pair(0, n-1));
    triangulate_DT<Edge_graph::Incident_facet_circulator>
      (P, Q, W, lambda, e_start, edge_graph, WC, visitor, false);


    if(W.get(0, n-1) == Weight::NOT_VALID()) {
#ifdef CGAL_HOLE_FILLING_VERBOSE
      #ifndef CGAL_TEST_SUITE
      CGAL_warning_msg(false, "Returning no output using Delaunay triangulation.\n Falling back to the general Triangulation framework.");
      #else
      std::cerr << "W: Returning no output using Delaunay triangulation.\n"
                << "Falling back to the general Triangulation framework.\n";
      #endif
#endif
      return Weight::NOT_VALID();
    }

    tracer(lambda, 0, n-1);
    return W.get(0, n-1);
  }
}; // End of Triangulate_hole_polyline_DT
#endif

/************************************************************************
 * Triangulate hole by using all search space
 ************************************************************************/
template<
  class Kernel,
  class Tracer,
  class WeightCalculator,
  class Visitor,
  template <class> class LookupTable
>
class Triangulate_hole_polyline {
public:
  typedef typename WeightCalculator::Weight  Weight;
  typedef typename Kernel::Point_3           Point_3;
  typedef std::vector<Point_3>               Polyline_3;

  Weight operator()(const Polyline_3& P,
                    const Polyline_3& Q,
                    Tracer& tracer,
                    const WeightCalculator& WC,
                    Visitor& visitor) const
  {
    CGAL_assertion(P.front() == P.back());
    CGAL_assertion(Q.empty() || (Q.front() == Q.back()));
    CGAL_assertion(Q.empty() || (P.size() == Q.size()));

    int n = static_cast<int>(P.size()) - 1;                       // because the first and last point are equal
    LookupTable<Weight> W(n,Weight::DEFAULT()); // do not forget that these default values are not changed for [i, i+1]
    LookupTable<int>    lambda(n,-1);

    triangulate_all(P, Q, WC, visitor, std::make_pair(0,n-1), W, lambda);

    if(W.get(0,n-1) == Weight::NOT_VALID() || n <= 2) {
#ifdef CGAL_HOLE_FILLING_VERBOSE
      #ifndef CGAL_TEST_SUITE
      CGAL_warning_msg(false, "Returning no output. No possible triangulation is found!");
      #else
      std::cerr << "W: Returning no output. No possible triangulation is found!\n";
      #endif
#endif
      return Weight::NOT_VALID();
    }

    tracer(lambda, 0, n-1);
    return W.get(0,n-1);
  }

  void triangulate_all(const Polyline_3& P,
                       const Polyline_3& Q,
                       const WeightCalculator& WC,
                       Visitor& visitor,
                       std::pair<int, int> range,
                       LookupTable<Weight>& W,
                       LookupTable<int>& lambda) const
  {
    int f = range.first, s = range.second;

    const int N = s * ( -3* f * (s - 1) + s * s - 1) /6;

    visitor.start_cubic_phase(N);
    for(int j = 2; j<= range.second; ++j) {              // determines range (2 - 3 - 4 )
      for(int i=range.first; i<= range.second-j; ++i) {  // iterates over ranges and find min triangulation in those ranges
        int k = i+j;                                     // like [0-2, 1-3, 2-4, ...], [0-3, 1-4, 2-5, ...]

        int m_min = -1;
        Weight w_min = Weight::NOT_VALID();
        // i is the range start (e.g. 1) k is the range end (e.g. 5) -> [1-5]. Now subdivide the region [1-5] with m -> 2,3,4
        for(int m = i+1; m<k; ++m) {
          visitor.cubic_step();
          // now the regions i-m and m-k might be valid(constructed) patches,
          if( W.get(i,m) == Weight::NOT_VALID() || W.get(m,k) == Weight::NOT_VALID() )
          { continue; }

          const Weight& w_imk = WC(P,Q,i,m,k, lambda);
          if(w_imk == Weight::NOT_VALID())
          { continue; }

          const Weight& w = W.get(i,m) + W.get(m,k) + w_imk;
          if(m_min == -1 || w < w_min) {
            w_min = w;
            m_min = m;
          }
        }

        // can be m_min = -1 and w_min = NOT_VALID which means no possible triangulation between i-k
        W.put(i,k,w_min);
        lambda.put(i,k, m_min);
      }
    }
    visitor.end_cubic_phase();
  }
};

#ifndef CGAL_HOLE_FILLING_DO_NOT_USE_CDT2

/************************************************************************
 * Triangulate hole by using a cdt_2
 ************************************************************************/
// /!\ points.first == points.last

template <typename Traits>
bool is_planar_2(
  const std::vector<typename Traits::Point_3>& points,
  const typename Traits::Vector_3& avg_normal,
  const typename Traits::FT max_squared_distance,
  const Traits& traits) {

  typedef typename Traits::FT FT;
  typedef typename Traits::Point_3 Point_3;
  typedef typename Traits::Plane_3 Plane_3;
  typedef typename Traits::Construct_projected_point_3 Projection_3;
  typedef typename Traits::Compute_squared_distance_3 Squared_distance_3;

  const Projection_3 projection_3 =
    traits.construct_projected_point_3_object();
  const Squared_distance_3 squared_distance_3 =
    traits.compute_squared_distance_3_object();

  const double n = static_cast<double>(points.size() - 1); // the first equals to the last
  if (n < 3) {
    return false; // cant be a plane!
  }

  // Compute centroid.
  const Point_3 centroid =
    CGAL::centroid(points.begin(), points.end() - 1);

  // Compute covariance matrix.
  FT xx = FT(0), yy = FT(0), zz = FT(0);
  FT xy = FT(0), xz = FT(0), yz = FT(0);
  for (const Point_3& p : points)
  {
    const FT dx = p.x() - centroid.x();
    const FT dy = p.y() - centroid.y();
    const FT dz = p.z() - centroid.z();
    xx += dx * dx; yy += dy * dy; zz += dz * dz;
    xy += dx * dy; xz += dx * dz; yz += dy * dz;
  }

  // Check the planarity.
  const FT x = yy * zz - yz * yz;
  const FT y = xx * zz - xz * xz;
  const FT z = xx * yy - xy * xy;
  FT maxv = -FT(1);
  maxv = (CGAL::max)(maxv, x);
  maxv = (CGAL::max)(maxv, y);
  maxv = (CGAL::max)(maxv, z);
  if (maxv <= FT(0)) {
    return false; // a good plane does not exist for sure!
  }

  // Here, avg_squared_distance is a little bit more tolerant than avg_distance^2.
  const Plane_3 plane = Plane_3(centroid, avg_normal);
  FT avg_squared_distance = FT(0);
  for (const Point_3& p : points)
  {
    const Point_3  q = projection_3(plane, p);
    avg_squared_distance += CGAL::abs(squared_distance_3(p, q));
  }
  avg_squared_distance /= FT(n);
  // std::cout << "avg squared distance: " << avg_squared_distance << std::endl;

  CGAL_assertion(max_squared_distance >= FT(0));
  if (avg_squared_distance > max_squared_distance) {
    return false; // the user distance criteria are not satisfied!
  }

  // std::cout << "The hole seems to be near planar." << std::endl;
  return true;
}

template <
  typename PointRange, // need size()
  typename Tracer,
  typename Visitor,
  typename Validity_checker,
  typename Traits
>
bool
triangulate_hole_polyline_with_cdt(const PointRange& points,
                                   Tracer& tracer,
                                   Visitor& visitor,
                                   const Validity_checker& is_valid,
                                   const Traits& traits,
                                   const typename Traits::FT max_squared_distance)
{
  typedef typename Traits::FT FT;
  typedef typename Traits::Point_3 Point_3;
  typedef typename Traits::Vector_3 Vector_3;
  typedef typename Traits::Collinear_3 Collinear_3;

  visitor.start_planar_phase();

  // Compute an average normal of the hole.
  const Collinear_3 collinear_3 =
    traits.collinear_3_object();

  std::vector<Point_3> P(std::begin(points), std::end(points));
  CGAL_assertion(P.size() >= 3);
  if (P.front() != P.back()) {
    P.push_back(P.front());
  }

  FT x = FT(0), y = FT(0), z = FT(0);
  int num_normals = 0;
  const Point_3& ref_point = P[0];
  const std::size_t size = P.size() - 1;
  for (std::size_t i = 1; i < size - 1; ++i) {
    const std::size_t ip = i + 1;

    const Point_3& p1 = ref_point; // 3 points, which form a triangle
    const Point_3& p2 = P[i];
    const Point_3& p3 = P[ip];

    // Skip in case we have collinear points.
    if (collinear_3(p1, p2, p3)) {
      continue;
    }

    // Computing the normal of a triangle.
    const Vector_3 n = CGAL::normal(p1, p2, p3);
    // If it is a positive normal ->
    if (
      ( n.x() >  FT(0)                                    ) ||
      ( n.x() == FT(0) && n.y() >  FT(0)                  ) ||
      ( n.x() == FT(0) && n.y() == FT(0) && n.z() > FT(0) )) {
      x += n.x(); y += n.y(); z += n.z();
    } else { // otherwise invert ->
      x -= n.x(); y -= n.y(); z -= n.z();
    }
    ++num_normals;
  }

  if (num_normals < 1) {
    // std::cerr << "WARNING: num normals, cdt 2 falls back to the original solution!" << std::endl;
    visitor.end_planar_phase(false);
    return false;
  }

  // Setting the final normal.
  FT ft_nn(num_normals);
  x /= ft_nn;
  y /= ft_nn;
  z /= ft_nn;
  const Vector_3 avg_normal = Vector_3(x, y, z);
  // std::cout << "avg normal: " << avg_normal << std::endl;

  if (avg_normal==NULL_VECTOR){
    visitor.end_planar_phase(false);
    return false;
  }
  // Checking the hole planarity.
  if (!is_planar_2(P, avg_normal, max_squared_distance, traits)) {
    // std::cerr << "WARNING: planarity, cdt 2 falls back to the original solution!" << std::endl;
        visitor.end_planar_phase(false);
    return false;
  }

  // Checking the hole simplicity.
  typedef CGAL::Projection_traits_3<Traits> P_traits;
  const P_traits p_traits(avg_normal);
  if (!is_simple_2(P.begin(), P.end() - 1, p_traits)) {
    // std::cerr << "WARNING: simplicity, cdt 2 falls back to the original solution!" << std::endl;
    visitor.end_planar_phase(false);
    return false;
  }

  Lookup_table_map<int> lambda(static_cast<int>(size), -1);

  // Create and fill the cdt_2.
  typedef CGAL::Triangulation_vertex_base_with_info_2<std::size_t, P_traits> Vb;
  typedef CGAL::Triangulation_face_base_with_info_2<bool, P_traits>          Fbi;
  typedef CGAL::Constrained_triangulation_face_base_2<P_traits, Fbi>         Fb;
  typedef CGAL::Triangulation_data_structure_2<Vb, Fb>                       TDS;
  // If the polygon is simple, there should be no intersection.
  typedef CGAL::No_constraint_intersection_tag                               Itag;
  typedef CGAL::Constrained_Delaunay_triangulation_2<P_traits, TDS, Itag>    CDT;
  P_traits cdt_traits(avg_normal);
  CDT cdt(cdt_traits);

  std::vector< std::pair<Point_3, std::size_t> > points_and_ids;
  points_and_ids.reserve(size);
  for (std::size_t i = 0; i < size; ++i) {
    points_and_ids.push_back(std::make_pair(P[i], i));
  }

  std::vector<typename CDT::Vertex_handle> vertices(size);
  cdt.insert(points_and_ids.begin(), points_and_ids.end());
  for (typename CDT::Vertex_handle v : cdt.finite_vertex_handles()) {
    vertices[v->info()] = v;
  }

  for (std::size_t i = 0; i < size; ++i) {
    const std::size_t ip = (i + 1) % size;
    if (vertices[i] != vertices[ip]) {
      cdt.insert_constraint(vertices[i], vertices[ip]);
    }
  }

  // Mark external faces.
  for (typename CDT::All_faces_iterator fit = cdt.all_faces_begin(),
    end = cdt.all_faces_end(); fit != end; ++fit) {
    fit->info() = false;
  }

  std::queue<typename CDT::Face_handle> face_queue;
  face_queue.push(cdt.infinite_vertex()->face());
  while (!face_queue.empty()) {

    typename CDT::Face_handle fh = face_queue.front();
    face_queue.pop();
    if (fh->info()) {
      continue;
    }

    fh->info() = true;
    for (int i = 0; i < 3; ++i) {
      if (!cdt.is_constrained(typename CDT::Edge(fh, i))) {
        face_queue.push(fh->neighbor(i));
      }
    }
  }

  if (cdt.dimension() != 2 || cdt.number_of_vertices() != size) {
    // std::cerr << "WARNING: dim + num vertices, cdt 2 falls back to the original solution!" << std::endl;
    visitor.end_planar_phase(false);
    return false;
  }

  // Fill the lambda.
  for (typename CDT::Finite_faces_iterator fit = cdt.finite_faces_begin(),
    end = cdt.finite_faces_end(); fit != end; ++fit) {
    if (!fit->info()) { // if it is not external

      std::vector<int> is(3);
      for (int i = 0; i < 3; ++i) {
        is[i] = static_cast<int>(fit->vertex(i)->info());
      }

      std::sort(is.begin(), is.end());
      lambda.put(is[0], is[2], is[1]);
      if (!is_valid(P, is[0], is[1], is[2])) {
        // std::cerr << "WARNING: validity, cdt 2 falls back to the original solution!" << std::endl;
        visitor.end_planar_phase(false);
        return false;
      }
    }
  }

  // Call the tracer. It correctly orients the patch faces.
  // std::cout << "CDT is being used!" << std::endl;
  tracer(lambda, 0, static_cast<int>(size) - 1);
  visitor.end_planar_phase(true);
  return true;
}

#endif

/*******************************************************************************
 * Internal entry point for both polyline and Polyhedron_3 triangulation functions
 ***********************************************************************************/
template <
  typename PointRange1,
  typename PointRange2,
  typename Tracer,
  typename WeightCalculator,
  typename Visitor,
  typename Kernel
>
typename WeightCalculator::Weight
triangulate_hole_polyline(const PointRange1& points,
                          const PointRange2& third_points,
                          Tracer& tracer,
                          const WeightCalculator& WC,
                          Visitor& visitor,
                          bool use_delaunay_triangulation,
                          bool skip_cubic_algorithm,
                          const Kernel&)
{
  CGAL_assertion(!points.empty());
  if (!use_delaunay_triangulation && skip_cubic_algorithm)
    return WeightCalculator::Weight::NOT_VALID();

  typedef Kernel        K;
  typedef typename K::Point_3    Point_3;
#ifndef CGAL_HOLE_FILLING_DO_NOT_USE_DT3
  typedef CGAL::internal::Triangulate_hole_polyline_DT<K, Tracer, WeightCalculator, Visitor> Fill_DT;
#else
  CGAL_USE(use_delaunay_triangulation);
#endif
  typedef CGAL::internal::Triangulate_hole_polyline<K, Tracer, WeightCalculator, Visitor>    Fill;

  std::vector<Point_3> P(boost::begin(points), boost::end(points));
  std::vector<Point_3> Q(boost::begin(third_points), boost::end(third_points));

  if(P.front() != P.back()){
    P.push_back(P.front());
    if( !Q.empty() && P.size() > Q.size()) {
      Q.push_back(Q.front());
    }
  }

  typename WeightCalculator::Weight w =
#ifndef CGAL_HOLE_FILLING_DO_NOT_USE_DT3
    use_delaunay_triangulation ? Fill_DT().operator()(P, Q, tracer, WC, visitor) :
#endif
    Fill().operator()(P, Q, tracer, WC, visitor);

#ifndef CGAL_HOLE_FILLING_DO_NOT_USE_DT3
  if(use_delaunay_triangulation &&
     w == WeightCalculator::Weight::NOT_VALID()
     &&!skip_cubic_algorithm)
  {
    w = Fill().operator()(P, Q, tracer, WC, visitor);
  }
#endif

#ifdef CGAL_PMP_HOLE_FILLING_DEBUG
  std::cerr << w << std::endl;
#endif
  return w;
}

} // namespace internal

} // namespace CGAL

#endif // CGAL_HOLE_FILLING_TRIANGULATE_HOLE_POLYLINE_H
