// Copyright (c) 2014  INRIA Sophia-Antipolis (France), INRIA Lorraine LORIA.
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
// Author(s)     : Fernando de Goes, Pierre Alliez, Ivo Vigan, Clément Jamin

#ifndef CGAL_RECONSTRUCTION_TRIANGULATION_2_H
#define CGAL_RECONSTRUCTION_TRIANGULATION_2_H

#include <CGAL/license/Optimal_transportation_reconstruction_2.h>


// local
#include <CGAL/OTR_2/Sample.h>
#include <CGAL/OTR_2/Reconstruction_edge_2.h>
#include <CGAL/OTR_2/Cost.h>
#include <CGAL/OTR_2/Reconstruction_vertex_base_2.h>
#include <CGAL/OTR_2/Reconstruction_face_base_2.h>

// CGAL
#include <CGAL/basic.h>
#include <CGAL/Random.h>
#include <CGAL/intersections.h>
#include <CGAL/Delaunay_triangulation_2.h>

// boost
#include <boost/multi_index/mem_fun.hpp>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/identity.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/container/deque.hpp>
#include <boost/optional.hpp>

// STL
#include <map>
#include <set>
#include <vector>
#include <queue>
#include <iostream>

#define CGAL_EPS   1e-15

namespace CGAL {
namespace OTR_2 {

/// \internal
///  The Reconstruction_triangulation_2 class
///  provides the reconstruction simplex as well as the transportation plan.
/// - Each vertex stores a normal vector.
/// - A vertex a Sample which got assigned to it by the transportation plan,
///   well as the corresponding relocated Point (of type Traits_::Point_2).
/// - In order to solve a linear system over the triangulation, a vertex may be constrained
///   or not (i.e. may contribute to the right or left member of the linear system),
///   and has a unique index.
/// The vertex class must derive from Reconstruction_vertex_base_3.
///
///  @param Traits_   The traits class
///  @param Tds       Model of TriangulationDataStructure_2.
///  The vertex class must derive from Reconstruction_vertex_base_2.
///  The face   class must derive from Reconstruction_face_base_2.
template<
  class Traits_,
  class Tds_ = Triangulation_data_structure_2<
  Reconstruction_vertex_base_2<Traits_>,
  Reconstruction_face_base_2<Traits_> > >
class Reconstruction_triangulation_2
  : public Delaunay_triangulation_2<Traits_, Tds_> 
{
public:

  typedef Delaunay_triangulation_2<Traits_, Tds_> Base;

  typedef typename Traits_::FT FT;
  typedef typename Traits_::Point_2 Point;
  typedef typename Traits_::Vector_2 Vector;
  typedef typename Traits_::Line_2 Line;
  typedef typename Traits_::Segment_2 Segment;
  typedef typename Traits_::Triangle_2 Triangle;

  typedef typename Base::Vertex Vertex;
  typedef typename Base::Vertex_handle Vertex_handle;
  typedef typename Base::Vertex_iterator Vertex_iterator;
  typedef typename Base::Vertex_circulator Vertex_circulator;
  typedef typename Base::Finite_vertices_iterator Finite_vertices_iterator;

  typedef typename Base::Edge Edge;
  typedef typename Base::Edge_circulator Edge_circulator;
  typedef typename Base::Finite_edges_iterator Finite_edges_iterator;

  typedef typename Base::Face_handle Face_handle;
  typedef typename Base::Face_circulator Face_circulator;
  typedef typename Base::Finite_faces_iterator Finite_faces_iterator;

  typedef std::map<Vertex_handle, Vertex_handle,
      less_Vertex_handle<Vertex_handle> > Vertex_handle_map;
  typedef std::map<Face_handle, Face_handle, 
                   less_Face_handle<Face_handle> > Face_handle_map;

  typedef std::set<Vertex_handle, 
                   less_Vertex_handle<Vertex_handle> > Vertex_handle_set;
  typedef std::set<Edge, less_Edge<Edge> > Edge_set;

  typedef std::vector<Edge> Edge_vector;

  typedef OTR_2::Cost<FT> Cost_;
  typedef OTR_2::Sample<Traits_> Sample_;
  typedef std::vector<Sample_*> Sample_vector;
  typedef typename Sample_vector::const_iterator Sample_vector_const_iterator;

  typedef OTR_2::Sample_with_priority<Sample_> PSample;
  typedef std::priority_queue<PSample, boost::container::deque<PSample>,
      OTR_2::greater_priority<PSample> > SQueue;

  typedef Reconstruction_edge_2<FT, Edge, 
                                Vertex_handle, Face_handle> Rec_edge_2;

  using Base::geom_traits;

  typedef boost::multi_index_container<
      Rec_edge_2,
      boost::multi_index::indexed_by<
      // sort by Rec_edge_2::operator<
      boost::multi_index::ordered_unique< boost::multi_index::identity<
      Rec_edge_2 > > ,
      // sort by Rec_edge_2::priority()
      boost::multi_index::ordered_non_unique<
      boost::multi_index::const_mem_fun<
      Rec_edge_2,const FT,&Rec_edge_2::priority> >
  >
  > MultiIndex;

  FT m_factor; // ghost vs solid
  mutable Random rng;

public:
  Reconstruction_triangulation_2(Traits_ traits = Traits_())
  : Base(traits), m_factor(1.)
  {
  }

  ~Reconstruction_triangulation_2() {
  }

  FT& ghost_factor() {
    return m_factor;
  }

  FT ghost_factor() const {
    return m_factor;
  }

  Edge random_finite_edge() const {
    std::size_t nbf = Base::number_of_faces();
    int offset = rng.get_int(0, static_cast<int>(nbf - 1));
    Finite_faces_iterator fi = Base::finite_faces_begin();
    for (int i = 0; i < offset; i++)
      fi++;
    Face_handle face = fi;
    int index = rng.get_int(0, 40) % 3;
    return Edge(face, index);
  }

  // ACCESS //

  Vertex_handle source_vertex(const Edge& edge) const {
    return edge.first->vertex(Base::ccw(edge.second));
  }

  Vertex_handle target_vertex(const Edge& edge) const {
    return edge.first->vertex(Base::cw(edge.second));
  }

  Vertex_handle opposite_vertex(const Edge& edge) const {
    return edge.first->vertex(edge.second);
  }

  bool is_pinned(const Edge& edge) const {
    Vertex_handle s = source_vertex(edge);
    if (s->pinned())
      return true;
    return false;
  }

  Edge twin_edge(const Edge& edge) const {
    Face_handle f = edge.first;
    Vertex_handle v = source_vertex(edge);
    Face_handle nf = f->neighbor(edge.second);
    return Edge(nf, Base::ccw(nf->index(v)));
  }

  Edge next_edge(const Edge& edge) const {
    Face_handle f = edge.first;
    int index = Base::ccw(edge.second);
    return Edge(f, index);
  }

  Edge prev_edge(const Edge& edge) const {
    Face_handle f = edge.first;
    int index = Base::cw(edge.second);
    return Edge(f, index);
  }

  FT get_length(const Edge& edge) const {
    Segment segment = get_segment(edge);
    return CGAL::approximate_sqrt(segment.squared_length());
  }

  Segment get_segment(const Edge& edge) const {
    const Point& ps = source_vertex(edge)->point();
    const Point& pt = target_vertex(edge)->point();
    return geom_traits().construct_segment_2_object()(ps, pt);
  }

  Triangle get_triangle(Face_handle face) const {
    Vertex_handle v0 = face->vertex(0);
    Vertex_handle v1 = face->vertex(1);
    Vertex_handle v2 = face->vertex(2);
    return geom_traits().construct_triangle_2_object()(
      v0->point(), v1->point(), v2->point());
  }

  // GET LINK //

  void get_vertices_from_edge_link(const Edge& edge,
      Vertex_handle_set& vertices) const {
    vertices.insert(opposite_vertex(edge));
    vertices.insert(opposite_vertex(twin_edge(edge)));
  }

  void get_vertices_from_vertex_link(Vertex_handle vertex,
      Vertex_handle_set& vertices) const {
    Vertex_circulator vcirc = Base::incident_vertices(vertex);
    Vertex_circulator vend = vcirc;
    CGAL_For_all(vcirc, vend)
    {
      Vertex_handle v = vcirc;
      vertices.insert(v);
    }
  }

  // boundary of star(vertex)
  // 'outward' chooses the orientation of the boundary
  void get_edges_from_star_minus_link(Vertex_handle vertex, Edge_vector& hull,
      bool outward = false) const {
    Face_circulator fcirc = Base::incident_faces(vertex);
    Face_circulator fend = fcirc;
    CGAL_For_all(fcirc, fend)
    {
      Face_handle face = fcirc;
      int index = face->index(vertex);
      Edge edge(face, index);
      if (outward)
        edge = twin_edge(edge);
      hull.push_back(edge);
    }
  }

  // ATTRIBUTES //

  bool is_ghost(const Edge& edge) const {
    return edge.first->ghost(edge.second);
  }

  int get_plan(const Edge& edge) const {
    return edge.first->plan(edge.second);
  }

  void set_plan(const Edge& edge, int simplex) const {
    edge.first->plan(edge.second) = simplex;
  }

  FT get_mass(const Edge& edge) const {
    return edge.first->mass(edge.second);
  }

  void set_mass(const Edge& edge, const FT mass) const {
    edge.first->mass(edge.second) = mass;
  }

  const Cost_& get_cost(const Edge& edge) const {
    return edge.first->cost(edge.second);
  }

  void set_vertex_cost(const Edge& edge, const Cost_& cost) const {
    edge.first->vertex_cost(edge.second) = cost;
  }

  void set_edge_cost(const Edge& edge, const Cost_& cost) const {
    edge.first->edge_cost(edge.second) = cost;
  }

  FT get_vertex_minus_edge_cost(const Edge& edge) const {
    const Cost_& vcost = edge.first->vertex_cost(edge.second);
    const Cost_& ecost = edge.first->edge_cost(edge.second);
    return vcost.finalize() - m_factor * ecost.finalize();
  }

  FT get_vertex_over_edge_cost(const Edge& edge) const {
    FT vvalue = edge.first->vertex_cost(edge.second).finalize();
    FT evalue = edge.first->edge_cost(edge.second).finalize();
    if (evalue == vvalue)
      return FT(1) / m_factor;
    return vvalue / (m_factor * evalue);
  }

  FT get_edge_relevance(const Edge& edge) const {
    FT M = get_mass(edge);
    if (M == FT(0))
      return FT(0);

    FT L = get_length(edge);
    FT cost = get_cost(edge).finalize();
    return M * L * L / cost;
  }

  FT get_density(const Edge& edge) const {
    FT length = get_length(edge);
    FT mass = get_mass(edge);
    return (mass / length);
  }

  unsigned int nb_samples(const Edge& edge) const {
    Edge twin = twin_edge(edge);
    return edge.first->samples(edge.second).size()
      + twin.first->samples(twin.second).size();
  }

  void collect_samples_from_edge(
    const Edge& edge, Sample_vector& samples) const 
  {
    const Sample_vector& edge_samples = edge.first->samples(edge.second);
    samples.insert(samples.end(), edge_samples.begin(), edge_samples.end());
  }

  void collect_samples_from_vertex(
    Vertex_handle vertex, Sample_vector& samples, bool cleanup) const 
  {
    Face_circulator fcirc = Base::incident_faces(vertex);
    Face_circulator fend = fcirc;
    CGAL_For_all(fcirc, fend)
    {
      Face_handle face = fcirc;
      int index = face->index(vertex);

      Edge edge(face, index);
      collect_samples_from_edge(edge, samples);

      Edge next = next_edge(edge);
      collect_samples_from_edge(next, samples);

      Edge prev = prev_edge(edge);
      collect_samples_from_edge(prev, samples);

      if (cleanup)
        face->clean_all_samples();
    }
    Sample_* sample = vertex->sample();
    if (sample)
      samples.push_back(sample);
    if (cleanup)
      vertex->set_sample(NULL);
  }

  void collect_all_samples(Sample_vector& samples) const {
    for (Finite_edges_iterator ei = Base::finite_edges_begin();
        ei != Base::finite_edges_end(); ++ei) {
      Edge edge = *ei;
      Edge twin = twin_edge(edge);
      collect_samples_from_edge(edge, samples);
      collect_samples_from_edge(twin, samples);
    }
    for (Finite_vertices_iterator vi = Base::finite_vertices_begin();
        vi != Base::finite_vertices_end(); ++vi) {
      Vertex_handle v = vi;
      Sample_* sample = v->sample();
      if (sample)
        samples.push_back(sample);
    }
  }

  void cleanup_assignments() {
    for (Finite_faces_iterator fi = Base::finite_faces_begin();
        fi != Base::finite_faces_end(); ++fi) {
      fi->clean_all_samples();
    }
    for (Finite_vertices_iterator vi = Base::finite_vertices_begin();
        vi != Base::finite_vertices_end(); ++vi) {
      vi->set_sample(NULL);
    }
  }

  // COST //

  Cost_ compute_total_cost() const {
    Cost_ sum;
    for (Finite_edges_iterator ei = Base::finite_edges_begin();
        ei != Base::finite_edges_end(); ++ei) {
      Edge edge = *ei;
      const Cost_& cost = get_cost(edge);
      sum.update_max(cost);
      sum.add(cost);
    }
    return sum;
  }

  Cost_ compute_cost_around_vertex(Vertex_handle vertex) const {
    Cost_ inner;
    Cost_ outer;
    Face_circulator fcirc = Base::incident_faces(vertex);
    Face_circulator fend = fcirc;
    CGAL_For_all(fcirc, fend)
    {
      Face_handle face = fcirc;
      int index = face->index(vertex);

      Edge edge(face, index);
      Cost_ cost = get_cost(edge);
      outer.update_max(cost);
      outer.add(cost);

      edge = next_edge(edge);
      cost = get_cost(edge);
      inner.update_max(cost);
      inner.add(cost);

      edge = next_edge(edge);
      cost = get_cost(edge);
      inner.update_max(cost);
      inner.add(cost);
    }
    inner.divide(2.0);

    Cost_ sum;
    sum.add(inner);
    sum.add(outer);
    sum.update_max(inner);
    sum.update_max(outer);
    return sum;
  }

  void reset_all_costs() {
    for (Finite_edges_iterator ei = Base::finite_edges_begin();
        ei != Base::finite_edges_end(); ++ei) {
      Edge edge = *ei;
      update_cost(edge);
    }
  }

  void update_cost(const Edge& edge) {
    compute_mass(edge);
    compute_edge_cost(edge);
    compute_vertex_cost(edge);
    select_plan(edge);
  }

  void compute_mass(const Edge& edge) const {
    FT mass = 0;

    typename Sample_vector::const_iterator it;
    const Sample_vector& samples0 = edge.first->samples(edge.second);
    for (it = samples0.begin(); it != samples0.end(); ++it) {
      Sample_* sample = *it;
      mass += sample->mass();
    }

    Edge twin = twin_edge(edge);
    const Sample_vector& samples1 = twin.first->samples(twin.second);
    for (it = samples1.begin(); it != samples1.end(); ++it) {
      Sample_* sample = *it;
      mass += sample->mass();
    }

    set_mass(edge, mass);
    set_mass(twin, mass);
  }

  void select_plan(const Edge& edge) const {
    // transport plan:
    // 0 - to vertex
    // 1 - to edge

    int plan = 0;
    FT diff = get_vertex_minus_edge_cost(edge);
    if (diff >= 0.0)
      plan = 1;

    Edge twin = twin_edge(edge);
    set_plan(edge, plan);
    set_plan(twin, plan);
  }

  void compute_edge_cost(const Edge& edge) const {
    SQueue squeue;
    FT M = get_mass(edge);
    FT L = get_length(edge);
    sort_samples_from_edge(edge, squeue);
    Cost_ cost = compute_cost_from_squeue(squeue, M, L);

    Edge twin = twin_edge(edge);
    set_edge_cost(edge, cost);
    set_edge_cost(twin, cost);
  }

  void sort_samples_from_edge(const Edge& edge, SQueue& squeue) const {
    typename Sample_vector::const_iterator it;
    const Sample_vector& samples0 = edge.first->samples(edge.second);
    for (it = samples0.begin(); it != samples0.end(); ++it) {
      Sample_* sample = *it;
      squeue.push(PSample(sample, sample->coordinate()));
    }

    Edge twin = twin_edge(edge);
    const Sample_vector& samples1 = twin.first->samples(twin.second);
    for (it = samples1.begin(); it != samples1.end(); ++it) {
      Sample_* sample = *it;
      squeue.push(PSample(sample, 1.0 - sample->coordinate()));
    }
  }

  Cost_ compute_cost_from_squeue(SQueue& squeue, const FT M, const FT L) const
  {
    if (squeue.empty())
      return Cost_();
    if (M == FT(0))
      return Cost_();

    Cost_ sum;
    FT start = 0;
    FT coef = L / M;
    while (!squeue.empty()) {
      PSample psample = squeue.top();
      squeue.pop();

      FT mass = psample.sample()->mass();
      FT coord = psample.priority() * L;
      FT bin = mass * coef;
      FT center = start + FT(0.5) * bin;
      FT pos = coord - center;

      FT norm2 = psample.sample()->distance2();
      FT tang2 = bin * bin / 12 + pos * pos;

      sum.add(Cost_(norm2, tang2), mass);
      sum.compute_max(norm2, tang2);

      start += bin;
    }
    return sum;
  }

  void compute_vertex_cost(const Edge& edge) const {
    Edge twin = twin_edge(edge);
    const Point& ps = source_vertex(edge)->point();
    const Point& pt = target_vertex(edge)->point();

    Sample_vector samples;
    collect_samples_from_edge(edge, samples);
    collect_samples_from_edge(twin, samples);

    Cost_ sum;
    for (Sample_vector_const_iterator it = samples.begin();
        it != samples.end(); ++it) {
      Sample_* sample = *it;
      FT mass = sample->mass();
      const Point& query = sample->point();

      FT Ds = geom_traits().compute_squared_distance_2_object()(query, ps);
      FT Dt = geom_traits().compute_squared_distance_2_object()(query, pt);
      FT dist2 = ((std::min))(Ds, Dt);

      FT norm2 = sample->distance2();
      FT tang2 = dist2 - norm2;

      sum.add(Cost_(norm2, tang2), mass);
      sum.compute_max(norm2, tang2);
    }
    set_vertex_cost(edge, sum);
    set_vertex_cost(twin, sum);
  }

  // SAMPLE //

  template<class Iterator> // value_type = Sample_*
  void assign_samples(Iterator begin, Iterator end) {
    for (Iterator it = begin; it != end; ++it) {
      Sample_* sample = *it;
      assign_sample(sample);
    }
  }

  template<class Iterator> // value_type = Sample_*
  void assign_samples_brute_force(Iterator begin, Iterator end) {
    for (Iterator it = begin; it != end; ++it) {
      Sample_* sample = *it;
      assign_sample_brute_force(sample);
    }
  }

  bool assign_sample(Sample_* sample) {
    const Point& point = sample->point();
    Face_handle face = Base::locate(point);

    if (face == Face_handle() || Base::is_infinite(face)) {
      //std::cout << "free bird" << std::endl;
      return false;
    }

    Vertex_handle vertex = find_nearest_vertex(point, face);
    if (vertex != Vertex_handle()) {
      assign_sample_to_vertex(sample, vertex);
      return true;
    }

    Edge edge = find_nearest_edge(point, face);
    assign_sample_to_edge(sample, edge);
    return true;
  }

  bool assign_sample_brute_force(Sample_* sample) {
    const Point& point = sample->point();
    Face_handle nearest_face = Face_handle();
    for (Finite_faces_iterator fi = Base::finite_faces_begin();
        fi != Base::finite_faces_end(); ++fi) {
      Face_handle face = fi;
      if (face_has_point(face, point)) {
        nearest_face = face;
        break;
      }
    }

    if (nearest_face == Face_handle()) {
      //std::cout << "free bird" << std::endl;
      return false;
    }

    Vertex_handle vertex = find_nearest_vertex(point, nearest_face);
    if (vertex != Vertex_handle()) {
      assign_sample_to_vertex(sample, vertex);
      return true;
    }

    Edge edge = find_nearest_edge(point, nearest_face);
    assign_sample_to_edge(sample, edge);
    return true;
  }

  bool face_has_point(Face_handle face, const Point& query) const {
    for (int i = 0; i < 3; ++i) {
      Edge edge(face, i);
      const Point& ps = source_vertex(edge)->point();
      const Point& pt = target_vertex(edge)->point();
      if (!compute_triangle_ccw(ps, pt, query))
        return false;
    }
    return true;
  }

  Vertex_handle find_nearest_vertex(const Point& point, Face_handle face) const
  {
    for (int i = 0; i < 3; ++i) {
      Vertex_handle vi = face->vertex(i);
      const Point& pi = vi->point();
      if (pi == point)
        return vi;
    }
    return Vertex_handle();
  }

  Edge find_nearest_edge(const Point& point, Face_handle face) const {
    Edge nearest(face, 0);
    FT min_dist2 = compute_distance2(point, get_segment(nearest));
    for (int i = 1; i < 3; ++i) {
      Edge edge(face, i);
      Segment segment = get_segment(edge);
      FT dist2 = compute_distance2(point, segment);
      if (dist2 < min_dist2) {
        min_dist2 = dist2;
        nearest = edge;
      }
    }

    return nearest;
  }

  void assign_sample_to_vertex(Sample_* sample, Vertex_handle vertex) const {
    /*if (vertex->sample()) {
      std::cout << "assign to vertex: vertex already has sample"
          << std::endl;
    }*/

    sample->distance2() = FT(0);
    sample->coordinate() = FT(0);
    vertex->set_sample(sample);
  }

  void assign_sample_to_edge(Sample_* sample, const Edge& edge) const {
    Segment segment = get_segment(edge);
    const Point& query = sample->point();
    sample->distance2() = compute_distance2(query, segment);
    sample->coordinate() = compute_coordinate(query, segment);
    edge.first->add_sample(edge.second, sample);
  }

  FT compute_distance2(const Point& query, const Segment& segment) const {

    if (geom_traits().orientation_2_object()(segment.source(), segment.target(), query) == COLLINEAR)
      return FT(0);

    Line line = geom_traits().construct_line_2_object()(segment);
    return geom_traits().compute_squared_distance_2_object()(query, line);
  }

  FT compute_coordinate(const Point& q, const Segment& segment) const {
    const Point& p0 = segment.source();
    const Point& p1 = segment.target();
    Vector p0p1 = geom_traits().construct_vector_2_object()(p0, p1);
    Vector p0q  = geom_traits().construct_vector_2_object()(p0, q);
    
    FT t = geom_traits().compute_scalar_product_2_object()(p0q, p0p1)
         / geom_traits().compute_scalar_product_2_object()(p0p1, p0p1);
    return t; // [0,1]
  }

  // SIGNED DISTANCE //

  // signed distance from line(a,b) to point t
  FT signed_distance(Vertex_handle a, Vertex_handle b,
      Vertex_handle t) const {
    const Point& pa = a->point();
    const Point& pb = b->point();
    const Point& pt = t->point();
    return compute_signed_distance(pa, pb, pt);
  }

  // signed distance from line(a,b) to point t
  FT compute_signed_distance(
    const Point& pa, const Point& pb, const Point& pt) const
  {
    if (pt == pa)
      return FT(0);
    if (pt == pb)
      return FT(0);
    if (pa == pb)
      return CGAL::approximate_sqrt(geom_traits().compute_squared_distance_2_object()(pa, pt));

    Vector vab = geom_traits().construct_vector_2_object()(pa, pb);
    // Normalize vab
    vab = geom_traits().construct_scaled_vector_2_object()(
      vab,
      FT(1) / CGAL::approximate_sqrt(geom_traits().compute_squared_length_2_object()(vab)));
    Vector vab90 = geom_traits().construct_vector_2_object()(-vab.y(), vab.x());
    Vector vat = geom_traits().construct_vector_2_object()(pa, pt);
    return geom_traits().compute_scalar_product_2_object()(vat, vab90);
  }

  // signed distance from t to the intersection of line(a,b) and line(t,s)
  // the pair::first is false if sign==-1 and true otherwise
  std::pair<bool,boost::optional<FT> >
  signed_distance_from_intersection(Vertex_handle a, Vertex_handle b,
      Vertex_handle t, Vertex_handle s) const {
    const Point& pa = a->point();
    const Point& pb = b->point();
    const Point& pt = t->point();
    const Point& ps = s->point();
    return compute_signed_distance_from_intersection(pa, pb, pt, ps);
  }

  // signed distance from t to the intersection of line(a,b) and line(t,s)
  // the pair::first is false if sign==-1 and true otherwise
  std::pair<bool,boost::optional<FT> >
  compute_signed_distance_from_intersection(
    const Point& pa, const Point& pb, const Point& pt, const Point& ps) const
  {
    FT Dabt = compute_signed_distance(pa, pb, pt);
    if (Dabt == FT(0))
      return std::make_pair(true,boost::make_optional(FT(0)));

    Line lab = geom_traits().construct_line_2_object()(
      pa, geom_traits().construct_vector_2_object()(pa, pb));
    Line lts = geom_traits().construct_line_2_object()(
      pt, geom_traits().construct_vector_2_object()(pt, ps));

    boost::optional<FT> Dqt;
    typename CGAL::cpp11::result_of<typename Traits_::Intersect_2(Line, Line)>::type
      result = intersection(lab, lts);
    if (result)
    {
      const Point* iq = boost::get<Point>(&(*result));
      if (iq)
        Dqt = CGAL::approximate_sqrt(geom_traits().compute_squared_distance_2_object()(*iq, pt));
    }

    return std::make_pair( (Dabt < FT(0) ? false : true) ,Dqt );
  }

  bool is_triangle_ccw(Vertex_handle a, Vertex_handle b, Vertex_handle c) const
  {
    const Point& pa = a->point();
    const Point& pb = b->point();
    const Point& pc = c->point();
    return compute_triangle_ccw(pa, pb, pc);
  }

  bool compute_triangle_ccw(
    const Point& pa, const Point& pb, const Point& pc) const
  {
    return geom_traits().orientation_2_object()(pa,pb,pc) != RIGHT_TURN;
  }

  // COMBINATORIAL TESTS //

  // (a,b) is cyclic if (a,b,c) and (a,c,b) exist
  bool is_edge_cyclic(const Edge& edge) const {
    Vertex_handle f = opposite_vertex(edge);
    Vertex_handle b = opposite_vertex(twin_edge(edge));
    return (f == b);
  }

  // b from (a,b) is cyclic if (a,b,c) and (b,a,c) exist
  bool is_target_cyclic(const Edge& edge) const {
    if (!is_edge_cyclic(edge))
      return false;

    Edge twin = twin_edge(edge);
    Edge prev = prev_edge(twin);
    Face_handle fp = prev.first->neighbor(prev.second);
    Face_handle ft = twin.first->neighbor(twin.second);
    return (fp == ft);
  }

  bool is_flippable(const Edge& edge) const {
    Edge twin = twin_edge(edge);
    if (Base::is_infinite(twin.first))
      return false;
    if (Base::is_infinite(edge.first))
      return false;

    Vertex_handle vs = source_vertex(edge);
    Vertex_handle vt = target_vertex(edge);
    Vertex_handle vf = opposite_vertex(edge);
    Vertex_handle vb = opposite_vertex(twin);

    return is_triangle_ccw(vs, vb, vf) && is_triangle_ccw(vt, vf, vb);
  }

  bool is_collapsible(const Edge& edge) const {
    return check_link_test(edge) && check_kernel_test(edge);
  }

  bool check_link_test(const Edge& edge) const {
    Vertex_handle s = source_vertex(edge);
    Vertex_handle t = target_vertex(edge);

    if (s == t)
      return false;
    typename Vertex_handle_set::const_iterator it;

    Vertex_handle_set svertices;
    get_vertices_from_vertex_link(s, svertices);

    Vertex_handle_set tvertices;
    get_vertices_from_vertex_link(t, tvertices);

    // link(s) inter link(t)
    Vertex_handle_set ivertices;
    for (it = svertices.begin(); it != svertices.end(); ++it) {
      Vertex_handle v = *it;
      if (tvertices.find(v) != tvertices.end())
        ivertices.insert(v);
    }

    Vertex_handle_set evertices;
    get_vertices_from_edge_link(edge, evertices);

    // link(edge) =? link(s) inter link(t)
    if (evertices.size() != ivertices.size())
      return false;

    for (it = evertices.begin(); it != evertices.end(); ++it) {
      Vertex_handle v = *it;
      if (ivertices.find(v) == ivertices.end())
        return false;
    }
    return true;
  }

  bool check_kernel_test(const Edge& edge) const {
    Vertex_handle s = source_vertex(edge);
    Vertex_handle t = target_vertex(edge);

    Edge_vector hull;
    hull.reserve(16);
    get_edges_from_star_minus_link(s, hull);
    return is_in_kernel(t->point(), hull.begin(), hull.end());
  }

  template<class Iterator> // value_type = Edge
  bool is_in_kernel(const Point& query, Iterator begin, Iterator end) const {
    for (Iterator it = begin; it != end; ++it) {
      Edge edge = *it;
      const Point& pa = source_vertex(edge)->point();
      const Point& pb = target_vertex(edge)->point();
      if (!compute_triangle_ccw(pa, pb, query))
        return false;
    }
    return true;
  }

  // COLLAPSE //

  // (s,a,b) + (s,b,c) -> (s,a,c) + (a,b,c)
  // st = (source,target) from 'make_collapsible'
  // return (a,c)
  Edge flip(const Edge& sb, Edge& st, int /*verbose*/ = 0) {
    Vertex_handle t = target_vertex(st);

    Edge sc = twin_edge(prev_edge(sb));
    Base::tds().flip(sb.first, sb.second);
    Edge ac = prev_edge(twin_edge(sc));

    Vertex_handle a = source_vertex(ac);
    if (a == t)
      st = prev_edge(ac);

    return ac;
  }

  void collapse(const Edge& edge, int /*verbose*/ = 0) {
    if (is_edge_cyclic(edge)) {
      collapse_cyclic_edge(edge);
      return;
    }

    Edge twin = twin_edge(edge);
    Base::tds().join_vertices(twin);
  }

  // (a,b,c) + (c,b,a) + (a,c,i) + (c,a,j) ->
  // (a,c,i) + (c,a,j)
  void collapse_cyclic_edge(const Edge& bc, int verbose = 0) {
    if (verbose > 1)
      std::cout << "collapse_cyclic_edge ... ";

    Edge cb = twin_edge(bc);
    Face_handle abc = bc.first;
    Face_handle cba = cb.first;

    Vertex_handle b = source_vertex(bc);
    Vertex_handle c = target_vertex(bc);
    Vertex_handle a = opposite_vertex(bc);

    Edge ac = twin_edge(next_edge(bc));
    Edge ca = twin_edge(prev_edge(cb));

    a->set_face(ac.first);
    c->set_face(ca.first);
    ac.first->set_neighbor(ac.second, ca.first);
    ca.first->set_neighbor(ca.second, ac.first);

    this->delete_face(abc);
    this->delete_face(cba);
    this->delete_vertex(b);

    if (verbose > 1)
      std::cout << "done" << std::endl;
  }

  void print_edge(Rec_edge_2 edge) const {
    int i = ((edge).edge()).second;
    Point a = ((edge).edge()).first->vertex((i+1)%3)->point();
    Point b = ((edge).edge()).first->vertex((i+2)%3)->point();
    std::cout <<"( " << (edge).priority()  <<  ") ( " << a << " , " << b << " )" << std::endl;
  }

  bool is_p_infinity(const std::pair<bool,boost::optional<FT> >& p) const
  {
    return p.first && p.second==boost::none;
  }

  bool is_m_infinity(const std::pair<bool,boost::optional<FT> >& p) const
  {
    return !p.first && p.second==boost::none;
  }

  bool is_infinity(const std::pair<bool,boost::optional<FT> >& p) const
  {
    return p.second==boost::none;
  }

  std::pair<bool,boost::optional<FT> > m_infinity() const
  {
    return std::pair<bool,boost::optional<FT> >(false,boost::optional<FT>());
  }

  template <class Iterator> // value_type = Edge
  bool make_collapsible(Edge& edge, Iterator begin, Iterator end, int verbose = 0)
  {
    Vertex_handle source = source_vertex(edge);
    Vertex_handle target = target_vertex(edge);

    MultiIndex multi_ind;
    for (Iterator it = begin; it != end; ++it)
    {
      Edge ab = twin_edge(*it);
      Vertex_handle a = source_vertex(ab);
      Vertex_handle b = target_vertex(ab);
      std::pair<bool,boost::optional<FT> > D = signed_distance_from_intersection(a, b, target, source);
      if (!D.first ) {
        CGAL_assertion(D.second!=boost::none);
        multi_ind.insert(Rec_edge_2(ab, *D.second));
      }
    }

    int nb_flips = 0;
    while (!multi_ind.empty())
    {
      Rec_edge_2 pedge = *(multi_ind.template get<1>()).begin();
      FT Dbc = pedge.priority();
      Edge bc = pedge.edge();
      (multi_ind.template get<0>()).erase(pedge);

      Edge sb = prev_edge(bc);
      Edge ab = prev_edge(twin_edge(sb));
      Edge sc = twin_edge(next_edge(bc));
      Edge cd = next_edge(sc);

      Vertex_handle a = source_vertex(ab);
      Vertex_handle b = source_vertex(bc);
      Vertex_handle c = target_vertex(bc);
      Vertex_handle d = target_vertex(cd);

      std::pair<bool,boost::optional<FT> > Dac=m_infinity();
      if (a != c && is_triangle_ccw(a, b, c))
        Dac = signed_distance_from_intersection(a, c, target, source);

      std::pair<bool,boost::optional<FT> > Dbd=m_infinity();
      if (b != d && is_triangle_ccw(b, c, d))
        Dbd = signed_distance_from_intersection(b, d, target, source);

      if ( is_m_infinity(Dac) && is_m_infinity(Dbd) )
      {
        if (verbose > 1)
          std::cerr << "--- No flips available ---"  << std::endl;
        return false;
      }

      FT value = Dbc <= 0 ? 1 : 2*Dbc; // value used if Dbd or Dac are +infinity
      if ( !is_infinity(Dac) )
      {
        if ( !is_infinity(Dbd))
          value = (std::max)(*Dac.second * (Dac.first?1:-1),
                             *Dbd.second * (Dbd.first?1:-1) );
        else
          if ( !is_p_infinity(Dbd) )
            value = *Dac.second * (Dac.first?1:-1);
      }
      else
        if ( !is_infinity(Dbd) && !is_p_infinity(Dac))
          value = *Dbd.second * (Dbd.first?1:-1);

      // if ( max(Dac,Dbd)+CGAL_EPS < Dbc )
      if (value + CGAL_EPS < Dbc)
      {
			/*
        std::cerr.precision(10);
        std::cerr << "--- Flip makes kernel worse ---" << std::endl;
        std::cerr << Dac << " or " << Dbd << " vs " << Dbc << std::endl;
        std::cerr << "a: " << a->point() << std::endl;
        std::cerr << "b: " << b->point() << std::endl;
        std::cerr << "c: " << c->point() << std::endl;
        std::cerr << "d: " << d->point() << std::endl;
        std::cerr << "t: " << target->point() << std::endl;
        std::cerr << "diff = " << Dbc - (std::max)(Dac, Dbd) << std::endl;
				*/
        return false;
      }

      // if (Dac > Dbd)
      if (is_p_infinity(Dac) || is_m_infinity(Dbd) ||
          (!is_p_infinity(Dbd) && !is_m_infinity(Dac)
           && *Dac.second * (Dac.first?1:-1) > *Dbd.second * (Dbd.first?1:-1)))
      {
        (multi_ind.template get<0>()).erase(Rec_edge_2(ab));

        Edge ac = flip(sb, edge, verbose);
        if (!Dac.first) {
          multi_ind.insert(Rec_edge_2(ac, - *Dac.second));
        }
      }
      else
      {
        (multi_ind.template get<0>()).erase(Rec_edge_2(cd));
        Edge bd = flip(sc, edge, verbose);
        if (!Dbd.first) {
          multi_ind.insert(Rec_edge_2(bd, - *Dbd.second));
        }
      }
      nb_flips++;
    }

    if (verbose > 1)
      std::cerr  << "Nb flips: "  << nb_flips << std::endl;

    return true;
  }
};

} } //namespace CGAL

#undef CGAL_EPS

#endif // CGAL_RECONSTRUCTION_TRIANGULATION_2_H
