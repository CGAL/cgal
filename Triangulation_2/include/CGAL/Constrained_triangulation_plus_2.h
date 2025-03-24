// Copyright (c) 2012  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Andreas Fabri, Mariette Yvinec

#ifndef CGAL_CONSTRAINED_TRIANGULATION_PLUS_2_H
#define CGAL_CONSTRAINED_TRIANGULATION_PLUS_2_H

#include <CGAL/license/Triangulation_2.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/Unique_hash_map.h>
#include <CGAL/assertions.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Triangulation_2/internal/Polyline_constraint_hierarchy_2.h>
#include <CGAL/Triangulation_2/internal/CTP2_subconstraint_graph.h>

#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_2/insert_constraints.h>

#include <boost/container/flat_set.hpp>

#include <array>
#include <type_traits>
#include <vector>

namespace CGAL {

// Comparison functor that compares two Vertex_handle.
// Used as 'Compare' functor for the constraint hierarchy.
template < class Tr >
class Ctp_2_compare_vertex_handles {
  const Tr* tr_p;

public:
  Ctp_2_compare_vertex_handles(const Tr* tr_p) : tr_p(tr_p) {}

  using Vertex_handle = typename Tr::Vertex_handle;

  bool operator()(const Vertex_handle& va,
                  const Vertex_handle& vb) const
  {
    return tr_p->compare_xy(va->point(), vb->point()) == SMALLER;
  }
}; // end class template Ctp_2_compare_vertex_handles

template <class Tr>
using Ctp_2_point_type = typename Tr::Geom_traits::Point_2;

template <class Tr>
using Ctp_2_hierarchy_type =
    Polyline_constraint_hierarchy_2<typename Tr::Vertex_handle,
                                    Ctp_2_compare_vertex_handles<Tr>,
                                    Ctp_2_point_type<Tr>>;

// Tr the base triangulation class
// Tr has to be Constrained or Constrained_Delaunay with Constrained_triangulation_plus_vertex_base

template < class Tr_>
class Constrained_triangulation_plus_2
  : public Tr_
  , protected Ctp_2_hierarchy_type<Tr_>
{
public:
  using Self = Constrained_triangulation_plus_2<Tr_>;
  using Base = Tr_;
  using Constraint_hierarchy = Ctp_2_hierarchy_type<Tr_>;
protected:
  const auto& hierarchy() const { return static_cast<const Constraint_hierarchy&>(*this); }
  auto& hierarchy() { return static_cast<Constraint_hierarchy&>(*this); }

private:
  using Tr = Tr_;


  template<class CDT>
  class Face_container
  {
    using Vertex_handle = typename CDT::Vertex_handle;
    using Face_handle = typename CDT::Face_handle;

    using Array = std::array<Vertex_handle, 3>;
    std::vector<Array> faces;
    CDT& cdt;

  public:
    using value_type = Face_handle;
    Face_container(CDT& cdt_ ) : cdt(cdt_) {}

    void push_back(Face_handle fh) {
      faces.push_back(Array{fh->vertex(0), fh->vertex(1), fh->vertex(2)});
    }

    template <class OutputIterator>
    void
    write_faces(OutputIterator out)
    {
      for(auto [v0, v1, v2] : make_range(faces.rbegin(), faces.rend())) {
        Face_handle fh;
        if(cdt.is_face(v0, v1, v2, fh)) {
          *out++ = fh;
        }
      }
    }
  };

public:
  // type aliases (aka type defs)
  using Triangulation = Tr;
  using Intersection_tag = typename Tr::Intersection_tag;

  // using-declarations of types or member functions from the two bases
  using Triangulation::vertices_begin;
  using Triangulation::vertices_end;
  using Triangulation::is_infinite;
  using Triangulation::number_of_vertices;

  using typename Triangulation::Edge;
  using typename Triangulation::Vertex;
#if defined(BOOST_MSVC) && (BOOST_MSVC < 1920)
  using Vertex_handle = typename Triangulation::Vertex_handle; // workaround for VC++ 19.16 (from MSVC 2017)
#else
  using typename Triangulation::Vertex_handle;
#endif
  using typename Triangulation::Face_handle;
  using typename Triangulation::Face_circulator;
  using typename Triangulation::Vertex_iterator;
  using typename Triangulation::Vertex_circulator;
  using typename Triangulation::Locate_type;
  using typename Triangulation::Line_face_circulator;
  using typename Triangulation::Geom_traits;
  using typename Triangulation::Constraint;
  using typename Triangulation::size_type;
  using typename Triangulation::List_edges;
  using typename Triangulation::List_faces;
  using typename Triangulation::List_vertices;
  using typename Triangulation::List_constraints;
  using typename Triangulation::Constrained_edges_iterator;

  using typename Constraint_hierarchy::Context;
  using typename Constraint_hierarchy::Context_iterator;
  using typename Constraint_hierarchy::Contexts;
  using typename Constraint_hierarchy::Constraint_iterator;
  using typename Constraint_hierarchy::Constraints;
  using typename Constraint_hierarchy::Subconstraint_iterator;
  using typename Constraint_hierarchy::Subconstraints;
  using typename Constraint_hierarchy::Subconstraint_and_contexts_iterator;
  using typename Constraint_hierarchy::Subconstraints_and_contexts;
  using typename Constraint_hierarchy::Constraint_id;
  using typename Constraint_hierarchy::Vertex_handle_compare;

#ifdef CGAL_CDT_2_DEBUG_INTERSECTIONS
  using Triangulation::display_vertex;
#endif // CGAL_CDT_2_DEBUG_INTERSECTIONS


  using Point = typename Triangulation::Geom_traits::Point_2;
  using Segment = typename Triangulation::Geom_traits::Segment_2;

  // Tag to mark the presence of a hierarchy of constraints
  using Constraint_hierarchy_tag = Tag_true;

  //Tag to distinguish Delaunay from regular triangulations
  using Weighted_tag = Tag_false;

  // Tag to distinguish periodic triangulations from others
  using Periodic_tag = Tag_false;

  // for user interface with the constraint hierarchy
  using Vertices_in_constraint_iterator = typename Constraint_hierarchy::Vertex_it;

  using Vertices_in_constraint = Iterator_range<Vertices_in_constraint_iterator>;

  using Points_in_constraint_iterator = typename Constraint_hierarchy::Point_it;
  using Points_in_constraint = Iterator_range<Points_in_constraint_iterator>;


  using Subconstraint = std::pair<Vertex_handle, Vertex_handle>;

  using Triangulation::geom_traits;
  using Triangulation::cw;
  using Triangulation::ccw;
  using Triangulation::incident_faces;

public:
  Constraint_hierarchy& hierarchy_ref()
  {
    return *this;
  }

  Constrained_triangulation_plus_2(const Geom_traits& gt=Geom_traits())
    : Triangulation(gt)
    , Constraint_hierarchy(Vertex_handle_compare(this))
  { }

  Constrained_triangulation_plus_2(const Constrained_triangulation_plus_2& ctp)
    : Constrained_triangulation_plus_2(ctp.geom_traits())
  { copy_triangulation(ctp);}

  Constrained_triangulation_plus_2(Constrained_triangulation_plus_2&&) = default;

  ~Constrained_triangulation_plus_2() override {}

  Constrained_triangulation_plus_2 & operator=(const Constrained_triangulation_plus_2& ctp)
  {
    copy_triangulation(ctp);
    return *this;
  }

  Constrained_triangulation_plus_2& operator=(Constrained_triangulation_plus_2&&) = default;

  template<class InputIterator>
  Constrained_triangulation_plus_2(InputIterator first,
                                   InputIterator last,
                                   const Geom_traits& gt=Geom_traits() )
    : Constrained_triangulation_plus_2(gt)
  {
    insert_constraints(first, last);
    CGAL_postcondition( this->is_valid() );
  }


  Constrained_triangulation_plus_2(const std::list<std::pair<Point,Point> > &constraints,
                                   const Geom_traits& gt=Geom_traits() )
    : Constrained_triangulation_plus_2(gt)
  {
    insert_constraints(constraints.begin(), constraints.end());
    CGAL_postcondition( this->is_valid() );
  }
  //Helping
  void clear() { Base::clear(); hierarchy().clear(); }
  void copy_triangulation(const Constrained_triangulation_plus_2 &ctp);
  void swap(Constrained_triangulation_plus_2 &ctp);

  // INSERTION
  Vertex_handle insert(const Point& a,
                       Face_handle start = Face_handle() );
  Vertex_handle insert(const Point& p,
                       Locate_type lt,
                       Face_handle loc, int li );

  Constraint_id insert_constraint(const Point& a, const Point& b)
  {
    Vertex_handle va= insert(a);
    // If the segment is "short" it is a good idea to start the next insertion
    // close to point a
    // Otherwise, to start here is as good as elsewhere
    Vertex_handle vb = insert(b, va->face());
    return insert_constraint(va, vb);
  }

  Constraint_id insert_constraint(const Constraint& c)
  {
    return insert_constraint(c.first, c.second);
  }

  Constraint_id insert_constraint(Vertex_handle va, Vertex_handle vb)
  {
#ifdef CGAL_CDT_2_DEBUG_INTERSECTIONS
  std::cerr << CGAL::internal::cdt_2_indent_level
            << "CT_plus_2::insert_constraint( " << display_vertex(va)
            << " , " << display_vertex(vb)
            << " )\n";
#endif // CGAL_CDT_2_DEBUG_INTERSECTIONS
    // protects against inserting a zero length constraint
    if(va == vb){
      return Constraint_id(nullptr);
    }
    // protects against inserting twice the same constraint
    Constraint_id cid = hierarchy().insert_constraint_old_API(va, vb);
    if (va != vb && (cid != Constraint_id(nullptr)) )  insert_subconstraint(va,vb);

    return cid;
  }

  template < class InputIterator>
  Constraint_id insert_constraint(InputIterator first, InputIterator last, bool close=false)
  {
    return insert_constraint_seq_impl(first, last, close);
  }

  template<typename Range>
  Constraint_id insert_constraint(const Range& r)
  {
    return insert_constraint_seq_impl(r.begin(), r.end(), false);
  }

  template < class PolygonTraits_2, class Container>
  Constraint_id insert_constraint(const Polygon_2<PolygonTraits_2,Container>& polygon)
  {
    return insert_constraint_seq_impl(polygon.vertices_begin(), polygon.vertices_end(), true);
  }
  /*
  template<typename InputIterator>
  size_type insert_constraints(InputIterator first, InputIterator last)
  {
    size_type n = 0;

    for(; first != last; ++first)
    {
      if(insert_constraint(*first))
        ++n;
    }
    return n;
  }
  */

  void split_subconstraint_graph_into_constraints(const std::function<bool(Vertex_handle)>& is_terminal
                                                  = std::function<bool(Vertex_handle)>())
  {
    internal::CTP2_graph_visitor<Self> visitor(*this);
    if (is_terminal)
      CGAL::split_graph_into_polylines (internal::CTP2_subconstraint_graph<Self>(*this), visitor,
                                        [&is_terminal](Vertex_handle vh,
                                                       const internal::CTP2_subconstraint_graph<Self>&) -> bool
                                        {
                                          return is_terminal(vh);
                                        });
    else
      CGAL::split_graph_into_polylines (internal::CTP2_subconstraint_graph<Self>(*this), visitor);
  }

  Vertex_handle push_back(const Point& p)
  {
    return insert(p);
  }

  Constraint_id push_back(const Constraint& c)
  {
    return insert_constraint(c.first, c.second);
  }

  // for backward compatibility
  // not const Point&, because otherwise VC6/7 messes it up with
  // the insert that takes an iterator range
  Constraint_id insert(Point a, Point b) { return insert_constraint(a, b); }
  Constraint_id insert(Vertex_handle va, Vertex_handle  vb) { return insert_constraint(va,vb); }



  template <class PointIterator, class IndicesIterator>
  std::size_t insert_constraints(PointIterator points_first,
                                 PointIterator points_beyond,
                                 IndicesIterator indices_first,
                                 IndicesIterator indices_beyond)
  {
    std::vector<Point> points(points_first, points_beyond);
    return internal::insert_constraints(*this,points, indices_first, indices_beyond);
  }


 template <class ConstraintIterator>
  std::size_t insert_constraints(ConstraintIterator first,
                                 ConstraintIterator beyond)
  {
    return internal::insert_constraints(*this,first,beyond);
  }


  Vertices_in_constraint_iterator
  insert_vertex_in_constraint(Constraint_id cid, Vertices_in_constraint_iterator pos,
                              Vertex_handle vh)
  {
    return insert_vertex_in_constraint(cid, pos, vh, Emptyset_iterator());
  }

  Vertices_in_constraint_iterator
  remove_vertex_from_constraint(Constraint_id cid, Vertices_in_constraint_iterator pos)
  {
    return remove_vertex_from_constraint(cid, pos, Emptyset_iterator());
  }



  // Removes pos from the constraint cid.
  // Returns the iterator to vertex that was just after pos (or end())
  // Writes the modified faces to out
  template <class OutputIterator>
  Vertices_in_constraint_iterator
  remove_vertex_from_constraint(Constraint_id cid, Vertices_in_constraint_iterator pos,
                                OutputIterator out)
  {
    if(pos == vertices_in_constraint_begin(cid)){
      // cid is [P, A, ..., B] -> split to aux=[P, A] and cid=[A...B]
      Constraint_id aux = hierarchy().split_head(cid, std::next(pos));
      remove_constraint(aux, out);
      return vertices_in_constraint_begin(cid);
    }

    if(pos == std::prev(vertices_in_constraint_end(cid))){
      // cid is [A, ..., B, P] -> split to cid=[A...B] and aux=[B,P]
      Constraint_id aux = hierarchy().split_tail(cid, std::prev(pos));
      remove_constraint(aux, out);
      return vertices_in_constraint_end(cid);
    }

    const auto second_vertex_of_cid = std::next(vertices_in_constraint_begin(cid));
    const auto next_to_last_vertex_of_cid = std::prev(vertices_in_constraint_end(cid), 2);

    Constraint_id head = nullptr, tail = nullptr;
    if(pos != second_vertex_of_cid){
      // cid is [A, ..., B, P, C, ..., D]
      // split to:
      //    head = [A...B] and,
      //    cid =  [B, P, C...D]
      // split off head
      head = hierarchy().split_head(cid, std::prev(pos));
    }
    if(pos != next_to_last_vertex_of_cid){
      // cid is now [B, P, C, ..., D]
      // split to:
      //    cid = [B, P, C] and,
      //    tail = [C...D]
      // split off tail
      tail = hierarchy().split_tail(cid,std::next(pos));
    }

    // now:
    //   cid is [B, P, C]
    //   head is null or [A...B]
    //   tail is null or [C...D]
    // Let create insert [C,D] and conditionally concatenate head and tail,
    // and return the iterator to C

    Vertex_handle b = *std::prev(pos);
    Vertex_handle c = *std::next(pos);

    Face_container<Constrained_triangulation_plus_2> fc(*this);
    Constraint_id bc = insert_constraint(b, c, std::back_inserter(fc));

    auto pos_before_c = std::prev(vertices_in_constraint_end(bc), 2);
    // `pos_before_c` is not necessarily == vertices_in_constraint_begin(bc)
    // there might have been intersecting constraints

    hierarchy().swap(cid, bc);
    remove_constraint(bc, std::back_inserter(fc)); // removes [B, P, C]
    // now cid is [B, C]

    if(head != nullptr){
      hierarchy().prepend(std::move(head), cid);
    }

    if(tail != nullptr){
      hierarchy().concatenate(cid, std::move(tail));
    }
    fc.write_faces(out);

    // we went one too far back because the last vertex `c` gets removed by concatenate/prepend
    return std::next(pos_before_c);
  }

  // Inserts vh before pos
  // Returns an iterator pointing on the newly inserted vertex
  // Writes the modified faces to out
  template <class OutputIterator>
  Vertices_in_constraint_iterator
  insert_vertex_in_constraint(Constraint_id cid, Vertices_in_constraint_iterator pos,
                              Vertex_handle v, OutputIterator out)
  {
    // Insertion before the first vertex
    if(pos == vertices_in_constraint_begin(cid)){
      //std::cout << "insertion before first vertex" << std::endl;
      Constraint_id head = insert_constraint(v, *pos, out);
      hierarchy().prepend(std::move(head), cid);
      return vertices_in_constraint_begin(cid);
    }

    // Insertion after the last vertex
    if(pos == vertices_in_constraint_end(cid)){
      //std::cout << "insertion after last vertex" << std::endl;
      Constraint_id tail = insert_constraint(*std::prev(pos), v, out);
      auto new_pos = std::prev(vertices_in_constraint_end(tail));
      hierarchy().concatenate(cid, std::move(tail));
      CGAL_assertion(v == *new_pos);
      return new_pos;
    }
    Vertices_in_constraint_iterator pred = std::prev(pos);
    Vertices_in_constraint_iterator latest_vertex = std::prev(vertices_in_constraint_end(cid));
    Vertex_handle a = *pred;
    Vertex_handle b = *pos;
    if(v == a || v == b){
      return pos;
    }

    // cid is [..., A, B, ...] and V=*v will be inserted between A and B

    Face_container<Constrained_triangulation_plus_2> fc(*this);
    Constraint_id a_v_b = insert_constraint(a, v, std::back_inserter(fc));
    Constraint_id aux = insert_constraint(v, b, std::back_inserter(fc));
    auto new_pos = vertices_in_constraint_begin(aux);
    concatenate(a_v_b, std::move(aux));

    // here:
    //   a_v_b is [A,.. V,.. B]
    //   aux is empty
    // and new_pos is the iterator to V
    CGAL_assertion(v == *new_pos);

    CGAL_assertion(std::distance(vertices_in_constraint_begin(a_v_b), new_pos) > 0 &&
                   std::distance(new_pos,   vertices_in_constraint_end(a_v_b)) > 0);
    // new_pos still points to something in a_v_b. In general a_v_b should only have three vertices,
    // but there might have been intersectiong constraints or vertices.

    const auto second_vertex_of_cid = std::next(vertices_in_constraint_begin(cid));
    // If the constraint consists only of a segment, and we want to insert
    // in the middle: cid is just the segment [A, B]
    if((pos == second_vertex_of_cid) && (second_vertex_of_cid == latest_vertex)){
      //std::cout << "insertion in constraint which is a segment" << std::endl;
      hierarchy().swap(cid, a_v_b);
      remove_constraint(a_v_b, std::back_inserter(fc));
      fc.write_faces(out);
      return new_pos;
    }

    Constraint_id head = nullptr, tail = nullptr;
    if(pos != second_vertex_of_cid){
      //std::cout << "split head" << std::endl;
      head = hierarchy().split_head(cid, pred);
      pred = vertices_in_constraint_begin(cid);
      pos = std::next(pred);
    }
    // head is now [..., A] or null
    //  cid is now [A, B, ...]
    CGAL_assertion(*pred == a);
    CGAL_assertion(*pos == b);
    if(pos != latest_vertex){
      //std::cout << "split tail" << std::endl;
      tail = hierarchy().split_tail(cid, pos);
    }
    // head is now [..., A] or null
    //  cid is now [A, B]
    // tail is now [B, ...] or null

    if(head != nullptr){
      //std::cout << "concatenate head" << std::endl;
      hierarchy().concatenate(head, std::move(a_v_b));
      hierarchy().swap(cid, head);
      remove_constraint(head, std::back_inserter(fc));
    } else {
      hierarchy().swap(cid, a_v_b);
      remove_constraint(a_v_b, std::back_inserter(fc));
    }
    // cid is now [..., A, V, B]
    // head is now null empty
    // a_v_b is now empty

    if(tail != nullptr){
      //std::cout << "concatenate tail" << std::endl;
      concatenate(cid, std::move(tail));
    }
    fc.write_faces(out);
    return new_pos;
  }

  template < class InputIterator, class OutputIterator>
  Constraint_id insert_constraint(InputIterator first, InputIterator last, OutputIterator out)
  {
    Face_handle hint;
    Face_container<Constrained_triangulation_plus_2> fc(*this);
    std::vector<Vertex_handle> vertices;
    for(;first!= last; first++){
      Vertex_handle vh = insert(*first, hint);
      hint = vh->face();
      // no duplicates
      if(vertices.empty() || (vertices.back() != vh)){
        vertices.push_back(vh);
      }
    }
    int n = vertices.size();
    if(n == 1){
      return nullptr;
    }
    Constraint_id ca = hierarchy().insert_constraint(vertices[0],vertices[1]);
    insert_subconstraint(vertices[0],vertices[1], std::back_inserter(fc));

    if(n>2){
      for(int j=1; j<n-1; j++){
        hierarchy().append_constraint(ca, vertices[j], vertices[j+1]);
        insert_subconstraint(vertices[j], vertices[j+1], std::back_inserter(fc));
      }
    }
    for(auto vh : vertices_in_constraint(ca)) {
      out = insert_incident_faces(vh, out);
    }
    //AF    vertices_in_constraint_begin(ca)->fixed() = true;
    // Vertices_in_constraint_iterator end = std::prev(vertices_in_constraint_end(ca));
    // end->fixed() = true;
    fc.write_faces(out);

    return ca;
  }


private:
  template < class InputIterator>
  Constraint_id insert_constraint_seq_impl(InputIterator first, InputIterator last, bool is_polygon)
  {
    Face_handle hint;
    std::vector<Vertex_handle> vertices;
    for(;first!= last; first++){
      Vertex_handle vh = insert(*first, hint);
      hint = vh->face();
      // no duplicates
      if(vertices.empty() || (vertices.back() != vh)){
        vertices.push_back(vh);
      }
    }
    if(is_polygon && (vertices.size()>1) && (vertices.front() != vertices.back())){
      vertices.push_back(vertices.front());
    }

    std::size_t n = vertices.size();
    if(n == 1){
      return Constraint_id{};
    }
    CGAL_assertion(n >= 2);

    Constraint_id ca = hierarchy().insert_constraint(vertices[0],vertices[1]);
    insert_subconstraint(vertices[0],vertices[1]);

    if(n>2){
      for(std::size_t j=1; j<n-1; j++){
        hierarchy().append_constraint(ca, vertices[j], vertices[j+1]);
        insert_subconstraint(vertices[j], vertices[j+1]);
      }
    }

    // fix first and last, one is redundant for is_polygon == true
    // vertices.front()->fixed() = true;
    // vertices.back()->fixed() = true;

    return ca;
  }

public:

  auto& file_output(std::ostream& os) const {
    os << static_cast<const Tr&>(*this);
    Unique_hash_map<Vertex_handle,int> V(0, number_of_vertices());
    int inum = 0;
    for(auto vh : this->finite_vertex_handles()) {
      V[vh] = inum++;
    }
    for(auto cid : constraints()) {
      os << cid.size();
      for(Vertex_handle vh : vertices_in_constraint(cid)) {
        os << " " << V[vh];
      }
      os << std::endl;
    }
    return os;
  }

  friend std::ostream& operator<<(std::ostream& os, const Constrained_triangulation_plus_2& ctp) {
    return ctp.file_output(os);
  }

  auto& file_input(std::istream& is) {
    is >> static_cast<Tr&>(*this);

    std::vector<Vertex_handle> vertices(number_of_vertices());
    auto [first, last] = this->finite_vertex_handles();
    std::copy(first, last, vertices.data());

    size_type n, id, id_next;
    while(is >> n) {
      is >> id >> id_next;
      Constraint_id cid = insert_constraint(vertices[id], vertices[id_next]);

      for(size_type i = 2; i < n; ++i) {
        id = id_next;
        is >> id_next;
        Constraint_id cid2 = insert_constraint(vertices[id], vertices[id_next]);
        cid = concatenate(cid, std::move(cid2));
      }
    }
    return is;
  }

  friend std::istream& operator>>(std::istream& is, Constrained_triangulation_plus_2& ctp) {
    return ctp.file_input(is);
  }

  template <class OutputIterator>
  typename Constrained_triangulation_plus_2<Tr>::Constraint_id
  insert_constraint(Vertex_handle va, Vertex_handle vb, OutputIterator out)
  {
    // protects against inserting a zero length constraint
    if(va == vb){
      return Constraint_id();
    }
    // protects against inserting twice the same constraint
    Constraint_id cid = hierarchy().insert_constraint(va, vb);
    if (va != vb && (cid != nullptr) )  insert_subconstraint(va,vb,out);

    for(auto vh : vertices_in_constraint(cid)) {
      out = insert_incident_faces(vh, out);
    }
    return cid;
  }

  Vertex_handle intersect(Face_handle f, int i,
                          Vertex_handle vaa,
                          Vertex_handle vbb) override;
  Vertex_handle intersect(Face_handle f, int i,
                          Vertex_handle vaa,
                          Vertex_handle vbb,
                          No_constraint_intersection_tag);
  Vertex_handle intersect(Face_handle f, int i,
                          Vertex_handle vaa,
                          Vertex_handle vbb,
                          No_constraint_intersection_requiring_constructions_tag);
  Vertex_handle intersect(Face_handle f, int i,
                          Vertex_handle vaa,
                          Vertex_handle vbb,
                          Exact_intersections_tag);
  Vertex_handle intersect(Face_handle f, int i,
                          Vertex_handle vaa,
                          Vertex_handle vbb,
                          Exact_predicates_tag);

  // REMOVAL

  template <class OutputIterator>
  void remove_constraint(Constraint_id cid, OutputIterator out)
  {
    std::vector<Vertex_handle> vertices(hierarchy().vertices_in_constraint_begin(cid),
                                        hierarchy().vertices_in_constraint_end(cid));

    hierarchy().remove_constraint(cid);
    for(auto it = vertices.begin(), succ = it; ++succ != vertices.end(); ++it){
      if(! is_subconstraint(*it, *succ)){ // this checks whether other constraints pass
        Face_handle fh;
        int i = -1;
        Triangulation::is_edge(*it, *succ, fh, i);
        CGAL_assertion(i != -1);
        Triangulation::remove_constrained_edge(fh,i, out); // this does also flipping if necessary.
      }
    }
  }
  void remove_constraint(Constraint_id cid)
  {
    remove_constraint(cid, Emptyset_iterator());
  }


  void simplify(Vertices_in_constraint_iterator v)
  {
    Vertices_in_constraint_iterator u = std::prev(v);
    Vertices_in_constraint_iterator w = std::next(v);
    bool unew = (*u != *w);
    hierarchy().simplify(u,v,w);

    Triangulation::remove_incident_constraints(*v);

    Triangulation::remove(*v);

    if(unew){
      Triangulation::insert_constraint(*u, *w);
    }
  }

  using Constraint_hierarchy::remove_points_without_corresponding_vertex;

  // CONCATENATE AND SPLIT

  // concatenates two constraints
  using Constraint_hierarchy::concatenate;

  // split a constraint in two constraints, so that vcit becomes the first
  // vertex of the new constraint
  // returns the new constraint
  Constraint_id
  split(Constraint_id first, Vertices_in_constraint_iterator vcit) {
    return hierarchy().split_tail(first, vcit);
  }

  // Query of the constraint hierarchy

  using Constraint_hierarchy::constraints_begin;
  using Constraint_hierarchy::constraints_end;
  using Constraint_hierarchy::constraints;
  using Constraint_hierarchy::subconstraints_begin;
  using Constraint_hierarchy::subconstraints_end;
  using Constraint_hierarchy::subconstraints;
  using Constraint_hierarchy::subconstraints_and_contexts_begin;
  using Constraint_hierarchy::subconstraints_and_contexts_end;
  using Constraint_hierarchy::subconstraints_and_contexts;
  using Constraint_hierarchy::context;
  using Constraint_hierarchy::number_of_enclosing_constraints;
  using Constraint_hierarchy::is_subconstraint;
  using Constraint_hierarchy::contexts_begin;
  using Constraint_hierarchy::contexts_end;
  using Constraint_hierarchy::contexts;
  using Constraint_hierarchy::vertices_in_constraint_begin;
  using Constraint_hierarchy::vertices_in_constraint_end;
  using Constraint_hierarchy::vertices_in_constraint;
  using Constraint_hierarchy::points_in_constraint_begin;
  using Constraint_hierarchy::points_in_constraint_end;

  Points_in_constraint points_in_constraint(Constraint_id cid) const
  {
    return Points_in_constraint(points_in_constraint_begin(cid), points_in_constraint_end(cid));
  }

  using Constraint_hierarchy::number_of_constraints;
  using Constraint_hierarchy::number_of_subconstraints;
  using Constraint_hierarchy::split_constraint;

protected:
  template <class OutputItertator>
  OutputItertator insert_incident_faces(Vertex_handle vh, OutputItertator out)
  {
    Face_circulator fc = incident_faces(vh), done = fc;
    if(fc != nullptr) {
      do {
        Face_handle fh = fc;
        *out++ = fh;
      } while(++fc != done);
    }
    return out;
  }

void
insert_subconstraint(Vertex_handle vaa,
                     Vertex_handle vbb)
{
  insert_subconstraint(vaa,vbb,Emptyset_iterator());
}




template <class OutputItertator>
void
insert_subconstraint(Vertex_handle vaa,
                     Vertex_handle vbb,
                     OutputItertator out)
  // insert the subconstraint [vaa vbb]
  // it will eventually be split into several subconstraints
{
#ifdef CGAL_CDT_2_DEBUG_INTERSECTIONS
  std::cerr << CGAL::internal::cdt_2_indent_level
            << "CT_plus_2::insert_subconstraint( " << display_vertex(vaa)
            << " , " << display_vertex(vbb)
            << " )\n";
  internal::Indentation_level::Exit_guard exit_guard = CGAL::internal::cdt_2_indent_level.open_new_scope();
  std::cerr << CGAL::internal::cdt_2_indent_level
            << "CT_plus_2::insert_constraint stack push [va, vb] ( " << display_vertex(vaa)
            << " , " << display_vertex(vbb)
            << " )\n";
#endif // CGAL_CDT_2_DEBUG_INTERSECTIONS
  std::stack<std::pair<Vertex_handle, Vertex_handle> > stack;
  stack.push(std::make_pair(vaa,vbb));

  while(! stack.empty()){
    auto [vaa, vbb] = stack.top();
    stack.pop();
    CGAL_precondition( vaa != vbb);
#ifdef CGAL_CDT_2_DEBUG_INTERSECTIONS
    std::cerr << CGAL::internal::cdt_2_indent_level
              << "CT_plus_2::insert_subconstraint, stack pop=( " << display_vertex(vaa)
              << " , " << display_vertex(vbb)
              << " ) remaining stack size: "
              << stack.size() << '\n';
    CGAL_assertion(this->is_valid());
#endif // CGAL_CDT_2_DEBUG_INTERSECTIONS
    Vertex_handle vi;

    Face_handle fr;
    int i;
    if(this->includes_edge(vaa,vbb,vi,fr,i)) {
#ifdef CGAL_CDT_2_DEBUG_INTERSECTIONS
    std::cerr << CGAL::internal::cdt_2_indent_level
              << "CT_plus_2::insert_subconstraint, the segment ( " << display_vertex(vaa)
              << " , " << display_vertex(vbb)
              << " ) is an edge with #"
              << vi->time_stamp() << "= " << vi->point()
              << '\n';
#endif // CGAL_CDT_2_DEBUG_INTERSECTIONS
      this->mark_constraint(fr,i);
      if (vi != vbb)  {
        hierarchy().split_constraint(vaa,vbb,vi);
#ifdef CGAL_CDT_2_DEBUG_INTERSECTIONS
  std::cerr << CGAL::internal::cdt_2_indent_level
            << "CT_plus_2::insert_constraint (includes_edge) stack push [vi, vbb] ( " << display_vertex(vi)
            << " , " << display_vertex(vbb)
            << " )\n";
#endif // CGAL_CDT_2_DEBUG_INTERSECTIONS
        stack.push(std::make_pair(vi,vbb));
      }
      continue;
    }

    List_faces intersected_faces;
    List_edges conflict_boundary_ab, conflict_boundary_ba;

    bool intersection  = this->find_intersected_faces(
                                                      vaa, vbb,
                                                      intersected_faces,
                                                      conflict_boundary_ab,
                                                      conflict_boundary_ba,
                                                      vi);

    if ( intersection) {
      if (vi != vaa && vi != vbb) {
        hierarchy().split_constraint(vaa,vbb,vi);
#ifdef CGAL_CDT_2_DEBUG_INTERSECTIONS
  std::cerr << CGAL::internal::cdt_2_indent_level
            << "CT_plus_2::insert_constraint stack push [vaa, vi] ( " << display_vertex(vaa)
            << " , " << display_vertex(vi)
            << " )\n";
  std::cerr << CGAL::internal::cdt_2_indent_level
            << "CT_plus_2::insert_constraint stack push [vi, vbb] ( " << display_vertex(vi)
            << " , " << display_vertex(vbb)
            << " )\n";
#endif // CGAL_CDT_2_DEBUG_INTERSECTIONS
        stack.push(std::make_pair(vaa,vi));
        stack.push(std::make_pair(vi,vbb));
      }
      else {
#ifdef CGAL_CDT_2_DEBUG_INTERSECTIONS
  std::cerr << CGAL::internal::cdt_2_indent_level
            << "CT_plus_2::insert_constraint stack push [vaa, vbb]( " << display_vertex(vaa)
            << " , " << display_vertex(vbb)
            << " )\n";
#endif // CGAL_CDT_2_DEBUG_INTERSECTIONS
        stack.push(std::make_pair(vaa,vbb));
      }

      continue;
    }


    //no intersection

    List_edges edges(conflict_boundary_ab);
    std::copy(conflict_boundary_ba.begin(), conflict_boundary_ba.end(), std::back_inserter(edges));

    // edges may contain mirror edges. They no longer exist after triangulate_hole
    // so we have to remove them before calling get_bounded_faces
    if(! edges.empty()){
      boost::container::flat_set<Face_handle> faces(intersected_faces.begin(), intersected_faces.end());
      for(typename List_edges::iterator it = edges.begin(); it!= edges.end();){
        if(faces.find(it->first) != faces.end()){
          typename List_edges::iterator it2 = it;
          ++it;
          edges.erase(it2);
        }else {
          ++it;
        }
      }
    }

    this->triangulate_hole(intersected_faces,
                           conflict_boundary_ab,
                           conflict_boundary_ba);

    this->get_bounded_faces(edges.begin(),
                            edges.end(),
                            out);

    if (vi != vbb) {
      hierarchy().split_constraint(vaa,vbb,vi);
      stack.push(std::make_pair(vi,vbb));
    }
  }
}



  //to debug
public:
  void print_hierarchy(std::ostream& os = std::cout) { hierarchy().print(os); }

  //template member functions
public:
  template < class InputIterator >
#if defined(_MSC_VER)
  std::ptrdiff_t insert(InputIterator first, InputIterator last, int i = 0)
#else
    std::ptrdiff_t insert(InputIterator first, InputIterator last)
#endif
  {
#if defined(_MSC_VER)
    CGAL_USE(i);
#endif
    size_type n = this->number_of_vertices();

    std::vector<Point> points (first, last);

    spatial_sort (points.begin(), points.end(), geom_traits());

    Face_handle hint;
    for (typename std::vector<Point>::const_iterator p = points.begin(), end = points.end();
            p != end; ++p)
        hint = insert (*p, hint)->face();

    return this->number_of_vertices() - n;
  }

};

template <class Tr>
void
Constrained_triangulation_plus_2<Tr>::
copy_triangulation(const Constrained_triangulation_plus_2 &ctp)
{
  Base::copy_triangulation(ctp);
  //the following assumes that the triangulation and its copy
  // iterate on their vertices in the same order
  std::map<Vertex_handle,Vertex_handle> vmap;
  Vertex_iterator vit = ctp.vertices_begin();
  Vertex_iterator vvit = this->vertices_begin();
  for( ; vit != ctp.vertices_end(); ++vit, ++vvit) {
    CGAL_assertion(vit->point() == vvit->point());
    vmap[vit] = vvit;
  }
  hierarchy().copy(ctp.hierarchy(), vmap);
}

template <class Tr>
void
Constrained_triangulation_plus_2<Tr>::
swap(Constrained_triangulation_plus_2 &ctp)
{
  Base::swap(ctp);
  hierarchy().swap(ctp.hierarchy());
}

template < class Tr >
inline
typename Constrained_triangulation_plus_2<Tr>::Vertex_handle
Constrained_triangulation_plus_2<Tr>::
insert(const Point& a, Face_handle start)
{
  Locate_type lt;
  int li;
  Face_handle loc = this->locate(a, lt, li, start);
  return insert(a,lt,loc,li);
}

template < class Tr>
typename Constrained_triangulation_plus_2<Tr>::Vertex_handle
Constrained_triangulation_plus_2<Tr>::
insert(const Point& a, Locate_type lt, Face_handle loc, int li)
{
  Vertex_handle v1, v2;
  bool insert_in_constrained_edge = false;

  if ( lt == Triangulation::EDGE && loc->is_constrained(li) )
  {
    if(std::is_same<typename Tr::Itag, No_constraint_intersection_tag>::value)
      throw typename Tr::Intersection_of_constraints_exception();

    insert_in_constrained_edge = true;
    v1=loc->vertex(ccw(li)); //endpoint of the constraint
    v2=loc->vertex(cw(li)); // endpoint of the constraint
  }

  Vertex_handle va = Triangulation::insert(a,lt,loc,li);
  // update the hierarchy
  if (insert_in_constrained_edge) {
#ifdef CGAL_CDT_2_DEBUG_INTERSECTIONS
    std::cerr << CGAL::internal::cdt_2_indent_level
              << "  CT_plus_2::insert(" << a << ") = #"
              << va->time_stamp()
              << "   insert in constrained edge:  #" << v1->time_stamp() << "= " << v1->point()
              << " , #" << v2->time_stamp() << "= " << v2->point()
              << std::endl;
#endif
    hierarchy().split_constraint(v1,v2,va);
  }
  return va;
}

template <class Tr>
typename Constrained_triangulation_plus_2<Tr>:: Vertex_handle
Constrained_triangulation_plus_2<Tr>::
intersect(Face_handle f, int i,
          Vertex_handle vaa,
          Vertex_handle vbb)
{
  return intersect(f, i, vaa, vbb, Intersection_tag());
}

template <class Tr>
typename Constrained_triangulation_plus_2<Tr>:: Vertex_handle
Constrained_triangulation_plus_2<Tr>::
intersect(Face_handle, int,
          Vertex_handle,
          Vertex_handle,
          No_constraint_intersection_tag)
{
  throw typename Tr::Intersection_of_constraints_exception();
  return Vertex_handle();
}

template <class Tr>
typename Constrained_triangulation_plus_2<Tr>:: Vertex_handle
Constrained_triangulation_plus_2<Tr>::
intersect(Face_handle, int,
          Vertex_handle,
          Vertex_handle,
          No_constraint_intersection_requiring_constructions_tag)
{
  throw typename Tr::Intersection_of_constraints_exception();
  return Vertex_handle();
}

template <class Tr>
typename Constrained_triangulation_plus_2<Tr>:: Vertex_handle
Constrained_triangulation_plus_2<Tr>::
intersect(Face_handle f, int i,
          Vertex_handle vaa,
          Vertex_handle vbb,
          Exact_intersections_tag)
// compute the intersection of the constraint edge (f,i)
// with the subconstraint (vaa,vbb) being inserted
// insert the intersection point
// (the  constraint edge (f,i) will be split in hierarchy by insert)
// and return the Vertex_handle of the new Vertex
{
  const Vertex_handle vcc = f->vertex(cw(i));
  const Vertex_handle vdd = f->vertex(ccw(i));
  const auto [vc, vd] = hierarchy().enclosing_constraint(vcc, vdd);
  const auto [va, vb] = hierarchy().enclosing_constraint(vaa, vbb);
  CGAL_assertion(vc != vd);
  CGAL_assertion(va != vb);

  const Point& pa = va->point();
  const Point& pb = vb->point();
  const Point& pc = vc->point();
  const Point& pd = vd->point();
#ifdef CGAL_CDT_2_DEBUG_INTERSECTIONS
  std::cerr << CGAL::internal::cdt_2_indent_level
            << "CT_plus_2::intersect segment ( " << display_vertex(va)
            << " , " << display_vertex(vb)
            << " ) with edge ( #"<< vc->time_stamp() << "= " << vc->point()
            << " , " << display_vertex(vd)
            << " , Exact_intersections_tag)\n";
#endif // CGAL_CDT_2_DEBUG_INTERSECTIONS
  Point pi(ORIGIN); // initialize although we are sure that it will be
                    // set by the intersection, but to quiet a warning
  Intersection_tag itag = Intersection_tag();
  CGAL_assertion_code( bool ok = )
  intersection(geom_traits(), pa, pb, pc, pd, pi, itag );
  CGAL_assertion(ok);

  Vertex_handle vi = insert(pi, Triangulation::EDGE, f, i);
#ifdef CGAL_CDT_2_DEBUG_INTERSECTIONS
  std::cerr << CGAL::internal::cdt_2_indent_level
            << "CT_plus_2::intersect, `vi` is ( " << display_vertex(vi)
            << " )\n";
#endif // CGAL_CDT_2_DEBUG_INTERSECTIONS
  return vi;
}

template <class Tr>
typename Constrained_triangulation_plus_2<Tr>::Vertex_handle
Constrained_triangulation_plus_2<Tr>::
intersect(Face_handle f, int i,
          Vertex_handle vaa,
          Vertex_handle vbb,
          Exact_predicates_tag itag)
{
  Vertex_handle  vcc, vdd;
  vcc = f->vertex(cw(i));
  vdd = f->vertex(ccw(i));

  const Point& pa = vaa->point();
  const Point& pb = vbb->point();
  const Point& pc = vcc->point();
  const Point& pd = vdd->point();
#ifdef CGAL_CDT_2_DEBUG_INTERSECTIONS
  std::cerr << CGAL::internal::cdt_2_indent_level
            << "CT_plus_2::intersect segment ( " << display_vertex(vaa)
            << " , " << display_vertex(vbb)
            << " ) with edge ( #"<< vcc->time_stamp() << "= " << vcc->point()
            << " , " << display_vertex(vdd)
            << " , Exact_predicates_tag)\n";
#endif // CGAL_CDT_2_DEBUG_INTERSECTIONS

  Vertex_handle vi = Triangulation::insert_intersection(
      f, i, vaa, vbb, vcc, vdd, pa, pb, pc, pd, itag);

  // vi == vc or vi == vd may happen even if intersection==true
  // due to approximate construction of the intersection
  if (vi != vcc && vi != vdd) {
    hierarchy().split_constraint(vcc,vdd,vi);
    insert_subconstraint(vcc,vi);
    insert_subconstraint(vi, vdd);
  }
  else {
    insert_subconstraint(vcc,vdd);
  }
  return vi;
}

} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif //CGAL_CONSTRAINED_TRIANGULATION_PLUS_2_H
