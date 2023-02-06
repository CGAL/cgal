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
#include <boost/tuple/tuple.hpp>

#include <CGAL/Default.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_2/insert_constraints.h>
#include <boost/container/flat_set.hpp>

#include <type_traits>

namespace CGAL {

// Comparison functor that compares two Vertex_handle.
// Used as 'Compare' functor for the constraint hierarchy.
template < class Tr >
class Pct2_vertex_handle_less_xy {
  const Tr* tr_p;

public:
  Pct2_vertex_handle_less_xy(const Tr* tr_p) : tr_p(tr_p) {}

  typedef typename Tr::Vertex_handle Vertex_handle;

  bool operator()(const Vertex_handle& va,
                  const Vertex_handle& vb) const
  {
    return tr_p->compare_xy(va->point(), vb->point()) == SMALLER;
  }
}; // end class template Pct2_vertex_handle_less_xy

// Tr the base triangulation class
// Tr has to be Constrained or Constrained_Delaunay with Constrained_triangulation_plus_vertex_base

template < class Tr_ = Default >
class Constrained_triangulation_plus_2
  : public
Default::Get< Tr_, Constrained_Delaunay_triangulation_2<
                      Exact_predicates_inexact_constructions_kernel
                      , Triangulation_data_structure_2<
                            Triangulation_vertex_base_2<Exact_predicates_inexact_constructions_kernel>
                          , Constrained_triangulation_face_base_2<Exact_predicates_inexact_constructions_kernel>
                          >
                      , CGAL::Exact_predicates_tag
                      > >::type
{
  typedef typename
  Default::Get< Tr_, Constrained_Delaunay_triangulation_2<
                  Exact_predicates_inexact_constructions_kernel
                  , Triangulation_data_structure_2<
                        Triangulation_vertex_base_2<Exact_predicates_inexact_constructions_kernel>
                      , Constrained_triangulation_face_base_2<Exact_predicates_inexact_constructions_kernel>
                      >
                  , CGAL::Exact_predicates_tag
                  > >::type Tr;


  template<class CDT>
  class Face_container
  {
    typedef typename CDT::Vertex_handle Vertex_handle;
    typedef typename CDT::Face_handle Face_handle;
  private:
    typedef boost::tuple<Vertex_handle, Vertex_handle, Vertex_handle> TFace;
    std::vector<TFace> faces;
    CDT& cdt;

  public:
    typedef Face_handle value_type;
    typedef Face_handle& reference;
    typedef const Face_handle& const_reference;

    Face_container(CDT& cdt_ ) : cdt(cdt_) {}

    void push_back(Face_handle fh)
    {
      faces.push_back(boost::make_tuple(fh->vertex(0),
                                        fh->vertex(1),
                                        fh->vertex(2)));
    }

    template <class OutputIterator>
    void
    write_faces(OutputIterator out)
    {
      for(typename std::vector<TFace>::reverse_iterator
            it = faces.rbegin(); it != faces.rend(); ++it) {
        Face_handle fh;
        if(cdt.is_face(boost::get<0>(*it), boost::get<1>(*it), boost::get<2>(*it), fh)){
          *out++ = fh;
        }
      }
    }
  };

public:
  typedef Tr                                   Triangulation;
  typedef typename Tr::Intersection_tag        Intersection_tag;
  typedef Constrained_triangulation_plus_2<Tr_> Self;
  typedef Tr                                   Base;


#ifndef CGAL_CFG_USING_BASE_MEMBER_BUG_2
  using Triangulation::vertices_begin;
  using Triangulation::vertices_end;
  using Triangulation::is_infinite;
  using Triangulation::number_of_vertices;
#endif
#ifdef CGAL_CDT_2_DEBUG_INTERSECTIONS
  using Triangulation::display_vertex;
#endif // CGAL_CDT_2_DEBUG_INTERSECTIONS


  typedef typename Triangulation::Edge             Edge;
  typedef typename Triangulation::Vertex           Vertex;
  typedef typename Triangulation::Vertex_handle    Vertex_handle;
  typedef typename Triangulation::Face_handle      Face_handle;
  typedef typename Triangulation::Face_circulator  Face_circulator ;
  typedef typename Triangulation::Vertex_iterator  Vertex_iterator;
  typedef typename Triangulation::Vertex_circulator  Vertex_circulator;
  typedef typename Triangulation::Locate_type      Locate_type;
  typedef typename Triangulation::Line_face_circulator Line_face_circulator;
  typedef typename Triangulation::Geom_traits      Geom_traits;
  typedef typename Geom_traits::Point_2            Point;
  typedef typename Geom_traits::Segment_2          Segment;
  typedef typename Triangulation::Constraint       Constraint;
  typedef typename Triangulation::size_type        size_type;

  typedef typename Triangulation::List_edges       List_edges;
  typedef typename Triangulation::List_faces       List_faces;
  typedef typename Triangulation::List_vertices    List_vertices;
  typedef typename Triangulation::List_constraints List_constraints;
  typedef typename Triangulation::Constrained_edges_iterator Constrained_edges_iterator;

  typedef Pct2_vertex_handle_less_xy<Self>         Vh_less_xy;
  typedef Polyline_constraint_hierarchy_2<Vertex_handle, Vh_less_xy, Point>
                                                   Constraint_hierarchy;
public:
  // Tag to mark the presence of a hierarchy of constraints
  typedef Tag_true                                 Constraint_hierarchy_tag;

  //Tag to distinguish Delaunay from regular triangulations
  typedef Tag_false                                Weighted_tag;

  // Tag to distinguish periodic triangulations from others
  typedef Tag_false                                Periodic_tag;

  // for user interface with the constraint hierarchy
  typedef typename Constraint_hierarchy::Vertex_it
                                            Vertices_in_constraint_iterator;

  typedef Iterator_range<Vertices_in_constraint_iterator> Vertices_in_constraint;

  typedef typename Constraint_hierarchy::Point_it
                                            Points_in_constraint_iterator;
  typedef Iterator_range<Points_in_constraint_iterator> Points_in_constraint;

  typedef typename Constraint_hierarchy::Context          Context;
  typedef typename Constraint_hierarchy::Context_iterator Context_iterator;
  typedef Iterator_range<Context_iterator>                Contexts;

  typedef typename Constraint_hierarchy::C_iterator   Constraint_iterator;
  typedef Iterator_range<Constraint_iterator> Constraints;

  typedef typename Constraint_hierarchy::Subconstraint_iterator  Subconstraint_iterator;
  typedef Iterator_range<Subconstraint_iterator> Subconstraints;

  typedef typename Constraint_hierarchy::Constraint_id Constraint_id;

  typedef std::pair<Vertex_handle, Vertex_handle> Subconstraint;

  using Triangulation::geom_traits;
  using Triangulation::cw;
  using Triangulation::ccw;
  using Triangulation::incident_faces;

protected:
  Constraint_hierarchy hierarchy;

public:
  Constraint_hierarchy& hierarchy_ref()
  {
    return hierarchy;
  }

  Constrained_triangulation_plus_2(const Geom_traits& gt=Geom_traits())
    : Triangulation(gt)
    , hierarchy(Vh_less_xy(this))
  { }

  Constrained_triangulation_plus_2(const Constrained_triangulation_plus_2& ctp)
    : Triangulation()
    , hierarchy(Vh_less_xy(this))
  { copy_triangulation(ctp);}

  Constrained_triangulation_plus_2(Constrained_triangulation_plus_2&&) = default;

  virtual ~Constrained_triangulation_plus_2() {}

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
     : Triangulation(gt)
     , hierarchy(Vh_less_xy(this))
  {
    insert_constraints(first, last);
    CGAL_postcondition( this->is_valid() );
  }


  Constrained_triangulation_plus_2(const std::list<std::pair<Point,Point> > &constraints,
                                   const Geom_traits& gt=Geom_traits() )
    : Triangulation(gt)
     , hierarchy(Vh_less_xy(this))
  {
    insert_constraints(constraints.begin(), constraints.end());
    CGAL_postcondition( this->is_valid() );
  }
  //Helping
  void clear() { Base::clear(); hierarchy.clear(); }
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
    Constraint_id cid = hierarchy.insert_constraint_old_API(va, vb);
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


  template <class OutputIterator>
  Vertices_in_constraint_iterator
  remove_vertex_from_constraint(Constraint_id cid, Vertices_in_constraint_iterator pos,
                                OutputIterator out)
  {
    if(pos == vertices_in_constraint_begin(cid)){
      ++pos;
      Constraint_id aux = hierarchy.split2(cid,pos);
      remove_constraint(aux, out);
      return pos;
    }

    Vertices_in_constraint_iterator it = vertices_in_constraint_end(cid);
    it--;
    if(pos == it){
      --pos;
      Constraint_id aux = hierarchy.split(cid, pos);
      remove_constraint(aux, out);
      return vertices_in_constraint_end(cid);
    }

    Vertices_in_constraint_iterator pp = pos;
    --pp;
    Vertex_handle a = *pp;
    pp = pos;
    ++pp;
    Vertex_handle b = *pp;
    --it;
    Vertices_in_constraint_iterator beg = vertices_in_constraint_begin(cid);
    ++beg;
    Face_container<Constrained_triangulation_plus_2> fc(*this);

    Constraint_id head = nullptr, tail = nullptr;
    if(pos != beg){
      // split off head
      --pos;
      head = hierarchy.split2(cid, pos);
      ++pos;
    }
    if(pos != it){
      // split off tail
      ++pos;
      tail = hierarchy.split(cid,pos);
    }

    Constraint_id aux = insert_constraint(a, b, std::back_inserter(fc));
    pos = vertices_in_constraint_end(aux);
    --pos;
    --pos; // this is not necessarily == vertices_in_constraint_begin(aux);
    hierarchy.swap(cid, aux);
    remove_constraint(aux, std::back_inserter(fc));

    if(head.vl_ptr()){
      hierarchy.concatenate2(head, cid);
    }

    if(tail.vl_ptr()){
      hierarchy.concatenate(cid, tail);
    }
    fc.write_faces(out);
    ++pos; // we went one too far back because the last vertex gets removed by concatenate
    return pos;
  }

  // Inserts vh before pos
  // Returns an iterator pointing on the newly inserted vertex
  // Writes the modified faces to out
  template <class OutputIterator>
  Vertices_in_constraint_iterator
  insert_vertex_in_constraint(Constraint_id cid, Vertices_in_constraint_iterator pos,
                              Vertex_handle vh, OutputIterator out)
  {
    // Insertion before the first vertex
    if(pos == vertices_in_constraint_begin(cid)){
      //std::cout << "insertion before first vertex" << std::endl;
      Constraint_id head = insert_constraint(vh, *pos, out);
      hierarchy.concatenate2(head, cid);
      return vertices_in_constraint_begin(cid);
    }

    // Insertion after the last vertex
    if(pos == vertices_in_constraint_end(cid)){
      //std::cout << "insertion after last vertex" << std::endl;
      pos--;
      Constraint_id tail = insert_constraint(*pos, vh, out);
      pos = vertices_in_constraint_end(tail);
      --pos;
      hierarchy.concatenate(cid, tail);
      return pos;
    }
    Vertex_handle b = *pos;
    --pos;
    Vertex_handle a = *pos;
    ++pos;
    Face_container<Constrained_triangulation_plus_2> fc(*this);
    Vertices_in_constraint_iterator beg = vertices_in_constraint_begin(cid), vcit;
    ++beg;
    vcit = beg;
    ++beg;
    // If the constraint consists only of a segment, and we want to insert
    // in the middle
    if((pos == vcit) && (beg == vertices_in_constraint_end(cid))){
      //std::cout << "insertion in constraint which is a segment" << std::endl;
      Constraint_id aux1 = insert_constraint(a, vh, std::back_inserter(fc));
      Constraint_id aux2 = insert_constraint(vh, b, std::back_inserter(fc));
      pos = vertices_in_constraint_begin(aux2);
      concatenate(aux1, aux2);
      hierarchy.swap(cid, aux1);
      remove_constraint(aux1, std::back_inserter(fc));
      fc.write_faces(out);
      return pos;

    }
    Constraint_id head = nullptr, tail = nullptr;
    Vertices_in_constraint_iterator bit = vertices_in_constraint_begin(cid);
    Vertices_in_constraint_iterator pred = pos;
    --pred;
    ++bit;
    if(pos != bit){
      //std::cout << "split head" << std::endl;
      head = split(cid, pred);
      std::swap(head,cid); // split2 does the job
      pred = vertices_in_constraint_begin(cid);
      pos = pred;
      ++pos;
    }
    Vertices_in_constraint_iterator eit = vertices_in_constraint_end(cid);
    --eit;
    if(pos != eit){
      //std::cout << "split tail" << std::endl;
      tail = split(cid, pos);
    }

    // make the new constraint
    Constraint_id aux1 = insert_constraint(a, vh, std::back_inserter(fc));
    Constraint_id aux2 = insert_constraint(vh, b, std::back_inserter(fc));
    pos = vertices_in_constraint_begin(aux2);
    concatenate(aux1, aux2);

    if(head.vl_ptr()){
      //std::cout << "concatenate head" << std::endl;
      remove_constraint(cid, std::back_inserter(fc));
      hierarchy.concatenate(head, aux1);
    } else {
      hierarchy.swap(cid, aux1);
      remove_constraint(aux1, std::back_inserter(fc));
      head = cid;
    }

    if(tail.vl_ptr()){
      //std::cout << "concatenate tail" << std::endl;
      concatenate(head, tail);
    }
    fc.write_faces(out);
    return pos;
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
    Constraint_id ca = hierarchy.insert_constraint(vertices[0],vertices[1]);
    insert_subconstraint(vertices[0],vertices[1], std::back_inserter(fc));

    if(n>2){
      for(int j=1; j<n-1; j++){
        hierarchy.append_constraint(ca, vertices[j], vertices[j+1]);
        insert_subconstraint(vertices[j], vertices[j+1], std::back_inserter(fc));
      }
    }
    for(Vertices_in_constraint_iterator vcit = vertices_in_constraint_begin(ca);
        vcit != vertices_in_constraint_end(ca);
        vcit++){
      insert_incident_faces(vcit, out);
    }
    //AF    vertices_in_constraint_begin(ca)->fixed() = true;
    // Vertices_in_constraint_iterator end = boost::prior(vertices_in_constraint_end(ca));
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
      return nullptr;
    }
    CGAL_assertion(n >= 2);

    Constraint_id ca = hierarchy.insert_constraint(vertices[0],vertices[1]);
    insert_subconstraint(vertices[0],vertices[1]);

    if(n>2){
      for(std::size_t j=1; j<n-1; j++){
        hierarchy.append_constraint(ca, vertices[j], vertices[j+1]);
        insert_subconstraint(vertices[j], vertices[j+1]);
      }
    }

    // fix first and last, one is redundant for is_polygon == true
    // vertices.front()->fixed() = true;
    // vertices.back()->fixed() = true;

    return ca;
  }

public:

  void
  file_output(std::ostream& os) const
  {
    os << static_cast<const Tr&>(*this);
    Unique_hash_map<Vertex_handle,int> V(0, number_of_vertices());
    int inum = 0;
    for(Vertex_iterator vit = vertices_begin(); vit != vertices_end() ; ++vit){
      if(! is_infinite(vit)){
        V[vit] = inum++;
      }
    }

    for(Constraint_iterator cit = constraints_begin(); cit != constraints_end(); ++cit){
      os << (*cit).second->all_size();
      for(Vertex_handle vh : vertices_in_constraint(*cit)){
         os << " " << V[vh];
       }
       os << std::endl;
    }
  }


  void file_input(std::istream& is)
  {

    is >> static_cast<Tr&>(*this);

    std::vector<Vertex_handle> V;
    V.reserve(number_of_vertices());
    for(Vertex_iterator vit = vertices_begin(); vit != vertices_end() ; ++vit){
      if(! is_infinite(vit)){
        V.push_back(vit);
      }
    }
    Constraint_id cid;
    int n, i0, i1;
    while(is >> n){
      is >> i0 >> i1;
      cid = insert_constraint(V[i0],V[i1]);

      for(int i = 2; i < n; i++){
        i0 = i1;
        is >> i1;
        Constraint_id cid2 = insert_constraint(V[i0],V[i1]);
        cid = concatenate(cid, cid2);
      }
    }
  }


  template <class OutputIterator>
  typename Constrained_triangulation_plus_2<Tr>::Constraint_id
  insert_constraint(Vertex_handle va, Vertex_handle vb, OutputIterator out)
  {
    // protects against inserting a zero length constraint
    if(va == vb){
    return Constraint_id(nullptr);
    }
    // protects against inserting twice the same constraint
    Constraint_id cid = hierarchy.insert_constraint(va, vb);
    if (va != vb && (cid != nullptr) )  insert_subconstraint(va,vb,out);

    for(Vertices_in_constraint_iterator vcit = vertices_in_constraint_begin(cid);
        vcit != vertices_in_constraint_end(cid);
        vcit++){
      insert_incident_faces(vcit, out);
    }
    return cid;
  }

  virtual Vertex_handle intersect(Face_handle f, int i,
                                  Vertex_handle vaa,
                                  Vertex_handle vbb);
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
    std::list<Vertex_handle> vertices(hierarchy.vertices_in_constraint_begin(cid),
                                      hierarchy.vertices_in_constraint_end(cid));

    hierarchy.remove_constraint(cid);
    for(typename std::list<Vertex_handle>::iterator it = vertices.begin(),
          succ = it;
        ++succ != vertices.end();
        ++it){
      if(! is_subconstraint(*it, *succ)){ // this checks whether other constraints pass
        Face_handle fh;
        int i;
        bool b = Triangulation::is_edge(*it, *succ, fh, i);
        CGAL_assume(b);
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
    Vertices_in_constraint_iterator u = boost::prior(v);
    Vertices_in_constraint_iterator w = boost::next(v);
    bool unew = (*u != *w);
    hierarchy.simplify(u,v,w);

    Triangulation::remove_incident_constraints(*v);

    Triangulation::remove(*v);

    if(unew){
      Triangulation::insert_constraint(*u, *w);
    }
  }

  std::size_t remove_points_without_corresponding_vertex(Constraint_id cid)
  {
    return hierarchy.remove_points_without_corresponding_vertex(cid);
  }
  std::size_t remove_points_without_corresponding_vertex()
  {
    return hierarchy.remove_points_without_corresponding_vertex();
  }


  // CONCATENATE AND SPLIT

  // concatenates two constraints
  Constraint_id
  concatenate(Constraint_id first, Constraint_id second);

  // split a constraint in two constraints, so that vcit becomes the first
  // vertex of the new constraint
  // returns the new constraint
  Constraint_id
  split(Constraint_id first, Vertices_in_constraint_iterator vcit);

  // Query of the constraint hierarchy
  Constraint_iterator constraints_begin() const;
  Constraint_iterator constraints_end()   const;
  Constraints constraints() const
  {
    return Constraints(constraints_begin(),constraints_end());
  }

  Subconstraint_iterator subconstraints_begin() const;
  Subconstraint_iterator subconstraints_end() const;

  Subconstraints subconstraints() const
  {
    return Subconstraints(subconstraints_begin(),subconstraints_end());
  }

  Context   context(Vertex_handle va, Vertex_handle vb); //AF: const;

  bool is_subconstraint(Vertex_handle va,
                        Vertex_handle vb);
  size_type number_of_enclosing_constraints(Vertex_handle va,
                                            Vertex_handle vb) const;
  Context_iterator   contexts_begin(Vertex_handle va,
                                    Vertex_handle vb) const;
  Context_iterator   contexts_end(Vertex_handle va,
                                  Vertex_handle vb) const;

  Contexts contexts(Vertex_handle va, Vertex_handle vb) const
  {
    return Contexts(contexts_begin(va,vb),contexts_end(va,vb));
  }

  Vertices_in_constraint_iterator vertices_in_constraint_begin(Constraint_id cid) const;
  Vertices_in_constraint_iterator vertices_in_constraint_end(Constraint_id cid) const;

  Vertices_in_constraint vertices_in_constraint(Constraint_id cid) const
  {
    return Vertices_in_constraint(vertices_in_constraint_begin(cid), vertices_in_constraint_end(cid));
  }

  Points_in_constraint_iterator points_in_constraint_begin(Constraint_id cid) const;
  Points_in_constraint_iterator points_in_constraint_end(Constraint_id cid) const ;

  Points_in_constraint points_in_constraint(Constraint_id cid) const
  {
    return Points_in_constraint(points_in_constraint_begin(cid), points_in_constraint_end(cid));
  }

  size_type number_of_constraints() {
    return static_cast<size_type> (hierarchy.number_of_constraints());}
  size_type number_of_subconstraints(){
    return static_cast<size_type> (hierarchy.number_of_subconstraints());}

  // public member, used by Mesh_2::Refine_edges
  void split_constraint(Vertex_handle v1, Vertex_handle v2,
                        Vertex_handle va) {
    hierarchy.split_constraint(v1,v2,va);
  }

protected:
  template <class OutputItertator>
  void insert_incident_faces(Vertices_in_constraint_iterator vcit, OutputItertator out)
  {
    Vertex_handle vh = *vcit;
    Face_circulator fc = incident_faces(vh), done = fc;
    Face_circulator null ;
    if ( fc != null )
    {
      do {
        Face_handle fh = fc;
        out = fh;
        out++;
        fc++;
      }while(fc != done);
    }
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
    boost::tie(vaa,vbb) = stack.top();
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
        hierarchy.split_constraint(vaa,vbb,vi);
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
        hierarchy.split_constraint(vaa,vbb,vi);
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
      hierarchy.split_constraint(vaa,vbb,vi);
      stack.push(std::make_pair(vi,vbb));
    }
  }
}



  //to debug
public:
  void print_hierarchy() { hierarchy.print(); }

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
  hierarchy.copy(ctp.hierarchy, vmap);
}

template <class Tr>
void
Constrained_triangulation_plus_2<Tr>::
swap(Constrained_triangulation_plus_2 &ctp)
{
  Base::swap(ctp);
  hierarchy.swap(ctp.hierarchy);
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
    hierarchy.split_constraint(v1,v2,va);
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
  Vertex_handle  vc, vd, va, vb;
  Vertex_handle  vcc, vdd;
  vcc = f->vertex(cw(i));
  vdd = f->vertex(ccw(i));
  CGAL_assertion_code( bool b1 = )
  hierarchy.enclosing_constraint(vcc,vdd,vc,vd);
  CGAL_assertion_code( bool b2 = )
  hierarchy.enclosing_constraint(vaa,vbb,va,vb);
  CGAL_assertion(b1);
  CGAL_assertion(b2);

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
    hierarchy.split_constraint(vcc,vdd,vi);
    insert_subconstraint(vcc,vi);
    insert_subconstraint(vi, vdd);
  }
  else {
    insert_subconstraint(vcc,vdd);
  }
  return vi;
}

  // CONCATENATE AND SPLIT

  // concatenates two constraints
template <class Tr>
typename Constrained_triangulation_plus_2<Tr>::Constraint_id
Constrained_triangulation_plus_2<Tr>::concatenate(Constraint_id first, Constraint_id second)
{
  return hierarchy.concatenate(first,second);
}

  // split a constraint in two constraints, so that vcit becomes the first
  // vertex of the new constraint
  // returns the new constraint
template <class Tr>
typename Constrained_triangulation_plus_2<Tr>::Constraint_id
Constrained_triangulation_plus_2<Tr>::split(Constraint_id first, Vertices_in_constraint_iterator vcit)
{
  return hierarchy.split(first, vcit);
}


template <class Tr>
std::ostream &
operator<<(std::ostream& os,
           const Constrained_triangulation_plus_2<Tr> &ct)
{
  ct.file_output(os);
  return os ;
}

template <class Tr>
std::istream &
operator>>(std::istream& is,
           Constrained_triangulation_plus_2<Tr> &ct)
{
  ct.file_input(is);
  return is ;
}

// Constraint Hierarchy Queries

template <class Tr>
inline
typename
Constrained_triangulation_plus_2<Tr>::Constraint_iterator
Constrained_triangulation_plus_2<Tr>::
constraints_begin() const
{
  return hierarchy.c_begin();
}

template <class Tr>
inline
typename
Constrained_triangulation_plus_2<Tr>::Constraint_iterator
Constrained_triangulation_plus_2<Tr>::
constraints_end() const
{
  return hierarchy.c_end();
}

template <class Tr>
inline
typename
Constrained_triangulation_plus_2<Tr>::Subconstraint_iterator
Constrained_triangulation_plus_2<Tr>::
subconstraints_begin() const
{
  return hierarchy.subconstraint_begin();
}

template <class Tr>
inline
typename
Constrained_triangulation_plus_2<Tr>::Subconstraint_iterator
Constrained_triangulation_plus_2<Tr>::
subconstraints_end() const
{
  return hierarchy.subconstraint_end();
}


template <class Tr>
inline
typename Constrained_triangulation_plus_2<Tr>::Context
Constrained_triangulation_plus_2<Tr>::
context(Vertex_handle va, Vertex_handle vb) // AF: const
{
  return hierarchy.context(va,vb);
}


template <class Tr>
inline
typename Constrained_triangulation_plus_2<Tr>::size_type
Constrained_triangulation_plus_2<Tr>::
number_of_enclosing_constraints(Vertex_handle va, Vertex_handle vb) const
{
 return static_cast<size_type>
   (hierarchy.number_of_enclosing_constraints(va,vb));
}

template <class Tr>
inline bool
Constrained_triangulation_plus_2<Tr>::
is_subconstraint(Vertex_handle va, Vertex_handle vb)
{
 return hierarchy.is_subconstrained_edge(va,vb);
}


template <class Tr>
inline
typename Constrained_triangulation_plus_2<Tr>::Context_iterator
Constrained_triangulation_plus_2<Tr>::
contexts_begin(Vertex_handle va, Vertex_handle vb) const
{
  return hierarchy.contexts_begin(va,vb);
}

template <class Tr>
inline
typename Constrained_triangulation_plus_2<Tr>::Context_iterator
Constrained_triangulation_plus_2<Tr>::
contexts_end(Vertex_handle va, Vertex_handle vb) const
{
  return hierarchy.contexts_end(va,vb);
}

template <class Tr>
inline
typename Constrained_triangulation_plus_2<Tr>::Vertices_in_constraint_iterator
Constrained_triangulation_plus_2<Tr>::
vertices_in_constraint_begin(Constraint_id cid) const
{
  return  hierarchy.vertices_in_constraint_begin(cid);
}

template <class Tr>
inline
typename Constrained_triangulation_plus_2<Tr>::Vertices_in_constraint_iterator
Constrained_triangulation_plus_2<Tr>::
vertices_in_constraint_end(Constraint_id cid) const
{
  return  hierarchy.vertices_in_constraint_end(cid);
}

template <class Tr>
inline
typename Constrained_triangulation_plus_2<Tr>::Points_in_constraint_iterator
Constrained_triangulation_plus_2<Tr>::
points_in_constraint_begin(Constraint_id cid) const
{
  return  hierarchy.points_in_constraint_begin(cid);
}

template <class Tr>
inline
typename Constrained_triangulation_plus_2<Tr>::Points_in_constraint_iterator
Constrained_triangulation_plus_2<Tr>::
points_in_constraint_end(Constraint_id cid) const
{
  return  hierarchy.points_in_constraint_end(cid);
}

} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif //CGAL_CONSTRAINED_TRIANGULATION_PLUS_2_H
