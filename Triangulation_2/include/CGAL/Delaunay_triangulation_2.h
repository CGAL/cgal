// Copyright (c) 1997  INRIA Sophia-Antipolis (France).
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
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Mariette Yvinec
//               : Olivier Devillers  (remove)
//               : Pedro de Castro (displacement)

#ifndef CGAL_DELAUNAY_TRIANGULATION_2_H
#define CGAL_DELAUNAY_TRIANGULATION_2_H

#include <CGAL/license/Triangulation_2.h>

#include <CGAL/Triangulation_2.h>
#include <CGAL/iterator.h>
#include <CGAL/Object.h>

#ifndef CGAL_TRIANGULATION_2_DONT_INSERT_RANGE_OF_POINTS_WITH_INFO
#include <CGAL/Spatial_sort_traits_adapter_2.h>
#include <CGAL/internal/info_check.h>
#include <CGAL/tss.h>
#include <boost/iterator/zip_iterator.hpp>
#include <boost/mpl/and.hpp>

#endif //CGAL_TRIANGULATION_2_DONT_INSERT_RANGE_OF_POINTS_WITH_INFO

namespace CGAL {

template < class Gt,
           class Tds = Triangulation_data_structure_2 <
                         Triangulation_vertex_base_2<Gt> > >
class Delaunay_triangulation_2
    : public Triangulation_2<Gt,Tds>
{
public:
  typedef Gt Geom_traits;
  typedef typename Geom_traits::Point_2       Point;
  typedef typename Geom_traits::Segment_2     Segment;
  typedef typename Geom_traits::Triangle_2    Triangle;

  typedef typename Geom_traits::Orientation_2 Orientation_2;
  typedef typename Geom_traits::Compare_x_2   Compare_x;
  typedef typename Geom_traits::Compare_y_2   Compare_y;
  typedef typename Geom_traits::Side_of_oriented_circle_2
                                              Side_of_oriented_circle;


  typedef Triangulation_2<Gt,Tds>                       Triangulation;
  typedef typename Triangulation::size_type             size_type;
  typedef typename Triangulation::Locate_type           Locate_type;
  typedef typename Triangulation::Face_handle           Face_handle;
  typedef typename Triangulation::Vertex_handle         Vertex_handle;
  typedef typename Triangulation::Edge                  Edge;
  typedef typename Triangulation::Edge_circulator       Edge_circulator;
  typedef typename Triangulation::Face_circulator       Face_circulator;
  typedef typename Triangulation::Vertex_circulator     Vertex_circulator;
  typedef typename Triangulation::Finite_edges_iterator Finite_edges_iterator;
  typedef typename Triangulation::Finite_faces_iterator Finite_faces_iterator;
  typedef typename Triangulation::Finite_vertices_iterator
                                                        Finite_vertices_iterator;
  typedef typename Triangulation::All_faces_iterator    All_faces_iterator;

  //Tag to distinguish Delaunay from regular triangulations
  typedef Tag_false                                     Weighted_tag;

  // Tag to distinguish periodic triangulations from others
  typedef Tag_false                                     Periodic_tag;

#ifndef CGAL_CFG_USING_BASE_MEMBER_BUG_2
  using Triangulation::side_of_oriented_circle;
  using Triangulation::circumcenter;
  using Triangulation::collinear_between;
  using Triangulation::test_dim_down;
  using Triangulation::make_hole;
  using Triangulation::fill_hole_delaunay;
  using Triangulation::delete_vertex;
#endif

 Delaunay_triangulation_2(const Gt& gt = Gt())
  : Triangulation_2<Gt,Tds>(gt) {}

  Delaunay_triangulation_2(
         const Delaunay_triangulation_2<Gt,Tds> &tr)
       : Triangulation_2<Gt,Tds>(tr)
  {   CGAL_triangulation_postcondition(is_valid());  }

 template <class InputIterator>
 Delaunay_triangulation_2(InputIterator first, InputIterator last,
                          const Gt& gt = Gt())
   : Triangulation_2<Gt,Tds>(gt)
 {
  insert(first,last);
 }

// CHECK -QUERY
  bool is_valid(bool verbose = false, int level = 0) const;

  Vertex_handle
  nearest_vertex(const Point& p, Face_handle f= Face_handle()) const;

  bool does_conflict(const Point &p, Face_handle fh) const;// deprecated
  bool test_conflict(const Point &p, Face_handle fh, bool strict = true) const;
  bool find_conflicts(const Point &p,               //deprecated
                      std::list<Face_handle>& conflicts,
                      Face_handle start= Face_handle()) const;

//  //template member functions, declared and defined at the end
//  template <class OutputItFaces, class OutputItBoundaryEdges>
//  std::pair<OutputItFaces,OutputItBoundaryEdges>
//  get_conflicts_and_boundary(const Point  &p,
//                             OutputItFaces fit,
//                             OutputItBoundaryEdges eit,
//                             Face_handle start) const;
//  template <class OutputItFaces>
//  OutputItFaces
//  get_conflicts (const Point  &p,
//                 OutputItFaces fit,
//                 Face_handle start) const;
//  template <class OutputItBoundaryEdges>
//  OutputItBoundaryEdges
//  get_boundary_of_conflicts(const Point  &p,
//                            OutputItBoundaryEdges eit,
//                            Face_handle start) const;

  // DUAL
  Point dual (Face_handle f) const;
  Object dual(const Edge &e) const;
  Object dual(const Edge_circulator& ec) const;
  Object dual(const Finite_edges_iterator& ei) const;

  //INSERTION-REMOVAL
  Vertex_handle insert(const Point  &p,
                       Face_handle start = Face_handle());
  Vertex_handle insert(const Point& p,
                       Locate_type lt,
                       Face_handle loc, int li);
  Vertex_handle push_back(const Point &p);

  void remove(Vertex_handle v);

  // DISPLACEMENT
  void restore_Delaunay(Vertex_handle v);

  Vertex_handle move_if_no_collision(Vertex_handle v, const Point &p);
  Vertex_handle move(Vertex_handle v, const Point &p);

protected: // some internal stuffs
  template <class OutputItFaces>
  Vertex_handle insert_and_give_new_faces(const Point  &p,
                                          OutputItFaces fit,
                                          Face_handle start = Face_handle());
  template <class OutputItFaces>
  Vertex_handle insert_and_give_new_faces(const Point& p,
                                          Locate_type lt,
                                          Face_handle loc, int li,
                                          OutputItFaces fit);

  template <class OutputItFaces>
  Vertex_handle move_if_no_collision_and_give_new_faces(Vertex_handle v,
                                                        const Point &p,
                                                        OutputItFaces fit);

  template <class OutputItFaces>
  void remove_and_give_new_faces(Vertex_handle v,
                                 OutputItFaces fit);

public:
  bool is_delaunay_after_displacement(Vertex_handle v,
                                      const Point &p) const;

#ifndef CGAL_CFG_USING_BASE_MEMBER_BUG_2
  using Triangulation::cw;
  using Triangulation::ccw;
  using Triangulation::geom_traits;
#endif

private:
#ifndef CGAL_DT2_USE_RECURSIVE_PROPAGATING_FLIP
  void non_recursive_propagating_flip(Face_handle f,int i);
  void propagating_flip(const Face_handle& f, int i, int depth=0);
#else
  void propagating_flip(const Face_handle& f,int i);
#endif

// auxilliary functions for remove
  void remove_degree_init(Vertex_handle v, std::vector<Face_handle> &f,
         std::vector<Vertex_handle> &w, std::vector<int> &i,int&d,int&maxd);
  void remove_degree_triangulate(Vertex_handle v, std::vector<Face_handle> &f,
         std::vector<Vertex_handle> &w, std::vector<int> &i,int d);
  void remove_degree_d(Vertex_handle v, std::vector<Face_handle> &f,
         std::vector<Vertex_handle> &w, std::vector<int> &i,int d);
  void remove_degree3(Vertex_handle v, std::vector<Face_handle> &f,
         std::vector<Vertex_handle> &w, std::vector<int> &i);
  void remove_degree4(Vertex_handle v, std::vector<Face_handle> &f,
         std::vector<Vertex_handle> &w, std::vector<int> &i);
  void remove_degree5(Vertex_handle v, std::vector<Face_handle> &f,
         std::vector<Vertex_handle> &w, std::vector<int> &i);
  void remove_degree5_star(Vertex_handle &v,
                           Face_handle &,Face_handle &,Face_handle &,
                           Face_handle &,Face_handle &,
                           Vertex_handle&,Vertex_handle&,Vertex_handle&,
                           Vertex_handle&,Vertex_handle&,
                           int,int,int,
                           int,int);
  void remove_degree6(Vertex_handle v, std::vector<Face_handle> &f,
                      std::vector<Vertex_handle> &w, std::vector<int> &i);
  void remove_degree6_star   (Vertex_handle &v,
                              Face_handle &,Face_handle &,Face_handle &,
                              Face_handle &,Face_handle &,Face_handle &,
                              Vertex_handle&,Vertex_handle&,Vertex_handle&,
                              Vertex_handle&,Vertex_handle&,Vertex_handle&,
                              int,int,int,
                              int,int,int);
  void remove_degree6_N      (Vertex_handle &v,
                              Face_handle &,Face_handle &,Face_handle &,
                              Face_handle &,Face_handle &,Face_handle &,
                              Vertex_handle&,Vertex_handle&,Vertex_handle&,
                              Vertex_handle&,Vertex_handle&,Vertex_handle&,
                              int,int,int,
                              int,int,int);
  void remove_degree6_antiN  (Vertex_handle &v,
                              Face_handle &,Face_handle &,Face_handle &,
                              Face_handle &,Face_handle &,Face_handle &,
                              Vertex_handle&,Vertex_handle&,Vertex_handle&,
                              Vertex_handle&,Vertex_handle&,Vertex_handle&,
                              int,int,int,
                              int,int,int);
  void remove_degree6_diamond(Vertex_handle &v,
                              Face_handle &,Face_handle &,Face_handle &,
                              Face_handle &,Face_handle &,Face_handle &,
                              Vertex_handle&,Vertex_handle&,Vertex_handle&,
                              Vertex_handle&,Vertex_handle&,Vertex_handle&,
                              int,int,int,
                              int,int,int);
  void remove_degree7(Vertex_handle v,std::vector<Face_handle> &f,
                      std::vector<Vertex_handle> &w, std::vector<int> &i);
  bool incircle(int x, int j, int, int l, std::vector<Face_handle> &f,
                std::vector<Vertex_handle> &w, std::vector<int> &i)
  {
    // k is supposed to be j+1 modulo degree, x is supposed to be finite
    //test if w[x] inside circle w[j]w[k]w[l] (f[j] has vertices w[j]w[k])
    // THE FOLLOWING LINE IS TO BE REMOVED. JUST THERE FOR STUPID PRECONDITION
    //if(geom_traits().orientation_2_object()(w[j]->point(),w[k]->point(),w[l]->point())!=POSITIVE) return true;
    f[j]->set_vertex(i[j], w[l]); // change vertex v for another one
    return (test_conflict(w[x]->point(), f[j]));
  }
  void rotate7(int j, std::vector<Vertex_handle> &w,
         std::vector<Face_handle> &f, std::vector<int> &i);
  void remove_degree7_star      (Vertex_handle&,int,std::vector<Face_handle> &f,
         std::vector<Vertex_handle> &w, std::vector<int> &i);
  void remove_degree7_zigzag    (Vertex_handle&,int,std::vector<Face_handle> &f,
         std::vector<Vertex_handle> &w, std::vector<int> &i);
  void remove_degree7_leftdelta (Vertex_handle&,int,std::vector<Face_handle> &f,
         std::vector<Vertex_handle> &w, std::vector<int> &i);
  void remove_degree7_rightdelta(Vertex_handle&,int,std::vector<Face_handle> &f,
         std::vector<Vertex_handle> &w, std::vector<int> &i);
  void remove_degree7_leftfan   (Vertex_handle&,int,std::vector<Face_handle> &f,
         std::vector<Vertex_handle> &w, std::vector<int> &i);
  void remove_degree7_rightfan  (Vertex_handle&,int,std::vector<Face_handle> &f,
         std::vector<Vertex_handle> &w, std::vector<int> &i);
// end of auxilliary functions for remove

  Vertex_handle nearest_vertex_2D(const Point& p, Face_handle f) const;
  Vertex_handle nearest_vertex_1D(const Point& p) const;

  void  look_nearest_neighbor(const Point& p,
                              Face_handle f,
                              int i,
                              Vertex_handle& nn) const;

public:
  template < class Stream>
  Stream& draw_dual(Stream & ps)
  {
    Finite_edges_iterator eit= this->finite_edges_begin();
    for (; eit != this->finite_edges_end(); ++eit) {
      Object o = dual(eit);
      typename Geom_traits::Line_2  l;
      typename Geom_traits::Ray_2   r;
      Segment s;
      if(CGAL::assign(s,o)) ps << s;
      if(CGAL::assign(r,o)) ps << r;
      if(CGAL::assign(l,o)) ps << l;
    }
    return ps;
  }

#ifndef CGAL_TRIANGULATION_2_DONT_INSERT_RANGE_OF_POINTS_WITH_INFO
  template < class InputIterator >
  std::ptrdiff_t
  insert(InputIterator first, InputIterator last,
         typename boost::enable_if<
         boost::is_convertible<
         typename std::iterator_traits<InputIterator>::value_type,
         Point
         >
         >::type* = NULL)
#else
  template < class InputIterator >
  std::ptrdiff_t
  insert(InputIterator first, InputIterator last)
#endif //CGAL_TRIANGULATION_2_DONT_INSERT_RANGE_OF_POINTS_WITH_INFO
  {
    size_type n = this->number_of_vertices();

    std::vector<Point> points (first, last);
    spatial_sort (points.begin(), points.end(), geom_traits());
    Face_handle f;
    for (typename std::vector<Point>::const_iterator p = points.begin(), end = points.end();
         p != end; ++p)
      f = insert (*p, f)->face();

    return this->number_of_vertices() - n;
  }

#ifndef CGAL_TRIANGULATION_2_DONT_INSERT_RANGE_OF_POINTS_WITH_INFO
private:
  //top stands for tuple-or-pair
  template <class Info>
  const Point& top_get_first(const std::pair<Point,Info>& pair) const { return pair.first; }
  template <class Info>
  const Info& top_get_second(const std::pair<Point,Info>& pair) const { return pair.second; }
  template <class Info>
  const Point& top_get_first(const boost::tuple<Point,Info>& tuple) const { return boost::get<0>(tuple); }
  template <class Info>
  const Info& top_get_second(const boost::tuple<Point,Info>& tuple) const { return boost::get<1>(tuple); }

  template <class Tuple_or_pair,class InputIterator>
  std::ptrdiff_t insert_with_info(InputIterator first,InputIterator last)
  {
    size_type n = this->number_of_vertices();
    std::vector<std::size_t> indices;
    std::vector<Point> points;
    std::vector<typename Tds::Vertex::Info> infos;
    std::size_t index=0;
    for (InputIterator it=first;it!=last;++it) {
      Tuple_or_pair value=*it;
      points.push_back(top_get_first(value));
      infos.push_back (top_get_second(value));
      indices.push_back(index++);
    }

    typedef typename Pointer_property_map<Point>::type Pmap;
    typedef Spatial_sort_traits_adapter_2<Geom_traits,Pmap> Search_traits;

    spatial_sort(indices.begin(), indices.end(),
                 Search_traits(make_property_map(points),geom_traits()));

    Vertex_handle v_hint;
    Face_handle hint;
    for (typename std::vector<std::size_t>::const_iterator
      it = indices.begin(), end = indices.end();
      it != end; ++it) {
      v_hint = insert(points[*it], hint);
      if(v_hint!=Vertex_handle()) {
        v_hint->info()=infos[*it];
        hint=v_hint->face();
      }
    }

    return this->number_of_vertices() - n;
  }

public:

  template < class InputIterator >
  std::ptrdiff_t
  insert(InputIterator first,
         InputIterator last,
         typename boost::enable_if<
           boost::is_convertible<
             typename std::iterator_traits<InputIterator>::value_type,
             std::pair<Point,typename internal::Info_check<typename Tds::Vertex>::type>
           > >::type* = NULL)
  {
    return insert_with_info< std::pair<Point,typename internal::Info_check<typename Tds::Vertex>::type> >(first,last);
  }

  template <class  InputIterator_1,class InputIterator_2>
  std::ptrdiff_t
  insert(boost::zip_iterator< boost::tuple<InputIterator_1,InputIterator_2> > first,
         boost::zip_iterator< boost::tuple<InputIterator_1,InputIterator_2> > last,
         typename boost::enable_if<
           boost::mpl::and_<
             boost::is_convertible< typename std::iterator_traits<InputIterator_1>::value_type, Point >,
             boost::is_convertible< typename std::iterator_traits<InputIterator_2>::value_type, typename internal::Info_check<typename Tds::Vertex>::type >
           >
         >::type* = NULL)
  {
    return insert_with_info< boost::tuple<Point,typename internal::Info_check<typename Tds::Vertex>::type> >(first,last);
  }
#endif //CGAL_TRIANGULATION_2_DONT_INSERT_RANGE_OF_POINTS_WITH_INFO

  template <class OutputItFaces, class OutputItBoundaryEdges>
  std::pair<OutputItFaces,OutputItBoundaryEdges>
  get_conflicts_and_boundary(const Point  &p,
                             OutputItFaces fit,
                             OutputItBoundaryEdges eit,
                             Face_handle start = Face_handle(),
                             bool strict = true) const
  {
    CGAL_triangulation_precondition(this->dimension() == 2);
    int li;
    Locate_type lt;
    Face_handle fh = this->locate(p,lt,li, start);
    switch(lt) {
      case Triangulation::OUTSIDE_AFFINE_HULL:
      case Triangulation::VERTEX:
        return std::make_pair(fit,eit);
      case Triangulation::FACE:
      case Triangulation::EDGE:
      case Triangulation::OUTSIDE_CONVEX_HULL:
        *fit++ = fh; //put fh in OutputItFaces
        std::pair<OutputItFaces,OutputItBoundaryEdges> pit = std::make_pair(fit,eit);
        pit = propagate_conflicts(p,fh,0,pit, strict);
        pit = propagate_conflicts(p,fh,1,pit, strict);
        pit = propagate_conflicts(p,fh,2,pit, strict);
        return pit;
    }
    CGAL_triangulation_assertion(false);
    return std::make_pair(fit,eit);
  }

  template <class OutputItFaces>
  OutputItFaces
  get_conflicts (const Point  &p,
                 OutputItFaces fit,
                 Face_handle start= Face_handle(),
                 bool strict = true) const
  {
    std::pair<OutputItFaces,Emptyset_iterator> pp =
      get_conflicts_and_boundary(p, fit, Emptyset_iterator(), start, strict);
    return pp.first;
  }

  template <class OutputItBoundaryEdges>
  OutputItBoundaryEdges
  get_boundary_of_conflicts(const Point  &p,
                            OutputItBoundaryEdges eit,
                            Face_handle start= Face_handle(),
                            bool strict = true) const
  {
    std::pair<Emptyset_iterator, OutputItBoundaryEdges> pp =
      get_conflicts_and_boundary(p, Emptyset_iterator(), eit, start, strict);
    return pp.second;
  }

private:
#ifdef CGAL_DT2_USE_RECURSIVE_PROPAGATE_CONFLICTS
  template <class OutputItFaces, class OutputItBoundaryEdges>
  std::pair<OutputItFaces,OutputItBoundaryEdges>
  propagate_conflicts (const Point &p,
                       const Face_handle fh,
                       const int i,
                       std::pair<OutputItFaces,OutputItBoundaryEdges>
                       pit,
                       bool strict = true) const
  {
    Face_handle fn = fh->neighbor(i);
    if(! test_conflict(p,fn,strict)) {
      *(pit.second)++ = Edge(fn, fn->index(fh));
    } else {
      *(pit.first)++ = fn;
      int j = fn->index(fh);
      pit = propagate_conflicts(p,fn,ccw(j),pit,strict);
      pit = propagate_conflicts(p,fn,cw(j), pit,strict);
    }
    return pit;
  }
#else // NO CGAL_DT2_USE_RECURSIVE_PROPAGATE_CONFLICTS
  template <class OutputItFaces, class OutputItBoundaryEdges>
  std::pair<OutputItFaces,OutputItBoundaryEdges>
  non_recursive_propagate_conflicts(const Point  &p,
                                    const Face_handle fh,
                                    const int i,
                                    std::pair<OutputItFaces,OutputItBoundaryEdges> pit,
                                    bool strict = true)  const
  {
    std::stack<std::pair<Face_handle, int> > stack;
    stack.push(std::make_pair(fh,i));
    while (!stack.empty())
    {
      const Face_handle fh=stack.top().first;
      const int i=stack.top().second;
      stack.pop();
      Face_handle fn = fh->neighbor(i);
      if(! test_conflict(p,fn,strict)) {
        *(pit.second)++ = Edge(fn, fn->index(fh));
      } else {
        *(pit.first)++ = fn;
        int j = fn->index(fh);
        stack.push(std::make_pair(fn,ccw(j)));
        stack.push(std::make_pair(fn,cw(j)));
      }
    }
    return pit;
  }

  template <class OutputItFaces, class OutputItBoundaryEdges>
  std::pair<OutputItFaces,OutputItBoundaryEdges>
  propagate_conflicts (const Point  &p,
                       const Face_handle fh,
                       const int i,
                       std::pair<OutputItFaces,OutputItBoundaryEdges>
                       pit,
                       bool strict = true,
                       int depth=0) const
  {
    if(depth == 100)
      return non_recursive_propagate_conflicts(p, fh, i, pit, strict);

    Face_handle fn = fh->neighbor(i);
    if(! test_conflict(p,fn,strict)) {
      *(pit.second)++ = Edge(fn, fn->index(fh));
    } else {
      *(pit.first)++ = fn;
      int j = fn->index(fh);
      pit = propagate_conflicts(p,fn,ccw(j),pit, strict, depth+1);
      pit = propagate_conflicts(p,fn,cw(j), pit, strict, depth+1);
    }
    return pit;
  }
#endif // NO CGAL_DT2_USE_RECURSIVE_PROPAGATE_CONFLICTS

protected:
  void restore_edges(Vertex_handle v)
  {
    std::list<Edge> edges;
    Face_circulator fc = this->incident_faces(v), done(fc);
    int degree = 0;

    do {
      if((++degree) > 3) break;
    } while(++fc != done);

    fc = this->incident_faces(v);
    done = fc;

    if(degree == 3) {
      do {
        int i = fc->index(v);
        edges.push_back(Edge(fc, i));
      } while(++fc != done);
    } else {
      do {
        int i = fc->index(v);
        edges.push_back(Edge(fc, i));
        edges.push_back(Edge(fc, this->cw(i)));
      } while(++fc != done);
    }

    while(!edges.empty()) {
      const Edge &e = edges.front();
      Face_handle f = e.first;
      int i = e.second;
      edges.pop_front();
      if(this->is_infinite(f->vertex(i))) continue;
      Face_handle fi = f->neighbor(i);
      int mi = this->_tds.mirror_index(f, i);
      Vertex_handle vm = this->_tds.mirror_vertex(f, i);
      if(this->is_infinite(vm)) continue;
      if(this->side_of_oriented_circle(f, vm->point(),true) == ON_POSITIVE_SIDE) {
        this->_tds.flip(f, i);
        edges.push_back(Edge(f, i));
        edges.push_back(Edge(f, this->cw(i)));
        edges.push_back(Edge(fi, this->cw(mi)));
        edges.push_back(Edge(fi, mi));
      }
    }
  }

  void restore_edges(Vertex_handle v, std::set<Face_handle> &faces)
  {
    typedef std::list<Edge> Edges_list;
    Edges_list edges;
    Face_circulator fc = this->incident_faces(v), done(fc);
    int degree = 0;

    do {
      if((++degree) > 3) break;
    } while(++fc != done);

    fc = this->incident_faces(v);
    done = fc;

    if(degree == 3) {
      do {
        int i = fc->index(v);
        edges.push_back(Edge(fc, i));
      } while(++fc != done);
    } else {
      do {
        int i = fc->index(v);
        edges.push_back(Edge(fc, i));
        edges.push_back(Edge(fc, this->cw(i)));
      } while(++fc != done);
    }

    while(!edges.empty()) {
      const Edge &e = edges.front();
      Face_handle f = e.first;
      int i = e.second;
      edges.pop_front();
      faces.insert(f);
      if(this->is_infinite(f->vertex(i))) continue;
      Face_handle fi = f->neighbor(i);
      int mi = this->_tds.mirror_index(f, i);
      Vertex_handle vm = this->_tds.mirror_vertex(f, i);
      if(this->is_infinite(vm)) continue;
      if(this->side_of_oriented_circle(f, vm->point()) == ON_POSITIVE_SIDE) {
        this->_tds.flip(f, i);
        edges.push_back(Edge(f, i));
        edges.push_back(Edge(f, this->cw(i)));
        edges.push_back(Edge(fi, this->cw(mi)));
        edges.push_back(Edge(fi, mi));
      }
    }
  }
};

template < class Gt, class Tds >
inline bool
Delaunay_triangulation_2<Gt,Tds>::
test_conflict(const Point  &p, Face_handle fh, bool strict) const
{
  if(! strict) {
    Oriented_side os = side_of_oriented_circle(fh,p,false);
    return os == ON_POSITIVE_SIDE;
  }
  // return true  if P is inside the circumcircle of fh
  // if fh is infinite, return true when p is in the positive
  // halfspace or on the boundary and in the  finite edge of fh
  Oriented_side os = side_of_oriented_circle(fh,p,true);
  if(os == ON_POSITIVE_SIDE)
    return true;

  if(os == ON_ORIENTED_BOUNDARY && this->is_infinite(fh)) {
    int i = fh->index(this->infinite_vertex());
    return collinear_between(fh->vertex(cw(i))->point(), p,
                             fh->vertex(ccw(i))->point());
  }

  return false;
}

template < class Gt, class Tds >
inline bool
Delaunay_triangulation_2<Gt,Tds>::
does_conflict(const Point  &p, Face_handle fh) const
{
  return test_conflict(p,fh);
}

template < class Gt, class Tds >
inline bool
Delaunay_triangulation_2<Gt,Tds>::
find_conflicts(const Point  &p,
               std::list<Face_handle>& conflicts,
               Face_handle start) const
{
  get_conflicts(p, std::back_inserter(conflicts), start);
  return (! conflicts.empty());
}

template < class Gt, class Tds >
bool
Delaunay_triangulation_2<Gt,Tds>::
is_valid(bool verbose, int level) const
{
  bool result = Triangulation_2<Gt,Tds>::is_valid(verbose, level);

  for(Finite_faces_iterator it = this->finite_faces_begin();
       it != this->finite_faces_end(); it++) {
    for(int i=0; i<3; i++) {
      if(! this->is_infinite(this->mirror_vertex(it,i))) {
        result = result && ON_POSITIVE_SIDE !=
            side_of_oriented_circle(it, this->mirror_vertex(it,i)->point(), false);
      }
      CGAL_triangulation_assertion(result);
    }
  }
  return result;
}

template < class Gt, class Tds >
typename Delaunay_triangulation_2<Gt,Tds>::Vertex_handle
Delaunay_triangulation_2<Gt,Tds>::
nearest_vertex(const Point  &p, Face_handle f) const
{
  switch (this->dimension()) {
  case 0:
    if(this->number_of_vertices() == 0) return Vertex_handle();
    return this->finite_vertex();
    //break;
  case 1:
    return nearest_vertex_1D(p);
    //break;
  case 2:
    return nearest_vertex_2D(p,f);
    //break;
  }
  return Vertex_handle();
}

template < class Gt, class Tds >
typename Delaunay_triangulation_2<Gt,Tds>::Vertex_handle
Delaunay_triangulation_2<Gt,Tds>::
nearest_vertex_2D(const Point& p, Face_handle f) const
{
  CGAL_triangulation_precondition(this->dimension() == 2);
  f = this->locate(p,f);

  typename Geom_traits::Compare_distance_2
    compare_distance =  this->geom_traits().compare_distance_2_object();
  Vertex_handle nn = !this->is_infinite(f->vertex(0)) ? f->vertex(0):f->vertex(1);
  if(!this->is_infinite(f->vertex(1)) && compare_distance(p,
              f->vertex(1)->point(),
              nn->point()) == SMALLER)
    nn = f->vertex(1);
  if(!this->is_infinite(f->vertex(2)) && compare_distance(p,
              f->vertex(2)->point(),
              nn->point()) == SMALLER)
    nn = f->vertex(2);

  look_nearest_neighbor(p,f,0,nn);
  look_nearest_neighbor(p,f,1,nn);
  look_nearest_neighbor(p,f,2,nn);
  return nn;
}

template < class Gt, class Tds >
typename Delaunay_triangulation_2<Gt,Tds>::Vertex_handle
Delaunay_triangulation_2<Gt,Tds>::
nearest_vertex_1D(const Point& p) const
{
  typename Geom_traits::Compare_distance_2
      compare_distance =  this->geom_traits().compare_distance_2_object();
  Vertex_handle nn;

  Finite_vertices_iterator vit=this->finite_vertices_begin();
  nn = vit;
  for (; vit != this->finite_vertices_end(); ++vit) {
    if(compare_distance(p, vit->point(), nn->point()) == SMALLER)
      nn = vit;
  }
  return nn;
}

template < class Gt, class Tds >
void
Delaunay_triangulation_2<Gt,Tds>::
look_nearest_neighbor(const Point& p,
                      Face_handle f,
                      int i,
                      Vertex_handle& nn) const
{
  Face_handle  ni=f->neighbor(i);
  if(ON_POSITIVE_SIDE != side_of_oriented_circle(ni,p,true))
    return;

  typename Geom_traits::Compare_distance_2
    compare_distance = this->geom_traits().compare_distance_2_object();
  i = ni->index(f);
  if(!this->is_infinite(ni->vertex(i)) &&
       compare_distance(p, ni->vertex(i)->point(), nn->point())  == SMALLER)
    nn = ni->vertex(i);

  // recursive exploration of triangles whose circumcircle contains p
  look_nearest_neighbor(p, ni, ccw(i), nn);
  look_nearest_neighbor(p, ni, cw(i), nn);
}

//DUALITY
template<class Gt, class Tds>
inline
typename Delaunay_triangulation_2<Gt,Tds>::Point
Delaunay_triangulation_2<Gt,Tds>::
dual(Face_handle f) const
{
  CGAL_triangulation_precondition(this->_tds.is_face(f));
  CGAL_triangulation_precondition(this->dimension()==2);
  return circumcenter(f);
}

template < class Gt, class Tds >
Object
Delaunay_triangulation_2<Gt,Tds>::
dual(const Edge &e) const
{
  CGAL_triangulation_precondition(this->_tds.is_edge(e.first,e.second));

  typedef typename Geom_traits::Line_2        Line;
  typedef typename Geom_traits::Ray_2         Ray;

  CGAL_triangulation_precondition (!this->is_infinite(e));
  if(this->dimension()== 1) {
    const Point& p = (e.first)->vertex(cw(e.second))->point();
    const Point& q = (e.first)->vertex(ccw(e.second))->point();
    Line l  = this->geom_traits().construct_bisector_2_object()(p,q);
    return make_object(l);
  }

  // dimension==2
  if((!this->is_infinite(e.first)) &&
      (!this->is_infinite(e.first->neighbor(e.second)))) {
    Segment s = this->geom_traits().construct_segment_2_object()
                          (dual(e.first),dual(e.first->neighbor(e.second)));
    return make_object(s);
  }

  // one of the adjacent faces is infinite
  Face_handle f; int i;
  if(this->is_infinite(e.first)) {
    f = e.first->neighbor(e.second); i=f->index(e.first);
  }
  else {
    f = e.first; i=e.second;
  }

  const Point& p = f->vertex(cw(i))->point();
  const Point& q = f->vertex(ccw(i))->point();
  Line l = this->geom_traits().construct_bisector_2_object()(p,q);
  Ray r = this->geom_traits().construct_ray_2_object()(dual(f), l);
  return make_object(r);
}

template < class Gt, class Tds >
inline Object
Delaunay_triangulation_2<Gt,Tds>::
dual(const Edge_circulator& ec) const
{
  return dual(*ec);
}

template < class Gt, class Tds >
inline Object
Delaunay_triangulation_2<Gt,Tds>::
dual(const Finite_edges_iterator& ei) const
{
  return dual(*ei);
}

///////////////////////////////////////////////////////////////
//  INSERT

template < class Gt, class Tds >
inline
typename Delaunay_triangulation_2<Gt,Tds>::Vertex_handle
Delaunay_triangulation_2<Gt,Tds>::
insert(const Point  &p, Face_handle start)
{
  Locate_type lt;
  int li;
  Face_handle loc = this->locate(p, lt, li, start);
  return insert(p, lt, loc, li);
}

template < class Gt, class Tds >
inline
typename Delaunay_triangulation_2<Gt,Tds>::Vertex_handle
Delaunay_triangulation_2<Gt,Tds>::
push_back(const Point &p)
{
  return insert(p);
}

template < class Gt, class Tds >
inline
typename Delaunay_triangulation_2<Gt,Tds>::Vertex_handle
Delaunay_triangulation_2<Gt,Tds>::
insert(const Point  &p, Locate_type lt, Face_handle loc, int li)
{
  Vertex_handle v = Triangulation_2<Gt,Tds>::insert(p,lt,loc,li);
  restore_Delaunay(v);
  return(v);
}

template < class Gt, class Tds >
template < class OutputItFaces >
inline
typename Delaunay_triangulation_2<Gt,Tds>::Vertex_handle
Delaunay_triangulation_2<Gt,Tds>::
insert_and_give_new_faces(const Point  &p,
                          OutputItFaces oif,
                          Face_handle start)
{
  Vertex_handle v = insert(p, start);
  int dimension = this->dimension();
  if(dimension == 2)
  {
    Face_circulator fc = this->incident_faces(v), done(fc);
    do {
      *oif++ = fc;
    } while(++fc != done);
  }
  else if(dimension == 1)
  {
    Face_handle c = v->face();
    *oif++ = c;
    *oif++ = c->neighbor((~(c->index(v)))&1);
  }
  else *oif++ = v->face(); // dimension == 0
  return v;
}

template < class Gt, class Tds >
template < class OutputItFaces >
inline
typename Delaunay_triangulation_2<Gt,Tds>::Vertex_handle
Delaunay_triangulation_2<Gt,Tds>::
insert_and_give_new_faces(const Point  &p,
                          Locate_type lt,
                          Face_handle loc, int li,
                          OutputItFaces oif)
{
  Vertex_handle v = insert(p, lt, loc, li);
  int dimension = this->dimension();
  if(dimension == 2)
  {
    Face_circulator fc = this->incident_faces(v), done(fc);
    do {
      *oif++ = fc;
    } while(++fc != done);
  }
  else if(dimension == 1)
  {
    Face_handle c = v->face();
    *oif++ = c;
    *oif++ = c->neighbor((~(c->index(v)))&1);
  }
  else *oif++ = v->face(); // dimension == 0
  return v;
}

template < class Gt, class Tds >
void
Delaunay_triangulation_2<Gt,Tds>::
restore_Delaunay(Vertex_handle v)
{
  if(this->dimension() <= 1) return;

  Face_handle f=v->face();
  Face_handle next;
  int i;
  Face_handle start(f);
  do {
    i = f->index(v);
    next = f->neighbor(ccw(i));  // turn ccw around v
    propagating_flip(f,i);
    f = next;
  } while(next != start);
  return;
}

#ifndef CGAL_DT2_USE_RECURSIVE_PROPAGATING_FLIP
template < class Gt, class Tds >
void
Delaunay_triangulation_2<Gt,Tds>::
non_recursive_propagating_flip(Face_handle f, int i)
{
  std::stack<Edge> edges;
  const Vertex_handle& vp = f->vertex(i);
  const Point& p = vp->point();
  edges.push(Edge(f,i));

  while(! edges.empty()) {
    const Edge& e = edges.top();
    f = e.first;
    i = e.second;
    const Face_handle& n = f->neighbor(i);

    if(ON_POSITIVE_SIDE != side_of_oriented_circle(n, p, true)) {
      edges.pop();
      continue;
    }

    this->flip(f, i);
    // As we haven't popped it, we don't have to push it
    edges.push(Edge(n,n->index(vp)));
  }
}

template < class Gt, class Tds >
void
Delaunay_triangulation_2<Gt,Tds>::
propagating_flip(const Face_handle& f,int i, int depth)
{
#ifdef CGAL_DT2_IMMEDIATELY_NON_RECURSIVE_PROPAGATING_FLIP
  non_recursive_propagating_flip(f,i);
#else
  int max_depth = 100;
  if(depth==max_depth) {
    non_recursive_propagating_flip(f,i);
    return;
  }
  Face_handle n = f->neighbor(i);

  if(ON_POSITIVE_SIDE !=
      side_of_oriented_circle(n, f->vertex(i)->point(), true)) {
    return;
  }

  this->flip(f, i);
  propagating_flip(f,i,depth+1);
  i = n->index(f->vertex(i));
  propagating_flip(n,i,depth+1);
#endif
}
#else
template < class Gt, class Tds >
void
Delaunay_triangulation_2<Gt,Tds>::
propagating_flip(const Face_handle& f,int i)
{
  Face_handle n = f->neighbor(i);

  if(ON_POSITIVE_SIDE != side_of_oriented_circle(n, f->vertex(i)->point(), true))
    return;

  this->flip(f, i);
  propagating_flip(f,i);
  i = n->index(f->vertex(i));
  propagating_flip(n,i);
}
#endif

///////////////////////////////////////////////////////////////
//  REMOVE    see INRIA RResearch Report 7104

template < class Gt, class Tds >
template <class OutputItFaces>
void
Delaunay_triangulation_2<Gt,Tds>::
remove_and_give_new_faces(Vertex_handle v, OutputItFaces fit)
{
  CGAL_triangulation_precondition(v != Vertex_handle());
  CGAL_triangulation_precondition(!this->is_infinite(v));

  if(this->number_of_vertices() == 1) this->remove_first(v);
  else if(this->number_of_vertices() == 2) this->remove_second(v);
  else if(this->dimension() == 1)
  {
    Point p = v->point();
    Triangulation::remove(v);
    *fit++ = this->locate(p);
  }
  else if(this->test_dim_down(v)) {
    this->_tds.remove_dim_down(v);
    for(All_faces_iterator afi = this-> all_faces_begin();
        afi != this->all_faces_end();
        afi++) *fit++ = afi;
  }
  else {
    CGAL_STATIC_THREAD_LOCAL_VARIABLE(int, maxd,30);
    CGAL_STATIC_THREAD_LOCAL_VARIABLE(std::vector<Face_handle>, f, maxd);
    CGAL_STATIC_THREAD_LOCAL_VARIABLE(std::vector<int>, i, maxd);
    CGAL_STATIC_THREAD_LOCAL_VARIABLE(std::vector<Vertex_handle>, w, maxd);

    int d;
    remove_degree_init(v,f,w,i,d,maxd);
    remove_degree_triangulate(v,f,w,i,d);
    this->delete_vertex(v);
    Face_circulator fc(v[0]),done;
    do *fit++ = fc++; while (fc!=done);
  }
  return;
}

template < class Gt, class Tds >
void
Delaunay_triangulation_2<Gt,Tds>::
remove(Vertex_handle v)
{
  int d;

  CGAL_triangulation_precondition(v != Vertex_handle());
  CGAL_triangulation_precondition(!this->is_infinite(v));

  if(this->dimension() <= 1) { Triangulation::remove(v); return; }

  CGAL_STATIC_THREAD_LOCAL_VARIABLE(int, maxd,30);
  CGAL_STATIC_THREAD_LOCAL_VARIABLE(std::vector<Face_handle>, f, maxd);
  CGAL_STATIC_THREAD_LOCAL_VARIABLE(std::vector<int>, i, maxd);
  CGAL_STATIC_THREAD_LOCAL_VARIABLE(std::vector<Vertex_handle>, w, maxd);

  remove_degree_init(v,f,w,i,d,maxd);
  if(d == 0) return; //  dim is going down
  remove_degree_triangulate(v,f,w,i,d);
  this->delete_vertex(v);
}

template < class Gt, class Tds >
void
Delaunay_triangulation_2<Gt,Tds>::
remove_degree_init(Vertex_handle v, std::vector<Face_handle> &f,
                   std::vector<Vertex_handle> &w, std::vector<int> &i,
                   int &d, int &maxd)
{
  f[0] = v->face();d=0;
  do{
    i[d] = f[d]->index(v);
    w[d] = f[d]->vertex(ccw(i[d]));
    if(this->is_infinite(w[d])) {
      f[0] = f[d]; i[0]=i[d]; w[0]=w[d];
      w[0]->set_face(f[0]->neighbor(i[0]));
      f[1] = f[0]->neighbor(ccw(i[0]));
      i[1] = f[1]->index(v);
      w[1] = f[1]->vertex(ccw(i[1]));
      if(this->is_infinite(f[1]->neighbor(i[1]))) {//otherwise dim remains 2
        if(this->test_dim_down(v)) {
          d=0;
          this->tds().remove_dim_down(v);
          return;
        }
      }
      d=1;
    }
    w[d]->set_face(f[d]->neighbor(i[d]));//do no longer bother about set_face
    ++d;
    if(d==maxd) { maxd *=2; f.resize(maxd); w.resize(maxd); i.resize(maxd);}
    f[d] = f[d-1]->neighbor(ccw(i[d-1]));
  } while(f[d]!=f[0]);
  // all vertices finite but possibly w[0]
}

template < class Gt, class Tds >
void
Delaunay_triangulation_2<Gt,Tds>::
remove_degree_triangulate(Vertex_handle v,
                          std::vector<Face_handle> &f,
                          std::vector<Vertex_handle> &w,
                          std::vector<int> &i,int d)
{
  switch (d) {
  case 3:
    remove_degree3(v,f,w,i);    break;
  case 4:
    remove_degree4(v,f,w,i);    break;
  case 5:
    remove_degree5(v,f,w,i);    break;
  case 6:
    remove_degree6(v,f,w,i);    break;
  case 7:
    remove_degree7(v,f,w,i);    break;
  default:
    remove_degree_d(v,f,w,i,d);    break;
  }
}

template < class Gt, class Tds >
void
Delaunay_triangulation_2<Gt,Tds>::
remove_degree_d(Vertex_handle v, std::vector<Face_handle> &,
                std::vector<Vertex_handle> &,
                std::vector<int> &,int)
{
  // removing a degree d vertex, (dim is not going down)
  // this is the old removal procedure that is used now only if d > 7

  std::list<Edge> hole;
  make_hole(v, hole);
  fill_hole_delaunay(hole);
  return;
}

template < class Gt, class Tds >
void
Delaunay_triangulation_2<Gt,Tds>::
remove_degree3(Vertex_handle, std::vector<Face_handle> &f,
               std::vector<Vertex_handle> &, std::vector<int> &i)
{
  // removing a degree 3 vertex
  // only w[0] can be infinite

  // modify the triangulation
  Face_handle nn = f[1]->neighbor(i[1]);
  this->tds().set_adjacency(f[0], ccw(i[0]), nn, nn->index(f[1]));
  nn = f[2]->neighbor(i[2]);
  this->tds().set_adjacency(f[0], cw(i[0]), nn, nn->index(f[2]));
  f[0]->set_vertex(i[0], f[1]->vertex(cw(i[1])));

  // clean container
  this->tds().delete_face(f[1]);
  this->tds().delete_face(f[2]);

  return;
}

template < class Gt, class Tds >
void
Delaunay_triangulation_2<Gt,Tds>::
remove_degree4(Vertex_handle, std::vector<Face_handle> &f,
               std::vector<Vertex_handle> &w, std::vector<int> &i)
{
  // removing a degree 4 vertex
  // only w[0] can be infinite

  Face_handle nn;
  // modify f[0] f[1] for incircle test
  f[0]->set_vertex(i[0], w[3]); //w0 w1 w3

  if(!test_conflict(w[2]->point(), f[0])) {
    // diagonal 1 3
    f[1]->set_vertex(i[1], w[3]); //w1 w2 w3
    nn = f[3]->neighbor(i[3]);
    this->tds().set_adjacency(f[0], cw(i[0]), nn, nn->index(f[3]));
    nn = f[2]->neighbor(i[2]);
    this->tds().set_adjacency(f[1], ccw(i[1]), nn, nn->index(f[2]));
    // clean container
    this->tds().delete_face(f[2]);
    this->tds().delete_face(f[3]);
  } else {
    // diagonal 0 2
    f[0]->set_vertex(i[0], w[2]); //w0 w1 w2
    f[3]->set_vertex(i[3], w[2]); //w3 w0 w2
    nn = f[1]->neighbor(i[1]);
    this->tds().set_adjacency(f[0], ccw(i[0]), nn, nn->index(f[1]));
    nn = f[2]->neighbor(i[2]);
    this->tds().set_adjacency(f[3], cw(i[3]), nn, nn->index(f[2]));
    // clean container
    this->tds().delete_face(f[1]);
    this->tds().delete_face(f[2]);
  }

  return;
}

template < class Gt, class Tds >
void
Delaunay_triangulation_2<Gt,Tds>::
remove_degree5(Vertex_handle v, std::vector<Face_handle> &f,
               std::vector<Vertex_handle> &w, std::vector<int> &i)
{
  // removing a degree 5 vertex
  // only w[0] can be infinite

  if(incircle(3,0,1,2,f,w,i)) {
    if(incircle(4,0,1,3,f,w,i)) {
      if(incircle(4,1,2,3,f,w,i)) {
        // star from 4
        remove_degree5_star(v,f[4],f[0],f[1],f[2],f[3],
                            w[4],w[0],w[1],w[2],w[3],
                            i[4],i[0],i[1],i[2],i[3]);
      } else {
        //star from 1
        remove_degree5_star(v,f[1],f[2],f[3],f[4],f[0],
                            w[1],w[2],w[3],w[4],w[0],
                            i[1],i[2],i[3],i[4],i[0]);
      }
    } else {
      // star from 3
      remove_degree5_star(v,f[3],f[4],f[0],f[1],f[2],
                          w[3],w[4],w[0],w[1],w[2],
                          i[3],i[4],i[0],i[1],i[2]);
    }
  } else {
    if(incircle(4,2,3,0,f,w,i)) {
      if(incircle(4,0,1,2,f,w,i)) {
        // star from 4
        remove_degree5_star(v,f[4],f[0],f[1],f[2],f[3],
                            w[4],w[0],w[1],w[2],w[3],
                            i[4],i[0],i[1],i[2],i[3]);
      } else {
        //star from 2
        remove_degree5_star(v,f[2],f[3],f[4],f[0],f[1],
                            w[2],w[3],w[4],w[0],w[1],
                            i[2],i[3],i[4],i[0],i[1]);
      }
    } else {
      // star from 0
      remove_degree5_star(v,f[0],f[1],f[2],f[3],f[4],
                          w[0],w[1],w[2],w[3],w[4],
                          i[0],i[1],i[2],i[3],i[4]);
    }
  }

  return;
}

template < class Gt, class Tds >
inline void
Delaunay_triangulation_2<Gt,Tds>::remove_degree5_star
(
 Vertex_handle &,
 Face_handle & f0, Face_handle & f1, Face_handle & f2,
 Face_handle & f3, Face_handle & f4,
 Vertex_handle &v0, Vertex_handle &, Vertex_handle &,
 Vertex_handle &, Vertex_handle &,
 int i0, int i1, int i2, int i3, int i4)
{ // removing a degree 5 vertex, staring from v0
  Face_handle nn;
  f1->set_vertex(i1, v0);  // f1 = v1v2v0
  f2->set_vertex(i2, v0);  // f2 = v2v3v0
  f3->set_vertex(i3, v0);  // f3 = v3v4v0
  nn = f0->neighbor(i0);
  this->tds().set_adjacency(f1, cw(i1), nn, nn->index(f0));
  nn = f4->neighbor(i4);
  this->tds().set_adjacency(f3, ccw(i3), nn, nn->index(f4));
  this->tds().delete_face(f0);
  this->tds().delete_face(f4);
}

template < class Gt, class Tds >
void
Delaunay_triangulation_2<Gt,Tds>::
remove_degree6(Vertex_handle v, std::vector<Face_handle> &f,
               std::vector<Vertex_handle> &w, std::vector<int> &i)
{
  // removing a degree 6 vertex
  // only w[0] can be infinite

  if(incircle(1,2,3,0,f,w,i)) {
    if(incircle(4,2,3,5,f,w,i)) {
      if(incircle(1,2,3,4,f,w,i)) {
        if(incircle(4,0,1,3,f,w,i)) {
          if(incircle(5,0,1,4,f,w,i)) {
            remove_degree6_star(v,f[1],f[2],f[3],f[4],f[5],f[0],
                                w[1],w[2],w[3],w[4],w[5],w[0],
                                i[1],i[2],i[3],i[4],i[5],i[0]);
          } else {
            remove_degree6_N(v,f[1],f[2],f[3],f[4],f[5],f[0],
                             w[1],w[2],w[3],w[4],w[5],w[0],
                             i[1],i[2],i[3],i[4],i[5],i[0]);
          }} else {
          remove_degree6_antiN(v,f[0],f[1],f[2],f[3],f[4],f[5],
                               w[0],w[1],w[2],w[3],w[4],w[5],
                               i[0],i[1],i[2],i[3],i[4],i[5]);
        }} else {
        if(incircle(5,1,2,4,f,w,i)) {
          remove_degree6_N(v,f[2],f[3],f[4],f[5],f[0],f[1],
                           w[2],w[3],w[4],w[5],w[0],w[1],
                           i[2],i[3],i[4],i[5],i[0],i[1]);
        } else {
          if(incircle(5,0,1,4,f,w,i)) {
            remove_degree6_antiN(v,f[1],f[2],f[3],f[4],f[5],f[0],
                                 w[1],w[2],w[3],w[4],w[5],w[0],
                                 i[1],i[2],i[3],i[4],i[5],i[0]);
          } else {
            remove_degree6_star(v,f[4],f[5],f[0],f[1],f[2],f[3],
                                w[4],w[5],w[0],w[1],w[2],w[3],
                                i[4],i[5],i[0],i[1],i[2],i[3]);
          }}}} else {
      if(incircle(1,2,3,5,f,w,i)) {
        if(incircle(1,3,4,5,f,w,i)) {
          if(incircle(4,0,1,3,f,w,i)) {
            if(incircle(5,0,1,4,f,w,i)) {
              remove_degree6_star(v,f[1],f[2],f[3],f[4],f[5],f[0],
                                  w[1],w[2],w[3],w[4],w[5],w[0],
                                  i[1],i[2],i[3],i[4],i[5],i[0]);
            } else {
              remove_degree6_N(v,f[1],f[2],f[3],f[4],f[5],f[0],
                               w[1],w[2],w[3],w[4],w[5],w[0],
                               i[1],i[2],i[3],i[4],i[5],i[0]);
            }} else {
            remove_degree6_antiN(v,f[0],f[1],f[2],f[3],f[4],f[5],
                                 w[0],w[1],w[2],w[3],w[4],w[5],
                                 i[0],i[1],i[2],i[3],i[4],i[5]);
          }} else {
          if(incircle(5,0,1,3,f,w,i)) {
            remove_degree6_diamond(v,f[1],f[2],f[3],f[4],f[5],f[0],
                                   w[1],w[2],w[3],w[4],w[5],w[0],
                                   i[1],i[2],i[3],i[4],i[5],i[0]);
          } else {
            if(incircle(4,5,0,3,f,w,i)) {
              remove_degree6_antiN(v,f[0],f[1],f[2],f[3],f[4],f[5],
                                   w[0],w[1],w[2],w[3],w[4],w[5],
                                   i[0],i[1],i[2],i[3],i[4],i[5]);
            } else {
              remove_degree6_star(v,f[3],f[4],f[5],f[0],f[1],f[2],
                                  w[3],w[4],w[5],w[0],w[1],w[2],
                                  i[3],i[4],i[5],i[0],i[1],i[2]);
            }}}} else {
        remove_degree6_star(v,f[5],f[0],f[1],f[2],f[3],f[4],
            w[5],w[0],w[1],w[2],w[3],w[4],
            i[5],i[0],i[1],i[2],i[3],i[4]);
      }}} else {
    if(incircle(4,2,3,5,f,w,i)) {
      if(incircle(4,2,3,0,f,w,i)) {
        if(incircle(4,0,1,2,f,w,i)) {
          if(incircle(4,1,2,5,f,w,i)) {
            if(incircle(4,0,1,5,f,w,i)) {
              remove_degree6_star(v,f[4],f[5],f[0],f[1],f[2],f[3],
                  w[4],w[5],w[0],w[1],w[2],w[3],
                  i[4],i[5],i[0],i[1],i[2],i[3]);
            } else {
              remove_degree6_antiN(v,f[1],f[2],f[3],f[4],f[5],f[0],
                  w[1],w[2],w[3],w[4],w[5],w[0],
                  i[1],i[2],i[3],i[4],i[5],i[0]);
            }} else {
            remove_degree6_N(v,f[2],f[3],f[4],f[5],f[0],f[1],
                w[2],w[3],w[4],w[5],w[0],w[1],
                i[2],i[3],i[4],i[5],i[0],i[1]);
          }} else {
          if(incircle(4,5,0,2,f,w,i)) {
            remove_degree6_diamond(v,f[0],f[1],f[2],f[3],f[4],f[5],
                                     w[0],w[1],w[2],w[3],w[4],w[5],
                                     i[0],i[1],i[2],i[3],i[4],i[5]);
          } else {
            if(incircle(5,0,1,2,f,w,i)) {
              remove_degree6_N(v,f[2],f[3],f[4],f[5],f[0],f[1],
                               w[2],w[3],w[4],w[5],w[0],w[1],
                               i[2],i[3],i[4],i[5],i[0],i[1]);
            } else {
              remove_degree6_star(v,f[2],f[3],f[4],f[5],f[0],f[1],
                                  w[2],w[3],w[4],w[5],w[0],w[1],
                                  i[2],i[3],i[4],i[5],i[0],i[1]);
            }}}} else {
        remove_degree6_star(v,f[0],f[1],f[2],f[3],f[4],f[5],
                            w[0],w[1],w[2],w[3],w[4],w[5],
                            i[0],i[1],i[2],i[3],i[4],i[5]);
      }} else {
      if(incircle(5,2,3,0,f,w,i)) {
        if(incircle(5,0,1,2,f,w,i)) {
          remove_degree6_star(v,f[5],f[0],f[1],f[2],f[3],f[4],
                              w[5],w[0],w[1],w[2],w[3],w[4],
                              i[5],i[0],i[1],i[2],i[3],i[4]);
        } else {
          remove_degree6_antiN(v,f[2],f[3],f[4],f[5],f[0],f[1],
                               w[2],w[3],w[4],w[5],w[0],w[1],
                               i[2],i[3],i[4],i[5],i[0],i[1]);
        }} else {
        if(incircle(4,5,0,3,f,w,i)) {
          remove_degree6_star(v,f[0],f[1],f[2],f[3],f[4],f[5],
                              w[0],w[1],w[2],w[3],w[4],w[5],
                              i[0],i[1],i[2],i[3],i[4],i[5]);
        } else {
          remove_degree6_N(v,f[0],f[1],f[2],f[3],f[4],f[5],
                           w[0],w[1],w[2],w[3],w[4],w[5],
                           i[0],i[1],i[2],i[3],i[4],i[5]);
        }
      }
    }
  }
}

template < class Gt, class Tds >
inline void
Delaunay_triangulation_2<Gt,Tds>::remove_degree6_star(
  Vertex_handle &,
  Face_handle & f0, Face_handle & f1, Face_handle & f2,
  Face_handle & f3, Face_handle & f4, Face_handle & f5,
  Vertex_handle &v0, Vertex_handle &, Vertex_handle &,
  Vertex_handle &, Vertex_handle &, Vertex_handle &,
  int i0, int i1, int i2, int i3, int i4, int i5)
{ // removing a degree 6 vertex, staring from v0
  Face_handle nn;
  f1->set_vertex(i1, v0);  // f1 = v1v2v0
  f2->set_vertex(i2, v0);  // f2 = v2v3v0
  f3->set_vertex(i3, v0);  // f3 = v3v4v0
  f4->set_vertex(i4, v0);  // f4 = v4v5v0
  nn = f0->neighbor(i0);
  this->tds().set_adjacency(f1, cw(i1), nn, nn->index(f0));
  nn = f5->neighbor(i5);
  this->tds().set_adjacency(f4, ccw(i4), nn, nn->index(f5));
  this->tds().delete_face(f0);
  this->tds().delete_face(f5);
}

template < class Gt, class Tds >
inline void
Delaunay_triangulation_2<Gt,Tds>::remove_degree6_N(
  Vertex_handle &,
  Face_handle & f0, Face_handle & f1, Face_handle & f2,
  Face_handle & f3, Face_handle & f4, Face_handle & f5,
  Vertex_handle &v0, Vertex_handle &, Vertex_handle &,
  Vertex_handle &v3, Vertex_handle &, Vertex_handle &,
  int i0, int i1, int i2, int i3, int i4, int i5)
{ // removing a degree 6 vertex, N configuration with diagonal v0v3
  Face_handle nn;
  f1->set_vertex(i1, v0);  // f1 = v1v2v0
  f2->set_vertex(i2, v0);  // f2 = v2v3v0
  f4->set_vertex(i4, v3);  // f4 = v4v5v3
  f5->set_vertex(i5, v3);  // f5 = v5v0v3
  nn = f0->neighbor(i0);
  this->tds().set_adjacency(f1, cw(i1), nn, nn->index(f0));
  nn = f3->neighbor(i3);
  this->tds().set_adjacency(f4, cw(i4), nn, nn->index(f3));
  this->tds().set_adjacency(f2, ccw(i2), f5, ccw(i5));
  this->tds().delete_face(f0);
  this->tds().delete_face(f3);
}

template < class Gt, class Tds >
inline void
Delaunay_triangulation_2<Gt,Tds>::remove_degree6_antiN(
  Vertex_handle &,
  Face_handle & f0, Face_handle & f1, Face_handle & f2,
  Face_handle & f3, Face_handle & f4, Face_handle & f5,
  Vertex_handle &v0, Vertex_handle &, Vertex_handle &,
  Vertex_handle &v3, Vertex_handle &, Vertex_handle &,
  int i0, int i1, int i2, int i3, int i4, int i5)
{ // removing a degree 6 vertex, antiN configuration with diagonal v0v3
  Face_handle nn;
  f0->set_vertex(i0, v3);  // f0 = v0v1v3
  f1->set_vertex(i1, v3);  // f1 = v1v2v3
  f3->set_vertex(i3, v0);  // f3 = v3v4v0
  f4->set_vertex(i4, v0);  // f4 = v4v5v0
  nn = f2->neighbor(i2);
  this->tds().set_adjacency(f1, ccw(i1), nn, nn->index(f2));
  nn = f5->neighbor(i5);
  this->tds().set_adjacency(f4, ccw(i4), nn, nn->index(f5));
  this->tds().set_adjacency(f0, cw(i0), f3, cw(i3));
  this->tds().delete_face(f2);
  this->tds().delete_face(f5);
}

template < class Gt, class Tds >
inline void
Delaunay_triangulation_2<Gt,Tds>::remove_degree6_diamond(
  Vertex_handle &,
  Face_handle & f0, Face_handle & f1, Face_handle & f2,
  Face_handle & f3, Face_handle & f4, Face_handle & f5,
  Vertex_handle &v0, Vertex_handle &, Vertex_handle &v2,
  Vertex_handle &, Vertex_handle &v4, Vertex_handle &,
  int i0, int i1, int i2, int i3, int i4, int i5)
{ // removing a degree 6 vertex, with chords v0v2 v2v4 v4v0
  Face_handle nn;
  f0->set_vertex(i0, v2);  // f0 = v0v1v2
  f2->set_vertex(i2, v4);  // f2 = v2v3v4
  f4->set_vertex(i4, v0);  // f4 = v4v5v0
  f1->set_vertex(i1, v4);
  f1->set_vertex(ccw(i1), v0);  // f1 = v0v2v4
  nn = f1->neighbor(i1);
  this->tds().set_adjacency(f0, ccw(i0), nn, nn->index(f1));
  nn = f3->neighbor(i3);
  this->tds().set_adjacency(f2, ccw(i2), nn, nn->index(f3));
  nn = f5->neighbor(i5);
  this->tds().set_adjacency(f4, ccw(i4), nn, nn->index(f5));
  this->tds().set_adjacency(f0, cw(i0), f1, i1);
  this->tds().set_adjacency(f4, cw(i4), f1, cw(i1));

  this->tds().delete_face(f3);
  this->tds().delete_face(f5);
}

template < class Gt, class Tds >
void
Delaunay_triangulation_2<Gt,Tds>::
remove_degree7(Vertex_handle v,std::vector<Face_handle> &f,
               std::vector<Vertex_handle> &w, std::vector<int> &i)
{
  // removing a degree 7 vertex
  // only w[0] can be infinite

  if(incircle(2,0,1,3,f,w,i)) { // sweeping from above
    if(incircle(2,3,4,0,f,w,i)) {
      if(incircle(5,3,4,6,f,w,i)) {
        if(incircle(5,3,4,2,f,w,i)) {
          if(incircle(6,2,3,5,f,w,i)) {
            if(incircle(6,0,1,2,f,w,i)) {
              remove_degree7_leftfan(v,6,f,w,i);
            } else {
              remove_degree7_zigzag(v,6,f,w,i);
            }} else {
            if(incircle(5,0,1,2,f,w,i)) {
              if(incircle(6,1,2,5,f,w,i)) {
                remove_degree7_zigzag(v,2,f,w,i);
              } else {
                if(incircle(6,0,1,5,f,w,i)) {
                  remove_degree7_rightfan(v,5,f,w,i);
                } else {
                  remove_degree7_star(v,5,f,w,i);
                }}} else {
              if(incircle(2,5,6,0,f,w,i)) {
                if(incircle(6,0,1,2,f,w,i)) {
                  remove_degree7_zigzag(v,2,f,w,i);
                } else {
                  remove_degree7_rightfan(v,2,f,w,i);
                }} else {
                remove_degree7_rightdelta(v,5,f,w,i);
              }}}} else {
          if(incircle(4,0,1,2,f,w,i)) {
            if(incircle(5,1,2,4,f,w,i)) {
              if(incircle(6,1,2,5,f,w,i)) {
                remove_degree7_leftfan(v,2,f,w,i);
              } else {
                if(incircle(6,0,1,5,f,w,i)) {
                  remove_degree7_zigzag(v,5,f,w,i);
                } else {
                  remove_degree7_leftfan(v,5,f,w,i);
                }}} else {
              if(incircle(5,0,1,4,f,w,i)) {
                if(incircle(6,0,1,5,f,w,i)) {
                  remove_degree7_rightfan(v,1,f,w,i);
                } else {
                  remove_degree7_zigzag(v,1,f,w,i);
                }} else {
                remove_degree7_rightfan(v,4,f,w,i);
              }}} else {
            if(incircle(2,4,5,0,f,w,i)) {
              if(incircle(5,0,1,2,f,w,i)) {
                if(incircle(6,1,2,5,f,w,i)) {
                  remove_degree7_leftfan(v,2,f,w,i);
                } else {
                  if(incircle(6,0,1,5,f,w,i)) {
                    remove_degree7_zigzag(v,5,f,w,i);
                  } else {
                    remove_degree7_leftfan(v,5,f,w,i);
                  }}} else {
                if(incircle(2,5,6,0,f,w,i)) {
                  if(incircle(6,0,1,2,f,w,i)) {
                    remove_degree7_leftfan(v,2,f,w,i);
                  } else {
                    remove_degree7_star(v,2,f,w,i);
                  }} else {
                  remove_degree7_leftdelta(v,2,f,w,i);
                }}} else {
              remove_degree7_rightdelta(v,0,f,w,i);
            }}}} else {
        if(incircle(6,3,4,2,f,w,i)) {
          if(incircle(6,0,1,2,f,w,i)) {
            remove_degree7_star(v,6,f,w,i);
          } else {
            remove_degree7_rightfan(v,6,f,w,i);
          }} else {
          if(incircle(4,0,1,2,f,w,i)) {
            if(incircle(2,4,5,6,f,w,i)) {
              if(incircle(5,1,2,4,f,w,i)) {
                if(incircle(6,1,2,5,f,w,i)) {
                  remove_degree7_leftfan(v,2,f,w,i);
                } else {
                  if(incircle(6,0,1,5,f,w,i)) {
                    remove_degree7_zigzag(v,5,f,w,i);
                  } else {
                    remove_degree7_leftfan(v,5,f,w,i);
                  }}} else {
                if(incircle(5,0,1,4,f,w,i)) {
                  if(incircle(6,0,1,5,f,w,i)) {
                    remove_degree7_rightfan(v,1,f,w,i);
                  } else {
                    remove_degree7_zigzag(v,1,f,w,i);
                  }} else{
                  remove_degree7_rightfan(v,4,f,w,i);
                }}} else {
              if(incircle(6,1,2,4,f,w,i)) {
                remove_degree7_leftdelta(v,6,f,w,i);
              } else {
                if(incircle(1,4,5,6,f,w,i)) {
                  if(incircle(1,4,5,0,f,w,i)) {
                    if(incircle(6,0,1,5,f,w,i)) {
                      remove_degree7_rightfan(v,1,f,w,i);
                    } else {
                      remove_degree7_zigzag(v,1,f,w,i);
                    }} else {
                    remove_degree7_rightfan(v,4,f,w,i);
                  }} else {
                  if(incircle(6,0,1,4,f,w,i)) {
                    remove_degree7_rightdelta(v,4,f,w,i);
                  } else {
                    if(incircle(6,4,5,0,f,w,i)) {
                      remove_degree7_star(v,4,f,w,i);
                    } else {
                      remove_degree7_rightfan(v,4,f,w,i);
                    }}}}}} else {
            if(incircle(2,4,5,6,f,w,i)) {
              if(incircle(2,4,5,0,f,w,i)) {
                if(incircle(5,0,1,2,f,w,i)) {
                  if(incircle(6,1,2,5,f,w,i)) {
                    remove_degree7_leftfan(v,2,f,w,i);
                  } else {
                    if(incircle(6,0,1,5,f,w,i)) {
                      remove_degree7_zigzag(v,5,f,w,i);
                    } else {
                      remove_degree7_leftfan(v,5,f,w,i);
                    }}} else {
                  if(incircle(2,5,6,0,f,w,i)) {
                    if(incircle(6,0,1,2,f,w,i)) {
                      remove_degree7_leftfan(v,2,f,w,i);
                    } else {
                      remove_degree7_star(v,2,f,w,i);
                    }} else {
                    remove_degree7_leftdelta(v,2,f,w,i);
                  }}} else {
                remove_degree7_rightdelta(v,0,f,w,i);
              }} else {
              if(incircle(2,6,0,4,f,w,i)) {
                if(incircle(6,0,1,2,f,w,i)) {
                  remove_degree7_leftdelta(v,6,f,w,i);
                } else {
                  remove_degree7_rightdelta(v,2,f,w,i);
                }} else {
                if(incircle(6,4,5,0,f,w,i)) {
                  remove_degree7_leftdelta(v,4,f,w,i);
                } else {
                  remove_degree7_rightdelta(v,0,f,w,i);
                }}}}}}} else{
      if(incircle(5,3,4,6,f,w,i)) {
        if(incircle(5,3,4,0,f,w,i)) {
          if(incircle(5,2,3,0,f,w,i)) {
            if(incircle(6,2,3,5,f,w,i)) {
              if(incircle(6,0,1,2,f,w,i)) {
                remove_degree7_leftfan(v,6,f,w,i);
              } else {
                remove_degree7_zigzag(v,6,f,w,i);
              }}else
              if(incircle(5,0,1,2,f,w,i)) {
                if(incircle(6,1,2,5,f,w,i)) {
                  remove_degree7_zigzag(v,2,f,w,i);
                } else {
                  if(incircle(6,0,1,5,f,w,i)) {
                    remove_degree7_rightfan(v,5,f,w,i);
                  } else {
                    remove_degree7_star(v,5,f,w,i);
                  }}} else {
                if(incircle(2,5,6,0,f,w,i)) {
                  if(incircle(6,0,1,2,f,w,i)) {
                    remove_degree7_zigzag(v,2,f,w,i);
                  } else {
                    remove_degree7_rightfan(v,2,f,w,i);
                  }} else {
                  remove_degree7_rightdelta(v,5,f,w,i);
                }}} else {
            if(incircle(3,5,6,0,f,w,i)) {
              if(incircle(6,2,3,0,f,w,i)) {
                if(incircle(6,0,1,2,f,w,i)) {
                  remove_degree7_leftfan(v,6,f,w,i);
                } else {
                  remove_degree7_zigzag(v,6,f,w,i);
                }} else {
                remove_degree7_leftfan(v,3,f,w,i);
              }} else {
              remove_degree7_leftdelta(v,0,f,w,i);
            }}} else {
          remove_degree7_star(v,0,f,w,i);
        }} else {
        if(incircle(6,3,4,0,f,w,i)) {
          if(incircle(6,2,3,0,f,w,i)) {
            if(incircle(6,0,1,2,f,w,i)) {
              remove_degree7_star(v,6,f,w,i);
            } else {
              remove_degree7_rightfan(v,6,f,w,i);
            }} else {
            remove_degree7_zigzag(v,3,f,w,i);
          }} else {
          if(incircle(6,4,5,0,f,w,i)) {
            remove_degree7_leftfan(v,0,f,w,i);
          } else {
            remove_degree7_star(v,0,f,w,i);
          }}}}} else {  //sweeping from below
    if(incircle(1,6,0,3,f,w,i)) {
      if(incircle(5,6,0,4,f,w,i)) {
        if(incircle(5,6,0,1,f,w,i)) {
          if(incircle(4,0,1,5,f,w,i)) {
            if(incircle(4,2,3,1,f,w,i)) {
              remove_degree7_rightfan(v,4,f,w,i);
            } else {
              remove_degree7_zigzag(v,4,f,w,i);
            }} else {
            if(incircle(5,2,3,1,f,w,i)) {
              if(incircle(4,1,2,5,f,w,i)) {
                remove_degree7_zigzag(v,1,f,w,i);
              } else {
                if(incircle(4,2,3,5,f,w,i)) {
                  remove_degree7_leftfan(v,5,f,w,i);
                } else {
                  remove_degree7_star(v,5,f,w,i);
                }}} else {
              if(incircle(1,4,5,3,f,w,i)) {
                if(incircle(4,2,3,1,f,w,i)) {
                  remove_degree7_zigzag(v,1,f,w,i);
                } else {
                  remove_degree7_leftfan(v,1,f,w,i);
                }} else {
                remove_degree7_leftdelta(v,5,f,w,i);
              }}}} else {
          if(incircle(6,2,3,1,f,w,i)) {
            if(incircle(5,1,2,6,f,w,i)) {
              if(incircle(4,1,2,5,f,w,i)) {
                remove_degree7_rightfan(v,1,f,w,i);
              } else {
                if(incircle(4,2,3,5,f,w,i)) {
                  remove_degree7_zigzag(v,5,f,w,i);
                } else {
                  remove_degree7_rightfan(v,5,f,w,i);
                }}} else {
              if(incircle(5,2,3,6,f,w,i)) {
                if(incircle(4,2,3,5,f,w,i)) {
                  remove_degree7_leftfan(v,2,f,w,i);
                } else {
                  remove_degree7_zigzag(v,2,f,w,i);
                }} else {
                remove_degree7_leftfan(v,6,f,w,i);
              }}} else {
            if(incircle(1,5,6,3,f,w,i)) {
              if(incircle(5,2,3,1,f,w,i)) {
                if(incircle(4,1,2,5,f,w,i)) {
                  remove_degree7_rightfan(v,1,f,w,i);
                } else {
                  if(incircle(4,2,3,5,f,w,i)) {
                    remove_degree7_zigzag(v,5,f,w,i);
                  } else {
                    remove_degree7_rightfan(v,5,f,w,i);
                  }}} else {
                if(incircle(1,4,5,3,f,w,i)) {
                  if(incircle(4,2,3,1,f,w,i)) {
                    remove_degree7_rightfan(v,1,f,w,i);
                  } else {
                    remove_degree7_star(v,1,f,w,i);
                  }} else {
                  remove_degree7_rightdelta(v,1,f,w,i);
                }}} else {
              remove_degree7_leftdelta(v,3,f,w,i);
            }}}} else {
        if(incircle(4,6,0,1,f,w,i)) {
          if(incircle(4,2,3,1,f,w,i)) {
            remove_degree7_star(v,4,f,w,i);
          } else {
            remove_degree7_leftfan(v,4,f,w,i);
          }} else {
          if(incircle(6,2,3,1,f,w,i)) {
            if(incircle(1,5,6,4,f,w,i)) {
              if(incircle(5,1,2,6,f,w,i)) {
                if(incircle(4,1,2,5,f,w,i)) {
                  remove_degree7_rightfan(v,1,f,w,i);
                } else {
                  if(incircle(4,2,3,5,f,w,i)) {
                    remove_degree7_zigzag(v,5,f,w,i);
                  } else {
                    remove_degree7_rightfan(v,5,f,w,i);
                  }}} else {
                if(incircle(5,2,3,6,f,w,i)) {
                  if(incircle(4,2,3,5,f,w,i)) {
                    remove_degree7_leftfan(v,2,f,w,i);
                  } else {
                    remove_degree7_zigzag(v,2,f,w,i);
                  }} else{
                  remove_degree7_leftfan(v,6,f,w,i);
                }}} else {
              if(incircle(4,1,2,6,f,w,i)) {
                remove_degree7_rightdelta(v,4,f,w,i);
              } else {
                if(incircle(2,5,6,4,f,w,i)) {
                  if(incircle(2,5,6,3,f,w,i)) {
                    if(incircle(4,2,3,5,f,w,i)) {
                      remove_degree7_leftfan(v,2,f,w,i);
                    } else {
                      remove_degree7_zigzag(v,2,f,w,i);
                    }} else {
                    remove_degree7_leftfan(v,6,f,w,i);
                  }} else {
                  if(incircle(4,2,3,6,f,w,i)) {
                    remove_degree7_leftdelta(v,6,f,w,i);
                  } else {
                    if(incircle(4,5,6,3,f,w,i)) {
                      remove_degree7_star(v,6,f,w,i);
                    } else {
                      remove_degree7_leftfan(v,6,f,w,i);
                    }}}}}} else {
            if(incircle(1,5,6,4,f,w,i)) {
              if(incircle(1,5,6,3,f,w,i)) {
                if(incircle(5,2,3,1,f,w,i)) {
                  if(incircle(4,1,2,5,f,w,i)) {
                    remove_degree7_rightfan(v,1,f,w,i);
                  } else {
                    if(incircle(4,2,3,5,f,w,i)) {
                      remove_degree7_zigzag(v,5,f,w,i);
                    } else {
                      remove_degree7_rightfan(v,5,f,w,i);
                    }}} else {
                  if(incircle(1,4,5,3,f,w,i)) {
                    if(incircle(4,2,3,1,f,w,i)) {
                      remove_degree7_rightfan(v,1,f,w,i);
                    } else {
                      remove_degree7_star(v,1,f,w,i);
                    }} else {
                    remove_degree7_rightdelta(v,1,f,w,i);
                  }}} else {
                remove_degree7_leftdelta(v,3,f,w,i);
              }} else {
              if(incircle(1,3,4,6,f,w,i)) {
                if(incircle(4,2,3,1,f,w,i)) {
                  remove_degree7_rightdelta(v,4,f,w,i);
                } else {
                  remove_degree7_leftdelta(v,1,f,w,i);
                }} else {
                if(incircle(4,5,6,3,f,w,i)) {
                  remove_degree7_rightdelta(v,6,f,w,i);
                } else {
                  remove_degree7_leftdelta(v,3,f,w,i);
                }}}}}}} else{
      if(incircle(5,6,0,4,f,w,i)) {
        if(incircle(5,6,0,3,f,w,i)) {
          if(incircle(5,0,1,3,f,w,i)) {
            if(incircle(4,0,1,5,f,w,i)) {
              if(incircle(4,2,3,1,f,w,i)) {
                remove_degree7_rightfan(v,4,f,w,i);
              } else {
                remove_degree7_zigzag(v,4,f,w,i);
              }}else
              if(incircle(5,2,3,1,f,w,i)) {
                if(incircle(4,1,2,5,f,w,i)) {
                  remove_degree7_zigzag(v,1,f,w,i);
                } else {
                  if(incircle(4,2,3,5,f,w,i)) {
                    remove_degree7_leftfan(v,5,f,w,i);
                  } else {
                    remove_degree7_star(v,5,f,w,i);
                  }}} else {
                if(incircle(1,4,5,3,f,w,i)) {
                  if(incircle(4,2,3,1,f,w,i)) {
                    remove_degree7_zigzag(v,1,f,w,i);
                  } else {
                    remove_degree7_leftfan(v,1,f,w,i);
                  }} else {
                  remove_degree7_leftdelta(v,5,f,w,i);
                }}} else {
            if(! incircle(3,4,5,0,f,w,i)) {
              if(incircle(4,0,1,3,f,w,i)) {
                if(incircle(4,2,3,1,f,w,i)) {
                  remove_degree7_rightfan(v,4,f,w,i);
                } else {
                  remove_degree7_zigzag(v,4,f,w,i);
                }} else {
                remove_degree7_rightfan(v,0,f,w,i);
              }} else {
              remove_degree7_rightdelta(v,3,f,w,i);
            }}} else {
          remove_degree7_star(v,3,f,w,i);
        }} else {
        if(incircle(4,6,0,3,f,w,i)) {
          if(incircle(4,0,1,3,f,w,i)) {
            if(incircle(4,2,3,1,f,w,i)) {
              remove_degree7_star(v,4,f,w,i);
            } else {
              remove_degree7_leftfan(v,4,f,w,i);
            }} else {
            remove_degree7_zigzag(v,0,f,w,i);
          }} else {
          if(incircle(4,5,6,3,f,w,i)) {
            remove_degree7_rightfan(v,3,f,w,i);
          } else {
            remove_degree7_star(v,3,f,w,i);
          }}}}}
}

template < class Gt, class Tds >
inline void
Delaunay_triangulation_2<Gt,Tds>::
rotate7(int j, std::vector<Vertex_handle> &w,
        std::vector<Face_handle> &f, std::vector<int> &i)
{
  if(j==0) return;
  Face_handle ff=f[0];
  int ii=i[0],k=0,kk=(6*j)%7;
  Vertex_handle ww=w[0];
  for (int jj=0; k!=kk; jj=k) { // 7 is prime
    k=(jj+j)%7;
    w[jj]=w[k]; f[jj]=f[k]; i[jj]=i[k];
  }
  w[kk]=ww;f[kk]=ff;i[kk]=ii;
}

template < class Gt, class Tds >
inline void
Delaunay_triangulation_2<Gt,Tds>::
remove_degree7_star   (Vertex_handle &, int j,
std::vector<Face_handle> &f, std::vector<Vertex_handle> &w, std::vector<int> &i)
{ // removing a degree 7 vertex, staring from w[j]
  rotate7(j,w,f,i);

  Face_handle nn;
  f[1]->set_vertex(i[1], w[0]);  // f1 = w1w2w0
  f[2]->set_vertex(i[2], w[0]);  // f2 = w2w3w0
  f[3]->set_vertex(i[3], w[0]);  // f3 = w3w4w0
  f[4]->set_vertex(i[4], w[0]);  // f4 = w4w5w0
  f[5]->set_vertex(i[5], w[0]);  // f5 = w5w6w0

  nn = f[0]->neighbor(i[0]);
  this->tds().set_adjacency(f[1], cw(i[1]), nn, nn->index(f[0]));
  nn = f[6]->neighbor(i[6]);
  this->tds().set_adjacency(f[5], ccw(i[5]), nn, nn->index(f[6]));
  this->tds().delete_face(f[0]);
  this->tds().delete_face(f[6]);
}
template < class Gt, class Tds >
inline void
Delaunay_triangulation_2<Gt,Tds>::
remove_degree7_zigzag (Vertex_handle &, int j,
 std::vector<Face_handle> &f,std::vector<Vertex_handle> &w, std::vector<int> &i)
{ // removing a degree 7 vertex, zigzag, w[j] = middle point

 rotate7(j,w,f,i);

  Face_handle nn;
  f[1]->set_vertex(i[1], w[3]);  // f1 = w1w2w3
  f[2]->set_vertex(ccw(i[2]), w[1]);
  f[2]->set_vertex(i[2], w[0]);  // f2 = w1w3w0
  f[3]->set_vertex(i[3], w[0]);  // f3 = w3w4w0
  f[4]->set_vertex(cw(i[4]), w[6]);
  f[4]->set_vertex(i[4], w[0]);  // f4 = w4w6w0
  f[5]->set_vertex(i[5], w[4]);  // f5 = w5w6w4

  nn = f[2]->neighbor(i[2]);
  this->tds().set_adjacency(f[1], ccw(i[1]), nn, nn->index(f[2]));
  nn = f[0]->neighbor(i[0]);
  this->tds().set_adjacency(f[2], cw(i[2]), nn, nn->index(f[0]));
  nn = f[6]->neighbor(i[6]);
  this->tds().set_adjacency(f[4], ccw(i[4]), nn, nn->index(f[6]));
  nn = f[4]->neighbor(i[4]);
  this->tds().set_adjacency(f[5], cw(i[5]), nn, nn->index(f[4]));
  this->tds().set_adjacency(f[1], cw(i[1]), f[2], i[2]);
  this->tds().set_adjacency(f[4], i[4], f[5], ccw(i[5]));

  this->tds().delete_face(f[0]);
  this->tds().delete_face(f[6]);
}
template < class Gt, class Tds >
inline void
Delaunay_triangulation_2<Gt,Tds>::
remove_degree7_leftdelta(Vertex_handle &, int j,
 std::vector<Face_handle> &f,std::vector<Vertex_handle> &w, std::vector<int> &i)
{ // removing a degree 7 vertex, left delta from w[j]
  rotate7(j,w,f,i);

  Face_handle nn;
  f[1]->set_vertex(i[1], w[0]);  // f1 = w1w2w0
  f[2]->set_vertex(i[2], w[0]);  // f2 = w2w3w0
  f[3]->set_vertex(cw(i[3]), w[5]);
  f[3]->set_vertex(i[3], w[0]);  // f3 = w3w5w0
  f[4]->set_vertex(i[4], w[3]);  // f4 = w4w5w3
  f[5]->set_vertex(i[5], w[0]);  // f5 = w5w6w0

  nn = f[0]->neighbor(i[0]);
  this->tds().set_adjacency(f[1], cw(i[1]), nn, nn->index(f[0]));
  nn = f[3]->neighbor(i[3]);
  this->tds().set_adjacency(f[4], cw(i[4]), nn, nn->index(f[3]));
  nn = f[6]->neighbor(i[6]);
  this->tds().set_adjacency(f[5], ccw(i[5]), nn, nn->index(f[6]));
  this->tds().set_adjacency(f[3], i[3], f[4], ccw(i[4]));
  this->tds().set_adjacency(f[3], ccw(i[3]), f[5], cw(i[5]));

  this->tds().delete_face(f[0]);
  this->tds().delete_face(f[6]);
}

template < class Gt, class Tds >
inline void
Delaunay_triangulation_2<Gt,Tds>::
remove_degree7_rightdelta(Vertex_handle &, int j,
 std::vector<Face_handle> &f,std::vector<Vertex_handle> &w, std::vector<int> &i)
{ // removing a degree 7 vertex, right delta from w[j]
  rotate7(j,w,f,i);

  Face_handle nn;
  f[1]->set_vertex(i[1], w[0]);  // f1 = w1w2w0
  f[2]->set_vertex(i[2], w[4]);  // f2 = w2w3w4
  f[3]->set_vertex(ccw(i[3]), w[2]);
  f[3]->set_vertex(i[3], w[0]);  // f3 = w2w4w0
  f[4]->set_vertex(i[4], w[0]);  // f4 = w4w5w0
  f[5]->set_vertex(i[5], w[0]);  // f5 = w5w6w0

  nn = f[0]->neighbor(i[0]);
  this->tds().set_adjacency(f[1], cw(i[1]), nn, nn->index(f[0]));
  nn = f[3]->neighbor(i[3]);
  this->tds().set_adjacency(f[2], ccw(i[2]), nn, nn->index(f[3]));
  nn = f[6]->neighbor(i[6]);
  this->tds().set_adjacency(f[5], ccw(i[5]), nn, nn->index(f[6]));
  this->tds().set_adjacency(f[1], ccw(i[1]), f[3], cw(i[3]));
  this->tds().set_adjacency(f[3], i[3], f[2], cw(i[2]));

  this->tds().delete_face(f[0]);
  this->tds().delete_face(f[6]);
}
template < class Gt, class Tds >
inline void
Delaunay_triangulation_2<Gt,Tds>::
remove_degree7_leftfan(Vertex_handle &, int j,
 std::vector<Face_handle> &f,std::vector<Vertex_handle> &w, std::vector<int> &i)
{ // removing a degree 7 vertex, left fan from w[j]
  rotate7(j,w,f,i);

  Face_handle nn;
  f[1]->set_vertex(i[1], w[0]);  // f1 = w1w2w0
  f[2]->set_vertex(i[2], w[0]);  // f2 = w2w3w0
  f[3]->set_vertex(i[3], w[0]);  // f3 = w3w4w0
  f[4]->set_vertex(i[4], w[6]);  // f4 = w4w5w6
  f[6]->set_vertex(i[6], w[4]);  // f6 = w6w0w4

  nn = f[0]->neighbor(i[0]);
  this->tds().set_adjacency(f[1], cw(i[1]), nn, nn->index(f[0]));
  nn = f[5]->neighbor(i[5]);
  this->tds().set_adjacency(f[4], ccw(i[4]), nn, nn->index(f[5]));
  this->tds().set_adjacency(f[3], ccw(i[3]), f[6], ccw(i[6]));
  this->tds().set_adjacency(f[6], cw(i[6]), f[4], cw(i[4]));

  this->tds().delete_face(f[0]);
  this->tds().delete_face(f[5]);
}

template < class Gt, class Tds >
inline void
Delaunay_triangulation_2<Gt,Tds>::
remove_degree7_rightfan(Vertex_handle &, int j,
 std::vector<Face_handle> &f,std::vector<Vertex_handle> &w, std::vector<int> &i)
{ // removing a degree 7 vertex, right fan from w[j]
  rotate7(j,w,f,i);

  Face_handle nn;
  f[0]->set_vertex(i[0], w[3]);  // f0 = w0w1w3
  f[2]->set_vertex(i[2], w[1]);  // f2 = w2w3w1
  f[3]->set_vertex(i[3], w[0]);  // f3 = w3w4w0
  f[4]->set_vertex(i[4], w[0]);  // f4 = w4w5w0
  f[5]->set_vertex(i[5], w[0]);  // f5 = w5w6w0

  nn = f[1]->neighbor(i[1]);
  this->tds().set_adjacency(f[2], cw(i[2]), nn, nn->index(f[1]));
  nn = f[6]->neighbor(i[6]);
  this->tds().set_adjacency(f[5], ccw(i[5]), nn, nn->index(f[6]));
  this->tds().set_adjacency(f[2], ccw(i[2]), f[0], ccw(i[0]));
  this->tds().set_adjacency(f[0], cw(i[0]), f[3], cw(i[3]));

  this->tds().delete_face(f[1]);
  this->tds().delete_face(f[6]);
}

///////////////////////////////////////////////////////////////
//  DISPLACEMENT

template <class Gt, class Tds >
typename Delaunay_triangulation_2<Gt,Tds>::Vertex_handle
Delaunay_triangulation_2<Gt,Tds>::
move_if_no_collision(Vertex_handle v, const Point &p)
{
  CGAL_triangulation_precondition(!this->is_infinite(v));
  if(v->point() == p) return v;
  const int dim = this->dimension();

  if(dim == 2) {
    Point ant = v->point();
    v->set_point(p);
    // This option optimizes only when most of the
    // displacements would not break the orientation
    // of the faces.. we will consider this as an a priori,
    // because otherwise it is pointless to just do
    // not rebuild from scratch.
    if(this->well_oriented(v)) {
      restore_edges(v);
      return v;
    }
    v->set_point(ant);
  }

  Locate_type lt;
  int li;
  Vertex_handle inserted;
  Face_handle loc = this->locate(p, lt, li, v->face());

  if(lt == Triangulation_2<Gt,Tds>::VERTEX)
    return loc->vertex(li);

  if(dim == 0) {
    v->point() = p;
    return v;
  }

  size_type n_vertices = this->tds().number_of_vertices();

  if((lt == Triangulation::OUTSIDE_AFFINE_HULL) &&
     (dim == 1) && (n_vertices == 3)) {
    v->point() = p;
    return v;
  }

  if((lt != Triangulation::OUTSIDE_AFFINE_HULL) && (dim == 1)) {
    if(loc->has_vertex(v)) {
      v->point() = p;
    } else {
      inserted = insert(p, lt, loc, li);
      Face_handle f = v->face();
      int i = f->index(v);
      if(i==0) {f = f->neighbor(1);}
      CGAL_triangulation_assertion(f->index(v) == 1);
      Face_handle g= f->neighbor(0);
      f->set_vertex(1, g->vertex(1));
      f->set_neighbor(0,g->neighbor(0));
      g->neighbor(0)->set_neighbor(1,f);
      g->vertex(1)->set_face(f);
      this->delete_face(g);
      Face_handle f_ins = inserted->face();
      i = f_ins->index(inserted);
      if(i==0) {f_ins = f_ins->neighbor(1);}
      CGAL_triangulation_assertion(f_ins->index(inserted) == 1);
      Face_handle g_ins = f_ins->neighbor(0);
      f_ins->set_vertex(1, v);
      g_ins->set_vertex(0, v);
      v->set_point(p);
      v->set_face(inserted->face());
      this->delete_vertex(inserted);
    }
    return v;
  }

  if((lt != Triangulation::OUTSIDE_AFFINE_HULL) && this->test_dim_down(v)) {
    // verify if p and two static vertices are collinear in this case
    int iinf = 0;
    Face_circulator finf = this->incident_faces(this->infinite_vertex()),
                    fdone(finf);
    do {
      if(!finf->has_vertex(v))
      {
        iinf = ~(finf->index(this->infinite_vertex()));
        break;
      }
    } while(++finf != fdone);
    if(this->orientation(finf->vertex(iinf&1)->point(),
                         finf->vertex(iinf&2)->point(),
                         p) == COLLINEAR)
    {
      v->point() = p;
      this->tds().dim_down(loc, loc->index(v));
      return v;
    }
  }

  inserted = insert(p, lt, loc, li);

  {
    int d;
    CGAL_STATIC_THREAD_LOCAL_VARIABLE(int, maxd,30);
    CGAL_STATIC_THREAD_LOCAL_VARIABLE(std::vector<Face_handle>, f, maxd);
    CGAL_STATIC_THREAD_LOCAL_VARIABLE(std::vector<int>, i, maxd);
    CGAL_STATIC_THREAD_LOCAL_VARIABLE(std::vector<Vertex_handle>, w, maxd);

    remove_degree_init(v,f,w,i,d,maxd);
    remove_degree_triangulate(v,f,w,i,d);
  }

  // fixing pointer
  Face_circulator fc = this->incident_faces(inserted), done(fc);
  std::vector<Face_handle> faces_pt;
  faces_pt.reserve(16);
  do { faces_pt.push_back(fc); } while(++fc != done);
  std::size_t ss = faces_pt.size();
  for(std::size_t k=0; k<ss; k++)
  {
    Face_handle f = faces_pt[k];
    int i = f->index(inserted);
    f->set_vertex(i, v);
  }

  v->set_point(p);
  v->set_face(inserted->face());

  this->delete_vertex(inserted);

  return v;
}

template <class Gt, class Tds >
typename Delaunay_triangulation_2<Gt,Tds>::Vertex_handle
Delaunay_triangulation_2<Gt,Tds>::
move(Vertex_handle v, const Point &p)
{
  CGAL_triangulation_precondition(!this->is_infinite(v));
  if(v->point() == p) return v;
  Vertex_handle w = move_if_no_collision(v,p);
  if(w != v) {
    remove(v);
    return w;
  }
  return v;
}

template <class Gt, class Tds >
bool
Delaunay_triangulation_2<Gt,Tds>::
is_delaunay_after_displacement(Vertex_handle v, const Point &p) const
{
  CGAL_triangulation_precondition(!this->is_infinite(v));
  CGAL_triangulation_precondition(this->dimension() == 2);
  CGAL_triangulation_precondition(!this->test_dim_down(v));
  if(v->point() == p) return true;
  Point ant = v->point();
  v->set_point(p);
  if(!this->well_oriented(v))
  {
    v->set_point(ant);
    return false;
  }
  std::list<Edge> edges;
  Face_circulator fc = this->incident_faces(v), done(fc);
  int degree = 0;
  do {
    if((++degree) > 3) break;
  } while(++fc != done);
  fc = this->incident_faces(v);
  done = fc;
  if(degree == 3) {
    do {
      int i = fc->index(v);
      edges.push_back(Edge(fc, i));
    } while(++fc != done);
  } else {
    do {
      int i = fc->index(v);
      edges.push_back(Edge(fc, i));
      edges.push_back(Edge(fc, this->cw(i)));
    } while(++fc != done);
  }
  while(!edges.empty()) {
    const Edge &e = edges.front();
    Face_handle f = e.first;
    int i = e.second;
    edges.pop_front();
    if(this->is_infinite(f->vertex(i))) continue;
    Vertex_handle vm = this->_tds.mirror_vertex(f, i);
    if(this->is_infinite(vm)) continue;
    if(this->side_of_oriented_circle(f, vm->point()) == ON_POSITIVE_SIDE) {
      v->set_point(ant);
      return false;
    }
  }
  v->set_point(ant);
  return true;
}

template <class Gt, class Tds >
template <class OutputItFaces>
typename Delaunay_triangulation_2<Gt,Tds>::Vertex_handle
Delaunay_triangulation_2<Gt,Tds>::
move_if_no_collision_and_give_new_faces(Vertex_handle v,
                                        const Point &p,
                                        OutputItFaces oif)
{
  CGAL_triangulation_precondition(!this->is_infinite(v));
  if(v->point() == p) return v;

  const int dim = this->dimension();

  if(dim == 2) {
    Point ant = v->point();
    v->set_point(p);
    // This option optimizes only when most of the
    // displacements would not break the orientation
    // of the faces.. we will consider this as an a priori,
    // because otherwise it is pointless to just do
    // not rebuild from scratch.
    if(well_oriented(v)) {
      std::set<Face_handle> faces_set;
      restore_edges(v, faces_set);
      for(typename std::set<Face_handle>::iterator ib = faces_set.begin(),
            iend = faces_set.end(); ib != iend; ib++) *oif++ = *ib;
      return v;
    }
    v->set_point(ant);
  }

  Locate_type lt;
  int li;
  Vertex_handle inserted;
  Face_handle loc = this->locate(p, lt, li, v->face());

  if(lt == Triangulation::VERTEX) return loc->vertex(li);

  if(dim == 0) {
    v->point() = p;
    return v;
  }

  size_type n_vertices = this->tds().number_of_vertices();

  if((lt == Triangulation::OUTSIDE_AFFINE_HULL) &&
     (dim == 1) && (n_vertices == 3)) {
    v->point() = p;
    for(All_faces_iterator afi = this-> all_faces_begin();
        afi != this->all_faces_end();
        afi++) *oif++ = afi;
    return v;
  }

  if((lt != Triangulation::OUTSIDE_AFFINE_HULL) && (dim == 1)) {
    if(loc->has_vertex(v)) {
      v->point() = p;
    } else {
      inserted = insert(p, lt, loc, li);
      Face_handle f = v->face();
      int i = f->index(v);
      if(i==0) {f = f->neighbor(1);}
      CGAL_triangulation_assertion(f->index(v) == 1);
      Face_handle g= f->neighbor(0);
      f->set_vertex(1, g->vertex(1));
      f->set_neighbor(0,g->neighbor(0));
      g->neighbor(0)->set_neighbor(1,f);
      g->vertex(1)->set_face(f);
      this->delete_face(g);
      *oif++ = f;
      Face_handle f_ins = inserted->face();
      i = f_ins->index(inserted);
      if(i==0) {f_ins = f_ins->neighbor(1);}
      CGAL_triangulation_assertion(f_ins->index(inserted) == 1);
      Face_handle g_ins = f_ins->neighbor(0);
      f_ins->set_vertex(1, v);
      g_ins->set_vertex(0, v);
      v->set_face(inserted->face());
      v->set_point(p);
      this->delete_vertex(inserted);
    }
    *oif++ = v->face();
    if(v->face()->neighbor(0)->has_vertex(v))
      *oif++ = v->face()->neighbor(0);
    if(v->face()->neighbor(1)->has_vertex(v))
      *oif++ = v->face()->neighbor(1);
    return v;
  }

  if((lt != Triangulation::OUTSIDE_AFFINE_HULL) && this->test_dim_down(v)) {
    // verify if p and two static vertices are collinear in this case
    int iinf;
    Face_circulator finf = incident_faces(this->infinite_vertex()),
      fdone(finf);
    do {
      if(!finf->has_vertex(v))
      {
        iinf = ~(finf->index(this->infinite_vertex()));
        break;
      }
    } while(++finf != fdone);
    if(this->orientation(finf->vertex(iinf&1)->point(),
                         finf->vertex(iinf&2)->point(),
                         p) == COLLINEAR)
    {
      v->point() = p;
      this->tds().dim_down(loc, loc->index(v));
      for(All_faces_iterator afi = this-> all_faces_begin();
          afi != this->all_faces_end();
          afi++) *oif++ = afi;
      return v;
    }
  }

  std::set<Face_handle> faces_set;
  inserted = Delaunay_triangulation_2<Gt,Tds>::insert(p, lt, loc, li);
  Face_circulator fc = this->incident_faces(inserted), done(fc);
  do { faces_set.insert(fc); } while(++fc != done);

  {
    CGAL_STATIC_THREAD_LOCAL_VARIABLE(int, maxd,30);
    CGAL_STATIC_THREAD_LOCAL_VARIABLE(std::vector<Face_handle>, f, maxd);
    CGAL_STATIC_THREAD_LOCAL_VARIABLE(std::vector<int>, i, maxd);
    CGAL_STATIC_THREAD_LOCAL_VARIABLE(std::vector<Vertex_handle>, w, maxd);

    int d;
    remove_degree_init(v,f,w,i,d,maxd);
    remove_degree_triangulate(v,f,w,i,d);
    this->delete_vertex(v);
    Face_circulator fc(v[0]),done;
    do *oif++ = fc++; while (fc!=done);
  }

  fc = this->incident_faces(inserted), done(fc);
  std::vector<Face_handle> faces_pt;
  faces_pt.reserve(16);
  do { faces_pt.push_back(fc); } while(++fc != done);
  int ss = faces_pt.size();
  for(int k=0; k<ss; k++)
  {
    Face_handle f = faces_pt[k];
    int i = f->index(inserted);
    f->set_vertex(i, v);
  }
  v->set_point(p);
  v->set_face(inserted->face());
  this->delete_vertex(inserted);

  for(typename std::set<Face_handle>::const_iterator ib = faces_set.begin(),
        iend = faces_set.end(); ib != iend; ib++)
    *oif++ = *ib;

  return v;
}

} //namespace CGAL

#endif // CGAL_DELAUNAY_TRIANGULATION_2_H
