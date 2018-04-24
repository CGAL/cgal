// Copyright (c) 1997-2013 INRIA Sophia-Antipolis (France).
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
// Author(s)     : Nico Kruithof <Nico@nghk.nl>

#ifndef CGAL_PERIODIC_2_DELAUNAY_TRIANGULATION_2_H
#define CGAL_PERIODIC_2_DELAUNAY_TRIANGULATION_2_H

#include <CGAL/license/Periodic_2_triangulation_2.h>

#include <CGAL/Periodic_2_triangulation_2.h>
#include <CGAL/iterator.h>
#include <CGAL/algorithm.h>

#ifndef CGAL_TRIANGULATION_2_DONT_INSERT_RANGE_OF_POINTS_WITH_INFO
#include <CGAL/Spatial_sort_traits_adapter_2.h>
#include <CGAL/internal/info_check.h>
#include <CGAL/tss.h>

#include <boost/iterator/zip_iterator.hpp>
#include <boost/mpl/and.hpp>
#endif //CGAL_TRIANGULATION_2_DONT_INSERT_RANGE_OF_POINTS_WITH_INFO

namespace CGAL
{

template <
  class Gt,
  class Tds = Triangulation_data_structure_2 <
    Periodic_2_triangulation_vertex_base_2<Gt>,
    Periodic_2_triangulation_face_base_2<Gt> > >
class Periodic_2_Delaunay_triangulation_2
  : public Periodic_2_triangulation_2<Gt, Tds>
{
  typedef Periodic_2_Delaunay_triangulation_2<Gt, Tds>          Self;
public:
  typedef Periodic_2_triangulation_2<Gt, Tds>                   Base;

public:
  typedef Tds                                  Triangulation_data_structure;
  typedef Gt                                   Geom_traits;

  typedef typename Gt::Periodic_2_offset_2     Offset;
  typedef typename Gt::Iso_rectangle_2         Iso_rectangle;
  typedef array<int, 2>                        Covering_sheets;

  typedef typename Gt::FT                      FT;
  typedef typename Gt::Point_2                 Point;
  typedef typename Gt::Segment_2               Segment;
  typedef typename Gt::Triangle_2              Triangle;

  typedef std::pair<Point, Offset>              Periodic_point;
  typedef array< std::pair<Point, Offset>, 2>   Periodic_segment;
  typedef array< std::pair<Point, Offset>, 3>   Periodic_triangle;
  typedef array< std::pair<Point, Offset>, 4>   Periodic_tetrahedron;

  typedef typename Base::size_type              size_type;
  typedef typename Base::Locate_type            Locate_type;
  typedef typename Base::Face_handle            Face_handle;
  typedef typename Base::Vertex_handle          Vertex_handle;
  typedef typename Base::Edge                   Edge;
  typedef typename Base::Edge_circulator        Edge_circulator;
  typedef typename Base::Face_circulator        Face_circulator;
  typedef typename Base::Vertex_circulator      Vertex_circulator;
  typedef typename Base::Finite_edges_iterator  Finite_edges_iterator;
  typedef typename Base::Finite_faces_iterator  Finite_faces_iterator;
  typedef typename Base::Finite_vertices_iterator Finite_vertices_iterator;
  typedef typename Base::All_faces_iterator     All_faces_iterator;

  typedef typename Base::Edge_iterator          Edge_iterator;
  typedef typename Base::Face_iterator          Face_iterator;
  typedef typename Base::Vertex_iterator        Vertex_iterator;

  typedef typename Base::Periodic_segment_iterator  Periodic_segment_iterator;
  typedef typename Base::Periodic_triangle_iterator Periodic_triangle_iterator;

  //Tag to distinguish Delaunay from regular triangulations
  typedef Tag_false                             Weighted_tag;

  // Tag to distinguish periodic triangulations from others
  typedef Tag_true                              Periodic_tag;

public:
#ifndef CGAL_CFG_USING_BASE_MEMBER_BUG_2
  using Base::empty;
  using Base::cw;
  using Base::ccw;
  using Base::create_face;
  using Base::insert_too_long_edge;
  using Base::locate;
  using Base::remove_degree_init;

  using Base::combine_offsets;
  using Base::get_offset;
  using Base::get_neighbor_offset;
  using Base::int_to_off;
  using Base::set_offsets;
  using Base::is_1_cover;
  using Base::number_of_sheets;

  using Base::dimension;
  using Base::domain;
  using Base::geom_traits;
  using Base::tds;
  using Base::is_infinite;
  using Base::number_of_vertices;
  using Base::faces_begin;
  using Base::finite_edges_begin;
  using Base::finite_edges_end;
  using Base::incident_faces;

  using Base::orientation;
  using Base::point;
  using Base::construct_segment;
#endif

  /// \name Constructors
  // \{
  /// Constructor
  Periodic_2_Delaunay_triangulation_2(const Iso_rectangle & domain = Iso_rectangle(0, 0, 1, 1),
                                      const Gt& gt = Gt())
    : Base(domain, gt) {}

  /// Copy constructor
  Periodic_2_Delaunay_triangulation_2(const Periodic_2_Delaunay_triangulation_2<Gt, Tds> &tr)
    : Base(tr)
  {
    CGAL_triangulation_postcondition( is_valid(true) );
  }

  /// Constructor with insertion of points
  template < class InputIterator >
  Periodic_2_Delaunay_triangulation_2(InputIterator first, InputIterator last,
                                      const Iso_rectangle & domain = Iso_rectangle(0, 0, 1, 1),
                                      const Gt& gt = Gt())
    : Periodic_2_triangulation_2<Gt, Tds>(domain, gt)
  {
    insert(first, last);
  }

  // \}

  /// \name Methods regarding the covering
  /// \{

  /// Checks whether the triangulation is a valid simplicial complex in the one cover.
  /// Uses an edge-length-criterion.
  bool is_extensible_triangulation_in_1_sheet_h1() const
  {
    if(!is_1_cover())
      return (this->_too_long_edge_counter == 0);

    FT longest_edge_squared_length(0);
    Segment s;

    for(Periodic_segment_iterator psit = this->periodic_segments_begin(Base::UNIQUE);
                                  psit != this->periodic_segments_end(Base::UNIQUE); ++psit)
    {
      s = construct_segment(*psit);
      longest_edge_squared_length = (std::max)(longest_edge_squared_length,
                                               s.squared_length());
    }
    return (longest_edge_squared_length < this->_edge_length_threshold);
  }

  /// Checks whether the triangulation is a valid simplicial complex in the one cover.
  /// Uses a criterion based on the maximal radius of the circumscribing circle.
  bool is_extensible_triangulation_in_1_sheet_h2() const
  {
    for(Periodic_triangle_iterator tit = this->periodic_triangles_begin(Base::UNIQUE);
                                   tit != this->periodic_triangles_end(Base::UNIQUE); ++tit)
    {
      Point cc = geom_traits().construct_circumcenter_2_object()(
                   tit->at(0).first, tit->at(1).first, tit->at(2).first,
                   tit->at(0).second, tit->at(1).second, tit->at(2).second);

      if (!(FT(16) * squared_distance(cc, point(tit->at(0))) <
            (domain().xmax() - domain().xmin()) * (domain().xmax() - domain().xmin())))
        return false;
    }
    return true;
  }

  // \}

  /// \name Insertion-Removal
  // \{
  Vertex_handle insert(const Point  &p,
                       Face_handle start = Face_handle() );
  Vertex_handle insert(const Point& p,
                       Locate_type lt,
                       Face_handle loc, int li );


  /// Inserts a point in the triangulation.
  Vertex_handle push_back(const Point &p);

#ifndef CGAL_TRIANGULATION_2_DONT_INSERT_RANGE_OF_POINTS_WITH_INFO
  template < class InputIterator >
  std::ptrdiff_t
  insert(InputIterator first, InputIterator last,
         bool is_large_point_set = true,
         typename boost::enable_if <
         boost::is_convertible <
         typename std::iterator_traits<InputIterator>::value_type,
         Point
         > >::type* = NULL)
#else
  template < class InputIterator >
  std::ptrdiff_t
  insert(InputIterator first, InputIterator last, bool is_large_point_set = true)
#endif //CGAL_TRIANGULATION_2_DONT_INSERT_RANGE_OF_POINTS_WITH_INFO 
  {
    if (first == last) return 0;

    size_type n = number_of_vertices();

    // The heuristic discards the existing triangulation so it can only be
    // applied to empty triangulations.
    if (n != 0) is_large_point_set = false;

    std::set<Vertex_handle> dummy_points;
    std::vector<Point> points(first, last);
    typename std::vector<Point>::iterator pbegin = points.begin();

    if (is_large_point_set)
      {
        std::vector<Vertex_handle> tmp_dummy_points = this->insert_dummy_points();
        std::copy(tmp_dummy_points.begin(), tmp_dummy_points.end(),
                  std::inserter(dummy_points, dummy_points.begin()));
      }
    else
      {
        CGAL::cpp98::random_shuffle (points.begin(), points.end());
        pbegin = points.begin();

        // The empty triangulation is a 1-cover by definition, insert at least one point
        insert(*pbegin);
        ++pbegin;
        while (!is_1_cover())
          {
            if (pbegin == points.end())
              return number_of_vertices() - n;
            insert(*pbegin);
            ++pbegin;
          }
      }

    CGAL_assertion(is_1_cover());

    // Insert the points
    spatial_sort (pbegin, points.end(), geom_traits());

    Face_handle f;
    Locate_type lt;
    int li;

    for (typename std::vector<Point>::const_iterator p = pbegin, end = points.end();
         p != end; ++p)
      {
        f = locate(*p, lt, li, f);

        if (lt == Base::VERTEX)
          {
            dummy_points.erase(f->vertex(li));
          }
        else
          {
            insert(*p, lt, f, li);
          }
      }

    if (is_large_point_set)
      {
        for (typename std::set<Vertex_handle>::const_iterator it = dummy_points.begin();
             it != dummy_points.end(); ++it)
          {
            remove(*it);
          }
      }

    return number_of_vertices() - n;
  }
#ifndef CGAL_TRIANGULATION_2_DONT_INSERT_RANGE_OF_POINTS_WITH_INFO
private:
  //top stands for tuple-or-pair
  template <class Info>
  const Point& top_get_first(const std::pair<Point, Info>& pair) const
  {
    return pair.first;
  }
  template <class Info>
  const Info& top_get_second(const std::pair<Point, Info>& pair) const
  {
    return pair.second;
  }
  template <class Info>
  const Point& top_get_first(const boost::tuple<Point, Info>& tuple) const
  {
    return boost::get<0>(tuple);
  }
  template <class Info>
  const Info& top_get_second(const boost::tuple<Point, Info>& tuple) const
  {
    return boost::get<1>(tuple);
  }

  template <class Tuple_or_pair, class InputIterator>
  std::ptrdiff_t insert_with_info(InputIterator first, InputIterator last, bool is_large_point_set)
  {
    if (first == last) return 0;

    std::vector<std::size_t> indices;
    std::vector<Point> points;
    std::vector<typename Tds::Vertex::Info> infos;
    std::size_t index = 0;
    for (InputIterator it = first; it != last; ++it)
      {
        Tuple_or_pair value = *it;
        points.push_back( top_get_first(value)  );
        infos.push_back ( top_get_second(value) );
        indices.push_back(index++);
      }

    typedef typename Pointer_property_map<Point>::type Pmap;
    typedef Spatial_sort_traits_adapter_2<Geom_traits, Pmap> Search_traits;

    size_type n = number_of_vertices();

    // The heuristic discards the existing triangulation so it can only be
    // applied to empty triangulations.
    if (n != 0) is_large_point_set = false;

    std::set<Vertex_handle> dummy_points;
    typename std::vector<std::size_t>::iterator pbegin = indices.begin();

    if (is_large_point_set)
      {
        std::vector<Vertex_handle> tmp_dummy_points = this->insert_dummy_points();
        std::copy(tmp_dummy_points.begin(), tmp_dummy_points.end(),
                  std::inserter(dummy_points, dummy_points.begin()));
      }
    else
      {
        CGAL::cpp98::random_shuffle(indices.begin(), indices.end());
        pbegin = indices.begin();

        Vertex_handle v_new;

        // The empty triangulation is a 1-cover by definition, insert at least one point
        v_new = insert(points[*pbegin]);
        v_new->info() = infos[*pbegin];
        ++pbegin;

        while (!is_1_cover())
          {
            if (pbegin == indices.end())
              return number_of_vertices() - n;
            v_new = insert(points[*pbegin]);
            v_new->info() = infos[*pbegin];
            ++pbegin;
          }
      }

    CGAL_assertion(is_1_cover());

    // Insert the points
    spatial_sort(indices.begin(),
                 indices.end(),
                 Search_traits(make_property_map(points), geom_traits()));

    Face_handle f;
    Locate_type lt;
    int li;

    Face_handle hint;
    for (typename std::vector<std::size_t>::const_iterator it = pbegin, end = indices.end();
         it != end; ++it)
      {
        f = locate(points[*it], lt, li, f);

        if (lt == Base::VERTEX)
          {
            // Always copy the info, it might be a dummy vertex
            f->vertex(li)->info() = infos[*it];
            dummy_points.erase(f->vertex(li));
          }
        else
          {
            Vertex_handle v_new = insert(points[*it], lt, f, li);
            v_new->info() = infos[*it];
          }
      }

    if (is_large_point_set)
      {
        for (typename std::set<Vertex_handle>::const_iterator it = dummy_points.begin();
             it != dummy_points.end(); ++it)
          {
            remove(*it);
          }
      }

    return number_of_vertices() - n;
  }


public:

  template < class InputIterator >
  std::ptrdiff_t
  insert( InputIterator first,
          InputIterator last,
          bool is_large_point_set = true,
          typename boost::enable_if <
          boost::is_convertible <
          typename std::iterator_traits<InputIterator>::value_type,
          std::pair<Point, typename internal::Info_check<typename Tds::Vertex>::type>
          > >::type* = NULL
        )
  {
    return insert_with_info< std::pair<Point, typename internal::Info_check<typename Tds::Vertex>::type> >(first, last, is_large_point_set);
  }

  template <class  InputIterator_1, class InputIterator_2>
  std::ptrdiff_t
  insert( boost::zip_iterator< boost::tuple<InputIterator_1, InputIterator_2> > first,
          boost::zip_iterator< boost::tuple<InputIterator_1, InputIterator_2> > last,
          bool is_large_point_set = true,
          typename boost::enable_if <
          boost::mpl::and_ <
          boost::is_convertible< typename std::iterator_traits<InputIterator_1>::value_type, Point >,
          boost::is_convertible< typename std::iterator_traits<InputIterator_2>::value_type, typename internal::Info_check<typename Tds::Vertex>::type >
          > >::type* = NULL)
  {
    return insert_with_info< boost::tuple<Point, typename internal::Info_check<typename Tds::Vertex>::type> >(first, last, is_large_point_set);
  }
#endif //CGAL_TRIANGULATION_2_DONT_INSERT_RANGE_OF_POINTS_WITH_INFO

  void  remove(Vertex_handle v );
  // \}

  /// \name Displacement
  // \{

  Vertex_handle move_if_no_collision(Vertex_handle v, const Point &p);
  Vertex_handle move_point(Vertex_handle v, const Point &p);
  // \}

  /// \name Check - Query
  // \{
  /// Returns the vertex closest to p, the point location will start from f
  Vertex_handle
  nearest_vertex(const Point& p, Face_handle f = Face_handle()) const;

  template <class OutputItFaces, class OutputItBoundaryEdges>
  std::pair<OutputItFaces, OutputItBoundaryEdges>
  get_conflicts_and_boundary(const Point  &p,
                             OutputItFaces fit,
                             OutputItBoundaryEdges eit,
                             Face_handle start = Face_handle()) const
  {
    CGAL_triangulation_precondition( dimension() == 2);
    int li;
    Locate_type lt;
    Face_handle fh = locate(p, lt, li, start);
    switch(lt)
    {
      default:
        break;
      case Base::VERTEX:
        return std::make_pair(fit, eit);
      case Base::FACE:
      case Base::EDGE:
      case Base::EMPTY:
        *fit++ = fh; //put fh in OutputItFaces
        std::pair<OutputItFaces, OutputItBoundaryEdges> pit = std::make_pair(fit, eit);
        pit = propagate_conflicts(p, fh, 0, pit);
        pit = propagate_conflicts(p, fh, 1, pit);
        pit = propagate_conflicts(p, fh, 2, pit);
        return pit;
    }
    CGAL_triangulation_assertion(false);
    return std::make_pair(fit, eit);
  }

  template <class OutputItFaces>
  OutputItFaces
  get_conflicts (const Point  &p,
                 OutputItFaces fit,
                 Face_handle start = Face_handle()) const
  {
    std::pair<OutputItFaces, Emptyset_iterator> pp =
      get_conflicts_and_boundary(p, fit, Emptyset_iterator(), start);
    return pp.first;
  }

  template <class OutputItBoundaryEdges>
  OutputItBoundaryEdges
  get_boundary_of_conflicts(const Point  &p,
                            OutputItBoundaryEdges eit,
                            Face_handle start = Face_handle()) const
  {
    std::pair<Emptyset_iterator, OutputItBoundaryEdges> pp =
      get_conflicts_and_boundary(p, Emptyset_iterator(), eit, start);
    return pp.second;
  }
  // \}

  /// Constructs the circumcenter of the face f, respects the offset
  Point circumcenter(Face_handle f) const
  {
    return construct_circumcenter(f->vertex(0)->point(),
                                  f->vertex(1)->point(),
                                  f->vertex(2)->point(),
                                  get_offset(f, 0),
                                  get_offset(f, 1),
                                  get_offset(f, 2));
  }
  Point construct_circumcenter(const Point &p1, const Point &p2, const Point &p3,
                               const Offset &o1, const Offset &o2, const Offset &o3) const
  {
    return geom_traits().construct_circumcenter_2_object()(p1, p2, p3, o1, o2, o3);
  }

  /// \name Dual
  // \{
  /// Returns the dual of f, which is the circumcenter of f.
  Point dual(Face_handle f) const;
  /// Returns the dual of e, which is always a segment in the periodic triangulation.
  Segment dual(const Edge &e) const ;
  /// Returns the dual of the edge pointed to by ec.
  Segment dual(const Edge_circulator& ec) const;
  /// Returns the dual of the edge pointed to by ei.
  Segment dual(const Edge_iterator& ei) const;

  template < class Stream>
  Stream& draw_dual(Stream & ps)
  {
    Finite_edges_iterator eit = finite_edges_begin();
    for (; eit != finite_edges_end(); ++eit)
      {
        ps << dual(eit);
      }
    return ps;
  }
  // \}

  /// \name Checking
  // \{
  bool is_valid(bool verbose = false, int level = 0) const;

  /// Checks whether f->vertex(i) lies outside the circumcircle of the face nb
  inline bool locally_Delaunay(const Face_handle &f, int i, const Face_handle &nb);
  // \}

private:
  /// Not in the documentation
  inline void restore_Delaunay(Vertex_handle v);

#ifndef CGAL_DT2_USE_RECURSIVE_PROPAGATING_FLIP
  void non_recursive_propagating_flip(Face_handle f, int i);
  void propagating_flip(const Face_handle& f, int i, int depth = 0);
#else
  void propagating_flip(const Face_handle& f, int i);
#endif

  // auxilliary functions for remove
  // returns false if we first need to convert to a 9-cover before the vertex can be removed
  bool remove_single_vertex(Vertex_handle v, const Offset &v_o);
  void remove_degree_triangulate(Vertex_handle v, std::vector<Face_handle> &f,
                                 std::vector<Vertex_handle> &w,
                                 std::vector<int> &i, int d);
  void remove_degree_triangulate(Vertex_handle v, std::vector<Face_handle> &f,
                                 std::vector<Vertex_handle> &w, std::vector<Offset> &offset_w,
                                 std::vector<int> &i, int d);
  void remove_degree_d(Vertex_handle v, std::vector<Face_handle> &f,
                       std::vector<Vertex_handle> &w,
                       std::vector<int> &i, int d);
  void remove_degree_d(Vertex_handle v, std::vector<Face_handle> &f,
                       std::vector<Vertex_handle> &w, std::vector<Offset> &offset_w,
                       std::vector<int> &i, int d);
  /// Assumes that all offsets are (0,0)
  void fill_hole_delaunay(std::list<Edge> & hole);
  /// Fill hole over a periodic boundary
  void fill_hole_delaunay(std::list<Edge> & hole,
                          std::map<Vertex_handle, Offset> &vertex_offsets);

  void remove_degree3(Vertex_handle v, std::vector<Face_handle> &f,
                      std::vector<Vertex_handle> &w,
                      std::vector<int> &i);
  void remove_degree3(Vertex_handle v, std::vector<Face_handle> &f,
                      std::vector<Vertex_handle> &w, std::vector<Offset> &o,
                      std::vector<int> &i);
  void remove_degree4(Vertex_handle v, std::vector<Face_handle> &f,
                      std::vector<Vertex_handle> &w, std::vector<Offset> &o,
                      std::vector<int> &i);
  void remove_degree4(Vertex_handle v, std::vector<Face_handle> &f,
                      std::vector<Vertex_handle> &w,
                      std::vector<int> &i);
  void remove_degree5(Vertex_handle v, std::vector<Face_handle> &f,
                      std::vector<Vertex_handle> &w,
                      std::vector<int> &i);
  void remove_degree5(Vertex_handle v, std::vector<Face_handle> &f,
                      std::vector<Vertex_handle> &w, std::vector<Offset> &o,
                      std::vector<int> &i);
  void remove_degree5_star   (Vertex_handle &v,
                              Face_handle & , Face_handle & , Face_handle & , Face_handle & , Face_handle & ,
                              Vertex_handle&, Vertex_handle&, Vertex_handle&, Vertex_handle&, Vertex_handle&,
                              int           , int           , int           , int           , int );
  void remove_degree5_star   (Vertex_handle &v,
                              Face_handle & , Face_handle & , Face_handle & , Face_handle & , Face_handle & ,
                              Vertex_handle&, Vertex_handle&, Vertex_handle&, Vertex_handle&, Vertex_handle&,
                              Offset&, Offset&, Offset&, Offset&, Offset&,
                              int           , int           , int           , int           , int );
  void remove_degree6(Vertex_handle v, std::vector<Face_handle> &f,
                      std::vector<Vertex_handle> &w,
                      std::vector<int> &i);
  void remove_degree6(Vertex_handle v, std::vector<Face_handle> &f,
                      std::vector<Vertex_handle> &w, std::vector<Offset> &o,
                      std::vector<int> &i);
  void remove_degree6_star   (Vertex_handle &v,
                              Face_handle & , Face_handle & , Face_handle & ,
                              Face_handle & , Face_handle & , Face_handle & ,
                              Vertex_handle&, Vertex_handle&, Vertex_handle&,
                              Vertex_handle&, Vertex_handle&, Vertex_handle&,
                              int           , int           , int           ,
                              int           , int           , int );
  void remove_degree6_star   (Vertex_handle &v,
                              Face_handle & , Face_handle & , Face_handle & ,
                              Face_handle & , Face_handle & , Face_handle & ,
                              Vertex_handle&, Vertex_handle&, Vertex_handle&,
                              Vertex_handle&, Vertex_handle&, Vertex_handle&,
                              Offset&, Offset&, Offset&,
                              Offset&, Offset&, Offset&,
                              int           , int           , int           ,
                              int           , int           , int );
  void remove_degree6_N      (Vertex_handle &v,
                              Face_handle & , Face_handle & , Face_handle & ,
                              Face_handle & , Face_handle & , Face_handle & ,
                              Vertex_handle&, Vertex_handle&, Vertex_handle&,
                              Vertex_handle&, Vertex_handle&, Vertex_handle&,
                              int           , int           , int           ,
                              int           , int           , int  );
  void remove_degree6_N      (Vertex_handle &v,
                              Face_handle & , Face_handle & , Face_handle & ,
                              Face_handle & , Face_handle & , Face_handle & ,
                              Vertex_handle&, Vertex_handle&, Vertex_handle&,
                              Vertex_handle&, Vertex_handle&, Vertex_handle&,
                              Offset&, Offset&, Offset&,
                              Offset&, Offset&, Offset&,
                              int           , int           , int           ,
                              int           , int           , int  );
  void remove_degree6_antiN  (Vertex_handle &v,
                              Face_handle & , Face_handle & , Face_handle & ,
                              Face_handle & , Face_handle & , Face_handle & ,
                              Vertex_handle&, Vertex_handle&, Vertex_handle&,
                              Vertex_handle&, Vertex_handle&, Vertex_handle&,
                              int           , int           , int           ,
                              int           , int           , int  );
  void remove_degree6_antiN  (Vertex_handle &v,
                              Face_handle & , Face_handle & , Face_handle & ,
                              Face_handle & , Face_handle & , Face_handle & ,
                              Vertex_handle&, Vertex_handle&, Vertex_handle&,
                              Vertex_handle&, Vertex_handle&, Vertex_handle&,
                              Offset&, Offset&, Offset&,
                              Offset&, Offset&, Offset&,
                              int           , int           , int           ,
                              int           , int           , int  );
  void remove_degree6_diamond(Vertex_handle &v,
                              Face_handle & , Face_handle & , Face_handle & ,
                              Face_handle & , Face_handle & , Face_handle & ,
                              Vertex_handle&, Vertex_handle&, Vertex_handle&,
                              Vertex_handle&, Vertex_handle&, Vertex_handle&,
                              int           , int           , int           ,
                              int           , int           , int  );
  void remove_degree6_diamond(Vertex_handle &v,
                              Face_handle & , Face_handle & , Face_handle & ,
                              Face_handle & , Face_handle & , Face_handle & ,
                              Vertex_handle&, Vertex_handle&, Vertex_handle&,
                              Vertex_handle&, Vertex_handle&, Vertex_handle&,
                              Offset&, Offset&, Offset&,
                              Offset&, Offset&, Offset&,
                              int           , int           , int           ,
                              int           , int           , int  );
  void remove_degree7(Vertex_handle v, std::vector<Face_handle> &f,
                      std::vector<Vertex_handle> &w,
                      std::vector<int> &i);
  void remove_degree7(Vertex_handle v, std::vector<Face_handle> &f,
                      std::vector<Vertex_handle> &w, std::vector<Offset> &o,
                      std::vector<int> &i);

  void rotate7(int j, std::vector<Vertex_handle> &w,
               std::vector<Face_handle> &f,
               std::vector<int> &i);
  void rotate7(int j, std::vector<Vertex_handle> &w,
               std::vector<Face_handle> &f, std::vector<Offset> &o,
               std::vector<int> &i);
  /// Returns whether the simplicity criterion is satisfied
  void get_offset_degree7(std::vector<Offset> &in_o, int out_o[]);
  void remove_degree7_star      (Vertex_handle&, int, std::vector<Face_handle> &f,
                                 std::vector<Vertex_handle> &w, std::vector<int> &i);
  void remove_degree7_star      (Vertex_handle&, int, std::vector<Face_handle> &f,
                                 std::vector<Vertex_handle> &w, std::vector<Offset> &o, std::vector<int> &i);
  void remove_degree7_zigzag    (Vertex_handle&, int, std::vector<Face_handle> &f,
                                 std::vector<Vertex_handle> &w, std::vector<int> &i);
  void remove_degree7_zigzag    (Vertex_handle&, int, std::vector<Face_handle> &f,
                                 std::vector<Vertex_handle> &w, std::vector<Offset> &o, std::vector<int> &i);
  void remove_degree7_leftdelta (Vertex_handle&, int, std::vector<Face_handle> &f,
                                 std::vector<Vertex_handle> &w, std::vector<int> &i);
  void remove_degree7_leftdelta (Vertex_handle&, int, std::vector<Face_handle> &f,
                                 std::vector<Vertex_handle> &w, std::vector<Offset> &o, std::vector<int> &i);
  void remove_degree7_rightdelta(Vertex_handle&, int, std::vector<Face_handle> &f,
                                 std::vector<Vertex_handle> &w, std::vector<int> &i);
  void remove_degree7_rightdelta(Vertex_handle&, int, std::vector<Face_handle> &f,
                                 std::vector<Vertex_handle> &w, std::vector<Offset> &o, std::vector<int> &i);
  void remove_degree7_leftfan   (Vertex_handle&, int, std::vector<Face_handle> &f,
                                 std::vector<Vertex_handle> &w, std::vector<int> &i);
  void remove_degree7_leftfan   (Vertex_handle&, int, std::vector<Face_handle> &f,
                                 std::vector<Vertex_handle> &w, std::vector<Offset> &o, std::vector<int> &i);
  void remove_degree7_rightfan  (Vertex_handle&, int, std::vector<Face_handle> &f,
                                 std::vector<Vertex_handle> &w, std::vector<int> &i);
  void remove_degree7_rightfan  (Vertex_handle&, int, std::vector<Face_handle> &f,
                                 std::vector<Vertex_handle> &w, std::vector<Offset> &o, std::vector<int> &i);

  /// Determines whether the point p lies on the (un-)bounded side of
  /// the circle through the vertices of f
  Oriented_side
  side_of_oriented_circle(Face_handle f,
                          const Point & p, bool perturb = false) const;
  /// Determines whether the point p lies on the (un-)bounded side of
  /// the circle through the points p0, p1 and p2
  Oriented_side
  side_of_oriented_circle(const Point &p0, const Point &p1, const Point &p2,
                          const Point &p, bool perturb) const;
  /// Determines whether the point (p,o) lies on the (un-)bounded side of
  /// the circle through the points (p0,o0), (p1,o1) and (p2,o2)
  Oriented_side
  side_of_oriented_circle(const Point &p0, const Point &p1, const Point &p2,
                          const Point &p, const Offset &o0, const Offset &o1, const Offset &o2,
                          const Offset &o, bool perturb) const;

  bool incircle(int x, int j, int k, int l, std::vector<Face_handle> &,
                std::vector<Vertex_handle> &w, std::vector<int> &)
  {

    return side_of_oriented_circle(w[j]->point(), w[k]->point(), w[l]->point(), w[x]->point(),
                                   true) ==  ON_POSITIVE_SIDE;
  }
  bool incircle(int x, int j, int k, int l, std::vector<Face_handle> &,
                std::vector<Vertex_handle> &w, std::vector<Offset> &o, std::vector<int> &)
  {

    return side_of_oriented_circle(w[j]->point(), w[k]->point(), w[l]->point(), w[x]->point(),
                                   o[j], o[k], o[l], o[x],
                                   true) ==  ON_POSITIVE_SIDE;
  }

// end of auxilliary functions for remove




  Vertex_handle nearest_vertex_2D(const Point& p, Face_handle f) const;

  void  look_nearest_neighbor(const Point& p,
                              Face_handle f,
                              int i,
                              Vertex_handle& nn) const;

  template <class OutputItFaces, class OutputItBoundaryEdges>
  std::pair<OutputItFaces, OutputItBoundaryEdges>
  propagate_conflicts (const Point  &p,
                       Face_handle fh,
                       int i,
                       std::pair<OutputItFaces, OutputItBoundaryEdges>
                       pit)  const
  {
    Face_handle fn = fh->neighbor(i);
    if (! test_conflict(p, fn))
      {
        *(pit.second)++ = Edge(fn, fn->index(fh));
      }
    else
      {
        *(pit.first)++ = fn;
        int j = fn->index(fh);
        pit = propagate_conflicts(p, fn, ccw(j), pit);
        pit = propagate_conflicts(p, fn, cw(j), pit);
      }
    return pit;
  }

  bool test_conflict(const Point &p, Face_handle fh) const
  {
    return side_of_oriented_circle(fh, p, true) ==  ON_POSITIVE_SIDE;
  }
};

template < class Gt, class Tds >
bool
Periodic_2_Delaunay_triangulation_2<Gt, Tds>::
is_valid(bool verbose, int level) const
{
  // Check the parent
  bool result = Periodic_2_triangulation_2<Gt, Tds>::is_valid(verbose, level);

  // Check in_sphere:
  if (dimension() == 2)
    {
      const Point *p[4];
      Offset off[4];
      for (Face_iterator fit = faces_begin();
           fit != this->faces_end(); ++fit)
        {
          for (int i = 0; i < 3; i++)
            {
              p[i] = &fit->vertex(i)->point();
              off[i] = get_offset(fit, i);
            }

          /// Check whether the vertices of the neighbor lie outside the circumcircle of the face
          for (int i = 0; i < 3; ++i)
            {
              p[3]   = &fit->vertex(i)->point();
              off[3] = combine_offsets(get_offset(fit, i), get_neighbor_offset(fit, i));

              result &= ON_POSITIVE_SIDE !=
                        side_of_oriented_circle(*p[0], *p[1], *p[2], *p[3],
                                                off[0], off[1], off[2], off[3],
                                                false);
              CGAL_triangulation_assertion(result);
            }
        }
    }

  return result;
}

template < class Gt, class Tds >
typename Periodic_2_Delaunay_triangulation_2<Gt, Tds>::Vertex_handle
Periodic_2_Delaunay_triangulation_2<Gt, Tds>::
nearest_vertex(const Point  &p, Face_handle f) const
{
  switch (dimension())
    {
    case 0:
      return Vertex_handle();
      //break;
    case 2:
      return nearest_vertex_2D(p, f);
      //break;
    }
  return Vertex_handle();
}

template < class Gt, class Tds >
typename Periodic_2_Delaunay_triangulation_2<Gt, Tds>::Vertex_handle
Periodic_2_Delaunay_triangulation_2<Gt, Tds>::
nearest_vertex_2D(const Point& p, Face_handle f) const
{
  CGAL_triangulation_precondition(dimension() == 2);
  f = locate(p, f);

  typename Geom_traits::Compare_distance_2 compare_distance =
    geom_traits().compare_distance_2_object();
  Vertex_handle nn =  f->vertex(0);
  if (compare_distance(p, f->vertex(1)->point(), nn->point()) == SMALLER)
    nn = f->vertex(1);
  if (compare_distance(p, f->vertex(2)->point(), nn->point()) == SMALLER)
    nn = f->vertex(2);

  look_nearest_neighbor(p, f, 0, nn);
  look_nearest_neighbor(p, f, 1, nn);
  look_nearest_neighbor(p, f, 2, nn);

  return nn;
}

template < class Gt, class Tds >
void
Periodic_2_Delaunay_triangulation_2<Gt, Tds>::
look_nearest_neighbor(const Point& p,
                      Face_handle f,
                      int i,
                      Vertex_handle& nn) const
{
  Face_handle  ni = f->neighbor(i);
  if ( this->side_of_oriented_circle(ni, p, true) != ON_POSITIVE_SIDE )
    return;

  typename Geom_traits::Compare_distance_2 compare_distance =
    geom_traits().compare_distance_2_object();
  i = ni->index(f);
  if (compare_distance(p, ni->vertex(i)->point(), nn->point()) == SMALLER)
    nn = ni->vertex(i);

  // recursive exploration of triangles whose circumcircle contains p
  look_nearest_neighbor(p, ni, ccw(i), nn);
  look_nearest_neighbor(p, ni, cw(i), nn);
}

//DUALITY
template<class Gt, class Tds>
inline
typename Periodic_2_Delaunay_triangulation_2<Gt, Tds>::Point
Periodic_2_Delaunay_triangulation_2<Gt, Tds>::
dual (Face_handle f) const
{
  CGAL_triangulation_precondition (dimension() == 2);
  return circumcenter(f);
}


template < class Gt, class Tds >
inline typename Gt::Segment_2
Periodic_2_Delaunay_triangulation_2<Gt, Tds>::
dual(const Edge &e) const
{
  // dimension==2
  Face_handle nb = e.first->neighbor(e.second);
  Point p0 = dual(e.first);
  Point p1 = dual(nb);
  Offset o = combine_offsets( Offset(), get_neighbor_offset(e.first, e.second));
  Segment s = construct_segment(p0, p1, o, Offset());

  return s;
}

template < class Gt, class Tds >
inline typename Gt::Segment_2
Periodic_2_Delaunay_triangulation_2<Gt, Tds>::
dual(const Edge_circulator& ec) const
{
  return dual(*ec);
}
template < class Gt, class Tds >
inline typename Gt::Segment_2
Periodic_2_Delaunay_triangulation_2<Gt, Tds>::
dual(const Edge_iterator& ei) const
{
  return dual(*ei);
}

///////////////////////////////////////////////////////////////
//  INSERT

template < class Gt, class Tds >
inline
typename Periodic_2_Delaunay_triangulation_2<Gt, Tds>::Vertex_handle
Periodic_2_Delaunay_triangulation_2<Gt, Tds>::
insert(const Point  &p,  Face_handle start)
{
  CGAL_triangulation_assertion((domain().xmin() <= p.x()) &&
                               (p.x() < domain().xmax()));
  CGAL_triangulation_assertion((domain().ymin() <= p.y()) &&
                               (p.y() < domain().ymax()));

  if (empty())
    {
      return this->insert_first(p);
    }

  if (start == Face_handle())
    {
      start = this->faces_begin();
    }

  Locate_type lt;
  int li;
  Face_handle loc = locate (p, lt, li, start);

  /// Call the insert function with the located simplex
  return insert(p, lt, loc, li);
}

template < class Gt, class Tds >
inline
typename Periodic_2_Delaunay_triangulation_2<Gt, Tds>::Vertex_handle
Periodic_2_Delaunay_triangulation_2<Gt, Tds>::
push_back(const Point &p)
{
  return insert(p);
}

template < class Gt, class Tds >
inline
typename Periodic_2_Delaunay_triangulation_2<Gt, Tds>::Vertex_handle
Periodic_2_Delaunay_triangulation_2<Gt, Tds>::
insert(const Point  &p, Locate_type lt, Face_handle loc, int li)
{
  Vertex_handle vh = Base::insert(p, lt, loc, li);

  if (lt != Base::VERTEX)
    {
      restore_Delaunay(vh);

      if (!is_1_cover())
        {
          typename Base::Virtual_vertex_reverse_map_it vertices_it =
            this->virtual_vertices_reverse().find(vh);
          CGAL_triangulation_assertion(vertices_it != this->virtual_vertices_reverse().end());
          const std::vector<Vertex_handle> &virtual_vertices = vertices_it->second;
          for (size_t i = 0; i < virtual_vertices.size(); ++i)
            {
              restore_Delaunay(virtual_vertices[i]);
            }

          this->try_to_convert_to_one_cover();
          if (is_1_cover())
            {
              CGAL_triangulation_assertion(is_valid());
            }
        }
    }

  return vh;
}

template < class Gt, class Tds >
void
Periodic_2_Delaunay_triangulation_2<Gt, Tds>::
restore_Delaunay(Vertex_handle v)
{
  Face_handle f = v->face();
  Face_handle next;
  int i;
  Face_handle start(f);
  do
    {
      i = f->index(v);
      next = f->neighbor(ccw(i));  // turn ccw around v
      propagating_flip(f, i);
      f = next;
    }
  while(next != start);
}


#ifndef CGAL_DT2_USE_RECURSIVE_PROPAGATING_FLIP
template < class Gt, class Tds >
void
Periodic_2_Delaunay_triangulation_2<Gt, Tds>::
non_recursive_propagating_flip(Face_handle f , int i)
{
  std::stack<Edge> edges;
  const Vertex_handle& vp = f->vertex(i);
  edges.push(Edge(f, i));

  while(! edges.empty())
    {
      const Edge& e = edges.top();
      f = e.first;
      i = e.second;
      const Face_handle& n = f->neighbor(i);

      if (locally_Delaunay(f, i, n))
        {
          edges.pop();
          continue;
        }
      this->flip_single_edge(f, i);
      // As we haven't popped it, we don't have to push it
      edges.push(Edge(n, n->index(vp)));
    }
}

template < class Gt, class Tds >
void
Periodic_2_Delaunay_triangulation_2<Gt, Tds>::
propagating_flip(const Face_handle& f, int i, int depth)
{
#ifdef CGAL_DT2_IMMEDIATELY_NON_RECURSIVE_PROPAGATING_FLIP
  non_recursive_propagating_flip(f, i);
#else
  int max_depth = 100;
  if(depth == max_depth)
    {
      non_recursive_propagating_flip(f, i);
      return;
    }
  Face_handle n = f->neighbor(i);

  if (locally_Delaunay(f, i, n))
    return;

  this->flip_single_edge(f, i);
  propagating_flip(f, i, depth + 1);
  i = n->index(f->vertex(i));
  propagating_flip(n, i, depth + 1);
#endif
}
#else
template < class Gt, class Tds >
void
Periodic_2_Delaunay_triangulation_2<Gt, Tds>::
propagating_flip(Face_handle& f, int i)
{
  Face_handle nb = f->neighbor(i);

  if (locally_Delaunay(f, nb))
    return;

  this->flip_single_edge(f, i);
  propagating_flip(f, i);
  i = nb->index(f->vertex(i));
  propagating_flip(nb, i);
}
#endif

template < class Gt, class Tds >
bool
Periodic_2_Delaunay_triangulation_2<Gt, Tds>::
locally_Delaunay(const Face_handle &f, int i, const Face_handle &nb)
{
  CGAL_BRANCH_PROFILER("locally_Delaunay(), simplicity check failures", tmp);

  bool simplicity_criterion = is_1_cover() && f->has_zero_offsets() && nb->has_zero_offsets();

  const Point *p[4];
  for (int index = 0; index < 3; ++index)
    {
      p[index]   = &nb->vertex(index)->point();
    }
  p[3]   = &f->vertex(i)->point();

  Oriented_side os;
  if (simplicity_criterion)
    {
      // No periodic offsets
      os = side_of_oriented_circle(*p[0], *p[1], *p[2], *p[3], true);
    }
  else
    {
      CGAL_BRANCH_PROFILER_BRANCH(tmp);

      Offset off[4];

      for (int index = 0; index < 3; ++index)
        {
          off[index] = get_offset(nb, index);
        }
      off[3] = combine_offsets(get_offset(f, i), get_neighbor_offset(f, i));

      os = side_of_oriented_circle(*p[0], *p[1], *p[2], *p[3],
                                   off[0], off[1], off[2], off[3], true);
    }

  return (ON_POSITIVE_SIDE != os);
}

///////////////////////////////////////////////////////////////
//  REMOVE    see INRIA RResearch Report 7104

template < class Gt, class Tds >
void
Periodic_2_Delaunay_triangulation_2<Gt, Tds>::
remove(Vertex_handle v)
{
  // Make sure we have the original vertex
  CGAL_assertion(v == this->get_original_vertex(v));

  CGAL_triangulation_precondition(v != Vertex_handle());
  CGAL_triangulation_precondition(dimension() == 2);

  if ( this->number_of_vertices() == 1)
    {
      // Last vertex
      Base::remove_first(v);
      return;
    }

  if (!remove_single_vertex(v, Offset()))
    {
      // The vertex was not removed as we need to revert to the 9-cover first
      this->convert_to_9_sheeted_covering();

      remove_single_vertex(v, Offset());
    }

  if (!is_1_cover())
    {
      CGAL_assertion(this->virtual_vertices_reverse().find(v) != this->virtual_vertices_reverse().end());

      const std::vector<Vertex_handle> &virtual_copies = this->virtual_vertices_reverse().find(v)->second;
      for (int i = 0; i < 8; ++i)
        {
          remove_single_vertex(virtual_copies[i], Offset((i + 1) / 3, (i + 1) % 3));
        }

      this->remove_from_virtual_copies(v);
    }
}

namespace internal{
namespace P2DT2{

template<class P2DT2>
struct Static_data{
  int maxd;
  std::vector<typename P2DT2::Face_handle> f;
  std::vector<int> i;
  std::vector<typename P2DT2::Vertex_handle> w;
  std::vector<typename P2DT2::Offset> offset_w;
  Static_data(int m)
    : maxd(m)
    , f(maxd)
    , i(maxd)
    , w(maxd)
    , offset_w(maxd)
  {}
};

} } //end of namespace internal::P2DT2

template < class Gt, class Tds >
bool
Periodic_2_Delaunay_triangulation_2<Gt, Tds>::
remove_single_vertex(Vertex_handle v, const Offset &v_o)
{
  typedef internal::P2DT2::
    Static_data< Periodic_2_Delaunay_triangulation_2<Gt, Tds> > Static_data;
  CGAL_STATIC_THREAD_LOCAL_VARIABLE(Static_data, sd, 30);

  int d;
  bool simplicity_criterion;

  if (remove_degree_init(v, v_o, sd.f, sd.w, sd.offset_w, sd.i, d, sd.maxd, simplicity_criterion))
    {
      if (is_1_cover())
        {
          // Don't delete if the hole is too big and the triangulation is a 1-cover
          return false;
        }
    }

  if (simplicity_criterion)
    remove_degree_triangulate(v, sd.f, sd.w, sd.i, d);
  else
    remove_degree_triangulate(v, sd.f, sd.w, sd.offset_w, sd.i, d);

  this->delete_vertex(v);

  return true;
}


template < class Gt, class Tds >
void
Periodic_2_Delaunay_triangulation_2<Gt, Tds>::
remove_degree_triangulate(Vertex_handle v,
                          std::vector<Face_handle> &f,
                          std::vector<Vertex_handle> &w,
                          std::vector<int> &i, int d)
{
  // degree: 3: 1%, 4: 9%, 5: 23%, 6: 35%, 7: 19%, r: 10%

  switch (d)
    {
    case 3:
      remove_degree3(v, f, w, i);
      break;
    case 4:
      remove_degree4(v, f, w, i);
      break;
    case 5:
      remove_degree5(v, f, w, i);
      break;
    case 6:
      remove_degree6(v, f, w, i);
      break;
    case 7:
      remove_degree7(v, f, w, i);
      break;
    default:
      remove_degree_d(v, f, w, i, d);
      break;
    }
}

template < class Gt, class Tds >
void
Periodic_2_Delaunay_triangulation_2<Gt, Tds>::
remove_degree_triangulate(Vertex_handle v,
                          std::vector<Face_handle> &f,
                          std::vector<Vertex_handle> &w,
                          std::vector<Offset> &offset_w,
                          std::vector<int> &i, int d)
{
  // degree: 3: 1%, 4: 9%, 5: 23%, 6: 35%, 7: 19%, r: 10%

  // Remove all the edges that are too long.
  // This only needs to be done when the simplicity condition is not
  // met because the simplicity condition implies is_1_cover(), hence
  // no too long edges.
  this->remove_too_long_edges_in_star(v);

  switch (d)
    {
    case 3:
      remove_degree3(v, f, w, offset_w, i);
      break;
    case 4:
      remove_degree4(v, f, w, offset_w, i);
      break;
    case 5:
      remove_degree5(v, f, w, offset_w, i);
      break;
    case 6:
      remove_degree6(v, f, w, offset_w, i);
      break;
    case 7:
      remove_degree7(v, f, w, offset_w, i);
      break;
    default:
      remove_degree_d(v, f, w, offset_w, i, d);
      break;
    }
}

template < class Gt, class Tds >
void
Periodic_2_Delaunay_triangulation_2<Gt, Tds>::
remove_degree_d(Vertex_handle v, std::vector<Face_handle> &,
                std::vector<Vertex_handle> &,
                std::vector<int> &, int)
{
  std::list<Edge> hole;
  this->make_hole(v, hole);

  fill_hole_delaunay(hole);

  return;
}
template < class Gt, class Tds >
void
Periodic_2_Delaunay_triangulation_2<Gt, Tds>::
remove_degree_d(Vertex_handle v, std::vector<Face_handle> &,
                std::vector<Vertex_handle> &w, std::vector<Offset> &offset_w,
                std::vector<int> &, int d)
{
  std::list<Edge> hole;
  this->make_hole(v, hole);

  std::map<Vertex_handle, Offset> vertex_offsets;
  for (int idx = 0; idx < d; ++idx)
    {
      vertex_offsets[w[idx]] = offset_w[idx];
    }

  fill_hole_delaunay(hole, vertex_offsets);

  return;
}

template < class Gt, class Tds >
void
Periodic_2_Delaunay_triangulation_2<Gt, Tds>::
remove_degree3(Vertex_handle, std::vector<Face_handle> &f,
               std::vector<Vertex_handle> &,
               std::vector<int> &i)
{
  // modify the triangulation
  Face_handle nn = f[1]->neighbor( i[1] );
  tds().set_adjacency(f[0], ccw(i[0]) , nn , nn->index(f[1])  );
  nn = f[2]->neighbor( i[2] );
  tds().set_adjacency(f[0], cw(i[0]) , nn , nn->index(f[2])  );
  f[0]->set_vertex  (            i[0] , f[1]->vertex( cw(i[1]) ) );

  // clean container
  tds().delete_face(f[1]);
  tds().delete_face(f[2]);

  this->set_offsets(f[0], 0, 0, 0);

  return;
}
template < class Gt, class Tds >
void
Periodic_2_Delaunay_triangulation_2<Gt, Tds>::
remove_degree3(Vertex_handle, std::vector<Face_handle> &f,
               std::vector<Vertex_handle> &, std::vector<Offset> &o,
               std::vector<int> &i)
{
  // modify the triangulation
  Face_handle nn = f[1]->neighbor( i[1] );
  tds().set_adjacency(f[0], ccw(i[0]) , nn , nn->index(f[1])  );
  nn = f[2]->neighbor( i[2] );
  tds().set_adjacency(f[0], cw(i[0]) , nn , nn->index(f[2])  );
  f[0]->set_vertex  (            i[0] , f[1]->vertex( cw(i[1]) ) );

  // clean container
  tds().delete_face(f[1]);
  tds().delete_face(f[2]);

  Offset oo[3];
  oo[    i[0] ] = o[2];
  oo[ cw(i[0])] = o[1];
  oo[ccw(i[0])] = o[0];

  if (oo[0].x() < 0 || oo[1].x() < 0 || oo[2].x() < 0)
    {
      Offset o(number_of_sheets()[0], 0);
      oo[0] += o;
      oo[1] += o;
      oo[2] += o;
    }
  if (oo[0].y() < 0 || oo[1].y() < 0 || oo[2].y() < 0)
    {
      Offset o(0, number_of_sheets()[1]);
      oo[0] += o;
      oo[1] += o;
      oo[2] += o;
    }
  this->set_offsets(f[0],
                    (oo[0].x() >= number_of_sheets()[0] ? 2 : 0) + (oo[0].y() >= number_of_sheets()[1] ? 1 : 0),
                    (oo[1].x() >= number_of_sheets()[0] ? 2 : 0) + (oo[1].y() >= number_of_sheets()[1] ? 1 : 0),
                    (oo[2].x() >= number_of_sheets()[0] ? 2 : 0) + (oo[2].y() >= number_of_sheets()[1] ? 1 : 0));

  return;
}

template < class Gt, class Tds >
void
Periodic_2_Delaunay_triangulation_2<Gt, Tds>::
remove_degree4(Vertex_handle, std::vector<Face_handle> &f,
               std::vector<Vertex_handle> &w,
               std::vector<int> &i )
{
  // removing a degree 4 vertex

  Face_handle nn;

  if ( !incircle(2, 0, 1, 3, f, w, i) )
    {
      // diagonal 1 3
      f[0]->set_vertex( i[0], w[3] ); //w0 w1 w3
      f[1]->set_vertex( i[1], w[3] ); //w1 w2 w3
      nn = f[3]->neighbor( i[3] );
      tds().set_adjacency(f[0], cw(i[0]) , nn , nn->index(f[3])  );
      nn = f[2]->neighbor( i[2] );
      tds().set_adjacency(f[1], ccw(i[1]) , nn , nn->index(f[2]) );
      // clean container
      tds().delete_face(f[2]);
      tds().delete_face(f[3]);

      f[0]->set_offsets(0, 0, 0);
      f[1]->set_offsets(0, 0, 0);

      insert_too_long_edge(f[0], ccw(i[0]));
    }
  else
    {
      // diagonal 0 2
      f[0]->set_vertex( i[0], w[2]); //w0 w1 w2
      f[3]->set_vertex( i[3], w[2]); //w3 w0 w2
      nn = f[1]->neighbor( i[1] );
      tds().set_adjacency(f[0], ccw(i[0]) , nn , nn->index(f[1])  );
      nn = f[2]->neighbor( i[2] );
      tds().set_adjacency(f[3], cw(i[3]) , nn , nn->index(f[2])  );
      // clean container
      tds().delete_face(f[1]);
      tds().delete_face(f[2]);

      f[0]->set_offsets(0, 0, 0);
      f[3]->set_offsets(0, 0, 0);
    }
}

template < class Gt, class Tds >
void
Periodic_2_Delaunay_triangulation_2<Gt, Tds>::
remove_degree4(Vertex_handle, std::vector<Face_handle> &f,
               std::vector<Vertex_handle> &w, std::vector<Offset> &o,
               std::vector<int> &i )
{
  // removing a degree 4 vertex

  Face_handle nn;

  int oo[4];
  if ((o[0] == o[1]) && (o[0] == o[2]) && (o[0] == o[3]))
    {
      for (int i = 0; i < 4; ++i) oo[i] = 0;
    }
  else
    {
      Covering_sheets cover = number_of_sheets();
      if ((o[0].x() < 0) || (o[1].x() < 0) || (o[2].x() < 0) || (o[3].x() < 0))
        for (int i = 0; i < 4; ++i)
          o[i] += Offset(cover[0], 0);

      if ((o[0].y() < 0) || (o[1].y() < 0) || (o[2].y() < 0) || (o[3].y() < 0))
        for (int i = 0; i < 4; ++i)
          o[i] += Offset(0, cover[1]);

      for (int i = 0; i < 4; ++i)
        {
          oo[i] = (o[i].x() >= cover[0] ? 2 : 0) + (o[i].y() >= cover[1] ? 1 : 0);
        }
    }

  if ( !incircle(2, 0, 1, 3, f, w, o, i) )
    {
      // diagonal 1 3
      f[0]->set_vertex( i[0], w[3] ); //w0 w1 w3
      f[1]->set_vertex( i[1], w[3] ); //w1 w2 w3
      nn = f[3]->neighbor( i[3] );
      tds().set_adjacency(f[0], cw(i[0]) , nn , nn->index(f[3])  );
      nn = f[2]->neighbor( i[2] );
      tds().set_adjacency(f[1], ccw(i[1]) , nn , nn->index(f[2]) );
      // clean container
      tds().delete_face(f[2]);
      tds().delete_face(f[3]);

      int o_face[3];
      o_face[i[0]]      = oo[3];
      o_face[ccw(i[0])] = oo[0];
      o_face[ cw(i[0])] = oo[1];
      this->set_offsets(f[0], o_face[0], o_face[1], o_face[2]);
      o_face[i[1]]      = oo[3];
      o_face[ccw(i[1])] = oo[1];
      o_face[ cw(i[1])] = oo[2];
      this->set_offsets(f[1], o_face[0], o_face[1], o_face[2]);

      insert_too_long_edge(f[0], ccw(i[0]));
    }
  else
    {
      // diagonal 0 2
      f[0]->set_vertex( i[0], w[2]); //w0 w1 w2
      f[3]->set_vertex( i[3], w[2]); //w3 w0 w2
      nn = f[1]->neighbor( i[1] );
      tds().set_adjacency(f[0], ccw(i[0]) , nn , nn->index(f[1])  );
      nn = f[2]->neighbor( i[2] );
      tds().set_adjacency(f[3], cw(i[3]) , nn , nn->index(f[2])  );
      // clean container
      tds().delete_face(f[1]);
      tds().delete_face(f[2]);

      int o_face[3];
      o_face[i[0]]      = oo[2];
      o_face[ccw(i[0])] = oo[0];
      o_face[ cw(i[0])] = oo[1];
      this->set_offsets(f[0], o_face[0], o_face[1], o_face[2]);
      o_face[i[3]]      = oo[2];
      o_face[ccw(i[3])] = oo[3];
      o_face[ cw(i[3])] = oo[0];
      this->set_offsets(f[3], o_face[0], o_face[1], o_face[2]);

      insert_too_long_edge(f[3], ccw(i[3]));
    }
}

template < class Gt, class Tds >
void
Periodic_2_Delaunay_triangulation_2<Gt, Tds>::
remove_degree5(Vertex_handle v, std::vector<Face_handle> &f,
               std::vector<Vertex_handle> &w,
               std::vector<int> &i)
{
  // removing a degree 5 vertex
  this->remove_too_long_edges_in_star(v);

  if (incircle(3, 0, 1, 2, f, w, i))
    {
      if (incircle(4, 0, 1, 3, f, w, i))
        {
          if (incircle(4, 1, 2, 3, f, w, i))
            {
              // star from 4
              remove_degree5_star(v, f[4], f[0], f[1], f[2], f[3],
                                  w[4], w[0], w[1], w[2], w[3],
                                  i[4], i[0], i[1], i[2], i[3]);
            }
          else
            {
              //star from 1
              remove_degree5_star(v, f[1], f[2], f[3], f[4], f[0],
                                  w[1], w[2], w[3], w[4], w[0],
                                  i[1], i[2], i[3], i[4], i[0]);


            }
        }
      else
        {
          // star from 3
          remove_degree5_star(v, f[3], f[4], f[0], f[1], f[2],
                              w[3], w[4], w[0], w[1], w[2],
                              i[3], i[4], i[0], i[1], i[2]);
        }
    }
  else
    {
      if (incircle(4, 2, 3, 0, f, w, i))
        {
          if (incircle(4, 0, 1, 2, f, w, i))
            {
              // star from 4
              remove_degree5_star(v, f[4], f[0], f[1], f[2], f[3],
                                  w[4], w[0], w[1], w[2], w[3],
                                  i[4], i[0], i[1], i[2], i[3]);
            }
          else
            {
              //star from 2
              remove_degree5_star(v, f[2], f[3], f[4], f[0], f[1],
                                  w[2], w[3], w[4], w[0], w[1],
                                  i[2], i[3], i[4], i[0], i[1]);
            }
        }
      else
        {
          // star from 0
          remove_degree5_star(v, f[0], f[1], f[2], f[3], f[4],
                              w[0], w[1], w[2], w[3], w[4],
                              i[0], i[1], i[2], i[3], i[4]);
        }
    }
  return;
}
template < class Gt, class Tds >
void
Periodic_2_Delaunay_triangulation_2<Gt, Tds>::
remove_degree5(Vertex_handle v, std::vector<Face_handle> &f,
               std::vector<Vertex_handle> &w, std::vector<Offset> &o,
               std::vector<int> &i)
{
  // removing a degree 5 vertex
  this->remove_too_long_edges_in_star(v);

  if (incircle(3, 0, 1, 2, f, w, o, i))
    {
      if (incircle(4, 0, 1, 3, f, w, o, i))
        {
          if (incircle(4, 1, 2, 3, f, w, o, i))
            {
              // star from 4
              remove_degree5_star(v, f[4], f[0], f[1], f[2], f[3],
                                  w[4], w[0], w[1], w[2], w[3],
                                  o[4], o[0], o[1], o[2], o[3],
                                  i[4], i[0], i[1], i[2], i[3]);
            }
          else
            {
              //star from 1
              remove_degree5_star(v, f[1], f[2], f[3], f[4], f[0],
                                  w[1], w[2], w[3], w[4], w[0],
                                  o[1], o[2], o[3], o[4], o[0],
                                  i[1], i[2], i[3], i[4], i[0]);


            }
        }
      else
        {
          // star from 3
          remove_degree5_star(v, f[3], f[4], f[0], f[1], f[2],
                              w[3], w[4], w[0], w[1], w[2],
                              o[3], o[4], o[0], o[1], o[2],
                              i[3], i[4], i[0], i[1], i[2]);
        }
    }
  else
    {
      if (incircle(4, 2, 3, 0, f, w, o, i))
        {
          if (incircle(4, 0, 1, 2, f, w, o, i))
            {
              // star from 4
              remove_degree5_star(v, f[4], f[0], f[1], f[2], f[3],
                                  w[4], w[0], w[1], w[2], w[3],
                                  o[4], o[0], o[1], o[2], o[3],
                                  i[4], i[0], i[1], i[2], i[3]);
            }
          else
            {
              //star from 2
              remove_degree5_star(v, f[2], f[3], f[4], f[0], f[1],
                                  w[2], w[3], w[4], w[0], w[1],
                                  o[2], o[3], o[4], o[0], o[1],
                                  i[2], i[3], i[4], i[0], i[1]);
            }
        }
      else
        {
          // star from 0
          remove_degree5_star(v, f[0], f[1], f[2], f[3], f[4],
                              w[0], w[1], w[2], w[3], w[4],
                              o[0], o[1], o[2], o[3], o[4],
                              i[0], i[1], i[2], i[3], i[4]);
        }
    }
  return;
}

template < class Gt, class Tds >
void
Periodic_2_Delaunay_triangulation_2<Gt, Tds>::remove_degree5_star
(Vertex_handle &,
 Face_handle &f0, Face_handle &f1, Face_handle &f2, Face_handle &f3, Face_handle &f4,
 Vertex_handle &v0, Vertex_handle &, Vertex_handle &, Vertex_handle &, Vertex_handle &,
 int i0, int i1, int i2, int i3, int i4 )
{
  // removing a degree 5 vertex, starring from v0
  Face_handle nn;
  f1->set_vertex( i1, v0) ;  // f1 = v1v2v0
  f2->set_vertex( i2, v0) ;  // f2 = v2v3v0
  f3->set_vertex( i3, v0) ;  // f3 = v3v4v0

  nn = f0->neighbor( i0 );
  tds().set_adjacency(f1, cw(i1) , nn , nn->index(f0) );
  nn = f4->neighbor( i4 );
  tds().set_adjacency(f3, ccw(i3) , nn , nn->index(f4) );
  tds().delete_face(f0);
  tds().delete_face(f4);

  f1->set_offsets(0, 0, 0);
  f2->set_offsets(0, 0, 0);
  f3->set_offsets(0, 0, 0);
}
template < class Gt, class Tds >
void
Periodic_2_Delaunay_triangulation_2<Gt, Tds>::remove_degree5_star
(Vertex_handle &,
 Face_handle &f0, Face_handle &f1, Face_handle &f2, Face_handle &f3, Face_handle &f4,
 Vertex_handle &v0, Vertex_handle &, Vertex_handle &, Vertex_handle &, Vertex_handle &,
 Offset &o0, Offset &o1, Offset &o2, Offset &o3, Offset &o4,
 int i0, int i1, int i2, int i3, int i4 )
{
  // removing a degree 5 vertex, starring from v0
  Face_handle nn;
  f1->set_vertex( i1, v0) ;  // f1 = v1v2v0
  f2->set_vertex( i2, v0) ;  // f2 = v2v3v0
  f3->set_vertex( i3, v0) ;  // f3 = v3v4v0

  nn = f0->neighbor( i0 );
  tds().set_adjacency(f1, cw(i1) , nn , nn->index(f0) );
  nn = f4->neighbor( i4 );
  tds().set_adjacency(f3, ccw(i3) , nn , nn->index(f4) );
  tds().delete_face(f0);
  tds().delete_face(f4);

  if (o0.x() < 0 || o1.x() < 0 || o2.x() < 0 || o3.x() < 0 || o4.x() < 0)
    {
      o0 += Offset(number_of_sheets()[0], 0);
      o1 += Offset(number_of_sheets()[0], 0);
      o2 += Offset(number_of_sheets()[0], 0);
      o3 += Offset(number_of_sheets()[0], 0);
      o4 += Offset(number_of_sheets()[0], 0);
    }
  if (o0.y() < 0 || o1.y() < 0 || o2.y() < 0 || o3.y() < 0 || o4.y() < 0)
    {
      o0 += Offset(0, number_of_sheets()[1]);
      o1 += Offset(0, number_of_sheets()[1]);
      o2 += Offset(0, number_of_sheets()[1]);
      o3 += Offset(0, number_of_sheets()[1]);
      o4 += Offset(0, number_of_sheets()[1]);
    }
  int oo0 = (o0.x() >= number_of_sheets()[0] ? 2 : 0) + (o0.y() >= number_of_sheets()[1] ? 1 : 0);
  int oo1 = (o1.x() >= number_of_sheets()[0] ? 2 : 0) + (o1.y() >= number_of_sheets()[1] ? 1 : 0);
  int oo2 = (o2.x() >= number_of_sheets()[0] ? 2 : 0) + (o2.y() >= number_of_sheets()[1] ? 1 : 0);
  int oo3 = (o3.x() >= number_of_sheets()[0] ? 2 : 0) + (o3.y() >= number_of_sheets()[1] ? 1 : 0);
  int oo4 = (o4.x() >= number_of_sheets()[0] ? 2 : 0) + (o4.y() >= number_of_sheets()[1] ? 1 : 0);

  int oo[3];
  oo[i1]      = oo0;
  oo[ccw(i1)] = oo1;
  oo[ cw(i1)] = oo2;
  this->set_offsets(f1, oo[0], oo[1], oo[2]);
  oo[i2]      = oo0;
  oo[ccw(i2)] = oo2;
  oo[ cw(i2)] = oo3;
  this->set_offsets(f2, oo[0], oo[1], oo[2]);
  oo[i3]      = oo0;
  oo[ccw(i3)] = oo3;
  oo[ cw(i3)] = oo4;
  this->set_offsets(f3, oo[0], oo[1], oo[2]);

  //insert_too_long_edges_in_star(f1->vertex(i1));
  insert_too_long_edge(f1, ccw(i1));
  insert_too_long_edge(f2, ccw(i2));
}

template < class Gt, class Tds >
void
Periodic_2_Delaunay_triangulation_2<Gt, Tds>::
remove_degree6(Vertex_handle v, std::vector<Face_handle> &f,
               std::vector<Vertex_handle> &w,
               std::vector<int> &i)
{
  if(incircle(1, 2, 3, 0, f, w, i))
    {
      if(incircle(4, 2, 3, 5, f, w, i))
        {
          if(incircle(1, 2, 3, 4, f, w, i))
            {
              if(incircle(4, 0, 1, 3, f, w, i))
                {
                  if(incircle(5, 0, 1, 4, f, w, i))
                    {
                      remove_degree6_star(v,
                                          f[1], f[2], f[3], f[4], f[5], f[0],
                                          w[1], w[2], w[3], w[4], w[5], w[0],
                                          i[1], i[2], i[3], i[4], i[5], i[0]);
                    }
                  else
                    {
                      remove_degree6_N(v,
                                       f[1], f[2], f[3], f[4], f[5], f[0],
                                       w[1], w[2], w[3], w[4], w[5], w[0],
                                       i[1], i[2], i[3], i[4], i[5], i[0]);
                    }
                }
              else
                {
                  remove_degree6_antiN(v,
                                       f[0], f[1], f[2], f[3], f[4], f[5],
                                       w[0], w[1], w[2], w[3], w[4], w[5],
                                       i[0], i[1], i[2], i[3], i[4], i[5]);
                }
            }
          else
            {
              if(incircle(5, 1, 2, 4, f, w, i))
                {
                  remove_degree6_N(v, f[2], f[3], f[4], f[5], f[0], f[1],
                                   w[2], w[3], w[4], w[5], w[0], w[1],
                                   i[2], i[3], i[4], i[5], i[0], i[1]);
                }
              else
                {
                  if(incircle(5, 0, 1, 4, f, w, i))
                    {
                      remove_degree6_antiN(v,
                                           f[1], f[2], f[3], f[4], f[5], f[0],
                                           w[1], w[2], w[3], w[4], w[5], w[0],
                                           i[1], i[2], i[3], i[4], i[5], i[0]);
                    }
                  else
                    {
                      remove_degree6_star(v, f[4], f[5], f[0], f[1], f[2], f[3],
                                          w[4], w[5], w[0], w[1], w[2], w[3],
                                          i[4], i[5], i[0], i[1], i[2], i[3]);
                    }
                }
            }
        }
      else
        {
          if(incircle(1, 2, 3, 5, f, w, i))
            {
              if(incircle(1, 3, 4, 5, f, w, i))
                {
                  if(incircle(4, 0, 1, 3, f, w, i))
                    {
                      if(incircle(5, 0, 1, 4, f, w, i))
                        {
                          remove_degree6_star(v, f[1], f[2], f[3], f[4], f[5], f[0],
                                              w[1], w[2], w[3], w[4], w[5], w[0],
                                              i[1], i[2], i[3], i[4], i[5], i[0]);
                        }
                      else
                        {
                          remove_degree6_N(v, f[1], f[2], f[3], f[4], f[5], f[0],
                                           w[1], w[2], w[3], w[4], w[5], w[0],
                                           i[1], i[2], i[3], i[4], i[5], i[0]);
                        }
                    }
                  else
                    {
                      remove_degree6_antiN(v, f[0], f[1], f[2], f[3], f[4], f[5],
                                           w[0], w[1], w[2], w[3], w[4], w[5],
                                           i[0], i[1], i[2], i[3], i[4], i[5]);
                    }
                }
              else
                {
                  if(incircle(5, 0, 1, 3, f, w, i))
                    {
                      remove_degree6_diamond(v, f[1], f[2], f[3], f[4], f[5], f[0],
                                             w[1], w[2], w[3], w[4], w[5], w[0],
                                             i[1], i[2], i[3], i[4], i[5], i[0]);
                    }
                  else
                    {
                      if(incircle(4, 5, 0, 3, f, w, i))
                        {
                          remove_degree6_antiN(v, f[0], f[1], f[2], f[3], f[4], f[5],
                                               w[0], w[1], w[2], w[3], w[4], w[5],
                                               i[0], i[1], i[2], i[3], i[4], i[5]);
                        }
                      else
                        {
                          remove_degree6_star(v, f[3], f[4], f[5], f[0], f[1], f[2],
                                              w[3], w[4], w[5], w[0], w[1], w[2],
                                              i[3], i[4], i[5], i[0], i[1], i[2]);
                        }
                    }
                }
            }
          else
            {
              remove_degree6_star(v, f[5], f[0], f[1], f[2], f[3], f[4],
                                  w[5], w[0], w[1], w[2], w[3], w[4],
                                  i[5], i[0], i[1], i[2], i[3], i[4]);
            }
        }
    }
  else
    {
      if(incircle(4, 2, 3, 5, f, w, i))
        {
          if(incircle(4, 2, 3, 0, f, w, i))
            {
              if(incircle(4, 0, 1, 2, f, w, i))
                {
                  if(incircle(4, 1, 2, 5, f, w, i))
                    {
                      if(incircle(4, 0, 1, 5, f, w, i))
                        {
                          remove_degree6_star(v, f[4], f[5], f[0], f[1], f[2], f[3],
                                              w[4], w[5], w[0], w[1], w[2], w[3],
                                              i[4], i[5], i[0], i[1], i[2], i[3]);
                        }
                      else
                        {
                          remove_degree6_antiN(v, f[1], f[2], f[3], f[4], f[5], f[0],
                                               w[1], w[2], w[3], w[4], w[5], w[0],
                                               i[1], i[2], i[3], i[4], i[5], i[0]);
                        }
                    }
                  else
                    {
                      remove_degree6_N(v, f[2], f[3], f[4], f[5], f[0], f[1],
                                       w[2], w[3], w[4], w[5], w[0], w[1],
                                       i[2], i[3], i[4], i[5], i[0], i[1]);
                    }
                }
              else
                {
                  if(incircle(4, 5, 0, 2, f, w, i))
                    {
                      remove_degree6_diamond(v, f[0], f[1], f[2], f[3], f[4], f[5],
                                             w[0], w[1], w[2], w[3], w[4], w[5],
                                             i[0], i[1], i[2], i[3], i[4], i[5]);
                    }
                  else
                    {
                      if(incircle(5, 0, 1, 2, f, w, i))
                        {
                          remove_degree6_N(v, f[2], f[3], f[4], f[5], f[0], f[1],
                                           w[2], w[3], w[4], w[5], w[0], w[1],
                                           i[2], i[3], i[4], i[5], i[0], i[1]);
                        }
                      else
                        {
                          remove_degree6_star(v, f[2], f[3], f[4], f[5], f[0], f[1],
                                              w[2], w[3], w[4], w[5], w[0], w[1],
                                              i[2], i[3], i[4], i[5], i[0], i[1]);
                        }
                    }
                }
            }
          else
            {
              remove_degree6_star(v, f[0], f[1], f[2], f[3], f[4], f[5],
                                  w[0], w[1], w[2], w[3], w[4], w[5],
                                  i[0], i[1], i[2], i[3], i[4], i[5]);
            }
        }
      else
        {
          if(incircle(5, 2, 3, 0, f, w, i))
            {
              if(incircle(5, 0, 1, 2, f, w, i))
                {
                  remove_degree6_star(v, f[5], f[0], f[1], f[2], f[3], f[4],
                                      w[5], w[0], w[1], w[2], w[3], w[4],
                                      i[5], i[0], i[1], i[2], i[3], i[4]);
                }
              else
                {
                  remove_degree6_antiN(v, f[2], f[3], f[4], f[5], f[0], f[1],
                                       w[2], w[3], w[4], w[5], w[0], w[1],
                                       i[2], i[3], i[4], i[5], i[0], i[1]);
                }
            }
          else
            {
              if(incircle(4, 5, 0, 3, f, w, i))
                {
                  remove_degree6_star(v, f[0], f[1], f[2], f[3], f[4], f[5],
                                      w[0], w[1], w[2], w[3], w[4], w[5],
                                      i[0], i[1], i[2], i[3], i[4], i[5]);
                }
              else
                {
                  remove_degree6_N(v, f[0], f[1], f[2], f[3], f[4], f[5],
                                   w[0], w[1], w[2], w[3], w[4], w[5],
                                   i[0], i[1], i[2], i[3], i[4], i[5]);
                }
            }
        }
    }
}

template < class Gt, class Tds >
void
Periodic_2_Delaunay_triangulation_2<Gt, Tds>::
remove_degree6(Vertex_handle v, std::vector<Face_handle> &f,
               std::vector<Vertex_handle> &w, std::vector<Offset> &o,
               std::vector<int> &i)
{
  // removing a degree 6 vertex
  this->remove_too_long_edges_in_star(v);

  if(incircle(1, 2, 3, 0, f, w, o, i))
    {
      if(incircle(4, 2, 3, 5, f, w, o, i))
        {
          if(incircle(1, 2, 3, 4, f, w, o, i))
            {
              if(incircle(4, 0, 1, 3, f, w, o, i))
                {
                  if(incircle(5, 0, 1, 4, f, w, o, i))
                    {
                      remove_degree6_star(v,
                                          f[1], f[2], f[3], f[4], f[5], f[0],
                                          w[1], w[2], w[3], w[4], w[5], w[0],
                                          o[1], o[2], o[3], o[4], o[5], o[0],
                                          i[1], i[2], i[3], i[4], i[5], i[0]);
                    }
                  else
                    {
                      remove_degree6_N(v,
                                       f[1], f[2], f[3], f[4], f[5], f[0],
                                       w[1], w[2], w[3], w[4], w[5], w[0],
                                       o[1], o[2], o[3], o[4], o[5], o[0],
                                       i[1], i[2], i[3], i[4], i[5], i[0]);
                    }
                }
              else
                {
                  remove_degree6_antiN(v,
                                       f[0], f[1], f[2], f[3], f[4], f[5],
                                       w[0], w[1], w[2], w[3], w[4], w[5],
                                       o[0], o[1], o[2], o[3], o[4], o[5],
                                       i[0], i[1], i[2], i[3], i[4], i[5]);
                }
            }
          else
            {
              if(incircle(5, 1, 2, 4, f, w, o, i))
                {
                  remove_degree6_N(v, f[2], f[3], f[4], f[5], f[0], f[1],
                                   w[2], w[3], w[4], w[5], w[0], w[1],
                                   o[2], o[3], o[4], o[5], o[0], o[1],
                                   i[2], i[3], i[4], i[5], i[0], i[1]);
                }
              else
                {
                  if(incircle(5, 0, 1, 4, f, w, o, i))
                    {
                      remove_degree6_antiN(v,
                                           f[1], f[2], f[3], f[4], f[5], f[0],
                                           w[1], w[2], w[3], w[4], w[5], w[0],
                                           o[1], o[2], o[3], o[4], o[5], o[0],
                                           i[1], i[2], i[3], i[4], i[5], i[0]);
                    }
                  else
                    {
                      remove_degree6_star(v, f[4], f[5], f[0], f[1], f[2], f[3],
                                          w[4], w[5], w[0], w[1], w[2], w[3],
                                          o[4], o[5], o[0], o[1], o[2], o[3],
                                          i[4], i[5], i[0], i[1], i[2], i[3]);
                    }
                }
            }
        }
      else
        {
          if(incircle(1, 2, 3, 5, f, w, o, i))
            {
              if(incircle(1, 3, 4, 5, f, w, o, i))
                {
                  if(incircle(4, 0, 1, 3, f, w, o, i))
                    {
                      if(incircle(5, 0, 1, 4, f, w, o, i))
                        {
                          remove_degree6_star(v, f[1], f[2], f[3], f[4], f[5], f[0],
                                              w[1], w[2], w[3], w[4], w[5], w[0],
                                              o[1], o[2], o[3], o[4], o[5], o[0],
                                              i[1], i[2], i[3], i[4], i[5], i[0]);
                        }
                      else
                        {
                          remove_degree6_N(v, f[1], f[2], f[3], f[4], f[5], f[0],
                                           w[1], w[2], w[3], w[4], w[5], w[0],
                                           o[1], o[2], o[3], o[4], o[5], o[0],
                                           i[1], i[2], i[3], i[4], i[5], i[0]);
                        }
                    }
                  else
                    {
                      remove_degree6_antiN(v, f[0], f[1], f[2], f[3], f[4], f[5],
                                           w[0], w[1], w[2], w[3], w[4], w[5],
                                           o[0], o[1], o[2], o[3], o[4], o[5],
                                           i[0], i[1], i[2], i[3], i[4], i[5]);
                    }
                }
              else
                {
                  if(incircle(5, 0, 1, 3, f, w, o, i))
                    {
                      remove_degree6_diamond(v, f[1], f[2], f[3], f[4], f[5], f[0],
                                             w[1], w[2], w[3], w[4], w[5], w[0],
                                             o[1], o[2], o[3], o[4], o[5], o[0],
                                             i[1], i[2], i[3], i[4], i[5], i[0]);
                    }
                  else
                    {
                      if(incircle(4, 5, 0, 3, f, w, o, i))
                        {
                          remove_degree6_antiN(v, f[0], f[1], f[2], f[3], f[4], f[5],
                                               w[0], w[1], w[2], w[3], w[4], w[5],
                                               o[0], o[1], o[2], o[3], o[4], o[5],
                                               i[0], i[1], i[2], i[3], i[4], i[5]);
                        }
                      else
                        {
                          remove_degree6_star(v, f[3], f[4], f[5], f[0], f[1], f[2],
                                              w[3], w[4], w[5], w[0], w[1], w[2],
                                              o[3], o[4], o[5], o[0], o[1], o[2],
                                              i[3], i[4], i[5], i[0], i[1], i[2]);
                        }
                    }
                }
            }
          else
            {
              remove_degree6_star(v, f[5], f[0], f[1], f[2], f[3], f[4],
                                  w[5], w[0], w[1], w[2], w[3], w[4],
                                  o[5], o[0], o[1], o[2], o[3], o[4],
                                  i[5], i[0], i[1], i[2], i[3], i[4]);
            }
        }
    }
  else
    {
      if(incircle(4, 2, 3, 5, f, w, o, i))
        {
          if(incircle(4, 2, 3, 0, f, w, o, i))
            {
              if(incircle(4, 0, 1, 2, f, w, o, i))
                {
                  if(incircle(4, 1, 2, 5, f, w, o, i))
                    {
                      if(incircle(4, 0, 1, 5, f, w, o, i))
                        {
                          remove_degree6_star(v, f[4], f[5], f[0], f[1], f[2], f[3],
                                              w[4], w[5], w[0], w[1], w[2], w[3],
                                              o[4], o[5], o[0], o[1], o[2], o[3],
                                              i[4], i[5], i[0], i[1], i[2], i[3]);
                        }
                      else
                        {
                          remove_degree6_antiN(v, f[1], f[2], f[3], f[4], f[5], f[0],
                                               w[1], w[2], w[3], w[4], w[5], w[0],
                                               o[1], o[2], o[3], o[4], o[5], o[0],
                                               i[1], i[2], i[3], i[4], i[5], i[0]);
                        }
                    }
                  else
                    {
                      remove_degree6_N(v, f[2], f[3], f[4], f[5], f[0], f[1],
                                       w[2], w[3], w[4], w[5], w[0], w[1],
                                       o[2], o[3], o[4], o[5], o[0], o[1],
                                       i[2], i[3], i[4], i[5], i[0], i[1]);
                    }
                }
              else
                {
                  if(incircle(4, 5, 0, 2, f, w, o, i))
                    {
                      remove_degree6_diamond(v, f[0], f[1], f[2], f[3], f[4], f[5],
                                             w[0], w[1], w[2], w[3], w[4], w[5],
                                             o[0], o[1], o[2], o[3], o[4], o[5],
                                             i[0], i[1], i[2], i[3], i[4], i[5]);
                    }
                  else
                    {
                      if(incircle(5, 0, 1, 2, f, w, o, i))
                        {
                          remove_degree6_N(v, f[2], f[3], f[4], f[5], f[0], f[1],
                                           w[2], w[3], w[4], w[5], w[0], w[1],
                                           o[2], o[3], o[4], o[5], o[0], o[1],
                                           i[2], i[3], i[4], i[5], i[0], i[1]);
                        }
                      else
                        {
                          remove_degree6_star(v, f[2], f[3], f[4], f[5], f[0], f[1],
                                              w[2], w[3], w[4], w[5], w[0], w[1],
                                              o[2], o[3], o[4], o[5], o[0], o[1],
                                              i[2], i[3], i[4], i[5], i[0], i[1]);
                        }
                    }
                }
            }
          else
            {
              remove_degree6_star(v, f[0], f[1], f[2], f[3], f[4], f[5],
                                  w[0], w[1], w[2], w[3], w[4], w[5],
                                  o[0], o[1], o[2], o[3], o[4], o[5],
                                  i[0], i[1], i[2], i[3], i[4], i[5]);
            }
        }
      else
        {
          if(incircle(5, 2, 3, 0, f, w, o, i))
            {
              if(incircle(5, 0, 1, 2, f, w, o, i))
                {
                  remove_degree6_star(v, f[5], f[0], f[1], f[2], f[3], f[4],
                                      w[5], w[0], w[1], w[2], w[3], w[4],
                                      o[5], o[0], o[1], o[2], o[3], o[4],
                                      i[5], i[0], i[1], i[2], i[3], i[4]);
                }
              else
                {
                  remove_degree6_antiN(v, f[2], f[3], f[4], f[5], f[0], f[1],
                                       w[2], w[3], w[4], w[5], w[0], w[1],
                                       o[2], o[3], o[4], o[5], o[0], o[1],
                                       i[2], i[3], i[4], i[5], i[0], i[1]);
                }
            }
          else
            {
              if(incircle(4, 5, 0, 3, f, w, o, i))
                {
                  remove_degree6_star(v, f[0], f[1], f[2], f[3], f[4], f[5],
                                      w[0], w[1], w[2], w[3], w[4], w[5],
                                      o[0], o[1], o[2], o[3], o[4], o[5],
                                      i[0], i[1], i[2], i[3], i[4], i[5]);
                }
              else
                {
                  remove_degree6_N(v, f[0], f[1], f[2], f[3], f[4], f[5],
                                   w[0], w[1], w[2], w[3], w[4], w[5],
                                   o[0], o[1], o[2], o[3], o[4], o[5],
                                   i[0], i[1], i[2], i[3], i[4], i[5]);
                }
            }
        }
    }
}

template < class Gt, class Tds >
void
Periodic_2_Delaunay_triangulation_2<Gt, Tds>::remove_degree6_star
(Vertex_handle &,
 Face_handle &  f0, Face_handle &  f1, Face_handle &  f2,
 Face_handle &  f3, Face_handle &  f4, Face_handle &  f5,
 Vertex_handle &v0, Vertex_handle &, Vertex_handle &,
 Vertex_handle &, Vertex_handle &, Vertex_handle &,
 int i0, int i1, int i2, int i3, int i4, int i5 )
{
  // removing a degree 6 vertex, staring from v0
  Face_handle nn;
  f1->set_vertex( i1, v0) ;  // f1 = v1v2v0
  f2->set_vertex( i2, v0) ;  // f2 = v2v3v0
  f3->set_vertex( i3, v0) ;  // f3 = v3v4v0
  f4->set_vertex( i4, v0) ;  // f4 = v4v5v0
  nn = f0->neighbor( i0 );
  tds().set_adjacency(f1, cw(i1), nn, nn->index(f0));
  nn = f5->neighbor( i5 );
  tds().set_adjacency(f4, ccw(i4), nn,  nn->index(f5));
  tds().delete_face(f0);
  tds().delete_face(f5);

  f1->set_offsets(0, 0, 0);
  f2->set_offsets(0, 0, 0);
  f3->set_offsets(0, 0, 0);
  f4->set_offsets(0, 0, 0);
}

template < class Gt, class Tds >
void
Periodic_2_Delaunay_triangulation_2<Gt, Tds>::remove_degree6_star
(Vertex_handle &v,
 Face_handle &  f0, Face_handle &  f1, Face_handle &  f2,
 Face_handle &  f3, Face_handle &  f4, Face_handle &  f5,
 Vertex_handle &v0, Vertex_handle &v1, Vertex_handle &v2,
 Vertex_handle &v3, Vertex_handle &v4, Vertex_handle &v5,
 Offset &o0, Offset &o1, Offset &o2,
 Offset &o3, Offset &o4, Offset &o5,
 int i0, int i1, int i2, int i3, int i4, int i5 )
{
  // removing a degree 6 vertex, staring from v0
  remove_degree6_star(v,
                      f0, f1, f2, f3, f4, f5,
                      v0, v1, v2, v3, v4, v5,
                      i0, i1, i2, i3, i4, i5);

  if (o0.x() < 0 || o1.x() < 0 || o2.x() < 0 || o3.x() < 0 || o4.x() < 0 || o5.x() < 0)
    {
      o0 += Offset(number_of_sheets()[0], 0);
      o1 += Offset(number_of_sheets()[0], 0);
      o2 += Offset(number_of_sheets()[0], 0);
      o3 += Offset(number_of_sheets()[0], 0);
      o4 += Offset(number_of_sheets()[0], 0);
      o5 += Offset(number_of_sheets()[0], 0);
    }
  if (o0.y() < 0 || o1.y() < 0 || o2.y() < 0 || o3.y() < 0 || o4.y() < 0 || o5.y() < 0)
    {
      o0 += Offset(0, number_of_sheets()[1]);
      o1 += Offset(0, number_of_sheets()[1]);
      o2 += Offset(0, number_of_sheets()[1]);
      o3 += Offset(0, number_of_sheets()[1]);
      o4 += Offset(0, number_of_sheets()[1]);
      o5 += Offset(0, number_of_sheets()[1]);
    }
  int oo0 = (o0.x() >= number_of_sheets()[0] ? 2 : 0) + (o0.y() >= number_of_sheets()[1] ? 1 : 0);
  int oo1 = (o1.x() >= number_of_sheets()[0] ? 2 : 0) + (o1.y() >= number_of_sheets()[1] ? 1 : 0);
  int oo2 = (o2.x() >= number_of_sheets()[0] ? 2 : 0) + (o2.y() >= number_of_sheets()[1] ? 1 : 0);
  int oo3 = (o3.x() >= number_of_sheets()[0] ? 2 : 0) + (o3.y() >= number_of_sheets()[1] ? 1 : 0);
  int oo4 = (o4.x() >= number_of_sheets()[0] ? 2 : 0) + (o4.y() >= number_of_sheets()[1] ? 1 : 0);
  int oo5 = (o5.x() >= number_of_sheets()[0] ? 2 : 0) + (o5.y() >= number_of_sheets()[1] ? 1 : 0);

  int oo[3];
  oo[i1]      = oo0;
  oo[ccw(i1)] = oo1;
  oo[ cw(i1)] = oo2;
  this->set_offsets(f1, oo[0], oo[1], oo[2]);
  oo[i2]      = oo0;
  oo[ccw(i2)] = oo2;
  oo[ cw(i2)] = oo3;
  this->set_offsets(f2, oo[0], oo[1], oo[2]);
  oo[i3]      = oo0;
  oo[ccw(i3)] = oo3;
  oo[ cw(i3)] = oo4;
  this->set_offsets(f3, oo[0], oo[1], oo[2]);
  oo[i4]      = oo0;
  oo[ccw(i4)] = oo4;
  oo[ cw(i4)] = oo5;
  this->set_offsets(f4, oo[0], oo[1], oo[2]);

  insert_too_long_edge(f1, ccw(i1));
  insert_too_long_edge(f2, ccw(i2));
  insert_too_long_edge(f3, ccw(i3));
}

template < class Gt, class Tds >
void
Periodic_2_Delaunay_triangulation_2<Gt, Tds>::remove_degree6_N
(
  Vertex_handle &,
  Face_handle &  f0, Face_handle &  f1, Face_handle &  f2,
  Face_handle &  f3, Face_handle &  f4, Face_handle &  f5,
  Vertex_handle &v0, Vertex_handle &, Vertex_handle &,
  Vertex_handle &v3, Vertex_handle &, Vertex_handle &,
  int i0, int i1, int i2, int i3, int i4, int i5 )
{
  // removing a degree 6 vertex, N configuration with diagonal v0v3
  Face_handle nn;
  f1->set_vertex( i1, v0) ;  // f1 = v1v2v0
  f2->set_vertex( i2, v0) ;  // f2 = v2v3v0
  f4->set_vertex( i4, v3) ;  // f4 = v4v5v3
  f5->set_vertex( i5, v3) ;  // f5 = v5v0v3
  nn = f0->neighbor( i0 );
  tds().set_adjacency(f1, cw(i1) , nn , nn->index(f0)  );
  nn = f3->neighbor( i3 );
  tds().set_adjacency(f4, cw(i4) , nn, nn->index(f3) );
  tds().set_adjacency(f2, ccw(i2) , f5 , ccw(i5)  );
  tds().delete_face(f0);
  tds().delete_face(f3);

  f1->set_offsets(0, 0, 0);
  f2->set_offsets(0, 0, 0);
  f4->set_offsets(0, 0, 0);
  f5->set_offsets(0, 0, 0);
}
template < class Gt, class Tds >
void
Periodic_2_Delaunay_triangulation_2<Gt, Tds>::remove_degree6_N
(
  Vertex_handle &v,
  Face_handle &  f0, Face_handle &  f1, Face_handle &  f2,
  Face_handle &  f3, Face_handle &  f4, Face_handle &  f5,
  Vertex_handle &v0, Vertex_handle &v1, Vertex_handle &v2,
  Vertex_handle &v3, Vertex_handle &v4, Vertex_handle &v5,
  Offset &o0, Offset &o1, Offset &o2,
  Offset &o3, Offset &o4, Offset &o5,
  int i0, int i1, int i2, int i3, int i4, int i5 )
{
  // removing a degree 6 vertex, N configuration with diagonal v0v3

  remove_degree6_N(v,
                   f0, f1, f2, f3, f4, f5,
                   v0, v1, v2, v3, v4, v5,
                   i0, i1, i2, i3, i4, i5);

  if (o0.x() < 0 || o1.x() < 0 || o2.x() < 0 || o3.x() < 0 || o4.x() < 0 || o5.x() < 0)
    {
      o0 += Offset(number_of_sheets()[0], 0);
      o1 += Offset(number_of_sheets()[0], 0);
      o2 += Offset(number_of_sheets()[0], 0);
      o3 += Offset(number_of_sheets()[0], 0);
      o4 += Offset(number_of_sheets()[0], 0);
      o5 += Offset(number_of_sheets()[0], 0);
    }
  if (o0.y() < 0 || o1.y() < 0 || o2.y() < 0 || o3.y() < 0 || o4.y() < 0 || o5.y() < 0)
    {
      o0 += Offset(0, number_of_sheets()[1]);
      o1 += Offset(0, number_of_sheets()[1]);
      o2 += Offset(0, number_of_sheets()[1]);
      o3 += Offset(0, number_of_sheets()[1]);
      o4 += Offset(0, number_of_sheets()[1]);
      o5 += Offset(0, number_of_sheets()[1]);
    }
  int oo0 = (o0.x() >= number_of_sheets()[0] ? 2 : 0) + (o0.y() >= number_of_sheets()[1] ? 1 : 0);
  int oo1 = (o1.x() >= number_of_sheets()[0] ? 2 : 0) + (o1.y() >= number_of_sheets()[1] ? 1 : 0);
  int oo2 = (o2.x() >= number_of_sheets()[0] ? 2 : 0) + (o2.y() >= number_of_sheets()[1] ? 1 : 0);
  int oo3 = (o3.x() >= number_of_sheets()[0] ? 2 : 0) + (o3.y() >= number_of_sheets()[1] ? 1 : 0);
  int oo4 = (o4.x() >= number_of_sheets()[0] ? 2 : 0) + (o4.y() >= number_of_sheets()[1] ? 1 : 0);
  int oo5 = (o5.x() >= number_of_sheets()[0] ? 2 : 0) + (o5.y() >= number_of_sheets()[1] ? 1 : 0);

  int oo[3];
  oo[i1]      = oo0;
  oo[ccw(i1)] = oo1;
  oo[ cw(i1)] = oo2;
  this->set_offsets(f1, oo[0], oo[1], oo[2]);
  oo[i2]      = oo0;
  oo[ccw(i2)] = oo2;
  oo[ cw(i2)] = oo3;
  this->set_offsets(f2, oo[0], oo[1], oo[2]);
  oo[i4]      = oo3;
  oo[ccw(i4)] = oo4;
  oo[ cw(i4)] = oo5;
  this->set_offsets(f4, oo[0], oo[1], oo[2]);
  oo[i5]      = oo3;
  oo[ccw(i5)] = oo5;
  oo[ cw(i5)] = oo0;
  this->set_offsets(f5, oo[0], oo[1], oo[2]);

  insert_too_long_edge(f1, ccw(i1));
  insert_too_long_edge(f2, ccw(i2));
  insert_too_long_edge(f4, ccw(i4));
  insert_too_long_edge(f5, ccw(i5));
}


template < class Gt, class Tds >
void
Periodic_2_Delaunay_triangulation_2<Gt, Tds>::remove_degree6_antiN
(
  Vertex_handle &,
  Face_handle &  f0, Face_handle &  f1, Face_handle &  f2,
  Face_handle &  f3, Face_handle &  f4, Face_handle &  f5,
  Vertex_handle &v0, Vertex_handle &, Vertex_handle &,
  Vertex_handle &v3, Vertex_handle &, Vertex_handle &,
  int i0, int i1, int i2, int i3, int i4, int i5 )
{
  // removing a degree 6 vertex, antiN configuration with diagonal v0v3
  Face_handle nn;
  f0->set_vertex( i0, v3) ;  // f0 = v0v1v3
  f1->set_vertex( i1, v3) ;  // f1 = v1v2v3
  f3->set_vertex( i3, v0) ;  // f3 = v3v4v0
  f4->set_vertex( i4, v0) ;  // f4 = v4v5v0
  nn = f2->neighbor( i2 );
  tds().set_adjacency(f1, ccw(i1) , nn , nn->index(f2)  );
  nn = f5->neighbor( i5 );
  tds().set_adjacency(f4, ccw(i4) , nn , nn->index(f5) );
  tds().set_adjacency(f0, cw(i0) , f3, cw(i3) );
  tds().delete_face(f2);
  tds().delete_face(f5);

  f0->set_offsets(0, 0, 0);
  f1->set_offsets(0, 0, 0);
  f3->set_offsets(0, 0, 0);
  f4->set_offsets(0, 0, 0);
}
template < class Gt, class Tds >
void
Periodic_2_Delaunay_triangulation_2<Gt, Tds>::remove_degree6_antiN
(
  Vertex_handle &v,
  Face_handle &  f0, Face_handle &  f1, Face_handle &  f2,
  Face_handle &  f3, Face_handle &  f4, Face_handle &  f5,
  Vertex_handle &v0, Vertex_handle &v1, Vertex_handle &v2,
  Vertex_handle &v3, Vertex_handle &v4, Vertex_handle &v5,
  Offset &o0, Offset &o1, Offset &o2,
  Offset &o3, Offset &o4, Offset &o5,
  int i0, int i1, int i2, int i3, int i4, int i5 )
{
  // removing a degree 6 vertex, antiN configuration with diagonal v0v3
  remove_degree6_antiN(v,
                       f0, f1, f2, f3, f4, f5,
                       v0, v1, v2, v3, v4, v5,
                       i0, i1, i2, i3, i4, i5);

  if (o0.x() < 0 || o1.x() < 0 || o2.x() < 0 || o3.x() < 0 || o4.x() < 0 || o5.x() < 0)
    {
      o0 += Offset(number_of_sheets()[0], 0);
      o1 += Offset(number_of_sheets()[0], 0);
      o2 += Offset(number_of_sheets()[0], 0);
      o3 += Offset(number_of_sheets()[0], 0);
      o4 += Offset(number_of_sheets()[0], 0);
      o5 += Offset(number_of_sheets()[0], 0);
    }
  if (o0.y() < 0 || o1.y() < 0 || o2.y() < 0 || o3.y() < 0 || o4.y() < 0 || o5.y() < 0)
    {
      o0 += Offset(0, number_of_sheets()[1]);
      o1 += Offset(0, number_of_sheets()[1]);
      o2 += Offset(0, number_of_sheets()[1]);
      o3 += Offset(0, number_of_sheets()[1]);
      o4 += Offset(0, number_of_sheets()[1]);
      o5 += Offset(0, number_of_sheets()[1]);
    }
  int oo0 = (o0.x() >= number_of_sheets()[0] ? 2 : 0) + (o0.y() >= number_of_sheets()[1] ? 1 : 0);
  int oo1 = (o1.x() >= number_of_sheets()[0] ? 2 : 0) + (o1.y() >= number_of_sheets()[1] ? 1 : 0);
  int oo2 = (o2.x() >= number_of_sheets()[0] ? 2 : 0) + (o2.y() >= number_of_sheets()[1] ? 1 : 0);
  int oo3 = (o3.x() >= number_of_sheets()[0] ? 2 : 0) + (o3.y() >= number_of_sheets()[1] ? 1 : 0);
  int oo4 = (o4.x() >= number_of_sheets()[0] ? 2 : 0) + (o4.y() >= number_of_sheets()[1] ? 1 : 0);
  int oo5 = (o5.x() >= number_of_sheets()[0] ? 2 : 0) + (o5.y() >= number_of_sheets()[1] ? 1 : 0);

  int oo[3];
  oo[i0]      = oo3;
  oo[ccw(i0)] = oo0;
  oo[ cw(i0)] = oo1;
  this->set_offsets(f0, oo[0], oo[1], oo[2]);
  oo[i1]      = oo3;
  oo[ccw(i1)] = oo1;
  oo[ cw(i1)] = oo2;
  this->set_offsets(f1, oo[0], oo[1], oo[2]);
  oo[i3]      = oo0;
  oo[ccw(i3)] = oo3;
  oo[ cw(i3)] = oo4;
  this->set_offsets(f3, oo[0], oo[1], oo[2]);
  oo[i4]      = oo0;
  oo[ccw(i4)] = oo4;
  oo[ cw(i4)] = oo5;
  this->set_offsets(f4, oo[0], oo[1], oo[2]);

  insert_too_long_edge(f0, ccw(i0));
  insert_too_long_edge(f0,  cw(i0));
  insert_too_long_edge(f3, ccw(i3));
}

template < class Gt, class Tds >
void
Periodic_2_Delaunay_triangulation_2<Gt, Tds>::remove_degree6_diamond
(
  Vertex_handle &,
  Face_handle &  f0, Face_handle &  f1, Face_handle &  f2,
  Face_handle &  f3, Face_handle &  f4, Face_handle &  f5,
  Vertex_handle &v0, Vertex_handle &, Vertex_handle &v2,
  Vertex_handle &, Vertex_handle &v4, Vertex_handle &,
  int i0, int i1, int i2, int i3, int i4, int i5 )
{
  // removing a degree 6 vertex, with chords v0v2 v2v4 v4v0
  Face_handle nn;
  f0->set_vertex( i0, v2) ;  // f0 = v0v1v2
  f2->set_vertex( i2, v4) ;  // f2 = v2v3v4
  f4->set_vertex( i4, v0) ;  // f4 = v4v5v0
  f1->set_vertex( i1, v4) ;
  f1->set_vertex( ccw(i1), v0) ;  // f1 = v0v2v4
  nn = f1->neighbor( i1 );
  tds().set_adjacency(f0, ccw(i0) , nn , nn->index(f1) );
  nn = f3->neighbor( i3 );
  tds().set_adjacency(f2, ccw(i2) , nn , nn->index(f3) );
  nn = f5->neighbor( i5 );
  tds().set_adjacency(f4, ccw(i4) , nn , nn->index(f5) );
  tds().set_adjacency(f0, cw(i0) , f1 , i1  );
  tds().set_adjacency(f4, cw(i4) , f1 , cw(i1) );

  tds().delete_face(f3);
  tds().delete_face(f5);


  f0->set_offsets(0, 0, 0);
  f1->set_offsets(0, 0, 0);
  f2->set_offsets(0, 0, 0);
  f4->set_offsets(0, 0, 0);
}
template < class Gt, class Tds >
void
Periodic_2_Delaunay_triangulation_2<Gt, Tds>::remove_degree6_diamond
(
  Vertex_handle &v,
  Face_handle &  f0, Face_handle &  f1, Face_handle &  f2,
  Face_handle &  f3, Face_handle &  f4, Face_handle &  f5,
  Vertex_handle &v0, Vertex_handle &v1, Vertex_handle &v2,
  Vertex_handle &v3, Vertex_handle &v4, Vertex_handle &v5,
  Offset &o0, Offset &o1, Offset &o2,
  Offset &o3, Offset &o4, Offset &o5,
  int i0, int i1, int i2, int i3, int i4, int i5 )
{
  // removing a degree 6 vertex, with chords v0v2 v2v4 v4v0
  remove_degree6_diamond(v,
                         f0, f1, f2, f3, f4, f5,
                         v0, v1, v2, v3, v4, v5,
                         i0, i1, i2, i3, i4, i5);

  if (o0.x() < 0 || o1.x() < 0 || o2.x() < 0 || o3.x() < 0 || o4.x() < 0 || o5.x() < 0)
    {
      o0 += Offset(number_of_sheets()[0], 0);
      o1 += Offset(number_of_sheets()[0], 0);
      o2 += Offset(number_of_sheets()[0], 0);
      o3 += Offset(number_of_sheets()[0], 0);
      o4 += Offset(number_of_sheets()[0], 0);
      o5 += Offset(number_of_sheets()[0], 0);
    }
  if (o0.y() < 0 || o1.y() < 0 || o2.y() < 0 || o3.y() < 0 || o4.y() < 0 || o5.y() < 0)
    {
      o0 += Offset(0, number_of_sheets()[1]);
      o1 += Offset(0, number_of_sheets()[1]);
      o2 += Offset(0, number_of_sheets()[1]);
      o3 += Offset(0, number_of_sheets()[1]);
      o4 += Offset(0, number_of_sheets()[1]);
      o5 += Offset(0, number_of_sheets()[1]);
    }
  int oo0 = (o0.x() >= number_of_sheets()[0] ? 2 : 0) + (o0.y() >= number_of_sheets()[1] ? 1 : 0);
  int oo1 = (o1.x() >= number_of_sheets()[0] ? 2 : 0) + (o1.y() >= number_of_sheets()[1] ? 1 : 0);
  int oo2 = (o2.x() >= number_of_sheets()[0] ? 2 : 0) + (o2.y() >= number_of_sheets()[1] ? 1 : 0);
  int oo3 = (o3.x() >= number_of_sheets()[0] ? 2 : 0) + (o3.y() >= number_of_sheets()[1] ? 1 : 0);
  int oo4 = (o4.x() >= number_of_sheets()[0] ? 2 : 0) + (o4.y() >= number_of_sheets()[1] ? 1 : 0);
  int oo5 = (o5.x() >= number_of_sheets()[0] ? 2 : 0) + (o5.y() >= number_of_sheets()[1] ? 1 : 0);

  int oo[3];
  oo[i0]      = oo2;
  oo[ccw(i0)] = oo0;
  oo[ cw(i0)] = oo1;
  this->set_offsets(f0, oo[0], oo[1], oo[2]);
  oo[i2]      = oo4;
  oo[ccw(i2)] = oo2;
  oo[ cw(i2)] = oo3;
  this->set_offsets(f2, oo[0], oo[1], oo[2]);
  oo[i4]      = oo0;
  oo[ccw(i4)] = oo4;
  oo[ cw(i4)] = oo5;
  this->set_offsets(f4, oo[0], oo[1], oo[2]);
  oo[i1]      = oo4;
  oo[ccw(i1)] = oo0;
  oo[ cw(i1)] = oo2;
  this->set_offsets(f1, oo[0], oo[1], oo[2]);

  insert_too_long_edge(f1, 0);
  insert_too_long_edge(f1, 1);
  insert_too_long_edge(f1, 2);
}


template < class Gt, class Tds >
void
Periodic_2_Delaunay_triangulation_2<Gt, Tds>::
remove_degree7(Vertex_handle v, std::vector<Face_handle> &f,
               std::vector<Vertex_handle> &w,
               std::vector<int> &i)
{
  // removing a degree 7 vertex

  if (incircle(2, 0, 1, 3, f, w, i)) // sweeping from above
    {
      if (incircle(2, 3, 4, 0, f, w, i))
        {
          if (incircle(5, 3, 4, 6, f, w, i))
            {
              if (incircle(5, 3, 4, 2, f, w, i))
                {
                  if (incircle(6, 2, 3, 5, f, w, i))
                    {
                      if (incircle(6, 0, 1, 2, f, w, i))
                        {
                          remove_degree7_leftfan(v,  6  , f, w, i);
                        }
                      else
                        {
                          remove_degree7_zigzag(v,  6  , f, w, i);
                        }
                    }
                  else
                    {
                      if (incircle(5, 0, 1, 2, f, w, i))
                        {
                          if (incircle(6, 1, 2, 5, f, w, i))
                            {
                              remove_degree7_zigzag(v, 2   , f, w, i);
                            }
                          else
                            {
                              if (incircle(6, 0, 1, 5, f, w, i))
                                {
                                  remove_degree7_rightfan(v, 5   , f, w, i);
                                }
                              else
                                {
                                  remove_degree7_star(v,  5  , f, w, i);
                                }
                            }
                        }
                      else
                        {
                          if (incircle(2, 5, 6, 0, f, w, i))
                            {
                              if (incircle(6, 0, 1, 2, f, w, i))
                                {
                                  remove_degree7_zigzag(v,  2  , f, w, i);
                                }
                              else
                                {
                                  remove_degree7_rightfan(v,  2  , f, w, i);
                                }
                            }
                          else
                            {
                              remove_degree7_rightdelta(v,  5  , f, w, i);
                            }
                        }
                    }
                }
              else
                {
                  if (incircle(4, 0, 1, 2, f, w, i))
                    {
                      if (incircle(5, 1, 2, 4, f, w, i))
                        {
                          if (incircle(6, 1, 2, 5, f, w, i))
                            {
                              remove_degree7_leftfan(v,  2  , f, w, i);
                            }
                          else
                            {
                              if (incircle(6, 0, 1, 5, f, w, i))
                                {
                                  remove_degree7_zigzag(v,  5  , f, w, i);
                                }
                              else
                                {
                                  remove_degree7_leftfan(v,  5  , f, w, i);
                                }
                            }
                        }
                      else
                        {
                          if (incircle(5, 0, 1, 4, f, w, i))
                            {
                              if (incircle(6, 0, 1, 5, f, w, i))
                                {
                                  remove_degree7_rightfan(v,  1  , f, w, i);
                                }
                              else
                                {
                                  remove_degree7_zigzag(v,  1  , f, w, i);
                                }
                            }
                          else
                            {
                              remove_degree7_rightfan(v,  4  , f, w, i);
                            }
                        }
                    }
                  else
                    {
                      if (incircle(2, 4, 5, 0, f, w, i))
                        {
                          if (incircle(5, 0, 1, 2, f, w, i))
                            {
                              if (incircle(6, 1, 2, 5, f, w, i))
                                {
                                  remove_degree7_leftfan(v,  2  , f, w, i);
                                }
                              else
                                {
                                  if (incircle(6, 0, 1, 5, f, w, i))
                                    {
                                      remove_degree7_zigzag(v,  5  , f, w, i);
                                    }
                                  else
                                    {
                                      remove_degree7_leftfan(v,  5  , f, w, i);
                                    }
                                }
                            }
                          else
                            {
                              if (incircle(2, 5, 6, 0, f, w, i))
                                {
                                  if (incircle(6, 0, 1, 2, f, w, i))
                                    {
                                      remove_degree7_leftfan(v,  2  , f, w, i);
                                    }
                                  else
                                    {
                                      remove_degree7_star(v,  2  , f, w, i);
                                    }
                                }
                              else
                                {
                                  remove_degree7_leftdelta(v,  2  , f, w, i);
                                }
                            }
                        }
                      else
                        {
                          remove_degree7_rightdelta(v,  0  , f, w, i);
                        }
                    }
                }
            }
          else
            {
              if (incircle(6, 3, 4, 2, f, w, i))
                {
                  if (incircle(6, 0, 1, 2, f, w, i))
                    {
                      remove_degree7_star(v,  6  , f, w, i);
                    }
                  else
                    {
                      remove_degree7_rightfan(v,  6  , f, w, i);
                    }
                }
              else
                {
                  if (incircle(4, 0, 1, 2, f, w, i))
                    {
                      if (incircle(2, 4, 5, 6, f, w, i))
                        {
                          if (incircle(5, 1, 2, 4, f, w, i))
                            {
                              if (incircle(6, 1, 2, 5, f, w, i))
                                {
                                  remove_degree7_leftfan(v,  2  , f, w, i);
                                }
                              else
                                {
                                  if (incircle(6, 0, 1, 5, f, w, i))
                                    {
                                      remove_degree7_zigzag(v,  5  , f, w, i);
                                    }
                                  else
                                    {
                                      remove_degree7_leftfan(v,  5  , f, w, i);
                                    }
                                }
                            }
                          else
                            {
                              if (incircle(5, 0, 1, 4, f, w, i))
                                {
                                  if (incircle(6, 0, 1, 5, f, w, i))
                                    {
                                      remove_degree7_rightfan(v,  1  , f, w, i);
                                    }
                                  else
                                    {
                                      remove_degree7_zigzag(v,  1  , f, w, i);
                                    }
                                }
                              else
                                {
                                  remove_degree7_rightfan(v,  4  , f, w, i);
                                }
                            }
                        }
                      else
                        {
                          if (incircle(6, 1, 2, 4, f, w, i))
                            {
                              remove_degree7_leftdelta(v,  6  , f, w, i);
                            }
                          else
                            {
                              if (incircle(1, 4, 5, 6, f, w, i))
                                {
                                  if (incircle(1, 4, 5, 0, f, w, i))
                                    {
                                      if (incircle(6, 0, 1, 5, f, w, i))
                                        {
                                          remove_degree7_rightfan(v,  1  , f, w, i);
                                        }
                                      else
                                        {
                                          remove_degree7_zigzag(v,  1  , f, w, i);
                                        }
                                    }
                                  else
                                    {
                                      remove_degree7_rightfan(v,  4  , f, w, i);
                                    }
                                }
                              else
                                {
                                  if (incircle(6, 0, 1, 4, f, w, i))
                                    {
                                      remove_degree7_rightdelta(v,  4  , f, w, i);
                                    }
                                  else
                                    {
                                      if (incircle(6, 4, 5, 0, f, w, i))
                                        {
                                          remove_degree7_star(v,  4  , f, w, i);
                                        }
                                      else
                                        {
                                          remove_degree7_rightfan(v,  4  , f, w, i);
                                        }
                                    }
                                }
                            }
                        }
                    }
                  else
                    {
                      if (incircle(2, 4, 5, 6, f, w, i))
                        {
                          if (incircle(2, 4, 5, 0, f, w, i))
                            {
                              if (incircle(5, 0, 1, 2, f, w, i))
                                {
                                  if (incircle(6, 1, 2, 5, f, w, i))
                                    {
                                      remove_degree7_leftfan(v,  2  , f, w, i);
                                    }
                                  else
                                    {
                                      if (incircle(6, 0, 1, 5, f, w, i))
                                        {
                                          remove_degree7_zigzag(v,  5  , f, w, i);
                                        }
                                      else
                                        {
                                          remove_degree7_leftfan(v,  5  , f, w, i);
                                        }
                                    }
                                }
                              else
                                {
                                  if (incircle(2, 5, 6, 0, f, w, i))
                                    {
                                      if (incircle(6, 0, 1, 2, f, w, i))
                                        {
                                          remove_degree7_leftfan(v,  2  , f, w, i);
                                        }
                                      else
                                        {
                                          remove_degree7_star(v,  2  , f, w, i);
                                        }
                                    }
                                  else
                                    {
                                      remove_degree7_leftdelta(v,  2  , f, w, i);
                                    }
                                }
                            }
                          else
                            {
                              remove_degree7_rightdelta(v,  0  , f, w, i);
                            }
                        }
                      else
                        {
                          if (incircle(2, 6, 0, 4, f, w, i))
                            {
                              if (incircle(6, 0, 1, 2, f, w, i))
                                {
                                  remove_degree7_leftdelta(v,  6  , f, w, i);
                                }
                              else
                                {
                                  remove_degree7_rightdelta(v,  2  , f, w, i);
                                }
                            }
                          else
                            {
                              if (incircle(6, 4, 5, 0, f, w, i))
                                {
                                  remove_degree7_leftdelta(v,  4  , f, w, i);
                                }
                              else
                                {
                                  remove_degree7_rightdelta(v,  0  , f, w, i);
                                }
                            }
                        }
                    }
                }
            }
        }
      else
        {
          if (incircle(5, 3, 4, 6, f, w, i))
            {
              if (incircle(5, 3, 4, 0, f, w, i))
                {
                  if (incircle(5, 2, 3, 0, f, w, i))
                    {
                      if (incircle(6, 2, 3, 5, f, w, i))
                        {
                          if (incircle(6, 0, 1, 2, f, w, i))
                            {
                              remove_degree7_leftfan(v,  6  , f, w, i);
                            }
                          else
                            {
                              remove_degree7_zigzag(v,  6  , f, w, i);
                            }
                        }
                      else if (incircle(5, 0, 1, 2, f, w, i))
                        {
                          if (incircle(6, 1, 2, 5, f, w, i))
                            {
                              remove_degree7_zigzag(v,  2  , f, w, i);
                            }
                          else
                            {
                              if (incircle(6, 0, 1, 5, f, w, i))
                                {
                                  remove_degree7_rightfan(v,  5  , f, w, i);
                                }
                              else
                                {
                                  remove_degree7_star(v,  5  , f, w, i);
                                }
                            }
                        }
                      else
                        {
                          if (incircle(2, 5, 6, 0, f, w, i))
                            {
                              if (incircle(6, 0, 1, 2, f, w, i))
                                {
                                  remove_degree7_zigzag(v,  2  , f, w, i);
                                }
                              else
                                {
                                  remove_degree7_rightfan(v,  2  , f, w, i);
                                }
                            }
                          else
                            {
                              remove_degree7_rightdelta(v,  5  , f, w, i);
                            }
                        }
                    }
                  else
                    {
                      if (incircle(3, 5, 6, 0, f, w, i))
                        {
                          if (incircle(6, 2, 3, 0, f, w, i))
                            {
                              if (incircle(6, 0, 1, 2, f, w, i))
                                {
                                  remove_degree7_leftfan(v,  6  , f, w, i);
                                }
                              else
                                {
                                  remove_degree7_zigzag(v,  6  , f, w, i);
                                }
                            }
                          else
                            {
                              remove_degree7_leftfan(v,  3  , f, w, i);
                            }
                        }
                      else
                        {
                          remove_degree7_leftdelta(v,  0  , f, w, i);
                        }
                    }
                }
              else
                {
                  remove_degree7_star(v,  0  , f, w, i);
                }
            }
          else
            {
              if (incircle(6, 3, 4, 0, f, w, i))
                {
                  if (incircle(6, 2, 3, 0, f, w, i))
                    {
                      if (incircle(6, 0, 1, 2, f, w, i))
                        {
                          remove_degree7_star(v,  6  , f, w, i);
                        }
                      else
                        {
                          remove_degree7_rightfan(v,  6  , f, w, i);
                        }
                    }
                  else
                    {
                      remove_degree7_zigzag(v,  3  , f, w, i);
                    }
                }
              else
                {
                  if (incircle(6, 4, 5, 0, f, w, i))
                    {
                      remove_degree7_leftfan(v,  0  , f, w, i);
                    }
                  else
                    {
                      remove_degree7_star(v,  0  , f, w, i);
                    }
                }
            }
        }
    }
  else    //sweeping from below
    {
      if (incircle(1, 6, 0, 3, f, w, i))
        {
          if (incircle(5, 6, 0, 4, f, w, i))
            {
              if (incircle(5, 6, 0, 1, f, w, i))
                {
                  if (incircle(4, 0, 1, 5, f, w, i))
                    {
                      if (incircle(4, 2, 3, 1, f, w, i))
                        {
                          remove_degree7_rightfan(v,  4  , f, w, i);
                        }
                      else
                        {
                          remove_degree7_zigzag(v,  4  , f, w, i);
                        }
                    }
                  else
                    {
                      if (incircle(5, 2, 3, 1, f, w, i))
                        {
                          if (incircle(4, 1, 2, 5, f, w, i))
                            {
                              remove_degree7_zigzag(v, 1   , f, w, i);
                            }
                          else
                            {
                              if (incircle(4, 2, 3, 5, f, w, i))
                                {
                                  remove_degree7_leftfan(v, 5   , f, w, i);
                                }
                              else
                                {
                                  remove_degree7_star(v,  5  , f, w, i);
                                }
                            }
                        }
                      else
                        {
                          if (incircle(1, 4, 5, 3, f, w, i))
                            {
                              if (incircle(4, 2, 3, 1, f, w, i))
                                {
                                  remove_degree7_zigzag(v,  1  , f, w, i);
                                }
                              else
                                {
                                  remove_degree7_leftfan(v,  1  , f, w, i);
                                }
                            }
                          else
                            {
                              remove_degree7_leftdelta(v,  5  , f, w, i);
                            }
                        }
                    }
                }
              else
                {
                  if (incircle(6, 2, 3, 1, f, w, i))
                    {
                      if (incircle(5, 1, 2, 6, f, w, i))
                        {
                          if (incircle(4, 1, 2, 5, f, w, i))
                            {
                              remove_degree7_rightfan(v,  1  , f, w, i);
                            }
                          else
                            {
                              if (incircle(4, 2, 3, 5, f, w, i))
                                {
                                  remove_degree7_zigzag(v,  5  , f, w, i);
                                }
                              else
                                {
                                  remove_degree7_rightfan(v,  5  , f, w, i);
                                }
                            }
                        }
                      else
                        {
                          if (incircle(5, 2, 3, 6, f, w, i))
                            {
                              if (incircle(4, 2, 3, 5, f, w, i))
                                {
                                  remove_degree7_leftfan(v,  2  , f, w, i);
                                }
                              else
                                {
                                  remove_degree7_zigzag(v,  2  , f, w, i);
                                }
                            }
                          else
                            {
                              remove_degree7_leftfan(v,  6  , f, w, i);
                            }
                        }
                    }
                  else
                    {
                      if (incircle(1, 5, 6, 3, f, w, i))
                        {
                          if (incircle(5, 2, 3, 1, f, w, i))
                            {
                              if (incircle(4, 1, 2, 5, f, w, i))
                                {
                                  remove_degree7_rightfan(v,  1  , f, w, i);
                                }
                              else
                                {
                                  if (incircle(4, 2, 3, 5, f, w, i))
                                    {
                                      remove_degree7_zigzag(v,  5  , f, w, i);
                                    }
                                  else
                                    {
                                      remove_degree7_rightfan(v,  5  , f, w, i);
                                    }
                                }
                            }
                          else
                            {
                              if (incircle(1, 4, 5, 3, f, w, i))
                                {
                                  if (incircle(4, 2, 3, 1, f, w, i))
                                    {
                                      remove_degree7_rightfan(v,  1  , f, w, i);
                                    }
                                  else
                                    {
                                      remove_degree7_star(v,  1  , f, w, i);
                                    }
                                }
                              else
                                {
                                  remove_degree7_rightdelta(v,  1  , f, w, i);
                                }
                            }
                        }
                      else
                        {
                          remove_degree7_leftdelta(v,  3  , f, w, i);
                        }
                    }
                }
            }
          else
            {
              if (incircle(4, 6, 0, 1, f, w, i))
                {
                  if (incircle(4, 2, 3, 1, f, w, i))
                    {
                      remove_degree7_star(v,  4  , f, w, i);
                    }
                  else
                    {
                      remove_degree7_leftfan(v,  4  , f, w, i);
                    }
                }
              else
                {
                  if (incircle(6, 2, 3, 1, f, w, i))
                    {
                      if (incircle(1, 5, 6, 4, f, w, i))
                        {
                          if (incircle(5, 1, 2, 6, f, w, i))
                            {
                              if (incircle(4, 1, 2, 5, f, w, i))
                                {
                                  remove_degree7_rightfan(v,  1  , f, w, i);
                                }
                              else
                                {
                                  if (incircle(4, 2, 3, 5, f, w, i))
                                    {
                                      remove_degree7_zigzag(v,  5  , f, w, i);
                                    }
                                  else
                                    {
                                      remove_degree7_rightfan(v,  5  , f, w, i);
                                    }
                                }
                            }
                          else
                            {
                              if (incircle(5, 2, 3, 6, f, w, i))
                                {
                                  if (incircle(4, 2, 3, 5, f, w, i))
                                    {
                                      remove_degree7_leftfan(v,  2  , f, w, i);
                                    }
                                  else
                                    {
                                      remove_degree7_zigzag(v,  2  , f, w, i);
                                    }
                                }
                              else
                                {
                                  remove_degree7_leftfan(v,  6  , f, w, i);
                                }
                            }
                        }
                      else
                        {
                          if (incircle(4, 1, 2, 6, f, w, i))
                            {
                              remove_degree7_rightdelta(v,  4  , f, w, i);
                            }
                          else
                            {
                              if (incircle(2, 5, 6, 4, f, w, i))
                                {
                                  if (incircle(2, 5, 6, 3, f, w, i))
                                    {
                                      if (incircle(4, 2, 3, 5, f, w, i))
                                        {
                                          remove_degree7_leftfan(v,  2  , f, w, i);
                                        }
                                      else
                                        {
                                          remove_degree7_zigzag(v,  2  , f, w, i);
                                        }
                                    }
                                  else
                                    {
                                      remove_degree7_leftfan(v,  6  , f, w, i);
                                    }
                                }
                              else
                                {
                                  if (incircle(4, 2, 3, 6, f, w, i))
                                    {
                                      remove_degree7_leftdelta(v,  6  , f, w, i);
                                    }
                                  else
                                    {
                                      if (incircle(4, 5, 6, 3, f, w, i))
                                        {
                                          remove_degree7_star(v,  6  , f, w, i);
                                        }
                                      else
                                        {
                                          remove_degree7_leftfan(v,  6  , f, w, i);
                                        }
                                    }
                                }
                            }
                        }
                    }
                  else
                    {
                      if (incircle(1, 5, 6, 4, f, w, i))
                        {
                          if (incircle(1, 5, 6, 3, f, w, i))
                            {
                              if (incircle(5, 2, 3, 1, f, w, i))
                                {
                                  if (incircle(4, 1, 2, 5, f, w, i))
                                    {
                                      remove_degree7_rightfan(v,  1  , f, w, i);
                                    }
                                  else
                                    {
                                      if (incircle(4, 2, 3, 5, f, w, i))
                                        {
                                          remove_degree7_zigzag(v,  5  , f, w, i);
                                        }
                                      else
                                        {
                                          remove_degree7_rightfan(v,  5  , f, w, i);
                                        }
                                    }
                                }
                              else
                                {
                                  if (incircle(1, 4, 5, 3, f, w, i))
                                    {
                                      if (incircle(4, 2, 3, 1, f, w, i))
                                        {
                                          remove_degree7_rightfan(v,  1  , f, w, i);
                                        }
                                      else
                                        {
                                          remove_degree7_star(v,  1  , f, w, i);
                                        }
                                    }
                                  else
                                    {
                                      remove_degree7_rightdelta(v,  1  , f, w, i);
                                    }
                                }
                            }
                          else
                            {
                              remove_degree7_leftdelta(v,  3  , f, w, i);
                            }
                        }
                      else
                        {
                          if (incircle(1, 3, 4, 6, f, w, i))
                            {
                              if (incircle(4, 2, 3, 1, f, w, i))
                                {
                                  remove_degree7_rightdelta(v,  4  , f, w, i);
                                }
                              else
                                {
                                  remove_degree7_leftdelta(v,  1  , f, w, i);
                                }
                            }
                          else
                            {
                              if (incircle(4, 5, 6, 3, f, w, i))
                                {
                                  remove_degree7_rightdelta(v,  6  , f, w, i);
                                }
                              else
                                {
                                  remove_degree7_leftdelta(v,  3  , f, w, i);
                                }
                            }
                        }
                    }
                }
            }
        }
      else
        {
          if (incircle(5, 6, 0, 4, f, w, i))
            {
              if (incircle(5, 6, 0, 3, f, w, i))
                {
                  if (incircle(5, 0, 1, 3, f, w, i))
                    {
                      if (incircle(4, 0, 1, 5, f, w, i))
                        {
                          if (incircle(4, 2, 3, 1, f, w, i))
                            {
                              remove_degree7_rightfan(v,  4  , f, w, i);
                            }
                          else
                            {
                              remove_degree7_zigzag(v,  4  , f, w, i);
                            }
                        }
                      else if (incircle(5, 2, 3, 1, f, w, i))
                        {
                          if (incircle(4, 1, 2, 5, f, w, i))
                            {
                              remove_degree7_zigzag(v,  1  , f, w, i);
                            }
                          else
                            {
                              if (incircle(4, 2, 3, 5, f, w, i))
                                {
                                  remove_degree7_leftfan(v,  5  , f, w, i);
                                }
                              else
                                {
                                  remove_degree7_star(v,  5  , f, w, i);
                                }
                            }
                        }
                      else
                        {
                          if (incircle(1, 4, 5, 3, f, w, i))
                            {
                              if (incircle(4, 2, 3, 1, f, w, i))
                                {
                                  remove_degree7_zigzag(v,  1  , f, w, i);
                                }
                              else
                                {
                                  remove_degree7_leftfan(v,  1  , f, w, i);
                                }
                            }
                          else
                            {
                              remove_degree7_leftdelta(v,  5  , f, w, i);
                            }
                        }
                    }
                  else
                    {
                      if (! incircle(3, 4, 5, 0, f, w, i))
                        {
                          if (incircle(4, 0, 1, 3, f, w, i))
                            {
                              if (incircle(4, 2, 3, 1, f, w, i))
                                {
                                  remove_degree7_rightfan(v,  4  , f, w, i);
                                }
                              else
                                {
                                  remove_degree7_zigzag(v,  4  , f, w, i);
                                }
                            }
                          else
                            {
                              remove_degree7_rightfan(v,  0  , f, w, i);
                            }
                        }
                      else
                        {
                          remove_degree7_rightdelta(v,  3  , f, w, i);
                        }
                    }
                }
              else
                {
                  remove_degree7_star(v,  3  , f, w, i);
                }
            }
          else
            {
              if (incircle(4, 6, 0, 3, f, w, i))
                {
                  if (incircle(4, 0, 1, 3, f, w, i))
                    {
                      if (incircle(4, 2, 3, 1, f, w, i))
                        {
                          remove_degree7_star(v,  4  , f, w, i);
                        }
                      else
                        {
                          remove_degree7_leftfan(v,  4  , f, w, i);
                        }
                    }
                  else
                    {
                      remove_degree7_zigzag(v,  0  , f, w, i);
                    }
                }
              else
                {
                  if (incircle(4, 5, 6, 3, f, w, i))
                    {
                      remove_degree7_rightfan(v,  3  , f, w, i);
                    }
                  else
                    {
                      remove_degree7_star(v,  3  , f, w, i);
                    }
                }
            }
        }
    }
}
template < class Gt, class Tds >
void
Periodic_2_Delaunay_triangulation_2<Gt, Tds>::
remove_degree7(Vertex_handle v, std::vector<Face_handle> &f,
               std::vector<Vertex_handle> &w, std::vector<Offset> &o,
               std::vector<int> &i)
{
  // removing a degree 7 vertex

  if (incircle(2, 0, 1, 3, f, w, o, i)) // sweeping from above
    {
      if (incircle(2, 3, 4, 0, f, w, o, i))
        {
          if (incircle(5, 3, 4, 6, f, w, o, i))
            {
              if (incircle(5, 3, 4, 2, f, w, o, i))
                {
                  if (incircle(6, 2, 3, 5, f, w, o, i))
                    {
                      if (incircle(6, 0, 1, 2, f, w, o, i))
                        {
                          remove_degree7_leftfan(v,  6  , f, w, o, i);
                        }
                      else
                        {
                          remove_degree7_zigzag(v,  6  , f, w, o, i);
                        }
                    }
                  else
                    {
                      if (incircle(5, 0, 1, 2, f, w, o, i))
                        {
                          if (incircle(6, 1, 2, 5, f, w, o, i))
                            {
                              remove_degree7_zigzag(v, 2   , f, w, o, i);
                            }
                          else
                            {
                              if (incircle(6, 0, 1, 5, f, w, o, i))
                                {
                                  remove_degree7_rightfan(v, 5   , f, w, o, i);
                                }
                              else
                                {
                                  remove_degree7_star(v,  5  , f, w, o, i);
                                }
                            }
                        }
                      else
                        {
                          if (incircle(2, 5, 6, 0, f, w, o, i))
                            {
                              if (incircle(6, 0, 1, 2, f, w, o, i))
                                {
                                  remove_degree7_zigzag(v,  2  , f, w, o, i);
                                }
                              else
                                {
                                  remove_degree7_rightfan(v,  2  , f, w, o, i);
                                }
                            }
                          else
                            {
                              remove_degree7_rightdelta(v,  5  , f, w, o, i);
                            }
                        }
                    }
                }
              else
                {
                  if (incircle(4, 0, 1, 2, f, w, o, i))
                    {
                      if (incircle(5, 1, 2, 4, f, w, o, i))
                        {
                          if (incircle(6, 1, 2, 5, f, w, o, i))
                            {
                              remove_degree7_leftfan(v,  2  , f, w, o, i);
                            }
                          else
                            {
                              if (incircle(6, 0, 1, 5, f, w, o, i))
                                {
                                  remove_degree7_zigzag(v,  5  , f, w, o, i);
                                }
                              else
                                {
                                  remove_degree7_leftfan(v,  5  , f, w, o, i);
                                }
                            }
                        }
                      else
                        {
                          if (incircle(5, 0, 1, 4, f, w, o, i))
                            {
                              if (incircle(6, 0, 1, 5, f, w, o, i))
                                {
                                  remove_degree7_rightfan(v,  1  , f, w, o, i);
                                }
                              else
                                {
                                  remove_degree7_zigzag(v,  1  , f, w, o, i);
                                }
                            }
                          else
                            {
                              remove_degree7_rightfan(v,  4  , f, w, o, i);
                            }
                        }
                    }
                  else
                    {
                      if (incircle(2, 4, 5, 0, f, w, o, i))
                        {
                          if (incircle(5, 0, 1, 2, f, w, o, i))
                            {
                              if (incircle(6, 1, 2, 5, f, w, o, i))
                                {
                                  remove_degree7_leftfan(v,  2  , f, w, o, i);
                                }
                              else
                                {
                                  if (incircle(6, 0, 1, 5, f, w, o, i))
                                    {
                                      remove_degree7_zigzag(v,  5  , f, w, o, i);
                                    }
                                  else
                                    {
                                      remove_degree7_leftfan(v,  5  , f, w, o, i);
                                    }
                                }
                            }
                          else
                            {
                              if (incircle(2, 5, 6, 0, f, w, o, i))
                                {
                                  if (incircle(6, 0, 1, 2, f, w, o, i))
                                    {
                                      remove_degree7_leftfan(v,  2  , f, w, o, i);
                                    }
                                  else
                                    {
                                      remove_degree7_star(v,  2  , f, w, o, i);
                                    }
                                }
                              else
                                {
                                  remove_degree7_leftdelta(v,  2  , f, w, o, i);
                                }
                            }
                        }
                      else
                        {
                          remove_degree7_rightdelta(v,  0  , f, w, o, i);
                        }
                    }
                }
            }
          else
            {
              if (incircle(6, 3, 4, 2, f, w, o, i))
                {
                  if (incircle(6, 0, 1, 2, f, w, o, i))
                    {
                      remove_degree7_star(v,  6  , f, w, o, i);
                    }
                  else
                    {
                      remove_degree7_rightfan(v,  6  , f, w, o, i);
                    }
                }
              else
                {
                  if (incircle(4, 0, 1, 2, f, w, o, i))
                    {
                      if (incircle(2, 4, 5, 6, f, w, o, i))
                        {
                          if (incircle(5, 1, 2, 4, f, w, o, i))
                            {
                              if (incircle(6, 1, 2, 5, f, w, o, i))
                                {
                                  remove_degree7_leftfan(v,  2  , f, w, o, i);
                                }
                              else
                                {
                                  if (incircle(6, 0, 1, 5, f, w, o, i))
                                    {
                                      remove_degree7_zigzag(v,  5  , f, w, o, i);
                                    }
                                  else
                                    {
                                      remove_degree7_leftfan(v,  5  , f, w, o, i);
                                    }
                                }
                            }
                          else
                            {
                              if (incircle(5, 0, 1, 4, f, w, o, i))
                                {
                                  if (incircle(6, 0, 1, 5, f, w, o, i))
                                    {
                                      remove_degree7_rightfan(v,  1  , f, w, o, i);
                                    }
                                  else
                                    {
                                      remove_degree7_zigzag(v,  1  , f, w, o, i);
                                    }
                                }
                              else
                                {
                                  remove_degree7_rightfan(v,  4  , f, w, o, i);
                                }
                            }
                        }
                      else
                        {
                          if (incircle(6, 1, 2, 4, f, w, o, i))
                            {
                              remove_degree7_leftdelta(v,  6  , f, w, o, i);
                            }
                          else
                            {
                              if (incircle(1, 4, 5, 6, f, w, o, i))
                                {
                                  if (incircle(1, 4, 5, 0, f, w, o, i))
                                    {
                                      if (incircle(6, 0, 1, 5, f, w, o, i))
                                        {
                                          remove_degree7_rightfan(v,  1  , f, w, o, i);
                                        }
                                      else
                                        {
                                          remove_degree7_zigzag(v,  1  , f, w, o, i);
                                        }
                                    }
                                  else
                                    {
                                      remove_degree7_rightfan(v,  4  , f, w, o, i);
                                    }
                                }
                              else
                                {
                                  if (incircle(6, 0, 1, 4, f, w, o, i))
                                    {
                                      remove_degree7_rightdelta(v,  4  , f, w, o, i);
                                    }
                                  else
                                    {
                                      if (incircle(6, 4, 5, 0, f, w, o, i))
                                        {
                                          remove_degree7_star(v,  4  , f, w, o, i);
                                        }
                                      else
                                        {
                                          remove_degree7_rightfan(v,  4  , f, w, o, i);
                                        }
                                    }
                                }
                            }
                        }
                    }
                  else
                    {
                      if (incircle(2, 4, 5, 6, f, w, o, i))
                        {
                          if (incircle(2, 4, 5, 0, f, w, o, i))
                            {
                              if (incircle(5, 0, 1, 2, f, w, o, i))
                                {
                                  if (incircle(6, 1, 2, 5, f, w, o, i))
                                    {
                                      remove_degree7_leftfan(v,  2  , f, w, o, i);
                                    }
                                  else
                                    {
                                      if (incircle(6, 0, 1, 5, f, w, o, i))
                                        {
                                          remove_degree7_zigzag(v,  5  , f, w, o, i);
                                        }
                                      else
                                        {
                                          remove_degree7_leftfan(v,  5  , f, w, o, i);
                                        }
                                    }
                                }
                              else
                                {
                                  if (incircle(2, 5, 6, 0, f, w, o, i))
                                    {
                                      if (incircle(6, 0, 1, 2, f, w, o, i))
                                        {
                                          remove_degree7_leftfan(v,  2  , f, w, o, i);
                                        }
                                      else
                                        {
                                          remove_degree7_star(v,  2  , f, w, o, i);
                                        }
                                    }
                                  else
                                    {
                                      remove_degree7_leftdelta(v,  2  , f, w, o, i);
                                    }
                                }
                            }
                          else
                            {
                              remove_degree7_rightdelta(v,  0  , f, w, o, i);
                            }
                        }
                      else
                        {
                          if (incircle(2, 6, 0, 4, f, w, o, i))
                            {
                              if (incircle(6, 0, 1, 2, f, w, o, i))
                                {
                                  remove_degree7_leftdelta(v,  6  , f, w, o, i);
                                }
                              else
                                {
                                  remove_degree7_rightdelta(v,  2  , f, w, o, i);
                                }
                            }
                          else
                            {
                              if (incircle(6, 4, 5, 0, f, w, o, i))
                                {
                                  remove_degree7_leftdelta(v,  4  , f, w, o, i);
                                }
                              else
                                {
                                  remove_degree7_rightdelta(v,  0  , f, w, o, i);
                                }
                            }
                        }
                    }
                }
            }
        }
      else
        {
          if (incircle(5, 3, 4, 6, f, w, o, i))
            {
              if (incircle(5, 3, 4, 0, f, w, o, i))
                {
                  if (incircle(5, 2, 3, 0, f, w, o, i))
                    {
                      if (incircle(6, 2, 3, 5, f, w, o, i))
                        {
                          if (incircle(6, 0, 1, 2, f, w, o, i))
                            {
                              remove_degree7_leftfan(v,  6  , f, w, o, i);
                            }
                          else
                            {
                              remove_degree7_zigzag(v,  6  , f, w, o, i);
                            }
                        }
                      else if (incircle(5, 0, 1, 2, f, w, o, i))
                        {
                          if (incircle(6, 1, 2, 5, f, w, o, i))
                            {
                              remove_degree7_zigzag(v,  2  , f, w, o, i);
                            }
                          else
                            {
                              if (incircle(6, 0, 1, 5, f, w, o, i))
                                {
                                  remove_degree7_rightfan(v,  5  , f, w, o, i);
                                }
                              else
                                {
                                  remove_degree7_star(v,  5  , f, w, o, i);
                                }
                            }
                        }
                      else
                        {
                          if (incircle(2, 5, 6, 0, f, w, o, i))
                            {
                              if (incircle(6, 0, 1, 2, f, w, o, i))
                                {
                                  remove_degree7_zigzag(v,  2  , f, w, o, i);
                                }
                              else
                                {
                                  remove_degree7_rightfan(v,  2  , f, w, o, i);
                                }
                            }
                          else
                            {
                              remove_degree7_rightdelta(v,  5  , f, w, o, i);
                            }
                        }
                    }
                  else
                    {
                      if (incircle(3, 5, 6, 0, f, w, o, i))
                        {
                          if (incircle(6, 2, 3, 0, f, w, o, i))
                            {
                              if (incircle(6, 0, 1, 2, f, w, o, i))
                                {
                                  remove_degree7_leftfan(v,  6  , f, w, o, i);
                                }
                              else
                                {
                                  remove_degree7_zigzag(v,  6  , f, w, o, i);
                                }
                            }
                          else
                            {
                              remove_degree7_leftfan(v,  3  , f, w, o, i);
                            }
                        }
                      else
                        {
                          remove_degree7_leftdelta(v,  0  , f, w, o, i);
                        }
                    }
                }
              else
                {
                  remove_degree7_star(v,  0  , f, w, o, i);
                }
            }
          else
            {
              if (incircle(6, 3, 4, 0, f, w, o, i))
                {
                  if (incircle(6, 2, 3, 0, f, w, o, i))
                    {
                      if (incircle(6, 0, 1, 2, f, w, o, i))
                        {
                          remove_degree7_star(v,  6  , f, w, o, i);
                        }
                      else
                        {
                          remove_degree7_rightfan(v,  6  , f, w, o, i);
                        }
                    }
                  else
                    {
                      remove_degree7_zigzag(v,  3  , f, w, o, i);
                    }
                }
              else
                {
                  if (incircle(6, 4, 5, 0, f, w, o, i))
                    {
                      remove_degree7_leftfan(v,  0  , f, w, o, i);
                    }
                  else
                    {
                      remove_degree7_star(v,  0  , f, w, o, i);
                    }
                }
            }
        }
    }
  else    //sweeping from below
    {
      if (incircle(1, 6, 0, 3, f, w, o, i))
        {
          if (incircle(5, 6, 0, 4, f, w, o, i))
            {
              if (incircle(5, 6, 0, 1, f, w, o, i))
                {
                  if (incircle(4, 0, 1, 5, f, w, o, i))
                    {
                      if (incircle(4, 2, 3, 1, f, w, o, i))
                        {
                          remove_degree7_rightfan(v,  4  , f, w, o, i);
                        }
                      else
                        {
                          remove_degree7_zigzag(v,  4  , f, w, o, i);
                        }
                    }
                  else
                    {
                      if (incircle(5, 2, 3, 1, f, w, o, i))
                        {
                          if (incircle(4, 1, 2, 5, f, w, o, i))
                            {
                              remove_degree7_zigzag(v, 1   , f, w, o, i);
                            }
                          else
                            {
                              if (incircle(4, 2, 3, 5, f, w, o, i))
                                {
                                  remove_degree7_leftfan(v, 5   , f, w, o, i);
                                }
                              else
                                {
                                  remove_degree7_star(v,  5  , f, w, o, i);
                                }
                            }
                        }
                      else
                        {
                          if (incircle(1, 4, 5, 3, f, w, o, i))
                            {
                              if (incircle(4, 2, 3, 1, f, w, o, i))
                                {
                                  remove_degree7_zigzag(v,  1  , f, w, o, i);
                                }
                              else
                                {
                                  remove_degree7_leftfan(v,  1  , f, w, o, i);
                                }
                            }
                          else
                            {
                              remove_degree7_leftdelta(v,  5  , f, w, o, i);
                            }
                        }
                    }
                }
              else
                {
                  if (incircle(6, 2, 3, 1, f, w, o, i))
                    {
                      if (incircle(5, 1, 2, 6, f, w, o, i))
                        {
                          if (incircle(4, 1, 2, 5, f, w, o, i))
                            {
                              remove_degree7_rightfan(v,  1  , f, w, o, i);
                            }
                          else
                            {
                              if (incircle(4, 2, 3, 5, f, w, o, i))
                                {
                                  remove_degree7_zigzag(v,  5  , f, w, o, i);
                                }
                              else
                                {
                                  remove_degree7_rightfan(v,  5  , f, w, o, i);
                                }
                            }
                        }
                      else
                        {
                          if (incircle(5, 2, 3, 6, f, w, o, i))
                            {
                              if (incircle(4, 2, 3, 5, f, w, o, i))
                                {
                                  remove_degree7_leftfan(v,  2  , f, w, o, i);
                                }
                              else
                                {
                                  remove_degree7_zigzag(v,  2  , f, w, o, i);
                                }
                            }
                          else
                            {
                              remove_degree7_leftfan(v,  6  , f, w, o, i);
                            }
                        }
                    }
                  else
                    {
                      if (incircle(1, 5, 6, 3, f, w, o, i))
                        {
                          if (incircle(5, 2, 3, 1, f, w, o, i))
                            {
                              if (incircle(4, 1, 2, 5, f, w, o, i))
                                {
                                  remove_degree7_rightfan(v,  1  , f, w, o, i);
                                }
                              else
                                {
                                  if (incircle(4, 2, 3, 5, f, w, o, i))
                                    {
                                      remove_degree7_zigzag(v,  5  , f, w, o, i);
                                    }
                                  else
                                    {
                                      remove_degree7_rightfan(v,  5  , f, w, o, i);
                                    }
                                }
                            }
                          else
                            {
                              if (incircle(1, 4, 5, 3, f, w, o, i))
                                {
                                  if (incircle(4, 2, 3, 1, f, w, o, i))
                                    {
                                      remove_degree7_rightfan(v,  1  , f, w, o, i);
                                    }
                                  else
                                    {
                                      remove_degree7_star(v,  1  , f, w, o, i);
                                    }
                                }
                              else
                                {
                                  remove_degree7_rightdelta(v,  1  , f, w, o, i);
                                }
                            }
                        }
                      else
                        {
                          remove_degree7_leftdelta(v,  3  , f, w, o, i);
                        }
                    }
                }
            }
          else
            {
              if (incircle(4, 6, 0, 1, f, w, o, i))
                {
                  if (incircle(4, 2, 3, 1, f, w, o, i))
                    {
                      remove_degree7_star(v,  4  , f, w, o, i);
                    }
                  else
                    {
                      remove_degree7_leftfan(v,  4  , f, w, o, i);
                    }
                }
              else
                {
                  if (incircle(6, 2, 3, 1, f, w, o, i))
                    {
                      if (incircle(1, 5, 6, 4, f, w, o, i))
                        {
                          if (incircle(5, 1, 2, 6, f, w, o, i))
                            {
                              if (incircle(4, 1, 2, 5, f, w, o, i))
                                {
                                  remove_degree7_rightfan(v,  1  , f, w, o, i);
                                }
                              else
                                {
                                  if (incircle(4, 2, 3, 5, f, w, o, i))
                                    {
                                      remove_degree7_zigzag(v,  5  , f, w, o, i);
                                    }
                                  else
                                    {
                                      remove_degree7_rightfan(v,  5  , f, w, o, i);
                                    }
                                }
                            }
                          else
                            {
                              if (incircle(5, 2, 3, 6, f, w, o, i))
                                {
                                  if (incircle(4, 2, 3, 5, f, w, o, i))
                                    {
                                      remove_degree7_leftfan(v,  2  , f, w, o, i);
                                    }
                                  else
                                    {
                                      remove_degree7_zigzag(v,  2  , f, w, o, i);
                                    }
                                }
                              else
                                {
                                  remove_degree7_leftfan(v,  6  , f, w, o, i);
                                }
                            }
                        }
                      else
                        {
                          if (incircle(4, 1, 2, 6, f, w, o, i))
                            {
                              remove_degree7_rightdelta(v,  4  , f, w, o, i);
                            }
                          else
                            {
                              if (incircle(2, 5, 6, 4, f, w, o, i))
                                {
                                  if (incircle(2, 5, 6, 3, f, w, o, i))
                                    {
                                      if (incircle(4, 2, 3, 5, f, w, o, i))
                                        {
                                          remove_degree7_leftfan(v,  2  , f, w, o, i);
                                        }
                                      else
                                        {
                                          remove_degree7_zigzag(v,  2  , f, w, o, i);
                                        }
                                    }
                                  else
                                    {
                                      remove_degree7_leftfan(v,  6  , f, w, o, i);
                                    }
                                }
                              else
                                {
                                  if (incircle(4, 2, 3, 6, f, w, o, i))
                                    {
                                      remove_degree7_leftdelta(v,  6  , f, w, o, i);
                                    }
                                  else
                                    {
                                      if (incircle(4, 5, 6, 3, f, w, o, i))
                                        {
                                          remove_degree7_star(v,  6  , f, w, o, i);
                                        }
                                      else
                                        {
                                          remove_degree7_leftfan(v,  6  , f, w, o, i);
                                        }
                                    }
                                }
                            }
                        }
                    }
                  else
                    {
                      if (incircle(1, 5, 6, 4, f, w, o, i))
                        {
                          if (incircle(1, 5, 6, 3, f, w, o, i))
                            {
                              if (incircle(5, 2, 3, 1, f, w, o, i))
                                {
                                  if (incircle(4, 1, 2, 5, f, w, o, i))
                                    {
                                      remove_degree7_rightfan(v,  1  , f, w, o, i);
                                    }
                                  else
                                    {
                                      if (incircle(4, 2, 3, 5, f, w, o, i))
                                        {
                                          remove_degree7_zigzag(v,  5  , f, w, o, i);
                                        }
                                      else
                                        {
                                          remove_degree7_rightfan(v,  5  , f, w, o, i);
                                        }
                                    }
                                }
                              else
                                {
                                  if (incircle(1, 4, 5, 3, f, w, o, i))
                                    {
                                      if (incircle(4, 2, 3, 1, f, w, o, i))
                                        {
                                          remove_degree7_rightfan(v,  1  , f, w, o, i);
                                        }
                                      else
                                        {
                                          remove_degree7_star(v,  1  , f, w, o, i);
                                        }
                                    }
                                  else
                                    {
                                      remove_degree7_rightdelta(v,  1  , f, w, o, i);
                                    }
                                }
                            }
                          else
                            {
                              remove_degree7_leftdelta(v,  3  , f, w, o, i);
                            }
                        }
                      else
                        {
                          if (incircle(1, 3, 4, 6, f, w, o, i))
                            {
                              if (incircle(4, 2, 3, 1, f, w, o, i))
                                {
                                  remove_degree7_rightdelta(v,  4  , f, w, o, i);
                                }
                              else
                                {
                                  remove_degree7_leftdelta(v,  1  , f, w, o, i);
                                }
                            }
                          else
                            {
                              if (incircle(4, 5, 6, 3, f, w, o, i))
                                {
                                  remove_degree7_rightdelta(v,  6  , f, w, o, i);
                                }
                              else
                                {
                                  remove_degree7_leftdelta(v,  3  , f, w, o, i);
                                }
                            }
                        }
                    }
                }
            }
        }
      else
        {
          if (incircle(5, 6, 0, 4, f, w, o, i))
            {
              if (incircle(5, 6, 0, 3, f, w, o, i))
                {
                  if (incircle(5, 0, 1, 3, f, w, o, i))
                    {
                      if (incircle(4, 0, 1, 5, f, w, o, i))
                        {
                          if (incircle(4, 2, 3, 1, f, w, o, i))
                            {
                              remove_degree7_rightfan(v,  4  , f, w, o, i);
                            }
                          else
                            {
                              remove_degree7_zigzag(v,  4  , f, w, o, i);
                            }
                        }
                      else if (incircle(5, 2, 3, 1, f, w, o, i))
                        {
                          if (incircle(4, 1, 2, 5, f, w, o, i))
                            {
                              remove_degree7_zigzag(v,  1  , f, w, o, i);
                            }
                          else
                            {
                              if (incircle(4, 2, 3, 5, f, w, o, i))
                                {
                                  remove_degree7_leftfan(v,  5  , f, w, o, i);
                                }
                              else
                                {
                                  remove_degree7_star(v,  5  , f, w, o, i);
                                }
                            }
                        }
                      else
                        {
                          if (incircle(1, 4, 5, 3, f, w, o, i))
                            {
                              if (incircle(4, 2, 3, 1, f, w, o, i))
                                {
                                  remove_degree7_zigzag(v,  1  , f, w, o, i);
                                }
                              else
                                {
                                  remove_degree7_leftfan(v,  1  , f, w, o, i);
                                }
                            }
                          else
                            {
                              remove_degree7_leftdelta(v,  5  , f, w, o, i);
                            }
                        }
                    }
                  else
                    {
                      if (! incircle(3, 4, 5, 0, f, w, o, i))
                        {
                          if (incircle(4, 0, 1, 3, f, w, o, i))
                            {
                              if (incircle(4, 2, 3, 1, f, w, o, i))
                                {
                                  remove_degree7_rightfan(v,  4  , f, w, o, i);
                                }
                              else
                                {
                                  remove_degree7_zigzag(v,  4  , f, w, o, i);
                                }
                            }
                          else
                            {
                              remove_degree7_rightfan(v,  0  , f, w, o, i);
                            }
                        }
                      else
                        {
                          remove_degree7_rightdelta(v,  3  , f, w, o, i);
                        }
                    }
                }
              else
                {
                  remove_degree7_star(v,  3  , f, w, o, i);
                }
            }
          else
            {
              if (incircle(4, 6, 0, 3, f, w, o, i))
                {
                  if (incircle(4, 0, 1, 3, f, w, o, i))
                    {
                      if (incircle(4, 2, 3, 1, f, w, o, i))
                        {
                          remove_degree7_star(v,  4  , f, w, o, i);
                        }
                      else
                        {
                          remove_degree7_leftfan(v,  4  , f, w, o, i);
                        }
                    }
                  else
                    {
                      remove_degree7_zigzag(v,  0  , f, w, o, i);
                    }
                }
              else
                {
                  if (incircle(4, 5, 6, 3, f, w, o, i))
                    {
                      remove_degree7_rightfan(v,  3  , f, w, o, i);
                    }
                  else
                    {
                      remove_degree7_star(v,  3  , f, w, o, i);
                    }
                }
            }
        }
    }
}



template < class Gt, class Tds >
void
Periodic_2_Delaunay_triangulation_2<Gt, Tds>::
rotate7(int j,  std::vector<Vertex_handle> &w,
        std::vector<Face_handle> &f, std::vector<int> &i)
{
  if (j == 0) return;
  Face_handle ff = f[0];
  int ii = i[0], k = 0, kk = (6 * j) % 7;
  Vertex_handle ww = w[0];
  for (int jj = 0; k != kk; jj = k) // 7 is prime
    {
      k = (jj + j) % 7;
      w[jj] = w[k];
      f[jj] = f[k];
      i[jj] = i[k];
    }
  w[kk] = ww;
  f[kk] = ff;
  i[kk] = ii;
}
template < class Gt, class Tds >
void
Periodic_2_Delaunay_triangulation_2<Gt, Tds>::
rotate7(int j,  std::vector<Vertex_handle> &w,
        std::vector<Face_handle> &f, std::vector<Offset> &o, std::vector<int> &i)
{
  if (j == 0) return;
  Face_handle ff = f[0];
  int ii = i[0], k = 0, kk = (6 * j) % 7;
  Vertex_handle ww = w[0];
  Offset oo = o[0];
  for (int jj = 0; k != kk; jj = k) // 7 is prime
    {
      k = (jj + j) % 7;
      w[jj] = w[k];
      f[jj] = f[k];
      o[jj] = o[k];
      i[jj] = i[k];
    }
  w[kk] = ww;
  f[kk] = ff;
  o[kk] = oo;
  i[kk] = ii;
}

template < class Gt, class Tds >
void
Periodic_2_Delaunay_triangulation_2<Gt, Tds>::
get_offset_degree7(std::vector<Offset> &in_o, int out_o[])
{
  bool add[2];

  add[0] = add[1] = false;

  for (int cnt = 0; cnt < 7; ++cnt)
    {
      add[0] |= in_o[cnt].x() < 0;
      add[1] |= in_o[cnt].y() < 0;
    }

  Covering_sheets c = number_of_sheets();
  if (add[0] || add[1])
    {
      const Offset oo = Offset(add[0] ? c[0] : 0, add[1] ? c[1] : 0);
      for (int i = 0; i < 7; ++i) in_o[i] += oo;
    }

  for (int cnt = 0; cnt < 7; ++cnt)
    out_o[cnt] = (in_o[cnt].x() >= c[0] ? 2 : 0) + (in_o[cnt].y() >= c[1] ? 1 : 0);
}

template < class Gt, class Tds >
void
Periodic_2_Delaunay_triangulation_2<Gt, Tds>::
remove_degree7_star   (Vertex_handle &, int j,
                       std::vector<Face_handle> &f, std::vector<Vertex_handle> &w, std::vector<int> &i)
{
  // removing a degree 7 vertex, staring from w[j]
  rotate7(j, w, f, i);

  Face_handle nn;
  f[1]->set_vertex( i[1], w[0]) ;  // f1 = w1w2w0
  f[2]->set_vertex( i[2], w[0]) ;  // f2 = w2w3w0
  f[3]->set_vertex( i[3], w[0]) ;  // f3 = w3w4w0
  f[4]->set_vertex( i[4], w[0]) ;  // f4 = w4w5w0
  f[5]->set_vertex( i[5], w[0]) ;  // f5 = w5w6w0

  nn = f[0]->neighbor( i[0] );
  tds().set_adjacency(f[1], cw(i[1]) , nn , nn->index(f[0])  );
  nn = f[6]->neighbor( i[6] );
  tds().set_adjacency(f[5], ccw(i[5]) , nn , nn->index(f[6]) );
  tds().delete_face(f[0]);
  tds().delete_face(f[6]);


  f[1]->set_offsets(0, 0, 0);
  f[2]->set_offsets(0, 0, 0);
  f[3]->set_offsets(0, 0, 0);
  f[4]->set_offsets(0, 0, 0);
  f[5]->set_offsets(0, 0, 0);
}
template < class Gt, class Tds >
void
Periodic_2_Delaunay_triangulation_2<Gt, Tds>::
remove_degree7_star   (Vertex_handle &v, int j,
                       std::vector<Face_handle> &f, std::vector<Vertex_handle> &w, std::vector<Offset> &o, std::vector<int> &i)
{
  // removing a degree 7 vertex, staring from w[j]

  // Rotate the offset as well
  rotate7(j, w, f, o, i);
  remove_degree7_star(v, /* !! */ 0, f, w, i);

  int oo[7];
  get_offset_degree7(o, oo);
  int o_face[3];
  int ii;
  ii = i[1];
  o_face[ii] = oo[0];
  o_face[ccw(ii)] = oo[1];
  o_face[ cw(ii)] = oo[2];
  this->set_offsets(f[1], o_face[0], o_face[1], o_face[2]);
  ii = i[2];
  o_face[ii] = oo[0];
  o_face[ccw(ii)] = oo[2];
  o_face[ cw(ii)] = oo[3];
  this->set_offsets(f[2], o_face[0], o_face[1], o_face[2]);
  ii = i[3];
  o_face[ii] = oo[0];
  o_face[ccw(ii)] = oo[3];
  o_face[ cw(ii)] = oo[4];
  this->set_offsets(f[3], o_face[0], o_face[1], o_face[2]);
  ii = i[4];
  o_face[ii] = oo[0];
  o_face[ccw(ii)] = oo[4];
  o_face[ cw(ii)] = oo[5];
  this->set_offsets(f[4], o_face[0], o_face[1], o_face[2]);
  ii = i[5];
  o_face[ii] = oo[0];
  o_face[ccw(ii)] = oo[5];
  o_face[ cw(ii)] = oo[6];
  this->set_offsets(f[5], o_face[0], o_face[1], o_face[2]);

  insert_too_long_edge(f[1], ccw(i[1]));
  insert_too_long_edge(f[2], ccw(i[2]));
  insert_too_long_edge(f[3], ccw(i[3]));
  insert_too_long_edge(f[4], ccw(i[4]));
}
template < class Gt, class Tds >
void
Periodic_2_Delaunay_triangulation_2<Gt, Tds>::
remove_degree7_zigzag (Vertex_handle &, int j,
                       std::vector<Face_handle> &f, std::vector<Vertex_handle> &w, std::vector<int> &i)
{
  // removing a degree 7 vertex, zigzag, w[j] = middle point
  rotate7(j, w, f, i);

  Face_handle nn;
  f[1]->set_vertex(    i[1] , w[3]) ;  // f1 = w1w2w3
  f[2]->set_vertex(ccw(i[2]), w[1]) ;
  f[2]->set_vertex(    i[2] , w[0]) ;  // f2 = w1w3w0
  f[3]->set_vertex(    i[3] , w[0]) ;  // f3 = w3w4w0
  f[4]->set_vertex( cw(i[4]), w[6]) ;
  f[4]->set_vertex(    i[4] , w[0]) ;  // f4 = w4w6w0
  f[5]->set_vertex(    i[5] , w[4]) ;  // f5 = w5w6w4

  nn = f[2]->neighbor( i[2] );
  tds().set_adjacency(f[1], ccw(i[1]) , nn, nn->index(f[2]) );
  nn = f[0]->neighbor( i[0] );
  tds().set_adjacency(f[2], cw(i[2]) , nn , nn->index(f[0]) );
  nn = f[6]->neighbor( i[6] );
  tds().set_adjacency(f[4], ccw(i[4]) , nn , nn->index(f[6])  );
  nn = f[4]->neighbor( i[4] );
  tds().set_adjacency(f[5], cw(i[5]) , nn , nn->index(f[4])  );
  tds().set_adjacency(f[1], cw(i[1]) , f[2] , i[2]   );
  tds().set_adjacency(f[4], i[4]  , f[5] , ccw(i[5])  );

  tds().delete_face(f[0]);
  tds().delete_face(f[6]);


  f[1]->set_offsets(0, 0, 0);
  f[2]->set_offsets(0, 0, 0);
  f[3]->set_offsets(0, 0, 0);
  f[4]->set_offsets(0, 0, 0);
  f[5]->set_offsets(0, 0, 0);
}
template < class Gt, class Tds >
void
Periodic_2_Delaunay_triangulation_2<Gt, Tds>::
remove_degree7_zigzag (Vertex_handle &v, int j,
                       std::vector<Face_handle> &f, std::vector<Vertex_handle> &w, std::vector<Offset> &o, std::vector<int> &i)
{
  // removing a degree 7 vertex, zigzag, w[j] = middle point

  // Rotate the offset as well
  rotate7(j, w, f, o, i);
  remove_degree7_zigzag(v, /* !! */ 0, f, w, i);

  int oo[7];
  get_offset_degree7(o, oo);
  int o_face[3];
  int ii;
  ii = i[1];
  o_face[ii] = oo[3];
  o_face[ccw(ii)] = oo[1];
  o_face[ cw(ii)] = oo[2];
  this->set_offsets(f[1], o_face[0], o_face[1], o_face[2]);
  ii = i[2];
  o_face[ii] = oo[0];
  o_face[ccw(ii)] = oo[1];
  o_face[ cw(ii)] = oo[3];
  this->set_offsets(f[2], o_face[0], o_face[1], o_face[2]);
  ii = i[3];
  o_face[ii] = oo[0];
  o_face[ccw(ii)] = oo[3];
  o_face[ cw(ii)] = oo[4];
  this->set_offsets(f[3], o_face[0], o_face[1], o_face[2]);
  ii = i[4];
  o_face[ii] = oo[0];
  o_face[ccw(ii)] = oo[4];
  o_face[ cw(ii)] = oo[6];
  this->set_offsets(f[4], o_face[0], o_face[1], o_face[2]);
  ii = i[5];
  o_face[ii] = oo[4];
  o_face[ccw(ii)] = oo[5];
  o_face[ cw(ii)] = oo[6];
  this->set_offsets(f[5], o_face[0], o_face[1], o_face[2]);

  insert_too_long_edge(f[1],  cw(i[1]));
  insert_too_long_edge(f[2], ccw(i[2]));
  insert_too_long_edge(f[3], ccw(i[3]));
  insert_too_long_edge(f[4],     i[4]);
}
template < class Gt, class Tds >
void
Periodic_2_Delaunay_triangulation_2<Gt, Tds>::
remove_degree7_leftdelta(Vertex_handle &, int j,
                         std::vector<Face_handle> &f, std::vector<Vertex_handle> &w, std::vector<int> &i)
{
  // removing a degree 7 vertex, left delta from w[j]
  rotate7(j, w, f, i);

  Face_handle nn;
  f[1]->set_vertex(    i[1] , w[0]) ;  // f1 = w1w2w0
  f[2]->set_vertex(    i[2] , w[0]) ;  // f2 = w2w3w0
  f[3]->set_vertex( cw(i[3]), w[5]) ;
  f[3]->set_vertex(    i[3] , w[0]) ;  // f3 = w3w5w0
  f[4]->set_vertex(    i[4] , w[3]) ;  // f4 = w4w5w3
  f[5]->set_vertex(    i[5] , w[0]) ;  // f5 = w5w6w0

  nn = f[0]->neighbor( i[0] );
  tds().set_adjacency(f[1], cw(i[1]) , nn , nn->index(f[0])  );
  nn = f[3]->neighbor( i[3] );
  tds().set_adjacency(f[4], cw(i[4]) , nn , nn->index(f[3]) );
  nn = f[6]->neighbor( i[6] );
  tds().set_adjacency(f[5], ccw(i[5]) , nn , nn->index(f[6])  );
  tds().set_adjacency(f[3], i[3]  , f[4] , ccw(i[4])  );
  tds().set_adjacency(f[3], ccw(i[3]) , f[5] ,  cw(i[5]) );

  tds().delete_face(f[0]);
  tds().delete_face(f[6]);

  f[1]->set_offsets(0, 0, 0);
  f[2]->set_offsets(0, 0, 0);
  f[3]->set_offsets(0, 0, 0);
  f[4]->set_offsets(0, 0, 0);
  f[5]->set_offsets(0, 0, 0);
}
template < class Gt, class Tds >
void
Periodic_2_Delaunay_triangulation_2<Gt, Tds>::
remove_degree7_leftdelta(Vertex_handle &v, int j,
                         std::vector<Face_handle> &f, std::vector<Vertex_handle> &w, std::vector<Offset> &o, std::vector<int> &i)
{
  // removing a degree 7 vertex, left delta from w[j]

  // Rotate the offset as well
  rotate7(j, w, f, o, i);
  remove_degree7_leftdelta(v, /* !! */ 0, f, w, i);

  int oo[7];
  get_offset_degree7(o, oo);

  int o_face[3];
  int ii;
  ii = i[1];
  o_face[ii] = oo[0];
  o_face[ccw(ii)] = oo[1];
  o_face[ cw(ii)] = oo[2];
  this->set_offsets(f[1], o_face[0], o_face[1], o_face[2]);
  ii = i[2];
  o_face[ii] = oo[0];
  o_face[ccw(ii)] = oo[2];
  o_face[ cw(ii)] = oo[3];
  this->set_offsets(f[2], o_face[0], o_face[1], o_face[2]);
  ii = i[3];
  o_face[ii] = oo[0];
  o_face[ccw(ii)] = oo[3];
  o_face[ cw(ii)] = oo[5];
  this->set_offsets(f[3], o_face[0], o_face[1], o_face[2]);
  ii = i[4];
  o_face[ii] = oo[3];
  o_face[ccw(ii)] = oo[4];
  o_face[ cw(ii)] = oo[5];
  this->set_offsets(f[4], o_face[0], o_face[1], o_face[2]);
  ii = i[5];
  o_face[ii] = oo[0];
  o_face[ccw(ii)] = oo[5];
  o_face[ cw(ii)] = oo[6];
  this->set_offsets(f[5], o_face[0], o_face[1], o_face[2]);

  insert_too_long_edge(f[1], ccw(i[1]));
  insert_too_long_edge(f[2], ccw(i[2]));
  insert_too_long_edge(f[3],     i[3]);
  insert_too_long_edge(f[3], ccw(i[3]));
}
template < class Gt, class Tds >
void
Periodic_2_Delaunay_triangulation_2<Gt, Tds>::
remove_degree7_rightdelta(Vertex_handle &, int j,
                          std::vector<Face_handle> &f, std::vector<Vertex_handle> &w, std::vector<int> &i)
{
  // removing a degree 7 vertex, right delta from w[j]
  rotate7(j, w, f, i);

  Face_handle nn;
  f[1]->set_vertex(    i[1] , w[0]) ;  // f1 = w1w2w0
  f[2]->set_vertex(    i[2] , w[4]) ;  // f2 = w2w3w4
  f[3]->set_vertex(ccw(i[3]), w[2]) ;
  f[3]->set_vertex(    i[3] , w[0]) ;  // f3 = w2w4w0
  f[4]->set_vertex(    i[4] , w[0]) ;  // f4 = w4w5w0
  f[5]->set_vertex(    i[5] , w[0]) ;  // f5 = w5w6w0

  nn = f[0]->neighbor( i[0] );
  tds().set_adjacency(f[1], cw(i[1]) , nn , nn->index(f[0]) );
  nn = f[3]->neighbor( i[3] );
  tds().set_adjacency(f[2], ccw(i[2]) , nn, nn->index(f[3]) );
  nn = f[6]->neighbor( i[6] );
  tds().set_adjacency(f[5], ccw(i[5]) , nn , nn->index(f[6]) );
  tds().set_adjacency(f[1], ccw(i[1]) , f[3], cw(i[3])  );
  tds().set_adjacency(f[3], i[3]  , f[2], cw(i[2]) );

  tds().delete_face(f[0]);
  tds().delete_face(f[6]);

  f[1]->set_offsets(0, 0, 0);
  f[2]->set_offsets(0, 0, 0);
  f[3]->set_offsets(0, 0, 0);
  f[4]->set_offsets(0, 0, 0);
  f[5]->set_offsets(0, 0, 0);
}
template < class Gt, class Tds >
void
Periodic_2_Delaunay_triangulation_2<Gt, Tds>::
remove_degree7_rightdelta(Vertex_handle &v, int j,
                          std::vector<Face_handle> &f, std::vector<Vertex_handle> &w, std::vector<Offset> &o, std::vector<int> &i)
{
  // removing a degree 7 vertex, right delta from w[j]

  // Rotate the offset as well
  rotate7(j, w, f, o, i);
  remove_degree7_rightdelta(v, /* !! */ 0, f, w, i);

  int oo[7];
  get_offset_degree7(o, oo);
  int o_face[3];
  int ii;
  ii = i[1];
  o_face[ii] = oo[0];
  o_face[ccw(ii)] = oo[1];
  o_face[ cw(ii)] = oo[2];
  this->set_offsets(f[1], o_face[0], o_face[1], o_face[2]);
  ii = i[2];
  o_face[ii] = oo[4];
  o_face[ccw(ii)] = oo[2];
  o_face[ cw(ii)] = oo[3];
  this->set_offsets(f[2], o_face[0], o_face[1], o_face[2]);
  ii = i[3];
  o_face[ii] = oo[0];
  o_face[ccw(ii)] = oo[2];
  o_face[ cw(ii)] = oo[4];
  this->set_offsets(f[3], o_face[0], o_face[1], o_face[2]);
  ii = i[4];
  o_face[ii] = oo[0];
  o_face[ccw(ii)] = oo[4];
  o_face[ cw(ii)] = oo[5];
  this->set_offsets(f[4], o_face[0], o_face[1], o_face[2]);
  ii = i[5];
  o_face[ii] = oo[0];
  o_face[ccw(ii)] = oo[5];
  o_face[ cw(ii)] = oo[6];
  this->set_offsets(f[5], o_face[0], o_face[1], o_face[2]);

  insert_too_long_edge(f[1], ccw(i[1]));
  insert_too_long_edge(f[2],  cw(i[2]));
  insert_too_long_edge(f[3], ccw(i[3]));
  insert_too_long_edge(f[4], ccw(i[4]));
}
template < class Gt, class Tds >
void
Periodic_2_Delaunay_triangulation_2<Gt, Tds>::
remove_degree7_leftfan(Vertex_handle &, int j,
                       std::vector<Face_handle> &f, std::vector<Vertex_handle> &w, std::vector<int> &i)
{
  // removing a degree 7 vertex, left fan from w[j]
  rotate7(j, w, f, i);

  Face_handle nn;
  f[1]->set_vertex(    i[1] , w[0]) ;  // f1 = w1w2w0
  f[2]->set_vertex(    i[2] , w[0]) ;  // f2 = w2w3w0
  f[3]->set_vertex(    i[3] , w[0]) ;  // f3 = w3w4w0
  f[4]->set_vertex(    i[4] , w[6]) ;  // f4 = w4w5w6
  f[6]->set_vertex(    i[6] , w[4]) ;  // f6 = w6w0w4

  nn = f[0]->neighbor( i[0] );
  tds().set_adjacency(f[1], cw(i[1]) , nn, nn->index(f[0]) );
  nn = f[5]->neighbor( i[5] );
  tds().set_adjacency(f[4], ccw(i[4]) , nn, nn->index(f[5]) );
  tds().set_adjacency(f[3], ccw(i[3]) , f[6], ccw(i[6]) );
  tds().set_adjacency(f[6], cw(i[6]) , f[4], cw(i[4]) );

  tds().delete_face(f[0]);
  tds().delete_face(f[5]);

  f[1]->set_offsets(0, 0, 0);
  f[2]->set_offsets(0, 0, 0);
  f[3]->set_offsets(0, 0, 0);
  f[4]->set_offsets(0, 0, 0);
  f[6]->set_offsets(0, 0, 0);
}
template < class Gt, class Tds >
void
Periodic_2_Delaunay_triangulation_2<Gt, Tds>::
remove_degree7_leftfan(Vertex_handle &v, int j,
                       std::vector<Face_handle> &f, std::vector<Vertex_handle> &w, std::vector<Offset> &o, std::vector<int> &i)
{
  // removing a degree 7 vertex, left fan from w[j]

  // Rotate the offset as well
  rotate7(j, w, f, o, i);
  remove_degree7_leftfan(v, /* !! */ 0, f, w, i);

  int oo[7];
  get_offset_degree7(o, oo);

  int o_face[3];
  int ii;
  ii = i[1];
  o_face[ii] = oo[0];
  o_face[ccw(ii)] = oo[1];
  o_face[ cw(ii)] = oo[2];
  this->set_offsets(f[1], o_face[0], o_face[1], o_face[2]);
  ii = i[2];
  o_face[ii] = oo[0];
  o_face[ccw(ii)] = oo[2];
  o_face[ cw(ii)] = oo[3];
  this->set_offsets(f[2], o_face[0], o_face[1], o_face[2]);
  ii = i[3];
  o_face[ii] = oo[0];
  o_face[ccw(ii)] = oo[3];
  o_face[ cw(ii)] = oo[4];
  this->set_offsets(f[3], o_face[0], o_face[1], o_face[2]);
  ii = i[4];
  o_face[ii] = oo[6];
  o_face[ccw(ii)] = oo[4];
  o_face[ cw(ii)] = oo[5];
  this->set_offsets(f[4], o_face[0], o_face[1], o_face[2]);
  ii = i[6];
  o_face[ii] = oo[4];
  o_face[ccw(ii)] = oo[6];
  o_face[ cw(ii)] = oo[0];
  this->set_offsets(f[6], o_face[0], o_face[1], o_face[2]);

  insert_too_long_edge(f[1], ccw(i[1]));
  insert_too_long_edge(f[2], ccw(i[2]));
  insert_too_long_edge(f[3], ccw(i[3]));
  insert_too_long_edge(f[4],  cw(i[4]));
}
template < class Gt, class Tds >
void
Periodic_2_Delaunay_triangulation_2<Gt, Tds>::
remove_degree7_rightfan(Vertex_handle &, int j,
                        std::vector<Face_handle> &f, std::vector<Vertex_handle> &w, std::vector<int> &i)
{
  // removing a degree 7 vertex, right fan from w[j]
  rotate7(j, w, f, i);

  Face_handle nn;
  f[0]->set_vertex(    i[0] , w[3]) ;  // f0 = w0w1w3
  f[2]->set_vertex(    i[2] , w[1]) ;  // f2 = w2w3w1
  f[3]->set_vertex(    i[3] , w[0]) ;  // f3 = w3w4w0
  f[4]->set_vertex(    i[4] , w[0]) ;  // f4 = w4w5w0
  f[5]->set_vertex(    i[5] , w[0]) ;  // f5 = w5w6w0

  nn = f[1]->neighbor( i[1] );
  tds().set_adjacency(f[2], cw(i[2]) , nn, nn->index(f[1]) );
  nn = f[6]->neighbor( i[6] );
  tds().set_adjacency(f[5], ccw(i[5]) , nn, nn->index(f[6]) );
  tds().set_adjacency(f[2], ccw(i[2]) , f[0], ccw(i[0])  );
  tds().set_adjacency(f[0], cw(i[0]) , f[3] , cw(i[3]) );

  tds().delete_face(f[1]);
  tds().delete_face(f[6]);

  f[0]->set_offsets(0, 0, 0);
  f[2]->set_offsets(0, 0, 0);
  f[3]->set_offsets(0, 0, 0);
  f[4]->set_offsets(0, 0, 0);
  f[5]->set_offsets(0, 0, 0);
}
template < class Gt, class Tds >
void
Periodic_2_Delaunay_triangulation_2<Gt, Tds>::
remove_degree7_rightfan(Vertex_handle &v, int j,
                        std::vector<Face_handle> &f, std::vector<Vertex_handle> &w, std::vector<Offset> &o, std::vector<int> &i)
{
  // removing a degree 7 vertex, right fan from w[j]

  // Rotate the offset as well
  rotate7(j, w, f, o, i);
  remove_degree7_rightfan(v, /* !! */ 0, f, w, i);

  int oo[7];
  get_offset_degree7(o, oo);
  int o_face[3];
  int ii;
  ii = i[0];
  o_face[ii] = oo[3];
  o_face[ccw(ii)] = oo[0];
  o_face[ cw(ii)] = oo[1];
  this->set_offsets(f[0], o_face[0], o_face[1], o_face[2]);
  ii = i[2];
  o_face[ii] = oo[1];
  o_face[ccw(ii)] = oo[2];
  o_face[ cw(ii)] = oo[3];
  this->set_offsets(f[2], o_face[0], o_face[1], o_face[2]);
  ii = i[3];
  o_face[ii] = oo[0];
  o_face[ccw(ii)] = oo[3];
  o_face[ cw(ii)] = oo[4];
  this->set_offsets(f[3], o_face[0], o_face[1], o_face[2]);
  ii = i[4];
  o_face[ii] = oo[0];
  o_face[ccw(ii)] = oo[4];
  o_face[ cw(ii)] = oo[5];
  this->set_offsets(f[4], o_face[0], o_face[1], o_face[2]);
  ii = i[5];
  o_face[ii] = oo[0];
  o_face[ccw(ii)] = oo[5];
  o_face[ cw(ii)] = oo[6];
  this->set_offsets(f[5], o_face[0], o_face[1], o_face[2]);

  insert_too_long_edge(f[0], ccw(i[0]));
  insert_too_long_edge(f[0],  cw(i[0]));
  insert_too_long_edge(f[3], ccw(i[3]));
  insert_too_long_edge(f[4], ccw(i[4]));
}

template<class Gt, class Tds>
Oriented_side
Periodic_2_Delaunay_triangulation_2<Gt, Tds>::
side_of_oriented_circle(const Point &p0, const Point &p1, const Point &p2,
                        const Point &p, bool perturb) const
{
  Oriented_side os = geom_traits().side_of_oriented_circle_2_object()(p0, p1, p2, p);
  if ((os != ON_ORIENTED_BOUNDARY) || (!perturb))
    return os;

  // We are now in a degenerate case => we do a symbolic perturbation.

  // We sort the points lexicographically.
  const Point * points[4] = { &p0, &p1, &p2, &p };
  std::sort(points, points + 4, typename Base::Perturbation_order(this));

  // We successively look whether the leading monomial, then 2nd monomial
  // of the determinant has non null coefficient.
  // 2 iterations are enough (cf paper)
  for (int i = 3; i > 0; --i)
    {
      if (points[i] == &p)
        return ON_NEGATIVE_SIDE; // since p0 p1 p2 are non collinear
      // and positively oriented
      Orientation o;
      if (points[i] == &p2 && (o = orientation(p0, p1, p)) != COLLINEAR)
        return Oriented_side(o);
      if (points[i] == &p1 && (o = orientation(p0, p, p2)) != COLLINEAR)
        return Oriented_side(o);
      if (points[i] == &p0 && (o = orientation(p, p1, p2)) != COLLINEAR)
        return Oriented_side(o);
    }
  CGAL_triangulation_assertion(false);
  return ON_NEGATIVE_SIDE;
}

template<class Gt, class Tds>
Oriented_side
Periodic_2_Delaunay_triangulation_2<Gt, Tds>::
side_of_oriented_circle(const Point &p0, const Point &p1, const Point &p2,
                        const Point &p,
                        const Offset &o0, const Offset &o1, const Offset &o2,
                        const Offset &o, bool perturb) const
{
  Oriented_side os = geom_traits().side_of_oriented_circle_2_object()(p0, p1, p2, p, o0, o1, o2, o);
  if ((os != ON_ORIENTED_BOUNDARY) || (!perturb))
    return os;

  // We are now in a degenerate case => we do a symbolic perturbation.
  // We sort the points lexicographically.
  Periodic_point pts[4] = { std::make_pair(p0, o0), std::make_pair(p1, o1),
                            std::make_pair(p2, o2), std::make_pair(p, o)
                          };
  const Periodic_point *points[4] = { &pts[0], &pts[1], &pts[2], &pts[3] };

  std::sort(points, points + 4, typename Base::Perturbation_order(this));

  // We successively look whether the leading monomial, then 2nd monomial
  // of the determinant has non null coefficient.
  // 2 iterations are enough (cf paper)
  for (int i = 3; i > 0; --i)
    {
      if (points[i] == &pts[3])
        return ON_NEGATIVE_SIDE; // since p0 p1 p2 are non collinear
      // and positively oriented
      Orientation orient;
      if ((points[i] == &pts[2]) && ((orient = orientation(p0, p1, p, o0, o1, o))
                                     != COLLINEAR))
        return Oriented_side(orient);
      if ((points[i] == &pts[1]) && ((orient = orientation(p0, p, p2, o0, o, o2))
                                     != COLLINEAR))
        return Oriented_side(orient);
      if ((points[i] == &pts[0]) && ((orient = orientation(p, p1, p2, o, o1, o2))
                                     != COLLINEAR))
        return Oriented_side(orient);
    }
  CGAL_triangulation_assertion(false);
  return ON_NEGATIVE_SIDE;
}

template<class Gt, class Tds>
Oriented_side
Periodic_2_Delaunay_triangulation_2<Gt, Tds>::
side_of_oriented_circle(Face_handle f, const Point & p, bool perturb) const
{
  Oriented_side os = ON_NEGATIVE_SIDE;

  int i = 0;
  // TODO: optimize which copies to check depending on the offsets in
  // the face.
  while (os == ON_NEGATIVE_SIDE && i < 4)
    {
      os = side_of_oriented_circle(f->vertex(0)->point(), f->vertex(1)->point(), f->vertex(2)->point(), p,
                                   get_offset(f, 0), get_offset(f, 1), get_offset(f, 2), combine_offsets(Offset(), int_to_off(i)),
                                   perturb);
      i++;
    }

  return os;
}

///////////////////////////////////////////////////////////////
//  DISPLACEMENT

template <class Gt, class Tds >
typename Periodic_2_Delaunay_triangulation_2<Gt, Tds>::Vertex_handle
Periodic_2_Delaunay_triangulation_2<Gt, Tds>::
move_if_no_collision(Vertex_handle v, const Point &p)
{
  Locate_type lt;
  int li;
  Vertex_handle inserted;
  Face_handle loc = locate(p, lt, li, v->face());

  if (lt == Base::VERTEX)
    return v;
  else
    /// This can be optimized by checking whether we can move v->point() to p
    return insert(p, lt, loc, li);
}

template <class Gt, class Tds >
typename Periodic_2_Delaunay_triangulation_2<Gt, Tds>::Vertex_handle
Periodic_2_Delaunay_triangulation_2<Gt, Tds>::
move_point(Vertex_handle v, const Point &p)
{
  if(v->point() == p) return v;
  Vertex_handle w = move_if_no_collision(v, p);
  if(w != v)
    {
      remove(v);
      return w;
    }
  return v;
}

template<class Gt, class Tds>
void Periodic_2_Delaunay_triangulation_2<Gt, Tds>::fill_hole_delaunay(std::list<Edge> & first_hole)
{
  typedef std::list<Edge> Hole;
  typedef std::list<Hole> Hole_list;

  Face_handle  f, ff, fn;
  int i, ii, in;
  Hole_list hole_list;
  Hole hole;

  hole_list.push_front(first_hole);

  while( ! hole_list.empty())
    {
      hole = hole_list.front();
      hole_list.pop_front();
      typename Hole::iterator hit = hole.begin();

      // if the hole has only three edges, create the triangle
      if (hole.size() == 3)
        {
          hit = hole.begin();
          f = (*hit).first;
          i = (*hit).second;
          ff = (* ++hit).first;
          ii = (*hit).second;
          fn = (* ++hit).first;
          in = (*hit).second;
          Face_handle newf = create_face(f, i, ff, ii, fn, in);
          newf->set_offsets(0, 0, 0);

          continue;
        }

      // else find an edge with two finite vertices
      // on the hole boundary
      // and the new triangle adjacent to that edge
      //  cut the hole and push it back

      // take the first neighboring face and pop it;
      ff = (hole.front()).first;
      ii = (hole.front()).second;
      hole.pop_front();

      Vertex_handle v0 = ff->vertex(cw(ii));
      Vertex_handle v1 = ff->vertex(ccw(ii));
      Vertex_handle v2 = Vertex_handle();
      Vertex_handle v3 = Vertex_handle();
      const Point& p0 = v0->point();
      const Point& p1 = v1->point();

      typename Hole::iterator hdone = hole.end();
      hit =  hole.begin();
      typename Hole::iterator cut_after(hit);

      // if tested vertex is c with respect to the vertex opposite
      // to NULL neighbor,
      // stop at the before last face;
      hdone--;
      while( hit != hdone)
        {
          fn = (*hit).first;
          in = (*hit).second;
          Vertex_handle vv = fn->vertex(ccw(in));

          const Point &p = vv->point();
          Orientation orient = orientation(p0, p1, p);

          if (orient == COUNTERCLOCKWISE)
            {
              if (v2 == Vertex_handle())
                {
                  v2 = vv;
                  v3 = vv;
                  cut_after = hit;
                }
              else
                {
                  Oriented_side side = side_of_oriented_circle(p0, p1, v3->point(), p, true);
                  if (side == ON_POSITIVE_SIDE)
                    {
                      v2 = vv;
                      v3 = vv;
                      cut_after = hit;
                    }
                }
            }

          ++hit;
        }

      // create new triangle and update adjacency relations
      Face_handle newf;

      //update the hole and push back in the Hole_List stack
      // if v2 belongs to the neighbor following or preceding *f
      // the hole remain a single hole
      // otherwise it is split in two holes

      fn = (hole.front()).first;
      in = (hole.front()).second;
      if (fn->has_vertex(v2, i) && i == fn->ccw(in))
        {
          newf = create_face(ff, ii, fn, in);

          newf->set_offsets(0, 0, 0);

          hole.pop_front();
          hole.push_front(Edge(newf, 1));
          hole_list.push_front(hole);
        }
      else
        {
          fn = (hole.back()).first;
          in = (hole.back()).second;
          if (fn->has_vertex(v2, i) && i == fn->cw(in))
            {
              newf = create_face(fn, in, ff, ii);
              newf->set_offsets(0, 0, 0);

              hole.pop_back();
              hole.push_back(Edge(newf, 1));
              hole_list.push_front(hole);
            }
          else
            {
              // split the hole in two holes
              CGAL_assertion(v2 != Vertex_handle());
              newf = create_face(ff, ii, v2);
              newf->set_offsets(0, 0, 0);

              Hole new_hole;
              ++cut_after;
              while( hole.begin() != cut_after )
                {
                  new_hole.push_back(hole.front());
                  hole.pop_front();
                }

              hole.push_front(Edge( newf, 1));
              new_hole.push_front(Edge( newf, 0));
              hole_list.push_front(hole);
              hole_list.push_front(new_hole);
            }
        }
    }
}

template<class Gt, class Tds>
void Periodic_2_Delaunay_triangulation_2<Gt, Tds>::fill_hole_delaunay(
  std::list<Edge> & first_hole,
  std::map<Vertex_handle, Offset> &vertex_offsets)
{
  typedef std::list<Edge> Hole;
  typedef std::list<Hole> Hole_list;

  Face_handle  f, ff, fn;
  int i, ii, in;
  Hole_list hole_list;
  Hole hole;

  hole_list.push_front(first_hole);

  while( ! hole_list.empty())
    {
      hole = hole_list.front();
      hole_list.pop_front();
      typename Hole::iterator hit = hole.begin();

      // if the hole has only three edges, create the triangle
      if (hole.size() == 3)
        {
          hit = hole.begin();
          f = (*hit).first;
          i = (*hit).second;
          ff = (* ++hit).first;
          ii = (*hit).second;
          fn = (* ++hit).first;
          in = (*hit).second;
          Face_handle newf = create_face(f, i, ff, ii, fn, in);
          Offset oo0(vertex_offsets[newf->vertex(0)]);
          Offset oo1(vertex_offsets[newf->vertex(1)]);
          Offset oo2(vertex_offsets[newf->vertex(2)]);
          if (oo0.x() < 0 || oo1.x() < 0 || oo2.x() < 0)
            {
              oo0 += Offset(number_of_sheets()[0], 0);
              oo1 += Offset(number_of_sheets()[0], 0);
              oo2 += Offset(number_of_sheets()[0], 0);
            }
          if (oo0.y() < 0 || oo1.y() < 0 || oo2.y() < 0)
            {
              oo0 += Offset(0, number_of_sheets()[1]);
              oo1 += Offset(0, number_of_sheets()[1]);
              oo2 += Offset(0, number_of_sheets()[1]);
            }
          set_offsets(newf,
                      (oo0.x() >= number_of_sheets()[0] ? 2 : 0) + (oo0.y() >= number_of_sheets()[1] ? 1 : 0),
                      (oo1.x() >= number_of_sheets()[0] ? 2 : 0) + (oo1.y() >= number_of_sheets()[1] ? 1 : 0),
                      (oo2.x() >= number_of_sheets()[0] ? 2 : 0) + (oo2.y() >= number_of_sheets()[1] ? 1 : 0));

          insert_too_long_edge(newf, 0);
          insert_too_long_edge(newf, 1);
          insert_too_long_edge(newf, 2);

          continue;
        }

      // else find an edge with two finite vertices
      // on the hole boundary
      // and the new triangle adjacent to that edge
      //  cut the hole and push it back

      // take the first neighboring face and pop it;
      ff = (hole.front()).first;
      ii = (hole.front()).second;
      hole.pop_front();

      Vertex_handle v0 = ff->vertex(cw(ii));
      Vertex_handle v1 = ff->vertex(ccw(ii));
      Vertex_handle v2 = Vertex_handle();
      Vertex_handle v3 = Vertex_handle();
      const Point& p0 = v0->point();
      const Point& p1 = v1->point();
      const Offset o0 = vertex_offsets[v0];
      const Offset o1 = vertex_offsets[v1];
      bool simplicity_criterion = (o0 == o1);

      typename Hole::iterator hdone = hole.end();
      hit =  hole.begin();
      typename Hole::iterator cut_after(hit);

      // if tested vertex is c with respect to the vertex opposite
      // to NULL neighbor,
      // stop at the before last face;
      hdone--;
      while( hit != hdone)
        {
          fn = (*hit).first;
          in = (*hit).second;
          Vertex_handle vv = fn->vertex(ccw(in));

          const Point &p = vv->point();
          CGAL_assertion(vertex_offsets.find(vv) != vertex_offsets.end());
          const Offset o = vertex_offsets[vv];
          Orientation orient;
          simplicity_criterion &= (o == o0);
          if (simplicity_criterion)
            orient = orientation(p0, p1, p);
          else
            orient = orientation(p0, p1, p, o0, o1, o);

          if (orient == COUNTERCLOCKWISE)
            {
              if (v2 == Vertex_handle())
                {
                  v2 = vv;
                  v3 = vv;
                  cut_after = hit;
                }
              else
                {
                  Offset o3 = vertex_offsets[v3];
                  Oriented_side side;
                  if (simplicity_criterion && (o3 == o0))
                    side = side_of_oriented_circle(p0, p1, v3->point(), p,
                                                   true);
                  else
                    side = side_of_oriented_circle(p0, p1, v3->point(), p,
                                                   o0, o1, o3, o,
                                                   true);
                  if (side == ON_POSITIVE_SIDE)
                    {
                      v2 = vv;
                      v3 = vv;
                      cut_after = hit;
                    }
                }
            }

          ++hit;
        }

      // create new triangle and update adjacency relations
      Face_handle newf;

      //update the hole and push back in the Hole_List stack
      // if v2 belongs to the neighbor following or preceding *f
      // the hole remain a single hole
      // otherwise it is split in two holes

      fn = (hole.front()).first;
      in = (hole.front()).second;
      if (fn->has_vertex(v2, i) && i == fn->ccw(in))
        {
          newf = create_face(ff, ii, fn, in);

          Offset oo0 = o0;
          Offset oo1 = o1;
          Offset oo2 = vertex_offsets[v2];
          if (oo0.x() < 0 || oo1.x() < 0 || oo2.x() < 0)
            {
              oo0 += Offset(number_of_sheets()[0], 0);
              oo1 += Offset(number_of_sheets()[0], 0);
              oo2 += Offset(number_of_sheets()[0], 0);
            }
          if (oo0.y() < 0 || oo1.y() < 0 || oo2.y() < 0)
            {
              oo0 += Offset(0, number_of_sheets()[1]);
              oo1 += Offset(0, number_of_sheets()[1]);
              oo2 += Offset(0, number_of_sheets()[1]);
            }
          set_offsets(newf,
                      (oo0.x() >= number_of_sheets()[0] ? 2 : 0) + (oo0.y() >= number_of_sheets()[1] ? 1 : 0),
                      (oo1.x() >= number_of_sheets()[0] ? 2 : 0) + (oo1.y() >= number_of_sheets()[1] ? 1 : 0),
                      (oo2.x() >= number_of_sheets()[0] ? 2 : 0) + (oo2.y() >= number_of_sheets()[1] ? 1 : 0));
          // set_offsets(newf, o0, o1, o2);
          insert_too_long_edge(newf, 0);
          insert_too_long_edge(newf, 1);

          hole.pop_front();
          hole.push_front(Edge(newf, 1));
          hole_list.push_front(hole);
        }
      else
        {
          fn = (hole.back()).first;
          in = (hole.back()).second;
          if (fn->has_vertex(v2, i) && i == fn->cw(in))
            {
              newf = create_face(fn, in, ff, ii);
              Offset oo0 = o0;
              Offset oo1 = o1;
              Offset oo2 = vertex_offsets[v2];
              if (oo0.x() < 0 || oo1.x() < 0 || oo2.x() < 0)
                {
                  oo0 += Offset(number_of_sheets()[0], 0);
                  oo1 += Offset(number_of_sheets()[0], 0);
                  oo2 += Offset(number_of_sheets()[0], 0);
                }
              if (oo0.y() < 0 || oo1.y() < 0 || oo2.y() < 0)
                {
                  oo0 += Offset(0, number_of_sheets()[1]);
                  oo1 += Offset(0, number_of_sheets()[1]);
                  oo2 += Offset(0, number_of_sheets()[1]);
                }
              set_offsets(newf,
                          (oo2.x() >= number_of_sheets()[0] ? 2 : 0) + (oo2.y() >= number_of_sheets()[1] ? 1 : 0),
                          (oo0.x() >= number_of_sheets()[0] ? 2 : 0) + (oo0.y() >= number_of_sheets()[1] ? 1 : 0),
                          (oo1.x() >= number_of_sheets()[0] ? 2 : 0) + (oo1.y() >= number_of_sheets()[1] ? 1 : 0));
              insert_too_long_edge(newf, 1);
              insert_too_long_edge(newf, 2);
              hole.pop_back();
              hole.push_back(Edge(newf, 1));
              hole_list.push_front(hole);
            }
          else
            {
              // split the hole in two holes
              CGAL_assertion(v2 != Vertex_handle());
              newf = create_face(ff, ii, v2);
              Offset oo0 = o0;
              Offset oo1 = o1;
              Offset oo2 = vertex_offsets[v2];
              if (oo0.x() < 0 || oo1.x() < 0 || oo2.x() < 0)
                {
                  oo0 += Offset(number_of_sheets()[0], 0);
                  oo1 += Offset(number_of_sheets()[0], 0);
                  oo2 += Offset(number_of_sheets()[0], 0);
                }
              if (oo0.y() < 0 || oo1.y() < 0 || oo2.y() < 0)
                {
                  oo0 += Offset(0, number_of_sheets()[1]);
                  oo1 += Offset(0, number_of_sheets()[1]);
                  oo2 += Offset(0, number_of_sheets()[1]);
                }
              set_offsets(newf,
                          (oo0.x() >= number_of_sheets()[0] ? 2 : 0) + (oo0.y() >= number_of_sheets()[1] ? 1 : 0),
                          (oo1.x() >= number_of_sheets()[0] ? 2 : 0) + (oo1.y() >= number_of_sheets()[1] ? 1 : 0),
                          (oo2.x() >= number_of_sheets()[0] ? 2 : 0) + (oo2.y() >= number_of_sheets()[1] ? 1 : 0));


              // set_offsets(newf, o0, o1, o2);
              insert_too_long_edge(newf, 0);
              insert_too_long_edge(newf, 1);

              Hole new_hole;
              ++cut_after;
              while( hole.begin() != cut_after )
                {
                  new_hole.push_back(hole.front());
                  hole.pop_front();
                }

              hole.push_front(Edge( newf, 1));
              new_hole.push_front(Edge( newf, 0));
              hole_list.push_front(hole);
              hole_list.push_front(new_hole);
            }
        }
    }
}
} //namespace CGAL

#endif // CGAL_PERIODIC_2_DELAUNAY_TRIANGULATION_2_H
