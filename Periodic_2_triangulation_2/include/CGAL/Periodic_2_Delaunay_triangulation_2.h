// Copyright (c) 1997-2013 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Nico Kruithof <Nico@nghk.nl>

#ifndef CGAL_PERIODIC_2_DELAUNAY_TRIANGULATION_2_H
#define CGAL_PERIODIC_2_DELAUNAY_TRIANGULATION_2_H

#include <CGAL/license/Periodic_2_triangulation_2.h>

#include <CGAL/Periodic_2_triangulation_2.h>
#include <CGAL/Periodic_2_triangulation_vertex_base_2.h>
#include <CGAL/Periodic_2_triangulation_face_base_2.h>

#include <CGAL/algorithm.h>
#include <CGAL/iterator.h>

#ifndef CGAL_TRIANGULATION_2_DONT_INSERT_RANGE_OF_POINTS_WITH_INFO
#include <CGAL/Spatial_sort_traits_adapter_2.h>
#include <CGAL/internal/info_check.h>
#include <CGAL/tss.h>

#include <boost/iterator/zip_iterator.hpp>
#include <boost/mpl/and.hpp>
#endif //CGAL_TRIANGULATION_2_DONT_INSERT_RANGE_OF_POINTS_WITH_INFO

namespace CGAL {

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
  typedef Tds                                                   Triangulation_data_structure;
  typedef Gt                                                    Geom_traits;

  typedef typename Gt::Periodic_2_offset_2                      Offset;
  typedef typename Gt::Domain                                   Domain;
  typedef std::array<int, 2>                                    Covering_sheets;

  typedef typename Gt::FT                                       FT;
  typedef typename Gt::Point_2                                  Point;
  typedef typename Gt::Segment_2                                Segment;
  typedef typename Gt::Triangle_2                               Triangle;

  typedef std::pair<Point, Offset>                              Periodic_point;
  typedef std::array< std::pair<Point, Offset>, 2>              Periodic_segment;
  typedef std::array< std::pair<Point, Offset>, 3>              Periodic_triangle;
  typedef std::array< std::pair<Point, Offset>, 4>              Periodic_tetrahedron;

  typedef typename Base::size_type                              size_type;
  typedef typename Base::Locate_type                            Locate_type;

  typedef typename Base::Vertex_handle                          Vertex_handle;
  typedef typename Base::Vertex_circulator                      Vertex_circulator;
  typedef typename Base::Vertex_iterator                        Vertex_iterator;
  typedef typename Base::Edge                                   Edge;
  typedef typename Base::Edge_circulator                        Edge_circulator;
  typedef typename Base::Edge_iterator                          Edge_iterator;
  typedef typename Base::Face_handle                            Face_handle;
  typedef typename Base::Face_circulator                        Face_circulator;
  typedef typename Base::Face_iterator                          Face_iterator;
  typedef typename Base::Finite_vertices_iterator               Finite_vertices_iterator;
  typedef typename Base::Finite_edges_iterator                  Finite_edges_iterator;
  typedef typename Base::Finite_faces_iterator                  Finite_faces_iterator;
  typedef typename Base::All_faces_iterator                     All_faces_iterator;

  typedef typename Base::Periodic_segment_iterator              Periodic_segment_iterator;
  typedef typename Base::Periodic_triangle_iterator             Periodic_triangle_iterator;

  //Tag to distinguish Delaunay from regular triangulations
  typedef Tag_false                                             Weighted_tag;

  // Tag to distinguish periodic triangulations from others
  typedef Tag_true                                              Periodic_tag;

private:
  struct Face_handle_hash
      : public CGAL::cpp98::unary_function<Face_handle, std::size_t>
  {
    std::size_t operator()(Face_handle fh) const
    {
      return boost::hash<typename Face_handle::pointer>()(&*fh);
    }
  };

  typedef std::unordered_set<Face_handle, Face_handle_hash>     Too_big_circumdisks_set;
  typedef typename Too_big_circumdisks_set::const_iterator      Too_big_circumdisks_set_it;

public:
#ifndef CGAL_CFG_USING_BASE_MEMBER_BUG_2
  using Base::empty;
  using Base::cw;
  using Base::ccw;
  using Base::locate;
  using Base::construct_point;
  using Base::periodic_point;

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
  using Base::number_of_vertices;
  using Base::vertices_begin;
  using Base::vertices_end;
  using Base::edges_begin;
  using Base::edges_end;
  using Base::faces_begin;
  using Base::faces_end;
  using Base::incident_faces;

  using Base::orientation;
  using Base::point;
  using Base::construct_segment;
#endif

public:
  /// Checks whether f->vertex(i) lies outside the circumcircle of the face nb
  inline bool locally_Delaunay(Face_handle f, int i, Face_handle nb);

  /// Determines whether the point p lies on the (un-)bounded side of
  /// the circle through the points p0, p1 and p2
  Oriented_side side_of_oriented_circle(const Point& p0, const Point& p1, const Point& p2, const Point& p,
                                        bool perturb) const;

  /// Determines whether the point (p,o) lies on the (un-)bounded side of
  /// the circle through the points (p0,o0), (p1,o1) and (p2,o2)
  Oriented_side side_of_oriented_circle(const Point& p0, const Point& p1, const Point& p2, const Point& p,
                                        const Offset& o0, const Offset& o1, const Offset& o2, const Offset& o,
                                        bool perturb) const;

  /// Determines whether the point p lies on the (un-)bounded side of
  /// the circle through the vertices of f
  Oriented_side side_of_oriented_circle(Face_handle f, const Point& p, bool perturb = false) const;

  Bounded_side side_of_circle(Face_handle f,
                              const Point& p, const Offset& o = Offset(),
                              bool perturb = false) const
  {
    return enum_cast<Bounded_side>(
             side_of_oriented_circle(
               f->vertex(0)->point(), f->vertex(1)->point(), f->vertex(2)->point(), p,
               this->get_offset(f, 0), this->get_offset(f, 1), this->get_offset(f, 2), o,
               perturb));
  }

  bool incircle(int x, int j, int k, int l,
                 std::vector<Face_handle>&,
                 std::vector<Vertex_handle>& w,
                 std::vector<int>&)
  {
    return _side_of_circle(w[j]->point(), w[k]->point(), w[l]->point(), w[x]->point(), true) ==  ON_POSITIVE_SIDE;
  }

  bool incircle(int x, int j, int k, int l,
                 std::vector<Face_handle> &,
                 std::vector<Vertex_handle> &w,
                 std::vector<Offset> &o,
                 std::vector<int> &)
  {
    return _side_of_circle(w[j]->point(), w[k]->point(), w[l]->point(), w[x]->point(),
                           o[j], o[k], o[l], o[x], true) ==  ON_POSITIVE_SIDE;
  }

  CGAL::Comparison_result
  compare_squared_circumradius_to_threshold(const Periodic_point& p0, const Periodic_point& p1,
                                            const Periodic_point& p2, const FT threshold) const
  {
    return geom_traits().compare_squared_radius_2_object()(p0.first, p1.first, p2.first,
                                                           p0.second, p1.second, p2.second,
                                                           threshold);
  }

  CGAL::Comparison_result
  compare_squared_circumradius_to_threshold(Face_handle face, const FT threshold) const
  {
    Periodic_point p0 = periodic_point(face, 0);
    Periodic_point p1 = periodic_point(face, 1);
    Periodic_point p2 = periodic_point(face, 2);

    return compare_squared_circumradius_to_threshold(p0, p1, p2, threshold);
  }

  /// Constructs the circumcenter of the face f, respects the offset
  Point circumcenter(Face_handle f) const
  {
    return construct_circumcenter(f->vertex(0)->point(), f->vertex(1)->point(), f->vertex(2)->point(),
                                  get_offset(f, 0), get_offset(f, 1), get_offset(f, 2));
  }
  Point construct_circumcenter(const Point& p1, const Point& p2, const Point& p3,
                               const Offset& o1, const Offset& o2, const Offset& o3) const
  {
    return geom_traits().construct_circumcenter_2_object()(p1, p2, p3, o1, o2, o3);
  }

public:
  /// \name Constructors
  Periodic_2_Delaunay_triangulation_2(const Gt& gt)
    : Base(gt)
  {
    update_cover_data_after_setting_domain();
  }

  Periodic_2_Delaunay_triangulation_2(const Domain& domain = Domain())
    : Periodic_2_Delaunay_triangulation_2(Gt(domain))
  { }

  /// Copy
  // @todo (can't be "= default" because some members are pointers)
  Periodic_2_Delaunay_triangulation_2(const Periodic_2_Delaunay_triangulation_2<Gt, Tds>& tr) = delete;
  Periodic_2_Delaunay_triangulation_2& operator=(const Periodic_2_Delaunay_triangulation_2&) = delete;

  /// Constructor with insertion of points
  template <class InputIterator>
  Periodic_2_Delaunay_triangulation_2(InputIterator first, InputIterator last,
                                      const Gt& gt = Gt())
    : Periodic_2_triangulation_2<Gt, Tds>(domain, gt)
  {
    insert(first, last);
  }

  void copy_multiple_covering(const Periodic_2_Delaunay_triangulation_2& tr);

  void swap(Periodic_2_Delaunay_triangulation_2& tr)
  {
    Base::swap(tr);

    std::swap(squared_circumradius_threshold, tr.squared_circumradius_threshold);
    std::swap(faces_with_too_big_circumdisk, tr.faces_with_too_big_circumdisk);
  }

  void clear()
  {
    Base::clear();
    faces_with_too_big_circumdisk.clear();
  }

  void update_cover_data_after_setting_domain()
  {
#ifndef CGAL_GENERIC_P2T2
    // the criterion is that the largest circumdisk must have a diameter smaller than c/2
    // (c being the square side), thus we need a squared circumdisk radius smaller than c*c/16
    squared_circumradius_threshold = (domain().xmax() - domain().xmin()) *
                                     (domain().xmax() - domain().xmin()) / FT(16);
#endif
  }

  void set_domain(const Domain& domain) override
  {
    clear();
    Base::set_domain(domain);
    update_cover_data_after_setting_domain();
  }

  /// \name Insertion-Removal
  Vertex_handle insert(const Point& p, Face_handle start = Face_handle())
  {
    Conflict_tester tester(p, this);
    Point_hider hider;
    Cover_manager cover_manager(*this);
    Vertex_handle vh = Base::insert_in_conflict(p, start, tester, hider, cover_manager);
    CGAL_assertion(vh != Vertex_handle());
    return vh;
  }

  Vertex_handle insert(const Point& p,
                       Locate_type lt,
                       Face_handle f,
                       int li)
  {
    Conflict_tester tester(p, this);
    Point_hider hider;
    Cover_manager cover_manager(*this);
    Vertex_handle vh = Base::insert_in_conflict(p, lt, f, li, tester, hider, cover_manager);
    CGAL_assertion(vh != Vertex_handle());
    return vh;
  }

  Vertex_handle push_back(const Point& p) { return insert(p); }

  /// Insertion with info

#if 0 // @todo need to introduce Periodic_2_Delaunay_triangulation_remove_traits_2, see 3D version

#ifndef CGAL_TRIANGULATION_2_DONT_INSERT_RANGE_OF_POINTS_WITH_INFO
  template <class InputIterator>
  std::ptrdiff_t
  insert(InputIterator first, InputIterator last,
         bool is_large_point_set = true,
         typename boost::enable_if<
                    boost::is_convertible<
                      typename std::iterator_traits<InputIterator>::value_type, Point> >::type* = nullptr)
#else
  template <class InputIterator>
  std::ptrdiff_t
  insert(InputIterator first, InputIterator last,
         bool is_large_point_set = true)
#endif //CGAL_TRIANGULATION_2_DONT_INSERT_RANGE_OF_POINTS_WITH_INFO
  {
    if(first == last)
      return 0;

    size_type n = number_of_vertices();

    // The heuristic discards the existing triangulation so it can only be
    // applied to empty triangulations.
    if(n != 0)
      is_large_point_set = false;

    std::vector<Point> points(first, last);
    std::vector<Vertex_handle> dummy_points;
    typename std::vector<Point>::iterator pbegin = points.begin();

    if(is_large_point_set)
    {
      std::vector<Vertex_handle> dummy_points = this->insert_dummy_points();
    }
    else
    {
      CGAL::cpp98::random_shuffle(points.begin(), points.end());
      pbegin = points.begin();

      for(;;)
      {
        if(pbegin == points.end())
          return number_of_vertices() - n;

        insert(*pbegin);
        ++pbegin;

        if(is_1_cover())
          break;
      }
    }

    CGAL_assertion(is_1_cover());

    // Organize the points
    spatial_sort(pbegin, points.end(), geom_traits());

    Face_handle hint;
    Conflict_tester tester(*pbegin, this);
    Point_hider hider;
    Cover_manager cover_manager(*this);

    // Actual insertion
    std::vector<Vertex_handle> double_vertices =
      Base::insert_in_conflict(points.begin(), points.end(), hint, tester, hider, cover_manager);

    CGAL_assertion_code(for(Vertex_handle vh : double_vertices))
    CGAL_assertion(vh != Vertex_handle());

    if(is_large_point_set)
    {
      typedef CGAL::Periodic_2_Delaunay_triangulation_remove_traits_2<Gt> P2removeT;
      typedef CGAL::Delaunay_triangulation_2<P2removeT> DT;
      typedef Vertex_remover<DT> Remover;

      P2removeT remove_traits(domain());
      DT dt(remove_traits);
      Remover remover(this, dt);
      Conflict_tester t(this);

      for(const Vertex_handle dummy_vh : dummy_points)
      {
        if(std::find(double_vertices.begin(), double_vertices.end(), dummy_vh) == double_vertices.end())
          Base::remove(dummy_vh, remover, t, cover_manager);
      }
    }

    return number_of_vertices() - n;
  }
#endif

#if 0
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
  std::ptrdiff_t insert_with_info(InputIterator first, InputIterator last,
                                  bool is_large_point_set)
  {
    if(first == last)
      return 0;

    std::vector<std::size_t> indices;
    std::vector<Point> points;
    std::vector<typename Tds::Vertex::Info> infos;
    std::size_t index = 0;
    for(InputIterator it = first; it != last; ++it)
    {
      Tuple_or_pair value = *it;
      points.push_back(top_get_first(value));
      infos.push_back(top_get_second(value));
      indices.push_back(index++);
    }

    typedef typename Pointer_property_map<Point>::type Pmap;
    typedef Spatial_sort_traits_adapter_2<Geom_traits, Pmap> Search_traits;

    size_type n = number_of_vertices();

    // The heuristic discards the existing triangulation so it can only be
    // applied to empty triangulations.
    if(n != 0)
      is_large_point_set = false;

    std::set<Vertex_handle> dummy_points;
    typename std::vector<std::size_t>::iterator ind_it = indices.begin();

    if(is_large_point_set)
    {
      std::vector<Vertex_handle> dummy_points_vector = this->insert_dummy_points();
      std::copy(dummy_points_vector.begin(), dummy_points_vector.end(),
                std::inserter(dummy_points, dummy_points.begin()));
    }
    else
    {
      CGAL::cpp98::random_shuffle(indices.begin(), indices.end());
      ind_it = indices.begin();

      for(;;)
      {
        if(ind_it == indices.end())
          return number_of_vertices() - n;

        Vertex_handle v_new = insert(points[*ind_it]);
        v_new->info() = infos[*ind_it];
        ++ind_it;

        if(is_1_cover())
          break;
      }
    }

    CGAL_assertion(is_1_cover());

    // Insert the points
    spatial_sort(indices.begin(), indices.end(),
                 Search_traits(make_property_map(points), geom_traits()));

    Face_handle f;
    for(typename std::vector<std::size_t>::const_iterator it=ind_it, end=indices.end();
        it != end; ++it)
    {
      Locate_type lt;
      int li, lj;
      Offset o;
      f = locate(points[*it], o, lt, li, lj, f);

      if(lt == Base::VERTEX)
      {
        // Always copy the info, it might be a dummy vertex
        f->vertex(li)->info() = infos[*it];
        dummy_points.erase(f->vertex(li));
      }
      else
      {
        Vertex_handle v_new = insert(points[*it], o, lt, f, li); // @fixme conflict_finder & stuff
        v_new->info() = infos[*it];
      }
    }

    if(is_large_point_set)
    {
      typedef CGAL::Periodic_2_Delaunay_triangulation_remove_traits_2<Gt> P2removeT;
      typedef CGAL::Delaunay_triangulation_2<P2removeT> DT;
      typedef Vertex_remover<DT> Remover;

      P2removeT remove_traits(domain());
      DT dt(remove_traits);
      Remover remover(this, dt);
      Conflict_tester t(this);

      for(const Vertex_handle dummy_vh : dummy_points)
      {
        if(std::find(double_vertices.begin(), double_vertices.end(), dummy_vh) == double_vertices.end())
          Base::remove(dummy_vh, remover, t, cover_manager);
      }
    }

    return number_of_vertices() - n;
  }

public:
  template <class InputIterator>
  std::ptrdiff_t
  insert(InputIterator first, InputIterator last,
         bool is_large_point_set = true,
         typename boost::enable_if<
                    boost::is_convertible<
                      typename std::iterator_traits<InputIterator>::value_type,
                      std::pair<Point, typename internal::Info_check<typename Tds::Vertex>::type> > >::type* = nullptr)
  {
    return insert_with_info<std::pair<Point, typename internal::Info_check<typename Tds::Vertex>::type> >(first, last, is_large_point_set);
  }

  template <class  InputIterator_1, class InputIterator_2>
  std::ptrdiff_t
  insert(boost::zip_iterator<boost::tuple<InputIterator_1, InputIterator_2> > first,
         boost::zip_iterator<boost::tuple<InputIterator_1, InputIterator_2> > last,
         bool is_large_point_set = true,
         typename boost::enable_if<
                    boost::mpl::and_<
                      boost::is_convertible<
                        typename std::iterator_traits<InputIterator_1>::value_type, Point>,
                      boost::is_convertible<
                        typename std::iterator_traits<InputIterator_2>::value_type,
                        typename internal::Info_check<typename Tds::Vertex>::type> > >::type* = nullptr)
  {
    return insert_with_info<boost::tuple<Point, typename internal::Info_check<typename Tds::Vertex>::type> >(first, last, is_large_point_set);
  }
#endif //CGAL_TRIANGULATION_2_DONT_INSERT_RANGE_OF_POINTS_WITH_INFO
#endif

  void remove(Vertex_handle v);

  template <typename InputIterator>
  std::ptrdiff_t remove(InputIterator first, InputIterator beyond)
  {
    std::size_t n = number_of_vertices();
    while(first != beyond)
      remove(*first++);

    return n - number_of_vertices();
  }

  /// \name Displacement

  Vertex_handle move_if_no_collision(Vertex_handle v, const Point& p);
  Vertex_handle move_point(Vertex_handle v, const Point& p);

  /// \name Conflict checking

  template <class OutputIteratorBoundaryEdges,
            class OutputIteratorFaces,
            class OutputIteratorInternalEdges>
  Triple<OutputIteratorBoundaryEdges, OutputIteratorFaces, OutputIteratorInternalEdges>
  find_conflicts(const Point& p,
                 Face_handle f,
                 OutputIteratorBoundaryEdges beit,
                 OutputIteratorFaces fit,
                 OutputIteratorInternalEdges ieit) const
  {
    CGAL_triangulation_precondition(number_of_vertices() != 0);

    std::vector<Edge> edges;
    edges.reserve(16);
    std::vector<Face_handle> faces;
    faces.reserve(16);

    Conflict_tester tester(p, this);

    Triple<typename std::back_insert_iterator<std::vector<Edge> >,
           typename std::back_insert_iterator<std::vector<Face_handle> >,
           OutputIteratorInternalEdges> tit =
        Base::find_conflicts(f, tester, make_triple(std::back_inserter(edges),
                                                    std::back_inserter(faces),
                                                    ieit));
    ieit = tit.third;

    // Reset the conflict flag on the boundary.
    for(const Edge& e : edges)
    {
      e.first->neighbor(e.second)->tds_data().clear();
      *beit++ = e;
    }

    // Reset the conflict flag in the conflict faces.
    for(Face_handle fh : faces)
    {
      fh->tds_data().clear();
      *fit++ = fh;
    }

    for(Vertex_handle vc : this->v_offsets)
      vc->clear_offset();
    this->v_offsets.clear();

    return make_triple(beit, fit, ieit);
  }

  template <class OutputIteratorBoundaryEdges,
            class OutputIteratorFaces>
  std::pair<OutputIteratorBoundaryEdges, OutputIteratorFaces>
  find_conflicts(const Point &p,
                 Face_handle f,
                 OutputIteratorBoundaryEdges beit,
                 OutputIteratorFaces fit) const
  {
    Triple<OutputIteratorBoundaryEdges, OutputIteratorFaces, Emptyset_iterator> t =
      find_conflicts(p, f, beit, fit, Emptyset_iterator());

    return std::make_pair(t.first, t.second);
  }

  template <class OutputIteratorFaces>
  OutputIteratorFaces
  find_conflicts(const Point &p,
                 Face_handle f,
                 OutputIteratorFaces fit) const
  {
    Triple<Emptyset_iterator, OutputIteratorFaces, Emptyset_iterator> t =
      find_conflicts(p, f, Emptyset_iterator(), fit, Emptyset_iterator());

    return t.second;
  }

public:
  class Point_hider
  {
  public:
    template <class InputIterator>
    inline void set_vertices(InputIterator, InputIterator) const { }

    inline void hide_point(Face_handle, const Point &) { }
    inline void hide(Point&, Face_handle) const { CGAL_triangulation_assertion(false); }
    inline void do_hide(const Point &, Face_handle) const { CGAL_triangulation_assertion(false); }

    template <class Conflict_tester>
    inline void hide_points(Vertex_handle, const Conflict_tester &) { }

    template <class Tester>
    inline bool replace_vertex(const Point&, Vertex_handle, const Tester&) const { return true; }
    inline Vertex_handle replace_vertex(Face_handle f, int index, const Point &) { return f->vertex(index); }
    inline void reinsert_vertices(Vertex_handle) { }
  };

  class Conflict_tester
  {
    // stores a pointer to the triangulation, a point, and an offset
    const Self* tr_ptr;
    Point p;
    mutable Offset o;

  public:
    /// Constructor
    Conflict_tester(const Self* tr_ptr) : tr_ptr(tr_ptr), p(Point()) { }
    Conflict_tester(const Point &pt, const Self* tr_ptr) : tr_ptr(tr_ptr), p(pt) { }

    /// returns true if the circumcircle of 'f' contains 'p'
    bool operator()(const Face_handle f, const Offset &off) const
    {
      return (tr_ptr->side_of_circle(f, p, tr_ptr->combine_offsets(o, off), true) == ON_BOUNDED_SIDE);
    }

    bool operator()(const Face_handle f, const Point& pt, const Offset &off) const
    {
      return (tr_ptr->side_of_circle(f, pt, o + off, true) == ON_BOUNDED_SIDE);
    }

    int compare_weight(Point, Point) const { return 0; }

    bool test_initial_face(Face_handle f, const Offset &off) const
    {
      if(!(operator()(f, off)))
        CGAL_triangulation_assertion(false);
      return true;
    }

    void set_point(const Point &_p) { p = _p; }
    const Point& point() const { return p; }
    void set_offset(const Offset &off) const { o = off; }
    const Offset& get_offset() const { return o; }
  };

  bool test_conflict(const Point& p, Face_handle fh) const
  {
    return side_of_oriented_circle(fh, p, true) ==  ON_POSITIVE_SIDE;
  }

  /// \name Check - Query

  Vertex_handle nearest_vertex_2D(const Point& p, Face_handle f) const;

  void look_nearest_neighbor(const Point& p,
                             Face_handle f,
                             int i,
                             Vertex_handle& nn) const;

  /// Returns the vertex closest to p, the point location will start from f
  Vertex_handle nearest_vertex(const Point& p, Face_handle f = Face_handle()) const;

public:
  /// \name Dual

  /// Returns the dual of f, which is the circumcenter of f.
  Point dual(Face_handle f) const;
  /// Returns the dual of e, which is always a segment in the periodic triangulation.
  Segment dual(const Edge& e) const ;
  /// Returns the dual of the edge pointed to by ec.
  Segment dual(const Edge_circulator& ec) const;
  /// Returns the dual of the edge pointed to by ei.
  Segment dual(const Edge_iterator& ei) const;

  template <class Stream>
  Stream& draw_dual(Stream& ps)
  {
    Edge_iterator eit = edges_begin(), eend = edges_end();
    for(; eit!=eend; ++eit)
      ps << dual(eit);

    return ps;
  }

public:
  /// \name Checking
  bool is_valid(bool verbose = false, int level = 0) const
  {
    // Check the parent
    bool result = Periodic_2_triangulation_2<Gt, Tds>::is_valid(verbose, level);

    // Check in_sphere:
    if(dimension() == 2)
    {
      const Point *p[4];
      Offset off[4];
      for(Face_iterator fit = faces_begin(); fit != this->faces_end(); ++fit)
      {
        for(int i=0; i<3; ++i)
        {
          p[i] = &fit->vertex(i)->point();
          off[i] = get_offset(fit, i);
        }

        /// Check whether the vertices of the neighbor lie outside the circumcircle of the face
        for(int i=0; i<3; ++i)
        {
          p[3] = &fit->vertex(i)->point();
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

private:
  /// \name Methods regarding the covering
  /// \{

  class Cover_manager
  {
    Self& tr;

  public:
    Cover_manager(Self& tr) : tr(tr) { }

    void create_initial_triangulation() { tr.create_initial_triangulation(); }

    template <class FaceIt>
    void insert_unsatisfying_elements(Vertex_handle v, const FaceIt begin, const FaceIt end) {
      tr.insert_faces_with_too_big_orthoball(v, begin, end);
    }

    template <class FaceIt>
    void delete_unsatisfying_elements(const FaceIt begin, const FaceIt end) {
      tr.delete_faces_with_too_big_circumdisk(begin, end);
    }

    bool can_be_converted_to_1_sheet() const { return tr.can_be_converted_to_1_sheet(); }

    bool update_cover_data_during_management(Face_handle new_fh,
                                             const std::vector<Face_handle>& new_faces)
    {
      return tr.update_cover_data_during_management(new_fh, new_faces);
    }
  };

public:
  void create_initial_triangulation()
  {
    CGAL_triangulation_assertion(faces_with_too_big_circumdisk.empty());

    for(Face_iterator iter=faces_begin(), end_iter = faces_end(); iter!=end_iter; ++iter)
      faces_with_too_big_circumdisk.insert(iter);
  }

  template <class FaceIt>
  void insert_faces_with_too_big_orthoball(Vertex_handle /*v*/, FaceIt begin, const FaceIt end)
  {
    for(; begin != end; ++begin)
      if(compare_squared_circumradius_to_threshold(*begin, squared_circumradius_threshold) != CGAL::SMALLER)
        faces_with_too_big_circumdisk.insert(*begin);
  }

  void insert_faces_with_too_big_orthoball(Face_iterator begin, Face_iterator end)
  {
    for(; begin != end; ++begin)
      if(compare_squared_circumradius_to_threshold(begin, squared_circumradius_threshold) != CGAL::SMALLER)
        faces_with_too_big_circumdisk.insert(begin);
  }

  template <class FaceIt>
  void delete_faces_with_too_big_circumdisk(FaceIt begin, const FaceIt end)
  {
    for(; begin != end; ++begin)
    {
      Too_big_circumdisks_set_it iter = faces_with_too_big_circumdisk.find(*begin);
      if(iter != faces_with_too_big_circumdisk.end())
        faces_with_too_big_circumdisk.erase(iter);
    }
  }

  bool can_be_converted_to_1_sheet() const { return faces_with_too_big_circumdisk.empty(); }

  // returns 'true/false' depending on whether the cover would (or has, if 'abort_if_cover_change'
  // is set to 'false') change.
  bool update_cover_data_during_management(Face_handle new_fh,
                                           const std::vector<Face_handle>& new_faces)
  {
    if(compare_squared_circumradius_to_threshold(new_fh, squared_circumradius_threshold) != CGAL::SMALLER)
    {
      if(is_1_cover())
      {
        // Whether we are changing the cover or simply aborting, we need to get rid of the new faces
        tds().delete_faces(new_faces.begin(), new_faces.end());
        return true;
      }
      else
      {
        faces_with_too_big_circumdisk.insert(new_fh);
      }
    }

    return false;
  }

  virtual void update_cover_data_after_converting_to_9_sheeted_covering()
  {
    for(Face_iterator iter = faces_begin(), end_iter = faces_end(); iter != end_iter; ++iter)
      if(compare_squared_circumradius_to_threshold(iter, squared_circumradius_threshold) != CGAL::SMALLER)
        faces_with_too_big_circumdisk.insert(iter);
  }

  virtual void clear_covering_data()
  {
    faces_with_too_big_circumdisk.clear();
  }

public:
  /// Checks whether the triangulation is a valid simplicial complex in the one cover.
  /// Uses an edge-length-criterion.
  bool is_extensible_triangulation_in_1_sheet_h1() const
  {
    if(!is_1_cover())
      return can_be_converted_to_1_sheet();

    return is_extensible_triangulation_in_1_sheet_h2();
  }

  /// Checks whether the triangulation is a valid simplicial complex in the one cover.
  /// Uses a criterion based on the maximal radius of the circumscribing circle.
  bool is_extensible_triangulation_in_1_sheet_h2() const
  {
    for(Periodic_triangle_iterator tit = this->periodic_triangles_begin(Base::UNIQUE);
        tit != this->periodic_triangles_end(Base::UNIQUE); ++tit)
    {
      if(compare_squared_circumradius_to_threshold(tit->at(0), tit->at(1), tit->at(2),
                                                   squared_circumradius_threshold) != CGAL::SMALLER)
        return false;
    }

    return true;
  }

private:
  /// This threshold should be chosen such that if all Delaunay balls have a squared radius smaller than this,
  /// we can be sure that there are no self-edges anymore.
  FT squared_circumradius_threshold;

  /// This container stores all the faces whose circumdisk squared radius is larger
  /// than the treshold `squared_circumradius_threshold`.
  Too_big_circumdisks_set faces_with_too_big_circumdisk;
};

template <class Gt, class Tds>
typename Periodic_2_Delaunay_triangulation_2<Gt, Tds>::Vertex_handle
Periodic_2_Delaunay_triangulation_2<Gt, Tds>::
nearest_vertex(const Point& p, Face_handle f) const
{
  switch(dimension())
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

template <class Gt, class Tds>
typename Periodic_2_Delaunay_triangulation_2<Gt, Tds>::Vertex_handle
Periodic_2_Delaunay_triangulation_2<Gt, Tds>::
nearest_vertex_2D(const Point& p, Face_handle f) const
{
  CGAL_triangulation_precondition(dimension() == 2);
  f = locate(p, f);

  typename Geom_traits::Compare_distance_2 compare_distance = geom_traits().compare_distance_2_object();

  Vertex_handle nn =  f->vertex(0);
  if(compare_distance(p, f->vertex(1)->point(), nn->point()) == SMALLER)
    nn = f->vertex(1);
  if(compare_distance(p, f->vertex(2)->point(), nn->point()) == SMALLER)
    nn = f->vertex(2);

  look_nearest_neighbor(p, f, 0, nn);
  look_nearest_neighbor(p, f, 1, nn);
  look_nearest_neighbor(p, f, 2, nn);

  return nn;
}

template <class Gt, class Tds>
void
Periodic_2_Delaunay_triangulation_2<Gt, Tds>::
look_nearest_neighbor(const Point& p,
                      Face_handle f,
                      int i,
                      Vertex_handle& nn) const
{
  Face_handle  ni = f->neighbor(i);
  if(this->side_of_oriented_circle(ni, p, true) != ON_POSITIVE_SIDE)
    return;

  typename Geom_traits::Compare_distance_2 compare_distance = geom_traits().compare_distance_2_object();
  i = ni->index(f);
  if(compare_distance(p, ni->vertex(i)->point(), nn->point()) == SMALLER)
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
dual(Face_handle f) const
{
  CGAL_triangulation_precondition(dimension() == 2);
  return circumcenter(f);
}

template <class Gt, class Tds>
inline typename Gt::Segment_2
Periodic_2_Delaunay_triangulation_2<Gt, Tds>::
dual(const Edge& e) const
{
  // dimension==2
  Face_handle nb = e.first->neighbor(e.second);
  Point p0 = dual(e.first);
  Point p1 = dual(nb);
  Offset o = combine_offsets(Offset(), get_neighbor_offset(e.first, e.second));
  Segment s = construct_segment(p0, p1, o, Offset());

  return s;
}

template <class Gt, class Tds>
inline typename Gt::Segment_2
Periodic_2_Delaunay_triangulation_2<Gt, Tds>::
dual(const Edge_circulator& ec) const
{
  return dual(*ec);
}

template <class Gt, class Tds>
inline typename Gt::Segment_2
Periodic_2_Delaunay_triangulation_2<Gt, Tds>::
dual(const Edge_iterator& ei) const
{
  return dual(*ei);
}

///////////////////////////////////////////////////////////////
//  INSERT

template <class Gt, class Tds>
bool
Periodic_2_Delaunay_triangulation_2<Gt, Tds>::
locally_Delaunay(Face_handle f, int i, Face_handle nb)
{
  CGAL_BRANCH_PROFILER("locally_Delaunay(), simplicity check failures", tmp);

  bool simplicity_criterion = is_1_cover() && f->has_zero_offsets() && nb->has_zero_offsets();

  const Point *p[4];
  for(int index = 0; index < 3; ++index)
    p[index]   = &nb->vertex(index)->point();

  p[3] = &f->vertex(i)->point();

  Oriented_side os;
  if(simplicity_criterion)
  {
    // No periodic offsets
    os = side_of_oriented_circle(*p[0], *p[1], *p[2], *p[3], true);
  }
  else
  {
    CGAL_BRANCH_PROFILER_BRANCH(tmp);

    Offset off[4];

    for(int index=0; index<3; ++index)
      off[index] = get_offset(nb, index);

    off[3] = combine_offsets(get_offset(f, i), get_neighbor_offset(f, i));

    os = side_of_oriented_circle(*p[0], *p[1], *p[2], *p[3],
                                 off[0], off[1], off[2], off[3], true);
  }

  return (ON_POSITIVE_SIDE != os);
}

///////////////////////////////////////////////////////////////
//  REMOVE    see INRIA RResearch Report 7104

template <class Gt, class Tds>
void
Periodic_2_Delaunay_triangulation_2<Gt, Tds>::
remove(Vertex_handle v)
{
  // @todo
  CGAL_assertion(false);
}

namespace P2DT2 {
namespace internal {

template<class P2DT2>
struct Static_data
{
  Static_data(int m) : maxd(m), f(maxd), i(maxd), w(maxd), offset_w(maxd) { }

  int maxd;
  std::vector<typename P2DT2::Face_handle> f;
  std::vector<int> i;
  std::vector<typename P2DT2::Vertex_handle> w;
  std::vector<typename P2DT2::Offset> offset_w;
};

} // namespace internal
} // namespace P2DT2

template<class Gt, class Tds>
Oriented_side
Periodic_2_Delaunay_triangulation_2<Gt, Tds>::
side_of_oriented_circle(const Point& p0, const Point& p1, const Point& p2, const Point& p, bool perturb) const
{
  Oriented_side os = geom_traits().side_of_oriented_circle_2_object()(p0, p1, p2, p);
  if((os != ON_ORIENTED_BOUNDARY) || (!perturb))
    return os;

  // We are now in a degenerate case => we do a symbolic perturbation.

  // We sort the points lexicographically.
  const Point * points[4] = { &p0, &p1, &p2, &p };
  std::sort(points, points + 4, typename Base::Perturbation_order(this));

  // We successively look whether the leading monomial, then 2nd monomial
  // of the determinant has non null coefficient.
  // 2 iterations are enough (cf paper)
  for(int i = 3; i > 0; --i)
  {
    if(points[i] == &p)
      return ON_NEGATIVE_SIDE; // since p0 p1 p2 are non collinear and positively oriented

    Orientation o;
    if(points[i] == &p2 && (o = orientation(p0, p1, p)) != COLLINEAR)
      return Oriented_side(o);
    if(points[i] == &p1 && (o = orientation(p0, p, p2)) != COLLINEAR)
      return Oriented_side(o);
    if(points[i] == &p0 && (o = orientation(p, p1, p2)) != COLLINEAR)
      return Oriented_side(o);
  }

  CGAL_triangulation_assertion(false);
  return ON_NEGATIVE_SIDE;
}

template<class Gt, class Tds>
Oriented_side
Periodic_2_Delaunay_triangulation_2<Gt, Tds>::
side_of_oriented_circle(const Point& p0, const Point& p1, const Point& p2, const Point& p,
                        const Offset& o0, const Offset& o1, const Offset& o2, const Offset& o,
                        bool perturb) const
{
  Oriented_side os = geom_traits().side_of_oriented_circle_2_object()(p0, p1, p2, p, o0, o1, o2, o);
  if((os != ON_ORIENTED_BOUNDARY) || (!perturb))
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
  for(int i = 3; i > 0; --i)
  {
    if(points[i] == &pts[3])
      return ON_NEGATIVE_SIDE; // since p0 p1 p2 are non collinear and positively oriented

    Orientation orient;
    if((points[i] == &pts[2]) && ((orient = orientation(p0, p1, p, o0, o1, o)) != COLLINEAR))
      return Oriented_side(orient);
    if((points[i] == &pts[1]) && ((orient = orientation(p0, p, p2, o0, o, o2)) != COLLINEAR))
      return Oriented_side(orient);
    if((points[i] == &pts[0]) && ((orient = orientation(p, p1, p2, o, o1, o2)) != COLLINEAR))
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
  while(os == ON_NEGATIVE_SIDE && i < 4)
  {
    os = side_of_oriented_circle(f->vertex(0)->point(), f->vertex(1)->point(), f->vertex(2)->point(), p,
                                 get_offset(f, 0), get_offset(f, 1), get_offset(f, 2), combine_offsets(Offset(), int_to_off(i)),
                                 perturb);
    ++i;
  }

  return os;
}

///////////////////////////////////////////////////////////////
//  DISPLACEMENT

template <class Gt, class Tds >
typename Periodic_2_Delaunay_triangulation_2<Gt, Tds>::Vertex_handle
Periodic_2_Delaunay_triangulation_2<Gt, Tds>::
move_if_no_collision(Vertex_handle v, const Point& p)
{
  Locate_type lt;
  int li;
  Vertex_handle inserted;
  Face_handle loc = locate(p, lt, li, v->face());

  if(lt == Base::VERTEX)
    return v;
  else
    /// This can be optimized by checking whether we can move v->point() to p
    return insert(p, lt, loc, li);
}

template <class Gt, class Tds >
typename Periodic_2_Delaunay_triangulation_2<Gt, Tds>::Vertex_handle
Periodic_2_Delaunay_triangulation_2<Gt, Tds>::
move_point(Vertex_handle v, const Point& p)
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

} // namespace CGAL

#endif // CGAL_PERIODIC_2_DELAUNAY_TRIANGULATION_2_H
