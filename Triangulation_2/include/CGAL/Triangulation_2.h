// Copyright (c) 1997-2010  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Olivier Devillers, Mariette Yvinec

#ifndef CGAL_TRIANGULATION_2_H
#define CGAL_TRIANGULATION_2_H

#include <CGAL/license/Triangulation_2.h>

#include <list>
#include <vector>
#include <map>
#include <algorithm>
#include <utility>
#include <iostream>

#include <CGAL/iterator.h>
#include <CGAL/function_objects.h>

#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_utils_2.h>

#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/Triangulation_vertex_base_2.h>
#include <CGAL/Triangulation_face_base_2.h>
#include <CGAL/Triangulation_line_face_circulator_2.h>
#include <CGAL/spatial_sort.h>
#include <CGAL/Spatial_sort_traits_adapter_2.h>

#include <CGAL/double.h>
#include <CGAL/internal/boost/function_property_map.hpp>

#include <boost/utility/result_of.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_smallint.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/iterator/transform_iterator.hpp>

#ifndef CGAL_NO_STRUCTURAL_FILTERING
#include <CGAL/internal/Static_filters/tools.h>
#include <CGAL/Triangulation_structural_filtering_traits.h>
#include <CGAL/determinant.h>
#endif // no CGAL_NO_STRUCTURAL_FILTERING

namespace CGAL {
template < class Gt, class Tds > class Triangulation_2;
template < class Gt, class Tds > std::istream& operator>>
  (std::istream& is, Triangulation_2<Gt,Tds> &tr);
template < class Gt, class Tds >  std::ostream& operator<<
  (std::ostream& os, const Triangulation_2<Gt,Tds> &tr);

#ifndef CGAL_NO_STRUCTURAL_FILTERING
namespace internal {
// structural filtering is performed only for EPIC
struct Structural_filtering_2_tag {};
struct No_structural_filtering_2_tag {};

template <bool filter>
struct Structural_filtering_selector_2 {
  typedef No_structural_filtering_2_tag  Tag;
};

template <>
struct Structural_filtering_selector_2<true> {
  typedef Structural_filtering_2_tag  Tag;
};
}
#endif // no CGAL_NO_STRUCTURAL_FILTERING

template < class Gt,
           class Tds = Triangulation_data_structure_2 <
                             Triangulation_vertex_base_2<Gt>,
                             Triangulation_face_base_2<Gt> > >
class Triangulation_2
  : public Triangulation_cw_ccw_2
{
  friend std::istream& operator>> <>
                (std::istream& is, Triangulation_2 &tr);
  typedef Triangulation_2<Gt,Tds>             Self;

public:
  typedef Tds                                 Triangulation_data_structure;
  typedef Gt                                  Geom_traits;

  // point types
  typedef typename Gt::Point_2                Point_2;
  typedef typename Tds::Vertex::Point         Point;

  typedef typename Geom_traits::Segment_2     Segment;
  typedef typename Geom_traits::Triangle_2    Triangle;
  typedef typename Geom_traits::Orientation_2 Orientation_2;
  typedef typename Geom_traits::Compare_x_2   Compare_x;
  typedef typename Geom_traits::Compare_y_2   Compare_y;

  typedef typename Tds::size_type              size_type;
  typedef typename Tds::difference_type        difference_type;

  typedef typename Tds::Vertex                 Vertex;
  typedef typename Tds::Face                   Face;
  typedef typename Tds::Edge                   Edge;
  typedef typename Tds::Vertex_handle          Vertex_handle;
  typedef typename Tds::Face_handle            Face_handle;

  typedef typename Tds::Face_circulator        Face_circulator;
  typedef typename Tds::Vertex_circulator      Vertex_circulator;
  typedef typename Tds::Edge_circulator        Edge_circulator;

  typedef typename Tds::Face_iterator          All_faces_iterator;
  typedef typename Tds::Edge_iterator          All_edges_iterator;
  typedef typename Tds::Halfedge_iterator      All_halfedges_iterator;
  typedef typename Tds::Vertex_iterator        All_vertices_iterator;

  class Perturbation_order
  {
    const Self *t;

  public:
    Perturbation_order(const Self *tr) : t(tr) {}

    bool operator()(const Point *p, const Point *q) const {
      return t->compare_xy(*p, *q) == SMALLER;
    }
  };

  friend class Perturbation_order;

  // This class is used to generate the Finite_*_iterators.
  class Infinite_tester
  {
    const Triangulation_2 *t;
  public:
    Infinite_tester() {}
    Infinite_tester(const Triangulation_2 *tr)	  : t(tr) {}

    bool operator()(const All_vertices_iterator & vit) const  {
      return t->is_infinite(vit);
    }
    bool operator()(const All_faces_iterator & fit ) const {
      return t->is_infinite(fit);
    }
    bool operator()(const All_edges_iterator & eit) const {
      return t->is_infinite(eit);
    }
  };

  //We derive in order to add a conversion to handle.
  class Finite_vertices_iterator :
    public Filter_iterator<All_vertices_iterator, Infinite_tester>
  {
    typedef Filter_iterator<All_vertices_iterator, Infinite_tester> Base;
    typedef Finite_vertices_iterator                          Self;
  public:
    Finite_vertices_iterator() : Base() {}
    Finite_vertices_iterator(const Base &b) : Base(b) {}
    Self & operator++() { Base::operator++(); return *this; }
    Self & operator--() { Base::operator--(); return *this; }
    Self operator++(int) { Self tmp(*this); ++(*this); return tmp; }
    Self operator--(int) { Self tmp(*this); --(*this); return tmp; }
    operator Vertex_handle() const { return Base::base(); }
  };

  class Finite_faces_iterator
    : public Filter_iterator<All_faces_iterator, Infinite_tester>
  {
    typedef Filter_iterator<All_faces_iterator, Infinite_tester> Base;
    typedef Finite_faces_iterator                           Self;
  public:
    Finite_faces_iterator() : Base() {}
    Finite_faces_iterator(const Base &b) : Base(b) {}
    Self & operator++() { Base::operator++(); return *this; }
    Self & operator--() { Base::operator--(); return *this; }
    Self operator++(int) { Self tmp(*this); ++(*this); return tmp; }
    Self operator--(int) { Self tmp(*this); --(*this); return tmp; }
    operator Face_handle() const { return Base::base(); }
  };

  typedef Filter_iterator<All_edges_iterator,
                           Infinite_tester>    Finite_edges_iterator;

  //for backward compatibility
  typedef Finite_faces_iterator                Face_iterator;
  typedef Finite_edges_iterator                Edge_iterator;
  typedef Finite_vertices_iterator             Vertex_iterator;

  typedef Triangulation_line_face_circulator_2<Self>  Line_face_circulator;

  // Auxiliary iterators for convenience
  // do not use default template argument to please VC++
  typedef Project_point<Vertex>                           Proj_point;

  typedef boost::transform_iterator<Proj_point,Finite_vertices_iterator> Point_iterator;
  typedef Point                value_type; // to have a back_inserter
  typedef const value_type&    const_reference;
  typedef value_type&          reference;

  enum Locate_type {VERTEX=0,
        EDGE, //1
        FACE, //2
        OUTSIDE_CONVEX_HULL, //3
        OUTSIDE_AFFINE_HULL}; //4

  //Tag to distinguish regular triangulations from others;
  typedef Tag_false  Weighted_tag;

protected:
  Gt _gt;
  Tds _tds;
  Vertex_handle _infinite_vertex;

public:
  // CONSTRUCTORS
  Triangulation_2(const Geom_traits& geom_traits=Geom_traits());
  Triangulation_2(const Triangulation_2<Gt,Tds> &tr);

  //Assignement
  Triangulation_2 &operator=(const Triangulation_2 &tr);

  //Helping
  void copy_triangulation(const Triangulation_2 &tr);
  void swap(Triangulation_2 &tr);
  void clear();

  //ACCESS FUNCTION
  const Geom_traits& geom_traits() const { return _gt;}
  const Tds & tds() const { return _tds;}
  Tds & tds() { return _tds;}
  int dimension() const { return _tds.dimension();}
  size_type number_of_vertices() const {return _tds.number_of_vertices() - 1;}
  size_type number_of_faces() const;
  Vertex_handle infinite_vertex() const;
  Vertex_handle finite_vertex() const;
  Face_handle infinite_face() const;
  Infinite_tester infinite_tester() const;

  //SETTING
  void set_infinite_vertex(const Vertex_handle& v) {_infinite_vertex=v;}

  // CHECKING
  bool is_valid(bool verbose = false, int level = 0) const;

  // TEST INFINITE FEATURES AND OTHER FEATURES
  bool is_infinite(Face_handle f) const;
  bool is_infinite(Vertex_handle v) const;
  bool is_infinite(Face_handle f, int i) const;
  bool is_infinite(const Edge& e) const;
  bool is_infinite(const Edge_circulator& ec) const;
  bool is_infinite(const All_edges_iterator& ei) const;
  bool is_edge(Vertex_handle va, Vertex_handle vb) const;
  bool is_edge(Vertex_handle va, Vertex_handle vb, Face_handle& fr,
               int & i) const;
  bool includes_edge(Vertex_handle va, Vertex_handle vb,
                     Vertex_handle& vbb,
                     Face_handle& fr, int & i) const;
  bool is_face(Vertex_handle v1, Vertex_handle v2, Vertex_handle v3) const;
  bool is_face(Vertex_handle v1, Vertex_handle v2, Vertex_handle v3,
               Face_handle &fr) const;

  // GEOMETRIC FEATURES AND CONSTRUCTION
  Point_2 construct_point(const Point& p) const;
  Triangle triangle(Face_handle f) const;
  Segment segment(Face_handle f, int i) const;
  Segment segment(const Edge& e) const;
  Segment segment(const Edge_circulator& ec) const;
  Segment segment(const All_edges_iterator& ei) const;
  Segment segment(const Finite_edges_iterator& ei) const;
  Point_2 circumcenter(Face_handle f) const;
  Point_2 circumcenter(const Point& p0,
                       const Point& p1,
                       const Point& p2) const;


  //MOVE - INSERTION - DELETION - Flip
public:
  void flip(Face_handle f, int i);

  Vertex_handle insert_first(const Point& p);
  Vertex_handle insert_second(const Point& p);
  Vertex_handle insert_in_edge(const Point& p, Face_handle f,int i);
  Vertex_handle insert_in_face(const Point& p, Face_handle f);
  Vertex_handle insert_outside_convex_hull(const Point& p, Face_handle f);
  Vertex_handle insert_outside_affine_hull(const Point& p);
  Vertex_handle insert(const Point &p, Face_handle start = Face_handle() );
  Vertex_handle insert(const Point& p,
                       Locate_type lt,
                       Face_handle loc, int li );
//  template < class InputIterator >
//  std::ptrdiff_t insert(InputIterator first, InputIterator last);
  Vertex_handle push_back(const Point& a);

  void remove_degree_3(Vertex_handle v, Face_handle f = Face_handle());
  void remove_first(Vertex_handle v);
  void remove_second(Vertex_handle v);
  void remove(Vertex_handle v);

  // MOVE
  Vertex_handle move_if_no_collision(Vertex_handle v, const Point &p);
  Vertex_handle move(Vertex_handle v, const Point &p);

protected: // some internal methods

  // INSERT, REMOVE, MOVE GIVING NEW FACES
  template <class OutputItFaces>
  Vertex_handle insert_and_give_new_faces(const Point &p,
                                          OutputItFaces fit,
                                          Face_handle start = Face_handle() );
  template <class OutputItFaces>
  Vertex_handle insert_and_give_new_faces(const Point& p,
                                          Locate_type lt,
                                          Face_handle loc, int li,
                                          OutputItFaces fit);

  template <class OutputItFaces>
  void remove_and_give_new_faces(Vertex_handle v,
                                 OutputItFaces fit);

  template <class OutputItFaces>
  Vertex_handle move_if_no_collision_and_give_new_faces(Vertex_handle v,
                                                        const Point &p,
                                                        OutputItFaces fit);

public:
  // POINT LOCATION
  Face_handle
  march_locate_1D(const Point& t, Locate_type& lt, int& li) const ;
  Face_handle
  march_locate_2D(Face_handle start,
                  const Point& t,
                  Locate_type& lt,
                  int& li) const;

  Face_handle
  march_locate_2D_LFC(Face_handle start,
                      const Point& t,
                      Locate_type& lt,
                      int& li) const;

  void
  compare_walks(const Point& p,
                Face_handle c1, Face_handle c2,
                Locate_type& lt1, Locate_type& lt2,
                int li1, int li2) const;

#ifdef CGAL_NO_STRUCTURAL_FILTERING
  Face_handle
  locate(const Point& p,
         Locate_type& lt,
         int& li,
         Face_handle start = Face_handle()) const;

  Face_handle
  locate(const Point &p, Face_handle start = Face_handle()) const;
#else  // no CGAL_NO_STRUCTURAL_FILTERING
#  ifndef CGAL_T2_STRUCTURAL_FILTERING_MAX_VISITED_CELLS
#    define CGAL_T2_STRUCTURAL_FILTERING_MAX_VISITED_CELLS 2500
#  endif // no CGAL_T2_STRUCTURAL_FILTERING_MAX_VISITED_CELLS

protected:
  Face_handle
  exact_locate(const Point& p,
               Locate_type& lt,
               int& li,
               Face_handle start) const;

  Face_handle
  generic_locate(const Point& p,
                 Locate_type& lt,
                 int& li,
                 Face_handle start,
                 internal::Structural_filtering_2_tag) const {
    return exact_locate(p, lt, li, inexact_locate(p, start));
  }

  Face_handle
  generic_locate(const Point& p,
                 Locate_type& lt,
                 int& li,
                 Face_handle start,
                 internal::No_structural_filtering_2_tag) const {
    return exact_locate(p, lt, li, start);
  }

  bool has_inexact_negative_orientation(const Point &p, const Point &q,
                                        const Point &r) const;

public:
  Face_handle
  inexact_locate(const Point& p,
                 Face_handle start = Face_handle(),
                 int max_num_cells =
                 CGAL_T2_STRUCTURAL_FILTERING_MAX_VISITED_CELLS) const;

  Face_handle
  locate(const Point & p,
         Locate_type & lt, int & li,
         Face_handle start = Face_handle()) const
  {
    typedef Triangulation_structural_filtering_traits<Geom_traits> TSFT;
    typedef typename internal::Structural_filtering_selector_2<
      TSFT::Use_structural_filtering_tag::value >::Tag Should_filter_tag;

    return generic_locate(p, lt, li, start, Should_filter_tag());
  }

  Face_handle
  locate(const Point & p, Face_handle start = Face_handle()) const
  {
    Locate_type lt;
    int li;
    return locate(p, lt, li, start);
  }

#endif // no CGAL_NO_STRUCTURAL_FILTERING

  //TRAVERSING : ITERATORS AND CIRCULATORS
  Finite_faces_iterator finite_faces_begin() const;
  Finite_faces_iterator finite_faces_end() const;
  Finite_vertices_iterator finite_vertices_begin() const;
  Finite_vertices_iterator finite_vertices_end() const;
  Finite_edges_iterator finite_edges_begin() const;
  Finite_edges_iterator finite_edges_end() const;
  Point_iterator points_begin() const;
  Point_iterator points_end() const;

  All_faces_iterator all_faces_begin() const;
  All_faces_iterator all_faces_end() const;
  All_vertices_iterator all_vertices_begin() const;
  All_vertices_iterator all_vertices_end() const;
  All_edges_iterator all_edges_begin() const;
  All_edges_iterator all_edges_end() const;
  All_halfedges_iterator all_halfedges_begin() const;
  All_halfedges_iterator all_halfedges_end() const;

  //for compatibility with previous versions
  Face_iterator faces_begin() const {return finite_faces_begin();}
  Face_iterator faces_end() const {return finite_faces_end();}
  Edge_iterator edges_begin() const {return finite_edges_begin();}
  Edge_iterator edges_end() const {return finite_edges_end();}
  Vertex_iterator vertices_begin() const {return finite_vertices_begin();}
  Vertex_iterator vertices_end() const {return finite_vertices_end();}

  Face_circulator incident_faces( Vertex_handle v,
                                  Face_handle f = Face_handle()) const;
  Vertex_circulator incident_vertices(Vertex_handle v,
                                      Face_handle f = Face_handle()) const;
  Edge_circulator incident_edges(Vertex_handle v,
                                 Face_handle f = Face_handle()) const;

  size_type degree(Vertex_handle v) const;

  Vertex_handle mirror_vertex(Face_handle f, int i) const;
  int mirror_index(Face_handle v, int i) const;
  Edge mirror_edge(Edge e) const;

  Line_face_circulator line_walk(const Point& p,
                                 const Point& q,
                                 Face_handle f = Face_handle()) const;

 // TO DEBUG
  void show_all() const;
  void show_vertex(Vertex_handle vh) const;
  void show_face( Face_handle fh) const;

  // IO
// template < class Stream >
// Stream& draw_triangulation(Stream& os) const;

  //PREDICATES
  Oriented_side
  oriented_side(const Point &p0, const Point &p1,
                const Point &p2, const Point &p) const;

  Bounded_side
  bounded_side(const Point &p0, const Point &p1,
               const Point &p2, const Point &p) const;

  Oriented_side
  oriented_side(Face_handle f, const Point &p) const;

  Oriented_side
  side_of_oriented_circle(const Point &p0, const Point &p1, const Point &p2,
                          const Point &p, bool perturb) const;

  Oriented_side
  side_of_oriented_circle(Face_handle f, const Point & p, bool perturb = false) const;

  bool
  collinear_between(const Point& p, const Point& q, const Point& r)
  const;

  Comparison_result compare_x(const Point& p, const Point& q) const;
  Comparison_result compare_xy(const Point& p, const Point& q) const;
  Comparison_result compare_y(const Point& p, const Point& q) const;
  bool               xy_equal(const Point& p, const Point& q) const;

  Orientation orientation(const Point& p,
                          const Point& q,
                          const Point& r) const;

protected:
  void remove_1D(Vertex_handle v);
  void remove_2D(Vertex_handle v);
  bool test_dim_down(Vertex_handle v) const;
  void fill_hole(Vertex_handle v, std::list<Edge> & hole);
  void fill_hole_delaunay(std::list<Edge> & hole);

  // output faces
  template <class OutputItFaces>
  void fill_hole(Vertex_handle v, std::list<Edge> & hole, OutputItFaces fit);

  template <class OutputItFaces>
  void fill_hole_delaunay(std::list<Edge> & hole, OutputItFaces fit);

  void make_hole(Vertex_handle v, std::list<Edge> & hole,
                 std::set<Face_handle> &faces_set);

public:
  void make_hole(Vertex_handle v, std::list<Edge> & hole);
//  template<class EdgeIt>
//  Vertex_handle star_hole( Point p,
//                           EdgeIt edge_begin,
//                           EdgeIt edge_end);

//  template<class EdgeIt, class FaceIt>
//  Vertex_handle star_hole( Point p,
//                           EdgeIt edge_begin,
//                           EdgeIt edge_end,
//                           FaceIt face_begin,
//                           FaceIt face_end);

  Face_handle create_face(Face_handle f1d, int i1,
                          Face_handle f2, int i2,
                          Face_handle f3, int i3);
  Face_handle create_face(Face_handle f1, int i1,
                          Face_handle f2, int i2);
  Face_handle create_face(Face_handle f, int i, Vertex_handle v);
  Face_handle create_face(Vertex_handle v1, Vertex_handle v2,Vertex_handle v3);
  Face_handle create_face(Vertex_handle v1, Vertex_handle v2,Vertex_handle v3,
                          Face_handle f1, Face_handle f2, Face_handle f3);
  Face_handle create_face();
  Face_handle create_face(Face_handle); //calls copy constructor of Face
  void delete_face(Face_handle f);
  void delete_vertex(Vertex_handle v);

  Vertex_handle collapse_edge(Edge e)
  {
    return _tds.collapse_edge(e);
  }

  Vertex_handle file_input(std::istream& is);
  void file_output(std::ostream& os) const;

private:
  Vertex_handle insert_outside_convex_hull_1(const Point& p, Face_handle f);
  Vertex_handle insert_outside_convex_hull_2(const Point& p, Face_handle f);

  // template members
public:
  template < class Stream >
  Stream& draw_triangulation(Stream& os) const
  {
    Finite_edges_iterator it = finite_edges_begin();
    for( ;it != finite_edges_end() ; ++it) {
      os << segment(it);
    }
    return os;
  }

template < class InputIterator >
std::ptrdiff_t insert(InputIterator first, InputIterator last)
{
  size_type n = number_of_vertices();

  std::vector<Point> points (first, last);

  typedef typename Geom_traits::Construct_point_2 Construct_point_2;
  typedef typename boost::result_of<const Construct_point_2(const Point&)>::type Ret;
  typedef CGAL::internal::boost_::function_property_map<Construct_point_2, Point, Ret> fpmap;
  typedef CGAL::Spatial_sort_traits_adapter_2<Geom_traits, fpmap> Search_traits_2;

  spatial_sort(points.begin(), points.end(),
               Search_traits_2(
                 CGAL::internal::boost_::make_function_property_map<Point, Ret, Construct_point_2>(
                   geom_traits().construct_point_2_object()), geom_traits()));

  Face_handle f;
  for (typename std::vector<Point>::const_iterator p = points.begin(), end = points.end();
          p != end; ++p)
      f = insert (*p, f)->face();

  return number_of_vertices() - n;
}

bool well_oriented(Vertex_handle v) const
{
  Face_circulator fc = incident_faces(v), done(fc);
  do {
    if(!is_infinite(fc)) {
      Vertex_handle v0 = fc->vertex(0);
      Vertex_handle v1 = fc->vertex(1);
      Vertex_handle v2 = fc->vertex(2);
      if(orientation(v0->point(),v1->point(),v2->point())
        != COUNTERCLOCKWISE) return false;
    }
  } while(++fc != done);
  return true;
}

bool from_convex_hull(Vertex_handle v) {
  CGAL_triangulation_precondition(!is_infinite(v));
  Vertex_circulator vc = incident_vertices(v), done(vc);
  do { if(is_infinite(vc)) return true; } while(++vc != done);
  return false;
}

public:
template<class EdgeIt>
Vertex_handle star_hole( const Point& p,
                         EdgeIt edge_begin,
                         EdgeIt edge_end) {
  std::list<Face_handle> empty_list;
  return star_hole(p,
                   edge_begin,
                   edge_end,
                   empty_list.begin(),
                   empty_list.end());
}

template<class EdgeIt, class FaceIt>
Vertex_handle star_hole( const Point& p,
                         EdgeIt edge_begin,
                         EdgeIt edge_end,
                         FaceIt face_begin,
                         FaceIt face_end) {
  Vertex_handle v = _tds.star_hole( edge_begin, edge_end,
                                    face_begin, face_end);
  v->set_point(p);
  return v;
}
};

// CONSTRUCTORS
template <class Gt, class Tds >
Triangulation_2<Gt, Tds>::
Triangulation_2(const Geom_traits& geom_traits)
  : _gt(geom_traits), _tds()
{
  _infinite_vertex = _tds.insert_first();
}

// copy constructor duplicates vertices and faces
template <class Gt, class Tds >
Triangulation_2<Gt, Tds>::
Triangulation_2(const Triangulation_2 &tr)
  : _gt(tr._gt)
{
  _infinite_vertex = _tds.copy_tds(tr._tds, tr.infinite_vertex());
}

//Assignement
template <class Gt, class Tds >
Triangulation_2<Gt, Tds> &
Triangulation_2<Gt, Tds>::
operator=(const Triangulation_2 &tr)
{
  copy_triangulation(tr);
  return *this;
}

// Helping functions

template <class Gt, class Tds >
void
Triangulation_2<Gt, Tds>::
copy_triangulation(const Triangulation_2 &tr)
{
  _tds.clear();
  _gt = tr._gt;
  _infinite_vertex = _tds.copy_tds(tr._tds, tr._infinite_vertex);
}

template <class Gt, class Tds >
void
Triangulation_2<Gt, Tds>::
swap(Triangulation_2 &tr)
{
  Vertex_handle v= _infinite_vertex;
  _infinite_vertex = tr._infinite_vertex;
  tr._infinite_vertex = v;

  _tds.swap(tr._tds);

  Geom_traits t = geom_traits();
  _gt = tr.geom_traits();
  tr._gt = t;
}

template <class Gt, class Tds >
void
Triangulation_2<Gt, Tds>::
clear()
{
  _tds.clear(); //detruit tous les sommets et toutes les faces
  _infinite_vertex = _tds.insert_first();
}

template <class Gt, class Tds >
typename Triangulation_2<Gt, Tds>::size_type
Triangulation_2<Gt, Tds>::
number_of_faces() const
{
  size_type count = _tds.number_of_faces();
  Face_circulator fc = incident_faces(infinite_vertex()), done(fc);
  if ( ! fc.is_empty() ) {
    do {
      --count; ++fc;
    } while (fc != done);
  }
  return count;
}

template <class Gt, class Tds >
inline
typename Triangulation_2<Gt,Tds>::Vertex_handle
Triangulation_2<Gt,Tds>::
infinite_vertex() const
{
  return _infinite_vertex;
}

template <class Gt, class Tds >
inline
typename Triangulation_2<Gt,Tds>::Vertex_handle
Triangulation_2<Gt,Tds>::
finite_vertex() const
{
  CGAL_triangulation_precondition (number_of_vertices() >= 1);
  return (finite_vertices_begin());
}

template <class Gt, class Tds >
inline
typename Triangulation_2<Gt,Tds>::Face_handle
Triangulation_2<Gt,Tds>::
infinite_face() const
{
  return infinite_vertex()->face();
}

template <class Gt, class Tds >
inline
typename Triangulation_2<Gt,Tds>::Infinite_tester
Triangulation_2<Gt,Tds>::
infinite_tester() const
{
  return Infinite_tester(this);
}

template <class Gt, class Tds >
bool
Triangulation_2<Gt, Tds>::
is_valid(bool verbose, int level) const
{
  bool result = _tds.is_valid(verbose, level);
  if (dimension() <= 0 || (dimension()==1 && number_of_vertices() == 2 ) )
    return result;

  if (dimension() == 1) {
    Finite_vertices_iterator it1 = finite_vertices_begin(),
                             it2(it1), it3(it1);
    ++it2;
    ++it3; ++it3;
    while( it3 != finite_vertices_end()) {
     Orientation s = orientation(it1->point(),
         it2->point(),
         it3->point());
     result = result && s == COLLINEAR ;
     CGAL_triangulation_assertion(result);
     ++it1 ; ++it2; ++it3;
    }
  }
  else { //dimension() == 2
    for(Finite_faces_iterator it=finite_faces_begin();
        it!=finite_faces_end(); it++) {
      CGAL_triangulation_assertion( ! is_infinite(it));
      Orientation s = orientation(it->vertex(0)->point(),
                                  it->vertex(1)->point(),
                                  it->vertex(2)->point());
      CGAL_triangulation_assertion( s == LEFT_TURN );
      result = result && ( s == LEFT_TURN );
    }

    Vertex_circulator start = incident_vertices(infinite_vertex());
    Vertex_circulator pc(start);
    Vertex_circulator qc(start); ++qc;
    Vertex_circulator rc(start); ++rc; ++rc;
    do {
      Orientation s = orientation(pc->point(),
                                  qc->point(),
                                  rc->point());
      CGAL_triangulation_assertion( s != LEFT_TURN );
      result = result && ( s != LEFT_TURN );
      ++pc ; ++qc ; ++rc;
    } while(pc != start);

    // check number of faces. This cannot be done by the Tds
    // which does not know the number of components nor the genus
    result = result && (number_of_faces() == 2*(number_of_vertices()+1)
                        - 4
                        - degree(infinite_vertex()));
    CGAL_triangulation_assertion( result);
  }
  return result;
}

template <class Gt, class Tds >
inline bool
Triangulation_2<Gt, Tds>::
is_infinite(Face_handle f) const
{
  return f->has_vertex(infinite_vertex());
}

template <class Gt, class Tds >
inline bool
Triangulation_2<Gt, Tds>::
is_infinite(Vertex_handle v) const
{
  return v == infinite_vertex();
}

template <class Gt, class Tds >
inline bool
Triangulation_2<Gt, Tds>::
is_infinite(Face_handle f, int i) const
{
  return is_infinite(f->vertex(ccw(i))) ||
         is_infinite(f->vertex(cw(i)));
}

template <class Gt, class Tds >
inline bool
Triangulation_2<Gt, Tds>::
is_infinite(const Edge& e) const
{
  return is_infinite(e.first,e.second);
}

template <class Gt, class Tds >
inline bool
Triangulation_2<Gt, Tds>::
is_infinite(const Edge_circulator& ec) const
{
  return is_infinite(*ec);
}

template <class Gt, class Tds >
inline bool
Triangulation_2<Gt, Tds>::
is_infinite(const All_edges_iterator& ei) const
{
  return is_infinite(*ei);
}

template <class Gt, class Tds >
inline bool
Triangulation_2<Gt, Tds>::
is_edge(Vertex_handle va, Vertex_handle vb) const
{
  return _tds.is_edge( va, vb);
}

template <class Gt, class Tds >
inline bool
Triangulation_2<Gt, Tds>::
is_edge(Vertex_handle va, Vertex_handle vb, Face_handle& fr, int & i) const
{
  return _tds.is_edge(va, vb, fr, i);
}

template <class Gt, class Tds >
bool
Triangulation_2<Gt, Tds>::
includes_edge(Vertex_handle va, Vertex_handle vb,
              Vertex_handle & vbb,
              Face_handle& fr, int & i) const
// returns true if the line segment ab contains an edge e of t
// incident to a, false otherwise
// if true, vbb becomes the vertex of e distinct from a
// fr is the face incident to e and e=(fr,i)
// fr is on the right side of a->b
{
  Vertex_handle v;
  Orientation orient;
  int indv;
  Edge_circulator ec = incident_edges(va), done(ec);
  if (ec != 0) {
    do {
      //find the index of the other vertex of *ec
      indv = 3 - ((*ec).first)->index(va) - (*ec).second ;
      v = ((*ec).first)->vertex(indv);
      if (!is_infinite(v)) {
        if (v==vb) {
          vbb = vb;
          fr=(*ec).first;
          i= (*ec).second;
          return true;
        }
        else {
          orient = orientation(va->point(),
                               vb->point(),
                               v->point());
          if((orient==COLLINEAR) &&
             (collinear_between (va->point(),
                                 v->point(),
                                 vb->point()))) {
            vbb = v;
            fr=(*ec).first;
            i= (*ec).second;
            return true;
          }
        }
      }
    } while (++ec != done);
  }
  return false;
}

template <class Gt, class Tds >
inline bool
Triangulation_2<Gt, Tds>::
is_face(Vertex_handle v1, Vertex_handle v2, Vertex_handle v3) const
{
  return _tds.is_face(v1, v2, v3);
}

template <class Gt, class Tds >
inline bool
Triangulation_2<Gt, Tds>::
is_face(Vertex_handle v1,
        Vertex_handle v2,
        Vertex_handle v3,
        Face_handle &fr) const
{
  return _tds.is_face(v1, v2, v3, fr);
}

template <class Gt, class Tds >
typename Triangulation_2<Gt, Tds>::Point_2
Triangulation_2<Gt, Tds>::
construct_point(const Point& p) const
{
  return geom_traits().construct_point_2_object()(p);
}

template <class Gt, class Tds >
typename Triangulation_2<Gt, Tds>::Triangle
Triangulation_2<Gt, Tds>::
triangle(Face_handle f) const
{
  CGAL_triangulation_precondition( ! is_infinite(f) );
  typename Gt::Construct_triangle_2
      construct_triangle = geom_traits().construct_triangle_2_object();
  return construct_triangle(construct_point(f->vertex(0)->point()),
                            construct_point(f->vertex(1)->point()),
                            construct_point(f->vertex(2)->point()));
}

template <class Gt, class Tds >
typename Triangulation_2<Gt, Tds>::Segment
Triangulation_2<Gt, Tds>::
segment(Face_handle f, int i) const
{
  CGAL_triangulation_precondition( ! is_infinite(f,i));
  typename Gt::Construct_segment_2
      construct_segment = geom_traits().construct_segment_2_object();
  return construct_segment(construct_point(f->vertex(ccw(i))->point()),
                           construct_point(f->vertex(cw(i))->point()));
}

template <class Gt, class Tds >
typename Triangulation_2<Gt, Tds>::Segment
Triangulation_2<Gt, Tds>::
segment(const Edge& e) const
{
  CGAL_triangulation_precondition(! is_infinite(e));
  typename Gt::Construct_segment_2
      construct_segment = geom_traits().construct_segment_2_object();
  return construct_segment(construct_point(e.first->vertex(ccw(e.second))->point()),
                           construct_point(e.first->vertex( cw(e.second))->point()));
}

template <class Gt, class Tds >
typename Triangulation_2<Gt, Tds>::Segment
Triangulation_2<Gt, Tds>::
segment(const Edge_circulator& ec) const
{
  return segment(*ec);
}

template <class Gt, class Tds >
typename Triangulation_2<Gt, Tds>::Segment
Triangulation_2<Gt, Tds>::
segment(const Finite_edges_iterator& ei) const
{
  return segment(*ei);
}

template <class Gt, class Tds >
typename Triangulation_2<Gt, Tds>::Segment
Triangulation_2<Gt, Tds>::
segment(const All_edges_iterator& ei) const
{
  return segment(*ei);
}

template <class Gt, class Tds >
void
Triangulation_2<Gt, Tds>::
flip(Face_handle f, int i)
{
  CGAL_triangulation_precondition ( f != Face_handle() );
  CGAL_triangulation_precondition (i == 0 || i == 1 || i == 2);
  CGAL_triangulation_precondition( dimension()==2);

  CGAL_triangulation_precondition( !is_infinite(f) &&
                                   !is_infinite(f->neighbor(i)) );
  CGAL_triangulation_precondition(
        orientation(f->vertex(i)->point(),
                    f->vertex(cw(i))->point(),
                    mirror_vertex(f,i)->point()) == RIGHT_TURN &&
        orientation(f->vertex(i)->point(),
                    f->vertex(ccw(i))->point(),
                    mirror_vertex(f,i)->point()) == LEFT_TURN);
  _tds.flip(f, i);
  return;
}

template <class Gt, class Tds >
typename Triangulation_2<Gt,Tds>::Vertex_handle
Triangulation_2<Gt,Tds>::
insert_first(const Point& p)
{
  CGAL_triangulation_precondition(number_of_vertices() == 0);
  Vertex_handle v = _tds.insert_second();
  v->set_point(p);
  return v;
}

template <class Gt, class Tds >
typename Triangulation_2<Gt,Tds>::Vertex_handle
Triangulation_2<Gt,Tds>::
insert_second(const Point& p)
{
  CGAL_triangulation_precondition(number_of_vertices() == 1);
   Vertex_handle v = _tds.insert_dim_up(infinite_vertex(), true);
   v->set_point(p);
   return v;
}

template <class Gt, class Tds >
typename Triangulation_2<Gt,Tds>::Vertex_handle
Triangulation_2<Gt,Tds>::
insert_in_edge(const Point& p, Face_handle f,int i)
{
 CGAL_triangulation_exactness_precondition(
        orientation(f->vertex(cw(i))->point(), p,
        f->vertex(ccw(i))->point()) == COLLINEAR &&
        collinear_between(f->vertex(cw(i))->point(), p,
        f->vertex(ccw(i))->point() ) );
  Vertex_handle v = _tds.insert_in_edge(f,i);
  v->set_point(p);
  return v;
}

template <class Gt, class Tds >
typename Triangulation_2<Gt,Tds>::Vertex_handle
Triangulation_2<Gt,Tds>::
insert_in_face(const Point& p, Face_handle f)
{
  CGAL_triangulation_precondition(oriented_side(f,p) == ON_POSITIVE_SIDE);
  Vertex_handle v= _tds.insert_in_face(f);
  v->set_point(p);
  return v;
}

template <class Gt, class Tds >
typename Triangulation_2<Gt,Tds>::Vertex_handle
Triangulation_2<Gt,Tds>::
insert_outside_convex_hull(const Point& p, Face_handle f)
{
  CGAL_triangulation_precondition(is_infinite(f) && dimension() >= 1);
  Vertex_handle v;
  if (dimension() == 1)
    v=insert_outside_convex_hull_1(p, f);
  else
    v=insert_outside_convex_hull_2(p, f);

  v->set_point(p);
  return v;
}

template <class Gt, class Tds >
typename Triangulation_2<Gt,Tds>::Vertex_handle
Triangulation_2<Gt,Tds>::
insert_outside_convex_hull_1(const Point& p, Face_handle f)
{
  CGAL_triangulation_precondition( is_infinite(f) && dimension()==1);
  CGAL_triangulation_precondition(
        orientation(mirror_vertex(f, f->index(infinite_vertex()))->point(),
                    f->vertex(1- f->index(infinite_vertex()))->point(),
                    p) == COLLINEAR &&
        collinear_between(mirror_vertex(f,f->index(infinite_vertex()))->point(),
                          f->vertex(1- f->index(infinite_vertex()))->point(),
                          p) );
  Vertex_handle v=_tds.insert_in_edge(f, 2);
  v->set_point(p);
  return v;
}

template <class Gt, class Tds >
typename Triangulation_2<Gt,Tds>::Vertex_handle
Triangulation_2<Gt,Tds>::
insert_outside_convex_hull_2(const Point& p, Face_handle f)
{
  CGAL_triangulation_precondition(is_infinite(f));

  int li = f->index(infinite_vertex());

  CGAL_triangulation_precondition(
        orientation(p,
                    f->vertex(ccw(li))->point(),
                    f->vertex(cw(li))->point()) == LEFT_TURN);

  std::list<Face_handle> ccwlist;
  std::list<Face_handle> cwlist;

  Face_circulator fc = incident_faces(infinite_vertex(), f);
  bool done = false;
  while(! done) {
    fc--;
    li = fc->index(infinite_vertex());
    const Point& q = fc->vertex(ccw(li))->point();
    const Point& r = fc->vertex(cw(li))->point();
    if(orientation(p,q,r) == LEFT_TURN ) { ccwlist.push_back(fc); }
    else {done=true;}
  }

  fc = incident_faces(infinite_vertex(), f);
  done = false;
  while(! done){
    fc++;
    li = fc->index(infinite_vertex());
    const Point& q = fc->vertex(ccw(li))->point();
    const Point& r = fc->vertex(cw(li))->point();
    if(orientation(p,q,r) == LEFT_TURN ) { cwlist.push_back(fc);}
    else {done=true;}
  }

  Vertex_handle v = _tds.insert_in_face(f);
  v->set_point(p);

  Face_handle fh;
  while ( ! ccwlist.empty()) {
    fh = ccwlist.front();
    li = ccw(fh->index(infinite_vertex()));
    _tds.flip( fh, li);
    ccwlist.pop_front();
  }

  while ( ! cwlist.empty()) {
    fh = cwlist.front();
    li = cw(fh->index(infinite_vertex()));
    _tds.flip( fh, li);
    cwlist.pop_front();
  }

  //reset infinite_vertex()->face();
  fc = incident_faces(v);
  while( ! is_infinite(fc)) { fc++;}
  infinite_vertex()->set_face(fc);

  return v;
}

template <class Gt, class Tds >
typename Triangulation_2<Gt,Tds>::Vertex_handle
Triangulation_2<Gt,Tds>::
insert_outside_affine_hull(const Point& p)
{
  CGAL_triangulation_precondition(dimension() < 2);
  bool conform = false;
  if (dimension() == 1) {
    Face_handle f = (*finite_edges_begin()).first;
    Orientation orient = orientation( f->vertex(0)->point(),
                                      f->vertex(1)->point(),
                                      p);
    CGAL_triangulation_precondition(orient != COLLINEAR);
    conform = ( orient == COUNTERCLOCKWISE);
  }
  Vertex_handle v = _tds.insert_dim_up( infinite_vertex(), conform);
  v->set_point(p);
  return v;
}

template <class Gt, class Tds >
typename Triangulation_2<Gt,Tds>::Vertex_handle
Triangulation_2<Gt,Tds>::
insert(const Point &p, Face_handle start)
{
  Locate_type lt;
  int li;
  Face_handle loc = locate (p, lt, li, start);
  return insert(p, lt, loc, li);
}

template <class Gt, class Tds >
typename Triangulation_2<Gt,Tds>::Vertex_handle
Triangulation_2<Gt,Tds>::
insert(const Point& p, Locate_type lt, Face_handle loc, int li)
  // insert a point p, whose localisation is known (lt, f, i)
{
  if(number_of_vertices() == 0) {
    return(insert_first(p));
  }

  if(number_of_vertices() == 1) {
    if (lt == VERTEX) return finite_vertex();
    else return(insert_second(p));
  }

  switch(lt) {
  case FACE:
    return insert_in_face(p,loc);
  case EDGE:
    return insert_in_edge(p,loc,li);
  case OUTSIDE_CONVEX_HULL:
    return insert_outside_convex_hull(p,loc);
  case OUTSIDE_AFFINE_HULL:
   return insert_outside_affine_hull(p);
  case VERTEX:
    return loc->vertex(li);
  }
  CGAL_triangulation_assertion(false); // locate step failed
  return Vertex_handle();
}


template <class Gt, class Tds >
inline
typename Triangulation_2<Gt,Tds>::Vertex_handle
Triangulation_2<Gt,Tds>::
push_back(const Point &p)
{
  return insert(p);
}

template <class Gt, class Tds >
inline void
Triangulation_2<Gt,Tds>::
remove_degree_3(Vertex_handle v, Face_handle f)
{
  if (f == Face_handle()) f=v->face();
  _tds.remove_degree_3(v, f);
  return;
}

template <class Gt, class Tds >
inline void
Triangulation_2<Gt,Tds>::
remove_first(Vertex_handle v)
{
  _tds.remove_second(v);
  return;
}

template <class Gt, class Tds >
inline void
Triangulation_2<Gt,Tds>::
remove_second(Vertex_handle v)
{
  _tds.remove_dim_down(v);
  return;
}

template <class Gt, class Tds >
void
Triangulation_2<Gt,Tds>::
remove(Vertex_handle v)
{
  CGAL_triangulation_precondition( v != Vertex_handle());
  CGAL_triangulation_precondition( !is_infinite(v));

  if  (number_of_vertices() == 1)
    remove_first(v);
  else if (number_of_vertices() == 2)
    remove_second(v);
  else if ( dimension() == 1)
    remove_1D(v);
  else
    remove_2D(v);
  return;
}

template <class Gt, class Tds >
inline void
Triangulation_2<Gt, Tds>::
remove_1D(Vertex_handle v)
{
  _tds.remove_1D(v);
}

template <class Gt, class Tds >
bool
Triangulation_2<Gt,Tds>::
test_dim_down(Vertex_handle v) const
{
  //test the dimensionality of the resulting triangulation
  //upon removing of vertex v
  //it goes down to 1 iff
  // 1) any finite face is incident to v
  // 2) all vertices are collinear
  CGAL_triangulation_precondition(dimension() == 2);
  bool dim1 = true;
  Finite_faces_iterator fit = finite_faces_begin();
  while (dim1==true && fit != finite_faces_end()) {
    dim1 = dim1 && fit->has_vertex(v);
    fit++;
  }
  Face_circulator fic = incident_faces(v);
  while (is_infinite(fic)) {++fic;}
  Face_circulator done(fic);
  Face_handle start(fic); int iv = start->index(v);
  const Point& p = start->vertex(cw(iv))->point();
  const Point& q = start->vertex(ccw(iv))->point();
  while ( dim1 && ++fic != done) {
    iv = fic->index(v);
    if (fic->vertex(ccw(iv)) != infinite_vertex()) {
      dim1 = dim1 &&
             orientation(p, q, fic->vertex(ccw(iv))->point()) == COLLINEAR;
    }
  }
  return dim1;
}

template <class Gt, class Tds >
void
Triangulation_2<Gt,Tds>::
remove_2D(Vertex_handle v)
{
  if (test_dim_down(v)) { _tds.remove_dim_down(v); }
  else {
    std::list<Edge> hole;
    make_hole(v, hole);
    fill_hole(v, hole);
    delete_vertex(v);
  }
  return;
}

template < class Gt, class Tds >
template < class OutputItFaces >
inline
typename Triangulation_2<Gt,Tds>::Vertex_handle
Triangulation_2<Gt,Tds>::
insert_and_give_new_faces(const Point &p,
                          OutputItFaces oif,
                          Face_handle start)
{
  Vertex_handle v = insert(p, start);
  int dimension = this->dimension();
  if(dimension == 2)
  {
    Face_circulator fc = incident_faces(v), done(fc);
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
typename Triangulation_2<Gt,Tds>::Vertex_handle
Triangulation_2<Gt,Tds>::
insert_and_give_new_faces(const Point &p,
                          Locate_type lt,
                          Face_handle loc, int li,
                          OutputItFaces oif)
{
  Vertex_handle v = insert(p, lt, loc, li);
  int dimension = this->dimension();
  if(dimension == 2)
  {
    Face_circulator fc = incident_faces(v), done(fc);
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
template <class OutputItFaces>
void
Triangulation_2<Gt,Tds>::
remove_and_give_new_faces(Vertex_handle v, OutputItFaces fit)
{
  CGAL_triangulation_precondition( v != Vertex_handle());
  CGAL_triangulation_precondition( !is_infinite(v));

  if(number_of_vertices() == 1) remove_first(v);
  else if(number_of_vertices() == 2) remove_second(v);
  else if( dimension() == 1)
  {
    Point p = v->point();
    remove(v);
    *fit++ = locate(p);
  }
  else if (test_dim_down(v)) {
    _tds.remove_dim_down(v);
    for(All_faces_iterator afi = tds().face_iterator_base_begin();
        afi != tds().face_iterator_base_begin();
        afi++) *fit++ = afi;
  }
  else {
    std::list<Edge> hole;
    make_hole(v, hole);
    fill_hole(v, hole, fit);
    delete_vertex(v);
  }
  return;
}

template <class Gt, class Tds >
void
Triangulation_2<Gt, Tds>::
make_hole ( Vertex_handle v, std::list<Edge> & hole)
{
  std::vector<Face_handle> to_delete;
  to_delete.reserve(16);

  Face_handle f, fn;
  int i, in ;
  Vertex_handle vv;

  Face_circulator fc = incident_faces(v);
  Face_circulator done(fc);
  do {
    f = fc; fc++;
    i = f->index(v);
    fn = f->neighbor(i);
    in = fn->index(f);
    vv = f->vertex(cw(i));
    vv->set_face(fn);
    vv = f->vertex(ccw(i));
    vv->set_face(fn);
    fn->set_neighbor(in, Face_handle());
    hole.push_back(Edge(fn,in));
    to_delete.push_back(f);
  } while(fc != done);

  std::size_t size = to_delete.size();
  for(std::size_t i=0; i<size; i++) {
    f = to_delete[i];
    delete_face(f);
  }
}

template <class Gt, class Tds >
void
Triangulation_2<Gt,Tds>::
make_hole(Vertex_handle v, std::list<Edge> & hole,
          std::set<Face_handle> &faces_set)
{
  std::vector<Face_handle> to_delete;
  to_delete.reserve(16);

  Face_handle f, fn;
  int i, in ;
  Vertex_handle vv;

  Face_circulator fc = incident_faces(v);
  Face_circulator done(fc);
  do {
    f = fc; fc++;
    i = f->index(v);
    fn = f->neighbor(i);
    in = fn->index(f);
    vv = f->vertex(cw(i));
    vv->set_face(fn);
    vv = f->vertex(ccw(i));
    vv->set_face(fn);
    fn->set_neighbor(in, Face_handle());
    hole.push_back(Edge(fn,in));
    to_delete.push_back(f);
  } while(fc != done);

  std::size_t size = to_delete.size();
  for(std::size_t i=0; i<size; i++) {
    f = to_delete[i];
    faces_set.erase(f);
    delete_face(f);
  }
}

template <class Gt, class Tds >
void
Triangulation_2<Gt, Tds>::
fill_hole(Vertex_handle v, std::list< Edge > & hole)
{
  // uses the fact that the hole is starshaped
  // with repect to v->point()
  typedef std::list<Edge> Hole;

  Face_handle ff, fn;
  int ii, in;
  Vertex_handle v0, v1, v2;
  Bounded_side side;

  //stack algorithm to create faces
  // create face v0,v1,v2
  //if v0,v1,v2 are finite vertices
  // and form a left_turn
  // and triangle v0v1v2 does not contain v->point()
  if( hole.size() != 3) {
    typename Hole::iterator hit = hole.begin();
    typename Hole::iterator next= hit;
    while( hit != hole.end() && hole.size() != 3) {
      ff = (*hit).first;
      ii = (*hit).second;
      v0 = ff->vertex(cw(ii));
      v1 = ff->vertex(ccw(ii));
      if( !is_infinite(v0) && !is_infinite(v1)) {
        next=hit; next++;
        if(next == hole.end()) next=hole.begin();
        fn = (*next).first;
        in = (*next).second;
        v2 = fn->vertex(ccw(in));
        if ( !is_infinite(v2) &&
             orientation(v0->point(), v1->point(), v2->point()) == LEFT_TURN ) {
          side = bounded_side(v0->point(),
                               v1->point(),
                               v2->point(),
                               v->point());

          if( side == ON_UNBOUNDED_SIDE ||
              (side == ON_BOUNDARY && orientation(v0->point(),
                                                  v->point(),
                                                  v2->point()) == COLLINEAR &&
               collinear_between(v0->point(),v->point(),v2->point()) ))
          {
            //create face
            Face_handle newf = create_face(ff,ii,fn,in);
            typename Hole::iterator tempo=hit;
            hit = hole.insert(hit,Edge(newf,1)); //push newf
            hole.erase(tempo); //erase ff
            hole.erase(next); //erase fn
            if (hit != hole.begin() ) --hit;
            continue;
          }
        }
      }
      ++hit;
    }
  }

  // either the hole has only three edges
  // or all its finite vertices are reflex or flat
  // except may be one vertex whose corresponding ear
  // includes the vertex being removed

  // deal with the last left_turn if any
  if(hole.size() != 3) {
    typename Hole::iterator hit=hole.begin();
    while(hit != hole.end()) {
      ff = (*hit).first;  ii = (*hit).second;
      hit++;
      if(hit != hole.end()) { fn = (*hit).first; in = (*hit).second;}
      else { fn = ((hole.front()).first); in = (hole.front()).second;}
      if ( !is_infinite(ff->vertex(cw(ii))) &&
           !is_infinite(fn->vertex(cw(in))) &&
           !is_infinite(fn->vertex(ccw(in))) &&
           orientation(ff->vertex(cw(ii))->point(),
                       fn->vertex(cw(in))->point(),
                       fn->vertex(ccw(in))->point()) == LEFT_TURN) {
        create_face(ff,ii,fn,in);
        break;
      }
    }
  }

  // deal with a reflex chain of convex hull edges
  if(hole.size() != 3) {
    // look for infinite vertex
    ff = (hole.front()).first;
    ii = (hole.front()).second;
    while ( ! is_infinite(ff->vertex(cw(ii)))){
      hole.push_back(hole.front());
      hole.pop_front();
      ff = (hole.front()).first;
      ii = (hole.front()).second;
    }
    //create faces
    while(hole.size() != 3){
      ff = (hole.front()).first;
      ii = (hole.front()).second;
      hole.pop_front();
      fn = (hole.front()).first;
      in = (hole.front()).second;
      hole.pop_front();
      Face_handle newf = create_face(ff,ii,fn,in);
      hole.push_front(Edge(newf,1));
    }
  }

  // now hole has three edges
  typename Hole::iterator hit;
  hit = hole.begin();
//  // I don't know why the following yelds a segmentation fault
//  create_face( (*hit).first, (*hit).second,
//               (* ++hit).first, (*hit).second,
//               (* ++hit).first, (*hit).second);
  ff = (*hit).first;      ii = (*hit).second;
  fn = (* ++hit).first;   in = (*hit).second;
  Face_handle f3 = (* ++hit).first;
  int i3 = (*hit).second;
  create_face(ff,ii,fn,in,f3,i3);
}

template < class Gt, class Tds >
template <class OutputItFaces>
void
Triangulation_2<Gt,Tds>::
fill_hole(Vertex_handle v, std::list<Edge> & hole, OutputItFaces fit)
{
  // uses the fact that the hole is starshaped
  // with repect to v->point()
  typedef std::list<Edge> Hole;

  Face_handle ff, fn;
  int ii , in;
  Vertex_handle v0, v1, v2;
  Bounded_side side;

  //stack algorithm to create faces
  // create face v0,v1,v2
  //if v0,v1,v2 are finite vertices
  // and form a left_turn
  // and triangle v0v1v2 does not contain v->point()
  if( hole.size() != 3) {
    typename Hole::iterator hit = hole.begin();
    typename Hole::iterator next= hit;
    while( hit != hole.end() && hole.size() != 3) {
      ff = (*hit).first;
      ii = (*hit).second;
      v0 = ff->vertex(cw(ii));
      v1 = ff->vertex(ccw(ii));
      if( !is_infinite(v0) && !is_infinite(v1)) {
        next=hit; next++;
        if(next == hole.end()) next=hole.begin();
        fn = (*next).first;
        in = (*next).second;
        v2 = fn->vertex(ccw(in));
        if ( !is_infinite(v2) &&
             orientation(v0->point(), v1->point(), v2->point()) == LEFT_TURN ) {
          side = bounded_side(v0->point(), v1->point(), v2->point(), v->point());
          if( side == ON_UNBOUNDED_SIDE ||
              (side == ON_BOUNDARY && orientation(v0->point(),
                                                  v->point(),
                                                  v2->point()) == COLLINEAR &&
               collinear_between(v0->point(),v->point(),v2->point()) ))
          {
            //create face
            Face_handle newf = create_face(ff,ii,fn,in);
            *fit++ = newf;
            typename Hole::iterator tempo=hit;
            hit = hole.insert(hit,Edge(newf,1)); //push newf
            hole.erase(tempo); //erase ff
            hole.erase(next); //erase fn
            if (hit != hole.begin() ) --hit;
            continue;
          }
        }
      }
      ++hit;
    }
  }

  // either the hole has only three edges
  // or all its finite vertices are reflex or flat
  // except may be one vertex whose corresponding ear
  // includes the vertex being removed

  // deal with the last left_turn if any
  if(hole.size() != 3) {
    typename Hole::iterator hit=hole.begin();
    while(hit != hole.end()) {
      ff = (*hit).first;
      ii = (*hit).second;
      hit++;
      if(hit != hole.end()) { fn = (*hit).first; in = (*hit).second;}
      else { fn = ((hole.front()).first); in = (hole.front()).second;}
      if ( !is_infinite(ff->vertex(cw(ii))) &&
           !is_infinite(fn->vertex(cw(in))) &&
           !is_infinite(fn->vertex(ccw(in))) &&
           orientation(ff->vertex(cw(ii))->point(),
                       fn->vertex(cw(in))->point(),
                       fn->vertex(ccw(in))->point()) == LEFT_TURN) {
        Face_handle  newf = create_face(ff,ii,fn,in);
        *fit++ = newf;
        break;
      }
    }
  }

  // deal with a reflex chain of convex hull edges
  if(hole.size() != 3) {
    // look for infinite vertex
    ff = (hole.front()).first;
    ii = (hole.front()).second;
    while ( ! is_infinite(ff->vertex(cw(ii)))){
      hole.push_back(hole.front());
      hole.pop_front();
      ff = (hole.front()).first;
      ii = (hole.front()).second;
    }
    //create faces
    while(hole.size() != 3){
      ff = (hole.front()).first;
      ii = (hole.front()).second;
      hole.pop_front();
      fn = (hole.front()).first;
      in = (hole.front()).second;
      hole.pop_front();
      Face_handle newf = create_face(ff,ii,fn,in);
      *fit++ = newf;
      hole.push_front(Edge(newf,1));
    }
  }

  // now hole has three edges
  typename Hole::iterator hit;
  hit = hole.begin();
//  // I don't know why the following yelds a segmentation fault
//  create_face( (*hit).first, (*hit).second,
//               (* ++hit).first, (*hit).second,
//               (* ++hit).first, (*hit).second);
  ff = (*hit).first;      ii = (*hit).second;
  fn = (* ++hit).first;   in = (*hit).second;
  Face_handle f3 = (* ++hit).first;
  int i3 = (*hit).second;
  Face_handle newf = create_face(ff,ii,fn,in,f3,i3);
  *fit++ = newf;
}

template <class Gt, class Tds >
void
Triangulation_2<Gt, Tds>::
fill_hole_delaunay(std::list<Edge> & first_hole)
{
  typedef std::list<Edge> Hole;
  typedef std::list<Hole> Hole_list;

  Face_handle f, ff, fn;
  int i, ii, in;
  Hole_list hole_list;

  hole_list.push_front(first_hole);

  while( ! hole_list.empty())
  {
    Hole& hole = hole_list.front();

    typename Hole::iterator hit = hole.begin();

    // if the hole has only three edges, create the triangle
    if (hole.size() == 3) {
      hit = hole.begin();
      f = (*hit).first;        i = (*hit).second;
      ff = (* ++hit).first;    ii = (*hit).second;
      fn = (* ++hit).first;    in = (*hit).second;
      create_face(f,i,ff,ii,fn,in);
      hole_list.pop_front();
      continue;
    }

    // else find an edge with two finite vertices
    // on the hole boundary
    // and the new triangle adjacent to that edge
    // cut the hole and push it back

    // first, ensure that a neighboring face
    // whose vertices on the hole boundary are finite
    // is the first of the hole
    bool finite= false;
    while (!finite){
      ff = (hole.front()).first;
      ii = (hole.front()).second;
      if ( is_infinite(ff->vertex(cw(ii))) ||
           is_infinite(ff->vertex(ccw(ii)))) {
        hole.push_back(hole.front());
        hole.pop_front();
      }
      else finite=true;
    }

    // take the first neighboring face and pop it;
    ff = (hole.front()).first;
    ii =(hole.front()).second;
    hole.pop_front();

    Vertex_handle v0 = ff->vertex(cw(ii));
    Vertex_handle v1 = ff->vertex(ccw(ii));
    Vertex_handle v2 = infinite_vertex();
    Vertex_handle v3;
    const Point& p0 = v0->point();
    const Point& p1 = v1->point();

    typename Hole::iterator hdone = hole.end();
    hit = hole.begin();
    typename Hole::iterator cut_after(hit);

    // if tested vertex is c with respect to the vertex opposite
    // to NULL neighbor,
    // stop at the before last face;
    hdone--;
    while( hit != hdone) {
      fn = (*hit).first;
      in = (*hit).second;
      Vertex_handle vv = fn->vertex(ccw(in));
      if (is_infinite(vv)) {
        if(is_infinite(v2)) cut_after = hit;
      }
      else {     // vv is a finite vertex
        const Point & p = vv->point();
        if (orientation(p0,p1,p) == COUNTERCLOCKWISE) {
          if (is_infinite(v2)) { v2=vv; v3=vv; cut_after=hit;}
          else{
            //
            if (this->side_of_oriented_circle(p0,p1,v3->point(),p,true) == ON_POSITIVE_SIDE){
              v2=vv; v3=vv; cut_after=hit;}
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
    if (fn->has_vertex(v2, i) && i == fn->ccw(in)) {
      newf = create_face(ff,ii,fn,in);
      hole.pop_front();
      hole.push_front(Edge( newf,1));
    }
    else{
      fn = (hole.back()).first;
      in = (hole.back()).second;
      if (fn->has_vertex(v2, i) && i== fn->cw(in)) {
        newf = create_face(fn,in,ff,ii);
        hole.pop_back();
        hole.push_back(Edge(newf,1));
      }
      else {
        // split the hole in two holes
        newf = create_face(ff,ii,v2);
        Hole new_hole;
        ++cut_after;
        while( hole.begin() != cut_after )
        {
          new_hole.push_back(hole.front());
          hole.pop_front();
        }

        hole.push_front(Edge( newf,1));
        new_hole.push_front(Edge( newf,0));
        hole_list.push_front(new_hole);
      }
    }
  }
}

template < class Gt, class Tds >
template <class OutputItFaces>
void
Triangulation_2<Gt,Tds>::
fill_hole_delaunay(std::list<Edge> & first_hole, OutputItFaces fit)
{
  typedef typename Gt::Side_of_oriented_circle_2 In_circle;
  typedef std::list<Edge>                        Hole;
  typedef std::list<Hole>                        Hole_list;

  In_circle in_circle = geom_traits().side_of_oriented_circle_2_object();

  Face_handle f, ff, fn;
  int i, ii, in;
  Hole_list hole_list;
  hole_list.push_front(first_hole);

  while(!hole_list.empty()) {
    Hole& hole = hole_list.front();
    typename Hole::iterator hit = hole.begin();

    if (hole.size() == 3) {
      hit = hole.begin();
      f = (*hit).first;        i = (*hit).second;
      ff = (* ++hit).first;    ii = (*hit).second;
      fn = (* ++hit).first;    in = (*hit).second;
      Face_handle newf = create_face(f,i,ff,ii,fn,in);
      *fit++ = newf;
      hole_list.pop_front();
      continue;
    }

    bool finite= false;
    while (!finite){
      ff = (hole.front()).first;
      ii = (hole.front()).second;
      if ( is_infinite(ff->vertex(cw(ii))) ||
     is_infinite(ff->vertex(ccw(ii)))) {
        hole.push_back(hole.front());
        hole.pop_front();
      } else finite=true;
    }

    ff = (hole.front()).first;
    ii =(hole.front()).second;
    hole.pop_front();

    Vertex_handle v0 = ff->vertex(cw(ii));
    Vertex_handle v1 = ff->vertex(ccw(ii));
    Vertex_handle v2 = infinite_vertex();
    const Point& p0 = v0->point();
    const Point& p1 = v1->point();

    typename Hole::iterator hdone = hole.end();
    hit = hole.begin();
    typename Hole::iterator cut_after(hit);

    hdone--;
    while( hit != hdone) {
      fn = (*hit).first;
      in = (*hit).second;
      Vertex_handle vv = fn->vertex(ccw(in));
      if (is_infinite(vv)) {
        if(is_infinite(v2)) cut_after = hit;
      } else {     // vv is a finite vertex
        const Point & p = vv->point();
        if (orientation(p0,p1,p) == CGAL::COUNTERCLOCKWISE) {
          if (is_infinite(v2)) { v2 = vv; cut_after = hit;}
          else{
            if (in_circle(p0,p1,v2->point(),p) == CGAL::ON_POSITIVE_SIDE){
              v2 = vv; cut_after = hit;
            }
          }
        }
      }
      ++hit;
    }

    Face_handle newf;

    fn = (hole.front()).first;
    in = (hole.front()).second;
    if (fn->has_vertex(v2, i) && i == fn->ccw(in)) {
      newf = create_face(ff,ii,fn,in);
      hole.pop_front();
      hole.push_front(Edge( newf,1));
    } else {
      fn = (hole.back()).first;
      in = (hole.back()).second;
      if (fn->has_vertex(v2, i) && i== fn->cw(in)) {
        newf = create_face(fn,in,ff,ii);
        hole.pop_back();
        hole.push_back(Edge(newf,1));
      } else {
        newf = create_face(ff,ii,v2);
        Hole new_hole;
        ++cut_after;
        while( hole.begin() != cut_after ) {
          new_hole.push_back(hole.front());
          hole.pop_front();
        }
        hole.push_front(Edge(newf, 1));
        new_hole.push_front(Edge(newf, 0));
        hole_list.push_front(new_hole);
      }
    }
    *fit++ = newf;
  }
}

template <class Gt, class Tds >
typename Triangulation_2<Gt,Tds>::Vertex_handle
Triangulation_2<Gt,Tds>::
move_if_no_collision(Vertex_handle v, const Point &p)
{
  CGAL_triangulation_precondition(!is_infinite(v));
  if(v->point() == p) return v;

  const int dim = dimension();

  Locate_type lt;
  int li;
  Face_handle loc = locate(p, lt, li, v->face());

  if(lt == VERTEX) return loc->vertex(li);

  if(dim == 0) {
    v->set_point(p);
    return v;
  }

  size_type n_vertices = tds().number_of_vertices();

  if((lt == OUTSIDE_AFFINE_HULL) && (dim == 1) && (n_vertices == 3)) {
    v->set_point(p);
    return v;
  }

  if((lt != OUTSIDE_AFFINE_HULL) && (dim == 1)) {
    if(loc->has_vertex(v)) {
      v->set_point(p);
    } else {
      Vertex_handle inserted = insert(p, lt, loc, li);
      Face_handle f = v->face();
      int i = f->index(v);
      if (i==0) {f = f->neighbor(1);}
      CGAL_triangulation_assertion(f->index(v) == 1);
      Face_handle g= f->neighbor(0);
      f->set_vertex(1, g->vertex(1));
      f->set_neighbor(0,g->neighbor(0));
      g->neighbor(0)->set_neighbor(1,f);
      g->vertex(1)->set_face(f);
      delete_face(g);
      Face_handle f_ins = inserted->face();
      i = f_ins->index(inserted);
      if (i==0) {f_ins = f_ins->neighbor(1);}
      CGAL_triangulation_assertion(f_ins->index(inserted) == 1);
      Face_handle g_ins = f_ins->neighbor(0);
      f_ins->set_vertex(1, v);
      g_ins->set_vertex(0, v);
      v->set_point(p);
      v->set_face(inserted->face());
      delete_vertex(inserted);
    }
    return v;
  }

  if((lt != OUTSIDE_AFFINE_HULL) && test_dim_down(v)) {
    // verify if p and two static vertices are collinear in this case
    int iinf = 0;
    Face_circulator finf = incident_faces(infinite_vertex()), fdone(finf);
    do {
      if(!finf->has_vertex(v))
      {
        iinf = ~(finf->index(infinite_vertex()));
        break;
      }
    } while(++finf != fdone);
    if(this->orientation(finf->vertex(iinf&1)->point(),
                         finf->vertex(iinf&2)->point(),
                         p) == COLLINEAR)
    {
      v->set_point(p);
      _tds.dim_down(loc, loc->index(v));
      return v;
    }
  }

  Vertex_handle inserted = insert(p, lt, loc, li);

  std::list<Edge> hole;
  make_hole(v, hole);
  fill_hole(v, hole);

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
  delete_vertex(inserted);

  return v;
}

template <class Gt, class Tds >
typename Triangulation_2<Gt,Tds>::Vertex_handle
Triangulation_2<Gt,Tds>::
move(Vertex_handle v, const Point &p)
{
  CGAL_triangulation_precondition(!is_infinite(v));
  if(v->point() == p) return v;
  Vertex_handle w = move_if_no_collision(v,p);
  if(w != v) {
    remove(v);
    return w;
  }
  return v;
}

template <class Gt, class Tds >
template <class OutputItFaces>
typename Triangulation_2<Gt,Tds>::Vertex_handle
Triangulation_2<Gt,Tds>::
move_if_no_collision_and_give_new_faces(Vertex_handle v,
                                            const Point &p,
                                        OutputItFaces oif)
{
  CGAL_triangulation_precondition(!is_infinite(v));
  if(v->point() == p) return v;
  const int dim = this->dimension();

  Locate_type lt;
  int li;
  Vertex_handle inserted;
  Face_handle loc = locate(p, lt, li, v->face());

  if(lt == VERTEX) return loc->vertex(li);

  if(dim == 0) {
    v->set_point(p);
    return v;
  }

  int n_vertices = tds().number_of_vertices();

  if((lt == OUTSIDE_AFFINE_HULL) && (dim == 1) && (n_vertices == 3)) {
    v->set_point(p);

    for(All_faces_iterator afi = tds().face_iterator_base_begin();
        afi != tds().face_iterator_base_begin();
        afi++) *oif++ = afi;

    return v;
  }

  if((lt != OUTSIDE_AFFINE_HULL) && (dim == 1)) {
    if(loc->has_vertex(v)) {
      v->set_point(p);
    } else {
      inserted = insert(p, lt, loc, li);
      Face_handle f = v->face();
      int i = f->index(v);
      if (i==0) {f = f->neighbor(1);}
      CGAL_triangulation_assertion(f->index(v) == 1);
      Face_handle g= f->neighbor(0);
      f->set_vertex(1, g->vertex(1));
      f->set_neighbor(0,g->neighbor(0));
      g->neighbor(0)->set_neighbor(1,f);
      g->vertex(1)->set_face(f);
      delete_face(g);
      *oif++ = f;
      Face_handle f_ins = inserted->face();
      i = f_ins->index(inserted);
      if (i==0) {f_ins = f_ins->neighbor(1);}
      CGAL_triangulation_assertion(f_ins->index(inserted) == 1);
      Face_handle g_ins = f_ins->neighbor(0);
      f_ins->set_vertex(1, v);
      g_ins->set_vertex(0, v);
      v->set_point(p);
      v->set_face(inserted->face());
      delete_vertex(inserted);
    }
    *oif++ = v->face();
    if(v->face()->neighbor(0)->has_vertex(v))
      *oif++ = v->face()->neighbor(0);
    if(v->face()->neighbor(1)->has_vertex(v))
      *oif++ = v->face()->neighbor(1);
    return v;
  }

  if((lt != OUTSIDE_AFFINE_HULL) && test_dim_down(v)) {
    // verify if p and two static vertices are collinear in this case
    int iinf;
    Face_circulator finf = incident_faces(infinite_vertex()), fdone(finf);
    do {
      if(!finf->has_vertex(v))
      {
        iinf = ~(finf->index(infinite_vertex()));
        break;
      }
    } while (++finf != fdone);
    if(this->orientation(finf->vertex(iinf&1)->point(),
                         finf->vertex(iinf&2)->point(),
                         p) == COLLINEAR)
    {
      v->set_point(p);
      _tds.dim_down(loc, loc->index(v));
      return v;
    }

    for(All_faces_iterator afi = tds().face_iterator_base_begin();
        afi != tds().face_iterator_base_begin();
        afi++)
      *oif++ = afi;
  }

  std::set<Face_handle> faces_set;
  inserted = insert(p, lt, loc, li);
  Face_circulator fc = incident_faces(inserted), done(fc);
  do { faces_set.insert(fc); } while (++fc != done);

  std::list<Edge> hole;
  make_hole(v, hole, faces_set);
  fill_hole(v, hole, oif);

  fc = this->incident_faces(inserted), done(fc);
  std::vector<Face_handle> faces_pt;
  faces_pt.reserve(16);
  do { faces_pt.push_back(fc); } while (++fc != done);
  int ss = faces_pt.size();
  for(int k=0; k<ss; k++)
  {
    Face_handle f = faces_pt[k];
    int i = f->index(inserted);
    f->set_vertex(i, v);
  }
  v->set_point(p);
  v->set_face(inserted->face());
  delete_vertex(inserted);

  for(typename std::set<Face_handle>::iterator ib = faces_set.begin(),
      iend = faces_set.end(); ib != iend; ib++)
    *oif++ = *ib;

  return v;
}

template <class Gt, class Tds >
inline
typename Triangulation_2<Gt, Tds>::Face_handle
Triangulation_2<Gt, Tds>::
create_face(Face_handle f1, int i1,
            Face_handle f2, int i2,
            Face_handle f3, int i3)
{
  return _tds.create_face(f1, i1, f2, i2, f3, i3);
}

template <class Gt, class Tds >
inline
typename Triangulation_2<Gt, Tds>::Face_handle
Triangulation_2<Gt, Tds>::
create_face(Face_handle f1, int i1,
            Face_handle f2, int i2)
{
  return _tds.create_face(f1, i1, f2, i2);
}

template <class Gt, class Tds >
inline
typename Triangulation_2<Gt, Tds>::Face_handle
Triangulation_2<Gt, Tds>::
create_face(Face_handle f, int i, Vertex_handle v)
{
  return _tds.create_face(f, i, v);
}

template <class Gt, class Tds >
inline
typename Triangulation_2<Gt, Tds>::Face_handle
Triangulation_2<Gt, Tds>::
create_face(Vertex_handle v1, Vertex_handle v2, Vertex_handle v3)
{
  return _tds.create_face(v1, v2, v3);
}

template <class Gt, class Tds >
inline
typename Triangulation_2<Gt, Tds>::Face_handle
Triangulation_2<Gt, Tds>::
create_face(Vertex_handle v1, Vertex_handle v2, Vertex_handle v3,
            Face_handle f1, Face_handle f2, Face_handle f3)
{
  return _tds.create_face(v1, v2, v3, f1, f2, f3);
}

template <class Gt, class Tds >
inline
typename Triangulation_2<Gt, Tds>::Face_handle
Triangulation_2<Gt, Tds>::
create_face(Face_handle fh)
{
  return _tds.create_face(fh);
}

template <class Gt, class Tds >
inline
typename Triangulation_2<Gt, Tds>::Face_handle
Triangulation_2<Gt, Tds>::
create_face()
{
  return _tds.create_face();
}

template <class Gt, class Tds >
inline
void
Triangulation_2<Gt, Tds>::
delete_face(Face_handle f)
{
  _tds.delete_face(f);
}

template <class Gt, class Tds >
inline
void
Triangulation_2<Gt, Tds>::
delete_vertex(Vertex_handle v)
{
  _tds.delete_vertex(v);
}

// POINT LOCATION
template <class Gt, class Tds >
typename Triangulation_2<Gt, Tds>::Face_handle
Triangulation_2<Gt, Tds>::
march_locate_1D(const Point& t,
                Locate_type& lt,
                int& li) const
{
  Face_handle ff = infinite_face();
  int iv = ff->index(infinite_vertex());
  Face_handle f = ff->neighbor(iv);
  Orientation pqt = orientation(f->vertex(0)->point(),
                                f->vertex(1)->point(),
                                t);
  if(pqt == RIGHT_TURN || pqt == LEFT_TURN) {
    lt = OUTSIDE_AFFINE_HULL;
    li = 4 ;// should not be used
    return Face_handle();
  }

  int i= f->index(ff);
  if (collinear_between(t,f->vertex(1-i)->point(),f->vertex(i)->point())) {
    lt = OUTSIDE_CONVEX_HULL;
    li = iv;
    return ff;
  }

  if( xy_equal(t,f->vertex(1-i)->point()) ){
    lt = VERTEX;
    li=1-i;
    return f;
  }

  ff = ff->neighbor(1-iv); //the other infinite face
  iv = ff->index(infinite_vertex());
  f = ff->neighbor(iv);
  i = f->index(ff);
  if (collinear_between(t,f->vertex(1-i)->point(),f->vertex(i)->point())) {
    lt = OUTSIDE_CONVEX_HULL;
    li = iv;
    return ff;
  }
  if( xy_equal(t,f->vertex(1-i)->point()) ){
    lt = VERTEX;
    li=1-i;
    return f;
  }

  Finite_edges_iterator eit= finite_edges_begin();
  Vertex_handle u,v;
  for( ; eit != finite_edges_end() ; eit++) {
    u = (*eit).first->vertex(0);
    v = (*eit).first->vertex(1);
    if(xy_equal(t,v->point())){
      lt = VERTEX;
      li = 1;
      return (*eit).first;
    }
    if(collinear_between(u->point(), t, v->point())){
      lt = EDGE;
      li = 2;
      return (*eit).first;
    }
  }
  CGAL_triangulation_assertion(false);
  return Face_handle();
}

template <class Gt, class Tds >
typename Triangulation_2<Gt, Tds>::Face_handle
Triangulation_2<Gt, Tds>::
march_locate_2D_LFC(Face_handle start,
                    const Point& t,
                    Locate_type& lt,
                    int& li) const
{
  //    CGAL_triangulation_precondition( ! is_infinite(start) );
  const Point& p = start->vertex(0)->point();
  const Point& q = start->vertex(1)->point();
  const Point& r = start->vertex(2)->point();
  if(xy_equal(t,p)) {
    lt = VERTEX;
    li = 0;
    return start;
  }

  Line_face_circulator lfc;

  Orientation o2 = orientation(p, q, t);
  Orientation o0 = orientation(q, r, t);
  Orientation o1 = orientation(r, p, t);
  if( (o2 == LEFT_TURN)&& (o1 == LEFT_TURN)) {
    lfc = Line_face_circulator(start, 0,
                               Line_face_circulator::vertex_edge,
                               this, p, t);
  } else if ( (o0 == LEFT_TURN)&& (o2 == LEFT_TURN)) {
    lfc = Line_face_circulator(start, 1,
                               Line_face_circulator::vertex_edge,
                               this, q, t);
  } else if ( (o1 == LEFT_TURN)&& (o0 == LEFT_TURN)) {
    lfc = Line_face_circulator(start, 2,
                               Line_face_circulator::vertex_edge,
                               this, r, t);
  } if( (o2 == RIGHT_TURN)&& (o1 == RIGHT_TURN)) {
    lfc = Line_face_circulator(start, 0,
                               Line_face_circulator::edge_vertex,
                               this, p, t);
  } else if ( (o0 == RIGHT_TURN)&& (o2 == RIGHT_TURN)) {
    lfc = Line_face_circulator(start, 1,
                               Line_face_circulator::edge_vertex,
                               this, q, t);
  } else if ( (o1 == RIGHT_TURN)&& (o0 == RIGHT_TURN)) {
    lfc = Line_face_circulator(start, 2,
                               Line_face_circulator::edge_vertex,
                               this, r, t);
  }else {
    lfc = Line_face_circulator(start->vertex(0), this, t);
  }
  if(lfc==0 || lfc.collinear_outside()){
    // point t lies outside or on the convex hull
    // we walk on the convex hull to find it out
    Face_circulator fc = incident_faces(infinite_vertex());
    Face_circulator done(fc);
    int ic = fc->index(infinite_vertex());
    if (xy_equal(t,fc->vertex(cw(ic))->point())){
      lt = VERTEX;
      li = cw(ic);
      return fc;
    }
    Orientation ori;
    do{ // walking ccw around convex hull
      ic = fc->index(infinite_vertex());
      if (xy_equal(t,fc->vertex(ccw(ic))->point())){
        lt = VERTEX;
        li = ccw(ic);
        return fc;
      }
      ori = orientation( fc->vertex(cw(ic))->point(),
                         fc->vertex(ccw(ic))->point(), t);
      if (ori == RIGHT_TURN) {
        lt = OUTSIDE_CONVEX_HULL;
        li = ic;
        return fc;
      }
      if (ori == COLLINEAR &&
          collinear_between(fc->vertex(cw(ic))->point(),
                            t,
                            fc->vertex(ccw(ic))->point()) ) {
        lt = EDGE;
        li = ic;
        return fc;
      }
    } while (--fc != done);
    //should not arrive there;
    CGAL_triangulation_assertion(fc != done);
  }

  while(! lfc.locate(t, lt, li) ){
    ++lfc;
  }
  return lfc;
}

template <class Gt, class Tds >
void
Triangulation_2<Gt, Tds>::
compare_walks(const Point& p,
              Face_handle c1, Face_handle c2,
              Locate_type& lt1, Locate_type& lt2,
              int li1, int li2) const
{
  bool b = true;
  b = b && (lt1 == lt2);
  if((lt1 == lt2) && (lt1 == VERTEX)) {
    b = b && ( c1->vertex(li1) == c2->vertex(li2) );
  } else if((lt1 == lt2) && (lt1 == EDGE)) {
    b = b && ((c1 == c2) || ( (c1->neighbor(li1) == c2) && (c2->neighbor(li2) == c1)));
  } else if((lt1 == lt2) && (lt1 == OUTSIDE_CONVEX_HULL)) {
    b = b && (is_infinite(c1) && is_infinite(c2));
  } else {
    b = b && (lt1 == lt2);
    b = b && (lt1 == FACE);
    b = b && (c1 == c2);
  }

  if ( c1 != c2) {
    std::cerr << "from compare_walks " << std::endl;
    std::cerr << "point " << p << std::endl;
    std::cerr << "locate 1 " << &*c1 << "\t" << lt1 << "\t" << li1
              << std::endl;
    std::cerr << "locate 2 " << &*c2 << "\t" << lt2 << "\t" << li2
              << std::endl;
    std::cerr << std::endl;
    show_face(c1);
    std::cerr << std::endl;
    show_face(c2);
    std::cerr << std::endl;
  }
  CGAL_triangulation_assertion(b);
}


#if 1
template <class Gt, class Tds >
typename Triangulation_2<Gt, Tds>::Face_handle
Triangulation_2<Gt, Tds>::
march_locate_2D(Face_handle c,
                const Point& t,
                Locate_type& lt,
                int& li) const
{
  CGAL_triangulation_assertion(! is_infinite(c));

  boost::rand48 rng;

  boost::uniform_smallint<> two(0, 1);
  boost::variate_generator<boost::rand48&, boost::uniform_smallint<> > coin(rng, two);

  Face_handle prev = Face_handle();
  bool first = true;
  while (1) {
    if ( is_infinite(c) ) {
      // c must contain t in its interior
      lt = OUTSIDE_CONVEX_HULL;
      li = c->index(infinite_vertex());
      return c;
    }

    // else c is finite
    // Instead of testing its edges in a random order we do the following
    // until we find a neighbor to go further:
    // As we come from prev we do not have to check the edge leading to prev
    // Now we flip a coin in order to decide if we start checking the
    // edge before or the edge after the edge leading to prev
    // We do loop unrolling in order to find out if this is faster.
    // In the very beginning we do not have a prev, but for the first step
    // we do not need randomness
    int left_first = coin()%2;

    const Point & p0 = c->vertex( 0 )->point();
    const Point & p1 = c->vertex( 1 )->point();
    const Point & p2 = c->vertex( 2 )->point();
    Orientation o0, o1, o2;

    if(first){
      prev = c;
      first = false;
      o0 = orientation(p0,p1,t);
      if ( o0 == NEGATIVE ) {
        c = c->neighbor( 2 );
        continue;
      }
      o1 = orientation(p1,p2,t);
      if ( o1 == NEGATIVE ) {
        c = c->neighbor( 0 );
        continue;
      }
      o2 = orientation(p2,p0,t);
      if ( o2 == NEGATIVE ) {
        c = c->neighbor( 1 );
        continue;
      }
    } else if(left_first){
      if(c->neighbor(0) == prev){
        prev = c;
        o0 = orientation(p0,p1,t);
        if ( o0 == NEGATIVE ) {
          c = c->neighbor( 2 );
          continue;
        }
        o2 = orientation(p2,p0,t);
        if ( o2 == NEGATIVE ) {
          c = c->neighbor( 1 );
          continue;
        }
        o1 = POSITIVE;
      } else if(c->neighbor(1) == prev){
        prev = c;
        o1 = orientation(p1,p2,t);
        if ( o1 == NEGATIVE ) {
          c = c->neighbor( 0 );
          continue;
        }
        o0 = orientation(p0,p1,t);
        if ( o0 == NEGATIVE ) {
          c = c->neighbor( 2 );
          continue;
        }
        o2 = POSITIVE;
      } else {
        prev = c;
        o2 = orientation(p2,p0,t);
        if ( o2 == NEGATIVE ) {
          c = c->neighbor( 1 );
          continue;
        }
        o1 = orientation(p1,p2,t);
        if ( o1 == NEGATIVE ) {
          c = c->neighbor( 0 );
          continue;
        }
        o0 = POSITIVE;
      }

    } else { // right_first
      if(c->neighbor(0) == prev){
        prev = c;
        o2 = orientation(p2,p0,t);
        if ( o2 == NEGATIVE ) {
          c = c->neighbor( 1 );
          continue;
        }
        o0 = orientation(p0,p1,t);
        if ( o0 == NEGATIVE ) {
          c = c->neighbor( 2 );
          continue;
        }
        o1 = POSITIVE;
      } else if(c->neighbor(1) == prev){
        prev = c;
        o0 = orientation(p0,p1,t);
        if ( o0 == NEGATIVE ) {
          c = c->neighbor( 2 );
          continue;
        }
        o1 = orientation(p1,p2,t);
        if ( o1 == NEGATIVE ) {
          c = c->neighbor( 0 );
          continue;
        }
        o2 = POSITIVE;
      } else {
        prev = c;
        o1 = orientation(p1,p2,t);
        if ( o1 == NEGATIVE ) {
          c = c->neighbor( 0 );
          continue;
        }
        o2 = orientation(p2,p0,t);
        if ( o2 == NEGATIVE ) {
          c = c->neighbor( 1 );
          continue;
        }
        o0 = POSITIVE;
      }
    }

    // now p is in c or on its boundary
    int sum = ( o0 == COLLINEAR )
              + ( o1 == COLLINEAR )
              + ( o2 == COLLINEAR );
    switch (sum) {
      case 0:
      {
        lt = FACE;
        li = 4;
        break;
      }
      case 1:
      {
        lt = EDGE;
        li = ( o0 == COLLINEAR ) ? 2 :
                                   ( o1 == COLLINEAR ) ? 0 :
                                                         1;
        break;
      }
      case 2:
      {
        lt = VERTEX;
        li = ( o0 != COLLINEAR ) ? 2 :
                                   ( o1 != COLLINEAR ) ? 0 :
                                                         1;
        break;
      }
    }
    return c;
  }
}
#else // not 1
template <class Gt, class Tds >
typename Triangulation_2<Gt, Tds>::Face_handle
Triangulation_2<Gt, Tds>::
march_locate_2D(Face_handle c,
                const Point& t,
                Locate_type& lt,
                int& li) const
{
  CGAL_triangulation_assertion(! is_infinite(c));

  boost::uniform_smallint<> three(0, 2);
  boost::variate_generator<boost::rand48&, boost::uniform_smallint<> > die3(rng, three);

  Face_handle prev = Face_handle();
  while (1) {
    if ( is_infinite(c) ) {
      // c must contain t in its interior
      lt = OUTSIDE_CONVEX_HULL;
      li = c->index(infinite_vertex());
      return c;
    }

    // else c is finite
    // we test its edges in a random order until we find a
    // neighbor to go further

    int i = die3();
    int ccwi = ccw(i);
    int cwi = cw(i);
    const Point & p0 = c->vertex( i )->point();
    const Point & p1 = c->vertex( ccwi )->point();
    Orientation o0, o1, o2;
    CGAL_triangulation_assertion(orientation(p0,p1,c->vertex( cwi )->point())==POSITIVE);
    if(c->neighbor(cwi) == prev){
      o0 = POSITIVE;
    } else {
      o0 = orientation(p0,p1,t);
      if ( o0 == NEGATIVE ) {
        prev = c;
        c = c->neighbor( cwi );
        continue;
      }
    }
    const Point & p2 = c->vertex( cwi )->point();
    if(c->neighbor(i) == prev){
      o1 = POSITIVE;
    } else {
      o1 = orientation(p1,p2,t);
      if ( o1 == NEGATIVE ) {
        prev = c;
        c = c->neighbor( i );
        continue;
      }
    }
    if(c->neighbor(ccwi) == prev){
      o2 = POSITIVE;
    } else {
      o2 = orientation(p2,p0,t);
      if ( o2 == NEGATIVE ) {
        prev = c;
        c = c->neighbor( ccwi );
        continue;
      }
    }

    // now p is in c or on its boundary
    int sum = ( o0 == COLLINEAR )
              + ( o1 == COLLINEAR )
              + ( o2 == COLLINEAR );
    switch (sum) {
      case 0:
      {
        lt = FACE;
        li = 4;
        break;
      }
      case 1:
      {
        lt = EDGE;
        li = ( o0 == COLLINEAR ) ? cwi :
                                   ( o1 == COLLINEAR ) ? i :
                                                         ccwi;
        break;
      }
      case 2:
      {
        lt = VERTEX;
        li = ( o0 != COLLINEAR ) ? cwi :
                                   ( o1 != COLLINEAR ) ? i :
                                                         ccwi;
        break;
      }
    }
    return c;
  }
}
#endif // not 1

template <class Gt, class Tds >
typename Triangulation_2<Gt, Tds>::Face_handle
Triangulation_2<Gt,Tds>::
#ifdef CGAL_NO_STRUCTURAL_FILTERING
locate(const Point& p,
       Locate_type& lt,
       int& li,
       Face_handle start) const
#else // no CGAL_NO_STRUCTURAL_FILTERING
exact_locate(const Point& p,
             Locate_type& lt,
             int& li,
             Face_handle start) const
#endif // no CGAL_NO_STRUCTURAL_FILTERING
{
  if (dimension() < 0) {
    lt = OUTSIDE_AFFINE_HULL;
    li = 4; // li should not be used in this case
    return Face_handle();
  }
  if( dimension() == 0) {
    // Do not use finite_vertex directly because there can be hidden vertices
    // (regular triangulations)
    if (xy_equal(p,finite_vertex()->face()->vertex(0)->point())){
      lt = VERTEX ;
    }
    else{
      lt = OUTSIDE_AFFINE_HULL;
    }
    li = 4; // li should not be used in this case
    return Face_handle();
  }
  if(dimension() == 1){
    return march_locate_1D(p, lt, li);
  }

  if(start == Face_handle()) {
    start = infinite_face()->
            neighbor(infinite_face()->index(infinite_vertex()));
  } else if(is_infinite(start)) {
    start = start->neighbor(start->index(infinite_vertex()));
  }

#if ( ! defined(CGAL_ZIG_ZAG_WALK)) && ( ! defined(CGAL_LFC_WALK))
#define CGAL_ZIG_ZAG_WALK
#endif

#ifdef CGAL_ZIG_ZAG_WALK
  Face_handle res1 = march_locate_2D(start, p, lt, li);
#endif
#ifdef CGAL_LFC_WALK
  Locate_type lt2;
  int li2;
  Face_handle res2 = march_locate_2D_LFC(start, p, lt2, li2);
#endif

#if defined(CGAL_ZIG_ZAG_WALK) && defined(CGAL_LFC_WALK)
  compare_walks(p,
                res1, res2,
                lt, lt2,
                li, li2);
#endif

#ifdef CGAL_ZIG_ZAG_WALK
  return res1;
#endif

#ifdef CGAL_LFC_WALK
  lt = lt2;
  li = li2;
  return res2;
#endif
}

#ifdef CGAL_NO_STRUCTURAL_FILTERING
template <class Gt, class Tds >
typename Triangulation_2<Gt, Tds>:: Face_handle
Triangulation_2<Gt, Tds>::
locate(const Point &p,
       Face_handle start) const
{
  Locate_type lt;
  int li;
  return locate(p, lt, li, start);
}
#else
template <class Gt, class Tds >
inline
typename Triangulation_2<Gt, Tds>::Face_handle
Triangulation_2<Gt, Tds>::
inexact_locate(const Point & t, Face_handle start, int n_of_turns) const
{
  if(dimension() < 2) return start;

  if(start == Face_handle()){
    start = infinite_face()->
      neighbor(infinite_face()->index(infinite_vertex()));
  } else if(is_infinite(start)){
    start = start->neighbor(start->index(infinite_vertex()));
  }

  Face_handle prev = Face_handle(), c = start;
  bool first = true;
  while (1) {

    if(!(n_of_turns--)) return c;

    if ( is_infinite(c) ) return c;

    const Point & p0 = c->vertex( 0 )->point();
    const Point & p1 = c->vertex( 1 )->point();
    const Point & p2 = c->vertex( 2 )->point();

    if(first) {
      prev = c;
      first = false;
      if(has_inexact_negative_orientation(p0,p1,t) ) {
        c = c->neighbor( 2 );
        continue;
      }
      if(has_inexact_negative_orientation(p1,p2,t) ) {
        c = c->neighbor( 0 );
        continue;
      }
      if (has_inexact_negative_orientation(p2,p0,t) ) {
        c = c->neighbor( 1 );
        continue;
      }
    } else {
      if(c->neighbor(0) == prev){
        prev = c;
        if (has_inexact_negative_orientation(p0,p1,t) ) {
          c = c->neighbor( 2 );
          continue;
        }
        if (has_inexact_negative_orientation(p2,p0,t) ) {
          c = c->neighbor( 1 );
          continue;
        }
      } else if(c->neighbor(1) == prev){
        prev = c;
        if (has_inexact_negative_orientation(p0,p1,t) ) {
          c = c->neighbor( 2 );
          continue;
        }
        if (has_inexact_negative_orientation(p1,p2,t) ) {
          c = c->neighbor( 0 );
          continue;
        }
      } else {
        prev = c;
        if (has_inexact_negative_orientation(p2,p0,t) ) {
          c = c->neighbor( 1 );
          continue;
        }
        if (has_inexact_negative_orientation(p1,p2,t) ) {
          c = c->neighbor( 0 );
          continue;
        }
      }
    }
    break;
  }
  return c;
}

template <class Gt, class Tds >
inline
bool
Triangulation_2<Gt, Tds>::
has_inexact_negative_orientation(const Point &p, const Point &q,
                                 const Point &r) const
{
  // So that this code works well with Lazy_kernel
  internal::Static_filters_predicates::Get_approx<Point> get_approx;

  const double px = to_double(get_approx(p).x());
  const double py = to_double(get_approx(p).y());
  const double qx = to_double(get_approx(q).x());
  const double qy = to_double(get_approx(q).y());
  const double rx = to_double(get_approx(r).x());
  const double ry = to_double(get_approx(r).y());

  const double pqx = qx - px;
  const double pqy = qy - py;
  const double prx = rx - px;
  const double pry = ry - py;

  return ( determinant(pqx, pqy, prx, pry) < 0);
}
#endif

template <class Gt, class Tds >
typename Triangulation_2<Gt, Tds>::Finite_faces_iterator
Triangulation_2<Gt, Tds>::
finite_faces_begin() const
{
  if ( dimension() < 2 )
    return finite_faces_end();
  return CGAL::filter_iterator( all_faces_end(),
                                Infinite_tester(this),
                                all_faces_begin() );
}

template <class Gt, class Tds >
typename Triangulation_2<Gt, Tds>::Finite_faces_iterator
Triangulation_2<Gt, Tds>::
finite_faces_end() const
{
  return CGAL::filter_iterator( all_faces_end(),
                                Infinite_tester(this) );
}

template <class Gt, class Tds >
typename Triangulation_2<Gt, Tds>::Finite_vertices_iterator
Triangulation_2<Gt, Tds>::
finite_vertices_begin() const
{
  if ( number_of_vertices() <= 0 )
    return finite_vertices_end();
  return CGAL::filter_iterator( all_vertices_end(),
                                Infinite_tester(this),
                                all_vertices_begin() );
}

template <class Gt, class Tds >
typename Triangulation_2<Gt, Tds>::Finite_vertices_iterator
Triangulation_2<Gt, Tds>::
finite_vertices_end() const
{
  return CGAL::filter_iterator(all_vertices_end(),
                               Infinite_tester(this));
}

template <class Gt, class Tds >
typename Triangulation_2<Gt, Tds>::Finite_edges_iterator
Triangulation_2<Gt, Tds>::
finite_edges_begin() const
{
  if ( dimension() < 1 )
    return finite_edges_end();
  return CGAL::filter_iterator( all_edges_end(),
                                infinite_tester(),
                                all_edges_begin());
}

template <class Gt, class Tds >
typename Triangulation_2<Gt, Tds>::Finite_edges_iterator
Triangulation_2<Gt, Tds>::
finite_edges_end() const
{
  return CGAL::filter_iterator(all_edges_end(),
                               infinite_tester() );
}

template <class Gt, class Tds >
typename Triangulation_2<Gt, Tds>::Point_iterator
Triangulation_2<Gt, Tds>::
points_begin() const
{
  return Point_iterator(finite_vertices_begin());
}

template <class Gt, class Tds >
typename Triangulation_2<Gt, Tds>::Point_iterator
Triangulation_2<Gt, Tds>::
points_end() const
{
  return Point_iterator(finite_vertices_end());
}

template <class Gt, class Tds >
typename Triangulation_2<Gt, Tds>::All_faces_iterator
Triangulation_2<Gt, Tds>::
all_faces_begin() const
{
  return _tds.faces_begin();
}

template <class Gt, class Tds >
typename Triangulation_2<Gt, Tds>::All_faces_iterator
Triangulation_2<Gt, Tds>::
all_faces_end() const
{
  return _tds.faces_end();;
}

template <class Gt, class Tds >
typename Triangulation_2<Gt, Tds>::All_vertices_iterator
Triangulation_2<Gt, Tds>::
all_vertices_begin() const
{
  return _tds.vertices_begin();
}

template <class Gt, class Tds >
typename Triangulation_2<Gt, Tds>::All_vertices_iterator
Triangulation_2<Gt, Tds>::
all_vertices_end() const
{
  return _tds.vertices_end();
}

template <class Gt, class Tds >
typename Triangulation_2<Gt, Tds>::All_edges_iterator
Triangulation_2<Gt, Tds>::
all_edges_begin() const
{
  return _tds.edges_begin();
}

template <class Gt, class Tds >
typename Triangulation_2<Gt, Tds>::All_edges_iterator
Triangulation_2<Gt, Tds>::
all_edges_end() const
{
  return _tds.edges_end();
}

template <class Gt, class Tds >
typename Triangulation_2<Gt, Tds>::All_halfedges_iterator
Triangulation_2<Gt, Tds>::
all_halfedges_begin() const
{
  return _tds.halfedges_begin();
}

template <class Gt, class Tds >
typename Triangulation_2<Gt, Tds>::All_halfedges_iterator
Triangulation_2<Gt, Tds>::
all_halfedges_end() const
{
  return _tds.halfedges_end();
}

template <class Gt, class Tds >
inline
typename Triangulation_2<Gt, Tds>::Face_circulator
Triangulation_2<Gt, Tds>::
incident_faces(Vertex_handle v, Face_handle f) const
{
  return _tds.incident_faces(v,f);
}

template <class Gt, class Tds >
inline
typename Triangulation_2<Gt, Tds>::Vertex_circulator
Triangulation_2<Gt, Tds>::
incident_vertices(Vertex_handle v, Face_handle f) const
{
  return _tds.incident_vertices(v,f);
}

template <class Gt, class Tds >
inline
typename Triangulation_2<Gt, Tds>::Edge_circulator
Triangulation_2<Gt, Tds>::
incident_edges(Vertex_handle v, Face_handle f) const
{
  return _tds.incident_edges(v,f);
}

template <class Gt, class Tds >
inline
typename Triangulation_2<Gt, Tds>::size_type
Triangulation_2<Gt, Tds>::
degree(Vertex_handle v) const
{
  return _tds.degree(v);
}

template <class Gt, class Tds >
inline
typename Triangulation_2<Gt, Tds>::Vertex_handle
Triangulation_2<Gt, Tds>::
mirror_vertex(Face_handle f, int i) const
{
  return _tds.mirror_vertex(f,i);
}

template <class Gt, class Tds >
inline
int
Triangulation_2<Gt, Tds>::
mirror_index(Face_handle f, int i) const
{
  return _tds.mirror_index(f,i);
}

template <class Gt, class Tds >
inline
typename Triangulation_2<Gt, Tds>::Edge
Triangulation_2<Gt, Tds>::
mirror_edge(const Edge e) const
{
  return _tds.mirror_edge(e);
}

template <class Gt, class Tds >
typename Triangulation_2<Gt, Tds>::Line_face_circulator
Triangulation_2<Gt, Tds>::
line_walk(const Point& p, const Point& q, Face_handle f) const
{
  CGAL_triangulation_precondition( (dimension() == 2) && ! xy_equal(p,q));
  Line_face_circulator lfc = (f == Face_handle())
                             ? Line_face_circulator(p, q, this)
                             : Line_face_circulator(p, q, f, this);

  // the following lines may be useless :
  //  Line_face_circulator(p,q...) returns either a null circulator
  //  or a pointer to a finite face (to be checked)
  if( (!lfc.is_empty()) && is_infinite( lfc )){
    do { ++lfc ;}
    while (is_infinite(lfc));
  }
  return lfc;
}

template <class Gt, class Tds >
Oriented_side
Triangulation_2<Gt, Tds>::
oriented_side(const Point &p0, const Point &p1,
        const Point &p2, const Point &p) const
{
  // return position of point p with respect to the oriented triangle p0p1p2
  // depends on the orientation of the vertices
  Bounded_side bs=bounded_side(p0,p1,p2,p);
  if (bs == ON_BOUNDARY) return ON_ORIENTED_BOUNDARY;
  Orientation ot = orientation(p0, p1, p2);
  if (bs == ON_BOUNDED_SIDE)
    return (ot == LEFT_TURN) ? ON_POSITIVE_SIDE : ON_NEGATIVE_SIDE;
  // bs == ON_UNBOUNDED_SIDE
  return (ot == LEFT_TURN) ? ON_NEGATIVE_SIDE : ON_POSITIVE_SIDE;
}

template <class Gt, class Tds >
Bounded_side
Triangulation_2<Gt, Tds>::
bounded_side(const Point &p0, const Point &p1,
             const Point &p2, const Point &p) const
{
  // return position of point p with respect to triangle p0p1p2
  CGAL_triangulation_precondition( orientation(p0, p1, p2) != COLLINEAR);
  Orientation o1 = orientation(p0, p1, p),
              o2 = orientation(p1, p2, p),
              o3 = orientation(p2, p0, p);

  if (o1 == COLLINEAR){
    if (o2 == COLLINEAR || o3 == COLLINEAR) return ON_BOUNDARY;
    if (collinear_between(p0, p, p1))        return ON_BOUNDARY;
    return ON_UNBOUNDED_SIDE;
  }

  if (o2 == COLLINEAR){
    if (o3 == COLLINEAR)                     return ON_BOUNDARY;
    if (collinear_between(p1, p, p2))        return ON_BOUNDARY;
    return ON_UNBOUNDED_SIDE;
  }

  if (o3 == COLLINEAR){
    if (collinear_between(p2, p, p0))        return ON_BOUNDARY;
    return ON_UNBOUNDED_SIDE;
  }

  // from here none ot, o1, o2 and o3 are known to be non null
    if (o1 == o2 && o2 == o3) return ON_BOUNDED_SIDE;
    return ON_UNBOUNDED_SIDE;
}

template <class Gt, class Tds >
Oriented_side
Triangulation_2<Gt, Tds>::
oriented_side(Face_handle f, const Point &p) const
{
  CGAL_triangulation_precondition ( dimension()==2);
  return oriented_side(f->vertex(0)->point(),
                       f->vertex(1)->point(),
                       f->vertex(2)->point(),
                       p);
}

template <class Gt, class Tds >
Oriented_side
Triangulation_2<Gt, Tds>::
side_of_oriented_circle(const Point &p0, const Point &p1, const Point &p2,
                        const Point &p, bool perturb) const
{
  //CGAL_triangulation_precondition( orientation(p0, p1, p2) == POSITIVE );
  // no reason for such precondition and it invalidates fast removal in Delaunay

  typename Gt::Side_of_oriented_circle_2 pred = geom_traits().side_of_oriented_circle_2_object();
  Oriented_side os = pred(construct_point(p0), construct_point(p1),
                          construct_point(p2), construct_point(p));

  if ((os != ON_ORIENTED_BOUNDARY) || (! perturb))
    return os;

  // We are now in a degenerate case => we do a symbolic perturbation.

  // We sort the points lexicographically.
  const Point * points[4] = {&p0, &p1, &p2, &p};
  std::sort(points, points+4, Perturbation_order(this) );

  // We successively look whether the leading monomial, then 2nd monomial
  // of the determinant has non null coefficient.
  // 2 iterations are enough if p0p1p2 is positive (cf paper)
  for (int i=3; i>0; --i) {
    if (points[i] == &p)
      return ON_NEGATIVE_SIDE; // since p0 p1 p2 are non collinear
    // and "conceptually" positively oriented
    Orientation o;
    if (points[i] == &p2 && (o = orientation(p0,p1,p)) != COLLINEAR )
      return Oriented_side(o);
    if (points[i] == &p1 && (o = orientation(p0,p,p2)) != COLLINEAR )
      return Oriented_side(o);
    if (points[i] == &p0 && (o = orientation(p,p1,p2)) != COLLINEAR )
      return Oriented_side(o);
  }
  // CGAL_triangulation_assertion(false);
  //no reason for such precondition and it invalidates fast removal in Delaunay
  return ON_NEGATIVE_SIDE;
}

template < class Gt, class Tds >
Oriented_side
Triangulation_2<Gt,Tds>::
side_of_oriented_circle(Face_handle f, const Point & p, bool perturb) const
{
  if ( ! is_infinite(f) ) {
    return this->side_of_oriented_circle(f->vertex(0)->point(),
                                         f->vertex(1)->point(),
                                         f->vertex(2)->point(),p, perturb);
  }

  int i = f->index(infinite_vertex());
  Orientation o = orientation(f->vertex(ccw(i))->point(),
                              f->vertex(cw(i))->point(),
                              p);
  return (o == NEGATIVE) ? ON_NEGATIVE_SIDE :
                           (o == POSITIVE) ? ON_POSITIVE_SIDE :
                                             ON_ORIENTED_BOUNDARY;
}

template <class Gt, class Tds >
bool
Triangulation_2<Gt, Tds>::
collinear_between(const Point& p, const Point& q, const Point& r) const
{
  // return true if point q is strictly between p and r
  // p,q and r are supposed to be collinear points
  Comparison_result c_pr = compare_x(p, r);
  Comparison_result c_pq;
  Comparison_result c_qr;
  if(c_pr == EQUAL) {
    //c_pr = compare_y(p, r);
    c_pq = compare_y(p, q);
    c_qr = compare_y(q, r);
  } else {
    c_pq = compare_x(p, q);
    c_qr = compare_x(q, r);
  }
  return ( (c_pq == SMALLER) && (c_qr == SMALLER) ) ||
         ( (c_pq == LARGER) && (c_qr == LARGER) );

}

template <class Gt, class Tds >
inline
Comparison_result
Triangulation_2<Gt, Tds>::
compare_x(const Point& p, const Point& q) const
{
  return geom_traits().compare_x_2_object()(construct_point(p),
                                            construct_point(q));
}

template <class Gt, class Tds >
inline
Comparison_result
Triangulation_2<Gt, Tds>::
compare_xy(const Point& p, const Point& q) const
{
  Comparison_result res = geom_traits().compare_x_2_object()(construct_point(p),
                                                             construct_point(q));
  if(res == EQUAL){
    return geom_traits().compare_y_2_object()(construct_point(p),
                                              construct_point(q));
  }
  return res;
}

template <class Gt, class Tds >
inline
Comparison_result
Triangulation_2<Gt, Tds>::
compare_y(const Point& p, const Point& q) const
{
  return geom_traits().compare_y_2_object()(construct_point(p),
                                            construct_point(q));
}

template <class Gt, class Tds >
inline
bool
Triangulation_2<Gt, Tds>::
xy_equal(const Point& p, const Point& q) const
{
  return compare_x(p,q)== EQUAL && compare_y(p,q)== EQUAL ;
}

template <class Gt, class Tds >
inline
Orientation
Triangulation_2<Gt, Tds>::
orientation(const Point& p, const Point& q,const Point& r ) const
{
  return geom_traits().orientation_2_object()(construct_point(p),
                                              construct_point(q),
                                              construct_point(r));
}

template<class Gt, class Tds>
inline
typename Triangulation_2<Gt,Tds>::Point_2
Triangulation_2<Gt,Tds>::
circumcenter(const Point& p0, const Point& p1, const Point& p2) const
{
  return
    geom_traits().construct_circumcenter_2_object()(construct_point(p0),
                                                    construct_point(p1),
                                                    construct_point(p2));
}

template <class Gt, class Tds >
typename Triangulation_2<Gt, Tds>::Point_2
Triangulation_2<Gt, Tds>::
circumcenter(Face_handle f) const
{
  CGAL_triangulation_precondition (dimension()==2);
  return circumcenter((f->vertex(0))->point(),
                      (f->vertex(1))->point(),
                      (f->vertex(2))->point());
}

template <class Gt, class Tds >
void
Triangulation_2<Gt, Tds>::
show_all() const
{
  std::cerr<< "PRINT THE COMPLETE TRIANGULATION :"<<std::endl;
  std::cerr << std::endl<<"====> "<< this;
  std::cerr <<  " dimension " << dimension() << std::endl;
  std::cerr << "nb of vertices " << number_of_vertices() << std::endl;

  if (dimension() < 1) return;
  if(dimension() == 1) {
    std::cerr<<" all edges "<<std::endl;
    All_edges_iterator aeit;
    for(aeit = all_edges_begin(); aeit != all_edges_end(); aeit++){
      show_face(aeit->first);
    }
    return;
  }

  std::cerr<<" finite faces "<<std::endl;
  Finite_faces_iterator fi;
  for(fi = finite_faces_begin(); fi != finite_faces_end(); fi++) {
    show_face(fi);
  }

  std::cerr <<" infinite faces "<<std::endl;
  All_faces_iterator afi;
  for(afi = all_faces_begin(); afi != all_faces_end(); afi++) {
    if(is_infinite(afi)) show_face(afi);
  }

  if (number_of_vertices()>1) {
    std::cerr << "print vertices of the regular triangulation"
        <<std::endl;
    All_vertices_iterator vi;
    for( vi = all_vertices_begin(); vi != all_vertices_end(); vi++){
      show_vertex(vi);
      std::cerr << "  / associated face: "
       << (void*)(&(*(vi->face())))<< std::endl;;
      }
      std::cerr<<std::endl;
  }
  return;
}

template <class Gt, class Tds >
void
Triangulation_2<Gt, Tds>::
show_vertex(Vertex_handle vh) const
{
  if(is_infinite(vh)) std::cerr << "inf \t";
  else std::cerr << vh->point() << "\t";
  return;
}

template <class Gt, class Tds >
void
Triangulation_2<Gt, Tds>::
show_face(Face_handle fh) const
{
  std::cerr << "face : "<<(void*)&(*fh)<<" => "<<std::endl;
  int i = fh->dimension();
  switch(i){
  case 0:
    std::cerr <<"point :" ; show_vertex(fh->vertex(0));
    std::cerr <<" / neighbor " << &(*(fh->neighbor(0)));
    std::cerr <<"[" ; show_vertex(fh->neighbor(0)->vertex(0));
    std::cerr <<"]"  << std::endl;
    break;
  case 1:
     std::cerr <<"point :" ; show_vertex(fh->vertex(0));
     std::cerr <<" / neighbor " << &(*(fh->neighbor(0)));
     std::cerr <<"[" ; show_vertex(fh->neighbor(0)->vertex(0));
     std::cerr <<"/" ; show_vertex(fh->neighbor(0)->vertex(1));
     std::cerr <<"]" <<std::endl;

     std::cerr <<"point :" ; show_vertex(fh->vertex(1));
     std::cerr <<" / neighbor " << &(*(fh->neighbor(1)));
     std::cerr <<"[" ; show_vertex(fh->neighbor(1)->vertex(0));
     std::cerr <<"/" ; show_vertex(fh->neighbor(1)->vertex(1));
     std::cerr <<"]" <<std::endl;
     break;
  case 2:
    std::cerr <<"point :" ; show_vertex(fh->vertex(0));
    std::cerr <<" / neighbor " << &(*(fh->neighbor(0)));
    std::cerr <<"[" ; show_vertex(fh->neighbor(0)->vertex(0));
    std::cerr <<"/" ; show_vertex(fh->neighbor(0)->vertex(1));
    std::cerr <<"/" ; show_vertex(fh->neighbor(0)->vertex(2));
    std::cerr <<"]" <<std::endl;

    std::cerr <<"point :" ; show_vertex(fh->vertex(1));
    std::cerr <<" / neighbor " << &(*(fh->neighbor(1)));
    std::cerr <<"[" ; show_vertex(fh->neighbor(1)->vertex(0));
    std::cerr <<"/" ; show_vertex(fh->neighbor(1)->vertex(1));
    std::cerr <<"/" ; show_vertex(fh->neighbor(1)->vertex(2));
    std::cerr <<"]" <<std::endl;

    std::cerr <<"point :" ; show_vertex(fh->vertex(2));
    std::cerr <<" / neighbor " << &(*(fh->neighbor(2)));
    std::cerr <<"[" ; show_vertex(fh->neighbor(2)->vertex(0));
    std::cerr <<"/" ; show_vertex(fh->neighbor(2)->vertex(1));
    std::cerr <<"/" ; show_vertex(fh->neighbor(2)->vertex(2));
    std::cerr <<"]" <<std::endl;
    break;
  }
  return;
}

template <class Gt, class Tds >
void
Triangulation_2<Gt, Tds>::
file_output(std::ostream& os) const
{
  _tds.file_output(os, infinite_vertex(), true);
}

template <class Gt, class Tds >
typename Triangulation_2<Gt, Tds>::Vertex_handle
Triangulation_2<Gt, Tds>::
file_input(std::istream& is)
{
  clear();
  Vertex_handle v= _tds.file_input(is, true);
  set_infinite_vertex(v);
  return v;
}

template <class Gt, class Tds >
std::ostream&
operator<<(std::ostream& os, const Triangulation_2<Gt, Tds> &tr)
{
  tr.file_output(os);
  return os ;
}

template < class Gt, class Tds >
std::istream&
operator>>(std::istream& is, Triangulation_2<Gt, Tds> &tr)
{
  tr.file_input(is);
  CGAL_triangulation_assertion(tr.is_valid());
  return is;
}

} //namespace CGAL

#endif //CGAL_TRIANGULATION_2_H
