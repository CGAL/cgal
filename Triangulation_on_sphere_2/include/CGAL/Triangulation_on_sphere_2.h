// Copyright (c) 1997, 2012, 2019 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mariette Yvinec,
//                 Claudia Werner
//                 Mael Rouxel-Labbé
//
#ifndef CGAL_TRIANGULATION_ON_SPHERE_2_H
#define CGAL_TRIANGULATION_ON_SPHERE_2_H

#include <CGAL/license/Triangulation_on_sphere_2.h>

#include <CGAL/assertions.h>
#include <CGAL/Triangulation_utils_2.h>
#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/Triangulation_on_sphere_vertex_base_2.h>
#include <CGAL/Triangulation_on_sphere_face_base_2.h>
#include <CGAL/Triangulation_on_sphere_2/IO/OFF.h>
#include <CGAL/Triangulation_on_sphere_2/internal/arc_on_sphere_2_subsampling.h> // included for convenience

#include <CGAL/iterator.h>
#include <CGAL/Iterator_project.h>
#include <CGAL/function_objects.h>
#include <CGAL/boost/iterator/transform_iterator.hpp>

#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_smallint.hpp>
#include <boost/random/variate_generator.hpp>

#include <algorithm>
#include <iostream>
#include <fstream>
#include <list>
#include <utility>

namespace CGAL {

template < class Gt, class Tds >
class Triangulation_on_sphere_2;

template < class Gt, class Tds >
std::istream& operator>>(std::istream& is, Triangulation_on_sphere_2<Gt, Tds>& tr);

template < class Gt, class Tds >
std::ostream& operator<<(std::ostream& os, const Triangulation_on_sphere_2<Gt, Tds>& tr);

// This class just provides some basic methods and cannot be used independently.
// No insertion or removal is implemented.
template <class Gt,
          class Tds = Triangulation_data_structure_2 <
                        Triangulation_on_sphere_vertex_base_2<Gt>,
                        Triangulation_on_sphere_face_base_2<Gt> > >
class Triangulation_on_sphere_2
  : public Triangulation_cw_ccw_2
{
  typedef Triangulation_on_sphere_2<Gt, Tds>      Self;

public:
  typedef Tds                                     Triangulation_data_structure;
  typedef Gt                                      Geom_traits;

  typedef typename Geom_traits::FT                FT;
  typedef typename Geom_traits::Point_3           Point_3;
  typedef typename Geom_traits::Point_on_sphere_2 Point;
  typedef typename Geom_traits::Segment_3         Segment_3;
  typedef typename Geom_traits::Triangle_3        Triangle_3;
  typedef typename Geom_traits::Arc_on_sphere_2   Arc_on_sphere_2;

  typedef typename Geom_traits::Construct_point_3 Construct_point_3;

  typedef typename Tds::size_type                 size_type;
  typedef typename Tds::difference_type           difference_type;
  typedef typename Tds::Vertex                    Vertex;
  typedef typename Tds::Vertex_handle             Vertex_handle;
  typedef typename Tds::Edge                      Edge;
  typedef typename Tds::Face                      Face;
  typedef typename Tds::Face_handle               Face_handle;

  enum Locate_type {VERTEX = 0,
                    EDGE, // 1
                    FACE, // 2
                    OUTSIDE_CONVEX_HULL, // 3
                    OUTSIDE_AFFINE_HULL, // 4
                    CONTOUR, // 5
                    NOT_ON_SPHERE, // 6
                    TOO_CLOSE}; // 7

protected:
  Gt _gt;
  Tds _tds;

public:
  //-----------------------Creation-----------------------------------------------------------------
  Triangulation_on_sphere_2(const Geom_traits& gt = Geom_traits()) : _gt(gt) { }
  Triangulation_on_sphere_2(const Point_3& center, const FT radius) : _gt(center, radius) { }
  Triangulation_on_sphere_2(const Triangulation_on_sphere_2<Gt, Tds>& tr);
  void clear();

public:
  // these functions clear the triangulation
  void set_center(const Point_3& c)
  {
    clear();
    _gt.set_center(c);
  }

  void set_radius(const FT radius)
  {
    clear();
    _gt.set_radius(radius);
  }

  void set_center_and_radius(const Point_3& c, const FT radius)
  {
    clear();
    _gt.set_center(c);
    _gt.set_radius(radius);
  }

  // Assignment
  void swap(Triangulation_on_sphere_2& tr);
  Triangulation_on_sphere_2& operator=(Triangulation_on_sphere_2 tr); // intentional copy

  // Destructor
  ~Triangulation_on_sphere_2() = default;

public:
  // Members
  const Geom_traits& geom_traits() const { return _gt; }
  const Tds& tds() const { return _tds; }
  Tds& tds() { return _tds; }

  int dimension() const { return _tds.dimension(); }

  size_type number_of_vertices() const { return _tds.number_of_vertices(); }
  size_type number_of_edges() const { return _tds.number_of_edges(); }
  size_type number_of_faces() const { return _tds.number_of_faces(); } // total number of faces (solid + ghost)
  size_type number_of_ghost_faces() const
  {
    std::size_t nb = 0;
    for(All_faces_iterator it=all_faces_begin(); it!=all_faces_end(); ++it)
      if(is_ghost(it))
        ++nb;

    return nb;
  }

  size_type number_of_solid_faces() const { return number_of_faces() - number_of_ghost_faces(); }

  // TDS access
  bool is_edge(Vertex_handle va, Vertex_handle vb) const { return _tds.is_edge(va, vb); }
  bool is_edge(Vertex_handle va, Vertex_handle vb, Face_handle& fr, int& i) const
  {
    return _tds.is_edge(va, vb, fr, i);
  }

  bool is_face(Vertex_handle v1, Vertex_handle v2, Vertex_handle v3) const
  {
    return _tds.is_face(v1, v2, v3);
  }

  bool is_face(Vertex_handle v1, Vertex_handle v2, Vertex_handle v3, Face_handle& fr) const
  {
    return _tds.is_face(v1, v2, v3, fr);
  }

  size_type degree(Vertex_handle v) const { return _tds.degree(v); }

  Vertex_handle mirror_vertex(Face_handle f, int i) const { return _tds.mirror_vertex(f, i); }
  int mirror_index(Face_handle v, int i) const { return _tds.mirror_index(v, i); }
  Edge mirror_edge(const Edge& e) const { return _tds.mirror_edge(e); }

  bool is_ghost(const Face_handle f) const { return f->is_ghost(); }
  bool is_ghost(const Face_handle f, const int i) const
  {
    CGAL_precondition(dimension() == 2);
    const bool g1 = is_ghost(f);
    const bool g2 = is_ghost(f->neighbor(i));
    return g1 && g2;
  }

  bool is_ghost(const Edge& e) const
  {
    return is_ghost(e.first, e.second);
  }

  bool is_contour(const Face_handle f, const int i) const
  {
    CGAL_precondition(dimension() == 2);
    const bool edge1 = is_ghost(f);
    const bool edge2 = is_ghost(f->neighbor(i));
    return edge1 != edge2; // xor
  }

  bool is_contour(const Edge& e) const
  {
    return is_contour(e.first, e.second);
  }

  template<class T>
  bool is_infinite(const T&, int = 0) const
  {
    return false;
  }

  //-----------------------TRAVERSING : ITERATORS AND CIRCULATORS-----------------------------------

  typedef typename Tds::Vertex_iterator           Vertices_iterator;
  typedef typename Tds::Vertex_iterator           All_vertices_iterator;
  typedef typename Tds::Edge_iterator             All_edges_iterator;
  typedef typename Tds::Face_iterator             All_faces_iterator;

  typedef typename Tds::Vertex_iterator           Finite_vertices_iterator;
  typedef typename Tds::Edge_iterator             Finite_edges_iterator;
  typedef typename Tds::Face_iterator             Finite_faces_iterator;

  typedef typename Tds::Vertex_circulator         Vertex_circulator;
  typedef typename Tds::Edge_circulator           Edge_circulator;
  typedef typename Tds::Face_circulator           Face_circulator;

  // This class is used to generate the solid iterators
  class Ghost_tester
  {
    const Triangulation_on_sphere_2& tr;

  public:
    Ghost_tester(const Triangulation_on_sphere_2& tr) : tr(tr) { }

    bool operator()(const All_faces_iterator fit) const { return tr.is_ghost(fit); }
    bool operator()(const All_edges_iterator eit) const { return tr.is_ghost(*eit); }
  };

  class Non_contour_tester
  {
    const Triangulation_on_sphere_2& tr;

  public:
    Non_contour_tester(const Triangulation_on_sphere_2& tr) : tr(tr) { }

    bool operator()(const All_edges_iterator eit) const
    {
      return !tr.is_contour(*eit);
    }
  };

  class Solid_faces_iterator
    : public Filter_iterator<All_faces_iterator, Ghost_tester>
  {
    typedef Filter_iterator<All_faces_iterator, Ghost_tester> Base;
    typedef Solid_faces_iterator                              Self;

  public:
    Solid_faces_iterator() : Base() { }
    Solid_faces_iterator(const Base& b) : Base(b) { }
    Self& operator++() { Base::operator++(); return *this; }
    Self& operator--() { Base::operator--(); return *this; }
    Self operator++(int) { Self tmp(*this); ++(*this); return tmp; }
    Self operator--(int) { Self tmp(*this); --(*this); return tmp; }
    operator const Face_handle&() const { return Base::base(); }
  };

  typedef Iterator_range<Prevent_deref<Vertices_iterator> >          Vertex_handles;
  typedef Iterator_range<All_edges_iterator>                         All_edges;
  typedef Iterator_range<Prevent_deref<All_faces_iterator> >         All_face_handles;

  // one solid and one ghost face incident to this edge
  typedef Filter_iterator<All_edges_iterator, Non_contour_tester>    Contour_edges_iterator;

  // solid edges: both incident faces are solid
  typedef Filter_iterator<All_edges_iterator, Ghost_tester>          Solid_edges_iterator;
  typedef Iterator_range<Solid_edges_iterator>                       Solid_edges;
  typedef Iterator_range<Prevent_deref<Solid_faces_iterator, const Face_handle&>> Solid_face_handles;

  typedef Project_point<Vertex>                                      Pt_proj;
  typedef boost::transform_iterator<Pt_proj, Vertices_iterator>      Point_iterator;
  typedef Iterator_range<Point_iterator>                             Points;

  // -- faces
  All_faces_iterator all_faces_begin() const { return _tds.faces_begin(); }
  All_faces_iterator all_faces_end() const { return _tds.faces_end(); }
  Finite_faces_iterator finite_faces_begin() const { return _tds.faces_begin(); }
  Finite_faces_iterator finite_faces_end() const { return _tds.faces_end(); }

  All_face_handles all_face_handles() const
  {
    return make_prevent_deref_range(all_faces_begin(), all_faces_end());
  }

  // -- edges
  Contour_edges_iterator contour_edges_begin() const
  {
    if(dimension() < 1)
      return contour_edges_end();

    return CGAL::filter_iterator(all_edges_end(), Non_contour_tester(*this), all_edges_begin());
  }

  Contour_edges_iterator contour_edges_end() const
  {
    return CGAL::filter_iterator(all_edges_end(), Non_contour_tester(*this));
  }

  Solid_faces_iterator solid_faces_begin() const
  {
    if(dimension() < 2)
      return solid_faces_end();

    return CGAL::filter_iterator(all_faces_end(), Ghost_tester(*this), all_faces_begin());
  }

  Solid_faces_iterator solid_faces_end() const
  {
    return CGAL::filter_iterator(all_faces_end(), Ghost_tester(*this));
  }

  Solid_face_handles solid_faces() const
  {
    return { solid_faces_begin(), solid_faces_end() };
  }

  Solid_edges_iterator solid_edges_begin() const
  {
    if(dimension() < 1)
      return solid_edges_end();

    return CGAL::filter_iterator(all_edges_end(), Ghost_tester(*this), all_edges_begin());
  }

  Solid_edges_iterator solid_edges_end() const
  {
    return CGAL::filter_iterator(all_edges_end(), Ghost_tester(*this));
  }

  Solid_edges solid_edges() const
  {
    return CGAL::make_range(solid_edges_begin(), solid_edges_end());
  }

  All_edges_iterator all_edges_begin() const { return _tds.edges_begin(); }
  All_edges_iterator all_edges_end() const { return _tds.edges_end(); }
  Finite_edges_iterator finite_edges_begin() const { return _tds.edges_begin(); }
  Finite_edges_iterator finite_edges_end() const { return _tds.edges_end(); }

  All_edges all_edges() const { return _tds.edges(); }

  // -- vertices
  Vertices_iterator vertices_begin() const { return _tds.vertices_begin(); }
  Vertices_iterator vertices_end() const { return _tds.vertices_end(); }
  All_vertices_iterator all_vertices_begin() const { return _tds.vertices_begin(); }
  All_vertices_iterator all_vertices_end() const { return _tds.vertices_end(); }
  Finite_vertices_iterator finite_vertices_begin() const { return _tds.vertices_begin(); }
  Finite_vertices_iterator finite_vertices_end() const { return _tds.vertices_end(); }

  Vertex_handles vertex_handles() const
  {
    return make_prevent_deref_range(vertices_begin(), vertices_end());
  }

  Point_iterator points_begin() const
  {
    return Point_iterator(vertices_begin());
  }

  Point_iterator points_end() const
  {
    return Point_iterator(vertices_end());
  }

  Points points() const
  {
    return Points(points_begin(), points_end());
  }

  Vertex_circulator incident_vertices(Vertex_handle v, Face_handle f = Face_handle()) const
  {
    return _tds.incident_vertices(v, f);
  }

  Edge_circulator incident_edges(Vertex_handle v, Face_handle f = Face_handle()) const
  {
    return _tds.incident_edges(v, f);
  }

  Face_circulator incident_faces(Vertex_handle v, Face_handle f = Face_handle()) const
  {
    return _tds.incident_faces(v, f);
  }

  //-----------------------GEOMETRY-----------------------------------------------------------------
  const Point& point(const Vertex_handle v) const { return v->point(); }
  const Point& point(const Face_handle f, const int i) const { return point(f->vertex(i)); }

  decltype(auto) construct_point(const Point& p) const
  {
    return geom_traits().construct_point_3_object()(p);
  }

  Segment_3 segment(const Face_handle f, int i) const
  {
    return geom_traits().construct_segment_3_object()(construct_point(point(f, ccw(i))),
                                                      construct_point(point(f, cw(i))));
  }

  Segment_3 segment(const Edge& e) const
  {
    return segment(e.first, e.second);
  }

  Segment_3 segment(const Edge_circulator& ec) const { return segment(*ec); }
  Segment_3 segment(const All_edges_iterator ei) const { return segment(*ei); }

  Arc_on_sphere_2 segment_on_sphere(const Face_handle f, int i) const
  {
    CGAL_precondition(!is_ghost(f, i));

    return geom_traits().construct_arc_on_sphere_2_object()(point(f, ccw(i)), point(f, cw(i)));
  }

  Arc_on_sphere_2 segment_on_sphere(const Edge& e) const
  {
    return segment_on_sphere(e.first, e.second);
  }

  Triangle_3 triangle(const Face_handle f) const
  {
    return geom_traits().construct_triangle_3_object()(construct_point(point(f, 0)),
                                                       construct_point(point(f, 1)),
                                                       construct_point(point(f, 2)));
  }

  //-----------------------PREDICATES---------------------------------------------------------------
  bool are_equal(const Point& p, const Point& q) const
  {
    return geom_traits().equal_on_sphere_2_object()(p, q);
  }

  Comparison_result compare(const Point& p, const Point& q) const
  {
    return geom_traits().compare_on_sphere_2_object()(p, q);
  }

  bool collinear_between(const Point& p, const Point& q, const Point& r) const
  {
    return geom_traits().collinear_are_strictly_ordered_on_great_circle_2_object()(p, q, r);
  }

  Orientation orientation(const Point& p, const Point& q, const Point& r, const Point& s) const
  {
    return geom_traits().side_of_oriented_circle_on_sphere_2_object()(p, q, r, s);
  }

  Orientation orientation(const Face_handle f, const Point& p) const
  {
    return orientation(point(f, 0), point(f, 1), point(f, 2), p);
  }

  Orientation orientation_on_sphere(const Point& p, const Point& q, const Point& r) const
  {
    return geom_traits().orientation_on_sphere_2_object()(p, q, r);
  }

  Orientation orientation_on_sphere(const Face_handle f) const
  {
    return orientation_on_sphere(point(f, 0), point(f, 1), point(f, 2));
  }

  //-----------------------LOCATION-----------------------------------------------------------------
  bool are_points_too_close(const Point& p, const Point& q, Locate_type& lt) const;
  void test_distance(const Point& p, Face_handle& f, Locate_type& lt, int& li) const;

  Face_handle march_locate_1D(const Point& t, Locate_type& lt, int& li) const ;
  Face_handle march_locate_2D(Face_handle f, const Point& t, Locate_type& lt, int& li) const;
  Face_handle locate(const Point& p, Locate_type& lt, int& li, Face_handle start) const;
  Face_handle locate(const Point& p, const Face_handle start) const;
  Face_handle locate_edge(const Point& p, Locate_type& lt, int& li, bool on_diametral_plane) const;

  //-----------------------DEBUG--------------------------------------------------------------------
  void check_neighboring_faces() const;
  bool is_plane() const;
  bool is_valid_vertex(Vertex_handle fh, bool verbose = false, int level = 0) const;
  bool is_valid(bool verbose = false, int level = 0) const;

  void show_all() const;
  void show_vertex(Vertex_handle vh) const;
  void show_face(Face_handle fh) const;

  // IN/OUT
  Vertex_handle file_input(std::istream& is);
  void file_output(std::ostream& os) const;
};

// copy constructor duplicates vertices and faces
template <typename Gt, typename Tds>
Triangulation_on_sphere_2<Gt, Tds>::
Triangulation_on_sphere_2(const Triangulation_on_sphere_2& other)
  : _gt(other.geom_traits())
{
  _tds.copy_tds(other.tds());
}

template <typename Gt, typename Tds>
void
Triangulation_on_sphere_2<Gt, Tds>::
clear()
{
  _tds.clear();
}

// Assignment
template <typename Gt, typename Tds>
void
Triangulation_on_sphere_2<Gt, Tds>::
swap(Triangulation_on_sphere_2& tr)
{
  using std::swap;
  _tds.swap(tr._tds);
  swap(_gt, tr._gt);
}

template <typename Gt, typename Tds>
Triangulation_on_sphere_2<Gt, Tds>&
Triangulation_on_sphere_2<Gt, Tds>::
operator=(Triangulation_on_sphere_2 tr) // intentional copy
{
  swap(tr);
  return *this;
}

//----------------------------------------POINT LOCATION---------------------------------------//

// tests whether the two points p and q are too close according to the lemma about hidden vertices.
template<class Gt, class Tds>
bool
Triangulation_on_sphere_2<Gt, Tds>::
are_points_too_close(const Point& p, const Point& q, Locate_type& lt) const
{
  if(geom_traits().are_points_too_close(p, q))
  {
    lt = TOO_CLOSE;
    return true;
  }

  return false;
}

/*
 calls are_points_too_close() for possible conflicts.
 If the point p is too close to an existing vertex, this vertex is returned.
 */
template <typename Gt, typename Tds>
void
Triangulation_on_sphere_2<Gt, Tds>::
test_distance(const Point& p,
              Face_handle& f,
              Locate_type& lt,
              int& li) const
{
  CGAL_precondition(dimension() >= -1);

  switch(dimension())
  {
    case -1: // 1 vertex
      CGAL_FALLTHROUGH;
    case 0: // 2 vertices
      CGAL_FALLTHROUGH;
    case 1: // 3+ coplanar vertices
    {
      for(Vertices_iterator vi=vertices_begin(); vi!=vertices_end(); ++vi)
      {
        if(are_points_too_close(point(vi), p, lt))
        {
          f = vi->face();
          li = f->index(vi);
          return;
        }
      }
      break;
    }
    case 2:
    {
      // @fixme this seems wrong, the face in conflict might not have the closest vertex

      const Vertex_handle v0 = f->vertex(0);
      if(are_points_too_close(point(v0), p, lt))
      {
        li = 0;
        return;
      }

      const Vertex_handle v1 = f->vertex(1);
      if(are_points_too_close(point(v1), p, lt))
      {
        li = 1;
        return;
      }

      const Vertex_handle v2 = f->vertex(2);
      if(are_points_too_close(point(v2), p, lt))
      {
        li = 2;
        return;
      }
      break;
    }
  }
}

/*
 * Location for degenerated cases: locates the conflicting edge in a 1-dimensional triangulation.
 * This method is used when the new point is coplanar with the existing vertices.
 * The Boolean 'on_diametral_plane' indicates whether the points are also coplanar with the
 * center of the sphere (true) or not (false).
 */
template <class Gt, class Tds>
typename Triangulation_on_sphere_2<Gt, Tds>::Face_handle
Triangulation_on_sphere_2<Gt, Tds>::
locate_edge(const Point& p,
            Locate_type& lt,
            int& li,
            const bool on_diametral_plane) const
{
  CGAL_precondition(dimension() == 1);

  Face_handle loc;
  if(on_diametral_plane)
  {
    All_edges_iterator eit;
    for(eit=all_edges_begin(); eit!=all_edges_end(); ++eit)
    {
      const Face_handle f = eit->first;
      const Face_handle fn = f->neighbor(0);
      const Point& q = point(fn, 1);

      // Check if *eit is an a "ghost" edge, that is if in the diametral plane, the center
      // of the sphere is to the right of (eit->source(); eit->target()).
      // Note that this doesn't check if p is actually "in conflict" with that edge,
      // so that edge is just kept in memory, and if no solid edge is found, then it has to be this edge
      if(collinear_between(point(f, 0), point(f, 1), q))
      {
        // The new point is on the same 3D plane as the existing vertices, but not on an existing edge
        loc = eit->first;
      }
      else
      {
        // not ghost, check the cone
        if(collinear_between(point(eit->first, 0), point(eit->first, 1), p))
        {
          loc = eit->first;
          lt = EDGE;
          li = 2;
          test_distance(p, loc, lt, li);
          return loc;
        }
      }
    }

    // Couldn't find a solid edge that is in conflict with p
    lt = OUTSIDE_CONVEX_HULL;
    li = 4;

    test_distance(p, loc, lt, li);
    return loc;
  }
  else // not coplanar with the center of the sphere
  {
    for(All_edges_iterator eit = all_edges_begin(); eit!=all_edges_end(); ++eit)
    {
      // The plane going through the center of the sphere, v1, and v2 cuts the non-diametral circle
      // in two, with the negative side of the plane corresponding to the shorter part, so
      // if p is coplanar with all the points on that circle and on the negative side, it splits
      // the edge v1-v2
      if(orientation_on_sphere(point(eit->first, 0), point(eit->first, 1), p) == ON_NEGATIVE_SIDE)
      {
        loc = eit->first;
        lt = EDGE;
        li = 2;
        test_distance(p, loc, lt, li);
        return loc;
      }
    }

    CGAL_assertion(false);
    return loc;
  }
}

template <typename Gt, typename Tds>
typename Triangulation_on_sphere_2<Gt, Tds>::Face_handle
Triangulation_on_sphere_2<Gt, Tds>::
march_locate_1D(const Point& p, Locate_type& lt, int& li) const
{
  CGAL_assertion(dimension() == 1);

  // Check if p is coplanar with the existing points first three points of the triangulation
  Face_handle f = all_edges_begin()->first;
  const Vertex_handle v1 = f->vertex(0);
  const Vertex_handle v2 = f->vertex(1);
  const Vertex_handle v3 = f->neighbor(0)->vertex(1);

  const Orientation orient = orientation(point(v1), point(v2), point(v3), p);
  if(orient != ON_ORIENTED_BOUNDARY)
  {
    lt = OUTSIDE_AFFINE_HULL;
    li = 4;
    test_distance(p, f, lt, li);
    return f;
  }

  // From then on, p is coplanar with all the triangulation's points

  // Check if p is coradial with one existing point
  Vertices_iterator vi;
  for(vi=vertices_begin(); vi!=vertices_end(); ++vi) // @todo turns insertion into O(n)
  {
    if(are_equal(point(vi), p))
    {
      lt = VERTEX;
      f = vi->face();
      li = f->index(vi); // could be simply '1'
      return f;
    }
  }

  // v1, v2, and v3 are coplanar so this is just checking if they are coplanar with the center of the sphere
  const Orientation pqr = orientation_on_sphere(point(v1), point(v2), point(v3));
  if(pqr == ON_ORIENTED_BOUNDARY)
    return locate_edge(p, lt, li, true /*on_diametral_plane*/);
  else
    return locate_edge(p, lt, li, false /*on_diametral_plane*/);
}

template <typename Gt, typename Tds>
typename Triangulation_on_sphere_2<Gt, Tds>::Face_handle
Triangulation_on_sphere_2<Gt, Tds>::
march_locate_2D(Face_handle f,
                const Point& t,
                Locate_type& lt,
                int& li) const
{
  CGAL_precondition(dimension() == 2);
  CGAL_precondition(!is_ghost(f));

  boost::rand48 rng;
  boost::uniform_smallint<> two(0, 1);
  boost::variate_generator<boost::rand48&, boost::uniform_smallint<> > coin(rng, two);

  Face_handle prev = Face_handle();
  bool first = true;

  for(;;)
  {
    if(is_ghost(f))
    {
      if(orientation(f, t) == ON_POSITIVE_SIDE) // conflict with the corresponding face
      {
        lt = OUTSIDE_CONVEX_HULL;
        li = 4;
        test_distance(t, f, lt, li);
        return f;
      }
      else // one of the neighbors has to be in conflict with 't'
      {
        Face_handle next = Face_handle();
        for(int i=0; i<=2; ++i)
        {
          next = f->neighbor(i);
          if(orientation(next, t) == ON_POSITIVE_SIDE)
          {
            lt = CONTOUR;
            li = 4;
            // @fixme same issue as the test_distance on the face in conflict in the generic case
            test_distance(t, next, lt, li);
            return next;
          }
        }

        CGAL_precondition(false);
      }
    }

    const Point& p0 = point(f, 0);
    const Point& p1 = point(f, 1);
    const Point& p2 = point(f, 2);

    // Instead of testing f's edges in a random order we do the following
    // until we find a neighbor to go further:
    // As we come from 'prev', we do not have to check the edge leading to 'prev'
    // Now, we flip a coin in order to decide if we start checking the
    // edge before or the edge after the edge leading to 'prev'
    // We do loop unrolling in order to find out if this is faster.
    // In the very beginning we do not have a prev, but for the first step
    // we do not need randomness
    int left_first = coin()%2;
    Orientation o0, o1, o2;

    /************************FIRST*************************/
    if(first)
    {
      prev = f;
      first = false;
      o0 = orientation_on_sphere(p0, p1, t); // classic march locate
      if(o0 == NEGATIVE)
      {
        f = f->neighbor(2);
        continue;
      }

      o1 = orientation_on_sphere(p1, p2, t);
      if(o1 == NEGATIVE)
      {
        f = f->neighbor(0);
        continue;
      }

      o2 = orientation_on_sphere(p2, p0, t);
      if(o2 == NEGATIVE)
      {
        f = f->neighbor(1);
        continue;
      }
    }
    //****LEFT****//
    else if(left_first)
    {
      if(f->neighbor(0) == prev)
      {
        prev = f;
        o0 = orientation_on_sphere(p0, p1, t);
        if(o0 == NEGATIVE)
        {
          f = f->neighbor(2);
          continue;
        }

        o2 = orientation_on_sphere(p2, p0, t);
        if(o2 == NEGATIVE)
        {
          f = f->neighbor(1);
          continue;
        }
        o1 = orientation_on_sphere(p1, p2, t);
      }
      else if(f->neighbor(1) == prev)
      {
        prev = f;
        o1 = orientation_on_sphere(p1, p2, t);
        if(o1 == NEGATIVE)
        {
          f = f->neighbor(0);
          continue;
        }

        o0 = orientation_on_sphere(p0, p1, t);
        if(o0 == NEGATIVE)
        {
          f = f->neighbor(2);
          continue;
        }

        o2 = orientation_on_sphere(p2, p0, t);
      }
      else
      {
        prev = f;
        o2 = orientation_on_sphere(p2, p0, t);
        if(o2 == NEGATIVE)
        {
          f = f->neighbor(1);
          continue;
        }

        o1 = orientation_on_sphere(p1, p2, t);
        if(o1 == NEGATIVE)
        {
          f = f->neighbor(0);
          continue;
        }

        o0 = orientation_on_sphere(p0, p1, t);
      }
    }
    else
    {
      if(f->neighbor(0) == prev)
      {
        prev = f;
        o2 = orientation_on_sphere(p2, p0, t);
        if(o2 == NEGATIVE)
        {
          f = f->neighbor(1);
          continue;
        }

        o0 = orientation_on_sphere(p0, p1, t);
        if(o0 == NEGATIVE)
        {
          f = f->neighbor(2);
          continue;
        }
        o1 = orientation_on_sphere(p1, p2, t);
      }
      else if(f->neighbor(1) == prev)
      {
        prev = f;
        o0 = orientation_on_sphere(p0, p1, t);
        if(o0 == NEGATIVE)
        {
          f = f->neighbor(2);
          continue;
        }

        o1 = orientation_on_sphere(p1, p2, t);
        if(o1 == NEGATIVE)
        {
          f = f->neighbor(0);
          continue;
        }

        o2 = orientation_on_sphere(p2, p0, t);
      }
      else
      {
        prev = f;
        o1 = orientation_on_sphere(p1, p2, t);
        if(o1 == NEGATIVE)
        {
          f = f->neighbor(0);
          continue;
        }

        o2 = orientation_on_sphere(p2, p0, t);
        if(o2 == NEGATIVE)
        {
          f = f->neighbor(1);
          continue;
        }

        o0 = orientation_on_sphere(p0, p1, t);
      }
    }

    //********FACE LOCATED************/
    int sum =  ((o0 == COLLINEAR) ? 1 : 0)
             + ((o1 == COLLINEAR) ? 1 : 0)
             + ((o2 == COLLINEAR) ? 1 : 0);

    switch(sum)
    {
      case 0:
      {
        lt = FACE;
        li = 4;
        break;
      }
      case 1:
      {
        lt = EDGE;
        li = (o0 == COLLINEAR) ? 2 : (o1 == COLLINEAR) ? 0 : 1;
        break;
      }
      case 2:
      {
        lt = VERTEX;
        li = (o0 != COLLINEAR) ? 2 : (o1 != COLLINEAR) ? 0 : 1;
        break;
      }
      default:
      {
        // impossible
        CGAL_assertion(false);
        return f;
      }
    }

    test_distance(t, f, lt, li);

    return f;
  }
}

// @todo implement new, faster walks
template <typename Gt, typename Tds>
typename Triangulation_on_sphere_2<Gt, Tds>::Face_handle
Triangulation_on_sphere_2<Gt, Tds>::
locate(const Point& p,
       Locate_type& lt,
       int& li,
       Face_handle start) const
{
  if(!geom_traits().is_on_sphere(p))
  {
    lt = NOT_ON_SPHERE;
    return Face_handle();
  }

  switch(dimension())
  {
    case -2: // empty triangulation
    {
      lt = OUTSIDE_AFFINE_HULL;
      li = 4;
      return Face_handle();
    }
    case -1: // 1 vertex
      CGAL_FALLTHROUGH;
    case 0: // 2 vertices
    {
      for(Vertices_iterator vi=vertices_begin(); vi!=vertices_end(); ++vi)
      {
        if(are_equal(p, point(vi)))
        {
          lt = VERTEX;
          Face_handle f = vi->face();
          li = f->index(vi);
          return f;
        }
      }

      lt = OUTSIDE_AFFINE_HULL;
      li = 4;

      Face_handle f;
      test_distance(p, f, lt, li);

      return f;
    }
    case 1: // 3+ coplanar vertices
    {
      return march_locate_1D(p, lt, li);
    }
  }

  // below is dimension() == 2
  if(start == Face_handle())
    start = all_faces_begin();

  if(is_ghost(start))
  {
    for(All_faces_iterator it=all_faces_begin(); it !=all_faces_end(); ++it)
    {
      if(!is_ghost(it))
      {
        start = it;
        break;
      }
    }
  }

#if(! defined(CGAL_ZIG_ZAG_WALK)) && (! defined(CGAL_LFC_WALK))
  #define CGAL_ZIG_ZAG_WALK
#endif

#ifdef CGAL_ZIG_ZAG_WALK
  Face_handle res1;
  res1 = march_locate_2D(start, p, lt, li);
#endif

#ifdef CGAL_LFC_WALK
  Locate_type lt2;
  int li2;
  Face_handle res2 = march_locate_2D_LFC(start, p, lt2, li2);
#endif

#if defined(CGAL_ZIG_ZAG_WALK) && defined(CGAL_LFC_WALK)
  compare_walks(p, res1, res2, lt, lt2, li, li2);
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

template <typename Gt, typename Tds>
typename Triangulation_on_sphere_2<Gt, Tds>:: Face_handle
Triangulation_on_sphere_2<Gt, Tds>::
locate(const Point& p,
       const Face_handle start) const
{
  Locate_type lt;
  int li;
  return locate(p, lt, li, start);
}

//--------------------------------------------DEBUG-------------------------------------------------

template <typename Gt, typename Tds>
void
Triangulation_on_sphere_2<Gt, Tds>::
show_all() const
{
  //Triangulation_2::show_all();
  std::cerr << "PRINTING COMPLETE TRIANGULATION:" << std::endl;
  std::cerr << std::endl<< "====> " << this;
  std::cerr << "dimension " << dimension() << std::endl;
  std::cerr << "nb of vertices " << number_of_vertices() << std::endl;

  if(dimension() < 1)
    return;

  if(dimension() == 1)
  {
    std::cerr << " all edges dim 1 " << std::endl;
    for(All_edges_iterator aeit=all_edges_begin(); aeit!=all_edges_end(); ++aeit)
    {
      show_face(aeit->first);
      std::cerr << "   ------------   " << std::endl;
    }

    return;
  }

  std::cerr << " faces " << std::endl;
  for(All_faces_iterator fi = all_faces_begin(); fi !=all_faces_end(); ++fi)
  {
    show_face(fi);
    std::cerr << "   ------------   " << std::endl;
  }

  if(number_of_vertices() > 1)
  {
    std::cerr << "print triangulation vertices:" << std::endl;
    for(Vertices_iterator vi=vertices_begin(); vi!=vertices_end(); ++vi)
    {
      show_vertex(vi);
      std::cerr << "  / associated face: " << (void*)(&(*(vi->face()))) << std::endl;;
    }

    std::cerr << std::endl;
  }
  return;
}

template <typename Gt, typename Tds>
void
Triangulation_on_sphere_2<Gt, Tds>::
show_vertex(Vertex_handle vh) const
{
  std::cerr << point(vh) << "\t";
  return;
}

template <typename Gt, typename Tds>
void
Triangulation_on_sphere_2<Gt, Tds>::
show_face(Face_handle fh) const
{
  std::cerr << "face : " << (void*)&(*fh) << " => " << std::endl;
  if(is_ghost(fh))
    std::cerr << "ghost " << std::endl;

  int i = fh->dimension();
  switch(i)
  {
    case 0:
      std::cerr << "point :" ; show_vertex(fh->vertex(0));
      std::cerr << " / neighbor " << &(*(fh->neighbor(0)));
      std::cerr << "[" ; show_vertex(fh->neighbor(0)->vertex(0));
      std::cerr << "]"  << std::endl;
      break;

    case 1:
      std::cerr << "point :" ; show_vertex(fh->vertex(0));
      std::cerr << " / neighbor " << &(*(fh->neighbor(0)));
      std::cerr << "[" ; show_vertex(fh->neighbor(0)->vertex(0));
      std::cerr << "/" ; show_vertex(fh->neighbor(0)->vertex(1));
      std::cerr << "]" << std::endl;

      std::cerr << "point :" ; show_vertex(fh->vertex(1));
      std::cerr << " / neighbor " << &(*(fh->neighbor(1)));
      std::cerr << "[" ; show_vertex(fh->neighbor(1)->vertex(0));
      std::cerr << "/" ; show_vertex(fh->neighbor(1)->vertex(1));
      std::cerr << "]" << std::endl;
      break;

    case 2:
      std::cerr << "point :" ; show_vertex(fh->vertex(0));
      std::cerr << " / neighbor " << &(*(fh->neighbor(0)));
      std::cerr << "[" ; show_vertex(fh->neighbor(0)->vertex(0));
      std::cerr << "/" ; show_vertex(fh->neighbor(0)->vertex(1));
      std::cerr << "/" ; show_vertex(fh->neighbor(0)->vertex(2));
      std::cerr << "]" << std::endl;

      std::cerr << "point :" ; show_vertex(fh->vertex(1));
      std::cerr << " / neighbor " << &(*(fh->neighbor(1)));
      std::cerr << "[" ; show_vertex(fh->neighbor(1)->vertex(0));
      std::cerr << "/" ; show_vertex(fh->neighbor(1)->vertex(1));
      std::cerr << "/" ; show_vertex(fh->neighbor(1)->vertex(2));
      std::cerr << "]" << std::endl;

      std::cerr << "point :" ; show_vertex(fh->vertex(2));
      std::cerr << " / neighbor " << &(*(fh->neighbor(2)));
      std::cerr << "[" ; show_vertex(fh->neighbor(2)->vertex(0));
      std::cerr << "/" ; show_vertex(fh->neighbor(2)->vertex(1));
      std::cerr << "/" ; show_vertex(fh->neighbor(2)->vertex(2));
      std::cerr << "]" << std::endl;
      break;
  }
}

// --------------------------CHECKING---------------------------

// checks whether neighboring faces are linked correctly to each other.
template <typename Gt, typename Tds>
void
Triangulation_on_sphere_2<Gt, Tds>::
check_neighboring_faces() const
{
  if(dimension() == 1)
  {
    for(All_faces_iterator eit=all_faces_begin(); eit!=all_faces_end(); ++eit)
    {
      CGAL_assertion_code(const Face_handle f1 = eit->neighbor(0);)
      CGAL_assertion_code(const Face_handle f2 = eit->neighbor(1);)
      CGAL_assertion(f1->has_neighbor(eit));
      CGAL_assertion(f2->has_neighbor(eit));
    }
  }

  for(All_faces_iterator eit=all_faces_begin(); eit!=all_faces_end(); ++eit)
  {
    CGAL_assertion_code(const Face_handle f1 = eit->neighbor(0);)
    CGAL_assertion_code(const Face_handle f2 = eit->neighbor(1);)
    CGAL_assertion_code(const Face_handle f3 = eit->neighbor(2);)
    CGAL_assertion(f1->has_neighbor(eit));
    CGAL_assertion(f2->has_neighbor(eit));
    CGAL_assertion(f3->has_neighbor(eit));
  }
}

// checks whether a given triangulation is plane (all points are coplanar)
template <typename Gt, typename Tds>
bool
Triangulation_on_sphere_2<Gt, Tds>::
is_plane() const
{
  if(number_of_vertices() <= 3)
    return true;

  bool is_plane = true;

  Vertices_iterator it1 = vertices_begin(), it2(it1), it3(it1), it4(it1);
  std::advance(it2, 1);
  std::advance(it3, 2);
  std::advance(it4, 3);

  while(it4 != vertices_end())
  {
    Orientation s = orientation(point(it1), point(it2), point(it3), point(it4));
    is_plane = is_plane && s == COPLANAR;

    if(!is_plane)
      return false;

    ++it1;
    ++it2;
    ++it3;
    ++it4;
  }

  return true;
}

template <typename Gt, typename Tds>
bool
Triangulation_on_sphere_2<Gt, Tds>::
is_valid_vertex(Vertex_handle vh, bool verbose, int /*level*/) const
{
  bool result = vh->face()->has_vertex(vh);
  if(!result)
  {
    if(verbose)
    {
      std::cerr << " from is_valid_vertex " << std::endl;
      std::cerr << "normal vertex " << &(*vh) << std::endl;
      show_vertex(vh);
      std::cerr << "\nvh_>face " << &*(vh->face()) << " " << std::endl;

      show_face(vh->face());
    }

    CGAL_assertion(false);
    return false;
  }

  return true;
}

template <typename Gt, typename Tds>
bool
Triangulation_on_sphere_2<Gt, Tds>::
is_valid(bool verbose,
         int level) const
{
  bool result = _tds.is_valid(verbose, level);
  if(dimension() <= 0 || (dimension() == 1 && number_of_vertices() == 2))
    return result;

  for(Vertices_iterator vit=vertices_begin(); vit!=vertices_end(); ++vit)
    result = result && is_valid_vertex(vit, verbose, level);

  if(dimension() == 1)
  {
    result = result && this->is_plane();
    CGAL_assertion(result);
  }
  else // dimension() == 2
  {
    for(All_faces_iterator it=all_faces_begin(); it!=all_faces_end(); ++it)
    {
      const Orientation s = orientation_on_sphere(point(it, 0), point(it, 1), point(it, 2));
      result = result && (s == LEFT_TURN || is_ghost(it));
      CGAL_assertion(result);
    }

    // check number of faces. This cannot be done by the TDS,
    // which does not know the number of components nor the genus
    result = result && (number_of_faces() == (2 * number_of_vertices() - 4));

    CGAL_assertion(result);
  }

  return result;
}

// ---------------------------------I/O--------------------------- //

template <typename Gt, typename Tds>
void
Triangulation_on_sphere_2<Gt, Tds>::
file_output(std::ostream& os) const
{
  os << _gt.center() << " " << _gt.radius() << "\n";
  _tds.file_output(os, Vertex_handle(), true);
}

template <typename Gt, typename Tds>
typename Triangulation_on_sphere_2<Gt, Tds>::Vertex_handle
Triangulation_on_sphere_2<Gt, Tds>::
file_input(std::istream& is)
{
  clear();

  Point_3 center;
  FT radius;
  is >> center >> radius;
  _gt.set_center(center);
  _gt.set_radius(radius);

  Vertex_handle v = _tds.file_input(is, false);
  CGAL_assertion(is_valid());
  return v;
}

template <typename Gt, typename Tds>
std::ostream&
operator<<(std::ostream& os, const Triangulation_on_sphere_2<Gt, Tds>& tr)
{
  tr.file_output(os);
  return os ;
}

template < class Gt, class Tds >
std::istream&
operator>>(std::istream& is, Triangulation_on_sphere_2<Gt, Tds>& tr)
{
  tr.file_input(is);
  return is;
}

} // namespace CGAL

#endif // CGAL_TRIANGULATION_ON_SPHERE_2_H
