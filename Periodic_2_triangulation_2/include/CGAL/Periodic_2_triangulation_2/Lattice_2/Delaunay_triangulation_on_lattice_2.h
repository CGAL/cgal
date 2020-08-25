// Copyright (c) 2020 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé
//                 Georg Osang

// See also publication:
//   Georg Osang, Mael Rouxel-Labbé, and Monique Teillaud
//   Generalizing CGAL Periodic Delaunay Triangulations
//   ESA 2020

#ifndef CGAL_P2T2_DELAUNAY_TRIANGULATION_ON_LATTICE_2_H
#define CGAL_P2T2_DELAUNAY_TRIANGULATION_ON_LATTICE_2_H

#include <CGAL/license/Periodic_2_triangulation_2.h>

#define CGAL_GENERIC_P2T2

#include <CGAL/Periodic_2_triangulation_2/Lattice_2/Triangulation_face_base_on_lattice_2.h>
#include <CGAL/Periodic_2_triangulation_2/Lattice_2/Triangulation_iterators_on_lattice_2.h>
#include <CGAL/Periodic_2_triangulation_2/Lattice_2/Triangulation_circulators_on_lattice_2.h>
#include <CGAL/Periodic_2_triangulation_2/Lattice_2/Lattice_2.h>
#include <CGAL/Periodic_2_triangulation_2/Square_flat_torus_2/Delaunay_triangulation_on_square_flat_torus_2.h>
#include <CGAL/Periodic_2_triangulation_vertex_base_2.h>
#include <CGAL/Periodic_2_triangulation_2/IO/periodic_2_triangulation_2_io.h>
#include <CGAL/draw_triangulation_2.h>

#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/Triangulation_utils_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/utility.h>

#include <utility>
#include <iostream>
#include <fstream>
#include <string>

namespace CGAL {

// TODO
// - Move p2_iterators to Square_flat_torus_2
// - Don't document Finite iterators / functions (but keep them)
// -

template <typename Gt_,
          typename Tds_ = CGAL::Default>
class Delaunay_triangulation_on_lattice_2
  : public Triangulation_cw_ccw_2
{
public:
  typedef Gt_                                                             Geom_traits;

  typedef typename Default::Get<Tds_,
                                Triangulation_data_structure_2<
                                  Periodic_2_triangulation_vertex_base_2<Gt_>,
                                  Triangulation_face_base_on_lattice_2<Gt_> > >::type
                                                                          Triangulation_data_structure;

private:
  typedef Geom_traits                                                     Gt;
  typedef Triangulation_data_structure                                    Tds;
  typedef Delaunay_triangulation_on_lattice_2<Gt, Tds>                    Self;

public:
  typedef Delaunay_triangulation_2<Gt, Tds>                               DT2;

  // Once the periodic triangulation fits on a single cover, we can (ab)use
  // the code initially meant for periodic triangulations on iso rectangles
  typedef Delaunay_triangulation_on_square_flat_torus_2<Gt, Tds>          P2DT2;

public:
  typedef typename Gt::Domain                                             Domain;
  typedef std::array<int, 2>                                              Covering_sheets;

  typedef typename Geom_traits::FT                                        FT;
  typedef typename Geom_traits::Segment_2                                 Segment;
  typedef typename Geom_traits::Vector_2                                  Vector;
  typedef typename Geom_traits::Triangle_2                                Triangle;

public:
  typedef typename Triangulation_data_structure::size_type                size_type;
  typedef typename Tds::Vertex                                            Vertex;
  typedef typename Tds::Vertex_handle                                     Vertex_handle;
  typedef typename Tds::Edge                                              Edge;
  typedef typename Tds::Face_handle                                       Face_handle;

  typedef typename Vertex::Point                                          Point;

  typedef typename Geom_traits::Periodic_2_offset_2                       Offset;
  typedef std::pair<Point, Offset>                                        Periodic_point;
  typedef array<std::pair<Point, Offset>, 2>                              Periodic_segment;
  typedef array<std::pair<Point, Offset>, 3>                              Periodic_triangle;
  typedef array<std::pair<Point, Offset>, 4>                              Periodic_tetrahedron;

  typedef Triangulation_ds_vertex_iterator_on_lattice_2<Self>             Vertex_iterator;
  typedef Triangulation_ds_edge_iterator_on_lattice_2<Self>               Edge_iterator;
  typedef Triangulation_ds_face_iterator_on_lattice_2<Self>               Face_iterator;

  typedef Triangulation_ds_vertex_circulator_on_lattice_2<Self>           Vertex_circulator;
  typedef Triangulation_ds_edge_circulator_on_lattice_2<Self>             Edge_circulator;
  typedef Triangulation_ds_face_circulator_on_lattice_2<Self>             Face_circulator;

private:
  typedef Construct_periodic_point<Self>                                  Construct_PP;
  typedef Construct_periodic_segment<Self>                                Construct_PS;
  typedef Construct_periodic_triangle<Self>                               Construct_PT;

public:
  typedef boost::transform_iterator<Construct_PP, Vertex_iterator>        Periodic_point_iterator;
  typedef boost::transform_iterator<Construct_PS, Vertex_iterator>        Periodic_segment_iterator;
  typedef boost::transform_iterator<Construct_PT, Vertex_iterator>        Periodic_triangle_iterator;

  //Tag to distinguish Delaunay from regular triangulations
  typedef Tag_false                                                       Weighted_tag;

  // Tag to distinguish periodic triangulations from others
  typedef Tag_true                                                        Periodic_tag;

  enum Locate_type
  {
    VERTEX = 0,
    EDGE, // 1
    FACE, // 2
    EMPTY , // 3
    OUTSIDE_CONVEX_HULL, // unused, for compatibility with Alpha shapes
    OUTSIDE_AFFINE_HULL // unused, for compatibility with Alpha shapes
  };

  enum Iterator_type
  {
    IT_DT2 = 0,
    IT_P2T2 // 1
  };

private:
  Geom_traits gt_;

  bool is_1_cover_;
  FT sq_circumradius_threshold;
  std::set<Face_handle> faces_with_too_big_circumradius;

  DT2 dt2;
  P2DT2 p2dt2;
  std::unordered_map<Vertex_handle, Vertex_handle> canonical_vertices;

public:
  /// Constructors
  Delaunay_triangulation_on_lattice_2(const Geom_traits& gt)
    :  gt_(gt), is_1_cover_(false), dt2(gt_), p2dt2(gt_)
  {
    sq_circumradius_threshold = 0.0625 * gt.get_domain().systole_sq_length();
  }

  Delaunay_triangulation_on_lattice_2(const Domain& lattice)
    : Delaunay_triangulation_on_lattice_2(Geom_traits(lattice))
  { }

  template <class InputIterator>
  Delaunay_triangulation_on_lattice_2(InputIterator first, InputIterator beyond,
                                      const Geom_traits& gt)
    : gt_(gt), is_1_cover_(false), dt2(gt_), p2dt2(gt_)
  {
    sq_circumradius_threshold = 0.0625 * gt.get_domain().systole_sq_length();
    insert(first, beyond);
  }

  template <class InputIterator>
  Delaunay_triangulation_on_lattice_2(InputIterator first, InputIterator beyond,
                                      const Domain& lattice)
    : Delaunay_triangulation_on_lattice_2(first, beyond, Geom_traits(lattice))
  { }

  // @todo (virtual vertices cannot be naively copied)
  Delaunay_triangulation_on_lattice_2& operator=(const Delaunay_triangulation_on_lattice_2&) = delete;
  void copy_triangulation(const Delaunay_triangulation_on_lattice_2 &tr) = delete;
  void swap(Delaunay_triangulation_on_lattice_2 &tr) = delete;

  void clear()
  {
    dt2.clear();
    p2dt2.clear();

    is_1_cover_ = false;
    canonical_vertices.clear();
    faces_with_too_big_circumradius.clear();
  }

  /// Access
  const Geom_traits& geom_traits() const { return gt_; }
  const Triangulation_data_structure& tds() const
  {
    if(is_1_cover())
      return p2dt2.tds();
    else
      return dt2.tds();
  }

  int dimension() const { return (number_of_vertices() == 0) ? -2 : 2; }

  // This triangulation, even when internally using a DT2 with duplicated points,
  // mimics a single-sheet periodic triangulation
  Covering_sheets number_of_sheets() const { return make_array(1, 1); }

public:
  bool is_1_cover() const { return is_1_cover_; }
  size_t too_big_faces() const { return faces_with_too_big_circumradius.size(); }

public:
  template<class T>
  bool is_infinite(const T&, int = 0, int = 0) const { return false; }

  /// Canonicity
  Offset compute_offset(const Edge& e) const
  {
    Face_handle f = e.first;
    int i = e.second;
    return compute_offset(f->vertex(dt2.cw(i)), f->vertex(dt2.ccw(i)));
  }

  Offset compute_offset(const Face_handle f, const int i) const
  {
    return compute_offset(f->vertex(dt2.cw(i)), f->vertex(dt2.ccw(i)));
  }

  // Offset of an edge represented by its two vertices
  Offset compute_offset(const Vertex_handle v1, const Vertex_handle v2) const
  {
    // @todo proper min implementation
    return min(v1->offset(), v2->offset());
  }

  Offset compute_offset(const Face_handle f) const
  {
    return min(f->vertex(0)->offset(), f->vertex(1)->offset(), f->vertex(2)->offset());
  }

  bool is_canonical(const Point& p) const
  {
    return gt_.get_domain().is_in_scaled_domain(p, 1);
  }

  bool is_canonical(const Vertex_handle v) const
  {
    if(is_1_cover())
      return true;

    return (v->offset() == Offset(0, 0));
  }

  bool is_canonical(const Edge& e) const
  {
    if(is_1_cover())
      return true;

    if(dt2.is_infinite(e.first))
      return false;

    return compute_offset(e) == Offset(0, 0);
  }

  bool is_canonical(const Face_handle f) const
  {
    if(is_1_cover())
      return true;

    if(dt2.is_infinite(f))
      return false;

    return compute_offset(f) == Offset(0, 0);
  }

  template <typename ForwardFaceIterator>
  void mark_canonical_faces(ForwardFaceIterator fit, ForwardFaceIterator beyond)
  {
    CGAL_precondition(!is_1_cover());
    for(; fit != beyond; ++fit)
      fit->set_canonical_flag(is_canonical(fit));
  }

  void mark_canonical_faces()
  {
    CGAL_precondition(!is_1_cover());
    return mark_canonical_faces(dt2.finite_faces_begin(), dt2.finite_faces_end());
  }

  void mark_canonical_faces(Vertex_handle v)
  {
    CGAL_precondition(!is_1_cover());
    CGAL_precondition(v != Vertex_handle());

    if(dt2.dimension() != 2)
      return;

    typename DT2::Face_circulator fc = dt2.incident_faces(v), done = fc;
    do
    {
      fc->set_canonical_flag(is_canonical(fc));
    }
    while(++fc != done);
  }

  void reset_all_canonicity()
  {
    CGAL_precondition(!is_1_cover());
    typename DT2::Finite_faces_iterator fit = dt2.faces_begin(), fend = dt2.faces_end();
    for(; fit!=fend; ++fit)
      fit->set_canonical_flag(false);
  }

  Vertex_handle canonical_vertex(const Vertex_handle v) const
  {
    return (is_canonical(v) ? v : canonical_vertices.at(v));
  }

  /// Canonicalization
  Point construct_canonical_point(const Point& p) const
  {
    return geom_traits().construct_canonical_point(p);
  }

  template <class InputIterator>
  std::vector<Point> construct_canonical_points(InputIterator first, InputIterator beyond) const
  {
    std::vector<Point> canonical_points;

    while(first != beyond)
    {
      const Point& p = *first++;
      canonical_points.push_back(construct_canonical_point(p));
    }

    return canonical_points;
  }

  /// Given a face having a vertex in the domain, and an offset such that
  /// at least one of the translated vertices has a vertex in the domain,
  /// return the handle of the translated face.
  ///
  /// Returns Face_handle() if the provided face is not a periodic face
  Face_handle find_translated_face(Face_handle f, const Offset& o) const
  {
    // The code commented out below does the same and is simpler,
    // but uses a construction and point location.
    //   Point translated_barycenter = gt_.get_domain().translate_by_offset(construct_barycenter(f), o);
    //   return dt2.locate(translated_barycenter);

    // Find a vertex whose translate is in the domain.
    bool vertex_found = false;
    int j=0;
    for(; j<3; ++j)
    {
      if(f->vertex(j)->offset() == -o)
      {
        vertex_found = true;
        break;
      }
    }

    CGAL_assertion_msg(vertex_found, "Invalid offset for translation of face.");
    Vertex_handle cv = canonical_vertices.at(f->vertex(j));
    Vertex_handle v_ccw = f->vertex(dt2.ccw(j));
    Vertex_handle v_cw = f->vertex(dt2.cw(j));

    // Scan through the incident faces and find the one that is
    // equivalent to f.
    typename DT2::Face_circulator fc = dt2.incident_faces(cv), done = fc;
    do
    {
      CGAL_assertion(!dt2.is_infinite(fc));

      int cj = fc->index(cv);
      Vertex_handle cv_ccw = fc->vertex(dt2.ccw(cj));
      Vertex_handle cv_cw = fc->vertex(dt2.cw(cj));

      if(canonical_vertices.at(cv_ccw) == canonical_vertices.at(v_ccw) &&
         cv_ccw->offset() == v_ccw->offset() + o &&
         canonical_vertices.at(cv_cw) == canonical_vertices.at(v_cw) &&
         cv_cw->offset() == v_cw->offset() + o)
      {
        return Face_handle(fc);
      }
    }
    while(++fc != done);

    // provided face is not a periodic face and doesn't have a copy with a vertex in the domain.
    return Face_handle();
  }

  /// Returns Face_handle() if the provided face is not a periodic face
  Face_handle get_canonical_face(Face_handle f) const
  {
    return find_translated_face(f, -compute_offset(f));
  }

  // @todo something smarter, this function is core to everything else
  // How could we keep "canonical neighbors" in memory instead of having
  // to find them later...?
  //
  // @cache this value?
  Face_handle neighbor(Face_handle f, int i) const
  {
    if(is_canonical(f->neighbor(i)))
      return f->neighbor(i);

    // Translate the face so that the corresponding edge is canonical.
    Offset edge_off = compute_offset(f, i);
    Face_handle f_trans = find_translated_face(f, -edge_off);
    CGAL_assertion(f_trans != Face_handle());

    // Find the vertex corresponding to vertex(i) in the original.
    bool vertex_found = false;
    int j = 0;
    for(; j<3; ++j)
    {
      if(canonical_vertices.at(f_trans->vertex(j)) == canonical_vertices.at(f->vertex(i)) &&
         f_trans->vertex(j)->offset() == f->vertex(i)->offset() - edge_off)
      {
        vertex_found = true;
        break;
      }
    }

    CGAL_assertion(vertex_found);

    // Get the neighbour in DT2 and check if it's canonical.
    Face_handle neighbor = f_trans->neighbor(j);
    if(is_canonical(neighbor))
      return neighbor;
    else // if not, find the canonical translate.
      return get_canonical_face(neighbor);
  }

  // -----------------------------------------------------------------------------------------------
  /// Returns the number of vertices. Counts all vertices that are
  /// representatives of the same point in the 1-cover as one vertex.
  size_type number_of_vertices() const
  {
    if(is_1_cover())
      return p2dt2.number_of_vertices();
    else
      return dt2.number_of_vertices() / 9; // @fixme wrong?
  }

  size_type number_of_edges() const
  {
    if(is_1_cover())
    {
      return p2dt2.number_of_edges();
    }
    else
    {
      // Exploiting Euler's formula that #f - #e + #v = 0 on the torus
      return number_of_faces() + number_of_vertices();
    }
  }

  size_type number_of_faces() const
  {
    if(is_1_cover())
    {
      return p2dt2.number_of_faces();
    }
    else
    {
      // @todo something better than this naive way (flag edge/face classes)
      size_type nf = 0;
      typename DT2::Finite_faces_iterator fit = dt2.finite_faces_begin(), fend = dt2.finite_faces_end();
      for(; fit!=fend; ++fit)
        if(is_canonical(fit))
          ++nf;

      return nf;
    }
  }

  size_type number_of_stored_vertices() const
  {
    if(is_1_cover())
      return p2dt2.number_of_vertices();
    else
      return dt2.number_of_vertices();
  }

  size_type number_of_stored_edges() const
  {
    if(is_1_cover())
      return p2dt2.number_of_edges();
    else
      return dt2.number_of_edges();
  }

  size_type number_of_stored_faces() const
  {
    if(is_1_cover())
      return p2dt2.number_of_faces();
    else
      return dt2.number_of_faces();
  }

  size_type number_of_finite_vertices() const { return number_of_vertices(); }
  size_type number_of_finite_edges() const { return number_of_edges(); }
  size_type number_of_finite_faces() const { return number_of_faces(); }

  /// Iterators and Circulators
  Vertex_iterator vertices_begin() const { return Vertex_iterator(this); }
  Vertex_iterator vertices_end() const { return Vertex_iterator(this, 1); }
  Edge_iterator edges_begin() const { return Edge_iterator(this); }
  Edge_iterator edges_end() const { return Edge_iterator(this, 1); }
  Face_iterator faces_begin() const { return Face_iterator(this); }
  Face_iterator faces_end() const { return Face_iterator(this, 1); }

  Periodic_point_iterator periodic_points_begin() const { return Periodic_point_iterator(vertices_begin()); }
  Periodic_point_iterator periodic_points_end() const { return Periodic_point_iterator(vertices_end()); }
  Periodic_segment_iterator periodic_segments_begin() const { return Periodic_segment_iterator(edges_begin()); }
  Periodic_segment_iterator periodic_segments_end() const { return Periodic_segment_iterator(edges_end()); }
  Periodic_triangle_iterator periodic_triangles_begin() const { return Periodic_triangle_iterator(faces_begin()); }
  Periodic_triangle_iterator periodic_triangles_end() const { return Periodic_triangle_iterator(faces_end()); }

  Vertex_circulator adjacent_vertices(Vertex_handle v, Face_handle f) const
  {
    if(!is_canonical(v))
      v = canonical_vertices.at(v);

    return Vertex_circulator(v, f, this);
  }

  Vertex_circulator adjacent_vertices(Vertex_handle v) const
  {
    return adjacent_vertices(v, Face_handle());
  }

  Edge_circulator incident_edges(Vertex_handle v, Face_handle f) const
  {
    if(!is_canonical(v))
      v = canonical_vertices.at(v);

    return Edge_circulator(v, f, this);
  }

  Edge_circulator incident_edges(Vertex_handle v) const
  {
    return incident_edges(v, Face_handle());
  }

  Face_circulator incident_faces(Vertex_handle v, Face_handle f) const
  {
    if(!is_canonical(v))
      v = canonical_vertices.at(v);

    return Face_circulator(v, f, this);
  }

  Face_circulator incident_faces(Vertex_handle v) const
  {
    return incident_faces(v, Face_handle());
  }

  // -----------------------------------------------------------------------------------------------

  /// Predicates
  FT compute_squared_circumradius(const Face_handle f) const
  {
    return geom_traits().compute_squared_radius_2_object()(dt2.point(f, 0),
                                                           dt2.point(f, 1),
                                                           dt2.point(f, 2));
  }

  /// Constructions
  Point construct_point(const Point& p, const Offset& off) const
  {
    return geom_traits().construct_point_2_object()(p, off);
  }

  Point construct_point(const Periodic_point& pp) const
  {
    return construct_point(pp.first, pp.second);
  }

  Segment construct_segment(const Point& p1, const Point& p2,
                            const Offset& o1, const Offset& o2) const
  {
    return geom_traits().construct_segment_2_object()(p1, p2, o1, o2);
  }

  Segment construct_segment(const Periodic_segment& seg) const
  {
    return construct_segment(seg[0].first, seg[1].first, seg[0].second, seg[1].second);
  }

  Triangle construct_triangle(const Point& p1, const Point& p2, const Point& p3,
                              const Offset& o1, const Offset& o2, const Offset& o3) const
  {
    return geom_traits().construct_triangle_2_object()(p1, p2, p3, o1, o2, o3);
  }

  Triangle construct_triangle(const Periodic_triangle& tri) const
  {
    return construct_triangle(tri[0].first, tri[1].first, tri[2].first,
                              tri[0].second, tri[1].second, tri[2].second);
  }

  // -----------------------------------------------------------------------------------------------
  Periodic_point periodic_point(const Vertex_handle v) const
  {
    if(is_1_cover())
      return p2dt2.periodic_point(v);
    else
      return Periodic_point(v->point().first, v->point().second);
  }

  Periodic_point periodic_point(const Face_handle f, int i) const
  {
    if(is_1_cover())
      return p2dt2.periodic_point(f, i);
    else
      return Periodic_point(f->vertex(i)->point().first, f->vertex(i)->point().second);
  }

  Periodic_segment periodic_segment(const Face_handle f, int i) const
  {
    if(is_1_cover())
      return p2dt2.periodic_segment(f, i);
    else
      return make_array(periodic_point(f, ccw(i)), periodic_point(f, cw(i)));
  }

  Periodic_segment periodic_segment(const Edge &e) const
  {
    return periodic_segment(e.first, e.second);
  }

  Periodic_triangle periodic_triangle(Face_handle f) const
  {
    return make_array(periodic_point(f, 0), periodic_point(f, 1), periodic_point(f, 2));
  }

  // -----------------------------------------------------------------------------------------------
  Point point(const Periodic_point& pp) const
  {
    return construct_point(pp.first, pp.second);
  }

  Point point(const Vertex_handle v) const
  {
    return point(periodic_point(v));
  }

  Point point(const Face_handle f, int i) const
  {
    return point(periodic_point(f, i));
  }

  Segment segment(const Periodic_segment& ps) const
  {
    return construct_segment(ps[0].first, ps[1].first, ps[0].second, ps[1].second);
  }

  Segment segment(Face_handle f, int i) const
  {
    return segment(periodic_segment(f, i));
  }

  Segment segment(const Edge& e) const
  {
    return segment(periodic_segment(e));
  }

  Segment segment(const Edge_circulator& ec) const { return segment(*ec); }
  Segment segment(const Edge_iterator& ei) const { return segment(*ei); }

  Triangle triangle(const Periodic_triangle& pt) const
  {
    return construct_triangle(pt[0].first, pt[1].first, pt[2].first,
                              pt[0].second, pt[1].second, pt[2].second);
  }

  Triangle triangle(Face_handle f) const
  {
    return triangle(periodic_triangle(f));
  }

  /// Locate functions
//  Face_handle locate(const Point& p, Offset& lo,
//                     Locate_type& lt, int& li,
//                     Face_handle start = Face_handle()) const
//  { }

  Face_handle locate(const Point& p,
                     Locate_type& lt, int& li,
                     Face_handle start = Face_handle()) const
  {
    if(is_1_cover())
      return p2dt2.locate(p, lt, li, start);
    else
      return dt2.locate(p, lt, li, start);
  }

  Face_handle locate(const Point& p,
                     Face_handle start = Face_handle()) const
  {
    if(is_1_cover())
      return p2dt2.locate(p, start);
    else
      return dt2.locate(p, start);
  }

  Vertex_handle nearest_vertex(const Point& p, Face_handle f = Face_handle())
  {
    // @todo
    CGAL_assertion(false);
  }

  /// Insertion and removal
  Vertex_handle insert_in_dt2(const Point& p)
  {
    CGAL_precondition(!is_1_cover_);
    CGAL_precondition(gt_.get_domain().is_in_scaled_domain(p));

    if(dt2.dimension() >= 2) // equivalent to !dt2.empty() since we insert duplicate vertices
    {
      // @todo avoid recomputing the conflict zone if possible (done also 'insert', sort of)
      std::vector<Face_handle> faces_in_conflict;
      dt2.get_conflicts(p, std::back_inserter(faces_in_conflict));
//      std::cout << faces_in_conflict.size() << " faces in conflict" << std::endl;

      size_t erased_faces = 0;
      for(Face_handle f : faces_in_conflict)
      {
        Face_handle cf = get_canonical_face(f);

        // We might get non-periodic "boundary" faces that don't have a canonical version
        if(cf != Face_handle())
        {
          // Some faces might appear multiple times in the conflict zone, but this is fine.
          erased_faces += faces_with_too_big_circumradius.erase(cf);
        }
      }
//      std::cout << "Faces with too big radius in conflict zone:" << erased_faces << std::endl;
    }

    Vertex_handle v = dt2.insert(p);

    CGAL_assertion(v != Vertex_handle());
    v->set_offset(Offset(0, 0));

    mark_canonical_faces(v);

    canonical_vertices[v] = v;
    for(const std::array<int, 2>& off : gt_.get_domain().overlapping_offsets)
    {
      // @fixme construction here means potential problems
      // solution is prob to template DT2 with a GT where Point := pair<Point, Offset>
      // and then do the full stack of filtering like for periodic traits
      // and then insert(make_pair(p, off))
      // ---> once again requires an exact lattice...
      const Vector off_v = gt_.construct_sum_of_vectors_2_object()(
                             gt_.construct_scaled_vector_2_object()(gt_.get_domain().basis()[0], off[0]),
                             gt_.construct_scaled_vector_2_object()(gt_.get_domain().basis()[1], off[1]));
      const Point off_p = p + off_v;

      if(gt_.get_domain().is_in_scaled_domain(off_p, 3))
      {
        Vertex_handle v_copy = dt2.insert(off_p);
        CGAL_assertion(v_copy != Vertex_handle());
        v_copy->set_offset(Offset(off[0], off[1]));

        canonical_vertices[v_copy] = v;

        mark_canonical_faces(v_copy);
      }
    }

    // Update the current maximum circumradius value
    typename DT2::Face_circulator fc = dt2.incident_faces(v), done(fc);
    do
    {
      CGAL_assertion(!dt2.is_infinite(fc));

      Face_handle cf = get_canonical_face(fc);
      CGAL_assertion(cf != Face_handle());

      const FT sq_cr = compute_squared_circumradius(cf);
      if(sq_cr > sq_circumradius_threshold)
        faces_with_too_big_circumradius.insert(cf);
    }
    while(++fc != done);

    if(faces_with_too_big_circumradius.empty())
       convert_to_1_cover();

    return v;
  }

  Vertex_handle insert_in_p2dt2(const Point& p)
  {
    CGAL_precondition(is_1_cover_);
    return p2dt2.insert(p);
  }

  Vertex_handle insert(const Point& p)
  {
    const Point cp = gt_.get_domain().construct_canonical_point(p);

    if(is_1_cover())
      return insert_in_p2dt2(cp);
    else
      return insert_in_dt2(cp);
  }

  Vertex_handle push_back(const Point& p)
  {
    return insert(p);
  }

  // @todo mirror it with other insert functions (e.g. spatial sorting)
  template <class InputIterator>
  void insert(InputIterator first, InputIterator beyond)
  {
    while(first != beyond)
      insert(*first++);

    CGAL_postcondition(dt2.dimension() == 2);
    CGAL_postcondition(dt2.is_valid());
  }

  void remove(Vertex_handle v)
  {
    // @todo
    CGAL_assertion(false);
  }

  Vertex_handle move_if_no_collision(Vertex_handle v, const Point& p)
  {
    Locate_type lt;
    int li;
    Vertex_handle inserted;
    Face_handle loc = locate(p, lt, li, v->face());

    if(lt == VERTEX)
      return v;
    else
      /// This can be optimized by checking whether we can move v->point() to p
      return insert(p, lt, loc, li);
  }

  Vertex_handle move_point(Vertex_handle v, const Point& p)
  {
    if(v->point() == p)
      return v;

    Vertex_handle w = move_if_no_collision(v, p);
    if(w != v)
    {
      remove(v);
      return w;
    }
    return v;
  }

  /// Covering

  // @todo two versions, one using the sufficient condition, one checking for real.
  // "For real" --> P3T3 is_embeddable_in_...
  // "For real" is implemented below.
  bool is_simplicial_complex() const
  {
    // Ensure that there is no edge between a vertex and itself
    typename DT2::Finite_faces_iterator fit = dt2.faces_begin(), fend = dt2.faces_end();
    for(; fit!=fend; ++fit)
    {
      if(!is_canonical(fit))
        continue;

      Vertex_handle v0 = canonical_vertex(fit->vertex(0));
      Vertex_handle v1 = canonical_vertex(fit->vertex(1));
      Vertex_handle v2 = canonical_vertex(fit->vertex(2));

      if(v0 == v1 || v0 == v2 || v1 == v2)
        return false;
    }

    // Ensure that there are no two edges between the same pair of vertices
    typename DT2::Finite_vertices_iterator vit = dt2.vertices_begin(),
                                           vend = dt2.vertices_end();
    for(; vit!=vend; ++vit)
    {
      if(!is_canonical(vit))
        continue;

      typename DT2::Vertex_circulator vc = dt2.incident_vertices(vit), done = vc;
      std::unordered_set<Vertex_handle> neighbours;
      do
      {
        Vertex_handle cv = canonical_vertex(vc);
        if(neighbours.find(cv) != neighbours.end())
          return false; // some neighbouring vertex appeared multiple times

        neighbours.insert(cv);
      }
      while(++vc != done);
    }

    return true;
  }

  void add_edge_to_incident_faces_map(const Face_handle f, int i,
                                      std::map<std::set<Vertex_handle>,
                                      std::vector<std::pair<Face_handle, int> > >& incident_faces_map)
  {
    typedef std::set<Vertex_handle>                                            Edge_vertices;
    typedef std::pair<Face_handle, int>                                        Incident_face;
    typedef std::map<Edge_vertices, std::vector<Incident_face> >               Incident_faces_map;

    // the opposite vertex of f in c is i
    Edge_vertices e;
    e.insert(f->vertex((i + 1) % 3));
    e.insert(f->vertex((i + 2) % 3));
    CGAL_precondition(e.size() == 2);

    Incident_face icf = std::make_pair(f, i);
    std::vector<Incident_face> vec;
    vec.push_back(icf);

    std::pair<typename Incident_faces_map::iterator, bool> is_insert_successful =
        incident_faces_map.insert(std::make_pair(e, vec));
    if(!is_insert_successful.second) // the entry already exists in the map
    {
      // a facet must have exactly two incident faces
      CGAL_assertion(is_insert_successful.first->second.size() == 1);
      is_insert_successful.first->second.push_back(icf);
    }
  }

  void convert_to_1_cover()
  {
    if(is_1_cover())
      return;

    std::cout << "Transition to Phase 2 after " << number_of_vertices() << " vertices" << std::endl;
    std::ofstream out("res_dt2.off");
    CGAL::write_DT2_to_OFF(out, dt2);
    out.close();

    p2dt2.clear();
    p2dt2.tds().set_dimension(2);

    std::unordered_map<Vertex_handle, Vertex_handle> vertex_correspondence_map;

    typename DT2::Finite_vertices_iterator vit = dt2.finite_vertices_begin(),
                                           vend = dt2.finite_vertices_end();
    for(; vit!=vend; ++vit)
    {
      if(!is_canonical(vit))
        continue;

      Vertex_handle v = p2dt2.tds().create_vertex();
      v->set_point(vit->point());
      vertex_correspondence_map[vit] = v;
    }

    // @todo array instead of vector
    typedef std::map<std::set<Vertex_handle>, // two vertices of an edge
                     std::vector<std::pair<Face_handle, int> > > Incident_faces_map;
    Incident_faces_map incident_faces_map;

    size_type cfc = 0;
    typename DT2::Finite_faces_iterator fit = dt2.finite_faces_begin(), fend = dt2.finite_faces_end();
    for(; fit!=fend; ++fit)
    {
      if(!is_canonical(fit))
        continue;

      ++cfc;

      Vertex_handle p2dt2_v0 = vertex_correspondence_map.at(canonical_vertex(fit->vertex(0)));
      Vertex_handle p2dt2_v1 = vertex_correspondence_map.at(canonical_vertex(fit->vertex(1)));
      Vertex_handle p2dt2_v2 = vertex_correspondence_map.at(canonical_vertex(fit->vertex(2)));

      Face_handle f = p2dt2.tds().create_face(p2dt2_v0, p2dt2_v1, p2dt2_v2);

      Offset min_off = CGAL::min(fit->vertex(0)->offset(),
                                 fit->vertex(1)->offset(),
                                 fit->vertex(2)->offset());

      f->set_offsets(fit->vertex(0)->offset() - min_off,
                      fit->vertex(1)->offset() - min_off,
                      fit->vertex(2)->offset() - min_off);

      add_edge_to_incident_faces_map(f, 0, incident_faces_map);
      add_edge_to_incident_faces_map(f, 1, incident_faces_map);
      add_edge_to_incident_faces_map(f, 2, incident_faces_map);

      // Set up incident face information
      for(int i=0; i<3; ++i)
      {
        if(f->vertex(i)->face() == Face_handle())
          f->vertex(i)->set_face(f);
      }
    }

    // Set up adjacencies
    typename Incident_faces_map::const_iterator ifit = incident_faces_map.begin();
    for(; ifit!=incident_faces_map.end(); ++ifit)
    {
      const std::vector<std::pair<Face_handle, int> >& adjacent_faces = ifit->second;
      CGAL_assertion(adjacent_faces.size() == 2);

      Face_handle f0 = adjacent_faces[0].first;
      int i0 = adjacent_faces[0].second;
      Face_handle f1 = adjacent_faces[1].first;
      int i1 = adjacent_faces[1].second;

      p2dt2.tds().set_adjacency(f0, i0, f1, i1);
    }

    std::ofstream out_p2t2("res_p2t2.off");
    CGAL::write_PD2T2_to_OFF(out_p2t2, p2dt2);
    out_p2t2.close();

    CGAL_postcondition(p2dt2.is_valid(true));

    is_1_cover_ = true;
  }

  /// Duals
  Point construct_barycenter(Face_handle f) const
  {
    Vector v0 = Vector(CGAL::ORIGIN, f->vertex(0)->point());
    Vector v1 = Vector(CGAL::ORIGIN, f->vertex(1)->point());
    Vector v2 = Vector(CGAL::ORIGIN, f->vertex(2)->point());

    Vector bcv = gt_.construct_scaled_vector_2_object()(
                   gt_.construct_sum_of_vectors_2_object()(
                     v0, gt_.construct_sum_of_vectors_2_object()(v1,v2)), FT(1)/3);

    return Point(bcv.x(), bcv.y());
  }

  /// is_valid and IO
  bool is_valid(Face_handle f, const bool verbose = false) const
  {
    if(is_1_cover())
      return true; // @todo is_valid(f) doesn't exist in (D)T2
    else
      return p2dt2.is_valid(f, verbose);
  }

  bool is_valid(const bool verbose = false)
  {
    if(is_1_cover())
    {
      if(!faces_with_too_big_circumradius.empty())
      {
        if(verbose)
        {
          std::cerr << "Error: periodic triangulation in one cover "
                    << "but there are faces with too big circumradius" << std::endl;
        }

        return false;
      }

      if(!p2dt2.is_valid())
        return false;
    }
    else
    {
      return dt2.is_valid();
    }
  }

  template <typename GT, typename TDS>
  friend std::ostream& operator<<(std::ostream& os, const Delaunay_triangulation_on_lattice_2<GT, TDS>& tr)
  {
    if(tr.is_1_cover())
      CGAL::write_PD2T2_to_OFF(os, tr.p2dt2);
    else
      CGAL::write_DT2_to_OFF(os, tr.dt2);

    return os;
  }
};

} // namespace CGAL

#endif // CGAL_P2T2_DELAUNAY_TRIANGULATION_ON_LATTICE_2_H
