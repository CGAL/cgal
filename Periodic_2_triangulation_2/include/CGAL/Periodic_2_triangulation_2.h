// Copyright (c) 1997-2020 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Nico Kruithof <Nico@nghk.nl>
//                 Mael Rouxel-Labbé

#ifndef CGAL_PERIODIC_2_TRIANGULATION_2_H
#define CGAL_PERIODIC_2_TRIANGULATION_2_H

#include <CGAL/license/Periodic_2_triangulation_2.h>

#include <CGAL/Periodic_2_triangulation_face_base_2.h>
#include <CGAL/Periodic_2_triangulation_vertex_base_2.h>
#include <CGAL/Periodic_2_triangulation_iterators_2.h>
#include <CGAL/periodic_2_triangulation_2_io.h>

#include <CGAL/basic.h>
#include <CGAL/iterator.h>
#include <CGAL/Iterator_project.h>
#include <CGAL/function_objects.h>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_utils_2.h>
#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/utility.h>

#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_smallint.hpp>
#include <boost/random/variate_generator.hpp>

#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <list>
#include <map>
#include <utility>
#include <vector>

namespace CGAL {

/// Periodic triangulation class.
/// Its main functionality is:
/// - Insertion of points
/// - Deletion of points
/// - Point location
template <class Gt,
          class Tds = Triangulation_data_structure_2 <
                        Periodic_2_triangulation_vertex_base_2<Gt>,
                        Periodic_2_triangulation_face_base_2<Gt> > >
class Periodic_2_triangulation_2
    : public Triangulation_cw_ccw_2
{
  typedef Periodic_2_triangulation_2<Gt, Tds>                               Self;

public:
  // Public types of Periodic_2_triangulation_2
  /// The triangulation data structure type
  typedef Tds                                                               Triangulation_data_structure;
  /// The traits class
  typedef Gt                                                                Geom_traits;

  /// The periodic offset type
  typedef typename Gt::Periodic_2_offset_2                                  Offset;
  /// The iso rectangle type
  typedef typename Gt::Domain                                               Domain;
  /// Integer tuple to store the number of sheets in each direction of space.
  typedef std::array<int, 2>                                                Covering_sheets;

  /// The point type
  typedef typename Gt::Point_2                                              Point;
  /// The vector type
  typedef typename Gt::Segment_2                                            Segment;
  /// The segment type
  typedef typename Gt::Vector_2                                             Vector;
  /// The triangle type
  typedef typename Gt::Triangle_2                                           Triangle;

  /// Represents a point-offset pair. The point in the pair lies in the original domain.
  typedef std::pair<Point, Offset>                                          Periodic_point;
  /// A pair of periodic points representing a segment in the periodic domain.
  typedef std::array<std::pair<Point, Offset>, 2>                           Periodic_segment;
  /// A triple of periodic points representing a triangle in the periodic domain.
  typedef std::array<std::pair<Point, Offset>, 3>                           Periodic_triangle;

  /// The vertex type
  typedef typename Tds::Vertex                                              Vertex;
  /// The face type
  typedef typename Tds::Face                                                Face;
  /// The edge type
  typedef typename Tds::Edge                                                Edge;

  /// Size type (an unsigned integral type)
  typedef typename Tds::size_type                                           size_type;
  /// Difference type (a signed integral type)
  typedef typename Tds::difference_type                                     difference_type;

  /// Handle to a vertex
  typedef typename Tds::Vertex_handle                                       Vertex_handle;
  /// Handle to a face
  typedef typename Tds::Face_handle                                         Face_handle;

  /// Iterator over the faces
  typedef typename Tds::Face_iterator                                       Face_iterator;
  /// Iterator over the edges
  typedef typename Tds::Edge_iterator                                       Edge_iterator;
  /// Iterator over the vertices
  typedef typename Tds::Vertex_iterator                                     Vertex_iterator;
  /// Iterator over the vertices whose corresponding points lie in the
  /// original domain, i.e. for each set of periodic copies the
  /// Unique_vertex_iterator iterates over exactly one representative.
  typedef Periodic_2_triangulation_unique_vertex_iterator_2<Self>           Unique_vertex_iterator;

  /// \name For compatibility with the Triangulation_2 class

  typedef Face_iterator                                                     Finite_faces_iterator;
  typedef Edge_iterator                                                     Finite_edges_iterator;
  typedef Vertex_iterator                                                   Finite_vertices_iterator;
  typedef Face_iterator                                                     All_faces_iterator;

  /// Circulator over all faces incident to a vertex
  typedef typename Tds::Face_circulator                                     Face_circulator;
  /// Circulator over all edges incident to a vertex
  typedef typename Tds::Edge_circulator                                     Edge_circulator;
  /// Circulator over all vertices incident to a vertex
  typedef typename Tds::Vertex_circulator                                   Vertex_circulator;

  /// \name Periodic iterator types

  /// Iterator over all periodic triangles
  typedef Periodic_2_triangulation_triangle_iterator_2<Self>                Periodic_triangle_iterator;
  /// Iterator over all periodic segments
  typedef Periodic_2_triangulation_segment_iterator_2<Self>                 Periodic_segment_iterator;
  /// Iterator over all periodic points
  typedef Periodic_2_triangulation_point_iterator_2<Self>                   Periodic_point_iterator;

  // Auxiliary iterators for convenience
  // do not use default template argument to please VC++
  /// Functor that returns the point given a vertex
  typedef Project_point<Vertex>                                             Proj_point;

  /// \name STL types

  /// value_type similar to stl containers
  typedef Point                                                             value_type; // to have a back_inserter
  /// const_reference similar to stl containers
  typedef const value_type&                                                 const_reference;
  /// reference similar to stl containers
  typedef value_type&                                                       reference;

  /// Tag to distinguish regular triangulations from others;
  typedef Tag_false                                                         Weighted_tag;

  /// Tag to distinguish periodic triangulations from others
  typedef Tag_true                                                          Periodic_tag;

protected:
  // Protected types of Periodic_2_triangulation_2
  typedef typename Gt::Orientation_2                                        Orientation_2;
  typedef typename Gt::Compare_x_2                                          Compare_x;
  typedef typename Gt::Compare_y_2                                          Compare_y;

  typedef typename Gt::FT                                                   FT;

  typedef std::pair<Vertex_handle, Offset>                                  Virtual_vertex;
  typedef std::map<Vertex_handle, Virtual_vertex>                           Virtual_vertex_map;
  typedef typename Virtual_vertex_map::const_iterator                       Virtual_vertex_map_it;

  /// Vector contains virtual copies with offset off:
  /// virtual copy with offset off is stored at position: i=3*off[0]+off[1]-1
  typedef std::map<Vertex_handle, std::vector<Vertex_handle> >              Virtual_vertex_reverse_map;
  typedef typename Virtual_vertex_reverse_map::const_iterator               Virtual_vertex_reverse_map_it;

public:
  /// \name Enumeration types

  /// Type determining how to iterate over the stored simplices in the triangulation
  enum Iterator_type
  {
    STORED = 0,
    UNIQUE, // 1
    STORED_COVER_DOMAIN, // 2
    UNIQUE_COVER_DOMAIN // 3
  };//3

  /// Return type of a point location query
  enum Locate_type
  {
    /// The query point lies on a vertex
    VERTEX = 0,
    /// The query point lies on an edge
    EDGE,
    /// The query point lies on a face
    FACE,
    /// The query point lies outside the affine hull of the triangulation,
    /// which is the case when the triangulation is empty.
    EMPTY,
    OUTSIDE_CONVEX_HULL, // unused, for compatibility with Alpha_shape_2
    OUTSIDE_AFFINE_HULL  // unused, for compatibility with Alpha_shape_2
  };

protected:
  /// \name Functors

  /// Functor for symbolically perturbing points
  class Perturbation_order
  {
    const Self *t;

  public:
    // Perturbation_order, public interface
    Perturbation_order(const Self *tr) : t(tr) { }

    bool operator()(const Point *p, const Point *q) const
    {
      return t->compare_xy(*p, *q) == SMALLER;
    }
    bool operator()(const Periodic_point *p, const Periodic_point *q) const
    {
      return t->compare_xy(p->first, q->first, p->second, q->second) == SMALLER;
    }
  };

private:
  class Finder;

  void copy_multiple_covering(const Periodic_2_triangulation_2 & tr);

public:
  /// \name Constructors

  /// Constructor
  Periodic_2_triangulation_2(const Gt& gt) : _gt(gt), _tds(), _cover(make_array(1, 1)) { }
  Periodic_2_triangulation_2(const Domain& domain) : Periodic_2_triangulation_2(Gt(domain)) { }

  /// Copy constructor
  Periodic_2_triangulation_2(const Periodic_2_triangulation_2<Gt, Tds>& tr) { copy_triangulation(tr); }

  /// Assignment
  Periodic_2_triangulation_2 &operator=(const Periodic_2_triangulation_2 &tr)
  {
    copy_triangulation(tr);
    return *this;
  }

  /// Copy the triangulation
  void copy_triangulation(const Periodic_2_triangulation_2 &tr);

  /// Swap two triangulations
  void swap(Periodic_2_triangulation_2 &tr);

  /// Clear the triangulation
  void clear()
  {
    _tds.clear();
    _tds.set_dimension(-2);

    v_offsets.clear();
    virtual_vertices.clear();
    virtual_vertices_reverse.clear();

    _cover = make_array(1, 1);
  }

  /// Changes the domain. Note that this function calls clear(), i.e.,
  /// it erases the existing triangulation.
  virtual void set_domain(const Domain& domain)
  {
    clear();
    _gt.set_domain(domain);
  }

  /// \name Access functions

  /// Returns the geometric traits used for the predicates and constructions.
  const Geom_traits& geom_traits() const { return _gt; }

  /// Returns the datastructure storing the triangulation.
  const Triangulation_data_structure & tds() const { return _tds; }

  /// Returns the datastructure storing the triangulation.
  Triangulation_data_structure & tds() { return _tds; }

  /// Returns the domain of the 1-sheeted cover.
  const Domain & domain() const { return _gt.get_domain(); }

  /// Returns the number of copies of the 1-sheeted cover stored in each of the principal directions.
  Covering_sheets number_of_sheets() const { return _cover; }

  /// Returns the dimension of the triangulation.
  int dimension() const
  {
    return _tds.dimension() == 2 ? 2 : 0;
  }

  /// \name Number of simplices

  /// Returns whether the triangulation is empty.
  bool empty() const { return _tds.dimension() < 2; }

  /// Returns the number of vertices. Counts all vertices that are
  /// representatives of the same point in the 1-cover as one vertex.
  size_type number_of_vertices() const
  {
    if(is_1_cover())
      return _tds.number_of_vertices();
    else
      return _tds.number_of_vertices() / 9;
  }

  /// Returns the number of edges. Counts all edges that are
  /// representatives of the same segment in the 1-cover as one edge.
  size_type number_of_edges() const
  {
    if(is_1_cover())
      return _tds.number_of_edges();
    else
      return _tds.number_of_edges() / 9;
  }

  /// Returns the number of faces. Counts all faces that are
  /// representatives of the same triangle in the 1-cover as one face.
  size_type number_of_faces() const
  {
    if(is_1_cover())
      return _tds.number_of_faces();
    else
      return _tds.number_of_faces() / 9;
  }

  /// Returns the number of vertices stored in the datastructure.
  size_type number_of_stored_vertices() const
  {
    return _tds.number_of_vertices();
  }

  /// Returns the number of edges stored in the datastructure.
  size_type number_of_stored_edges() const
  {
    return _tds.number_of_edges();
  }

  /// Returns the number of faces stored in the datastructure.
  size_type number_of_stored_faces() const
  {
    return _tds.number_of_faces();
  }

  /// \name Undocumented functions, needed by the geometric iterators

  /// [Undoc] Combines two offsets, where the first offset is defined by the
  /// virtual vertex and the second by the face.
  Offset combine_offsets(const Offset& o_c, const Offset& o_t) const
  {
    Offset o_ct(_cover[0] * o_t.x(), _cover[1] * o_t.y());
    return o_c + o_ct;
  }

  /// [Undoc] Returns the offset of nb==h->neighbor(i) with respect to fh.
  /// Get the offset between the origins of the internal offset coordinate
  /// systems of two neighboring faces with respect from ch to nb.
  ///
  /// - Find two corresponding vertices from each face
  /// - Return the difference of their offsets.
  ///
  Offset get_neighbor_offset(Face_handle fh, int i, Face_handle nb, int j) const
  {
    // Redundance in the signature
    CGAL_triangulation_precondition(fh->neighbor(i) == nb);
    CGAL_triangulation_precondition(nb->neighbor(j) == fh);
    CGAL_triangulation_precondition(fh->vertex(cw(i)) == nb->vertex(ccw(j)));

    return int_to_off(nb->offset(ccw(j))) - int_to_off(fh->offset(cw(i)));
  }

  /// [Undoc] Returns the offset of nb==fh->neighbor(i) with respect to fh.
  /// Get the offset between the origins of the internal offset coordinate
  /// systems of two neighboring faces with respect from fh to nb.
  ///
  /// - Find two corresponding vertices from each face
  /// - Return the difference of their offsets.
  ///
  Offset get_neighbor_offset(Face_handle fh, int i) const
  {
    Face_handle nb = fh->neighbor(i);
    return get_neighbor_offset(fh, i, nb, nb->index(fh));
  }

  /// [Undoc] returns the combined offset of the vertex
  /// (if we are not on the 1-cover) and the offset defined by the face.
  Offset get_offset(Face_handle f, int i) const
  {
    if(is_1_cover())
      return int_to_off(f->offset(i));

    Virtual_vertex_map_it it = virtual_vertices.find(f->vertex(i));
    if(it != virtual_vertices.end())
      return combine_offsets(it->second.second, int_to_off(f->offset(i)));
    else
      return combine_offsets(Offset(), int_to_off(f->offset(i)));
  }

  /// [Undoc] Returns the offset of the vertex if we are not on the 1-cover.
  Offset get_offset(Vertex_handle vh) const
  {
    if(is_1_cover())
      return Offset();
    Virtual_vertex_map_it it = virtual_vertices.find(vh);
    if(it != virtual_vertices.end())
      return it->second.second;
    else
      return Offset();
  }

  /// Converts an offset to a bit pattern where bit1==offx and bit0==offy.
#ifndef CGAL_GENERIC_P2T2
  int off_to_int(const Offset& off) const
  {
    CGAL_triangulation_assertion(off.x() == 0 || off.x() == 1);
    CGAL_triangulation_assertion(off.y() == 0 || off.y() == 1);
    int i = ((off.x() & 1) << 1) + (off.y() & 1);
    return i;
  }

  /// Creates an offset from a bit pattern.
  Offset int_to_off(int i) const
  {
    return Offset((i >> 1) & 1, i & 1);
  }
#else // CGAL_P2T2_USE_COMPACT_OFFSET
  Offset off_to_int(const Offset& off) { return off; }
  Offset int_to_off(const Offset& off) const { return off; }
#endif // CGAL_P2T2_USE_COMPACT_OFFSET

  /// Periodic functions

  /// These functions give the pair (vertex, offset) that corresponds
  /// to the i-th vertex of face f. The vertex returned is not a virtual copy.
  void get_vertex(Face_handle f, int i, Vertex_handle &vh, Offset &off) const;

  /// These functions give the pair (vertex, offset) that corresponds
  /// to the i-th vertex of vertex vh. The vertex returned is not a virtual copy.
  void get_vertex(Vertex_handle vh_i, Vertex_handle &vh, Offset &off) const;

  /// Returns the face containing the three vertices defined by vh[0], vh1[1] and vh[2].
  inline Face_handle get_face(const Vertex_handle* vh) const;

#ifndef CGAL_GENERIC_P2T2
  /// Assigns the offsets to the vertices of the face f, and makes the offset minimal in each direction.
  void set_offsets(Face_handle f, int o0, int o1, int o2)
  {
    int off0[2] = { (o0 >> 1) & 1, (o0 & 1) };
    int off1[2] = { (o1 >> 1) & 1, (o1 & 1) };
    int off2[2] = { (o2 >> 1) & 1, (o2 & 1) };

    // Make sure that there is at least one zero offset in every direction
    for(int i=0; i<2; ++i)
    {
      int min_off = (std::min)((std::min)(off0[i], off1[i]), off2[i]);
      if(min_off != 0)
      {
        off0[i] -= min_off;
        off1[i] -= min_off;
        off2[i] -= min_off;
      }
    }
    o0 = ((off0[0] & 1) << 1) + (off0[1] & 1);
    o1 = ((off1[0] & 1) << 1) + (off1[1] & 1);
    o2 = ((off2[0] & 1) << 1) + (off2[1] & 1);
    f->set_offsets(o0, o1, o2);
  }
#endif // CGAL_P2T2_USE_COMPACT_OFFSET

  /// Assigns the offsets to the vertices of the face f, and makes the offset minimal in each direction.
  template<class Offset>
  void set_offsets(Face_handle f, const Offset& o0, const Offset& o1, const Offset& o2)
  {
    int off0[2] = { o0.x(), o0.y() };
    int off1[2] = { o1.x(), o1.y() };
    int off2[2] = { o2.x(), o2.y() };
    for(int i=0; i<2; ++i)
    {
      int min_off = (std::min)((std::min)(off0[i], off1[i]), off2[i]);
      if(min_off != 0)
      {
        off0[i] -= min_off;
        off1[i] -= min_off;
        off2[i] -= min_off;
      }
    }

    CGAL_triangulation_assertion((std::min)((std::min)(off0[0], off1[0]), off2[0]) == 0);
    CGAL_triangulation_assertion((std::min)((std::min)(off0[1], off1[1]), off2[1]) == 0);
    CGAL_triangulation_assertion((0 <= off0[0]) && (off0[0] < 2));
    CGAL_triangulation_assertion((0 <= off1[0]) && (off1[0] < 2));
    CGAL_triangulation_assertion((0 <= off2[0]) && (off2[0] < 2));
    CGAL_triangulation_assertion((0 <= off0[1]) && (off0[1] < 2));
    CGAL_triangulation_assertion((0 <= off1[1]) && (off1[1] < 2));
    CGAL_triangulation_assertion((0 <= off2[1]) && (off2[1] < 2));

#ifndef CGAL_GENERIC_P2T2
    int o0i = ((off0[0] & 1) << 1) + (off0[1] & 1);
    int o1i = ((off1[0] & 1) << 1) + (off1[1] & 1);
    int o2i = ((off2[0] & 1) << 1) + (off2[1] & 1);
    f->set_offsets(o0i, o1i, o2i);
#else // CGAL_P2T2_USE_COMPACT_OFFSET
    // @gp2t2 correct to take min above?
    f->set_offsets(Offset(off0[0], off0[1]), Offset(off1[0], off1[1]), Offset(off2[0], off2[1]));
#endif // CGAL_P2T2_USE_COMPACT_OFFSET
  }

  /// \name Geometric access functions

  /// Returns the periodic point given by vertex v. If t is
  /// represented in the 1-sheeted covering space, the offset is
  /// always zero. Otherwise v can correspond to a periodic copy
  /// outside domain of an input point.
  Periodic_point periodic_point(Vertex_handle v) const
  {
    return Periodic_point(v->point(), get_offset(v));
  }

  /// If t is represented in the 1-sheeted covering space, this
  /// function returns the periodic point given by the i-th vertex of
  /// face f, that is the point in the original domain and the offset
  /// of the vertex in f. If t is represented in the 9-sheeted
  /// covering space, this offset is possibly added to another offset
  /// determining the periodic copy.
  /// \pre i == {0,1,2}
  Periodic_point periodic_point(Face_handle f, int i) const
  {
    return Periodic_point(f->vertex(i)->point(), get_offset(f, i));
  }

  /// Returns the periodic segment formed by the two point-offset
  /// pairs corresponding to the two vertices of edge (f,i).
  /// \pre i == {0,1,2}
  Periodic_segment periodic_segment(Face_handle f, int i) const
  {
    CGAL_triangulation_precondition(number_of_vertices() != 0);
    CGAL_triangulation_precondition(i >= 0 && i <= 2);

    return make_array(periodic_point(f, ccw(i)),
                      periodic_point(f, cw(i)));
  }

  /// Same as the previous method for edge e.
  Periodic_segment periodic_segment(const Edge &e) const
  {
    return periodic_segment(e.first, e.second);
  }

  /// Returns the periodic triangle formed by the three point-offset
  /// pairs corresponding to the three vertices of face `f`.
  Periodic_triangle periodic_triangle(Face_handle f) const
  {
    return make_array(periodic_point(f, 0), periodic_point(f, 1), periodic_point(f, 2));
  }

  /// Converts the Periodic_point pp (point-offset pair) to the corresponding
  /// Point in \f$R^2\f$.
  Point point(const Periodic_point& pp) const { return construct_point(pp.first, pp.second); }
  Point point(Vertex_handle v) const { return point(periodic_point(v)); }
  Point point(Face_handle fh, int i) const { return point(periodic_point(fh, i)); }

  /// Converts the Periodic_segment ps to a Segment in \f$R^2\f$.
  Segment segment(const Periodic_segment &ps) const
  {
    return construct_segment(ps[0].first, ps[1].first, ps[0].second, ps[1].second);
  }
  /// Converts the Periodic_triangle pt to a Triagle in \f$R^2\f$.
  Triangle triangle(const Periodic_triangle &pt) const
  {
    return construct_triangle(pt[0].first, pt[1].first, pt[2].first,
                              pt[0].second, pt[1].second, pt[2].second);
  }

  /// Constructs the segment associated with the edge (f,i), respects the offset
  Segment segment(Face_handle f, int i) const { return segment(periodic_segment(f, i)); }

  /// Constructs the segment associated with the edge e, respects the offset
  Segment segment(const Edge& e) const { return segment(periodic_segment(e)); }

  /// Constructs the segment associated with the edge ec, respects the offset
  Segment segment(const Edge_circulator& ec) const
  {
    return segment(periodic_segment(ec->first, ec->second));
  }

  /// Constructs the segment associated with the edge ei, respects the offset
  Segment segment(const Edge_iterator& ei) const
  {
    return segment(periodic_segment(ei->first, ei->second));
  }

  /// Constructs the triangle associated with the face f, respects the offset
  Triangle triangle(Face_handle f) const { return triangle(periodic_triangle(f)); }

  /// \name Queries on simplices

  /// Returns false, no infinite simplices in the periodic triangulation
  template<class T>
  bool is_infinite(const T&, int = 0) const { return false; }

  size_type degree(Vertex_handle v) const { return _tds.degree(v); }

  bool is_edge(Vertex_handle va, Vertex_handle vb) const
  {
    return _tds.is_edge(va, vb);
  }
  bool is_edge(Vertex_handle va, Vertex_handle vb, Face_handle& fr, int&  i) const
  {
    return _tds.is_edge(va, vb, fr, i);
  }
  bool is_face(Vertex_handle v1, Vertex_handle v2, Vertex_handle v3) const
  {
    return _tds.is_face(v1, v2, v3);
  }
  bool is_face(Vertex_handle v1, Vertex_handle v2, Vertex_handle v3,
               Face_handle &fr) const
  {
    return _tds.is_face(v1, v2, v3, fr);
  }

  /// \name Traversal of the Triangulation

  /// Iterator over all stored vertices. Starts at an arbitrary vertex.
  /// Returns vertices_end() if t.number_of_vertices()=0.
  Vertex_iterator vertices_begin() const  { return _tds.vertices_begin(); }
  /// Past the end Vertex_iterator.
  Vertex_iterator vertices_end() const { return _tds.vertices_end(); }

  /// Iterator over all stored edges. Starts at an arbitrary edge.
  /// Returns edges_end() if t.number_of_vertices()=0.
  Edge_iterator edges_begin() const { return _tds.edges_begin(); }
  /// Past the end Edge_iterator.
  Edge_iterator edges_end() const {  return _tds.edges_end(); }

  /// Iterator over all stored faces. Starts at an arbitrary face.
  /// Returns faces_end() if t.number_of_vertices()=0.
  Face_iterator faces_begin() const { return _tds.faces_begin(); }
  /// Past the end Face_iterator.
  Face_iterator faces_end() const { return _tds.faces_end(); }

  /// Iterator over all stored vertices. Starts at an arbitrary vertex.
  /// Returns vertices_end() if t.number_of_vertices()=0.
  Vertex_iterator finite_vertices_begin() const { return _tds.vertices_begin(); }
  /// Past the end Vertex_iterator.
  Vertex_iterator finite_vertices_end() const { return _tds.vertices_end(); }

  /// Iterator over all stored edges. Starts at an arbitrary edge.
  /// Returns edges_end() if t.number_of_vertices()=0.
  Edge_iterator finite_edges_begin() const { return _tds.edges_begin(); }
  /// Past the end Edge_iterator.
  Edge_iterator finite_edges_end() const { return _tds.edges_end(); }

  /// Iterator over all stored faces. Starts at an arbitrary face.
  /// Returns faces_end() if t.number_of_vertices()=0.
  Face_iterator finite_faces_begin() const { return _tds.faces_begin(); }
  /// Past the end Face_iterator.
  Face_iterator finite_faces_end() const { return _tds.faces_end(); }

  /// Iterator over all stored vertices. Starts at an arbitrary vertex.
  /// Returns vertices_end() if t.number_of_vertices()=0.
  Vertex_iterator all_vertices_begin() const { return _tds.vertices_begin(); }
  /// Past the end Vertex_iterator.
  Vertex_iterator all_vertices_end() const { return _tds.vertices_end(); }

  /// Iterator over all stored edges. Starts at an arbitrary edge.
  /// Returns edges_end() if t.number_of_vertices()=0.
  Edge_iterator all_edges_begin() const { return _tds.edges_begin(); }
  /// Past the end Edge_iterator.
  Edge_iterator all_edges_end() const { return _tds.edges_end(); }

  /// Iterator over all stored faces. Starts at an arbitrary face.
  /// Returns faces_end() if t.number_of_vertices()=0.
  Face_iterator all_faces_begin() const { return _tds.faces_begin(); }
  /// Past the end Face_iterator.
  Face_iterator all_faces_end() const { return _tds.faces_end(); }

  /// begin iterator over the non-virtual vertices
  Unique_vertex_iterator unique_vertices_begin() const
  {
    return CGAL::filter_iterator(vertices_end(),
                                 Periodic_2_triangulation_2_internal::Domain_tester<Self>(this),
                                 vertices_begin());
  }
  /// past-the-end iterator over the non-virtual vertices
  Unique_vertex_iterator unique_vertices_end() const
  {
    return CGAL::filter_iterator(vertices_end(),
                                 Periodic_2_triangulation_2_internal::Domain_tester<Self>(this));
  }

  /// \name Geometric iterators

  /// Start iterator over the points
  Periodic_point_iterator periodic_points_begin(Iterator_type it = STORED) const
  {
    return Periodic_point_iterator(this, it);
  }
  /// Past-the-end iterator over the points
  Periodic_point_iterator periodic_points_end(Iterator_type it = STORED) const
  {
    return Periodic_point_iterator(this, 1, it);
  }

  /// Start iterator over the segments
  Periodic_segment_iterator periodic_segments_begin(Iterator_type it = STORED) const
  {
    return Periodic_segment_iterator(this, it);
  }
  /// Past-the-end iterator over the segments
  Periodic_segment_iterator periodic_segments_end(Iterator_type it = STORED) const
  {
    return Periodic_segment_iterator(this, 1, it);
  }

  /// Start iterator over the triangles
  Periodic_triangle_iterator periodic_triangles_begin(Iterator_type it = STORED) const
  {
    return Periodic_triangle_iterator(this, it);
  }
  /// Past-the-end iterator over the triangles
  Periodic_triangle_iterator periodic_triangles_end(Iterator_type it = STORED) const
  {
    return Periodic_triangle_iterator(this, 1, it);
  }

  /// \name Incident simplices
  Face_circulator incident_faces(Vertex_handle v, Face_handle f = Face_handle()) const
  {
    return _tds.incident_faces(v, f);
  }
  Edge_circulator incident_edges(Vertex_handle v, Face_handle f = Face_handle()) const
  {
    return _tds.incident_edges(v, f);
  }
  Vertex_circulator adjacent_vertices(Vertex_handle v, Face_handle f = Face_handle()) const
  {
    return _tds.incident_vertices(v, f);
  }

  /// \name Traversal between adjacent faces

  Vertex_handle mirror_vertex(Face_handle f, int i) const { return _tds.mirror_vertex(f, i); }
  int mirror_index(Face_handle f, int i) const { return _tds.mirror_index(f, i); }

  /// \name Predicates and Constructions

  /// Determines whether the point q lies strictly between the points p and r
  /// p,q and r are supposed to be collinear points
  bool collinear_between(const Point& p, const Point& q, const Point& r) const;

  /// Determines whether the point (q,o_q) lies strictly between the points (p,o_p) and (r,o_r)
  /// (q,o_q), (p,o_p) and (r,o_r) are supposed to be collinear points
  bool collinear_between(const Point& p, const Point& q, const Point& r,
                         const Offset& o_p, const Offset& o_q, const Offset& o_r) const;

  /// Compares the x-coordinates of p and q
  Comparison_result compare_x(const Point& p, const Point& q) const
  {
    return geom_traits().compare_x_2_object()(p, q);
  }
  /// Compares the x-coordinates of (p,o_p) and (q,o_q)
  Comparison_result compare_x(const Point& p, const Point& q,
                              const Offset& o_p, const Offset& o_q) const
  {
    return geom_traits().compare_x_2_object()(p, q, o_p, o_q);
  }

  /// Compares (p,o_p) and (q,o_q) lexicographically
  Comparison_result compare_xy(const Point& p, const Point& q) const
  {
    Comparison_result res = geom_traits().compare_x_2_object()(p, q);
    if(res == EQUAL)
      return geom_traits().compare_y_2_object()(p, q);

    return res;
  }
  /// Compares p and q lexicographically
  Comparison_result compare_xy(const Point& p, const Point& q,
                               const Offset& o_p, const Offset& o_q) const
  {
    Comparison_result res = geom_traits().compare_x_2_object()(p, q, o_p, o_q);
    if(res == EQUAL)
      return geom_traits().compare_y_2_object()(p, q, o_p, o_q);

    return res;
  }

  /// Compares the y-coordinates of p and q
  Comparison_result compare_y(const Point& p, const Point& q) const
  {
    return geom_traits().compare_y_2_object()(p, q);
  }

  /// Compares the y-coordinates of (p,o_p) and (q,o_q)
  Comparison_result compare_y(const Point& p, const Point& q,
                              const Offset& o_p, const Offset& o_q) const
  {
    return geom_traits().compare_y_2_object()(p, q, o_p, o_q);
  }

  /// Checks for equality of p and q
  bool xy_equal(const Point& p, const Point& q) const { return compare_xy(p, q) == EQUAL; }

  /// Returns the orientation of p,r,q
  Orientation orientation(const Point& p, const Point& q, const Point& r) const
  {
    return geom_traits().orientation_2_object()(p, q, r);
  }
  /// Returns the orientation of (p,o0), (q,o1), (r,o2)
  Orientation orientation(const Point& p, const Point& q, const Point& r,
                          const Offset& o_p, const Offset& o_q, const Offset& o_r) const
  {
    return geom_traits().orientation_2_object()(p, q, r, o_p, o_q, o_r);
  }

  /// Determines whether the point p lies on the (un-)bounded side of the triangle (p0,p1,p2)
  Bounded_side bounded_side(const Point& p0, const Point& p1, const Point& p2, const Point& p) const;
  /// Determines whether the point (p,o) lies on the (un-)bounded side of the triangle ((p0,o0),(p1,o1),(p2,o2))
  Bounded_side bounded_side(const Point& p0, const Point& p1, const Point& p2, const Point& p,
                            const Offset& o0, const Offset& o1, const Offset& o2, const Offset& o) const;

  Oriented_side oriented_side(Face_handle f, const Point& p, const Offset& o) const;

  /// Testing where the point (p,off) lies w.r.t. the face f
  Bounded_side side_of_face(const Point& p, const Offset& off, Face_handle f, Locate_type &lt, int& li) const;
  /// Testing where the point (p,off) lies w.r.t. the face f
  Bounded_side side_of_face(const Point& p, Face_handle f, Locate_type &lt, int& li) const
  {
    return side_of_face(p, Offset(), f, lt, li);
  }

  /// \name Wrapping the traits
  Point construct_point(const Point& p, const Offset& o) const
  {
    return geom_traits().construct_point_2_object()(p, o);
  }
  Point construct_point(const Periodic_point& pp) const
  {
    return construct_point(pp.first, pp.second);
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
  Segment construct_segment(const Point& p1, const Point& p2,
                            const Offset& o1, const Offset& o2) const
  {
    return geom_traits().construct_segment_2_object()(p1, p2, o1, o2);
  }
  Segment construct_segment(const Periodic_segment& seg) const
  {
    return construct_segment(seg[0].first, seg[1].first, seg[0].second, seg[1].second);
  }

  // @todo structural filtering
  /// \name Point location

  template<class Conflict_tester>
  Offset get_location_offset(const Conflict_tester& tester, Face_handle f, bool& found) const
  {
    CGAL_triangulation_precondition(number_of_vertices() != 0);

    if(f->has_zero_offsets() && tester(f, Offset()))
    {
      found = true;
      return Offset();
    }
    else
    {
#ifndef CGAL_GENERIC_P2T2
      // Main idea seems to just test all possibilities.
      int cumm_off = f->offset(0) | f->offset(1) | f->offset(2);
      for(int i=0; i<4; ++i)
      {
        // Check if for each bit position: bit of i <= bit of cumm_off
        if(((cumm_off | (~i)) & 3) == 3)
        {
          const Offset o = int_to_off(i);
          if(tester(f, o))
          {
            found = true;
            return o;
          }
        }
      }
#else // CGAL_P2T2_USE_COMPACT_OFFSET
      int min_off_x = f->offset(2).x(),
          max_off_x = f->offset(2).x(),
          min_off_y = f->offset(2).y(),
          max_off_y = f->offset(2).y();
      for(int i=0; i<2; ++i)
      {
        min_off_x = (std::min)(min_off_x, f->offset(i).x());
        max_off_x = (std::max)(max_off_x, f->offset(i).x());
        min_off_y = (std::min)(min_off_y, f->offset(i).y());
        max_off_y = (std::max)(max_off_y, f->offset(i).y());
      }

      if(min_off_x == max_off_x && min_off_y == max_off_y)
      {
        CGAL_triangulation_assertion(tester(f, Offset(min_off_x, min_off_y)));
        found = true;
        return Offset(min_off_x, min_off_y);
      }
      else
      {
        for(int i=min_off_x; i<=max_off_x; ++i)
        {
          for(int j=min_off_y; j<=max_off_y; ++j)
          {
            const Offset oij(i, j);
            if(tester(f, oij))
            {
              found = true;
              return oij;
            }
          }
        }
      }
#endif // CGAL_P2T2_USE_COMPACT_OFFSET
    }

    CGAL_triangulation_assertion(false);
    return Offset();
  }

  // In this overload, we know we must find the offset
  template<class Conflict_tester>
  Offset get_location_offset(const Conflict_tester& tester, Face_handle f, const Point& p) const
  {
    bool found = false;
    Offset o = get_location_offset(tester, f, found);
    CGAL_triangulation_assertion(found);
    return o;
  }

  /// Do a remembering heuristic walk to locate point (p,o)
  Face_handle march_locate_2D(Face_handle f, const Point& p, const Offset& o, Offset& lo,
                              Locate_type& lt, int& li) const;

  /// Checks whether the result of two point location queries are equivalent.
  bool compare_walks(const Point& p, Face_handle f1, Face_handle f2,
                     Locate_type& lt1, Locate_type& lt2, int li1, int li2) const;

  Face_handle locate(const Point& p, const Offset& o_p, Offset& lo,
                     Locate_type& lt, int& li,
                     Face_handle start = Face_handle()) const
  {
    if(dimension() <= 0)
    {
      lt = EMPTY;
      li = 4;
      return Face_handle();
    }

    // Triangulation is not empty
    if(start == Face_handle())
      start = faces_begin();

    return march_locate_2D(start, p, o_p, lo, lt, li);
  }

  /// Locates the simplex containing the point represented by p and o.
  ///
  /// The type of the simplex is stored in lt.
  /// The simplex containing the point is returned using lt and li.
  /// The Face_handle start is the start point of the heuristic walk.
  Face_handle locate(const Point& p, const Offset& o,
                     Locate_type& lt, int& li,
                     Face_handle start = Face_handle()) const
  {
    Offset unused_off_query;
    return locate(p, Offset(), unused_off_query, lt, li, start);
  }

  /// Wrapper function for locate if the offset is omitted.
  Face_handle locate(const Point& p,
                     Locate_type& lt,
                     int& li,
                     Face_handle start = Face_handle()) const
  {
    return locate(p, Offset(), lt, li, start);
  }

  /// Wrapper function for locate if only the requested point is given.
  Face_handle locate(const Point& p, Face_handle start = Face_handle()) const
  {
    Locate_type lt;
    int li;
    return locate(p, Offset(), lt, li, start);
  }

  /// \name Insertion
#ifndef CGAL_GENERIC_P2T2
  Vertex_handle create_initial_triangulation(const Point& p);
#endif

  /// Inserts a point in the triangulation
  /// \param p the point to be inserted
  /// \param start the start face for point location
  /// \return The new vertex handle or an existing Vertex_handle if p was inserted before
  Vertex_handle insert(const Point& p, Face_handle start = Face_handle());

  /// Inserts a point in the triangulation
  /// \pre The point has been located in the triangulation
  Vertex_handle insert(const Point& p, Locate_type lt, Face_handle loc, int li);

  /// Insert a point in the triangulation
  Vertex_handle push_back(const Point& p) { return insert(p); }

  /// Inserts p in the face f and sets the offsets of the newly created faces
  /// Insert periodic copies in all periodic copies of the domain
  Vertex_handle insert_in_face(const Point& p, Face_handle f);

  template < class Conflict_tester, class Point_hider, class CoverManager >
  Vertex_handle periodic_insert(const Point& p, const Offset& o,
                                Locate_type lt, Face_handle f,
                                const Conflict_tester& tester,
                                Point_hider& hider,
                                CoverManager& cover_manager,
                                Vertex_handle vh = Vertex_handle())
  {
#ifdef CGAL_DEBUG_P2T2
    std::cout << "Periodic insert, Offset: " << o << std::endl;
#endif

    CGAL_triangulation_precondition(number_of_vertices() != 0);

    CGAL_triangulation_assertion_code(Locate_type lt_assert; int i_assert;)
    CGAL_triangulation_assertion(side_of_face(tester.point(), o, f, lt_assert, i_assert) != ON_UNBOUNDED_SIDE);

    tester.set_offset(o);

    // Choose the periodic copy of tester.point() that is in conflict with c.
    bool found = false;
    Offset current_off = get_location_offset(tester, f, found);
#ifdef CGAL_DEBUG_P2T2
    std::cout << "get_location_offset: " << current_off << std::endl;
#endif

    CGAL_triangulation_assertion(side_of_face(tester.point(), combine_offsets(o, current_off),
                                              f, lt_assert, i_assert)
                                 != ON_UNBOUNDED_SIDE);

    // If the new point is not in conflict with its face, it is hidden.
    if(!found || !tester.test_initial_face(f, current_off))
    {
      hider.hide_point(f, p);
      return Vertex_handle();
    }

    // Ok, we really insert the point now.
    // First, find the conflict region.
    std::vector<Edge> edges;
    edges.reserve(16);
    std::vector<Face_handle> faces;
    faces.reserve(16); // @todo what should those values be

    find_conflicts(f, current_off, tester, make_triple(std::back_inserter(edges),
                                                       std::back_inserter(faces),
                                                       Emptyset_iterator()));

    // Reset the conflict flag on the boundary.
    //
    // @todo technically the TDS insert_in_hole should do that, but TDS_2 doesn't have tds_data usage
    // because T2 doesn't use it (yet?)
    for(typename std::vector<Edge>::iterator eit=edges.begin(); eit!=edges.end(); ++eit)
      eit->first->neighbor(eit->second)->tds_data().clear();

#ifdef CGAL_DEBUG_P2T2
    std::cout << faces.size() << " faces in conflict" << std::endl;
    for(const auto fh : faces)
    {
      std::cout << "fh in conflict: " << &*fh << std::endl;
      std::cout << point(fh, 0) << " 0 "
                << point(fh, 1) << " 0 "
                << point(fh, 2) << " 0 "
                << point(fh, 0) << " 0" << std::endl;
    }
#endif

    // The points that are hidden by the conflicting faces,
    // as they will be deleted during the insertion.
    hider.set_vertices(faces.begin(), faces.end());

    if(!is_1_cover())
      cover_manager.delete_unsatisfying_elements(faces.begin(), faces.end());

    // Compute the star and put it into the data structure.
    Vertex_handle v = tds().create_vertex();
    _tds.insert_in_hole(v, edges.front(), faces.begin(), faces.end()); // @todo same optimization as in 3D?
    v->set_point(p);

    // Store the new faces from the star in nbs.
    // @todo this could be done within the _insert_in_hole without losing any
    // time because each face is visited in any case.
    std::vector<Face_handle> nbs;
    Face_circulator fc = incident_faces(v), done(fc);
    do { nbs.push_back(fc); }
    while(++fc != done);

    // For all neighbors of the newly added vertex v: fetch their offsets from
    // the tester and reset them in the triangulation data structure.
    for(Face_handle fh : nbs)
    {
      Offset off[3];
      for(int i=0; i<3; ++i)
        off[i] = fh->vertex(i)->offset();

      set_offsets(fh, off[0], off[1], off[2]);
    }

    for(Vertex_handle cv : v_offsets)
      cv->clear_offset();
    v_offsets.clear();

    if(vh != Vertex_handle())
    {
  //    CGAL_triangulation_assertion(virtual_vertices.find(v) == virtual_vertices.end());
      virtual_vertices[v] = Virtual_vertex(vh, o);
      virtual_vertices_reverse[vh].push_back(v);
    }

    if(!is_1_cover())
      cover_manager.insert_unsatisfying_elements(v, nbs.begin(), nbs.end());

    // Store the hidden points in their new faces.
    hider.reinsert_vertices(v);
    return v;
  }

  // COMMON INSERTION for DELAUNAY and REGULAR TRIANGULATION
  template <class Conflict_tester, class Point_hider, class CoverManager>
  Vertex_handle insert_in_conflict(const Point& p,
                                   Locate_type lt,
                                   Face_handle f, int li,
                                   const Conflict_tester& tester,
                                   Point_hider& hider,
                                   CoverManager& cover_manager)
  {
#ifdef CGAL_DEBUG_P2T2
    std::cout << "Insert: " << p << std::endl;
#endif

#ifndef CGAL_GENERIC_P2T2
    // @todo enable this check for generic, with a kind of pt.is_canonical(p);
    CGAL_triangulation_assertion((domain().xmin() <= p.x()) && (p.x() < domain().xmax()));
    CGAL_triangulation_assertion((domain().ymin() <= p.y()) && (p.y() < domain().ymax()));

    if(number_of_vertices() == 0)
    {
      Vertex_handle vh = create_initial_triangulation(p);
      cover_manager.create_initial_triangulation();

      return vh;
    }
#endif

    if((lt == VERTEX) && (tester.compare_weight(f->vertex(li)->point(), p) == 0))
      return f->vertex(li);

    Vertex_handle vstart;
    if(!is_1_cover())
    {
      Virtual_vertex_map_it vvmit = virtual_vertices.find(f->vertex(0));
      if(vvmit == virtual_vertices.end())
        vstart = f->vertex(0);
      else
        vstart = vvmit->second.first;

      CGAL_triangulation_assertion(virtual_vertices.find(vstart) == virtual_vertices.end());
      CGAL_triangulation_assertion(virtual_vertices_reverse.find(vstart) != virtual_vertices_reverse.end());
    }

    CGAL_triangulation_assertion(number_of_vertices() != 0);
    CGAL_triangulation_expensive_assertion(is_valid());

    Vertex_handle vh = periodic_insert(p, Offset(), lt, f, tester, hider, cover_manager);

    if(is_1_cover())
      return vh;

    virtual_vertices_reverse[vh] = std::vector<Vertex_handle>();
    Offset lo;

    // insert copies
    for(int i=0; i<_cover[0]; ++i)
    {
      for(int j=0; j<_cover[1]; ++j)
      {
        if((i != 0) || (j != 0))
        {
          f = locate(p, Offset(i,j), lo, lt, li, Face_handle());
          periodic_insert(p, Offset(i,j), lt, f, tester, hider, cover_manager, vh);
        }
      }
    }

    CGAL_triangulation_assertion(is_valid(true));

    // We will never be there because generic P2T2 is always single cover
    // but 'convert_to_1_sheeted_covering' is hidden
#ifndef CGAL_GENERIC_P2T2
    // Fall back to 1-cover if the criterion that the largest circumradius is short-enough
    if(cover_manager.can_be_converted_to_1_sheet())
    {
      std::cout << "Can convert!!" << std::endl;
      std::cin.get();

      convert_to_1_sheeted_covering();
      CGAL_triangulation_expensive_assertion(is_valid());
    }
#endif

    return vh;
  }

  template <class Conflict_tester, class Point_hider, class CoverManager>
  Vertex_handle insert_in_conflict(const Point& p,
                                   Face_handle start,
                                   const Conflict_tester& tester,
                                   Point_hider& hider,
                                   CoverManager& cover_manager)
  {
    Locate_type lt;
    int li = 0;
    Offset lo;
    Face_handle f = locate(p, Offset(), lo, lt, li, start);
    return insert_in_conflict(p, lt, f, li, tester, hider, cover_manager);
  }

  template <class InputIterator, class Conflict_tester, class Point_hider, class CoverManager>
  std::vector<Vertex_handle>
  insert_in_conflict(InputIterator begin, InputIterator end,
                     Face_handle start,
                     Conflict_tester& tester,
                     Point_hider& hider,
                     CoverManager& cover_manager)
  {
    Vertex_handle new_vertex;
    std::vector<Vertex_handle> double_vertices;
    Locate_type lt = Locate_type();
    int li=0;

    CGAL_triangulation_assertion_code(Locate_type lta;)
    CGAL_triangulation_assertion_code(int ia = 0;)

    Face_handle hint;
    while(begin != end)
    {
      tester.set_point(*begin);
      Offset lo;
      hint = locate(*begin, Offset(), lo, lt, li, start);

      CGAL_triangulation_assertion_code(if(number_of_vertices() != 0) {);
        CGAL_triangulation_assertion(side_of_face(*begin, Offset(), hint, lta, ia) != ON_UNBOUNDED_SIDE);
        CGAL_triangulation_assertion(lta == lt);
        CGAL_triangulation_assertion(ia == li);
      CGAL_triangulation_assertion_code(});

      new_vertex = insert_in_conflict(*begin, lt, hint, li, tester, hider, cover_manager);
      if(lt == VERTEX)
        double_vertices.push_back(new_vertex);

      if(new_vertex != Vertex_handle())
        start = new_vertex->face();

      ++begin;
    }

    return double_vertices;
  }

  template <class Conflict_test,
            class OutputIteratorBoundaryEdges,
            class OutputIteratorFaces,
            class OutputIteratorInternalEdges>
  Triple<OutputIteratorBoundaryEdges,
         OutputIteratorFaces,
         OutputIteratorInternalEdges>
  find_conflicts(Face_handle d,
                 const Offset& current_off,
                 const Conflict_test& tester,
                 Triple<OutputIteratorBoundaryEdges,
                        OutputIteratorFaces,
                        OutputIteratorInternalEdges> it) const
  {
    CGAL_triangulation_precondition(number_of_vertices() != 0);
    CGAL_triangulation_precondition(tester(d, current_off));

//    for(auto vit=all_vertices_begin(); vit!=all_vertices_end(); ++vit)
//      vit->clear_offset();

    std::stack<std::pair<Face_handle, Offset> > face_stack;
    face_stack.emplace(d, current_off);
    d->tds_data().mark_in_conflict();
    *it.second++ = d;

    do
    {
      Face_handle f = face_stack.top().first;

      Offset current_off2 = face_stack.top().second;
      face_stack.pop();

      for(int i=0; i<3; ++i)
      {
        Face_handle test = f->neighbor(i);

#ifdef CGAL_DEBUG_P2T2
        std::cout << "neighbor: " << std::endl;
        std::cout << point(test, 0) << " 0 "
                  << point(test, 1) << " 0 "
                  << point(test, 2) << " 0 "
                  << point(test, 0) << " 0" << std::endl;
#endif

        if(test->tds_data().is_in_conflict())
        {
          if(f < test)
            *it.third++ = Edge(f, i); // Internal edge

#ifdef CGAL_DEBUG_P2T2
          std::cout << "already in conflict" << std::endl;
#endif
          continue; // test was already in conflict.
        }

        if(test->tds_data().is_clear())
        {
          Offset o_test = current_off2 + get_neighbor_offset(f, i);
#ifdef CGAL_DEBUG_P2T2
          std::cout << "o_test: " << current_off2 << " + " << get_neighbor_offset(f, i) << " = " << o_test << std::endl;
#endif
          if(tester(test, o_test))
          {
            if(f < test)
              *it.third++ = Edge(f, i); // Internal edge

#ifdef CGAL_DEBUG_P2T2
            std::cout << "In conflict too!" << std::endl;
#endif
            face_stack.emplace(test, o_test);
            test->tds_data().mark_in_conflict();
            *it.second++ = test;
            continue;
          }

#ifdef CGAL_DEBUG_P2T2
          std::cout << "Not in conflict" << std::endl;
#endif

          test->tds_data().mark_on_boundary(); // test is on the boundary.
        }

        *it.first++ = Edge(f, i);
        for(int j=0; j<3; ++j)
        {
          if(j == i)
            continue;

#ifdef CGAL_DEBUG_P2T2
          std::cout << "marking " << f->vertex(j)->point() << " (" << point(f, j)
                    << ") with offset " << int_to_off(f->offset(j)) << " - " << current_off2
                    << " = " << int_to_off(f->offset(j)) - current_off2 << std::endl;
#endif

          if(!f->vertex(j)->get_offset_flag())
          {
#ifdef CGAL_DEBUG_P2T2
            std::cout << "NEW MARK" << std::endl;
#endif
            f->vertex(j)->set_offset(int_to_off(f->offset(j)) - current_off2);
            v_offsets.push_back(f->vertex(j));
          }
        }
      }
    }
    while(!face_stack.empty());

    return it;
  }

  template <class Conflict_test,
            class OutputIteratorBoundaryEdges,
            class OutputIteratorFaces,
            class OutputIteratorInternalEdges>
  Triple<OutputIteratorBoundaryEdges,
         OutputIteratorFaces,
         OutputIteratorInternalEdges>
  find_conflicts(Face_handle f,
                 const Conflict_test& tester,
                 Triple<OutputIteratorBoundaryEdges,
                        OutputIteratorFaces,
                        OutputIteratorInternalEdges> it) const
  {
    bool b = false;
    Offset off = get_location_offset(tester, f, b);
    if(b)
      return find_conflicts(f, off, tester, it);

    CGAL_triangulation_assertion(false); // @tmp only true if we are in insert
    return it;
  }

private:
  /// Inserts (p,o) in the face f and sets the offsets of the newly created faces
  /// Doesn't insert periodic copies
  Vertex_handle insert_in_face(const Point& p, const Offset& o,
                               Face_handle f,
                               Vertex_handle vh);

  /// Remove a vertex without removing it's possible periodic copies.
  /// Helper functions
  void remove_degree_3_single_copy(Vertex_handle vh);

protected:
  /// Inserts a point with an offset in the triangulation
  /// \pre The point has been located in the triangulation
  Vertex_handle insert(const Point& p, const Offset& o, Locate_type lt,
                       Face_handle loc, int li, Vertex_handle vh);

public:
  /// \name Methods regarding the covering

  /// [Undoc] Returns the non-virtual copy of the vertex.
  Vertex_handle get_original_vertex(Vertex_handle vh) const
  {
    if(is_1_cover())
      return vh;

    Virtual_vertex_map_it it = virtual_vertices.find(vh);
    if(it != virtual_vertices.end())
      return it->second.first;
    else
      return vh;
  }

  /// Tests whether a vertex is a periodic copy of a vertex in the 3-cover.
  bool is_virtual(Vertex_handle v)
  {
    if(is_1_cover())
      return false;
    return (virtual_vertices.find(v) != virtual_vertices.end());
  }

  const std::vector<Vertex_handle>& periodic_copies(const Vertex_handle v) const
  {
    CGAL_triangulation_precondition(number_of_sheets() != make_array(1, 1));
    CGAL_triangulation_precondition(virtual_vertices.find(v) == virtual_vertices.end());
    CGAL_triangulation_assertion(virtual_vertices_reverse.find(v) != virtual_vertices_reverse.end());
    return virtual_vertices_reverse.find(v)->second;
  }

  /// \name Miscellaneous

  /// Checks if the face is valid.
  bool is_valid(Face_handle fh, bool verbose = false, int level = 0) const;
  /// Checks if the triangulation is valid.
  bool is_valid(bool verbose = false, int level = 0) const;

  bool well_oriented(Vertex_handle v) const
  {
    Face_circulator fc = incident_faces(v), done(fc);
    do
    {
      Orientation o;

      Vertex_handle v0 = fc->vertex(0);
      Vertex_handle v1 = fc->vertex(1);
      Vertex_handle v2 = fc->vertex(2);
      if(fc->has_zero_offsets())
      {
        o = orientation(v0->point(), v1->point(), v2->point());
      }
      else
      {
        Offset off0 = get_offset(fc, 0);
        Offset off1 = get_offset(fc, 1);
        Offset off2 = get_offset(fc, 2);
        o = orientation(v0->point(), v1->point(), v2->point(), off0, off1, off2);
      }

      if(o != COUNTERCLOCKWISE)
        return false;
    }
    while(++fc != done);

    return true;
  }

  /** \name Checking helpers */
  /// calls has_self_edges for every face of the triangulation
  bool has_self_edges() const
  {
    Face_iterator it;
    for(it = all_faces_begin(); it != all_faces_end(); ++it)
      if(has_self_edges(it))
        return true;

    return false;
  }

  bool has_self_edges(Face_handle fh) const
  {
    CGAL_triangulation_assertion((fh->vertex(0) != fh->vertex(1)) || (fh->offset(0) != fh->offset(1)));
    CGAL_triangulation_assertion((fh->vertex(0) != fh->vertex(2)) || (fh->offset(0) != fh->offset(2)));
    CGAL_triangulation_assertion((fh->vertex(1) != fh->vertex(2)) || (fh->offset(1) != fh->offset(2)));
    return ((fh->vertex(0) == fh->vertex(1)) ||
            (fh->vertex(0) == fh->vertex(2)) ||
            (fh->vertex(1) == fh->vertex(2)));
  }

#ifdef MACRO_THAT_DOESNT_EXIT_TO_MAKE_GP2T2_WORK
  /// Serialize the triangulation to an output stream
  std::ostream& save(std::ostream& os) const;

  /// Deserialize the triangulation from an input stream
  std::istream& load(std::istream& is);
#endif

  template<class Stream>
  Stream& draw_triangulation(Stream& os) const
  {
    Edge_iterator it = edges_begin();
    for(; it != edges_end(); ++it)
      os << segment(it);

    return os;
  }

protected:
  std::vector<Vertex_handle> insert_dummy_points();

public:
  /// [Undoc] Returns whether the stored triangulation covers a 1-cover.
  bool is_1_cover() const { return (_cover[0] == 1) && (_cover[1] == 1); }

  /// Checks whether the triangulation is a valid simplicial complex in the one cover.
  bool is_triangulation_in_1_sheet() const;

#ifndef CGAL_GENERIC_P2T2
  /// Convert a 9 sheeted cover (used for sparse triangulations) to a single sheeted cover.
  /// \pre !is_1_cover();
  void convert_to_1_sheeted_covering();

  /// Convert a single sheeted cover (used for dense triangulations) to a 9 sheeted cover.
  /// \pre is_1_cover();
  void convert_to_9_sheeted_covering();
#endif

  virtual void clear_covering_data() { CGAL_triangulation_assertion(false); }
  virtual void update_cover_data_after_converting_to_9_sheeted_covering() { CGAL_triangulation_assertion(false); }

  virtual ~Periodic_2_triangulation_2() { }

protected:
  inline void try_to_convert_to_one_cover()
  {
    if(!is_1_cover() && is_triangulation_in_1_sheet())
    {
      CGAL_triangulation_expensive_assertion(is_valid());
      this->convert_to_1_sheeted_covering();
      CGAL_triangulation_expensive_assertion(is_valid());
    }
  }

protected:
  /// \name Triangulation data members

  /// Geometric traits
  Gt _gt;
  /// Triangulation data structure
  Tds _tds;

private:
  // Private data of Periodic_2_triangulation_2
  /// \name Periodic members

  /// Determines if we currently compute in 3-cover or 1-cover.
  Covering_sheets _cover;

private:
  /// map of offsets for periodic copies of vertices
  Virtual_vertex_map virtual_vertices;
  /// map of a non-virtual vertex to its virtual copies
  Virtual_vertex_reverse_map virtual_vertices_reverse;

protected:
  /// v_offsets temporarily stores all the vertices on the border of a conflict region.
  mutable std::vector<Vertex_handle> v_offsets;
};

// Helping functions
template <class Gt, class Tds>
class Periodic_2_triangulation_2<Gt, Tds>::Finder
{
  const Self* _t;
  const Point&  _p;
public:
  Finder(const Self* t, const Point& p) : _t(t), _p(p) { }

  bool operator()(const Vertex_handle v) { return _t->xy_equal(v->point(), _p); }
};

template <class Gt, class Tds>
inline void
Periodic_2_triangulation_2<Gt, Tds>::
copy_multiple_covering(const Periodic_2_triangulation_2<Gt, Tds>& tr)
{
  // Write the respective offsets in the vertices to make them
  // automatically copy with the tds.
  for(Vertex_iterator vit = tr.vertices_begin(); vit != tr.vertices_end(); ++vit)
    vit->set_offset(tr.get_offset(vit));

  // copy the tds
  _tds = tr.tds();

  // make a list of all vertices that belong to the original
  // domain and initialize the basic structure of
  // virtual_vertices_reverse
  std::list<Vertex_handle> vlist;
  for(Vertex_iterator vit = vertices_begin(); vit != vertices_end(); ++vit)
  {
    if(vit->offset() == Offset())
    {
      vlist.push_back(vit);
      virtual_vertices_reverse.insert(std::make_pair(vit, std::vector<Vertex_handle>(8)));
      CGAL_triangulation_assertion(virtual_vertices_reverse.find(vit)->second.size() == 8);
    }
  }

  // Iterate over all vertices that are not in the original domain
  // and construct the respective entries to virtual_vertices and
  // virtual_vertices_reverse
  for(Vertex_iterator vit2 = vertices_begin(); vit2 != vertices_end(); ++vit2)
  {
    if(vit2->offset() != Offset())
    {
      //TODO: use some binding, maybe boost instead of the Finder.
      typename std::list<Vertex_handle>::iterator vlist_it
          = std::find_if(vlist.begin(), vlist.end(), Finder(this, vit2->point()));
      Offset off = vit2->offset();
      virtual_vertices.insert(std::make_pair(vit2, std::make_pair(*vlist_it, off)));
      virtual_vertices_reverse.find(*vlist_it)->second[3 * off[0] + off[1] - 1] = vit2;
      CGAL_triangulation_assertion(get_offset(vit2) == off);
    }
  }

  // Cleanup vertex offsets
  for(Vertex_iterator vit = vertices_begin(); vit != vertices_end(); ++vit)
    vit->clear_offset();
  for(Vertex_iterator vit = tr.vertices_begin(); vit != tr.vertices_end(); ++vit)
    vit->clear_offset();
}

template<class Gt, class Tds>
void
Periodic_2_triangulation_2<Gt, Tds>::
copy_triangulation(const Periodic_2_triangulation_2 &tr)
{
  _tds.clear();
  _gt = tr.geom_traits();
  _cover = tr._cover;

  if(tr.is_1_cover())
    _tds = tr.tds();
  else
    copy_multiple_covering(tr);

  CGAL_triangulation_expensive_postcondition(*this == tr);
}

template<class Gt, class Tds>
void
Periodic_2_triangulation_2<Gt, Tds>::
swap(Periodic_2_triangulation_2 &tr)
{
  _tds.swap(tr._tds);

  Geom_traits t = geom_traits();
  _gt = tr.geom_traits();
  tr._gt = t;

  std::swap(tr._cover, _cover);
  std::swap(tr.virtual_vertices, virtual_vertices);
  std::swap(tr.virtual_vertices_reverse, virtual_vertices_reverse);
}

/**
 * - fh->offset(i) is an bit tuple encapsulated in an integer. Each bit
 *   represents the offset in one direction --> 2-cover!
 * - int_to_off(int) decodes this again.
 * - Finally the offset vector is multiplied by cover.
 *   So if we are working in 3-cover we translate it to the neighboring
 *   3-cover and not only to the neighboring domain.
 */
template<class Gt, class Tds>
inline void
Periodic_2_triangulation_2<Gt, Tds>::
get_vertex(Face_handle fh, int i, Vertex_handle &vh, Offset &off) const
{
  off = combine_offsets(Offset(), int_to_off(fh->offset(i)));
  vh = fh->vertex(i);

  if(is_1_cover())
    return;

  Vertex_handle vh_i = vh;
  get_vertex(vh_i, vh, off);
}

template<class Gt, class Tds>
inline void
Periodic_2_triangulation_2<Gt, Tds>::
get_vertex(Vertex_handle vh_i, Vertex_handle &vh, Offset &off) const
{
  Virtual_vertex_map_it it = virtual_vertices.find(vh_i);

  if(it == virtual_vertices.end())
  {
    // if vh_i is not contained in virtual_vertices, then it is in the
    // original domain.
    vh = vh_i;
    CGAL_triangulation_assertion(vh != Vertex_handle());
  }
  else
  {
    // otherwise it has to be looked up as well as its offset.
    vh = it->second.first;
    off += it->second.second;
  }
}

/** Find the Face that consists of the three given vertices
 *
 *  Iterates over all faces and compare the three vertices of each face
 *  with the three vertices in vh.
 */
template<class Gt, class Tds>
inline typename Periodic_2_triangulation_2<Gt, Tds>::Face_handle
Periodic_2_triangulation_2<Gt, Tds>::
get_face(const Vertex_handle* vh) const
{
  bool contains_v[2];
  Face_circulator fc = incident_faces(vh[2]);
  Face_circulator done(fc);
  do
  {
    CGAL_triangulation_assertion(fc->vertex(0) == vh[2] || fc->vertex(1) == vh[2] || fc->vertex(2) == vh[2]);
    for(int j=0; j<2; ++j)
      contains_v[j] = (fc->vertex(0) == vh[j]) || (fc->vertex(1) == vh[j]) || (fc->vertex(2) == vh[j]);

    if(contains_v[0] && contains_v[1])
      return fc;
  }
  while(++fc != done);

  CGAL_triangulation_assertion(false);
  return Face_handle();
}

template<class Gt, class Tds>
Bounded_side
Periodic_2_triangulation_2<Gt, Tds>::
side_of_face(const Point& q, const Offset& off, Face_handle f, Locate_type& lt, int& li) const
{
  CGAL_triangulation_precondition(number_of_vertices() != 0);

  Orientation o0, o1, o2;
  o0 = o1 = o2 = ZERO;

  if(f->has_zero_offsets() && is_1_cover())
  {
    CGAL_triangulation_assertion(off == Offset());

    const Point& p0 = f->vertex(0)->point();
    const Point& p1 = f->vertex(1)->point();
    const Point& p2 = f->vertex(2)->point();

    if(((o0 = orientation(q, p1, p2)) == RIGHT_TURN) ||
       ((o1 = orientation(p0, q, p2)) == RIGHT_TURN) ||
       ((o2 = orientation(p0, p1, q)) == RIGHT_TURN))
    {
      return ON_UNBOUNDED_SIDE;
    }
  }
  else // Special case for the periodic space.
  {
    Offset off_q;
    Offset offs[3];
    const Point *p[3];
    for(int i=0; i<3; ++i)
    {
      p[i] = &(f->vertex(i)->point());
      offs[i] = get_offset(f, i);
    }

    CGAL_triangulation_assertion(orientation(*p[0], *p[1], *p[2],
                                             offs[0], offs[1], offs[2]) == POSITIVE);

    bool found = false;

#ifndef CGAL_GENERIC_P2T2
    int cumm_off = f->offset(0) | f->offset(1) | f->offset(2);
    for(int i = 0; (i < 4) && (!found); ++i)
    {
      if((cumm_off | ((~i) & 3)) == 3)
      {
        o0 = o1 = o2 = RIGHT_TURN;
        off_q = combine_offsets(off, int_to_off(i));

        if(((o0 = orientation(q, *p[1], *p[2], off_q, offs[1], offs[2])) != RIGHT_TURN) &&
           ((o1 = orientation(*p[0], q, *p[2], offs[0], off_q, offs[2])) != RIGHT_TURN) &&
           ((o2 = orientation(*p[0], *p[1], q, offs[0], offs[1], off_q)) != RIGHT_TURN))
        {
          found = true;
        }
      }
    }
#else // CGAL_P2T2_USE_COMPACT_OFFSET
    int min_off_x = f->offset(2).x(),
        max_off_x = f->offset(2).x(),
        min_off_y = f->offset(2).y(),
        max_off_y = f->offset(2).y();
    for(int i=0; i<2; ++i)
    {
      min_off_x = (std::min)(min_off_x, f->offset(i).x());
      max_off_x = (std::max)(max_off_x, f->offset(i).x());
      min_off_y = (std::min)(min_off_y, f->offset(i).y());
      max_off_y = (std::max)(max_off_y, f->offset(i).y());
    }

    for(int i=min_off_x; i<=max_off_x && (!found); ++i)
    {
      for(int j=min_off_y; j<=max_off_y && (!found); ++j)
      {
        o0 = o1 = o2 = RIGHT_TURN;
        off_q = combine_offsets(off, Offset(i, j));

        if(((o0 = orientation(q, *p[1], *p[2], off_q, offs[1], offs[2])) != RIGHT_TURN) &&
           ((o1 = orientation(*p[0], q, *p[2], offs[0], off_q, offs[2])) != RIGHT_TURN) &&
           ((o2 = orientation(*p[0], *p[1], q, offs[0], offs[1], off_q)) != RIGHT_TURN))
        {
          found = true;
          // we also quit the double loop now, so that the oi's are the correct ones
        }
      }
    }
#endif // CGAL_P2T2_USE_COMPACT_OFFSET

    if(!found)
      return ON_UNBOUNDED_SIDE;
  }

  // now all the oi's are >=0
  // sum gives the number of faces p lies on
  int sum = ((o0 == ZERO) ? 1 : 0) + ((o1 == ZERO) ? 1 : 0) + ((o2 == ZERO) ? 1 : 0);

  switch(sum)
  {
    case 0:
    {
      lt = FACE;
      li = 4;
      return ON_BOUNDED_SIDE;
    }
    case 1:
    {
      lt = EDGE;
      // i = index such that q lies on edge (f,li)
      li = (o0 == ZERO) ? 0 : (o1 == ZERO) ? 1 : 2;
      return ON_BOUNDARY;
    }
    case 2:
    {
      lt = VERTEX;
      // i = index such that q lies on vertex li
      li = (o0 != ZERO) ? 0 : (o1 != ZERO) ? 1 : 2;
      return ON_BOUNDARY;
    }
    default:
    {
      // impossible : cannot be on 3 edges for a real triangle
      CGAL_triangulation_assertion(false);
      return ON_BOUNDARY;
    }
  }
}

template<class Gt, class Tds>
Oriented_side
Periodic_2_triangulation_2<Gt, Tds>::
oriented_side(Face_handle f, const Point& p, const Offset& o) const
{
  Point& p0 = f->vertex(0)->point();
  Point& p1 = f->vertex(1)->point();
  Point& p2 = f->vertex(2)->point();

  if(f->has_zero_offsets() && is_1_cover())
  {
    CGAL_precondition(o == Offset());

    // return position of point p with respect to the oriented triangle p0p1p2
    // the orientation of the vertices is assumed to be counter clockwise
    CGAL_triangulation_assertion(orientation(p0, p1, p2) == LEFT_TURN);

    Bounded_side bs = bounded_side(p0, p1, p2, p);
    switch(bs)
    {
      case ON_BOUNDARY:
        return ON_ORIENTED_BOUNDARY;
      case ON_BOUNDED_SIDE:
        return ON_POSITIVE_SIDE;
      case ON_UNBOUNDED_SIDE:
        return ON_NEGATIVE_SIDE;
    }
  }
  else // Special case for the periodic space.
  {
    Offset off_q;
    Offset off0 = get_offset(f, 0);
    Offset off1 = get_offset(f, 1);
    Offset off2 = get_offset(f, 2);

    // return position of point p with respect to the oriented triangle p0p1p2
    // the orientation of the vertices is assumed to be counter clockwise
    CGAL_triangulation_assertion(orientation(p0, p1, p2, off0, off1, off2) == LEFT_TURN);

    Bounded_side bs;

#ifndef CGAL_GENERIC_P2T2
    int cumm_off = f->offset(0) | f->offset(1) | f->offset(2);
    for(int i=0; i<4; ++i)
    {
      if((cumm_off | ((~i) & 3)) == 3)
      {
        off_q = combine_offsets(o, int_to_off(i));
        bs = bounded_side(p0, p1, p2, p, off0, off1, off2, off_q);
        if(bs != ON_UNBOUNDED_SIDE)
        {
          return (bs == ON_BOUNDARY ? ON_ORIENTED_BOUNDARY : ON_POSITIVE_SIDE);
        }
      }
    }
#else // CGAL_P2T2_USE_COMPACT_OFFSET
    int min_off_x = f->offset(2).x(),
        max_off_x = f->offset(2).x(),
        min_off_y = f->offset(2).y(),
        max_off_y = f->offset(2).y();
    for(int i=0; i<2; ++i)
    {
      min_off_x = (std::min)(min_off_x, f->offset(i).x());
      max_off_x = (std::max)(max_off_x, f->offset(i).x());
      min_off_y = (std::min)(min_off_y, f->offset(i).y());
      max_off_y = (std::max)(max_off_y, f->offset(i).y());
    }

    for(int i=min_off_x; i<=max_off_x; ++i)
    {
      for(int j=min_off_y; j<=max_off_y; ++j)
      {
        off_q = combine_offsets(o, Offset(i, j));
        bs = bounded_side(p0, p1, p2, p, off0, off1, off2, off_q);
        if(bs != ON_UNBOUNDED_SIDE)
        {
          return (bs == ON_BOUNDARY ? ON_ORIENTED_BOUNDARY : ON_POSITIVE_SIDE);
        }
      }
    }
#endif // CGAL_P2T2_USE_COMPACT_OFFSET

    return ON_NEGATIVE_SIDE;
  }

  CGAL_triangulation_assertion(false);
  return ON_NEGATIVE_SIDE;
}

template<class Gt, class Tds>
bool
Periodic_2_triangulation_2<Gt, Tds>::
collinear_between(const Point& p, const Point& q, const Point& r) const
{
  // return true if point q is strictly between p and r
  // p,q and r are supposed to be collinear points
  Comparison_result c_pr = compare_x(p, r);
  Comparison_result c_pq;
  Comparison_result c_qr;
  if(c_pr == EQUAL)
  {
    c_pq = compare_y(p, q);
    c_qr = compare_y(q, r);
  }
  else
  {
    c_pq = compare_x(p, q);
    c_qr = compare_x(q, r);
  }

  return ((c_pq == SMALLER) && (c_qr == SMALLER)) ||
         ((c_pq == LARGER)  && (c_qr == LARGER));
}

template<class Gt, class Tds>
bool Periodic_2_triangulation_2<Gt, Tds>::
collinear_between(const Point& p, const Point& q, const Point& r,
                  const Offset& o_p, const Offset& o_q, const Offset& o_r) const
{
  // return true if point q is strictly between p and r
  // p,q and r are supposed to be collinear points
  Comparison_result c_pr = compare_x(p, r, o_p, o_r);
  Comparison_result c_pq;
  Comparison_result c_qr;
  if(c_pr == EQUAL)
  {
    c_pq = compare_y(p, q, o_p, o_q);
    c_qr = compare_y(q, r, o_q, o_r);
  }
  else
  {
    c_pq = compare_x(p, q, o_p, o_q);
    c_qr = compare_x(q, r, o_q, o_r);
  }
  return (((c_pq == SMALLER) && (c_qr == SMALLER)) ||
          ((c_pq == LARGER)  && (c_qr == LARGER)));
}

template<class Gt, class Tds>
Bounded_side
Periodic_2_triangulation_2<Gt, Tds>::
bounded_side(const Point& p0, const Point& p1, const Point& p2, const Point& p) const
{
  // return position of point p with respect to triangle p0p1p2
  CGAL_triangulation_precondition(orientation(p0, p1, p2) != COLLINEAR);

  Orientation o1 = orientation(p0, p1, p);
  Orientation o2 = orientation(p1, p2, p);
  Orientation o3 = orientation(p2, p0, p);

  if(o1 == COLLINEAR)
  {
    if(o2 == COLLINEAR || o3 == COLLINEAR)
      return ON_BOUNDARY;
    if(collinear_between(p0, p, p1))
      return ON_BOUNDARY;

    return ON_UNBOUNDED_SIDE;
  }

  if(o2 == COLLINEAR)
  {
    if(o3 == COLLINEAR)
      return ON_BOUNDARY;

    if(collinear_between(p1, p, p2))
      return ON_BOUNDARY;

    return ON_UNBOUNDED_SIDE;
  }

  if(o3 == COLLINEAR)
  {
    if(collinear_between(p2, p, p0))
      return ON_BOUNDARY;

    return ON_UNBOUNDED_SIDE;
  }

  // from here none ot, o1, o2 and o3 are known to be non null
  if(o1 == o2 && o2 == o3)
    return ON_BOUNDED_SIDE;

  return ON_UNBOUNDED_SIDE;
}

template<class Gt, class Tds>
Bounded_side Periodic_2_triangulation_2<Gt, Tds>::
bounded_side(const Point& p0, const Point& p1, const Point& p2, const Point& p,
             const Offset& o0, const Offset& o1, const Offset& o2, const Offset& o) const
{
  // return position of point p with respect to triangle p0p1p2
  CGAL_triangulation_precondition(orientation(p0, p1, p2, o0, o1, o2) != COLLINEAR);
  Orientation orient1 = orientation(p0, p1, p, o0, o1, o);
  Orientation orient2 = orientation(p1, p2, p, o1, o2, o);
  Orientation orient3 = orientation(p2, p0, p, o2, o0, o);

  if(orient1 == COLLINEAR)
  {
    if(orient2 == COLLINEAR || orient3 == COLLINEAR)
      return ON_BOUNDARY;

    if(collinear_between(p0, p, p1, o0, o, o1))
      return ON_BOUNDARY;

    return ON_UNBOUNDED_SIDE;
  }

  if(orient2 == COLLINEAR)
  {
    if(orient3 == COLLINEAR)
      return ON_BOUNDARY;

    if(collinear_between(p1, p, p2, o1, o, o2))
      return ON_BOUNDARY;

    return ON_UNBOUNDED_SIDE;
  }

  if(orient3 == COLLINEAR)
  {
    if(collinear_between(p2, p, p0, o2, o, o0))
      return ON_BOUNDARY;

    return ON_UNBOUNDED_SIDE;
  }

  // from here none ot, o1, o2 and o3 are known to be non null
  if(orient1 == orient2 && orient2 == orient3)
    return ON_BOUNDED_SIDE;

  return ON_UNBOUNDED_SIDE;
}

#ifndef CGAL_GENERIC_P2T2
template<class Gt, class Tds>
typename Periodic_2_triangulation_2<Gt, Tds>::Vertex_handle
Periodic_2_triangulation_2<Gt, Tds>::
create_initial_triangulation(const Point& p)
{
  CGAL_triangulation_assertion(empty());
  // The empty triangulation has a single sheeted cover
  _cover = make_array(3, 3);

  /// Virtual vertices, one per periodic domain
  Vertex_handle vir_vertices[3][3];

  /// Virtual faces, two per periodic domain
  Face_handle faces[3][3][2];

  // Initialise vertices:
  vir_vertices[0][0] = _tds.create_vertex();
  vir_vertices[0][0]->set_point(p);
  virtual_vertices_reverse[vir_vertices[0][0]] = std::vector<Vertex_handle>();
  for(int i=0; i<_cover[0]; ++i)
  {
    for(int j=0; j<_cover[1]; ++j)
    {
      if((i != 0) || (j != 0))
      {
        // Initialise virtual vertices out of the domain for debugging
        vir_vertices[i][j] = _tds.create_vertex();
        vir_vertices[i][j]->set_point(p); //+Offset(i,j));
        virtual_vertices[vir_vertices[i][j]] = Virtual_vertex(vir_vertices[0][0], Offset(i, j));
        virtual_vertices_reverse[vir_vertices[0][0]].push_back(vir_vertices[i][j]);
      }
    }
  }

  // Create faces:
  for(int i=0; i<_cover[0]; ++i)
    for(int j=0; j<_cover[1]; ++j)
      for(int f = 0; f < 2; f++)
        faces[i][j][f] = _tds.create_face(); // f faces per 'rectangle'

  // table containing the vertex information
  // index to the right vertex: [number of faces][vertex][offset]
  int vertex_ind[2][3][2] = { { { 0, 0 }, { 1, 1 }, { 0, 1 } },
                              { { 0, 0 }, { 1, 0 }, { 1, 1 } } };
  // Table containing the neighbor information
  // [number of faces][neighbor][offset,face]
  int neighb_ind[2][3][3] = { { { 0, 1, 1 }, { -1, 0, 1 }, { 0, 0, 1 } },
                              { { 1, 0, 0 }, { 0, 0, 0 }, { 0, -1, 0 } } };

  for(int i=0; i<_cover[0]; ++i)
  {
    for(int j=0; j<_cover[1]; ++j)
    {
      int offset = ((i == _cover[0] - 1 ? 2 : 0) | (j == _cover[1] - 1 ? 1 : 0));
      for(int f=0; f<2; ++f)
      {
        faces[i][j][f]->set_vertices(vir_vertices[(i + vertex_ind[f][0][0])
            % _cover[0]][(j + vertex_ind[f][0][1]) % _cover[1]],
            vir_vertices[(i + vertex_ind[f][1][0]) % _cover[0]][(
              j + vertex_ind[f][1][1]) % _cover[1]], vir_vertices[(
                i + vertex_ind[f][2][0]) % _cover[0]][(j + vertex_ind[f][2][1]) % _cover[1]]);

        set_offsets(faces[i][j][f], offset & (vertex_ind[f][0][0] * 2
            + vertex_ind[f][0][1] * 1), offset & (vertex_ind[f][1][0] * 2
            + vertex_ind[f][1][1] * 1), offset & (vertex_ind[f][2][0] * 2
            + vertex_ind[f][2][1] * 1));
        faces[i][j][f]->set_neighbors(faces[(i + _cover[0]
            + neighb_ind[f][0][0]) % _cover[0]][(j + _cover[1]
            + neighb_ind[f][0][1]) % _cover[1]][neighb_ind[f][0][2]], faces[(
              i + _cover[0] + neighb_ind[f][1][0]) % _cover[0]][(j + _cover[1]
            + neighb_ind[f][1][1]) % _cover[1]][neighb_ind[f][1][2]], faces[(
              i + _cover[0] + neighb_ind[f][2][0]) % _cover[0]][(j + _cover[1]
            + neighb_ind[f][2][1]) % _cover[1]][neighb_ind[f][2][2]]);
      }
    }
  }

  // set pointers from the vertices to incident faces.
  for(int i=0; i<_cover[0]; ++i)
    for(int j=0; j<_cover[1]; ++j)
      vir_vertices[i][j]->set_face(faces[i][j][0]);

  _tds.set_dimension(2);

  return vir_vertices[0][0];
}
#endif

template<class Gt, class Tds>
typename Periodic_2_triangulation_2<Gt, Tds>::Vertex_handle
Periodic_2_triangulation_2<Gt, Tds>::
insert_in_face(const Point& p, const Offset& o, Face_handle f, Vertex_handle vh)
{
  CGAL_triangulation_assertion(f != Face_handle());
  CGAL_triangulation_assertion(number_of_vertices() != 0);
  CGAL_triangulation_assertion((!is_1_cover()) || (o == Offset()));

  const bool simplicity_criterion = f->has_zero_offsets() && o.is_zero();

  Offset current_off;

  // Save the neighbors and the offsets
  Face_handle nb[3];
  int nb_index[3];
  Offset offsets[3];
  CGAL_triangulation_assertion_code(Vertex_handle vertices[3];)

  if(!simplicity_criterion)
  {
    // Choose the periodic copy of tester.point() that is inside f.
    current_off = get_location_offset(f, p, o); // @fixme

// @tmp restore after oriented side has been converted to generic P2T2
//    CGAL_triangulation_assertion(oriented_side(f, p, combine_offsets(o, current_off)) != ON_NEGATIVE_SIDE);

    for(int i=0; i<3; ++i)
    {
      nb[i] = f->neighbor(i);
      nb_index[i] = nb[i]->index(f);
      offsets[i] = f->offset(i);
      CGAL_triangulation_assertion_code(vertices[i] = f->vertex(i););
    }
  }

  // Insert the new vertex
  Vertex_handle v = _tds.insert_in_face(f);
  v->set_point(p);

  if(!simplicity_criterion)
  {
    // Update the offsets
    Offset v_offset = off_to_int(current_off);
    Offset new_offsets[3];
    for(int i=0; i<3; ++i)
    {
      Face_handle new_face = nb[i]->neighbor(nb_index[i]);
      int v_index = new_face->index(v);

      CGAL_triangulation_assertion(new_face->vertex(ccw(v_index)) == vertices[ccw(i)]);
      CGAL_triangulation_assertion(new_face->vertex(cw(v_index)) == vertices[cw(i)]);

      new_offsets[v_index] = v_offset;
      new_offsets[ccw(v_index)] = offsets[ccw(i)];
      new_offsets[cw(v_index)] = offsets[cw(i)];
      set_offsets(new_face, new_offsets[0], new_offsets[1], new_offsets[2]);
    }
  }

  if(!is_1_cover())
  {
    // update the book-keeping in case of a periodic copy
    if(vh != Vertex_handle())
    {
      virtual_vertices[v] = Virtual_vertex(vh, o);
      virtual_vertices_reverse[vh].push_back(v);
    }
  }

  return v;
}

template<class Gt, class Tds>
typename Periodic_2_triangulation_2<Gt, Tds>::Vertex_handle
Periodic_2_triangulation_2<Gt, Tds>::
insert(const Point& p, Face_handle start)
{
  if(number_of_stored_vertices() == 0)
    return insert_first(p);

  if(start == Face_handle())
    start = faces_begin();

  Locate_type lt;
  int li;
  Face_handle loc = locate(p, lt, li, start);

  if(start != Face_handle())
    CGAL_triangulation_assertion(start->vertex(0) != Vertex_handle());

  return insert(p, lt, loc, li);
}

template<class Gt, class Tds>
typename Periodic_2_triangulation_2<Gt, Tds>::Vertex_handle
Periodic_2_triangulation_2<Gt, Tds>::
insert(const Point& p, Locate_type lt, Face_handle loc, int li)
{
  // @tmp
//  if(number_of_stored_vertices() == 0)
//    return insert_first(p);

  // vstart is a vertex incident to the Face_handle start that will be used as
  // for creating a start point for the virtual vertices.
  // We use the virtual copies of a vertex incident to loc.
  Vertex_handle vstart;
  if(!is_1_cover())
  {
    Virtual_vertex_map_it vvmit = virtual_vertices.find(loc->vertex(0));
    if(vvmit == virtual_vertices.end())
      vstart = loc->vertex(0);
    else
      vstart = vvmit->second.first;

    // vstart should be non-virtual, but should have virtual copies
    CGAL_triangulation_assertion(virtual_vertices.find(vstart) == virtual_vertices.end());
    CGAL_triangulation_assertion(virtual_vertices_reverse.find(vstart) != virtual_vertices_reverse.end());
  }

  Vertex_handle vh = insert(p, Offset(), lt, loc, li, Vertex_handle());

  // Don't add periodic copies if we are on the 1-cover
  if(is_1_cover())
    return vh;

  // Don't continue if the point lies on a vertex as this will break the
  // start_vertices array below.
  if(lt == VERTEX)
    return vh;

  const std::vector<Vertex_handle>& start_vertices = virtual_vertices_reverse.find(vstart)->second;
  CGAL_triangulation_assertion(start_vertices.size() == size_t(number_of_sheets()[0] * number_of_sheets()[1] - 1));

  for(int i=0; i<number_of_sheets()[0]; ++i)
  {
    for(int j=0; j<number_of_sheets()[1]; ++j)
    {
      if((i != 0) || (j != 0))
      {
        loc = start_vertices[i * 3 + j - 1]->face();
        Offset offset(i, j);

        loc = locate(p, offset, lt, li, loc);

        insert(p, offset, lt, loc, li, vh);
      }
    }
  }

  return vh;
}

// insert a point p, whose localization is known (lt, f, i)
template<class Gt, class Tds>
typename Periodic_2_triangulation_2<Gt, Tds>::Vertex_handle
Periodic_2_triangulation_2<Gt, Tds>::
insert(const Point& p, const Offset& o, Locate_type lt, Face_handle loc, int li, Vertex_handle vh)
{
  Vertex_handle result;
  switch(lt)
  {
    case FACE:
    {
      result = insert_in_face(p, o, loc, vh);
      break;
    }
    case EDGE:
    {
      CGAL_triangulation_assertion(false);
      break;
    }
    case VERTEX:
    {
      // The vertex is a special case, we can return immediately
      CGAL_triangulation_assertion(vh == Vertex_handle());
      return loc->vertex(li);
    }
    case EMPTY:
    {
      result = insert_first(p);
      break;
    }
    default:
    {
      CGAL_triangulation_assertion(false); // locate step failed
      return Vertex_handle();
    }
  }

  if(!is_1_cover() && (vh == Vertex_handle()))
    virtual_vertices_reverse[result] = std::vector<Vertex_handle>();

  return result;
}

template<class Gt, class Tds>
bool Periodic_2_triangulation_2<Gt, Tds>::
compare_walks(const Point& p,
              Face_handle f1, Face_handle f2,
              Locate_type& lt1, Locate_type& lt2,
              int li1, int li2) const
{
  bool b = true;
  b &= (lt1 == lt2);
  if((lt1 == lt2) && (lt1 == VERTEX))
  {
    b &= (f1->vertex(li1) == f2->vertex(li2));
  }
  else if((lt1 == lt2) && (lt1 == EDGE))
  {
    b &= ((f1 == f2) || ((f1->neighbor(li1) == f2) && (f2->neighbor(li2) == f1)));
  }
  else if((lt1 == lt2) && (lt1 == EMPTY))
  {
    // Skip
  }
  else
  {
    b &= (lt1 == lt2);
    b &= (lt1 == FACE);
    b &= (f1 == f2);
  }

  if(!b)
  {
    std::cout << "from compare_walks " << std::endl;
    std::cout << "point " << p << std::endl;
    std::cout << "locate 1 " << &*f1 << "\t" << lt1 << "\t" << li1 << std::endl;
    std::cout << "locate 2 " << &*f2 << "\t" << lt2 << "\t" << li2 << std::endl;
    std::cout << std::endl;
    show_face(f1);
    std::cout << std::endl;
    show_face(f2);
    std::cout << std::endl;
  }

  CGAL_triangulation_assertion(b);
  return b;
}

template<class Gt, class Tds>
typename Periodic_2_triangulation_2<Gt, Tds>::Face_handle
Periodic_2_triangulation_2<Gt, Tds>::
march_locate_2D(Face_handle f, const Point& query, const Offset& o_p, Offset& off_query,
                Locate_type& lt, int& li) const
{
  CGAL_triangulation_assertion(!empty());

#ifdef CGAL_DEBUG_P2T2
  std::cout << "                   ###########" << std::endl;
  std::cout << "March locate: " << query << " op: " << o_p << " " << construct_point(query, o_p) << std::endl;
  std::cout << "p2t2: " << number_of_vertices() << " " << number_of_faces() << std::endl;
#endif

  off_query = o_p;

  // Random generator
  boost::rand48 rng;
  boost::uniform_smallint<> two(0, 1);
  boost::variate_generator<boost::rand48&, boost::uniform_smallint<> > coin(rng, two);

  // Give the point the best start-offset possible
#ifndef CGAL_GENERIC_P2T2 // @fixme is it ok to ignore?
  if(is_1_cover() && !f->has_zero_offsets())
  {
    int cumm_off = f->offset(0) | f->offset(1) | f->offset(2);
    if(((cumm_off & 2) == 2) && (FT(2) * query.x() < (domain().xmax() + domain().xmin())))
      off_query += Offset(1, 0);
    if(((cumm_off & 1) == 1) && (FT(2) * query.y() < (domain().ymax() + domain().ymin())))
      off_query += Offset(0, 1);
  }
#endif

#ifdef CGAL_DEBUG_P2T2
  std::cout << "Off query post initial opti: " << off_query << std::endl;
#endif

  Face_handle prev = Face_handle();
  int prev_index = 0;
  Offset off[3];
  Orientation o[3];
  for(;;)
  {
#ifdef CGAL_DEBUG_P2T2
    std::cout << "------" << std::endl;
    std::cout << "current face: " << f->vertex(0)->point() << " "
                                  << f->vertex(1)->point() << " "
                                  << f->vertex(2)->point() << std::endl;
    std::cout << "current offs: " << int_to_off(f->offset(0)) << " "
                                  << int_to_off(f->offset(1)) << " "
                                  << int_to_off(f->offset(2)) << std::endl;
    std::cout << "face points: " << point(f, 0) << " 0 "
                                 << point(f, 1) << " 0 "
                                 << point(f, 2) << " 0 "
                                 << point(f, 0) << " 0" << std::endl;
    std::cout << "prev_index: " << prev_index << std::endl;
    std::cout << "previously on the negative side of " << point(f, ccw(prev_index)) << " 0 " << point(f, cw(prev_index)) << " 0" << std::endl;
    std::cout << "off_query: " << off_query << std::endl;
#endif

    // Instead of testing its edges in a random order we do the following
    // until we find a neighbor to go further:
    // As we come from prev we do not have to check the edge leading to prev
    // Now we flip a coin in order to decide if we start checking the
    // edge before or the edge after the edge leading to prev
    int left_first = coin() % 2;
    std::cout << "left_first: " << left_first << std::endl;

    bool simplicity_criterion = f->has_zero_offsets() && off_query.is_null() && is_1_cover();
    std::cout << "simplicity_criterion: " << simplicity_criterion << std::endl;

    const Point *p[3] =
    {
      &f->vertex(0)->point(),
      &f->vertex(1)->point(),
      &f->vertex(2)->point()
    };

    // Get the offsets
    if(!simplicity_criterion)
    {
      if(!is_1_cover())
      {
        // Just fetch the vertices of c as points with offsets
        for(int i=0; i<3; ++i)
          off[i] = get_offset(f, i);
      }
      else
      {
        // We are on the one cover and on the boundary between domains
        // Hence, we need to check predicates with offsets
        for(int i=0; i<3; ++i)
          off[i] = int_to_off(f->offset(i));
      }
    }

    if(prev == Face_handle())
    {
      prev = f;
      // First step, also check the prev_index
      if(simplicity_criterion)
      {
        o[ccw(prev_index)] = orientation(*p[ccw(prev_index)], *p[cw(prev_index)], query);
      }
      else
      {
        o[ccw(prev_index)] = orientation(*p[ccw(prev_index)], *p[cw(prev_index)], query,
                                         off[ccw(prev_index)], off[cw(prev_index)], off_query);
      }

      std::cout << "First; right turn? " << std::boolalpha << (o[ccw(prev_index)] == RIGHT_TURN) << std::endl;
      std::cout << "side of " << point(f, ccw(prev_index)) << " 0 " << point(f, cw(prev_index)) << " 0" << std::endl;
      if(o[ccw(prev_index)] == RIGHT_TURN)
      {
        // This assignment is already done: prev = f
        f = f->neighbor(prev_index);
        int new_index = f->index(prev);
        if(!(simplicity_criterion && f->has_zero_offsets()))
          off_query = combine_offsets(off_query, get_neighbor_offset(prev, prev_index, f, new_index));

        prev_index = new_index;
        continue;
      }
    }
    else
    {
      o[ccw(prev_index)] = LEFT_TURN;
    }

    if(left_first)
    {
      if(simplicity_criterion)
      {
        std::cout << "simple orientation: " << *p[prev_index] << " " << *p[ccw(prev_index)] << " " << query << std::endl;
        o[prev_index] = orientation(*p[prev_index], *p[ccw(prev_index)], query);
      }
      else
      {
        o[prev_index] = orientation(*p[prev_index], *p[ccw(prev_index)], query,
                                    off[prev_index], off[ccw(prev_index)], off_query);
      }

      std::cout << "left first; right turn? " << std::boolalpha << (o[prev_index] == RIGHT_TURN) << std::endl;
      std::cout << "side of " << point(f, prev_index) << " 0 " << point(f, ccw(prev_index)) << " 0" << std::endl;
      if(o[prev_index] == RIGHT_TURN)
      {
        prev = f;
        f = f->neighbor(cw(prev_index));
        int new_index = f->index(prev);
        if(!(simplicity_criterion && f->has_zero_offsets()))
          off_query = combine_offsets(off_query, get_neighbor_offset(prev, cw(prev_index), f, new_index));

        prev_index = new_index;
        continue;
      }
    }

    {
      // Do right side
      if(simplicity_criterion)
      {
        o[cw(prev_index)] = orientation(*p[cw(prev_index)], *p[prev_index], query);
      }
      else
      {
        o[cw(prev_index)] = orientation(*p[cw(prev_index)], *p[prev_index], query,
                                        off[cw(prev_index)], off[prev_index], off_query);
      }

      std::cout << "right; right turn? " << std::boolalpha << (o[cw(prev_index)] == RIGHT_TURN) << std::endl;
      std::cout << "side of " << point(f, cw(prev_index)) << " 0 " << point(f, prev_index) << " 0" << std::endl;
      if(o[cw(prev_index)] == RIGHT_TURN)
      {
        prev = f;
        f = f->neighbor(ccw(prev_index));
        int new_index = f->index(prev);
        if(!(simplicity_criterion && f->has_zero_offsets()))
          off_query = combine_offsets(off_query, get_neighbor_offset(prev, ccw(prev_index), f, new_index));

        prev_index = new_index;
        continue;
      }
    }

    if(!left_first)
    {
      if(simplicity_criterion)
      {
        o[prev_index] = orientation(*p[prev_index], *p[ccw(prev_index)], query);
      }
      else
      {
        o[prev_index] = orientation(*p[prev_index], *p[ccw(prev_index)], query,
                                    off[prev_index], off[ccw(prev_index)], off_query);
      }

      std::cout << "left second; right turn? " << std::boolalpha << (o[prev_index] == RIGHT_TURN) << std::endl;
      std::cout << "side of " << point(f, ccw(prev_index)) << " 0 " << point(f, cw(prev_index)) << " 0" << std::endl;
      if(o[prev_index] == RIGHT_TURN)
      {
        prev = f;
        f = f->neighbor(cw(prev_index));
        int new_index = f->index(prev);
        if(!(simplicity_criterion && f->has_zero_offsets()))
          off_query = combine_offsets(off_query, get_neighbor_offset(prev, cw(prev_index), f, new_index));

        prev_index = new_index;
        continue;
      }
    }

#ifdef CGAL_DEBUG_P2T2
    std::cout << "within..." << std::endl;
#endif

    // now p is in c or on its boundary
    int sum = (o[0] == COLLINEAR) + (o[1] == COLLINEAR) + (o[2] == COLLINEAR);
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
        li = (o[0] == COLLINEAR) ? 2 : (o[1] == COLLINEAR) ? 0 : 1;
        break;
      }
      case 2:
      {
        lt = VERTEX;
        li = (o[0] != COLLINEAR) ? 2 : (o[1] != COLLINEAR) ? 0 : 1;
        break;
      }
    }
    return f;
  }
}

#ifndef CGAL_GENERIC_P2T2
/** Delete each redundant face and the not anymore needed data
 *  structures.
 *
 *  This function consists of four iterations over all faces and one
 *  iteration over all vertices:
 *  -# Face iteration: mark all faces that are to delete
 *  -# Face iteration: redirect neighbors of remaining faces
 *  -# Face iteration: redirect vertices of remaining faces
 *  -# Face iteration: delete all faces marked in the 1. iteration
 *  -# Vertex iteration: delete all vertices outside the original domain.
 */
template<class Gt, class Tds>
void
Periodic_2_triangulation_2<Gt, Tds>::
convert_to_1_sheeted_covering()
{
  // ###################################################################
  // ### First face iteration ##########################################
  // ###################################################################
  {
    if(is_1_cover())
      return;

    CGAL_triangulation_expensive_assertion(is_triangulation_in_1_sheet());

    bool to_delete, has_simplifiable_offset;
    Virtual_vertex_map_it vvmit;
    // First iteration over all faces: Mark the faces that are to delete.
    // Faces are to delete if they cannot be translated anymore in the
    // direction of one of the axes without yielding negative offsets.
    for(Face_iterator it = all_faces_begin(); it != all_faces_end(); ++it)
    {
      to_delete = false;
      // for all directions in 2D Space
      for(int j=0; j<2; ++j)
      {
        has_simplifiable_offset = true;
        // for all vertices of face it
        for(int i=0; i<3; ++i)
        {
          vvmit = virtual_vertices.find(it->vertex(i));
          if(vvmit == virtual_vertices.end())
          {
            // if it->vertex(i) lies inside the original domain:

            // the face cannot be moved any more because if we did, then
            // it->vertex(i) will get at least one negative offset.
            has_simplifiable_offset = false;
          }
          else
          {
            // if it->vertex(i) lies outside the original domain:

            // The face can certainly be deleted if the offset contains a 2
            to_delete |= (vvmit->second.second[j] == 2);
            // The face can be moved into one direction only if the offset of
            // all for vertices is >=1 for this direction. Since we already
            // tested for 2 it is sufficient to test here for 1.
            has_simplifiable_offset &= (vvmit->second.second[j] == 1);
          }
        }
        // if the offset can be simplified, i.e. the face can be moved, then
        // it can be deleted.
        if(has_simplifiable_offset)
          to_delete = true;
      }
      // Mark all faces that are to delete. They cannot be deleted yet,
      // because neighboring information still needs to be extracted.
      it->set_additional_flag(to_delete ? 1 : 0);
    }
  }

  // ###################################################################
  // ### Second face iteration #########################################
  // ###################################################################
  {
    Vertex_handle vert[3], nbv[3];
    Offset off[3];
    Face_handle nb, new_neighbor;
    std::vector<Triple<Face_handle, int, Face_handle> > new_neighbor_relations;

    // Second iteration over all faces: redirect neighbors where necessary
    for(Face_iterator it = all_faces_begin(); it != all_faces_end(); ++it)
    {
      // Skip all faces that are to delete.
      if(it->get_additional_flag() == 1)
        continue;

      // Redirect neighbors: Only neighbors that are marked by the
      // additional_flag have to be substituted by one of their periodic
      // copies. The unmarked neighbors stay the same.
      for(int i=0; i<3; ++i)
      {
        if(it->neighbor(i)->get_additional_flag() != 1)
          continue;

        nb = it->neighbor(i);

        for(int j=0; j<3; ++j)
        {
          off[j] = Offset();
          get_vertex(nb, j, vert[j], off[j]);
        }

        int x, y;
        x = (std::min)((std::min)(off[0][0], off[1][0]), off[2][0]);
        y = (std::min)((std::min)(off[0][1], off[1][1]), off[2][1]);

        // The vector from nb to the "original" periodic copy of nb, that is
        // the copy that will not be deleted.
        Offset difference_offset(x, y);
        CGAL_triangulation_assertion(!difference_offset.is_null());

        // We now have to find the "original" periodic copy of nb from
        // its vertices. Therefore, we first have to find the vertices.
        for(int j=0; j<3; ++j)
        {
          CGAL_triangulation_assertion((off[j] - difference_offset)[0] >= 0);
          CGAL_triangulation_assertion((off[j] - difference_offset)[1] >= 0);
          CGAL_triangulation_assertion((off[j] - difference_offset)[0] < 3);
          CGAL_triangulation_assertion((off[j] - difference_offset)[1] < 3);

          // find the Vertex_handles of the vertices of the "original"
          // periodic copy of nb. If the vertex is inside the original
          // domain, there is nothing to do
          if((off[j] - difference_offset).is_null())
          {
            nbv[j] = vert[j];
            // If the vertex is outside the original domain, we have to search
            // in virtual_vertices in the "wrong" direction. That means, we
            // cannot use virtual_vertices.find but have to use
            // virtual_vertices_reverse.
          }
          else
          {
            Offset nbo = off[j] - difference_offset;
            nbv[j] = virtual_vertices_reverse.find(vert[j]) ->second[nbo[0]
                * 3 + nbo[1] - 1];
          }
        }

        // Find the new neighbor by its 4 vertices
        new_neighbor = get_face(nbv);

        // Store the new neighbor relation. This cannot be applied yet because
        // it would disturb the functioning of get_face(...)
        new_neighbor_relations.push_back(make_triple(it, i, new_neighbor));
      }
    }

    // Apply the new neighbor relations now.
    for(unsigned int i=0; i< new_neighbor_relations.size(); ++i)
    {
      new_neighbor_relations[i].first->set_neighbor(
            new_neighbor_relations[i].second, new_neighbor_relations[i].third);
    }
  }

  // ###################################################################
  // ### Third face iteration ##########################################
  // ###################################################################
  {
    Vertex_handle vert[3];
    Offset off[3];
    // Third iteration over all faces: redirect vertices where necessary
    for(Face_iterator it = all_faces_begin(); it != all_faces_end(); ++it)
    {
      // Skip all faces that are marked to delete
      if(it->get_additional_flag() == 1)
        continue;
      // Find the corresponding vertices of it in the original domain
      // and set them as new vertices of it.
      for(int i=0; i<3; ++i)
      {
        off[i] = Offset();
        get_vertex(it, i, vert[i], off[i]);
        it->set_vertex(i, vert[i]);

        // redirect also the face pointer of the vertex.
        it->vertex(i)->set_face(it);
      }
      // Set the offsets.
      set_offsets(it, off[0], off[1], off[2]);
      CGAL_triangulation_assertion(int_to_off(it->offset(0)) == off[0]);
      CGAL_triangulation_assertion(int_to_off(it->offset(1)) == off[1]);
      CGAL_triangulation_assertion(int_to_off(it->offset(2)) == off[2]);
    }
  }

  // ###################################################################
  // ### Fourth face iteration #########################################
  // ###################################################################
  {
    // Delete the marked faces.
    std::vector<Face_handle> faces_to_delete;
    for(Face_iterator fit = all_faces_begin(); fit != all_faces_end(); ++fit)
    {
      if(fit->get_additional_flag() == 1)
        faces_to_delete.push_back(fit);
    }
    for(typename std::vector<Face_handle>::iterator it =
        faces_to_delete.begin(); it != faces_to_delete.end(); ++it)
    {
      _tds.delete_face(*it);
    }
  }

  // ###################################################################
  // ### Vertex iteration ##############################################
  // ###################################################################
  {
    // Delete all the vertices in virtual_vertices, that is all vertices
    // outside the original domain.
    std::vector<Vertex_handle> vertices_to_delete;
    for(Vertex_iterator vit = all_vertices_begin(); vit != all_vertices_end(); ++vit)
    {
      if(virtual_vertices.count(vit) != 0)
      {
        CGAL_triangulation_assertion(virtual_vertices.count(vit) == 1);
        vertices_to_delete.push_back(vit);
      }
    }

    for(typename std::vector<Vertex_handle>::iterator it =
        vertices_to_delete.begin(); it != vertices_to_delete.end(); ++it)
    {
      _tds.delete_vertex(*it);
    }
  }

  _cover = make_array(1, 1);
  virtual_vertices.clear();
  virtual_vertices_reverse.clear();

  clear_covering_data();

  CGAL_triangulation_assertion(is_1_cover());
}

template<class Gt, class Tds>
void
Periodic_2_triangulation_2<Gt, Tds>::
convert_to_9_sheeted_covering()
{
  if(_cover == make_array(3, 3))
    return;

  CGAL_triangulation_precondition(is_1_cover());

  // Create 9 copies of each vertex and write virtual_vertices and
  // virtual_vertices_reverse
  std::list<Vertex_handle> original_vertices;
  // try to use std::copy instead of the following loop.
  for(Vertex_iterator vit = vertices_begin(); vit != vertices_end(); ++vit)
    original_vertices.push_back(vit);

  for(typename std::list<Vertex_handle>::iterator vit = original_vertices.begin();
      vit != original_vertices.end(); ++vit)
  {
    Vertex_handle v_cp;
    std::vector<Vertex_handle> copies;
    for(int i=0; i<3; ++i)
    {
      for(int j=0; j<3; ++j)
      {
        if(i == 0 && j == 0)
          continue;
        v_cp = _tds.create_vertex(*vit);
        copies.push_back(v_cp);
        virtual_vertices.insert(std::make_pair(v_cp, std::make_pair(*vit,
                                                                     Offset(i, j))));
      }
    }

    virtual_vertices_reverse.insert(std::make_pair(*vit, copies));
  }

  // Create 9 copies of each face from the respective copies of the
  // vertices and write virtual_faces and virtual_faces_reverse.
  typedef std::map<Face_handle, std::pair<Face_handle, Offset> >  Virtual_face_map;
  typedef std::map<Face_handle, std::vector<Face_handle> >        Virtual_face_reverse_map;
  typedef typename Virtual_face_reverse_map::const_iterator       VCRMIT;

  Virtual_face_map virtual_faces;
  Virtual_face_reverse_map virtual_faces_reverse;

  std::list<Face_handle> original_faces;
  for(Face_iterator fit = faces_begin(); fit != faces_end(); ++fit)
    original_faces.push_back(fit);

  // Store vertex offsets in a separate data structure
  std::list<Offset> off_v;
  for(typename std::list<Vertex_handle>::iterator vit =
      original_vertices.begin(); vit != original_vertices.end(); ++vit)
  {
    Face_handle fff = (*vit)->face();
    int v_index = fff->index(*vit);
    off_v.push_back(int_to_off(fff->offset(v_index)));
  }

  // Store neighboring offsets in a separate data structure
  std::list<std::array<Offset, 3> > off_nb;
  for(typename std::list<Face_handle>::iterator fit = original_faces.begin();
      fit != original_faces.end(); ++fit)
  {
    std::array<Offset, 3> off_nb_f;
    for(int i=0; i<3; ++i)
    {
      Face_handle fff = *fit;
      Face_handle nnn = fff->neighbor(i);
      off_nb_f[i] = get_neighbor_offset(fff, i, nnn, nnn->index(fff));
    }
    off_nb.push_back(off_nb_f);
  }

  // Create copies of faces
  for(typename std::list<Face_handle>::iterator fit = original_faces.begin();
      fit != original_faces.end(); ++fit)
  {
    Face_handle c_cp;
    Vertex_handle v0, v1, v2;
    std::vector<Face_handle> copies;
    Virtual_vertex_reverse_map_it vvrmit[3];

    Offset vvoff[3];
    for(int i=0; i<3; ++i)
    {
      vvrmit[i] = virtual_vertices_reverse.find((*fit)->vertex(i));
      CGAL_triangulation_assertion(vvrmit[i] != virtual_vertices_reverse.end());
      vvoff[i] = int_to_off((*fit)->offset(i));
    }

    Vertex_handle vvh[3];
    for(int n = 0; n < 8; ++n) // iterate over faces
    {
      for(int i=0; i<3; ++i) // iterate over vertices of the face
      {
        // Decomposition of n into an offset (nx,ny):
        // nx = ((n+1)/3)%3, ny = (n+1)%3
        int o_i = ((n + 1) / 3 + vvoff[i].x() + 3) % 3;
        int o_j = ((n + 1) + vvoff[i].y() + 3) % 3;
        int n_c = 3 * o_i + o_j - 1;

        CGAL_triangulation_assertion(n_c >= -1);

        if(n_c == -1)
          vvh[i] = (*fit)->vertex(i);
        else
          vvh[i] = vvrmit[i]->second[n_c];
      }

      c_cp = _tds.create_face(vvh[0], vvh[1], vvh[2]);
      copies.push_back(c_cp);
    }

    virtual_faces_reverse.insert(std::make_pair(*fit, copies));
  }

  // Set new vertices of boundary faces of the original domain.
  for(typename std::list<Face_handle>::iterator fit = original_faces.begin();
      fit != original_faces.end(); ++fit)
  {
    for(int i=0; i<3; ++i)
    {
      Virtual_vertex_reverse_map_it vvrmit = virtual_vertices_reverse.find((*fit)->vertex(i));
      CGAL_triangulation_assertion(vvrmit != virtual_vertices_reverse.end());
      Offset vvoff = int_to_off((*fit)->offset(i));
      if(!vvoff.is_null())
      {
        int n_f = 3 * vvoff.x() + vvoff.y() - 1;
        CGAL_triangulation_assertion(n_f >= 0);
        CGAL_triangulation_assertion(static_cast<unsigned int>(n_f) < vvrmit->second.size());
        (*fit)->set_vertex(i, vvrmit->second[n_f]);
      }
    }
  }

  // Set neighboring relations of face copies
  typename std::list<std::array<Offset, 3> >::iterator oit = off_nb.begin();
  for(typename std::list<Face_handle>::iterator fit = original_faces.begin();
      fit != original_faces.end(); ++fit, ++oit)
  {
    CGAL_triangulation_assertion(oit != off_nb.end());
    VCRMIT c_cp = virtual_faces_reverse.find(*fit);
    CGAL_triangulation_assertion(c_cp != virtual_faces_reverse.end());
    for(int i=0; i<3; ++i)
    {
      Face_handle fit_nb = (*fit)->neighbor(i);
      VCRMIT c_cp_nb = virtual_faces_reverse.find(fit_nb);
      CGAL_triangulation_assertion(c_cp_nb != virtual_faces_reverse.end());
      Offset nboff = (*oit)[i];
      for(int n = 0; n < 8; ++n)
      {
        int n_nb;
        if(nboff.is_null())
        {
          n_nb = n;
        }
        else
        {
          int o_i = ((n + 1) / 3 - nboff.x() + 3) % 3;
          int o_j = (n + 1 - nboff.y() + 3) % 3;
          n_nb = 3 * o_i + o_j - 1;
        }

        if(n_nb == -1)
        {
          CGAL_triangulation_assertion(fit_nb->has_vertex(c_cp->second[n]->vertex(ccw(i))));
          CGAL_triangulation_assertion(fit_nb->has_vertex(c_cp->second[n]->vertex(cw(i))));
          c_cp->second[n]->set_neighbor(i, fit_nb);
        }
        else
        {
          CGAL_triangulation_assertion(n_nb >= 0);
          CGAL_triangulation_assertion(static_cast<unsigned int>(n_nb) <= c_cp_nb->second.size());
          CGAL_triangulation_assertion(c_cp_nb->second[n_nb]->has_vertex(c_cp->second[n]->vertex(ccw(i))));
          CGAL_triangulation_assertion(c_cp_nb->second[n_nb]->has_vertex(c_cp->second[n]->vertex(cw(i))));
          c_cp->second[n]->set_neighbor(i, c_cp_nb->second[n_nb]);
        }
      }
    }
  }

  // Set neighboring relations of original faces
  oit = off_nb.begin();
  for(typename std::list<Face_handle>::iterator fit = original_faces.begin();
      fit != original_faces.end(); ++fit, ++oit)
  {
    CGAL_triangulation_assertion(oit != off_nb.end());
    for(int i=0; i<3; ++i)
    {
      Offset nboff = (*oit)[i];
      if(!nboff.is_null())
      {
        Face_handle fit_nb = (*fit)->neighbor(i);
        VCRMIT c_cp_nb = virtual_faces_reverse.find(fit_nb);
        CGAL_triangulation_assertion(c_cp_nb != virtual_faces_reverse.end());
        int o_i = (3 - nboff.x()) % 3;
        int o_j = (3 - nboff.y()) % 3;
        int n_nb = 3 * o_i + o_j - 1;

        CGAL_triangulation_assertion(n_nb >= 0);
        CGAL_triangulation_assertion(static_cast<unsigned int>(n_nb) <= c_cp_nb->second.size());
        CGAL_triangulation_assertion(c_cp_nb->second[n_nb]->has_vertex((*fit)->vertex(ccw(i))));
        CGAL_triangulation_assertion(c_cp_nb->second[n_nb]->has_vertex((*fit)->vertex(cw(i))));

        (*fit)->set_neighbor(i, c_cp_nb->second[n_nb]);
      }
    }
  }

  // Set incident faces
  for(Face_iterator fit = faces_begin(); fit != faces_end(); ++fit)
    for(int i=0; i<3; ++i)
      fit->vertex(i)->set_face(fit);

  // Set offsets where necessary
  for(typename std::list<Face_handle>::iterator fit = original_faces.begin();
      fit != original_faces.end(); ++fit)
  {
    VCRMIT c_cp = virtual_faces_reverse.find(*fit);
    CGAL_triangulation_assertion(c_cp != virtual_faces_reverse.end());
    Offset off[3];

    for(int i=0; i<3; ++i)
      off[i] = int_to_off((*fit)->offset(i));

    if(off[0].is_null() && off[1].is_null() && off[2].is_null())
      continue;

    for(int n = 0; n < 8; ++n)
    {
      Offset off_cp[4];
      int o_i = ((n + 1) / 3) % 3;
      int o_j = (n + 1) % 3;
      if(o_i != 2 && o_j != 2)
        continue;

      for(int i=0; i<3; ++i)
      {
        off_cp[i] = Offset((o_i == 2) ? off[i].x() : 0, (o_j == 2) ? off[i].y() : 0);
        CGAL_triangulation_assertion(off_cp[i].x() == 0 || off_cp[i].x() == 1);
        CGAL_triangulation_assertion(off_cp[i].y() == 0 || off_cp[i].y() == 1);
      }

      set_offsets(c_cp->second[n], off_cp[0], off_cp[1], off_cp[2]);
    }
  }

  // Iterate over all original faces and reset offsets.
  for(typename std::list<Face_handle>::iterator fit = original_faces.begin();
      fit != original_faces.end(); ++fit)
  {
    //This statement does not seem to have any effect
    set_offsets(*fit, 0, 0, 0);
    CGAL_triangulation_assertion((*fit)->offset(0) == 0);
    CGAL_triangulation_assertion((*fit)->offset(1) == 0);
    CGAL_triangulation_assertion((*fit)->offset(2) == 0);
  }

  _cover = make_array(3, 3);

  update_cover_data_after_converting_to_9_sheeted_covering();

  CGAL_triangulation_expensive_assertion(is_valid());
  CGAL_triangulation_assertion(!is_1_cover());
}
#endif // CGAL_GENERIC_P2T2


template<class Gt, class Tds>
bool
Periodic_2_triangulation_2<Gt, Tds>::
is_valid(Face_handle fh, bool verbose, int /*level*/) const
{
  bool result = true;

  const Point *p[3];
  Offset off[3];

  int xmin, xmax, ymin, ymax;
  xmin = ymin = 3;
  xmax = ymax = 0;
  for(int i=0; i<3; ++i)
  {
    Offset o = get_offset(fh, i);
    xmin = (std::min)(xmin, o[0]);
    xmax = (std::max)(xmax, o[0]);
    ymin = (std::min)(ymin, o[1]);
    ymax = (std::max)(ymax, o[1]);

    p[i] = &(fh->vertex(i)->point());
    off[i] = get_offset(fh, i);
  }

#ifndef CGAL_GENERIC_P2T2
  // Should at most cross 1 border in each direction
  result &= (xmax - xmin <= 1);
  result &= (ymax - ymin <= 1);
  if(!result)
  {
    if(verbose)
    {
      std::cout << "Periodic_2_triangulation_2: impossible offsets:" << "\n"
                << *p[0] << " \t" << off[0] << "\n"
                << *p[1] << " \t" << off[1] << "\n"
                << *p[2] << " \t" << off[2] << std::endl;
      std::cout << construct_point(*p[0], off[0]) << " 0 "
                << construct_point(*p[1], off[1]) << " 0 "
                << construct_point(*p[2], off[2]) << " 0" << std::endl;

      std::cout << "min/max: " << xmin << "," << xmax << " " << ymin << "," << ymax << std::endl;
      for(int i=0; i<3; ++i)
      {
        Offset o = get_offset(fh, i);
        std::cout << "Offset: " << o << std::endl;
      }

      std::cout << std::endl;
    }
    CGAL_triangulation_assertion(false);
  }
#endif

  if(orientation(*p[0], *p[1], *p[2], off[0], off[1], off[2]) != POSITIVE)
  {
    if(verbose)
    {
      std::cout << "Periodic_2_triangulation_2: wrong orientation:" << "\n"
                << *p[0] << " \t" << off[0] << "\n"
                << *p[1] << " \t" << off[1] << "\n"
                << *p[2] << " \t" << off[2] << std::endl;
      std::cout << construct_point(*p[0], off[0]) << " 0 "
                << construct_point(*p[1], off[1]) << " 0 "
                << construct_point(*p[2], off[2]) << " 0" << std::endl;
    }

    result = false;
  }

  return result;
}

template<class Gt, class Tds>
bool
Periodic_2_triangulation_2<Gt, Tds>::
is_valid(bool verbose, int level) const
{
  bool result = _tds.is_valid(verbose, level);
  CGAL_triangulation_assertion(result);

  for(Face_iterator fit = faces_begin(); fit != faces_end(); ++fit)
  {
    result = is_valid(fit, verbose, level);
    CGAL_triangulation_assertion(result);
  }

  // Check for the right number of simplices
  int copies = number_of_sheets()[0] * number_of_sheets()[1];
  result &= (number_of_stored_vertices() == copies * number_of_vertices());
  result &= (number_of_stored_edges() == copies * number_of_edges());
  result &= (number_of_stored_faces() == copies * number_of_faces());
  CGAL_triangulation_assertion(result);

  // check number of euler characteristic. This cannot be done by the Tds
  // which does not know the genus
  result &= (number_of_stored_vertices() - number_of_stored_edges() + number_of_stored_faces() == 0);
  CGAL_triangulation_assertion(result);

  result &= !has_self_edges();
  CGAL_triangulation_assertion(result);

  // Edges should not be longer than 1 periodicity
  for(Face_iterator fit = faces_begin(); fit != faces_end(); ++fit)
    result &= is_valid(fit, verbose, level);

  CGAL_triangulation_assertion(result);

  result &= is_1_cover() == virtual_vertices.empty();
  result &= is_1_cover() == virtual_vertices_reverse.empty();
  result &= (virtual_vertices.size() ==
             (number_of_sheets()[0] * number_of_sheets()[1] - 1) * virtual_vertices_reverse.size());
  CGAL_triangulation_assertion(result);

  for(Virtual_vertex_map_it it = virtual_vertices.begin(); it != virtual_vertices.end(); ++it)
  {
    const Vertex_handle &copy = it->first;
    const Vertex_handle &orig = it->second.first;
    const Offset& off = it->second.second;
    size_t index = number_of_sheets()[0] * off[0] + off[1] - 1;
    Virtual_vertex_reverse_map_it rev_it = virtual_vertices_reverse.find(orig);
    if(rev_it != virtual_vertices_reverse.end())
    {
      if(index < rev_it->second.size())
        result &= (rev_it->second[index] == copy);
      else
        result &= false;
    }
    else
    {
      result &= false;
    }
  }
  CGAL_triangulation_assertion(result);

  for(Virtual_vertex_reverse_map_it it = virtual_vertices_reverse.begin(); it != virtual_vertices_reverse.end(); ++it)
  {
    const std::vector<Vertex_handle>& copies = it->second;
    result &= copies.size() == 8;
    for(size_t i=0; i< copies.size(); ++i)
    {
      Virtual_vertex_map_it copy_it = virtual_vertices.find(copies[i]);
      if(copy_it != virtual_vertices.end())
        result &= copy_it->second.first == it->first;
      else
        result &= false;
    }
  }

  return result;
}

#ifdef MACRO_THAT_DOESNT_EXIT_TO_MAKE_GP2T2_WORK
template<class Gt, class Tds>
std::ostream&
Periodic_2_triangulation_2<Gt, Tds>::
save(std::ostream& os) const
{
  // writes :
  // the number of vertices
  // the domain as four coordinates: xmin ymin ymax zmax
  // the current covering that guarantees the triangulation to be a
  //     simplicial complex
  // the non combinatorial information on vertices (points in case of 1-sheeted
  //     covering, point-offset pairs otherwise)
  //     ALL PERIODIC COPIES OF ONE VERTEX MUST BE STORED CONSECUTIVELY
  // the number of faces
  // the faces by the indices of their vertices in the preceding list
  // of vertices, plus the non combinatorial information on each face
  // the neighbors of each face by their index in the preceding list of faces

  // outputs dimension, domain and number of vertices
  Covering_sheets cover = number_of_sheets();
  size_type n = number_of_vertices();

  if(is_ascii(os))
  {
    os << domain() << std::endl
       << cover[0] << " " << cover[1] << std::endl
       << n*cover[0]*cover[1] << std::endl;
  }
  else
  {
    os << domain();
    write(os, cover[0]);
    write(os, cover[1]);
    write(os, n * cover[0]*cover[1]);
  }

  std::cout << "Line:" << __LINE__ << " cover[0]:" << cover[0] << " cover[1]:" << cover[1] << " n*c0*f1:" << (n * cover[0]*cover[1]) << std::endl;
  std::cout << "save, #Vertices: " << n << std::endl;

  if(n == 0)
    return os;

  // write the vertices
  Unique_hash_map<Vertex_handle, std::size_t > V;
  std::size_t i = 0;
  if(is_1_cover())
  {
    for(Vertex_iterator it = vertices_begin(); it != vertices_end(); ++it)
    {
      V[it] = i++;
      os << it->point();
      if(is_ascii(os))
        os << std::endl;
    }
  }
  else
  {
    Virtual_vertex_map_it vit, vvit;
    std::vector<Vertex_handle> vv;
    for(Vertex_iterator it = vertices_begin(); it != vertices_end(); ++it)
    {
      vit = virtual_vertices.find(it);
      if(vit != virtual_vertices.end())
        continue;

      V[it] = i++;
      if(is_ascii(os))
        os << it->point() << std::endl << Offset(0, 0) << std::endl;
      else
        os << it->point() << Offset(0, 0);

      CGAL_triangulation_assertion(virtual_vertices_reverse.find(it) != virtual_vertices_reverse.end());
      vv = virtual_vertices_reverse.find(it)->second;
      CGAL_triangulation_assertion(vv.size() == 8);

      for(std::size_t j=0; j<vv.size(); ++j)
      {
        vvit = virtual_vertices.find(vv[j]);
        CGAL_triangulation_assertion(vvit != virtual_vertices.end());
        V[vv[j]] = i++;
        if(is_ascii(os))
          os << vv[j]->point() << std::endl << vvit->second.second << std::endl;
        else
          os << vv[j]->point() << vvit->second.second;
      }
    }
  }
  CGAL_triangulation_postcondition(i == _cover[0]*_cover[1]*n);

  Unique_hash_map<Face_handle, std::size_t> F;
  int inum = 0;
  // asks the tds for the combinatorial information
  // vertices of the faces
  size_type m = _tds.number_of_faces();
  if(is_ascii(os))
    os << std::endl << m << std::endl;
  else
    write(os, m);

  std::cout << "save, #Faces: " << m << std::endl;

  for(Face_iterator ib = faces_begin(); ib != faces_end(); ++ib)
  {
    F[ib] = inum++;
    for(int j=0; j<3; ++j)
    {
      if(is_ascii(os))
        os << V[ib->vertex(j)] << " ";
      else
        write(os, V[ib->vertex(j)]);
    }
    os << *ib;

    if(is_ascii(os))
      os << "\n";
  }

  if(is_ascii(os))
    os << "\n";

  std::cout << "save, face check: " << inum << " == " << m << std::endl;
  CGAL_triangulation_assertion(m == (size_type)inum);

  // neighbor pointers of the  faces
  for(Face_iterator it = faces_begin(); it != faces_end(); ++it)
  {
    for(int j=0; j<3; ++j)
    {
      CGAL_triangulation_assertion(F.is_defined(it->neighbor(j)));
      if(is_ascii(os))
        os << F[it->neighbor(j)] << " ";
      else
        write(os, F[it->neighbor(j)]);
    }

    if(is_ascii(os))
      os << "\n";
  }

  // write offsets
  //for(unsigned int i=0; i<number_of_faces(); ++i) {
  for(Face_iterator it = faces_begin(); it != faces_end(); ++it)
  {
    //Face_handle ch = std::find(faces_begin(), faces_end(), i);
    Face_handle ch(it);
    for(int j=0; j<3; ++j)
    {
      if(is_ascii(os))
      {
        os << ch->offset(j);
        if(j == 3)
          os << std::endl;
        else
          os << ' ';
      }
      else
      {
        write(os, ch->offset(j));
      }
    }
  }

  // write the non combinatorial information on the faces
  // using the << operator of Face
  // works because the iterator of the tds traverses the faces in the
  // same order as the iterator of the triangulation
  if(number_of_vertices() != 0)
  {
    for(Face_iterator it = faces_begin(); it != faces_end(); ++it)
    {
      os << *it; // other information
      if(is_ascii(os))
        os << std::endl;
    }
  }

  return os;
}

template<class Gt, class Tds>
std::istream&
Periodic_2_triangulation_2<Gt, Tds>::
load(std::istream& is)
{
  // reads
  // the current covering that guarantees the triangulation to be a
  //     simplicial complex
  // the number of vertices
  // the non combinatorial information on vertices (points in case of 1-sheeted
  //     covering, point-offset pairs otherwise)
  //     ALL PERIODIC COPIES OF ONE VERTEX MUST BE STORED CONSECUTIVELY
  // the number of faces
  // the faces by the indices of their vertices in the preceding list
  // of vertices, plus the non combinatorial information on each face
  // the neighbors of each face by their index in the preceding list of face
  CGAL_triangulation_precondition(is.good());

  clear();

  Iso_rectangle domain(0, 0, 1, 1);
  int cx = 0, cy = 0;
  size_type n = 0;

  if(is_ascii(is))
  {
    is >> domain;
    is >> cx >> cy >> n;
  }
  else
  {
    is >> domain;
    read(is, cx);
    read(is, cy);
    read(is, n);
  }
  std::cout << "Line:" << __LINE__ << " cx:" << cx << " cy:" << cy << " n:" << n << std::endl;

  CGAL_triangulation_assertion((n / (cx * cy))*cx*cy == n);

  tds().set_dimension((n == 0 ? -2 : 2));
  set_domain(domain);
  _cover = make_array(cx, cy);

  if(n == 0)
    return is;

  std::map< std::size_t, Vertex_handle > V;

  if(cx == 1 && cy == 1)
  {
    Point p;
    for(std::size_t i=0; i< n; ++i)
    {
      V[i] = tds().create_vertex();
      is >> p;
      V[i]->set_point(p);
    }
  }
  else
  {
    Vertex_handle v, w;
    std::vector<Vertex_handle> vv;
    Offset off;
    Point p;
    for(std::size_t i=0; i< n; ++i)
    {
      v = tds().create_vertex();
      V[i] = v;
      is >> p >> off;
      V[i]->set_point(p);
      vv.clear();
      for(int j = 1; j < cx * cy; ++j)
      {
        i++;
        w = tds().create_vertex();
        V[i] = w;
        is >> p >> off;
        V[i]->set_point(p);
        vv.push_back(w);
        virtual_vertices[w] = std::make_pair(v, off);
      }
      virtual_vertices_reverse[v] = vv;
    }
  }

  // Creation of the faces
  std::size_t index;
  size_type m;
  if(is_ascii(is))
    is >> m;
  else
    read(is, m);

  std::vector<Face_handle> F(m);
  std::cout << "load, #Faces: " << m << std::endl;
  {
    for(size_t i=0; i< m; ++i)
    {
      F[i] = _tds.create_face();
      for(int j=0; j<3; ++j)
      {
        if(is_ascii(is))
          is >> index;
        else
          read(is, index);

        CGAL_triangulation_assertion(index < V.size());
        F[i]->set_vertex(j, V[index]);

        // The face pointer of vertices is set too often,
        // but otherwise we had to use one more map
        V[index]->set_face(F[i]);
      }

      // read in non combinatorial info of the face
      is >> *(F[i]);
    }
  }

  // Setting the neighbor pointers
  {
    for(size_t i=0; i< m; ++i)
    {
      for(int j=0; j<3; ++j)
      {
        if(is_ascii(is))
          is >> index;
        else
          read(is, index);

        if(index >= F.size())
        {
          std::cout << __FILE__ << ", " << __FUNCTION__ << ", l:" << __LINE__ << "  f="
                    << i << "<" << m << ", index=" << j << " nb=" << index << " #F=" << F.size()
                    << std::endl;
        }
        CGAL_triangulation_assertion(i < F.size());
        CGAL_triangulation_assertion(index < F.size());
        F[i]->set_neighbor(j, F[index]);
      }
    }
  }

  // read offsets
  int off[3] = {0, 0, 0};
  for(std::size_t j=0; j<m; ++j)
  {
    if(is_ascii(is))
    {
      is >> off[0] >> off[1] >> off[2];
    }
    else
    {
      read(is, off[0]);
      read(is, off[1]);
      read(is, off[2]);
    }
    set_offsets(F[j], off[0], off[1], off[2]);
  }

  // read potential other information
  for(std::size_t j=0; j<m; ++j)
    is >> *(F[j]);

  CGAL_triangulation_expensive_assertion(is_valid());
  return is;
}

namespace internal {

/// Internal function used by operator==.
//TODO: introduce offsets
template <class GT, class Tds1, class Tds2>
bool
test_next(const Periodic_2_triangulation_2<GT, Tds1>& t1,
          const Periodic_2_triangulation_2<GT, Tds2>& t2,
          typename Periodic_2_triangulation_2<GT, Tds1>::Face_handle f1,
          typename Periodic_2_triangulation_2<GT, Tds2>::Face_handle f2,
          std::map<typename Periodic_2_triangulation_2<GT, Tds1>::Face_handle,
          typename Periodic_2_triangulation_2<GT, Tds2>::Face_handle >& Cmap,
          std::map<typename Periodic_2_triangulation_2<GT, Tds1>::Vertex_handle,
          typename Periodic_2_triangulation_2<GT, Tds2>::Vertex_handle >& Vmap)
{
  // This function tests and registers the 4 neighbors of f1/f2,
  // and recursively calls itself over them.
  // Returns false if an inequality has been found.

  // Precondition: f1, f2 have been registered as well as their 4 vertices.
  CGAL_triangulation_precondition(t1.number_of_vertices() != 0);
  CGAL_triangulation_precondition(Cmap[f1] == f2);
  CGAL_triangulation_precondition(Vmap.find(f1->vertex(0)) != Vmap.end());
  CGAL_triangulation_precondition(Vmap.find(f1->vertex(1)) != Vmap.end());
  CGAL_triangulation_precondition(Vmap.find(f1->vertex(2)) != Vmap.end());

  typedef Periodic_2_triangulation_2<GT, Tds1>                               Tr1;
  typedef Periodic_2_triangulation_2<GT, Tds2>                               Tr2;
  typedef typename Tr1::Vertex_handle                                        Vertex_handle1;
  typedef typename Tr1::Face_handle                                          Face_handle1;
  typedef typename Tr2::Vertex_handle                                        Vertex_handle2;
  typedef typename Tr2::Face_handle                                          Face_handle2;
  typedef typename std::map<Face_handle1, Face_handle2>::const_iterator      Cit;
  typedef typename std::map<Vertex_handle1, Vertex_handle2 >::const_iterator Vit;

  for(int i=0; i<= 2; ++i)
  {
    Face_handle1 n1 = f1->neighbor(i);
    Cit cit = Cmap.find(n1);
    Vertex_handle1 v1 = f1->vertex(i);
    Vertex_handle2 v2 = Vmap[v1];
    Face_handle2 n2 = f2->neighbor(f2->index(v2));
    if(cit != Cmap.end())
    {
      // n1 was already registered.
      if(cit->second != n2)
        return false;

      continue;
    }

    // n1 has not yet been registered.
    // We check that the new vertices match geometrically.
    // And we register them.
    Vertex_handle1 vn1 = n1->vertex(n1->index(f1));
    Vertex_handle2 vn2 = n2->vertex(n2->index(f2));
    Vit vit = Vmap.find(vn1);
    if(vit != Vmap.end())
    {
      // vn1 already registered
      if(vit->second != vn2)
        return false;
    }
    else
    {
      if(t1.geom_traits().compare_xy_2_object()(vn1->point(), vn2->point()) != 0)
        return false;

      // We register vn1/vn2.
      Vmap.insert(std::make_pair(vn1, vn2));
    }

    // We register n1/n2.
    Cmap.insert(std::make_pair(n1, n2));
    // We recurse on n1/n2.
    if(!test_next(t1, t2, n1, n2, Cmap, Vmap))
      return false;
  }

  return true;
}

} // namespace internal

template<class Gt, class Tds>
std::istream&
operator>>(std::istream& is, Periodic_2_triangulation_2<Gt, Tds>& tr)
{
  return tr.load(is);
}

template<class Gt, class Tds>
std::ostream&
operator<<(std::ostream& os, Periodic_2_triangulation_2<Gt, Tds>& tr)
{
  return tr.save(os);
}

template < class GT, class Tds1, class Tds2  >
bool
operator==(const Periodic_2_triangulation_2<GT, Tds1>& t1,
           const Periodic_2_triangulation_2<GT, Tds2>& t2)
{
  typedef typename Periodic_2_triangulation_2<GT, Tds1>::Vertex_handle        Vertex_handle1;
  typedef typename Periodic_2_triangulation_2<GT, Tds1>::Face_handle          Face_handle1;
  typedef typename Periodic_2_triangulation_2<GT, Tds2>::Vertex_handle        Vertex_handle2;
  typedef typename Periodic_2_triangulation_2<GT, Tds2>::Vertex_handle        Vertex_iterator2;
  typedef typename Periodic_2_triangulation_2<GT, Tds2>::Face_handle          Face_handle2;
  typedef typename Periodic_2_triangulation_2<GT, Tds2>::Face_circulator      Face_circulator2;

  typedef typename Periodic_2_triangulation_2<GT, Tds1>::Point                Point;
  typedef typename Periodic_2_triangulation_2<GT, Tds1>::Offset               Offset;

  // Some quick checks.
  if(t1.domain() != t2.domain() || t1.number_of_sheets() != t2.number_of_sheets())
    return false;

  if(t1.number_of_vertices() != t2.number_of_vertices() || t1.number_of_faces() != t2.number_of_faces())
    return false;

  // Special case for empty triangulations
  if(t1.number_of_vertices() == 0)
    return true;

  // We will store the mapping between the 2 triangulations vertices and
  // faces in 2 maps.
  std::map<Vertex_handle1, Vertex_handle2> Vmap;
  std::map<Face_handle1, Face_handle2> Cmap;

  // find a common point
  Vertex_handle1 v1 = static_cast<Vertex_handle1>(t1.vertices_begin());
  Vertex_handle2 iv2;
  for(Vertex_iterator2 vit2 = t2.vertices_begin(); vit2 != t2.vertices_end(); ++vit2)
  {
    if(t1.compare_xy(vit2->point(), v1->point(), t2.get_offset(vit2), t1.get_offset(v1)) != EQUAL)
      continue;

    iv2 = static_cast<Vertex_handle2>(vit2);
    break;
  }

  if(iv2 == Vertex_handle2())
    return false;

  Vmap.insert(std::make_pair(v1, iv2));

  // We pick one face of t1, and try to match it against the
  // faces of t2.
  Face_handle1 c = v1->face();
  Vertex_handle1 v2 = f->vertex(t1.cw(f->index(v1)));
  Vertex_handle1 v3 = f->vertex(t1.ccw(f->index(v1)));
  Point p2 = v2->point();
  Point p3 = v3->point();
  Offset o2 = t1.get_offset(v2);
  Offset o3 = t1.get_offset(v3);

  Face_circulator2 fc = t2.incident_faces(iv2), done(fc);
  do
  {
    int inf = fc->index(iv2);

    if(t1.compare_xy(p2, fc->vertex((inf + 1) % 3)->point(),
                     o2, t2.get_offset(fc->vertex((inf + 1) % 3))) == EQUAL)
      Vmap.insert(std::make_pair(v2, fc->vertex((inf + 1) % 3)));
    else if(t1.compare_xy(p2, fc->vertex((inf + 2) % 3)->point(),
                          o2, t2.get_offset(fc->vertex((inf + 2) % 3))) == EQUAL)
      Vmap.insert(std::make_pair(v2, fc->vertex((inf + 2) % 3)));
    else
      continue; // None matched v2.

    if(t1.compare_xy(p3, fc->vertex((inf + 1) % 3)->point(),
                     o3, t2.get_offset(fc->vertex((inf + 1) % 3))) == EQUAL)
      Vmap.insert(std::make_pair(v3, fc->vertex((inf + 1) % 3)));
    else if(t1.compare_xy(p3, fc->vertex((inf + 2) % 3)->point(),
                          o3, t2.get_offset(fc->vertex((inf + 2) % 3))) == EQUAL)
      Vmap.insert(std::make_pair(v3, fc->vertex((inf + 2) % 3)));
    else
      continue; // None matched v3.

    // Found it !
    Cmap.insert(std::make_pair(c, fc));
    break;
  }
  while(++fc != done);

  if(Cmap.size() == 0)
    return false;

  // We now have one face, we need to propagate recursively.
  return internal::test_next(t1, t2, Cmap.begin()->first, Cmap.begin()->second, Cmap, Vmap);
}

template < class GT, class Tds1, class Tds2 >
inline
bool
operator!=(const Periodic_2_triangulation_2<GT, Tds1>& t1,
           const Periodic_2_triangulation_2<GT, Tds2>& t2)
{
  return ! (t1 == t2);
}
#endif // MACRO_THAT_DOESNT_EXIT_TO_MAKE_GP2T2_WORK

#define CGAL_INCLUDE_FROM_PERIODIC_2_TRIANGULATION_2_H
#include <CGAL/Periodic_2_triangulation_dummy_12.h>
#undef CGAL_INCLUDE_FROM_PERIODIC_2_TRIANGULATION_2_H

} //namespace CGAL

#endif //CGAL_PERIODIC_2_TRIANGULATION_2_H
