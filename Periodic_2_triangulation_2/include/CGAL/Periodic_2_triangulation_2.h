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
//
// Author(s)     : Nico Kruithof <Nico@nghk.nl>

#ifndef CGAL_PERIODIC_2_TRIANGULATION_2_H
#define CGAL_PERIODIC_2_TRIANGULATION_2_H

#include <CGAL/basic.h>

#include <list>
#include <vector>
#include <map>
#include <algorithm>
#include <utility>
#include <iostream>

#include <CGAL/iterator.h>
#include <CGAL/Iterator_project.h>
#include <CGAL/function_objects.h>

#include <CGAL/triangulation_assertions.h>
#include <CGAL/Triangulation_utils_2.h>

#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/Triangulation_vertex_base_2.h>
#include <CGAL/Triangulation_face_base_2.h>
#include <CGAL/Periodic_2_triangulation_face_base_2.h>
#include <CGAL/Periodic_2_triangulation_vertex_base_2.h>
#include <CGAL/Periodic_2_triangulation_iterators_2.h>
#include <CGAL/spatial_sort.h>

#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_smallint.hpp>
#include <boost/random/variate_generator.hpp>

#include <CGAL/utility.h>

namespace CGAL
{

/// Periodic triangulation class.
/// Its main functionality is:
/// - Insertion of points
/// - Deletion of points
/// - Point location
template < class Gt,
         class Tds = Triangulation_data_structure_2 <
         Periodic_2_triangulation_vertex_base_2<Gt>,
         Periodic_2_triangulation_face_base_2<Gt> > >
class Periodic_2_triangulation_2: public Triangulation_cw_ccw_2
{
  typedef Periodic_2_triangulation_2<Gt, Tds> Self;

public:
  // Public types of Periodic_2_triangulation_2
  /// The triangulation data structure type
  typedef Tds Triangulation_data_structure;
  /// The traits class
  typedef Gt Geom_traits;

  /// The periodic offset type
  typedef typename Gt::Periodic_2_offset_2 Offset;
  /// The iso rectangle type
  typedef typename Gt::Iso_rectangle_2 Iso_rectangle;
  /// Integer tuple to store the number of sheets in each direction of space.
  typedef array<int, 2> Covering_sheets;

  /// The point type
  typedef typename Gt::Point_2 Point;
  /// The segment type
  typedef typename Gt::Segment_2 Segment;
  /// The triangle type
  typedef typename Gt::Triangle_2 Triangle;

  /// Represents a point-offset pair. The point in the pair lies in the original domain.
  typedef std::pair<Point, Offset> Periodic_point;
  /// A pair of periodic points representing a segment in the periodic domain.
  typedef array<std::pair<Point, Offset>, 2> Periodic_segment;
  /// A triple of periodic points representing a triangle in the periodic domain.
  typedef array<std::pair<Point, Offset>, 3> Periodic_triangle;

  /// The vertex type
  typedef typename Tds::Vertex Vertex;
  /// The face type
  typedef typename Tds::Face Face;
  /// The edge type
  typedef typename Tds::Edge Edge;

  /// Size type (an unsigned integral type)
  typedef typename Tds::size_type size_type;
  /// Difference type (a signed integral type)
  typedef typename Tds::difference_type difference_type;

  /// Handle to a vertex
  typedef typename Tds::Vertex_handle Vertex_handle;
  /// Handle to a face
  typedef typename Tds::Face_handle Face_handle;

  /// Iterator over the faces
  typedef typename Tds::Face_iterator Face_iterator;
  /// Iterator over the edges
  typedef typename Tds::Edge_iterator Edge_iterator;
  /// Iterator over the vertices
  typedef typename Tds::Vertex_iterator Vertex_iterator;
  /// Iterator over the vertices whose corresponding points lie in the
  /// original domain, i.e. for each set of periodic copies the
  /// Unique_vertex_iterator iterates over exactly one representative.
  typedef Periodic_2_triangulation_unique_vertex_iterator_2<Self>
  Unique_vertex_iterator;

  /// \name For compatibility with the Triangulation_2 class
  // \{
  typedef Face_iterator Finite_faces_iterator;
  typedef Edge_iterator Finite_edges_iterator;
  typedef Vertex_iterator Finite_vertices_iterator;
  typedef Face_iterator All_faces_iterator;
  // \}

  /// Circulator over all faces incident to a vertex
  typedef typename Tds::Face_circulator Face_circulator;
  /// Circulator over all edges incident to a vertex
  typedef typename Tds::Edge_circulator Edge_circulator;
  /// Circulator over all vertices incident to a vertex
  typedef typename Tds::Vertex_circulator Vertex_circulator;

  /// \name Periodic iterator types
  //\{
  /// Iterator over all periodic triangles
  typedef Periodic_2_triangulation_triangle_iterator_2<Self>
  Periodic_triangle_iterator;
  /// Iterator over all periodic segments
  typedef Periodic_2_triangulation_segment_iterator_2<Self>
  Periodic_segment_iterator;
  /// Iterator over all periodic points
  typedef Periodic_2_triangulation_point_iterator_2<Self>
  Periodic_point_iterator;
  //\}

  /// \name Enumeration types
  //\{

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
    EMPTY
  };
  //\}

  // Auxiliary iterators for convenience
  // do not use default template argument to please VC++
  /// Functor that returns the point given a vertex
  typedef Project_point<Vertex> Proj_point;

  /// \name STL types
  // \{
  /// value_type similar to stl containers
  typedef Point value_type; // to have a back_inserter
  /// const_reference similar to stl containers
  typedef const value_type& const_reference;
  /// reference similar to stl containers
  typedef value_type& reference;
  // \}


  /// Tag to distinguish Regular triangulations from others;
  typedef Tag_false Weighted_tag;

protected:
  // Protected types of Periodic_2_triangulation_2
  typedef typename Gt::Orientation_2 Orientation_2;
  typedef typename Gt::Compare_x_2 Compare_x;
  typedef typename Gt::Compare_y_2 Compare_y;


  typedef typename Gt::FT FT;
  typedef std::pair<Vertex_handle, Offset> Virtual_vertex;
  typedef std::map<Vertex_handle, Virtual_vertex> Virtual_vertex_map;
  typedef typename Virtual_vertex_map::const_iterator Virtual_vertex_map_it;

  /// Vector is contains virtual copies with offset off:
  /// virtual copy with offset off is stored at position: i=3*off[0]+off[1]-1
  typedef std::map<Vertex_handle, std::vector<Vertex_handle> >
  Virtual_vertex_reverse_map;
  typedef typename Virtual_vertex_reverse_map::const_iterator
  Virtual_vertex_reverse_map_it;
  typedef std::map<Vertex_handle, std::list<Vertex_handle> > Too_long_edges_map;
  typedef typename Too_long_edges_map::const_iterator Too_long_edges_map_it;

  /// \name Functors
  // \{
  /// Functor for symbolically perturbing points
  class Perturbation_order
  {
    const Self *t;

  public:
    // Perturbation_order, public interface
    Perturbation_order(const Self *tr) :
      t(tr)
    {
    }

    bool operator()(const Point *p, const Point *q) const
    {
      return t->compare_xy(*p, *q) == SMALLER;
    }
    bool operator()(const Periodic_point *p, const Periodic_point *q) const
    {
      return t->compare_xy(p->first, q->first, p->second, q->second) == SMALLER;
    }
  };
  // \}
  friend class Perturbation_order;

private:
  // Copy constructor helpers
  class Finder;
  void copy_multiple_covering(const Periodic_2_triangulation_2 & tr);

public:
  // Public functions of Periodic_2_triangulation_2
  /// \name Constructors
  //\{
  /// Constructor
  Periodic_2_triangulation_2(
    const Iso_rectangle &domain = Iso_rectangle(0, 0, 1, 1),
    const Geom_traits &geom_traits = Geom_traits());

  /// Copy constructor
  Periodic_2_triangulation_2(const Periodic_2_triangulation_2<Gt, Tds> &tr);
  /// Assignment
  Periodic_2_triangulation_2 &operator=(const Periodic_2_triangulation_2 &tr);

  /// Copy the triangulation
  void copy_triangulation(const Periodic_2_triangulation_2 &tr);
  /// Swap two triangulations
  void swap(Periodic_2_triangulation_2 &tr);
  /// Clear the triangulation
  void clear();

  /// Serialize the triangulation to an output stream
  std::ostream& save(std::ostream& os) const;

  /// Deserialize the triangulation from an input stream
  std::istream& load(std::istream& is);

  //\}

  /// \name Access functions
  //\{

  /// Returns the geometric traits used for the predicates and constructions.
  const Geom_traits& geom_traits() const
  {
    return _gt;
  }
  /// Returns the datastructure storing the triangulation.
  const Triangulation_data_structure & tds() const
  {
    return _tds;
  }
  /// Returns the datastructure storing the triangulation.
  Triangulation_data_structure & tds()
  {
    return _tds;
  }
  /// Returns the domain of the 1-sheeted cover.
  const Iso_rectangle & domain() const
  {
    return _domain;
  }
  /// Returns the number of copies of the 1-sheeted cover stored in each of
  /// the principal directions.
  Covering_sheets number_of_sheets() const
  {
    return _cover;
  }
  /// Returns the dimension of the triangulation.
  int dimension() const
  {
    return _tds.dimension() == 2 ? 2 : 0;
  }

  //\}

  /// \name Number of simplices
  //\{
  /// Returns whether the triangulation is empty.
  bool empty() const
  {
    return _tds.dimension() < 2;
  }

  /// Returns the number of vertices. Counts all vertices that are
  /// representatives of the same point in the 1-cover as one vertex.
  size_type number_of_vertices() const
  {
    if (is_1_cover())
      return _tds.number_of_vertices();
    else
      return _tds.number_of_vertices() / 9;
  }

  /// Returns the number of edges. Counts all edges that are
  /// representatives of the same segment in the 1-cover as one edge.
  size_type number_of_edges() const
  {
    if (is_1_cover())
      return _tds.number_of_edges();
    else
      return _tds.number_of_edges() / 9;
  }
  /// Returns the number of faces. Counts all faces that are
  /// representatives of the same triangle in the 1-cover as one face.
  size_type number_of_faces() const
  {
    if (is_1_cover())
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
  //\}


  /// \name Methods regarding the covering
  /// \{

  /// Checks whether the triangulation is a valid simplicial complex in the one cover.
  /// Uses an edge-length-criterion.
  bool is_extensible_triangulation_in_1_sheet_h1() const;

  /// Checks whether the triangulation is a valid simplicial complex in the one cover.
  /// Uses a criterion based on the maximal radius of the circumscribing circle.
  bool is_extensible_triangulation_in_1_sheet_h2() const;

  /// Checks whether the triangulation is a valid simplicial complex in the one cover.
  bool is_triangulation_in_1_sheet() const;

  /// Convert a 9 sheeted cover (used for sparse triangulations) to a single sheeted cover.
  /// \pre !is_1_cover();
  void convert_to_1_sheeted_covering();
  /// Convert a single sheeted cover (used for dense triangulations) to a 9 sheeted cover.
  /// \pre is_1_cover();
  void convert_to_9_sheeted_covering();
  // \}

  /// \name Geometric access functions
  // \{

  /// Returns the periodic point given by vertex v. If t is
  /// represented in the 1-sheeted covering space, the offset is
  /// always zero. Otherwise v can correspond to a periodic copy
  /// outside domain of an input point.
  Periodic_point periodic_point(const Vertex_handle &v) const
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
  Periodic_point periodic_point(const Face_handle &f, int i) const
  {
    return Periodic_point(f->vertex(i)->point(), get_offset(f, i));
  }

  /// Returns the periodic segment formed by the two point-offset
  /// pairs corresponding to the two vertices of edge (f,i).
  /// \pre i == {0,1,2}
  Periodic_segment periodic_segment(const Face_handle &f, int i) const
  {
    CGAL_triangulation_precondition( number_of_vertices() != 0 );
    CGAL_triangulation_precondition( i >= 0 && i <= 2);
    return make_array(periodic_point(f, ccw(i)),
                      periodic_point(f, cw(i)));
  }

  /// Same as the previous method for edge e.
  Periodic_segment periodic_segment(const Edge &e) const
  {
    return periodic_segment(e.first, e.second);
  }

  /// Returns the periodic triangle formed by the three point-offset
  /// pairs corresponding to the three vertices of facet f.
  Periodic_triangle periodic_triangle(const Face_handle &f) const
  {
    return make_array(periodic_point(f, 0),
                      periodic_point(f, 1),
                      periodic_point(f, 2));
  }

  /// Converts the Periodic_point pp (point-offset pair) to the corresponding
  /// Point in \f$R^2\f$.
  Point point(const Periodic_point & pp) const
  {
    return construct_point(pp.first, pp.second);
  }
  Point point(const Vertex_handle &v) const
  {
    return point(periodic_point(v));
  }
  Point point(const Face_handle &fh, int i) const
  {
    return point(periodic_point(fh, i));
  }
  /// Converts the Periodic_segment ps to a Segment in \f$R^2\f$.
  Segment segment(const Periodic_segment &ps) const
  {
    return construct_segment(ps[0].first, ps[1].first, ps[0].second, ps[1].second);
  }
  /// Converts the Periodic_triangle pt to a Triagle in \f$R^2\f$.
  Triangle triangle(const Periodic_triangle &pt) const
  {
    Triangle triang = construct_triangle(pt[0].first, pt[1].first, pt[2].first,
                                         pt[0].second, pt[1].second, pt[2].second);
    return triang;
  }

  /// Constructs the segment associated with the edge (f,i), respects the offset
  Segment segment(Face_handle f, int i) const
  {
    return segment(periodic_segment(f, i));
  }
  /// Constructs the segment associated with the edge e, respects the offset
  Segment segment(const Edge& e) const
  {
    return segment(periodic_segment(e));
  }
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
  Triangle triangle(Face_handle f) const
  {
    return triangle(periodic_triangle(f));
  }
  //\}

  Point move_in_domain(const Point &p)
  {
    typename Gt::FT x = p.x();
    typename Gt::FT y = p.y();

    while (x < _domain.xmin()) x += _domain.xmax() - _domain.xmin();
    while (x >= _domain.xmax()) x -= _domain.xmax() - _domain.xmin();

    while (y < _domain.ymin()) y += _domain.ymax() - _domain.ymin();
    while (y >= _domain.ymax()) y -= _domain.ymax() - _domain.ymin();

    return Point(x, y);
  }

  /// \name Queries on simplices
  // \{
  bool is_edge(Vertex_handle va, Vertex_handle vb) const
  {
    return _tds.is_edge(va, vb);
  }
  bool is_edge(Vertex_handle va, Vertex_handle vb, Face_handle& fr, int & i) const
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
  // \}

  /// \name Queries
  // \{

  /// Wrapper function for locate if only the requested point is given.
  Face_handle locate(const Point &p, Face_handle start = Face_handle()) const
  {
    Locate_type lt;
    int li;
    return locate(p, Offset(), lt, li, start);
  }

  /// Wrapper function for locate if the offset is omitted.
  Face_handle locate(const Point& p, Locate_type& lt, int& li,
                     Face_handle start = Face_handle()) const
  {
    return locate(p, Offset(), lt, li, start);
  }

  /// Returns the oriented side of the point p with respect to the
  /// triangle defined by the face f
  Oriented_side oriented_side(Face_handle f, const Point& p) const
  {
    return oriented_side(f, p, Offset());
  }

  // \}


  /// \name Traversal of the Triangulation
  // \{
  /// Iterator over all stored vertices. Starts at an arbitrary vertex.
  /// Returns vertices_end() if t.number_of_vertices()=0.
  Vertex_iterator vertices_begin() const
  {
    return _tds.vertices_begin();
  }
  /// Past the end Vertex_iterator.
  Vertex_iterator vertices_end() const
  {
    return _tds.vertices_end();
  }
  /// Iterator over all stored edges. Starts at an arbitrary edge.
  /// Returns edges_end() if t.number_of_vertices()=0.
  Edge_iterator edges_begin() const
  {
    return _tds.edges_begin();
  }
  /// Past the end Edge_iterator.
  Edge_iterator edges_end() const
  {
    return _tds.edges_end();
  }
  /// Iterator over all stored faces. Starts at an arbitrary face.
  /// Returns faces_end() if t.number_of_vertices()=0.
  Face_iterator faces_begin() const
  {
    return _tds.faces_begin();
  }
  /// Past the end Face_iterator.
  Face_iterator faces_end() const
  {
    return _tds.faces_end();
  }

  /// Iterator over all stored vertices. Starts at an arbitrary vertex.
  /// Returns vertices_end() if t.number_of_vertices()=0.
  Vertex_iterator finite_vertices_begin() const
  {
    return _tds.vertices_begin();
  }
  /// Past the end Vertex_iterator.
  Vertex_iterator finite_vertices_end() const
  {
    return _tds.vertices_end();
  }
  /// Iterator over all stored edges. Starts at an arbitrary edge.
  /// Returns edges_end() if t.number_of_vertices()=0.
  Edge_iterator finite_edges_begin() const
  {
    return _tds.edges_begin();
  }
  /// Past the end Edge_iterator.
  Edge_iterator finite_edges_end() const
  {
    return _tds.edges_end();
  }
  /// Iterator over all stored faces. Starts at an arbitrary face.
  /// Returns faces_end() if t.number_of_vertices()=0.
  Face_iterator finite_faces_begin() const
  {
    return _tds.faces_begin();
  }
  /// Past the end Face_iterator.
  Face_iterator finite_faces_end() const
  {
    return _tds.faces_end();
  }

  /// Iterator over all stored vertices. Starts at an arbitrary vertex.
  /// Returns vertices_end() if t.number_of_vertices()=0.
  Vertex_iterator all_vertices_begin() const
  {
    return _tds.vertices_begin();
  }
  /// Past the end Vertex_iterator.
  Vertex_iterator all_vertices_end() const
  {
    return _tds.vertices_end();
  }
  /// Iterator over all stored edges. Starts at an arbitrary edge.
  /// Returns edges_end() if t.number_of_vertices()=0.
  Edge_iterator all_edges_begin() const
  {
    return _tds.edges_begin();
  }
  /// Past the end Edge_iterator.
  Edge_iterator all_edges_end() const
  {
    return _tds.edges_end();
  }
  /// Iterator over all stored faces. Starts at an arbitrary face.
  /// Returns faces_end() if t.number_of_vertices()=0.
  Face_iterator all_faces_begin() const
  {
    return _tds.faces_begin();
  }
  /// Past the end Face_iterator.
  Face_iterator all_faces_end() const
  {
    return _tds.faces_end();
  }

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

  // \}
  /// \name Geometric iterators
  //\{

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
  //\}

  /// \name Incident simplices
  // \{

  Face_circulator incident_faces(Vertex_handle v, Face_handle f = Face_handle()) const
  {
    return _tds.incident_faces(v, f);
  }
  Edge_circulator incident_edges(Vertex_handle v, Face_handle f = Face_handle()) const
  {
    return _tds.incident_edges(v, f);
  }
/*   Vertex_circulator incident_vertices(Vertex_handle v, Face_handle f = */
/*                                         Face_handle()) const */
/*   { */
/*     bool DEPRECATED_USE_ADJACENT_VERTICES; */
/*     return adjacent_vertices(v, f); */
/*   } */
  Vertex_circulator adjacent_vertices(Vertex_handle v, Face_handle f =
                                        Face_handle()) const
  {
    return _tds.incident_vertices(v, f);
  }
  // \}

  /// \name Traversal between adjacent faces
  // \{

  Vertex_handle mirror_vertex(Face_handle f, int i) const
  {
    return _tds.mirror_vertex(f, i);
  }
  int mirror_index(Face_handle f, int i) const
  {
    return _tds.mirror_index(f, i);
  }
  //\}


  /// \name Modifiers
  // \{

  /// Flips the edge and all periodic copies
  void flip(Face_handle f, int i);

  /// Inserts a point in the triangulation
  /// \param p the point to be inserted
  /// \param start the start face for point location
  /// \return The new vertex handle or an existing Vertex_handle if p was inserted before
  Vertex_handle insert(const Point &p, Face_handle start = Face_handle());

  /// Inserts a point in the triangulation
  /// \pre The point has been located in the triangulation
  Vertex_handle insert(const Point& p, Locate_type lt, Face_handle loc, int li);

  /// Insert a point in the triangulation
  Vertex_handle push_back(const Point& p)
  {
    return insert(p);
  }

  // \}

  /// \name Advanced modifiers
  //\{

  /// Insert the first vertex in the triangulation and creates the 9-cover.
  Vertex_handle insert_first(const Point& p);
  /// Inserts p in the face f and sets the offsets of the newly created faces
  /// Insert periodic copies in all periodic copies of the domain
  Vertex_handle insert_in_face(const Point& p, Face_handle f);
  /// Inserts (p,o) in the edge (f,i) and sets the offsets of the newly created faces
  /// Insert periodic copies in all periodic copies of the domain
  Vertex_handle insert_in_edge(const Point& p, Face_handle f, int i);

  /// Remove a degree 3 vertex from a 2D triangulation
  void remove_degree_3(Vertex_handle v);

  bool remove_degree_init(Vertex_handle v, const Offset &v_o,
                          std::vector<Face_handle> &f,
                          std::vector<Vertex_handle> &w, std::vector<Offset> &offset_w,
                          std::vector<int> &i, int &d, int &maxd, bool &simplicity_criterion);

  /// Remove a vertex from a 2D triangulation with number_of_vertices() == 1
  void remove_first(Vertex_handle v);

  /// Changes the domain. Note that this function calls clear(), i.e.,
  /// it erases the existing triangulation.
  void set_domain(const Iso_rectangle &domain)
  {
    clear();
    _domain = domain;
    _gt.set_domain(_domain);
    _edge_length_threshold =
      FT(0.166) * (_domain.xmax() - _domain.xmin()) * (_domain.xmax() - _domain.xmin());
  }
  //\}


  /// \name Point location

  /// Do a remembering heuristic walk to locate point (p,o)
  Face_handle
  march_locate_2D(Face_handle f, const Point& p, const Offset& o,
                  Locate_type& lt, int& li) const;

  /// Checks whether the result of two point location queries are equivalent.
  bool
  compare_walks(const Point& p, Face_handle c1, Face_handle c2,
                Locate_type& lt1, Locate_type& lt2, int li1, int li2) const;

  /// Testing where the point (p,off) lies w.r.t. the face f
  Bounded_side side_of_face(const Point &p, const Offset &off, Face_handle f,
                            Locate_type &lt, int &li) const;
  /// Testing where the point (p,off) lies w.r.t. the face f
  Bounded_side side_of_face(const Point &p, Face_handle f, Locate_type &lt,
                            int &li) const
  {
    return side_of_face(p, Offset(), f, lt, li);
  }
  //\}

  /// \name Predicates and Constructions
  //\{
  /// Determines whether the point p lies on the (un-)bounded side of
  /// the triangle (p0,p1,p2)
  Bounded_side
  bounded_side(const Point &p0, const Point &p1, const Point &p2, const Point &p) const;

  /// Determines whether the point (p,o) lies on the (un-)bounded side of
  /// the triangle ((p0,o0),(p1,o1),(p2,o2))
  Bounded_side
  bounded_side(const Point &p0, const Point &p1, const Point &p2, const Point &p,
               const Offset &o0, const Offset &o1, const Offset &o2, const Offset &o) const;

  /// Determines whether the point q lies strictly between the points p and r
  /// p,q and r are supposed to be collinear points
  bool
  collinear_between(const Point& p, const Point& q, const Point& r) const;

  /// Determines whether the point (q,o_q) lies strictly between the points (p,o_p) and (r,o_r)
  /// (q,o_q), (p,o_p) and (r,o_r) are supposed to be collinear points
  bool
  collinear_between(const Point& p, const Point& q, const Point& r,
                    const Offset& o_p, const Offset& o_q, const Offset& o_r) const;

  /// Compares the x-coordinates of p and q
  Comparison_result compare_x(const Point& p, const Point& q) const;
  /// Compares the x-coordinates of (p,o_p) and (q,o_q)
  Comparison_result compare_x(const Point& p, const Point& q,
                              const Offset &o_p, const Offset &o_q) const;

  /// Compares p and q lexicographically
  Comparison_result compare_xy(const Point& p, const Point& q,
                               const Offset &o_p, const Offset &o_q) const;
  /// Compares (p,o_p) and (q,o_q) lexicographically
  Comparison_result compare_xy(const Point& p, const Point& q) const;
  /// Compares the x-coordinates of p and q
  Comparison_result compare_y(const Point& p, const Point& q) const;
  /// Compares the x-coordinates of (p,o_p) and (q,o_q)
  Comparison_result compare_y(const Point& p, const Point& q,
                              const Offset &o_p, const Offset &o_q) const;
  /// Checks for equality of p and q
  bool xy_equal(const Point& p, const Point& q) const;
  /// Returns the orientation of p1,p2,p3
  Orientation orientation(const Point& p1, const Point& p2, const Point& p3) const;
  /// Returns the orientation of (p1,o1), (p2,o2), (p3,o3)
  Orientation orientation(const Point& p1, const Point& p2, const Point& p3,
                          const Offset& o1, const Offset& o2, const Offset& o3) const;

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
  //\}


  /// \name Miscellaneous
  //\{

  /// Returns whether the union of the faces f and f->neighbor(i) form
  /// a convex quadrilateral.
  bool flippable(Face_handle f, int i);

  size_type degree(Vertex_handle v) const
  {
    return _tds.degree(v);
  }

  /// Checks if the triangulation is valid.
  bool is_valid(bool verbose = false, int level = 0) const;
  /// Checks if the face is valid.
  bool is_valid(Face_handle fh, bool verbose = false, int level = 0) const;

  //\}


  /// \name Undocumented functions, needed by the geometric iterators
  // \{
  /// [Undoc] Returns whether the stored triangulation covers a 1-cover.
  bool is_1_cover() const
  {
    return (_cover[0] == 1) && (_cover[1] == 1);
  }

  /// [Undoc] Combines two offsets, where the first offset is defined by the
  /// virtual vertex and the second by the face.
  Offset combine_offsets(const Offset& o_c, const Offset& o_t) const
  {
    Offset o_ct(_cover[0] * o_t.x(), _cover[1] * o_t.y());
    return o_c + o_ct;
  }
  /// [Undoc] Returns the offset of nb==ch->neighbor(i) with respect to ch.
  /// Get the offset between the origins of the internal offset coordinate
  /// systems of two neighboring faces with respect from ch to nb.
  ///
  /// - Find two corresponding vertices from each face
  /// - Return the difference of their offsets.
  ///
  Offset get_neighbor_offset(Face_handle fh, int i) const
  {
    Face_handle nb = fh->neighbor(i);
    return get_neighbor_offset(fh, i, nb, nb->index(fh));
  }
  /// [Undoc] Returns the offset of nb==ch->neighbor(i) with respect to ch.
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
  /// [Undoc] returns the combined offset of the vertex
  /// (if we are not on the 1-cover) and the offset defined by the face.
  Offset get_offset(Face_handle f, int i) const
  {
    if (is_1_cover())
      return int_to_off(f->offset(i));

    Virtual_vertex_map_it it = _virtual_vertices.find(f->vertex(i));
    if (it != _virtual_vertices.end())
      return combine_offsets(it->second.second, int_to_off(f->offset(i)));
    else
      return combine_offsets(Offset(), int_to_off(f->offset(i)));
  }
  /// [Undoc] Returns the offset of the vertex if we are not on the 1-cover.
  Offset get_offset(Vertex_handle vh) const
  {
    if (is_1_cover())
      return Offset();
    Virtual_vertex_map_it it = _virtual_vertices.find(vh);
    if (it != _virtual_vertices.end())
      return it->second.second;
    else
      return Offset();
  }

  /// Converts an offset to a bit pattern where bit1==offx and bit0==offy.
  int off_to_int(const Offset & off) const
  {
    CGAL_triangulation_assertion( off.x() == 0 || off.x() == 1 );
    CGAL_triangulation_assertion( off.y() == 0 || off.y() == 1 );
    int i = ((off.x() & 1) << 1) + (off.y() & 1);
    return i;
  }
  /// Creates an offset from a bit pattern.
  Offset int_to_off(int i) const
  {
    return Offset((i >> 1) & 1, i & 1);
  }

  // \}
  // Protected functions of Periodic_2_triangulation_2
  /// Const accessor to the virtual vertices reverse map,
  /// used to optimize point location for periodic copies.
  const Virtual_vertex_reverse_map &virtual_vertices_reverse() const
  {
    return _virtual_vertices_reverse;
  }

  /// [Undoc] Returns the non-virtual copy of the vertex.
  Vertex_handle get_original_vertex(Vertex_handle vh) const
  {
    if (is_1_cover())
      return vh;
    Virtual_vertex_map_it it = _virtual_vertices.find(vh);
    if (it != _virtual_vertices.end())
      return it->second.first;
    else
      return vh;
  }


  /// Tests whether a vertex is a periodic copy of a vertex in the 3-cover.
  bool is_virtual(Vertex_handle v)
  {
    if (is_1_cover())
      return false;
    return (_virtual_vertices.find(v) != _virtual_vertices.end());
  }

  const std::vector<Vertex_handle>& periodic_copies(const Vertex_handle v) const
  {
    CGAL_triangulation_precondition(number_of_sheets() != make_array(1, 1) );
    CGAL_triangulation_precondition(_virtual_vertices.find(v) == _virtual_vertices.end());
    CGAL_triangulation_assertion(_virtual_vertices_reverse.find(v) != _virtual_vertices_reverse.end());
    return _virtual_vertices_reverse.find(v)->second;
  }

  template<class Stream>
  Stream& draw_triangulation(Stream& os) const
  {
    Edge_iterator it = edges_begin();
    for (; it != edges_end(); ++it)
      {
        os << segment(it);
      }
    return os;
  }
protected:
  std::vector<Vertex_handle> insert_dummy_points();


  inline void try_to_convert_to_one_cover()
  {
    // Fall back to 1-cover if the criterion that the longest edge is shorter
    // than sqrt(0.166) is fulfilled.
    if (_too_long_edge_counter == 0 && !is_1_cover())
      {
        CGAL_triangulation_expensive_assertion( is_valid() );
        this->convert_to_1_sheeted_covering();
        CGAL_triangulation_expensive_assertion( is_valid() );
      }
  }

protected:
  // Protected functions of Periodic_2_triangulation_2

  /// Inserts a point with an offset in the triangulation
  /// \pre The point has been located in the triangulation
  Vertex_handle insert(const Point& p, const Offset& o, Locate_type lt,
                       Face_handle loc, int li, Vertex_handle vh);

  /// \name Helper functions for queries
  // \{

  /// Locates the simplex containing the point represented by p and o.
  ///
  /// The type of the simplex is stored in lt.
  /// The simplex containing the point is returned using lt and li.
  /// The Face_handle start is the start point of the heuristic walk.
  Face_handle
  locate(const Point& p, const Offset &o, Locate_type& lt, int& li,
         Face_handle start = Face_handle()) const;
  /// Returns the oriented side of the point (p,o) with respect to the
  /// triangle defined by the face f
  Oriented_side
  oriented_side(Face_handle f, const Point& p, const Offset &o) const;
  // \}

  /// \name Insertion helpers
  //\{
  /// Insert too long edges in the star of v
  void insert_too_long_edges_in_star(Vertex_handle v);

  /// Insert too long edge
  void insert_too_long_edge(Face_handle f, int i);

  /// Remove too long edges in the star of v
  void remove_too_long_edges_in_star(Vertex_handle v);

  /// Removes an edge if it is too long
  void remove_too_long_edge(Face_handle f, int i);

  /// Check whether an edge is too long
  bool edge_is_too_long(const Point &p1, const Point &p2) const;

  /// Flips the edge, no periodic copies are flipped
  void flip_single_edge(Face_handle f, int i);

  /// Remove a vertex from the virtual copies maps
  /// Used when a Delaunay vertex is removed
  void remove_from_virtual_copies(Vertex_handle v);
  //\}

  /// \name Wrapping the traits
  //\{
  Point construct_point(const Point& p, const Offset &o) const
  {
    return geom_traits().construct_point_2_object()(p, o);
  }
  Point construct_point(const Periodic_point& pp) const
  {
    return construct_point(pp.first, pp.second);
  }
  Triangle construct_triangle(const Point &p1, const Point &p2,
                              const Point &p3, const Offset &o1, const Offset &o2, const Offset &o3) const
  {
    return geom_traits().construct_triangle_2_object()(p1, p2, p3, o1, o2, o3);
  }
  Triangle construct_triangle(const Periodic_triangle& tri) const
  {
    return construct_triangle(tri[0].first, tri[1].first, tri[2].first,
                              tri[0].second, tri[1].second, tri[2].second);
  }
  Segment construct_segment(const Point &p1, const Point &p2, const Offset &o1,
                            const Offset &o2) const
  {
    return geom_traits().construct_segment_2_object()(p1, p2, o1, o2);
  }
  Segment construct_segment(const Periodic_segment& seg) const
  {
    return construct_segment(seg[0].first, seg[1].first, seg[0].second,
                             seg[1].second);
  }
  //\}

  /// Test whether removing vertex v decreases the dimension of the triangulation.
  bool test_dim_down(Vertex_handle /*v*/) const
  {
    //test the dimensionality of the resulting triangulation
    //upon removing of vertex v
    return number_of_vertices() == 1;
  }
  void make_hole(Vertex_handle v, std::list<Edge> & hole);


  Face_handle create_face(Face_handle f1, int i1, Face_handle f2, int i2, Face_handle f3, int i3);
  Face_handle create_face(Face_handle f1, int i1, Face_handle f2, int i2);
  Face_handle create_face(Face_handle f, int i, Vertex_handle v);
  Face_handle create_face(Vertex_handle v1, Vertex_handle v2, Vertex_handle v3);
  Face_handle create_face(Vertex_handle v1, Vertex_handle v2, Vertex_handle v3, Face_handle f1, Face_handle f2, Face_handle f3);
  Face_handle create_face();
  Face_handle create_face(Face_handle); //calls copy constructor of Face
  void delete_face(Face_handle f);
  void delete_vertex(Vertex_handle v);

  // template members

  bool well_oriented(Vertex_handle v) const
  {
    Face_circulator fc = incident_faces(v), done(fc);
    do
      {
        Orientation o;

        Vertex_handle v0 = fc->vertex(0);
        Vertex_handle v1 = fc->vertex(1);
        Vertex_handle v2 = fc->vertex(2);
        if (fc->has_zero_offsets())
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
        if (o != COUNTERCLOCKWISE) return false;
      }
    while (++fc != done);
    return true;
  }

  template<class EdgeIt>
  Vertex_handle star_hole(const Point& p, EdgeIt edge_begin, EdgeIt edge_end)
  {
    std::list<Face_handle> empty_list;
    return star_hole(p, edge_begin, edge_end, empty_list.begin(),
                     empty_list.end());
  }

  template<class EdgeIt, class FaceIt>
  Vertex_handle star_hole(const Point& p, EdgeIt edge_begin, EdgeIt edge_end,
                          FaceIt face_begin, FaceIt face_end)
  {
    CGAL_assertion(is_1_cover());

    Vertex_handle v = _tds.star_hole(edge_begin, edge_end, face_begin, face_end);
    v->set_point(p);
    return v;
  }

  /// Periodic functions
  //\{

  /// These functions give the pair (vertex, offset) that corresponds
  /// to the i-th vertex of face f. The vertex returned is not a virtual copy.
  void get_vertex(Face_handle f, int i, Vertex_handle &vh, Offset &off) const;
  /// These functions give the pair (vertex, offset) that corresponds
  /// to the i-th vertex of vertex vh. The vertex returned is not a virtual copy.
  void get_vertex(Vertex_handle vh_i, Vertex_handle &vh, Offset &off) const;
  /// Returns the face containing the three vertices defined by vh[0], vh1[1] and vh[2].
  inline Face_handle get_face(const Vertex_handle* vh) const;
  /// Constructs a list of too long edges in the triangulation.
  int find_too_long_edges(
    std::map<Vertex_handle, std::list<Vertex_handle> >& edges) const;

  /// Returns the offset such that (p, o) lies on the bounded side of the face f.
  Offset get_location_offset(Face_handle f, const Point &p, const Offset &o) const
  {
    CGAL_triangulation_precondition( number_of_vertices() != 0 );

    if (is_1_cover() && f->has_zero_offsets())
      {
        // default case:
        return Offset();
      }
    else
      {
        int cumm_off = f->offset(0) | f->offset(1) | f->offset(2);
        // Special case for the periodic space.
        // Fetch vertices and respective offsets of c from _virtual_vertices
        const Point *pts[3];
        Offset off[3];
        for (int i = 0; i < 3; i++)
          {
            pts[i] = &(f->vertex(i)->point());
            off[i] = get_offset(f, i);
          }

        // Main idea seems to just test all possibilities.
        for (int i = 0; i < 4; i++)
          {
            if (((cumm_off | (~i)) & 3) == 3)
              {
                if (bounded_side(*pts[0], *pts[1], *pts[2], p, off[0], off[1],
                                 off[2], combine_offsets(o, int_to_off(i))) != ON_UNBOUNDED_SIDE)
                  {
                    return int_to_off(i);
                  }
              }
          }
      }
    CGAL_assertion(false);
    return Offset();
  }

  /// Assigns the offsets to the vertices of the face f, and makes the offset minimal in each direction.
  void set_offsets(Face_handle f, int o0, int o1, int o2)
  {
    int off0[2] = { (o0 >> 1) & 1, (o0 & 1) };
    int off1[2] = { (o1 >> 1) & 1, (o1 & 1) };
    int off2[2] = { (o2 >> 1) & 1, (o2 & 1) };
    // Make sure that there is at least one zero offset in every direction
    for (int i = 0; i < 2; i++)
      {
        int min_off = (std::min)((std::min)(off0[i], off1[i]), off2[i]);
        if (min_off != 0)
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

  /// Assigns the offsets to the vertices of the face f, and makes the offset minimal in each direction.
  template<class Offset>
  void set_offsets(Face_handle f, const Offset &o0, const Offset &o1, const Offset &o2)
  {
    int off0[2] = { o0.x(), o0.y() };
    int off1[2] = { o1.x(), o1.y() };
    int off2[2] = { o2.x(), o2.y() };
    for (int i = 0; i < 2; i++)
      {
        int min_off = (std::min)((std::min)(off0[i], off1[i]), off2[i]);
        if (min_off != 0)
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

    int o0i = ((off0[0] & 1) << 1) + (off0[1] & 1);
    int o1i = ((off1[0] & 1) << 1) + (off1[1] & 1);
    int o2i = ((off2[0] & 1) << 1) + (off2[1] & 1);
    f->set_offsets(o0i, o1i, o2i);
  }
  //\}

  /// Checks the too_long_edges bookkeeping
  bool is_valid_too_long_edges(bool verbose = false, int level = 0) const;

  /** @name Checking helpers */ //@{
  /// calls has_self_edges for every face of the triangulation
  bool has_self_edges() const
  {
    Face_iterator it;
    for ( it = all_faces_begin(); it != all_faces_end(); ++it )
      if (has_self_edges(it)) return true;
    return false;
  }
  bool has_self_edges(Face_handle fh) const
  {
    CGAL_triangulation_assertion((fh->vertex(0) != fh->vertex(1)) ||
                                 (fh->offset(0) != fh->offset(1)));
    CGAL_triangulation_assertion((fh->vertex(0) != fh->vertex(2)) ||
                                 (fh->offset(0) != fh->offset(2)));
    CGAL_triangulation_assertion((fh->vertex(1) != fh->vertex(2)) ||
                                 (fh->offset(1) != fh->offset(2)));
    return ((fh->vertex(0) == fh->vertex(1)) ||
            (fh->vertex(0) == fh->vertex(2)) ||
            (fh->vertex(1) == fh->vertex(2)));
  }

  //@}


protected:
  // Protected data of Periodic_2_triangulation_2
  /// \name Triangulation data members
  // \{

  /// Geometric traits
  Gt _gt;
  /// Triangulation data structure
  Tds _tds;
  // \}

  /// Returns false, no infinite simplices in the periodic triangulation
  template <class T>
  inline bool is_infinite(T) const
  {
    return false;
  }

private:
  /// Inserts (p,o) in the face f and sets the offsets of the newly created faces
  /// Doesn't insert periodic copies
  Vertex_handle insert_in_face(const Point& p, const Offset &o,
                               Face_handle f,
                               Vertex_handle vh);
  /// Inserts (p,o) in the edge (f,i) and sets the offsets of the newly created faces
  /// Doesn't insert periodic copies
  Vertex_handle insert_in_edge(const Point& p, const Offset &o,
                               Face_handle f, int i,
                               Vertex_handle vh);

  /// Remove a vertex without removing it's possible periodic copies.
  /// Helper functions
  void remove_degree_3_single_copy(Vertex_handle vh);

  // Private data of Periodic_2_triangulation_2
  /// \name Periodic members
  //\{

  /// Determines if we currently compute in 3-cover or 1-cover.
  Covering_sheets _cover;

  /// The domain
  Iso_rectangle _domain;

  /// This threshold should be chosen such that if all edges are shorter,
  /// we can be sure that there are no self-edges anymore.
  FT _edge_length_threshold;

  /// This adjacency list stores all edges that are longer than
  /// edge_length_threshold.
  Too_long_edges_map _too_long_edges;
  /// Number of edges that are too long
  size_t _too_long_edge_counter;

  /// map of offsets for periodic copies of vertices
  Virtual_vertex_map _virtual_vertices;
  /// map of a non-virtual vertex to its virtual copies
  Virtual_vertex_reverse_map _virtual_vertices_reverse;
  //\}
}; // class Periodic_2_triangulation_2

// CONSTRUCTORS
template<class Gt, class Tds>
Periodic_2_triangulation_2<Gt, Tds>::Periodic_2_triangulation_2(
  const Iso_rectangle & domain, const Geom_traits& geom_traits)
  : _gt(geom_traits), _tds()
  , _cover(make_array(1, 1))
  , _domain(domain)
  , _too_long_edge_counter(0)
{
  CGAL_triangulation_precondition(_domain.xmax() - _domain.xmin() ==
                                  _domain.ymax() - _domain.ymin());
  set_domain(_domain);
}

// copy constructor duplicates vertices and faces
template<class Gt, class Tds>
Periodic_2_triangulation_2<Gt, Tds>::Periodic_2_triangulation_2(const Periodic_2_triangulation_2 &tr)
{
  copy_triangulation(tr);
}

//Assignment
template<class Gt, class Tds>
Periodic_2_triangulation_2<Gt, Tds> &
Periodic_2_triangulation_2<Gt, Tds>::operator=(
  const Periodic_2_triangulation_2 &tr)
{
  copy_triangulation(tr);
  return *this;
}

// Helping functions
template < class GT, class Tds >
class Periodic_2_triangulation_2<GT, Tds>::Finder
{
  const Self* _t;
  const Point & _p;
public:
  Finder(const Self* t, const Point &p) : _t(t), _p(p) {}
  bool operator()(const Vertex_handle v)
  {
    return _t->xy_equal(v->point(), _p);
  }
};

template < class GT, class Tds >
inline void
Periodic_2_triangulation_2<GT, Tds>::
copy_multiple_covering(const Periodic_2_triangulation_2<GT, Tds> & tr)
{
  // Write the respective offsets in the vertices to make them
  // automatically copy with the tds.
  for (Vertex_iterator vit = tr.vertices_begin() ;
       vit != tr.vertices_end() ; ++vit)
    {
      vit->set_offset(tr.get_offset(vit));
    }
  // copy the tds
  _tds = tr.tds();
  // make a list of all vertices that belong to the original
  // domain and initialize the basic structure of
  // virtual_vertices_reverse
  std::list<Vertex_handle> vlist;
  for (Vertex_iterator vit = vertices_begin() ;
       vit != vertices_end() ; ++vit)
    {
      if (vit->offset() == Offset())
        {
          vlist.push_back(vit);
          _virtual_vertices_reverse.insert(
            std::make_pair(vit, std::vector<Vertex_handle>(8)));
          CGAL_triangulation_assertion(_virtual_vertices_reverse.find(vit)
                                       ->second.size() == 8);
        }
    }
  // Iterate over all vertices that are not in the original domain
  // and construct the respective entries to virtual_vertices and
  // virtual_vertices_reverse
  for (Vertex_iterator vit2 = vertices_begin() ;
       vit2 != vertices_end() ; ++vit2)
    {
      if (vit2->offset() != Offset())
        {
          //TODO: use some binding, maybe boost instead of the Finder.
          typename std::list<Vertex_handle>::iterator vlist_it
          = std::find_if(vlist.begin(), vlist.end(),
                         Finder(this, vit2->point()));
          Offset off = vit2->offset();
          _virtual_vertices.insert(std::make_pair(vit2,
                                                  std::make_pair(*vlist_it, off)));
          _virtual_vertices_reverse.find(*vlist_it)
          ->second[3 * off[0] + off[1] - 1] = vit2;
          CGAL_triangulation_assertion(get_offset(vit2) == off);
        }
    }
  // Cleanup vertex offsets
  for (Vertex_iterator vit = vertices_begin() ;
       vit != vertices_end() ; ++vit)
    vit->clear_offset();
  for (Vertex_iterator vit = tr.vertices_begin() ;
       vit != tr.vertices_end() ; ++vit)
    vit->clear_offset();
  // Build up the too_long_edges container
  _too_long_edge_counter = 0;
  _too_long_edges.clear();
  for (Vertex_iterator vit = vertices_begin() ;
       vit != vertices_end() ; ++vit)
    _too_long_edges[vit] = std::list<Vertex_handle>();
  std::pair<Vertex_handle, Vertex_handle> edge_to_add;
  Point p1, p2;
  int i, j;
  for (Edge_iterator eit = edges_begin() ;
       eit != edges_end() ; ++eit)
    {
      if (&*(eit->first->vertex(cw(eit->second)))
          < &*(eit->first->vertex(ccw(eit->second))))
        {
          i = cw(eit->second);
          j = ccw(eit->second);
        }
      else
        {
          i = ccw(eit->second);
          j = cw(eit->second);
        }
      edge_to_add = std::make_pair(eit->first->vertex(i),
                                   eit->first->vertex(j));
      p1 = construct_point(eit->first->vertex(i)->point(),
                           get_offset(eit->first, i));
      p2 = construct_point(eit->first->vertex(j)->point(),
                           get_offset(eit->first, j));
      Vertex_handle v_no = eit->first->vertex(i);
      if (squared_distance(p1, p2) > _edge_length_threshold)
        {
          CGAL_triangulation_assertion(
            find(_too_long_edges[v_no].begin(),
                 _too_long_edges[v_no].end(),
                 edge_to_add.second) == _too_long_edges[v_no].end());
          _too_long_edges[v_no].push_back(edge_to_add.second);
          _too_long_edge_counter++;
        }
    }
}

template<class Gt, class Tds>
void Periodic_2_triangulation_2<Gt, Tds>::copy_triangulation(
  const Periodic_2_triangulation_2 &tr)
{
  _tds.clear();
  _gt = tr._gt;
  _cover = tr._cover;
  _domain = tr._domain;
  _edge_length_threshold = tr._edge_length_threshold;
  _too_long_edge_counter = tr._too_long_edge_counter;
  if (tr.is_1_cover())
    {
      _tds = tr.tds();
    }
  else
    {
      copy_multiple_covering(tr);
    }
  CGAL_assertion(_too_long_edge_counter == tr._too_long_edge_counter);
  CGAL_triangulation_expensive_postcondition(*this == tr);
}

template<class Gt, class Tds>
void Periodic_2_triangulation_2<Gt, Tds>::swap(Periodic_2_triangulation_2 &tr)
{
  _tds.swap(tr._tds);

  Geom_traits t = geom_traits();
  _gt = tr.geom_traits();
  tr._gt = t;

  std::swap(tr._cover, _cover);
  std::swap(tr._domain, _domain);

  std::swap(tr._edge_length_threshold, _edge_length_threshold);
  std::swap(tr._too_long_edges, _too_long_edges);
  std::swap(tr._too_long_edge_counter, _too_long_edge_counter);

  std::swap(tr._virtual_vertices, _virtual_vertices);
  std::swap(tr._virtual_vertices_reverse, _virtual_vertices_reverse);
}

template<class Gt, class Tds>
void Periodic_2_triangulation_2<Gt, Tds>::clear()
{
  _tds.clear();
  _tds.set_dimension(-2);

  _too_long_edges.clear();
  _too_long_edge_counter = 0;

  _virtual_vertices.clear();
  _virtual_vertices_reverse.clear();

  _cover = make_array(1, 1);
}

template<class Gt, class Tds>
bool Periodic_2_triangulation_2<Gt, Tds>::is_valid(Face_handle fh, bool /*verbose*/, int /*level*/) const
{
  bool result = true;

  int xmin, xmax, ymin, ymax;
  xmin = ymin = 3;
  xmax = ymax = 0;
  for (int i = 0; i < 3; ++i)
    {
      Offset o = get_offset(fh, i);
      xmin = (std::min)(xmin, o[0]);
      xmax = (std::max)(xmax, o[0]);
      ymin = (std::min)(ymin, o[1]);
      ymax = (std::max)(ymax, o[1]);
    }
  // Should at most cross 1 border in each direction
  result &= (xmax - xmin <= 1);
  result &= (ymax - ymin <= 1);
  if (!result)
    {
      std::cerr << "min/max: " << xmin << "," << xmax << " " << ymin << "," << ymax << std::endl;
      for (int i = 0; i < 3; ++i)
        {
          Offset o = get_offset(fh, i);
          std::cerr << "Offset: " << o << std::endl;
        }
      std::cerr << std::endl;
      CGAL_triangulation_assertion(false);
    }

  return result;
}

template<class Gt, class Tds>
bool Periodic_2_triangulation_2<Gt, Tds>::is_valid(bool verbose, int level) const
{
  bool result = _tds.is_valid(verbose, level);
  CGAL_triangulation_assertion(result);

  if (dimension() == 2)
    {
      // Check positive orientation:
      const Point *p[3];
      Offset off[3];
      for (Face_iterator fit = faces_begin(); fit != faces_end(); ++fit)
        {
          for (int i = 0; i < 3; i++)
            {
              p[i] = &fit->vertex(i)->point();
              off[i] = get_offset(fit, i);
            }

          if (orientation(*p[0], *p[1], *p[2], off[0], off[1], off[2]) != POSITIVE)
            {
              if (verbose)
                {
                  std::cerr
                      << "Periodic_2_triangulation_2: wrong orientation:" << "\n"
                      << *p[0] << " \t" << off[0] << "\n"
                      << *p[1] << " \t" << off[1] << "\n"
                      << *p[2] << " \t" << off[2] << std::endl;
                }
              result = false;
            }
        }
    }
  CGAL_triangulation_assertion(result);

  // Check for the right number of simplices
  int copies = number_of_sheets()[0] * number_of_sheets()[1];
  result &= (number_of_stored_vertices() == copies * number_of_vertices());
  result &= (number_of_stored_edges() == copies * number_of_edges());
  result &= (number_of_stored_faces() == copies * number_of_faces());
  CGAL_triangulation_assertion(result);

  // check number of euler characteristic. This cannot be done by the Tds
  // which does not know the genus
  result &= (number_of_stored_vertices() - number_of_stored_edges()
             + number_of_stored_faces() == 0);
  CGAL_triangulation_assertion(result);

  result &= !has_self_edges();
  CGAL_triangulation_assertion(result);

  // Edges should not be longer than 1 periodicity
  for (Face_iterator fit = faces_begin(); fit != faces_end(); ++fit)
    {
      result &= is_valid(fit, verbose, level);
    }
  CGAL_triangulation_assertion(result);

  result &= is_1_cover() == _virtual_vertices.empty();
  result &= is_1_cover() == _virtual_vertices_reverse.empty();
  result &= (_virtual_vertices.size() == (number_of_sheets()[0]
                                          * number_of_sheets()[1] - 1) * _virtual_vertices_reverse.size());
  CGAL_triangulation_assertion(result);

  for (Virtual_vertex_map_it it = _virtual_vertices.begin(); it
       != _virtual_vertices.end(); ++it)
    {
      const Vertex_handle &copy = it->first;
      const Vertex_handle &orig = it->second.first;
      const Offset &off = it->second.second;
      size_t index = number_of_sheets()[0] * off[0] + off[1] - 1;
      Virtual_vertex_reverse_map_it rev_it = _virtual_vertices_reverse.find(orig);
      if (rev_it != _virtual_vertices_reverse.end())
        {
          if (index < rev_it->second.size())
            {
              result &= (rev_it->second[index] == copy);
            }
          else
            {
              result &= false;
            }
        }
      else
        {
          result &= false;
        }
    }
  CGAL_triangulation_assertion(result);

  for (Virtual_vertex_reverse_map_it it = _virtual_vertices_reverse.begin(); it
       != _virtual_vertices_reverse.end(); ++it)
    {
      const std::vector<Vertex_handle> &copies = it->second;
      result &= copies.size() == 8;
      for (size_t i = 0; i < copies.size(); ++i)
        {
          Virtual_vertex_map_it copy_it = _virtual_vertices.find(copies[i]);
          if (copy_it != _virtual_vertices.end())
            {
              result &= copy_it->second.first == it->first;
            }
          else
            {
              result &= false;
            }
        }
    }

  // Check the too_long_edges administration
  result &= is_valid_too_long_edges(verbose, level);

  return result;
}

template<class Gt, class Tds>
bool Periodic_2_triangulation_2<Gt, Tds>::is_valid_too_long_edges(bool verbose, int /*level*/) const
{
  bool result = true;

  result &= is_1_cover() == _too_long_edges.empty();
  CGAL_triangulation_assertion(result);
  size_t too_long_edges = 0;
  for (Too_long_edges_map_it it = _too_long_edges.begin(); it
       != _too_long_edges.end(); ++it)
    {
      too_long_edges += it->second.size();
    }
  CGAL_triangulation_assertion(result);
  if (_too_long_edge_counter != too_long_edges)
    {
      if (verbose) std::cout << "Too long edge counter is incorrect: " << _too_long_edge_counter << " != " << too_long_edges << std::endl;
      result = false;
    }
  CGAL_triangulation_assertion(result);

  /// Expensive check whether the right too long edges are in the list
  if (is_1_cover())
    {
      for (Edge_iterator eit = edges_begin(); eit != edges_end(); ++eit)
        {
          Vertex_handle vh1 = eit->first->vertex(ccw(eit->second));
          Vertex_handle vh2 = eit->first->vertex(cw(eit->second));
          Point p1 = construct_point(vh1->point(), get_offset(eit->first, ccw(eit->second)));
          Point p2 = construct_point(vh2->point(), get_offset(eit->first, cw(eit->second)));
          result &= (!edge_is_too_long(p1, p2));
        }
      CGAL_triangulation_assertion(result);
    }
  else
    {
      too_long_edges = 0;
      for (Edge_iterator eit = edges_begin(); eit != edges_end(); ++eit)
        {
          Vertex_handle vh1 = eit->first->vertex(ccw(eit->second));
          Vertex_handle vh2 = eit->first->vertex(cw(eit->second));
          Point p1 = construct_point(vh1->point(),
                                     get_offset(eit->first, ccw(eit->second)));
          Point p2 = construct_point(vh2->point(),
                                     get_offset(eit->first, cw(eit->second)));

          if (&*vh2 < &*vh1)
            std::swap(vh1, vh2);
          CGAL_triangulation_assertion(&*vh1 < &*vh2);

          bool too_long = edge_is_too_long(p1, p2);
          if (too_long != edge_is_too_long(p2, p1))
            {
              if (verbose) std::cout << "Long edge criterion not symmetric c(v1,v2) != c(v2,v1)" << std::endl;
              result = false;
            }
          CGAL_triangulation_assertion(result);

          Too_long_edges_map_it it = _too_long_edges.find(vh1);
          if (it == _too_long_edges.end())
            {
              if (too_long)
                {
                  if (verbose) std::cout << "1. Too long edge not in the datastructure" << std::endl;
                  result = false;
                }
              result &= !too_long;
              CGAL_triangulation_assertion(result);
            }
          else
            {
              typename std::list<Vertex_handle>::const_iterator it2 = find(it->second.begin(), it->second.end(), vh2);
              if (too_long)
                {
                  too_long_edges++;
                  if (it2 == it->second.end())
                    {
                      if (verbose) std::cout << "2. Too long edge not in the datastructure" << std::endl;
                      result = false;
                    }
                  CGAL_triangulation_assertion(result);
                }
              else
                {
                  if (it2 != it->second.end())
                    {
                      if (verbose) std::cout << "Edge is not too long, but contained in the datastructure" << std::endl;
                      result = false;
                    }
                  CGAL_triangulation_assertion(result);
                }
            }
        }

      if (_too_long_edge_counter != too_long_edges)
        {
          if (verbose)
            std::cout << "Counts do not match: " << _too_long_edge_counter << " != " << too_long_edges << std::endl;
          result = false;
        }
      CGAL_triangulation_assertion(result);
    }

  return result;
}

template<class Gt, class Tds>
bool Periodic_2_triangulation_2<Gt, Tds>::flippable(Face_handle f, int i)
{
  Face_handle nb = f->neighbor(i);
  int j = nb->index(f);

  const Point *p[4];

  p[0] = &f->vertex(i)->point();      // i
  p[1] = &nb->vertex(j)->point();     // opposite
  p[2] = &f->vertex(ccw(i))->point(); // ccw
  p[3] = &f->vertex(cw(i))->point();  // cw

  if (is_1_cover() && f->has_zero_offsets() && nb->has_zero_offsets())
    {
      // if (orientation(*p[0], *p[1], *p[2]) != RIGHT_TURN)
      //   return false;
      // if (orientation(*p[0], *p[1], *p[3]) != LEFT_TURN)
      //   return false;
      if (orientation(*p[0], *p[1], *p[2]) == LEFT_TURN)
        return false;
      if (orientation(*p[0], *p[1], *p[3]) == RIGHT_TURN)
        return false;
    }
  else
    {
      Offset off[4];
      off[0] = get_offset(f, i);
      off[1] = combine_offsets(get_offset(nb, j), get_neighbor_offset(nb, j, f, i));
      off[2] = get_offset(f, ccw(i));
      off[3] = get_offset(f, cw(i));

      // if (orientation(*p[0], *p[1], *p[2], off[0], off[1], off[2]) != RIGHT_TURN)
      //   return false;
      // if (orientation(*p[0], *p[1], *p[3], off[0], off[1], off[3]) != LEFT_TURN)
      //   return false;
      if (orientation(*p[0], *p[1], *p[2], off[0], off[1], off[2]) == LEFT_TURN)
        return false;
      if (orientation(*p[0], *p[1], *p[3], off[0], off[1], off[3]) == RIGHT_TURN)
        return false;
    }

  return true;
}

template<class Gt, class Tds>
void Periodic_2_triangulation_2<Gt, Tds>::flip(Face_handle f, int i)
{
  if (is_1_cover())
    {
      flip_single_edge(f, i);
      return;
    }

  Vertex_handle vh1 = f->vertex( cw(i));
  Vertex_handle vh2 = f->vertex(ccw(i));
  Virtual_vertex_map_it it_vh1 = _virtual_vertices.find(vh1);
  Virtual_vertex_map_it it_vh2 = _virtual_vertices.find(vh2);

  Offset vh1_offset, vh2_offset;
  if (it_vh1 != _virtual_vertices.end())
    {
      vh1 = it_vh1->second.first;
      vh1_offset = it_vh1->second.second;
    }
  if (it_vh2 != _virtual_vertices.end())
    {
      vh2 = it_vh2->second.first;
      vh2_offset = it_vh2->second.second;
    }

  CGAL_triangulation_assertion(
    virtual_vertices_reverse().find(vh1) != virtual_vertices_reverse().end());
  CGAL_triangulation_assertion(
    virtual_vertices_reverse().find(vh2) != virtual_vertices_reverse().end());
  const std::vector<Vertex_handle> &v1s =
    virtual_vertices_reverse().find(vh1)->second;
  const std::vector<Vertex_handle> &v2s =
    virtual_vertices_reverse().find(vh2)->second;

  CGAL_assertion(v1s.size() == 8);
  CGAL_assertion(v1s.size() == v2s.size());

  Face_handle fh;
  int index=0;
  Vertex_handle vh1_copy, vh2_copy;

  // Virtual copies
  for (int x = 0; x < 3; ++x)
    {
      for (int y = 0; y < 3; ++y)
        {
          int i1 = 3 * ((x + vh1_offset.x()) % 3) + ((y + vh1_offset.y()) % 3);
          int i2 = 3 * ((x + vh2_offset.x()) % 3) + ((y + vh2_offset.y()) % 3);

          if (i1 == 0)
            vh1_copy = vh1;
          else
            vh1_copy = v1s[i1 - 1];
          if (i2 == 0)
            vh2_copy = vh2;
          else
            vh2_copy = v2s[i2 - 1];

          bool found = is_edge(vh1_copy, vh2_copy, fh, index);
	  CGAL_USE(found);
          CGAL_assertion(found);
	  if (found)
	    flip_single_edge(fh, index);
        }
    }

  try_to_convert_to_one_cover();
}
template<class Gt, class Tds>
void Periodic_2_triangulation_2<Gt, Tds>::flip_single_edge(Face_handle f, int i)
{
  CGAL_triangulation_precondition(f != Face_handle());
  CGAL_triangulation_precondition(i == 0 || i == 1 || i == 2);
  CGAL_triangulation_precondition(dimension() == 2);

  CGAL_triangulation_precondition(flippable(f, i));

  if (!is_1_cover())
    remove_too_long_edge(f, i);

  Face_handle nb = f->neighbor(i);
  if (f->has_zero_offsets() && nb->has_zero_offsets())
    {
      _tds.flip(f, i);

      if (!is_1_cover())
        insert_too_long_edge(f, ccw(i));

      return;
    }

  int nb_index = nb->index(f);
  int offsets[4];
  offsets[0] = f->offset(i);
  offsets[1] = f->offset(cw(i));
  offsets[2] = f->offset(ccw(i));
  offsets[3] = nb->offset(nb_index);

  // Move the offsets of f and nb in the same space by correcting for nb_offset
  Offset nb_offset = get_neighbor_offset(f, i, nb, nb_index);
  if (nb_offset.x() != 0)
    {
      if (nb_offset.x() == 1)
        {
          CGAL_assertion(((offsets[0] & 2) | (offsets[1] & 2) | (offsets[2] & 2)) == 0);
          offsets[0] |= 2;
          offsets[1] |= 2;
          offsets[2] |= 2;
        }
      else
        {
          CGAL_triangulation_assertion(nb_offset.x() == -1);
          CGAL_assertion((offsets[3] & 2) == 0);
          offsets[3] |= 2;
        }
    }
  if (nb_offset.y() != 0)
    {
      if (nb_offset.y() == 1)
        {
          CGAL_assertion(((offsets[0] & 1) | (offsets[1] & 1) | (offsets[2] & 1)) == 0);
          offsets[0] |= 1;
          offsets[1] |= 1;
          offsets[2] |= 1;
        }
      else
        {
          CGAL_triangulation_assertion(nb_offset.y() == -1);
          CGAL_assertion((offsets[3] & 1) == 0);
          offsets[3] |= 1;
        }
    }
  CGAL_assertion((offsets[0] & offsets[1] & offsets[2] & offsets[3]) == 0);
  CGAL_triangulation_assertion_code(Vertex_handle vh = f->vertex(i));
  CGAL_triangulation_assertion_code(Vertex_handle vh_ccw = f->vertex(ccw(i)));
  _tds.flip(f, i);
  // Combinatorial checks
  CGAL_triangulation_assertion(vh == f->vertex(i));
  CGAL_triangulation_assertion(vh_ccw == f->vertex(ccw(i)));
  CGAL_triangulation_assertion(f->vertex(i) == nb->vertex(cw(nb_index)));
  CGAL_triangulation_assertion(f->vertex(cw(i)) == nb->vertex(nb_index));

  // Restore the offsets
  int new_off[3];
  // For face f
  new_off[i] = offsets[0];
  new_off[ccw(i)] = offsets[2];
  new_off[cw(i)] = offsets[3];
  set_offsets(f, new_off[0], new_off[1], new_off[2]);
  // For face nb
  new_off[nb_index] = offsets[3];
  new_off[ccw(nb_index)] = offsets[1];
  new_off[cw(nb_index)] = offsets[0];
  set_offsets(nb, new_off[0], new_off[1], new_off[2]);

  if (!is_1_cover())
    insert_too_long_edge(f, ccw(i));
}

template<class Gt, class Tds>
void
Periodic_2_triangulation_2<Gt, Tds>::remove_from_virtual_copies(Vertex_handle v)
{
  typename Virtual_vertex_reverse_map::iterator rev_it = _virtual_vertices_reverse.find(v);
  CGAL_triangulation_assertion(rev_it != _virtual_vertices_reverse.end());

  const std::vector<Vertex_handle> &virtual_copies = rev_it->second;
  for (size_t i = 0; i < virtual_copies.size(); ++i)
    {
      _virtual_vertices.erase(virtual_copies[i]);
    }
  _virtual_vertices_reverse.erase(rev_it);
}

template<class Gt, class Tds>
typename Periodic_2_triangulation_2<Gt, Tds>::Vertex_handle Periodic_2_triangulation_2 <
Gt, Tds >::insert_first(const Point& p)
{
  CGAL_assertion(empty());
  // The empty triangulation has a single sheeted cover
  _cover = make_array(3, 3);

  /// Virtual vertices, one per periodic domain
  Vertex_handle vir_vertices[3][3];
  /// Virtual faces, two per periodic domain
  Face_handle faces[3][3][2];

  // Initialise vertices:
  vir_vertices[0][0] = _tds.create_vertex();
  vir_vertices[0][0]->set_point(p);
  _virtual_vertices_reverse[vir_vertices[0][0]] = std::vector<Vertex_handle>();
  for (int i = 0; i < _cover[0]; i++)
    {
      for (int j = 0; j < _cover[1]; j++)
        {
          if ((i != 0) || (j != 0))
            {
              // Initialise virtual vertices out of the domain for debugging
              vir_vertices[i][j] = _tds.create_vertex();
              vir_vertices[i][j]->set_point(p); //+Offset(i,j));
              _virtual_vertices[vir_vertices[i][j]] = Virtual_vertex(
                  vir_vertices[0][0], Offset(i, j));
              _virtual_vertices_reverse[vir_vertices[0][0]].push_back(
                vir_vertices[i][j]);
            }
        }
    }

  // Create faces:
  for (int i = 0; i < _cover[0]; i++)
    {
      for (int j = 0; j < _cover[1]; j++)
        {
          for (int f = 0; f < 2; f++)
            {
              // f faces per 'rectangle'
              faces[i][j][f] = _tds.create_face();
            }
        }
    }

  // table containing the vertex information
  // index to the right vertex: [number of faces][vertex][offset]
  int vertex_ind[2][3][2] = { { { 0, 0 }, { 1, 1 }, { 0, 1 } }, { { 0, 0 }, {
        1, 0
      }, { 1, 1 }
    }
  };
  // Table containing the neighbor information
  // [number of faces][neighbor][offset,face]
  int neighb_ind[2][3][3] = { { { 0, 1, 1 }, { -1, 0, 1 }, { 0, 0, 1 } }, { {
        1, 0, 0
      }, { 0, 0, 0 }, { 0, -1, 0 }
    }
  };
  for (int i = 0; i < _cover[0]; i++)
    {
      for (int j = 0; j < _cover[1]; j++)
        {
          int offset =
            ((i == _cover[0] - 1 ? 2 : 0) | (j == _cover[1] - 1 ? 1 : 0));
          for (int f = 0; f < 2; f++)
            {
              faces[i][j][f]->set_vertices(vir_vertices[(i + vertex_ind[f][0][0])
                                           % _cover[0]][(j + vertex_ind[f][0][1]) % _cover[1]],
                                           vir_vertices[(i + vertex_ind[f][1][0]) % _cover[0]][(j
                                               + vertex_ind[f][1][1]) % _cover[1]], vir_vertices[(i
                                                   + vertex_ind[f][2][0]) % _cover[0]][(j + vertex_ind[f][2][1])
                                                       % _cover[1]]);
              set_offsets(faces[i][j][f], offset & (vertex_ind[f][0][0] * 2
                                                    + vertex_ind[f][0][1] * 1), offset & (vertex_ind[f][1][0] * 2
                                                        + vertex_ind[f][1][1] * 1), offset & (vertex_ind[f][2][0] * 2
                                                            + vertex_ind[f][2][1] * 1));
              faces[i][j][f]->set_neighbors(faces[(i + _cover[0]
                                                   + neighb_ind[f][0][0]) % _cover[0]][(j + _cover[1]
                                                       + neighb_ind[f][0][1]) % _cover[1]][neighb_ind[f][0][2]], faces[(i
                                                           + _cover[0] + neighb_ind[f][1][0]) % _cover[0]][(j + _cover[1]
                                                               + neighb_ind[f][1][1]) % _cover[1]][neighb_ind[f][1][2]], faces[(i
                                                                   + _cover[0] + neighb_ind[f][2][0]) % _cover[0]][(j + _cover[1]
                                                                       + neighb_ind[f][2][1]) % _cover[1]][neighb_ind[f][2][2]]);
            }
        }
    }
  // set pointers from the vertices to incident faces.
  for (int i = 0; i < _cover[0]; i++)
    {
      for (int j = 0; j < _cover[1]; j++)
        {
          vir_vertices[i][j]->set_face(faces[i][j][0]);
        }
    }

  _tds.set_dimension(2);

  // create the base for too_long_edges;
  CGAL_triangulation_assertion(_too_long_edges.empty() );
  CGAL_triangulation_assertion(_too_long_edge_counter == 0);

  // Insert all vertices as the first vertex in the _too_long_edges list
  int k = 0;
  std::list<Vertex_handle> empty_list;
  for (Vertex_iterator vit = vertices_begin(); vit != vertices_end(); ++vit)
    {
      _too_long_edges[vit] = empty_list;
      k++;
    }

  // Insert all edges as all edges are too long
  _too_long_edge_counter = 0;
  for (Edge_iterator eit = edges_begin(); eit != edges_end(); eit++)
    {
      Vertex_handle vh1 = eit->first->vertex(ccw(eit->second));
      Vertex_handle vh2 = eit->first->vertex(cw(eit->second));
      if (&*vh1 < &*vh2)
        {
          _too_long_edges[vh1].push_back(vh2);
        }
      else
        {
          _too_long_edges[vh2].push_back(vh1);
        }
      _too_long_edge_counter++;
    }

  return vir_vertices[0][0];
}

template<class Gt, class Tds>
typename Periodic_2_triangulation_2<Gt, Tds>::Vertex_handle
Periodic_2_triangulation_2<Gt, Tds>::insert_in_edge(const Point& p,
    Face_handle f, int i)
{
  return insert(p, EDGE, f, i);
}
template<class Gt, class Tds>
typename Periodic_2_triangulation_2<Gt, Tds>::Vertex_handle
Periodic_2_triangulation_2<Gt, Tds>::insert_in_edge(const Point& p, const Offset &o,
    Face_handle f, int i,
    Vertex_handle vh)
{
  // Insert in edge calls an insert_in_face and a flip.
  // Therefore there is no need to update the too_long_edges bookkeeping directly.

  CGAL_triangulation_assertion(number_of_vertices() != 0);
  CGAL_triangulation_assertion((!is_1_cover()) || (o == Offset()));

  // Backup of the neighbor and its relative offset
  Face_handle nb = f->neighbor(i);
  int j = nb->index(f);
  CGAL_triangulation_assertion_code(Offset current_offset = get_location_offset(f, p, o));
  CGAL_triangulation_assertion
  (orientation(f->vertex(cw(i))->point(), p, f->vertex(ccw(i))->point(),
               get_offset(f, cw(i)), combine_offsets(o, current_offset), get_offset(f, ccw(i))) == COLLINEAR &&
   collinear_between(f->vertex(cw(i))->point(), p, f->vertex(ccw(i))->point(),
                     get_offset(f, cw(i)), combine_offsets(o, current_offset), get_offset(f, ccw(i))) );

  /// Insert in the face and flip an edge
  Vertex_handle v = insert_in_face(p, o, f, vh);
  flip_single_edge(nb, j);

  return v;
}

template<class Gt, class Tds>
typename Periodic_2_triangulation_2<Gt, Tds>::Vertex_handle
Periodic_2_triangulation_2<Gt, Tds>::insert_in_face(const Point& p, Face_handle f)
{
  return insert(p, FACE, f, 0);
}
template<class Gt, class Tds>
typename Periodic_2_triangulation_2<Gt, Tds>::Vertex_handle
Periodic_2_triangulation_2<Gt, Tds>::insert_in_face(const Point& p, const Offset &o,
    Face_handle f,
    Vertex_handle vh)
{
  CGAL_triangulation_assertion(f != Face_handle());
  CGAL_triangulation_assertion(number_of_vertices() != 0);
  CGAL_triangulation_assertion((!is_1_cover()) || (o == Offset()));

  const bool simplicity_criterion = f->has_zero_offsets() && o.is_zero();


  Offset current_off;

  // Save the neighbors and the offsets
  Face_handle nb[3];
  int nb_index[3];
  int offsets[3];
  CGAL_triangulation_assertion_code( Vertex_handle vertices[3]; )

  if (!simplicity_criterion)
    {
      // Choose the periodic copy of tester.point() that is inside c.
      current_off = get_location_offset(f, p, o);

      CGAL_triangulation_assertion(oriented_side(f, p, combine_offsets(o, current_off)) != ON_NEGATIVE_SIDE);

      for (int i = 0; i < 3; ++i)
        {
          nb[i] = f->neighbor(i);
          nb_index[i] = nb[i]->index(f);
          offsets[i] = f->offset(i);
          CGAL_triangulation_assertion_code( vertices[i] = f->vertex(i); );
        }
    }

  // Insert the new vertex
  Vertex_handle v = _tds.insert_in_face(f);
  v->set_point(p);

  if (!simplicity_criterion)
    {
      // Update the offsets
      int v_offset = off_to_int(current_off);
      int new_offsets[3];
      for (int i = 0; i < 3; ++i)
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

  if (!is_1_cover())
    {
      // update the book-keeping in case of a periodic copy
      if (vh != Vertex_handle())
        {
          _virtual_vertices[v] = Virtual_vertex(vh, o);
          _virtual_vertices_reverse[vh].push_back(v);
        }

      insert_too_long_edges_in_star(v);
    }

  return v;
}

template<class Gt, class Tds>
typename Periodic_2_triangulation_2<Gt, Tds>::Vertex_handle
Periodic_2_triangulation_2<Gt, Tds>::insert(const Point &p, Face_handle start)
{
  CGAL_triangulation_assertion((_domain.xmin() <= p.x()) &&
                               (p.x() < _domain.xmax()));
  CGAL_triangulation_assertion((_domain.ymin() <= p.y()) &&
                               (p.y() < _domain.ymax()));

  if (number_of_stored_vertices() == 0)
    {
      return insert_first(p);
    }

  if (start == Face_handle())
    {
      start = faces_begin();
    }

  Locate_type lt;
  int li;
  Face_handle loc = locate(p, lt, li, start);

  if (start != Face_handle())
    {
      CGAL_assertion(start->vertex(0) != Vertex_handle());
    }
  return insert(p, lt, loc, li);
}

template<class Gt, class Tds>
typename Periodic_2_triangulation_2<Gt, Tds>::Vertex_handle
Periodic_2_triangulation_2<Gt, Tds>::insert(const Point& p,
    Locate_type lt, Face_handle loc, int li)
{
  if (number_of_stored_vertices() == 0)
    {
      return insert_first(p);
    }

  // vstart is a vertex incident to the Face_handle start that will be used as
  // for creating a start point for the virtual vertices.
  // We use the virtual copies of a vertex incident to loc.
  Vertex_handle vstart;
  if (!is_1_cover())
    {
      Virtual_vertex_map_it vvmit = _virtual_vertices.find(loc->vertex(0));
      if (vvmit == _virtual_vertices.end())
        {
          vstart = loc->vertex(0);
        }
      else
        {
          vstart = vvmit->second.first;
        }

      // vstart should be non-virtual, but should have virtual copies
      CGAL_triangulation_assertion(_virtual_vertices.find(vstart)
                                   == _virtual_vertices.end());
      CGAL_triangulation_assertion(_virtual_vertices_reverse.find(vstart)
                                   != _virtual_vertices_reverse.end());
    }

  Vertex_handle vh = insert(p, Offset(), lt, loc, li, Vertex_handle());

  // Don't add periodic copies if we are on the 1-cover
  if (is_1_cover())
    return vh;

  // Don't continue if the point lies on a vertex as this will break the
  // start_vertices array below.
  if (lt == VERTEX)
    return vh;

  const std::vector<Vertex_handle> &start_vertices =
    _virtual_vertices_reverse.find(vstart)->second;
  CGAL_assertion(start_vertices.size() == size_t(number_of_sheets()[0] * number_of_sheets()[1] - 1));

  for (int i = 0; i < number_of_sheets()[0]; i++)
    {
      for (int j = 0; j < number_of_sheets()[1]; j++)
        {
          if ((i != 0) || (j != 0))
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

template<class Gt, class Tds>
typename Periodic_2_triangulation_2<Gt, Tds>::Vertex_handle Periodic_2_triangulation_2 <
Gt, Tds >::insert(const Point& p, const Offset &o, Locate_type lt,
                  Face_handle loc, int li, Vertex_handle vh)
// insert a point p, whose localization is known (lt, f, i)
{
  Vertex_handle result;
  switch (lt)
    {
    case FACE:
    {
      result = insert_in_face(p, o, loc, vh);
      break;
    }
    case EDGE:
    {
      result = insert_in_edge(p, o, loc, li, vh);
      break;
    }
    case VERTEX:
    {
      // The vertex is a special case, we can return immediately
      CGAL_assertion(vh == Vertex_handle());
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

  if (!is_1_cover() && (vh == Vertex_handle()))
    {
      _virtual_vertices_reverse[result] = std::vector<Vertex_handle>();
    }

  return result;
}

template<class Gt, class Tds>
inline void Periodic_2_triangulation_2<Gt, Tds>::remove_degree_3(Vertex_handle v)
{
  CGAL_assertion(number_of_vertices() > 1);
  CGAL_assertion(degree(v) == 3);

  if (is_1_cover())
    {
      remove_degree_3_single_copy(v);
      return;
    }

  {
    Virtual_vertex_map_it it = _virtual_vertices.find(v);
    if (it != _virtual_vertices.end())
      {
        v = it->second.first;
      }
  }

  remove_too_long_edges_in_star(v);

  typename Virtual_vertex_reverse_map::iterator reverse_it =
    _virtual_vertices_reverse.find(v);
  CGAL_assertion(reverse_it != _virtual_vertices_reverse.end());

  const std::vector<Vertex_handle> &virtual_copies = reverse_it->second;
  for (typename std::vector<Vertex_handle>::const_iterator it = virtual_copies.begin();
       it != virtual_copies.end(); ++it)
    {
      _virtual_vertices.erase(*it);
      remove_degree_3_single_copy(*it);
    }

  _virtual_vertices_reverse.erase(reverse_it);
  remove_degree_3_single_copy(v);

}

template<class Gt, class Tds>
inline void Periodic_2_triangulation_2<Gt, Tds>::remove_degree_3_single_copy(Vertex_handle vh)
{
  Face_handle f = vh->face();
  int i = ccw(f->index(vh));
  Face_handle f2 = f->neighbor(i);
  int j = f2->index(f);
  // Get the offsets in ccw order
  Offset off[3];
  off[i]      = get_offset(f, i);
  off[ccw(i)] = get_offset(f, ccw(i));
  off[cw(i)]  = combine_offsets(get_offset(f2, j), get_neighbor_offset(f2, j, f, i));
  if (off[0].x() < 0 || off[1].x() < 0 || off[2].x() < 0)
    {
      Offset o(number_of_sheets()[0], 0);
      off[0] += o;
      off[1] += o;
      off[2] += o;
    }
  if (off[0].y() < 0 || off[1].y() < 0 || off[2].y() < 0)
    {
      Offset o(0, number_of_sheets()[1]);
      off[0] += o;
      off[1] += o;
      off[2] += o;
    }

  // Remove the vertex, keep face f
  _tds.remove_degree_3(vh, f);

  // Reset the offsets
  set_offsets(f,
              (off[0].x() >= number_of_sheets()[0] ? 2 : 0) + (off[0].y() >= number_of_sheets()[1] ? 1 : 0),
              (off[1].x() >= number_of_sheets()[0] ? 2 : 0) + (off[1].y() >= number_of_sheets()[1] ? 1 : 0),
              (off[2].x() >= number_of_sheets()[0] ? 2 : 0) + (off[2].y() >= number_of_sheets()[1] ? 1 : 0));
}

template<class Gt, class Tds>
inline void Periodic_2_triangulation_2<Gt, Tds>::remove_first(Vertex_handle)
{
  CGAL_assertion(number_of_vertices() == 1);
  clear();
  return;
}


template < class Gt, class Tds >
bool
Periodic_2_triangulation_2<Gt, Tds>::
remove_degree_init(Vertex_handle v, const Offset &v_o,
                   std::vector<Face_handle> &f,
                   std::vector<Vertex_handle> &w,
                   std::vector<Offset> &offset_w,
                   std::vector<int> &i,
                   int &d, int &maxd,
                   bool &simplicity_criterion)
{
  Bbox_2 bbox = v->point().bbox();
  simplicity_criterion = is_1_cover();

  f[0] = v->face();
  d = 0;

  do
    {
      i[d] = f[d]->index(v);
      w[d] = f[d]->vertex( ccw(i[d]) );
      offset_w[d] = get_offset(f[d], ccw(i[d])) - get_offset(f[d], i[d]) + v_o;
      w[d]->set_face( f[d]->neighbor(i[d])); // do no longer bother about set_face
      simplicity_criterion &= (offset_w[d] == offset_w[0]);

      bbox = bbox + this->construct_point(w[d]->point(), offset_w[d]).bbox();

      ++d;
      if ( d == maxd)
        {
          maxd *= 2;
          f.resize(maxd);
          w.resize(maxd);
          offset_w.resize(maxd);
          i.resize(maxd);
        }

      f[d] = f[d - 1]->neighbor( ccw(i[d - 1]) );

    }
  while(f[d] != f[0]);

  return is_1_cover() &&
         this->edge_is_too_long(Point(bbox.xmin(), bbox.ymin()), Point(bbox.xmax(), bbox.ymax()));
}

template<class Gt, class Tds>
void Periodic_2_triangulation_2<Gt, Tds>::make_hole(Vertex_handle v, std::list<Edge> & hole)
{
  remove_too_long_edges_in_star(v);

  std::list<Face_handle> to_delete;

  Face_handle f, fn;
  int i, in;
  Vertex_handle vv;

  Face_circulator fc = incident_faces(v);
  Face_circulator done(fc);
  do
    {
      f = fc;
      fc++;
      i = f->index(v);
      fn = f->neighbor(i);
      in = fn->index(f);
      vv = f->vertex(cw(i));
      if (vv->face() == f)
        vv->set_face(fn);
      vv = f->vertex(ccw(i));
      if (vv->face() == f)
        vv->set_face(fn);
      fn->set_neighbor(in, Face_handle());
      hole.push_back(Edge(fn, in));
      to_delete.push_back(f);
    }
  while (fc != done);

  while (!to_delete.empty())
    {
      delete_face(to_delete.front());
      to_delete.pop_front();
    }
  return;
}

template<class Gt, class Tds>
inline typename Periodic_2_triangulation_2<Gt, Tds>::Face_handle Periodic_2_triangulation_2 <
Gt, Tds >::create_face(Face_handle f1, int i1, Face_handle f2, int i2,
                       Face_handle f3, int i3)
{
  return _tds.create_face(f1, i1, f2, i2, f3, i3);
}

template<class Gt, class Tds>
inline typename Periodic_2_triangulation_2<Gt, Tds>::Face_handle Periodic_2_triangulation_2 <
Gt, Tds >::create_face(Face_handle f1, int i1, Face_handle f2, int i2)
{
  return _tds.create_face(f1, i1, f2, i2);
}

template<class Gt, class Tds>
inline typename Periodic_2_triangulation_2<Gt, Tds>::Face_handle Periodic_2_triangulation_2 <
Gt, Tds >::create_face(Face_handle f, int i, Vertex_handle v)
{
  return _tds.create_face(f, i, v);
}

template<class Gt, class Tds>
inline typename Periodic_2_triangulation_2<Gt, Tds>::Face_handle Periodic_2_triangulation_2 <
Gt, Tds >::create_face(Vertex_handle v1, Vertex_handle v2, Vertex_handle v3)
{
  return _tds.create_face(v1, v2, v3);
}

template<class Gt, class Tds>
inline typename Periodic_2_triangulation_2<Gt, Tds>::Face_handle Periodic_2_triangulation_2 <
Gt, Tds >::create_face(Vertex_handle v1, Vertex_handle v2, Vertex_handle v3,
                       Face_handle f1, Face_handle f2, Face_handle f3)
{
  return _tds.create_face(v1, v2, v3, f1, f2, f3);
}

template<class Gt, class Tds>
inline typename Periodic_2_triangulation_2<Gt, Tds>::Face_handle Periodic_2_triangulation_2 <
Gt, Tds >::create_face(Face_handle fh)
{
  return _tds.create_face(fh);
}

template<class Gt, class Tds>
inline typename Periodic_2_triangulation_2<Gt, Tds>::Face_handle Periodic_2_triangulation_2 <
Gt, Tds >::create_face()
{
  return _tds.create_face();
}

template<class Gt, class Tds>
inline
void Periodic_2_triangulation_2<Gt, Tds>::delete_face(Face_handle f)
{
  _tds.delete_face(f);
}

template<class Gt, class Tds>
inline
void Periodic_2_triangulation_2<Gt, Tds>::delete_vertex(Vertex_handle v)
{
  _tds.delete_vertex(v);
}

template<class Gt, class Tds>
bool Periodic_2_triangulation_2<Gt, Tds>::compare_walks(const Point& p,
    Face_handle c1, Face_handle c2, Locate_type& lt1, Locate_type& lt2,
    int li1, int li2) const
{
  bool b = true;
  b &= (lt1 == lt2);
  if ((lt1 == lt2) && (lt1 == VERTEX))
    {
      b &= (c1->vertex(li1) == c2->vertex(li2));
    }
  else if ((lt1 == lt2) && (lt1 == EDGE))
    {
      b &= ((c1 == c2)
            || ((c1->neighbor(li1) == c2) && (c2->neighbor(li2) == c1)));
    }
  else if ((lt1 == lt2) && (lt1 == EMPTY))
    {
      // Skip
    }
  else
    {
      b &= (lt1 == lt2);
      b &= (lt1 == FACE);
      b &= (c1 == c2);
    }

  if (!b)
    {
      std::cerr << "from compare_walks " << std::endl;
      std::cerr << "point " << p << std::endl;
      std::cerr << "locate 1 " << &*c1 << "\t" << lt1 << "\t" << li1 << std::endl;
      std::cerr << "locate 2 " << &*c2 << "\t" << lt2 << "\t" << li2 << std::endl;
      std::cerr << std::endl;
      show_face(c1);
      std::cerr << std::endl;
      show_face(c2);
      std::cerr << std::endl;
    }

  CGAL_triangulation_assertion(b);
  return b;
}

template<class Gt, class Tds>
typename Periodic_2_triangulation_2<Gt, Tds>::Face_handle
Periodic_2_triangulation_2<Gt, Tds>::
march_locate_2D(Face_handle f, const Point& query,
                const Offset& o_p, Locate_type& lt, int& li) const
{
  CGAL_assertion(!empty());

  Offset off_query = o_p;

  // Random generator
  boost::rand48 rng;
  boost::uniform_smallint<> two(0, 1);
  boost::variate_generator<boost::rand48&, boost::uniform_smallint<> > coin(rng, two);

  // Give the point the best start-offset possible
  if (is_1_cover() && !f->has_zero_offsets())
    {
      int cumm_off = f->offset(0) | f->offset(1) | f->offset(2);
      if (((cumm_off & 2) == 2) &&
          (FT(2) * query.x() < (_domain.xmax() + _domain.xmin())))
        off_query += Offset(1, 0);
      if (((cumm_off & 1) == 1) &&
          (FT(2) * query.y() < (_domain.ymax() + _domain.ymin())))
        off_query += Offset(0, 1);
    }

  Face_handle prev = Face_handle();
  int prev_index = 0;
  Offset off[3];
  Orientation o[3];
  while (1)
    {
      // Instead of testing its edges in a random order we do the following
      // until we find a neighbor to go further:
      // As we come from prev we do not have to check the edge leading to prev
      // Now we flip a coin in order to decide if we start checking the
      // edge before or the edge after the edge leading to prev
      int left_first = coin() % 2;

      bool simplicity_criterion =
        f->has_zero_offsets() && off_query.is_null() && is_1_cover();

      const Point *p[3] =
      {
        &f->vertex(0)->point(),
        &f->vertex(1)->point(),
        &f->vertex(2)->point()
      };

      // Get the offsets
      if (!simplicity_criterion)
        {
          if (!is_1_cover())
            {
              // Just fetch the vertices of c as points with offsets
              for (int i = 0; i < 3; i++)
                {
                  off[i] = get_offset(f, i);
                }
            }
          else
            {
              // We are on the one cover and on the boundary between domains
              // Hence, we need to check predicates with offsets
              for (int i = 0; i < 3; i++)
                {
                  off[i] = int_to_off(f->offset(i));
                }
            }
        }

      if (prev == Face_handle())
        {
          prev = f;
          // First step, also check the prev_index
          if (simplicity_criterion)
            {
              o[ccw(prev_index)] =
                orientation(*p[ccw(prev_index)], *p[cw(prev_index)], query);
            }
          else
            {
              o[ccw(prev_index)] =
                orientation(*p[ccw(prev_index)], *p[cw(prev_index)], query,
                            off[ccw(prev_index)], off[cw(prev_index)], off_query);
            }
          if (o[ccw(prev_index)] == NEGATIVE)
            {
              // This assignment is already done: prev = f
              f = f->neighbor(prev_index);
              int new_index = f->index(prev);
              if (!(simplicity_criterion && f->has_zero_offsets()))
                off_query = combine_offsets(off_query,
                                            get_neighbor_offset(prev, prev_index,
                                                f, new_index));
              prev_index = new_index;
              continue;
            }
        }
      else
        {
          o[ccw(prev_index)] = POSITIVE;
        }

      if (left_first)
        {
          if (simplicity_criterion)
            {
              o[prev_index] =
                orientation(*p[prev_index], *p[ccw(prev_index)], query);
            }
          else
            {
              o[prev_index] =
                orientation(*p[prev_index], *p[ccw(prev_index)], query,
                            off[prev_index], off[ccw(prev_index)], off_query);
            }
          if (o[prev_index] == NEGATIVE)
            {
              prev = f;
              f = f->neighbor(cw(prev_index));
              int new_index = f->index(prev);
              if (!(simplicity_criterion && f->has_zero_offsets()))
                off_query = combine_offsets(off_query,
                                            get_neighbor_offset(prev, cw(prev_index), f, new_index));
              prev_index = new_index;
              continue;
            }
        }
      {
        // Do right side
        if (simplicity_criterion)
          {
            o[cw(prev_index)] =
              orientation(*p[cw(prev_index)], *p[prev_index], query);
          }
        else
          {
            o[cw(prev_index)] =
              orientation(*p[cw(prev_index)], *p[prev_index], query,
                          off[cw(prev_index)], off[prev_index], off_query);
          }
        if (o[cw(prev_index)] == NEGATIVE)
          {
            prev = f;
            f = f->neighbor(ccw(prev_index));
            int new_index = f->index(prev);
            if (!(simplicity_criterion && f->has_zero_offsets()))
              off_query = combine_offsets(off_query,
                                          get_neighbor_offset(prev, ccw(prev_index), f, new_index));
            prev_index = new_index;
            continue;
          }
      }
      if (!left_first)
        {
          if (simplicity_criterion)
            {
              o[prev_index] = orientation(*p[prev_index], *p[ccw(prev_index)], query);
            }
          else
            {
              o[prev_index] = orientation(*p[prev_index], *p[ccw(prev_index)], query,
                                          off[prev_index], off[ccw(prev_index)], off_query);
            }
          if (o[prev_index] == NEGATIVE)
            {
              prev = f;
              f = f->neighbor(cw(prev_index));
              int new_index = f->index(prev);
              if (!(simplicity_criterion && f->has_zero_offsets()))
                off_query = combine_offsets(off_query,
                                            get_neighbor_offset(prev, cw(prev_index), f, new_index));
              prev_index = new_index;
              continue;
            }
        }

      // now p is in c or on its boundary
      int sum = (o[0] == COLLINEAR) + (o[1] == COLLINEAR) + (o[2] == COLLINEAR);
      switch (sum)
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

template<class Gt, class Tds>
typename Periodic_2_triangulation_2<Gt, Tds>::Face_handle Periodic_2_triangulation_2 <
Gt, Tds >::locate(const Point& p, const Offset &o, Locate_type& lt, int& li,
                  Face_handle start) const
{
  CGAL_triangulation_assertion((_domain.xmin() <= p.x()) &&
                               (p.x() < _domain.xmax()));
  CGAL_triangulation_assertion((_domain.ymin() <= p.y()) &&
                               (p.y() < _domain.ymax()));

  if (dimension() <= 0)
    {
      lt = EMPTY;
      li = 4;
      return Face_handle();
    }

  // Triangulation is not empty
  if (start == Face_handle())
    {
      start = faces_begin();
    }

  return march_locate_2D(start, p, o, lt, li);
}

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
void Periodic_2_triangulation_2<Gt, Tds>::convert_to_1_sheeted_covering()
{
  // ###################################################################
  // ### First face iteration ##########################################
  // ###################################################################
  {
    if (is_1_cover())
      return;
    CGAL_triangulation_expensive_assertion(is_triangulation_in_1_sheet());

    bool to_delete, has_simplifiable_offset;
    Virtual_vertex_map_it vvmit;
    // First iteration over all faces: Mark the faces that are to delete.
    // Faces are to delete if they cannot be translated anymore in the
    // direction of one of the axes without yielding negative offsets.
    for (Face_iterator it = all_faces_begin(); it != all_faces_end(); ++it)
      {
        to_delete = false;
        // for all directions in 2D Space
        for (int j = 0; j < 2; j++)
          {
            has_simplifiable_offset = true;
            // for all vertices of face it
            for (int i = 0; i < 3; i++)
              {
                vvmit = _virtual_vertices.find(it->vertex(i));
                if (vvmit == _virtual_vertices.end())
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
            if (has_simplifiable_offset)
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
    for (Face_iterator it = all_faces_begin(); it != all_faces_end(); ++it)
      {
        // Skip all faces that are to delete.
        if (it->get_additional_flag() == 1)
          continue;

        // Redirect neighbors: Only neighbors that are marked by the
        // additional_flag have to be substituted by one of their periodic
        // copies. The unmarked neighbors stay the same.
        for (int i = 0; i < 3; i++)
          {
            if (it->neighbor(i)->get_additional_flag() != 1)
              continue;

            nb = it->neighbor(i);

            for (int j = 0; j < 3; j++)
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
            CGAL_triangulation_assertion( !difference_offset.is_null() );

            // We now have to find the "original" periodic copy of nb from
            // its vertices. Therefore, we first have to find the vertices.
            for (int j = 0; j < 3; j++)
              {
                CGAL_triangulation_assertion( (off[j] - difference_offset)[0] >= 0);
                CGAL_triangulation_assertion( (off[j] - difference_offset)[1] >= 0);
                CGAL_triangulation_assertion( (off[j] - difference_offset)[0] < 3);
                CGAL_triangulation_assertion( (off[j] - difference_offset)[1] < 3);

                // find the Vertex_handles of the vertices of the "original"
                // periodic copy of nb. If the vertex is inside the original
                // domain, there is nothing to do
                if ((off[j] - difference_offset).is_null())
                  {
                    nbv[j] = vert[j];
                    // If the vertex is outside the original domain, we have to search
                    // in _virtual_vertices in the "wrong" direction. That means, we
                    // cannot use _virtual_vertices.find but have to use
                    // _virtual_vertices_reverse.
                  }
                else
                  {
                    Offset nbo = off[j] - difference_offset;
                    nbv[j] = _virtual_vertices_reverse.find(vert[j]) ->second[nbo[0]
                             * 3 + nbo[1] - 1];
                  }
              }
            // Find the new neighbor by its 4 vertices
            new_neighbor = get_face(nbv);

            // Store the new neighbor relation. This cannot be applied yet because
            // it would disturb the functioning of get_face( ... )
            new_neighbor_relations.push_back(make_triple(it, i, new_neighbor));
          }
      }
    // Apply the new neighbor relations now.
    for (unsigned int i = 0; i < new_neighbor_relations.size(); i++)
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
    for (Face_iterator it = all_faces_begin(); it != all_faces_end(); ++it)
      {
        // Skip all faces that are marked to delete
        if (it->get_additional_flag() == 1)
          continue;
        // Find the corresponding vertices of it in the original domain
        // and set them as new vertices of it.
        for (int i = 0; i < 3; i++)
          {
            off[i] = Offset();
            get_vertex(it, i, vert[i], off[i]);
            it->set_vertex(i, vert[i]);
            CGAL_triangulation_assertion(vert[i]->point()[0] < _domain.xmax());
            CGAL_triangulation_assertion(vert[i]->point()[1] < _domain.ymax());
            CGAL_triangulation_assertion(vert[i]->point()[0] >= _domain.xmin());
            CGAL_triangulation_assertion(vert[i]->point()[1] >= _domain.ymin());

            // redirect also the face pointer of the vertex.
            it->vertex(i)->set_face(it);
          }
        // Set the offsets.
        set_offsets(it, off[0], off[1], off[2]);
        CGAL_triangulation_assertion( int_to_off(it->offset(0)) == off[0] );
        CGAL_triangulation_assertion( int_to_off(it->offset(1)) == off[1] );
        CGAL_triangulation_assertion( int_to_off(it->offset(2)) == off[2] );
      }
  }

  // ###################################################################
  // ### Fourth face iteration #########################################
  // ###################################################################
  {
    // Delete the marked faces.
    std::vector<Face_handle> faces_to_delete;
    for (Face_iterator fit = all_faces_begin(); fit != all_faces_end(); ++fit)
      {
        if (fit->get_additional_flag() == 1)
          faces_to_delete.push_back(fit);
      }
    for (typename std::vector<Face_handle>::iterator it =
           faces_to_delete.begin(); it != faces_to_delete.end(); ++it)
      {
        _tds.delete_face(*it);
      }
  }

  // ###################################################################
  // ### Vertex iteration ##############################################
  // ###################################################################
  {
    // Delete all the vertices in _virtual_vertices, that is all vertices
    // outside the original domain.
    std::vector<Vertex_handle> vertices_to_delete;
    for (Vertex_iterator vit = all_vertices_begin(); vit != all_vertices_end(); ++vit)
      {
        if (_virtual_vertices.count(vit) != 0)
          {
            CGAL_triangulation_assertion( _virtual_vertices.count( vit ) == 1 );
            vertices_to_delete.push_back(vit);
          }
      }
    for (typename std::vector<Vertex_handle>::iterator it =
           vertices_to_delete.begin(); it != vertices_to_delete.end(); ++it)
      {
        _tds.delete_vertex(*it);
      }
  }

  _cover = make_array(1, 1);
  _virtual_vertices.clear();
  _virtual_vertices_reverse.clear();
  _too_long_edge_counter = 0;
  _too_long_edges.clear();

  CGAL_triangulation_assertion(is_1_cover());
}

template<class Gt, class Tds>
void Periodic_2_triangulation_2<Gt, Tds>::convert_to_9_sheeted_covering()
{
  if (_cover == make_array(3, 3))
    return;
  CGAL_triangulation_precondition(is_1_cover());

  // Create 9 copies of each vertex and write virtual_vertices and
  // virtual_vertices_reverse
  std::list<Vertex_handle> original_vertices;
  // try to use std::copy instead of the following loop.
  for (Vertex_iterator vit = vertices_begin(); vit != vertices_end(); ++vit)
    original_vertices.push_back(vit);
  for (typename std::list<Vertex_handle>::iterator vit =
         original_vertices.begin(); vit != original_vertices.end(); ++vit)
    {
      Vertex_handle v_cp;
      std::vector<Vertex_handle> copies;
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
          {
            if (i == 0 && j == 0)
              continue;
            v_cp = _tds.create_vertex(*vit);
            copies.push_back(v_cp);
            _virtual_vertices.insert(std::make_pair(v_cp, std::make_pair(*vit,
                                                    Offset(i, j))));
          }
      _virtual_vertices_reverse.insert(std::make_pair(*vit, copies));
    }

  // Create 9 copies of each face from the respective copies of the
  // vertices and write virtual_faces and virtual_faces_reverse.
  typedef std::map<Face_handle, std::pair<Face_handle, Offset> >
  Virtual_face_map;
  typedef std::map<Face_handle, std::vector<Face_handle> >
  Virtual_face_reverse_map;
  typedef typename Virtual_face_reverse_map::const_iterator VCRMIT;

  Virtual_face_map virtual_faces;
  Virtual_face_reverse_map virtual_faces_reverse;

  std::list<Face_handle> original_faces;
  for (Face_iterator fit = faces_begin(); fit != faces_end(); ++fit)
    original_faces.push_back(fit);

  // Store vertex offsets in a separate data structure
  std::list<Offset> off_v;
  for (typename std::list<Vertex_handle>::iterator vit =
         original_vertices.begin(); vit != original_vertices.end(); ++vit)
    {
      Face_handle ccc = (*vit)->face();
      int v_index = ccc->index(*vit);
      off_v.push_back(int_to_off(ccc->offset(v_index)));
    }

  // Store neighboring offsets in a separate data structure
  std::list<array<Offset, 3> > off_nb;
  for (typename std::list<Face_handle>::iterator fit = original_faces.begin(); fit
       != original_faces.end(); ++fit)
    {
      array<Offset, 3> off_nb_f;
      for (int i = 0; i < 3; i++)
        {
          Face_handle fff = *fit;
          Face_handle nnn = fff->neighbor(i);
          off_nb_f[i] = get_neighbor_offset(fff, i, nnn, nnn->index(fff));
        }
      off_nb.push_back(off_nb_f);
    }

  // Create copies of faces
  for (typename std::list<Face_handle>::iterator fit = original_faces.begin(); fit
       != original_faces.end(); ++fit)
    {
      Face_handle c_cp;
      Vertex_handle v0, v1, v2;
      std::vector<Face_handle> copies;
      Virtual_vertex_reverse_map_it vvrmit[3];
      Offset vvoff[3];
      for (int i = 0; i < 3; i++)
        {
          vvrmit[i] = _virtual_vertices_reverse.find((*fit)->vertex(i));
          CGAL_triangulation_assertion(
            vvrmit[i] != _virtual_vertices_reverse.end());
          vvoff[i] = int_to_off((*fit)->offset(i));
        }
      Vertex_handle vvh[3];
      for (int n = 0; n < 8; n++)   // iterate over faces
        {
          for (int i = 0; i < 3; i++)   // iterate over vertices of the face
            {
              // Decomposition of n into an offset (nx,ny):
              // nx = ((n+1)/3)%3, ny = (n+1)%3
              int o_i = ((n + 1) / 3 + vvoff[i].x() + 3) % 3;
              int o_j = ((n + 1) + vvoff[i].y() + 3) % 3;
              int n_c = 3 * o_i + o_j - 1;
              CGAL_triangulation_assertion(n_c >= -1);
              if (n_c == -1)
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
  for (typename std::list<Face_handle>::iterator fit = original_faces.begin(); fit
       != original_faces.end(); ++fit)
    {
      for (int i = 0; i < 3; i++)
        {
          Virtual_vertex_reverse_map_it vvrmit = _virtual_vertices_reverse.find(
              (*fit)->vertex(i));
          CGAL_triangulation_assertion(vvrmit != _virtual_vertices_reverse.end());
          Offset vvoff = int_to_off((*fit)->offset(i));
          if (!vvoff.is_null())
            {
              int n_f = 3 * vvoff.x() + vvoff.y() - 1;
              CGAL_triangulation_assertion(n_f >= 0);
              CGAL_triangulation_assertion(static_cast<unsigned int>(n_f) < vvrmit->second.size());
              (*fit)->set_vertex(i, vvrmit->second[n_f]);
            }
        }
    }

  // Set neighboring relations of face copies
  typename std::list<array<Offset, 3> >::iterator oit = off_nb.begin();
  for (typename std::list<Face_handle>::iterator fit = original_faces.begin(); fit
       != original_faces.end(); ++fit, ++oit)
    {
      CGAL_triangulation_assertion( oit != off_nb.end() );
      VCRMIT c_cp = virtual_faces_reverse.find(*fit);
      CGAL_triangulation_assertion(c_cp != virtual_faces_reverse.end());
      for (int i = 0; i < 3; i++)
        {
          Face_handle fit_nb = (*fit)->neighbor(i);
          VCRMIT c_cp_nb = virtual_faces_reverse.find(fit_nb);
          CGAL_triangulation_assertion(c_cp_nb != virtual_faces_reverse.end());
          Offset nboff = (*oit)[i];
          for (int n = 0; n < 8; n++)
            {
              int n_nb;
              if (nboff.is_null())
                n_nb = n;
              else
                {
                  int o_i = ((n + 1) / 3 - nboff.x() + 3) % 3;
                  int o_j = (n + 1 - nboff.y() + 3) % 3;
                  n_nb = 3 * o_i + o_j - 1;
                }
              if (n_nb == -1)
                {
                  CGAL_triangulation_assertion(fit_nb->has_vertex(c_cp->second[n]->vertex(ccw(i))) );
                  CGAL_triangulation_assertion(fit_nb->has_vertex(c_cp->second[n]->vertex( cw(i))) );
                  c_cp->second[n]->set_neighbor(i, fit_nb);
                }
              else
                {
                  CGAL_triangulation_assertion(n_nb >= 0);
                  CGAL_triangulation_assertion(static_cast<unsigned int>(n_nb) <= c_cp_nb->second.size());
                  CGAL_triangulation_assertion(c_cp_nb->second[n_nb]->has_vertex(c_cp->second[n]->vertex(ccw(i))) );
                  CGAL_triangulation_assertion(c_cp_nb->second[n_nb]->has_vertex(c_cp->second[n]->vertex( cw(i))) );
                  c_cp->second[n]->set_neighbor(i, c_cp_nb->second[n_nb]);
                }
            }
        }
    }

  // Set neighboring relations of original faces
  oit = off_nb.begin();
  for (typename std::list<Face_handle>::iterator fit = original_faces.begin(); fit
       != original_faces.end(); ++fit, ++oit)
    {
      CGAL_triangulation_assertion( oit != off_nb.end() );
      for (int i = 0; i < 3; i++)
        {
          Offset nboff = (*oit)[i];
          if (!nboff.is_null())
            {
              Face_handle fit_nb = (*fit)->neighbor(i);
              VCRMIT c_cp_nb = virtual_faces_reverse.find(fit_nb);
              CGAL_triangulation_assertion(c_cp_nb != virtual_faces_reverse.end());
              int o_i = (3 - nboff.x()) % 3;
              int o_j = (3 - nboff.y()) % 3;
              int n_nb = 3 * o_i + o_j - 1;
              CGAL_triangulation_assertion(n_nb >= 0);
              CGAL_triangulation_assertion(static_cast<unsigned int>(n_nb) <= c_cp_nb->second.size());
              CGAL_triangulation_assertion(c_cp_nb->second[n_nb]->has_vertex((*fit)->vertex(ccw(i))) );
              CGAL_triangulation_assertion(c_cp_nb->second[n_nb]->has_vertex((*fit)->vertex( cw(i))) );
              (*fit)->set_neighbor(i, c_cp_nb->second[n_nb]);
            }
        }
    }

  // Set incident faces
  for (Face_iterator fit = faces_begin(); fit != faces_end(); ++fit)
    {
      for (int i = 0; i < 3; i++)
        {
          fit->vertex(i)->set_face(fit);
        }
    }

  // Set offsets where necessary
  for (typename std::list<Face_handle>::iterator fit = original_faces.begin(); fit
       != original_faces.end(); ++fit)
    {
      VCRMIT c_cp = virtual_faces_reverse.find(*fit);
      CGAL_triangulation_assertion( c_cp != virtual_faces_reverse.end());
      Offset off[3];
      for (int i = 0; i < 3; i++)
        off[i] = int_to_off((*fit)->offset(i));
      if (off[0].is_null() && off[1].is_null() && off[2].is_null())
        continue;
      for (int n = 0; n < 8; n++)
        {
          Offset off_cp[4];
          int o_i = ((n + 1) / 3) % 3;
          int o_j = (n + 1) % 3;
          if (o_i != 2 && o_j != 2)
            continue;
          for (int i = 0; i < 3; i++)
            {
              off_cp[i] = Offset((o_i == 2) ? off[i].x() : 0, (o_j == 2) ? off[i].y()
                                 : 0);
              CGAL_triangulation_assertion(off_cp[i].x() == 0 || off_cp[i].x() == 1);
              CGAL_triangulation_assertion(off_cp[i].y() == 0 || off_cp[i].y() == 1);
            }
          set_offsets(c_cp->second[n], off_cp[0], off_cp[1], off_cp[2]);
        }
    }

  // Iterate over all original faces and reset offsets.
  for (typename std::list<Face_handle>::iterator fit = original_faces.begin(); fit
       != original_faces.end(); ++fit)
    {
      //This statement does not seem to have any effect
      set_offsets(*fit, 0, 0, 0);
      CGAL_triangulation_assertion((*fit)->offset(0) == 0);
      CGAL_triangulation_assertion((*fit)->offset(1) == 0);
      CGAL_triangulation_assertion((*fit)->offset(2) == 0);
    }

  _cover = make_array(3, 3);

  // Set up too long edges data structure
  int i = 0;
  for (Vertex_iterator vit = vertices_begin(); vit != vertices_end(); ++vit)
    {
      _too_long_edges[vit] = std::list<Vertex_handle>();
      ++i;
    }
  _too_long_edge_counter = find_too_long_edges(_too_long_edges);

  CGAL_triangulation_expensive_assertion(is_valid());

  CGAL_triangulation_assertion(!is_1_cover());
}

// iterate over all edges and store the ones that are longer than
// edge_length_threshold in edges. Return the number of too long edges.
template<class GT, class Tds>
inline int Periodic_2_triangulation_2<GT, Tds>::find_too_long_edges(std::map <
    Vertex_handle, std::list<Vertex_handle> > & edges) const
{
  Point p1, p2;
  int counter = 0;
  Vertex_handle v_no, vh;
  for (Edge_iterator eit = edges_begin(); eit != edges_end(); eit++)
    {
      p1 = construct_point(eit->first->vertex(ccw(eit->second))->point(),
                           get_offset(eit->first, ccw(eit->second)));
      p2 = construct_point(eit->first->vertex(cw(eit->second))->point(),
                           get_offset(eit->first, cw(eit->second)));
      if (edge_is_too_long(p1, p2))
        {
          if (&*(eit->first->vertex(ccw(eit->second))) < &*(eit->first->vertex(cw(
                eit->second))))
            {
              v_no = eit->first->vertex(ccw(eit->second));
              vh = eit->first->vertex(cw(eit->second));
            }
          else
            {
              v_no = eit->first->vertex(cw(eit->second));
              vh = eit->first->vertex(ccw(eit->second));
            }
          edges[v_no].push_back(vh);
          counter++;
        }
    }
  return counter;
}

/**
 * - fh->offset(i) is an bit tuple encapsulated in an integer. Each bit
 *   represents the offset in one direction --> 2-cover!
 * - int_to_off(int) decodes this again.
 * - Finally the offset vector is multiplied by cover.
 *   So if we are working in 3-cover we translate it to the neighboring
 *   3-cover and not only to the neighboring domain.
 */
template<class GT, class Tds>
inline void Periodic_2_triangulation_2<GT, Tds>::get_vertex(Face_handle fh,
    int i, Vertex_handle &vh, Offset &off) const
{
  off = combine_offsets(Offset(), int_to_off(fh->offset(i)));
  vh = fh->vertex(i);

  if (is_1_cover())
    return;
  Vertex_handle vh_i = vh;
  get_vertex(vh_i, vh, off);
  return;
}

template<class GT, class Tds>
inline void Periodic_2_triangulation_2<GT, Tds>::get_vertex(Vertex_handle vh_i,
    Vertex_handle &vh, Offset &off) const
{
  Virtual_vertex_map_it it = _virtual_vertices.find(vh_i);

  if (it == _virtual_vertices.end())
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
      CGAL_triangulation_assertion(vh->point().x() < _domain.xmax());
      CGAL_triangulation_assertion(vh->point().y() < _domain.ymax());
      CGAL_triangulation_assertion(vh->point().x() >= _domain.xmin());
      CGAL_triangulation_assertion(vh->point().y() >= _domain.ymin());
    }
}

/** Find the Face that consists of the three given vertices
 *
 *  Iterates over all faces and compare the three vertices of each face
 *  with the three vertices in vh.
 */
template<class GT, class Tds>
inline typename Periodic_2_triangulation_2<GT, Tds>::Face_handle Periodic_2_triangulation_2 <
GT, Tds >::get_face(const Vertex_handle* vh) const
{
  bool contains_v[2];
  Face_circulator fc = incident_faces(vh[2]);
  Face_circulator done(fc);
  do
    {
      CGAL_triangulation_assertion(
        fc->vertex(0) == vh[2] ||
        fc->vertex(1) == vh[2] ||
        fc->vertex(2) == vh[2]);
      for (int j = 0; j < 2; j++)
        {
          contains_v[j] = (fc->vertex(0) == vh[j]) || (fc->vertex(1) == vh[j])
                          || (fc->vertex(2) == vh[j]);
        }
      if (contains_v[0] && contains_v[1])
        {
          return fc;
        }
    }
  while (++fc != done);

  CGAL_triangulation_assertion(false);
  return Face_handle();
}

template<class Gt, class Tds>
Bounded_side Periodic_2_triangulation_2<Gt, Tds>::side_of_face(const Point &q,
    const Offset &off, Face_handle f, Locate_type &lt, int &li) const
{
  CGAL_triangulation_precondition(number_of_vertices() != 0);

  Orientation o0, o1, o2;
  o0 = o1 = o2 = ZERO;
  int cumm_off = f->offset(0) | f->offset(1) | f->offset(2);

  if ((cumm_off == 0) && is_1_cover())
    {
      CGAL_triangulation_assertion(off == Offset());

      const Point &p0 = f->vertex(0)->point();
      const Point &p1 = f->vertex(1)->point();
      const Point &p2 = f->vertex(2)->point();

      if (((o0 = orientation(q, p1, p2)) == NEGATIVE) || ((o1 = orientation(p0,
          q, p2)) == NEGATIVE) || ((o2 = orientation(p0, p1, q)) == NEGATIVE))
        {
          return ON_UNBOUNDED_SIDE;
        }
    }
  else     // Special case for the periodic space.
    {
      Offset off_q;
      Offset offs[3];
      const Point *p[3];
      for (int i = 0; i < 3; i++)
        {
          p[i] = &(f->vertex(i)->point());
          offs[i] = get_offset(f, i);
        }
      CGAL_triangulation_assertion(orientation(*p[0], *p[1], *p[2],
                                   offs[0], offs[1], offs[2]) == POSITIVE);
      bool found = false;
      for (int i = 0; (i < 4) && (!found); i++)
        {
          if ((cumm_off | ((~i) & 3)) == 3)
            {
              o0 = o1 = o2 = NEGATIVE;
              off_q = combine_offsets(off, int_to_off(i));

              if (((o0 = orientation(q, *p[1], *p[2], off_q, offs[1], offs[2]))
                   != NEGATIVE) && ((o1 = orientation(*p[0], q, *p[2], offs[0], off_q,
                                                      offs[2])) != NEGATIVE) && ((o2 = orientation(*p[0], *p[1], q,
                                                          offs[0], offs[1], off_q)) != NEGATIVE))
                {
                  found = true;
                }
            }
        }
      if (!found)
        return ON_UNBOUNDED_SIDE;
    }

  // now all the oi's are >=0
  // sum gives the number of facets p lies on
  int sum = ((o0 == ZERO) ? 1 : 0) + ((o1 == ZERO) ? 1 : 0) + ((o2 == ZERO) ? 1
            : 0);

  switch (sum)
    {
    case 0:
    {
      lt = FACE;
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
Oriented_side Periodic_2_triangulation_2<Gt, Tds>::oriented_side(Face_handle f,
    const Point& p, const Offset &o) const
{
  Point &p0 = f->vertex(0)->point();
  Point &p1 = f->vertex(1)->point();
  Point &p2 = f->vertex(2)->point();

  int cumm_off = f->offset(0) | f->offset(1) | f->offset(2);

  if ((cumm_off == 0) && is_1_cover())
    {
      CGAL_precondition(o == Offset());

      // return position of point p with respect to the oriented triangle p0p1p2
      // the orientation of the vertices is assumed to be counter clockwise
      CGAL_assertion(orientation(p0, p1, p2) == LEFT_TURN);

      Bounded_side bs = bounded_side(p0, p1, p2, p);
      switch (bs)
        {
        case ON_BOUNDARY:
          return ON_ORIENTED_BOUNDARY;
        case ON_BOUNDED_SIDE:
          return ON_POSITIVE_SIDE;
        case ON_UNBOUNDED_SIDE:
          return ON_NEGATIVE_SIDE;
        }
    }
  else     // Special case for the periodic space.
    {
      Offset off_q;
      Offset off0 = get_offset(f, 0);
      Offset off1 = get_offset(f, 1);
      Offset off2 = get_offset(f, 2);

      // return position of point p with respect to the oriented triangle p0p1p2
      // the orientation of the vertices is assumed to be counter clockwise
      CGAL_assertion(orientation(p0, p1, p2, off0, off1, off2) == LEFT_TURN);

      Bounded_side bs;
      for (int i = 0; (i <= 7); i++)
        {
          if ((cumm_off | ((~i) & 3)) == 3)
            {
              off_q = combine_offsets(o, int_to_off(i));
              bs = bounded_side(p0, p1, p2, p, off0, off1, off2, off_q);
              if (bs != ON_UNBOUNDED_SIDE)
                {
                  return (bs == ON_BOUNDARY ? ON_ORIENTED_BOUNDARY : ON_POSITIVE_SIDE);
                }
            }
        }

      return ON_NEGATIVE_SIDE;
    }

  CGAL_assertion(false);
  return ON_NEGATIVE_SIDE;
}

template<class Gt, class Tds>
Bounded_side Periodic_2_triangulation_2<Gt, Tds>::bounded_side(const Point &p0, const Point &p1, const Point &p2, const Point &p) const
{

  // return position of point p with respect to triangle p0p1p2
  CGAL_triangulation_precondition( orientation(p0, p1, p2) != COLLINEAR);

  Orientation o1 = orientation(p0, p1, p);
  Orientation o2 = orientation(p1, p2, p);
  Orientation o3 = orientation(p2, p0, p);

  if (o1 == COLLINEAR)
    {
      if (o2 == COLLINEAR || o3 == COLLINEAR)
        return ON_BOUNDARY;
      if (collinear_between(p0, p, p1))
        return ON_BOUNDARY;
      return ON_UNBOUNDED_SIDE;
    }

  if (o2 == COLLINEAR)
    {
      if (o3 == COLLINEAR)
        return ON_BOUNDARY;
      if (collinear_between(p1, p, p2))
        return ON_BOUNDARY;
      return ON_UNBOUNDED_SIDE;
    }

  if (o3 == COLLINEAR)
    {
      if (collinear_between(p2, p, p0))
        return ON_BOUNDARY;
      return ON_UNBOUNDED_SIDE;
    }

  // from here none ot, o1, o2 and o3 are known to be non null
  if (o1 == o2 && o2 == o3)
    return ON_BOUNDED_SIDE;
  return ON_UNBOUNDED_SIDE;
}

template<class Gt, class Tds>
Bounded_side Periodic_2_triangulation_2<Gt, Tds>::bounded_side(const Point &p0,
    const Point &p1, const Point &p2, const Point &p, const Offset &o0,
    const Offset &o1, const Offset &o2, const Offset &o) const
{
  // return position of point p with respect to triangle p0p1p2
  CGAL_triangulation_precondition( orientation(p0, p1, p2, o0, o1, o2) != COLLINEAR);
  Orientation orient1 = orientation(p0, p1, p, o0, o1, o);
  Orientation orient2 = orientation(p1, p2, p, o1, o2, o);
  Orientation orient3 = orientation(p2, p0, p, o2, o0, o);

  if (orient1 == COLLINEAR)
    {
      if (orient2 == COLLINEAR || orient3 == COLLINEAR)
        return ON_BOUNDARY;
      if (collinear_between(p0, p, p1, o0, o, o1))
        return ON_BOUNDARY;
      return ON_UNBOUNDED_SIDE;
    }

  if (orient2 == COLLINEAR)
    {
      if (orient3 == COLLINEAR)
        return ON_BOUNDARY;
      if (collinear_between(p1, p, p2, o1, o, o2))
        return ON_BOUNDARY;
      return ON_UNBOUNDED_SIDE;
    }

  if (orient3 == COLLINEAR)
    {
      if (collinear_between(p2, p, p0, o2, o, o0))
        return ON_BOUNDARY;
      return ON_UNBOUNDED_SIDE;
    }

  // from here none ot, o1, o2 and o3 are known to be non null
  if (orient1 == orient2 && orient2 == orient3)
    return ON_BOUNDED_SIDE;
  return ON_UNBOUNDED_SIDE;
}


template<class Gt, class Tds>
bool Periodic_2_triangulation_2<Gt, Tds>::collinear_between(const Point& p,
    const Point& q, const Point& r) const
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
  return ( (c_pq == SMALLER) && (c_qr == SMALLER) ) ||
         ( (c_pq == LARGER)  && (c_qr == LARGER) );

}

template<class Gt, class Tds>
bool Periodic_2_triangulation_2<Gt, Tds>::collinear_between(const Point& p,
    const Point& q, const Point& r, const Offset& o_p, const Offset& o_q,
    const Offset& o_r) const
{
  // return true if point q is strictly between p and r
  // p,q and r are supposed to be collinear points
  Comparison_result c_pr = compare_x(p, r, o_p, o_r);
  Comparison_result c_pq;
  Comparison_result c_qr;
  if (c_pr == EQUAL)
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
inline Comparison_result Periodic_2_triangulation_2<Gt, Tds>::compare_x(
  const Point& p, const Point& q) const
{
  return geom_traits().compare_x_2_object()(p, q);
}

template<class Gt, class Tds>
inline Comparison_result Periodic_2_triangulation_2<Gt, Tds>::compare_x(
  const Point& p, const Point& q, const Offset &o_p, const Offset &o_q) const
{
  return geom_traits().compare_x_2_object()(p, q, o_p, o_q);
}

template<class Gt, class Tds>
inline Comparison_result Periodic_2_triangulation_2<Gt, Tds>::compare_xy(
  const Point& p, const Point& q) const
{
  Comparison_result res = geom_traits().compare_x_2_object()(p, q);
  if (res == EQUAL)
    {
      return geom_traits().compare_y_2_object()(p, q);
    }
  return res;
}

template<class Gt, class Tds>
inline Comparison_result Periodic_2_triangulation_2<Gt, Tds>::compare_xy(
  const Point& p, const Point& q, const Offset &o_p, const Offset &o_q) const
{
  Comparison_result res = geom_traits().compare_x_2_object()(p, q, o_p, o_q);
  if (res == EQUAL)
    {
      return geom_traits().compare_y_2_object()(p, q, o_p, o_q);
    }
  return res;
}

template<class Gt, class Tds>
inline Comparison_result Periodic_2_triangulation_2<Gt, Tds>::compare_y(
  const Point& p, const Point& q) const
{
  return geom_traits().compare_y_2_object()(p, q);
}

template<class Gt, class Tds>
inline Comparison_result Periodic_2_triangulation_2<Gt, Tds>::compare_y(
  const Point& p, const Point& q, const Offset &o_p, const Offset &o_q) const
{
  return geom_traits().compare_y_2_object()(p, q, o_p, o_q);
}

template<class Gt, class Tds>
inline
bool Periodic_2_triangulation_2<Gt, Tds>::xy_equal(const Point& p,
    const Point& q) const
{
  return compare_xy(p, q) == EQUAL;
}

template<class Gt, class Tds>
inline Orientation Periodic_2_triangulation_2<Gt, Tds>::orientation(
  const Point& p0, const Point& p1, const Point& p2) const
{
  return geom_traits().orientation_2_object()(p0, p1, p2);
}
template<class Gt, class Tds>
inline Orientation Periodic_2_triangulation_2<Gt, Tds>::orientation(
  const Point& p0, const Point& p1, const Point& p2, const Offset& o0,
  const Offset& o1, const Offset& o2) const
{
  return geom_traits().orientation_2_object()(p0, p1, p2, o0, o1, o2);
}


template<class Gt, class Tds>
Oriented_side Periodic_2_triangulation_2<Gt, Tds>::side_of_oriented_circle(
  const Point &p0, const Point &p1, const Point &p2, const Point &p,
  bool perturb) const
{
  Oriented_side os = geom_traits().side_of_oriented_circle_2_object()(p0, p1, p2, p);
  if ((os != ON_ORIENTED_BOUNDARY) || (!perturb))
    return os;

  // We are now in a degenerate case => we do a symbolic perturbation.

  // We sort the points lexicographically.
  const Point * points[4] = { &p0, &p1, &p2, &p };
  std::sort(points, points + 4, Perturbation_order(this));

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
Oriented_side Periodic_2_triangulation_2<Gt, Tds>::side_of_oriented_circle(
  const Point &p0, const Point &p1, const Point &p2, const Point &p,
  const Offset &o0, const Offset &o1, const Offset &o2, const Offset &o,
  bool perturb) const
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

  std::sort(points, points + 4, Perturbation_order(this));

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
Oriented_side Periodic_2_triangulation_2<Gt, Tds>::side_of_oriented_circle(
  Face_handle f, const Point & p, bool perturb) const
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


template<class Gt, class Tds>
void Periodic_2_triangulation_2<Gt, Tds>::insert_too_long_edges_in_star(Vertex_handle vh)
{
  // Insert the too long edges in the star of vh
  Face_handle f = vh->face();
  Face_handle f_start = f;

  do
    {
      int i = ccw(f->index(vh));

      insert_too_long_edge(f, i);

      // Proceed to the next face
      f = f->neighbor(i);
    }
  while (f != f_start);
}

template<class Gt, class Tds>
void Periodic_2_triangulation_2<Gt, Tds>::insert_too_long_edge(Face_handle f, int i)
{
  Vertex_handle vh1 = f->vertex(ccw(i));
  Vertex_handle vh2 = f->vertex(cw(i));
  CGAL_assertion(vh1 != Vertex_handle());
  CGAL_assertion(vh2 != Vertex_handle());
  Point p1 = construct_point(vh1->point(), get_offset(f, ccw(i)));
  Point p2 = construct_point(vh2->point(), get_offset(f, cw(i)));

  if (&*vh1 < &*vh2)
    {
      if (edge_is_too_long(p1, p2) &&
          (find(_too_long_edges[vh1].begin(), _too_long_edges[vh1].end(), vh2) == _too_long_edges[vh1].end()))
        {
          _too_long_edges[vh1].push_back(vh2);
          _too_long_edge_counter++;
        }
    }
  else
    {
      CGAL_triangulation_precondition(&*vh2 < &*vh1);
      if (edge_is_too_long(p2, p1) &&
          (find(_too_long_edges[vh2].begin(), _too_long_edges[vh2].end(), vh1) == _too_long_edges[vh2].end()))
        {
          _too_long_edges[vh2].push_back(vh1);
          _too_long_edge_counter++;
        }
    }
}

template<class Gt, class Tds>
void Periodic_2_triangulation_2<Gt, Tds>::remove_too_long_edges_in_star(
  Vertex_handle vh)
{
  if (is_1_cover())
    return;

  // Insert the too long edges in the star of vh
  Face_handle f = vh->face();
  Face_handle f_start = f;

  do
    {
      int i = f->index(vh);
      int i2 = ccw(i);
      Vertex_handle vh2 = f->vertex(i2);

      // Point corresponding to v
      Point p1 = construct_point(vh->point(), get_offset(f, f->index(vh)));
      // Point corresponding to the other vertex
      Point p2 = construct_point(vh2->point(), get_offset(f, i2));

      if (&*vh < &*vh2)
        {
          if (edge_is_too_long(p1, p2) &&
              (find(_too_long_edges[vh].begin(), _too_long_edges[vh].end(), vh2) !=
               _too_long_edges[vh].end()))
            {
              _too_long_edges[vh].remove(vh2);
              _too_long_edge_counter--;
            }
        }
      else
        {
          CGAL_triangulation_precondition(&*vh2 < &*vh);
          if (edge_is_too_long(p1, p2) &&
              (find(_too_long_edges[vh2].begin(), _too_long_edges[vh2].end(), vh) !=
               _too_long_edges[vh2].end()))
            {
              _too_long_edges[vh2].remove(vh);
              _too_long_edge_counter--;
            }
        }

      // Proceed to the next face
      f = f->neighbor(i2);
    }
  while (f != f_start);
}

template<class Gt, class Tds>
void Periodic_2_triangulation_2<Gt, Tds>::remove_too_long_edge(Face_handle f,
    int i)
{
  Vertex_handle vh1 = f->vertex(cw(i));
  Vertex_handle vh2 = f->vertex(ccw(i));
  Point p1 = construct_point(vh1->point(), get_offset(f, cw(i)));
  Point p2 = construct_point(vh2->point(), get_offset(f, ccw(i)));
  if (edge_is_too_long(p1, p2))
    {
      if (&*vh1 < &*vh2)
        {
          typename std::list<Vertex_handle>::iterator it = find(
                _too_long_edges[vh1].begin(), _too_long_edges[vh1].end(), vh2);
          if (it != _too_long_edges[vh1].end())
            {
              _too_long_edges[vh1].erase(it);
              _too_long_edge_counter--;
            }
        }
      else
        {
          typename std::list<Vertex_handle>::iterator it = find(
                _too_long_edges[vh2].begin(), _too_long_edges[vh2].end(), vh1);
          if (it != _too_long_edges[vh2].end())
            {
              _too_long_edges[vh2].erase(it);
              _too_long_edge_counter--;
            }
        }
    }
}

template<class Gt, class Tds>
bool Periodic_2_triangulation_2<Gt, Tds>::edge_is_too_long(const Point &p1,
    const Point &p2) const
{
  return squared_distance(p1, p2) > _edge_length_threshold;
}

template<class GT, class Tds>
inline bool Periodic_2_triangulation_2<GT, Tds>::is_extensible_triangulation_in_1_sheet_h1() const
{
  if (!is_1_cover())
    {
      if (_too_long_edge_counter == 0)
        return true;
      else
        return false;
    }
  else
    {
      typename Geom_traits::FT longest_edge_squared_length(0);
      Segment s;
      for (Periodic_segment_iterator psit = periodic_segments_begin(UNIQUE); psit
           != periodic_segments_end(UNIQUE); ++psit)
        {
          s = construct_segment(*psit);
          longest_edge_squared_length = (std::max)(longest_edge_squared_length,
                                        s.squared_length());
        }
      return (longest_edge_squared_length < _edge_length_threshold);
    }
}

template<class GT, class Tds>
inline bool Periodic_2_triangulation_2<GT, Tds>::is_extensible_triangulation_in_1_sheet_h2() const
{
  typedef typename Geom_traits::Construct_circumcenter_2 Construct_circumcenter;
  typedef typename Geom_traits::FT FT;
  Construct_circumcenter construct_circumcenter =
    _gt.construct_circumcenter_2_object();
  for (Periodic_triangle_iterator tit = periodic_triangles_begin(UNIQUE); tit
       != periodic_triangles_end(UNIQUE); ++tit)
    {
      Point cc = construct_circumcenter(tit->at(0).first, tit->at(1).first,
                                        tit->at(2).first, tit->at(0).second, tit->at(1).second,
                                        tit->at(2).second);

      if (!(FT(16) * squared_distance(cc, point(tit->at(0))) < (_domain.xmax()
            - _domain.xmin()) * (_domain.xmax() - _domain.xmin())))
        return false;
    }
  return true;
}

template<class GT, class Tds>
inline bool Periodic_2_triangulation_2<GT, Tds>::is_triangulation_in_1_sheet() const
{
  if (is_1_cover())
    return true;
  for (Vertex_iterator vit = vertices_begin(); vit != vertices_end(); ++vit)
    {
      if (_virtual_vertices.find(vit) == _virtual_vertices.end())
        continue;
      std::set<Vertex_handle> nb_v_odom;
      Vertex_handle vh;
      Offset off;
      Vertex_circulator vcir = adjacent_vertices(vit);
      Vertex_circulator vstart = vcir;
      size_t degree = 0;
      do
        {
          get_vertex(vcir, vh, off);
          nb_v_odom.insert(vh);
          degree++;
        }
      while (++vcir != vstart);
      if (degree != nb_v_odom.size())
        return false;
    }
  return true;
}

template<class Gt, class Tds>
std::ostream&
Periodic_2_triangulation_2<Gt, Tds>::save(std::ostream& os) const
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


  if (is_ascii(os))
    os << domain() << std::endl
       << cover[0] << " " << cover[1] << std::endl
       << n*cover[0]*cover[1] << std::endl;
  else
    {
      os << domain();
      write(os, cover[0]);
      write(os, cover[1]);
      write(os, n * cover[0]*cover[1]);
    }
  std::cout << "Line:" << __LINE__ << " cover[0]:" << cover[0] << " cover[1]:" << cover[1] << " n*c0*c1:" << (n * cover[0]*cover[1]) << std::endl;

  std::cout << "save, #Vertices: " << n << std::endl;

  if (n == 0)
    return os;

  // write the vertices
  Unique_hash_map<Vertex_handle, std::size_t > V;
  std::size_t i = 0;
  if (is_1_cover())
    {
      for (Vertex_iterator it = vertices_begin(); it != vertices_end(); ++it)
        {
          V[it] = i++;
          os << it->point();
          if (is_ascii(os))
            os << std::endl;
        }
    }
  else
    {
      Virtual_vertex_map_it vit, vvit;
      std::vector<Vertex_handle> vv;
      for (Vertex_iterator it = vertices_begin(); it != vertices_end(); ++it)
        {
          vit = _virtual_vertices.find(it);
          if (vit != _virtual_vertices.end()) continue;
          V[it] = i++;
          if (is_ascii(os))
            os << it->point() << std::endl
               << Offset(0, 0) << std::endl;
          else
            os << it->point() << Offset(0, 0);
          CGAL_triangulation_assertion(_virtual_vertices_reverse.find(it)
                                       != _virtual_vertices_reverse.end());
          vv = _virtual_vertices_reverse.find(it)->second;
          CGAL_triangulation_assertion(vv.size() == 8);
          for (std::size_t j = 0; j < vv.size(); j++)
            {
              vvit = _virtual_vertices.find(vv[j]);
              CGAL_triangulation_assertion(vvit != _virtual_vertices.end());
              V[vv[j]] = i++;
              if (is_ascii(os))
                os << vv[j]->point() << std::endl
                   << vvit->second.second << std::endl;
              else os << vv[j]->point() << vvit->second.second;
            }
        }
    }
  CGAL_triangulation_postcondition(i == _cover[0]*_cover[1]*n);

  Unique_hash_map<Face_handle, std::size_t> F;
  int inum = 0;
  // asks the tds for the combinatorial information
  // vertices of the faces
  size_type m = _tds.number_of_faces();
  if (is_ascii(os)) os << std::endl << m << std::endl;
  else write(os, m);
  std::cout << "save, #Faces: " << m << std::endl;

  for( Face_iterator ib = faces_begin();
       ib != faces_end(); ++ib)
    {
      F[ib] = inum++;
      for(int j = 0; j < 3 ; ++j)
        {
          if(is_ascii(os)) os << V[ib->vertex(j)] << " ";
          else write(os, V[ib->vertex(j)]);
        }
      os << *ib ;
      if(is_ascii(os)) os << "\n";
    }
  if(is_ascii(os)) os << "\n";

  std::cout << "save, face check: " << inum << " == " << m << std::endl;
  CGAL_assertion(m == (size_type)inum);

  // neighbor pointers of the  faces
  for( Face_iterator it = faces_begin();
       it != faces_end(); ++it)
    {
      for(int j = 0; j < 3; ++j)
        {
          CGAL_assertion(F.is_defined(it->neighbor(j)));
          if(is_ascii(os))  os << F[it->neighbor(j)] << " ";
          else write(os, F[it->neighbor(j)]);
        }
      if(is_ascii(os)) os << "\n";
    }

  // write offsets
  //for (unsigned int i=0 ; i<number_of_faces() ; i++) {
  for (Face_iterator it = faces_begin(); it != faces_end(); ++it)
    {
      //Face_handle ch = std::find(faces_begin(), faces_end(), i);
      Face_handle ch(it);
      for (int j = 0; j < 3; j++)
        {
          if(is_ascii(os))
            {
              os << ch->offset(j);
              if ( j == 3 )
                os << std::endl;
              else
                os << ' ';
            }
          else write(os, ch->offset(j));
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
Periodic_2_triangulation_2<Gt, Tds>::load(std::istream& is)
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

  if (is_ascii(is))
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
  _domain = domain;
  _gt.set_domain(domain);
  _cover = make_array(cx, cy);

  if ( n == 0 ) return is;

  std::map< std::size_t, Vertex_handle > V;

  if (cx == 1 && cy == 1)
    {
      Point p;
      for (std::size_t i = 0; i < n; i++)
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
      for (std::size_t i = 0; i < n; i++)
        {
          v = tds().create_vertex();
          V[i] = v;
          is >> p >> off;
          V[i]->set_point(p);
          vv.clear();
          for (int j = 1; j < cx * cy; j++)
            {
              i++;
              w = tds().create_vertex();
              V[i] = w;
              is >> p >> off;
              V[i]->set_point(p);
              vv.push_back(w);
              _virtual_vertices[w] = std::make_pair(v, off);
            }
          _virtual_vertices_reverse[v] = vv;
        }
    }

  // Creation of the faces
  std::size_t index;
  size_type m;
  if (is_ascii(is)) is >> m;
  else read(is, m);
  std::vector<Face_handle> F(m);
  std::cout << "load, #Faces: " << m << std::endl;
  {
    for(size_t i = 0; i < m; ++i)
      {
        F[i] = _tds.create_face() ;
        for(int j = 0; j < 3 ; ++j)
          {
            if (is_ascii(is)) is >> index;
            else read(is, index);
            CGAL_assertion(index < V.size());
            F[i]->set_vertex(j, V[index]);
            // The face pointer of vertices is set too often,
            // but otherwise we had to use one more map
            V[index]->set_face(F[i]);
          }
        // read in non combinatorial info of the face
        is >> *(F[i]) ;
      }
  }

  // Setting the neighbor pointers
  {
    for(size_t i = 0; i < m; ++i)
      {
        for(int j = 0; j < 3; ++j)
          {
            if (is_ascii(is)) is >> index;
            else read(is, index);
            if (index >= F.size()) {
              std::cout << __FILE__ << ", " << __FUNCTION__ << ", l:" << __LINE__ << "  f="
                        << i << "<" << m << ", index=" << j << " nb=" << index << " #F=" << F.size()
                        << std::endl;
            }
            CGAL_assertion(i < F.size());
            CGAL_assertion(index < F.size());
            F[i]->set_neighbor(j, F[index]);
          }
      }
  }

  // read offsets
  int off[3] = {0, 0, 0};
  for (std::size_t j = 0 ; j < m; j++)
    {
      if (is_ascii(is))
        is >> off[0] >> off[1] >> off[2];
      else
        {
          read(is, off[0]);
          read(is, off[1]);
          read(is, off[2]);
        }
      set_offsets(F[j], off[0], off[1], off[2]);
    }

  // read potential other information
  for (std::size_t j = 0 ; j < m; j++)
    is >> *(F[j]);

  int i = 0;
  for (Vertex_iterator vi = vertices_begin();
       vi != vertices_end(); ++vi)
    {
      _too_long_edges[vi] = std::list<Vertex_handle>();
      ++i;
    }

  _edge_length_threshold = FT(0.166) * (_domain.xmax() - _domain.xmin())
                           * (_domain.xmax() - _domain.xmin());
  _too_long_edge_counter = find_too_long_edges(_too_long_edges);

  CGAL_triangulation_expensive_assertion( is_valid() );
  return is;
}


namespace internal
{

/// Internal function used by operator==.
//TODO: introduce offsets
template <class GT, class Tds1, class Tds2>
bool
test_next(const Periodic_2_triangulation_2<GT, Tds1> &t1,
          const Periodic_2_triangulation_2<GT, Tds2> &t2,
          typename Periodic_2_triangulation_2<GT, Tds1>::Face_handle c1,
          typename Periodic_2_triangulation_2<GT, Tds2>::Face_handle c2,
          std::map < typename Periodic_2_triangulation_2<GT, Tds1>::Face_handle,
          typename Periodic_2_triangulation_2<GT, Tds2>::Face_handle > &Cmap,
          std::map < typename Periodic_2_triangulation_2<GT, Tds1>::Vertex_handle,
          typename Periodic_2_triangulation_2<GT, Tds2>::Vertex_handle > &Vmap)
{
  // This function tests and registers the 4 neighbors of c1/c2,
  // and recursively calls itself over them.
  // Returns false if an inequality has been found.

  // Precondition: c1, c2 have been registered as well as their 4 vertices.
  CGAL_triangulation_precondition(t1.number_of_vertices() != 0);
  CGAL_triangulation_precondition(Cmap[c1] == c2);
  CGAL_triangulation_precondition(Vmap.find(c1->vertex(0)) != Vmap.end());
  CGAL_triangulation_precondition(Vmap.find(c1->vertex(1)) != Vmap.end());
  CGAL_triangulation_precondition(Vmap.find(c1->vertex(2)) != Vmap.end());

  typedef Periodic_2_triangulation_2<GT, Tds1> Tr1;
  typedef Periodic_2_triangulation_2<GT, Tds2> Tr2;
  typedef typename Tr1::Vertex_handle  Vertex_handle1;
  typedef typename Tr1::Face_handle    Face_handle1;
  typedef typename Tr2::Vertex_handle  Vertex_handle2;
  typedef typename Tr2::Face_handle    Face_handle2;
  typedef typename std::map<Face_handle1, Face_handle2>::const_iterator  Cit;
  typedef typename std::map < Vertex_handle1,
          Vertex_handle2 >::const_iterator Vit;

  for (int i = 0; i <= 2; ++i)
    {
      Face_handle1 n1 = c1->neighbor(i);
      Cit cit = Cmap.find(n1);
      Vertex_handle1 v1 = c1->vertex(i);
      Vertex_handle2 v2 = Vmap[v1];
      Face_handle2 n2 = c2->neighbor(c2->index(v2));
      if (cit != Cmap.end())
        {
          // n1 was already registered.
          if (cit->second != n2)
            return false;
          continue;
        }
      // n1 has not yet been registered.
      // We check that the new vertices match geometrically.
      // And we register them.
      Vertex_handle1 vn1 = n1->vertex(n1->index(c1));
      Vertex_handle2 vn2 = n2->vertex(n2->index(c2));
      Vit vit = Vmap.find(vn1);
      if (vit != Vmap.end())
        {
          // vn1 already registered
          if (vit->second != vn2)
            return false;
        }
      else
        {
          if (t1.geom_traits().compare_xy_2_object()(vn1->point(),
              vn2->point()) != 0)
            return false;

          // We register vn1/vn2.
          Vmap.insert(std::make_pair(vn1, vn2));
        }

      // We register n1/n2.
      Cmap.insert(std::make_pair(n1, n2));
      // We recurse on n1/n2.
      if (!test_next(t1, t2, n1, n2, Cmap, Vmap))
        return false;
    }

  return true;
}

} // namespace internal


template<class Gt, class Tds>
std::istream&
operator>>(std::istream& is, Periodic_2_triangulation_2<Gt, Tds> &tr)
{
  return tr.load(is);
}
template<class Gt, class Tds>
std::ostream&
operator<<(std::ostream& os, Periodic_2_triangulation_2<Gt, Tds> &tr)
{
  return tr.save(os);
}

template < class GT, class Tds1, class Tds2  >
bool
operator==(const Periodic_2_triangulation_2<GT, Tds1> &t1,
           const Periodic_2_triangulation_2<GT, Tds2> &t2)
{
  typedef typename Periodic_2_triangulation_2<GT, Tds1>::Vertex_handle
  Vertex_handle1;
  typedef typename Periodic_2_triangulation_2<GT, Tds1>::Face_handle
  Face_handle1;
  typedef typename Periodic_2_triangulation_2<GT, Tds2>::Vertex_handle
  Vertex_handle2;
  typedef typename Periodic_2_triangulation_2<GT, Tds2>::Vertex_handle
  Vertex_iterator2;
  typedef typename Periodic_2_triangulation_2<GT, Tds2>::Face_handle
  Face_handle2;
  typedef typename Periodic_2_triangulation_2<GT, Tds2>::Face_circulator
  Face_circulator2;

  typedef typename Periodic_2_triangulation_2<GT, Tds1>::Point      Point;
  typedef typename Periodic_2_triangulation_2<GT, Tds1>::Offset     Offset;

  // Some quick checks.
  if (   t1.domain()           != t2.domain()
         || t1.number_of_sheets() != t2.number_of_sheets())
    return false;

  if (   t1.number_of_vertices() != t2.number_of_vertices()
         || t1.number_of_faces() != t2.number_of_faces())
    return false;

  // Special case for empty triangulations
  if (t1.number_of_vertices() == 0)
    return true;

  // We will store the mapping between the 2 triangulations vertices and
  // faces in 2 maps.
  std::map<Vertex_handle1, Vertex_handle2> Vmap;
  std::map<Face_handle1, Face_handle2> Cmap;

  // find a common point
  Vertex_handle1 v1 = static_cast<Vertex_handle1>(t1.vertices_begin());
  Vertex_handle2 iv2;
  for (Vertex_iterator2 vit2 = t2.vertices_begin() ;
       vit2 != t2.vertices_end(); ++vit2)
    {
      if (t1.compare_xy(vit2->point(), v1->point(),
                        t2.get_offset(vit2), t1.get_offset(v1)) != EQUAL)
        continue;
      iv2 = static_cast<Vertex_handle2>(vit2);
      break;
    }
  if (iv2 == Vertex_handle2())
    return false;
  Vmap.insert(std::make_pair(v1, iv2));

  // We pick one face of t1, and try to match it against the
  // faces of t2.
  Face_handle1 c = v1->face();
  Vertex_handle1 v2 = c->vertex(t1.cw(c->index(v1)));
  Vertex_handle1 v3 = c->vertex(t1.ccw(c->index(v1)));
  Point p2 = v2->point();
  Point p3 = v3->point();
  Offset o2 = t1.get_offset(v2);
  Offset o3 = t1.get_offset(v3);

  Face_circulator2 fc = t2.incident_faces(iv2), done(fc);
  do
    {
      int inf = fc->index(iv2);

      if (t1.compare_xy(p2, fc->vertex((inf + 1) % 3)->point(),
                        o2, t2.get_offset(fc->vertex((inf + 1) % 3))) == EQUAL)
        Vmap.insert(std::make_pair(v2, fc->vertex((inf + 1) % 3)));
      else if (t1.compare_xy(p2, fc->vertex((inf + 2) % 3)->point(),
                             o2, t2.get_offset(fc->vertex((inf + 2) % 3))) == EQUAL)
        Vmap.insert(std::make_pair(v2, fc->vertex((inf + 2) % 3)));
      else
        continue; // None matched v2.

      if (t1.compare_xy(p3, fc->vertex((inf + 1) % 3)->point(),
                        o3, t2.get_offset(fc->vertex((inf + 1) % 3))) == EQUAL)
        Vmap.insert(std::make_pair(v3, fc->vertex((inf + 1) % 3)));
      else if (t1.compare_xy(p3, fc->vertex((inf + 2) % 3)->point(),
                             o3, t2.get_offset(fc->vertex((inf + 2) % 3))) == EQUAL)
        Vmap.insert(std::make_pair(v3, fc->vertex((inf + 2) % 3)));
      else
        continue; // None matched v3.

      // Found it !
      Cmap.insert(std::make_pair(c, fc));
      break;
    }
  while (++fc != done);

  if (Cmap.size() == 0)
    return false;

  // We now have one face, we need to propagate recursively.
  return internal::test_next(t1, t2,
                             Cmap.begin()->first, Cmap.begin()->second, Cmap, Vmap);
}

template < class GT, class Tds1, class Tds2 >
inline
bool
operator!=(const Periodic_2_triangulation_2<GT, Tds1> &t1,
           const Periodic_2_triangulation_2<GT, Tds2> &t2)
{
  return ! (t1 == t2);
}

#include <CGAL/Periodic_2_triangulation_dummy_12.h>

} //namespace CGAL


#endif //CGAL_PERIODIC_2_TRIANGULATION_2_H
