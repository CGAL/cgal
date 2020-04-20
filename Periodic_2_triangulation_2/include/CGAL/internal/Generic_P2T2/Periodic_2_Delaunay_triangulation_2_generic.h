// Copyright (c) 2019-2020 XXXXX
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Mael Rouxel-Labbé
//                 Georg Osang

#ifndef CGAL_PERIODIC_2_DELAUNAY_TRIANGULATION_2_GENERIC_H
#define CGAL_PERIODIC_2_DELAUNAY_TRIANGULATION_2_GENERIC_H

#include <CGAL/license/Periodic_2_triangulation_2.h>

#include <CGAL/internal/Generic_P2T2/Periodic_2_triangulation_vertex_base_2_generic.h>
#include <CGAL/internal/Generic_P2T2/Periodic_2_triangulation_face_base_2_generic.h>
#include <CGAL/internal/Generic_P2T2/Periodic_2_triangulation_iterators_2_generic.h>

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Periodic_2_Delaunay_triangulation_2.h>
#include <CGAL/Periodic_2_Delaunay_triangulation_traits_2.h>
#include <CGAL/Lattice_2.h>
#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/draw_triangulation_2.h>

#include <CGAL/utility.h>

#include <utility>
#include <iostream>
#include <fstream>
#include <string>

// @todo steal the get_vertex(Face_handle, int, Vertex& /*canonical*/, Offset) from P2T2_GT

namespace CGAL {

template <typename GT,
          typename TDS>
class Periodic_2_Delaunay_triangulation_2_generic
{
  typedef Periodic_2_Delaunay_triangulation_2_generic<GT, TDS>            Self;

public:
  typedef Delaunay_triangulation_2<GT, TDS>                               DT2;
  typedef Periodic_2_Delaunay_triangulation_2<GT, TDS>                    P2DT2;

public:
  typedef GT                                                              Geom_traits;
  typedef TDS                                                             Triangulation_data_structure;

  typedef typename TDS::Vertex                                            Vertex;
  typedef typename TDS::Vertex_handle                                     Vertex_handle;
  typedef typename TDS::Edge                                              Edge;
  typedef typename TDS::Face_handle                                       Face_handle;

  typedef typename Vertex::Point                                          Point;
  typedef Periodic_2_triangulation_point_iterator_2_generic<Self>         Periodic_point_iterator;

public:
  typedef typename GT::Domain                                             Lattice;

public:
  typedef typename Geom_traits::FT                                        FT;
  typedef typename Geom_traits::Segment_2                                 Segment;
  typedef typename Geom_traits::Vector_2                                  Vector;
  typedef typename Geom_traits::Triangle_2                                Triangle;

  typedef typename P2DT2::Offset                                          Offset;
  typedef std::pair<Point, Offset>                                        Periodic_point;
  typedef array<std::pair<Point, Offset>, 2>                              Periodic_segment;
  typedef array<std::pair<Point, Offset>, 3>                              Periodic_triangle;
  typedef array<std::pair<Point, Offset>, 4>                              Periodic_tetrahedron;

  typedef typename Triangulation_data_structure::size_type                size_type;

  //Tag to distinguish Delaunay from regular triangulations
  typedef Tag_false                                                       Weighted_tag;

  // Tag to distinguish periodic triangulations from others
  typedef Tag_true                                                        Periodic_tag;

  enum Locate_type
  {
    VERTEX = 0,
    EDGE, //1
    FACE, //2
    EMPTY , //4
    OUTSIDE_CONVEX_HULL, // unused, for compatibility with Alpha shapes
    OUTSIDE_AFFINE_HULL // unused, for compatibility with Alpha shapes
  };

  enum Iterator_type
  {
    IT_DT2 = 0,
    IT_P2T2 // 1
  };

  /// Constructors
  Periodic_2_Delaunay_triangulation_2_generic(const GT& gt)
    : is_1_cover_(false), gt_(gt), dt2(gt_), p2dt2(gt_)
  {
    sq_circumradius_threshold = 0.0625 * gt.get_domain().systole_sq_length();
  }

  Periodic_2_Delaunay_triangulation_2_generic(const Lattice& lattice)
    : Periodic_2_Delaunay_triangulation_2_generic(GT(lattice))
  { }

  template <class InputIterator>
  Periodic_2_Delaunay_triangulation_2_generic(InputIterator first, InputIterator beyond,
                                              const GT& gt)
    : is_1_cover_(false), gt_(gt), dt2(gt_), p2dt2(gt_)
  {
    sq_circumradius_threshold = 0.0625 * gt.get_domain().systole_sq_length();
    insert(first, beyond);
  }

  template <class InputIterator>
  Periodic_2_Delaunay_triangulation_2_generic(InputIterator first, InputIterator beyond,
                                              const Lattice& lattice)
    : Periodic_2_Delaunay_triangulation_2_generic(first, beyond, GT(lattice))
  { }

  // @todo (virtual vertices cannot be naively copied)
  Periodic_2_Delaunay_triangulation_2_generic& operator=(const Periodic_2_Delaunay_triangulation_2_generic&) = delete;

  // @todo
  //  void copy_triangulation(const Periodic_2_Delaunay_triangulation_2_generic &tr) { }
  //  void swap(Periodic_2_Delaunay_triangulation_2_generic &tr) { }

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
    if(is_1_cover_)
      return p2dt2.tds();
    else
      return dt2.tds();
  }

  bool is_1_cover() const { return is_1_cover_; }

  size_t too_big_faces() { return faces_with_too_big_circumradius.size(); }

  /// Returns the number of vertices. Counts all vertices that are
  /// representatives of the same point in the 1-cover as one vertex.
  size_type number_of_vertices() const
  {
    if(is_1_cover_)
      return p2dt2.number_of_vertices();
    else
      return dt2.number_of_vertices() / 9;
  }

  size_type number_of_edges() const
  {
    if(is_1_cover_)
    {
      return p2dt2.number_of_edges();
    }
    else
    {
      // Exploiting Euler's formula that #f - #e + #v = 0 on the torus
      return number_of_faces() + number_of_vertices();

      //   // Alternative naive implementation
      //   size_type ne = 0;
      //   typename DT2::Finite_edges_iterator eit = dt2.finite_edges_begin(),
      //                                       eend = dt2.finite_edges_end();
      //   for(; eit!=eend; ++eit)
      //     if(is_canonical(*eit))
      //       ++ne;

      //   return ne;
    }
  }

  size_type number_of_faces() const
  {
    if(is_1_cover_)
    {
      return p2dt2.number_of_faces();
    }
    else
    {
      // @todo something better than this naive way (flag edge/face classes)
      size_type nf = 0;
      typename DT2::Finite_faces_iterator fit = dt2.finite_faces_begin(),
                                          fend = dt2.finite_faces_end();
      for(; fit!=fend; ++fit)
        if(is_canonical(fit))
          ++nf;

      return nf;
    }
  }

  /// Start iterator over the points
  Periodic_point_iterator periodic_points_begin() const
  {
    if(is_1_cover_)
      return Periodic_point_iterator(this, IT_P2T2);
    else
      return Periodic_point_iterator(this, IT_DT2);
  }

  /// Past-the-end iterator over the points
  Periodic_point_iterator periodic_points_end() const
  {
    if(is_1_cover_)
      return Periodic_point_iterator(this, 1, IT_P2T2);
    else
      return Periodic_point_iterator(this, 1, IT_DT2);
  }

  //   Vertices_iterator vertices_begin() const { }
  //   Vertices_iterator vertices_end() const { }
  //   Edges_iterator edges_begin() const { }
  //   Edges_iterator edges_end() const { }
  //   Faces_iterator faces_begin() const { }
  //   Faces_iterator faces_end() const { }

  // (and finite versions, e.g. Finite_vertices_iterator finite_vertices_iterator, for other packages)

  //   Vertex_circulator adjacent_vertices(Vertex_handle vh) const { }
  //   Edge_circulator incident_edges(Vertex_handle vh) const { }
  //   Face_circulator incident_faces(Vertex_handle vh) const { }

  int dimension() const { return (number_of_vertices() == 0) ? -2 : 2; }

  // @todo two versions, one using the sufficient condition, one checking for real.
  // "For real" --> P3T3 is_embeddable_in_...
  // "For real" is implemented below.
  bool is_simplicial_complex() const
  {
    // Ensure there is no edge between a vertex and itself.
    typename DT2::Finite_faces_iterator fit = dt2.faces_begin(),
                                        fend = dt2.faces_end();
    for(; fit!=fend; ++fit)
    {
      if(!is_canonical(fit))
        continue;

      Vertex_handle vh0 = canonical_vertex(fit->vertex(0));
      Vertex_handle vh1 = canonical_vertex(fit->vertex(1));
      Vertex_handle vh2 = canonical_vertex(fit->vertex(2));

      if(vh0 == vh1 || vh0 == vh2 || vh1 == vh2)
        return false;
    }

    // Ensure there are no two edges between the same pair of vertices.
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
          return false; // Some neighbouring vertex appeared multiple times
        neighbours.insert(cv);
      }
      while(++vc != done);
    }

    return true;
  }

  /// Predicates
  FT compute_squared_circumradius(const Face_handle fh) const
  {
    return geom_traits().compute_squared_radius_2_object()(dt2.point(fh, 0),
                                                           dt2.point(fh, 1),
                                                           dt2.point(fh, 2));
  }

  /// Constructions
  // @todo
  Point construct_point(const Point& /*p*/, const Offset& /*off*/) const { CGAL_assertion(false); }
  std::pair<Point /*canonical point*/, Offset> periodic_point(const Point& /*p*/) const { CGAL_assertion(false); }

  /// Canonicalization
  Point construct_canonical_point(const Point& p) const
  {
    return gt_.get_domain().construct_canonical_point(p);
  }

  /// Canonicity functions
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

  bool is_canonical(const Point& p) const
  {
    return gt_.get_domain().is_in_scaled_domain(p, 1);
  }

  Vertex_handle canonical_vertex(const Vertex_handle vh) const
  {
    return (is_canonical(vh) ? vh : canonical_vertices.at(vh));
  }

  bool is_canonical(const Vertex_handle vh) const
  {
    return (vh->offset() == Offset(0, 0));
  }

  Offset compute_offset(const Edge& e) const
  {
    Face_handle fh = e.first;
    int i = e.second;
    return compute_offset(fh->vertex(dt2.cw(i)), fh->vertex(dt2.ccw(i)));
  }

  bool is_canonical(const Edge e) const
  {
    if(dt2.is_infinite(e.first))
      return false;

    return compute_offset(e) == Offset(0, 0);
  }

  // Offset of an edge represented by its two vertices
  Offset compute_offset(const Vertex_handle vh1, const Vertex_handle vh2) const
  {
    // @todo proper min implementation
    return min(vh1->offset(), vh2->offset());
  }

  // Canonicity of an edge represented by its two vertices
  bool is_canonical(const Vertex_handle vh1, const Vertex_handle vh2) const
  {
    if(dt2.is_infinite(vh1) || dt2.is_infinite(vh2))
      return false;

    return compute_offset(vh1, vh2) == Offset(0, 0);
  }

  Offset compute_offset(const Face_handle fh) const
  {
    return min(fh->vertex(0)->offset(), fh->vertex(1)->offset(), fh->vertex(2)->offset());
  }

  bool is_canonical(const Face_handle fh) const
  {
    if(dt2.is_infinite(fh))
      return false;

    return compute_offset(fh) == Offset(0, 0);
  }

  template <typename ForwardFaceIterator>
  void mark_canonical_faces(ForwardFaceIterator fit, ForwardFaceIterator beyond)
  {
    while(fit != beyond)
    {
      fit->set_canonical_flag(is_canonical(fit++));
    }
  }

  void mark_canonical_faces()
  {
    return mark_canonical_faces(dt2.finite_faces_begin(), dt2.finite_faces_end());
  }

  void mark_canonical_faces(Vertex_handle vh)
  {
    CGAL_precondition(vh != Vertex_handle());

    if(dt2.dimension() != 2)
      return;

    typename DT2::Face_circulator fc = dt2.incident_faces(vh), done = fc;
    do
    {
      fc->set_canonical_flag(is_canonical(fc));
    }
    while(++fc != done);
  }

  void reset_all_canonicity()
  {
    typename DT2::Finite_faces_iterator fit = dt2.faces_begin(),
                                        fend = dt2.faces_end();
    for(; fit!=fend; ++fit)
      fit->set_canonical_flag(false);
  }

  // Returns Face_handle() if the provided face is not a periodic face
  Face_handle get_canonical_face(Face_handle fh) const
  {
    std::cout << "Getting canonical face of: " << dt2.point(fh, 0) << " "
                                               << dt2.point(fh, 1) << " "
                                               << dt2.point(fh, 2) << std::endl;

    return find_translated_face(fh, -compute_offset(fh));
  }

  /// Low level functions to mascarade the DT2 as a periodic triangulation
  Point construct_barycenter(Face_handle fh) const
  {
    Vector v0 = Vector(CGAL::ORIGIN, fh->vertex(0)->point());
    Vector v1 = Vector(CGAL::ORIGIN, fh->vertex(1)->point());
    Vector v2 = Vector(CGAL::ORIGIN, fh->vertex(2)->point());

    Vector bcv = gt_.construct_scaled_vector_2_object()(
                   gt_.construct_sum_of_vectors_2_object()(
                     v0, gt_.construct_sum_of_vectors_2_object()(v1,v2)), FT(1)/3);

    return Point(bcv.x(), bcv.y());
  }

  // Given a face having a vertex in the domain, and an offset such that
  // at least one of the translated vertices has a vertex in the domain,
  // return the handle of the translated face.
  //
  // Returns Face_handle() if the provided face is not a periodic face
  Face_handle find_translated_face(Face_handle fh, const Offset& o) const
  {
    // The code commented out below does the same and is simpler,
    // but uses a construction and point location.
    //    Point translated_barycenter = gt_.get_domain().translate_by_offset(
    //        construct_barycenter(fh), o);
    //    return dt2.locate(translated_barycenter);

    // Find a vertex whose translate is in the domain.
    bool vertex_found = false;
    int j=0;
    for(; j<3; ++j)
    {
      if(fh->vertex(j)->offset() == -o)
      {
        vertex_found = true;
        break;
      }
    }

    CGAL_assertion_msg(vertex_found, "Invalid offset for translation of face.");
    Vertex_handle cv = canonical_vertices.at(fh->vertex(j));
    Vertex_handle v_ccw = fh->vertex(dt2.ccw(j));
    Vertex_handle v_cw = fh->vertex(dt2.cw(j));

    // Scan through the incident faces and find the one that is
    // equivalent to fh.
    typename DT2::Face_circulator fc = dt2.incident_faces(cv), done = fc;
    do
    {
      if(dt2.is_infinite(fc)) // shouldn't ever happen
        continue;

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

  // @todo something smarter, this function is core to everything else
  // How could we keep "canonical neighbors" in memory instead of having
  // to find them later...?
  //
  // @cache this value?
  Face_handle neighbor(Face_handle fh, int i) const
  {
    if(is_canonical(fh->neighbor(i)))
      return fh->neighbor(i);

    // Translate the face so that the corresponding edge is canonical.
    Offset edge_off = compute_offset(std::make_pair(fh, i));
    Face_handle fh_trans = find_translated_face(fh, -edge_off);
    CGAL_assertion(fh_trans != Face_handle());

    // Find the vertex corresponding to vertex(i) in the original.
    bool vertex_found = false;
    int j=0;
    for(; j<3; ++j) {
      if(canonical_vertices.at(fh_trans->vertex(j)) == canonical_vertices.at(fh->vertex(i))
         && fh_trans->vertex(j)->offset() == fh->vertex(i)->offset() - edge_off) {
        vertex_found = true;
        break;
      }
    }
    CGAL_assertion(vertex_found);

    // Get the neighbour in DT2 and check if it's canonical.
    Face_handle neighbor = fh_trans->neighbor(j);
    if(is_canonical(neighbor))
      return neighbor;
    else // if not, find the canonical translate.
      return get_canonical_face(neighbor);
  }

  /// Iterators and Circulators
  // @todo this should return a (custom) circulator
  std::set<Face_handle> incident_faces(Vertex_handle vh) const
  {
    std::set<Face_handle> ifhs;

    // Need to start from the canonical vertex to make sure there is no "weird" (to be defined properly) faces
    if(!is_canonical(vh))
      vh = canonical_vertices.at(vh);

    typename DT2::Face_circulator tds_fc = dt2.incident_faces(vh), done = tds_fc;
    do
    {
      ifhs.insert(get_canonical_face(tds_fc));

      // @todo when proper periodic traits are used to construct dt2 (that is,
      // insert (p, off) instead of constructing the offsetted point,
      // then the vertices will store the _canonical_ (geometric) point,
      // and we don't have to use the canonical_vertices[] but rather we can
      // just compare fh->vertex(0)->point() and vh->point()
      //
      // "Real" offsetted point is then done via a 'construct_point(cp, off)' function
      // in the triangulation class

      // The alternative is to use the neighbourhood relation of the faces, and
      // traverse along edges around the central vertex, in the process possibly
      // shifting the center vertex around to stay canonical.
      // Incomplete code of this in a previous commit.
    }
    while(++tds_fc != done);

    return ifhs;
  }

  /// Locate functions
  //  Face_handle locate(const Point& p, Offset& lo,
  //                     Locate_type& lt, int& li,
  //                     Face_handle start = Face_handle()) const
  //  { }

  //  Face_handle locate(const Point& p,
  //                     Locate_type& lt, int& li,
  //                     Face_handle start = Face_handle()) const
  //  {
  //    Offset lo;
  //    return locate(p, lo, lt, li, start);
  //  }

  //  Face_handle locate(const Point& p,
  //                     Face_handle start = Face_handle()) const
  //  {
  //    Locate_type lt;
  //    int li;
  //    return locate(p, lt, li, start);
  //  }

  /// Insertion and removal
  Vertex_handle insert(const Point& p)
  {
    if(is_1_cover_)
      return insert_in_p2dt2(p);
    else
      return insert_in_dt2(p);
  }

  Vertex_handle insert_in_dt2(const Point& p)
  {
    const Point cp = gt_.get_domain().construct_canonical_point(p);
    std::cout << "Insert (DT2): " << p << " canonical: " << cp << std::endl;

    std::cout << dt2.number_of_vertices() << " vertices" << std::endl;
    if(dt2.dimension() >= 2) // equivalent to !dt2.empty() since we insert duplicate vertices
    {
      // @todo avoid recomputing the conflict zone if possible (done also 'insert', sort of)
      std::vector<Face_handle> faces_in_conflict;
      dt2.get_conflicts(cp, std::back_inserter(faces_in_conflict));
      std::cout << faces_in_conflict.size() << " faces in conflict" << std::endl;

      size_t erased_faces = 0;
      for(Face_handle fh : faces_in_conflict)
      {
        Face_handle cfh = get_canonical_face(fh);

        // We might get non-periodic "boundary" faces that don't have a canonical version
        if(cfh != Face_handle())
        {
          // Some faces might appear multiple times in the conflict zone, but this is fine.
          erased_faces += faces_with_too_big_circumradius.erase(cfh);
        }
      }
      std::cout << "Faces with too big radius in conflict zone:" << erased_faces << std::endl;
    }

    Vertex_handle vh = dt2.insert(cp);

    CGAL_assertion(vh != Vertex_handle());
    vh->set_offset(Offset(0, 0));

    mark_canonical_faces(vh);

    canonical_vertices[vh] = vh;
    for(const std::vector<int>& off : overlapping_offsets)
    {
      // @fixme Insert without constructions (cf P3T3)
      const Vector off_v = gt_.construct_sum_of_vectors_2_object()(
                             gt_.construct_scaled_vector_2_object()(gt_.get_domain().basis()[0], off[0]),
                             gt_.construct_scaled_vector_2_object()(gt_.get_domain().basis()[1], off[1]));

      const Point off_p = cp + off_v;
      if(gt_.get_domain().is_in_scaled_domain(off_p, 3))
      {
        Vertex_handle vh_copy = dt2.insert(off_p);
        CGAL_assertion(vh_copy != Vertex_handle());
        vh_copy->set_offset(Offset(off[0], off[1]));

        canonical_vertices[vh_copy] = vh;

        mark_canonical_faces(vh_copy);
      }
    }

    // Update the current maximum circumradius value
    std::cout << "Gather faces w/ too big circumradius" << std::endl;
    typename DT2::Face_circulator fc = dt2.incident_faces(vh), done(fc);
    do
    {
      CGAL_assertion(!dt2.is_infinite(fc));

      Face_handle cfh = get_canonical_face(fc);
      CGAL_assertion(cfh != Face_handle());

      const FT sq_cr = compute_squared_circumradius(cfh);
      std::cout << " sq_cr: " << sq_cr << " sys:" << sq_circumradius_threshold << std::endl;
      if(sq_cr > sq_circumradius_threshold)
        faces_with_too_big_circumradius.insert(cfh);
    }
    while(++fc != done);

    std::cout << faces_with_too_big_circumradius.size() << " faces with too big sq_cr" << std::endl;

    if(faces_with_too_big_circumradius.empty())
       convert_to_1_cover();

    return vh;
  }

  Vertex_handle insert_in_p2dt2(const Point& p)
  {
    const Point cp = gt_.get_domain().construct_canonical_point(p);
    std::cout << "Insert (P2DT2): " << p << " canonical: " << cp << std::endl;

    return p2dt2.insert(p);
  }

  // @todo mirror it with other insert functions (e.g. spatial sorting)
  template <class InputIterator>
  void insert(InputIterator first, InputIterator beyond)
  {
    std::size_t np = std::distance(first, beyond);

    while(first != beyond)
      insert(*first++);

    CGAL_postcondition(dt2.dimension() == 2);
    CGAL_postcondition(dt2.number_of_vertices() == 9 * np);
    CGAL_postcondition(dt2.is_valid());
  }

  /// Covering conversion
  void add_edge_to_incident_faces_map(const Face_handle fh, int i,
                                      std::map<std::set<Vertex_handle>,
                                      std::vector<std::pair<Face_handle, int> > >& incident_faces_map)
  {
    typedef std::set<Vertex_handle>                                            Edge_vertices;
    typedef std::pair<Face_handle, int>                                        Incident_face;
    typedef std::map<Edge_vertices, std::vector<Incident_face> >               Incident_faces_map;

    // the opposite vertex of f in c is i
    Edge_vertices e;
    e.insert(fh->vertex((i + 1) % 3));
    e.insert(fh->vertex((i + 2) % 3));
    CGAL_precondition(e.size() == 2);

    Incident_face icf = std::make_pair(fh, i);
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
    if(is_1_cover_)
      return;

    std::cout << "Converting..." << std::endl;

    is_1_cover_ = true;

    p2dt2.clear();
    p2dt2.tds().set_dimension(2);

    std::unordered_map<Vertex_handle, Vertex_handle> vertex_correspondence_map;

    typename DT2::Finite_vertices_iterator vit = dt2.finite_vertices_begin(),
                                           vend = dt2.finite_vertices_end();
    for(; vit!=vend; ++vit)
    {
      if(!is_canonical(vit))
        continue;

      Vertex_handle vh = p2dt2.tds().create_vertex();
      vh->set_point(vit->point());
      vertex_correspondence_map[vit] = vh;
    }

    // @todo array instead of vector
    typedef std::map<std::set<Vertex_handle>, // two vertices of an edge
                     std::vector<std::pair<Face_handle, int> > > Incident_faces_map;
    Incident_faces_map incident_faces_map;

    size_type cfc = 0;
    typename DT2::Finite_faces_iterator fit = dt2.finite_faces_begin(),
                                        fend = dt2.finite_faces_end();
    for(; fit!=fend; ++fit)
    {
      if(!is_canonical(fit))
        continue;

      ++cfc;

      Vertex_handle p2dt2_vh0 = vertex_correspondence_map.at(canonical_vertex(fit->vertex(0)));
      Vertex_handle p2dt2_vh1 = vertex_correspondence_map.at(canonical_vertex(fit->vertex(1)));
      Vertex_handle p2dt2_vh2 = vertex_correspondence_map.at(canonical_vertex(fit->vertex(2)));

      Face_handle fh = p2dt2.tds().create_face(p2dt2_vh0, p2dt2_vh1, p2dt2_vh2);

      Offset min_off = CGAL::min(fit->vertex(0)->offset(),
                                 fit->vertex(1)->offset(),
                                 fit->vertex(2)->offset());

      fh->set_offsets(fit->vertex(0)->offset() - min_off,
                      fit->vertex(1)->offset() - min_off,
                      fit->vertex(2)->offset() - min_off);

      add_edge_to_incident_faces_map(fh, 0, incident_faces_map);
      add_edge_to_incident_faces_map(fh, 1, incident_faces_map);
      add_edge_to_incident_faces_map(fh, 2, incident_faces_map);

      // Set up incident face information
      for(int i=0; i<3; ++i)
      {
        if(fh->vertex(i)->face() == Face_handle())
          fh->vertex(i)->set_face(fh);
      }
    }

    std::cout << dt2.number_of_vertices() / 9 << " canonical vertices in dt2" << std::endl;
    std::cout << cfc << " canonical faces in dt2" << std::endl;

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

    std::cout << p2dt2.number_of_vertices() << " vertices in p2dt2" << std::endl;
    std::cout << p2dt2.number_of_faces() << " canonical faces in p2dt2" << std::endl;

    std::ofstream out_dt2("dt2_post_convert.off");
    CGAL::draw_t2(out_dt2, p2dt2);
    out_dt2.close();
    std::ofstream out_pc("p2dt2_post_convert.off");
    CGAL::write_P2T2_to_OFF(out_pc, p2dt2);
    out_pc.close();

    CGAL_postcondition(p2dt2.is_valid(true));
  }

public: // @tmp, shouldn't be exposed
  Geom_traits gt_;

  DT2 dt2;
  P2DT2 p2dt2;

private:
  bool is_1_cover_;

  std::unordered_map<Vertex_handle, Vertex_handle> canonical_vertices;

  FT sq_circumradius_threshold;
  std::set<Face_handle> faces_with_too_big_circumradius;

  // A list of those offsets such that the domain translated along the offset
  // overlaps the scaled domain.
  // @todo: factor this into the lattice.
  const std::vector<std::vector<int> > overlapping_offsets =
  {
    // entirely contained in scaled domains
    {-1, -1}, {0, 1}, {1, 0}, {-1, 0}, {0, -1}, {1, 1},
    // intersecting the scaled domain
    {-1, -2}, {1, 2}, {-2, -1}, {2, 1}, {-1, 1}, {1, -1}
  };
};

} // namespace CGAL

#endif // CGAL_PERIODIC_2_DELAUNAY_TRIANGULATION_2_GENERIC_H
