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

#ifndef CGAL_PERIODIC_3_DELAUNAY_TRIANGULATION_3_GENERIC_H
#define CGAL_PERIODIC_3_DELAUNAY_TRIANGULATION_3_GENERIC_H

#include <CGAL/license/Periodic_3_triangulation_3.h>

#include <CGAL/internal/Generic_P3T3/Periodic_3_triangulation_vertex_base_3_generic.h>
#include <CGAL/internal/Generic_P3T3/Periodic_3_triangulation_cell_base_3_generic.h>
// #include <CGAL/internal/Generic_P3T3/Periodic_3_triangulation_iterators_3_generic.h>

#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Periodic_3_Delaunay_triangulation_3.h>
#include <CGAL/Periodic_3_Delaunay_triangulation_traits_3.h>
#include <CGAL/Lattice_3.h>
#include <CGAL/Triangulation_data_structure_3.h>

#include <CGAL/utility.h>

#include <utility>
#include <iostream>
#include <fstream>
#include <string>

// @todo steal the get_vertex(Cell_handle, int, Vertex& /*canonical*/, Offset) from P2T2_GT

namespace CGAL {

template <typename GT,
          typename TDS>
class Periodic_3_Delaunay_triangulation_3_generic
{
  typedef Periodic_3_Delaunay_triangulation_3_generic<GT, TDS>            Self;

public:
  typedef Delaunay_triangulation_3<GT, TDS>                               DT3;
  typedef Periodic_3_Delaunay_triangulation_3<GT, TDS>                    P3DT3;

public:
  typedef GT                                                              Geom_traits;
  typedef TDS                                                             Triangulation_data_structure;

  typedef typename TDS::Vertex                                            Vertex;
  typedef typename TDS::Vertex_handle                                     Vertex_handle;
  typedef typename TDS::Edge                                              Edge;
  typedef typename TDS::Facet                                             Facet;
  typedef typename TDS::Cell_handle                                       Cell_handle;

  typedef typename Vertex::Point                                          Point;
  // typedef Periodic_3_triangulation_point_iterator_3_generic<Self>         Periodic_point_iterator;

public:
  typedef typename GT::Domain                                             Lattice;

public:
  typedef typename Geom_traits::FT                                        FT;
  typedef typename Geom_traits::Segment_3                                 Segment;
  typedef typename Geom_traits::Vector_3                                  Vector;
  typedef typename Geom_traits::Triangle_3                                Triangle;

  typedef typename P3DT3::Offset                                          Offset;
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
    FACET, //2
    CELL, //4
    OUTSIDE_CONVEX_HULL, // unused, for compatibility with Alpha shapes
    OUTSIDE_AFFINE_HULL // unused, for compatibility with Alpha shapes
  };

  // enum Iterator_type
  // {
  //   IT_DT3 = 0,
  //   IT_P3T3 // 1
  // };

  /// Constructors
  Periodic_3_Delaunay_triangulation_3_generic(const GT& gt)
    : is_1_cover_(false), gt_(gt), dt3(gt_), p3dt3(gt_)
  {
    sq_circumradius_threshold = 0.0625 * gt.get_domain().systole_sq_length();
  }

  Periodic_3_Delaunay_triangulation_3_generic(const Lattice& lattice)
    : Periodic_3_Delaunay_triangulation_3_generic(GT(lattice))
  { }

  template <class InputIterator>
  Periodic_3_Delaunay_triangulation_3_generic(InputIterator first, InputIterator beyond,
                                              const GT& gt)
    : is_1_cover_(false), gt_(gt), dt3(gt_), p3dt3(gt_)
  {
    sq_circumradius_threshold = 0.0625 * gt.get_domain().systole_sq_length();
    insert(first, beyond);
  }

  template <class InputIterator>
  Periodic_3_Delaunay_triangulation_3_generic(InputIterator first, InputIterator beyond,
                                              const Lattice& lattice)
    : Periodic_3_Delaunay_triangulation_3_generic(first, beyond, GT(lattice))
  { }

  // @todo (virtual vertices cannot be naively copied)
  Periodic_3_Delaunay_triangulation_3_generic& operator=(const Periodic_3_Delaunay_triangulation_3_generic&) = delete;

  // @todo
  //  void copy_triangulation(const Periodic_3_Delaunay_triangulation_3_generic &tr) { }
  //  void swap(Periodic_3_Delaunay_triangulation_3_generic &tr) { }

  void clear()
  {
    dt3.clear();
    p3dt3.clear();

    is_1_cover_ = false;
    canonical_vertices.clear();
    cells_with_too_big_circumradius.clear();
  }

  /// Access
  const Geom_traits& geom_traits() const { return gt_; }
  const Triangulation_data_structure& tds() const
  {
    if(is_1_cover_)
      return p3dt3.tds();
    else
      return dt3.tds();
  }

  bool is_1_cover() const { return is_1_cover_; }

  size_t too_big_cells() { return cells_with_too_big_circumradius.size(); }

  /// Returns the number of vertices. Counts all vertices that are
  /// representatives of the same point in the 1-cover as one vertex.
  size_type number_of_vertices() const
  {
    if(is_1_cover_)
      return p3dt3.number_of_vertices();
    else
      return dt3.number_of_vertices() / 27;
  }

  // size_type number_of_edges() const
  // {
  //   if(is_1_cover_)
  //   {
  //     return p3dt3.number_of_edges();
  //   }
  //   else
  //   {
  //     // @todo
  //     return -1;
  //   }
  // }

  // size_type number_of_facets() const
  // {
  //   if(is_1_cover_)
  //   {
  //     return p3dt3.number_of_facets();
  //   }
  //   else
  //   {
  //     // @todo
  //     return -1;
  //   }
  // }

  size_type number_of_cells() const
  {
    if(is_1_cover_)
    {
      return p3dt3.number_of_cells();
    }
    else
    {
      // @todo something better than this naive way (flag edge/face classes)
      size_type nf = 0;
      typename DT3::Finite_cells_iterator fit = dt3.finite_cells_begin(),
                                          fend = dt3.finite_cells_end();
      for(; fit!=fend; ++fit)
        if(is_canonical(fit))
          ++nf;

      return nf;
    }
  }

  // /// Start iterator over the points
  // Periodic_point_iterator periodic_points_begin() const
  // {
  //   if(is_1_cover_)
  //     return Periodic_point_iterator(this, IT_P2T2);
  //   else
  //     return Periodic_point_iterator(this, IT_DT3);
  // }

  // /// Past-the-end iterator over the points
  // Periodic_point_iterator periodic_points_end() const
  // {
  //   if(is_1_cover_)
  //     return Periodic_point_iterator(this, 1, IT_P2T2);
  //   else
  //     return Periodic_point_iterator(this, 1, IT_DT3);
  // }

  //   Vertices_iterator vertices_begin() const { }
  //   Vertices_iterator vertices_end() const { }
  //   Edges_iterator edges_begin() const { }
  //   Edges_iterator edges_end() const { }
  //   Faces_iterator faces_begin() const { }
  //   Faces_iterator faces_end() const { }

  // (and finite versions, e.g. Finite_vertices_iterator finite_vertices_iterator, for other packages)

  //   Vertex_circulator adjacent_vertices(Vertex_handle vh) const { }
  //   Edge_circulator incident_edges(Vertex_handle vh) const { }
  //   Cell_circulator incident_cells(Vertex_handle vh) const { }

  int dimension() const { return (number_of_vertices() == 0) ? -2 : 3; }

  // @todo two versions, one using the sufficient condition, one checking for real.
  // "For real" --> P3T3 is_embeddable_in_...
  // "For real" is implemented below.
  bool is_simplicial_complex() const
  {
    // Ensure there is no edge between a vertex and itself.
    typename DT3::Finite_cells_iterator fit = dt3.cells_begin(),
                                        fend = dt3.cells_end();
    for(; fit!=fend; ++fit)
    {
      if(!is_canonical(fit))
        continue;

      Vertex_handle vh0 = canonical_vertex(fit->vertex(0));
      Vertex_handle vh1 = canonical_vertex(fit->vertex(1));
      Vertex_handle vh2 = canonical_vertex(fit->vertex(2));
      Vertex_handle vh3 = canonical_vertex(fit->vertex(3));

      if(vh0 == vh1 || vh0 == vh2 || vh0 == vh3 || vh1 == vh2 || vh1 == vh3 || vh2 == vh3)
        return false;
    }

    // Ensure there are no two edges between the same pair of vertices.
    typename DT3::Finite_vertices_iterator vit = dt3.vertices_begin(),
                                           vend = dt3.vertices_end();
    for(; vit!=vend; ++vit)
    {
      if(!is_canonical(vit))
        continue;

      typename DT3::Vertex_circulator vc = dt3.incident_vertices(vit), done = vc;
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
  FT compute_squared_circumradius(const Cell_handle fh) const
  {
    return geom_traits().compute_squared_radius_3_object()(dt3.point(fh, 0),
                                                           dt3.point(fh, 1),
                                                           dt3.point(fh, 2),
                                                           dt3.point(fh, 3));
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
    return (vh->offset() == Offset(0, 0, 0));
  }

  Offset compute_offset(const Edge& e) const
  {
    Cell_handle fh = e.first;
    int i = e.second;
    int j = e.third;
    return compute_offset(fh->vertex(i), fh->vertex(j));
  }

  bool is_canonical(const Edge e) const
  {
    if(dt3.is_infinite(e.first))
      return false;

    return compute_offset(e) == Offset(0, 0, 0);
  }

  Offset compute_offset(const Facet& e) const
  {
    Cell_handle fh = e.first;
    int i = e.second;
    return compute_offset(fh->vertex((i+1)%4), fh->vertex((i+2)%4), fh->vertex((i+3)%4));
  }

  bool is_canonical(const Facet e) const
  {
    if(dt3.is_infinite(e.first))
      return false;

    return compute_offset(e) == Offset(0, 0, 0);
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
    if(dt3.is_infinite(vh1) || dt3.is_infinite(vh2))
      return false;

    return compute_offset(vh1, vh2) == Offset(0, 0, 0);
  }

  // Offset of an facet represented by its three vertices
  Offset compute_offset(const Vertex_handle vh1, const Vertex_handle vh2, const Vertex_handle vh3) const
  {
    // @todo proper min implementation
    return min(vh1->offset(), vh2->offset(), vh3->offset());
  }

  // Canonicity of an facet represented by its three vertices
  bool is_canonical(const Vertex_handle vh1, const Vertex_handle vh2, const Vertex_handle vh3) const
  {
    if(dt3.is_infinite(vh1) || dt3.is_infinite(vh2) || dt3.is_infinite(vh3))
      return false;

    return compute_offset(vh1, vh2) == Offset(0, 0, 0);
  }

  Offset compute_offset(const Cell_handle fh) const
  {
    return min(fh->vertex(0)->offset(), fh->vertex(1)->offset(), fh->vertex(2)->offset(), fh->vertex(3)->offset());
  }

  bool is_canonical(const Cell_handle fh) const
  {
    if(dt3.is_infinite(fh))
      return false;

    return compute_offset(fh) == Offset(0, 0, 0);
  }

  template <typename ForwardCellIterator>
  void mark_canonical_cells(ForwardCellIterator fit, ForwardCellIterator beyond)
  {
    while(fit != beyond)
    {
      fit->set_canonical_flag(is_canonical(fit++));
    }
  }

  void mark_canonical_cells()
  {
    return mark_canonical_cells(dt3.finite_cells_begin(), dt3.finite_cells_end());
  }

  void mark_canonical_cells(Vertex_handle vh)
  {
    CGAL_precondition(vh != Vertex_handle());

    if(dt3.dimension() != 3)
      return;

    std::vector<Cell_handle> incident_chs;
    dt3.incident_cells(vh, std::back_inserter(incident_chs));
    for(Cell_handle ch : incident_chs)
      ch->set_canonical_flag(is_canonical(ch));
  }

  void reset_all_canonicity()
  {
    typename DT3::Finite_cells_iterator fit = dt3.cells_begin(),
                                        fend = dt3.cells_end();
    for(; fit!=fend; ++fit)
      fit->set_canonical_flag(false);
  }

  // Returns Cell_handle() if the provided cell is not a periodic cell
  Cell_handle get_canonical_cell(Cell_handle fh) const
  {
    std::cout << "Getting canonical cell of: " << dt3.point(fh, 0) << " "
                                               << dt3.point(fh, 1) << " "
                                               << dt3.point(fh, 2) << " "
                                               << dt3.point(fh, 3) << std::endl;

    return find_translated_cell(fh, -compute_offset(fh));
  }

  /// Low level functions to mascarade the DT3 as a periodic triangulation
  Point construct_barycenter(Cell_handle fh) const
  {
    Vector v0 = Vector(CGAL::ORIGIN, fh->vertex(0)->point());
    Vector v1 = Vector(CGAL::ORIGIN, fh->vertex(1)->point());
    Vector v2 = Vector(CGAL::ORIGIN, fh->vertex(2)->point());
    Vector v3 = Vector(CGAL::ORIGIN, fh->vertex(3)->point());

    Vector bcv = gt_.construct_scaled_vector_3_object()(
                   gt_.construct_sum_of_vectors_3_object()(v0,
                     gt_.construct_sum_of_vectors_3_object()(v1,
                       gt_.construct_sum_of_vectors_3_object()(v2,v3))), FT(1)/3);

    return Point(bcv.x(), bcv.y(), bcv.z());
  }

  // Given a cell having a vertex in the domain, and an offset such that
  // at least one of the translated vertices has a vertex in the domain,
  // return the handle of the translated cell.
  //
  // Returns Cell_handle() if the provided cell is not a periodic cell
  Cell_handle find_translated_cell(Cell_handle ch, const Offset& o) const
  {
    // The code commented out below does the same and is simpler,
    // but uses a construction and point location.
    //    Point translated_barycenter = gt_.get_domain().translate_by_offset(
    //        construct_barycenter(fh), o);
    //    return dt3.locate(translated_barycenter);

    // Find a vertex whose translate is in the domain.
    bool vertex_found = false;
    int j=0;
    for(; j<4; ++j)
    {
      if(ch->vertex(j)->offset() == -o)
      {
        vertex_found = true;
        break;
      }
    }

    CGAL_assertion_msg(vertex_found, "Invalid offset for translation of cell.");
    Vertex_handle cv = canonical_vertices.at(ch->vertex(j));

    // Scan through the incident cells and find the one that is
    // equivalent to fh.
    std::vector<Cell_handle> incident_chs;
    dt3.incident_cells(cv, std::back_inserter(incident_chs));
    CGAL_assertion(!incident_chs.empty());

    for(Cell_handle cch : incident_chs)
    {
      CGAL_assertion(!dt3.is_infinite(cch));

      int cj = cch->index(cv);

      bool found = false;
      for(int i=1; i<4; ++i)
      {
        found = false;
        Vertex_handle vi = ch->vertex((j+i) % 4);
        for(int ci=1; ci<4; ++ci)
        {
            Vertex_handle cvi = cch->vertex((cj+ci) % 4);
            if(canonical_vertices.at(vi) == canonical_vertices.at(cvi) && cvi->offset() == vi->offset() + o)
            {
              found = true;
              break;
            }
        }
        if(!found)
          break;
      }

      if(!found)
        continue;

      return cch;
    }

    // provided cell is not a periodic cell and doesn't have a copy with a vertex in the domain.
    return Cell_handle();
  }

  // @todo something smarter, this function is core to everything else
  // How could we keep "canonical neighbors" in memory instead of having
  // to find them later...?
  //
  // @cache this value?
  Cell_handle neighbor(Cell_handle fh, int i) const
  {
    if(is_canonical(fh->neighbor(i)))
      return fh->neighbor(i);

    // Translate the cell so that the corresponding edge is canonical.
    Offset facet_off = compute_offset(std::make_pair(fh, i));
    Cell_handle fh_trans = find_translated_cell(fh, -facet_off);
    CGAL_assertion(fh_trans != Cell_handle());

    // Find the vertex corresponding to vertex(i) in the original.
    bool vertex_found = false;
    int j=0;
    for(; j<4; ++j) {
      if(canonical_vertices.at(fh_trans->vertex(j)) == canonical_vertices.at(fh->vertex(i))
         && fh_trans->vertex(j)->offset() == fh->vertex(i)->offset() - facet_off) {
        vertex_found = true;
        break;
      }
    }
    CGAL_assertion(vertex_found);

    // Get the neighbour in DT3 and check if it's canonical.
    Cell_handle neighbor = fh_trans->neighbor(j);
    if(is_canonical(neighbor))
      return neighbor;
    else // if not, find the canonical translate.
      return get_canonical_cell(neighbor);
  }

  /// Iterators and Circulators
  // @todo this should return a (custom) circulator
  std::set<Cell_handle> incident_cells(Vertex_handle vh) const
  {
    std::set<Cell_handle> ifhs;

    // Need to start from the canonical vertex to make sure there is no "weird" (to be defined properly) cells
    if(!is_canonical(vh))
      vh = canonical_vertices.at(vh);

    typename DT3::Cell_circulator tds_fc = dt3.incident_cells(vh), done = tds_fc;
    do
    {
      ifhs.insert(get_canonical_cell(tds_fc));

      // @todo when proper periodic traits are used to construct dt3 (that is,
      // insert (p, off) instead of constructing the offsetted point,
      // then the vertices will store the _canonical_ (geometric) point,
      // and we don't have to use the canonical_vertices[] but rather we can
      // just compare fh->vertex(0)->point() and vh->point()
      //
      // "Real" offsetted point is then done via a 'construct_point(cp, off)' function
      // in the triangulation class

      // The alternative is to use the neighbourhood relation of the cells, and
      // traverse along edges around the central vertex, in the process possibly
      // shifting the center vertex around to stay canonical.
      // Incomplete code of this in a previous commit.
    }
    while(++tds_fc != done);

    return ifhs;
  }

  /// Locate functions
  //  Cell_handle locate(const Point& p, Offset& lo,
  //                     Locate_type& lt, int& li,
  //                     Cell_handle start = Cell_handle()) const
  //  { }

  //  Cell_handle locate(const Point& p,
  //                     Locate_type& lt, int& li,
  //                     Cell_handle start = Cell_handle()) const
  //  {
  //    Offset lo;
  //    return locate(p, lo, lt, li, start);
  //  }

  //  Cell_handle locate(const Point& p,
  //                     Cell_handle start = Cell_handle()) const
  //  {
  //    Locate_type lt;
  //    int li;
  //    return locate(p, lt, li, start);
  //  }

  /// Insertion and removal
  Vertex_handle insert(const Point& p)
  {
    if(is_1_cover_)
      return insert_in_p3dt3(p);
    else
      return insert_in_dt3(p);
  }

  Vertex_handle insert_in_dt3(const Point& p)
  {
    const Point cp = gt_.get_domain().construct_canonical_point(p);
    std::cout << "Insert (DT3): " << p << " canonical: " << cp << std::endl;

    std::cout << dt3.number_of_vertices() << " vertices (before insertion)" << std::endl;
    if(dt3.dimension() >= 3) // equivalent to !dt3.empty() since we insert duplicate vertices
    {
      // @todo avoid recomputing the conflict zone if possible (done also 'insert', sort of)
      Cell_handle c = dt3.locate(p);

      std::vector<Cell_handle> cells_in_conflict;
      dt3.find_conflicts(cp, c, CGAL::Emptyset_iterator(), std::back_inserter(cells_in_conflict));
      std::cout << cells_in_conflict.size() << " cells in conflict" << std::endl;

      size_t erased_cells = 0;  // @tmp
      for(Cell_handle ch : cells_in_conflict)
      {
        Cell_handle cch = get_canonical_cell(ch);

        // We might get non-periodic "boundary" cells that don't have a canonical version
        if(cch != Cell_handle())
        {
          // Some cells might appear multiple times in the conflict zone, but this is fine.
          erased_cells += cells_with_too_big_circumradius.erase(cch);
        }
      }
      std::cout << "Cells with too big radius in conflict zone:" << erased_cells << std::endl;
    }

    Vertex_handle vh = dt3.insert(cp);

    CGAL_assertion(vh != Vertex_handle());
    vh->set_offset(Offset(0, 0, 0));

    mark_canonical_cells(vh);

    canonical_vertices[vh] = vh;
    for(const std::vector<int>& off : overlapping_offsets)
    {
      // @fixme Insert without constructions (cf P3T3); factor into lattice class
      const Vector off_v = gt_.construct_sum_of_vectors_3_object()(
                             gt_.construct_scaled_vector_3_object()(gt_.get_domain().basis()[0], off[0]),
                             gt_.construct_sum_of_vectors_3_object()(
                               gt_.construct_scaled_vector_3_object()(gt_.get_domain().basis()[1], off[1]),
                               gt_.construct_scaled_vector_3_object()(gt_.get_domain().basis()[2], off[2])));

      const Point off_p = cp + off_v;
      if(gt_.get_domain().is_in_scaled_domain(off_p, 3))
      {
        Vertex_handle vh_copy = dt3.insert(off_p);
        CGAL_assertion(vh_copy != Vertex_handle());
        vh_copy->set_offset(Offset(off[0], off[1], off[2]));

        canonical_vertices[vh_copy] = vh;

        mark_canonical_cells(vh_copy);
      }
    }

    // Update the current maximum circumradius value
    std::cout << "Gather cells w/ too big circumradius" << std::endl;
    std::vector<Cell_handle> incident_chs;
    dt3.incident_cells(vh, std::back_inserter(incident_chs));

    for(Cell_handle ch : incident_chs)
    {
      CGAL_assertion(!dt3.is_infinite(ch));

      Cell_handle cfh = get_canonical_cell(ch);
      CGAL_assertion(cfh != Cell_handle());

      const FT sq_cr = compute_squared_circumradius(cfh);
      std::cout << " sq_cr: " << sq_cr << " sys:" << sq_circumradius_threshold << std::endl;
      if(sq_cr > sq_circumradius_threshold)
        cells_with_too_big_circumradius.insert(cfh);
    }

    std::cout << cells_with_too_big_circumradius.size() << " cells with too big sq_cr" << std::endl;

    if(cells_with_too_big_circumradius.empty())
       convert_to_1_cover();

    return vh;
  }

  Vertex_handle insert_in_p3dt3(const Point& p)
  {
    const Point cp = gt_.get_domain().construct_canonical_point(p);
    std::cout << "Insert (P3DT3): " << p << " canonical: " << cp << std::endl;

    return p3dt3.insert(p);
  }

  // @todo mirror it with other insert functions (e.g. spatial sorting)
  template <class InputIterator>
  void insert(InputIterator first, InputIterator beyond)
  {
    std::size_t np = std::distance(first, beyond);

    while(first != beyond)
      insert(*first++);

    CGAL_postcondition(dt3.dimension() == 3);
    CGAL_postcondition(dt3.number_of_vertices() == 27 * np);
    CGAL_postcondition(dt3.is_valid());
  }

  /// Covering conversion
  void add_facet_to_incident_cells_map(const Cell_handle fh, int i,
                                      std::map<std::set<Vertex_handle>,
                                      std::vector<std::pair<Cell_handle, int> > >& incident_cells_map)
  {
    typedef std::set<Vertex_handle>                                            Facet_vertices;
    typedef std::pair<Cell_handle, int>                                        Incident_cell;
    typedef std::map<Facet_vertices, std::vector<Incident_cell> >              Incident_cells_map;

    // the opposite vertex of f in c is i
    Facet_vertices f;
    f.insert(fh->vertex((i + 1) % 4));
    f.insert(fh->vertex((i + 2) % 4));
    f.insert(fh->vertex((i + 3) % 4));
    CGAL_precondition(f.size() == 3);

    Incident_cell icf = std::make_pair(fh, i);
    std::vector<Incident_cell> vec;
    vec.push_back(icf);

    std::pair<typename Incident_cells_map::iterator, bool> is_insert_successful =
        incident_cells_map.insert(std::make_pair(f, vec));
    if(!is_insert_successful.second) // the entry already exists in the map
    {
      // a facet must have exactly two incident cells
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

    p3dt3.clear();
    p3dt3.tds().set_dimension(3);

    std::unordered_map<Vertex_handle, Vertex_handle> vertex_correspondence_map;

    typename DT3::Finite_vertices_iterator vit = dt3.finite_vertices_begin(),
                                           vend = dt3.finite_vertices_end();
    for(; vit!=vend; ++vit)
    {
      if(!is_canonical(vit))
        continue;

      Vertex_handle vh = p3dt3.tds().create_vertex();
      vh->set_point(vit->point());
      vertex_correspondence_map[vit] = vh;
    }

    // @todo array instead of vector
    typedef std::map<std::set<Vertex_handle>, // three vertices of a facet
                     std::vector<std::pair<Cell_handle, int> > > Incident_cells_map;
    Incident_cells_map incident_cells_map;

    size_type cfc = 0;
    typename DT3::Finite_cells_iterator fit = dt3.finite_cells_begin(),
                                        fend = dt3.finite_cells_end();
    for(; fit!=fend; ++fit)
    {
      if(!is_canonical(fit))
        continue;

      ++cfc;

      Vertex_handle p3dt3_vh0 = vertex_correspondence_map.at(canonical_vertex(fit->vertex(0)));
      Vertex_handle p3dt3_vh1 = vertex_correspondence_map.at(canonical_vertex(fit->vertex(1)));
      Vertex_handle p3dt3_vh2 = vertex_correspondence_map.at(canonical_vertex(fit->vertex(2)));
      Vertex_handle p3dt3_vh3 = vertex_correspondence_map.at(canonical_vertex(fit->vertex(3)));

      Cell_handle fh = p3dt3.tds().create_cell(p3dt3_vh0, p3dt3_vh1, p3dt3_vh2, p3dt3_vh3);

      Offset min_off = CGAL::min(fit->vertex(0)->offset(),
                                 fit->vertex(1)->offset(),
                                 fit->vertex(2)->offset(),
                                 fit->vertex(3)->offset());

      fh->set_offsets(fit->vertex(0)->offset() - min_off,
                      fit->vertex(1)->offset() - min_off,
                      fit->vertex(2)->offset() - min_off,
                      fit->vertex(3)->offset() - min_off);

      add_facet_to_incident_cells_map(fh, 0, incident_cells_map);
      add_facet_to_incident_cells_map(fh, 1, incident_cells_map);
      add_facet_to_incident_cells_map(fh, 2, incident_cells_map);
      add_facet_to_incident_cells_map(fh, 3, incident_cells_map);

      // Set up incident cell information
      for(int i=0; i<4; ++i)
      {
        if(fh->vertex(i)->cell() == Cell_handle())
          fh->vertex(i)->set_cell(fh);
      }
    }

    std::cout << dt3.number_of_vertices() / 27 << " canonical vertices in dt3" << std::endl;
    std::cout << cfc << " canonical cells in dt3" << std::endl;

    // Set up adjacencies
    typename Incident_cells_map::const_iterator ifit = incident_cells_map.begin();
    for(; ifit!=incident_cells_map.end(); ++ifit)
    {
      const std::vector<std::pair<Cell_handle, int> >& adjacent_cells = ifit->second;
      CGAL_assertion(adjacent_cells.size() == 2);

      Cell_handle f0 = adjacent_cells[0].first;
      int i0 = adjacent_cells[0].second;
      Cell_handle f1 = adjacent_cells[1].first;
      int i1 = adjacent_cells[1].second;

      p3dt3.tds().set_adjacency(f0, i0, f1, i1);
    }

    std::cout << p3dt3.number_of_vertices() << " vertices in p3dt3" << std::endl;
    std::cout << p3dt3.number_of_cells() << " canonical cells in p3dt3" << std::endl;

    // std::ofstream out_dt3("dt3_post_convert.off");
    // CGAL::draw_t2(out_dt3, p3dt3);
    // out_dt3.close();
    // std::ofstream out_pc("p3dt3_post_convert.off");
    // CGAL::write_P2T2_to_OFF(out_pc, p3dt3);
    // out_pc.close();

    CGAL_postcondition(p3dt3.is_valid(true));
  }

public: // @tmp, shouldn't be exposed
  Geom_traits gt_;

  DT3 dt3;
  P3DT3 p3dt3;

private:
  bool is_1_cover_;

  std::unordered_map<Vertex_handle, Vertex_handle> canonical_vertices;

  FT sq_circumradius_threshold;
  std::set<Cell_handle> cells_with_too_big_circumradius;

  // A list of those offsets such that the domain translated along the offset
  // overlaps the scaled domain.
  // @todo: factor this into the lattice.
  // @todo: let lattice filter out non-overlapping ones
  // @todo: prove sufficiency
  const std::vector<std::vector<int> > overlapping_offsets =
  {
    // offsets that are entirely contained within the scaled domain
    {-1, -1, -1}, {0, 0, 1}, {0, 1, 0}, {1, 0, 0}, {1, 1, 0}, {1, 0, 1}, {0, -1, -1},
    {0, 1, 1}, {-1, 0, -1}, {-1, -1, 0}, {1, 1, 1}, {0, 0, -1}, {0, -1, 0}, {-1, 0, 0},
    // offsets that have a guaranteed intersection with the scaled domain
    {1, 2, 0}, {1, 0, 2}, {-1, -2, -2}, {2, 1, 0}, {2, 0, 1}, {1, -1, -1}, {0, 1, 2},
    {-2, -1, -2}, {0, 2, 1}, {-1, 1, -1}, {-2, -2, -1}, {-1, -1, 1}, {1, 1, 2},
    {-1, -1, -2}, {1, 2, 1}, {0, 1, -1}, {-1, -2, -1}, {0, -1, 1}, {2, 1, 1},
    {1, 0, -1}, {1, -1, 0}, {-2, -1, -1}, {-1, 0, 1}, {-1, 1, 0}, {1, 2, 2},
    {-1, 0, -2}, {-1, -2, 0}, {2, 1, 2}, {0, -1, -2}, {2, 2, 1}, {1, 1, -1},
    {0, -2, -1}, {1, -1, 1}, {-2, -1, 0}, {-2, 0, -1}, {-1, 1, 1},
    // offsets that might have an intersection with the scaled domain (6 of them will)
    {3, 2, 1}, {2, 1, -1}, {3, 1, 2}, {2, -1, 1}, {1, -1, -2}, {1, -2, -1},
    {2, 3, 1}, {1, 2, -1}, {1, 3, 2}, {-1, 2, 1}, {-1, 1, -2}, {-2, 1, -1},
    {2, 1, 3}, {1, -1, 2}, {1, 2, 3}, {-1, 1, 2}, {-1, -2, 1}, {-2, -1, 1},
    {-1, -2, -3}, {-1, -3, -2}, {-2, -1, -3}, {-3, -1, -2}, {-2, -3, -1}, {-3, -2, -1}
  };
};

} // namespace CGAL

#endif // CGAL_PERIODIC_3_DELAUNAY_TRIANGULATION_3_GENERIC_H
