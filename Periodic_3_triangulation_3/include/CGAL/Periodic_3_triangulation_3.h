// Copyright (c) 1999-2003,2006-2009,2014-2015,2017   INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//                 Sylvain Pion <Sylvain.Pion@sophia.inria.fr>
//                 Nico Kruithof <Nico.Kruithof@sophia.inria.fr>
//                 Manuel Caroli <Manuel.Caroli@sophia.inria.fr>
//                 Aymeric Pelle <Aymeric.Pelle@sophia.inria.fr>
//                 Mael Rouxel-Labbé

#ifndef CGAL_PERIODIC_3_TRIANGULATION_3_H
#define CGAL_PERIODIC_3_TRIANGULATION_3_H

#include <CGAL/license/Periodic_3_triangulation_3.h>

#include <CGAL/basic.h>

#include <CGAL/internal/Periodic_3_triangulation_iterators_3.h>
#include <CGAL/Periodic_3_triangulation_ds_cell_base_3.h>
#include <CGAL/Periodic_3_triangulation_ds_vertex_base_3.h>
#include <CGAL/Periodic_3_triangulation_traits_3.h>
#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Triangulation_cell_base_3.h>
#include <CGAL/Triangulation_vertex_base_3.h>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/internal/canonicalize_helper.h>

#include <CGAL/array.h>
#include <CGAL/internal/Exact_type_selector.h>
#include <CGAL/NT_converter.h>
#include <CGAL/Unique_hash_map.h>
#include <CGAL/use.h>

#ifndef CGAL_NO_STRUCTURAL_FILTERING
#include <CGAL/internal/Static_filters/tools.h>
#include <CGAL/Triangulation_structural_filtering_traits.h>
#include <CGAL/determinant.h>
#endif // no CGAL_NO_STRUCTURAL_FILTERING

#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_smallint.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/utility/result_of.hpp>
#include <boost/unordered_map.hpp>

#include <iostream>
#include <algorithm>
#include <cmath>
#include <functional>
#include <list>
#include <utility>

namespace CGAL {

template < class GT, class TDS > class Periodic_3_triangulation_3;

template < class GT, class TDS > std::istream& operator>>
    (std::istream& is, Periodic_3_triangulation_3<GT,TDS>& tr);
template < class GT, class TDS > std::ostream& operator<<
    (std::ostream& os, const Periodic_3_triangulation_3<GT,TDS>& tr);

#ifndef CGAL_NO_STRUCTURAL_FILTERING
namespace internal {
// structural filtering is performed only for EPIC
struct Periodic_structural_filtering_3_tag {};
struct No_periodic_structural_filtering_3_tag {};

template <bool filter>
struct Periodic_structural_filtering_selector_3
{
#ifdef FORCE_STRUCTURAL_FILTERING
  typedef Periodic_structural_filtering_3_tag  Tag;
#else
  typedef No_periodic_structural_filtering_3_tag  Tag;
#endif
};

template <>
struct Periodic_structural_filtering_selector_3<true>
{
  typedef Periodic_structural_filtering_3_tag  Tag;
};
}
#endif // no CGAL_NO_STRUCTURAL_FILTERING

/**\class Periodic_3_triangulation_3
 *
 * \brief Implements functionality for computing in periodic space.
 *
 * There are several things that are special to computing in $\mathbb{T}^3$
 * such as
 * - periodicity --> offsets
 * - no infinite vertex
 * - no degenerate dimensions
 * All functions that are affected can be found in this class. When it is
 * necessary to provide different implementations for Delaunay and regular
 * triangulation, we work with visitors.
 */

template < class GT,
           class TDS = Triangulation_data_structure_3 <
                         Triangulation_vertex_base_3<GT,
                           Periodic_3_triangulation_ds_vertex_base_3<> >,
                         Triangulation_cell_base_3<GT,
                           Periodic_3_triangulation_ds_cell_base_3<> > > >
class Periodic_3_triangulation_3
  : public Triangulation_utils_3
{
  friend std::istream& operator>> <>
  (std::istream& is, Periodic_3_triangulation_3<GT, TDS>& tr);

  typedef Periodic_3_triangulation_3<GT,TDS>   Self;

public:
  typedef GT                                   Geometric_traits;
  typedef TDS                                  Triangulation_data_structure;

  typedef typename GT::Periodic_3_offset_3     Offset;
  typedef typename GT::Iso_cuboid_3            Iso_cuboid;
  typedef std::array<int, 3>           Covering_sheets;

  // point types
  typedef typename TDS::Vertex::Point          Point;
  typedef typename GT::Point_3                 Point_3;
  typedef std::pair<Point, Offset>             Periodic_point;
  typedef std::pair<Point_3, Offset>           Periodic_point_3;

  typedef typename GT::Segment_3               Segment;
  typedef typename GT::Triangle_3              Triangle;
  typedef typename GT::Tetrahedron_3           Tetrahedron;

  typedef std::array<Periodic_point, 2> Periodic_segment;
  typedef std::array<Periodic_point, 3> Periodic_triangle;
  typedef std::array<Periodic_point, 4> Periodic_tetrahedron;

  typedef std::array<Periodic_point_3, 2> Periodic_segment_3;
  typedef std::array<Periodic_point_3, 3> Periodic_triangle_3;
  typedef std::array<Periodic_point_3, 4> Periodic_tetrahedron_3;

  typedef typename TDS::Vertex                 Vertex;
  typedef typename TDS::Cell                   Cell;
  typedef typename TDS::Facet                  Facet;
  typedef typename TDS::Edge                   Edge;

  typedef typename TDS::Vertex_handle          Vertex_handle;
  typedef typename TDS::Cell_handle            Cell_handle;

  typedef typename TDS::size_type              size_type;
  typedef typename TDS::difference_type        difference_type;

  typedef typename TDS::Cell_iterator          Cell_iterator;
  typedef typename TDS::Facet_iterator         Facet_iterator;
  typedef typename TDS::Edge_iterator          Edge_iterator;
  typedef typename TDS::Vertex_iterator        Vertex_iterator;

  typedef typename TDS::Cell_circulator        Cell_circulator;
  typedef typename TDS::Facet_circulator       Facet_circulator;

  typedef Cell_iterator                        All_cells_iterator;
  typedef Facet_iterator                       All_facets_iterator;
  typedef Edge_iterator                        All_edges_iterator;
  typedef Vertex_iterator                      All_vertices_iterator;
  typedef Periodic_3_triangulation_unique_vertex_iterator_3<Self>
                                               Unique_vertex_iterator;

private:
  typedef typename GT::FT                      FT;
  typedef std::pair< Vertex_handle, Offset >   Virtual_vertex;
  typedef boost::unordered_map<Vertex_handle, Virtual_vertex>
                                               Virtual_vertex_map;
  typedef typename Virtual_vertex_map::const_iterator
                                               Virtual_vertex_map_it;
  typedef boost::unordered_map<Vertex_handle, std::vector<Vertex_handle > >
                                               Virtual_vertex_reverse_map;
  typedef typename Virtual_vertex_reverse_map::const_iterator
                                               Virtual_vertex_reverse_map_it;
  typedef Triple< Vertex_handle, Vertex_handle, Vertex_handle >
                                               Vertex_triple;

public:
  typedef Periodic_3_triangulation_tetrahedron_iterator_3<Self>
                                               Periodic_tetrahedron_iterator;
  typedef Periodic_3_triangulation_triangle_iterator_3<Self>
                                               Periodic_triangle_iterator;
  typedef Periodic_3_triangulation_segment_iterator_3<Self>
                                               Periodic_segment_iterator;
  typedef Periodic_3_triangulation_point_iterator_3<Self>
                                               Periodic_point_iterator;

  typedef Point                                value_type;
  typedef const value_type&                    const_reference;

  //Tag to distinguish regular triangulations from others
  typedef Tag_false                            Weighted_tag;

  // Tag to distinguish periodic triangulations from others
  typedef Tag_true                             Periodic_tag;

public:
  enum Iterator_type {
    STORED = 0,
    UNIQUE, //1
    STORED_COVER_DOMAIN, //2
    UNIQUE_COVER_DOMAIN };//3

  enum Locate_type {
    VERTEX = 0,
    EDGE, //1
    FACET, //2
    CELL, //3
    EMPTY , //4
    OUTSIDE_CONVEX_HULL, // unused, for compatibility with Alpha_shape_3
    OUTSIDE_AFFINE_HULL }; // unused, for compatibility with Alpha_shape_3

  // unused and undocumented types and functions required to be
  // compatible with Alpha_shape_3 / Periodic_3_mesh_3
public:
  typedef Cell_iterator                        Finite_cells_iterator;
  typedef Facet_iterator                       Finite_facets_iterator;
  typedef Edge_iterator                        Finite_edges_iterator;
  typedef Vertex_iterator                      Finite_vertices_iterator;

  int dimension() const { return (number_of_vertices() == 0) ? -2 : 3; }

  template<class T>
  bool is_infinite(const T&, int = 0, int = 0) const { return false; }

  size_type number_of_finite_cells() const { return number_of_cells(); }
  size_type number_of_finite_facets() const { return number_of_facets(); }
  size_type number_of_finite_edges() const { return number_of_edges(); }
  size_type number_of_finite_vertices() const { return number_of_vertices(); }

private:
  Geometric_traits  _gt;
  Triangulation_data_structure _tds;

protected:
  /// map of offsets for periodic copies of vertices
  Virtual_vertex_map virtual_vertices;
  Virtual_vertex_reverse_map  virtual_vertices_reverse;

protected:
  /// v_offsets temporarily stores all the vertices on the border of a
  /// conflict region.
  mutable std::vector<Vertex_handle> v_offsets;

  /// Determines if we currently compute in 3-cover or 1-cover.
  Covering_sheets _cover;

public:
  /** @name Creation */
  Periodic_3_triangulation_3(const Iso_cuboid& domain = Iso_cuboid(0,0,0, 1,1,1),
                             const Geometric_traits& gt = Geometric_traits())
    : _gt(gt), _tds()
  {
    typedef typename internal::Exact_field_selector<FT>::Type EFT;
    typedef NT_converter<FT,EFT> NTC;
    CGAL_USE_TYPE(NTC);
    CGAL_triangulation_precondition_code( NTC ntc; )
    CGAL_triangulation_precondition(ntc(domain.xmax())-ntc(domain.xmin())
                                    == ntc(domain.ymax())-ntc(domain.ymin()));
    CGAL_triangulation_precondition(ntc(domain.ymax())-ntc(domain.ymin())
                                    == ntc(domain.zmax())-ntc(domain.zmin()));
    CGAL_triangulation_precondition(ntc(domain.zmax())-ntc(domain.zmin())
                                    == ntc(domain.xmax())-ntc(domain.xmin()));

    _gt.set_domain(domain);
    _cover = CGAL::make_array(3,3,3);
    init_tds();
  }

protected:
  // Copy constructor helpers
  class Finder;

public:
  // Copy constructor duplicates vertices and cells
  Periodic_3_triangulation_3(const Periodic_3_triangulation_3& tr)
    : _gt(tr.geom_traits()),
      _cover(tr.number_of_sheets())
  {
    if(is_1_cover())
      tds() = tr.tds();
    else
      copy_multiple_covering(tr);

    CGAL_triangulation_expensive_postcondition(*this == tr);
    CGAL_triangulation_expensive_postcondition(is_valid());
  }

  virtual ~Periodic_3_triangulation_3() {}

  void copy_multiple_covering(const Periodic_3_triangulation_3& tr)
  {
    // Write the respective offsets in the vertices to make them
    // automatically copy with the tds.
    for(Vertex_iterator vit = tr.vertices_begin(); vit != tr.vertices_end(); ++vit)
      vit->set_offset(tr.get_offset(vit));

    // copy the tds
    tds() = tr.tds();

    // make a list of all vertices that belong to the original
    // domain and initialize the basic structure of
    // virtual_vertices_reverse
    std::list<Vertex_handle> vlist;
    for(Vertex_iterator vit = vertices_begin(); vit != vertices_end(); ++vit)
    {
      if(vit->offset() == Offset())
      {
        vlist.push_back(vit);
        virtual_vertices_reverse.insert(std::make_pair(vit,std::vector<Vertex_handle>(26)));
        CGAL_triangulation_assertion(virtual_vertices_reverse.find(vit)->second.size() == 26);
      }
    }

    // Iterate over all vertices that are not in the original domain
    // and construct the respective entries to virtual_vertices and
    // virtual_vertices_reverse
    for(Vertex_iterator vit2 = vertices_begin(); vit2 != vertices_end(); ++vit2)
    {
      if(vit2->offset() != Offset())
      {
        typename std::list<Vertex_handle>::iterator vlist_it
          = std::find_if(vlist.begin(), vlist.end(), Finder(this,vit2->point()));
        Offset off = vit2->offset();
        virtual_vertices.insert(std::make_pair(vit2, std::make_pair(*vlist_it,off)));
        virtual_vertices_reverse.find(*vlist_it)->second[9*off[0]+3*off[1]+off[2]-1] = vit2;
        CGAL_triangulation_assertion(get_offset(vit2) == off);
      }
    }

    // Cleanup vertex offsets
    for(Vertex_iterator vit = vertices_begin(); vit != vertices_end(); ++vit)
      vit->clear_offset();
    for(Vertex_iterator vit = tr.vertices_begin(); vit != tr.vertices_end(); ++vit)
      vit->clear_offset();
  }

  /** @name Assignment */
  Periodic_3_triangulation_3& operator=(Periodic_3_triangulation_3 tr)
  {
    swap(tr);
    return *this;
  }

  void swap(Periodic_3_triangulation_3& tr)
  {
    std::swap(tr._gt, _gt);
    _tds.swap(tr._tds);
    std::swap(virtual_vertices,tr.virtual_vertices);
    std::swap(virtual_vertices_reverse,tr.virtual_vertices_reverse);
    std::swap(_cover, tr._cover);
  }

  virtual void clear_covering_data() { }

  /// Clears the triangulation and initializes it again.
  void clear()
  {
    _tds.clear();
    clear_covering_data();
    virtual_vertices.clear();
    virtual_vertices_reverse.clear();
    v_offsets.clear();

    _cover = CGAL::make_array(3,3,3);
    init_tds();
  }

private:
  /// Initializes the triangulation data structure
  void init_tds()
  {
    _tds.set_dimension(-2);
    v_offsets.reserve(48);
  }

public:
  /** @name Access functions */
  const Geometric_traits& geom_traits() const { return _gt; }
  const TDS& tds() const { return _tds; }
  TDS& tds() { return _tds; }

  bool is_parallel() const { return false; }

  virtual void gather_cell_hidden_points(const Cell_handle /*cit*/, std::vector<Point>& /* hidden_points*/) { }
  virtual void reinsert_hidden_points_after_converting_to_1_sheeted(const std::vector<Point>& /* hidden_points*/) { }

  const Iso_cuboid& domain() const { return _gt.get_domain(); }

  virtual void update_cover_data_after_setting_domain() {}

  void set_domain(const Iso_cuboid& domain)
  {
    clear();
    _gt.set_domain(domain);
    update_cover_data_after_setting_domain();
  }

  const Covering_sheets& number_of_sheets() const { return _cover; }
  const std::pair<Vertex_handle, Offset> original_vertex(const Vertex_handle v) const
  {
    return (virtual_vertices.find(v) == virtual_vertices.end()) ?
      std::make_pair(v,Offset()) : virtual_vertices.find(v)->second;
  }
  const std::vector<Vertex_handle>& periodic_copies(const Vertex_handle v) const
  {
    CGAL_triangulation_precondition(number_of_sheets() != CGAL::make_array(1,1,1));
    CGAL_triangulation_precondition(virtual_vertices.find(v) == virtual_vertices.end());
    CGAL_triangulation_assertion(
        virtual_vertices_reverse.find(v) != virtual_vertices_reverse.end());
    return virtual_vertices_reverse.find(v)->second;
  }

  bool is_triangulation_in_1_sheet() const;

  void convert_to_1_sheeted_covering();
  virtual void update_cover_data_after_converting_to_27_sheeted_covering() { }
  void convert_to_27_sheeted_covering();

  size_type number_of_cells() const {
    if(is_1_cover()) return _tds.number_of_cells();
    else return _tds.number_of_cells()/27;
  }
  size_type number_of_facets() const {
    if(is_1_cover()) return _tds.number_of_facets();
    else return _tds.number_of_facets()/27;
  }
  size_type number_of_edges() const {
    if(is_1_cover()) return _tds.number_of_edges();
    else return _tds.number_of_edges()/27;
  }
  size_type number_of_vertices() const {
    if(is_1_cover()) return _tds.number_of_vertices();
    else return _tds.number_of_vertices()/27;
  }

  size_type number_of_stored_cells() const {
    return _tds.number_of_cells();
  }
  size_type number_of_stored_facets() const {
    return _tds.number_of_facets();
  }
  size_type number_of_stored_edges() const {
    return _tds.number_of_edges();
  }
  size_type number_of_stored_vertices() const {
    return _tds.number_of_vertices();
  }

public:
  bool is_1_cover() const
  {
    bool flag;
    flag = ((_cover[0] == 1) && (_cover[1] == 1) && (_cover[2] == 1));
    return flag;
  }

  void set_cover(const Covering_sheets& cover)
  {
    _cover = cover;
  }

public:
  bool is_virtual(Vertex_handle v)
  {
    if(is_1_cover())
      return false;
    return (virtual_vertices.find(v) != virtual_vertices.end());
  }

public:
  // Offset converters
  int off_to_int(const Offset& off) const
  {
    CGAL_triangulation_assertion( off.x() == 0 || off.x() == 1 );
    CGAL_triangulation_assertion( off.y() == 0 || off.y() == 1 );
    CGAL_triangulation_assertion( off.z() == 0 || off.z() == 1 );
    int i = ((off.x()&1)<<2) + ((off.y()&1)<<1) + ((off.z()&1));
    return i;
  }

  Offset int_to_off(int i) const
  {
    return Offset((i>>2)&1,(i>>1)&1,i&1);
  }

  void set_offsets(Cell_handle c, int o0,int o1,int o2,int o3)
  {
    int off0[3] = {(o0>>2)&1,(o0>>1)&1,(o0&1)};
    int off1[3] = {(o1>>2)&1,(o1>>1)&1,(o1&1)};
    int off2[3] = {(o2>>2)&1,(o2>>1)&1,(o2&1)};
    int off3[3] = {(o3>>2)&1,(o3>>1)&1,(o3&1)};
    for(int i=0; i<3; i++) {
      int min_off = (std::min)((std::min)(off0[i],off1[i]),
                               (std::min)(off2[i],off3[i]));
      if(min_off != 0) {
        off0[i] -= min_off; off1[i] -= min_off;
        off2[i] -= min_off; off3[i] -= min_off;
      }
    }
    o0 = ((off0[0]&1)<<2)+((off0[1]&1)<<1)+(off0[2]&1);
    o1 = ((off1[0]&1)<<2)+((off1[1]&1)<<1)+(off1[2]&1);
    o2 = ((off2[0]&1)<<2)+((off2[1]&1)<<1)+(off2[2]&1);
    o3 = ((off3[0]&1)<<2)+((off3[1]&1)<<1)+(off3[2]&1);
    c->set_offsets(o0,o1,o2,o3);
  }

  template <class Offset>
  void set_offsets(Cell_handle c, Offset o0,Offset o1,Offset o2,Offset o3)
  {
    int off0[3] = {o0.x(),o0.y(),o0.z()};
    int off1[3] = {o1.x(),o1.y(),o1.z()};
    int off2[3] = {o2.x(),o2.y(),o2.z()};
    int off3[3] = {o3.x(),o3.y(),o3.z()};
    for(int i=0; i<3; i++) {
      int min_off = (std::min)((std::min)(off0[i],off1[i]),
                               (std::min)(off2[i],off3[i]));
      if(min_off != 0) {
        off0[i] -= min_off; off1[i] -= min_off;
        off2[i] -= min_off; off3[i] -= min_off;
      }
    }

    CGAL_triangulation_assertion((std::min)((std::min)(off0[0],off1[0]),
                                            (std::min)(off2[0],off3[0])) == 0);
    CGAL_triangulation_assertion((std::min)((std::min)(off0[1],off1[1]),
                                            (std::min)(off2[1],off3[1])) == 0);
    CGAL_triangulation_assertion((std::min)((std::min)(off0[2],off1[2]),
                                            (std::min)(off2[2],off3[2])) == 0);
    CGAL_triangulation_assertion((0 <= off0[0]) && (off0[0] < 2));
    CGAL_triangulation_assertion((0 <= off1[0]) && (off1[0] < 2));
    CGAL_triangulation_assertion((0 <= off2[0]) && (off2[0] < 2));
    CGAL_triangulation_assertion((0 <= off3[0]) && (off3[0] < 2));
    CGAL_triangulation_assertion((0 <= off0[1]) && (off0[1] < 2));
    CGAL_triangulation_assertion((0 <= off1[1]) && (off1[1] < 2));
    CGAL_triangulation_assertion((0 <= off2[1]) && (off2[1] < 2));
    CGAL_triangulation_assertion((0 <= off3[1]) && (off3[1] < 2));
    CGAL_triangulation_assertion((0 <= off0[2]) && (off0[2] < 2));
    CGAL_triangulation_assertion((0 <= off1[2]) && (off1[2] < 2));
    CGAL_triangulation_assertion((0 <= off2[2]) && (off2[2] < 2));
    CGAL_triangulation_assertion((0 <= off3[2]) && (off3[2] < 2));

    int o0i = ((off0[0]&1)<<2)+((off0[1]&1)<<1)+(off0[2]&1);
    int o1i = ((off1[0]&1)<<2)+((off1[1]&1)<<1)+(off1[2]&1);
    int o2i = ((off2[0]&1)<<2)+((off2[1]&1)<<1)+(off2[2]&1);
    int o3i = ((off3[0]&1)<<2)+((off3[1]&1)<<1)+(off3[2]&1);
    c->set_offsets(o0i,o1i,o2i,o3i);
  }

public:
  /** @name Wrapping the traits */
  // Note that calling functors with "construct_point(p), offset" and not
  // construct_point(p, offset) is done on purpose. Indeed, construct_point(p)
  // is not a real construction: it's either the identity (when `Point` is `Point_3`)
  // or getting the bare point within a weighted point. However, construct_point(p, offset)
  // is a real construction. When using filters, construct_point() must be done
  // within the filtered predicates with the appropriate construct_point_3_object.

  template<typename P> // can be Point or Point_3
  Comparison_result compare_xyz(const P& p1, const P& p2) const {
    return geom_traits().compare_xyz_3_object()(construct_point(p1),
                                                construct_point(p2));
  }
  template<typename P> // can be Point or Point_3
  Comparison_result compare_xyz(const P& p1, const P& p2,
                                const Offset& o1, const Offset& o2) const {
    return geom_traits().compare_xyz_3_object()(
             construct_point(p1), construct_point(p2), o1, o2);
  }
  template<typename P> // can be Point or Point_3
  Comparison_result compare_xyz(const std::pair<P, Offset>& pp1,
                                const std::pair<P, Offset>& pp2) const {
    return geom_traits().compare_xyz_3_object()(
             construct_point(pp1.first), construct_point(pp2.first), pp1.second, pp2.second);
  }

  template<typename P> // can be Point or Point_3
  Orientation orientation(const P& p1, const P& p2, const P& p3, const P& p4) const {
    return geom_traits().orientation_3_object()(construct_point(p1), construct_point(p2),
                                                construct_point(p3), construct_point(p4));
  }
  template<typename P> // can be Point or Point_3
  Orientation orientation(const P& p1, const P& p2, const P& p3, const P& p4,
                          const Offset& o1, const Offset& o2,
                          const Offset& o3, const Offset& o4) const {
    return geom_traits().orientation_3_object()(
             construct_point(p1), construct_point(p2), construct_point(p3), construct_point(p4),
             o1, o2, o3, o4);
  }
  template<typename P>  // can be Point or Point_3
  Orientation orientation(const std::pair<P, Offset>& pp1, const std::pair<P, Offset>& pp2,
                          const std::pair<P, Offset>& pp3, const std::pair<P, Offset>& pp4) const {
    return geom_traits().orientation_3_object()(
             construct_point(pp1.first), construct_point(pp2.first),
             construct_point(pp3.first), construct_point(pp4.first),
             pp1.second, pp2.second, pp3.second, pp4.second);
  }

  template<typename P>  // can be Point or Point_3
  bool equal(const P& p1, const P& p2) const {
    return compare_xyz(construct_point(p1), construct_point(p2)) == EQUAL;
  }
  template<typename P> // can be Point or Point_3
  bool equal(const P& p1, const P& p2, const Offset& o1, const Offset& o2) const {
    return compare_xyz(construct_point(p1), construct_point(p2), o1, o2) == EQUAL;
  }
  template<typename P>  // can be Point or Point_3
  bool equal(const std::pair<P, Offset>& pp1, const std::pair<P, Offset>& pp2) const {
    return compare_xyz(construct_point(pp1.first), construct_point(pp1.second),
                       pp2.first, pp2.second) == EQUAL;
  }

  template<typename P>  // can be Point or Point_3
  bool coplanar(const P& p1, const P& p2, const P& p3, const P& p4) const {
    return orientation(construct_point(p1), construct_point(p2),
                       construct_point(p3), construct_point(p4)) == COPLANAR;
  }
  template<typename P> // can be Point or Point_3
  bool coplanar(const P& p1, const P& p2, const P& p3, const P& p4,
                const Offset& o1, const Offset& o2, const Offset& o3, const Offset& o4) const {
    return orientation(construct_point(p1), construct_point(p2),
                       construct_point(p3), construct_point(p4),
                       o1, o2, o3, o4) == COPLANAR;
  }
  template<typename P>  // can be Point or Point_3
  bool coplanar(const std::pair<P, Offset>& pp1, const std::pair<P, Offset>& pp2,
                const std::pair<P, Offset>& pp3, const std::pair<P, Offset>& pp4) const {
    return orientation(construct_point(pp1.first), construct_point(pp2.first),
                       construct_point(pp3.first), construct_point(pp4.first),
                       pp1.second, pp2.second, pp3.second, pp4.second) == COPLANAR;
  }

public:
  /** @name Geometric access functions */

  // *************************************************************************
  // -*-*-*-*-*-*-*-*-*-*-*-*-*-* POINT -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
  // *************************************************************************

  // Note that construct_point()-like functions have return type Point_3,
  // but       point()          -like functions have return type Point

  template<typename P> // can be Point or Point_3
  typename boost::result_of<const typename GT::Construct_point_3(const P&)>::type
  construct_point(const P& p) const {
    return geom_traits().construct_point_3_object()(p);
  }

  template<typename P> // can be Point or Point_3
  Point_3 construct_point(const P& p, const Offset& o) const {
    return geom_traits().construct_point_3_object()(p, o);
  }
  template<typename P> // can be Point or Point_3
  Point_3 construct_point(const std::pair<P, Offset>& pp) const {
    return construct_point(pp.first, pp.second);
  }

  Periodic_point_3 construct_periodic_point(const Point_3& p,
                                            bool& had_to_use_exact) const
  {
    // The function is a different file to be able to be used where there is
    // no triangulation (namely, the domains of Periodic_3_mesh_3).
    return ::CGAL::P3T3::internal::construct_periodic_point(p, had_to_use_exact, geom_traits());
  }

  Periodic_point_3 construct_periodic_point(const Point_3& p) const
  {
    bool useless = false;
    return construct_periodic_point(p, useless);
  }

  // ---------------------------------------------------------------------------
  // The following functions return objects of type Point and Periodic_point,
  // _not_ Point_3 and Periodic_point_3.
  // They are templated by `construct_point` to distingush between Delaunay and
  // regular triangulations
  // ---------------------------------------------------------------------------

  template <class ConstructPoint>
  Point point(const Periodic_point& pp, ConstructPoint cp) const
  {
    return cp(pp.first, pp.second);
  }

  // The following functions return the "real" position in space (unconstrained
  // to the fundamental domain) of the vertices v and c->vertex(idx),
  // respectively
  template <class ConstructPoint>
  Point point(Vertex_handle v, ConstructPoint cp) const {
    return point(periodic_point(v), cp);
  }

  template <class ConstructPoint>
  Point point(Cell_handle c, int idx, ConstructPoint cp) const
  {
//    if(is_1_cover())
    return point(periodic_point(c, idx), cp);

    Offset vec_off[4];
    for(int i=0; i<4; i++)
      vec_off[i] = periodic_point(c,i).second;

    int ox = vec_off[0].x();
    int oy = vec_off[0].y();
    int oz = vec_off[0].z();
    for(int i=1; i<4; i++) {
      ox = (std::min)(ox, vec_off[i].x());
      oy = (std::min)(oy, vec_off[i].y());
      oz = (std::min)(oz, vec_off[i].z());
    }
    Offset diff_off(-ox, -oy, -oz);
    if(diff_off.is_null())
      return point(periodic_point(c, idx), cp);

    for(unsigned int i=0; i<4; i++)
      vec_off[i] += diff_off;
    Vertex_handle canonic_vh[4];
    for(int i=0; i<4; i++) {
      Virtual_vertex_map_it vvmit = virtual_vertices.find(c->vertex(i));
      Vertex_handle orig_vh;
      if(vvmit == virtual_vertices.end())
        orig_vh = c->vertex(i);
      else
        orig_vh = vvmit->second.first;

      if(vec_off[i].is_null()) {
        canonic_vh[i] = orig_vh;
      }
      else {
        CGAL_assertion(virtual_vertices_reverse.find(orig_vh)
                       != virtual_vertices_reverse.end());
        canonic_vh[i] = virtual_vertices_reverse.find(orig_vh)->
                          second[9*vec_off[i][0]+3*vec_off[i][1]+vec_off[i][2]-1];
      }
    }

    std::vector<Cell_handle> cells;
    incident_cells(canonic_vh[0], std::back_inserter(cells));
    for(unsigned int i=0; i<cells.size(); i++) {
      CGAL_assertion(cells[i]->has_vertex(canonic_vh[0]));
      if(cells[i]->has_vertex(canonic_vh[1])
          && cells[i]->has_vertex(canonic_vh[2])
          && cells[i]->has_vertex(canonic_vh[3]) )
        return point(periodic_point(cells[i], cells[i]->index(canonic_vh[idx])), cp);
    }
    CGAL_assertion(false);
    return Point();
  }

  Periodic_point periodic_point(const Vertex_handle v) const
  {
    if(is_1_cover())
      return std::make_pair(v->point(), Offset(0,0,0));

    Virtual_vertex_map_it it = virtual_vertices.find(v);
    if(it == virtual_vertices.end()) {
      // if v is not contained in virtual_vertices, then it is in the
      // original domain.
      return std::make_pair(v->point(), Offset(0,0,0));
    } else {
      // otherwise it has to be looked up as well as its offset.
      return std::make_pair(it->second.first->point(), it->second.second);
    }
  }

  Periodic_point periodic_point(const Cell_handle c, int i) const
  {
    if(is_1_cover())
      return std::make_pair(c->vertex(i)->point(), int_to_off(c->offset(i)));

    Virtual_vertex_map_it it = virtual_vertices.find(c->vertex(i));
    if(it == virtual_vertices.end()) {
      // if c->vertex(i) is not contained in virtual_vertices, then it
      // is in the original domain.
      return std::make_pair(c->vertex(i)->point(),
                            combine_offsets(Offset(),int_to_off(c->offset(i))) );
    } else {
      // otherwise it has to be looked up as well as its offset.
      return std::make_pair(it->second.first->point(),
                            combine_offsets(it->second.second, int_to_off(c->offset(i))) );
    }
  }

  // *************************************************************************
  // -*-*-*-*-*-*-*-*-*-*-*-*-*-* SEGMENT -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
  // *************************************************************************

  // Note that an object of type `Segment` is constructed from two objects of type `Point_3`
  // while an object of type `Periodic_segment` is constructed from two objects of type `Point`

  template<typename P> // can be Point or Point_3
  Segment construct_segment(const P& p1, const P& p2) const {
    return geom_traits().construct_segment_3_object()(construct_point(p1),
                                                      construct_point(p2));
  }
  template<typename P> // can be Point or Point_3
  Segment construct_segment(const P& p1, const P& p2,
                            const Offset& o1, const Offset& o2) const {
    return geom_traits().construct_segment_3_object()(construct_point(p1, o1),
                                                      construct_point(p2, o2));
  }
  template<typename PS> // can be Periodic_segment or Periodic_segment_3
  Segment construct_segment(const PS& ps) const {
    return geom_traits().construct_segment_3_object()(
             construct_point(ps[0].first, ps[0].second),
             construct_point(ps[1].first, ps[1].second));
  }
  template<typename P> // can be Point or Point_3
  Periodic_segment construct_periodic_segment(const P& p1, const P& p2,
                                              const Offset& o1 = Offset(),
                                              const Offset& o2 = Offset()) const {
    return CGAL::make_array(std::make_pair(construct_point(p1), o1),
                            std::make_pair(construct_point(p2), o2));
  }

  Periodic_segment periodic_segment(const Cell_handle c, int i, int j) const {
    CGAL_triangulation_precondition( i != j );
    CGAL_triangulation_precondition( number_of_vertices() != 0 );
    CGAL_triangulation_precondition( i >= 0 && i <= 3 && j >= 0 && j <= 3 );
    return CGAL::make_array(std::make_pair(c->vertex(i)->point(), get_offset(c,i)),
                            std::make_pair(c->vertex(j)->point(), get_offset(c,j)));
  }
  Periodic_segment periodic_segment(const Edge& e) const {
    return periodic_segment(e.first,e.second,e.third);
  }
  Periodic_segment periodic_segment(const Cell_handle c, Offset offset, int i, int j) const {
    Periodic_segment result = periodic_segment(c,i,j);
    offset.x() *= _cover[0];
    offset.y() *= _cover[1];
    offset.z() *= _cover[2];
    result[0].second += offset;
    result[1].second += offset;
    return result;
  }

  // *************************************************************************
  // -*-*-*-*-*-*-*-*-*-*-*-*-*-* TRIANGLE -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
  // *************************************************************************

  // Note that an object of type `Triangle` is constructed from three objects of type `Point_3`
  // while an object of type `Periodic_triangle` is constructed from three objects of type `Point`

  template<typename P> // can be Point or Point_3
  Triangle construct_triangle(const P& p1, const P& p2, const P& p3) const {
    return geom_traits().construct_triangle_3_object()(construct_point(p1),
                                                       construct_point(p2),
                                                       construct_point(p3));
  }
  template<typename P> // can be Point or Point_3
  Triangle construct_triangle(const P& p1, const P& p2, const P& p3,
                              const Offset& o1, const Offset& o2, const Offset& o3) const {
    return geom_traits().construct_triangle_3_object()(construct_point(p1, o1),
                                                       construct_point(p2, o2),
                                                       construct_point(p3, o3));
  }
  template<typename PT> // can be Periodic_triangle or Periodic_triangle_3
  Triangle construct_triangle(const PT& tri) const {
    return geom_traits().construct_triangle_3_object()(
             construct_point(tri[0].first, tri[0].second),
             construct_point(tri[1].first, tri[1].second),
             construct_point(tri[2].first, tri[2].second));
  }
  template<typename P> // can be Point or Point_3
  Periodic_triangle construct_periodic_triangle(const P& p1, const P& p2,
                                                const P& p3,
                                                const Offset& o1 = Offset(),
                                                const Offset& o2 = Offset(),
                                                const Offset& o3 = Offset()) const {
    return CGAL::make_array(std::make_pair(construct_point(p1), o1),
                            std::make_pair(construct_point(p2), o2),
                            std::make_pair(construct_point(p3), o3));
  }

  Periodic_triangle periodic_triangle(const Cell_handle c, int i) const {
    CGAL_triangulation_precondition( number_of_vertices() != 0 );
    CGAL_triangulation_precondition( i >= 0 && i <= 3 );
    if( (i&1)==0 )
      return CGAL::make_array(std::make_pair(c->vertex((i+2)&3)->point(), get_offset(c,(i+2)&3)),
                              std::make_pair(c->vertex((i+1)&3)->point(), get_offset(c,(i+1)&3)),
                              std::make_pair(c->vertex((i+3)&3)->point(), get_offset(c,(i+3)&3)) );

    return CGAL::make_array(std::make_pair(c->vertex((i+1)&3)->point(), get_offset(c,(i+1)&3)),
                            std::make_pair(c->vertex((i+2)&3)->point(), get_offset(c,(i+2)&3)),
                            std::make_pair(c->vertex((i+3)&3)->point(), get_offset(c,(i+3)&3)) );
  }
  Periodic_triangle periodic_triangle(const Facet& f) const {
    return periodic_triangle(f.first, f.second);
  }
  Periodic_triangle periodic_triangle(const Cell_handle c, Offset offset, int i) const {
    Periodic_triangle result = periodic_triangle(c,i);
    offset.x() *= _cover[0];
    offset.y() *= _cover[1];
    offset.z() *= _cover[2];
    result[0].second += offset;
    result[1].second += offset;
    result[2].second += offset;
    return result;
  }

  // *************************************************************************
  // -*-*-*-*-*-*-*-*-*-*-*-*-*-* TETRAHEDRON -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
  // *************************************************************************

  // Note that an object of type `Tetrahedron` is constructed from four objects of type `Point_3`
  // while an object of type `Periodic_tetrahedron` is constructed from four objects of type `Point`

  template<typename P> // can be Point or Point_3
  Tetrahedron construct_tetrahedron(const P& p1, const P& p2,
                                    const P& p3, const P& p4) const {
    return geom_traits().construct_tetrahedron_3_object()(construct_point(p1),
                                                          construct_point(p2),
                                                          construct_point(p3),
                                                          construct_point(p4));
  }
  template<typename P> // can be Point or Point_3
  Tetrahedron construct_tetrahedron(const P& p1, const P& p2,
                                    const P& p3, const P& p4,
                                    const Offset& o1, const Offset& o2,
                                    const Offset& o3, const Offset& o4) const {
    return geom_traits().construct_tetrahedron_3_object()(
             construct_point(p1, o1), construct_point(p2, o2),
             construct_point(p3, o3), construct_point(p4, o4));
  }
  template<typename PT> // can be Periodic_tetrahedron or Periodic_tetrahedron_3
  Tetrahedron construct_tetrahedron(const PT& tet) const {
    return geom_traits().construct_tetrahedron_3_object()(
             construct_point(tet[0].first, tet[0].second),
             construct_point(tet[1].first, tet[1].second),
             construct_point(tet[2].first, tet[2].second),
             construct_point(tet[3].first, tet[3].second));
  }
  template<typename P> // can be Point or Point_3
  Periodic_tetrahedron construct_periodic_tetrahedron(const P& p1, const P& p2,
                                                      const P& p3, const P& p4,
                                                      const Offset& o1 = Offset(),
                                                      const Offset& o2 = Offset(),
                                                      const Offset& o3 = Offset(),
                                                      const Offset& o4 = Offset()) const {
    return CGAL::make_array(std::make_pair(construct_point(p1), o1),
                            std::make_pair(construct_point(p2), o2),
                            std::make_pair(construct_point(p3), o3),
                            std::make_pair(construct_point(p4), o4));
  }

  Periodic_tetrahedron periodic_tetrahedron(const Cell_handle c) const
  {
    CGAL_triangulation_precondition( number_of_vertices() != 0 );
    return CGAL::make_array(std::make_pair(c->vertex(0)->point(), get_offset(c,0)),
                            std::make_pair(c->vertex(1)->point(), get_offset(c,1)),
                            std::make_pair(c->vertex(2)->point(), get_offset(c,2)),
                            std::make_pair(c->vertex(3)->point(), get_offset(c,3)) );
  }
  Periodic_tetrahedron periodic_tetrahedron(const Cell_handle c, Offset offset) const
  {
    Periodic_tetrahedron result = periodic_tetrahedron(c);
    offset.x() *= _cover[0];
    offset.y() *= _cover[1];
    offset.z() *= _cover[2];
    result[0].second += offset;
    result[1].second += offset;
    result[2].second += offset;
    result[3].second += offset;
    return result;
  }

  // end of geometric functions

public:
  /** @name Queries */
  bool is_vertex(const Point& p, Vertex_handle& v) const;

  bool is_vertex(Vertex_handle v) const {
    return _tds.is_vertex(v);
  }

  bool is_edge(Vertex_handle u, Vertex_handle v, Cell_handle& c, int& i, int& j) const {
    return _tds.is_edge(u, v, c, i, j);
  }

  bool is_edge(Vertex_handle u, const Offset& off_u,
               Vertex_handle v, const Offset& off_v,
               Cell_handle& c, int& i, int& j) const
  {
    if(!_tds.is_edge(u,v,c,i,j))
      return false;

    if((get_offset(c,i) == off_u) && (get_offset(c,j) == off_v))
      return true;

    // it might be that different cells containing (u,v) yield
    // different offsets, which forces us to test for all possibilities.
    else {
      Cell_circulator ccirc = incident_cells(c,i,j,c);
      while(++ccirc != c) {
        i = ccirc->index(u);
        j = ccirc->index(v);
        if((get_offset(ccirc,i) == off_u) && (get_offset(ccirc,j) == off_v)) {
          c = ccirc;
          return true;
        }
      }
      return false;
    }
  }

  bool is_facet(Vertex_handle u, Vertex_handle v, Vertex_handle w,
                Cell_handle& c, int& i, int& j, int& k) const {
    return _tds.is_facet(u, v, w, c, i, j, k);
  }

  bool is_facet(Vertex_handle u, const Offset& off_u,
                Vertex_handle v, const Offset& off_v,
                Vertex_handle w, const Offset& off_w,
                Cell_handle& c, int& i, int& j, int& k) const
  {
    if(!_tds.is_facet(u,v,w,c,i,j,k)) return false;
    if((get_offset(c,i) == off_u)
        && (get_offset(c,j) == off_v)
        && (get_offset(c,k) == off_w) )
      return true;
    // it might be that c and c->neighbor(l) yield different offsets
    // which forces us to test for both possibilities.
    else {
      int l = 6-i-j-k;
      c = c->neighbor(l);
      i = c->index(u);
      j = c->index(v);
      k = c->index(w);
      return ((get_offset(c,i) == off_u)
              && (get_offset(c,j) == off_v)
              && (get_offset(c,k) == off_w) );
    }
  }

  bool is_cell(Cell_handle c) const {
    return _tds.is_cell(c);
  }

  bool is_cell(Vertex_handle u, Vertex_handle v,
      Vertex_handle w, Vertex_handle t,
      Cell_handle& c, int& i, int& j, int& k, int& l) const {
    return _tds.is_cell(u, v, w, t, c, i, j, k, l);
  }

  bool is_cell(Vertex_handle u, Vertex_handle v, Vertex_handle w,
      Vertex_handle t, Cell_handle& c) const
  {
    int i,j,k,l;
    return _tds.is_cell(u, v, w, t, c, i, j, k, l);
  }

  bool is_cell(Vertex_handle u, const Offset& off_u,
               Vertex_handle v, const Offset& off_v,
               Vertex_handle w, const Offset& off_w,
               Vertex_handle t, const Offset& off_t,
               Cell_handle& c, int& i, int& j, int& k, int& l) const
  {
    if(!_tds.is_cell(u,v,w,t,c,i,j,k,l)) return false;
    return ((get_offset(c,i) == off_u)
            && (get_offset(c,j) == off_v)
            && (get_offset(c,k) == off_w)
            && (get_offset(c,l) == off_t) );
    return false;
  }

  bool is_cell(Vertex_handle u, const Offset& off_u,
               Vertex_handle v, const Offset& off_v,
               Vertex_handle w, const Offset& off_w,
               Vertex_handle t, const Offset& off_t,
               Cell_handle& c) const
  {
    int i, j, k, l;
    return is_cell(u,off_u,v,off_v,w,off_w,t,off_t,c,i,j,k,l);
  }

  bool has_vertex(const Facet& f, Vertex_handle v, int& j) const {
    return _tds.has_vertex(f.first, f.second, v, j);
  }
  bool has_vertex(Cell_handle c, int i, Vertex_handle v, int& j) const {
    return _tds.has_vertex(c, i, v, j);
  }
  bool has_vertex(const Facet& f, Vertex_handle v) const {
    return _tds.has_vertex(f.first, f.second, v);
  }
  bool has_vertex(Cell_handle c, int i, Vertex_handle v) const {
    return _tds.has_vertex(c, i, v);
  }

  bool are_equal(Cell_handle c, int i, Cell_handle n, int j) const {
    return _tds.are_equal(c, i, n, j);
  }
  bool are_equal(const Facet& f, const Facet& g) const {
    return _tds.are_equal(f.first, f.second, g.first, g.second);
  }
  bool are_equal(const Facet& f, Cell_handle n, int j) const {
    return _tds.are_equal(f.first, f.second, n, j);
  }

#ifdef CGAL_NO_STRUCTURAL_FILTERING
  Cell_handle
  periodic_locate(const Point& p, const Offset& o_p, Offset& lo,
                  Locate_type& lt, int& li, int& lj,
                  Cell_handle start = Cell_handle()) const;
#else // no CGAL_NO_STRUCTURAL_FILTERING
#  ifndef CGAL_P3T3_STRUCTURAL_FILTERING_MAX_VISITED_CELLS
#    define CGAL_P3T3_STRUCTURAL_FILTERING_MAX_VISITED_CELLS 2500
#  endif // no CGAL_P3T3_STRUCTURAL_FILTERING_MAX_VISITED_CELLS

public:
  Cell_handle
  inexact_periodic_locate(const Point& p, const Offset& o_p,
                 Cell_handle start = Cell_handle(),
                 int max_num_cells = CGAL_P3T3_STRUCTURAL_FILTERING_MAX_VISITED_CELLS) const;
protected:
  Cell_handle
  exact_periodic_locate(const Point& p, const Offset& o_p, Offset& lo,
                        Locate_type& lt,
                        int& li, int& lj,
                        Cell_handle start) const;

  Cell_handle
  generic_periodic_locate(const Point& p, const Offset& o_p, Offset& lo,
                          Locate_type& lt,
                          int& li, int& lj,
                          Cell_handle start,
                          internal::Periodic_structural_filtering_3_tag) const {
    return exact_periodic_locate(p, o_p, lo, lt, li, lj, inexact_periodic_locate(p, o_p, start));
  }

  Cell_handle
  generic_periodic_locate(const Point& p, const Offset& o_p, Offset& lo,
                          Locate_type& lt,
                          int& li, int& lj,
                          Cell_handle start,
                          internal::No_periodic_structural_filtering_3_tag) const {
    return exact_periodic_locate(p, o_p, lo, lt, li, lj, start);
  }

  Orientation
  inexact_orientation(const Point& p, const Point& q,
                      const Point& r, const Point& s,
                      const Offset& o_p = Offset(), const Offset& o_q = Offset(),
                      const Offset& o_r = Offset(), const Offset& o_s = Offset()) const
  {
    const Point_3 bp = construct_point(p, o_p);
    const Point_3 bq = construct_point(q, o_q);
    const Point_3 br = construct_point(r, o_r);
    const Point_3 bs = construct_point(s, o_s);

    // So that this code works well with Lazy_kernel
    internal::Static_filters_predicates::Get_approx<Point_3> get_approx;
    const double px = to_double(get_approx(bp).x());
    const double py = to_double(get_approx(bp).y());
    const double pz = to_double(get_approx(bp).z());
    const double qx = to_double(get_approx(bq).x());
    const double qy = to_double(get_approx(bq).y());
    const double qz = to_double(get_approx(bq).z());
    const double rx = to_double(get_approx(br).x());
    const double ry = to_double(get_approx(br).y());
    const double rz = to_double(get_approx(br).z());
    const double sx = to_double(get_approx(bs).x());
    const double sy = to_double(get_approx(bs).y());
    const double sz = to_double(get_approx(bs).z());

    const double pqx = qx - px;
    const double pqy = qy - py;
    const double pqz = qz - pz;
    const double prx = rx - px;
    const double pry = ry - py;
    const double prz = rz - pz;
    const double psx = sx - px;
    const double psy = sy - py;
    const double psz = sz - pz;

    const double det = determinant(pqx, pqy, pqz,
                                   prx, pry, prz,
                                   psx, psy, psz);
    if(det > 0) return POSITIVE;
    if(det < 0) return NEGATIVE;
    return ZERO;
  }

public:
  Cell_handle
  periodic_locate(const Point& p, const Offset& o_p, Offset& lo,
                  Locate_type& lt, int& li, int& lj,
                  Cell_handle start = Cell_handle()) const
  {
    typedef Triangulation_structural_filtering_traits<Geometric_traits> TSFT;
    typedef typename internal::Periodic_structural_filtering_selector_3<
      TSFT::Use_structural_filtering_tag::value >::Tag Should_filter_tag;

    return generic_periodic_locate(p, o_p, lo, lt, li, lj, start, Should_filter_tag());
  }

  Cell_handle
  inexact_locate(const Point& p,
                 Cell_handle start = Cell_handle(),
                 int max_num_cells = CGAL_P3T3_STRUCTURAL_FILTERING_MAX_VISITED_CELLS) const
  {
    return inexact_periodic_locate(p, Offset(), start, max_num_cells);
  }
#endif // no CGAL_NO_STRUCTURAL_FILTERING

protected:
  /** @name Location helpers */
//  Cell_handle periodic_locate(const Point& p, const Offset& o_p,
//    Locate_type& lt, int& li, int& lj, Cell_handle start) const;

  Bounded_side side_of_cell(const Point& p, const Offset& off, Cell_handle c,
                            Locate_type& lt, int& i, int& j) const;

public:
  /** @name Point Location */
  /** Wrapper function for locate if only the request point is given.
    */
  Cell_handle locate(const Point& p, Cell_handle start = Cell_handle()) const
  {
    Locate_type lt;
    int li, lj;
    return locate( p, lt, li, lj, start);
  }

  /** Wrapper function calling locate with an empty offset if there was no
    * offset given.
    */
  Cell_handle locate(const Point& p,
                     Locate_type& lt, int& li, int& lj,
                     Cell_handle start = Cell_handle()) const
  {
    Offset lo;
    return locate( p, lo, lt, li, lj, start);
  }

  Cell_handle locate(const Point& p, Offset& lo,
                     Cell_handle start = Cell_handle()) const
  {
    Locate_type lt;
    int li, lj;
    return locate( p, lo, lt, li, lj, start);
  }

  Cell_handle locate(const Point& p, Offset& lo,
                     Locate_type& lt, int& li, int& lj,
                     Cell_handle start = Cell_handle()) const
  {
    Cell_handle ch = periodic_locate(p, Offset(), lo, lt, li, lj, start);
    for(unsigned i = 0; i < 3; ++i)
      if(lo[i] >= 1)
        lo[i] = -1;
    return ch;
  }

  Bounded_side side_of_cell(const Point& p, Cell_handle c,
                            Locate_type& lt, int& i, int& j) const
  {
    if(number_of_vertices() == 0) {
      lt = EMPTY;
      return ON_UNBOUNDED_SIDE;
    }
    return side_of_cell(p,Offset(),c,lt,i,j);
  }

private:
  /** @name Insertion helpers */
  template < class Conflict_tester, class Point_hider, class CoverManager >
  Vertex_handle periodic_insert(const Point& p, const Offset& o, Locate_type lt,
                                Cell_handle c, const Conflict_tester& tester,
                                Point_hider& hider, CoverManager& cover_manager,
                                Vertex_handle vh = Vertex_handle());

  template <class Point_iterator, class Offset_iterator>
  void periodic_sort(Point_iterator /*p_begin*/, Point_iterator /*p_end*/,
                     Offset_iterator /*o_begin*/, Offset_iterator /*o_end*/) const {
    std::cout << "Periodic_sort not yet implemented" << std::endl;
  }

  Vertex_handle create_initial_triangulation(const Point& p);

public:
  std::vector<Vertex_handle> insert_dummy_points();

protected:
  // this is needed for compatibility reasons
  template <class Conflict_test, class OutputIteratorBoundaryFacets,
            class OutputIteratorCells, class OutputIteratorInternalFacets>
  Triple<OutputIteratorBoundaryFacets, OutputIteratorCells,
         OutputIteratorInternalFacets>
  find_conflicts(Cell_handle c,
                 const Conflict_test& tester,
                 Triple<OutputIteratorBoundaryFacets, OutputIteratorCells,
                 OutputIteratorInternalFacets> it) const
  {
    bool b = false;
    Offset off = get_location_offset(tester, c, b);
    if(b)
      return find_conflicts(c,off,tester,it);
    return it;
  }

  template <class Conflict_test, class OutputIteratorBoundaryFacets,
            class OutputIteratorCells, class OutputIteratorInternalFacets>
  Triple<OutputIteratorBoundaryFacets, OutputIteratorCells,
         OutputIteratorInternalFacets>
  find_conflicts(Cell_handle c, const Offset& current_off,
                 const Conflict_test& tester,
                 Triple<OutputIteratorBoundaryFacets,
                        OutputIteratorCells,
                        OutputIteratorInternalFacets> it) const;

protected:

  // COMMON INSERTION for DELAUNAY and REGULAR TRIANGULATION
  template < class Conflict_tester, class Point_hider, class CoverManager >
  Vertex_handle insert_in_conflict(const Point& p, Cell_handle start,
                                   const Conflict_tester& tester,
                                   Point_hider& hider,
                                   CoverManager& cover_manager)
  {
    Locate_type lt = Locate_type();
    int li=0, lj=0;
    Offset lo;
    Cell_handle c = periodic_locate(p, Offset(), lo, lt, li, lj, start);
    return insert_in_conflict(p,lt,c,li,lj,tester,hider, cover_manager);
  }

  template < class Conflict_tester, class Point_hider, class CoverManager >
  Vertex_handle insert_in_conflict(const Point& p, Locate_type lt,
                                   Cell_handle c, int li, int lj,
                                   const Conflict_tester& tester,
                                   Point_hider& hider,
                                   CoverManager& cover_manager);

  template < class InputIterator, class Conflict_tester,
             class Point_hider, class CoverManager>
  std::vector<Vertex_handle> insert_in_conflict(
      InputIterator begin, InputIterator end, Cell_handle start,
      Conflict_tester& tester, Point_hider& hider, CoverManager& cover_manager)
  {
    Vertex_handle new_vertex;
    std::vector<Vertex_handle> double_vertices;
    Locate_type lt = Locate_type();
    int li=0, lj=0;
    CGAL_triangulation_assertion_code( Locate_type lta = Locate_type(); )
    CGAL_triangulation_assertion_code( int ia = 0; )
    CGAL_triangulation_assertion_code( int ja = 0; )
    Cell_handle hint;
    while(begin!=end) {
      tester.set_point(*begin);
      Offset lo;
      hint = periodic_locate(*begin, Offset(), lo, lt, li, lj, start);
      CGAL_triangulation_assertion_code( if(number_of_vertices() != 0) { );
      CGAL_triangulation_assertion(side_of_cell(
        *begin,Offset(), hint, lta, ia, ja) != ON_UNBOUNDED_SIDE);
        CGAL_triangulation_assertion(lta == lt);
        CGAL_triangulation_assertion(ia == li);
        CGAL_triangulation_assertion(ja == lj);
        CGAL_triangulation_assertion_code( }
      );

      new_vertex = insert_in_conflict(*begin,lt,hint,li,lj,tester,hider, cover_manager);
      if(lt == VERTEX)
        double_vertices.push_back(new_vertex);
      if(new_vertex != Vertex_handle())
        start = new_vertex->cell();
      begin++;
    }
    return double_vertices;
  }

private:
  /** @name Removal helpers */
  Vertex_triple make_vertex_triple(const Facet& f) const
  {
    Cell_handle ch = f.first;
    int i = f.second;
    return Vertex_triple(ch->vertex(vertex_triple_index(i,0)),
                         ch->vertex(vertex_triple_index(i,1)),
                         ch->vertex(vertex_triple_index(i,2)));
  }

  void make_canonical(Vertex_triple& t) const;

  void make_hole(Vertex_handle v, std::map<Vertex_triple,Facet>& outer_map,
                 std::vector<Cell_handle>& hole);

protected:
  /** @name Removal */
  template < class PointRemover, class CoverManager >
  bool periodic_remove(Vertex_handle v, PointRemover& remover, CoverManager& cover_manager,
                       const bool abort_if_cover_change = false);

  template < class PointRemover, class CT, class CoverManager >
  void remove(Vertex_handle v, PointRemover& remover, CT& ct, CoverManager& cover_manager);

  void delete_vertex(Vertex_handle vertex_handle)
  {
    tds().delete_vertex(vertex_handle);

    if(!is_1_cover())
    {
      typename Virtual_vertex_map::iterator iter = this->virtual_vertices.find(vertex_handle);
      if(iter != this->virtual_vertices.end())
      {
        Vertex_handle vh = iter->second.first;
        this->virtual_vertices.erase(iter);

        typename Virtual_vertex_reverse_map::iterator origin_it = this->virtual_vertices_reverse.find(vh);
        std::vector<Vertex_handle>& copies = origin_it->second;
        typename std::vector<Vertex_handle>::iterator copy_iter = std::find(copies.begin(), copies.end(), vertex_handle);
        CGAL_triangulation_assertion(copy_iter != copies.end());
        copies.erase(copy_iter);
        if(copies.empty())
          virtual_vertices_reverse.erase(origin_it);
      }
      return;
    }

    CGAL_triangulation_assertion(this->virtual_vertices.find(vertex_handle) == this->virtual_vertices.end());
    CGAL_triangulation_assertion(this->virtual_vertices_reverse.find(vertex_handle) == this->virtual_vertices_reverse.end());
  }

public:
  /** @name Traversal */
  Cell_iterator cells_begin() const {
    return _tds.cells_begin();
  }
  Cell_iterator cells_end() const {
    return _tds.cells_end();
  }

  Vertex_iterator vertices_begin() const {
    return _tds.vertices_begin();
  }
  Vertex_iterator vertices_end() const {
    return _tds.vertices_end();
  }

  Edge_iterator edges_begin() const {
    return _tds.edges_begin();
  }
  Edge_iterator edges_end() const {
    return _tds.edges_end();
  }

  Facet_iterator facets_begin() const {
    return _tds.facets_begin();
  }
  Facet_iterator facets_end() const {
    return _tds.facets_end();
  }

  Cell_iterator finite_cells_begin() const {
    return _tds.cells_begin();
  }
  Cell_iterator finite_cells_end() const {
    return _tds.cells_end();
  }

  Vertex_iterator finite_vertices_begin() const {
    return _tds.vertices_begin();
  }
  Vertex_iterator finite_vertices_end() const {
    return _tds.vertices_end();
  }

  Edge_iterator finite_edges_begin() const {
    return _tds.edges_begin();
  }
  Edge_iterator finite_edges_end() const {
    return _tds.edges_end();
  }

  Facet_iterator finite_facets_begin() const {
    return _tds.facets_begin();
  }
  Facet_iterator finite_facets_end() const {
    return _tds.facets_end();
  }

  All_cells_iterator all_cells_begin() const {
    return _tds.cells_begin();
  }
  All_cells_iterator all_cells_end() const {
    return _tds.cells_end();
  }

  All_vertices_iterator all_vertices_begin() const {
    return _tds.vertices_begin();
  }
  All_vertices_iterator all_vertices_end() const {
    return _tds.vertices_end();
  }

  All_edges_iterator all_edges_begin() const {
    return _tds.edges_begin();
  }
  All_edges_iterator all_edges_end() const {
    return _tds.edges_end();
  }

  All_facets_iterator all_facets_begin() const {
    return _tds.facets_begin();
  }
  All_facets_iterator all_facets_end() const {
    return _tds.facets_end();
  }

  Unique_vertex_iterator unique_vertices_begin() const {
    return CGAL::filter_iterator(vertices_end(), Domain_tester<Self>(this),
                                 vertices_begin());
  }
  Unique_vertex_iterator unique_vertices_end() const {
    return CGAL::filter_iterator(vertices_end(), Domain_tester<Self>(this));
  }

  // Geometric iterators
  Periodic_tetrahedron_iterator periodic_tetrahedra_begin(
      Iterator_type it = STORED) const {
    return Periodic_tetrahedron_iterator(this, it);
  }
  Periodic_tetrahedron_iterator periodic_tetrahedra_end(
      Iterator_type it = STORED) const {
    return Periodic_tetrahedron_iterator(this, 1, it);
  }

  Periodic_triangle_iterator periodic_triangles_begin(
      Iterator_type it = STORED) const {
    return Periodic_triangle_iterator(this, it);
  }
  Periodic_triangle_iterator periodic_triangles_end(
      Iterator_type it = STORED) const {
    return Periodic_triangle_iterator(this, 1, it);
  }

  Periodic_segment_iterator periodic_segments_begin(
      Iterator_type it = STORED) const {
    return Periodic_segment_iterator(this, it);
  }
  Periodic_segment_iterator periodic_segments_end(
      Iterator_type it = STORED) const {
    return Periodic_segment_iterator(this, 1, it);
  }

  Periodic_point_iterator periodic_points_begin(
      Iterator_type it = STORED) const {
    return Periodic_point_iterator(this, it);
  }
  Periodic_point_iterator periodic_points_end(
      Iterator_type it = STORED) const  {
    return Periodic_point_iterator(this, 1, it);
  }

  // Circulators
  Cell_circulator incident_cells(const Edge& e) const {
    return _tds.incident_cells(e);
  }
  Cell_circulator incident_cells(Cell_handle c, int i, int j) const {
    return _tds.incident_cells(c, i, j);
  }
  Cell_circulator incident_cells(const Edge& e, Cell_handle start) const {
    return _tds.incident_cells(e, start);
  }
  Cell_circulator incident_cells(Cell_handle c, int i, int j,
      Cell_handle start) const {
    return _tds.incident_cells(c, i, j, start);
  }

  Facet_circulator incident_facets(const Edge& e) const {
    return _tds.incident_facets(e);
  }
  Facet_circulator incident_facets(Cell_handle c, int i, int j) const {
    return _tds.incident_facets(c, i, j);
  }
  Facet_circulator incident_facets(const Edge& e, const Facet& start) const {
    return _tds.incident_facets(e, start);
  }
  Facet_circulator incident_facets(Cell_handle c, int i, int j,
      const Facet& start) const {
    return _tds.incident_facets(c, i, j, start);
  }
  Facet_circulator incident_facets(const Edge& e,
      Cell_handle start, int f) const {
    return _tds.incident_facets(e, start, f);
  }
  Facet_circulator incident_facets(Cell_handle c, int i, int j,
      Cell_handle start, int f) const {
    return _tds.incident_facets(c, i, j, start, f);
  }

  // around a vertex
  template <class OutputIterator>
  OutputIterator incident_cells(Vertex_handle v, OutputIterator cells) const {
    return _tds.incident_cells(v, cells);
  }

  template <class OutputIterator>
  OutputIterator incident_facets(Vertex_handle v, OutputIterator facets) const {
    return _tds.incident_facets(v, facets);
  }

  template <class OutputIterator>
  OutputIterator incident_edges(Vertex_handle v, OutputIterator edges) const {
    return _tds.incident_edges(v, edges);
  }

  template <class OutputIterator>
  OutputIterator adjacent_vertices(Vertex_handle v, OutputIterator vertices) const {
    return _tds.adjacent_vertices(v, vertices);
  }

  //deprecated, don't use anymore
  template <class OutputIterator>
  OutputIterator incident_vertices(Vertex_handle v, OutputIterator vertices) const {
    return _tds.adjacent_vertices(v, vertices);
  }

  size_type degree(Vertex_handle v) const {
    return _tds.degree(v);
  }

  // Functions forwarded from TDS.
  int mirror_index(Cell_handle c, int i) const {
    return _tds.mirror_index(c, i);
  }

  Vertex_handle mirror_vertex(Cell_handle c, int i) const {
    return _tds.mirror_vertex(c, i);
  }

  Facet mirror_facet(Facet f) const {
    return _tds.mirror_facet(f);
  }

private:
  /** @name Checking helpers */
  /// calls has_self_edges for every cell of the triangulation
  bool has_self_edges() const
  {
    Cell_iterator it;
    for( it = all_cells_begin(); it != all_cells_end(); ++it )
      if(has_self_edges(it)) return true;
    return false;
  }

  bool has_self_edges(Cell_handle c) const;

public:
  /** @name Checking */
  bool is_valid(bool verbose = false, int level = 0) const;
  bool is_valid(Cell_handle c, bool verbose = false, int level = 0) const;

protected:
  template <class ConflictTester>
  bool is_valid_conflict(ConflictTester& tester, bool verbose = false,
                         int level = 0) const;

public:
  /** @name Functors */
  class Perturbation_order
  {
    const Self *t;

  public:
    Perturbation_order(const Self *tr)
      : t(tr) {}

    template<typename P>
    bool operator()(const P* p, const P* q) const {
      return t->compare_xyz(*p, *q) == SMALLER;
    }
  };

public:
  // undocumented access functions
  Offset get_offset(Cell_handle ch, int i) const
  {
    if(is_1_cover())
      return int_to_off(ch->offset(i));

    Virtual_vertex_map_it it = virtual_vertices.find(ch->vertex(i));
    if(it != virtual_vertices.end())
      return combine_offsets(it->second.second, int_to_off(ch->offset(i)));
    else
      return combine_offsets(Offset(), int_to_off(ch->offset(i)));
  }

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

  Offset combine_offsets(const Offset& o_c, const Offset& o_t) const
  {
    Offset o_ct(_cover[0]*o_t.x(), _cover[1]*o_t.y(), _cover[2]*o_t.z());
    return o_c + o_ct;
  }

  // These functions give the pair (vertex, offset) that corresponds to the
  // i-th vertex of cell ch or vertex vh, respectively.
  void get_vertex(Cell_handle ch, int i, Vertex_handle& vh, Offset& off) const;
  void get_vertex(Vertex_handle vh_i, Vertex_handle& vh, Offset& off) const;

protected:
  // Auxiliary functions
  Cell_handle get_cell(const Vertex_handle* vh) const;

  template<class Conflict_tester>
  Offset get_location_offset(const Conflict_tester& tester,
                             Cell_handle c) const;

  template<class Conflict_tester>
  Offset get_location_offset(const Conflict_tester& tester,
                             Cell_handle c, bool& found) const;

  Offset neighbor_offset(Cell_handle ch, int i, Cell_handle nb) const;

public:
  Offset neighbor_offset(Cell_handle ch, int i) const
  {
    return neighbor_offset(ch, i, ch->neighbor(i));
  }

protected:
  /** @name Friends */
  friend std::istream& operator>> <>
      (std::istream& is, Periodic_3_triangulation_3<GT,TDS>& tr);
  friend std::ostream& operator<< <>
      (std::ostream& os, const Periodic_3_triangulation_3<GT,TDS>& tr);

protected:
  template <class ConstructCircumcenter>
  Periodic_point_3 periodic_circumcenter(Cell_handle c,
                                         ConstructCircumcenter construct_circumcenter) const
  {
    CGAL_triangulation_precondition(c != Cell_handle());

    Point_3 p = construct_circumcenter(c->vertex(0)->point(), c->vertex(1)->point(),
                                       c->vertex(2)->point(), c->vertex(3)->point(),
                                       get_offset(c, 0), get_offset(c, 1),
                                       get_offset(c, 2), get_offset(c, 3));

    return construct_periodic_point(p);
  }

public:
  bool is_canonical(const Facet& f) const
  {
    if(number_of_sheets() == CGAL::make_array(1,1,1))
      return true;

    Offset cell_off0 = int_to_off(f.first->offset((f.second+1)&3));
    Offset cell_off1 = int_to_off(f.first->offset((f.second+2)&3));
    Offset cell_off2 = int_to_off(f.first->offset((f.second+3)&3));
    Offset diff_off((cell_off0.x() == 1
                     && cell_off1.x() == 1
                     && cell_off2.x() == 1) ? -1 : 0,
                    (cell_off0.y() == 1
                     && cell_off1.y() == 1
                     && cell_off2.y() == 1) ? -1 : 0,
                    (cell_off0.z() == 1
                     && cell_off1.z() == 1
                     && cell_off2.z() == 1) ? -1 : 0);
    Offset off0 = combine_offsets(get_offset(f.first, (f.second+1)&3), diff_off);
    Offset off1 = combine_offsets(get_offset(f.first, (f.second+2)&3), diff_off);
    Offset off2 = combine_offsets(get_offset(f.first, (f.second+3)&3), diff_off);

    // If there is one offset with entries larger than 1 then we are
    // talking about a vertex that is too far away from the original
    // domain to belong to a canonical triangle.
    if(off0.x() > 1) return false;
    if(off0.y() > 1) return false;
    if(off0.z() > 1) return false;
    if(off1.x() > 1) return false;
    if(off1.y() > 1) return false;
    if(off1.z() > 1) return false;
    if(off2.x() > 1) return false;
    if(off2.y() > 1) return false;
    if(off2.z() > 1) return false;

    // If there is one direction of space for which all offsets are
    // non-zero then the edge is not canonical because we can
    // take the copy closer towards the origin in that direction.
    int offx = off0.x() & off1.x() & off2.x();
    int offy = off0.y() & off1.y() & off2.y();
    int offz = off0.z() & off1.z() & off2.z();

    return (offx == 0 && offy == 0 && offz == 0);
  }

protected:
  template <class ConstructCircumcenter>
  bool canonical_dual_segment(Cell_handle c, int i, Periodic_segment_3& ps,
                              ConstructCircumcenter construct_circumcenter) const
  {
    CGAL_triangulation_precondition(c != Cell_handle());
    Offset off = neighbor_offset(c,i,c->neighbor(i));
    Periodic_point_3 p1 = periodic_circumcenter(c, construct_circumcenter);
    Periodic_point_3 p2 = periodic_circumcenter(c->neighbor(i), construct_circumcenter);
    Offset o1 = -p1.second;
    Offset o2 = combine_offsets(-p2.second,-off);
    Offset cumm_off((std::min)(o1.x(),o2.x()),
                    (std::min)(o1.y(),o2.y()),
                    (std::min)(o1.z(),o2.z()));
    const Periodic_point_3 pp1 = std::make_pair(construct_point(p1), o1-cumm_off);
    const Periodic_point_3 pp2 = std::make_pair(construct_point(p2), o2-cumm_off);
    ps = CGAL::make_array(pp1, pp2);
    return (cumm_off == Offset(0,0,0));
  }

  template <class OutputIterator, class ConstructCircumcenter>
  OutputIterator dual(Cell_handle c, int i, int j,
                      OutputIterator points,
                      ConstructCircumcenter construct_circumcenter) const
  {
    Cell_circulator cstart = incident_cells(c, i, j);

    Offset offv = periodic_point(c,i).second;
    Vertex_handle v = c->vertex(i);

    Cell_circulator ccit = cstart;
    do {
      Point_3 dual_orig = periodic_circumcenter(ccit, construct_circumcenter).first;
      int idx = ccit->index(v);
      Offset off = periodic_point(ccit,idx).second;
      Point_3 dual = construct_point(std::make_pair(dual_orig, -off+offv));
      *points++ = dual;
      ++ccit;
    } while(ccit != cstart);

    return points;
  }

  template <class OutputIterator, class ConstructCircumcenter>
  OutputIterator dual(Vertex_handle v, OutputIterator points,
                      ConstructCircumcenter construct_circumcenter) const
  {
    std::vector<Cell_handle> cells;
    incident_cells(v,std::back_inserter(cells));

    for(unsigned int i=0; i<cells.size(); i++) {
      Point_3 dual_orig = periodic_circumcenter(cells[i], construct_circumcenter).first;
      int idx = cells[i]->index(v);
      Offset off = periodic_point(cells[i],idx).second;
      Point_3 dual = construct_point(std::make_pair(dual_orig, -off));
      *points++ = dual;
    }

    return points;
  }

  template <class Stream, class ConstructCircumcenter>
  Stream& draw_dual(Stream& os, ConstructCircumcenter construct_circumcenter) const
  {
    CGAL_triangulation_assertion_code( unsigned int i = 0; )
    for(Facet_iterator fit = facets_begin(), end = facets_end();
         fit != end; ++fit) {
      if(!is_canonical(*fit))
        continue;

      Periodic_segment_3 pso;
      canonical_dual_segment(fit->first, fit->second, pso, construct_circumcenter);
      Segment so = construct_segment(pso);
      CGAL_triangulation_assertion_code ( ++i; )
      os << so.source() << '\n';
      os << so.target() << '\n';
    }
    CGAL_triangulation_assertion( i == number_of_facets() );

    return os;
  }

  /// Volume computations

  // Note: Polygon area computation requires to evaluate square roots
  // and thus cannot be done without changing the Traits concept.

  template <class ConstructCircumcenter>
  FT dual_volume(Vertex_handle v, ConstructCircumcenter construct_circumcenter) const
  {
    std::list<Edge> edges;
    incident_edges(v, std::back_inserter(edges));

    FT vol(0);
    for(typename std::list<Edge>::iterator eit = edges.begin();
         eit != edges.end(); ++eit) {

      // compute the dual of the edge *eit but handle the translations
      // with respect to the dual of v. That is why we cannot use one
      // of the existing dual functions here.
      Facet_circulator fstart = incident_facets(*eit);
      Facet_circulator fcit = fstart;
      std::vector<Point_3> pts;
      do {
        // TODO: possible speed-up by caching the circumcenters
        Point_3 dual_orig = periodic_circumcenter(fcit->first, construct_circumcenter).first;
        int idx = fcit->first->index(v);
        Offset off = periodic_point(fcit->first,idx).second;
        pts.push_back(construct_point(std::make_pair(dual_orig,-off)));
        ++fcit;
      } while(fcit != fstart);

      Point_3 orig(0,0,0);
      for(unsigned int i=1; i<pts.size()-1; i++)
        vol += Tetrahedron(orig,pts[0],pts[i],pts[i+1]).volume();
    }
    return vol;
  }

  /// Centroid computations

  // Note: Centroid computation for polygons requires to evaluate
  // square roots and thus cannot be done without changing the
  // Traits concept.

  template <class ConstructCircumcenter>
  Point_3 dual_centroid(Vertex_handle v, ConstructCircumcenter construct_circumcenter) const
  {
    std::list<Edge> edges;
    incident_edges(v, std::back_inserter(edges));

    FT vol(0);
    FT x(0), y(0), z(0);
    for(typename std::list<Edge>::iterator eit = edges.begin();
         eit != edges.end(); ++eit) {

      // compute the dual of the edge *eit but handle the translations
      // with respect to the dual of v. That is why we cannot use one
      // of the existing dual functions here.
      Facet_circulator fstart = incident_facets(*eit);
      Facet_circulator fcit = fstart;
      std::vector<Point_3> pts;
      do {
        // TODO: possible speed-up by caching the circumcenters
        Point_3 dual_orig = periodic_circumcenter(fcit->first, construct_circumcenter).first;
        int idx = fcit->first->index(v);
        Offset off = periodic_point(fcit->first,idx).second;
        pts.push_back(construct_point(std::make_pair(dual_orig,-off)));
        ++fcit;
      } while(fcit != fstart);

      Point_3 orig(0,0,0);
      FT tetvol;
      for(unsigned int i=1; i<pts.size()-1; i++) {
        tetvol = Tetrahedron(orig,pts[0],pts[i],pts[i+1]).volume();
        x += (pts[0].x() + pts[i].x() + pts[i+1].x()) * tetvol;
        y += (pts[0].y() + pts[i].y() + pts[i+1].y()) * tetvol;
        z += (pts[0].z() + pts[i].z() + pts[i+1].z()) * tetvol;
        vol += tetvol;
      }
    }
    x /= ( 4 * vol );
    y /= ( 4 * vol );
    z /= ( 4 * vol );

    Iso_cuboid d = domain();
    x = (x < d.xmin() ? x+d.xmax()-d.xmin()
  : (x >= d.xmax() ? x-d.xmax()+d.xmin() : x));
    y = (y < d.ymin() ? y+d.ymax()-d.ymin()
  : (y >= d.ymax() ? y-d.ymax()+d.ymin() : y));
    z = (z < d.zmin() ? z+d.zmax()-d.zmin()
  : (z >= d.zmax() ? z-d.zmax()+d.zmin() : z));

    CGAL_triangulation_postcondition((d.xmin()<=x)&&(x<d.xmax()));
    CGAL_triangulation_postcondition((d.ymin()<=y)&&(y<d.ymax()));
    CGAL_triangulation_postcondition((d.zmin()<=z)&&(z<d.zmax()));

    return Point_3(x,y,z);
  }
};

template < class GT, class TDS >
inline bool
Periodic_3_triangulation_3<GT,TDS>::
is_triangulation_in_1_sheet() const
{
  if(is_1_cover())
    return true;

  for(Vertex_iterator vit = vertices_begin(); vit != vertices_end(); ++vit) {
    if(virtual_vertices.find(vit) == virtual_vertices.end())
      continue;

    std::vector<Vertex_handle> nb_v;
    std::set<Vertex_handle> nb_v_odom;
    Vertex_handle vh;
    Offset off;
    adjacent_vertices(vit, std::back_inserter(nb_v));
    for(unsigned int i=0; i<nb_v.size(); i++) {
      get_vertex(nb_v[i],vh,off);
      nb_v_odom.insert(vh);
    }
    if(nb_v.size() != nb_v_odom.size())
      return false;
  }
  return true;
}

template < class GT, class TDS >
inline bool
Periodic_3_triangulation_3<GT,TDS>::
is_vertex(const Point& p, Vertex_handle& v) const
{
  Locate_type lt;
  int li, lj;
  Cell_handle c = locate( p, lt, li, lj );
  if( lt != VERTEX )
    return false;
  v = c->vertex(li);
  return true;
}

template < class GT, class TDS >
inline void
Periodic_3_triangulation_3<GT,TDS>::
make_canonical(Vertex_triple& t) const
{
  int i = (t.first < t.second) ? 0 : 1;
  if(i==0) {
    i = (t.first < t.third) ? 0 : 2;
  } else {
    i = (t.second < t.third) ? 1 : 2;
  }
  Vertex_handle tmp;
  switch(i) {
  case 0: return;
  case 1:
    tmp = t.first;
    t.first = t.second;
    t.second = t.third;
    t.third = tmp;
    return;
  default:
    tmp = t.first;
    t.first = t.third;
    t.third = t.second;
    t.second = tmp;
  }
}

/** Assumes a point, an offset, and a cell to start from.
  * Gives the locate type and the simplex (cell and indices) containing p.
  *
  * Performs a remembering stochastic walk if the triangulation is not empty.
  * After the walk the type of the simplex containing p is determined.
  *
  * returns the cell p lies in
  * starts at cell "start"
  * returns a cell Cell_handel if lt == CELL
  * returns a facet (Cell_handle,li) if lt == FACET
  * returns an edge (Cell_handle,li,lj) if lt == EDGE
  * returns a vertex (Cell_handle,li) if lt == VERTEX
  */
template < class GT, class TDS >
inline typename Periodic_3_triangulation_3<GT,TDS>::Cell_handle
Periodic_3_triangulation_3<GT,TDS>::
#ifdef CGAL_NO_STRUCTURAL_FILTERING
periodic_locate
#else
exact_periodic_locate
#endif
(const Point& p, const Offset& o_p, Offset& lo,
 Locate_type& lt, int& li, int& lj, Cell_handle start) const
{
  int cumm_off = 0;
  Offset off_query = o_p;
  if(number_of_vertices() == 0) {
    lo = Offset();
    lt = EMPTY;
    return Cell_handle();
  }
  CGAL_triangulation_assertion(number_of_vertices() != 0);

  if(start == Cell_handle()) {
    start = cells_begin();
  }

  cumm_off = start->offset(0) | start->offset(1)
           | start->offset(2) | start->offset(3);
  if(is_1_cover() && cumm_off != 0) {
    if(((cumm_off & 4) == 4) && (FT(2)*p.x()<(domain().xmax()+domain().xmin())))
      off_query += Offset(1,0,0);
    if(((cumm_off & 2) == 2) && (FT(2)*p.y()<(domain().ymax()+domain().ymin())))
      off_query += Offset(0,1,0);
    if(((cumm_off & 1) == 1) && (FT(2)*p.z()<(domain().zmax()+domain().zmin())))
      off_query += Offset(0,0,1);
  }

  CGAL_triangulation_postcondition(start!=Cell_handle());
  CGAL_triangulation_assertion(start->neighbor(0)->neighbor(
      start->neighbor(0)->index(start))==start);
  CGAL_triangulation_assertion(start->neighbor(1)->neighbor(
      start->neighbor(1)->index(start))==start);
  CGAL_triangulation_assertion(start->neighbor(2)->neighbor(
      start->neighbor(2)->index(start))==start);
  CGAL_triangulation_assertion(start->neighbor(3)->neighbor(
      start->neighbor(3)->index(start))==start);

  // We implement the remembering visibility/stochastic walk.

  // Remembers the previous cell to avoid useless orientation tests.
  Cell_handle previous = Cell_handle();
  Cell_handle c = start;

  // Stores the results of the 4 orientation tests.  It will be used
  // at the end to decide if p lies on a face/edge/vertex/interior.
  Orientation o[4];

  boost::rand48 rng;
  boost::uniform_smallint<> four(0, 3);
  boost::variate_generator<boost::rand48&, boost::uniform_smallint<> > die4(rng, four);

  // Now treat the cell c.
try_next_cell:
  // For the remembering stochastic walk,
  // we need to start trying with a random index :
  int i = die4();
  // For the remembering visibility walk (Delaunay only), we don't :
  // int i = 0;

  cumm_off = c->offset(0) | c->offset(1) | c->offset(2) | c->offset(3);

  bool simplicity_criterion = (cumm_off == 0) && (off_query.is_null());

  // We know that the 4 vertices of c are positively oriented.
  // So, in order to test if p is seen outside from one of c's facets,
  // we just replace the corresponding point by p in the orientation
  // test.  We do this using the arrays below.

  Offset off[4];
  const Point* pts[4] = { &(c->vertex(0)->point()),
                          &(c->vertex(1)->point()),
                          &(c->vertex(2)->point()),
                          &(c->vertex(3)->point()) };

  if(!simplicity_criterion && is_1_cover() ) {
    for(int i=0; i<4; i++) {
      off[i] = int_to_off(c->offset(i));
    }
  }

  if(!is_1_cover()) {
    // Just fetch the vertices of c as points with offsets
    for(int i=0; i<4; i++) {
      pts[i] = &(c->vertex(i)->point());
      off[i] = get_offset(c,i);
    }
  }

  for(int j=0; j != 4; ++j, i = (i+1)&3) {
    Cell_handle next = c->neighbor(i);
    if(previous == next) {
      o[i] = POSITIVE;
      continue;
    }

    CGAL_triangulation_assertion(next->neighbor(next->index(c)) == c);

    // We temporarily put p at i's place in pts.
    const Point* backup = pts[i];
    pts[i] = &p;

    if(simplicity_criterion && is_1_cover() ) {
      o[i] = orientation(*pts[0], *pts[1], *pts[2], *pts[3]);

      if( o[i] != NEGATIVE ) {
        pts[i] = backup;
        continue;
      }
    }
    else {
      Offset backup_off;

      backup_off = off[i];
      off[i] = off_query;
      o[i] = orientation(*pts[0], *pts[1], *pts[2], *pts[3],
          off[0], off[1], off[2], off[3]);

      if( o[i] != NEGATIVE ) {
        pts[i] = backup;
        off[i] = backup_off;
        continue;
      }
    }

    // Test whether we need to adapt the offset of the query point.
    // This means, if we get out of the current cover.
    off_query = combine_offsets(off_query, neighbor_offset(c,i,next));
    previous = c;
    c = next;
    goto try_next_cell;
  }


  // Ok, now we have found the cell. It remains to find the dimension of the
  // intersected simplex.
  // now p is in c or on its boundary
  int sum = ( o[0] == COPLANAR )
          + ( o[1] == COPLANAR )
          + ( o[2] == COPLANAR )
          + ( o[3] == COPLANAR );
  switch(sum) {
  case 0:
      lt = CELL;
      break;
  case 1:
      lt = FACET;
      li = ( o[0] == COPLANAR ) ? 0 :
             ( o[1] == COPLANAR ) ? 1 :
               ( o[2] == COPLANAR ) ? 2 : 3;
      break;
  case 2:
    lt = EDGE;
    li = ( o[0] != COPLANAR ) ? 0 :
           ( o[1] != COPLANAR ) ? 1 : 2;
    lj = ( o[li+1] != COPLANAR ) ? li+1 :
           ( o[li+2] != COPLANAR ) ? li+2 : li+3;
    break;
  case 3:
    lt = VERTEX;
    li = ( o[0] != COPLANAR ) ? 0 :
           ( o[1] != COPLANAR ) ? 1 :
             ( o[2] != COPLANAR ) ? 2 : 3;
    break;
  default:
    // the point cannot lie on four facets
    CGAL_triangulation_assertion(false);
  }
  lo = off_query;

  return c;
}


#ifndef CGAL_NO_STRUCTURAL_FILTERING
template < class GT, class TDS >
typename Periodic_3_triangulation_3<GT,TDS>::Cell_handle
Periodic_3_triangulation_3<GT,TDS>::
inexact_periodic_locate(const Point& p, const Offset& o_p,
                        Cell_handle start,
                        int n_of_turns) const
{
  int cumm_off = 0;
  Offset off_query = o_p;
  if(number_of_vertices() == 0) {
    return Cell_handle();
  }
  CGAL_triangulation_assertion(number_of_vertices() != 0);

  if(start == Cell_handle()) {
    start = cells_begin();
  }

  cumm_off = start->offset(0) | start->offset(1)
          | start->offset(2) | start->offset(3);
  if(is_1_cover() && cumm_off != 0) {
    if(((cumm_off & 4) == 4) && (FT(2)*p.x()<(domain().xmax()+domain().xmin())))
      off_query += Offset(1,0,0);
    if(((cumm_off & 2) == 2) && (FT(2)*p.y()<(domain().ymax()+domain().ymin())))
      off_query += Offset(0,1,0);
    if(((cumm_off & 1) == 1) && (FT(2)*p.z()<(domain().zmax()+domain().zmin())))
      off_query += Offset(0,0,1);
  }

  CGAL_triangulation_postcondition(start!=Cell_handle());
  CGAL_triangulation_assertion(start->neighbor(0)->neighbor(
      start->neighbor(0)->index(start))==start);
  CGAL_triangulation_assertion(start->neighbor(1)->neighbor(
      start->neighbor(1)->index(start))==start);
  CGAL_triangulation_assertion(start->neighbor(2)->neighbor(
      start->neighbor(2)->index(start))==start);
  CGAL_triangulation_assertion(start->neighbor(3)->neighbor(
      start->neighbor(3)->index(start))==start);

  // We implement the remembering visibility/stochastic walk.

  // Remembers the previous cell to avoid useless orientation tests.
  Cell_handle previous = Cell_handle();
  Cell_handle c = start;

  // Now treat the cell c.
try_next_cell:
  --n_of_turns;
  cumm_off = c->offset(0) | c->offset(1) | c->offset(2) | c->offset(3);

  bool simplicity_criterion = (cumm_off == 0) && (off_query.is_null());

  // We know that the 4 vertices of c are positively oriented.
  // So, in order to test if p is seen outside from one of c's facets,
  // we just replace the corresponding point by p in the orientation
  // test.  We do this using the arrays below.

  Offset off[4];
  const Point* pts[4] = { &(c->vertex(0)->point()),
                          &(c->vertex(1)->point()),
                          &(c->vertex(2)->point()),
                          &(c->vertex(3)->point()) };

  if(!simplicity_criterion && is_1_cover() ) {
    for(int i=0; i<4; i++) {
      off[i] = int_to_off(c->offset(i));
    }
  }

  if(!is_1_cover()) {
    // Just fetch the vertices of c as points with offsets
    for(int i=0; i<4; i++) {
      pts[i] = &(c->vertex(i)->point());
      off[i] = get_offset(c,i);
    }
  }

  for(int i=0; i != 4; ++i) {
    Cell_handle next = c->neighbor(i);
    if(previous == next) {
      continue;
    }

    // We temporarily put p at i's place in pts.
    const Point* backup = pts[i];
    pts[i] = &p;

    if(simplicity_criterion && is_1_cover() ) {
      if( inexact_orientation(*pts[0], *pts[1], *pts[2], *pts[3]) != NEGATIVE ) {
        pts[i] = backup;
        continue;
      }
    }
    else {
      Offset backup_off;

      backup_off = off[i];
      off[i] = off_query;

      if( inexact_orientation(*pts[0], *pts[1], *pts[2], *pts[3],
          off[0], off[1], off[2], off[3]) != NEGATIVE ) {
        pts[i] = backup;
        off[i] = backup_off;
        continue;
      }
    }

    // Test whether we need to adapt the offset of the query point.
    // This means, if we get out of the current cover.
    off_query = combine_offsets(off_query, neighbor_offset(c,i,next));
    previous = c;
    c = next;
    if(n_of_turns)
      goto try_next_cell;
  }

  return c;
}
#endif


/**
 * returns
 * ON_BOUNDED_SIDE if (q, off) inside the cell
 * ON_BOUNDARY if (q, off) on the boundary of the cell
 * ON_UNBOUNDED_SIDE if (q, off) lies outside the cell
 *
 * lt has a meaning only when ON_BOUNDED_SIDE or ON_BOUNDARY
 */
template < class GT, class TDS >
inline Bounded_side Periodic_3_triangulation_3<GT,TDS>::side_of_cell(
    const Point& q, const Offset& off, Cell_handle c,
    Locate_type& lt, int& i, int& j) const
{
  CGAL_triangulation_precondition( number_of_vertices() != 0 );

  Orientation o0,o1,o2,o3;
  o0 = o1 = o2 = o3 = ZERO;

  int cumm_off = c->offset(0) | c->offset(1) | c->offset(2) | c->offset(3);
  if((cumm_off == 0) && (is_1_cover())) {
    CGAL_triangulation_assertion(off == Offset());
    const Point& p0  = c->vertex(0)->point();
    const Point& p1  = c->vertex(1)->point();
    const Point& p2  = c->vertex(2)->point();
    const Point& p3  = c->vertex(3)->point();

    if(((o0 = orientation(q ,p1,p2,p3)) == NEGATIVE) ||
        ((o1 = orientation(p0,q ,p2,p3)) == NEGATIVE) ||
        ((o2 = orientation(p0,p1,q ,p3)) == NEGATIVE) ||
        ((o3 = orientation(p0,p1,p2,q )) == NEGATIVE) ) {
      return ON_UNBOUNDED_SIDE;
    }
  } else { // Special case for the periodic space.
    Offset off_q;
    Offset offs[4];
    const Point *p[4];
    for(int i=0; i<4; i++) {
      p[i] = &(c->vertex(i)->point());
      offs[i] = get_offset(c,i);
    }
    CGAL_triangulation_assertion(orientation(*p[0], *p[1], *p[2], *p[3],
                               offs[0], offs[1], offs[2], offs[3]) == POSITIVE);
    bool found=false;
    for(int i=0; (i<8)&&(!found); i++) {
      if((cumm_off | ((~i)&7)) == 7) {
        o0 = o1 = o2 = o3 = NEGATIVE;
        off_q = combine_offsets(off, int_to_off(i));

        if(((o0 = orientation(q,  *p[1],  *p[2],  *p[3],
                               off_q  ,offs[1],offs[2],offs[3])) != NEGATIVE)&&
            ((o1 = orientation(*p[0],      q,  *p[2],  *p[3],
                               offs[0],  off_q,offs[2],offs[3])) != NEGATIVE)&&
            ((o2 = orientation(*p[0],  *p[1],      q,  *p[3],
                               offs[0],offs[1],  off_q,offs[3])) != NEGATIVE)&&
            ((o3 = orientation(*p[0],  *p[1],  *p[2],      q,
                               offs[0],offs[1],offs[2],  off_q)) != NEGATIVE)) {
          found = true;
        }
      }
    }
    if(!found) return ON_UNBOUNDED_SIDE;
  }

  // now all the oi's are >=0
  // sum gives the number of facets p lies on
  int sum = ( (o0 == ZERO) ? 1 : 0 )
            + ( (o1 == ZERO) ? 1 : 0 )
            + ( (o2 == ZERO) ? 1 : 0 )
            + ( (o3 == ZERO) ? 1 : 0 );

  switch(sum) {
  case 0:
    {
      lt = CELL;
      return ON_BOUNDED_SIDE;
    }
  case 1:
    {
      lt = FACET;
      // i = index such that p lies on facet(i)
      i = ( o0 == ZERO ) ? 0 :
            ( o1 == ZERO ) ? 1 :
              ( o2 == ZERO ) ? 2 : 3;
      return ON_BOUNDARY;
    }
  case 2:
    {
      lt = EDGE;
      // i = smallest index such that p does not lie on facet(i)
      // i must be < 3 since p lies on 2 facets
      i = ( o0 == POSITIVE ) ? 0 :
            ( o1 == POSITIVE ) ? 1 : 2;
      // j = larger index such that p not on facet(j)
      // j must be > 0 since p lies on 2 facets
      j = ( o3 == POSITIVE ) ? 3 :
            ( o2 == POSITIVE ) ? 2 : 1;
      return ON_BOUNDARY;
    }
  case 3:
    {
      lt = VERTEX;
      // i = index such that p does not lie on facet(i)
      i = ( o0 == POSITIVE ) ? 0 :
            ( o1 == POSITIVE ) ? 1 :
              ( o2 == POSITIVE ) ? 2 : 3;
      return ON_BOUNDARY;
    }
  default:
    {
      // impossible : cannot be on 4 facets for a real tetrahedron
      CGAL_triangulation_assertion(false);
      return ON_BOUNDARY;
    }
  }
} // side_of_cell

/*! \brief Insert point.
*
* Inserts the point p into the triangulation. It assumes that
* the cell c containing p is already known.
*
* Implementation:
* - some precondition checking
* - find and mark conflicting cells --> find_conflicts
* Conflicting cells are stored in the vector cells.
* - backup hidden points
* - Delete the edges of the marked cells from too_long_edges
* - Insert the new vertex in the hole obtained by removing the
*   conflicting cells (star-approach) --> _tds._insert_in_hole
* - find out about offsets
* - Insert the newly added edges that are "too long"
*   to too_long_edges
* - reinsert hidden points
*/
template < class GT, class TDS >
template < class Conflict_tester, class Point_hider, class CoverManager >
inline typename Periodic_3_triangulation_3<GT,TDS>::Vertex_handle
Periodic_3_triangulation_3<GT,TDS>::periodic_insert(
    const Point& p, const Offset& o,
    Locate_type /*lt*/, Cell_handle c, const Conflict_tester& tester,
    Point_hider& hider, CoverManager& cover_manager, Vertex_handle vh)
{
  Vertex_handle v;
  CGAL_triangulation_precondition(number_of_vertices() != 0);
  CGAL_triangulation_assertion_code(
      Locate_type lt_assert; int i_assert; int j_assert;);
  CGAL_triangulation_precondition(side_of_cell(tester.point(),o, c,
      lt_assert, i_assert, j_assert) != ON_UNBOUNDED_SIDE);

  tester.set_offset(o);

  // Choose the periodic copy of tester.point() that is in conflict with c.
  bool found = false;
  Offset current_off = get_location_offset(tester, c, found);

  CGAL_triangulation_assertion(side_of_cell(tester.point(),
      combine_offsets(o,current_off),c,lt_assert,i_assert,j_assert)
      != ON_UNBOUNDED_SIDE);

  // If the new point is not in conflict with its cell, it is hidden.
  if(!found || !tester.test_initial_cell(c, current_off)) {
    hider.hide_point(c,p);
    return Vertex_handle();
  }

  // Ok, we really insert the point now.
  // First, find the conflict region.
  std::vector<Cell_handle> cells;
  cells.reserve(32);

  Facet facet;

  find_conflicts(c, current_off, tester,
                 make_triple(Oneset_iterator<Facet>(facet),
                             std::back_inserter(cells),
                             Emptyset_iterator()));

  // Remember the points that are hidden by the conflicting cells,
  // as they will be deleted during the insertion.
  hider.set_vertices(cells.begin(), cells.end());

  if(!is_1_cover())
    cover_manager.delete_unsatisfying_elements(cells.begin(), cells.end());

  // Insertion. Warning: facets[0].first MUST be in conflict!
  // Compute the star and put it into the data structure.
  // Store the new cells from the star in nbs.
  v = _tds._insert_in_hole(cells.begin(), cells.end(), facet.first, facet.second);
  v->set_point(p);

  //TODO: this could be done within the _insert_in_hole without losing any
  //time because each cell is visited in any case.
  //- Do timings to argue to modify _insert_in_conflicts if need be
  //- Find the modified _insert_in_hole in the branch svn history of TDS
  std::vector<Cell_handle> nbs;
  incident_cells(v, std::back_inserter(nbs));
  // For all neighbors of the newly added vertex v: fetch their offsets from
  // the tester and reset them in the triangulation data structure.
  for(typename std::vector<Cell_handle>::iterator cit = nbs.begin();
      cit != nbs.end(); cit++) {
    Offset off[4];
    for(int i=0; i<4; i++) {
      off[i] = (*cit)->vertex(i)->offset();
    }
    set_offsets(*cit, off[0], off[1], off[2], off[3]);
  }

  for(typename std::vector<Vertex_handle>::iterator voit = v_offsets.begin();
      voit != v_offsets.end(); ++voit) {
    (*voit)->clear_offset();
  }
  v_offsets.clear();

  if(vh != Vertex_handle()) {
//    CGAL_triangulation_assertion(virtual_vertices.find(v) == virtual_vertices.end());
    virtual_vertices[v] = Virtual_vertex(vh,o);
    virtual_vertices_reverse[vh].push_back(v);
  }

  if(!is_1_cover())
    cover_manager.insert_unsatisfying_elements(v, nbs.begin(), nbs.end());

  // Store the hidden points in their new cells.
  hider.reinsert_vertices(v);
  return v;
}

/** Inserts the first point to a triangulation.
  *
  * With inserting the first point the 3-sheeted covering is constructed.
  * So first, the 27 vertices are inserted and are added to virtual_vertices
  * Then 6*27 cells are created.
  * Then all links are set.
 */
template < class GT, class TDS >
inline typename Periodic_3_triangulation_3<GT,TDS>::Vertex_handle
Periodic_3_triangulation_3<GT,TDS>::create_initial_triangulation(const Point& p)
{
  /// Virtual vertices, one per periodic instance
  Vertex_handle vir_vertices[3][3][3];
  /// Virtual cells, 6 per periodic instance
  Cell_handle cells[3][3][3][6];

  // Initialise vertices:
  vir_vertices[0][0][0] = _tds.create_vertex();
  vir_vertices[0][0][0]->set_point(p);
  virtual_vertices_reverse[vir_vertices[0][0][0]] = std::vector<Vertex_handle>();
  for(int i=0; i<_cover[0]; i++) {
    for(int j=0; j<_cover[1]; j++) {
      for(int k=0; k<_cover[2]; k++) {
        if((i!=0)||(j!=0)||(k!=0)) {
          // Initialise virtual vertices out of the domain for debugging
          vir_vertices[i][j][k] =
            _tds.create_vertex();
          vir_vertices[i][j][k]->set_point(p); //+Offset(i,j,k));
          virtual_vertices[vir_vertices[i][j][k]] =
            Virtual_vertex(vir_vertices[0][0][0], Offset(i,j,k));
          virtual_vertices_reverse[vir_vertices[0][0][0]].push_back(
            vir_vertices[i][j][k]);
        }
        CGAL_triangulation_assertion(vir_vertices[i][j][k] != Vertex_handle());
        CGAL_triangulation_assertion(vir_vertices[0][0][0]->point() == p);
      }
    }
  }

  // Create cells:
  for(int i=0; i<_cover[0]; i++) {
    for(int j=0; j<_cover[1]; j++) {
      for(int k=0; k<_cover[2]; k++) {
        for(int l=0; l<6; l++) {
          // 6 cells per 'cube'
          cells[i][j][k][l] = _tds.create_cell();
          for(int n=0; n<4; n++)
            CGAL_triangulation_assertion(cells[i][j][k][l] != Cell_handle());
        }
      }
    }
  }
  // set vertex and neighbor information
  // index to the right vertex: [number of cells][vertex][offset]
  int vertex_ind[6][4][3] = {
    { {0, 0, 0},  {0, 1, 0},  {0, 0, 1},  {1, 0, 0} },
    { {1, 1, 0},  {0, 1, 1},  {1, 0, 1},  {1, 1, 1} },
    { {1, 0, 0},  {0, 1, 1},  {0, 1, 0},  {0, 0, 1} },
    { {1, 0, 0},  {0, 1, 1},  {0, 0, 1},  {1, 0, 1} },
    { {1, 0, 0},  {0, 1, 1},  {1, 0, 1},  {1, 1, 0} },
    { {1, 0, 0},  {0, 1, 1},  {1, 1, 0},  {0, 1, 0} }
  };
  int neighb_ind[6][4][4] = {
    { { 0, 0, 0, 2},  { 0,-1, 0, 5},  { 0, 0,-1, 3},  {-1, 0, 0, 4} },
    { { 0, 0, 1, 5},  { 1, 0, 0, 2},  { 0, 1, 0, 3},  { 0, 0, 0, 4} },
    { {-1, 0, 0, 1},  { 0, 0, 0, 0},  { 0, 0, 0, 3},  { 0, 0, 0, 5} },
    { { 0, 0, 1, 0},  { 0,-1, 0, 1},  { 0, 0, 0, 4},  { 0, 0, 0, 2} },
    { { 0, 0, 0, 1},  { 1, 0, 0, 0},  { 0, 0, 0, 5},  { 0, 0, 0, 3} },
    { { 0, 1, 0, 0},  { 0, 0,-1, 1},  { 0, 0, 0, 2},  { 0, 0, 0, 4} }
  };
  for(int i=0; i<_cover[0]; i++) {
    for(int j=0; j<_cover[1]; j++) {
      for(int k=0; k<_cover[2]; k++) {
        int offset = (i==_cover[0]-1 ? 4 : 0) |
                     (j==_cover[1]-1 ? 2 : 0) |
                     (k==_cover[2]-1 ? 1 : 0);
        for(int l=0; l<6; l++) {
          // cell 0:
          cells[i][j][k][l]->set_vertices(
              vir_vertices
                [(i+vertex_ind[l][0][0])%_cover[0]]
                [(j+vertex_ind[l][0][1])%_cover[1]]
                [(k+vertex_ind[l][0][2])%_cover[2]],
              vir_vertices
                [(i+vertex_ind[l][1][0])%_cover[0]]
                [(j+vertex_ind[l][1][1])%_cover[1]]
                [(k+vertex_ind[l][1][2])%_cover[2]],
              vir_vertices
                [(i+vertex_ind[l][2][0])%_cover[0]]
                [(j+vertex_ind[l][2][1])%_cover[1]]
                [(k+vertex_ind[l][2][2])%_cover[2]],
              vir_vertices
                [(i+vertex_ind[l][3][0])%_cover[0]]
                [(j+vertex_ind[l][3][1])%_cover[1]]
                [(k+vertex_ind[l][3][2])%_cover[2]]);
          set_offsets(cells[i][j][k][l],
              offset & (vertex_ind[l][0][0]*4 +
                        vertex_ind[l][0][1]*2 +
                        vertex_ind[l][0][2]*1),
              offset & (vertex_ind[l][1][0]*4 +
                        vertex_ind[l][1][1]*2 +
                        vertex_ind[l][1][2]*1),
              offset & (vertex_ind[l][2][0]*4 +
                        vertex_ind[l][2][1]*2 +
                        vertex_ind[l][2][2]*1),
              offset & (vertex_ind[l][3][0]*4 +
                        vertex_ind[l][3][1]*2 +
                        vertex_ind[l][3][2]*1));
          cells[i][j][k][l]->set_neighbors(
              cells [(i+_cover[0]+neighb_ind[l][0][0])%_cover[0]]
                    [(j+_cover[1]+neighb_ind[l][0][1])%_cover[1]]
                    [(k+_cover[2]+neighb_ind[l][0][2])%_cover[2]]
                    [         neighb_ind[l][0][3]       ],
              cells [(i+_cover[0]+neighb_ind[l][1][0])%_cover[0]]
                    [(j+_cover[1]+neighb_ind[l][1][1])%_cover[1]]
                    [(k+_cover[2]+neighb_ind[l][1][2])%_cover[2]]
                    [         neighb_ind[l][1][3]       ],
              cells [(i+_cover[0]+neighb_ind[l][2][0])%_cover[0]]
                    [(j+_cover[1]+neighb_ind[l][2][1])%_cover[1]]
                    [(k+_cover[2]+neighb_ind[l][2][2])%_cover[2]]
                    [         neighb_ind[l][2][3]       ],
              cells [(i+_cover[0]+neighb_ind[l][3][0])%_cover[0]]
                    [(j+_cover[1]+neighb_ind[l][3][1])%_cover[1]]
                    [(k+_cover[2]+neighb_ind[l][3][2])%_cover[2]]
                    [         neighb_ind[l][3][3]       ]
          );
        }
      }
    }
  }
  // set pointers from the vertices to incident cells.
  for(int i=0; i<_cover[0]; i++) {
    for(int j=0; j<_cover[1]; j++) {
      for(int k=0; k<_cover[2]; k++) {
        vir_vertices[i][j][k]->set_cell(cells[i][j][k][0]);
      }
    }
  }

  _tds.set_dimension(3);

  return vir_vertices[0][0][0];
}
#define CGAL_INCLUDE_FROM_PERIODIC_3_TRIANGULATION_3_H
#include <CGAL/internal/Periodic_3_triangulation_dummy_36.h>
#undef CGAL_INCLUDE_FROM_PERIODIC_3_TRIANGULATION_3_H

/** finds all cells that are in conflict with the currently added point
  * (stored in tester).
  *
  * The result will be a hole of which the following data is returned:
  * - boundary facets
  * - cells
  * - internal facets.
  *
  * c is the current cell, which must be in conflict.
  * tester is the function object that tests if a cell is in conflict.
  */
template <class GT, class TDS>
template <class Conflict_test,
          class OutputIteratorBoundaryFacets,
          class OutputIteratorCells,
          class OutputIteratorInternalFacets>
Triple<OutputIteratorBoundaryFacets,
       OutputIteratorCells,
       OutputIteratorInternalFacets>
Periodic_3_triangulation_3<GT,TDS>::
find_conflicts(Cell_handle d, const Offset& current_off,
               const Conflict_test& tester,
               Triple<OutputIteratorBoundaryFacets,
               OutputIteratorCells,
               OutputIteratorInternalFacets> it) const
{
  CGAL_triangulation_precondition( number_of_vertices() != 0 );
  CGAL_triangulation_precondition( tester(d, current_off) );

  std::stack<std::pair<Cell_handle, Offset> > cell_stack;
  cell_stack.push(std::make_pair(d,current_off));
  d->tds_data().mark_in_conflict();
  *it.second++ = d;

  do {
    Cell_handle c = cell_stack.top().first;
    Offset current_off2 = cell_stack.top().second;
    cell_stack.pop();

    for(int i=0; i< 4; ++i) {
      Cell_handle test = c->neighbor(i);
      if(test->tds_data().is_in_conflict()) {
        if(c < test) {
          *it.third++ = Facet(c, i); // Internal facet.
        }
        continue; // test was already in conflict.
      }
      if(test->tds_data().is_clear()) {
        Offset o_test = current_off2 + neighbor_offset(c, i, test);
        if(tester(test,o_test)) {
          if(c < test)
            *it.third++ = Facet(c, i); // Internal facet.

          cell_stack.push(std::make_pair(test,o_test));
          test->tds_data().mark_in_conflict();
          *it.second++ = test;
          continue;
        }
        test->tds_data().mark_on_boundary(); // test is on the boundary.
      }
      *it.first++ = Facet(c, i);
      for(int j = 0; j<4; j++) {
        if(j==i) continue;
        if(!c->vertex(j)->get_offset_flag()) {
          c->vertex(j)->set_offset(int_to_off(c->offset(j))-current_off2);
          v_offsets.push_back(c->vertex(j));
        }
      }
    }
  } while(!cell_stack.empty());
  return it;
}

/*! \brief Insert point into triangulation.
 *
 * Inserts the point p into the triangulation. It expects
 * - a cell to start the point location
 * - a testing function to determine cells in conflict
 * - a testing function to determine if a vertex is hidden.
 *
 * Implementation:
 * - If the triangulation is empty call a special function
 * (create_initial_triangulation) to construct the basic
 * triangulation.
 * - Run point location to get the cell c containing point p.
 * - Call periodic_insert to insert p into the 3-cover.
 * - Also insert the eight periodic copies of p.
 */
template < class GT, class TDS >
template < class Conflict_tester, class Point_hider, class CoverManager >
inline typename Periodic_3_triangulation_3<GT,TDS>::Vertex_handle
Periodic_3_triangulation_3<GT,TDS>::insert_in_conflict(const Point& p,
    Locate_type lt, Cell_handle c, int li, int lj,
    const Conflict_tester& tester, Point_hider& hider, CoverManager& cover_manager) {
  CGAL_triangulation_assertion((domain().xmin() <= p.x())
                               && (p.x() < domain().xmax()));
  CGAL_triangulation_assertion((domain().ymin() <= p.y())
                               && (p.y() < domain().ymax()));
  CGAL_triangulation_assertion((domain().zmin() <= p.z())
                               && (p.z() < domain().zmax()));

  if(number_of_vertices() == 0) {
    Vertex_handle vh = create_initial_triangulation(p);
    cover_manager.create_initial_triangulation();
    return vh;
  }

  if((lt == VERTEX) &&
      (tester.compare_weight(c->vertex(li)->point(),p) == 0) ) {
    return c->vertex(li);
  }

  Vertex_handle vstart;
  if(!is_1_cover()) {
    Virtual_vertex_map_it vvmit = virtual_vertices.find(c->vertex(0));
    if(vvmit == virtual_vertices.end())
      vstart = c->vertex(0);
    else
      vstart = vvmit->second.first;
    CGAL_triangulation_assertion(virtual_vertices.find(vstart)
                                 == virtual_vertices.end());
    CGAL_triangulation_assertion(virtual_vertices_reverse.find(vstart)
                                 != virtual_vertices_reverse.end());
  }

  CGAL_triangulation_assertion( number_of_vertices() != 0 );
  CGAL_triangulation_expensive_assertion(is_valid());
  hider.set_original_cube(true);
  Vertex_handle vh = periodic_insert(p, Offset(), lt, c, tester, hider, cover_manager);
  if(is_1_cover()) {
    return vh;
  }

  hider.set_original_cube(false);

  virtual_vertices_reverse[vh] = std::vector<Vertex_handle>();
  Offset lo;
  // insert 26 periodic copies
  for(int i=0; i<_cover[0]; i++) {
    for(int j=0; j<_cover[1]; j++) {
      for(int k=0; k<_cover[2]; k++) {
        if((i!=0)||(j!=0)||(k!=0)) {
          c = periodic_locate(p, Offset(i,j,k), lo, lt, li, lj, Cell_handle());
          periodic_insert(p, Offset(i,j,k), lt, c, tester, hider,cover_manager,vh);
        }
      }
    }
  }
  CGAL_triangulation_expensive_assertion(is_valid());

  // Fall back to 1-cover if the criterion that the longest edge is shorter
  // than sqrt(0.166) is fulfilled.
  if( cover_manager.can_be_converted_to_1_sheet() ) {
    CGAL_triangulation_expensive_assertion(is_valid());
    convert_to_1_sheeted_covering();
    CGAL_triangulation_expensive_assertion( is_valid() );
  }
  return vh;
}

/// tests if two vertices of one cell are just periodic copies of each other
template < class GT, class TDS >
inline bool Periodic_3_triangulation_3<GT,TDS>::has_self_edges(Cell_handle c) const {
  CGAL_triangulation_assertion((c->vertex(0) != c->vertex(1)) ||
                               (c->offset(0) != c->offset(1)));
  CGAL_triangulation_assertion((c->vertex(0) != c->vertex(2)) ||
                               (c->offset(0) != c->offset(2)));
  CGAL_triangulation_assertion((c->vertex(0) != c->vertex(3)) ||
                               (c->offset(0) != c->offset(3)));
  CGAL_triangulation_assertion((c->vertex(1) != c->vertex(2)) ||
                               (c->offset(1) != c->offset(2)));
  CGAL_triangulation_assertion((c->vertex(1) != c->vertex(3)) ||
                               (c->offset(1) != c->offset(3)));
  CGAL_triangulation_assertion((c->vertex(2) != c->vertex(3)) ||
                               (c->offset(2) != c->offset(3)));
  return ((c->vertex(0) == c->vertex(1)) ||
          (c->vertex(0) == c->vertex(2)) ||
          (c->vertex(0) == c->vertex(3)) ||
      (c->vertex(1) == c->vertex(2)) ||
      (c->vertex(1) == c->vertex(3)) ||
      (c->vertex(2) == c->vertex(3)));
}

/*! \brief Tests if the triangulation is valid.
 *
 * A triangulation is valid if
 * - A cell is not its own neighbor.
 * - A cell has no two equal neighbors
 * - A cell has no two equal vertex-offset pairs
 * - A cell is positively oriented.
 * - The point of a neighbor of cell c that does not belong to c is not inside
 *   the circumcircle of c.
 */
template < class GT, class TDS >
bool
Periodic_3_triangulation_3<GT,TDS>::
is_valid(bool verbose, int level) const
{
  if(!is_1_cover())
  {
    for(Virtual_vertex_reverse_map_it iter = virtual_vertices_reverse.begin(), end_iter = virtual_vertices_reverse.end();
         iter != end_iter;
         ++iter)
    {
      for(typename Virtual_vertex_reverse_map::mapped_type::const_iterator iter_2 = iter->second.begin(),
           end_iter_2 = iter->second.end();
           iter_2 != end_iter_2;
           ++iter_2)
      {
        CGAL_triangulation_assertion(virtual_vertices.find(*iter_2) != virtual_vertices.end());
        CGAL_triangulation_assertion(virtual_vertices.at(*iter_2).first == iter->first);
      }
    }
  }

  bool error = false;
  for(Cell_iterator cit = cells_begin();
       cit != cells_end(); ++cit) {
    for(int i=0; i<4; i++) {
      CGAL_triangulation_assertion(cit != cit->neighbor(i));
      for(int j=i+1; j<4; j++) {
        CGAL_triangulation_assertion(cit->neighbor(i) != cit->neighbor(j));
        CGAL_triangulation_assertion(cit->vertex(i) != cit->vertex(j));
      }
    }
    // Check positive orientation:
    const Point *p[4]; Offset off[4];
    for(int i=0; i<4; i++) {
      p[i] = &cit->vertex(i)->point();
      off[i] = get_offset(cit,i);
    }
    if(orientation(*p[0], *p[1], *p[2], *p[3],
                   off[0], off[1], off[2], off[3]) != POSITIVE) {
      if(verbose) {
        std::cerr<<"Periodic_3_triangulation_3: wrong orientation:"<<std::endl;
        std::cerr<<off[0]<<'\t'<<*p[0]<<'\n'
                         <<off[1]<<'\t'<<*p[1]<<'\n'
                         <<off[2]<<'\t'<<*p[2]<<'\n'
                         <<off[3]<<'\t'<<*p[3]<<std::endl;
      }
      error = true;
    }
  }

  if(!has_self_edges()) {
    if(! _tds.is_valid(verbose, level) ) {
      return false;
    }
  }

  return !error;
}

template < class GT, class TDS >
bool Periodic_3_triangulation_3<GT,TDS>::is_valid(Cell_handle ch,
                                                  bool verbose, int level) const
{
  if( ! _tds.is_valid(ch,verbose,level) )
    return false;

  bool error = false;
  const Point *p[4]; Offset off[4];
  for(int i=0; i<4; i++) {
    p[i] = &ch->vertex(i)->point();
    off[i] = get_offset(ch,i);
  }

  if(orientation(*p[0], *p[1], *p[2], *p[3],
                  off[0], off[1], off[2], off[3]) != POSITIVE) {
    error = true;
  }

  return !error;
}

template < class GT, class TDS >
template < class ConflictTester >
bool Periodic_3_triangulation_3<GT,TDS>::
is_valid_conflict(ConflictTester& tester, bool verbose, int level) const
{
  Cell_iterator it;
  for( it = cells_begin(); it != cells_end(); ++it ) {
    is_valid(it, verbose, level);
    for(int i=0; i<4; i++ ) {
      Offset o_nb = neighbor_offset(it, i, it->neighbor(i));
      Offset o_vt = get_offset(it->neighbor(i), it->neighbor(i)->index(it));
      if(tester(it,
                it->neighbor(i)->vertex(it->neighbor(i)->index(it))->point(),
                o_vt - o_nb)) {
        if(verbose) {
          std::cerr << "non-empty sphere:\n"
                    << "Point[0]: " << it->vertex(0)->point()
                    << " Off: " << int_to_off(it->offset(0)) << "\n"
                    << "Point[1]: " << it->vertex(1)->point()
                    << " Off: " << int_to_off(it->offset(1)) << "\n"
                    << "Point[2]: " << it->vertex(2)->point()
                    << " Off: " << int_to_off(it->offset(2)) << "\n"
                    << "Point[3]: " << it->vertex(3)->point()
                    << " Off: " << int_to_off(it->offset(3)) << "\n"
                    << "Foreign point: "
                    << it->neighbor(i)->vertex(it->neighbor(i)->index(it))->point()
                    << "\t Off: " << o_vt - o_nb << std::endl;
        }
        return false;
      }
    }
  }
  return true;
}

template < class GT, class TDS >
inline void Periodic_3_triangulation_3<GT,TDS>::make_hole(Vertex_handle v,
                                      std::map<Vertex_triple,Facet>& outer_map,
                                      std::vector<Cell_handle>& hole)
{
  incident_cells(v, std::back_inserter(hole));

  for(typename std::vector<Cell_handle>::iterator cit = hole.begin();
       cit != hole.end(); ++cit) {
    int indv = (*cit)->index(v);
    Cell_handle opp_cit = (*cit)->neighbor( indv );
    Facet f(opp_cit, opp_cit->index(*cit));
    Vertex_triple vt = make_vertex_triple(f);
    make_canonical(vt);
    outer_map[vt] = f;
    for(int i=0; i<4; i++)
      if( i != indv )
        (*cit)->vertex(i)->set_cell(opp_cit);
  }
}

/*! \brief Remove a vertex from the triangulation.
 *
 * Removes vertex v from the triangulation.
 */
template < class GT, class TDS >
template < class PointRemover, class Conflict_tester, class CoverManager>
inline void
Periodic_3_triangulation_3<GT,TDS>::
remove(Vertex_handle v, PointRemover& r, Conflict_tester& t, CoverManager& cover_manager)
{
  CGAL_expensive_precondition(is_vertex(v));

  if(!is_1_cover())
  {
    std::vector<Vertex_handle> vhrem;
    if(number_of_vertices() == 1) {
      clear();
      return;
    }
    Virtual_vertex_map_it vvmit = virtual_vertices.find(v);
    if(vvmit != virtual_vertices.end()) v = vvmit->second.first;
    CGAL_triangulation_assertion(virtual_vertices_reverse.find(v)
                                 != virtual_vertices_reverse.end());
    vhrem = virtual_vertices_reverse.find(v)->second;
    virtual_vertices_reverse.erase(v);
    CGAL_triangulation_assertion(vhrem.size()==26);
    for(int i=0; i<26; i++) {
      periodic_remove(vhrem[i],r, cover_manager);
      virtual_vertices.erase(vhrem[i]);
      CGAL_triangulation_expensive_assertion(is_valid());
    }
    periodic_remove(v,r, cover_manager);
  }
  else
  {
    periodic_remove(v, r, cover_manager);
    if(!is_1_cover())
      remove(v, r, t, cover_manager);
  }
}

/*! \brief Remove a vertex from the triangulation.
 *
 * Removes vertex v from the triangulation.
 * It expects a reference to an instance of a PointRemover.
 *
 * Implementation:
 * - Compute the hole, that is, all cells incident to v. Cells outside of
 *   this hole are not affected by the deletion of v.
 * - Triangulate the hole. This is done computing the triangulation
 *   in Euclidean space for the points on the border of the hole.
 * - Sew this triangulation into the hole.
 * - Test for all newly added edges, whether they are shorter than the
 *   edge_length_threshold. If not, convert to 3-cover.
 */
template < class GT, class TDS >
template < class PointRemover, class CoverManager >
inline bool
Periodic_3_triangulation_3<GT,TDS>::
periodic_remove(Vertex_handle v, PointRemover& remover, CoverManager& cover_manager,
                const bool abort_if_cover_change)
{
  // Construct the set of vertex triples on the boundary
  // with the facet just behind
  typedef std::map<Vertex_triple,Facet>           Vertex_triple_Facet_map;
  typedef PointRemover                            Point_remover;
  typedef typename Point_remover::CellE_handle    CellE_handle;
  typedef typename Point_remover::VertexE_handle  VertexE_handle;
  typedef typename Point_remover::FacetE          FacetE;
  typedef typename Point_remover::VertexE_triple  VertexE_triple;
  typedef typename Point_remover::Finite_cellsE_iterator
      Finite_cellsE_iterator;
  typedef typename Point_remover::Vertex_triple_FacetE_map
      Vertex_triple_FacetE_map;

  // First compute the hole and its boundary vertices.
  std::vector<Cell_handle> hole;
  hole.reserve(64);
  Vertex_triple_Facet_map outer_map;
  Vertex_triple_FacetE_map inner_map;

  make_hole(v, outer_map, hole);

  CGAL_triangulation_assertion(outer_map.size()==hole.size());

  if(!is_1_cover()) {
    cover_manager.delete_unsatisfying_elements(hole.begin(), hole.end());
  }

  // Build up the map between Vertices on the boundary and offsets
  // collect all vertices on the boundary
  std::vector<Vertex_handle> vertices;
  vertices.reserve(64);

  // The set is needed to ensure that each vertex is inserted only once.
  std::set<Vertex_handle> tmp_vertices;
  // The map connects vertices to offsets in the hole
  std::map<Vertex_handle, Offset> vh_off_map;

  for(typename std::vector<Cell_handle>::iterator cit = hole.begin();
      cit != hole.end(); ++cit)
  {
    // Put all incident vertices in tmp_vertices.
    for(int j=0; j<4; ++j) {
      if((*cit)->vertex(j) != v) {
        tmp_vertices.insert((*cit)->vertex(j));
        vh_off_map[(*cit)->vertex(j)] = int_to_off((*cit)->offset(j))
                                        - int_to_off((*cit)->offset((*cit)->index(v)));
      }
    }
  }

  // Now output the vertices.
  std::copy(tmp_vertices.begin(), tmp_vertices.end(), std::back_inserter(vertices));

  // create a Delaunay/regular triangulation of the points on the boundary
  // in Euclidean space and make a map from the vertices in remover.tmp
  // towards the vertices in *this

  Unique_hash_map<VertexE_handle,Vertex_handle> vmap;
  CellE_handle ch;
  remover.tmp.clear();

  for(unsigned int i=0; i < vertices.size(); i++) {
    typedef typename Point_remover::Triangulation_R3::Point TRPoint;
    CGAL_triangulation_assertion(
          get_offset(vertices[i]) + combine_offsets(Offset(), vh_off_map[vertices[i]])
       == combine_offsets(get_offset(vertices[i]),vh_off_map[vertices[i]]));
    TRPoint trp = std::make_pair(vertices[i]->point(),
                                 combine_offsets( get_offset(vertices[i]), vh_off_map[vertices[i]]) );
    VertexE_handle vh = remover.tmp.insert(trp, ch);
    vmap[vh] = vertices[i];
    CGAL_triangulation_assertion(vmap.is_defined(vh));
  }
  CGAL_triangulation_assertion(remover.tmp.number_of_vertices() != 0);

  // Construct the set of vertex triples of tmp
  // We reorient the vertex triple so that it matches those from outer_map
  // Also note that we use the vertices of *this, not of tmp
  for(Finite_cellsE_iterator it = remover.tmp.finite_cells_begin();
      it != remover.tmp.finite_cells_end();
      ++it) {
    VertexE_triple vt_aux;
    for(int i=0; i < 4; i++) {
      FacetE f = std::pair<CellE_handle,int>(it,i);
      vt_aux = VertexE_triple(
          f.first->vertex(vertex_triple_index(f.second,0)),
          f.first->vertex(vertex_triple_index(f.second,1)),
          f.first->vertex(vertex_triple_index(f.second,2)));
      if(vmap.is_defined(vt_aux.first)
          && vmap.is_defined(vt_aux.second)
          && vmap.is_defined(vt_aux.third) ) {
        Vertex_triple vt(vmap[vt_aux.first],vmap[vt_aux.third],
            vmap[vt_aux.second]);
        make_canonical(vt);
        inner_map[vt]= f;
      }
    }
  }

  // A structure for storing the new neighboring relations
  typedef boost::tuple<Cell_handle, int, Cell_handle> Neighbor_relation;
  std::vector<Neighbor_relation> nr_vec;
  std::vector<Cell_handle> new_cells;

  // Grow inside the hole, by extending the surface
  while(! outer_map.empty()) {
    typename Vertex_triple_Facet_map::iterator oit = outer_map.begin();

    typename Vertex_triple_Facet_map::value_type o_vt_f_pair = *oit;
    Cell_handle o_ch = o_vt_f_pair.second.first;
    unsigned int o_i = o_vt_f_pair.second.second;

    typename Vertex_triple_FacetE_map::iterator iit =
        inner_map.find(o_vt_f_pair.first);

    CGAL_triangulation_assertion(iit != inner_map.end());
    typename Vertex_triple_FacetE_map::value_type i_vt_f_pair = *iit;
    CellE_handle i_ch = i_vt_f_pair.second.first;
    unsigned int i_i = i_vt_f_pair.second.second;

    // create a new cell to glue to the outer surface
    Cell_handle new_ch = _tds.create_cell();
    new_cells.push_back(new_ch);
    new_ch->set_vertices(vmap[i_ch->vertex(0)], vmap[i_ch->vertex(1)],
                         vmap[i_ch->vertex(2)], vmap[i_ch->vertex(3)]);
    set_offsets(new_ch, vh_off_map[vmap[i_ch->vertex(0)]],
                        vh_off_map[vmap[i_ch->vertex(1)]],
                        vh_off_map[vmap[i_ch->vertex(2)]],
                        vh_off_map[vmap[i_ch->vertex(3)]]);

    // Update the edge length management
    if(cover_manager.update_cover_data_during_management(new_ch, new_cells,
                                                         abort_if_cover_change))
    {
      CGAL_triangulation_expensive_postcondition(_tds.is_valid());
      return false; // removing would cause / has caused a cover change
    }

    // The neighboring relation needs to be stored temporarily in
    // nr_vec. It cannot be applied directly because then we could not
    // easily cancel the removing process if a cell is encountered
    // that does not obey the edge-length criterion.
    nr_vec.push_back(boost::make_tuple(o_ch,o_i,new_ch));
    nr_vec.push_back(boost::make_tuple(new_ch,i_i,o_ch));

    // for the other faces check, if they can also be glued
    for(unsigned int i = 0; i < 4; i++) {
      if(i != i_i) {
        Facet f = std::pair<Cell_handle,int>(new_ch,i);
        Vertex_triple vt = make_vertex_triple(f);
        make_canonical(vt);
        std::swap(vt.second,vt.third);
        typename Vertex_triple_Facet_map::iterator oit2 = outer_map.find(vt);
        if(oit2 == outer_map.end()) {
          std::swap(vt.second,vt.third);
          outer_map[vt]= f;
        } else {
        // glue the faces
          typename Vertex_triple_Facet_map::value_type o_vt_f_pair2 = *oit2;
          Cell_handle o_ch2 = o_vt_f_pair2.second.first;
          int o_i2 = o_vt_f_pair2.second.second;
          nr_vec.push_back(boost::make_tuple(o_ch2,o_i2,new_ch));
          nr_vec.push_back(boost::make_tuple(new_ch,i,o_ch2));
          outer_map.erase(oit2);
        }
      }
    }
    outer_map.erase(oit);
  }

  // finally set the neighboring relations
  for(unsigned int i=0; i<nr_vec.size(); i++) {
    nr_vec[i].template get<0>()->set_neighbor(nr_vec[i].template get<1>(),nr_vec[i].template get<2>());
  }

  // Output the hidden points.
  for(typename std::vector<Cell_handle>::iterator
      hi = hole.begin(), hend = hole.end(); hi != hend; ++hi)
  {
    remover.add_hidden_points(*hi);
  }

  _tds.delete_vertex(v);
  _tds.delete_cells(hole.begin(), hole.end());
  CGAL_triangulation_expensive_assertion(is_valid());
  return true; // sucessfully removed the vertex
}

// ############################################################################
// ############################################################################
// ############################################################################

/** \brief Delete each redundant cell and the not anymore needed data
 *  structures.
 *
 *  This function consists of four iterations over all cells and one
 *  iteration over all vertices:
 *  1. cell iteration: mark all cells that are to delete
 *  2. cell iteration: redirect neighbors of remaining cells
 *  3. cell iteration: redirect vertices of remaining cells
 *  4. cell iteration: delete all cells marked in the 1. iteration
 *  Vertex iteration: delete all vertices outside the original domain.
 */
template < class GT, class TDS >
inline void
Periodic_3_triangulation_3<GT,TDS>::convert_to_1_sheeted_covering()
{
  // ###################################################################
  // ### First cell iteration ##########################################
  // ###################################################################
  {
    if(is_1_cover())
      return;

    bool to_delete, has_simplifiable_offset;
    Virtual_vertex_map_it vvmit;
    // First iteration over all cells: Mark the cells that are to be deleted.
    // Cells will be deleted if they cannot be translated anymore in the
    // direction of one of the axes without yielding negative offsets.
    for( Cell_iterator it = all_cells_begin();
    it != all_cells_end(); ++it ) {
      to_delete = false;
      // for all directions in 3D Space
      for( int j=0; j<3; j++ ) {
        has_simplifiable_offset = true;
        // for all vertices of cell it
        for( int i=0; i<4; i++ ) {
          vvmit = virtual_vertices.find(it->vertex(i));
          // if it->vertex(i) lies inside the original domain:
          if(vvmit == virtual_vertices.end()) {
            // the cell cannot be moved any more because if we did, then
            // it->vertex(i) will get at least one negative offset.
            has_simplifiable_offset = false;
            // if it->vertex(i) lies outside the original domain:
          } else {
            // The cell can certainly be deleted if the offset contains a 2
            to_delete = to_delete || (vvmit->second.second[j] == 2);
            // The cell can be moved into one direction only if the offset of
            // all for vertices is >=1 for this direction. Since we already
            // tested for 2 it is sufficient to test here for 1.
            has_simplifiable_offset = has_simplifiable_offset
                && (vvmit->second.second[j] == 1);
          }
        }
        // if the offset can be simplified, i.e. the cell can be moved, then
        // it can be deleted.
        if(has_simplifiable_offset)
          to_delete = true;
      }
      // Mark all cells that are to delete. They cannot be deleted yet,
      // because neighboring information still needs to be extracted.
      if(to_delete) {
        it->set_additional_flag(1);
      }
    }
  }

  // ###################################################################
  // ### Second cell iteration #########################################
  // ###################################################################
  {
    Vertex_handle vert[4], nbv[4];
    Offset off[4];
    Cell_handle nb, new_neighbor;
    std::vector<Triple<Cell_handle, int, Cell_handle> > new_neighbor_relations;

    // Second iteration over all cells: redirect neighbors where necessary
    for(Cell_iterator it = all_cells_begin();
        it != all_cells_end(); ++it) {
      // Skip all cells that are to delete.
      if(it->get_additional_flag() == 1)
        continue;

      // Redirect neighbors: Only neighbors that are marked by the
      // additional_flag have to be substituted by one of their periodic
      // copies. The unmarked neighbors stay the same.
      for( int i = 0; i < 4; i++ ) {
        if( it->neighbor(i)->get_additional_flag() != 1 )
          continue;

        nb = it->neighbor(i);

        for( int j = 0; j < 4; j++ ) {
          off[j] = Offset();
          get_vertex( nb, j, vert[j], off[j]);
        }
        int x,y,z;
        x = (std::min) ( (std::min) ( off[0][0], off[1][0] ),
            (std::min) ( off[2][0], off[3][0] ) );
        y = (std::min) ( (std::min) ( off[0][1], off[1][1] ),
            (std::min) ( off[2][1], off[3][1] ) );
        z = (std::min) ( (std::min) ( off[0][2], off[1][2] ),
            (std::min) ( off[2][2], off[3][2] ) );

        // The vector from nb to the "original" periodic copy of nb, that is
        // the copy that will not be deleted.
        Offset difference_offset(x,y,z);
        CGAL_triangulation_assertion( !difference_offset.is_null() );

        // We now have to find the "original" periodic copy of nb from
        // its vertices. Therefore, we first have to find the vertices.
        for( int j = 0; j < 4; j++ ) {
          CGAL_triangulation_assertion( (off[j]-difference_offset)[0] >= 0);
          CGAL_triangulation_assertion( (off[j]-difference_offset)[1] >= 0);
          CGAL_triangulation_assertion( (off[j]-difference_offset)[2] >= 0);
          CGAL_triangulation_assertion( (off[j]-difference_offset)[0] < 3);
          CGAL_triangulation_assertion( (off[j]-difference_offset)[1] < 3);
          CGAL_triangulation_assertion( (off[j]-difference_offset)[2] < 3);

          // find the Vertex_handles of the vertices of the "original"
          // periodic copy of nb. If the vertex is inside the original
          // domain, there is nothing to do
          if( (off[j]-difference_offset).is_null() ) {
            nbv[j] = vert[j];
            // If the vertex is outside the original domain, we have to search
            // in virtual_vertices in the "wrong" direction. That means, we
            // cannot use virtual_vertices.find but have to use
            // virtual_vertices_reverse.
          } else {
            Offset nbo = off[j]-difference_offset;
            nbv[j] = virtual_vertices_reverse.find(vert[j])
              ->second[nbo[0]*9+nbo[1]*3+nbo[2]-1];
          }
        }
        // Find the new neighbor by its 4 vertices
        new_neighbor = get_cell( nbv );

        // Store the new neighbor relation. This cannot be applied yet because
        // it would disturb the functioning of get_cell( ... )
        new_neighbor_relations.push_back(make_triple(it, i, new_neighbor));
      }
    }
    // Apply the new neighbor relations now.
    for(unsigned int i=0; i<new_neighbor_relations.size(); i++) {
      new_neighbor_relations[i].first->set_neighbor(
          new_neighbor_relations[i].second,
          new_neighbor_relations[i].third);
    }
  }

  // ###################################################################
  // ### Third cell iteration ##########################################
  // ###################################################################
  {
    Vertex_handle vert[4];
    Offset off[4];
    // Third iteration over all cells: redirect vertices where necessary
    for(Cell_iterator it = all_cells_begin();
        it != all_cells_end(); ++it) {
      // Skip all cells that are marked to delete
      if(it->get_additional_flag() == 1)
        continue;

      // Find the corresponding vertices of it in the original domain
      // and set them as new vertices of it.
      for( int i = 0; i < 4; i++ ) {
        off[i] = Offset();
        get_vertex( it, i, vert[i], off[i]);
        it->set_vertex( i, vert[i]);
        CGAL_triangulation_assertion(vert[i]->point()[0] < domain().xmax());
        CGAL_triangulation_assertion(vert[i]->point()[1] < domain().ymax());
        CGAL_triangulation_assertion(vert[i]->point()[2] < domain().zmax());
        CGAL_triangulation_assertion(vert[i]->point()[0] >= domain().xmin());
        CGAL_triangulation_assertion(vert[i]->point()[1] >= domain().ymin());
        CGAL_triangulation_assertion(vert[i]->point()[2] >= domain().zmin());

        // redirect also the cell pointer of the vertex.
        it->vertex(i)->set_cell(it);
      }
      // Set the offsets.
      set_offsets(it, off[0], off[1], off[2], off[3] );
      CGAL_triangulation_assertion( int_to_off(it->offset(0)) == off[0] );
      CGAL_triangulation_assertion( int_to_off(it->offset(1)) == off[1] );
      CGAL_triangulation_assertion( int_to_off(it->offset(2)) == off[2] );
      CGAL_triangulation_assertion( int_to_off(it->offset(3)) == off[3] );
    }
  }

  // ###################################################################
  // ### Fourth cell iteration #########################################
  // ###################################################################
  std::vector<Point> hidden_points;
  {
    // Delete the marked cells.
    std::vector<Cell_handle> cells_to_delete;
    for( Cell_iterator cit = all_cells_begin();
                        cit != all_cells_end(); ++cit ) {
      if( cit->get_additional_flag() == 1 )
      {
        gather_cell_hidden_points(cit, hidden_points);
        cells_to_delete.push_back( cit );
      }
    }
    _tds.delete_cells(cells_to_delete.begin(), cells_to_delete.end());
  }

  // ###################################################################
  // ### Vertex iteration ##############################################
  // ###################################################################
  {
    // Delete all the vertices in virtual_vertices, that is all vertices
    // outside the original domain.
    std::vector<Vertex_handle> vertices_to_delete;
    for( Vertex_iterator vit = all_vertices_begin();
          vit != all_vertices_end(); ++vit ) {
      if( virtual_vertices.count( vit ) != 0 ) {
        CGAL_triangulation_assertion( virtual_vertices.count( vit ) == 1 );
        vertices_to_delete.push_back( vit );
      }
    }
    _tds.delete_vertices(vertices_to_delete.begin(), vertices_to_delete.end());
  }
  _cover = CGAL::make_array(1,1,1);
  virtual_vertices.clear();
  virtual_vertices_reverse.clear();
  reinsert_hidden_points_after_converting_to_1_sheeted(hidden_points);
}

template < class GT, class TDS >
inline void
Periodic_3_triangulation_3<GT,TDS>::convert_to_27_sheeted_covering()
{
  if(_cover == CGAL::make_array(3,3,3))
    return;

  CGAL_triangulation_precondition(is_1_cover());

  // Create 27 copies of each vertex and write virtual_vertices and
  // virtual_vertices_reverse
  std::list<Vertex_handle> original_vertices;
  // try to use std::copy instead of the following loop.
  for(Vertex_iterator vit = vertices_begin(); vit != vertices_end(); ++vit)
    original_vertices.push_back(vit);
  for(typename std::list<Vertex_handle>::iterator vit
       = original_vertices.begin(); vit != original_vertices.end(); ++vit) {
    Vertex_handle v_cp;
    std::vector<Vertex_handle> copies;
    for(int i=0; i<3; i++)
      for(int j=0; j<3; j++)
        for(int k=0; k<3; k++) {
          if(i==0 && j==0 && k==0)
            continue;
          v_cp = _tds.create_vertex(*vit);
          copies.push_back(v_cp);
          virtual_vertices.insert(std::make_pair(v_cp,
                                                 std::make_pair(*vit,Offset(i,j,k))));
        }
    virtual_vertices_reverse.insert(std::make_pair(*vit,copies));
  }

  // Create 27 copies of each cell from the respective copies of the
  // vertices and write virtual_cells and virtual_cells_reverse.
  typedef std::map<Cell_handle, std::vector<Cell_handle > >
      Virtual_cell_reverse_map;
  typedef typename Virtual_cell_reverse_map::const_iterator VCRMIT;

  Virtual_cell_reverse_map virtual_cells_reverse;

  std::list<Cell_handle> original_cells;
  for(Cell_iterator cit = cells_begin(); cit != cells_end(); ++cit)
    original_cells.push_back(cit);

  // Store vertex offsets in a separate data structure
  std::list< Offset > off_v;
  for(typename std::list<Vertex_handle>::iterator vit
       = original_vertices.begin(); vit != original_vertices.end(); ++vit) {
    Cell_handle ccc = (*vit)->cell();
    int v_index = ccc->index(*vit);
    off_v.push_back(int_to_off(ccc->offset(v_index)));
  }

  // Store neighboring offsets in a separate data structure
  std::list<std::array<Offset,4> > off_nb;
  for(typename std::list<Cell_handle>::iterator cit = original_cells.begin();
       cit != original_cells.end(); ++cit) {
    std::array<Offset,4> off_nb_c;
    for(int i=0; i<4; i++) {
      Cell_handle ccc = *cit;
      Cell_handle nnn = ccc->neighbor(i);
      off_nb_c[i] = neighbor_offset(ccc,i,nnn);
    }
    off_nb.push_back(off_nb_c);
  }

  // Create copies of cells
  for(typename std::list<Cell_handle>::iterator cit = original_cells.begin();
      cit != original_cells.end(); ++cit) {
    Cell_handle c_cp;
    std::vector<Cell_handle> copies;
    Virtual_vertex_reverse_map_it vvrmit[4];
    Offset vvoff[4];
    for(int i=0; i<4; i++) {
      vvrmit[i] = virtual_vertices_reverse.find((*cit)->vertex(i));
      CGAL_triangulation_assertion(
            vvrmit[i] != virtual_vertices_reverse.end());
      vvoff[i] = int_to_off((*cit)->offset(i));
    }
    Vertex_handle vvh[4];
    for(int n=0; n<26; n++) {
      for(int i=0; i<4; i++) {
        // Decomposition of n into an offset (nx,ny,nz):
        // nx = (n+1)/9, ny = ((n+1)/3)%3, nz = (n+1)%3
        int o_i = ((n+1)/9+vvoff[i].x()+3)%3;
        int o_j = ((n+1)/3+vvoff[i].y()+3)%3;
        int o_k = ((n+1)+vvoff[i].z()+3)%3;
        int n_c = 9*o_i+3*o_j+o_k-1;
        CGAL_triangulation_assertion(n_c >= -1);
        if(n_c == -1) vvh[i] = (*cit)->vertex(i);
        else           vvh[i] = vvrmit[i]->second[n_c];
      }
      c_cp = _tds.create_cell(vvh[0], vvh[1], vvh[2], vvh[3]);
      copies.push_back(c_cp);
    }
    virtual_cells_reverse.insert(std::make_pair(*cit,copies));
  }

  // Set new vertices of boundary cells of the original domain.
  for(typename std::list<Cell_handle>::iterator cit = original_cells.begin();
       cit != original_cells.end(); ++cit) {
    for(int i=0; i<4; i++) {
      Virtual_vertex_reverse_map_it vvrmit
          = virtual_vertices_reverse.find((*cit)->vertex(i));
      CGAL_triangulation_assertion(vvrmit != virtual_vertices_reverse.end());
      Offset vvoff = int_to_off((*cit)->offset(i));
      if(!vvoff.is_null()) {
        int n_c = 9*vvoff.x()+3*vvoff.y()+vvoff.z()-1;
        CGAL_triangulation_assertion(n_c >= 0);
        CGAL_triangulation_assertion(static_cast<unsigned int>(n_c)
                                     < vvrmit->second.size());
        (*cit)->set_vertex(i,vvrmit->second[n_c]);
      }
    }
  }

  // Set neighboring relations of cell copies
  typename std::list<std::array<Offset,4> >::iterator oit = off_nb.begin();
  for(typename std::list<Cell_handle>::iterator cit = original_cells.begin();
       cit != original_cells.end(); ++cit, ++oit) {
    CGAL_triangulation_assertion( oit != off_nb.end() );
    VCRMIT c_cp = virtual_cells_reverse.find(*cit);
    CGAL_triangulation_assertion(c_cp != virtual_cells_reverse.end());
    for(int i=0; i<4; i++) {
      Cell_handle cit_nb = (*cit)->neighbor(i);
      VCRMIT c_cp_nb = virtual_cells_reverse.find(cit_nb);
      CGAL_triangulation_assertion(c_cp_nb != virtual_cells_reverse.end());
      Offset nboff = (*oit)[i];
      for(int n=0; n<26; n++) {
        int n_nb;
        if(nboff.is_null()) n_nb = n;
        else {
          int o_i = ((n+1)/9-nboff.x()+3)%3;
          int o_j = ((n+1)/3-nboff.y()+3)%3;
          int o_k = (n+1-nboff.z()+3)%3;
          n_nb = 9*o_i+3*o_j+o_k-1;
        }
        if(n_nb == -1) {
          CGAL_triangulation_assertion(cit_nb->has_vertex(
                                         c_cp->second[n]->vertex((i+1)%4)) );
          CGAL_triangulation_assertion(cit_nb->has_vertex(
                                         c_cp->second[n]->vertex((i+2)%4)) );
          CGAL_triangulation_assertion(cit_nb->has_vertex(
                                         c_cp->second[n]->vertex((i+3)%4)) );
          c_cp->second[n]->set_neighbor(i,cit_nb);
        }
        else {
          CGAL_triangulation_assertion(n_nb >= 0);
          CGAL_triangulation_assertion(static_cast<unsigned int>(n_nb)
                                       <= c_cp_nb->second.size());
          CGAL_triangulation_assertion(c_cp_nb->second[n_nb]
                                       ->has_vertex(c_cp->second[n]->vertex((i+1)%4)) );
          CGAL_triangulation_assertion(c_cp_nb->second[n_nb]
                                       ->has_vertex(c_cp->second[n]->vertex((i+2)%4)) );
          CGAL_triangulation_assertion(c_cp_nb->second[n_nb]
                                       ->has_vertex(c_cp->second[n]->vertex((i+3)%4)) );
          c_cp->second[n]->set_neighbor(i,c_cp_nb->second[n_nb]);
        }
      }
    }
  }

  // Set neighboring relations of original cells
  oit = off_nb.begin();
  for(typename std::list<Cell_handle>::iterator cit = original_cells.begin();
       cit != original_cells.end(); ++cit, ++oit) {
    CGAL_triangulation_assertion( oit != off_nb.end() );
    for(int i=0; i<4; i++) {
      Offset nboff = (*oit)[i];
      if(!nboff.is_null()) {
        Cell_handle cit_nb = (*cit)->neighbor(i);
        VCRMIT c_cp_nb = virtual_cells_reverse.find(cit_nb);
        CGAL_triangulation_assertion(c_cp_nb != virtual_cells_reverse.end());
        int o_i = (3-nboff.x())%3;
        int o_j = (3-nboff.y())%3;
        int o_k = (3-nboff.z())%3;
        int n_nb = 9*o_i+3*o_j+o_k-1;
        CGAL_triangulation_assertion(n_nb >= 0);
        CGAL_triangulation_assertion(static_cast<unsigned int>(n_nb)
                                     <= c_cp_nb->second.size());
        CGAL_triangulation_assertion(c_cp_nb->second[n_nb]
                                     ->has_vertex((*cit)->vertex((i+1)%4)) );
        CGAL_triangulation_assertion(c_cp_nb->second[n_nb]
                                     ->has_vertex((*cit)->vertex((i+2)%4)) );
        CGAL_triangulation_assertion(c_cp_nb->second[n_nb]
                                     ->has_vertex((*cit)->vertex((i+3)%4)) );
        (*cit)->set_neighbor(i,c_cp_nb->second[n_nb]);
      }
    }
  }

  // Set incident cells
  for(Cell_iterator cit = cells_begin(); cit != cells_end(); ++cit) {
    for(int i=0; i<4; i++) {
      cit->vertex(i)->set_cell(cit);
    }
  }

  // Set offsets where necessary
  for(typename std::list<Cell_handle>::iterator cit = original_cells.begin();
       cit != original_cells.end(); ++cit) {
    VCRMIT c_cp = virtual_cells_reverse.find(*cit);
    CGAL_triangulation_assertion( c_cp != virtual_cells_reverse.end());
    Offset off[4];
    for(int i=0; i<4; i++)
      off[i] = int_to_off((*cit)->offset(i));

    if(off[0].is_null() && off[1].is_null() &&
        off[2].is_null() && off[3].is_null())
      continue;

    for(int n=0; n<26; n++) {
      Offset off_cp[4];
      int o_i = (n+1)/9;
      int o_j = ((n+1)/3)%3;
      int o_k = (n+1)%3;
      if(o_i!=2 && o_j!=2 && o_k !=2)
        continue;

      for(int i=0; i<4; i++) {
        off_cp[i] = Offset((o_i==2)?off[i].x():0,
                           (o_j==2)?off[i].y():0,
                           (o_k==2)?off[i].z():0);
        CGAL_triangulation_assertion(off_cp[i].x() == 0 || off_cp[i].x() == 1);
        CGAL_triangulation_assertion(off_cp[i].y() == 0 || off_cp[i].y() == 1);
        CGAL_triangulation_assertion(off_cp[i].z() == 0 || off_cp[i].z() == 1);
      }
      set_offsets(c_cp->second[n],off_cp[0],off_cp[1],off_cp[2],off_cp[3]);
    }
  }

  // Iterate over all original cells and reset offsets.
  for(typename std::list<Cell_handle>::iterator cit = original_cells.begin();
       cit != original_cells.end(); ++cit) {
    //This statement does not seem to have any effect
    set_offsets(*cit, 0,0,0,0);
    CGAL_triangulation_assertion((*cit)->offset(0) == 0);
    CGAL_triangulation_assertion((*cit)->offset(1) == 0);
    CGAL_triangulation_assertion((*cit)->offset(2) == 0);
    CGAL_triangulation_assertion((*cit)->offset(3) == 0);
  }

  _cover = CGAL::make_array(3,3,3);
  CGAL_triangulation_expensive_assertion(is_valid());

  update_cover_data_after_converting_to_27_sheeted_covering();
}

template < class GT, class TDS >
class Periodic_3_triangulation_3<GT,TDS>::Finder
{
  const Self* _t;
  const Point& _p;

public:
  Finder(const Self* t, const Point &p) : _t(t), _p(p) {}
  bool operator()(const Vertex_handle v) {
    return _t->equal(v->point(), _p);
  }
};

/** Find the cell that consists of the four given vertices
 *
 *  Iterates over all cells and compare the four vertices of each cell
 *  with the four vertices in vh.
 */
template < class GT, class TDS >
inline typename Periodic_3_triangulation_3<GT,TDS>::Cell_handle
Periodic_3_triangulation_3<GT,TDS>::get_cell(const Vertex_handle* vh) const
{
  bool contains_v[4];
  std::vector<Cell_handle> cells;
  incident_cells(vh[3],std::back_inserter(cells));
  for( typename std::vector<Cell_handle>::iterator it = cells.begin();
       it != cells.end(); it++ ) {
    CGAL_triangulation_assertion(
          (*it)->vertex(0) == vh[3] || (*it)->vertex(1) == vh[3]
        ||(*it)->vertex(2) == vh[3] || (*it)->vertex(3) == vh[3]);
    for( int j=0; j<3; j++ ) {
      contains_v[j] = false;
      contains_v[j] = ( (*it)->vertex(0) == vh[j] )
          || ( (*it)->vertex(1) == vh[j] )
          || ( (*it)->vertex(2) == vh[j] )
          || ( (*it)->vertex(3) == vh[j] );
    }
    if(contains_v[0] && contains_v[1] && contains_v[2]) {
      return (*it);
    }
  }
  CGAL_triangulation_assertion(false);
  return Cell_handle();
}

/*! \brief Get the offset of tester.point() such that
 * this point is in conflict with c w.r.t tester.get_offset().
 *
 * Implementation: Just try all eight possibilities.
 */
template < class GT, class TDS >
template < class Conflict_tester >
inline typename Periodic_3_triangulation_3<GT,TDS>::Offset
Periodic_3_triangulation_3<GT,TDS>::get_location_offset(
    const Conflict_tester& tester, Cell_handle c) const
{
  CGAL_triangulation_precondition( number_of_vertices() != 0 );

  int cumm_off = c->offset(0) | c->offset(1) | c->offset(2) | c->offset(3);
  if(cumm_off == 0) {
    // default case:
    CGAL_triangulation_assertion(tester(c, Offset()));
    return Offset();
  } else {
    // Main idea seems to just test all possibilities.
    for(int i=0; i<8; i++) {
      if(((cumm_off | (~i))&7) == 7) {
        if(tester(c,int_to_off(i))) {
          return int_to_off(i);
        }
      }
    }
  }
  CGAL_triangulation_assertion(false);
  return Offset();
}

template < class GT, class TDS >
template < class Conflict_tester >
inline typename Periodic_3_triangulation_3<GT,TDS>::Offset
Periodic_3_triangulation_3<GT,TDS>::get_location_offset(
    const Conflict_tester& tester, Cell_handle c, bool& found) const
{
  CGAL_triangulation_precondition( number_of_vertices() != 0 );

  found = false;

  int cumm_off = c->offset(0) | c->offset(1) | c->offset(2) | c->offset(3);
  if(cumm_off == 0 && tester(c, Offset())) {
    // default case:
    found = true;
    return Offset();
  } else {
    // Main idea seems to just test all possibilities.
    for(int i=0; i<8; i++) {
      if(((cumm_off | (~i))&7) == 7) {
        if(tester(c,int_to_off(i))) {
          found = true;
          return int_to_off(i);
        }
      }
    }
  }

  return Offset();
}

/** Get the offset between the origins of the internal offset coordinate
  * systems of two neighboring cells with respect from ch to nb.
  *
  * - Find two corresponding vertices from each cell
  * - Return the difference of their offsets.
  */
template < class GT, class TDS >
inline typename Periodic_3_triangulation_3<GT,TDS>::Offset
Periodic_3_triangulation_3<GT,TDS>::neighbor_offset(
    Cell_handle ch, int i, Cell_handle nb) const
{
  // Redundance in the signature!
  CGAL_triangulation_precondition(ch->neighbor(i) == nb);
  CGAL_triangulation_precondition(nb->neighbor(nb->index(ch)) == ch);

  Vertex_handle vertex_ch;
  int index_ch, index_nb;
  // ensure that vertex_ch \in nb and vertex_nb \in ch
  index_ch = (i==0? 1 : 0);
  vertex_ch = ch->vertex(index_ch);
  index_nb = nb->index(vertex_ch);

  return int_to_off(nb->offset(index_nb)) - int_to_off(ch->offset(index_ch));
}

/**
 * - ch->offset(i) is an bit triple encapsulated in an integer. Each bit
 *   represents the offset in one direction --> 2-cover!
 * - it_to_off(int) decodes this again.
 * - Finally the offset vector is multiplied by cover.
 *   So if we are working in 3-cover we translate it to the neighboring
 *   3-cover and not only to the neighboring domain.
 */
template < class GT, class TDS >
inline void Periodic_3_triangulation_3<GT, TDS>::get_vertex(
    Cell_handle ch, int i, Vertex_handle& vh, Offset& off) const
{
  off = combine_offsets(Offset(),int_to_off(ch->offset(i)));
  vh = ch->vertex(i);

  if(is_1_cover())
    return;

  Vertex_handle vh_i = vh;
  get_vertex(vh_i, vh, off);
  return;
}

template < class GT, class TDS >
inline void Periodic_3_triangulation_3<GT, TDS>::get_vertex(
    Vertex_handle vh_i, Vertex_handle& vh, Offset& off) const
{
  Virtual_vertex_map_it it = virtual_vertices.find(vh_i);

  if(it == virtual_vertices.end()) {
    // if 'vh_i' is not contained in virtual_vertices, then it is in
    // the original domain.
    vh = vh_i;
    CGAL_triangulation_assertion(vh != Vertex_handle());
  } else {
    // otherwise it has to be looked up as well as its offset.
    vh = it->second.first;
    off += it->second.second;
    CGAL_triangulation_assertion(vh->point().x() < domain().xmax());
    CGAL_triangulation_assertion(vh->point().y() < domain().ymax());
    CGAL_triangulation_assertion(vh->point().z() < domain().zmax());
    CGAL_triangulation_assertion(vh->point().x() >= domain().xmin());
    CGAL_triangulation_assertion(vh->point().y() >= domain().ymin());
    CGAL_triangulation_assertion(vh->point().z() >= domain().zmin());
  }
}

template < class GT, class TDS >
std::istream&
operator>> (std::istream& is, Periodic_3_triangulation_3<GT,TDS>& tr)
  // reads
  // the current covering that guarantees the triangulation to be a
  //     simplicial complex
  // the number of vertices
  // the non combinatorial information on vertices (points in case of 1-sheeted
  //     covering, point-offset pairs otherwise)
  //     ALL PERIODIC COPIES OF ONE VERTEX MUST BE STORED CONSECUTIVELY
  // the number of cells
  // the cells by the indices of their vertices in the preceding list
  // of vertices, plus the non combinatorial information on each cell
  // the neighbors of each cell by their index in the preceding list of cells
{
  CGAL_triangulation_precondition(is.good());

  typedef Periodic_3_triangulation_3<GT,TDS>            Triangulation;
  typedef typename Triangulation::size_type             size_type;
  typedef typename Triangulation::Vertex_handle         Vertex_handle;
  typedef typename Triangulation::Cell_handle           Cell_handle;
  typedef typename Triangulation::Offset                Offset;
  typedef typename Triangulation::Iso_cuboid            Iso_cuboid;

  tr.clear();

  Iso_cuboid domain(0,0,0,1,1,1);
  int cx=0, cy=0, cz=0;
  size_type n=0;

  if(is_ascii(is)) {
    is >> domain;
    is >> cx >> cy >> cz;
    is >> n;
  }
  else {
    is >> domain;
    read(is,cx);
    read(is,cy);
    read(is,cz);
    read(is,n);
  }

  CGAL_triangulation_assertion((n/(cx*cy*cz))*cx*cy*cz == n);

  tr.tds().set_dimension((n==0?-2:3));
  tr._gt.set_domain(domain);
  tr._cover = CGAL::make_array(cx,cy,cz);

  if( n==0 ) return is;

  std::vector< Vertex_handle > V(n);

  if(cx==1 && cy==1 && cz==1) {
    for(std::size_t i=0; i < n; i++) {
      V[i] = tr.tds().create_vertex();
      is >> *V[i];
    }
  } else {
    Vertex_handle v,w;
    std::vector<Vertex_handle> vv;
    Offset off;
    for(std::size_t i=0; i < n; i++) {
      v = tr.tds().create_vertex();
      V[i] = v;
      is >> *V[i] >> off;
      vv.clear();
      for(int j=1; j<cx*cy*cz; j++) {
        i++;
        w = tr.tds().create_vertex();
        V[i] = w;
        is >> *V[i] >> off;
        vv.push_back(w);
        tr.virtual_vertices[w]=std::make_pair(v,off);
      }
      tr.virtual_vertices_reverse[v]=vv;
    }
  }

  std::vector< Cell_handle > C;
  std::size_t m;
  tr.tds().read_cells(is, V, m, C);

  // read offsets
  int off[4] = {0,0,0,0};
  for(std::size_t j=0; j < m; j++) {
    if(is_ascii(is))
      is >> off[0] >> off[1] >> off[2] >> off[3];
    else {
      read(is,off[0]);
      read(is,off[1]);
      read(is,off[2]);
      read(is,off[3]);
    }
    tr.set_offsets(C[j],off[0],off[1],off[2],off[3]);
  }

  // read potential other information
  for(std::size_t j=0; j < m; j++)
    is >> *(C[j]);

  CGAL_triangulation_expensive_assertion( tr.is_valid() );
  return is;
}

template < class GT, class TDS >
std::ostream&
operator<< (std::ostream& os,const Periodic_3_triangulation_3<GT,TDS>& tr)
// writes :
// the number of vertices
// the domain as six coordinates: xmin ymin zmin xmax ymax zmax
// the current covering that guarantees the triangulation to be a
//     simplicial complex
// the non combinatorial information on vertices (points in case of 1-sheeted
//     covering, point-offset pairs otherwise)
//     ALL PERIODIC COPIES OF ONE VERTEX MUST BE STORED CONSECUTIVELY
// the number of cells
// the cells by the indices of their vertices in the preceding list
// of vertices, plus the non combinatorial information on each cell
// the neighbors of each cell by their index in the preceding list of cells
{
  typedef Periodic_3_triangulation_3<GT,TDS>       Triangulation;
  typedef typename Triangulation::size_type        size_type;
  typedef typename Triangulation::Vertex_handle    Vertex_handle;
  typedef typename Triangulation::Vertex_iterator  Vertex_iterator;
  typedef typename Triangulation::Cell_handle      Cell_handle;
  typedef typename Triangulation::Cell_iterator    Cell_iterator;
  typedef typename Triangulation::Covering_sheets  Covering_sheets;
  typedef typename Triangulation::Offset           Offset;
  typedef typename Triangulation::Virtual_vertex_map_it Virtual_vertex_map_it;
  typedef typename Triangulation::Iso_cuboid       Iso_cuboid;

  // outputs dimension, domain and number of vertices
  Iso_cuboid domain = tr.domain();
  Covering_sheets cover = tr.number_of_sheets();
  size_type n = tr.number_of_vertices();

  if(is_ascii(os))
    os << domain << std::endl
       << cover[0] << " " << cover[1] << " " << cover[2] << std::endl
       << n*cover[0]*cover[1]*cover[2] << std::endl;
  else {
    os << domain;
    write(os,cover[0]);
    write(os,cover[1]);
    write(os,cover[2]);
    write(os,n*cover[0]*cover[1]*cover[2]);
  }

  if(n == 0)
    return os;

  // write the vertices
  Unique_hash_map<Vertex_handle, std::size_t > V;
  std::size_t i=0;
  if(tr.is_1_cover()) {
    for(Vertex_iterator it=tr.vertices_begin(); it!=tr.vertices_end(); ++it) {
      V[it] = i++;
      os << it->point();
      if(is_ascii(os))
        os << std::endl;
    }
  } else {
    Virtual_vertex_map_it vit, vvit;
    std::vector<Vertex_handle> vv;
    for(Vertex_iterator it=tr.vertices_begin(); it!=tr.vertices_end(); ++it) {
      vit = tr.virtual_vertices.find(it);
      if(vit != tr.virtual_vertices.end())
        continue;

      V[it]=i++;
      if(is_ascii(os))
        os << it->point() << std::endl
           << Offset(0,0,0) << std::endl;
      else os << it->point() << Offset(0,0,0);
      CGAL_triangulation_assertion(tr.virtual_vertices_reverse.find(it)
          != tr.virtual_vertices_reverse.end());
      vv = tr.virtual_vertices_reverse.find(it)->second;
      CGAL_triangulation_assertion(vv.size() == 26);
      for(std::size_t j=0; j<vv.size(); j++) {
        vvit = tr.virtual_vertices.find(vv[j]);
        CGAL_triangulation_assertion(vvit != tr.virtual_vertices.end());
        V[vv[j]] = i++;
        if(is_ascii(os))
          os << vv[j]->point() << std::endl
             << vvit->second.second << std::endl;
        else os << vv[j]->point() << vvit->second.second;
      }
    }
  }
  CGAL_triangulation_postcondition(i == tr.number_of_sheets()[0] *
                                        tr.number_of_sheets()[1] *
                                        tr.number_of_sheets()[2] * n);

  // asks the tds for the combinatorial information
  tr.tds().print_cells(os, V);

  // write offsets
  //for(unsigned int i=0; i<tr.number_of_cells(); i++) {
  for(Cell_iterator it=tr.cells_begin(); it!=tr.cells_end(); ++it) {
    //Cell_handle ch = std::find(tr.cells_begin(), tr.cells_end(), i);
    Cell_handle ch(it);
    for(int j=0; j<4; j++) {
      if(is_ascii(os)) {
        os << ch->offset(j);
        if( j==3 )
          os << std::endl;
        else
          os << ' ';
      }
      else write(os,ch->offset(j));
    }
  }

  // write the non combinatorial information on the cells
  // using the << operator of Cell
  // works because the iterator of the tds traverses the cells in the
  // same order as the iterator of the triangulation
  if(tr.number_of_vertices() != 0) {
      for(Cell_iterator it=tr.cells_begin(); it != tr.cells_end(); ++it) {
    os << *it; // other information
    if(is_ascii(os))
      os << std::endl;
    }
  }
  return os;
}

namespace internal {

/// Internal function used by operator==.
// This function tests and registers the 4 neighbors of c1/c2,
// and performs a bfs traversal
// Returns false if an inequality has been found.
//TODO: introduce offsets
template <class GT, class TDS1, class TDS2>
bool
test_next(const Periodic_3_triangulation_3<GT, TDS1>& t1,
          const Periodic_3_triangulation_3<GT, TDS2>& t2,
          typename Periodic_3_triangulation_3<GT, TDS1>::Cell_handle c1,
          typename Periodic_3_triangulation_3<GT, TDS2>::Cell_handle c2,
          std::map<typename Periodic_3_triangulation_3<GT, TDS1>::Cell_handle,
          typename Periodic_3_triangulation_3<GT, TDS2>::Cell_handle>& Cmap,
          std::map<typename Periodic_3_triangulation_3<GT, TDS1>::Vertex_handle,
          typename Periodic_3_triangulation_3<GT, TDS2>::Vertex_handle>& Vmap)
{
  typedef Periodic_3_triangulation_3<GT, TDS1> Tr1;
  typedef Periodic_3_triangulation_3<GT, TDS2> Tr2;
  typedef typename Tr1::Vertex_handle  Vertex_handle1;
  typedef typename Tr1::Cell_handle    Cell_handle1;
  typedef typename Tr2::Vertex_handle  Vertex_handle2;
  typedef typename Tr2::Cell_handle    Cell_handle2;
  typedef typename std::map<Cell_handle1, Cell_handle2>::const_iterator  Cit;
  typedef typename std::map<Vertex_handle1, Vertex_handle2>::const_iterator Vit;

  std::vector<std::pair<Cell_handle1, Cell_handle2> > queue;
  queue.push_back(std::make_pair(c1,c2));

  while(! queue.empty()) {
    boost::tie(c1,c2) = queue.back();
    queue.pop_back();

    // Precondition: c1, c2 have been registered as well as their 4 vertices.
    CGAL_triangulation_precondition(t1.number_of_vertices() != 0);
    CGAL_triangulation_precondition(Cmap[c1] == c2);
    CGAL_triangulation_precondition(Vmap.find(c1->vertex(0)) != Vmap.end());
    CGAL_triangulation_precondition(Vmap.find(c1->vertex(1)) != Vmap.end());
    CGAL_triangulation_precondition(Vmap.find(c1->vertex(2)) != Vmap.end());
    CGAL_triangulation_precondition(Vmap.find(c1->vertex(3)) != Vmap.end());

    for(int i=0; i <= 3; ++i) {
      Cell_handle1 n1 = c1->neighbor(i);
      Cit cit = Cmap.find(n1);
      Vertex_handle1 v1 = c1->vertex(i);
      Vertex_handle2 v2 = Vmap[v1];
      Cell_handle2 n2 = c2->neighbor(c2->index(v2));
      if(cit != Cmap.end()) {
        // n1 was already registered.
        if(cit->second != n2)
          return false;
        continue;
      }
      // n1 has not yet been registered.
      // We check that the new vertices match geometrically.
      // And we register them.
      Vertex_handle1 vn1 = n1->vertex(n1->index(c1));
      Vertex_handle2 vn2 = n2->vertex(n2->index(c2));
      Vit vit = Vmap.find(vn1);
      if(vit != Vmap.end()) {
        // vn1 already registered
        if(vit->second != vn2)
          return false;
      }
      else {
        if(t1.geom_traits().compare_xyz_3_object()(t1.construct_point(vn1->point()),
                                                    t2.construct_point(vn2->point())) != 0)
          return false;

        // We register vn1/vn2.
        Vmap.insert(std::make_pair(vn1, vn2));
      }

      // We register n1/n2.
      Cmap.insert(std::make_pair(n1, n2));
      queue.push_back(std::make_pair(n1, n2));
    }
  }
  return true;
}

} // namespace internal


template < class GT, class TDS1, class TDS2  >
bool
operator==(const Periodic_3_triangulation_3<GT,TDS1>& t1,
           const Periodic_3_triangulation_3<GT,TDS2>& t2)
{
  typedef typename Periodic_3_triangulation_3<GT,TDS1>::Vertex_handle
      Vertex_handle1;
  typedef typename Periodic_3_triangulation_3<GT,TDS1>::Cell_handle
      Cell_handle1;
  typedef typename Periodic_3_triangulation_3<GT,TDS2>::Vertex_handle
      Vertex_handle2;
  typedef typename Periodic_3_triangulation_3<GT,TDS2>::Vertex_handle
      Vertex_iterator2;
  typedef typename Periodic_3_triangulation_3<GT,TDS2>::Cell_handle
      Cell_handle2;

  typedef typename Periodic_3_triangulation_3<GT,TDS1>::Point      Point;
  typedef typename Periodic_3_triangulation_3<GT,TDS1>::Offset     Offset;

  // typedef typename Periodic_3_triangulation_3<GT,TDS1>
  //     ::Geometric_traits::Compare_xyz_3                       Compare_xyz_3;
  // Compare_xyz_3 cmp1 = t1.geom_traits().compare_xyz_3_object();
  // Compare_xyz_3 cmp2 = t2.geom_traits().compare_xyz_3_object();

  // Some quick checks.
  if(   t1.domain()           != t2.domain()
      || t1.number_of_sheets() != t2.number_of_sheets())
    return false;

  if(   t1.number_of_vertices() != t2.number_of_vertices()
      || t1.number_of_cells() != t2.number_of_cells())
    return false;

  // Special case for empty triangulations
  if(t1.number_of_vertices() == 0)
    return true;

  // We will store the mapping between the 2 triangulations vertices and
  // cells in 2 maps.
  std::map<Vertex_handle1, Vertex_handle2> Vmap;
  std::map<Cell_handle1, Cell_handle2> Cmap;

  // find a common point
  Vertex_handle1 v1 = static_cast<Vertex_handle1>(t1.vertices_begin());
  Vertex_handle2 iv2;
  for(Vertex_iterator2 vit2 = t2.vertices_begin();
      vit2 != t2.vertices_end(); ++vit2) {
    if(!t1.equal(vit2->point(), v1->point(),
                  t2.get_offset(vit2), t1.get_offset(v1)))
      continue;
    iv2 = static_cast<Vertex_handle2>(vit2);
    break;
  }
  if(iv2 == Vertex_handle2())
    return false;
  Vmap.insert(std::make_pair(v1, iv2));

  // We pick one cell of t1, and try to match it against the
  // cells of t2.
  Cell_handle1 c = v1->cell();
  Vertex_handle1 v2 = c->vertex((c->index(v1)+1)%4);
  Vertex_handle1 v3 = c->vertex((c->index(v1)+2)%4);
  Vertex_handle1 v4 = c->vertex((c->index(v1)+3)%4);
  Point p2 = v2->point();
  Point p3 = v3->point();
  Point p4 = v4->point();
  Offset o2 = t1.get_offset(v2);
  Offset o3 = t1.get_offset(v3);
  Offset o4 = t1.get_offset(v4);

  std::vector<Cell_handle2> ics;
  t2.incident_cells(iv2, std::back_inserter(ics));
  for(typename std::vector<Cell_handle2>::const_iterator cit = ics.begin();
                                                          cit != ics.end(); ++cit) {
    int inf = (*cit)->index(iv2);

    if(t1.equal(p2, (*cit)->vertex((inf+1)%4)->point(),
                 o2, t2.get_offset((*cit)->vertex((inf+1)%4))))
      Vmap.insert(std::make_pair(v2, (*cit)->vertex((inf+1)%4)));
    else if(t1.equal(p2, (*cit)->vertex((inf+2)%4)->point(),
                      o2, t2.get_offset((*cit)->vertex((inf+2)%4))))
      Vmap.insert(std::make_pair(v2, (*cit)->vertex((inf+2)%4)));
    else if(t1.equal(p2, (*cit)->vertex((inf+3)%4)->point(),
                      o2, t2.get_offset((*cit)->vertex((inf+3)%4))))
      Vmap.insert(std::make_pair(v2, (*cit)->vertex((inf+3)%4)));
    else
      continue; // None matched v2.

    if(t1.equal(p3, (*cit)->vertex((inf+1)%4)->point(),
                 o3, t2.get_offset((*cit)->vertex((inf+1)%4))))
      Vmap.insert(std::make_pair(v3, (*cit)->vertex((inf+1)%4)));
    else if(t1.equal(p3, (*cit)->vertex((inf+2)%4)->point(),
                      o3, t2.get_offset((*cit)->vertex((inf+2)%4))))
      Vmap.insert(std::make_pair(v3, (*cit)->vertex((inf+2)%4)));
    else if(t1.equal(p3, (*cit)->vertex((inf+3)%4)->point(),
                      o3, t2.get_offset((*cit)->vertex((inf+3)%4))))
      Vmap.insert(std::make_pair(v3, (*cit)->vertex((inf+3)%4)));
    else
      continue; // None matched v3.

    if(t1.equal(p4, (*cit)->vertex((inf+1)%4)->point(),
                 o4, t2.get_offset((*cit)->vertex((inf+1)%4))))
      Vmap.insert(std::make_pair(v4,(*cit)->vertex((inf+1)%4)));
    else if(t1.equal(p4, (*cit)->vertex((inf+2)%4)->point(),
                      o4, t2.get_offset((*cit)->vertex((inf+2)%4))))
      Vmap.insert(std::make_pair(v4,(*cit)->vertex((inf+2)%4)));
    else if(t1.equal(p4, (*cit)->vertex((inf+3)%4)->point(),
                      o4, t2.get_offset((*cit)->vertex((inf+3)%4))))
      Vmap.insert(std::make_pair(v4,(*cit)->vertex((inf+3)%4)));
    else
      continue; // None matched v4.

    // Found it !
    Cmap.insert(std::make_pair(c, *cit));
    break;
  }

  if(Cmap.size() == 0)
    return false;

  // We now have one cell, we need to compare in a bfs graph traversal
  return internal::test_next(t1, t2, Cmap.begin()->first, Cmap.begin()->second, Cmap, Vmap);
}

template < class GT, class TDS1, class TDS2 >
inline
bool
operator!=(const Periodic_3_triangulation_3<GT,TDS1>& t1,
    const Periodic_3_triangulation_3<GT,TDS2>& t2)
{
  return ! (t1 == t2);
}

} // namespace CGAL

#endif // CGAL_PERIODIC_3_TRIANGULATION_3_H
