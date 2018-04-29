// Copyright (c) 1999-2004,2006-2009,2013-2015,2017  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Monique Teillaud <Monique.Teillaud@inria.fr>
//                 Aymeric Pelle <Aymeric.Pelle@sophia.inria.fr>
//                 Mael Rouxel-Labbé

#ifndef CGAL_PERIODIC_3_REGULAR_TRIANGULATION_3_H
#define CGAL_PERIODIC_3_REGULAR_TRIANGULATION_3_H

#include <CGAL/license/Periodic_3_triangulation_3.h>

// Needed by remove to fill the hole.
#include <CGAL/internal/Periodic_3_regular_triangulation_remove_traits_3.h>

#include <CGAL/Periodic_3_triangulation_3.h>
#include <CGAL/Regular_triangulation_vertex_base_3.h>
#include <CGAL/Regular_triangulation_cell_base_3.h>
#include <CGAL/internal/Has_nested_type_Bare_point.h>

#include <CGAL/enum.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/spatial_sort.h>
#include <CGAL/utility.h>
#include <CGAL/internal/boost/function_property_map.hpp>

#include <boost/mpl/if.hpp>
#include <boost/mpl/identity.hpp>
#include <boost/unordered_set.hpp>

#include <cstdlib>
#include <functional>
#include <iterator>
#include <list>
#include <map>
#include <vector>
#include <utility>

namespace CGAL
{
template < class Gt,
           class Tds =
             Triangulation_data_structure_3 <
               Regular_triangulation_vertex_base_3<Gt,
                 Periodic_3_triangulation_ds_vertex_base_3<> >,
               Regular_triangulation_cell_base_3<Gt,
                 Periodic_3_triangulation_ds_cell_base_3<> >
             >
         >
class Periodic_3_regular_triangulation_3
  : public Periodic_3_triangulation_3<Gt, Tds>
{
  typedef Periodic_3_regular_triangulation_3<Gt, Tds>             Self;

public:
  typedef Periodic_3_triangulation_3<Gt, Tds>      Tr_Base;

  typedef Gt                                       Geometric_traits;
  typedef Geometric_traits                         Geom_traits;
  typedef Tds                                      Triangulation_data_structure;

  typedef typename Gt::FT                          FT;

  typedef typename Tr_Base::Periodic_segment_3             Periodic_segment_3;
  typedef typename Tr_Base::Periodic_triangle_3            Periodic_triangle_3;
  typedef typename Tr_Base::Periodic_tetrahedron_3         Periodic_tetrahedron_3;
  typedef typename Tr_Base::Periodic_segment               Periodic_segment;
  typedef typename Tr_Base::Periodic_triangle              Periodic_triangle;
  typedef typename Tr_Base::Periodic_tetrahedron           Periodic_tetrahedron;
  typedef typename Tr_Base::Periodic_tetrahedron_iterator  Periodic_tetrahedron_iterator;

  typedef typename Tr_Base::Vertex_handle          Vertex_handle;
  typedef typename Tr_Base::Cell_handle            Cell_handle;

  typedef typename Tr_Base::Vertex                 Vertex;
  typedef typename Tr_Base::Edge                   Edge;
  typedef typename Tr_Base::Facet                  Facet;
  typedef typename Tr_Base::Cell                   Cell;

  typedef typename Tr_Base::Vertex_iterator        Vertex_iterator;
  typedef typename Tr_Base::Edge_iterator          Edge_iterator;
  typedef typename Tr_Base::Facet_iterator         Facet_iterator;
  typedef typename Tr_Base::Facet_circulator       Facet_circulator;
  typedef typename Tr_Base::Cell_iterator          Cell_iterator;
  typedef typename Tr_Base::Cell_circulator        Cell_circulator;

  typedef typename Tr_Base::All_vertices_iterator  All_vertices_iterator;
  typedef typename Tr_Base::All_edges_iterator     All_edges_iterator;
  typedef typename Tr_Base::All_facets_iterator    All_facets_iterator;
  typedef typename Tr_Base::All_cells_iterator     All_cells_iterator;

  typedef typename Tr_Base::size_type              size_type;
  typedef typename Tr_Base::difference_type        difference_type;
  typedef typename Tr_Base::Locate_type            Locate_type;

  typedef typename Tr_Base::Offset                 Offset;
  typedef typename Tr_Base::Iso_cuboid             Iso_cuboid;
  typedef typename Tr_Base::Covering_sheets        Covering_sheets;

  // Traits are not supposed to define "Bare_point", but checking for backward compatibility
  typedef typename boost::mpl::eval_if_c<
    internal::Has_nested_type_Bare_point<Gt>::value,
    typename internal::Bare_point_type<Gt>,
    boost::mpl::identity<typename Gt::Point_3>
  >::type                                          Bare_point;
  typedef typename Gt::Weighted_point_3            Weighted_point;

  typedef std::pair<Bare_point, Offset>            Periodic_bare_point;
  typedef std::pair<Weighted_point, Offset>        Periodic_weighted_point;

  typedef typename Gt::Segment_3                   Segment;
  typedef typename Gt::Triangle_3                  Triangle;
  typedef typename Gt::Tetrahedron_3               Tetrahedron;

  typedef typename Gt::Object_3                    Object;

  //Tag to distinguish Delaunay from Regular triangulations
  typedef Tag_true                                 Weighted_tag;

  // Tag to distinguish periodic triangulations from others
  typedef Tag_true                                 Periodic_tag;

#ifndef CGAL_CFG_USING_BASE_MEMBER_BUG_2
  using Tr_Base::cw;
  using Tr_Base::ccw;
  using Tr_Base::domain;
  using Tr_Base::geom_traits;
  using Tr_Base::int_to_off;
  using Tr_Base::is_1_cover;
  using Tr_Base::is_virtual;
  using Tr_Base::mirror_vertex;
  using Tr_Base::number_of_vertices;
  using Tr_Base::number_of_edges;
  using Tr_Base::number_of_facets;
  using Tr_Base::number_of_cells;
  using Tr_Base::number_of_sheets;
  using Tr_Base::next_around_edge;
  using Tr_Base::orientation;
  using Tr_Base::point;
  using Tr_Base::set_offsets;
  using Tr_Base::swap;
  using Tr_Base::tds;
  using Tr_Base::vertex_triple_index;
  using Tr_Base::vertices_begin;
  using Tr_Base::vertices_end;
  using Tr_Base::edges_begin;
  using Tr_Base::edges_end;
  using Tr_Base::facets_begin;
  using Tr_Base::facets_end;
  using Tr_Base::cells_begin;
  using Tr_Base::cells_end;
#endif

  // For strict-ansi compliance
  using Tr_Base::combine_offsets;
  using Tr_Base::draw_dual;
  using Tr_Base::find_conflicts;
  using Tr_Base::get_offset;
  using Tr_Base::get_original_vertex;
  using Tr_Base::get_location_offset;
  using Tr_Base::adjacent_vertices;
  using Tr_Base::incident_edges;
  using Tr_Base::incident_facets;
  using Tr_Base::incident_cells;
  using Tr_Base::is_valid_conflict;
  using Tr_Base::locate;
  using Tr_Base::neighbor_offset;
  using Tr_Base::periodic_point;
#ifndef CGAL_NO_STRUCTURAL_FILTERING
  using Tr_Base::inexact_locate;
#endif

private:
  struct Cell_handle_hash
    : public CGAL::unary_function<Cell_handle, std::size_t>
  {
    std::size_t operator()(const Cell_handle& ch) const {
      return boost::hash<typename Cell_handle::pointer>()(&*ch);
    }
  };

  boost::unordered_set<Cell_handle, Cell_handle_hash> cells_with_too_big_orthoball;

  class Cover_manager
  {
    Periodic_3_regular_triangulation_3& tr;

public:
    Cover_manager (Periodic_3_regular_triangulation_3& tr)
      : tr(tr)
    { }

    void create_initial_triangulation() {
      tr.create_initial_triangulation();
    }

    template <class CellIt>
    void delete_unsatisfying_elements(const CellIt begin, const CellIt end) {
      tr.delete_cells_with_too_big_orthoball(begin, end);
    }

    template <class CellIt>
    void insert_unsatisfying_elements(Vertex_handle v, const CellIt begin, const CellIt end) {
      tr.insert_cells_with_too_big_orthoball(v, begin, end);
    }

    bool can_be_converted_to_1_sheet () const {
      return tr.can_be_converted_to_1_sheet();
    }

    bool update_cover_data_during_management (Cell_handle new_ch, const std::vector<Cell_handle>& new_cells) {
      return tr.update_cover_data_during_management(new_ch, new_cells);
    }
  };

public:
  /** @name Creation */ //@{
  Periodic_3_regular_triangulation_3 (const Iso_cuboid& domain = Iso_cuboid(0, 0, 0, 1, 1, 1),
                                      const Geometric_traits& gt = Geometric_traits())
    : Tr_Base(domain, gt)
  { }

  // copy constructor duplicates vertices and cells
  Periodic_3_regular_triangulation_3 (const Periodic_3_regular_triangulation_3& tr)
    : Tr_Base(tr)
  {
    if(is_1_cover())
      tds() = tr.tds();
    else
      this->copy_multiple_covering(tr);

    CGAL_triangulation_expensive_postcondition(*this == tr);
    CGAL_triangulation_expensive_postcondition( is_valid() );
  }

  template < typename InputIterator >
  Periodic_3_regular_triangulation_3(InputIterator first, InputIterator last,
                                     const Iso_cuboid& domain = Iso_cuboid(0,0,0,1,1,1),
                                     const Geometric_traits& gt = Geometric_traits(),
                                     bool is_large_point_set = false)
    : Tr_Base(domain, gt)
  {
    insert(first, last, is_large_point_set);
  }

  void copy_multiple_covering(const Periodic_3_regular_triangulation_3& tr)
  {
    // Write the respective offsets in the vertices to make them
    // automatically copy with the tds.
    for(Vertex_iterator vit = tr.vertices_begin(); vit != tr.vertices_end(); ++vit) {
      vit->set_offset(tr.get_offset(vit));
    }

    // copy the tds
    tds() = tr.tds();

    // make a list of all vertices that belong to the original
    // domain and initialize the basic structure of
    // virtual_vertices_reverse
    std::list<Vertex_handle> vlist;
    for(Vertex_iterator vit = vertices_begin(); vit != vertices_end(); ++vit) {
      if(vit->offset() == Offset()) {
        vlist.push_back(vit);
        this->virtual_vertices_reverse.insert(
          std::make_pair(vit,std::vector<Vertex_handle>(26)));
        CGAL_triangulation_assertion(this->virtual_vertices_reverse.find(vit)
                                       ->second.size() == 26);
      }
    }

    // Iterate over all vertices that are not in the original domain
    // and construct the respective entries to virtual_vertices and
    // virtual_vertices_reverse
    for(Vertex_iterator vit2 = vertices_begin();
         vit2 != vertices_end(); ++vit2) {
      if(vit2->offset() != Offset()) {
        typename std::list<Vertex_handle>::iterator vlist_it
            = std::find_if(vlist.begin(), vlist.end(),
                           typename Tr_Base::Finder(this,vit2->point()));
        Offset off = vit2->offset();
        this->virtual_vertices.insert(std::make_pair(vit2,
                                        std::make_pair(*vlist_it,off)));
        this->virtual_vertices_reverse.find(*vlist_it)
          ->second[9*off[0]+3*off[1]+off[2]-1]=vit2;
        CGAL_triangulation_assertion(get_offset(vit2) == off);
      }
    }

    // Cleanup vertex offsets
    for(Vertex_iterator vit = vertices_begin();
         vit != vertices_end(); ++vit)
      vit->clear_offset();
    for(Vertex_iterator vit = tr.vertices_begin();
         vit != tr.vertices_end(); ++vit)
      vit->clear_offset();

    insert_cells_with_too_big_orthoball(tr.cells_begin(), tr.cells_end());
  }

  Periodic_3_regular_triangulation_3 operator= (Periodic_3_regular_triangulation_3 tr)
  {
    tr.swap(*this);
    return *this;
  }

  void swap(Periodic_3_regular_triangulation_3&tr)
  {
    std::swap(cells_with_too_big_orthoball, tr.cells_with_too_big_orthoball);
    Tr_Base::swap(tr);
  }

  void create_initial_triangulation()
  {
    CGAL_triangulation_assertion( cells_with_too_big_orthoball.empty() );

    for(Cell_iterator iter = cells_begin(), end_iter = cells_end(); iter != end_iter; ++iter)
      cells_with_too_big_orthoball.insert(iter);
  }

  template <class CellIt>
  void delete_cells_with_too_big_orthoball(CellIt begin, const CellIt end)
  {
    for(; begin != end; ++begin)
    {
      typename boost::unordered_set<Cell_handle>::iterator iter = cells_with_too_big_orthoball.find(*begin);
      if(iter != cells_with_too_big_orthoball.end())
      {
        cells_with_too_big_orthoball.erase(iter);
      }
    }
  }

  CGAL::Comparison_result compare_orthsphere_radius_to_threshold (
      const Periodic_weighted_point& p0, const Periodic_weighted_point& p1,
      const Periodic_weighted_point& p2, const Periodic_weighted_point& p3,
      const FT threshold) const
  {
    return geom_traits().compare_weighted_squared_radius_3_object()(
             p0.first,  p1.first,  p2.first,  p3.first,
             p0.second, p1.second, p2.second, p3.second,
             threshold);
  }

  CGAL::Comparison_result compare_orthsphere_radius_to_threshold (Cell_handle cell,
                                                                  const FT threshold) const
  {
    Periodic_weighted_point p0 = periodic_point(cell, 0);
    Periodic_weighted_point p1 = periodic_point(cell, 1);
    Periodic_weighted_point p2 = periodic_point(cell, 2);
    Periodic_weighted_point p3 = periodic_point(cell, 3);

    return compare_orthsphere_radius_to_threshold(p0, p1, p2, p3, threshold);
  }

  template <class CellIt>
  void insert_cells_with_too_big_orthoball(Vertex_handle /*v*/, CellIt begin, const CellIt end)
  {
    FT threshold = FT(0.015625) * (domain().xmax()-domain().xmin()) * (domain().xmax()-domain().xmin());
    for(; begin != end; ++begin)
    {
      if(compare_orthsphere_radius_to_threshold(*begin, threshold) != CGAL::SMALLER)
      {
        cells_with_too_big_orthoball.insert(*begin);
      }
    }
  }

  void insert_cells_with_too_big_orthoball(Cell_iterator begin, Cell_iterator end)
  {
    FT threshold = FT(0.015625) * (domain().xmax()-domain().xmin()) * (domain().xmax()-domain().xmin());
    for(; begin != end; ++begin)
    {
      if(compare_orthsphere_radius_to_threshold(begin, threshold) != CGAL::SMALLER)
      {
        cells_with_too_big_orthoball.insert(begin);
      }
    }
  }

  bool can_be_converted_to_1_sheet () const
  {
    return cells_with_too_big_orthoball.empty();
  }

  bool update_cover_data_during_management (Cell_handle new_ch, const std::vector<Cell_handle>& new_cells)
  {
    bool result = false;
    FT threshold = FT(0.015625) * (domain().xmax() - domain().xmin()) * (domain().xmax() - domain().xmin());

    if(compare_orthsphere_radius_to_threshold(new_ch, threshold) != CGAL::SMALLER)
    {
      if(is_1_cover())
      {
        tds().delete_cells(new_cells.begin(), new_cells.end());
        this->convert_to_27_sheeted_covering();
        result = true;
      }
      else
      {
        cells_with_too_big_orthoball.insert(new_ch);
      }
    }

    return result;
  }

  virtual void update_cover_data_after_converting_to_27_sheeted_covering ()
  {
    FT threshold = FT(0.015625) * (domain().xmax()-domain().xmin()) * (domain().xmax() - domain().xmin());
    for(Cell_iterator iter = cells_begin(), end_iter = cells_end(); iter != end_iter; ++iter)
    {
      if(compare_orthsphere_radius_to_threshold(iter, threshold) != CGAL::SMALLER)
      {
        cells_with_too_big_orthoball.insert(iter);
      }
    }
  }

  // in the regular setting, there is nothing to do since we use `cells_with_too_big_orthoball`
  // instead of an edge length threshold
  virtual void update_cover_data_after_setting_domain () {}

  // the function below is used in `convert_to_1_sheeted_covering()` of P3T3
  // but only makes sense for regular triangulations.
  virtual void gather_cell_hidden_points(const Cell_handle cit,
                                         std::vector<Weighted_point>& hidden_points)
  {
    std::copy(cit->hidden_points_begin(), cit->hidden_points_end(),
              std::back_inserter(hidden_points));
  }

  virtual void reinsert_hidden_points_after_converting_to_1_sheeted(const std::vector<Weighted_point>& hidden_points)
  {
    typename std::vector<Weighted_point>::const_iterator wpv_it = hidden_points.begin();
    typename std::vector<Weighted_point>::const_iterator wpv_end = hidden_points.end();

    while(wpv_it != wpv_end)
      insert(*wpv_it++);
  }

  /** @name Insertion */ //@{
  Vertex_handle insert(const Weighted_point& point,
                       Cell_handle start = Cell_handle())
  {
    Conflict_tester tester(point, this);
    Point_hider hider(this);
    Cover_manager cover_manager(*this);
    CGAL_triangulation_precondition(point.weight() >= 0);
    CGAL_triangulation_precondition_msg
    (
      point.weight() < ( FT(0.015625) * (domain().xmax()-domain().xmin()) * (domain().xmax()-domain().xmin()) ),
     "point.weight() < 1/64 * domain_size * domain_size"
    );
    return Tr_Base::insert_in_conflict(point, start, tester, hider, cover_manager);
  }

  Vertex_handle insert(const Weighted_point& point,
                       Locate_type lt, Cell_handle c,
                       int li, int lj)
  {
    Conflict_tester tester(point, this);
    Point_hider hider(this);
    Cover_manager cover_manager(*this);
    CGAL_triangulation_precondition(point.weight() >= 0);
    CGAL_triangulation_precondition_msg
    (
      point.weight() < ( FT(0.015625) * (domain().xmax()-domain().xmin()) * (domain().xmax()-domain().xmin()) ),
      "point.weight() < 1/64 * domain_size * domain_size"
    );
    return Tr_Base::insert_in_conflict(point,lt,c,li,lj, tester,hider,cover_manager);
  }

  template < class InputIterator >
  std::ptrdiff_t insert(InputIterator first, InputIterator last,
                        bool is_large_point_set = false)
  {
    if(first == last)
      return 0;

    CGAL_triangulation_precondition_code
    (
      bool precondition_is_satisfied = true;
      FT upper_bound = FT(0.015625) * (domain().xmax()-domain().xmin()) * (domain().xmax()-domain().xmin());
      for(InputIterator pc_first = first, pc_last = last; pc_first != pc_last; ++pc_first)
      {
        if(pc_first->weight() < FT(0) || pc_first->weight() >= upper_bound)
        {
          precondition_is_satisfied = false;
          break;
        }
      }
    )

    CGAL_triangulation_precondition_msg
    (
      precondition_is_satisfied,
      "0 <= point.weight() < 1/64 * domain_size * domain_size"
    );

    size_type n = number_of_vertices();
    // The heuristic discards the existing triangulation so it can only be
    // applied to empty triangulations.
    if(n != 0)
      is_large_point_set = false;

    std::vector<Weighted_point> points(first, last);
    CGAL::cpp98::random_shuffle(points.begin(), points.end());
    Cell_handle hint;
    std::vector<Vertex_handle> dummy_points_vhs, double_vertices;
    std::vector<Weighted_point> dummy_points;
    typename std::vector<Weighted_point>::iterator pbegin = points.begin();
    if(is_large_point_set)
    {
      dummy_points_vhs = insert_dummy_points();
      dummy_points.reserve(dummy_points_vhs.size());
      for(typename std::vector<Vertex_handle>::iterator iter = dummy_points_vhs.begin(), end_iter = dummy_points_vhs.end(); iter != end_iter; ++iter)
        dummy_points.push_back((*iter)->point());
    }
    else
    {
      while(!is_1_cover())
      {
        insert(*pbegin);
        ++pbegin;
        if(pbegin == points.end())
          return number_of_vertices() - n;
      }
    }

    // Spatial sorting can only be applied to bare points, so we need an adaptor
    typedef typename Geom_traits::Construct_point_3 Construct_point_3;
    typedef typename boost::result_of<const Construct_point_3(const Weighted_point&)>::type Ret;
    typedef CGAL::internal::boost_::function_property_map<Construct_point_3, Weighted_point, Ret> fpmap;
    typedef CGAL::Spatial_sort_traits_adapter_3<Geom_traits, fpmap> Search_traits_3;

    spatial_sort(pbegin, points.end(),
                 Search_traits_3(
                   CGAL::internal::boost_::make_function_property_map<Weighted_point, Ret, Construct_point_3>(
                       geom_traits().construct_point_3_object()), geom_traits()));

    Conflict_tester tester(*pbegin, this);
    Point_hider hider(this);
    Cover_manager cover_manager(*this);
    double_vertices = Tr_Base::insert_in_conflict(pbegin, points.end(), hint, tester, hider, cover_manager);

    if(is_large_point_set)
    {
      for(unsigned int i = 0; i < dummy_points_vhs.size(); ++i)
      {
        bool is_hidden = false;
        for(Cell_iterator iter = this->cells_begin(); iter != this->cells_end(); ++iter)
        {
          typename Cell::Point_iterator it = std::find(iter->hidden_points_begin(), iter->hidden_points_end(), dummy_points[i]);
          if(it != iter->hidden_points_end())
          {
            is_hidden = true;
            iter->unhide_point(it);
          }
        }
        if(!is_hidden)
          if(std::find(double_vertices.begin(), double_vertices.end(), dummy_points_vhs[i]) == double_vertices.end())
            remove(dummy_points_vhs[i]);
      }
    }

    return number_of_vertices() - n;
  }
//@}

  void remove(Vertex_handle v)
  {
    typedef CGAL::Periodic_3_regular_triangulation_remove_traits_3< Gt > P3removeT;
    typedef CGAL::Regular_triangulation_3< P3removeT > Euclidean_triangulation;
    typedef Vertex_remover< Euclidean_triangulation > Remover;
    P3removeT remove_traits(domain());
    Euclidean_triangulation tmp(remove_traits);
    Remover remover(this, tmp);
    Conflict_tester ct(this);
    Cover_manager cover_manager(*this);

    Tr_Base::remove(v, remover, ct, cover_manager);

    // Re-insert the points that v was hiding.
    for(typename Remover::Hidden_points_iterator hi = remover.hidden_points_begin();
        hi != remover.hidden_points_end(); ++hi)
    {
      insert(*hi);
    }

    CGAL_triangulation_expensive_assertion(is_valid());
  }

public:
   /** @name Wrapping the traits */ //@{
  bool less_power_distance (const Bare_point &p, const Weighted_point &q, const Weighted_point &r)  const
  {
    return geom_traits().compare_power_distance_3_object()(p, q, r) == SMALLER;
  }

  bool less_power_distance (const Bare_point &p, const Weighted_point &q, const Weighted_point &r,
                            const Offset &o1, const Offset &o2, const Offset &o3)  const
  {
    return geom_traits().compare_power_distance_3_object()(p, q, r, o1, o2, o3) == SMALLER;
  }

  Bare_point construct_weighted_circumcenter(const Weighted_point &p, const Weighted_point &q,
                                             const Weighted_point &r) const
  {
    return geom_traits().construct_weighted_circumcenter_3_object()(p,q,r);
  }
  Bare_point construct_weighted_circumcenter(const Weighted_point &p, const Weighted_point &q,
                                             const Weighted_point &r,
                                             const Offset& o1, const Offset& o2,
                                             const Offset& o3) const
  {
    return geom_traits().construct_weighted_circumcenter_3_object()(p,q,r,o1,o2,o3);
  }
  Bare_point construct_weighted_circumcenter (const Weighted_point &p, const Weighted_point &q,
                                              const Weighted_point &r, const Weighted_point &s) const
  {
    return geom_traits().construct_weighted_circumcenter_3_object()(p,q,r,s);
  }
  Bare_point construct_weighted_circumcenter (const Weighted_point &p, const Weighted_point &q,
                                              const Weighted_point &r, const Weighted_point &s,
                                              const Offset& o1, const Offset& o2,
                                              const Offset& o3, const Offset& o4) const
  {
    return geom_traits().construct_weighted_circumcenter_3_object()(p,q,r,s, o1,o2,o3,o4);
  }

  Comparison_result compare_distance(const Bare_point &p, const Weighted_point &q,
                                     const Weighted_point &r) const {
    return geom_traits().compare_power_distance_3_object()(p, q, r);
  }
  Comparison_result compare_distance(const Bare_point& p, const Weighted_point& q,
                                     const Weighted_point& r,
                                     const Offset &o_p, const Offset &o_q,
                                     const Offset &o_r) const {
    return geom_traits().compare_power_distance_3_object()(p, q, r, o_p, o_q, o_r);
  }

  Oriented_side power_test(const Weighted_point &p, const Weighted_point &q) const
  {
    CGAL_triangulation_precondition(this->equal(p, q));
    return geom_traits().power_side_of_oriented_power_sphere_3_object()(p, q);
  }
  Oriented_side power_test(const Weighted_point &p, const Weighted_point &q,
                           const Weighted_point &r, const Weighted_point &s,
                           const Weighted_point &t,
                           const Offset &o_p, const Offset &o_q,
                           const Offset &o_r, const Offset &o_s,
                           const Offset &o_t) const
  {
    return geom_traits().power_side_of_oriented_power_sphere_3_object()(
             p, q, r, s, t, o_p, o_q, o_r, o_s, o_t);
  }

  Oriented_side side_of_oriented_power_sphere(const Weighted_point &p, const Weighted_point &q,
                                              const Weighted_point &r, const Weighted_point &s,
                                              const Weighted_point &t,
                                              const Offset &o_p, const Offset &o_q,
                                              const Offset &o_r, const Offset &o_s,
                                              const Offset &o_t) const
  {
    return power_test(p,q,r,s,t,o_p,o_q,o_r,o_s,o_t);
  }

  Bounded_side side_of_power_sphere(const Cell_handle& c, const Weighted_point& p,
                                    const Offset & offset = Offset(),
                                    bool perturb = false) const
  {
    Bounded_side bs = ON_UNBOUNDED_SIDE;
    int i = 0;
    // TODO: optimize which copies to check depending on the offsets in
    // the cell.
    while(bs == ON_UNBOUNDED_SIDE && i<8) {
      bs = _side_of_power_sphere(c,p,combine_offsets(offset,int_to_off(i)),perturb);
      i++;
    }
    return bs;
  }

  Bounded_side _side_of_power_sphere(const Cell_handle& c, const Weighted_point& p,
                                     const Offset & offset = Offset(),
                                     bool perturb = false) const;

  Weighted_point construct_weighted_point(const Weighted_point& p, const Offset &o) const
  {
    return geom_traits().construct_weighted_point_3_object()(p, o);
  }
  Weighted_point construct_weighted_point(const Periodic_weighted_point& pp) const
  {
    return construct_weighted_point(pp.first, pp.second);
  }

public:
  /** @name Geometric access functions */
  /// @{
  Weighted_point point(const Periodic_weighted_point& pp) const
  {
    return point(pp, geom_traits().construct_weighted_point_3_object());
  }

  // The following functions return the "real" position in space (unrestrained
  // to the original periodic domain) of the vertices v and c->vertex(idx),
  // respectively

  Weighted_point point(Vertex_handle v) const
  {
    return point(v, geom_traits().construct_weighted_point_3_object());
  }

  Weighted_point point(Cell_handle c, int idx) const
  {
    return point(c, idx, geom_traits().construct_weighted_point_3_object());
  }

  // end of geometric functions
  /// @}

#define CGAL_INCLUDE_FROM_PERIODIC_3_REGULAR_TRIANGULATION_3_H
#include <CGAL/internal/Periodic_3_regular_triangulation_dummy_288.h>
#undef CGAL_INCLUDE_FROM_PERIODIC_3_REGULAR_TRIANGULATION_3_H

  Vertex_handle nearest_power_vertex(const Bare_point& p, Cell_handle start) const
  {
    CGAL_triangulation_precondition(p.x() < domain().xmax());
    CGAL_triangulation_precondition(p.y() < domain().ymax());
    CGAL_triangulation_precondition(p.z() < domain().zmax());
    CGAL_triangulation_precondition(p.x() >= domain().xmin());
    CGAL_triangulation_precondition(p.y() >= domain().ymin());
    CGAL_triangulation_precondition(p.z() >= domain().zmin());

    if(number_of_vertices() == 0)
      return Vertex_handle();

    typename Gt::Construct_weighted_point_3 p2wp =
        geom_traits().construct_weighted_point_3_object();

    Locate_type lt;
    int li, lj;
    Offset query_offset;
    Cell_handle c = locate(p2wp(p), query_offset, lt, li, lj, start);

#ifdef CGAL_PERIODIC_DEBUG_NEAREST_POWER_VERTEX
    std::cout << "nearest power vertex: " << p << std::endl;
    std::cout << "vertices: " << number_of_vertices() << std::endl;
    std::cout << "stored vertices: " << this->number_of_stored_vertices() << std::endl;
    std::cout << "Locate: " << p << std::endl;
    std::cout << "Cell: " << &*c << std::endl;
    std::cout << this->point(c, 0) << std::endl;
    std::cout << this->point(c, 1) << std::endl;
    std::cout << this->point(c, 2) << std::endl;
    std::cout << this->point(c, 3) << std::endl;
    std::cout << "offset query: " << query_offset << std::endl;
    std::cout << "bounded side : " << geom_traits().bounded_side_3_object()(
                   Tetrahedron(this->point(c, 0).point(), this->point(c, 1).point(),
                               this->point(c, 2).point(), this->point(c, 3).point()), p) << std::endl;
    std::cout << "power side: " << geom_traits().power_side_of_bounded_power_sphere_3_object()(
                   this->point(c, 0), this->point(c, 1),
                   this->point(c, 2), this->point(c, 3), p2wp(p)) << std::endl;
    std::cout << "power distance: " << geom_traits().compute_power_distance_to_power_sphere_3_object()(
                   this->point(c, 0), this->point(c, 1),
                   this->point(c, 2), this->point(c, 3), p2wp(p)) << std::endl;
#endif

    // - start with the closest vertex from the located cell.
    // - repeatedly take the nearest of its incident vertices if any
    // - if not, we're done.

    // Take the opposite because periodic_locate() returns the offset such that
    // cell + offset contains 'p' but here we need to move 'p'
    query_offset = this->combine_offsets(Offset(), -query_offset);

    Vertex_handle nearest = nearest_vertex_in_cell(c, p, query_offset);
    Offset offset_of_nearest = get_min_dist_offset(p, query_offset, nearest);

#ifdef CGAL_PERIODIC_DEBUG_NEAREST_POWER_VERTEX
    std::cout << "nearest vertex in cell : " << &*nearest << " : " << nearest->point() << std::endl;
    std::cout << "offset_of_nearest: " << offset_of_nearest << std::endl;
#endif

    std::vector<Vertex_handle> vs;
    vs.reserve(32);
    while(true)
    {
      Vertex_handle tmp = nearest;
#ifdef CGAL_PERIODIC_DEBUG_NEAREST_POWER_VERTEX
      std::cout << "tmp set to : " << &*nearest << " : " << nearest->point()
                << " || offset: " << nearest->offset() << std::endl;
#endif

      adjacent_vertices(nearest, std::back_inserter(vs));
      for(typename std::vector<Vertex_handle>::const_iterator vsit = vs.begin();
          vsit != vs.end(); ++vsit)
      {
        // Can happen in 27-sheeted triangulations composed of few points
        if((*vsit)->point() == nearest->point())
          continue;

        const Offset min_dist_offset = get_min_dist_offset(p, query_offset, *vsit);
        if(compare_distance(p, (*vsit)->point(), tmp->point(),
                            query_offset, min_dist_offset, offset_of_nearest) == SMALLER)
        {
          tmp = *vsit;
          offset_of_nearest = min_dist_offset;
#ifdef CGAL_PERIODIC_DEBUG_NEAREST_POWER_VERTEX
          std::cout << " Closer adjacent vertex: " << &*tmp << " : " << tmp->point()
                    << " || offset " << offset_of_nearest << std::endl;
#endif
        }
      }

      if(tmp == nearest)
        break;

      vs.clear();
      nearest = tmp;
    }

    return get_original_vertex(nearest);
  }

  bool is_Gabriel(Cell_handle c, int i) const
  {
    CGAL_triangulation_precondition(number_of_vertices() != 0);
    typename Geom_traits::Power_side_of_bounded_power_sphere_3
      side_of_bounded_power_sphere =
      geom_traits().power_side_of_bounded_power_sphere_3_object();

    if(side_of_bounded_power_sphere (
         c->vertex(vertex_triple_index(i,0))->point(),
         c->vertex(vertex_triple_index(i,1))->point(),
         c->vertex(vertex_triple_index(i,2))->point(),
         c->vertex(i)->point(),
         get_offset(c,vertex_triple_index(i,0)),
         get_offset(c,vertex_triple_index(i,1)),
         get_offset(c,vertex_triple_index(i,2)),
         get_offset(c,i) ) == ON_BOUNDED_SIDE )
      return false;

    Cell_handle neighbor = c->neighbor(i);
    int in = neighbor->index(c);

    if(side_of_bounded_power_sphere(
         neighbor->vertex(vertex_triple_index(in,0))->point(),
         neighbor->vertex(vertex_triple_index(in,1))->point(),
         neighbor->vertex(vertex_triple_index(in,2))->point(),
         neighbor->vertex(in)->point(),
         get_offset(neighbor,vertex_triple_index(in,0)),
         get_offset(neighbor,vertex_triple_index(in,1)),
         get_offset(neighbor,vertex_triple_index(in,2)),
         get_offset(neighbor, in) ) == ON_BOUNDED_SIDE )
      return false;

    return true;
  }

  bool is_Gabriel(Cell_handle c, int i, int j) const
  {
    CGAL_triangulation_precondition(number_of_vertices() != 0);
    typename Geom_traits::Power_side_of_bounded_power_sphere_3
      side_of_bounded_power_sphere =
      geom_traits().power_side_of_bounded_power_sphere_3_object();

    Facet_circulator fcirc = incident_facets(c,i,j), fdone(fcirc);
    Vertex_handle v1 = c->vertex(i);
    Vertex_handle v2 = c->vertex(j);
    do {
      // test whether the vertex of cc opposite to *fcirc
      // is inside the sphere defined by the edge e = (s, i,j)
      // It is necessary to fetch the offsets from the current cell.
      Cell_handle cc = fcirc->first;
      int i1 = cc->index(v1);
      int i2 = cc->index(v2);
      int i3 = fcirc->second;
      Offset off1 = get_offset(cc, i1);
      Offset off2 = get_offset(cc, i2);
      Offset off3 = get_offset(cc, i3);
      if(side_of_bounded_power_sphere(v1->point(), v2->point(), cc->vertex(i3)->point(),
                                      off1, off2, off3) == ON_BOUNDED_SIDE)
        return false;
    } while(++fcirc != fdone);
    return true;
  }

  bool is_Gabriel(const Facet& f) const
  {
    return is_Gabriel(f.first, f.second);
  }

  bool is_Gabriel(const Edge& e) const
  {
    return is_Gabriel(e.first, e.second, e.third);
  }

  bool is_Gabriel(Vertex_handle v) const
  {
    typename Geom_traits::Power_side_of_bounded_power_sphere_3
      side_of_bounded_orthogonal_sphere =
      geom_traits().power_side_of_bounded_power_sphere_3_object();

    const Bare_point& bp = geom_traits().construct_point_3_object()(v->point());

    Vertex_handle nearest_v = nearest_power_vertex(bp, v->cell());

    // Need to find the offset such that power distance v->point()
    // to nearest_v->point() is minimum
    Offset nearest_v_off = get_min_dist_offset_general(bp, nearest_v);

    return (side_of_bounded_orthogonal_sphere(
              v->point(), nearest_v->point(), Offset(), nearest_v_off) != CGAL::ON_BOUNDED_SIDE);
  }

  Offset get_min_dist_offset(const Bare_point& p, const Offset & o,
                             const Vertex_handle vh) const {
    Offset mdo = get_offset(vh);
    Offset min_off = Offset(0,0,0);
    min_off = (compare_distance(p,vh->point(),vh->point(),
      o,combine_offsets(mdo,min_off),combine_offsets(mdo,Offset(0,0,1)))
        == SMALLER ? min_off : Offset(0,0,1) );
    min_off = (compare_distance(p,vh->point(),vh->point(),
      o,combine_offsets(mdo,min_off),combine_offsets(mdo,Offset(0,1,0)))
        == SMALLER ? min_off : Offset(0,1,0) );
    min_off = (compare_distance(p,vh->point(),vh->point(),
      o,combine_offsets(mdo,min_off),combine_offsets(mdo,Offset(0,1,1)))
        == SMALLER ? min_off : Offset(0,1,1) );
    min_off = (compare_distance(p,vh->point(),vh->point(),
      o,combine_offsets(mdo,min_off),combine_offsets(mdo,Offset(1,0,0)))
        == SMALLER ? min_off : Offset(1,0,0) );
    min_off = (compare_distance(p,vh->point(),vh->point(),
      o,combine_offsets(mdo,min_off),combine_offsets(mdo,Offset(1,0,1)))
        == SMALLER ? min_off : Offset(1,0,1) );
    min_off = (compare_distance(p,vh->point(),vh->point(),
      o,combine_offsets(mdo,min_off),combine_offsets(mdo,Offset(1,1,0)))
        == SMALLER ? min_off : Offset(1,1,0) );
    min_off = (compare_distance(p,vh->point(),vh->point(),
      o,combine_offsets(mdo,min_off),combine_offsets(mdo,Offset(1,1,1)))
        == SMALLER ? min_off : Offset(1,1,1) );
    return combine_offsets(mdo,min_off);
  }

  // In this version, `p` must be in the domain, and we check all possible
  // offsets to find the minimum
  Offset get_min_dist_offset_general(const Bare_point& p, const Vertex_handle vh) const
  {
    CGAL_triangulation_precondition(p.x() < domain().xmax());
    CGAL_triangulation_precondition(p.y() < domain().ymax());
    CGAL_triangulation_precondition(p.z() < domain().zmax());
    CGAL_triangulation_precondition(p.x() >= domain().xmin());
    CGAL_triangulation_precondition(p.y() >= domain().ymin());
    CGAL_triangulation_precondition(p.z() >= domain().zmin());

    Offset min_off = Offset(0,0,0);

    for(int i=-1; i<=1; ++i) {
      for(int j=-1; j<=1; ++j) {
        for(int k=-1; k<=1; ++k)
        {
          if(i==0 && j==0 && k==0)
            continue;

          Offset loc_off(i, j, k);
          if(compare_distance(p, vh->point(), vh->point(), Offset(), min_off, loc_off) == LARGER)
            min_off = loc_off;
        }
      }
    }

    return min_off;
  }

  Vertex_handle nearest_vertex_in_cell(const Cell_handle& c, const Bare_point& p,
                                       const Offset & o) const
  {
    CGAL_triangulation_precondition(number_of_vertices() != 0);
    Vertex_handle nearest = c->vertex(0);
    for(int i=1; i<4; i++) {
      nearest = (compare_distance(p,nearest->point(),c->vertex(i)->point(),
        o,get_offset(c,c->index(nearest)),get_offset(c,i)) == SMALLER) ?
        nearest : c->vertex(i);
    }
    return nearest;
  }

  size_type number_of_hidden_points () const
  {
    size_type count = 0;
    for(Cell_iterator iter = cells_begin(), end_iter = cells_end(); iter != end_iter; ++iter)
      count += std::distance(iter->hidden_points_begin(), iter->hidden_points_end());
    return count;
  }

  bool is_valid(bool verbose = false, int level = 0) const;
  bool is_valid(Cell_handle c, bool verbose = false, int level = 0) const;

protected:
  // Protected, because inheritors(e.g. periodic triangulation for meshing)
  // of the class Periodic_3_regular_triangulation_3 use this class
  class Conflict_tester;

private:
  class Point_hider;

#ifndef CGAL_CFG_OUTOFLINE_TEMPLATE_MEMBER_DEFINITION_BUG
  template <class TriangulationR3> struct Vertex_remover;
#else
  template <class TriangulationR3>
  struct Vertex_remover
  {
    typedef TriangulationR3      Triangulation_R3;

    typedef typename std::vector<Weighted_point>::iterator Hidden_points_iterator;

    typedef Triple < Vertex_handle, Vertex_handle, Vertex_handle > Vertex_triple;

    typedef typename Triangulation_R3::Triangulation_data_structure TDSE;
    typedef typename Triangulation_R3::Cell_handle        CellE_handle;
    typedef typename Triangulation_R3::Vertex_handle      VertexE_handle;
    typedef typename Triangulation_R3::Facet              FacetE;
    typedef typename Triangulation_R3::Finite_cells_iterator Finite_cellsE_iterator;

    typedef Triple< VertexE_handle, VertexE_handle, VertexE_handle >
    VertexE_triple;

    typedef std::map<Vertex_triple,Facet> Vertex_triple_Facet_map;
    typedef std::map<Vertex_triple, FacetE> Vertex_triple_FacetE_map;
    typedef typename Vertex_triple_FacetE_map::iterator
    Vertex_triple_FacetE_map_it;

    Vertex_remover(const Self *t, Triangulation_R3 &tmp_) : _t(t),tmp(tmp_) {}

    const Self *_t;
    Triangulation_R3 &tmp;

    void add_hidden_points(Cell_handle ch) {
      std::copy(ch->hidden_points_begin(), ch->hidden_points_end(),
    std::back_inserter(hidden));
    }

    Hidden_points_iterator hidden_points_begin() {
      return hidden.begin();
    }
    Hidden_points_iterator hidden_points_end() {
      return hidden.end();
    }
    //private:
    // The removal of v may un-hide some points,
    // Space functions output them.
    std::vector<Weighted_point> hidden;
  };
#endif //CGAL_CFG_OUTOFLINE_TEMPLATE_MEMBER_DEFINITION_BUG

public:
  Periodic_bare_point periodic_weighted_circumcenter(Cell_handle c) const {
    return Tr_Base::periodic_circumcenter(c,
      geom_traits().construct_weighted_circumcenter_3_object());
  }

  /** @name Voronoi diagram */ //@{
  Bare_point dual(Cell_handle c) const {
    return Tr_Base::construct_point(periodic_weighted_circumcenter(c));
  }

  bool canonical_dual_segment(Cell_handle c, int i, Periodic_segment_3& ps) const {
    return Tr_Base::canonical_dual_segment(c, i, ps,
      geom_traits().construct_weighted_circumcenter_3_object());
  }

  Periodic_segment_3 dual(const Facet & f) const {
    return dual( f.first, f.second );
  }

  Periodic_segment_3 dual(Cell_handle c, int i) const{
    Periodic_segment_3 ps;
    canonical_dual_segment(c,i,ps);
    return ps;
  }

  template <class OutputIterator>
  OutputIterator dual(const Edge & e, OutputIterator points) const {
    return Tr_Base::dual(e.first, e.second, e.third, points);
  }

  template <class OutputIterator>
  OutputIterator dual(Cell_handle c, int i, int j,
      OutputIterator points) const {
    Tr_Base::dual(c, i, j, points, geom_traits().construct_weighted_circumcenter_3_object());
    return points;
  }

  template <class OutputIterator>
  OutputIterator dual(Vertex_handle v, OutputIterator points) const {
    Tr_Base::dual(v, points, geom_traits().construct_weighted_circumcenter_3_object());
    return points;
  }

  template <class Stream>
  Stream& draw_dual(Stream& os) const {
    return Tr_Base::draw_dual(os, geom_traits().construct_weighted_circumcenter_3_object());
  }

  /// Volume computations
  FT dual_volume(Vertex_handle v) const {
    return Tr_Base::dual_volume(v, geom_traits().construct_weighted_circumcenter_3_object());
  }

  /// Centroid computations
  Bare_point dual_centroid(Vertex_handle v) const {
    return Tr_Base::dual_centroid(
             v, geom_traits().construct_weighted_circumcenter_3_object());
  }
//@}

  template <class OutputIteratorBoundaryFacets, class OutputIteratorCells>
  std::pair<OutputIteratorBoundaryFacets, OutputIteratorCells>
  find_conflicts(const Weighted_point &p, Cell_handle c,
                 OutputIteratorBoundaryFacets bfit, OutputIteratorCells cit) const
  {
    Triple<OutputIteratorBoundaryFacets,OutputIteratorCells,Emptyset_iterator>
      t = find_conflicts(p, c, bfit, cit, Emptyset_iterator());
    return std::make_pair(t.first, t.second);
  }

  template <class OutputIteratorBoundaryFacets, class OutputIteratorCells,
            class OutputIteratorInternalFacets>
  Triple<OutputIteratorBoundaryFacets, OutputIteratorCells,
         OutputIteratorInternalFacets>
  find_conflicts(const Weighted_point &p, Cell_handle c,
                 OutputIteratorBoundaryFacets bfit, OutputIteratorCells cit,
                 OutputIteratorInternalFacets ifit) const;

  /// Returns the vertices on the boundary of the conflict hole.
  template <class OutputIterator>
  OutputIterator vertices_in_conflict(const Weighted_point&p, Cell_handle c,
                                      OutputIterator res) const;

  inline bool
  is_extensible_triangulation_in_1_sheet_h1() const
  {
    if(!is_1_cover())
      return can_be_converted_to_1_sheet();
    return is_extensible_triangulation_in_1_sheet_h2();
  }

  inline bool
  is_extensible_triangulation_in_1_sheet_h2() const
  {
    FT threshold = FT(0.015625) * (domain().xmax()-domain().xmin()) * (domain().xmax()-domain().xmin());

    for(Periodic_tetrahedron_iterator tit = this->periodic_tetrahedra_begin(Tr_Base::UNIQUE);
        tit != this->periodic_tetrahedra_end(Tr_Base::UNIQUE); ++tit)
    {
      if(compare_orthsphere_radius_to_threshold(tit->at(0), tit->at(1),
                                                 tit->at(2), tit->at(3),
                                                 threshold) != CGAL::SMALLER)
        return false;
    }
    return true;
  }
};

template < class Gt, class Tds >
template <class OutputIterator>
OutputIterator
Periodic_3_regular_triangulation_3<Gt,Tds>::vertices_in_conflict(
    const Weighted_point&p, Cell_handle c, OutputIterator res) const
{
  if(number_of_vertices() == 0)
    return res;

  // Get the facets on the boundary of the hole.
  std::vector<Facet> facets;
  find_conflicts(p, c, std::back_inserter(facets), Emptyset_iterator());

  // Then extract uniquely the vertices.
  std::set<Vertex_handle> vertices;
  for(typename std::vector<Facet>::const_iterator i = facets.begin();
       i != facets.end(); ++i) {
    vertices.insert(i->first->vertex((i->second+1)&3));
    vertices.insert(i->first->vertex((i->second+2)&3));
    vertices.insert(i->first->vertex((i->second+3)&3));
  }

  return std::copy(vertices.begin(), vertices.end(), res);
}

template < class Gt, class Tds >
template <class OutputIteratorBoundaryFacets, class OutputIteratorCells,
          class OutputIteratorInternalFacets>
Triple<OutputIteratorBoundaryFacets, OutputIteratorCells,
       OutputIteratorInternalFacets>
Periodic_3_regular_triangulation_3<Gt,Tds>::find_conflicts(
    const Weighted_point& p, Cell_handle c,
    OutputIteratorBoundaryFacets bfit,
    OutputIteratorCells cit,
    OutputIteratorInternalFacets ifit) const
{
  CGAL_triangulation_precondition(number_of_vertices() != 0);

  std::vector<Facet> facets;
  facets.reserve(64);
  std::vector<Cell_handle> cells;
  cells.reserve(32);

  Conflict_tester tester(p, this);
  Triple<typename std::back_insert_iterator<std::vector<Facet> >,
         typename std::back_insert_iterator<std::vector<Cell_handle> >,
         OutputIteratorInternalFacets> tit =
           Tr_Base::find_conflicts(c, tester,
                                   make_triple(std::back_inserter(facets),
                                               std::back_inserter(cells),
                                               ifit));
  ifit = tit.third;

  // Reset the conflict flag on the boundary.
  for(typename std::vector<Facet>::iterator fit=facets.begin();
  fit != facets.end(); ++fit) {
    fit->first->neighbor(fit->second)->tds_data().clear();
    *bfit++ = *fit;
  }

  // Reset the conflict flag in the conflict cells.
  for(typename std::vector<Cell_handle>::iterator ccit=cells.begin();
      ccit != cells.end(); ++ccit) {
    (*ccit)->tds_data().clear();
    *cit++ = *ccit;
  }

  for(typename std::vector<Vertex_handle>::iterator
      voit = this->v_offsets.begin(); voit != this->v_offsets.end(); ++voit) {
    (*voit)->clear_offset();
  }

  this->v_offsets.clear();

  return make_triple(bfit, cit, ifit);
}

template < class Gt, class Tds >
Bounded_side Periodic_3_regular_triangulation_3<Gt,Tds>::
_side_of_power_sphere(const Cell_handle &c, const Weighted_point &q,
                      const Offset &offset, bool perturb ) const
{
  Weighted_point p0 = c->vertex(0)->point(),
                 p1 = c->vertex(1)->point(),
                 p2 = c->vertex(2)->point(),
                 p3 = c->vertex(3)->point();
  Offset o0 = this->get_offset(c,0),
         o1 = this->get_offset(c,1),
         o2 = this->get_offset(c,2),
         o3 = this->get_offset(c,3),
         oq = offset;

  CGAL_triangulation_precondition( orientation(p0, p1, p2, p3, o0, o1, o2, o3) == POSITIVE );

  Oriented_side os = ON_NEGATIVE_SIDE;
  os = side_of_oriented_power_sphere(p0, p1, p2, p3, q, o0, o1, o2, o3, oq);

  if(os != ON_ORIENTED_BOUNDARY || !perturb)
    return (Bounded_side) os;

  // We are now in a degenerate case => we do a symbolic perturbation.
  // We sort the points lexicographically.
  Periodic_weighted_point pts[5] = {std::make_pair(p0,o0), std::make_pair(p1,o1),
                                    std::make_pair(p2,o2), std::make_pair(p3,o3),
                                    std::make_pair(q,oq)};
  const Periodic_weighted_point *points[5] ={&pts[0],&pts[1],&pts[2],&pts[3],&pts[4]};

  std::sort(points, points+5, typename Tr_Base::Perturbation_order(this));

  // We successively look whether the leading monomial, then 2nd monomial
  // of the determinant has non null coefficient.
  for(int i=4; i>1; --i) {
    if(points[i] == &pts[4]) {
      CGAL_triangulation_assertion(orientation(p0, p1, p2, p3, o0, o1, o2, o3)
          == POSITIVE);
      // since p0 p1 p2 p3 are non coplanar and positively oriented
      return ON_UNBOUNDED_SIDE;
    }
    Orientation o;
    if(points[i] == &pts[3] &&
        (o = orientation(p0, p1, p2, q, o0, o1, o2, oq)) != COPLANAR ) {
      return (Bounded_side) o;
    }
    if(points[i] == &pts[2] &&
        (o = orientation(p0, p1, q, p3, o0, o1, oq, o3)) != COPLANAR ) {
      return (Bounded_side) o;
    }
    if(points[i] == &pts[1] &&
        (o = orientation(p0, q, p2, p3, o0, oq, o2, o3)) != COPLANAR ) {
      return (Bounded_side) o;
    }
    if(points[i] == &pts[0] &&
        (o = orientation(q, p1, p2 ,p3, oq, o1, o2, o3)) != COPLANAR ) {
      return (Bounded_side) o;
    }
  }

  CGAL_triangulation_assertion(false);
  return ON_UNBOUNDED_SIDE;
}

template < class Gt, class Tds >
bool
Periodic_3_regular_triangulation_3<Gt,Tds>::
is_valid(bool verbose, int level) const
{
  if(!Tr_Base::is_valid(verbose, level)) {
    if(verbose)
      std::cerr << "Regular: invalid base" << std::endl;
    return false;
  }

  Conflict_tester tester(this);
  if(!is_valid_conflict(tester, verbose, level)) {
    if(verbose)
      std::cerr << "Regular: conflict problems" << std::endl;
    return false;
  }

  if(verbose)
    std::cerr << "Regular valid triangulation" << std::endl;
  return true;
}

template < class GT, class TDS >
bool
Periodic_3_regular_triangulation_3<GT,TDS>::
is_valid(Cell_handle ch, bool verbose, int level) const
{
  bool error = false;
  if(!Tr_Base::is_valid(ch, verbose, level)) {
    error = true;
    if(verbose) {
      std::cerr << "geometrically invalid cell" << std::endl;
      for(int i=0; i<4; i++ )
        std::cerr << ch->vertex(i)->point() << ", ";
      std::cerr << std::endl;
    }
  }

  for(Vertex_iterator vit = vertices_begin(); vit != vertices_end(); ++ vit) {
    const Periodic_weighted_point& pwp = periodic_point(vit);
    for(int i=-1; i<=1; i++) {
      for(int j=-1; j<=1; j++) {
        for(int k=-1; k<=1; k++) {
          const Periodic_weighted_point& ofpwp = std::make_pair(pwp.first,
                                                                pwp.second + Offset(i,j,k));
          if(periodic_point(ch,0) == ofpwp
              || periodic_point(ch,1) == ofpwp
              || periodic_point(ch,2) == ofpwp
              || periodic_point(ch,3) == ofpwp)
            continue;
          if(_side_of_power_sphere(ch, ofpwp.first, ofpwp.second, true)
              != ON_UNBOUNDED_SIDE) {
            error = true;
            if(verbose) {
              std::cerr << "Regular invalid cell" << std::endl;
              for(int i=0; i<4; i++) {
                Periodic_weighted_point pp = periodic_point(ch, i);
                std::cerr <<"("<<pp.first <<","<<pp.second<< "), ";
              }
              std::cerr << std::endl;
            }
          }
        }
      }
    }
  }
  return !error;
}

template < class GT, class Tds >
class Periodic_3_regular_triangulation_3<GT,Tds>::Conflict_tester
{
  // stores a pointer to the triangulation,
  // a point, and an offset
  const Self* t;
  Weighted_point p;
  // stores the offset of a point in 27-cover
  mutable Offset o;

public:
  /// Constructor
  Conflict_tester(const Self* _t) : t(_t), p(Weighted_point()) { }
  Conflict_tester(const Weighted_point& pt, const Self *_t) : t(_t), p(pt) { }

  /** The functor
    *
    * gives true if the circumcircle of c contains p
    */
  bool operator()(const Cell_handle c, const Offset& off) const {
    return (t->_side_of_power_sphere(c, p, t->combine_offsets(o, off), true)
             == ON_BOUNDED_SIDE);
  }

  bool operator()(const Cell_handle c, const Weighted_point& pt,
                  const Offset& off) const {
    return (t->_side_of_power_sphere(c, pt, o + off, true) == ON_BOUNDED_SIDE);
  }

  int compare_weight(const Weighted_point& p, const Weighted_point& q) const {
    return t->power_test(p, q);
  }

  bool test_initial_cell(Cell_handle c, const Offset& off) const {
    return (operator()(c, off));
  }

  void set_point(const Weighted_point& _p) {
    p = _p;
  }

  void set_offset(const Offset& off) const {
    o = off;
  }

  const Offset& get_offset() const {
    return o;
  }

  const Weighted_point& point() const {
    return p;
  }

};

template < class GT, class Tds>
class Periodic_3_regular_triangulation_3<GT,Tds>::Point_hider
{
  Self* t;
  mutable std::vector<Vertex_handle> vertices;
  mutable std::vector<Weighted_point> hidden_points;
  mutable bool is_original_cube;

public:
  Point_hider(Self* tr) : t(tr), is_original_cube(false) { }

  void set_original_cube (bool b) const {
    is_original_cube = b;
  }

  template <class InputIterator>
  inline void set_vertices(InputIterator start, InputIterator end) const
  {
    while(start != end) {
      std::copy((*start)->hidden_points_begin(),
                (*start)->hidden_points_end(),
                std::back_inserter(hidden_points));

      for(int i=0; i<=3; i++) {
        Vertex_handle v = (*start)->vertex(i);
        if(v->cell() != Cell_handle()) {
          vertices.push_back(v);
          v->set_cell(Cell_handle());
        }
      }
      start ++;
    }
  }

  inline void reinsert_vertices(Vertex_handle v)
  {
    Locate_type lt = Locate_type();
    int li=0, lj=0;

    Cell_handle hc = v->cell();
    for(typename std::vector<Vertex_handle>::iterator
        vi = vertices.begin(); vi != vertices.end(); ++vi) {
      if((*vi)->cell() != Cell_handle())
        continue;
      if(is_original_cube) {
        hc = t->locate((*vi)->point(), lt, li, lj, hc);
        hc->hide_point((*vi)->point());
      }
      t->delete_vertex(*vi);
    }
    vertices.clear();
      for(typename std::vector<Weighted_point>::iterator
          hp = hidden_points.begin(); hp != hidden_points.end(); ++hp) {
        hc = t->locate(*hp, lt, li, lj, hc);
        hc->hide_point(*hp);
      }
      hidden_points.clear();
  }

  inline Vertex_handle replace_vertex(Cell_handle c, int index, const Weighted_point& p)
  {
    Vertex_handle v = c->vertex(index);
    c->hide_point(v->point());
    v->set_point(p);
    return v;
  }

  inline void hide_point(Cell_handle c, const Weighted_point& p)
  {
    if(is_original_cube)
      c->hide_point(p);
  }

//  inline void hide(Weighted_point&, Cell_handle ) const  // useless?
//  {
//    CGAL_triangulation_assertion(false);
//  }
//
//  inline void do_hide(const Weighted_point&, Cell_handle ) const // useless?
//  {
//    CGAL_triangulation_assertion(false);
//  }

//  template < class Tester >
//  inline bool replace_vertex(const Weighted_point&, Vertex_handle, const Tester&) const // useless?
//  {
//    return true;
//  }
//
//  template <class Conflict_tester>
//  inline void hide_points(Vertex_handle,
//                          const Conflict_tester &)
//  {
//    // No points to hide in the Delaunay triangulation.
//  }
};

#ifndef CGAL_CFG_OUTOFLINE_TEMPLATE_MEMBER_DEFINITION_BUG
template <class GT, class Tds>
template <class TriangulationR3>
struct Periodic_3_regular_triangulation_3<GT,Tds>::Vertex_remover
{
  typedef TriangulationR3                                Triangulation_R3;

  typedef typename std::vector<Weighted_point>::iterator Hidden_points_iterator;

  typedef Triple < Vertex_handle, Vertex_handle, Vertex_handle > Vertex_triple;

  typedef typename Triangulation_R3::Triangulation_data_structure TDSE;
  typedef typename Triangulation_R3::Cell_handle                  CellE_handle;
  typedef typename Triangulation_R3::Vertex_handle                VertexE_handle;
  typedef typename Triangulation_R3::Facet                        FacetE;
  typedef typename Triangulation_R3::Finite_cells_iterator        Finite_cellsE_iterator;

  typedef Triple<VertexE_handle, VertexE_handle, VertexE_handle>  VertexE_triple;

  typedef std::map<Vertex_triple,Facet> Vertex_triple_Facet_map;
  typedef std::map<Vertex_triple, FacetE> Vertex_triple_FacetE_map;
  typedef typename Vertex_triple_FacetE_map::iterator
  Vertex_triple_FacetE_map_it;

  Vertex_remover(const Self* t, Triangulation_R3& tmp_) : _t(t),tmp(tmp_) {}

  const Self* _t;
  Triangulation_R3 &tmp;

  void add_hidden_points(Cell_handle ch) {
    std::copy(ch->hidden_points_begin(), ch->hidden_points_end(),
              std::back_inserter(hidden));
  }

  Hidden_points_iterator hidden_points_begin() {
    return hidden.begin();
  }
  Hidden_points_iterator hidden_points_end() {
    return hidden.end();
  }
  private:
  // The removal of v may un-hide some points,
  // Space functions output them.
  std::vector<Weighted_point> hidden;
};
#endif //CGAL_CFG_OUTOFLINE_TEMPLATE_MEMBER_DEFINITION_BUG

template < class GT, class TDS >
std::istream &
operator>> (std::istream& is, Periodic_3_regular_triangulation_3<GT, TDS>& tr)
{
  typedef Periodic_3_regular_triangulation_3<GT,TDS>   P3RT3;
  typedef typename P3RT3::Tr_Base                      Tr_Base;

  is >> static_cast<Tr_Base&>(tr);

  tr.insert_cells_with_too_big_orthoball(tr.cells_begin(), tr.cells_end());

  CGAL_triangulation_expensive_assertion( tr.is_valid() );
  return is;
}

} // namespace CGAL

#endif
