// Copyright (c) 1999-2004,2006-2009,2013-2015   INRIA Sophia-Antipolis (France).
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
//
// Author(s)     : Monique Teillaud <Monique.Teillaud@inria.fr>
//                 Aymeric Pelle <Aymeric.Pelle@sophia.inria.fr>
#ifndef CGAL_PERIODIC_3_REGULAR_TRIANGULATION_3_H
#define CGAL_PERIODIC_3_REGULAR_TRIANGULATION_3_H

#include <CGAL/Periodic_3_triangulation_3.h>
#include <CGAL/spatial_sort.h>

// Needed by remove to fill the hole.
#include <CGAL/Periodic_3_regular_triangulation_remove_traits_3.h>
#include <CGAL/Regular_triangulation_3.h>

#include <boost/unordered_set.hpp>


namespace CGAL
{
template < class Gt,
		   class Tds = Triangulation_data_structure_3 < Triangulation_vertex_base_3<Gt, Periodic_3_triangulation_ds_vertex_base_3<> >,
		                                                Regular_triangulation_cell_base_3<Gt, Periodic_3_triangulation_ds_cell_base_3<> > >
		 >
class Periodic_3_regular_triangulation_3 : public Periodic_3_triangulation_3<Gt,Tds>
{
	typedef Periodic_3_regular_triangulation_3<Gt,Tds>  Self;

public:
	typedef Periodic_3_triangulation_3<Gt,Tds>  Base;

	typedef Gt                                    Geometric_traits;
	typedef Tds                                   Triangulation_data_structure;

	typedef Geometric_traits                      Geom_traits;
	typedef typename Gt::FT                       FT;

	typedef typename Gt::Weighted_point_3         Weighted_point;
	typedef typename Gt::Bare_point               Bare_point;
	typedef typename Gt::Segment_3                Segment;
	typedef typename Gt::Triangle_3               Triangle;
	typedef typename Gt::Tetrahedron_3            Tetrahedron;

	typedef typename Base::Periodic_point         Periodic_point;
	typedef typename Base::Periodic_segment       Periodic_segment;
	typedef typename Base::Periodic_triangle      Periodic_triangle;
	typedef typename Base::Periodic_tetrahedron   Periodic_tetrahedron;

  typedef typename Base::Periodic_tetrahedron_iterator  Periodic_tetrahedron_iterator;

	typedef typename Base::Cell_handle            Cell_handle;
	typedef typename Base::Vertex_handle          Vertex_handle;

	typedef typename Base::Cell                   Cell;
	typedef typename Base::Vertex                 Vertex;
	typedef typename Base::Facet                  Facet;
	typedef typename Base::Edge                   Edge;

	typedef typename Base::Cell_circulator        Cell_circulator;
	typedef typename Base::Facet_circulator       Facet_circulator;
	typedef typename Base::Cell_iterator          Cell_iterator;
	typedef typename Base::Facet_iterator         Facet_iterator;
	typedef typename Base::Edge_iterator          Edge_iterator;
	typedef typename Base::Vertex_iterator        Vertex_iterator;

	typedef typename Base::All_cells_iterator     All_cells_iterator;
	typedef typename Base::All_facets_iterator    All_facets_iterator;
	typedef typename Base::All_edges_iterator     All_edges_iterator;
	typedef typename Base::All_vertices_iterator  All_vertices_iterator;

  typedef typename Base::size_type              size_type;
	typedef typename Base::difference_type        difference_type;

	typedef typename Base::Locate_type            Locate_type;
	typedef typename Base::Iterator_type          Iterator_type;

	typedef typename Base::Offset                 Offset;
	typedef typename Base::Iso_cuboid             Iso_cuboid;
	typedef typename Base::Covering_sheets        Covering_sheets;

#ifndef CGAL_CFG_USING_BASE_MEMBER_BUG_2
using Base::cw;
using Base::ccw;
using Base::domain;
using Base::geom_traits;
using Base::int_to_off;
using Base::number_of_sheets;
using Base::number_of_vertices;
using Base::number_of_edges;
using Base::number_of_facets;
using Base::number_of_cells;
using Base::cells_begin;
using Base::cells_end;
using Base::vertices_begin;
using Base::vertices_end;
using Base::facets_begin;
using Base::facets_end;
using Base::tds;
using Base::next_around_edge;
using Base::vertex_triple_index;
using Base::mirror_vertex;
using Base::orientation;
using Base::swap;
using Base::is_1_cover;
using Base::is_virtual;
using Base::point;
using Base::set_offsets;
#endif

// For strict-ansi compliance
using Base::adjacent_vertices;
using Base::combine_offsets;
using Base::get_offset;
using Base::get_original_vertex;
using Base::get_location_offset;
using Base::neighbor_offset;
using Base::incident_edges;
using Base::incident_facets;
using Base::incident_cells;
using Base::is_valid_conflict;
using Base::locate;
using Base::periodic_point;
using Base::segment;

private:
  struct Cell_handle_hash : public std::unary_function<Cell_handle, std::size_t>
  {
    std::size_t operator()(const Cell_handle& ch) const { return boost::hash<typename Cell_handle::pointer>()(&*ch); }
  };
  boost::unordered_set<Cell_handle, Cell_handle_hash> cells_with_too_big_orthoball;

  class Cover_manager
  {
    Periodic_3_regular_triangulation_3& tr;

  public:
    Cover_manager (Periodic_3_regular_triangulation_3& tr)
  : tr(tr)
  {}

    void create_initial_triangulation()
    {
      tr.create_initial_triangulation();
    }

    template <class CellIt>
    void delete_too_long_edges(const CellIt begin, const CellIt end)
    {
      tr.delete_too_long_edges(begin, end);
    }

    template <class CellIt>
    void insert_too_long_edges(Vertex_handle v, const CellIt begin, const CellIt end)
    {
      tr.insert_too_long_edges(v, begin, end);
    }

    bool can_be_converted_to_1_sheet () const
    {
      return tr.can_be_converted_to_1_sheet();
    }

    bool update_cover_data_during_management (Cell_handle new_ch, const std::vector<Cell_handle>& new_cells)
    {
      return tr.update_cover_data_during_management(new_ch, new_cells);
    }
  };

public:
	/** @name Creation */ //@{
	Periodic_3_regular_triangulation_3 (const Iso_cuboid& domain = Iso_cuboid(0, 0, 0, 1, 1, 1),
			                            const Geometric_traits& gt = Geometric_traits())
	: Base(domain, gt)
	{
	}

  // copy constructor duplicates vertices and cells
  Periodic_3_regular_triangulation_3 (const Periodic_3_regular_triangulation_3& tr)
  : Base(tr)
  {
    if (is_1_cover()) {
      tds() = tr.tds();
    } else {
      this->copy_multiple_covering(tr);
    }
    CGAL_triangulation_expensive_postcondition(*this == tr);
    CGAL_triangulation_expensive_postcondition( is_valid() );
  }

  void copy_multiple_covering(const Periodic_3_regular_triangulation_3& tr) {
    // Write the respective offsets in the vertices to make them
    // automatically copy with the tds.
    for (Vertex_iterator vit = tr.vertices_begin() ;
         vit != tr.vertices_end() ; ++vit) {
      vit->set_offset(tr.get_offset(vit));
    }
    // copy the tds
    tds() = tr.tds();
    // make a list of all vertices that belong to the original
    // domain and initialize the basic structure of
    // virtual_vertices_reverse
    std::list<Vertex_handle> vlist;
    for (Vertex_iterator vit = vertices_begin() ;
         vit != vertices_end() ; ++vit) {
      if (vit->offset() == Offset()) {
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
    for (Vertex_iterator vit2 = vertices_begin() ;
         vit2 != vertices_end() ; ++vit2) {
      if (vit2->offset() != Offset()) {
        //TODO: use some binding, maybe boost instead of the Finder.
        typename std::list<Vertex_handle>::iterator vlist_it
    = std::find_if(vlist.begin(), vlist.end(),
             typename Base::Finder(this,vit2->point()));
        Offset off = vit2->offset();
        this->virtual_vertices.insert(std::make_pair(vit2,
                 std::make_pair(*vlist_it,off)));
        this->virtual_vertices_reverse.find(*vlist_it)
    ->second[9*off[0]+3*off[1]+off[2]-1]=vit2;
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

    insert_too_long_edges(tr.cells_begin(), tr.cells_end());
  }

  template < typename InputIterator >
  Periodic_3_regular_triangulation_3(InputIterator first, InputIterator last,
      const Iso_cuboid& domain = Iso_cuboid(0,0,0,1,1,1),
      const Geometric_traits& gt = Geometric_traits(),
      bool is_large_point_set = false)
    : Base(domain, gt)
  {
    insert(first, last, is_large_point_set);
  }

  Periodic_3_regular_triangulation_3 operator= (Periodic_3_regular_triangulation_3 tr)
  {
    tr.swap(*this);
    return *this;
  }

  void swap(Periodic_3_regular_triangulation_3&tr) {
    std::swap(cells_with_too_big_orthoball,tr.cells_with_too_big_orthoball);
    Base::swap(tr);
  }

  void create_initial_triangulation()
  {
    CGAL_triangulation_assertion( cells_with_too_big_orthoball.empty() );

    for (Cell_iterator iter = cells_begin(), end_iter = cells_end(); iter != end_iter; ++iter)
      cells_with_too_big_orthoball.insert(iter);
  }

  template <class CellIt>
  void delete_too_long_edges(CellIt begin, const CellIt end)
  {
    for (; begin != end; ++begin)
    {
      typename boost::unordered_set<Cell_handle>::iterator iter = cells_with_too_big_orthoball.find(*begin);
      if (iter != cells_with_too_big_orthoball.end())
      {
        cells_with_too_big_orthoball.erase(iter);
      }
    }
  }

  FT squared_orthoball_radius (const Periodic_point& p0, const Periodic_point& p1, const Periodic_point& p2, const Periodic_point& p3) const
  {
    typename Geometric_traits::Construct_weighted_circumcenter_3 construct_weighted_circumcenter_3
                                          = geom_traits().construct_weighted_circumcenter_3_object();

    Bare_point weighted_circumcenter = construct_weighted_circumcenter_3(
        p0.first,  p1.first,  p2.first,  p3.first,
        p0.second, p1.second, p2.second, p3.second);
    Weighted_point pt = point(p0);
    FT ao_2 = squared_distance(static_cast<const Bare_point&>(pt), weighted_circumcenter);
    FT io_2 = ao_2 - pt.weight();
    return io_2;
  }

  FT squared_orthoball_radius (Cell_handle cell)
  {
    Periodic_point p0 = periodic_point(cell, 0);
    Periodic_point p1 = periodic_point(cell, 1);
    Periodic_point p2 = periodic_point(cell, 2);
    Periodic_point p3 = periodic_point(cell, 3);

    return squared_orthoball_radius(p0, p1, p2, p3);
  }

  template <class CellIt>
  void insert_too_long_edges(Vertex_handle /*v*/, CellIt begin, const CellIt end)
  {
    FT threshold = FT(0.015625) * (domain().xmax()-domain().xmin()) * (domain().xmax()-domain().xmin());
    for (; begin != end; ++begin)
    {
      if (squared_orthoball_radius(*begin) >= threshold)
      {
        cells_with_too_big_orthoball.insert(*begin);
      }
    }
  }

  void insert_too_long_edges(Cell_iterator begin, Cell_iterator end)
  {
    FT threshold = FT(0.015625) * (domain().xmax()-domain().xmin()) * (domain().xmax()-domain().xmin());
    for (; begin != end; ++begin)
    {
      if (squared_orthoball_radius(begin) >= threshold)
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

    if (squared_orthoball_radius(new_ch) >= threshold)
    {
      if (is_1_cover())
      {
        tds().delete_cells(new_cells.begin(), new_cells.end());
        this->convert_to_27_sheeted_covering();
        result = true;
      }
      else
        cells_with_too_big_orthoball.insert(new_ch);
    }

    return result;
  }

  virtual void update_cover_data_after_converting_to_27_sheeted_covering ()
  {
    FT threshold = FT(0.015625) * (domain().xmax()-domain().xmin()) * (domain().xmax() - domain().xmin());
    for (Cell_iterator iter = cells_begin(), end_iter = cells_end(); iter != end_iter; ++iter)
    {
      if (squared_orthoball_radius(iter) >= threshold)
      {
        cells_with_too_big_orthoball.insert(iter);
      }
    }
  }

  virtual void update_cover_data_after_setting_domain () {}

  virtual void reinsert_hidden_points_after_converting_to_1_sheeted (std::vector<Weighted_point>& hidden_points)
  {
    while (hidden_points.size())
    {
      insert(hidden_points.back());
      hidden_points.pop_back();
    }
  }

  /** @name Insertion */ //@{
   Vertex_handle insert(const Weighted_point& point, Cell_handle start = Cell_handle()) {
     Conflict_tester tester(point, this);
     Point_hider hider(this);
     Cover_manager cover_manager(*this);
     CGAL_triangulation_precondition(point.weight() >= 0);
     CGAL_triangulation_precondition_msg(point.weight() < ( FT(0.015625) * (domain().xmax()-domain().xmin()) * (domain().xmax()-domain().xmin()) ),
         "point.weight() < 1/64 * domain_size * domain_size");
     return Base::insert_in_conflict(point, start, tester, hider, cover_manager);
   }

   Vertex_handle insert(const Weighted_point& point, Locate_type lt, Cell_handle c,
        int li, int lj) {
      Conflict_tester tester(point, this);
      Point_hider hider(this);
      Cover_manager cover_manager(*this);
      CGAL_triangulation_precondition(point.weight() >= 0);
      CGAL_triangulation_precondition_msg(point.weight() < ( FT(0.015625) * (domain().xmax()-domain().xmin()) * (domain().xmax()-domain().xmin()) ),
          "point.weight() < 1/64 * domain_size * domain_size");
      return Base::insert_in_conflict(point,lt,c,li,lj, tester,hider,cover_manager);
    }

   template < class InputIterator >
   std::ptrdiff_t insert(InputIterator first, InputIterator last,
       bool is_large_point_set = false)
  {
    if (first == last)
      return 0;

    CGAL_triangulation_precondition_code
    (
        bool precondition_is_satisfied = true;
        FT upper_bound = FT(0.015625) * (domain().xmax()-domain().xmin()) * (domain().xmax()-domain().xmin());
        for (InputIterator pc_first = first, pc_last = last; pc_first != pc_last; ++pc_first)
          if (pc_first->weight() < FT(0) || pc_first->weight() >= upper_bound)
          {
            precondition_is_satisfied = false;
            break;
          }
    )
    CGAL_triangulation_precondition_msg(precondition_is_satisfied,
        "0 <= point.weight() < 1/64 * domain_size * domain_size");

    size_type n = number_of_vertices();
    // The heuristic discards the existing triangulation so it can only be
    // applied to empty triangulations.
    if (n != 0)
      is_large_point_set = false;

    std::vector<Weighted_point> points(first, last);
    std::random_shuffle(points.begin(), points.end());
    Cell_handle hint;
    std::vector<Vertex_handle> dummy_points_vhs, double_vertices;
    std::vector<Weighted_point> dummy_points;
    typename std::vector<Weighted_point>::iterator pbegin = points.begin();
    if (is_large_point_set)
    {
      dummy_points_vhs = insert_dummy_points();
      dummy_points.reserve(dummy_points_vhs.size());
      for (typename std::vector<Vertex_handle>::iterator iter = dummy_points_vhs.begin(), end_iter = dummy_points_vhs.end(); iter != end_iter; ++iter)
        dummy_points.push_back((*iter)->point());
    }
    else
      while (!is_1_cover())
      {
        insert(*pbegin);
        ++pbegin;
        if (pbegin == points.end())
          return number_of_vertices() - n;
      }

    // Use Geom_traits::K for efficiency: spatial_sort creates a lot
    // of copies of the traits but does not need the domain that is
    // stored in it.
    spatial_sort(pbegin, points.end(), typename Geom_traits::K());

    Conflict_tester tester(*pbegin, this);
    Point_hider hider(this);
    Cover_manager cover_manager(*this);
    double_vertices = Base::insert_in_conflict(pbegin, points.end(), hint, tester, hider, cover_manager);

    if (is_large_point_set)
    {
      for (unsigned int i = 0; i < dummy_points_vhs.size(); ++i)
      {
        bool is_hidden = false;
        for (Cell_iterator iter = this->cells_begin(); iter != this->cells_end(); ++iter)
        {
          typename Cell::Point_iterator it = std::find(iter->hidden_points_begin(), iter->hidden_points_end(), dummy_points[i]);
          if (it != iter->hidden_points_end())
          {
            is_hidden = true;
            iter->unhide_point(it);
          }
        }
        if (!is_hidden)
          if (std::find(double_vertices.begin(), double_vertices.end(), dummy_points_vhs[i]) == double_vertices.end())
            remove(dummy_points_vhs[i]);
      }
    }

    return number_of_vertices() - n;
  }
   //@}

   void remove(Vertex_handle v)
   {
     typedef CGAL::Periodic_3_regular_triangulation_remove_traits_3< Gt > P3removeT;
     typedef CGAL::Regular_triangulation_3< P3removeT >
       Euclidean_triangulation;
     typedef Vertex_remover< Euclidean_triangulation > Remover;
     P3removeT remove_traits(domain());
     Euclidean_triangulation tmp(remove_traits);
     Remover remover(this, tmp);
     Conflict_tester ct(this);
     Cover_manager cover_manager(*this);

     Base::remove(v, remover, ct, cover_manager);

     // Re-insert the points that v was hiding.
     for (typename Remover::Hidden_points_iterator hi = remover.hidden_points_begin();
          hi != remover.hidden_points_end();
          ++hi)
     {
         insert(*hi);
     }

     CGAL_triangulation_expensive_assertion(is_valid());
   }

protected:
	bool less_power_distance (const Bare_point &p, const Weighted_point &q, const Weighted_point &r)  const
	{
		return geom_traits().compare_power_distance_3_object()(p, q, r) == SMALLER;
	}

	bool less_power_distance (const Bare_point &p, const Weighted_point &q, const Weighted_point &r,
			                  const Offset &o1, const Offset &o2, const Offset &o3)  const
	{
		return geom_traits().compare_power_distance_3_object()(p, q, r, o1, o2, o3) == SMALLER;
	}

	Bare_point construct_weighted_circumcenter (const Weighted_point &p, const Weighted_point &q, const Weighted_point &r, const Weighted_point &s) const
	{
		return geom_traits().construct_weighted_circumcenter_3_object()(p,q,r,s);
	}

	Bare_point construct_weighted_circumcenter (const Weighted_point &p, const Weighted_point &q, const Weighted_point &r, const Weighted_point &s,
			                                    const Offset& o1, const Offset& o2, const Offset& o3, const Offset& o4) const
	{
		return geom_traits().construct_weighted_circumcenter_3_object()(p,q,r,s, o1,o2,o3,o4);
	}

	Bare_point construct_weighted_circumcenter(const Weighted_point &p, const Weighted_point &q, const Weighted_point &r) const
	{
		return geom_traits().construct_weighted_circumcenter_3_object()(p,q,r);
	}

	Bare_point construct_weighted_circumcenter(const Weighted_point &p, const Weighted_point &q, const Weighted_point &r,
			                                   const Offset& o1, const Offset& o2, const Offset& o3) const
	{
		return geom_traits().construct_weighted_circumcenter_3_object()(p,q,r);
	}

//protected:
public:
	std::vector<Vertex_handle> insert_dummy_points ()
	{
	  this->clear();

	  std::vector<Vertex_handle> vertices;
	  vertices.reserve(6*6*8);
	  Cell_handle cells[1728];

    FT domain_x(domain().xmax() - domain().xmin());
    FT domain_y(domain().ymax() - domain().ymin());
    FT domain_z(domain().zmax() - domain().zmin());

	  for (unsigned i = 0; i < 6; ++i)
	    for (unsigned j = 0; j < 6; ++j)
	      for (unsigned k = 0; k < 8; ++k)
	      {
	        Vertex_handle vh = tds().create_vertex();
	        vertices.push_back(vh);

	        FT x = (FT(i) * domain_x / FT(6)) + domain().xmin();
	        if (k % 2)
	          x += FT(1) * domain_x / FT(12);
	        FT y = (FT(j) * domain_y / FT(6)) + domain().ymin();
	        if (k % 2)
	          y += FT(1) * domain_y / FT(12);
	        FT z = (FT(k) * domain_z / FT(8)) + domain().zmin();

	        vh->set_point(Weighted_point(Bare_point(x, y, z), 0));
	      }

	  for (Cell_handle* iter = cells, * end_iter = cells + 1728; iter != end_iter; ++iter)
	    *iter = tds().create_cell();

	  cells[0]->set_vertices(vertices[31],vertices[279],vertices[32],vertices[39]);
	  cells[0]->set_neighbors(cells[1220],cells[495],cells[1015],cells[1489]);
	  set_offsets(cells[0],4,0,5,4);
	  cells[1]->set_vertices(vertices[48],vertices[1],vertices[0],vertices[41]);
	  cells[1]->set_neighbors(cells[1402],cells[567],cells[123],cells[295]);
	  set_offsets(cells[1],2,2,2,0);
	  cells[2]->set_vertices(vertices[23],vertices[24],vertices[72],vertices[31]);
	  cells[2]->set_neighbors(cells[499],cells[428],cells[1261],cells[195]);
	  set_offsets(cells[2],0,1,1,0);
	  cells[3]->set_vertices(vertices[59],vertices[11],vertices[51],vertices[60]);
	  cells[3]->set_neighbors(cells[221],cells[7],cells[298],cells[424]);
	  set_offsets(cells[3],0,0,0,0);
	  cells[4]->set_vertices(vertices[226],vertices[169],vertices[218],vertices[178]);
	  cells[4]->set_neighbors(cells[959],cells[1337],cells[47],cells[873]);
	  set_offsets(cells[4],0,0,0,0);
	  cells[5]->set_vertices(vertices[62],vertices[22],vertices[15],vertices[70]);
	  cells[5]->set_neighbors(cells[447],cells[341],cells[288],cells[18]);
	  set_offsets(cells[5],0,0,0,0);
	  cells[6]->set_vertices(vertices[68],vertices[116],vertices[108],vertices[59]);
	  cells[6]->set_neighbors(cells[478],cells[498],cells[549],cells[193]);
	  set_offsets(cells[6],0,0,0,0);
	  cells[7]->set_vertices(vertices[60],vertices[59],vertices[108],vertices[51]);
	  cells[7]->set_neighbors(cells[210],cells[658],cells[3],cells[498]);
	  set_offsets(cells[7],0,0,0,0);
	  cells[8]->set_vertices(vertices[54],vertices[6],vertices[7],vertices[14]);
	  cells[8]->set_neighbors(cells[266],cells[171],cells[62],cells[103]);
	  set_offsets(cells[8],0,0,0,0);
	  cells[9]->set_vertices(vertices[3],vertices[251],vertices[11],vertices[10]);
	  cells[9]->set_neighbors(cells[912],cells[185],cells[945],cells[803]);
	  set_offsets(cells[9],4,0,4,4);
	  cells[10]->set_vertices(vertices[7],vertices[8],vertices[0],vertices[48]);
	  cells[10]->set_neighbors(cells[295],cells[27],cells[188],cells[963]);
	  set_offsets(cells[10],0,1,1,1);
	  cells[11]->set_vertices(vertices[68],vertices[20],vertices[11],vertices[60]);
	  cells[11]->set_neighbors(cells[179],cells[298],cells[239],cells[144]);
	  set_offsets(cells[11],0,0,0,0);
	  cells[12]->set_vertices(vertices[5],vertices[54],vertices[6],vertices[45]);
	  cells[12]->set_neighbors(cells[175],cells[1425],cells[574],cells[62]);
	  set_offsets(cells[12],2,2,2,0);
	  cells[13]->set_vertices(vertices[55],vertices[7],vertices[56],vertices[15]);
	  cells[13]->set_neighbors(cells[112],cells[314],cells[66],cells[230]);
	  set_offsets(cells[13],0,0,1,0);
	  cells[14]->set_vertices(vertices[91],vertices[138],vertices[90],vertices[50]);
	  cells[14]->set_neighbors(cells[821],cells[156],cells[696],cells[820]);
	  set_offsets(cells[14],0,0,0,2);
	  cells[15]->set_vertices(vertices[3],vertices[11],vertices[51],vertices[58]);
	  cells[15]->set_neighbors(cells[424],cells[249],cells[185],cells[221]);
	  set_offsets(cells[15],0,0,0,0);
	  cells[16]->set_vertices(vertices[176],vertices[168],vertices[216],vertices[169]);
	  cells[16]->set_neighbors(cells[994],cells[254],cells[571],cells[1314]);
	  set_offsets(cells[16],0,0,0,0);
	  cells[17]->set_vertices(vertices[57],vertices[58],vertices[106],vertices[66]);
	  cells[17]->set_neighbors(cells[96],cells[679],cells[263],cells[431]);
	  set_offsets(cells[17],0,0,0,0);
	  cells[18]->set_vertices(vertices[62],vertices[14],vertices[15],vertices[22]);
	  cells[18]->set_neighbors(cells[1494],cells[5],cells[79],cells[189]);
	  set_offsets(cells[18],0,0,0,0);
	  cells[19]->set_vertices(vertices[73],vertices[65],vertices[120],vertices[113]);
	  cells[19]->set_neighbors(cells[272],cells[766],cells[618],cells[738]);
	  set_offsets(cells[19],0,0,0,0);
	  cells[20]->set_vertices(vertices[122],vertices[114],vertices[115],vertices[67]);
	  cells[20]->set_neighbors(cells[742],cells[379],cells[397],cells[517]);
	  set_offsets(cells[20],0,0,0,0);
	  cells[21]->set_vertices(vertices[62],vertices[63],vertices[15],vertices[55]);
	  cells[21]->set_neighbors(cells[314],cells[66],cells[487],cells[341]);
	  set_offsets(cells[21],0,0,0,0);
	  cells[22]->set_vertices(vertices[68],vertices[76],vertices[21],vertices[69]);
	  cells[22]->set_neighbors(cells[300],cells[488],cells[258],cells[481]);
	  set_offsets(cells[22],0,0,0,0);
	  cells[23]->set_vertices(vertices[102],vertices[110],vertices[62],vertices[55]);
	  cells[23]->set_neighbors(cells[487],cells[371],cells[657],cells[627]);
	  set_offsets(cells[23],0,0,0,0);
	  cells[24]->set_vertices(vertices[174],vertices[166],vertices[167],vertices[119]);
	  cells[24]->set_neighbors(cells[1037],cells[869],cells[785],cells[1075]);
	  set_offsets(cells[24],0,0,0,0);
	  cells[25]->set_vertices(vertices[161],vertices[153],vertices[208],vertices[201]);
	  cells[25]->set_neighbors(cells[967],cells[138],cells[1027],cells[1217]);
	  set_offsets(cells[25],0,0,0,0);
	  cells[26]->set_vertices(vertices[29],vertices[30],vertices[78],vertices[38]);
	  cells[26]->set_neighbors(cells[464],cells[359],cells[1583],cells[404]);
	  set_offsets(cells[26],0,0,0,0);
	  cells[27]->set_vertices(vertices[47],vertices[0],vertices[48],vertices[7]);
	  cells[27]->set_neighbors(cells[10],cells[609],cells[1353],cells[363]);
	  set_offsets(cells[27],0,3,3,2);
	  cells[28]->set_vertices(vertices[18],vertices[10],vertices[9],vertices[58]);
	  cells[28]->set_neighbors(cells[39],cells[75],cells[106],cells[1006]);
	  set_offsets(cells[28],0,0,0,0);
	  cells[29]->set_vertices(vertices[166],vertices[109],vertices[158],vertices[118]);
	  cells[29]->set_neighbors(cells[739],cells[690],cells[224],cells[182]);
	  set_offsets(cells[29],0,0,0,0);
	  cells[30]->set_vertices(vertices[34],vertices[74],vertices[26],vertices[25]);
	  cells[30]->set_neighbors(cells[242],cells[1550],cells[520],cells[186]);
	  set_offsets(cells[30],0,0,0,0);
	  cells[31]->set_vertices(vertices[3],vertices[4],vertices[52],vertices[12]);
	  cells[31]->set_neighbors(cells[337],cells[120],cells[1457],cells[331]);
	  set_offsets(cells[31],0,0,0,0);
	  cells[32]->set_vertices(vertices[66],vertices[19],vertices[67],vertices[74]);
	  cells[32]->set_neighbors(cells[432],cells[63],cells[484],cells[470]);
	  set_offsets(cells[32],0,0,0,0);
	  cells[33]->set_vertices(vertices[64],vertices[24],vertices[17],vertices[72]);
	  cells[33]->set_neighbors(cells[287],cells[213],cells[195],cells[87]);
	  set_offsets(cells[33],0,0,0,0);
	  cells[34]->set_vertices(vertices[248],vertices[193],vertices[241],vertices[201]);
	  cells[34]->set_neighbors(cells[1502],cells[1509],cells[1481],cells[301]);
	  set_offsets(cells[34],0,0,0,0);
	  cells[35]->set_vertices(vertices[197],vertices[204],vertices[149],vertices[157]);
	  cells[35]->set_neighbors(cells[1190],cells[558],cells[443],cells[168]);
	  set_offsets(cells[35],0,0,0,0);
	  cells[36]->set_vertices(vertices[17],vertices[66],vertices[9],vertices[57]);
	  cells[36]->set_neighbors(cells[263],cells[270],cells[338],cells[232]);
	  set_offsets(cells[36],0,0,0,0);
	  cells[37]->set_vertices(vertices[204],vertices[196],vertices[149],vertices[156]);
	  cells[37]->set_neighbors(cells[852],cells[1190],cells[976],cells[168]);
	  set_offsets(cells[37],0,0,0,0);
	  cells[38]->set_vertices(vertices[31],vertices[80],vertices[72],vertices[79]);
	  cells[38]->set_neighbors(cells[684],cells[401],cells[251],cells[508]);
	  set_offsets(cells[38],0,1,1,0);
	  cells[39]->set_vertices(vertices[1],vertices[10],vertices[58],vertices[9]);
	  cells[39]->set_neighbors(cells[28],cells[166],cells[1507],cells[190]);
	  set_offsets(cells[39],0,0,0,0);
	  cells[40]->set_vertices(vertices[74],vertices[27],vertices[75],vertices[82]);
	  cells[40]->set_neighbors(cells[514],cells[697],cells[482],cells[515]);
	  set_offsets(cells[40],0,0,0,0);
	  cells[41]->set_vertices(vertices[43],vertices[92],vertices[44],vertices[35]);
	  cells[41]->set_neighbors(cells[143],cells[1004],cells[493],cells[523]);
	  set_offsets(cells[41],0,0,0,0);
	  cells[42]->set_vertices(vertices[60],vertices[52],vertices[5],vertices[12]);
	  cells[42]->set_neighbors(cells[337],cells[206],cells[120],cells[149]);
	  set_offsets(cells[42],0,0,0,0);
	  cells[43]->set_vertices(vertices[160],vertices[105],vertices[152],vertices[153]);
	  cells[43]->set_neighbors(cells[951],cells[1112],cells[955],cells[678]);
	  set_offsets(cells[43],0,0,0,0);
	  cells[44]->set_vertices(vertices[258],vertices[249],vertices[18],vertices[10]);
	  cells[44]->set_neighbors(cells[1006],cells[749],cells[1510],cells[145]);
	  set_offsets(cells[44],0,0,4,4);
	  cells[45]->set_vertices(vertices[116],vertices[69],vertices[68],vertices[61]);
	  cells[45]->set_neighbors(cells[488],cells[193],cells[317],cells[258]);
	  set_offsets(cells[45],0,0,0,0);
	  cells[46]->set_vertices(vertices[21],vertices[70],vertices[69],vertices[78]);
	  cells[46]->set_neighbors(cells[292],cells[452],cells[380],cells[201]);
	  set_offsets(cells[46],0,0,0,0);
	  cells[47]->set_vertices(vertices[177],vertices[169],vertices[226],vertices[178]);
	  cells[47]->set_neighbors(cells[4],cells[1039],cells[646],cells[1186]);
	  set_offsets(cells[47],0,0,0,0);
	  cells[48]->set_vertices(vertices[201],vertices[153],vertices[193],vertices[202]);
	  cells[48]->set_neighbors(cells[1149],cells[1501],cells[811],cells[1229]);
	  set_offsets(cells[48],0,0,0,0);
	  cells[49]->set_vertices(vertices[1],vertices[50],vertices[49],vertices[58]);
	  cells[49]->set_neighbors(cells[320],cells[166],cells[190],cells[583]);
	  set_offsets(cells[49],0,0,0,0);
	  cells[50]->set_vertices(vertices[93],vertices[140],vertices[92],vertices[52]);
	  cells[50]->set_neighbors(cells[349],cells[597],cells[788],cells[812]);
	  set_offsets(cells[50],0,0,0,2);
	  cells[51]->set_vertices(vertices[108],vertices[53],vertices[101],vertices[61]);
	  cells[51]->set_neighbors(cells[281],cells[706],cells[247],cells[382]);
	  set_offsets(cells[51],0,0,0,0);
	  cells[52]->set_vertices(vertices[41],vertices[90],vertices[2],vertices[42]);
	  cells[52]->set_neighbors(cells[477],cells[1606],cells[102],cells[564]);
	  set_offsets(cells[52],0,0,2,0);
	  cells[53]->set_vertices(vertices[1],vertices[2],vertices[50],vertices[10]);
	  cells[53]->set_neighbors(cells[72],cells[190],cells[1449],cells[328]);
	  set_offsets(cells[53],0,0,0,0);
	  cells[54]->set_vertices(vertices[98],vertices[51],vertices[91],vertices[139]);
	  cells[54]->set_neighbors(cells[823],cells[147],cells[872],cells[187]);
	  set_offsets(cells[54],2,2,0,0);
	  cells[55]->set_vertices(vertices[87],vertices[86],vertices[134],vertices[79]);
	  cells[55]->set_neighbors(cells[414],cells[808],cells[578],cells[218]);
	  set_offsets(cells[55],0,0,0,0);
	  cells[56]->set_vertices(vertices[172],vertices[115],vertices[164],vertices[124]);
	  cells[56]->set_neighbors(cells[970],cells[984],cells[845],cells[772]);
	  set_offsets(cells[56],0,0,0,0);
	  cells[57]->set_vertices(vertices[163],vertices[155],vertices[212],vertices[164]);
	  cells[57]->set_neighbors(cells[59],cells[833],cells[703],cells[441]);
	  set_offsets(cells[57],0,0,0,0);
	  cells[58]->set_vertices(vertices[285],vertices[277],vertices[44],vertices[37]);
	  cells[58]->set_neighbors(cells[1193],cells[1354],cells[1172],cells[1710]);
	  set_offsets(cells[58],0,0,4,4);
	  cells[59]->set_vertices(vertices[212],vertices[155],vertices[204],vertices[164]);
	  cells[59]->set_neighbors(cells[1030],cells[1258],cells[57],cells[1281]);
	  set_offsets(cells[59],0,0,0,0);
	  cells[60]->set_vertices(vertices[81],vertices[73],vertices[82],vertices[33]);
	  cells[60]->set_neighbors(cells[465],cells[290],cells[556],cells[601]);
	  set_offsets(cells[60],0,0,0,0);
	  cells[61]->set_vertices(vertices[13],vertices[61],vertices[53],vertices[60]);
	  cells[61]->set_neighbors(cells[247],cells[113],cells[274],cells[296]);
	  set_offsets(cells[61],0,0,0,0);
	  cells[62]->set_vertices(vertices[5],vertices[6],vertices[54],vertices[14]);
	  cells[62]->set_neighbors(cells[8],cells[321],cells[1417],cells[12]);
	  set_offsets(cells[62],0,0,0,0);
	  cells[63]->set_vertices(vertices[114],vertices[66],vertices[67],vertices[74]);
	  cells[63]->set_neighbors(cells[32],cells[397],cells[240],cells[82]);
	  set_offsets(cells[63],0,0,0,0);
	  cells[64]->set_vertices(vertices[100],vertices[108],vertices[99],vertices[51]);
	  cells[64]->set_neighbors(cells[210],cells[791],cells[658],cells[894]);
	  set_offsets(cells[64],0,0,0,0);
	  cells[65]->set_vertices(vertices[185],vertices[146],vertices[186],vertices[234]);
	  cells[65]->set_neighbors(cells[245],cells[1173],cells[1414],cells[1125]);
	  set_offsets(cells[65],0,2,0,0);
	  cells[66]->set_vertices(vertices[62],vertices[7],vertices[55],vertices[15]);
	  cells[66]->set_neighbors(cells[13],cells[21],cells[189],cells[114]);
	  set_offsets(cells[66],0,0,0,0);
	  cells[67]->set_vertices(vertices[144],vertices[96],vertices[137],vertices[97]);
	  cells[67]->set_neighbors(cells[861],cells[1129],cells[674],cells[1072]);
	  set_offsets(cells[67],2,2,0,2);
	  cells[68]->set_vertices(vertices[83],vertices[82],vertices[130],vertices[75]);
	  cells[68]->set_neighbors(cells[814],cells[781],cells[566],cells[388]);
	  set_offsets(cells[68],0,0,0,0);
	  cells[69]->set_vertices(vertices[175],vertices[216],vertices[167],vertices[215]);
	  cells[69]->set_neighbors(cells[919],cells[810],cells[837],cells[1313]);
	  set_offsets(cells[69],0,1,0,0);
	  cells[70]->set_vertices(vertices[48],vertices[88],vertices[89],vertices[41]);
	  cells[70]->set_neighbors(cells[378],cells[123],cells[567],cells[787]);
	  set_offsets(cells[70],2,0,0,0);
	  cells[71]->set_vertices(vertices[29],vertices[78],vertices[69],vertices[77]);
	  cells[71]->set_neighbors(cells[561],cells[528],cells[279],cells[452]);
	  set_offsets(cells[71],0,0,0,0);
	  cells[72]->set_vertices(vertices[50],vertices[2],vertices[3],vertices[10]);
	  cells[72]->set_neighbors(cells[1423],cells[306],cells[53],cells[374]);
	  set_offsets(cells[72],0,0,0,0);
	  cells[73]->set_vertices(vertices[153],vertices[105],vertices[145],vertices[154]);
	  cells[73]->set_neighbors(cells[557],cells[1151],cells[969],cells[951]);
	  set_offsets(cells[73],0,0,0,0);
	  cells[74]->set_vertices(vertices[33],vertices[82],vertices[34],vertices[25]);
	  cells[74]->set_neighbors(cells[520],cells[948],cells[465],cells[449]);
	  set_offsets(cells[74],0,0,0,0);
	  cells[75]->set_vertices(vertices[66],vertices[18],vertices[9],vertices[58]);
	  cells[75]->set_neighbors(cells[28],cells[263],cells[202],cells[232]);
	  set_offsets(cells[75],0,0,0,0);
	  cells[76]->set_vertices(vertices[173],vertices[220],vertices[213],vertices[165]);
	  cells[76]->set_neighbors(cells[1295],cells[1060],cells[1302],cells[1047]);
	  set_offsets(cells[76],0,0,0,0);
	  cells[77]->set_vertices(vertices[149],vertices[198],vertices[197],vertices[206]);
	  cells[77]->set_neighbors(cells[1213],cells[558],cells[151],cells[1182]);
	  set_offsets(cells[77],0,0,0,0);
	  cells[78]->set_vertices(vertices[130],vertices[138],vertices[129],vertices[81]);
	  cells[78]->set_neighbors(cells[381],cells[98],cells[591],cells[1017]);
	  set_offsets(cells[78],0,0,0,0);
	  cells[79]->set_vertices(vertices[13],vertices[14],vertices[62],vertices[22]);
	  cells[79]->set_neighbors(cells[18],cells[288],cells[1528],cells[278]);
	  set_offsets(cells[79],0,0,0,0);
	  cells[80]->set_vertices(vertices[75],vertices[76],vertices[124],vertices[84]);
	  cells[80]->set_neighbors(cells[545],cells[89],cells[458],cells[624]);
	  set_offsets(cells[80],0,0,0,0);
	  cells[81]->set_vertices(vertices[220],vertices[163],vertices[212],vertices[172]);
	  cells[81]->set_neighbors(cells[833],cells[1301],cells[233],cells[1062]);
	  set_offsets(cells[81],0,0,0,0);
	  cells[82]->set_vertices(vertices[114],vertices[66],vertices[59],vertices[67]);
	  cells[82]->set_neighbors(cells[470],cells[469],cells[63],cells[429]);
	  set_offsets(cells[82],0,0,0,0);
	  cells[83]->set_vertices(vertices[273],vertices[265],vertices[34],vertices[274]);
	  cells[83]->set_neighbors(cells[1563],cells[1691],cells[930],cells[1659]);
	  set_offsets(cells[83],0,0,4,0);
	  cells[84]->set_vertices(vertices[235],vertices[194],vertices[242],vertices[282]);
	  cells[84]->set_neighbors(cells[226],cells[1696],cells[1694],cells[1002]);
	  set_offsets(cells[84],0,2,2,0);
	  cells[85]->set_vertices(vertices[225],vertices[217],vertices[274],vertices[226]);
	  cells[85]->set_neighbors(cells[1178],cells[1655],cells[1187],cells[1330]);
	  set_offsets(cells[85],0,0,0,0);
	  cells[86]->set_vertices(vertices[76],vertices[28],vertices[19],vertices[68]);
	  cells[86]->set_neighbors(cells[271],cells[480],cells[481],cells[339]);
	  set_offsets(cells[86],0,0,0,0);
	  cells[87]->set_vertices(vertices[64],vertices[16],vertices[17],vertices[24]);
	  cells[87]->set_neighbors(cells[1553],cells[33],cells[297],cells[170]);
	  set_offsets(cells[87],0,0,0,0);
	  cells[88]->set_vertices(vertices[128],vertices[120],vertices[79],vertices[127]);
	  cells[88]->set_neighbors(cells[654],cells[737],cells[1038],cells[813]);
	  set_offsets(cells[88],1,1,0,0);
	  cells[89]->set_vertices(vertices[75],vertices[84],vertices[124],vertices[132]);
	  cells[89]->set_neighbors(cells[365],cells[418],cells[633],cells[80]);
	  set_offsets(cells[89],0,0,0,0);
	  cells[90]->set_vertices(vertices[25],vertices[74],vertices[17],vertices[65]);
	  cells[90]->set_neighbors(cells[486],cells[466],cells[358],cells[242]);
	  set_offsets(cells[90],0,0,0,0);
	  cells[91]->set_vertices(vertices[125],vertices[172],vertices[165],vertices[117]);
	  cells[91]->set_neighbors(cells[712],cells[982],cells[459],cells[953]);
	  set_offsets(cells[91],0,0,0,0);
	  cells[92]->set_vertices(vertices[72],vertices[24],vertices[25],vertices[32]);
	  cells[92]->set_neighbors(cells[1610],cells[437],cells[499],cells[287]);
	  set_offsets(cells[92],0,0,0,0);
	  cells[93]->set_vertices(vertices[47],vertices[40],vertices[88],vertices[0]);
	  cells[93]->set_neighbors(cells[235],cells[363],cells[1615],cells[550]);
	  set_offsets(cells[93],0,1,1,3);
	  cells[94]->set_vertices(vertices[73],vertices[122],vertices[113],vertices[121]);
	  cells[94]->set_neighbors(cells[792],cells[766],cells[115],cells[618]);
	  set_offsets(cells[94],0,0,0,0);
	  cells[95]->set_vertices(vertices[260],vertices[212],vertices[252],vertices[205]);
	  cells[95]->set_neighbors(cells[1243],cells[1465],cells[964],cells[1024]);
	  set_offsets(cells[95],0,0,0,0);
	  cells[96]->set_vertices(vertices[106],vertices[58],vertices[59],vertices[66]);
	  cells[96]->set_neighbors(cells[451],cells[429],cells[17],cells[681]);
	  set_offsets(cells[96],0,0,0,0);
	  cells[97]->set_vertices(vertices[1],vertices[56],vertices[49],vertices[48]);
	  cells[97]->set_neighbors(cells[604],cells[584],cells[335],cells[131]);
	  set_offsets(cells[97],0,0,0,0);
	  cells[98]->set_vertices(vertices[129],vertices[130],vertices[81],vertices[121]);
	  cells[98]->set_neighbors(cells[764],cells[617],cells[1065],cells[78]);
	  set_offsets(cells[98],0,0,0,0);
	  cells[99]->set_vertices(vertices[36],vertices[76],vertices[28],vertices[27]);
	  cells[99]->set_neighbors(cells[339],cells[1574],cells[220],cells[107]);
	  set_offsets(cells[99],0,0,0,0);
	  cells[100]->set_vertices(vertices[155],vertices[202],vertices[147],vertices[154]);
	  cells[100]->set_neighbors(cells[1191],cells[961],cells[1179],cells[1034]);
	  set_offsets(cells[100],0,0,0,0);
	  cells[101]->set_vertices(vertices[56],vertices[8],vertices[9],vertices[16]);
	  cells[101]->set_neighbors(cells[1505],cells[199],cells[385],cells[336]);
	  set_offsets(cells[101],0,0,0,0);
	  cells[102]->set_vertices(vertices[41],vertices[90],vertices[42],vertices[33]);
	  cells[102]->set_neighbors(cells[303],cells[1183],cells[560],cells[52]);
	  set_offsets(cells[102],0,0,0,0);
	  cells[103]->set_vertices(vertices[54],vertices[7],vertices[6],vertices[47]);
	  cells[103]->set_neighbors(cells[234],cells[426],cells[612],cells[8]);
	  set_offsets(cells[103],2,2,2,0);
	  cells[104]->set_vertices(vertices[28],vertices[20],vertices[68],vertices[21]);
	  cells[104]->set_neighbors(cells[434],cells[481],cells[1251],cells[271]);
	  set_offsets(cells[104],0,0,0,0);
	  cells[105]->set_vertices(vertices[3],vertices[12],vertices[60],vertices[11]);
	  cells[105]->set_neighbors(cells[179],cells[221],cells[803],cells[120]);
	  set_offsets(cells[105],0,0,0,0);
	  cells[106]->set_vertices(vertices[18],vertices[10],vertices[58],vertices[11]);
	  cells[106]->set_neighbors(cells[185],cells[202],cells[912],cells[28]);
	  set_offsets(cells[106],0,0,0,0);
	  cells[107]->set_vertices(vertices[36],vertices[29],vertices[28],vertices[76]);
	  cells[107]->set_neighbors(cells[439],cells[99],cells[502],cells[1527]);
	  set_offsets(cells[107],0,0,0,0);
	  cells[108]->set_vertices(vertices[5],vertices[93],vertices[45],vertices[52]);
	  cells[108]->set_neighbors(cells[597],cells[276],cells[599],cells[574]);
	  set_offsets(cells[108],2,0,0,2);
	  cells[109]->set_vertices(vertices[133],vertices[134],vertices[85],vertices[125]);
	  cells[109]->set_neighbors(cells[462],cells[828],cells[590],cells[389]);
	  set_offsets(cells[109],0,0,0,0);
	  cells[110]->set_vertices(vertices[253],vertices[252],vertices[245],vertices[205]);
	  cells[110]->set_neighbors(cells[1469],cells[1305],cells[1465],cells[1484]);
	  set_offsets(cells[110],0,0,0,0);
	  cells[111]->set_vertices(vertices[96],vertices[48],vertices[89],vertices[49]);
	  cells[111]->set_neighbors(cells[584],cells[858],cells[604],cells[350]);
	  set_offsets(cells[111],2,2,0,2);
	  cells[112]->set_vertices(vertices[7],vertices[8],vertices[56],vertices[15]);
	  cells[112]->set_neighbors(cells[385],cells[13],cells[1304],cells[188]);
	  set_offsets(cells[112],0,1,1,0);
	  cells[113]->set_vertices(vertices[5],vertices[13],vertices[53],vertices[60]);
	  cells[113]->set_neighbors(cells[61],cells[149],cells[206],cells[192]);
	  set_offsets(cells[113],0,0,0,0);
	  cells[114]->set_vertices(vertices[54],vertices[7],vertices[55],vertices[62]);
	  cells[114]->set_neighbors(cells[66],cells[371],cells[171],cells[611]);
	  set_offsets(cells[114],0,0,0,0);
	  cells[115]->set_vertices(vertices[73],vertices[122],vertices[121],vertices[130]);
	  cells[115]->set_neighbors(cells[521],cells[764],cells[778],cells[94]);
	  set_offsets(cells[115],0,0,0,0);
	  cells[116]->set_vertices(vertices[148],vertices[141],vertices[188],vertices[189]);
	  cells[116]->set_neighbors(cells[926],cells[1383],cells[390],cells[332]);
	  set_offsets(cells[116],2,0,0,0);
	  cells[117]->set_vertices(vertices[167],vertices[208],vertices[159],vertices[207]);
	  cells[117]->set_neighbors(cells[1141],cells[1161],cells[1204],cells[1268]);
	  set_offsets(cells[117],0,1,0,0);
	  cells[118]->set_vertices(vertices[82],vertices[35],vertices[90],vertices[42]);
	  cells[118]->set_neighbors(cells[313],cells[303],cells[497],cells[293]);
	  set_offsets(cells[118],0,0,0,0);
	  cells[119]->set_vertices(vertices[37],vertices[29],vertices[36],vertices[84]);
	  cells[119]->set_neighbors(cells[502],cells[368],cells[527],cells[1629]);
	  set_offsets(cells[119],0,0,0,0);
	  cells[120]->set_vertices(vertices[3],vertices[12],vertices[52],vertices[60]);
	  cells[120]->set_neighbors(cells[42],cells[223],cells[105],cells[31]);
	  set_offsets(cells[120],0,0,0,0);
	  cells[121]->set_vertices(vertices[55],vertices[104],vertices[63],vertices[56]);
	  cells[121]->set_neighbors(cells[209],cells[314],cells[150],cells[330]);
	  set_offsets(cells[121],0,1,0,1);
	  cells[122]->set_vertices(vertices[70],vertices[23],vertices[15],vertices[63]);
	  cells[122]->set_neighbors(cells[135],cells[341],cells[148],cells[447]);
	  set_offsets(cells[122],0,0,0,0);
	  cells[123]->set_vertices(vertices[1],vertices[89],vertices[41],vertices[48]);
	  cells[123]->set_neighbors(cells[70],cells[1],cells[584],cells[289]);
	  set_offsets(cells[123],2,0,0,2);
	  cells[124]->set_vertices(vertices[12],vertices[243],vertices[252],vertices[251]);
	  cells[124]->set_neighbors(cells[1257],cells[1576],cells[1468],cells[1379]);
	  set_offsets(cells[124],4,0,0,0);
	  cells[125]->set_vertices(vertices[280],vertices[225],vertices[233],vertices[232]);
	  cells[125]->set_neighbors(cells[1407],cells[1642],cells[1487],cells[1648]);
	  set_offsets(cells[125],0,0,0,0);
	  cells[126]->set_vertices(vertices[196],vertices[148],vertices[236],vertices[189]);
	  cells[126]->set_neighbors(cells[1383],cells[238],cells[862],cells[768]);
	  set_offsets(cells[126],2,2,0,0);
	  cells[127]->set_vertices(vertices[221],vertices[213],vertices[261],vertices[270]);
	  cells[127]->set_neighbors(cells[1585],cells[1536],cells[1586],cells[891]);
	  set_offsets(cells[127],0,0,0,0);
	  cells[128]->set_vertices(vertices[28],vertices[268],vertices[20],vertices[261]);
	  cells[128]->set_neighbors(cells[217],cells[1251],cells[1571],cells[1298]);
	  set_offsets(cells[128],4,0,4,0);
	  cells[129]->set_vertices(vertices[84],vertices[37],vertices[85],vertices[92]);
	  cells[129]->set_neighbors(cells[513],cells[632],cells[524],cells[568]);
	  set_offsets(cells[129],0,0,0,0);
	  cells[130]->set_vertices(vertices[99],vertices[148],vertices[147],vertices[156]);
	  cells[130]->set_neighbors(cells[1192],cells[153],cells[720],cells[1099]);
	  set_offsets(cells[130],0,0,0,0);
	  cells[131]->set_vertices(vertices[1],vertices[9],vertices[49],vertices[56]);
	  cells[131]->set_neighbors(cells[412],cells[97],cells[336],cells[166]);
	  set_offsets(cells[131],0,0,0,0);
	  cells[132]->set_vertices(vertices[103],vertices[55],vertices[96],vertices[104]);
	  cells[132]->set_neighbors(cells[150],cells[638],cells[330],cells[838]);
	  set_offsets(cells[132],0,0,1,1);
	  cells[133]->set_vertices(vertices[236],vertices[188],vertices[179],vertices[228]);
	  cells[133]->set_neighbors(cells[925],cells[269],cells[327],cells[1392]);
	  set_offsets(cells[133],0,0,0,0);
	  cells[134]->set_vertices(vertices[199],vertices[192],vertices[240],vertices[200]);
	  cells[134]->set_neighbors(cells[1143],cells[1498],cells[1221],cells[1444]);
	  set_offsets(cells[134],0,1,1,1);
	  cells[135]->set_vertices(vertices[15],vertices[23],vertices[64],vertices[63]);
	  cells[135]->set_neighbors(cells[169],cells[141],cells[122],cells[453]);
	  set_offsets(cells[135],0,0,1,0);
	  cells[136]->set_vertices(vertices[21],vertices[261],vertices[269],vertices[28]);
	  cells[136]->set_neighbors(cells[1571],cells[1464],cells[1251],cells[1209]);
	  set_offsets(cells[136],4,0,0,4);
	  cells[137]->set_vertices(vertices[159],vertices[200],vertices[199],vertices[207]);
	  cells[137]->set_neighbors(cells[1450],cells[1205],cells[1141],cells[1222]);
	  set_offsets(cells[137],0,1,0,0);
	  cells[138]->set_vertices(vertices[209],vertices[161],vertices[208],vertices[201]);
	  cells[138]->set_neighbors(cells[25],cells[176],cells[960],cells[1273]);
	  set_offsets(cells[138],0,0,0,0);
	  cells[139]->set_vertices(vertices[121],vertices[113],vertices[161],vertices[170]);
	  cells[139]->set_neighbors(cells[954],cells[444],cells[792],cells[992]);
	  set_offsets(cells[139],0,0,0,0);
	  cells[140]->set_vertices(vertices[135],vertices[142],vertices[134],vertices[87]);
	  cells[140]->set_neighbors(cells[850],cells[822],cells[848],cells[1106]);
	  set_offsets(cells[140],0,0,0,0);
	  cells[141]->set_vertices(vertices[15],vertices[64],vertices[56],vertices[63]);
	  cells[141]->set_neighbors(cells[209],cells[314],cells[135],cells[184]);
	  set_offsets(cells[141],0,1,1,0);
	  cells[142]->set_vertices(vertices[94],vertices[86],vertices[85],vertices[134]);
	  cells[142]->set_neighbors(cells[311],cells[880],cells[218],cells[407]);
	  set_offsets(cells[142],0,0,0,0);
	  cells[143]->set_vertices(vertices[92],vertices[84],vertices[44],vertices[35]);
	  cells[143]->set_neighbors(cells[538],cells[41],cells[302],cells[524]);
	  set_offsets(cells[143],0,0,0,0);
	  cells[144]->set_vertices(vertices[19],vertices[20],vertices[11],vertices[68]);
	  cells[144]->set_neighbors(cells[11],cells[433],cells[271],cells[1567]);
	  set_offsets(cells[144],0,0,0,0);
	  cells[145]->set_vertices(vertices[257],vertices[249],vertices[18],vertices[258]);
	  cells[145]->set_neighbors(cells[44],cells[1558],cells[1378],cells[1556]);
	  set_offsets(cells[145],0,0,4,0);
	  cells[146]->set_vertices(vertices[206],vertices[198],vertices[151],vertices[158]);
	  cells[146]->set_neighbors(cells[641],cells[1111],cells[151],cells[1240]);
	  set_offsets(cells[146],0,0,0,0);
	  cells[147]->set_vertices(vertices[98],vertices[91],vertices[138],vertices[139]);
	  cells[147]->set_neighbors(cells[376],cells[413],cells[54],cells[696]);
	  set_offsets(cells[147],2,0,0,0);
	  cells[148]->set_vertices(vertices[70],vertices[71],vertices[23],vertices[63]);
	  cells[148]->set_neighbors(cells[169],cells[122],cells[294],cells[355]);
	  set_offsets(cells[148],0,0,0,0);
	  cells[149]->set_vertices(vertices[5],vertices[60],vertices[53],vertices[52]);
	  cells[149]->set_neighbors(cells[650],cells[599],cells[42],cells[113]);
	  set_offsets(cells[149],0,0,0,0);
	  cells[150]->set_vertices(vertices[55],vertices[96],vertices[104],vertices[56]);
	  cells[150]->set_neighbors(cells[637],cells[121],cells[200],cells[132]);
	  set_offsets(cells[150],0,1,1,1);
	  cells[151]->set_vertices(vertices[149],vertices[158],vertices[198],vertices[206]);
	  cells[151]->set_neighbors(cells[146],cells[77],cells[692],cells[643]);
	  set_offsets(cells[151],0,0,0,0);
	  cells[152]->set_vertices(vertices[91],vertices[90],vertices[83],vertices[43]);
	  cells[152]->set_neighbors(cells[329],cells[211],cells[156],cells[820]);
	  set_offsets(cells[152],0,0,0,0);
	  cells[153]->set_vertices(vertices[99],vertices[156],vertices[147],vertices[107]);
	  cells[153]->set_neighbors(cells[843],cells[704],cells[907],cells[130]);
	  set_offsets(cells[153],0,0,0,0);
	  cells[154]->set_vertices(vertices[116],vertices[67],vertices[107],vertices[59]);
	  cells[154]->set_neighbors(cells[469],cells[478],cells[549],cells[286]);
	  set_offsets(cells[154],0,0,0,0);
	  cells[155]->set_vertices(vertices[180],vertices[220],vertices[173],vertices[172]);
	  cells[155]->set_neighbors(cells[1302],cells[993],cells[870],cells[1345]);
	  set_offsets(cells[155],0,0,0,0);
	  cells[156]->set_vertices(vertices[50],vertices[90],vertices[91],vertices[43]);
	  cells[156]->set_neighbors(cells[152],cells[519],cells[536],cells[14]);
	  set_offsets(cells[156],2,0,0,0);
	  cells[157]->set_vertices(vertices[224],vertices[216],vertices[264],vertices[217]);
	  cells[157]->set_neighbors(cells[1446],cells[1315],cells[1320],cells[1605]);
	  set_offsets(cells[157],0,0,0,0);
	  cells[158]->set_vertices(vertices[161],vertices[153],vertices[210],vertices[162]);
	  cells[158]->set_neighbors(cells[208],cells[1008],cells[746],cells[1027]);
	  set_offsets(cells[158],0,0,0,0);
	  cells[159]->set_vertices(vertices[127],vertices[135],vertices[182],vertices[134]);
	  cells[159]->set_neighbors(cells[1106],cells[780],cells[822],cells[1076]);
	  set_offsets(cells[159],0,0,0,0);
	  cells[160]->set_vertices(vertices[80],vertices[32],vertices[33],vertices[40]);
	  cells[160]->set_neighbors(cells[1439],cells[562],cells[539],cells[402]);
	  set_offsets(cells[160],0,0,0,0);
	  cells[161]->set_vertices(vertices[138],vertices[89],vertices[137],vertices[129]);
	  cells[161]->set_neighbors(cells[815],cells[652],cells[381],cells[856]);
	  set_offsets(cells[161],0,0,0,0);
	  cells[162]->set_vertices(vertices[100],vertices[51],vertices[91],vertices[52]);
	  cells[162]->set_neighbors(cells[592],cells[492],cells[372],cells[823]);
	  set_offsets(cells[162],2,2,0,2);
	  cells[163]->set_vertices(vertices[39],vertices[40],vertices[80],vertices[88]);
	  cells[163]->set_neighbors(cells[562],cells[546],cells[550],cells[539]);
	  set_offsets(cells[163],0,1,1,1);
	  cells[164]->set_vertices(vertices[157],vertices[109],vertices[149],vertices[158]);
	  cells[164]->set_neighbors(cells[691],cells[692],cells[182],cells[971]);
	  set_offsets(cells[164],0,0,0,0);
	  cells[165]->set_vertices(vertices[169],vertices[161],vertices[218],vertices[170]);
	  cells[165]->set_neighbors(cells[248],cells[959],cells[444],cells[1049]);
	  set_offsets(cells[165],0,0,0,0);
	  cells[166]->set_vertices(vertices[1],vertices[58],vertices[49],vertices[9]);
	  cells[166]->set_neighbors(cells[177],cells[131],cells[39],cells[49]);
	  set_offsets(cells[166],0,0,0,0);
	  cells[167]->set_vertices(vertices[92],vertices[52],vertices[91],vertices[43]);
	  cells[167]->set_neighbors(cells[563],cells[211],cells[525],cells[349]);
	  set_offsets(cells[167],0,2,0,0);
	  cells[168]->set_vertices(vertices[197],vertices[196],vertices[149],vertices[204]);
	  cells[168]->set_neighbors(cells[37],cells[35],cells[1244],cells[1199]);
	  set_offsets(cells[168],0,0,0,0);
	  cells[169]->set_vertices(vertices[23],vertices[71],vertices[64],vertices[63]);
	  cells[169]->set_neighbors(cells[725],cells[135],cells[148],cells[203]);
	  set_offsets(cells[169],0,0,1,0);
	  cells[170]->set_vertices(vertices[64],vertices[16],vertices[9],vertices[17]);
	  cells[170]->set_neighbors(cells[1554],cells[270],cells[87],cells[199]);
	  set_offsets(cells[170],0,0,0,0);
	  cells[171]->set_vertices(vertices[54],vertices[14],vertices[7],vertices[62]);
	  cells[171]->set_neighbors(cells[189],cells[114],cells[321],cells[8]);
	  set_offsets(cells[171],0,0,0,0);
	  cells[172]->set_vertices(vertices[222],vertices[213],vertices[262],vertices[214]);
	  cells[172]->set_neighbors(cells[1290],cells[1349],cells[1338],cells[1057]);
	  set_offsets(cells[172],0,0,0,0);
	  cells[173]->set_vertices(vertices[247],vertices[254],vertices[199],vertices[207]);
	  cells[173]->set_neighbors(cells[1249],cells[1441],cells[1291],cells[769]);
	  set_offsets(cells[173],0,0,0,0);
	  cells[174]->set_vertices(vertices[128],vertices[136],vertices[135],vertices[87]);
	  cells[174]->set_neighbors(cells[531],cells[758],cells[273],cells[782]);
	  set_offsets(cells[174],1,1,0,0);
	  cells[175]->set_vertices(vertices[45],vertices[54],vertices[6],vertices[94]);
	  cells[175]->set_neighbors(cells[426],cells[581],cells[255],cells[12]);
	  set_offsets(cells[175],0,2,2,0);
	  cells[176]->set_vertices(vertices[209],vertices[201],vertices[208],vertices[256]);
	  cells[176]->set_neighbors(cells[940],cells[1226],cells[1542],cells[138]);
	  set_offsets(cells[176],0,0,0,0);
	  cells[177]->set_vertices(vertices[58],vertices[57],vertices[49],vertices[9]);
	  cells[177]->set_neighbors(cells[412],cells[166],cells[263],cells[431]);
	  set_offsets(cells[177],0,0,0,0);
	  cells[178]->set_vertices(vertices[13],vertices[22],vertices[70],vertices[21]);
	  cells[178]->set_neighbors(cells[406],cells[367],cells[1580],cells[288]);
	  set_offsets(cells[178],0,0,0,0);
	  cells[179]->set_vertices(vertices[20],vertices[12],vertices[11],vertices[60]);
	  cells[179]->set_neighbors(cells[105],cells[11],cells[319],cells[183]);
	  set_offsets(cells[179],0,0,0,0);
	  cells[180]->set_vertices(vertices[112],vertices[72],vertices[64],vertices[65]);
	  cells[180]->set_neighbors(cells[213],cells[569],cells[456],cells[724]);
	  set_offsets(cells[180],0,0,0,0);
	  cells[181]->set_vertices(vertices[67],vertices[76],vertices[68],vertices[116]);
	  cells[181]->set_neighbors(cells[258],cells[549],cells[619],cells[480]);
	  set_offsets(cells[181],0,0,0,0);
	  cells[182]->set_vertices(vertices[166],vertices[109],vertices[157],vertices[158]);
	  cells[182]->set_neighbors(cells[164],cells[1041],cells[29],cells[268]);
	  set_offsets(cells[182],0,0,0,0);
	  cells[183]->set_vertices(vertices[251],vertices[12],vertices[11],vertices[20]);
	  cells[183]->set_neighbors(cells[179],cells[1566],cells[1253],cells[803]);
	  set_offsets(cells[183],0,4,4,4);
	  cells[184]->set_vertices(vertices[15],vertices[16],vertices[56],vertices[64]);
	  cells[184]->set_neighbors(cells[199],cells[141],cells[453],cells[385]);
	  set_offsets(cells[184],0,1,1,1);
	  cells[185]->set_vertices(vertices[10],vertices[11],vertices[3],vertices[58]);
	  cells[185]->set_neighbors(cells[15],cells[306],cells[106],cells[9]);
	  set_offsets(cells[185],0,0,0,0);
	  cells[186]->set_vertices(vertices[34],vertices[27],vertices[26],vertices[74]);
	  cells[186]->set_neighbors(cells[285],cells[30],cells[482],cells[1565]);
	  set_offsets(cells[186],0,0,0,0);
	  cells[187]->set_vertices(vertices[51],vertices[98],vertices[91],vertices[50]);
	  cells[187]->set_neighbors(cells[696],cells[473],cells[448],cells[54]);
	  set_offsets(cells[187],2,2,0,2);
	  cells[188]->set_vertices(vertices[7],vertices[8],vertices[48],vertices[56]);
	  cells[188]->set_neighbors(cells[335],cells[230],cells[112],cells[10]);
	  set_offsets(cells[188],0,1,1,1);
	  cells[189]->set_vertices(vertices[62],vertices[14],vertices[7],vertices[15]);
	  cells[189]->set_neighbors(cells[1495],cells[66],cells[18],cells[171]);
	  set_offsets(cells[189],0,0,0,0);
	  cells[190]->set_vertices(vertices[1],vertices[10],vertices[50],vertices[58]);
	  cells[190]->set_neighbors(cells[306],cells[49],cells[39],cells[53]);
	  set_offsets(cells[190],0,0,0,0);
	  cells[191]->set_vertices(vertices[117],vertices[164],vertices[157],vertices[109]);
	  cells[191]->set_neighbors(cells[647],cells[268],cells[642],cells[1018]);
	  set_offsets(cells[191],0,0,0,0);
	  cells[192]->set_vertices(vertices[5],vertices[62],vertices[53],vertices[13]);
	  cells[192]->set_neighbors(cells[296],cells[113],cells[278],cells[207]);
	  set_offsets(cells[192],0,0,0,0);
	  cells[193]->set_vertices(vertices[108],vertices[116],vertices[68],vertices[61]);
	  cells[193]->set_neighbors(cells[45],cells[479],cells[700],cells[6]);
	  set_offsets(cells[193],0,0,0,0);
	  cells[194]->set_vertices(vertices[181],vertices[173],vertices[230],vertices[182]);
	  cells[194]->set_neighbors(cells[244],cells[1165],cells[986],cells[1109]);
	  set_offsets(cells[194],0,0,0,0);
	  cells[195]->set_vertices(vertices[23],vertices[24],vertices[64],vertices[72]);
	  cells[195]->set_neighbors(cells[33],cells[203],cells[2],cells[297]);
	  set_offsets(cells[195],0,1,1,1);
	  cells[196]->set_vertices(vertices[127],vertices[174],vertices[119],vertices[126]);
	  cells[196]->set_neighbors(cells[785],cells[759],cells[710],cells[869]);
	  set_offsets(cells[196],0,0,0,0);
	  cells[197]->set_vertices(vertices[6],vertices[245],vertices[246],vertices[254]);
	  cells[197]->set_neighbors(cells[1026],cells[1485],cells[1466],cells[1056]);
	  set_offsets(cells[197],4,0,0,0);
	  cells[198]->set_vertices(vertices[81],vertices[90],vertices[89],vertices[41]);
	  cells[198]->set_neighbors(cells[399],cells[378],cells[560],cells[789]);
	  set_offsets(cells[198],0,0,0,0);
	  cells[199]->set_vertices(vertices[56],vertices[16],vertices[9],vertices[64]);
	  cells[199]->set_neighbors(cells[170],cells[420],cells[184],cells[101]);
	  set_offsets(cells[199],0,0,0,0);
	  cells[200]->set_vertices(vertices[55],vertices[56],vertices[48],vertices[96]);
	  cells[200]->set_neighbors(cells[604],cells[576],cells[150],cells[230]);
	  set_offsets(cells[200],0,1,1,1);
	  cells[201]->set_vertices(vertices[21],vertices[70],vertices[61],vertices[69]);
	  cells[201]->set_neighbors(cells[722],cells[488],cells[46],cells[367]);
	  set_offsets(cells[201],0,0,0,0);
	  cells[202]->set_vertices(vertices[66],vertices[58],vertices[11],vertices[18]);
	  cells[202]->set_neighbors(cells[106],cells[423],cells[75],cells[451]);
	  set_offsets(cells[202],0,0,0,0);
	  cells[203]->set_vertices(vertices[23],vertices[72],vertices[64],vertices[71]);
	  cells[203]->set_neighbors(cells[724],cells[169],cells[428],cells[195]);
	  set_offsets(cells[203],0,1,1,0);
	  cells[204]->set_vertices(vertices[61],vertices[118],vertices[70],vertices[110]);
	  cells[204]->set_neighbors(cells[629],cells[677],cells[284],cells[722]);
	  set_offsets(cells[204],0,0,0,0);
	  cells[205]->set_vertices(vertices[76],vertices[75],vertices[67],vertices[27]);
	  cells[205]->set_neighbors(cells[515],cells[283],cells[458],cells[624]);
	  set_offsets(cells[205],0,0,0,0);
	  cells[206]->set_vertices(vertices[13],vertices[60],vertices[5],vertices[12]);
	  cells[206]->set_neighbors(cells[42],cells[1532],cells[319],cells[113]);
	  set_offsets(cells[206],0,0,0,0);
	  cells[207]->set_vertices(vertices[5],vertices[54],vertices[53],vertices[62]);
	  cells[207]->set_neighbors(cells[411],cells[192],cells[321],cells[600]);
	  set_offsets(cells[207],0,0,0,0);
	  cells[208]->set_vertices(vertices[210],vertices[153],vertices[202],vertices[162]);
	  cells[208]->set_neighbors(cells[783],cells[250],cells[158],cells[811]);
	  set_offsets(cells[208],0,0,0,0);
	  cells[209]->set_vertices(vertices[63],vertices[64],vertices[56],vertices[104]);
	  cells[209]->set_neighbors(cells[215],cells[121],cells[726],cells[141]);
	  set_offsets(cells[209],0,1,1,1);
	  cells[210]->set_vertices(vertices[108],vertices[59],vertices[99],vertices[51]);
	  cells[210]->set_neighbors(cells[635],cells[64],cells[7],cells[630]);
	  set_offsets(cells[210],0,0,0,0);
	  cells[211]->set_vertices(vertices[83],vertices[92],vertices[91],vertices[43]);
	  cells[211]->set_neighbors(cells[167],cells[152],cells[493],cells[416]);
	  set_offsets(cells[211],0,0,0,0);
	  cells[212]->set_vertices(vertices[57],vertices[106],vertices[105],vertices[114]);
	  cells[212]->set_neighbors(cells[950],cells[450],cells[679],cells[304]);
	  set_offsets(cells[212],0,0,0,0);
	  cells[213]->set_vertices(vertices[64],vertices[72],vertices[17],vertices[65]);
	  cells[213]->set_neighbors(cells[466],cells[454],cells[180],cells[33]);
	  set_offsets(cells[213],0,0,0,0);
	  cells[214]->set_vertices(vertices[266],vertices[258],vertices[18],vertices[259]);
	  cells[214]->set_neighbors(cells[421],cells[1121],cells[1564],cells[1558]);
	  set_offsets(cells[214],0,0,4,0);
	  cells[215]->set_vertices(vertices[104],vertices[64],vertices[56],vertices[57]);
	  cells[215]->set_neighbors(cells[420],cells[395],cells[509],cells[209]);
	  set_offsets(cells[215],0,0,0,0);
	  cells[216]->set_vertices(vertices[168],vertices[161],vertices[208],vertices[216]);
	  cells[216]->set_neighbors(cells[1273],cells[1312],cells[994],cells[1087]);
	  set_offsets(cells[216],0,0,0,0);
	  cells[217]->set_vertices(vertices[268],vertices[260],vertices[20],vertices[261]);
	  cells[217]->set_neighbors(cells[1188],cells[128],cells[522],cells[1023]);
	  set_offsets(cells[217],0,0,4,0);
	  cells[218]->set_vertices(vertices[94],vertices[86],vertices[134],vertices[87]);
	  cells[218]->set_neighbors(cells[55],cells[850],cells[596],cells[142]);
	  set_offsets(cells[218],0,0,0,0);
	  cells[219]->set_vertices(vertices[45],vertices[94],vertices[46],vertices[37]);
	  cells[219]->set_neighbors(cells[555],cells[1219],cells[573],cells[581]);
	  set_offsets(cells[219],0,0,0,0);
	  cells[220]->set_vertices(vertices[84],vertices[76],vertices[36],vertices[27]);
	  cells[220]->set_neighbors(cells[99],cells[503],cells[458],cells[502]);
	  set_offsets(cells[220],0,0,0,0);
	  cells[221]->set_vertices(vertices[3],vertices[60],vertices[51],vertices[11]);
	  cells[221]->set_neighbors(cells[3],cells[15],cells[105],cells[223]);
	  set_offsets(cells[221],0,0,0,0);
	  cells[222]->set_vertices(vertices[14],vertices[253],vertices[262],vertices[22]);
	  cells[222]->set_neighbors(cells[1411],cells[875],cells[1528],cells[1224]);
	  set_offsets(cells[222],4,0,0,4);
	  cells[223]->set_vertices(vertices[3],vertices[52],vertices[51],vertices[60]);
	  cells[223]->set_neighbors(cells[372],cells[221],cells[120],cells[592]);
	  set_offsets(cells[223],0,0,0,0);
	  cells[224]->set_vertices(vertices[117],vertices[109],vertices[166],vertices[118]);
	  cells[224]->set_neighbors(cells[29],cells[844],cells[580],cells[268]);
	  set_offsets(cells[224],0,0,0,0);
	  cells[225]->set_vertices(vertices[118],vertices[110],vertices[111],vertices[63]);
	  cells[225]->set_neighbors(cells[717],cells[430],cells[629],cells[354]);
	  set_offsets(cells[225],0,0,0,0);
	  cells[226]->set_vertices(vertices[242],vertices[194],vertices[233],vertices[282]);
	  cells[226]->set_neighbors(cells[1617],cells[1679],cells[84],cells[1455]);
	  set_offsets(cells[226],2,2,0,0);
	  cells[227]->set_vertices(vertices[112],vertices[57],vertices[104],vertices[105]);
	  cells[227]->set_neighbors(cells[683],cells[842],cells[510],cells[509]);
	  set_offsets(cells[227],0,0,0,0);
	  cells[228]->set_vertices(vertices[120],vertices[80],vertices[72],vertices[73]);
	  cells[228]->set_neighbors(cells[305],cells[738],cells[735],cells[684]);
	  set_offsets(cells[228],0,0,0,0);
	  cells[229]->set_vertices(vertices[0],vertices[240],vertices[241],vertices[248]);
	  cells[229]->set_neighbors(cells[301],cells[975],cells[1496],cells[1684]);
	  set_offsets(cells[229],4,0,0,0);
	  cells[230]->set_vertices(vertices[55],vertices[7],vertices[48],vertices[56]);
	  cells[230]->set_neighbors(cells[188],cells[200],cells[13],cells[610]);
	  set_offsets(cells[230],0,0,1,1);
	  cells[231]->set_vertices(vertices[65],vertices[114],vertices[113],vertices[122]);
	  cells[231]->set_neighbors(cells[748],cells[618],cells[353],cells[504]);
	  set_offsets(cells[231],0,0,0,0);
	  cells[232]->set_vertices(vertices[17],vertices[18],vertices[9],vertices[66]);
	  cells[232]->set_neighbors(cells[75],cells[36],cells[361],cells[1555]);
	  set_offsets(cells[232],0,0,0,0);
	  cells[233]->set_vertices(vertices[171],vertices[163],vertices[220],vertices[172]);
	  cells[233]->set_neighbors(cells[81],cells[870],cells[917],cells[1322]);
	  set_offsets(cells[233],0,0,0,0);
	  cells[234]->set_vertices(vertices[6],vertices[247],vertices[47],vertices[7]);
	  cells[234]->set_neighbors(cells[1353],cells[103],cells[266],cells[1724]);
	  set_offsets(cells[234],6,2,4,6);
	  cells[235]->set_vertices(vertices[88],vertices[40],vertices[41],vertices[0]);
	  cells[235]->set_neighbors(cells[1686],cells[567],cells[93],cells[507]);
	  set_offsets(cells[235],0,0,0,2);
	  cells[236]->set_vertices(vertices[4],vertices[245],vertices[45],vertices[5]);
	  cells[236]->set_neighbors(cells[1425],cells[276],cells[966],cells[1207]);
	  set_offsets(cells[236],6,2,4,6);
	  cells[237]->set_vertices(vertices[244],vertices[4],vertices[44],vertices[283]);
	  cells[237]->set_neighbors(cells[1364],cells[324],cells[1283],cells[1711]);
	  set_offsets(cells[237],2,6,4,0);
	  cells[238]->set_vertices(vertices[237],vertices[196],vertices[236],vertices[189]);
	  cells[238]->set_neighbors(cells[126],cells[1427],cells[1429],cells[1706]);
	  set_offsets(cells[238],0,2,0,0);
	  cells[239]->set_vertices(vertices[13],vertices[60],vertices[20],vertices[68]);
	  cells[239]->set_neighbors(cells[11],cells[434],cells[274],cells[319]);
	  set_offsets(cells[239],0,0,0,0);
	  cells[240]->set_vertices(vertices[65],vertices[66],vertices[114],vertices[74]);
	  cells[240]->set_neighbors(cells[63],cells[353],cells[486],cells[257]);
	  set_offsets(cells[240],0,0,0,0);
	  cells[241]->set_vertices(vertices[159],vertices[152],vertices[111],vertices[151]);
	  cells[241]->set_neighbors(cells[898],cells[585],cells[941],cells[767]);
	  set_offsets(cells[241],0,1,0,0);
	  cells[242]->set_vertices(vertices[25],vertices[26],vertices[17],vertices[74]);
	  cells[242]->set_neighbors(cells[485],cells[90],cells[30],cells[1612]);
	  set_offsets(cells[242],0,0,0,0);
	  cells[243]->set_vertices(vertices[13],vertices[21],vertices[61],vertices[68]);
	  cells[243]->set_neighbors(cells[488],cells[274],cells[434],cells[367]);
	  set_offsets(cells[243],0,0,0,0);
	  cells[244]->set_vertices(vertices[230],vertices[173],vertices[222],vertices[182]);
	  cells[244]->set_neighbors(cells[1071],cells[1357],cells[194],cells[1116]);
	  set_offsets(cells[244],0,0,0,0);
	  cells[245]->set_vertices(vertices[234],vertices[146],vertices[186],vertices[187]);
	  cells[245]->set_neighbors(cells[1138],cells[1171],cells[1115],cells[65]);
	  set_offsets(cells[245],0,2,0,0);
	  cells[246]->set_vertices(vertices[173],vertices[165],vertices[222],vertices[174]);
	  cells[246]->set_neighbors(cells[256],cells[1071],cells[540],cells[1060]);
	  set_offsets(cells[246],0,0,0,0);
	  cells[247]->set_vertices(vertices[60],vertices[61],vertices[53],vertices[108]);
	  cells[247]->set_neighbors(cells[51],cells[438],cells[479],cells[61]);
	  set_offsets(cells[247],0,0,0,0);
	  cells[248]->set_vertices(vertices[218],vertices[161],vertices[210],vertices[170]);
	  cells[248]->set_neighbors(cells[1008],cells[1289],cells[165],cells[979]);
	  set_offsets(cells[248],0,0,0,0);
	  cells[249]->set_vertices(vertices[3],vertices[58],vertices[51],vertices[50]);
	  cells[249]->set_neighbors(cells[448],cells[473],cells[306],cells[15]);
	  set_offsets(cells[249],0,0,0,0);
	  cells[250]->set_vertices(vertices[210],vertices[202],vertices[155],vertices[162]);
	  cells[250]->set_neighbors(cells[1179],cells[1010],cells[208],cells[1236]);
	  set_offsets(cells[250],0,0,0,0);
	  cells[251]->set_vertices(vertices[31],vertices[39],vertices[80],vertices[79]);
	  cells[251]->set_neighbors(cells[577],cells[38],cells[387],cells[495]);
	  set_offsets(cells[251],0,0,1,0);
	  cells[252]->set_vertices(vertices[282],vertices[235],vertices[234],vertices[227]);
	  cells[252]->set_neighbors(cells[1418],cells[1650],cells[1452],cells[1694]);
	  set_offsets(cells[252],0,0,0,0);
	  cells[253]->set_vertices(vertices[158],vertices[150],vertices[103],vertices[110]);
	  cells[253]->set_neighbors(cells[932],cells[702],cells[740],cells[670]);
	  set_offsets(cells[253],0,0,0,0);
	  cells[254]->set_vertices(vertices[176],vertices[169],vertices[216],vertices[224]);
	  cells[254]->set_neighbors(cells[1320],cells[1360],cells[1358],cells[16]);
	  set_offsets(cells[254],0,0,0,0);
	  cells[255]->set_vertices(vertices[94],vertices[54],vertices[93],vertices[45]);
	  cells[255]->set_neighbors(cells[574],cells[532],cells[175],cells[881]);
	  set_offsets(cells[255],0,2,0,0);
	  cells[256]->set_vertices(vertices[222],vertices[165],vertices[214],vertices[174]);
	  cells[256]->set_neighbors(cells[918],cells[1311],cells[246],cells[1338]);
	  set_offsets(cells[256],0,0,0,0);
	  cells[257]->set_vertices(vertices[57],vertices[66],vertices[114],vertices[65]);
	  cells[257]->set_neighbors(cells[240],cells[450],cells[338],cells[679]);
	  set_offsets(cells[257],0,0,0,0);
	  cells[258]->set_vertices(vertices[116],vertices[76],vertices[68],vertices[69]);
	  cells[258]->set_neighbors(cells[22],cells[45],cells[375],cells[181]);
	  set_offsets(cells[258],0,0,0,0);
	  cells[259]->set_vertices(vertices[228],vertices[171],vertices[220],vertices[180]);
	  cells[259]->set_neighbors(cells[870],cells[1345],cells[261],cells[962]);
	  set_offsets(cells[259],0,0,0,0);
	  cells[260]->set_vertices(vertices[58],vertices[106],vertices[98],vertices[49]);
	  cells[260]->set_neighbors(cells[483],cells[320],cells[431],cells[312]);
	  set_offsets(cells[260],0,0,0,0);
	  cells[261]->set_vertices(vertices[179],vertices[171],vertices[228],vertices[180]);
	  cells[261]->set_neighbors(cells[259],cells[925],cells[1095],cells[1122]);
	  set_offsets(cells[261],0,0,0,0);
	  cells[262]->set_vertices(vertices[268],vertices[220],vertices[260],vertices[213]);
	  cells[262]->set_neighbors(cells[1284],cells[522],cells[1381],cells[1365]);
	  set_offsets(cells[262],0,0,0,0);
	  cells[263]->set_vertices(vertices[66],vertices[58],vertices[9],vertices[57]);
	  cells[263]->set_neighbors(cells[177],cells[36],cells[17],cells[75]);
	  set_offsets(cells[263],0,0,0,0);
	  cells[264]->set_vertices(vertices[270],vertices[262],vertices[22],vertices[263]);
	  cells[264]->set_neighbors(cells[1025],cells[1202],cells[779],cells[1533]);
	  set_offsets(cells[264],0,0,4,0);
	  cells[265]->set_vertices(vertices[212],vertices[204],vertices[205],vertices[157]);
	  cells[265]->set_neighbors(cells[443],cells[1255],cells[1258],cells[1243]);
	  set_offsets(cells[265],0,0,0,0);
	  cells[266]->set_vertices(vertices[6],vertices[14],vertices[247],vertices[7]);
	  cells[266]->set_neighbors(cells[708],cells[234],cells[8],cells[1239]);
	  set_offsets(cells[266],4,4,0,4);
	  cells[267]->set_vertices(vertices[275],vertices[44],vertices[36],vertices[35]);
	  cells[267]->set_neighbors(cells[538],cells[1424],cells[1666],cells[1705]);
	  set_offsets(cells[267],0,4,4,4);
	  cells[268]->set_vertices(vertices[117],vertices[109],vertices[157],vertices[166]);
	  cells[268]->set_neighbors(cells[182],cells[980],cells[224],cells[191]);
	  set_offsets(cells[268],0,0,0,0);
	  cells[269]->set_vertices(vertices[227],vertices[236],vertices[179],vertices[228]);
	  cells[269]->set_neighbors(cells[133],cells[1201],cells[1362],cells[1200]);
	  set_offsets(cells[269],0,0,0,0);
	  cells[270]->set_vertices(vertices[64],vertices[17],vertices[9],vertices[57]);
	  cells[270]->set_neighbors(cells[36],cells[420],cells[454],cells[170]);
	  set_offsets(cells[270],0,0,0,0);
	  cells[271]->set_vertices(vertices[28],vertices[20],vertices[19],vertices[68]);
	  cells[271]->set_neighbors(cells[144],cells[86],cells[104],cells[1467]);
	  set_offsets(cells[271],0,0,0,0);
	  cells[272]->set_vertices(vertices[120],vertices[65],vertices[112],vertices[113]);
	  cells[272]->set_neighbors(cells[729],cells[949],cells[19],cells[456]);
	  set_offsets(cells[272],0,0,0,0);
	  cells[273]->set_vertices(vertices[87],vertices[88],vertices[128],vertices[136]);
	  cells[273]->set_neighbors(cells[425],cells[174],cells[806],cells[745]);
	  set_offsets(cells[273],0,1,1,1);
	  cells[274]->set_vertices(vertices[13],vertices[68],vertices[61],vertices[60]);
	  cells[274]->set_neighbors(cells[479],cells[61],cells[239],cells[243]);
	  set_offsets(cells[274],0,0,0,0);
	  cells[275]->set_vertices(vertices[86],vertices[37],vertices[38],vertices[46]);
	  cells[275]->set_neighbors(cells[1474],cells[315],cells[555],cells[553]);
	  set_offsets(cells[275],0,0,0,0);
	  cells[276]->set_vertices(vertices[52],vertices[5],vertices[4],vertices[45]);
	  cells[276]->set_neighbors(cells[236],cells[476],cells[108],cells[337]);
	  set_offsets(cells[276],2,2,2,0);
	  cells[277]->set_vertices(vertices[225],vertices[217],vertices[224],vertices[272]);
	  cells[277]->set_neighbors(cells[1315],cells[1598],cells[1646],cells[1174]);
	  set_offsets(cells[277],0,0,0,0);
	  cells[278]->set_vertices(vertices[5],vertices[14],vertices[62],vertices[13]);
	  cells[278]->set_neighbors(cells[79],cells[192],cells[1529],cells[321]);
	  set_offsets(cells[278],0,0,0,0);
	  cells[279]->set_vertices(vertices[29],vertices[78],vertices[77],vertices[86]);
	  cells[279]->set_neighbors(cells[442],cells[518],cells[359],cells[71]);
	  set_offsets(cells[279],0,0,0,0);
	  cells[280]->set_vertices(vertices[195],vertices[242],vertices[250],vertices[243]);
	  cells[280]->set_neighbors(cells[1009],cells[1242],cells[1380],cells[1513]);
	  set_offsets(cells[280],0,0,0,0);
	  cells[281]->set_vertices(vertices[53],vertices[110],vertices[101],vertices[61]);
	  cells[281]->set_neighbors(cells[639],cells[51],cells[676],cells[427]);
	  set_offsets(cells[281],0,0,0,0);
	  cells[282]->set_vertices(vertices[114],vertices[106],vertices[107],vertices[59]);
	  cells[282]->set_neighbors(cells[649],cells[469],cells[429],cells[501]);
	  set_offsets(cells[282],0,0,0,0);
	  cells[283]->set_vertices(vertices[27],vertices[76],vertices[19],vertices[67]);
	  cells[283]->set_neighbors(cells[480],cells[432],cells[205],cells[339]);
	  set_offsets(cells[283],0,0,0,0);
	  cells[284]->set_vertices(vertices[61],vertices[118],vertices[110],vertices[109]);
	  cells[284]->set_neighbors(cells[739],cells[639],cells[547],cells[204]);
	  set_offsets(cells[284],0,0,0,0);
	  cells[285]->set_vertices(vertices[26],vertices[27],vertices[19],vertices[74]);
	  cells[285]->set_neighbors(cells[432],cells[484],cells[186],cells[1515]);
	  set_offsets(cells[285],0,0,0,0);
	  cells[286]->set_vertices(vertices[67],vertices[115],vertices[116],vertices[107]);
	  cells[286]->set_neighbors(cells[802],cells[154],cells[742],cells[663]);
	  set_offsets(cells[286],0,0,0,0);
	  cells[287]->set_vertices(vertices[72],vertices[24],vertices[17],vertices[25]);
	  cells[287]->set_neighbors(cells[1611],cells[466],cells[92],cells[33]);
	  set_offsets(cells[287],0,0,0,0);
	  cells[288]->set_vertices(vertices[13],vertices[22],vertices[62],vertices[70]);
	  cells[288]->set_neighbors(cells[5],cells[383],cells[178],cells[79]);
	  set_offsets(cells[288],0,0,0,0);
	  cells[289]->set_vertices(vertices[1],vertices[89],vertices[50],vertices[41]);
	  cells[289]->set_neighbors(cells[399],cells[328],cells[123],cells[583]);
	  set_offsets(cells[289],2,0,2,0);
	  cells[290]->set_vertices(vertices[33],vertices[82],vertices[81],vertices[90]);
	  cells[290]->set_neighbors(cells[794],cells[560],cells[303],cells[60]);
	  set_offsets(cells[290],0,0,0,0);
	  cells[291]->set_vertices(vertices[282],vertices[42],vertices[281],vertices[242]);
	  cells[291]->set_neighbors(cells[1232],cells[1679],cells[1697],cells[1689]);
	  set_offsets(cells[291],0,4,0,2);
	  cells[292]->set_vertices(vertices[69],vertices[78],vertices[70],vertices[118]);
	  cells[292]->set_neighbors(cells[723],cells[722],cells[640],cells[46]);
	  set_offsets(cells[292],0,0,0,0);
	  cells[293]->set_vertices(vertices[82],vertices[35],vertices[83],vertices[90]);
	  cells[293]->set_neighbors(cells[329],cells[388],cells[118],cells[566]);
	  set_offsets(cells[293],0,0,0,0);
	  cells[294]->set_vertices(vertices[118],vertices[71],vertices[70],vertices[63]);
	  cells[294]->set_neighbors(cells[148],cells[629],cells[430],cells[723]);
	  set_offsets(cells[294],0,0,0,0);
	  cells[295]->set_vertices(vertices[1],vertices[8],vertices[48],vertices[0]);
	  cells[295]->set_neighbors(cells[10],cells[1],cells[1101],cells[335]);
	  set_offsets(cells[295],0,0,0,0);
	  cells[296]->set_vertices(vertices[13],vertices[62],vertices[53],vertices[61]);
	  cells[296]->set_neighbors(cells[676],cells[61],cells[383],cells[192]);
	  set_offsets(cells[296],0,0,0,0);
	  cells[297]->set_vertices(vertices[23],vertices[16],vertices[64],vertices[24]);
	  cells[297]->set_neighbors(cells[87],cells[195],cells[939],cells[453]);
	  set_offsets(cells[297],0,1,1,1);
	  cells[298]->set_vertices(vertices[68],vertices[60],vertices[11],vertices[59]);
	  cells[298]->set_neighbors(cells[3],cells[433],cells[498],cells[11]);
	  set_offsets(cells[298],0,0,0,0);
	  cells[299]->set_vertices(vertices[136],vertices[89],vertices[81],vertices[129]);
	  cells[299]->set_neighbors(cells[381],cells[357],cells[815],cells[765]);
	  set_offsets(cells[299],0,0,0,0);
	  cells[300]->set_vertices(vertices[76],vertices[29],vertices[21],vertices[69]);
	  cells[300]->set_neighbors(cells[452],cells[22],cells[528],cells[439]);
	  set_offsets(cells[300],0,0,0,0);
	  cells[301]->set_vertices(vertices[240],vertices[193],vertices[241],vertices[248]);
	  cells[301]->set_neighbors(cells[34],cells[229],cells[1432],cells[1683]);
	  set_offsets(cells[301],0,0,0,0);
	  cells[302]->set_vertices(vertices[35],vertices[84],vertices[83],vertices[92]);
	  cells[302]->set_neighbors(cells[800],cells[493],cells[143],cells[537]);
	  set_offsets(cells[302],0,0,0,0);
	  cells[303]->set_vertices(vertices[90],vertices[82],vertices[42],vertices[33]);
	  cells[303]->set_neighbors(cells[449],cells[102],cells[290],cells[118]);
	  set_offsets(cells[303],0,0,0,0);
	  cells[304]->set_vertices(vertices[57],vertices[106],vertices[97],vertices[105]);
	  cells[304]->set_neighbors(cells[366],cells[683],cells[212],cells[352]);
	  set_offsets(cells[304],0,0,0,0);
	  cells[305]->set_vertices(vertices[72],vertices[80],vertices[25],vertices[73]);
	  cells[305]->set_neighbors(cells[467],cells[505],cells[228],cells[437]);
	  set_offsets(cells[305],0,0,0,0);
	  cells[306]->set_vertices(vertices[10],vertices[58],vertices[3],vertices[50]);
	  cells[306]->set_neighbors(cells[249],cells[72],cells[190],cells[185]);
	  set_offsets(cells[306],0,0,0,0);
	  cells[307]->set_vertices(vertices[227],vertices[234],vertices[179],vertices[187]);
	  cells[307]->set_neighbors(cells[1171],cells[1200],cells[1418],cells[1368]);
	  set_offsets(cells[307],0,0,0,0);
	  cells[308]->set_vertices(vertices[164],vertices[156],vertices[109],vertices[116]);
	  cells[308]->set_neighbors(cells[751],cells[642],cells[761],cells[647]);
	  set_offsets(cells[308],0,0,0,0);
	  cells[309]->set_vertices(vertices[147],vertices[196],vertices[195],vertices[204]);
	  cells[309]->set_neighbors(cells[952],cells[1012],cells[976],cells[1376]);
	  set_offsets(cells[309],0,0,0,0);
	  cells[310]->set_vertices(vertices[35],vertices[27],vertices[34],vertices[82]);
	  cells[310]->set_neighbors(cells[482],cells[497],cells[514],cells[1660]);
	  set_offsets(cells[310],0,0,0,0);
	  cells[311]->set_vertices(vertices[77],vertices[85],vertices[86],vertices[134]);
	  cells[311]->set_neighbors(cells[142],cells[777],cells[462],cells[496]);
	  set_offsets(cells[311],0,0,0,0);
	  cells[312]->set_vertices(vertices[98],vertices[58],vertices[51],vertices[106]);
	  cells[312]->set_neighbors(cells[681],cells[511],cells[260],cells[448]);
	  set_offsets(cells[312],0,0,0,0);
	  cells[313]->set_vertices(vertices[90],vertices[35],vertices[43],vertices[42]);
	  cells[313]->set_neighbors(cells[1139],cells[477],cells[118],cells[329]);
	  set_offsets(cells[313],0,0,0,0);
	  cells[314]->set_vertices(vertices[15],vertices[63],vertices[56],vertices[55]);
	  cells[314]->set_neighbors(cells[121],cells[13],cells[21],cells[141]);
	  set_offsets(cells[314],0,0,1,0);
	  cells[315]->set_vertices(vertices[39],vertices[86],vertices[38],vertices[46]);
	  cells[315]->set_neighbors(cells[275],cells[1671],cells[500],cells[552]);
	  set_offsets(cells[315],0,0,0,0);
	  cells[316]->set_vertices(vertices[184],vertices[96],vertices[143],vertices[136]);
	  cells[316]->set_neighbors(cells[671],cells[816],cells[685],cells[1120]);
	  set_offsets(cells[316],1,3,0,1);
	  cells[317]->set_vertices(vertices[69],vertices[116],vertices[109],vertices[61]);
	  cells[317]->set_neighbors(cells[700],cells[547],cells[45],cells[752]);
	  set_offsets(cells[317],0,0,0,0);
	  cells[318]->set_vertices(vertices[66],vertices[11],vertices[59],vertices[19]);
	  cells[318]->set_neighbors(cells[433],cells[470],cells[423],cells[451]);
	  set_offsets(cells[318],0,0,0,0);
	  cells[319]->set_vertices(vertices[20],vertices[12],vertices[60],vertices[13]);
	  cells[319]->set_neighbors(cells[206],cells[239],cells[1279],cells[179]);
	  set_offsets(cells[319],0,0,0,0);
	  cells[320]->set_vertices(vertices[50],vertices[58],vertices[98],vertices[49]);
	  cells[320]->set_neighbors(cells[260],cells[631],cells[49],cells[448]);
	  set_offsets(cells[320],0,0,0,0);
	  cells[321]->set_vertices(vertices[5],vertices[14],vertices[54],vertices[62]);
	  cells[321]->set_neighbors(cells[171],cells[207],cells[278],cells[62]);
	  set_offsets(cells[321],0,0,0,0);
	  cells[322]->set_vertices(vertices[277],vertices[269],vertices[38],vertices[278]);
	  cells[322]->set_neighbors(cells[346],cells[1719],cells[1472],cells[1181]);
	  set_offsets(cells[322],0,0,4,0);
	  cells[323]->set_vertices(vertices[191],vertices[238],vertices[231],vertices[183]);
	  cells[323]->set_neighbors(cells[972],cells[1401],cells[325],cells[1434]);
	  set_offsets(cells[323],0,0,0,0);
	  cells[324]->set_vertices(vertices[284],vertices[244],vertices[44],vertices[283]);
	  cells[324]->set_neighbors(cells[237],cells[1623],cells[1701],cells[1709]);
	  set_offsets(cells[324],0,2,4,0);
	  cells[325]->set_vertices(vertices[183],vertices[191],vertices[238],vertices[190]);
	  cells[325]->set_neighbors(cells[973],cells[419],cells[1155],cells[323]);
	  set_offsets(cells[325],0,0,0,0);
	  cells[326]->set_vertices(vertices[227],vertices[219],vertices[276],vertices[228]);
	  cells[326]->set_neighbors(cells[1374],cells[1362],cells[1201],cells[1625]);
	  set_offsets(cells[326],0,0,0,0);
	  cells[327]->set_vertices(vertices[228],vertices[236],vertices[188],vertices[181]);
	  cells[327]->set_neighbors(cells[1114],cells[1325],cells[1184],cells[133]);
	  set_offsets(cells[327],0,0,0,0);
	  cells[328]->set_vertices(vertices[1],vertices[50],vertices[2],vertices[41]);
	  cells[328]->set_neighbors(cells[564],cells[1317],cells[289],cells[53]);
	  set_offsets(cells[328],2,2,2,0);
	  cells[329]->set_vertices(vertices[90],vertices[35],vertices[83],vertices[43]);
	  cells[329]->set_neighbors(cells[493],cells[152],cells[313],cells[293]);
	  set_offsets(cells[329],0,0,0,0);
	  cells[330]->set_vertices(vertices[103],vertices[55],vertices[104],vertices[63]);
	  cells[330]->set_neighbors(cells[121],cells[698],cells[699],cells[132]);
	  set_offsets(cells[330],0,0,1,0);
	  cells[331]->set_vertices(vertices[3],vertices[52],vertices[4],vertices[43]);
	  cells[331]->set_neighbors(cells[525],cells[890],cells[563],cells[31]);
	  set_offsets(cells[331],2,2,2,0);
	  cells[332]->set_vertices(vertices[141],vertices[148],vertices[188],vertices[100]);
	  cells[332]->set_neighbors(cells[1135],cells[1113],cells[707],cells[116]);
	  set_offsets(cells[332],0,2,0,2);
	  cells[333]->set_vertices(vertices[146],vertices[106],vertices[154],vertices[97]);
	  cells[333]->set_neighbors(cells[366],cells[457],cells[626],cells[754]);
	  set_offsets(cells[333],0,0,0,0);
	  cells[334]->set_vertices(vertices[70],vertices[30],vertices[23],vertices[78]);
	  cells[334]->set_neighbors(cells[405],cells[355],cells[380],cells[396]);
	  set_offsets(cells[334],0,0,0,0);
	  cells[335]->set_vertices(vertices[48],vertices[8],vertices[1],vertices[56]);
	  cells[335]->set_neighbors(cells[336],cells[97],cells[188],cells[295]);
	  set_offsets(cells[335],0,0,0,0);
	  cells[336]->set_vertices(vertices[56],vertices[8],vertices[1],vertices[9]);
	  cells[336]->set_neighbors(cells[1506],cells[131],cells[101],cells[335]);
	  set_offsets(cells[336],0,0,0,0);
	  cells[337]->set_vertices(vertices[52],vertices[4],vertices[5],vertices[12]);
	  cells[337]->set_neighbors(cells[966],cells[42],cells[31],cells[276]);
	  set_offsets(cells[337],0,0,0,0);
	  cells[338]->set_vertices(vertices[65],vertices[17],vertices[57],vertices[66]);
	  cells[338]->set_neighbors(cells[36],cells[257],cells[486],cells[454]);
	  set_offsets(cells[338],0,0,0,0);
	  cells[339]->set_vertices(vertices[27],vertices[28],vertices[19],vertices[76]);
	  cells[339]->set_neighbors(cells[86],cells[283],cells[99],cells[1195]);
	  set_offsets(cells[339],0,0,0,0);
	  cells[340]->set_vertices(vertices[138],vertices[91],vertices[83],vertices[131]);
	  cells[340]->set_neighbors(cells[835],cells[829],cells[376],cells[820]);
	  set_offsets(cells[340],0,0,0,0);
	  cells[341]->set_vertices(vertices[62],vertices[70],vertices[15],vertices[63]);
	  cells[341]->set_neighbors(cells[122],cells[21],cells[373],cells[5]);
	  set_offsets(cells[341],0,0,0,0);
	  cells[342]->set_vertices(vertices[250],vertices[241],vertices[201],vertices[249]);
	  cells[342]->set_neighbors(cells[1509],cells[1562],cells[1422],cells[1502]);
	  set_offsets(cells[342],0,0,0,0);
	  cells[343]->set_vertices(vertices[244],vertices[197],vertices[285],vertices[245]);
	  cells[343]->set_neighbors(cells[1669],cells[1630],cells[1522],cells[1712]);
	  set_offsets(cells[343],2,2,0,2);
	  cells[344]->set_vertices(vertices[282],vertices[233],vertices[225],vertices[234]);
	  cells[344]->set_neighbors(cells[1415],cells[1100],cells[1617],cells[1693]);
	  set_offsets(cells[344],0,0,0,0);
	  cells[345]->set_vertices(vertices[239],vertices[151],vertices[191],vertices[192]);
	  cells[345]->set_neighbors(cells[909],cells[613],cells[1355],cells[1436]);
	  set_offsets(cells[345],0,2,0,3);
	  cells[346]->set_vertices(vertices[38],vertices[269],vertices[30],vertices[278]);
	  cells[346]->set_neighbors(cells[911],cells[1581],cells[322],cells[1583]);
	  set_offsets(cells[346],4,0,4,0);
	  cells[347]->set_vertices(vertices[237],vertices[189],vertices[229],vertices[238]);
	  cells[347]->set_neighbors(cells[1163],cells[1627],cells[1363],cells[1427]);
	  set_offsets(cells[347],0,0,0,0);
	  cells[348]->set_vertices(vertices[230],vertices[190],vertices[238],vertices[181]);
	  cells[348]->set_neighbors(cells[921],cells[1168],cells[1165],cells[419]);
	  set_offsets(cells[348],0,0,0,0);
	  cells[349]->set_vertices(vertices[140],vertices[91],vertices[92],vertices[52]);
	  cells[349]->set_neighbors(cells[167],cells[50],cells[492],cells[416]);
	  set_offsets(cells[349],0,0,0,2);
	  cells[350]->set_vertices(vertices[96],vertices[136],vertices[89],vertices[48]);
	  cells[350]->set_neighbors(cells[787],cells[111],cells[839],cells[857]);
	  set_offsets(cells[350],2,0,0,2);
	  cells[351]->set_vertices(vertices[142],vertices[54],vertices[94],vertices[95]);
	  cells[351]->set_neighbors(cells[605],cells[849],cells[675],cells[881]);
	  set_offsets(cells[351],0,2,0,0);
	  cells[352]->set_vertices(vertices[106],vertices[57],vertices[97],vertices[49]);
	  cells[352]->set_neighbors(cells[588],cells[483],cells[431],cells[304]);
	  set_offsets(cells[352],0,0,0,0);
	  cells[353]->set_vertices(vertices[65],vertices[74],vertices[114],vertices[122]);
	  cells[353]->set_neighbors(cells[397],cells[231],cells[616],cells[240]);
	  set_offsets(cells[353],0,0,0,0);
	  cells[354]->set_vertices(vertices[118],vertices[158],vertices[111],vertices[110]);
	  cells[354]->set_neighbors(cells[702],cells[225],cells[739],cells[690]);
	  set_offsets(cells[354],0,0,0,0);
	  cells[355]->set_vertices(vertices[70],vertices[78],vertices[23],vertices[71]);
	  cells[355]->set_neighbors(cells[356],cells[148],cells[723],cells[334]);
	  set_offsets(cells[355],0,0,0,0);
	  cells[356]->set_vertices(vertices[78],vertices[31],vertices[23],vertices[71]);
	  cells[356]->set_neighbors(cells[428],cells[355],cells[393],cells[405]);
	  set_offsets(cells[356],0,0,0,0);
	  cells[357]->set_vertices(vertices[128],vertices[136],vertices[81],vertices[129]);
	  cells[357]->set_neighbors(cells[299],cells[617],cells[1073],cells[425]);
	  set_offsets(cells[357],0,0,0,0);
	  cells[358]->set_vertices(vertices[74],vertices[73],vertices[65],vertices[25]);
	  cells[358]->set_neighbors(cells[505],cells[90],cells[377],cells[616]);
	  set_offsets(cells[358],0,0,0,0);
	  cells[359]->set_vertices(vertices[29],vertices[38],vertices[78],vertices[86]);
	  cells[359]->set_neighbors(cells[551],cells[279],cells[553],cells[26]);
	  set_offsets(cells[359],0,0,0,0);
	  cells[360]->set_vertices(vertices[119],vertices[160],vertices[111],vertices[159]);
	  cells[360]->set_neighbors(cells[767],cells[750],cells[1036],cells[762]);
	  set_offsets(cells[360],0,1,0,0);
	  cells[361]->set_vertices(vertices[26],vertices[18],vertices[17],vertices[66]);
	  cells[361]->set_neighbors(cells[232],cells[485],cells[370],cells[559]);
	  set_offsets(cells[361],0,0,0,0);
	  cells[362]->set_vertices(vertices[118],vertices[126],vertices[78],vertices[71]);
	  cells[362]->set_neighbors(cells[544],cells[723],cells[755],cells[640]);
	  set_offsets(cells[362],0,0,0,0);
	  cells[363]->set_vertices(vertices[47],vertices[88],vertices[48],vertices[0]);
	  cells[363]->set_neighbors(cells[567],cells[27],cells[93],cells[608]);
	  set_offsets(cells[363],0,1,3,3);
	  cells[364]->set_vertices(vertices[211],vertices[210],vertices[203],vertices[163]);
	  cells[364]->set_neighbors(cells[968],cells[905],cells[1282],cells[1048]);
	  set_offsets(cells[364],0,0,0,0);
	  cells[365]->set_vertices(vertices[132],vertices[84],vertices[124],vertices[77]);
	  cells[365]->set_neighbors(cells[545],cells[472],cells[579],cells[89]);
	  set_offsets(cells[365],0,0,0,0);
	  cells[366]->set_vertices(vertices[154],vertices[106],vertices[105],vertices[97]);
	  cells[366]->set_neighbors(cells[304],cells[557],cells[333],cells[950]);
	  set_offsets(cells[366],0,0,0,0);
	  cells[367]->set_vertices(vertices[13],vertices[70],vertices[61],vertices[21]);
	  cells[367]->set_neighbors(cells[201],cells[243],cells[178],cells[383]);
	  set_offsets(cells[367],0,0,0,0);
	  cells[368]->set_vertices(vertices[44],vertices[37],vertices[36],vertices[84]);
	  cells[368]->set_neighbors(cells[119],cells[538],cells[524],cells[1193]);
	  set_offsets(cells[368],0,0,0,0);
	  cells[369]->set_vertices(vertices[132],vertices[140],vertices[83],vertices[92]);
	  cells[369]->set_neighbors(cells[416],cells[800],cells[871],cells[743]);
	  set_offsets(cells[369],0,0,0,0);
	  cells[370]->set_vertices(vertices[26],vertices[18],vertices[66],vertices[19]);
	  cells[370]->set_neighbors(cells[423],cells[484],cells[1569],cells[361]);
	  set_offsets(cells[370],0,0,0,0);
	  cells[371]->set_vertices(vertices[102],vertices[54],vertices[55],vertices[62]);
	  cells[371]->set_neighbors(cells[114],cells[23],cells[411],cells[669]);
	  set_offsets(cells[371],0,0,0,0);
	  cells[372]->set_vertices(vertices[52],vertices[60],vertices[100],vertices[51]);
	  cells[372]->set_neighbors(cells[658],cells[162],cells[223],cells[650]);
	  set_offsets(cells[372],0,0,0,0);
	  cells[373]->set_vertices(vertices[110],vertices[70],vertices[62],vertices[63]);
	  cells[373]->set_neighbors(cells[341],cells[487],cells[629],cells[677]);
	  set_offsets(cells[373],0,0,0,0);
	  cells[374]->set_vertices(vertices[50],vertices[3],vertices[2],vertices[43]);
	  cells[374]->set_neighbors(cells[1228],cells[536],cells[519],cells[72]);
	  set_offsets(cells[374],2,2,2,0);
	  cells[375]->set_vertices(vertices[116],vertices[76],vertices[69],vertices[124]);
	  cells[375]->set_neighbors(cells[666],cells[410],cells[619],cells[258]);
	  set_offsets(cells[375],0,0,0,0);
	  cells[376]->set_vertices(vertices[139],vertices[91],vertices[138],vertices[131]);
	  cells[376]->set_neighbors(cells[340],cells[1126],cells[836],cells[147]);
	  set_offsets(cells[376],0,0,0,0);
	  cells[377]->set_vertices(vertices[82],vertices[73],vertices[74],vertices[25]);
	  cells[377]->set_neighbors(cells[358],cells[520],cells[465],cells[491]);
	  set_offsets(cells[377],0,0,0,0);
	  cells[378]->set_vertices(vertices[89],vertices[88],vertices[81],vertices[41]);
	  cells[378]->set_neighbors(cells[403],cells[198],cells[70],cells[765]);
	  set_offsets(cells[378],0,0,0,0);
	  cells[379]->set_vertices(vertices[75],vertices[122],vertices[115],vertices[67]);
	  cells[379]->set_neighbors(cells[20],cells[494],cells[644],cells[474]);
	  set_offsets(cells[379],0,0,0,0);
	  cells[380]->set_vertices(vertices[21],vertices[30],vertices[70],vertices[78]);
	  cells[380]->set_neighbors(cells[334],cells[46],cells[404],cells[406]);
	  set_offsets(cells[380],0,0,0,0);
	  cells[381]->set_vertices(vertices[138],vertices[89],vertices[129],vertices[81]);
	  cells[381]->set_neighbors(cells[299],cells[78],cells[789],cells[161]);
	  set_offsets(cells[381],0,0,0,0);
	  cells[382]->set_vertices(vertices[100],vertices[53],vertices[101],vertices[108]);
	  cells[382]->set_neighbors(cells[51],cells[734],cells[438],cells[887]);
	  set_offsets(cells[382],0,0,0,0);
	  cells[383]->set_vertices(vertices[13],vertices[62],vertices[61],vertices[70]);
	  cells[383]->set_neighbors(cells[677],cells[367],cells[288],cells[296]);
	  set_offsets(cells[383],0,0,0,0);
	  cells[384]->set_vertices(vertices[94],vertices[39],vertices[47],vertices[46]);
	  cells[384]->set_neighbors(cells[1021],cells[582],cells[500],cells[530]);
	  set_offsets(cells[384],0,0,0,0);
	  cells[385]->set_vertices(vertices[15],vertices[8],vertices[56],vertices[16]);
	  cells[385]->set_neighbors(cells[101],cells[184],cells[587],cells[112]);
	  set_offsets(cells[385],0,1,1,1);
	  cells[386]->set_vertices(vertices[59],vertices[19],vertices[68],vertices[67]);
	  cells[386]->set_neighbors(cells[480],cells[549],cells[470],cells[433]);
	  set_offsets(cells[386],0,0,0,0);
	  cells[387]->set_vertices(vertices[86],vertices[39],vertices[31],vertices[79]);
	  cells[387]->set_neighbors(cells[251],cells[435],cells[578],cells[552]);
	  set_offsets(cells[387],0,0,0,0);
	  cells[388]->set_vertices(vertices[90],vertices[82],vertices[130],vertices[83]);
	  cells[388]->set_neighbors(cells[68],cells[776],cells[293],cells[794]);
	  set_offsets(cells[388],0,0,0,0);
	  cells[389]->set_vertices(vertices[142],vertices[85],vertices[133],vertices[134]);
	  cells[389]->set_neighbors(cells[109],cells[1105],cells[880],cells[834]);
	  set_offsets(cells[389],0,0,0,0);
	  cells[390]->set_vertices(vertices[148],vertices[101],vertices[141],vertices[189]);
	  cells[390]->set_neighbors(cells[1104],cells[116],cells[1146],cells[707]);
	  set_offsets(cells[390],2,2,0,0);
	  cells[391]->set_vertices(vertices[211],vertices[258],vertices[251],vertices[203]);
	  cells[391]->set_neighbors(cells[1549],cells[1252],cells[1048],cells[1520]);
	  set_offsets(cells[391],0,0,0,0);
	  cells[392]->set_vertices(vertices[31],vertices[30],vertices[23],vertices[271]);
	  cells[392]->set_neighbors(cells[1593],cells[1261],cells[1395],cells[405]);
	  set_offsets(cells[392],4,4,4,0);
	  cells[393]->set_vertices(vertices[78],vertices[79],vertices[31],vertices[71]);
	  cells[393]->set_neighbors(cells[401],cells[356],cells[544],cells[435]);
	  set_offsets(cells[393],0,0,0,0);
	  cells[394]->set_vertices(vertices[81],vertices[73],vertices[128],vertices[121]);
	  cells[394]->set_neighbors(cells[659],cells[617],cells[764],cells[774]);
	  set_offsets(cells[394],0,0,0,0);
	  cells[395]->set_vertices(vertices[104],vertices[49],vertices[57],vertices[56]);
	  cells[395]->set_neighbors(cells[412],cells[215],cells[637],cells[588]);
	  set_offsets(cells[395],0,0,0,0);
	  cells[396]->set_vertices(vertices[70],vertices[22],vertices[23],vertices[30]);
	  cells[396]->set_neighbors(cells[974],cells[334],cells[406],cells[447]);
	  set_offsets(cells[396],0,0,0,0);
	  cells[397]->set_vertices(vertices[114],vertices[74],vertices[67],vertices[122]);
	  cells[397]->set_neighbors(cells[644],cells[20],cells[353],cells[63]);
	  set_offsets(cells[397],0,0,0,0);
	  cells[398]->set_vertices(vertices[101],vertices[110],vertices[158],vertices[109]);
	  cells[398]->set_neighbors(cells[739],cells[691],cells[639],cells[740]);
	  set_offsets(cells[398],0,0,0,0);
	  cells[399]->set_vertices(vertices[90],vertices[50],vertices[89],vertices[41]);
	  cells[399]->set_neighbors(cells[289],cells[198],cells[564],cells[821]);
	  set_offsets(cells[399],0,2,0,0);
	  cells[400]->set_vertices(vertices[77],vertices[125],vertices[126],vertices[117]);
	  cells[400]->set_neighbors(cells[933],cells[655],cells[796],cells[799]);
	  set_offsets(cells[400],0,0,0,0);
	  cells[401]->set_vertices(vertices[31],vertices[79],vertices[72],vertices[71]);
	  cells[401]->set_neighbors(cells[732],cells[428],cells[393],cells[38]);
	  set_offsets(cells[401],0,0,1,0);
	  cells[402]->set_vertices(vertices[80],vertices[32],vertices[25],vertices[33]);
	  cells[402]->set_neighbors(cells[1412],cells[467],cells[160],cells[437]);
	  set_offsets(cells[402],0,0,0,0);
	  cells[403]->set_vertices(vertices[88],vertices[33],vertices[81],vertices[41]);
	  cells[403]->set_neighbors(cells[560],cells[378],cells[507],cells[534]);
	  set_offsets(cells[403],0,0,0,0);
	  cells[404]->set_vertices(vertices[21],vertices[30],vertices[78],vertices[29]);
	  cells[404]->set_neighbors(cells[26],cells[452],cells[1584],cells[380]);
	  set_offsets(cells[404],0,0,0,0);
	  cells[405]->set_vertices(vertices[78],vertices[30],vertices[23],vertices[31]);
	  cells[405]->set_neighbors(cells[392],cells[356],cells[464],cells[334]);
	  set_offsets(cells[405],0,0,0,0);
	  cells[406]->set_vertices(vertices[21],vertices[22],vertices[70],vertices[30]);
	  cells[406]->set_neighbors(cells[396],cells[380],cells[1535],cells[178]);
	  set_offsets(cells[406],0,0,0,0);
	  cells[407]->set_vertices(vertices[37],vertices[86],vertices[85],vertices[94]);
	  cells[407]->set_neighbors(cells[142],cells[573],cells[555],cells[496]);
	  set_offsets(cells[407],0,0,0,0);
	  cells[408]->set_vertices(vertices[178],vertices[121],vertices[170],vertices[130]);
	  cells[408]->set_neighbors(cells[521],cells[1051],cells[1065],cells[634]);
	  set_offsets(cells[408],0,0,0,0);
	  cells[409]->set_vertices(vertices[158],vertices[103],vertices[151],vertices[111]);
	  cells[409]->set_neighbors(cells[898],cells[585],cells[702],cells[670]);
	  set_offsets(cells[409],0,0,0,0);
	  cells[410]->set_vertices(vertices[124],vertices[116],vertices[117],vertices[69]);
	  cells[410]->set_neighbors(cells[752],cells[721],cells[375],cells[693]);
	  set_offsets(cells[410],0,0,0,0);
	  cells[411]->set_vertices(vertices[53],vertices[54],vertices[102],vertices[62]);
	  cells[411]->set_neighbors(cells[371],cells[627],cells[207],cells[668]);
	  set_offsets(cells[411],0,0,0,0);
	  cells[412]->set_vertices(vertices[56],vertices[57],vertices[9],vertices[49]);
	  cells[412]->set_neighbors(cells[177],cells[131],cells[395],cells[420]);
	  set_offsets(cells[412],0,0,0,0);
	  cells[413]->set_vertices(vertices[139],vertices[186],vertices[138],vertices[98]);
	  cells[413]->set_neighbors(cells[688],cells[147],cells[1081],cells[1126]);
	  set_offsets(cells[413],0,0,0,2);
	  cells[414]->set_vertices(vertices[134],vertices[86],vertices[126],vertices[79]);
	  cells[414]->set_neighbors(cells[440],cells[807],cells[55],cells[777]);
	  set_offsets(cells[414],0,0,0,0);
	  cells[415]->set_vertices(vertices[113],vertices[105],vertices[162],vertices[114]);
	  cells[415]->set_neighbors(cells[793],cells[748],cells[504],cells[512]);
	  set_offsets(cells[415],0,0,0,0);
	  cells[416]->set_vertices(vertices[140],vertices[91],vertices[83],vertices[92]);
	  cells[416]->set_neighbors(cells[211],cells[369],cells[349],cells[835]);
	  set_offsets(cells[416],0,0,0,0);
	  cells[417]->set_vertices(vertices[148],vertices[101],vertices[149],vertices[156]);
	  cells[417]->set_neighbors(cells[747],cells[852],cells[471],cells[1146]);
	  set_offsets(cells[417],0,0,0,0);
	  cells[418]->set_vertices(vertices[75],vertices[124],vertices[123],vertices[132]);
	  cells[418]->set_neighbors(cells[572],cells[801],cells[89],cells[548]);
	  set_offsets(cells[418],0,0,0,0);
	  cells[419]->set_vertices(vertices[183],vertices[238],vertices[230],vertices[190]);
	  cells[419]->set_neighbors(cells[348],cells[1164],cells[325],cells[972]);
	  set_offsets(cells[419],0,0,0,0);
	  cells[420]->set_vertices(vertices[56],vertices[64],vertices[9],vertices[57]);
	  cells[420]->set_neighbors(cells[270],cells[412],cells[215],cells[199]);
	  set_offsets(cells[420],0,0,0,0);
	  cells[421]->set_vertices(vertices[259],vertices[258],vertices[18],vertices[251]);
	  cells[421]->set_neighbors(cells[749],cells[1518],cells[1520],cells[214]);
	  set_offsets(cells[421],0,0,4,0);
	  cells[422]->set_vertices(vertices[182],vertices[174],vertices[175],vertices[127]);
	  cells[422]->set_neighbors(cells[1077],cells[1076],cells[780],cells[847]);
	  set_offsets(cells[422],0,0,0,0);
	  cells[423]->set_vertices(vertices[19],vertices[66],vertices[11],vertices[18]);
	  cells[423]->set_neighbors(cells[202],cells[1568],cells[370],cells[318]);
	  set_offsets(cells[423],0,0,0,0);
	  cells[424]->set_vertices(vertices[58],vertices[11],vertices[51],vertices[59]);
	  cells[424]->set_neighbors(cells[3],cells[681],cells[451],cells[15]);
	  set_offsets(cells[424],0,0,0,0);
	  cells[425]->set_vertices(vertices[128],vertices[88],vertices[81],vertices[136]);
	  cells[425]->set_neighbors(cells[765],cells[357],cells[273],cells[773]);
	  set_offsets(cells[425],0,0,0,0);
	  cells[426]->set_vertices(vertices[54],vertices[47],vertices[6],vertices[94]);
	  cells[426]->set_neighbors(cells[582],cells[175],cells[605],cells[103]);
	  set_offsets(cells[426],2,0,2,0);
	  cells[427]->set_vertices(vertices[53],vertices[102],vertices[101],vertices[110]);
	  cells[427]->set_neighbors(cells[831],cells[281],cells[627],cells[888]);
	  set_offsets(cells[427],0,0,0,0);
	  cells[428]->set_vertices(vertices[23],vertices[31],vertices[72],vertices[71]);
	  cells[428]->set_neighbors(cells[401],cells[203],cells[356],cells[2]);
	  set_offsets(cells[428],0,0,1,0);
	  cells[429]->set_vertices(vertices[106],vertices[66],vertices[59],vertices[114]);
	  cells[429]->set_neighbors(cells[82],cells[282],cells[679],cells[96]);
	  set_offsets(cells[429],0,0,0,0);
	  cells[430]->set_vertices(vertices[71],vertices[118],vertices[111],vertices[63]);
	  cells[430]->set_neighbors(cells[225],cells[682],cells[294],cells[716]);
	  set_offsets(cells[430],0,0,0,0);
	  cells[431]->set_vertices(vertices[58],vertices[57],vertices[106],vertices[49]);
	  cells[431]->set_neighbors(cells[352],cells[260],cells[177],cells[17]);
	  set_offsets(cells[431],0,0,0,0);
	  cells[432]->set_vertices(vertices[74],vertices[19],vertices[67],vertices[27]);
	  cells[432]->set_neighbors(cells[283],cells[515],cells[285],cells[32]);
	  set_offsets(cells[432],0,0,0,0);
	  cells[433]->set_vertices(vertices[19],vertices[68],vertices[11],vertices[59]);
	  cells[433]->set_neighbors(cells[298],cells[318],cells[386],cells[144]);
	  set_offsets(cells[433],0,0,0,0);
	  cells[434]->set_vertices(vertices[13],vertices[68],vertices[20],vertices[21]);
	  cells[434]->set_neighbors(cells[104],cells[1573],cells[243],cells[239]);
	  set_offsets(cells[434],0,0,0,0);
	  cells[435]->set_vertices(vertices[78],vertices[86],vertices[31],vertices[79]);
	  cells[435]->set_neighbors(cells[387],cells[393],cells[440],cells[551]);
	  set_offsets(cells[435],0,0,0,0);
	  cells[436]->set_vertices(vertices[112],vertices[152],vertices[111],vertices[160]);
	  cells[436]->set_neighbors(cells[767],cells[762],cells[678],cells[680]);
	  set_offsets(cells[436],1,1,0,1);
	  cells[437]->set_vertices(vertices[72],vertices[32],vertices[25],vertices[80]);
	  cells[437]->set_neighbors(cells[402],cells[305],cells[508],cells[92]);
	  set_offsets(cells[437],0,0,0,0);
	  cells[438]->set_vertices(vertices[60],vertices[108],vertices[53],vertices[100]);
	  cells[438]->set_neighbors(cells[382],cells[650],cells[658],cells[247]);
	  set_offsets(cells[438],0,0,0,0);
	  cells[439]->set_vertices(vertices[28],vertices[29],vertices[21],vertices[76]);
	  cells[439]->set_neighbors(cells[300],cells[481],cells[107],cells[1464]);
	  set_offsets(cells[439],0,0,0,0);
	  cells[440]->set_vertices(vertices[126],vertices[86],vertices[78],vertices[79]);
	  cells[440]->set_neighbors(cells[435],cells[544],cells[414],cells[442]);
	  set_offsets(cells[440],0,0,0,0);
	  cells[441]->set_vertices(vertices[163],vertices[155],vertices[203],vertices[212]);
	  cells[441]->set_neighbors(cells[1281],cells[905],cells[57],cells[968]);
	  set_offsets(cells[441],0,0,0,0);
	  cells[442]->set_vertices(vertices[77],vertices[86],vertices[78],vertices[126]);
	  cells[442]->set_neighbors(cells[440],cells[561],cells[777],cells[279]);
	  set_offsets(cells[442],0,0,0,0);
	  cells[443]->set_vertices(vertices[197],vertices[204],vertices[157],vertices[205]);
	  cells[443]->set_neighbors(cells[265],cells[460],cells[1470],cells[35]);
	  set_offsets(cells[443],0,0,0,0);
	  cells[444]->set_vertices(vertices[169],vertices[121],vertices[161],vertices[170]);
	  cells[444]->set_neighbors(cells[139],cells[165],cells[634],cells[1044]);
	  set_offsets(cells[444],0,0,0,0);
	  cells[445]->set_vertices(vertices[129],vertices[137],vertices[184],vertices[136]);
	  cells[445]->set_neighbors(cells[685],cells[965],cells[815],cells[853]);
	  set_offsets(cells[445],0,0,0,0);
	  cells[446]->set_vertices(vertices[79],vertices[126],vertices[119],vertices[71]);
	  cells[446]->set_neighbors(cells[755],cells[728],cells[544],cells[759]);
	  set_offsets(cells[446],0,0,0,0);
	  cells[447]->set_vertices(vertices[70],vertices[22],vertices[15],vertices[23]);
	  cells[447]->set_neighbors(cells[1525],cells[122],cells[396],cells[5]);
	  set_offsets(cells[447],0,0,0,0);
	  cells[448]->set_vertices(vertices[50],vertices[58],vertices[51],vertices[98]);
	  cells[448]->set_neighbors(cells[312],cells[187],cells[320],cells[249]);
	  set_offsets(cells[448],0,0,0,0);
	  cells[449]->set_vertices(vertices[42],vertices[82],vertices[34],vertices[33]);
	  cells[449]->set_neighbors(cells[74],cells[1438],cells[303],cells[497]);
	  set_offsets(cells[449],0,0,0,0);
	  cells[450]->set_vertices(vertices[57],vertices[114],vertices[105],vertices[65]);
	  cells[450]->set_neighbors(cells[504],cells[510],cells[257],cells[212]);
	  set_offsets(cells[450],0,0,0,0);
	  cells[451]->set_vertices(vertices[58],vertices[11],vertices[59],vertices[66]);
	  cells[451]->set_neighbors(cells[318],cells[96],cells[202],cells[424]);
	  set_offsets(cells[451],0,0,0,0);
	  cells[452]->set_vertices(vertices[21],vertices[78],vertices[69],vertices[29]);
	  cells[452]->set_neighbors(cells[71],cells[300],cells[404],cells[46]);
	  set_offsets(cells[452],0,0,0,0);
	  cells[453]->set_vertices(vertices[15],vertices[16],vertices[64],vertices[23]);
	  cells[453]->set_neighbors(cells[297],cells[135],cells[1235],cells[184]);
	  set_offsets(cells[453],0,1,1,0);
	  cells[454]->set_vertices(vertices[64],vertices[65],vertices[17],vertices[57]);
	  cells[454]->set_neighbors(cells[338],cells[270],cells[569],cells[213]);
	  set_offsets(cells[454],0,0,0,0);
	  cells[455]->set_vertices(vertices[151],vertices[192],vertices[200],vertices[152]);
	  cells[455]->set_neighbors(cells[1079],cells[941],cells[893],cells[1221]);
	  set_offsets(cells[455],0,1,1,1);
	  cells[456]->set_vertices(vertices[65],vertices[72],vertices[120],vertices[112]);
	  cells[456]->set_neighbors(cells[730],cells[272],cells[180],cells[738]);
	  set_offsets(cells[456],0,0,0,0);
	  cells[457]->set_vertices(vertices[145],vertices[146],vertices[154],vertices[97]);
	  cells[457]->set_neighbors(cells[333],cells[557],cells[1091],cells[719]);
	  set_offsets(cells[457],0,0,0,0);
	  cells[458]->set_vertices(vertices[84],vertices[75],vertices[76],vertices[27]);
	  cells[458]->set_neighbors(cells[205],cells[220],cells[565],cells[80]);
	  set_offsets(cells[458],0,0,0,0);
	  cells[459]->set_vertices(vertices[125],vertices[172],vertices[117],vertices[124]);
	  cells[459]->set_neighbors(cells[984],cells[796],cells[928],cells[91]);
	  set_offsets(cells[459],0,0,0,0);
	  cells[460]->set_vertices(vertices[205],vertices[157],vertices[197],vertices[206]);
	  cells[460]->set_neighbors(cells[558],cells[1225],cells[1294],cells[443]);
	  set_offsets(cells[460],0,0,0,0);
	  cells[461]->set_vertices(vertices[147],vertices[99],vertices[154],vertices[146]);
	  cells[461]->set_neighbors(cells[754],cells[922],cells[1136],cells[704]);
	  set_offsets(cells[461],0,0,0,0);
	  cells[462]->set_vertices(vertices[77],vertices[85],vertices[134],vertices[125]);
	  cells[462]->set_neighbors(cells[109],cells[799],cells[795],cells[311]);
	  set_offsets(cells[462],0,0,0,0);
	  cells[463]->set_vertices(vertices[174],vertices[117],vertices[166],vertices[126]);
	  cells[463]->set_neighbors(cells[844],cells[785],cells[933],cells[1042]);
	  set_offsets(cells[463],0,0,0,0);
	  cells[464]->set_vertices(vertices[78],vertices[30],vertices[31],vertices[38]);
	  cells[464]->set_neighbors(cells[1395],cells[551],cells[26],cells[405]);
	  set_offsets(cells[464],0,0,0,0);
	  cells[465]->set_vertices(vertices[33],vertices[73],vertices[82],vertices[25]);
	  cells[465]->set_neighbors(cells[377],cells[74],cells[467],cells[60]);
	  set_offsets(cells[465],0,0,0,0);
	  cells[466]->set_vertices(vertices[72],vertices[25],vertices[17],vertices[65]);
	  cells[466]->set_neighbors(cells[90],cells[213],cells[505],cells[287]);
	  set_offsets(cells[466],0,0,0,0);
	  cells[467]->set_vertices(vertices[80],vertices[33],vertices[25],vertices[73]);
	  cells[467]->set_neighbors(cells[465],cells[305],cells[556],cells[402]);
	  set_offsets(cells[467],0,0,0,0);
	  cells[468]->set_vertices(vertices[262],vertices[214],vertices[254],vertices[207]);
	  cells[468]->set_neighbors(cells[1493],cells[1538],cells[1350],cells[1483]);
	  set_offsets(cells[468],0,0,0,0);
	  cells[469]->set_vertices(vertices[67],vertices[114],vertices[107],vertices[59]);
	  cells[469]->set_neighbors(cells[282],cells[154],cells[82],cells[742]);
	  set_offsets(cells[469],0,0,0,0);
	  cells[470]->set_vertices(vertices[66],vertices[19],vertices[59],vertices[67]);
	  cells[470]->set_neighbors(cells[386],cells[82],cells[32],cells[318]);
	  set_offsets(cells[470],0,0,0,0);
	  cells[471]->set_vertices(vertices[156],vertices[148],vertices[101],vertices[108]);
	  cells[471]->set_neighbors(cells[734],cells[705],cells[720],cells[417]);
	  set_offsets(cells[471],0,0,0,0);
	  cells[472]->set_vertices(vertices[132],vertices[124],vertices[125],vertices[77]);
	  cells[472]->set_neighbors(cells[796],cells[795],cells[365],cells[928]);
	  set_offsets(cells[472],0,0,0,0);
	  cells[473]->set_vertices(vertices[3],vertices[51],vertices[91],vertices[50]);
	  cells[473]->set_neighbors(cells[187],cells[519],cells[249],cells[592]);
	  set_offsets(cells[473],2,2,0,2);
	  cells[474]->set_vertices(vertices[123],vertices[122],vertices[115],vertices[75]);
	  cells[474]->set_neighbors(cells[379],cells[548],cells[594],cells[924]);
	  set_offsets(cells[474],0,0,0,0);
	  cells[475]->set_vertices(vertices[112],vertices[104],vertices[63],vertices[111]);
	  cells[475]->set_neighbors(cells[698],cells[682],cells[680],cells[726]);
	  set_offsets(cells[475],1,1,0,0);
	  cells[476]->set_vertices(vertices[52],vertices[45],vertices[4],vertices[92]);
	  cells[476]->set_neighbors(cells[575],cells[525],cells[597],cells[276]);
	  set_offsets(cells[476],2,0,2,0);
	  cells[477]->set_vertices(vertices[90],vertices[43],vertices[2],vertices[42]);
	  cells[477]->set_neighbors(cells[931],cells[52],cells[313],cells[536]);
	  set_offsets(cells[477],0,0,2,0);
	  cells[478]->set_vertices(vertices[108],vertices[116],vertices[107],vertices[59]);
	  cells[478]->set_neighbors(cells[154],cells[630],cells[6],cells[756]);
	  set_offsets(cells[478],0,0,0,0);
	  cells[479]->set_vertices(vertices[60],vertices[68],vertices[61],vertices[108]);
	  cells[479]->set_neighbors(cells[193],cells[247],cells[498],cells[274]);
	  set_offsets(cells[479],0,0,0,0);
	  cells[480]->set_vertices(vertices[76],vertices[68],vertices[19],vertices[67]);
	  cells[480]->set_neighbors(cells[386],cells[283],cells[181],cells[86]);
	  set_offsets(cells[480],0,0,0,0);
	  cells[481]->set_vertices(vertices[28],vertices[76],vertices[21],vertices[68]);
	  cells[481]->set_neighbors(cells[22],cells[104],cells[86],cells[439]);
	  set_offsets(cells[481],0,0,0,0);
	  cells[482]->set_vertices(vertices[82],vertices[27],vertices[34],vertices[74]);
	  cells[482]->set_neighbors(cells[186],cells[520],cells[40],cells[310]);
	  set_offsets(cells[482],0,0,0,0);
	  cells[483]->set_vertices(vertices[98],vertices[106],vertices[97],vertices[49]);
	  cells[483]->set_neighbors(cells[352],cells[860],cells[260],cells[626]);
	  set_offsets(cells[483],0,0,0,0);
	  cells[484]->set_vertices(vertices[26],vertices[74],vertices[19],vertices[66]);
	  cells[484]->set_neighbors(cells[32],cells[370],cells[485],cells[285]);
	  set_offsets(cells[484],0,0,0,0);
	  cells[485]->set_vertices(vertices[74],vertices[26],vertices[17],vertices[66]);
	  cells[485]->set_neighbors(cells[361],cells[486],cells[484],cells[242]);
	  set_offsets(cells[485],0,0,0,0);
	  cells[486]->set_vertices(vertices[74],vertices[66],vertices[17],vertices[65]);
	  cells[486]->set_neighbors(cells[338],cells[90],cells[240],cells[485]);
	  set_offsets(cells[486],0,0,0,0);
	  cells[487]->set_vertices(vertices[110],vertices[63],vertices[62],vertices[55]);
	  cells[487]->set_neighbors(cells[21],cells[23],cells[699],cells[373]);
	  set_offsets(cells[487],0,0,0,0);
	  cells[488]->set_vertices(vertices[68],vertices[69],vertices[21],vertices[61]);
	  cells[488]->set_neighbors(cells[201],cells[243],cells[45],cells[22]);
	  set_offsets(cells[488],0,0,0,0);
	  cells[489]->set_vertices(vertices[144],vertices[96],vertices[103],vertices[143]);
	  cells[489]->set_neighbors(cells[838],cells[892],cells[1120],cells[638]);
	  set_offsets(cells[489],3,3,2,0);
	  cells[490]->set_vertices(vertices[180],vertices[188],vertices[133],vertices[181]);
	  cells[490]->set_neighbors(cells[825],cells[1107],cells[1325],cells[1123]);
	  set_offsets(cells[490],0,0,0,0);
	  cells[491]->set_vertices(vertices[73],vertices[74],vertices[122],vertices[82]);
	  cells[491]->set_neighbors(cells[697],cells[778],cells[377],cells[616]);
	  set_offsets(cells[491],0,0,0,0);
	  cells[492]->set_vertices(vertices[100],vertices[91],vertices[140],vertices[52]);
	  cells[492]->set_neighbors(cells[349],cells[788],cells[162],cells[827]);
	  set_offsets(cells[492],2,0,0,2);
	  cells[493]->set_vertices(vertices[35],vertices[92],vertices[83],vertices[43]);
	  cells[493]->set_neighbors(cells[211],cells[329],cells[41],cells[302]);
	  set_offsets(cells[493],0,0,0,0);
	  cells[494]->set_vertices(vertices[67],vertices[75],vertices[124],vertices[115]);
	  cells[494]->set_neighbors(cells[548],cells[663],cells[379],cells[624]);
	  set_offsets(cells[494],0,0,0,0);
	  cells[495]->set_vertices(vertices[31],vertices[32],vertices[80],vertices[39]);
	  cells[495]->set_neighbors(cells[539],cells[251],cells[0],cells[508]);
	  set_offsets(cells[495],0,1,1,0);
	  cells[496]->set_vertices(vertices[77],vertices[37],vertices[86],vertices[85]);
	  cells[496]->set_neighbors(cells[407],cells[311],cells[568],cells[518]);
	  set_offsets(cells[496],0,0,0,0);
	  cells[497]->set_vertices(vertices[42],vertices[35],vertices[34],vertices[82]);
	  cells[497]->set_neighbors(cells[310],cells[449],cells[118],cells[1618]);
	  set_offsets(cells[497],0,0,0,0);
	  cells[498]->set_vertices(vertices[60],vertices[68],vertices[108],vertices[59]);
	  cells[498]->set_neighbors(cells[6],cells[7],cells[298],cells[479]);
	  set_offsets(cells[498],0,0,0,0);
	  cells[499]->set_vertices(vertices[31],vertices[24],vertices[72],vertices[32]);
	  cells[499]->set_neighbors(cells[92],cells[508],cells[1413],cells[2]);
	  set_offsets(cells[499],0,1,1,1);
	  cells[500]->set_vertices(vertices[86],vertices[39],vertices[94],vertices[46]);
	  cells[500]->set_neighbors(cells[384],cells[555],cells[315],cells[596]);
	  set_offsets(cells[500],0,0,0,0);
	  cells[501]->set_vertices(vertices[114],vertices[154],vertices[107],vertices[106]);
	  cells[501]->set_neighbors(cells[541],cells[282],cells[950],cells[622]);
	  set_offsets(cells[501],0,0,0,0);
	  cells[502]->set_vertices(vertices[84],vertices[29],vertices[36],vertices[76]);
	  cells[502]->set_neighbors(cells[107],cells[220],cells[526],cells[119]);
	  set_offsets(cells[502],0,0,0,0);
	  cells[503]->set_vertices(vertices[35],vertices[84],vertices[36],vertices[27]);
	  cells[503]->set_neighbors(cells[220],cells[1424],cells[565],cells[538]);
	  set_offsets(cells[503],0,0,0,0);
	  cells[504]->set_vertices(vertices[65],vertices[114],vertices[105],vertices[113]);
	  cells[504]->set_neighbors(cells[415],cells[729],cells[231],cells[450]);
	  set_offsets(cells[504],0,0,0,0);
	  cells[505]->set_vertices(vertices[72],vertices[73],vertices[25],vertices[65]);
	  cells[505]->set_neighbors(cells[358],cells[466],cells[738],cells[305]);
	  set_offsets(cells[505],0,0,0,0);
	  cells[506]->set_vertices(vertices[166],vertices[158],vertices[159],vertices[111]);
	  cells[506]->set_neighbors(cells[585],cells[750],cells[690],cells[1031]);
	  set_offsets(cells[506],0,0,0,0);
	  cells[507]->set_vertices(vertices[88],vertices[40],vertices[33],vertices[41]);
	  cells[507]->set_neighbors(cells[1653],cells[403],cells[235],cells[562]);
	  set_offsets(cells[507],0,0,0,0);
	  cells[508]->set_vertices(vertices[31],vertices[32],vertices[72],vertices[80]);
	  cells[508]->set_neighbors(cells[437],cells[38],cells[495],cells[499]);
	  set_offsets(cells[508],0,1,1,1);
	  cells[509]->set_vertices(vertices[104],vertices[57],vertices[112],vertices[64]);
	  cells[509]->set_neighbors(cells[569],cells[726],cells[215],cells[227]);
	  set_offsets(cells[509],0,0,0,0);
	  cells[510]->set_vertices(vertices[65],vertices[57],vertices[112],vertices[105]);
	  cells[510]->set_neighbors(cells[227],cells[729],cells[450],cells[569]);
	  set_offsets(cells[510],0,0,0,0);
	  cells[511]->set_vertices(vertices[99],vertices[51],vertices[106],vertices[98]);
	  cells[511]->set_neighbors(cells[312],cells[595],cells[872],cells[635]);
	  set_offsets(cells[511],0,0,0,0);
	  cells[512]->set_vertices(vertices[113],vertices[105],vertices[153],vertices[162]);
	  cells[512]->set_neighbors(cells[969],cells[746],cells[415],cells[955]);
	  set_offsets(cells[512],0,0,0,0);
	  cells[513]->set_vertices(vertices[92],vertices[37],vertices[85],vertices[45]);
	  cells[513]->set_neighbors(cells[573],cells[598],cells[535],cells[129]);
	  set_offsets(cells[513],0,0,0,0);
	  cells[514]->set_vertices(vertices[82],vertices[27],vertices[75],vertices[35]);
	  cells[514]->set_neighbors(cells[565],cells[566],cells[310],cells[40]);
	  set_offsets(cells[514],0,0,0,0);
	  cells[515]->set_vertices(vertices[74],vertices[27],vertices[67],vertices[75]);
	  cells[515]->set_neighbors(cells[205],cells[644],cells[40],cells[432]);
	  set_offsets(cells[515],0,0,0,0);
	  cells[516]->set_vertices(vertices[129],vertices[128],vertices[176],vertices[121]);
	  cells[516]->set_neighbors(cells[990],cells[713],cells[617],cells[1073]);
	  set_offsets(cells[516],0,0,0,0);
	  cells[517]->set_vertices(vertices[122],vertices[162],vertices[115],vertices[114]);
	  cells[517]->set_neighbors(cells[790],cells[20],cells[748],cells[841]);
	  set_offsets(cells[517],0,0,0,0);
	  cells[518]->set_vertices(vertices[29],vertices[86],vertices[77],vertices[37]);
	  cells[518]->set_neighbors(cells[496],cells[527],cells[553],cells[279]);
	  set_offsets(cells[518],0,0,0,0);
	  cells[519]->set_vertices(vertices[3],vertices[91],vertices[43],vertices[50]);
	  cells[519]->set_neighbors(cells[156],cells[374],cells[473],cells[563]);
	  set_offsets(cells[519],2,0,0,2);
	  cells[520]->set_vertices(vertices[82],vertices[74],vertices[34],vertices[25]);
	  cells[520]->set_neighbors(cells[30],cells[74],cells[377],cells[482]);
	  set_offsets(cells[520],0,0,0,0);
	  cells[521]->set_vertices(vertices[130],vertices[121],vertices[170],vertices[122]);
	  cells[521]->set_neighbors(cells[792],cells[854],cells[115],cells[408]);
	  set_offsets(cells[521],0,0,0,0);
	  cells[522]->set_vertices(vertices[268],vertices[260],vertices[261],vertices[213]);
	  cells[522]->set_neighbors(cells[1344],cells[891],cells[262],cells[217]);
	  set_offsets(cells[522],0,0,0,0);
	  cells[523]->set_vertices(vertices[43],vertices[92],vertices[4],vertices[44]);
	  cells[523]->set_neighbors(cells[575],cells[1364],cells[41],cells[525]);
	  set_offsets(cells[523],0,0,2,0);
	  cells[524]->set_vertices(vertices[84],vertices[37],vertices[92],vertices[44]);
	  cells[524]->set_neighbors(cells[535],cells[143],cells[368],cells[129]);
	  set_offsets(cells[524],0,0,0,0);
	  cells[525]->set_vertices(vertices[43],vertices[52],vertices[4],vertices[92]);
	  cells[525]->set_neighbors(cells[476],cells[523],cells[167],cells[331]);
	  set_offsets(cells[525],0,2,2,0);
	  cells[526]->set_vertices(vertices[76],vertices[29],vertices[77],vertices[84]);
	  cells[526]->set_neighbors(cells[527],cells[545],cells[502],cells[528]);
	  set_offsets(cells[526],0,0,0,0);
	  cells[527]->set_vertices(vertices[84],vertices[29],vertices[77],vertices[37]);
	  cells[527]->set_neighbors(cells[518],cells[568],cells[119],cells[526]);
	  set_offsets(cells[527],0,0,0,0);
	  cells[528]->set_vertices(vertices[76],vertices[29],vertices[69],vertices[77]);
	  cells[528]->set_neighbors(cells[71],cells[666],cells[526],cells[300]);
	  set_offsets(cells[528],0,0,0,0);
	  cells[529]->set_vertices(vertices[274],vertices[226],vertices[266],vertices[219]);
	  cells[529]->set_neighbors(cells[1373],cells[1152],cells[916],cells[1178]);
	  set_offsets(cells[529],0,0,0,0);
	  cells[530]->set_vertices(vertices[94],vertices[47],vertices[39],vertices[87]);
	  cells[530]->set_neighbors(cells[554],cells[596],cells[606],cells[384]);
	  set_offsets(cells[530],0,0,0,0);
	  cells[531]->set_vertices(vertices[136],vertices[95],vertices[135],vertices[87]);
	  cells[531]->set_neighbors(cells[848],cells[174],cells[806],cells[672]);
	  set_offsets(cells[531],1,0,0,0);
	  cells[532]->set_vertices(vertices[85],vertices[94],vertices[93],vertices[45]);
	  cells[532]->set_neighbors(cells[255],cells[598],cells[573],cells[667]);
	  set_offsets(cells[532],0,0,0,0);
	  cells[533]->set_vertices(vertices[119],vertices[112],vertices[71],vertices[111]);
	  cells[533]->set_neighbors(cells[682],cells[716],cells[762],cells[731]);
	  set_offsets(cells[533],0,1,0,0);
	  cells[534]->set_vertices(vertices[80],vertices[33],vertices[81],vertices[88]);
	  cells[534]->set_neighbors(cells[403],cells[773],cells[562],cells[556]);
	  set_offsets(cells[534],0,0,0,0);
	  cells[535]->set_vertices(vertices[92],vertices[37],vertices[45],vertices[44]);
	  cells[535]->set_neighbors(cells[1354],cells[575],cells[524],cells[513]);
	  set_offsets(cells[535],0,0,0,0);
	  cells[536]->set_vertices(vertices[50],vertices[43],vertices[2],vertices[90]);
	  cells[536]->set_neighbors(cells[477],cells[564],cells[156],cells[374]);
	  set_offsets(cells[536],2,0,2,0);
	  cells[537]->set_vertices(vertices[83],vertices[75],vertices[84],vertices[35]);
	  cells[537]->set_neighbors(cells[565],cells[302],cells[566],cells[633]);
	  set_offsets(cells[537],0,0,0,0);
	  cells[538]->set_vertices(vertices[44],vertices[84],vertices[36],vertices[35]);
	  cells[538]->set_neighbors(cells[503],cells[267],cells[143],cells[368]);
	  set_offsets(cells[538],0,0,0,0);
	  cells[539]->set_vertices(vertices[39],vertices[32],vertices[80],vertices[40]);
	  cells[539]->set_neighbors(cells[160],cells[163],cells[1220],cells[495]);
	  set_offsets(cells[539],0,1,1,1);
	  cells[540]->set_vertices(vertices[173],vertices[125],vertices[165],vertices[174]);
	  cells[540]->set_neighbors(cells[982],cells[246],cells[589],cells[953]);
	  set_offsets(cells[540],0,0,0,0);
	  cells[541]->set_vertices(vertices[107],vertices[154],vertices[99],vertices[106]);
	  cells[541]->set_neighbors(cells[754],cells[649],cells[501],cells[704]);
	  set_offsets(cells[541],0,0,0,0);
	  cells[542]->set_vertices(vertices[176],vertices[184],vertices[129],vertices[177]);
	  cells[542]->set_neighbors(cells[853],cells[1043],cells[1084],cells[965]);
	  set_offsets(cells[542],0,0,0,0);
	  cells[543]->set_vertices(vertices[140],vertices[93],vertices[85],vertices[133]);
	  cells[543]->set_neighbors(cells[834],cells[625],cells[884],cells[812]);
	  set_offsets(cells[543],0,0,0,0);
	  cells[544]->set_vertices(vertices[126],vertices[79],vertices[78],vertices[71]);
	  cells[544]->set_neighbors(cells[393],cells[362],cells[446],cells[440]);
	  set_offsets(cells[544],0,0,0,0);
	  cells[545]->set_vertices(vertices[124],vertices[76],vertices[77],vertices[84]);
	  cells[545]->set_neighbors(cells[526],cells[365],cells[80],cells[666]);
	  set_offsets(cells[545],0,0,0,0);
	  cells[546]->set_vertices(vertices[87],vertices[88],vertices[39],vertices[80]);
	  cells[546]->set_neighbors(cells[163],cells[577],cells[745],cells[554]);
	  set_offsets(cells[546],0,1,0,1);
	  cells[547]->set_vertices(vertices[61],vertices[69],vertices[118],vertices[109]);
	  cells[547]->set_neighbors(cells[580],cells[284],cells[317],cells[722]);
	  set_offsets(cells[547],0,0,0,0);
	  cells[548]->set_vertices(vertices[75],vertices[124],vertices[115],vertices[123]);
	  cells[548]->set_neighbors(cells[845],cells[474],cells[418],cells[494]);
	  set_offsets(cells[548],0,0,0,0);
	  cells[549]->set_vertices(vertices[68],vertices[67],vertices[116],vertices[59]);
	  cells[549]->set_neighbors(cells[154],cells[6],cells[386],cells[181]);
	  set_offsets(cells[549],0,0,0,0);
	  cells[550]->set_vertices(vertices[39],vertices[40],vertices[88],vertices[47]);
	  cells[550]->set_neighbors(cells[93],cells[554],cells[1595],cells[163]);
	  set_offsets(cells[550],0,1,1,0);
	  cells[551]->set_vertices(vertices[31],vertices[78],vertices[38],vertices[86]);
	  cells[551]->set_neighbors(cells[359],cells[552],cells[435],cells[464]);
	  set_offsets(cells[551],0,0,0,0);
	  cells[552]->set_vertices(vertices[31],vertices[86],vertices[38],vertices[39]);
	  cells[552]->set_neighbors(cells[315],cells[1015],cells[387],cells[551]);
	  set_offsets(cells[552],0,0,0,0);
	  cells[553]->set_vertices(vertices[29],vertices[38],vertices[86],vertices[37]);
	  cells[553]->set_neighbors(cells[275],cells[518],cells[1246],cells[359]);
	  set_offsets(cells[553],0,0,0,0);
	  cells[554]->set_vertices(vertices[87],vertices[47],vertices[39],vertices[88]);
	  cells[554]->set_neighbors(cells[550],cells[546],cells[607],cells[530]);
	  set_offsets(cells[554],0,0,0,1);
	  cells[555]->set_vertices(vertices[94],vertices[86],vertices[46],vertices[37]);
	  cells[555]->set_neighbors(cells[275],cells[219],cells[407],cells[500]);
	  set_offsets(cells[555],0,0,0,0);
	  cells[556]->set_vertices(vertices[80],vertices[33],vertices[73],vertices[81]);
	  cells[556]->set_neighbors(cells[60],cells[774],cells[534],cells[467]);
	  set_offsets(cells[556],0,0,0,0);
	  cells[557]->set_vertices(vertices[145],vertices[154],vertices[105],vertices[97]);
	  cells[557]->set_neighbors(cells[366],cells[906],cells[457],cells[73]);
	  set_offsets(cells[557],0,0,0,0);
	  cells[558]->set_vertices(vertices[149],vertices[206],vertices[197],vertices[157]);
	  cells[558]->set_neighbors(cells[460],cells[35],cells[692],cells[77]);
	  set_offsets(cells[558],0,0,0,0);
	  cells[559]->set_vertices(vertices[257],vertices[18],vertices[17],vertices[26]);
	  cells[559]->set_neighbors(cells[361],cells[1613],cells[1560],cells[1555]);
	  set_offsets(cells[559],0,4,4,4);
	  cells[560]->set_vertices(vertices[33],vertices[90],vertices[81],vertices[41]);
	  cells[560]->set_neighbors(cells[198],cells[403],cells[102],cells[290]);
	  set_offsets(cells[560],0,0,0,0);
	  cells[561]->set_vertices(vertices[69],vertices[77],vertices[78],vertices[126]);
	  cells[561]->set_neighbors(cells[442],cells[640],cells[655],cells[71]);
	  set_offsets(cells[561],0,0,0,0);
	  cells[562]->set_vertices(vertices[80],vertices[40],vertices[33],vertices[88]);
	  cells[562]->set_neighbors(cells[507],cells[534],cells[163],cells[160]);
	  set_offsets(cells[562],0,0,0,0);
	  cells[563]->set_vertices(vertices[3],vertices[91],vertices[52],vertices[43]);
	  cells[563]->set_neighbors(cells[167],cells[331],cells[519],cells[592]);
	  set_offsets(cells[563],2,0,2,0);
	  cells[564]->set_vertices(vertices[41],vertices[50],vertices[2],vertices[90]);
	  cells[564]->set_neighbors(cells[536],cells[52],cells[399],cells[328]);
	  set_offsets(cells[564],0,2,2,0);
	  cells[565]->set_vertices(vertices[35],vertices[75],vertices[84],vertices[27]);
	  cells[565]->set_neighbors(cells[458],cells[503],cells[514],cells[537]);
	  set_offsets(cells[565],0,0,0,0);
	  cells[566]->set_vertices(vertices[82],vertices[35],vertices[75],vertices[83]);
	  cells[566]->set_neighbors(cells[537],cells[68],cells[293],cells[514]);
	  set_offsets(cells[566],0,0,0,0);
	  cells[567]->set_vertices(vertices[48],vertices[88],vertices[41],vertices[0]);
	  cells[567]->set_neighbors(cells[235],cells[1],cells[363],cells[70]);
	  set_offsets(cells[567],2,0,0,2);
	  cells[568]->set_vertices(vertices[84],vertices[37],vertices[77],vertices[85]);
	  cells[568]->set_neighbors(cells[496],cells[579],cells[129],cells[527]);
	  set_offsets(cells[568],0,0,0,0);
	  cells[569]->set_vertices(vertices[112],vertices[57],vertices[65],vertices[64]);
	  cells[569]->set_neighbors(cells[454],cells[180],cells[509],cells[510]);
	  set_offsets(cells[569],0,0,0,0);
	  cells[570]->set_vertices(vertices[265],vertices[257],vertices[26],vertices[266]);
	  cells[570]->set_neighbors(cells[1560],cells[1385],cells[1059],cells[1613]);
	  set_offsets(cells[570],0,0,4,0);
	  cells[571]->set_vertices(vertices[176],vertices[121],vertices[168],vertices[169]);
	  cells[571]->set_neighbors(cells[1044],cells[16],cells[713],cells[990]);
	  set_offsets(cells[571],0,0,0,0);
	  cells[572]->set_vertices(vertices[132],vertices[123],vertices[172],vertices[124]);
	  cells[572]->set_neighbors(cells[845],cells[928],cells[418],cells[660]);
	  set_offsets(cells[572],0,0,0,0);
	  cells[573]->set_vertices(vertices[37],vertices[94],vertices[85],vertices[45]);
	  cells[573]->set_neighbors(cells[532],cells[513],cells[219],cells[407]);
	  set_offsets(cells[573],0,0,0,0);
	  cells[574]->set_vertices(vertices[5],vertices[93],vertices[54],vertices[45]);
	  cells[574]->set_neighbors(cells[255],cells[12],cells[108],cells[600]);
	  set_offsets(cells[574],2,0,2,0);
	  cells[575]->set_vertices(vertices[92],vertices[45],vertices[4],vertices[44]);
	  cells[575]->set_neighbors(cells[1471],cells[523],cells[535],cells[476]);
	  set_offsets(cells[575],0,0,2,0);
	  cells[576]->set_vertices(vertices[95],vertices[48],vertices[96],vertices[55]);
	  cells[576]->set_neighbors(cells[200],cells[899],cells[610],cells[839]);
	  set_offsets(cells[576],0,3,3,2);
	  cells[577]->set_vertices(vertices[79],vertices[87],vertices[39],vertices[80]);
	  cells[577]->set_neighbors(cells[546],cells[251],cells[661],cells[578]);
	  set_offsets(cells[577],0,0,0,1);
	  cells[578]->set_vertices(vertices[86],vertices[87],vertices[39],vertices[79]);
	  cells[578]->set_neighbors(cells[577],cells[387],cells[55],cells[596]);
	  set_offsets(cells[578],0,0,0,0);
	  cells[579]->set_vertices(vertices[85],vertices[84],vertices[132],vertices[77]);
	  cells[579]->set_neighbors(cells[365],cells[795],cells[568],cells[632]);
	  set_offsets(cells[579],0,0,0,0);
	  cells[580]->set_vertices(vertices[69],vertices[117],vertices[118],vertices[109]);
	  cells[580]->set_neighbors(cells[224],cells[547],cells[752],cells[603]);
	  set_offsets(cells[580],0,0,0,0);
	  cells[581]->set_vertices(vertices[45],vertices[94],vertices[6],vertices[46]);
	  cells[581]->set_neighbors(cells[582],cells[1628],cells[219],cells[175]);
	  set_offsets(cells[581],0,0,2,0);
	  cells[582]->set_vertices(vertices[94],vertices[47],vertices[6],vertices[46]);
	  cells[582]->set_neighbors(cells[1488],cells[581],cells[384],cells[426]);
	  set_offsets(cells[582],0,0,2,0);
	  cells[583]->set_vertices(vertices[1],vertices[49],vertices[50],vertices[89]);
	  cells[583]->set_neighbors(cells[631],cells[289],cells[584],cells[49]);
	  set_offsets(cells[583],2,2,2,0);
	  cells[584]->set_vertices(vertices[1],vertices[49],vertices[89],vertices[48]);
	  cells[584]->set_neighbors(cells[111],cells[123],cells[97],cells[583]);
	  set_offsets(cells[584],2,2,0,2);
	  cells[585]->set_vertices(vertices[159],vertices[158],vertices[151],vertices[111]);
	  cells[585]->set_neighbors(cells[409],cells[241],cells[506],cells[1111]);
	  set_offsets(cells[585],0,0,0,0);
	  cells[586]->set_vertices(vertices[281],vertices[40],vertices[33],vertices[273]);
	  cells[586]->set_neighbors(cells[1439],cells[628],cells[1651],cells[1653]);
	  set_offsets(cells[586],0,4,4,0);
	  cells[587]->set_vertices(vertices[255],vertices[16],vertices[8],vertices[15]);
	  cells[587]->set_neighbors(cells[385],cells[1304],cells[1491],cells[1266]);
	  set_offsets(cells[587],0,5,5,4);
	  cells[588]->set_vertices(vertices[97],vertices[49],vertices[57],vertices[104]);
	  cells[588]->set_neighbors(cells[395],cells[683],cells[636],cells[352]);
	  set_offsets(cells[588],0,0,0,0);
	  cells[589]->set_vertices(vertices[182],vertices[125],vertices[173],vertices[174]);
	  cells[589]->set_neighbors(cells[540],cells[1071],cells[1070],cells[846]);
	  set_offsets(cells[589],0,0,0,0);
	  cells[590]->set_vertices(vertices[133],vertices[125],vertices[182],vertices[134]);
	  cells[590]->set_neighbors(cells[1070],cells[1105],cells[109],cells[846]);
	  set_offsets(cells[590],0,0,0,0);
	  cells[591]->set_vertices(vertices[130],vertices[138],vertices[81],vertices[90]);
	  cells[591]->set_neighbors(cells[789],cells[794],cells[776],cells[78]);
	  set_offsets(cells[591],0,0,0,0);
	  cells[592]->set_vertices(vertices[3],vertices[51],vertices[52],vertices[91]);
	  cells[592]->set_neighbors(cells[162],cells[563],cells[473],cells[223]);
	  set_offsets(cells[592],2,2,2,0);
	  cells[593]->set_vertices(vertices[104],vertices[97],vertices[152],vertices[105]);
	  cells[593]->set_neighbors(cells[906],cells[842],cells[683],cells[797]);
	  set_offsets(cells[593],0,0,0,0);
	  cells[594]->set_vertices(vertices[130],vertices[122],vertices[123],vertices[75]);
	  cells[594]->set_neighbors(cells[474],cells[781],cells[814],cells[854]);
	  set_offsets(cells[594],0,0,0,0);
	  cells[595]->set_vertices(vertices[99],vertices[98],vertices[106],vertices[146]);
	  cells[595]->set_neighbors(cells[626],cells[754],cells[824],cells[511]);
	  set_offsets(cells[595],0,0,0,0);
	  cells[596]->set_vertices(vertices[86],vertices[94],vertices[39],vertices[87]);
	  cells[596]->set_neighbors(cells[530],cells[578],cells[218],cells[500]);
	  set_offsets(cells[596],0,0,0,0);
	  cells[597]->set_vertices(vertices[52],vertices[92],vertices[93],vertices[45]);
	  cells[597]->set_neighbors(cells[598],cells[108],cells[476],cells[50]);
	  set_offsets(cells[597],2,0,0,0);
	  cells[598]->set_vertices(vertices[93],vertices[92],vertices[85],vertices[45]);
	  cells[598]->set_neighbors(cells[513],cells[532],cells[597],cells[812]);
	  set_offsets(cells[598],0,0,0,0);
	  cells[599]->set_vertices(vertices[5],vertices[53],vertices[93],vertices[52]);
	  cells[599]->set_neighbors(cells[651],cells[108],cells[149],cells[600]);
	  set_offsets(cells[599],2,2,0,2);
	  cells[600]->set_vertices(vertices[5],vertices[53],vertices[54],vertices[93]);
	  cells[600]->set_neighbors(cells[668],cells[574],cells[599],cells[207]);
	  set_offsets(cells[600],2,2,2,0);
	  cells[601]->set_vertices(vertices[73],vertices[82],vertices[130],vertices[81]);
	  cells[601]->set_neighbors(cells[794],cells[764],cells[60],cells[778]);
	  set_offsets(cells[601],0,0,0,0);
	  cells[602]->set_vertices(vertices[185],vertices[145],vertices[146],vertices[194]);
	  cells[602]->set_neighbors(cells[719],cells[1414],cells[1408],cells[1091]);
	  set_offsets(cells[602],0,2,2,2);
	  cells[603]->set_vertices(vertices[69],vertices[126],vertices[118],vertices[117]);
	  cells[603]->set_neighbors(cells[844],cells[580],cells[655],cells[640]);
	  set_offsets(cells[603],0,0,0,0);
	  cells[604]->set_vertices(vertices[96],vertices[56],vertices[48],vertices[49]);
	  cells[604]->set_neighbors(cells[97],cells[111],cells[637],cells[200]);
	  set_offsets(cells[604],0,0,0,0);
	  cells[605]->set_vertices(vertices[54],vertices[47],vertices[94],vertices[95]);
	  cells[605]->set_neighbors(cells[606],cells[351],cells[612],cells[426]);
	  set_offsets(cells[605],2,0,0,0);
	  cells[606]->set_vertices(vertices[95],vertices[47],vertices[94],vertices[87]);
	  cells[606]->set_neighbors(cells[530],cells[849],cells[607],cells[605]);
	  set_offsets(cells[606],0,0,0,0);
	  cells[607]->set_vertices(vertices[87],vertices[47],vertices[88],vertices[95]);
	  cells[607]->set_neighbors(cells[608],cells[806],cells[606],cells[554]);
	  set_offsets(cells[607],0,0,1,0);
	  cells[608]->set_vertices(vertices[95],vertices[47],vertices[88],vertices[48]);
	  cells[608]->set_neighbors(cells[363],cells[727],cells[609],cells[607]);
	  set_offsets(cells[608],0,0,1,3);
	  cells[609]->set_vertices(vertices[95],vertices[7],vertices[47],vertices[48]);
	  cells[609]->set_neighbors(cells[27],cells[608],cells[610],cells[612]);
	  set_offsets(cells[609],0,2,0,3);
	  cells[610]->set_vertices(vertices[55],vertices[7],vertices[95],vertices[48]);
	  cells[610]->set_neighbors(cells[609],cells[576],cells[230],cells[611]);
	  set_offsets(cells[610],2,2,0,3);
	  cells[611]->set_vertices(vertices[54],vertices[7],vertices[95],vertices[55]);
	  cells[611]->set_neighbors(cells[610],cells[669],cells[114],cells[612]);
	  set_offsets(cells[611],2,2,0,2);
	  cells[612]->set_vertices(vertices[54],vertices[7],vertices[47],vertices[95]);
	  cells[612]->set_neighbors(cells[609],cells[605],cells[611],cells[103]);
	  set_offsets(cells[612],2,2,0,0);
	  cells[613]->set_vertices(vertices[239],vertices[191],vertices[232],vertices[192]);
	  cells[613]->set_neighbors(cells[1259],cells[1394],cells[345],cells[1437]);
	  set_offsets(cells[613],0,0,1,3);
	  cells[614]->set_vertices(vertices[145],vertices[97],vertices[152],vertices[144]);
	  cells[614]->set_neighbors(cells[797],cells[686],cells[1130],cells[906]);
	  set_offsets(cells[614],0,0,0,0);
	  cells[615]->set_vertices(vertices[119],vertices[166],vertices[111],vertices[118]);
	  cells[615]->set_neighbors(cells[690],cells[716],cells[1013],cells[750]);
	  set_offsets(cells[615],0,0,0,0);
	  cells[616]->set_vertices(vertices[65],vertices[74],vertices[122],vertices[73]);
	  cells[616]->set_neighbors(cells[491],cells[618],cells[358],cells[353]);
	  set_offsets(cells[616],0,0,0,0);
	  cells[617]->set_vertices(vertices[128],vertices[129],vertices[81],vertices[121]);
	  cells[617]->set_neighbors(cells[98],cells[394],cells[516],cells[357]);
	  set_offsets(cells[617],0,0,0,0);
	  cells[618]->set_vertices(vertices[65],vertices[122],vertices[113],vertices[73]);
	  cells[618]->set_neighbors(cells[94],cells[19],cells[616],cells[231]);
	  set_offsets(cells[618],0,0,0,0);
	  cells[619]->set_vertices(vertices[67],vertices[124],vertices[76],vertices[116]);
	  cells[619]->set_neighbors(cells[375],cells[181],cells[663],cells[624]);
	  set_offsets(cells[619],0,0,0,0);
	  cells[620]->set_vertices(vertices[188],vertices[139],vertices[131],vertices[140]);
	  cells[620]->set_neighbors(cells[836],cells[621],cells[1096],cells[1061]);
	  set_offsets(cells[620],0,0,0,0);
	  cells[621]->set_vertices(vertices[180],vertices[188],vertices[131],vertices[140]);
	  cells[621]->set_neighbors(cells[620],cells[1068],cells[1123],cells[896]);
	  set_offsets(cells[621],0,0,0,0);
	  cells[622]->set_vertices(vertices[162],vertices[154],vertices[107],vertices[114]);
	  cells[622]->set_neighbors(cells[501],cells[790],cells[793],cells[770]);
	  set_offsets(cells[622],0,0,0,0);
	  cells[623]->set_vertices(vertices[286],vertices[237],vertices[277],vertices[229]);
	  cells[623]->set_neighbors(cells[1003],cells[1480],cells[1627],cells[1635]);
	  set_offsets(cells[623],0,0,0,0);
	  cells[624]->set_vertices(vertices[67],vertices[75],vertices[76],vertices[124]);
	  cells[624]->set_neighbors(cells[80],cells[619],cells[494],cells[205]);
	  set_offsets(cells[624],0,0,0,0);
	  cells[625]->set_vertices(vertices[132],vertices[140],vertices[85],vertices[133]);
	  cells[625]->set_neighbors(cells[543],cells[828],cells[927],cells[871]);
	  set_offsets(cells[625],0,0,0,0);
	  cells[626]->set_vertices(vertices[146],vertices[98],vertices[106],vertices[97]);
	  cells[626]->set_neighbors(cells[483],cells[333],cells[915],cells[595]);
	  set_offsets(cells[626],0,0,0,0);
	  cells[627]->set_vertices(vertices[53],vertices[62],vertices[102],vertices[110]);
	  cells[627]->set_neighbors(cells[23],cells[427],cells[676],cells[411]);
	  set_offsets(cells[627],0,0,0,0);
	  cells[628]->set_vertices(vertices[281],vertices[273],vertices[33],vertices[42]);
	  cells[628]->set_neighbors(cells[1438],cells[1183],cells[1689],cells[586]);
	  set_offsets(cells[628],0,0,4,4);
	  cells[629]->set_vertices(vertices[110],vertices[118],vertices[70],vertices[63]);
	  cells[629]->set_neighbors(cells[294],cells[373],cells[225],cells[204]);
	  set_offsets(cells[629],0,0,0,0);
	  cells[630]->set_vertices(vertices[108],vertices[107],vertices[99],vertices[59]);
	  cells[630]->set_neighbors(cells[649],cells[210],cells[478],cells[907]);
	  set_offsets(cells[630],0,0,0,0);
	  cells[631]->set_vertices(vertices[98],vertices[49],vertices[89],vertices[50]);
	  cells[631]->set_neighbors(cells[583],cells[855],cells[320],cells[859]);
	  set_offsets(cells[631],2,2,0,2);
	  cells[632]->set_vertices(vertices[92],vertices[84],vertices[132],vertices[85]);
	  cells[632]->set_neighbors(cells[579],cells[871],cells[129],cells[800]);
	  set_offsets(cells[632],0,0,0,0);
	  cells[633]->set_vertices(vertices[75],vertices[84],vertices[132],vertices[83]);
	  cells[633]->set_neighbors(cells[800],cells[801],cells[537],cells[89]);
	  set_offsets(cells[633],0,0,0,0);
	  cells[634]->set_vertices(vertices[178],vertices[121],vertices[169],vertices[170]);
	  cells[634]->set_neighbors(cells[444],cells[959],cells[408],cells[1092]);
	  set_offsets(cells[634],0,0,0,0);
	  cells[635]->set_vertices(vertices[99],vertices[51],vertices[59],vertices[106]);
	  cells[635]->set_neighbors(cells[681],cells[649],cells[511],cells[210]);
	  set_offsets(cells[635],0,0,0,0);
	  cells[636]->set_vertices(vertices[97],vertices[49],vertices[104],vertices[96]);
	  cells[636]->set_neighbors(cells[637],cells[674],cells[861],cells[588]);
	  set_offsets(cells[636],0,0,0,0);
	  cells[637]->set_vertices(vertices[96],vertices[49],vertices[104],vertices[56]);
	  cells[637]->set_neighbors(cells[395],cells[150],cells[604],cells[636]);
	  set_offsets(cells[637],0,0,0,0);
	  cells[638]->set_vertices(vertices[103],vertices[104],vertices[96],vertices[144]);
	  cells[638]->set_neighbors(cells[674],cells[489],cells[908],cells[132]);
	  set_offsets(cells[638],0,1,1,1);
	  cells[639]->set_vertices(vertices[61],vertices[109],vertices[110],vertices[101]);
	  cells[639]->set_neighbors(cells[398],cells[281],cells[706],cells[284]);
	  set_offsets(cells[639],0,0,0,0);
	  cells[640]->set_vertices(vertices[69],vertices[126],vertices[78],vertices[118]);
	  cells[640]->set_neighbors(cells[362],cells[292],cells[603],cells[561]);
	  set_offsets(cells[640],0,0,0,0);
	  cells[641]->set_vertices(vertices[198],vertices[150],vertices[151],vertices[158]);
	  cells[641]->set_neighbors(cells[670],cells[146],cells[643],cells[1215]);
	  set_offsets(cells[641],0,0,0,0);
	  cells[642]->set_vertices(vertices[117],vertices[164],vertices[109],vertices[116]);
	  cells[642]->set_neighbors(cells[308],cells[752],cells[693],cells[191]);
	  set_offsets(cells[642],0,0,0,0);
	  cells[643]->set_vertices(vertices[149],vertices[150],vertices[198],vertices[158]);
	  cells[643]->set_neighbors(cells[641],cells[151],cells[942],cells[1216]);
	  set_offsets(cells[643],0,0,0,0);
	  cells[644]->set_vertices(vertices[122],vertices[74],vertices[67],vertices[75]);
	  cells[644]->set_neighbors(cells[515],cells[379],cells[697],cells[397]);
	  set_offsets(cells[644],0,0,0,0);
	  cells[645]->set_vertices(vertices[115],vertices[162],vertices[155],vertices[107]);
	  cells[645]->set_neighbors(cells[770],cells[715],cells[790],cells[1001]);
	  set_offsets(cells[645],0,0,0,0);
	  cells[646]->set_vertices(vertices[129],vertices[177],vertices[178],vertices[169]);
	  cells[646]->set_neighbors(cells[47],cells[1092],cells[1043],cells[946]);
	  set_offsets(cells[646],0,0,0,0);
	  cells[647]->set_vertices(vertices[164],vertices[156],vertices[157],vertices[109]);
	  cells[647]->set_neighbors(cells[971],cells[191],cells[308],cells[1014]);
	  set_offsets(cells[647],0,0,0,0);
	  cells[648]->set_vertices(vertices[148],vertices[99],vertices[139],vertices[100]);
	  cells[648]->set_neighbors(cells[791],cells[1135],cells[894],cells[786]);
	  set_offsets(cells[648],2,2,0,2);
	  cells[649]->set_vertices(vertices[107],vertices[106],vertices[99],vertices[59]);
	  cells[649]->set_neighbors(cells[635],cells[630],cells[282],cells[541]);
	  set_offsets(cells[649],0,0,0,0);
	  cells[650]->set_vertices(vertices[52],vertices[60],vertices[53],vertices[100]);
	  cells[650]->set_neighbors(cells[438],cells[651],cells[372],cells[149]);
	  set_offsets(cells[650],0,0,0,0);
	  cells[651]->set_vertices(vertices[53],vertices[100],vertices[93],vertices[52]);
	  cells[651]->set_neighbors(cells[788],cells[599],cells[650],cells[886]);
	  set_offsets(cells[651],2,2,0,2);
	  cells[652]->set_vertices(vertices[186],vertices[137],vertices[129],vertices[138]);
	  cells[652]->set_neighbors(cells[161],cells[656],cells[688],cells[878]);
	  set_offsets(cells[652],0,0,0,0);
	  cells[653]->set_vertices(vertices[186],vertices[139],vertices[131],vertices[179]);
	  cells[653]->set_neighbors(cells[1061],cells[817],cells[1098],cells[1126]);
	  set_offsets(cells[653],0,0,0,0);
	  cells[654]->set_vertices(vertices[127],vertices[120],vertices[79],vertices[119]);
	  cells[654]->set_neighbors(cells[728],cells[759],cells[1019],cells[88]);
	  set_offsets(cells[654],0,1,0,0);
	  cells[655]->set_vertices(vertices[69],vertices[77],vertices[126],vertices[117]);
	  cells[655]->set_neighbors(cells[400],cells[603],cells[721],cells[561]);
	  set_offsets(cells[655],0,0,0,0);
	  cells[656]->set_vertices(vertices[178],vertices[186],vertices[129],vertices[138]);
	  cells[656]->set_neighbors(cells[652],cells[1017],cells[1127],cells[946]);
	  set_offsets(cells[656],0,0,0,0);
	  cells[657]->set_vertices(vertices[102],vertices[55],vertices[103],vertices[110]);
	  cells[657]->set_neighbors(cells[699],cells[932],cells[23],cells[903]);
	  set_offsets(cells[657],0,0,0,0);
	  cells[658]->set_vertices(vertices[60],vertices[108],vertices[100],vertices[51]);
	  cells[658]->set_neighbors(cells[64],cells[372],cells[7],cells[438]);
	  set_offsets(cells[658],0,0,0,0);
	  cells[659]->set_vertices(vertices[128],vertices[73],vertices[120],vertices[121]);
	  cells[659]->set_neighbors(cells[766],cells[733],cells[394],cells[735]);
	  set_offsets(cells[659],0,0,0,0);
	  cells[660]->set_vertices(vertices[180],vertices[123],vertices[172],vertices[132]);
	  cells[660]->set_neighbors(cells[572],cells[711],cells[1069],cells[1063]);
	  set_offsets(cells[660],0,0,0,0);
	  cells[661]->set_vertices(vertices[79],vertices[80],vertices[128],vertices[87]);
	  cells[661]->set_neighbors(cells[745],cells[737],cells[577],cells[813]);
	  set_offsets(cells[661],0,1,1,0);
	  cells[662]->set_vertices(vertices[120],vertices[113],vertices[160],vertices[168]);
	  cells[662]->set_neighbors(cells[744],cells[987],cells[991],cells[949]);
	  set_offsets(cells[662],0,0,0,0);
	  cells[663]->set_vertices(vertices[67],vertices[124],vertices[116],vertices[115]);
	  cells[663]->set_neighbors(cells[970],cells[286],cells[494],cells[619]);
	  set_offsets(cells[663],0,0,0,0);
	  cells[664]->set_vertices(vertices[276],vertices[267],vertices[28],vertices[268]);
	  cells[664]->set_neighbors(cells[1526],cells[1631],cells[1624],cells[1045]);
	  set_offsets(cells[664],0,0,4,0);
	  cells[665]->set_vertices(vertices[195],vertices[202],vertices[203],vertices[250]);
	  cells[665]->set_neighbors(cells[947],cells[1242],cells[1513],cells[1237]);
	  set_offsets(cells[665],0,0,0,0);
	  cells[666]->set_vertices(vertices[124],vertices[76],vertices[69],vertices[77]);
	  cells[666]->set_neighbors(cells[528],cells[721],cells[545],cells[375]);
	  set_offsets(cells[666],0,0,0,0);
	  cells[667]->set_vertices(vertices[93],vertices[94],vertices[85],vertices[142]);
	  cells[667]->set_neighbors(cells[880],cells[834],cells[881],cells[532]);
	  set_offsets(cells[667],0,0,0,0);
	  cells[668]->set_vertices(vertices[53],vertices[54],vertices[93],vertices[102]);
	  cells[668]->set_neighbors(cells[798],cells[885],cells[411],cells[600]);
	  set_offsets(cells[668],2,2,0,2);
	  cells[669]->set_vertices(vertices[102],vertices[54],vertices[95],vertices[55]);
	  cells[669]->set_neighbors(cells[611],cells[900],cells[371],cells[675]);
	  set_offsets(cells[669],2,2,0,2);
	  cells[670]->set_vertices(vertices[150],vertices[103],vertices[151],vertices[158]);
	  cells[670]->set_neighbors(cells[409],cells[641],cells[253],cells[1157]);
	  set_offsets(cells[670],0,0,0,0);
	  cells[671]->set_vertices(vertices[96],vertices[95],vertices[143],vertices[136]);
	  cells[671]->set_neighbors(cells[672],cells[316],cells[839],cells[899]);
	  set_offsets(cells[671],3,0,0,1);
	  cells[672]->set_vertices(vertices[143],vertices[95],vertices[135],vertices[136]);
	  cells[672]->set_neighbors(cells[531],cells[816],cells[671],cells[902]);
	  set_offsets(cells[672],0,0,0,1);
	  cells[673]->set_vertices(vertices[230],vertices[221],vertices[270],vertices[222]);
	  cells[673]->set_neighbors(cells[1586],cells[1276],cells[1116],cells[1303]);
	  set_offsets(cells[673],0,0,0,0);
	  cells[674]->set_vertices(vertices[144],vertices[104],vertices[96],vertices[97]);
	  cells[674]->set_neighbors(cells[636],cells[67],cells[797],cells[638]);
	  set_offsets(cells[674],0,0,0,0);
	  cells[675]->set_vertices(vertices[102],vertices[54],vertices[142],vertices[95]);
	  cells[675]->set_neighbors(cells[351],cells[901],cells[669],cells[798]);
	  set_offsets(cells[675],2,2,0,0);
	  cells[676]->set_vertices(vertices[53],vertices[62],vertices[110],vertices[61]);
	  cells[676]->set_neighbors(cells[677],cells[281],cells[296],cells[627]);
	  set_offsets(cells[676],0,0,0,0);
	  cells[677]->set_vertices(vertices[61],vertices[70],vertices[62],vertices[110]);
	  cells[677]->set_neighbors(cells[373],cells[676],cells[204],cells[383]);
	  set_offsets(cells[677],0,0,0,0);
	  cells[678]->set_vertices(vertices[112],vertices[105],vertices[152],vertices[160]);
	  cells[678]->set_neighbors(cells[43],cells[436],cells[956],cells[842]);
	  set_offsets(cells[678],0,0,0,0);
	  cells[679]->set_vertices(vertices[57],vertices[66],vertices[106],vertices[114]);
	  cells[679]->set_neighbors(cells[429],cells[212],cells[257],cells[17]);
	  set_offsets(cells[679],0,0,0,0);
	  cells[680]->set_vertices(vertices[152],vertices[112],vertices[111],vertices[104]);
	  cells[680]->set_neighbors(cells[475],cells[867],cells[842],cells[436]);
	  set_offsets(cells[680],1,1,0,1);
	  cells[681]->set_vertices(vertices[106],vertices[58],vertices[51],vertices[59]);
	  cells[681]->set_neighbors(cells[424],cells[635],cells[96],cells[312]);
	  set_offsets(cells[681],0,0,0,0);
	  cells[682]->set_vertices(vertices[71],vertices[112],vertices[63],vertices[111]);
	  cells[682]->set_neighbors(cells[475],cells[430],cells[533],cells[725]);
	  set_offsets(cells[682],0,1,0,0);
	  cells[683]->set_vertices(vertices[105],vertices[57],vertices[104],vertices[97]);
	  cells[683]->set_neighbors(cells[588],cells[593],cells[304],cells[227]);
	  set_offsets(cells[683],0,0,0,0);
	  cells[684]->set_vertices(vertices[79],vertices[80],vertices[72],vertices[120]);
	  cells[684]->set_neighbors(cells[228],cells[732],cells[813],cells[38]);
	  set_offsets(cells[684],0,1,1,1);
	  cells[685]->set_vertices(vertices[137],vertices[96],vertices[184],vertices[136]);
	  cells[685]->set_neighbors(cells[316],cells[445],cells[857],cells[1072]);
	  set_offsets(cells[685],0,2,0,0);
	  cells[686]->set_vertices(vertices[192],vertices[144],vertices[145],vertices[152]);
	  cells[686]->set_neighbors(cells[614],cells[1079],cells[893],cells[914]);
	  set_offsets(cells[686],0,0,0,0);
	  cells[687]->set_vertices(vertices[123],vertices[170],vertices[163],vertices[115]);
	  cells[687]->set_neighbors(cells[804],cells[826],cells[924],cells[1055]);
	  set_offsets(cells[687],0,0,0,0);
	  cells[688]->set_vertices(vertices[186],vertices[137],vertices[138],vertices[98]);
	  cells[688]->set_neighbors(cells[856],cells[413],cells[1124],cells[652]);
	  set_offsets(cells[688],0,0,0,2);
	  cells[689]->set_vertices(vertices[176],vertices[128],vertices[127],vertices[168]);
	  cells[689]->set_neighbors(cells[1038],cells[977],cells[990],cells[1035]);
	  set_offsets(cells[689],1,1,0,1);
	  cells[690]->set_vertices(vertices[166],vertices[158],vertices[111],vertices[118]);
	  cells[690]->set_neighbors(cells[354],cells[615],cells[29],cells[506]);
	  set_offsets(cells[690],0,0,0,0);
	  cells[691]->set_vertices(vertices[101],vertices[158],vertices[149],vertices[109]);
	  cells[691]->set_neighbors(cells[164],cells[747],cells[398],cells[942]);
	  set_offsets(cells[691],0,0,0,0);
	  cells[692]->set_vertices(vertices[149],vertices[158],vertices[206],vertices[157]);
	  cells[692]->set_neighbors(cells[1041],cells[558],cells[164],cells[151]);
	  set_offsets(cells[692],0,0,0,0);
	  cells[693]->set_vertices(vertices[124],vertices[164],vertices[117],vertices[116]);
	  cells[693]->set_neighbors(cells[642],cells[410],cells[970],cells[984]);
	  set_offsets(cells[693],0,0,0,0);
	  cells[694]->set_vertices(vertices[190],vertices[133],vertices[181],vertices[182]);
	  cells[694]->set_neighbors(cells[986],cells[1165],cells[935],cells[934]);
	  set_offsets(cells[694],0,0,0,0);
	  cells[695]->set_vertices(vertices[133],vertices[188],vertices[140],vertices[141]);
	  cells[695]->set_neighbors(cells[1113],cells[884],cells[825],cells[1123]);
	  set_offsets(cells[695],0,0,0,0);
	  cells[696]->set_vertices(vertices[91],vertices[98],vertices[138],vertices[50]);
	  cells[696]->set_neighbors(cells[855],cells[14],cells[187],cells[147]);
	  set_offsets(cells[696],0,2,0,2);
	  cells[697]->set_vertices(vertices[122],vertices[74],vertices[75],vertices[82]);
	  cells[697]->set_neighbors(cells[40],cells[814],cells[491],cells[644]);
	  set_offsets(cells[697],0,0,0,0);
	  cells[698]->set_vertices(vertices[111],vertices[104],vertices[63],vertices[103]);
	  cells[698]->set_neighbors(cells[330],cells[717],cells[867],cells[475]);
	  set_offsets(cells[698],0,1,0,0);
	  cells[699]->set_vertices(vertices[110],vertices[55],vertices[103],vertices[63]);
	  cells[699]->set_neighbors(cells[330],cells[717],cells[487],cells[657]);
	  set_offsets(cells[699],0,0,0,0);
	  cells[700]->set_vertices(vertices[116],vertices[108],vertices[109],vertices[61]);
	  cells[700]->set_neighbors(cells[706],cells[317],cells[193],cells[751]);
	  set_offsets(cells[700],0,0,0,0);
	  cells[701]->set_vertices(vertices[161],vertices[113],vertices[160],vertices[153]);
	  cells[701]->set_neighbors(cells[955],cells[1217],cells[746],cells[744]);
	  set_offsets(cells[701],0,0,0,0);
	  cells[702]->set_vertices(vertices[111],vertices[158],vertices[103],vertices[110]);
	  cells[702]->set_neighbors(cells[253],cells[717],cells[354],cells[409]);
	  set_offsets(cells[702],0,0,0,0);
	  cells[703]->set_vertices(vertices[163],vertices[115],vertices[155],vertices[164]);
	  cells[703]->set_neighbors(cells[715],cells[57],cells[772],cells[1001]);
	  set_offsets(cells[703],0,0,0,0);
	  cells[704]->set_vertices(vertices[147],vertices[99],vertices[107],vertices[154]);
	  cells[704]->set_neighbors(cells[541],cells[961],cells[461],cells[153]);
	  set_offsets(cells[704],0,0,0,0);
	  cells[705]->set_vertices(vertices[109],vertices[156],vertices[101],vertices[108]);
	  cells[705]->set_neighbors(cells[471],cells[706],cells[751],cells[747]);
	  set_offsets(cells[705],0,0,0,0);
	  cells[706]->set_vertices(vertices[109],vertices[108],vertices[101],vertices[61]);
	  cells[706]->set_neighbors(cells[51],cells[639],cells[700],cells[705]);
	  set_offsets(cells[706],0,0,0,0);
	  cells[707]->set_vertices(vertices[101],vertices[148],vertices[141],vertices[100]);
	  cells[707]->set_neighbors(cells[332],cells[887],cells[734],cells[390]);
	  set_offsets(cells[707],2,2,0,2);
	  cells[708]->set_vertices(vertices[255],vertices[14],vertices[7],vertices[247]);
	  cells[708]->set_neighbors(cells[266],cells[1263],cells[944],cells[1495]);
	  set_offsets(cells[708],0,4,4,0);
	  cells[709]->set_vertices(vertices[231],vertices[230],vertices[278],vertices[223]);
	  cells[709]->set_neighbors(cells[1636],cells[1643],cells[1297],cells[1675]);
	  set_offsets(cells[709],0,0,0,0);
	  cells[710]->set_vertices(vertices[134],vertices[174],vertices[127],vertices[126]);
	  cells[710]->set_neighbors(cells[196],cells[807],cells[736],cells[780]);
	  set_offsets(cells[710],0,0,0,0);
	  cells[711]->set_vertices(vertices[172],vertices[180],vertices[132],vertices[125]);
	  cells[711]->set_neighbors(cells[1064],cells[928],cells[993],cells[660]);
	  set_offsets(cells[711],0,0,0,0);
	  cells[712]->set_vertices(vertices[172],vertices[164],vertices[165],vertices[117]);
	  cells[712]->set_neighbors(cells[1018],cells[91],cells[984],cells[983]);
	  set_offsets(cells[712],0,0,0,0);
	  cells[713]->set_vertices(vertices[129],vertices[121],vertices[176],vertices[169]);
	  cells[713]->set_neighbors(cells[571],cells[1043],cells[1092],cells[516]);
	  set_offsets(cells[713],0,0,0,0);
	  cells[714]->set_vertices(vertices[229],vertices[236],vertices[276],vertices[228]);
	  cells[714]->set_neighbors(cells[1362],cells[1664],cells[1184],cells[1665]);
	  set_offsets(cells[714],0,0,0,0);
	  cells[715]->set_vertices(vertices[115],vertices[107],vertices[155],vertices[164]);
	  cells[715]->set_neighbors(cells[866],cells[703],cells[802],cells[645]);
	  set_offsets(cells[715],0,0,0,0);
	  cells[716]->set_vertices(vertices[119],vertices[118],vertices[111],vertices[71]);
	  cells[716]->set_neighbors(cells[430],cells[533],cells[755],cells[615]);
	  set_offsets(cells[716],0,0,0,0);
	  cells[717]->set_vertices(vertices[111],vertices[110],vertices[103],vertices[63]);
	  cells[717]->set_neighbors(cells[699],cells[698],cells[225],cells[702]);
	  set_offsets(cells[717],0,0,0,0);
	  cells[718]->set_vertices(vertices[192],vertices[145],vertices[193],vertices[200]);
	  cells[718]->set_neighbors(cells[1142],cells[1143],cells[1079],cells[1361]);
	  set_offsets(cells[718],0,0,0,0);
	  cells[719]->set_vertices(vertices[194],vertices[146],vertices[154],vertices[145]);
	  cells[719]->set_neighbors(cells[457],cells[1150],cells[602],cells[922]);
	  set_offsets(cells[719],0,0,0,0);
	  cells[720]->set_vertices(vertices[99],vertices[108],vertices[148],vertices[156]);
	  cells[720]->set_neighbors(cells[471],cells[130],cells[907],cells[894]);
	  set_offsets(cells[720],0,0,0,0);
	  cells[721]->set_vertices(vertices[77],vertices[124],vertices[117],vertices[69]);
	  cells[721]->set_neighbors(cells[410],cells[655],cells[666],cells[796]);
	  set_offsets(cells[721],0,0,0,0);
	  cells[722]->set_vertices(vertices[61],vertices[69],vertices[70],vertices[118]);
	  cells[722]->set_neighbors(cells[292],cells[204],cells[547],cells[201]);
	  set_offsets(cells[722],0,0,0,0);
	  cells[723]->set_vertices(vertices[118],vertices[78],vertices[70],vertices[71]);
	  cells[723]->set_neighbors(cells[355],cells[294],cells[362],cells[292]);
	  set_offsets(cells[723],0,0,0,0);
	  cells[724]->set_vertices(vertices[71],vertices[72],vertices[64],vertices[112]);
	  cells[724]->set_neighbors(cells[180],cells[725],cells[730],cells[203]);
	  set_offsets(cells[724],0,1,1,1);
	  cells[725]->set_vertices(vertices[63],vertices[112],vertices[71],vertices[64]);
	  cells[725]->set_neighbors(cells[724],cells[169],cells[726],cells[682]);
	  set_offsets(cells[725],0,1,0,1);
	  cells[726]->set_vertices(vertices[63],vertices[104],vertices[112],vertices[64]);
	  cells[726]->set_neighbors(cells[509],cells[725],cells[209],cells[475]);
	  set_offsets(cells[726],0,1,1,1);
	  cells[727]->set_vertices(vertices[95],vertices[88],vertices[136],vertices[48]);
	  cells[727]->set_neighbors(cells[787],cells[839],cells[608],cells[806]);
	  set_offsets(cells[727],0,1,1,3);
	  cells[728]->set_vertices(vertices[79],vertices[120],vertices[71],vertices[119]);
	  cells[728]->set_neighbors(cells[731],cells[446],cells[654],cells[732]);
	  set_offsets(cells[728],0,1,0,0);
	  cells[729]->set_vertices(vertices[113],vertices[65],vertices[112],vertices[105]);
	  cells[729]->set_neighbors(cells[510],cells[956],cells[504],cells[272]);
	  set_offsets(cells[729],0,0,0,0);
	  cells[730]->set_vertices(vertices[112],vertices[72],vertices[120],vertices[71]);
	  cells[730]->set_neighbors(cells[732],cells[731],cells[724],cells[456]);
	  set_offsets(cells[730],1,1,1,0);
	  cells[731]->set_vertices(vertices[120],vertices[112],vertices[71],vertices[119]);
	  cells[731]->set_neighbors(cells[533],cells[728],cells[978],cells[730]);
	  set_offsets(cells[731],1,1,0,0);
	  cells[732]->set_vertices(vertices[120],vertices[72],vertices[79],vertices[71]);
	  cells[732]->set_neighbors(cells[401],cells[728],cells[730],cells[684]);
	  set_offsets(cells[732],1,1,0,0);
	  cells[733]->set_vertices(vertices[121],vertices[128],vertices[168],vertices[120]);
	  cells[733]->set_neighbors(cells[1038],cells[991],cells[659],cells[990]);
	  set_offsets(cells[733],0,0,0,0);
	  cells[734]->set_vertices(vertices[148],vertices[100],vertices[101],vertices[108]);
	  cells[734]->set_neighbors(cells[382],cells[471],cells[894],cells[707]);
	  set_offsets(cells[734],0,0,0,0);
	  cells[735]->set_vertices(vertices[120],vertices[80],vertices[73],vertices[128]);
	  cells[735]->set_neighbors(cells[774],cells[659],cells[813],cells[228]);
	  set_offsets(cells[735],0,0,0,0);
	  cells[736]->set_vertices(vertices[134],vertices[125],vertices[174],vertices[126]);
	  cells[736]->set_neighbors(cells[933],cells[710],cells[799],cells[1070]);
	  set_offsets(cells[736],0,0,0,0);
	  cells[737]->set_vertices(vertices[87],vertices[128],vertices[79],vertices[127]);
	  cells[737]->set_neighbors(cells[88],cells[808],cells[758],cells[661]);
	  set_offsets(cells[737],0,1,0,0);
	  cells[738]->set_vertices(vertices[65],vertices[72],vertices[73],vertices[120]);
	  cells[738]->set_neighbors(cells[228],cells[19],cells[456],cells[505]);
	  set_offsets(cells[738],0,0,0,0);
	  cells[739]->set_vertices(vertices[118],vertices[109],vertices[158],vertices[110]);
	  cells[739]->set_neighbors(cells[398],cells[354],cells[284],cells[29]);
	  set_offsets(cells[739],0,0,0,0);
	  cells[740]->set_vertices(vertices[101],vertices[110],vertices[150],vertices[158]);
	  cells[740]->set_neighbors(cells[253],cells[942],cells[398],cells[831]);
	  set_offsets(cells[740],0,0,0,0);
	  cells[741]->set_vertices(vertices[168],vertices[160],vertices[119],vertices[167]);
	  cells[741]->set_neighbors(cells[1036],cells[809],cells[1269],cells[987]);
	  set_offsets(cells[741],1,1,0,0);
	  cells[742]->set_vertices(vertices[115],vertices[114],vertices[107],vertices[67]);
	  cells[742]->set_neighbors(cells[469],cells[286],cells[20],cells[790]);
	  set_offsets(cells[742],0,0,0,0);
	  cells[743]->set_vertices(vertices[132],vertices[140],vertices[131],vertices[83]);
	  cells[743]->set_neighbors(cells[835],cells[763],cells[369],cells[1068]);
	  set_offsets(cells[743],0,0,0,0);
	  cells[744]->set_vertices(vertices[168],vertices[113],vertices[160],vertices[161]);
	  cells[744]->set_neighbors(cells[701],cells[1087],cells[992],cells[662]);
	  set_offsets(cells[744],0,0,0,0);
	  cells[745]->set_vertices(vertices[87],vertices[80],vertices[128],vertices[88]);
	  cells[745]->set_neighbors(cells[773],cells[273],cells[546],cells[661]);
	  set_offsets(cells[745],0,1,1,1);
	  cells[746]->set_vertices(vertices[161],vertices[113],vertices[153],vertices[162]);
	  cells[746]->set_neighbors(cells[512],cells[158],cells[954],cells[701]);
	  set_offsets(cells[746],0,0,0,0);
	  cells[747]->set_vertices(vertices[156],vertices[101],vertices[149],vertices[109]);
	  cells[747]->set_neighbors(cells[691],cells[971],cells[705],cells[417]);
	  set_offsets(cells[747],0,0,0,0);
	  cells[748]->set_vertices(vertices[122],vertices[113],vertices[162],vertices[114]);
	  cells[748]->set_neighbors(cells[415],cells[517],cells[231],cells[818]);
	  set_offsets(cells[748],0,0,0,0);
	  cells[749]->set_vertices(vertices[251],vertices[258],vertices[18],vertices[10]);
	  cells[749]->set_neighbors(cells[44],cells[912],cells[1277],cells[421]);
	  set_offsets(cells[749],0,0,4,4);
	  cells[750]->set_vertices(vertices[119],vertices[166],vertices[159],vertices[111]);
	  cells[750]->set_neighbors(cells[506],cells[360],cells[615],cells[1037]);
	  set_offsets(cells[750],0,0,0,0);
	  cells[751]->set_vertices(vertices[116],vertices[156],vertices[109],vertices[108]);
	  cells[751]->set_neighbors(cells[705],cells[700],cells[756],cells[308]);
	  set_offsets(cells[751],0,0,0,0);
	  cells[752]->set_vertices(vertices[117],vertices[116],vertices[109],vertices[69]);
	  cells[752]->set_neighbors(cells[317],cells[580],cells[410],cells[642]);
	  set_offsets(cells[752],0,0,0,0);
	  cells[753]->set_vertices(vertices[183],vertices[224],vertices[175],vertices[223]);
	  cells[753]->set_neighbors(cells[1085],cells[851],cells[981],cells[1359]);
	  set_offsets(cells[753],0,1,0,0);
	  cells[754]->set_vertices(vertices[154],vertices[146],vertices[99],vertices[106]);
	  cells[754]->set_neighbors(cells[595],cells[541],cells[333],cells[461]);
	  set_offsets(cells[754],0,0,0,0);
	  cells[755]->set_vertices(vertices[126],vertices[118],vertices[119],vertices[71]);
	  cells[755]->set_neighbors(cells[716],cells[446],cells[362],cells[1013]);
	  set_offsets(cells[755],0,0,0,0);
	  cells[756]->set_vertices(vertices[116],vertices[107],vertices[156],vertices[108]);
	  cells[756]->set_neighbors(cells[907],cells[751],cells[478],cells[761]);
	  set_offsets(cells[756],0,0,0,0);
	  cells[757]->set_vertices(vertices[150],vertices[102],vertices[190],vertices[143]);
	  cells[757]->set_neighbors(cells[1118],cells[1156],cells[938],cells[895]);
	  set_offsets(cells[757],2,2,0,0);
	  cells[758]->set_vertices(vertices[128],vertices[135],vertices[127],vertices[87]);
	  cells[758]->set_neighbors(cells[822],cells[737],cells[174],cells[1035]);
	  set_offsets(cells[758],1,0,0,0);
	  cells[759]->set_vertices(vertices[127],vertices[126],vertices[119],vertices[79]);
	  cells[759]->set_neighbors(cells[446],cells[654],cells[807],cells[196]);
	  set_offsets(cells[759],0,0,0,0);
	  cells[760]->set_vertices(vertices[183],vertices[190],vertices[182],vertices[135]);
	  cells[760]->set_neighbors(cells[1083],cells[996],cells[920],cells[1164]);
	  set_offsets(cells[760],0,0,0,0);
	  cells[761]->set_vertices(vertices[164],vertices[107],vertices[156],vertices[116]);
	  cells[761]->set_neighbors(cells[756],cells[308],cells[802],cells[866]);
	  set_offsets(cells[761],0,0,0,0);
	  cells[762]->set_vertices(vertices[112],vertices[160],vertices[111],vertices[119]);
	  cells[762]->set_neighbors(cells[360],cells[533],cells[978],cells[436]);
	  set_offsets(cells[762],1,1,0,0);
	  cells[763]->set_vertices(vertices[131],vertices[132],vertices[83],vertices[123]);
	  cells[763]->set_neighbors(cells[801],cells[819],cells[1069],cells[743]);
	  set_offsets(cells[763],0,0,0,0);
	  cells[764]->set_vertices(vertices[73],vertices[130],vertices[121],vertices[81]);
	  cells[764]->set_neighbors(cells[98],cells[394],cells[601],cells[115]);
	  set_offsets(cells[764],0,0,0,0);
	  cells[765]->set_vertices(vertices[136],vertices[88],vertices[81],vertices[89]);
	  cells[765]->set_neighbors(cells[378],cells[299],cells[787],cells[425]);
	  set_offsets(cells[765],0,0,0,0);
	  cells[766]->set_vertices(vertices[121],vertices[73],vertices[120],vertices[113]);
	  cells[766]->set_neighbors(cells[19],cells[991],cells[94],cells[659]);
	  set_offsets(cells[766],0,0,0,0);
	  cells[767]->set_vertices(vertices[160],vertices[152],vertices[111],vertices[159]);
	  cells[767]->set_neighbors(cells[241],cells[360],cells[1214],cells[436]);
	  set_offsets(cells[767],1,1,0,0);
	  cells[768]->set_vertices(vertices[187],vertices[148],vertices[236],vertices[196]);
	  cells[768]->set_neighbors(cells[126],cells[1196],cells[1133],cells[1327]);
	  set_offsets(cells[768],0,2,0,2);
	  cells[769]->set_vertices(vertices[247],vertices[246],vertices[199],vertices[254]);
	  cells[769]->set_neighbors(cells[1300],cells[173],cells[1485],cells[1722]);
	  set_offsets(cells[769],0,0,0,0);
	  cells[770]->set_vertices(vertices[162],vertices[154],vertices[155],vertices[107]);
	  cells[770]->set_neighbors(cells[961],cells[645],cells[622],cells[1179]);
	  set_offsets(cells[770],0,0,0,0);
	  cells[771]->set_vertices(vertices[178],vertices[186],vertices[179],vertices[226]);
	  cells[771]->set_neighbors(cells[1384],cells[1332],cells[1039],cells[817]);
	  set_offsets(cells[771],0,0,0,0);
	  cells[772]->set_vertices(vertices[172],vertices[115],vertices[163],vertices[164]);
	  cells[772]->set_neighbors(cells[703],cells[833],cells[56],cells[826]);
	  set_offsets(cells[772],0,0,0,0);
	  cells[773]->set_vertices(vertices[128],vertices[80],vertices[81],vertices[88]);
	  cells[773]->set_neighbors(cells[534],cells[425],cells[745],cells[774]);
	  set_offsets(cells[773],0,0,0,0);
	  cells[774]->set_vertices(vertices[128],vertices[80],vertices[73],vertices[81]);
	  cells[774]->set_neighbors(cells[556],cells[394],cells[773],cells[735]);
	  set_offsets(cells[774],0,0,0,0);
	  cells[775]->set_vertices(vertices[170],vertices[210],vertices[163],vertices[162]);
	  cells[775]->set_neighbors(cells[1010],cells[804],cells[1008],cells[1289]);
	  set_offsets(cells[775],0,0,0,0);
	  cells[776]->set_vertices(vertices[83],vertices[130],vertices[90],vertices[138]);
	  cells[776]->set_neighbors(cells[591],cells[820],cells[829],cells[388]);
	  set_offsets(cells[776],0,0,0,0);
	  cells[777]->set_vertices(vertices[77],vertices[134],vertices[86],vertices[126]);
	  cells[777]->set_neighbors(cells[414],cells[442],cells[799],cells[311]);
	  set_offsets(cells[777],0,0,0,0);
	  cells[778]->set_vertices(vertices[73],vertices[82],vertices[122],vertices[130]);
	  cells[778]->set_neighbors(cells[814],cells[115],cells[601],cells[491]);
	  set_offsets(cells[778],0,0,0,0);
	  cells[779]->set_vertices(vertices[270],vertices[262],vertices[263],vertices[215]);
	  cells[779]->set_neighbors(cells[1308],cells[1296],cells[1591],cells[264]);
	  set_offsets(cells[779],0,0,0,0);
	  cells[780]->set_vertices(vertices[127],vertices[182],vertices[174],vertices[134]);
	  cells[780]->set_neighbors(cells[1070],cells[710],cells[159],cells[422]);
	  set_offsets(cells[780],0,0,0,0);
	  cells[781]->set_vertices(vertices[83],vertices[130],vertices[123],vertices[75]);
	  cells[781]->set_neighbors(cells[594],cells[801],cells[68],cells[819]);
	  set_offsets(cells[781],0,0,0,0);
	  cells[782]->set_vertices(vertices[176],vertices[136],vertices[135],vertices[128]);
	  cells[782]->set_neighbors(cells[174],cells[1035],cells[1073],cells[1074]);
	  set_offsets(cells[782],1,1,0,1);
	  cells[783]->set_vertices(vertices[162],vertices[153],vertices[202],vertices[154]);
	  cells[783]->set_neighbors(cells[1151],cells[1179],cells[969],cells[208]);
	  set_offsets(cells[783],0,0,0,0);
	  cells[784]->set_vertices(vertices[211],vertices[251],vertices[259],vertices[260]);
	  cells[784]->set_neighbors(cells[1280],cells[1331],cells[1252],cells[1520]);
	  set_offsets(cells[784],0,0,0,0);
	  cells[785]->set_vertices(vertices[174],vertices[166],vertices[119],vertices[126]);
	  cells[785]->set_neighbors(cells[1013],cells[196],cells[463],cells[24]);
	  set_offsets(cells[785],0,0,0,0);
	  cells[786]->set_vertices(vertices[187],vertices[99],vertices[139],vertices[148]);
	  cells[786]->set_neighbors(cells[648],cells[1134],cells[1099],cells[1137]);
	  set_offsets(cells[786],0,2,0,2);
	  cells[787]->set_vertices(vertices[136],vertices[88],vertices[89],vertices[48]);
	  cells[787]->set_neighbors(cells[70],cells[350],cells[727],cells[765]);
	  set_offsets(cells[787],0,0,0,2);
	  cells[788]->set_vertices(vertices[93],vertices[100],vertices[140],vertices[52]);
	  cells[788]->set_neighbors(cells[492],cells[50],cells[651],cells[883]);
	  set_offsets(cells[788],0,2,0,2);
	  cells[789]->set_vertices(vertices[138],vertices[89],vertices[81],vertices[90]);
	  cells[789]->set_neighbors(cells[198],cells[591],cells[821],cells[381]);
	  set_offsets(cells[789],0,0,0,0);
	  cells[790]->set_vertices(vertices[115],vertices[162],vertices[107],vertices[114]);
	  cells[790]->set_neighbors(cells[622],cells[742],cells[517],cells[645]);
	  set_offsets(cells[790],0,0,0,0);
	  cells[791]->set_vertices(vertices[99],vertices[51],vertices[139],vertices[100]);
	  cells[791]->set_neighbors(cells[823],cells[648],cells[64],cells[872]);
	  set_offsets(cells[791],2,2,0,2);
	  cells[792]->set_vertices(vertices[121],vertices[113],vertices[170],vertices[122]);
	  cells[792]->set_neighbors(cells[818],cells[521],cells[94],cells[139]);
	  set_offsets(cells[792],0,0,0,0);
	  cells[793]->set_vertices(vertices[162],vertices[105],vertices[154],vertices[114]);
	  cells[793]->set_neighbors(cells[950],cells[622],cells[415],cells[969]);
	  set_offsets(cells[793],0,0,0,0);
	  cells[794]->set_vertices(vertices[90],vertices[82],vertices[81],vertices[130]);
	  cells[794]->set_neighbors(cells[601],cells[591],cells[388],cells[290]);
	  set_offsets(cells[794],0,0,0,0);
	  cells[795]->set_vertices(vertices[85],vertices[132],vertices[125],vertices[77]);
	  cells[795]->set_neighbors(cells[472],cells[462],cells[579],cells[828]);
	  set_offsets(cells[795],0,0,0,0);
	  cells[796]->set_vertices(vertices[125],vertices[124],vertices[117],vertices[77]);
	  cells[796]->set_neighbors(cells[721],cells[400],cells[472],cells[459]);
	  set_offsets(cells[796],0,0,0,0);
	  cells[797]->set_vertices(vertices[104],vertices[97],vertices[144],vertices[152]);
	  cells[797]->set_neighbors(cells[614],cells[908],cells[593],cells[674]);
	  set_offsets(cells[797],0,0,0,0);
	  cells[798]->set_vertices(vertices[93],vertices[54],vertices[142],vertices[102]);
	  cells[798]->set_neighbors(cells[675],cells[879],cells[668],cells[881]);
	  set_offsets(cells[798],0,2,0,2);
	  cells[799]->set_vertices(vertices[77],vertices[134],vertices[126],vertices[125]);
	  cells[799]->set_neighbors(cells[736],cells[400],cells[462],cells[777]);
	  set_offsets(cells[799],0,0,0,0);
	  cells[800]->set_vertices(vertices[92],vertices[84],vertices[83],vertices[132]);
	  cells[800]->set_neighbors(cells[633],cells[369],cells[632],cells[302]);
	  set_offsets(cells[800],0,0,0,0);
	  cells[801]->set_vertices(vertices[75],vertices[132],vertices[123],vertices[83]);
	  cells[801]->set_neighbors(cells[763],cells[781],cells[633],cells[418]);
	  set_offsets(cells[801],0,0,0,0);
	  cells[802]->set_vertices(vertices[115],vertices[107],vertices[164],vertices[116]);
	  cells[802]->set_neighbors(cells[761],cells[970],cells[286],cells[715]);
	  set_offsets(cells[802],0,0,0,0);
	  cells[803]->set_vertices(vertices[3],vertices[12],vertices[11],vertices[251]);
	  cells[803]->set_neighbors(cells[183],cells[9],cells[1468],cells[105]);
	  set_offsets(cells[803],4,4,4,0);
	  cells[804]->set_vertices(vertices[170],vertices[162],vertices[163],vertices[115]);
	  cells[804]->set_neighbors(cells[1001],cells[687],cells[841],cells[775]);
	  set_offsets(cells[804],0,0,0,0);
	  cells[805]->set_vertices(vertices[141],vertices[142],vertices[133],vertices[190]);
	  cells[805]->set_neighbors(cells[935],cells[934],cells[1046],cells[882]);
	  set_offsets(cells[805],0,0,0,0);
	  cells[806]->set_vertices(vertices[87],vertices[88],vertices[136],vertices[95]);
	  cells[806]->set_neighbors(cells[727],cells[531],cells[607],cells[273]);
	  set_offsets(cells[806],0,1,1,0);
	  cells[807]->set_vertices(vertices[134],vertices[126],vertices[127],vertices[79]);
	  cells[807]->set_neighbors(cells[759],cells[808],cells[414],cells[710]);
	  set_offsets(cells[807],0,0,0,0);
	  cells[808]->set_vertices(vertices[87],vertices[134],vertices[127],vertices[79]);
	  cells[808]->set_neighbors(cells[807],cells[737],cells[55],cells[822]);
	  set_offsets(cells[808],0,0,0,0);
	  cells[809]->set_vertices(vertices[127],vertices[168],vertices[119],vertices[167]);
	  cells[809]->set_neighbors(cells[741],cells[869],cells[1078],cells[1019]);
	  set_offsets(cells[809],0,1,0,0);
	  cells[810]->set_vertices(vertices[175],vertices[222],vertices[215],vertices[167]);
	  cells[810]->set_neighbors(cells[1306],cells[69],cells[1310],cells[989]);
	  set_offsets(cells[810],0,0,0,0);
	  cells[811]->set_vertices(vertices[210],vertices[153],vertices[201],vertices[202]);
	  cells[811]->set_neighbors(cells[48],cells[1500],cells[208],cells[1027]);
	  set_offsets(cells[811],0,0,0,0);
	  cells[812]->set_vertices(vertices[85],vertices[140],vertices[92],vertices[93]);
	  cells[812]->set_neighbors(cells[50],cells[598],cells[543],cells[871]);
	  set_offsets(cells[812],0,0,0,0);
	  cells[813]->set_vertices(vertices[79],vertices[80],vertices[120],vertices[128]);
	  cells[813]->set_neighbors(cells[735],cells[88],cells[661],cells[684]);
	  set_offsets(cells[813],0,1,1,1);
	  cells[814]->set_vertices(vertices[130],vertices[82],vertices[122],vertices[75]);
	  cells[814]->set_neighbors(cells[697],cells[594],cells[68],cells[778]);
	  set_offsets(cells[814],0,0,0,0);
	  cells[815]->set_vertices(vertices[137],vertices[89],vertices[136],vertices[129]);
	  cells[815]->set_neighbors(cells[299],cells[445],cells[161],cells[857]);
	  set_offsets(cells[815],0,0,0,0);
	  cells[816]->set_vertices(vertices[184],vertices[143],vertices[135],vertices[136]);
	  cells[816]->set_neighbors(cells[672],cells[1074],cells[316],cells[904]);
	  set_offsets(cells[816],1,0,0,1);
	  cells[817]->set_vertices(vertices[178],vertices[186],vertices[131],vertices[179]);
	  cells[817]->set_neighbors(cells[653],cells[1097],cells[771],cells[1127]);
	  set_offsets(cells[817],0,0,0,0);
	  cells[818]->set_vertices(vertices[170],vertices[113],vertices[162],vertices[122]);
	  cells[818]->set_neighbors(cells[748],cells[841],cells[792],cells[954]);
	  set_offsets(cells[818],0,0,0,0);
	  cells[819]->set_vertices(vertices[130],vertices[131],vertices[83],vertices[123]);
	  cells[819]->set_neighbors(cells[763],cells[781],cells[1050],cells[829]);
	  set_offsets(cells[819],0,0,0,0);
	  cells[820]->set_vertices(vertices[83],vertices[138],vertices[90],vertices[91]);
	  cells[820]->set_neighbors(cells[14],cells[152],cells[340],cells[776]);
	  set_offsets(cells[820],0,0,0,0);
	  cells[821]->set_vertices(vertices[138],vertices[89],vertices[90],vertices[50]);
	  cells[821]->set_neighbors(cells[399],cells[14],cells[855],cells[789]);
	  set_offsets(cells[821],0,0,0,2);
	  cells[822]->set_vertices(vertices[127],vertices[135],vertices[134],vertices[87]);
	  cells[822]->set_neighbors(cells[140],cells[808],cells[758],cells[159]);
	  set_offsets(cells[822],0,0,0,0);
	  cells[823]->set_vertices(vertices[139],vertices[51],vertices[91],vertices[100]);
	  cells[823]->set_neighbors(cells[162],cells[827],cells[791],cells[54]);
	  set_offsets(cells[823],0,2,0,2);
	  cells[824]->set_vertices(vertices[99],vertices[146],vertices[139],vertices[98]);
	  cells[824]->set_neighbors(cells[1081],cells[872],cells[595],cells[1137]);
	  set_offsets(cells[824],2,2,0,2);
	  cells[825]->set_vertices(vertices[188],vertices[141],vertices[133],vertices[181]);
	  cells[825]->set_neighbors(cells[934],cells[490],cells[926],cells[695]);
	  set_offsets(cells[825],0,0,0,0);
	  cells[826]->set_vertices(vertices[123],vertices[115],vertices[163],vertices[172]);
	  cells[826]->set_neighbors(cells[772],cells[917],cells[845],cells[687]);
	  set_offsets(cells[826],0,0,0,0);
	  cells[827]->set_vertices(vertices[100],vertices[91],vertices[139],vertices[140]);
	  cells[827]->set_neighbors(cells[836],cells[1096],cells[492],cells[823]);
	  set_offsets(cells[827],2,0,0,0);
	  cells[828]->set_vertices(vertices[132],vertices[133],vertices[85],vertices[125]);
	  cells[828]->set_neighbors(cells[109],cells[795],cells[1064],cells[625]);
	  set_offsets(cells[828],0,0,0,0);
	  cells[829]->set_vertices(vertices[130],vertices[138],vertices[83],vertices[131]);
	  cells[829]->set_neighbors(cells[340],cells[819],cells[999],cells[776]);
	  set_offsets(cells[829],0,0,0,0);
	  cells[830]->set_vertices(vertices[229],vertices[236],vertices[181],vertices[189]);
	  cells[830]->set_neighbors(cells[1114],cells[1163],cells[1427],cells[1184]);
	  set_offsets(cells[830],0,0,0,0);
	  cells[831]->set_vertices(vertices[101],vertices[102],vertices[150],vertices[110]);
	  cells[831]->set_neighbors(cells[932],cells[740],cells[427],cells[937]);
	  set_offsets(cells[831],0,0,0,0);
	  cells[832]->set_vertices(vertices[220],vertices[211],vertices[260],vertices[212]);
	  cells[832]->set_neighbors(cells[1336],cells[1284],cells[1062],cells[1365]);
	  set_offsets(cells[832],0,0,0,0);
	  cells[833]->set_vertices(vertices[172],vertices[163],vertices[212],vertices[164]);
	  cells[833]->set_neighbors(cells[57],cells[983],cells[772],cells[81]);
	  set_offsets(cells[833],0,0,0,0);
	  cells[834]->set_vertices(vertices[93],vertices[85],vertices[133],vertices[142]);
	  cells[834]->set_neighbors(cells[389],cells[882],cells[667],cells[543]);
	  set_offsets(cells[834],0,0,0,0);
	  cells[835]->set_vertices(vertices[140],vertices[91],vertices[131],vertices[83]);
	  cells[835]->set_neighbors(cells[340],cells[743],cells[416],cells[836]);
	  set_offsets(cells[835],0,0,0,0);
	  cells[836]->set_vertices(vertices[140],vertices[91],vertices[139],vertices[131]);
	  cells[836]->set_neighbors(cells[376],cells[620],cells[835],cells[827]);
	  set_offsets(cells[836],0,0,0,0);
	  cells[837]->set_vertices(vertices[223],vertices[216],vertices[175],vertices[215]);
	  cells[837]->set_neighbors(cells[69],cells[989],cells[1604],cells[1085]);
	  set_offsets(cells[837],0,1,0,0);
	  cells[838]->set_vertices(vertices[103],vertices[55],vertices[143],vertices[96]);
	  cells[838]->set_neighbors(cells[899],cells[489],cells[132],cells[903]);
	  set_offsets(cells[838],2,2,0,3);
	  cells[839]->set_vertices(vertices[95],vertices[136],vertices[96],vertices[48]);
	  cells[839]->set_neighbors(cells[350],cells[576],cells[727],cells[671]);
	  set_offsets(cells[839],0,1,3,3);
	  cells[840]->set_vertices(vertices[187],vertices[194],vertices[146],vertices[147]);
	  cells[840]->set_neighbors(cells[922],cells[1136],cells[1420],cells[1115]);
	  set_offsets(cells[840],0,2,2,2);
	  cells[841]->set_vertices(vertices[170],vertices[162],vertices[115],vertices[122]);
	  cells[841]->set_neighbors(cells[517],cells[924],cells[818],cells[804]);
	  set_offsets(cells[841],0,0,0,0);
	  cells[842]->set_vertices(vertices[105],vertices[112],vertices[152],vertices[104]);
	  cells[842]->set_neighbors(cells[680],cells[593],cells[227],cells[678]);
	  set_offsets(cells[842],0,0,0,0);
	  cells[843]->set_vertices(vertices[155],vertices[107],vertices[147],vertices[156]);
	  cells[843]->set_neighbors(cells[153],cells[1029],cells[866],cells[961]);
	  set_offsets(cells[843],0,0,0,0);
	  cells[844]->set_vertices(vertices[126],vertices[117],vertices[166],vertices[118]);
	  cells[844]->set_neighbors(cells[224],cells[1013],cells[603],cells[463]);
	  set_offsets(cells[844],0,0,0,0);
	  cells[845]->set_vertices(vertices[123],vertices[115],vertices[172],vertices[124]);
	  cells[845]->set_neighbors(cells[56],cells[572],cells[548],cells[826]);
	  set_offsets(cells[845],0,0,0,0);
	  cells[846]->set_vertices(vertices[133],vertices[125],vertices[173],vertices[182]);
	  cells[846]->set_neighbors(cells[589],cells[986],cells[590],cells[1016]);
	  set_offsets(cells[846],0,0,0,0);
	  cells[847]->set_vertices(vertices[182],vertices[222],vertices[175],vertices[174]);
	  cells[847]->set_neighbors(cells[1310],cells[422],cells[1071],cells[1357]);
	  set_offsets(cells[847],0,0,0,0);
	  cells[848]->set_vertices(vertices[135],vertices[95],vertices[142],vertices[87]);
	  cells[848]->set_neighbors(cells[849],cells[140],cells[531],cells[902]);
	  set_offsets(cells[848],0,0,0,0);
	  cells[849]->set_vertices(vertices[142],vertices[95],vertices[94],vertices[87]);
	  cells[849]->set_neighbors(cells[606],cells[850],cells[848],cells[351]);
	  set_offsets(cells[849],0,0,0,0);
	  cells[850]->set_vertices(vertices[134],vertices[142],vertices[94],vertices[87]);
	  cells[850]->set_neighbors(cells[849],cells[218],cells[140],cells[880]);
	  set_offsets(cells[850],0,0,0,0);
	  cells[851]->set_vertices(vertices[183],vertices[230],vertices[223],vertices[175]);
	  cells[851]->set_neighbors(cells[1351],cells[753],cells[1356],cells[1297]);
	  set_offsets(cells[851],0,0,0,0);
	  cells[852]->set_vertices(vertices[196],vertices[148],vertices[149],vertices[156]);
	  cells[852]->set_neighbors(cells[417],cells[37],cells[1192],cells[862]);
	  set_offsets(cells[852],0,0,0,0);
	  cells[853]->set_vertices(vertices[184],vertices[137],vertices[129],vertices[177]);
	  cells[853]->set_neighbors(cells[878],cells[542],cells[1090],cells[445]);
	  set_offsets(cells[853],0,0,0,0);
	  cells[854]->set_vertices(vertices[130],vertices[170],vertices[123],vertices[122]);
	  cells[854]->set_neighbors(cells[924],cells[594],cells[521],cells[1051]);
	  set_offsets(cells[854],0,0,0,0);
	  cells[855]->set_vertices(vertices[98],vertices[89],vertices[138],vertices[50]);
	  cells[855]->set_neighbors(cells[821],cells[696],cells[631],cells[856]);
	  set_offsets(cells[855],2,0,0,2);
	  cells[856]->set_vertices(vertices[98],vertices[89],vertices[137],vertices[138]);
	  cells[856]->set_neighbors(cells[161],cells[688],cells[855],cells[859]);
	  set_offsets(cells[856],2,0,0,0);
	  cells[857]->set_vertices(vertices[96],vertices[89],vertices[136],vertices[137]);
	  cells[857]->set_neighbors(cells[815],cells[685],cells[858],cells[350]);
	  set_offsets(cells[857],2,0,0,0);
	  cells[858]->set_vertices(vertices[96],vertices[49],vertices[89],vertices[137]);
	  cells[858]->set_neighbors(cells[859],cells[857],cells[861],cells[111]);
	  set_offsets(cells[858],2,2,0,0);
	  cells[859]->set_vertices(vertices[137],vertices[49],vertices[89],vertices[98]);
	  cells[859]->set_neighbors(cells[631],cells[856],cells[860],cells[858]);
	  set_offsets(cells[859],0,2,0,2);
	  cells[860]->set_vertices(vertices[97],vertices[49],vertices[137],vertices[98]);
	  cells[860]->set_neighbors(cells[859],cells[915],cells[483],cells[861]);
	  set_offsets(cells[860],2,2,0,2);
	  cells[861]->set_vertices(vertices[96],vertices[49],vertices[137],vertices[97]);
	  cells[861]->set_neighbors(cells[860],cells[67],cells[636],cells[858]);
	  set_offsets(cells[861],2,2,0,2);
	  cells[862]->set_vertices(vertices[189],vertices[196],vertices[148],vertices[149]);
	  cells[862]->set_neighbors(cells[852],cells[1146],cells[1429],cells[126]);
	  set_offsets(cells[862],0,2,2,2);
	  cells[863]->set_vertices(vertices[229],vertices[181],vertices[221],vertices[230]);
	  cells[863]->set_neighbors(cells[1109],cells[1166],cells[1168],cells[1375]);
	  set_offsets(cells[863],0,0,0,0);
	  cells[864]->set_vertices(vertices[194],vertices[147],vertices[235],vertices[195]);
	  cells[864]->set_neighbors(cells[1376],cells[1002],cells[1227],cells[1420]);
	  set_offsets(cells[864],2,2,0,2);
	  cells[865]->set_vertices(vertices[165],vertices[157],vertices[205],vertices[214]);
	  cells[865]->set_neighbors(cells[1294],cells[1250],cells[985],cells[1255]);
	  set_offsets(cells[865],0,0,0,0);
	  cells[866]->set_vertices(vertices[164],vertices[107],vertices[155],vertices[156]);
	  cells[866]->set_neighbors(cells[843],cells[1030],cells[761],cells[715]);
	  set_offsets(cells[866],0,0,0,0);
	  cells[867]->set_vertices(vertices[104],vertices[152],vertices[103],vertices[111]);
	  cells[867]->set_neighbors(cells[898],cells[698],cells[680],cells[908]);
	  set_offsets(cells[867],1,1,0,0);
	  cells[868]->set_vertices(vertices[151],vertices[103],vertices[144],vertices[152]);
	  cells[868]->set_neighbors(cells[908],cells[893],cells[898],cells[1033]);
	  set_offsets(cells[868],0,0,1,1);
	  cells[869]->set_vertices(vertices[127],vertices[174],vertices[167],vertices[119]);
	  cells[869]->set_neighbors(cells[24],cells[809],cells[196],cells[1077]);
	  set_offsets(cells[869],0,0,0,0);
	  cells[870]->set_vertices(vertices[180],vertices[171],vertices[220],vertices[172]);
	  cells[870]->set_neighbors(cells[233],cells[155],cells[1063],cells[259]);
	  set_offsets(cells[870],0,0,0,0);
	  cells[871]->set_vertices(vertices[85],vertices[132],vertices[92],vertices[140]);
	  cells[871]->set_neighbors(cells[369],cells[812],cells[625],cells[632]);
	  set_offsets(cells[871],0,0,0,0);
	  cells[872]->set_vertices(vertices[98],vertices[51],vertices[139],vertices[99]);
	  cells[872]->set_neighbors(cells[791],cells[824],cells[511],cells[54]);
	  set_offsets(cells[872],2,2,0,2);
	  cells[873]->set_vertices(vertices[226],vertices[169],vertices[217],vertices[218]);
	  cells[873]->set_neighbors(cells[1307],cells[1131],cells[4],cells[1186]);
	  set_offsets(cells[873],0,0,0,0);
	  cells[874]->set_vertices(vertices[34],vertices[275],vertices[267],vertices[27]);
	  cells[874]->set_neighbors(cells[1389],cells[1565],cells[1660],cells[1524]);
	  set_offsets(cells[874],4,0,0,4);
	  cells[875]->set_vertices(vertices[14],vertices[262],vertices[255],vertices[22]);
	  cells[875]->set_neighbors(cells[1025],cells[1494],cells[222],cells[1234]);
	  set_offsets(cells[875],4,0,0,4);
	  cells[876]->set_vertices(vertices[258],vertices[201],vertices[250],vertices[210]);
	  cells[876]->set_neighbors(cells[1500],cells[1285],cells[1288],cells[1562]);
	  set_offsets(cells[876],0,0,0,0);
	  cells[877]->set_vertices(vertices[186],vertices[137],vertices[185],vertices[177]);
	  cells[877]->set_neighbors(cells[1090],cells[1173],cells[878],cells[1125]);
	  set_offsets(cells[877],0,0,0,0);
	  cells[878]->set_vertices(vertices[186],vertices[137],vertices[177],vertices[129]);
	  cells[878]->set_neighbors(cells[853],cells[946],cells[652],cells[877]);
	  set_offsets(cells[878],0,0,0,0);
	  cells[879]->set_vertices(vertices[93],vertices[142],vertices[141],vertices[102]);
	  cells[879]->set_neighbors(cells[1046],cells[885],cells[798],cells[882]);
	  set_offsets(cells[879],0,0,0,2);
	  cells[880]->set_vertices(vertices[142],vertices[94],vertices[85],vertices[134]);
	  cells[880]->set_neighbors(cells[142],cells[389],cells[850],cells[667]);
	  set_offsets(cells[880],0,0,0,0);
	  cells[881]->set_vertices(vertices[93],vertices[54],vertices[94],vertices[142]);
	  cells[881]->set_neighbors(cells[351],cells[667],cells[798],cells[255]);
	  set_offsets(cells[881],0,2,0,0);
	  cells[882]->set_vertices(vertices[93],vertices[133],vertices[141],vertices[142]);
	  cells[882]->set_neighbors(cells[805],cells[879],cells[834],cells[884]);
	  set_offsets(cells[882],0,0,0,0);
	  cells[883]->set_vertices(vertices[100],vertices[93],vertices[140],vertices[141]);
	  cells[883]->set_neighbors(cells[884],cells[1113],cells[886],cells[788]);
	  set_offsets(cells[883],2,0,0,0);
	  cells[884]->set_vertices(vertices[141],vertices[93],vertices[140],vertices[133]);
	  cells[884]->set_neighbors(cells[543],cells[695],cells[882],cells[883]);
	  set_offsets(cells[884],0,0,0,0);
	  cells[885]->set_vertices(vertices[53],vertices[102],vertices[93],vertices[141]);
	  cells[885]->set_neighbors(cells[879],cells[886],cells[888],cells[668]);
	  set_offsets(cells[885],2,2,0,0);
	  cells[886]->set_vertices(vertices[53],vertices[141],vertices[93],vertices[100]);
	  cells[886]->set_neighbors(cells[883],cells[651],cells[887],cells[885]);
	  set_offsets(cells[886],2,0,0,2);
	  cells[887]->set_vertices(vertices[53],vertices[101],vertices[141],vertices[100]);
	  cells[887]->set_neighbors(cells[707],cells[886],cells[382],cells[888]);
	  set_offsets(cells[887],2,2,0,2);
	  cells[888]->set_vertices(vertices[53],vertices[102],vertices[141],vertices[101]);
	  cells[888]->set_neighbors(cells[937],cells[887],cells[427],cells[885]);
	  set_offsets(cells[888],2,2,0,2);
	  cells[889]->set_vertices(vertices[225],vertices[234],vertices[177],vertices[226]);
	  cells[889]->set_neighbors(cells[1185],cells[1187],cells[1655],cells[1416]);
	  set_offsets(cells[889],0,0,0,0);
	  cells[890]->set_vertices(vertices[243],vertices[4],vertices[43],vertices[3]);
	  cells[890]->set_neighbors(cells[331],cells[1228],cells[1457],cells[1456]);
	  set_offsets(cells[890],2,6,4,6);
	  cells[891]->set_vertices(vertices[221],vertices[268],vertices[261],vertices[213]);
	  cells[891]->set_neighbors(cells[522],cells[127],cells[1381],cells[1633]);
	  set_offsets(cells[891],0,0,0,0);
	  cells[892]->set_vertices(vertices[191],vertices[103],vertices[143],vertices[144]);
	  cells[892]->set_neighbors(cells[489],cells[1153],cells[1033],cells[1158]);
	  set_offsets(cells[892],0,2,0,3);
	  cells[893]->set_vertices(vertices[151],vertices[144],vertices[192],vertices[152]);
	  cells[893]->set_neighbors(cells[686],cells[455],cells[868],cells[909]);
	  set_offsets(cells[893],0,1,1,1);
	  cells[894]->set_vertices(vertices[99],vertices[100],vertices[148],vertices[108]);
	  cells[894]->set_neighbors(cells[734],cells[720],cells[64],cells[648]);
	  set_offsets(cells[894],0,0,0,0);
	  cells[895]->set_vertices(vertices[141],vertices[102],vertices[190],vertices[150]);
	  cells[895]->set_neighbors(cells[757],cells[1145],cells[937],cells[1046]);
	  set_offsets(cells[895],0,2,0,2);
	  cells[896]->set_vertices(vertices[180],vertices[188],vertices[179],vertices[131]);
	  cells[896]->set_neighbors(cells[1061],cells[1095],cells[621],cells[925]);
	  set_offsets(cells[896],0,0,0,0);
	  cells[897]->set_vertices(vertices[188],vertices[139],vertices[187],vertices[179]);
	  cells[897]->set_neighbors(cells[1098],cells[1392],cells[1061],cells[1134]);
	  set_offsets(cells[897],0,0,0,0);
	  cells[898]->set_vertices(vertices[151],vertices[103],vertices[152],vertices[111]);
	  cells[898]->set_neighbors(cells[867],cells[241],cells[409],cells[868]);
	  set_offsets(cells[898],0,0,1,0);
	  cells[899]->set_vertices(vertices[143],vertices[55],vertices[95],vertices[96]);
	  cells[899]->set_neighbors(cells[576],cells[671],cells[838],cells[900]);
	  set_offsets(cells[899],0,2,0,3);
	  cells[900]->set_vertices(vertices[102],vertices[55],vertices[95],vertices[143]);
	  cells[900]->set_neighbors(cells[899],cells[901],cells[903],cells[669]);
	  set_offsets(cells[900],2,2,0,0);
	  cells[901]->set_vertices(vertices[143],vertices[95],vertices[102],vertices[142]);
	  cells[901]->set_neighbors(cells[675],cells[1118],cells[902],cells[900]);
	  set_offsets(cells[901],0,0,2,0);
	  cells[902]->set_vertices(vertices[135],vertices[95],vertices[143],vertices[142]);
	  cells[902]->set_neighbors(cells[901],cells[1082],cells[848],cells[672]);
	  set_offsets(cells[902],0,0,0,0);
	  cells[903]->set_vertices(vertices[102],vertices[55],vertices[143],vertices[103]);
	  cells[903]->set_neighbors(cells[838],cells[938],cells[657],cells[900]);
	  set_offsets(cells[903],2,2,0,2);
	  cells[904]->set_vertices(vertices[184],vertices[143],vertices[183],vertices[135]);
	  cells[904]->set_neighbors(cells[920],cells[998],cells[816],cells[1154]);
	  set_offsets(cells[904],1,0,0,0);
	  cells[905]->set_vertices(vertices[211],vertices[163],vertices[203],vertices[212]);
	  cells[905]->set_neighbors(cells[441],cells[1336],cells[1062],cells[364]);
	  set_offsets(cells[905],0,0,0,0);
	  cells[906]->set_vertices(vertices[145],vertices[97],vertices[105],vertices[152]);
	  cells[906]->set_neighbors(cells[593],cells[951],cells[614],cells[557]);
	  set_offsets(cells[906],0,0,0,0);
	  cells[907]->set_vertices(vertices[99],vertices[108],vertices[156],vertices[107]);
	  cells[907]->set_neighbors(cells[756],cells[153],cells[630],cells[720]);
	  set_offsets(cells[907],0,0,0,0);
	  cells[908]->set_vertices(vertices[104],vertices[144],vertices[103],vertices[152]);
	  cells[908]->set_neighbors(cells[868],cells[867],cells[797],cells[638]);
	  set_offsets(cells[908],1,1,0,1);
	  cells[909]->set_vertices(vertices[192],vertices[144],vertices[151],vertices[191]);
	  cells[909]->set_neighbors(cells[1033],cells[345],cells[1259],cells[893]);
	  set_offsets(cells[909],3,3,2,0);
	  cells[910]->set_vertices(vertices[32],vertices[265],vertices[272],vertices[24]);
	  cells[910]->set_neighbors(cells[1318],cells[1492],cells[1610],cells[1657]);
	  set_offsets(cells[910],4,0,0,4);
	  cells[911]->set_vertices(vertices[278],vertices[269],vertices[30],vertices[270]);
	  cells[911]->set_neighbors(cells[1390],cells[1589],cells[1479],cells[346]);
	  set_offsets(cells[911],0,0,4,0);
	  cells[912]->set_vertices(vertices[251],vertices[18],vertices[11],vertices[10]);
	  cells[912]->set_neighbors(cells[106],cells[9],cells[749],cells[1518]);
	  set_offsets(cells[912],0,4,4,4);
	  cells[913]->set_vertices(vertices[232],vertices[184],vertices[224],vertices[177]);
	  cells[913]->set_neighbors(cells[1084],cells[1396],cells[1347],cells[1398]);
	  set_offsets(cells[913],0,0,0,0);
	  cells[914]->set_vertices(vertices[185],vertices[144],vertices[145],vertices[192]);
	  cells[914]->set_neighbors(cells[686],cells[1409],cells[1053],cells[1130]);
	  set_offsets(cells[914],0,2,2,2);
	  cells[915]->set_vertices(vertices[146],vertices[97],vertices[137],vertices[98]);
	  cells[915]->set_neighbors(cells[860],cells[1124],cells[626],cells[1089]);
	  set_offsets(cells[915],2,2,0,2);
	  cells[916]->set_vertices(vertices[227],vertices[226],vertices[274],vertices[219]);
	  cells[916]->set_neighbors(cells[529],cells[1440],cells[1372],cells[1654]);
	  set_offsets(cells[916],0,0,0,0);
	  cells[917]->set_vertices(vertices[171],vertices[123],vertices[163],vertices[172]);
	  cells[917]->set_neighbors(cells[826],cells[233],cells[1063],cells[1055]);
	  set_offsets(cells[917],0,0,0,0);
	  cells[918]->set_vertices(vertices[174],vertices[165],vertices[214],vertices[166]);
	  cells[918]->set_neighbors(cells[985],cells[1075],cells[1042],cells[256]);
	  set_offsets(cells[918],0,0,0,0);
	  cells[919]->set_vertices(vertices[216],vertices[208],vertices[167],vertices[215]);
	  cells[919]->set_neighbors(cells[1204],cells[69],cells[1545],cells[1312]);
	  set_offsets(cells[919],1,1,0,0);
	  cells[920]->set_vertices(vertices[183],vertices[143],vertices[190],vertices[135]);
	  cells[920]->set_neighbors(cells[1082],cells[760],cells[904],cells[1155]);
	  set_offsets(cells[920],0,0,0,0);
	  cells[921]->set_vertices(vertices[238],vertices[190],vertices[189],vertices[181]);
	  cells[921]->set_neighbors(cells[1067],cells[1163],cells[348],cells[1167]);
	  set_offsets(cells[921],0,0,0,0);
	  cells[922]->set_vertices(vertices[147],vertices[146],vertices[154],vertices[194]);
	  cells[922]->set_neighbors(cells[719],cells[1191],cells[840],cells[461]);
	  set_offsets(cells[922],0,0,0,0);
	  cells[923]->set_vertices(vertices[178],vertices[170],vertices[171],vertices[123]);
	  cells[923]->set_neighbors(cells[1055],cells[1054],cells[1051],cells[1080]);
	  set_offsets(cells[923],0,0,0,0);
	  cells[924]->set_vertices(vertices[123],vertices[170],vertices[115],vertices[122]);
	  cells[924]->set_neighbors(cells[841],cells[474],cells[854],cells[687]);
	  set_offsets(cells[924],0,0,0,0);
	  cells[925]->set_vertices(vertices[180],vertices[188],vertices[228],vertices[179]);
	  cells[925]->set_neighbors(cells[133],cells[261],cells[896],cells[1325]);
	  set_offsets(cells[925],0,0,0,0);
	  cells[926]->set_vertices(vertices[189],vertices[141],vertices[188],vertices[181]);
	  cells[926]->set_neighbors(cells[825],cells[1114],cells[1067],cells[116]);
	  set_offsets(cells[926],0,0,0,0);
	  cells[927]->set_vertices(vertices[180],vertices[140],vertices[132],vertices[133]);
	  cells[927]->set_neighbors(cells[625],cells[1064],cells[1123],cells[1068]);
	  set_offsets(cells[927],0,0,0,0);
	  cells[928]->set_vertices(vertices[132],vertices[172],vertices[125],vertices[124]);
	  cells[928]->set_neighbors(cells[459],cells[472],cells[572],cells[711]);
	  set_offsets(cells[928],0,0,0,0);
	  cells[929]->set_vertices(vertices[4],vertices[243],vertices[244],vertices[252]);
	  cells[929]->set_neighbors(cells[1453],cells[1523],cells[1379],cells[1283]);
	  set_offsets(cells[929],4,0,0,0);
	  cells[930]->set_vertices(vertices[225],vertices[265],vertices[273],vertices[274]);
	  cells[930]->set_neighbors(cells[83],cells[1692],cells[1330],cells[1656]);
	  set_offsets(cells[930],0,0,0,0);
	  cells[931]->set_vertices(vertices[2],vertices[283],vertices[42],vertices[43]);
	  cells[931]->set_neighbors(cells[1139],cells[477],cells[1445],cells[1699]);
	  set_offsets(cells[931],6,0,4,4);
	  cells[932]->set_vertices(vertices[150],vertices[102],vertices[103],vertices[110]);
	  cells[932]->set_neighbors(cells[657],cells[253],cells[831],cells[938]);
	  set_offsets(cells[932],0,0,0,0);
	  cells[933]->set_vertices(vertices[125],vertices[117],vertices[174],vertices[126]);
	  cells[933]->set_neighbors(cells[463],cells[736],cells[400],cells[982]);
	  set_offsets(cells[933],0,0,0,0);
	  cells[934]->set_vertices(vertices[141],vertices[133],vertices[181],vertices[190]);
	  cells[934]->set_neighbors(cells[694],cells[1067],cells[805],cells[825]);
	  set_offsets(cells[934],0,0,0,0);
	  cells[935]->set_vertices(vertices[190],vertices[142],vertices[133],vertices[182]);
	  cells[935]->set_neighbors(cells[1105],cells[694],cells[1083],cells[805]);
	  set_offsets(cells[935],0,0,0,0);
	  cells[936]->set_vertices(vertices[231],vertices[278],vertices[279],vertices[271]);
	  cells[936]->set_neighbors(cells[1245],cells[1093],cells[1643],cells[1715]);
	  set_offsets(cells[936],0,0,0,0);
	  cells[937]->set_vertices(vertices[101],vertices[102],vertices[141],vertices[150]);
	  cells[937]->set_neighbors(cells[895],cells[1104],cells[831],cells[888]);
	  set_offsets(cells[937],2,2,0,2);
	  cells[938]->set_vertices(vertices[150],vertices[102],vertices[143],vertices[103]);
	  cells[938]->set_neighbors(cells[903],cells[1158],cells[932],cells[757]);
	  set_offsets(cells[938],2,2,0,2);
	  cells[939]->set_vertices(vertices[23],vertices[263],vertices[16],vertices[24]);
	  cells[939]->set_neighbors(cells[1233],cells[297],cells[1590],cells[1235]);
	  set_offsets(cells[939],4,0,5,5);
	  cells[940]->set_vertices(vertices[256],vertices[201],vertices[208],vertices[248]);
	  cells[940]->set_neighbors(cells[1210],cells[1508],cells[1543],cells[176]);
	  set_offsets(cells[940],0,0,0,0);
	  cells[941]->set_vertices(vertices[151],vertices[200],vertices[159],vertices[152]);
	  cells[941]->set_neighbors(cells[1214],cells[241],cells[455],cells[1222]);
	  set_offsets(cells[941],0,1,0,1);
	  cells[942]->set_vertices(vertices[101],vertices[150],vertices[149],vertices[158]);
	  cells[942]->set_neighbors(cells[643],cells[691],cells[740],cells[1147]);
	  set_offsets(cells[942],0,0,0,0);
	  cells[943]->set_vertices(vertices[228],vertices[219],vertices[268],vertices[220]);
	  cells[943]->set_neighbors(cells[1386],cells[1366],cells[962],cells[1374]);
	  set_offsets(cells[943],0,0,0,0);
	  cells[944]->set_vertices(vertices[14],vertices[254],vertices[247],vertices[255]);
	  cells[944]->set_neighbors(cells[1291],cells[708],cells[1234],cells[1239]);
	  set_offsets(cells[944],4,0,0,0);
	  cells[945]->set_vertices(vertices[3],vertices[243],vertices[251],vertices[10]);
	  cells[945]->set_neighbors(cells[1007],cells[9],cells[1423],cells[1468]);
	  set_offsets(cells[945],4,0,0,4);
	  cells[946]->set_vertices(vertices[178],vertices[186],vertices[177],vertices[129]);
	  cells[946]->set_neighbors(cells[878],cells[646],cells[656],cells[1039]);
	  set_offsets(cells[946],0,0,0,0);
	  cells[947]->set_vertices(vertices[210],vertices[250],vertices[203],vertices[202]);
	  cells[947]->set_neighbors(cells[665],cells[1236],cells[1500],cells[1285]);
	  set_offsets(cells[947],0,0,0,0);
	  cells[948]->set_vertices(vertices[273],vertices[34],vertices[25],vertices[33]);
	  cells[948]->set_neighbors(cells[74],cells[1412],cells[1438],cells[1659]);
	  set_offsets(cells[948],0,4,4,4);
	  cells[949]->set_vertices(vertices[113],vertices[120],vertices[160],vertices[112]);
	  cells[949]->set_neighbors(cells[978],cells[956],cells[272],cells[662]);
	  set_offsets(cells[949],0,0,0,0);
	  cells[950]->set_vertices(vertices[114],vertices[105],vertices[154],vertices[106]);
	  cells[950]->set_neighbors(cells[366],cells[501],cells[212],cells[793]);
	  set_offsets(cells[950],0,0,0,0);
	  cells[951]->set_vertices(vertices[153],vertices[105],vertices[152],vertices[145]);
	  cells[951]->set_neighbors(cells[906],cells[1102],cells[73],cells[43]);
	  set_offsets(cells[951],0,0,0,0);
	  cells[952]->set_vertices(vertices[244],vertices[195],vertices[196],vertices[204]);
	  cells[952]->set_neighbors(cells[309],cells[1244],cells[1377],cells[1458]);
	  set_offsets(cells[952],0,0,0,0);
	  cells[953]->set_vertices(vertices[173],vertices[172],vertices[165],vertices[125]);
	  cells[953]->set_neighbors(cells[91],cells[540],cells[993],cells[1302]);
	  set_offsets(cells[953],0,0,0,0);
	  cells[954]->set_vertices(vertices[170],vertices[113],vertices[161],vertices[162]);
	  cells[954]->set_neighbors(cells[746],cells[1008],cells[818],cells[139]);
	  set_offsets(cells[954],0,0,0,0);
	  cells[955]->set_vertices(vertices[113],vertices[105],vertices[160],vertices[153]);
	  cells[955]->set_neighbors(cells[43],cells[701],cells[512],cells[956]);
	  set_offsets(cells[955],0,0,0,0);
	  cells[956]->set_vertices(vertices[112],vertices[105],vertices[160],vertices[113]);
	  cells[956]->set_neighbors(cells[955],cells[949],cells[729],cells[678]);
	  set_offsets(cells[956],0,0,0,0);
	  cells[957]->set_vertices(vertices[279],vertices[287],vertices[239],vertices[280]);
	  cells[957]->set_neighbors(cells[1676],cells[1247],cells[1180],cells[1725]);
	  set_offsets(cells[957],0,0,0,1);
	  cells[958]->set_vertices(vertices[219],vertices[266],vertices[259],vertices[211]);
	  cells[958]->set_neighbors(cells[1564],cells[1575],cells[1334],cells[1622]);
	  set_offsets(cells[958],0,0,0,0);
	  cells[959]->set_vertices(vertices[178],vertices[169],vertices[218],vertices[170]);
	  cells[959]->set_neighbors(cells[165],cells[1080],cells[634],cells[4]);
	  set_offsets(cells[959],0,0,0,0);
	  cells[960]->set_vertices(vertices[209],vertices[161],vertices[201],vertices[210]);
	  cells[960]->set_neighbors(cells[1027],cells[1288],cells[979],cells[138]);
	  set_offsets(cells[960],0,0,0,0);
	  cells[961]->set_vertices(vertices[155],vertices[154],vertices[147],vertices[107]);
	  cells[961]->set_neighbors(cells[704],cells[843],cells[770],cells[100]);
	  set_offsets(cells[961],0,0,0,0);
	  cells[962]->set_vertices(vertices[228],vertices[171],vertices[219],vertices[220]);
	  cells[962]->set_neighbors(cells[1272],cells[943],cells[259],cells[1122]);
	  set_offsets(cells[962],0,0,0,0);
	  cells[963]->set_vertices(vertices[247],vertices[8],vertices[0],vertices[7]);
	  cells[963]->set_neighbors(cells[10],cells[1353],cells[1263],cells[1442]);
	  set_offsets(cells[963],0,5,5,4);
	  cells[964]->set_vertices(vertices[213],vertices[212],vertices[260],vertices[205]);
	  cells[964]->set_neighbors(cells[95],cells[1319],cells[1241],cells[1284]);
	  set_offsets(cells[964],0,0,0,0);
	  cells[965]->set_vertices(vertices[129],vertices[184],vertices[176],vertices[136]);
	  cells[965]->set_neighbors(cells[1074],cells[1073],cells[445],cells[542]);
	  set_offsets(cells[965],0,0,0,0);
	  cells[966]->set_vertices(vertices[4],vertices[12],vertices[245],vertices[5]);
	  cells[966]->set_neighbors(cells[1531],cells[236],cells[337],cells[1175]);
	  set_offsets(cells[966],4,4,0,4);
	  cells[967]->set_vertices(vertices[208],vertices[153],vertices[200],vertices[201]);
	  cells[967]->set_neighbors(cells[1229],cells[1210],cells[25],cells[997]);
	  set_offsets(cells[967],0,0,0,0);
	  cells[968]->set_vertices(vertices[163],vertices[210],vertices[203],vertices[155]);
	  cells[968]->set_neighbors(cells[1236],cells[441],cells[1010],cells[364]);
	  set_offsets(cells[968],0,0,0,0);
	  cells[969]->set_vertices(vertices[162],vertices[105],vertices[153],vertices[154]);
	  cells[969]->set_neighbors(cells[73],cells[783],cells[793],cells[512]);
	  set_offsets(cells[969],0,0,0,0);
	  cells[970]->set_vertices(vertices[124],vertices[115],vertices[164],vertices[116]);
	  cells[970]->set_neighbors(cells[802],cells[693],cells[663],cells[56]);
	  set_offsets(cells[970],0,0,0,0);
	  cells[971]->set_vertices(vertices[157],vertices[156],vertices[149],vertices[109]);
	  cells[971]->set_neighbors(cells[747],cells[164],cells[647],cells[1190]);
	  set_offsets(cells[971],0,0,0,0);
	  cells[972]->set_vertices(vertices[238],vertices[230],vertices[231],vertices[183]);
	  cells[972]->set_neighbors(cells[1297],cells[323],cells[419],cells[1675]);
	  set_offsets(cells[972],0,0,0,0);
	  cells[973]->set_vertices(vertices[191],vertices[150],vertices[238],vertices[190]);
	  cells[973]->set_neighbors(cells[1167],cells[325],cells[1156],cells[1431]);
	  set_offsets(cells[973],0,2,0,0);
	  cells[974]->set_vertices(vertices[23],vertices[22],vertices[263],vertices[30]);
	  cells[974]->set_neighbors(cells[1202],cells[1593],cells[396],cells[1525]);
	  set_offsets(cells[974],4,4,0,4);
	  cells[975]->set_vertices(vertices[0],vertices[248],vertices[241],vertices[8]);
	  cells[975]->set_neighbors(cells[1503],cells[1101],cells[1442],cells[229]);
	  set_offsets(cells[975],4,0,0,4);
	  cells[976]->set_vertices(vertices[147],vertices[156],vertices[196],vertices[204]);
	  cells[976]->set_neighbors(cells[37],cells[309],cells[1029],cells[1192]);
	  set_offsets(cells[976],0,0,0,0);
	  cells[977]->set_vertices(vertices[176],vertices[168],vertices[127],vertices[175]);
	  cells[977]->set_neighbors(cells[1078],cells[1032],cells[1314],cells[689]);
	  set_offsets(cells[977],1,1,0,0);
	  cells[978]->set_vertices(vertices[160],vertices[120],vertices[119],vertices[112]);
	  cells[978]->set_neighbors(cells[731],cells[762],cells[949],cells[987]);
	  set_offsets(cells[978],1,1,0,1);
	  cells[979]->set_vertices(vertices[218],vertices[161],vertices[209],vertices[210]);
	  cells[979]->set_neighbors(cells[960],cells[1517],cells[248],cells[1049]);
	  set_offsets(cells[979],0,0,0,0);
	  cells[980]->set_vertices(vertices[165],vertices[117],vertices[157],vertices[166]);
	  cells[980]->set_neighbors(cells[268],cells[985],cells[1042],cells[1018]);
	  set_offsets(cells[980],0,0,0,0);
	  cells[981]->set_vertices(vertices[231],vertices[224],vertices[183],vertices[223]);
	  cells[981]->set_neighbors(cells[753],cells[1297],cells[995],cells[1393]);
	  set_offsets(cells[981],0,1,0,0);
	  cells[982]->set_vertices(vertices[125],vertices[117],vertices[165],vertices[174]);
	  cells[982]->set_neighbors(cells[1042],cells[540],cells[933],cells[91]);
	  set_offsets(cells[982],0,0,0,0);
	  cells[983]->set_vertices(vertices[172],vertices[212],vertices[165],vertices[164]);
	  cells[983]->set_neighbors(cells[1254],cells[712],cells[833],cells[1301]);
	  set_offsets(cells[983],0,0,0,0);
	  cells[984]->set_vertices(vertices[172],vertices[164],vertices[117],vertices[124]);
	  cells[984]->set_neighbors(cells[693],cells[459],cells[56],cells[712]);
	  set_offsets(cells[984],0,0,0,0);
	  cells[985]->set_vertices(vertices[165],vertices[157],vertices[214],vertices[166]);
	  cells[985]->set_neighbors(cells[1040],cells[918],cells[980],cells[865]);
	  set_offsets(cells[985],0,0,0,0);
	  cells[986]->set_vertices(vertices[173],vertices[181],vertices[133],vertices[182]);
	  cells[986]->set_neighbors(cells[694],cells[846],cells[194],cells[1107]);
	  set_offsets(cells[986],0,0,0,0);
	  cells[987]->set_vertices(vertices[120],vertices[160],vertices[119],vertices[168]);
	  cells[987]->set_neighbors(cells[741],cells[1019],cells[662],cells[978]);
	  set_offsets(cells[987],1,1,0,1);
	  cells[988]->set_vertices(vertices[225],vertices[185],vertices[232],vertices[177]);
	  cells[988]->set_neighbors(cells[1347],cells[1396],cells[1416],cells[1407]);
	  set_offsets(cells[988],0,0,0,0);
	  cells[989]->set_vertices(vertices[223],vertices[222],vertices[215],vertices[175]);
	  cells[989]->set_neighbors(cells[810],cells[837],cells[1351],cells[1592]);
	  set_offsets(cells[989],0,0,0,0);
	  cells[990]->set_vertices(vertices[176],vertices[128],vertices[168],vertices[121]);
	  cells[990]->set_neighbors(cells[733],cells[571],cells[516],cells[689]);
	  set_offsets(cells[990],0,0,0,0);
	  cells[991]->set_vertices(vertices[120],vertices[113],vertices[168],vertices[121]);
	  cells[991]->set_neighbors(cells[992],cells[733],cells[766],cells[662]);
	  set_offsets(cells[991],0,0,0,0);
	  cells[992]->set_vertices(vertices[121],vertices[113],vertices[168],vertices[161]);
	  cells[992]->set_neighbors(cells[744],cells[1044],cells[139],cells[991]);
	  set_offsets(cells[992],0,0,0,0);
	  cells[993]->set_vertices(vertices[180],vertices[172],vertices[173],vertices[125]);
	  cells[993]->set_neighbors(cells[953],cells[1016],cells[711],cells[155]);
	  set_offsets(cells[993],0,0,0,0);
	  cells[994]->set_vertices(vertices[168],vertices[161],vertices[216],vertices[169]);
	  cells[994]->set_neighbors(cells[1274],cells[16],cells[1044],cells[216]);
	  set_offsets(cells[994],0,0,0,0);
	  cells[995]->set_vertices(vertices[231],vertices[224],vertices[223],vertices[272]);
	  cells[995]->set_neighbors(cells[1170],cells[1299],cells[1410],cells[981]);
	  set_offsets(cells[995],0,1,0,1);
	  cells[996]->set_vertices(vertices[175],vertices[183],vertices[182],vertices[135]);
	  cells[996]->set_neighbors(cells[760],cells[1076],cells[1117],cells[1356]);
	  set_offsets(cells[996],0,0,0,0);
	  cells[997]->set_vertices(vertices[160],vertices[153],vertices[200],vertices[208]);
	  cells[997]->set_neighbors(cells[967],cells[1267],cells[1217],cells[1112]);
	  set_offsets(cells[997],0,0,0,0);
	  cells[998]->set_vertices(vertices[176],vertices[184],vertices[183],vertices[135]);
	  cells[998]->set_neighbors(cells[904],cells[1117],cells[1074],cells[1352]);
	  set_offsets(cells[998],1,1,0,0);
	  cells[999]->set_vertices(vertices[178],vertices[138],vertices[130],vertices[131]);
	  cells[999]->set_neighbors(cells[829],cells[1050],cells[1127],cells[1017]);
	  set_offsets(cells[999],0,0,0,0);
	  cells[1000]->set_vertices(vertices[193],vertices[194],vertices[202],vertices[145]);
	  cells[1000]->set_neighbors(cells[1150],cells[1149],cells[1371],cells[1454]);
	  set_offsets(cells[1000],0,0,0,0);
	  cells[1001]->set_vertices(vertices[163],vertices[162],vertices[155],vertices[115]);
	  cells[1001]->set_neighbors(cells[645],cells[703],cells[804],cells[1010]);
	  set_offsets(cells[1001],0,0,0,0);
	  cells[1002]->set_vertices(vertices[242],vertices[194],vertices[235],vertices[195]);
	  cells[1002]->set_neighbors(cells[864],cells[1700],cells[1132],cells[84]);
	  set_offsets(cells[1002],2,2,0,2);
	  cells[1003]->set_vertices(vertices[277],vertices[237],vertices[284],vertices[229]);
	  cells[1003]->set_neighbors(cells[1463],cells[1661],cells[623],cells[1707]);
	  set_offsets(cells[1003],0,0,0,0);
	  cells[1004]->set_vertices(vertices[43],vertices[44],vertices[283],vertices[35]);
	  cells[1004]->set_neighbors(cells[1666],cells[1139],cells[41],cells[1364]);
	  set_offsets(cells[1004],4,4,0,4);
	  cells[1005]->set_vertices(vertices[235],vertices[187],vertices[227],vertices[236]);
	  cells[1005]->set_neighbors(cells[1200],cells[1459],cells[1196],cells[1418]);
	  set_offsets(cells[1005],0,0,0,0);
	  cells[1006]->set_vertices(vertices[249],vertices[10],vertices[9],vertices[18]);
	  cells[1006]->set_neighbors(cells[28],cells[1556],cells[44],cells[1507]);
	  set_offsets(cells[1006],0,4,4,4);
	  cells[1007]->set_vertices(vertices[243],vertices[250],vertices[251],vertices[10]);
	  cells[1007]->set_neighbors(cells[1277],cells[945],cells[1223],cells[1519]);
	  set_offsets(cells[1007],0,0,0,4);
	  cells[1008]->set_vertices(vertices[170],vertices[161],vertices[210],vertices[162]);
	  cells[1008]->set_neighbors(cells[158],cells[775],cells[954],cells[248]);
	  set_offsets(cells[1008],0,0,0,0);
	  cells[1009]->set_vertices(vertices[2],vertices[242],vertices[243],vertices[250]);
	  cells[1009]->set_neighbors(cells[280],cells[1223],cells[1511],cells[1619]);
	  set_offsets(cells[1009],4,0,0,0);
	  cells[1010]->set_vertices(vertices[163],vertices[210],vertices[155],vertices[162]);
	  cells[1010]->set_neighbors(cells[250],cells[1001],cells[775],cells[968]);
	  set_offsets(cells[1010],0,0,0,0);
	  cells[1011]->set_vertices(vertices[203],vertices[155],vertices[195],vertices[204]);
	  cells[1011]->set_neighbors(cells[1012],cells[1159],cells[1281],cells[1237]);
	  set_offsets(cells[1011],0,0,0,0);
	  cells[1012]->set_vertices(vertices[147],vertices[204],vertices[195],vertices[155]);
	  cells[1012]->set_neighbors(cells[1011],cells[1034],cells[1029],cells[309]);
	  set_offsets(cells[1012],0,0,0,0);
	  cells[1013]->set_vertices(vertices[126],vertices[166],vertices[119],vertices[118]);
	  cells[1013]->set_neighbors(cells[615],cells[755],cells[844],cells[785]);
	  set_offsets(cells[1013],0,0,0,0);
	  cells[1014]->set_vertices(vertices[164],vertices[204],vertices[157],vertices[156]);
	  cells[1014]->set_neighbors(cells[1190],cells[647],cells[1030],cells[1258]);
	  set_offsets(cells[1014],0,0,0,0);
	  cells[1015]->set_vertices(vertices[31],vertices[38],vertices[279],vertices[39]);
	  cells[1015]->set_neighbors(cells[1671],cells[0],cells[552],cells[1238]);
	  set_offsets(cells[1015],4,4,0,4);
	  cells[1016]->set_vertices(vertices[133],vertices[180],vertices[173],vertices[125]);
	  cells[1016]->set_neighbors(cells[993],cells[846],cells[1064],cells[1107]);
	  set_offsets(cells[1016],0,0,0,0);
	  cells[1017]->set_vertices(vertices[129],vertices[138],vertices[130],vertices[178]);
	  cells[1017]->set_neighbors(cells[999],cells[1065],cells[656],cells[78]);
	  set_offsets(cells[1017],0,0,0,0);
	  cells[1018]->set_vertices(vertices[165],vertices[164],vertices[157],vertices[117]);
	  cells[1018]->set_neighbors(cells[191],cells[980],cells[712],cells[1254]);
	  set_offsets(cells[1018],0,0,0,0);
	  cells[1019]->set_vertices(vertices[120],vertices[168],vertices[119],vertices[127]);
	  cells[1019]->set_neighbors(cells[809],cells[654],cells[1038],cells[987]);
	  set_offsets(cells[1019],1,1,0,0);
	  cells[1020]->set_vertices(vertices[260],vertices[203],vertices[251],vertices[252]);
	  cells[1020]->set_neighbors(cells[1257],cells[1576],cells[1024],cells[1252]);
	  set_offsets(cells[1020],0,0,0,0);
	  cells[1021]->set_vertices(vertices[39],vertices[46],vertices[287],vertices[47]);
	  cells[1021]->set_neighbors(cells[1488],cells[1595],cells[384],cells[1727]);
	  set_offsets(cells[1021],4,4,0,4);
	  cells[1022]->set_vertices(vertices[286],vertices[198],vertices[239],vertices[246]);
	  cells[1022]->set_neighbors(cells[1486],cells[1634],cells[1478],cells[1473]);
	  set_offsets(cells[1022],0,2,0,2);
	  cells[1023]->set_vertices(vertices[268],vertices[259],vertices[20],vertices[260]);
	  cells[1023]->set_neighbors(cells[1280],cells[217],cells[1331],cells[1298]);
	  set_offsets(cells[1023],0,0,4,0);
	  cells[1024]->set_vertices(vertices[260],vertices[203],vertices[252],vertices[212]);
	  cells[1024]->set_neighbors(cells[1066],cells[95],cells[1336],cells[1020]);
	  set_offsets(cells[1024],0,0,0,0);
	  cells[1025]->set_vertices(vertices[263],vertices[262],vertices[22],vertices[255]);
	  cells[1025]->set_neighbors(cells[875],cells[1537],cells[1308],cells[264]);
	  set_offsets(cells[1025],0,0,4,0);
	  cells[1026]->set_vertices(vertices[245],vertices[197],vertices[246],vertices[254]);
	  cells[1026]->set_neighbors(cells[1212],cells[197],cells[1265],cells[1669]);
	  set_offsets(cells[1026],0,0,0,0);
	  cells[1027]->set_vertices(vertices[161],vertices[153],vertices[201],vertices[210]);
	  cells[1027]->set_neighbors(cells[811],cells[960],cells[158],cells[25]);
	  set_offsets(cells[1027],0,0,0,0);
	  cells[1028]->set_vertices(vertices[199],vertices[246],vertices[198],vertices[206]);
	  cells[1028]->set_neighbors(cells[1213],cells[1240],cells[1300],cells[1486]);
	  set_offsets(cells[1028],0,0,0,0);
	  cells[1029]->set_vertices(vertices[147],vertices[156],vertices[204],vertices[155]);
	  cells[1029]->set_neighbors(cells[1030],cells[1012],cells[843],cells[976]);
	  set_offsets(cells[1029],0,0,0,0);
	  cells[1030]->set_vertices(vertices[164],vertices[155],vertices[204],vertices[156]);
	  cells[1030]->set_neighbors(cells[1029],cells[1014],cells[866],cells[59]);
	  set_offsets(cells[1030],0,0,0,0);
	  cells[1031]->set_vertices(vertices[166],vertices[206],vertices[159],vertices[158]);
	  cells[1031]->set_neighbors(cells[1111],cells[506],cells[1041],cells[1264]);
	  set_offsets(cells[1031],0,0,0,0);
	  cells[1032]->set_vertices(vertices[135],vertices[176],vertices[127],vertices[175]);
	  cells[1032]->set_neighbors(cells[977],cells[1076],cells[1117],cells[1035]);
	  set_offsets(cells[1032],0,1,0,0);
	  cells[1033]->set_vertices(vertices[151],vertices[103],vertices[191],vertices[144]);
	  cells[1033]->set_neighbors(cells[892],cells[909],cells[868],cells[1157]);
	  set_offsets(cells[1033],2,2,0,3);
	  cells[1034]->set_vertices(vertices[147],vertices[202],vertices[155],vertices[195]);
	  cells[1034]->set_neighbors(cells[1237],cells[1012],cells[1227],cells[100]);
	  set_offsets(cells[1034],0,0,0,0);
	  cells[1035]->set_vertices(vertices[135],vertices[128],vertices[127],vertices[176]);
	  cells[1035]->set_neighbors(cells[689],cells[1032],cells[782],cells[758]);
	  set_offsets(cells[1035],0,1,0,1);
	  cells[1036]->set_vertices(vertices[167],vertices[160],vertices[119],vertices[159]);
	  cells[1036]->set_neighbors(cells[360],cells[1037],cells[1268],cells[741]);
	  set_offsets(cells[1036],0,1,0,0);
	  cells[1037]->set_vertices(vertices[167],vertices[166],vertices[159],vertices[119]);
	  cells[1037]->set_neighbors(cells[750],cells[1036],cells[24],cells[1160]);
	  set_offsets(cells[1037],0,0,0,0);
	  cells[1038]->set_vertices(vertices[168],vertices[128],vertices[127],vertices[120]);
	  cells[1038]->set_neighbors(cells[88],cells[1019],cells[733],cells[689]);
	  set_offsets(cells[1038],1,1,0,1);
	  cells[1039]->set_vertices(vertices[178],vertices[186],vertices[226],vertices[177]);
	  cells[1039]->set_neighbors(cells[1185],cells[47],cells[946],cells[771]);
	  set_offsets(cells[1039],0,0,0,0);
	  cells[1040]->set_vertices(vertices[214],vertices[157],vertices[206],vertices[166]);
	  cells[1040]->set_neighbors(cells[1041],cells[1264],cells[985],cells[1294]);
	  set_offsets(cells[1040],0,0,0,0);
	  cells[1041]->set_vertices(vertices[166],vertices[157],vertices[206],vertices[158]);
	  cells[1041]->set_neighbors(cells[692],cells[1031],cells[182],cells[1040]);
	  set_offsets(cells[1041],0,0,0,0);
	  cells[1042]->set_vertices(vertices[174],vertices[117],vertices[165],vertices[166]);
	  cells[1042]->set_neighbors(cells[980],cells[918],cells[463],cells[982]);
	  set_offsets(cells[1042],0,0,0,0);
	  cells[1043]->set_vertices(vertices[176],vertices[177],vertices[129],vertices[169]);
	  cells[1043]->set_neighbors(cells[646],cells[713],cells[1358],cells[542]);
	  set_offsets(cells[1043],0,0,0,0);
	  cells[1044]->set_vertices(vertices[169],vertices[121],vertices[168],vertices[161]);
	  cells[1044]->set_neighbors(cells[992],cells[994],cells[444],cells[571]);
	  set_offsets(cells[1044],0,0,0,0);
	  cells[1045]->set_vertices(vertices[36],vertices[267],vertices[28],vertices[276]);
	  cells[1045]->set_neighbors(cells[664],cells[1632],cells[1382],cells[1574]);
	  set_offsets(cells[1045],4,0,4,0);
	  cells[1046]->set_vertices(vertices[141],vertices[102],vertices[142],vertices[190]);
	  cells[1046]->set_neighbors(cells[1118],cells[805],cells[895],cells[879]);
	  set_offsets(cells[1046],0,2,0,0);
	  cells[1047]->set_vertices(vertices[221],vertices[220],vertices[213],vertices[173]);
	  cells[1047]->set_neighbors(cells[76],cells[1103],cells[1339],cells[1381]);
	  set_offsets(cells[1047],0,0,0,0);
	  cells[1048]->set_vertices(vertices[211],vertices[210],vertices[258],vertices[203]);
	  cells[1048]->set_neighbors(cells[1285],cells[391],cells[364],cells[1516]);
	  set_offsets(cells[1048],0,0,0,0);
	  cells[1049]->set_vertices(vertices[169],vertices[161],vertices[209],vertices[218]);
	  cells[1049]->set_neighbors(cells[979],cells[1307],cells[165],cells[1274]);
	  set_offsets(cells[1049],0,0,0,0);
	  cells[1050]->set_vertices(vertices[178],vertices[131],vertices[130],vertices[123]);
	  cells[1050]->set_neighbors(cells[819],cells[1051],cells[1054],cells[999]);
	  set_offsets(cells[1050],0,0,0,0);
	  cells[1051]->set_vertices(vertices[170],vertices[178],vertices[130],vertices[123]);
	  cells[1051]->set_neighbors(cells[1050],cells[854],cells[923],cells[408]);
	  set_offsets(cells[1051],0,0,0,0);
	  cells[1052]->set_vertices(vertices[240],vertices[233],vertices[192],vertices[280]);
	  cells[1052]->set_neighbors(cells[1642],cells[1678],cells[1681],cells[1443]);
	  set_offsets(cells[1052],2,0,2,0);
	  cells[1053]->set_vertices(vertices[185],vertices[144],vertices[192],vertices[232]);
	  cells[1053]->set_neighbors(cells[1259],cells[1406],cells[1086],cells[914]);
	  set_offsets(cells[1053],0,2,2,0);
	  cells[1054]->set_vertices(vertices[131],vertices[178],vertices[171],vertices[123]);
	  cells[1054]->set_neighbors(cells[923],cells[1094],cells[1050],cells[1097]);
	  set_offsets(cells[1054],0,0,0,0);
	  cells[1055]->set_vertices(vertices[171],vertices[170],vertices[163],vertices[123]);
	  cells[1055]->set_neighbors(cells[687],cells[917],cells[923],cells[1286]);
	  set_offsets(cells[1055],0,0,0,0);
	  cells[1056]->set_vertices(vertices[246],vertices[6],vertices[285],vertices[245]);
	  cells[1056]->set_neighbors(cells[1197],cells[1669],cells[197],cells[1340]);
	  set_offsets(cells[1056],2,6,0,2);
	  cells[1057]->set_vertices(vertices[270],vertices[213],vertices[262],vertices[222]);
	  cells[1057]->set_neighbors(cells[172],cells[1591],cells[1586],cells[1585]);
	  set_offsets(cells[1057],0,0,0,0);
	  cells[1058]->set_vertices(vertices[259],vertices[28],vertices[19],vertices[267]);
	  cells[1058]->set_neighbors(cells[1195],cells[1370],cells[1526],cells[1467]);
	  set_offsets(cells[1058],0,4,4,0);
	  cells[1059]->set_vertices(vertices[217],vertices[257],vertices[265],vertices[266]);
	  cells[1059]->set_neighbors(cells[570],cells[1367],cells[1514],cells[1608]);
	  set_offsets(cells[1059],0,0,0,0);
	  cells[1060]->set_vertices(vertices[173],vertices[165],vertices[213],vertices[222]);
	  cells[1060]->set_neighbors(cells[1338],cells[1103],cells[246],cells[76]);
	  set_offsets(cells[1060],0,0,0,0);
	  cells[1061]->set_vertices(vertices[188],vertices[139],vertices[179],vertices[131]);
	  cells[1061]->set_neighbors(cells[653],cells[896],cells[620],cells[897]);
	  set_offsets(cells[1061],0,0,0,0);
	  cells[1062]->set_vertices(vertices[220],vertices[163],vertices[211],vertices[212]);
	  cells[1062]->set_neighbors(cells[905],cells[832],cells[81],cells[1322]);
	  set_offsets(cells[1062],0,0,0,0);
	  cells[1063]->set_vertices(vertices[180],vertices[123],vertices[171],vertices[172]);
	  cells[1063]->set_neighbors(cells[917],cells[870],cells[660],cells[1094]);
	  set_offsets(cells[1063],0,0,0,0);
	  cells[1064]->set_vertices(vertices[180],vertices[133],vertices[132],vertices[125]);
	  cells[1064]->set_neighbors(cells[828],cells[711],cells[1016],cells[927]);
	  set_offsets(cells[1064],0,0,0,0);
	  cells[1065]->set_vertices(vertices[129],vertices[121],vertices[178],vertices[130]);
	  cells[1065]->set_neighbors(cells[408],cells[1017],cells[98],cells[1092]);
	  set_offsets(cells[1065],0,0,0,0);
	  cells[1066]->set_vertices(vertices[212],vertices[203],vertices[252],vertices[204]);
	  cells[1066]->set_neighbors(cells[1159],cells[1243],cells[1281],cells[1024]);
	  set_offsets(cells[1066],0,0,0,0);
	  cells[1067]->set_vertices(vertices[141],vertices[181],vertices[189],vertices[190]);
	  cells[1067]->set_neighbors(cells[921],cells[1145],cells[934],cells[926]);
	  set_offsets(cells[1067],0,0,0,0);
	  cells[1068]->set_vertices(vertices[131],vertices[140],vertices[132],vertices[180]);
	  cells[1068]->set_neighbors(cells[927],cells[1069],cells[621],cells[743]);
	  set_offsets(cells[1068],0,0,0,0);
	  cells[1069]->set_vertices(vertices[131],vertices[123],vertices[180],vertices[132]);
	  cells[1069]->set_neighbors(cells[660],cells[1068],cells[763],cells[1094]);
	  set_offsets(cells[1069],0,0,0,0);
	  cells[1070]->set_vertices(vertices[182],vertices[125],vertices[174],vertices[134]);
	  cells[1070]->set_neighbors(cells[736],cells[780],cells[590],cells[589]);
	  set_offsets(cells[1070],0,0,0,0);
	  cells[1071]->set_vertices(vertices[182],vertices[173],vertices[222],vertices[174]);
	  cells[1071]->set_neighbors(cells[246],cells[847],cells[589],cells[244]);
	  set_offsets(cells[1071],0,0,0,0);
	  cells[1072]->set_vertices(vertices[137],vertices[96],vertices[144],vertices[184]);
	  cells[1072]->set_neighbors(cells[1120],cells[1128],cells[685],cells[67]);
	  set_offsets(cells[1072],0,2,2,0);
	  cells[1073]->set_vertices(vertices[129],vertices[136],vertices[176],vertices[128]);
	  cells[1073]->set_neighbors(cells[782],cells[516],cells[357],cells[965]);
	  set_offsets(cells[1073],0,0,0,0);
	  cells[1074]->set_vertices(vertices[176],vertices[184],vertices[135],vertices[136]);
	  cells[1074]->set_neighbors(cells[816],cells[782],cells[965],cells[998]);
	  set_offsets(cells[1074],1,1,0,1);
	  cells[1075]->set_vertices(vertices[174],vertices[214],vertices[167],vertices[166]);
	  cells[1075]->set_neighbors(cells[1160],cells[24],cells[918],cells[1311]);
	  set_offsets(cells[1075],0,0,0,0);
	  cells[1076]->set_vertices(vertices[135],vertices[182],vertices[175],vertices[127]);
	  cells[1076]->set_neighbors(cells[422],cells[1032],cells[159],cells[996]);
	  set_offsets(cells[1076],0,0,0,0);
	  cells[1077]->set_vertices(vertices[175],vertices[174],vertices[167],vertices[127]);
	  cells[1077]->set_neighbors(cells[869],cells[1078],cells[422],cells[1310]);
	  set_offsets(cells[1077],0,0,0,0);
	  cells[1078]->set_vertices(vertices[175],vertices[168],vertices[127],vertices[167]);
	  cells[1078]->set_neighbors(cells[809],cells[1077],cells[1313],cells[977]);
	  set_offsets(cells[1078],0,1,0,0);
	  cells[1079]->set_vertices(vertices[192],vertices[145],vertices[200],vertices[152]);
	  cells[1079]->set_neighbors(cells[1102],cells[455],cells[686],cells[718]);
	  set_offsets(cells[1079],0,0,0,0);
	  cells[1080]->set_vertices(vertices[178],vertices[218],vertices[171],vertices[170]);
	  cells[1080]->set_neighbors(cells[1286],cells[923],cells[959],cells[1337]);
	  set_offsets(cells[1080],0,0,0,0);
	  cells[1081]->set_vertices(vertices[139],vertices[146],vertices[186],vertices[98]);
	  cells[1081]->set_neighbors(cells[1124],cells[413],cells[824],cells[1138]);
	  set_offsets(cells[1081],0,2,0,2);
	  cells[1082]->set_vertices(vertices[190],vertices[143],vertices[142],vertices[135]);
	  cells[1082]->set_neighbors(cells[902],cells[1083],cells[920],cells[1118]);
	  set_offsets(cells[1082],0,0,0,0);
	  cells[1083]->set_vertices(vertices[182],vertices[190],vertices[142],vertices[135]);
	  cells[1083]->set_neighbors(cells[1082],cells[1106],cells[760],cells[935]);
	  set_offsets(cells[1083],0,0,0,0);
	  cells[1084]->set_vertices(vertices[177],vertices[176],vertices[184],vertices[224]);
	  cells[1084]->set_neighbors(cells[1352],cells[913],cells[1358],cells[542]);
	  set_offsets(cells[1084],0,0,0,0);
	  cells[1085]->set_vertices(vertices[224],vertices[216],vertices[175],vertices[223]);
	  cells[1085]->set_neighbors(cells[837],cells[753],cells[1605],cells[1360]);
	  set_offsets(cells[1085],1,1,0,0);
	  cells[1086]->set_vertices(vertices[185],vertices[144],vertices[232],vertices[184]);
	  cells[1086]->set_neighbors(cells[1397],cells[1347],cells[1128],cells[1053]);
	  set_offsets(cells[1086],0,2,0,0);
	  cells[1087]->set_vertices(vertices[168],vertices[160],vertices[208],vertices[161]);
	  cells[1087]->set_neighbors(cells[1217],cells[216],cells[744],cells[1269]);
	  set_offsets(cells[1087],0,0,0,0);
	  cells[1088]->set_vertices(vertices[198],vertices[151],vertices[239],vertices[199]);
	  cells[1088]->set_neighbors(cells[1355],cells[1486],cells[1240],cells[1436]);
	  set_offsets(cells[1088],2,2,0,2);
	  cells[1089]->set_vertices(vertices[185],vertices[97],vertices[137],vertices[146]);
	  cells[1089]->set_neighbors(cells[915],cells[1125],cells[1091],cells[1129]);
	  set_offsets(cells[1089],0,2,0,2);
	  cells[1090]->set_vertices(vertices[185],vertices[137],vertices[184],vertices[177]);
	  cells[1090]->set_neighbors(cells[853],cells[1347],cells[877],cells[1128]);
	  set_offsets(cells[1090],0,0,0,0);
	  cells[1091]->set_vertices(vertices[145],vertices[97],vertices[185],vertices[146]);
	  cells[1091]->set_neighbors(cells[1089],cells[602],cells[457],cells[1130]);
	  set_offsets(cells[1091],2,2,0,2);
	  cells[1092]->set_vertices(vertices[129],vertices[121],vertices[169],vertices[178]);
	  cells[1092]->set_neighbors(cells[634],cells[646],cells[1065],cells[713]);
	  set_offsets(cells[1092],0,0,0,0);
	  cells[1093]->set_vertices(vertices[231],vertices[272],vertices[271],vertices[279]);
	  cells[1093]->set_neighbors(cells[1477],cells[936],cells[1292],cells[1299]);
	  set_offsets(cells[1093],0,1,0,0);
	  cells[1094]->set_vertices(vertices[131],vertices[123],vertices[171],vertices[180]);
	  cells[1094]->set_neighbors(cells[1063],cells[1095],cells[1069],cells[1054]);
	  set_offsets(cells[1094],0,0,0,0);
	  cells[1095]->set_vertices(vertices[131],vertices[179],vertices[180],vertices[171]);
	  cells[1095]->set_neighbors(cells[261],cells[1094],cells[1097],cells[896]);
	  set_offsets(cells[1095],0,0,0,0);
	  cells[1096]->set_vertices(vertices[188],vertices[139],vertices[140],vertices[100]);
	  cells[1096]->set_neighbors(cells[827],cells[1113],cells[1135],cells[620]);
	  set_offsets(cells[1096],0,0,0,2);
	  cells[1097]->set_vertices(vertices[178],vertices[179],vertices[131],vertices[171]);
	  cells[1097]->set_neighbors(cells[1095],cells[1054],cells[1332],cells[817]);
	  set_offsets(cells[1097],0,0,0,0);
	  cells[1098]->set_vertices(vertices[187],vertices[139],vertices[186],vertices[179]);
	  cells[1098]->set_neighbors(cells[653],cells[1171],cells[897],cells[1138]);
	  set_offsets(cells[1098],0,0,0,0);
	  cells[1099]->set_vertices(vertices[147],vertices[99],vertices[187],vertices[148]);
	  cells[1099]->set_neighbors(cells[786],cells[1133],cells[130],cells[1136]);
	  set_offsets(cells[1099],2,2,0,2);
	  cells[1100]->set_vertices(vertices[274],vertices[282],vertices[225],vertices[234]);
	  cells[1100]->set_neighbors(cells[344],cells[1655],cells[1650],cells[1692]);
	  set_offsets(cells[1100],0,0,0,0);
	  cells[1101]->set_vertices(vertices[0],vertices[8],vertices[241],vertices[1]);
	  cells[1101]->set_neighbors(cells[1447],cells[1402],cells[295],cells[975]);
	  set_offsets(cells[1101],4,4,0,4);
	  cells[1102]->set_vertices(vertices[200],vertices[145],vertices[153],vertices[152]);
	  cells[1102]->set_neighbors(cells[951],cells[1112],cells[1079],cells[1142]);
	  set_offsets(cells[1102],0,0,0,0);
	  cells[1103]->set_vertices(vertices[221],vertices[173],vertices[213],vertices[222]);
	  cells[1103]->set_neighbors(cells[1060],cells[1586],cells[1116],cells[1047]);
	  set_offsets(cells[1103],0,0,0,0);
	  cells[1104]->set_vertices(vertices[101],vertices[150],vertices[141],vertices[189]);
	  cells[1104]->set_neighbors(cells[1145],cells[390],cells[1147],cells[937]);
	  set_offsets(cells[1104],2,2,0,0);
	  cells[1105]->set_vertices(vertices[182],vertices[142],vertices[133],vertices[134]);
	  cells[1105]->set_neighbors(cells[389],cells[590],cells[1106],cells[935]);
	  set_offsets(cells[1105],0,0,0,0);
	  cells[1106]->set_vertices(vertices[135],vertices[142],vertices[182],vertices[134]);
	  cells[1106]->set_neighbors(cells[1105],cells[159],cells[140],cells[1083]);
	  set_offsets(cells[1106],0,0,0,0);
	  cells[1107]->set_vertices(vertices[180],vertices[181],vertices[133],vertices[173]);
	  cells[1107]->set_neighbors(cells[986],cells[1016],cells[1342],cells[490]);
	  set_offsets(cells[1107],0,0,0,0);
	  cells[1108]->set_vertices(vertices[247],vertices[248],vertices[255],vertices[207]);
	  cells[1108]->set_neighbors(cells[1504],cells[1291],cells[1441],cells[1534]);
	  set_offsets(cells[1108],0,1,0,0);
	  cells[1109]->set_vertices(vertices[181],vertices[173],vertices[221],vertices[230]);
	  cells[1109]->set_neighbors(cells[1116],cells[863],cells[194],cells[1343]);
	  set_offsets(cells[1109],0,0,0,0);
	  cells[1110]->set_vertices(vertices[266],vertices[209],vertices[258],vertices[218]);
	  cells[1110]->set_neighbors(cells[1517],cells[1335],cells[1369],cells[1559]);
	  set_offsets(cells[1110],0,0,0,0);
	  cells[1111]->set_vertices(vertices[159],vertices[206],vertices[151],vertices[158]);
	  cells[1111]->set_neighbors(cells[146],cells[585],cells[1031],cells[1189]);
	  set_offsets(cells[1111],0,0,0,0);
	  cells[1112]->set_vertices(vertices[160],vertices[152],vertices[200],vertices[153]);
	  cells[1112]->set_neighbors(cells[1102],cells[997],cells[43],cells[1214]);
	  set_offsets(cells[1112],0,0,0,0);
	  cells[1113]->set_vertices(vertices[141],vertices[188],vertices[140],vertices[100]);
	  cells[1113]->set_neighbors(cells[1096],cells[883],cells[332],cells[695]);
	  set_offsets(cells[1113],0,0,0,2);
	  cells[1114]->set_vertices(vertices[236],vertices[189],vertices[188],vertices[181]);
	  cells[1114]->set_neighbors(cells[926],cells[327],cells[830],cells[1383]);
	  set_offsets(cells[1114],0,0,0,0);
	  cells[1115]->set_vertices(vertices[194],vertices[146],vertices[234],vertices[187]);
	  cells[1115]->set_neighbors(cells[245],cells[1419],cells[840],cells[1414]);
	  set_offsets(cells[1115],2,2,0,0);
	  cells[1116]->set_vertices(vertices[230],vertices[173],vertices[221],vertices[222]);
	  cells[1116]->set_neighbors(cells[1103],cells[673],cells[244],cells[1109]);
	  set_offsets(cells[1116],0,0,0,0);
	  cells[1117]->set_vertices(vertices[176],vertices[183],vertices[175],vertices[135]);
	  cells[1117]->set_neighbors(cells[996],cells[1032],cells[998],cells[1359]);
	  set_offsets(cells[1117],1,0,0,0);
	  cells[1118]->set_vertices(vertices[190],vertices[102],vertices[142],vertices[143]);
	  cells[1118]->set_neighbors(cells[901],cells[1082],cells[757],cells[1046]);
	  set_offsets(cells[1118],0,2,0,0);
	  cells[1119]->set_vertices(vertices[230],vertices[238],vertices[278],vertices[229]);
	  cells[1119]->set_neighbors(cells[1673],cells[1166],cells[1168],cells[1675]);
	  set_offsets(cells[1119],0,0,0,0);
	  cells[1120]->set_vertices(vertices[144],vertices[96],vertices[143],vertices[184]);
	  cells[1120]->set_neighbors(cells[316],cells[1153],cells[1072],cells[489]);
	  set_offsets(cells[1120],3,3,0,1);
	  cells[1121]->set_vertices(vertices[26],vertices[266],vertices[18],vertices[259]);
	  cells[1121]->set_neighbors(cells[214],cells[1569],cells[1176],cells[1560]);
	  set_offsets(cells[1121],4,0,4,0);
	  cells[1122]->set_vertices(vertices[179],vertices[171],vertices[219],vertices[228]);
	  cells[1122]->set_neighbors(cells[962],cells[1201],cells[261],cells[1333]);
	  set_offsets(cells[1122],0,0,0,0);
	  cells[1123]->set_vertices(vertices[133],vertices[180],vertices[140],vertices[188]);
	  cells[1123]->set_neighbors(cells[621],cells[695],cells[490],cells[927]);
	  set_offsets(cells[1123],0,0,0,0);
	  cells[1124]->set_vertices(vertices[146],vertices[137],vertices[186],vertices[98]);
	  cells[1124]->set_neighbors(cells[688],cells[1081],cells[915],cells[1125]);
	  set_offsets(cells[1124],2,0,0,2);
	  cells[1125]->set_vertices(vertices[146],vertices[137],vertices[185],vertices[186]);
	  cells[1125]->set_neighbors(cells[877],cells[65],cells[1124],cells[1089]);
	  set_offsets(cells[1125],2,0,0,0);
	  cells[1126]->set_vertices(vertices[131],vertices[186],vertices[138],vertices[139]);
	  cells[1126]->set_neighbors(cells[413],cells[376],cells[653],cells[1127]);
	  set_offsets(cells[1126],0,0,0,0);
	  cells[1127]->set_vertices(vertices[131],vertices[178],vertices[138],vertices[186]);
	  cells[1127]->set_neighbors(cells[656],cells[1126],cells[817],cells[999]);
	  set_offsets(cells[1127],0,0,0,0);
	  cells[1128]->set_vertices(vertices[144],vertices[137],vertices[184],vertices[185]);
	  cells[1128]->set_neighbors(cells[1090],cells[1086],cells[1129],cells[1072]);
	  set_offsets(cells[1128],2,0,0,0);
	  cells[1129]->set_vertices(vertices[144],vertices[97],vertices[137],vertices[185]);
	  cells[1129]->set_neighbors(cells[1089],cells[1128],cells[1130],cells[67]);
	  set_offsets(cells[1129],2,2,0,0);
	  cells[1130]->set_vertices(vertices[144],vertices[97],vertices[185],vertices[145]);
	  cells[1130]->set_neighbors(cells[1091],cells[914],cells[614],cells[1129]);
	  set_offsets(cells[1130],2,2,0,2);
	  cells[1131]->set_vertices(vertices[226],vertices[217],vertices[266],vertices[218]);
	  cells[1131]->set_neighbors(cells[1369],cells[1373],cells[873],cells[1178]);
	  set_offsets(cells[1131],0,0,0,0);
	  cells[1132]->set_vertices(vertices[194],vertices[242],vertices[202],vertices[195]);
	  cells[1132]->set_neighbors(cells[1513],cells[1227],cells[1002],cells[1454]);
	  set_offsets(cells[1132],0,0,0,0);
	  cells[1133]->set_vertices(vertices[187],vertices[147],vertices[148],vertices[196]);
	  cells[1133]->set_neighbors(cells[1192],cells[768],cells[1421],cells[1099]);
	  set_offsets(cells[1133],0,2,2,2);
	  cells[1134]->set_vertices(vertices[148],vertices[139],vertices[187],vertices[188]);
	  cells[1134]->set_neighbors(cells[897],cells[1327],cells[1135],cells[786]);
	  set_offsets(cells[1134],2,0,0,0);
	  cells[1135]->set_vertices(vertices[148],vertices[139],vertices[188],vertices[100]);
	  cells[1135]->set_neighbors(cells[1096],cells[332],cells[648],cells[1134]);
	  set_offsets(cells[1135],2,0,0,2);
	  cells[1136]->set_vertices(vertices[146],vertices[99],vertices[187],vertices[147]);
	  cells[1136]->set_neighbors(cells[1099],cells[840],cells[461],cells[1137]);
	  set_offsets(cells[1136],2,2,0,2);
	  cells[1137]->set_vertices(vertices[146],vertices[99],vertices[139],vertices[187]);
	  cells[1137]->set_neighbors(cells[786],cells[1138],cells[1136],cells[824]);
	  set_offsets(cells[1137],2,2,0,0);
	  cells[1138]->set_vertices(vertices[146],vertices[139],vertices[186],vertices[187]);
	  cells[1138]->set_neighbors(cells[1098],cells[245],cells[1137],cells[1081]);
	  set_offsets(cells[1138],2,0,0,0);
	  cells[1139]->set_vertices(vertices[43],vertices[283],vertices[42],vertices[35]);
	  cells[1139]->set_neighbors(cells[1341],cells[313],cells[1004],cells[931]);
	  set_offsets(cells[1139],4,0,4,4);
	  cells[1140]->set_vertices(vertices[229],vertices[276],vertices[269],vertices[221]);
	  cells[1140]->set_neighbors(cells[1460],cells[1670],cells[1664],cells[1638]);
	  set_offsets(cells[1140],0,0,0,0);
	  cells[1141]->set_vertices(vertices[208],vertices[200],vertices[159],vertices[207]);
	  cells[1141]->set_neighbors(cells[137],cells[117],cells[1451],cells[1267]);
	  set_offsets(cells[1141],1,1,0,0);
	  cells[1142]->set_vertices(vertices[200],vertices[145],vertices[193],vertices[153]);
	  cells[1142]->set_neighbors(cells[1149],cells[1229],cells[1102],cells[718]);
	  set_offsets(cells[1142],0,0,0,0);
	  cells[1143]->set_vertices(vertices[240],vertices[192],vertices[193],vertices[200]);
	  cells[1143]->set_neighbors(cells[718],cells[1432],cells[134],cells[1443]);
	  set_offsets(cells[1143],0,0,0,0);
	  cells[1144]->set_vertices(vertices[30],vertices[261],vertices[22],vertices[270]);
	  cells[1144]->set_neighbors(cells[1533],cells[1202],cells[1390],cells[1535]);
	  set_offsets(cells[1144],4,0,4,0);
	  cells[1145]->set_vertices(vertices[141],vertices[190],vertices[189],vertices[150]);
	  cells[1145]->set_neighbors(cells[1167],cells[1104],cells[895],cells[1067]);
	  set_offsets(cells[1145],0,0,0,2);
	  cells[1146]->set_vertices(vertices[148],vertices[101],vertices[189],vertices[149]);
	  cells[1146]->set_neighbors(cells[1147],cells[862],cells[417],cells[390]);
	  set_offsets(cells[1146],2,2,0,2);
	  cells[1147]->set_vertices(vertices[101],vertices[150],vertices[189],vertices[149]);
	  cells[1147]->set_neighbors(cells[1216],cells[1146],cells[942],cells[1104]);
	  set_offsets(cells[1147],2,2,0,2);
	  cells[1148]->set_vertices(vertices[237],vertices[246],vertices[286],vertices[285]);
	  cells[1148]->set_neighbors(cells[1208],cells[1635],cells[1668],cells[1478]);
	  set_offsets(cells[1148],0,2,0,0);
	  cells[1149]->set_vertices(vertices[193],vertices[202],vertices[153],vertices[145]);
	  cells[1149]->set_neighbors(cells[1151],cells[1142],cells[1000],cells[48]);
	  set_offsets(cells[1149],0,0,0,0);
	  cells[1150]->set_vertices(vertices[194],vertices[154],vertices[202],vertices[145]);
	  cells[1150]->set_neighbors(cells[1151],cells[1000],cells[719],cells[1191]);
	  set_offsets(cells[1150],0,0,0,0);
	  cells[1151]->set_vertices(vertices[202],vertices[154],vertices[153],vertices[145]);
	  cells[1151]->set_neighbors(cells[73],cells[1149],cells[1150],cells[783]);
	  set_offsets(cells[1151],0,0,0,0);
	  cells[1152]->set_vertices(vertices[274],vertices[266],vertices[267],vertices[219]);
	  cells[1152]->set_neighbors(cells[1622],cells[1440],cells[529],cells[1620]);
	  set_offsets(cells[1152],0,0,0,0);
	  cells[1153]->set_vertices(vertices[144],vertices[143],vertices[191],vertices[184]);
	  cells[1153]->set_neighbors(cells[1154],cells[1397],cells[1120],cells[892]);
	  set_offsets(cells[1153],3,0,0,1);
	  cells[1154]->set_vertices(vertices[191],vertices[143],vertices[183],vertices[184]);
	  cells[1154]->set_neighbors(cells[904],cells[1309],cells[1153],cells[1155]);
	  set_offsets(cells[1154],0,0,0,1);
	  cells[1155]->set_vertices(vertices[183],vertices[143],vertices[191],vertices[190]);
	  cells[1155]->set_neighbors(cells[1156],cells[325],cells[920],cells[1154]);
	  set_offsets(cells[1155],0,0,0,0);
	  cells[1156]->set_vertices(vertices[191],vertices[143],vertices[150],vertices[190]);
	  cells[1156]->set_neighbors(cells[757],cells[973],cells[1155],cells[1158]);
	  set_offsets(cells[1156],0,0,2,0);
	  cells[1157]->set_vertices(vertices[150],vertices[103],vertices[191],vertices[151]);
	  cells[1157]->set_neighbors(cells[1033],cells[1215],cells[670],cells[1158]);
	  set_offsets(cells[1157],2,2,0,2);
	  cells[1158]->set_vertices(vertices[150],vertices[103],vertices[143],vertices[191]);
	  cells[1158]->set_neighbors(cells[892],cells[1156],cells[1157],cells[938]);
	  set_offsets(cells[1158],2,2,0,0);
	  cells[1159]->set_vertices(vertices[252],vertices[195],vertices[204],vertices[203]);
	  cells[1159]->set_neighbors(cells[1011],cells[1066],cells[1256],cells[1377]);
	  set_offsets(cells[1159],0,0,0,0);
	  cells[1160]->set_vertices(vertices[167],vertices[214],vertices[159],vertices[166]);
	  cells[1160]->set_neighbors(cells[1264],cells[1037],cells[1075],cells[1161]);
	  set_offsets(cells[1160],0,0,0,0);
	  cells[1161]->set_vertices(vertices[167],vertices[214],vertices[207],vertices[159]);
	  cells[1161]->set_neighbors(cells[1260],cells[117],cells[1160],cells[1203]);
	  set_offsets(cells[1161],0,0,0,0);
	  cells[1162]->set_vertices(vertices[286],vertices[287],vertices[46],vertices[246]);
	  cells[1162]->set_neighbors(cells[1400],cells[1208],cells[1634],cells[1726]);
	  set_offsets(cells[1162],0,0,4,2);
	  cells[1163]->set_vertices(vertices[238],vertices[189],vertices[229],vertices[181]);
	  cells[1163]->set_neighbors(cells[830],cells[1168],cells[921],cells[347]);
	  set_offsets(cells[1163],0,0,0,0);
	  cells[1164]->set_vertices(vertices[190],vertices[230],vertices[183],vertices[182]);
	  cells[1164]->set_neighbors(cells[1356],cells[760],cells[1165],cells[419]);
	  set_offsets(cells[1164],0,0,0,0);
	  cells[1165]->set_vertices(vertices[190],vertices[181],vertices[230],vertices[182]);
	  cells[1165]->set_neighbors(cells[194],cells[1164],cells[694],cells[348]);
	  set_offsets(cells[1165],0,0,0,0);
	  cells[1166]->set_vertices(vertices[229],vertices[221],vertices[278],vertices[230]);
	  cells[1166]->set_neighbors(cells[1303],cells[1119],cells[863],cells[1670]);
	  set_offsets(cells[1166],0,0,0,0);
	  cells[1167]->set_vertices(vertices[238],vertices[150],vertices[189],vertices[190]);
	  cells[1167]->set_neighbors(cells[1145],cells[921],cells[973],cells[1430]);
	  set_offsets(cells[1167],0,2,0,0);
	  cells[1168]->set_vertices(vertices[230],vertices[238],vertices[229],vertices[181]);
	  cells[1168]->set_neighbors(cells[1163],cells[863],cells[348],cells[1119]);
	  set_offsets(cells[1168],0,0,0,0);
	  cells[1169]->set_vertices(vertices[263],vertices[256],vertices[255],vertices[16]);
	  cells[1169]->set_neighbors(cells[1266],cells[1491],cells[1316],cells[1582]);
	  set_offsets(cells[1169],0,1,0,5);
	  cells[1170]->set_vertices(vertices[272],vertices[224],vertices[223],vertices[264]);
	  cells[1170]->set_neighbors(cells[1605],cells[1348],cells[1315],cells[995]);
	  set_offsets(cells[1170],1,1,0,1);
	  cells[1171]->set_vertices(vertices[234],vertices[187],vertices[186],vertices[179]);
	  cells[1171]->set_neighbors(cells[1098],cells[1384],cells[307],cells[245]);
	  set_offsets(cells[1171],0,0,0,0);
	  cells[1172]->set_vertices(vertices[285],vertices[46],vertices[277],vertices[37]);
	  cells[1172]->set_neighbors(cells[1474],cells[58],cells[1219],cells[1672]);
	  set_offsets(cells[1172],0,4,0,4);
	  cells[1173]->set_vertices(vertices[185],vertices[186],vertices[177],vertices[234]);
	  cells[1173]->set_neighbors(cells[1185],cells[1416],cells[65],cells[877]);
	  set_offsets(cells[1173],0,0,0,0);
	  cells[1174]->set_vertices(vertices[225],vertices[177],vertices[224],vertices[217]);
	  cells[1174]->set_neighbors(cells[1321],cells[277],cells[1187],cells[1396]);
	  set_offsets(cells[1174],0,0,0,0);
	  cells[1175]->set_vertices(vertices[4],vertices[252],vertices[245],vertices[12]);
	  cells[1175]->set_neighbors(cells[1484],cells[966],cells[1379],cells[1523]);
	  set_offsets(cells[1175],4,0,0,4);
	  cells[1176]->set_vertices(vertices[267],vertices[266],vertices[26],vertices[259]);
	  cells[1176]->set_neighbors(cells[1121],cells[1370],cells[1622],cells[1620]);
	  set_offsets(cells[1176],0,0,4,0);
	  cells[1177]->set_vertices(vertices[14],vertices[245],vertices[254],vertices[253]);
	  cells[1177]->set_neighbors(cells[1305],cells[1224],cells[1530],cells[1466]);
	  set_offsets(cells[1177],4,0,0,0);
	  cells[1178]->set_vertices(vertices[274],vertices[217],vertices[266],vertices[226]);
	  cells[1178]->set_neighbors(cells[1131],cells[529],cells[85],cells[1367]);
	  set_offsets(cells[1178],0,0,0,0);
	  cells[1179]->set_vertices(vertices[162],vertices[202],vertices[155],vertices[154]);
	  cells[1179]->set_neighbors(cells[100],cells[770],cells[783],cells[250]);
	  set_offsets(cells[1179],0,0,0,0);
	  cells[1180]->set_vertices(vertices[40],vertices[287],vertices[279],vertices[280]);
	  cells[1180]->set_neighbors(cells[957],cells[1541],cells[1405],cells[1594]);
	  set_offsets(cells[1180],5,0,0,1);
	  cells[1181]->set_vertices(vertices[269],vertices[277],vertices[38],vertices[29]);
	  cells[1181]->set_neighbors(cells[1246],cells[1583],cells[1230],cells[322]);
	  set_offsets(cells[1181],0,0,4,4);
	  cells[1182]->set_vertices(vertices[237],vertices[149],vertices[198],vertices[197]);
	  cells[1182]->set_neighbors(cells[77],cells[1206],cells[1199],cells[1428]);
	  set_offsets(cells[1182],0,2,2,2);
	  cells[1183]->set_vertices(vertices[41],vertices[281],vertices[33],vertices[42]);
	  cells[1183]->set_neighbors(cells[628],cells[102],cells[1606],cells[1653]);
	  set_offsets(cells[1183],4,0,4,4);
	  cells[1184]->set_vertices(vertices[229],vertices[228],vertices[181],vertices[236]);
	  cells[1184]->set_neighbors(cells[327],cells[830],cells[714],cells[1375]);
	  set_offsets(cells[1184],0,0,0,0);
	  cells[1185]->set_vertices(vertices[234],vertices[186],vertices[177],vertices[226]);
	  cells[1185]->set_neighbors(cells[1039],cells[889],cells[1384],cells[1173]);
	  set_offsets(cells[1185],0,0,0,0);
	  cells[1186]->set_vertices(vertices[177],vertices[169],vertices[217],vertices[226]);
	  cells[1186]->set_neighbors(cells[873],cells[1187],cells[47],cells[1321]);
	  set_offsets(cells[1186],0,0,0,0);
	  cells[1187]->set_vertices(vertices[226],vertices[225],vertices[217],vertices[177]);
	  cells[1187]->set_neighbors(cells[1174],cells[1186],cells[889],cells[85]);
	  set_offsets(cells[1187],0,0,0,0);
	  cells[1188]->set_vertices(vertices[261],vertices[260],vertices[20],vertices[253]);
	  cells[1188]->set_neighbors(cells[1293],cells[1572],cells[1344],cells[217]);
	  set_offsets(cells[1188],0,0,4,0);
	  cells[1189]->set_vertices(vertices[199],vertices[206],vertices[151],vertices[159]);
	  cells[1189]->set_neighbors(cells[1111],cells[1222],cells[1205],cells[1240]);
	  set_offsets(cells[1189],0,0,0,0);
	  cells[1190]->set_vertices(vertices[157],vertices[204],vertices[149],vertices[156]);
	  cells[1190]->set_neighbors(cells[37],cells[971],cells[1014],cells[35]);
	  set_offsets(cells[1190],0,0,0,0);
	  cells[1191]->set_vertices(vertices[202],vertices[194],vertices[147],vertices[154]);
	  cells[1191]->set_neighbors(cells[922],cells[100],cells[1150],cells[1227]);
	  set_offsets(cells[1191],0,0,0,0);
	  cells[1192]->set_vertices(vertices[147],vertices[148],vertices[196],vertices[156]);
	  cells[1192]->set_neighbors(cells[852],cells[976],cells[130],cells[1133]);
	  set_offsets(cells[1192],0,0,0,0);
	  cells[1193]->set_vertices(vertices[44],vertices[277],vertices[36],vertices[37]);
	  cells[1193]->set_neighbors(cells[1629],cells[368],cells[58],cells[1198]);
	  set_offsets(cells[1193],4,0,4,4);
	  cells[1194]->set_vertices(vertices[276],vertices[284],vertices[227],vertices[236]);
	  cells[1194]->set_neighbors(cells[1459],cells[1362],cells[1665],cells[1703]);
	  set_offsets(cells[1194],0,0,0,0);
	  cells[1195]->set_vertices(vertices[267],vertices[28],vertices[19],vertices[27]);
	  cells[1195]->set_neighbors(cells[339],cells[1515],cells[1574],cells[1058]);
	  set_offsets(cells[1195],0,4,4,4);
	  cells[1196]->set_vertices(vertices[235],vertices[187],vertices[236],vertices[196]);
	  cells[1196]->set_neighbors(cells[768],cells[1462],cells[1421],cells[1005]);
	  set_offsets(cells[1196],0,0,0,2);
	  cells[1197]->set_vertices(vertices[245],vertices[6],vertices[285],vertices[45]);
	  cells[1197]->set_neighbors(cells[1628],cells[1207],cells[1425],cells[1056]);
	  set_offsets(cells[1197],2,6,0,4);
	  cells[1198]->set_vertices(vertices[44],vertices[284],vertices[36],vertices[277]);
	  cells[1198]->set_neighbors(cells[1616],cells[1193],cells[1710],cells[1705]);
	  set_offsets(cells[1198],4,0,4,0);
	  cells[1199]->set_vertices(vertices[196],vertices[149],vertices[237],vertices[197]);
	  cells[1199]->set_neighbors(cells[1182],cells[1426],cells[168],cells[1429]);
	  set_offsets(cells[1199],2,2,0,2);
	  cells[1200]->set_vertices(vertices[227],vertices[187],vertices[179],vertices[236]);
	  cells[1200]->set_neighbors(cells[1392],cells[269],cells[1005],cells[307]);
	  set_offsets(cells[1200],0,0,0,0);
	  cells[1201]->set_vertices(vertices[228],vertices[227],vertices[219],vertices[179]);
	  cells[1201]->set_neighbors(cells[1372],cells[1122],cells[269],cells[326]);
	  set_offsets(cells[1201],0,0,0,0);
	  cells[1202]->set_vertices(vertices[30],vertices[270],vertices[22],vertices[263]);
	  cells[1202]->set_neighbors(cells[264],cells[974],cells[1388],cells[1144]);
	  set_offsets(cells[1202],4,0,4,0);
	  cells[1203]->set_vertices(vertices[215],vertices[214],vertices[207],vertices[167]);
	  cells[1203]->set_neighbors(cells[1161],cells[1204],cells[1306],cells[1350]);
	  set_offsets(cells[1203],0,0,0,0);
	  cells[1204]->set_vertices(vertices[215],vertices[208],vertices[167],vertices[207]);
	  cells[1204]->set_neighbors(cells[117],cells[1203],cells[1546],cells[919]);
	  set_offsets(cells[1204],0,1,0,0);
	  cells[1205]->set_vertices(vertices[199],vertices[206],vertices[159],vertices[207]);
	  cells[1205]->set_neighbors(cells[1260],cells[137],cells[1249],cells[1189]);
	  set_offsets(cells[1205],0,0,0,0);
	  cells[1206]->set_vertices(vertices[246],vertices[197],vertices[237],vertices[198]);
	  cells[1206]->set_neighbors(cells[1182],cells[1478],cells[1213],cells[1668]);
	  set_offsets(cells[1206],2,2,0,2);
	  cells[1207]->set_vertices(vertices[4],vertices[245],vertices[285],vertices[45]);
	  cells[1207]->set_neighbors(cells[1197],cells[1471],cells[236],cells[1630]);
	  set_offsets(cells[1207],6,2,0,4);
	  cells[1208]->set_vertices(vertices[286],vertices[246],vertices[46],vertices[285]);
	  cells[1208]->set_neighbors(cells[1340],cells[1672],cells[1148],cells[1162]);
	  set_offsets(cells[1208],0,2,4,0);
	  cells[1209]->set_vertices(vertices[21],vertices[30],vertices[269],vertices[261]);
	  cells[1209]->set_neighbors(cells[1390],cells[136],cells[1535],cells[1584]);
	  set_offsets(cells[1209],4,4,0,0);
	  cells[1210]->set_vertices(vertices[208],vertices[200],vertices[248],vertices[201]);
	  cells[1210]->set_neighbors(cells[1481],cells[940],cells[967],cells[1451]);
	  set_offsets(cells[1210],0,0,0,0);
	  cells[1211]->set_vertices(vertices[214],vertices[205],vertices[254],vertices[206]);
	  cells[1211]->set_neighbors(cells[1225],cells[1493],cells[1294],cells[1483]);
	  set_offsets(cells[1211],0,0,0,0);
	  cells[1212]->set_vertices(vertices[246],vertices[197],vertices[206],vertices[254]);
	  cells[1212]->set_neighbors(cells[1225],cells[1300],cells[1026],cells[1213]);
	  set_offsets(cells[1212],0,0,0,0);
	  cells[1213]->set_vertices(vertices[246],vertices[197],vertices[198],vertices[206]);
	  cells[1213]->set_neighbors(cells[77],cells[1028],cells[1212],cells[1206]);
	  set_offsets(cells[1213],0,0,0,0);
	  cells[1214]->set_vertices(vertices[160],vertices[152],vertices[159],vertices[200]);
	  cells[1214]->set_neighbors(cells[941],cells[1267],cells[1112],cells[767]);
	  set_offsets(cells[1214],1,1,0,1);
	  cells[1215]->set_vertices(vertices[191],vertices[198],vertices[150],vertices[151]);
	  cells[1215]->set_neighbors(cells[641],cells[1157],cells[1436],cells[1431]);
	  set_offsets(cells[1215],0,2,2,2);
	  cells[1216]->set_vertices(vertices[189],vertices[149],vertices[150],vertices[198]);
	  cells[1216]->set_neighbors(cells[643],cells[1430],cells[1428],cells[1147]);
	  set_offsets(cells[1216],0,2,2,2);
	  cells[1217]->set_vertices(vertices[160],vertices[153],vertices[208],vertices[161]);
	  cells[1217]->set_neighbors(cells[25],cells[1087],cells[701],cells[997]);
	  set_offsets(cells[1217],0,0,0,0);
	  cells[1218]->set_vertices(vertices[264],vertices[209],vertices[216],vertices[256]);
	  cells[1218]->set_neighbors(cells[1226],cells[1601],cells[1599],cells[1446]);
	  set_offsets(cells[1218],0,0,0,0);
	  cells[1219]->set_vertices(vertices[45],vertices[46],vertices[285],vertices[37]);
	  cells[1219]->set_neighbors(cells[1172],cells[1354],cells[219],cells[1628]);
	  set_offsets(cells[1219],4,4,0,4);
	  cells[1220]->set_vertices(vertices[39],vertices[279],vertices[32],vertices[40]);
	  cells[1220]->set_neighbors(cells[1541],cells[539],cells[1594],cells[0]);
	  set_offsets(cells[1220],4,0,5,5);
	  cells[1221]->set_vertices(vertices[151],vertices[192],vertices[199],vertices[200]);
	  cells[1221]->set_neighbors(cells[134],cells[1222],cells[455],cells[1355]);
	  set_offsets(cells[1221],0,1,0,1);
	  cells[1222]->set_vertices(vertices[151],vertices[200],vertices[199],vertices[159]);
	  cells[1222]->set_neighbors(cells[137],cells[1189],cells[941],cells[1221]);
	  set_offsets(cells[1222],0,1,0,0);
	  cells[1223]->set_vertices(vertices[2],vertices[250],vertices[243],vertices[10]);
	  cells[1223]->set_neighbors(cells[1007],cells[1423],cells[1387],cells[1009]);
	  set_offsets(cells[1223],4,0,0,4);
	  cells[1224]->set_vertices(vertices[14],vertices[253],vertices[254],vertices[262]);
	  cells[1224]->set_neighbors(cells[1482],cells[1234],cells[222],cells[1177]);
	  set_offsets(cells[1224],4,0,0,0);
	  cells[1225]->set_vertices(vertices[254],vertices[197],vertices[206],vertices[205]);
	  cells[1225]->set_neighbors(cells[460],cells[1211],cells[1265],cells[1212]);
	  set_offsets(cells[1225],0,0,0,0);
	  cells[1226]->set_vertices(vertices[216],vertices[208],vertices[256],vertices[209]);
	  cells[1226]->set_neighbors(cells[176],cells[1218],cells[1273],cells[1545]);
	  set_offsets(cells[1226],0,0,0,0);
	  cells[1227]->set_vertices(vertices[147],vertices[194],vertices[202],vertices[195]);
	  cells[1227]->set_neighbors(cells[1132],cells[1034],cells[864],cells[1191]);
	  set_offsets(cells[1227],0,0,0,0);
	  cells[1228]->set_vertices(vertices[2],vertices[243],vertices[43],vertices[3]);
	  cells[1228]->set_neighbors(cells[890],cells[374],cells[1423],cells[1445]);
	  set_offsets(cells[1228],6,2,4,6);
	  cells[1229]->set_vertices(vertices[200],vertices[153],vertices[193],vertices[201]);
	  cells[1229]->set_neighbors(cells[48],cells[1481],cells[967],cells[1142]);
	  set_offsets(cells[1229],0,0,0,0);
	  cells[1230]->set_vertices(vertices[36],vertices[277],vertices[269],vertices[29]);
	  cells[1230]->set_neighbors(cells[1181],cells[1527],cells[1629],cells[1639]);
	  set_offsets(cells[1230],4,0,0,4);
	  cells[1231]->set_vertices(vertices[241],vertices[281],vertices[41],vertices[2]);
	  cells[1231]->set_neighbors(cells[1606],cells[1317],cells[1540],cells[1685]);
	  set_offsets(cells[1231],2,0,4,6);
	  cells[1232]->set_vertices(vertices[242],vertices[42],vertices[281],vertices[2]);
	  cells[1232]->set_neighbors(cells[1606],cells[1540],cells[1699],cells[291]);
	  set_offsets(cells[1232],2,4,0,6);
	  cells[1233]->set_vertices(vertices[24],vertices[264],vertices[263],vertices[16]);
	  cells[1233]->set_neighbors(cells[1316],cells[939],cells[1324],cells[1248]);
	  set_offsets(cells[1233],5,1,0,5);
	  cells[1234]->set_vertices(vertices[14],vertices[254],vertices[255],vertices[262]);
	  cells[1234]->set_neighbors(cells[1538],cells[875],cells[1224],cells[944]);
	  set_offsets(cells[1234],4,0,0,0);
	  cells[1235]->set_vertices(vertices[23],vertices[263],vertices[15],vertices[16]);
	  cells[1235]->set_neighbors(cells[1491],cells[453],cells[939],cells[1525]);
	  set_offsets(cells[1235],4,0,4,5);
	  cells[1236]->set_vertices(vertices[210],vertices[202],vertices[203],vertices[155]);
	  cells[1236]->set_neighbors(cells[1237],cells[968],cells[250],cells[947]);
	  set_offsets(cells[1236],0,0,0,0);
	  cells[1237]->set_vertices(vertices[155],vertices[202],vertices[203],vertices[195]);
	  cells[1237]->set_neighbors(cells[665],cells[1011],cells[1034],cells[1236]);
	  set_offsets(cells[1237],0,0,0,0);
	  cells[1238]->set_vertices(vertices[31],vertices[38],vertices[271],vertices[279]);
	  cells[1238]->set_neighbors(cells[1245],cells[1489],cells[1015],cells[1395]);
	  set_offsets(cells[1238],4,4,0,0);
	  cells[1239]->set_vertices(vertices[6],vertices[254],vertices[247],vertices[14]);
	  cells[1239]->set_neighbors(cells[944],cells[266],cells[1466],cells[1485]);
	  set_offsets(cells[1239],4,0,0,4);
	  cells[1240]->set_vertices(vertices[199],vertices[198],vertices[151],vertices[206]);
	  cells[1240]->set_neighbors(cells[146],cells[1189],cells[1028],cells[1088]);
	  set_offsets(cells[1240],0,0,0,0);
	  cells[1241]->set_vertices(vertices[213],vertices[212],vertices[205],vertices[165]);
	  cells[1241]->set_neighbors(cells[1255],cells[1250],cells[1295],cells[964]);
	  set_offsets(cells[1241],0,0,0,0);
	  cells[1242]->set_vertices(vertices[195],vertices[250],vertices[203],vertices[243]);
	  cells[1242]->set_neighbors(cells[1519],cells[1256],cells[280],cells[665]);
	  set_offsets(cells[1242],0,0,0,0);
	  cells[1243]->set_vertices(vertices[212],vertices[252],vertices[205],vertices[204]);
	  cells[1243]->set_neighbors(cells[1470],cells[265],cells[1066],cells[95]);
	  set_offsets(cells[1243],0,0,0,0);
	  cells[1244]->set_vertices(vertices[197],vertices[244],vertices[196],vertices[204]);
	  cells[1244]->set_neighbors(cells[952],cells[168],cells[1521],cells[1426]);
	  set_offsets(cells[1244],0,0,0,0);
	  cells[1245]->set_vertices(vertices[279],vertices[278],vertices[38],vertices[271]);
	  cells[1245]->set_neighbors(cells[1581],cells[1238],cells[936],cells[1716]);
	  set_offsets(cells[1245],0,0,4,0);
	  cells[1246]->set_vertices(vertices[277],vertices[37],vertices[38],vertices[29]);
	  cells[1246]->set_neighbors(cells[553],cells[1181],cells[1629],cells[1474]);
	  set_offsets(cells[1246],0,4,4,4);
	  cells[1247]->set_vertices(vertices[279],vertices[280],vertices[239],vertices[231]);
	  cells[1247]->set_neighbors(cells[1597],cells[1714],cells[1292],cells[957]);
	  set_offsets(cells[1247],0,1,0,0);
	  cells[1248]->set_vertices(vertices[271],vertices[264],vertices[263],vertices[24]);
	  cells[1248]->set_neighbors(cells[1233],cells[1590],cells[1547],cells[1637]);
	  set_offsets(cells[1248],0,1,0,5);
	  cells[1249]->set_vertices(vertices[254],vertices[206],vertices[199],vertices[207]);
	  cells[1249]->set_neighbors(cells[1205],cells[173],cells[1493],cells[1300]);
	  set_offsets(cells[1249],0,0,0,0);
	  cells[1250]->set_vertices(vertices[213],vertices[165],vertices[205],vertices[214]);
	  cells[1250]->set_neighbors(cells[865],cells[1290],cells[1338],cells[1241]);
	  set_offsets(cells[1250],0,0,0,0);
	  cells[1251]->set_vertices(vertices[28],vertices[20],vertices[21],vertices[261]);
	  cells[1251]->set_neighbors(cells[1573],cells[136],cells[128],cells[104]);
	  set_offsets(cells[1251],4,4,4,0);
	  cells[1252]->set_vertices(vertices[211],vertices[203],vertices[251],vertices[260]);
	  cells[1252]->set_neighbors(cells[1020],cells[784],cells[1336],cells[391]);
	  set_offsets(cells[1252],0,0,0,0);
	  cells[1253]->set_vertices(vertices[12],vertices[251],vertices[260],vertices[20]);
	  cells[1253]->set_neighbors(cells[1280],cells[1293],cells[183],cells[1576]);
	  set_offsets(cells[1253],4,0,0,4);
	  cells[1254]->set_vertices(vertices[165],vertices[212],vertices[157],vertices[164]);
	  cells[1254]->set_neighbors(cells[1258],cells[1018],cells[983],cells[1255]);
	  set_offsets(cells[1254],0,0,0,0);
	  cells[1255]->set_vertices(vertices[165],vertices[212],vertices[205],vertices[157]);
	  cells[1255]->set_neighbors(cells[265],cells[865],cells[1254],cells[1241]);
	  set_offsets(cells[1255],0,0,0,0);
	  cells[1256]->set_vertices(vertices[243],vertices[195],vertices[252],vertices[203]);
	  cells[1256]->set_neighbors(cells[1159],cells[1257],cells[1242],cells[1453]);
	  set_offsets(cells[1256],0,0,0,0);
	  cells[1257]->set_vertices(vertices[251],vertices[243],vertices[252],vertices[203]);
	  cells[1257]->set_neighbors(cells[1256],cells[1020],cells[1519],cells[124]);
	  set_offsets(cells[1257],0,0,0,0);
	  cells[1258]->set_vertices(vertices[212],vertices[204],vertices[157],vertices[164]);
	  cells[1258]->set_neighbors(cells[1014],cells[1254],cells[59],cells[265]);
	  set_offsets(cells[1258],0,0,0,0);
	  cells[1259]->set_vertices(vertices[192],vertices[144],vertices[191],vertices[232]);
	  cells[1259]->set_neighbors(cells[1397],cells[613],cells[1053],cells[909]);
	  set_offsets(cells[1259],3,3,0,1);
	  cells[1260]->set_vertices(vertices[214],vertices[206],vertices[207],vertices[159]);
	  cells[1260]->set_neighbors(cells[1205],cells[1161],cells[1264],cells[1493]);
	  set_offsets(cells[1260],0,0,0,0);
	  cells[1261]->set_vertices(vertices[31],vertices[271],vertices[23],vertices[24]);
	  cells[1261]->set_neighbors(cells[1590],cells[2],cells[1413],cells[392]);
	  set_offsets(cells[1261],4,0,4,5);
	  cells[1262]->set_vertices(vertices[219],vertices[259],vertices[267],vertices[268]);
	  cells[1262]->set_neighbors(cells[1526],cells[1624],cells[1575],cells[1622]);
	  set_offsets(cells[1262],0,0,0,0);
	  cells[1263]->set_vertices(vertices[255],vertices[247],vertices[7],vertices[8]);
	  cells[1263]->set_neighbors(cells[963],cells[1304],cells[1534],cells[708]);
	  set_offsets(cells[1263],0,0,4,5);
	  cells[1264]->set_vertices(vertices[214],vertices[206],vertices[159],vertices[166]);
	  cells[1264]->set_neighbors(cells[1031],cells[1160],cells[1040],cells[1260]);
	  set_offsets(cells[1264],0,0,0,0);
	  cells[1265]->set_vertices(vertices[245],vertices[197],vertices[254],vertices[205]);
	  cells[1265]->set_neighbors(cells[1225],cells[1305],cells[1469],cells[1026]);
	  set_offsets(cells[1265],0,0,0,0);
	  cells[1266]->set_vertices(vertices[255],vertices[256],vertices[8],vertices[16]);
	  cells[1266]->set_neighbors(cells[1278],cells[587],cells[1169],cells[1548]);
	  set_offsets(cells[1266],0,1,5,5);
	  cells[1267]->set_vertices(vertices[160],vertices[200],vertices[159],vertices[208]);
	  cells[1267]->set_neighbors(cells[1141],cells[1268],cells[997],cells[1214]);
	  set_offsets(cells[1267],1,1,0,1);
	  cells[1268]->set_vertices(vertices[160],vertices[208],vertices[159],vertices[167]);
	  cells[1268]->set_neighbors(cells[117],cells[1036],cells[1269],cells[1267]);
	  set_offsets(cells[1268],1,1,0,0);
	  cells[1269]->set_vertices(vertices[168],vertices[160],vertices[167],vertices[208]);
	  cells[1269]->set_neighbors(cells[1268],cells[1312],cells[1087],cells[741]);
	  set_offsets(cells[1269],1,1,0,1);
	  cells[1270]->set_vertices(vertices[217],vertices[169],vertices[216],vertices[209]);
	  cells[1270]->set_neighbors(cells[1274],cells[1446],cells[1307],cells[1320]);
	  set_offsets(cells[1270],0,0,0,0);
	  cells[1271]->set_vertices(vertices[264],vertices[257],vertices[256],vertices[16]);
	  cells[1271]->set_neighbors(cells[1552],cells[1316],cells[1324],cells[1599]);
	  set_offsets(cells[1271],0,0,0,4);
	  cells[1272]->set_vertices(vertices[219],vertices[171],vertices[211],vertices[220]);
	  cells[1272]->set_neighbors(cells[1322],cells[1386],cells[962],cells[1275]);
	  set_offsets(cells[1272],0,0,0,0);
	  cells[1273]->set_vertices(vertices[216],vertices[161],vertices[208],vertices[209]);
	  cells[1273]->set_neighbors(cells[138],cells[1226],cells[1274],cells[216]);
	  set_offsets(cells[1273],0,0,0,0);
	  cells[1274]->set_vertices(vertices[169],vertices[161],vertices[216],vertices[209]);
	  cells[1274]->set_neighbors(cells[1273],cells[1270],cells[1049],cells[994]);
	  set_offsets(cells[1274],0,0,0,0);
	  cells[1275]->set_vertices(vertices[219],vertices[218],vertices[211],vertices[171]);
	  cells[1275]->set_neighbors(cells[1287],cells[1272],cells[1326],cells[1334]);
	  set_offsets(cells[1275],0,0,0,0);
	  cells[1276]->set_vertices(vertices[230],vertices[270],vertices[223],vertices[222]);
	  cells[1276]->set_neighbors(cells[1592],cells[1351],cells[673],cells[1636]);
	  set_offsets(cells[1276],0,0,0,0);
	  cells[1277]->set_vertices(vertices[251],vertices[250],vertices[258],vertices[10]);
	  cells[1277]->set_neighbors(cells[1510],cells[749],cells[1007],cells[1549]);
	  set_offsets(cells[1277],0,0,0,4);
	  cells[1278]->set_vertices(vertices[256],vertices[249],vertices[8],vertices[16]);
	  cells[1278]->set_neighbors(cells[1505],cells[1266],cells[1552],cells[1544]);
	  set_offsets(cells[1278],0,0,4,4);
	  cells[1279]->set_vertices(vertices[12],vertices[20],vertices[253],vertices[13]);
	  cells[1279]->set_neighbors(cells[1572],cells[1532],cells[319],cells[1293]);
	  set_offsets(cells[1279],4,4,0,4);
	  cells[1280]->set_vertices(vertices[259],vertices[251],vertices[20],vertices[260]);
	  cells[1280]->set_neighbors(cells[1253],cells[1023],cells[784],cells[1566]);
	  set_offsets(cells[1280],0,0,4,0);
	  cells[1281]->set_vertices(vertices[212],vertices[155],vertices[203],vertices[204]);
	  cells[1281]->set_neighbors(cells[1011],cells[1066],cells[59],cells[441]);
	  set_offsets(cells[1281],0,0,0,0);
	  cells[1282]->set_vertices(vertices[218],vertices[210],vertices[211],vertices[163]);
	  cells[1282]->set_neighbors(cells[364],cells[1287],cells[1289],cells[1516]);
	  set_offsets(cells[1282],0,0,0,0);
	  cells[1283]->set_vertices(vertices[244],vertices[4],vertices[283],vertices[243]);
	  cells[1283]->set_neighbors(cells[1456],cells[1663],cells[929],cells[237]);
	  set_offsets(cells[1283],2,6,0,2);
	  cells[1284]->set_vertices(vertices[220],vertices[260],vertices[213],vertices[212]);
	  cells[1284]->set_neighbors(cells[964],cells[1295],cells[832],cells[262]);
	  set_offsets(cells[1284],0,0,0,0);
	  cells[1285]->set_vertices(vertices[258],vertices[210],vertices[250],vertices[203]);
	  cells[1285]->set_neighbors(cells[947],cells[1549],cells[1048],cells[876]);
	  set_offsets(cells[1285],0,0,0,0);
	  cells[1286]->set_vertices(vertices[171],vertices[218],vertices[163],vertices[170]);
	  cells[1286]->set_neighbors(cells[1289],cells[1055],cells[1080],cells[1287]);
	  set_offsets(cells[1286],0,0,0,0);
	  cells[1287]->set_vertices(vertices[171],vertices[218],vertices[211],vertices[163]);
	  cells[1287]->set_neighbors(cells[1282],cells[1322],cells[1286],cells[1275]);
	  set_offsets(cells[1287],0,0,0,0);
	  cells[1288]->set_vertices(vertices[209],vertices[201],vertices[258],vertices[210]);
	  cells[1288]->set_neighbors(cells[876],cells[1517],cells[960],cells[1561]);
	  set_offsets(cells[1288],0,0,0,0);
	  cells[1289]->set_vertices(vertices[218],vertices[210],vertices[163],vertices[170]);
	  cells[1289]->set_neighbors(cells[775],cells[1286],cells[248],cells[1282]);
	  set_offsets(cells[1289],0,0,0,0);
	  cells[1290]->set_vertices(vertices[213],vertices[205],vertices[262],vertices[214]);
	  cells[1290]->set_neighbors(cells[1483],cells[172],cells[1250],cells[1588]);
	  set_offsets(cells[1290],0,0,0,0);
	  cells[1291]->set_vertices(vertices[255],vertices[254],vertices[247],vertices[207]);
	  cells[1291]->set_neighbors(cells[173],cells[1108],cells[1538],cells[944]);
	  set_offsets(cells[1291],0,0,0,0);
	  cells[1292]->set_vertices(vertices[279],vertices[272],vertices[280],vertices[231]);
	  cells[1292]->set_neighbors(cells[1596],cells[1247],cells[1093],cells[1607]);
	  set_offsets(cells[1292],0,1,1,0);
	  cells[1293]->set_vertices(vertices[12],vertices[260],vertices[253],vertices[20]);
	  cells[1293]->set_neighbors(cells[1188],cells[1279],cells[1253],cells[1577]);
	  set_offsets(cells[1293],4,0,0,4);
	  cells[1294]->set_vertices(vertices[214],vertices[157],vertices[205],vertices[206]);
	  cells[1294]->set_neighbors(cells[460],cells[1211],cells[1040],cells[865]);
	  set_offsets(cells[1294],0,0,0,0);
	  cells[1295]->set_vertices(vertices[220],vertices[212],vertices[213],vertices[165]);
	  cells[1295]->set_neighbors(cells[1241],cells[76],cells[1301],cells[1284]);
	  set_offsets(cells[1295],0,0,0,0);
	  cells[1296]->set_vertices(vertices[223],vertices[270],vertices[263],vertices[215]);
	  cells[1296]->set_neighbors(cells[779],cells[1603],cells[1592],cells[1391]);
	  set_offsets(cells[1296],0,0,0,0);
	  cells[1297]->set_vertices(vertices[231],vertices[230],vertices[223],vertices[183]);
	  cells[1297]->set_neighbors(cells[851],cells[981],cells[972],cells[709]);
	  set_offsets(cells[1297],0,0,0,0);
	  cells[1298]->set_vertices(vertices[28],vertices[259],vertices[20],vertices[268]);
	  cells[1298]->set_neighbors(cells[1023],cells[128],cells[1526],cells[1467]);
	  set_offsets(cells[1298],4,0,4,0);
	  cells[1299]->set_vertices(vertices[231],vertices[272],vertices[223],vertices[271]);
	  cells[1299]->set_neighbors(cells[1348],cells[1643],cells[1093],cells[995]);
	  set_offsets(cells[1299],0,1,0,0);
	  cells[1300]->set_vertices(vertices[246],vertices[206],vertices[199],vertices[254]);
	  cells[1300]->set_neighbors(cells[1249],cells[769],cells[1212],cells[1028]);
	  set_offsets(cells[1300],0,0,0,0);
	  cells[1301]->set_vertices(vertices[220],vertices[212],vertices[165],vertices[172]);
	  cells[1301]->set_neighbors(cells[983],cells[1302],cells[81],cells[1295]);
	  set_offsets(cells[1301],0,0,0,0);
	  cells[1302]->set_vertices(vertices[173],vertices[220],vertices[165],vertices[172]);
	  cells[1302]->set_neighbors(cells[1301],cells[953],cells[155],cells[76]);
	  set_offsets(cells[1302],0,0,0,0);
	  cells[1303]->set_vertices(vertices[278],vertices[221],vertices[270],vertices[230]);
	  cells[1303]->set_neighbors(cells[673],cells[1636],cells[1166],cells[1479]);
	  set_offsets(cells[1303],0,0,0,0);
	  cells[1304]->set_vertices(vertices[15],vertices[255],vertices[7],vertices[8]);
	  cells[1304]->set_neighbors(cells[1263],cells[112],cells[587],cells[1495]);
	  set_offsets(cells[1304],4,0,4,5);
	  cells[1305]->set_vertices(vertices[253],vertices[245],vertices[254],vertices[205]);
	  cells[1305]->set_neighbors(cells[1265],cells[1482],cells[110],cells[1177]);
	  set_offsets(cells[1305],0,0,0,0);
	  cells[1306]->set_vertices(vertices[222],vertices[214],vertices[215],vertices[167]);
	  cells[1306]->set_neighbors(cells[1203],cells[810],cells[1311],cells[1349]);
	  set_offsets(cells[1306],0,0,0,0);
	  cells[1307]->set_vertices(vertices[217],vertices[169],vertices[209],vertices[218]);
	  cells[1307]->set_neighbors(cells[1049],cells[1369],cells[873],cells[1270]);
	  set_offsets(cells[1307],0,0,0,0);
	  cells[1308]->set_vertices(vertices[215],vertices[262],vertices[263],vertices[255]);
	  cells[1308]->set_neighbors(cells[1025],cells[1582],cells[1490],cells[779]);
	  set_offsets(cells[1308],0,0,0,0);
	  cells[1309]->set_vertices(vertices[191],vertices[184],vertices[183],vertices[232]);
	  cells[1309]->set_neighbors(cells[1398],cells[1401],cells[1397],cells[1154]);
	  set_offsets(cells[1309],0,1,0,1);
	  cells[1310]->set_vertices(vertices[175],vertices[222],vertices[167],vertices[174]);
	  cells[1310]->set_neighbors(cells[1311],cells[1077],cells[847],cells[810]);
	  set_offsets(cells[1310],0,0,0,0);
	  cells[1311]->set_vertices(vertices[222],vertices[214],vertices[167],vertices[174]);
	  cells[1311]->set_neighbors(cells[1075],cells[1310],cells[256],cells[1306]);
	  set_offsets(cells[1311],0,0,0,0);
	  cells[1312]->set_vertices(vertices[168],vertices[208],vertices[167],vertices[216]);
	  cells[1312]->set_neighbors(cells[919],cells[1313],cells[216],cells[1269]);
	  set_offsets(cells[1312],1,1,0,1);
	  cells[1313]->set_vertices(vertices[168],vertices[216],vertices[167],vertices[175]);
	  cells[1313]->set_neighbors(cells[69],cells[1078],cells[1314],cells[1312]);
	  set_offsets(cells[1313],1,1,0,0);
	  cells[1314]->set_vertices(vertices[176],vertices[168],vertices[175],vertices[216]);
	  cells[1314]->set_neighbors(cells[1313],cells[1360],cells[16],cells[977]);
	  set_offsets(cells[1314],1,1,0,1);
	  cells[1315]->set_vertices(vertices[272],vertices[217],vertices[224],vertices[264]);
	  cells[1315]->set_neighbors(cells[157],cells[1170],cells[1647],cells[277]);
	  set_offsets(cells[1315],0,0,0,0);
	  cells[1316]->set_vertices(vertices[264],vertices[256],vertices[263],vertices[16]);
	  cells[1316]->set_neighbors(cells[1169],cells[1233],cells[1271],cells[1602]);
	  set_offsets(cells[1316],1,1,0,5);
	  cells[1317]->set_vertices(vertices[41],vertices[2],vertices[1],vertices[241]);
	  cells[1317]->set_neighbors(cells[1449],cells[1402],cells[1231],cells[328]);
	  set_offsets(cells[1317],4,6,6,2);
	  cells[1318]->set_vertices(vertices[272],vertices[265],vertices[264],vertices[24]);
	  cells[1318]->set_neighbors(cells[1609],cells[1547],cells[910],cells[1647]);
	  set_offsets(cells[1318],0,0,0,4);
	  cells[1319]->set_vertices(vertices[213],vertices[260],vertices[253],vertices[205]);
	  cells[1319]->set_neighbors(cells[1465],cells[1588],cells[964],cells[1344]);
	  set_offsets(cells[1319],0,0,0,0);
	  cells[1320]->set_vertices(vertices[224],vertices[169],vertices[216],vertices[217]);
	  cells[1320]->set_neighbors(cells[1270],cells[157],cells[1321],cells[254]);
	  set_offsets(cells[1320],0,0,0,0);
	  cells[1321]->set_vertices(vertices[177],vertices[169],vertices[224],vertices[217]);
	  cells[1321]->set_neighbors(cells[1320],cells[1174],cells[1186],cells[1358]);
	  set_offsets(cells[1321],0,0,0,0);
	  cells[1322]->set_vertices(vertices[171],vertices[163],vertices[211],vertices[220]);
	  cells[1322]->set_neighbors(cells[1062],cells[1272],cells[233],cells[1287]);
	  set_offsets(cells[1322],0,0,0,0);
	  cells[1323]->set_vertices(vertices[233],vertices[185],vertices[234],vertices[194]);
	  cells[1323]->set_neighbors(cells[1414],cells[1617],cells[1408],cells[1415]);
	  set_offsets(cells[1323],0,0,0,2);
	  cells[1324]->set_vertices(vertices[24],vertices[257],vertices[264],vertices[16]);
	  cells[1324]->set_neighbors(cells[1271],cells[1233],cells[1553],cells[1609]);
	  set_offsets(cells[1324],4,0,0,4);
	  cells[1325]->set_vertices(vertices[180],vertices[188],vertices[181],vertices[228]);
	  cells[1325]->set_neighbors(cells[327],cells[1342],cells[925],cells[490]);
	  set_offsets(cells[1325],0,0,0,0);
	  cells[1326]->set_vertices(vertices[226],vertices[218],vertices[219],vertices[171]);
	  cells[1326]->set_neighbors(cells[1275],cells[1333],cells[1337],cells[1373]);
	  set_offsets(cells[1326],0,0,0,0);
	  cells[1327]->set_vertices(vertices[187],vertices[148],vertices[188],vertices[236]);
	  cells[1327]->set_neighbors(cells[1383],cells[1392],cells[768],cells[1134]);
	  set_offsets(cells[1327],0,2,0,0);
	  cells[1328]->set_vertices(vertices[42],vertices[282],vertices[34],vertices[275]);
	  cells[1328]->set_neighbors(cells[1680],cells[1618],cells[1698],cells[1690]);
	  set_offsets(cells[1328],4,0,4,0);
	  cells[1329]->set_vertices(vertices[233],vertices[273],vertices[281],vertices[282]);
	  cells[1329]->set_neighbors(cells[1689],cells[1679],cells[1693],cells[1652]);
	  set_offsets(cells[1329],0,0,0,0);
	  cells[1330]->set_vertices(vertices[225],vertices[217],vertices[265],vertices[274]);
	  cells[1330]->set_neighbors(cells[1367],cells[930],cells[85],cells[1646]);
	  set_offsets(cells[1330],0,0,0,0);
	  cells[1331]->set_vertices(vertices[268],vertices[211],vertices[259],vertices[260]);
	  cells[1331]->set_neighbors(cells[784],cells[1023],cells[1365],cells[1575]);
	  set_offsets(cells[1331],0,0,0,0);
	  cells[1332]->set_vertices(vertices[178],vertices[179],vertices[171],vertices[226]);
	  cells[1332]->set_neighbors(cells[1333],cells[1337],cells[771],cells[1097]);
	  set_offsets(cells[1332],0,0,0,0);
	  cells[1333]->set_vertices(vertices[179],vertices[226],vertices[219],vertices[171]);
	  cells[1333]->set_neighbors(cells[1326],cells[1122],cells[1332],cells[1372]);
	  set_offsets(cells[1333],0,0,0,0);
	  cells[1334]->set_vertices(vertices[219],vertices[218],vertices[266],vertices[211]);
	  cells[1334]->set_neighbors(cells[1335],cells[958],cells[1275],cells[1373]);
	  set_offsets(cells[1334],0,0,0,0);
	  cells[1335]->set_vertices(vertices[266],vertices[218],vertices[258],vertices[211]);
	  cells[1335]->set_neighbors(cells[1516],cells[1564],cells[1334],cells[1110]);
	  set_offsets(cells[1335],0,0,0,0);
	  cells[1336]->set_vertices(vertices[211],vertices[203],vertices[260],vertices[212]);
	  cells[1336]->set_neighbors(cells[1024],cells[832],cells[905],cells[1252]);
	  set_offsets(cells[1336],0,0,0,0);
	  cells[1337]->set_vertices(vertices[178],vertices[226],vertices[171],vertices[218]);
	  cells[1337]->set_neighbors(cells[1326],cells[1080],cells[4],cells[1332]);
	  set_offsets(cells[1337],0,0,0,0);
	  cells[1338]->set_vertices(vertices[222],vertices[165],vertices[213],vertices[214]);
	  cells[1338]->set_neighbors(cells[1250],cells[172],cells[256],cells[1060]);
	  set_offsets(cells[1338],0,0,0,0);
	  cells[1339]->set_vertices(vertices[228],vertices[220],vertices[221],vertices[173]);
	  cells[1339]->set_neighbors(cells[1047],cells[1343],cells[1345],cells[1366]);
	  set_offsets(cells[1339],0,0,0,0);
	  cells[1340]->set_vertices(vertices[246],vertices[6],vertices[46],vertices[285]);
	  cells[1340]->set_neighbors(cells[1628],cells[1208],cells[1056],cells[1400]);
	  set_offsets(cells[1340],2,6,4,0);
	  cells[1341]->set_vertices(vertices[283],vertices[275],vertices[42],vertices[35]);
	  cells[1341]->set_neighbors(cells[1618],cells[1139],cells[1666],cells[1698]);
	  set_offsets(cells[1341],0,0,4,4);
	  cells[1342]->set_vertices(vertices[180],vertices[181],vertices[173],vertices[228]);
	  cells[1342]->set_neighbors(cells[1343],cells[1345],cells[1325],cells[1107]);
	  set_offsets(cells[1342],0,0,0,0);
	  cells[1343]->set_vertices(vertices[181],vertices[228],vertices[221],vertices[173]);
	  cells[1343]->set_neighbors(cells[1339],cells[1109],cells[1342],cells[1375]);
	  set_offsets(cells[1343],0,0,0,0);
	  cells[1344]->set_vertices(vertices[213],vertices[260],vertices[261],vertices[253]);
	  cells[1344]->set_neighbors(cells[1188],cells[1587],cells[1319],cells[522]);
	  set_offsets(cells[1344],0,0,0,0);
	  cells[1345]->set_vertices(vertices[180],vertices[228],vertices[173],vertices[220]);
	  cells[1345]->set_neighbors(cells[1339],cells[155],cells[259],cells[1342]);
	  set_offsets(cells[1345],0,0,0,0);
	  cells[1346]->set_vertices(vertices[280],vertices[273],vertices[272],vertices[32]);
	  cells[1346]->set_neighbors(cells[1657],cells[1607],cells[1645],cells[1433]);
	  set_offsets(cells[1346],0,0,0,4);
	  cells[1347]->set_vertices(vertices[185],vertices[184],vertices[232],vertices[177]);
	  cells[1347]->set_neighbors(cells[913],cells[988],cells[1090],cells[1086]);
	  set_offsets(cells[1347],0,0,0,0);
	  cells[1348]->set_vertices(vertices[272],vertices[264],vertices[223],vertices[271]);
	  cells[1348]->set_neighbors(cells[1637],cells[1299],cells[1547],cells[1170]);
	  set_offsets(cells[1348],1,1,0,0);
	  cells[1349]->set_vertices(vertices[222],vertices[262],vertices[215],vertices[214]);
	  cells[1349]->set_neighbors(cells[1350],cells[1306],cells[172],cells[1591]);
	  set_offsets(cells[1349],0,0,0,0);
	  cells[1350]->set_vertices(vertices[215],vertices[214],vertices[262],vertices[207]);
	  cells[1350]->set_neighbors(cells[468],cells[1490],cells[1203],cells[1349]);
	  set_offsets(cells[1350],0,0,0,0);
	  cells[1351]->set_vertices(vertices[230],vertices[222],vertices[223],vertices[175]);
	  cells[1351]->set_neighbors(cells[989],cells[851],cells[1357],cells[1276]);
	  set_offsets(cells[1351],0,0,0,0);
	  cells[1352]->set_vertices(vertices[224],vertices[176],vertices[184],vertices[183]);
	  cells[1352]->set_neighbors(cells[998],cells[1398],cells[1359],cells[1084]);
	  set_offsets(cells[1352],1,1,1,0);
	  cells[1353]->set_vertices(vertices[247],vertices[0],vertices[47],vertices[7]);
	  cells[1353]->set_neighbors(cells[27],cells[234],cells[963],cells[1713]);
	  set_offsets(cells[1353],2,7,4,6);
	  cells[1354]->set_vertices(vertices[45],vertices[285],vertices[44],vertices[37]);
	  cells[1354]->set_neighbors(cells[58],cells[535],cells[1219],cells[1471]);
	  set_offsets(cells[1354],4,0,4,4);
	  cells[1355]->set_vertices(vertices[199],vertices[151],vertices[239],vertices[192]);
	  cells[1355]->set_neighbors(cells[345],cells[1444],cells[1221],cells[1088]);
	  set_offsets(cells[1355],2,2,0,3);
	  cells[1356]->set_vertices(vertices[183],vertices[230],vertices[175],vertices[182]);
	  cells[1356]->set_neighbors(cells[1357],cells[996],cells[1164],cells[851]);
	  set_offsets(cells[1356],0,0,0,0);
	  cells[1357]->set_vertices(vertices[230],vertices[222],vertices[175],vertices[182]);
	  cells[1357]->set_neighbors(cells[847],cells[1356],cells[244],cells[1351]);
	  set_offsets(cells[1357],0,0,0,0);
	  cells[1358]->set_vertices(vertices[176],vertices[169],vertices[224],vertices[177]);
	  cells[1358]->set_neighbors(cells[1321],cells[1084],cells[1043],cells[254]);
	  set_offsets(cells[1358],0,0,0,0);
	  cells[1359]->set_vertices(vertices[176],vertices[224],vertices[175],vertices[183]);
	  cells[1359]->set_neighbors(cells[753],cells[1117],cells[1352],cells[1360]);
	  set_offsets(cells[1359],1,1,0,0);
	  cells[1360]->set_vertices(vertices[176],vertices[216],vertices[175],vertices[224]);
	  cells[1360]->set_neighbors(cells[1085],cells[1359],cells[254],cells[1314]);
	  set_offsets(cells[1360],1,1,0,1);
	  cells[1361]->set_vertices(vertices[192],vertices[145],vertices[233],vertices[193]);
	  cells[1361]->set_neighbors(cells[1371],cells[1443],cells[718],cells[1409]);
	  set_offsets(cells[1361],2,2,0,2);
	  cells[1362]->set_vertices(vertices[276],vertices[236],vertices[227],vertices[228]);
	  cells[1362]->set_neighbors(cells[269],cells[326],cells[714],cells[1194]);
	  set_offsets(cells[1362],0,0,0,0);
	  cells[1363]->set_vertices(vertices[237],vertices[189],vertices[238],vertices[198]);
	  cells[1363]->set_neighbors(cells[1430],cells[1641],cells[1428],cells[347]);
	  set_offsets(cells[1363],0,0,0,2);
	  cells[1364]->set_vertices(vertices[283],vertices[4],vertices[44],vertices[43]);
	  cells[1364]->set_neighbors(cells[523],cells[1004],cells[1456],cells[237]);
	  set_offsets(cells[1364],0,6,4,4);
	  cells[1365]->set_vertices(vertices[268],vertices[211],vertices[260],vertices[220]);
	  cells[1365]->set_neighbors(cells[832],cells[262],cells[1386],cells[1331]);
	  set_offsets(cells[1365],0,0,0,0);
	  cells[1366]->set_vertices(vertices[228],vertices[268],vertices[221],vertices[220]);
	  cells[1366]->set_neighbors(cells[1381],cells[1339],cells[943],cells[1461]);
	  set_offsets(cells[1366],0,0,0,0);
	  cells[1367]->set_vertices(vertices[274],vertices[217],vertices[265],vertices[266]);
	  cells[1367]->set_neighbors(cells[1059],cells[1385],cells[1178],cells[1330]);
	  set_offsets(cells[1367],0,0,0,0);
	  cells[1368]->set_vertices(vertices[227],vertices[226],vertices[179],vertices[234]);
	  cells[1368]->set_neighbors(cells[1384],cells[307],cells[1654],cells[1372]);
	  set_offsets(cells[1368],0,0,0,0);
	  cells[1369]->set_vertices(vertices[217],vertices[209],vertices[266],vertices[218]);
	  cells[1369]->set_neighbors(cells[1110],cells[1131],cells[1307],cells[1514]);
	  set_offsets(cells[1369],0,0,0,0);
	  cells[1370]->set_vertices(vertices[26],vertices[259],vertices[19],vertices[267]);
	  cells[1370]->set_neighbors(cells[1058],cells[1515],cells[1176],cells[1569]);
	  set_offsets(cells[1370],4,0,4,0);
	  cells[1371]->set_vertices(vertices[233],vertices[145],vertices[194],vertices[193]);
	  cells[1371]->set_neighbors(cells[1000],cells[1455],cells[1361],cells[1408]);
	  set_offsets(cells[1371],0,2,2,2);
	  cells[1372]->set_vertices(vertices[227],vertices[226],vertices[219],vertices[179]);
	  cells[1372]->set_neighbors(cells[1333],cells[1201],cells[1368],cells[916]);
	  set_offsets(cells[1372],0,0,0,0);
	  cells[1373]->set_vertices(vertices[226],vertices[266],vertices[219],vertices[218]);
	  cells[1373]->set_neighbors(cells[1334],cells[1326],cells[1131],cells[529]);
	  set_offsets(cells[1373],0,0,0,0);
	  cells[1374]->set_vertices(vertices[276],vertices[219],vertices[268],vertices[228]);
	  cells[1374]->set_neighbors(cells[943],cells[1461],cells[326],cells[1624]);
	  set_offsets(cells[1374],0,0,0,0);
	  cells[1375]->set_vertices(vertices[229],vertices[228],vertices[221],vertices[181]);
	  cells[1375]->set_neighbors(cells[1343],cells[863],cells[1184],cells[1664]);
	  set_offsets(cells[1375],0,0,0,0);
	  cells[1376]->set_vertices(vertices[235],vertices[147],vertices[196],vertices[195]);
	  cells[1376]->set_neighbors(cells[309],cells[1458],cells[864],cells[1421]);
	  set_offsets(cells[1376],0,2,2,2);
	  cells[1377]->set_vertices(vertices[244],vertices[195],vertices[204],vertices[252]);
	  cells[1377]->set_neighbors(cells[1159],cells[1521],cells[1453],cells[952]);
	  set_offsets(cells[1377],0,0,0,0);
	  cells[1378]->set_vertices(vertices[209],vertices[249],vertices[257],vertices[258]);
	  cells[1378]->set_neighbors(cells[145],cells[1559],cells[1561],cells[1551]);
	  set_offsets(cells[1378],0,0,0,0);
	  cells[1379]->set_vertices(vertices[4],vertices[243],vertices[252],vertices[12]);
	  cells[1379]->set_neighbors(cells[124],cells[1175],cells[1457],cells[929]);
	  set_offsets(cells[1379],4,0,0,4);
	  cells[1380]->set_vertices(vertices[242],vertices[195],vertices[283],vertices[243]);
	  cells[1380]->set_neighbors(cells[1663],cells[1619],cells[280],cells[1700]);
	  set_offsets(cells[1380],2,2,0,2);
	  cells[1381]->set_vertices(vertices[221],vertices[220],vertices[268],vertices[213]);
	  cells[1381]->set_neighbors(cells[262],cells[891],cells[1047],cells[1366]);
	  set_offsets(cells[1381],0,0,0,0);
	  cells[1382]->set_vertices(vertices[275],vertices[267],vertices[36],vertices[276]);
	  cells[1382]->set_neighbors(cells[1045],cells[1704],cells[1626],cells[1389]);
	  set_offsets(cells[1382],0,0,4,0);
	  cells[1383]->set_vertices(vertices[236],vertices[148],vertices[188],vertices[189]);
	  cells[1383]->set_neighbors(cells[116],cells[1114],cells[126],cells[1327]);
	  set_offsets(cells[1383],0,2,0,0);
	  cells[1384]->set_vertices(vertices[226],vertices[234],vertices[186],vertices[179]);
	  cells[1384]->set_neighbors(cells[1171],cells[771],cells[1368],cells[1185]);
	  set_offsets(cells[1384],0,0,0,0);
	  cells[1385]->set_vertices(vertices[274],vertices[265],vertices[26],vertices[266]);
	  cells[1385]->set_neighbors(cells[570],cells[1620],cells[1367],cells[1563]);
	  set_offsets(cells[1385],0,0,4,0);
	  cells[1386]->set_vertices(vertices[219],vertices[211],vertices[268],vertices[220]);
	  cells[1386]->set_neighbors(cells[1365],cells[943],cells[1272],cells[1575]);
	  set_offsets(cells[1386],0,0,0,0);
	  cells[1387]->set_vertices(vertices[2],vertices[241],vertices[250],vertices[10]);
	  cells[1387]->set_neighbors(cells[1422],cells[1223],cells[1449],cells[1511]);
	  set_offsets(cells[1387],4,0,0,4);
	  cells[1388]->set_vertices(vertices[271],vertices[270],vertices[30],vertices[263]);
	  cells[1388]->set_neighbors(cells[1202],cells[1593],cells[1391],cells[1589]);
	  set_offsets(cells[1388],0,0,4,0);
	  cells[1389]->set_vertices(vertices[267],vertices[275],vertices[36],vertices[27]);
	  cells[1389]->set_neighbors(cells[1424],cells[1574],cells[874],cells[1382]);
	  set_offsets(cells[1389],0,0,4,4);
	  cells[1390]->set_vertices(vertices[269],vertices[261],vertices[30],vertices[270]);
	  cells[1390]->set_neighbors(cells[1144],cells[911],cells[1536],cells[1209]);
	  set_offsets(cells[1390],0,0,4,0);
	  cells[1391]->set_vertices(vertices[223],vertices[270],vertices[271],vertices[263]);
	  cells[1391]->set_neighbors(cells[1388],cells[1637],cells[1296],cells[1644]);
	  set_offsets(cells[1391],0,0,0,0);
	  cells[1392]->set_vertices(vertices[187],vertices[188],vertices[179],vertices[236]);
	  cells[1392]->set_neighbors(cells[133],cells[1200],cells[1327],cells[897]);
	  set_offsets(cells[1392],0,0,0,0);
	  cells[1393]->set_vertices(vertices[231],vertices[232],vertices[183],vertices[224]);
	  cells[1393]->set_neighbors(cells[1398],cells[981],cells[1410],cells[1401]);
	  set_offsets(cells[1393],0,1,0,1);
	  cells[1394]->set_vertices(vertices[239],vertices[280],vertices[192],vertices[232]);
	  cells[1394]->set_neighbors(cells[1642],cells[613],cells[1597],cells[1678]);
	  set_offsets(cells[1394],0,1,3,1);
	  cells[1395]->set_vertices(vertices[31],vertices[30],vertices[271],vertices[38]);
	  cells[1395]->set_neighbors(cells[1581],cells[1238],cells[464],cells[392]);
	  set_offsets(cells[1395],4,4,0,4);
	  cells[1396]->set_vertices(vertices[225],vertices[232],vertices[224],vertices[177]);
	  cells[1396]->set_neighbors(cells[913],cells[1174],cells[988],cells[1598]);
	  set_offsets(cells[1396],0,0,0,0);
	  cells[1397]->set_vertices(vertices[232],vertices[144],vertices[191],vertices[184]);
	  cells[1397]->set_neighbors(cells[1153],cells[1309],cells[1086],cells[1259]);
	  set_offsets(cells[1397],1,3,0,1);
	  cells[1398]->set_vertices(vertices[232],vertices[184],vertices[183],vertices[224]);
	  cells[1398]->set_neighbors(cells[1352],cells[1393],cells[913],cells[1309]);
	  set_offsets(cells[1398],1,1,0,1);
	  cells[1399]->set_vertices(vertices[246],vertices[199],vertices[239],vertices[287]);
	  cells[1399]->set_neighbors(cells[1640],cells[1634],cells[1722],cells[1486]);
	  set_offsets(cells[1399],2,2,0,0);
	  cells[1400]->set_vertices(vertices[6],vertices[246],vertices[46],vertices[287]);
	  cells[1400]->set_neighbors(cells[1162],cells[1488],cells[1723],cells[1340]);
	  set_offsets(cells[1400],6,2,4,0);
	  cells[1401]->set_vertices(vertices[231],vertices[191],vertices[183],vertices[232]);
	  cells[1401]->set_neighbors(cells[1309],cells[1393],cells[1437],cells[323]);
	  set_offsets(cells[1401],0,0,0,1);
	  cells[1402]->set_vertices(vertices[41],vertices[241],vertices[1],vertices[0]);
	  cells[1402]->set_neighbors(cells[1101],cells[1],cells[1685],cells[1317]);
	  set_offsets(cells[1402],4,2,6,6);
	  cells[1403]->set_vertices(vertices[242],vertices[193],vertices[250],vertices[241]);
	  cells[1403]->set_neighbors(cells[1502],cells[1511],cells[1404],cells[1512]);
	  set_offsets(cells[1403],0,0,0,0);
	  cells[1404]->set_vertices(vertices[241],vertices[193],vertices[281],vertices[242]);
	  cells[1404]->set_neighbors(cells[1499],cells[1540],cells[1403],cells[1683]);
	  set_offsets(cells[1404],2,2,0,2);
	  cells[1405]->set_vertices(vertices[40],vertices[240],vertices[287],vertices[280]);
	  cells[1405]->set_neighbors(cells[1676],cells[1180],cells[1688],cells[1677]);
	  set_offsets(cells[1405],5,3,0,1);
	  cells[1406]->set_vertices(vertices[233],vertices[185],vertices[192],vertices[232]);
	  cells[1406]->set_neighbors(cells[1053],cells[1642],cells[1407],cells[1409]);
	  set_offsets(cells[1406],0,0,2,0);
	  cells[1407]->set_vertices(vertices[225],vertices[185],vertices[233],vertices[232]);
	  cells[1407]->set_neighbors(cells[1406],cells[125],cells[988],cells[1415]);
	  set_offsets(cells[1407],0,0,0,0);
	  cells[1408]->set_vertices(vertices[185],vertices[145],vertices[194],vertices[233]);
	  cells[1408]->set_neighbors(cells[1371],cells[1323],cells[1409],cells[602]);
	  set_offsets(cells[1408],0,2,2,0);
	  cells[1409]->set_vertices(vertices[192],vertices[145],vertices[185],vertices[233]);
	  cells[1409]->set_neighbors(cells[1408],cells[1406],cells[1361],cells[914]);
	  set_offsets(cells[1409],2,2,0,0);
	  cells[1410]->set_vertices(vertices[231],vertices[272],vertices[232],vertices[224]);
	  cells[1410]->set_neighbors(cells[1598],cells[1393],cells[995],cells[1596]);
	  set_offsets(cells[1410],0,1,1,1);
	  cells[1411]->set_vertices(vertices[261],vertices[253],vertices[22],vertices[262]);
	  cells[1411]->set_neighbors(cells[222],cells[1533],cells[1587],cells[1579]);
	  set_offsets(cells[1411],0,0,4,0);
	  cells[1412]->set_vertices(vertices[32],vertices[273],vertices[25],vertices[33]);
	  cells[1412]->set_neighbors(cells[948],cells[402],cells[1439],cells[1658]);
	  set_offsets(cells[1412],4,0,4,4);
	  cells[1413]->set_vertices(vertices[31],vertices[271],vertices[24],vertices[32]);
	  cells[1413]->set_neighbors(cells[1492],cells[499],cells[1489],cells[1261]);
	  set_offsets(cells[1413],4,0,5,5);
	  cells[1414]->set_vertices(vertices[185],vertices[146],vertices[234],vertices[194]);
	  cells[1414]->set_neighbors(cells[1115],cells[1323],cells[602],cells[65]);
	  set_offsets(cells[1414],0,2,0,2);
	  cells[1415]->set_vertices(vertices[233],vertices[185],vertices[225],vertices[234]);
	  cells[1415]->set_neighbors(cells[1416],cells[344],cells[1323],cells[1407]);
	  set_offsets(cells[1415],0,0,0,0);
	  cells[1416]->set_vertices(vertices[225],vertices[185],vertices[177],vertices[234]);
	  cells[1416]->set_neighbors(cells[1173],cells[889],cells[1415],cells[988]);
	  set_offsets(cells[1416],0,0,0,0);
	  cells[1417]->set_vertices(vertices[14],vertices[245],vertices[5],vertices[6]);
	  cells[1417]->set_neighbors(cells[1425],cells[62],cells[1466],cells[1530]);
	  set_offsets(cells[1417],4,0,4,4);
	  cells[1418]->set_vertices(vertices[227],vertices[235],vertices[234],vertices[187]);
	  cells[1418]->set_neighbors(cells[1419],cells[307],cells[1005],cells[252]);
	  set_offsets(cells[1418],0,0,0,0);
	  cells[1419]->set_vertices(vertices[235],vertices[194],vertices[234],vertices[187]);
	  cells[1419]->set_neighbors(cells[1115],cells[1418],cells[1420],cells[1694]);
	  set_offsets(cells[1419],0,2,0,0);
	  cells[1420]->set_vertices(vertices[194],vertices[147],vertices[187],vertices[235]);
	  cells[1420]->set_neighbors(cells[1421],cells[1419],cells[864],cells[840]);
	  set_offsets(cells[1420],2,2,0,0);
	  cells[1421]->set_vertices(vertices[187],vertices[147],vertices[196],vertices[235]);
	  cells[1421]->set_neighbors(cells[1376],cells[1196],cells[1420],cells[1133]);
	  set_offsets(cells[1421],0,2,2,0);
	  cells[1422]->set_vertices(vertices[250],vertices[241],vertices[249],vertices[10]);
	  cells[1422]->set_neighbors(cells[1448],cells[1510],cells[1387],cells[342]);
	  set_offsets(cells[1422],0,0,0,4);
	  cells[1423]->set_vertices(vertices[2],vertices[10],vertices[243],vertices[3]);
	  cells[1423]->set_neighbors(cells[945],cells[1228],cells[72],cells[1223]);
	  set_offsets(cells[1423],4,4,0,4);
	  cells[1424]->set_vertices(vertices[275],vertices[35],vertices[36],vertices[27]);
	  cells[1424]->set_neighbors(cells[503],cells[1389],cells[1660],cells[267]);
	  set_offsets(cells[1424],0,4,4,4);
	  cells[1425]->set_vertices(vertices[245],vertices[6],vertices[45],vertices[5]);
	  cells[1425]->set_neighbors(cells[12],cells[236],cells[1417],cells[1197]);
	  set_offsets(cells[1425],2,6,4,6);
	  cells[1426]->set_vertices(vertices[244],vertices[196],vertices[237],vertices[197]);
	  cells[1426]->set_neighbors(cells[1199],cells[1712],cells[1244],cells[1475]);
	  set_offsets(cells[1426],2,2,0,2);
	  cells[1427]->set_vertices(vertices[229],vertices[237],vertices[236],vertices[189]);
	  cells[1427]->set_neighbors(cells[238],cells[830],cells[347],cells[1463]);
	  set_offsets(cells[1427],0,0,0,0);
	  cells[1428]->set_vertices(vertices[189],vertices[149],vertices[198],vertices[237]);
	  cells[1428]->set_neighbors(cells[1182],cells[1363],cells[1429],cells[1216]);
	  set_offsets(cells[1428],0,2,2,0);
	  cells[1429]->set_vertices(vertices[196],vertices[149],vertices[189],vertices[237]);
	  cells[1429]->set_neighbors(cells[1428],cells[238],cells[1199],cells[862]);
	  set_offsets(cells[1429],2,2,0,0);
	  cells[1430]->set_vertices(vertices[198],vertices[150],vertices[189],vertices[238]);
	  cells[1430]->set_neighbors(cells[1167],cells[1363],cells[1431],cells[1216]);
	  set_offsets(cells[1430],2,2,0,0);
	  cells[1431]->set_vertices(vertices[191],vertices[150],vertices[198],vertices[238]);
	  cells[1431]->set_neighbors(cells[1430],cells[1435],cells[973],cells[1215]);
	  set_offsets(cells[1431],0,2,2,0);
	  cells[1432]->set_vertices(vertices[200],vertices[193],vertices[240],vertices[248]);
	  cells[1432]->set_neighbors(cells[301],cells[1498],cells[1481],cells[1143]);
	  set_offsets(cells[1432],0,0,0,0);
	  cells[1433]->set_vertices(vertices[273],vertices[225],vertices[280],vertices[272]);
	  cells[1433]->set_neighbors(cells[1487],cells[1346],cells[1656],cells[1648]);
	  set_offsets(cells[1433],0,0,0,0);
	  cells[1434]->set_vertices(vertices[191],vertices[239],vertices[231],vertices[238]);
	  cells[1434]->set_neighbors(cells[1674],cells[323],cells[1435],cells[1437]);
	  set_offsets(cells[1434],0,0,0,0);
	  cells[1435]->set_vertices(vertices[191],vertices[198],vertices[239],vertices[238]);
	  cells[1435]->set_neighbors(cells[1473],cells[1434],cells[1431],cells[1436]);
	  set_offsets(cells[1435],0,2,0,0);
	  cells[1436]->set_vertices(vertices[198],vertices[151],vertices[191],vertices[239]);
	  cells[1436]->set_neighbors(cells[345],cells[1435],cells[1088],cells[1215]);
	  set_offsets(cells[1436],2,2,0,0);
	  cells[1437]->set_vertices(vertices[231],vertices[191],vertices[232],vertices[239]);
	  cells[1437]->set_neighbors(cells[613],cells[1597],cells[1434],cells[1401]);
	  set_offsets(cells[1437],0,0,1,0);
	  cells[1438]->set_vertices(vertices[273],vertices[34],vertices[33],vertices[42]);
	  cells[1438]->set_neighbors(cells[449],cells[628],cells[1690],cells[948]);
	  set_offsets(cells[1438],0,4,4,4);
	  cells[1439]->set_vertices(vertices[40],vertices[32],vertices[33],vertices[273]);
	  cells[1439]->set_neighbors(cells[1412],cells[586],cells[1645],cells[160]);
	  set_offsets(cells[1439],4,4,4,0);
	  cells[1440]->set_vertices(vertices[227],vertices[274],vertices[267],vertices[219]);
	  cells[1440]->set_neighbors(cells[1152],cells[1625],cells[916],cells[1578]);
	  set_offsets(cells[1440],0,0,0,0);
	  cells[1441]->set_vertices(vertices[199],vertices[248],vertices[247],vertices[207]);
	  cells[1441]->set_neighbors(cells[1108],cells[173],cells[1450],cells[1497]);
	  set_offsets(cells[1441],0,1,0,0);
	  cells[1442]->set_vertices(vertices[0],vertices[247],vertices[248],vertices[8]);
	  cells[1442]->set_neighbors(cells[1534],cells[975],cells[963],cells[1496]);
	  set_offsets(cells[1442],5,0,1,5);
	  cells[1443]->set_vertices(vertices[240],vertices[192],vertices[233],vertices[193]);
	  cells[1443]->set_neighbors(cells[1361],cells[1682],cells[1143],cells[1052]);
	  set_offsets(cells[1443],2,2,0,2);
	  cells[1444]->set_vertices(vertices[199],vertices[192],vertices[239],vertices[240]);
	  cells[1444]->set_neighbors(cells[1678],cells[1640],cells[134],cells[1355]);
	  set_offsets(cells[1444],2,3,0,3);
	  cells[1445]->set_vertices(vertices[2],vertices[243],vertices[283],vertices[43]);
	  cells[1445]->set_neighbors(cells[1456],cells[931],cells[1228],cells[1619]);
	  set_offsets(cells[1445],6,2,0,4);
	  cells[1446]->set_vertices(vertices[217],vertices[209],vertices[216],vertices[264]);
	  cells[1446]->set_neighbors(cells[1218],cells[157],cells[1600],cells[1270]);
	  set_offsets(cells[1446],0,0,0,0);
	  cells[1447]->set_vertices(vertices[1],vertices[241],vertices[249],vertices[8]);
	  cells[1447]->set_neighbors(cells[1503],cells[1506],cells[1101],cells[1448]);
	  set_offsets(cells[1447],4,0,0,4);
	  cells[1448]->set_vertices(vertices[1],vertices[10],vertices[249],vertices[241]);
	  cells[1448]->set_neighbors(cells[1422],cells[1447],cells[1449],cells[1507]);
	  set_offsets(cells[1448],4,4,0,0);
	  cells[1449]->set_vertices(vertices[1],vertices[2],vertices[10],vertices[241]);
	  cells[1449]->set_neighbors(cells[1387],cells[1448],cells[1317],cells[53]);
	  set_offsets(cells[1449],4,4,4,0);
	  cells[1450]->set_vertices(vertices[199],vertices[200],vertices[248],vertices[207]);
	  cells[1450]->set_neighbors(cells[1451],cells[1441],cells[137],cells[1498]);
	  set_offsets(cells[1450],0,1,1,0);
	  cells[1451]->set_vertices(vertices[208],vertices[200],vertices[207],vertices[248]);
	  cells[1451]->set_neighbors(cells[1450],cells[1508],cells[1210],cells[1141]);
	  set_offsets(cells[1451],1,1,0,1);
	  cells[1452]->set_vertices(vertices[275],vertices[235],vertices[282],vertices[227]);
	  cells[1452]->set_neighbors(cells[252],cells[1649],cells[1702],cells[1695]);
	  set_offsets(cells[1452],0,0,0,0);
	  cells[1453]->set_vertices(vertices[243],vertices[195],vertices[244],vertices[252]);
	  cells[1453]->set_neighbors(cells[1377],cells[929],cells[1256],cells[1663]);
	  set_offsets(cells[1453],0,0,0,0);
	  cells[1454]->set_vertices(vertices[194],vertices[193],vertices[202],vertices[242]);
	  cells[1454]->set_neighbors(cells[1512],cells[1132],cells[1455],cells[1000]);
	  set_offsets(cells[1454],0,0,0,0);
	  cells[1455]->set_vertices(vertices[242],vertices[193],vertices[233],vertices[194]);
	  cells[1455]->set_neighbors(cells[1371],cells[226],cells[1454],cells[1499]);
	  set_offsets(cells[1455],2,2,0,2);
	  cells[1456]->set_vertices(vertices[243],vertices[4],vertices[283],vertices[43]);
	  cells[1456]->set_neighbors(cells[1364],cells[1445],cells[890],cells[1283]);
	  set_offsets(cells[1456],2,6,0,4);
	  cells[1457]->set_vertices(vertices[3],vertices[4],vertices[12],vertices[243]);
	  cells[1457]->set_neighbors(cells[1379],cells[1468],cells[890],cells[31]);
	  set_offsets(cells[1457],4,4,4,0);
	  cells[1458]->set_vertices(vertices[244],vertices[195],vertices[235],vertices[196]);
	  cells[1458]->set_neighbors(cells[1376],cells[1570],cells[952],cells[1662]);
	  set_offsets(cells[1458],2,2,0,2);
	  cells[1459]->set_vertices(vertices[284],vertices[235],vertices[227],vertices[236]);
	  cells[1459]->set_neighbors(cells[1005],cells[1194],cells[1462],cells[1702]);
	  set_offsets(cells[1459],0,0,0,0);
	  cells[1460]->set_vertices(vertices[276],vertices[268],vertices[269],vertices[221]);
	  cells[1460]->set_neighbors(cells[1633],cells[1140],cells[1461],cells[1631]);
	  set_offsets(cells[1460],0,0,0,0);
	  cells[1461]->set_vertices(vertices[276],vertices[228],vertices[268],vertices[221]);
	  cells[1461]->set_neighbors(cells[1366],cells[1460],cells[1664],cells[1374]);
	  set_offsets(cells[1461],0,0,0,0);
	  cells[1462]->set_vertices(vertices[284],vertices[196],vertices[235],vertices[236]);
	  cells[1462]->set_neighbors(cells[1196],cells[1459],cells[1706],cells[1570]);
	  set_offsets(cells[1462],0,2,0,0);
	  cells[1463]->set_vertices(vertices[284],vertices[237],vertices[236],vertices[229]);
	  cells[1463]->set_neighbors(cells[1427],cells[1665],cells[1003],cells[1706]);
	  set_offsets(cells[1463],0,0,0,0);
	  cells[1464]->set_vertices(vertices[21],vertices[269],vertices[29],vertices[28]);
	  cells[1464]->set_neighbors(cells[1527],cells[439],cells[136],cells[1584]);
	  set_offsets(cells[1464],4,0,4,4);
	  cells[1465]->set_vertices(vertices[260],vertices[252],vertices[253],vertices[205]);
	  cells[1465]->set_neighbors(cells[110],cells[1319],cells[95],cells[1577]);
	  set_offsets(cells[1465],0,0,0,0);
	  cells[1466]->set_vertices(vertices[6],vertices[245],vertices[254],vertices[14]);
	  cells[1466]->set_neighbors(cells[1177],cells[1239],cells[1417],cells[197]);
	  set_offsets(cells[1466],4,0,0,4);
	  cells[1467]->set_vertices(vertices[259],vertices[20],vertices[19],vertices[28]);
	  cells[1467]->set_neighbors(cells[271],cells[1058],cells[1298],cells[1567]);
	  set_offsets(cells[1467],0,4,4,4);
	  cells[1468]->set_vertices(vertices[3],vertices[12],vertices[251],vertices[243]);
	  cells[1468]->set_neighbors(cells[124],cells[945],cells[1457],cells[803]);
	  set_offsets(cells[1468],4,4,0,0);
	  cells[1469]->set_vertices(vertices[245],vertices[252],vertices[197],vertices[205]);
	  cells[1469]->set_neighbors(cells[1470],cells[1265],cells[110],cells[1522]);
	  set_offsets(cells[1469],0,0,0,0);
	  cells[1470]->set_vertices(vertices[252],vertices[204],vertices[197],vertices[205]);
	  cells[1470]->set_neighbors(cells[443],cells[1469],cells[1243],cells[1521]);
	  set_offsets(cells[1470],0,0,0,0);
	  cells[1471]->set_vertices(vertices[4],vertices[285],vertices[44],vertices[45]);
	  cells[1471]->set_neighbors(cells[1354],cells[575],cells[1207],cells[1711]);
	  set_offsets(cells[1471],6,0,4,4);
	  cells[1472]->set_vertices(vertices[229],vertices[269],vertices[277],vertices[278]);
	  cells[1472]->set_neighbors(cells[322],cells[1480],cells[1670],cells[1638]);
	  set_offsets(cells[1472],0,0,0,0);
	  cells[1473]->set_vertices(vertices[238],vertices[198],vertices[239],vertices[286]);
	  cells[1473]->set_neighbors(cells[1022],cells[1674],cells[1641],cells[1435]);
	  set_offsets(cells[1473],0,2,0,0);
	  cells[1474]->set_vertices(vertices[277],vertices[46],vertices[38],vertices[37]);
	  cells[1474]->set_neighbors(cells[275],cells[1246],cells[1172],cells[1718]);
	  set_offsets(cells[1474],0,4,4,4);
	  cells[1475]->set_vertices(vertices[237],vertices[196],vertices[244],vertices[284]);
	  cells[1475]->set_neighbors(cells[1570],cells[1708],cells[1706],cells[1426]);
	  set_offsets(cells[1475],0,2,2,0);
	  cells[1476]->set_vertices(vertices[231],vertices[286],vertices[278],vertices[238]);
	  cells[1476]->set_neighbors(cells[1673],cells[1675],cells[1674],cells[1715]);
	  set_offsets(cells[1476],0,0,0,0);
	  cells[1477]->set_vertices(vertices[279],vertices[272],vertices[271],vertices[32]);
	  cells[1477]->set_neighbors(cells[1492],cells[1489],cells[1607],cells[1093]);
	  set_offsets(cells[1477],0,1,0,5);
	  cells[1478]->set_vertices(vertices[286],vertices[198],vertices[246],vertices[237]);
	  cells[1478]->set_neighbors(cells[1206],cells[1148],cells[1641],cells[1022]);
	  set_offsets(cells[1478],0,2,2,0);
	  cells[1479]->set_vertices(vertices[278],vertices[221],vertices[269],vertices[270]);
	  cells[1479]->set_neighbors(cells[1536],cells[911],cells[1303],cells[1670]);
	  set_offsets(cells[1479],0,0,0,0);
	  cells[1480]->set_vertices(vertices[278],vertices[286],vertices[277],vertices[229]);
	  cells[1480]->set_neighbors(cells[623],cells[1472],cells[1673],cells[1719]);
	  set_offsets(cells[1480],0,0,0,0);
	  cells[1481]->set_vertices(vertices[200],vertices[193],vertices[248],vertices[201]);
	  cells[1481]->set_neighbors(cells[34],cells[1210],cells[1229],cells[1432]);
	  set_offsets(cells[1481],0,0,0,0);
	  cells[1482]->set_vertices(vertices[262],vertices[205],vertices[253],vertices[254]);
	  cells[1482]->set_neighbors(cells[1305],cells[1224],cells[1483],cells[1588]);
	  set_offsets(cells[1482],0,0,0,0);
	  cells[1483]->set_vertices(vertices[262],vertices[205],vertices[254],vertices[214]);
	  cells[1483]->set_neighbors(cells[1211],cells[468],cells[1290],cells[1482]);
	  set_offsets(cells[1483],0,0,0,0);
	  cells[1484]->set_vertices(vertices[12],vertices[252],vertices[245],vertices[253]);
	  cells[1484]->set_neighbors(cells[110],cells[1531],cells[1577],cells[1175]);
	  set_offsets(cells[1484],4,0,0,0);
	  cells[1485]->set_vertices(vertices[6],vertices[246],vertices[247],vertices[254]);
	  cells[1485]->set_neighbors(cells[769],cells[1239],cells[197],cells[1723]);
	  set_offsets(cells[1485],4,0,0,0);
	  cells[1486]->set_vertices(vertices[246],vertices[198],vertices[239],vertices[199]);
	  cells[1486]->set_neighbors(cells[1088],cells[1399],cells[1028],cells[1022]);
	  set_offsets(cells[1486],2,2,0,2);
	  cells[1487]->set_vertices(vertices[272],vertices[225],vertices[280],vertices[232]);
	  cells[1487]->set_neighbors(cells[125],cells[1596],cells[1598],cells[1433]);
	  set_offsets(cells[1487],0,0,0,0);
	  cells[1488]->set_vertices(vertices[47],vertices[46],vertices[287],vertices[6]);
	  cells[1488]->set_neighbors(cells[1400],cells[1724],cells[582],cells[1021]);
	  set_offsets(cells[1488],4,4,0,6);
	  cells[1489]->set_vertices(vertices[31],vertices[271],vertices[32],vertices[279]);
	  cells[1489]->set_neighbors(cells[1477],cells[0],cells[1238],cells[1413]);
	  set_offsets(cells[1489],4,0,5,0);
	  cells[1490]->set_vertices(vertices[215],vertices[262],vertices[255],vertices[207]);
	  cells[1490]->set_neighbors(cells[1538],cells[1539],cells[1350],cells[1308]);
	  set_offsets(cells[1490],0,0,0,0);
	  cells[1491]->set_vertices(vertices[263],vertices[255],vertices[15],vertices[16]);
	  cells[1491]->set_neighbors(cells[587],cells[1235],cells[1169],cells[1537]);
	  set_offsets(cells[1491],0,0,4,5);
	  cells[1492]->set_vertices(vertices[32],vertices[272],vertices[271],vertices[24]);
	  cells[1492]->set_neighbors(cells[1547],cells[1413],cells[910],cells[1477]);
	  set_offsets(cells[1492],5,1,0,5);
	  cells[1493]->set_vertices(vertices[214],vertices[254],vertices[207],vertices[206]);
	  cells[1493]->set_neighbors(cells[1249],cells[1260],cells[1211],cells[468]);
	  set_offsets(cells[1493],0,0,0,0);
	  cells[1494]->set_vertices(vertices[14],vertices[22],vertices[255],vertices[15]);
	  cells[1494]->set_neighbors(cells[1537],cells[1495],cells[18],cells[875]);
	  set_offsets(cells[1494],4,4,0,4);
	  cells[1495]->set_vertices(vertices[15],vertices[14],vertices[7],vertices[255]);
	  cells[1495]->set_neighbors(cells[708],cells[1304],cells[1494],cells[189]);
	  set_offsets(cells[1495],4,4,4,0);
	  cells[1496]->set_vertices(vertices[0],vertices[247],vertices[240],vertices[248]);
	  cells[1496]->set_neighbors(cells[1497],cells[229],cells[1442],cells[1721]);
	  set_offsets(cells[1496],5,0,1,1);
	  cells[1497]->set_vertices(vertices[199],vertices[240],vertices[247],vertices[248]);
	  cells[1497]->set_neighbors(cells[1496],cells[1441],cells[1498],cells[1720]);
	  set_offsets(cells[1497],0,1,0,1);
	  cells[1498]->set_vertices(vertices[199],vertices[200],vertices[240],vertices[248]);
	  cells[1498]->set_neighbors(cells[1432],cells[1497],cells[1450],cells[134]);
	  set_offsets(cells[1498],0,1,1,1);
	  cells[1499]->set_vertices(vertices[281],vertices[193],vertices[233],vertices[242]);
	  cells[1499]->set_neighbors(cells[1455],cells[1679],cells[1404],cells[1682]);
	  set_offsets(cells[1499],0,2,0,2);
	  cells[1500]->set_vertices(vertices[210],vertices[201],vertices[250],vertices[202]);
	  cells[1500]->set_neighbors(cells[1501],cells[947],cells[811],cells[876]);
	  set_offsets(cells[1500],0,0,0,0);
	  cells[1501]->set_vertices(vertices[202],vertices[193],vertices[201],vertices[250]);
	  cells[1501]->set_neighbors(cells[1502],cells[1500],cells[1512],cells[48]);
	  set_offsets(cells[1501],0,0,0,0);
	  cells[1502]->set_vertices(vertices[250],vertices[193],vertices[201],vertices[241]);
	  cells[1502]->set_neighbors(cells[34],cells[342],cells[1403],cells[1501]);
	  set_offsets(cells[1502],0,0,0,0);
	  cells[1503]->set_vertices(vertices[248],vertices[241],vertices[8],vertices[249]);
	  cells[1503]->set_neighbors(cells[1447],cells[1544],cells[1509],cells[975]);
	  set_offsets(cells[1503],0,0,4,0);
	  cells[1504]->set_vertices(vertices[256],vertices[248],vertices[207],vertices[255]);
	  cells[1504]->set_neighbors(cells[1108],cells[1539],cells[1548],cells[1508]);
	  set_offsets(cells[1504],1,1,0,0);
	  cells[1505]->set_vertices(vertices[249],vertices[16],vertices[9],vertices[8]);
	  cells[1505]->set_neighbors(cells[101],cells[1506],cells[1278],cells[1557]);
	  set_offsets(cells[1505],0,4,4,4);
	  cells[1506]->set_vertices(vertices[1],vertices[249],vertices[9],vertices[8]);
	  cells[1506]->set_neighbors(cells[1505],cells[336],cells[1447],cells[1507]);
	  set_offsets(cells[1506],4,0,4,4);
	  cells[1507]->set_vertices(vertices[1],vertices[10],vertices[9],vertices[249]);
	  cells[1507]->set_neighbors(cells[1006],cells[1506],cells[1448],cells[39]);
	  set_offsets(cells[1507],4,4,4,0);
	  cells[1508]->set_vertices(vertices[256],vertices[208],vertices[207],vertices[248]);
	  cells[1508]->set_neighbors(cells[1451],cells[1504],cells[940],cells[1546]);
	  set_offsets(cells[1508],1,1,0,1);
	  cells[1509]->set_vertices(vertices[248],vertices[241],vertices[249],vertices[201]);
	  cells[1509]->set_neighbors(cells[342],cells[1543],cells[34],cells[1503]);
	  set_offsets(cells[1509],0,0,0,0);
	  cells[1510]->set_vertices(vertices[250],vertices[249],vertices[258],vertices[10]);
	  cells[1510]->set_neighbors(cells[44],cells[1277],cells[1422],cells[1562]);
	  set_offsets(cells[1510],0,0,0,4);
	  cells[1511]->set_vertices(vertices[2],vertices[241],vertices[242],vertices[250]);
	  cells[1511]->set_neighbors(cells[1403],cells[1009],cells[1387],cells[1540]);
	  set_offsets(cells[1511],4,0,0,0);
	  cells[1512]->set_vertices(vertices[202],vertices[193],vertices[250],vertices[242]);
	  cells[1512]->set_neighbors(cells[1403],cells[1513],cells[1454],cells[1501]);
	  set_offsets(cells[1512],0,0,0,0);
	  cells[1513]->set_vertices(vertices[195],vertices[202],vertices[250],vertices[242]);
	  cells[1513]->set_neighbors(cells[1512],cells[280],cells[1132],cells[665]);
	  set_offsets(cells[1513],0,0,0,0);
	  cells[1514]->set_vertices(vertices[217],vertices[209],vertices[257],vertices[266]);
	  cells[1514]->set_neighbors(cells[1559],cells[1059],cells[1369],cells[1600]);
	  set_offsets(cells[1514],0,0,0,0);
	  cells[1515]->set_vertices(vertices[26],vertices[267],vertices[19],vertices[27]);
	  cells[1515]->set_neighbors(cells[1195],cells[285],cells[1565],cells[1370]);
	  set_offsets(cells[1515],4,0,4,4);
	  cells[1516]->set_vertices(vertices[218],vertices[258],vertices[211],vertices[210]);
	  cells[1516]->set_neighbors(cells[1048],cells[1282],cells[1517],cells[1335]);
	  set_offsets(cells[1516],0,0,0,0);
	  cells[1517]->set_vertices(vertices[218],vertices[209],vertices[258],vertices[210]);
	  cells[1517]->set_neighbors(cells[1288],cells[1516],cells[979],cells[1110]);
	  set_offsets(cells[1517],0,0,0,0);
	  cells[1518]->set_vertices(vertices[18],vertices[251],vertices[11],vertices[259]);
	  cells[1518]->set_neighbors(cells[1566],cells[1568],cells[421],cells[912]);
	  set_offsets(cells[1518],4,0,4,0);
	  cells[1519]->set_vertices(vertices[243],vertices[250],vertices[203],vertices[251]);
	  cells[1519]->set_neighbors(cells[1549],cells[1257],cells[1007],cells[1242]);
	  set_offsets(cells[1519],0,0,0,0);
	  cells[1520]->set_vertices(vertices[211],vertices[258],vertices[259],vertices[251]);
	  cells[1520]->set_neighbors(cells[421],cells[784],cells[391],cells[1564]);
	  set_offsets(cells[1520],0,0,0,0);
	  cells[1521]->set_vertices(vertices[244],vertices[204],vertices[197],vertices[252]);
	  cells[1521]->set_neighbors(cells[1470],cells[1522],cells[1377],cells[1244]);
	  set_offsets(cells[1521],0,0,0,0);
	  cells[1522]->set_vertices(vertices[245],vertices[244],vertices[197],vertices[252]);
	  cells[1522]->set_neighbors(cells[1521],cells[1469],cells[1523],cells[343]);
	  set_offsets(cells[1522],0,0,0,0);
	  cells[1523]->set_vertices(vertices[4],vertices[244],vertices[245],vertices[252]);
	  cells[1523]->set_neighbors(cells[1522],cells[1175],cells[929],cells[1630]);
	  set_offsets(cells[1523],4,0,0,0);
	  cells[1524]->set_vertices(vertices[275],vertices[274],vertices[34],vertices[267]);
	  cells[1524]->set_neighbors(cells[1621],cells[874],cells[1578],cells[1680]);
	  set_offsets(cells[1524],0,0,4,0);
	  cells[1525]->set_vertices(vertices[23],vertices[22],vertices[15],vertices[263]);
	  cells[1525]->set_neighbors(cells[1537],cells[1235],cells[974],cells[447]);
	  set_offsets(cells[1525],4,4,4,0);
	  cells[1526]->set_vertices(vertices[267],vertices[259],vertices[28],vertices[268]);
	  cells[1526]->set_neighbors(cells[1298],cells[664],cells[1262],cells[1058]);
	  set_offsets(cells[1526],0,0,4,0);
	  cells[1527]->set_vertices(vertices[36],vertices[269],vertices[28],vertices[29]);
	  cells[1527]->set_neighbors(cells[1464],cells[107],cells[1230],cells[1632]);
	  set_offsets(cells[1527],4,0,4,4);
	  cells[1528]->set_vertices(vertices[13],vertices[14],vertices[22],vertices[253]);
	  cells[1528]->set_neighbors(cells[222],cells[1579],cells[1529],cells[79]);
	  set_offsets(cells[1528],4,4,4,0);
	  cells[1529]->set_vertices(vertices[13],vertices[253],vertices[5],vertices[14]);
	  cells[1529]->set_neighbors(cells[1530],cells[278],cells[1528],cells[1532]);
	  set_offsets(cells[1529],4,0,4,4);
	  cells[1530]->set_vertices(vertices[253],vertices[245],vertices[5],vertices[14]);
	  cells[1530]->set_neighbors(cells[1417],cells[1529],cells[1177],cells[1531]);
	  set_offsets(cells[1530],0,0,4,4);
	  cells[1531]->set_vertices(vertices[253],vertices[12],vertices[5],vertices[245]);
	  cells[1531]->set_neighbors(cells[966],cells[1530],cells[1484],cells[1532]);
	  set_offsets(cells[1531],0,4,4,0);
	  cells[1532]->set_vertices(vertices[13],vertices[12],vertices[5],vertices[253]);
	  cells[1532]->set_neighbors(cells[1531],cells[1529],cells[1279],cells[206]);
	  set_offsets(cells[1532],4,4,4,0);
	  cells[1533]->set_vertices(vertices[270],vertices[261],vertices[22],vertices[262]);
	  cells[1533]->set_neighbors(cells[1411],cells[264],cells[1585],cells[1144]);
	  set_offsets(cells[1533],0,0,4,0);
	  cells[1534]->set_vertices(vertices[247],vertices[248],vertices[8],vertices[255]);
	  cells[1534]->set_neighbors(cells[1548],cells[1263],cells[1108],cells[1442]);
	  set_offsets(cells[1534],0,1,5,0);
	  cells[1535]->set_vertices(vertices[21],vertices[22],vertices[30],vertices[261]);
	  cells[1535]->set_neighbors(cells[1144],cells[1209],cells[1580],cells[406]);
	  set_offsets(cells[1535],4,4,4,0);
	  cells[1536]->set_vertices(vertices[221],vertices[261],vertices[269],vertices[270]);
	  cells[1536]->set_neighbors(cells[1390],cells[1479],cells[127],cells[1633]);
	  set_offsets(cells[1536],0,0,0,0);
	  cells[1537]->set_vertices(vertices[263],vertices[22],vertices[15],vertices[255]);
	  cells[1537]->set_neighbors(cells[1494],cells[1491],cells[1025],cells[1525]);
	  set_offsets(cells[1537],0,4,4,0);
	  cells[1538]->set_vertices(vertices[262],vertices[254],vertices[255],vertices[207]);
	  cells[1538]->set_neighbors(cells[1291],cells[1490],cells[468],cells[1234]);
	  set_offsets(cells[1538],0,0,0,0);
	  cells[1539]->set_vertices(vertices[215],vertices[256],vertices[207],vertices[255]);
	  cells[1539]->set_neighbors(cells[1504],cells[1490],cells[1582],cells[1546]);
	  set_offsets(cells[1539],0,1,0,0);
	  cells[1540]->set_vertices(vertices[242],vertices[281],vertices[241],vertices[2]);
	  cells[1540]->set_neighbors(cells[1231],cells[1511],cells[1232],cells[1404]);
	  set_offsets(cells[1540],2,0,2,6);
	  cells[1541]->set_vertices(vertices[40],vertices[280],vertices[279],vertices[32]);
	  cells[1541]->set_neighbors(cells[1607],cells[1220],cells[1645],cells[1180]);
	  set_offsets(cells[1541],5,1,0,5);
	  cells[1542]->set_vertices(vertices[209],vertices[201],vertices[256],vertices[249]);
	  cells[1542]->set_neighbors(cells[1543],cells[1551],cells[1561],cells[176]);
	  set_offsets(cells[1542],0,0,0,0);
	  cells[1543]->set_vertices(vertices[256],vertices[201],vertices[248],vertices[249]);
	  cells[1543]->set_neighbors(cells[1509],cells[1544],cells[1542],cells[940]);
	  set_offsets(cells[1543],0,0,0,0);
	  cells[1544]->set_vertices(vertices[248],vertices[249],vertices[8],vertices[256]);
	  cells[1544]->set_neighbors(cells[1278],cells[1548],cells[1543],cells[1503]);
	  set_offsets(cells[1544],0,0,4,0);
	  cells[1545]->set_vertices(vertices[216],vertices[208],vertices[215],vertices[256]);
	  cells[1545]->set_neighbors(cells[1546],cells[1601],cells[1226],cells[919]);
	  set_offsets(cells[1545],1,1,0,1);
	  cells[1546]->set_vertices(vertices[215],vertices[208],vertices[207],vertices[256]);
	  cells[1546]->set_neighbors(cells[1508],cells[1539],cells[1545],cells[1204]);
	  set_offsets(cells[1546],0,1,0,1);
	  cells[1547]->set_vertices(vertices[272],vertices[264],vertices[271],vertices[24]);
	  cells[1547]->set_neighbors(cells[1248],cells[1492],cells[1318],cells[1348]);
	  set_offsets(cells[1547],1,1,0,5);
	  cells[1548]->set_vertices(vertices[255],vertices[248],vertices[8],vertices[256]);
	  cells[1548]->set_neighbors(cells[1544],cells[1266],cells[1504],cells[1534]);
	  set_offsets(cells[1548],0,1,5,1);
	  cells[1549]->set_vertices(vertices[258],vertices[250],vertices[251],vertices[203]);
	  cells[1549]->set_neighbors(cells[1519],cells[391],cells[1285],cells[1277]);
	  set_offsets(cells[1549],0,0,0,0);
	  cells[1550]->set_vertices(vertices[265],vertices[26],vertices[25],vertices[34]);
	  cells[1550]->set_neighbors(cells[30],cells[1659],cells[1563],cells[1612]);
	  set_offsets(cells[1550],0,4,4,4);
	  cells[1551]->set_vertices(vertices[209],vertices[249],vertices[256],vertices[257]);
	  cells[1551]->set_neighbors(cells[1552],cells[1599],cells[1378],cells[1542]);
	  set_offsets(cells[1551],0,0,0,0);
	  cells[1552]->set_vertices(vertices[257],vertices[249],vertices[256],vertices[16]);
	  cells[1552]->set_neighbors(cells[1278],cells[1271],cells[1557],cells[1551]);
	  set_offsets(cells[1552],0,0,0,4);
	  cells[1553]->set_vertices(vertices[24],vertices[16],vertices[17],vertices[257]);
	  cells[1553]->set_neighbors(cells[1554],cells[1614],cells[1324],cells[87]);
	  set_offsets(cells[1553],4,4,4,0);
	  cells[1554]->set_vertices(vertices[16],vertices[257],vertices[9],vertices[17]);
	  cells[1554]->set_neighbors(cells[1555],cells[170],cells[1553],cells[1557]);
	  set_offsets(cells[1554],4,0,4,4);
	  cells[1555]->set_vertices(vertices[257],vertices[18],vertices[9],vertices[17]);
	  cells[1555]->set_neighbors(cells[232],cells[1554],cells[559],cells[1556]);
	  set_offsets(cells[1555],0,4,4,4);
	  cells[1556]->set_vertices(vertices[249],vertices[18],vertices[9],vertices[257]);
	  cells[1556]->set_neighbors(cells[1555],cells[1557],cells[145],cells[1006]);
	  set_offsets(cells[1556],0,4,4,0);
	  cells[1557]->set_vertices(vertices[16],vertices[249],vertices[9],vertices[257]);
	  cells[1557]->set_neighbors(cells[1556],cells[1554],cells[1552],cells[1505]);
	  set_offsets(cells[1557],4,0,4,0);
	  cells[1558]->set_vertices(vertices[266],vertices[257],vertices[18],vertices[258]);
	  cells[1558]->set_neighbors(cells[145],cells[214],cells[1559],cells[1560]);
	  set_offsets(cells[1558],0,0,4,0);
	  cells[1559]->set_vertices(vertices[266],vertices[209],vertices[257],vertices[258]);
	  cells[1559]->set_neighbors(cells[1378],cells[1558],cells[1110],cells[1514]);
	  set_offsets(cells[1559],0,0,0,0);
	  cells[1560]->set_vertices(vertices[26],vertices[257],vertices[18],vertices[266]);
	  cells[1560]->set_neighbors(cells[1558],cells[1121],cells[570],cells[559]);
	  set_offsets(cells[1560],4,0,4,0);
	  cells[1561]->set_vertices(vertices[209],vertices[201],vertices[249],vertices[258]);
	  cells[1561]->set_neighbors(cells[1562],cells[1378],cells[1288],cells[1542]);
	  set_offsets(cells[1561],0,0,0,0);
	  cells[1562]->set_vertices(vertices[258],vertices[201],vertices[249],vertices[250]);
	  cells[1562]->set_neighbors(cells[342],cells[1510],cells[876],cells[1561]);
	  set_offsets(cells[1562],0,0,0,0);
	  cells[1563]->set_vertices(vertices[34],vertices[265],vertices[26],vertices[274]);
	  cells[1563]->set_neighbors(cells[1385],cells[1621],cells[83],cells[1550]);
	  set_offsets(cells[1563],4,0,4,0);
	  cells[1564]->set_vertices(vertices[266],vertices[258],vertices[259],vertices[211]);
	  cells[1564]->set_neighbors(cells[1520],cells[958],cells[1335],cells[214]);
	  set_offsets(cells[1564],0,0,0,0);
	  cells[1565]->set_vertices(vertices[34],vertices[267],vertices[26],vertices[27]);
	  cells[1565]->set_neighbors(cells[1515],cells[186],cells[874],cells[1621]);
	  set_offsets(cells[1565],4,0,4,4);
	  cells[1566]->set_vertices(vertices[251],vertices[20],vertices[11],vertices[259]);
	  cells[1566]->set_neighbors(cells[1567],cells[1518],cells[1280],cells[183]);
	  set_offsets(cells[1566],0,4,4,0);
	  cells[1567]->set_vertices(vertices[259],vertices[20],vertices[11],vertices[19]);
	  cells[1567]->set_neighbors(cells[144],cells[1568],cells[1467],cells[1566]);
	  set_offsets(cells[1567],0,4,4,4);
	  cells[1568]->set_vertices(vertices[18],vertices[259],vertices[11],vertices[19]);
	  cells[1568]->set_neighbors(cells[1567],cells[423],cells[1569],cells[1518]);
	  set_offsets(cells[1568],4,0,4,4);
	  cells[1569]->set_vertices(vertices[26],vertices[18],vertices[19],vertices[259]);
	  cells[1569]->set_neighbors(cells[1568],cells[1370],cells[1121],cells[370]);
	  set_offsets(cells[1569],4,4,4,0);
	  cells[1570]->set_vertices(vertices[244],vertices[196],vertices[235],vertices[284]);
	  cells[1570]->set_neighbors(cells[1462],cells[1701],cells[1475],cells[1458]);
	  set_offsets(cells[1570],2,2,0,0);
	  cells[1571]->set_vertices(vertices[269],vertices[268],vertices[28],vertices[261]);
	  cells[1571]->set_neighbors(cells[128],cells[136],cells[1633],cells[1631]);
	  set_offsets(cells[1571],0,0,4,0);
	  cells[1572]->set_vertices(vertices[13],vertices[253],vertices[261],vertices[20]);
	  cells[1572]->set_neighbors(cells[1188],cells[1573],cells[1279],cells[1579]);
	  set_offsets(cells[1572],4,0,0,4);
	  cells[1573]->set_vertices(vertices[13],vertices[261],vertices[21],vertices[20]);
	  cells[1573]->set_neighbors(cells[1251],cells[434],cells[1572],cells[1580]);
	  set_offsets(cells[1573],4,0,4,4);
	  cells[1574]->set_vertices(vertices[267],vertices[36],vertices[28],vertices[27]);
	  cells[1574]->set_neighbors(cells[99],cells[1195],cells[1389],cells[1045]);
	  set_offsets(cells[1574],0,4,4,4);
	  cells[1575]->set_vertices(vertices[219],vertices[211],vertices[259],vertices[268]);
	  cells[1575]->set_neighbors(cells[1331],cells[1262],cells[1386],cells[958]);
	  set_offsets(cells[1575],0,0,0,0);
	  cells[1576]->set_vertices(vertices[12],vertices[251],vertices[252],vertices[260]);
	  cells[1576]->set_neighbors(cells[1020],cells[1577],cells[1253],cells[124]);
	  set_offsets(cells[1576],4,0,0,0);
	  cells[1577]->set_vertices(vertices[12],vertices[252],vertices[253],vertices[260]);
	  cells[1577]->set_neighbors(cells[1465],cells[1293],cells[1576],cells[1484]);
	  set_offsets(cells[1577],4,0,0,0);
	  cells[1578]->set_vertices(vertices[227],vertices[274],vertices[275],vertices[267]);
	  cells[1578]->set_neighbors(cells[1524],cells[1626],cells[1440],cells[1649]);
	  set_offsets(cells[1578],0,0,0,0);
	  cells[1579]->set_vertices(vertices[13],vertices[22],vertices[261],vertices[253]);
	  cells[1579]->set_neighbors(cells[1411],cells[1572],cells[1528],cells[1580]);
	  set_offsets(cells[1579],4,4,0,0);
	  cells[1580]->set_vertices(vertices[13],vertices[22],vertices[21],vertices[261]);
	  cells[1580]->set_neighbors(cells[1535],cells[1573],cells[1579],cells[178]);
	  set_offsets(cells[1580],4,4,4,0);
	  cells[1581]->set_vertices(vertices[38],vertices[278],vertices[30],vertices[271]);
	  cells[1581]->set_neighbors(cells[1589],cells[1395],cells[1245],cells[346]);
	  set_offsets(cells[1581],4,0,4,0);
	  cells[1582]->set_vertices(vertices[215],vertices[256],vertices[255],vertices[263]);
	  cells[1582]->set_neighbors(cells[1169],cells[1308],cells[1602],cells[1539]);
	  set_offsets(cells[1582],0,1,0,0);
	  cells[1583]->set_vertices(vertices[269],vertices[38],vertices[30],vertices[29]);
	  cells[1583]->set_neighbors(cells[26],cells[1584],cells[1181],cells[346]);
	  set_offsets(cells[1583],0,4,4,4);
	  cells[1584]->set_vertices(vertices[21],vertices[30],vertices[29],vertices[269]);
	  cells[1584]->set_neighbors(cells[1583],cells[1464],cells[1209],cells[404]);
	  set_offsets(cells[1584],4,4,4,0);
	  cells[1585]->set_vertices(vertices[270],vertices[213],vertices[261],vertices[262]);
	  cells[1585]->set_neighbors(cells[1587],cells[1533],cells[1057],cells[127]);
	  set_offsets(cells[1585],0,0,0,0);
	  cells[1586]->set_vertices(vertices[221],vertices[213],vertices[270],vertices[222]);
	  cells[1586]->set_neighbors(cells[1057],cells[673],cells[1103],cells[127]);
	  set_offsets(cells[1586],0,0,0,0);
	  cells[1587]->set_vertices(vertices[213],vertices[253],vertices[261],vertices[262]);
	  cells[1587]->set_neighbors(cells[1411],cells[1585],cells[1588],cells[1344]);
	  set_offsets(cells[1587],0,0,0,0);
	  cells[1588]->set_vertices(vertices[213],vertices[205],vertices[253],vertices[262]);
	  cells[1588]->set_neighbors(cells[1482],cells[1587],cells[1290],cells[1319]);
	  set_offsets(cells[1588],0,0,0,0);
	  cells[1589]->set_vertices(vertices[278],vertices[270],vertices[30],vertices[271]);
	  cells[1589]->set_neighbors(cells[1388],cells[1581],cells[1644],cells[911]);
	  set_offsets(cells[1589],0,0,4,0);
	  cells[1590]->set_vertices(vertices[271],vertices[263],vertices[23],vertices[24]);
	  cells[1590]->set_neighbors(cells[939],cells[1261],cells[1248],cells[1593]);
	  set_offsets(cells[1590],0,0,4,5);
	  cells[1591]->set_vertices(vertices[270],vertices[222],vertices[262],vertices[215]);
	  cells[1591]->set_neighbors(cells[1349],cells[779],cells[1592],cells[1057]);
	  set_offsets(cells[1591],0,0,0,0);
	  cells[1592]->set_vertices(vertices[223],vertices[222],vertices[270],vertices[215]);
	  cells[1592]->set_neighbors(cells[1591],cells[1296],cells[989],cells[1276]);
	  set_offsets(cells[1592],0,0,0,0);
	  cells[1593]->set_vertices(vertices[271],vertices[30],vertices[23],vertices[263]);
	  cells[1593]->set_neighbors(cells[974],cells[1590],cells[1388],cells[392]);
	  set_offsets(cells[1593],0,4,4,0);
	  cells[1594]->set_vertices(vertices[39],vertices[279],vertices[40],vertices[287]);
	  cells[1594]->set_neighbors(cells[1180],cells[1595],cells[1727],cells[1220]);
	  set_offsets(cells[1594],4,0,5,0);
	  cells[1595]->set_vertices(vertices[39],vertices[287],vertices[40],vertices[47]);
	  cells[1595]->set_neighbors(cells[1615],cells[550],cells[1021],cells[1594]);
	  set_offsets(cells[1595],4,0,5,4);
	  cells[1596]->set_vertices(vertices[272],vertices[232],vertices[280],vertices[231]);
	  cells[1596]->set_neighbors(cells[1597],cells[1292],cells[1410],cells[1487]);
	  set_offsets(cells[1596],1,1,1,0);
	  cells[1597]->set_vertices(vertices[280],vertices[232],vertices[239],vertices[231]);
	  cells[1597]->set_neighbors(cells[1437],cells[1247],cells[1596],cells[1394]);
	  set_offsets(cells[1597],1,1,0,0);
	  cells[1598]->set_vertices(vertices[272],vertices[225],vertices[232],vertices[224]);
	  cells[1598]->set_neighbors(cells[1396],cells[1410],cells[277],cells[1487]);
	  set_offsets(cells[1598],0,0,0,0);
	  cells[1599]->set_vertices(vertices[264],vertices[209],vertices[256],vertices[257]);
	  cells[1599]->set_neighbors(cells[1551],cells[1271],cells[1600],cells[1218]);
	  set_offsets(cells[1599],0,0,0,0);
	  cells[1600]->set_vertices(vertices[217],vertices[209],vertices[264],vertices[257]);
	  cells[1600]->set_neighbors(cells[1599],cells[1608],cells[1514],cells[1446]);
	  set_offsets(cells[1600],0,0,0,0);
	  cells[1601]->set_vertices(vertices[264],vertices[216],vertices[215],vertices[256]);
	  cells[1601]->set_neighbors(cells[1545],cells[1602],cells[1218],cells[1604]);
	  set_offsets(cells[1601],1,1,0,1);
	  cells[1602]->set_vertices(vertices[264],vertices[256],vertices[215],vertices[263]);
	  cells[1602]->set_neighbors(cells[1582],cells[1603],cells[1316],cells[1601]);
	  set_offsets(cells[1602],1,1,0,0);
	  cells[1603]->set_vertices(vertices[223],vertices[264],vertices[215],vertices[263]);
	  cells[1603]->set_neighbors(cells[1602],cells[1296],cells[1637],cells[1604]);
	  set_offsets(cells[1603],0,1,0,0);
	  cells[1604]->set_vertices(vertices[223],vertices[216],vertices[215],vertices[264]);
	  cells[1604]->set_neighbors(cells[1601],cells[1603],cells[1605],cells[837]);
	  set_offsets(cells[1604],0,1,0,1);
	  cells[1605]->set_vertices(vertices[224],vertices[216],vertices[223],vertices[264]);
	  cells[1605]->set_neighbors(cells[1604],cells[1170],cells[157],cells[1085]);
	  set_offsets(cells[1605],1,1,0,1);
	  cells[1606]->set_vertices(vertices[281],vertices[42],vertices[41],vertices[2]);
	  cells[1606]->set_neighbors(cells[52],cells[1231],cells[1232],cells[1183]);
	  set_offsets(cells[1606],0,4,4,6);
	  cells[1607]->set_vertices(vertices[280],vertices[272],vertices[279],vertices[32]);
	  cells[1607]->set_neighbors(cells[1477],cells[1541],cells[1346],cells[1292]);
	  set_offsets(cells[1607],1,1,0,5);
	  cells[1608]->set_vertices(vertices[217],vertices[257],vertices[264],vertices[265]);
	  cells[1608]->set_neighbors(cells[1609],cells[1647],cells[1059],cells[1600]);
	  set_offsets(cells[1608],0,0,0,0);
	  cells[1609]->set_vertices(vertices[265],vertices[257],vertices[264],vertices[24]);
	  cells[1609]->set_neighbors(cells[1324],cells[1318],cells[1614],cells[1608]);
	  set_offsets(cells[1609],0,0,0,4);
	  cells[1610]->set_vertices(vertices[32],vertices[24],vertices[25],vertices[265]);
	  cells[1610]->set_neighbors(cells[1611],cells[1658],cells[910],cells[92]);
	  set_offsets(cells[1610],4,4,4,0);
	  cells[1611]->set_vertices(vertices[24],vertices[265],vertices[17],vertices[25]);
	  cells[1611]->set_neighbors(cells[1612],cells[287],cells[1610],cells[1614]);
	  set_offsets(cells[1611],4,0,4,4);
	  cells[1612]->set_vertices(vertices[265],vertices[26],vertices[17],vertices[25]);
	  cells[1612]->set_neighbors(cells[242],cells[1611],cells[1550],cells[1613]);
	  set_offsets(cells[1612],0,4,4,4);
	  cells[1613]->set_vertices(vertices[257],vertices[26],vertices[17],vertices[265]);
	  cells[1613]->set_neighbors(cells[1612],cells[1614],cells[570],cells[559]);
	  set_offsets(cells[1613],0,4,4,0);
	  cells[1614]->set_vertices(vertices[24],vertices[257],vertices[17],vertices[265]);
	  cells[1614]->set_neighbors(cells[1613],cells[1611],cells[1609],cells[1553]);
	  set_offsets(cells[1614],4,0,4,0);
	  cells[1615]->set_vertices(vertices[287],vertices[40],vertices[47],vertices[0]);
	  cells[1615]->set_neighbors(cells[93],cells[1713],cells[1677],cells[1595]);
	  set_offsets(cells[1615],0,5,4,7);
	  cells[1616]->set_vertices(vertices[284],vertices[276],vertices[36],vertices[277]);
	  cells[1616]->set_neighbors(cells[1639],cells[1198],cells[1661],cells[1704]);
	  set_offsets(cells[1616],0,0,4,0);
	  cells[1617]->set_vertices(vertices[282],vertices[194],vertices[233],vertices[234]);
	  cells[1617]->set_neighbors(cells[1323],cells[344],cells[1694],cells[226]);
	  set_offsets(cells[1617],0,2,0,0);
	  cells[1618]->set_vertices(vertices[42],vertices[275],vertices[34],vertices[35]);
	  cells[1618]->set_neighbors(cells[1660],cells[497],cells[1341],cells[1328]);
	  set_offsets(cells[1618],4,0,4,4);
	  cells[1619]->set_vertices(vertices[2],vertices[242],vertices[283],vertices[243]);
	  cells[1619]->set_neighbors(cells[1380],cells[1445],cells[1009],cells[1699]);
	  set_offsets(cells[1619],6,2,0,2);
	  cells[1620]->set_vertices(vertices[274],vertices[266],vertices[26],vertices[267]);
	  cells[1620]->set_neighbors(cells[1176],cells[1621],cells[1152],cells[1385]);
	  set_offsets(cells[1620],0,0,4,0);
	  cells[1621]->set_vertices(vertices[34],vertices[274],vertices[26],vertices[267]);
	  cells[1621]->set_neighbors(cells[1620],cells[1565],cells[1524],cells[1563]);
	  set_offsets(cells[1621],4,0,4,0);
	  cells[1622]->set_vertices(vertices[219],vertices[266],vertices[267],vertices[259]);
	  cells[1622]->set_neighbors(cells[1176],cells[1262],cells[958],cells[1152]);
	  set_offsets(cells[1622],0,0,0,0);
	  cells[1623]->set_vertices(vertices[275],vertices[284],vertices[44],vertices[283]);
	  cells[1623]->set_neighbors(cells[324],cells[1666],cells[1667],cells[1705]);
	  set_offsets(cells[1623],0,0,4,0);
	  cells[1624]->set_vertices(vertices[276],vertices[219],vertices[267],vertices[268]);
	  cells[1624]->set_neighbors(cells[1262],cells[664],cells[1374],cells[1625]);
	  set_offsets(cells[1624],0,0,0,0);
	  cells[1625]->set_vertices(vertices[227],vertices[219],vertices[267],vertices[276]);
	  cells[1625]->set_neighbors(cells[1624],cells[1626],cells[326],cells[1440]);
	  set_offsets(cells[1625],0,0,0,0);
	  cells[1626]->set_vertices(vertices[227],vertices[267],vertices[275],vertices[276]);
	  cells[1626]->set_neighbors(cells[1382],cells[1703],cells[1625],cells[1578]);
	  set_offsets(cells[1626],0,0,0,0);
	  cells[1627]->set_vertices(vertices[238],vertices[237],vertices[286],vertices[229]);
	  cells[1627]->set_neighbors(cells[623],cells[1673],cells[347],cells[1641]);
	  set_offsets(cells[1627],0,0,0,0);
	  cells[1628]->set_vertices(vertices[285],vertices[6],vertices[46],vertices[45]);
	  cells[1628]->set_neighbors(cells[581],cells[1219],cells[1197],cells[1340]);
	  set_offsets(cells[1628],0,6,4,4);
	  cells[1629]->set_vertices(vertices[36],vertices[37],vertices[277],vertices[29]);
	  cells[1629]->set_neighbors(cells[1246],cells[1230],cells[119],cells[1193]);
	  set_offsets(cells[1629],4,4,0,4);
	  cells[1630]->set_vertices(vertices[4],vertices[244],vertices[285],vertices[245]);
	  cells[1630]->set_neighbors(cells[343],cells[1207],cells[1523],cells[1711]);
	  set_offsets(cells[1630],6,2,0,2);
	  cells[1631]->set_vertices(vertices[276],vertices[268],vertices[28],vertices[269]);
	  cells[1631]->set_neighbors(cells[1571],cells[1632],cells[1460],cells[664]);
	  set_offsets(cells[1631],0,0,4,0);
	  cells[1632]->set_vertices(vertices[36],vertices[276],vertices[28],vertices[269]);
	  cells[1632]->set_neighbors(cells[1631],cells[1527],cells[1639],cells[1045]);
	  set_offsets(cells[1632],4,0,4,0);
	  cells[1633]->set_vertices(vertices[221],vertices[268],vertices[269],vertices[261]);
	  cells[1633]->set_neighbors(cells[1571],cells[1536],cells[891],cells[1460]);
	  set_offsets(cells[1633],0,0,0,0);
	  cells[1634]->set_vertices(vertices[286],vertices[239],vertices[287],vertices[246]);
	  cells[1634]->set_neighbors(cells[1399],cells[1162],cells[1022],cells[1725]);
	  set_offsets(cells[1634],0,0,0,2);
	  cells[1635]->set_vertices(vertices[237],vertices[286],vertices[277],vertices[285]);
	  cells[1635]->set_neighbors(cells[1672],cells[1707],cells[1148],cells[623]);
	  set_offsets(cells[1635],0,0,0,0);
	  cells[1636]->set_vertices(vertices[278],vertices[230],vertices[270],vertices[223]);
	  cells[1636]->set_neighbors(cells[1276],cells[1644],cells[709],cells[1303]);
	  set_offsets(cells[1636],0,0,0,0);
	  cells[1637]->set_vertices(vertices[223],vertices[264],vertices[263],vertices[271]);
	  cells[1637]->set_neighbors(cells[1248],cells[1391],cells[1348],cells[1603]);
	  set_offsets(cells[1637],0,1,0,0);
	  cells[1638]->set_vertices(vertices[229],vertices[276],vertices[277],vertices[269]);
	  cells[1638]->set_neighbors(cells[1639],cells[1472],cells[1140],cells[1661]);
	  set_offsets(cells[1638],0,0,0,0);
	  cells[1639]->set_vertices(vertices[277],vertices[276],vertices[36],vertices[269]);
	  cells[1639]->set_neighbors(cells[1632],cells[1230],cells[1638],cells[1616]);
	  set_offsets(cells[1639],0,0,4,0);
	  cells[1640]->set_vertices(vertices[287],vertices[199],vertices[239],vertices[240]);
	  cells[1640]->set_neighbors(cells[1444],cells[1676],cells[1720],cells[1399]);
	  set_offsets(cells[1640],0,2,0,3);
	  cells[1641]->set_vertices(vertices[238],vertices[198],vertices[286],vertices[237]);
	  cells[1641]->set_neighbors(cells[1478],cells[1627],cells[1363],cells[1473]);
	  set_offsets(cells[1641],0,2,0,0);
	  cells[1642]->set_vertices(vertices[280],vertices[233],vertices[192],vertices[232]);
	  cells[1642]->set_neighbors(cells[1406],cells[1394],cells[125],cells[1052]);
	  set_offsets(cells[1642],0,0,2,0);
	  cells[1643]->set_vertices(vertices[231],vertices[278],vertices[271],vertices[223]);
	  cells[1643]->set_neighbors(cells[1644],cells[1299],cells[709],cells[936]);
	  set_offsets(cells[1643],0,0,0,0);
	  cells[1644]->set_vertices(vertices[278],vertices[270],vertices[271],vertices[223]);
	  cells[1644]->set_neighbors(cells[1391],cells[1643],cells[1636],cells[1589]);
	  set_offsets(cells[1644],0,0,0,0);
	  cells[1645]->set_vertices(vertices[40],vertices[273],vertices[280],vertices[32]);
	  cells[1645]->set_neighbors(cells[1346],cells[1541],cells[1439],cells[1651]);
	  set_offsets(cells[1645],4,0,0,4);
	  cells[1646]->set_vertices(vertices[225],vertices[217],vertices[272],vertices[265]);
	  cells[1646]->set_neighbors(cells[1647],cells[1656],cells[1330],cells[277]);
	  set_offsets(cells[1646],0,0,0,0);
	  cells[1647]->set_vertices(vertices[272],vertices[217],vertices[264],vertices[265]);
	  cells[1647]->set_neighbors(cells[1608],cells[1318],cells[1646],cells[1315]);
	  set_offsets(cells[1647],0,0,0,0);
	  cells[1648]->set_vertices(vertices[273],vertices[225],vertices[233],vertices[280]);
	  cells[1648]->set_neighbors(cells[125],cells[1652],cells[1433],cells[1693]);
	  set_offsets(cells[1648],0,0,0,0);
	  cells[1649]->set_vertices(vertices[275],vertices[282],vertices[274],vertices[227]);
	  cells[1649]->set_neighbors(cells[1650],cells[1578],cells[1452],cells[1680]);
	  set_offsets(cells[1649],0,0,0,0);
	  cells[1650]->set_vertices(vertices[274],vertices[282],vertices[234],vertices[227]);
	  cells[1650]->set_neighbors(cells[252],cells[1654],cells[1649],cells[1100]);
	  set_offsets(cells[1650],0,0,0,0);
	  cells[1651]->set_vertices(vertices[40],vertices[280],vertices[273],vertices[281]);
	  cells[1651]->set_neighbors(cells[1652],cells[586],cells[1688],cells[1645]);
	  set_offsets(cells[1651],4,0,0,0);
	  cells[1652]->set_vertices(vertices[273],vertices[280],vertices[233],vertices[281]);
	  cells[1652]->set_neighbors(cells[1681],cells[1329],cells[1651],cells[1648]);
	  set_offsets(cells[1652],0,0,0,0);
	  cells[1653]->set_vertices(vertices[41],vertices[40],vertices[33],vertices[281]);
	  cells[1653]->set_neighbors(cells[586],cells[1183],cells[1686],cells[507]);
	  set_offsets(cells[1653],4,4,4,0);
	  cells[1654]->set_vertices(vertices[227],vertices[234],vertices[274],vertices[226]);
	  cells[1654]->set_neighbors(cells[1655],cells[916],cells[1368],cells[1650]);
	  set_offsets(cells[1654],0,0,0,0);
	  cells[1655]->set_vertices(vertices[274],vertices[234],vertices[225],vertices[226]);
	  cells[1655]->set_neighbors(cells[889],cells[85],cells[1654],cells[1100]);
	  set_offsets(cells[1655],0,0,0,0);
	  cells[1656]->set_vertices(vertices[225],vertices[265],vertices[272],vertices[273]);
	  cells[1656]->set_neighbors(cells[1657],cells[1433],cells[930],cells[1646]);
	  set_offsets(cells[1656],0,0,0,0);
	  cells[1657]->set_vertices(vertices[273],vertices[265],vertices[272],vertices[32]);
	  cells[1657]->set_neighbors(cells[910],cells[1346],cells[1658],cells[1656]);
	  set_offsets(cells[1657],0,0,0,4);
	  cells[1658]->set_vertices(vertices[32],vertices[265],vertices[25],vertices[273]);
	  cells[1658]->set_neighbors(cells[1659],cells[1412],cells[1657],cells[1610]);
	  set_offsets(cells[1658],4,0,4,0);
	  cells[1659]->set_vertices(vertices[265],vertices[34],vertices[25],vertices[273]);
	  cells[1659]->set_neighbors(cells[948],cells[1658],cells[83],cells[1550]);
	  set_offsets(cells[1659],0,4,4,0);
	  cells[1660]->set_vertices(vertices[34],vertices[35],vertices[275],vertices[27]);
	  cells[1660]->set_neighbors(cells[1424],cells[874],cells[310],cells[1618]);
	  set_offsets(cells[1660],4,4,0,4);
	  cells[1661]->set_vertices(vertices[277],vertices[284],vertices[276],vertices[229]);
	  cells[1661]->set_neighbors(cells[1665],cells[1638],cells[1003],cells[1616]);
	  set_offsets(cells[1661],0,0,0,0);
	  cells[1662]->set_vertices(vertices[283],vertices[195],vertices[235],vertices[244]);
	  cells[1662]->set_neighbors(cells[1458],cells[1701],cells[1663],cells[1700]);
	  set_offsets(cells[1662],0,2,0,2);
	  cells[1663]->set_vertices(vertices[243],vertices[195],vertices[283],vertices[244]);
	  cells[1663]->set_neighbors(cells[1662],cells[1283],cells[1453],cells[1380]);
	  set_offsets(cells[1663],2,2,0,2);
	  cells[1664]->set_vertices(vertices[229],vertices[228],vertices[276],vertices[221]);
	  cells[1664]->set_neighbors(cells[1461],cells[1140],cells[1375],cells[714]);
	  set_offsets(cells[1664],0,0,0,0);
	  cells[1665]->set_vertices(vertices[276],vertices[284],vertices[236],vertices[229]);
	  cells[1665]->set_neighbors(cells[1463],cells[714],cells[1661],cells[1194]);
	  set_offsets(cells[1665],0,0,0,0);
	  cells[1666]->set_vertices(vertices[283],vertices[44],vertices[275],vertices[35]);
	  cells[1666]->set_neighbors(cells[267],cells[1341],cells[1004],cells[1623]);
	  set_offsets(cells[1666],0,4,0,4);
	  cells[1667]->set_vertices(vertices[235],vertices[284],vertices[275],vertices[283]);
	  cells[1667]->set_neighbors(cells[1623],cells[1695],cells[1701],cells[1702]);
	  set_offsets(cells[1667],0,0,0,0);
	  cells[1668]->set_vertices(vertices[285],vertices[197],vertices[237],vertices[246]);
	  cells[1668]->set_neighbors(cells[1206],cells[1148],cells[1669],cells[1712]);
	  set_offsets(cells[1668],0,2,0,2);
	  cells[1669]->set_vertices(vertices[245],vertices[197],vertices[285],vertices[246]);
	  cells[1669]->set_neighbors(cells[1668],cells[1056],cells[1026],cells[343]);
	  set_offsets(cells[1669],2,2,0,2);
	  cells[1670]->set_vertices(vertices[229],vertices[221],vertices[269],vertices[278]);
	  cells[1670]->set_neighbors(cells[1479],cells[1472],cells[1166],cells[1140]);
	  set_offsets(cells[1670],0,0,0,0);
	  cells[1671]->set_vertices(vertices[39],vertices[38],vertices[279],vertices[46]);
	  cells[1671]->set_neighbors(cells[1717],cells[1727],cells[315],cells[1015]);
	  set_offsets(cells[1671],4,4,0,4);
	  cells[1672]->set_vertices(vertices[277],vertices[286],vertices[46],vertices[285]);
	  cells[1672]->set_neighbors(cells[1208],cells[1172],cells[1635],cells[1718]);
	  set_offsets(cells[1672],0,0,4,0);
	  cells[1673]->set_vertices(vertices[238],vertices[286],vertices[278],vertices[229]);
	  cells[1673]->set_neighbors(cells[1480],cells[1119],cells[1627],cells[1476]);
	  set_offsets(cells[1673],0,0,0,0);
	  cells[1674]->set_vertices(vertices[231],vertices[239],vertices[286],vertices[238]);
	  cells[1674]->set_neighbors(cells[1473],cells[1476],cells[1434],cells[1714]);
	  set_offsets(cells[1674],0,0,0,0);
	  cells[1675]->set_vertices(vertices[230],vertices[238],vertices[231],vertices[278]);
	  cells[1675]->set_neighbors(cells[1476],cells[709],cells[1119],cells[972]);
	  set_offsets(cells[1675],0,0,0,0);
	  cells[1676]->set_vertices(vertices[287],vertices[240],vertices[239],vertices[280]);
	  cells[1676]->set_neighbors(cells[1678],cells[957],cells[1405],cells[1640]);
	  set_offsets(cells[1676],0,3,0,1);
	  cells[1677]->set_vertices(vertices[0],vertices[287],vertices[40],vertices[240]);
	  cells[1677]->set_neighbors(cells[1405],cells[1687],cells[1721],cells[1615]);
	  set_offsets(cells[1677],7,0,5,3);
	  cells[1678]->set_vertices(vertices[239],vertices[240],vertices[192],vertices[280]);
	  cells[1678]->set_neighbors(cells[1052],cells[1394],cells[1676],cells[1444]);
	  set_offsets(cells[1678],0,3,3,1);
	  cells[1679]->set_vertices(vertices[233],vertices[282],vertices[281],vertices[242]);
	  cells[1679]->set_neighbors(cells[291],cells[1499],cells[226],cells[1329]);
	  set_offsets(cells[1679],0,0,0,2);
	  cells[1680]->set_vertices(vertices[282],vertices[274],vertices[34],vertices[275]);
	  cells[1680]->set_neighbors(cells[1524],cells[1328],cells[1649],cells[1691]);
	  set_offsets(cells[1680],0,0,4,0);
	  cells[1681]->set_vertices(vertices[281],vertices[280],vertices[233],vertices[240]);
	  cells[1681]->set_neighbors(cells[1052],cells[1682],cells[1688],cells[1652]);
	  set_offsets(cells[1681],0,0,0,2);
	  cells[1682]->set_vertices(vertices[240],vertices[193],vertices[233],vertices[281]);
	  cells[1682]->set_neighbors(cells[1499],cells[1681],cells[1683],cells[1443]);
	  set_offsets(cells[1682],2,2,0,0);
	  cells[1683]->set_vertices(vertices[240],vertices[193],vertices[281],vertices[241]);
	  cells[1683]->set_neighbors(cells[1404],cells[1684],cells[301],cells[1682]);
	  set_offsets(cells[1683],2,2,0,2);
	  cells[1684]->set_vertices(vertices[0],vertices[240],vertices[281],vertices[241]);
	  cells[1684]->set_neighbors(cells[1683],cells[1685],cells[229],cells[1687]);
	  set_offsets(cells[1684],6,2,0,2);
	  cells[1685]->set_vertices(vertices[0],vertices[281],vertices[41],vertices[241]);
	  cells[1685]->set_neighbors(cells[1231],cells[1402],cells[1684],cells[1686]);
	  set_offsets(cells[1685],6,0,4,2);
	  cells[1686]->set_vertices(vertices[0],vertices[40],vertices[41],vertices[281]);
	  cells[1686]->set_neighbors(cells[1653],cells[1685],cells[1687],cells[235]);
	  set_offsets(cells[1686],6,4,4,0);
	  cells[1687]->set_vertices(vertices[0],vertices[240],vertices[40],vertices[281]);
	  cells[1687]->set_neighbors(cells[1688],cells[1686],cells[1684],cells[1677]);
	  set_offsets(cells[1687],6,2,4,0);
	  cells[1688]->set_vertices(vertices[40],vertices[280],vertices[281],vertices[240]);
	  cells[1688]->set_neighbors(cells[1681],cells[1687],cells[1405],cells[1651]);
	  set_offsets(cells[1688],4,0,0,2);
	  cells[1689]->set_vertices(vertices[273],vertices[42],vertices[281],vertices[282]);
	  cells[1689]->set_neighbors(cells[291],cells[1329],cells[1690],cells[628]);
	  set_offsets(cells[1689],0,4,0,0);
	  cells[1690]->set_vertices(vertices[42],vertices[273],vertices[34],vertices[282]);
	  cells[1690]->set_neighbors(cells[1691],cells[1328],cells[1689],cells[1438]);
	  set_offsets(cells[1690],4,0,4,0);
	  cells[1691]->set_vertices(vertices[282],vertices[273],vertices[34],vertices[274]);
	  cells[1691]->set_neighbors(cells[83],cells[1680],cells[1692],cells[1690]);
	  set_offsets(cells[1691],0,0,4,0);
	  cells[1692]->set_vertices(vertices[273],vertices[282],vertices[225],vertices[274]);
	  cells[1692]->set_neighbors(cells[1100],cells[930],cells[1691],cells[1693]);
	  set_offsets(cells[1692],0,0,0,0);
	  cells[1693]->set_vertices(vertices[273],vertices[233],vertices[225],vertices[282]);
	  cells[1693]->set_neighbors(cells[344],cells[1692],cells[1329],cells[1648]);
	  set_offsets(cells[1693],0,0,0,0);
	  cells[1694]->set_vertices(vertices[235],vertices[194],vertices[282],vertices[234]);
	  cells[1694]->set_neighbors(cells[1617],cells[252],cells[1419],cells[84]);
	  set_offsets(cells[1694],0,2,0,0);
	  cells[1695]->set_vertices(vertices[275],vertices[235],vertices[283],vertices[282]);
	  cells[1695]->set_neighbors(cells[1696],cells[1698],cells[1452],cells[1667]);
	  set_offsets(cells[1695],0,0,0,0);
	  cells[1696]->set_vertices(vertices[283],vertices[235],vertices[242],vertices[282]);
	  cells[1696]->set_neighbors(cells[84],cells[1697],cells[1695],cells[1700]);
	  set_offsets(cells[1696],0,0,2,0);
	  cells[1697]->set_vertices(vertices[42],vertices[283],vertices[242],vertices[282]);
	  cells[1697]->set_neighbors(cells[1696],cells[291],cells[1698],cells[1699]);
	  set_offsets(cells[1697],4,0,2,0);
	  cells[1698]->set_vertices(vertices[42],vertices[275],vertices[283],vertices[282]);
	  cells[1698]->set_neighbors(cells[1695],cells[1697],cells[1328],cells[1341]);
	  set_offsets(cells[1698],4,0,0,0);
	  cells[1699]->set_vertices(vertices[2],vertices[242],vertices[42],vertices[283]);
	  cells[1699]->set_neighbors(cells[1697],cells[931],cells[1619],cells[1232]);
	  set_offsets(cells[1699],6,2,4,0);
	  cells[1700]->set_vertices(vertices[242],vertices[195],vertices[235],vertices[283]);
	  cells[1700]->set_neighbors(cells[1662],cells[1696],cells[1380],cells[1002]);
	  set_offsets(cells[1700],2,2,0,0);
	  cells[1701]->set_vertices(vertices[235],vertices[244],vertices[284],vertices[283]);
	  cells[1701]->set_neighbors(cells[324],cells[1667],cells[1662],cells[1570]);
	  set_offsets(cells[1701],0,2,0,0);
	  cells[1702]->set_vertices(vertices[275],vertices[235],vertices[227],vertices[284]);
	  cells[1702]->set_neighbors(cells[1459],cells[1703],cells[1667],cells[1452]);
	  set_offsets(cells[1702],0,0,0,0);
	  cells[1703]->set_vertices(vertices[275],vertices[284],vertices[227],vertices[276]);
	  cells[1703]->set_neighbors(cells[1194],cells[1626],cells[1704],cells[1702]);
	  set_offsets(cells[1703],0,0,0,0);
	  cells[1704]->set_vertices(vertices[284],vertices[275],vertices[36],vertices[276]);
	  cells[1704]->set_neighbors(cells[1382],cells[1616],cells[1703],cells[1705]);
	  set_offsets(cells[1704],0,0,4,0);
	  cells[1705]->set_vertices(vertices[44],vertices[275],vertices[36],vertices[284]);
	  cells[1705]->set_neighbors(cells[1704],cells[1198],cells[1623],cells[267]);
	  set_offsets(cells[1705],4,0,4,0);
	  cells[1706]->set_vertices(vertices[237],vertices[196],vertices[284],vertices[236]);
	  cells[1706]->set_neighbors(cells[1462],cells[1463],cells[238],cells[1475]);
	  set_offsets(cells[1706],0,2,0,0);
	  cells[1707]->set_vertices(vertices[277],vertices[237],vertices[285],vertices[284]);
	  cells[1707]->set_neighbors(cells[1708],cells[1710],cells[1003],cells[1635]);
	  set_offsets(cells[1707],0,0,0,0);
	  cells[1708]->set_vertices(vertices[285],vertices[237],vertices[244],vertices[284]);
	  cells[1708]->set_neighbors(cells[1475],cells[1709],cells[1707],cells[1712]);
	  set_offsets(cells[1708],0,0,2,0);
	  cells[1709]->set_vertices(vertices[44],vertices[285],vertices[244],vertices[284]);
	  cells[1709]->set_neighbors(cells[1708],cells[324],cells[1710],cells[1711]);
	  set_offsets(cells[1709],4,0,2,0);
	  cells[1710]->set_vertices(vertices[44],vertices[277],vertices[285],vertices[284]);
	  cells[1710]->set_neighbors(cells[1707],cells[1709],cells[1198],cells[58]);
	  set_offsets(cells[1710],4,0,0,0);
	  cells[1711]->set_vertices(vertices[4],vertices[244],vertices[44],vertices[285]);
	  cells[1711]->set_neighbors(cells[1709],cells[1471],cells[1630],cells[237]);
	  set_offsets(cells[1711],6,2,4,0);
	  cells[1712]->set_vertices(vertices[244],vertices[197],vertices[237],vertices[285]);
	  cells[1712]->set_neighbors(cells[1668],cells[1708],cells[343],cells[1426]);
	  set_offsets(cells[1712],2,2,0,0);
	  cells[1713]->set_vertices(vertices[247],vertices[287],vertices[47],vertices[0]);
	  cells[1713]->set_neighbors(cells[1615],cells[1353],cells[1721],cells[1724]);
	  set_offsets(cells[1713],2,0,4,7);
	  cells[1714]->set_vertices(vertices[231],vertices[239],vertices[279],vertices[286]);
	  cells[1714]->set_neighbors(cells[1725],cells[1715],cells[1674],cells[1247]);
	  set_offsets(cells[1714],0,0,0,0);
	  cells[1715]->set_vertices(vertices[231],vertices[286],vertices[279],vertices[278]);
	  cells[1715]->set_neighbors(cells[1716],cells[936],cells[1476],cells[1714]);
	  set_offsets(cells[1715],0,0,0,0);
	  cells[1716]->set_vertices(vertices[286],vertices[278],vertices[38],vertices[279]);
	  cells[1716]->set_neighbors(cells[1245],cells[1717],cells[1715],cells[1719]);
	  set_offsets(cells[1716],0,0,4,0);
	  cells[1717]->set_vertices(vertices[46],vertices[286],vertices[38],vertices[279]);
	  cells[1717]->set_neighbors(cells[1716],cells[1671],cells[1726],cells[1718]);
	  set_offsets(cells[1717],4,0,4,0);
	  cells[1718]->set_vertices(vertices[46],vertices[277],vertices[38],vertices[286]);
	  cells[1718]->set_neighbors(cells[1719],cells[1717],cells[1672],cells[1474]);
	  set_offsets(cells[1718],4,0,4,0);
	  cells[1719]->set_vertices(vertices[286],vertices[277],vertices[38],vertices[278]);
	  cells[1719]->set_neighbors(cells[322],cells[1716],cells[1480],cells[1718]);
	  set_offsets(cells[1719],0,0,4,0);
	  cells[1720]->set_vertices(vertices[247],vertices[199],vertices[287],vertices[240]);
	  cells[1720]->set_neighbors(cells[1640],cells[1721],cells[1497],cells[1722]);
	  set_offsets(cells[1720],2,2,0,3);
	  cells[1721]->set_vertices(vertices[0],vertices[247],vertices[287],vertices[240]);
	  cells[1721]->set_neighbors(cells[1720],cells[1677],cells[1496],cells[1713]);
	  set_offsets(cells[1721],7,2,0,3);
	  cells[1722]->set_vertices(vertices[246],vertices[199],vertices[287],vertices[247]);
	  cells[1722]->set_neighbors(cells[1720],cells[1723],cells[769],cells[1399]);
	  set_offsets(cells[1722],2,2,0,2);
	  cells[1723]->set_vertices(vertices[6],vertices[246],vertices[287],vertices[247]);
	  cells[1723]->set_neighbors(cells[1722],cells[1724],cells[1485],cells[1400]);
	  set_offsets(cells[1723],6,2,0,2);
	  cells[1724]->set_vertices(vertices[47],vertices[287],vertices[247],vertices[6]);
	  cells[1724]->set_neighbors(cells[1723],cells[234],cells[1488],cells[1713]);
	  set_offsets(cells[1724],4,0,2,6);
	  cells[1725]->set_vertices(vertices[286],vertices[239],vertices[279],vertices[287]);
	  cells[1725]->set_neighbors(cells[957],cells[1726],cells[1634],cells[1714]);
	  set_offsets(cells[1725],0,0,0,0);
	  cells[1726]->set_vertices(vertices[286],vertices[279],vertices[46],vertices[287]);
	  cells[1726]->set_neighbors(cells[1727],cells[1162],cells[1725],cells[1717]);
	  set_offsets(cells[1726],0,0,4,0);
	  cells[1727]->set_vertices(vertices[39],vertices[46],vertices[279],vertices[287]);
	  cells[1727]->set_neighbors(cells[1726],cells[1594],cells[1021],cells[1671]);
	  set_offsets(cells[1727],4,4,0,0);
	  vertices[0]->set_cell(cells[1721]);
	  vertices[1]->set_cell(cells[1507]);
	  vertices[2]->set_cell(cells[1699]);
	  vertices[3]->set_cell(cells[1468]);
	  vertices[4]->set_cell(cells[1711]);
	  vertices[5]->set_cell(cells[1532]);
	  vertices[6]->set_cell(cells[1724]);
	  vertices[7]->set_cell(cells[1495]);
	  vertices[8]->set_cell(cells[1548]);
	  vertices[9]->set_cell(cells[1557]);
	  vertices[10]->set_cell(cells[1510]);
	  vertices[11]->set_cell(cells[1568]);
	  vertices[12]->set_cell(cells[1577]);
	  vertices[13]->set_cell(cells[1580]);
	  vertices[14]->set_cell(cells[1530]);
	  vertices[15]->set_cell(cells[1537]);
	  vertices[16]->set_cell(cells[1557]);
	  vertices[17]->set_cell(cells[1614]);
	  vertices[18]->set_cell(cells[1569]);
	  vertices[19]->set_cell(cells[1569]);
	  vertices[20]->set_cell(cells[1573]);
	  vertices[21]->set_cell(cells[1584]);
	  vertices[22]->set_cell(cells[1580]);
	  vertices[23]->set_cell(cells[1593]);
	  vertices[24]->set_cell(cells[1614]);
	  vertices[25]->set_cell(cells[1659]);
	  vertices[26]->set_cell(cells[1621]);
	  vertices[27]->set_cell(cells[1660]);
	  vertices[28]->set_cell(cells[1632]);
	  vertices[29]->set_cell(cells[1629]);
	  vertices[30]->set_cell(cells[1593]);
	  vertices[31]->set_cell(cells[1489]);
	  vertices[32]->set_cell(cells[1658]);
	  vertices[33]->set_cell(cells[1653]);
	  vertices[34]->set_cell(cells[1691]);
	  vertices[35]->set_cell(cells[1666]);
	  vertices[36]->set_cell(cells[1705]);
	  vertices[37]->set_cell(cells[1629]);
	  vertices[38]->set_cell(cells[1719]);
	  vertices[39]->set_cell(cells[1727]);
	  vertices[40]->set_cell(cells[1688]);
	  vertices[41]->set_cell(cells[1686]);
	  vertices[42]->set_cell(cells[1699]);
	  vertices[43]->set_cell(cells[1456]);
	  vertices[44]->set_cell(cells[1711]);
	  vertices[45]->set_cell(cells[1628]);
	  vertices[46]->set_cell(cells[1727]);
	  vertices[47]->set_cell(cells[1724]);
	  vertices[48]->set_cell(cells[839]);
	  vertices[49]->set_cell(cells[861]);
	  vertices[50]->set_cell(cells[855]);
	  vertices[51]->set_cell(cells[872]);
	  vertices[52]->set_cell(cells[788]);
	  vertices[53]->set_cell(cells[888]);
	  vertices[54]->set_cell(cells[881]);
	  vertices[55]->set_cell(cells[903]);
	  vertices[56]->set_cell(cells[637]);
	  vertices[57]->set_cell(cells[683]);
	  vertices[58]->set_cell(cells[681]);
	  vertices[59]->set_cell(cells[681]);
	  vertices[60]->set_cell(cells[658]);
	  vertices[61]->set_cell(cells[722]);
	  vertices[62]->set_cell(cells[677]);
	  vertices[63]->set_cell(cells[726]);
	  vertices[64]->set_cell(cells[726]);
	  vertices[65]->set_cell(cells[738]);
	  vertices[66]->set_cell(cells[679]);
	  vertices[67]->set_cell(cells[742]);
	  vertices[68]->set_cell(cells[549]);
	  vertices[69]->set_cell(cells[752]);
	  vertices[70]->set_cell(cells[723]);
	  vertices[71]->set_cell(cells[755]);
	  vertices[72]->set_cell(cells[738]);
	  vertices[73]->set_cell(cells[778]);
	  vertices[74]->set_cell(cells[697]);
	  vertices[75]->set_cell(cells[814]);
	  vertices[76]->set_cell(cells[666]);
	  vertices[77]->set_cell(cells[799]);
	  vertices[78]->set_cell(cells[723]);
	  vertices[79]->set_cell(cells[813]);
	  vertices[80]->set_cell(cells[813]);
	  vertices[81]->set_cell(cells[794]);
	  vertices[82]->set_cell(cells[814]);
	  vertices[83]->set_cell(cells[835]);
	  vertices[84]->set_cell(cells[800]);
	  vertices[85]->set_cell(cells[880]);
	  vertices[86]->set_cell(cells[777]);
	  vertices[87]->set_cell(cells[850]);
	  vertices[88]->set_cell(cells[806]);
	  vertices[89]->set_cell(cells[859]);
	  vertices[90]->set_cell(cells[821]);
	  vertices[91]->set_cell(cells[836]);
	  vertices[92]->set_cell(cells[871]);
	  vertices[93]->set_cell(cells[886]);
	  vertices[94]->set_cell(cells[881]);
	  vertices[95]->set_cell(cells[902]);
	  vertices[96]->set_cell(cells[1120]);
	  vertices[97]->set_cell(cells[1130]);
	  vertices[98]->set_cell(cells[1124]);
	  vertices[99]->set_cell(cells[1137]);
	  vertices[100]->set_cell(cells[1135]);
	  vertices[101]->set_cell(cells[1147]);
	  vertices[102]->set_cell(cells[1118]);
	  vertices[103]->set_cell(cells[1158]);
	  vertices[104]->set_cell(cells[908]);
	  vertices[105]->set_cell(cells[969]);
	  vertices[106]->set_cell(cells[950]);
	  vertices[107]->set_cell(cells[961]);
	  vertices[108]->set_cell(cells[907]);
	  vertices[109]->set_cell(cells[971]);
	  vertices[110]->set_cell(cells[932]);
	  vertices[111]->set_cell(cells[898]);
	  vertices[112]->set_cell(cells[978]);
	  vertices[113]->set_cell(cells[992]);
	  vertices[114]->set_cell(cells[950]);
	  vertices[115]->set_cell(cells[1001]);
	  vertices[116]->set_cell(cells[970]);
	  vertices[117]->set_cell(cells[1042]);
	  vertices[118]->set_cell(cells[1013]);
	  vertices[119]->set_cell(cells[1037]);
	  vertices[120]->set_cell(cells[1038]);
	  vertices[121]->set_cell(cells[1092]);
	  vertices[122]->set_cell(cells[924]);
	  vertices[123]->set_cell(cells[1094]);
	  vertices[124]->set_cell(cells[984]);
	  vertices[125]->set_cell(cells[1070]);
	  vertices[126]->set_cell(cells[1013]);
	  vertices[127]->set_cell(cells[1078]);
	  vertices[128]->set_cell(cells[1073]);
	  vertices[129]->set_cell(cells[1092]);
	  vertices[130]->set_cell(cells[1065]);
	  vertices[131]->set_cell(cells[1127]);
	  vertices[132]->set_cell(cells[1069]);
	  vertices[133]->set_cell(cells[1123]);
	  vertices[134]->set_cell(cells[1106]);
	  vertices[135]->set_cell(cells[1117]);
	  vertices[136]->set_cell(cells[1074]);
	  vertices[137]->set_cell(cells[1129]);
	  vertices[138]->set_cell(cells[1127]);
	  vertices[139]->set_cell(cells[1138]);
	  vertices[140]->set_cell(cells[1123]);
	  vertices[141]->set_cell(cells[1145]);
	  vertices[142]->set_cell(cells[1118]);
	  vertices[143]->set_cell(cells[1158]);
	  vertices[144]->set_cell(cells[1397]);
	  vertices[145]->set_cell(cells[1409]);
	  vertices[146]->set_cell(cells[1414]);
	  vertices[147]->set_cell(cells[1421]);
	  vertices[148]->set_cell(cells[1383]);
	  vertices[149]->set_cell(cells[1429]);
	  vertices[150]->set_cell(cells[1431]);
	  vertices[151]->set_cell(cells[1436]);
	  vertices[152]->set_cell(cells[1214]);
	  vertices[153]->set_cell(cells[1229]);
	  vertices[154]->set_cell(cells[1191]);
	  vertices[155]->set_cell(cells[1281]);
	  vertices[156]->set_cell(cells[1192]);
	  vertices[157]->set_cell(cells[1294]);
	  vertices[158]->set_cell(cells[1111]);
	  vertices[159]->set_cell(cells[1268]);
	  vertices[160]->set_cell(cells[1269]);
	  vertices[161]->set_cell(cells[1274]);
	  vertices[162]->set_cell(cells[1179]);
	  vertices[163]->set_cell(cells[1322]);
	  vertices[164]->set_cell(cells[1258]);
	  vertices[165]->set_cell(cells[1338]);
	  vertices[166]->set_cell(cells[1264]);
	  vertices[167]->set_cell(cells[1313]);
	  vertices[168]->set_cell(cells[1314]);
	  vertices[169]->set_cell(cells[1358]);
	  vertices[170]->set_cell(cells[1289]);
	  vertices[171]->set_cell(cells[1337]);
	  vertices[172]->set_cell(cells[1302]);
	  vertices[173]->set_cell(cells[1345]);
	  vertices[174]->set_cell(cells[1311]);
	  vertices[175]->set_cell(cells[1360]);
	  vertices[176]->set_cell(cells[1360]);
	  vertices[177]->set_cell(cells[1416]);
	  vertices[178]->set_cell(cells[1337]);
	  vertices[179]->set_cell(cells[1392]);
	  vertices[180]->set_cell(cells[1345]);
	  vertices[181]->set_cell(cells[1375]);
	  vertices[182]->set_cell(cells[1357]);
	  vertices[183]->set_cell(cells[1401]);
	  vertices[184]->set_cell(cells[1398]);
	  vertices[185]->set_cell(cells[1416]);
	  vertices[186]->set_cell(cells[1384]);
	  vertices[187]->set_cell(cells[1421]);
	  vertices[188]->set_cell(cells[1392]);
	  vertices[189]->set_cell(cells[1430]);
	  vertices[190]->set_cell(cells[1167]);
	  vertices[191]->set_cell(cells[1437]);
	  vertices[192]->set_cell(cells[1678]);
	  vertices[193]->set_cell(cells[1683]);
	  vertices[194]->set_cell(cells[1694]);
	  vertices[195]->set_cell(cells[1700]);
	  vertices[196]->set_cell(cells[1706]);
	  vertices[197]->set_cell(cells[1712]);
	  vertices[198]->set_cell(cells[1641]);
	  vertices[199]->set_cell(cells[1722]);
	  vertices[200]->set_cell(cells[1498]);
	  vertices[201]->set_cell(cells[1562]);
	  vertices[202]->set_cell(cells[1513]);
	  vertices[203]->set_cell(cells[1549]);
	  vertices[204]->set_cell(cells[1521]);
	  vertices[205]->set_cell(cells[1588]);
	  vertices[206]->set_cell(cells[1493]);
	  vertices[207]->set_cell(cells[1546]);
	  vertices[208]->set_cell(cells[1546]);
	  vertices[209]->set_cell(cells[1600]);
	  vertices[210]->set_cell(cells[1517]);
	  vertices[211]->set_cell(cells[1575]);
	  vertices[212]->set_cell(cells[1336]);
	  vertices[213]->set_cell(cells[1588]);
	  vertices[214]->set_cell(cells[1493]);
	  vertices[215]->set_cell(cells[1604]);
	  vertices[216]->set_cell(cells[1605]);
	  vertices[217]->set_cell(cells[1647]);
	  vertices[218]->set_cell(cells[1517]);
	  vertices[219]->set_cell(cells[1625]);
	  vertices[220]->set_cell(cells[1386]);
	  vertices[221]->set_cell(cells[1670]);
	  vertices[222]->set_cell(cells[1592]);
	  vertices[223]->set_cell(cells[1644]);
	  vertices[224]->set_cell(cells[1605]);
	  vertices[225]->set_cell(cells[1693]);
	  vertices[226]->set_cell(cells[1655]);
	  vertices[227]->set_cell(cells[1703]);
	  vertices[228]->set_cell(cells[1664]);
	  vertices[229]->set_cell(cells[1673]);
	  vertices[230]->set_cell(cells[1675]);
	  vertices[231]->set_cell(cells[1715]);
	  vertices[232]->set_cell(cells[1642]);
	  vertices[233]->set_cell(cells[1693]);
	  vertices[234]->set_cell(cells[1694]);
	  vertices[235]->set_cell(cells[1702]);
	  vertices[236]->set_cell(cells[1706]);
	  vertices[237]->set_cell(cells[1712]);
	  vertices[238]->set_cell(cells[1675]);
	  vertices[239]->set_cell(cells[1725]);
	  vertices[240]->set_cell(cells[1721]);
	  vertices[241]->set_cell(cells[1685]);
	  vertices[242]->set_cell(cells[1700]);
	  vertices[243]->set_cell(cells[1663]);
	  vertices[244]->set_cell(cells[1712]);
	  vertices[245]->set_cell(cells[1669]);
	  vertices[246]->set_cell(cells[1723]);
	  vertices[247]->set_cell(cells[1724]);
	  vertices[248]->set_cell(cells[1548]);
	  vertices[249]->set_cell(cells[1562]);
	  vertices[250]->set_cell(cells[1562]);
	  vertices[251]->set_cell(cells[1576]);
	  vertices[252]->set_cell(cells[1577]);
	  vertices[253]->set_cell(cells[1588]);
	  vertices[254]->set_cell(cells[1538]);
	  vertices[255]->set_cell(cells[1582]);
	  vertices[256]->set_cell(cells[1602]);
	  vertices[257]->set_cell(cells[1614]);
	  vertices[258]->set_cell(cells[1564]);
	  vertices[259]->set_cell(cells[1622]);
	  vertices[260]->set_cell(cells[1577]);
	  vertices[261]->set_cell(cells[1633]);
	  vertices[262]->set_cell(cells[1591]);
	  vertices[263]->set_cell(cells[1637]);
	  vertices[264]->set_cell(cells[1647]);
	  vertices[265]->set_cell(cells[1659]);
	  vertices[266]->set_cell(cells[1622]);
	  vertices[267]->set_cell(cells[1626]);
	  vertices[268]->set_cell(cells[1633]);
	  vertices[269]->set_cell(cells[1670]);
	  vertices[270]->set_cell(cells[1644]);
	  vertices[271]->set_cell(cells[1644]);
	  vertices[272]->set_cell(cells[1657]);
	  vertices[273]->set_cell(cells[1693]);
	  vertices[274]->set_cell(cells[1692]);
	  vertices[275]->set_cell(cells[1705]);
	  vertices[276]->set_cell(cells[1704]);
	  vertices[277]->set_cell(cells[1719]);
	  vertices[278]->set_cell(cells[1719]);
	  vertices[279]->set_cell(cells[1727]);
	  vertices[280]->set_cell(cells[1688]);
	  vertices[281]->set_cell(cells[1689]);
	  vertices[282]->set_cell(cells[1698]);
	  vertices[283]->set_cell(cells[1701]);
	  vertices[284]->set_cell(cells[1710]);
	  vertices[285]->set_cell(cells[1712]);
	  vertices[286]->set_cell(cells[1726]);
	  vertices[287]->set_cell(cells[1727]);

	  tds().set_dimension(3);
	  this->set_cover(make_array(1,1,1));

	  return vertices;
	}

  Oriented_side
  power_test(const Weighted_point &p, const Weighted_point &q) const
  {
      CGAL_triangulation_precondition(this->equal(p, q));
      return geom_traits().power_test_3_object()(p, q);
  }

  Oriented_side
  power_test(const Weighted_point &p, const Weighted_point &q,
       const Weighted_point &r, const Weighted_point &s,
       const Weighted_point &t, const Offset &o_p,
       const Offset &o_q, const Offset &o_r, const Offset &o_s,
       const Offset &o_t) const
  {
      return geom_traits().power_test_3_object()(p, q, r, s, t, o_p, o_q, o_r, o_s, o_t);
  }

  Oriented_side side_of_oriented_power_sphere(const Weighted_point &p, const Weighted_point &q,
      const Weighted_point &r, const Weighted_point &s, const Weighted_point &t, const Offset &o_p,
      const Offset &o_q, const Offset &o_r, const Offset &o_s,
      const Offset &o_t) const
  {
    return power_test(p,q,r,s,t,o_p,o_q,o_r,o_s,o_t);
  }

  Bounded_side side_of_power_sphere(const Cell_handle& c, const Weighted_point& p,
      const Offset & offset = Offset(), bool perturb = false) const{
    Bounded_side bs = ON_UNBOUNDED_SIDE;
    int i=0;
    // TODO: optimize which copies to check depending on the offsets in
    // the cell.
    while (bs == ON_UNBOUNDED_SIDE && i<8) {
      bs= _side_of_power_sphere(c,p,combine_offsets(offset,int_to_off(i)),perturb);
      i++;
    }
    return bs;
  }

  Vertex_handle nearest_power_vertex(const Weighted_point& p, Cell_handle start) const
  {
    if (number_of_vertices() == 0)
      return Vertex_handle();

    Locate_type lt;
    int li, lj;
    Cell_handle c = locate(p, lt, li, lj, start);
    if (lt == Base::VERTEX)
      return c->vertex(li);
    const Conflict_tester tester(p, this);
    Offset o = combine_offsets(Offset(), get_location_offset(tester, c));

    // - start with the closest vertex from the located cell.
    // - repeatedly take the nearest of its incident vertices if any
    // - if not, we're done.
    Vertex_handle nearest = nearest_vertex_in_cell(c, p, o);
    std::vector<Vertex_handle> vs;
    vs.reserve(32);
    while (true)
    {
      Vertex_handle tmp = nearest;
      Offset tmp_off = get_min_dist_offset(p, o, tmp);
      adjacent_vertices(nearest, std::back_inserter(vs));
      for (typename std::vector<Vertex_handle>::const_iterator vsit = vs.begin(); vsit != vs.end(); ++vsit)
        tmp =
            (compare_distance(p, tmp->point(), (*vsit)->point(), o, tmp_off, get_min_dist_offset(p, o, *vsit)) == SMALLER) ?
                tmp : *vsit;
      if (tmp == nearest)
        break;
      vs.clear();
      nearest = tmp;
    }

    return get_original_vertex(nearest);
  }

  Offset get_min_dist_offset(const Weighted_point & p, const Offset & o, const Vertex_handle vh) const {
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

  Vertex_handle nearest_vertex_in_cell(const Cell_handle& c, const Weighted_point & p, const Offset & o) const {
    CGAL_triangulation_precondition(number_of_vertices() != 0);
    Vertex_handle nearest = c->vertex(0);
    for (int i=1 ; i<4 ; i++) {
      nearest = (compare_distance(p,nearest->point(),c->vertex(i)->point(),
        o,get_offset(c,c->index(nearest)),get_offset(c,i)) == SMALLER) ?
        nearest : c->vertex(i);
    }
    return nearest;
  }

  Comparison_result compare_distance(const Weighted_point &p, const Weighted_point &q,
      const Weighted_point &r) const {
      return geom_traits().compare_power_distance_3_object()(p, q, r);
  }
  Comparison_result compare_distance(const Weighted_point& p, const Weighted_point& q,
      const Weighted_point& r, const Offset &o_p, const Offset &o_q,
      const Offset &o_r) const {
    return geom_traits().compare_power_distance_3_object()(p, q, r, o_p, o_q, o_r);
  }

  Bounded_side _side_of_power_sphere(const Cell_handle& c, const Weighted_point& p,
      const Offset & offset = Offset(), bool perturb = false) const;

  size_type number_of_hidden_points () const
  {
    size_type count = 0;
    for (Cell_iterator iter = cells_begin(), end_iter = cells_end(); iter != end_iter; ++iter)
      count += std::distance(iter->hidden_points_begin(), iter->hidden_points_end());
    return count;
  }

  bool is_valid(bool verbose = false, int level = 0) const;
  bool is_valid(Cell_handle c, bool verbose = false, int level = 0) const;

protected:
  // Protected, because inheritors(e.g. periodic triangulation for meshing)
  // of the class Periodic_3_Delaunay_triangulation_3 use this class
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
  Periodic_point periodic_weighted_circumcenter(Cell_handle c) const {
    return Base::periodic_circumcenter(c, geom_traits().construct_weighted_circumcenter_3_object());
  }

  /** @name Voronoi diagram */ //@{
  Bare_point dual(Cell_handle c) const {
    return point(periodic_weighted_circumcenter(c));
  }

  bool canonical_dual_segment(Cell_handle c, int i, Periodic_segment& ps) const {
    return Base::canonical_dual_segment(c, i, ps, geom_traits().construct_weighted_circumcenter_3_object());
  }

  Periodic_segment dual(const Facet & f) const {
    return dual( f.first, f.second );
  }
  Periodic_segment dual(Cell_handle c, int i) const{
    Periodic_segment ps;
    canonical_dual_segment(c,i,ps);
    return ps;
  }

  template <class OutputIterator>
  OutputIterator dual(const Edge & e, OutputIterator points) const {
    return dual(e.first, e.second, e.third, points);
  }

  template <class OutputIterator>
  OutputIterator dual(Cell_handle c, int i, int j,
      OutputIterator points) const {
    Base::dual(c, i, j, points, geom_traits().construct_weighted_circumcenter_3_object());
    return points;
  }

  template <class OutputIterator>
  OutputIterator dual(Vertex_handle v, OutputIterator points) const {
    Base::dual(v, points, geom_traits().construct_weighted_circumcenter_3_object());
    return points;
  }

  template <class Stream>
  Stream& draw_dual(Stream& os) const {
    return Base::draw_dual(os, geom_traits().construct_weighted_circumcenter_3_object());
  }

  /// Volume computations

  FT dual_volume(Vertex_handle v) const {
    return Base::dual_volume(v, geom_traits().construct_weighted_circumcenter_3_object());
  }

  /// Centroid computations

  Bare_point dual_centroid(Vertex_handle v) const {
    return Base::dual_centroid(v, geom_traits().construct_weighted_circumcenter_3_object());
  }
//@}

  template <class OutputIteratorBoundaryFacets, class OutputIteratorCells>
  std::pair<OutputIteratorBoundaryFacets, OutputIteratorCells>
  find_conflicts(const Weighted_point &p, Cell_handle c,
      OutputIteratorBoundaryFacets bfit, OutputIteratorCells cit) const {
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
  is_extensible_triangulation_in_1_sheet_h1() const {
    if (!is_1_cover())
      return can_be_converted_to_1_sheet();
    return is_extensible_triangulation_in_1_sheet_h2();
  }

  inline bool
  is_extensible_triangulation_in_1_sheet_h2() const
  {
    typedef typename Geometric_traits::Construct_weighted_circumcenter_3 Construct_weighted_circumcenter_3;
    typedef typename Geometric_traits::FT FT;

    FT threshold = FT(0.015625) * (domain().xmax()-domain().xmin()) * (domain().xmax()-domain().xmin());
    Construct_weighted_circumcenter_3 construct_weighted_circumcenter = geom_traits().construct_weighted_circumcenter_3_object();

    for (Periodic_tetrahedron_iterator tit = this->periodic_tetrahedra_begin(Base::UNIQUE);
         tit != this->periodic_tetrahedra_end(Base::UNIQUE);
         ++tit)
    {
      Bare_point cc = construct_weighted_circumcenter(tit->at(0).first, tit->at(1).first, tit->at(2).first, tit->at(3).first,
                                        tit->at(0).second, tit->at(1).second, tit->at(2).second, tit->at(3).second);

      if (squared_orthoball_radius(tit->at(0), tit->at(1), tit->at(2), tit->at(3)) >= threshold)
        return false;
    }
    return true;
  }
};

template < class Gt, class Tds >
template <class OutputIterator>
OutputIterator
Periodic_3_regular_triangulation_3<Gt,Tds>::vertices_in_conflict(
    const Weighted_point&p, Cell_handle c, OutputIterator res) const {
  if (number_of_vertices() == 0) return res;

  // Get the facets on the boundary of the hole.
  std::vector<Facet> facets;
  find_conflicts(p, c, std::back_inserter(facets), Emptyset_iterator());

  // Then extract uniquely the vertices.
  std::set<Vertex_handle> vertices;
  for (typename std::vector<Facet>::const_iterator i = facets.begin();
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
Periodic_3_regular_triangulation_3<Gt,Tds>::find_conflicts( const Weighted_point
&p,
    Cell_handle c, OutputIteratorBoundaryFacets bfit,
    OutputIteratorCells cit, OutputIteratorInternalFacets ifit) const {
  CGAL_triangulation_precondition(number_of_vertices() != 0);

  std::vector<Facet> facets;
  facets.reserve(64);
  std::vector<Cell_handle> cells;
  cells.reserve(32);

  Conflict_tester tester(p, this);
  Triple<typename std::back_insert_iterator<std::vector<Facet> >,
         typename std::back_insert_iterator<std::vector<Cell_handle> >,
         OutputIteratorInternalFacets> tit = Base::find_conflicts(c, tester,
              make_triple(std::back_inserter(facets),
                      std::back_inserter(cells), ifit));
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

  for (typename std::vector<Vertex_handle>::iterator
   voit = this->v_offsets.begin();
       voit != this->v_offsets.end() ; ++voit) {
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
  os= side_of_oriented_power_sphere(p0, p1, p2, p3, q, o0, o1, o2, o3, oq);

  if (os != ON_ORIENTED_BOUNDARY || !perturb)
    return (Bounded_side) os;

  //We are now in a degenerate case => we do a symbolic perturbation.
  // We sort the points lexicographically.
  Periodic_point pts[5] = {std::make_pair(p0,o0), std::make_pair(p1,o1),
         std::make_pair(p2,o2), std::make_pair(p3,o3),
         std::make_pair(q,oq)};
  const Periodic_point *points[5] ={&pts[0],&pts[1],&pts[2],&pts[3],&pts[4]};

  std::sort(points, points+5,
      typename Base::template Perturbation_order< typename Gt::Compare_xyz_3 >(geom_traits().compare_xyz_3_object()));

  // We successively look whether the leading monomial, then 2nd monomial
  // of the determinant has non null coefficient.
  for (int i=4; i>1; --i) {
    if (points[i] == &pts[4]) {
      CGAL_triangulation_assertion(orientation(p0, p1, p2, p3, o0, o1, o2, o3)
          == POSITIVE);
      // since p0 p1 p2 p3 are non coplanar and positively oriented
      return ON_UNBOUNDED_SIDE;
    }
    Orientation o;
    if (points[i] == &pts[3] &&
        (o = orientation(p0, p1, p2, q, o0, o1, o2, oq)) != COPLANAR ) {
      return (Bounded_side) o;
    }
    if (points[i] == &pts[2] &&
        (o = orientation(p0, p1, q, p3, o0, o1, oq, o3)) != COPLANAR ) {
      return (Bounded_side) o;
    }
    if (points[i] == &pts[1] &&
        (o = orientation(p0, q, p2, p3, o0, oq, o2, o3)) != COPLANAR ) {
      return (Bounded_side) o;
    }
    if (points[i] == &pts[0] &&
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
  if (!Base::is_valid(verbose, level)) {
    if (verbose)
      std::cerr << "Regular: invalid base" << std::endl;
    return false;
  }

  Conflict_tester tester(this);
  if (!is_valid_conflict(tester, verbose, level)) {
    if (verbose)
      std::cerr << "Regular: conflict problems" << std::endl;
    return false;
  }

  if (verbose)
    std::cerr << "Regular valid triangulation" << std::endl;
  return true;
}

template < class GT, class TDS >
bool
Periodic_3_regular_triangulation_3<GT,TDS>::
is_valid(Cell_handle ch, bool verbose, int level) const {
  bool error = false;
  if (!Base::is_valid(ch, verbose, level)) {
    error = true;
    if (verbose) {
      std::cerr << "geometrically invalid cell" << std::endl;
      for (int i=0; i<4; i++ )
  std::cerr << ch->vertex(i)->point() << ", ";
      std::cerr << std::endl;
    }
  }
  for (Vertex_iterator vit = vertices_begin(); vit != vertices_end(); ++ vit) {
    for (int i=-1; i<=1; i++)
      for (int j=-1; j<=1; j++)
  for (int k=-1; k<=1; k++) {
    if (periodic_point(ch,0) == std::make_pair(periodic_point(vit).first,
      periodic_point(vit).second+Offset(i,j,k))
    || periodic_point(ch,1) == std::make_pair(periodic_point(vit).first,
      periodic_point(vit).second+Offset(i,j,k))
          || periodic_point(ch,2) == std::make_pair(periodic_point(vit).first,
      periodic_point(vit).second+Offset(i,j,k))
    || periodic_point(ch,3) == std::make_pair(periodic_point(vit).first,
                  periodic_point(vit).second+Offset(i,j,k)) )
      continue;
    if (_side_of_power_sphere(ch, periodic_point(vit).first,
      periodic_point(vit).second+Offset(i,j,k),true)
        != ON_UNBOUNDED_SIDE) {
      error = true;
      if (verbose) {
        std::cerr << "Regular invalid cell" << std::endl;
        for (int i=0; i<4; i++ ) {
    Periodic_point pp = periodic_point(ch,i);
    std::cerr <<"("<<pp.first <<","<<pp.second<< "), ";
        }
        std::cerr << std::endl;
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
  const Self *t;
  Weighted_point p;
  // stores the offset of a point in 27-cover
  mutable Offset o;

public:
  /// Constructor
  Conflict_tester(const Self *_t) : t(_t), p(Weighted_point()) {}
  Conflict_tester(const Weighted_point &pt, const Self *_t) : t(_t), p(pt) { }

  /** The functor
    *
    * gives true if the circumcircle of c contains p
    */
  bool operator()(const Cell_handle c, const Offset &off) const {
    return (t->_side_of_power_sphere(c, p, t->combine_offsets(o, off), true)
        == ON_BOUNDED_SIDE);
  }

  bool operator()(const Cell_handle c, const Weighted_point& pt,
      const Offset &off) const {
    return (t->_side_of_power_sphere(c, pt, o + off, true) == ON_BOUNDED_SIDE);
  }

  int compare_weight(const Weighted_point& p, const Weighted_point& q) const
  {
    return t->power_test(p, q);
  }

  bool test_initial_cell(Cell_handle c, const Offset &off) const
  {
    return (operator()(c, off));
  }

  void set_point(const Weighted_point &_p) {
    p = _p;
  }

  void set_offset(const Offset &off) const {
    o = off;
  }

  const Offset &get_offset() const {
    return o;
  }

  const Weighted_point &point() const {
    return p;
  }

};

template < class GT, class Tds>
class Periodic_3_regular_triangulation_3<GT,Tds>::Point_hider
{
  Self *t;
  mutable std::vector<Vertex_handle> vertices;
  mutable std::vector<Weighted_point> hidden_points;
  mutable bool is_original_cube;

public:
  Point_hider(Self *tr) : t(tr), is_original_cube(false) {}

  void set_original_cube (bool b) const {
    is_original_cube = b;
  }

  template <class InputIterator>
  inline void set_vertices(InputIterator start, InputIterator end) const
  {
    while (start != end) {
        std::copy((*start)->hidden_points_begin(),
            (*start)->hidden_points_end(),
            std::back_inserter(hidden_points));

      for (int i=0; i<=3; i++) {
        Vertex_handle v = (*start)->vertex(i);
        if (v->cell() != Cell_handle()) {
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
    for (typename std::vector<Vertex_handle>::iterator
        vi = vertices.begin(); vi != vertices.end(); ++vi) {
      if ((*vi)->cell() != Cell_handle()) continue;
      if (is_original_cube)
      {
        hc = t->locate((*vi)->point(), lt, li, lj, hc);
        hc->hide_point((*vi)->point());
      }
      t->delete_vertex(*vi);
    }
    vertices.clear();
      for (typename std::vector<Weighted_point>::iterator
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
    if (is_original_cube)
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
//      const Conflict_tester &)
//  {
//    // No points to hide in the Delaunay triangulation.
//  }
};

#ifndef CGAL_CFG_OUTOFLINE_TEMPLATE_MEMBER_DEFINITION_BUG
template <class GT, class Tds>
template <class TriangulationR3>
struct Periodic_3_regular_triangulation_3<GT,Tds>::Vertex_remover
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
  private:
  // The removal of v may un-hide some points,
  // Space functions output them.
  std::vector<Weighted_point> hidden;
};
#endif //CGAL_CFG_OUTOFLINE_TEMPLATE_MEMBER_DEFINITION_BUG

template < class GT, class TDS >
std::istream &
operator>> (std::istream& is, Periodic_3_regular_triangulation_3<GT,TDS> &tr)
{
  typedef Periodic_3_regular_triangulation_3<GT,TDS>   P3RT3;
  typedef typename P3RT3::Base                          Base;
  typedef typename P3RT3::Vertex_iterator               Vertex_iterator;
  typedef typename GT::FT FT;
  typedef typename P3RT3::Vertex_handle                 Vertex_handle;

  is >> static_cast<Base&>(tr);

  tr.insert_too_long_edges(tr.cells_begin(), tr.cells_end());

  CGAL_triangulation_expensive_assertion( tr.is_valid() );
  return is;
}
}// namespace CGAL

#endif
