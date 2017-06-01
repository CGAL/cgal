// Copyright (c) 1999-2004,2006-2009, 2017  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//                 Sylvain Pion <Sylvain.Pion@sophia.inria.fr>
//                 Andreas Fabri <Andreas.Fabri@sophia.inria.fr>
//                 Nico Kruithof <Nico.Kruithof@sophia.inria.fr>
//                 Manuel Caroli <Manuel.Caroli@sophia.inria.fr>

#ifndef CGAL_PERIODIC_3_DELAUNAY_TRIANGULATION_3_H
#define CGAL_PERIODIC_3_DELAUNAY_TRIANGULATION_3_H

#include <CGAL/license/Periodic_3_triangulation_3.h>

#include <CGAL/Periodic_3_triangulation_3.h>
#include <CGAL/spatial_sort.h>

// Needed by remove to fill the hole.
#include <CGAL/internal/Periodic_3_Delaunay_triangulation_remove_traits_3.h>
#include <CGAL/Delaunay_triangulation_3.h>

#include <algorithm>
#include <iostream>
#include <vector>
#include <utility>

namespace CGAL {

template < class GT, class TDS > class Periodic_3_Delaunay_triangulation_3;
template < class GT, class TDS > std::istream& operator>>
    (std::istream& is, Periodic_3_Delaunay_triangulation_3<GT,TDS> &tr);

template < class Gt,
           class Tds = Triangulation_data_structure_3 <
                         Triangulation_vertex_base_3<
                           Gt, Periodic_3_triangulation_ds_vertex_base_3<>
                         >,
                         Triangulation_cell_base_3<
                           Gt, Periodic_3_triangulation_ds_cell_base_3<>
                         >
                       >
         >
class Periodic_3_Delaunay_triangulation_3 :
  public Periodic_3_triangulation_3<Gt, Tds>
{
  friend std::istream& operator>> <>
  (std::istream& is, Periodic_3_Delaunay_triangulation_3<Gt, Tds> &tr);

  typedef Periodic_3_Delaunay_triangulation_3<Gt,Tds>          Self;
public:
  typedef Periodic_3_triangulation_3<Gt,Tds>                   Base;

public:
  /** @name Template parameter types */ //@{
  typedef Gt                                    Geometric_traits;
  typedef Tds                                   Triangulation_data_structure;
  //@}

  ///Compatibility typedef:
  typedef Geometric_traits                      Geom_traits;
  typedef typename Gt::FT                       FT;

  typedef typename Gt::Point_3                  Point;
  typedef typename Gt::Segment_3                Segment;
  typedef typename Gt::Triangle_3               Triangle;
  typedef typename Gt::Tetrahedron_3            Tetrahedron;

  typedef typename Base::Periodic_point         Periodic_point;
  typedef typename Base::Periodic_segment       Periodic_segment;
  typedef typename Base::Periodic_triangle      Periodic_triangle;
  typedef typename Base::Periodic_tetrahedron   Periodic_tetrahedron;

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
  typedef typename Base::Periodic_segment_iterator  Periodic_segment_iterator;
  typedef typename Base::Periodic_tetrahedron_iterator  Periodic_tetrahedron_iterator;
  //@}

#ifndef CGAL_CFG_USING_BASE_MEMBER_BUG_2
  using Base::cw;
  using Base::ccw;
  using Base::domain;
  using Base::geom_traits;
  using Base::insert_dummy_points;
  using Base::int_to_off;
  using Base::is_1_cover;
  using Base::is_virtual;
  using Base::number_of_sheets;
  using Base::number_of_vertices;
  using Base::number_of_edges;
  using Base::number_of_facets;
  using Base::number_of_cells;
  using Base::next_around_edge;
  using Base::vertex_triple_index;
  using Base::mirror_vertex;
  using Base::orientation;
  using Base::point;
  using Base::swap;
  using Base::tds;
  using Base::vertices_begin;
  using Base::vertices_end;
  using Base::edges_begin;
  using Base::edges_end;
  using Base::facets_begin;
  using Base::facets_end;
  using Base::cells_begin;
  using Base::cells_end;
#endif

  // For strict-ansi compliance
  using Base::adjacent_vertices;
  using Base::combine_offsets;
  using Base::construct_point;
  using Base::convert_to_27_sheeted_covering;
  using Base::draw_dual;
  using Base::incident_edges;
  using Base::incident_facets;
  using Base::incident_cells;
  using Base::is_valid_conflict;
  using Base::get_offset;
  using Base::get_original_vertex;
  using Base::get_location_offset;
  using Base::locate;
  using Base::neighbor_offset;
  using Base::periodic_point;

private:
  /// This threshold should be chosen such that if all edges are shorter,
  /// we can be sure that there are no self-edges anymore.
  FT edge_length_threshold;

  /// This adjacency list stores all edges that are longer than
  /// edge_length_threshold.
  std::map< Vertex_handle, std::list<Vertex_handle> > too_long_edges;
  unsigned int too_long_edge_counter;

  class Cover_manager
  {
    Periodic_3_Delaunay_triangulation_3& tr;

  public:
    Cover_manager (Periodic_3_Delaunay_triangulation_3& tr)
    : tr(tr)
    {}

    void create_initial_triangulation()
    {
      tr.create_initial_triangulation();
    }

    template <class CellIt>
    void delete_unsatisfying_elements(const CellIt begin, const CellIt end)
    {
      tr.delete_too_long_edges(begin, end);
    }

    template <class CellIt>
    void insert_unsatisfying_elements(Vertex_handle v, const CellIt begin, const CellIt end)
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
  Periodic_3_Delaunay_triangulation_3(const Iso_cuboid& domain = Iso_cuboid(0,0,0,1,1,1),
                                      const Geometric_traits& gt = Geometric_traits())
    : Base(domain, gt), too_long_edge_counter(0)
  {
    edge_length_threshold = FT(0.166) * (domain.xmax()-domain.xmin())
                                      * (domain.xmax()-domain.xmin());
  }

  // copy constructor duplicates vertices and cells
  Periodic_3_Delaunay_triangulation_3(const Periodic_3_Delaunay_triangulation_3& tr)
    : Base(tr), edge_length_threshold(tr.edge_length_threshold)
  {
    if(is_1_cover()) {
      tds() = tr.tds();
    } else {
      this->copy_multiple_covering(tr);
    }
    CGAL_triangulation_expensive_postcondition(*this == tr);
    CGAL_triangulation_expensive_postcondition( is_valid() );
  }

  template < typename InputIterator >
  Periodic_3_Delaunay_triangulation_3(InputIterator first, InputIterator last,
                                      const Iso_cuboid& domain = Iso_cuboid(0,0,0,1,1,1),
                                      const Geometric_traits& gt = Geometric_traits())
    : Base(domain, gt), too_long_edge_counter(0)
  {
    edge_length_threshold = FT(0.166) * (domain.xmax()-domain.xmin())
                                      * (domain.xmax()-domain.xmin());
    insert(first, last);
  }

  Periodic_3_Delaunay_triangulation_3 operator=(Periodic_3_Delaunay_triangulation_3 tr)
  {
    tr.swap(*this);
    return *this;
  }

  void copy_multiple_covering(const Periodic_3_Delaunay_triangulation_3 & tr);

  void swap(Periodic_3_Delaunay_triangulation_3&tr)
  {
    std::swap(edge_length_threshold,tr.edge_length_threshold);
    std::swap(too_long_edges,tr.too_long_edges);
    std::swap(too_long_edge_counter,tr.too_long_edge_counter);
    Base::swap(tr);
  }

  virtual void clear_covering_data()
  {
    too_long_edges.clear();
    too_long_edge_counter = 0;
  }

  virtual void update_cover_data_after_setting_domain ()
  {
    edge_length_threshold = FT(0.166) * (domain().xmax()-domain().xmin())
                                      * (domain().xmax()-domain().xmin());
  }

  virtual void update_cover_data_after_converting_to_27_sheeted_covering()
  {
    // Set up too long edges data structure
    int i=0;
    for(Vertex_iterator vit = vertices_begin(); vit != vertices_end(); ++vit) {
      too_long_edges[vit] = std::list<Vertex_handle>();
      ++i;
    }
    too_long_edge_counter = this->find_too_long_edges(too_long_edges);
  }

  bool is_extensible_triangulation_in_1_sheet_h1() const;
  bool is_extensible_triangulation_in_1_sheet_h2() const;

  // iterate over all edges and store the ones that are longer than
  // edge_length_threshold in edges. Return the number of too long edges.
  int find_too_long_edges(std::map<Vertex_handle, std::list<Vertex_handle> >& edges) const
  {
    Point p1, p2;
    int counter = 0;
    Vertex_handle v_no,vh;
    for(Edge_iterator eit = edges_begin(); eit != edges_end(); eit++) {
      p1 = construct_point(eit->first->vertex(eit->second)->point(),
                           get_offset(eit->first, eit->second));
      p2 = construct_point(eit->first->vertex(eit->third)->point(),
                           get_offset(eit->first, eit->third));
      if(squared_distance(p1,p2) > edge_length_threshold) {
        if(&*(eit->first->vertex(eit->second)) < &*(eit->first->vertex(eit->third))) {
          v_no = eit->first->vertex(eit->second);
          vh = eit->first->vertex(eit->third);
        } else {
          v_no = eit->first->vertex(eit->third);
          vh = eit->first->vertex(eit->second);
        }
        edges[v_no].push_back(vh);
        counter++;
      }
    }
    return counter;
  }
  //@}

  void create_initial_triangulation()
  {
    // create the base for too_long_edges;
    CGAL_triangulation_assertion( too_long_edges.empty() );
    CGAL_triangulation_assertion(too_long_edge_counter == 0);

    for(Vertex_iterator vit = vertices_begin(); vit !=vertices_end(); ++vit )
      too_long_edges[vit] = std::list<Vertex_handle>();;

    std::vector<Cell_handle> temp_inc_cells;
    for(Vertex_iterator vit = vertices_begin(); vit !=vertices_end(); ++vit ) {
      temp_inc_cells.clear();
      incident_cells(vit, std::back_inserter(temp_inc_cells));
      for(unsigned int i=0; i<temp_inc_cells.size(); i++) {
        int k = temp_inc_cells[i]->index(vit);
        for(int j=0; j<4; j++) {
          if(j==k) continue;
          if(&*vit > &*(temp_inc_cells[i]->vertex(j))) continue;
          if((find(too_long_edges[vit].begin(),
                   too_long_edges[vit].end(),
                   temp_inc_cells[i]->vertex(j)) ==
              too_long_edges[vit].end())) {
            too_long_edges[vit].push_back(temp_inc_cells[i]->vertex(j));
            too_long_edge_counter++;
          }
        }
      }
    }
  }

  template <class CellIt>
  void delete_too_long_edges(const CellIt begin, const CellIt end)
  {
    std::pair< Vertex_handle, Vertex_handle > edge_to_delete, edge_to_delete2;
    typename std::list< Vertex_handle >::iterator sit;
    // Iterate over all cells that are in the star. That means that those cells
    // are going to be deleted. Therefore, all of them have to be deleted from
    // too_long_edges, if they are contained in it.
    for(CellIt it = begin; it != end; ++it) {
      for(int j=0; j<4; j++) {
        for(int k=0; k<4; k++) {
          if(&*((*it)->vertex(j)) < &*((*it)->vertex(k))) {
            edge_to_delete = std::make_pair((*it)->vertex(j),(*it)->vertex(k));
          } else {
            edge_to_delete = std::make_pair((*it)->vertex(k),(*it)->vertex(j));
          }
          Vertex_handle v_no = edge_to_delete.first;
          sit = std::find(too_long_edges[v_no].begin(),
              too_long_edges[v_no].end(),
              edge_to_delete.second);
          if(sit != too_long_edges[v_no].end()) {
            too_long_edges[v_no].erase(sit);
            too_long_edge_counter--;
          }
        }
      }
    }
  }

  template <class CellIt>
  void insert_too_long_edges(Vertex_handle v, const CellIt begin, const CellIt end)
  {
    CGAL_triangulation_precondition(number_of_vertices() != 0);
    // add newly added edges to too_long_edges, if necessary.
    Point p1,p2;
    std::pair< Vertex_handle, Vertex_handle > edge_to_add;
    std::list<Vertex_handle> empty_list;
    too_long_edges[v] = empty_list;
    // Iterate over all cells of the new star.
    for(CellIt it = begin; it != end; ++it) {
      // Consider all possible vertex pairs.
      for(int k=0; k<4; k++) {
        for(int j=0; j<4; j++) {
          if(j==k) continue;
          if(&*((*it)->vertex(j)) > &*((*it)->vertex(k))) continue;
          // make the offsets canonical (wrt. to some notion)
          // add to too_long_edges, if not yet added and if "too long"
          CGAL_triangulation_precondition(
                &*((*it)->vertex(j))< &*((*it)->vertex(k)));

          edge_to_add = std::make_pair((*it)->vertex(j), (*it)->vertex(k));

          p1 = construct_point((*it)->vertex(j)->point(), get_offset(*it, j));
          p2 = construct_point((*it)->vertex(k)->point(), get_offset(*it, k));

          if((squared_distance(p1,p2) > edge_length_threshold)
             && (find(too_long_edges[(*it)->vertex(j)].begin(),
                      too_long_edges[(*it)->vertex(j)].end(),
                      edge_to_add.second)
                 == too_long_edges[(*it)->vertex(j)].end())) {
            too_long_edges[(*it)->vertex(j)].push_back(edge_to_add.second);
            too_long_edge_counter++;
          }
        }
      }
    }
  }

  bool can_be_converted_to_1_sheet() const
  {
    return too_long_edge_counter == 0;
  }

  bool update_cover_data_during_management(Cell_handle new_ch,
                                           const std::vector<Cell_handle>& new_cells)
  {
    for(int i=0; i < 4; i++) {
      for(int j=0; j < 4; j++) {
        if(j==i) continue;
        if(&*(new_ch->vertex(i)) > &*(new_ch->vertex(j))) continue;

        Point p1 = construct_point(new_ch->vertex(i)->point(),
                                   get_offset(new_ch, i));
        Point p2 = construct_point(new_ch->vertex(j)->point(),
                                   get_offset(new_ch, j));
        Vertex_handle v_no = new_ch->vertex(i);

        if(squared_distance(p1, p2) > edge_length_threshold) {
          // If the cell does not fulfill the edge-length criterion
          // revert all changes to the triangulation and transform it
          // to a triangulation in the needed covering space.
          if(is_1_cover()) {
            tds().delete_cells(new_cells.begin(), new_cells.end());
            convert_to_27_sheeted_covering();
            return true;
          }
          else if(find(too_long_edges[v_no].begin(),
                       too_long_edges[v_no].end(),
                       new_ch->vertex(j))
                  == too_long_edges[v_no].end()) {
            too_long_edges[v_no].push_back(new_ch->vertex(j));
            too_long_edge_counter++;
          }
        }
      }
    }
    return false;
  }

  /** @name Insertion */ //@{
  Vertex_handle insert(const Point& p, Cell_handle start = Cell_handle())
  {
    Conflict_tester tester(p, this);
    Point_hider hider;
    Cover_manager cover_manager(*this);
    return Base::insert_in_conflict(p, start, tester, hider, cover_manager);
  }

  Vertex_handle insert(const Point& p, Locate_type lt, Cell_handle c,
                       int li, int lj)
  {
    Conflict_tester tester(p, this);
    Point_hider hider;
    Cover_manager cover_manager(*this);
    return Base::insert_in_conflict(p,lt,c,li,lj, tester,hider, cover_manager);
  }

  template < class InputIterator >
  std::ptrdiff_t insert(InputIterator first, InputIterator last,
                        bool is_large_point_set = false)
  {
    if(first == last) return 0;
    size_type n = number_of_vertices();
    // The heuristic discards the existing triangulation so it can only be
    // applied to empty triangulations.
    if(n!=0)
      is_large_point_set = false;

    std::vector<Point> points(first, last);
    std::random_shuffle (points.begin(), points.end());
    Cell_handle hint;
    std::vector<Vertex_handle> dummy_points, double_vertices;
    typename std::vector<Point>::iterator pbegin = points.begin();
    if(is_large_point_set)
      dummy_points = insert_dummy_points();
    else while(!is_1_cover()) {
      insert(*pbegin);
      ++pbegin;
      if(pbegin == points.end())
        return number_of_vertices() - n;
    }

    spatial_sort (pbegin, points.end(), this->geom_traits());

    Conflict_tester tester(*pbegin,this);
    Point_hider hider;
    Cover_manager cover_manager(*this);
    double_vertices = Base::insert_in_conflict(
      points.begin(), points.end(), hint, tester, hider, cover_manager);

    if(is_large_point_set) {
      typedef CGAL::Periodic_3_Delaunay_triangulation_remove_traits_3<Gt> P3removeT;
      typedef CGAL::Delaunay_triangulation_3< P3removeT > DT;
      typedef Vertex_remover< DT > Remover;
      P3removeT remove_traits(domain());
      DT dt(remove_traits);
      Remover remover(this,dt);
      Conflict_tester t(this);
      for(unsigned int i=0; i<dummy_points.size(); i++) {
        if(std::find(double_vertices.begin(), double_vertices.end(),
                      dummy_points[i]) == double_vertices.end())
          Base::remove(dummy_points[i],remover,t, cover_manager);
      }
    }

    return number_of_vertices() - n;
  }
  //@}

  /** @name Point moving */ //@{
  // @todo should be deprecated and a function move() should be introduced
  // see what is done in /Triangulation_3
  // Also need to introduce move() for periodic regular triangulations
  Vertex_handle move_point(Vertex_handle v, const Point& p);
  //@}

public:
  /** @name Removal */ //@{
  void remove(Vertex_handle v);

  template < typename InputIterator >
  std::ptrdiff_t remove(InputIterator first, InputIterator beyond)
  {
    std::size_t n = number_of_vertices();
    while(first != beyond) {
      remove(*first);
      ++first;
    }
    return n - number_of_vertices();
  }
  //@}

public:
  /** @name Wrapping the traits */ //@{
  Oriented_side side_of_oriented_sphere(const Point& p, const Point& q,
                                        const Point& r, const Point& s,
                                        const Point& t) const {
    return geom_traits().side_of_oriented_sphere_3_object()(p,q,r,s,t);
  }
  Oriented_side side_of_oriented_sphere(const Point& p, const Point& q,
                                        const Point& r, const Point& s,
                                        const Point& t,
                                        const Offset& o_p, const Offset& o_q,
                                        const Offset& o_r, const Offset& o_s,
                                        const Offset& o_t) const {
    return geom_traits().side_of_oriented_sphere_3_object()(
          p,q,r,s,t,o_p,o_q,o_r,o_s,o_t);
  }
  Comparison_result compare_distance(const Point& p, const Point& q,
                                     const Point& r) const {
      return geom_traits().compare_distance_3_object()(p, q, r);
  }
  Comparison_result compare_distance(const Point& p, const Point& q,
                                     const Point& r,
                                     const Offset& o_p, const Offset& o_q,
                                     const Offset& o_r) const {
    return geom_traits().compare_distance_3_object()(p, q, r, o_p, o_q, o_r);
  }
  //@}

private:
  /** @name Query helpers */ //@{
  Bounded_side _side_of_sphere(const Cell_handle& c, const Point& p,
      const Offset & offset = Offset(), bool perturb = false) const;

  Offset get_min_dist_offset(const Point& p, const Offset& o,
                             const Vertex_handle vh) const;
  //@}

public:
  /** @name Queries */ //@{
  Bounded_side side_of_sphere(const Cell_handle& c, const Point& p,
      const Offset & offset = Offset(), bool perturb = false) const{
    Bounded_side bs = ON_UNBOUNDED_SIDE;
    int i=0;
    // TODO: optimize which copies to check depending on the offsets in
    // the cell.
    while(bs == ON_UNBOUNDED_SIDE && i<8) {
      bs= _side_of_sphere(c,p,combine_offsets(offset,int_to_off(i)),perturb);
      i++;
    }
    return bs;
  }

  Vertex_handle nearest_vertex(const Point& p, Cell_handle c = Cell_handle()) const;
  Vertex_handle nearest_vertex_in_cell(const Cell_handle& c, const Point & p,
                                       const Offset & offset = Offset()) const;

  /// Undocumented wrapper for find_conflicts.
  template <class OutputIteratorBoundaryFacets, class OutputIteratorCells>
  std::pair<OutputIteratorBoundaryFacets, OutputIteratorCells>
  find_conflicts(const Point &p, Cell_handle c,
                 OutputIteratorBoundaryFacets bfit,
                 OutputIteratorCells cit) const
  {
    Triple<OutputIteratorBoundaryFacets,OutputIteratorCells,Emptyset_iterator>
    t = find_conflicts(p, c, bfit, cit, Emptyset_iterator());
    return std::make_pair(t.first, t.second);
  }

  template <class OutputIteratorBoundaryFacets, class OutputIteratorCells,
            class OutputIteratorInternalFacets>
  Triple<OutputIteratorBoundaryFacets, OutputIteratorCells,
         OutputIteratorInternalFacets>
  find_conflicts(const Point &p, Cell_handle c,
                 OutputIteratorBoundaryFacets bfit, OutputIteratorCells cit,
                 OutputIteratorInternalFacets ifit) const;

  /// Returns the vertices on the boundary of the conflict hole.
  template <class OutputIterator>
  OutputIterator vertices_in_conflict(const Point&p, Cell_handle c,
                                      OutputIterator res) const;

  bool is_Gabriel(const Cell_handle c, int i) const;
  bool is_Gabriel(const Cell_handle c, int i, int j) const;
  bool is_Gabriel(const Facet& f)const {
    return is_Gabriel(f.first, f.second);
  }
  bool is_Gabriel(const Edge& e) const {
    return is_Gabriel(e.first, e.second, e.third);
  }
  //@}

private:
  /** @name Voronoi diagram helpers */ //@{
  bool is_canonical(const Periodic_segment &ps) const
  {
    if(number_of_sheets() == make_array(1,1,1)) return true;
    Offset o0 = ps.at(0).second;
    Offset o1 = ps.at(1).second;
    Offset cumm_off((std::min)(o0.x(),o1.x()),(std::min)(o0.y(),o1.y()),
                    (std::min)(o0.z(),o1.z()));
    return (cumm_off == Offset(0,0,0));
  }
  //@}

public:
  /** @name Geometric access functions */
  /// @{

  Point point(const Periodic_point& pp) const
  {
    return point(pp, geom_traits().construct_point_3_object());
  }

  // The following functions return the "real" position in space (unrestrained
  // to the original periodic domain) of the vertices v and c->vertex(idx),
  // respectively

  Point point(Vertex_handle v) const
  {
    return point(v, geom_traits().construct_point_3_object());
  }

  Point point(Cell_handle c, int idx) const
  {
    return point(c, idx, geom_traits().construct_point_3_object());
  }

  // end of geometric functions
  /// @}

  Periodic_point periodic_circumcenter(Cell_handle c) const {
    return Base::periodic_circumcenter(c, geom_traits().construct_circumcenter_3_object());
  }

public:
  /** @name Voronoi diagram */ //@{
  Point dual(Cell_handle c) const {
    return point(periodic_circumcenter(c));
  }

  bool canonical_dual_segment(Cell_handle c, int i, Periodic_segment& ps) const {
    return Base::canonical_dual_segment(c, i, ps, geom_traits().construct_circumcenter_3_object());
  }
  Periodic_segment dual(const Facet & f) const {
    return dual( f.first, f.second );
  }
  Periodic_segment dual(Cell_handle c, int i) const
  {
    Periodic_segment ps;
    canonical_dual_segment(c,i,ps);
    return ps;
  }

  template <class OutputIterator>
  OutputIterator dual(const Edge & e, OutputIterator points) const {
    return dual(e.first, e.second, e.third, points);
  }

  template <class OutputIterator>
  OutputIterator dual(Cell_handle c, int i, int j, OutputIterator points) const {
    Base::dual(c, i, j, points, geom_traits().construct_circumcenter_3_object());
    return points;
  }

  template <class OutputIterator>
  OutputIterator dual(Vertex_handle v, OutputIterator points) const {
    Base::dual(v, points, geom_traits().construct_circumcenter_3_object());
    return points;
  }

  template <class Stream>
  Stream& draw_dual(Stream& os) const {
    return Base::draw_dual(os, geom_traits().construct_circumcenter_3_object());
  }

  /// Volume computations
  FT dual_volume(Vertex_handle v) const {
    return Base::dual_volume(v, geom_traits().construct_circumcenter_3_object());
  }

  /// Centroid computations
  Point dual_centroid(Vertex_handle v) const {
    return Base::dual_centroid(v, geom_traits().construct_circumcenter_3_object());
  }
  //@}

  /** @name Checking */ //@{
  bool is_valid(bool verbose = false, int level = 0) const;
  bool is_valid(Cell_handle c, bool verbose = false, int level = 0) const;
  //@}

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

    typedef typename std::vector<Point>::iterator Hidden_points_iterator;

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

    void add_hidden_points(Cell_handle) {
      std::copy (hidden_points_begin(), hidden_points_end(),
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
    std::vector<Point> hidden;
  };
#endif //CGAL_CFG_OUTOFLINE_TEMPLATE_MEMBER_DEFINITION_BUG

  // unused and undocumented types and functions required to be
  // compatible to Alpha_shape_3
public:
  typedef Cell_iterator   Finite_cells_iterator;
  typedef Facet_iterator  Finite_facets_iterator;
  typedef Edge_iterator   Finite_edges_iterator;
  typedef Vertex_iterator Finite_vertices_iterator;

  int dimension() const { return 3; }
  template < class T >
  bool is_infinite(const T&, int = 0, int = 0) const { return false; }
  size_type number_of_finite_cells() const { return number_of_cells(); }
  size_type number_of_finite_facets() const { return number_of_facets(); }
  size_type number_of_finite_edges() const { return number_of_edges(); }
  size_type number_of_finite_vertices() const { return number_of_vertices(); }
};

template < class GT, class Tds >
typename Periodic_3_Delaunay_triangulation_3<GT,Tds>::Vertex_handle
Periodic_3_Delaunay_triangulation_3<GT,Tds>::nearest_vertex(const Point& p,
                                                            Cell_handle start) const
{
  if(number_of_vertices() == 0)
    return Vertex_handle();

  Locate_type lt;
  int li, lj;
  Cell_handle c = locate(p, lt, li, lj, start);
  if(lt == Base::VERTEX) return c->vertex(li);
  const Conflict_tester tester(p, this);
  Offset o = combine_offsets(Offset(),get_location_offset(tester,c));

  // - start with the closest vertex from the located cell.
  // - repeatedly take the nearest of its incident vertices if any
  // - if not, we're done.
  Vertex_handle nearest = nearest_vertex_in_cell(c, p, o);
  std::vector<Vertex_handle> vs;
  vs.reserve(32);
  while(true) {
    Vertex_handle tmp = nearest;
    Offset tmp_off = get_min_dist_offset(p,o,tmp);
    adjacent_vertices(nearest, std::back_inserter(vs));
    for(typename std::vector<Vertex_handle>::const_iterator
         vsit = vs.begin(); vsit != vs.end(); ++vsit)
      tmp = (compare_distance(p,tmp->point(),(*vsit)->point(),
                              o,tmp_off,get_min_dist_offset(p,o,*vsit))
             == SMALLER) ? tmp : *vsit;
    if(tmp == nearest)
      break;
    vs.clear();
    nearest = tmp;
  }

  return get_original_vertex(nearest);
}

// just trying the eight possibilities
template < class GT, class Tds >
typename Periodic_3_Delaunay_triangulation_3<GT,Tds>::Offset
Periodic_3_Delaunay_triangulation_3<GT,Tds>::get_min_dist_offset(
    const Point& p, const Offset& o, const Vertex_handle vh) const
{
  Offset mdo = get_offset(vh);
  Offset min_off = Offset(0,0,0);
  min_off = (compare_distance(p, vh->point(), vh->point(),
                              o, combine_offsets(mdo,min_off),
                                 combine_offsets(mdo,Offset(0,0,1)))
             == SMALLER ? min_off : Offset(0,0,1) );
  min_off = (compare_distance(p, vh->point(), vh->point(),
                              o, combine_offsets(mdo,min_off),
                                 combine_offsets(mdo,Offset(0,1,0)))
             == SMALLER ? min_off : Offset(0,1,0) );
  min_off = (compare_distance(p, vh->point(), vh->point(),
                              o, combine_offsets(mdo,min_off),
                                 combine_offsets(mdo,Offset(0,1,1)))
             == SMALLER ? min_off : Offset(0,1,1) );
  min_off = (compare_distance(p, vh->point(), vh->point(),
                              o, combine_offsets(mdo,min_off),
                                 combine_offsets(mdo,Offset(1,0,0)))
             == SMALLER ? min_off : Offset(1,0,0) );
  min_off = (compare_distance(p, vh->point(), vh->point(),
                              o, combine_offsets(mdo,min_off),
                                 combine_offsets(mdo,Offset(1,0,1)))
             == SMALLER ? min_off : Offset(1,0,1) );
  min_off = (compare_distance(p, vh->point(), vh->point(),
                              o, combine_offsets(mdo,min_off),
                                 combine_offsets(mdo,Offset(1,1,0)))
             == SMALLER ? min_off : Offset(1,1,0) );
  min_off = (compare_distance(p, vh->point(), vh->point(),
                              o, combine_offsets(mdo,min_off),
                                 combine_offsets(mdo,Offset(1,1,1)))
             == SMALLER ? min_off : Offset(1,1,1) );
  return combine_offsets(mdo,min_off);
}

/// Returns the finite vertex of the cell c which is the closest to p.
template < class GT, class Tds >
typename Periodic_3_Delaunay_triangulation_3<GT,Tds>::Vertex_handle
Periodic_3_Delaunay_triangulation_3<GT,Tds>::nearest_vertex_in_cell(
    const Cell_handle& c, const Point& p, const Offset& o) const
{
  CGAL_triangulation_precondition(number_of_vertices() != 0);
  Vertex_handle nearest = c->vertex(0);
  for(int i=1; i<4; i++) {
    nearest = (compare_distance(p,nearest->point(),c->vertex(i)->point(),
                                o, get_offset(c,c->index(nearest)),
                                   get_offset(c,i)) == SMALLER) ?
                nearest : c->vertex(i);
  }
  return nearest;
}

// ############################################################################

// TODO: reintroduce the commented lines.
template < class Gt, class Tds >
typename Periodic_3_Delaunay_triangulation_3<Gt,Tds>::Vertex_handle
Periodic_3_Delaunay_triangulation_3<Gt,Tds>::
move_point(Vertex_handle v, const Point& p)
{
  CGAL_triangulation_expensive_precondition(is_vertex(v));
  // Remember an incident vertex to restart
  // the point location after the removal.
  // Cell_handle c = v->cell();
  //Vertex_handle old_neighbor = c->vertex(c->index(v) == 0 ? 1 : 0);
  //  CGAL_triangulation_assertion(old_neighbor != v);

  remove(v);

  if(number_of_vertices() == 0)
    return insert(p);
  return insert(p);//, old_neighbor->cell());
}

template < class Gt, class Tds >
void Periodic_3_Delaunay_triangulation_3<Gt,Tds>::remove(Vertex_handle v)
{
  typedef CGAL::Periodic_3_Delaunay_triangulation_remove_traits_3<Gt> P3removeT;
  typedef CGAL::Delaunay_triangulation_3< P3removeT >
    Euclidean_triangulation;
  typedef Vertex_remover< Euclidean_triangulation > Remover;
  P3removeT remove_traits(domain());
  Euclidean_triangulation tmp(remove_traits);
  Remover remover(this, tmp);
  Conflict_tester ct(this);
  Cover_manager cover_manager(*this);

  Base::remove(v, remover, ct, cover_manager);
  CGAL_triangulation_expensive_assertion(is_valid());
}

template < class Gt, class Tds >
template <class OutputIteratorBoundaryFacets, class OutputIteratorCells,
          class OutputIteratorInternalFacets>
Triple<OutputIteratorBoundaryFacets, OutputIteratorCells,
       OutputIteratorInternalFacets>
Periodic_3_Delaunay_triangulation_3<Gt,Tds>::find_conflicts(
    const Point& p, Cell_handle c, OutputIteratorBoundaryFacets bfit,
    OutputIteratorCells cit, OutputIteratorInternalFacets ifit) const
{
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

  for(typename std::vector<Vertex_handle>::iterator
      voit = this->v_offsets.begin();
      voit != this->v_offsets.end(); ++voit) {
    (*voit)->clear_offset();
  }
  this->v_offsets.clear();

  return make_triple(bfit, cit, ifit);
}

template < class Gt, class Tds >
template <class OutputIterator>
OutputIterator
Periodic_3_Delaunay_triangulation_3<Gt,Tds>::vertices_in_conflict(
    const Point& p, Cell_handle c, OutputIterator res) const
{
  if(number_of_vertices() == 0) return res;

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
Bounded_side Periodic_3_Delaunay_triangulation_3<Gt,Tds>::
_side_of_sphere(const Cell_handle& c, const Point& q,
                const Offset& offset, bool perturb ) const
{
  Point p0 = c->vertex(0)->point(),
        p1 = c->vertex(1)->point(),
        p2 = c->vertex(2)->point(),
        p3 = c->vertex(3)->point();
  Offset o0 = this->get_offset(c,0),
        o1 = this->get_offset(c,1),
        o2 = this->get_offset(c,2),
        o3 = this->get_offset(c,3),
        oq = offset;

  Oriented_side os = ON_NEGATIVE_SIDE;
  os= side_of_oriented_sphere(p0, p1, p2, p3, q, o0, o1, o2, o3, oq);

  if(os != ON_ORIENTED_BOUNDARY || !perturb)
    return (Bounded_side) os;

  //We are now in a degenerate case => we do a symbolic perturbation.
  // We sort the points lexicographically.
  Periodic_point pts[5] = {std::make_pair(p0,o0), std::make_pair(p1,o1),
                           std::make_pair(p2,o2), std::make_pair(p3,o3),
                           std::make_pair(q,oq)};
  const Periodic_point *points[5] ={&pts[0],&pts[1],&pts[2],&pts[3],&pts[4]};

  std::sort(points, points+5, typename Base::Perturbation_order(this));

  // We successively look whether the leading monomial, then 2nd monomial
  // of the determinant has non null coefficient.
  // 2 iterations are enough (cf paper)
  for(int i=4; i>2; --i) {
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
bool Periodic_3_Delaunay_triangulation_3<Gt,Tds>::
is_Gabriel(const Cell_handle c, int i) const
{
  CGAL_triangulation_precondition(number_of_vertices() != 0);
  typename Geom_traits::Side_of_bounded_sphere_3
    side_of_bounded_sphere =
    geom_traits().side_of_bounded_sphere_3_object();

  if(side_of_bounded_sphere (
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

  if(side_of_bounded_sphere(
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

template < class Gt, class Tds >
bool Periodic_3_Delaunay_triangulation_3<Gt,Tds>::
is_Gabriel(const Cell_handle c, int i, int j) const
{
  typename Geom_traits::Side_of_bounded_sphere_3
    side_of_bounded_sphere =
    geom_traits().side_of_bounded_sphere_3_object();

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
    Offset off1 = int_to_off(cc->offset(i1));
    Offset off2 = int_to_off(cc->offset(i2));
    Offset off3 = int_to_off(cc->offset(i3));
    if(side_of_bounded_sphere(
          v1->point(), v2->point(), cc->vertex(fcirc->second)->point(),
          off1, off2, off3) == ON_BOUNDED_SIDE ) return false;
  } while(++fcirc != fdone);
  return true;
}

template < class Gt, class Tds >
bool
Periodic_3_Delaunay_triangulation_3<Gt,Tds>::
is_valid(bool verbose, int level) const
{
  if(!Base::is_valid(verbose, level)) {
    if(verbose)
      std::cerr << "Delaunay: invalid base" << std::endl;
    return false;
  }

  Conflict_tester tester(this);
  if(!is_valid_conflict(tester, verbose, level)) {
    if(verbose)
      std::cerr << "Delaunay: conflict problems" << std::endl;
    return false;
  }

  if(verbose)
    std::cerr << "Delaunay valid triangulation" << std::endl;
  return true;
}

template < class GT, class TDS >
bool
Periodic_3_Delaunay_triangulation_3<GT,TDS>::
is_valid(Cell_handle ch, bool verbose, int level) const
{
  bool error = false;
  if(!Base::is_valid(ch, verbose, level)) {
    error = true;
    if(verbose) {
      std::cerr << "geometrically invalid cell" << std::endl;
      for(int i=0; i<4; i++)
        std::cerr << ch->vertex(i)->point() << ", ";
      std::cerr << std::endl;
    }
  }

  for(Vertex_iterator vit = vertices_begin(); vit != vertices_end(); ++ vit) {
    for(int i=-1; i<=1; i++) {
      for(int j=-1; j<=1; j++) {
        for(int k=-1; k<=1; k++) {
          if(periodic_point(ch,0) == std::make_pair(periodic_point(vit).first,
                                                     periodic_point(vit).second+Offset(i,j,k))
              || periodic_point(ch,1) == std::make_pair(periodic_point(vit).first,
                                                        periodic_point(vit).second+Offset(i,j,k))
              || periodic_point(ch,2) == std::make_pair(periodic_point(vit).first,
                                                        periodic_point(vit).second+Offset(i,j,k))
              || periodic_point(ch,3) == std::make_pair(periodic_point(vit).first,
                                                        periodic_point(vit).second+Offset(i,j,k)) )
            continue;
          if(_side_of_sphere(ch, periodic_point(vit).first,
                              periodic_point(vit).second+Offset(i,j,k),true)
              != ON_UNBOUNDED_SIDE) {
            error = true;
            if(verbose) {
              std::cerr << "Delaunay invalid cell" << std::endl;
              for(int i=0; i<4; i++) {
                Periodic_point pp = periodic_point(ch,i);
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
class Periodic_3_Delaunay_triangulation_3<GT,Tds>::Conflict_tester
{
  // stores a pointer to the triangulation,
  // a point, and an offset
  const Self *t;
  Point p;
  // stores the offset of a point in 27-cover
  mutable Offset o;

public:
  /// Constructor
  Conflict_tester(const Self *_t) : t(_t), p(Point()) {}
  Conflict_tester(const Point &pt, const Self *_t) : t(_t), p(pt) { }

  /** The functor
    *
    * gives true if the circumcircle of c contains p
    */
  bool operator()(const Cell_handle c, const Offset &off) const {
    return (t->_side_of_sphere(c, p, t->combine_offsets(o, off), true)
        == ON_BOUNDED_SIDE);
  }

  bool operator()(const Cell_handle c, const Point& pt,
      const Offset &off) const {
    return (t->_side_of_sphere(c, pt, o + off, true) == ON_BOUNDED_SIDE);
  }

  int compare_weight(Point, Point) const
  {
    return 0;
  }

  bool test_initial_cell(Cell_handle c, const Offset &off) const
  {
    if(!(operator()(c, off)))
      CGAL_triangulation_assertion(false);
    return true;
  }

  void set_point(const Point &_p) {
    p = _p;
  }

  void set_offset(const Offset &off) const {
    o = off;
  }

  const Offset &get_offset() const {
    return o;
  }

  const Point &point() const {
    return p;
  }

};

template < class GT, class Tds>
class Periodic_3_Delaunay_triangulation_3<GT,Tds>::Point_hider
{
public:
  Point_hider() { }

  void set_original_cube (bool) const { }

  template <class InputIterator>
  inline void set_vertices(InputIterator, InputIterator) const { }
  inline void reinsert_vertices(Vertex_handle ) { }
  inline Vertex_handle replace_vertex(Cell_handle c, int index, const Point &) {
    return c->vertex(index);
  }
  inline void hide_point(Cell_handle, const Point &) { }

  inline void hide(Point &, Cell_handle ) const {
    CGAL_triangulation_assertion(false);
  }

  inline void do_hide(const Point &, Cell_handle ) const {
    CGAL_triangulation_assertion(false);
  }
  template < class Tester >
  inline bool replace_vertex(const Point &, Vertex_handle ,
                             const Tester &) const {
    return true;
  }
  template <class Conflict_tester>
  inline void hide_points(Vertex_handle, const Conflict_tester &) {
    // No points to hide in the Delaunay triangulation.
  }
};

#ifndef CGAL_CFG_OUTOFLINE_TEMPLATE_MEMBER_DEFINITION_BUG
template <class GT, class Tds>
template <class TriangulationR3>
struct Periodic_3_Delaunay_triangulation_3<GT,Tds>::Vertex_remover
{
  typedef TriangulationR3      Triangulation_R3;

  typedef typename std::vector<Point>::iterator Hidden_points_iterator;

  typedef Triple < Vertex_handle, Vertex_handle, Vertex_handle > Vertex_triple;

  typedef typename Triangulation_R3::Triangulation_data_structure TDSE;
  typedef typename Triangulation_R3::Cell_handle        CellE_handle;
  typedef typename Triangulation_R3::Vertex_handle      VertexE_handle;
  typedef typename Triangulation_R3::Facet              FacetE;
  typedef typename Triangulation_R3::Finite_cells_iterator Finite_cellsE_iterator;

  typedef Triple< VertexE_handle, VertexE_handle, VertexE_handle > VertexE_triple;

  typedef std::map<Vertex_triple, Facet>  Vertex_triple_Facet_map;
  typedef std::map<Vertex_triple, FacetE> Vertex_triple_FacetE_map;
  typedef typename Vertex_triple_FacetE_map::iterator Vertex_triple_FacetE_map_it;

  Vertex_remover(const Self *t, Triangulation_R3 &tmp_) : _t(t),tmp(tmp_) {}

  const Self *_t;
  Triangulation_R3 &tmp;

  void add_hidden_points(Cell_handle) {
    std::copy(hidden_points_begin(), hidden_points_end(),
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
  std::vector<Point> hidden;
};
#endif //CGAL_CFG_OUTOFLINE_TEMPLATE_MEMBER_DEFINITION_BUG

template < class GT, class TDS >
inline bool
Periodic_3_Delaunay_triangulation_3<GT,TDS>::
is_extensible_triangulation_in_1_sheet_h1() const
{
  if(!is_1_cover()) {
    if(too_long_edge_counter == 0) return true;
    else return false;
  } else {
    typename Geometric_traits::FT longest_edge_squared_length(0);
    Segment s;
    for(Periodic_segment_iterator psit = this->periodic_segments_begin(Base::UNIQUE);
        psit != this->periodic_segments_end(Base::UNIQUE); ++psit) {
      s = this->construct_segment(*psit);
      longest_edge_squared_length = (std::max)(longest_edge_squared_length,
    s.squared_length());
    }
    return (longest_edge_squared_length < edge_length_threshold);
  }
}

template < class GT, class TDS >
inline bool
Periodic_3_Delaunay_triangulation_3<GT,TDS>::
is_extensible_triangulation_in_1_sheet_h2() const
{
  typedef typename Geometric_traits::Construct_circumcenter_3
    Construct_circumcenter;
  typedef typename Geometric_traits::FT FT;
  Construct_circumcenter construct_circumcenter
    = geom_traits().construct_circumcenter_3_object();
  for(Periodic_tetrahedron_iterator tit = this->periodic_tetrahedra_begin(Base::UNIQUE);
       tit != this->periodic_tetrahedra_end(Base::UNIQUE); ++tit) {
    Point cc = construct_circumcenter(
  tit->at(0).first, tit->at(1).first,
  tit->at(2).first, tit->at(3).first,
  tit->at(0).second, tit->at(1).second,
  tit->at(2).second, tit->at(3).second);

    if( !(FT(16)*squared_distance(cc,point(tit->at(0)))
      < (domain().xmax()-domain().xmin())*(domain().xmax()-domain().xmin())) )
      return false;
  }
  return true;
}

template < class GT, class TDS >
inline void
Periodic_3_Delaunay_triangulation_3<GT,TDS>::
copy_multiple_covering(const Periodic_3_Delaunay_triangulation_3<GT,TDS> & tr)
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
  for(Vertex_iterator vit2 = vertices_begin(); vit2 != vertices_end(); ++vit2) {
    if(vit2->offset() != Offset()) {
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
  for(Vertex_iterator vit = vertices_begin(); vit != vertices_end(); ++vit)
    vit->clear_offset();
  for(Vertex_iterator vit = tr.vertices_begin(); vit != tr.vertices_end(); ++vit)
    vit->clear_offset();

  // Build up the too_long_edges container
  too_long_edge_counter = 0;
  too_long_edges.clear();
  for(Vertex_iterator vit = vertices_begin(); vit != vertices_end(); ++vit)
    too_long_edges[vit] = std::list<Vertex_handle>();

  std::pair<Vertex_handle, Vertex_handle> edge_to_add;
  Point p1, p2;
  int i, j;
  for(Edge_iterator eit = edges_begin(); eit != edges_end(); ++eit) {
    if(&*(eit->first->vertex(eit->second)) < &*(eit->first->vertex(eit->third))) {
      i = eit->second; j = eit->third;
    } else {
      i = eit->third; j = eit->second;
    }
    edge_to_add = std::make_pair(eit->first->vertex(i),
                                 eit->first->vertex(j));
    p1 = construct_point(eit->first->vertex(i)->point(),
                         get_offset(eit->first, i));
    p2 = construct_point(eit->first->vertex(j)->point(),
                         get_offset(eit->first, j));
    Vertex_handle v_no = eit->first->vertex(i);
    if(squared_distance(p1,p2) > edge_length_threshold) {
      CGAL_triangulation_assertion(
            find(too_long_edges[v_no].begin(),
                 too_long_edges[v_no].end(),
                 edge_to_add.second) == too_long_edges[v_no].end());
      too_long_edges[v_no].push_back(edge_to_add.second);
      too_long_edge_counter++;
    }
  }
}

template < class GT, class TDS >
std::istream &
operator>> (std::istream& is, Periodic_3_Delaunay_triangulation_3<GT,TDS> &tr)
{
  typedef Periodic_3_Delaunay_triangulation_3<GT,TDS>   P3DT3;
  typedef typename P3DT3::Base                          Base;
  typedef typename P3DT3::Vertex_iterator               Vertex_iterator;
  typedef typename GT::FT FT;
  typedef typename P3DT3::Vertex_handle                 Vertex_handle;

  is >> static_cast<Base&>(tr);

  int i = 0;
  for(Vertex_iterator vi = tr.vertices_begin(); vi != tr.vertices_end(); ++vi) {
    tr.too_long_edges[vi]=std::list<Vertex_handle>();
    ++i;
  }

  tr.edge_length_threshold = FT(0.166) * (tr.domain().xmax()-tr.domain().xmin())
                                       * (tr.domain().xmax()-tr.domain().xmin());
  tr.too_long_edge_counter = tr.find_too_long_edges(tr.too_long_edges);

  CGAL_triangulation_expensive_assertion( tr.is_valid() );
  return is;
}

} // namespace CGAL

#endif // CGAL_PERIODIC_3_DELAUNAY_TRIANGULATION_3_H
