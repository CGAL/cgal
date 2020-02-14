// Copyright (c) 2019 GeometryFactory (France).
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
//
// Author(s)     : Jane Tournois

#ifndef CGAL_INTERNAL_COLLAPSE_SHORT_EDGES_H
#define CGAL_INTERNAL_COLLAPSE_SHORT_EDGES_H

#include <boost/bimap.hpp>
#include <boost/bimap/set_of.hpp>
#include <boost/bimap/multiset_of.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/array.hpp>
#include <boost/unordered_set.hpp>

#include <CGAL/Triangulation_incremental_builder_3.h>

#include <vector>
#include <algorithm>

#include <CGAL/Tetrahedral_remeshing/internal/tetrahedral_remeshing_helpers.h>


namespace CGAL
{
namespace Tetrahedral_remeshing
{
namespace internal
{
  enum Edge_type     { FEATURE, BOUNDARY, INSIDE, MIXTE,
                       NO_COLLAPSE, INVALID, IMAGINARY, MIXTE_IMAGINARY, HULL_EDGE };
  enum Collapse_type { TO_MIDPOINT, TO_V0, TO_V1, IMPOSSIBLE };
  enum Result_type   { VALID,
                       V_PROBLEM, C_PROBLEM, E_PROBLEM,
                       TOPOLOGICAL_PROBLEM, ORIENTATION_PROBLEM, SHARED_NEIGHBOR_PROBLEM };

  template<typename C3t3, typename Visitor>
  class CollapseTriangulation
  {
    typedef typename C3t3::Triangulation                        Tr;
    typedef typename C3t3::Edge                                 Edge;
    typedef typename C3t3::Cell_handle                          Cell_handle;
    typedef typename C3t3::Vertex_handle                        Vertex_handle;
    typedef typename C3t3::Subdomain_index                      Subdomain_index;
    typedef typename C3t3::Triangulation::Point                 Point_3;
    typedef typename C3t3::Triangulation::Geom_traits::Vector_3 Vector_3;

    typedef CGAL::Triangulation_incremental_builder_3<Tr> Builder;

  public:
    CollapseTriangulation(C3t3& c3t3,
                          const Edge& edge,
                          Collapse_type _collapse_type,
                          Visitor& visitor)
    {
      v0_init = edge.first->vertex(edge.second);
      v1_init = edge.first->vertex(edge.third);

      std::vector<Vertex_handle> vertices_to_insert;
      c3t3.triangulation().finite_incident_vertices(v0_init,
                             std::back_inserter(vertices_to_insert));
      vertices_to_insert.push_back(v0_init);
      c3t3.triangulation().finite_incident_vertices(v1_init,
                             std::back_inserter(vertices_to_insert));

      // create incremental builder
      Builder builder(triangulation, true);
      builder.begin_triangulation(3);

      collapse_type = _collapse_type;

      //To add the vertices only once
      for (Vertex_handle vh : vertices_to_insert)
      {
        if (v2v.left.find(vh) == v2v.left.end())
        {
          Vertex_handle new_vh = builder.add_vertex();
          new_vh->set_point(vh->point());
          new_vh->set_dimension(vh->in_dimension());

          v2v.left.insert(std::make_pair(vh, new_vh));
        }
      }

      std::vector<Cell_handle> cells_to_insert;
      c3t3.triangulation().finite_incident_cells(v0_init, std::back_inserter(cells_to_insert));
      c3t3.triangulation().finite_incident_cells(v1_init, std::back_inserter(cells_to_insert));

      //To add the cells only once
      for (Cell_handle ch : cells_to_insert)
      {
        if (c2c.left.find(ch) == c2c.left.end())
        {
          Cell_handle new_ch = builder.add_cell(v2v.left.at(ch->vertex(0)), v2v.left.at(ch->vertex(1)),
                                                v2v.left.at(ch->vertex(2)), v2v.left.at(ch->vertex(3)));
          new_ch->set_subdomain_index(ch->subdomain_index());
          visitor.after_add_cell(ch, new_ch);

          c2c.left.insert(std::make_pair(ch, new_ch));
        }
      }

      // finished
      builder.end_triangulation();
    }

    void update()
    {
      vh0 = v2v.left.at(v0_init);
      vh1 = v2v.left.at(v1_init);

      Cell_handle ch;
      int i0, i1;
      not_an_edge = true;
      if (triangulation.is_edge(vh0, vh1, ch, i0, i1))
      {
        edge = Edge(ch, i0, i1);
        not_an_edge = false;
      }

      to_remove.clear();
      sharing_neighbor.clear();

      typedef typename Tr::Cell_circulator Cell_circulator;
      Cell_circulator circ = triangulation.incident_cells(edge);
      Cell_circulator done = circ;
      do
      {
        to_remove[circ] = true;
        if (circ->neighbor(circ->index(vh0))->has_neighbor(circ->neighbor(circ->index(vh1))))
        {
          sharing_neighbor[circ->neighbor(circ->index(vh0))] = true;
          sharing_neighbor[circ->neighbor(circ->index(vh1))] = true;
        }
      } while (++circ != done);

      collapsed = false;
    }

    Result_type collapse()
    {
      if (not_an_edge)
      {
        std::cout << "LocalTriangulation::Not an edge..." << std::endl;
        return E_PROBLEM;
      }
      else
      {
        Vector_3 v0_new_pos = vec(vh0->point());

        if (collapse_type == TO_MIDPOINT){
          v0_new_pos = v0_new_pos + 0.5 * Vector_3(point(vh0->point()), point(vh1->point()));
        }
        else if (collapse_type == TO_V1){
          v0_new_pos = vec(point(vh1->point()));
        }

        boost::unordered_set<Cell_handle> invalid_cells;

        typedef typename Tr::Cell_circulator Cell_circulator;
        Cell_circulator circ = triangulation.incident_cells(edge);
        Cell_circulator done = circ;

        std::vector<Cell_handle> cells_to_remove;

        //Update the vertex before removing it
        std::vector<Cell_handle> find_incident;
        triangulation.incident_cells(vh0, std::back_inserter(find_incident));

        std::vector<Cell_handle> cells_to_update;
        triangulation.incident_cells(vh1, std::back_inserter(cells_to_update));

//        Result_type r = VALID;
        do
        {
          int v0_id = circ->index(vh0);
          int v1_id = circ->index(vh1);

          Cell_handle n0_ch = circ->neighbor(v0_id);
          Cell_handle n1_ch = circ->neighbor(v1_id);

          int ch_id_in_n0 = n0_ch->index(circ);
          int ch_id_in_n1 = n1_ch->index(circ);

//          if (n0_ch->has_neighbor(n1_ch))
//            r = SHARED_NEIGHBOR_PROBLEM;

          //Update neighbors before removing cell
          n0_ch->set_neighbor(ch_id_in_n0, n1_ch);
          n1_ch->set_neighbor(ch_id_in_n1, n0_ch);

          Subdomain_index si_n0 = n0_ch->subdomain_index();
          Subdomain_index si_n1 = n1_ch->subdomain_index();
          Subdomain_index si = circ->subdomain_index();

          if (si_n0 != si && si_n1 != si)
            return TOPOLOGICAL_PROBLEM;

          if ( triangulation.is_infinite(n0_ch->vertex(ch_id_in_n0))
            && triangulation.is_infinite(n1_ch->vertex(ch_id_in_n1)))
            return TOPOLOGICAL_PROBLEM;

          if ( triangulation.is_infinite(n0_ch)
            && triangulation.is_infinite(n1_ch)
            && !triangulation.is_infinite(circ))
            return TOPOLOGICAL_PROBLEM;

          cells_to_remove.push_back(circ);

          invalid_cells.insert(circ);

        } while (++circ != done);


        vh0->set_point(Point_3(v0_new_pos.x(), v0_new_pos.y(), v0_new_pos.z()));
        vh1->set_point(Point_3(v0_new_pos.x(), v0_new_pos.y(), v0_new_pos.z()));

        Vertex_handle infinite_vertex = triangulation.infinite_vertex();

        bool v0_updated = false;
        for (unsigned int i = 0; i < find_incident.size(); i++)
        {
          const Cell_handle ch = find_incident[i];
          if (invalid_cells.find(ch) == invalid_cells.end()) //valid cell
          {
            if (triangulation.is_infinite(ch))
              infinite_vertex->set_cell(ch);
            else {
              vh0->set_cell(ch);
              v0_updated = true;
            }
          }
        }

        //Update the vertex before removing it
        for (unsigned int i = 0; i < cells_to_update.size(); i++)
        {
          Cell_handle & ch = cells_to_update[i];

          if (invalid_cells.find(ch) == invalid_cells.end()) //valid cell
          {
            ch->set_vertex(ch->index(vh1), vh0);

            if (triangulation.is_infinite(ch))
              infinite_vertex->set_cell(ch);
            else {
              if (!v0_updated) {
                vh0->set_cell(ch);
                v0_updated = true;
              }
            }
          }
        }

        if (!v0_updated){
          std::cout << "CollapseTriangulation::PB i cell not valid!!!" << std::endl;
          return V_PROBLEM;
        }
        triangulation.tds().delete_vertex(vh1);

        //Removing cells
        for (unsigned int i = 0; i < cells_to_remove.size(); i++){
          triangulation.tds().delete_cell(cells_to_remove[i]);
        }

        typedef typename Tr::Finite_cells_iterator Finite_cells_iterator;
        for (Finite_cells_iterator cit = triangulation.finite_cells_begin();
             cit != triangulation.finite_cells_end(); ++cit)
        {
          if (!is_well_oriented(triangulation, cit))
            return ORIENTATION_PROBLEM;
        }

        typedef typename Tr::Cell_iterator Cell_iterator;
        for (Cell_iterator cit = triangulation.cells_begin();
             cit != triangulation.cells_end(); ++cit)
        {
          if (!triangulation.tds().is_valid(cit, true))
            return C_PROBLEM;
        }

        typedef typename Tr::Vertex_iterator Vertex_iterator;
        for (Vertex_iterator vit = triangulation.vertices_begin();
             vit != triangulation.vertices_end(); ++vit)
        {
          if (!triangulation.tds().is_valid(vit, true))
            return V_PROBLEM;
        }

        //int si_nb_vh0 = nb_incident_subdomains(vh0, c3t3);
        //int si_nb_vh1 = nb_incident_subdomains(vh1, c3t3);
        //int vertices_subdomain_nb_vh0 = std::max(si_nb_vh0, si_nb_vh1);
        //bool is_on_hull_vh0 = is_on_convex_hull(vh0, c3t3) || is_on_convex_hull(vh1, c3t3);

        //if( is_valid_for_domains() )
        return VALID;

       // return TOPOLOGICAL_PROBLEM;
      }
    }

  protected:
    Tr triangulation;
    boost::bimap<Vertex_handle, Vertex_handle> v2v;/*vertex of main tr - vertex of collapse tr*/
    boost::bimap<Cell_handle, Cell_handle>     c2c;/*cell of main tr - cell of collapse tr*/

    boost::unordered_map<Cell_handle, bool> to_remove; //default is false
    boost::unordered_map<Cell_handle, bool> sharing_neighbor;//default is false

    Collapse_type collapse_type;

    Vertex_handle v0_init;
    Vertex_handle v1_init;

    Vertex_handle vh0;
    Vertex_handle vh1;

    Edge edge;

    bool collapsed;
    bool not_an_edge;
  };

  template<typename C3t3, typename CellSelector>
  bool topology_test(const typename C3t3::Edge& edge,
                     const C3t3& c3t3,
                     const CellSelector& cell_selector)
  {
    typedef typename C3t3::Vertex_handle Vertex_handle;
    typedef typename C3t3::Cell_handle   Cell_handle;
    typedef typename C3t3::Edge          Edge;
    typedef typename C3t3::Facet         Facet;
    typedef typename C3t3::Triangulation::Facet_circulator Facet_circulator;

    const Vertex_handle v0 = edge.first->vertex(edge.second);
    const Vertex_handle v1 = edge.first->vertex(edge.third);

    // the "topology test" checks that :
    // no incident non-boundary facet has 3 boundary edges
    // no incident boundary facet has 3 feature edges

    Facet_circulator fcirc = c3t3.triangulation().incident_facets(edge);
    Facet_circulator fdone = fcirc;
    do
    {
      if (c3t3.triangulation().is_infinite(fcirc->first))
        continue;

      const Facet& f = *fcirc;
      if (is_boundary(c3t3, f, cell_selector))
        //boundary : check that facet does not have 3 feature edges
      {
        //Get the ids of the opposite vertices
        for (int i = 1; i < 4; i++)
        {
          Vertex_handle vi = f.first->vertex((f.second + i) % 4);
          if (vi != v0 && vi != v1 && nb_incident_subdomains(vi, c3t3) > 1)
          {
            if (is_edge_in_complex(v0, vi, c3t3)
              && is_edge_in_complex(v1, vi, c3t3))
              return false;
          }
        }
      }
      else //non-boundary : check that facet does not have 3 boundary edges
      {
        const Cell_handle circ = f.first;
        const int i = f.second;
        if ( is_boundary(c3t3, Edge(circ, (i + 1) % 4, (i + 2) % 4), cell_selector)
          && is_boundary(c3t3, Edge(circ, (i + 2) % 4, (i + 3) % 4), cell_selector)
          && is_boundary(c3t3, Edge(circ, (i + 3) % 4, (i + 1) % 4), cell_selector))
          return false;
      }
    } while (++fcirc != fdone);

    return true;
  }

  template<typename C3t3, typename CellSelector>
  Collapse_type get_collapse_type(const typename C3t3::Edge& edge,
                                  const C3t3& c3t3,
                                  CellSelector cell_selector)
  {
    bool update_v0 = false;
    bool update_v1 = false;
    get_edge_info(edge, update_v0, update_v1, c3t3, cell_selector);

    if (update_v0 && update_v1) return TO_MIDPOINT;
    else if (update_v0)         return TO_V1;
    else if (update_v1)         return TO_V0;
    else                        return IMPOSSIBLE;
  }

  template<typename C3t3>
  Edge_type get_edge_type(const typename C3t3::Edge& edge,
                          const C3t3& c3t3)
  {
    typedef typename C3t3::Vertex_handle Vertex_handle;
    typedef typename C3t3::Triangulation::Cell_circulator Cell_circulator;
    typedef typename C3t3::Subdomain_index Subdomain_index;

    const Vertex_handle & v0 = edge.first->vertex(edge.second);
    const Vertex_handle & v1 = edge.first->vertex(edge.third);

    int dim0 = c3t3.in_dimension(v0);
    int dim1 = c3t3.in_dimension(v1);

    bool is_v0_on_hull = is_on_convex_hull(v0, c3t3);
    bool is_v1_on_hull = is_on_convex_hull(v1, c3t3);

    if (dim0 == 3 && dim1 == 3)
    {
      if (is_v0_on_hull && is_v1_on_hull)
      {
        Cell_circulator circ = c3t3.triangulation().incident_cells(edge);
        Cell_circulator done = circ;
        do
        {
          if (c3t3.triangulation().is_infinite(circ))
            return HULL_EDGE;
        }
        while (++circ != done);
        return NO_COLLAPSE;
      }
      else if (is_v0_on_hull || is_v1_on_hull)
      {
        return MIXTE_IMAGINARY;
      }
      return INSIDE;
    }

    if (dim0 == 2 && dim1 == 2)
    {
      Cell_circulator circ = c3t3.triangulation().incident_cells(edge);
      Cell_circulator done = circ;

      std::vector<Subdomain_index> indices;
      do
      {
        Subdomain_index current_si = circ->subdomain_index();

        if (std::find(indices.begin(), indices.end(), current_si) == indices.end()){
          indices.push_back(current_si);
        }

        Subdomain_index si_n = circ->neighbor(circ->index(v0))->subdomain_index();
        if (si_n ==
          circ->neighbor(circ->index(v1))->subdomain_index() && si_n != current_si){
          return NO_COLLAPSE;
        }

      }
      while (++circ != done);

      std::size_t nb_si_v0 = nb_incident_subdomains(v0, c3t3);
      std::size_t nb_si_v1 = nb_incident_subdomains(v1, c3t3);

      if (indices.size() >= (std::min)(nb_si_v0, nb_si_v1)){
        return BOUNDARY;
      }

      return NO_COLLAPSE;
    }

    if (dim0 == 3 && dim1 == 2)
    {
      if (is_v0_on_hull)
        return NO_COLLAPSE;
      return MIXTE;
    }

    if (dim1 == 3 && dim0 == 2)
    {
      if (is_v1_on_hull)
        return NO_COLLAPSE;
      return MIXTE;
    }

    //std::cerr << "ERROR : get_edge_type did not return anything valid!" << std::endl;
    return NO_COLLAPSE;
  }

  template<typename C3t3>
  bool is_valid_collapse(const typename C3t3::Edge& edge,
                         const C3t3& c3t3)
  {
    typedef typename C3t3::Vertex_handle Vertex_handle;
    typedef typename C3t3::Cell_handle   Cell_handle;
    typedef typename C3t3::Triangulation::Cell_circulator Cell_circulator;

    const Vertex_handle v0 = edge.first->vertex(edge.second);
    const Vertex_handle v1 = edge.first->vertex(edge.third);

    Cell_circulator circ = c3t3.triangulation().incident_cells(edge);
    Cell_circulator done = circ;
    do
    {
      int v0_id = circ->index(v0);
      int v1_id = circ->index(v1);

      Cell_handle n0_ch = circ->neighbor(v0_id);
      Cell_handle n1_ch = circ->neighbor(v1_id);

      if ( n0_ch->has_vertex(v0)
        || n1_ch->has_vertex(v1)
        || n0_ch->has_neighbor(n1_ch))
        return false;
    }
    while (++circ != done);

    return true;
  }

  template<typename C3t3, typename CellSelector>
  bool is_valid_collapse(const typename C3t3::Edge& edge,
                         const Collapse_type& collapse_type,
                         const typename C3t3::Triangulation::Point& new_pos,
                         const C3t3& c3t3,
                         const bool /*protect_boundaries*/,
                         CellSelector cell_selector)
  {
    typedef typename C3t3::Vertex_handle        Vertex_handle;
    typedef typename C3t3::Cell_handle          Cell_handle;
    typedef typename C3t3::Triangulation::Point Point;

    Vertex_handle v0 = edge.first->vertex(edge.second);
    Vertex_handle v1 = edge.first->vertex(edge.third);

    ////about protection of boundaries
    //if (protect_boundaries)
    //{
    //  if (c3t3.is_in_complex(edge)
    //   || helpers::is_boundary(c3t3, edge, cell_selector))
    //    return false;
    //}
    //we need to check that surfaces are not broken anyhow
    bool v0_boundary = is_boundary_vertex(v0, c3t3, cell_selector);
    bool v1_boundary = is_boundary_vertex(v1, c3t3, cell_selector);
    if (collapse_type == TO_V0 && v1_boundary && !v0_boundary)
      return false;
    if (collapse_type == TO_V1 && v0_boundary && !v1_boundary)
      return false;
    if (collapse_type == TO_MIDPOINT && (v0_boundary ^ v1_boundary))//both or none to allow collapse
      return false;
    
    std::vector<Cell_handle> cells_to_check;
    if (collapse_type == TO_V1 || collapse_type == TO_MIDPOINT)
    {
      c3t3.triangulation().finite_incident_cells(v0,
                             std::back_inserter(cells_to_check));

      for (std::size_t i = 0; i < cells_to_check.size(); i++)
      {
        const Cell_handle& ch = cells_to_check[i];
        if (!ch->has_vertex(v1))
        {
          //check orientation
          boost::array<Point, 4> pts = { ch->vertex(0)->point(),
                                         ch->vertex(1)->point(),
                                         ch->vertex(2)->point(),
                                         ch->vertex(3)->point()};
          pts[ch->index(v0)] = new_pos;
          if (CGAL::orientation(point(pts[0]), point(pts[1]),
                                point(pts[2]), point(pts[3])) != CGAL::POSITIVE)
            return false;
        }
      }
      cells_to_check.clear();
    }
    else if (collapse_type == TO_V0 || collapse_type == TO_MIDPOINT)
    {
      c3t3.triangulation().finite_incident_cells(v1,
                             std::back_inserter(cells_to_check));

      for (std::size_t i = 0; i < cells_to_check.size(); i++)
      {
        const Cell_handle& ch = cells_to_check[i];
        if (!ch->has_vertex(v0))
        {
          //check orientation
          //check orientation
          boost::array<Point, 4> pts = { ch->vertex(0)->point(),
                                         ch->vertex(1)->point(),
                                         ch->vertex(2)->point(),
                                         ch->vertex(3)->point() };
          pts[ch->index(v1)] = new_pos;
          if (CGAL::orientation(point(pts[0]), point(pts[1]),
                                point(pts[2]), point(pts[3])) != CGAL::POSITIVE)
            return false;
        }
      }
      cells_to_check.clear();
    }

    return is_valid_collapse(edge, c3t3);
  }

  template<typename C3t3, typename SqLengthMap>
  bool are_edge_lengths_valid(const typename C3t3::Vertex_handle v1,
                              const typename C3t3::Vertex_handle v2,
                              const C3t3& c3t3,
                              const typename C3t3::Triangulation::Point& new_pos,
                              SqLengthMap& edges_sqlength,
                              const typename C3t3::Triangulation::Geom_traits::FT& sqhigh,
                              const bool /* adaptive */ = false)
  {
    //SqLengthMap::key_type is Vertex_handle
    //SqLengthMap::value_type is double
    typedef typename C3t3::Triangulation::Geom_traits::FT FT;
    typedef typename C3t3::Edge                           Edge;
    typedef typename C3t3::Vertex_handle                  Vertex_handle;

    std::vector<Edge> inc_edges;
    c3t3.triangulation().finite_incident_edges(v1,
                           std::back_inserter(inc_edges));

    for (std::size_t i = 0; i < inc_edges.size(); i++)
    {
      const Edge& ei = inc_edges[i];

      Vertex_handle ivh = ei.first->vertex(ei.second);
      if (ivh == v1)
        ivh = ei.first->vertex(ei.third);

      if (v2 != ivh && edges_sqlength.find(ivh) == edges_sqlength.end())
      {
        FT sqlen_i = CGAL::squared_distance(new_pos, ivh->point());

        //if (adaptive){
        //  if (is_boundary_edge(ei) || is_hull_edge(ei)){
        //    if (sqlen_i > split_length)
        //      return false;
        //  }
        //  else if (sqlen_i > 4.*getAimedLength(ei, aimed_length) / 3.){// && is_in_complex(ei)  ){
        //    return false;
        //  }
        //}
        //else {
        if (sqlen_i > sqhigh) {
          return false;
        }
        //}

        edges_sqlength[ivh] = sqlen_i;
      }
    }

    return true;
  }

  template<typename C3t3, typename SqLengthMap>
  bool are_edge_lengths_valid(const typename C3t3::Edge& edge,
                              const C3t3& c3t3,
                              const Collapse_type& collapse_type,
                              const typename C3t3::Triangulation::Point& new_pos,
                              SqLengthMap& edges_sqlength,
                              const typename C3t3::Triangulation::Geom_traits::FT& sqhigh,
                              const bool adaptive = false)
  {
    //SqLengthMap::key_type is Vertex_handle
    //SqLengthMap::value_type is double

    typedef typename C3t3::Vertex_handle Vertex_handle;

    edges_sqlength.clear();
    Vertex_handle v0 = edge.first->vertex(edge.second);
    Vertex_handle v1 = edge.first->vertex(edge.third);

    if (collapse_type == TO_V1 || collapse_type == TO_MIDPOINT)
    {
      if (!are_edge_lengths_valid(v0, v1, c3t3, new_pos,
                                  edges_sqlength, sqhigh, adaptive))
        return false;
    }
    else if (collapse_type == TO_V0 || collapse_type == TO_MIDPOINT)
    {
      if (!are_edge_lengths_valid(v1, v0, c3t3, new_pos,
                                  edges_sqlength, sqhigh, adaptive))
        return false;
    }
    return true;
  }

  template<typename C3t3>
  void merge_surface_patch_indices(typename C3t3::Facet& f1,
                                   typename C3t3::Facet& f2,
                                   C3t3& c3t3)
  {
    const bool in_cx_f1 = c3t3.is_in_complex(f1);
    const bool in_cx_f2 = c3t3.is_in_complex(f2);

    if (in_cx_f1 && !in_cx_f2)
    {
      typename C3t3::Surface_patch_index patch = c3t3.surface_patch_index(f1);
      c3t3.remove_from_complex(f1);
      c3t3.add_to_complex(f1, patch);
      c3t3.add_to_complex(f2, patch);
    }
    else if (in_cx_f2 && !in_cx_f1)
    {
      typename C3t3::Surface_patch_index patch = c3t3.surface_patch_index(f2);
      c3t3.remove_from_complex(f2);
      c3t3.add_to_complex(f1, patch);
      c3t3.add_to_complex(f2, patch);
    }
    else
    {
      CGAL_assertion(
        //f1 and f2 are not both in complex
        !(in_cx_f1 && in_cx_f2)
        // unless they are on the same surface
        || c3t3.surface_patch_index(f1) == c3t3.surface_patch_index(f2));
    }
  }

  template<typename C3t3>
  typename C3t3::Vertex_handle
  collapse(const typename C3t3::Cell_handle ch,
           const int to, const int from,
           C3t3& c3t3)
  {
    typedef typename C3t3::Triangulation Tr;
    typedef typename C3t3::Vertex_handle Vertex_handle;
    typedef typename C3t3::Cell_handle   Cell_handle;
    typedef typename C3t3::Facet         Facet;
    typedef typename Tr::Cell_circulator Cell_circulator;

    Tr& tr = c3t3.triangulation();

    Vertex_handle vh0 = ch->vertex(to);
    Vertex_handle vh1 = ch->vertex(from);

    std::vector<Cell_handle> cells_to_remove;

    //Update the vertex before removing it
    std::vector<Cell_handle> find_incident;
    tr.incident_cells(vh0, std::back_inserter(find_incident));

    std::vector<Cell_handle> cells_to_update;
    tr.incident_cells(vh1, std::back_inserter(cells_to_update));

//    if (vh1->in_dimension() == 2 && c3t3.is_in_complex(vh1))
//      std::cout << "Collapsing a feature vertex!!!!!!" << std::endl;

    boost::unordered_set<Cell_handle> invalid_cells;
    bool valid = true;
    Cell_circulator circ = tr.incident_cells(ch, to, from);
    Cell_circulator done = circ;
    do
    {
      const int v0_id = circ->index(vh0);
      const int v1_id = circ->index(vh1);

      Cell_handle n0_ch = circ->neighbor(v0_id);
      Cell_handle n1_ch = circ->neighbor(v1_id);

      const int ch_id_in_n0 = n0_ch->index(circ);
      const int ch_id_in_n1 = n1_ch->index(circ);

      //Merge surface patch indices
      merge_surface_patch_indices(Facet(n0_ch, ch_id_in_n0),
                                  Facet(n1_ch, ch_id_in_n1),
                                  c3t3);

      //Update neighbors before removing cell
      n0_ch->set_neighbor(ch_id_in_n0, n1_ch);
      n1_ch->set_neighbor(ch_id_in_n1, n0_ch);

      //Update vertices cell pointer
      //if( !triangulation.is_infinite( n0_ch ) )
      int nb_on_boundary_n0 = 0;
      for (int i = 0; i < 3; i++)
      {
        int vid = Tr::vertex_triple_index(ch_id_in_n0, i);
        n0_ch->vertex(vid)->set_cell(n0_ch);
        if (c3t3.in_dimension(n0_ch->vertex(vid)))
          nb_on_boundary_n0++;
      }
      //else
      int nb_on_boundary_n1 = 0;
      for (int i = 0; i < 3; i++)
      {
        int vid = Tr::vertex_triple_index(ch_id_in_n1, i);
        n1_ch->vertex(vid)->set_cell(n1_ch);
        if (c3t3.in_dimension(n1_ch->vertex(vid)))
          nb_on_boundary_n1++;
      }

      if ( tr.is_infinite(n0_ch->vertex(ch_id_in_n0))
        && tr.is_infinite(n1_ch->vertex(ch_id_in_n1)))
        return Vertex_handle();

      cells_to_remove.push_back(circ);

      invalid_cells.insert(circ);

    } while (++circ != done);

    const Vertex_handle infinite_vertex = tr.infinite_vertex();

    bool v0_updated = false;
    for (const Cell_handle ch : find_incident)
    {
      if (invalid_cells.find(ch) == invalid_cells.end())//valid cell
      {
        if (tr.is_infinite(ch))
          infinite_vertex->set_cell(ch);
        //else {
        vh0->set_cell(ch);
        v0_updated = true;
        //}
      }
    }

    //Update the vertex before removing it
    for (const Cell_handle ch : cells_to_update)
    {
      if (invalid_cells.find(ch) == invalid_cells.end()) //valid cell
      {
        ch->set_vertex(ch->index(vh1), vh0);

        if (tr.is_infinite(ch))
          infinite_vertex->set_cell(ch);
        //else {
        if (!v0_updated) {
          vh0->set_cell(ch);
          v0_updated = true;
        }
        //}
      }
    }

    if (!v0_updated)
      std::cout << "PB i cell not valid!!!" << std::endl;

    // Delete vertex
    c3t3.triangulation().tds().delete_vertex(vh1);

    // Delete cells
    for (Cell_handle cell_to_remove : cells_to_remove)
    {
      if (cell_to_remove->subdomain_index() > 0)
        c3t3.remove_from_complex(cell_to_remove);
      c3t3.triangulation().tds().delete_cell(cell_to_remove);
    }

    if (!valid){
      std::cout << "Global triangulation collapse bug!!" << std::endl;
      return Vertex_handle();
    }

    return vh0;
  }


  template<typename C3t3>
  typename C3t3::Vertex_handle collapse(typename C3t3::Edge& edge,
                                        const Collapse_type& collapse_type,
                                        C3t3& c3t3)
  {
    typedef typename C3t3::Vertex_handle Vertex_handle;
    typedef typename C3t3::Triangulation::Point Point_3;

    Vertex_handle vh0 = edge.first->vertex(edge.second);
    Vertex_handle vh1 = edge.first->vertex(edge.third);

    int dim_vh0 = c3t3.in_dimension(vh0);
    int dim_vh1 = c3t3.in_dimension(vh1);

    Vertex_handle vh = Vertex_handle();

    Point_3 p0 = vh0->point();
    Point_3 p1 = vh1->point();

    //Collapse at mid point
    if (collapse_type == TO_MIDPOINT)
    {
      Point_3 new_position(CGAL::midpoint(point(vh0->point()), point(vh1->point())));
      vh0->set_point(new_position);
      vh1->set_point(new_position);

      vh = collapse(edge.first, edge.second, edge.third, c3t3);
      c3t3.set_dimension(vh, (std::min)(dim_vh0, dim_vh1));
    }
    else //Collapse at vertex
    {
      if (collapse_type == TO_V1)
      {
        vh0->set_point(p1);
        vh = collapse(edge.first, edge.third, edge.second, c3t3);
        c3t3.set_dimension(vh, (std::min)(dim_vh0, dim_vh1));
      }
      else //Collapse at v0
      {
        if (collapse_type == TO_V0)
        {
          vh1->set_point(p0);
          vh = collapse(edge.first, edge.second, edge.third, c3t3);
          c3t3.set_dimension(vh, (std::min)(dim_vh0, dim_vh1));
        }
        else
          CGAL_assertion(false);
      }
    }
    return vh;
  }

  template<typename C3t3, typename CellSelector, typename Visitor>
  typename C3t3::Vertex_handle collapse_edge(typename C3t3::Edge& edge,
                     C3t3& c3t3,
                     const typename C3t3::Triangulation::Geom_traits::FT& sqhigh,
                     const bool protect_boundaries,
                     CellSelector cell_selector,
                     Visitor& visitor)
  {
    typedef typename C3t3::Triangulation   Tr;
    typedef typename Tr::Point             Point;
    typedef typename Tr::Vertex_handle     Vertex_handle;

    Vertex_handle v0 = edge.first->vertex(edge.second);
    Vertex_handle v1 = edge.first->vertex(edge.third);

    Collapse_type collapse_type = get_collapse_type(edge, c3t3, cell_selector);
    Edge_type edge_type = get_edge_type(edge, c3t3);

    if (collapse_type != IMPOSSIBLE && edge_type != NO_COLLAPSE)
    {
      Point new_pos;
      switch(collapse_type)
      {
      case TO_V0:
        new_pos = v0->point(); break;
      case TO_V1:
        new_pos = v1->point(); break;
      default:
        CGAL_assertion(collapse_type == TO_MIDPOINT);
        new_pos = Point(CGAL::midpoint(point(v0->point()), point(v1->point())));
      }

      boost::unordered_map<Vertex_handle, double> edges_sqlength_after_collapse;
      if (is_valid_collapse(edge, collapse_type, new_pos, c3t3,
                            protect_boundaries, cell_selector))
      {
        if (are_edge_lengths_valid(edge, c3t3, collapse_type, new_pos,
                                   edges_sqlength_after_collapse, sqhigh
                                   /*, adaptive = false*/))
        {
          CollapseTriangulation<C3t3, Visitor> local_tri(c3t3, edge, collapse_type, visitor);
          local_tri.update();

          Result_type res = local_tri.collapse();
          if (res == VALID)
            return collapse(edge, collapse_type, c3t3);
        }
      }
    }

    return Vertex_handle();
  }

  template<typename C3T3, typename CellSelector>
  bool can_be_collapsed(const typename C3T3::Edge& e,
                        const C3T3& c3t3,
                        const bool protect_boundaries,
                        CellSelector cell_selector)
  {
    if (is_outside(e, c3t3, cell_selector))
      return false;

    if (protect_boundaries)
    {
      if (c3t3.is_in_complex(e))
        return false;
      else if (is_boundary(c3t3, e, cell_selector))
        return false;

#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
      if (!is_internal(e, c3t3, cell_selector))
      {
        std::cerr << "e is not inside!?" << std::endl;
        typename C3T3::Vertex_handle v1 = e.first->vertex(e.second);
        typename C3T3::Vertex_handle v2 = e.first->vertex(e.third);
        std::cerr << v1->point() << " " << v2->point() << std::endl;
      }
#endif

      CGAL_assertion(is_internal(e, c3t3, cell_selector));
      return true;
    }
    else
    {
      return true;
    }
  }

  template<typename C3T3, typename CellSelector, typename Visitor>
  void collapse_short_edges(C3T3& c3t3,
    const typename C3T3::Triangulation::Geom_traits::FT& low,
    const typename C3T3::Triangulation::Geom_traits::FT& high,
    const bool protect_boundaries,
    CellSelector cell_selector,
    Visitor& visitor)
  {
    typedef typename C3T3::Triangulation       T3;
    typedef typename T3::Cell_handle           Cell_handle;
    typedef typename T3::Edge                  Edge;
    typedef typename T3::Finite_edges_iterator Finite_edges_iterator;
    typedef typename T3::Vertex_handle         Vertex_handle;
    typedef typename std::pair<Vertex_handle, Vertex_handle> Edge_vv;

    typedef typename T3::Geom_traits     Gt;
    typedef typename T3::Geom_traits::FT FT;
    typedef boost::bimap<
      boost::bimaps::set_of<Edge_vv>,
      boost::bimaps::multiset_of<FT, std::less<FT> > >  Boost_bimap;
    typedef typename Boost_bimap::value_type            short_edge;

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
    std::cout << "Collapse short edges (" << low << ", " << high << ")...";
    std::cout.flush();
    std::size_t nb_collapses = 0;
#endif
    const FT sq_low = low*low;
    const FT sq_high = high*high;

    //collect long edges
    T3& tr = c3t3.triangulation();
    Boost_bimap short_edges;
    for (Finite_edges_iterator eit = tr.finite_edges_begin();
         eit != tr.finite_edges_end(); ++eit)
    {
      const Edge& e = *eit;
      if (!can_be_collapsed(e, c3t3, protect_boundaries, cell_selector))
        continue;

      typename Gt::Compute_squared_length_3 sql
        = tr.geom_traits().compute_squared_length_3_object();
      FT sqlen = sql(tr.segment(e));
      if (sqlen < sq_low)
        short_edges.insert(short_edge(make_vertex_pair<T3>(e), sqlen));
    }

#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
    debug::dump_edges(short_edges, "short_edges.polylines.txt");

    std::ofstream short_success("short_collapse_success.polylines.txt");
    std::ofstream short_fail("short_collapse_fail.polylines.txt");
    std::ofstream short_cancel("short_collapse_canceled.polylines.txt");
#endif

    while(!short_edges.empty())
    {
      //the edge with shortest length
      typename Boost_bimap::right_map::iterator eit = short_edges.right.begin();
      Edge_vv e = eit->second;
      short_edges.right.erase(eit);

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE_PROGRESS
      std::cout << "\rCollapse... (" << short_edges.left.size() << " short edges, ";
      std::cout << nb_collapses << " collapses)";
      std::cout.flush();
#endif
      Cell_handle cell;
      int i1, i2;
      if ( tr.tds().is_vertex(e.first)
        && tr.tds().is_vertex(e.second)
        && tr.tds().is_edge(e.first, e.second, cell, i1, i2)
        && tr.segment(Edge(cell, i1, i2)).squared_length() < sq_low)
      {
#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
        const typename T3::Point p1 = e.first->point();
        const typename T3::Point p2 = e.second->point();
#endif

        Edge edge(cell, i1, i2);

        if (!can_be_collapsed(edge, c3t3, protect_boundaries, cell_selector))
        {
#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
            short_cancel << "2 " << point(p1) << " " << point(p2) << std::endl;
#endif
          continue;
        }
#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
        Vertex_handle vh =
#endif
        collapse_edge(edge, c3t3, sq_high,
                      protect_boundaries, cell_selector, 
                      visitor);

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
        if (vh != Vertex_handle())
          ++nb_collapses;
#endif
#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
        if (vh != Vertex_handle())
          short_success << "2 " << point(p1) << " " << point(p2) << std::endl;
        else
          short_fail << "2 " << point(p1) << " " << point(p2) << std::endl;
#endif
      }
    }//end loop on short_edges
#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
    short_success.close();
    short_fail.close();
#endif

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
    std::cout << " done (" << nb_collapses << " collapses)." << std::endl;
#endif
  }
}
}
}

#endif // CGAL_INTERNAL_COLLAPSE_SHORT_EDGES_H
