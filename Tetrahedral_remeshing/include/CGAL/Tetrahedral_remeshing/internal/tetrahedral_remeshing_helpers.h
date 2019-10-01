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

#ifndef CGAL_INTERNAL_TET_REMESHING_HELPERS_H
#define CGAL_INTERNAL_TET_REMESHING_HELPERS_H

#include <utility>

#include <CGAL/Point_3.h>
#include <CGAL/Weighted_point_3.h>
#include <CGAL/Vector_3.h>

namespace CGAL
{
namespace Tetrahedral_remeshing
{
  enum Subdomain_relation { EQUAL, DIFFERENT, INCLUDED, INCLUDES };
  enum Sliver_removal_result { INVALID_ORIENTATION, INVALID_CELL, INVALID_VERTEX,
    NOT_FLIPPABLE, EDGE_PROBLEM, VALID_FLIP, NO_BEST_CONFIGURATION, EXISTING_EDGE };

  template<typename K>
  CGAL::Point_3<K> point(const CGAL::Point_3<K>& p)
  {
    return p;
  }
  template<typename K>
  CGAL::Point_3<K> point(const CGAL::Weighted_point_3<K>& wp)
  {
    typename K::Construct_point_3 pt = K().construct_point_3_object();
    return pt(wp);
  }

  template <typename K>
  CGAL::Vector_3<K> vec(const CGAL::Point_3<K>& p)
  {
    typename K::Construct_vector_3 v = K().construct_vector_3_object();
    return v(CGAL::ORIGIN, p);
  }
  template <typename K>
  CGAL::Vector_3<K> vec(const CGAL::Weighted_point_3<K>& wp)
  {
    return vec(point(wp));
  }

  
  const int indices_table[4][3] = { { 3, 1, 2 },
                                    { 3, 2, 0 },
                                    { 3, 0, 1 },
                                    { 2, 1, 0 } };

  int indices(const int& i, const int& j)
  {
    CGAL_assertion(i >= 0 && i < 4);
    CGAL_assertion(j >= 0 && j < 3);
    return indices_table[i][j];
  }

  template<typename Gt>
  typename Gt::FT dihedral_angle(const CGAL::Point_3<Gt>& p,
    const CGAL::Point_3<Gt>& q,
    const CGAL::Point_3<Gt>& r,
    const CGAL::Point_3<Gt>& s)
  {
    return Gt().compute_approximate_dihedral_angle_3_object()(p, q, r, s);
  }

  template<typename Gt>
  typename Gt::FT min_dihedral_angle(const CGAL::Point_3<Gt>& p,
    const CGAL::Point_3<Gt>& q,
    const CGAL::Point_3<Gt>& r,
    const CGAL::Point_3<Gt>& s)
  {
    typedef typename Gt::FT FT;
    FT a = CGAL::abs(dihedral_angle(p, q, r, s));
    FT min_dh = a;

    a = CGAL::abs(dihedral_angle(p, r, q, s));
    min_dh = (std::min)(a, min_dh);

    a = CGAL::abs(dihedral_angle(p, s, q, r));
    min_dh = (std::min)(a, min_dh);

    a = CGAL::abs(dihedral_angle(q, r, p, s));
    min_dh = (std::min)(a, min_dh);

    a = CGAL::abs(dihedral_angle(q, s, p, r));
    min_dh = (std::min)(a, min_dh);

    a = CGAL::abs(dihedral_angle(r, s, p, q));
    min_dh = (std::min)(a, min_dh);

    return min_dh;
  }

  template<typename Gt, typename VertexHandle>
  typename Gt::FT min_dihedral_angle(VertexHandle v0,
    VertexHandle v1,
    VertexHandle v2,
    VertexHandle v3)
  {
    return min_dihedral_angle(point(v0->point()),
      point(v1->point()),
      point(v2->point()),
      point(v3->point()));
  }

  template<typename Gt, typename CellHandle>
  typename Gt::FT min_dihedral_angle(CellHandle c)
  {
    return min_dihedral_angle(point(c->vertex(0)->point()),
      point(c->vertex(1)->point()),
      point(c->vertex(2)->point()),
      point(c->vertex(3)->point()));
  }

  template<typename Tr>
  std::pair<typename Tr::Vertex_handle, typename Tr::Vertex_handle>
    make_vertex_pair(const typename Tr::Edge& e)
  {
    typedef typename Tr::Vertex_handle Vertex_handle;
    Vertex_handle v1 = e.first->vertex(e.second);
    Vertex_handle v2 = e.first->vertex(e.third);
    if (v2 < v1) std::swap(v1, v2);

    return std::make_pair(v1, v2);
  }

  template<typename Vh>
  std::pair<Vh, Vh> make_vertex_pair(const Vh v1, const Vh v2)
  {
    if (v2 < v1) return std::make_pair(v2, v1);
    else         return std::make_pair(v1, v2);
  }

  template<typename Vh>
  CGAL::Triple<Vh, Vh, Vh> make_vertex_triple(const Vh vh0, const Vh vh1, const Vh vh2)
  {
    CGAL::Triple<Vh, Vh, Vh> ft(vh0, vh1, vh2);
    if (ft.template get<1>() < ft.template get<0>()) std::swap(ft.template get<0>(), ft.template get<1>());
    if (ft.template get<2>() < ft.template get<1>()) std::swap(ft.template get<1>(), ft.template get<2>());
    if (ft.template get<1>() < ft.template get<0>()) std::swap(ft.template get<0>(), ft.template get<1>());
    return ft;
  }

  template<typename VertexHandle>
  bool is_on_feature(const VertexHandle v)
  {
    return (v->in_dimension() == 1);
  }

  template<typename CellHandle>
  CGAL::Orientation orientation(const CellHandle ch)
  {
    return CGAL::orientation(point(ch->vertex(0)->point()),
                             point(ch->vertex(1)->point()),
                             point(ch->vertex(2)->point()),
                             point(ch->vertex(3)->point()));
  }

  template<typename CellHandle>
  bool is_well_oriented(const CellHandle ch)
  {
    return CGAL::POSITIVE == orientation(ch);
  }

  template<typename VertexHandle>
  bool is_well_oriented(const VertexHandle v0, const VertexHandle v1,
    const VertexHandle v2, const VertexHandle v3)
  {
    return CGAL::POSITIVE == CGAL::orientation(point(v0->point()),
                                               point(v1->point()),
                                               point(v2->point()),
                                               point(v3->point()));
  }


  template<typename C3T3, typename CellSelector>
  bool is_boundary(const C3T3& c3t3,
    const typename C3T3::Triangulation::Edge& e,
    CellSelector cell_selector)
  {
    typedef typename C3T3::Triangulation   Tr;
    typedef typename Tr::Facet_circulator  Facet_circulator;
    typedef typename Tr::Facet             Facet;

    Facet_circulator fcirc = c3t3.triangulation().incident_facets(e);
    Facet_circulator fend = fcirc;
    std::vector<Facet> boundary_facets;

    do
    {
      Facet f = *fcirc;
      if (c3t3.is_in_complex(f))
        return true;
      else if (cell_selector(f.first) // XOR
        ^ cell_selector(f.first->neighbor(f.second)))
        return true;
      else if (c3t3.triangulation().is_infinite(f) //XOR
        ^ c3t3.triangulation().is_infinite(f.first->neighbor(f.second)))
        return true;

      ++fcirc;
    } while (fcirc != fend);

    return false;
  }

  template<typename C3t3, typename CellSelector>
  bool is_boundary_edge(const typename C3t3::Vertex_handle& v0,
    const typename C3t3::Vertex_handle& v1,
    const C3t3& c3t3,
    CellSelector cell_selector)
  {
    typedef typename C3t3::Edge        Edge;
    typedef typename C3t3::Cell_handle Cell_handle;

    Cell_handle cell;
    int i0, i1;
    if (c3t3.triangulation().tds().is_edge(v0, v1, cell, i0, i1))
      return is_boundary(c3t3, Edge(cell, i0, i1), cell_selector);
    else
      return false;
  }

  template<typename C3t3, typename CellSelector>
  bool is_boundary_vertex(const typename C3t3::Vertex_handle& v,
    const C3t3& c3t3,
    CellSelector cell_selector)
  {
    typedef typename C3t3::Facet Facet;
    std::vector<Facet> facets;
    c3t3.triangulation().incident_facets(v, std::back_inserter(facets));

    BOOST_FOREACH(Facet f, facets)
    {
      if (c3t3.is_in_complex(f))
        return true;
      if (cell_selector(f.first) ^ cell_selector(f.first->neighbor(f.second)))
        return true;
    }
    return false;
  }

  template<typename C3t3, typename CellSelector>
  bool is_edge_in_complex(const typename C3t3::Vertex_handle& v0,
    const typename C3t3::Vertex_handle& v1,
    const C3t3& c3t3,
    CellSelector /*cell_selector*/)
  {
    typedef typename C3t3::Edge        Edge;
    typedef typename C3t3::Cell_handle Cell_handle;

    Cell_handle cell;
    int i0, i1;
    if (c3t3.triangulation().tds().is_edge(v0, v1, cell, i0, i1))
      return c3t3.is_in_complex(Edge(cell, i0, i1));
    else
      return false;
  }

  template<typename C3t3, typename CellSelector>
  bool topology_test(const typename C3t3::Edge& edge,
    const C3t3& c3t3,
    CellSelector cell_selector)
  {
    typedef typename C3t3::Vertex_handle Vertex_handle;
    typedef typename C3t3::Triangulation::Facet_circulator Facet_circulator;
    typedef typename C3t3::Subdomain_index Subdomain_index;

    Vertex_handle v0 = edge.first->vertex(edge.second);
    Vertex_handle v1 = edge.first->vertex(edge.third);

    Facet_circulator fcirc = c3t3.triangulation().incident_facets(edge);
    Facet_circulator fdone = fcirc;
    do
    {
      if (c3t3.triangulation().is_infinite(fcirc->first))
        continue;

      Subdomain_index si_circ = fcirc->first->subdomain_index();
      Subdomain_index si_neigh = fcirc->first->neighbor(fcirc->second)->subdomain_index();
      if (si_circ == si_neigh)
      {
        //Get the ids of the opposite vertices
        for (int i = 1; i < 4; i++)
        {
          Vertex_handle vi = fcirc->first->vertex((fcirc->second + i) % 4);
          if (vi != v0 && vi != v1 && nb_incident_subdomains(vi, c3t3) > 1)
          {
            if (is_edge_in_complex(v0, vi, c3t3, cell_selector)
              && is_edge_in_complex(v1, vi, c3t3, cell_selector))
              return false;
          }
        }
      }
    } while (++fcirc != fdone);

    return true;
  }

  template<typename C3t3>
  Subdomain_relation compare_subdomains(typename C3t3::Vertex_handle v0,
    typename C3t3::Vertex_handle v1,
    const C3t3& c3t3)
  {
    typedef typename C3t3::Subdomain_index Subdomain_index;

    std::vector<Subdomain_index> subdomains_v0;
    incident_subdomains(v0, c3t3, std::back_inserter(subdomains_v0));
    std::sort(subdomains_v0.begin(), subdomains_v0.end());

    std::vector<Subdomain_index> subdomains_v1;
    incident_subdomains(v1, c3t3, std::back_inserter(subdomains_v1));
    std::sort(subdomains_v1.begin(), subdomains_v1.end());

    if (subdomains_v0.size() == subdomains_v1.size())
    {
      for (unsigned int i = 0; i < subdomains_v0.size(); i++)
        if (subdomains_v0[i] != subdomains_v1[i])
          return DIFFERENT;
      return EQUAL;
    }
    else
    {
      std::vector<Subdomain_index>
        intersection((std::min)(subdomains_v0.size(), subdomains_v1.size()), -1);
      typename std::vector<Subdomain_index>::iterator
        end_it = std::set_intersection(subdomains_v0.begin(), subdomains_v0.end(),
          subdomains_v1.begin(), subdomains_v1.end(),
          intersection.begin());
      std::ptrdiff_t intersection_size = (end_it - intersection.begin());

      if (subdomains_v0.size() > subdomains_v1.size()
        && intersection_size == std::ptrdiff_t(subdomains_v1.size()))
      {
        return INCLUDES;
      }
      else if (intersection_size == std::ptrdiff_t(subdomains_v0.size())) {
        return INCLUDED;
      }
    }
    return DIFFERENT;
  }



  template<typename C3t3, typename CellSelector>
  void get_edge_info(const typename C3t3::Edge& edge,
    bool& update_v0,
    bool& update_v1,
    const C3t3& c3t3,
    CellSelector cell_selector)
  {
    typedef typename C3t3::Vertex_handle Vertex_handle;

    Vertex_handle v0 = edge.first->vertex(edge.second);
    Vertex_handle v1 = edge.first->vertex(edge.third);

    int dim0 = c3t3.in_dimension(v0);
    int dim1 = c3t3.in_dimension(v1);

    std::size_t nb_si_v0 = nb_incident_subdomains(v0, c3t3);
    std::size_t nb_si_v1 = nb_incident_subdomains(v1, c3t3);

    update_v0 = false;
    update_v1 = false;

    bool is_v0_on_hull = is_on_hull(v0, c3t3);
    bool is_v1_on_hull = is_on_hull(v1, c3t3);

    //Same type imaginary or inside vertices
    if (dim0 == 3 && dim1 == 3)
    {
      if (is_v0_on_hull && is_v1_on_hull)//both endvertices are on hull
      {
        if (is_on_hull(edge, c3t3)) //edge also is on hull
        {
          update_v0 = true;
          update_v1 = true;
        }
      }
      else
      {
        if (!is_v0_on_hull) //v0 not on hull
          update_v0 = true;
        if (!is_v1_on_hull) //v1 not on hull
          update_v1 = true;
      }
      return;
    }
    //Feature edge case
    if (nb_si_v0 > 2 && nb_si_v1 > 2)
    {
      if (c3t3.is_in_complex(edge))
      {
        if (!topology_test(edge, c3t3, cell_selector))
          return;

        if (nb_si_v0 > nb_si_v1) {
          update_v1 = true;
        }
        else if (nb_si_v1 > nb_si_v0) {
          update_v0 = true;
        }
        else {
          update_v0 = true;
          update_v1 = true;
        }
      }
      return;
    }

    if (dim0 == 2 && dim1 == 2)
    {
      if (is_boundary(c3t3, edge, cell_selector))
      {
        if (!topology_test(edge, c3t3, cell_selector))
          return;
        Subdomain_relation subdomain_rel = compare_subdomains(v0, v1, c3t3);

        //Vertices on the same surface
        if (subdomain_rel == INCLUDES) {
          update_v1 = true;
        }
        else if (subdomain_rel == INCLUDED) {
          update_v0 = true;
        }
        else if (subdomain_rel == EQUAL)
        {
          if (c3t3.number_of_edges() == 0)
          {
            update_v0 = true;
            update_v1 = true;
          }
          else
          {
            bool v0_on_feature = is_on_feature(v0);
            bool v1_on_feature = is_on_feature(v1);

            if (v0_on_feature && v1_on_feature) {
              if (c3t3.is_in_complex(edge)) {
                if (!c3t3.is_in_complex(v0))
                  update_v0 = true;
                if (!c3t3.is_in_complex(v1))
                  update_v1 = true;
              }
            }
            else {
              if (!v0_on_feature) {
                update_v0 = true;
              }
              if (!v1_on_feature) {
                update_v1 = true;
              }
            }
          }
        }
      }

      return;
    }
    //In the case of mixte edges
    if (dim0 == 2 && dim1 == 3 && !is_v1_on_hull) {
      update_v1 = true;
      return;
    }

    if (dim1 == 2 && dim0 == 3 && !is_v0_on_hull) {
      update_v0 = true;
      return;
    }
  }

  template<typename C3t3, typename OutputIterator>
  OutputIterator incident_subdomains(const typename C3t3::Vertex_handle v,
    const C3t3& c3t3,
    OutputIterator oit)
  {
    typedef typename C3t3::Triangulation::Cell_handle Cell_handle;
    std::vector<Cell_handle> cells;
    c3t3.triangulation().incident_cells(v, std::back_inserter(cells));

    for (std::size_t i = 0; i < cells.size(); ++i)
      *oit++ = cells[i]->subdomain_index();

    return oit;
  }

  template<typename C3t3, typename OutputIterator>
  OutputIterator incident_subdomains(const typename C3t3::Edge& e,
    const C3t3& c3t3,
    OutputIterator oit)
  {
    typedef typename C3t3::Triangulation::Cell_circulator Cell_circulator;

    Cell_circulator circ = c3t3.triangulation().incident_cells(e);
    Cell_circulator end = circ;
    do
    {
      *oit++ = circ->subdomain_index();
    } while (++circ != end);

    return oit;
  }

  template<typename C3t3>
  std::size_t nb_incident_subdomains(const typename C3t3::Vertex_handle v,
    const C3t3& c3t3)
  {
    typedef typename C3t3::Subdomain_index Subdomain_index;

    boost::unordered_set<Subdomain_index> indices;
    incident_subdomains(v, c3t3, std::inserter(indices, indices.begin()));

    return indices.size();
  }

  template<typename C3t3>
  std::size_t nb_incident_subdomains(const typename C3t3::Edge& e,
    const C3t3& c3t3)
  {
    typedef typename C3t3::Subdomain_index Subdomain_index;

    boost::unordered_set<Subdomain_index> indices;
    incident_subdomains(e, c3t3, std::inserter(indices, indices.begin()));

    return indices.size();
  }

  template<typename C3t3>
  std::size_t nb_incident_complex_edges(const typename C3t3::Vertex_handle v,
    const C3t3& c3t3)
  {
    typedef typename C3t3::Edge Edge;
    boost::unordered_set<Edge> edges;
    c3t3.triangulation().incident_edges(v,
      std::inserter(edges, edges.begin()));

    std::size_t count = 0;
    for (typename boost::unordered_set<Edge>::iterator eit = edges.begin();
      eit != edges.end();
      ++eit)
    {
      if (c3t3.is_in_complex(*eit))
        ++count;
    }
    return count;
  }


  template<typename C3t3>
  bool is_feature(const typename C3t3::Vertex_handle v,
    const typename C3t3::Vertex_handle neighbor,
    const C3t3& c3t3)
  {
    typename C3t3::Cell_handle ch;
    int i0, i1;
    if (c3t3.triangulation().is_edge(v, neighbor, ch, i0, i1))
    {
      typename C3t3::Edge edge(ch, i0, i1);
      return c3t3.is_in_complex(edge);
    }
    return false;
  }

  template<typename C3t3>
  bool is_feature(const typename C3t3::Vertex_handle v, const C3t3& c3t3)
  {
    typedef typename C3t3::Edge Edge;

    if (nb_incident_subdomains(v, c3t3) > 2)
    {
      std::vector<Edge> edges;
      c3t3.triangulation().finite_incident_edges(v, std::back_inserter(edges));

      int feature_count = 0;
      BOOST_FOREACH(Edge ei, edges)
      {
        if (c3t3.is_in_complex(ei))
        {
          feature_count++;
          if (feature_count >= 3)
            return true;
        }
      }
    }
    else if (c3t3.number_of_corners() > 0)
    {
      return c3t3.is_in_complex(v);
    }
    return false;
  }

  /**
  * returns true iff `v` is on the outer hull of c3t3.triangulation()
  * i.e. finite and incident to at least one infinite cell
  */
  template<typename C3t3>
  bool is_on_hull(const typename C3t3::Vertex_handle v,
                  const C3t3& c3t3)
  {
    if (v == c3t3.triangulation().infinite_vertex())
      return true;

    //on hull == incident to infinite cell
    typedef typename C3t3::Triangulation::Cell_handle Cell_handle;

    std::vector<Cell_handle> cells;
    c3t3.triangulation().incident_cells(v, std::back_inserter(cells));
    for (std::size_t i = 0; i < cells.size(); ++i)
    {
      if (c3t3.triangulation().is_infinite(cells[i]))
        return true;
    }
    return false;
  }

  template<typename C3t3>
  bool is_on_domain_hull(const typename C3t3::Vertex_handle v,
    const C3t3& c3t3,
    const typename C3t3::Subdomain_index& imaginary_index)
  {
    if (v == c3t3.triangulation().infinite_vertex())
      return false;

    on hull == incident to infinite cell
    typedef typename C3t3::Triangulation::Cell_handle Cell_handle;

    bool met_inside_cell = false;
    bool met_outside_cell = false;

    std::vector<Cell_handle> cells;
    c3t3.triangulation().incident_cells(v, std::back_inserter(cells));
    for (std::size_t i = 0; i < cells.size(); ++i)
    {
      if (c3t3.triangulation().is_infinite(cells[i])
        || !c3t3.is_in_complex(cells[i])
        || cells[i]->subdomain_index() == imaginary_index)
        met_outside_cell = true;
      else
        met_inside_cell = true;

      if (met_inside_cell && met_outside_cell)
        return true;
    }
    return false;
  }

  /**
  * returns true iff `edge` is on the outer hull
  * of c3t3.triangulation()
  * i.e. finite and incident to at least one infinite cell
  */
  template<typename C3t3>
  bool is_on_hull(const typename C3t3::Edge & edge,
    const C3t3& c3t3)
  {
    typedef typename C3t3::Triangulation::Cell_circulator Cell_circulator;
    Cell_circulator circ = c3t3.triangulation().incident_cells(edge);
    Cell_circulator done = circ;
    do
    {
      if (c3t3.triangulation().is_infinite(circ))
        return true;
    } while (++circ != done);

    return false;
  }

  template<typename C3t3>
  bool is_on_domain_hull(const typename C3t3::Edge & edge,
    const C3t3& c3t3,
    const typename C3t3::Subdomain_index& imaginary_index)
  {
    typedef typename C3t3::Triangulation::Cell_circulator Cell_circulator;

    bool met_inside_cell = false;
    bool met_outside_cell = false;

    Cell_circulator circ = c3t3.triangulation().incident_cells(edge);
    Cell_circulator done = circ;
    do
    {
      if (c3t3.triangulation().is_infinite(circ)
        || !c3t3.is_in_complex(circ)
        || circ->subdomain_index() == imaginary_index)
        met_outside_cell = true;
      else
        met_inside_cell = true;

      if (met_inside_cell && met_outside_cell)
        return true;
    } while (++circ != done);

    return false;
  }

  template<typename C3t3>
  bool is_imaginary(const typename C3t3::Vertex_handle v,
    const C3t3& c3t3,
    const typename C3t3::Subdomain_index& imaginary_index)
  {
    typedef typename C3t3::Triangulation::Cell_handle Cell_handle;

    std::vector<Cell_handle> cells;
    c3t3.triangulation().incident_cells(v, std::back_inserter(cells));

    BOOST_FOREACH(Cell_handle c, cells)
    {
      if (c->subdomain_index() != imaginary_index)
        return false;
    }
    return true;
  }

  /**
  * returns true off edge is fully imaginary
  * i.e. if all its incident cells are not in the complex,
  * and have their subdomain index == imaginary_index
  */
  template<typename C3t3>
  bool is_imaginary(const typename C3t3::Edge & edge,
    const C3t3& c3t3,
    const typename C3t3::Subdomain_index& imaginary_index)
  {
    typedef typename C3t3::Triangulation::Cell_circulator Cell_circulator;
    Cell_circulator circ = c3t3.triangulation().incident_cells(edge);
    Cell_circulator done = circ;
    do
    {
      if (c3t3.is_in_complex(circ)
        && circ->subdomain_index() != imaginary_index)
        return false;
    } while (++circ != done);

    return    true;
  }

  template<typename C3t3, typename CellSelector>
  bool is_outside(const typename C3t3::Edge & edge,
    const C3t3& c3t3,
    const typename C3t3::Subdomain_index& imaginary_index,
    CellSelector cell_selector)
  {
    typedef typename C3t3::Triangulation::Cell_circulator Cell_circulator;
    Cell_circulator circ = c3t3.triangulation().incident_cells(edge);
    Cell_circulator done = circ;
    do
    {
      // is cell infinite?
      if (c3t3.triangulation().is_infinite(circ))
        continue;
      // is cell imaginary?
      if (c3t3.is_in_complex(circ) && circ->subdomain_index() == imaginary_index)
        continue;
      // circ does not belong to the selection
      if (!cell_selector(circ))
        continue;

      // none of the above conditions was met
      return false;
    } while (circ != done);

    return true; //all cells have met the loop conditions
  }

  template<typename C3t3, typename CellSelector>
  bool is_selected(const typename C3t3::Vertex_handle v,
    const C3t3& c3t3,
    CellSelector cell_selector)
  {
    typedef typename C3t3::Triangulation::Cell_handle Cell_handle;

    std::vector<Cell_handle> cells;
    c3t3.triangulation().incident_cells(v, std::back_inserter(cells));

    BOOST_FOREACH(Cell_handle c, cells)
    {
      if (!cell_selector(c))
        return false;
    }
    return true;
  }

  template<typename C3t3, typename CellSelector>
  bool is_inside(const typename C3t3::Edge& edge,
    const C3t3& c3t3,
    const typename C3t3::Subdomain_index& imaginary_index,
    CellSelector cell_selector)
  {
    typedef typename C3t3::Triangulation::Cell_circulator Cell_circulator;
    Cell_circulator circ = c3t3.triangulation().incident_cells(edge);
    Cell_circulator done = circ;

    const typename C3t3::Subdomain_index si = circ->subdomain_index();
    if (si == imaginary_index || !c3t3.is_in_complex(circ))
      return false;
    do
    {
      if (c3t3.triangulation().is_infinite(circ))
        return false;
      if (si != circ->subdomain_index())
        return false;
      if (!cell_selector(circ))
        return false;
    } while (++circ != done);

    return true;
  }

  template<typename Tr, typename K>
  bool is_convex(const Tr& tr,
    const CGAL::Iso_cuboid_3<K>& bbox,
    typename Tr::Facet& facet)
  {
    typedef typename Tr::Cell_handle            Cell_handle;
    typedef typename Tr::Vertex_handle          Vertex_handle;
    typedef typename Tr::Facet                  Facet;
    typedef typename Tr::Finite_facets_iterator Finite_facets_iterator;

    for (Finite_facets_iterator fit = tr.finite_facets_begin();
      fit != tr.finite_facets_end(); ++fit)
    {
      Facet f = *fit;
      Facet mf = tr.mirror_facet(f);
      if (!tr.is_infinite(f.first) && !tr.is_infinite(mf.first))
        continue;

      if (tr.is_infinite(mf.first))
        f = mf;
      CGAL_assertion(tr.is_infinite(f.first));

      boost::array<Vertex_handle, 3> vs;
      for (int i = 0; i < 3; ++i)
        vs[i] = f.first->vertex((f.second + i + 1) % 4);
      if (f.second % 2 == 0)
        std::swap(vs[0], vs[1]);

      Cell_handle fin_c = f.first->neighbor(f.second);
      Vertex_handle v4 = fin_c->vertex(fin_c->index(f.first));

      CGAL_assertion(!tr.is_infinite(fin_c));
      CGAL_assertion(!f.first->has_vertex(v4));

      CGAL_assertion(CGAL::NEGATIVE
        == CGAL::orientation(vs[0]->point(), vs[1]->point(),
          vs[2]->point(), v4->point()));

      for (int i = 1; i < 4; ++i)
      {
        nfi is neighbor of f on convex hull
        Cell_handle ni = f.first->neighbor((f.second + i) % 4);
        CGAL_assertion(tr.is_infinite(ni));

        collect points
        Vertex_handle v3 = ni->vertex(ni->index(f.first));
        CGAL_assertion(v3 != vs[0] && v3 != vs[1] && v3 != vs[2]
          && v3 != tr.infinite_vertex());
        CGAL_assertion(!f.first->has_vertex(v3));

        CGAL::Orientation o2 = CGAL::orientation(vs[0]->point(),
          vs[1]->point(), vs[2]->point(), v3->point());
        if (o2 == CGAL::POSITIVE)
        {
          facet = f;

          if (!bbox.is_degenerate()
            && bbox.has_on_boundary(vs[0]->point())
            && bbox.has_on_boundary(vs[1]->point())
            && bbox.has_on_boundary(vs[2]->point()))
          {
            facet = Facet(ni, ni->index(tr.infinite_vertex()));
          }

          return false;
        }
      }
    }
    return true;
  }

  template<typename Tr>
  bool is_convex(const Tr& tr)
  {
    typename Tr::Facet f;
    typename Tr::Geom_traits::Iso_cuboid_3 bb(CGAL::ORIGIN, CGAL::ORIGIN);
    return is_convex(tr, bb, f);
  }

  template<typename Facet, typename Gt>
  typename Gt::Vector_3 normal(const Facet& f, const Gt& gt)
  {
    namespace PMP = CGAL::Polygon_mesh_processing;
    typedef typename Gt::Vector_3 Vector;
    typedef typename Gt::Point_3  Point;

    Point p0 = point(f.first->vertex((f.second + 1) % 4)->point());
    Point p1 = point(f.first->vertex((f.second + 2) % 4)->point());
    const Point& p2 = point(f.first->vertex((f.second + 3) % 4)->point());

    if (f.second % 2 == 0)//equivalent to the commented orientation test
      std::swap(p0, p1);

    Vector n = PMP::internal::triangle_normal(p0, p1, p2, gt);

    if (!typename Gt::Equal_3()(n, CGAL::NULL_VECTOR))
      PMP::internal::normalize(n, gt);

    return n;
  }

  template<typename C3t3, typename CellSelector, typename OutputIterator>
  OutputIterator get_inside_edges(const C3t3& c3t3,
    const typename C3t3::Subdomain_index& imaginary_index,
    CellSelector cell_selector,
    OutputIterator oit)/*holds pairs of Vertex_handles*/
  {
    for (typename C3t3::Triangulation::Finite_edges_iterator
      eit = c3t3.triangulation().finite_edges_begin();
      eit != c3t3.triangulation().finite_edges_end();
      ++eit)
    {
      const typename C3t3::Edge& e = *eit;
//        if ( !c3t3.is_in_complex(e)
//          && !is_boundary_edge(e, c3t3)
//          && !is_on_hull(e, c3t3)
//          && !is_imaginary(e, c3t3, imaginary_index))
      if (is_inside(e, c3t3, imaginary_index, cell_selector))
      {
        *oit++ = make_vertex_pair<typename C3t3::Triangulation>(e);
      }
    }
    return oit;
  }


  namespace debug
  {
    // forward-declaration
    template<typename Tr, typename CellRange>
    void dump_cells(const CellRange& cells, const char* filename);

    template <typename Bimap>
    void dump_edges(const Bimap& edges, const char* filename)
    {
      std::ofstream ofs(filename);
      ofs.precision(17);

      BOOST_FOREACH(typename Bimap::left_const_reference it, edges.left)
      {
        ofs << "2 " << it.first.first->point()
          << " " << it.first.second->point() << std::endl;
      }
      ofs.close();
    }

    template<typename Facet, typename OutputStream>
    void dump_facet(const Facet& f, OutputStream& os)
    {
      os << "4 ";
      os << f.first->vertex((f.second + 1) % 4)->point() << " "
        << f.first->vertex((f.second + 2) % 4)->point() << " "
        << f.first->vertex((f.second + 3) % 4)->point() << " "
        << f.first->vertex((f.second + 1) % 4)->point();
      os << std::endl;
    }

    template<typename FacetRange>
    void dump_facets(const FacetRange& facets, const char* filename)
    {
      std::ofstream os(filename);
      for (typename FacetRange::const_iterator fit = facets.begin();
        fit != facets.end(); ++fit)
      {
        typename FacetRange::value_type f = *fit;
        dump_facet(f, os);
      }
    }

    template<typename CellRange>
    void dump_polylines(const CellRange& cells, const char* filename)
    {
      std::ofstream ofs(filename);
      if (!ofs) return;

      for (typename CellRange::const_iterator it = cells.begin();
        it != cells.end(); ++it)
      {
        for (int i = 0; i < 4; ++i)
          dump_facet(std::make_pair(*it, i), ofs);
      }
      ofs.close();
    }

    template<typename Tr>
    void dump_surface_off(const Tr& tr, const char* filename)
    {
      typedef typename Tr::Vertex_handle              Vertex_handle;
      typedef typename Tr::Cell_handle                Cell_handle;
      typedef typename Tr::Finite_facets_iterator     Finite_facets_iterator;
      typedef boost::bimap<Vertex_handle, int>                   Bimap_t;
      typedef typename Bimap_t::left_map::value_type             value_type;

      //collect vertices
      Bimap_t vertices;
      std::size_t nbf = 0;
      int index = 0;
      for (Finite_facets_iterator fit = tr.finite_facets_begin();
        fit != tr.finite_facets_end(); ++fit)
      {
        Cell_handle c = fit->first;
        int i = fit->second;
        if (tr.is_infinite(c) || tr.is_infinite(c->neighbor(i)))
        {
          nbf++;
          for (int j = 1; j < 4; ++j)
          {
            Vertex_handle vij = c->vertex((i + j) % 4);
            if (vertices.left.find(vij) == vertices.left.end())
              vertices.left.insert(value_type(vij, index++));
          }
        }
      }

      //write header
      std::ofstream ofs(filename);
      ofs.precision(17);
      ofs << "OFF" << std::endl;
      ofs << vertices.left.size() << " " << nbf << " 0" << std::endl << std::endl;

      // write vertices
      for (typename Bimap_t::right_iterator vit = vertices.right.begin();
        vit != vertices.right.end(); ++vit)
      {
        ofs << vit->second->point() << std::endl;
      }

      //write facets
      std::size_t nbf_print = 0;
      for (Finite_facets_iterator fit = tr.finite_facets_begin();
        fit != tr.finite_facets_end(); ++fit)
      {
        Cell_handle c = fit->first;
        int i = fit->second;
        if (tr.is_infinite(c) || tr.is_infinite(c->neighbor(i)))
        {
          ofs << "3  " << vertices.left.at(c->vertex((i + 1) % 4)) << " "
            << vertices.left.at(c->vertex((i + 2) % 4)) << " "
            << vertices.left.at(c->vertex((i + 3) % 4)) << std::endl;
          ++nbf_print;
        }
      }
      CGAL_assertion(nbf == nbf_print);

      ofs.close();
    }

    template<typename Tr>
    void dump_cells_off(const Tr& tr, const char* filename)
    {
      typedef typename Tr::Vertex_handle              Vertex_handle;
      typedef typename Tr::Cell_handle                Cell_handle;
      typedef typename Tr::Finite_facets_iterator     Finite_facets_iterator;
      typedef typename Tr::Finite_vertices_iterator   Finite_vertices_iterator;
      typedef boost::bimap<Vertex_handle, int>                   Bimap_t;
      typedef typename Bimap_t::left_map::value_type             value_type;

      //write header
      std::ofstream ofs(filename);
      ofs.precision(17);
      ofs << "OFF" << std::endl;
      ofs << tr.number_of_vertices()
        << " " << tr.number_of_finite_facets() << " 0" << std::endl << std::endl;

      //collect and write vertices
      Bimap_t vertices;
      int index = 0;
      for (Finite_vertices_iterator vit = tr.finite_vertices_begin();
        vit != tr.finite_vertices_end(); ++vit)
      {
        vertices.left.insert(value_type(vit, index++));
        ofs << vit->point().x() << " "
          << vit->point().y() << " "
          << vit->point().z() << std::endl;
      }

      //write facets
      for (Finite_facets_iterator fit = tr.finite_facets_begin();
        fit != tr.finite_facets_end(); ++fit)
      {
        Cell_handle c = fit->first;
        int i = fit->second;
        ofs << "3  " << vertices.left.at(c->vertex((i + 1) % 4)) << " "
          << vertices.left.at(c->vertex((i + 2) % 4)) << " "
          << vertices.left.at(c->vertex((i + 3) % 4)) << std::endl;
      }
      ofs.close();
    }

    template<typename Tr, typename CellRange, typename IndexRange>
    void dump_cells(const CellRange& cells,
      const IndexRange& indices,
      const char* filename)
    {
      typedef typename Tr::Vertex_handle                         Vertex_handle;
      typedef typename Tr::Point                                 Point;
      typedef boost::bimap<Vertex_handle, int>                   Bimap_t;
      typedef typename Bimap_t::left_map::value_type             value_type;

      CGAL_assertion(indices.empty() || cells.size() == indices.size());

      //collect vertices
      Bimap_t vertices;
      int index = 1;
      for (typename CellRange::const_iterator cit = cells.begin();
        cit != cells.end();
        ++cit)
      {
        for (int i = 0; i < 4; ++i)
        {
          Vertex_handle vi = (*cit)->vertex(i);
          if (vertices.left.find(vi) == vertices.left.end())
            vertices.left.insert(value_type(vi, index++));
        }
      }

      //write cells
      std::ofstream ofs(filename);
      ofs.precision(17);
      ofs << "MeshVersionFormatted 1" << std::endl;
      ofs << "Dimension 3" << std::endl;
      ofs << "Vertices" << std::endl << vertices.size() << std::endl;
      for (typename Bimap_t::right_const_iterator vit = vertices.right.begin();
        vit != vertices.right.end();
        ++vit)
      {
        const Point& p = vit->second->point();
        ofs << p.x() << " " << p.y() << " " << p.z() << " 2" << std::endl;
      }
      ofs << "Tetrahedra " << std::endl << cells.size() << std::endl;
      typename IndexRange::const_iterator iit = indices.begin();
      for (typename CellRange::const_iterator cit = cells.begin();
        cit != cells.end();
        ++cit)
      {
        ofs << vertices.left.at((*cit)->vertex(0))
          << " " << vertices.left.at((*cit)->vertex(1))
          << " " << vertices.left.at((*cit)->vertex(2))
          << " " << vertices.left.at((*cit)->vertex(3));

        if (iit == indices.end())
          ofs << " 1" << std::endl;
        else
        {
          // std::cerr << "Cell #" << (cit - cells.begin())
          //           << " has original index " << *iit << std::endl;
          ofs << " " << (*iit) << std::endl;
          ++iit;
        }
      }
      ofs << "End" << std::endl;
      ofs.close();
    }

    template<typename Tr, typename CellRange>
    void dump_cells(const CellRange& cells, const char* filename)
    {
      std::vector<int> indices;
      dump_cells<Tr>(cells, indices, filename);
    }

    template<typename Tr>
    void dump_cells_in_complex(const Tr& tr, const char* filename)
    {
      std::vector<typename Tr::Cell_handle> cells;
      std::vector<typename Tr::Cell::Subdomain_index> indices;

      for (typename Tr::Finite_cells_iterator cit = tr.finite_cells_begin();
        cit != tr.finite_cells_end(); ++cit)
      {
        if (cit->subdomain_index() > 0)
        {
          cells.push_back(cit);
          indices.push_back(cit->subdomain_index());
        }
      }
      dump_cells<Tr>(cells, indices, filename);
    }

    template<typename C3t3>
    void dump_facets_in_complex(const C3t3& c3t3, const char* filename)
    {
      typedef typename C3t3::Triangulation              Tr;
      typedef typename Tr::Vertex_handle                Vertex_handle;
      typedef typename Tr::Cell_handle                  Cell_handle;
      typedef typename C3t3::Facets_in_complex_iterator Facets_in_complex_iterator;
      typedef boost::bimap<Vertex_handle, int>          Bimap_t;
      typedef typename Bimap_t::left_map::value_type  value_type;

      //collect vertices
      Bimap_t vertices;
      std::size_t nbf = 0;
      int index = 0;
      for (Facets_in_complex_iterator fit = c3t3.facets_in_complex_begin();
        fit != c3t3.facets_in_complex_end(); ++fit)
      {
        Cell_handle c = fit->first;
        int i = fit->second;

        nbf++;
        for (int j = 1; j < 4; ++j)
        {
          Vertex_handle vij = c->vertex((i + j) % 4);
          if (vertices.left.find(vij) == vertices.left.end())
            vertices.left.insert(value_type(vij, index++));
        }
      }

      //write header
      std::ofstream ofs(filename);
      ofs.precision(17);
      ofs << "OFF" << std::endl;
      ofs << vertices.left.size() << " " << nbf << " 0" << std::endl << std::endl;

      // write vertices
      for (typename Bimap_t::right_iterator vit = vertices.right.begin();
        vit != vertices.right.end(); ++vit)
      {
        ofs << vit->second->point() << std::endl;
      }

      //write facets
      std::size_t nbf_print = 0;
      for (Facets_in_complex_iterator fit = c3t3.facets_in_complex_begin();
        fit != c3t3.facets_in_complex_end(); ++fit)
      {
        Cell_handle c = fit->first;
        int i = fit->second;
        ofs << "3  " << vertices.left.at(c->vertex((i + 1) % 4)) << " "
          << vertices.left.at(c->vertex((i + 2) % 4)) << " "
          << vertices.left.at(c->vertex((i + 3) % 4)) << std::endl;
        ++nbf_print;
      }
      CGAL_assertion(nbf == nbf_print);

      ofs.close();
    }

    template<typename C3T3>
    void dump_edges_in_complex(const C3T3& c3t3, const char* filename)
    {
      std::ofstream ofs(filename);
      ofs.precision(17);
      for (typename C3T3::Edges_in_complex_iterator eit = c3t3.edges_in_complex_begin();
        eit != c3t3.edges_in_complex_end(); ++eit)
      {
        const typename C3T3::Edge& e = *eit;
        ofs << "2 "
          << e.first->vertex(e.second)->point() << " "
          << e.first->vertex(e.third)->point() << "\n";
      }
      ofs.close();
    }

    template<typename Tr>
    void dump_vertices_by_dimension(const Tr& tr, const char* prefix)
    {
      typedef typename Tr::Vertex_handle Vertex_handle;
      std::vector< std::vector<Vertex_handle> > vertices_per_dimension(4);

      for (typename Tr::Finite_vertices_iterator
        vit = tr.finite_vertices_begin();
        vit != tr.finite_vertices_end();
        ++vit)
      {
        //vertices_per_dimension[vit->info()].push_back(vit);
        vertices_per_dimension[vit->in_dimension()].push_back(vit);
      }

      for (int i = 0; i < 4; ++i)
      {
        //dimension is i
        const std::vector<Vertex_handle>& vertices_di = vertices_per_dimension[i];

        std::cout << "Dimension " << i << " : " << vertices_di.size() << std::endl;

        std::ostringstream oss;
        oss << prefix << "_dimension_" << i << ".off";

        std::ofstream ofs(oss.str());
        ofs.precision(17);
        ofs << "OFF" << std::endl;
        ofs << vertices_di.size() << " 0 0" << std::endl << std::endl;

        for (std::size_t j = 0; j < vertices_di.size(); ++j)
        {
          ofs << vertices_di[j]->point() << std::endl;
        }

        ofs.close();
      }
    }

    template<typename Tr>
    void dump_triangulation_cells(const Tr& tr, const char* filename)
    {
      std::vector<typename Tr::Cell_handle> cells(tr.number_of_finite_cells());
      int i = 0;
      for (typename Tr::Finite_cells_iterator cit = tr.finite_cells_begin();
        cit != tr.finite_cells_end(); ++cit)
      {
        cells[i++] = cit;
      }
      dump_cells<Tr>(cells, filename);
    }

    template<typename Tr>
    void dump_without_imaginary(const Tr& tr, const char* filename,
      const int imaginary_index)
    {
      std::vector<typename Tr::Cell_handle> cells;
      std::vector<std::ptrdiff_t> indices;

      for (typename Tr::Finite_cells_iterator cit = tr.finite_cells_begin();
        cit != tr.finite_cells_end(); ++cit)
      {
        if (cit->subdomain_index() > 0
          && cit->subdomain_index() != imaginary_index)
        {
          cells.push_back(cit);
          indices.push_back(1);
        }
      }
      dump_cells<Tr>(cells, indices, filename);
    }

    //template<typename VertexPairsSet>
    //void dump_edges(const VertexPairsSet& edges, const char* filename)
    //{
    //  std::ofstream ofs(filename);
    //  BOOST_FOREACH(typename VertexPairsSet::key_type vp, edges)
    //  {
    //    ofs << "2 " << vp.first->point()
    //      << " " << vp.second->point() << std::endl;
    //  }
    //  ofs.close();
    //}

  } //namespace debug
 } //namespace Tetrahedral_remeshing
} //namespace CGAL

#endif //CGAL_INTERNAL_TET_REMESHING_HELPERS_H
