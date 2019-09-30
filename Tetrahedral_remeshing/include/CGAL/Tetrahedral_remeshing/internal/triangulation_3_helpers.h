// Copyright (c) 2018 GeometryFactory (France).
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
// Author(s)     : Jane Tournois
//
//******************************************************************************
//
//******************************************************************************

#ifndef CGAL_TRIANGULATION_3_HELPERS_H
#define CGAL_TRIANGULATION_3_HELPERS_H

#include <CGAL/Iso_cuboid_3.h>
#include <CGAL/Point_3.h>
#include <CGAL/Weighted_point_3.h>
#include <CGAL/Vector_3.h>
#include <CGAL/utility.h>
#include <CGAL/Origin.h>

#include <CGAL/Polygon_mesh_processing/compute_normal.h>

#include <boost/bimap.hpp>
#include <boost/unordered_set.hpp>

#include <iostream>
#include <fstream>

#include <vector>

namespace CGAL
{
  namespace debug
  {
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

  } // end namespace debug (in ::CGAL)

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
    return Gt().compute_approximate_dihedral_angle_3_object()(p,q,r,s);
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
    return CGAL::POSITIVE == CGAL::orientation(point(v0->point()), point(v1->point()),
                                               point(v2->point()), point(v3->point()));
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
    }
    while (++circ != end);

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
    else if(c3t3.number_of_corners() > 0)
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

    //on hull == incident to infinite cell
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
      if ( c3t3.is_in_complex(circ)
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
      //circ does not belong to the selection
      if (!cell_selector(circ))
        continue;

      //none of the above conditions was met
      return false;
    }
    while (circ != done);

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
    if (si == imaginary_index || !c3t3.is_in_complex(circ) )
      return false;
    do
    {
      if (c3t3.triangulation().is_infinite(circ))
        return false;
      if (si != circ->subdomain_index())
        return false;
      if (!cell_selector(circ))
        return false;
    }
    while (++circ != done);

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
        //nfi is neighbor of f on convex hull
        Cell_handle ni = f.first->neighbor((f.second + i) % 4);
        CGAL_assertion(tr.is_infinite(ni));

        //collect points
        Vertex_handle v3 = ni->vertex(ni->index(f.first));
        CGAL_assertion(   v3 != vs[0] && v3 != vs[1] && v3 != vs[2]
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

    //if (CGAL::POSITIVE != CGAL::orientation(p0, p1, p2, p3))
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
      //if ( !c3t3.is_in_complex(e)
      //  && !is_boundary_edge(e, c3t3)
      //  && !is_on_hull(e, c3t3)
      //  && !is_imaginary(e, c3t3, imaginary_index))
      if (is_inside(e, c3t3, imaginary_index, cell_selector))
      {
        *oit++ = make_vertex_pair<typename C3t3::Triangulation>(e);
      }
    }
    return oit;
  }

namespace internal
{
  template<typename Tr, typename IncidentFacetsMap>
  void add_to_incidence_map(const typename Tr::Cell_handle& c,
                            const int& index,
                            IncidentFacetsMap& incidence_map,
                            const bool verbose = false)
  {
    CGAL_USE(verbose);
    typedef typename IncidentFacetsMap::key_type Vertex_set;
    Vertex_set vertices;
    vertices.insert(c->vertex((index + 1) % 4));
    vertices.insert(c->vertex((index + 2) % 4));
    vertices.insert(c->vertex((index + 3) % 4));
    CGAL_assertion(vertices.size() == 3);

    typename IncidentFacetsMap::iterator it = incidence_map.find(vertices);
    if (it == incidence_map.end())
    {
      std::vector<typename Tr::Facet> facets(1);
      facets[0] = typename Tr::Facet(c, index);
      incidence_map.insert(std::make_pair(vertices, facets));
    }
    else
    {
      it->second.push_back(typename Tr::Facet(c, index));

      CGAL_assertion(it->second.size() == 2);
      CGAL_assertion(it->second[0] != it->second[1]);
    }
  }

  template<typename Tr>
  typename Tr::Cell_handle
    create_neighbor_infinite_cell(const typename Tr::Cell_handle c,
                                  const int i,
                                  Tr& tr)
  {
    CGAL_assertion(!tr.is_infinite(c));
    CGAL_assertion_code(std::size_t nbc = tr.number_of_cells());

    typedef typename Tr::Cell_handle Cell_handle;
    Cell_handle opp_c;
    // the infinite cell that we are creating needs to be well oriented
    if (i == 0 || i == 2)
    {
      opp_c = create_cell(c->vertex((i + 3) % 4),
        tr.infinite_vertex(),
        c->vertex((i + 1) % 4),
        c->vertex((i + 2) % 4), tr);
    }
    else
    {
      opp_c = create_cell(tr.infinite_vertex(),
        c->vertex((i + 1) % 4),
        c->vertex((i + 2) % 4),
        c->vertex((i + 3) % 4), tr);
    }
    tr.infinite_vertex()->set_cell(opp_c);

    CGAL_assertion(nbc + 1 == tr.number_of_cells());
    return opp_c;
  }


}//end namespace internal

}//end namespace CGAL

#endif //CGAL_TRIANGULATION_3_HELPERS_H
