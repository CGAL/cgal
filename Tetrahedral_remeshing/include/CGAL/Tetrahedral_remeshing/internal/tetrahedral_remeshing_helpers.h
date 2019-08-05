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

namespace CGAL
{
namespace Tetrahedral_remeshing
{
  enum Subdomain_relation { EQUAL, DIFFERENT, INCLUDED, INCLUDES };
  enum Sliver_removal_result { INVALID_ORIENTATION, INVALID_CELL, INVALID_VERTEX,
    NOT_FLIPPABLE, EDGE_PROBLEM, VALID_FLIP, NO_BEST_CONFIGURATION, EXISTING_EDGE };

namespace helpers
{
  template<typename SubdomainIndex>
  std::pair<SubdomainIndex, SubdomainIndex>
  make_surface_patch_index(const SubdomainIndex& s1, const SubdomainIndex& s2)
  {
    CGAL_assertion(s1 != s2);
    if (s1 < s2)  
      return std::make_pair(s1, s2);
    else
      return std::make_pair(s2, s1);
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
            if ( is_edge_in_complex(v0, vi, c3t3, cell_selector)
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

        if (nb_si_v0 > nb_si_v1){
          update_v1 = true;
        }
        else if (nb_si_v1 > nb_si_v0){
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
        if (subdomain_rel == INCLUDES){
          update_v1 = true;
        }
        else if (subdomain_rel == INCLUDED){
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

            if (v0_on_feature && v1_on_feature){
              if (c3t3.is_in_complex(edge)){
                if (!c3t3.is_in_complex(v0))
                  update_v0 = true;
                if (!c3t3.is_in_complex(v1))
                  update_v1 = true;
              }
            }
            else {
              if (!v0_on_feature){
                update_v0 = true;
              }
              if (!v1_on_feature){
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


  template<typename C3T3>
  void print_subdomain_indices(const C3T3& c3t3)
  {
    typedef typename C3T3::Triangulation       Tr;
    typedef typename Tr::Finite_cells_iterator Finite_cells_iterator;

    std::cout << "SUBDOMAINS : " << std::endl;
    unsigned int line_id = 0;
    for (Finite_cells_iterator cit = c3t3.triangulation().finite_cells_begin();
         cit != c3t3.triangulation().finite_cells_end();
         ++cit, ++line_id)
    {
      if (line_id % 10 == 0)
        std::cout << std::endl;
      std::cout << "\t" << cit->subdomain_index();
    }

  }

#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
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
#endif


}
}
}

#endif //CGAL_INTERNAL_TET_REMESHING_HELPERS_H
