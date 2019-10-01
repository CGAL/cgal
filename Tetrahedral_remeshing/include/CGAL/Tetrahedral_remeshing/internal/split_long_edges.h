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

#ifndef CGAL_INTERNAL_SPLIT_LONG_EDGES_H
#define CGAL_INTERNAL_SPLIT_LONG_EDGES_H

#include <boost/bimap.hpp>
#include <boost/bimap/set_of.hpp>
#include <boost/bimap/multiset_of.hpp>
#include <boost/unordered_map.hpp>

#include <CGAL/Tetrahedral_remeshing/internal/tetrahedral_remeshing_helpers.h>

#include <functional>
#include <utility>

namespace CGAL
{
namespace Tetrahedral_remeshing
{
namespace internal
{
  template<typename C3t3>
  typename C3t3::Vertex_handle split_edge(const typename C3t3::Edge& e,
                                          C3t3& c3t3)
  {
    typedef typename C3t3::Triangulation   Tr;
    typedef typename C3t3::Subdomain_index Subdomain_index;
    typedef typename Tr::Point             Point;
    typedef typename Tr::Facet             Facet;
    typedef typename Tr::Vertex_handle     Vertex_handle;
    typedef typename Tr::Cell_handle       Cell_handle;
    typedef typename Tr::Cell_circulator   Cell_circulator;
    typedef typename Tr::Cell::Info        Cell_info;

    Tr& tr = c3t3.triangulation();
    Vertex_handle v1 = e.first->vertex(e.second);
    Vertex_handle v2 = e.first->vertex(e.third);

    //backup subdomain info of incident cells before making changes
    short dimension = (c3t3.is_in_complex(e)) ? 1 : 3;
    boost::unordered_map<Facet, std::pair<Subdomain_index, Cell_info> > info;

    Cell_circulator circ = tr.incident_cells(e);
    Cell_circulator end = circ;
    Subdomain_index prev = c3t3.subdomain_index(circ);
    Subdomain_index curr = prev;
    do
    {
      //keys are the opposite facets to the ones not containing e,
      //because they will not be modified
      Facet opp_facet = tr.mirror_facet(Facet(circ, circ->index(v1)));
      info.insert(std::make_pair(opp_facet,
          std::make_pair(c3t3.subdomain_index(circ), circ->info())));

      opp_facet = tr.mirror_facet(Facet(circ, circ->index(v2)));
      info.insert(std::make_pair(opp_facet,
          std::make_pair(c3t3.subdomain_index(circ), circ->info())));

      ++circ;
      prev = curr;
      curr = c3t3.subdomain_index(circ);

      if (prev != curr && dimension == 3)
        dimension = 2;
    } while (circ != end);

    // insert midpoint
    Vertex_handle new_v = tr.tds().insert_in_edge(e);
    const Point m(CGAL::midpoint(point(v1->point()), point(v2->point())));
    new_v->set_point(m);

    // update dimension
    c3t3.set_dimension(new_v, dimension);

    // update c3t3
    std::vector<Cell_handle> new_cells;
    tr.incident_cells(new_v, std::back_inserter(new_cells));
    for (std::size_t i = 0; i < new_cells.size(); ++i)
    {
      Cell_handle nci = new_cells[i];
      Facet fi(nci, nci->index(new_v));
      Subdomain_index n_index = info.at(tr.mirror_facet(fi)).first;
      c3t3.set_subdomain_index(nci, n_index);
      nci->info() = info.at(tr.mirror_facet(fi)).second;
    }

    return new_v;
  }

  template<typename C3T3, typename CellSelector>
  bool can_be_split(const typename C3T3::Edge& e,
                    const C3T3& c3t3,
                    const bool protect_boundaries,
                    const typename C3T3::Subdomain_index& imaginary_index,
                    CellSelector cell_selector)
  {
    if (is_outside(e, c3t3, imaginary_index, cell_selector))
      return false;
    if (is_imaginary(e, c3t3, imaginary_index))
      return false;

#ifdef CGAL_LIMITED_APERTURE_EDGE_SELECTION
    if (CGAL::helpers::is_on_the_outer_box(e, c3t3, imaginary_index))
      return true;
#endif

    if (protect_boundaries)
    {
      if (c3t3.is_in_complex(e))
        return false;
      else if (helpers::is_boundary(c3t3, e, cell_selector))
        return false;

#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
      if (!is_inside(e, c3t3, imaginary_index, cell_selector))
      {
        std::cerr << "e is not inside!?" << std::endl;
        typename C3T3::Vertex_handle v1 = e.first->vertex(e.second);
        typename C3T3::Vertex_handle v2 = e.first->vertex(e.third);
        std::cerr << v1->point() << " " << v2->point() << std::endl;
      }
#endif

      CGAL_assertion(is_inside(e, c3t3, imaginary_index, cell_selector));
      return true;
    }
    else
    {
      return true;
    }
  }

  template<typename C3T3, typename CellSelector>
  void split_long_edges(C3T3& c3t3,
    const typename C3T3::Triangulation::Geom_traits::FT& high,
    const bool protect_boundaries,
    const typename C3T3::Subdomain_index& imaginary_index,
    CellSelector cell_selector)
  {
    typedef typename C3T3::Triangulation       T3;
    typedef typename T3::Cell_handle           Cell_handle;
    typedef typename T3::Edge                  Edge;
    typedef typename T3::Finite_edges_iterator Finite_edges_iterator;
    typedef typename T3::Vertex_handle         Vertex_handle;
    typedef typename std::pair<Vertex_handle, Vertex_handle> Edge_vv;

    typedef typename T3::Geom_traits::FT FT;
    typedef boost::bimap<
      boost::bimaps::set_of<Edge_vv>,
      boost::bimaps::multiset_of<FT, std::greater<FT> > >  Boost_bimap;
    typedef typename Boost_bimap::value_type               long_edge;

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
    std::cout << "Split long edges (" << high << ")...";
    std::cout.flush();
    std::size_t nb_splits = 0;
#endif
    const FT sq_high = high*high;

    //collect long edges
    T3& tr = c3t3.triangulation();
    Boost_bimap long_edges;
    for (Finite_edges_iterator eit = tr.finite_edges_begin();
         eit != tr.finite_edges_end(); ++eit)
    {
      Edge e = *eit;
      if (!can_be_split(e, c3t3, protect_boundaries, imaginary_index, cell_selector))
        continue;

      FT sqlen = tr.segment(e).squared_length();
      if (sqlen > sq_high)
        long_edges.insert(long_edge(make_vertex_pair<T3>(e), sqlen));
    }

#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
    debug::dump_edges(long_edges, "long_edges.polylines.txt");

    std::ofstream ofs("midpoints.off");
    ofs << "OFF" << std::endl;
    ofs << long_edges.size() << " 0 0" << std::endl;
#endif
    while(!long_edges.empty())
    {
      //the edge with longest length
      typename Boost_bimap::right_map::iterator eit = long_edges.right.begin();
      Edge_vv e = eit->second;
#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE_PROGRESS
      const double sqlen = eit->first;
#endif
      long_edges.right.erase(eit);

      Cell_handle cell;
      int i1, i2;
      if ( tr.tds().is_edge(e.first, e.second, cell, i1, i2))
      {
        Edge edge(cell, i1, i2);

        //check that splittability has not changed
        if (!can_be_split(edge, c3t3, protect_boundaries, imaginary_index, cell_selector))
          continue;

        Vertex_handle vh = split_edge(edge, c3t3);
        //CGAL_assertion(tr.is_valid(true));

#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
        ofs << vh->point() << std::endl;
#endif
        if (vh != Vertex_handle())
        {
#if  defined(CGAL_TETRAHEDRAL_REMESHING_VERBOSE_PROGRESS) \
  || defined(CGAL_TETRAHEDRAL_REMESHING_VERBOSE)
          ++nb_splits;
#endif
          ////insert newly created edges if needed
          //std::vector<Edge> new_edges;
          //tr.incident_edges(vh, std::back_inserter(new_edges));
          //
          //for (std::size_t i = 0; i < new_edges.size(); ++i)
          //{
          //  const Edge& ei = new_edges[i];
          //  Segment seg(ei.first->vertex(ei.second)->point(),
          //              ei.first->vertex(ei.third)->point());
          //
          //  const FT sqlen_i = seg.squared_length();
          //  if (sqlen_i > sq_high)
          //    long_edges.insert(long_edge(make_vertex_pair<T3>(ei), sqlen_i));
          //}
        }

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE_PROGRESS
        std::cout << "\rSplit (" << high << ")... ("
          << long_edges.left.size() << " long edges, "
          << "length  = " << std::sqrt(sqlen) << ", "
          << std::sqrt(CGAL::squared_distance(point(e.first->point()),
                                              point(e.second->point()))) << ", "
          << nb_splits << " splits)";
        std::cout.flush();
#endif
      }
    }//end loop on long_edges

#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
    if(ofs.is_open())
      ofs.close();
#endif

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
    std::cout << " done (" << nb_splits << " splits)." << std::endl;
#endif
    }
}
}
}

#endif // CGAL_INTERNAL_SPLIT_LONG_EDGES_H
