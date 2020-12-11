// Copyright (c) 2020 GeometryFactory (France) and Telecom Paris (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Jane Tournois, Noura Faraj, Jean-Marc Thiery, Tamy Boubekeur

#ifndef CGAL_INTERNAL_SPLIT_LONG_EDGES_H
#define CGAL_INTERNAL_SPLIT_LONG_EDGES_H

#include <CGAL/license/Tetrahedral_remeshing.h>

#include <boost/bimap.hpp>
#include <boost/bimap/set_of.hpp>
#include <boost/bimap/multiset_of.hpp>
#include <boost/unordered_map.hpp>
#include <boost/container/small_vector.hpp>

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
  typedef typename C3t3::Triangulation       Tr;
  typedef typename C3t3::Subdomain_index     Subdomain_index;
  typedef typename C3t3::Surface_patch_index Surface_patch_index;
  typedef typename Tr::Geom_traits::Point_3 Point;
  typedef typename Tr::Facet                Facet;
  typedef typename Tr::Vertex_handle        Vertex_handle;
  typedef typename Tr::Cell_handle          Cell_handle;
  typedef typename Tr::Cell_circulator      Cell_circulator;

  Tr& tr = c3t3.triangulation();
  const Vertex_handle v1 = e.first->vertex(e.second);
  const Vertex_handle v2 = e.first->vertex(e.third);

  const Point m = tr.geom_traits().construct_midpoint_3_object()
    (point(v1->point()), point(v2->point()));

  //backup subdomain info of incident cells before making changes
  short dimension = 0;
  if(c3t3.is_in_complex(e))
    dimension = 1;
  else
  {
    const std::size_t nb_patches = nb_incident_surface_patches(e, c3t3);
    if(nb_patches == 1)
      dimension = 2;
    else if(nb_patches == 0)
      dimension = 3;
    else
      CGAL_assertion(false);//e should be in complex
  }
  CGAL_assertion(dimension > 0);

  boost::unordered_map<Facet, Subdomain_index> cells_info;
  boost::unordered_map<Facet, std::pair<Vertex_handle, Surface_patch_index> > facets_info;

  // check orientation and collect incident cells to avoid circulating twice
  boost::container::small_vector<Cell_handle, 30> inc_cells;
  Cell_circulator circ = tr.incident_cells(e);
  Cell_circulator end = circ;
  do
  {
    inc_cells.push_back(circ);
    if (tr.is_infinite(circ))
    {
      ++circ;
      continue;
    }

    //1st half-cell
    std::array<Point, 4> pts = { point(circ->vertex(0)->point()),
                                 point(circ->vertex(1)->point()),
                                 point(circ->vertex(2)->point()),
                                 point(circ->vertex(3)->point()) };
    const int i1 = circ->index(v1);
    const Point p1 = pts[i1];
    pts[i1] = m;
    if(CGAL::orientation(pts[0], pts[1], pts[2], pts[3]) != CGAL::POSITIVE)
      return Vertex_handle();

    //2nd half-cell
    pts[i1] = p1;
    pts[circ->index(v2)] = m;
    if (CGAL::orientation(pts[0], pts[1], pts[2], pts[3]) != CGAL::POSITIVE)
      return Vertex_handle();

    ++circ;
  }
  while (circ != end);

  for(Cell_handle c : inc_cells)
  {
    const int index_v1 = c->index(v1);
    const int index_v2 = c->index(v2);

    //keys are the opposite facets to the ones not containing e,
    //because they will not be modified
    const Subdomain_index subdomain = c3t3.subdomain_index(c);
    const Facet opp_facet1 = tr.mirror_facet(Facet(c, index_v1));
    const Facet opp_facet2 = tr.mirror_facet(Facet(c, index_v2));

    // volume data
    cells_info.insert(std::make_pair(opp_facet1, subdomain));
    cells_info.insert(std::make_pair(opp_facet2, subdomain));
    if (c3t3.is_in_complex(c))
      c3t3.remove_from_complex(c);

    // surface data for facets of the cells to be split
    const int findex = CGAL::Triangulation_utils_3::next_around_edge(index_v1, index_v2);
    Surface_patch_index patch = c3t3.surface_patch_index(c, findex);
    Vertex_handle opp_vertex = c->vertex(findex);
    facets_info.insert(std::make_pair(opp_facet1,
                                      std::make_pair(opp_vertex, patch)));
    facets_info.insert(std::make_pair(opp_facet2,
                                      std::make_pair(opp_vertex, patch)));

    if(c3t3.is_in_complex(c, findex))
      c3t3.remove_from_complex(c, findex);
  }

  // insert midpoint
  Vertex_handle new_v = tr.tds().insert_in_edge(e);
  new_v->set_point(typename Tr::Point(m));
  new_v->set_dimension(dimension);

  // update c3t3 with subdomain and surface patch indices
  std::vector<Cell_handle> new_cells;
  tr.incident_cells(new_v, std::back_inserter(new_cells));
  for (Cell_handle new_cell : new_cells)
  {
    const Facet fi(new_cell, new_cell->index(new_v));
    const Facet mfi = tr.mirror_facet(fi);

    //get subdomain info back
    CGAL_assertion(cells_info.find(mfi) != cells_info.end());
    Subdomain_index n_index = cells_info.at(mfi);
    if (Subdomain_index() != n_index)
      c3t3.add_to_complex(new_cell, n_index);
    else
      new_cell->set_subdomain_index(Subdomain_index());

    // get surface info back
    CGAL_assertion(facets_info.find(mfi) != facets_info.end());
    const std::pair<Vertex_handle, Surface_patch_index> v_and_opp_patch = facets_info.at(mfi);

    // facet opposite to new_v (status wrt c3t3 is unchanged)
    new_cell->set_surface_patch_index(new_cell->index(new_v),
                                      mfi.first->surface_patch_index(mfi.second));

    // new half-facet (added or not to c3t3 depending on the stored surface patch index)
    if (Surface_patch_index() == v_and_opp_patch.second)
      new_cell->set_surface_patch_index(new_cell->index(v_and_opp_patch.first),
                                        Surface_patch_index());
    else
      c3t3.add_to_complex(new_cell,
                          new_cell->index(v_and_opp_patch.first),
                          v_and_opp_patch.second);

    // newly created internal facet
    for (int i = 0; i < 4; ++i)
    {
      const Vertex_handle vi = new_cell->vertex(i);
      if (vi == v1 || vi == v2)
      {
        new_cell->set_surface_patch_index(i, Surface_patch_index());
        break;
      }
    }

    //the 4th facet (new_v, v_and_opp_patch.first, v1 or v2)
    // will have its patch tagged from the other side, if needed
  }

  set_index(new_v, c3t3);

  return new_v;
}

template<typename C3T3, typename CellSelector>
bool can_be_split(const typename C3T3::Edge& e,
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
void split_long_edges(C3T3& c3t3,
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
    if (!can_be_split(e, c3t3, protect_boundaries, cell_selector))
      continue;

    typename Gt::Compute_squared_length_3 sql
      = tr.geom_traits().compute_squared_length_3_object();
    FT sqlen = sql(tr.segment(e));
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
      if (!can_be_split(edge, c3t3, protect_boundaries, cell_selector))
        continue;

      visitor.before_split(tr, edge);
      Vertex_handle vh = split_edge(edge, c3t3);

      if(vh != Vertex_handle())
        visitor.after_split(tr, vh);

#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
      if (vh != Vertex_handle())
        ofs << vh->point() << std::endl;
#endif

#if  defined(CGAL_TETRAHEDRAL_REMESHING_VERBOSE_PROGRESS) \
|| defined(CGAL_TETRAHEDRAL_REMESHING_VERBOSE)
      if (vh != Vertex_handle())
        ++nb_splits;
#endif

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE_PROGRESS
      std::cout << "\rSplit (" << high << ")... ("
                << long_edges.left.size() << " long edges, "
                << "length  = " << std::sqrt(sqlen) << ", "
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

} // internal
} // Tetrahedral_remeshing
} // CGAL

#endif // CGAL_INTERNAL_SPLIT_LONG_EDGES_H
