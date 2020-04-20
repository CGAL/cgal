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
// Author(s)     : Jane Tournois, Noura Faraj

#ifndef CGAL_INTERNAL_COLLAPSE_SHORT_EDGES_H
#define CGAL_INTERNAL_COLLAPSE_SHORT_EDGES_H

#include <CGAL/license/Tetrahedral_remeshing.h>

#include <boost/bimap.hpp>
#include <boost/bimap/set_of.hpp>
#include <boost/bimap/multiset_of.hpp>
#include <boost/array.hpp>
#include <boost/unordered_set.hpp>
#include <boost/container/small_vector.hpp>

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
      for (const Cell_handle ch : find_incident)
      {
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
      for (Cell_handle ch : cells_to_update)
      {
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
      for (Cell_handle ch : cells_to_remove){
        triangulation.tds().delete_cell(ch);
      }

      for (Cell_handle cit : triangulation.finite_cell_handles())
      {
        if (!is_well_oriented(triangulation, cit))
          return ORIENTATION_PROBLEM;
        if (!triangulation.tds().is_valid(cit, true))
          return C_PROBLEM;
      }

      for (Vertex_handle vit : triangulation.finite_vertex_handles())
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

//template<typename C3t3>
//Edge_type get_edge_type(const typename C3t3::Edge& edge,
//                        const C3t3& c3t3)
//{
//  typedef typename C3t3::Vertex_handle Vertex_handle;
//  typedef typename C3t3::Triangulation::Cell_circulator Cell_circulator;
//  typedef typename C3t3::Subdomain_index Subdomain_index;

//  const Vertex_handle & v0 = edge.first->vertex(edge.second);
//  const Vertex_handle & v1 = edge.first->vertex(edge.third);

//  const int dim0 = c3t3.in_dimension(v0);
//  const int dim1 = c3t3.in_dimension(v1);

//  const bool is_v0_on_hull = is_on_convex_hull(v0, c3t3);
//  const bool is_v1_on_hull = is_on_convex_hull(v1, c3t3);

//  if (c3t3.is_in_complex(edge))
//    return FEATURE;

//  else if (dim0 == 3 && dim1 == 3)
//    return INSIDE;

//  else if (dim0 == 2 && dim1 == 2)
//  {
//    Cell_circulator circ = c3t3.triangulation().incident_cells(edge);
//    Cell_circulator done = circ;

//    std::vector<Subdomain_index> indices;
//    do
//    {
//      Subdomain_index current_si = circ->subdomain_index();

//      if (std::find(indices.begin(), indices.end(), current_si) == indices.end()) {
//        indices.push_back(current_si);
//      }

//      Subdomain_index si_n0 = circ->neighbor(circ->index(v0))->subdomain_index();
//      Subdomain_index si_n1 = circ->neighbor(circ->index(v1))->subdomain_index();
//      if (si_n0 == si_n1 && si_n0 != current_si)
//        return NO_COLLAPSE;

//    } while (++circ != done);

//    const std::size_t nb_si_v0 = nb_incident_subdomains(v0, c3t3);
//    const std::size_t nb_si_v1 = nb_incident_subdomains(v1, c3t3);

//    if (indices.size() >= (std::min)(nb_si_v0, nb_si_v1)) {
//      return BOUNDARY;
//    }
//  }

//  //std::cerr << "ERROR : get_edge_type did not return anything valid!" << std::endl;
//  return NO_COLLAPSE;
//}

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

    if (n0_ch->has_vertex(v0)
        || n1_ch->has_vertex(v1)
        || n0_ch->has_neighbor(n1_ch))
    {
#ifdef CGAL_DEBUG_TET_REMESHING_IN_PLUGIN
      if (c3t3.is_in_complex(edge))
        ++nb_invalid_collapse_short;
#endif
      return false;
    }
  }
  while (++circ != done);

  return true;
}

template<typename C3t3>
bool is_valid_collapse(const typename C3t3::Edge& edge,
                       const Collapse_type& collapse_type,
                       const typename C3t3::Triangulation::Point& new_pos,
                       const C3t3& c3t3)
{
  typedef typename C3t3::Vertex_handle        Vertex_handle;
  typedef typename C3t3::Cell_handle          Cell_handle;
  typedef typename C3t3::Triangulation::Point Point;

  const Vertex_handle v0 = edge.first->vertex(edge.second);
  const Vertex_handle v1 = edge.first->vertex(edge.third);

#ifdef CGAL_DEBUG_TET_REMESHING_IN_PLUGIN
  const bool in_cx = c3t3.is_in_complex(edge);
  if (in_cx)
  {
    if (collapse_type == TO_MIDPOINT)
      nb_test_midpoint++;
    else if (collapse_type == TO_V1)
      nb_test_v1++;
    else
      nb_test_v0++;
  }
#endif

  if (collapse_type == TO_V1 || collapse_type == TO_MIDPOINT)
  {
    std::vector<Cell_handle> cells_to_check;
    c3t3.triangulation().finite_incident_cells(v0,
        std::back_inserter(cells_to_check));

    for (const Cell_handle ch : cells_to_check)
    {
      if (!ch->has_vertex(v1))
      {
        //check orientation
        boost::array<Point, 4> pts = { ch->vertex(0)->point(),
                                       ch->vertex(1)->point(),
                                       ch->vertex(2)->point(),
                                       ch->vertex(3)->point()};
        pts[ch->index(v0)] = new_pos;
        if (CGAL::orientation(point(pts[0]), point(pts[1]), point(pts[2]), point(pts[3]))
            != CGAL::POSITIVE)
        {
#ifdef CGAL_DEBUG_TET_REMESHING_IN_PLUGIN
          if (in_cx)
          {
            if (collapse_type == TO_MIDPOINT)
              nb_orientation_midpoint++;
            else
              nb_orientation_v1++;
          }
#endif
          return false;
        }
      }
    }
  }
  if (collapse_type == TO_V0 || collapse_type == TO_MIDPOINT)
  {
    std::vector<Cell_handle> cells_to_check;
    c3t3.triangulation().finite_incident_cells(v1,
        std::back_inserter(cells_to_check));

    for (const Cell_handle ch : cells_to_check)
    {
      if (!ch->has_vertex(v0))
      {
        //check orientation
        boost::array<Point, 4> pts = { ch->vertex(0)->point(),
                                       ch->vertex(1)->point(),
                                       ch->vertex(2)->point(),
                                       ch->vertex(3)->point() };
        pts[ch->index(v1)] = new_pos;
        if (CGAL::orientation(point(pts[0]), point(pts[1]), point(pts[2]), point(pts[3]))
            != CGAL::POSITIVE)
        {
#ifdef CGAL_DEBUG_TET_REMESHING_IN_PLUGIN
          if (in_cx)
          {
            if (collapse_type == TO_MIDPOINT)
              nb_orientation_midpoint++;
            else
              nb_orientation_v0++;
          }
#endif
          return false;
        }
      }
    }
  }

  return is_valid_collapse(edge, c3t3);
}

template<typename C3t3, typename CellSelector>
bool are_edge_lengths_valid(const typename C3t3::Edge& edge,
                            const C3t3& c3t3,
                            const typename C3t3::Triangulation::Point& new_pos,
                            const typename C3t3::Triangulation::Geom_traits::FT& sqhigh,
                            const CellSelector& cell_selector,
                            const bool /* adaptive */ = false)
{
  //SqLengthMap::key_type is Vertex_handle
  //SqLengthMap::value_type is double
  typedef typename C3t3::Triangulation::Geom_traits::FT FT;
  typedef typename C3t3::Edge                           Edge;
  typedef typename C3t3::Vertex_handle                  Vertex_handle;

  const Vertex_handle v1 = edge.first->vertex(edge.second);
  const Vertex_handle v2 = edge.first->vertex(edge.third);

  boost::unordered_map<Vertex_handle, FT> edges_sqlength_after_collapse;

  std::vector<Edge> inc_edges;
  c3t3.triangulation().finite_incident_edges(v1,
      std::back_inserter(inc_edges));
  c3t3.triangulation().finite_incident_edges(v2,
      std::back_inserter(inc_edges));

  for (const Edge& ei : inc_edges)
  {
    if (is_outside(ei, c3t3, cell_selector))
      continue;

    Vertex_handle vh = ei.first->vertex(ei.second);
    if (vh == v1 || vh == v2)
      vh = ei.first->vertex(ei.third);
    if (vh == v1 || vh == v2)
      continue;

    if (edges_sqlength_after_collapse.find(vh) == edges_sqlength_after_collapse.end())
    {
      const FT sqlen = CGAL::squared_distance(new_pos, point(vh->point()));

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

      if (sqlen > sqhigh) {
        return false;
      }
      //}
      edges_sqlength_after_collapse[vh] = sqlen;
    }
  }

  return true;
}

template<typename C3t3>
void merge_surface_patch_indices(const typename C3t3::Facet& f1,
                                 const typename C3t3::Facet& f2,
                                 C3t3& c3t3)
{
  const bool in_cx_f1 = c3t3.is_in_complex(f1);
  const bool in_cx_f2 = c3t3.is_in_complex(f2);

  if (in_cx_f1 && !in_cx_f2)
  {
    typename C3t3::Surface_patch_index patch = c3t3.surface_patch_index(f1);
    f2.first->set_surface_patch_index(f2.second, patch);
  }
  else if (in_cx_f2 && !in_cx_f1)
  {
    typename C3t3::Surface_patch_index patch = c3t3.surface_patch_index(f2);
    f1.first->set_surface_patch_index(f1.second, patch);
  }
  else if(in_cx_f1 && in_cx_f2)
  {
    CGAL_assertion(c3t3.surface_patch_index(f1) == c3t3.surface_patch_index(f2));

    typename C3t3::Surface_patch_index patch = c3t3.surface_patch_index(f2);
    c3t3.remove_from_complex(f2);
    f2.first->set_surface_patch_index(f2.second, patch);
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


  //Update the vertex before removing it
  std::vector<Cell_handle> find_incident;
  tr.incident_cells(vh0, std::back_inserter(find_incident));

  std::vector<Cell_handle> cells_to_update;
  tr.incident_cells(vh1, std::back_inserter(cells_to_update));

  boost::container::small_vector<Cell_handle, 30> inc_cells;
  Cell_circulator circ = tr.incident_cells(ch, to, from);
  Cell_circulator done = circ;
  do
  {
    for (int i = 0; i < 4; ++i)
    {
      const Vertex_handle vi = circ->vertex(i);
      if (vi != vh0 && vi != vh1)
      {
        const Facet fi(circ, i);
        if (c3t3.is_in_complex(fi))
          c3t3.remove_from_complex(fi);
      }
    }
    inc_cells.push_back(circ);
  }
  while (++circ != done);

  bool valid = true;
  std::vector<Cell_handle> cells_to_remove;
  boost::unordered_set<Cell_handle> invalid_cells;

  for(const Cell_handle c : inc_cells)
  {
    const int v0_id = c->index(vh0);
    const int v1_id = c->index(vh1);

    Cell_handle n0_ch = c->neighbor(v0_id);
    Cell_handle n1_ch = c->neighbor(v1_id);

    const int ch_id_in_n0 = n0_ch->index(c);
    const int ch_id_in_n1 = n1_ch->index(c);

    //Merge surface patch indices
    merge_surface_patch_indices(Facet(n0_ch, ch_id_in_n0),
                                Facet(n1_ch, ch_id_in_n1),
                                c3t3);

    //Update neighbors before removing cell
    n0_ch->set_neighbor(ch_id_in_n0, n1_ch);
    n1_ch->set_neighbor(ch_id_in_n1, n0_ch);

    //Update vertices cell pointer
    for (int i = 0; i < 3; i++)
    {
      int vid = Tr::vertex_triple_index(ch_id_in_n0, i);
      n0_ch->vertex(vid)->set_cell(n0_ch);
    }
    for (int i = 0; i < 3; i++)
    {
      int vid = Tr::vertex_triple_index(ch_id_in_n1, i);
      n1_ch->vertex(vid)->set_cell(n1_ch);
    }

    if (tr.is_infinite(n0_ch->vertex(ch_id_in_n0))
      && tr.is_infinite(n1_ch->vertex(ch_id_in_n1)))
    {
      std::cout << "Collapse infinite issue!" << std::endl;
      return Vertex_handle();
    }
    cells_to_remove.push_back(c);

    invalid_cells.insert(c);
  }

  const Vertex_handle infinite_vertex = tr.infinite_vertex();

  bool v0_updated = false;
  for (const Cell_handle c : find_incident)
  {
    if (invalid_cells.find(c) == invalid_cells.end())//valid cell
    {
      if (tr.is_infinite(c))
        infinite_vertex->set_cell(c);
      //else {
      vh0->set_cell(c);
      v0_updated = true;
      //}
    }
  }

  // update complex edges
  const std::array<std::array<int, 2>, 6> edges
    = { { 0,1, 0,2, 0,3, 1,2, 1,3, 2,3 } }; //vertex indices in cells
  const Vertex_handle vkept = vh0;
  const Vertex_handle vdeleted = vh1;
  for (const Cell_handle c : cells_to_update)
  {
    for (const std::array<int, 2>& ei : edges)
    {
      Vertex_handle eiv0 = c->vertex(ei[0]);
      Vertex_handle eiv1 = c->vertex(ei[1]);
      if (eiv1 == vdeleted && eiv0 != vkept) //replace eiv1 by vkept
      {
        if (c3t3.is_in_complex(eiv0, eiv1))
        {
          c3t3.add_to_complex(eiv0, vkept, c3t3.curve_index(eiv0, eiv1));
          c3t3.remove_from_complex(eiv0, eiv1);
        }
      }
      else if (eiv0 == vdeleted && eiv1 != vkept) //replace eiv0 by vkept
      {
        if (c3t3.is_in_complex(eiv0, eiv1))
        {
          c3t3.add_to_complex(vkept, eiv1, c3t3.curve_index(eiv0, eiv1));
          c3t3.remove_from_complex(eiv0, eiv1);
        }
      }
    }
  }

  // update complex facets

  //Update the vertex before removing it
  for (const Cell_handle c : cells_to_update)
  {
    if (invalid_cells.find(c) == invalid_cells.end()) //valid cell
    {
      c->set_vertex(c->index(vh1), vh0);

      if (tr.is_infinite(c))
        infinite_vertex->set_cell(c);
      //else {
      if (!v0_updated) {
        vh0->set_cell(c);
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
    // remove cell
    if (c3t3.is_in_complex(cell_to_remove))
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

  const int dim_vh0 = c3t3.in_dimension(vh0);
  const int dim_vh1 = c3t3.in_dimension(vh1);

  Vertex_handle vh = Vertex_handle();

  const Point_3 p0 = vh0->point();
  const Point_3 p1 = vh1->point();

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
    const bool /* protect_boundaries */,
    CellSelector cell_selector,
    Visitor& visitor)
{
  typedef typename C3t3::Triangulation   Tr;
  typedef typename Tr::Point             Point;
  typedef typename Tr::Vertex_handle     Vertex_handle;

  const Vertex_handle v0 = edge.first->vertex(edge.second);
  const Vertex_handle v1 = edge.first->vertex(edge.third);

  Collapse_type collapse_type = get_collapse_type(edge, c3t3, cell_selector);

#ifdef CGAL_DEBUG_TET_REMESHING_IN_PLUGIN
  const bool in_cx = c3t3.is_in_complex(edge);
  if (in_cx && collapse_type == IMPOSSIBLE)
    nb_impossible++;
#endif

  if (collapse_type == IMPOSSIBLE)
    return Vertex_handle();

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

  if (!is_valid_collapse(edge, collapse_type, new_pos, c3t3))
  {
#ifdef TET_REMESHING_COLLAPSE_FALLBACK_EXPERIMENTS
    if (collapse_type == TO_MIDPOINT)
    {
      // with TO_MIDPOINT, we are authorized to test TO_V0 and TO_V1
      if (is_valid_collapse(edge, TO_V0, v0->point(), c3t3))
      {
        collapse_type = TO_V0;
        new_pos = v0->point();
      }
      else if (is_valid_collapse(edge, TO_V1, v1->point(), c3t3))
      {
        collapse_type = TO_V1;
        new_pos = v1->point();
      }
      else
      {
#ifdef CGAL_DEBUG_TET_REMESHING_IN_PLUGIN
        if (in_cx)
          nb_invalid_collapse++;
#endif
        return Vertex_handle();
      }
    }
    else
#endif //TET_REMESHING_COLLAPSE_FALLBACK_EXPERIMENTS
    {
#ifdef CGAL_DEBUG_TET_REMESHING_IN_PLUGIN
      if (in_cx)
        nb_invalid_collapse++;
#endif
      return Vertex_handle();
    }
  }

  if (are_edge_lengths_valid(edge, c3t3, new_pos, sqhigh, cell_selector/*, adaptive = false*/))
  {
    CollapseTriangulation<C3t3, Visitor> local_tri(c3t3, edge, collapse_type, visitor);
    local_tri.update();

    Result_type res = local_tri.collapse();
    if (res == VALID)
    {
#ifdef CGAL_DEBUG_TET_REMESHING_IN_PLUGIN
      if (in_cx)
        nb_valid_collapse++;
#endif
      return collapse(edge, collapse_type, c3t3);
    }
  }
#ifdef CGAL_DEBUG_TET_REMESHING_IN_PLUGIN
  else if (in_cx)
    nb_invalid_lengths++;
#endif
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

  T3& tr = c3t3.triangulation();
  typename Gt::Compute_squared_length_3 sql
    = tr.geom_traits().compute_squared_length_3_object();

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
  std::cout << "Collapse short edges (" << low << ", " << high << ")...";
  std::cout.flush();
  std::size_t nb_collapses = 0;
#endif
  const FT sq_low = low*low;
  const FT sq_high = high*high;

  //collect long edges
  Boost_bimap short_edges;
  for (Finite_edges_iterator eit = tr.finite_edges_begin();
       eit != tr.finite_edges_end(); ++eit)
  {
    const Edge& e = *eit;
    if (!can_be_collapsed(e, c3t3, protect_boundaries, cell_selector))
      continue;

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
    FT sqlen = eit->first;
    std::cout << "\rCollapse... (" << short_edges.left.size() << " short edges, ";
    std::cout << std::sqrt(sqlen) << ", ";
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

      Vertex_handle vh = collapse_edge(edge, c3t3, sq_high,
                                       protect_boundaries, cell_selector,
                                       visitor);
      if (vh != Vertex_handle())
      {
        std::vector<Edge> incident_short;
        c3t3.triangulation().finite_incident_edges(vh,
            std::back_inserter(incident_short));
        for (const Edge& eshort : incident_short)
        {
          if (!can_be_collapsed(eshort, c3t3, protect_boundaries, cell_selector))
            continue;

          const FT sqlen = sql(tr.segment(eshort));
          if (sqlen < sq_low)
            short_edges.insert(short_edge(make_vertex_pair<T3>(eshort), sqlen));
        }

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
        ++nb_collapses;
#endif
      }
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
