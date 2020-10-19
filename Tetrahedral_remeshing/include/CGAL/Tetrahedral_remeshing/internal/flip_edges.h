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

#ifndef CGAL_INTERNAL_FLIP_EDGES_H
#define CGAL_INTERNAL_FLIP_EDGES_H

#include <CGAL/license/Tetrahedral_remeshing.h>

#include <CGAL/Triangulation_utils_3.h>
#include <CGAL/utility.h>

#include <CGAL/Tetrahedral_remeshing/internal/tetrahedral_remeshing_helpers.h>

#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>
#include <boost/container/small_vector.hpp>
#include <boost/optional.hpp>

#include <limits>
#include <queue>

namespace CGAL
{
namespace Tetrahedral_remeshing
{
namespace internal
{
enum Flip_Criterion{ MIN_ANGLE_BASED, AVERAGE_ANGLE_BASED,
                     VALENCE_BASED, VALENCE_MIN_DH_BASED };

//template<typename VertexPair>
//void flip_inside_edges(std::vector<VertexPair>&)
//{
//  //TODO
//}

//outer_mirror_facets contains the set of facets of the outer hull
//of the set of cells modified by the flip operation,
//"seen from" outside
//i.e. for each facet f among those, f.first has not been modified by flip
template<typename C3t3, typename CellSet, typename FacetSet>
void update_c3t3_facets(C3t3& c3t3,
                        const CellSet& cells_to_update,
                        const FacetSet& outer_mirror_facets)
{
  typedef typename C3t3::Facet       Facet;
  typedef typename C3t3::Cell_handle Cell_handle;
  typedef typename C3t3::Surface_patch_index Surface_patch_index;

  for (Cell_handle c : cells_to_update)
  {
    //their subdomain indices have not been modified because we kept the same cells
    //surface patch indices need to be fixed though
    for (int i = 0; i < 4; ++i)
    {
      const Facet f(c, i);
      const Facet mf = c3t3.triangulation().mirror_facet(f);
      if (outer_mirror_facets.find(mf) != outer_mirror_facets.end())
      {
        //we are on the border of the modified zone, c3t3 info is valid outside,
        //on mirror facet
        const typename C3t3::Surface_patch_index patch = c3t3.surface_patch_index(mf);
        if (c3t3.is_in_complex(mf))
          f.first->set_surface_patch_index(f.second, patch);
        else
          f.first->set_surface_patch_index(f.second, Surface_patch_index());
      }
      else
      {
        //we are inside the modified zone, c3t3 info is not valid anymore
        if (c3t3.is_in_complex(f) || c3t3.is_in_complex(mf))
        {
          f.first->set_surface_patch_index(f.second, Surface_patch_index());
          mf.first->set_surface_patch_index(mf.second, Surface_patch_index());
        }
      }
    }
  }
}

template<typename C3t3, typename IncCellsVectorMap>
Sliver_removal_result flip_3_to_2(typename C3t3::Edge& edge,
                                  C3t3& c3t3,
                                  const std::vector<typename C3t3::Vertex_handle>& vertices_around_edge,
                                  const Flip_Criterion& criterion,
                                  IncCellsVectorMap& inc_cells)
{
  typedef typename C3t3::Triangulation Tr;
  typedef typename C3t3::Facet         Facet;
  typedef typename C3t3::Vertex_handle Vertex_handle;
  typedef typename C3t3::Cell_handle   Cell_handle;
  typedef typename Tr::Cell_circulator Cell_circulator;
  typedef typename Tr::Geom_traits     Gt;
  typedef typename Gt::FT              FT;

  //Edge to face flip
  Tr& tr = c3t3.triangulation();

  Cell_circulator circ = tr.incident_cells(edge);
  Cell_circulator done = circ;

  Vertex_handle vh0 = edge.first->vertex(edge.second);
  Vertex_handle vh1 = edge.first->vertex(edge.third);

  //Select 2 cells to keep and update and one to remove
  Cell_handle ch0 = Cell_handle(circ++);
  Cell_handle ch1 = Cell_handle(circ++);
  Cell_handle cell_to_remove = Cell_handle(circ++);
  if (circ != done)
  {
    std::cout << "Wrong flip function" << std::endl;
    return NOT_FLIPPABLE;
  }

  //Check structural validity
  Cell_handle c;
  int i0, i1, i3;
  if (tr.is_facet(vertices_around_edge[0], vertices_around_edge[1], vertices_around_edge[2],
                  c, i0, i1, i3))
    return NOT_FLIPPABLE;

  //Check topological validity
  const typename C3t3::Subdomain_index subdomain = ch0->subdomain_index();
  if ( subdomain != ch1->subdomain_index()
       || subdomain != cell_to_remove->subdomain_index()
       || ch1->subdomain_index() != cell_to_remove->subdomain_index())
    return NOT_FLIPPABLE;

  Vertex_handle vh2;
  Vertex_handle vh3;

  for (int i = 0; i < 3; ++i){
    if (!ch0->has_vertex(vertices_around_edge[i]))
      vh2 = vertices_around_edge[i];
    else if (!ch1->has_vertex(vertices_around_edge[i]))
      vh3 = vertices_around_edge[i];
  }

  int vh0_id = ch0->index(vh0);
  int vh1_id = ch1->index(vh1);

  //Check if flip valid
  if (!is_well_oriented(tr, vh2,
                        ch0->vertex(indices(vh0_id, 0)),
                        ch0->vertex(indices(vh0_id, 1)),
                        ch0->vertex(indices(vh0_id, 2)))
      || !is_well_oriented(tr, vh3,
                           ch1->vertex(indices(vh1_id, 0)),
                           ch1->vertex(indices(vh1_id, 1)),
                           ch1->vertex(indices(vh1_id, 2))))
    return NOT_FLIPPABLE;

  ///********************VALIDITY CHECK***************************/
  //double curr_min_dh;
  //bool check_validity = false;
  //std::vector<typename Tr::Tetrahedron> pre_sliver_Removal_cells;
  //if (check_validity){
  //  pre_sliver_Removal_cells.clear();
  //  pre_sliver_Removal_cells.push_back(K::Tetrahedron_3(ch0->vertex(0)->point(), ch0->vertex(1)->point(), ch0->vertex(2)->point(), ch0->vertex(3)->point()));
  //  pre_sliver_Removal_cells.push_back(K::Tetrahedron_3(ch1->vertex(0)->point(), ch1->vertex(1)->point(), ch1->vertex(2)->point(), ch1->vertex(3)->point()));
  //  pre_sliver_Removal_cells.push_back(K::Tetrahedron_3(cell_to_remove->vertex(0)->point(), cell_to_remove->vertex(1)->point(),
  //    cell_to_remove->vertex(2)->point(), cell_to_remove->vertex(3)->point()));

  //  curr_min_dh = min_dihedral_angle<Gt>(ch0);
  //  curr_min_dh = std::min(curr_min_dh, min_dihedral_angle<Gt>(ch1));
  //  curr_min_dh = std::min(curr_min_dh, min_dihedral_angle<Gt>(cell_to_remove));

  //  pre_sliver_Removal_vertices.clear();
  //  for (int i = 0; i < vertices_around_edge.size(); ++i){
  //    pre_sliver_Removal_vertices.push_back(Point_3(vertices_around_edge[i]->point()));
  //  }

  //  previous_edges.clear();
  //  previous_edges.push_back(std::make_pair(vh0->point(), vh1->point()));
  //}
  /*************************************************************/


  if (criterion == MIN_ANGLE_BASED)
  {
    //Current worst dihedral angle
    Dihedral_angle_cosine curr_max_cosdh = max_cos_dihedral_angle(tr, ch0);
    curr_max_cosdh = (std::max)(curr_max_cosdh, max_cos_dihedral_angle(tr, ch1));
    curr_max_cosdh = (std::max)(curr_max_cosdh, max_cos_dihedral_angle(tr, cell_to_remove));

    //Result worst dihedral angle
    if (curr_max_cosdh < max_cos_dihedral_angle(tr, vh2,
                                         ch0->vertex(indices(vh0_id, 0)),
                                         ch0->vertex(indices(vh0_id, 1)),
                                         ch0->vertex(indices(vh0_id, 2)))
        || curr_max_cosdh < max_cos_dihedral_angle(tr, vh3,
                                         ch1->vertex(indices(vh1_id, 0)),
                                         ch1->vertex(indices(vh1_id, 1)),
                                         ch1->vertex(indices(vh1_id, 2))))
      return NO_BEST_CONFIGURATION;
  }
  else if (criterion == AVERAGE_ANGLE_BASED)
  {
    //Current worst dihedral angle
    double average_min_dh = min_dihedral_angle(tr, ch0);
    average_min_dh += min_dihedral_angle(tr, ch1);
    average_min_dh += min_dihedral_angle(tr, cell_to_remove);

    average_min_dh /= 3.;

    FT new_average_min_dh = 0.5 *
                            (min_dihedral_angle(tr, vh2, ch0->vertex(indices(vh0_id, 0)),
                                ch0->vertex(indices(vh0_id, 1)),
                                ch0->vertex(indices(vh0_id, 2)))
                           + min_dihedral_angle(tr, vh3, ch1->vertex(indices(vh1_id, 0)),
                                 ch1->vertex(indices(vh1_id, 1)),
                                 ch1->vertex(indices(vh1_id, 2))));
    //Result worst dihedral angle
    if (average_min_dh > new_average_min_dh)
      return NO_BEST_CONFIGURATION;
  }

  //Keep the facets
  typedef CGAL::Triple<Vertex_handle, Vertex_handle, Vertex_handle> Facet_vvv;
  typedef boost::unordered_map<Facet_vvv, std::size_t> FaceMapIndex;
  boost::unordered_set<Facet> outer_mirror_facets;

  FaceMapIndex facet_map_indices;
  std::vector<Facet> mirror_facets;
  circ = Cell_circulator(done);
  do
  {
    // facet opposite to vh0
    int curr_vh0_id = circ->index(vh0);
    Facet n_vh0_facet = tr.mirror_facet(Facet(circ, curr_vh0_id));

    outer_mirror_facets.insert(n_vh0_facet);

    Facet_vvv face0 = make_vertex_triple(circ->vertex(indices(curr_vh0_id, 0)),
                                         circ->vertex(indices(curr_vh0_id, 1)),
                                         circ->vertex(indices(curr_vh0_id, 2)));

    typename FaceMapIndex::iterator it = facet_map_indices.find(face0);
    if (it == facet_map_indices.end())
    {
      facet_map_indices[face0] = mirror_facets.size();
      mirror_facets.push_back(n_vh0_facet);
    }

    // facet opposite to vh1
    int curr_vh1_id = circ->index(vh1);
    Facet n_vh1_facet = tr.mirror_facet(Facet(circ, curr_vh1_id));

    outer_mirror_facets.insert(n_vh1_facet);

    Facet_vvv face1 = make_vertex_triple(circ->vertex(indices(curr_vh1_id, 0)),
                                         circ->vertex(indices(curr_vh1_id, 1)),
                                         circ->vertex(indices(curr_vh1_id, 2)));
    it = facet_map_indices.find(face1);
    if (it == facet_map_indices.end())
    {
      facet_map_indices[face1] = mirror_facets.size();
      mirror_facets.push_back(n_vh1_facet);
    }
  }
  while (++circ != done);

  /*
  c3t3.remove_from_complex( ch0 );
  c3t3.remove_from_complex( ch1 );
  c3t3.remove_from_complex( cell_to_remove );

  tr.flip(edge);

  for( int i = 0 ; i < facets.size() ; i ++ ){
  Cell_handle new_cell = facets[i].first->neighbor( facets[i].second );
  c3t3.add_to_complex( new_cell, si );
  }
  */

  //Update cells
  ch0->set_vertex(vh0_id, vh2);
  ch1->set_vertex(vh1_id, vh3);

  // "New" cells are not created, only modified/updated
  std::vector<Cell_handle> cells_to_update;
  cells_to_update.push_back(ch0);
  cells_to_update.push_back(ch1);

  //Update adjacencies and vertices' cells
  for (Cell_handle ch : cells_to_update)
  {
    for (int v = 0; v < 4; ++v)
    {
      Facet_vvv face = make_vertex_triple(ch->vertex(indices(v, 0)),
                                          ch->vertex(indices(v, 1)),
                                          ch->vertex(indices(v, 2)));
      typename FaceMapIndex::iterator it = facet_map_indices.find(face);
      if (it == facet_map_indices.end())
      {
        facet_map_indices[face] = mirror_facets.size();
        mirror_facets.push_back(Facet(ch, v));
      }
      else
      {
        Facet mirror_facet = mirror_facets[it->second];

        //Update neighbor
        mirror_facet.first->set_neighbor(mirror_facet.second, ch);
        ch->set_neighbor(v, mirror_facet.first);
      }
      ch->vertex(v)->set_cell(ch);

      inc_cells[ch->vertex(v)] = boost::none;
      ch->reset_cache_validity();
    }
  }

  // Update c3t3
  update_c3t3_facets(c3t3, cells_to_update, outer_mirror_facets);

  c3t3.remove_from_complex(cell_to_remove);
  tr.tds().delete_cell(cell_to_remove);

  /********************VALIDITY CHECK***************************/
  //if (check_validity)
  //{
  //  post_sliver_Removal_cells.clear();
  //  post_sliver_Removal_cells.push_back(ch0);
  //  post_sliver_Removal_cells.push_back(ch1);

  //  double new_min_dh = min_dihedral_angle<Gt>(ch0);
  //  new_min_dh = std::min(new_min_dh, min_dihedral_angle<Gt>(ch1));

  //  post_sliver_Removal_vertices.clear();
  //  post_sliver_Removal_vertices.push_back(vh2);
  //  post_sliver_Removal_vertices.push_back(vh3);

  //  if (!is_well_oriented(ch0))
  //    return INVALID_ORIENTATION;
  //  if (!is_well_oriented(ch1))
  //    return INVALID_ORIENTATION;
  //  if (!tr.is_valid(ch0))
  //    return INVALID_CELL;
  //  if (!tr.is_valid(ch1))
  //    return INVALID_CELL;

  //  for (int i = 0; i < 4; ++i){
  //    if (!tr.is_valid(ch0->neighbor(i)))
  //      return INVALID_CELL;
  //    if (!tr.is_valid(ch1->neighbor(i)))
  //      return INVALID_CELL;
  //    if (!tr.tds().is_valid(ch0->vertex(i)))
  //      return INVALID_VERTEX;
  //    if (!tr.tds().is_valid(ch1->vertex(i))){
  //      return INVALID_VERTEX;
  //    }
  //  }

  //  if ((curr_min_dh - new_min_dh) > 0.01){
  //    std::cout << "Three_to_two_flip::Flip not improving the quality: " << curr_min_dh << " to " << new_min_dh << std::endl;
  //    return INVALID_CELL;
  //  }
  //}
  /***********************************************************/

  return VALID_FLIP;
}

template<typename C3t3, typename CandidatesQueue>
void find_best_flip_to_improve_dh(C3t3& c3t3,
                                  typename C3t3::Edge& edge,
                                  typename C3t3::Vertex_handle vh2,
                                  typename C3t3::Vertex_handle vh3,
                                  CandidatesQueue& candidates,
                                  const Dihedral_angle_cosine& curr_max_cos_dh,
                                  bool is_sliver_well_oriented = true,
                                  int e_id = 0)
{
  typedef typename C3t3::Triangulation  Tr;
  typedef typename C3t3::Vertex_handle  Vertex_handle;
  typedef typename C3t3::Cell_handle    Cell_handle;
  typedef typename C3t3::Facet          Facet;
  typedef typename Tr::Facet_circulator Facet_circulator;
  typedef typename Tr::Cell_circulator  Cell_circulator;

  // std::cout << "find_best_flip_to_improve_dh boundary " << std::endl;
  Tr& tr = c3t3.triangulation();

  Vertex_handle vh0 = edge.first->vertex(edge.second);
  Vertex_handle vh1 = edge.first->vertex(edge.third);

  Facet_circulator curr_fcirc = tr.incident_facets(edge);
  Facet_circulator curr_fdone = curr_fcirc;

  //Only keep the possible flips
  std::vector<Vertex_handle> opposite_vertices;
  int nb_cells_around_edge = 0;
  do
  {
    Vertex_handle vh;
    //Get the ids of the opposite vertices
    for (int i = 0; i < 3; ++i)
    {
      Vertex_handle curr_vertex = curr_fcirc->first->vertex(indices(curr_fcirc->second, i));
      if ( curr_vertex != vh0
           && curr_vertex != vh1
           && (curr_vertex == vh2 || curr_vertex == vh3))
      {
        vh = curr_vertex;
        Facet_circulator facet_circulator(curr_fcirc);
        Facet_circulator facet_done(curr_fcirc);

        facet_done--;
        facet_circulator++;
        facet_circulator++;

        bool is_edge = false;
        do
        {
          //Get the ids of the opposite vertices
          for (int j = 0; j < 3; ++j)
          {
            Vertex_handle curr = facet_circulator->first->vertex(
                                          indices(facet_circulator->second, j));
            if (curr != vh0  && curr != vh1)
            {
              Cell_handle ch;
              int i0, i1;
              if (tr.is_edge(curr, vh, ch, i0, i1))
                is_edge = true;
            }
          }
        } while (++facet_circulator != facet_done);

        if (!is_edge && !tr.is_infinite(vh))
          opposite_vertices.push_back(vh);
      }
    }
    nb_cells_around_edge++;
  }
  while (++curr_fcirc != curr_fdone);

  if (nb_cells_around_edge < 4)
    return;

  //Facets that will be used to create new cells i.e. all the facets opposite to vh1 and don't have vh
  //Facets that will be used to update cells i.e. all the facets opposite to vh0 will be set to vh: facet.first->set_vertex( facet.second, vh )

  Cell_circulator cell_circulator = tr.incident_cells(edge);
  Cell_circulator done = cell_circulator;

  boost::container::small_vector<Facet, 60> facets;
  for (std::size_t i = 0; i < opposite_vertices.size(); ++i)
  {
    Vertex_handle vh = opposite_vertices[i];
    bool keep = true;

    do
    {
      //Store it if it do not have vh
      if (!cell_circulator->has_vertex(vh))
      {
        //Facets opposite to vh0
        Facet facet_vh0(cell_circulator, cell_circulator->index(vh0));

        //Facets opposite to vh1
        Facet facet_vh1(cell_circulator, cell_circulator->index(vh1));

        facets.push_back(facet_vh1);
        facets.push_back(facet_vh0);
      }
    } while (++cell_circulator != done);


    Dihedral_angle_cosine max_flip_cos_dh(CGAL::NEGATIVE, 1., 1.);
    for (const Facet& fi : facets)
    {
      if (!tr.is_infinite(fi.first))
      {
        if (is_well_oriented(tr, vh, fi.first->vertex(indices(fi.second, 0)),
                             fi.first->vertex(indices(fi.second, 1)),
                             fi.first->vertex(indices(fi.second, 2))))
        {
          max_flip_cos_dh = (std::max)(
            max_flip_cos_dh,
            max_cos_dihedral_angle(tr, vh, fi.first->vertex(indices(fi.second, 0)),
                                           fi.first->vertex(indices(fi.second, 1)),
                                           fi.first->vertex(indices(fi.second, 2))));
        }
        else
        {
          keep = false;
          break;
        }

        if (max_flip_cos_dh.is_one())//it will not get worse than 1.
        {
          keep = false;
          break;
        }
      }
    }
    facets.clear();

    if (keep && (max_flip_cos_dh < curr_max_cos_dh  || !is_sliver_well_oriented))
    {
      //std::cout << "vh " << vh->info() <<" old " << curr_max_cos_dh << " min " << min_flip_tan_dh << std::endl;
      candidates.push(std::make_pair(max_flip_cos_dh, std::make_pair(vh, e_id)));
    }
  }
}

template<typename Vertex_handle, typename CellVector, typename Cell_handle>
bool is_edge_uv(Vertex_handle u,
                Vertex_handle v,
                const CellVector& cells_incident_to_u,
                Cell_handle& cell,
                int& i,
                int& j)
{
  if (u == v)
    return false;

  for (typename CellVector::value_type c : cells_incident_to_u)
  {
    if (c->has_vertex(v, j))
    {
      cell = c;
      i = cell->index(u);
      return true;
    }
  }
  return false;
}

template<typename Vertex_handle, typename CellVector>
bool is_edge_uv(Vertex_handle u,
                Vertex_handle v,
                const CellVector& cells_incident_to_u)
{
  typename CellVector::value_type c;
  int i, j;
  return is_edge_uv(u, v, cells_incident_to_u, c, i, j);
}

template<typename C3t3, typename CandidatesQueue,
         typename IncCellsVectorMap>
void find_best_flip_to_improve_dh(C3t3& c3t3,
                                  typename C3t3::Edge& edge,
                                  CandidatesQueue& candidates,
                                  const Dihedral_angle_cosine& curr_max_cosdh,
                                  IncCellsVectorMap& inc_cells,
                                  bool is_sliver_well_oriented = true,
                                  int e_id = 0)
{
  typedef typename C3t3::Triangulation  Tr;
  typedef typename C3t3::Vertex_handle  Vertex_handle;
  typedef typename C3t3::Cell_handle    Cell_handle;
  typedef typename C3t3::Facet          Facet;
  typedef typename Tr::Facet_circulator Facet_circulator;
  typedef typename Tr::Cell_circulator  Cell_circulator;

  Tr& tr = c3t3.triangulation();

  Vertex_handle vh0 = edge.first->vertex(edge.second);
  Vertex_handle vh1 = edge.first->vertex(edge.third);

  Facet_circulator curr_fcirc = tr.incident_facets(edge);
  Facet_circulator curr_fdone = curr_fcirc;

  //Only keep the possible flips
  std::vector<Vertex_handle> opposite_vertices;
  int nb_cells_around_edge = 0;
  do
  {
    Vertex_handle vh;
    //Get the ids of the opposite vertices
    for (int i = 0; i < 3; ++i)
    {
      Vertex_handle curr_vertex = curr_fcirc->first->vertex(
                                    indices(curr_fcirc->second, i));
      if (curr_vertex != vh0 && curr_vertex != vh1)
      {
        vh = curr_vertex;
        break;
      }
    }

    if(tr.is_infinite(vh))
      continue;

    boost::optional<boost::container::small_vector<Cell_handle, 64>>& o_inc_vh = inc_cells[vh];
    if (o_inc_vh == boost::none)
    {
      boost::container::small_vector<Cell_handle, 64> inc_vec;
      tr.incident_cells(vh, std::back_inserter(inc_vec));
      o_inc_vh = inc_vec;
    }

    Facet_circulator facet_circulator = curr_fcirc;
    Facet_circulator facet_done = curr_fcirc;

    facet_done--;
    facet_circulator++;
    facet_circulator++;
    bool is_edge = false;
    do
    {
      //Get the ids of the opposite vertices
      for (int i = 0; i < 3; ++i)
      {
        Vertex_handle curr_vertex = facet_circulator->first->vertex(
                                      indices(facet_circulator->second, i));
        if (curr_vertex != vh0  && curr_vertex != vh1)
        {
          if (is_edge_uv(vh, curr_vertex, boost::get(o_inc_vh)))
          {
            is_edge = true;
            break;
          }
        }
      }
    } while (++facet_circulator != facet_done);

    if (!is_edge)
      opposite_vertices.push_back(vh);

    nb_cells_around_edge++;
  }
  while (++curr_fcirc != curr_fdone);
  if (nb_cells_around_edge < 4)
    return;

  //Facets that will be used to create new cells
  //    i.e. all the facets opposite to vh1 and don't have vh
  //Facets that will be used to update cells
  //    i.e. all the facets opposite to vh0 will be set to vh:
  //    facet.first->set_vertex( facet.second, vh )

  Cell_circulator cell_circulator = tr.incident_cells(edge);
  Cell_circulator done = cell_circulator;

  boost::container::small_vector<Facet, 60> facets;
  for (Vertex_handle vh : opposite_vertices)
  {
    bool keep = true;
    do
    {
      //Store it if it do not have vh
      if (!cell_circulator->has_vertex(vh))
      {
        //Facets opposite to vh0
        Facet facet_vh0(cell_circulator, cell_circulator->index(vh0));

        //Facets opposite to vh1
        Facet facet_vh1(cell_circulator, cell_circulator->index(vh1));

        facets.push_back(facet_vh1);
        facets.push_back(facet_vh0);
      }
    }
    while (++cell_circulator != done);

    Dihedral_angle_cosine max_flip_cos_dh(CGAL::NEGATIVE, 1., 1.);
    for (const Facet& fi : facets)
    {
      if (!tr.is_infinite(fi.first))
      {
        if (is_well_oriented(tr, vh, fi.first->vertex(indices(fi.second, 0)),
                             fi.first->vertex(indices(fi.second, 1)),
                             fi.first->vertex(indices(fi.second, 2))))
        {
          max_flip_cos_dh = (std::max)(max_flip_cos_dh,
            max_cos_dihedral_angle(tr, vh, fi.first->vertex(indices(fi.second, 0)),
                                           fi.first->vertex(indices(fi.second, 1)),
                                           fi.first->vertex(indices(fi.second, 2))));
        }
        else
        {
          keep = false;
          break;
        }

        if (max_flip_cos_dh.is_one())//it will not get worse than 1.
        {
          keep = false;
          break;
        }
      }
    }
    facets.clear();

    if (keep && (max_flip_cos_dh < curr_max_cosdh || !is_sliver_well_oriented))
    {
      //std::cout << "vh " << vh->info() <<" old " << curr_max_cosdh << " min " << min_flip_tan_dh << std::endl;
      candidates.push(std::make_pair(max_flip_cos_dh, std::make_pair(vh, e_id)));
    }
  }
}

template<typename C3t3, typename IncCellsVectorMap, typename Visitor>
Sliver_removal_result flip_n_to_m(C3t3& c3t3,
                                  typename C3t3::Edge& edge,
                                  typename C3t3::Vertex_handle vh,
                                  IncCellsVectorMap& inc_cells,
                                  Visitor& visitor,
                                  bool check_validity = false)
{
  CGAL_USE(check_validity);
  // std::cout << "n_to_m_flip::start" << std::endl;
  typedef typename C3t3::Triangulation  Tr;
  typedef typename C3t3::Vertex_handle  Vertex_handle;
  typedef typename C3t3::Cell_handle    Cell_handle;
  typedef typename C3t3::Facet          Facet;
  typedef typename Tr::Facet_circulator Facet_circulator;
  typedef typename Tr::Cell_circulator  Cell_circulator;

  Tr& tr = c3t3.triangulation();

  Vertex_handle vh0 = edge.first->vertex(edge.second);
  Vertex_handle vh1 = edge.first->vertex(edge.third);

  //This vertex will have its valence augmenting a lot,
  //TODO take the best one

  //TODO!!!! Check that the created edges do not exist!!!

  Facet_circulator facet_circulator = tr.incident_facets(edge);
  Facet_circulator done_facet_circulator = facet_circulator;
  bool look_for_vh_iterator = true;
  do
  {
    facet_circulator++;

    //Get the ids of the opposite vertices
    for (int i = 0; i < 3; ++i)
    {
      if (facet_circulator->first->vertex(indices(facet_circulator->second, i)) == vh)
        look_for_vh_iterator = false;
    }
  } while (facet_circulator != done_facet_circulator && look_for_vh_iterator);

  if (look_for_vh_iterator) {
    std::cout << "Vertex not an opposite of the edge!!" << std::endl;
    return NOT_FLIPPABLE;
  }

  Facet_circulator facet_done(facet_circulator);
  facet_done--;
  facet_circulator++;
  facet_circulator++;

  boost::optional<boost::container::small_vector<Cell_handle, 64>>& o_inc_vh = inc_cells[vh];
  if (o_inc_vh == boost::none)
  {
    boost::container::small_vector<Cell_handle, 64> inc_vec;
    tr.incident_cells(vh, std::back_inserter(inc_vec));
    o_inc_vh = inc_vec;
  }

  do
  {
    //Get the ids of the opposite vertices
    for (int i = 0; i < 3; ++i)
    {
      Vertex_handle curr_vertex = facet_circulator->first->vertex(
                                    indices(facet_circulator->second, i));
      if (curr_vertex != vh0  && curr_vertex != vh1)
      {
        if (is_edge_uv(vh, curr_vertex, boost::get(o_inc_vh)))
          return NOT_FLIPPABLE;
      }
    }
  } while (++facet_circulator != facet_done);


  boost::container::small_vector<Cell_handle, 20> to_remove;

  //Neighbors that will need to be updated after flip
  boost::unordered_set<Facet> neighbor_facets;

  //Facets that will be used to create new cells
  // i.e. all the facets opposite to vh1 and don't have vh
  std::vector<Facet> facets_for_new_cells;

  //Facets that will be used to update cells
  // i.e. all the facets opposite to vh0 will be set to vh :
  // facet.first->set_vertex( facet.second, vh )
  std::vector<Facet> facets_for_updated_cells;

  Cell_circulator cell_circulator = tr.incident_cells(edge);
  Cell_circulator done = cell_circulator;
  do
  {
    //Facets opposite to vh0
    Facet facet_vh0(cell_circulator, cell_circulator->index(vh0));
    neighbor_facets.insert(tr.mirror_facet(facet_vh0));

    //Facets opposite to vh1
    Facet facet_vh1(cell_circulator, cell_circulator->index(vh1));
    neighbor_facets.insert(tr.mirror_facet(facet_vh1));

    //Store it if it do not have vh
    if (cell_circulator->has_vertex(vh)) {
      to_remove.push_back(cell_circulator);
    }
    else
    {
      facets_for_new_cells.push_back(facet_vh1);
      facets_for_updated_cells.push_back(facet_vh0);
    }
    //
    //        if( ! is_well_oriented( cell_circulator ) )
    //            return WRONG;
  }
  while (++cell_circulator != done);

  //Check that the result will be valid
  for (const Facet& fi : facets_for_new_cells)
  {
    if ( !tr.is_infinite(fi.first)
         && !is_well_oriented(tr, vh, fi.first->vertex(indices(fi.second, 0)),
                              fi.first->vertex(indices(fi.second, 1)),
                              fi.first->vertex(indices(fi.second, 2))))
      return NOT_FLIPPABLE;
  }
  for (const Facet& fi : facets_for_updated_cells)
  {
    if ( !tr.is_infinite(fi.first)
         && !is_well_oriented(tr, vh, fi.first->vertex(indices(fi.second, 0)),
                              fi.first->vertex(indices(fi.second, 1)),
                              fi.first->vertex(indices(fi.second, 2))))
      return NOT_FLIPPABLE;
  }

  ///********************VALIDITY CHECK***************************/
  //double current_min_dh = DBL_MAX;

  //if (check_validity){

  //  pre_sliver_Removal_cells.clear();
  //  do{
  //    pre_sliver_Removal_cells.push_back(K::Tetrahedron_3(cell_circulator->vertex(0)->point(), cell_circulator->vertex(1)->point(),
  //      cell_circulator->vertex(2)->point(), cell_circulator->vertex(3)->point()));

  //    if (!tr.is_infinite(cell_circulator))
  //      current_min_dh = std::min(current_min_dh, min_dihedral_angle(cell_circulator));
  //  } while (++cell_circulator != done);

  //  pre_sliver_Removal_vertices.clear();
  //  pre_sliver_Removal_vertices.push_back(vh->point());

  //  previous_edges.clear();
  //  previous_edges.push_back(std::make_pair(vh0->point(), vh1->point()));
  //}
  ///*************************************************************/

  //Subdomain index?
  typename C3t3::Subdomain_index subdomain = to_remove[0]->subdomain_index();
  visitor.before_flip(to_remove[0]);

#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
  for (std::size_t i = 1; i < to_remove.size(); ++i)
    CGAL_assertion(subdomain == to_remove[i]->subdomain_index());
#endif

  std::vector<Cell_handle> cells_to_update;

  //Create new cells
  for (const Facet& fi : facets_for_new_cells)
  {
    Cell_handle new_cell = tr.tds().create_cell();

    for (int v = 0; v < 4; v++){
      new_cell->set_vertex(v, fi.first->vertex(v));
    }

    new_cell->set_vertex(fi.second, vh);

    c3t3.add_to_complex(new_cell, subdomain);
    visitor.after_flip(new_cell);
    cells_to_update.push_back(new_cell);
  }

  //Update_existing cells
  for (const Facet& fi : facets_for_updated_cells)
  {
    fi.first->set_vertex(fi.second, vh);
    cells_to_update.push_back(fi.first);
  }

  typedef CGAL::Triple<Vertex_handle, Vertex_handle, Vertex_handle> Facet_vvv;
  typedef boost::unordered_map<Facet_vvv, std::size_t> FaceMapIndex;

  FaceMapIndex facet_map_indices;
  std::vector<Facet> facets;

  for (const Facet& f : neighbor_facets)
  {
    Cell_handle ch = f.first;
    int v = f.second;

    Facet_vvv face = make_vertex_triple(ch->vertex(indices(v,0)),
                                        ch->vertex(indices(v,1)),
                                        ch->vertex(indices(v,2)));
    typename FaceMapIndex::iterator it = facet_map_indices.find(face);
    if (it == facet_map_indices.end())
    {
      facet_map_indices[face] = facets.size();
      facets.push_back(Facet(ch, v));
    }
  }

  //Update adjacencies and vertices cells
  for (Cell_handle ch : cells_to_update)
  {
    for (int v = 0; v < 4; v++)
    {
      Facet_vvv face = make_vertex_triple(ch->vertex(indices(v,0)),
                                          ch->vertex(indices(v,1)),
                                          ch->vertex(indices(v,2)));
      typename FaceMapIndex::iterator it = facet_map_indices.find(face);
      if (it == facet_map_indices.end())
      {
        facet_map_indices[face] = facets.size();
        facets.push_back(Facet(ch, v));
      }
      else
      {
        Facet facet = facets[it->second];

        //Update neighbor
        facet.first->set_neighbor(facet.second, ch);
        ch->set_neighbor(v, facet.first);
      }
      ch->vertex(v)->set_cell(ch);

      inc_cells[ch->vertex(v)] = boost::none;
      ch->reset_cache_validity();
    }
  }

  // Update c3t3
  update_c3t3_facets(c3t3, cells_to_update, neighbor_facets);

  //Remove cells
  for (Cell_handle ch : to_remove)
  {
    c3t3.remove_from_complex(ch);
    tr.tds().delete_cell(ch);
  }

  ///********************VALIDITY CHECK***************************/
  //if (check_validity){

  //  double new_min_dh = DBL_MAX;

  //  post_sliver_Removal_cells.clear();
  //  for (unsigned int i = 0; i < cells_to_update.size(); ++i){
  //    post_sliver_Removal_cells.push_back(cells_to_update[i]);

  //    if (!tr.is_infinite(cells_to_update[i]))
  //      new_min_dh = std::min(new_min_dh, min_dihedral_angle(cells_to_update[i]));
  //  }

  //  post_sliver_Removal_vertices.clear();
  //  for (unsigned int i = 0; i < vertices_around_edge.size(); ++i){
  //    post_sliver_Removal_vertices.push_back(vertices_around_edge[i]);
  //  }

  //  current_edges.clear();
  //  for (unsigned int i = 0; i < vertices_around_edge.size(); ++i){
  //    current_edges.push_back(std::make_pair(vertices_around_edge[i]->point(), vh->point()));
  //  }


  //  for (unsigned int i = 0; i < cells_to_update.size(); ++i){
  //    if (!tr.is_valid(cells_to_update[i]))
  //      return INVALID_CELL;

  //    for (int v = 0; v < 4; v++){
  //      if (!tr.is_valid(cells_to_update[i]->neighbor(v)))
  //        return INVALID_CELL;

  //      if (!tr.tds().is_valid(cells_to_update[i]->vertex(v)))
  //        return INVALID_VERTEX;

  //    }
  //  }

  //  if ((current_min_dh - new_min_dh) > 0.01){
  //    std::cout << pre_sliver_Removal_cells.size() << " to " << post_sliver_Removal_cells.size() << " flip not improving the quality: " <<
  //      current_min_dh << " to " << new_min_dh << std::endl;
  //    return INVALID_CELL;
  //  }

  //}
  ///***********************************************************/

  // std::cout << "n_to_m_flip::end with success" << std::endl;

  return VALID_FLIP;
}


template<typename C3t3, typename IncCellsVectorMap, typename Visitor>
Sliver_removal_result flip_n_to_m(typename C3t3::Edge& edge,
                                  C3t3& c3t3,
                                  const std::vector<typename C3t3::Vertex_handle>& boundary_vertices,
                                  const Flip_Criterion& criterion,
                                  IncCellsVectorMap& inc_cells,
                                  Visitor& visitor)
{
  typedef typename C3t3::Vertex_handle Vertex_handle;
  typedef typename C3t3::Triangulation::Cell_circulator Cell_circulator;
  typename C3t3::Triangulation& tr = c3t3.triangulation();

  Sliver_removal_result result = NOT_FLIPPABLE;

  typedef std::pair<Dihedral_angle_cosine, std::pair<Vertex_handle, int> > CosAngle_and_vertex;

  //std::cout << "n_to_m_flip " << boundary_vertices.size() << std::endl;
  if (criterion == MIN_ANGLE_BASED)
  {
    std::priority_queue<CosAngle_and_vertex,
                        std::vector<CosAngle_and_vertex>,
                        std::greater<CosAngle_and_vertex>
                      > candidates;

    Cell_circulator circ = c3t3.triangulation().incident_cells(edge);
    Cell_circulator done = circ;

    Dihedral_angle_cosine curr_max_cosdh = max_cos_dihedral_angle(tr, circ++);
    while (circ != done)
    {
      curr_max_cosdh = (std::max)(curr_max_cosdh, max_cos_dihedral_angle(tr, circ++));
    }
    if (boundary_vertices.size() == 2)
      find_best_flip_to_improve_dh(c3t3, edge, boundary_vertices[0], boundary_vertices[1],
                                   candidates, curr_max_cosdh);
    else
      find_best_flip_to_improve_dh(c3t3, edge, candidates, curr_max_cosdh,
                                   inc_cells);

    bool flip_performed = false;
    while (!candidates.empty() && !flip_performed)
    {
      CosAngle_and_vertex curr_cost_vpair = candidates.top();
      candidates.pop();

//      std::cout << "\tcurrent   cos = " << curr_max_cosdh
//        << "\t angle = " << std::acos(curr_max_cosdh) * 180./CGAL_PI << std::endl;
//      std::cout << "\tcandidate cos = " << curr_cost_vpair.first
//        << "\t angle = " << std::acos(curr_cost_vpair.first) * 180./CGAL_PI << std::endl;
//      std::cout << std::endl;

      if (curr_max_cosdh <= curr_cost_vpair.first)
        return NO_BEST_CONFIGURATION;

      result = flip_n_to_m(c3t3, edge, curr_cost_vpair.second.first, inc_cells, visitor);

      if (result != NOT_FLIPPABLE)
        flip_performed = true;
    }
  }

  return result;
}

template<typename C3t3, typename IncCellsVectorMap, typename Visitor>
Sliver_removal_result find_best_flip(typename C3t3::Edge& edge,
                                     C3t3& c3t3,
                                     const Flip_Criterion& criterion,
                                     IncCellsVectorMap& inc_cells,
                                     Visitor& visitor)
{
  typedef typename C3t3::Triangulation        Tr;
  typedef typename C3t3::Vertex_handle        Vertex_handle;
  typedef typename Tr::Facet_circulator       Facet_circulator;

  Tr& tr = c3t3.triangulation();

  const Vertex_handle v0 = edge.first->vertex(edge.second);
  const Vertex_handle v1 = edge.first->vertex(edge.third);

  Facet_circulator circ = tr.incident_facets(edge);
  Facet_circulator done = circ;

  //Identify the vertices around this edge
  boost::unordered_set<Vertex_handle> vertices_around_edge;
  bool boundary_edge = false;
  bool hull_edge = false;

  boost::unordered_set<Vertex_handle> boundary_vertices;
//  boost::unordered_set<Vertex_handle> hull_vertices;
  do
  {
    //Get the ids of the opposite vertices
    for (int i = 0; i < 3; ++i)
    {
      Vertex_handle vi = circ->first->vertex(indices(circ->second, i));
      if (vi != v0 && vi != v1)
      {
        vertices_around_edge.insert(vi);

        if ( circ->first->subdomain_index()
             != circ->first->neighbor(circ->second)->subdomain_index())
        {
          boundary_edge = true;
          boundary_vertices.insert(vi);
        }

        if ( tr.is_infinite(circ->first)
             != tr.is_infinite(circ->first->neighbor(circ->second)))
        {
          hull_edge = true;
          //hull_vertices.insert(vi);
        }
      }
    }
  }
  while (++circ != done);


  //Check if not feature edge
  if (boundary_vertices.size() > 2)
    return NOT_FLIPPABLE;

  // perform flip when possible
  Sliver_removal_result res = NOT_FLIPPABLE;
  if (vertices_around_edge.size() == 3)
  {
    if (!boundary_edge && !hull_edge)
    {
      std::vector<Vertex_handle> vertices;
      vertices.insert(vertices.end(), vertices_around_edge.begin(), vertices_around_edge.end());
      res = flip_3_to_2(edge, c3t3, vertices, criterion, inc_cells);
    }
  }
  else
  {
    //TODO fix for hull edges
    // if( hull_edge )
    //    return n_to_m_flip( edge, hull_vertices, flip_criterion, check_validity );
    if (!hull_edge)
    {
      std::vector<Vertex_handle> vertices;
      vertices.insert(vertices.end(), boundary_vertices.begin(), boundary_vertices.end());
      res = flip_n_to_m(edge, c3t3, vertices, criterion, inc_cells, visitor);
      //return n_to_m_flip(edge, boundary_vertices, flip_criterion);
    }
  }

  return res;
}


template<typename VertexPair, typename C3t3, typename Visitor>
std::size_t flip_all_edges(const std::vector<VertexPair>& edges,
                           C3t3& c3t3,
                           const Flip_Criterion& criterion,
                           Visitor& visitor)
{
  typedef typename C3t3::Triangulation Tr;
  typedef typename Tr::Cell_handle   Cell_handle;
  typedef typename Tr::Vertex_handle Vertex_handle;
  typedef typename Tr::Edge          Edge;

  Tr& tr = c3t3.triangulation();

  std::unordered_map<Vertex_handle,
    boost::optional<boost::container::small_vector<Cell_handle, 64> > > inc_cells;

  std::size_t count = 0;
  for (const VertexPair& vp : edges)
  {
    boost::optional<boost::container::small_vector<Cell_handle, 64>>&
      o_inc_vh = inc_cells[vp.first];
    if (o_inc_vh == boost::none)
    {
      boost::container::small_vector<Cell_handle, 64> inc_vec;
      tr.incident_cells(vp.first, std::back_inserter(inc_vec));
      o_inc_vh = inc_vec;
    }

    Cell_handle ch;
    int i0, i1;
    if (is_edge_uv(vp.first, vp.second, boost::get(o_inc_vh), ch, i0, i1))
    {
      Edge edge(ch, i0, i1);

      Sliver_removal_result res = find_best_flip(edge, c3t3, criterion, inc_cells, visitor);
      if (res == INVALID_CELL || res == INVALID_VERTEX || res == INVALID_ORIENTATION)
      {
        std::cout << "FLIP PROBLEM!!!!" << std::endl;
        return count;
      }
      if (res == VALID_FLIP)
      {
        ++count;
#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE_PROGRESS
        std::cout << "\rFlip... (";
        std::cout << count << " flips)";
        std::cout.flush();
#endif
      }
    }
  }

  for(Cell_handle c : c3t3.triangulation().finite_cell_handles())
    c->reset_cache_validity();

  return count;
}

template<typename C3T3, typename CellSelector, typename Visitor>
void flip_edges(C3T3& c3t3,
                const bool protect_boundaries,
                CellSelector cell_selector,
                Visitor& visitor)
{
  CGAL_USE(protect_boundaries);
  typedef typename C3T3::Triangulation       T3;
  typedef typename T3::Vertex_handle         Vertex_handle;
  typedef typename std::pair<Vertex_handle, Vertex_handle> Edge_vv;

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
  std::cout << "Flip edges...";
  std::cout.flush();
  std::size_t nb_flips = 0;
#endif

  //const Flip_Criterion criterion = VALENCE_MIN_DH_BASED;

  //collect long edges

  //compute vertices normals map?

  // typedef typename C3T3::Surface_patch_index Surface_patch_index;
  // typedef boost::unordered_map<Surface_patch_index, unsigned int> Spi_map;
  //if (!protect_boundaries)
  //{
  //  std::cout << "\tBoundary flips" << std::endl;
  //  //Boundary flip
  //  std::vector<Spi_map> boundary_vertices_valences;
  //  std::vector<Edge_vv> boundary_edges;

  //  collectBoundaryEdges(boundary_edges);

  //  computeVerticesValences(boundary_vertices_valences);

  //  if (criterion == VALENCE_BASED)
  //    flipBoundaryEdges(boundary_edges, boundary_vertices_valences, VALENCE_BASED);
  //  else
  //    flipBoundaryEdges(boundary_edges, boundary_vertices_valences, MIN_ANGLE_BASED);
  //}

  std::vector<Edge_vv> inside_edges;
  get_internal_edges(c3t3,
                     cell_selector,
                     std::back_inserter(inside_edges));

  //if (criterion == VALENCE_BASED)
  //  flip_inside_edges(inside_edges);
  //else
  //{
#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
  nb_flips =
#endif
    flip_all_edges(inside_edges, c3t3, MIN_ANGLE_BASED, visitor);
  //}

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
  std::cout << " done (" << nb_flips << " flips)." << std::endl;
#endif
}

}//namespace internal
}//namespace Tetrahedral_remeshing
}//namespace CGAL

#endif // CGAL_INTERNAL_FLIP_EDGES_H
