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

#include <CGAL/Tetrahedral_remeshing/internal/Elementary_operation.h>
#include <CGAL/Tetrahedral_remeshing/internal/tetrahedral_remeshing_helpers.h>

#include <boost/container/small_vector.hpp>
#include <boost/functional/hash.hpp>

#ifdef CGAL_LINKED_WITH_TBB
#include <tbb/concurrent_unordered_map.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#endif

#include <unordered_map>
#include <unordered_set>
#include <algorithm>
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

template<typename C3t3, typename IncCellsVectorMap, typename CellSelector>
Sliver_removal_result flip_3_to_2(typename C3t3::Edge& edge,
                                  C3t3& c3t3,
                                  const boost::container::small_vector<typename C3t3::Vertex_handle, 16>& vertices_around_edge,
                                  const Flip_Criterion& criterion,
                                  IncCellsVectorMap& inc_cells,
                                  CellSelector& cell_selector)
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
  // Only the answer is used here, so the threadsafe walk -- which returns no
  // cell -- serves this site exactly.
  if (is_facet_tagged(tr, vertices_around_edge[0], vertices_around_edge[1],
                      vertices_around_edge[2]))
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
  typedef std::unordered_map<Facet_vvv, std::size_t> FaceMapIndex;
  std::unordered_set<Facet, boost::hash<Facet>> outer_mirror_facets;

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

      inc_cells[ch->vertex(v)].clear();
      ch->reset_cache_validity();
    }
  }

  // Update c3t3
  update_c3t3_facets(c3t3, cells_to_update, outer_mirror_facets);

  treat_before_delete(cell_to_remove, cell_selector, c3t3);
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
              if (is_edge_tagged(tr, curr, vh))
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

  for (std::size_t i = 0; i < opposite_vertices.size(); ++i)
  {
    Vertex_handle vh = opposite_vertices[i];
    bool keep = true;

    boost::container::small_vector<Facet, 60> facets;
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
      if (!tr.is_infinite(fi.first) && c3t3.is_in_complex(fi.first))
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

        // the worst angle of the flip only ever grows from here, so once it has
        // reached the one the edge already has, this vertex cannot be kept -
        // unless the cells are inverted, where it is kept whatever it measures
        if (is_sliver_well_oriented && !(max_flip_cos_dh < curr_max_cos_dh))
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

  //Collect the vertex opposite to the edge in each facet around it, in
  //circulation order. The chord test below used to re-circulate the facets
  //around the edge once per such vertex to enumerate the other apices; those
  //are the same vertices collected here, so the tests can be indexed on the
  //ring instead. is_edge_uv is read-only, so stopping at the first chord
  //found yields the same verdict as running the ring to its end.
  boost::container::small_vector<Vertex_handle, 32> ring_apices;
  do
  {
    //Get the id of the opposite vertex
    for (int i = 0; i < 3; ++i)
    {
      Vertex_handle curr_vertex = curr_fcirc->first->vertex(
                                    indices(curr_fcirc->second, i));
      if (curr_vertex != vh0 && curr_vertex != vh1)
      {
        ring_apices.push_back(curr_vertex);
        break;
      }
    }
  }
  while (++curr_fcirc != curr_fdone);

  //Only keep the possible flips
  std::vector<Vertex_handle> opposite_vertices;
  int nb_cells_around_edge = 0;
  const int n_apices = static_cast<int>(ring_apices.size());

  //Gather the star of every finite apex first. An infinite apex has no star
  //and is never a candidate, but it can still be the far end of a chord, so
  //it keeps its place on the ring and is marked by a null star here.
  using Star = boost::container::small_vector<Cell_handle, 64>;
  boost::container::small_vector<Star*, 32> apex_star(n_apices, nullptr);
  for (int p = 0; p < n_apices; ++p)
  {
    const Vertex_handle vh = ring_apices[p];

    if(tr.is_infinite(vh))
      continue;

    Star& o_inc_vh = inc_cells[vh];
    if (o_inc_vh.empty())
      incident_cells_tagged(tr, vh, std::back_inserter(o_inc_vh));

    apex_star[p] = &o_inc_vh;
    nb_cells_around_edge++;
  }

  //a chord is an edge joining two apices that are not neighbors on the ring
  //(positions p-1 and p+1).
  //
  //Each unordered pair is settled ONCE here. The test used to sit inside the
  //loop above, running j from p+2 to p+n_apices-2 for every apex in turn, so
  //every pair was asked about twice -- once from each of its two ends -- and
  //each ask walks one apex's whole star looking for a cell that holds the
  //other. A chord is found for only 3.6% of the pairs, so the second ask
  //almost always repeats the first ask's full walk.
  //
  //The relation is symmetric and so is the exclusion it causes: a chord
  //between p and j disqualifies both of them. Asking once from p and marking
  //both ends therefore leaves `opposite_vertices` exactly as it was, with the
  //same members in the same ascending order of p. A pair whose finite ends
  //are both already disqualified is skipped, which is what the old `break`
  //bought, and the old `break` never suppressed a pair the other end still
  //needed: that end asked about it in its own turn.
  boost::container::small_vector<char, 32> chorded(n_apices, 0);
  for (int p = 0; p < n_apices; ++p)
  {
    for (int j = p + 2; j < n_apices && j <= p + n_apices - 2; ++j)
    {
      const bool ask_from_p = (apex_star[p] != nullptr) && !chorded[p];
      const bool ask_from_j = (apex_star[j] != nullptr) && !chorded[j];
      if (!ask_from_p && !ask_from_j)
        continue;

      const int u = ask_from_p ? p : j;
      const int w = ask_from_p ? j : p;
      if (is_edge_uv(ring_apices[u], ring_apices[w], *apex_star[u]))
      {
        if (apex_star[p] != nullptr) chorded[p] = 1;
        if (apex_star[j] != nullptr) chorded[j] = 1;
      }
    }
  }

  for (int p = 0; p < n_apices; ++p)
  {
    if (apex_star[p] != nullptr && !chorded[p])
      opposite_vertices.push_back(ring_apices[p]);
  }

  if (nb_cells_around_edge < 4)
    return;

  //Facets that will be used to create new cells
  //    i.e. all the facets opposite to vh1 and don't have vh
  //Facets that will be used to update cells
  //    i.e. all the facets opposite to vh0 will be set to vh:
  //    facet.first->set_vertex( facet.second, vh )

  // The facets an apex is judged on are produced and judged in ONE pass.
  //
  // They used to be collected into a `small_vector<Facet, 60>` by a full turn
  // of the ring, and only then evaluated -- and the evaluation abandons the
  // apex on its FIRST facet most of the time, because the three tests below
  // (inverted cell, a worst angle of one, an angle no better than the edge
  // already has) are the common case, not the exception. Every facet the turn
  // collected past that point was collected for nothing.
  //
  // Same facets, in the same order -- the ring is circulated from the same
  // cell for every apex, and a cell still yields the facet opposite `vh1`
  // before the one opposite `vh0` -- and the same three exits, so `keep` and
  // `max_flip_cos_dh` are what they were. `max_flip_cos_dh` is a max over the
  // facets, and it is only READ when no exit was taken, i.e. when the fold ran
  // over the whole ring either way.
  for (Vertex_handle vh : opposite_vertices)
  {
    bool keep = true;
    Dihedral_angle_cosine max_flip_cos_dh(CGAL::NEGATIVE, 1., 1.);

    Cell_circulator cell_circulator = tr.incident_cells(edge);
    Cell_circulator done = cell_circulator;
    do
    {
      //Cells that do not have vh are the ones the flip rewrites
      if (!cell_circulator->has_vertex(vh))
      {
        //Facets opposite to vh1, then opposite to vh0
        const Facet ring_facets[2]
          = { Facet(cell_circulator, cell_circulator->index(vh1)),
              Facet(cell_circulator, cell_circulator->index(vh0)) };

        for (const Facet& fi : ring_facets)
        {
          if (tr.is_infinite(fi.first))
            continue;

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

          // the worst angle of the flip only ever grows from here, so once it has
          // reached the one the edge already has, this vertex cannot be kept -
          // unless the cells are inverted, where it is kept whatever it measures
          if (is_sliver_well_oriented && !(max_flip_cos_dh < curr_max_cosdh))
          {
            keep = false;
            break;
          }
        }

        if (!keep)
          break;
      }
    }
    while (++cell_circulator != done);

    if (keep && (max_flip_cos_dh < curr_max_cosdh || !is_sliver_well_oriented))
    {
      //std::cout << "vh " << vh->info() <<" old " << curr_max_cosdh << " min " << min_flip_tan_dh << std::endl;
      candidates.push(std::make_pair(max_flip_cos_dh, std::make_pair(vh, e_id)));
    }
  }
}

template<typename C3t3,
         typename IncCellsVectorMap,
         typename CellSelector,
         typename Visitor>
Sliver_removal_result flip_n_to_m(C3t3& c3t3,
                                  typename C3t3::Edge& edge,
                                  typename C3t3::Vertex_handle vh,
                                  IncCellsVectorMap& inc_cells,
                                  CellSelector& cell_selector,
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

  const Vertex_handle vh0 = edge.first->vertex(edge.second);
  const Vertex_handle vh1 = edge.first->vertex(edge.third);

  //This vertex will have its valence augmenting a lot,
  //TODO take the best one

  //TODO!!!! Check that the created edges do not exist!!!

  boost::container::small_vector<Facet, 2> facets_in_complex;

  Facet_circulator facet_circulator = tr.incident_facets(edge);
  Facet_circulator done_facet_circulator = facet_circulator;
  bool look_for_vh_iterator = true;
  do
  {
    if (c3t3.is_in_complex(*facet_circulator))
    {
      facets_in_complex.push_back(*facet_circulator);
    }

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

  boost::container::small_vector<Cell_handle, 64>& o_inc_vh = inc_cells[vh];
  if (o_inc_vh.empty())
    incident_cells_tagged(tr, vh, std::back_inserter(o_inc_vh));

  do
  {
    //Get the ids of the opposite vertices
    for (int i = 0; i < 3; ++i)
    {
      Vertex_handle curr_vertex = facet_circulator->first->vertex(
                                    indices(facet_circulator->second, i));
      if (curr_vertex != vh0  && curr_vertex != vh1)
      {
        if (is_edge_uv(vh, curr_vertex, o_inc_vh))
          return NOT_FLIPPABLE;
      }
    }
  } while (++facet_circulator != facet_done);


  boost::container::small_vector<Cell_handle, 20> to_remove;

  //Neighbors that will need to be updated after flip
  std::unordered_set<Facet, boost::hash<Facet>> neighbor_facets;

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

  //Surface
  for (const Facet& f : facets_in_complex)
    c3t3.remove_from_complex(f);

  //Subdomain index
  typedef typename C3t3::Subdomain_index Subdomain_index;
  const Subdomain_index subdomain = to_remove[0]->subdomain_index();
  bool selected = get(cell_selector, to_remove[0]);
  visitor.before_flip(to_remove[0]);

  std::vector<Cell_handle> cells_to_update;

  //Create new cells
  for (const Facet& fi : facets_for_new_cells)
  {
    Cell_handle new_cell = tr.tds().create_cell();

    for (int v = 0; v < 4; v++){
      new_cell->set_vertex(v, fi.first->vertex(v));
    }

    new_cell->set_vertex(fi.second, vh);

    treat_new_cell(new_cell, subdomain, cell_selector, selected, c3t3);

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
  typedef std::unordered_map<Facet_vvv, std::size_t> FaceMapIndex;

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

      inc_cells[ch->vertex(v)].clear();
    }
    ch->reset_cache_validity();
  }

  // Update c3t3
  update_c3t3_facets(c3t3, cells_to_update, neighbor_facets);

  //Remove cells
  for (Cell_handle ch : to_remove)
  {
    treat_before_delete(ch, cell_selector, c3t3);
    ch->reset_cache_validity();
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


template<typename C3t3, typename IncCellsVectorMap, typename CellSelector, typename Visitor>
Sliver_removal_result flip_n_to_m(typename C3t3::Edge& edge,
                                  C3t3& c3t3,
                                  const std::vector<typename C3t3::Vertex_handle>& boundary_vertices,
                                  const Flip_Criterion& criterion,
                                  IncCellsVectorMap& inc_cells,
                                  CellSelector& cell_selector,
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
    do
    {
      curr_max_cosdh = (std::max)(curr_max_cosdh, max_cos_dihedral_angle(tr, circ));
    } while (++circ != done);

    if (boundary_vertices.size() == 2)
      find_best_flip_to_improve_dh(c3t3, edge, boundary_vertices[0], boundary_vertices[1],
                                   candidates, curr_max_cosdh);
    else
      find_best_flip_to_improve_dh(c3t3, edge, candidates, curr_max_cosdh, inc_cells);

    bool flip_performed = false;
    while (!candidates.empty() && !flip_performed)
    {
      CosAngle_and_vertex curr_cost_vpair = candidates.top();
      candidates.pop();

//      std::cout << "\tcurrent   cos = " << curr_max_cosdh.value()
//        << "\t angle = " << std::acos(curr_max_cosdh.value()) * 180./CGAL_PI << std::endl;
//      std::cout << "\tcandidate cos = " << curr_cost_vpair.first.value()
//        << "\t angle = " << std::acos(curr_cost_vpair.first.value()) * 180./CGAL_PI << std::endl;
//      std::cout << std::endl;

      if (curr_max_cosdh <= curr_cost_vpair.first)
        return NO_BEST_CONFIGURATION;

      result = flip_n_to_m(c3t3, edge, curr_cost_vpair.second.first, inc_cells,
                           cell_selector, visitor);

      if (result != NOT_FLIPPABLE)
        flip_performed = true;
    }
  }

  return result;
}

template<typename C3t3, typename IncCellsVectorMap, typename CellSelector, typename Visitor>
Sliver_removal_result find_best_flip(typename C3t3::Edge& edge,
                                     C3t3& c3t3,
                                     const Flip_Criterion& criterion,
                                     IncCellsVectorMap& inc_cells,
                                     CellSelector& cell_selector,
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

  //Identify the vertices around this edge. The ring of apices around an edge
  //holds a handful of distinct vertices, so they are kept inline and scanned
  //rather than hashed. flip_3_to_2 does not depend on their order : it picks
  //vh2/vh3 by testing each vertex against ch0/ch1 individually, and is_facet
  //is symmetric in its three vertices.
  boost::container::small_vector<Vertex_handle, 16> vertices_around_edge;
  bool boundary_edge = false;
  bool hull_edge = false;

  std::unordered_set<Vertex_handle> boundary_vertices;
//  std::unordered_set<Vertex_handle> hull_vertices;
  do
  {
    //Get the ids of the opposite vertices
    for (int i = 0; i < 3; ++i)
    {
      Vertex_handle vi = circ->first->vertex(indices(circ->second, i));
      if (vi != v0 && vi != v1)
      {
        if (std::find(vertices_around_edge.begin(), vertices_around_edge.end(), vi)
            == vertices_around_edge.end())
          vertices_around_edge.push_back(vi);

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
      res = flip_3_to_2(edge, c3t3, vertices_around_edge, criterion, inc_cells, cell_selector);
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
      res = flip_n_to_m(edge, c3t3, vertices, criterion, inc_cells, cell_selector, visitor);
      //return n_to_m_flip(edge, boundary_vertices, flip_criterion);
    }
  }

  return res;
}


template<typename VertexPair, typename C3t3,
         typename IncidentCellsVectorMap, typename CellSelector, typename Visitor>
std::size_t flip_all_edges(const std::vector<VertexPair>& edges,
                           C3t3& c3t3,
                           IncidentCellsVectorMap& inc_cells,
                           const Flip_Criterion& criterion,
                           CellSelector& cell_selector,
                           Visitor& visitor)
{
  typedef typename C3t3::Triangulation Tr;
  typedef typename Tr::Cell_handle   Cell_handle;
  typedef typename Tr::Edge          Edge;

  Tr& tr = c3t3.triangulation();

  std::size_t count = 0;
  for (const VertexPair& vp : edges)
  {
    boost::container::small_vector<Cell_handle, 64>& o_inc_vh = inc_cells[vp.first];
    if (o_inc_vh.empty())
      tr.incident_cells(vp.first, std::back_inserter(o_inc_vh));

    Cell_handle ch;
    int i0, i1;
    if (is_edge_uv(vp.first, vp.second, o_inc_vh, ch, i0, i1))
    {
      Edge edge(ch, i0, i1);

      Sliver_removal_result res
        = find_best_flip(edge, c3t3, criterion, inc_cells, cell_selector, visitor);
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

  return count;
}

/**
* Per boundary vertex, the distinct subdomain indices of the cells around it.
*
* ONLY THE COUNT is ever read -- `for_each_boundary_valence_patch()` asks
* whether a vertex touches more than one subdomain -- so the per-vertex set is
* a short inline array scanned linearly, not a hash set. A boundary vertex
* touches one or two subdomains almost everywhere, and a `std::unordered_set`
* per vertex paid a node allocation to hold that.
*
* The container is concurrent under `Parallel_tag` because the parallel build
* fills it one vertex per thread.
*/
template<typename C3T3>
using Vertex_subdomains_map_t = Concurrency_selected_container_t<
  typename C3T3::Triangulation::Concurrency_tag,
  boost::unordered_map<typename C3T3::Vertex_handle,
                       boost::container::small_vector<typename C3T3::Subdomain_index, 4> >,
#ifdef CGAL_LINKED_WITH_TBB
  tbb::concurrent_unordered_map<typename C3T3::Vertex_handle,
                                boost::container::small_vector<typename C3T3::Subdomain_index, 4>,
                                boost::hash<typename C3T3::Vertex_handle> >
#else
  boost::unordered_map<typename C3T3::Vertex_handle,
                       boost::container::small_vector<typename C3T3::Subdomain_index, 4> >
#endif
  >;

// Records `si` for `v` if it is not there already.
template<typename Subdomains, typename Subdomain_index>
void record_subdomain(Subdomains& subdomains, const Subdomain_index& si)
{
  for (const Subdomain_index& s : subdomains)
    if (s == si)
      return;
  subdomains.push_back(si);
}

template<typename C3t3>
void collect_subdomains_on_boundary(const C3t3& c3t3,
  Vertex_subdomains_map_t<C3t3>& vertices_subdomain_indices)
{
  for (auto c : c3t3.triangulation().all_cell_handles())
  {
    for (auto v : c3t3.triangulation().vertices(c))
    {
      const int dim = v->in_dimension();
      if(dim >= 0 && dim < 3)
        record_subdomain(vertices_subdomain_indices[v], c->subdomain_index());
    }
  }
}

#ifdef CGAL_LINKED_WITH_TBB
/**
* The two whole-mesh scans the boundary flip's preprocessing opens with: the
* boundary edges, and the subdomain indices incident to each boundary vertex.
* Both are read-only walks of the same container in front of a phase that then
* runs on every thread, so on the parallel path they run over ONE shared cell
* snapshot -- walked once here where the serial arm walks it twice, once for
* the edge iterator and once for `all_cell_handles()`.
*
* The valence reduction that follows them stays serial; it is a per-boundary-edge
* loop, not a whole-mesh one.
*/
template<typename C3T3, typename CellSelector>
void collect_boundary_edges_and_subdomains_parallel(
  const C3T3& c3t3,
  const CellSelector& cell_selector,
  std::vector<typename C3T3::Edge>& boundary_edges,
  Vertex_subdomains_map_t<C3T3>& vertices_subdomain_indices)
{
  using Edge = typename C3T3::Edge;
  const typename C3T3::Triangulation& tr = c3t3.triangulation();

  using Cell_handle = typename C3T3::Cell_handle;
  using Vertex_handle = typename C3T3::Vertex_handle;
  using Subdomain_index = typename C3T3::Subdomain_index;

  const std::vector<Cell_handle> cells = gather_all_cells(tr);

  // Each thread unions into its own map and the maps are unioned afterwards:
  // the result is a SET per vertex, so it does not depend on the order the
  // cells were visited in, and only its size is ever read.
  //
  // Building it BY VERTEX instead -- one owner per vertex, no merge at all --
  // was tried and is slower: it trades this cell scan for a star walk per
  // boundary vertex, and the pointer chasing costs more than the merge it
  // removes (-1.116% instructions but +2.82% wall on 1146193_cdt_0.5). See
  // `rejected/bsubmap_pervertex`.
  using Vertex_subdomains_map
    = boost::unordered_map<Vertex_handle,
                           boost::container::small_vector<Subdomain_index, 4> >;
  tbb::enumerable_thread_specific<Vertex_subdomains_map> tl_vsi;
  tbb::parallel_for(tbb::blocked_range<std::size_t>(0, cells.size()),
    [&](const tbb::blocked_range<std::size_t>& range)
    {
      Vertex_subdomains_map& local = tl_vsi.local();
      for (std::size_t ci = range.begin(); ci != range.end(); ++ci)
      {
        const Cell_handle c = cells[ci];
        for (auto v : tr.vertices(c))
        {
          const int dim = v->in_dimension();
          if (dim >= 0 && dim < 3)
            record_subdomain(local[v], c->subdomain_index());
        }
      }
    });

  for (const Vertex_subdomains_map& local : tl_vsi)
    for (const auto& vsi : local)
    {
      auto& s = vertices_subdomain_indices[vsi.first];
      for (const Subdomain_index& si : vsi.second)
        record_subdomain(s, si);
    }

  // Same per-edge test and the same edge set as the serial walk; the scan
  // applies the smallest-incident-cell ownership rule the finite edge
  // iterator applies, in the same cell and slot order, so `boundary_edges`
  // comes back in the order the serial walk fills it -- and the valence loop
  // and the candidate list below it are built from that order.
  boundary_edges = parallel_collect_from_finite_edges<Edge>(tr, cells,
    [&](const Edge& e, std::vector<Edge>& out)
    {
      if (is_boundary(c3t3, e, cell_selector))
        out.push_back(e);
    });
}

/**
* The per-boundary-edge valence reduction, split in two so that the parallel
* half needs no per-thread map and no reduction afterwards.
*
* The cost of the serial loop below is the facet circulator and the
* `is_in_complex()` test it runs per incident facet; the map increments
* themselves are two hash lookups per qualifying facet. `bflip_valence_parallel`
* (rejected 2026-09-16) reduced per-thread maps and paid for the merge with
* instructions the serial loop never spends -- +0.7 to +2.4% on the small
* configs. Here the parallel pass only walks the circulators, writing the
* patches it would have counted into a slot it owns, and the map is then filled
* by a serial pass that does exactly the increments the serial loop does, in
* the same order. Nothing is merged, and nothing is allocated per thread.
*
* Counting is commutative, so edge order never mattered for the result; it is
* kept anyway because the serial pass walks `boundary_edges` as it is.
*/
template<typename C3T3, typename PatchOutput>
void for_each_boundary_valence_patch(
  const C3T3& c3t3,
  const typename C3T3::Edge& e,
  const Vertex_subdomains_map_t<C3T3>& vertices_subdomain_indices,
  PatchOutput out)
{
  using Surface_patch_index = typename C3T3::Surface_patch_index;
  using Facet_circulator = typename C3T3::Triangulation::Facet_circulator;

  const typename C3T3::Triangulation& tr = c3t3.triangulation();
  const auto v0 = e.first->vertex(e.second);
  const auto v1 = e.first->vertex(e.third);

  // `find()` where the serial loop uses `operator[]`: the entry it would
  // insert is empty, so a missing vertex has size 0 either way, and reading
  // through `find()` leaves the map untouched by a parallel pass.
  const auto it0 = vertices_subdomain_indices.find(v0);
  const auto it1 = vertices_subdomain_indices.find(v1);
  const std::size_t n0 = (it0 == vertices_subdomain_indices.end()) ? 0 : it0->second.size();
  const std::size_t n1 = (it1 == vertices_subdomain_indices.end()) ? 0 : it1->second.size();

  Facet_circulator facet_circulator = tr.incident_facets(e);
  Facet_circulator done(facet_circulator);

  //In case of feature edge
  if (n0 > 2 && n1 > 2)
  {
    do
    {
      if (c3t3.is_in_complex(*facet_circulator))
        out(c3t3.surface_patch_index(*facet_circulator));
    } while (++facet_circulator != done);
  }
  //Normal surface edge, or non-manifold edge on dangling facet
  else
  {
    Surface_patch_index first_patch = Surface_patch_index();
    do
    {
      if (c3t3.is_in_complex(*facet_circulator))
      {
        const Surface_patch_index surfi = c3t3.surface_patch_index(*facet_circulator);
        if (first_patch == Surface_patch_index())
          first_patch = surfi;
        else if (first_patch == surfi)
          continue;

        out(surfi);
      }
    } while (++facet_circulator != done);
  }
}

template<typename C3T3, typename BoundaryValencesMap>
void compute_boundary_vertices_valences_parallel(
  const C3T3& c3t3,
  const std::vector<typename C3T3::Edge>& boundary_edges,
  BoundaryValencesMap& boundary_vertices_valences,
  const Vertex_subdomains_map_t<C3T3>& vertices_subdomain_indices)
{
  using Surface_patch_index = typename C3T3::Surface_patch_index;

  // One slot per edge, written by the one thread that owns that edge. Four
  // patches cover a manifold surface edge (two facets) with room to spare; an
  // edge that needs more sets `overflow` and is walked again by the serial
  // pass, so the buffer's size is a performance choice and never a
  // correctness one.
  static constexpr unsigned int inline_patches = 4;
  struct Edge_patches
  {
    std::array<Surface_patch_index, inline_patches> patches;
    unsigned char n = 0;
    bool overflow = false;
  };

  std::vector<Edge_patches> per_edge(boundary_edges.size());

  tbb::parallel_for(tbb::blocked_range<std::size_t>(0, boundary_edges.size()),
    [&](const tbb::blocked_range<std::size_t>& range)
    {
      for (std::size_t ei = range.begin(); ei != range.end(); ++ei)
      {
        Edge_patches& slot = per_edge[ei];
        for_each_boundary_valence_patch(c3t3, boundary_edges[ei],
          vertices_subdomain_indices,
          [&slot](const Surface_patch_index& surfi)
          {
            if (slot.n < inline_patches)
              slot.patches[slot.n++] = surfi;
            else
              slot.overflow = true;
          });
      }
    });

  for (std::size_t ei = 0; ei < boundary_edges.size(); ++ei)
  {
    const typename C3T3::Edge& e = boundary_edges[ei];
    const auto v0 = e.first->vertex(e.second);
    const auto v1 = e.first->vertex(e.third);
    const Edge_patches& slot = per_edge[ei];

    if (slot.overflow)
    {
      for_each_boundary_valence_patch(c3t3, e, vertices_subdomain_indices,
        [&](const Surface_patch_index& surfi)
        {
          boundary_vertices_valences[v0][surfi]++;
          boundary_vertices_valences[v1][surfi]++;
        });
      continue;
    }

    // One outer lookup per endpoint, not one per patch. This map is the
    // concurrent one -- the routine runs only on the parallel path -- and its
    // references stay valid across insertions.
    if (slot.n != 0)
    {
      auto& m0 = boundary_vertices_valences[v0];
      auto& m1 = boundary_vertices_valences[v1];
      for (unsigned char k = 0; k < slot.n; ++k)
      {
        m0[slot.patches[k]]++;
        m1[slot.patches[k]]++;
      }
    }
  }
}

#endif // CGAL_LINKED_WITH_TBB

template<typename C3T3, typename CellSelector, typename BoundaryValencesMap>
void collectBoundaryEdgesAndComputeVerticesValences(
  const C3T3& c3t3,
  const CellSelector& cell_selector,
  std::vector<typename C3T3::Edge>& boundary_edges,
  BoundaryValencesMap& boundary_vertices_valences,
  Vertex_subdomains_map_t<C3T3>& vertices_subdomain_indices)
{
  typedef typename C3T3::Surface_patch_index Surface_patch_index;
  typedef typename C3T3::Vertex_handle       Vertex_handle;
  typedef typename C3T3::Edge                Edge;
  typedef typename C3T3::Triangulation::Facet_circulator Facet_circulator;

  const typename C3T3::Triangulation& tr = c3t3.triangulation();

  boundary_edges.clear();
  boundary_vertices_valences.clear();

#ifdef CGAL_LINKED_WITH_TBB
  if constexpr (is_parallel_triangulation<typename C3T3::Triangulation>())
  {
    collect_boundary_edges_and_subdomains_parallel(c3t3, cell_selector,
                                                   boundary_edges,
                                                   vertices_subdomain_indices);
  }
  else
#endif
  {
    for (const Edge& e : tr.finite_edges())
    {
      if (is_boundary(c3t3, e, cell_selector))
        boundary_edges.push_back(e);
    }

    // collect incident subdomain indices at vertices
    collect_subdomains_on_boundary(c3t3, vertices_subdomain_indices);
  }

#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
  CGAL::Tetrahedral_remeshing::debug::dump_edges(boundary_edges,
                                                 "boundary_edges.polylines.txt");
#endif

#ifdef CGAL_LINKED_WITH_TBB
  if constexpr (is_parallel_triangulation<typename C3T3::Triangulation>())
  {
    compute_boundary_vertices_valences_parallel(c3t3, boundary_edges,
                                                boundary_vertices_valences,
                                                vertices_subdomain_indices);
    return;
  }
#endif

  for (const Edge& e : boundary_edges)
  {
    const Vertex_handle v0 = e.first->vertex(e.second);
    const Vertex_handle v1 = e.first->vertex(e.third);

    const std::size_t n0 = vertices_subdomain_indices[v0].size();
    const std::size_t n1 = vertices_subdomain_indices[v1].size();

    //In case of feature edge
    if (n0 > 2 && n1 > 2)
    {
      Facet_circulator facet_circulator = tr.incident_facets(e);
      Facet_circulator done(facet_circulator);
      do
      {
        if (c3t3.is_in_complex(*facet_circulator))
        {
          Surface_patch_index surfi = c3t3.surface_patch_index(*facet_circulator);
          boundary_vertices_valences[v0][surfi]++;
          boundary_vertices_valences[v1][surfi]++;
        }
      } while (++facet_circulator != done);
    }
    //Normal surface edge, or non-manifold edge on dangling facet
    else
    {
      Facet_circulator facet_circulator = tr.incident_facets(e);
      Facet_circulator done(facet_circulator);
      Surface_patch_index first_patch = Surface_patch_index();
      do
      {
        if (c3t3.is_in_complex(*facet_circulator))
        {
          Surface_patch_index surfi = c3t3.surface_patch_index(*facet_circulator);
          if (first_patch == Surface_patch_index())
            first_patch = surfi;
          else if (first_patch == surfi)
            continue;

          boundary_vertices_valences[v0][surfi]++;
          boundary_vertices_valences[v1][surfi]++;
        }
      } while (++facet_circulator != done);
    }
  }
}

template<typename C3T3, typename IncCellsVector, typename Visitor>
Sliver_removal_result flip_n_to_m_on_surface(typename C3T3::Edge& edge,
    C3T3& c3t3,
    typename C3T3::Vertex_handle v0i,//v0 of new edge that will replace edge
    typename C3T3::Vertex_handle v1i,//v1 of new edge that will replace edge
    const IncCellsVector& cells_around_edge,
    Flip_Criterion /*flip_criterion*/,
    Visitor& /*visitor*/)
{
  typedef typename C3T3::Vertex_handle Vertex_handle;
  typedef typename C3T3::Cell_handle   Cell_handle;

  typename C3T3::Triangulation& tr = c3t3.triangulation();

  const Vertex_handle u = edge.first->vertex(edge.second);
  const Vertex_handle v = edge.first->vertex(edge.third);

  typedef std::pair<int, int> IndInCell;
  std::map<Cell_handle, IndInCell> indices;
  for (Cell_handle c : cells_around_edge)
  {
    indices[c] = std::make_pair(c->index(u), c->index(v));
  }

  for (Cell_handle c : cells_around_edge)
  {
    int i = indices[c].first;
    int j = indices[c].second;
    c->set_vertex(i, v0i);
    c->set_vertex(j, v1i);

    if (!is_well_oriented(tr, c))
    {
      c->set_vertex(j, v0i);
      c->set_vertex(i, v1i);
      if (!is_well_oriented(tr, c))
      {
        //rollback all changes
        for (Cell_handle cc : cells_around_edge)
        {
          const int ii = indices[cc].first;
          if (cc->vertex(ii) != u)
          {
            cc->set_vertex(ii, u);
            cc->set_vertex(indices[cc].second, v);
          }
        }
        return NOT_FLIPPABLE;
      }
    }
  }

  for (Cell_handle c : cells_around_edge)
    c->reset_cache_validity();

  return VALID_FLIP;
}

//v0i and v1i are the vertices opposite to `edge`
//on facets of the surface
template<typename C3T3, typename IncCellsVectorMap,
         typename Visitor>
Sliver_removal_result flip_on_surface(C3T3& c3t3,
    typename C3T3::Edge& edge,
    const typename C3T3::Vertex_handle v0i,//v0 of new edge that will replace edge
    const typename C3T3::Vertex_handle v1i,//v1 of new edge that will replace edge
    IncCellsVectorMap& inc_cells,
    Flip_Criterion flip_criterion,
    Visitor& visitor)
{
  typedef typename C3T3::Triangulation Tr;
  typedef typename Tr::Cell_handle     Cell_handle;
  typedef typename Tr::Vertex_handle   Vertex_handle;
  typedef typename Tr::Cell_circulator Cell_circulator;

  Tr& tr = c3t3.triangulation();
  Cell_circulator circ = tr.incident_cells(edge);
  Cell_circulator done(circ);

  std::vector<Cell_handle> cells_around_edge;
  do
  {
    cells_around_edge.push_back(circ);
  } while (++circ != done);

  if (cells_around_edge.size() != 4)
  {
    if (cells_around_edge.size() > 4){
//////      if (flip_criterion == VALENCE_BASED){
//////        return find_best_n_m_flip(edge, vh0_index, vh1_index);
//////      }
//////      else {
//        std::vector<Vertex_handle> boundary_vertices;
//        boundary_vertices.push_back(v0i);
//        boundary_vertices.push_back(v1i);
#ifdef CGAL_FLIP_ON_SURFACE_DISABLE_NM_FLIP
        return NOT_FLIPPABLE;
#else
        return flip_n_to_m_on_surface(edge, c3t3, v0i, v1i,
                                      cells_around_edge, flip_criterion,
                                      visitor);
#endif
//////    }
    }
    else
      return NOT_FLIPPABLE;
  }

#ifdef CGAL_FLIP_ON_SURFACE_DISABLE_44_FLIP
  return NOT_FLIPPABLE;
#endif

  inc_cells[edge.first->vertex(edge.second)].clear();
  inc_cells[edge.first->vertex(edge.third)].clear();

  Cell_handle ch0, ch1, ch2, ch3;
  ch0 = cells_around_edge[0];
  ch1 = cells_around_edge[1];
  ch2 = cells_around_edge[2];
  ch3 = cells_around_edge[3];

  Dihedral_angle_cosine curr_max_cosdh = max_cos_dihedral_angle(tr, ch0);
  for (int i = 1; i < 4; ++i)
    curr_max_cosdh = (std::max)(curr_max_cosdh,
                                max_cos_dihedral_angle(tr, cells_around_edge[i]));

  Vertex_handle vh0, vh1, vh2, vh3, vh4, vh5;

  int ivh4_in_ch0 = ch0->index(ch1);
  vh4 = ch0->vertex(ivh4_in_ch0);

  int ivh2_in_ch0 = ch0->index(ch3);
  vh2 = ch0->vertex(ivh2_in_ch0);

  vh5 = ch1->vertex(ch1->index(ch0));
  vh0 = ch2->vertex(ch2->index(ch1));

  for (int j = 0; j < 3; j++){
    if (indices(ivh4_in_ch0, j) == ivh2_in_ch0){
      int j1 = (j + 1) % 3;
      int j2 = (j + 2) % 3;
      vh1 = ch0->vertex(indices(ivh4_in_ch0, j1));
      vh3 = ch0->vertex(indices(ivh4_in_ch0, j2));
      break;
    }
  }

  bool planar_flip;
  if ((vh0 == v0i && vh2 == v1i) || (vh2 == v0i && vh0 == v1i))
    planar_flip = true;
  else if ((vh4 == v0i && vh5 == v1i) || (vh5 == v0i && vh4 == v1i))
    planar_flip = false;
  else
    return NOT_FLIPPABLE;

  typedef typename C3T3::Facet Facet;
  typedef typename C3T3::Surface_patch_index Surface_patch_index;

  if (planar_flip)
  {
#ifdef CGAL_FLIP_ON_SURFACE_DISABLE_PLANAR_44_FLIP
    return NOT_FLIPPABLE;
#endif
    Surface_patch_index patch = c3t3.surface_patch_index(ch0, ch0->index(vh4));
    CGAL_assertion(patch != Surface_patch_index());
    CGAL_assertion(c3t3.is_in_complex(ch0, ch0->index(vh4)));
    c3t3.remove_from_complex(ch0, ch0->index(vh4));
    CGAL_assertion(c3t3.is_in_complex(ch3, ch3->index(vh4)));
    c3t3.remove_from_complex(ch3, ch3->index(vh4));

    boost::unordered_map<Facet, Surface_patch_index> opposite_facet_in_complex;
    for (Cell_handle chi : cells_around_edge)
    {
      Facet f1(chi, chi->index(vh1));
      Facet f2(chi, chi->index(vh3));

      if (c3t3.is_in_complex(f1))
      {
        Surface_patch_index spi = c3t3.surface_patch_index(f1);
        opposite_facet_in_complex[c3t3.triangulation().mirror_facet(f1)] = spi;
        c3t3.remove_from_complex(f1);
      }
      if (c3t3.is_in_complex(f2))
      {
        Surface_patch_index spi = c3t3.surface_patch_index(f2);
        opposite_facet_in_complex[c3t3.triangulation().mirror_facet(f2)] = spi;
        c3t3.remove_from_complex(f2);
      }
    }

    Cell_handle n_ch3_vh1 = ch3->neighbor(ch3->index(vh1));
    Cell_handle n_ch0_vh3 = ch0->neighbor(ch0->index(vh3));

    Cell_handle n_ch2_vh1 = ch2->neighbor(ch2->index(vh1));
    Cell_handle n_ch1_vh3 = ch1->neighbor(ch1->index(vh3));

    ch3->set_vertex(ch3->index(vh3), vh2);
    ch0->set_vertex(ch0->index(vh1), vh0);
    ch2->set_vertex(ch2->index(vh3), vh2);
    ch1->set_vertex(ch1->index(vh1), vh0);

    Sliver_removal_result db = VALID_FLIP;
    if (!is_well_oriented(tr, ch0)
      || !is_well_oriented(tr, ch1)
      || !is_well_oriented(tr, ch2)
      || !is_well_oriented(tr, ch3))
      db = NOT_FLIPPABLE;
    else if (curr_max_cosdh < max_cos_dihedral_angle(tr, ch0, false)
      || curr_max_cosdh < max_cos_dihedral_angle(tr, ch1, false)
      || curr_max_cosdh < max_cos_dihedral_angle(tr, ch2, false)
      || curr_max_cosdh < max_cos_dihedral_angle(tr, ch3, false))
      db = NO_BEST_CONFIGURATION;

    if(db != VALID_FLIP)
    {
      ch3->set_vertex(ch3->index(vh2), vh3);
      ch0->set_vertex(ch0->index(vh0), vh1);
      ch2->set_vertex(ch2->index(vh2), vh3);
      ch1->set_vertex(ch1->index(vh0), vh1);

      c3t3.add_to_complex(ch0, ch0->index(vh4), patch);
      c3t3.add_to_complex(ch3, ch3->index(vh4), patch);

      for (Cell_handle chi : cells_around_edge)
      {
        Facet f1(chi, chi->index(vh1));
        Facet f2(chi, chi->index(vh3));

        auto it = opposite_facet_in_complex.find(c3t3.triangulation().mirror_facet(f1));
        if (it != opposite_facet_in_complex.end())
          c3t3.add_to_complex(f1, it->second);

        it = opposite_facet_in_complex.find(c3t3.triangulation().mirror_facet(f2));
        if (it != opposite_facet_in_complex.end())
          c3t3.add_to_complex(f2, it->second);
      }

      return db;
    }

    //Top cells 2-2 flip
    ch3->set_neighbor(ch3->index(vh1), ch0);
    ch3->set_neighbor(ch3->index(vh0), n_ch0_vh3);
    n_ch0_vh3->set_neighbor(n_ch0_vh3->index(ch0), ch3);

    ch0->set_neighbor(ch0->index(vh3), ch3);
    ch0->set_neighbor(ch0->index(vh2), n_ch3_vh1);
    n_ch3_vh1->set_neighbor(n_ch3_vh1->index(ch3), ch0);

    //Bottom cells 2-2 flip
    ch2->set_neighbor(ch2->index(vh1), ch1);
    ch2->set_neighbor(ch2->index(vh0), n_ch1_vh3);
    n_ch1_vh3->set_neighbor(n_ch1_vh3->index(ch1), ch2);

    ch1->set_neighbor(ch1->index(vh3), ch2);
    ch1->set_neighbor(ch1->index(vh2), n_ch2_vh1);
    n_ch2_vh1->set_neighbor(n_ch2_vh1->index(ch2), ch1);

    for (Cell_handle ci : cells_around_edge)
    {
      for (int j = 0; j < 4; j++)
        ci->vertex(j)->set_cell(ci);
    }

    c3t3.add_to_complex(ch0, ch0->index(vh4), patch);
    c3t3.add_to_complex(ch3, ch3->index(vh4), patch);

    for (Cell_handle chi : cells_around_edge)
    {
      Facet f1(chi, chi->index(vh0));
      Facet f2(chi, chi->index(vh2));

      auto it = opposite_facet_in_complex.find(c3t3.triangulation().mirror_facet(f1));
      if (it != opposite_facet_in_complex.end())
        c3t3.add_to_complex(f1, it->second);

      it = opposite_facet_in_complex.find(c3t3.triangulation().mirror_facet(f2));
      if (it != opposite_facet_in_complex.end())
        c3t3.add_to_complex(f2, it->second);
    }

    for(Cell_handle c : cells_around_edge)
      c->reset_cache_validity();

    return db;
  }
  else //Non planar flip
  {
#ifdef CGAL_FLIP_ON_SURFACE_DISABLE_NON_PLANAR_44_FLIP
    return NOT_FLIPPABLE;
#endif
    typename C3T3::Surface_patch_index patch = c3t3.surface_patch_index(ch0, ch0->index(vh2));
    CGAL_assertion(patch != typename C3T3::Surface_patch_index());

    CGAL_assertion(c3t3.is_in_complex(ch0, ch0->index(vh2)));
    c3t3.remove_from_complex(ch0, ch0->index(vh2));
    CGAL_assertion(c3t3.is_in_complex(ch1, ch1->index(vh2)));
    c3t3.remove_from_complex(ch1, ch1->index(vh2));

    boost::unordered_map<Facet, Surface_patch_index> opposite_facet_in_complex;
    for (Cell_handle chi : cells_around_edge)
    {
      Facet f1(chi, chi->index(vh1));
      Facet f2(chi, chi->index(vh3));

      if (c3t3.is_in_complex(f1))
      {
        Surface_patch_index spi = c3t3.surface_patch_index(f1);
        opposite_facet_in_complex[c3t3.triangulation().mirror_facet(f1)] = spi;
        c3t3.remove_from_complex(f1);
      }
      if (c3t3.is_in_complex(f2))
      {
        Surface_patch_index spi = c3t3.surface_patch_index(f2);
        opposite_facet_in_complex[c3t3.triangulation().mirror_facet(f2)] = spi;
        c3t3.remove_from_complex(f2);
      }
    }

    // Top Flip
    ch3->set_vertex(ch3->index(vh1), vh5);
    ch2->set_vertex(ch2->index(vh3), vh4);
    ch0->set_vertex(ch0->index(vh1), vh5);
    ch1->set_vertex(ch1->index(vh3), vh4);

    Sliver_removal_result db = VALID_FLIP;
    if (!is_well_oriented(tr, ch0)
      || !is_well_oriented(tr, ch1)
      || !is_well_oriented(tr, ch2)
      || !is_well_oriented(tr, ch3))
      db = NOT_FLIPPABLE;
    else if (curr_max_cosdh < max_cos_dihedral_angle(tr, ch0, false)
      || curr_max_cosdh < max_cos_dihedral_angle(tr, ch1, false)
      || curr_max_cosdh < max_cos_dihedral_angle(tr, ch2, false)
      || curr_max_cosdh < max_cos_dihedral_angle(tr, ch3, false))
      db = NO_BEST_CONFIGURATION;

    if (db == NOT_FLIPPABLE || db == NO_BEST_CONFIGURATION)
    {
      ch3->set_vertex(ch3->index(vh5), vh1);
      ch2->set_vertex(ch2->index(vh4), vh3);
      ch0->set_vertex(ch0->index(vh5), vh1);
      ch1->set_vertex(ch1->index(vh4), vh3);

      c3t3.add_to_complex(ch0, ch0->index(vh2), patch);
      c3t3.add_to_complex(ch1, ch1->index(vh2), patch);

      for (Cell_handle chi : cells_around_edge)
      {
        Facet f1(chi, chi->index(vh1));
        Facet f2(chi, chi->index(vh3));

        auto it = opposite_facet_in_complex.find(c3t3.triangulation().mirror_facet(f1));
        if (it != opposite_facet_in_complex.end())
          c3t3.add_to_complex(f1, it->second);

        it = opposite_facet_in_complex.find(c3t3.triangulation().mirror_facet(f2));
        if (it != opposite_facet_in_complex.end())
          c3t3.add_to_complex(f2, it->second);
      }

      return db;
    }

    //Left cells 2-2 flip
    Cell_handle n_ch3_vh3 = ch3->neighbor(ch3->index(vh3));
    Cell_handle n_ch2_vh1 = ch2->neighbor(ch2->index(vh1));

    ch3->set_neighbor(ch3->index(vh3), ch2);
    ch3->set_neighbor(ch3->index(vh4), n_ch2_vh1);
    n_ch2_vh1->set_neighbor(n_ch2_vh1->index(ch2), ch3);

    ch2->set_neighbor(ch2->index(vh1), ch3);
    ch2->set_neighbor(ch2->index(vh5), n_ch3_vh3);
    n_ch3_vh3->set_neighbor(n_ch3_vh3->index(ch3), ch2);

    //Right cells 2-2 flip
    Cell_handle n_ch0_vh3 = ch0->neighbor(ch0->index(vh3));
    Cell_handle n_ch1_vh1 = ch1->neighbor(ch1->index(vh1));

    ch0->set_neighbor(ch0->index(vh3), ch1);
    ch0->set_neighbor(ch0->index(vh4), n_ch1_vh1);
    n_ch1_vh1->set_neighbor(n_ch1_vh1->index(ch1), ch0);

    ch1->set_neighbor(ch1->index(vh1), ch0);
    ch1->set_neighbor(ch1->index(vh5), n_ch0_vh3);
    n_ch0_vh3->set_neighbor(n_ch0_vh3->index(ch0), ch1);

    for (const Cell_handle ce : cells_around_edge)
    {
      for (int j = 0; j < 4; j++)
        ce->vertex(j)->set_cell(ce);
    }

    c3t3.add_to_complex(ch0, ch0->index(vh2), patch);
    c3t3.add_to_complex(ch1, ch1->index(vh2), patch);

    for (Cell_handle chi : cells_around_edge)
    {
      Facet f1(chi, chi->index(vh4));
      Facet f2(chi, chi->index(vh5));

      auto it = opposite_facet_in_complex.find(c3t3.triangulation().mirror_facet(f1));
      if (it != opposite_facet_in_complex.end())
        c3t3.add_to_complex(f1, it->second);

      it = opposite_facet_in_complex.find(c3t3.triangulation().mirror_facet(f2));
      if (it != opposite_facet_in_complex.end())
        c3t3.add_to_complex(f2, it->second);
    }

    for (Cell_handle c : cells_around_edge)
      c->reset_cache_validity();

    return VALID_FLIP;
  }

  return NOT_FLIPPABLE;
}

template<typename C3T3, typename SurfaceIndexMapMap,
         typename IncidentCellsVectorMap, typename Flip_Criterion,
         typename CellSelector, typename Visitor>
std::size_t flipBoundaryEdges(C3T3& c3t3,
                              const std::vector<typename C3T3::Edge>& boundary_edges,
                              SurfaceIndexMapMap& boundary_vertices_valences,
                              IncidentCellsVectorMap& inc_cells,
                              const Flip_Criterion& flip_criterion,
                              CellSelector& cell_selector,
                              Visitor& visitor)
{
  typedef typename C3T3::Vertex_handle Vertex_handle;
  typedef typename C3T3::Cell_handle Cell_handle;
  typedef typename C3T3::Facet Facet;
  typedef typename C3T3::Edge Edge;
  typedef typename C3T3::Triangulation Tr;
  typedef std::pair<Vertex_handle, Vertex_handle> Edge_vv;

  std::size_t nb_success = 0;

  Tr& tr = c3t3.triangulation();

  std::vector<Edge_vv> candidate_edges_for_flip;
  for(const Edge& e : boundary_edges) {
    if(!c3t3.is_in_complex(e))
      candidate_edges_for_flip.push_back(make_vertex_pair(e));
  }

  for(const auto& [vh0, vh1] : candidate_edges_for_flip) {
    boost::container::small_vector<Cell_handle, 64>& inc_vh0 = inc_cells[vh0];
    if(inc_vh0.empty())
      tr.incident_cells(vh0, std::back_inserter(inc_vh0));

    Cell_handle c;
    int i, j;
    if(!is_edge_uv(vh0, vh1, inc_vh0, c, i, j))
      continue;

    Edge edge(c, i, j);
    std::vector<Facet> boundary_facets;
    const bool on_boundary = is_boundary_edge(edge, c3t3, cell_selector, boundary_facets);

    //    if (on_boundary && boundary_facets.empty())
    //    {
    //      std::cerr << vh0->point().point() << "\t " << vh1->point().point() << std::endl;
    //      bool b = is_boundary_edge(vh0, vh1, c3t3, cell_selector);
    //      CGAL::Tetrahedral_remeshing::debug::dump_c3t3(c3t3, "dump_c3t3_about_boundary_");
    //      CGAL::Tetrahedral_remeshing::debug::dump_facets_in_complex(c3t3, "dump_facets_about_boundary_.off");
    //      CGAL::Tetrahedral_remeshing::debug::dump_facets_from_selection(
    //        c3t3, cell_selector, "dump_facets_from_selection_.off");
    //      std::cerr << "valid = " << tr.tds().is_valid(true) << std::endl;
    //      std::cerr << "boundary = " << b << std::endl;
    //      CGAL_assertion(on_boundary);
    //    }
    //    else if (on_boundary && boundary_facets.size() != 2)
    //    {
    //      std::cerr << vh0->point().point() << "\t " << vh1->point().point() << std::endl;
    //      CGAL::Tetrahedral_remeshing::debug::dump_c3t3(c3t3, "dump_c3t3_about_boundary_");
    //      CGAL::Tetrahedral_remeshing::debug::dump_facets(boundary_facets, "dump_boundary_facets.polylines.txt");
    //      std::vector<Facet> dummy_facets;
    //      bool b = is_boundary_edge(edge, c3t3, cell_selector, dummy_facets, true/**/);
    //      std::cerr << "boundary = " << b << std::endl;
    //    }

    if(!on_boundary || boundary_facets.size() != 2)
      continue;
    CGAL_assertion(boundary_facets.size() == 2);

    if(flip_surface_edge(c3t3, edge, boundary_facets, boundary_vertices_valences,
                          inc_cells, flip_criterion, visitor))
    {
      ++nb_success;
      CGAL_expensive_assertion(tr.tds().is_valid());
    }
  }
  return nb_success;
}

//todo : write this function
template <typename C3t3,
          typename BV_valences,
          typename IncidentCellsVectorMap,
          typename Flip_criterion,
          typename Visitor>
bool flip_surface_edge(C3t3& c3t3,
                        typename C3t3::Edge& edge,
                        const std::vector<typename C3t3::Facet>& boundary_facets,
                        BV_valences& boundary_vertices_valences,
                        IncidentCellsVectorMap& inc_cells,
                        const Flip_criterion& flip_criterion,
                        Visitor& visitor)
{
  using Vertex_handle       = typename C3t3::Vertex_handle;
  using Facet               = typename C3t3::Facet;
  using Cell_handle         = typename C3t3::Cell_handle;
  using Surface_patch_index = typename C3t3::Surface_patch_index;

  auto& tr = c3t3.triangulation();

    const Facet& f0 = boundary_facets[0];
    const Facet& f1 = boundary_facets[1];

    const Vertex_handle vh0 = edge.first->vertex(edge.second);
    const Vertex_handle vh1 = edge.first->vertex(edge.third);

    // find 3rd and 4th vertices to flip on surface
    const Vertex_handle vh2 = third_vertex(f0, vh0, vh1, tr);
    const Vertex_handle vh3 = third_vertex(f1, vh0, vh1, tr);

    CGAL_expensive_assertion(debug::check_facets(vh0, vh1, vh2, vh3, c3t3));

    // Not tds().is_edge(): that walk MARKS tds_data() on every cell it
    // visits, and it visits the stars of vh2 and vh3, which no lock zone
    // covers -- the same reason the split and the collapse keep it off the
    // parallel path. Both vertices are corners of cells of the zone, so the
    // non-marking walk reads a star this thread holds the apex of.
    bool vh2_vh3_share_an_edge;
    if constexpr (is_parallel_triangulation<typename C3t3::Triangulation>())
      vh2_vh3_share_an_edge = tr.is_edge_threadsafe(vh2, vh3);
    else
      vh2_vh3_share_an_edge = tr.tds().is_edge(vh2, vh3);

    if (!vh2_vh3_share_an_edge) // most-likely to happen early exit
    {
      const Surface_patch_index surfi = c3t3.surface_patch_index(boundary_facets[0]);

      int v0 = boundary_vertices_valences.at(vh0)[surfi];
      int v1 = boundary_vertices_valences.at(vh1)[surfi];
      int v2 = boundary_vertices_valences.at(vh2)[surfi];
      int v3 = boundary_vertices_valences.at(vh3)[surfi];

      if(v0 < 2 || v1 < 2 || v2 < 2 || v3 < 2)
        return false;

      int m0 = (boundary_vertices_valences.at(vh0).size() > 1 ? 4 : 6);
      int m1 = (boundary_vertices_valences.at(vh1).size() > 1 ? 4 : 6);
      int m2 = (boundary_vertices_valences.at(vh2).size() > 1 ? 4 : 6);
      int m3 = (boundary_vertices_valences.at(vh3).size() > 1 ? 4 : 6);

      int initial_cost = (v0 - m0)*(v0 - m0)
                       + (v1 - m1)*(v1 - m1)
                       + (v2 - m2)*(v2 - m2)
                       + (v3 - m3)*(v3 - m3);
      v0--;
      v1--;
      v2++;
      v3++;

      int final_cost = (v0 - m0)*(v0 - m0)
                     + (v1 - m1)*(v1 - m1)
                     + (v2 - m2)*(v2 - m2)
                     + (v3 - m3)*(v3 - m3);
      if (initial_cost > final_cost)
      {
        CGAL_expensive_assertion_code(std::size_t nbf =
          std::distance(c3t3.facets_in_complex_begin(),
                        c3t3.facets_in_complex_end()));
        CGAL_expensive_assertion_code(std::size_t nbe =
          std::distance(c3t3.edges_in_complex_begin(),
                        c3t3.edges_in_complex_end()));

        Sliver_removal_result db = flip_on_surface(c3t3, edge, vh2, vh3,
                                                   inc_cells,
                                                   flip_criterion,
                                                   visitor);
        if (db == VALID_FLIP)
        {
          CGAL_expensive_assertion(tr.tds().is_edge(vh2, vh3));
          Cell_handle c;
          int li, lj, lk;
          CGAL_expensive_assertion_code(bool b =)
          tr.tds().is_facet(vh2, vh3, vh0, c, li, lj, lk);
          CGAL_expensive_assertion(b);
          c3t3.add_to_complex(c, (6 - li - lj - lk), surfi);

          CGAL_expensive_assertion_code(b = )
          tr.tds().is_facet(vh2, vh3, vh1, c, li, lj, lk);
          CGAL_expensive_assertion(b);
          c3t3.add_to_complex(c, (6 - li - lj - lk), surfi);

          CGAL_expensive_assertion_code(std::size_t nbf_post =
            std::distance(c3t3.facets_in_complex_begin(),
                          c3t3.facets_in_complex_end()));
          CGAL_expensive_assertion(nbf == nbf_post);
          CGAL_expensive_assertion_code(std::size_t nbe_post =
            std::distance(c3t3.edges_in_complex_begin(),
                          c3t3.edges_in_complex_end()));
          CGAL_expensive_assertion(nbe == nbe_post);

          boundary_vertices_valences[vh0][surfi]--;
          boundary_vertices_valences[vh1][surfi]--;
          boundary_vertices_valences[vh2][surfi]++;
          boundary_vertices_valences[vh3][surfi]++;

          return true;
        }
      }
    }
  return false;
}

// Shared state for the internal and boundary edge-flip operations: the cell
// selector, the visitor, and the incident-cells cache used by find_best_flip
// and flip_on_surface.
template <typename C3t3,
          typename CellSelector,
          typename Visitor>
class Edge_flip_operation_base
{
protected:
  using Tr            = typename C3t3::Triangulation;
  using Cell_handle   = typename Tr::Cell_handle;
  using Vertex_handle = typename Tr::Vertex_handle;
  using Edge          = typename Tr::Edge;
  using Facet         = typename Tr::Facet;
  using Cells_vector  = boost::container::small_vector<Cell_handle, 64>;
  using Concurrency_tag = typename Tr::Concurrency_tag;

  // The cache is written by every thread that flips, so it must tolerate
  // concurrent insertion when remeshing in parallel. Each *value* is only ever
  // touched by the thread holding the lock zone of the vertex it belongs to.
#ifdef CGAL_LINKED_WITH_TBB
  using Incident_cells_map = Concurrency_selected_container_t<
      Concurrency_tag,
      std::unordered_map<Vertex_handle, Cells_vector>,
      tbb::concurrent_unordered_map<Vertex_handle, Cells_vector> >;
#else
  using Incident_cells_map = std::unordered_map<Vertex_handle, Cells_vector>;
#endif

  CellSelector& m_cell_selector;
  Visitor& m_visitor;
  // Shared across the internal and boundary flip passes, exactly as the former
  // flip_edges() shared a single inc_cells map between flip_all_edges() and
  // flipBoundaryEdges().
  Incident_cells_map& inc_cells;

  Edge_flip_operation_base(CellSelector& cell_selector,
                        Visitor& visitor,
                        Incident_cells_map& incident_cells)
      : m_cell_selector(cell_selector)
      , m_visitor(visitor)
      , inc_cells(incident_cells) {}
};

// ---------------------------------------------------------------------------
// The minimum viable lock zone of an internal flip.
//
// The zone taken by default is the union of the two endpoint stars. A flip
// only ever writes the ring cells -- the cells around the flipped edge,
// rewritten in full -- and their mirror cells, the cells across the ring's
// outer facets, which it re-stitches through set_neighbor() and
// update_c3t3_facets(). Ring and mirrors are a strict subset of the two stars:
// the ring IS star(v0) n star(v1), and a mirror cell still contains an
// endpoint, so it is a star cell too.
//
// `lock_zone()` below runs only under `Parallel_tag`. The `Located_edge`
// it hands to `execute_operation()` lives in tetrahedral_remeshing_helpers.h,
// because the split's lock zone does the same thing.

// Flip of internal (non-boundary) edges. Mirrors the former flip_all_edges():
// reset the cell caches, collect the internal edges, then run find_best_flip
// on each in turn, sharing the incident-cells cache.
template <typename C3t3, typename CellSelector, typename Visitor>
class Internal_edge_flip_operation
    : public Edge_flip_operation_base<C3t3, CellSelector, Visitor>,
      public Elementary_operation<C3t3,
                                 std::pair<typename C3t3::Vertex_handle, typename C3t3::Vertex_handle>,
                                 std::vector<std::pair<typename C3t3::Vertex_handle, typename C3t3::Vertex_handle>>>
{
  using BaseClass = Edge_flip_operation_base<C3t3, CellSelector, Visitor>;
  using typename BaseClass::Cell_handle;
  using typename BaseClass::Vertex_handle;
  using typename BaseClass::Edge;
  using typename BaseClass::Cells_vector;
  using BaseClass::m_cell_selector;
  using BaseClass::m_visitor;
  using BaseClass::inc_cells;

public:
  using Incident_cells_map = typename BaseClass::Incident_cells_map;
  using Edge_vv = std::pair<Vertex_handle, Vertex_handle>;
  using Base_operation = Elementary_operation<C3t3, Edge_vv, std::vector<Edge_vv>>;
  using Element_type = typename Base_operation::Element_type;
  static_assert(std::is_same_v<Element_type, Edge_vv>, "Element_type must be Edge_vv");
  using Element_range = typename Base_operation::Element_range;

  Internal_edge_flip_operation(CellSelector& cell_selector,
                            Visitor& visitor,
                            Incident_cells_map& incident_cells)
      : BaseClass(cell_selector, visitor, incident_cells) {}

  Element_range get_elements(const C3t3& c3t3) const override
  {
#ifdef CGAL_LINKED_WITH_TBB
    if constexpr (is_parallel)
    {
      // Both halves of this function -- the cache reset over the complex cells
      // and the edge collection -- run on one thread in front of a phase that
      // then runs on all of them, and together they are the largest serial
      // section left in the parallel run. They share ONE cell snapshot, so the
      // compact container is walked once here where the serial arm walks it
      // twice: once for `cells_in_complex()` and once for the edge iterator.
      const typename C3t3::Triangulation& tr = c3t3.triangulation();
      const std::vector<Cell_handle> cells = gather_all_cells(tr);

      tbb::parallel_for(tbb::blocked_range<std::size_t>(0, cells.size()),
        [&](const tbb::blocked_range<std::size_t>& range)
        {
          for (std::size_t ci = range.begin(); ci != range.end(); ++ci)
          {
            if (c3t3.is_in_complex(cells[ci]))
              cells[ci]->reset_cache_validity();//we will use sliver_value
                                                //to store the cos_dihedral_angle
          }
        });

      // Same per-edge test and the same edge set as get_internal_edges(); only
      // the traversal reaching each edge exactly once differs, and it applies
      // the smallest-incident-cell ownership rule the serial edge iterator
      // applies, in the same cell and slot order -- so the candidate order is
      // unchanged too. The per-edge test below is spelled exactly as
      // get_internal_edges() spells it.
      return parallel_collect_from_finite_edges<Edge_vv>(tr, cells,
        [&](const Edge& e, std::vector<Edge_vv>& out)
        {
          if (is_internal(e, c3t3, m_cell_selector))
            out.push_back(make_vertex_pair(e));
        });
    }
#endif

    for (auto c : c3t3.cells_in_complex())
      c->reset_cache_validity();//we will use sliver_value
                                //to store the cos_dihedral_angle

    std::vector<Edge_vv> inside_edges;
    get_internal_edges(c3t3, m_cell_selector, std::back_inserter(inside_edges));
    return inside_edges;
  }

  bool execute_operation(const Element_type& vp, C3t3& c3t3) override
  {
    Cell_handle ch;
    int i0, i1;
    if constexpr (is_parallel)
    {
      // lock_zone() ran first and located the edge in order to lock its ring.
      // No match means it found no cell carrying the edge, i.e. this pair is
      // no longer an edge and there is nothing to flip.
      const auto& located = last_located_edge<Vertex_handle, Cell_handle>();
      if (!located.matches(vp.first, vp.second))
        return false;
      ch = located.c;
      i0 = located.i0;
      i1 = located.i1;
    }
    else
    {
      Cells_vector& o_inc_vh = inc_cells[vp.first];
      if (o_inc_vh.empty())
        c3t3.triangulation().incident_cells(vp.first, std::back_inserter(o_inc_vh));

      if (!is_edge_uv(vp.first, vp.second, o_inc_vh, ch, i0, i1))
        return false;
    }

    Edge edge(ch, i0, i1);
    const Sliver_removal_result res
      = find_best_flip(edge, c3t3, MIN_ANGLE_BASED, inc_cells, m_cell_selector, m_visitor);
    return (res == VALID_FLIP);
  }

  // `lock_zone()` is called by the parallel executor alone; the sequential one
  // goes straight to `execute_operation()`, which must then locate the edge
  // itself.
  static constexpr bool is_parallel
    = std::is_convertible_v<typename BaseClass::Concurrency_tag, CGAL::Parallel_tag>;

  /**
  * Locks everything `execute_operation()` may write, and nothing else.
  *
  * A flip rewrites the RING -- the cells around the flipped edge -- in full,
  * and re-stitches the MIRROR cells across the ring's outer facets through
  * `set_neighbor()` and `update_c3t3_facets()`. Cells created by the flip have
  * only ring vertices and one ring apex as corners. Nothing else is written:
  * the rest of the two endpoint stars is read, never modified.
  *
  * A cell is protected by holding its four vertices, so the zone is the vertex
  * union of the ring and its mirrors. That is a strict subset of the two
  * endpoint stars -- the ring IS star(v0) n star(v1), and a mirror cell still
  * contains an endpoint -- and it is reached here in O(ring) rather than
  * O(star): one early-exit search for a cell carrying the edge, then a
  * `Cell_circulator`, which rotates through neighbour pointers and so needs
  * neither a visited set nor marks.
  *
  * The two endpoints are locked first for two reasons: `v->cell()` is only
  * dependable once `v` is held, and holding them is what makes the search safe
  * to read -- every cell it visits contains `v0`, and a thread writing such a
  * cell would have to hold all four of its vertices, `v0` among them.
  *
  * Returning false leaves partial locks behind; the executor releases them all
  * before retrying.
  */
  bool lock_zone(const Element_type& vp, const C3t3& c3t3) const
  {
    typedef typename C3t3::Triangulation::Cell_circulator Cell_circulator;
    const typename C3t3::Triangulation& tr = c3t3.triangulation();

    last_located_edge<Vertex_handle, Cell_handle>().clear();

    if (!tr.try_lock_vertex(vp.first)
     || !tr.try_lock_vertex(vp.second))
      return false;

    Cell_handle edge_cell;
    const Vertex_handle other = vp.second;
    if (!tr.find_first_incident_cell_threadsafe(vp.first,
          [other](const Cell_handle c) { return c->has_vertex(other); },
          edge_cell))
      return true; // no longer an edge; execute_operation() will decline it

    const int i0 = edge_cell->index(vp.first);
    const int i1 = edge_cell->index(vp.second);
    const Edge edge(edge_cell, i0, i1);

    // One vertex per cell, not four. Both endpoints are already held, and a
    // ring cell shares three of its four vertices with the ring cell before
    // it -- consecutive cells of the circulator are neighbours -- so only the
    // fourth can still be unheld. A cell across an outer facet shares that
    // facet's three vertices with the ring cell, which is held in full by the
    // time it is reached, so the same is true of it. The first ring cell has
    // no predecessor and is locked in full.
    Cell_circulator circ = tr.incident_cells(edge);
    const Cell_circulator done = circ;
    Cell_handle previous_ring_cell;
    do
    {
      const Cell_handle c = circ;
      if (previous_ring_cell == Cell_handle())
      {
        if (!tr.try_lock_cell(c))
          return false;
      }
      else if (!tr.try_lock_vertex(c->vertex(c->index(previous_ring_cell))))
        return false;

      // the two cells across the ring cell's outer facets
      const Cell_handle m0 = c->neighbor(c->index(vp.first));
      const Cell_handle m1 = c->neighbor(c->index(vp.second));
      if (!tr.try_lock_vertex(m0->vertex(m0->index(c)))
       || !tr.try_lock_vertex(m1->vertex(m1->index(c))))
        return false;

      previous_ring_cell = c;
    }
    while (++circ != done);

    last_located_edge<Vertex_handle, Cell_handle>()
      .set(vp.first, vp.second, edge_cell, i0, i1);
    return true;
  }

  // no edge is preferred over another: shuffling spreads the threads out
  static constexpr bool requires_ordered_processing = false;

  std::string operation_name() const override { return "Flip edges (internal)"; }
};

// Flip of boundary edges. Mirrors the former flipBoundaryEdges(): compute the
// per-vertex boundary valences once, then, for each boundary edge, flip on the
// surface when it lowers the valence cost, updating the valences accordingly.
template <typename C3t3, typename CellSelector, typename Visitor>
class Boundary_edge_flip_operation
    : public Edge_flip_operation_base<C3t3, CellSelector, Visitor>,
      public Elementary_operation<C3t3,
                                 std::pair<typename C3t3::Vertex_handle, typename C3t3::Vertex_handle>,
                                 std::vector<std::pair<typename C3t3::Vertex_handle, typename C3t3::Vertex_handle>>>
{
  using BaseClass = Edge_flip_operation_base<C3t3, CellSelector, Visitor>;
  using typename BaseClass::Cell_handle;
  using typename BaseClass::Vertex_handle;
  using typename BaseClass::Edge;
  using typename BaseClass::Facet;
  using typename BaseClass::Cells_vector;
  using BaseClass::m_cell_selector;
  using BaseClass::m_visitor;
  using BaseClass::inc_cells;

  using Subdomain_index = typename C3t3::Subdomain_index;
  using Surface_patch_index = typename C3t3::Surface_patch_index;
  /**
  * The per-patch valences of ONE boundary vertex.
  *
  * A boundary vertex lies on a single surface patch almost everywhere and on a
  * handful at worst, so a `boost::unordered_map` here pays a hash and a node
  * allocation to hold, typically, one pair -- and it is asked ~4 times per
  * boundary edge while the valences are built, which is the serial half of the
  * boundary flip's preprocessing (0.608 s of a 41 s run on `fine_0.35`, at
  * 0.91x). A linear scan over an inline array answers the same three
  * questions: `operator[]` to read or bump a patch's valence, and `size() > 1`
  * to ask whether the vertex is on more than one patch.
  *
  * Same values, same insertion-on-miss semantics. The container is never
  * iterated, so no order is observable.
  */
  struct Spi_map
  {
    boost::container::small_vector<std::pair<Surface_patch_index, unsigned int>, 2> m_v;

    unsigned int& operator[](const Surface_patch_index& si)
    {
      for (auto& p : m_v)
        if (p.first == si)
          return p.second;
      m_v.emplace_back(si, 0u);
      return m_v.back().second;
    }
    std::size_t size() const { return m_v.size(); }
  };

  // Filled once by get_elements() for every boundary vertex. A flip then
  // decrements two valences and increments two others; all four vertices lie
  // in the stars locked by lock_zone(), so only the map's structure needs to
  // be safe against concurrent insertion.
#ifdef CGAL_LINKED_WITH_TBB
  using Boundary_valences_map = Concurrency_selected_container_t<
      typename BaseClass::Concurrency_tag,
      boost::unordered_map<Vertex_handle, Spi_map>,
      tbb::concurrent_unordered_map<Vertex_handle, Spi_map> >;
#else
  using Boundary_valences_map = boost::unordered_map<Vertex_handle, Spi_map>;
#endif

  mutable Boundary_valences_map m_boundary_vertices_valences;

public:
  using Incident_cells_map = typename BaseClass::Incident_cells_map;
  using Edge_vv = std::pair<Vertex_handle, Vertex_handle>;
  using Base_operation = Elementary_operation<C3t3, Edge_vv, std::vector<Edge_vv>>;
  using Element_type = typename Base_operation::Element_type;
  static_assert(std::is_same_v<Element_type, Edge_vv>, "Element_type must be Edge_vv");
  using Element_range = typename Base_operation::Element_range;

  Boundary_edge_flip_operation(CellSelector& cell_selector, Visitor& visitor, Incident_cells_map& incident_cells)
      : BaseClass(cell_selector, visitor, incident_cells) {}

  Element_range get_elements(const C3t3& c3t3) const override
  {
    std::vector<Edge> boundary_edges;
    Vertex_subdomains_map_t<C3t3> vertices_subdomain_indices;
    m_boundary_vertices_valences.clear();
    collectBoundaryEdgesAndComputeVerticesValences(c3t3,
                                                   m_cell_selector,
                                                   boundary_edges,
                                                   m_boundary_vertices_valences,
                                                   vertices_subdomain_indices);

#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
    if (!debug::are_cell_orientations_valid(c3t3.triangulation()))
      std::cerr << "ERROR in ORIENTATION" << std::endl;
#endif

    std::vector<Edge_vv> candidate_edges_for_flip;
    for (const Edge& e : boundary_edges)
    {
      if (!c3t3.is_in_complex(e))
        candidate_edges_for_flip.push_back(make_vertex_pair(e));
    }
    return candidate_edges_for_flip;
  }

  // `lock_zone()` is called by the parallel executor alone; the sequential one
  // goes straight to `execute_operation()`, which must then locate the edge
  // itself.
  static constexpr bool is_parallel
    = std::is_convertible_v<typename BaseClass::Concurrency_tag, CGAL::Parallel_tag>;

  bool execute_operation(const Element_type& vp, C3t3& c3t3) override
  {
    const Vertex_handle vh0 = vp.first;
    const Vertex_handle vh1 = vp.second;
    typename C3t3::Triangulation& tr = c3t3.triangulation();

    Cell_handle c;
    int i, j;
    if constexpr (is_parallel)
    {
      // lock_zone() ran first and located the edge in order to lock its ring.
      // No match means it found no cell carrying the edge, i.e. this pair is
      // no longer an edge and there is nothing to flip.
      const auto& located = last_located_edge<Vertex_handle, Cell_handle>();
      if (!located.matches(vh0, vh1))
        return false;
      c = located.c;
      i = located.i0;
      j = located.i1;
      CGAL_USE(tr);
    }
    else
    {
      Cells_vector& inc_vh0 = inc_cells[vh0];
      if (inc_vh0.empty())
        tr.incident_cells(vh0, std::back_inserter(inc_vh0));

      if (!is_edge_uv(vh0, vh1, inc_vh0, c, i, j))
        return false;
    }

    Edge edge(c, i, j);
    std::vector<Facet> boundary_facets;
    const bool on_boundary = is_boundary_edge(edge, c3t3, m_cell_selector, boundary_facets);

    if (!on_boundary || boundary_facets.size() != 2)
      return false;

    return flip_surface_edge(c3t3, edge,
                              boundary_facets,
                              m_boundary_vertices_valences,
                              inc_cells,
                              MIN_ANGLE_BASED,
                              m_visitor);
  }

  /**
  * The same minimal zone as the internal flip -- the ring and the cells across
  * its outer facets -- located the same way. See the comment there for why a
  * ring cell is protected by holding its four vertices, and why the endpoints
  * have to be held before the search may be read.
  *
  * A surface flip writes no more than an internal one does: the 4-4 case
  * rewrites the four ring cells and re-stitches the four cells across their
  * outer facets, `flip_n_to_m_on_surface()` does the same for a longer ring,
  * and `set_cell()` is called only on vertices of ring cells, which the ring's
  * own locks cover. The valences of the four vertices are updated in place,
  * and all four are corners of ring cells. The two reads that leave the ring
  * -- whether `vh2` and `vh3` already share an edge, and `is_boundary_edge()`
  * -- either stay on the ring and its mirrors or walk a star whose apex this
  * zone holds, which is why the first of them had to stop marking.
  *
  * The two whole endpoint stars this used to take were never the write set.
  * They were what `execute_operation()` searched for the edge, and it now
  * takes the located edge instead. Taking them cost O(star) per candidate, on
  * a mesh where one star reaches 8004 cells, plus a copy of each into
  * `inc_cells`.
  */
  bool lock_zone(const Element_type& vp, const C3t3& c3t3) const
  {
    typedef typename C3t3::Triangulation::Cell_circulator Cell_circulator;
    const typename C3t3::Triangulation& tr = c3t3.triangulation();

    last_located_edge<Vertex_handle, Cell_handle>().clear();

    if (!tr.try_lock_vertex(vp.first)
     || !tr.try_lock_vertex(vp.second))
      return false;

    Cell_handle edge_cell;
    const Vertex_handle other = vp.second;
    if (!tr.find_first_incident_cell_threadsafe(vp.first,
          [other](const Cell_handle c) { return c->has_vertex(other); },
          edge_cell))
      return true; // no longer an edge; execute_operation() will decline it

    const int i0 = edge_cell->index(vp.first);
    const int i1 = edge_cell->index(vp.second);
    const Edge edge(edge_cell, i0, i1);

    // One vertex per cell, not four. Both endpoints are already held, and a
    // ring cell shares three of its four vertices with the ring cell before
    // it -- consecutive cells of the circulator are neighbours -- so only the
    // fourth can still be unheld. A cell across an outer facet shares that
    // facet's three vertices with the ring cell, which is held in full by the
    // time it is reached, so the same is true of it. The first ring cell has
    // no predecessor and is locked in full.
    Cell_circulator circ = tr.incident_cells(edge);
    const Cell_circulator done = circ;
    Cell_handle previous_ring_cell;
    do
    {
      const Cell_handle c = circ;
      if (previous_ring_cell == Cell_handle())
      {
        if (!tr.try_lock_cell(c))
          return false;
      }
      else if (!tr.try_lock_vertex(c->vertex(c->index(previous_ring_cell))))
        return false;

      // the two cells across the ring cell's outer facets
      const Cell_handle m0 = c->neighbor(c->index(vp.first));
      const Cell_handle m1 = c->neighbor(c->index(vp.second));
      if (!tr.try_lock_vertex(m0->vertex(m0->index(c)))
       || !tr.try_lock_vertex(m1->vertex(m1->index(c))))
        return false;

      previous_ring_cell = c;
    }
    while (++circ != done);

    last_located_edge<Vertex_handle, Cell_handle>()
      .set(vp.first, vp.second, edge_cell, i0, i1);
    return true;
  }

  static constexpr bool requires_ordered_processing = false;

  std::string operation_name() const override { return "Flip edges (boundary)"; }
};

}//namespace internal
}//namespace Tetrahedral_remeshing
}//namespace CGAL

#endif // CGAL_INTERNAL_FLIP_EDGES_H
