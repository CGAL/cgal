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

#ifndef CGAL_INTERNAL_COLLAPSE_SHORT_EDGES_H
#define CGAL_INTERNAL_COLLAPSE_SHORT_EDGES_H

#include <CGAL/license/Tetrahedral_remeshing.h>

#include <boost/bimap.hpp>
#include <boost/bimap/set_of.hpp>
#include <boost/bimap/multiset_of.hpp>
#include <boost/container/small_vector.hpp>
#include <boost/functional/hash.hpp>

#include <vector>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <unordered_set>

#include <CGAL/SMDS_3/tet_soup_to_c3t3.h>
#include <CGAL/utility.h>
#include <CGAL/Tetrahedral_remeshing/internal/tetrahedral_remeshing_helpers.h>

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
#include <CGAL/Real_timer.h>
#endif

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
                     ANGLE_PROBLEM,
                     TOPOLOGICAL_PROBLEM, ORIENTATION_PROBLEM, SHARED_NEIGHBOR_PROBLEM };

template<typename C3t3>
class CollapseTriangulation
{
  typedef typename C3t3::Triangulation                        Tr;
  typedef typename C3t3::Edge                                 Edge;
  typedef typename C3t3::Cell_handle                          Cell_handle;
  typedef typename C3t3::Vertex_handle                        Vertex_handle;
  typedef typename C3t3::Subdomain_index                      Subdomain_index;
  typedef typename C3t3::Surface_patch_index                  Surface_patch_index;
  typedef typename C3t3::Triangulation::Point                 Point_3;
  typedef typename C3t3::Triangulation::Geom_traits::Vector_3 Vector_3;

public:
  CollapseTriangulation(const Edge& e,
                        const std::unordered_set<Cell_handle>& cells_to_insert,
                        Collapse_type _collapse_type)
    : collapse_type(_collapse_type)
    , v0_init(e.first->vertex(e.second))
    , v1_init(e.first->vertex(e.third))
  {
    typedef std::array<int, 3> Facet;
    typedef std::array<int, 4> Tet;

    std::unordered_set<Vertex_handle> vertices_to_insert;
    for (Cell_handle ch : cells_to_insert)
    {
      for(int i = 0; i < 4; ++i)
        vertices_to_insert.insert(ch->vertex(i));
    }

    CGAL_expensive_assertion(vertices_to_insert.end()
      != std::find(vertices_to_insert.begin(), vertices_to_insert.end(), v1_init));
    CGAL_expensive_assertion(vertices_to_insert.end()
      != std::find(vertices_to_insert.begin(), vertices_to_insert.end(), v0_init));

    std::unordered_map<Vertex_handle, int> v2i;/*vertex of main tr - vertex of collapse tr*/

    //To add the vertices only once
    std::vector<Point_3> points;
    int index = 0;
    for (Vertex_handle vh : vertices_to_insert)
    {
      if (v2i.find(vh) == v2i.end())
      {
        points.push_back(vh->point());
        v2i.insert(std::make_pair(vh, index++));
      }
    }

    std::vector<Tet> finite_cells;
    std::vector<int> subdomains;
    for (Cell_handle ch : cells_to_insert)
    {
      finite_cells.push_back( { v2i.at(ch->vertex(0)),
                                v2i.at(ch->vertex(1)),
                                v2i.at(ch->vertex(2)),
                                v2i.at(ch->vertex(3)) } );
      subdomains.push_back(ch->subdomain_index());
    }

    // finished
    std::vector<Vertex_handle> new_vertices;
    std::map<Facet, typename C3t3::Surface_patch_index> border_facets;
    if (CGAL::SMDS_3::build_triangulation_impl(
            triangulation, points, finite_cells, subdomains, border_facets,
            new_vertices, /*verbose*/ false,
            /*replace_domain_0*/ false,
            /*allow_non_manifold*/false))
    {
      CGAL_expensive_assertion(triangulation.tds().is_valid());
      CGAL_assertion(triangulation.infinite_vertex() == new_vertices[0]);

      // update()
      vh0 = new_vertices[v2i.at(v0_init) + 1];
      vh1 = new_vertices[v2i.at(v1_init) + 1];

      Cell_handle ch;
      int i0, i1;
      not_an_edge = true;
      CGAL_assertion(triangulation.tds().is_vertex(vh0));
      CGAL_assertion(triangulation.tds().is_vertex(vh1));
      if (triangulation.is_edge(vh0, vh1, ch, i0, i1))
      {
        edge = Edge(ch, i0, i1);
        not_an_edge = false;
      }
    }
#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
    else
      std::cout << "Warning : CollapseTriangulation is not valid!" << std::endl;
#endif
  }

  Result_type collapse()
  {
    if (not_an_edge)
    {
#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
      std::cout << "CollapseTriangulation::Not an edge..." << std::endl;
#endif
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

      std::unordered_set<Cell_handle> invalid_cells;

      typedef typename Tr::Cell_circulator Cell_circulator;
      Cell_circulator circ = triangulation.incident_cells(edge);
      Cell_circulator done = circ;

      std::vector<Cell_handle> cells_to_remove;

      //Update the vertex before removing it
      std::vector<Cell_handle> find_incident;
      triangulation.incident_cells(vh0, std::back_inserter(find_incident));

      std::vector<Cell_handle> cells_to_update;
      triangulation.incident_cells(vh1, std::back_inserter(cells_to_update));

      do
      {
        int v0_id = circ->index(vh0);
        int v1_id = circ->index(vh1);

        Cell_handle n0_ch = circ->neighbor(v0_id);
        Cell_handle n1_ch = circ->neighbor(v1_id);

        int ch_id_in_n0 = n0_ch->index(circ);
        int ch_id_in_n1 = n1_ch->index(circ);

        if (n0_ch->has_neighbor(n1_ch))
          return SHARED_NEIGHBOR_PROBLEM;

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

      // compute and keep worst angle
      Dihedral_angle_cosine curr_max_cos
        = (std::max)(max_cos_dihedral_angle_in_range(triangulation, cells_to_remove, false),
                     max_cos_dihedral_angle_in_range(triangulation, cells_to_update, false));


      vh0->set_point(Point_3(v0_new_pos.x(), v0_new_pos.y(), v0_new_pos.z()));
      vh1->set_point(Point_3(v0_new_pos.x(), v0_new_pos.y(), v0_new_pos.z()));

      Vertex_handle infinite_vertex = triangulation.infinite_vertex();

      bool v0_updated = false;
      for (const Cell_handle& ch : find_incident)
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

      // check validity of cells
      for (Cell_handle cit : triangulation.finite_cell_handles())
      {
        if (!is_well_oriented(triangulation, cit))
          return ORIENTATION_PROBLEM;
        if (!triangulation.tds().is_valid(cit, true))
          return C_PROBLEM;
      }

      // check angles
      for (Cell_handle cit : triangulation.finite_cell_handles())
      {
        auto max_cos_after_collapse = max_cos_dihedral_angle(triangulation, cit, false);
        if (      curr_max_cos < max_cos_after_collapse  // angles decreased
         && acceptable_max_cos < max_cos_after_collapse) // && angles go below acceptable bound
          return ANGLE_PROBLEM;
      }

      // check validity of vertices
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

  const Collapse_type collapse_type;

  const Vertex_handle v0_init;
  const Vertex_handle v1_init;

  Vertex_handle vh0;
  Vertex_handle vh1;

  Edge edge;

  bool not_an_edge;

  const Dihedral_angle_cosine acceptable_max_cos{0.995}; // 0.995 cos <=> 5.7 degrees
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

    for (const Cell_handle& ch : cells_to_check)
    {
      if (!ch->has_vertex(v1))
      {
        //check orientation
        std::array<Point, 4> pts = { ch->vertex(0)->point(),
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

    for (const Cell_handle& ch : cells_to_check)
    {
      if (!ch->has_vertex(v0))
      {
        //check orientation
        std::array<Point, 4> pts = { ch->vertex(0)->point(),
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

template<typename Facet, typename Vh>
bool facet_has_edge(const Facet& f, const Vh v0, const Vh v1)
{
  std::array<std::array<int, 2>, 3> edges = {{ {{1,2}}, {{2,3}}, {{3,1}} }};

  for (int i = 0; i < 3; ++i)
  {
    const std::array<int, 2>& ei = edges[i];
    if ( f.first->vertex((f.second + ei[0]) % 4) == v0
      && f.first->vertex((f.second + ei[1]) % 4) == v1)
      return true;
    if ( f.first->vertex((f.second + ei[0]) % 4) == v1
      && f.first->vertex((f.second + ei[1]) % 4) == v0)
      return true;
  }
  return false;
}

template<typename C3t3, typename CellSelector>
bool collapse_preserves_surface_star(const typename C3t3::Edge& edge,
                                     const C3t3& c3t3,
                                     const typename C3t3::Triangulation::Point& new_pos,
                                     const CellSelector& cell_selector)
{
  typedef typename C3t3::Triangulation       Tr;
  typedef typename C3t3::Vertex_handle       Vertex_handle;
  typedef typename C3t3::Facet               Facet;
  typedef typename Tr::Geom_traits::Vector_3 Vector_3;
  typedef typename Tr::Geom_traits::Point_3  Point_3;

  const Tr& tr = c3t3.triangulation();

  const Vertex_handle v0 = edge.first->vertex(edge.second);
  const Vertex_handle v1 = edge.first->vertex(edge.third);
  if (c3t3.in_dimension(v0) != 2 || c3t3.in_dimension(v1) != 2)
    return true;//other cases should not be treated here

  typename Tr::Geom_traits gt = c3t3.triangulation().geom_traits();
  typename Tr::Geom_traits::Construct_opposite_vector_3
    opp = gt.construct_opposite_vector_3_object();
  typename Tr::Geom_traits::Compute_scalar_product_3
    product = gt.compute_scalar_product_3_object();
  typename Tr::Geom_traits::Construct_normal_3
    normal = gt.construct_normal_3_object();

  std::unordered_set<Facet, boost::hash<Facet>> facets;
  tr.finite_incident_facets(v0, std::inserter(facets, facets.end()));
  tr.finite_incident_facets(v1, std::inserter(facets, facets.end()));

// note : checking a 2nd ring of facets does not change the result
//  std::unordered_set<Facet, boost::hash<Facet>> ring2;
//  for (const Facet& f : facets)
//  {
//    for (int i = 1; i < 4; ++i)
//    {
//      Vertex_handle vi = f.first->vertex((f.second + i) % 4);
//      tr.finite_incident_facets(vi, std::inserter(ring2, ring2.end()));
//    }
//  }
//  facets.insert(ring2.begin(), ring2.end());

  Vector_3 reference_normal = CGAL::NULL_VECTOR;
  //Point_3 reference_c;
  for (const Facet& f : facets)
  {
    if (!is_boundary(c3t3, f, cell_selector))
      continue;
    if (facet_has_edge(f, v0, v1))
      continue; //this facet will collapse if collapse happens

    std::array<Point_3, 3> pts = {{ point(f.first->vertex((f.second + 1) % 4)->point()),
                                    point(f.first->vertex((f.second + 2) % 4)->point()),
                                    point(f.first->vertex((f.second + 3) % 4)->point()) }};
    if(f.second % 2 == 0)
      std::swap(pts[0], pts[1]);

    Vector_3 n_before_collapse = normal(pts[0], pts[1], pts[2]);

    const Facet& mf = tr.mirror_facet(f);
    bool do_opp = false;
    if (  c3t3.triangulation().is_infinite(mf.first)
      ||  c3t3.subdomain_index(mf.first) < c3t3.subdomain_index(f.first))
    {
      n_before_collapse = opp(n_before_collapse);
      do_opp = true;
    }

    if (reference_normal == CGAL::NULL_VECTOR)
    {
      //reference_c = CGAL::centroid(pts[0], pts[1], pts[2]);
      reference_normal = n_before_collapse;
    }

    // check after move
    for (int i = 0; i < 3; ++i)
    {
      const Vertex_handle vi = f.first->vertex((f.second + i + 1) % 4);
      if (vi == v0 || vi == v1)
      {
        if (f.second % 2 == 0)
        {
          if(i == 0)      pts[1] = point(new_pos);
          else if(i == 1) pts[0] = point(new_pos);
          else            pts[2] = point(new_pos);
        }
        else
          pts[i] = point(new_pos);
        break;
      }
    }

    Vector_3 n_after_collapse = normal(pts[0], pts[1], pts[2]);
    if(do_opp)
      n_after_collapse = opp(n_after_collapse);

    const double dotref = product(reference_normal, n_after_collapse);
    if(dotref < 0)
      return false;
    const double dot = product(n_before_collapse, n_after_collapse);
    if(dot < 0)
      return false;

//    if (dot * dotref < 0)
//    {
//      std::cout << "\ncollapse edge : ";
//      std::cout << point(v0->point()) << " " << point(v1->point()) << std::endl;
//      const auto vs = c3t3.triangulation().vertices(f);
//      const std::array<Point_3, 3> ps = { {point(vs[0]->point()),
//                                           point(vs[1]->point()),
//                                           point(vs[2]->point())} };
//      std::cout << "facet : ";
//      std::cout << ps[0] << " " << ps[1] << " " << ps[2] << std::endl;
//      const Point_3 c = CGAL::centroid(ps[0], ps[1], ps[2]);
//      std::cout << "n_before_collapse ";
//      std::cout << c << " " << (c + n_before_collapse) << std::endl;
//      std::cout << "n_after_collapse  ";
//      std::cout << c << " " << (c + n_after_collapse) << std::endl;
//      std::cout << "reference_normal  ";
//      std::cout << reference_c << " " << (reference_c + reference_normal) << std::endl;
//      std::cout << std::endl;
//    }
//    if (dotref < 0 || dot < 0)
//      return false;
  }

  return true;
}

template<typename C3t3, typename Sizing, typename CellSelector>
bool are_edge_lengths_valid(const typename C3t3::Edge& edge,
                            const C3t3& c3t3,
                            const Collapse_type& collapse_type,
                            const typename C3t3::Triangulation::Point& new_pos,
                            const Sizing& sizing,
                            const CellSelector& cell_selector)
{
  //SqLengthMap::key_type is Vertex_handle
  //SqLengthMap::value_type is double
  typedef typename C3t3::Triangulation::Geom_traits::FT FT;
  typedef typename C3t3::Edge                           Edge;
  typedef typename C3t3::Vertex_handle                  Vertex_handle;

  const Vertex_handle v0 = edge.first->vertex(edge.second);
  const Vertex_handle v1 = edge.first->vertex(edge.third);

  std::vector<Edge> inc_edges;
  if (collapse_type == TO_V1 || collapse_type == TO_MIDPOINT)
    c3t3.triangulation().finite_incident_edges(v0,
      std::back_inserter(inc_edges));
  if (collapse_type == TO_V0 || collapse_type == TO_MIDPOINT)
    c3t3.triangulation().finite_incident_edges(v1,
      std::back_inserter(inc_edges));

  typename C3t3::Index new_index{};
  int new_dim = -1;
  FT sizing_at_new_pos = FT(0);
  if (collapse_type == TO_V0)
  {
    new_index = c3t3.index(v0);
    new_dim = c3t3.in_dimension(v0);
    sizing_at_new_pos = sizing_at_vertex(v0, sizing, c3t3, cell_selector);
  }
  else if (collapse_type == TO_V1)
  {
    new_index = c3t3.index(v1);
    new_dim = c3t3.in_dimension(v1);
    sizing_at_new_pos = sizing_at_vertex(v1, sizing, c3t3, cell_selector);
  }
  else if (collapse_type == TO_MIDPOINT)
  {
    new_index = max_dimension_index(v0, v1);
    new_dim = (std::max)(v0->in_dimension(), v1->in_dimension());
#ifdef CGAL_AVERAGE_SIZING_AFTER_COLLAPSE
    sizing_at_new_pos = sizing_at_midpoint(edge, new_dim, new_index, sizing, c3t3, cell_selector);
#else
    sizing_at_new_pos = sizing(point(new_pos), new_dim, new_index);
    if(new_dim < 3 && sizing_at_new_pos == FT(0))
      sizing_at_new_pos = max_sizing_in_incident_cells(edge, sizing, c3t3, cell_selector);
#endif
  }
  else
    CGAL_assertion(false);

  std::unordered_map<Vertex_handle, FT> edges_sqlength_after_collapse;
  for (const Edge& ei : inc_edges)
  {
    if (is_outside(ei, c3t3, cell_selector))
      continue;

    Vertex_handle vh = ei.first->vertex(ei.second);
    if (vh == v0 || vh == v1)
      vh = ei.first->vertex(ei.third);
    if (vh == v0 || vh == v1)
      continue;

    FT sqlen = FT(0);
    if (edges_sqlength_after_collapse.find(vh) == edges_sqlength_after_collapse.end())
    {
      sqlen = CGAL::squared_distance(point(new_pos), point(vh->point()));
      edges_sqlength_after_collapse[vh] = sqlen;
    }
    else
      sqlen = edges_sqlength_after_collapse[vh];

    const FT sizing_at_vh = sizing_at_vertex(vh, sizing, c3t3, cell_selector);
    const FT sqhigh
        = CGAL::square(FT(4) / FT(3)) * (std::max)(CGAL::square(sizing_at_vh),
                                                   CGAL::square(sizing_at_new_pos));
    if (sqlen > sqhigh)
      return false;

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
      //
      //}
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

template<typename C3t3, typename CellSelector, typename ShortEdgesBimap>
typename C3t3::Vertex_handle
collapse(const typename C3t3::Cell_handle ch,
         const int to, const int from,
         CellSelector& cell_selector,
         C3t3& c3t3,
         ShortEdgesBimap& short_edges)
{
  typedef typename C3t3::Triangulation Tr;
  typedef typename C3t3::Vertex_handle Vertex_handle;
  typedef typename C3t3::Cell_handle   Cell_handle;
  typedef typename C3t3::Facet         Facet;
  typedef typename Tr::Cell_circulator Cell_circulator;

  Tr& tr = c3t3.triangulation();

  Vertex_handle vkept = ch->vertex(to);
  const Vertex_handle vdeleted = ch->vertex(from);

  //Update the vertex before removing it
  std::vector<Cell_handle> incident_to_vkept;
  tr.incident_cells(vkept, std::back_inserter(incident_to_vkept));

  std::vector<Cell_handle> incident_to_vdeleted;
  tr.incident_cells(vdeleted, std::back_inserter(incident_to_vdeleted));

  boost::container::small_vector<Cell_handle, 30> incident_to_edge;
  Cell_circulator circ = tr.incident_cells(ch, to, from);
  Cell_circulator done = circ;
  do
  {
    for (int i = 0; i < 4; ++i)
    {
      const Vertex_handle vi = circ->vertex(i);
      if (vi != vkept && vi != vdeleted)
      {
        const Facet fi(circ, i);
        if (c3t3.is_in_complex(fi))
          c3t3.remove_from_complex(fi);
      }
    }
    incident_to_edge.push_back(circ);
  }
  while (++circ != done);

  if(c3t3.is_in_complex(ch->vertex(from), ch->vertex(to)))
    c3t3.remove_from_complex(ch->vertex(from), ch->vertex(to));

  std::vector<Cell_handle> cells_to_remove;
  std::unordered_set<Cell_handle> invalid_cells;

  for(const Cell_handle& c : incident_to_edge)
  {
    const int v0_id = c->index(vkept);
    const int v1_id = c->index(vdeleted);

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
  for (const Cell_handle& c : incident_to_vkept)
  {
    if (invalid_cells.find(c) == invalid_cells.end())//valid cell
    {
      if (tr.is_infinite(c))
        infinite_vertex->set_cell(c);
      //else {
      vkept->set_cell(c);
      v0_updated = true;
      //}
    }
  }

  // update complex edges
  for (const Cell_handle& c : incident_to_vdeleted)
  {
    for (const auto& ei : cell_edges(c, tr))
    {
      remove_from_bimap(ei, short_edges);

      const Vertex_handle eiv0 = c->vertex(ei.second);
      const Vertex_handle eiv1 = c->vertex(ei.third);
      if (eiv1 == vdeleted && eiv0 != vkept) //replace eiv1 by vkept
      {
        if (c3t3.is_in_complex(eiv0, eiv1))
        {
          if (!c3t3.is_in_complex(eiv0, vkept))
            c3t3.add_to_complex(eiv0, vkept, c3t3.curve_index(eiv0, eiv1));
          c3t3.remove_from_complex(eiv0, eiv1);
        }
      }
      else if (eiv0 == vdeleted && eiv1 != vkept) //replace eiv0 by vkept
      {
        if (c3t3.is_in_complex(eiv0, eiv1))
        {
          if (!c3t3.is_in_complex(vkept, eiv1))
            c3t3.add_to_complex(vkept, eiv1, c3t3.curve_index(eiv0, eiv1));
          c3t3.remove_from_complex(eiv0, eiv1);
        }
      }
    }
  }

  //Update the vertex before removing it
  for (const Cell_handle& c : incident_to_vdeleted)
  {
    if (invalid_cells.find(c) == invalid_cells.end()) //valid cell
    {
      c->set_vertex(c->index(vdeleted), vkept);

      if (tr.is_infinite(c))
        infinite_vertex->set_cell(c);
      //else {
      if (!v0_updated) {
        vkept->set_cell(c);
        v0_updated = true;
      }
      //}
    }
  }

  if (!v0_updated)
    std::cout << "PB i cell not valid!!!" << std::endl;

  // Delete vertex
  c3t3.triangulation().tds().delete_vertex(vdeleted);

  // Delete cells
  for (Cell_handle cell_to_remove : cells_to_remove)
  {
    // remove cell
    treat_before_delete(cell_to_remove, cell_selector, c3t3);
    c3t3.triangulation().tds().delete_cell(cell_to_remove);
  }

  return vkept;
}


template<typename C3t3, typename CellSelector, typename ShortEdgesBimap>
typename C3t3::Vertex_handle collapse(typename C3t3::Edge& edge,
                                      const Collapse_type& collapse_type,
                                      CellSelector& cell_selector,
                                      C3t3& c3t3,
                                      ShortEdgesBimap& short_edges)
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

    vh = collapse(edge.first, edge.second, edge.third, cell_selector, c3t3, short_edges);
    c3t3.set_dimension(vh, (std::min)(dim_vh0, dim_vh1));
  }
  else //Collapse at vertex
  {
    if (collapse_type == TO_V1)
    {
      vh0->set_point(p1);
      vh = collapse(edge.first, edge.third, edge.second, cell_selector, c3t3, short_edges);
      c3t3.set_dimension(vh, (std::min)(dim_vh0, dim_vh1));
    }
    else //Collapse at v0
    {
      if (collapse_type == TO_V0)
      {
        vh1->set_point(p0);
        vh = collapse(edge.first, edge.second, edge.third, cell_selector, c3t3, short_edges);
        c3t3.set_dimension(vh, (std::min)(dim_vh0, dim_vh1));
      }
      else
        CGAL_assertion(false);
    }
  }
  return vh;
}

template<typename C3t3>
bool is_cells_set_manifold(const C3t3&,
    std::unordered_set<typename C3t3::Cell_handle>& cells)
{
  typedef typename C3t3::Cell_handle Cell_handle;
  typedef typename C3t3::Vertex_handle Vh;
  typedef std::array<Vh, 3> FV;
  typedef std::pair<Vh, Vh> EV;

  std::unordered_map<FV, int, boost::hash<FV>> facets;
  for (Cell_handle c : cells)
  {
    for (int i = 0; i < 4; ++i)
    {
      const FV fvi = make_vertex_array(c->vertex((i + 1) % 4),
        c->vertex((i + 2) % 4),
        c->vertex((i + 3) % 4));
      typename std::unordered_map<FV, int, boost::hash<FV>>::iterator fit = facets.find(fvi);
      if (fit == facets.end())
        facets.insert(std::make_pair(fvi, 1));
      else
        fit->second++;
    }
  }

  std::unordered_map<EV, int, boost::hash<EV>> edges;
  for (const auto& fvv : facets)
  {
    if (fvv.second != 1)
      continue;

    for (int i = 0; i < 3; ++i)
    {
      const EV evi = make_vertex_pair(fvv.first[i], fvv.first[(i + 1) % 3]);
      typename std::unordered_map<EV, int, boost::hash<EV>>::iterator eit = edges.find(evi);
      if (eit == edges.end())
        edges.insert(std::make_pair(evi, 1));
      else
        eit->second++;
    }
  }

  for (const auto& evv : edges)
    if (evv.second != 2)
      return false;

  return true;
}

template<typename C3t3,
         typename Sizing,
         typename CellSelector,
         typename ShortEdgesBimap,
         typename Visitor>
typename C3t3::Vertex_handle collapse_edge(typename C3t3::Edge& edge,
    C3t3& c3t3,
    const Sizing& sizing,
    const bool /* protect_boundaries */,
    CellSelector cell_selector,
    ShortEdgesBimap& short_edges,
    Visitor& )
{
  typedef typename C3t3::Triangulation   Tr;
  typedef typename Tr::Point             Point;
  typedef typename Tr::Vertex_handle     Vertex_handle;
  typedef typename Tr::Cell_handle       Cell_handle;

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
    {
#ifdef CGAL_DEBUG_TET_REMESHING_IN_PLUGIN
      if (in_cx)
        nb_invalid_collapse++;
#endif
      return Vertex_handle();
    }
  }

  if (are_edge_lengths_valid(edge, c3t3, collapse_type, new_pos, sizing, cell_selector)
    && collapse_preserves_surface_star(edge, c3t3, new_pos, cell_selector))
  {
    CGAL_expensive_assertion(c3t3.triangulation().tds().is_edge(
                       edge.first->vertex(edge.second),
                       edge.first->vertex(edge.third)));

    Vertex_handle v0_init = edge.first->vertex(edge.second);
    Vertex_handle v1_init = edge.first->vertex(edge.third);

    std::unordered_set<Cell_handle> cells_to_insert;
    c3t3.triangulation().finite_incident_cells(v0_init,
      std::inserter(cells_to_insert, cells_to_insert.end()));
    c3t3.triangulation().finite_incident_cells(v1_init,
      std::inserter(cells_to_insert, cells_to_insert.end()));

    if(!is_cells_set_manifold(c3t3, cells_to_insert))
      return Vertex_handle();

    CollapseTriangulation<C3t3> local_tri(edge, cells_to_insert, collapse_type);

    Result_type res = local_tri.collapse();
    if (res == VALID)
    {
#ifdef CGAL_DEBUG_TET_REMESHING_IN_PLUGIN
      if (in_cx)
        nb_valid_collapse++;
#endif
      return collapse(edge, collapse_type, cell_selector, c3t3, short_edges);
    }
  }
#ifdef CGAL_DEBUG_TET_REMESHING_IN_PLUGIN
  else if (in_cx)
    nb_invalid_lengths++;
#endif
  return Vertex_handle();
}

template<typename C3T3, typename CellSelector>
auto can_be_collapsed(const typename C3T3::Edge& e,
                      const C3T3& c3t3,
                      const bool protect_boundaries,
                      CellSelector cell_selector)
{
  struct Collapsible
  {
    bool can_be_collapsed;
    bool on_boundary;
  };

  const bool boundary = c3t3.is_in_complex(e)
                     || is_boundary(c3t3, e, cell_selector);

  if (protect_boundaries)
  {
    if (boundary)
      return Collapsible{false, boundary};

    return Collapsible{ is_internal(e, c3t3, cell_selector), boundary};
  }
  else
  {
    return Collapsible{ is_selected(e, c3t3.triangulation(), cell_selector),
                        boundary };
  }
}

template<typename C3T3,
         typename Sizing,
         typename CellSelector,
         typename Visitor>
void collapse_short_edges(C3T3& c3t3,
                          const Sizing& sizing,
                          const bool protect_boundaries,
                          CellSelector cell_selector,
                          Visitor& visitor)
{
  typedef typename C3T3::Triangulation       T3;
  typedef typename T3::Edge                  Edge;
  typedef typename T3::Vertex_handle         Vertex_handle;

  typedef typename T3::Geom_traits::FT FT;
  typedef boost::bimap<
        boost::bimaps::set_of<Edge, Compare_edges<Edge> >,
        boost::bimaps::multiset_of<FT, std::less<FT> > >  Boost_bimap;
  typedef typename Boost_bimap::value_type            short_edge;

  T3& tr = c3t3.triangulation();

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
  std::cout << "Collapse short edges...";
  std::cout.flush();
  std::size_t nb_collapses = 0;
  CGAL::Real_timer timer;
  timer.start();
#endif

  //collect long edges
  Boost_bimap short_edges;
  for (const Edge& e : tr.finite_edges())
  {
    auto [collapsible, boundary] = can_be_collapsed(e, c3t3, protect_boundaries, cell_selector);
    if (!collapsible)
      continue;

    const auto sqlen = is_too_short(e, boundary, sizing, c3t3, cell_selector);
    if(sqlen != std::nullopt)
      short_edges.insert(short_edge(e, sqlen.value()));
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
    Edge e = eit->second;

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE_PROGRESS
    FT sqlen = eit->first;
    std::cout << "\rCollapse... (" << short_edges.left.size() << " short edges, ";
    std::cout << std::sqrt(sqlen) << ", ";
    std::cout << nb_collapses << " collapses)";
    std::cout.flush();
#endif

    short_edges.right.erase(eit);

    CGAL_expensive_assertion_code(const bool bd = is_boundary_edge(e));
    CGAL_expensive_assertion(!!is_too_short(e, bd, sizing, c3t3, cell_selector));
    CGAL_expensive_assertion(can_be_collapsed(e, c3t3, protect_boundaries, cell_selector));

#ifdef CGAL_TETRAHEDRAL_REMESHING_DEBUG
    const auto p1 = e.first->vertex(e.second)->point();
    const auto p2 = e.first->vertex(e.third)->point();
#endif

    Vertex_handle vh = collapse_edge(e, c3t3, sizing,
                                     protect_boundaries, cell_selector,
                                     short_edges,
                                     visitor);
    if (vh != Vertex_handle())
    {
      std::vector<Edge> incident_short;
      c3t3.triangulation().finite_incident_edges(vh,
          std::back_inserter(incident_short));
      for (const Edge& eshort : incident_short)
      {
        const auto [collapsible, boundary]
          = can_be_collapsed(eshort, c3t3, protect_boundaries, cell_selector);
        if (!collapsible)
          continue;

        const auto sqlen = is_too_short(eshort, boundary, sizing, c3t3, cell_selector);
        update_bimap(eshort, short_edges, sqlen);
      }

#ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
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
  timer.stop();
  std::cout << " done (" << nb_collapses << " collapses, in "
            << timer.time() << " seconds)." << std::endl;
#endif
}
}
}
}

#endif // CGAL_INTERNAL_COLLAPSE_SHORT_EDGES_H
