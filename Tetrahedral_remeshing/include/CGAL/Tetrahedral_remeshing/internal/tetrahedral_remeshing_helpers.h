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

#ifndef CGAL_INTERNAL_TET_REMESHING_HELPERS_H
#define CGAL_INTERNAL_TET_REMESHING_HELPERS_H

#include <CGAL/license/Tetrahedral_remeshing.h>

#include <utility>
#include <array>
#include <iterator>

#include <CGAL/Point_3.h>
#include <CGAL/Weighted_point_3.h>
#include <CGAL/Vector_3.h>
#include <CGAL/utility.h>

#include <CGAL/IO/File_binary_mesh_3.h>

#include <boost/container/flat_set.hpp>
#include <boost/container/small_vector.hpp>

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
  CGAL_assertion(i < 4 && j < 3);
  if(i < 4 && j < 3)
    return indices_table[i][j];
  else
    return -1;
}

template<typename Gt, typename Point>
typename Gt::FT dihedral_angle(const Point& p,
                               const Point& q,
                               const Point& r,
                               const Point& s,
                               const Gt& gt)
{
  return gt.compute_approximate_dihedral_angle_3_object()(p, q, r, s);
}

template<typename Point, typename Geom_traits>
typename Geom_traits::FT min_dihedral_angle(const Point& p,
                                            const Point& q,
                                            const Point& r,
                                            const Point& s,
                                            const Geom_traits& gt)
{
  typedef typename Geom_traits::FT FT;
  FT a = CGAL::abs(dihedral_angle(p, q, r, s, gt));
  FT min_dh = a;

  a = CGAL::abs(dihedral_angle(p, r, s, q, gt));
  min_dh = (std::min)(a, min_dh);

  a = CGAL::abs(dihedral_angle(p, s, q, r, gt));
  min_dh = (std::min)(a, min_dh);

  a = CGAL::abs(dihedral_angle(q, r, p, s, gt));
  min_dh = (std::min)(a, min_dh);

  a = CGAL::abs(dihedral_angle(q, s, r, p, gt));
  min_dh = (std::min)(a, min_dh);

  a = CGAL::abs(dihedral_angle(r, s, p, q, gt));
  min_dh = (std::min)(a, min_dh);

  return min_dh;
}

template<typename Tr>
typename Tr::Geom_traits::FT min_dihedral_angle(const Tr& tr,
                                                const typename Tr::Vertex_handle v0,
                                                const typename Tr::Vertex_handle v1,
                                                const typename Tr::Vertex_handle v2,
                                                const typename Tr::Vertex_handle v3)
{
  return min_dihedral_angle(point(v0->point()),
                            point(v1->point()),
                            point(v2->point()),
                            point(v3->point()),
                            tr.geom_traits());
}

template<typename Tr>
typename Tr::Geom_traits::FT min_dihedral_angle(const Tr& tr,
                                                const typename Tr::Cell_handle c)
{
  return min_dihedral_angle(tr,
                            c->vertex(0),
                            c->vertex(1),
                            c->vertex(2),
                            c->vertex(3));
}

struct Dihedral_angle_cosine
{
  CGAL::Sign m_sgn;
  double m_sq_num;
  double m_sq_den;

  Dihedral_angle_cosine(const CGAL::Sign& sgn, const double& sq_num, const double& sq_den)
    : m_sgn(sgn)
    , m_sq_num(sq_num)
    , m_sq_den(sq_den)
  {}

  bool is_one() const
  {
    return m_sgn == CGAL::POSITIVE && m_sq_num == m_sq_den;
  }
  double signed_square_value() const
  {
    switch(m_sgn)
    {
    case CGAL::POSITIVE:
      return m_sq_num / m_sq_den;
    case CGAL::ZERO:
      return 0.;
    default:
      CGAL_assertion(m_sgn == CGAL::NEGATIVE);
      return -1. * m_sq_num / m_sq_den;
    };
  }

  friend bool operator<(const Dihedral_angle_cosine& l,
                        const Dihedral_angle_cosine& r)
  {
    //if numerators have different signs
    if (l.m_sgn == CGAL::NEGATIVE && r.m_sgn != CGAL::NEGATIVE)
      return true;

    else if (l.m_sgn == CGAL::POSITIVE && r.m_sgn != CGAL::POSITIVE)
      return false;

    else if (l.m_sgn == CGAL::ZERO)
      return (r.m_sgn == CGAL::POSITIVE);

    //else numerators have the same sign
    else if (l.m_sgn == CGAL::POSITIVE) //both angles are in [0; PI/2[
    {
      CGAL_assertion(r.m_sgn == CGAL::POSITIVE);

      return  (l.m_sq_num * r.m_sq_den < r.m_sq_num* l.m_sq_den);
    }
    else //both angles are in [PI/2; PI]
    {
      CGAL_assertion(l.m_sgn != CGAL::POSITIVE);
      CGAL_assertion(r.m_sgn != CGAL::POSITIVE);

      return (l.m_sq_num * r.m_sq_den >= r.m_sq_num* l.m_sq_den);
    }
  }

  friend bool operator<=(const Dihedral_angle_cosine& l,
                         const Dihedral_angle_cosine& r)
  {
    if(l < r)
      return true;
    else
      return l.m_sgn == r.m_sgn
         &&  l.m_sq_num * r.m_sq_den == r.m_sq_num * l.m_sq_den;
  }
};

template<typename Gt>
Dihedral_angle_cosine cos_dihedral_angle(const typename Gt::Point_3& i,
                                         const typename Gt::Point_3& j,
                                         const typename Gt::Point_3& k,
                                         const typename Gt::Point_3& l,
                                         const Gt& gt)
{
  CGAL_assertion(CGAL::orientation(i, j, k, l) != CGAL::NEGATIVE);

  typename Gt::Construct_vector_3 vector = gt.construct_vector_3_object();
  typename Gt::Construct_cross_product_vector_3 cross_product =
    gt.construct_cross_product_vector_3_object();
  typename Gt::Compute_scalar_product_3 scalar_product =
    gt.compute_scalar_product_3_object();

  typedef typename Gt::FT FT;
  typedef typename Gt::Vector_3 Vector_3;
  const Vector_3 ij = vector(i, j);
  const Vector_3 ik = vector(i, k);
  const Vector_3 il = vector(i, l);

  const Vector_3 ijik = cross_product(ij, ik);
  if(CGAL::NULL_VECTOR == ijik)
    return Dihedral_angle_cosine(CGAL::POSITIVE, 1.,1.);

  const Vector_3 ilij = cross_product(il, ij);
  if (CGAL::NULL_VECTOR == ilij)
    return Dihedral_angle_cosine(CGAL::POSITIVE, 1.,1.);

  const FT num = scalar_product(ijik, ilij);
  if(num == 0.)
    return Dihedral_angle_cosine(CGAL::ZERO, 0.,1.);

  const double sqden = CGAL::to_double(
    scalar_product(ijik, ijik) * scalar_product(ilij, ilij));

  return Dihedral_angle_cosine(CGAL::sign(num), CGAL::square(num), sqden);
}

template<typename Point, typename Geom_traits>
Dihedral_angle_cosine max_cos_dihedral_angle(const Point& p,
                                             const Point& q,
                                             const Point& r,
                                             const Point& s,
                                             const Geom_traits& gt)
{
  Dihedral_angle_cosine a = cos_dihedral_angle(p, q, r, s, gt);
  Dihedral_angle_cosine max_cos_dh = a;
  if(max_cos_dh.is_one()) return max_cos_dh;

  a = cos_dihedral_angle(p, r, s, q, gt);
  if(max_cos_dh < a) max_cos_dh = a;
  if (max_cos_dh.is_one()) return max_cos_dh;

  a = cos_dihedral_angle(p, s, q, r, gt);
  if (max_cos_dh < a) max_cos_dh = a;
  if (max_cos_dh.is_one()) return max_cos_dh;

  a = cos_dihedral_angle(q, r, p, s, gt);
  if (max_cos_dh < a) max_cos_dh = a;
  if (max_cos_dh.is_one()) return max_cos_dh;

  a = cos_dihedral_angle(q, s, r, p, gt);
  if (max_cos_dh < a) max_cos_dh = a;
  if (max_cos_dh.is_one()) return max_cos_dh;

  a = cos_dihedral_angle(r, s, p, q, gt);
  if (max_cos_dh < a) max_cos_dh = a;

  return max_cos_dh;
}

template<typename Tr>
Dihedral_angle_cosine max_cos_dihedral_angle(const Tr& tr,
                                             const typename Tr::Vertex_handle v0,
                                             const typename Tr::Vertex_handle v1,
                                             const typename Tr::Vertex_handle v2,
                                             const typename Tr::Vertex_handle v3)
{
  return max_cos_dihedral_angle(point(v0->point()),
                                point(v1->point()),
                                point(v2->point()),
                                point(v3->point()),
                                tr.geom_traits());
}

template<typename Tr>
Dihedral_angle_cosine max_cos_dihedral_angle(const Tr& tr,
                                             const typename Tr::Cell_handle c)
{
  if (c->is_cache_valid())
    return Dihedral_angle_cosine(CGAL::sign(c->sliver_value()),
                                 CGAL::abs(c->sliver_value()), 1.);

  Dihedral_angle_cosine cos_dh = max_cos_dihedral_angle(tr,
                                                        c->vertex(0),
                                                        c->vertex(1),
                                                        c->vertex(2),
                                                        c->vertex(3));
  c->set_sliver_value(cos_dh.signed_square_value());
  return cos_dh;
}

template<typename C3t3>
bool is_peelable(const C3t3& c3t3,
                 const typename C3t3::Cell_handle ch,
                 std::array<bool, 4>& facets_on_surface)
{
  typedef typename C3t3::Triangulation::Geom_traits::FT FT;
  typedef typename C3t3::Facet                          Facet;

  if(!c3t3.is_in_complex(ch))
    return false;

  bool on_surface = false;
  for (int i = 0; i < 4; ++i)
  {
    facets_on_surface[i] = !c3t3.is_in_complex(ch->neighbor(i));
    on_surface = on_surface || facets_on_surface[i];
  }
  if(!on_surface)
    return false;

  FT area_on_surface = 0.;
  FT area_inside = 0.;
  for (int i = 0; i < 4; ++i)
  {
    Facet f(ch, i);
    const FT facet_area = CGAL::approximate_sqrt(c3t3.triangulation().triangle(f).squared_area());
    if(facets_on_surface[i])
      area_on_surface += facet_area;
    else
      area_inside += facet_area;
  }

  return (area_inside < 1.5 * area_on_surface);
}

template<typename Tr>
typename Tr::Geom_traits::Vector_3 facet_normal(const Tr& tr,
                                                const typename Tr::Facet& f)
{
  const typename Tr::Geom_traits gt = tr.geom_traits();
  typename Tr::Geom_traits::Construct_normal_3 cn
    = gt.construct_normal_3_object();
  return cn(point(f.first->vertex((f.second + 1) % 4)->point()),
            point(f.first->vertex((f.second + 2) % 4)->point()),
            point(f.first->vertex((f.second + 3) % 4)->point()));
}

template<typename Vh>
std::pair<Vh, Vh> make_vertex_pair(const Vh v1, const Vh v2)
{
  if (v2 < v1) return std::make_pair(v2, v1);
  else         return std::make_pair(v1, v2);
}

template<typename Tr>
std::pair<typename Tr::Vertex_handle, typename Tr::Vertex_handle>
make_vertex_pair(const typename Tr::Edge& e)
{
  typedef typename Tr::Vertex_handle Vertex_handle;
  Vertex_handle v1 = e.first->vertex(e.second);
  Vertex_handle v2 = e.first->vertex(e.third);
  return make_vertex_pair(v1, v2);
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

template<typename Vh>
std::array<Vh, 3> make_vertex_array(const Vh vh0, const Vh vh1, const Vh vh2)
{
  std::array<Vh, 3> ft = { {vh0, vh1, vh2} };
  if (ft[1] < ft[0]) std::swap(ft[0], ft[1]);
  if (ft[2] < ft[1]) std::swap(ft[1], ft[2]);
  if (ft[1] < ft[0]) std::swap(ft[0], ft[1]);
  return ft;
}

template<typename Facet>
Facet canonical_facet(const Facet& f)
{
  const typename Facet::first_type c = f.first;
  const int i = f.second;
  const typename Facet::first_type c2 = c->neighbor(i);
  return (c2 < c) ? std::make_pair(c2, c2->index(c)) : std::make_pair(c, i);
}

template<typename VertexHandle>
bool is_on_feature(const VertexHandle v)
{
  return (v->in_dimension() == 1 || v->in_dimension() == 0);
}

template<typename Tr>
bool is_well_oriented(const Tr& tr,
                      const typename Tr::Vertex_handle v0,
                      const typename Tr::Vertex_handle v1,
                      const typename Tr::Vertex_handle v2,
                      const typename Tr::Vertex_handle v3)
{
  return CGAL::POSITIVE == tr.geom_traits().orientation_3_object()(
           point(v0->point()),
           point(v1->point()),
           point(v2->point()),
           point(v3->point()));
}

template<typename Tr>
bool is_well_oriented(const Tr& tr, const typename Tr::Cell_handle ch)
{
  return is_well_oriented(tr, ch->vertex(0), ch->vertex(1),
                          ch->vertex(2), ch->vertex(3));
}

template<typename C3T3, typename CellSelector>
bool is_boundary(const C3T3& c3t3,
                 const typename C3T3::Facet& f,
                 const CellSelector& cell_selector)
{
  return c3t3.is_in_complex(f)
         || cell_selector(f.first) != cell_selector(f.first->neighbor(f.second));
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

  do
  {
    const Facet& f = *fcirc;
    if (is_boundary(c3t3, f, cell_selector))
      return true;
  }
  while (++fcirc != fend);

  return false;
}

template<typename C3t3, typename CellSelector>
bool is_boundary_edge(const typename C3t3::Vertex_handle& v0,
                      const typename C3t3::Vertex_handle& v1,
                      const C3t3& c3t3,
                      const CellSelector& cell_selector)
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

  for(const Facet& f : facets)
  {
    if (c3t3.is_in_complex(f))
      return true;
    if (cell_selector(f.first) ^ cell_selector(f.first->neighbor(f.second)))
      return true;
  }
  return false;
}

template<typename C3t3>
typename C3t3::Surface_patch_index surface_patch_index(const typename C3t3::Vertex_handle v,
    const C3t3& c3t3)
{
  typedef typename C3t3::Surface_patch_index Surface_patch_index;
  typedef typename C3t3::Facet Facet;
  std::vector<Facet> facets;
  c3t3.triangulation().incident_facets(v, std::back_inserter(facets));

  for(const Facet& f : facets)
  {
    if (c3t3.is_in_complex(f))
      return c3t3.surface_patch_index(f);
  }
  return Surface_patch_index();
}

template<typename C3t3>
void set_index(typename C3t3::Vertex_handle v, const C3t3& c3t3)
{
  switch (v->in_dimension())
  {
  case 3:
    v->set_index(v->cell()->subdomain_index());
    break;
  case 2:
    v->set_index(surface_patch_index(v, c3t3));
    break;
  case 1:
    v->set_index(typename C3t3::Curve_index(1));
    break;
  case 0:
    v->set_index(boost::get<typename C3t3::Corner_index>(v->index()));
    break;
  default:
    CGAL_assertion(false);
  }
}

template<typename C3t3>
bool is_edge_in_complex(const typename C3t3::Vertex_handle& v0,
                        const typename C3t3::Vertex_handle& v1,
                        const C3t3& c3t3)
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

template<typename C3t3, typename OutputIterator>
OutputIterator incident_surface_patches(const typename C3t3::Edge& e,
                                        const C3t3& c3t3,
                                        OutputIterator oit)
{
  typedef typename C3t3::Triangulation::Facet_circulator Facet_circulator;
  typedef typename C3t3::Triangulation::Facet Facet;

  Facet_circulator circ = c3t3.triangulation().incident_facets(e);
  Facet_circulator end = circ;
  do
  {
    const Facet& f = *circ;
    if(c3t3.is_in_complex(f))
      *oit++ = c3t3.surface_patch_index(f);
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

template <typename C3t3>
std::size_t nb_incident_surface_patches(const typename C3t3::Edge& e,
                                        const C3t3& c3t3)
{
  typedef typename C3t3::Surface_patch_index Surface_patch_index;

  boost::unordered_set<Surface_patch_index> indices;
  incident_surface_patches(e, c3t3, std::inserter(indices, indices.begin()));

  return indices.size();
}

template<typename C3t3>
std::size_t nb_incident_complex_edges(const typename C3t3::Vertex_handle v,
                                      const C3t3& c3t3)
{
  typedef typename C3t3::Edge Edge;
  boost::unordered_set<Edge> edges;
  c3t3.triangulation().finite_incident_edges(v, std::inserter(edges, edges.begin()));

  std::size_t count = 0;
  for (const Edge& e : edges)
  {
    if (c3t3.is_in_complex(e))
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

  if (c3t3.number_of_corners() > 0)
  {
    return c3t3.is_in_complex(v);
  }
  else if (nb_incident_subdomains(v, c3t3) > 3)
  {
    std::vector<Edge> edges;
    c3t3.triangulation().finite_incident_edges(v, std::back_inserter(edges));

    int feature_count = 0;
    for(const Edge& ei : edges)
    {
      if (c3t3.is_in_complex(ei))
      {
        feature_count++;
        if (feature_count >= 3)
          return true;
      }
    }
  }
  return false;
}

/**
* returns true iff `v` is on the outer hull of c3t3.triangulation()
* i.e. finite and incident to at least one infinite cell
*/
template<typename C3t3>
bool is_on_convex_hull(const typename C3t3::Vertex_handle v,
                       const C3t3& c3t3)
{
  if (v == c3t3.triangulation().infinite_vertex())
    return true;

  //on hull == incident to infinite cell
  typedef typename C3t3::Triangulation::Cell_handle Cell_handle;

  std::vector<Cell_handle> cells;
  c3t3.triangulation().incident_cells(v, std::back_inserter(cells));
  for (Cell_handle ci : cells)
  {
    if (c3t3.triangulation().is_infinite(ci))
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
bool is_on_convex_hull(const typename C3t3::Edge & edge,
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

template<typename C3t3, typename CellSelector>
bool is_outside(const typename C3t3::Edge & edge,
                const C3t3& c3t3,
                CellSelector cell_selector)
{
  typedef typename C3t3::Triangulation::Cell_circulator Cell_circulator;
  Cell_circulator circ = c3t3.triangulation().incident_cells(edge);
  Cell_circulator done = circ;
  do
  {
    // is cell in complex?
    if (c3t3.is_in_complex(circ))
      return false;
    // does circ belong to the selection?
    if (cell_selector(circ))
      return false;

    ++circ;
  } while (circ != done);

  return true; //all incident cells are outside or infinite
}

template<typename C3t3, typename CellSelector>
bool is_selected(const typename C3t3::Vertex_handle v,
                 const C3t3& c3t3,
                 CellSelector cell_selector)
{
  typedef typename C3t3::Triangulation::Cell_handle Cell_handle;

  std::vector<Cell_handle> cells;
  c3t3.triangulation().incident_cells(v, std::back_inserter(cells));

  for(Cell_handle c : cells)
  {
    if (cell_selector(c))
      return true;
  }
  return false;
}

template<typename C3t3, typename CellSelector>
bool is_internal(const typename C3t3::Edge& edge,
                 const C3t3& c3t3,
                 CellSelector cell_selector)
{
  const typename C3t3::Vertex_handle vs = edge.first->vertex(edge.second);
  const typename C3t3::Vertex_handle vt = edge.first->vertex(edge.third);

  typedef typename C3t3::Triangulation::Cell_circulator Cell_circulator;
  Cell_circulator circ = c3t3.triangulation().incident_cells(edge);
  Cell_circulator done = circ;

  const typename C3t3::Subdomain_index si = circ->subdomain_index();
  do
  {
    if (c3t3.triangulation().is_infinite(circ))
      return false;
    if (si != circ->subdomain_index())
      return false;
    if (!cell_selector(circ))
      return false;
    if (c3t3.is_in_complex(
          circ,
          CGAL::Triangulation_utils_3::next_around_edge(circ->index(vs), circ->index(vt))))
      return false;
  } while (++circ != done);

  return true;
}

template<typename Gt>
void normalize(typename Gt::Vector_3& v, const Gt& gt)
{
  typedef typename Gt::FT       FT;

  const FT norm = CGAL::approximate_sqrt(gt.compute_squared_length_3_object()(v));
  if (norm != FT(0))
    v = gt.construct_divided_vector_3_object()(v, norm);
}

template<typename Facet, typename Gt>
typename Gt::Vector_3 normal(const Facet& f, const Gt& gt)
{
  typedef typename Gt::Vector_3 Vector_3;
  typedef typename Gt::Point_3  Point;
  typedef typename Gt::FT       FT;

  Point p0 = point(f.first->vertex((f.second + 1) % 4)->point());
  Point p1 = point(f.first->vertex((f.second + 2) % 4)->point());
  const Point& p2 = point(f.first->vertex((f.second + 3) % 4)->point());

  if (f.second % 2 == 0)//equivalent to the commented orientation test
    std::swap(p0, p1);

  Vector_3 n = gt.construct_cross_product_vector_3_object()(
                 gt.construct_vector_3_object()(p1, p2),
                 gt.construct_vector_3_object()(p1, p0));

  //cross-product(AB, AC)'s norm is the area of the parallelogram
  //formed by these 2 vectors.
  //the triangle's area is half of it
  return gt.construct_scaled_vector_3_object()(n, FT(1) / FT(2));
}

template<typename C3t3, typename CellSelector, typename OutputIterator>
OutputIterator get_internal_edges(const C3t3& c3t3,
                                  CellSelector cell_selector,
                                  OutputIterator oit)/*holds Edges*/
{
  for (typename C3t3::Triangulation::Finite_edges_iterator
       eit = c3t3.triangulation().finite_edges_begin();
       eit != c3t3.triangulation().finite_edges_end();
       ++eit)
  {
    const typename C3t3::Edge& e = *eit;
    if (is_internal(e, c3t3, cell_selector))
    {
      *oit++ = make_vertex_pair<typename C3t3::Triangulation>(e);
    }
  }
  return oit;
}

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
      if (is_boundary(c3t3, Edge(circ, (i + 1) % 4, (i + 2) % 4), cell_selector)
          && is_boundary(c3t3, Edge(circ, (i + 2) % 4, (i + 3) % 4), cell_selector)
          && is_boundary(c3t3, Edge(circ, (i + 3) % 4, (i + 1) % 4), cell_selector))
        return false;
    }
  } while (++fcirc != fdone);

  return true;
}

template<typename C3t3>
Subdomain_relation compare_subdomains(const typename C3t3::Vertex_handle v0,
                                      const typename C3t3::Vertex_handle v1,
                                      const C3t3& c3t3)
{
  typedef typename C3t3::Subdomain_index Subdomain_index;
  typedef boost::container::flat_set<Subdomain_index,
    std::less<Subdomain_index>,
    boost::container::small_vector<Subdomain_index, 30> > Set_of_subdomains;

  Set_of_subdomains subdomains_v0;
  incident_subdomains(v0, c3t3,
    std::inserter(subdomains_v0, subdomains_v0.begin()));

  Set_of_subdomains subdomains_v1;
  incident_subdomains(v1, c3t3,
    std::inserter(subdomains_v1, subdomains_v1.begin()));

  if (subdomains_v0.size() == subdomains_v1.size())
  {
    if(std::equal(subdomains_v0.begin(), subdomains_v0.end(), subdomains_v1.begin()))
      return EQUAL;
    else
      return DIFFERENT;
  }
  else
  {
    boost::container::small_vector<Subdomain_index, 30>
      intersection((std::min)(subdomains_v0.size(), subdomains_v1.size()), -1);
    typename boost::container::small_vector<Subdomain_index, 30>::iterator
    end_it = std::set_intersection(subdomains_v0.begin(), subdomains_v0.end(),
                                   subdomains_v1.begin(), subdomains_v1.end(),
                                   intersection.begin());
    std::ptrdiff_t intersection_size =
      std::distance(intersection.begin(), end_it);

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
                   const CellSelector& cell_selector)
{
  typedef typename C3t3::Vertex_handle Vertex_handle;

  update_v0 = false;
  update_v1 = false;

  const Vertex_handle v0 = edge.first->vertex(edge.second);
  const Vertex_handle v1 = edge.first->vertex(edge.third);

  const int dim0 = c3t3.in_dimension(v0);
  const int dim1 = c3t3.in_dimension(v1);

  if (dim0 == 3)
  {
    CGAL_assertion(!is_on_convex_hull(v0, c3t3));
    update_v0 = true;
    if (dim1 == 3)
    {
      CGAL_assertion(!is_on_convex_hull(v1, c3t3));
      update_v1 = true;
      return;
    }
    else // dim1 is 2, 1, or 0
      return;
  }
  else if (dim1 == 3)
  {
    update_v1 = true;
    return;
  }

  // from now on, all cases lie on surfaces, or between surfaces
  CGAL_assertion(dim0 != 3 && dim1 != 3);

  //feature edges and feature vertices
  if (dim0 < 2 || dim1 < 2)
  {
    if (!topology_test(edge, c3t3, cell_selector))
    {
#ifdef CGAL_DEBUG_TET_REMESHING_IN_PLUGIN
      nb_topology_test++;
#endif
      return;
    }

    if (c3t3.is_in_complex(edge))
    {
      const std::size_t nb_si_v0 = nb_incident_subdomains(v0, c3t3);
      const std::size_t nb_si_v1 = nb_incident_subdomains(v1, c3t3);

      if (nb_si_v0 > nb_si_v1) {
        if (!c3t3.is_in_complex(v1))
          update_v1 = true;
      }
      else if (nb_si_v1 > nb_si_v0) {
        if (!c3t3.is_in_complex(v0))
          update_v0 = true;
      }
      else {
        if (!c3t3.is_in_complex(v0))
          update_v0 = true;
        if (!c3t3.is_in_complex(v1))
          update_v1 = true;
      }
    }
    else
    {
      if (dim0 == 2 && is_boundary_edge(v0, v1, c3t3, cell_selector))
      {
        update_v0 = true;
        return;
      }
      else if(dim1 == 2 && is_boundary_edge(v0, v1, c3t3, cell_selector))
      {
        update_v1 = true;
        return;
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
          const bool v0_on_feature = is_on_feature(v0);
          const bool v1_on_feature = is_on_feature(v1);

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
  }
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

  for(typename Bimap::left_const_reference it : edges.left)
  {
    ofs << "2 " << point(it.first.first->point())
        << " " << point(it.first.second->point()) << std::endl;
  }
  ofs.close();
}

template<typename Facet, typename OutputStream>
void dump_facet(const Facet& f, OutputStream& os)
{
  os << "4 ";
  os << point(f.first->vertex((f.second + 1) % 4)->point()) << " "
     << point(f.first->vertex((f.second + 2) % 4)->point()) << " "
     << point(f.first->vertex((f.second + 3) % 4)->point()) << " "
     << point(f.first->vertex((f.second + 1) % 4)->point());
  os << std::endl;
}

template<typename FacetRange>
void dump_facets(const FacetRange& facets, const char* filename)
{
  std::ofstream os(filename);
  for (typename FacetRange::value_type f : facets)
  {
    dump_facet(f, os);
  }
  os.close();
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

template<typename C3t3>
void check_surface_patch_indices(const C3t3& c3t3)
{
  typedef typename C3t3::Vertex_handle Vertex_handle;
  for (Vertex_handle v : c3t3.triangulation().finite_vertex_handles())
  {
    if (v->in_dimension() != 2)
      continue;
    CGAL_assertion(surface_patch_index(v, c3t3) != typename C3t3::Surface_patch_index());
  }
}

template<typename Tr>
bool are_cell_orientations_valid(const Tr& tr)
{
  typedef typename Tr::Geom_traits::Point_3 Point_3;
  typedef typename Tr::Facet                Facet;

  std::set<Facet> facets;
  for (const typename Tr::Cell_handle ch : tr.finite_cell_handles())
  {
    const Point_3& p0 = point(ch->vertex(0)->point());
    const Point_3& p1 = point(ch->vertex(1)->point());
    const Point_3& p2 = point(ch->vertex(2)->point());
    const Point_3& p3 = point(ch->vertex(3)->point());

    const CGAL::Orientation o = CGAL::orientation(p0, p1, p2, p3);
    if (o != CGAL::POSITIVE)
    {
      facets.insert(canonical_facet(Facet(ch, 0)));
      facets.insert(canonical_facet(Facet(ch, 1)));
      facets.insert(canonical_facet(Facet(ch, 2)));
      facets.insert(canonical_facet(Facet(ch, 3)));
    }
  }
  if (!facets.empty())
  {
    std::cerr << "Warning : there are inverted cells!\n"
              << "\tSee cells_with_negative_volume.polylines.txt" << std::endl;
    dump_facets(facets, "cells_with_negative_volume.polylines.txt");
  }
  return facets.empty();
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
    ofs << point(vit->second->point()) << std::endl;
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

template<typename CellRange, typename Tr>
void dump_cells_off(const CellRange& cells, const Tr& /*tr*/, const char* filename)
{
  typedef typename Tr::Vertex_handle              Vertex_handle;
  typedef typename Tr::Cell_handle                Cell_handle;
  typedef boost::bimap<Vertex_handle, int>        Bimap_t;
  typedef typename Bimap_t::left_map::value_type  value_type;

  Bimap_t vertices;
  int index = 0;
  boost::unordered_set<std::array<Vertex_handle, 3> > facets;

  for (Cell_handle c : cells)
  {
    //collect vertices
    for (int i = 0; i < 4; ++i)
    {
      Vertex_handle vi = c->vertex(i);
      if (vertices.left.find(c->vertex(i)) == vertices.left.end())
        vertices.left.insert(value_type(vi, index++));
    }
    //collect facets
    for (int i = 0; i < 4; ++i)
    {
      //if (tr.is_infinite(c->neighbor(i)))
      {
        std::array<Vertex_handle, 3> fi = make_vertex_array(c->vertex((i + 1) % 4),
          c->vertex((i + 2) % 4),
          c->vertex((i + 3) % 4));
          facets.insert(fi);
      }
    }
  }

  //write header
  std::ofstream ofs(filename);
  ofs.precision(17);
  ofs << "OFF" << std::endl;
  ofs << vertices.size() << " " << facets.size() << " 0" << std::endl << std::endl;

  for(const typename Bimap_t::right_map::value_type& v : vertices.right)
    ofs << v.second->point().x() << " "
        << v.second->point().y() << " "
        << v.second->point().z() << std::endl;

  for(const std::array<Vertex_handle, 3>& f : facets)
    ofs << "3  " << vertices.left.at(f[0]) << " "
                 << vertices.left.at(f[1]) << " "
                 << vertices.left.at(f[2]) << std::endl;

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

template<typename CellRange>
void dump_cells_polylines(const CellRange& cells, const char* filename)
{
  std::ofstream ofs(filename);
  ofs.precision(17);
  for (auto c : cells)
  {
    ofs << "2 " << point(c->vertex(0)->point()) << " "
                << point(c->vertex(1)->point()) <<std::endl;
    ofs << "2 " << point(c->vertex(0)->point()) << " "
                << point(c->vertex(2)->point()) << std::endl;
    ofs << "2 " << point(c->vertex(0)->point()) << " "
                << point(c->vertex(3)->point()) << std::endl;
    ofs << "2 " << point(c->vertex(1)->point()) << " "
                << point(c->vertex(2)->point()) << std::endl;
    ofs << "2 " << point(c->vertex(1)->point()) << " "
                << point(c->vertex(3)->point()) << std::endl;
    ofs << "2 " << point(c->vertex(2)->point()) << " "
                << point(c->vertex(3)->point()) << std::endl;
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
    ofs << point(vit->second->point()) << std::endl;
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
        << point(e.first->vertex(e.second)->point()) << " "
        << point(e.first->vertex(e.third)->point()) << "\n";
  }
  ofs.close();
}

template<typename Tr, typename CellSelector>
void dump_cells_with_small_dihedral_angle(const Tr& tr,
    const double angle_bound,
    CellSelector cell_select,
    const char* filename)
{
  typedef typename Tr::Cell_handle           Cell_handle;
  typedef typename Tr::Cell::Subdomain_index Subdomain_index;
  std::vector<Cell_handle>     cells;
  std::vector<Subdomain_index> indices;

  for (typename Tr::Finite_cells_iterator cit = tr.finite_cells_begin();
       cit != tr.finite_cells_end(); ++cit)
  {
    Cell_handle c = cit;
    if (c->subdomain_index() != Subdomain_index() && cell_select(c))
    {
      double dh = min_dihedral_angle(tr, c);
      if (dh < angle_bound)
      {
        cells.push_back(c);
        indices.push_back(c->subdomain_index());
      }
    }
  }
  std::cout << "bad cells : " << cells.size() << std::endl;
  dump_cells<Tr>(cells, indices, filename);
  dump_cells_off(cells, tr, "bad_cells.off");
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
    if (vit->in_dimension() == -1)
      continue;//far point
    CGAL_assertion(vit->in_dimension() >= 0 && vit->in_dimension() < 4);

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

    for (Vertex_handle vj : vertices_di)
    {
      ofs << point(vj->point()) << std::endl;
    }

    ofs.close();
  }
}

template<typename Tr>
void dump_triangulation_cells(const Tr& tr, const char* filename)
{
  std::vector<typename Tr::Cell_handle> cells(tr.number_of_finite_cells());
  std::vector<typename Tr::Cell::Subdomain_index> indices(tr.number_of_finite_cells());
  int i = 0;
  for (typename Tr::Finite_cells_iterator cit = tr.finite_cells_begin();
       cit != tr.finite_cells_end(); ++cit)
  {
    cells[i] = cit;
    indices[i++] = cit->subdomain_index();
  }
  dump_cells<Tr>(cells, indices, filename);
}

template<typename C3t3>
void dump_binary(const C3t3& c3t3, const char* filename)
{
  std::ofstream os(filename, std::ios::binary | std::ios::out);
  CGAL::Mesh_3::save_binary_file(os, c3t3);
  os.close();
}

template<typename C3t3>
void dump_medit(const C3t3& c3t3, const char* filename)
{
  std::ofstream os(filename, std::ios::out);
  c3t3.output_to_medit(os, true, true);
  os.close();
}

template<typename C3t3>
void dump_c3t3(const C3t3& c3t3, const char* filename_no_extension)
{
  std::string filename_medit(filename_no_extension);
  filename_medit.append(".mesh");
  dump_medit(c3t3, filename_medit.c_str());

  std::string filename_binary(filename_no_extension);
  filename_binary.append(".binary.cgal");
  dump_binary(c3t3, filename_binary.c_str());
}


} //namespace debug
} //namespace Tetrahedral_remeshing
} //namespace CGAL

#endif //CGAL_INTERNAL_TET_REMESHING_HELPERS_H
