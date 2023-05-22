// Copyright (c) 2020 GeometryFactory (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Jane Tournois

#ifndef CGAL_PMP_REMESHING_ADAPTIVE_SIZING_FIELD_H
#define CGAL_PMP_REMESHING_ADAPTIVE_SIZING_FIELD_H

#include <CGAL/license/Polygon_mesh_processing/meshing_hole_filling.h>

#include <CGAL/Polygon_mesh_processing/internal/Isotropic_remeshing/Sizing_field.h>

#include <CGAL/Polygon_mesh_processing/interpolated_corrected_curvatures.h>

#include <CGAL/number_utils.h>

namespace CGAL
{
namespace Polygon_mesh_processing
{
template <class PolygonMesh>
class Adaptive_sizing_field : public CGAL::Sizing_field<PolygonMesh>
{
private:
  typedef CGAL::Sizing_field<PolygonMesh> Base;

public:
  typedef typename Base::K          K;
  typedef typename Base::FT         FT;
  typedef typename Base::Point_3    Point_3;
  typedef typename Base::halfedge_descriptor halfedge_descriptor;
  typedef typename Base::vertex_descriptor   vertex_descriptor;

  typedef typename CGAL::dynamic_vertex_property_t<FT> Vertex_property_tag;
  typedef typename boost::property_map<PolygonMesh,
                                       Vertex_property_tag>::type VertexSizingMap;

  //todo ip: set a property map that can calculate curvature in one go. I think I'm generating constant maps (without put)
  // try 1
  typedef Principal_curvatures_and_directions<K> Principal_curvatures;
//  typedef Constant_property_map<vertex_descriptor, Principal_curvatures> Vertex_curvature_map;

  // try 2
  typedef Constant_property_map<vertex_descriptor, Principal_curvatures_and_directions<K>> Default_principal_map;
  typedef typename internal_np::Lookup_named_param_def<internal_np::vertex_principal_curvatures_and_directions_map_t,
                                                       parameters::Default_named_parameters,
                                                       Default_principal_map>::type
                                                         Vertex_curvature_map;

    Adaptive_sizing_field(const double tol
                        , const std::pair<FT, FT>& edge_len_min_max
                        , PolygonMesh& pmesh)
    : tol(tol)
    , m_sq_short(CGAL::square(edge_len_min_max.first))
    , m_sq_long( CGAL::square(edge_len_min_max.second))
    , m_pmesh(pmesh)
  {
    m_vertex_sizing_map = get(Vertex_property_tag(), m_pmesh);
  }

private:
  FT sqlength(const vertex_descriptor& va,
              const vertex_descriptor& vb) const
  {
    typename boost::property_map<PolygonMesh, CGAL::vertex_point_t>::const_type
      vpmap = get(CGAL::vertex_point, m_pmesh);
    return FT(CGAL::squared_distance(get(vpmap, va), get(vpmap, vb)));
  }

  FT sqlength(const halfedge_descriptor& h) const
  {
    return sqlength(target(h, m_pmesh), source(h, m_pmesh));
  }

public:
  void calc_sizing_map()
  {
#ifdef CGAL_PMP_REMESHING_VERBOSE
    int oversize  = 0;
    int undersize = 0;
    int insize    = 0;
    std::cout << "Calculating sizing field..." << std::endl;
#endif

    //todo ip: how to make this work?
//    Vertex_curvature_map vertex_curvature_map;
//    interpolated_corrected_principal_curvatures_and_directions(m_pmesh
//                                                               , vertex_curvature_map);

    // calculate square vertex sizing field (L(x_i))^2 from curvature field
    for(vertex_descriptor v : vertices(m_pmesh))
    {
//      auto vertex_curv = get(vertex_curvature_map, v); //todo ip: how to make this work?
      //todo ip: temp solution
      const Principal_curvatures vertex_curv = interpolated_corrected_principal_curvatures_and_directions_one_vertex(m_pmesh, v);
      const FT max_absolute_curv = std::max(std::abs(vertex_curv.max_curvature), std::abs(vertex_curv.min_curvature));
      const FT vertex_size_sq = 6 * tol / max_absolute_curv - 3 * CGAL::square(tol);
      if (vertex_size_sq > m_sq_long)
      {
        put(m_vertex_sizing_map, v, m_sq_long);
#ifdef CGAL_PMP_REMESHING_VERBOSE
        ++oversize;
#endif
      }
      else if (vertex_size_sq < m_sq_short)
      {
        put(m_vertex_sizing_map, v, m_sq_short);
#ifdef CGAL_PMP_REMESHING_VERBOSE
        ++undersize;
#endif
      }
      else
      {
        put(m_vertex_sizing_map, v, vertex_size_sq);
#ifdef CGAL_PMP_REMESHING_VERBOSE
        ++insize;
#endif
      }
    }
#ifdef CGAL_PMP_REMESHING_VERBOSE
    std::cout << " done (" << insize    << " from curvature, "
                          << oversize  << " set to max, "
                          << undersize << " set to min)" << std::endl;
#endif
  }

  boost::optional<FT> is_too_long(const halfedge_descriptor& h) const
  {
    const FT sqlen = sqlength(h);
    FT sqtarg_len = std::min(get(m_vertex_sizing_map, source(h, m_pmesh)),
                             get(m_vertex_sizing_map, target(h, m_pmesh)));
    CGAL_assertion(get(m_vertex_sizing_map, source(h, m_pmesh)));
    CGAL_assertion(get(m_vertex_sizing_map, target(h, m_pmesh)));
    if(sqlen > sqtarg_len)
      return sqlen;
    else
      return boost::none;
  }

  boost::optional<FT> is_too_long(const vertex_descriptor& va,
                                  const vertex_descriptor& vb) const
  {
    const FT sqlen = sqlength(va, vb);
    FT sqtarg_len = std::min(get(m_vertex_sizing_map, va),
                             get(m_vertex_sizing_map, vb));
    CGAL_assertion(get(m_vertex_sizing_map, va));
    CGAL_assertion(get(m_vertex_sizing_map, vb));
    if (sqlen > 16./9. * sqtarg_len)
      return sqlen;
    else
      return boost::none;
  }

  boost::optional<FT> is_too_short(const halfedge_descriptor& h) const
  {
    const FT sqlen = sqlength(h);
    FT sqtarg_len = std::min(get(m_vertex_sizing_map, source(h, m_pmesh)),
                             get(m_vertex_sizing_map, target(h, m_pmesh)));
    CGAL_assertion(get(m_vertex_sizing_map, source(h, m_pmesh)));
    CGAL_assertion(get(m_vertex_sizing_map, target(h, m_pmesh)));
    if (sqlen < 16./25. * sqtarg_len)
      return sqlen;
    else
      return boost::none;
  }

  virtual Point_3 split_placement(const halfedge_descriptor& h) const
  {
    typename boost::property_map<PolygonMesh, CGAL::vertex_point_t>::const_type
      vpmap = get(CGAL::vertex_point, m_pmesh);
    return CGAL::midpoint(get(vpmap, target(h, m_pmesh)),
                          get(vpmap, source(h, m_pmesh)));
  }

  void update_sizing_map(const vertex_descriptor& v)
  {
    // calculating it as the average of two vertices on other ends
    // of halfedges as updating is done during an edge split
    int i = 0;
    FT vertex_size_sq = 0;
    CGAL_assertion(CGAL::halfedges_around_target(v, m_pmesh) == 2);
    for (halfedge_descriptor ha: CGAL::halfedges_around_target(v, m_pmesh))
    {
      vertex_size_sq += get(m_vertex_sizing_map, source(ha, m_pmesh));
      ++i;
    }
    vertex_size_sq /= i;

    put(m_vertex_sizing_map, v, vertex_size_sq);
  }

  //todo ip: is_protected_constraint_too_long() from PR

private:
  const FT tol;
  const FT m_sq_short;
  const FT m_sq_long;
  PolygonMesh& m_pmesh;
  VertexSizingMap m_vertex_sizing_map;
};

}//end namespace Polygon_mesh_processing
}//end namespace CGAL

#endif //CGAL_PMP_REMESHING_ADAPTIVE_SIZING_FIELD_H
