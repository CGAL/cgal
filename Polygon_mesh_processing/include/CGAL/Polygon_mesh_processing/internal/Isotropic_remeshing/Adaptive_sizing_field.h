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
  typedef typename Base::FT         FT;
  typedef typename Base::Point_3    Point_3;
  typedef typename Base::halfedge_descriptor halfedge_descriptor;
  typedef typename Base::vertex_descriptor   vertex_descriptor;
  typedef typename CGAL::dynamic_vertex_property_t<FT>            Vertex_property_tag;
  typedef typename boost::property_map<PolygonMesh,
                                       Vertex_property_tag>::type VertexSizingMap;

    Adaptive_sizing_field(const std::pair<FT, FT>& edge_len_min_max
                        , PolygonMesh& pmesh)
    : m_sq_short( CGAL::square(edge_len_min_max.first))
    , m_sq_long(  CGAL::square(edge_len_min_max.second))
    , m_pmesh(pmesh)
  {
      //todo ip: initialize sizing map with default values
      //todo ip: might end up using directly the property map of the curvature calculation (if mutable)?
      vertex_sizing_map_ = get(Vertex_property_tag(), m_pmesh);
      for(vertex_descriptor v : vertices(m_pmesh)){
          put(vertex_sizing_map_, v, m_sq_long);
      }
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
      //todo ip
      // calculate curvature

      // loop over curvature property field and calculate the target mesh size for a vertex
      // don't forget to store squared length

  }

  boost::optional<FT> is_too_long(const halfedge_descriptor& h) const
  {
    const FT sqlen = sqlength(h);
    FT sqtarg_len = std::min(get(vertex_sizing_map_, source(h, m_pmesh)),
                             get(vertex_sizing_map_, target(h, m_pmesh)));
    CGAL_assertion(get(vertex_sizing_map_, source(h, m_pmesh)));
    CGAL_assertion(get(vertex_sizing_map_, target(h, m_pmesh)));
    if(sqlen > sqtarg_len)
      return sqlen;
    else
      return boost::none;
  }

  boost::optional<FT> is_too_long(const vertex_descriptor& va,
                                  const vertex_descriptor& vb) const
  {
    const FT sqlen = sqlength(va, vb);
    FT sqtarg_len = std::min(get(vertex_sizing_map_, va),
                             get(vertex_sizing_map_, vb));
    CGAL_assertion(get(vertex_sizing_map_, va));
    CGAL_assertion(get(vertex_sizing_map_, vb));
    if (sqlen > sqtarg_len)
      return sqlen;
    else
      return boost::none;
  }

  boost::optional<FT> is_too_short(const halfedge_descriptor& h) const
  {
    const FT sqlen = sqlength(h);
    FT sqtarg_len = std::min(get(vertex_sizing_map_, source(h, m_pmesh)),
                             get(vertex_sizing_map_, target(h, m_pmesh)));
    CGAL_assertion(get(vertex_sizing_map_, source(h, m_pmesh)));
    CGAL_assertion(get(vertex_sizing_map_, target(h, m_pmesh)));
    if (sqlen < sqtarg_len)
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

  void update_sizing_map(const vertex_descriptor& vnew)
  {
    //todo ip: calculate curvature for the vertex
    //dummy
    put(vertex_sizing_map_, vnew, m_sq_short);
  }

private:
  FT m_sq_short;
  FT m_sq_long;
  PolygonMesh& m_pmesh;
  VertexSizingMap vertex_sizing_map_;
};

}//end namespace Polygon_mesh_processing
}//end namespace CGAL

#endif //CGAL_PMP_REMESHING_ADAPTIVE_SIZING_FIELD_H
