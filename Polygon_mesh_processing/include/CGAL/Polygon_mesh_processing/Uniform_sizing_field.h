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

#ifndef CGAL_PMP_REMESHING_UNIFORM_SIZING_FIELD_H
#define CGAL_PMP_REMESHING_UNIFORM_SIZING_FIELD_H

#include <CGAL/license/Polygon_mesh_processing/meshing_hole_filling.h>

#include <CGAL/Polygon_mesh_processing/internal/Isotropic_remeshing/Sizing_field.h>

#include <CGAL/number_utils.h>

namespace CGAL
{
namespace Polygon_mesh_processing
{
template <class PolygonMesh,
          class VPMap =  typename boost::property_map<PolygonMesh, CGAL::vertex_point_t>::const_type>
class Uniform_sizing_field : public CGAL::Sizing_field<PolygonMesh, VPMap>
{
private:
  typedef CGAL::Sizing_field<PolygonMesh, VPMap> Base;

public:
  typedef typename Base::FT         FT;
  typedef typename Base::Point_3    Point_3;
  typedef typename Base::halfedge_descriptor halfedge_descriptor;
  typedef typename Base::vertex_descriptor   vertex_descriptor;

  Uniform_sizing_field<PolygonMesh, VPMap>(const FT size, const VPMap& vpmap)
    : m_sq_short( CGAL::square(4./5. * size))
    , m_sq_long(  CGAL::square(4./3. * size))
    , m_vpmap(vpmap)
  {}

  Uniform_sizing_field<PolygonMesh>(const FT& size, const PolygonMesh& pmesh)
    : m_sq_short( CGAL::square(4./5. * size))
    , m_sq_long(  CGAL::square(4./3. * size))
    , m_vpmap(get(CGAL::vertex_point, pmesh))
  {}

private:
  FT sqlength(const vertex_descriptor va,
              const vertex_descriptor vb) const
  {
    return FT(CGAL::squared_distance(get(m_vpmap, va), get(m_vpmap, vb)));
  }

  FT sqlength(const halfedge_descriptor& h, const PolygonMesh& pmesh) const
  {
    return sqlength(target(h, pmesh), source(h, pmesh));
  }

public:
  boost::optional<FT> is_too_long(const halfedge_descriptor h, const PolygonMesh& pmesh) const
  {
    const FT sqlen = sqlength(h, pmesh);
    if(sqlen > m_sq_long)
      return sqlen;
    else
      return boost::none;
  }

  boost::optional<FT> is_too_long(const vertex_descriptor va,
                                  const vertex_descriptor vb) const
  {
    const FT sqlen = sqlength(va, vb);
    if (sqlen > m_sq_long)
      return sqlen;
    else
      return boost::none;
  }

  boost::optional<FT> is_too_short(const halfedge_descriptor h, const PolygonMesh& pmesh) const
  {
    const FT sqlen = sqlength(h, pmesh);
    if (sqlen < m_sq_long)
      return sqlen;
    else
      return boost::none;
  }

  virtual Point_3 split_placement(const halfedge_descriptor h, const PolygonMesh& pmesh) const
  {
    return CGAL::midpoint(get(m_vpmap, target(h, pmesh)),
                          get(m_vpmap, source(h, pmesh)));
  }

private:
  FT m_sq_short;
  FT m_sq_long;
  const VPMap m_vpmap;
};

}//end namespace Polygon_mesh_processing
}//end namespace CGAL

#endif //CGAL_PMP_REMESHING_UNIFORM_SIZING_FIELD_H
