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

#include <CGAL/Polygon_mesh_processing/Sizing_field_base.h>

#include <CGAL/number_utils.h>

namespace CGAL
{
namespace Polygon_mesh_processing
{
/*!
* \ingroup PMP_meshing_grp
* provides a set of instructions for isotropic remeshing to achieve uniform
* mesh edge lengths.
*
* Edges longer than the 4/3 of the target edge length will be split in half, while
* edges shorter than the 4/5 of the target edge length will be collapsed.
*
* \cgalModels{PMPSizingField}
*
* \sa `isotropic_remeshing()`
* \sa `Adaptive_sizing_field`
*
* @tparam PolygonMesh model of `MutableFaceGraph` that
*         has an internal property map for `CGAL::vertex_point_t`.
* @tparam VPMap a property map associating points to the vertices of `pmesh`.
*         It is a a class model of `ReadWritePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
*          as key type and `%Point_3` as value type. Default is `boost::get(CGAL::vertex_point, pmesh)`.
*/
template <class PolygonMesh,
          class VPMap =  typename boost::property_map<PolygonMesh, CGAL::vertex_point_t>::const_type>
class Uniform_sizing_field
#ifndef DOXYGEN_RUNNING
  : public Sizing_field_base<PolygonMesh, VPMap>
#endif
{
private:
  typedef Sizing_field_base<PolygonMesh, VPMap> Base;

public:
  typedef typename Base::FT         FT;
  typedef typename Base::Point_3    Point_3;
  typedef typename Base::halfedge_descriptor halfedge_descriptor;
  typedef typename Base::vertex_descriptor   vertex_descriptor;

  /// \name Creation
  /// @{

  /*!
  * returns an object to serve as criterion for uniform edge lengths.
  * \param size is the target edge length for the isotropic remeshing. If set to 0,
  *        the criterion for edge length is ignored and edges are neither split nor collapsed.
  * \param vpmap is the input vertex point map that associates points to the vertices of
  *        an input mesh.
  */
  Uniform_sizing_field<PolygonMesh, VPMap>(const FT size, const VPMap& vpmap)
    : m_sq_short( CGAL::square(4./5. * size))
    , m_sq_long(  CGAL::square(4./3. * size))
    , m_vpmap(vpmap)
  {}

  /*!
  * returns an object to serve as criterion for uniform edge lengths. It calls the first
  * constructor using default values for the vertex point map of the input polygon mesh.
  *
  * @param size is the target edge length for the isotropic remeshing. If set to 0,
  *        the criterion for edge length is ignored and edges are neither split nor collapsed.
  * @param pmesh a polygon mesh with triangulated surface patches to be remeshed. The default
  *        vertex point map of pmesh is used to construct the class.
  */
  Uniform_sizing_field<PolygonMesh, VPMap>(const FT size, const PolygonMesh& pmesh)
    : Uniform_sizing_field(size, get(CGAL::vertex_point, pmesh))
  {}

  /// @}

private:
  FT sqlength(const vertex_descriptor va,
              const vertex_descriptor vb) const
  {
    return FT(squared_distance(get(m_vpmap, va), get(m_vpmap, vb)));
  }

  FT sqlength(const halfedge_descriptor& h, const PolygonMesh& pmesh) const
  {
    return sqlength(target(h, pmesh), source(h, pmesh));
  }

public:
  std::optional<FT> is_too_long(const halfedge_descriptor h, const PolygonMesh& pmesh) const
  {
    const FT sqlen = sqlength(h, pmesh);
    if(sqlen > m_sq_long)
      return sqlen;
    else
      return std::nullopt;
  }

  std::optional<FT> is_too_long(const vertex_descriptor va, const vertex_descriptor vb) const
  {
    const FT sqlen = sqlength(va, vb);
    if (sqlen > m_sq_long)
      return sqlen;
    else
      return std::nullopt;
  }

  std::optional<FT> is_too_short(const halfedge_descriptor h, const PolygonMesh& pmesh) const
  {
    const FT sqlen = sqlength(h, pmesh);
    if (sqlen < m_sq_short)
      return sqlen;
    else
      return std::nullopt;
  }

  Point_3 split_placement(const halfedge_descriptor h, const PolygonMesh& pmesh) const
  {
    return midpoint(get(m_vpmap, target(h, pmesh)),
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
