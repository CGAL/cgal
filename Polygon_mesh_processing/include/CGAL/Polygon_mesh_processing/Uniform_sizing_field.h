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

#include <CGAL/Polygon_mesh_processing/internal/Sizing_field_base.h>

#include <CGAL/number_utils.h>

namespace CGAL
{
namespace Polygon_mesh_processing
{
/*!
* \ingroup PMP_meshing_grp
* a sizing field describing a uniform target edge length for
* `CGAL::Polygon_mesh_processing::isotropic_remeshing()`.
*
* Edges longer than 4/3 of the target edge length will be split in half, while
* edges shorter than 4/5 of the target edge length will be collapsed.
*
* \cgalModels{PMPSizingField}
*
* \sa `isotropic_remeshing()`
* \sa `Adaptive_sizing_field`
*
* @tparam PolygonMesh model of `MutableFaceGraph` that
*         has an internal property map for `CGAL::vertex_point_t`.
* @tparam VPMap property map associating points to the vertices of `pmesh`,
*         model of `ReadWritePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
*         as key type and `%Point_3` as value type. Default is `boost::get(CGAL::vertex_point, pmesh)`.
*/
template <class PolygonMesh,
          class VPMap =  typename boost::property_map<PolygonMesh, CGAL::vertex_point_t>::const_type>
class Uniform_sizing_field
#ifndef DOXYGEN_RUNNING
: public internal::Sizing_field_base<PolygonMesh, VPMap>
#endif
{
private:
  typedef internal::Sizing_field_base<PolygonMesh, VPMap> Base;

public:
  typedef typename Base::FT         FT;
  typedef typename Base::Point_3    Point_3;
  typedef typename Base::halfedge_descriptor halfedge_descriptor;
  typedef typename Base::vertex_descriptor   vertex_descriptor;

  /// \name Creation
  /// @{

  /*!
  * Constructor.
  * \param size the target edge length for isotropic remeshing. If set to 0,
  *        the criterion for edge length is ignored and edges are neither split nor collapsed.
  * \param vpmap is the input vertex point map that associates points to the vertices of
  *        the input mesh.
  */
  Uniform_sizing_field(const FT size, const VPMap& vpmap)
    : m_size(size)
    , m_sq_short( CGAL::square(4./5. * size))
    , m_sq_long(  CGAL::square(4./3. * size))
    , m_vpmap(vpmap)
  {}

  /*!
  * Constructor using internal vertex point map of the input polygon mesh.
  *
  * @param size the target edge length for isotropic remeshing. If set to 0,
  *        the criterion for edge length is ignored and edges are neither split nor collapsed.
  * @param pmesh a polygon mesh with triangulated surface patches to be remeshed. The default
  *        vertex point map of `pmesh` is used to construct the class.
  */
  Uniform_sizing_field(const FT size, const PolygonMesh& pmesh)
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
  FT at(const vertex_descriptor /* v */) const
  {
    return m_size;
  }

  std::optional<FT> is_too_long(const vertex_descriptor va, const vertex_descriptor vb) const
  {
    const FT sqlen = sqlength(va, vb);
    if (sqlen > m_sq_long)
      //no need to return the ratio for the uniform field
      return sqlen;
    else
      return std::nullopt;
  }

  std::optional<FT> is_too_short(const halfedge_descriptor h, const PolygonMesh& pmesh) const
  {
    const FT sqlen = sqlength(h, pmesh);
    if (sqlen < m_sq_short)
      //no need to return the ratio for the uniform field
      return sqlen;
    else
      return std::nullopt;
  }

  Point_3 split_placement(const halfedge_descriptor h, const PolygonMesh& pmesh) const
  {
    return midpoint(get(m_vpmap, target(h, pmesh)),
                    get(m_vpmap, source(h, pmesh)));
  }

  void update(const vertex_descriptor /* v */, const PolygonMesh& /* pmesh */)
  {}

private:
  const FT m_size;
  const FT m_sq_short;
  const FT m_sq_long;
  const VPMap m_vpmap;
};

}//end namespace Polygon_mesh_processing
}//end namespace CGAL

#endif //CGAL_PMP_REMESHING_UNIFORM_SIZING_FIELD_H
