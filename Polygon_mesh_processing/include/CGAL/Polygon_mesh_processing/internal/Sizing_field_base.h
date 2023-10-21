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

#ifndef CGAL_PMP_REMESHING_SIZING_FIELD_H
#define CGAL_PMP_REMESHING_SIZING_FIELD_H

#include <CGAL/license/Polygon_mesh_processing/meshing_hole_filling.h>

#include <CGAL/Kernel_traits.h>
#include <CGAL/boost/graph/properties.h>

#include <boost/optional.hpp>

namespace CGAL
{
namespace Polygon_mesh_processing
{
namespace internal
{
/*!
* \ingroup PMP_meshing_grp
* pure virtual class serving as a base for sizing field classes used in isotropic
* remeshing.
*
* \cgalModels{PMPSizingField}
*
* \sa `isotropic_remeshing()`
* \sa `Uniform_sizing_field`
* \sa `Adaptive_sizing_field`
*
* @tparam PolygonMesh model of `MutableFaceGraph` that
*         has an internal property map for `CGAL::vertex_point_t`.
* @tparam VPMap property map associating points to the vertices of `pmesh`,
*         model of `ReadWritePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
*         as key type and `%Point_3` as value type. Default is `boost::get(CGAL::vertex_point, pmesh)`.
*/
template <class PolygonMesh,
          class VPMap = typename boost::property_map<PolygonMesh, CGAL::vertex_point_t>::const_type>
class Sizing_field_base
{
private:
  typedef PolygonMesh                     PM;
  typedef typename boost::property_traits<VPMap>::value_type Point;

public:
  typedef typename CGAL::Kernel_traits<Point>::Kernel K;
  typedef typename boost::graph_traits<PM>::face_descriptor     face_descriptor;
  typedef typename boost::graph_traits<PM>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<PM>::vertex_descriptor   vertex_descriptor;
  typedef Point                                                 Point_3;
  typedef typename K::FT                                        FT;

public:
  virtual FT at(const vertex_descriptor v) const = 0;
  virtual std::optional<FT> is_too_long(const vertex_descriptor va,
                                        const vertex_descriptor vb) const = 0;
  virtual std::optional<FT> is_too_short(const halfedge_descriptor h,
                                         const PolygonMesh& pmesh) const = 0;
  virtual Point_3 split_placement(const halfedge_descriptor h, const PolygonMesh& pmesh) const = 0;
  virtual void update(const vertex_descriptor v, const PolygonMesh& pmesh) = 0;

};

}//end namespace internal
}//end namespace Polygon_mesh_processing
}//end namespace CGAL

#endif //CGAL_PMP_REMESHING_SIZING_FIELD_H
