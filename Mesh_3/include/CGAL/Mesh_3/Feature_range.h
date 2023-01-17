// Copyright (c) 2023 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Jane Tournois
//
//******************************************************************************
//
//******************************************************************************

#ifndef CGAL_MESH_3_FEATURE_RANGE_H
#define CGAL_MESH_3_FEATURE_RANGE_H

#include <CGAL/license/Mesh_3.h>

#include <CGAL/Mesh_3/polylines_to_protect.h>

#include <vector>

namespace CGAL
{
namespace Mesh_3
{
/*!
* \ingroup PkgMesh3FeatureDetection
*
* Functor for feature addition
*/
template<typename PolylineRange>
struct Feature_range
{
private:
  const PolylineRange& polylines_;

public:
  Feature_range(const PolylineRange& polylines)
    : polylines_(polylines)
  {}

public:
  /*!
  *
  * \tparam Mesh_domain class model of `MeshDomainWithFeatures_3`
  * \param domain the mesh domain to be enriched with polyline features
  */
  template<typename Mesh_domain>
  std::vector<std::vector<typename Mesh_domain::Point_3>> operator()(Mesh_domain& domain)
  {
    using Point = typename Mesh_domain::Point_3;

    std::vector<std::vector<Point>> polyline_graph;
    CGAL::polylines_to_protect(polyline_graph, std::begin(polylines_), std::end(polylines_));

    return polyline_graph;
  }

  const PolylineRange& polylines() { return polylines_; }
};

template<typename PolylineRange>
auto feature_range(const PolylineRange range)
{
  return Feature_range<PolylineRange>(range);
}

}//end namespace Mesh_3
}//end namespace CGAL


#endif //CGAL_MESH_3_FEATURE_RANGE_H
