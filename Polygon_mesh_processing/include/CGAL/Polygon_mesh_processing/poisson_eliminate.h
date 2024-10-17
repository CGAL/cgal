// Copyright (c) 2024 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Andreas Fabri

#ifndef CGAL_POLYGON_MESH_PROCESSING_POISSON_ELIMINATE_H
#define CGAL_POLYGON_MESH_PROCESSING_POISSON_ELIMINATE_H

#include <CGAL/license/Polygon_mesh_processing/distance.h>


#include <cyVector.h>
#include <cySampleElim.h>

#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/Polygon_mesh_processing/measure.h>


namespace CGAL {
namespace Polygon_mesh_processing {

template <class TriangleMesh, class OutputIterator, class NamedParameters = parameters::Default_named_parameters>
poisson_eliminate(const TriangleMesh& sm, OutputIterator out, const NamedParameters& np = parameters::default_values())
{
  typedef typename GetGeomTraits<TriangleMesh, NamedParameters>::type             GeomTraits;
  typedef typename GeomTraits::Point_3                                            Point_3;

  Bbox_3 bb = bbox_3(sm.points().begin(), sm.points().end());
  cy::Vec3d bl(bb.xmin(), bb.ymin(), bb.zmin());
  cy::Vec3d tr(bb.xmax(), bb.ymax(), bb.zmax());

  std::vector<cy::Vec3d> inputPoints, outputPoints;
  std::vector<Point_3> points;

  // @todo write with a transform_iterator directly into inputPoints
  sample_triangle_mesh(sm,
                       std::back_inserter(points),
                       CGAL::parameters::number_of_points_on_faces(2* num_vertices(sm))
                         .do_sample_vertices(false)
                         .do_sample_edges(false));
  double area = CGAL::Polygon_mesh_processing::area(sm);

  for(int i = 0; i < points.size(); ++i){
    inputPoints.push_back(cy::Vec3d(to_double(points[i].x()), to_double(points[i].y()), to_double(points[i].z())));
  }

  outputPoints.resize(num_vertices(sm)/2);

  cy::WeightedSampleElimination< cy::Vec3d, double, 3, int > wse;
  wse.SetBoundsMin(bl);
  wse.SetBoundsMax(tr);
  bool isProgressive = true;

  double d_max = 2 * wse.GetMaxPoissonDiskRadius( 2, outputPoints.size(), area );

  wse.Eliminate( inputPoints.data(), inputPoints.size(),
                 outputPoints.data(), outputPoints.size(),
                 isProgressive,
                 d_max, 2 );

  for (const cy::Vec3d& p : outputPoints){
    *out++ = Point_3(p.x, p.y, p.z);
  }
}

} // namespace Polygon_mesh_processing
} // namespace CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_POISSON_ELIMINATE_H