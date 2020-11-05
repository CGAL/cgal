// Copyright (c) 2018 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Konstantinos Katrioplas (konst.katrioplas@gmail.com)

#ifndef CGAL_POLYGON_MESH_PROCESSING_INTERNAL_SMOOTHING_EVALUATION_H
#define CGAL_POLYGON_MESH_PROCESSING_INTERNAL_SMOOTHING_EVALUATION_H

#include <CGAL/license/Polygon_mesh_processing/meshing_hole_filling.h>

#include <CGAL/Polygon_mesh_processing/measure.h>

#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>

#include <set>
#include <vector>

namespace CGAL {
namespace Polygon_mesh_processing {
namespace internal {

template<typename PolygonMesh, typename GeomTraits>
class Quality_evaluator
{
  typedef typename GeomTraits::Point_3                                          Point;
  typedef typename GeomTraits::Vector_3                                         Vector;
  typedef typename GeomTraits::Line_3                                           Line;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor        halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor            face_descriptor;
  typedef typename boost::property_map<PolygonMesh, CGAL::vertex_point_t>::type VertexPointMap;
  typedef typename boost::property_traits<VertexPointMap>::reference            Point_ref;

public:
  Quality_evaluator(PolygonMesh& pmesh,
                    const GeomTraits& traits)
    : mesh_(pmesh), traits_(traits)
  {
    std::size_t number_of_triangles = faces(mesh_).size();
    angles_.reserve(number_of_triangles * 3);
    areas_.reserve(number_of_triangles);
    aspect_ratios_.reserve(number_of_triangles);
    vpmap_ = get(CGAL::vertex_point, mesh_);
  }

  void gather_angles()
  {
    for(halfedge_descriptor hi : halfedges(mesh_))
    {
      const Point_ref a = get(vpmap_, source(hi, mesh_));
      const Point_ref b = get(vpmap_, target(hi, mesh_));
      const Point_ref c = get(vpmap_, target(next(hi, mesh_), mesh_));

      angles_.push_back(traits_.compute_approximate_angle_3_object()(a, b, c));
    }

#ifdef CGAL_PMP_SMOOTHING_DEBUG
    std::cout << "angles_ size = " << angles_.size() << std::endl;
#endif
  }

  template <typename Stream>
  void extract_angles(Stream& output)
  {
    for(unsigned int i=0; i!=angles_.size(); ++i)
      output << angles_[i] << std::endl;
    output.close();
  }

  void measure_areas()
  {
    for(face_descriptor f : faces(mesh_))
      areas_.push_back(face_area(f, mesh_));

#ifdef CGAL_PMP_SMOOTHING_DEBUG
    std::cout << "areas_ size = " << areas_.size() << std::endl;
#endif
  }

  template <typename Stream>
  void extract_areas(Stream& output)
  {
    for(unsigned int i=0; i!=areas_.size(); ++i)
      output << areas_[i] << std::endl;
    output.close();
  }

  void calc_aspect_ratios()
  {
    for(face_descriptor f : faces(mesh_))
      aspect_ratios_.push_back(CGAL::Polygon_mesh_processing::face_aspect_ratio(f, mesh_));

#ifdef CGAL_PMP_SMOOTHING_DEBUG
    std::cout << "aspect_ratios_ size = " << aspect_ratios_.size() << std::endl;
#endif
  }

  template <typename Stream>
  void extract_aspect_ratios(Stream& output)
  {
    for(unsigned int i=0; i!=aspect_ratios_.size(); ++i)
      output << aspect_ratios_[i] << std::endl;
    output.close();
  }

private:
  PolygonMesh& mesh_;
  const GeomTraits& traits_;
  VertexPointMap vpmap_;

  std::vector<double> angles_;
  std::vector<double> areas_;
  std::vector<double> aspect_ratios_;
};

} // namespace internal
} // namespace Polygon_mesh_processing
} // namespace CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_INTERNAL_SMOOTHING_EVALUATION_H
