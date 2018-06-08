// Copyright (c) 2018 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Konstantinos Katrioplas (konst.katrioplas@gmail.com)

#ifndef CGAL_POLYGON_MESH_PROCESSING_INTERNAL_SMOOTHING_EVALUATION_H
#define CGAL_POLYGON_MESH_PROCESSING_INTERNAL_SMOOTHING_EVALUATION_H

#include <CGAL/license/Polygon_mesh_processing/meshing_hole_filling.h>

#include <vector>
#include <set>
#include <boost/graph/graph_traits.hpp>
#include <boost/foreach.hpp>
#include <boost/property_map/property_map.hpp>
#include <CGAL/Polygon_mesh_processing/measure.h>

namespace CGAL {
namespace Polygon_mesh_processing {
namespace internal {

template<typename PolygonMesh, typename GeomTraits>
class Quality_evaluator
{
  typedef typename GeomTraits::Point_3 Point;
  typedef typename GeomTraits::Vector_3 Vector;
  typedef typename GeomTraits::Line_3 Line;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor face_descriptor;
  typedef typename boost::property_map<PolygonMesh, CGAL::vertex_point_t>::type VertexPointMap;

public:
  Quality_evaluator(PolygonMesh& pmesh) : mesh_(pmesh)
  {
    std::size_t number_of_triangles = faces(mesh_).size();
    angles_.reserve(number_of_triangles * 3);
    areas_.reserve(number_of_triangles);
    aspect_ratios_.reserve(number_of_triangles);
    vpmap_ = get(CGAL::vertex_point, mesh_);
  }

  void gather_angles()
  {
    double rad_to_deg = 180. / CGAL_PI;

    BOOST_FOREACH(halfedge_descriptor hi, halfedges(mesh_))
    {
      Point a = get(vpmap_, source(hi, mesh_));
      Point b = get(vpmap_, target(hi, mesh_));
      Point c = get(vpmap_, target(next(hi, mesh_), mesh_));
      Vector ba(b, a);
      Vector bc(b, c);
      double cos_angle = CGAL::to_double((ba * bc))
        / CGAL::approximate_sqrt(ba.squared_length() * bc.squared_length());
      angles_.push_back(std::acos(cos_angle) * rad_to_deg);
    }
    #ifdef CGAL_PMP_SMOOTHING_VERBOSE
    std::cout<<"angles_ size= "<<angles_.size()<<std::endl;
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
    BOOST_FOREACH(face_descriptor f, faces(mesh_))
      areas_.push_back(face_area(f, mesh_));

    #ifdef CGAL_PMP_SMOOTHING_VERBOSE
    std::cout<<"areas_ size= "<<areas_.size()<<std::endl;
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
    BOOST_FOREACH(face_descriptor f, faces(mesh_))
    {
      halfedge_descriptor h = halfedge(f, mesh_);
      Point points[3];
      points[0] = get(vpmap_, target(h, mesh_));
      points[1] = get(vpmap_, target(next(h, mesh_), mesh_));
      points[2] = get(vpmap_, target(next(next(h, mesh_), mesh_), mesh_));

      double min_alt = std::numeric_limits<double>::infinity();
      double longest_edge = 0;

      for(int i=0; i<3; ++i)
      {
        double alt = CGAL::approximate_sqrt(CGAL::squared_distance(points[(0+i)%3], Line(points[(1+i)%3], points[(2+i)%3])));
        double edge =  CGAL::approximate_sqrt(CGAL::squared_distance(points[(1+i)%3], points[(2+i)%3]));
        if(alt < min_alt) { min_alt = alt; }
        if(edge > longest_edge) { longest_edge = edge; }
      }
      aspect_ratios_.push_back(longest_edge / min_alt);
    }
    #ifdef CGAL_PMP_SMOOTHING_VERBOSE
    std::cout<<"aspect_ratios_ size= "<<aspect_ratios_.size()<<std::endl;
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
  VertexPointMap vpmap_;
  std::vector<double> angles_;
  std::vector<double> areas_;
  std::vector<double> aspect_ratios_;

};

} //namespaces
}
}


#endif // CGAL_POLYGON_MESH_PROCESSING_INTERNAL_SMOOTHING_EVALUATION_H
