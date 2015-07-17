// Copyright (c) 2015 GeometryFactory (France).
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
//
//
// Author(s)     : Sebastien Loriot

#ifndef CGAL_POLYGON_MESH_PROCESSING_DISTANCE_H
#define CGAL_POLYGON_MESH_PROCESSING_DISTANCE_H

#include <algorithm>
#include <cmath>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>

/// \cond SKIP_IN_MANUAL

namespace CGAL{
namespace Polygon_mesh_processing {
namespace internal{

template <class Kernel, class OutputIterator>
OutputIterator
triangle_grid_sampling(const typename Kernel::Point_3& p0,
                const typename Kernel::Point_3& p1,
                const typename Kernel::Point_3& p2,
                double distance,
                OutputIterator out)
{
  const double d_p0p1 = std::sqrt( CGAL::squared_distance(p0, p1) );
  const double d_p0p2 = std::sqrt( CGAL::squared_distance(p0, p2) );

  const double n = (std::max)(std::ceil( d_p0p1 / distance ),
                                   std::ceil( d_p0p2 / distance ));

  for (double i=1; i<n; ++i)
    for (double j=1; j<n-i; ++j)
    {
      const double c0=(1-(i+j)/n), c1=i/n, c2=j/n;
      *out++=typename Kernel::Point_3(
              p0.x()*c0+p1.x()*c1+p2.x()*c2,
              p0.y()*c0+p1.y()*c1+p2.y()*c2,
              p0.z()*c0+p1.z()*c1+p2.z()*c2
            );
    }
  return out;
}

template <class Kernel, class OutputIterator>
OutputIterator
triangle_grid_sampling(const typename Kernel::Triangle_3& t, double distance, OutputIterator out)
{
  return triangle_grid_sampling<Kernel>(t[0], t[1], t[2], distance, out);
}

} //end of namespace internal

template <class Kernel, class TriangleRange, class OutputIterator>
OutputIterator
sample_triangles(const TriangleRange& triangles, double distance, OutputIterator out)
{
  typedef typename Kernel::Point_3 Point_3;
  typedef typename Kernel::Vector_3 Vector_3;
  typedef typename Kernel::Triangle_3 Triangle_3;

  std::set< std::pair<Point_3, Point_3> > sampled_edges;

  // sample edges but skip endpoints
  BOOST_FOREACH(const Triangle_3& t, triangles)
  {
    for (int i=0;i<3; ++i)
    {
      const Point_3& p0=t[i];
      const Point_3& p1=t[(i+1)%3];
      if ( sampled_edges.insert(CGAL::make_sorted_pair(p0, p1)).second )
      {
        const double d_p0p1 = std::sqrt( CGAL::squared_distance(p0, p1) );

        const double nb_pts = std::ceil( d_p0p1 / distance );
        const Vector_3 step_vec = (p1 - p0) / nb_pts;
        for (double i=1; i<nb_pts; ++i)
        {
          *out++=p0 + step_vec * i;
        }
      }
    }
  }

  // sample triangles
  BOOST_FOREACH(const Triangle_3& t, triangles)
    out=internal::triangle_grid_sampling<Kernel>(t, distance, out);

  //add endpoints
  std::set< Point_3 > endpoints;
  BOOST_FOREACH(const Triangle_3& t, triangles)
    for(int i=0; i<3; ++i)
    {
      if ( endpoints.insert(t[i]).second ) *out++=t[i];
    }
  return out;
}

template <class Kernel, class TriangleRange1, class TriangleRange2>
double approximated_Hausdorff_distance(
  const TriangleRange1& triangles_1,
  const TriangleRange2& triangles_2,
  double targeted_precision
)
{
  std::vector<typename Kernel::Point_3> sample_points;
  sample_triangles<Kernel>( triangles_1,
                            targeted_precision,
                            std::back_inserter(sample_points) );
  // std::ofstream out("/tmp/samples.xyz");
  // std::copy(sample_points.begin(), sample_points.end(), std::ostream_iterator<typename Kernel::Point_3>(out,"\n"));

  typedef typename TriangleRange2::const_iterator Triangle_iterator;
  typedef AABB_triangle_primitive<Kernel, Triangle_iterator> Primitive;
  typedef AABB_traits<Kernel, Primitive> Traits;
  typedef AABB_tree< Traits > Tree;

  Tree tree(triangles_2.begin(), triangles_2.end());
  tree.accelerate_distance_queries();

  double hdist = 0;
  BOOST_FOREACH(const typename Kernel::Point_3& pt, sample_points)
  {
    double d = std::sqrt( tree.squared_distance(pt) );
    // if ( d > 1e-1 ) std::cout << pt << "\n";
    if (d>hdist) hdist=d;
  }

  return hdist;
}

template <class Kernel, class TriangleRange1, class TriangleRange2>
double approximated_symmetric_Hausdorff_distance(
  const TriangleRange1& triangles_1,
  const TriangleRange2& triangles_2,
  double targeted_precision
)
{
  return (std::max)(
    approximated_Hausdorff_distance<Kernel>(triangles_1, triangles_2, targeted_precision),
    approximated_Hausdorff_distance<Kernel>(triangles_1, triangles_2, targeted_precision)
  );
}


///\todo add approximated_Hausdorff_distance(FaceGraph, FaceGraph)
///\todo add approximated_Hausdorff_distance(TriangleRange, FaceGraph)
///\todo add approximated_Hausdorff_distance(FaceGraph, TriangleRange)
///\todo add approximated_Hausdorff_distance(PointRange, FaceGraph)
///\todo add approximated_Hausdorff_distance(PointRange, TriangleRange)
///\todo add barycentric_triangle_sampling(Triangle, PointPerAreaUnit)
///\todo add random_triangle_sampling(Triangle, PointPerAreaUnit)
/// \todo maybe we should use a sampling using epec to avoid being too far from the sampled triangle
  



} } // end of namespace CGAL::Polygon_mesh_processing

/// \endcond

#endif //CGAL_POLYGON_MESH_PROCESSING_DISTANCE_H
