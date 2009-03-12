// Author: Laurent Saboret

#ifndef CGAL_NORMAL_ORIENTATION_WRT_CAMERAS_H
#define CGAL_NORMAL_ORIENTATION_WRT_CAMERAS_H

#include <CGAL/Orientable_normal_3.h>

#include <iterator>


/// Orient a 3D point's normal w.r.t. the position of cameras
/// that reconstructed the point by photogrammetry.
template < class Gt, ///< Geometric traits class.
           class OrientableNormal_3,
           class InputIterator
>
void
normal_orientation_wrt_cameras(const typename Gt::Point_3& p, ///< 3D point position
                               OrientableNormal_3& normal, ///< 3D point normal (in and out)
                               InputIterator first_camera,  ///< 3D point cameras
                               InputIterator beyond_camera) ///< 3D point cameras
{
    typedef typename Gt::FT       FT;
    typedef typename Gt::Point_3  Point;
    typedef OrientableNormal_3    Normal;
    typedef typename Gt::Vector_3 Vector;

    Vector n = normal;

    //                                        ->
    //                                        cp     ->
    // Search for the camera position c / | ------ * n | is maximum
    //                                      ||cp||
    FT max_dot_product = 0;
    Point max_camera;
    for(InputIterator c = first_camera; c != beyond_camera; c++)
    {
    Vector cp = p - *c;
    FT norm_cp = std::sqrt(cp * cp);
    if (norm_cp != 0.0)
    {
      FT dot = (cp * n) / norm_cp;
      if (dot < 0)
        dot = -dot;

      if (max_dot_product < dot)
      {
         max_dot_product = dot;
         max_camera = *c;
      }
    }
    }

    //        ->             ->
    // Orient n backwards to cp
    if (max_dot_product > 0)
    {
      Vector cp = p - max_camera;
      FT dot = (cp * n);
      if (dot > 0)
         n = -n;
      normal = Normal(n, true /* oriented */);
    }
    else // if failure
    {
      normal = Normal(n, false /* non oriented */);
    }
}


/// Orient the normals of the [first, beyond) range of vertices
/// w.r.t. the position of cameras
/// that reconstructed the points by photogrammetry.
///
/// @commentheading Preconditions:
/// - VertexIterator is a model of ForwardIterator.
/// - VertexPointMap is a model of boost::readable_property_map.
/// - VertexNormalMap is a model of boost::lvalue_property_map.
/// - VertexCamerasMap is a model of boost::readable_property_map.
/// - Normals must be unit vectors.
template < class VertexIterator, class VertexPointMap, class VertexNormalMap, class VertexCamerasMap >
void
normal_orientation_wrt_cameras(VertexIterator first, ///< range of vertices
                               VertexIterator beyond, ///< range of vertices
                               VertexPointMap vertex_point_map, ///< property map VertexIterator -> Point_3
                               VertexNormalMap vertex_normal_map, ///< property map VertexIterator -> Normal (in and out)
                               VertexCamerasMap vertex_cameras_map) ///< property map VertexIterator -> pair of camera iterators
{
  typedef typename std::iterator_traits<VertexIterator>::value_type Vertex_type;
  typedef typename Vertex_type::Geom_traits Gt;

  // iterate over input points and orient normals
  for (VertexIterator it = first; it != beyond; it++)
  {
      normal_orientation_wrt_cameras<Gt>(get(vertex_point_map,it),
                                      vertex_normal_map[it],
                                      get(vertex_cameras_map,it).first, get(vertex_cameras_map,it).second);
  }
}


#endif // CGAL_NORMAL_ORIENTATION_WRT_CAMERAS_H

