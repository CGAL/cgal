// Author: Laurent Saboret

#ifndef ORIENT_NORMALS_WRT_CAMERAS_H
#define ORIENT_NORMALS_WRT_CAMERAS_H

#include <CGAL/Orientable_normal_3.h>
#include <CGAL/point_set_processing_assertions.h>

#include <boost/property_map.hpp>

#include <iterator>


// ----------------------------------------------------------------------------
// Private section
// ----------------------------------------------------------------------------
namespace CGALi {


/// Orient a 3D point's normal w.r.t. the position of cameras
/// that reconstructed the point by photogrammetry.
template < class Gt, ///< Geometric traits class.
           class OrientableNormal_3,
           class InputIterator
>
void
orient_normal_wrt_cameras(const typename Gt::Point_3& p, ///< 3D point position
                          OrientableNormal_3& normal, ///< 3D point normal (in and out)
                          InputIterator first_camera,  ///< Iterator over first 3D point's camera
                          InputIterator beyond_camera) ///< Past-the-end iterator
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


} /* namespace CGALi */


// ----------------------------------------------------------------------------
// Public section
// ----------------------------------------------------------------------------


/// Orient the normals of the [first, beyond) range of vertices
/// w.r.t. the position of cameras
/// that reconstructed the points by photogrammetry.
///
/// @commentheading Preconditions:
/// - VertexIterator is a model of ForwardIterator.
/// - VertexPointMap is a model of boost::readable_property_map with a value_type model of Kernel::Point_3.
/// - VertexNormalMap is a model of boost::lvalue_property_map with a value_type model of OrientableNormal_3.
/// - Normals must be unit vectors.
///
/// @return the number of un-oriented normals.

template<class VertexIterator, class VertexPointMap, class VertexNormalMap, class VertexCamerasMap>
unsigned int
orient_normals_wrt_cameras(
    VertexIterator first, ///< iterator over first input vertex
    VertexIterator beyond, ///< past-the-end iterator
    VertexPointMap vertex_point_map, ///< property map VertexIterator -> Point_3
    VertexNormalMap vertex_normal_map, ///< property map VertexIterator -> Normal (in and out)
    VertexCamerasMap vertex_cameras_map) ///< property map VertexIterator -> pair of camera iterators
{
    CGAL_TRACE("Call orient_normals_wrt_cameras()\n");

    typedef typename std::iterator_traits<VertexIterator>::value_type Vertex_type;
    typedef typename Vertex_type::Geom_traits Gt;

    // iterate over input points and orient normals
    for (VertexIterator it = first; it != beyond; it++)
    {
        CGALi::orient_normal_wrt_cameras<Gt>(get(vertex_point_map,it),
                                             vertex_normal_map[it],
                                             get(vertex_cameras_map,it).first, get(vertex_cameras_map,it).second);
    }

    // Count un-oriented normals
    unsigned int unoriented_normals = 0;
    for (VertexIterator it = first; it != beyond; it++)
        if ( ! vertex_normal_map[it].is_oriented() )
          unoriented_normals++;
    CGAL_TRACE("  => %u normals are unoriented\n", unoriented_normals);

    CGAL_TRACE("End of orient_normals_wrt_cameras()\n");

    return unoriented_normals;
}


#endif // ORIENT_NORMALS_WRT_CAMERAS_H

