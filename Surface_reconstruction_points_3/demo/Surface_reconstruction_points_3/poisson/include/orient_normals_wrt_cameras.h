// Author: Laurent Saboret

#ifndef ORIENT_NORMALS_WRT_CAMERAS_H
#define ORIENT_NORMALS_WRT_CAMERAS_H

#include "Gyroviz_point_3.h"

#include <CGAL/point_set_processing_assertions.h>

#include <iterator>


// ----------------------------------------------------------------------------
// Private section
// ----------------------------------------------------------------------------
namespace CGALi {


/// Orient a 3D point's normal w.r.t. the position of cameras
/// that reconstructed the point by photogrammetry.
//
// @return true if the orientation is robust.
template < class Gt, ///< Geometric traits class.
           class InputIterator
>
bool
orient_normal_wrt_cameras(const typename Gt::Point_3& p, ///< 3D point position
                          typename Gt::Vector_3& normal, ///< 3D point normal (in and out)
                          InputIterator first_camera,  ///< Iterator over first 3D point's camera
                          InputIterator beyond_camera) ///< Past-the-end iterator
{
    typedef typename Gt::FT       FT;
    typedef typename Gt::Point_3  Point;
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
    if (max_dot_product > 0 &&
        std::abs(max_dot_product) > 0.2588) // oriented iff angle < ~75 degrees
    {
      Vector cp = p - max_camera;
      FT dot = (cp * n);
      if (dot > 0)
         n = -n;
      normal = n;
      return true /* oriented */;
    }
    else // if failure
    {
      normal = n;
      return false /* non oriented */;
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
/// This method modifies the order of input points, and returns 
/// an iterator over the first point with an unoriented normal (see erase-remove idiom).
/// Warning: this method should not be called on sorted containers.
///
/// @commentheading Precondition: normals must be unit vectors.
///
/// @commentheading Template Parameters:
/// @param ForwardIterator value_type must be convertible to Gyroviz_point_3<Kernel>.
/// @param Kernel Geometric traits class. It can be omitted and deduced automatically from the iterator type.
///
/// @return iterator over the first point with an unoriented normal.

// This variant requires the kernel.
template <typename ForwardIterator,
          typename Kernel
>
ForwardIterator
orient_normals_wrt_cameras(
           ForwardIterator first,   ///< iterator over the first input/output point.
           ForwardIterator beyond,  ///< past-the-end iterator.
           const Kernel& kernel)    ///< geometric traits.
{
    CGAL_TRACE("Call orient_normals_wrt_cameras()\n");

    typedef typename std::iterator_traits<ForwardIterator>::value_type Gyroviz_point;

    // Iterate over input points and orient normals.
    // Copy points with robust normal orientation to oriented_points[], the others to unoriented_points[].
    std::deque<Gyroviz_point> oriented_points, unoriented_points;
    for (ForwardIterator it = first; it != beyond; it++)
    {
        bool oriented = CGALi::orient_normal_wrt_cameras<Kernel>(
                                              *it,
                                              it->normal(),
                                              it->cameras_begin(), it->cameras_end());
        if (oriented)
          oriented_points.push_back(*it);
        else
          unoriented_points.push_back(*it);
    }

    // Replace [first, beyond) range by the content of oriented_points[], then unoriented_points[].
    ForwardIterator first_unoriented_point =
      std::copy(oriented_points.begin(), oriented_points.end(), first);
    std::copy(unoriented_points.begin(), unoriented_points.end(), first_unoriented_point);

    CGAL_TRACE("  => %u normals are unoriented\n", unoriented_points.size());
    CGAL_TRACE("End of orient_normals_wrt_cameras()\n");

    return first_unoriented_point;
}

/// @cond SKIP_IN_MANUAL
// This variant deduces the kernel from the iterator type.
template <typename ForwardIterator>
ForwardIterator
orient_normals_wrt_cameras(
       ForwardIterator first,       ///< iterator over the first input/output point
       ForwardIterator beyond)      ///< past-the-end iterator
{
    typedef typename std::iterator_traits<ForwardIterator>::value_type Value_type;
    typedef typename CGAL::Kernel_traits<Value_type>::Kernel Kernel;
    return orient_normals_wrt_cameras(first,beyond,Kernel());
}
/// @endcond


#endif // ORIENT_NORMALS_WRT_CAMERAS_H

