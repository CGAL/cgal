// Author : Nader Salman

#ifndef REMOVE_OUTLIERS_WRT_CAMERA_CONE_ANGLE_3_H
#define REMOVE_OUTLIERS_WRT_CAMERA_CONE_ANGLE_3_H

#include <CGAL/basic.h>
#include "Gyroviz_point_3.h"

#include <iterator>
#include <vector>
#include <map>


/// Compute greatest camera angle for a single point.
///
/// @heading Parameters:
/// @param Kernel Geometric traits class.
/// @param InputIterator value_type is Point_3.
///
/// @return computed greatest camera angle (radians).
template < typename Kernel, typename CameraInputIterator >
double compute_greatest_camera_angle_3(const typename Kernel::Point_3& position, 
                                       CameraInputIterator first_camera, 
                                       CameraInputIterator beyond_camera)
{
    // geometric types
    typedef typename Kernel::FT FT;
    typedef typename Kernel::Point_3 Point;
    typedef typename Kernel::Vector_3 Vector;

    std::vector<Point> cameras;
    std::copy(first_camera, beyond_camera, std::back_inserter(cameras));

    // give a score to each vertex: the score will help detecting outliers			  
    FT greatest_camera_angle=0, v1_v2, n_v1, n_v2, intermediate_score;
    Vector v1, v2;

    for(int i=0; i<cameras.size()-1; ++i)
    {
        for(int j=i+1; j<cameras.size(); ++j)
        {
            v1 = cameras[i] - position;
            v2 = cameras[j] - position;
            n_v1  = sqrt(v1.squared_length());
            n_v2  = sqrt(v2.squared_length()); 
            v1_v2 = v1 * v2; // scalar product
            intermediate_score = acos(v1_v2/(n_v1*n_v2));

            if(intermediate_score > greatest_camera_angle)
                greatest_camera_angle = intermediate_score;
        }
    }
    return greatest_camera_angle;
}



/// Remove vertices / cameras cone's angle < min_cameras_cone_angle.
/// This variant requires the kernel.
///
/// Precondition: min_cameras_cone_angle >= 0.
///
/// @heading Parameters:
/// @param InputIterator value_type is Gyroviz_point.
/// @param OutputIterator value_type is Gyroviz_point.
/// @param Kernel Geometric traits class.
///
/// @return past-the-end output iterator.
template <typename InputIterator,
typename OutputIterator,
typename Kernel
>
OutputIterator
remove_outliers_wrt_camera_cone_angle_3(InputIterator first,            ///< input points 
                                        InputIterator beyond,
                                        OutputIterator output,          ///< output points
                                        double min_cameras_cone_angle,  ///< min angle of camera's cone (radians)
                                        const Kernel& /*kernel*/)
{
    // geometric types
    typedef typename Kernel::FT                FT;
    typedef typename Kernel::Point_3           Point;
    typedef Gyroviz_point_3<Kernel>            Gyroviz_point_3;
    

  
    // precondition: at least one element in the container.
    // to fix: should have at least three distinct points
    // but this is costly to check
    CGAL_precondition(first != beyond);

    // precondition: at least 0
    CGAL_precondition(min_cameras_cone_angle >= 0);

    // iterate over input points
    for(InputIterator point_it = first; point_it != beyond; point_it++)
    {
        FT greatest_camera_angle = compute_greatest_camera_angle_3<Kernel>(*point_it, point_it->cameras_begin(), point_it->cameras_end());
        if (greatest_camera_angle >= min_cameras_cone_angle)
            *output++ = *point_it;
    }

    return output;
}

/// Remove vertices / cameras cone's angle < min_cameras_cone_angle.
/// This function is mutating the input point set.
/// This variant requires the kernel.
///
/// Precondition: min_cameras_cone_angle >= 0.
///
/// @heading Parameters:
/// @param ForwardIterator value_type is Point_3.
/// @param Kernel Geometric traits class.
template <typename ForwardIterator,
typename Kernel
>
void
remove_outliers_wrt_camera_cone_angle_3(ForwardIterator first,          ///< input/output points
                                        ForwardIterator beyond,
                                        double min_cameras_cone_angle,  ///< min angle of camera's cone (radians)
                                        const Kernel& /*kernel*/)
{
    CGAL_precondition(false); // nyi
}


/// Remove vertices / cameras cone's angle < min_cameras_cone_angle.
/// This variant deduces the kernel from iterator types.
///
/// Precondition: min_cameras_cone_angle >= 0.
///
/// @heading Parameters:
/// @param InputIterator value_type is Point_3.
/// @param OutputIterator value_type is Point_3.
///
/// @return past-the-end output iterator.
template <typename InputIterator,
typename OutputIterator
>
OutputIterator
remove_outliers_wrt_camera_cone_angle_3(InputIterator first,            ///< input points
                                        InputIterator beyond,
                                        OutputIterator output,          ///< output points
                                        double min_cameras_cone_angle)  ///< min angle of camera's cone (radians)
{
    typedef typename std::iterator_traits<InputIterator>::value_type Value_type;
    typedef typename CGAL::Kernel_traits<Value_type>::Kernel Kernel;
    return remove_outliers_wrt_camera_cone_angle_3(first,beyond,output,min_cameras_cone_angle,Kernel());
}

/// Remove vertices / cameras cone's angle < min_cameras_cone_angle.
/// This function is mutating the input point set.
/// This variant deduces the kernel from iterator types.
///
/// Precondition: min_cameras_cone_angle >= 0.
///
/// @heading Parameters:
/// @param ForwardIterator value_type is Point_3.
template <typename ForwardIterator>
void
remove_outliers_wrt_camera_cone_angle_3(ForwardIterator first,          ///< input/output points
                                        ForwardIterator beyond,
                                        double min_cameras_cone_angle)  ///< min angle of camera's cone (radians)
{
    typedef typename std::iterator_traits<ForwardIterator>::value_type Value_type;
    typedef typename CGAL::Kernel_traits<Value_type>::Kernel Kernel;
    remove_outliers_wrt_camera_cone_angle_3(first,beyond,min_cameras_cone_angle,Kernel());
}


#endif // REMOVE_OUTLIERS_WRT_CAMERA_CONE_ANGLE_3_H

