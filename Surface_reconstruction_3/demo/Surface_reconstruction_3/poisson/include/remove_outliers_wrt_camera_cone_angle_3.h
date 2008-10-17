// Author : Nader Salman

#ifndef REMOVE_OUTLIERS_WRT_CAMERA_CONE_ANGLE_3_H
#define REMOVE_OUTLIERS_WRT_CAMERA_CONE_ANGLE_3_H

#include "Gyroviz_point_3.h"

#include <iterator>
#include <vector>
#include <map>


/// Compute cameras' cone angle for a single point.
///
/// @heading Parameters:
/// @param Kernel Geometric traits class.
///
/// @return computed greatest camera angle (radians).
template < typename Kernel >
double compute_greatest_camera_angle_3(const Gyroviz_point_3<Kernel>& gpt)
{
    // geometric types
    typedef typename Kernel::FT FT;
    typedef typename Kernel::Point_3 Point;
    typedef typename Kernel::Vector_3 Vector;
    
    assert(gpt.cameras_begin() != gpt.cameras_end());

    std::vector<Point> cameras(gpt.cameras_begin(), gpt.cameras_end());

    // give a score to each vertex: the score will help detecting outliers			  
    FT greatest_camera_angle=0, v1_v2, n_v1, n_v2, intermediate_score;
    Vector v1, v2;
    for(int i=0; i<cameras.size()-1; ++i)
    {
        for(int j=i+1; j<cameras.size(); ++j)
        {
            v1 = cameras[i] - gpt;
            v2 = cameras[j] - gpt;
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

/// Utility class for remove_outliers_wrt_camera_cone_angle_3():
/// Predicate that indicates if cameras cone's angle >= min_camera_cone_angle.
template <typename Kernel>
struct Is_cameras_cone_angle_greater_equal
{
    typedef typename Kernel::FT FT;

private:
    double m_min_camera_cone_angle;

public:
    Is_cameras_cone_angle_greater_equal(FT min_camera_cone_angle) ///< min angle of camera's cone (radians)
      : m_min_camera_cone_angle (min_camera_cone_angle) {}

    bool operator() (const Gyroviz_point_3<Kernel>& gpt) const 
    {
        FT greatest_camera_angle = compute_greatest_camera_angle_3<Kernel>(gpt);
        return (greatest_camera_angle >= m_min_camera_cone_angle);
    }
};

/// Remove points / cameras cone's angle < min_camera_cone_angle.
/// This variant requires the kernel.
///
/// Precondition: min_camera_cone_angle >= 0.
///
/// @heading Parameters:
/// @param InputIterator value_type must be Gyroviz_point_3.
/// @param OutputIterator value_type must be Gyroviz_point_3.
/// @param Kernel Geometric traits class.
///
/// @return past-the-end output iterator.
template <typename InputIterator,
          typename OutputIterator,
          typename Kernel
>
OutputIterator
remove_outliers_wrt_camera_cone_angle_3(InputIterator first,           ///< input points 
                                        InputIterator beyond,
                                        OutputIterator output,         ///< output points
                                        double min_camera_cone_angle,  ///< min angle of camera's cone (radians)
                                        const Kernel& /*kernel*/)
{
    // geometric types
    typedef typename Kernel::FT                FT;
    typedef typename Kernel::Point_3           Point;
    typedef Gyroviz_point_3<Kernel>            Gyroviz_point_3;
    
    // precondition: at least 0
    CGAL_precondition(min_camera_cone_angle >= 0);

    // iterate over input points and output points / camera angle >= min_camera_cone_angle
    for(InputIterator point_it = first; point_it != beyond; point_it++)
    {
        FT greatest_camera_angle = compute_greatest_camera_angle_3<Kernel>(*point_it);
        if (greatest_camera_angle >= min_camera_cone_angle)
            *output++ = *point_it;
    }
    return output;
}

/// Remove points / cameras cone's angle < min_camera_cone_angle.
/// This function is mutating the input point set.
/// This variant requires the kernel.
///
/// Warning: 
/// This method modifies the order of points, thus
//  should not be called on sorted containers.
///
/// Precondition: min_camera_cone_angle >= 0.
///
/// @heading Parameters:
/// @param ForwardIterator value_type must be Gyroviz_point_3.
/// @param Kernel Geometric traits class.
///
/// @return First iterator to remove (see erase-remove idiom).
template <typename ForwardIterator,
          typename Kernel
>
ForwardIterator
remove_outliers_wrt_camera_cone_angle_3(ForwardIterator first,         ///< input/output points
                                        ForwardIterator beyond,
                                        double min_camera_cone_angle,  ///< min angle of camera's cone (radians)
                                        const Kernel& /*kernel*/)
{
    // geometric types
    typedef typename Kernel::FT                FT;
    typedef typename Kernel::Point_3           Point;
    typedef Gyroviz_point_3<Kernel>            Gyroviz_point_3;
    
    // precondition: at least 0
    CGAL_precondition(min_camera_cone_angle >= 0);

    ForwardIterator first_iterator_to_remove =
      std::partition(first, beyond, Is_cameras_cone_angle_greater_equal<Kernel>(min_camera_cone_angle));

    return first_iterator_to_remove;
}

/// Remove points / cameras cone's angle < min_camera_cone_angle.
/// This variant deduces the kernel from iterator types.
///
/// Precondition: min_camera_cone_angle >= 0.
///
/// @heading Parameters:
/// @param InputIterator value_type must be Gyroviz_point_3.
/// @param OutputIterator value_type must be Gyroviz_point_3.
///
/// @return past-the-end output iterator.
template <typename InputIterator,
          typename OutputIterator
>
OutputIterator
remove_outliers_wrt_camera_cone_angle_3(InputIterator first,           ///< input points
                                        InputIterator beyond,
                                        OutputIterator output,         ///< output points
                                        double min_camera_cone_angle)  ///< min angle of camera's cone (radians)
{
    typedef typename std::iterator_traits<InputIterator>::value_type Value_type;
    typedef typename CGAL::Kernel_traits<Value_type>::Kernel Kernel;
    return remove_outliers_wrt_camera_cone_angle_3(first,beyond,output,min_camera_cone_angle,Kernel());
}

/// Remove points / cameras cone's angle < min_camera_cone_angle.
/// This function is mutating the input point set.
/// This variant deduces the kernel from iterator types.
///
/// Warning: 
/// This method modifies the order of points, thus
//  should not be called on sorted containers.
///
/// Precondition: min_camera_cone_angle >= 0.
///
/// @heading Parameters:
/// @param ForwardIterator value_type must be Gyroviz_point_3.
///
/// @return First iterator to remove (see erase-remove idiom).
template <typename ForwardIterator>
ForwardIterator
remove_outliers_wrt_camera_cone_angle_3(ForwardIterator first,         ///< input/output points
                                        ForwardIterator beyond,
                                        double min_camera_cone_angle)  ///< min angle of camera's cone (radians)
{
    typedef typename std::iterator_traits<ForwardIterator>::value_type Value_type;
    typedef typename CGAL::Kernel_traits<Value_type>::Kernel Kernel;
    return remove_outliers_wrt_camera_cone_angle_3(first,beyond,min_camera_cone_angle,Kernel());
}


#endif // REMOVE_OUTLIERS_WRT_CAMERA_CONE_ANGLE_3_H

