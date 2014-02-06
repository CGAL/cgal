#ifndef CGAL_SURFACE_MESH_SEGMENTATION_DISK_SAMPLERS_H
// Copyright (c) 2014  GeometryFactory Sarl (France).
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


#define CGAL_SURFACE_MESH_SEGMENTATION_DISK_SAMPLERS_H

/// @cond CGAL_DOCUMENT_INTERNAL

/**
 * @file Disk_samplers.h
 * @brief This file contains 3 sampling methods, which can be used as a template parameter for CGAL::internal::SDF_calculation.
 */
#include <cmath>
#include <CGAL/number_type_basic.h>
#include <CGAL/assertions.h>
#include <limits>

#define CGAL_ANGLE_ST_DEV_DIVIDER 3.0

namespace CGAL
{

namespace internal
{

///////////////////////////////////////////////////////// /////////////////////////////////////////////////////////
//                     *                               // //                     *                               //
//                            *                        // //                            *                        //
//                  *     *         *    *             // //                                                     //
//           *                  *                      // //                   *     *       *     *             //
//                *     *                              // //            *                                        //
//                          *       *   *    *         // //                  *          *                       //
//      *    *    *    *       *                       // //      *      *         *            *    *           //
//                         *         *   *       *     // //                          *    *               *     //
//       *    *      *                                 // //         *      *   *  *    *  *    *                //
//               *            *    *          *        // //                       *  *              *           //
//         *            *   *           *       *      // //            *       *    *** *    *          *       //
//     *       *                *    *      *          // //      *          *    *  *     *      *              //
//                  *    *                             // //                     *     *  *   *                  //
//         *                   *        *    *         // //           *    *        *               *           //
//             *    *              *                   // //                      *    *        *                //
//                        *              *             // //                 *             *           *         //
//           *        *       *    *          *        // //            *            *                           //
//               *                                     // //                              *     *                //
//                         *     *      *              // //                     *                               //
//               *    *                                // //                *          *                         //
//                           *       *                 // //                                   *                 //
///////////////////////////////////////////////////////// /////////////////////////////////////////////////////////
// Uniform //                                  // Custom power (biased to center) //
/**
 * @brief Uses Vogel's method to sample points from unit-disk.
 *
 * @tparam Tuple should have a constructor which takes 3 double parameters.
 * @tparam uniform indicates how sampling occurs (uniform or biased to center)
 * @see Disk_samplers.h, SDF_calculation
 */
template<class Tuple, bool uniform = false>
class Vogel_disk_sampling
{
public:
  /**
   * Samples points from unit-disk.
   * @param number_of_points number of points to be picked
   * @param[out] out_it sampled points from disk, each point is tuple of:
   *   - coordinate-x
   *   - coordinate-y
   *   - weight (proportional to distance from disk center)
   */
  template<class OutputIterator>
  void operator()(std::size_t number_of_points,
                  OutputIterator out_it) const {
    const double golden_ratio = 3.0 - std::sqrt(5.0);

    if(uniform) {
      for(std::size_t i = 0; i < number_of_points; ++i) {
        double Q = i * golden_ratio * CGAL_PI;
        double R = std::sqrt(static_cast<double>(i) / number_of_points);
        double weight =  exp(-0.5 * (std::pow(R / CGAL_ANGLE_ST_DEV_DIVIDER, 2)));
        *out_it++ = Tuple(R * cos(Q), R * sin(Q), weight);
      }
    } else {
      const double custom_power =
        1.0; // it determines how much sampling is biased to center
      // for uniform result one can use 0.5 (i.e. sqrt)
      for(std::size_t i = 0; i < number_of_points; ++i) {
        double Q = i * golden_ratio * CGAL_PI;
        double R = std::pow(static_cast<double>(i) / number_of_points, custom_power);
        // use uniform weigths, since we already give importance to locations that are close to center.
        *out_it++ = Tuple(R * cos(Q), R * sin(Q), 1.0);
      }
    }
  }
};

/////////////////////////////////////////////////////////
//                          *                          //
//                          *                          //
//                          *                          //
//          *                               *          //
//            *             *             *            //
//              *           *           *              //
//                *                   *                //
//                  *       *       *                  //
//                    *     *     *                    //
//                      *       *                      //
//                        * * *                        //
//    * *  *  *  *  * *  *     *  * *  *  *  *  * *    //
//                        * * *                        //
//                      *       *                      //
//                    *     *     *                    //
//                  *       *       *                  //
//                *                   *                //
//              *           *           *              //
//            *             *             *            //
//          *                               *          //
//                          *                          //
//                          *                          //
//                          *                          //
/////////////////////////////////////////////////////////
/**
 * @brief Uses polar mapping to sample points from unit-disk.
 *
 * @tparam Tuple should have a constructor which takes 3 double parameters.
 */
template<class Tuple>
class Polar_disk_sampling
{
public:
  /**
   * Samples points from unit-disk.
   * @param number_of_points number of points to be picked
   * @param[out] out_it sampled points from disk, each point is tuple of:
   *   - coordinate-x
   *   - coordinate-y
   *   - weight (proportional to distance from disk center)
   *
   * Note: returned samples size = \f$ \lfloor \sqrt {number\_of\_points} \rfloor ^ 2 \f$
   */
  template<class OutputIterator>
  void operator()(std::size_t number_of_points,
                  OutputIterator out_it) const {
    const std::size_t number_of_points_sqrt = static_cast<std::size_t>(std::sqrt(
          static_cast<double>(number_of_points)));
    // use cone_angle / 3 as one standard deviation while weighting.

    for(std::size_t i = 1; i <= number_of_points_sqrt; ++i)
      for(std::size_t j = 1; j <= number_of_points_sqrt; ++j) {
        double w1 = static_cast<double>(i) / number_of_points_sqrt;
        double w2 = static_cast<double>(j) / number_of_points_sqrt;
        double R = w1;
        double Q = 2 * w2 * CGAL_PI;

        double weight = exp(-0.5 * (pow(R / CGAL_ANGLE_ST_DEV_DIVIDER, 2)));
        *out_it++ = Tuple(R * cos(Q), R * sin(Q), weight);
      }
  }
};

/////////////////////////////////////////////////////////
//                   *   *     *   *                   //
//              *                       *              //
//                                                     //
//          *        *   *     *   *        *          //
//       *                                     *       //
//               *                     *               //
//                        *   *                        //
//     *      *      *             *      *      *     //
//                                                     //
//    *     *      *      *   *      *      *     *    //
//                                                     //
//    *     *      *      *   *      *      *     *    //
//                                                     //
//     *      *      *             *      *      *     //
//                        *   *                        //
//               *                     *               //
//       *                                     *       //
//          *        *   *     *   *        *          //
//                                                     //
//              *                       *              //
//                   *   *     *   *                   //
/////////////////////////////////////////////////////////
/**
 * @brief Uses concentric mapping to sample points from unit-disk.
 *
 * @tparam Tuple should have a constructor which takes 3 double parameters.
 */
template<class Tuple>
class Concentric_disk_sampling
{
public:
  /**
   * Samples points from unit-disk.
   * @param number_of_points number of points to be picked
   * @param cone_angle opening angle of cone (might be necessary for weighting)
   * @param[out] out_it sampled points from disk, each point is tuple of:
   *   - coordinate-x
   *   - coordinate-y
   *   - weight (proportional to angle between cone-normal)
   *
   * Note: returned samples size = \f$ \lfloor \sqrt {number\_of\_points} \rfloor ^ 2 \f$
   */
  template<class OutputIterator>
  void operator()(std::size_t number_of_points,
                  OutputIterator out_it) const {
    const std::size_t number_of_points_sqrt = static_cast<std::size_t>(std::sqrt(
          static_cast<double>(number_of_points)));
    const double fraction = (number_of_points_sqrt == 1) ? 0.0
                            : 2.0 / (number_of_points_sqrt -1);

    for(std::size_t i = 0; i < number_of_points_sqrt; ++i)
      for(std::size_t j = 0; j < number_of_points_sqrt; ++j) {
        double w1 = -1 + i * fraction;
        double w2 = -1 + j * fraction;
        double R, Q;
        if(w1 == 0 && w2 == 0) {
          R = 0;
          Q = 0;
        } else if(w1 > -w2) {
          if(w1 > w2) {
            R = w1;
            Q = (w2 / w1);
          } else        {
            R = w2;
            Q = (2 - w1 / w2);
          }
        } else {
          if(w1 < w2) {
            R = -w1;
            Q = (4 + w2 / w1);
          } else        {
            R = -w2;
            Q = (6 - w1 / w2);
          }
        }
        Q *= (CGAL_PI / 4.0);

        double weight = exp(-0.5 * (pow(R / CGAL_ANGLE_ST_DEV_DIVIDER, 2)));
        *out_it++ = Tuple(R * cos(Q), R * sin(Q), weight);
      }
  }
};

}//namespace internal
/// @endcond
}//namespace CGAL
#undef CGAL_ANGLE_ST_DEV_DIVIDER
#endif //CGAL_SURFACE_MESH_SEGMENTATION_DISK_SAMPLERS_H
