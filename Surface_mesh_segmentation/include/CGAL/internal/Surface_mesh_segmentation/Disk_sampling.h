#ifndef CGAL_SURFACE_MESH_SEGMENTATION_DISK_SAMPLING_H
#define CGAL_SURFACE_MESH_SEGMENTATION_DISK_SAMPLING_H
/**
 * @file Disk_sampling.h
 * This file contains 3 sampling methods, which can be used as a template parameter for `SDF_calculation<Polyhedron, DiskSampling>`.
 */
#include <CGAL/number_type_basic.h>

#include <boost/tuple/tuple.hpp>
#include <vector>
#include <cmath>

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
 * @see Disk_sampling.h
 */
class Vogel_disk_sampling
{
protected:
  typedef boost::tuple<double, double, double> Disk_sample;
  typedef std::vector<Disk_sample>             Disk_samples_list;
public:
  /**
   * Samples points from unit-disk.
   * @param number_of_points number of points to be picked
   * @param cone_angle opening angle of cone (might be necessary for weighting)
   * @param[out] samples sampled points from disk, each point is tuple of:
   *   - get<0> = coordinate-x
   *   - get<1> = coordinate-y
   *   - get<2> = weight (proportional to angle between cone-normal)
   */
  void operator()(int number_of_points, double cone_angle,
                  Disk_samples_list& samples) const {
    const double length_of_normal = 1.0 / tan(cone_angle / 2.0);
    const double angle_st_dev = cone_angle / CGAL_ANGLE_ST_DEV_DIVIDER;
    const double golden_ratio = 3.0 - std::sqrt(5.0);

#if 0
    for(int i = 0; i < number_of_points; ++i) {
      double Q = i * golden_ratio * CGAL_PI;
      double R = std::sqrt(static_cast<double>(i) / number_of_points);
      double angle = atan(R / length_of_normal);
      double weight =  exp(-0.5 * (std::pow(angle / angle_st_dev, 2))); // weight
      samples.push_back(Disk_sample(R * cos(Q), R * sin(Q), weight));
    }
#else
    const double custom_power = 8.0 / 8.0;
    for(int i = 0; i < number_of_points; ++i) {
      double Q = i * golden_ratio * CGAL_PI;
      double R = std::pow(static_cast<double>(i) / number_of_points, custom_power);
      // use uniform weigths, since we already give importance to locations that are close to center.
      samples.push_back(Disk_sample(R * cos(Q), R * sin(Q), 1.0));
    }
#endif
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
 */
class Polar_disk_sampling
{
protected:
  typedef boost::tuple<double, double, double> Disk_sample;
  typedef std::vector<Disk_sample>             Disk_samples_list;
public:
  /**
   * Samples points from unit-disk.
   * @param number_of_points number of points to be picked
   * @param cone_angle opening angle of cone (might be necessary for weighting)
   * @param[out] samples sampled points from disk, each point is tuple of:
   *   - get<0> = coordinate-x
   *   - get<1> = coordinate-y
   *   - get<2> = weight (proportional to angle between cone-normal)
   *
   * NOTE:
   * Returned samples size = floor( sqrt (number_of_points) ) ^ 2
   */
  void operator()(int number_of_points, double cone_angle,
                  Disk_samples_list& samples) const {
    const int number_of_points_sqrt = static_cast<int>(std::sqrt(
                                        static_cast<double>(number_of_points)));
    const double length_of_normal = 1.0 / tan(cone_angle / 2.0);
    // use cone_angle / 3 as one standard deviation while weighting.
    const double angle_st_dev = cone_angle / CGAL_ANGLE_ST_DEV_DIVIDER;

    for(int i = 1; i <= number_of_points_sqrt; ++i)
      for(int j = 1; j <= number_of_points_sqrt; ++j) {
        double w1 = static_cast<double>(i) / number_of_points_sqrt;
        double w2 = static_cast<double>(j) / number_of_points_sqrt;
        double R = w1;
        double Q = 2 * w2 * CGAL_PI;
        double angle = atan(R / length_of_normal);
        double weight = exp(-0.5 * (pow(angle / angle_st_dev, 2)));
        samples.push_back(Disk_sample(R * cos(Q), R * sin(Q), weight));
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
 */
class Concentric_disk_sampling
{
protected:
  typedef boost::tuple<double, double, double> Disk_sample;
  typedef std::vector<Disk_sample>             Disk_samples_list;
public:
  /**
   * Samples points from unit-disk.
   * @param number_of_points number of points to be picked
   * @param cone_angle opening angle of cone (might be necessary for weighting)
   * @param[out] samples sampled points from disk, each point is tuple of:
   *   - get<0> = coordinate-x
   *   - get<1> = coordinate-y
   *   - get<2> = weight (proportional to angle between cone-normal)
   *
   * NOTE:
   * Returned samples size = floor( sqrt (number_of_points) ) ^ 2
   */
  void operator()(int number_of_points, double cone_angle,
                  Disk_samples_list& samples) const {
    const int number_of_points_sqrt = static_cast<int>(std::sqrt(
                                        static_cast<double>(number_of_points)));
    const double length_of_normal = 1.0 / tan(cone_angle / 2.0);
    const double fraction = 2.0 / (number_of_points_sqrt -1);
    // use cone_angle / 3 as one standard deviation while weighting.
    const double angle_st_dev = cone_angle / CGAL_ANGLE_ST_DEV_DIVIDER;

    for(int i = 0; i < number_of_points_sqrt; ++i)
      for(int j = 0; j < number_of_points_sqrt; ++j) {
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
        double angle = atan(R / length_of_normal);
        double weight = exp(-0.5 * (pow(angle / angle_st_dev, 2)));
        samples.push_back(Disk_sample(R * cos(Q), R * sin(Q), weight));
      }
  }
};

}//namespace internal
}//namespace CGAL
#undef CGAL_ANGLE_ST_DEV_DIVIDER
#endif //CGAL_SURFACE_MESH_SEGMENTATION_DISK_SAMPLING_H