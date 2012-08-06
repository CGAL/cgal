#ifndef CGAL_SURFACE_MESH_SEGMENTATION_DISK_SAMPLING_H
#define CGAL_SURFACE_MESH_SEGMENTATION_DISK_SAMPLING_H

#include <boost/tuple/tuple.hpp>
#include <vector>

#define CGAL_ANGLE_ST_DEV_DIVIDER 2.0

namespace CGAL
{
namespace internal
{
/**
 * @brief Uses Vogel's method to sample points from unit-disk.
 */
class Vogel_disk_sampling
{
protected:
  typedef boost::tuple<double, double, double> Disk_sample;
  typedef std::vector<Disk_sample>             Disk_samples_list;
public:
  Vogel_disk_sampling() { }
  /**
   * Samples points from unit-disk.
   * @param number_of_points number of points to be picked
   * @param cone_angle opening angle of cone (might be necessary for weighting)
   * @param[out] samples sampled points from disk, each point is tuple of:
   *   - get<0> = coordinate-x
   *   - get<1> = coordinate-y
   *   - get<2> = weight (proportional to distance from origin (another words angle between cone-normal)).
   */
  void sample(int number_of_points, double cone_angle,
              Disk_samples_list& samples) {
    const double length_of_normal = 1.0 / tan(cone_angle / 2.0);
    const double angle_st_dev = cone_angle / CGAL_ANGLE_ST_DEV_DIVIDER;
    const double golden_ratio = 3.0 - std::sqrt(5.0);

#if 0
    for(int i = 0; i < number_of_points; ++i) {
      double Q = i * golden_ratio * CGAL_PI;
      double R = std::sqrt(static_cast<double>(i) / number_of_points);
      double angle = atan(R / length_of_normal);
      angle =  exp(-0.5 * (std::pow(angle / angle_st_dev, 2))); // weight
      samples.push_back(Disk_sample(R * cos(Q), R * sin(Q), angle));
    }
#else
    double custom_power = 8.0 / 8.0;
    for(int i = 0; i < number_of_points; ++i) {
      double Q = i * golden_ratio * CGAL_PI;
      double R = std::pow(static_cast<double>(i) / number_of_points, custom_power);
      samples.push_back(Disk_sample(R * cos(Q), R * sin(Q), 1.0));
    }
#endif
  }
};
}//namespace internal
}//namespace CGAL
#undef CGAL_ANGLE_ST_DEV_DIVIDER
#endif //CGAL_SURFACE_MESH_SEGMENTATION_DISK_SAMPLING_H