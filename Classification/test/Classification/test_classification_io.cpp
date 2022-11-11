#if defined (_MSC_VER) && !defined (_WIN64)
#pragma warning(disable:4244) // boost::number_distance::distance()
                              // converts 64 to 32 bits integers
#endif

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Classification.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/Random.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
typedef CGAL::Point_set_3<Point> Point_set;
typedef Point_set::Point_map Point_map;

typedef Kernel::Iso_cuboid_3 Iso_cuboid_3;

namespace Classification = CGAL::Classification;

typedef Classification::Label_handle                                            Label_handle;
typedef Classification::Feature_handle                                          Feature_handle;
typedef Classification::Label_set                                               Label_set;
typedef Classification::Feature_set                                             Feature_set;

typedef Classification::ETHZ::Random_forest_classifier                           Classifier;

typedef Classification::Planimetric_grid<Kernel, Point_set, Point_map>        Planimetric_grid;
typedef Classification::Point_set_neighborhood<Kernel, Point_set, Point_map>  Neighborhood;
typedef Classification::Local_eigen_analysis                                    Local_eigen_analysis;

typedef Classification::Feature::Distance_to_plane<Point_set, Point_map>      Distance_to_plane;
typedef Classification::Feature::Elevation<Kernel, Point_set, Point_map>      Elevation;

int main (int, char**)
{
  Point_set points;

  for (std::size_t i = 0; i < 1000; ++ i)
    points.insert (Point (CGAL::get_default_random().get_double(),
                          CGAL::get_default_random().get_double(),
                          CGAL::get_default_random().get_double()));

  Iso_cuboid_3 bbox = CGAL::bounding_box (points.points().begin(), points.points().end());

  float grid_resolution = 0.34f;
  float radius_dtm = 15.0f;

  Planimetric_grid grid (points, points.point_map(), bbox, grid_resolution);
  Neighborhood neighborhood (points, points.point_map());
  Local_eigen_analysis eigen
    = Local_eigen_analysis::create_from_point_set
    (points, points.point_map(), neighborhood.k_neighbor_query(6));

  Feature_set features;
  Feature_handle distance_to_plane = features.add<Distance_to_plane> (points, points.point_map(), eigen);
  Feature_handle elevation = features.add<Elevation> (points, points.point_map(), grid,
                                                      radius_dtm);

  Label_set labels;

  std::vector<int> training_set (points.size(), -1);
  for (std::size_t i = 0; i < 3; ++ i)
  {
    std::ostringstream oss;
    oss << "label_" << i;
    Label_handle lh = labels.add(oss.str().c_str());

    for (std::size_t j = 0; j < 100; ++ j)
      training_set[std::size_t(CGAL::get_default_random().get_int(0, int(training_set.size())))] = int(i);
  }

  Classifier classifier (labels, features);
  classifier.train (training_set);

  std::ofstream outf ("output_config.gz", std::ios::binary);
  outf.precision(18);
  classifier.save_configuration(outf);
  outf.close();

  Classifier classifier2 (labels, features);
  std::ifstream inf ("output_config.gz", std::ios::binary);
  classifier2.load_configuration(inf);

  Classifier classifier3 (classifier, features);

  std::vector<std::size_t> label_indices (points.size());
  std::vector<std::size_t> label_indices_2 (points.size());
  std::vector<std::size_t> label_indices_3 (points.size());

  Classification::classify<CGAL::Sequential_tag> (points, labels, classifier, label_indices);
  Classification::classify<CGAL::Sequential_tag> (points, labels, classifier2, label_indices_2);
  Classification::classify<CGAL::Sequential_tag> (points, labels, classifier3, label_indices_3);

  assert (label_indices == label_indices_2);
  assert (label_indices == label_indices_3);

  return EXIT_SUCCESS;
}
