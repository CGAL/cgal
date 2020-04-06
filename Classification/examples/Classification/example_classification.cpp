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
#include <CGAL/bounding_box.h>
#include <CGAL/tags.h>
#include <CGAL/IO/read_ply_points.h>
#include <CGAL/IO/write_ply_points.h>

#include <CGAL/Real_timer.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;
typedef Kernel::Iso_cuboid_3 Iso_cuboid_3;
typedef std::vector<Point> Point_range;
typedef CGAL::Identity_property_map<Point> Pmap;

namespace Classification = CGAL::Classification;

typedef Classification::Sum_of_weighted_features_classifier Classifier;

typedef Classification::Planimetric_grid<Kernel, Point_range, Pmap>             Planimetric_grid;
typedef Classification::Point_set_neighborhood<Kernel, Point_range, Pmap>       Neighborhood;
typedef Classification::Local_eigen_analysis                                    Local_eigen_analysis;

typedef Classification::Label_handle                                            Label_handle;
typedef Classification::Feature_handle                                          Feature_handle;
typedef Classification::Label_set                                               Label_set;
typedef Classification::Feature_set                                             Feature_set;

typedef Classification::Feature::Distance_to_plane<Point_range, Pmap>           Distance_to_plane;
typedef Classification::Feature::Elevation<Kernel, Point_range, Pmap>           Elevation;
typedef Classification::Feature::Vertical_dispersion<Kernel, Point_range, Pmap> Dispersion;


///////////////////////////////////////////////////////////////////
//! [Analysis]

int main (int argc, char** argv)
{
  std::string filename (argc > 1 ? argv[1] : "data/b9.ply");
  std::ifstream in (filename.c_str());
  std::vector<Point> pts;

  std::cerr << "Reading input" << std::endl;
  if (!in
      || !(CGAL::read_ply_points (in, std::back_inserter (pts))))
  {
    std::cerr << "Error: cannot read " << filename << std::endl;
    return EXIT_FAILURE;
  }

  float grid_resolution = 0.34f;
  unsigned int number_of_neighbors = 6;

  std::cerr << "Computing useful structures" << std::endl;

  Iso_cuboid_3 bbox = CGAL::bounding_box (pts.begin(), pts.end());

  Planimetric_grid grid (pts, Pmap(), bbox, grid_resolution);
  Neighborhood neighborhood (pts, Pmap());
  Local_eigen_analysis eigen
    = Local_eigen_analysis::create_from_point_set
    (pts, Pmap(), neighborhood.k_neighbor_query(number_of_neighbors));

  //! [Analysis]
  ///////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////
  //! [Features]

  float radius_neighbors = 1.7f;
  float radius_dtm = 15.0f;

  std::cerr << "Computing features" << std::endl;
  Feature_set features;

  features.begin_parallel_additions(); // No effect in sequential mode

  Feature_handle distance_to_plane = features.add<Distance_to_plane> (pts, Pmap(), eigen);
  Feature_handle dispersion = features.add<Dispersion> (pts, Pmap(), grid,
                                                        radius_neighbors);
  Feature_handle elevation = features.add<Elevation> (pts, Pmap(), grid,
                                                      radius_dtm);

  features.end_parallel_additions(); // No effect in sequential mode

  //! [Features]
  ///////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////
  //! [Labels]

  Label_set labels;

  // Init name only
  Label_handle ground = labels.add ("ground");

  // Init name and color
  Label_handle vegetation = labels.add ("vegetation", CGAL::Color(0,255,0));

  // Init name, Color and standard index (here, ASPRS building index)
  Label_handle roof = labels.add ("roof", CGAL::Color (255, 0, 0), 6);

  //! [Labels]
  ///////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////
  //! [Weights]

  std::cerr << "Setting weights" << std::endl;
  Classifier classifier (labels, features);
  classifier.set_weight (distance_to_plane, 6.75e-2f);
  classifier.set_weight (dispersion, 5.45e-1f);
  classifier.set_weight (elevation, 1.47e1f);

  std::cerr << "Setting effects" << std::endl;
  classifier.set_effect (ground, distance_to_plane, Classifier::NEUTRAL);
  classifier.set_effect (ground, dispersion, Classifier::NEUTRAL);
  classifier.set_effect (ground, elevation, Classifier::PENALIZING);

  classifier.set_effect (vegetation, distance_to_plane,  Classifier::FAVORING);
  classifier.set_effect (vegetation, dispersion, Classifier::FAVORING);
  classifier.set_effect (vegetation, elevation, Classifier::NEUTRAL);

  classifier.set_effect (roof, distance_to_plane,  Classifier::NEUTRAL);
  classifier.set_effect (roof, dispersion, Classifier::NEUTRAL);
  classifier.set_effect (roof, elevation, Classifier::FAVORING);

  //! [Weights]
  ///////////////////////////////////////////////////////////////////

  // Run classification
  std::cerr << "Classifying" << std::endl;

  ///////////////////////////////////////////////////////////////////
  //! [Classify]
  std::vector<int> label_indices (pts.size(), -1);

  CGAL::Real_timer t;
  t.start();
  Classification::classify<CGAL::Parallel_if_available_tag> (pts, labels, classifier, label_indices);
  t.stop();
  std::cerr << "Raw classification performed in " << t.time() << " second(s)" << std::endl;
  t.reset();
  //! [Classify]
  ///////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////
  //! [Smoothing]
  t.start();
  Classification::classify_with_local_smoothing<CGAL::Parallel_if_available_tag>
    (pts, Pmap(), labels, classifier,
     neighborhood.sphere_neighbor_query(radius_neighbors),
     label_indices);
  t.stop();
  std::cerr << "Classification with local smoothing performed in " << t.time() << " second(s)" << std::endl;
  t.reset();
  //! [Smoothing]
  ///////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////
  //! [Graph_cut]
  t.start();
  Classification::classify_with_graphcut<CGAL::Parallel_if_available_tag>
    (pts, Pmap(), labels, classifier,
     neighborhood.k_neighbor_query(12),
     0.2f, 4, label_indices);
  t.stop();
  std::cerr << "Classification with graphcut performed in " << t.time() << " second(s)" << std::endl;
  //! [Graph_cut]
  ///////////////////////////////////////////////////////////////////

  // Save the output in a colored PLY format

  std::vector<unsigned char> red, green, blue;
  red.reserve(pts.size());
  green.reserve(pts.size());
  blue.reserve(pts.size());

  for (std::size_t i = 0; i < pts.size(); ++ i)
  {
    Label_handle label = labels[std::size_t(label_indices[i])];
    unsigned r = 0, g = 0, b = 0;
    if (label == ground)
    {
      r = 245; g = 180; b = 0;
    }
    else if (label == vegetation)
    {
      r = 0; g = 255; b = 27;
    }
    else if (label == roof)
    {
      r = 255; g = 0; b = 170;
    }
    red.push_back(r);
    green.push_back(g);
    blue.push_back(b);
  }

  std::ofstream f ("classification.ply");

  CGAL::write_ply_points_with_properties
    (f, CGAL::make_range (boost::counting_iterator<std::size_t>(0),
                          boost::counting_iterator<std::size_t>(pts.size())),
     CGAL::make_ply_point_writer (CGAL::make_property_map(pts)),
     std::make_pair(CGAL::make_property_map(red), CGAL::PLY_property<unsigned char>("red")),
     std::make_pair(CGAL::make_property_map(green), CGAL::PLY_property<unsigned char>("green")),
     std::make_pair(CGAL::make_property_map(blue), CGAL::PLY_property<unsigned char>("blue")));


  std::cerr << "All done" << std::endl;
  return EXIT_SUCCESS;
}
