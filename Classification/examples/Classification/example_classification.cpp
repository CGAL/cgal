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

#include <CGAL/IO/read_ply_points.h>

#include <CGAL/Real_timer.h>

#ifdef CGAL_LINKED_WITH_TBB
typedef CGAL::Parallel_tag Concurrency_tag;
#else
typedef CGAL::Sequential_tag Concurrency_tag;
#endif

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

#ifdef CGAL_LINKED_WITH_TBB
  features.begin_parallel_additions();
#endif
  
  Feature_handle distance_to_plane = features.add<Distance_to_plane> (pts, Pmap(), eigen);
  Feature_handle dispersion = features.add<Dispersion> (pts, Pmap(), grid,
                                                        radius_neighbors);
  Feature_handle elevation = features.add<Elevation> (pts, Pmap(), grid,
                                                      radius_dtm);

#ifdef CGAL_LINKED_WITH_TBB
  features.end_parallel_additions();
#endif
  
  //! [Features]
  ///////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////
  //! [Labels]

  Label_set labels;
  Label_handle ground = labels.add ("ground");
  Label_handle vegetation = labels.add ("vegetation");
  Label_handle roof = labels.add ("roof");

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
  Classification::classify<Concurrency_tag> (pts, labels, classifier, label_indices);
  t.stop();
  std::cerr << "Raw classification performed in " << t.time() << " second(s)" << std::endl;
  t.reset();
  //! [Classify]
  ///////////////////////////////////////////////////////////////////
  
  ///////////////////////////////////////////////////////////////////
  //! [Smoothing]
  t.start();
  Classification::classify_with_local_smoothing<Concurrency_tag>
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
  Classification::classify_with_graphcut<Concurrency_tag>
    (pts, Pmap(), labels, classifier,
     neighborhood.k_neighbor_query(12),
     0.2f, 4, label_indices);
  t.stop();
  std::cerr << "Classification with graphcut performed in " << t.time() << " second(s)" << std::endl;
  //! [Graph_cut]
  ///////////////////////////////////////////////////////////////////
  
  // Save the output in a colored PLY format

  std::ofstream f ("classification.ply");
  f << "ply" << std::endl
    << "format ascii 1.0" << std::endl
    << "element vertex " << pts.size() << std::endl
    << "property float x" << std::endl
    << "property float y" << std::endl
    << "property float z" << std::endl
    << "property uchar red" << std::endl
    << "property uchar green" << std::endl
    << "property uchar blue" << std::endl
    << "end_header" << std::endl;
  
  for (std::size_t i = 0; i < pts.size(); ++ i)
  {
    f << pts[i] << " ";
      
    Label_handle label = labels[std::size_t(label_indices[i])];
    if (label == ground)
      f << "245 180 0" << std::endl;
    else if (label == vegetation)
      f << "0 255 27" << std::endl;
    else if (label == roof)
      f << "255 0 170" << std::endl;
    else
    {
      f << "0 0 0" << std::endl;
      std::cerr << "Error: unknown classification label" << std::endl;
    }
  }
  
  std::cerr << "All done" << std::endl;
  return EXIT_SUCCESS;
}
