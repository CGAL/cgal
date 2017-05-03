#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Classification.h>
#include <CGAL/bounding_box.h>

#include <CGAL/IO/read_ply_points.h>

#include <CGAL/Real_timer.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;
typedef Kernel::Iso_cuboid_3 Iso_cuboid_3;
typedef std::vector<Point> Point_range;
typedef CGAL::Identity_property_map<Point> Pmap;

namespace Classification = CGAL::Classification;

typedef Classification::Sum_of_weighted_features_classifier Classification_classifier;

typedef Classification::Planimetric_grid<Kernel, Point_range, Pmap>             Planimetric_grid;
typedef Classification::Point_set_neighborhood<Kernel, Point_range, Pmap>       Neighborhood;
typedef Classification::Local_eigen_analysis                                    Local_eigen_analysis;

typedef Classification::Label_handle                                            Label_handle;
typedef Classification::Feature_handle                                          Feature_handle;
typedef Classification::Label_set                                               Label_set;
typedef Classification::Feature_set                                             Feature_set;

typedef Classification::Feature::Distance_to_plane<Point_range, Pmap>           Distance_to_plane;
typedef Classification::Feature::Linearity                                      Linearity;
typedef Classification::Feature::Omnivariance                                   Omnivariance;
typedef Classification::Feature::Planarity                                      Planarity;
typedef Classification::Feature::Surface_variation                              Surface_variation;
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

  float grid_resolution = 0.34;
  float radius_neighbors = 1.7;
  float radius_dtm = 15.0;

  std::cerr << "Computing useful structures" << std::endl;

  Iso_cuboid_3 bbox = CGAL::bounding_box (pts.begin(), pts.end());

  Planimetric_grid grid (pts, Pmap(), bbox, grid_resolution);
  Neighborhood neighborhood (pts, Pmap());
  Local_eigen_analysis eigen (pts, Pmap(), neighborhood.k_neighbor_query(6));

  //! [Analysis]
  ///////////////////////////////////////////////////////////////////
  
  ///////////////////////////////////////////////////////////////////
  //! [Features]

  std::cerr << "Computing features" << std::endl;
  Feature_set features;
  Feature_handle distance_to_plane = features.add<Distance_to_plane> (pts, Pmap(), eigen);
  Feature_handle linearity = features.add<Linearity> (pts, eigen);
  Feature_handle omnivariance = features.add<Omnivariance> (pts, eigen);
  Feature_handle planarity = features.add<Planarity> (pts, eigen);
  Feature_handle surface_variation = features.add<Surface_variation> (pts, eigen);
  Feature_handle dispersion = features.add<Dispersion> (pts, Pmap(), grid,
                                                        grid_resolution,
                                                        radius_neighbors);
  Feature_handle elevation = features.add<Elevation> (pts, Pmap(), grid,
                                                      grid_resolution,
                                                      radius_dtm);

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
  Classification_classifier classifier (labels, features);
  classifier.set_weight (distance_to_plane, 6.75e-2);
  classifier.set_weight (linearity, 1.19);
  classifier.set_weight (omnivariance, 1.34e-1);
  classifier.set_weight (planarity, 7.32e-1);
  classifier.set_weight (surface_variation, 1.36e-1);
  classifier.set_weight (dispersion, 5.45e-1);
  classifier.set_weight (elevation, 1.47e1);
  
  std::cerr << "Setting effects" << std::endl;
  classifier.set_effect (ground, distance_to_plane, Classification_classifier::NEUTRAL);
  classifier.set_effect (ground, linearity,  Classification_classifier::PENALIZING);
  classifier.set_effect (ground, omnivariance, Classification_classifier::NEUTRAL);
  classifier.set_effect (ground, planarity, Classification_classifier::FAVORING);
  classifier.set_effect (ground, surface_variation, Classification_classifier::PENALIZING);
  classifier.set_effect (ground, dispersion, Classification_classifier::NEUTRAL);
  classifier.set_effect (ground, elevation, Classification_classifier::PENALIZING);
  
  classifier.set_effect (vegetation, distance_to_plane,  Classification_classifier::FAVORING);
  classifier.set_effect (vegetation, linearity,  Classification_classifier::NEUTRAL);
  classifier.set_effect (vegetation, omnivariance, Classification_classifier::FAVORING);
  classifier.set_effect (vegetation, planarity, Classification_classifier::NEUTRAL);
  classifier.set_effect (vegetation, surface_variation, Classification_classifier::NEUTRAL);
  classifier.set_effect (vegetation, dispersion, Classification_classifier::FAVORING);
  classifier.set_effect (vegetation, elevation, Classification_classifier::NEUTRAL);

  classifier.set_effect (roof, distance_to_plane,  Classification_classifier::NEUTRAL);
  classifier.set_effect (roof, linearity,  Classification_classifier::PENALIZING);
  classifier.set_effect (roof, omnivariance, Classification_classifier::FAVORING);
  classifier.set_effect (roof, planarity, Classification_classifier::FAVORING);
  classifier.set_effect (roof, surface_variation, Classification_classifier::PENALIZING);
  classifier.set_effect (roof, dispersion, Classification_classifier::NEUTRAL);
  classifier.set_effect (roof, elevation, Classification_classifier::FAVORING);

  //! [Weights]
  ///////////////////////////////////////////////////////////////////

  // Run classification
  std::cerr << "Classifying" << std::endl;

  ///////////////////////////////////////////////////////////////////
  //! [Classify]
  std::vector<std::size_t> label_indices;
    
  CGAL::Real_timer t;
  t.start();
  Classification::classify<CGAL::Parallel_tag> (pts, labels, classifier, label_indices);
  t.stop();
  std::cerr << "Raw classification performed in " << t.time() << " second(s)" << std::endl;
  t.reset();
  //! [Classify]
  ///////////////////////////////////////////////////////////////////
  
  ///////////////////////////////////////////////////////////////////
  //! [Smoothing]
  t.start();
  Classification::classify_with_local_smoothing<CGAL::Parallel_tag>
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
  Classification::classify_with_graphcut<CGAL::Sequential_tag>
    (pts, Pmap(), labels, classifier,
     neighborhood.k_neighbor_query(12),
     0.2, 4, label_indices);
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
      
    Label_handle label = labels[label_indices[i]];
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
