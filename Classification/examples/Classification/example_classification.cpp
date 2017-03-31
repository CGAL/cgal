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

namespace Classif = CGAL::Classification;

typedef Classif::Sum_of_weighted_features_predicate Classification_predicate;

typedef Classif::Planimetric_grid<Kernel, Point_range, Pmap>             Planimetric_grid;
typedef Classif::Point_set_neighborhood<Kernel, Point_range, Pmap>       Neighborhood;
typedef Classif::Local_eigen_analysis<Kernel, Point_range, Pmap>         Local_eigen_analysis;

typedef Classif::Label_handle                                            Label_handle;
typedef Classif::Feature_handle                                          Feature_handle;
typedef Classif::Label_set                                               Label_set;
typedef Classif::Feature_set                                             Feature_set;

typedef Classif::Feature::Distance_to_plane<Kernel, Point_range, Pmap>   Distance_to_plane;
typedef Classif::Feature::Linearity<Kernel, Point_range, Pmap>           Linearity;
typedef Classif::Feature::Omnivariance<Kernel, Point_range, Pmap>        Omnivariance;
typedef Classif::Feature::Planarity<Kernel, Point_range, Pmap>           Planarity;
typedef Classif::Feature::Surface_variation<Kernel, Point_range, Pmap>   Surface_variation;
typedef Classif::Feature::Elevation<Kernel, Point_range, Pmap>           Elevation;
typedef Classif::Feature::Vertical_dispersion<Kernel, Point_range, Pmap> Dispersion;


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

  double grid_resolution = 0.34;
  double radius_neighbors = 1.7;
  double radius_dtm = 15.0;

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
  Feature_handle d2p = features.add<Distance_to_plane> (pts, Pmap(), eigen);
  Feature_handle lin = features.add<Linearity> (pts, eigen);
  Feature_handle omni = features.add<Omnivariance> (pts, eigen);
  Feature_handle plan = features.add<Planarity> (pts, eigen);
  Feature_handle surf = features.add<Surface_variation> (pts, eigen);
  Feature_handle disp = features.add<Dispersion> (pts, Pmap(), grid,
                                                  grid_resolution,
                                                  radius_neighbors);
  Feature_handle elev = features.add<Elevation> (pts, Pmap(), grid,
                                                 grid_resolution,
                                                 radius_dtm);

  //! [Features]
  ///////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////
  //! [Labels]

  Label_set labels;
  Label_handle ground = labels.add ("ground");
  Label_handle vege = labels.add ("vegetation");
  Label_handle roof = labels.add ("roof");

  //! [Labels]
  ///////////////////////////////////////////////////////////////////
  
  ///////////////////////////////////////////////////////////////////
  //! [Weights]

  std::cerr << "Setting weights" << std::endl;
  Classification_predicate predicate (labels, features);
  predicate.set_weight (d2p, 6.75e-2);
  predicate.set_weight (lin, 1.19);
  predicate.set_weight (omni, 1.34e-1);
  predicate.set_weight (plan, 7.32e-1);
  predicate.set_weight (surf, 1.36e-1);
  predicate.set_weight (disp, 5.45e-1);
  predicate.set_weight (elev, 1.47e1);
  
  std::cerr << "Setting effects" << std::endl;
  predicate.set_effect (ground, d2p, Classification_predicate::NEUTRAL);
  predicate.set_effect (ground, lin,  Classification_predicate::PENALIZING);
  predicate.set_effect (ground, omni, Classification_predicate::NEUTRAL);
  predicate.set_effect (ground, plan, Classification_predicate::FAVORING);
  predicate.set_effect (ground, surf, Classification_predicate::PENALIZING);
  predicate.set_effect (ground, disp, Classification_predicate::NEUTRAL);
  predicate.set_effect (ground, elev, Classification_predicate::PENALIZING);
  
  predicate.set_effect (vege, d2p,  Classification_predicate::FAVORING);
  predicate.set_effect (vege, lin,  Classification_predicate::NEUTRAL);
  predicate.set_effect (vege, omni, Classification_predicate::FAVORING);
  predicate.set_effect (vege, plan, Classification_predicate::NEUTRAL);
  predicate.set_effect (vege, surf, Classification_predicate::NEUTRAL);
  predicate.set_effect (vege, disp, Classification_predicate::FAVORING);
  predicate.set_effect (vege, elev, Classification_predicate::NEUTRAL);

  predicate.set_effect (roof, d2p,  Classification_predicate::NEUTRAL);
  predicate.set_effect (roof, lin,  Classification_predicate::PENALIZING);
  predicate.set_effect (roof, omni, Classification_predicate::FAVORING);
  predicate.set_effect (roof, plan, Classification_predicate::FAVORING);
  predicate.set_effect (roof, surf, Classification_predicate::PENALIZING);
  predicate.set_effect (roof, disp, Classification_predicate::NEUTRAL);
  predicate.set_effect (roof, elev, Classification_predicate::FAVORING);

  //! [Weights]
  ///////////////////////////////////////////////////////////////////

  // Run classification
  std::cerr << "Classifying" << std::endl;

  ///////////////////////////////////////////////////////////////////
  //! [Classify]
  std::vector<std::size_t> label_indices;
    
  CGAL::Real_timer t;
  t.start();
  Classif::classify<CGAL::Parallel_tag> (pts, labels, predicate, label_indices);
  t.stop();
  std::cerr << "Raw classification performed in " << t.time() << " second(s)" << std::endl;
  t.reset();
  //! [Classify]
  ///////////////////////////////////////////////////////////////////
  
  ///////////////////////////////////////////////////////////////////
  //! [Smoothing]
  t.start();
  Classif::classify_with_local_smoothing<CGAL::Parallel_tag>
    (pts, Pmap(), labels, predicate,
     neighborhood.range_neighbor_query(radius_neighbors),
     label_indices);
  t.stop();
  std::cerr << "Classification with local smoothing performed in " << t.time() << " second(s)" << std::endl;
  t.reset();
  //! [Smoothing]
  ///////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////
  //! [Graph_cut]
  t.start();
  Classif::classify_with_graphcut<CGAL::Sequential_tag>
    (pts, Pmap(), labels, predicate,
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
    else if (label == vege)
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
