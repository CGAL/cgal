#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Classifier.h>
#include <CGAL/Classification/Point_set_neighborhood.h>
#include <CGAL/Classification/Planimetric_grid.h>
#include <CGAL/Classification/Feature_base.h>
#include <CGAL/Classification/Feature/Eigen.h>
#include <CGAL/Classification/Feature/Distance_to_plane.h>
#include <CGAL/Classification/Feature/Vertical_dispersion.h>
#include <CGAL/Classification/Feature/Elevation.h>

#include <CGAL/IO/read_ply_points.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;
typedef Kernel::Iso_cuboid_3 Iso_cuboid_3;
typedef std::vector<Point> Point_range;
typedef CGAL::Identity_property_map<Point> Pmap;

typedef CGAL::Classifier<Point_range, Pmap> Classifier;

typedef CGAL::Classification::Planimetric_grid<Kernel, Point_range, Pmap>             Planimetric_grid;
typedef CGAL::Classification::Point_set_neighborhood<Kernel, Point_range, Pmap>       Neighborhood;
typedef CGAL::Classification::Local_eigen_analysis<Kernel, Point_range, Pmap>         Local_eigen_analysis;

typedef CGAL::Classification::Label_handle                                            Label_handle;
typedef CGAL::Classification::Feature_handle                                          Feature_handle;

typedef CGAL::Classification::Feature::Distance_to_plane<Kernel, Point_range, Pmap>   Distance_to_plane;
typedef CGAL::Classification::Feature::Linearity<Kernel, Point_range, Pmap>           Linearity;
typedef CGAL::Classification::Feature::Omnivariance<Kernel, Point_range, Pmap>        Omnivariance;
typedef CGAL::Classification::Feature::Planarity<Kernel, Point_range, Pmap>           Planarity;
typedef CGAL::Classification::Feature::Surface_variation<Kernel, Point_range, Pmap>   Surface_variation;
typedef CGAL::Classification::Feature::Elevation<Kernel, Point_range, Pmap>           Elevation;
typedef CGAL::Classification::Feature::Vertical_dispersion<Kernel, Point_range, Pmap> Dispersion;


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

  Classifier classifier (pts, Pmap());

  //! [Analysis]
  ///////////////////////////////////////////////////////////////////
  
  ///////////////////////////////////////////////////////////////////
  //! [Features]

  std::cerr << "Computing features" << std::endl;
  Feature_handle d2p = classifier.add_feature<Distance_to_plane> (Pmap(), eigen);
  Feature_handle lin = classifier.add_feature<Linearity> (eigen);
  Feature_handle omni = classifier.add_feature<Omnivariance> (eigen);
  Feature_handle plan = classifier.add_feature<Planarity> (eigen);
  Feature_handle surf = classifier.add_feature<Surface_variation> (eigen);
  Feature_handle disp = classifier.add_feature<Dispersion> (Pmap(), grid,
                                                                grid_resolution,
                                                                radius_neighbors);
  Feature_handle elev = classifier.add_feature<Elevation> (Pmap(), grid,
                                                               grid_resolution,
                                                               radius_dtm);
  
  std::cerr << "Setting weights" << std::endl;
  d2p->set_weight(6.75e-2);
  lin->set_weight(1.19);
  omni->set_weight(1.34e-1);
  plan->set_weight(7.32e-1);
  surf->set_weight(1.36e-1);
  disp->set_weight(5.45e-1);
  elev->set_weight(1.47e1);

  
  //! [Features]
  ///////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////
  //! [Classification Labels]
  
  std::cerr << "Setting up labels" << std::endl;
  
  // Create label and define how features affect them
  Label_handle ground = classifier.add_label ("ground");
  ground->set_feature_effect (d2p,  CGAL::Classification::Feature::NEUTRAL);
  ground->set_feature_effect (lin,  CGAL::Classification::Feature::PENALIZING);
  ground->set_feature_effect (omni, CGAL::Classification::Feature::NEUTRAL);
  ground->set_feature_effect (plan, CGAL::Classification::Feature::FAVORING);
  ground->set_feature_effect (surf, CGAL::Classification::Feature::PENALIZING);
  ground->set_feature_effect (disp, CGAL::Classification::Feature::NEUTRAL);
  ground->set_feature_effect (elev, CGAL::Classification::Feature::PENALIZING);

  Label_handle vege = classifier.add_label ("vegetation");
  vege->set_feature_effect (d2p,  CGAL::Classification::Feature::FAVORING);
  vege->set_feature_effect (lin,  CGAL::Classification::Feature::NEUTRAL);
  vege->set_feature_effect (omni, CGAL::Classification::Feature::FAVORING);
  vege->set_feature_effect (plan, CGAL::Classification::Feature::NEUTRAL);
  vege->set_feature_effect (surf, CGAL::Classification::Feature::NEUTRAL);
  vege->set_feature_effect (disp, CGAL::Classification::Feature::FAVORING);
  vege->set_feature_effect (elev, CGAL::Classification::Feature::NEUTRAL);
  
  Label_handle roof = classifier.add_label ("roof");
  roof->set_feature_effect (d2p,  CGAL::Classification::Feature::NEUTRAL);
  roof->set_feature_effect (lin,  CGAL::Classification::Feature::PENALIZING);
  roof->set_feature_effect (omni, CGAL::Classification::Feature::FAVORING);
  roof->set_feature_effect (plan, CGAL::Classification::Feature::FAVORING);
  roof->set_feature_effect (surf, CGAL::Classification::Feature::PENALIZING);
  roof->set_feature_effect (disp, CGAL::Classification::Feature::NEUTRAL);
  roof->set_feature_effect (elev, CGAL::Classification::Feature::FAVORING);

  //! [Classification Labels]
  ///////////////////////////////////////////////////////////////////

  // Run classification
  classifier.run_with_graphcut (neighborhood.k_neighbor_query(12), 0.2);
  
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
      
      Label_handle label = classifier.label_of (i);
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
