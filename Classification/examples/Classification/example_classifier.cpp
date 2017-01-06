#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Classifier.h>
#include <CGAL/Classification/Point_set_neighborhood.h>
#include <CGAL/Classification/Planimetric_grid.h>
#include <CGAL/Classification/Attribute_base.h>
#include <CGAL/Classification/Attribute/Eigen.h>
#include <CGAL/Classification/Attribute/Distance_to_plane.h>
#include <CGAL/Classification/Attribute/Vertical_dispersion.h>
#include <CGAL/Classification/Attribute/Elevation.h>

#include <CGAL/IO/read_ply_points.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;
typedef Kernel::Iso_cuboid_3 Iso_cuboid_3;
typedef std::vector<Point> Point_range;
typedef CGAL::Identity_property_map<Point> Pmap;

typedef CGAL::Classifier<Point_range, Pmap> Classifier;

typedef CGAL::Classification::Planimetric_grid<Kernel, Point_range, Pmap>       Planimetric_grid;
typedef CGAL::Classification::Point_set_neighborhood<Kernel, Point_range, Pmap> Neighborhood;
typedef CGAL::Classification::Local_eigen_analysis<Kernel, Point_range, Pmap>   Local_eigen_analysis;

typedef CGAL::Classification::Type_handle                                           Type_handle;
typedef CGAL::Classification::Attribute_handle                                      Attribute_handle;

typedef CGAL::Classification::Attribute::Distance_to_plane<Kernel, Point_range, Pmap>   Distance_to_plane;
typedef CGAL::Classification::Attribute::Linearity<Kernel, Point_range, Pmap>           Linearity;
typedef CGAL::Classification::Attribute::Omnivariance<Kernel, Point_range, Pmap>        Omnivariance;
typedef CGAL::Classification::Attribute::Planarity<Kernel, Point_range, Pmap>           Planarity;
typedef CGAL::Classification::Attribute::Surface_variation<Kernel, Point_range, Pmap>   Surface_variation;
typedef CGAL::Classification::Attribute::Elevation<Kernel, Point_range, Pmap>           Elevation;
typedef CGAL::Classification::Attribute::Vertical_dispersion<Kernel, Point_range, Pmap> Dispersion;


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
  //! [Attributes]

  Classifier psc (pts, Pmap());
  std::cerr << "Computing attributes" << std::endl;
  Attribute_handle d2p = psc.add_attribute<Distance_to_plane> (Pmap(), eigen);
  Attribute_handle lin = psc.add_attribute<Linearity> (eigen);
  Attribute_handle omni = psc.add_attribute<Omnivariance> (eigen);
  Attribute_handle plan = psc.add_attribute<Planarity> (eigen);
  Attribute_handle surf = psc.add_attribute<Surface_variation> (eigen);
  Attribute_handle disp = psc.add_attribute<Dispersion> (Pmap(), grid,
                                                         grid_resolution,
                                                         radius_neighbors);
  Attribute_handle elev = psc.add_attribute<Elevation> (Pmap(), grid,
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

  
  //! [Attributes]
  ///////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////
  //! [Classification Types]
  
  std::cerr << "Setting up classification types" << std::endl;
  
  // Create classification type and define how attributes affect them
  Type_handle ground = psc.add_classification_type ("ground");
  ground->set_attribute_effect (d2p,  CGAL::Classification::Attribute::NEUTRAL);
  ground->set_attribute_effect (lin,  CGAL::Classification::Attribute::PENALIZING);
  ground->set_attribute_effect (omni, CGAL::Classification::Attribute::NEUTRAL);
  ground->set_attribute_effect (plan, CGAL::Classification::Attribute::FAVORING);
  ground->set_attribute_effect (surf, CGAL::Classification::Attribute::PENALIZING);
  ground->set_attribute_effect (disp, CGAL::Classification::Attribute::NEUTRAL);
  ground->set_attribute_effect (elev, CGAL::Classification::Attribute::PENALIZING);

  Type_handle vege = psc.add_classification_type ("vegetation");
  vege->set_attribute_effect (d2p,  CGAL::Classification::Attribute::FAVORING);
  vege->set_attribute_effect (lin,  CGAL::Classification::Attribute::NEUTRAL);
  vege->set_attribute_effect (omni, CGAL::Classification::Attribute::FAVORING);
  vege->set_attribute_effect (plan, CGAL::Classification::Attribute::NEUTRAL);
  vege->set_attribute_effect (surf, CGAL::Classification::Attribute::NEUTRAL);
  vege->set_attribute_effect (disp, CGAL::Classification::Attribute::FAVORING);
  vege->set_attribute_effect (elev, CGAL::Classification::Attribute::NEUTRAL);
  
  Type_handle roof = psc.add_classification_type ("roof");
  roof->set_attribute_effect (d2p,  CGAL::Classification::Attribute::NEUTRAL);
  roof->set_attribute_effect (lin,  CGAL::Classification::Attribute::PENALIZING);
  roof->set_attribute_effect (omni, CGAL::Classification::Attribute::FAVORING);
  roof->set_attribute_effect (plan, CGAL::Classification::Attribute::FAVORING);
  roof->set_attribute_effect (surf, CGAL::Classification::Attribute::PENALIZING);
  roof->set_attribute_effect (disp, CGAL::Classification::Attribute::NEUTRAL);
  roof->set_attribute_effect (elev, CGAL::Classification::Attribute::FAVORING);

  //! [Classification Types]
  ///////////////////////////////////////////////////////////////////

  // Run classification
  psc.run_with_graphcut (neighborhood.k_neighbor_query(12), 0.2);
  
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
      
      Type_handle type = psc.classification_type_of (i);
      if (type == ground)
        f << "245 180 0" << std::endl;
      else if (type == vege)
        f << "0 255 27" << std::endl;
      else if (type == roof)
        f << "255 0 170" << std::endl;
      else
        {
          f << "0 0 0" << std::endl;
          std::cerr << "Error: unknown classification type" << std::endl;
        }
    }
  
  std::cerr << "All done" << std::endl;
  return EXIT_SUCCESS;
}
