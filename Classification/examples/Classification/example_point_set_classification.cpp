#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Classifier.h>
#include <CGAL/Classification/Point_set_neighborhood.h>
#include <CGAL/Classification/Planimetric_grid.h>
#include <CGAL/Classification/Attribute.h>
#include <CGAL/Classification/Attributes_eigen.h>
#include <CGAL/Classification/Attribute_distance_to_plane.h>
#include <CGAL/Classification/Attribute_vertical_dispersion.h>
#include <CGAL/Classification/Attribute_elevation.h>

#include <CGAL/IO/read_ply_points.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;
typedef Kernel::Iso_cuboid_3 Iso_cuboid_3;
typedef std::vector<Point>::iterator Iterator;
typedef CGAL::Identity_property_map<Point> Pmap;

typedef CGAL::Classifier<Iterator, Pmap> Classifier;

typedef CGAL::Classification::Planimetric_grid<Kernel, Iterator, Pmap>       Planimetric_grid;
typedef CGAL::Classification::Point_set_neighborhood<Kernel, Iterator, Pmap> Neighborhood;
typedef CGAL::Classification::Local_eigen_analysis<Kernel, Iterator, Pmap>   Local_eigen_analysis;

typedef CGAL::Classification::Type_handle                                           Type_handle;
typedef CGAL::Classification::Attribute_handle                                      Attribute_handle;

typedef CGAL::Classification::Attribute_distance_to_plane<Kernel, Iterator, Pmap>   Distance_to_plane;
typedef CGAL::Classification::Attribute_linearity<Kernel, Iterator, Pmap>           Linearity;
typedef CGAL::Classification::Attribute_omnivariance<Kernel, Iterator, Pmap>        Omnivariance;
typedef CGAL::Classification::Attribute_planarity<Kernel, Iterator, Pmap>           Planarity;
typedef CGAL::Classification::Attribute_surface_variation<Kernel, Iterator, Pmap>   Surface_variation;
typedef CGAL::Classification::Attribute_elevation<Kernel, Iterator, Pmap>           Elevation;
typedef CGAL::Classification::Attribute_vertical_dispersion<Kernel, Iterator, Pmap> Dispersion;


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

  Planimetric_grid grid (pts.begin(), pts.end(), Pmap(), bbox, grid_resolution);
  Neighborhood neighborhood (pts.begin(), pts.end(), Pmap());
  double garbage;
  Local_eigen_analysis eigen (pts.begin(), pts.end(), Pmap(), neighborhood.k_neighbor_query(6), garbage);
  
  //! [Analysis]
  ///////////////////////////////////////////////////////////////////
  
  ///////////////////////////////////////////////////////////////////
  //! [Attributes]

  std::cerr << "Computing attributes" << std::endl;
  Attribute_handle d2p (new Distance_to_plane (pts.begin(), pts.end(), Pmap(), eigen));
  Attribute_handle lin (new Linearity (pts.begin(), pts.end(), eigen));
  Attribute_handle omni (new Omnivariance (pts.begin(), pts.end(), eigen));
  Attribute_handle plan (new Planarity (pts.begin(), pts.end(), eigen));
  Attribute_handle surf (new Surface_variation (pts.begin(), pts.end(), eigen));
  Attribute_handle disp (new Dispersion (pts.begin(), pts.end(), Pmap(), grid,
                                         grid_resolution,
                                         radius_neighbors));
  Attribute_handle elev (new Elevation (pts.begin(), pts.end(), Pmap(), grid,
                                        grid_resolution,
                                        radius_dtm));
  
  std::cerr << "Setting weights" << std::endl;
  d2p->weight  = 6.75e-2;
  lin->weight  = 1.19;
  omni->weight = 1.34e-1;
  plan->weight = 7.32e-1;
  surf->weight = 1.36e-1;
  disp->weight = 5.45e-1;
  elev->weight = 1.47e1;

  
  // Add attributes to classification object
  Classifier psc (pts.begin (), pts.end(), Pmap());
  psc.add_attribute (d2p);
  psc.add_attribute (lin);
  psc.add_attribute (omni);
  psc.add_attribute (plan);
  psc.add_attribute (surf);
  psc.add_attribute (disp);
  psc.add_attribute (elev);

  //! [Attributes]
  ///////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////
  //! [Classification Types]
  
  std::cerr << "Setting up classification types" << std::endl;
  
  // Create classification type and define how attributes affect them
  Type_handle ground = psc.add_classification_type ("ground");
  ground->set_attribute_effect (d2p,  CGAL::Classification::Type::NEUTRAL_ATT);
  ground->set_attribute_effect (lin,  CGAL::Classification::Type::PENALIZED_ATT);
  ground->set_attribute_effect (omni, CGAL::Classification::Type::NEUTRAL_ATT);
  ground->set_attribute_effect (plan, CGAL::Classification::Type::FAVORED_ATT);
  ground->set_attribute_effect (surf, CGAL::Classification::Type::PENALIZED_ATT);
  ground->set_attribute_effect (disp, CGAL::Classification::Type::NEUTRAL_ATT);
  ground->set_attribute_effect (elev, CGAL::Classification::Type::PENALIZED_ATT);

  Type_handle vege = psc.add_classification_type ("vegetation");
  vege->set_attribute_effect (d2p,  CGAL::Classification::Type::FAVORED_ATT);
  vege->set_attribute_effect (lin,  CGAL::Classification::Type::NEUTRAL_ATT);
  vege->set_attribute_effect (omni, CGAL::Classification::Type::FAVORED_ATT);
  vege->set_attribute_effect (plan, CGAL::Classification::Type::NEUTRAL_ATT);
  vege->set_attribute_effect (surf, CGAL::Classification::Type::NEUTRAL_ATT);
  vege->set_attribute_effect (disp, CGAL::Classification::Type::FAVORED_ATT);
  vege->set_attribute_effect (elev, CGAL::Classification::Type::NEUTRAL_ATT);
  
  Type_handle roof = psc.add_classification_type ("roof");
  roof->set_attribute_effect (d2p,  CGAL::Classification::Type::NEUTRAL_ATT);
  roof->set_attribute_effect (lin,  CGAL::Classification::Type::PENALIZED_ATT);
  roof->set_attribute_effect (omni, CGAL::Classification::Type::FAVORED_ATT);
  roof->set_attribute_effect (plan, CGAL::Classification::Type::FAVORED_ATT);
  roof->set_attribute_effect (surf, CGAL::Classification::Type::PENALIZED_ATT);
  roof->set_attribute_effect (disp, CGAL::Classification::Type::NEUTRAL_ATT);
  roof->set_attribute_effect (elev, CGAL::Classification::Type::FAVORED_ATT);

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
