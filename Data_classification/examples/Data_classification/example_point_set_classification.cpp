#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Point_set_classification.h>
#include <CGAL/Data_classification/Planimetric_grid.h>
#include <CGAL/Data_classification/Segmentation_attribute_vertical_dispersion.h>
#include <CGAL/Data_classification/Segmentation_attribute_elevation.h>
#include <CGAL/Data_classification/Segmentation_attribute_verticality.h>
#include <CGAL/Data_classification/Segmentation_attribute_distance_to_plane.h>
#include <CGAL/IO/read_ply_points.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;
typedef Kernel::Iso_cuboid_3 Iso_cuboid_3;
typedef std::vector<Point>::iterator Iterator;
typedef CGAL::Identity_property_map<Point> Pmap;

typedef CGAL::Point_set_classification<Kernel, Iterator, Pmap> Classification;

typedef CGAL::Data_classification::Planimetric_grid<Kernel, Iterator, Pmap>      Planimetric_grid;
typedef CGAL::Data_classification::Neighborhood<Kernel, Iterator, Pmap>          Neighborhood;
typedef CGAL::Data_classification::Local_eigen_analysis<Kernel, Iterator, Pmap>  Local_eigen_analysis;
typedef CGAL::Segmentation_attribute_vertical_dispersion<Kernel, Iterator, Pmap> Dispersion;
typedef CGAL::Segmentation_attribute_elevation<Kernel, Iterator, Pmap>           Elevation;
typedef CGAL::Segmentation_attribute_verticality<Kernel, Iterator, Pmap>         Verticality;
typedef CGAL::Segmentation_attribute_distance_to_plane<Kernel, Iterator, Pmap>   Distance_to_plane;


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

  double grid_resolution = 0.5;
  double radius_neighbors = 1.5;
  double radius_dtm = 12.5;

  std::cerr << "Computing useful structures" << std::endl;

  Iso_cuboid_3 bbox = CGAL::bounding_box (pts.begin(), pts.end());
  Planimetric_grid grid (pts.begin(), pts.end(), Pmap(), bbox, grid_resolution);
  Neighborhood neighborhood (pts.begin(), pts.end(), Pmap());
  Local_eigen_analysis eigen (pts.begin(), pts.end(), Pmap(), neighborhood, radius_neighbors);
  
  std::cerr << "Computing attributes" << std::endl;
  Dispersion disp (pts.begin(), pts.end(), Pmap(), grid,
                   grid_resolution,
                   radius_neighbors,
                   1.78); // Weight
  
  Elevation elev (pts.begin(), pts.end(), Pmap(), bbox, grid,
                  grid_resolution,
                  radius_neighbors,
                  radius_dtm,
                  2.86); // Weight
  
  Verticality verti (pts.begin(), pts.end(), eigen,
                     3.70); // Weight
  
  Distance_to_plane d2p (pts.begin(), pts.end(), Pmap(), eigen,
                         0.0016);

  Classification psc (pts.begin (), pts.end(), Pmap());
  
  // // Add attributes to PSC
  psc.add_segmentation_attribute (&disp);
  psc.add_segmentation_attribute (&elev);
  psc.add_segmentation_attribute (&verti);
  psc.add_segmentation_attribute (&d2p);

  // // Create classification type and define how attributes affect them
  CGAL::Classification_type ground ("ground");
  ground.set_attribute_effect (&disp, CGAL::Classification_type::PENALIZED_ATT);
  ground.set_attribute_effect (&elev, CGAL::Classification_type::PENALIZED_ATT);
  ground.set_attribute_effect (&verti, CGAL::Classification_type::PENALIZED_ATT);
  ground.set_attribute_effect (&d2p, CGAL::Classification_type::PENALIZED_ATT);

  CGAL::Classification_type vege ("vegetation");
  vege.set_attribute_effect (&disp, CGAL::Classification_type::FAVORED_ATT);
  vege.set_attribute_effect (&elev, CGAL::Classification_type::FAVORED_ATT);
  vege.set_attribute_effect (&verti, CGAL::Classification_type::NEUTRAL_ATT);
  vege.set_attribute_effect (&d2p, CGAL::Classification_type::FAVORED_ATT);
  
  CGAL::Classification_type roof ("roof");
  roof.set_attribute_effect (&disp, CGAL::Classification_type::PENALIZED_ATT);
  roof.set_attribute_effect (&elev, CGAL::Classification_type::FAVORED_ATT);
  roof.set_attribute_effect (&verti, CGAL::Classification_type::NEUTRAL_ATT);
  roof.set_attribute_effect (&d2p, CGAL::Classification_type::NEUTRAL_ATT);

  // Add types to PSC
  psc.add_classification_type (&vege);
  psc.add_classification_type (&ground);
  psc.add_classification_type (&roof);

  // Run classification
  psc.run_with_graphcut (neighborhood, 0.5);

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
      
      CGAL::Classification_type* type = psc.classification_type_of (i);
      if (type == &ground)
        f << "245 180 0" << std::endl;
      else if (type == &vege)
        f << "0 255 27" << std::endl;
      else if (type == &roof)
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
