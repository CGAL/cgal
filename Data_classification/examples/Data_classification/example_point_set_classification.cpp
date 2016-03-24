#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Point_set_classification.h>
#include <CGAL/IO/read_ply_points.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef CGAL::Point_set_classification<Kernel> PSC;

typedef CGAL::Scatter_segmentation_attribute<Kernel> Scatter;
typedef CGAL::Elevation_segmentation_attribute<Kernel> Elevation;
typedef CGAL::Distance_to_plane_segmentation_attribute<Kernel> Planarity;

typedef CGAL::Vegetation_classification_type<Kernel> Vegetation;
typedef CGAL::Ground_classification_type<Kernel> Ground;
typedef CGAL::Roof_classification_type<Kernel> Roof;

typedef Kernel::Point_3 Point;

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

  std::cerr << "Initializing" << std::endl;
  PSC psc (pts.begin (), pts.end(), 1.5);
  psc.initialization();

  std::cerr << "Computing attributes" << std::endl;
  Scatter scat (psc, 0.1);
  Planarity plan (psc, 3.5);
  Elevation elev (psc, 16);

  Vegetation vege (scat, plan, elev);
  Ground ground (scat, plan, elev);
  Roof roof (scat, plan, elev);

  psc.segmentation_classes.push_back (&vege);
  psc.segmentation_classes.push_back (&ground);
  psc.segmentation_classes.push_back (&roof);

  psc.point_cloud_classification(2);
  
  std::cerr << "All done" << std::endl;
  return EXIT_SUCCESS;
}
