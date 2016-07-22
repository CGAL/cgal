#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Point_set_classification.h>
#include <CGAL/IO/read_ply_points.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef CGAL::Point_set_classification<Kernel> PSC;

typedef CGAL::Segmentation_attribute_scatter<Kernel> Scatter;
typedef CGAL::Segmentation_attribute_elevation<Kernel> Elevation;
typedef CGAL::Segmentation_attribute_horizontality<Kernel> Horizontality;
typedef CGAL::Segmentation_attribute_nonplanarity<Kernel> NonPlanarity;

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
  PSC psc (pts.begin (), pts.end(), 0.8);
  psc.initialization();

  std::cerr << "Computing attributes" << std::endl;
  // Attributes with user-defined weights
  Scatter scat (psc, 0.26);
  Elevation elev (psc, 0.08);
  Horizontality horiz (psc, 0.13);
  NonPlanarity plan (psc, 0.72);

  // Add attributes to PSC
  psc.add_segmentation_attribute (&scat);
  psc.add_segmentation_attribute (&elev);
  psc.add_segmentation_attribute (&horiz);
  psc.add_segmentation_attribute (&plan);

  // Create classification type and define how attributes affect them
  CGAL::Classification_type ground ("ground");
  ground.set_attribute_effect (&scat, CGAL::Classification_type::PENALIZED_ATT);
  ground.set_attribute_effect (&elev, CGAL::Classification_type::PENALIZED_ATT);
  ground.set_attribute_effect (&horiz, CGAL::Classification_type::PENALIZED_ATT);
  ground.set_attribute_effect (&plan, CGAL::Classification_type::PENALIZED_ATT);

  CGAL::Classification_type vege ("vegetation");
  vege.set_attribute_effect (&scat, CGAL::Classification_type::FAVORED_ATT);
  vege.set_attribute_effect (&elev, CGAL::Classification_type::NEUTRAL_ATT);
  vege.set_attribute_effect (&horiz, CGAL::Classification_type::NEUTRAL_ATT);
  vege.set_attribute_effect (&plan, CGAL::Classification_type::PENALIZED_ATT);
  
  CGAL::Classification_type roof ("roof");
  roof.set_attribute_effect (&scat, CGAL::Classification_type::NEUTRAL_ATT);
  roof.set_attribute_effect (&elev, CGAL::Classification_type::FAVORED_ATT);
  roof.set_attribute_effect (&horiz, CGAL::Classification_type::NEUTRAL_ATT);
  roof.set_attribute_effect (&plan, CGAL::Classification_type::NEUTRAL_ATT);

  // Add types to PSC
  psc.add_classification_type (&vege);
  psc.add_classification_type (&ground);
  psc.add_classification_type (&roof);

  psc.classify (1); // Run with method=1 (global regularization with graphcut)

  // Recover output
  std::vector<Point> pts_ground, pts_vege, pts_roof;
  for (std::size_t i = 0; i < pts.size(); ++ i)
    {
      Classification_type* type = psc.classification_type_of (i);
      switch (type)
        {
        case &ground:
          pts_ground.push_back (i);
          break;
        case &vege:
          pts_vege.push_back (i);
          break;
        case &roof:
          pts_roof.push_back (i);
          break;
        default:
          std::cerr << "Error: unknown classification type" << std::endl;
        }
    }
  
  std::cerr << "All done" << std::endl;
  return EXIT_SUCCESS;
}
