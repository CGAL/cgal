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

#include <CGAL/IO/read_ply_points.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;
typedef Kernel::Iso_cuboid_3 Iso_cuboid_3;
typedef std::vector<Point> Point_range;
typedef CGAL::Identity_property_map<Point> Pmap;

typedef CGAL::Classifier<Point_range, Pmap> Classifier;

typedef CGAL::Classification::Point_set_neighborhood<Kernel, Point_range, Pmap>       Neighborhood;
typedef CGAL::Classification::Local_eigen_analysis<Kernel, Point_range, Pmap>         Local_eigen_analysis;

typedef CGAL::Classification::Label_handle                                            Label_handle;
typedef CGAL::Classification::Feature_handle                                          Feature_handle;

typedef CGAL::Classification::Feature::Sphericity<Kernel, Point_range, Pmap>          Sphericity;


// User-defined feature that identifies a specific area of the 3D
// space. This feature takes value 1 for points that lie inside the
// area and 0 for the others.
class My_feature : public CGAL::Classification::Feature_base
{
  const Point_range& range;
  double xmin, xmax, ymin, ymax;
public:
  My_feature (const Point_range& range, // constructor should start with item range
              double xmin, double xmax, double ymin, double ymax)
    : range (range), xmin(xmin), xmax(xmax), ymin(ymin), ymax(ymax)
  { }

  double value (std::size_t pt_index)
  {
    if (xmin < range[pt_index].x() && range[pt_index].x() < xmax &&
        ymin < range[pt_index].y() && range[pt_index].y() < ymax)
      return 1.;
    else
      return 0.;
  }
    
};

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

  Neighborhood neighborhood (pts, Pmap());
  Local_eigen_analysis eigen (pts, Pmap(), neighborhood.k_neighbor_query(6));

  Classifier classifier (pts, Pmap());

  std::cerr << "Computing features" << std::endl;
  Feature_handle sphericity = classifier.add_feature<Sphericity> (eigen);

  // Feature that identifies points whose x coordinate is between -20
  // and 20 and whose y coordinate is between -15 and 15
  Feature_handle my_feature = classifier.add_feature<My_feature> (-20., 20., -15., 15.);
                                                     
  std::cerr << "Setting weights" << std::endl;
  sphericity->set_weight(0.5);
  my_feature->set_weight(0.25);

  std::cerr << "Setting up labels" << std::endl;
  Label_handle a = classifier.add_label ("label_A");
  a->set_feature_effect (sphericity,  CGAL::Classification::Feature::FAVORING);
  a->set_feature_effect (my_feature,  CGAL::Classification::Feature::FAVORING);

  Label_handle b = classifier.add_label ("label_B");
  b->set_feature_effect (sphericity,  CGAL::Classification::Feature::PENALIZING);
  b->set_feature_effect (my_feature,  CGAL::Classification::Feature::PENALIZING);

  classifier.run_with_graphcut (neighborhood.k_neighbor_query(12), 0.2);

  std::cerr << "All done" << std::endl;
  return EXIT_SUCCESS;
}
