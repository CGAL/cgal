#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Classification.h>
#include <CGAL/IO/read_ply_points.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;
typedef Kernel::Iso_cuboid_3 Iso_cuboid_3;
typedef std::vector<Point> Point_range;
typedef CGAL::Identity_property_map<Point> Pmap;

namespace Classif = CGAL::Classification;

typedef Classif::Sum_of_weighted_features_predicate                      Classification_predicate;

typedef Classif::Point_set_neighborhood<Kernel, Point_range, Pmap>       Neighborhood;
typedef Classif::Local_eigen_analysis<Kernel, Point_range, Pmap>         Local_eigen_analysis;

typedef Classif::Label_handle                                            Label_handle;
typedef Classif::Feature_handle                                          Feature_handle;
typedef Classif::Label_set                                               Label_set;
typedef Classif::Feature_set                                             Feature_set;

typedef Classif::Feature::Sphericity<Kernel, Point_range, Pmap>          Sphericity;


// User-defined feature that identifies a specific area of the 3D
// space. This feature takes value 1 for points that lie inside the
// area and 0 for the others.
class My_feature : public CGAL::Classification::Feature_base
{
  const Point_range& range;
  double xmin, xmax, ymin, ymax;
public:
  My_feature (const Point_range& range,
              double xmin, double xmax, double ymin, double ymax)
    : range (range), xmin(xmin), xmax(xmax), ymin(ymin), ymax(ymax)
  {
    this->set_name ("my_feature");
  }

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

  Label_set labels;
  Label_handle a = labels.add ("label_A");
  Label_handle b = labels.add ("label_B");
  
  std::cerr << "Computing features" << std::endl;
  Feature_set features;
  
  Feature_handle sphericity = features.add<Sphericity> (pts, eigen);

  // Feature that identifies points whose x coordinate is between -20
  // and 20 and whose y coordinate is between -15 and 15
  Feature_handle my_feature = features.add<My_feature> (pts, -20., 20., -15., 15.);

  Classification_predicate predicate (labels, features);
  
  std::cerr << "Setting weights" << std::endl;
  predicate.set_weight(sphericity, 0.5);
  predicate.set_weight(my_feature, 0.25);

  std::cerr << "Setting up labels" << std::endl;
  predicate.set_effect (a, sphericity, Classification_predicate::FAVORING);
  predicate.set_effect (a, my_feature, Classification_predicate::FAVORING);
  predicate.set_effect (b, sphericity, Classification_predicate::PENALIZING);
  predicate.set_effect (b, my_feature, Classification_predicate::PENALIZING);

  std::vector<std::size_t> label_indices;
  Classif::classify_with_graphcut<CGAL::Sequential_tag>
    (pts, Pmap(), Pmap(), labels, predicate,
     neighborhood.k_neighbor_query(12),
     0.5, 1, label_indices);

  std::cerr << "All done" << std::endl;
  return EXIT_SUCCESS;
}
