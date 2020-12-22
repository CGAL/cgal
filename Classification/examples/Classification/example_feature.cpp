#if defined (_MSC_VER) && !defined (_WIN64)
#pragma warning(disable:4244) // boost::number_distance::distance()
                              // converts 64 to 32 bits integers
#endif

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

namespace Classification = CGAL::Classification;

typedef Classification::Sum_of_weighted_features_classifier                      Classifier;

typedef Classification::Point_set_neighborhood<Kernel, Point_range, Pmap>       Neighborhood;
typedef Classification::Local_eigen_analysis                                    Local_eigen_analysis;

typedef Classification::Label_handle                                            Label_handle;
typedef Classification::Feature_handle                                          Feature_handle;
typedef Classification::Label_set                                               Label_set;
typedef Classification::Feature_set                                             Feature_set;

typedef Classification::Feature::Verticality<Kernel>                            Verticality;

///////////////////////////////////////////////////////////////////
//! [Feature]

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

  float value (std::size_t pt_index)
  {
    if (xmin < range[pt_index].x() && range[pt_index].x() < xmax &&
        ymin < range[pt_index].y() && range[pt_index].y() < ymax)
      return 1.f;
    else
      return 0.f;
  }

};

//! [Feature]
///////////////////////////////////////////////////////////////////

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
  Local_eigen_analysis eigen
    = Local_eigen_analysis::create_from_point_set (pts, Pmap(), neighborhood.k_neighbor_query(6));

  Label_set labels;
  Label_handle a = labels.add ("label_A");
  Label_handle b = labels.add ("label_B");

  ///////////////////////////////////////////////////////////////////
  //! [Addition]

  std::cerr << "Computing features" << std::endl;
  Feature_set features;

  // Feature that identifies points whose x coordinate is between -20
  // and 20 and whose y coordinate is between -15 and 15
  Feature_handle my_feature = features.add<My_feature> (pts, -20., 20., -15., 15.);

  //! [Addition]
  ///////////////////////////////////////////////////////////////////

  Feature_handle verticality = features.add<Verticality> (pts, eigen);

  Classifier classifier (labels, features);

  std::cerr << "Setting weights" << std::endl;
  classifier.set_weight(verticality, 0.5);
  classifier.set_weight(my_feature, 0.25);

  std::cerr << "Setting up labels" << std::endl;
  classifier.set_effect (a, verticality, Classifier::FAVORING);
  classifier.set_effect (a, my_feature, Classifier::FAVORING);
  classifier.set_effect (b, verticality, Classifier::PENALIZING);
  classifier.set_effect (b, my_feature, Classifier::PENALIZING);

  std::cerr << "Classifying" << std::endl;
  std::vector<int> label_indices(pts.size(), -1);
  Classification::classify_with_graphcut<CGAL::Parallel_if_available_tag>
    (pts, Pmap(), labels, classifier,
     neighborhood.k_neighbor_query(12),
     0.5, 1, label_indices);

  std::cerr << "All done" << std::endl;
  return EXIT_SUCCESS;
}
