#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Random.h>
#include <CGAL/tuple.h>
#include <CGAL/Point_set_classifier.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
typedef CGAL::cpp11::array<unsigned char, 3> Color;
typedef Kernel::Iso_cuboid_3 Iso_cuboid_3;

typedef boost::tuple<Point, Vector, Color, std::size_t> Point_with_info;
typedef std::vector<Point_with_info> Point_range;

typedef CGAL::Nth_of_tuple_property_map<0, Point_with_info> Pmap;
typedef CGAL::Nth_of_tuple_property_map<1, Point_with_info> Nmap;
typedef CGAL::Nth_of_tuple_property_map<2, Point_with_info> Cmap;
typedef CGAL::Nth_of_tuple_property_map<3, Point_with_info> Emap;

typedef CGAL::Point_set_classifier<Kernel, Point_range, Pmap> PSC;

typedef CGAL::Classification::Type_handle                                           Type_handle;
typedef CGAL::Classification::Attribute_handle                                      Attribute_handle;

typedef CGAL::Classification::Attribute::Effect Attribute_effect;

int main (int, char**)
{
  Point_range range;
  range.reserve(10000);
  for (std::size_t i = 0; i < 10000; ++ i)
    range.push_back (boost::make_tuple
                     (Point (CGAL::get_default_random().get_double(),
                             CGAL::get_default_random().get_double(),
                             CGAL::get_default_random().get_double()),
                      Vector (CGAL::get_default_random().get_double(),
                              CGAL::get_default_random().get_double(),
                              CGAL::get_default_random().get_double()),
                      CGAL::make_array ((unsigned char)(CGAL::get_default_random().get_int(0, 255)),
                                        (unsigned char)(CGAL::get_default_random().get_int(0, 255)),
                                        (unsigned char)(CGAL::get_default_random().get_int(0, 255))),
                      std::size_t(CGAL::get_default_random().get_int(0, 10))));

  PSC classifier (range, Pmap());
  classifier.generate_attributes (5, Nmap(), Cmap(), Emap());
  assert (classifier.number_of_scales() == 5);

  std::vector<Type_handle> types;
  types.push_back (classifier.add_classification_type ("type1"));
  types.push_back (classifier.add_classification_type ("type2"));
  types.push_back (classifier.add_classification_type ("type3"));
  types.push_back (classifier.add_classification_type ("type4"));
  types.push_back (classifier.add_classification_type ("type5"));

  for (std::size_t i = 0; i < 100; ++ i)
    classifier.set_inlier (types[CGAL::get_default_random().get_int(0, types.size())], i);

  classifier.train(500);


  std::ofstream fconfig ("config.xml");
  classifier.save_configuration (fconfig);
  fconfig.close();
  
  std::ofstream fclassif ("classif.ply");
  classifier.write_classification_to_ply(fclassif);
  fclassif.close();
  
  classifier.clear();

  std::ifstream fconfig2 ("config.xml");
  assert (classifier.load_configuration (fconfig2, Nmap(), Cmap(), Emap()));
  
  
  return EXIT_SUCCESS;
}
