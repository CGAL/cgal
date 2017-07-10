#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Random.h>
#include <CGAL/tuple.h>
#include <CGAL/Classifier.h>
#include <CGAL/Classification/Local_eigen_analysis.h>
#include <CGAL/Classification/Point_set_neighborhood.h>
#include <CGAL/Classification/Planimetric_grid.h>
#include <CGAL/Classification/Attribute_base.h>
#include <CGAL/Classification/Attribute/Distance_to_plane.h>
#include <CGAL/Classification/Attribute/Echo_scatter.h>
#include <CGAL/Classification/Attribute/Eigen.h>
#include <CGAL/Classification/Attribute/Elevation.h>
#include <CGAL/Classification/Attribute/Hsv.h>
#include <CGAL/Classification/Attribute/Vertical_dispersion.h>
#include <CGAL/Classification/Attribute/Verticality.h>

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

typedef CGAL::Classifier<Point_range, Pmap> Classifier;

typedef CGAL::Classification::Planimetric_grid<Kernel, Point_range, Pmap>       Planimetric_grid;
typedef CGAL::Classification::Point_set_neighborhood<Kernel, Point_range, Pmap> Neighborhood;
typedef CGAL::Classification::Local_eigen_analysis<Kernel, Point_range, Pmap>   Local_eigen_analysis;

typedef CGAL::Classification::Type_handle                                           Type_handle;
typedef CGAL::Classification::Attribute_handle                                      Attribute_handle;

typedef CGAL::Classification::Attribute::Anisotropy
<Kernel, Point_range, Pmap> Anisotropy;
typedef CGAL::Classification::Attribute::Distance_to_plane
<Kernel, Point_range, Pmap> Distance_to_plane;
typedef CGAL::Classification::Attribute::Eigentropy
<Kernel, Point_range, Pmap> Eigentropy;
typedef CGAL::Classification::Attribute::Elevation
<Kernel, Point_range, Pmap>                    Elevation;
typedef CGAL::Classification::Attribute::Linearity
<Kernel, Point_range, Pmap> Linearity;
typedef CGAL::Classification::Attribute::Omnivariance
<Kernel, Point_range, Pmap> Omnivariance;
typedef CGAL::Classification::Attribute::Planarity
<Kernel, Point_range, Pmap> Planarity;
typedef CGAL::Classification::Attribute::Sphericity
<Kernel, Point_range, Pmap> Sphericity;
typedef CGAL::Classification::Attribute::Sum_eigenvalues
<Kernel, Point_range, Pmap> Sum_eigen;
typedef CGAL::Classification::Attribute::Surface_variation
<Kernel, Point_range, Pmap> Surface_variation;
typedef CGAL::Classification::Attribute::Vertical_dispersion
<Kernel, Point_range, Pmap>                    Dispersion;
typedef CGAL::Classification::Attribute::Verticality
<Kernel, Point_range, Pmap> Verticality;
typedef CGAL::Classification::RGB_Color RGB_Color;
typedef CGAL::Classification::Attribute::Hsv
<Kernel, Point_range, Cmap> Hsv;
typedef CGAL::Classification::Attribute::Hsv
<Kernel, Point_range, Cmap> Hsv;
typedef CGAL::Classification::Attribute::Echo_scatter
<Kernel, Point_range, Pmap, Emap> Echo_scatter;

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

  double grid_resolution = 0.34;
  double radius_neighbors = 1.7;
  double radius_dtm = 15.0;
  
  Classifier classifier (range, Pmap());
  Iso_cuboid_3 bbox = CGAL::bounding_box (boost::make_transform_iterator
                                          (range.begin(),
                                           CGAL::Property_map_to_unary_function<Pmap>()),
                                          boost::make_transform_iterator
                                          (range.end(),
                                           CGAL::Property_map_to_unary_function<Pmap>()));

  Planimetric_grid grid (range, Pmap(), bbox, grid_resolution);
  Neighborhood neighborhood (range, Pmap());
  Local_eigen_analysis eigen (range, Pmap(), neighborhood.k_neighbor_query(6));

  std::vector<Attribute_handle> att;
  
  att.push_back (classifier.add_attribute<Anisotropy> (eigen));
  att.push_back (classifier.add_attribute<Dispersion> (Pmap(), grid,
                                                       grid_resolution,
                                                       radius_neighbors));
  att.push_back (classifier.add_attribute<Distance_to_plane> (Pmap(), eigen));
  att.push_back (classifier.add_attribute<Echo_scatter> (Emap(), grid, grid_resolution, radius_neighbors));
  att.push_back (classifier.add_attribute<Eigentropy> (eigen));
  att.push_back (classifier.add_attribute<Elevation> (Pmap(), grid,
                                                      grid_resolution,
                                                      radius_dtm));
  att.push_back (classifier.add_attribute<Hsv> (Cmap(), 0, 180, 90));
  att.push_back (classifier.add_attribute<Hsv> (Cmap(), 1, 50, 25));
  att.push_back (classifier.add_attribute<Hsv> (Cmap(), 2, 50, 25));
  att.push_back (classifier.add_attribute<Linearity> (eigen));
  att.push_back (classifier.add_attribute<Omnivariance> (eigen));
  att.push_back (classifier.add_attribute<Planarity> (eigen));
  att.push_back (classifier.add_attribute<Sphericity> (eigen));
  att.push_back (classifier.add_attribute<Sum_eigen> (eigen));
  att.push_back (classifier.add_attribute<Surface_variation> (eigen));
  att.push_back (classifier.add_attribute<Verticality> (eigen));
  att.push_back (classifier.add_attribute<Verticality> (Nmap()));
  att.push_back (classifier.add_attribute<Omnivariance> (eigen));
  assert (classifier.remove_attribute (att.back()));
  att.pop_back();
  assert (att.size() == classifier.number_of_attributes());
  
  std::vector<Type_handle> types;
  types.push_back (classifier.add_classification_type ("type1"));
  types.push_back (classifier.add_classification_type ("type2"));
  types.push_back (classifier.add_classification_type ("type3"));
  types.push_back (classifier.add_classification_type ("type4"));
  types.push_back (classifier.add_classification_type ("type5"));
  types.push_back (classifier.add_classification_type ("type6"));
  assert (classifier.remove_classification_type (types.back()));
  types.pop_back();
  assert (types.size() == classifier.number_of_classification_types());

  for (std::size_t i = 0; i < types.size(); ++ i)
    for (std::size_t j = 0; j < att.size(); ++ j)
      {
        att[j]->set_weight (CGAL::get_default_random().get_double());
        types[i]->set_attribute_effect
          (att[j], (Attribute_effect)(CGAL::get_default_random().get_int(0,3)));
      }

  classifier.run();
  classifier.run_with_local_smoothing(neighborhood.k_neighbor_query(6));
  classifier.run_with_graphcut(neighborhood.k_neighbor_query(12), 0.5);

  std::vector<Type_handle> output;
  output.reserve (10000);
  for (std::size_t i = 0; i < range.size(); ++ i)
    {
      output.push_back (classifier.classification_type_of(i));
      assert (0 <= classifier.confidence_of(i));
    }
  

  classifier.clear_classification_types();
  classifier.clear_attributes();
  
  return EXIT_SUCCESS;
}
