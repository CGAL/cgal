#include <CGAL/Generalized_map.h>
#include <CGAL/Cell_attribute.h>
#include <vector>
#include <algorithm>
#include <string>
#include <cstdlib>

struct Myitem
{
  using Use_index=CGAL::Tag_true; // use indices
  using Index_type=std::uint16_t; // 16 bits
  template<class GMap>
  struct Dart_wrapper
  {
    typedef CGAL::Cell_attribute<GMap> Vol_attrib;
    typedef std::tuple<void,void,void,Vol_attrib> Attributes;
  };
};

using GMap_3=CGAL::Generalized_map<3,Myitem>;
using Dart_descriptor=GMap_3::Dart_descriptor;

int main()
{
  GMap_3 gm;

  // Create 2 hexahedra.
  Dart_descriptor d1=gm.make_combinatorial_hexahedron();
  Dart_descriptor d2=gm.make_combinatorial_hexahedron();
  gm.sew<3>(d1, d2); // 3-Sew the two hexahedra along one facet.

  // Create two 3-attributes and associated them to darts.
  gm.set_attribute<3>(d1, gm.create_attribute<3>());
  gm.set_attribute<3>(d2, gm.create_attribute<3>());

  // Associate a vector of size_t to darts
  std::vector<std::size_t> array_for_darts(gm.upper_bound_on_dart_ids());
  std::generate(array_for_darts.begin(), array_for_darts.end(), std::rand);

  // Associate a vector of string to 3-attributes
  std::vector<std::string> array_for_vols(gm.upper_bound_on_attribute_ids<3>());
  std::size_t i=0;
  for(auto& e: array_for_vols)
  { e="vol"+std::to_string(i++); }

  std::cout<<"Value in array for darts d1 and d2: "
           <<array_for_darts[d1]<<" and "<<array_for_darts[d2]<<std::endl;

  std::cout<<"Value in array for volumes of dart d1 and d2: "
           <<array_for_vols[gm.attribute<3>(d1)]<<" and "
           <<array_for_vols[gm.attribute<3>(d2)]<<std::endl;

  return EXIT_SUCCESS;
}
