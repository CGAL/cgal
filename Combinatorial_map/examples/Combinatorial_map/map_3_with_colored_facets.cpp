#include <CGAL/Combinatorial_map.h>
#include <CGAL/Cell_attribute.h>
#include <iostream>
#include <algorithm>
#include <cstdlib>

struct Sum_functor
{
  template<class Cell_attribute>
  void operator()(Cell_attribute& ca1,Cell_attribute& ca2)
  { ca1.info()=ca1.info()+ca2.info(); }
};
struct Divide_by_two_functor
{
  template<class Cell_attribute>
  void operator()(Cell_attribute& ca1,Cell_attribute& ca2)
  {
    ca1.info()=(ca1.info()/2);
    ca2.info()=(ca1.info());
  }
};

struct Myitem
{
  template<class CMap>
  struct Dart_wrapper
  {
    typedef CGAL::Cell_attribute<CMap, int, CGAL::Tag_true,
     Sum_functor, Divide_by_two_functor> Facet_attribute;
    typedef std::tuple<void,void,Facet_attribute> Attributes;
  };
};

typedef CGAL::Combinatorial_map<3,Myitem> CMap_3;
typedef CMap_3::Dart_handle               Dart_handle;

int main()
{
  CMap_3 cm;

  // Create 2 hexahedra.
  Dart_handle dh1 = cm.make_combinatorial_hexahedron();
  Dart_handle dh2 = cm.make_combinatorial_hexahedron();

  // 1) Create all 2-attributes and associated them to darts.
  for (CMap_3::Dart_range::iterator
       it=cm.darts().begin(), itend=cm.darts().end();
       it!=itend; ++it)
  {
    if ( cm.attribute<2>(it)==NULL )
      cm.set_attribute<2>(it, cm.create_attribute<2>());
  }

  // 2) Set the color of all facets of the first hexahedron to 7.
  for (CMap_3::One_dart_per_incident_cell_range<2, 3>::iterator
       it=cm.one_dart_per_incident_cell<2,3>(dh1).begin(),
       itend=cm.one_dart_per_incident_cell<2,3>(dh1).end(); it!=itend; ++it)
  { cm.info<2>(it)=7; }

  // 3) Set the color of all facets of the second hexahedron to 13.
  for (CMap_3::One_dart_per_incident_cell_range<2, 3>::iterator it=
       cm.one_dart_per_incident_cell<2,3>(dh2).begin(),
       itend=cm.one_dart_per_incident_cell<2,3>(dh2).end(); it!=itend; ++it)
  { cm.info<2>(it)=13; }

  // 4) 3-Sew the two hexahedra along one facet.
  cm.sew<3>(dh1, dh2);

  // 5) Display all the values of 2-attributes.
  for (CMap_3::Attribute_range<2>::type::iterator
       it=cm.attributes<2>().begin(), itend=cm.attributes<2>().end();
       it!=itend; ++it)
  {
    std::cout<<cm.info_of_attribute<2>(it)<<"; ";
  }
  std::cout<<std::endl;

  // 6) Insert a vertex in the facet between the two hexahedra.
  cm.insert_cell_0_in_cell_2(dh2);

  // 7) Display all the values of 2-attributes.
  for (CMap_3::Attribute_range<2>::type::iterator
       it=cm.attributes<2>().begin(), itend=cm.attributes<2>().end();
       it!=itend; ++it)
  {
    std::cout<<cm.info_of_attribute<2>(it)<<"; ";
  }
  std::cout<<std::endl;
  cm.display_characteristics(std::cout);
  std::cout<<", valid="<<cm.is_valid()<<std::endl;

  return EXIT_SUCCESS;
}
