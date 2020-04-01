#include <CGAL/Generalized_map.h>
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
  template<class GMap>
  struct Dart_wrapper
  {
    typedef CGAL::Cell_attribute<GMap, int, CGAL::Tag_true,
     Sum_functor, Divide_by_two_functor> Facet_attribute;
    typedef std::tuple<void,void,Facet_attribute> Attributes;
  };
};

typedef CGAL::Generalized_map<3,Myitem> GMap_3;
typedef GMap_3::Dart_handle             Dart_handle;

int main()
{
  GMap_3 gm;

  // Create 2 hexahedra.
  Dart_handle dh1 = gm.make_combinatorial_hexahedron();
  Dart_handle dh2 = gm.make_combinatorial_hexahedron();

  // 1) Create all 2-attributes and associated them to darts.
  for (GMap_3::Dart_range::iterator
         it=gm.darts().begin(), itend=gm.darts().end();
       it!=itend; ++it)
  {
    if ( gm.attribute<2>(it)==NULL )
      gm.set_attribute<2>(it, gm.create_attribute<2>());
  }

  // 2) Set the color of all facets of the first hexahedron to 7.
  for (GMap_3::One_dart_per_incident_cell_range<2, 3>::iterator
         it=gm.one_dart_per_incident_cell<2,3>(dh1).begin(),
         itend=gm.one_dart_per_incident_cell<2,3>(dh1).end(); it!=itend; ++it)
  { gm.info<2>(it)=7; }

  // 3) Set the color of all facets of the second hexahedron to 13.
  for (GMap_3::One_dart_per_incident_cell_range<2, 3>::iterator it=
         gm.one_dart_per_incident_cell<2,3>(dh2).begin(),
         itend=gm.one_dart_per_incident_cell<2,3>(dh2).end(); it!=itend; ++it)
  { gm.info<2>(it)=13; }

  // 4) 3-Sew the two hexahedra along one facet.
  gm.sew<3>(dh1, dh2);

  // 5) Display all the values of 2-attributes.
  for (GMap_3::Attribute_range<2>::type::iterator
       it=gm.attributes<2>().begin(), itend=gm.attributes<2>().end();
       it!=itend; ++it)
  {
    std::cout<<gm.info_of_attribute<2>(it)<<"; ";
  }
  std::cout<<std::endl;

  // 6) Insert a vertex in the facet between the two hexahedra.
  gm.insert_cell_0_in_cell_2(dh2);

  // 7) Display all the values of 2-attributes.
  for (GMap_3::Attribute_range<2>::type::iterator
         it=gm.attributes<2>().begin(), itend=gm.attributes<2>().end();
       it!=itend; ++it)
  {
    std::cout<<gm.info_of_attribute<2>(it)<<"; ";
  }
  std::cout<<std::endl;
  gm.display_characteristics(std::cout);
  std::cout<<", valid="<<gm.is_valid()<<std::endl;

  return EXIT_SUCCESS;
}
