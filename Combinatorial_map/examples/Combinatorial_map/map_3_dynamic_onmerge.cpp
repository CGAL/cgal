#include <CGAL/Combinatorial_map.h>
#include <CGAL/Combinatorial_map_constructors.h>
#include <CGAL/Combinatorial_map_operations.h>
#include <CGAL/Cell_attribute.h>
#include <iostream>
#include <algorithm>
#include <cstdlib>

struct Myitem
{
  template<class CMap>
  struct Dart_wrapper
  {
    typedef CGAL::Dart<3, CMap> Dart;
    typedef CGAL::Cell_attribute<CMap, int> Facet_attribute;
    typedef CGAL::cpp11::tuple<void,void,Facet_attribute> Attributes;
  };
};

typedef CGAL::Combinatorial_map<3,Myitem> CMap_3;
typedef CMap_3::Dart_handle               Dart_handle;
typedef CMap_3::Attribute_type<2>::type   Facet_attribute;

struct Merge_functor
{
  void operator()(Facet_attribute& ca1, Facet_attribute& ca2)
  { ca1.info()=ca1.info()+ca2.info(); }
};

struct Split_functor
{
  Split_functor(CMap_3& amap) : mmap(amap)
  {}

  void operator()(Facet_attribute& ca1, Facet_attribute& ca2)
  {
    set_color_of_facet(ca1.dart());
    set_color_of_facet(ca2.dart());
  }

private:

  void set_color_of_facet(CMap_3::Dart_handle dh)
  {
    int nb=0;
    int sum=0;
    for (CMap_3::One_dart_per_incident_cell_range<2, 3>::iterator it=
           mmap.one_dart_per_incident_cell<2,3>(dh).begin(),
           itend=mmap.one_dart_per_incident_cell<2,3>(dh).end();
         it!=itend; ++it, ++nb)
    { sum+=it->attribute<2>()->info(); }

    dh->attribute<2>()->info()=(sum/nb);
  }


  CMap_3& mmap;
};

void display_map_and_2attributs(CMap_3& cm)
{
  for (CMap_3::Attribute_range<2>::type::iterator
       it=cm.attributes<2>().begin(), itend=cm.attributes<2>().end();
       it!=itend; ++it)
  {
    std::cout<<it->info()<<"; ";
  }
  std::cout<<std::endl;
  cm.display_characteristics(std::cout);
  std::cout<<", valid="<<cm.is_valid()<<std::endl;
}

int main()
{
  CMap_3 cm;

  // Create 2 hexahedra.
  Dart_handle dh1 = CGAL::make_combinatorial_hexahedron(cm);
  Dart_handle dh2 = CGAL::make_combinatorial_hexahedron(cm);

  // 1) Create all 2-attributes and associated them to darts.
  for (CMap_3::Dart_range::iterator
       it=cm.darts().begin(), itend=cm.darts().end();
       it!=itend; ++it)
  {
    if ( it->attribute<2>()==NULL )
      cm.set_attribute<2>(it, cm.create_attribute<2>());
  }

  // 2) Set the color of all facets of the first hexahedron to 7.
  for (CMap_3::One_dart_per_incident_cell_range<2, 3>::iterator
       it=cm.one_dart_per_incident_cell<2,3>(dh1).begin(),
       itend=cm.one_dart_per_incident_cell<2,3>(dh1).end(); it!=itend; ++it)
  { it->attribute<2>()->info()=7; }

  // 3) Set the color of all facets of the second hexahedron to 13.
  for (CMap_3::One_dart_per_incident_cell_range<2, 3>::iterator it=
       cm.one_dart_per_incident_cell<2,3>(dh2).begin(),
       itend=cm.one_dart_per_incident_cell<2,3>(dh2).end(); it!=itend; ++it)
  { it->attribute<2>()->info()=13; }

  cm.onsplit_functor<2>()=Split_functor(cm);
  cm.onmerge_functor<2>()=Merge_functor();

  // 4) 3-Sew the two hexahedra along one facet.
  cm.sew<3>(dh1, dh2);

  // 5) Display all the values of 2-attributes.
  display_map_and_2attributs(cm);

  // 6) Insert a vertex in the facet between the two hexahedra.
  Dart_handle resdart=CGAL::insert_cell_0_in_cell_2(cm, dh2);

  // 7) Display all the values of 2-attributes.
  display_map_and_2attributs(cm);

  // Now there is no dynamic functor called when two 2-attributes are merged
  cm.onmerge_functor<2>()=boost::function<void(Facet_attribute&,
                                               Facet_attribute&)>();

  CGAL::remove_cell<CMap_3, 1>(cm, resdart);
  display_map_and_2attributs(cm);

  return EXIT_SUCCESS;
}
