#include <CGAL/Combinatorial_map.h>
#include <CGAL/Combinatorial_map_constructors.h>
#include <CGAL/Combinatorial_map_operations.h>
#include <CGAL/Cell_attribute.h>
#include <iostream>
#include <cstdlib>

// My item class: no functor is associated with Face_attribute.
struct Myitem
{
  template<class CMap>
  struct Dart_wrapper
  {
    typedef CGAL::Dart<3, CMap> Dart;
    typedef CGAL::Cell_attribute<CMap, int> Face_attribute;
    typedef CGAL::cpp11::tuple<void,void,Face_attribute> Attributes;
  };
};

// Definition of my combinatorial map.
typedef CGAL::Combinatorial_map<3,Myitem> CMap_3;
typedef CMap_3::Dart_handle               Dart_handle;
typedef CMap_3::Attribute_type<2>::type   Face_attribute;

// Functor called when two faces are merged.
struct Merge_functor
{
  // operator() automatically called before a merge.
  void operator()(Face_attribute& ca1, Face_attribute& ca2)
  {
    ca1.info()=ca1.info()+ca2.info();
    std::cout<<"After on merge faces: info of face1="<<ca1.info()
             <<", info of face2="<<ca2.info()<<std::endl;
  }
};

// Functor called when one face is split in two.
struct Split_functor
{
  Split_functor(CMap_3& amap) : mmap(amap)
  {}

  // operator() automatically called after a split.
  void operator()(Face_attribute& ca1, Face_attribute& ca2)
  {
    set_color_of_face(ca1.dart());
    set_color_of_face(ca2.dart());
    std::cout<<"After on split faces: info of face1="<<ca1.info()
             <<", info of face2="<<ca2.info()<<std::endl;
  }

private:

  // The info of a face is the mean of the info of all its neighboors faces.
  void set_color_of_face(CMap_3::Dart_handle dh)
  {
    int nb=0;
    int sum=0;
    for (CMap_3::Dart_of_orbit_range<1>::iterator
           it=mmap.darts_of_orbit<1>(dh).begin(),
           itend=mmap.darts_of_orbit<1>(dh).end();
         it!=itend; ++it, ++nb)
    { sum+=mmap.info<2>(mmap.beta<2>(it)); }

    mmap.info<2>(dh)=(sum/nb);
  }

  CMap_3& mmap;
};

// Function allowing to display all the 2-attributes, and the characteristics
// of a given combinatorial map.
void display_map_and_2attributes(CMap_3& cm)
{
  for (CMap_3::Attribute_range<2>::type::iterator
       it=cm.attributes<2>().begin(), itend=cm.attributes<2>().end();
       it!=itend; ++it)
  {
    std::cout<<cm.info_of_attribute<2>(it)<<"; ";
  }
  std::cout<<std::endl;
  cm.display_characteristics(std::cout);
  std::cout<<", valid="<<cm.is_valid()<<std::endl;
}

int main()
{
  CMap_3 cm;

  // 0) Create 2 hexahedra.
  Dart_handle dh1 = CGAL::make_combinatorial_hexahedron(cm);
  Dart_handle dh2 = CGAL::make_combinatorial_hexahedron(cm);

  // 1) Create 2-attributes of the first hexahedron, info()==7.
  for (CMap_3::One_dart_per_incident_cell_range<2, 3>::iterator
       it=cm.one_dart_per_incident_cell<2,3>(dh1).begin(),
       itend=cm.one_dart_per_incident_cell<2,3>(dh1).end(); it!=itend; ++it)
  { cm.set_attribute<2>(it, cm.create_attribute<2>(7)); }

  // 2) Create 2-attributes of the second hexahedron, info()==13.
  for (CMap_3::One_dart_per_incident_cell_range<2, 3>::iterator it=
       cm.one_dart_per_incident_cell<2,3>(dh2).begin(),
       itend=cm.one_dart_per_incident_cell<2,3>(dh2).end(); it!=itend; ++it)
  { cm.set_attribute<2>(it, cm.create_attribute<2>(13)); }

  // 3) Set the onsplit and onmerge functors.
  cm.onsplit_functor<2>()=Split_functor(cm);
  cm.onmerge_functor<2>()=Merge_functor();

  // 4) 3-Sew the two hexahedra along one face. This calls 1 onmerge.
  cm.sew<3>(dh1, dh2);

  // 5) Display all the values of 2-attributes.
  display_map_and_2attributes(cm);

  // 6) Insert a vertex in the face between the two hexahedra.
  //    This calls 4 onsplit.
  Dart_handle resdart=CGAL::insert_cell_0_in_cell_2(cm, dh2);

  // 7) Display all the values of 2-attributes.
  display_map_and_2attributes(cm);

  // 8) "Remove" the dynamic onmerge functor.
  cm.onmerge_functor<2>()=boost::function<void(Face_attribute&,
                                               Face_attribute&)>();

  // 9) Remove one edge: this merges two faces, however no dynamic
  //    functor is called (because it was removed).
  CGAL::remove_cell<CMap_3, 1>(cm, resdart);

  // 10) Display all the values of 2-attributes.
  display_map_and_2attributes(cm);

  return EXIT_SUCCESS;
}
