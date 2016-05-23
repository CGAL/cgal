#include <CGAL/Generalized_map.h>
#include <CGAL/Cell_attribute.h>
#include <iostream>
#include <cstdlib>

// My item class: no functor is associated with Face_attribute.
struct Myitem
{
  template<class GMap>
  struct Dart_wrapper
  {
    typedef CGAL::GMap_dart<3, GMap> Dart;
    typedef CGAL::Cell_attribute<GMap, int> Face_attribute;
    typedef CGAL::cpp11::tuple<void,void,Face_attribute> Attributes;
  };
};

// Definition of my combinatorial map.
typedef CGAL::Generalized_map<3,Myitem> GMap_3;
typedef GMap_3::Dart_handle             Dart_handle;
typedef GMap_3::Attribute_type<2>::type Face_attribute;

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
  Split_functor(GMap_3& amap) : mmap(amap)
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
  void set_color_of_face(GMap_3::Dart_handle dh)
  {
    int nb=0;
    int sum=0;
    for (GMap_3::Dart_of_cell_range<1>::iterator
           it=mmap.darts_of_cell<1>(dh).begin(),
           itend=mmap.darts_of_cell<1>(dh).end();
         it!=itend; ++it, ++nb)
    { sum+=mmap.info<2>(mmap.alpha<2>(it)); }

    mmap.info<2>(dh)=(sum/nb);
  }

  GMap_3& mmap;
};

// Function allowing to display all the 2-attributes, and the characteristics
// of a given combinatorial map.
void display_map_and_2attributes(GMap_3& gm)
{
  for (GMap_3::Attribute_range<2>::type::iterator
       it=gm.attributes<2>().begin(), itend=gm.attributes<2>().end();
       it!=itend; ++it)
  {
    std::cout<<gm.info_of_attribute<2>(it)<<"; ";
  }
  std::cout<<std::endl;
  gm.display_characteristics(std::cout);
  std::cout<<", valid="<<gm.is_valid()<<std::endl;
}

int main()
{
  GMap_3 gm;

  // 0) Create 2 hexahedra.
  Dart_handle dh1 = gm.make_combinatorial_hexahedron();
  Dart_handle dh2 = gm.make_combinatorial_hexahedron();

  // 1) Create 2-attributes of the first hexahedron, info()==7.
  for (GMap_3::One_dart_per_incident_cell_range<2, 3>::iterator
       it=gm.one_dart_per_incident_cell<2,3>(dh1).begin(),
       itend=gm.one_dart_per_incident_cell<2,3>(dh1).end(); it!=itend; ++it)
  { gm.set_attribute<2>(it, gm.create_attribute<2>(7)); }

  // 2) Create 2-attributes of the second hexahedron, info()==13.
  for (GMap_3::One_dart_per_incident_cell_range<2, 3>::iterator it=
       gm.one_dart_per_incident_cell<2,3>(dh2).begin(),
       itend=gm.one_dart_per_incident_cell<2,3>(dh2).end(); it!=itend; ++it)
  { gm.set_attribute<2>(it, gm.create_attribute<2>(13)); }

  // 3) Set the onsplit and onmerge functors.
  gm.onsplit_functor<2>()=Split_functor(gm);
  gm.onmerge_functor<2>()=Merge_functor();

  // 4) 3-Sew the two hexahedra along one face. This calls 1 onmerge.
  gm.sew<3>(dh1, dh2);

  // 5) Display all the values of 2-attributes.
  display_map_and_2attributes(gm);

  // 6) Insert a vertex in the face between the two hexahedra.
  //    This calls 4 onsplit.
  // TODO Dart_handle resdart=gm.insert_cell_0_in_cell_2(dh2);

  // 7) Display all the values of 2-attributes.
  display_map_and_2attributes(gm);

  // 8) "Remove" the dynamic onmerge functor.
  gm.onmerge_functor<2>()=boost::function<void(Face_attribute&,
                                               Face_attribute&)>();

  // 9) Remove one edge: this merges two faces, however no dynamic
  //    functor is called (because it was removed).
  // TODO  gm.remove_cell<1>(resdart);

  // 10) Display all the values of 2-attributes.
  display_map_and_2attributes(gm);

  return EXIT_SUCCESS;
}
