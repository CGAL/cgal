#include <CGAL/Generalized_map.h>
#include <CGAL/Cell_attribute.h>
#include <iostream>
#include <cstdlib>

// My item class: no static functor is associated with Face_attribute.
struct Myitem
{
  template<class GMap>
  struct Dart_wrapper
  {
    typedef CGAL::Cell_attribute<GMap, double> Face_attribute; // A weight
    typedef std::tuple<void,void,Face_attribute> Attributes;
  };
};

// Definition of my generalized map.
typedef CGAL::Generalized_map<3,Myitem> GMap_3;
typedef GMap_3::Dart_handle             Dart_handle;
typedef GMap_3::Attribute_type<2>::type Face_attribute;

// Functor called when two faces are merged.
struct Merge_functor
{
  // operator() automatically called before a merge.
  void operator()(Face_attribute& ca1, Face_attribute& ca2)
  {
    ca1.info()=ca1.info()+ca2.info(); // Update can be done incrementally.
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
     // We need to reinitalize the weight of the two faces
    GMap_3::size_type nb1=mmap.darts_of_cell<2>(ca1.dart()).size();
    GMap_3::size_type nb2=mmap.darts_of_cell<2>(ca2.dart()).size();
    mmap.info<2>(ca1.dart())*=(double(nb1)/(nb1+nb2));
    mmap.info<2>(ca2.dart())*=(double(nb2)/(nb1+nb2));
    std::cout<<"After on split faces: info of face1="<<ca1.info()
             <<", info of face2="<<ca2.info()<<std::endl;
  }

private:
  GMap_3& mmap;
};

// Function allowing to display all the 2-attributes, and the characteristics
// of a given combinatorial map.
void display_map_and_2attributes(GMap_3& gm)
{
  for (GMap_3::Attribute_range<2>::type::iterator
       it=gm.attributes<2>().begin(), itend=gm.attributes<2>().end();
       it!=itend; ++it)
  { std::cout<<gm.info_of_attribute<2>(it)<<"; "; }
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

  // 1) Create and initialize 2-attributes.
  for (GMap_3::One_dart_per_cell_range<2>::iterator
       it=gm.one_dart_per_cell<2>().begin(),
       itend=gm.one_dart_per_cell<2>().end(); it!=itend; ++it)
  {
    gm.set_attribute<2>(it, gm.create_attribute<2>(1));
  }

  // 2) Set the onsplit and onmerge functors.
  gm.onsplit_functor<2>()=Split_functor(gm);
  gm.onmerge_functor<2>()=Merge_functor();

  // 3) 3-Sew the two hexahedra along one face. This calls 1 onmerge.
  gm.sew<3>(dh1, dh2);

  // 4) Display all the values of 2-attributes.
  display_map_and_2attributes(gm);

  // 5) Insert a vertex in the face between the two hexahedra.
  //    This calls 3 onsplit.
  Dart_handle resdart=gm.insert_cell_0_in_cell_2(dh2);

  // 6) Display all the values of 2-attributes.
  display_map_and_2attributes(gm);

  // 7) "Remove" the dynamic onmerge functor.
  gm.onmerge_functor<2>()=boost::function<void(Face_attribute&,
                                               Face_attribute&)>();

  // 8) Remove one edge: this merges two faces, however no dynamic
  //    functor is called (because it was removed).
  gm.remove_cell<1>(resdart);

  // 9) Display all the values of 2-attributes.
  display_map_and_2attributes(gm);

  return EXIT_SUCCESS;
}
