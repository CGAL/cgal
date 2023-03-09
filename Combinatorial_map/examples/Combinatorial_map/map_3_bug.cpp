#include <CGAL/Combinatorial_map.h>
#include <CGAL/Cell_attribute.h>
#include <vector>
#include <algorithm>
#include <string>
#include <cstdlib>

struct MyInfo
{
  MyInfo() :data(1)
  {}

  MyInfo(int i) :data(i)
  {}

  int data;
};

struct Myitem1
{
  using Use_index=CGAL::Tag_true; // use indices
  using Index_type=std::uint16_t; // 16 bits
  template<class CMap>
  struct Dart_wrapper
  {
    typedef CGAL::Cell_attribute<CMap, MyInfo> attrib;
    typedef std::tuple<attrib> Attributes;
  };
};

struct Myitem2
{
  template<class CMap>
  struct Dart_wrapper
  {
    typedef CGAL::Cell_attribute<CMap, MyInfo> attrib;
    typedef std::tuple<attrib> Attributes;
  };
};

using CMap1=CGAL::Combinatorial_map<3,Myitem1>;
using CMap2=CGAL::Combinatorial_map<3,Myitem2>;

#define NB 100
int main()
{
  CMap1 cm1;
  CMap2 cm2;

  CMap1::Attribute_descriptor<0>::type a1=cm1.create_attribute<0>(2);
  for(std::size_t i=0; i<NB; ++i)
  {
    CMap1::Attribute_descriptor<0>::type a2=cm1.create_attribute<0>(cm1.info_of_attribute<0>(a1));
    if(cm1.info_of_attribute<0>(a1).data!=cm1.info_of_attribute<0>(a2).data)
    { std::cout<<"ERROR1: "<<cm1.info_of_attribute<0>(a1).data<<" and "<<cm1.info_of_attribute<0>(a2).data<<" ("<<i<<")"<<std::endl; }
  }
  
  CMap2::Attribute_descriptor<0>::type a3=cm2.create_attribute<0>(3);
  for(std::size_t i=0; i<NB; ++i)
  {
    CMap2::Attribute_descriptor<0>::type a4=cm2.create_attribute<0>(cm2.info_of_attribute<0>(a3));
    if(cm2.info_of_attribute<0>(a3).data!=cm2.info_of_attribute<0>(a4).data)
    { std::cout<<"ERROR2: "<<cm2.info_of_attribute<0>(a3).data<<" and "<<cm2.info_of_attribute<0>(a4).data<<" ("<<i<<")"<<std::endl; }
  }

  return EXIT_SUCCESS;
}
