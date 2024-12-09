#include <CGAL/Generalized_map.h>
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
  template<class GMap>
  struct Dart_wrapper
  {
    typedef CGAL::Cell_attribute<GMap, MyInfo> attrib;
    typedef std::tuple<void, void, attrib> Attributes;
  };
};

struct Myitem2
{
  template<class GMap>
  struct Dart_wrapper
  {
    typedef CGAL::Cell_attribute<GMap, MyInfo> attrib;
    typedef std::tuple<void, void, attrib> Attributes;
  };
};

using GMap1=CGAL::Generalized_map<3,Myitem1>;
using GMap2=CGAL::Generalized_map<3,Myitem2>;

#define NB 1000
template<typename GMap>
bool test(const std::string& s)
{
  bool res=true;
  GMap m;
  // 1) create a face and one attribute.
  typename GMap::Dart_descriptor dd=m.make_combinatorial_polygon(4);
  m.template set_attribute<2>(dd, m.template create_attribute<2>(2));
  // 2) Split this face NB times => will create new 2-attributes for new faces
  for(std::size_t i=0; i<NB; ++i)
  {
    typename GMap::Dart_descriptor
      newd=m.insert_cell_1_in_cell_2(dd, m.template alpha<0,1>(dd));
    if(m.template attribute<2>(newd)==GMap::null_descriptor)
    {
      std::cout<<"ERROR1: "<<s<<": "
               <<"attribute<2>(newd)==GMap::null_descriptor"<<std::endl;
      res=false;
    }
    else if(m.template info<2>(newd).data!=2)
    {
      std::cout<<"ERROR2: "<<s<<": "<<m.template info<2>(newd).data<<std::endl;
      res=false;
    }

    newd=m.template opposite<2>(newd);
    if(m.template attribute<2>(newd)==GMap::null_descriptor)
    {
      std::cout<<"ERROR3: "<<s<<": "
               <<"attribute<2>(newd)==GMap::null_descriptor"<<std::endl;
      res=false;
    }
    else if(m.template info<2>(newd).data!=2)
    {
      std::cout<<"ERROR4: "<<s<<": "<<m.template info<2>(newd).data<<std::endl;
      res=false;
    }
  }
  return res;
}

int main()
{
  if(!test<GMap1>("GMap1") || !test<GMap2>("GMap2"))
  { return EXIT_FAILURE; }

  return EXIT_SUCCESS;
}
