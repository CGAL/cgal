#include <cstdlib>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <CGAL/Combinatorial_map.h>
#include <CGAL/Combinatorial_map_constructors.h>
#include <CGAL/Combinatorial_map_operations.h>
#include <CGAL/Combinatorial_map_save_load.h>
#include <CGAL/Cell_attribute.h>
#include <CGAL/Linear_cell_complex.h>
#include <CGAL/Linear_cell_complex_constructors.h>
#include <CGAL/Linear_cell_complex_operations.h>
#include <cstdlib>
//******************************************************************************
// Functor used to display all the attributes of a given dimension.
template<typename CMap>
struct My_functor_display_attrib
{
  template <unsigned int i>
  static void run(const CMap& amap)
  {
    std::cout<<i<<"-attributes: ";

    typename CMap::template Attribute_range<i>::type::const_iterator
      it_attrib    = amap.template attributes<i>().begin(),
      itend_attrib = amap.template attributes<i>().end();

    for (; it_attrib!=itend_attrib; ++it_attrib)
    {
      std::cout<<it_attrib->info()<<"  ";
    }
    std::cout<<std::endl;
  }
};

// Display the map, i.e. all its characteristics, then all its attributes.
template < class CMap >
void display_map(const CMap& amap)
{
  amap.display_characteristics(std::cout)
    <<", valid="<<amap.is_valid()<<std::endl;
  
  CMap::Helper::template Foreach_enabled_attributes
    <My_functor_display_attrib<CMap> >::run(amap);
}

//******************************************************************************
// Functor used to display all the attributes of a given dimension.
template<typename CMap, unsigned int i, typename Info=typename CMap::template Attribute_type<i>::type::Info>
struct My_functor_display_one_attrib_lcc
{
  static void run(const CMap& amap)
  {
    std::cout<<i<<"-attributes: ";

    typename CMap::template Attribute_range<i>::type::const_iterator
      it_attrib    = amap.template attributes<i>().begin(),
      itend_attrib = amap.template attributes<i>().end();

    for (; it_attrib!=itend_attrib; ++it_attrib)
    {
      std::cout<<it_attrib->info()<<"  ";
    }
    std::cout<<std::endl;
  }
};
template<typename CMap, unsigned int i>
struct My_functor_display_one_attrib_lcc<CMap, i, void>
{
  static void run(const CMap& amap)
  {
    std::cout<<i<<"-attributes: ";

    std::cout<<amap.template attributes<i>().size()<<" attributes without info."
             <<std::endl;
  }
};
template<typename CMap>
struct My_functor_display_one_attrib_lcc<CMap, 0, void>
{
  static void run(const CMap& amap)
  {
    std::cout<<"0-attributes: ";

    typename CMap::template Attribute_range<0>::type::const_iterator
      it_attrib    = amap.template attributes<0>().begin(),
      itend_attrib = amap.template attributes<0>().end();

    for (; it_attrib!=itend_attrib; ++it_attrib)
    {
      std::cout<<"("<<it_attrib->point()<<")  ";
    }
    std::cout<<std::endl;
  }
};
template<typename CMap, typename Info>
struct My_functor_display_one_attrib_lcc<CMap, 0, Info>
{
  static void run(const CMap& amap)
  {
    std::cout<<"0-attributes: ";

    typename CMap::template Attribute_range<0>::type::const_iterator
      it_attrib    = amap.template attributes<0>().begin(),
      itend_attrib = amap.template attributes<0>().end();

    for (; it_attrib!=itend_attrib; ++it_attrib)
    {
      std::cout<<"("<<it_attrib->point()<<" "<<it_attrib->info()<<")  ";
    }
    std::cout<<std::endl;
  }
};

template<typename CMap>
struct My_functor_display_attrib_lcc
{
  template <unsigned int i>
  static void run(const CMap& amap)
  { My_functor_display_one_attrib_lcc<CMap, i>::run(amap); }
};

// Display the map, i.e. all its characteristics, then all its attributes.
template < class CMap >
void display_lcc(const CMap& amap)
{
  amap.display_characteristics(std::cout)
    <<", valid="<<amap.is_valid()<<std::endl;
  
  CMap::Helper::template Foreach_enabled_attributes
    <My_functor_display_attrib_lcc<CMap> >::run(amap);
}

//******************************************************************************
/* Configuration example 1
   Map containing attributes defined with basic types. */
struct MesAttributs1
{
  template < class CMap>
  struct Dart_wrapper
  {
    typedef CGAL::Dart<3, CMap> Dart;
    typedef CGAL::Cell_attribute<CMap, int> T1;
    typedef CGAL::Cell_attribute<CMap, float> T2;
    typedef CGAL::cpp0x::tuple<void,T1,T2> Attributes;
  };
};

typedef CGAL::Combinatorial_map<3, MesAttributs1> CMap_3a;
typedef CMap_3a::Dart_handle Dart_handle_a;

void example1()
{
  CMap_3a cm;

  // 1) create two tetrahedra 3-sewn
  Dart_handle_a dh1 = CGAL::make_combinatorial_tetrahedron(cm);
  Dart_handle_a dh2 = CGAL::make_combinatorial_tetrahedron(cm);
  cm.sew<3>(dh1,dh2);

  // 2) make attributes
  for (CMap_3a::Dart_range::iterator it(cm.darts().begin()),
         itend(cm.darts().end()); it!=itend; ++it)
  {
    if ( it->attribute<1>()==NULL )
    {
      cm.set_attribute<1>(it, cm.create_attribute<1>());
      it->attribute<1>()->info()=cm.number_of_attributes<1>();
    }
    if( it->attribute<2>()==NULL)
    {
      cm.set_attribute<2>(it, cm.create_attribute<2>());
      it->attribute<2>()->info()=(cm.number_of_attributes<2>()/2.0);
    }
  }
  
  std::cout<<"************** example 1 **************"<<std::endl;
  display_map(cm);
  
  // 3) save
  std::ofstream output("output_example_1.xml");
  output<<cm;
}

//******************************************************************************
/* Configuration example 2
   Map containing custom attributes, without overload write_cmap_attribute_node.
   In this case, custom attributes are not saved. */
struct ACustomType
{
  int anint;
  float afloat;
  ACustomType(int ai=0, float af=0.0) : anint(ai), afloat(af)
  {}
  friend std::ostream& operator<<(std::ostream& os, const ACustomType& a)
  { return os<<"["<<a.anint<<", "<<a.afloat<<"]"; }
};

struct MesAttributs2
{
  template < class CMap>
  struct Dart_wrapper
  {
    typedef CGAL::Dart<3, CMap> Dart;
    typedef CGAL::Cell_attribute<CMap, int> T1;
    typedef CGAL::Cell_attribute<CMap, ACustomType> T2;
    typedef CGAL::cpp0x::tuple<void,T1,T2> Attributes;
  };
};

typedef CGAL::Combinatorial_map<3, MesAttributs2> CMap_3b;
typedef CMap_3b::Dart_handle Dart_handle_b;

void example2()
{
  CMap_3b cm;

  // 1) create two tetrahedra 3-sewn
  Dart_handle_b dh1 = CGAL::make_combinatorial_tetrahedron(cm);
  Dart_handle_b dh2 = CGAL::make_combinatorial_tetrahedron(cm);
  cm.sew<3>(dh1,dh2);

  // 2) make attributes
  for (CMap_3b::Dart_range::iterator it(cm.darts().begin()),
         itend(cm.darts().end()); it!=itend; ++it)
  {
    if ( it->attribute<1>()==NULL )
    {
      cm.set_attribute<1>(it, cm.create_attribute<1>());
      it->attribute<1>()->info()=cm.number_of_attributes<1>();
    }
    if( it->attribute<2>()==NULL)
    {
      cm.set_attribute<2>(it, cm.create_attribute<2>());
      it->attribute<2>()->info()=ACustomType((int)cm.number_of_attributes<2>(),
                                             cm.number_of_attributes<2>()/2.0);
    }
  }

  std::cout<<"************** example 2 **************"<<std::endl;
  display_map(cm);

  // 3) save
  std::ofstream output("output_example_2.xml");
  output<<cm;
}

//******************************************************************************
/* Configuration example 3
   Map containing custom attributes, and overloading write_cmap_attribute_node.
   Here custom attributes are also saved. */
struct ACustomType2
{
  int anint;
  float afloat;
  ACustomType2(int ai=0, float af=0.0) : anint(ai), afloat(af)
  {}
  friend std::ostream& operator<<(std::ostream& os, const ACustomType2& a)
  { return os<<"["<<a.anint<<", "<<a.afloat<<"]"; }
};

namespace CGAL {

// Definition of function allowing to save custon information.
template<>
void write_cmap_attribute_node<ACustomType2>(boost::property_tree::ptree & node,
                                            const ACustomType2& arg)
{
  boost::property_tree::ptree & nValue = node.add("v","");
  nValue.add("v1",arg.anint);
  nValue.add("v2",arg.afloat);
}

}

struct MesAttributs3
{
  template < class CMap>
  struct Dart_wrapper
  {
    typedef CGAL::Dart<3, CMap> Dart;
    typedef CGAL::Cell_attribute<CMap, int> T1;
    typedef CGAL::Cell_attribute<CMap, ACustomType2> T2;
    typedef CGAL::cpp0x::tuple<void,T1,T2> Attributes;
  };
};

typedef CGAL::Combinatorial_map<3, MesAttributs3> CMap_3c;
typedef CMap_3c::Dart_handle Dart_handle_c;
typedef CMap_3c::Dart_const_handle Dart_const_handle_c;

namespace CGAL
{
// Definition of function allowing to same custom darts.
template<>
void write_cmap_dart_node(boost::property_tree::ptree & node,
                          Dart_const_handle_c dh)
{
  static int i=0;
  boost::property_tree::ptree & nValue = node.add("v","");
  nValue.add("myn",i++);
}

}


void example3()
{
  CMap_3c cm;

  // 1) create two tetrahedra 3-sewn
  Dart_handle_c dh1 = CGAL::make_combinatorial_tetrahedron(cm);
  Dart_handle_c dh2 = CGAL::make_combinatorial_tetrahedron(cm);
  cm.sew<3>(dh1,dh2);

  // 2) make attributes
  for (CMap_3c::Dart_range::iterator it(cm.darts().begin()),
         itend(cm.darts().end()); it!=itend; ++it)
  {
    if ( it->attribute<1>()==NULL )
    {
      cm.set_attribute<1>(it, cm.create_attribute<1>());
      it->attribute<1>()->info()=cm.number_of_attributes<1>();
    }
    if( it->attribute<2>()==NULL)
    {
      cm.set_attribute<2>(it, cm.create_attribute<2>());
      it->attribute<2>()->info()=ACustomType2(cm.number_of_attributes<2>(),
                                              cm.number_of_attributes<2>()/2.0);
    }
  }

  std::cout<<"************** example 3 **************"<<std::endl;
  display_map(cm);

  // 3) save
  std::ofstream output("output_example_3.xml");
  output<<cm;
}

//******************************************************************************
/* Configuration example 4
   Linear cell complex without attributes*/
typedef CGAL::Linear_cell_complex<3> LCC_1;

void example4()
{
  LCC_1 lcc;
  LCC_1::Dart_handle d1 = lcc.make_tetrahedron(LCC_1::Point(-1, 0, 0),
                                               LCC_1::Point(0, 2, 0), 
                                               LCC_1::Point(1, 0, 0),
                                               LCC_1::Point(1, 1, 2));
  LCC_1::Dart_handle d2 = lcc.make_tetrahedron(LCC_1::Point(0, 2, -1),
                                               LCC_1::Point(-1, 0, -1),
                                               LCC_1::Point(1, 0, -1),
                                               LCC_1::Point(1, 1, -3));
  lcc.sew<3>(d1, d2);

  std::cout<<"************** example 4 **************"<<std::endl;
  display_lcc(lcc);

  // 3) save  
  std::ofstream output("output_example_4.xml");
  output<<lcc;
}

//******************************************************************************
/* Configuration example 5
   Linear cell complex with attributes*/
struct Myitem7
{
  template<class Refs>
  struct Dart_wrapper
  {
    typedef CGAL::Dart<3, Refs > Dart;    
    typedef CGAL::Cell_attribute_with_point< Refs, int, CGAL::Tag_true> 
    Vertex_attribute;
    typedef CGAL::Cell_attribute<Refs, int> T1;    
    typedef CGAL::cpp0x::tuple<Vertex_attribute, T1> Attributes;
  };
};
typedef CGAL::Linear_cell_complex_traits
<3, CGAL::Exact_predicates_inexact_constructions_kernel> Traits;
typedef CGAL::Linear_cell_complex<3,3,Traits,Myitem7>     LCC_2;

void example5()
{
  LCC_2 lcc;
  LCC_2::Dart_handle d1 = lcc.make_tetrahedron(LCC_2::Point(-1, 0, 0),
                                               LCC_2::Point(0, 2, 0), 
                                               LCC_2::Point(1, 0, 0),
                                               LCC_2::Point(1, 1, 2));
  LCC_2::Dart_handle d2 = lcc.make_tetrahedron(LCC_2::Point(0, 2, -1),
                                               LCC_2::Point(-1, 0, -1),
                                               LCC_2::Point(1, 0, -1),
                                               LCC_2::Point(1, 1, -3));
  lcc.sew<3>(d1, d2);

  // 2) make attributes
  for (LCC_2::Dart_range::iterator it(lcc.darts().begin()),
         itend(lcc.darts().end()); it!=itend; ++it)
  {
    it->attribute<0>()->info()=rand()%10;
    if ( it->attribute<1>()==NULL )
    {
      lcc.set_attribute<1>(it, lcc.create_attribute<1>());
      it->attribute<1>()->info()=lcc.number_of_attributes<1>();
    }
  }
  
  std::cout<<"************** example 5 **************"<<std::endl;
  display_lcc(lcc);

  // 3) save  
  std::ofstream output("output_example_5.xml");
  output<<lcc;
}

//******************************************************************************
//==================================== main ====================================
int main()
{
  example1();
  example2();
  example3();
  example4();
  example5();
  
  return EXIT_SUCCESS;
}
