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

bool example1(const char* filename)
{
  std::cout<<"************** example 1 **************"<<std::endl;
  CMap_3a cm;
  std::ifstream input(filename);
  if (!input) return false;
  input>>cm;
  display_map(cm);
  return true;
}

//******************************************************************************
/* Configuration example 2
   Map containing custom attributes, without overload read_cmap_attribute_node.
   In this case, custom attributes are not loaded. Use custom load dart. */
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

namespace CGAL
{
// Definition of function allowing to same custom darts.
template<>
void read_cmap_dart_node<CMap_3b::Dart_handle>
(const boost::property_tree::ptree::value_type & v,
 CMap_3b::Dart_handle dh)
{
  std::cout<<"Read dart "<<v.second.get<int>("myn")<<std::endl;
  }
}


bool example2(const char* filename)
{
  std::cout<<"************** example 2 **************"<<std::endl;
  CMap_3b cm;
  std::ifstream input(filename);
  if (!input) return false;
  input>>cm;
  display_map(cm);
  return true;
}

//******************************************************************************
/* Configuration example 3
   Map containing custom attributes, and overloading read_cmap_attribute_node.
   Here custom attributes are also loaded. */
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

template<>
void read_cmap_attribute_node
(const boost::property_tree::ptree::value_type &v,ACustomType2 &val)
{
  val.anint  = v.second.get<int>("v1");
  val.afloat = v.second.get<float>("v2");
    
  /* Example showing how to iterate through all the son of the node v.
    BOOST_FOREACH( const boost::property_tree::ptree::value_type &v2, v.second)
    {
      if(v2.first=="v1")
      {
        val.anint=boost::lexical_cast< int >(v2.second.data());
      }
      else if (v2.first=="v2")
      {
        val.afloat=boost::lexical_cast< float >(v2.second.data());
      }
    }
  */
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

bool example3(const char* filename)
{
  std::cout<<"************** example 3 **************"<<std::endl;
  CMap_3c cm;
  std::ifstream input(filename);
  if (!input) return false;
  input>>cm;
  display_map(cm);
  return true;
}

//******************************************************************************
/* Configuration example 4
   Map containing custom attributes, and overloading read_cmap_attribute_node.
   But dimension of the map is 2 (instead of 3 in previous examples). */
struct MesAttributs4
{
  template < class CMap>
  struct Dart_wrapper
  {
    typedef CGAL::Dart<2, CMap> Dart;
    typedef CGAL::Cell_attribute<CMap, int> T1;
    typedef CGAL::Cell_attribute<CMap, int> T2;
    typedef CGAL::cpp0x::tuple<void,T1,T2> Attributes;
  };
};

typedef CGAL::Combinatorial_map<2, MesAttributs4> CMap_2;

bool example4(const char* filename)
{
  std::cout<<"************** example 4 **************"<<std::endl;
  CMap_2 cm;
  std::ifstream input(filename);
  if (!input) return false;
  input>>cm;
  display_map(cm);
  return true;
}

//******************************************************************************
/* Configuration example 5
   Map containing custom attributes, and overloading read_cmap_attribute_node.
   But dimension of the map is 5. */
struct MesAttributs5
{
  template < class CMap>
  struct Dart_wrapper
  {
    typedef CGAL::Dart<5, CMap> Dart;
    typedef CGAL::Cell_attribute<CMap, int> T1;
    typedef CGAL::Cell_attribute<CMap, ACustomType2> T2;
    typedef CGAL::cpp0x::tuple<void,T1,T2> Attributes;
  };
};

typedef CGAL::Combinatorial_map<5, MesAttributs5> CMap_5;

bool example5(const char* filename)
{
  std::cout<<"************** example 5 **************"<<std::endl;
  CMap_5 cm;
  std::ifstream input(filename);
  if (!input) return false;
  input>>cm;
  display_map(cm);
  return true;
}

//******************************************************************************
/* Configuration example 6
   Linear cell complex without attributes*/
typedef CGAL::Linear_cell_complex<3> LCC_1;

bool example6(const char* filename)
{
  std::cout<<"************** example 6 **************"<<std::endl;
  LCC_1 cm;
  std::ifstream input(filename);
  if (!input) return false;
  input>>cm;
  display_lcc(cm);
  return true;
}

//******************************************************************************
/* Configuration example 7
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

bool example7(const char* filename)
{
  std::cout<<"************** example 7 **************"<<std::endl;
  LCC_2 cm;
  std::ifstream input(filename);
  if (!input) return false;
  input>>cm;
  display_lcc(cm);
  return true;
}

//******************************************************************************
//==================================== main ====================================
int main(int argc, char* argv[])
{
  if(argc!=2)
  {
    std::cout<<"Usage: a.out filename"<<std::endl
             <<"   filename being an xml combinatorial map file.\n";
    return EXIT_FAILURE;
  }

  if ( !example1(argv[1]) || !example2(argv[1]) || !example3(argv[1]) ||
       !example4(argv[1]) || !example5(argv[1]) || !example6(argv[1]) ||
       !example7(argv[1]))
  {
    std::cout<<"Error reading "<<argv[1]<<std::endl;
    return EXIT_FAILURE;
  }  

  return EXIT_SUCCESS;
}
