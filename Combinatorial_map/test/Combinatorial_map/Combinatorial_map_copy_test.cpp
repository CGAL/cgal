#include <CGAL/Combinatorial_map.h>
#include <CGAL/Cell_attribute.h>
#include "Combinatorial_map_test_iterators.h"
#include <CGAL/HalfedgeDS_default.h>
#if CGAL_USE_OPENMESH
#  include <OpenMesh/Core/IO/MeshIO.hh>
#  include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#  include <CGAL/boost/graph/graph_traits_TriMesh_ArrayKernelT.h>
  typedef OpenMesh::TriMesh_ArrayKernelT</* MyTraits*/> OpenMesh_mesh;
#endif // CGAL_USE_OPENMESH

#include <iostream>
#include <fstream>
#include <cstdlib>

using namespace std;

struct Map_2_dart_items
{
  /// Dart_wrapper defines the type of darts used.
  template < class Refs >
  struct Dart_wrapper
  {
    typedef CGAL::Cell_attribute< Refs, int, CGAL::Tag_true > Int_attrib;
    typedef CGAL::Cell_attribute< Refs, double, CGAL::Tag_true > Double_attrib;

    typedef std::tuple<Double_attrib, void, Double_attrib> Attributes;
  };
};

struct Map_2_dart_max_items_3
{
  /// Dart_wrapper defines the type of darts used.
  template < class Refs >
  struct Dart_wrapper
  {
    typedef int Dart_info;

    typedef CGAL::Cell_attribute< Refs, int, CGAL::Tag_true > Int_attrib;
    typedef CGAL::Cell_attribute< Refs, double, CGAL::Tag_true > Double_attrib;

    typedef std::tuple<Int_attrib, Int_attrib,
          Double_attrib> Attributes;
  };
};

struct Map_3_dart_items_3
{
  /// Dart_wrapper defines the type of darts used.
  template < class Refs >
  struct Dart_wrapper
  {
    typedef CGAL::Cell_attribute< Refs, int, CGAL::Tag_true > Int_attrib;
    typedef CGAL::Cell_attribute< Refs, double, CGAL::Tag_true > Double_attrib;

    typedef std::tuple<Double_attrib, void,
          Int_attrib, Double_attrib> Attributes;
  };
};

struct Map_3_dart_max_items_3
{
  /// Dart_wrapper defines the type of darts used.
  template < class Refs >
  struct Dart_wrapper
  {
    typedef int* Dart_info;

    typedef CGAL::Cell_attribute< Refs, int, CGAL::Tag_true > Int_attrib;
    typedef CGAL::Cell_attribute< Refs, double, CGAL::Tag_true > Double_attrib;

    typedef std::tuple<Int_attrib, Int_attrib,
          Int_attrib, Double_attrib> Attributes;
  };
};

struct MonInfo
{
  MonInfo(long long int i=0) : mnb(i==0?rand():i),
                               ptr(reinterpret_cast<char*>(this))
  {}
  long long int mnb;
  std::string s;
  char *ptr;

  bool operator==(const MonInfo& info) const
  { return mnb==info.mnb && s==info.s && ptr==info.ptr; }
};

class Another_map_3_dart_items_3
{
public:
  /// Dart_wrapper defines the type of darts used.
  template < class Refs >
  struct Dart_wrapper
  {
    typedef MonInfo Dart_info;

    typedef CGAL::Cell_attribute< Refs, int > Int_attrib;

    typedef std::tuple<Int_attrib, void, Int_attrib> Attributes;
  };
};

struct Map_dart_items_4
{
  template < class Refs >
  struct Dart_wrapper
  {
    typedef CGAL::Cell_attribute< Refs, int > Int_attrib;
    typedef CGAL::Cell_attribute< Refs, double > Double_attrib;

    typedef std::tuple<Int_attrib, void,
          Int_attrib, void, Int_attrib>
    Attributes;
  };
};

struct Map_dart_max_items_4
{
  template < class Refs >
  struct Dart_wrapper
  {
    typedef char Dart_info;

    typedef CGAL::Cell_attribute< Refs, int > Int_attrib;
    typedef CGAL::Cell_attribute< Refs, double > Double_attrib;

    typedef std::tuple<Int_attrib, Int_attrib,
          Int_attrib, Int_attrib, Double_attrib>
    Attributes;
  };
};

// noinfo, void, void, void
typedef CGAL::Combinatorial_map<2, CGAL::Generic_map_min_items > Map1;

// noinfo, double, void, double
typedef CGAL::Combinatorial_map<2, Map_2_dart_items > Map2;

// info=int, int, int, double
typedef CGAL::Combinatorial_map<2, Map_2_dart_max_items_3> Map3;

// noinfo, void, void, void, void
typedef CGAL::Combinatorial_map<3, CGAL::Generic_map_min_items > Map4;

// noinfo, double, void, int, double
typedef CGAL::Combinatorial_map<3, Map_3_dart_items_3> Map5;

// info=int*, int, int, int, double
typedef CGAL::Combinatorial_map<3, Map_3_dart_max_items_3> Map6;

// info=MonInfo, int, void, int, void
typedef CGAL::Combinatorial_map<3, Another_map_3_dart_items_3> Map7;

// noinfo, int, void, int, void, int
typedef CGAL::Combinatorial_map<4, Map_dart_items_4> Map8;

// info=char, int, int, int, int, double
typedef CGAL::Combinatorial_map<4, Map_dart_max_items_4> Map9;

struct Traits { typedef int Point_2; typedef int Point_3; };
typedef CGAL::HalfedgeDS_default<Traits> HDS;

////////////////////////////////////////////////////////////////////////////////
template<typename Map, unsigned int i, typename Attr=typename Map::
         template Attribute_type<i>::type>
struct CreateAttributes
{
  static void run(Map& map)
  {
    long long int nb=0;
    for(typename Map::Dart_range::iterator it=map.darts().begin(),
        itend=map.darts().end(); it!=itend; ++it)
    {
      if ( map.template attribute<i>(it)==map.null_handle )
      {
        map.template set_attribute<i>
          (it, map.template create_attribute<i>
           (typename Map::template Attribute_type<i>::type::Info(nb)));
        ++nb;
      }
    }
  }
};

template<typename Map, unsigned int i>
struct CreateAttributes<Map, i, CGAL::Void>
{
  static void run(Map&)
  {
  }
};

template<typename Map, typename Info=typename Map::Dart_info>
struct InitDartInfo
{
  static void run(Map& map)
  {
    long long int nb=0;
    for(typename Map::Dart_range::iterator it=map.darts().begin(),
        itend=map.darts().end(); it!=itend; ++it)
    {
      map.info(it)=Info(nb);
    }
  }
};

template<typename Map>
struct InitDartInfo<Map, CGAL::Void>
{
  static void run(Map&)
  {}
};

template<typename Map, unsigned int i, typename Attr=typename Map::
         template Attribute_type<i>::type>
struct DisplayNumberOfAttribs
{
  static void run(Map& amap)
  {
    std::cout<<i<<"-attributes="<<amap.template number_of_attributes<i>()<<"  ";
  }
};
template<typename Map, unsigned int i>
struct DisplayNumberOfAttribs<Map,i,CGAL::Void>
{
  static void run(Map&)
  {}
};

template<typename Map, unsigned int i, typename Attr=typename Map::
         template Attribute_type<i>::type>
struct DisplayAttribs
{
  static void run(Map& amap)
  {
    std::cout<<i<<"-attributes: ";
    for ( typename Map::template Attribute_range<i>::type::iterator
          it=amap.template attributes<i>().begin(),
          itend=amap.template attributes<i>().end();
          it!=itend; ++it )
    {
      std::cout<<amap.template info<i>(it)<<"; ";
    }
    std::cout<<std::endl;
  }
};
template<typename Map, unsigned int i>
struct DisplayAttribs<Map,i,CGAL::Void>
{
  static void run(Map&)
  {}
};

template<typename Map>
void displayAllAttribs2D(Map& amap, const char* c)
{
  std::cout<<c;
  DisplayAttribs<Map,0>::run(amap);
  DisplayAttribs<Map,1>::run(amap);
  DisplayAttribs<Map,2>::run(amap);
}

template<typename Map>
void displayAllAttribs3D(Map& amap, const char* c)
{
  std::cout<<c;
  DisplayAttribs<Map,0>::run(amap);
  DisplayAttribs<Map,1>::run(amap);
  DisplayAttribs<Map,2>::run(amap);
  DisplayAttribs<Map,3>::run(amap);
}

template<typename Map>
void displayAllAttribs4D(Map& amap, const char* c)
{
  std::cout<<c;
  DisplayAttribs<Map,0>::run(amap);
  DisplayAttribs<Map,1>::run(amap);
  DisplayAttribs<Map,2>::run(amap);
  DisplayAttribs<Map,3>::run(amap);
  DisplayAttribs<Map,4>::run(amap);
}

template<typename Map>
void create2Dmap(Map& map)
{
  for ( int i=0; i<15; ++i )
    map.make_combinatorial_hexahedron();
  CreateAttributes<Map,0>::run(map);
  CreateAttributes<Map,1>::run(map);
  CreateAttributes<Map,2>::run(map);
  CGAL_assertion ( map.is_valid() );
}
template<typename Map>
void create3Dmap(Map& map)
{
  for ( int i=0; i<15; ++i )
    map.make_combinatorial_hexahedron();

  for ( int i=0; i<20; ++i )
  {
    typename Map::Dart_handle d1=map.darts().begin();
    while ( !map.template is_free<3>(d1) ) ++d1;
    typename Map::Dart_handle d2=map.darts().begin();
    while ( !map.template is_sewable<3>(d1, d2) ) ++d2;
    map.template sew<3>(d1,d2);
  }
  InitDartInfo<Map>::run(map);
  CreateAttributes<Map,0>::run(map);
  CreateAttributes<Map,1>::run(map);
  CreateAttributes<Map,2>::run(map);
  CreateAttributes<Map,3>::run(map);
}
template<typename Map>
void create4Dmap(Map& map)
{
  for ( int i=0; i<45; ++i )
    map.make_combinatorial_hexahedron();

  for ( int i=0; i<40; ++i )
  {
    typename Map::Dart_handle d1=map.darts().begin();
    while ( !map.template is_free<3>(d1) ) ++d1;
    typename Map::Dart_handle d2=map.darts().begin();
    while ( !map.template is_sewable<3>(d1, d2) ) ++d2;
    map.template sew<3>(d1,d2);
  }

  for ( int i=0; i<20; ++i )
  {
    typename Map::Dart_handle d1=map.darts().begin();
    while ( !map.template is_free<4>(d1) ) ++d1;
    typename Map::Dart_handle d2=map.darts().begin();
    while ( !map.template is_sewable<4>(d1, d2) ) ++d2;
    map.template sew<4>(d1,d2);
  }
  InitDartInfo<Map>::run(map);
  CreateAttributes<Map,0>::run(map);
  CreateAttributes<Map,1>::run(map);
  CreateAttributes<Map,2>::run(map);
  CreateAttributes<Map,3>::run(map);
  CreateAttributes<Map,4>::run(map);
  CGAL_assertion ( map.is_valid() );
}

bool testCopy()
{
  Map1 map1; create2Dmap(map1);
  Map2 map2; create2Dmap(map2);
  Map3 map3; create2Dmap(map3);

  Map4 map4; create3Dmap(map4);
  Map5 map5; create3Dmap(map5);
  Map6 map6; create3Dmap(map6);
  Map7 map7; create3Dmap(map7);

  Map8 map8; create4Dmap(map8);
  Map9 map9; create4Dmap(map9);

  // Test iterators
  if ( !test_iterators_2(map1) )
  {  assert(false); return false; }
  if ( !test_iterators_3(map4) )
  {  assert(false); return false; }
  if ( !test_iterators_4(map8) )
  {  assert(false); return false; }

  // First copy of same types
  {
  Map1 map1p(map1);
  if ( !map1p.is_valid() || !map1.is_isomorphic_to(map1p) )
  { assert(false); return false; }
  Map2 map2p(map2);
  if ( !map2p.is_valid() || !map2.is_isomorphic_to(map2p) )
  { assert(false); return false; }
  Map3 map3p(map3);
  if ( !map3p.is_valid() || !map3.is_isomorphic_to(map3p) )
  { assert(false); return false; }
  Map4 map4p(map4);
  if ( !map4p.is_valid() || !map4.is_isomorphic_to(map4p) )
  { assert(false); return false; }
  Map5 map5p(map5);
  if ( !map5p.is_valid() || !map5.is_isomorphic_to(map5p) )
  { assert(false); return false; }
  Map6 map6p(map6);
  if ( !map6p.is_valid() || !map6.is_isomorphic_to(map6p) )
  { assert(false); return false; }
  Map7 map7p(map7);
  if ( !map7p.is_valid() || !map7.is_isomorphic_to(map7p) )
  { assert(false); return false; }
  Map8 map8p(map8);
  if ( !map8p.is_valid() || !map8.is_isomorphic_to(map8p) )
  { assert(false); return false; }
  Map9 map9p(map9);
  if ( !map9p.is_valid() || !map9.is_isomorphic_to(map9p) )
  { assert(false); return false; }
  }

  // Second copy of same dimensions but different attributes
  // Maps are still isomorphic but no same attributes
  {
    // 2D
    Map2 map1p(map1); assert(map1p.is_valid());
    if ( !map1.is_isomorphic_to(map1p) ) { assert(false); return false; }

    Map3 map1t(map1); assert(map1t.is_valid());
    if ( !map1.is_isomorphic_to(map1t, false) ) { assert(false); return false; }

    if ( !map1p.is_isomorphic_to(map1t, false) )
    { assert(false); return false; }

    Map1 map2p(map2); assert(map2p.is_valid());
    if ( map2.is_isomorphic_to(map2p) ) { assert(false); return false; }
    if ( !map2.is_isomorphic_to(map2p, false, false) )
    { assert(false); return false; }

    Map3 map2t(map2); assert(map2t.is_valid());
    if ( map2.is_isomorphic_to(map2t) ) { assert(false); return false; }
    if ( !map2.is_isomorphic_to(map2t, false, false) )
    { assert(false); return false; }

    if ( map2p.is_isomorphic_to(map2t) ) { assert(false); return false; }
    if ( !map2p.is_isomorphic_to(map2t, false, false) )
    { assert(false); return false; }

    Map1 map3p(map3); assert(map3p.is_valid());
    if ( map3.is_isomorphic_to(map3p) ) { assert(false); return false; }
    if ( !map3.is_isomorphic_to(map3p, false, false) )
    { assert(false); return false; }

    Map2 map3t(map3); assert(map3t.is_valid());
    if ( map3.is_isomorphic_to(map3t) ) { assert(false); return false; }
    if ( !map3.is_isomorphic_to(map3t, false, false) )
    { assert(false); return false; }

    if ( map3p.is_isomorphic_to(map3t) ) { assert(false); return false; }
    if ( !map3p.is_isomorphic_to(map3t, false, false) )
    { assert(false); return false; }

    assert( map1.is_isomorphic_to(map1p)==map1p.is_isomorphic_to(map1) );
    assert( map1.is_isomorphic_to(map1t)==map1t.is_isomorphic_to(map1) );
    assert( map2.is_isomorphic_to(map2p)==map2p.is_isomorphic_to(map2) );
    assert( map2.is_isomorphic_to(map2t)==map2t.is_isomorphic_to(map2) );
    assert( map3.is_isomorphic_to(map3p)==map3p.is_isomorphic_to(map3) );
    assert( map3.is_isomorphic_to(map3t)==map3t.is_isomorphic_to(map3) );

    // 3D
    Map4 map5a(map5); assert(map5a.is_valid());
    if ( map5.is_isomorphic_to(map5a) ) { assert(false); return false; }
    if ( !map5.is_isomorphic_to(map5a, false, false) )
    { assert(false); return false; }

    Map6 map5b(map5); assert(map5b.is_valid());
    if ( map5.is_isomorphic_to(map5b) ) { assert(false); return false; }
    if ( !map5.is_isomorphic_to(map5b, false, false) )
    { assert(false); return false; }
    assert( map5b.number_of_attributes<0>()==0 &&
            map5b.number_of_attributes<1>()==0 &&
            map5b.number_of_attributes<2>()==map5.number_of_attributes<2>() &&
            map5b.number_of_attributes<3>()==map5.number_of_attributes<3>() );

    Map7 map5c(map5); assert(map5c.is_valid());
    if ( map5.is_isomorphic_to(map5c) ) { assert(false); return false; }
    if ( !map5.is_isomorphic_to(map5c, false, false) )
    { assert(false); return false; }
    assert( map5c.number_of_attributes<0>()==0 &&
            map5c.number_of_attributes<2>()==map5.number_of_attributes<2>() );

    assert( map5.is_isomorphic_to(map5a)==map5a.is_isomorphic_to(map5) );
    assert( map5.is_isomorphic_to(map5b)==map5b.is_isomorphic_to(map5) );
    assert( map5.is_isomorphic_to(map5c)==map5c.is_isomorphic_to(map5) );

    // 4D
    Map8 map9a(map9); assert(map9a.is_valid());
    if ( map9.is_isomorphic_to(map9a) ) { assert(false); return false; }
    if ( !map9.is_isomorphic_to(map9a, false, false) )
    { assert(false); return false; }
    assert( map9a.number_of_attributes<0>()==map9.number_of_attributes<0>() &&
            map9a.number_of_attributes<2>()==map9.number_of_attributes<2>() &&
            map9a.number_of_attributes<4>()==0 );
    assert( map9a.is_isomorphic_to(map9)==map9.is_isomorphic_to(map9a) );

    Map9 map8a(map8); assert(map8a.is_valid());
    if ( map8.is_isomorphic_to(map8a) ) { assert(false); return false; }
    if ( !map8.is_isomorphic_to(map8a, false, false) )
    { assert(false); return false; }
    assert( map8a.number_of_attributes<0>()==map8.number_of_attributes<0>() &&
            map8a.number_of_attributes<1>()==0 &&
            map8a.number_of_attributes<2>()==map8.number_of_attributes<2>() &&
            map8a.number_of_attributes<3>()==0 &&
            map8a.number_of_attributes<4>()==0 );
    assert( map8a.is_isomorphic_to(map8)==map8.is_isomorphic_to(map8a) );

  }

  // Third copy of different dimensions and different attributes
  {
    Map5 map2a(map2); assert(map2a.is_valid());
    if ( map2a.is_isomorphic_to(map2) ) { assert(false); return false; }
    if ( !map2a.is_isomorphic_to(map2, false, false) )
    { assert(false); return false; }
    assert( map2a.number_of_attributes<0>()==map2.number_of_attributes<0>() &&
            map2a.number_of_attributes<2>()==0 &&
            map2a.number_of_attributes<3>()==0 );
    assert( map2a.is_isomorphic_to(map2)==map2.is_isomorphic_to(map2a) );

    Map2 map5a(map5); assert(map5a.is_valid());
    if ( map5a.is_isomorphic_to(map5) ) { assert(false); return false; }
    assert( map5a.number_of_attributes<0>()>=map5.number_of_attributes<0>() &&
            map5a.number_of_attributes<2>()==0 );

    Map5 map9a(map9); assert(map9a.is_valid());
    if ( map9a.is_isomorphic_to(map9) ) { assert(false); return false; }
    assert( map9a.number_of_attributes<0>()==0 &&
            map9a.number_of_attributes<2>()>=map9.number_of_attributes<2>() &&
            map9a.number_of_attributes<3>()==0 );
    assert( map9a.is_isomorphic_to(map9)==map9.is_isomorphic_to(map9a) );

    CGAL::Cast_converter_cmap_attributes<Map9,Map5,0> c0;
    CGAL::Default_converter_cmap_attributes<Map9,Map5,1> c1;
    CGAL::Default_converter_cmap_attributes<Map9,Map5,2> c2;
    CGAL::Cast_converter_cmap_attributes<Map9,Map5,3> c3;

    std::tuple<CGAL::Cast_converter_cmap_attributes<Map9,Map5,0>,
        CGAL::Default_converter_cmap_attributes<Map9,Map5,1>,
        CGAL::Default_converter_cmap_attributes<Map9,Map5,2>,
        CGAL::Cast_converter_cmap_attributes<Map9,Map5,3> > myconverters
        (c0, c1, c2, c3);

    Map5 map9b(map9, myconverters); assert(map9a.is_valid());
    if ( map9b.is_isomorphic_to(map9) ) { assert(false); return false; }
    assert( map9b.number_of_attributes<0>()>=map9.number_of_attributes<0>() &&
            map9b.number_of_attributes<2>()>=map9.number_of_attributes<2>() &&
            map9b.number_of_attributes<3>()>=map9.number_of_attributes<3>() );
    assert( map9b.is_isomorphic_to(map9)==map9.is_isomorphic_to(map9b) );

    CGAL::Cast_converter_cmap_attributes<Map5,Map9,0> cb0;
    CGAL::Default_converter_cmap_attributes<Map5,Map9,1> cb1;
    CGAL::Default_converter_cmap_attributes<Map5,Map9,2> cb2;
    CGAL::Cast_converter_cmap_attributes<Map5,Map9,3> cb3;

    std::tuple<CGAL::Cast_converter_cmap_attributes<Map5,Map9,0>,
        CGAL::Default_converter_cmap_attributes<Map5,Map9,1>,
        CGAL::Default_converter_cmap_attributes<Map5,Map9,2>,
        CGAL::Cast_converter_cmap_attributes<Map5,Map9,3> > myconverters2
        (cb0, cb1, cb2, cb3);

    Map9 map5b(map5, myconverters2); assert(map5b.is_valid());
    if ( map5b.is_isomorphic_to(map5) ) { assert(false); return false; }
    if ( !map5b.is_isomorphic_to(map5, false, false) )
    { assert(false); return false; }
    assert( map5b.number_of_attributes<0>()==map5.number_of_attributes<0>() &&
            map5b.number_of_attributes<2>()==map5.number_of_attributes<2>() &&
            map5b.number_of_attributes<3>()==map5.number_of_attributes<3>() );
    assert( map5b.is_isomorphic_to(map5)==map5.is_isomorphic_to(map5b) );
  }

  { // Test reverse orientation
    Map9 map9b(map9);
    if ( !map9.is_isomorphic_to(map9b) ) { assert(false); return false; }

    map9b.reverse_orientation();
    if ( map9.is_isomorphic_to(map9b) ) { assert(false); return false; }

    map9b.reverse_orientation();
    if ( !map9.is_isomorphic_to(map9b) ) { assert(false); return false; }

    map9b.reverse_orientation_connected_component(map9b.first_dart());
    if ( map9.is_isomorphic_to(map9b) ) { assert(false); return false; }

    map9b.reverse_orientation_connected_component(map9b.first_dart());
    if ( !map9.is_isomorphic_to(map9b) ) { assert(false); return false; }
  }

  /*    displayAllAttribs4D(map9, "map9******************\n");
      displayAllAttribs3D(map9b, "map9b******************\n");*/

  /*   std::cout<<map9b.number_of_attributes<0>()<<"  "
             <<map9.number_of_attributes<0>()<<"  "
            <<map9b.number_of_attributes<2>()<<"  "
           <<map9.number_of_attributes<2>()<<"  "
          <<map9b.number_of_attributes<3>()<<"  "
         <<map9.number_of_attributes<3>()<<std::endl;*/

  //map5.display_characteristics(std::cout)<<std::endl;
  //map5a.display_characteristics(std::cout)<<std::endl;

  return true;
}

bool testImportFromHalfedge()
{
  bool res=true;

  HDS hds;
  CGAL::HalfedgeDS_decorator<HDS> decorator(hds);
  decorator.create_loop();
  decorator.create_segment();

  Map1 map1; map1.import_from_halfedge_graph(hds);
  Map2 map2; map2.import_from_halfedge_graph(hds);
  Map3 map3; map3.import_from_halfedge_graph(hds);

#if CGAL_USE_OPENMESH
  {
    // test the compilation, with an empty mesh
    OpenMesh_mesh hds;
    Map1 map1; map1.import_from_halfedge_graph(hds);
    Map2 map2; map2.import_from_halfedge_graph(hds);
    Map3 map3; map3.import_from_halfedge_graph(hds);
  }
#endif // CGAL_USE_OPENMESH
  return res; // TODO compare number of darts/cells
}

int main()
{
  std::cout<<"Combinatorial map copy test (v1)."<<std::flush;

  if ( !testCopy() )
  {
    std::cout<<" Failed."<<std::endl;
    return EXIT_FAILURE;
  }

  std::cout<<" Success."<<std::endl;
  return EXIT_SUCCESS;
}
