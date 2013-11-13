#include <CGAL/Linear_cell_complex.h>
#include <CGAL/Cell_attribute_with_point.h>
#include <CGAL/Linear_cell_complex_constructors.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <iostream>
#include <fstream>

using namespace std;

struct Map_2_dart_items
{
  /// Dart_wrapper defines the type of darts used.
  template < class Refs >
  struct Dart_wrapper
  {
    typedef CGAL::Dart< 2, Refs > Dart;

    typedef CGAL::Cell_attribute< Refs, int > Int_attrib;
    typedef CGAL::Cell_attribute< Refs, double > Double_attrib;
    typedef CGAL::Cell_attribute_with_point< Refs, double > Double_attrib_wp;

    typedef CGAL::cpp11::tuple<Double_attrib_wp, void, Double_attrib> Attributes;
  };
};

struct Map_2_dart_max_items_3
{
  /// Dart_wrapper defines the type of darts used.
  template < class Refs >
  struct Dart_wrapper
  {
    typedef CGAL::Dart< 2, Refs > Dart;

    typedef CGAL::Cell_attribute_with_point< Refs, int > Int_attrib_wp;
    typedef CGAL::Cell_attribute< Refs, int > Int_attrib;
    typedef CGAL::Cell_attribute< Refs, double > Double_attrib;

    typedef CGAL::cpp11::tuple<Int_attrib_wp, Int_attrib,
          Double_attrib> Attributes;
  };
};

struct Map_3_dart_items_3
{
  /// Dart_wrapper defines the type of darts used.
  template < class Refs >
  struct Dart_wrapper
  {
    typedef CGAL::Dart< 3, Refs > Dart;

    typedef CGAL::Cell_attribute< Refs, int > Int_attrib;
    typedef CGAL::Cell_attribute< Refs, double > Double_attrib;
    typedef CGAL::Cell_attribute_with_point< Refs, double > Double_attrib_wp;

    typedef CGAL::cpp11::tuple<Double_attrib_wp, void,
          Int_attrib, Double_attrib> Attributes;
  };
};

struct Map_3_dart_max_items_3
{
  /// Dart_wrapper defines the type of darts used.
  template < class Refs >
  struct Dart_wrapper
  {
    typedef CGAL::Dart< 3, Refs > Dart;

    typedef CGAL::Cell_attribute_with_point< Refs, int > Int_attrib_wp;
    typedef CGAL::Cell_attribute< Refs, int > Int_attrib;
    typedef CGAL::Cell_attribute< Refs, double > Double_attrib;

    typedef CGAL::cpp11::tuple<Int_attrib_wp, Int_attrib,
          Int_attrib, Double_attrib> Attributes;
  };
};

class Another_map_3_dart_items_3
{
public:
  /// Dart_wrapper defines the type of darts used.
  template < class Refs >
  struct Dart_wrapper
  {
    typedef CGAL::Dart< 3, Refs > Dart;

    typedef CGAL::Cell_attribute_with_point< Refs, int > Int_attrib_wp;
    typedef CGAL::Cell_attribute< Refs, int > Int_attrib;

    typedef CGAL::cpp11::tuple<Int_attrib_wp, void, Int_attrib> Attributes;
  };
};

struct Map_dart_items_4
{
  template < class Refs >
  struct Dart_wrapper
  {
    typedef CGAL::Dart< 4, Refs > Dart;

    typedef CGAL::Cell_attribute_with_point< Refs, int > Int_attrib_wp;
    typedef CGAL::Cell_attribute< Refs, int > Int_attrib;
    typedef CGAL::Cell_attribute< Refs, double > Double_attrib;

    typedef CGAL::cpp11::tuple<Int_attrib_wp, void,
          Int_attrib, void, Int_attrib>
    Attributes;
  };
};

struct Map_dart_max_items_4
{
  template < class Refs >
  struct Dart_wrapper
  {
    typedef CGAL::Dart< 4, Refs > Dart;

    typedef CGAL::Cell_attribute_with_point< Refs, int > Int_attrib_wp;
    typedef CGAL::Cell_attribute< Refs, int > Int_attrib;
    typedef CGAL::Cell_attribute< Refs, double > Double_attrib;

    typedef CGAL::cpp11::tuple<Int_attrib_wp, Int_attrib,
          Int_attrib, Int_attrib, Double_attrib>
    Attributes;
  };
};

typedef CGAL::Linear_cell_complex_traits
<3, CGAL::Exact_predicates_inexact_constructions_kernel> Traits3_a;

typedef CGAL::Linear_cell_complex_traits
<3, CGAL::Exact_predicates_exact_constructions_kernel> Traits3_b;

typedef CGAL::Linear_cell_complex_traits<4> Traits4_a;

// Point_3, void, void
typedef CGAL::Linear_cell_complex<2,3, Traits3_a, CGAL::Linear_cell_complex_min_items<2> > Map1;

// Point_3+double, void, double
typedef CGAL::Linear_cell_complex<2,3, Traits3_a, Map_2_dart_items > Map2;

// Point_3+int, int, double
typedef CGAL::Linear_cell_complex<2,3, Traits3_b, Map_2_dart_max_items_3> Map3;

// Point_3, void, void, void
typedef CGAL::Linear_cell_complex<3,3, Traits3_a, CGAL::Linear_cell_complex_min_items<3> > Map4;

// Point_3+double, void, int, double
typedef CGAL::Linear_cell_complex<3,3, Traits3_a, Map_3_dart_items_3> Map5;

// Point_3+int, int, int, double
typedef CGAL::Linear_cell_complex<3,3, Traits3_b, Map_3_dart_max_items_3> Map6;

// Point_3+int, void, int, void
typedef CGAL::Linear_cell_complex<3,3, Traits3_b, Another_map_3_dart_items_3> Map7;

// Point_4+int, void, int, void, int
typedef CGAL::Linear_cell_complex<4,4, Traits4_a, Map_dart_items_4> Map8;

// Point_4+int, int, int, int, double
typedef CGAL::Linear_cell_complex<4,4, Traits4_a, Map_dart_max_items_4> Map9;

struct Converter_map9_points_into_map5_points
{
  Map5::Attribute_handle<0>::type operator()
  (const Map9& map1, Map5& map2, Map9::Dart_const_handle dh1,
   Map5::Dart_handle dh2) const
  {
    CGAL_assertion( map1.attribute<0>(dh1)!=map1.null_handle );

    Map5::Attribute_handle<0>::type res = map2.attribute<0>(dh2);
    if ( res==map2.null_handle )
    {
      res = map2.create_attribute<0>();
    }

    const Map9::Point & p = map1.point(dh1);
    map2.point_of_vertex_attribute(res) = Map5::Point(p[0],p[1],p[2]);
    return res;
  }
};


/*
template<typename Map>
typename Map::Dart_handle getRandomDart(Map& map)
{
  int nb = rand()%map.number_of_darts();
  typename Map::Dart_range::iterator it=map.darts().begin();
  for ( int i=0; i<nb; ++i, ++it )
  {}
  return it;
}
*/

template<typename Map, int i, typename Info=
         typename Map::template Attribute_type<i>::type::Info>
struct SetInfoIfNonVoid
{
  static void run(Map& map, typename Map::template Attribute_handle<i>::type attr,
                  int nb)
  {
    map.template info_of_attribute<i>(attr)=nb;
  }
};
template<typename Map, int i>
struct SetInfoIfNonVoid<Map, i, void>
{
  static void run(Map&, typename Map::template Attribute_handle<i>::type,
                  int)
  {}
};

template<typename Map, unsigned int i, typename Attr=typename Map::
         template Attribute_type<i>::type>
struct CreateAttributes
{
  static void run(Map& map)
  {
    int nb=0;
    for(typename Map::Dart_range::iterator it=map.darts().begin(),
        itend=map.darts().end(); it!=itend; ++it)
    {
      if ( map.template attribute<i>(it)==map.null_handle )
      {
        map.template set_attribute<i>(it, map.template create_attribute<i>());
        SetInfoIfNonVoid<Map, i>::run(map, map.template attribute<i>(it), ++nb);
      }
    }
  }
};

template<typename Map, typename Attr>
struct CreateAttributes<Map, 0, Attr>
{
  static void run(Map& amap)
  {
    int nb=0;
    for ( typename Map::template Attribute_range<0>::type::iterator
          it=amap.template attributes<0>().begin(),
          itend=amap.template attributes<0>().end(); it!=itend; ++it )
      SetInfoIfNonVoid<Map, 0>::run(amap, it, ++nb);
  }
};

template<typename Map, unsigned int i>
struct CreateAttributes<Map, i, CGAL::Void>
{
  static void run(Map&)
  {}
};

template<typename Map>
struct CreateAttributes<Map, 0, CGAL::Void>
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
      std::cout<<it->info()<<"; ";
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

template<typename Map,typename Attr,typename Info=typename Attr::Info>
struct DisplayVertexAttrib
{
  static void run(Map& amap)
  {
    std::cout<<"0-attributes: ";
    for ( typename Map::template Attribute_range<0>::type::iterator
          it=amap.template attributes<0>().begin(),
          itend=amap.template attributes<0>().end();
          it!=itend; ++it )
    {
      std::cout<<it->info()<<"; ";
    }
    std::cout<<std::endl;
  }
};
template<typename Map,typename Attr>
struct DisplayVertexAttrib<Map,Attr,void>
{
  static void run(Map&)
  {}
};

template<typename Map, typename Attr>
struct DisplayAttribs<Map,0,Attr>
{
  static void run(Map& amap)
  { DisplayVertexAttrib<Map,Attr>::run(amap); }
};

template<typename Map>
void displayAllAttribs2D(Map& amap, const char* c)
{
  std::cout<<c;
  DisplayAttribs<Map,0>::run(amap);
  DisplayAttribs<Map,1>::run(amap);
  DisplayAttribs<Map,2>::run(amap);

  std::cout<<"Points: ";
  for ( typename Map::template Attribute_range<0>::type::iterator
        it=amap.template attributes<0>().begin(),
        itend=amap.template attributes<0>().end();
        it!=itend; ++it )
  {
    std::cout<<it->point()<<"; ";
  }
  std::cout<<std::endl;
}

template<typename Map>
void displayAllAttribs3D(Map& amap, const char* c)
{
  std::cout<<c;
  DisplayAttribs<Map,0>::run(amap);
  DisplayAttribs<Map,1>::run(amap);
  DisplayAttribs<Map,2>::run(amap);
  DisplayAttribs<Map,3>::run(amap);

  std::cout<<"Points: ";
  for ( typename Map::template Attribute_range<0>::type::iterator
        it=amap.template attributes<0>().begin(),
        itend=amap.template attributes<0>().end();
        it!=itend; ++it )
  {
    std::cout<<it->point()<<"; ";
  }
  std::cout<<std::endl;
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

  std::cout<<"Points: ";
  for ( typename Map::template Attribute_range<0>::type::iterator
        it=amap.template attributes<0>().begin(),
        itend=amap.template attributes<0>().end();
        it!=itend; ++it )
  {
    std::cout<<it->point()<<"; ";
  }
  std::cout<<std::endl;

}

template<typename Map>
void create2Dmap(Map& map)
{
  for ( int i=0; i<15; ++i )
    map.make_tetrahedron(typename Map::Point(i, 0, 0),typename Map::Point(i, 2, 0),
                         typename Map::Point(i+1, 0, 0),typename Map::Point(i+1, 1, 2));
  CreateAttributes<Map,0>::run(map);
  CreateAttributes<Map,1>::run(map);
  CreateAttributes<Map,2>::run(map);
  CGAL_assertion ( map.is_valid() );
}
template<typename Map>
void create3Dmap(Map& map)
{
  for ( int i=0; i<15; ++i )
    map.make_tetrahedron(typename Map::Point(i, 0, 0),typename Map::Point(i, 2, 0),
                         typename Map::Point(i+1, 0, 0),typename Map::Point(i+1, 1, 2));

  for ( int i=0; i<20; ++i )
  {
    typename Map::Dart_handle d1=map.darts().begin();
    while ( !map.template is_free<3>(d1) ) ++d1;
    typename Map::Dart_handle d2=map.darts().begin();
    while ( !map.template is_sewable<3>(d1, d2) ) ++d2;
    map.template sew<3>(d1,d2);
  }
  CreateAttributes<Map,0>::run(map);
  CreateAttributes<Map,1>::run(map);
  CreateAttributes<Map,2>::run(map);
  CreateAttributes<Map,3>::run(map);
  CGAL_assertion ( map.is_valid() );
}

template<typename LCC>
typename LCC::Point apoint(typename LCC::FT x, typename LCC::FT y,
                           typename LCC::FT z, typename LCC::FT t)
{
  std::vector<typename LCC::FT> tab;
  tab.push_back(x); tab.push_back(y);
  tab.push_back(z); tab.push_back(t);
  typename LCC::Point p(4,tab.begin(),tab.end());
  return p;
}

template<typename Map>
void create4Dmap(Map& map)
{
  for ( int i=0; i<45; ++i )
    map.make_tetrahedron(apoint<Map>(i, 0, 0, 0),apoint<Map>(i, 2, 0, 0),
                         apoint<Map>(i+1, 0, 0, 0),apoint<Map>(i+1, 1, 2, 0));

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
    if ( map1.is_isomorphic_to(map1p) ) { assert(false); return false; }
    if ( !map1.is_isomorphic_to(map1p, false) ) { assert(false); return false; }

    Map3 map1t(map1); assert(map1t.is_valid());
    if ( map1.is_isomorphic_to(map1t) ) { assert(false); return false; }
    if ( !map1.is_isomorphic_to(map1t, false) ) { assert(false); return false; }

    if ( map1p.is_isomorphic_to(map1t) ) { assert(false); return false; }
    if ( !map1p.is_isomorphic_to(map1t, false) ) { assert(false); return false; }

    Map1 map2p(map2); assert(map2p.is_valid());
    if ( map2.is_isomorphic_to(map2p) ) { assert(false); return false; }
    if ( !map2.is_isomorphic_to(map2p, false) ) { assert(false); return false; }

    Map3 map2t(map2); assert(map2t.is_valid());
    if ( map2.is_isomorphic_to(map2t) ) { assert(false); return false; }
    if ( !map2.is_isomorphic_to(map2t, false) ) { assert(false); return false; }

    if ( map2p.is_isomorphic_to(map2t) ) { assert(false); return false; }
    if ( !map2p.is_isomorphic_to(map2t, false) ) { assert(false); return false; }

    Map1 map3p(map3); assert(map3p.is_valid());
    if ( map3.is_isomorphic_to(map3p) ) { assert(false); return false; }
    if ( !map3.is_isomorphic_to(map3p, false) ) { assert(false); return false; }

    Map2 map3t(map3); assert(map3t.is_valid());
    if ( map3.is_isomorphic_to(map3t) ) { assert(false); return false; }
    if ( !map3.is_isomorphic_to(map3t, false) ) { assert(false); return false; }

    if ( map3p.is_isomorphic_to(map3t) ) { assert(false); return false; }
    if ( !map3p.is_isomorphic_to(map3t, false) ) { assert(false); return false; }

    assert( map1.is_isomorphic_to(map1p)==map1p.is_isomorphic_to(map1) );
    assert( map1.is_isomorphic_to(map1t)==map1t.is_isomorphic_to(map1) );
    assert( map2.is_isomorphic_to(map2p)==map2p.is_isomorphic_to(map2) );
    assert( map2.is_isomorphic_to(map2t)==map2t.is_isomorphic_to(map2) );
    assert( map3.is_isomorphic_to(map3p)==map3p.is_isomorphic_to(map3) );
    assert( map3.is_isomorphic_to(map3t)==map3t.is_isomorphic_to(map3) );

    // 3D
    Map4 map5a(map5); assert(map5a.is_valid());
    if ( map5.is_isomorphic_to(map5a) ) { assert(false); return false; }
    if ( !map5.is_isomorphic_to(map5a, false) ) { assert(false); return false; }

    Map6 map5b(map5); assert(map5b.is_valid());
    if ( map5.is_isomorphic_to(map5b) ) { assert(false); return false; }
    if ( !map5.is_isomorphic_to(map5b, false) ) { assert(false); return false; }
    assert( map5b.number_of_attributes<0>()==map5.number_of_attributes<0>() &&
            map5b.number_of_attributes<1>()==0 &&
            map5b.number_of_attributes<2>()==map5.number_of_attributes<2>() &&
            map5b.number_of_attributes<3>()==map5.number_of_attributes<3>() );

    Map7 map5c(map5); assert(map5c.is_valid());
    if ( map5.is_isomorphic_to(map5c) ) { assert(false); return false; }
    if ( !map5.is_isomorphic_to(map5c, false) ) { assert(false); return false; }
    assert( map5c.number_of_attributes<0>()==map5.number_of_attributes<0>() &&
            map5c.number_of_attributes<2>()==map5.number_of_attributes<2>() );

    assert( map5.is_isomorphic_to(map5a)==map5a.is_isomorphic_to(map5) );
    assert( map5.is_isomorphic_to(map5b)==map5b.is_isomorphic_to(map5) );
    assert( map5.is_isomorphic_to(map5c)==map5c.is_isomorphic_to(map5) );

    // 4D
    Map8 map9a(map9); assert(map9a.is_valid());
    if ( map9.is_isomorphic_to(map9a) ) { assert(false); return false; }
    if ( !map9.is_isomorphic_to(map9a, false) ) { assert(false); return false; }
    assert( map9a.number_of_attributes<0>()==map9.number_of_attributes<0>() &&
            map9a.number_of_attributes<2>()==map9.number_of_attributes<2>() &&
            map9a.number_of_attributes<4>()==0 );
    assert( map9a.is_isomorphic_to(map9)==map9.is_isomorphic_to(map9a) );

    Map9 map8a(map8); assert(map8a.is_valid());
    if ( map8.is_isomorphic_to(map8a) ) { assert(false); return false; }
    if ( !map8.is_isomorphic_to(map8a, false) ) { assert(false); return false; }
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
    if ( !map2a.is_isomorphic_to(map2, false) ) { assert(false); return false; }
    assert( map2a.number_of_attributes<0>()==map2.number_of_attributes<0>() &&
            map2a.number_of_attributes<2>()==0 &&
            map2a.number_of_attributes<3>()==0 );
    assert( map2a.is_isomorphic_to(map2)==map2.is_isomorphic_to(map2a) );

    Map2 map5a(map5); assert(map5a.is_valid());
    if ( map5a.is_isomorphic_to(map5) ) { assert(false); return false; }
    assert( map5a.number_of_attributes<0>()==map2.number_of_attributes<0>() &&
            map5a.number_of_attributes<2>()==0 );

    Map5 map9a(map9); assert(map9a.is_valid());
    if ( map9a.is_isomorphic_to(map9) ) { assert(false); return false; }
    assert( map9a.number_of_attributes<0>()>=map9.number_of_attributes<0>() &&
            map9a.number_of_attributes<2>()>=map9.number_of_attributes<2>() &&
            map9a.number_of_attributes<3>()==0 );
    assert( map9a.is_isomorphic_to(map9)==map9.is_isomorphic_to(map9a) );

    CGAL::Cast_converter_cmap_attributes<Map9,Map5,0> c0;
    CGAL::Default_converter_cmap_attributes<Map9,Map5,1> c1;
    CGAL::Default_converter_cmap_attributes<Map9,Map5,2> c2;
    CGAL::Cast_converter_cmap_attributes<Map9,Map5,3> c3;

    CGAL::cpp11::tuple<CGAL::Cast_converter_cmap_attributes<Map9,Map5,0>,
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

    Converter_map9_points_into_map5_points mypointconverter;

    Map5 map9c(map9, myconverters, mypointconverter); assert(map9a.is_valid());
    if ( map9c.is_isomorphic_to(map9) ) { assert(false); return false; }
    assert( map9c.number_of_attributes<0>()>=map9.number_of_attributes<0>() &&
            map9c.number_of_attributes<2>()>=map9.number_of_attributes<2>() &&
            map9c.number_of_attributes<3>()>=map9.number_of_attributes<3>() );
    assert( map9c.is_isomorphic_to(map9)==map9.is_isomorphic_to(map9c) );

    CGAL::Cast_converter_cmap_attributes<Map5,Map9,0> cb0;
    CGAL::Default_converter_cmap_attributes<Map5,Map9,1> cb1;
    CGAL::Default_converter_cmap_attributes<Map5,Map9,2> cb2;
    CGAL::Cast_converter_cmap_attributes<Map5,Map9,3> cb3;

    CGAL::cpp11::tuple<CGAL::Cast_converter_cmap_attributes<Map5,Map9,0>,
        CGAL::Default_converter_cmap_attributes<Map5,Map9,1>,
        CGAL::Default_converter_cmap_attributes<Map5,Map9,2>,
        CGAL::Cast_converter_cmap_attributes<Map5,Map9,3> > myconverters2
        (cb0, cb1, cb2, cb3);

    Map9 map5b(map5, myconverters2); assert(map5b.is_valid());
    if ( map5b.is_isomorphic_to(map5) ) { assert(false); return false; }
    if ( !map5b.is_isomorphic_to(map5, false) ) { assert(false); return false; }
    assert( map5b.number_of_attributes<0>()==map5.number_of_attributes<0>() &&
            map5b.number_of_attributes<2>()==map5.number_of_attributes<2>() &&
            map5b.number_of_attributes<3>()==map5.number_of_attributes<3>() );
    assert( map5b.is_isomorphic_to(map5)==map5.is_isomorphic_to(map5b) );
  }

  /*map2.display_characteristics(std::cout)<<std::endl;
    map2a.display_characteristics(std::cout)<<std::endl;
    displayAllAttribs2D(mapXX, "mapXX******************\n");
    displayAllAttribs2D(mapYY, "mapYY******************\n");*/
  /*map5.display_characteristics(std::cout)<<std::endl;
    map5b.display_characteristics(std::cout)<<std::endl;
    displayAllAttribs3D(map5, "map5******************\n");
    displayAllAttribs3D(map5b, "map5b******************\n");*/

  return true;
}

int main()
{
  std::cout<<"Linear cell complex copy test (v1)."<<std::flush;

  if ( !testCopy() )
  {
    std::cout<<" Failed."<<std::endl;
    return EXIT_FAILURE;
  }

  std::cout<<" Success."<<std::endl;
  return EXIT_SUCCESS;
}
