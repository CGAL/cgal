#include <CGAL/Combinatorial_map.h>
#include <CGAL/Cell_attribute.h>

#include "Face_graph_wrapper_test.h"
#include "Combinatorial_map_2_test.h"
#include "Combinatorial_map_3_test.h"
#include <cstdlib>

struct Min_items: public CGAL::Generic_map_min_items
{
#ifdef USE_COMPACT_CONTAINER_WITH_INDEX
  typedef CGAL::Tag_true Use_index;
#endif
};

struct Map_2_dart_items
{
#ifdef USE_COMPACT_CONTAINER_WITH_INDEX
  typedef CGAL::Tag_true Use_index;
#endif

  /// Dart_wrapper defines the type of darts used.
  template < class Refs >
  struct Dart_wrapper
  {
    typedef void Dart_info;

    typedef CGAL::Cell_attribute< Refs, int, CGAL::Tag_true > Int_attrib;
    typedef CGAL::Cell_attribute< Refs, double, CGAL::Tag_true > Double_attrib;

    typedef std::tuple<Double_attrib, void, Double_attrib> Attributes;
  };
};

struct Map_2_dart_max_items_3
{
#ifdef USE_COMPACT_CONTAINER_WITH_INDEX
  typedef CGAL::Tag_true Use_index;
#endif

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
#ifdef USE_COMPACT_CONTAINER_WITH_INDEX
  typedef CGAL::Tag_true Use_index;
#endif

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
#ifdef USE_COMPACT_CONTAINER_WITH_INDEX
  typedef CGAL::Tag_true Use_index;
#endif

  /// Dart_wrapper defines the type of darts used.
  template < class Refs >
  struct Dart_wrapper
  {
    typedef double Dart_info;

    typedef CGAL::Cell_attribute< Refs, int, CGAL::Tag_true > Int_attrib;
    typedef CGAL::Cell_attribute< Refs, double, CGAL::Tag_true > Double_attrib;

    typedef std::tuple<Int_attrib, Int_attrib,
          Int_attrib, Double_attrib> Attributes;
  };
};

struct MonInfo
{
  MonInfo(int i=0) : mnb(i==0?rand():i), ptr(reinterpret_cast<char*>(this))
  {}
  int mnb;
  std::string s;
  char *ptr;

  bool operator==(const MonInfo& info) const
  { return mnb==info.mnb && s==info.s && ptr==info.ptr; }
};

struct Another_map_3_dart_items_3
{
#ifdef USE_COMPACT_CONTAINER_WITH_INDEX
  typedef CGAL::Tag_true Use_index;
#endif

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
#ifdef USE_COMPACT_CONTAINER_WITH_INDEX
  typedef CGAL::Tag_true Use_index;
#endif

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
#ifdef USE_COMPACT_CONTAINER_WITH_INDEX
  typedef CGAL::Tag_true Use_index;
#endif

  template < class Refs >
  struct Dart_wrapper
  {
    typedef char* Dart_info;

    typedef CGAL::Cell_attribute< Refs, int > Int_attrib;
    typedef CGAL::Cell_attribute< Refs, double > Double_attrib;

    typedef std::tuple<Int_attrib, Int_attrib,
          Int_attrib, Double_attrib, Double_attrib>
    Attributes;
  };
};

typedef CGAL::Combinatorial_map<2, Min_items> Map1;

typedef CGAL::Combinatorial_map<2, Map_2_dart_items> Map2;

typedef CGAL::Combinatorial_map<2, Map_2_dart_max_items_3> Map3;

typedef CGAL::Combinatorial_map<3, Min_items> Map4;

typedef CGAL::Combinatorial_map<3, Map_3_dart_items_3> Map5;

typedef CGAL::Combinatorial_map<3, Map_3_dart_max_items_3> Map6;

typedef CGAL::Combinatorial_map<3, Another_map_3_dart_items_3> Map7;

typedef CGAL::Combinatorial_map<4, Map_dart_items_4> Map8;

typedef CGAL::Combinatorial_map<4, Map_dart_max_items_4> Map9;

int main()
{
  if (!test_face_graph_wrapper())
  {
    std::cout<<"ERROR during test_face_graph_wrapper."<<std::endl;
    return EXIT_FAILURE;
  }

  if ( !test_get_new_mark<Map1>() )
  {
    std::cout<<"ERROR during test_get_new_mark."<<std::endl;
    return EXIT_FAILURE;
  }

  if ( !test2D<Map1>() )
  {
    std::cout<<"ERROR during test2D<Map1>."<<std::endl;
    return EXIT_FAILURE;
  }

  if ( !test2D<Map2>() )
  {
    std::cout<<"ERROR during test2D<Map2>."<<std::endl;
    return EXIT_FAILURE;
  }

  if ( !test2D<Map3>() )
  {
    std::cout<<"ERROR during test2D<Map3>."<<std::endl;
    return EXIT_FAILURE;
  }

  if ( !test3D<Map4>() )
  {
    std::cout<<"ERROR during test3D<Map4>."<<std::endl;
    return EXIT_FAILURE;
  }

  if ( !test3D<Map5>() )
  {
    std::cout<<"ERROR during test3D<Map5>."<<std::endl;
    return EXIT_FAILURE;
  }

  if ( !test3D<Map6>() )
  {
    std::cout<<"ERROR during test3D<Map6>."<<std::endl;
    return EXIT_FAILURE;
  }

  if ( !test3D<Map7>() )
  {
    std::cout<<"ERROR during test3D<Map7>."<<std::endl;
    return EXIT_FAILURE;
  }

  if ( !test3D<Map8>() )
  {
    std::cout<<"ERROR during test3D<Map8>."<<std::endl;
    return EXIT_FAILURE;
  }

  if ( !test3D<Map9>() )
  {
    std::cout<<"ERROR during test3D<Map9>."<<std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
