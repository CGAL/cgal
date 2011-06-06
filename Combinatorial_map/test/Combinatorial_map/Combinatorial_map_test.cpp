#include <CGAL/Combinatorial_map.h>
#include <CGAL/Cell_attribute.h>

#include "Combinatorial_map_2_test.h"
#include "Combinatorial_map_3_test.h"

struct Map_2_dart_max_items_3
{
  /// Dart_wrapper defines the type of darts used.
  template < class Refs >
  struct Dart_wrapper
  {
    typedef CGAL::Dart< 2, Refs > Dart;

    typedef CGAL::Cell_attribute< Refs, int > Int_attrib;
    typedef CGAL::Cell_attribute< Refs, double > Double_attrib;
    
    typedef CGAL::cpp0x::tuple<Int_attrib, Int_attrib, 
			       Double_attrib> Attributes;
  };
};

struct Map_3_dart_max_items_3
{
  /// Dart_wrapper defines the type of darts used.
  template < class Refs >
  struct Dart_wrapper
  {
    typedef CGAL::Dart< 3, Refs > Dart;

    typedef CGAL::Cell_attribute< Refs, int > Int_attrib;
    typedef CGAL::Cell_attribute< Refs, double > Double_attrib;
    
    typedef CGAL::cpp0x::tuple<Int_attrib, Int_attrib, 
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

    typedef CGAL::Cell_attribute< Refs, int > Int_attrib;
    
    typedef CGAL::cpp0x::tuple<Int_attrib, void, Int_attrib> 
    Attributes;
  };
};

struct Map_dart_max_items_4
{
  template < class Refs >
  struct Dart_wrapper
  {
    typedef CGAL::Dart< 4, Refs > Dart;

    typedef CGAL::Cell_attribute< Refs, int > Int_attrib;
    typedef CGAL::Cell_attribute< Refs, double > Double_attrib;
    
    typedef CGAL::cpp0x::tuple<Int_attrib, Int_attrib, 
			       Int_attrib, Double_attrib, Double_attrib> 
    Attributes;
  };
};



int main()
{
  typedef CGAL::Combinatorial_map<2,
    CGAL::Combinatorial_map_min_items<2> > Map1;
  test2D<Map1>();
  
  typedef CGAL::Combinatorial_map<2, Map_2_dart_max_items_3 > Map2;
  test2D<Map2>();

  typedef CGAL::Combinatorial_map<2, Map_2_dart_max_items_3> Map3;
  test2D<Map3>();

  typedef CGAL::Combinatorial_map<3,
    CGAL::Combinatorial_map_min_items<3> > Map4;
  test3D<Map4>();
  
  typedef CGAL::Combinatorial_map<3, Map_3_dart_max_items_3> Map5;
  test3D<Map5>();

  typedef CGAL::Combinatorial_map<3, Map_3_dart_max_items_3> Map6;
  test3D<Map6>();

  typedef CGAL::Combinatorial_map<3, Another_map_3_dart_items_3> Map7;
  test3D<Map7>();

  typedef CGAL::Combinatorial_map<3, Another_map_3_dart_items_3> Map8;
  test3D<Map8>();

  typedef CGAL::Combinatorial_map<4, Map_dart_max_items_4> Map9;
  test3D<Map9>();

  return EXIT_SUCCESS;
}
