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
    
    typedef CGAL::cpp11::tuple<Int_attrib, Int_attrib, 
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
    
    typedef CGAL::cpp11::tuple<Int_attrib, Int_attrib, 
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
    
    typedef CGAL::cpp11::tuple<Int_attrib, void, Int_attrib> 
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
    
    typedef CGAL::cpp11::tuple<Int_attrib, Int_attrib, 
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
  if ( !test3D<Map4>() )
  {
    std::cout<<"ERROR during Test_LCC_4<LCC4b>."<<std::endl;
    return EXIT_FAILURE;
  }

  typedef CGAL::Combinatorial_map<3, Map_3_dart_max_items_3> Map5;
  if ( !test3D<Map5>() )
  {
    std::cout<<"ERROR during test3D<Map5>."<<std::endl;
    return EXIT_FAILURE;
  }

  typedef CGAL::Combinatorial_map<3, Map_3_dart_max_items_3> Map6;
  if ( !test3D<Map6>() )
  {
    std::cout<<"ERROR during test3D<Map6>."<<std::endl;
    return EXIT_FAILURE;
  }

  typedef CGAL::Combinatorial_map<3, Another_map_3_dart_items_3> Map7;
  if ( !test3D<Map7>() )
  {
    std::cout<<"ERROR during test3D<Map7>."<<std::endl;
    return EXIT_FAILURE;
  }

  typedef CGAL::Combinatorial_map<3, Another_map_3_dart_items_3> Map8;
  if ( !test3D<Map8>() )
  {
    std::cout<<"ERROR during test3D<Map8>."<<std::endl;
    return EXIT_FAILURE;
  }

  typedef CGAL::Combinatorial_map<4, Map_dart_max_items_4> Map9;
  if ( !test3D<Map9>() )
  {
    std::cout<<"ERROR during test3D<Map9>."<<std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
