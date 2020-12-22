#include <CGAL/Combinatorial_map.h>
#include <CGAL/Cell_attribute.h>
#include <CGAL/Face_graph_wrapper.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include "Combinatorial_map_2_test.h"
#include "Combinatorial_map_3_test.h"
#include <cstdlib>

struct Map_2_dart_items
{
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
    typedef char* Dart_info;

    typedef CGAL::Cell_attribute< Refs, int > Int_attrib;
    typedef CGAL::Cell_attribute< Refs, double > Double_attrib;

    typedef std::tuple<Int_attrib, Int_attrib,
          Int_attrib, Double_attrib, Double_attrib>
    Attributes;
  };
};

typedef CGAL::Combinatorial_map<2, CGAL::Generic_map_min_items > Map1;

typedef CGAL::Combinatorial_map<2, Map_2_dart_items > Map2;

typedef CGAL::Combinatorial_map<2, Map_2_dart_max_items_3> Map3;

typedef CGAL::Combinatorial_map<3, CGAL::Generic_map_min_items > Map4;

typedef CGAL::Combinatorial_map<3, Map_3_dart_items_3> Map5;

typedef CGAL::Combinatorial_map<3, Map_3_dart_max_items_3> Map6;

typedef CGAL::Combinatorial_map<3, Another_map_3_dart_items_3> Map7;

typedef CGAL::Combinatorial_map<4, Map_dart_items_4> Map8;

typedef CGAL::Combinatorial_map<4, Map_dart_max_items_4> Map9;

bool test_get_new_mark()
{
  cout << "***************************** TEST GET_NEW_MARK:"
       << endl;

  Map1 map;

  Map1::size_type marks[Map1::NB_MARKS];
  for (Map1::size_type i=0; i<Map1::NB_MARKS; ++i)
  {
    try
    {
      marks[i] = map.get_new_mark();
    }
    catch (Map1::Exception_no_more_available_mark)
    {
      std::cerr<<"No more free mark, exit."<<std::endl;
      return false;
    }
  }

  cout << "Creation of NB_MARK marks: OK" << endl;

  bool res = false;
  Map1::size_type mark=0;
  try
  {
    mark = map.get_new_mark();
  }
  catch (Map1::Exception_no_more_available_mark)
  {
    std::cout<<"The creation of an additional mark throw an exception: OK"<<std::endl;
    res = true;
  }

  if ( !res )
  {
      std::cerr<<"PB we can reserve NB_MARK+1 !! mark, exit."<<std::endl;
      map.free_mark(mark); // This is never supposed to occur.
      return false;
  }

  for (Map1::size_type i=0; i<Map1::NB_MARKS; ++i)
  {
    map.free_mark(marks[i]);
  }

  cout << "***************************** TEST GET_NEW_MARK DONE" << endl;

  return true;
}

bool test_face_graph_wrapper()
{
  bool res=true;

  typedef CGAL::Surface_mesh<CGAL::Simple_cartesian<double>::Point_3> SMesh;
  SMesh m;
  std::ifstream in1("data/head.off");
  if (in1.fail())
  {
    std::cout<<"Error: impossible to open 'data/head.off'"<<std::endl;
    return false;
  }
  in1>>m;

  CGAL::Face_graph_wrapper<SMesh> fgw1(m);
  std::vector<unsigned int> cells=fgw1.count_all_cells();
  if (cells[0]!=1487 || cells[1]!=4406 || cells[2]!=2921 ||
      fgw1.number_of_darts()!=8812)
  {
    std::cout<<"Error: incorrect number of cells in test_face_graph_wrapper "
             <<"for Surface_mesh: "
             <<cells[0]<<", "<<cells[1]<<", "<<cells[2]<<", "<<fgw1.number_of_darts()
             <<std::endl;
    res=false;
  }

  typedef CGAL::Polyhedron_3<CGAL::Simple_cartesian<double> > Polyhedron;
  Polyhedron p;
  std::ifstream in2("data/head.off");
  if (in2.fail())
  {
    std::cout<<"Error: impossible to open 'data/head.off'"<<std::endl;
    return false;
  }
  in2>>p;
  CGAL::Face_graph_wrapper<Polyhedron> fgw2(p);
  cells=fgw2.count_all_cells();
  if (cells[0]!=1487 || cells[1]!=4406 || cells[2]!=2921 ||
      fgw2.number_of_darts()!=8812)
  {
    std::cout<<"Error: incorrect number of cells in test_face_graph_wrapper "
             <<"for Polyhedron."
             <<cells[0]<<", "<<cells[1]<<", "<<cells[2]<<", "<<fgw2.number_of_darts()
             <<std::endl;
    res=false;
  }

  return res;
}

int main()
{
  if (!test_face_graph_wrapper())
  {
    std::cout<<"ERROR during test_face_graph_wrapper."<<std::endl;
    return EXIT_FAILURE;
  }

  if ( !test_get_new_mark() )
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
