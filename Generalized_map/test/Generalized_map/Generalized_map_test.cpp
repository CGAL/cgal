#include <CGAL/Generalized_map.h>
#include <CGAL/Cell_attribute.h>

#include "Generalized_map_2_test.h"
#include "Generalized_map_3_test.h"
#include "Generalized_map_4_test.h"

struct Sum_functor
{
  template<class Cell_attribute>
  void operator()(Cell_attribute& ca1,Cell_attribute& ca2)
  { ca1.info()=ca1.info()+ca2.info(); }
};
struct Divide_by_two_functor
{
  template<class Cell_attribute>
  void operator()(Cell_attribute& ca1,Cell_attribute& ca2)
  {
    ca1.info()=(ca1.info()/2);
    ca2.info()=(ca1.info());
  }
};

struct Myitems_2
{
  template <class GMAP>
  struct Dart_wrapper
  {
    typedef void Dart_info;

    typedef CGAL::Cell_attribute<GMAP,int,CGAL::Tag_true,Sum_functor,
                                 Divide_by_two_functor> myattrib;
    typedef std::tuple<myattrib, myattrib, myattrib>
    Attributes;
  };
};

struct Myitems_2c
{
  template <class GMAP>
  struct Dart_wrapper
  {
    typedef int Dart_info;

    typedef CGAL::Cell_attribute<GMAP,int,CGAL::Tag_false,Sum_functor,
                                 Divide_by_two_functor> myattrib;
    typedef std::tuple<myattrib, myattrib, myattrib>
    Attributes;
  };
};
struct Myitems_3
{
  template <class GMAP>
  struct Dart_wrapper
  {
    typedef double* Dart_info;

    typedef CGAL::Cell_attribute<GMAP,int,CGAL::Tag_true,Sum_functor,
                                 Divide_by_two_functor> myattrib;
    typedef std::tuple<myattrib, myattrib, myattrib, myattrib>
    Attributes;
  };
};

struct Myitems_3c
{
  template <class GMAP>
  struct Dart_wrapper
  {
    typedef void Dart_info;

    typedef CGAL::Cell_attribute<GMAP,int,CGAL::Tag_false,Sum_functor,
                                 Divide_by_two_functor> myattrib;
    typedef std::tuple<myattrib, myattrib, myattrib, myattrib>
    Attributes;
  };
};

struct Myitems_4
{
  template <class GMAP>
  struct Dart_wrapper
  {
    typedef double Dart_info;

    typedef CGAL::Cell_attribute<GMAP,int,CGAL::Tag_true,Sum_functor,
                                 Divide_by_two_functor> myattrib;
    typedef std::tuple<myattrib, myattrib, myattrib, myattrib,myattrib>
    Attributes;
  };
};

struct Myitems_4c
{
  template <class GMAP>
  struct Dart_wrapper
  {
    typedef char* Dart_info;

    typedef CGAL::Cell_attribute<GMAP,int,CGAL::Tag_false,Sum_functor,
                                 Divide_by_two_functor> myattrib;
    typedef std::tuple<myattrib, myattrib, myattrib, myattrib,myattrib>
    Attributes;
  };
};

int main()
{
  std::cout<<"Generalized_map start test (v1)."<<std::flush;

  trace_display_msg("\ntest_GMAP_2<GMAP2>");
  typedef CGAL::Generalized_map<2> GMAP2;
  if ( !test_GMAP_2<GMAP2>() )
  {
    std::cout<<" Error during Test_GMAP_2<GMAP2>."<<std::endl;
    return EXIT_FAILURE;
  }

  trace_display_msg("test_GMAP_3<GMAP3>");
  typedef CGAL::Generalized_map<3> GMAP3;
  if ( !test_GMAP_3<GMAP3>() )
  {
    std::cout<<" Error during Test_GMAP_3<GMAP3>."<<std::endl;
    return EXIT_FAILURE;
  }

  trace_display_msg("test_GMAP_4<GMAP4>");
  typedef CGAL::Generalized_map<4> GMAP4;
  if ( !test_GMAP_4<GMAP4>() )
  {
    std::cout<<" Error during Test_GMAP_4<GMAP4>."<<std::endl;
    return EXIT_FAILURE;
  }

  trace_display_msg("test_GMAP_2<GMAP2b>");
  typedef CGAL::Generalized_map<2,Myitems_2> GMAP2b;
  if ( !test_GMAP_2<GMAP2b>() )
  {
    std::cout<<" Error during Test_GMAP_2<GMAP2b>."<<std::endl;
    return EXIT_FAILURE;
  }

  trace_display_msg("test_GMAP_2<GMAP2c>");
  typedef CGAL::Generalized_map<2,Myitems_2c> GMAP2c;
  if ( !test_GMAP_2<GMAP2c>() )
  {
    std::cout<<" Error during Test_GMAP_2<GMAP2c>."<<std::endl;
    return EXIT_FAILURE;
  }

  trace_display_msg("test_GMAP_3<GMAP3b>");
  typedef CGAL::Generalized_map<3,Myitems_3> GMAP3b;
  if ( !test_GMAP_3<GMAP3b>() )
  {
    std::cout<<" Error during Test_GMAP_3<GMAP3b>."<<std::endl;
    return EXIT_FAILURE;
  }

  trace_display_msg("test_GMAP_3<GMAP3c>");
  typedef CGAL::Generalized_map<3,Myitems_3c> GMAP3c;
  if ( !test_GMAP_3<GMAP3c>() )
  {
    std::cout<<" Error during Test_GMAP_3<GMAP3c>."<<std::endl;
    return EXIT_FAILURE;
  }

  trace_display_msg("test_GMAP_4<GMAP4b>");
  typedef CGAL::Generalized_map<4,Myitems_4> GMAP4b;
  if ( !test_GMAP_4<GMAP4b>() )
  {
    std::cout<<" Error during Test_GMAP_4<GMAP4b>."<<std::endl;
    return EXIT_FAILURE;
  }

  trace_display_msg("test_GMAP_4<GMAP4c>");
  typedef CGAL::Generalized_map<4,Myitems_4c> GMAP4c;
  if ( !test_GMAP_4<GMAP4c>() )
  {
    std::cout<<" Error during Test_GMAP_4<GMAP4c>."<<std::endl;
    return EXIT_FAILURE;
  }

  std::cout<<" Success."<<std::endl;
  return EXIT_SUCCESS;
}
