#include "Linear_cell_complex_2_test.h"
#include "Linear_cell_complex_3_test.h"
#include "Linear_cell_complex_4_test.h"

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
  template <class LCC>
  struct Dart_wrapper
  {
    typedef CGAL::Dart<2, LCC> Dart;

    typedef CGAL::Cell_attribute_with_point<LCC,int,CGAL::Tag_true,Sum_functor,
                                            Divide_by_two_functor> myattrib;
    typedef CGAL::cpp11::tuple<myattrib, myattrib, myattrib>
    Attributes;
  };
};

struct Myitems_2c
{
  template <class LCC>
  struct Dart_wrapper
  {
    typedef CGAL::Dart<2, LCC> Dart;

    typedef CGAL::Cell_attribute_with_point<LCC,int,CGAL::Tag_false,Sum_functor,
                                            Divide_by_two_functor> myattrib;
    typedef CGAL::cpp11::tuple<myattrib, myattrib, myattrib>
    Attributes;
  };
};
struct Myitems_3
{
  template <class LCC>
  struct Dart_wrapper
  {
    typedef CGAL::Dart<3, LCC> Dart;

    typedef CGAL::Cell_attribute_with_point<LCC,int,CGAL::Tag_true,Sum_functor,
                                            Divide_by_two_functor> myattrib;
    typedef CGAL::cpp11::tuple<myattrib, myattrib, myattrib, myattrib>
    Attributes;
  };
};

struct Myitems_3c
{
  template <class LCC>
  struct Dart_wrapper
  {
    typedef CGAL::Dart<3, LCC> Dart;

    typedef CGAL::Cell_attribute_with_point<LCC,int,CGAL::Tag_false,Sum_functor,
                                            Divide_by_two_functor> myattrib;
    typedef CGAL::cpp11::tuple<myattrib, myattrib, myattrib, myattrib>
    Attributes;
  };
};

struct Myitems_4
{
  template <class LCC>
  struct Dart_wrapper
  {
    typedef CGAL::Dart<4, LCC> Dart;

    typedef CGAL::Cell_attribute_with_point<LCC,int,CGAL::Tag_true,Sum_functor,
                                            Divide_by_two_functor> myattrib;
    typedef CGAL::cpp11::tuple<myattrib, myattrib, myattrib, myattrib,myattrib>
    Attributes;
  };
};

struct Myitems_4c
{
  template <class LCC>
  struct Dart_wrapper
  {
    typedef CGAL::Dart<4, LCC> Dart;

    typedef CGAL::Cell_attribute_with_point<LCC,int,CGAL::Tag_false,Sum_functor,
                                            Divide_by_two_functor> myattrib;
    typedef CGAL::cpp11::tuple<myattrib, myattrib, myattrib, myattrib,myattrib>
    Attributes;
  };
};

int main()
{
  std::cout<<"Linear_cell_complex start test (v1)."<<std::flush;

  trace_display_msg("\ntest_LCC_2<LCC2>");
  typedef CGAL::Linear_cell_complex<2> LCC2;
  if ( !test_LCC_2<LCC2>() )
  {
    std::cout<<" Error during Test_LCC_2<LCC2>."<<std::endl;
    return EXIT_FAILURE;
  }
  
  trace_display_msg("test_LCC_3<LCC3>");
  typedef CGAL::Linear_cell_complex<3> LCC3;
  if ( !test_LCC_3<LCC3>() )
  {
    std::cout<<" Error during Test_LCC_3<LCC3>."<<std::endl;
    return EXIT_FAILURE;
  }
  
  trace_display_msg("test_LCC_4<LCC4>");
  typedef CGAL::Linear_cell_complex<4> LCC4;
  if ( !test_LCC_4<LCC4>() )
  {
    std::cout<<" Error during Test_LCC_4<LCC4>."<<std::endl;
    return EXIT_FAILURE;
  }
  
  trace_display_msg("test_LCC_2<LCC2b>");
  typedef CGAL::Linear_cell_complex<2,2,
                                    CGAL::Linear_cell_complex_traits<2>,
                                    Myitems_2> LCC2b;
  if ( !test_LCC_2<LCC2b>() )
  {
    std::cout<<" Error during Test_LCC_2<LCC2b>."<<std::endl;
    return EXIT_FAILURE;
  }
  
  trace_display_msg("test_LCC_2<LCC2c>");
  typedef CGAL::Linear_cell_complex<2,2,
                                    CGAL::Linear_cell_complex_traits<2>,
                                    Myitems_2c> LCC2c;
  if ( !test_LCC_2<LCC2c>() )
  {
    std::cout<<" Error during Test_LCC_2<LCC2c>."<<std::endl;
    return EXIT_FAILURE;
  }

  trace_display_msg("test_LCC_3<LCC3b>");
  typedef CGAL::Linear_cell_complex<3,3,
                                    CGAL::Linear_cell_complex_traits<3>,
                                    Myitems_3> LCC3b;
  if ( !test_LCC_3<LCC3b>() )
  {
    std::cout<<" Error during Test_LCC_3<LCC3b>."<<std::endl;
    return EXIT_FAILURE;
  }

  trace_display_msg("test_LCC_3<LCC3c>");
  typedef CGAL::Linear_cell_complex<3,3,
                                    CGAL::Linear_cell_complex_traits<3>,
                                    Myitems_3c> LCC3c;
  if ( !test_LCC_3<LCC3c>() )
  {
    std::cout<<" Error during Test_LCC_3<LCC3c>."<<std::endl;
    return EXIT_FAILURE;
  }
  
  trace_display_msg("test_LCC_4<LCC4b>");
  typedef CGAL::Linear_cell_complex<4,4,
                                    CGAL::Linear_cell_complex_traits<4>,
                                    Myitems_4> LCC4b;
  if ( !test_LCC_4<LCC4b>() )
  {
    std::cout<<" Error during Test_LCC_4<LCC4b>."<<std::endl;
    return EXIT_FAILURE;
  }

  trace_display_msg("test_LCC_4<LCC4c>");
  typedef CGAL::Linear_cell_complex<4,4,
                                    CGAL::Linear_cell_complex_traits<4>,
                                    Myitems_4c> LCC4c;
  if ( !test_LCC_4<LCC4c>() )
  {
    std::cout<<" Error during Test_LCC_4<LCC4c>."<<std::endl;
    return EXIT_FAILURE;
  }

  std::cout<<" Success."<<std::endl;
  return EXIT_SUCCESS;
}
