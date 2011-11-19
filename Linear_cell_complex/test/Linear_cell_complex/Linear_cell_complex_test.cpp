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
    typedef CGAL::cpp0x::tuple<myattrib, myattrib, myattrib>
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
    typedef CGAL::cpp0x::tuple<myattrib, myattrib, myattrib, myattrib>
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
    typedef CGAL::cpp0x::tuple<myattrib, myattrib, myattrib, myattrib,myattrib>
    Attributes;
  };
};


int main()
{
  typedef CGAL::Linear_cell_complex<2> LCC2;
  if ( !test_LCC_2<LCC2>() )
  {
    std::cout<<"ERROR during Test_LCC_2<LCC2>."<<std::endl;
    return EXIT_FAILURE;
  }
  
  typedef CGAL::Linear_cell_complex<3> LCC3;
  if ( !test_LCC_3<LCC3>() )
  {
    std::cout<<"ERROR during Test_LCC_3<LCC3>."<<std::endl;
    return EXIT_FAILURE;
  }
  
  typedef CGAL::Linear_cell_complex<4> LCC4;
  if ( !test_LCC_4<LCC4>() )
  {
    std::cout<<"ERROR during Test_LCC_4<LCC4>."<<std::endl;
    return EXIT_FAILURE;
  }
  
  typedef CGAL::Linear_cell_complex<2,2,
                                    CGAL::Linear_cell_complex_traits<2>,
                                    Myitems_2> LCC2b;
  if ( !test_LCC_2<LCC2b>() )
  {
    std::cout<<"ERROR during Test_LCC_2<LCC2b>."<<std::endl;
    return EXIT_FAILURE;
  }
  
  typedef CGAL::Linear_cell_complex<3,3,
                                    CGAL::Linear_cell_complex_traits<3>,
                                    Myitems_3> LCC3b;
  if ( !test_LCC_3<LCC3b>() )
  {
    std::cout<<"ERROR during Test_LCC_3<LCC3b>."<<std::endl;
    return EXIT_FAILURE;
  }
  
  typedef CGAL::Linear_cell_complex<4,4,
                                    CGAL::Linear_cell_complex_traits<4>,
                                    Myitems_4> LCC4b;
  if ( !test_LCC_4<LCC4b>() )
  {
    std::cout<<"ERROR during Test_LCC_4<LCC4b>."<<std::endl;
    return EXIT_FAILURE;
  }
  
  return EXIT_SUCCESS;
}
