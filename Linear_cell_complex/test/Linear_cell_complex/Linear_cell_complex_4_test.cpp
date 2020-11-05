#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex_for_generalized_map.h>

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

template<typename LCC>
struct Myattrib : public  CGAL::Cell_attribute_with_point
    <LCC,int,CGAL::Tag_true,Sum_functor, Divide_by_two_functor>
{
  typedef CGAL::Cell_attribute_with_point
  <LCC,int,CGAL::Tag_true,Sum_functor, Divide_by_two_functor>  Base;

  Myattrib()
  {}

  Myattrib(const typename LCC::Point& p): Base(p, 0)
  {}
};

struct MonInfo
{
  MonInfo(int i=0) : mnb(i==0?rand():i), ptr(reinterpret_cast<char*>(this))
  {}

  bool operator==(const MonInfo& info) const
  { return mnb==info.mnb && s==info.s && ptr==info.ptr; }

  int mnb;
  std::string s;
  char *ptr;
};

struct Myitems_4b
{
  template <class LCC>
  struct Dart_wrapper
  {
    typedef MonInfo Dart_info;

    typedef Myattrib<LCC> myattrib;
    typedef std::tuple<myattrib, myattrib, myattrib, myattrib,myattrib>
    Attributes;
  };
};

struct Myitems_4c
{
  template <class LCC>
  struct Dart_wrapper
  {
    typedef Myattrib<LCC> myattrib;
    typedef std::tuple<myattrib, myattrib, myattrib, myattrib,myattrib>
    Attributes;
  };
};

int main()
{
  std::cout<<"Linear_cell_complex_4_test start (v1)."<<std::flush;

  // ****************** TEST FOR CMAP ******************
  trace_display_msg("test_LCC_4<LCC4>");
  typedef CGAL::Linear_cell_complex_for_combinatorial_map<4> LCC4;
  if ( !test_LCC_4<LCC4>() )
  {
    std::cout<<" Error during Test_LCC_4<LCC4>."<<std::endl;
    return EXIT_FAILURE;
  }

  trace_display_msg("test_LCC_4<LCC4b>");
  typedef CGAL::Linear_cell_complex_for_combinatorial_map<4,4,
                                    CGAL::Linear_cell_complex_traits<4>,
                                    Myitems_4b> LCC4b;
  if ( !test_LCC_4<LCC4b>() )
  {
    std::cout<<" Error during Test_LCC_4<LCC4b>."<<std::endl;
    return EXIT_FAILURE;
  }

  trace_display_msg("test_LCC_4<LCC4c>");
  typedef CGAL::Linear_cell_complex_for_combinatorial_map<4,4,
                                    CGAL::Linear_cell_complex_traits<4>,
                                    Myitems_4c> LCC4c;
  if ( !test_LCC_4<LCC4c>() )
  {
    std::cout<<" Error during Test_LCC_4<LCC4c>."<<std::endl;
    return EXIT_FAILURE;
  }

  // ****************** TEST FOR GMAP ******************
  trace_display_msg("test_LCC_4<GLCC4>");
  typedef CGAL::Linear_cell_complex_for_generalized_map<4> GLCC4;
  if ( !test_LCC_4<GLCC4>() )
  {
    std::cout<<" Error during Test_LCC_4<GLCC4>."<<std::endl;
    return EXIT_FAILURE;
  }

  typedef CGAL::Linear_cell_complex_for_generalized_map<4,4,
                                    CGAL::Linear_cell_complex_traits<4>,
                                    Myitems_4b> GLCC4b;
  trace_display_msg("test_LCC_4<GLCC4b>");
  if ( !test_LCC_4<GLCC4b>() )
  {
    std::cout<<" Error during Test_LCC_4<GLCC4b>."<<std::endl;
    return EXIT_FAILURE;
  }

  trace_display_msg("test_LCC_4<GLCC4c>");
  typedef CGAL::Linear_cell_complex_for_generalized_map<4,4,
                                    CGAL::Linear_cell_complex_traits<4>,
                                    Myitems_4c> GLCC4c;
  if ( !test_LCC_4<GLCC4c>() )
  {
    std::cout<<" Error during Test_LCC_4<GLCC4c>."<<std::endl;
    return EXIT_FAILURE;
  }

  std::cout<<" Success."<<std::endl;
  return EXIT_SUCCESS;
}
