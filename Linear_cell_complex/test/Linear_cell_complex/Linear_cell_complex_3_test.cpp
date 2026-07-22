#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex_for_generalized_map.h>

#include "Linear_cell_complex_3_test.h"

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

  Myattrib(const typename LCC::Point& p=CGAL::ORIGIN): Base(p, 0)
  {}
};

struct MonInfo
{
  MonInfo(long long int i=0) : mnb(i==0?rand():static_cast<int>(i)),
                               ptr(reinterpret_cast<char*>(this))
  {}

  bool operator==(const MonInfo& info) const
  { return mnb==info.mnb && s==info.s && ptr==info.ptr; }

  friend std::ostream& operator<<(std::ostream& os, const MonInfo& i)
  {
    os<<"("<<i.mnb<<", "<<i.s<<", "<<reinterpret_cast<std::size_t>(i.ptr)<<")";
    return os;
  }

  int mnb;
  std::string s;
  char *ptr;
};

struct Min_items: public CGAL::Linear_cell_complex_min_items
{
#ifdef USE_COMPACT_CONTAINER_WITH_INDEX
  typedef CGAL::Tag_true Use_index;
#endif
};

struct Myitems_3b
{
#ifdef USE_COMPACT_CONTAINER_WITH_INDEX
  typedef CGAL::Tag_true Use_index;
  using Index_type=std::uint32_t;
#endif
  template <class LCC>
  struct Dart_wrapper
  {
    typedef MonInfo Dart_info;
    typedef Myattrib<LCC> myattrib;
    typedef std::tuple<myattrib, myattrib, myattrib, myattrib>
    Attributes;
  };
};

struct Myitems_3c
{
#ifdef USE_COMPACT_CONTAINER_WITH_INDEX
  typedef CGAL::Tag_true Use_index;
#endif
  template <class LCC>
  struct Dart_wrapper
  {
    typedef Myattrib<LCC> myattrib;
    typedef std::tuple<myattrib, myattrib, myattrib, myattrib>
    Attributes;
  };
};

int main()
{
  std::cout<<"Linear_cell_complex_3_test start (v1)."<<std::flush;

  // ****************** TEST FOR CMAP ******************
  trace_display_msg("test_LCC_3<LCC3>");
  typedef CGAL::Linear_cell_complex_for_combinatorial_map<3,3,
      CGAL::Linear_cell_complex_traits<3>, Min_items> LCC3;
  if ( !test_LCC_3<LCC3>() )
  {
    std::cout<<" Error during Test_LCC_3<LCC3>."<<std::endl;
    return EXIT_FAILURE;
  }

  trace_display_msg("test_LCC_3<LCC3b>");
  typedef CGAL::Linear_cell_complex_for_combinatorial_map<3,3,
                                    CGAL::Linear_cell_complex_traits<3>,
                                    Myitems_3b> LCC3b;
  if ( !test_LCC_3<LCC3b>() )
  {
    std::cout<<" Error during Test_LCC_3<LCC3b>."<<std::endl;
    return EXIT_FAILURE;
  }

  trace_display_msg("test_LCC_3<LCC3c>");
  typedef CGAL::Linear_cell_complex_for_combinatorial_map<3,3,
                                    CGAL::Linear_cell_complex_traits<3>,
                                    Myitems_3c> LCC3c;
  if ( !test_LCC_3<LCC3c>() )
  {
    std::cout<<" Error during Test_LCC_3<LCC3c>."<<std::endl;
    return EXIT_FAILURE;
  }

  // ****************** TEST FOR GMAP ******************
  trace_display_msg("test_LCC_3<GLCC3>");
  typedef CGAL::Linear_cell_complex_for_generalized_map<3,3,
      CGAL::Linear_cell_complex_traits<3>, Min_items> GLCC3;
  if ( !test_LCC_3<GLCC3>() )
  {
    std::cout<<" Error during Test_LCC_3<GLCC3>."<<std::endl;
    return EXIT_FAILURE;
  }

  trace_display_msg("test_LCC_3<GLCC3b>");
  typedef CGAL::Linear_cell_complex_for_generalized_map<3,3,
                                    CGAL::Linear_cell_complex_traits<3>,
                                    Myitems_3b> GLCC3b;
  if ( !test_LCC_3<GLCC3b>() )
  {
    std::cout<<" Error during Test_LCC_3<GLCC3b>."<<std::endl;
    return EXIT_FAILURE;
  }

  trace_display_msg("test_LCC_3<GLCC3c>");
  typedef CGAL::Linear_cell_complex_for_generalized_map<3,3,
                                    CGAL::Linear_cell_complex_traits<3>,
                                    Myitems_3c> GLCC3c;
  if ( !test_LCC_3<GLCC3c>() )
  {
    std::cout<<" Error during Test_LCC_3<GLCC3c>."<<std::endl;
    return EXIT_FAILURE;
  }

  std::cout<<" Success."<<std::endl;
  return EXIT_SUCCESS;
}
