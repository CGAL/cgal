#include <CGAL/use.h>

template <typename T>
struct Foo {
  typedef typename T::fake fake; 
  // Fake typedef to check that Foo<double> is not instantiated.
};

int test_use_type()
{
  typedef Foo<double> type;

  // If the following line is commented, g++-4.8 -Wall displays that
  // warning:
  //   typedef ‘type’ locally defined but not used [-Wunused-local-typedefs]
  CGAL_USE_TYPE(type);

  return 0;
}

int test_use() {
  int unused;

  // If the following line is commented, g++-4.8 -Wall displays that
  // warning:
  //   unused variable ‘unused’ [-Wunused-variable]
  CGAL_USE(unused);
  return 0;
}

int main() {
  return test_use_type() + test_use();
}
