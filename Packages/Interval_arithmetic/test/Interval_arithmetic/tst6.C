
// This is only to test whether the compiler is able to emit a warning when a
// variable (Filter) is used before being initialized.

#include <CGAL/basic.h>
#include <CGAL/Arithmetic_filter.h>

template < class NT >
void test (const NT &)
{
  NT a;
  NT b=1;
  NT c;
  c = b+a; // a is used but not initialized.
}

int main()
{
  test (CGAL::Interval_nt_advanced());
  test (CGAL::Filtered_exact <double, double> ());
  return 0;
}
