
// This is only to test whether the compiler is able to emit a warning when a
// variable (Filter) is used before being initialized.

#include <CGAL/basic.h>
#include <CGAL/Arithmetic_filter.h>

// Just to look at how good the CGAL_IA_FORCE_TO_DOUBLE macro is compiled.
double force2mem(const double a)
{
  return CGAL_IA_FORCE_TO_DOUBLE(a);
}

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
