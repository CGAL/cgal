// This is only to test whether the compiler is able to emit a warning when a
// variable (Filter) is used before being initialized.

#include <CGAL/Filtered_exact.h>

void empty_handler(const char*, const char*, const char*, int, const char *)
{
  // Do nothing.
}

CGAL::Interval_nt_advanced add( const CGAL::Interval_nt_advanced & a,
				const CGAL::Interval_nt_advanced & b)
{
    CGAL::Interval_nt_advanced c,d;
    c = a+b;
    d = a+b;
    return c+d;
}

template < class NT >
void test (const NT &)
{
  CGAL::Failure_behaviour backup = CGAL::set_error_behaviour(CGAL::CONTINUE);
  CGAL::Failure_function prev    = CGAL::set_error_handler(empty_handler);
  NT a;
  NT b=1;
  NT c;
  c = b+a; // a is used but not initialized.
  CGAL::set_error_handler(prev);
  CGAL::set_error_behaviour(backup);
}

int main()
{
  test (CGAL::Interval_nt_advanced());
  test (CGAL::Filtered_exact <double, double> ());
  return 0;
}
