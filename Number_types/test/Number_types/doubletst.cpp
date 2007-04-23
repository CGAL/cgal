#include <CGAL/basic.h>

void test_is_integer()
{
  std::cout << "Testing is_integer(double)" << std::endl;
  assert (! CGAL::is_integer(0.5));
  assert (! CGAL::is_integer(0.25));
  assert (! CGAL::is_integer(0.1));
  assert (! CGAL::is_integer(1e-100));
  assert (! CGAL::is_integer(1.5));
  assert (! CGAL::is_integer(15.1));
  assert (! CGAL::is_integer(0.5));
  assert (! CGAL::is_integer(0.25));
  assert (! CGAL::is_integer(0.1));
  assert (! CGAL::is_integer(1e-100));
  assert (! CGAL::is_integer(1.5));
  assert (! CGAL::is_integer(15.1));

  assert (CGAL::is_integer(0));
  assert (CGAL::is_integer(1));
  assert (CGAL::is_integer(2));
  assert (CGAL::is_integer(1e100));
  assert (CGAL::is_integer(0));
  assert (CGAL::is_integer(1));
  assert (CGAL::is_integer(2));
  assert (CGAL::is_integer(1e100));
}

int main()
{
    double zero = 0.0;
    double posnormal = 1.3;
    double negnormal = -1.0;
    double nan = zero/zero;
    double posinf = posnormal/zero;
    double neginf = negnormal/zero;
    
    if (!CGAL:: is_valid(zero))
	return 1;
    if (!CGAL_NTS is_finite(zero))
	return 1;
    if (!CGAL:: is_valid(posnormal))
	return 1;
    if (!CGAL_NTS is_finite(posnormal))
	return 1;
    if (!CGAL:: is_valid(negnormal))
	return 1;
    if (!CGAL_NTS is_finite(negnormal))
	return 1;
    if (CGAL:: is_valid(nan))
	return 1;
    if (CGAL_NTS is_finite(nan))
	return 1;
    if (!CGAL:: is_valid(posinf))
	return 1;
    if (CGAL_NTS is_finite(posinf))
	return 1;
    if (!CGAL:: is_valid(neginf))
	return 1;
    if (CGAL_NTS is_finite(neginf))
	return 1;

    test_is_integer();

    return 0;
}
