#include <mpfr.h>
#include <CGAL/MpfiInterval.h>

void test_mpfr () {
	mpfr_t one, two;
	mpfr_inits (one, two, NULL);
	mpfr_set_si (one, 1, GMP_RNDN);
	mpfr_set_si (two, 2, GMP_RNDN);

	CGAL::MpfiInterval interval1 (1);
	CGAL::MpfiInterval interval2 (0.5, 1.5);

	std::cout << "interval1 = " << interval1 << "\tis point: " << interval1.is_point () << std::endl;
	std::cout << "interval2 = " << interval2 << "\tis point: " << interval2.is_point () << std::endl;

	std::cout << "interval1==one: " << (interval1 == one) << std::endl;
	std::cout << "one==interval1: " << (one == interval1) << std::endl;

	std::cout << "interval1<=one: ";
	try { std::cout << (interval1 <= one) << std::endl;
	} catch (CGAL::comparison_overlap_exn &o) { std::cout << o.what () << std::endl; }

	std::cout << "interval1<one: ";
	try { std::cout << (interval1 < one) << std::endl;
	} catch (CGAL::comparison_overlap_exn &o) { std::cout << o.what () << std::endl; }

	std::cout << "interval1==two: " << (interval1 == two) << std::endl;

	std::cout << "interval2==one: ";
	try { std::cout << (interval2 == one) << std::endl;
	} catch (CGAL::comparison_overlap_exn &o) { std::cout << o.what () << std::endl; }

	std::cout << "interval2==two: ";
	try { std::cout << (interval2 == two) << std::endl;
	} catch (CGAL::comparison_overlap_exn &o) { std::cout << o.what () << std::endl; }
}

int main () {
	test_mpfr ();
	return 0;
}
