#include <CGAL/MpfiInterval.h>

void test_int () {
	int one = 1;
	int two = 2;

	CGAL::MpfiInterval interval1 (1);
	CGAL::MpfiInterval interval2 (0.5, 1.5);

	std::cout << "interval1 = " << interval1 << "\tis point: " << interval1.is_point () << std::endl;
	std::cout << "interval2 = " << interval2 << "\tis point: " << interval2.is_point () << std::endl;

	std::cout << "interval1==one: " << (interval1 == one) << std::endl;

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
	test_int ();
	return 0;
}
