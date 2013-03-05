#include <CGAL/basic.h>
#include <CGAL/Interval_nt.h>

#if defined(BOOST_MSVC)
#  pragma warning(disable:4723)
#endif

int main()
{
    long double zero = 0.0;
    long double posnormal = 1.3l;
    long double negnormal = -1.0;
    long double nan = zero/zero;
    long double posinf = posnormal/zero;
    long double neginf = negnormal/zero;

    if (!CGAL::is_valid(zero))
	return 1;
    if (!CGAL_NTS is_finite(zero))
	return 1;
    if (!CGAL::is_valid(posnormal))
	return 1;
    if (!CGAL_NTS is_finite(posnormal))
	return 1;
    if (!CGAL::is_valid(negnormal))
	return 1;
    if (!CGAL_NTS is_finite(negnormal))
	return 1;
    if (CGAL::is_valid(nan))
	return 1;
    if (CGAL_NTS is_finite(nan))
	return 1;
    if (!CGAL::is_valid(posinf))
	return 1;
    if (CGAL_NTS is_finite(posinf))
	return 1;
    if (!CGAL::is_valid(neginf))
	return 1;
    if (CGAL_NTS is_finite(neginf))
	return 1;

    if (sizeof(double) != sizeof(long double)) {
        CGAL::Interval_nt<> ia = CGAL_NTS to_interval(posnormal);
        std::cout.precision(30);
        std::cout << posnormal << std::endl;
        std::cout << ia << std::endl;
        if (ia.inf() == ia.sup())
	    return 1;
        ia = CGAL_NTS to_interval(-posnormal);
        std::cout.precision(30);
        std::cout << -posnormal << std::endl;
        std::cout << ia << std::endl;
        if (ia.inf() == ia.sup())
	    return 1;
    }

    return 0;
}
