#include <CGAL/basic.h>

int main()
{
    float zero = 0.0;
    float posnormal = 1.3;
    float negnormal = -1.0;
    float nan = zero/zero;
    float posinf = posnormal/zero;
    float neginf = negnormal/zero;
    
    if (!CGAL::is_valid(zero))
	return 1;
    if (!CGAL::is_finite(zero))
	return 1;
    if (!CGAL::is_valid(posnormal))
	return 1;
    if (!CGAL::is_finite(posnormal))
	return 1;
    if (!CGAL::is_valid(negnormal))
	return 1;
    if (!CGAL::is_finite(negnormal))
	return 1;
    if (CGAL::is_valid(nan))
	return 1;
    if (CGAL::is_finite(nan))
	return 1;
    if (!CGAL::is_valid(posinf))
	return 1;
    if (CGAL::is_finite(posinf))
	return 1;
    if (!CGAL::is_valid(neginf))
	return 1;
    if (CGAL::is_finite(neginf))
	return 1;
    return 0;
}
