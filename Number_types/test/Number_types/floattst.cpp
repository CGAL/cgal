#include <CGAL/config.h>
#include <CGAL/float.h>
#include <CGAL/number_utils.h>

#if defined(BOOST_MSVC)
#  pragma warning(disable:4723)
#endif

int main()
{
    float zero = 0.0;
    float posnormal = (float)1.3;
    float negnormal = -1.0;
    float nan = zero/zero;
    float posinf = posnormal/zero;
    float neginf = negnormal/zero;

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
    return 0;
}
