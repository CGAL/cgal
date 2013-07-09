/*
 * Author: Francisc Bungiu 
 * E-mail: fbungiu@gmail.com
 */

#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Simple_visibility_2.h>
#include <CGAL/test_model_methods.h>
#include <CGAL/test_utils.h>
#include <CGAL/Simple_visibility_2.h>

int main() {
    typedef CGAL::Gmpq                                              Number_type;
    typedef CGAL::Cartesian<Number_type> 							Kernel;
    typedef CGAL::Arr_segment_traits_2<Kernel> 						Traits_2;
    typedef Traits_2::Point_2										Point_2;
    typedef Traits_2::X_monotone_curve_2							Segment_2;
    typedef CGAL::Arrangement_2<Traits_2>							Arrangement_2;

    Point_2 p1(0, 0), p2(8, 0), p3(8, 8), p4(0, 8);

    Arrangement_2 arr;
    Segment_2 segments[4];

    segments[0] = Segment_2(p1, p2);
    segments[1] = Segment_2(p2, p3);
    segments[2] = Segment_2(p3, p4);
    segments[3] = Segment_2(p4, p1);

    Point_2 q_point(4, 4);

    CGAL::insert(arr, &segments[0], &segments[4]);

    CGAL::Visibility_2::Simple_visibility_2<Arrangement_2> visibility(arr);
    assert(true == (CGAL::test_is_attached<CGAL::Visibility_2::Simple_visibility_2<Arrangement_2> >(visibility)));
    visibility.detach();
    assert(false == (CGAL::test_is_attached<CGAL::Visibility_2::Simple_visibility_2<Arrangement_2> >(visibility)));
    visibility.attach(arr);
    assert(true == (CGAL::test_is_attached<CGAL::Visibility_2::Simple_visibility_2<Arrangement_2> >(visibility)));

    assert(true == (CGAL::test_are_equal<Arrangement_2>(arr, arr)));

    return 0;
}