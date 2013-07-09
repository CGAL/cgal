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
    Point_2 r1(0, 0), r2(5, 1), r3(7,9), r4(10, 15);

    Arrangement_2 arr, arr1;
    Segment_2 segments_p[4], segments_r[4];

    segments_p[0] = Segment_2(p1, p2);
    segments_p[1] = Segment_2(p2, p3);
    segments_p[2] = Segment_2(p3, p4);
    segments_p[3] = Segment_2(p4, p1);

    segments_r[0] = Segment_2(r1, r2);
    segments_r[1] = Segment_2(r2, r3);
    segments_r[2] = Segment_2(r3, r4);
    segments_r[3] = Segment_2(r4, r1);

    Point_2 q_point(4, 4);

    CGAL::insert(arr, &segments_p[0], &segments_p[4]);
    CGAL::insert(arr1, &segments_r[0], &segments_r[2]);

    CGAL::Visibility_2::Simple_visibility_2<Arrangement_2> visibility(arr);
    assert(true == (CGAL::test_is_attached<CGAL::Visibility_2::Simple_visibility_2<Arrangement_2> >(visibility)));
    visibility.detach();
    assert(false == (CGAL::test_is_attached<CGAL::Visibility_2::Simple_visibility_2<Arrangement_2> >(visibility)));
    visibility.attach(arr);
    assert(true == (CGAL::test_is_attached<CGAL::Visibility_2::Simple_visibility_2<Arrangement_2> >(visibility)));

    assert(false == (CGAL::test_are_equal<Arrangement_2>(arr, arr1)));

    return 0;
}