/*
 * Author: Francisc Bungiu 
 * E-mail: fbungiu@gmail.com
 */

#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Simple_visibility_2.h>
#include <CGAL/test_model_methods.h>
#include <CGAL/test_utils.h>
#include <CGAL/test_simple_polygons.h>
#include <CGAL/Simple_visibility_2.h>

#include <iostream>
#include <fstream>

int main() {
{
    typedef CGAL::Gmpq                                              Number_type;
    typedef CGAL::Cartesian<Number_type> 							Kernel;
    typedef CGAL::Arr_segment_traits_2<Kernel> 						Traits_2;
    typedef Traits_2::Point_2										Point_2;
    typedef Traits_2::X_monotone_curve_2							Segment_2;
    typedef CGAL::Arrangement_2<Traits_2>							Arrangement_2;

     // First read arrangement 
    Arrangement_2 arr;
    std::ifstream input("./data/simple_polygon_test_case_1.in");
    CGAL::create_arrangement_from_input<Arrangement_2>(arr, input);
    CGAL::Visibility_2::Simple_visibility_2<Arrangement_2> visibility;
    assert(false == (CGAL::test_is_attached<CGAL::Visibility_2::Simple_visibility_2<Arrangement_2> >(visibility)));
    visibility.attach(arr);
    assert(true == (CGAL::test_is_attached<CGAL::Visibility_2::Simple_visibility_2<Arrangement_2> >(visibility)));
    visibility.detach();
    assert(false == (CGAL::test_is_attached<CGAL::Visibility_2::Simple_visibility_2<Arrangement_2> >(visibility)));
    visibility.attach(arr);
    assert(true == (CGAL::test_are_equal<Arrangement_2>(arr, visibility.arr())));  

    // Run test cases from https://cgal.geometryfactory.com/CGAL/Members/wiki/Visibility/TestCases
    assert(true == (CGAL::simple_polygon_test_case_1<CGAL::Visibility_2::Simple_visibility_2<Arrangement_2>, Arrangement_2> ()));
}
{
    typedef CGAL::Exact_predicates_exact_constructions_kernel       Kernel;
    typedef CGAL::Arr_segment_traits_2<Kernel>                      Traits_2;
    typedef Traits_2::Point_2                                       Point_2;
    typedef Traits_2::X_monotone_curve_2                            Segment_2;
    typedef CGAL::Arrangement_2<Traits_2>                           Arrangement_2;

    // First read arrangement 
    Arrangement_2 arr;
    std::ifstream input("./data/simple_polygon_test_case_1.in");
    CGAL::create_arrangement_from_input<Arrangement_2>(arr, input);
    CGAL::Visibility_2::Simple_visibility_2<Arrangement_2> visibility;
    assert(false == (CGAL::test_is_attached<CGAL::Visibility_2::Simple_visibility_2<Arrangement_2> >(visibility)));
    visibility.attach(arr);
    assert(true == (CGAL::test_is_attached<CGAL::Visibility_2::Simple_visibility_2<Arrangement_2> >(visibility)));
    visibility.detach();
    assert(false == (CGAL::test_is_attached<CGAL::Visibility_2::Simple_visibility_2<Arrangement_2> >(visibility)));
    visibility.attach(arr);
    assert(true == (CGAL::test_are_equal<Arrangement_2>(arr, visibility.arr())));  

    // Run test cases from https://cgal.geometryfactory.com/CGAL/Members/wiki/Visibility/TestCases
    assert(true == (CGAL::simple_polygon_test_case_1<CGAL::Visibility_2::Simple_visibility_2<Arrangement_2>, Arrangement_2> ()));
}
    return 0;
}