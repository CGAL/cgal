/*
 * Author: Kan Huang
 * E-mail: huangkandiy@gmail.com
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



#include <iostream>
#include <fstream>
#include <string>

int main() {
{
    typedef CGAL::Gmpq                                              Number_type;
    typedef CGAL::Cartesian<Number_type> 							Kernel;
    typedef CGAL::Arr_segment_traits_2<Kernel> 						Traits_2;
    typedef Traits_2::Point_2										Point_2;
    typedef Traits_2::X_monotone_curve_2							Segment_2;
    typedef CGAL::Arrangement_2<Traits_2>							Arrangement_2;



    std::cout<<"Kernel Cartesian<Gmpq>"<<std::endl;
    for (int i=1; i!=6; i++) {
        std::cout<<"Test "<<i<<" begins"<<std::endl;
        std::string input_arr_file, ans_file;
        input_arr_file = "data/Arrangement_Test/in" + std::string(i);
        ans_file = "data/Arrangement_Test/non_regular_out" + std::string(i);
        std::ifstream input(input_arr_file);
        std::ifstream ans_input(ans_file);
        Arrangement_2 arr_in, arr_out;
        CGAL::create_arrangement_from_file<Arrangement_2>(arr_in, input);
        CGAL::create_arrangement_from_file<Arrangement_2>(arr_out, ans_input);

        CGAL::Visibility_2::Naive_visibility_2<Arrangement_2> visibility;


    }


    assert(true == (CGAL::test_are_equal<Arrangement_2>(arr, visibility.arr())));  


    // Run test cases from https://cgal.geometryfactory.com/CGAL/Members/wiki/Visibility/TestCases
    assert(true == (CGAL::simple_polygon_test_case_1<CGAL::Visibility_2::Simple_visibility_2<Arrangement_2>, Arrangement_2> ()));
    assert(true == (CGAL::simple_polygon_test_case_2<CGAL::Visibility_2::Simple_visibility_2<Arrangement_2>, Arrangement_2> ()));


    typedef CGAL::Exact_predicates_exact_constructions_kernel       Kernel_b;
    typedef CGAL::Arr_segment_traits_2<Kernel>                      Traits_2_b;
    typedef Traits_2::Point_2                                       Point_2_b;
    typedef Traits_2::X_monotone_curve_2                            Segment_2_b;
    typedef CGAL::Arrangement_2<Traits_2>                           Arrangement_2_b;


    // First read arrangement 
    Arrangement_2 arr;


    std::ifstream input("./data/simple_polygon_test_case_1.in");
    CGAL::create_arrangement_from_file<Arrangement_2>(arr, input);

    assert(false == (CGAL::test_is_attached<CGAL::Visibility_2::Simple_visibility_2<Arrangement_2> >(visibility)));
    visibility.attach(arr);
    assert(true == (CGAL::test_is_attached<CGAL::Visibility_2::Simple_visibility_2<Arrangement_2> >(visibility)));
    visibility.detach();
    assert(false == (CGAL::test_is_attached<CGAL::Visibility_2::Simple_visibility_2<Arrangement_2> >(visibility)));
    visibility.attach(arr);
    assert(true == (CGAL::test_are_equal<Arrangement_2>(arr, visibility.arr())));  

    // Run test cases from https://cgal.geometryfactory.com/CGAL/Members/wiki/Visibility/TestCases
    assert(true == (CGAL::simple_polygon_test_case_1<CGAL::Visibility_2::Simple_visibility_2<Arrangement_2>, Arrangement_2> ()));
    assert(true == (CGAL::simple_polygon_test_case_2<CGAL::Visibility_2::Simple_visibility_2<Arrangement_2>, Arrangement_2> ()));

    return 0;
}
