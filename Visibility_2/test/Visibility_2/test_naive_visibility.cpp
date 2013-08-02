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
//#include <CGAL/Rotational_sweep_visibility_2.h>
#include <CGAL/test_model_methods.h>
#include <CGAL/test_utils.h>




#include <iostream>
#include <fstream>
#include <string>

int main() {
    /*
    int case_number = 1;
    //test kernel Cartesian<Gmpq>
    {
        typedef CGAL::Gmpq                                              Number_type;
        typedef CGAL::Cartesian<Number_type> 							Kernel;
        typedef CGAL::Arr_segment_traits_2<Kernel> 						Traits_2;
        typedef Traits_2::Point_2										Point_2;
        typedef Traits_2::X_monotone_curve_2							Segment_2;
        typedef CGAL::Arrangement_2<Traits_2>							Arrangement_2;

        std::cout<<"Kernel: Cartesian<Gmpq>"<<std::endl;
        for (int i=1; i <= case_number; i++) {
            std::cout<<"Test "<<i<<" begins"<<std::endl;
            std::string input_arr_file("data/Arrangement_Test/in");
            input_arr_file += number2string(i);
            std::string ans_file("data/Arrangement_Test/non_regular_out");
            ans_file += number2string(i);

            std::ifstream input(input_arr_file.c_str());
//            std::cout<<"read file over\n";
            std::ifstream ans_input(ans_file.c_str());
//            std::cout<<"read another file\n";
            Arrangement_2 arr_in, arr_ans, arr_vb;
            CGAL::create_arrangement_from_file<Arrangement_2>(arr_in, input);
            CGAL::create_arrangement_from_file<Arrangement_2>(arr_ans, ans_input);

            CGAL::Visibility_2::Naive_visibility_2<Arrangement_2, CGAL::Tag_false> vb(arr_in);
            typename Arrangement_2::Face_const_handle fit = arr_in.faces_begin();
            do {
                if (!fit->is_unbounded()) break;
            } while (++fit != arr_in.faces_end());
            vb.visibility_region(Point_2(0, 0), fit, arr_vb);
            CGAL::Visibility_2::print_arrangement(arr_vb);
            CGAL::Visibility_2::print_arrangement(arr_ans);
            std::cout<<(true == (CGAL::test_are_equal<Arrangement_2>(arr_ans, arr_vb)))<<std::endl;
            std::cout<<CGAL::compare_arr_by_edges(arr_ans, arr_vb)<<std::endl;
        }
    }
    //test kernel Exact_predicates_exact_constructions_kernel
//    {
//        typedef CGAL::Exact_predicates_exact_constructions_kernel       Kernel;
//        typedef CGAL::Arr_segment_traits_2<Kernel>                      Traits_2;
//        typedef Traits_2::Point_2                                       Point_2;
//        typedef Traits_2::X_monotone_curve_2                            Segment_2;
//        typedef CGAL::Arrangement_2<Traits_2>                           Arrangement_2;

//        std::cout<<"Kernel: Exact_predicates_exact_constructions"<<std::endl;
//        for (int i=1; i <= case_number; i++) {
//            std::cout<<"Test "<<i<<" begins"<<std::endl;
//            std::string input_arr_file("data/Arrangement_Test/in");
//            input_arr_file += number2string(i);
//            std::string ans_file("data/Arrangement_Test/non_regular_out");
//            ans_file += number2string(i);

//            std::ifstream input(input_arr_file.c_str());
//            std::ifstream ans_input(ans_file.c_str());
//            Arrangement_2 arr_in, arr_ans, arr_vb;
//            CGAL::create_arrangement_from_file<Arrangement_2>(arr_in, input);
//            CGAL::create_arrangement_from_file<Arrangement_2>(arr_ans, ans_input);

//            CGAL::Visibility_2::Naive_visibility_2<Arrangement_2, CGAL::Tag_false> vb(arr_in);
//            typename Arrangement_2::Face_const_handle fit = arr_in.faces_begin();
//            do {
//                if (!fit->is_unbounded()) break;
//            } while (++fit != arr_in.faces_end());
//            vb.visibility_region(Point_2(0, 0), fit, arr_vb);
//            std::cout<<(true == (CGAL::test_are_equal<Arrangement_2>(arr_ans, arr_vb)))<<std::endl;
//        }
//    }*/
    return 0;
}
