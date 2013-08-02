#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Simple_visibility_2.h>
#include <CGAL/Point_2.h>
#include <CGAL/Segment_2.h>
#include <istream>

typedef CGAL::Exact_predicates_exact_constructions_kernel               Kernel;
typedef Kernel::Point_2                                                 Point_2;
typedef Kernel::Segment_2                                               Segment_2;
typedef CGAL::Arr_segment_traits_2<Kernel>                              Traits_2;
typedef CGAL::Arrangement_2<Traits_2>                                   Arrangement_2;

int main() {
    //create environment
    Point_2 p1(0, 4), p2(0, 0), p3(3, 2), p4(4, 0), p5(4, 4), p6(1, 2);
    Segment_2 s[6];
    s[0] = Segment_2(p1, p2);
    s[1] = Segment_2(p2, p3);
    s[2] = Segment_2(p3, p4);
    s[3] = Segment_2(p4, p5);
    s[4] = Segment_2(p5, p6);
    s[5] = Segment_2(p6, p1);
    Arrangement_2 env;
    CGAL::insert_curves(env, &s[0], &s[6]);
    //visibility query
    Point_2 query_point(0.5, 2);
    Arrangement_2::Face_const_handle face = ++env.faces_begin();

    Arrangement_2 non_regular_output;
    CGAL::Visibility_2::Simple_visibility_2<Arrangement_2, CGAL::Tag_false> non_regular_visibility(env);
    non_regular_visibility.visibility_region(query_point, face, non_regular_output);
    std::cout<<"Non-regularized visibility region of p has "<<non_regular_output.number_of_vertices()<<" vertices."<<std::endl;

    Arrangement_2 regular_output;
    CGAL::Visibility_2::Simple_visibility_2<Arrangement_2, CGAL::Tag_true> regular_visibility(env);
    regular_visibility.visibility_region(query_point, face, regular_output);
    std::cout<<"Regularized visibility region of p has "<<regular_output.number_of_vertices()<<" vertices."<<std::endl;

}

