#include <CGAL/tags.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Simple_visibility_2.h>
#include <vector>

typedef CGAL::Exact_predicates_exact_constructions_kernel               Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>                              Traits_2;
typedef CGAL::Arrangement_2<Traits_2>                                   Arrangement_2;
typedef Arrangement_2::Point_2                                          Point_2;
typedef Arrangement_2::X_monotone_curve_2                               Segment_2;


int main() {
    Point_2 p(0.5, 2);
    Arrangement_2 env;
    int vertice_x[6] = {0, 0, 3, 4, 4, 1};
    int vertice_y[6] = {4, 0, 2, 0, 4, 2};
    std::vector<Point_2> points;
    for (int i = 0; i != 6; i++) {
        points.push_back(Point_2(vertice_x[i], vertice_y[i]));
    }
    for (int i = 0; i != 5; i++) {
        CGAL::insert(env, Segment_2(points[i], points[i+1]));
    }
    Arrangement_2 output;
    CGAL::Visibility_2::Simple_visibility_2<Arrangement_2, CGAL::Tag_false> simple_visibility(env);
    simple_visibility.visibility_region(p, env.faces_begin(), output);
}

