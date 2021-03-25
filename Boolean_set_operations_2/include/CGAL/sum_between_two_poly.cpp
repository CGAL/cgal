#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/minkowski_sum_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Point_2.h>
#include <CGAL/Vector_2.h>
#include <CGAL/Polygon_2_algorithms.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Aff_transformation_2.h>
#include <cmath>
#include <iostream>


typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef CGAL::Polygon_2<Kernel>                           Polygon_2;
typedef CGAL::Polygon_with_holes_2<Kernel>                Polygon_with_holes_2;
typedef CGAL::Point_2<Kernel>                             Point_2;
typedef Kernel::FT                                        FT;
typedef CGAL::Aff_transformation_2<Kernel>                Transformation;

Transformation rotate(CGAL::ROTATION,sin(M_PI),cos(M_PI));

Kernel kernel;
FT squared_distance_from_origin(Polygon_with_holes_2 sum) {
    CGAL::Origin origin;
    auto comp_sqr_distance = kernel.compute_squared_distance_2_object();
    auto n = sum.outer_boundary().size();
    const auto& edge = sum.outer_boundary().edge(1);
    auto sqr_distance = comp_sqr_distance(origin, edge);
    for (auto i = 1; i < n; ++i) {
        const auto& edge= sum.outer_boundary().edge(i);
        auto tmp = comp_sqr_distance(origin, edge);
        if (tmp < sqr_distance) sqr_distance = tmp;
    }
    return sqr_distance;
}
int main()
{
    Polygon_2 P, Q;
    //test cases
    P.push_back (Point_2 (1, 0));
    P.push_back (Point_2 (2, 0));
    P.push_back (Point_2 (2, 1));
    P.push_back (Point_2 (1, 1));
    Q.push_back (Point_2 (1, 3));
    Q.push_back (Point_2 (2, 3 ));
    Q.push_back (Point_2 (1, 5));
    Polygon_2 Q_Reflection;
    for (auto &i : Q) {
        //Reflection
        Q_Reflection.push_back(rotate(i));
    }
    Polygon_with_holes_2 sum = CGAL::minkowski_sum_2(P, Q_Reflection);
    Kernel traits;
    FT x;
    switch (CGAL::bounded_side_2(sum.outer_boundary().begin(),sum.outer_boundary().end(),Point_2(0.0,0.0), traits))
    {
        case CGAL::ON_BOUNDED_SIDE :
            x=0;
            break;
        case CGAL::ON_BOUNDARY:
        case CGAL::ON_UNBOUNDED_SIDE:
            x=squared_distance_from_origin(sum);
            break;
    }
    std::cout<<x<<std::endl;
}

