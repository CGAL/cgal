//
// Created by kabir on 11/04/21.
//
#include "distance_between_two_polygons.h"
#include <CGAL/Polygon_2.h>
#include <CGAL/Point_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <iostream>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef CGAL::Polygon_2<Kernel> Polygon_2;
typedef CGAL::Point_2<Kernel> Point_2;

int main()
{
    Polygon_2 P, Q;
    //test cases
    P.push_back (Point_2 (1, 0));
    P.push_back (Point_2 (2, 0));
    P.push_back (Point_2 (2, 1));
    P.push_back (Point_2 (1, 1));
    Q.push_back (Point_2 (1, 2.5));
    Q.push_back (Point_2 (2, 3 ));
    Q.push_back (Point_2 (1, 5));
    std::cout<<CGAL::distance_between_two_polygons(P,Q);
}


