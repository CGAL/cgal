#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Maximal_empty_spheres/maximal_empty_spheres.h>
#include <iostream>
#include <vector>

using Circle_2 = CGAL::Exact_predicates_inexact_constructions_kernel::Circle_2;
using Point_2 = CGAL::Exact_predicates_inexact_constructions_kernel::Point_2;

int main(){
    std::vector<Circle_2> circles, result;
    std::ifstream in("data/2D/circles.csv");
    double x, y, r;
    while(in >> x){

        in.ignore(10,','); in >> y;  in.ignore(10,',');  in >> r;
        circles.emplace_back(Point_2(x, y), r*r);
        std::cout << "Circle: " << circles.back() << std::endl;

    }
    CGAL::maximal_empty_spheres_2(circles, std::back_inserter(result));
    return 0;
}
