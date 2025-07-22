#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Maximal_empty_spheres/maximal_empty_spheres.h>
#include <iostream>
#include <vector>

using Sphere_3 = CGAL::Exact_predicates_inexact_constructions_kernel::Sphere_3;
using Point_3 = CGAL::Exact_predicates_inexact_constructions_kernel::Point_3;

int main(){
    std::vector<Sphere_3> spheres, result;
    std::ifstream in("data/3D/spheres.csv");
    double x, y, z, r;
    while(in >> x){

        in.ignore(10,','); in >> y;  in.ignore(10,','); in >> z; in.ignore(10,','); in >> r;
        spheres.emplace_back(Point_3(x, y, z), r*r);
        std::cout << "Sphere: " << spheres.back() << std::endl;

    }
    CGAL::maximal_empty_spheres(spheres, std::back_inserter(result));
    return 0;
}
