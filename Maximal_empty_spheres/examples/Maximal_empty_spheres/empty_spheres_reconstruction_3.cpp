#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Maximal_empty_spheres/maximal_empty_spheres.h>
#include <iostream>
#include <vector>

// using Sphere_3 = CGAL::Exact_predicates_inexact_constructions_kernel::Sphere_3;
// using Point_3 = CGAL::Exact_predicates_inexact_constructions_kernel::Point_3;

int main(){
    std::ifstream in("data/3D/spheres.csv");
    std::vector<Eigen::RowVector4d> input_spheres;

    double x, y, z, r;
    while(in >> x){

        in.ignore(10,','); in >> y;  in.ignore(10,','); in >> z; in.ignore(10,','); in >> r;
        if (r >= 0) input_spheres.emplace_back(Eigen::RowVector4d(x,y,z,fabs(r)));
        // std::cout << "Sphere: " << x << y << z << r << std::endl;
    }

    std::cout << "Read " << input_spheres.size() << " spheres" << std::endl;

    Eigen::MatrixXd G(input_spheres.size(), 4);
    for (int i=0; i<input_spheres.size(); i++) G.row(i) = input_spheres[i];

    Eigen::MatrixXi contact_indices;
    Eigen::MatrixXd res;
    CGAL::maximal_empty_spheres<CGAL::Dimension_tag<3>>(G, res, &contact_indices);

    std::cout << "Res.shape:             " << res.rows() << ", " << res.cols() << std::endl;
    std::cout << "contact_indices.shape: " << contact_indices.rows() << ", " << contact_indices.cols() << std::endl;
    std::cout << contact_indices << std::endl;

    /*
    std::vector<Sphere_3> spheres, result;
    std::ifstream in("data/3D/spheres.csv");
    double x, y, z, r;
    while(in >> x){

        in.ignore(10,','); in >> y;  in.ignore(10,','); in >> z; in.ignore(10,','); in >> r;
        spheres.emplace_back(Point_3(x, y, z), r*r);
        std::cout << "Sphere: " << spheres.back() << std::endl;

    }
    CGAL::maximal_empty_spheres(spheres, std::back_inserter(result));
    */
    return 0;
}
