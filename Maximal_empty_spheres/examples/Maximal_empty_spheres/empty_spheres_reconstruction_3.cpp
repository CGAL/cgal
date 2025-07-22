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
        // I only filter out the positive spheres here as a first demo. The proper version uses positive / negative spheres separately.
        // TODO: implement the proper version
        if (r >= 0) input_spheres.emplace_back(Eigen::RowVector4d(x,y,z,fabs(r)));
        // std::cout << "Sphere: " << x << y << z << r << std::endl;
    }

    std::cout << "Read " << input_spheres.size() << " spheres" << std::endl;

    Eigen::MatrixXd G(input_spheres.size(), 4);
    for (int i=0; i<input_spheres.size(); i++) G.row(i) = input_spheres[i];

    Eigen::MatrixXi contact_indices;
    Eigen::MatrixXd res;
    CGAL::maximal_empty_spheres<CGAL::Dimension_tag<3>>(G, res, &contact_indices);

    // std::cout << res << std::endl;
    // std::cout << "contact_indices:" << std::endl;
    // std::cout << contact_indices << std::endl;

    // contact_indices contains the indices G and contains the spheres that each spheres of res is adjacent to.
    // calculate the contact sphere with largest (absolute) radius for each sphere in G:
    // (careful: results spheres all have negative radius)
    // TODO: filter out the contact spheres outside of a bounding box
    Eigen::VectorXi contact_point_indices = Eigen::VectorXi::Constant(G.rows(),-1);
    Eigen::VectorXd contact_point_radii   = Eigen::VectorXd::Constant(G.rows(),-1.);
    for (int i=0; i<contact_indices.rows(); i++){
        double r= fabs(res(i,3)); // contact spheres have negative radius, use abs value here
        for (int j=0; j<contact_indices.cols(); j++){
            int n = contact_indices(i,j);
            if ( (contact_point_indices(n) < 0) || ((contact_point_indices(n) >= 0) && (contact_point_radii(n) < r)) ){
                contact_point_indices(n) = i;
                contact_point_radii(n)   = r;
            } 
        }
    }

    std::cout << "CP Indices: " << std::endl;
    std::cout << contact_point_indices.transpose() << std::endl;

    std::cout << "CP Radii: " << std::endl;
    std::cout << contact_point_radii.transpose() << std::endl;

    // TODO: calculate the contact points and normals
    // TODO: integrate Poisson Reconstruction

    return 0;
}
