#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Maximal_empty_spheres/maximal_empty_spheres.h>
#include <iostream>
#include <vector>

#include "polyscope/polyscope.h"
#include "polyscope/point_cloud.h"
#include "polyscope/surface_mesh.h"

// using Sphere_3 = CGAL::Exact_predicates_inexact_constructions_kernel::Sphere_3;
// using Point_3 = CGAL::Exact_predicates_inexact_constructions_kernel::Point_3;


// ------------------- PSR Code -----------------
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/IO/facets_in_complex_2_to_triangle_mesh.h>
#include <CGAL/Poisson_reconstruction_function.h>
#include <CGAL/property_map.h>
#include <CGAL/IO/read_points.h>
#include <CGAL/compute_average_spacing.h>
#include <CGAL/Polygon_mesh_processing/distance.h>
#include <boost/iterator/transform_iterator.hpp>
#include <vector>
#include <fstream>

// https://doc.cgal.org/latest/Poisson_surface_reconstruction_3/Poisson_surface_reconstruction_3_2poisson_reconstruction_example_8cpp-example.html
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
typedef std::pair<Point, Vector> Point_with_normal;
typedef CGAL::First_of_pair_property_map<Point_with_normal> Point_map;
typedef CGAL::Second_of_pair_property_map<Point_with_normal> Normal_map;
typedef Kernel::Sphere_3 Sphere;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef CGAL::Poisson_reconstruction_function<Kernel> Poisson_reconstruction_function;
typedef CGAL::Surface_mesh_default_triangulation_3 STr;
typedef CGAL::Surface_mesh_complex_2_in_triangulation_3<STr> C2t3;
typedef CGAL::Implicit_surface_3<Kernel, Poisson_reconstruction_function> Surface_3;

int poisson_reconstruct(const std::vector<Point_with_normal> &points){

    // more less https://doc.cgal.org/latest/Poisson_surface_reconstruction_3/Poisson_surface_reconstruction_3_2poisson_reconstruction_example_8cpp-example.html

    // Poisson options
    FT sm_angle = 20.0; // Min triangle angle in degrees.
    FT sm_radius = 30; // Max triangle size w.r.t. point set average spacing.
    FT sm_distance = 0.375; // Surface Approximation error w.r.t. point set average spacing.
    
    Poisson_reconstruction_function function(points.begin(), points.end(), Point_map(), Normal_map());
    // Computes the Poisson indicator function f()
    // at each vertex of the triangulation.
    if ( ! function.compute_implicit_function() )
      return EXIT_FAILURE;
    // Computes average spacing
    FT average_spacing = CGAL::compute_average_spacing<CGAL::Sequential_tag>
      (points, 6 /* knn = 1 ring */,
       CGAL::parameters::point_map (Point_map()));
    // Gets one point inside the implicit surface
    // and computes implicit function bounding sphere radius.
    Point inner_point = function.get_inner_point();
    Sphere bsphere = function.bounding_sphere();
    FT radius = std::sqrt(bsphere.squared_radius());
    // Defines the implicit surface: requires defining a
    // conservative bounding sphere centered at inner point.
    FT sm_sphere_radius = 5.0 * radius;
    FT sm_dichotomy_error = sm_distance*average_spacing/1000.0; // Dichotomy error must be << sm_distance
    Surface_3 surface(function,
                      Sphere(inner_point,sm_sphere_radius*sm_sphere_radius),
                      sm_dichotomy_error/sm_sphere_radius);
    // Defines surface mesh generation criteria
    CGAL::Surface_mesh_default_criteria_3<STr> criteria(sm_angle,  // Min triangle angle (degrees)
                                                        sm_radius*average_spacing,  // Max triangle size
                                                        sm_distance*average_spacing); // Approximation error
    // Generates surface mesh with manifold option
    STr tr; // 3D Delaunay triangulation for surface mesh generation
    C2t3 c2t3(tr); // 2D complex in 3D Delaunay triangulation
    CGAL::make_surface_mesh(c2t3,                                 // reconstructed mesh
                            surface,                              // implicit surface
                            criteria,                             // meshing criteria
                            CGAL::Manifold_with_boundary_tag());  // require manifold mesh
    
    if(tr.number_of_vertices() == 0)
      return EXIT_FAILURE;

    

    // saves reconstructed surface mesh
    std::ofstream out("tmp.off");
    Polyhedron output_mesh;
    CGAL::facets_in_complex_2_to_triangle_mesh(c2t3, output_mesh);

    std::cout << output_mesh.size_of_vertices() << std::endl;
    // std::cout << output_mesh<< std::endl;
    out << output_mesh;

    /*

    // computes the approximation error of the reconstruction
    double max_dist =
      CGAL::Polygon_mesh_processing::approximate_max_distance_to_point_set
      (output_mesh,
       CGAL::make_range (boost::make_transform_iterator
                         (points.begin(), CGAL::Property_map_to_unary_function<Point_map>()),
                         boost::make_transform_iterator
                         (points.end(), CGAL::Property_map_to_unary_function<Point_map>())),
       4000);
    std::cout << "Max distance to point_set: " << max_dist << std::endl;

    */

    return EXIT_SUCCESS;
}

// ----------------------------------------------------


int main(){
    std::ifstream in("data/3D/spheres.csv");
    std::vector<Eigen::RowVector4d> input_spheres;

    bool filter_contact_spheres_bbx=true;

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

    Eigen::RowVectorXd bbxmin = G.block(0,0,G.rows(),3).colwise().minCoeff();
    Eigen::RowVectorXd bbxmax = G.block(0,0,G.rows(),3).colwise().maxCoeff();

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
            if ((!filter_contact_spheres_bbx)
                ||
                (  (res.block(n,0,1,3).array() >  bbxmin.array()).all()
                && (res.block(n,0,1,3).array() <= bbxmax.array()).all())){
                if ( (contact_point_indices(n) < 0) || ((contact_point_indices(n) >= 0) && (contact_point_radii(n) < r)) ){
                    contact_point_indices(n) = i;
                    contact_point_radii(n)   = r;
                } 
            }
        }
    }

    std::cout << "CP Indices: " << std::endl;
    std::cout << contact_point_indices.transpose() << std::endl;

    std::cout << "CP Radii: " << std::endl;
    std::cout << contact_point_radii.transpose() << std::endl;

    // calculate the contact points and normals
    std::vector<Point_with_normal> Pwns;
    for (int i=0; i<G.rows(); i++){
        if (contact_point_indices[i] >= 0){
            Point Csdf = Point(G(i,0),G(i,1),G(i,2));
            double rsdf = G(i,3);

            Point Ccontact = Point(
                        res(contact_point_indices[i],0),
                        res(contact_point_indices[i],1),
                        res(contact_point_indices[i],2)
                        );
            Vector D = Ccontact-Csdf;
            D /= sqrt(D.squared_length());
            Point  P = Csdf + fabs(rsdf)*D;
            Vector N = (rsdf >= 0)? D: -D;
            Pwns.emplace_back(P,N);
        }
    }

    std::cout << "PSR from " << Pwns.size() << " points with normals" << std::endl;
    poisson_reconstruct(Pwns);

    // POLYSCOPE DEBUGGING
    polyscope::init(); 
    polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::None;
    polyscope::SlicePlane* psPlane = polyscope::addSceneSlicePlane();
    psPlane->setDrawPlane(false);
    psPlane->setDrawWidget(true);

    auto pc_centers = polyscope::registerPointCloud("SDF Spheres", G.block(0,0,G.rows(),3));
    auto q = pc_centers->addScalarQuantity("SDF radius", G.col(3).array().abs()); // add the quantity
    pc_centers->setPointRadiusQuantity(q,false); // set the quantity as the radius

    auto res_centers = polyscope::registerPointCloud("Solutions", res.block(0,0,res.rows(),3));
    res_centers->setEnabled(false);
    auto res_r = res_centers->addScalarQuantity("SDF radius", res.col(3).array().abs()); // add the quantity
    res_centers->setPointRadiusQuantity(res_r,false); // set the quantity as the radius

    Eigen::MatrixXd P_(Pwns.size(),3), N_(Pwns.size(),3);
    for (int i=0; i<Pwns.size(); i++){
        for (int d=0; d<3; d++) {
            P_(i,d) = Pwns[i].first[d];
            N_(i,d) = Pwns[i].second[d];
        }
    }

    auto pc_contactpoints = polyscope::registerPointCloud("Contact Points", P_);
    pc_contactpoints->addVectorQuantity("Normals", N_);

    polyscope::show();

    return 0;
}
