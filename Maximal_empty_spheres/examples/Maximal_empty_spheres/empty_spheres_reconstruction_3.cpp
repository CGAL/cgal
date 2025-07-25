#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Maximal_empty_spheres/maximal_empty_spheres.h>
#include <CGAL/Dimension.h>

#include <iostream>
#include <fstream>
#include <vector>

#include <Eigen/Core>
// ------------------- PSR Code -----------------
#include <CGAL/poisson_surface_reconstruction.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/IO/polygon_mesh_io.h>
#include <CGAL/compute_average_spacing.h>
#include <CGAL/Timer.h>


// https://cgal.geometryfactory.com/CGAL/doc/master/Poisson_surface_reconstruction_3/Poisson_surface_reconstruction_3_2poisson_reconstruction_function_8cpp-example.htmltypedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
typedef std::pair<Point, Vector> Point_with_normal;
typedef CGAL::First_of_pair_property_map<Point_with_normal> Point_map;
typedef CGAL::Second_of_pair_property_map<Point_with_normal> Normal_map;
typedef CGAL::Surface_mesh<Point> Mesh;

void contact_points(const Eigen::MatrixXd &G, std::vector<Point_with_normal> &Pwns, Eigen::MatrixXd *bbxl=NULL){
    Eigen::MatrixXi contact_indices;
    Eigen::MatrixXd res;
    Eigen::MatrixXd G_ = G;
    G_.col(3) = G.col(3).array().abs();
    CGAL::maximal_empty_spheres<CGAL::Dimension_tag<3>>(G_, res, &contact_indices);

    // std::cout << res << std::endl;
    // std::cout << "contact_indices:" << std::endl;
    // std::cout << contact_indices << std::endl;

    // contact_indices contains the indices G and contains the spheres that each spheres of res is adjacent to.
    // calculate the contact sphere with largest (absolute) radius for each sphere in G:
    // (careful: results spheres all have negative radius)
    Eigen::VectorXi contact_point_indices = Eigen::VectorXi::Constant(G.rows(),-1);
    Eigen::VectorXd contact_point_radii   = Eigen::VectorXd::Constant(G.rows(),-1.);
    for (int i=0; i<contact_indices.rows(); i++){
        double r = fabs(res(i,3)); // contact spheres have negative radius, use abs value here
        if ((!bbxl)
            ||
            (  (res.block(i,0,1,3).array() >  bbxl->row(0).array()).all()
            && (res.block(i,0,1,3).array() <= bbxl->row(1).array()).all())){

            for (int j=0; j<contact_indices.cols(); j++){
                int n = contact_indices(i,j);
                if ( (contact_point_indices(n) < 0) || ((contact_point_indices(n) >= 0) && (contact_point_radii(n) < r)) ){
                    contact_point_indices(n) = i;
                    contact_point_radii(n)   = r;
                }
            }
        }
    }

    // std::cout << "CP Indices: " << std::endl;
    // std::cout << contact_point_indices.transpose() << std::endl;

    // std::cout << "CP Radii: " << std::endl;
    // std::cout << contact_point_radii.transpose() << std::endl;

    // calculate the contact points and normals
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
            double vl = sqrt(D.squared_length());
            // std::cout << "vl-(r+rc): " << vl - (fabs(rsdf)+fabs(res(contact_point_indices[i],3))) << std::endl;
            D /= vl;
            // std::cout << "D.squared_lenght(): " << D.squared_length() << std::endl;
            Point  P = Csdf + fabs(rsdf)*D;
            Vector N = (rsdf >= 0)? -D: D;
            Pwns.emplace_back(P,N);
        }
    }

}

int main(int argc, char** argv) {
    const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("data/3D/spheres.csv");;
    std::ifstream in(filename);
    // std::ifstream in("data/3D/koala.obj_gridvals_30_regular.csv");

    CGAL::Timer timer;
    std::vector<Eigen::RowVector4d> input_spheres;

    bool filter_contact_spheres_bbx=true;

    double x, y, z, r;
    while(in >> x){

        in.ignore(10,','); in >> y;  in.ignore(10,','); in >> z; in.ignore(10,','); in >> r;
        // I only filter out the positive spheres here as a first demo. The proper version uses positive / negative spheres separately.
        // TODO: implement the proper version
        input_spheres.emplace_back(Eigen::RowVector4d(x,y,z,r));
        // std::cout << "Sphere: " << x << y << z << r << std::endl;
    }

    std::cout << "Read " << input_spheres.size() << " spheres" << std::endl;

    Eigen::MatrixXd G(input_spheres.size(), 4);
    for (int i=0; i<input_spheres.size(); i++) G.row(i) = input_spheres[i];

    Eigen::MatrixXd bbxl(2,3);
    bbxl.row(0) = G.block(0,0,G.rows(),3).colwise().minCoeff();
    bbxl.row(1) = G.block(0,0,G.rows(),3).colwise().maxCoeff();

    // std::cout << "... pos/neg subsets: " << std::endl;

    int nnrc = (G.col(3).array() >= 0.).array().count();
    // std::cout << ( G.col(2).array() >= 0) << std::endl;
    // std::cout << "nnrc: " << nnrc << ", nrc: " << G.rows()-nnrc << std::endl;
    Eigen::MatrixXd Gp(nnrc,4);
    Eigen::MatrixXd Gn(G.rows()-nnrc,4);
    int np=0;
    int nn=0;
    for (int i=0; i<input_spheres.size(); i++){
        if (input_spheres[i](3) >= 0){
            Gp.row(np++) = input_spheres[i];
        } else {
            Gn.row(nn++) = input_spheres[i];
        }
    }

    std::cout << "... main calculation... " << std::endl;
    timer.start();
    std::vector<Point_with_normal> Pwns;
    contact_points(Gp, Pwns, (filter_contact_spheres_bbx)?&bbxl:NULL);
    contact_points(Gn, Pwns, (filter_contact_spheres_bbx)?&bbxl:NULL);

    std::cout << "Computing " << Pwns.size() << " contact points with normals in " << timer.time() << " sec." << std::endl;
    std::cout << "PSR from " << Pwns.size() << " points with normals" << std::endl;
    timer.reset();

    Mesh output_mesh;

    double average_spacing = CGAL::compute_average_spacing<CGAL::Sequential_tag>
      (Pwns, 6, CGAL::parameters::point_map(CGAL::First_of_pair_property_map<Point_with_normal>()));

    std::cout << "Average spacing: " << average_spacing << std::endl;

    // These are the defaults
    double sm_angle = 20.0;
    double sm_radius = 30.0;
    double sm_distance = 0.375;


    if (CGAL::poisson_surface_reconstruction_delaunay
      (Pwns.begin(), Pwns.end(),
       CGAL::First_of_pair_property_map<Point_with_normal>(),
       CGAL::Second_of_pair_property_map<Point_with_normal>(),
       output_mesh, average_spacing, sm_angle, sm_radius, sm_distance))
    {
       CGAL::IO::write_polygon_mesh("tmp.off", output_mesh,
                                    CGAL::parameters::stream_precision(17));
    } else {
       std::cerr << "Poisson reconstruction failed." << std::endl;
       return EXIT_FAILURE;
    }
    std::cout << "Poisson reconstruction done in " << timer.time() << " sec." << std::endl;
    return 0;
}
