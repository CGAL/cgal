#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Maximal_empty_spheres/contact_points_from_signed_distances.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <filesystem>
#include <string>

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
typedef Kernel::Sphere_3 Sphere;
typedef std::pair<Point, Vector> Point_with_normal;
typedef CGAL::First_of_pair_property_map<Point_with_normal> Point_map;
typedef CGAL::Second_of_pair_property_map<Point_with_normal> Normal_map;
typedef CGAL::Surface_mesh<Point> Mesh;



int main(int argc, char** argv) {

    bool filter_contact_spheres_with_bbox = true;
    int debug_level = 0;

    const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("data/3D/spheres.csv");;
    std::ifstream in(filename);

    CGAL::Timer timer;
    std::vector<std::pair<Sphere,int>> input_spheres;
    std::vector<Point_with_normal> pwns;

    double x, y, z, r;
    while(in >> x){
        in.ignore(10,','); in >> y;  in.ignore(10,','); in >> z; in.ignore(10,','); in >> r;
        int inside = r < 0 ? -1 : 1; // positive radius means outside, negative radius means inside
        input_spheres.emplace_back(std::make_pair(Sphere(Point(x,y,z),CGAL::square(r)),inside));
    }
    std::cout << "Read " << input_spheres.size() << " spheres" << std::endl;

   timer.start();
    CGAL::contact_points_from_signed_distances(input_spheres, std::back_inserter(pwns), filter_contact_spheres_with_bbox, debug_level);

    std::cout << "Computed " << pwns.size() << " contact points with normals in " << timer.time() << " sec." << std::endl;

    Mesh output_mesh;
    double average_spacing = CGAL::compute_average_spacing<CGAL::Sequential_tag>
      (pwns, 6, CGAL::parameters::point_map(CGAL::First_of_pair_property_map<Point_with_normal>()));

    // These are the defaults
    double sm_angle = 20.0;
    double sm_radius = 30.0;
    double sm_distance = 0.375;


    {
      std::string output_filename = std::filesystem::path(filename).filename().replace_extension(".pwn").string();
      std::cout << "Writing contact points to " << output_filename << std::endl;
      std::ofstream out(output_filename);
      for(const auto& pwn : pwns) {
        out << pwn.first << " " << pwn.second << std::endl;
      }
    }

    std::cout << "Poisson reconstruction with average spacing: " << average_spacing << std::endl;

    timer.reset();
    if (CGAL::poisson_surface_reconstruction_delaunay
      (pwns.begin(), pwns.end(),
       CGAL::First_of_pair_property_map<Point_with_normal>(),
       CGAL::Second_of_pair_property_map<Point_with_normal>(),
       output_mesh, average_spacing, sm_angle, sm_radius, sm_distance))
    {
       timer.stop();
       std::cout << "done in " << timer.time() << " sec." << std::endl;
       std::string output_filename = std::filesystem::path(filename).filename().replace_extension(".off").string();
       std::cout << "writing output mesh to " << output_filename << std::endl;
       CGAL::IO::write_polygon_mesh(output_filename, output_mesh,
                                    CGAL::parameters::stream_precision(17));
    } else {
       std::cerr << "Poisson reconstruction failed." << std::endl;
       return EXIT_FAILURE;
    }
    return 0;
}
