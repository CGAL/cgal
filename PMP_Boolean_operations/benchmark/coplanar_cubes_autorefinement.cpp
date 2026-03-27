#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/autorefinement.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/transform.h>
#include <CGAL/Real_timer.h>

#include <boost/container/small_vector.hpp>

#include <iostream>
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel K;
typedef CGAL::Surface_mesh<K::Point_3>                      Mesh;
typedef K::Point_3 Point_3;

typedef CGAL::Simple_cartesian<double> Cartesian;
typedef Cartesian::Point_3 Double_Point_3;

namespace PMP = CGAL::Polygon_mesh_processing;
namespace params = CGAL::parameters;

struct Sphere_function {
    double radius;
    Sphere_function(double r) : radius(r) {}
    Kernel::FT operator()(const Kernel::Point_3& p) const {
        return p.x()*p.x() + p.y()*p.y() + p.z()*p.z() - radius*radius;
    }
};

//Thanks Roberto!
K::Aff_transformation_3
random_rotation(CGAL::Random &gen)
{
  double a=gen.get_double(0,2*CGAL_PI);

  double ca = cos(a);
  double sa = sin(a);

  K::Aff_transformation_3 aff(1, 0, 0,
                              0, ca,-sa,
                              0, sa, ca);
  std::cout << "Rotation by " << a << std::endl;
  return aff;
}

int main(int argc, char** argv)
{
  Mesh cube;
  std::vector<Point_3> points;
  std::vector<boost::container::small_vector<std::size_t, 3>> faces;

  CGAL::make_hexahedron(
      Point_3(-1,-1,-1),
      Point_3(1,-1,-1),
      Point_3(1,1,-1),
      Point_3(-1,1,-1),
      Point_3(-1,1,1),
      Point_3(-1,-1,1),
      Point_3(1,-1,1),
      Point_3(1,1,1),
      cube,
      CGAL::parameters::do_not_triangulate_faces(false)
      );

  std::cout << "Iterative intersection of rotative cubes with snapping" << std::endl;

  int i=0;
  CGAL::Random random_gen = argc == 1 ? CGAL::get_default_random() : CGAL::Random(std::stoi(argv[1]));

  Mesh inter=cube;
  while(true)
  {
    std::cout << "Iteration " << i << std::endl;

    CGAL::Real_timer t;
    t.start();

    std::cout << "Add a randomly rotated cube to the scene" << std::endl;
    Mesh rotated_cube=cube;
    PMP::transform(random_rotation(random_gen), rotated_cube);

    std::cout << "compute_intersection" << std::endl;
    bool OK = PMP::corefine_and_compute_intersection(inter, rotated_cube, inter);

    if(!OK){
      std::cout << "No manifold, stop experiment" << std::endl;
      exit(0);
    }

    points.clear();
    faces.clear();
    PMP::polygon_mesh_to_polygon_soup(inter, points, faces);

    std::cout << "Snapped the points on double" << std::endl;
    bool success=PMP::autorefine_triangle_soup(points, faces, CGAL::parameters::apply_iterative_snap_rounding(true).erase_all_duplicates(true).concurrency_tag(CGAL::Parallel_if_available_tag()));
    t.stop();
    if(!success){
      std::cout << "Rounding failed" << std::endl;
      exit(0);
    }

    //dump model every 100 iterations
    if(i%100==0){
      std::cout << "Dump model" << std::endl;
      std::vector<Double_Point_3> double_points;
      for(auto &p: points)
        double_points.emplace_back(CGAL::to_double(p.x()),CGAL::to_double(p.y()),CGAL::to_double(p.z()));
      std::ofstream outfile("cubes_"+std::to_string(i)+".off");
      outfile.precision(17);
      outfile << "OFF\n" << points.size() << " " << faces.size() << " 0\n";
      for(auto p : points)
      outfile << p.x() << " " << p.y() << " " << p.z() << std::endl;

      for(auto &t : faces)
          outfile << "3" << " " << t[0] << " " << t[1] << " " << t[2] << std::endl;
      outfile.close();//
    }

    std::cout << "#points = " << points.size() << " and #triangles = " << faces.size() << " in " << t.time() << " sec.\n\n" << std::endl;

    inter.clear();
    PMP::polygon_soup_to_polygon_mesh(points, faces, inter);
    ++i;
  }

  return 0;
}

