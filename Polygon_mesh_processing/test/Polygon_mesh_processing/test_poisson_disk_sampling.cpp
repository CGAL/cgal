#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/distance.h>
#ifdef CGAL_USE_BSURF
#include <CGAL/Polygon_mesh_processing/Bsurf/locally_shortest_path.h>
#endif

#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/stitch_borders.h>
#include <iostream>
#include <string>
#include <vector>

#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits_3.h>

#include <CGAL/Polygon_mesh_processing/remesh.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel   K;
typedef K::Point_3                                            Point;

typedef CGAL::Surface_mesh<Point>                             Mesh;
typedef boost::graph_traits<Mesh>::vertex_descriptor          vertex_descriptor;
typedef boost::graph_traits<Mesh>::face_descriptor            face_descriptor;

typedef CGAL::AABB_face_graph_triangle_primitive<Mesh>        AABB_face_graph_primitive;
typedef CGAL::AABB_traits_3<K, AABB_face_graph_primitive>     AABB_face_graph_traits;

namespace PMP = CGAL::Polygon_mesh_processing;

typedef PMP::Face_location<Mesh, double>                      Face_location;
#ifdef CGAL_USE_BSURF
typedef PMP::Edge_location<Mesh, double>                      Edge_location;
#endif

#ifdef CGAL_USE_BSURF
double geodesiceApproximation(const Face_location source, const Face_location target, Mesh mesh)
{
  PMP::Dual_geodesic_solver<double> solver;
  PMP::init_geodesic_dual_solver(solver, mesh);
  std::vector<Edge_location> edge_locations;

  CGAL::Polygon_mesh_processing::locally_shortest_path<double>(source, target, mesh, edge_locations, solver);

  return PMP::path_length<K>(edge_locations,source,target,mesh);
}
#endif

int main(int argc, char* argv[])
{
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/eight.off");

  Mesh mesh;
  if(!PMP::IO::read_polygon_mesh(filename, mesh))
  {
    std::cerr << "Invalid input." << std::endl;
    return 1;
  }

  const double sampling_radius = (argc > 2) ? std::atof(argv[2]) : 0.07;

  const std::size_t number_of_darts = (argc > 3) ? std::atof(argv[3]) : 30;

  CGAL::Random random_seed = (argc > 4) ? std::atof(argv[4]) : CGAL::get_default_random();

  std::vector<Point> euclidean_points;
  PMP::sample_triangle_mesh(mesh,
                            std::back_inserter(euclidean_points),
                            CGAL::parameters::use_poisson_disk_sampling_euclidean(true).
                                              sampling_radius(sampling_radius).number_of_darts(number_of_darts).random_seed(random_seed));

  std::ofstream out1("euclidean-sampling.xyz");
  out1 << std::setprecision(17);
  std::copy(euclidean_points.begin(), euclidean_points.end(), std::ostream_iterator<Point>(out1, "\n"));

  std::cout << euclidean_points.size() << std::endl;

 //test if euclidean_points are far enough apart.

  bool result=true;
  double d;
  double c;

  for(std::size_t i=0; i<euclidean_points.size(); i++){
      c = 100;
      for(std::size_t j = 0; j<euclidean_points.size(); j++){

        d = sqrt(CGAL::squared_distance(euclidean_points[i],euclidean_points[j]));
        //no point too close
        if(d< sampling_radius && d != 0){
            result = false;
        }
        // no point too far
        if(d < c && d !=0){
            c = d;
        }
      }
      if(c  > 2*sampling_radius){
          result = false;
      }
  }
  std::cout << "The Euclidean test result is: " << result << std::endl;

  Mesh square;
  Mesh::Halfedge_index h = CGAL::make_quad(Point(0,0,0), Point(1,0,0),
                                           Point(1,1,0),Point(0,1,0),
                                           square);
  CGAL::Euler::split_face(h, next(next(h,square), square), square);
  PMP::isotropic_remeshing(faces(square), 0.05, square);
  std::cout << "Square now has " << faces(square).size() << " faces" << std::endl;
  PMP::stitch_borders(square);
  CGAL::IO::write_polygon_mesh("square.off", square, CGAL::parameters::stream_precision(17));
  std::vector<Point> square_points;

  PMP::sample_triangle_mesh(square,
                              std::back_inserter(square_points),
                              CGAL::parameters::use_poisson_disk_sampling_euclidean(true).
                                                sampling_radius(sampling_radius).number_of_darts(number_of_darts).random_seed(random_seed));

  std::cout << "There are: " << square_points.size() << " points in the square." << std::endl;

  std::ofstream out("square.xyz");
  out << std::setprecision(17);
  std::copy(square_points.begin(), square_points.end(), std::ostream_iterator<Point>(out, "\n"));

 /*
  //

#ifdef TEST_PURPOSE
  CGAL::AABB_tree<AABB_face_graph_traits> tree;
  PMP::build_AABB_tree(mesh, tree);
  Face_location query_location_target = PMP::locate_with_AABB_tree(euclidean_points[0], tree, mesh);
  Face_location query_location_source = PMP::locate_with_AABB_tree(euclidean_points[1], tree, mesh);

  std::cout << "Euclidean: " << sqrt(CGAL::squared_distance(euclidean_points[0],euclidean_points[1])) << std::endl;
  std::cout << "Geodesic: " << geodesiceApproximation(query_location_source, query_location_target, mesh) << std::endl;

    //Test geodesic
  //why can't I use geodesic?????

  result=true;
  std::vector<Point> geodesic_points;
  PMP::sample_triangle_mesh(mesh,
                              std::back_inserter(geodesic_points),
                              CGAL::parameters::use_poisson_disk_sampling_geodesic(true).
                                                sampling_radius(sampling_radius).number_of_darts(number_of_darts).random_seed(random_seed));
  std::ofstream out("geodesic-sampling.xyz");
  out << std::setprecision(17);
  std::copy(geodesic_points.begin(), geodesic_points.end(), std::ostream_iterator<Point>(out, "\n"));
  std::cout << geodesic_points.size() << std::endl;

  for(std::size_t i=0; i<geodesic_points.size(); i++){
    c = 100;
    Face_location query_location_source = PMP::locate_with_AABB_tree(geodesic_points[i], tree, mesh);
    for(std::size_t j = 0; j<geodesic_points.size(); j++){
        Face_location query_location_target = PMP::locate_with_AABB_tree(geodesic_points[j], tree, mesh);
        d = geodesiceApproximation(query_location_source, query_location_target, mesh);
        //no point too close
        if(d< sampling_radius && d != 0){
            result = false;
        }
        // no point too far
        if(d < c && d !=0){
            c = d;
        }
    }
    if(c  > 2*sampling_radius){
        result = false;
    }
}
    std::cout << "The geodesic test result is: " << result << std::endl;
#endif


 */
  return 0;
}
