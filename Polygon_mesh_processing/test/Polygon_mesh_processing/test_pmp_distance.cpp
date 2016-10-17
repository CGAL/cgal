#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Timer.h>


#include <CGAL/boost/graph/property_maps.h>

#include <CGAL/Polygon_mesh_processing/distance.h>

#include <fstream>
#include <ostream>

// typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Simple_cartesian<double> K;
typedef CGAL::Surface_mesh<K::Point_3> Mesh;

bool general_tests(const Mesh& m1,
                   const Mesh& m2 )
{

    std::cout << "Symetric distance between meshes (sequential) "
              << CGAL::Polygon_mesh_processing::approximated_symmetric_Hausdorff_distance<CGAL::Sequential_tag>(m1,m2,4000)
              << "\n";

    std::vector<K::Point_3> points;

    points.push_back(K::Point_3(0,0,0));
    points.push_back(K::Point_3(0,1,0));
    points.push_back(K::Point_3(1,0,0));
    points.push_back(K::Point_3(0,0,1));
    points.push_back(K::Point_3(0,2,0));

    Mesh m;
    CGAL::make_tetrahedron(points[0],
                           points[1],
                           points[2],
                           points[3],
                           m);

    std::cout << "Max distance to point set "
              << CGAL::Polygon_mesh_processing::max_distance_to_point_set(m,points,4000)
              << "\n";

    std::cout << "Max distance to triangle mesh (sequential) "
              << CGAL::Polygon_mesh_processing::max_distance_to_triangle_mesh<CGAL::Sequential_tag>(points,m)
              << "\n";

}

int main(int argc, char** argv)
{
  Mesh m1,m2;

  std::ifstream input(argv[1]);
  input >> m1;
  input.close();

  input.open(argv[2]);
  input >> m2;

  std::cout << "First mesh has " << num_faces(m1) << " faces\n";
  std::cout << "Second mesh has " << num_faces(m2) << " faces\n";

  CGAL::Timer time;
  time.start();
  std::cout << "Distance between meshes (parallel) "
            << CGAL::Polygon_mesh_processing::approximated_Hausdorff_distance<CGAL::Parallel_tag, K>(m1,m2,0.001, get(CGAL::vertex_point, m1), get(CGAL::vertex_point, m2), CGAL::Polygon_mesh_processing::GRID)
            << "\n";
  time.stop();
  std::cout << "done in " << time.time() << "s.\n";

  time.reset();
  time.start();
  std::cout << "Distance between meshes (sequential) "
            << CGAL::Polygon_mesh_processing::approximated_Hausdorff_distance<CGAL::Sequential_tag, K>(m1,m2,0.001, get(CGAL::vertex_point, m1), get(CGAL::vertex_point, m2), CGAL::Polygon_mesh_processing::GRID)
            << "\n";
  time.stop();
  std::cout << "done in " << time.time() << "s.\n";

  general_tests(m1,m2);
}


