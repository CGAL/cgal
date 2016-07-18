#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Timer.h>


#include <CGAL/boost/graph/property_maps.h>

#include <fstream>
#include <ostream>

// typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Simple_cartesian<double> K;
typedef CGAL::Surface_mesh<K::Point_3> Mesh;


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

  std::vector<K::Triangle_3> t1;
  t1.reserve(num_faces(m1));
  CGAL::Triangle_from_face_descriptor_map<Mesh> map1(&m1);
  BOOST_FOREACH(Mesh::Face_index f, faces(m1))
    t1.push_back(get(map1,f));

  std::vector<K::Triangle_3> t2;
  t2.reserve(num_faces(m2));
  CGAL::Triangle_from_face_descriptor_map<Mesh> map2(&m2);
  BOOST_FOREACH(Mesh::Face_index f, faces(m2))
    t2.push_back(get(map2,f));


  CGAL::Timer time;
  time.start();
  std::cout << "Distance between meshes (parallel)"
            << CGAL::Polygon_mesh_processing::approximated_Hausdorff_distance<CGAL::Parallel_tag, K>(t1,t2,0.001)
            << "\n";
  time.stop();
  std::cout << "done in " << time.time() << "s.\n";

  time.reset();
  time.start();
  std::cout << "Distance between meshes (sequential)"
            << CGAL::Polygon_mesh_processing::approximated_Hausdorff_distance<CGAL::Sequential_tag, K>(t1,t2,0.001)
            << "\n";
  time.stop();
  std::cout << "done in " << time.time() << "s.\n";
}


