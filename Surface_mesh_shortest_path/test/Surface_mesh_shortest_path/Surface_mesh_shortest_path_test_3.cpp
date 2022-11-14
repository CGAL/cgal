#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>

#include <CGAL/Surface_mesh_shortest_path/Surface_mesh_shortest_path_traits.h>
#include <CGAL/Surface_mesh_shortest_path/Surface_mesh_shortest_path.h>
#include <CGAL/Surface_mesh_shortest_path/internal/misc_functions.h>

#include <CGAL/boost/graph/iterator.h>

#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_tree.h>

#include <CGAL/Random.h>

#include <CGAL/test_util.h>
#include <iostream>
#include <fstream>
#include <utility>
#include <cmath>

#include "check.h"

int main(int argc, char* argv[])
{
  typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
  typedef CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3> Polyhedron_3;
  typedef CGAL::Surface_mesh_shortest_path_traits<Kernel, Polyhedron_3> Traits;
  typedef Traits::Barycentric_coordinates Barycentric_coordinates;
  typedef Traits::FT FT;
  typedef Traits::Point_3 Point_3;
  typedef Traits::Triangle_3 Triangle_3;
  typedef boost::graph_traits<Polyhedron_3> Graph_traits;
  typedef Graph_traits::face_descriptor face_descriptor;
  typedef Graph_traits::face_iterator face_iterator;
  typedef CGAL::Surface_mesh_shortest_path<Traits> Surface_mesh_shortest_path;
  typedef boost::property_map<Polyhedron_3, CGAL::vertex_point_t>::const_type VPM;
  typedef boost::property_map<Polyhedron_3, CGAL::face_index_t>::const_type FIM;

  typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron_3, VPM> AABB_face_graph_primitive;
  typedef CGAL::AABB_traits<Kernel, AABB_face_graph_primitive> AABB_face_graph_traits;

  Traits traits;

  Traits::Construct_barycenter_3 construct_barycenter_3(traits.construct_barycenter_3_object());

  std::string mesh(argv[1]);

  int randSeed = 6008991;

  if (argc > 2)
  {
    randSeed = std::atoi(argv[2]);
  }

  CGAL::Random random(randSeed);

  Polyhedron_3 polyhedron;

  std::ifstream in(mesh.c_str());

  in >> polyhedron;

  in.close();

  CGAL::set_halfedgeds_items_id(polyhedron);

  Surface_mesh_shortest_path shortestPaths(polyhedron, traits);

  face_iterator facesBegin, facesEnd;
  boost::tie(facesBegin, facesEnd) = faces(polyhedron);

  std::vector<face_descriptor> facesList;

  for (face_iterator facesCurrent = facesBegin; facesCurrent != facesEnd; ++facesCurrent)
  {
    facesList.push_back(*facesCurrent);
  }

  size_t numTrials = 30;

  FIM faceIndexMap(get(boost::face_index, polyhedron));
  VPM vertexPointMap(get(CGAL::vertex_point, polyhedron));

  for (size_t i = 0; i < numTrials; ++i)
  {
    size_t faceIndex = random.get_int(0, static_cast<int>(facesList.size()));
    face_descriptor face = facesList[faceIndex];

    Triangle_3 faceTriangle = CGAL::Surface_mesh_shortest_paths_3::internal::triangle_from_halfedge<Triangle_3, Polyhedron_3, VPM>(halfedge(face, polyhedron), polyhedron, vertexPointMap);

    Barycentric_coordinates location = CGAL::test::random_coordinates<Traits>(random);

    Point_3 location3d = construct_barycenter_3(faceTriangle[0], location[0], faceTriangle[1], location[1], faceTriangle[2], location[2]);

    Surface_mesh_shortest_path::Face_location faceLocation = shortestPaths.locate<AABB_face_graph_traits>(location3d);

    CHECK_EQUAL(faceIndexMap[face], faceIndexMap[faceLocation.first]);
    CHECK_CLOSE(location[0], faceLocation.second[0], FT(0.0001));
    CHECK_CLOSE(location[1], faceLocation.second[1], FT(0.0001));
    CHECK_CLOSE(location[2], faceLocation.second[2], FT(0.0001));
  }

  return 0;
}
