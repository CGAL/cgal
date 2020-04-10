#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <iostream>
#include <fstream>
#include <CGAL/IO/OFF_reader.h>
#include <CGAL/boost/graph/Euler_operations.h>

#include <CGAL/Real_timer.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> Mesh;

int main(int argc, char** argv)
{
  {
  std::cout << "Reading from stream\n";
  Mesh m;
  CGAL::Real_timer timer;
  timer.start();
  std::ifstream in((argc>1) ? argv[1] : "data/genus3.off");
  in >> m;
  timer.stop();
  std::cout << "  is_valid? " << CGAL::is_valid_polygon_mesh(m) << "\n";
  std::cout << "Total time: " << timer.time() << std::endl << std::endl;
  }
////////////////////////////////
  {
  std::cout << "Reading from soup + iterative add_face\n";
  Mesh m;
  CGAL::Real_timer timer;
  timer.start();
  std::ifstream in((argc>1) ? argv[1] : "data/blobby.off");
  std::vector<Kernel::Point_3> points;
  std::vector<std::array<Mesh::Vertex_index, 3> > triangles;
  CGAL::read_OFF(in, points, triangles);
  std::cout << "  Read soup: " << timer.time() << std::endl;
  m.reserve(static_cast<unsigned int>(points.size()),
            static_cast<unsigned int>(3*triangles.size()/2),
            static_cast<unsigned int>(triangles.size()));
  for (const Kernel::Point_3& pt : points)
    m.add_vertex(pt);
  CGAL::Real_timer subtimer;
  subtimer.start();
  for (const std::array<Mesh::Vertex_index, 3>& t : triangles)
    CGAL::Euler::add_face(t, m);
  subtimer.stop();
  timer.stop();
  std::cout << "  is_valid? " << CGAL::is_valid_polygon_mesh(m) << "\n";
  std::cout << "  time for iterative add_face: " << subtimer.time() << std::endl;
  std::cout << "Total time: " << timer.time() << std::endl << std::endl;
  }
////////////////////////////////
  {
  std::cout << "Reading from soup + add_faces\n";
  Mesh m;
  CGAL::Real_timer timer;
  timer.start();
  std::ifstream in((argc>1) ? argv[1] : "data/blobby.off");
  std::vector<Kernel::Point_3> points;
  std::vector<std::array<Mesh::Vertex_index, 3> > triangles;
  CGAL::read_OFF(in, points, triangles);
  std::cout << "  Read soup: " << timer.time() << std::endl;
  m.reserve(static_cast<unsigned int>(points.size()),
            static_cast<unsigned int>(3*triangles.size()/2),
            static_cast<unsigned int>(triangles.size()));
  for (const Kernel::Point_3& pt : points)
    m.add_vertex(pt);
  CGAL::Real_timer subtimer;
  subtimer.start();
  CGAL::Euler::add_faces(triangles, m);
  subtimer.stop();
  timer.stop();
  std::cout << "  is_valid? " << CGAL::is_valid_polygon_mesh(m) << "\n";
  std::cout << "  time for add_faces: " << subtimer.time() << std::endl;
  std::cout << "Total time: " << timer.time() << std::endl;
  }

  return 0;
}
