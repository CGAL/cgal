#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/Real_timer.h>
#include <CGAL/IO/polygon_soup_io.h>

#include <iostream>
#include <fstream>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> Mesh;

void convert_to_vertex_triples(
  const std::vector<std::array<std::size_t, 3> >& faces_ids,
        std::vector<std::array<Mesh::Vertex_index, 3> >& triangles)
{
  triangles.reserve(faces_ids.size());
  for (const std::array<std::size_t, 3>& a : faces_ids)
    triangles.push_back(
      CGAL::make_array( Mesh::Vertex_index(static_cast<Mesh::size_type>(a[0])),
                        Mesh::Vertex_index(static_cast<Mesh::size_type>(a[1])),
                        Mesh::Vertex_index(static_cast<Mesh::size_type>(a[2])) ) );
}
int main(int argc, char** argv)
{
  {
  std::cout << "Reading from stream\n";
  CGAL::Real_timer timer;
  timer.start();

  Mesh m;
  const char* filename = (argc>1) ? argv[1] : "data/genus3.off";
  CGAL::IO::read_polygon_mesh(filename, m);

  std::cout << "  is_valid? " << CGAL::is_valid_polygon_mesh(m) << "\n";
  std::cout << "Total time: " << timer.time() << std::endl << std::endl;
  }

////////////////////////////////

  {
  std::cout << "Reading from soup + iterative add_face\n";

  CGAL::Real_timer timer;
  timer.start();

  const char* filename = (argc>1) ? argv[1] : "data/blobby.off";
  std::vector<Kernel::Point_3> points;
  std::vector<std::array<std::size_t, 3> > faces_ids;
  CGAL::IO::read_polygon_soup(filename, points, faces_ids);
  std::cout << "  Read soup: " << timer.time() << std::endl;

  std::vector<std::array<Mesh::Vertex_index, 3> > triangles;
  convert_to_vertex_triples(faces_ids, triangles);

  Mesh m;
  m.reserve(static_cast<Mesh::size_type>(points.size()),
            static_cast<Mesh::size_type>(3*triangles.size()/2),
            static_cast<Mesh::size_type>(triangles.size()));
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

  CGAL::Real_timer timer;
  timer.start();

  const char* filename = (argc>1) ? argv[1] : "data/blobby.off";

  std::vector<Kernel::Point_3> points;
  std::vector<std::array<std::size_t, 3> > faces_ids;
  CGAL::IO::read_polygon_soup(filename, points, faces_ids);
  std::cout << "  Read soup: " << timer.time() << std::endl;

  std::vector<std::array<Mesh::Vertex_index, 3> > triangles;
  convert_to_vertex_triples(faces_ids, triangles);

  Mesh m;
  m.reserve(static_cast<Mesh::size_type>(points.size()),
            static_cast<Mesh::size_type>(3*triangles.size()/2),
            static_cast<Mesh::size_type>(triangles.size()));
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
