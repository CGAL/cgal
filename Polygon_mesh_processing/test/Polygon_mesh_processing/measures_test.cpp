#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/properties_Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/measure.h>

#include <iostream>
#include <fstream>
#include <iterator>
#include <list>

#include <boost/foreach.hpp>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point;

typedef CGAL::Polyhedron_3<K>       Polyhedron;
typedef CGAL::Surface_mesh<Point>   Surface_mesh;

namespace PMP = CGAL::Polygon_mesh_processing;

template<typename Mesh>
void test(const Mesh& pmesh)
{
  typedef boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;
  typedef boost::graph_traits<Mesh>::face_descriptor     face_descriptor;

  halfedge_descriptor border_he;
  BOOST_FOREACH(halfedge_descriptor h, halfedges(pmesh))
  {
    if (is_border(h, pmesh))
    {
      border_he = h;
      break;
    }
  }
  double border_l = PMP::border_length(border_he, pmesh);
  std::cout << "length of hole border = " << border_l << std::endl;

  std::list<face_descriptor> patch;
  BOOST_FOREACH(halfedge_descriptor h, halfedges(pmesh))
  {
    if (is_border(h, pmesh) || is_border(opposite(h, pmesh), pmesh))
      continue;
    else
    {
      double face_area = PMP::area(face(h, pmesh), pmesh);
      std::cout << "face area = " << face_area << std::endl;
 
      patch.push_back(face(h, pmesh));
      patch.push_back(face(opposite(h, pmesh), pmesh));
      patch.push_back(face(opposite(next(h, pmesh), pmesh), pmesh));
      patch.push_back(face(opposite(prev(h, pmesh), pmesh), pmesh));
      break;
    }
  }
  double patch_area = PMP::area(patch, pmesh);
  std::cout << "patch area = " << patch_area << std::endl;

  double mesh_area = PMP::area(pmesh);
  std::cout << "mesh area = " << mesh_area << std::endl;
}

void test_polyhedron(const char* filename)
{
  //run test for a Polyhedron
  Polyhedron poly; // file should contain oriented polyhedron
  std::ifstream input(filename);

  if (!input || !(input >> poly))
  {
    std::cerr << "Error: cannot read Polyhedron : " << filename << "\n";
    assert(!poly.empty());
    assert(false);
    return;
  }

  test(poly);
}

void test_surface_mesh(const char* filename)
{
  Surface_mesh sm;
  std::ifstream input(filename);

  if (!input || !(input >> sm))
  {
    std::cerr << "Error: cannot read Surface mesh : " << filename << "\n";
    assert(sm.number_of_vertices() > 0);
    assert(false);
    return;
  }

  test(sm);
}

int main(int argc, char* argv[])
{
  const char* filename = (argc > 1) ? argv[1] : "data/mech-holes-shark.off";

  test_polyhedron(filename);
  test_surface_mesh(filename);

  std::cerr << "All done." << std::endl;
  return 0;
}
