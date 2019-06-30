#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
//#include <CGAL/Exact_integer.h>
//#include <CGAL/Extended_homogeneous.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/draw_nef_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/convert_nef_polyhedron_to_polygon_mesh.h>

#include <fstream>
#include <iostream>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
//typedef CGAL::Exact_integer NT;
//typedef CGAL::Extended_homogeneous<NT> Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> Surface_mesh;
typedef CGAL::Nef_polyhedron_3<Kernel> Nef_polyhedron;
typedef Nef_polyhedron::Plane_3 Plane_3;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;

int main(int argc, char *argv[])
{
  Polyhedron P;
  std::ifstream ifs((argc > 1) ? argv[1] : "data/cross.off");
  ifs >> P;
  std::cout << "============== Polyhedron =============" << std::endl
            << std::endl;
  std::cout << "Number of vertices: " << P.size_of_vertices() << std::endl;
  std::cout << "Number of halfedges: " << P.size_of_halfedges() << std::endl;
  std::cout << "Number of facets: " << P.size_of_facets() << std::endl
            << std::endl;

  Nef_polyhedron N(P);

  std::cout << "================= Nef ==================" << std::endl
            << std::endl;
  std::cout << "Number of volumes: " << N.number_of_volumes() << std::endl;
  std::cout << "Number of vertices: " << N.number_of_vertices() << std::endl;
  std::cout << "Number of edges: " << N.number_of_edges() << std::endl;
  std::cout << "Number of facets: " << N.number_of_facets() << std::endl
            << std::endl;
  std::cout << "Number of halfedges: " << N.number_of_halfedges() << std::endl;
  std::cout << "Number of halffacets: " << N.number_of_halffacets() << std::endl
            << std::endl;
  std::cout << "Number of sfaces: " << N.number_of_sfaces() << std::endl;
  std::cout << "Number of shalfedges : " << N.number_of_shalfedges()
            << std::endl;
  std::cout << "Number of shalfloops: " << N.number_of_shalfloops() << std::endl
            << std::endl;

  // output the result into a Surface_mesh
  Surface_mesh output;
  CGAL::convert_nef_polyhedron_to_polygon_mesh(N, output);
  std::ofstream out;
  out.open("out.off");
  out << output;
  out.close();

  CGAL::draw(N);
}
