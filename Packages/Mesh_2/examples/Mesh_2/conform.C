// file examples/Mesh_2/conform.C
#include <CGAL/basic.h>
#include <iostream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Conforming_Delaunay_triangulation_2.h>

struct K : public CGAL::Exact_predicates_inexact_constructions_kernel {};
typedef CGAL::Constrained_Delaunay_triangulation_2<K> CDT;
typedef CDT::Vertex_handle Vertex_handle;
typedef CDT::Point Point;

int main(int, char**)
{
  CDT cdt;

  Vertex_handle va = cdt.insert(Point(-4,0));
  Vertex_handle vb = cdt.insert(Point(0,-1));
  Vertex_handle vc = cdt.insert(Point(4,0));
  Vertex_handle vd = cdt.insert(Point(0,1));
  cdt.insert(Point(2, 0.6));

  cdt.insert_constraint(va, vb);
  cdt.insert_constraint(vb, vc);
  cdt.insert_constraint(vc, vd);
  cdt.insert_constraint(vd, va);

  std::cout << "Number of vertices: " << cdt.number_of_vertices() 
	    << std::endl;

  std::cout << "Making conforming Delaunay..." << std::endl;
  CGAL::make_conforming_Delaunay_2(cdt);

  std::cout << "Number of vertices: " << cdt.number_of_vertices() 
	    << std::endl;

  std::cout << "Making conforming Gabriel..." << std::endl;
  CGAL::make_conforming_Gabriel_2(cdt);

  std::cout << "Number of vertices: " << cdt.number_of_vertices() 
	    << std::endl;
}
