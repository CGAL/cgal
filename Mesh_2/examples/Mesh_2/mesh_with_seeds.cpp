#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>

#include <iostream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_2<K> Vb;
typedef CGAL::Delaunay_mesh_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds> CDT;
typedef CGAL::Delaunay_mesh_size_criteria_2<CDT> Criteria;

typedef CDT::Vertex_handle Vertex_handle;
typedef CDT::Point Point;

int main()
{
  CDT cdt;
  Vertex_handle va = cdt.insert(Point(2,0));
  Vertex_handle vb = cdt.insert(Point(0,2));
  Vertex_handle vc = cdt.insert(Point(-2,0));
  Vertex_handle vd = cdt.insert(Point(0,-2));

  cdt.insert_constraint(va, vb);
  cdt.insert_constraint(vb, vc);
  cdt.insert_constraint(vc, vd);
  cdt.insert_constraint(vd, va);

  va = cdt.insert(Point(3,3));
  vb = cdt.insert(Point(-3,3));
  vc = cdt.insert(Point(-3,-3));
  vd = cdt.insert(Point(3,0-3));

  cdt.insert_constraint(va, vb);
  cdt.insert_constraint(vb, vc);
  cdt.insert_constraint(vc, vd);
  cdt.insert_constraint(vd, va);

  std::list<Point> list_of_seeds;

  list_of_seeds.push_back(Point(0, 0));

  std::cout << "Number of vertices: " << cdt.number_of_vertices() << std::endl;

  std::cout << "Meshing the domain..." << std::endl;
  CGAL::refine_Delaunay_mesh_2(cdt, list_of_seeds.begin(), list_of_seeds.end(),
                               Criteria());

  std::cout << "Number of vertices: " << cdt.number_of_vertices() << std::endl;
  std::cout << "Number of finite faces: " << cdt.number_of_faces() << std::endl;
  int mesh_faces_counter = 0;
  for(CDT::Finite_faces_iterator fit = cdt.finite_faces_begin();
      fit != cdt.finite_faces_end(); ++fit)
  {
    if(fit->is_in_domain()) ++mesh_faces_counter;
  }
  std::cout << "Number of faces in the mesh domain: " << mesh_faces_counter << std::endl;
}
