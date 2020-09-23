#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>

#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef bool Vertex_info;
typedef bool Face_info;

typedef CGAL::Triangulation_vertex_base_with_info_2<Vertex_info,K> Vb;
typedef CGAL::Triangulation_face_base_with_info_2<Face_info,K> Fb_w_i;
typedef CGAL::Constrained_triangulation_face_base_2<K,Fb_w_i> C_fb_w_i;
typedef CGAL::Delaunay_mesh_face_base_2<K,C_fb_w_i> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb> Tds;
typedef CGAL::Constrained_Delaunay_triangulation_2<K,Tds> CDT;

typedef CDT::Point Point;
typedef CDT::Vertex_handle Vertex_handle;
typedef CDT::Finite_vertices_iterator Finite_vertices_iterator;
typedef CDT::Finite_faces_iterator Finite_faces_iterator;
typedef CDT::Face_circulator Face_circulator;

int main()
{
  CDT cdt;

  cdt.insert(Point(0,0));
  cdt.insert(Point(1,1));
  cdt.insert(Point(-1,1));
  cdt.insert(Point(-1,-1));

  for (Finite_faces_iterator fit=cdt.finite_faces_begin();
       fit != cdt.finite_faces_end();
       ++fit)
    fit->info() = true;
}
