// Bug reported by a user against CGAL-4.0.2
// A single triangle cannot be meshed!
//
// The bug report is in the message:
//   Date: Wed, 25 Jul 2012 02:20:13 -0700 (PDT)
//   From: gerth <gert.hutter@gmail.com>
//   To: cgal-discuss@inria.fr
//   Subject: [cgal-discuss] Delauny_Mesher gives a CGAL::Assertion_exception
//
// The fault was a buggy assertion added in the following revision:
//   -----------------------------------------------------------------------
//   r66359 | lrineau | 2011-11-16 18:58:22 +0100 (Wed, 16 Nov 2011) | 3 lines
//
//   Add an assertion that checks that a point constructed as the midpoint of an
//   edge is located either on the edge or inside one of the two incident faces.
//
//   -----------------------------------------------------------------------
// And the assertion message was:
//   terminate called after throwing an instance of 'CGAL::Assertion_exception'
//     what():  CGAL ERROR: assertion violation!
//   Expr: zone.locate_type != Tr::FACE || zone.fh == f || zone.fh == n
//   File: /home/lrineau/wc/Mesh_2/test/Mesh_2/../../include/CGAL/Mesh_2/Refine_edges.h
//   Line: 430
//   Explanation: Your data set contains at least a vertex that is very close to 
//     a constrained edge! Mesh_2 cannot mesh that sort of data set.

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

  Vertex_handle va = cdt.insert(Point(0,0));
  Vertex_handle vb = cdt.insert(Point(-3.5240345245692879,-12.4));
  Vertex_handle vc = cdt.insert(Point(-9.3417057867627591,2.6237986627966166));

  cdt.insert_constraint(va, vb);
  cdt.insert_constraint(vb, vc);
  cdt.insert_constraint(vc, va);

  std::cout << "Number of vertices: " << cdt.number_of_vertices() << std::endl;

  std::cout << "Meshing the triangulation..." << std::endl;
  CGAL::refine_Delaunay_mesh_2(cdt, Criteria(0.125, 1.0));

  std::cout << "Number of vertices: " << cdt.number_of_vertices() << std::endl;
}
