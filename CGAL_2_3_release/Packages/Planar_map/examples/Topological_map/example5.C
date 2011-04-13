// examples/Topological_map/example5.C
// -----------------------------------
#include <CGAL/basic.h>
#include <iostream>

#include <CGAL/Topological_map_bases.h>
#include <CGAL/Pm_default_dcel.h>
#include <CGAL/Topological_map.h>


typedef CGAL::Pm_dcel<CGAL::Tpm_vertex_base,
                     CGAL::Tpm_halfedge_base,
                     CGAL::Tpm_face_base> Dcel;  

typedef CGAL::Topological_map<Dcel> Tpm;

typedef  Tpm::Halfedge_handle Halfedge_handle;
typedef  Tpm::Vertex_handle   Vertex_handle;
typedef  Tpm::Face_handle     Face_handle;


int main() {

  Tpm t;

  Face_handle uf=t.unbounded_face();

  std::cout << "inserting edge e1 in face interior ..." ;
  Halfedge_handle e1 = t.insert_in_face_interior(uf);
  CGAL_assertion(t.is_valid());
  std::cout << "map is valid." << std::endl;

  std::cout << "inserting edge e2 from target vertex of e1 ..." ;
  Halfedge_handle e2=t.insert_from_vertex(e1);
  CGAL_assertion(t.is_valid());
  std::cout <<  "map is valid." << std::endl;

  std::cout << "inserting edge e3 between target vertices of e2 and e1->twin() ...";
  Halfedge_handle e3=t.insert_at_vertices(e2,e1->twin());
  CGAL_assertion(t.is_valid());
  std::cout << "map is valid." << std::endl;


  std::cout << "splitting edge ...";
  Halfedge_handle e4=t.split_edge(e3);
  CGAL_assertion(t.is_valid());
  std::cout << "map is valid." << std::endl;

  std::cout << "merging the split edge ...";
  Halfedge_handle e5=t.merge_edge(e4,e4->next_halfedge());
  CGAL_assertion(t.is_valid());
  std::cout << "map is valid." << std::endl;

  std::cout << "removing edge ...";
#ifndef CGAL_NO_ASSERTIONS // in order to avoid warnings
  Face_handle f=
#endif
    t.remove_edge(e5); 
  CGAL_assertion(t.is_valid());
  std::cout << "map is valid." << std::endl;

  CGAL_assertion(f==uf); //asserting the returned face is the unbounded face
  


  return 0;
}
