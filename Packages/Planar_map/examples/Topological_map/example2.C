//example2

#include <CGAL/basic.h>
#include <iostream>

#include <CGAL/Topological_map_bases.h>
#include <CGAL/Pm_default_dcel.h>
#include <CGAL/Topological_map.h>

using namespace CGAL;

class Face_with_info : public Tpm_face_base {
  int inf;
public:
  Face_with_info() : Tpm_face_base(), inf(0) {}

  int info() {return inf;}
  void set_info(int i) {inf=i;}
};



typedef Pm_dcel<Tpm_vertex_base,
                     Tpm_halfedge_base,
                     Face_with_info > Dcel;  

typedef Topological_map<Dcel> Tpm;

typedef  Tpm::Halfedge_handle Halfedge_handle;
typedef  Tpm::Vertex_handle   Vertex_handle;
typedef  Tpm::Face_handle     Face_handle;


int main() {
  
  Tpm t;

  Face_handle uf=t.unbounded_face();

  std::cout << "inserting e1 in face interior..." << std::endl;
  Halfedge_handle e1 = t.insert_in_face_interior(uf);
  CGAL_assertion(t.is_valid());

  std::cout << "inserting e2 from vertex..." << std::endl;
  Halfedge_handle e2=t.insert_from_vertex(e1);
  CGAL_assertion(t.is_valid());

  std::cout << "inserting e3 between vertices of e2 and e1->twin()..." << std::endl;
  Halfedge_handle e3=t.insert_at_vertices(e2,e1->twin());
  CGAL_assertion(t.is_valid());
  
  std::cout << "\nsetting info of the new face to 10..." << std::endl;
  Face_handle nf=e3->face() ;
  nf->set_info(10);

  std::cout << "\nunbounded face info = " << uf->info() << std::endl;
  std::cout << "new face info = " << nf->info() << std::endl;

  return 0;
}
