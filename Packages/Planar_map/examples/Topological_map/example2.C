// examples/Topological_map/example2.C
// -----------------------------------

#include <CGAL/basic.h>
#include <iostream>

#include <CGAL/Topological_map_bases.h>
#include <CGAL/Pm_default_dcel.h>
#include <CGAL/Topological_map.h>

class Face_with_info : public CGAL::Tpm_face_base {
  int inf;
public:
  Face_with_info() : CGAL::Tpm_face_base(), inf(0) {}

  int info() { return inf; }
  void set_info(int i) { inf = i; }
};

typedef CGAL::Pm_dcel<CGAL::Tpm_vertex_base,
                      CGAL::Tpm_halfedge_base,
                      Face_with_info > Dcel;  

typedef CGAL::Topological_map<Dcel> Tpm;

typedef Tpm::Halfedge_handle Halfedge_handle;
typedef Tpm::Vertex_handle   Vertex_handle;
typedef Tpm::Face_handle     Face_handle;

int main()
{
  
  Tpm t;

  Face_handle uf = t.unbounded_face();

  Halfedge_handle e1 = t.insert_in_face_interior(uf);
  CGAL_assertion(t.is_valid());
  std::cout << "Edge e1 inserted in unbounded face interior" << std::endl;

  Halfedge_handle e2 = t.insert_from_vertex(e1);
  CGAL_assertion(t.is_valid());
  std::cout << "Edge e2 inserted from target vertex of e1" << std::endl;

  Halfedge_handle e3 = t.insert_at_vertices(e2, e1->twin());
  CGAL_assertion(t.is_valid());
  std::cout << "Edge e3 inserted between target vertices of e2 and "
            << "twin of e1" << std::endl;
  
  std::cout << std::endl <<"Setting info of the new face to 10" << std::endl;
  Face_handle nf = e3->face();
  nf->set_info(10);

  std::cout << "Unbounded face info = " << uf->info() << std::endl;
  std::cout << "New face info = " << nf->info() << std::endl;

  return 0;
}
