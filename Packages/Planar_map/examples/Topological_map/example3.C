//example3

#include <CGAL/basic.h>
#include <iostream>

#include <CGAL/Pm_default_dcel.h>
#include <CGAL/Topological_map.h>

typedef int  Point;
typedef int  Curve;

typedef CGAL::Pm_dcel<CGAL::Pm_vertex_base<Point>,
                     CGAL::Pm_halfedge_base<Curve>,
                     CGAL::Pm_face_base > Dcel;  

typedef CGAL::Topological_map<Dcel> Tpm;

typedef  Tpm::Halfedge_handle Halfedge_handle;
typedef  Tpm::Vertex_handle   Vertex_handle;
typedef  Tpm::Face_handle     Face_handle;


int main() {
  
  Point p1=1,p2=2,p3=3;
  Curve c1=1,c2=2,c3=3;

  Tpm t;

  Face_handle uf=t.unbounded_face();

  std::cout << "inserting curve " << c1 << " in face interior" << std::endl;
  Halfedge_handle e1 = t.insert_in_face_interior(uf);

  e1->set_curve(c1);e1->twin()->set_curve(c1);
  e1->source()->set_point(p1);
  e1->target()->set_point(p2);

  CGAL_assertion(t.is_valid());

  std::cout << "inserting curve " << c2 << " from vertex" << std::endl;
  Halfedge_handle e2=t.insert_from_vertex(e1);

  e2->set_curve(c2);
  e2->twin()->set_curve(c2);
  e2->target()->set_point(p3);

  CGAL_assertion(t.is_valid());

  std::cout << "inserting curve " << c3 << " between vertices " ;
  std::cout << e2->target()->point() << " and " ;
  std::cout << e1->twin()->target()->point() << std::endl;
  Halfedge_handle e3=t.insert_at_vertices(e2,e1->twin());

  e3->set_curve(c3);
  e3->twin()->set_curve(c3);
  
  CGAL_assertion(t.is_valid());

  Face_handle nf=e3->face() ;

  //checking that the new face has the right order of vertices ...,3,1,...
  Tpm::Ccb_halfedge_circulator cc=nf->outer_ccb();
  std::cout << "\norder of vertices in new face : ";
  do {
    std::cout << cc->source()->point() << " " ;
    ++cc;
    } while (cc != nf->outer_ccb());
  std::cout << std::endl;
    

  return 0;
}
