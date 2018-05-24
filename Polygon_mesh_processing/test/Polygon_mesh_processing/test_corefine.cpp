#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>

#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3> Surface_mesh;
typedef CGAL::Polyhedron_3<K> Polyhedron_3;

template <class TriangleMesh>
struct My_new_face_visitor
{
  typedef boost::graph_traits<TriangleMesh> GT;
  typedef typename GT::face_descriptor face_descriptor;

  void before_subface_creations(face_descriptor /*f_old*/,TriangleMesh&){}
  void after_subface_creations(TriangleMesh&){++(*i);}
  void before_subface_created(TriangleMesh&){}
  void after_subface_created(face_descriptor /*f_new*/,TriangleMesh&){}

  My_new_face_visitor()
    : i (new int(0) )
  {}

  boost::shared_ptr<int> i;
};

void test(const char* f1, const char* f2)
{
  std::cout << "Corefining " << f1
            << " and " << f2 << "\n";

  std::cout << "  with Surface_mesh\n";
  Surface_mesh sm1, sm2;
  std::ifstream input(f1);
  assert(input);
  input >> sm1;
  input.close();
  input.open(f2);
  assert(input);
  input >> sm2;
  input.close();
  My_new_face_visitor<Surface_mesh> sm_v;

  CGAL::Polygon_mesh_processing::corefine(sm1, sm2,
    CGAL::Polygon_mesh_processing::parameters::new_face_visitor(sm_v));

  assert(sm1.is_valid());
  assert(sm2.is_valid());
  assert(*(sm_v.i) != 0);

  std::cout << "  with Polyhedron_3\n";
  Polyhedron_3 P, Q;
  input.open(f1);
  assert(input);
  input >> P;
  input.close();
  input.open(f2);
  assert(input);
  input >> Q;
  My_new_face_visitor<Polyhedron_3> sm_p;

  CGAL::Polygon_mesh_processing::corefine(P, Q,
    CGAL::Polygon_mesh_processing::parameters::new_face_visitor(sm_p));
  assert(*(sm_p.i) != 0);

  assert(*(sm_v.i) == *(sm_p.i));
  assert(P.is_valid());
  assert(Q.is_valid());
}
int main(int argc, char** argv)
{
  for(int i=0; i< (argc-1)/2;++i)
  {
    test(argv[2*i+1], argv[2*(i+1)]);
    test(argv[2*(i+1)], argv[2*i+1]);
  }
}
