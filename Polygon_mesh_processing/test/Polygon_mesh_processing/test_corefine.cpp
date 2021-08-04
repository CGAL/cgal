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
struct My_visitor :
  public CGAL::Polygon_mesh_processing::Corefinement::Default_visitor<TriangleMesh>
{
  void after_subface_creations(TriangleMesh&){++(*i);}

  My_visitor()
    : i (new int(0) )
  {}

  std::shared_ptr<int> i;
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
  My_visitor<Surface_mesh> sm_v;

  std::size_t nb_v_before = num_vertices(sm1) + num_vertices(sm2);
  CGAL::Polygon_mesh_processing::corefine(sm1, sm2,
    CGAL::Polygon_mesh_processing::parameters::visitor(sm_v));
  std::size_t nb_v_after = num_vertices(sm1) + num_vertices(sm2);

  assert(sm1.is_valid());
  assert(sm2.is_valid());
  assert((*(sm_v.i) != 0) == (nb_v_before!=nb_v_after));

  std::cout << "  with Polyhedron_3\n";
  Polyhedron_3 P, Q;
  input.open(f1);
  assert(input);
  input >> P;
  input.close();
  input.open(f2);
  assert(input);
  input >> Q;
  My_visitor<Polyhedron_3> sm_p;

  nb_v_before = num_vertices(P) + num_vertices(Q);
  CGAL::Polygon_mesh_processing::corefine(P, Q,
    CGAL::Polygon_mesh_processing::parameters::visitor(sm_p));
  nb_v_after = num_vertices(P) + num_vertices(Q);

  assert((*(sm_p.i) != 0)  == (nb_v_before!=nb_v_after));

  assert(*(sm_v.i) == *(sm_p.i));
  assert(P.is_valid());
  assert(Q.is_valid());
}

void test_no_modifications(const char* f1, const char* f2)
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
  My_visitor<Surface_mesh> sm_v;

  std::size_t nb_v_before1 = vertices(sm1).size();
  std::size_t nb_v_before2 = vertices(sm2).size();

  CGAL::Polygon_mesh_processing::corefine(sm1, sm2,
    CGAL::parameters::visitor(sm_v),
    CGAL::parameters::do_not_modify(true));

  std::size_t nb_v_after1 = vertices(sm1).size();
  std::size_t nb_v_after2 = vertices(sm2).size();

  assert(sm1.is_valid());
  assert(sm2.is_valid());
  assert(nb_v_after2==nb_v_before2);

  assert((*(sm_v.i) != 0) == (nb_v_before1!=nb_v_after1));

  std::cout << "  with Polyhedron_3\n";
  Polyhedron_3 P, Q;
  input.open(f1);
  assert(input);
  input >> P;
  input.close();
  input.open(f2);
  assert(input);
  input >> Q;
  My_visitor<Polyhedron_3> sm_p;

  nb_v_before1 = vertices(P).size();
  nb_v_before2 = vertices(Q).size();
  CGAL::Polygon_mesh_processing::corefine(P, Q,
    CGAL::parameters::visitor(sm_p).do_not_modify(true));
  nb_v_after1 = vertices(P).size();
  nb_v_after2 = vertices(Q).size();

  assert(nb_v_after1==nb_v_before1);
  assert((*(sm_p.i) != 0)  == (nb_v_before2!=nb_v_after2));

  assert(P.is_valid());
  assert(Q.is_valid());
}

int main(int argc, char** argv)
{
  for(int i=0; i< (argc-1)/2;++i)
  {
    test(argv[2*i+1], argv[2*(i+1)]);
    test(argv[2*(i+1)], argv[2*i+1]);
    test_no_modifications(argv[2*(i+1)], argv[2*i+1]);
  }
}
