#include<CGAL/Gmpz.h>
#include<CGAL/Homogeneous.h>
#include<CGAL/Nef_polyhedron_3.h>
#include<CGAL/IO/Nef_polyhedron_iostream_3.h>
#include<CGAL/Polyhedron_3.h>
#include<CGAL/convex_decomposition_3.h>
#include<CGAL/Nef_3/Nary_union.h>
#include<CGAL/Nef_3/SNC_indexed_items.h>
#include<CGAL/Minkowski_sum_3/Gaussian_map.h>
#include<CGAL/Minkowski_sum_3/Gaussian_map_to_nef_3.h>

typedef CGAL::Gmpz RT;
typedef CGAL::Homogeneous<RT> Kernel;
typedef CGAL::Nef_polyhedron_3<Kernel, CGAL::SNC_indexed_items> Nef_polyhedron_3;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron_3;
typedef CGAL::Nary_union<Nef_polyhedron_3> Nary_union;
typedef Nef_polyhedron_3::SFace_const_iterator SFace_const_iterator;

void test_convex_parts(Nef_polyhedron_3& N)
{
  //  Nef_polyhedron_3 N(C);
  CGAL::convex_decomposition_3(N);
  std::cerr << N;
  std::list<Nef_polyhedron_3> convex_parts;
  Nef_polyhedron_3::Volume_const_iterator ci;
  for(ci = ++N.volumes_begin(); ci != N.volumes_end(); ++ci) {
    if(ci->mark()) {
      Polyhedron_3 P;
      N.convert_inner_shell_to_polyhedron(ci->shells_begin(), P);
      Nef_polyhedron_3 tmp0(P), tmp1;

      CGAL::Gaussian_map<Kernel,Nef_polyhedron_3> G(N, ci);
      CGAL::Gaussian_map_to_nef_3<Nef_polyhedron_3> 
	Convertor(G);
      tmp1.delegate(Convertor, true);

      CGAL_assertion(tmp1.is_valid());      
      CGAL_assertion(tmp1.closure().symmetric_difference(tmp0).is_empty());
    }
  }
}

int main(int argc, char* argv[])
{
  Nef_polyhedron_3 N;
  if(argc == 1) {
    std::ifstream in("star.nef3");
    in >> N;
  } else {
    std::ifstream in(argv[1]);
    in >> N;
  }
  test_convex_parts(N);
}
