#include <CGAL/Exact_integer.h>
#include<CGAL/Homogeneous.h>
#include<CGAL/Nef_polyhedron_3.h>
#include<CGAL/IO/Nef_polyhedron_iostream_3.h>
#include<CGAL/Polyhedron_3.h>
#include<CGAL/convex_decomposition_3.h>
#include<CGAL/Nef_nary_union_3.h>
#include<CGAL/Nef_3/SNC_indexed_items.h>
#include<fstream>

typedef CGAL::Exact_integer  RT;
typedef CGAL::Homogeneous<RT> Kernel;
typedef CGAL::Nef_polyhedron_3<Kernel> Nef_polyhedron_3;
//typedef CGAL::Nef_polyhedron_3<Kernel, CGAL::SNC_indexed_items> Nef_polyhedron_3;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron_3;
typedef CGAL::Nef_nary_union_3<Nef_polyhedron_3> Nary_union;

void check_decomposition(Nef_polyhedron_3& N)
{
  CGAL::convex_decomposition_3(N);
  std::list<Nef_polyhedron_3> convex_parts;
  Nef_polyhedron_3::Volume_const_iterator ci;
  for(ci = ++N.volumes_begin(); ci != N.volumes_end(); ++ci) {
    if(ci->mark()) {
      Polyhedron_3 P;
      N.convert_inner_shell_to_polyhedron(ci->shells_begin(), P);
      Nef_polyhedron_3 tmp(P);
      convex_parts.push_back(tmp);
    }
  }

  std::list<Nef_polyhedron_3>::iterator pi0, pi1;
  for(pi0 = convex_parts.begin(); pi0 != convex_parts.end(); ++pi0) {
    pi1 = pi0;
    for(++pi1; pi1 != convex_parts.end(); ++pi1) {
      Nef_polyhedron_3 tmp = pi0->intersection(*pi1).interior();
      assert(tmp.is_empty());
    }
  }

  Nary_union nu;
  for(pi0 = convex_parts.begin(); pi0 != convex_parts.end(); ++pi0)
    nu.add_polyhedron(*pi0);
  assert(nu.get_union().symmetric_difference(N).is_empty());
}

int main()
{
  Nef_polyhedron_3 N;
  std::cin >> N;
  check_decomposition(N);
}
