#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include <CGAL/Nef_3/SNC_indexed_items.h>
#include <CGAL/convex_decomposition_3.h>
#include <list>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron_3;
typedef CGAL::Nef_polyhedron_3<Kernel, CGAL::SNC_indexed_items> Nef_polyhedron_3;
typedef Nef_polyhedron_3::Volume_const_iterator Volume_const_iterator;

std::size_t run(std::string path)
{
  Polyhedron_3 input;
  std::ifstream(path) >> input;

  Nef_polyhedron_3 N(input);

  CGAL::convex_decomposition_3(N);
  std::list<Polyhedron_3> convex_parts;

  Volume_const_iterator ci = ++N.volumes_begin();
  for( ; ci != N.volumes_end(); ++ci) {
    if(ci->mark()) {
      Polyhedron_3 P;
      N.convert_inner_shell_to_polyhedron(ci->shells_begin(), P);
      convex_parts.push_back(P);
    }
  }

//  int i=0;
  for (const Polyhedron_3& P : convex_parts)
  {
//    std::ofstream("out_"+std::to_string(i++)+".off") << std::setprecision(17) << P;
    assert(P.size_of_vertices()!=0);
  }

  return convex_parts.size();
}

int main()
{
  std::size_t val = run("data/in1.off");
  assert(val==9);
  val = run("data/in2.off");
  assert(val==10);
  val = run("data/in3.off");
  assert(val==13);
  val = run("data/in4.off");
  assert(val==17);
}
