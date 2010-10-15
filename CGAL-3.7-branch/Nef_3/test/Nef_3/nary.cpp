
#define CGAL_NEF3_SORT_OUTPUT 1

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include <CGAL/Nef_3/SNC_indexed_items.h>
#include <CGAL/Nef_nary_union_3.h>
#include <CGAL/Nef_nary_intersection_3.h>
#include <fstream>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef CGAL::Nef_polyhedron_3<Kernel, CGAL::SNC_indexed_items> Nef_polyhedron;
typedef CGAL::Nef_nary_union_3<Nef_polyhedron> Union;
typedef CGAL::Nef_nary_intersection_3<Nef_polyhedron> Intersection;
typedef Kernel::Vector_3 Vector_3;
typedef Kernel::Aff_transformation_3 Aff_transformation_3;

int main()
{
  Nef_polyhedron N;
  std::ifstream in("data/cube.nef3.SH");
  in >> N;

  Nef_polyhedron C0, C1, C2;
  C0.transform(Aff_transformation_3(CGAL::TRANSLATION, Vector_3(1, 0, 0, 100)));
  C1.transform(Aff_transformation_3(CGAL::TRANSLATION, Vector_3(0, 1, 0, 100)));
  C2.transform(Aff_transformation_3(CGAL::TRANSLATION, Vector_3(0, 0, 1, 100)));

  Union u;
  u.add_polyhedron(C0);
  u.add_polyhedron(C1);
  u.add_polyhedron(C2);
  u.get_union();
  u.add_polyhedron(N);

  Intersection i;
  i.add_polyhedron(C0);
  i.add_polyhedron(C1);
  i.add_polyhedron(C2);
  i.get_intersection();
  i.add_polyhedron(N);

  return 0;
}
