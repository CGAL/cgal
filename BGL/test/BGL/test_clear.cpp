#include "test_Prefix.h"

#include <CGAL/boost/graph/helpers.h>

template<typename Mesh>
void test() {
  const std::string fname = "data/7_faces_triangle.off";
  Mesh m;
  if(!read_a_mesh(m, fname)) {
    std::cout << "Error reading file: " << fname << std::endl;
  }

  CGAL::clear(m);
  assert(num_vertices(m) == 0);
  assert(num_faces(m) == 0);
  assert(num_edges(m) == 0);
  assert(is_valid(m));
}

int main()
{
  test<SM>();
  test<Polyhedron>();
#if defined(CGAL_USE_OPENMESH)
  test<OMesh>();
#endif
  return 0;
}
