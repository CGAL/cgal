#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex_vtk_io.h>

#include <cstdlib>
#include <fstream>
#include <cassert>
#include <cstdio>

typedef CGAL::Linear_cell_complex_for_combinatorial_map<3, 3> LCC;
typedef LCC::Point Point;

int main()
{
  LCC lcc_out;

  // Create a hexahedron
  auto d = lcc_out.make_hexahedron(
    lcc_out.create_vertex_attribute(Point(0, 0, 0)),
    lcc_out.create_vertex_attribute(Point(1, 0, 0)),
    lcc_out.create_vertex_attribute(Point(1, 1, 0)),
    lcc_out.create_vertex_attribute(Point(0, 1, 0)),
    lcc_out.create_vertex_attribute(Point(0, 0, 1)),
    lcc_out.create_vertex_attribute(Point(1, 0, 1)),
    lcc_out.create_vertex_attribute(Point(1, 1, 1)),
    lcc_out.create_vertex_attribute(Point(0, 1, 1))
  );

  std::vector<float> v_scalars = {1,2,3,4,5,6,7,8};
  std::vector<float> vol_scalars = {42.0f};

  const char* fname = "tmp_test_lcc_vtk.vtk";

  assert(CGAL::write_vtk(lcc_out, fname, &v_scalars, &vol_scalars));

  LCC lcc_in;
  std::vector<float> v_scalars_in, vol_scalars_in;

  assert(CGAL::read_vtk(lcc_in, fname, &v_scalars_in, &vol_scalars_in));

  assert(lcc_in.number_of_vertex_attributes() == lcc_out.number_of_vertex_attributes());
  assert(lcc_in.template number_of_cells<3>() == lcc_out.template number_of_cells<3>());
  assert(v_scalars_in.size() == v_scalars.size());
  assert(vol_scalars_in.size() == vol_scalars.size());

  std::remove(fname);

  return EXIT_SUCCESS;
}