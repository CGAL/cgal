#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex_vtk_io.h>
#include <cassert>
#include <string>
#include <unordered_map>
#include <vector>

typedef CGAL::Linear_cell_complex_for_combinatorial_map<3, 3> LCC;

int main() {
  LCC lcc1, lcc2;
  std::vector<float> v_scalars1, vol_scalars1;
  std::vector<float> v_scalars2, vol_scalars2;

  const std::string input_file = "data/beam-with-mixed-cells.vtk";
  assert(CGAL::read_lcc_from_vtk(lcc1, input_file.c_str(), &v_scalars1, &vol_scalars1));

  // Build index maps
  std::unordered_map<LCC::Vertex_attribute_const_descriptor, std::size_t> vertex_indices;
  std::size_t idx = 0;
  for(auto itv = lcc1.vertex_attributes().begin(), itvend = lcc1.vertex_attributes().end(); itv != itvend; ++itv)
    vertex_indices[itv] = idx++;
  std::unordered_map<LCC::Dart_const_descriptor, std::size_t> volume_indices;
  idx = 0;
  for(auto itvol = lcc1.one_dart_per_cell<3>().begin(), itvolend = lcc1.one_dart_per_cell<3>().end(); itvol != itvolend;
      ++itvol)
    volume_indices[itvol] = idx++;

  const char* tmp_file = "tmp_test_lcc_vtk.vtk";
  assert(CGAL::write_lcc_to_vtk(
    lcc1, tmp_file,
    [&v_scalars1, &vertex_indices](const LCC& lcc, LCC::Dart_const_descriptor d) {
      return v_scalars1[vertex_indices.at(lcc.attribute<0>(d))];
    },
    [&vol_scalars1, &volume_indices](const LCC& lcc, LCC::Dart_const_descriptor d) {
      return vol_scalars1[volume_indices.at(d)];
    }));
  assert(CGAL::read_lcc_from_vtk(lcc2, tmp_file, &v_scalars2, &vol_scalars2));

  assert(lcc1.is_isomorphic_to(lcc2, false, true, true));

  std::remove(tmp_file);

  return EXIT_SUCCESS;
}