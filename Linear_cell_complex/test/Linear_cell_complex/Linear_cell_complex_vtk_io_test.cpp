#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex/IO/VTK.h>
#include <cassert>
#include <vector>
#include <cstdlib>

typedef CGAL::Linear_cell_complex_for_combinatorial_map<3, 3> LCC;

int main() {
  LCC lcc1, lcc2;
  std::vector<float>        vertex_scalars1, vertex_scalars2;
  std::vector<std::size_t>  volume_scalars1, volume_scalars2;

  const char* filename = "data/beam-with-mixed-cells.vtk";

  
  bool ok_read1 = CGAL::IO::read_VTK<LCC,float,std::size_t>(lcc1, filename,
                                                            &vertex_scalars1, &volume_scalars1);
  assert(ok_read1);

  
  std::size_t nb_vertices = 0;
  for(auto itv = lcc1.vertex_attributes().begin(), itvend = lcc1.vertex_attributes().end(); itv != itvend; ++itv)
    ++nb_vertices;
  if(vertex_scalars1.size() != nb_vertices) {
    vertex_scalars1.resize(nb_vertices);
    for(std::size_t i=0;i<nb_vertices;++i)
      vertex_scalars1[i] = static_cast<float>(i);
  }

  
  std::size_t nb_volumes = 0;
  for(auto itvol = lcc1.one_dart_per_cell<3>().begin(),
           itvolend = lcc1.one_dart_per_cell<3>().end(); itvol != itvolend; ++itvol)
    ++nb_volumes;

  if(volume_scalars1.size() != nb_volumes) {
    volume_scalars1.clear();
    volume_scalars1.reserve(nb_volumes);
    for(auto itvol = lcc1.one_dart_per_cell<3>().begin(),
             itvolend = lcc1.one_dart_per_cell<3>().end(); itvol != itvolend; ++itvol) {
      auto d = lcc1.dart_descriptor(*itvol);
      std::size_t nbv = lcc1.template one_dart_per_incident_cell<0,3>(d).size();
      volume_scalars1.push_back(nbv);
    }
  }

  
  bool ok_write = CGAL::IO::write_VTK<LCC,float,std::size_t>(lcc1, filename,
                                                             &vertex_scalars1, &volume_scalars1);
  assert(ok_write);

  
  bool ok_read2 = CGAL::IO::read_VTK<LCC,float,std::size_t>(lcc2, filename,
                                                            &vertex_scalars2, &volume_scalars2);
  assert(ok_read2);

  assert(lcc1.is_isomorphic_to(lcc2, false, true, true));
  assert(vertex_scalars1.size() == vertex_scalars2.size());
  assert(volume_scalars1.size() == volume_scalars2.size());
  for(std::size_t i=0;i<vertex_scalars1.size();++i)
    assert(vertex_scalars1[i] == vertex_scalars2[i]);
  for(std::size_t i=0;i<volume_scalars1.size();++i)
    assert(volume_scalars1[i] == volume_scalars2[i]);

  return EXIT_SUCCESS;
}
