#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Tetrahedral_remeshing/Remeshing_triangulation_3.h>
#include <CGAL/IO/File_medit.h>
#include <CGAL/Iso_cuboid_3.h>

#include <CGAL/tags.h>

#include <string>
#include <iostream>
#include <fstream>

int main()
{
  using K = CGAL::Exact_predicates_inexact_constructions_kernel;
  using Tr = CGAL::Tetrahedral_remeshing::Remeshing_triangulation_3<K>;

  // Open file elephant
  std::string filename = CGAL::data_file_path("meshes/elephant.mesh");
  std::ifstream in(filename, std::ios_base::in);
  if(!in) {
    std::cerr << "Error! Cannot open file " << filename << std::endl;
    return 1;
  }
  Tr tr;
  CGAL::IO::read_MEDIT(in, tr);
  assert(tr.is_valid());

  std::ofstream os("elephant_out.mesh");
  CGAL::IO::write_MEDIT(os, tr,
    CGAL::parameters::all_vertices(false).all_cells(true));
  os.close();

  CGAL::Bbox_3 bb(tr.finite_vertices_begin()->point().bbox());
  for (auto v : tr.finite_vertex_handles())
    bb = bb + v->point().bbox();

  CGAL::Iso_cuboid_3<K> isocuboid(bb);
  for (int i = 0; i < 8; ++i)
    tr.insert(isocuboid.vertex(i));

  std::ofstream os2("elephant_out_not_all_vertices.mesh");
  CGAL::IO::write_MEDIT(os2, tr,
    CGAL::parameters::all_vertices(false));
  os2.close();

  Tr tr2;
  std::ifstream is2("elephant_out_not_all_vertices.mesh");
  CGAL::IO::read_MEDIT(is2, tr2);;
  is2.close();
  assert(tr2.is_valid());

  return EXIT_SUCCESS;
}
