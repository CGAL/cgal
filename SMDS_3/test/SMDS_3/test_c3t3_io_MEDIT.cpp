#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Tetrahedral_remeshing/Remeshing_triangulation_3.h>
#include <CGAL/IO/File_medit.h>
#include <CGAL/Iso_cuboid_3.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
#include <CGAL/make_mesh_3.h>

#include <CGAL/tags.h>

#include <string>
#include <iostream>
#include <fstream>

int test_MEDIT_with_features()
{
  using K = CGAL::Exact_predicates_inexact_constructions_kernel;
  using Polyhedron = CGAL::Mesh_polyhedron_3<K>::type;
  using Mesh_domain = CGAL::Polyhedral_mesh_domain_with_features_3<K>;
  using Tr = CGAL::Mesh_triangulation_3<Mesh_domain>::type;
  using C3t3 = CGAL::Mesh_complex_3_in_triangulation_3<Tr, Mesh_domain::Corner_index, Mesh_domain::Curve_index>;
  using Mesh_criteria = CGAL::Mesh_criteria_3<Tr>;

  const std::string fname = CGAL::data_file_path("meshes/fandisk.off");
  std::ifstream input(fname);
  Polyhedron polyhedron;
  input >> polyhedron;
  if(input.fail()) {
    std::cerr << "Error: Cannot read file " << fname << std::endl;
    return EXIT_FAILURE;
  }

  if(!CGAL::is_triangle_mesh(polyhedron)) {
    std::cerr << "Input geometry is not triangulated." << std::endl;
    return EXIT_FAILURE;
  }

  Mesh_domain domain(polyhedron);
  domain.detect_features();
  Mesh_criteria criteria(CGAL::parameters::edge_size(0.025)
                                          .facet_distance(0.005)
                                          .cell_size(0.05));
  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria);

  const std::size_t nb_facets = c3t3.number_of_facets();
  const std::size_t nb_edges = c3t3.number_of_edges();
  const std::size_t nb_corners = c3t3.number_of_corners();

  std::ofstream os("fandisk_out.mesh");
  CGAL::IO::write_MEDIT(os, c3t3,
                        CGAL::parameters::all_vertices(true).all_cells(true));
  os.close();

  std::ifstream is("fandisk_out.mesh");
  C3t3 c3t3_in;
  CGAL::IO::read_MEDIT(is, c3t3_in);

  assert(nb_facets == c3t3_in.number_of_facets());
  assert(nb_edges == c3t3_in.number_of_edges());
  assert(nb_corners == c3t3_in.number_of_corners());

  return EXIT_SUCCESS;
}

int test()
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

  // non convex
  std::ofstream os3("elephant_out_not_all_cells.mesh");
  CGAL::IO::write_MEDIT(os3, tr, CGAL::parameters::all_cells(false));
  os3.close();

  Tr tr3;
  std::ifstream is3("elephant_out_not_all_cells.mesh");
  CGAL::IO::read_MEDIT(is3, tr3);
  is3.close();
  assert(tr3.is_valid());

  return EXIT_SUCCESS;
}

int main()
{
  if(test_MEDIT_with_features() != EXIT_SUCCESS)
    return EXIT_FAILURE;

  return test();
}
