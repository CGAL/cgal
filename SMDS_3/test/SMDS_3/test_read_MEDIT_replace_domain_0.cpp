#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Tetrahedral_remeshing/Remeshing_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/IO/File_medit.h>

#include <cassert>
#include <set>
#include <sstream>
#include <string>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Tr = CGAL::Tetrahedral_remeshing::Remeshing_triangulation_3<K>;
using C3t3 = CGAL::Mesh_complex_3_in_triangulation_3<Tr>;

std::string mesh_data(bool cgal_mesh)
{
  return std::string("MeshVersionFormatted 1\nDimension 3\n") +
    (cgal_mesh ? "# CGAL::Mesh_complex_3_in_triangulation_3\n" : "") +
    "Vertices\n5\n"
    "0 0 0 0\n1 0 0 0\n0 1 0 0\n0 0 1 0\n0 0 -1 0\n"
    "Tetrahedra\n2\n1 2 3 4 0\n1 3 2 5 7\nEnd\n";
}

void check_subdomains(const Tr& tr, bool replaced)
{
  assert(tr.is_valid());
  assert(tr.number_of_vertices() == 5);
  assert(tr.number_of_finite_cells() == 2);

  std::set<int> subdomains;
  for(auto cell : tr.finite_cell_handles())
    subdomains.insert(cell->subdomain_index());
  assert(subdomains == (replaced ? std::set<int>{7, 8} : std::set<int>{0, 7}));

  for(auto cell : tr.all_cell_handles())
    if(tr.is_infinite(cell))
      assert(cell->subdomain_index() == 0);
}

int main()
{
  for(bool cgal_mesh : {false, true})
  {
    Tr tr;
    std::istringstream input(mesh_data(cgal_mesh));
    bool read = CGAL::IO::read_MEDIT(input, tr);
    assert(read);
    check_subdomains(tr, false);

    C3t3 c3t3;
    std::istringstream complex_input(mesh_data(cgal_mesh));
    read = CGAL::IO::read_MEDIT(complex_input, c3t3);
    assert(read);
    check_subdomains(c3t3.triangulation(), false);
    assert(c3t3.number_of_cells_in_complex() == 1);

    for(bool replace : {false, true})
    {
      const bool replaced = replace && !cgal_mesh;
      Tr tr_with_parameter;
      std::istringstream parameter_input(mesh_data(cgal_mesh));
      read = CGAL::IO::read_MEDIT(parameter_input, tr_with_parameter,
                                 CGAL::parameters::replace_domain_0(replace));
      assert(read);
      check_subdomains(tr_with_parameter, replaced);

      C3t3 complex_with_parameter;
      std::istringstream complex_parameter_input(mesh_data(cgal_mesh));
      read = CGAL::IO::read_MEDIT(complex_parameter_input, complex_with_parameter,
                                 CGAL::parameters::replace_domain_0(replace));
      assert(read);
      check_subdomains(complex_with_parameter.triangulation(), replaced);
      assert(complex_with_parameter.number_of_cells_in_complex() == (replaced ? 2 : 1));
    }
  }
}
