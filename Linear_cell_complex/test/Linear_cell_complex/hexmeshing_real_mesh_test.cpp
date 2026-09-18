#include <CGAL/generate_hexahedral_mesh_using_two_refinement.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/config.h>
#include <string>

bool render_two_refinement(const std::string& file, int cube_cells_per_dim,
                           int nb_levels, bool trim, bool smooth)
{
  CGAL::Polyhedron_3<CGAL::Exact_predicates_inexact_constructions_kernel> poly;
  CGAL::IO::read_polygon_mesh(CGAL::data_file_path("meshes/" + file), poly);

  auto lcc=CGAL::generate_hexahedral_mesh_using_two_refinement
    (poly, cube_cells_per_dim, nb_levels,
     CGAL::parameters::use_trimming(trim).use_smoothing(smooth));

  bool res=lcc.is_valid();
  if(!res)
  { std::cout<<"[ERROR] lcc is not valid."<<std::endl; }

  std::size_t nb_non_hex=0;
  for(auto it=lcc.one_dart_per_cell<3>().begin(),
        itend=lcc.one_dart_per_cell<3>().end(); it!=itend; ++it)
  {
    if(!lcc.is_volume_combinatorial_hexahedron(it))
    { ++nb_non_hex; }
  }

  if(nb_non_hex>0)
  {
    std::cout<<"[ERROR] "<<nb_non_hex<<" volumes are not hexahedron."<<std::endl;
    res=false;
  }

  return res;
}

int main()
{
  bool res=true;

  if(!render_two_refinement("bunny00.off", 16, 0, true, true))
  {
    std::cout<<"[ERROR] in hexmeshing_real_mesh_test: test1."<<std::endl;
    res=false;
  }

  if(!render_two_refinement("dragon_res3.off", 16, 0, false, true))
  {
    std::cout<<"[ERROR] in hexmeshing_real_mesh_test: test2."<<std::endl;
    res=false;
  }

  if(!render_two_refinement("bunny00.off", 10, 1, false, false))
  {
    std::cout<<"[ERROR] in hexmeshing_real_mesh_test: test3."<<std::endl;
    res=false;
  }

  if(!render_two_refinement("dragon_res3.off", 10, 1, true, false))
  {
    std::cout<<"[ERROR] in hexmeshing_real_mesh_test: test4."<<std::endl;
    res=false;
  }

  if(!render_two_refinement("bunny00.off", 6, 2, true, true))
  {
    std::cout<<"[ERROR] in hexmeshing_real_mesh_test: test5."<<std::endl;
    res=false;
  }

  if(!render_two_refinement("dragon_res2.off", 6, 2, false, false))
  {
    std::cout<<"[ERROR] in hexmeshing_real_mesh_test: test6."<<std::endl;
    res=false;
  }

  if(!res)
  {  return EXIT_FAILURE; }

  return EXIT_SUCCESS;
}
