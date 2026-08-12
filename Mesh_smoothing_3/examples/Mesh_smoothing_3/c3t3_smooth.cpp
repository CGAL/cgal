#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Triangulation_3.h>
#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Simplicial_mesh_cell_base_3.h>
#include <CGAL/Simplicial_mesh_vertex_base_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>

#include <CGAL/tetrahedral_remeshing.h>

#include <CGAL/tags.h>

#include <CGAL/IO/File_medit.h>
#include <fstream>

#include <CGAL/Mesh_smoothing_3/boundary_aware_mesh_smoothing.h>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;

using Subdomain_index = int;
using Surface_patch_index = unsigned char;
using Curve_index = char;
using Corner_index = short;

using Cb = CGAL::Simplicial_mesh_cell_base_3<K, Subdomain_index, Surface_patch_index>;
using Vb = CGAL::Simplicial_mesh_vertex_base_3<K, Subdomain_index, Surface_patch_index,
                                                  Curve_index, Corner_index>;

using Tds = CGAL::Triangulation_data_structure_3<Vb, Cb, CGAL::Sequential_tag>;
using Triangulation = CGAL::Triangulation_3<K, Tds>;

using C3t3 = CGAL::Mesh_complex_3_in_triangulation_3<Triangulation>;


int main(int argc, char* argv[])
{
    std::string filename = (argc > 1) ? std::string(argv[1])
                        : "../data/mambo_m3.mesh";

    C3t3 c3t3;

    std::ifstream is(filename, std::ios_base::in);
    Triangulation tr;
    if(!CGAL::IO::read_MEDIT(is,tr))
    {
        std::cerr << "Failed to read" << std::endl;
        return EXIT_FAILURE;
    }
    c3t3.set_triangulation(std::move(tr));

    std::ofstream os("c3t3_initial.mesh");
    CGAL::IO::write_MEDIT(os, c3t3.triangulation(), CGAL::parameters::all_vertices(true));
    os.close();

    CGAL::boundary_aware_mesh_smoothing(c3t3, CGAL::parameters::verbose(true).number_of_iterations(100));

    std::ofstream os2("c3t3_smoothed.mesh");
    CGAL::IO::write_MEDIT(os2, c3t3.triangulation(), CGAL::parameters::all_vertices(true));
    os2.close();

    std::cout << "Done" << std::endl;
    return EXIT_SUCCESS;

}




