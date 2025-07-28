/*!
  \ingroup PkgLinearCellComplexExamples
  \brief Minimal example: read a `.3map` file, compute per-volume vertex count, and write to VTK.

  This example loads a 3D linear cell complex from a `.3map` file using 
  \cgal function `CGAL::load_combinatorial_map()`. It computes for each 
  3-cell (volume) the number of incident vertices (0-cells), stores these values 
  in a `std::vector<std::size_t>`, and writes the result to a `.vtk` file with 
  `CGAL::write_lcc_to_vtk()`, using the computed values as scalars for each volume.
*/

#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex_vtk_io.h>
#include <CGAL/Combinatorial_map_save_load.h>
#include <vector>
#include <fstream>
#include <cstdlib>

typedef CGAL::Linear_cell_complex_for_combinatorial_map<3> LCC;

int main()
{
    LCC lcc; CGAL::load_combinatorial_map("data/beam-with-mixed-cells.3map", lcc);
    std::vector<std::size_t> v;
    for(auto it = lcc.template one_dart_per_cell<3>().begin(), itend = lcc.template one_dart_per_cell<3>().end(); it != itend; ++it)
        v.push_back(std::distance(lcc.template one_dart_per_incident_cell<0,3>(lcc.dart_descriptor(*it)).begin(), lcc.template one_dart_per_incident_cell<0,3>(lcc.dart_descriptor(*it)).end()));
    CGAL::write_lcc_to_vtk(lcc, "data/beam-with-mixed-cells.vtk", nullptr, &v);
    return EXIT_SUCCESS;
}