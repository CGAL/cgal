/*!
  \ingroup PkgLinearCellComplexExamples
  \brief Minimal example: read a `.3map` file, compute per-volume vertex count, and write to VTK.

  This example loads a 3D linear cell complex from a `.3map` file using
  \cgal function `CGAL::load_combinatorial_map()`. It computes for each
  3-cell (volume) the number of incident vertices (0-cells), stores these values
  in a `std::vector<std::size_t>`, and writes the result to a `.vtk` file with
  `CGAL::IO::write_VTK()`, using the computed values as scalars for each volume.
*/

#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex/IO/VTK.h>
#include <CGAL/Combinatorial_map_save_load.h>
#include <vector>
#include <cstdlib>

int main()
{
  CGAL::Linear_cell_complex_for_combinatorial_map<3> lcc;
  CGAL::load_combinatorial_map("data/beam-with-mixed-cells.3map", lcc);

    // Compute per-volume vertex count
    std::vector<std::size_t> volume_scalars;
    for(auto it=lcc.template one_dart_per_cell<3>().begin(),
         itend=lcc.template one_dart_per_cell<3>().end(); it!=itend; ++it)
    {
      std::size_t nbv=lcc.template one_dart_per_incident_cell<0,3>(it).size();
      volume_scalars.push_back(nbv);
    }

    if(!CGAL::IO::write_VTK(lcc, "beam-with-mixed-cells.vtk", nullptr,
                            &volume_scalars))
    { return EXIT_FAILURE; }

    return EXIT_SUCCESS;
}
