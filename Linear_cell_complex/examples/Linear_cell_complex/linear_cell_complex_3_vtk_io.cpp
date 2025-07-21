/*!
  \ingroup PkgLinearCellComplexExamples

  \brief Example showing how to read and write a 3D linear cell complex from/to a VTK file.

  This example loads a volumetric mesh from a `.vtk` file into a `Linear_cell_complex_for_combinatorial_map<3>`,
  optionally reads scalar fields (vertex/volume), and writes the modified structure back to a `.vtk` file.

  \cgalFeature{Linear_cell_complex_vtk_io}

  \cgalExampleOutput{
  Loaded LCC from data/lcc_input.vtk
   - 80 vertices
   - 52 volumes
  Wrote LCC to data/lcc_output.vtk
  }
*/


#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex_traits.h>
#include <CGAL/Linear_cell_complex_vtk_io.h>
#include <CGAL/IO/io.h>
#include <iostream>
#include <vector>

int main()
{
  using LCC = CGAL::Linear_cell_complex_for_combinatorial_map<3, 3, CGAL::Linear_cell_complex_traits<3>>;
  LCC lcc;

  std::vector<float> vertex_scalars, volume_scalars;

   const char* input_filename = "data/lcc_input.vtk";
  if (!CGAL::read_vtk(lcc, input_filename, &vertex_scalars, &volume_scalars)) {
    std::cerr << "Failed to read: " << input_filename << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "Loaded LCC from " << input_filename << std::endl;
  std::cout << " - " << lcc.number_of_vertex_attributes() << " vertices" << std::endl;
  std::cout << " - " << lcc.template one_dart_per_cell<3>().size() << " volumes" << std::endl;

  const char* output_filename = "data/lcc_output.vtk";
  if (!CGAL::write_vtk(lcc, output_filename, &vertex_scalars, &volume_scalars)) {
    std::cerr << "Failed to write: " << output_filename << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "Wrote LCC to " << output_filename << std::endl;
  return EXIT_SUCCESS;
}