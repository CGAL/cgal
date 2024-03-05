#include <CGAL/Simple_cartesian.h>

#include <CGAL/Isosurfacing_3/Cartesian_grid_3.h>
#include <CGAL/Isosurfacing_3/dual_contouring_3.h>
#include <CGAL/Isosurfacing_3/Dual_contouring_domain_3.h>
#include <CGAL/Isosurfacing_3/Finite_difference_gradient_3.h>
#include <CGAL/Isosurfacing_3/Interpolated_discrete_values_3.h>
#include <CGAL/Isosurfacing_3/marching_cubes_3.h>
#include <CGAL/Isosurfacing_3/Marching_cubes_domain_3.h>

#include <CGAL/Image_3.h>
#include <CGAL/IO/read_vtk_image_data.h>
#include <CGAL/Isosurfacing_3/IO/Image_3.h>
#include <CGAL/IO/polygon_soup_io.h>

#include <vtkNew.h>
#include <vtkImageData.h>
#include <vtkMetaImageReader.h>
#include <vtkXMLImageDataReader.h>
#include <vtkTIFFReader.h>
#include <vtkNrrdReader.h>
#include <vtkMINCImageReader.h>

#include <iostream>
#include <vector>

using Kernel = CGAL::Simple_cartesian<double>;
using FT = typename Kernel::FT;
using Point = typename Kernel::Point_3;

using Grid = CGAL::Isosurfacing::Cartesian_grid_3<Kernel>;
using Values = CGAL::Isosurfacing::Interpolated_discrete_values_3<Grid>;
using Gradients = CGAL::Isosurfacing::Finite_difference_gradient_3<Kernel>;

using Point_range = std::vector<Point>;
using Polygon_range = std::vector<std::vector<std::size_t> >;

namespace IS = CGAL::Isosurfacing;

void run_marching_cubes(const Grid& grid,
                        const FT isovalue,
                        const Values& values)
{
  using Domain = IS::Marching_cubes_domain_3<Grid, Values, IS::Linear_interpolation_edge_intersection>;

  std::cout << "\n ---- " << std::endl;
  std::cout << "Running Marching Cubes with isovalue = " << isovalue << std::endl;

  // fill up values

  // create a domain from the grid
  Domain domain { grid, values };

  // prepare collections for the output indexed soup
  Point_range points;
  Polygon_range triangles;

  // execute marching cubes
  IS::marching_cubes<CGAL::Parallel_if_available_tag>(domain, isovalue, points, triangles);

  std::cout << "Output #vertices: " << points.size() << std::endl;
  std::cout << "Output #triangles: " << triangles.size() << std::endl;

  // save output indexed mesh to a file, in the OFF format
  CGAL::IO::write_polygon_soup("marching_cubes_vtk_image.off", points, triangles);
}

void run_dual_contouring(const Grid& grid,
                         const FT isovalue,
                         const Values& values)
{
  using Domain = IS::Dual_contouring_domain_3<Grid, Values, Gradients, IS::Linear_interpolation_edge_intersection>;

  std::cout << "\n ---- " << std::endl;
  std::cout << "Running Dual Contouring with isovalue = " << isovalue << std::endl;

  // fill up values and gradients
  const FT step = CGAL::approximate_sqrt(grid.spacing().squared_length()) * 0.01; // finite difference step
  Gradients gradients { values, step };
  Domain domain { grid, values, gradients };

  Point_range points;
  Polygon_range triangles;

  // run dual contouring isosurfacing
  IS::dual_contouring<CGAL::Parallel_if_available_tag>(domain, isovalue, points, triangles);

  std::cout << "Output #vertices: " << points.size() << std::endl;
  std::cout << "Output #triangles: " << triangles.size() << std::endl;
  CGAL::IO::write_polygon_soup("dual_contouring_vtk_image.off", points, triangles);
}

template <typename VtkReader>
void run(const char* filename,
         const FT isovalue)
{
  vtkNew<VtkReader> reader;
  reader->SetFileName(filename);
  reader->Update();
  CGAL::Image_3 image = CGAL::IO::read_vtk_image_data(reader->GetOutput());

  // convert image to a Cartesian grid
  Grid grid;
  Values values { grid }; // 'values' keeps a reference to the grid
  if(!IS::IO::read_Image_3(image, grid, values))
  {
    std::cerr << "Error: Cannot convert image to Cartesian grid" << std::endl;
    return;
  }

  std::cout << "Span: " << grid.span() << std::endl;
  std::cout << "Cell dimensions: " << grid.spacing()[0] << " " << grid.spacing()[1] << " " << grid.spacing()[2] << std::endl;
  std::cout << "Cell #: " << grid.xdim() << ", " << grid.ydim() << ", " << grid.zdim() << std::endl;

  run_marching_cubes(grid, isovalue, values);

  run_dual_contouring(grid, isovalue, values);
}

int main(int argc, char* argv[])
{
  if(argc == 1)
  {
    std::cerr << "Usage: " << argv[0] << " <vtk image> <isovalue = 0>" << std::endl;
    return EXIT_FAILURE;
  }

  const char* filename = argv[1];
  const FT isovalue = (argc > 2) ? std::stod(argv[2]) : 0;

  const std::string ext = CGAL::IO::internal::get_file_extension(filename);
  if(ext == "mhd" || ext == "mha")
    run<vtkMetaImageReader>(filename, isovalue);
  else if(ext == "vti")
    run<vtkXMLImageDataReader>(filename, isovalue);
  else if(ext == "tif")
    run<vtkTIFFReader>(filename, isovalue);
  else if(ext == "nrrd")
    run<vtkNrrdReader>(filename, isovalue);
  else if(ext == "mnc")
    run<vtkMINCImageReader>(filename, isovalue);
  else
  {
    std::cerr << "Error: Unsupported file format" << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
