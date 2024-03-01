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

int main(int argc, char* argv[])
{
  if(argc == 1)
  {
    std::cerr << "Usage: " << argv[0] << " <vtk image> <isovalue = 0>" << std::endl;
    return EXIT_FAILURE;
  }

  const std::string filename = argv[1];
  const FT isovalue = (argc > 2) ? std::stod(argv[2]) : 0;

  CGAL::Image_3 image;

  std::string ext = CGAL::IO::internal::get_file_extension(filename);
  if(ext == "mhd" || ext == "mha")
  {
    // load volumetric image from a file
    vtkNew<vtkMetaImageReader> reader;
    reader->SetFileName(argv[1]);
    reader->Update();
    image = CGAL::IO::read_vtk_image_data(reader->GetOutput());
  }
  else if(ext == "vti")
  {
    vtkNew<vtkXMLImageDataReader> reader;
    reader->SetFileName(argv[1]);
    reader->Update();
    image = CGAL::IO::read_vtk_image_data(reader->GetOutput());
  }
  else if(ext == "tif")
  {
    vtkNew<vtkTIFFReader> reader;
    reader->SetFileName(argv[1]);
    reader->Update();
    image = CGAL::IO::read_vtk_image_data(reader->GetOutput());
  }
  else
  {
    std::cerr << "Error: Unsupported file format" << std::endl;
    return EXIT_FAILURE;
  }

  // convert image to a Cartesian grid
  Grid grid;
  Values values { grid }; // 'values' keeps a reference to the grid
  if(!IS::IO::read_Image_3(image, grid, values))
  {
    std::cerr << "Error: Cannot convert image to Cartesian grid" << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "Span: " << grid.span() << std::endl;
  std::cout << "Cell dimensions: " << grid.spacing()[0] << " " << grid.spacing()[1] << " " << grid.spacing()[2] << std::endl;
  std::cout << "Cell #: " << grid.xdim() << ", " << grid.ydim() << ", " << grid.zdim() << std::endl;

  for (std::size_t i=0; i<grid.xdim(); ++i)
    for (std::size_t j=0; j<grid.ydim(); ++j)
      for (std::size_t k=0; k<grid.zdim(); ++k)
        values(i, j, k) = - values(i, j, k); // inside out

  run_marching_cubes(grid, isovalue, values);

  run_dual_contouring(grid, isovalue, values);

  return EXIT_SUCCESS;
}
