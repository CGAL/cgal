#define CGAL_DEBUG_CLIPPING

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/clip.h>

#include <CGAL/IO/OFF_reader.h>
#include <CGAL/IO/STL_reader.h>

#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel             K;
// typedef CGAL::Exact_predicates_exact_constructions_kernel               K;

typedef CGAL::Surface_mesh<K::Point_3>                                  Mesh;

typedef K::Point_3                                                      Point_3;
typedef K::Plane_3                                                      Plane_3;

namespace PMP = CGAL::Polygon_mesh_processing;

template <typename K, typename Mesh>
bool read_mesh(const char* filename, Mesh& sm)
{
  typedef typename K::Point_3                                   Point;

  std::ifstream in(filename, std::ios::binary);
  if(!in.good())
  {
    std::cerr << "Error: can't read file: " << filename << std::endl;
    std::exit(1);
  }

  std::string fn(filename);
  if(fn.substr(fn.find_last_of(".") + 1) == "stl")
  {
    std::vector<Point> points;
    std::vector<std::vector<int> > faces;
    if(!CGAL::read_STL(in, points, faces))
    {
      std::cerr << "Error: cannot read STL mesh\n";
      return false;
    }

    std::cout << "Cleaning polygon soup..." << std::endl;
    PMP::repair_polygon_soup(points, faces);

    if(!CGAL::Polygon_mesh_processing::orient_polygon_soup(points, faces))
      std::cerr << "W: File does not describe a polygon mesh" << std::endl;

    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, faces, sm);
  }
  else if(fn.substr(fn.find_last_of(".") + 1) == "off")
  {
    if(!(in >> sm))
    {
      std::cerr << "Cannot read .OFF as polygon mesh, reading .OFF as a polygon soup...\n";

      in.clear(); // clear fail and eof bits
      in.seekg(0, std::ios::beg); // back at the start

      std::vector<Point> points;
      std::vector<std::vector<int> > faces;
      if(!CGAL::read_OFF(in, points, faces))
      {
        std::cerr << "Error: cannot read OFF mesh as a polygon soup\n";
        return false;
      }

      std::cout << "Cleaning polygon soup..." << std::endl;
      PMP::repair_polygon_soup(points, faces);

      if(!CGAL::Polygon_mesh_processing::orient_polygon_soup(points, faces))
        std::cerr << "W: File does not describe a polygon mesh" << std::endl;

      CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, faces, sm);
    }
  }
  else
  {
    std::cerr << "Unknown file type" << std::endl;
    return false;
  }

  if(!CGAL::is_triangle_mesh(sm))
  {
    std::cerr << "Input geometry is not triangulated." << std::endl;
    return false;
  }

  return true;
}

int main(int argc, char** argv)
{
  const char* filename = (argc > 1) ? argv[1] : "data/mannequin-devil.off";
  Mesh mesh;
  if(!read_mesh<K>(filename, mesh))
    return EXIT_FAILURE;

  std::cout << "Input mesh: " << vertices(mesh).size() << " nv " << faces(mesh).size() << " nf " << std::endl;

#if 0//def CGAL_USE_CLIPPER_PLANE
  Plane_3 pl(Point_3(0,0,-600), Point_3(1,0,-600), Point_3(0,1,-600));
  PMP::clip(mesh, pl,
            CGAL::parameters::geom_traits(K())
                             .input_might_have_self_intersection(true));
#else
  filename = (argc > 2) ? argv[2] : "data/sphere_clipper.off";
  Mesh clipper;
  if(!read_mesh<K>(filename, clipper))
    return EXIT_FAILURE;

  std::cout << "Clipper mesh: " << num_vertices(clipper) << " nv " << num_faces(clipper) << " nf" << std::endl;

  if(PMP::is_outward_oriented(clipper))
    PMP::reverse_face_orientations(clipper);

  PMP::clip(mesh, clipper, CGAL::parameters::input_might_have_self_intersection(true));
#endif

  std::ofstream out("clipped.off");
  out.precision(17);
  out << mesh;

  std::cout << "Clipped!" << std::endl;

  return EXIT_SUCCESS;
}
