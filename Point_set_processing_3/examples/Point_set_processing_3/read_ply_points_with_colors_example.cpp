#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/property_map.h>
#include <CGAL/IO/read_ply_points.h>

#include <utility>
#include <vector>
#include <fstream>

// types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;

// Point with normal vector stored as a std::pair.
typedef std::pair<Point, Vector> Pwn;

// Color is red/green/blue array
typedef CGAL::cpp11::array<unsigned char, 3> Color;


class My_ply_interpreter
{
  std::vector<Pwn>& points;
  std::vector<Color>& colors;
    
public:
  My_ply_interpreter (std::vector<Pwn>& points,
                      std::vector<Color>& colors)
    : points (points), colors (colors)
  { }

  bool is_applicable (const std::vector<CGAL::Ply_read_number*>& readers)
  {
    std::size_t nb_property_found = 0;

    // Check if all properties are found in PLY input header
    for (std::size_t i = 0; i < readers.size (); ++ i)
      if (readers[i]->name () == "x" || // Point
          readers[i]->name () == "y" ||
          readers[i]->name () == "z" ||
          readers[i]->name () == "nx" || // Normal
          readers[i]->name () == "ny" ||
          readers[i]->name () == "nz" ||
          readers[i]->name () == "red" || // Color
          readers[i]->name () == "green" ||
          readers[i]->name () == "blue")
        nb_property_found ++;

    return (nb_property_found == 9);
  }
      
  void operator() (const std::vector<CGAL::Ply_read_number*>& readers)
  {
    FT x, y, z, nx, ny, nz;
    Color c = {{ 0, 0, 0 }};
    
    for (std::size_t i = 0; i < readers.size (); ++ i)
      if (readers[i]->name () == "x")
        readers[i]->assign (x);
      else if (readers[i]->name () == "y")
        readers[i]->assign (y);
      else if (readers[i]->name () == "z")
        readers[i]->assign (z);
      else if (readers[i]->name () == "nx")
        readers[i]->assign (nx);
      else if (readers[i]->name () == "ny")
        readers[i]->assign (ny);
      else if (readers[i]->name () == "nz")
        readers[i]->assign (nz);
      else if (readers[i]->name () == "red")
        readers[i]->assign (c[0]);
      else if (readers[i]->name () == "green")
        readers[i]->assign (c[1]);
      else if (readers[i]->name () == "blue")
        readers[i]->assign (c[2]);

    points.push_back (std::make_pair (Point (x, y, z), Vector (nx, ny, nz)));
    colors.push_back (c);
  }

};



int main(int argc, char*argv[])
{
  const char* fname = (argc>1)?argv[1]:"data/colors.ply";
    // Reads a .ply point set file with normal vectors and colors

  std::vector<Pwn> points; // store points with normals
  std::vector<Color> colors; // store colors in separate container

  My_ply_interpreter interpreter(points, colors); // init interpreter
  
  std::ifstream in(fname);
  if (!in ||
      !CGAL::read_ply_custom_points (in, interpreter, Kernel()))
    {
      std::cerr << "Error: cannot read file " << fname << std::endl;
      return EXIT_FAILURE;
    }

  for (std::size_t i = 0; i < points.size (); ++ i)
    if (colors[i][0] == 255 && colors[i][1] == 0 && colors[i][2] == 0)
      std::cerr << "Point " << points[i].first << " is red." << std::endl;
    else if (colors[i][0] == 0 && colors[i][1] == 255 && colors[i][2] == 0)
      std::cerr << "Point " << points[i].first << " is green." << std::endl;
    else if (colors[i][0] == 0 && colors[i][1] == 0 && colors[i][2] == 255)
      std::cerr << "Point " << points[i].first << " is blue." << std::endl;
  
  return EXIT_SUCCESS;
}
