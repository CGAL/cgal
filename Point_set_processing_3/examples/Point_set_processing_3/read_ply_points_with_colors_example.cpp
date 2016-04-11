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

// Custom interpreter that reads points, normals and colors and stores
// them in the appropriate container
class My_ply_interpreter
{
  std::vector<Pwn>& points;
  std::vector<Color>& colors;
    
public:
  My_ply_interpreter (std::vector<Pwn>& points,
                      std::vector<Color>& colors)
    : points (points), colors (colors)
  { }

  // Init and test if input file contains the right properties
  bool is_applicable (CGAL::Ply_reader& reader)
  {
    return reader.does_tag_exist<FT> ("x")
      && reader.does_tag_exist<FT> ("y")
      && reader.does_tag_exist<FT> ("z")
      && reader.does_tag_exist<FT> ("nx")
      && reader.does_tag_exist<FT> ("ny")
      && reader.does_tag_exist<FT> ("nz")
      && reader.does_tag_exist<unsigned char> ("red")
      && reader.does_tag_exist<unsigned char> ("green")
      && reader.does_tag_exist<unsigned char> ("blue");
  }

  // Describes how to process one line (= one point object)
  void process_line (CGAL::Ply_reader& reader)
  {
    FT x = (FT)0., y = (FT)0., z = (FT)0.,
      nx = (FT)0., ny = (FT)0., nz = (FT)0.;
    Color c = {{ 0, 0, 0 }};

    reader.assign (x, "x");
    reader.assign (y, "y");
    reader.assign (z, "z");
    reader.assign (nx, "nx");
    reader.assign (ny, "ny");
    reader.assign (nz, "nz");
    reader.assign (c[0], "red");
    reader.assign (c[1], "green");
    reader.assign (c[2], "blue");

    points.push_back (std::make_pair (Point (x, y, z), Vector (nx, ny, nz)));
    colors.push_back (c);
  }

};



int main(int argc, char*argv[])
{
  const char* fname = (argc>1) ? argv[1] : "data/colors.ply";
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

  // Display points with pure r/g/b colors
  for (std::size_t i = 0; i < points.size (); ++ i)
    if (colors[i][0] == 255 && colors[i][1] == 0 && colors[i][2] == 0)
      std::cerr << "Point " << points[i].first << " is red." << std::endl;
    else if (colors[i][0] == 0 && colors[i][1] == 255 && colors[i][2] == 0)
      std::cerr << "Point " << points[i].first << " is green." << std::endl;
    else if (colors[i][0] == 0 && colors[i][1] == 0 && colors[i][2] == 255)
      std::cerr << "Point " << points[i].first << " is blue." << std::endl;
  
  return EXIT_SUCCESS;
}
