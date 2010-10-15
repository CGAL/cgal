#include "debug_macros.h"
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Gray_level_image_3.h>
#include <CGAL/Implicit_surface_3.h>

#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>

#include "implicit_functions.h"
#include "parameters.h"

#include <iostream>
#include <fstream>
#include <cassert>


/////////////// Types ///////////////

typedef CGAL::Surface_mesh_default_triangulation_3 Tr;
typedef CGAL::Surface_mesh_complex_2_in_triangulation_3<Tr> C2t3;

typedef Tr::Geom_traits Gt;

typedef Gt::Point_3 Point_3;
typedef Gt::Sphere_3 Sphere_3;
typedef Gt::FT FT;

typedef CGAL::Implicit_function_wrapper<FT, Point_3> Implicit_function_wrapper;

typedef CGAL::Implicit_surface_3<Gt, Implicit_function_wrapper> Surface;

typedef CGAL::Gray_level_image_3<FT, Point_3> Gray_image;

Gray_image* isosurface = 0;

double generic_inrimage_function(double x, double y, double z)
{
  assert(isosurface != 0);

  return (*isosurface)(Point_3(x, y, z));
}

/// Global variables
std::ostream *out = 0;
std::string filename = std::string();
std::string function_name = "sphere";
bool output_to_file = false;
char* argv0 = "";

void usage(std::string error = "")
{
  if( error != "" )
    std:: cerr << "Error: " << error << std::endl;
  std::cerr << "Usage:\n  "
            << argv0
            << " [-f function_name]"
            << " [output_file.mesh|-]\n"
            << "If output_file.mesh is '-', outputs to standard out.\n"
            << "-f  define the implicite function to use\n";
  for(String_options::iterator it = string_options.begin();
      it != string_options.end();
      ++it)
    std::cerr << "--" << it->first << " default value is \""
	      << it->second << "\".\n";
  for(Double_options::iterator it = double_options.begin();
      it != double_options.end();
      ++it)
    std::cerr << "--" << it->first << " default value is "
	      << it->second << ".\n";
  exit(EXIT_FAILURE);
}

void parse_argv(int argc, char** argv, int extra_args = 0)
{
  if (argc >=(2 + extra_args))
    {
      std::string arg = argv[1+extra_args];
      if( arg == "-h" || arg == "--help")
        usage();
      else if( arg == "-f" )
        {
          if( argc < (3 + extra_args) )
            usage("-f must be followed by a function name!");
          function_name = argv[2 + extra_args];
          parse_argv(argc, argv, extra_args + 2);
        }
      else if( arg.substr(0, 2) == "--" )
	{
	  Double_options::iterator opt_it =
	    double_options.find(arg.substr(2, arg.length()-2));
	  if( opt_it != double_options.end() )
	    {
	      if( argc < (3 + extra_args) )
		usage((arg + " must be followed by a double!").c_str());
	      std::stringstream s;
	      double val;
	      s << argv[extra_args + 2];
	      s >> val;
	      if( !s )
		usage(("Bad double after " + arg + "!").c_str());
	      opt_it->second = val;
	      parse_argv(argc, argv, extra_args + 2);
	    }
	  else
          {
            String_options::iterator opt_it =
                string_options.find(arg.substr(2, arg.length()-2));
            if( opt_it != string_options.end() )
            {
              if( argc < (3 + extra_args) )
                usage((arg + " must be followed by a string!").c_str());
              std::string s = argv[extra_args + 2];
              opt_it->second = s;
              parse_argv(argc, argv, extra_args + 2);
            }
            else
              usage(("Invalid option " + arg).c_str());
          }
	}
      else
	{
	  output_to_file = true;
	  filename = argv[1+extra_args];
	  parse_argv(argc, argv, extra_args + 1);
	}
    }
}

/////////////// Main function ///////////////

int main(int argc, char **argv) {
  argv0 = argv[0];

  init_parameters();
  functions["generic_inrimage"] = &generic_inrimage_function;

  parse_argv(argc, argv);

  std::string inrimage = string_options["inrimage"];
  if( inrimage != "")
    {
      function_name = "generic_inrimage";
      isosurface = new Gray_image(inrimage.c_str(),
                                  double_options["isovalue"]);
    }

  // Function
  FT sphere_radius = double_options["enclosing_sphere_radius"];

  Surface surface(functions[function_name],
                  Sphere_3(Point_3(double_options["center_x"],
                                   double_options["center_y"],
                                   double_options["center_z"]),
                           sphere_radius*sphere_radius),
                  double_options["precision"]);

  // 2D-complex in 3D-Delaunay triangulation
  Tr tr;
  C2t3 c2t3(tr);

  // Initial point sample
  std::string read_initial_points = string_options["read_initial_points"];
  if( read_initial_points != "")
  {
    std::ifstream in( read_initial_points.c_str() );
    int n;
    in >> n;
    if(in) 
    {
    if(out)
      std::cerr << "Reading initial points to file \"" 
                << read_initial_points << "\"...\n";
      double_options["number_of_initial_points"] = 0;
      while( in.good() && !in.eof() )
      {
	Point_3 p;
	if(in >> p)
        {
          tr.insert(p);
          --n;
        }
      }
      CGAL_assertion( n == 0 );
    }
    else
    {
      std::cerr << "Cannot read initial points from file \""
                << read_initial_points << "\"!\n";
    }
  }
  if(double_options["number_of_initial_points"] != 0)
  {
    const int number_of_initial_points =
      static_cast<int>(double_options["number_of_initial_points"]);

    std::vector<Point_3> initial_point_sample;
    initial_point_sample.reserve(number_of_initial_points);

    typedef CGAL::Surface_mesh_traits_generator_3<Surface>::type Oracle;
    Oracle::Construct_initial_points initial_points =
      Oracle().construct_initial_points_object();
    
    initial_points(surface, std::back_inserter(initial_point_sample),
                   number_of_initial_points);

    std::ofstream out(string_options["dump_of_initial_points"].c_str());
    if(out)
      std::cerr << "Dumping initial points to file \"" 
                << string_options["dump_of_initial_points"] << "\"...\n";

    out << initial_point_sample.size() << "\n";
    for(std::vector<Point_3>::const_iterator it =
	  initial_point_sample.begin();
	it != initial_point_sample.end();
	++it)
      out << *it <<"\n";

    tr.insert (initial_point_sample.begin(), initial_point_sample.end());
  }

  // defining meshing criteria
  CGAL::Surface_mesh_default_criteria_3<Tr>
    criteria(double_options["angle_bound"],   // angular bound
             double_options["radius_bound"],  // radius bound
             double_options["distance_bound"]); //distance bound

  std::cerr << "Number of initial points: "
            << tr.number_of_vertices() << std::endl;

  // Surface meshing
  CGAL::make_surface_mesh(c2t3, surface, criteria,
                          CGAL::Manifold_with_boundary_tag(),
                          0); // zero extra initial points

  std::cerr << "\nNumber of points after refine_surface(): "
            << tr.number_of_vertices() << std::endl;

  std::cerr << "\nWriting " << filename.c_str() << "...\n";
  std::ofstream out(filename.c_str());
  if( !out )
    std::cerr << "Cannot open " << filename.c_str() << "!\n";
  else
  {
    CGAL::output_surface_facets_to_off(out, c2t3);
    std::cerr << "done.\n";
  }
}
