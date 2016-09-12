#include "debug_macros.h"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Robust_circumcenter_traits_3.h>
#include <CGAL/Delaunay_triangulation_3.h>

#include <CGAL/Polyhedral_surface_3.h>

#include <CGAL/Surface_mesher/Standard_criteria.h>

// #define CGAL_C2T3_USE_POLYHEDRON
#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>

#include <CGAL/Surface_mesher/Point_surface_indices_oracle_visitor.h>

#include <CGAL/make_surface_mesh.h>
#include <CGAL/Surface_mesh_default_edges_criteria_3.h>

#include "parameters.h"

#include <CGAL/Timer.h>

#include <iostream>
#include <fstream>

#include <sstream>

#include <boost/tuple/tuple.hpp> // boost::tie

using boost::tie;

/////////////// Types ///////////////

struct K2 : public CGAL::Exact_predicates_inexact_constructions_kernel {};
typedef CGAL::Robust_circumcenter_traits_3<K2>  K;
typedef CGAL::Complex_2_in_triangulation_vertex_base_3<K> Vb;
typedef CGAL::Surface_mesh_cell_base_3<K> Cb;
typedef CGAL::Triangulation_data_structure_3<Vb, Cb> Tds;
typedef CGAL::Delaunay_triangulation_3<K, Tds> Tr;

typedef K::Point_3 Point_3;
typedef K::Sphere_3 Sphere_3;
typedef K::FT FT;

typedef CGAL::Polyhedral_surface_3<K, CGAL::Surface_mesher::Has_edges> Surface;

struct Edge_info {
  bool is_cached;
  Point_3 lineic_center;
};

typedef CGAL::Complex_2_in_triangulation_3<Tr, Edge_info> C2t3;

typedef CGAL::Surface_mesher::Polyhedral_oracle<Surface> Surface_mesh_traits;

typedef CGAL::Surface_mesher::Refine_criterion<Tr> Criterion;
typedef CGAL::Surface_mesher::Standard_criteria <Criterion > Criteria;
typedef CGAL::Surface_mesh_default_edges_criteria_3<Tr> Edges_criteria;

typedef Surface_mesh_traits::Construct_initial_points Initial_points;

typedef CGAL::Simple_cartesian<double> Simple_kernel;
typedef Simple_kernel::Iso_rectangle_2 Rectangle_2;
typedef Simple_kernel::Segment_2 Segment_2;
typedef Simple_kernel::Point_2 Point_2;

/// Global variables
std::ostream *out = 0;
std::string filename = std::string();
std::string function_name = "";
char* argv0 = "";

void usage(std::string error = "")
{
  if( error != "" )
    std:: cerr << "Error: " << error << std::endl;
  std::cerr << "Usage:\n  "
            << argv0
            << " -f function_name"
            << " [output_file.off|-]\n"
            << "If output_file.off is '-', outputs to standard out.\n"
            << "-f define the OFF file to remesh.\n";
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

std::pair<std::ostream*, const bool>
open_file_for_writing(std::string filename,
                      std::string display_string = "Writing to ")
{
  if( filename != "")
  {
    if( filename == "-" )
    {
      std::cerr << display_string << "standard out...\n";
      return std::make_pair(&std::cout, false);
    }
    else
    {
      std::ofstream* result = new std::ofstream(filename.c_str());
      if( *result )
      {
        std::cerr << display_string << "file " << filename << "...\n";
        return std::make_pair(result, true);
      }
      else
      {
        delete result;
        std::cerr << "Error: cannot create " << filename << "\n";
        usage();
        return std::pair<std::ostream*, bool>(0, false);
      }
    }
  }
  else
    return std::pair<std::ostream*, bool>(0, false);
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
	  filename = argv[1+extra_args];
	  parse_argv(argc, argv, extra_args + 1);
	}
    }
}

/////////////// Main function ///////////////

int main(int argc, char **argv) {
  argv0 = argv[0];

  usage_ptr = &usage;

  init_parameters();

  parse_argv(argc, argv);

  if( function_name == "" )
    usage("Empty input file name");

  std::ifstream surface_ifs(function_name.c_str());
  Surface surface(surface_ifs,
		  FT(get_double_option("sharp_edge_cosine_bound")));
  surface_ifs.close();

  std::cerr << "Surface bounding box: " << surface.bbox() << "\n";

  // 2D-complex in 3D-Delaunay triangulation
  Tr tr;
  C2t3 c2t3(tr);

  CGAL::Timer timer;

  bool need_delete = false;
  std::ostream* out = 0;

  Surface_mesh_traits surface_mesh_traits;

  // Initial point sample
  std::string read_initial_points = get_string_option("read_initial_points");
  if( read_initial_points != "")
  {
    std::ifstream in( read_initial_points.c_str() );
    int n;
    in >> n;
    CGAL_assertion(in);
    while( !in.eof() )
      {
	Point_3 p;
	if(in >> p)
	  {
	    tr.insert(p);
	    --n;
	  }
      }
    CGAL_assertion( n == 0 );
    double_options["number_of_initial_points"] = 0;
  }
  else
  {
    const int number_of_initial_points =
      static_cast<int>(get_double_option("number_of_initial_points"));

    std::vector<Point_3> initial_point_sample;
    initial_point_sample.reserve(number_of_initial_points);

    Initial_points get_initial_points =
      surface_mesh_traits.construct_initial_points_object();

    get_initial_points(surface,
                       std::back_inserter(initial_point_sample),
                       number_of_initial_points);

    tie(out, need_delete) =
      open_file_for_writing(get_string_option("dump_of_initial_points"),
                            "Writing initial points to ");
    if( out )
    {
      *out << initial_point_sample.size() << "\n";
      for(std::vector<Point_3>::const_iterator it =
            initial_point_sample.begin();
          it != initial_point_sample.end();
          ++it)
        *out << *it <<"\n";
      if(need_delete)
        delete out;
    }
    tr.insert (initial_point_sample.begin(), initial_point_sample.end());
  }

  std::cerr << "Inserting bbox points... ";
  {
    const double
      xm = (surface.bbox().xmax()-surface.bbox().xmin())/2,
      ym = (surface.bbox().ymax()-surface.bbox().ymin())/2,
      zm = (surface.bbox().zmax()-surface.bbox().zmin())/2;
    

    K::Iso_cuboid_3 bbox(surface.bbox().xmin()-xm,
			 surface.bbox().ymin()-ym,
			 surface.bbox().zmin()-zm,
			 surface.bbox().xmax()+xm,
			 surface.bbox().ymax()+ym,
			 surface.bbox().zmax()+zm);
    for(int i =0; i < 8; ++i)
      tr.insert(bbox.vertex(i));
  }
  std::cerr << "done\n";

  // Meshing criteria
  CGAL::Surface_mesher::Curvature_size_criterion<Tr>
    curvature_size_criterion (get_double_option("distance_bound"));
  CGAL::Surface_mesher::Uniform_size_criterion<Tr>
    uniform_size_criterion (get_double_option("radius_bound"));
  CGAL::Surface_mesher::Aspect_ratio_criterion<Tr>
    aspect_ratio_criterion (get_double_option("angle_bound"));

  std::vector<Criterion*> criterion_vector;
  criterion_vector.push_back(&aspect_ratio_criterion);
  criterion_vector.push_back(&uniform_size_criterion);
  criterion_vector.push_back(&curvature_size_criterion);
  Criteria criteria (criterion_vector);

  Edges_criteria edges_criteria(get_double_option("edges_radius_bound"),
				get_double_option("edges_distance_bound"));
				

  std::cerr << "\nInitial number of points: " << tr.number_of_vertices()
            << std::endl;


  // Surface meshing

  timer.start();
  CGAL::make_piecewise_smooth_surface_mesh(c2t3,
					   surface,
					   criteria,
					   edges_criteria,
					   CGAL::Non_manifold_tag(),
// 					   CGAL::Manifold_with_boundary_tag(),
					   0);
  timer.stop();
  std::cerr << ::boost::format("\nFinal number of points: %1%"
			       "\nFinal number of facets: %2%"
			       "\nFinal number of edges: %3%"
			       "\nFinal number of marked edges: %4%"
			       "\nTotal time: %5%\n")
    % tr.number_of_vertices()
    % c2t3.number_of_facets()
    % c2t3.number_of_edges()
    % c2t3.number_of_marked_edges()
    % timer.time();

  tie(out, need_delete) =
    open_file_for_writing(filename,
                          "Writing finale surface off to ");
  if( out )
  {
    CGAL::output_surface_facets_to_off(*out, c2t3);
    if(need_delete)
      delete out;
  }
  std::cerr << " done\n";

#ifdef CGAL_SURFACE_MESHER_TEST_OPTIONS
  check_all_options_have_been_used();
#endif
}
