#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Complex_2_in_triangulation_vertex_base_3.h>
#include <CGAL/Complex_2_in_triangulation_surface_mesh_cell_base_3.h>
#include <CGAL/Mesh_3/Complex_2_in_triangulation_cell_base_3.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Regular_triangulation_cell_base_3.h>
#include <CGAL/Implicit_surfaces_mesher_3.h>
#include <CGAL/Surface_mesher/Criteria/Standard_criteria.h>
#include <CGAL/Mesh_3/Facet_on_surface_criterion.h>

#include <CGAL/Mesh_3/Slivers_exuder.h>

#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>
#include <CGAL/IO/File_medit.h>

#include <CGAL/Surface_mesher/Oracles/Implicit_oracle.h>

#include <CGAL/Mesh_criteria_3.h>

#include <CGAL/Point_traits.h>
#include <CGAL/Point_with_surface_index.h>

#include "implicit_function.h"
#include "parameters.h"

#include "debug.h"

#include <iostream>
#include <fstream>

#ifdef CGAL_USE_QT
// for distribution
#  include <CGAL/IO/Qt_widget.h>
#  include <qapplication.h> // needed by QPainter :-(
#  include <qpixmap.h>  // qt drawing to pixmap
#endif // CGAL_USE_QT

#include <algorithm> // std::max_element()
#include <sstream>

/////////////// Types /////////////// 

struct K : public CGAL::Exact_predicates_inexact_constructions_kernel {};
typedef CGAL::Regular_triangulation_euclidean_traits_3<K> Regular_traits;
typedef CGAL::Point_with_surface_index_geom_traits<Regular_traits> My_traits;
typedef CGAL::Triangulation_vertex_base_with_info_3<bool, My_traits> Vb1;
typedef CGAL::Complex_2_in_triangulation_vertex_base_3<My_traits, Vb1> Vb;
typedef CGAL::Regular_triangulation_cell_base_3<My_traits> Cb1;
typedef CGAL::Complex_2_in_triangulation_cell_base_3<My_traits, Cb1> Cb2;
typedef CGAL::Complex_2_in_triangulation_surface_mesh_cell_base_3<My_traits, Cb2> Cb3;
typedef CGAL::Mesh_3::Complex_2_in_triangulation_cell_base_3<My_traits, Cb3> Cb;
typedef CGAL::Triangulation_data_structure_3<Vb, Cb> Tds;
typedef CGAL::Regular_triangulation_3<My_traits, Tds> Tr;
typedef Function <K::FT> Func;
typedef CGAL::Surface_mesher::Implicit_oracle<My_traits, Func> Impl_oracle;
typedef CGAL::Surface_mesher::Refine_criterion<Tr> Criterion;
typedef CGAL::Surface_mesher::Standard_criteria <Criterion > Criteria;
typedef CGAL::Mesh_criteria_3<Tr> Tets_criteria;
typedef CGAL::Implicit_surfaces_mesher_3<Tr, Impl_oracle,
                                         Criteria,
                                         Tets_criteria> Mesher;

typedef CGAL::Simple_cartesian<double> Simple_kernel;
typedef Simple_kernel::Iso_rectangle_2 Rectangle_2;
typedef Simple_kernel::Segment_2 Segment_2;
typedef Simple_kernel::Point_2 Point_2;

typedef My_traits::Point_3 Point_3;

typedef CGAL::Point_traits<Point_3> Point_traits;
typedef Point_traits::Bare_point Bare_point_3;

typedef enum { RADIUS_RATIO, ANGLE, MIN_ANGLE} Distribution_type;

/// Global variables 
typedef std::map<std::string, std::string> String_options;

String_options string_options;
std::ostream *out = 0;
std::string filename = std::string();
std::string function_name = "sphere";
bool output_to_file = false;
bool dump_distribution = false;
std::string distribution_filename;
int distribution_x = 200;
int distribution_y = 100;
int distribution_size = 50;
Distribution_type distribution_type = RADIUS_RATIO;

void init_string_options()
{
  string_options["distribution_after_filename"] = "";
  string_options["mesh_after_filename"] = "";
  string_options["initial_surface_off"] = "";
  string_options["initial_surface_ghs"] = "";
  string_options["surface_off"] = "";
  string_options["surface_ghs"] = "";
  string_options["slivers_off"] = "";
  string_options["read_initial_points"] = "";
  string_options["inrimage"] = "";
}

void usage(char *argv0, std::string error = "")
{
  if( error != "" )
    std:: cerr << "Error: " << error << std::endl;
  std::cerr << "Usage:\n  " 
            << argv0
            << " [-d <output_distribution_pixmap.png>"
            << "[-s x y] [-n N] [-t (ANGLE|RADIUS_RATIO|MIN_ANGLE)] ]"
            << " [-f function_name]"
            << " [output_file.mesh|-]\n"
            << "If output_file.mesh is '-', outputs to standard out.\n"
            << "-d  Output distribution to file"
            << " output_distribution_pixmap.png\n"
            << "-s  define pixmap size (x, y), default is (200,100)\n" 
            << "-n  define size of the distribution, default is 50\n"
            << "-t  define the type of distribution, "
            << "default is RADIUS_RATIO\n"
            << "-f  define the implicite function to use\n";
  for(String_options::iterator it = string_options.begin();
      it != string_options.end();
      ++it)
    std::cerr << "--" << it->first << " default value is "
	      << it->second << ".\n";
  for(Double_options::iterator it = double_options.begin();
      it != double_options.end();
      ++it)
    std::cerr << "--" << it->first << " default value is "
	      << it->second << ".\n";
  exit(1);
}

void parse_argv(int argc, char** argv, int extra_args = 0)
{
  if (argc >=(2 + extra_args))
    {
      std::string arg = argv[1+extra_args];
      if( arg == "-h" || arg == "--help")
        usage(argv[0]);
      else if( arg == "-d" )
	{
	  dump_distribution = true;
          if( argc < 3)
            usage(argv[0], "-d must be followed by a filename!");
          distribution_filename = argv[2 + extra_args];
          parse_argv(argc, argv, extra_args + 2);
        }
      else if( arg == "-t" )
	{
          if( argc < (3+extra_args) )
            usage(argv[0], "-t must be followed by a ANGLE or"
		  "RADIUS_RATIO or MIN_ANGLE!");
          std::string arg2 = argv[2 + extra_args];
          if( arg2 == "ANGLE" ) distribution_type = ANGLE;
          else if( arg2 == "RADIUS_RATIO" ) distribution_type = RADIUS_RATIO;
          else if( arg2 == "MIN_ANGLE" ) distribution_type = MIN_ANGLE;
          else 
            usage(argv[0], "Bad type. Should be ANGLE or RADIUS_RATIO"
		  " or MIN_ANGLE!");
          parse_argv(argc, argv, extra_args + 2);
        }
      else if( arg == "-s" )
        {
          if( argc < (4+extra_args) )
            usage(argv[0], "-s must be followed by x y!");
          {
            std::stringstream s;
            s << argv[2+extra_args];
            s >> distribution_x;
            if( !s )
              usage(argv[0], "Bad integer x!");
          }
          {
            std::stringstream s;
            s << argv[3+extra_args];
            s >> distribution_y;
            if( !s )
              usage(argv[0], "Bad integer y!");
          }
          parse_argv(argc, argv, extra_args + 3);
        }
      else if( arg == "-n" )
        {
          if( argc < (3+extra_args) )
            usage(argv[0], "-n must be followed by an integer!");
          std::stringstream s;
          s << argv[2+extra_args];
          s >> distribution_size;
          if( !s )
            usage(argv[0], "Bad integer N!");
          parse_argv(argc, argv, extra_args + 2);
        }
      else if( arg == "-f" )
        {
          if( argc < (3 + extra_args) )
            usage(argv[0], "-f must be followed by a function name!");
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
		usage(argv[0],
		      (arg + " must be followed by a double!").c_str());
	      std::stringstream s;
	      double val;
	      s << argv[extra_args + 2];
	      s >> val;
	      if( !s )
		usage(argv[0], ("Bad double after " + arg + "!").c_str());
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
                usage(argv[0],
                      (arg + " must be followed by a string!").c_str());
              std::string s = argv[extra_args + 2];
              opt_it->second = s;
              parse_argv(argc, argv, extra_args + 2);
            }
            else
              usage(argv[0], ("Invalid option: " + arg + "!").c_str());
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

#ifdef CGAL_USE_QT
template <typename Triangulation>
void output_distribution_to_png(Triangulation& tr,
                                const int number_of_classes = 100,
                                std::string filename,
                                Distribution_type type = RADIUS_RATIO)
{
  typedef Triangulation Tr;

  std::vector<double> qualities;
  double max_quality = 0;
  
  for(typename Tr::Finite_cells_iterator it = tr.finite_cells_begin();
      it != tr.finite_cells_end();
      it++)
    {
      if(it->is_in_domain())
        {
          typename Tr::Tetrahedron t = tr.tetrahedron(it);
          switch( type )
            {
            case ANGLE:
              for(int i = 0; i < 4; ++i)
                for(int j = i + 1; j < 4; j++)
                  qualities.push_back(dihedral_angle(t, i, j));
              break;
	    case MIN_ANGLE:
	      qualities.push_back(min_dihedral_angle(t));
	      break;
            default: // RADIUS_RATIO
              double q = CGAL::to_double(radius_ratio(t));
              qualities.push_back(q);
              //              max_quality = std::max(max_quality, q);
            }
        }
    }

  if( type == ANGLE ) max_quality = 180.;
  else if ( type == MIN_ANGLE ) max_quality = 90.;
  else max_quality = 1.;
  
  const int number_of_cells = qualities.size();

  std::vector<int> distribution(number_of_classes);
  
  for(int j=0;j<number_of_cells;j++)
    { // This block is (c) Pierre Alliez 2005
      int index = number_of_classes-1;
      // saturate highest value to last bin

      if(qualities[j] < max_quality)
      {
        double dratio = qualities[j]/max_quality;
        index = static_cast<int>(dratio*(double)number_of_classes);
      }
      distribution[index]++;
    }

  const int max_occurrence = *std::max_element(distribution.begin(), 
                                               distribution.end());

  CGAL::Qt_widget *widget = new CGAL::Qt_widget();
  qApp->setMainWidget(widget);
  widget->resize(distribution_x, distribution_y);
  widget->set_window(-0.2, 1.2, -0.2, 1.2, true); // x_min, x_max,
                                                  // y_min, y_max.
  widget->show();
  
  widget->lock();
  *widget << CGAL::FillColor(CGAL::Color(200, 200, 200))
          << CGAL::Color(200, 200, 200)
          << Rectangle_2(Point_2(0, 0), Point_2(1,1));
  
  if( number_of_classes == 0 ) return;
  const double width = 1.0 / number_of_classes;

  *widget << CGAL::FillColor(CGAL::BLACK);
//   *widget << Segment_2(Point_2(0., 0.), Point_2(1., 0.));
  for(int k=0;k<number_of_classes;k++)
    if(distribution[k]>0)
      {
        double height = ( distribution[k]+0. ) / max_occurrence;
        *widget << CGAL::BLACK 
		<< Rectangle_2(Point_2(k*width, 0),
			       Point_2((k+1)*width, height));
      }
    else
      *widget << CGAL::RED << Segment_2(Point_2(k*width, 0),
					Point_2((k+1)*width, 0));
  
  widget->unlock();
  if( widget->get_pixmap().save( QString(filename.c_str()),
                                 "PNG") )
    std::cerr << "Distribution saved to file " << filename
              << std::endl;
  else
    {
      std::cerr << "Error: cannot save distribution to file "
                << filename << std::endl; 
      exit(1);
    }
  qApp->exec();
}
#endif // CGAL_USE_QT

/////////////// Main function /////////////// 

int main(int argc, char **argv) {

#ifdef CGAL_USE_QT
  QApplication app (argc, argv);
#endif // CGAL_USE_QT
  init_parameters();
  
  init_string_options();

  parse_argv(argc, argv);

  init_functions();
  
  std::string inrimage = string_options["inrimage"];
  if( inrimage != "")
    {
      function_name = "generic_inrimage";
      isosurface = new CGAL::Inrimage_isosurface(inrimage.c_str(),
                                                 double_options["isovalue"]);
    }
  
  // Function
  Func F(functions[function_name]);

  // Impl_oracle (NB: parity oracle is toggled)
  Impl_oracle O (F, 
		 Point_traits().point(
                   Bare_point_3(double_options["center_x"],
                                double_options["center_y"],
                                double_options["center_z"])),
		 double_options["enclosing_sphere_radius"],
		 double_options["precision"],
		 bipolar_oracle);

  // 2D-complex in 3D-Delaunay triangulation
  Tr T;

  // Initial point sample
  std::string read_initial_points = string_options["read_initial_points"];
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
	  T.insert(p);
	  --n;
	}
    }
    CGAL_assertion( n == 0 );
  }
  else
  {
    Impl_oracle::Points initial_point_sample = 
      O.random_points (static_cast<int>(
         double_options["number_of_initial_points"]));
    
    {
      std::ofstream out("dump_of_initial_points");
      
      out << initial_point_sample.size() << "\n";
      for(Impl_oracle::Points::const_iterator it = 
	    initial_point_sample.begin();
	  it != initial_point_sample.end();
	  ++it)
	out << *it << "\n";
    }

    T.insert (initial_point_sample.begin(), initial_point_sample.end());
  }


  // Meshing criteria
  CGAL::Surface_mesher::Curvature_size_criterion<Tr> 
    c_s_crit (double_options["curvature_bound"]);
  CGAL::Surface_mesher::Uniform_size_criterion<Tr>
    u_s_crit (double_options["size_bound"]);
  CGAL::Surface_mesher::Aspect_ratio_criterion<Tr>
    a_r_crit (double_options["facets_aspect_ratio_bound"]);
  CGAL::Mesh_3::Facet_on_surface_criterion<Tr>
    facet_on_surf;

  std::vector<Criterion*> crit_vect;
  crit_vect.push_back (&c_s_crit);
  crit_vect.push_back (&u_s_crit);
  crit_vect.push_back(&a_r_crit);
  crit_vect.push_back(&facet_on_surf);
  Criteria C (crit_vect);

  Tets_criteria tets_criteria(double_options["tets_aspect_ratio_bound"],
			      double_options["tets_size_bound"]);

  std::cerr << "Initial number of points: " << T.number_of_vertices() 
            << std::endl;
  // Surface meshing
  Mesher mesher (T, O, C, tets_criteria);
  mesher.refine_surface();
  {
    std::string dump_initial_surface_filename = 
      string_options["initial_surface_off"];
    if( dump_initial_surface_filename != "" )
      {
	std::ofstream dump(dump_initial_surface_filename.c_str());
	if( dump )
	  {
	    std::cerr << "Writing initial surface to off file "
		      << dump_initial_surface_filename << "..." << std::endl;
	    output_oriented_surface_facets_to_off(dump, T);
	  }
	else
	  usage(argv[0], ("Error: cannot create " + 
			  dump_initial_surface_filename).c_str());
      }
  }
  {
    std::string dump_initial_surface_filename = 
      string_options["initial_surface_ghs"];
    if( dump_initial_surface_filename != "" )
      {
	std::ofstream dump_points((dump_initial_surface_filename +
				   ".points").c_str());
	std::ofstream dump_faces((dump_initial_surface_filename +
				  ".faces").c_str());
	if( dump_points && dump_faces )
	  {
	    std::cerr << "Writing initial surface to ghs file "
		      << dump_initial_surface_filename << "..." << std::endl;
	    output_surface_facets_to_ghs(dump_points, dump_faces, T);
	  }
	else
	  usage(argv[0], ("Error: cannot create " + 
			  dump_initial_surface_filename).c_str());
      }
  }
  mesher.refine_mesh();
//   int i = 100;
//   while(!mesher.done())
//     {
//       std::stringstream s;
//       s << "out." << i++ << ".mesh";
//       std::cerr << s.str() << std::endl;
//       std::ofstream os3(s.str().c_str());
//       output_pslg_to_medit(os3, T);
//       mesher.step_by_step();
//     }
//   exit(0);
  std::cerr << "Final number of points: " << T.number_of_vertices() 
            << std::endl;

  if( output_to_file )
    {
      std::ofstream file;
      if( filename == "-" ) {
	std::cerr << "Writing to standard out..." << std::endl;
	out = &std::cout;
      }
      else {
	file.open(filename.c_str());
      
	if( file ) {
	  std::cerr << "Writing to file " << filename << "..." << std::endl;
	  out = &file;
	}
	else {
	  std::cerr << "Error: cannot open " << filename << std::endl;
	  usage(argv[0]);
	}
      }
    
      // Output
      output_pslg_to_medit(*out, T);
    }

#ifdef CGAL_USE_QT
  if( dump_distribution )
    output_distribution_to_png(T, distribution_size, 
			       distribution_filename, distribution_type);
#endif
  
  CGAL::Mesh_3::Slivers_exuder<Tr> exuder(T, double_options["pumping_bound"]);
  int number_of_pump = static_cast<int>(double_options["number_of_pump"]);
  for(int i = 0; i < number_of_pump; ++i)
  {
    exuder.init();
    exuder.pump_vertices();
  }
  
#ifdef CGAL_USE_QT
  std::string distribution_after_filename = 
    string_options["distribution_after_filename"];
  if( distribution_after_filename != "" )
    output_distribution_to_png(T, distribution_size,
			       distribution_after_filename, distribution_type);
#endif

  std::string mesh_after_filename = string_options["mesh_after_filename"];
  if( mesh_after_filename != "" )
  {
    std::ofstream file(mesh_after_filename.c_str());
    if( file ) {
      std::cerr << "Writing to file " << mesh_after_filename
		<< "..." << std::endl;
      output_pslg_to_medit(file, T);
    }
    else
      usage(argv[0] , ("Error: cannot create " + mesh_after_filename).c_str());
  }

  {
    std::string dump_final_surface_filename = string_options["surface_off"];
    if( dump_final_surface_filename != "" )
      {
	std::ofstream dump(dump_final_surface_filename.c_str());
	if( dump ) {
	  std::cerr << "Writing final surface to off file "
		    << dump_final_surface_filename << "..." << std::endl;
	  output_oriented_surface_facets_to_off(dump, T);
	}
	else
	  usage(argv[0], ("Error: cannot create " + 
			  dump_final_surface_filename).c_str());
      }
  }
  {
    std::string dump_final_surface_filename = string_options["surface_ghs"];
    if( dump_final_surface_filename != "" )
      {
	std::ofstream dump_points((dump_final_surface_filename +
				   ".points").c_str());
	std::ofstream dump_faces((dump_final_surface_filename + 
				  ".faces").c_str());
	if( dump_points && dump_faces ) {
	  std::cerr << "Writing final surface to ghs file "
		    << dump_final_surface_filename << "..." << std::endl;
	  output_surface_facets_to_ghs(dump_points, dump_faces, T);
	}
	else
	  usage(argv[0], ("Error: cannot create " + 
			  dump_final_surface_filename).c_str());
      }
  }
  
  std::string dump_slivers_filename = string_options["slivers_off"];
  if( dump_slivers_filename != "" )
  {
    std::ofstream dump(dump_slivers_filename.c_str());
    if( dump ) {
      std::cerr << "Writing slivers to off file "
		<< dump_slivers_filename << "..." << std::endl;
      output_slivers_to_off(dump, T, double_options["sliver_test"]);
    }
    else
      usage(argv[0], ("Error: cannot create " + 
		      dump_slivers_filename).c_str());
  }

  std::cerr << " done\n";
}
