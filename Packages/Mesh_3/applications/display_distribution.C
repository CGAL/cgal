#include <CGAL/Simple_cartesian.h>
#include <iostream>
#include <fstream>
#include <sstream>

#include <CGAL/Mesh_3/Slivers_exuder.h>

#ifndef CGAL_USE_QT
int main()
{
  std::cerr << "This tool requires Qt.\n";
  exit(1);
}
#else
#include <CGAL/IO/Qt_widget.h>
#include <qpainter.h>
#include <qapplication.h>

struct K : public CGAL::Simple_cartesian<double> {};
typedef K::Point_3 Point_3;
typedef K::Tetrahedron_3 Tetrahedron_3;
typedef K::Point_2 Point_2;
typedef K::Segment_2 Segment_2;
typedef CGAL::Iso_rectangle_2<K> Rectangle_2;

/// Global variables 
typedef std::map<std::string, std::string> String_options;
typedef std::map<std::string, double> Double_options;

String_options string_options;
Double_options double_options;

void init_options()
{
  string_options["tets"] = "";
  string_options["noboite"] = "";
  string_options["mesh"] = "";
  double_options["echelle"] = 1;
}

void usage(char *argv0, std::string error = "")
{
  if( error != "" )
    std:: cerr << "Error: " << error << std::endl;
  std::cerr << "Usage:\n  " 
            << argv0
            << " [(--tets|--noboite|--mesh) file]"
            << std::endl;
  exit(1);
}

void parse_argv(int argc, char** argv, int extra_args = 0)
{
  if (argc >=(2 + extra_args))
    {
      std::string arg = argv[1+extra_args];
      if( arg.substr(0, 2) == "--" )
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
    }
} // end parse_argv

void output_distribution_to_png(std::vector<double>& elements,
                                double max,
                                const int number_of_classes,
                                std::string filename)
{
  const int number_of_cells = elements.size();

  std::vector<int> distribution(number_of_classes);
  
  for(int j=0;j<number_of_cells;j++)
    { // This block is (c) Pierre Alliez 2005
      int index = number_of_classes-1;
      // saturate highest value to last bin

      if(elements[j] < max)
	{
	  double dratio = elements[j]/max;
	  index = static_cast<int>(dratio*(double)number_of_classes);
	}
      distribution[index]++;
    }

//   const int max_occurrence = *std::max_element(distribution.begin(), 
// 					       distribution.end());

  CGAL::Qt_widget *widget = new CGAL::Qt_widget();
  qApp->setMainWidget(widget);
  widget->resize(400, 400);
  widget->set_window(0, 1, 0, 1, true); // x_min, x_max,
                                        // y_min, y_max.
  widget->show();
  
  widget->lock();
  *widget << CGAL::FillColor(CGAL::Color(200, 200, 200))
	  << CGAL::Color(200, 200, 200)
	  << Rectangle_2(Point_2(0, 0), Point_2(1,1));
  
  if( number_of_classes == 0 ) return;
  const double width = 1.0 / number_of_classes;

  const double echelle = double_options["echelle"];

  *widget << CGAL::FillColor(CGAL::BLACK);
  //   *widget << Segment_2(Point_2(0., 0.), Point_2(1., 0.));
  for(int k=0;k<number_of_classes;k++)
    if(distribution[k]>0)
      {
	double height = ( (distribution[k]+0.)/number_of_cells ) * echelle;
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

bool failed(const char* error)
{
  std::cerr << error << std::endl;
  return false;
}

bool read_tets(std::vector<double>& elements, std::istream& in)
{
  // read header
  int nb_vertices = 0;
  int nb_tets = 0;
  std::string head;

  in >> nb_vertices >> head;
  if( !in || head != "vertices" )
    return false;

  in >> nb_tets >> head;
  if( !in || head != "tets" )
    return false;

  std::vector<Point_3> points;
  points.reserve(nb_vertices);

  // read points
  for(int i=0;i<nb_vertices;i++)
    {
      float x,y,z;
      in >> x >> y >> z;
      if( !in )
	return false;
      points.push_back(Point_3(x,y,z));
    }

  // read tets
  for(int i=0;i<nb_tets;i++)
    {
      int dummy, i0,i1,i2,i3;
      in >> dummy >> i0 >> i1 >> i2 >> i3;
      if( dummy != 4 || !in )
	return false;
      Tetrahedron_3 t = Tetrahedron_3(points[i0],
				      points[i1],
				      points[i2],
				      points[i3]);
      elements.push_back(radius_ratio(t));
    }
  return true;
}

bool read_mesh(std::vector<double>& elements, std::istream& in)
{
  // Header.
  std::string head;
  int version;
  in >> head >> version;
  if( head != "MeshVersionFormatted" ||
      version != 1 ||
      ! in)
    return failed("read_mesh: bad version");
  
  int dimension;
  in >> head >> dimension;
  if( head != "Dimension" ||
      dimension!= 3 ||
      ! in)
    return failed("read_mesh: bad dimension");

  // Vertices
  int number_of_vertices;
  in >> head >> number_of_vertices;
  if( head != "Vertices" || ! in )
    return failed("read_mesh: bad file (missing Vertices)");
  std::vector<Point_3> points;
  points.reserve(number_of_vertices);
  for(int i = 0; i < number_of_vertices; ++i)
  {
    int dummy_i;
    in >> points[i] >> dummy_i;
    if( !in )
      return failed("read_mesh: bad file (reading of vertices)");
  }

  // Facets
  int number_of_facets_on_surface;
  in >> head >> number_of_facets_on_surface;
  if( !in || head != "Triangles" )
    return failed("read_mesh: bad file (missing Triangles)");
  for(int i = 0; i < 4 * number_of_facets_on_surface; ++i)
  {
    double dummy;
    in >> dummy;
  }
  
  // Tetrahedra
  int number_of_cells;
  in >> head >> number_of_cells;
  if( !in || head != "Tetrahedra")
    return failed("read_mesh: bad file (missing Tetrahedra)");
  for(int i = 0; i < number_of_cells; ++i)
  {
    int i0, i1, i2, i3, dummy;
    in >> i0 >> i1 >> i2 >> i3 >> dummy;
    if( !in )
      return failed("read_mesh: bad file (reading of cells)");
    const Tetrahedron_3 t = Tetrahedron_3(points[i0-1],
                                          points[i1-1],
                                          points[i2-1],
                                          points[i3-1]);
    elements.push_back(radius_ratio(t));
  }
  in >> head;
  if ( !in || head != "End")
    return failed("read_mesh: bad file (missing End)");
  else 
    return true;
}

bool read_noboite(std::vector<double>& elements, std::istream& in)
{
  int nb_vertices = 0;
  int nb_tets = 0;
  int nb_input_points;
  int dummy;
  
  in >> nb_tets >> nb_vertices >> nb_input_points 
     >> dummy
     >> dummy
     >> dummy
     >> dummy
     >> dummy
     >> dummy
     >> dummy
     >> dummy
     >> dummy
     >> dummy
     >> dummy
     >> dummy
     >> dummy
     >> dummy;

  if( !in ) 
    return false;
  
  std::vector<Point_3> points;
  points.reserve(nb_vertices);
  
  elements.clear();
  elements.reserve(nb_tets);
  
  // read tets
  std::vector<int> tets;
  tets.reserve(4 * nb_tets);
  for(int i=0;i<nb_tets;i++)
  {
    int i0,i1,i2,i3;
    in >> i0 >> i1 >> i2 >> i3;
    if( !in )
      return false;
    tets.push_back(i0-1);
    tets.push_back(i1-1);
    tets.push_back(i2-1);
    tets.push_back(i3-1);
  }

    
  // read points
  for(int i=0;i<nb_vertices;i++)
  {
    double x,y,z;
    in >> x >> y >> z;
    if( !in ) 
      return false;
    points.push_back(Point_3(x,y,z));
  }
  
  // compute elements
  for(int i = 0; i < 4 * nb_tets; i += 4)
  {
    const Tetrahedron_3 t = Tetrahedron_3(points[tets[i]],
                                          points[tets[i+1]],
                                          points[tets[i+2]],
                                          points[tets[i+3]]);
    elements.push_back(radius_ratio(t));
  }
  return true;
}

int main(int argc, char** argv)
{
  QApplication app(argc, argv);
  init_options();
  parse_argv(argc, argv, 0);
  
  bool ghs = false;
  bool tets = false;
  bool mesh = false;
  std::string filename = "";
  std::string ghs_filename = string_options["noboite"];
  std::string tets_filename = string_options["tets"];
  std::string mesh_filename = string_options["mesh"];
  
  if(ghs_filename != "")
  {
    ghs = true;
    filename = ghs_filename;
  }
  if(mesh_filename != "")
  {
    mesh = true;
    filename = mesh_filename;
  }
  if(tets_filename != "")
  {
    tets = true;
    filename = tets_filename;
  }
  
  std::vector<double> elements;
  std::ifstream in_file(filename.c_str());
  if(tets)
    tets = read_tets(elements, in_file);
  if(mesh)
    mesh = read_mesh(elements, in_file);
  if(ghs)
    ghs = read_noboite(elements, in_file);
  std::stringstream png_name;
  png_name << filename << "-echelle-" << double_options["echelle"]
	   << ".png";
  if(tets || mesh || ghs)
    output_distribution_to_png(elements, 1., 100, png_name.str());
  else
    usage(argv[0], "Error: cannot read file " + filename + "!");
}
#endif // CGAL_USE_QT
