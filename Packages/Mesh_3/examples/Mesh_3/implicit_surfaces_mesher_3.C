#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Complex_2_in_triangulation_surface_mesh_cell_base_3.h>
#include <CGAL/Mesh_3/Complex_2_in_triangulation_cell_base_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Implicit_surfaces_mesher_3.h>
#include <CGAL/Chew_4_surfaces/Criteria/Standard_criteria.h>
#include <CGAL/IO/Complex_2_in_triangulation_3_file_writer.h>

#include <CGAL/IO/File_medit.h>

#include <CGAL/Chew_4_surfaces/Oracles/Implicit_oracle.h>

#include <CGAL/Mesh_criteria_3.h>

#include "implicit_function.h"
#include "parameters.h"

#include "debug.h"

#include <iostream>
#include <fstream>

// for distribution
#include <CGAL/IO/Qt_widget.h>
#include <qapplication.h> // needed by QPainter :-(
#include <algorithm> // std::max_element()
#include <qpixmap.h>  // qt drawing to pixmap
#include <sstream>

/////////////// Types /////////////// 

struct K : public CGAL::Exact_predicates_inexact_constructions_kernel {};
typedef CGAL::Triangulation_vertex_base_with_info_3<bool, K> Vb;
typedef CGAL::Complex_2_in_triangulation_surface_mesh_cell_base_3<K> CCb;
typedef CGAL::Mesh_3::Complex_2_in_triangulation_cell_base_3<K, CCb> Cb;
typedef CGAL::Triangulation_data_structure_3<Vb, Cb> Tds;
typedef CGAL::Delaunay_triangulation_3<K, Tds> Del;
typedef Function <K::FT> Func;
typedef CGAL::Chew_4_surfaces::Implicit_oracle<K, Func> Oracle;
typedef CGAL::Chew_4_surfaces::Refine_criterion<Del> Criterion;
typedef CGAL::Chew_4_surfaces::Standard_criteria <Criterion > Criteria;
typedef CGAL::Mesh_criteria_3<Del> Tets_criteria;
typedef CGAL::Implicit_surfaces_mesher_3<Del, Oracle,
                                         Criteria,
                                         Tets_criteria> Mesher;

typedef CGAL::Simple_cartesian<double> Simple_kernel;
typedef Simple_kernel::Iso_rectangle_2 Rectangle;
typedef Simple_kernel::Segment_2 Segment;
typedef Simple_kernel::Point_2 Point;

/// Global variables 
std::ostream *out = 0;
std::string filename = std::string();
bool output_to_file = false;
bool dump_distribution = false;
std::string distribution_filename;
int distribution_x = 200;
int distribution_y = 100;
int distribution_size = 50;

void usage(char *argv0, std::string error = "")
{
  if( error != "" )
    std:: cerr << error << std::endl;
  std::cerr << "Usage:\n  " 
            << argv0
            << " [-d <output_distribution_pixmap.png> [-s x y] [-n N] ]"
            << " [output_file.mesh|-]\n"
            << "If output_file.mesh is '-', outputs to standard out.\n"
            << "-d  Output distribution to file"
            << " output_distribution_pixmap.png\n"
            << "-s  define pixmap size (x, y), default is (200,100)\n" 
            << "-n  define size of the distribution, default is 50\n"
            << std::endl;
  exit(1);
}

void parse_argv(int argc, char** argv)
{
  int extra_args= 0;
  if (argc >=2)
    {
      std::string arg = argv[1];
      if( arg == "-d" )
	{
	  dump_distribution = true;
          if( argc < 3)
            usage(argv[0]);
          distribution_filename = argv[2];
          extra_args += 2;
        }
      if( argc >= (2+extra_args) )
        {
          arg = argv[1+extra_args];
          if( arg == "-s" )
            {
              if( argc < (4+extra_args) )
                usage(argv[0], "-s x y!");
              {
                std::stringstream s;
                s << argv[2+extra_args];
                s >> distribution_x;
                if( !s )
                  usage(argv[0], "x!");
              }
              {
                std::stringstream s;
                s << argv[3+extra_args];
                s >> distribution_y;
                if( !s )
                  usage(argv[0], "y!");
              }
              extra_args += 3;
            }
        }
      if( argc >= (2+extra_args) )
        {
          arg = argv[1+extra_args];
          if( arg == "-n" )
            {
              if( argc < (3+extra_args) )
                usage(argv[0]);
              std::stringstream s;
              s << argv[2+extra_args];
              s >> distribution_size;
              if( !s )
                usage(argv[0]);
              extra_args += 2;
            }
        }
    }

  if (argc >= (2 + extra_args)) {
    output_to_file = true;
    filename = argv[1+extra_args];
  }
}

template <typename K>
double
area(const CGAL::Tetrahedron_3<K>& t)
{// This function is (c) Pierre Alliez 2005
  typedef typename K::Triangle_3 Triangle;
    
  Triangle t1 = Triangle(t[0],t[1],t[2]);
  Triangle t2 = Triangle(t[0],t[1],t[3]);
  Triangle t3 = Triangle(t[2],t[1],t[3]);
  Triangle t4 = Triangle(t[2],t[0],t[3]);
  double a1 = std::sqrt(CGAL_NTS to_double(t1.squared_area()));
  double a2 = std::sqrt(CGAL_NTS to_double(t2.squared_area()));
  double a3 = std::sqrt(CGAL_NTS to_double(t3.squared_area()));
  double a4 = std::sqrt(CGAL_NTS to_double(t4.squared_area()));
  return a1 + a2 + a3 + a4;
}

template <typename K>
double
circumradius(const CGAL::Tetrahedron_3<K>& t)
{
  typename K::Point_3 center = circumcenter( t.vertex(0),
                                             t.vertex(1),
                                             t.vertex(2),
                                             t.vertex(3));
  return CGAL::sqrt(CGAL::to_double( squared_distance( center, t.vertex(0))));
}

template <typename K>
double
radius_ratio(const typename CGAL::Tetrahedron_3<K>& t)
{ // This function is (c) Pierre Alliez 2005
  typename K::FT inradius = 3 * std::abs(t.volume()) / area(t);
  double circ = circumradius(t);
  if(circ == 0)
    return 0;
  else
    return (3 * CGAL::to_double(inradius) / circ);
}

template <typename Triangulation>
void output_distribution_to_png(Triangulation& tr,
                                const int number_of_classes = 100)
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
          double q = CGAL::to_double(radius_ratio(t));
          qualities.push_back(q);
          max_quality = std::max(max_quality, q);
        }
    }

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

//   QPixmap* pixmap = new QPixmap( distribution_x, 
//                                  distribution_y);
  CGAL::Qt_widget *widget = new CGAL::Qt_widget();
  qApp->setMainWidget(widget);
  widget->resize(distribution_x, distribution_y);
  widget->set_window(-0.2, 1.2, -0.2, 1.2); // x_min, x_max, y_min, y_max.
  widget->show();
  
  widget->lock();
  *widget << CGAL::RED << Segment(Point(0.,0.), Point(1., 0.))
          << CGAL::BLACK;
  
  if( number_of_classes == 0 ) return;
  const double width = 1.0 / number_of_classes;

  *widget << CGAL::FillColor(CGAL::BLACK);
  *widget << Segment(Point(0., 0.), Point(1., 0.));
  for(int k=0;k<number_of_classes;k++)
    {
      double height = ( distribution[k]+0. ) / max_occurrence;
      *widget << Rectangle(Point(k*width, 0),
                           Point((k+1)*width, height));
    }
  std::cerr << "Max: " << max_quality << std::endl;
  
  widget->unlock();
  qApp->exec();
//   if( pixmap->save( QString(distribution_filename.c_str()), "PNG") )
//     std::cerr << "Distribution saved to file " << distribution_filename
//               << std::endl;
//   else
//     {
//       std::cerr << "Error: cannot save distribution to file "
//                 << distribution_filename << std::endl; 
//       exit(1);
//     }
}


/////////////// Main function /////////////// 

int main(int argc, char **argv) {

  QApplication app (argc, argv);
  
  parse_argv(argc, argv);

  // Function
  Func F;

  // Oracle (NB: parity oracle is toggled)
  Oracle O (F, K::Point_3 (0,0,0),
	    enclosing_sphere_radius,
	    precision,
	    bipolar_oracle);

  // 2D-complex in 3D-Delaunay triangulation
  Del T;

  // Initial point sample
  Oracle::Points initial_point_sample = 
    O.random_points (number_of_initial_points);
  T.insert (initial_point_sample.begin(), initial_point_sample.end());
  
  // Meshing criteria
  CGAL::Chew_4_surfaces::Curvature_size_criterion<Del> 
    c_s_crit (curvature_bound);
//   CGAL::Chew_4_surfaces::Uniform_size_criterion<Del>
//     u_s_crit (size_bound);
  CGAL::Chew_4_surfaces::Aspect_ratio_criterion<Del>
    a_r_crit (aspect_ratio_bound);

  std::vector<Criterion*> crit_vect;
  crit_vect.push_back (&c_s_crit);
//   crit_vect.push_back (&u_s_crit);
//   crit_vect.push_back(&a_r_crit);
  Criteria C (crit_vect);

  Tets_criteria tets_criteria(tets_aspect_ratio_bound, tets_size_bound);

  std::cerr << "Initial number of points: " << T.number_of_vertices() 
            << std::endl;
  
  // Surface meshing
  Mesher mesher (T, O, C, tets_criteria);
  mesher.refine_surface();
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
  if( dump_distribution )
    output_distribution_to_png(T, distribution_size);
  std::cerr << " done\n";
}
