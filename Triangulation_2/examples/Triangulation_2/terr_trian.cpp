// Delaunay Triangulation of a set of 3D points in the xy-plane.
// (Terrain triangulation)

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/IO/Verbose_ostream.h>
#include <CGAL/IO/OFF.h>
#include <CGAL/boost/graph/IO/OFF.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/boost/graph/graph_traits_Delaunay_triangulation_2.h>
#include <CGAL/boost/graph/graph_traits_Triangulation_2.h>

#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Point_3<Kernel>  Point_3;

typedef CGAL::Projection_traits_xy_3<Kernel>               Gt;

typedef  CGAL::Triangulation_2<Gt>                 Triangulation;
typedef  CGAL::Delaunay_triangulation_2<Gt>        Delaunay_triangulation;

bool  verbose      = false;
bool  delaunay     = false;
bool  incr         = false;

// main function with standard unix commandline arguments
// ------------------------------------------------------
int main( int argc, char **argv) {
  int n = 0; // number of filenames
  char *filename[2];
  bool help = false;
  for (int i = 1; i < argc; i++) { // check commandline options
    if ( strcmp( "-v", argv[i]) == 0)
      verbose = true;
    else if ( strcmp( "-delaunay", argv[i]) == 0)
      delaunay = true;
    else if ( strcmp( "-incr", argv[i]) == 0)
      incr = true;
    else if ( (strcmp( "-h", argv[i]) == 0) ||
              (strcmp( "-help", argv[i]) == 0))
      help = true;
    else if ( n < 2 ) {
      filename[ n++] = argv[i];
    } else {
      ++n;
      break;
    }
  }
  if ((n > 2) || help) {
    if ( ! help)
      std::cerr << "Error: in parameter list" << std::endl;
    std::cerr << "Usage: " << argv[0] << " [<options>] [<infile> [<outfile>]]"
         << std::endl;
    std::cerr << "       Terrain triangulation in the xy-plane." << std::endl;
    std::cerr << "       -delaunay  Delaunay triangulation (triangulation otherwise)." << std::endl;
    std::cerr << "       -incr      Incremental insertion (no flips)." << std::endl;
    std::cerr << "       -v         verbose." << std::endl;
    exit( ! help);
  }

  CGAL::Verbose_ostream vout( verbose);
  vout << argv[0] << ": verbosity on." << std::endl;

  const char*  iname = "cin";
  std::istream* p_in  = &std::cin;
  std::ifstream in;
  if ( n > 0) {
      in.open( filename[0]);
      p_in = &in;
      iname = filename[0];
  }
  if ( !*p_in) {
      std::cerr << argv[0] << ": error: cannot open file '" << iname
       << "' for reading." << std::endl;
      exit(1);
  }

  CGAL::File_scanner_OFF scanner( * p_in, true);
  if ( !*p_in)
      exit(1);

  const char*  oname = "cout";
  std::ostream* p_out = &std::cout;
  std::ofstream out;
  if ( n > 1) {
      out.open( filename[1]);
      p_out = &out;
      oname = filename[1];
  }
  if ( !*p_out) {
      std::cerr << argv[0] << ": error: cannot open file '"<< oname
           << "' for writing." << std::endl;
      exit(1);
  }

  if ( delaunay)
  {
    Delaunay_triangulation triang;
    if (incr)
    {
      vout << "Scanning and triangulating ... " << std::flush;
      for ( std::size_t j = 0; j < scanner.size_of_vertices(); j++) {
          double x, y, z;
          scanner.scan_vertex( x, y, z);
          Point_3 p( x, y, z);
          triang.insert( p);
      }
    }
    else
    {
      vout << "Scanning ... " << std::flush;
      std::vector<Point_3> points;
      for ( std::size_t j = 0; j < scanner.size_of_vertices(); j++) {
          double x, y, z;
          scanner.scan_vertex( x, y, z);
          points.push_back(Point_3(x, y, z));
      }
      vout << "done." << std::endl;
      vout << "Triangulating ... " << std::flush;
      triang.insert(points.begin(), points.end());
    }
    vout << " done." << std::endl;
    vout << "write_triangulation(" << oname << ") ... " << std::flush;
    CGAL::IO::write_OFF(*p_out, triang);
    vout << "done." << std::endl;

  } else {
    Triangulation triang;
    if(incr)
    {
      vout << "Scanning and triangulating ... " << std::flush;
      for ( std::size_t j = 0; j < scanner.size_of_vertices(); j++) {
          double x, y, z;
          scanner.scan_vertex( x, y, z);
          Point_3 p( x, y, z);
          triang.insert( p);
      }
    }
    else
    {
      vout << "Scanning ... " << std::flush;
      std::vector<Point_3> points;
      for ( std::size_t j = 0; j < scanner.size_of_vertices(); j++) {
          double x, y, z;
          scanner.scan_vertex( x, y, z);
          points.push_back(Point_3( x, y, z));
      }
      vout << "done." << std::endl;
      vout << "Triangulating ... " << std::flush;
      triang.insert(points.begin(), points.end());
    }
    vout << "done." << std::endl;
    vout << "write_triangulation(" << oname << ") ... " << std::flush;
    CGAL::IO::write_OFF(*p_out, triang);
    vout << "done." << std::endl;
  }
  if ( !*p_in) {
      std::cerr << argv[0] << " read error: while reading file '"<< iname << "'."
                << std::endl;
      std::exit(1);
  }
  if ( !*p_out) {
      std::cerr << argv[0] << " write error: while writing file '"<< oname << "'."
                << std::endl;
      exit(1);
  }

  return 0;
}
