#include <iostream.h>
#include <fstream.h>

#include <assert.h>

#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>

#include "CGAL/Stream_lines_2.h"
#include "CGAL/Euler_integrator_2.h"
#include "CGAL/Runge_kutta_integrator_2.h"
#include "CGAL/Triangular_field_2.h"

typedef double coord_type;
typedef CGAL::Cartesian<coord_type> K1;
typedef  struct K : public K1 {};
typedef CGAL::Triangular_field_2<K> Field;
typedef CGAL::Euler_integrator_2<Field> Euler_integrator;
typedef CGAL::Runge_kutta_integrator_2<Euler_integrator, Field> Runge_kutta_integrator;
typedef CGAL::Stream_lines_2<Field, Runge_kutta_integrator> Stl;
typedef CGAL::Stream_lines_2<Field, Runge_kutta_integrator>::Stream_line_iterator_2 stl_iterator;

int main(int argc, char **argv)
{
  Euler_integrator euler_integrator(1);
  Runge_kutta_integrator runge_kutta_integrator(euler_integrator);
  
  ifstream infile(argv[1], ios::in);
  Field triangular_field(infile);
  infile.close();
  
  /* the placement of streamlines */  
  std::cout << "processing... " << endl;
  double dSep = (long double) atof(argv[3]);
  double dRat = (long double) atof(argv[4]);
  Stl Stream_lines(triangular_field, runge_kutta_integrator,dSep,dRat);
  std::cout << "success" << endl;
  
  stl_iterator begin_iterator = Stream_lines.begin();
  stl_iterator end_iterator = Stream_lines.end();
  
  char * cName;
  cName = argv[2];
  strcat(cName, ".stl");
  std::cout << cName << endl;
  ofstream fw(cName,ios::out);
  Stream_lines.print_stream_lines(fw);
}
