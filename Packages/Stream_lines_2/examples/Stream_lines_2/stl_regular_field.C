#include <iostream>
#include <fstream>

#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>

#include <CGAL/Stream_lines_2.h>
#include <CGAL/Euler_integrator_2.h>
#include <CGAL/Runge_kutta_integrator_2.h>
#include <CGAL/Regular_grid_2.h>
#include <CGAL/Triangular_field_2.h>

typedef double coord_type;
typedef CGAL::Cartesian<coord_type> K;
typedef CGAL::Regular_grid_2<K> Field;
typedef CGAL::Euler_integrator_2<Field> Euler_integrator;
typedef CGAL::Runge_kutta_integrator_2<Euler_integrator, Field> Runge_kutta_integrator;
typedef CGAL::Stream_lines_2<Field, Runge_kutta_integrator> Stl;
typedef CGAL::Stream_lines_2<Field, Runge_kutta_integrator>::Stream_line_iterator_2 Stl_iterator;

int main(int argc, char **argv)
{
  Euler_integrator euler_integrator(1);  
  Runge_kutta_integrator runge_kutta_integrator(euler_integrator);
  
  std::ifstream infile("data/vnoise.vec.cin", std::ios::in);
  double iXSize, iYSize;
  iXSize = iYSize = 512;
  Field regular_grid_2(infile, iXSize, iYSize);
  infile.close();
  
  /* the placement of streamlines */  
  std::cout << "processing...\n";
  double dSep = 3.5;
  double dRat = 1.6;
  Stl Stream_lines(regular_grid_2, runge_kutta_integrator,dSep,dRat);
  std::cout << "placement generated\n";
  
  Stl_iterator begin_iterator = Stream_lines.begin();
  Stl_iterator end_iterator = Stream_lines.end();
  
	/*writing streamlines to streamlines_on_regular_grid.stl */
  std::cout << "streamlines_on_regular_grid.stl\n";
  std::ofstream fw("streamlines_on_regular_grid.stl",std::ios::out);
  Stream_lines.print_stream_lines(fw);
  
}
