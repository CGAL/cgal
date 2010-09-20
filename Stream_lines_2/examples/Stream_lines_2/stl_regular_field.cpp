#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Stream_lines_2.h>
#include <CGAL/Runge_kutta_integrator_2.h>
#include <CGAL/Regular_grid_2.h>
#include <CGAL/Triangular_field_2.h>

#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
typedef CGAL::Regular_grid_2<K> Field;
typedef CGAL::Runge_kutta_integrator_2<Field> Runge_kutta_integrator;
typedef CGAL::Stream_lines_2<Field, Runge_kutta_integrator> Strl;
typedef Strl::Point_iterator_2 Point_iterator;
typedef Strl::Stream_line_iterator_2 Strl_iterator;
typedef Strl::Point_2 Point_2;
typedef Strl::Vector_2 Vector_2;

int main()
{
  Runge_kutta_integrator runge_kutta_integrator;

  /*data.vec.cin is an ASCII file containing the vector field values*/
  std::ifstream infile("data/vnoise.vec.cin", std::ios::in);
  double iXSize, iYSize;
  unsigned int x_samples, y_samples;
  iXSize = iYSize = 512;
  infile >> x_samples;
  infile >> y_samples;
  Field regular_grid_2(x_samples, y_samples, iXSize, iYSize);
  /*fill the grid with the appropreate values*/
  for (unsigned int i=0;i<x_samples;i++)
    for (unsigned int j=0;j<y_samples;j++)
      {
        double xval, yval;
        infile >> xval;
        infile >> yval;
        regular_grid_2.set_field(i, j, Vector_2(xval, yval));
      }
  infile.close();

  /* the placement of streamlines */
  std::cout << "processing...\n";
  double dSep = 3.5;
  double dRat = 1.6;
  Strl Stream_lines(regular_grid_2, runge_kutta_integrator,dSep,dRat);
  std::cout << "placement generated\n";

  /*writing streamlines to streamlines_on_regular_grid_1.stl */
  std::ofstream fw("streamlines_on_regular_grid_1.stl",std::ios::out);
  fw << Stream_lines.number_of_lines() << "\n";
  for(Strl_iterator sit = Stream_lines.begin(); sit != Stream_lines.end(); sit++)
    {
      fw << "\n";
      for(Point_iterator pit = sit->first; pit != sit->second; pit++){
	Point_2 p = *pit;
	fw << p.x() << " " << p.y() << "\n";
      }
    }

  fw.close();

}
