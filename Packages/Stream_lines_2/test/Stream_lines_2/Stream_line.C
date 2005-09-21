#include <iostream>
#include <fstream>

#include <CGAL/Cartesian.h>

#include <CGAL/Stream_lines_2.h>
#include <CGAL/Euler_integrator_2.h>
#include <CGAL/Runge_kutta_integrator_2.h>
#include <CGAL/Regular_grid_2.h>
#include <CGAL/Triangular_field_2.h>

typedef double coord_type;
typedef CGAL::Cartesian<coord_type> K;
typedef CGAL::Regular_grid_2<K> Field;
typedef CGAL::Runge_kutta_integrator_2<Field> Runge_kutta_integrator;
typedef CGAL::Stream_lines_2<Field, Runge_kutta_integrator> Stl;
typedef Stl::Stream_line_iterator_2 Stl_iterator;

int main()
{
	Runge_kutta_integrator runge_kutta_integrator;
	for (int i=1;i<=23;i++)
	if (i!=11 && i!=12 && i!=17 && i!=18 && i!=20 && i!=21)
	{
		char name[80];
		sprintf(name, "data/%d.vec.cin", i);
		char namer[80];
		sprintf(namer, "data/%d.stl", i);
		std::ifstream infile(name, std::ios::in);
		double iXSize, iYSize;
		iXSize = iYSize = 512;
		Field regular_grid_2(infile, iXSize, iYSize);
		infile.close();
		/* the placement of streamlines */
		double dSep = 3.5;
		double dRat = 1.6;
		std::cout << "processing... (" << name << ")\n";
		Stl Stream_lines(regular_grid_2, runge_kutta_integrator,dSep,dRat,1);
		std::cout << "placement generated\n";
	}
	std::cout << "success\n";
}
