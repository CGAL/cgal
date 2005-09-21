#include <iostream>
#include <fstream>

#include <CGAL/Cartesian.h>

#include <CGAL/Stream_lines_2.h>
#include <CGAL/Euler_integrator_2.h>
#include <CGAL/Runge_kutta_integrator_2.h>
#include <CGAL/Triangular_field_2.h>

typedef double coord_type;
typedef CGAL::Cartesian<coord_type> K;
typedef CGAL::Triangular_field_2<K> Field;
typedef CGAL::Runge_kutta_integrator_2<Field> Runge_kutta_integrator;
typedef CGAL::Stream_lines_2<Field, Runge_kutta_integrator> Stl;
typedef Stl::Stream_line_iterator_2 stl_iterator;

int main()
{
	Runge_kutta_integrator runge_kutta_integrator(1);

	std::ifstream infile("data/irregular_data.tri.cin", std::ios::in);
	Field triangular_field(infile);
	infile.close();

	/* the placement of streamlines */
	std::cout << "processing...\n";
	double dSep = 30.0;
	double dRat = 1.6;
	Stl Stream_lines(triangular_field, runge_kutta_integrator,dSep,dRat);
	std::cout << "placement generated\n";

	/*writing streamlines to streamlines.stl */
	std::cout << "streamlines.stl\n";
	std::ofstream fw("streamlines.stl",std::ios::out);
	Stream_lines.print_stream_lines(fw);
}
