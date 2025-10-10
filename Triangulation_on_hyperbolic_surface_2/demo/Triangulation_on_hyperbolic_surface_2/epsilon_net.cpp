#include "window.h"

#include <CGAL/Exact_rational.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Circular_kernel_2.h>
#include <CGAL/Algebraic_kernel_for_circles_2_2.h>
#include <CGAL/Hyperbolic_Delaunay_triangulation_CK_traits_2.h>
#include <CGAL/Hyperbolic_surface_traits_2.h>

#include <CGAL/Delaunay_triangulation_on_hyperbolic_surface_2.h>
#include <CGAL/Hyperbolic_fundamental_domain_factory_2.h>

#include <iostream>
#include <fstream>
#include <CGAL/Triangulation_on_hyperbolic_surface_2_IO.h>

#include <CGAL/Timer.h>

// typedef CGAL::Exact_rational																					NumberType;
typedef CGAL::Gmpq																								NumberType;
typedef CGAL::Circular_kernel_2<CGAL::Simple_cartesian<NumberType>,CGAL::Algebraic_kernel_for_circles_2_2<NumberType>> Kernel;
typedef CGAL::Hyperbolic_Delaunay_triangulation_CK_traits_2<Kernel>                                             ParentTraits;
typedef CGAL::Hyperbolic_surface_traits_2<ParentTraits>                                                        	Traits;
typedef typename Traits::Hyperbolic_point_2                                                                     Point;

typedef CGAL::Hyperbolic_fundamental_domain_2<Traits>                                                           Domain;
typedef CGAL::Hyperbolic_isometry_2<Traits>                                                                     Isometry;
typedef CGAL::Hyperbolic_fundamental_domain_factory_2<Traits>                                                   Factory;

typedef CGAL::Triangulation_on_hyperbolic_surface_2<Traits, CGAL::Delaunay_triangulation_attributes<Traits>>    Base;
typedef CGAL::Delaunay_triangulation_on_hyperbolic_surface_2<Traits>                                            Delaunay_triangulation;
typedef typename Delaunay_triangulation::Anchor                                                                 Anchor;


double eps = 0.1;
int seed = time(NULL);
int p = 0;

// void parse_command_line(int argc, char ** argv) {
// 	for (int i = 1; i < argc; i++) {
// 		if (!strcmp(argv[i], "--epsilon")) {
// 			i++;
// 			eps = std::stod(argv[i]);
// 		} else if (!strcmp(argv[i], "--seed")) {
// 			i++;
// 			seed = atoi(argv[i]);
// 		} else if(!strcmp(argv[i], "--precision")) {
// 			i++;
// 			p = atoi(argv[i]);
// 		} else {
// 			std::cout << "ERROR: Unknown option: " << argv[i] << std::endl;
// 			exit (1);
// 		}
// 	}
// }

int main(int argc, char *argv[])
{
	// parse_command_line(argc, argv);
	if (argc > 1) {
		eps = std::stod(argv[1]);
	}

	if (argc > 2) {
		p = atoi(argv[2]);
	}

	// 1. GENERATE THE INPUT
	Domain domain;
	if (argc <= 3) {
		seed = time(NULL);
		std::cout << "Generating surface with random seed " << seed << "..." << std::endl;
	} else {
		seed = atoi(argv[3]);
	}
	Factory factory;
	if (seed >= 0) {
		std::cout << "Generating surface with seed " << seed << "..." << std::endl;
		domain = factory.make_hyperbolic_fundamental_domain_g2(seed);
	} else {
		switch (seed) {
			// case -1:
			// 	std::cout << "Loading the surface with a very small systole..." << std::endl;
			// 	Base t;
			// 	std::ifstream("/home/clanuel/Documents/camille/cgal_camille/benchmarks/input_triangulations/dt_thin_surface.txt") >> t;
			// 	dt = Delaunay_triangulation(t);
			case -5:
				std::cout << "Loading the genus 5 surface..." << std::endl;
				std::ifstream("/home/clanuel/Documents/camille/cgal_camille/benchmarks/input_domains/FMgenus5.txt") >> domain;
				break;
			case -7:
				std::cout << "Loading the genus 7 surface..." << std::endl;
				std::ifstream("/home/clanuel/Documents/camille/cgal_camille/benchmarks/input_domains/FM-genus-7.txt") >> domain;
				break;
			case -31:
				std::cout << "Loading the genus 3 surface (1)..." << std::endl;
				std::ifstream("/home/clanuel/Documents/camille/cgal_camille/benchmarks/input_domains/FM-genus-3.1.txt") >> domain;
				break;
			case -32:
				std::cout << "Loading the genus 3 surface (2)..." << std::endl;
				std::ifstream("/home/clanuel/Documents/camille/cgal_camille/benchmarks/input_domains/FM-genus-3.2.txt") >> domain;
				break;
			case -33:
				std::cout << "Loading the genus 3 surface (3)..." << std::endl;
				std::ifstream("/home/clanuel/Documents/camille/cgal_camille/benchmarks/input_domains/FM-genus-3.3.txt") >> domain;
				break;
			default:
				exit (1);
		}
	}
	Delaunay_triangulation dt = Delaunay_triangulation(domain);

	// 2. GET A VERTEX
	// So that if you run the demo on a same surface but with different values of epsilon,
	// the drawing will be centered at the same vertex and it will look similar.
	Point v0 = dt.anchor().vertices[0];

	// 3. COMPUTE EPSILON-NET and display useful info
	if constexpr(!std::is_same<NumberType, CGAL::Gmpq>::value) {
		std::cout << "WARNING: Not using the CGAL::Gmpq number type. Precision will be ignored and to_double approximation will be used instead." << std::endl;
	}
	std::cout << "Computing a " << eps << "-net with floating-point precision " << p*53 << "..." << std::endl;
	CGAL::Timer timer;
	timer.start();
	std::cout << "Is epsilon-net? " << dt.epsilon_net(eps, p) << std::endl;
	timer.stop();
	std::cout << "Done in " << timer.time() << " seconds." << std::endl;
	dt.combinatorial_map().display_characteristics(std::cout) << std::endl;
	// std::cout << dt.is_epsilon_covering(eps) << dt.is_epsilon_packing(eps) << std::endl;
	// std::cout << dt.shortest_edge() << std::endl;
	// std::cout << dt.shortest_loop() << std::endl;

	// 4. SET THE FIRST ANCHOR OF THE DRAWING
	Anchor anchor = dt.locate(v0);
	int index = 0;
	for (int i = 0; i < 3; i++) {
		if (v0 == anchor.vertices[i]) {
			index = i;
		}
	}

	Anchor start = Anchor();
	start.dart = anchor.dart;
	for (int i = 0; i < 3; i++) {
		start.vertices[i] = anchor.vertices[(i + index) % 3];
		if (i < index) {
			start.dart = dt.Base::ccw(start.dart);
		}
	}

	// 5. DRAW
	QApplication app(argc, argv);
	app.setApplicationName("Demo: epsilon-net");
	DemoWindow window;
	window.item().draw_triangulation(dt, start);
	window.show();
	QStringList args = app.arguments();
	args.removeAt(0);
	return app.exec();
}
