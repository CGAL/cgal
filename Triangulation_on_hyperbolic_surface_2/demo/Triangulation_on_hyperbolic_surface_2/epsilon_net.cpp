#include "window.h"

#include <CGAL/Delaunay_triangulation_on_hyperbolic_surface_2.h>
#include <CGAL/Hyperbolic_Dirichlet_domain_2.h>
#include <CGAL/Hyperbolic_fundamental_domain_factory_2.h>

#include <iostream>
#include <fstream>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_rational.h>

#include <CGAL/Algebraic_kernel_for_circles_2_2.h>
#include <CGAL/Circular_kernel_2.h>
#include <CGAL/Circular_kernel_2/Intersection_traits.h>
#include <CGAL/Hyperbolic_Delaunay_triangulation_CK_traits_2.h>

#include <CGAL/Timer.h>

typedef CGAL::Exact_rational		NumberType;
typedef CGAL::Circular_kernel_2<CGAL::Simple_cartesian<NumberType>,CGAL::Algebraic_kernel_for_circles_2_2<NumberType>> Kernel;
typedef CGAL::Hyperbolic_Delaunay_triangulation_CK_traits_2<Kernel>                                             ParentTraits;
typedef CGAL::Hyperbolic_surface_traits_2<ParentTraits>                                                        	Traits;
typedef typename Traits::Complex                                                                                ComplexNumber;
typedef typename Traits::Hyperbolic_point_2                                                                     Point;

typedef CGAL::Hyperbolic_fundamental_domain_2<Traits>                                                           Domain;
typedef CGAL::Hyperbolic_isometry_2<Traits>                                                                     Isometry;
typedef CGAL::Hyperbolic_fundamental_domain_factory_2<Traits>                                                   Factory;

typedef CGAL::Triangulation_on_hyperbolic_surface_2<Traits, CGAL::Delaunay_triangulation_attributes<Traits>>    Base;
typedef CGAL::Delaunay_triangulation_on_hyperbolic_surface_2<Traits>                                            Delaunay_triangulation;
typedef typename Delaunay_triangulation::Anchor                                                                 Anchor;

int main(int argc, char *argv[])
{
	// 1. GENERATE THE INPUT
	Domain domain;
	int seed;
	if (argc <= 2) {
		seed = time(NULL);
		std::cout << "Generating surface with random seed " << seed << "..." << std::endl;
	} else {
		seed = atoi(argv[2]);
		std::cout << "Generating surface with seed " << seed << "..." << std::endl;
	}
	Factory factory;
	domain = factory.make_hyperbolic_fundamental_domain_g2(seed);

	// 2. GET A VERTEX
	// So that if you run the demo on a same surface but with different values of epsilon,
	// the drawing will be centered at the same vertex and it will look similar.
	Delaunay_triangulation my_triangulation = Delaunay_triangulation(domain);
	my_triangulation.make_Delaunay();
	Anchor & my_anchor = my_triangulation.anchor();
	Point v0 = my_anchor.vertices[0];

	// 3. COMPUTE EPSILON-NET and display useful info
	double eps = 0.1;
	if (argc > 1) {
		eps = std::stod(argv[1]);
	}
	
	std::cout << "Computing a " << eps << "-net..." << std::endl;
	CGAL::Timer timer;
	timer.start();
	std::cout << my_triangulation.epsilon_net(eps) << std::endl;
	timer.stop();
	std::cout << "Done in " << timer.time() << " seconds." << std::endl;
	my_triangulation.combinatorial_map().display_characteristics(std::cout) << std::endl;

	// 4. SET THE FIRST ANCHOR OF THE DRAWING
	Anchor anch = std::get < 0 > (my_triangulation.locate_visibility_walk(v0));
	int index = 0;
	for (int i = 0; i < 3; i++) {
		if (v0 == anch.vertices[i]) {
			index = i;
		}
	}

	Anchor start = Anchor();
	start.dart = anch.dart;
	for (int i = 0; i < 3; i++) {
		start.vertices[i] = anch.vertices[(i + index) % 3];
		if (i < index) {
			start.dart = my_triangulation.Base::ccw(start.dart);
		}
	}

	// 5. DRAW
	QApplication app(argc, argv);
	app.setApplicationName("Demo: epsilon-net");
	DemoWindow window;
	window.item().draw_triangulation(my_triangulation, start);
	window.show();
	QStringList args = app.arguments();
	args.removeAt(0);
	return app.exec();
}
