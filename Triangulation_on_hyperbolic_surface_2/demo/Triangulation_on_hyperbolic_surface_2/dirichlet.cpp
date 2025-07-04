#include "window.h"

#include <CGAL/Delaunay_triangulation_on_hyperbolic_surface_2.h>
#include <CGAL/Hyperbolic_Dirichlet_domain_2.h>
#include <CGAL/Hyperbolic_fundamental_domain_factory_2.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_rational.h>

#include <CGAL/Algebraic_kernel_for_circles_2_2.h>
#include <CGAL/Circular_kernel_2.h>
#include <CGAL/Circular_kernel_2/Intersection_traits.h>
#include <CGAL/Hyperbolic_Delaunay_triangulation_CK_traits_2.h>

typedef CGAL::Exact_rational		NumberType;
typedef CGAL::Circular_kernel_2<CGAL::Simple_cartesian<NumberType>,CGAL::Algebraic_kernel_for_circles_2_2<NumberType>> Kernel;
typedef CGAL::Hyperbolic_Delaunay_triangulation_CK_traits_2<Kernel>                                             ParentTraits;
typedef CGAL::Hyperbolic_surface_traits_2<ParentTraits>                                                        	Traits;
typedef typename Traits::Complex                                    											Complex;
typedef typename Traits::Hyperbolic_point_2                                                                     Point;

typedef CGAL::Hyperbolic_fundamental_domain_2<Traits>                                                           Domain;
typedef CGAL::Hyperbolic_isometry_2<Traits>                                                                     Isometry;
typedef CGAL::Hyperbolic_fundamental_domain_factory_2<Traits>                                                   Factory;

typedef CGAL::Triangulation_on_hyperbolic_surface_2<Traits, CGAL::Delaunay_triangulation_attributes<Traits>>    Base;
typedef CGAL::Delaunay_triangulation_on_hyperbolic_surface_2<Traits>                                            Delaunay_triangulation;
typedef typename Delaunay_triangulation::CMap 																	CMap;
typedef typename Delaunay_triangulation::Anchor                                                                 Anchor;

int main(int argc, char **argv)
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
	Delaunay_triangulation dt = Delaunay_triangulation(domain);
	Point v0 = dt.anchor().vertices[0];

	// 3. DRAW DIRICHLET DOMAIN
	QApplication app(argc, argv);
	app.setApplicationName("Demo: Dirichlet domain and epsilon-net");
	DemoWindow window;
	window.item().draw_dirichlet(domain);

	// 3. COMPUTE EPSILON-NET and display useful info
	double eps = 0.1;
	if (argc > 1) {
		eps = std::stod(argv[1]);
	}

	std::cout << "Computing a " << eps << "-net..." << std::endl;
	std::cout << dt.epsilon_net(eps) << std::endl;
	CMap & cmap = dt.combinatorial_map();
	cmap.display_characteristics(std::cout) << std::endl;

	// // 4. SET THE FIRST ANCHOR OF THE DRAWING
	// // To center the drawing s.t. v0 is translated at 0
	Anchor anchor = dt.locate(v0);
	int index = 0;
	for (int i = 0; i < 3; i++) {
		if (v0 == anchor.vertices[i]) {
			index = i;
		}
	}

	Isometry center_v0 = CGAL::hyperbolic_translation < Traits > (v0);
	Anchor start = Anchor();
	start.dart = anchor.dart;
	for (int i = 0; i < 3; i++) {
		start.vertices[i] = center_v0.evaluate(anchor.vertices[(i + index) % 3]);
		if (i < index) {
			start.dart = dt.Base::ccw(start.dart);
		}
	}

	// 5. DRAW the triangulation
	// using a BFS algo to explore the triangles
	std::queue < Anchor > bfs_queue;
	std::vector < Anchor > to_draw;

	size_t in_queue = cmap.get_new_mark();  // mark darts of triangles with an anchor in the queue
	cmap.unmark_all(in_queue);
	bfs_queue.push(start);
	cmap.mark(start.dart, in_queue);
	cmap.mark(dt.Base::ccw(start.dart), in_queue);
	cmap.mark(dt.Base::cw(start.dart), in_queue);

	while (!bfs_queue.empty()) {
		Anchor & current = bfs_queue.front();
		to_draw.push_back(current);
		auto invader = current.dart;
		for (int i = 0; i < 3; i++) {
			auto invaded = dt.Base::opposite(invader);
			if (!cmap.is_marked(invaded, in_queue)) {
				Complex cross_ratio = dt.Base::get_cross_ratio(invader);
				Point & c = current.vertices[i % 3];
				Point & a = current.vertices[(i + 1) % 3];
				Point & b = current.vertices[(i + 2) % 3];
				Point d =
				    dt.Base::fourth_point_from_cross_ratio(a, b, c, cross_ratio);
				bfs_queue.push(Anchor(invaded, a, c, d));
				cmap.mark(invaded, in_queue);
				cmap.mark(dt.Base::ccw(invaded), in_queue);
				cmap.mark(dt.Base::cw(invaded), in_queue);
			}
			invader = dt.Base::ccw(invader);
		}
		bfs_queue.pop();
	}
	cmap.free_mark(in_queue);

	window.item().draw_triangles(to_draw);

	// 	// 4. SET THE FIRST ANCHOR OF THE DRAWING
	// Anchor anchor = dt.locate(v0);
	// int index = 0;
	// for (int i = 0; i < 3; i++) {
	// 	if (v0 == anchor.vertices[i]) {
	// 		index = i;
	// 	}
	// }

	// Anchor start = Anchor();
	// start.dart = anchor.dart;
	// for (int i = 0; i < 3; i++) {
	// 	start.vertices[i] = anchor.vertices[(i + index) % 3];
	// 	if (i < index) {
	// 		start.dart = dt.Base::ccw(start.dart);
	// 	}
	// }

	// window.item().draw_triangulation(dt, start);
	window.show();
	QStringList args = app.arguments();
	args.removeAt(0);
	return app.exec();
}