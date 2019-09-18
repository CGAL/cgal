#include <CGAL/Frechet_distance.h>
#include <CGAL/Cartesian.h>

#include <vector>

int main()
{
	using Kernel = CGAL::Cartesian<double>;
	using NT = Kernel::FT;
	using Point = Kernel::Point_2;
	using Points = std::vector<Point>;

	Points curve1 = { Point(0., 0.), Point(1., 1.) };
	Points curve2 = { Point(0., 1.), Point(1., 2.) };
	NT distance1 = 2.;
	NT distance2 = .5;

	if (!continuous_Frechet_distance_less_than(curve1, curve2, distance1)) {
		return EXIT_FAILURE;
	}
	if (continuous_Frechet_distance_less_than(curve1, curve2, distance2)) {
		return EXIT_FAILURE;
	}

	auto dist = continuous_Frechet_distance(curve1, curve2);
	if (dist > 1.01 || dist < 0.99) {
		return EXIT_FAILURE;
	}


	return EXIT_SUCCESS;
}
