#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Convex_hull_3/dual/halfspace_intersection_3.h>
#include <CGAL/Convex_hull_3/dual/halfspace_intersection_with_constructions_3.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>

#include <cmath>
#include <set>

using Point = std::tuple <signed long, signed long, signed long>;
using K = CGAL::Exact_predicates_exact_constructions_kernel;

constexpr signed char L = 24;
constexpr signed long R = 1000000000000l;

namespace PMP=CGAL::Polygon_mesh_processing;

int main() {
	double memo[4*L];
	for (signed char i = 0; i < L; i++)
		memo[i] = cos(2*M_PI*i/L);
	for (signed char i = 1; i < 4; i++)
		for (signed char j = 0; j < L; j++)
			memo[i*L + j] = memo[j];
	const double* const cm = memo + 2*L;

	std::set <Point> set;
	for (signed char a = 0; a < L; a++)
		for (signed char b = 0; b < L; b++)
			for (signed char c = 0; c < L; c++)
				for (signed char d = 0; d < L; d++) {
					const double x = 2*(cm[a + b] + cm[b + c] + cm[c + d] + cm[d + c - a] + cm[c - a + d - b] + cm[d - b - a]);
					const double y = 2*(cm[a    ] +             cm[c    ] +                 cm[c - a        ]                );
					const double z = 2*(            cm[b    ] +             cm[d        ] +                     cm[d - b    ]);

					set.insert(Point(round(R*x), round(R*y), round(R*z)));
				}

	std::vector <CGAL::Plane_3 <K> > in;
	for (const auto& p:  set)
		in.push_back(K::Plane_3(-get <0> (p), -get <1> (p), -get <2> (p), -R));

  {
	CGAL::Polyhedron_3 <K> hull;
	CGAL::halfspace_intersection_3(in.begin(), in.end(), hull, CGAL::ORIGIN);
  assert(PMP::triangulate_faces(hull));
  assert(!PMP::does_self_intersect(hull));
  }
  {
  //CGAL::Polyhedron_3 <K> hull;
  //CGAL::halfspace_intersection_with_constructions_3(in.begin(), in.end(), hull, CGAL::ORIGIN);
  //assert(PMP::triangulate_faces(hull));
  //assert(!PMP::does_self_intersect(hull));
  }
}
