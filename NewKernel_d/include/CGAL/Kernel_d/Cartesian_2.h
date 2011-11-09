#ifndef CGAL_KD_CARTESIAN_2_H
#define CGAL_KD_CARTESIAN_2_H
#include <CGAL/functor_tags.h>
#include <CGAL/predicates/sign_of_determinant.h>

namespace CGAL {
namespace CartesianDKernelFunctors {
template <class R_>
struct Orientation_of_points_2 : private Store_kernel<R_> {
	CGAL_FUNCTOR_INIT_STORE(Orientation_of_points_2);
	typedef typename R_::template Type<Point_tag>::type Point;
	typedef typename R_::Orientation result_type;
	typedef typename R_::FT FT;
	template <class Iter>
	result_type operator() (Iter f, Iter e) const {
		typename R_::template Functor<Compute_cartesian_coordinate_tag>::type c(this->kernel());
		Point const& A = *f;
		FT const& ax = c(A, 0);
		FT const& ay = c(A, 1);
		Point const& B = *++f;
		Point const& C = *++f;
		CGAL_assertion(++f == e);
		return sign_of_determinant(
				c(B, 0) - ax, c(B, 1) - ay,
				c(C, 0) - ax, c(C, 1) - ay);
	}
};
}
}

#endif
