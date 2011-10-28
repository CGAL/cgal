#ifndef CGAL_KD_CARTESIAN_PER_DIM_H
#define CGAL_KD_CARTESIAN_PER_DIM_H
#include <CGAL/functor_tags.h>
#include <CGAL/Dimension.h>
#include <CGAL/predicates/sign_of_determinant.h>
#include <CGAL/Kernel_d/Cartesian_2.h>

namespace CGAL {
template <class Dim_, class R_, class Derived_>
struct Cartesian_per_dimension : public R_ {};

// TODO: we want to put the actual functors in some other file. Maybe name the
// classes Orientation_2 to make it easy to typedef them all.
template <class R_, class Derived_>
struct Cartesian_per_dimension<Dimension_tag<2>,R_,Derived_> : public R_ {
	typedef R_ Kernel_base;
	template<class F, class D=void> struct Functor
		: Kernel_base::template Functor<F,D>
	{ };

#define CGAL_override_in_dim2(X) \
	template<class D> struct Functor<X##_tag,D> { \
		typedef CartesianDKernelFunctors::X##_2<Derived_> type; \
	}
	CGAL_override_in_dim2(Orientation_of_points);
};
}

#endif
