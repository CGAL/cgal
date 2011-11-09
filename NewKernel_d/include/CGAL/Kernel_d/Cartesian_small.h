#ifndef CGAL_KD_CARTESIAN_SMALL_H
#define CGAL_KD_CARTESIAN_SMALL_H
#include <CGAL/functor_tags.h>
#include <CGAL/predicates/sign_of_determinant.h>

// WARNING: not written yet, I just wanted to stick the BOOST stuff somewhere not to lose it.

namespace CGAL {
namespace CartesianDKernelFunctors {
template <class R_>
struct Orientation_of_points_small : private Store_kernel<R_> {
	CGAL_FUNCTOR_INIT_STORE(Orientation_of_points_small);
	typedef typename R_::template Type<Point_tag>::type Point;
	typedef typename R_::Orientation result_type;
	typedef typename R_::FT FT;
#if 0
	template<class Iter>
	result_type operator()(Iter f, Iter e)const{
		typename R_::template Functor<Compute_cartesian_coordinate_tag>::type c(this->kernel());
		Point const&A=*f;
		FT const& ax=c(A,0);
		FT const& ay=c(A,1);
		Point const&B=*++f;
		Point const&C=*++f;
		CGAL_assertion(++f==e);
		return sign_of_determinant(c(B,0)-ax,c(B,1)-ay,c(C,0)-ax,c(C,1)-ay);
	}
#endif

// TODO: iterate outside the class to have a single declaration of the right dimension.
#define VAR(Z,J,I) m(I,J)=c(p##I,J)-c(x,J);
#define VAR2(Z,I,N) BOOST_PP_ENUM(N,VAR,I)
#define CODE(Z,N,_) \
        result_type operator()(Point const&x, BOOST_PP_ENUM_PARAMS(N,Point const&p)) const { \
                typename R::template Functor<Compute_cartesian_coordinate_tag>::type c(this->kernel()); \
                Matrix m(N,N); \
                BOOST_PP_REPEAT(N,VAR2,N) \
                return sign_of_determinant(BOOST_PP_ENUM(N,VAR2,N)); \
        }

BOOST_PP_REPEAT_FROM_TO(2, 10, CODE, _ )
#undef CODE
#undef VAR2
#undef VAR
#endif

};
}
}

#endif
