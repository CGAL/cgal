#ifndef CGAL_KERNEL_D_CARTESIAN_COMPLETE_H
#define CGAL_KERNEL_D_CARTESIAN_COMPLETE_H

#include <CGAL/Kernel_d/function_objects_cartesian.h>
#include <CGAL/Kernel_d/Segmentd.h>
#include <CGAL/Kernel_d/Cartesian_per_dimension.h>
#include <boost/mpl/if.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/mpl/has_xxx.hpp>


namespace CGAL {
namespace internal {
BOOST_MPL_HAS_XXX_TRAIT_DEF(Segment)
BOOST_MPL_HAS_XXX_TRAIT_DEF(Ray)
BOOST_MPL_HAS_XXX_TRAIT_DEF(Direction)
BOOST_MPL_HAS_XXX_TRAIT_DEF(Line)
}

template<class R_,class Derived_> struct Cartesian_define_all_functors_
: public R_
{
	typedef R_ Kernel_base;
	template<class F, class D=void> struct Functor
		: Kernel_base::template Functor<F,D>
	{ };
#define CGAL_Kernel_cons2(F,f) \
	template<class D> struct Functor <F##_tag,D> { \
		typedef CartesianDKernelFunctors::F<Derived_> type; \
	};
	//type f() const { return type(); }
#define CGAL_Kernel_obj2(X,Y) \
	template<class D> struct Functor <Construct_ttag<X##_tag>,D> { \
		typedef CartesianDKernelFunctors::Construct_##Y<Derived_> type; \
	};
#define CGAL_Kernel_pred(F,f) CGAL_Kernel_cons2(F,f)
#define CGAL_Kernel_comp2(F,f) CGAL_Kernel_cons2(F,f)

#include <CGAL/Kernel_d/interface_macros.h>

};

template<class R_,class Derived_> struct Cartesian_define_all_functors
: public Cartesian_per_dimension<typename R_::Default_ambient_dimension,Cartesian_define_all_functors_<R_,Derived_>,Derived_> {};

template<class R_,bool force_=false> struct Cartesian_complete_types 
: public R_
{
	typedef R_ Kernel_base;
	template <class T,class=void> struct Type : Kernel_base::template Type<T> {};
#define CGAL_Kernel_obj2(X,Y) \
	template <class D> struct Type<X##_tag,D> { \
		typedef typename Kernel_base::template Type<X##_tag> B_; \
		typedef typename boost::mpl::if_c<force_||!internal::has_type<B_>::value,Wrap_type<X##Cd<R_> >,B_>::type::type type; \
	};
#include <CGAL/Kernel_d/interface_macros.h>
};


template<class R_,bool force_=false,class Derived_=R_> struct Cartesian_complete_constructors 
: public R_
{
	typedef R_ Kernel_base;
	template<class F,class D=void,class=typename map_functor_type<F>::type> struct Functor :
		R_::template Functor<F,D> {};
	template<class F,class D> struct Functor<F,D,Construct_tag> {
		typedef typename Kernel_base::template Functor<F>::type Base_functor;
		typedef typename boost::mpl::if_c<force_||boost::is_same<Base_functor,Null_functor>::value,
			typename Cartesian_define_all_functors<R_,Derived_>::template Functor<F>::type,
			Base_functor>::type type;
	};

};

template<class R_,bool force_=false,class Derived_=R_> struct Cartesian_complete_predicates 
: public R_
{
	typedef R_ Kernel_base;
	template<class F,class D=void,class=typename map_functor_type<F>::type> struct Functor :
		R_::template Functor<F,D> {};
	template<class F,class D> struct Functor<F,D,Predicate_tag> {
		typedef typename Kernel_base::template Functor<F>::type Base_functor;
		typedef typename boost::mpl::if_c<force_||boost::is_same<Base_functor,Null_functor>::value,
			typename Cartesian_define_all_functors<R_,Derived_>::template Functor<F>::type,
			Base_functor>::type type;
	};

};

template<class R_,bool force_=false,class Derived_=R_> struct Cartesian_complete_computes 
: public R_
{
	typedef R_ Kernel_base;
	template<class F,class D=void,class=typename map_functor_type<F>::type> struct Functor :
		R_::template Functor<F,D> {};
	template<class F,class D> struct Functor<F,D,Compute_tag> {
		typedef typename Kernel_base::template Functor<F>::type Base_functor;
		typedef typename boost::mpl::if_c<force_||boost::is_same<Base_functor,Null_functor>::value,
			typename Cartesian_define_all_functors<R_,Derived_>::template Functor<F>::type,
			Base_functor>::type type;
	};

};

}

#endif // CGAL_KERNEL_D_CARTESIAN_COMPLETE_H
