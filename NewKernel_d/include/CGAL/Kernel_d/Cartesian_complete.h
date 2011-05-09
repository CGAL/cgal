#ifndef CGAL_KERNEL_D_CARTESIAN_COMPLETE_H
#define CGAL_KERNEL_D_CARTESIAN_COMPLETE_H

#include <CGAL/Kernel_d/function_objects_cartesian.h>
#include <CGAL/Kernel_d/Segmentd.h>
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

template<class R_,bool force_=false> struct Cartesian_complete_types 
: public R_
{
	typedef R_ Kernel_base;
#define CGAL_Kernel_obj2(X) typedef typename boost::mpl::if_c<force_||!internal::has_##X<R_>::value,X##Cd<R_>,R_>::type::X X;
#include <CGAL/Kernel_d/interface_macros.h>
};


template<class R_,bool force_=false,class Derived_=R_> struct Cartesian_complete_constructors 
: public R_
{
	typedef R_ Kernel_base;
	template<class F> struct Construct {
		typedef typename Kernel_base::template Construct<F>::type Base_functor;
		typedef typename boost::mpl::if_<boost::is_same<Base_functor,Null_functor>,
			typename Cartesian_complete_constructors<R_,true,Derived_>::template Construct<F>::type,
			Base_functor>::type type;
	};

};

template<class R_,class Derived_> struct Cartesian_complete_constructors <R_,true,Derived_>
: public R_
{
	typedef R_ Kernel_base;
	template<class F, int i=0> struct Construct
		: Kernel_base::template Construct<F>
	{ };
#define CGAL_Kernel_cons2(F,f) \
	template<int i> struct Construct <F##_tag,i> { \
		typedef CartesianDKernelFunctors::F<Derived_> type; \
		type f() const { return type(); } \
	};
#include <CGAL/Kernel_d/interface_macros.h>

};


template<class R_,bool force_=false,class Derived_=R_> struct Cartesian_complete_predicates 
: public R_
{
	typedef R_ Kernel_base;
	template<class F> struct Predicate {
		typedef typename Kernel_base::template Predicate<F>::type Base_functor;
		typedef typename boost::mpl::if_<boost::is_same<Base_functor,Null_functor>,
			typename Cartesian_complete_predicates<R_,true,Derived_>::template Predicate<F>::type,
			Base_functor>::type type;
	};

};

template<class R_,class Derived_> struct Cartesian_complete_predicates <R_,true,Derived_>
: public R_
{
	typedef R_ Kernel_base;
	template<class F, int i=0> struct Predicate
		: Kernel_base::template Predicate<F>
	{ };
#define CGAL_Kernel_pred(F,f) \
	template<int i> struct Predicate <F##_tag,i> { \
		typedef CartesianDKernelFunctors::F<Derived_> type; \
		type f() const { return type(); } \
	};
#include <CGAL/Kernel_d/interface_macros.h>

};


template<class R_,bool force_=false,class Derived_=R_> struct Cartesian_complete_computes 
: public R_
{
	typedef R_ Kernel_base;
	template<class F> struct Compute {
		typedef typename Kernel_base::template Compute<F>::type Base_functor;
		typedef typename boost::mpl::if_<boost::is_same<Base_functor,Null_functor>,
			typename Cartesian_complete_computes<R_,true,Derived_>::template Compute<F>::type,
			Base_functor>::type type;
	};

};

template<class R_,class Derived_> struct Cartesian_complete_computes <R_,true,Derived_>
: public R_
{
	typedef R_ Kernel_base;
	template<class F, int i=0> struct Compute
		: Kernel_base::template Compute<F>
	{ };
#define CGAL_Kernel_comp2(F,f) \
	template<int i> struct Compute <F##_tag,i> { \
		typedef CartesianDKernelFunctors::F<Derived_> type; \
		type f() const { return type(); } \
	};
#include <CGAL/Kernel_d/interface_macros.h>

};
}

#endif // CGAL_KERNEL_D_CARTESIAN_COMPLETE_H
